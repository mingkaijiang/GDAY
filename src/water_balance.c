#include "water_balance.h"
#include "nrutil.h"

void initialise_soils_day(control *c, fluxes *f, params *p, state *s) {
    /* Initialise soil water state & parameters  */

    double *fsoil_top = NULL, *fsoil_root = NULL;
    int     i;

    /* site params not known, so derive them based on Cosby et al */

    if (c->calc_sw_params) {
        fsoil_top = get_soil_fracs(p->topsoil_type);
        fsoil_root = get_soil_fracs(p->rootsoil_type);

        /* top soil */
        calc_soil_params(fsoil_top, &p->theta_fc_topsoil, &p->theta_wp_topsoil,
                         &p->theta_sp_topsoil, &p->b_topsoil,
                         &p->psi_sat_topsoil);

        /* Plant available water in top soil (mm) */
        p->wcapac_topsoil = p->topsoil_depth  *\
                            (p->theta_fc_topsoil - p->theta_wp_topsoil);
        /* Root zone */
        calc_soil_params(fsoil_root, &p->theta_fc_root, &p->theta_wp_root,
                         &p->theta_sp_root, &p->b_root, &p->psi_sat_root);

        /* Plant available water in rooting zone (mm) */
        p->wcapac_root = p->rooting_depth * \
                            (p->theta_fc_root - p->theta_wp_root);
    }


    /* calculate Landsberg and Waring SW modifier parameters if not
       specified by the user based on a site calibration */
    if (p->ctheta_topsoil < -900.0 && p->ntheta_topsoil  < -900.0 &&
        p->ctheta_root < -900.0 && p->ntheta_root < -900.0) {
        get_soil_params(p->topsoil_type, &p->ctheta_topsoil, &p->ntheta_topsoil);
        get_soil_params(p->rootsoil_type, &p->ctheta_root, &p->ntheta_root);
    }
    /*
    printf("%f\n", p->wcapac_topsoil);
    printf("%f\n\n", p->wcapac_root);

    printf("%f\n", p->ctheta_topsoil);
    printf("%f\n", p->ntheta_topsoil);
    printf("%f\n", p->ctheta_root);
    printf("%f\n", p->ntheta_root);
    printf("%f\n", p->rooting_depth);

    exit(1); */

    free(fsoil_top);
    free(fsoil_root);

    return;
}


double calc_net_radiation(params *p, double sw_rad, double tair) {

    double net_rad, net_lw;

    /* Net loss of long-wave radn, Monteith & Unsworth '90, pg 52, eqn 4.17 */
    net_lw = 107.0 - 0.3 * tair;            /* W m-2 */

    /* Net radiation recieved by a surf, Monteith & Unsw '90, pg 54 eqn 4.21
        - note the minus net_lw is correct as eqn 4.17 is reversed in
          eqn 4.21, i.e Lu-Ld vs. Ld-Lu
        - NB: this formula only really holds for cloudless skies!
        - Bounding to zero, as we can't have negative soil evaporation, but you
          can have negative net radiation.
        - units: W m-2
    */
    net_rad = MAX(0.0, (1.0 - p->albedo) * sw_rad - net_lw);

    return (net_rad);
}

double calc_stomatal_conductance(params *p, state *s, double vpd, double Ca,
                                 double gpp) {
    /*
        Calculate stomatal conductance using Belinda's model. For the medlyn
        model this is already in conductance to CO2, so the 1.6 from the
        corrigendum to Medlyn et al 2011 is missing here

    References:
    -----------
    For conversion factor for conductance see...
    * Jones (1992) Plants and microclimate, pg 56 + Appendix 3
    * Diaz et al (2007) Forest Ecology and Management, 244, 32-40.

    Stomatal Model:
    * Medlyn et al. (2011) Global Change Biology, 17, 2134-2144.
    **Note** Corrigendum -> Global Change Biology, 18, 3476.

    Parameters:
    -----------
    g1 : float
        slope
    wtfac : float
        water availability scaler [0,1]
    vpd : float
        vapour pressure deficit (Pa)
    Ca : float
        atmospheric co2 [umol mol-1]
    gpp : float
        photosynthesis at the canopy scale (umol m-2 s-1)

    Returns:
    --------
    gs : float
        stomatal conductance (mol CO2 m-2 s-1)
    */
    double gs_over_a, g1, g0, gsc;

    g1 = p->g1 * s->wtfac_root;
    g0 = 0.0; /* p->g0; */
    gs_over_a = (1.0 + g1 / sqrt(vpd * PA_2_KPA)) / Ca;
    gsc = MAX(g0, g0 + gs_over_a * gpp);

    /* mol m-2 s-1 */
    return (gsc);

}


double canopy_boundary_layer_conduct(params *p, double canht, double wind,
                                     double press, double tair) {
    /*  Canopy boundary layer conductance, ga (from Jones 1992 p 68)


    Notes:
    ------
    'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
    3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
    Drake et al, 2010, 17, pg. 1526.

    References:
    ------------
    * Jones 1992, pg. 67-8.
    * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
      of what is in Monteith (ga = 1 / ra)
    * Allen et al. (1989) pg. 651.
    * Gash et al. (1999) Ag forest met, 94, 149-158.

    Parameters:
    -----------
    params : p
        parameters structure
    canht : float
        canopy height (m)
    wind : float
        wind speed (m s-1)
    press : float
        atmospheric pressure (Pa)
    tair : float
        air temperature (deg C)

    Returns:
    --------
    ga : float
        canopy boundary layer conductance (mol m-2 s-1)
    */

    /* z0m roughness length governing momentum transfer [m] */
    double z0m, z0h, d, arg1, arg2, arg3, ga, cmolar;
    double vk = 0.41;

    /* Convert from mm s-1 to mol m-2 s-1 */
    cmolar = press / (RGAS * (tair + DEG_TO_KELVIN));

    /* roughness length for momentum */
    z0m = p->dz0v_dh * canht;

    /*
       z0h roughness length governing transfer of heat and vapour [m]
      *Heat tranfer typically less efficent than momentum transfer. There is
       a lot of variability in values quoted for the ratio of these two...
       JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
       and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
       the default I am following Monteith and Unsworth, by setting the
       ratio to be 1, the code below is identical to that on page 249,
       eqn 15.7
    */
    z0h = p->z0h_z0m * z0m;

    /* zero plan displacement height [m] */
    d = p->displace_ratio * canht;

    arg1 = (vk * vk) * wind;
    arg2 = log((canht - d) / z0m);
    arg3 = log((canht - d) / z0h);

    ga = (arg1 / (arg2 * arg3)) * cmolar;

    return (ga);
}

double calc_slope_of_sat_vapour_pressure_curve(double tair) {
    /*
        Constant slope in Penman-Monteith equation

        Parameters:
        -----------
        tavg : float
            average daytime temperature

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [Pa K-1]

    */
    double arg1, arg2, slope;

    /* Const slope in Penman-Monteith equation  (Pa K-1) */
    arg1 = calc_sat_water_vapour_press(tair + 0.1);
    arg2 = calc_sat_water_vapour_press(tair);
    slope = (arg1 - arg2) / 0.1;

    return (slope);
}


double calc_pyschrometric_constant(double press, double lambda) {
    /* Psychrometric constant ratio of specific heat of moist air at
    a constant pressure to latent heat of vaporisation.

    Parameters:
    -----------
    press : float
        air pressure (Pa)
    lambda : float
         latent heat of water vaporization (J mol-1)


    Returns:
    --------
    gamma : float
        pyschrometric constant [Pa K-1]

    */
    return ( CP * MASS_AIR * press / lambda );

}

double calc_latent_heat_of_vapourisation(double tair) {
    /*   Latent heat of water vapour at air temperature

        Returns:
        -----------
        lambda : float
            latent heat of water vaporization [J mol-1]
    */
    return ( (H2OLV0 - 2.365E3 * tair) * H2OMW );

}

double calc_sat_water_vapour_press(double tac) {
    /*
        Calculate saturated water vapour pressure (Pa) at
        temperature TAC (Celsius). From Jones 1992 p 110 (note error in
        a - wrong units)
    */
    return (613.75 * exp(17.502 * tac / (240.97 + tac)));
}



double *get_soil_fracs(char *soil_type) {
    /*
     * Based on Table 2 in Cosby et al 1984, page 2.
     * Fractions of silt, sand and clay (in that order)
     */
    double *fsoil = malloc(3 * sizeof(double));

    if (strcmp(soil_type, "sand") == 0) {
        fsoil[0] = 0.05;
        fsoil[1] = 0.92;
        fsoil[2] = 0.03;
    } else if (strcmp(soil_type, "loamy_sand") == 0) {
        fsoil[0] = 0.12;
        fsoil[1] = 0.82;
        fsoil[2] = 0.06;
    } else if (strcmp(soil_type, "sandy_loam") == 0) {
        fsoil[0] = 0.32;
        fsoil[1] = 0.58;
        fsoil[2] = 0.1;
    } else if (strcmp(soil_type, "loam") == 0) {
        fsoil[0] = 0.39;
        fsoil[1] = 0.43;
        fsoil[2] = 0.18;
    } else if (strcmp(soil_type, "silty_loam") == 0) {
        fsoil[0] = 0.7;
        fsoil[1] = 0.17;
        fsoil[2] = 0.13;
    } else if (strcmp(soil_type, "sandy_clay_loam") == 0) {
        fsoil[0] = 0.15;
        fsoil[1] = 0.58;
        fsoil[2] = 0.27;
    } else if (strcmp(soil_type, "clay_loam") == 0) {
        fsoil[0] = 0.34;
        fsoil[1] = 0.32;
        fsoil[2] = 0.34;
    } else if (strcmp(soil_type, "silty_clay_loam") == 0) {
        fsoil[0] = 0.56;
        fsoil[1] = 0.1;
        fsoil[2] = 0.34;
    } else if (strcmp(soil_type, "sandy_clay") == 0) {
        fsoil[0] = 0.06;
        fsoil[1] = 0.52;
        fsoil[2] = 0.42;
    } else if (strcmp(soil_type, "silty_clay") == 0) {
        fsoil[0] = 0.47;
        fsoil[1] = 0.06;
        fsoil[2] = 0.47;
    } else if (strcmp(soil_type, "clay") == 0) {
        fsoil[0] = 0.2;
        fsoil[1] = 0.22;
        fsoil[2] = 0.58;
    } else {
        prog_error("Could not understand soil type", __LINE__);
    }

    return (fsoil);
}

void get_soil_params(char *soil_type, double *c_theta, double *n_theta) {
    /* For a given soil type, get the parameters for the soil
    moisture availability based on Landsberg and Waring, with updated
    parameters from Landsberg and Sands (2011), pg 190, Table 7.1

    Table also has values from Saxton for soil texture, perhaps makes more
    sense to use those than Cosby? Investigate?

    Reference
    ---------
    * Landsberg and Sands (2011) Physiological ecology of forest production.
    * Landsberg and Waring (1997) Forest Ecology & Management, 95, 209-228.
    * Clapp & Hornberger (1978) Water Resources Research, 14, 601–604.
    */

    if (strcmp(soil_type, "clay") == 0) {
        *c_theta = 0.4;
        *n_theta = 3.0;
    } else if (strcmp(soil_type, "clay_loam") == 0) {
        *c_theta = 0.5;
        *n_theta = 5.0;
    } else if (strcmp(soil_type, "loam") == 0) {
        *c_theta = 0.55;
        *n_theta = 6.0;
    } else if (strcmp(soil_type, "loamy_sand") == 0) {
        *c_theta = 0.65;
        *n_theta = 8.0;
    } else if (strcmp(soil_type, "sand") == 0) {
        *c_theta = 0.7;
        *n_theta = 9.0;
    } else if (strcmp(soil_type, "sandy_clay") == 0) {
        *c_theta = 0.45;
        *n_theta = 4.0;
    } else if (strcmp(soil_type, "sandy_clay_loam") == 0) {
        *c_theta = 0.525;
        *n_theta = 5.5;
    } else if (strcmp(soil_type, "sandy_loam") == 0) {
        *c_theta = 0.6;
        *n_theta = 7.0;
    } else if (strcmp(soil_type, "silt") == 0) {
        *c_theta = 0.625;
        *n_theta = 7.5;
    } else if (strcmp(soil_type, "silty_clay") == 0) {
        *c_theta = 0.425;
        *n_theta = 3.5;
    } else if (strcmp(soil_type, "silty_clay_loam") == 0) {
        *c_theta = 0.475;
        *n_theta = 4.5;
    } else if (strcmp(soil_type, "silty_loam") == 0) {
        *c_theta = 0.575;
        *n_theta = 6.5;
    } else {
        prog_error("There are no parameters for your soil type", __LINE__);
    }

    return;
}

void calc_soil_params(double *fsoil, double *theta_fc, double *theta_wp,
                      double *theta_sp, double *b, double *psi_sat_mpa) {
    /* Cosby parameters for use within the Clapp Hornberger soil hydraulics
    scheme are calculated based on the texture components of the soil.

    NB: Cosby et al were ambiguous in their paper as to what log base to
    use.  The correct implementation is base 10, as below.

    Parameters:
    ----------
    fsoil : list
        fraction of silt, sand, and clay (in that order

    Returns:
    --------
    theta_fc : float
        volumetric soil water concentration at field capacity
    theta_wp : float
        volumetric soil water concentration at the wilting point

    */
    /* soil suction of 3.364m and 152.9m, or equivalent of -0.033 & -1.5 MPa */
    double pressure_head_wilt = -152.9;
    double pressure_head_crit = -3.364;
    double psi_sat;

    /* *Note* subtle unit change to be consistent with fractions as opposed
     * to percentages of sand, silt, clay, e.g. I've changed the slope in
     * the "b" Clapp paramter from 0.157 to 15.7
     *
     * Also Cosby is unclear about which log base were used. 'Generally' now
     * assumed that logarithms to the base 10
     */

    /* Clapp Hornberger exponent [-] */
    *b = 3.1 + 15.7 * fsoil[CLAY] - 0.3 * fsoil[SAND];

    /*
     * soil matric potential at saturation, taking inverse of log (base10)
     * units = m
     */
    psi_sat = CM_2_M * -(pow(10.0, (1.54 - 0.95 * fsoil[SAND] +\
              0.63 * fsoil[SILT])));
    *psi_sat_mpa = psi_sat * METER_OF_HEAD_TO_MPA;

    /* volumetric soil moisture concentrations at the saturation point */
    *theta_sp = 0.505 - 0.037 * fsoil[CLAY] - 0.142 * fsoil[SAND];

    /*
     * volumetric soil moisture concentrations at the wilting point
     * assumed to equal suction of -1.5 MPa or a depth of water of 152.9 m
     */
    *theta_wp = *theta_sp * pow((psi_sat / pressure_head_wilt), (1.0 / *b));

    /*
     * volumetric soil moisture concentrations at field capacity assumed to
     * equal a suction of -0.0033 MPa or a depth of water of 3.364 m
     */
    *theta_fc = *theta_sp * pow((psi_sat / pressure_head_crit), (1.0 / *b));

    return;

}

void calculate_soil_water_fac(control *c, params *p, state *s) {
    /* Estimate a relative water availability factor [0..1]

    A drying soil results in physiological stress that can induce stomatal
    closure and reduce transpiration. Further, N mineralisation depends on
    top soil moisture.

    s->qs = 0.2 in SDGVM

    References:
    -----------
    * Landsberg and Waring (1997) Forest Ecology and Management, 95, 209-228.
      See --> Figure 2.
    * Egea et al. (2011) Agricultural Forest Meteorology, 151, 1370-1384.

    But similarly see:
    * van Genuchten (1981) Soil Sci. Soc. Am. J, 44, 892--898.
    * Wang and Leuning (1998) Ag Forest Met, 91, 89-111.

    * Pepper et al. (2008) Functional Change Biology, 35, 493-508

    Returns:
    --------
    wtfac_topsoil : float
        water availability factor for the top soil [0,1]
    wtfac_root : float
        water availability factor for the root zone [0,1]
    */

    double moisture_ratio_topsoil, moisture_ratio_root;
    double b, sf, psi_f;
    /*double psi_swp_topsoil;*/

    if (c->water_balance == BUCKET && c->sw_stress_model == 0) {
        /* JULES type model, see Egea et al. (2011) */
        s->wtfac_topsoil = calc_beta(s->pawater_topsoil, p->topsoil_depth,
                                     p->theta_fc_topsoil, p->theta_wp_topsoil,
                                     p->qs);

        s->wtfac_root = calc_beta(s->pawater_root, p->rooting_depth,
                                     p->theta_fc_root, p->theta_wp_root,
                                     p->qs);

    } else if (c->water_balance == BUCKET && c->sw_stress_model == 1) {
        /* Landsberg and Waring, (1997) */
        moisture_ratio_topsoil = s->pawater_topsoil / p->wcapac_topsoil;
        moisture_ratio_root = s->pawater_root / p->wcapac_root;

        s->wtfac_topsoil = calc_sw_modifier(moisture_ratio_topsoil,
                                            p->ctheta_topsoil,
                                            p->ntheta_topsoil);
        s->wtfac_root = calc_sw_modifier(moisture_ratio_root, p->ctheta_root,
                                         p->ntheta_root);

    } else if (c->water_balance == BUCKET && c->sw_stress_model == 2) {
        /*
            Zhou et al.(2013) Agricultural & Forest Met. 182–183, 204–214
            Assuming that overnight 􏰀pre-dawn leaf water potential =
            pre-dawn soil water potential.
        */
        //fprintf(stderr, "Zhou model not implemented\n");
        //exit(EXIT_FAILURE);

        // Hardwiring this for testing. Values taken from Table, 1 in
        // De Kauwe et al. 2015, Biogeosciences
        b = 0.82;
        sf = 1.9;
        psi_f = -1.85;

        s->wtfac_topsoil = exp(b * s->predawn_swp);
        s->wtfac_root = exp(b * s->predawn_swp);

        //s->wtfac_topsoil_ns = (1.0 + exp(sf * psi_f)) / \
        //                      (1.0 + exp(sf * (psi_f - s->predawn_swp)));
        //s->wtfac_root_ns = (1.0 + exp(sf * psi_f)) / \
        //                      (1.0 + exp(sf * (psi_f - s->predawn_swp)));

        /*
        s->wtfac_topsoil = exp(p->g1_b * s->psi_s_topsoil);
        s->wtfac_root = exp(p->g1_b * s->psi_s_root);

        ! SW modifier for Vcmax (non-stomatal limitation)
        s->wtfac_topsoil_ns = (1.0 + exp(p->vcmax_sf * p->vcmax_psi_f)) / \
                              (1.0 + exp(p->vcmax_sf * \
                                        (p->vcmax_psi_f - s->psi_s_topsoil)));
        s->wtfac_root_ns = (1.0 + exp(p->vcmax_sf * p->vcmax_psi_f)) / \
                            (1.0 + exp(p->vcmax_sf * \
                                        (p->vcmax_psi_f - s->psi_s_root)));
        */
    }
    return;
}

double calc_beta(double paw, double depth, double fc, double wp,
                 double exponent) {
    /*
        Soil water modifier, standard JULES/CABLE type approach

        equation 16 in Egea

        Note: we don't need to subtract the wp in the denominator here
              because our plant available water (paw) isn't bounded by
              the wilting point, it reaches zero

        Reference:
        ----------
        * Egea et al. (2011) Agricultural and Forest Meteorology.
    */

    double beta, theta;

    theta = paw / depth;
    beta = pow(theta / (fc - wp), exponent);
    if (beta > fc) {
        beta = 1.0;
    } else if (beta <= wp) {
        beta = 0.0;
    }

    return (beta);
}

double calc_sw_modifier(double theta, double c_theta, double n_theta) {
    /*
        Soil water modifier, equation 2 in Landsberg and Waring.
        Note: "The values of c_theta and n_theta are, nevertheless, chosen
              without specific empirical justification" :)

        Reference:
        ----------
        * Landsberg and Waring (1997) Forest Ecology and Management 95, 209-228.
    */
    return (1.0  / (1.0 + pow(((1.0 - theta) / c_theta), n_theta)));
}



void update_daily_water_struct(fluxes *f, double day_soil_evap,
                               double day_transpiration, double day_et,
                               double day_interception, double day_thoughfall,
                               double day_canopy_evap,
                               double day_runoff) {

    /* add half hour fluxes to day total store */
    f->soil_evap = day_soil_evap;
    f->transpiration = day_transpiration;
    f->et = day_et;
    f->interception = day_interception;
    f->throughfall = day_thoughfall;
    f->canopy_evap = day_canopy_evap;
    f->runoff = day_runoff;

    return;
}
