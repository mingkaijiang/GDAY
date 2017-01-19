#include "water_balance.h"
#include "nrutil.h"

void initialise_soils_day(control *c, fluxes *f, params *p, state *s) {
    /* Initialise soil water state & parameters  */

    double *fsoil_top = NULL, *fsoil_root = NULL;
    int     i;

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
         s->wtfac_topsoil = exp(b * s->predawn_swp);
        s->wtfac_root = exp(b * s->predawn_swp);
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
