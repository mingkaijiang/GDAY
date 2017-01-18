/* ============================================================================
* Photosynthesis - C3 & C4
*
* see below
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   14.01.2016
*
* =========================================================================== */
#include "photosynthesis.h"

void simple_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s) {
    /* 
    Modifies mate_C3_photosynthesis using a simplier approach 
    
    */
    
    // fprintf(stderr, "npp in simple photo 1 %f\n", f->npp);
    double lue_avg, conv;
    
    /* Covert PAR units (umol PAR MJ-1) */
    conv = MJ_TO_J * J_2_UMOL;
    m->par *= conv;
    
    /* lue in umol C umol-1 PAR */
    lue_avg = lue_simplified(p, s, m->Ca);
    
    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
      f->apar = 0.0;
    else
      f->apar = m->par * s->fipar;
    
    /* convert umol m-2 d-1 -> gC m-2 d-1 */
    conv = UMOL_TO_MOL * MOL_C_TO_GRAMS_C;
    
    if (s->lai > 0.0) {
      /* calculation for npp */
      f->gpp_gCm2 = lue_avg * f->apar * conv; // / (1.0 - exp(-p->kext * s->lai));
    } else {
      f->gpp_gCm2 = 0.0;
    }
    
    // fprintf(stderr, "npp_gCm2 = %f\n", f->npp_gCm2);
    // fprintf(stderr, "apar = %f\n", f->apar);
    // fprintf(stderr, "par = %f\n", m->par);
    // fprintf(stderr, "fipar = %f\n", s->fipar);
    // fprintf(stderr, "lai = %f\n", s->lai);
    
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    
    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * G_AS_TONNES / M2_AS_HA;
    f->npp = f->npp_gCm2 * G_AS_TONNES / M2_AS_HA;
    
    /* save apar in MJ m-2 d-1 */
    f->apar *= UMOL_2_JOL * J_TO_MJ;
    
    return;
  
}


void mate_C3_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s,
                            double daylen, double ncontent, double pcontent) {
    /*

    MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    pm[noon to sunset]).

    MATE is connected to G'DAY via LAI and leaf N, P content.

    Plant autotrophic respiration is calculated via carbon-use efficiency
    (CUE=NPP/GPP).

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

    Rubisco kinetic parameter values are from:
    * Bernacchi et al. (2001) PCE, 24, 253-259.
    * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */
    double N0, P0, gamma_star_am,
           gamma_star_pm, Km_am, Km_pm, jmax_am, jmax_pm, vcmax_am, vcmax_pm,
           ci_am, ci_pm, alpha_am, alpha_pm, ac_am, ac_pm, aj_am, aj_pm,
           ap_am, ap_pm,
           asat_am, asat_pm, lue_am, lue_pm, lue_avg, conv;
    double mt = p->measurement_temp + DEG_TO_KELVIN;

    /* Calculate mate params & account for temperature dependencies */
    N0 = calculate_top_of_canopy_n(p, s, ncontent);   //Unit: g N m-2;

    if (c->pcycle == TRUE) {
        P0 = calculate_top_of_canopy_p(p, s, pcontent);   //Unit: g P m-2
    } else {
        P0 = 0.0;
    }
    
    gamma_star_am = calculate_co2_compensation_point(p, m->Tk_am, mt);
    gamma_star_pm = calculate_co2_compensation_point(p, m->Tk_pm, mt);

    Km_am = calculate_michaelis_menten_parameter(p, m->Tk_am, mt);
    Km_pm = calculate_michaelis_menten_parameter(p, m->Tk_pm, mt);

    if (c->pcycle == TRUE) {
        calculate_jmax_and_vcmax_with_p(c, p, s, m->Tk_am, N0, P0, &jmax_am,
                                        &vcmax_am, mt);
        calculate_jmax_and_vcmax_with_p(c, p, s, m->Tk_pm, N0, P0, &jmax_pm,
                                        &vcmax_pm, mt);
    } else {
        calculate_jmax_and_vcmax(c, p, s, m->Tk_am, N0, &jmax_am,
                                 &vcmax_am, mt);
        calculate_jmax_and_vcmax(c, p, s, m->Tk_pm, N0, &jmax_pm,
                                 &vcmax_pm, mt);
    }

    ci_am = calculate_ci(c, p, s, m->vpd_am, m->Ca);
    ci_pm = calculate_ci(c, p, s, m->vpd_pm, m->Ca);

    /* quantum efficiency calculated for C3 plants */
    alpha_am = calculate_quantum_efficiency(p, ci_am, gamma_star_am);
    alpha_pm = calculate_quantum_efficiency(p, ci_pm, gamma_star_pm);

    /* Rubisco carboxylation limited rate of photosynthesis */
    ac_am = assim(ci_am, gamma_star_am, vcmax_am, Km_am);
    ac_pm = assim(ci_pm, gamma_star_pm, vcmax_pm, Km_pm);

    /* Light-limited rate of photosynthesis allowed by RuBP regeneration */
    aj_am = assim(ci_am, gamma_star_am, jmax_am/4.0, 2.0*gamma_star_am);
    aj_pm = assim(ci_pm, gamma_star_pm, jmax_pm/4.0, 2.0*gamma_star_pm);

    if (c->pcycle == TRUE) {
        if (c->triose_p == TRUE) {
            /* Triose-phosphates limited rate of photosynthesis */
            ap_am = assim_p(P0);
            ap_pm = assim_p(P0);

            /* light-saturated photosynthesis rate at the top of the canopy */
            asat_am = MIN(aj_am, MIN(ac_am, ap_am));
            asat_pm = MIN(aj_pm, MIN(ac_pm, ap_pm));
        } else {
            asat_am = MIN(aj_am, ac_am);
            asat_pm = MIN(aj_pm, ac_pm);
        }
    } else {
        asat_am = MIN(aj_am, ac_am);
        asat_pm = MIN(aj_pm, ac_pm);
    }

    /* Covert PAR units (umol PAR MJ-1) */
    conv = MJ_TO_J * J_2_UMOL;
    m->par *= conv;
    
    /* LUE (umol C umol-1 PAR) ; note conversion in epsilon */
    lue_am = epsilon(p, asat_am, m->par, alpha_am, daylen);
    lue_pm = epsilon(p, asat_pm, m->par, alpha_pm, daylen);

    /* use average to simulate canopy photosynthesis */
    lue_avg = (lue_am + lue_pm) / 2.0;
    
    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
        f->apar = 0.0;
    else
        f->apar = m->par * s->fipar;

    /* convert umol m-2 d-1 -> gC m-2 d-1 */
    conv = UMOL_TO_MOL * MOL_C_TO_GRAMS_C;
    
    f->gpp_gCm2 = f->apar * lue_avg * conv;
    f->gpp_am = (f->apar / 2.0) * lue_am * conv;
    f->gpp_pm = (f->apar / 2.0) * lue_pm * conv;
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    
    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * G_AS_TONNES / M2_AS_HA;
    f->npp = f->npp_gCm2 * G_AS_TONNES / M2_AS_HA;
    
    // fprintf(stderr, "npp_gCm2 = %f\n", f->npp_gCm2);
    // fprintf(stderr, "lue_avg = %f\n", lue_avg);
    // fprintf(stderr, "apar = %f\n", f->apar);
    // fprintf(stderr, "par = %f\n", m->par);
    // fprintf(stderr, "fipar = %f\n", s->fipar);
    // fprintf(stderr, "conv = %f\n", conv);
    // fprintf(stderr, "kn = %f\n", p->kn);
    // fprintf(stderr, "lai = %f\n", s->lai);
    // fprintf(stderr, "npp mate %f\n", f->npp);

    /* save apar in MJ m-2 d-1 */
    f->apar *= UMOL_2_JOL * J_TO_MJ;

    return;
}

double calculate_top_of_canopy_n(params *p, state *s, double ncontent)  {

    /*
    Calculate the canopy N at the top of the canopy (g N m-2), N0.
    Assuming an exponentially decreasing N distribution within the canopy:

    Note: swapped kext with kn;

    Returns:
    -------
    N0 : float (g N m-2)
        Top of the canopy N

    References:
    -----------
    * Chen et al 93, Oecologia, 93,63-69.

    */
    double N0;

    if (s->lai > 0.0) {
        /* calculation for canopy N content at the top of the canopy */
        N0 = ncontent * p->kn / (1.0 - exp(-p->kn * s->lai));
    } else {
        N0 = 0.0;
    }

    return (N0);
}

double calculate_top_of_canopy_p(params *p, state *s, double pcontent)  {

    /*
    Calculate the canopy P at the top of the canopy (g P m-2), P0.
    Assuming an exponentially decreasing P distribution within the canopy:

    Note: swapped kext with kp;

    Returns:
    -------
    P0 : float (g P m-2)
    Top of the canopy P

    */
    double P0;

    if (s->lai > 0.0) {
        /* calculation for canopy P content at the top of the canopy */
        P0 = pcontent * p->kp / (1.0 - exp(-p->kp * s->lai));
    } else {
        P0 = 0.0;
    }

    return (P0);
}

double calculate_co2_compensation_point(params *p, double Tk, double mt) {
    /*
        CO2 compensation point in the absence of mitochondrial respiration
        Rate of photosynthesis matches the rate of respiration and the net CO2
        assimilation is zero.

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)

        Returns:
        -------
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
    */
    return (arrh(mt, p->gamstar25, p->eag, Tk));
}

double arrh(double mt, double k25, double Ea, double Tk) {
    /*
        Temperature dependence of kinetic parameters is described by an
        Arrhenius function

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
    */
    return (k25 * exp((Ea * (Tk - mt)) / (mt * RGAS * Tk)));
}


double peaked_arrh(double mt, double k25, double Ea, double Tk, double deltaS,
                   double Hd) {
    /*
        Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.

    */
    double arg1, arg2, arg3;

    arg1 = arrh(mt, k25, Ea, Tk);
    arg2 = 1.0 + exp((mt * deltaS - Hd) / (mt * RGAS));
    arg3 = 1.0 + exp((Tk * deltaS - Hd) / (Tk * RGAS));


    return (arg1 * arg2 / arg3);
}

double calculate_michaelis_menten_parameter(params *p, double Tk, double mt) {
    /*
        Effective Michaelis-Menten coefficent of Rubisco activity

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)

        Returns:
        -------
        Km : float
            Effective Michaelis-Menten constant for Rubisco catalytic activity

        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Bernacchi et al. (2001) PCE, 24, 253-259.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */

    double Kc, Ko;

    /* Michaelis-Menten coefficents for carboxylation by Rubisco */
    Kc = arrh(mt, p->kc25, p->eac, Tk);

    /* Michaelis-Menten coefficents for oxygenation by Rubisco */
    Ko = arrh(mt, p->ko25, p->eao, Tk);

    /* return effective Michaelis-Menten coefficient for CO2 */
    return ( Kc * (1.0 + p->oi / Ko) ) ;

}

void calculate_jmax_and_vcmax(control *c, params *p, state *s, double Tk,
                              double N0, double *jmax, double *vcmax,
                              double mt) {
    /*
        Calculate the maximum RuBP regeneration rate for light-saturated
        leaves at the top of the canopy (Jmax) and the maximum rate of
        rubisco-mediated carboxylation at the top of the canopy (Vcmax).

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)
        N0 : float
            leaf N   (g N m-2)


        Returns:
        --------
        jmax : float (umol/m2/sec)
            the maximum rate of electron transport at 25 degC
        vcmax : float (umol/m2/sec)
            the maximum rate of electron transport at 25 degC
    */
    double jmax25, vcmax25;
    double conv;

    *vcmax = 0.0;
    *jmax = 0.0;

    if (c->modeljm == 0) {
        *jmax = p->jmax;
        *vcmax = p->vcmax;
    } else if (c->modeljm == 1) {
        /* the maximum rate of electron transport at 25 degC */
        jmax25 = p->jmaxna * N0 + p->jmaxnb;

        /* this response is well-behaved for TLEAF < 0.0 */
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                            p->delsj, p->edj);

        /* the maximum rate of electron transport at 25 degC */
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;

        *vcmax = arrh(mt, vcmax25, p->eav, Tk);

    } else if (c->modeljm == 2) {
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = arrh(mt, vcmax25, p->eav, Tk);

        jmax25 = p->jv_slope * vcmax25 - p->jv_intercept;
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk, p->delsj,
                               p->edj);

    } else if (c->modeljm == 3) {
        /* the maximum rate of electron transport at 25 degC */
        jmax25 = p->jmax;

        /* this response is well-behaved for TLEAF < 0.0 */
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                            p->delsj, p->edj);

        /* the maximum rate of electron transport at 25 degC */
        vcmax25 = p->vcmax;
        *vcmax = arrh(mt, vcmax25, p->eav, Tk);

    }

    /* reduce photosynthetic capacity with moisture stress */
    *jmax *= s->wtfac_root;
    *vcmax *= s->wtfac_root;
    /*  Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
    adj_for_low_temp(*(&jmax), Tk);
    adj_for_low_temp(*(&vcmax), Tk);

    return;

}

void calculate_jmax_and_vcmax_with_p(control *c, params *p, state *s, double Tk,
                              double N0, double P0, double *jmax, double *vcmax,
                              double mt) {
  /*
  Calculate the maximum RuBP regeneration rate for light-saturated
  leaves at the top of the canopy (Jmax) and the maximum rate of
  rubisco-mediated carboxylation at the top of the canopy (Vcmax).

  Parameters:
  ----------
  Tk : float
  air temperature (Kelvin)
  N0 : float
  leaf N   (g N m-2)
  P0 : float
  leaf P   (g P m-2)

  Returns:
  --------
  jmax : float (umol/m2/sec)
  the maximum rate of electron transport at 25 degC
  vcmax : float (umol/m2/sec)
  the maximum rate of electron transport at 25 degC
  */
  double jmax25, vcmax25;
  double jmax25p, jmax25n;
  double vcmax25p, vcmax25n;

  *vcmax = 0.0;
  *jmax = 0.0;

  if (c->modeljm == 0) {
    *jmax = p->jmax;
    *vcmax = p->vcmax;
  } else if (c->modeljm == 1) {
    /* the maximum rate of electron transport at 25 degC */
    jmax25n = p->jmaxna * pow(N0, p->jmaxnb);

    /* P limitation on jmax */
    jmax25p = p->jmaxpa * pow(P0, p->jmaxpb);

    jmax25 = MIN(jmax25n, jmax25p);

    /* this response is well-behaved for TLEAF < 0.0 */
    *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                        p->delsj, p->edj);

    /* the maximum rate of electron transport at 25 degC */
    vcmax25n = p->vcmaxna * pow(N0, p->vcmaxnb);

    /* P limitation on jmax */
    vcmax25p = p->vcmaxpa * pow(P0, p->vcmaxpb);

    vcmax25 = MIN(vcmax25n, vcmax25p);

    *vcmax = arrh(mt, vcmax25, p->eav, Tk);

  } else if (c->modeljm == 2) {
    vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
    *vcmax = arrh(mt, vcmax25, p->eav, Tk);

    jmax25 = p->jv_slope * vcmax25 - p->jv_intercept;
    *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk, p->delsj,
                        p->edj);

  } else if (c->modeljm == 3) {
    /* the maximum rate of electron transport at 25 degC */
    jmax25 = p->jmax;

    /* this response is well-behaved for TLEAF < 0.0 */
    *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                        p->delsj, p->edj);

    /* the maximum rate of electron transport at 25 degC */
    vcmax25 = p->vcmax;
    *vcmax = arrh(mt, vcmax25, p->eav, Tk);

  }

  /* reduce photosynthetic capacity with moisture stress */
  *jmax *= s->wtfac_root;
  *vcmax *= s->wtfac_root;
  /*  Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
  adj_for_low_temp(*(&jmax), Tk);
  adj_for_low_temp(*(&vcmax), Tk);
  
  // fprintf(stderr, "jmax25n %f\n", jmax25n);
  // fprintf(stderr, "jmax25p %f\n", jmax25p);
  // fprintf(stderr, "vcmax25n %f\n", vcmax25n);
  // fprintf(stderr, "vcmax25p %f\n", vcmax25p);

  return;

}


void adj_for_low_temp(double *param, double Tk) {
    /*
    Function allowing Jmax/Vcmax to be forced linearly to zero at low T

    Parameters:
    ----------
    Tk : float
        air temperature (Kelvin)
    */
    double lower_bound = 0.0;
    double upper_bound = 10.0;
    double Tc;

    Tc = Tk - DEG_TO_KELVIN;

    if (Tc < lower_bound)
        *param = 0.0;
    else if (Tc < upper_bound)
        *param *= (Tc - lower_bound) / (upper_bound - lower_bound);

    return;
}

double calculate_ci(control *c, params *p, state *s, double vpd, double Ca) {
    /*
    Calculate the intercellular (Ci) concentration

    Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
    A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.

    Parameters:
    ----------
    vpd : float
        vapour pressure deficit [Pa]
    Ca : float
        ambient co2 concentration

    Returns:
    -------
    ci:ca : float
        ratio of intercellular to atmospheric CO2 concentration

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    */

    double g1w, cica, ci=0.0;

    if (c->gs_model == MEDLYN) {
        g1w = p->g1 * s->wtfac_root;
        cica = g1w / (g1w + sqrt(vpd * PA_2_KPA));
        ci = cica * Ca;
    } else {
        prog_error("Only Belindas gs model is implemented", __LINE__);
    }

    return (ci);
}

double calculate_quantum_efficiency(params *p, double ci, double gamma_star) {
    /*

    Quantum efficiency for AM/PM periods replacing Sands 1996
    temperature dependancy function with eqn. from Medlyn, 2000 which is
    based on McMurtrie and Wang 1993.

    Parameters:
    ----------
    ci : float
        intercellular CO2 concentration.
    gamma_star : float [am/pm]
        CO2 compensation point in the abscence of mitochondrial respiration

    Returns:
    -------
    alpha : float
        Quantum efficiency

    References:
    -----------
    * Medlyn et al. (2000) Can. J. For. Res, 30, 873-888
    * McMurtrie and Wang (1993) PCE, 16, 1-13.

    */
    return (assim(ci, gamma_star, p->alpha_j/4.0, 2.0*gamma_star));
}

double assim(double ci, double gamma_star, double a1, double a2) {
    /*
    Morning and afternoon calcultion of photosynthesis with the
    limitation defined by the variables passed as a1 and a2, i.e. if we
    are calculating vcmax or jmax limited.

    Parameters:
    ----------
    ci : float
        intercellular CO2 concentration.
    gamma_star : float
        CO2 compensation point in the abscence of mitochondrial respiration
    a1 : float
        variable depends on whether the calculation is light or rubisco
        limited.
    a2 : float
        variable depends on whether the calculation is light or rubisco
        limited.

    Returns:
    -------
    assimilation_rate : float
        assimilation rate assuming either light or rubisco limitation.
    */
    if (ci < gamma_star)
        return (0.0);
    else
        return (a1 * (ci - gamma_star) / (a2 + ci));

}

double assim_p(double P0) {
    //
    // Calculate photosynthesis assimilation rate based on P limitation
    // Ref: Ellsworth et al., 2015, Plant, Cell and Environment, able 2
    //
    // Returns:
    // -------
    // ap: float
    // assimilation rate assuming P limitation
    //
    double ap, conv, tp;

    conv = MMOL_2_MOL * 31.0;
    tp = 6.51 + 14.64 * (P0 * conv);
    ap = 3.0 * tp;

    return(ap);
}

double lue_simplified(params *p, state *s, double co2) {
    /*
     * New LUE function replacing epsilon function for a simplified calculation
     * of LUE
     * 
     * Calculations:
     *   CaResp = 1.632 * (co2-60.9) / (co2+121.8)    ##RCO2
     *   Nresp = min(df/Nref, 1)                      ##Rate-limiting effect of low N
     *   return(LUE0 * CaResp * Nresp)
     * 
     * Parameters:
     * Nref: leaf N:C for saturation of photosynthesis   
     * LUE0: maximum LUE in kg C GJ-1
     */
    double lue, CaResp, Nresp, conv, nref;
  
    CaResp = 1.632 * (co2 - 60.9) / (co2 + 121.8);
    Nresp = MIN(s->shootnc / p->nref, 1);
    
    /* converting unit for lue0 from kg C GJ-1 to umol C umol -1 PAR */
    conv = (KG_AS_G / MOL_C_TO_GRAMS_C * MOL_TO_UMOL) / (J_2_UMOL * GJ_TO_J);
      
    lue = p->lue0 * conv * CaResp * Nresp;
    
    // fprintf(stderr, "co2 %f\n", co2);
    // fprintf(stderr, "CaResp %f\n", CaResp);
    // fprintf(stderr, "Nresp %f\n", Nresp);
    // fprintf(stderr, "conv %f\n", conv);
    // fprintf(stderr, "shootnc %f\n", s->shootnc);
    // fprintf(stderr, "shootc %f\n", s->shoot);
    // fprintf(stderr, "shootn %f\n", s->shootn);
    // fprintf(stderr, "plantn %f\n", s->plantn);
    // fprintf(stderr, "soiln %f\n", s->soiln);
    
    return (lue);
}

double epsilon(params *p, double asat, double par, double alpha,
               double daylen) {
    /*
    Canopy scale LUE using method from Sands 1995, 1996.

    Sands derived daily canopy LUE from Asat by modelling the light response
    of photosysnthesis as a non-rectangular hyperbola with a curvature
    (theta) and a quantum efficiency (alpha).

    Assumptions of the approach are:
     - horizontally uniform canopy
     - PAR varies sinusoidally during daylight hours
     - extinction coefficient is constant all day
     - Asat and incident radiation decline through the canopy following
       Beer's Law.
     - leaf transmission is assumed to be zero.

    * Numerical integration of "g" is simplified to 6 intervals.

    Parameters:
    ----------
    asat : float
        Light-saturated photosynthetic rate at the top of the canopy
    par : float
        photosyntetically active radiation (umol m-2 d-1)
    theta : float
        curvature of photosynthetic light response curve
    alpha : float
        quantum yield of photosynthesis (mol mol-1)

    Returns:
    -------
    lue : float
        integrated light use efficiency over the canopy (umol C umol-1 PAR)

    Notes:
    ------
    NB. I've removed solar irradiance to PAR conversion. Sands had
    gamma = 2000000 to convert from SW radiation in MJ m-2 day-1 to
    umol PAR on the basis that 1 MJ m-2 = 2.08 mol m-2 & mol to umol = 1E6.
    We are passing PAR in umol m-2 d-1, thus avoiding the above.

    References:
    -----------
    See assumptions above...
    * Sands, P. J. (1995) Australian Journal of Plant Physiology,
      22, 601-14.

    */
    double delta, q, integral_g, sinx, arg1, arg2, arg3, lue, h;
    int i;

    /* subintervals scalar, i.e. 6 intervals */
    delta = 0.16666666667;

    /* number of seconds of daylight */
    h = daylen * SECS_IN_HOUR;

    if (asat > 0.0) {
        /* normalised daily irradiance */
        q = M_PI * p->kext * alpha * par / (2.0 * h * asat);
        integral_g = 0.0;
        for (i = 1; i < 13; i+=2) {
            sinx = sin(M_PI * i / 24.);
            arg1 = sinx;
            arg2 = 1.0 + q * sinx;
            arg3 = (sqrt(pow((1.0 + q * sinx), 2) - 4.0 * p->theta * q * sinx));
            integral_g += arg1 / (arg2 + arg3);
        }
        integral_g *= delta;
        lue = alpha * integral_g * M_PI;
    } else {
        lue = 0.0;
    }

    return (lue);
}

