/* ============================================================================
* Calls photosynthesis model, water balance and evolves aboveground plant
* C & Nstate. Pools recieve C through allocation of accumulated photosynthate
* and N from both soil uptake and retranslocation within the plant. Key feedback
* through soil N mineralisation and plant N uptake
*
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   17.02.2015
*
* =========================================================================== */
#include "plant_growth.h"

void calc_day_growth(control *c, fluxes *f, met_arrays *ma,
                     met *m, nrutil *nr, params *p, state *s, double day_length,
                     int doy, double fdecay, double rdecay)
{
   double dummy=0.0;
   double nitfac, pitfac, npitfac;
   double ncbnew, nccnew, ncwimm, ncwnew;
   double pcbnew, pccnew, pcwimm, pcwnew;
   int    recalc_wb;


    /* calculate daily GPP/NPP, respiration and update water balance */
    carbon_daily_production(c, f, m, p, s, day_length);

    // leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
    //foliage in young stand, and leaf P:C as a fraction of Pcmaxyoung
    nitfac = MIN(1.0, s->shootnc / p->ncmaxf);
    pitfac = MIN(1.0, s->shootpc / p->pcmaxf);
    
    // fprintf(stderr, "shootnc %f\n", s->shootnc);
    // fprintf(stderr, "shootpc %f\n", s->shootpc);
    // fprintf(stderr, "shootc %f\n", s->shoot);
    // fprintf(stderr, "shootn %f\n", s->shootn);
    // fprintf(stderr, "shootp %f\n", s->shootp);

    
    /* checking for pcycle control parameter */
    if (c->pcycle == TRUE) {
        npitfac = MIN(nitfac, pitfac);
    } else {
        npitfac = nitfac;
    }

    /* figure out the C allocation fractions */
    /* daily allocation ...*/
    calc_carbon_allocation_fracs(c, f, p, s, npitfac);

    /* Distribute new C, N and P through the system */
    carbon_allocation(c, f, p, s, npitfac, doy);

    calculate_cnp_wood_ratios(c, p, s, npitfac, nitfac, pitfac,
                              &ncbnew, &nccnew, &ncwimm,
                              &ncwnew, &pcbnew, &pccnew, &pcwimm,
                              &pcwnew);
                              
    recalc_wb = np_allocation(c, f, p, s, ncbnew, nccnew, ncwimm, ncwnew,
                              pcbnew, pccnew, pcwimm, pcwnew,
                              fdecay, rdecay, doy);
    
    update_plant_state(c, f, p, s, fdecay, rdecay, doy);

    precision_control(f, s);

    return;
}

void carbon_daily_production(control *c, fluxes *f, met *m, params *p, state *s,
                             double daylen) {
    /* Calculate GPP, NPP and plant respiration at the daily timestep

    Parameters:
    -----------
    daylen : float
        daytime length (hrs)

    References:
    -----------
    * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.
    */
    double fc, leafn, leafp, ncontent, pcontent;
  
  // fprintf(stderr, "npp in carbon_daily_production 1 %f\n", f->npp);
  
    if (s->lai > 0.0) {
        /* average leaf nitrogen content (g N m-2 leaf) */
        leafn = (s->shootnc * p->cfracts / p->sla * KG_AS_G);
        /* average leaf phosphorus content (g P m-2 leaf) */
        leafp = (s->shootpc * p->cfracts / p->sla * KG_AS_G);

        /* total nitrogen content of the canopy */
        ncontent = leafn * s->lai;
        /* total phosphorus content of the canopy */
        pcontent = leafp * s->lai;

    } else {
        ncontent = 0.0;
        pcontent = 0.0;
    }

    /* When canopy is not closed, canopy light interception is reduced
        - calculate the fractional ground cover */
    if (s->lai < p->lai_closed) {
        /* discontinuous canopies */
        fc = s->lai / p->lai_closed;
    } else {
        fc = 1.0;
    }

    /* fIPAR - the fraction of intercepted PAR = IPAR/PAR incident at the
       top of the canopy, accounting for partial closure based on Jackson
       and Palmer (1979). */
    if (s->lai > 0.0)
        s->fipar = ((1.0 - exp(-p->kext * s->lai / fc)) * fc);
    else
        s->fipar = 0.0;
    
    /* Estimate photosynthesis */
    simple_photosynthesis(c, f, m, p, s);

    /* Calculate plant respiration */

    f->auto_resp = f->gpp * p->cue;
    
    /* Calculate NPP */
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;
    
    return;
}

void calculate_cnp_wood_ratios(control *c, params *p, state *s,
                               double npitfac, double nitfac, double pitfac,
                               double *ncbnew, double *nccnew,
                               double *ncwimm, double *ncwnew,
                               double *pcbnew, double *pccnew,
                               double *pcwimm, double *pcwnew) {
    /* Estimate the N:C and P:C ratio in the branch and stem. Option to vary
    the N:C and P:C ratio of the stem following Jeffreys (1999) or keep it a fixed
    fraction

    Parameters:
    -----------
    npitfac: float
       min of nitfac and pitfac;
    nitfac: float
       leaf N:C as a fraction of Ncmaxyoung;
    pitfac : float
       leaf P:C as a fraction of Pcmaxyoung;

    Returns:
    --------
    ncbnew : float
        N:C ratio of branch
    nccnew : double
        N:C ratio of coarse root
    ncwimm : float
        N:C ratio of immobile stem
    ncwnew : float
        N:C ratio of mobile stem
    pcbnew : float
        P:C ratio of branch
    pccnew : double
        P:C ratio of coarse root
    pcwimm : float
        P:C ratio of immobile stem
    pcwnew : float
        P:C ratio of mobile stem

    References:
    ----------
    * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
      with modelled implications for forest productivity under elevated
      atmospheric carbon dioxide. PhD.
    */

    /* calculate N:C ratios */
    if (nitfac < npitfac) {
        /* n:c ratio of new branch wood*/
        *ncbnew = p->ncbnew + nitfac * (p->ncbnew - p->ncbnewz);

        /* n:c ratio of coarse root */
        *nccnew = p->nccnew + nitfac * (p->nccnew - p->nccnewz);

        /* fixed N:C in the stemwood */
        if (c->fixed_stem_nc) {
            /* n:c ratio of stemwood - immobile pool and new ring */
            *ncwimm = p->ncwimm + nitfac * (p->ncwimm - p->ncwimmz);

            /* New stem ring N:C at critical leaf N:C (mobile) */
            *ncwnew = p->ncwnew + nitfac * (p->ncwnew - p->ncwnewz);

        // vary stem N:C based on reln with foliage, see Jeffreys PhD thesis.
        // Jeffreys 1999 showed that N:C ratio of new wood increases with
        // foliar N:C ratio,modelled here based on evidence as a linear
        // function.
        } else {
            *ncwimm = MAX(0.0, (0.0282 * s->shootnc + 0.000234) * p->fhw);

            /* New stem ring N:C at critical leaf N:C (mobile) */
            *ncwnew = MAX(0.0, 0.162 * s->shootnc - 0.00143);
        }
    } else {
        /* n:c ratio of new branch wood*/
        *ncbnew = p->ncbnew + npitfac * (p->ncbnew - p->ncbnewz);

        /* n:c ratio of coarse root */
        *nccnew = p->nccnew + npitfac * (p->nccnew - p->nccnewz);

        /* fixed N:C in the stemwood */
        if (c->fixed_stem_nc) {
            /* n:c ratio of stemwood - immobile pool and new ring */
            *ncwimm = p->ncwimm + npitfac * (p->ncwimm - p->ncwimmz);

            /* New stem ring N:C at critical leaf N:C (mobile) */
            *ncwnew = p->ncwnew + npitfac * (p->ncwnew - p->ncwnewz);

        // vary stem N:C based on reln with foliage, see Jeffreys PhD thesis.
        // Jeffreys 1999 showed that N:C ratio of new wood increases with
        // foliar N:C ratio,modelled here based on evidence as a linear
        // function.
        } else {
            *ncwimm = MAX(0.0, (0.0282 * s->shootnc + 0.000234) * p->fhw);

            /* New stem ring N:C at critical leaf N:C (mobile) */
            *ncwnew = MAX(0.0, 0.162 * s->shootnc - 0.00143);
        }
    }

    /* calculate P:C ratios */
    if (pitfac < npitfac) {
        /* p:c ratio of new branch wood*/
        *pcbnew = p->pcbnew + pitfac * (p->pcbnew - p->pcbnewz);

        /* p:c ratio of coarse root */
        *pccnew = p->pccnew + pitfac * (p->pccnew - p->pccnewz);

        /* fixed P:C in the stemwood */
        if (c->fixed_stem_pc) {
            /* p:c ratio of stemwood - immobile pool and new ring */
            *pcwimm = p->pcwimm + pitfac * (p->pcwimm - p->pcwimmz);

            /* New stem ring P:C at critical leaf P:C (mobile) */
            *pcwnew = p->pcwnew + pitfac * (p->pcwnew - p->pcwnewz);

            /* vary stem P:C based on reln with foliage,
            equation based on data from Attiwill 1978 - 1980 paper series. */
        } else {
            *pcwimm = MAX(0.0, -0.0016 * s->shootpc + 0.000003);

            // New stem ring P:C at critical leaf P:C (mobile),
            // equation based on data from Attiwill 1978 - 1980 paper series
            *pcwnew = MAX(0.0, -0.0022 * s->shootpc + 0.000009);
        }
    } else {
        /* p:c ratio of new branch wood*/
        *pcbnew = p->pcbnew + npitfac * (p->pcbnew - p->pcbnewz);

        /* p:c ratio of coarse root */
        *pccnew = p->pccnew + npitfac * (p->pccnew - p->pccnewz);

        /* fixed P:C in the stemwood */
        if (c->fixed_stem_pc) {
            /* p:c ratio of stemwood - immobile pool and new ring */
            *pcwimm = p->pcwimm + npitfac * (p->pcwimm - p->pcwimmz);

            /* New stem ring P:C at critical leaf P:C (mobile) */
            *pcwnew = p->pcwnew + npitfac * (p->pcwnew - p->pcwnewz);

        /* vary stem P:C based on reln with foliage,
        equation based on data from Attiwill 1978 - 1980 paper series. */
        } else {
            *pcwimm = MAX(0.0, -0.0016 * s->shootpc + 0.000003);

            /* New stem ring P:C at critical leaf P:C (mobile),
            equation based on data from Attiwill 1978 - 1980 paper series */
            *pcwnew = MAX(0.0, -0.0022 * s->shootpc + 0.000009);
        }
    }

    return;
}

int np_allocation(control *c, fluxes *f, params *p, state *s, double ncbnew,
                  double nccnew, double ncwimm, double ncwnew, double pcbnew,
                  double pccnew, double pcwimm, double pcwnew,double fdecay,
                  double rdecay, int doy) {
    /*
        Nitrogen and phosphorus distribution - allocate available N and
        P (mineral) through system. N and P is first allocated to the woody
        component, surplus N and P is then allocated to the shoot and roots
        with flexible ratios.

        References:
        -----------
        McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.

        Parameters:
        -----------
        ncbnew : float
            N:C ratio of branch
        ncwimm : float
            N:C ratio of immobile stem
        ncwnew : float
            N:C ratio of mobile stem
        pcbnew : float
            P:C ratio of branch
        pcwimm : float
            P:C ratio of immobile stem
        pcwnew : float
            P:C ratio of mobile stem
        fdecay : float
            foliage decay rate
        rdecay : float
            fine root decay rate
    */

    int    recalc_wb;
    double nsupply, psupply, rtot, ntot, ptot, arg;
    double depth_guess = 1.0;

    /* default is we don't need to recalculate the water balance,
       however if we cut back on NPP due to available N and P below then we do
       need to do this */
    recalc_wb = FALSE;

    /* N and P retranslocated proportion from dying plant tissue and stored within
       the plant */
    f->retrans = nitrogen_retrans(c, f, p, s, fdecay, rdecay, doy);
    f->retransp = phosphorus_retrans(c, f, p, s, fdecay, rdecay, doy);
    f->nuptake = calculate_nuptake(c, p, s);
    f->puptake = calculate_puptake(c, p, s, f);
    
    // fprintf(stderr, "nuptake %f\n", f->nuptake);
    // fprintf(stderr, "puptake %f\n", f->puptake);

    /* Mineralised nitrogen lost from the system by volatilisation/leaching */
    f->nloss = p->rateloss * s->inorgn;

    /* Mineralised P lost from the system by leaching */
    f->ploss = p->prateloss * s->inorgavlp;

    
    // fprintf(stderr, "nloss %f\n", f->nloss);
    // fprintf(stderr, "ploss %f\n", f->ploss);
    

    /* total nitrogen/phosphorus to allocate */
    ntot = MAX(0.0, f->nuptake + f->retrans);
    ptot = MAX(0.0, f->puptake + f->retransp);
    
    /* allocate N to pools with fixed N:C ratios */
    /* N flux into new ring (immobile component -> structural components) */
    f->npstemimm = f->npp * f->alstem * ncwimm;
    
    /* N flux into new ring (mobile component -> can be retrans for new
    woody tissue) */
    f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
    f->npbranch = f->npp * f->albranch * ncbnew;

    /* allocate P to pools with fixed P:C ratios */
    f->ppstemimm = f->npp * f->alstem * pcwimm;
    f->ppstemmob = f->npp * f->alstem * (pcwnew - pcwimm);
    f->ppbranch = f->npp * f->albranch * pcbnew;
    
    /* If we have allocated more N than we have avail, cut back C prodn */
    arg = f->npstemimm + f->npstemmob + f->npbranch;
    if (arg > ntot && c->fixleafnc == FALSE && c->ncycle) {
      recalc_wb = cut_back_production(c, f, p, s, ntot, ncbnew, nccnew,
                                      ncwimm, ncwnew, doy);
    }
    
    /* If we have allocated more P than we have avail, cut back C prodn */
    arg = f->ppstemimm + f->ppstemmob + f->ppbranch;
    if (arg > ptot && c->fixleafpc == FALSE && c->pcycle) {
      recalc_wb = cut_back_production(c, f, p, s, ptot, pcbnew, pccnew,
                                      pcwimm, pcwnew, doy);
    }
    
    /* Nitrogen reallocation to flexible-ratio pools */
    ntot -= f->npbranch + f->npstemimm + f->npstemmob;
    ntot = MAX(0.0, ntot);
    
    /* allocate remaining N to flexible-ratio pools */
    f->npleaf = ntot * f->alleaf / (f->alleaf + f->alroot * p->ncrfac);
    f->nproot = ntot - f->npleaf;
    
    /* Phosphorus reallocation to flexible-ratio pools */
    ptot -= f->ppbranch + f->ppstemimm + f->ppstemmob;
    ptot = MAX(0.0, ptot);
    
    /* allocate remaining P to flexible-ratio pools */
    f->ppleaf = ptot * f->alleaf / (f->alleaf + f->alroot * p->pcrfac);
    f->pproot = ptot - f->ppleaf;
    

    return (recalc_wb);
}




int cut_back_production(control *c, fluxes *f, params *p, state *s,
                        double tot, double xcbnew, double xccnew,
                        double xcwimm, double xcwnew, int doy) {

    double lai_inc, conv;
    double pcbnew, pccnew, pcwimm, pcwnew;
    /* default is we don't need to recalculate the water balance,
       however if we cut back on NPP due to available N and P below then we do
       need to do this */
    int recalc_wb = FALSE;

    /* Need to readjust the LAI for the reduced growth as this will
       have already been increased. First we need to figure out how
       much we have increased LAI by, important it is done here
       before cpleaf is reduced! */
    if (float_eq(s->shoot, 0.0)) {
        lai_inc = 0.0;
    } else {
        lai_inc = (f->cpleaf *
                   (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
                   f->deadleaves * s->lai / s->shoot);
    }

    f->npp *= tot / (f->npstemimm + f->npstemmob + \
                      f->npbranch);

    /* need to adjust growth values accordingly as well */
    f->cpleaf = f->npp * f->alleaf;
    f->cproot = f->npp * f->alroot;
    f->cpbranch = f->npp * f->albranch;
    f->cpstem = f->npp * f->alstem;



    if (c->pcycle) {
        f->ppbranch = f->npp * f->albranch * xcbnew;
        f->ppstemimm = f->npp * f->alstem * xcwimm;
        f->ppstemmob = f->npp * f->alstem * (xcwnew - xcwimm);
    } else {
        f->npbranch = f->npp * f->albranch * xcbnew;
        f->npstemimm = f->npp * f->alstem * xcwimm;
        f->npstemmob = f->npp * f->alstem * (xcwnew - xcwimm);
    }

    /* Also need to recalculate GPP and thus Ra and return a flag
       so that we know to recalculate the water balance. */
    f->gpp = f->npp / p->cue;
    conv = G_AS_TONNES / M2_AS_HA;
    f->gpp_gCm2 = f->gpp / conv;
    
    // following two lines commented out for simplified version */
    f->gpp_am = f->gpp_gCm2 / 2.0;
    f->gpp_pm = f->gpp_gCm2 / 2.0;

    /* New respiration flux */
    f->auto_resp =  f->gpp - f->npp;
    recalc_wb = TRUE;

    /* Now reduce LAI for down-regulated growth. */
    /* update leaf area [m2 m-2] */
    if (float_eq(s->shoot, 0.0)) {
      s->lai = 0.0;
    } else {
      s->lai -= lai_inc;
      s->lai += (f->cpleaf *
        (p->sla * M2_AS_HA / \
        (KG_AS_TONNES * p->cfracts)) -
        f->deadleaves * s->lai / s->shoot);
    }

    return (recalc_wb);
}

double calculate_growth_stress_limitation(params *p, state *s, control *c) {
    //
    // Calculate level of stress due to nitrogen, phosphorus or water
    // availability
    //
    double nlim, plim, current_limitation, nutrient_lim;
    double nc_opt = 0.04;
    double pc_opt = 0.004;

    /* N limitation based on leaf NC ratio */
    if (s->shootnc < p->nf_min) {
        nlim = 0.0;
    } else if (s->shootnc < nc_opt && s->shootnc > p->nf_min) {
        nlim = 1.0 - ((nc_opt - s->shootnc) / (nc_opt - p->nf_min));
    } else {
        nlim = 1.0;
    }

    /*
     * Limitation by nutrients or water. Water constraint is
     * implicit, in that, water stress results in an increase of root mass,
     * which are assumed to spread horizontally within the rooting zone.
     * So in effect, building additional root mass doesnt alleviate the
     * water limitation within the model. However, it does more
     * accurately reflect an increase in root C production at a water
     * limited site. This implementation is also consistent with other
     * approaches, e.g. LPJ. In fact I dont see much evidence for models
     * that have a flexible bucket depth. Minimum constraint is limited to
     * 0.1, following Zaehle et al. 2010 (supp), eqn 18.
     */
    current_limitation = MAX(0.1, nlim);

    if(c->pcycle == TRUE) {
        /* P limitation based on leaf PC ratio */
        if (s->shootpc < p->pf_min) {
            plim = 0.0;
        } else if (s->shootpc < pc_opt && s->shootpc > p->pf_min) {
            plim = 1.0 - ((pc_opt - s->shootpc) / (pc_opt - p->pf_min));
        } else {
            plim = 1.0;
        }
        nutrient_lim = MIN(nlim, plim);
        current_limitation = MAX(0.1,nutrient_lim);
    }
    
    /*fprintf(stderr, "plim %f\n", plim);
    fprintf(stderr, "nlim %f\n", nlim);
    fprintf(stderr, "shootnc %f\n", s->shootnc);
    fprintf(stderr, "shootpc %f\n", s->shootpc);
    fprintf(stderr, "shootc %f\n", s->shoot);
    fprintf(stderr, "shootn %f\n", s->shootn);
    fprintf(stderr, "shootp %f\n", s->shootp);
    */
    return (current_limitation);
}


void calc_carbon_allocation_fracs(control *c, fluxes *f, params *p, state *s,
                                  double npitfac) {
    /* Carbon allocation fractions to move photosynthate through the plant.

    Parameters:
    -----------
    npitfac : float
        the smallest value of leaf N:C as a fraction of 'Ncmaxf' (max 1.0) &
        leaf P:C as a fraction of "Pcmaxf" (max 1.0)

    Returns:
    --------
    alleaf : float
        allocation fraction for shoot
    alroot : float
        allocation fraction for fine roots
    albranch : float
        allocation fraction for branches
    alstem : float
        allocation fraction for stem

    References:
    -----------
    Corbeels, M. et al (2005) Ecological Modelling, 187, 449-474.
    McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
    */
    double total_alloc;

    f->alleaf = p->c_alloc_f;
    
    f->albranch = p->c_alloc_b;

    f->alroot = p->c_alloc_r;
    
    f->alstem = p->c_alloc_s;
    
    /* Total allocation should be one, if not print warning */
    total_alloc = f->alroot + f->alleaf + f->albranch + f->alstem;
    if (total_alloc > 1.0+EPSILON) {
        fprintf(stderr, "Allocation fracs > 1: %.13f\n", total_alloc);
        exit(EXIT_FAILURE);
    }

    return;
}

void carbon_allocation(control *c, fluxes *f, params *p, state *s,
                       double npitfac, int doy) {
    /* C distribution - allocate available C through system

    Parameters:
    -----------
    npitfac : float
        leaf N:C as a fraction of 'Ncmaxf' (max 1.0)
    */
    double days_left;
    f->cpleaf = f->npp * f->alleaf;
    f->cproot = f->npp * f->alroot;
    f->cpbranch = f->npp * f->albranch;
    f->cpstem = f->npp * f->alstem;

    /* evaluate SLA of new foliage accounting for variation in SLA
       with tree and leaf age (Sands and Landsberg, 2002). Assume
       SLA of new foliage is linearly related to leaf N:C ratio
       via nitfac. Based on date from two E.globulus stands in SW Aus, see
       Corbeels et al (2005) Ecological Modelling, 187, 449-474.
       (m2 onesided/kg DW)
     This needs to be updated to consider P effect
    */
    p->sla = p->slazero + npitfac * (p->slamax - p->slazero);

    /* update leaf area [m2 m-2] */
    if (float_eq(s->shoot, 0.0)) {
      s->lai = 0.0;
    } else {
      s->lai += (f->cpleaf *
        (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
        f->deadleaves * s->lai / s->shoot);
    }

    return;
}

void update_plant_state(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
    /*
    Daily change in C content

    Parameters:
    -----------
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate

    */
    double ncmaxf, ncmaxr;
    double extrasn, extrarn;   /* extra_s_n and extra_r_n - extra shoot/root n uptake */
    double pcmaxf, pcmaxr;
    double extrasp, extrarp;   /* extra_s_p and extra_r_p - extra shoot/root p uptake */

    /*
    ** Carbon pools
    */
    s->shoot += f->cpleaf - f->deadleaves;
    s->root += f->cproot - f->deadroots;
    s->branch += f->cpbranch - f->deadbranch;
    s->stem += f->cpstem - f->deadstems;

    //fprintf(stderr, "cproot %f\n", f->cproot);
    //fprintf(stderr, "deadroots %f\n", f->deadroots);

    /* annoying but can't see an easier way with the code as it is.
       If we are modelling grases, i.e. no stem them without this
       the sapwood will end up being reduced to a silly number as
       deadsapwood will keep being removed from the pool, even though there
       is no wood. */
    if (float_eq(s->stem, 0.01)) {
        s->sapwood = 0.01;
    } else if (s->stem < 0.01) {
        s->sapwood = 0.01;
    } else {
        s->sapwood += f->cpstem - f->deadsapwood;
    }


    /*
    ** Nitrogen and Phosphorus pools
    */
    s->shootn += f->npleaf - fdecay * s->shootn;
    s->shootp += f->ppleaf - fdecay * s->shootp;

    s->branchn += f->npbranch - p->bdecay * s->branchn;
    s->rootn += f->nproot - rdecay * s->rootn;
    s->stemnimm += f->npstemimm - p->wdecay * s->stemnimm;
    s->stemnmob += (f->npstemmob - p->wdecay * s->stemnmob - p->retransmob *
                    s->stemnmob);
    s->stemn = s->stemnimm + s->stemnmob;

    s->branchp += f->ppbranch - p->bdecay * s->branchp;

    s->rootp += f->pproot - rdecay * s->rootp;

    s->stempimm += f->ppstemimm - p->wdecay * s->stempimm;

    s->stempmob += (f->ppstemmob - p->wdecay * s->stempmob - p->retransmob *
                    s->stempmob);

    s->stemp = s->stempimm + s->stempmob;

    /*
     =============================
     Enforce maximum N:C and P:C ratios.
     =============================
     */
    
    /* If foliage or root N/C exceeds its max, then N uptake is cut back
    Similarly, of foliage or root P/C exceeds max, then P uptake is cut back */
    
    /* maximum leaf n:c and p:c ratios is function of stand age*/
    ncmaxf = p->ncmaxf;
    pcmaxf = p->pcmaxf;
    
    extrasn = 0.0;
    
    if (s->lai > 0.0) {
      
      if (s->shootn > (s->shoot * ncmaxf)) {
        extrasn = s->shootn - s->shoot * ncmaxf;
        
        /* Ensure N uptake cannot be reduced below zero. */
        if (extrasn >  f->nuptake)
          extrasn = f->nuptake;
        
        s->shootn -= extrasn;
        //f->nuptake -= extrasn;
      }
    }
    
    extrasp = 0.0;
    if (s->lai > 0.0) {
      
      if (s->shootp > (s->shoot * pcmaxf)) {
        extrasp = s->shootp - s->shoot * pcmaxf;
        
        /* Ensure P uptake cannot be reduced below zero. */
        if (extrasp >  f->puptake)
          extrasp = f->puptake;
        
        s->shootp -= extrasp;
        //f->puptake -= extrasp;
      }
    }
    
    /* if root N:C ratio exceeds its max, then nitrogen uptake is cut
    back. n.b. new ring n/c max is already set because it is related
    to leaf n:c */
    
    /* max root n:c */
    ncmaxr = ncmaxf * p->ncrfac;
    extrarn = 0.0;
    if (s->rootn > (s->root * ncmaxr)) {
      extrarn = s->rootn - s->root * ncmaxr;
      
      /* Ensure N uptake cannot be reduced below zero. */
      if ((extrasn + extrarn) > f->nuptake)
        extrarn = f->nuptake - extrasn;
      
      s->rootn -= extrarn;
      f->nuptake -= (extrarn+extrasn);
    }
    
    /* max root p:c */
    pcmaxr = pcmaxf * p->pcrfac;
    extrarp = 0.0;
    if (s->rootp > (s->root * pcmaxr)) {
      extrarp = s->rootp - s->root * pcmaxr;
      
      /* Ensure P uptake cannot be reduced below zero. */
      if ((extrasp + extrarp) > f->puptake)
        extrarp = f->puptake - extrasp;
      
      s->rootp -= extrarp;
      f->puptake -= (extrarp + extrasp);
    }
    

    return;
}

void precision_control(fluxes *f, state *s) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-10;

    /* C, N & P state variables */
    if (s->shoot < tolerance) {
        f->deadleaves += s->shoot;
        f->deadleafn += s->shootn;
        f->deadleafp += s->shootp;
        s->shoot = 0.0;
        s->shootn = 0.0;
        s->shootp = 0.0;
    }

    if (s->branch < tolerance) {
        f->deadbranch += s->branch;
        f->deadbranchn += s->branchn;
        f->deadbranchp += s->branchp;
        s->branch = 0.0;
        s->branchn = 0.0;
        s->branchp = 0.0;
    }

    if (s->root < tolerance) {
        f->deadrootn += s->rootn;
        f->deadrootp += s->rootp;
        f->deadroots += s->root;
        s->root = 0.0;
        s->rootn = 0.0;
        s->rootp = 0.0;
    }

    /* Not setting these to zero as this just leads to errors with desert
       regrowth...instead seeding them to a small value with a CN~25 and CP~300. */

    if (s->stem < tolerance) {
        f->deadstems += s->stem;
        f->deadstemn += s->stemn;
        f->deadstemp += s->stemp;
        s->stem = 0.001;
        s->stemn = 0.00004;
        s->stemp = 0.000003;
        s->stemnimm = 0.00004;
        s->stemnmob = 0.0;
        s->stempimm = 0.000003;
        s->stempmob = 0.0;
    }

    /* need separate one as this will become very small if there is no
       mobile stem N/P */
    if (s->stemnmob < tolerance) {
        f->deadstemn += s->stemnmob;
        s->stemnmob = 0.0;
    }

    if (s->stemnimm < tolerance) {
        f->deadstemn += s->stemnimm;
        s->stemnimm = 0.00004;
    }

    if (s->stempmob < tolerance) {
      f->deadstemp += s->stempmob;
      s->stempmob = 0.0;
    }

    if (s->stempimm < tolerance) {
      f->deadstemp += s->stempimm;
      s->stempimm = 0.000003;

    }

    return;
}


double nitrogen_retrans(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
    /* Nitrogen retranslocated from senesced plant matter.
    Constant rate of n translocated from mobile pool

    Parameters:
    -----------
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate

    Returns:
    --------
    N retrans : float
        N retranslocated plant matter

    */
    double leafretransn, rootretransn, branchretransn,stemretransn;

    leafretransn = p->fretrans * fdecay * s->shootn;

    rootretransn = p->rretrans * rdecay * s->rootn;
    branchretransn = p->bretrans * p->bdecay * s->branchn;
    stemretransn = (p->wretrans * p->wdecay * s->stemnmob + p->retransmob *
                    s->stemnmob);

    /* store for NCEAS output */
    f->leafretransn = leafretransn;

    return (leafretransn + rootretransn + branchretransn +
            stemretransn);
}

double phosphorus_retrans(control *c, fluxes *f, params *p, state *s,
                          double fdecay, double rdecay, int doy) {
    /*
        Phosphorus retranslocated from senesced plant matter.
        Constant rate of p translocated from mobile pool

        Parameters:
        -----------
        fdecay : float
        foliage decay rate
        rdecay : float
        fine root decay rate

        Returns:
        --------
        P retrans : float
        P retranslocated plant matter
    */
    double leafretransp, rootretransp, branchretransp, stemretransp;

    leafretransp = p->fretransp * fdecay * s->shootp;
    

    rootretransp = p->rretrans * rdecay * s->rootp;
    branchretransp = p->bretrans * p->bdecay * s->branchp;
    stemretransp = (p->wretrans * p->wdecay * s->stempmob + p->retransmob *
                    s->stempmob);

    /* store for NCEAS output */
    f->leafretransp = leafretransp;

    return (leafretransp + rootretransp + branchretransp +
            stemretransp);
}

double calculate_nuptake(control *c, params *p, state *s) {
    /*
        N uptake depends on the rate at which soil mineral N is made
        available to the plants.

        Returns:
        --------
        nuptake : float
            N uptake

        References:
        -----------
        * Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.
        * Raich et al. 1991, Ecological Applications, 1, 399-429.

    */
    double nuptake, U0, Kr;

    if (c->nuptake_model == 0) {
        /* Constant N uptake */
        nuptake = p->nuptakez;
    } else if (c->nuptake_model == 1) {
        /* evaluate nuptake : proportional to dynamic inorganic N pool */
        nuptake = p->rateuptake * s->inorgn;
    } else if (c->nuptake_model == 2) {
        /* N uptake is a saturating function on root biomass following
           Dewar and McMurtrie, 1996. */

        /* supply rate of available mineral N */
        U0 = p->rateuptake * s->inorgn;
        Kr = p->kr;
        nuptake = MAX(U0 * s->root / (s->root + Kr), 0.0);
        
        // fprintf(stderr, "U0 %f\n", U0);
        // fprintf(stderr, "root %f\n", s->root);
        // fprintf(stderr, "root+Kr %f\n", s->root+Kr);

    } else {
        fprintf(stderr, "Unknown N uptake option\n");
        exit(EXIT_FAILURE);
    }

    return (nuptake);
}

double calculate_puptake(control *c, params *p, state *s, fluxes *f) {
    /*
        P uptake depends on the rate at which soil mineral P is made
        available to the plants.

        Returns:
        --------
        puptake : float
        P uptake
    */
    double puptake, U0, Kr;

    if (c->puptake_model == 0) {
        /* Constant P uptake */
        puptake = p->puptakez;
    } else if (c->puptake_model == 1) {
        // evaluate puptake : proportional to lab P pool that is
        // available to plant uptake
        puptake = p->prateuptake * s->inorgavlp;
    } else if (c->puptake_model == 2) {
        /* P uptake is a saturating function on root biomass, as N */

        /* supply rate of available mineral P */
        U0 = p->prateuptake * s->inorgavlp;


        Kr = p->krp;
        puptake = MAX(U0 * s->root / (s->root + Kr), 0.0);
    } else {
        fprintf(stderr, "Unknown P uptake option\n");
        exit(EXIT_FAILURE);
    }

    return (puptake);
}

