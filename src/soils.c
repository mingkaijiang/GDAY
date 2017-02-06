/* ============================================================================
* Soil C, N and P flows into 4 litter pools (structural and metabolic, both
* above and belowground) and 3 SOM pools (Active, slow and passive). Soil P
* flows into 5 inorganic pools (parent, lab, sorb, ssorb, and occluded).
* In essence the CENTURY model.
*
* Active pool -> soil microbes & microbial products, turnover time of mths-yrs.
* Slow pool -> resistant plant material, turnover time of 20-50 yrs.
* Passive pool -> very resistant to decomp, turnover time of > 400 yrs.
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   14.08.2016
*
* =========================================================================== */
#include "soils.h"

void calculate_csoil_flows(control *c, fluxes *f, params *p, state *s,
                           double tsoil) {
    double lnleaf, lnroot;
    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    
    f->tfac_soil_decomp = calc_soil_temp_factor(tsoil);

    /* calculate model decay rates */
    calculate_decay_rates(c, f, p, s);

    /*
     * plant litter inputs to the metabolic and structural pools determined
     * by ratio of lignin/N ratio
     */
    lnleaf = calc_ligin_nratio_leaves(c, f, p);
    lnroot = calc_ligin_nratio_fine_roots(c, f, p);
    
    p->fmleaf = metafract(lnleaf);
    p->fmroot = metafract(lnroot);
    
    /* input from faeces */
    partition_plant_litter(c, f, p);
    cfluxes_from_structural_pool(f, p, s);
    cfluxes_from_metabolic_pool(f, p, s);
    cfluxes_from_active_pool(f, p, s, frac_microb_resp);
    cfluxes_from_slow_pool(f, p, s);
    cfluxes_from_passive_pool(f, p, s);
    calculate_soil_respiration(c, f, p, s);
    
    /* update the C pools */
    calculate_cpools(c, f, s, p);

    /* calculate NEP */
    f->nep = f->npp - f->hetero_resp;

    /* save fluxes for NCEAS output */
    f->co2_rel_from_surf_struct_litter = f->co2_to_air[0];
    f->co2_rel_from_soil_struct_litter = f->co2_to_air[1];
    f->co2_rel_from_surf_metab_litter = f->co2_to_air[2];
    f->co2_rel_from_soil_metab_litter = f->co2_to_air[3];
    f->co2_rel_from_active_pool = f->co2_to_air[4];
    f->co2_rel_from_slow_pool = f->co2_to_air[5];
    f->co2_rel_from_passive_pool = f->co2_to_air[6];

    return;
}

void calculate_decay_rates(control *c, fluxes *f, params *p, state *s) {
    /* Model decay rates - decomposition rates have a strong temperature
    and moisture dependency. Note same temperature is assumed for all 3
    SOM pools, found by Knorr et al (2005) to be untrue. N mineralisation
    depends on top soil moisture (most variable) (Connell et al. 1995)

    References:
    -----------
    Knorr et al. (2005) Nature, 433, 298-301.
    Connell et al. (1995) Biol. Fert. Soils, 20, 213-220.
    */
    double soil_text, lignin_cont_leaf, lignin_cont_root;

    /* abiotic decomposition factor - impact of soil moisture
       and soil temperature on microbial activity */

    /*  Effect of soil texture (silt + clay content) on active SOM turnover
        -> higher turnover for sandy soils */
    soil_text = 1.0 - (0.75 * p->finesoil);

    /* Impact of lignin content */
    lignin_cont_leaf = exp(-3.0 * p->ligshoot);
    lignin_cont_root = exp(-3.0 * p->ligroot);

    /* decay rate of surface structural pool */
    p->decayrate[0] = p->kdec1 * lignin_cont_leaf * f->tfac_soil_decomp;

    /* decay rate of surface metabolic pool */
    p->decayrate[1] = p->kdec2 * f->tfac_soil_decomp;

    /* decay rate of soil structural pool */
    p->decayrate[2] = p->kdec3 * lignin_cont_root * f->tfac_soil_decomp;

    /* decay rate of soil metabolic pool */
    p->decayrate[3] = p->kdec4 * f->tfac_soil_decomp;

    /* decay rate of active pool */
    p->decayrate[4] = p->kdec5 * soil_text * f->tfac_soil_decomp;

    /* decay rate of slow pool */
    p->decayrate[5] = p->kdec6 * f->tfac_soil_decomp;

    /* decay rate of passive pool */
    p->decayrate[6] = p->kdec7 * f->tfac_soil_decomp;

    return;
}

double calc_soil_temp_factor(double tsoil) {
    /*
    Soil-temperature activity factor (A9). Fit to Parton's fig 2a

    Parameters:
    -----------
    tsoil : double
        soil temperature (deg C)

    Returns:
    --------
    tfac : double
        soil temperature factor [degC]

    */
    double tfac;
    if (tsoil > 0.0) {
        tfac = MAX(0.0, 0.0326 + 0.00351 * pow(tsoil, 1.652) - \
                        pow((tsoil / 41.748), 7.19));
    } else {
        /* negative number cannot be raised to a fractional power
           number would need to be complex */
        tfac = 0.0;
    }

    return (tfac);
}


double calc_ligin_nratio_leaves(control *c, fluxes *f, params *p) {
    /* Estimate Lignin/N ratio, as this dictates the how plant litter is
    seperated between metabolic and structural pools.

    Returns:
    --------
    lnleaf : float
        lignin:N ratio of leaf

    */
    double lnleaf, nc_leaf_litter;

    nc_leaf_litter = ratio_of_litternc_to_live_leafnc(c, f, p);

    if (float_eq(nc_leaf_litter, 0.0)) {
        /* catch divide by zero if we have no leaves */
        lnleaf = 0.0;
    } else {
        lnleaf = p->ligshoot / p->cfracts / nc_leaf_litter;
    }

    return (lnleaf);

}

double calc_ligin_nratio_fine_roots(control *c, fluxes *f, params *p) {
    /* Estimate Lignin/N ratio, as this dictates the how plant litter is
    seperated between metabolic and structural pools.

    Returns:
    --------
    lnroot : float
        lignin:N ratio of fine root
    */
    double lnroot, nc_root_litter;

    nc_root_litter = ratio_of_litternc_to_live_rootnc(c, f, p);

    if (float_eq(nc_root_litter, 0.0)) {
        /* catch divide by zero if we have no roots */
        lnroot = 0.0;
    } else {
        lnroot = p->ligroot / p->cfracts / nc_root_litter;
    }

    return (lnroot);
}

double ratio_of_litternc_to_live_leafnc(control *c, fluxes *f, params *p) {
    /* ratio of litter N:C to live leaf N:C

    Returns:
    --------
    nc_leaf_litter : float
        N:C ratio of litter to foliage

    */
    double nc_leaf_litter;

    if (c->use_eff_nc){
        nc_leaf_litter = p->liteffnc * (1.0 - p->fretransn);
    } else {
        if (float_eq(f->deadleaves, 0.0)){
            nc_leaf_litter = 0.0;
        } else {
            nc_leaf_litter = f->deadleafn / f->deadleaves;
        }
    }
    
    return (nc_leaf_litter);
}

double ratio_of_litternc_to_live_rootnc(control *c, fluxes *f, params *p) {
    /* ratio of litter N:C to live root N:C

    Returns:
    --------
    nc_root_litter : float
        N:C ratio of litter to live root

    */
    double nc_root_litter;

    if (c->use_eff_nc){
        nc_root_litter = p->liteffnc * p->ncrfac;
    } else {
        if (float_eq(f->deadroots, 0.0)){
            nc_root_litter = 0.0;
        } else{
            nc_root_litter = f->deadrootn / f->deadroots;
        }
    }

    return (nc_root_litter);
}

double metafract(double lig2n) {
    /* Calculate what fraction of the litter will be partitioned to the
    metabolic pool which is given by the lignin:N ratio.

    Parameters:
    -----------
    lig2n : float
        lignin to N ratio

    Returns:
    --------
    metabolic fraction : float
        partitioned fraction to metabolic pool [must be positive]
    */

    /* Original implementation based on Parton et al. */
    return (MAX(0.0, 0.85 - (0.018 * lig2n)));
}


void partition_plant_litter(control *c, fluxes *f, params *p) {
    /* Partition litter from the plant (surface) and roots into metabolic
    and structural pools  */

    double leaf_material, wood_material;
    /*
     * Surface (leaves, branches, stem) Litter
     */

    /* ...to the structural pool*/
    leaf_material = f->deadleaves * (1.0 - p->fmleaf);
    wood_material = f->deadbranch + f->deadstems;
    f->surf_struct_litter = leaf_material + wood_material;        // annual 

    /* ...to the metabolic pool */
    f->surf_metab_litter = f->deadleaves * p->fmleaf;             // annual

    /*
    ** Root Litter
    */

    /* ...to the structural pool */
    f->soil_struct_litter = f->deadroots * (1.0 - p->fmroot);    // annual

    /* ...to the metabolic pool */
    f->soil_metab_litter = f->deadroots * p->fmroot;             // annual
    
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
    

    return;
}

void cfluxes_from_structural_pool(fluxes *f, params *p, state *s) {

    /* Send structural c fluxes to other SOM pools */
    
    double structout_surf = s->structsurf * p->decayrate[0];       // annual
    double structout_soil = s->structsoil * p->decayrate[2];       // annual

    /* C flux surface structural pool -> slow pool */
    f->surf_struct_to_slow = structout_surf * p->ligshoot * 0.7;               // annual

    /* C flux surface structural pool -> active pool */
    f->surf_struct_to_active = structout_surf * (1.0 - p->ligshoot) * 0.55;   // annual

    /* C flux soil structural pool -> slow pool */
    f->soil_struct_to_slow = structout_soil * p->ligroot * 0.7;               // annual

    /* soil structural pool -> active pool */
    f->soil_struct_to_active = structout_soil * (1.0 - p->ligroot) * 0.45;    // annual

    /* Respiration fluxes */

    /* CO2 lost during transfer of structural C to the slow pool */
    f->co2_to_air[0] = (structout_surf *
                        (p->ligshoot * 0.3 + (1.0 - p->ligshoot) * 0.45));   // annual

    /* CO2 lost during transfer structural C  to the active pool */
    f->co2_to_air[1] = (structout_soil *
                        (p->ligroot * 0.3 + (1.0 - p->ligroot) * 0.55));    // annual

    return;
}

void cfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {
    /* Send C from metabolic pools to other SOM pools */

    /* C flux surface metabolic pool -> active pool */
    f->surf_metab_to_active = s->metabsurf * p->decayrate[1] * 0.45;

    /* C flux soil metabolic pool  -> active pool */
    f->soil_metab_to_active = s->metabsoil * p->decayrate[3] * 0.45;

    /* Respiration fluxes */
    f->co2_to_air[2] = s->metabsurf * p->decayrate[1] * 0.55;
    f->co2_to_air[3] = s->metabsoil * p->decayrate[3] * 0.55;
    
    return;
}

void cfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {
    /* Send C fluxes from active pool to other SOM pools */
    
    double activeout = s->activesoil * p->decayrate[4];

    /* C flux active pool -> slow pool */
    f->active_to_slow = activeout * (1.0 - frac_microb_resp - 0.004);

    /* C flux active pool -> passive pool */
    f->active_to_passive = activeout * 0.004;

    /* Respiration fluxes */
    f->co2_to_air[4] = activeout * frac_microb_resp;

    return;
}

void cfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
    /* Send C fluxes from slow pool to other SOM pools */

    double slowout = s->slowsoil * p->decayrate[5];

    /* C flux slow pool -> active pool */
    f->slow_to_active = slowout * 0.42;

    /* slow pool -> passive pool */
    f->slow_to_passive = slowout * 0.03;

    /* Respiration fluxes */
    f->co2_to_air[5] = slowout * 0.55;

    return;
}

void cfluxes_from_passive_pool(fluxes *f, params *p, state *s) {

    /* C flux passive pool -> active pool */
    f->passive_to_active = s->passivesoil * p->decayrate[6] * 0.45;

    /* Respiration fluxes */
    f->co2_to_air[6] = s->passivesoil * p->decayrate[6] * 0.55;

    return;
}

void calculate_soil_respiration(control *c, fluxes *f, params *p, state *s) {
    /* calculate the total soil respiration (heterotrophic) flux, i.e.
    the amount of CO2 released back to the atmosphere */

    /* total CO2 production */
    f->hetero_resp = (f->co2_to_air[0] + f->co2_to_air[1] + f->co2_to_air[2] +
                      f->co2_to_air[3] + f->co2_to_air[4] + f->co2_to_air[5] +
                      f->co2_to_air[6]);
    
    return;
}

void calculate_cpools(control *c, fluxes *f, state *s, params *p) {
    /* Calculate new soil carbon pools. */

    /* Update pools */
    s->structsurf += (f->surf_struct_litter -
                     (f->surf_struct_to_slow + f->surf_struct_to_active +
                      f->co2_to_air[0]));

    s->structsoil += (f->soil_struct_litter -
                     (f->soil_struct_to_slow + f->soil_struct_to_active +
                      f->co2_to_air[1]));
    
    s->metabsurf += (f->surf_metab_litter -
                     (f->surf_metab_to_active + f->co2_to_air[2]));

    s->metabsoil += (f->soil_metab_litter -
                     (f->soil_metab_to_active + f->co2_to_air[3]));
    
    /* store the C SOM fluxes for Nitrogen/Phosphorus calculations */
    f->c_into_active = (f->surf_struct_to_active + f->soil_struct_to_active +
                        f->surf_metab_to_active + f->soil_metab_to_active +
                        f->slow_to_active + f->passive_to_active);

    f->c_into_slow = (f->surf_struct_to_slow + f->soil_struct_to_slow +
                      f->active_to_slow);

    f->c_into_passive = f->active_to_passive + f->slow_to_passive;

    s->activesoil += (f->c_into_active -
                      (f->active_to_slow + f->active_to_passive +
                       f->co2_to_air[4]));

    s->slowsoil +=  (f->c_into_slow -
                    (f->slow_to_active + f->slow_to_passive +
                     f->co2_to_air[5]));

    s->passivesoil += (f->c_into_passive -
                        (f->passive_to_active + f->co2_to_air[6]));
    
    precision_control_soil_c(f, s, p);
    
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
    
    
    return;
}

void precision_control_soil_c(fluxes *f, state *s, params *p) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-08, excess;
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    
    /* C & N state variables */
    if (s->metabsurf < tolerance) {
        excess = s->metabsurf;
        f->surf_metab_to_active = excess * 0.45;
        f->co2_to_air[2] = excess * 0.55;
        s->metabsurf = 0.0;
    }

    if (s->metabsoil < tolerance) {
        excess = s->metabsoil;
        f->soil_metab_to_active = excess * 0.45;
        f->co2_to_air[3] = excess * 0.55;
        s->metabsoil = 0.0;
    }

    return;
}


void calculate_nsoil_flows(control *c, fluxes *f, params *p, state *s) {

    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double nsurf, nsoil, active_nc_slope, slow_nc_slope, passive_nc_slope;

    n_inputs_from_plant_litter(f, p, &nsurf, &nsoil);
    partition_plant_litter_n(c, f, p, nsurf, nsoil);

    /* SOM nitrogen effluxes.  These are assumed to have the source n:c
       ratio prior to the increase of N:C due to co2 evolution. */
    nfluxes_from_structural_pools(f, p, s);
    nfluxes_from_metabolic_pool(f, p, s);
    nfluxes_from_active_pool(f, p, s, frac_microb_resp);
    nfluxes_from_slow_pool(f, p, s);
    nfluxes_from_passive_pool(f, p, s);

    /* gross N mineralisation */
    calculate_n_mineralisation(c, f);
    
    /* calculate N immobilisation */
    calculate_n_immobilisation(f, p, s, &(f->nimmob), &active_nc_slope,
                               &slow_nc_slope, &passive_nc_slope);
    
    /* calculate N net mineralisation */
    calc_n_net_mineralisation(c, f);
    
    /* Update model soil N pools */
    calculate_npools(c, f, p, s, active_nc_slope, slow_nc_slope,
                     passive_nc_slope);

    return;
}

void n_inputs_from_plant_litter(fluxes *f, params *p, double *nsurf,
                              double *nsoil) {
    /* inputs from plant litter.

    surface and soil pools are independent. Structural input flux n:c can
    be either constant or a fixed fraction of metabolic input flux.

    Returns:
    --------
    nsurf : float
        N input from surface pool
    nsoil : float
        N input from soil pool
    */

    /* surface and soil inputs (faeces n goes to abovgrd litter pools) */
    *nsurf = f->deadleafn + f->deadbranchn + f->deadstemn;
    *nsoil = f->deadrootn;

    return;
}

void partition_plant_litter_n(control *c, fluxes *f, params *p, double nsurf,
                              double nsoil) {
    /* Partition litter N from the plant (surface) and roots into metabolic
    and structural pools

    Parameters:
    -----------
    nsurf : float
        N input from surface pool
    nsoil : float
        N input from soil pool
    */

    double c_surf_struct_litter, c_soil_struct_litter;

    /* constant structural input n:c as per century */

    /* n flux -> surface structural pool */
    f->n_surf_struct_litter = f->surf_struct_litter / p->structcn;
    
    /* n flux -> soil structural pool */
    f->n_soil_struct_litter = f->soil_struct_litter / p->structcn;
    
    /* if not enough N for structural, all available N goes to structural */
    if (f->n_surf_struct_litter > nsurf)
      f->n_surf_struct_litter = nsurf;
    if (f->n_soil_struct_litter > nsoil)
      f->n_soil_struct_litter = nsoil;


    /* remaining N goes to metabolic pools */
    f->n_surf_metab_litter = nsurf - f->n_surf_struct_litter;
    f->n_soil_metab_litter = nsoil - f->n_soil_struct_litter;

    return;
}

void nfluxes_from_structural_pools(fluxes *f, params *p, state *s) {
    /* from structural pool */
    double sigwt;
    double structout_surf = s->structsurfn * p->decayrate[0];
    double structout_soil = s->structsoiln * p->decayrate[2];

    sigwt = structout_surf / (p->ligshoot * 0.7 + (1.0 - p->ligshoot) * 0.55);

    /* N flux from surface structural pool -> slow pool */
    f->n_surf_struct_to_slow = sigwt * p->ligshoot * 0.7;

    /* N flux surface structural pool -> active pool */
    f->n_surf_struct_to_active = sigwt * (1.0 - p->ligshoot) * 0.55;

    sigwt = structout_soil / (p->ligroot * 0.7 + (1. - p->ligroot) * 0.45);


    /* N flux from soil structural pool -> slow pool */
    f->n_soil_struct_to_slow = sigwt * p->ligroot * 0.7;

    /* N flux from soil structural pool -> active pool */
    f->n_soil_struct_to_active = sigwt * (1.0 - p->ligroot) * 0.45;

    return;
}

void nfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {

    /* N flux surface metabolic pool -> active pool */
    f->n_surf_metab_to_active = s->metabsurfn * p->decayrate[1];

    /* N flux soil metabolic pool  -> active pool */
    f->n_soil_metab_to_active = s->metabsoiln * p->decayrate[3];
    
    return;
}

void nfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {

    double activeout, sigwt;
    /* N fluxes from active pool */
    activeout = s->activesoiln * p->decayrate[4];
    sigwt = activeout / (1.0 - frac_microb_resp);

    /* N flux active pool -> slow pool */
    f->n_active_to_slow = sigwt * (1.0 - frac_microb_resp - 0.004);

    /* N flux active pool -> passive pool */
    f->n_active_to_passive = sigwt * 0.004;
    
    return;
}

void nfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
    /* N fluxes from slow pools */

    double slowout = s->slowsoiln * p->decayrate[5];
    double sigwt = slowout / 0.45;

    /* C flux slow pool -> active pool */
    f->n_slow_to_active = sigwt * 0.42;

    /* slow pool -> passive pool */
    f->n_slow_to_passive = sigwt * 0.03;

    return;
}

void nfluxes_from_passive_pool(fluxes *f, params *p, state *s) {
    /* N fluxes from passive pool */

    /* C flux passive pool -> active pool */
    f->n_passive_to_active = s->passivesoiln * p->decayrate[6];

    return;
}

void calculate_n_mineralisation(control *c, fluxes *f) {
    /* N gross mineralisation rate is given by the excess of N outflows
    over inflows. Nitrogen mineralisation is the process by which organic
    N is converted to plant available inorganic N, i.e. microbes decompose
    organic N from organic matter to ammonia (NH3) and ammonium (NH4),
    called ammonification.

    Returns:
    --------
    value : float
        Gross N mineralisation
    */
    f->ngross =  (f->n_surf_struct_to_slow + f->n_surf_struct_to_active +
                  f->n_soil_struct_to_slow + f->n_soil_struct_to_active +
                  f->n_surf_metab_to_active + f->n_soil_metab_to_active +
                  f->n_active_to_slow + f->n_active_to_passive +
                  f->n_slow_to_active + f->n_slow_to_passive +
                  f->n_passive_to_active);
  
  /* add diagnostic statement if needed */
  if (c->diagnosis) {
  }
  

    return;
}

void calculate_n_immobilisation(fluxes *f, params *p, state *s, double *nimmob,
                                double *active_nc_slope, double *slow_nc_slope,
                                double *passive_nc_slope) {
    /* N immobilised in new soil organic matter, the reverse of
    mineralisation. Micro-organisms in the soil compete with plants for N.
    Immobilisation is the process by which nitrate and ammonium are taken up
    by the soil organisms and thus become unavailable to the plant
    (->organic N).

    When C:N ratio is high the microorganisms need more nitrogen from
    the soil to decompose the carbon in organic materials. This N will be
    immobilised until these microorganisms die and the nitrogen is
    released.

    General equation for new soil N:C ratio vs Nmin, expressed as linear
    equation passing through point Nmin0, actncmin (etc). Values can be
    Nmin0=0, Actnc0=Actncmin

    if Nmin < Nmincrit:
        New soil N:C = soil N:C (when Nmin=0) + slope * Nmin

    if Nmin > Nmincrit
        New soil N:C = max soil N:C

    NB N:C ratio of new passive SOM can change even if assume Passiveconst

    Returns:
    --------
    nimob : float
        N immobilsed
    */
    double nmin, arg1, arg2, arg3, numer1, numer2, denom;

    /* N:C new SOM - active, slow and passive */
    *active_nc_slope = calculate_nc_slope(p, p->actncmax, p->actncmin);
    *slow_nc_slope = calculate_nc_slope(p, p->slowncmax, p->slowncmin);
    *passive_nc_slope = calculate_nc_slope(p, p->passncmax, p->passncmin);

    /* convert units */
    nmin = p->nmin0 / M2_AS_HA * G_AS_TONNES;
    
    arg1 = (p->passncmin - *passive_nc_slope * nmin) * f->c_into_passive;
    arg2 = (p->slowncmin - *slow_nc_slope * nmin) * f->c_into_slow;
    arg3 = f->c_into_active * (p->actncmin - *active_nc_slope * nmin);
    numer1 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * p->passncmax;
    arg2 = f->c_into_slow * p->slowncmax;
    arg3 = f->c_into_active * p->actncmax;
    numer2 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * *passive_nc_slope;
    arg2 = f->c_into_slow * *slow_nc_slope;
    arg3 = f->c_into_active * *active_nc_slope;
    denom = arg1 + arg2 + arg3;

    /* evaluate N immobilisation in new SOM */
    *nimmob = numer1 + denom * s->inorgn;
    if (*nimmob > numer2)
        *nimmob = numer2;
    
    return;
}


void calc_n_net_mineralisation(control *c, fluxes *f) {
    /* N Net mineralisation from microbial activity */
    f->nmineralisation = f->ngross - f->nimmob + f->nlittrelease;
  
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
  
    
    return;
}

double calculate_nc_slope(params *p, double ncmax, double ncmin) {
    /* Returns N:C ratio of the mineral pool slope

    based on fig 4 of Parton et al 1993. Standard slow pool C:N is different
    to the values in Parton. Bill quotes 12-20, whereas McMurtrie et al '01
    use 10-40.

    Parameters
    ----------
    ncmax : float
        SOM pools maximum N:C
    ncmin: float
        SOM pools minimum N:C

    Returns:
    --------
    value : float
        SOM pool N:C ratio
    */
    double arg1, arg2, conv;

    arg1 = ncmax - ncmin;
    arg2 = p->nmincrit - p->nmin0;
    conv = M2_AS_HA / G_AS_TONNES;

    return (arg1 / arg2 * conv);
}

void calculate_npools(control *c, fluxes *f, params *p, state *s,
                      double active_nc_slope, double slow_nc_slope,
                      double passive_nc_slope) {
    /*
    Update N pools in the soil

    Parameters
    ----------
    active_nc_slope : float
        active NC slope
    slow_nc_slope: float
        slow NC slope
    passive_nc_slope : float
        passive NC slope

    */
    double n_into_active, n_out_of_active, n_into_slow, n_out_of_slow,
           n_into_passive, n_out_of_passive, arg, active_nc, fixn, slow_nc,
           pass_nc;

    /*
        net N release implied by separation of litter into structural
        & metabolic. The following pools only fix or release N at their
        limiting n:c values.
    */

    /* N released or fixed from the N inorganic pool is incremented with
       each call to nc_limit and stored in f->nlittrelease */
    f->nlittrelease = 0.0;
    
    s->structsurfn += (f->n_surf_struct_litter -
                        (f->n_surf_struct_to_slow +
                         f->n_surf_struct_to_active));

    s->structsoiln += (f->n_soil_struct_litter -
                       (f->n_soil_struct_to_slow + f->n_soil_struct_to_active));

    s->structsurfn += nc_limit(f, s->structsurf, s->structsurfn,
                               1.0/p->structcn, 1.0/p->structcn);
    
    s->structsoiln += nc_limit(f, s->structsoil, s->structsoiln,
                               1.0/p->structcn, 1.0/p->structcn);
    
    s->metabsurfn += f->n_surf_metab_litter - f->n_surf_metab_to_active;
    
    s->metabsurfn += nc_limit(f, s->metabsurf, s->metabsurfn,1.0/25.0, 1.0/10.0);
    
    s->metabsoiln += (f->n_soil_metab_litter - f->n_soil_metab_to_active);
    s->metabsoiln += nc_limit(f, s->metabsoil, s->metabsoiln, 1.0/25.0,
                              1.0/10.0);
    
    /* When nothing is being added to the metabolic pools, there is the
       potential scenario with the way the model works for tiny bits to be
       removed with each timestep. Effectively with time this value which is
       zero can end up becoming zero but to a silly decimal place */
    precision_control_soil_n(f, s, p);

    /* Update SOM pools */
    n_into_active = (f->n_surf_struct_to_active + f->n_soil_struct_to_active +
                     f->n_surf_metab_to_active + f->n_soil_metab_to_active +
                     f->n_slow_to_active + f->n_passive_to_active);

    n_out_of_active = f->n_active_to_slow + f->n_active_to_passive;

    n_into_slow = (f->n_surf_struct_to_slow + f->n_soil_struct_to_slow +
                   f->n_active_to_slow);

    n_out_of_slow = f->n_slow_to_active + f->n_slow_to_passive;
    n_into_passive = f->n_active_to_passive + f->n_slow_to_passive;
    n_out_of_passive = f->n_passive_to_active;

    /* N:C of the SOM pools increases linearly btw prescribed min and max
       values as the Nconc of the soil increases. */
    arg = s->inorgn - p->nmin0 / M2_AS_HA * G_AS_TONNES;

    /* active */
    active_nc = p->actncmin + active_nc_slope * arg;
    if (active_nc > p->actncmax)
        active_nc = p->actncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_active, n_into_active, active_nc);
    s->activesoiln += n_into_active + fixn - n_out_of_active;

    /* slow */
    slow_nc = p->slowncmin + slow_nc_slope * arg;
    if (slow_nc > p->slowncmax)
        slow_nc = p->slowncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_slow, n_into_slow, slow_nc);
    s->slowsoiln += n_into_slow + fixn - n_out_of_slow;

    /* passive, update passive pool only if passiveconst=0 */
    pass_nc = p->passncmin + passive_nc_slope * arg;
    if (pass_nc > p->passncmax)
        pass_nc = p->passncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_passive, n_into_passive, pass_nc);
    s->passivesoiln += n_into_passive + fixn - n_out_of_passive;
    
    /* Daily increment of soil inorganic N pool, diff btw in and effluxes
       (grazer urine n goes directly into inorganic pool) nb inorgn may be
       unstable if rateuptake is large */
    s->inorgn += f->ninflow + f->nmineralisation - f->nloss - f->nuptake;  
    
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
    
    

    return;
}


double nc_limit(fluxes *f, double cpool, double npool, double ncmin,
                double ncmax) {
    /* Release N to 'Inorgn' pool or fix N from 'Inorgn', in order to keep
    the  N:C ratio of a litter pool within the range 'ncmin' to 'ncmax'.

    Parameters:
    -----------
    cpool : float
        various C pool (state)
    npool : float
        various N pool (state)
    ncmin : float
        maximum N:C ratio
    ncmax : float
        minimum N:C ratio

    Returns:
    --------
    fix/rel : float
        amount of N to be added/released from the inorganic pool

    */
    double rel, fix;
    double nmax = cpool * ncmax;
    double nmin = cpool * ncmin;

    if (npool > nmax) {
        /* release */
        rel = npool - nmax;
        f->nlittrelease += rel;
        return (-rel);
    } else if (npool < nmin) {
        /* fix */
        fix = nmin - npool;
    
        if (f->nlittrelease < fix) {
            fix = f->nlittrelease;
           }
      
        f->nlittrelease -= fix;

        return (fix);
    } else {
        return (0.0);
    }
}

double nc_flux(double cflux, double nflux, double nc_ratio) {
    /*
    Release N to Inorganic pool or fix N from the Inorganic pool in order
    to normalise the N:C ratio of a net flux

    Parameters:
    -----------
    cflux : float
        C flux into SOM pool
    nflux : float
        N flux into SOM pool
    nc_ratio : float
        preferred N:C ratio

    Returns:
        fix : float
        Returns the amount of N required to be fixed
    */

    return (cflux * nc_ratio) - nflux;
}


void precision_control_soil_n(fluxes *f, state *s, params *p) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-08, excess;

    if (s->metabsurfn < tolerance) {
        excess = s->metabsurfn;
        f->n_surf_metab_to_active = excess;
        s->metabsurfn = 0.0;
    }

    if (s->metabsoiln < tolerance) {
        excess = s->metabsoiln;
        f->n_soil_metab_to_active = excess;
        s->metabsoiln = 0.0;
    }

    return;
}

void calculate_psoil_flows(control *c, fluxes *f, params *p, state *s) {

    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double psurf, psoil, active_pc_slope, slow_pc_slope, passive_pc_slope;

    p_inputs_from_plant_litter(f, p, &psurf, &psoil);
    partition_plant_litter_p(c, f, p, psurf, psoil);

    /* SOM phosphorus effluxes.  These are assumed to have the source p:c
    ratio prior to the increase of P:C due to co2 evolution. */
    pfluxes_from_structural_pools(f, p, s);
    pfluxes_from_metabolic_pool(f, p, s);
    pfluxes_from_active_pool(f, p, s, frac_microb_resp);
    pfluxes_from_slow_pool(f, p, s);
    pfluxes_from_passive_pool(f, p, s);

    /* calculate P parent influxe to mineral P */
    calculate_p_parent_fluxes(c, f, p, s);

    /* gross P mineralisation */
    calculate_p_mineralisation(f);

    /* calculate P immobilisation */
    calculate_p_immobilisation(f, p, s, &(f->pimmob), &active_pc_slope,
                             &slow_pc_slope, &passive_pc_slope);

    /* calculate P net mineralisation*/
    calc_p_net_mineralisation(c, f);

    /* SIM phosphorus dynamics */
    calculate_p_ssorb_to_avl(s, f, p, c);
    calculate_p_ssorb_to_occ(s, f, p);
    calculate_p_avl_to_ssorb(s, f, p);

    /* Update model soil P pools */
    calculate_ppools(c, f, p, s, active_pc_slope, slow_pc_slope,
                     passive_pc_slope);
    
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
  
    return;
}

void p_inputs_from_plant_litter(fluxes *f, params *p, double *psurf,
                                double *psoil) {
    /* inputs from plant litter.

    surface and soil pools are independent. Structural input flux p:c can
    be either constant or a fixed fraction of metabolic input flux.

    Returns:
    --------
    psurf : float
    P input from surface pool
    psoil : float
    P input from soil pool
    */

    /* surface and soil inputs (faeces p goes to abovgrd litter pools) */
    *psurf = f->deadleafp + f->deadbranchp + f->deadstemp;
    *psoil = f->deadrootp;

    return;
}

void partition_plant_litter_p(control *c, fluxes *f, params *p, double psurf,
                              double psoil) {
    /* Partition litter P from the plant (surface) and roots into metabolic
    and structural pools

    Parameters:
    -----------
    psurf : float
    P input from surface pool
    psoil : float
    P input from soil pool
    */

    double c_surf_struct_litter, c_soil_struct_litter;

    /* constant structural input p:c as per century */

    /* p flux -> surface structural pool */
    f->p_surf_struct_litter = f->surf_struct_litter / p->structcp;
    
    /* p flux -> soil structural pool */
    f->p_soil_struct_litter = f->soil_struct_litter / p->structcp;
    
    /* if not enough P for structural, all available P goes to structural */
    if (f->p_surf_struct_litter > psurf)
      f->p_surf_struct_litter = psurf;
    if (f->p_soil_struct_litter > psoil)
      f->p_soil_struct_litter = psoil;
    

    /* remaining P goes to metabolic pools */
    f->p_surf_metab_litter = psurf - f->p_surf_struct_litter;
    f->p_soil_metab_litter = psoil - f->p_soil_struct_litter;

    return;
}

void pfluxes_from_structural_pools(fluxes *f, params *p, state *s) {
    /* from structural pool */
    double sigwt;
    double structout_surf = s->structsurfp * p->decayrate[0];
    double structout_soil = s->structsoilp * p->decayrate[2];

    sigwt = structout_surf / (p->ligshoot * 0.7 + (1.0 - p->ligshoot) * 0.55);

    /* P flux from surface structural pool -> slow pool */
    f->p_surf_struct_to_slow = sigwt * p->ligshoot * 0.7;

    /* P flux surface structural pool -> active pool */
    f->p_surf_struct_to_active = sigwt * (1.0 - p->ligshoot) * 0.55;

    sigwt = structout_soil / (p->ligroot * 0.7 + (1. - p->ligroot) * 0.45);


    /* P flux from soil structural pool -> slow pool */
    f->p_soil_struct_to_slow = sigwt * p->ligroot * 0.7;

    /* N flux from soil structural pool -> active pool */
    f->p_soil_struct_to_active = sigwt * (1.0 - p->ligroot) * 0.45;

    return;
}

void pfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {

    /* P flux surface metabolic pool -> active pool */
    f->p_surf_metab_to_active = s->metabsurfp * p->decayrate[1];

    /* P flux soil metabolic pool  -> active pool */
    f->p_soil_metab_to_active = s->metabsoilp * p->decayrate[3];

    return;
}


void pfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {

    double activeout, sigwt;
    /* P fluxes from active pool */
    activeout = s->activesoilp * p->decayrate[4];
    sigwt = activeout / (1.0 - frac_microb_resp);

    /* P flux active pool -> slow pool */
    f->p_active_to_slow = sigwt * (1.0 - frac_microb_resp - 0.004);

    /* P flux active pool -> passive pool */
    f->p_active_to_passive = sigwt * 0.004;

    return;
}

void pfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
    /* P fluxes from slow pools */

    double slowout = s->slowsoilp * p->decayrate[5];
    double sigwt = slowout / 0.45;

    /* P flux slow pool -> active pool */
    f->p_slow_to_active = sigwt * 0.42;

    /* slow pool -> passive pool */
    f->p_slow_to_passive = sigwt * 0.03;

    return;
}

void pfluxes_from_passive_pool(fluxes *f, params *p, state *s) {
    /* P fluxes from passive pool */

    /* P flux passive pool -> active pool */
    f->p_passive_to_active = s->passivesoilp * p->decayrate[6];

    return;
}

void calculate_p_parent_fluxes(control *c, fluxes *f, params *p, state *s) {
    /*
        Calculate weathering of parent P materials, i.e.
        the fluxes enterring into mineral P pool;

        Fluxes in = out so that parent P pool is a constant pool;
    */

    /* parent material weathering */
    f->p_par_to_avl = p->p_rate_par_weather * s->inorgparp;
  
  /* add diagnostic statement if needed */
  if (c->diagnosis) {
  }
  

    return;
}

void calculate_p_mineralisation(fluxes *f) {
    /* P gross mineralisation rate is given by the excess of P outflows
    over inflows. P mineralisation is the process by which organic P is
    converted to plant available inorganic P, i.e. the microbial conversion
    of organic P to H2PO4- or HPO42- forms of plant available P known as
    orthophosphates.

    Returns:
    --------
    value : float
    Gross P mineralisation
    Unit: t/ha/d
    */
    f->pgross =  (f->p_surf_struct_to_slow + f->p_surf_struct_to_active +
                  f->p_soil_struct_to_slow + f->p_soil_struct_to_active +
                  f->p_surf_metab_to_active + f->p_soil_metab_to_active +
                  f->p_active_to_slow + f->p_active_to_passive +
                  f->p_slow_to_active + f->p_slow_to_passive +
                  f->p_passive_to_active);

    return;
}

void calculate_p_immobilisation(fluxes *f, params *p, state *s, double *pimmob,
                                double *active_pc_slope, double *slow_pc_slope,
                                double *passive_pc_slope) {
    /* P immobilised in new soil organic matter, the reverse of
    mineralisation. Micro-organisms in the soil compete with plants for P.
    Immobilisation occurs when plant available P forms are consumed by microbes, turning
    the P into organic P forms that are not available to plants.

    When C:P ratio is high the microorganisms need more P from
    the soil to decompose the carbon in organic materials. This P will be
    immobilised until these microorganisms die and the P is
    released.

    General equation for new soil P:C ratio vs Pmin, expressed as linear
    equation passing through point Pmin0, actpcmin (etc). Values can be
    Pmin0=0, Actpc0=Actpcmin

    if Pmin < Pmincrit:
    New soil P:C = soil P:C (when Pmin=0) + slope * Pmin

    if Pmin > Pmincrit
    New soil P:C = max soil P:C

    NB P:C ratio of new passive SOM can change even if assume Passiveconst

    Returns:
    --------
    pimob : float
    P immobilsed
    */
    double pmin, arg1, arg2, arg3, numer1, numer2, denom;

    /* P:C new SOM - active, slow and passive */
    *active_pc_slope = calculate_pc_slope(p, p->actpcmax, p->actpcmin);
    *slow_pc_slope = calculate_pc_slope(p, p->slowpcmax, p->slowpcmin);
    *passive_pc_slope = calculate_pc_slope(p, p->passpcmax, p->passpcmin);

    /* convert units */
    pmin = p->pmin0 / M2_AS_HA * G_AS_TONNES;

    arg1 = (p->passpcmin - *passive_pc_slope * pmin) * f->c_into_passive;
    arg2 = (p->slowpcmin - *slow_pc_slope * pmin) * f->c_into_slow;
    arg3 = f->c_into_active * (p->actpcmin - *active_pc_slope * pmin);
    numer1 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * p->passpcmax;
    arg2 = f->c_into_slow * p->slowpcmax;
    arg3 = f->c_into_active * p->actpcmax;
    numer2 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * *passive_pc_slope;
    arg2 = f->c_into_slow * *slow_pc_slope;
    arg3 = f->c_into_active * *active_pc_slope;
    denom = arg1 + arg2 + arg3;

    /* evaluate P immobilisation in new SOM */
    *pimmob = numer1 + denom * s->inorgavlp;
    if (*pimmob > numer2)
        *pimmob = numer2;

    return;
}

void calc_p_net_mineralisation(control *c, fluxes *f) {
    /*
        P Net mineralisation from microbial activity
    */
    f->pmineralisation = f->pgross - f->pimmob + f->plittrelease;
  
  /* add diagnostic statement if needed */
  if (c->diagnosis) {
  }
  

    return;
}

double calculate_pc_slope(params *p, double pcmax, double pcmin) {
    /* Returns P:C ratio of the mineral pool slope
    Need to check back for good relationships; Currently using olde NC relationship
    based on Parton et al., 1993

    Parameters
    ----------
    pcmax : float
    SOM pools maximum P:C
    pcmin: float
    SOM pools minimum P:C

    Returns:
    --------
    value : float
    SOM pool P:C ratio
    */
    double arg1, arg2, conv;

    arg1 = pcmax - pcmin;
    arg2 = p->pmincrit - p->pmin0;
    conv = M2_AS_HA / G_AS_TONNES;

    return (arg1 / arg2 * conv);
}


void calculate_p_avl_to_ssorb(state *s, fluxes *f, params *p) {
  
  /* P flux from sorbed pool to strongly sorbed P pool */
    f->p_avl_to_ssorb = p->k1 * s->inorgavlp;
  
  return;
}

void calculate_p_ssorb_to_avl(state *s, fluxes *f, params *p, control *c) {
    /*
        calculate P transfer from strongly sorbed P pool to
        sorbed P pool;

    */
    
    if (s->inorgssorbp > 0.0) {
        f->p_ssorb_to_avl = p->k2 * s->inorgssorbp;
    } else {
        f->p_ssorb_to_avl = 0.0;
    }

    return;
}

void calculate_p_ssorb_to_occ(state *s, fluxes *f, params *p) {

    /* P flux from strongly sorbed pool to occluded P pool */
    if (s->inorgssorbp > 0.0) {
        f->p_ssorb_to_occ = p->k3 * s->inorgssorbp;
    } else {
        f->p_ssorb_to_occ = 0.0;
    }

    return;
}


void calculate_ppools(control *c, fluxes *f, params *p, state *s,
                      double active_pc_slope, double slow_pc_slope,
                      double passive_pc_slope) {
    /*
        Update P pools in the soil

        Parameters
        ----------
        active_pc_slope : float
            active PC slope
        slow_pc_slope: float
            slow PC slope
        passive_pc_slope : float
            passive PC slope

    */

    double p_into_active, p_out_of_active, p_into_slow, p_out_of_slow,
           p_into_passive, p_out_of_passive, arg, active_pc, fixp, slow_pc,
           pass_pc;

    double tot_avl_in, tot_avl_out, net_parent;

    /*
        net P release implied by separation of litter into structural
        & metabolic. The following pools only fix or release P at their
        limiting p:c values.
    */

    /* P released or fixed from the P inorganic labile pool is incremented with
    each call to pc_limit and stored in f->plittrelease */
    f->plittrelease = 0.0;

    s->structsurfp += (f->p_surf_struct_litter -
                      (f->p_surf_struct_to_slow +
                       f->p_surf_struct_to_active));

    s->structsoilp += (f->p_soil_struct_litter -
                      (f->p_soil_struct_to_slow +
                       f->p_soil_struct_to_active));

    s->structsurfp += pc_limit(f, s->structsurf, s->structsurfp,
                               1.0/p->structcp, 1.0/p->structcp);
    s->structsoilp += pc_limit(f, s->structsoil, s->structsoilp,
                               1.0/p->structcp, 1.0/p->structcp);


    /* pcmin & pcmax from Parton 1989 fig 2 */
    s->metabsurfp += f->p_surf_metab_litter - f->p_surf_metab_to_active;
    s->metabsurfp += pc_limit(f, s->metabsurf, s->metabsurfp,
                              1.0/150.0, 1.0/80.0);

    /* pcmin & pcmax from Parton 1989 fig 2 */
    s->metabsoilp += (f->p_soil_metab_litter - f->p_soil_metab_to_active);
    s->metabsoilp += pc_limit(f, s->metabsoil, s->metabsoilp,
                              1.0/150.0, 1.0/80.0);


    /* When nothing is being added to the metabolic pools, there is the
    potential scenario with the way the model works for tiny bits to be
    removed with each timestep. Effectively with time this value which is
    zero can end up becoming zero but to a silly decimal place */
    precision_control_soil_p(f, s, p);

    /* Update SOM pools */
    p_into_active = (f->p_surf_struct_to_active + f->p_soil_struct_to_active +
    f->p_surf_metab_to_active + f->p_soil_metab_to_active +
    f->p_slow_to_active + f->p_passive_to_active);

    p_out_of_active = f->p_active_to_slow + f->p_active_to_passive;

    p_into_slow = (f->p_surf_struct_to_slow + f->p_soil_struct_to_slow +
                   f->p_active_to_slow);

    p_out_of_slow = f->p_slow_to_active + f->p_slow_to_passive;
    p_into_passive = f->p_active_to_passive + f->p_slow_to_passive;
    p_out_of_passive = f->p_passive_to_active;

    // P:C of the SOM pools increases linearly btw prescribed min and max
    // values as the Pconc of the soil increases.
    arg = s->inorgavlp - p->pmin0 / M2_AS_HA * G_AS_TONNES;

    /* active */
    active_pc = p->actpcmin + active_pc_slope * arg;
    if (active_pc > p->actpcmax)
        active_pc = p->actpcmax;

    // release P to Inorganic labile pool or fix P from the Inorganic pool in order
    // to normalise the P:C ratio of a net flux
    fixp = pc_flux(f->c_into_active, p_into_active, active_pc);
    s->activesoilp += p_into_active + fixp - p_out_of_active;

    /* slow */
    slow_pc = p->slowpcmin + slow_pc_slope * arg;
    if (slow_pc > p->slowpcmax)
        slow_pc = p->slowpcmax;

    /* release P to Inorganic pool or fix P from the Inorganic pool in order
    to normalise the P:C ratio of a net flux */
    fixp = pc_flux(f->c_into_slow, p_into_slow, slow_pc);
    s->slowsoilp += p_into_slow + fixp - p_out_of_slow;

    /* passive, update passive pool only if passiveconst=0 */
    pass_pc = p->passpcmin + passive_pc_slope * arg;
    if (pass_pc > p->passpcmax)
        pass_pc = p->passpcmax;

    /* release P to Inorganic pool or fix P from the Inorganic pool in order
    to normalise the P:C ratio of a net flux */
    fixp = pc_flux(f->c_into_passive, p_into_passive, pass_pc);
    s->passivesoilp += p_into_passive + fixp - p_out_of_passive;

    /* Daily increment of soil inorganic available P pool (lab + sorb) */
    tot_avl_in = f->p_par_to_avl + f->pmineralisation + f->p_ssorb_to_avl;
    
    if (s->inorgavlp > 0) {
      tot_avl_out = f->puptake + f->ploss + f->p_avl_to_ssorb;
    } else {
      f->puptake = 0.0;
      f->ploss = 0.0;
      f->p_avl_to_ssorb = 0.0;
      tot_avl_out = f->puptake + f->ploss + f->p_avl_to_ssorb;
    }
    
    s->inorgavlp += tot_avl_in - tot_avl_out;

    /* Daily increment of soil inorganic secondary P pool (strongly sorbed) */
    s->inorgssorbp += f->p_avl_to_ssorb - f->p_ssorb_to_occ - f->p_ssorb_to_avl;

    /* Daily increment of soil inorganic occluded P pool */
    s->inorgoccp += f->p_ssorb_to_occ;

    /* Daily increment of soil inorganic parent P pool */
    s->inorgparp += f->p_atm_dep - f->p_par_to_avl;   
    
    /* add diagnostic statement if needed */
    if (c->diagnosis) {
    }
    
    
    return;
}


double pc_limit(fluxes *f, double cpool, double ppool, double pcmin,
                double pcmax) {
    /*
        Release P to 'Inorgavlp' pool or fix P from 'Inorgavlp', in order to keep
        the  P:C ratio of a litter pool within the range 'pcmin' to 'pcmax'.

        Parameters:
        -----------
        cpool : float
            various C pool (state)
        ppool : float
            various P pool (state)
        pcmin : float
            min P:C ratio
        pcmax : float
            max P:C ratio

        Returns:
        --------
        fix/rel : float
            amount of P to be added/released from the inorganic pool

    */
    double rel, fix;
    double pmax = cpool * pcmax;
    double pmin = cpool * pcmin;

    if (ppool > pmax) {
        /* release */
        rel = ppool - pmax;
        f->plittrelease += rel;
        return (-rel);
    } else if (ppool < pmin) {
        /* fix */
        fix = pmin - ppool;
        f->plittrelease -= fix;
        return (fix);
    } else {
        return (0.0);
    }
}

double pc_flux(double cflux, double pflux, double pc_ratio) {
    /*
        Release P to Inorganic pool or fix P from the Inorganic pool in order
        to normalise the P:C ratio of a net flux

        Parameters:
        -----------
        cflux : float
            C flux into SOM pool
        pflux : float
            P flux into SOM pool
        pc_ratio : float
            preferred P:C ratio

        Returns:
            fix : float
            Returns the amount of P required to be fixed
    */

    return (cflux * pc_ratio) - pflux;
}


void precision_control_soil_p(fluxes *f, state *s, params *p) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-08, excess;

    if (s->metabsurfp < tolerance) {
        excess = s->metabsurfp;
        f->p_surf_metab_to_active = excess;
        s->metabsurfp = 0.0;
    }

    if (s->metabsoilp < tolerance) {
        excess = s->metabsoilp;
        f->p_soil_metab_to_active = excess;
        s->metabsoilp = 0.0;
    }

    return;
}
