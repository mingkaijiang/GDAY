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
