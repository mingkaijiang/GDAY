/* ============================================================================
* Calculate C and N litter production
*
* Litter production for each pool is assumed to be proportional to biomass
* pool size.
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
#include "litter_production.h"

void calculate_litterfall(control *c, fluxes *f, params *p, state *s,
                          int doy, double *fdecay, double *rdecay) {

    double  ncflit, ncrlit;
    double  pcflit, pcrlit;

    /* litter N:C ratios, roots and shoot */
    ncflit = s->shootnc * (1.0 - p->fretrans);
    ncrlit = s->rootnc;
    
    /* litter P:C ratios, roots and shoot */
    pcflit = s->shootpc * (1.0 - p->fretransp);
    pcrlit = s->rootpc;

    /* C litter production */
    f->deadroots = *rdecay * s->root;
    f->deadstems = p->wdecay * s->stem;
    f->deadbranch = p->bdecay * s->branch;
    f->deadsapwood = (p->wdecay + p->sapturnover) * s->sapwood;
    
    // fprintf(stderr, "f->deadroots %f\n", f->deadroots);
    // fprintf(stderr, "rdecay %f\n", *rdecay);
    // fprintf(stderr, "s->root %f\n", s->root);

    f->deadleaves = *fdecay * s->shoot;
    
    /* N litter production */
    f->deadleafn = f->deadleaves * ncflit;

    /* P litter production */
    f->deadleafp = f->deadleaves * pcflit;

    /* Assuming fraction is retranslocated before senescence, i.e. a fracion
       of nutrients is stored within the plant */
    f->deadrootn = f->deadroots * ncrlit;

    f->deadrootp = f->deadroots * pcrlit;

    /* N in stemwood litter - only mobile n is retranslocated */
    f->deadstemn = p->wdecay * (s->stemnimm + s->stemnmob);

    /* P in stemwood litter - only mobile p is retranslocated */
    f->deadstemp = p->wdecay * (s->stempimm + s->stempmob);
        
    return;

}
