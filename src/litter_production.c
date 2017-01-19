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
    ncrlit = s->rootnc * (1.0 - p->rretrans);
    
    /* litter P:C ratios, roots and shoot */
    pcflit = s->shootpc * (1.0 - p->fretransp);
    pcrlit = s->rootpc * (1.0 - p->rretrans);

    /* C litter production */
    f->deadroots = *rdecay * s->root;
    f->deadcroots = p->crdecay * s->croot;
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
    f->deadcrootn = p->crdecay * s->crootn * (1.0 - p->cretrans);
    f->deadbranchn = p->bdecay * s->branchn * (1.0 - p->bretrans);

    f->deadrootp = f->deadroots * pcrlit;
    f->deadcrootp = p->crdecay * s->crootp * (1.0 - p->cretrans);
    f->deadbranchp = p->bdecay * s->branchp * (1.0 - p->bretrans);

    /* N in stemwood litter - only mobile n is retranslocated */
    f->deadstemn = p->wdecay * (s->stemnimm + s->stemnmob * \
                    (1.0 - p->wretrans));

    /* P in stemwood litter - only mobile p is retranslocated */
    f->deadstemp = p->wdecay * (s->stempimm + s->stempmob * \
                   (1.0 - p->wretrans));
    
    f->ceaten = 0.0;
    f->neaten = 0.0;
    f->peaten = 0.0;
        
    return;

}
