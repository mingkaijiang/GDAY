#ifndef LITTER_H
#define LITTER_H

#include "gday.h"

/* litter stuff */
float  decay_in_dry_soils(double, double, params *, state *);
void   calculate_litterfall(control *, fluxes *, params *, state *, int,
                            double *, double *);

#endif /* LITTER */
