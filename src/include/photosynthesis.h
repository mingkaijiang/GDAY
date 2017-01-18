#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"

/* Daily funcs */
void simple_photosynthesis(control *, fluxes *, met *, params *, state *);

double  arrh(double, double, double, double);

double  lue_simplified(params *, state *, double co2);


#endif /* PHOTOSYNTHESIS */
