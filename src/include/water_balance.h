#ifndef WATER_BALANCE_H
#define WATER_BALANCE_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"

void    initialise_soils_day(control *, fluxes *, params *, state *);

double  calc_sat_water_vapour_press(double);
double  calc_sw_modifier(double, double, double);

double *get_soil_fracs(char *);
double  calc_beta(double, double, double, double, double);
void    get_soil_params(char *, double *, double *);
void    calc_soil_params(double *, double *, double *,
                        double *, double *, double *);
void    calculate_soil_water_fac(control *, params *, state *);

#endif /* WATER_BALANCE */
