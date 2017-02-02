#ifndef PLANT_GROWTH_H
#define PLANT_GROWTH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "utilities.h"
#include "photosynthesis.h"

/* C stuff */
void    calc_annual_growth(control *, fluxes *, met *,
                        nrutil *, params *, state *);
void    carbon_allocation(control *, fluxes *, params *, state *, double);
void    calc_carbon_allocation_fracs(control *c, fluxes *, params *, state *,
                                     double);
void    update_plant_state(control *, fluxes *, params *, state *);
void    precision_control(fluxes *, state *);
void    carbon_annual_production(control *, fluxes *, met *m, params *, state *);

void    calculate_cnp_wood_ratios(control *c, params *, state *, double, double,
                                  double, double *, 
                                  double *, double *,double *, double *,
                                  double *);

/* N stuff */
void    np_allocation(control *c, fluxes *, params *, state *, double,
                     double, double, double, double,
                     double);
double calculate_nuptake(control *, params *, state *);

double nitrogen_retrans(control *, fluxes *, params *, state *);


/* P stuff */
double calculate_growth_stress_limitation(params *, state *, control *);
double calculate_puptake(control *, params *, state *, fluxes *);
double phosphorus_retrans(control *, fluxes *, params *, state *);

void cut_back_production(control *, fluxes *, params *, state *, double,
                        double, double, double);

#endif /* PLANT_GROWTH */
