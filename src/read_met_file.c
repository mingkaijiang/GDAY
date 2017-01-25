#include "read_met_file.h"



void read_daily_met_data_simple(char **argv, control *c, met *m, params *p, fluxes *f)
{
    double current_yr = -999.9;
    
    /* unpack met forcing */
    m->year = 1;
    m->co2 = p->co2_in;
    m->par = p->I0/365.25;
    
    m->ndep = p->ndep_in;
    m->nfix = p->nfix_in;
    m->pdep = p->pdep_in;
    m->tsoil = p->tsoil_in;
    
    /* Build an array of the unique years as we loop over the input file */
    if (current_yr != m->year) {
      p->num_years++;
      current_yr = m->year;
    }

    return;
}

