#include "read_met_file.h"



void read_daily_met_data_simple(char **argv, control *c, met *m, params *p, fluxes *f)
{
    int    file_len = 0;
    double current_yr = -999.9;
  
    /* work out how big the file is */
    file_len = 1;
    
    c->total_num_days = file_len;
    
    c->num_years = 0;
    
    /* unpack met forcing */
    m->year = 17;
    m->co2 = p->co2_in;
    m->par = p->I0/365.25;    // 3 GJ m-2 yr-1 to MJ m-2 d-1;
    
    m->ndep = p->ndep_in;
    m->nfix = p->nfix_in;
    m->pdep = p->pdep_in;
    m->tsoil = p->tsoil_in;
    
    /* Build an array of the unique years as we loop over the input file */
    if (current_yr != m->year) {
      c->num_years++;
      current_yr = m->year;
    }

    return;
}

