#include "read_met_file.h"



void read_annual_met_data_simple(char **argv, control *c, met *m, params *p)
{
    double current_yr = -999.9;
    
    /* unpack met forcing */
    m->year = 1;
    m->co2 = p->co2_in;
    m->par = p->I0;
    
    m->ndep = p->ndep_in;
    m->nfix = p->nfix_in;
    m->pdep = p->pdep_in;
    m->tsoil = p->tsoil_in;
    
    return;
}

