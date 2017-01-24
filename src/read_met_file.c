#include "read_met_file.h"



void read_daily_met_data_simple(char **argv, control *c, met_arrays *ma, params *p)
{
  int    file_len = 0;
  int    i = 0;
  double current_yr = -999.9;
  
  /* work out how big the file is */
  file_len = 1;

  c->total_num_days = file_len;
  
  /* allocate memory for meteorological arrays */
  if ((ma->year = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for year array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->tsoil = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for tsoil array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->co2 = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for co2 array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->ndep = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for ndep array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->nfix = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for nfix array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->pdep = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for pdep array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->par = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for par array\n");
    exit(EXIT_FAILURE);
  }
  
  
  i = 0;
  c->num_years = 0;
  
  /* unpack met forcing */
  ma->year[i] = 17;
  ma->co2[i] = p->co2_in;
  ma->par[i] = p->I0/365.25;    // 3 GJ m-2 yr-1 to MJ m-2 d-1;
  
  ma->ndep[i] = p->ndep_in;
  ma->nfix[i] = p->nfix_in;
  ma->pdep[i] = p->pdep_in;
  ma->tsoil[i] = p->tsoil_in;
 
    /* Build an array of the unique years as we loop over the input file */
    if (current_yr != ma->year[i]) {
      c->num_years++;
      current_yr = ma->year[i];
    }
    i++;
    
  return;
}

