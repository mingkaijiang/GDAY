#include "read_met_file.h"

void read_daily_met_data(char **argv, control *c, met_arrays *ma)
{
    FILE  *fp;
    char   line[STRING_LENGTH];
    int    file_len = 0;
    int    i = 0;
    int    nvars = 9;
    int    skipped_lines = 0;
    double current_yr = -999.9;

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open daily Met file %s for read\n",
                c->met_fname);
		exit(EXIT_FAILURE);
	 }

    /* work out how big the file is */
    file_len = 0;
//    while (fgets(line, STRING_LENGTH, fp) != NULL) {
//        /* ignore comment line */
//        if (*line == '#')
//            continue;
//        file_len++;
//    }
//    rewind(fp);
    c->total_num_days = file_len;

    /* allocate memory for meteorological arrays */
    if ((ma->year = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for year array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->prjday = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for prjday array\n");
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

    if ((ma->par_am = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par_pm = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par_pm array\n");
		exit(EXIT_FAILURE);
    }


    i = 0;
    c->num_years = 0;
    skipped_lines = 0;
    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#') {
            skipped_lines++;
            continue;
        }

        if (sscanf(line, "%lf,%lf,\
                          %lf,\
                          %lf,%lf,\
                          %lf,%lf,\
                          %lf,%lf",\
                          &(ma->year[i]), &(ma->prjday[i]), \
                          &(ma->tsoil[i]), \
                          &(ma->co2[i]), &(ma->ndep[i]), \
                          &(ma->nfix[i]),  &(ma->pdep[i]), \
                          &(ma->par_am[i]), &(ma->par_pm[i])) != nvars) {
            fprintf(stderr, "%s: badly formatted input in met file on line %d %d\n", \
                    *argv, (int)i+1+skipped_lines, nvars);
            exit(EXIT_FAILURE);
        }

        /* Build an array of the unique years as we loop over the input file */
        if (current_yr != ma->year[i]) {
            c->num_years++;
            current_yr = ma->year[i];
        }
        i++;
    }
    fclose(fp);
    return;
}



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
  
  if ((ma->prjday = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for prjday array\n");
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
  
  if ((ma->par_am = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for par_am array\n");
    exit(EXIT_FAILURE);
  }
  
  if ((ma->par_pm = (double *)calloc(file_len, sizeof(double))) == NULL) {
    fprintf(stderr,"Error allocating space for par_pm array\n");
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

