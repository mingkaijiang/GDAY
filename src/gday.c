/* ============================================================================
* Generic Decomposition And Yield (GDAY) model.
*
* Simple version for quasi-equilibrium analysis
* 
* This version runs at annual timestep by using met forcing data averaged annually
* and fluxes summed at the end of each year.
*
* Paramaeter descriptions are in gday.h
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe & Mingkai Jiang
*
* DATE:
*   18.01.2017
*
* =========================================================================== */

#include "gday.h"

int main(int argc, char **argv)
{
    int error = 0;

    /*
     * Setup structures, initialise stuff, e.g. zero fluxes.
     */
    control *c;
    fluxes *f;
    met *m;
    params *p;
    state *s;
    nrutil *nr;

    c = (control *)malloc(sizeof(control));
    if (c == NULL) {
        fprintf(stderr, "control structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    f = (fluxes *)malloc(sizeof(fluxes));
    if (f == NULL) {
    	fprintf(stderr, "fluxes structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    m = (met *)malloc(sizeof(met));
    if (m == NULL) {
    	fprintf(stderr, "met structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    p = (params *)malloc(sizeof(params));
    if (p == NULL) {
    	fprintf(stderr, "params structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    s = (state *)malloc(sizeof(state));
    if (s == NULL) {
    	fprintf(stderr, "state structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    nr = (nrutil *)malloc(sizeof(nrutil));
    if (nr == NULL) {
        fprintf(stderr, "nrutil structure: Not allocated enough memory!\n");
        exit(EXIT_FAILURE);
    } 

    initialise_control(c);
    initialise_params(p);
    initialise_fluxes(f);
    initialise_state(s);
    initialise_nrutil(nr);
    
    clparser(argc, argv, c);
    /*
     * Read .ini parameter file and meterological data
     */
    error = parse_ini_file(c, p, s);
    if (error != 0) {
        prog_error("Error reading .INI file on line", __LINE__);
    }
    strcpy(c->git_code_ver, build_git_sha);
    if (c->PRINT_GIT) {
        fprintf(stderr, "\n%s\n", c->git_code_ver);
        exit(EXIT_FAILURE);
    }
    
    /* read met data */
    read_daily_met_data_simple(argv, c, m, p, f);
    
    /* model runs */
    if (c->spin_up) {
        spin_up_pools(c, f, m, p, s, nr);
    } else {
        run_sim(c, f, m, p, s, nr);
    }

    /* clean up */
    fclose(c->ofp);

    fclose(c->ifp);
    if (c->output_ascii == FALSE) {
        fclose(c->ofp_hdr);
    }

    free(c);
    free(m);
    free(p);
    free(s);
    free(f);

    exit(EXIT_SUCCESS);
}

void run_sim(control *c, fluxes *f,  met *m, 
             params *p, state *s, nrutil *nr){

    int    nyr, doy, window_size, i, dummy = 0;
    int    fire_found = FALSE;;

    double fdecay, rdecay, current_limitation, npitfac, year;

    /* Setup output file */

    if (c->print_options == DAILY && c->spin_up == FALSE) {
        /* Daily outputs */
        open_output_file(c, c->out_fname, &(c->ofp));

        if (c->output_ascii) {
            write_output_header(c, p, &(c->ofp));
        } else {
            open_output_file(c, c->out_fname_hdr, &(c->ofp_hdr));
            write_output_header(c, p, &(c->ofp_hdr));
        }
    } else if (c->print_options == END && c->spin_up == FALSE) {
        /* Final state + param file */
        open_output_file(c, c->out_param_fname, &(c->ofp));
    }

    day_end_calculations(c, p, s);

    s->lai = MAX(0.01, (p->sla * M2_AS_HA / KG_AS_TONNES /
                          p->cfracts * s->shoot));

    /* ====================== **
    **   Y E A R    L O O P   **
    ** ====================== */
    for (nyr = 0; nyr < p->num_years; nyr++) {

      fprintf(stderr, "nyr %d\n", nyr);
      
      unpack_met_data_simple(f, m, p);
      
      // fprintf(stderr, "ninflow %f\n", f->ninflow);
      
      fdecay = p->fdecay;
      rdecay = p->rdecay;
      
      // fprintf(stderr, "fdecay %f\n", p->fdecay);
      
      calculate_litterfall(c, f, p, s, &fdecay, &rdecay);
      
      calc_day_growth(c, f, m, nr, p, s,
                      fdecay, rdecay);
      
      calculate_csoil_flows(c, f, p, s, m->tsoil);
      calculate_nsoil_flows(c, f, p, s);
      
      if (c->pcycle == TRUE) {
        calculate_psoil_flows(c, f, p, s);
      }
      
      /* update stress SMA */
      current_limitation = calculate_growth_stress_limitation(p, s, c);
      
      /* Turn off all N calculations */
      if (c->ncycle == FALSE)
        reset_all_n_pools_and_fluxes(f, s);
      
      /* Turn off all P calculations */
      if (c->pcycle == FALSE)
        reset_all_p_pools_and_fluxes(f, s);
      
      /* calculate C:N ratios and increment annual flux sum */
      day_end_calculations(c, p, s);
      
      if (c->print_options == DAILY && c->spin_up == FALSE) {
        if(c->output_ascii)
          write_daily_outputs_ascii(c, f, s, year);
        else
          write_daily_outputs_binary(c, f, s, year);
      }
        
    }
    /* ========================= **
    **   E N D   O F   Y E A R   **
    ** ========================= */

    if (c->print_options == END && c->spin_up == FALSE) {
        write_final_state(c, p, s);
    }

    return;


}

void spin_up_pools(control *c, fluxes *f, met *m,
                   params *p, state *s, nrutil *nr){
    /* Spin up model plant & soil pools to equilibrium.

    - Examine sequences of 50 years and check if C pools are changing
      by more than 0.005 units per 1000 yrs.

    References:
    ----------
    Adapted from...
    * Murty, D and McMurtrie, R. E. (2000) Ecological Modelling, 134,
      185-205, specifically page 196.
    */
    double tol_c = 5E-03;
    double tol_n = 5E-04;
    double tol_p = 5E-04;
    double prev_plantc = 99999.9;
    double prev_soilc = 99999.9;
    double prev_plantn = 99999.9;
    double prev_soiln = 99999.9;
    double prev_plantp = 99999.9;
    double prev_soilp = 99999.9;
    int i, cntrl_flag;
    
    /* check for convergences in units of kg/m2 */
    double conv = TONNES_HA_2_KG_M2;

    /* Final state + param file */
    open_output_file(c, c->out_param_fname, &(c->ofp));

    fprintf(stderr, "Spinning up the model...\n");
    while (TRUE) {
        if (fabs((prev_plantc) - (s->plantc)) < tol_c &&
            fabs((prev_soilc) - (s->soilc)) < tol_c &&
            fabs((prev_plantn) - (s->plantn)) < tol_n &&
            fabs((prev_soiln) - (s->soiln)) < tol_n &&
            fabs((prev_plantp) - (s->plantp)) < tol_p &&
            fabs((prev_soilp) - (s->soilp)) < tol_p) {
            break;
        } else {
            prev_plantc = s->plantc;
            prev_soilc = s->soilc;
            prev_plantn = s->plantn;
            prev_soiln = s->soiln;
            prev_plantp = s->plantp;
            prev_soilp = s->soilp;

            /* 1700 years (17 yrs x 100 cycles) */
            for (i = 0; i < 100; i++) {
                run_sim(c, f, m, p, s, nr); /* run GDAY */
            }
            if (c->pcycle) {
                /* Have we reached a steady state? */
                fprintf(stderr,
                        "Spinup: Plant C - %f, Soil C - %f, Soil N - %f, Soil P - %f\n",
                        s->plantc, s->soilc, s->soiln, s->soilp);
            } else if (c->ncycle) {
              /* Have we reached a steady state? */
              fprintf(stderr,
                      "Spinup: Plant C - %f, Soil C - %f, Soil N - %f\n",
                      s->plantc, s->soilc, s->soiln);
            } else {
              /* Have we reached a steady state? */
              fprintf(stderr,
                      "Spinup: Plant C - %f, Soil C - %f\n",
                      s->plantc, s->soilc);
            }
        }
    }
    write_final_state(c, p, s);

    return;
}

void clparser(int argc, char **argv, control *c) {
    int i;

    for (i = 1; i < argc; i++) {
        if (*argv[i] == '-') {
            if (!strncasecmp(argv[i], "-p", 2)) {
			    strcpy(c->cfg_fname, argv[++i]);
            } else if (!strncasecmp(argv[i], "-s", 2)) {
                c->spin_up = TRUE;
            } else if (!strncasecmp(argv[i], "-ver", 4)) {
                c->PRINT_GIT = TRUE;
            } else if (!strncasecmp(argv[i], "-u", 2) ||
                       !strncasecmp(argv[i], "-h", 2)) {
                usage(argv);
                exit(EXIT_FAILURE);
            } else {
                fprintf(stderr, "%s: unknown argument on command line: %s\n",
                               argv[0], argv[i]);
                usage(argv);
                exit(EXIT_FAILURE);
            }
        }
    }
    return;
}


void usage(char **argv) {
    fprintf(stderr, "\n========\n");
    fprintf(stderr, " USAGE:\n");
    fprintf(stderr, "========\n");
    fprintf(stderr, "%s [options]\n", argv[0]);
    fprintf(stderr, "\n\nExpected input file is a .ini/.cfg style param file, passed with the -p flag .\n");
    fprintf(stderr, "\nThe options are:\n");
    fprintf(stderr, "\n++General options:\n" );
    fprintf(stderr, "[-ver          \t] Print the git hash tag.]\n");
    fprintf(stderr, "[-p       fname\t] Location of parameter file (.ini/.cfg).]\n");
    fprintf(stderr, "[-s            \t] Spin-up GDAY, when it the model is finished it will print the final state to the param file.]\n");
    fprintf(stderr, "\n++Print this message:\n" );
    fprintf(stderr, "[-u/-h         \t] usage/help]\n");

    return;
}

void reset_all_n_pools_and_fluxes(fluxes *f, state *s) {
    /*
        If the N-Cycle is turned off the way I am implementing this is to
        do all the calculations and then reset everything at the end. This is
        a waste of resources but saves on multiple IF statements.
    */

    /*
    ** State
    */
    s->shootn = 0.0;
    s->rootn = 0.0;
    s->branchn = 0.0;
    s->stemnimm = 0.0;
    s->stemnmob = 0.0;
    s->structsurfn = 0.0;
    s->metabsurfn = 0.0;
    s->structsoiln = 0.0;
    s->metabsoiln = 0.0;
    s->activesoiln = 0.0;
    s->slowsoiln = 0.0;
    s->passivesoiln = 0.0;
    s->inorgn = 0.0;
    s->stemn = 0.0;
    s->stemnimm = 0.0;
    s->stemnmob = 0.0;

    /*
    ** Fluxes
    */
    f->nuptake = 0.0;
    f->nloss = 0.0;
    f->npassive = 0.0;
    f->ngross = 0.0;
    f->nimmob = 0.0;
    f->nlittrelease = 0.0;
    f->nmineralisation = 0.0;
    f->npleaf = 0.0;
    f->nproot = 0.0;
    f->npbranch = 0.0;
    f->npstemimm = 0.0;
    f->npstemmob = 0.0;
    f->deadleafn = 0.0;
    f->deadrootn = 0.0;
    f->deadbranchn = 0.0;
    f->deadstemn = 0.0;
    f->leafretransn = 0.0;
    f->n_surf_struct_litter = 0.0;
    f->n_surf_metab_litter = 0.0;
    f->n_soil_struct_litter = 0.0;
    f->n_soil_metab_litter = 0.0;
    f->n_surf_struct_to_slow = 0.0;
    f->n_soil_struct_to_slow = 0.0;
    f->n_surf_struct_to_active = 0.0;
    f->n_soil_struct_to_active = 0.0;
    f->n_surf_metab_to_active = 0.0;
    f->n_surf_metab_to_active = 0.0;
    f->n_active_to_slow = 0.0;
    f->n_active_to_passive = 0.0;
    f->n_slow_to_active = 0.0;
    f->n_slow_to_passive = 0.0;
    f->n_passive_to_active = 0.0;

    return;
}

void reset_all_p_pools_and_fluxes(fluxes *f, state *s) {
    /*
        If the P-Cycle is turned off the way I am implementing this is to
        do all the calculations and then reset everything at the end. This is
        a waste of resources but saves on multiple IF statements.
    */

    /*
    ** State
    */
    s->shootp = 0.0;
    s->rootp = 0.0;
    s->branchp = 0.0;
    s->stempimm = 0.0;
    s->stempmob = 0.0;
    s->structsurfp = 0.0;
    s->metabsurfp = 0.0;
    s->structsoilp = 0.0;
    s->metabsoilp = 0.0;
    s->activesoilp = 0.0;
    s->slowsoilp = 0.0;
    s->passivesoilp = 0.0;
    s->inorgp = 0.0;
    s->inorgavlp = 0.0;
    s->inorgssorbp = 0.0;
    s->inorgoccp = 0.0;
    s->inorgparp = 0.0;
    s->stemp = 0.0;
    s->stempimm = 0.0;
    s->stempmob = 0.0;

    /*
    ** Fluxes
    */
    f->puptake = 0.0;
    f->ploss = 0.0;
    f->ppassive = 0.0;
    f->pgross = 0.0;
    f->pimmob = 0.0;
    f->plittrelease = 0.0;
    f->pmineralisation = 0.0;
    f->ppleaf = 0.0;
    f->pproot = 0.0;
    f->ppbranch = 0.0;
    f->ppstemimm = 0.0;
    f->ppstemmob = 0.0;
    f->deadleafp = 0.0;
    f->deadrootp = 0.0;
    f->deadbranchp = 0.0;
    f->deadstemp = 0.0;
    f->leafretransp = 0.0;
    f->p_surf_struct_litter = 0.0;
    f->p_surf_metab_litter = 0.0;
    f->p_soil_struct_litter = 0.0;
    f->p_soil_metab_litter = 0.0;
    f->p_surf_struct_to_slow = 0.0;
    f->p_soil_struct_to_slow = 0.0;
    f->p_surf_struct_to_active = 0.0;
    f->p_soil_struct_to_active = 0.0;
    f->p_surf_metab_to_active = 0.0;
    f->p_surf_metab_to_active = 0.0;
    f->p_active_to_slow = 0.0;
    f->p_active_to_passive = 0.0;
    f->p_slow_to_active = 0.0;
    f->p_slow_to_passive = 0.0;
    f->p_slow_biochemical = 0.0;
    f->p_passive_to_active = 0.0;
    f->p_avl_to_ssorb = 0.0;
    f->p_ssorb_to_avl = 0.0;
    f->p_ssorb_to_occ = 0.0;
    f->p_par_to_avl = 0.0;
    f->p_atm_dep = 0.0;

    return;
}

void day_end_calculations(control *c, params *p, state *s) {
    /* Calculate derived values from state variables.

    Parameters:
    -----------
    day : integer
        day of simulation

    INIT : logical
        logical defining whether it is the first day of the simulation
    */
    
    /* update N:C and P:C of plant pool */
    if (float_eq(s->shoot, 0.0)) {
        s->shootnc = 0.0;
        s->shootpc = 0.0;
    } else {
        s->shootnc = s->shootn / s->shoot;
        s->shootpc = s->shootp / s->shoot;
    }

    /* Explicitly set the shoot N:C */
    if (c->ncycle == FALSE)
        s->shootnc = p->prescribed_leaf_NC;

    if (c->pcycle == FALSE)
        s->shootpc = p->prescribed_leaf_PC;

    if (float_eq(s->root, 0.0)) {
        s->rootnc = 0.0;
        s->rootpc = 0.0;
    } else {
        s->rootnc = MAX(0.0, s->rootn / s->root);
        s->rootpc = MAX(0.0, s->rootp / s->root);
    }

    /* total plant, soil & litter nitrogen */
    s->soiln = s->inorgn + s->activesoiln + s->slowsoiln + s->passivesoiln;
    s->litternag = s->structsurfn + s->metabsurfn;
    s->litternbg = s->structsoiln + s->metabsoiln;
    s->littern = s->litternag + s->litternbg;
    s->plantn = s->shootn + s->rootn + s->branchn + s->stemn;
    s->totaln = s->plantn + s->littern + s->soiln;

    /* total plant, soil & litter phosphorus */
    s->inorgp = s->inorgavlp + s->inorgssorbp + s->inorgoccp + s->inorgparp;
    s->soilp = s->inorgp + s->activesoilp + s->slowsoilp + s->passivesoilp;
    s->litterpag = s->structsurfp + s->metabsurfp;
    s->litterpbg = s->structsoilp + s->metabsoilp;
    s->litterp = s->litterpag + s->litterpbg;
    s->plantp = s->shootp + s->rootp + s->branchp + s->stemp;
    s->totalp = s->plantp + s->litterp + s->soilp;

    /* total plant, soil, litter and system carbon */
    s->soilc = s->activesoil + s->slowsoil + s->passivesoil;
    s->littercag = s->structsurf + s->metabsurf;
    s->littercbg = s->structsoil + s->metabsoil;
    s->litterc = s->littercag + s->littercbg;
    s->plantc = s->root + s->shoot + s->stem + s->branch;
    s->totalc = s->soilc + s->litterc + s->plantc;

    /* optional constant passive pool */
    if (c->passiveconst) {
        s->passivesoil = p->passivesoilz;
        s->passivesoiln = p->passivesoilnz;
        s->passivesoilp = p->passivesoilpz;
    }
    
    return;
}

void unpack_met_data_simple(fluxes *f, met *m, params *p) {
  
  /* unpack met forcing */
  m->Ca = p->co2_in;
  //m->par = p->I0/365.25;  
  m->par = p->I0;    
  
  m->ndep = p->ndep_in;
  m->nfix = p->nfix_in;
  m->pdep = p->pdep_in;
  m->tsoil = p->tsoil_in;
  
  //f->ninflow = (m->ndep + m->nfix) / 365.25;
  //f->p_atm_dep = m->pdep / 365.25;
  
  f->ninflow = (m->ndep + m->nfix);
  f->p_atm_dep = m->pdep;
  
  return;
}
