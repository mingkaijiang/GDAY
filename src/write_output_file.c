/* ============================================================================
* Print output file (ascii/binary)
*
*
*
* NOTES:
*   Currently I have only implemented the ascii version
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   25.02.2015
*
* =========================================================================== */
#include "write_output_file.h"


void open_output_file(control *c, char *fname, FILE **fp) {
    *fp = fopen(fname, "w");
    if (*fp == NULL)
        prog_error("Error opening output file for write on line", __LINE__);
}

void write_output_header(control *c, params *p, FILE **fp) {
    /*
        Write annual state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    int ncols = 86;  /* change with number of variables ? total count below is 93 */
    int nrows = p->num_years;

    /* Git version */
    fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "year,");

    /*
    ** STATE
    */
    /* plant */
    fprintf(*fp, "shoot,lai,branch,stem,root,");
    fprintf(*fp, "shootn,branchn,stemn,rootn,");
    fprintf(*fp, "shootp,branchp,stemp,rootp,");

    /* belowground */
    fprintf(*fp, "soilc,soiln,soilp,inorgn,");
    fprintf(*fp, "inorgp,inorgavlp,inorgssorbp,inorgoccp,inorgparp,");
    fprintf(*fp, "litterc,littercag,littercbg,litternag,litternbg,");
    fprintf(*fp, "litterpag,litterpbg,");
    fprintf(*fp, "activesoil,slowsoil,passivesoil,");
    fprintf(*fp, "activesoiln,slowsoiln,passivesoiln,activesoilp,slowsoilp,passivesoilp,");

    /*
    ** FLUXES
    */

    /* litter */
    fprintf(*fp, "deadleaves,deadbranch,deadstems,deadroots,");
    fprintf(*fp, "deadleafn,deadbranchn,deadstemn,deadrootn,");
    fprintf(*fp, "deadleafp,deadbranchp,deadstemp,deadrootp,");


    /* C fluxes */
    fprintf(*fp, "nep,gpp,npp,hetero_resp,auto_resp,apar,");

    /* C, N and P growth */
    fprintf(*fp, "cpleaf,cpbranch,cpstem,cproot,");
    fprintf(*fp, "npleaf,npbranch,npstemimm,npstemmob,nproot,");
    fprintf(*fp, "ppleaf,ppbranch,ppstemimm,ppstemmob,pproot,");


    /* N stuff */
    fprintf(*fp, "nuptake,ngross,nmineralisation,nloss,");

    /* P stuff */
    fprintf(*fp, "puptake,pgross,pmineralisation,ploss,");

    /* traceability stuff */
    fprintf(*fp, "tfac_soil_decomp,c_into_active,c_into_slow,");
    fprintf(*fp, "c_into_passive,active_to_slow,active_to_passive,");
    fprintf(*fp, "slow_to_active,slow_to_passive,passive_to_active,");
    fprintf(*fp, "co2_rel_from_surf_struct_litter,");
    fprintf(*fp, "co2_rel_from_soil_struct_litter,");
    fprintf(*fp, "co2_rel_from_surf_metab_litter,");
    fprintf(*fp, "co2_rel_from_soil_metab_litter,");
    fprintf(*fp, "co2_rel_from_active_pool,");
    fprintf(*fp, "co2_rel_from_slow_pool,");
    fprintf(*fp, "co2_rel_from_passive_pool,");

    /* Misc */
    fprintf(*fp, "leafretransn,");
    fprintf(*fp, "leafretransp\n");


    if (c->output_ascii == FALSE) {
        fprintf(*fp, "nrows=%d\n", nrows);
        fprintf(*fp, "ncols=%d\n", ncols);
    }
    return;
}

void write_annual_outputs_ascii(control *c, fluxes *f, state *s, int year) {
    /*
        Write annual state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */


    /* time stuff */
    fprintf(c->ofp, "%.10f,", (double)year);

    /*
    ** STATE

    */

    /* plant */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->shoot,s->lai,s->branch,s->stem,s->root);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    s->shootn,s->branchn,s->stemn,s->rootn);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    s->shootp,s->branchp,s->stemp,s->rootp);

    /* belowground */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    s->soilc,s->soiln,s->soilp,s->inorgn);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->inorgp,s->inorgavlp,s->inorgssorbp,
                    s->inorgoccp,s->inorgparp);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->litterc,s->littercag,s->littercbg,s->litternag,s->litternbg);

    fprintf(c->ofp, "%.10f,%.10f,",
                    s->litterpag,s->litterpbg);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    s->activesoil,s->slowsoil,s->passivesoil);

    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->activesoiln,s->slowsoiln,s->passivesoiln,
                    s->activesoilp,s->slowsoilp,s->passivesoilp);
    /*
    ** FLUXES
    */

    /* litter */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->deadleaves,f->deadbranch,f->deadstems,f->deadroots);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->deadleafn,f->deadbranchn,f->deadstemn,f->deadrootn);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->deadleafp,f->deadbranchp,f->deadstemp,f->deadrootp);

    /* C fluxes */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->nep,f->gpp,f->npp,f->hetero_resp,f->auto_resp,
                    f->apar);

    /* C N and P growth */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->cpleaf,f->cpbranch,f->cpstem,f->cproot);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->npleaf,f->npbranch,f->npstemimm,f->npstemmob,f->nproot);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->ppleaf,f->ppbranch,f->ppstemimm,f->ppstemmob,f->pproot);

    /* N stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->nuptake,f->ngross,f->nmineralisation,f->nloss);

    /* P stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->puptake,f->pgross,f->pmineralisation,f->ploss);


    /* traceability stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->tfac_soil_decomp,f->c_into_active,f->c_into_slow);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->c_into_passive,f->active_to_slow,f->active_to_passive);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->slow_to_active,f->slow_to_passive,f->passive_to_active);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_active_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_slow_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_passive_pool);

    /* Misc */
    fprintf(c->ofp, "%.10f,", f->leafretransn);
    fprintf(c->ofp, "%.10f\n", f->leafretransp);


    return;
}

void write_annual_outputs_binary(control *c, fluxes *f, state *s, int year) {
    /*
        Write annual state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    double temp;

    /* time stuff */
    temp = (double)year;
    fwrite(&temp, sizeof(double), 1, c->ofp);


    /* plant */
    fwrite(&(s->shoot), sizeof(double), 1, c->ofp);
    fwrite(&(s->lai), sizeof(double), 1, c->ofp);
    fwrite(&(s->branch), sizeof(double), 1, c->ofp);
    fwrite(&(s->stem), sizeof(double), 1, c->ofp);
    fwrite(&(s->root), sizeof(double), 1, c->ofp);

    /* C fluxes */
    fwrite(&(f->npp), sizeof(double), 1, c->ofp);

    return;
}


int write_final_state(control *c, params *p, state *s)
{
    /*
    Write the final state to the input param file so we can easily restart
    the model. This function copies the input param file with the exception
    of anything in the git hash and the state which it replaces with the updated
    stuff.

    */

    char line[STRING_LENGTH];
    char saved_line[STRING_LENGTH];
    char section[STRING_LENGTH] = "";
    char prev_name[STRING_LENGTH] = "";
    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;
    int line_number = 0;
    int match = FALSE;

    while (fgets(line, sizeof(line), c->ifp) != NULL) {
        strcpy(saved_line, line);
        line_number++;
        start = lskip(rstrip(line));
        if (*start == ';' || *start == '#') {
            /* Per Python ConfigParser, allow '#' comments at start of line */
        }
        else if (*start == '[') {
            /* A "[section]" line */
            end = find_char_or_comment(start + 1, ']');
            if (*end == ']') {
                *end = '\0';
                strncpy0(section, start + 1, sizeof(section));
                *prev_name = '\0';

            }
            else if (!error) {
                /* No ']' found on section line */
                error = line_number;

            }
        }
        else if (*start && *start != ';') {
            /* Not a comment, must be a name[=:]value pair */
            end = find_char_or_comment(start, '=');
            if (*end != '=') {
                end = find_char_or_comment(start, ':');
            }
            if (*end == '=' || *end == ':') {
                *end = '\0';
                name = rstrip(start);
                value = lskip(end + 1);
                end = find_char_or_comment(value, '\0');
                if (*end == ';')
                    *end = '\0';
                rstrip(value);

                /* Valid name[=:]value pair found, call handler */
                strncpy0(prev_name, name, sizeof(prev_name));

                if (!ohandler(section, name, value, c, p, s, &match) && !error)
                    error = line_number;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = line_number;
                break;
            }
        }
        if (match == FALSE)
            fprintf(c->ofp, "%s", saved_line);
        else
            match = FALSE; /* reset match flag */
    }
    return error;

}


int ohandler(char *section, char *name, char *value, control *c, params *p,
             state *s, int *match)
{
    /*
    Search for matches of the git and state values and where found write the
    current state values to the output parameter file.

    - also added previous ncd as this potential can be changed internally
    */

    #define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    /*
    ** GIT
    */
    if (MATCH("git", "git_hash")) {
        fprintf(c->ofp, "git_hash = %s\n", c->git_code_ver);
        *match = TRUE;
    }

    /*
    ** STATE
    */

    if (MATCH("state", "activesoil")) {
        fprintf(c->ofp, "activesoil = %.10f\n", s->activesoil);
        *match = TRUE;
    } else if (MATCH("state", "activesoiln")) {
        fprintf(c->ofp, "activesoiln = %.10f\n", s->activesoiln);
        *match = TRUE;
    } else if (MATCH("state", "activesoilp")) {
        fprintf(c->ofp, "activesoilp = %.10f\n", s->activesoilp);
        *match = TRUE;
    } else if (MATCH("state", "branch")) {
        fprintf(c->ofp, "branch = %.10f\n", s->branch);
        *match = TRUE;
    } else if (MATCH("state", "branchn")) {
        fprintf(c->ofp, "branchn = %.10f\n", s->branchn);
        *match = TRUE;
    } else if (MATCH("state", "branchp")) {
        fprintf(c->ofp, "branchp = %.10f\n", s->branchp);
        *match = TRUE;
    } else if (MATCH("state", "inorgn")) {
        fprintf(c->ofp, "inorgn = %.10f\n", s->inorgn);
        *match = TRUE;
    } else if (MATCH("state", "inorgp")) {
        fprintf(c->ofp, "inorgp = %.10f\n", s->inorgp);
        *match = TRUE;
    } else if (MATCH("state", "inorgavlp")) {
        fprintf(c->ofp, "inorgavlp = %.10f\n", s->inorgavlp);
        *match = TRUE;
    } else if (MATCH("state", "inorgssorbp")) {
        fprintf(c->ofp, "inorgssorbp = %.10f\n", s->inorgssorbp);
        *match = TRUE;
    } else if (MATCH("state", "inorgoccp")) {
        fprintf(c->ofp, "inorgoccp = %.10f\n", s->inorgoccp);
        *match = TRUE;
    } else if (MATCH("state", "inorgparp")) {
        fprintf(c->ofp, "inorgparp = %.10f\n", s->inorgparp);
        *match = TRUE;
    } else if (MATCH("state", "lai")) {
        fprintf(c->ofp, "lai = %.10f\n", s->lai);
        *match = TRUE;
    } else if (MATCH("state", "metabsoil")) {
        fprintf(c->ofp, "metabsoil = %.10f\n", s->metabsoil);
        *match = TRUE;
    } else if (MATCH("state", "metabsoiln")) {
        fprintf(c->ofp, "metabsoiln = %.10f\n", s->metabsoiln);
        *match = TRUE;
    } else if (MATCH("state", "metabsoilp")) {
        fprintf(c->ofp, "metabsoilp = %.10f\n", s->metabsoilp);
        *match = TRUE;
    } else if (MATCH("state", "metabsurf")) {
        fprintf(c->ofp, "metabsurf = %.10f\n", s->metabsurf);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfn")) {
        fprintf(c->ofp, "metabsurfn = %.10f\n", s->metabsurfn);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfp")) {
        fprintf(c->ofp, "metabsurfp = %.10f\n", s->metabsurfp);
        *match = TRUE;
    } else if (MATCH("state", "passivesoil")) {
        fprintf(c->ofp, "passivesoil = %.10f\n", s->passivesoil);
        *match = TRUE;
    } else if (MATCH("state", "passivesoiln")) {
        fprintf(c->ofp, "passivesoiln = %.10f\n", s->passivesoiln);
        *match = TRUE;
    } else if (MATCH("state", "passivesoilp")) {
        fprintf(c->ofp, "passivesoilp = %.10f\n", s->passivesoilp);
        *match = TRUE;
    } else if (MATCH("state", "root")) {
        fprintf(c->ofp, "root = %.10f\n", s->root);
        *match = TRUE;
    } else if (MATCH("state", "rootn")) {
        fprintf(c->ofp, "rootn = %.10f\n", s->rootn);
        *match = TRUE;
    } else if (MATCH("state", "rootp")) {
        fprintf(c->ofp, "rootp = %.10f\n", s->rootp);
        *match = TRUE;
    } else if (MATCH("state", "shoot")) {
        fprintf(c->ofp, "shoot = %.10f\n", s->shoot);
        *match = TRUE;
    } else if (MATCH("state", "shootn")) {
        fprintf(c->ofp, "shootn = %.10f\n", s->shootn);
        *match = TRUE;
    } else if (MATCH("state", "shootp")) {
        fprintf(c->ofp, "shootp = %.10f\n", s->shootp);
        *match = TRUE;
    } else if (MATCH("state", "sla")) {
        fprintf(c->ofp, "sla = %.10f\n", s->sla);
        *match = TRUE;
    } else if (MATCH("state", "slowsoil")) {
        fprintf(c->ofp, "slowsoil = %.10f\n", s->slowsoil);
        *match = TRUE;
    } else if (MATCH("state", "slowsoiln")) {
        fprintf(c->ofp, "slowsoiln = %.10f\n", s->slowsoiln);
        *match = TRUE;
    } else if (MATCH("state", "slowsoilp")) {
        fprintf(c->ofp, "slowsoilp = %.10f\n", s->slowsoilp);
        *match = TRUE;
    } else if (MATCH("state", "stem")) {
        fprintf(c->ofp, "stem = %.10f\n", s->stem);
        *match = TRUE;
    } else if (MATCH("state", "stemn")) {
        fprintf(c->ofp, "stemn = %.10f\n", s->stemn);
        *match = TRUE;
    } else if (MATCH("state", "stemnimm")) {
        fprintf(c->ofp, "stemnimm = %.10f\n", s->stemnimm);
        *match = TRUE;
    } else if (MATCH("state", "stemnmob")) {
        fprintf(c->ofp, "stemnmob = %.10f\n", s->stemnmob);
        *match = TRUE;
    } else if (MATCH("state", "stemp")) {
        fprintf(c->ofp, "stemp = %.10f\n", s->stemp);
        *match = TRUE;
    } else if (MATCH("state", "stempimm")) {
        fprintf(c->ofp, "stempimm = %.10f\n", s->stempimm);
        *match = TRUE;
    } else if (MATCH("state", "stempmob")) {
        fprintf(c->ofp, "stempmob = %.10f\n", s->stempmob);
        *match = TRUE;
    } else if (MATCH("state", "structsoil")) {
        fprintf(c->ofp, "structsoil = %.10f\n", s->structsoil);
        *match = TRUE;
    } else if (MATCH("state", "structsoiln")) {
        fprintf(c->ofp, "structsoiln = %.10f\n", s->structsoiln);
        *match = TRUE;
    } else if (MATCH("state", "structsoilp")) {
        fprintf(c->ofp, "structsoilp = %.10f\n", s->structsoilp);
        *match = TRUE;
    } else if (MATCH("state", "structsurf")) {
        fprintf(c->ofp, "structsurf = %.10f\n", s->structsurf);
        *match = TRUE;
    } else if (MATCH("state", "structsurfn")) {
        fprintf(c->ofp, "structsurfn = %.10f\n", s->structsurfn);
        *match = TRUE;
    } else if (MATCH("state", "structsurfp")) {
        fprintf(c->ofp, "structsurfp = %.10f\n", s->structsurfp);
        *match = TRUE;
    }

    return (1);
}
