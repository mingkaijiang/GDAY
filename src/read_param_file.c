/* Read the .ini file into the various structures. */

#include "read_param_file.h"

int parse_ini_file(control *c, params *p, state *s) {
    /*

    Loop through the file which is passed on the standard in, and break
    the file up into the relevant sections...

    */

    char line[STRING_LENGTH];
    char section[STRING_LENGTH] = "";
    char prev_name[STRING_LENGTH] = "";
    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;
    int line_number = 0;

    if ((c->ifp = fopen(c->cfg_fname, "r")) == NULL){
        prog_error("Error opening output file for write on line", __LINE__);
    }

    while (fgets(line, sizeof(line), c->ifp) != NULL) {
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

                if (!handler(section, name, value, c, p, s) && !error)
                    error = line_number;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = line_number;
                break;
            }
        }
    }

    if (c->print_options == END) {
        /* we need to re-read this file to dump the final state */
        rewind(c->ifp);
    }

    return error;


}



int handler(char *section, char *name, char *value, control *c,
            params *p, state *s)
{
    /*

    Assigns the values from the .INI file straight into the various
    structures

    */
    char *temp = value;

    #define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    /*
    ** GIT
    */
    if (MATCH("git", "git_hash")) {
        strcpy(c->git_hash, temp);
    }

    /*
    ** FILES
    */
    if (MATCH("files", "cfg_fname")) {
        strcpy(c->cfg_fname, temp);
    } else if (MATCH("files", "met_fname")) {
      strcpy(c->met_fname, temp);
    } else if (MATCH("files", "out_fname")) {
        strcpy(c->out_fname, temp);
    } else if (MATCH("files", "out_param_fname")) {
        strcpy(c->out_param_fname, temp);
    }

    /*
    ** CONTROL
    */
    if (MATCH("control", "fixed_stem_nc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixed_stem_nc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->fixed_stem_nc = TRUE;
        else {
            fprintf(stderr, "Unknown fixed_stem_nc option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "fixed_stem_pc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixed_stem_pc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
                 strcmp(temp, "TRUE") == 0 ||
                 strcmp(temp, "true") == 0)
             c->fixed_stem_pc = TRUE;
      else {
          fprintf(stderr, "Unknown fixed_stem_pc option: %s\n", temp);
          exit(EXIT_FAILURE);
      }
    } else if (MATCH("control", "diagnosis")) {
      if (strcmp(temp, "False") == 0 ||
          strcmp(temp, "FALSE") == 0 ||
          strcmp(temp, "false") == 0)
        c->diagnosis = FALSE;
      else if (strcmp(temp, "True") == 0 ||
               strcmp(temp, "TRUE") == 0 ||
               strcmp(temp, "true") == 0)
        c->diagnosis = TRUE;
      else {
        fprintf(stderr, "Unknown diagnosis option: %s\n", temp);
        exit(EXIT_FAILURE);
      }
    } else if (MATCH("control", "fixleafnc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixleafnc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->fixleafnc = TRUE;
        else {
            fprintf(stderr, "Unknown fixleafnc option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "fixleafpc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixleafpc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
               strcmp(temp, "TRUE") == 0 ||
               strcmp(temp, "true") == 0)
               c->fixleafpc = TRUE;
        else {
            fprintf(stderr, "Unknown fixleafpc option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "ncycle")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->ncycle = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->ncycle = TRUE;
        else {
            fprintf(stderr, "Unknown ncycle option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "pcycle")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->pcycle = FALSE;
        else if (strcmp(temp, "True") == 0 ||
                 strcmp(temp, "TRUE") == 0 ||
                 strcmp(temp, "true") == 0)
            c->pcycle = TRUE;
        else {
            fprintf(stderr, "Unknown pcycle option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "nuptake_model")) {
        c->nuptake_model = atoi(value);
    } else if (MATCH("control", "puptake_model")) {
        c->puptake_model = atoi(value);
    } else if (MATCH("control", "print_options")) {
        if (strcmp(temp, "Annual") == 0 ||
            strcmp(temp, "ANNUAL") == 0 ||
            strcmp(temp, "annual") == 0)
            c->print_options = ANNUAL;
        else if (strcmp(temp, "End") == 0 ||
            strcmp(temp, "END") == 0 ||
            strcmp(temp, "end") == 0)
            c->print_options = END;
        else {
            fprintf(stderr, "Unknown print option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "use_eff_nc")) {
        c->use_eff_nc = atoi(value);
    }


    /*
    ** State
    */
    if (MATCH("state", "activesoil")) {
        s->activesoil = atof(value);
    } else if (MATCH("state", "activesoiln")) {
        s->activesoiln = atof(value);
    } else if (MATCH("state", "activesoilp")) {
        s->activesoilp = atof(value);
    } else if (MATCH("state", "branch")) {
        s->branch = atof(value);
    } else if (MATCH("state", "branchn")) {
        s->branchn = atof(value);
    } else if (MATCH("state", "branchp")) {
        s->branchp = atof(value);
    } else if (MATCH("state", "inorgn")) {
        s->inorgn = atof(value);
    } else if (MATCH("state", "inorgp")) {
        s->inorgp = atof(value);
    } else if (MATCH("state", "inorgavlp")) {
        s->inorgavlp = atof(value);
    } else if (MATCH("state", "inorgssorbp")) {
        s->inorgssorbp = atof(value);
    } else if (MATCH("state", "inorgoccp")) {
        s->inorgoccp = atof(value);
    } else if (MATCH("state", "inorgparp")) {
        s->inorgparp = atof(value);
    } else if (MATCH("state", "lai")) {
        s->lai = atof(value);
    } else if (MATCH("state", "metabsoil")) {
        s->metabsoil = atof(value);
    } else if (MATCH("state", "metabsoiln")) {
        s->metabsoiln = atof(value);
    } else if (MATCH("state", "metabsoilp")) {
        s->metabsoilp = atof(value);
    } else if (MATCH("state", "metabsurf")) {
        s->metabsurf = atof(value);
    } else if (MATCH("state", "metabsurfn")) {
        s->metabsurfn = atof(value);
    } else if (MATCH("state", "metabsurfp")) {
        s->metabsurfp = atof(value);
    } else if (MATCH("state", "passivesoil")) {
        s->passivesoil = atof(value);
    } else if (MATCH("state", "passivesoiln")) {
        s->passivesoiln = atof(value);
    } else if (MATCH("state", "passivesoilp")) {
        s->passivesoilp = atof(value);
    } else if (MATCH("state", "root")) {
        s->root = atof(value);
    } else if (MATCH("state", "rootn")) {
        s->rootn = atof(value);
    } else if (MATCH("state", "rootp")) {
        s->rootp = atof(value);
    } else if (MATCH("state", "shoot")) {
        s->shoot = atof(value);
    } else if (MATCH("state", "shootn")) {
        s->shootn = atof(value);
    } else if (MATCH("state", "shootp")) {
        s->shootp = atof(value);
    } else if (MATCH("state", "sla")) {
        s->sla = atof(value);
    } else if (MATCH("state", "slowsoil")) {
        s->slowsoil = atof(value);
    } else if (MATCH("state", "slowsoiln")) {
        s->slowsoiln = atof(value);
    } else if (MATCH("state", "slowsoilp")) {
        s->slowsoilp = atof(value);
    } else if (MATCH("state", "stem")) {
        s->stem = atof(value);
    } else if (MATCH("state", "stemn")) {
        s->stemn = atof(value);
    } else if (MATCH("state", "stemnimm")) {
        s->stemnimm = atof(value);
    } else if (MATCH("state", "stemnmob")) {
        s->stemnmob = atof(value);
    } else if (MATCH("state", "stemp")) {
        s->stemp = atof(value);
    } else if (MATCH("state", "stempimm")) {
        s->stempimm = atof(value);
    } else if (MATCH("state", "stempmob")) {
        s->stempmob = atof(value);
    } else if (MATCH("state", "structsoil")) {
        s->structsoil = atof(value);
    } else if (MATCH("state", "structsoiln")) {
        s->structsoiln = atof(value);
    } else if (MATCH("state", "structsoilp")) {
        s->structsoilp = atof(value);
    } else if (MATCH("state", "structsurf")) {
        s->structsurf = atof(value);
    } else if (MATCH("state", "structsurfn")) {
        s->structsurfn = atof(value);
    } else if (MATCH("state", "structsurfp")) {
        s->structsurfp = atof(value);
    }

    /* Params */
    if (MATCH("params", "actncmax")) {
        p->actncmax = atof(value);
    } else if (MATCH("params", "actncmin")) {
        p->actncmin = atof(value);
    } else if (MATCH("params", "actpcmax")) {
        p->actpcmax = atof(value);
    } else if (MATCH("params", "actpcmin")) {
        p->actpcmin = atof(value);
    } else if (MATCH("params", "bdecay")) {
        p->bdecay = atof(value);
    } else if (MATCH("params", "bretrans")) {
      p->bretrans = atof(value);
    } else if (MATCH("params", "c_alloc_b")) {
        p->c_alloc_b = atof(value);
    } else if (MATCH("params", "c_alloc_f")) {
        p->c_alloc_f = atof(value);
    } else if (MATCH("params", "c_alloc_r")) {
        p->c_alloc_r = atof(value);
    } else if (MATCH("params", "c_alloc_s")) {
        p->c_alloc_s = atof(value);
    } else if (MATCH("params", "cfracts")) {
        p->cfracts = atof(value);
    } else if (MATCH("params", "cue")) {
        p->cue = atof(value);
    } else if (MATCH("params", "co2_in")) {
      p->co2_in = atof(value);
    } else if (MATCH("params", "fdecay")) {
        p->fdecay = atof(value);
    } else if (MATCH("params", "finesoil")) {
        p->finesoil = atof(value);
    } else if (MATCH("params", "fretransn")) {
        p->fretransn = atof(value);
    } else if (MATCH("params", "fretransp")) {
        p->fretransp = atof(value);
    } else if (MATCH("params", "I0")) {
      p->I0 = atof(value);
    } else if (MATCH("params", "k1")) {
      p->k1 = atof(value);
    } else if (MATCH("params", "k2")) {
      p->k2 = atof(value);
    } else if (MATCH("params", "k3")) {
      p->k3 = atof(value);
    } else if (MATCH("params", "kdec1")) {
        p->kdec1 = atof(value);
    } else if (MATCH("params", "kdec2")) {
        p->kdec2 = atof(value);
    } else if (MATCH("params", "kdec3")) {
        p->kdec3 = atof(value);
    } else if (MATCH("params", "kdec4")) {
        p->kdec4 = atof(value);
    } else if (MATCH("params", "kdec5")) {
        p->kdec5 = atof(value);
    } else if (MATCH("params", "kdec6")) {
        p->kdec6 = atof(value);
    } else if (MATCH("params", "kdec7")) {
        p->kdec7 = atof(value);
    } else if (MATCH("params", "kr")) {
        p->kr = atof(value);
    } else if (MATCH("params", "krp")) {
        p->krp = atof(value);
    } else if (MATCH("params", "ligroot")) {
        p->ligroot = atof(value);
    } else if (MATCH("params", "ligshoot")) {
        p->ligshoot = atof(value);
    } else if (MATCH("params", "liteffnc")) {
        p->liteffnc = atof(value);
    } else if (MATCH("params", "ncbnewz")) {
        p->ncbnewz = atof(value);
    } else if (MATCH("params", "ndep_in")) {
      p->ndep_in = atof(value);
    } else if (MATCH("params", "nfix_in")) {
      p->nfix_in = atof(value);
    } else if (MATCH("params", "lue0")) {
      p->lue0 = atof(value);
    } else if (MATCH("params", "pcbnewz")) {
        p->pcbnewz = atof(value);
    } else if (MATCH("params", "ncmaxf")) {
        p->ncmaxf = atof(value);
    } else if (MATCH("params", "nref")) {
      p->nref = atof(value);
    } else if (MATCH("params", "ncmaxr")) {
        p->ncmaxr = atof(value);
    } else if (MATCH("params", "pcmaxf")) {
        p->pcmaxf = atof(value);
    } else if (MATCH("params", "pcmaxr")) {
        p->pcmaxr = atof(value);
    } else if (MATCH("params", "ncrfac")) {
        p->ncrfac = atof(value);
    } else if (MATCH("params", "pcrfac")) {
        p->pcrfac = atof(value);
    } else if (MATCH("params", "pdep_in")) {
      p->pdep_in = atof(value);
    } else if (MATCH("params", "ncwimmz")) {
        p->ncwimmz = atof(value);
    } else if (MATCH("params", "ncwnewz")) {
        p->ncwnewz = atof(value);
    } else if (MATCH("params", "pcwimmz")) {
        p->pcwimmz = atof(value);
    } else if (MATCH("params", "pcwnewz")) {
        p->pcwnewz = atof(value);
    } else if (MATCH("params", "nf_min")) {
        p->nf_min = atof(value);
    } else if (MATCH("params", "pf_min")) {
        p->pf_min = atof(value);
    } else if (MATCH("params", "nmin0")) {
        p->nmin0 = atof(value);
    } else if (MATCH("params", "nmincrit")) {
        p->nmincrit = atof(value);
    } else if (MATCH("params", "nuptakez")) {
        p->nuptakez = atof(value);
    } else if (MATCH("params", "pmin0")) {
        p->pmin0 = atof(value);
    } else if (MATCH("params", "pmincrit")) {
        p->pmincrit = atof(value);
    } else if (MATCH("params", "puptakez")) {
      p->puptakez = atof(value);
    } else if (MATCH("params", "p_rate_par_weather")) {
        p->p_rate_par_weather = atof(value);
    } else if (MATCH("params", "passivesoilz")) {
        p->passivesoilz = atof(value);
    } else if (MATCH("params", "passivesoilnz")) {
        p->passivesoilnz = atof(value);
    } else if (MATCH("params", "passivesoilpz")) {
        p->passivesoilpz = atof(value);
    } else if (MATCH("params", "passncmax")) {
        p->passncmax = atof(value);
    } else if (MATCH("params", "passncmin")) {
        p->passncmin = atof(value);
    } else if (MATCH("params", "passpcmax")) {
        p->passpcmax = atof(value);
    } else if (MATCH("params", "passpcmin")) {
        p->passpcmin = atof(value);
    } else if (MATCH("params", "prescribed_leaf_NC")) {
        p->prescribed_leaf_NC = atof(value);
    } else if (MATCH("params", "prescribed_leaf_PC")) {
        p->prescribed_leaf_PC = atof(value);
    } else if (MATCH("params", "prime_y")) {
        p->prime_y = atof(value);
    } else if (MATCH("params", "prime_z")) {
        p->prime_z = atof(value);
    } else if (MATCH("params", "rateloss")) {
        p->rateloss = atof(value);
    } else if (MATCH("params", "prateloss")) {
        p->prateloss = atof(value);
    } else if (MATCH("params", "rateuptake")) {
        p->rateuptake = atof(value);
    } else if (MATCH("params", "prateuptake")) {
        p->prateuptake = atof(value);
    } else if (MATCH("params", "rdecay")) {
        p->rdecay = atof(value);
    } else if (MATCH("params", "rretrans")) {
      p->rretrans = atof(value);
    } else if (MATCH("params", "retransmob")) {
      p->retransmob = atof(value);
    } else if (MATCH("params", "sla")) {
        p->sla = atof(value);
    } else if (MATCH("params", "slowncmax")) {
        p->slowncmax = atof(value);
    } else if (MATCH("params", "slowncmin")) {
        p->slowncmin = atof(value);
    } else if (MATCH("params", "slowpcmax")) {
        p->slowpcmax = atof(value);
    } else if (MATCH("params", "slowpcmin")) {
        p->slowpcmin = atof(value);
    } else if (MATCH("params", "structcn")) {
        p->structcn = atof(value);
    } else if (MATCH("params", "structcp")) {
        p->structcp = atof(value);
    } else if (MATCH("params", "tsoil_in")) {
      p->tsoil_in = atof(value);
    } else if (MATCH("params", "wdecay")) {
        p->wdecay = atof(value);
    } else if (MATCH("params", "wretrans")) {
        p->wretrans = atof(value);
    }

    return (1);
}
