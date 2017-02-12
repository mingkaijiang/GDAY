#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "gday.h"

typedef struct {
    FILE *ifp;
    FILE *ofp;
    FILE *ofp_sd;
    char  cfg_fname[STRING_LENGTH];
    char  met_fname[STRING_LENGTH];
    char  out_fname[STRING_LENGTH];
    char  out_param_fname[STRING_LENGTH];
    char  git_hash[STRING_LENGTH];
    int   diagnosis;
    int   fixed_stem_nc;
    int   fixed_stem_pc;
    int   fixleafnc;
    int   fixleafpc;
    int   ncycle;
    int   pcycle;
    int   num_years;
    int   nuptake_model;
    int   puptake_model;
    int   print_options;
    int   use_eff_nc;
    int   num_months;
    int   total_num_months;
    char  git_code_ver[STRING_LENGTH];
    int   spin_up;
    int   PRINT_GIT;
    long  month_idx;

} control;


typedef struct {
    double activesoil;                  /* active C som pool (t/ha) */
    double activesoiln;                 /* active N som pool (t/ha) */
    double activesoilp;                 /* active P som pool (t/ha) */
    double inorgn;                      /* Inorganic soil N pool - dynamic (t/ha) */
    double inorgp;                      /* Inorganic soil P pool - dynamic (t/ha) */
    double inorgavlp;                   /* Inorganic soil P pool - available mineral P = lab + sorbed (t/ha) */
    double inorgssorbp;                 /* Inorganic soil P pool - strongly sorbed P (t/ha) */
    double inorgoccp;                   /* Inorganic soil P pool - occluded P (t/ha) */
    double inorgparp;                   /* Inorganic soil P pool - parent P (t/ha) */
    double lai;                         /* leaf area index m2 (leaf) m-2 (ground) */
    double fipar;
    double metabsoil;                   /* metabolic soil c (t/ha) */
    double metabsoiln;                  /* metabolic soil n (t/ha) */
    double metabsoilp;                  /* metabolic soil p (t/ha) */
    double metabsurf;                   /* metabolic surface c (t/ha) */
    double metabsurfn;                  /* metabolic surface n (t/ha) */
    double metabsurfp;                  /* metabolic surface p (t/ha) */
    double passivesoil;                 /* passive C som pool (t/ha) */
    double passivesoiln;                /* passive N som pool (t/ha) */
    double passivesoilp;                /* passive P som pool (t/ha) */
    double prev_sma;
    double root;                        /* root c (t/ha) */
    double rootn;                       /* root n (t/ha) */
    double rootp;                       /* root p (t/ha) */
    double shoot;                       /* shoot c (t/ha) */
    double shootn;                      /* shoot n (t/ha) */
    double shootp;                      /* shoot p (t/ha) */
    double slowsoil;                    /* slow C som pool (t/ha) */
    double slowsoiln;                   /* slow N som pool (t/ha) */
    double slowsoilp;                   /* slow P som pool (t/ha) */
    double stem;
    double stemn;                       /* Stem N (t/ha) = stemnimm + stemnmob */
    double stemnimm;
    double stemnmob;
    double stemp;                       /* Stem P (t/ha) = stempimm + stempmob */
    double stempimm;
    double stempmob;
    double structsoil;                  /* soil structural c (t/ha) */
    double structsoiln;                 /* soil structural n (t/ha) */
    double structsoilp;                 /* soil structural p (t/ha) */
    double structsurf;                  /* surface structural c (t/ha) */
    double structsurfn;                 /* surface structural n (t/ha) */
    double structsurfp;                 /* surface structural p (t/ha) */
    double shootnc;                     /* shoot nc ratio */
    double rootnc;                      /* root pn ratio */
    double shootpc;                     /* shoot pc ratio */
    double rootpc;                      /* root pc ratio */
    double c_to_alloc_shoot;
    double n_to_alloc_shoot;
    double p_to_alloc_shoot;
    double c_to_alloc_root;
    double n_to_alloc_root;
    double p_to_alloc_root;
    double c_to_alloc_stem;
    double n_to_alloc_stemmob;
    double n_to_alloc_stemimm;
    double p_to_alloc_stemmob;
    double p_to_alloc_stemimm;
    double anpp;                    /* aboveground NPP */
    double litterc;                 /* litter carbon */
    double littern;                 /* litter nitrogen */
    double litterp;                 /* litter phosphorus */
    double littercbg;               /* litter C belowground */
    double littercag;               /* litter C aboveground */
    double litternag;               /* litter N aboveground */
    double litternbg;               /* litter N belowground */
    double litterpag;               /* litter P aboveground */
    double litterpbg;               /* litter P belowground */
    double plantc;                  /* plant C */
    double plantn;                  /* plant N */
    double plantp;                  /* plant P */
    double totalc;                  /* total C */
    double totaln;                  /* total N */
    double totalp;                  /* total P */
    double soilc;                   /* Soil C */
    double soiln;                   /* soil N */
    double soilp;                   /* soil P */

} state;

typedef struct {
    double actncmax;                        /* Active pool (=1/3) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double actncmin;                        /* Active pool (=1/15) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double actpcmax;                        /* Active pool (=1/30) P:C ratio of new SOM - maximum [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993) */
    double actpcmin;                        /* Active pool (=1/80) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993) */
    double c_alloc_f;                       /* allocation to leaves at leaf n_crit and p_crit. */
    double c_alloc_r;                       /* allocation to roots at root n_crit and p_crit. */
    double c_alloc_s;                       /* allocation to stem at zero stem n/c and p/c. */
    double cfracts;                         /* carbon fraction of dry biomass */
    double co2_in;                          /* annual version co2 concentration ppm */
    double cue;                             /* carbon use efficiency, or the ratio of NPP to GPP */
    double decayrate[7];
    double fdecay;                          /* foliage turnover rate (1/yr) */
    double finesoil;                        /* clay+silt fraction */
    double fmleaf;
    double fmroot;
    double fretransn;                        /* foliage n retranslocation fraction - 46-57% in young E. globulus trees - see Corbeels et al 2005 ecological modelling 187, pg 463. Roughly 50% from review Aerts '96 */
    double fretransp;                       /* foliage p retranslocation fraction - 39.5-69 in Southern US FACE site - Finzi et al. 2001 Ecology  */
    double I0;                              /* annual version radiation MJ/m2/yr */
    double k1;                              /* P transfer rate coefficient from labile to secondary inorganic P pool */
    double k2;                              /* P transfer rate coefficient from secondary inorganic P to labile P */
    double k3;                              /* P transfer rate coefficient from secondary inorganic to occluded P pool */
    double kdec1;                           /* surface structural decay rate (1/yr) */
    double kdec2;                           /* surface metabolic decay rate (1/yr) */
    double kdec3;                           /* soil structural decay rate (1/yr) */
    double kdec4;                           /* soil metabolic decay rate(1/yr) */
    double kdec5;                           /* active pool decay rate (1/yr) */
    double kdec6;                           /* slow pool decay rate (1/yr) */
    double kdec7;                           /* passive pool decay rate (1/yr) */
    double kext;                            /* extinction coefficient for light  */
    double kr;                              /* N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha) see Silvia's PhD, Dewar and McM, 96. */
    double krp;                             /* P uptake coefficent */
    double ligroot;                         /* lignin-to-biomass ratio in root litter; Values from White et al. = 0.22  - Value in Smith et al. 2013 = 0.16, note subtly difference in eqn C9. */
    double ligshoot;                        /* lignin-to-biomass ratio in leaf litter; Values from White et al. DBF = 0.18; ENF = 0.24l; GRASS = 0.09; Shrub = 0.15 - Value in smith et al. 2013 = 0.2, note subtly difference in eqn C9. */
    double liteffnc;
    double lue0;                            /* maximum LUE in kg C GJ-1 */
    double ncmaxf;                          /* max N:C ratio of foliage in old stand, if the same as young=no effect */
    double nref;                            /* leaf nc for saturation of photosynthesis */
    double ncmaxr;                          /* max N:C ratio of roots */
    double ncrfac;                          /* N:C of fine root prodn / N:C of leaf prodn */
    double ncwimmz;                         /* N alloc param: Immobile stem N C at zero leaf N C */
    double ncwnewz;                         /* N alloc param: New stem ring N:C at zero leaf N:C (mobile) */
    double ndep_in;                         /* annual version ndep t/ha/yr */
    double nfix_in;                         /* annual version nfix t/ha/yr */
    double nf_min;                          /* leaf N:C minimum N concentration which allows productivity */
    double nmin0;                           /* mineral N pool corresponding to Actnc0,etc (g/m2) */
    double nmincrit;                        /* Critical mineral N pool at max soil N:C (g/m2) (Parton et al 1993, McMurtrie et al 2001). */
    double nuptakez;                        /* constant N uptake per year (1/yr) */
    double p_rate_par_weather;              /* parent P material weathering rate [yr-1] */
    double passivesoilnz;
    double passivesoilpz;
    double passivesoilz;
    double passncmax;                       /* Passive pool (=1/7) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double passncmin;                       /* Passive pool (=1/10) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double passpcmax;                       /* Passive pool (=1/20) P:C ratio of new SOM - maximum [units: gP/gC] */
    double passpcmin;                       /* Passive pool (=1/200) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC] */
    double pcmaxf;                          /* max P:C ratio of foliage in old stand, if the same as young=no effect */
    double pcmaxr;                          /* max P:C ratio of roots */
    double pcrfac;                          /* P:C of fine root prodp / P:C c of leaf prodp */
    double pcwimmz;                         /* P alloc param: Immobile stem P C at zero leaf P C */
    double pcwnewz;                         /* P alloc param: New stem ring P:C at zero leaf P:C (mobile) */
    double pdep_in;                         /* Annual version p deposition t/ha/yr */
    double pf_min;                          /* leaf P:C minimum P concentration which allows productivity */
    double pmin0;                           /* mineral P pool corresponding to Actpc0,etc (g/m2) */
    double pmincrit;                        /* Critical mineral P pool at max soil P:C (g/m2) */
    double prateloss;                       /* Rate of P loss from mineral P pool (/yr), Ref Wang et al., 2007, GB1018 */
    double prateuptake;                     /* Rate of P uptake from mineral P pool (/yr), guess value */
    double prescribed_leaf_NC;              /* If the N-Cycle is switched off this needs to be set, e.g. 0.03 */
    double prescribed_leaf_PC;              /* If the P-Cycle is switched off this needs to be set, e.g. 0.00249 */
    double puptakez;                        /* constant P uptake per year (1/yr) */
    double rateloss;                        /* Rate of N loss from mineral N pool (/yr) */
    double rateuptake;                      /* Rate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf Seems to correspond to very low NPP values */
    double rdecay;                          /* root turnover rate (1/yr) */
    double rretrans;                        /* root retranslocation coefficient */
    double retransmob;                      /* mobilized wood retranslocation coefficient */
    double sla;                             /* specific leaf area (m2 one-sided/kg DW) */
    double slowncmax;                       /* Slow pool (=1/15) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double slowncmin;                       /* Slow pool (=1/40) N:C of new SOM - when Nmin=Nmin0" [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double slowpcmax;                       /* Slow pool (=1/90) P:C ratio of new SOM - maximum [units: gP/gC]. */
    double slowpcmin;                       /* Slow pool (=1/200) P:C of new SOM - when Pmin=Pmin0" [units: gP/gC]. */
    double structcn;                        /* C:N ratio of structural bit of litter input */
    double structcp;                        /* C:P ratio of structural bit of litter input, Ref Attiwill 1980, Aus. J. Bot. 28, 199-222 Table 9 sum of branch, stem, sap and heartwood; */
    double tsoil_in;                        /* annual version tsoil [degree C] */
    double wdecay;                          /* wood turnover rate (1/yr) */
    double wretrans;                        /* wood retranslocation coefficient */
    double prime_y;
    double prime_z;

} params;

typedef struct {

    double *year;
    double *prjmonth;
    double *par;
    double *tsoil;
    double *co2;
    double *ndep;
    double *nfix;       
    double *pdep;


} met_arrays;


typedef struct {

    double par;
    double Ca;
    double ndep;
    double nfix;       
    double pdep;
    double tsoil;
    double year;
    double co2;


} met;



typedef struct {
    double gpp_gCm2;
    double npp_gCm2;
    double gpp_am;
    double gpp_pm;
    double gpp;
    double npp;
    double npp_photo;
    double gpp_photo;
    double nep;
    double auto_resp;
    double hetero_resp;
    double retransn;         /* plnat n retranslocation */
    double retransp;        /* plant p retranslocation */
    double apar;

    /* n */
    double nuptake;         /* n plant uptake rate */
    double nloss;           /* n loss by leaching and volatilisation */
    double npassive;        /* n passive -> active */
    double ngross;          /* N gross mineralisation */
    double nimmob;          /* N immobilisation in SOM */
    double nlittrelease;    /* N rel litter = struct + metab */
    double activelossf;     /* frac of active C -> CO2 */
    double nmineralisation;

    /* p */
    double puptake;         /* P plant uptake rate */
    double ploss;           /* P loss by leaching */
    double ppassive;        /* P passive -> active */
    double pgross;          /* P gross mineralisation */
    double pimmob;          /* P immobilisation in SOM */
    double plittrelease;    /* P rel litter = struct + metab */
    double pmineralisation; /* P net mineralised */

    /* daily C production */
    double cpleaf;
    double cproot;
    double cpstem;

    /* daily N production */
    double npleaf;
    double nproot;
    double npstemimm;
    double npstemmob;

    /* daily P production */
    double ppleaf;
    double pproot;
    double ppstemimm;
    double ppstemmob;

    /* dying stuff */
    double deadleaves;      /* Leaf litter C production (t/ha/yr) */
    double deadroots;       /* Root litter C production (t/ha/yr) */
    double deadstems;       /* Stem litter C production (t/ha/yr) */
    double deadleafn;       /* Leaf litter N production (t/ha/yr) */
    double deadrootn;       /* Root litter N production (t/ha/yr) */
    double deadstemn;       /* Stem litter N production (t/ha/yr) */
    double deadleafp;       /* Leaf litter P production (t/ha/yr) */
    double deadrootp;       /* Root litter P production (t/ha/yr) */
    double deadstemp;       /* Stem litter P production (t/ha/yr) */

    /* retranslocation */
    double leafretransn;    /* N retranslocation leaf */
    double leafretransp;    /* P version of leafretransn */
    double rootretransn;
    double rootretransp;
    double stemretransn;
    double stemretransp;


    /* C, N & P Surface litter */
    double surf_struct_litter;
    double surf_metab_litter;
    double n_surf_struct_litter;
    double n_surf_metab_litter;
    double p_surf_struct_litter;
    double p_surf_metab_litter;

    /* C, N & P Root Litter */
    double soil_struct_litter;
    double soil_metab_litter;
    double n_soil_struct_litter;
    double n_soil_metab_litter;
    double p_soil_struct_litter;
    double p_soil_metab_litter;


    /* C, N & P litter fluxes to slow pool */
    double surf_struct_to_slow;
    double soil_struct_to_slow;
    double n_surf_struct_to_slow;
    double n_soil_struct_to_slow;
    double p_surf_struct_to_slow;
    double p_soil_struct_to_slow;

    /* C, N & P litter fluxes to active pool */
    double surf_struct_to_active;
    double soil_struct_to_active;
    double n_surf_struct_to_active;
    double n_soil_struct_to_active;
    double p_surf_struct_to_active;
    double p_soil_struct_to_active;

    /* Metabolic fluxes to active pool */
    double surf_metab_to_active;
    double soil_metab_to_active;
    double n_surf_metab_to_active;
    double n_soil_metab_to_active;
    double p_surf_metab_to_active;
    double p_soil_metab_to_active;

    /* fluxes out of active pool */
    double active_to_slow;
    double active_to_passive;
    double n_active_to_slow;
    double n_active_to_passive;
    double p_active_to_slow;
    double p_active_to_passive;

    /* Fluxes from slow pool */
    double slow_to_active;
    double slow_to_passive;
    double n_slow_to_active;
    double n_slow_to_passive;
    double p_slow_to_active;
    double p_slow_to_passive;

    /* fluxes from passive pool */
    double passive_to_active;
    double n_passive_to_active;
    double p_passive_to_active;

    /* C, N & P source fluxes from the active, slow and passive pools */
    double c_into_active;
    double c_into_slow;
    double c_into_passive;

    /* inorganic P flux exchanges */
    double p_avl_to_ssorb;
    double p_ssorb_to_avl;
    double p_ssorb_to_occ;
    double p_par_to_avl;
    double p_atm_dep;


    /* CO2 flows to the air */
    double co2_to_air[7];

    /* C allocated fracs  */
    double alleaf;             /* allocation to leaf */
    double alroot;             /* allocation to fine root */
    double alstem;             /* allocation to stems */

    /* Misc stuff */
    double tfac_soil_decomp;
    double co2_rel_from_surf_struct_litter;
    double co2_rel_from_soil_struct_litter;
    double co2_rel_from_surf_metab_litter;
    double co2_rel_from_soil_metab_litter;
    double co2_rel_from_active_pool;
    double co2_rel_from_slow_pool;
    double co2_rel_from_passive_pool;

    double ninflow;

} fluxes;


typedef struct {
    double  *ystart;
    double   *yscal;
    double   *y;
    double   *dydx;
    double   *xp;
	  double  **yp;
    int       N;
    int       kmax;

    double   *ak2;
    double   *ak3;
    double   *ak4;
    double   *ak5;
    double   *ak6;
    double   *ytemp;
    double   *yerr;



} nrutil;

#endif
