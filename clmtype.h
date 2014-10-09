#ifndef CNM_HEADER
    #define CNM_HEADER

#define SHR_CONST_TKFRZ     273.15

/*****************************************************************************
 * Begin definition of conservation check structures
 ****************************************************************************/
/*
 * carbon balance structure
 */
typedef carbon_balance_type
{
   double begcb;        /* carbon mass, beginning of time step (gC/m**2) */
   double endcb;        /* carbon mass, end of time step (gC/m**2) */
   double errcb;        /* carbon balance error for the timestep (gC/m**2) */
} carbon_balance_type;

/*
 * nitrogen balance structure
 */
typedef nitrogen_balance_type
{
   double begnb;        /* nitrogen mass, beginning of time step (gN/m**2) */
   double endnb;        /* nitrogen mass, end of time step (gN/m**2) */
   double errnb;        /* nitrogen balance error for the timestep (gN/m**2) */
} nitrogen_balance_type;

/*****************************************************************************
 * End definition of conservation check structures
 ****************************************************************************/

/*****************************************************************************
 * Begin definition of structures defined at the pft_type level
 ****************************************************************************/
/*
 * physical state variables structure
 */
typedef struct pstate_type
{
   int frac_veg_nosno;          /* fraction of vegetation not covered by snow (0 OR 1) [-] */
   int frac_veg_nosno_alb;          /* fraction of vegetation not covered by snow (0 OR 1) [-] */
   double emv;      /* vegetation emissivity */
   double z0mv;     /* roughness length over vegetation, momentum [m] */
   double z0hv;     /* roughness length over vegetation, sensible heat [m] */
   double z0qv;     /* roughness length over vegetation, latent heat [m] */
   double *rootfr;       /* fraction of roots in each soil layer  (nlevgrnd) */
   double *rootr;        /* effective fraction of roots in each soil layer  (nlevgrnd) */
   double *rresis;       /* root resistance by layer (0-1)  (nlevgrnd) */
   double dewmx;        /* Maximum allowed dew [mm] */
   double rssun;        /* sunlit stomatal resistance (s/m) */
   double rssha;        /* shaded stomatal resistance (s/m) */
   double laisun;       /* sunlit projected leaf area index */
   double laisha;       /* shaded projected leaf area index */
   double btran;        /* transpiration wetness factor (0 to 1) */
   double fsun;     /* sunlit fraction of canopy */
   double tlai;     /* one-sided leaf area index, no burying by snow */
   double tsai;     /* one-sided stem area index, no burying by snow */
   double elai;     /* one-sided leaf area index with burying by snow */
   double esai;     /* one-sided stem area index with burying by snow */
   double fwet;     /* fraction of canopy that is wet (0 to 1) */
   double fdry;     /* fraction of foliage that is green and dry [-] (new) */
   double dt_veg;       /* change in t_veg, last iteration (Kelvin) */
   double htop;     /* canopy top (m) */
   double hbot;     /* canopy bottom (m) */
   double z0m;      /* momentum roughness length (m) */
   double displa;       /* displacement height (m) */
   double albd;     /* surface albedo (direct)                       (numrad) */
   double albi;     /* surface albedo (indirect)                      (numrad) */
   double fabd;     /* flux absorbed by veg per unit direct flux     (numrad) */
   double fabi;     /* flux absorbed by veg per unit diffuse flux    (numrad) */
   double ftdd;     /* down direct flux below veg per unit dir flx   (numrad) */
   double ftid;     /* down diffuse flux below veg per unit dir flx  (numrad) */
   double ftii;     /* down diffuse flux below veg per unit dif flx  (numrad) */
   double u10;      /* 10-m wind (m/s) (for dust model) */
   double u10_clm;      /* 10-m wind (m/s) */
   double va;       /* atmospheric wind speed plus convective velocity (m/s) */
   double ram1;     /* aerodynamical resistance (s/m) */
   double fv;       /* friction velocity (m/s) (for dust model) */
   double forc_hgt_u_pft;       /* wind forcing height (10m+z0m+d) (m) */
   double forc_hgt_t_pft;       /* temperature forcing height (10m+z0m+d) (m) */
   double forc_hgt_q_pft;       /* specific humidity forcing height (10m+z0m+d) (m) */
   /* Variables for prognostic crop model */
   double hdidx;        /*  cold hardening index? */
   double cumvd;        /*  cumulative vernalization d?ependence? */
   double htmx;     /*  max hgt attained by a crop during yr (m) */
   double vf;       /*  vernalization factor for cereal */
   double gddmaturity;      /*  growing degree days (gdd) needed to harvest (ddays) */
   double gdd0;     /*  growing degree-days base  0C from planting  (ddays) */
   double gdd8;     /*  growing degree-days base  8C from planting  (ddays) */
   double gdd10;        /*  growing degree-days base 10C from planting  (ddays) */
   double gdd020;       /*  20-year average of gdd0                     (ddays) */
   double gdd820;       /*  20-year average of gdd8                     (ddays) */
   double gdd1020;      /*  20-year average of gdd10                    (ddays) */
   double gddplant;     /*  accum gdd past planting date for crop       (ddays) */
   double gddtsoi;      /*  growing degree-days from planting (top two soil layers) (ddays) */
   double huileaf;      /*  heat unit index needed from planting to leaf emergence */
   double huigrain;     /*  heat unit index needed to reach vegetative maturity */
   double aleafi;       /*  saved leaf allocation coefficient from phase 2 */
   double astemi;       /*  saved stem allocation coefficient from phase 2 */
   double aleaf;        /*  leaf allocation coefficient */
   double astem;        /*  stem allocation coefficient */
   int croplive;            /*  Flag, true if planted, not harvested*/
   int cropplant;           /*  Flag, true if planted*/
   int harvdate;            /*  harvest date*/
                            /* cropplant and harvdate could be 2D to facilitate rotation */
   int idop;            /*  date of planting*/
   int peaklai;         /*  1: max allowed lai; 0: not at max*/
   double vds;      /* deposition velocity term (m/s) (for dry dep SO4, NH4NO3) */
   /* new variables for CN code */
   double slasun;       /* specific leaf area for sunlit canopy, projected area basis (m^2/gC) */
   double slasha;       /* specific leaf area for shaded canopy, projected area basis (m^2/gC) */
   double lncsun;       /* leaf N concentration per unit projected LAI (gN leaf/m^2) */
   double lncsha;       /* leaf N concentration per unit projected LAI (gN leaf/m^2) */
   double vcmxsun;      /* sunlit leaf Vcmax (umolCO2/m^2/s) */
   double vcmxsha;      /* shaded leaf Vcmax (umolCO2/m^2/s) */
   double gdir;     /* leaf projection in solar direction (0 to 1) */
   double omega;        /* fraction of intercepted radiation that is scattered (0 to 1) */
   double eff_kid;      /* effective extinction coefficient for indirect from direct */
   double eff_kii;      /* effective extinction coefficient for indirect from indirect */
   double sun_faid;     /* fraction sun canopy absorbed indirect from direct */
   double sun_faii;     /* fraction sun canopy absorbed indirect from indirect */
   double sha_faid;     /* fraction shade canopy absorbed indirect from direct */
   double sha_faii;     /* fraction shade canopy absorbed indirect from indirect */
   double cisun;        /* sunlit intracellular CO2 (Pa) */
   double cisha;        /* shaded intracellular CO2 (Pa) */
   double alphapsnsun;      /* sunlit 13c fractionation ([]) */
   double alphapsnsha;      /* shaded 13c fractionation ([]) */
   double sandfrac;     /*  sand fraction */
   double clayfrac;     /*  clay fraction */
   /* for dry deposition of chemical tracers */
   double mlaidiff;     /*  difference between lai month one and month two */
   double rb1;      /*  aerodynamical resistance (s/m) */
   double annlai;       /*  12 months of monthly lai from input data set   */

   double *bsw2;     /* Clapp and Hornberger "b" for CN code */
   double *psisat;       /* soil water potential at saturation for CN code (MPa) */
   double *vwcsat;       /* volumetric water content at saturation for CN code (m3/m3) */
   double decl;     /*  solar declination angle (radians) */
   double coszen;       /* cosine of solar zenith angle */
   double *soilpsi;      /* soil water potential in each soil layer (MPa) */
   double fpi;      /* fraction of potential immobilization (no units) */
   double fpg;      /* fraction of potential gpp (no units) */
   double annsum_counter;       /* seconds since last annual accumulator turnover */
   double annsum_npp;      /* annual sum of NPP, averaged from pft-level (gC/m2/yr) */
   double annavg_t2m;      /* annual average of 2m air temperature, averaged from pft-level (K) */
   double *watfc;        /* volumetric soil water at field capacity (nlevsoi) */
    
} pstate_type;

/*
 * ecophysiological constants structure
 */
typedef struct epc_type
{
   int noveg;           /* value for not vegetated*/
   int tree;            /* tree or not?*/
   double smpso;        /* soil water potential at full stomatal opening (mm) */
   double smpsc;        /* soil water potential at full stomatal closure (mm) */
   double fnitr;        /* foliage nitrogen limitation factor (-) */
   double foln;     /* foliage nitrogen (%) */
   double dleaf;        /* characteristic leaf dimension (m) */
   double c3psn;        /* photosynthetic pathway: 0. = c4, 1. = c3 */
   double mp;       /* slope of conductance-to-photosynthesis relationship */
   double qe25;     /* quantum efficiency at 25C (umol CO2 / umol photon) */
   double xl;       /* leaf/stem orientation index */
   double *rhol;     /* leaf reflectance: 1=vis, 2=nir   (numrad) */
   double *rhos;     /* stem reflectance: 1=vis, 2=nir   (numrad) */
   double *taul;     /* leaf transmittance: 1=vis, 2=nir (numrad) */
   double *taus;     /* stem transmittance: 1=vis, 2=nir (numrad) */
   double z0mr;     /* ratio of momentum roughness length to canopy top height (-) */
   double displar;      /* ratio of displacement height to canopy top height (-) */
   double roota_par;        /* CLM rooting distribution parameter [1/m] */
   double rootb_par;        /* CLM rooting distribution parameter [1/m] */
   double sla;      /* specific leaf area [m2 leaf g-1 carbon] */
   /* new variables for CN code */
   double dwood;        /* wood density (gC/m3) */
   double slatop;       /* specific leaf area at top of canopy, projected area basis [m^2/gC] */
   double dsladlai;     /* dSLA/dLAI, projected area basis [m^2/gC] */
   double leafcn;       /* leaf C:N (gC/gN) */
   double flnr;     /* fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf) */
   double woody;        /* binary flag for woody lifeform (1=woody, 0=not woody) */
   double lflitcn;      /* leaf litter C:N (gC/gN) */
   double frootcn;      /* fine root C:N (gC/gN) */
   double livewdcn;     /* live wood (phloem and ray parenchyma) C:N (gC/gN) */
   double deadwdcn;     /* dead wood (xylem and heartwood) C:N (gC/gN) */
   double graincn;      /* grain C:N (gC/gN) for prognostic crop model */
   double froot_leaf;       /* allocation parameter: new fine root C per new leaf C (gC/gC) */
   double stem_leaf;        /* allocation parameter: new stem c per new leaf C (gC/gC) */
   double croot_stem;       /* allocation parameter: new coarse root C per new stem C (gC/gC) */
   double flivewd;      /* allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units) */
   double fcur;     /* allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage */
   double lf_flab;      /* leaf litter labile fraction */
   double lf_fcel;      /* leaf litter cellulose fraction */
   double lf_flig;      /* leaf litter lignin fraction */
   double fr_flab;      /* fine root litter labile fraction */
   double fr_fcel;      /* fine root litter cellulose fraction */
   double fr_flig;      /* fine root litter lignin fraction */
   double leaf_long;        /* leaf longevity (yrs) */
   double evergreen;        /* binary flag for evergreen leaf habit (0 or 1) */
   double stress_decid;     /* binary flag for stress-deciduous leaf habit (0 or 1) */
   double season_decid;     /* binary flag for seasonal-deciduous leaf habit (0 or 1) */
   /* new variables for fire code */
   double resist;       /* resistance to fire (no units) */
} epc_type;

/*****************************************************************************
 * pft ecophysiological variables structure
 ****************************************************************************/
typedef struct epv_type
{
   double dormant_flag;     /* dormancy flag */
   double days_active;      /* number of days since last dormancy */
   double onset_flag;       /* onset flag */
   double onset_counter;        /* onset days counter */
   double onset_gddflag;        /* onset flag for growing degree day sum */
   double onset_fdd;        /* onset freezing degree days counter */
   double onset_gdd;        /* onset growing degree days */
   double onset_swi;        /* onset soil water index */
   double offset_flag;      /* offset flag */
   double offset_counter;       /* offset days counter */
   double offset_fdd;       /* offset freezing degree days counter */
   double offset_swi;       /* offset soil water index */
   double lgsf;     /* long growing season factor [0-1] */
   double bglfr;        /* background litterfall rate (1/s) */
   double bgtr;     /* background transfer growth rate (1/s) */
   double dayl;     /* daylength (seconds) */
   double prev_dayl;        /* daylength from previous timestep (seconds) */
   double annavg_t2m;       /* annual average 2m air temperature (K) */
   double tempavg_t2m;      /* temporary average 2m air temperature (K) */
   double gpp;      /* GPP flux before downregulation (gC/m2/s) */
   double availc;       /* C flux available for allocation (gC/m2/s) */
   double xsmrpool_recover;     /* C flux assigned to recovery of negative cpool (gC/m2/s) */
   double xsmrpool_c13ratio;        /* C13/C(12+13) ratio for xsmrpool (proportion) */
   double alloc_pnow;       /* fraction of current allocation to display as new growth (DIM) */
   double c_allometry;      /* C allocation index (DIM) */
   double n_allometry;      /* N allocation index (DIM) */
   double plant_ndemand;        /* N flux required to support initial GPP (gN/m2/s) */
   double tempsum_potential_gpp;    /* temporary annual sum of potential GPP */
   double annsum_potential_gpp;     /* annual sum of potential GPP */
   double tempmax_retransn;     /* temporary annual max of retranslocated N pool (gN/m2) */
   double annmax_retransn;      /* annual max of retranslocated N pool (gN/m2) */
   double avail_retransn;       /* N flux available from retranslocation pool (gN/m2/s) */
   double plant_nalloc;     /* total allocated N flux (gN/m2/s) */
   double plant_calloc;     /* total allocated C flux (gC/m2/s) */
   double excess_cflux;     /* C flux not allocated due to downregulation (gC/m2/s) */
   double downreg;      /* fractional reduction in GPP due to N limitation (DIM) */
   double prev_leafc_to_litter;     /* previous timestep leaf C litterfall flux (gC/m2/s) */
   double prev_frootc_to_litter;    /* previous timestep froot C litterfall flux (gC/m2/s) */
   double tempsum_npp;      /* temporary annual sum of NPP (gC/m2/yr) */
   double annsum_npp;       /* annual sum of NPP (gC/m2/yr) */
   double tempsum_litfall;      /* temporary annual sum of litfall (gC/m2/yr) */
   double annsum_litfall;       /* annual sum of litfall (gC/m2/yr) */
   double rc13_canair;      /* C13O2/C12O2 in canopy air */
   double rc13_psnsun;      /* C13O2/C12O2 in sunlit canopy psn flux */
   double rc13_psnsha;      /* C13O2/C12O2 in shaded canopy psn flux */
} epv_type;

/*
 * pft carbon state variables structure
 */
typedef struct cstate_type
{
   double leafcmax;     /*  (gC/m2) ann max leaf C */
   /* variables for prognostic crop model */
   double grainc;       /*  (gC/m2) grain C */
   double grainc_storage;       /*  (gC/m2) grain C storage */
   double grainc_xfer;      /*  (gC/m2) grain C transfer */

   double leafc;        /*  (gC/m2) leaf C */
   double leafc_storage;        /*  (gC/m2) leaf C storage */
   double leafc_xfer;       /*  (gC/m2) leaf C transfer */
   double frootc;       /*  (gC/m2) fine root C */
   double frootc_storage;       /*  (gC/m2) fine root C storage */
   double frootc_xfer;      /*  (gC/m2) fine root C transfer */
   double livestemc;        /*  (gC/m2) live stem C */
   double livestemc_storage;        /*  (gC/m2) live stem C storage */
   double livestemc_xfer;       /*  (gC/m2) live stem C transfer */
   double deadstemc;        /*  (gC/m2) dead stem C */
   double deadstemc_storage;        /*  (gC/m2) dead stem C storage */
   double deadstemc_xfer;       /*  (gC/m2) dead stem C transfer */
   double livecrootc;       /*  (gC/m2) live coarse root C */
   double livecrootc_storage;       /*  (gC/m2) live coarse root C storage */
   double livecrootc_xfer;      /*  (gC/m2) live coarse root C transfer */
   double deadcrootc;       /*  (gC/m2) dead coarse root C */
   double deadcrootc_storage;       /*  (gC/m2) dead coarse root C storage */
   double deadcrootc_xfer;      /*  (gC/m2) dead coarse root C transfer */
   double gresp_storage;        /*  (gC/m2) growth respiration storage */
   double gresp_xfer;       /*  (gC/m2) growth respiration transfer */
   double cpool;        /*  (gC/m2) temporary photosynthate C pool */
   double xsmrpool;     /*  (gC/m2) abstract C pool to meet excess MR demand */
   double pft_ctrunc;       /*  (gC/m2) pft-level sink for C truncation */
   /* summary (diagnostic) state variables, not involved in mass balance */
   double dispvegc;     /*  (gC/m2) displayed veg carbon, excluding storage and cpool */
   double storvegc;     /*  (gC/m2) stored vegetation carbon, excluding cpool */
   double totvegc;      /*  (gC/m2) total vegetation carbon, excluding cpool */
   double totpftc;      /*  (gC/m2) total pft-level carbon, including cpool */
   double woodc;        /*  (gC/m2) wood C */
} cstate_type;

/*
 * pft nitrogen state variables structure
 */
typedef struct nstate_type
{
    /* variables for prognostic crop model */
   double grainn;       /*  (gN/m2) grain N  */
   double grainn_storage;       /*  (gN/m2) grain N storage */
   double grainn_xfer;      /*  (gN/m2) grain N transfer */

   double leafn;        /*  (gN/m2) leaf N  */
   double leafn_storage;        /*  (gN/m2) leaf N storage */
   double leafn_xfer;       /*  (gN/m2) leaf N transfer */
   double frootn;       /*  (gN/m2) fine root N */
   double frootn_storage;       /*  (gN/m2) fine root N storage */
   double frootn_xfer;      /*  (gN/m2) fine root N transfer */
   double livestemn;        /*  (gN/m2) live stem N */
   double livestemn_storage;        /*  (gN/m2) live stem N storage */
   double livestemn_xfer;       /*  (gN/m2) live stem N transfer */
   double deadstemn;        /*  (gN/m2) dead stem N */
   double deadstemn_storage;        /*  (gN/m2) dead stem N storage */
   double deadstemn_xfer;       /*  (gN/m2) dead stem N transfer */
   double livecrootn;       /*  (gN/m2) live coarse root N */
   double livecrootn_storage;       /*  (gN/m2) live coarse root N storage */
   double livecrootn_xfer;      /*  (gN/m2) live coarse root N transfer */
   double deadcrootn;       /*  (gN/m2) dead coarse root N */
   double deadcrootn_storage;       /*  (gN/m2) dead coarse root N storage */
   double deadcrootn_xfer;      /*  (gN/m2) dead coarse root N transfer */
   double retransn;     /*  (gN/m2) plant pool of retranslocated N */
   double npool;        /*  (gN/m2) temporary plant N pool */
   double pft_ntrunc;       /*  (gN/m2) pft-level sink for N truncation */
   /* summary (diagnostic) state variables, not involved in mass balance */
   double dispvegn;     /*  (gN/m2) displayed veg nitrogen, excluding storage */
   double storvegn;     /*  (gN/m2) stored vegetation nitrogen */
   double totvegn;      /*  (gN/m2) total vegetation nitrogen */
   double totpftn;      /*  (gN/m2) total pft-level nitrogen */

   double cwdn;     /*  (gN/m2) coarse woody debris N */
   double litr1n;       /*  (gN/m2) litter labile N */
   double litr2n;       /*  (gN/m2) litter cellulose N */
   double litr3n;       /*  (gN/m2) litter lignin N */
   double soil1n;       /*  (gN/m2) soil organic matter N (fast pool) */
   double soil2n;       /*  (gN/m2) soil organic matter N (medium pool) */
   double soil3n;       /*  (gN/m2) soil orgainc matter N (slow pool) */
   double soil4n;       /*  (gN/m2) soil orgainc matter N (slowest pool) */
   double sminn;        /*  (gN/m2) soil mineral N */
   double col_ntrunc;       /*  (gN/m2) column-level sink for N truncation */

   double totlitn;      /*  (gN/m2) total litter nitrogen */
   double totsomn;      /*  (gN/m2) total soil organic matter nitrogen */
   double totecosysn;       /*  (gN/m2) total ecosystem nitrogen, incl veg  */
   double totcoln;      /*  (gN/m2) total column nitrogen, incl veg */
} nstate_type;

/*
 * pft carbon flux variables structure
 */
typedef struct cflux_type
{
   double psnsun;       /* sunlit leaf photosynthesis (umol CO2 /m**2/ s) */
   double psnsha;       /* shaded leaf photosynthesis (umol CO2 /m**2/ s) */
   double fpsn;     /* photosynthesis (umol CO2 /m**2 /s) */
   double fco2;     /* net CO2 flux (umol CO2 /m**2 /s) [+ = to atm] */
   /* new variables for CN code */
   /* gap mortality fluxes */
   double m_leafc_to_litter;        /*  leaf C mortality (gC/m2/s) */
   double m_leafc_storage_to_litter;        /*  leaf C storage mortality (gC/m2/s) */
   double m_leafc_xfer_to_litter;       /*  leaf C transfer mortality (gC/m2/s) */
   double m_frootc_to_litter;       /*  fine root C mortality (gC/m2/s) */
   double m_frootc_storage_to_litter;       /*  fine root C storage mortality (gC/m2/s) */
   double m_frootc_xfer_to_litter;      /*  fine root C transfer mortality (gC/m2/s) */
   double m_livestemc_to_litter;        /*  live stem C mortality (gC/m2/s) */
   double m_livestemc_storage_to_litter;        /*  live stem C storage mortality (gC/m2/s) */
   double m_livestemc_xfer_to_litter;       /*  live stem C transfer mortality (gC/m2/s) */
   double m_deadstemc_to_litter;        /*  dead stem C mortality (gC/m2/s) */
   double m_deadstemc_storage_to_litter;        /*  dead stem C storage mortality (gC/m2/s) */
   double m_deadstemc_xfer_to_litter;       /*  dead stem C transfer mortality (gC/m2/s) */
   double m_livecrootc_to_litter;       /*  live coarse root C mortality (gC/m2/s) */
   double m_livecrootc_storage_to_litter;       /*  live coarse root C storage mortality (gC/m2/s) */
   double m_livecrootc_xfer_to_litter;      /*  live coarse root C transfer mortality (gC/m2/s) */
   double m_deadcrootc_to_litter;       /*  dead coarse root C mortality (gC/m2/s) */
   double m_deadcrootc_storage_to_litter;       /*  dead coarse root C storage mortality (gC/m2/s) */
   double m_deadcrootc_xfer_to_litter;      /*  dead coarse root C transfer mortality (gC/m2/s) */
   double m_gresp_storage_to_litter;        /*  growth respiration storage mortality (gC/m2/s) */
   double m_gresp_xfer_to_litter;       /*  growth respiration transfer mortality (gC/m2/s) */
   /* harvest mortality fluxes */
   double hrv_leafc_to_litter;      /*  leaf C harvest mortality (gC/m2/s) */
   double hrv_leafc_storage_to_litter;      /*  leaf C storage harvest mortality (gC/m2/s) */
   double hrv_leafc_xfer_to_litter;     /*  leaf C transfer harvest mortality (gC/m2/s) */
   double hrv_frootc_to_litter;     /*  fine root C harvest mortality (gC/m2/s) */
   double hrv_frootc_storage_to_litter;     /*  fine root C storage harvest mortality (gC/m2/s) */
   double hrv_frootc_xfer_to_litter;        /*  fine root C transfer harvest mortality (gC/m2/s) */
   double hrv_livestemc_to_litter;      /*  live stem C harvest mortality (gC/m2/s) */
   double hrv_livestemc_storage_to_litter;      /*  live stem C storage harvest mortality (gC/m2/s) */
   double hrv_livestemc_xfer_to_litter;     /*  live stem C transfer harvest mortality (gC/m2/s) */
   double hrv_deadstemc_to_prod10c;     /*  dead stem C harvest to 10-year product pool (gC/m2/s) */
   double hrv_deadstemc_to_prod100c;        /*  dead stem C harvest to 100-year product pool (gC/m2/s) */
   double hrv_deadstemc_storage_to_litter;      /*  dead stem C storage harvest mortality (gC/m2/s) */
   double hrv_deadstemc_xfer_to_litter;     /*  dead stem C transfer harvest mortality (gC/m2/s) */
   double hrv_livecrootc_to_litter;     /*  live coarse root C harvest mortality (gC/m2/s) */
   double hrv_livecrootc_storage_to_litter;     /*  live coarse root C storage harvest mortality (gC/m2/s) */
   double hrv_livecrootc_xfer_to_litter;        /*  live coarse root C transfer harvest mortality (gC/m2/s) */
   double hrv_deadcrootc_to_litter;     /*  dead coarse root C harvest mortality (gC/m2/s) */
   double hrv_deadcrootc_storage_to_litter;     /*  dead coarse root C storage harvest mortality (gC/m2/s) */
   double hrv_deadcrootc_xfer_to_litter;        /*  dead coarse root C transfer harvest mortality (gC/m2/s) */
   double hrv_gresp_storage_to_litter;      /*  growth respiration storage harvest mortality (gC/m2/s) */
   double hrv_gresp_xfer_to_litter;     /*  growth respiration transfer harvest mortality (gC/m2/s) */
   double hrv_xsmrpool_to_atm;      /*  excess MR pool harvest mortality (gC/m2/s) */
   /* PFT-level fire fluxes */
   double m_leafc_to_fire;      /*  leaf C fire loss (gC/m2/s) */
   double m_leafc_storage_to_fire;      /*  leaf C storage fire loss (gC/m2/s) */
   double m_leafc_xfer_to_fire;     /*  leaf C transfer fire loss (gC/m2/s) */
   double m_frootc_to_fire;     /*  fine root C fire loss (gC/m2/s) */
   double m_frootc_storage_to_fire;     /*  fine root C storage fire loss (gC/m2/s) */
   double m_frootc_xfer_to_fire;        /*  fine root C transfer fire loss (gC/m2/s) */
   double m_livestemc_to_fire;      /*  live stem C fire loss (gC/m2/s) */
   double m_livestemc_storage_to_fire;      /*  live stem C storage fire loss (gC/m2/s) */
   double m_livestemc_xfer_to_fire;     /*  live stem C transfer fire loss (gC/m2/s) */
   double m_deadstemc_to_fire;      /*  dead stem C fire loss (gC/m2/s) */
   double m_deadstemc_to_litter_fire;       /*  dead stem C fire mortality to litter (gC/m2/s) */
   double m_deadstemc_storage_to_fire;      /*  dead stem C storage fire loss (gC/m2/s) */
   double m_deadstemc_xfer_to_fire;     /*  dead stem C transfer fire loss (gC/m2/s) */
   double m_livecrootc_to_fire;     /*  live coarse root C fire loss (gC/m2/s) */
   double m_livecrootc_storage_to_fire;     /*  live coarse root C storage fire loss (gC/m2/s) */
   double m_livecrootc_xfer_to_fire;        /*  live coarse root C transfer fire loss (gC/m2/s) */
   double m_deadcrootc_to_fire;     /*  dead coarse root C fire loss (gC/m2/s) */
   double m_deadcrootc_to_litter_fire;      /*  dead coarse root C fire mortality to litter (gC/m2/s) */
   double m_deadcrootc_storage_to_fire;     /*  dead coarse root C storage fire loss (gC/m2/s) */
   double m_deadcrootc_xfer_to_fire;        /*  dead coarse root C transfer fire loss (gC/m2/s) */
   double m_gresp_storage_to_fire;      /*  growth respiration storage fire loss (gC/m2/s) */
   double m_gresp_xfer_to_fire;     /*  growth respiration transfer fire loss (gC/m2/s) */
   /* phenology fluxes from transfer pools */
   double grainc_xfer_to_grainc;        /*  grain C growth from storage for prognostic crop(gC/m2/s) */
   double leafc_xfer_to_leafc;      /*  leaf C growth from storage (gC/m2/s) */
   double frootc_xfer_to_frootc;        /*  fine root C growth from storage (gC/m2/s) */
   double livestemc_xfer_to_livestemc;      /*  live stem C growth from storage (gC/m2/s) */
   double deadstemc_xfer_to_deadstemc;      /*  dead stem C growth from storage (gC/m2/s) */
   double livecrootc_xfer_to_livecrootc;        /*  live coarse root C growth from storage (gC/m2/s) */
   double deadcrootc_xfer_to_deadcrootc;        /*  dead coarse root C growth from storage (gC/m2/s) */
   /* leaf and fine root litterfall */
   double leafc_to_litter;      /*  leaf C litterfall (gC/m2/s) */
   double frootc_to_litter;     /*  fine root C litterfall (gC/m2/s) */
   double livestemc_to_litter;      /*  live stem C litterfall (gC/m2/s) */
   double grainc_to_food;       /*  grain C to food for prognostic crop(gC/m2/s) */
   /* maintenance respiration fluxes */
   double leaf_mr;      /*  leaf maintenance respiration (gC/m2/s) */
   double froot_mr;     /*  fine root maintenance respiration (gC/m2/s) */
   double livestem_mr;      /*  live stem maintenance respiration (gC/m2/s) */
   double livecroot_mr;     /*  live coarse root maintenance respiration (gC/m2/s) */
   double leaf_curmr;       /*  leaf maintenance respiration from current GPP (gC/m2/s) */
   double froot_curmr;      /*  fine root maintenance respiration from current GPP (gC/m2/s) */
   double livestem_curmr;       /*  live stem maintenance respiration from current GPP (gC/m2/s) */
   double livecroot_curmr;      /*  live coarse root maintenance respiration from current GPP (gC/m2/s) */
   double leaf_xsmr;        /*  leaf maintenance respiration from storage (gC/m2/s) */
   double froot_xsmr;       /*  fine root maintenance respiration from storage (gC/m2/s) */
   double livestem_xsmr;        /*  live stem maintenance respiration from storage (gC/m2/s) */
   double livecroot_xsmr;       /*  live coarse root maintenance respiration from storage (gC/m2/s) */
   /* photosynthesis fluxes */
   double psnsun_to_cpool;      /*  C fixation from sunlit canopy (gC/m2/s) */
   double psnshade_to_cpool;        /*  C fixation from shaded canopy (gC/m2/s) */
   /* allocation fluxes, from current GPP */
   double cpool_to_xsmrpool;        /*  allocation to maintenance respiration storage pool (gC/m2/s) */
   double cpool_to_grainc;      /*  allocation to grain C for prognostic crop(gC/m2/s) */
   double cpool_to_grainc_storage;      /*  allocation to grain C storage for prognostic crop(gC/m2/s) */
   double cpool_to_leafc;       /*  allocation to leaf C (gC/m2/s) */
   double cpool_to_leafc_storage;       /*  allocation to leaf C storage (gC/m2/s) */
   double cpool_to_frootc;      /*  allocation to fine root C (gC/m2/s) */
   double cpool_to_frootc_storage;      /*  allocation to fine root C storage (gC/m2/s) */
   double cpool_to_livestemc;       /*  allocation to live stem C (gC/m2/s) */
   double cpool_to_livestemc_storage;       /*  allocation to live stem C storage (gC/m2/s) */
   double cpool_to_deadstemc;       /*  allocation to dead stem C (gC/m2/s) */
   double cpool_to_deadstemc_storage;       /*  allocation to dead stem C storage (gC/m2/s) */
   double cpool_to_livecrootc;      /*  allocation to live coarse root C (gC/m2/s) */
   double cpool_to_livecrootc_storage;      /*  allocation to live coarse root C storage (gC/m2/s) */
   double cpool_to_deadcrootc;      /*  allocation to dead coarse root C (gC/m2/s) */
   double cpool_to_deadcrootc_storage;      /*  allocation to dead coarse root C storage (gC/m2/s) */
   double cpool_to_gresp_storage;       /*  allocation to growth respiration storage (gC/m2/s) */
   /* growth respiration fluxes */
   double xsmrpool_to_atm;      /*  excess MR pool harvest mortality (gC/m2/s) */
   double cpool_leaf_gr;        /*  leaf growth respiration (gC/m2/s) */
   double cpool_leaf_storage_gr;        /*  leaf growth respiration to storage (gC/m2/s) */
   double transfer_leaf_gr;     /*  leaf growth respiration from storage (gC/m2/s) */
   double cpool_froot_gr;       /*  fine root growth respiration (gC/m2/s) */
   double cpool_froot_storage_gr;       /*  fine root  growth respiration to storage (gC/m2/s) */
   double transfer_froot_gr;        /*  fine root  growth respiration from storage (gC/m2/s) */
   double cpool_livestem_gr;        /*  live stem growth respiration (gC/m2/s) */
   double cpool_livestem_storage_gr;        /*  live stem growth respiration to storage (gC/m2/s) */
   double transfer_livestem_gr;     /*  live stem growth respiration from storage (gC/m2/s) */
   double cpool_deadstem_gr;        /*  dead stem growth respiration (gC/m2/s) */
   double cpool_deadstem_storage_gr;        /*  dead stem growth respiration to storage (gC/m2/s) */
   double transfer_deadstem_gr;     /*  dead stem growth respiration from storage (gC/m2/s) */
   double cpool_livecroot_gr;       /*  live coarse root growth respiration (gC/m2/s) */
   double cpool_livecroot_storage_gr;       /*  live coarse root growth respiration to storage (gC/m2/s) */
   double transfer_livecroot_gr;        /*  live coarse root growth respiration from storage (gC/m2/s) */
   double cpool_deadcroot_gr;       /*  dead coarse root growth respiration (gC/m2/s) */
   double cpool_deadcroot_storage_gr;       /*  dead coarse root growth respiration to storage (gC/m2/s) */
   double transfer_deadcroot_gr;        /*  dead coarse root growth respiration from storage (gC/m2/s) */
   /* growth respiration for prognostic crop model */
   double cpool_grain_gr;       /*  grain growth respiration (gC/m2/s) */
   double cpool_grain_storage_gr;       /*  grain growth respiration to storage (gC/m2/s) */
   double transfer_grain_gr;        /*  grain growth respiration from storage (gC/m2/s) */
   /* annual turnover of storage to transfer pools */
   double grainc_storage_to_xfer;       /*  grain C shift storage to transfer for prognostic crop model (gC/m2/s) */
   double leafc_storage_to_xfer;        /*  leaf C shift storage to transfer (gC/m2/s) */
   double frootc_storage_to_xfer;       /*  fine root C shift storage to transfer (gC/m2/s) */
   double livestemc_storage_to_xfer;        /*  live stem C shift storage to transfer (gC/m2/s) */
   double deadstemc_storage_to_xfer;        /*  dead stem C shift storage to transfer (gC/m2/s) */
   double livecrootc_storage_to_xfer;       /*  live coarse root C shift storage to transfer (gC/m2/s) */
   double deadcrootc_storage_to_xfer;       /*  dead coarse root C shift storage to transfer (gC/m2/s) */
   double gresp_storage_to_xfer;        /*  growth respiration shift storage to transfer (gC/m2/s) */
   /* turnover of livewood to deadwood */
   double livestemc_to_deadstemc;       /*  live stem C turnover (gC/m2/s) */
   double livecrootc_to_deadcrootc;     /*  live coarse root C turnover (gC/m2/s) */
   /* summary (diagnostic) flux variables, not involved in mass balance */
   double gpp;      /*  (gC/m2/s) gross primary production  */
   double mr;       /*  (gC/m2/s) maintenance respiration */
   double current_gr;       /*  (gC/m2/s) growth resp for new growth displayed in this timestep */
   double transfer_gr;      /*  (gC/m2/s) growth resp for transfer growth displayed in this timestep */
   double storage_gr;       /*  (gC/m2/s) growth resp for growth sent to storage for later display */
   double gr;       /*  (gC/m2/s) total growth respiration */
   double ar;       /*  (gC/m2/s) autotrophic respiration (MR + GR) */
   double rr;       /*  (gC/m2/s) root respiration (fine root MR + total root GR) */
   double npp;      /*  (gC/m2/s) net primary production */
   double agnpp;        /*  (gC/m2/s) aboveground NPP */
   double bgnpp;        /*  (gC/m2/s) belowground NPP */
   double litfall;      /*  (gC/m2/s) litterfall (leaves and fine roots) */
   double vegfire;      /*  (gC/m2/s) pft-level fire loss (obsolete, mark for removal) */
   double wood_harvestc;        /*  (gC/m2/s) pft-level wood harvest (to product pools) */
   double pft_cinputs;      /*  (gC/m2/s) pft-level carbon inputs (for balance checking) */
   double pft_coutputs;     /*  (gC/m2/s) pft-level carbon outputs (for balance checking) */
   /* CLAMP summary (diagnostic) variables, not involved in mass balance */
   double frootc_alloc;     /*  (gC/m2/s) pft-level fine root C alloc */
   double frootc_loss;      /*  (gC/m2/s) pft-level fine root C loss */
   double leafc_alloc;      /*  (gC/m2/s) pft-level leaf C alloc */
   double leafc_loss;       /*  (gC/m2/s) pft-level leaf C loss */
   double woodc_alloc;      /*  (gC/m2/s) pft-level wood C alloc */
   double woodc_loss;       /*  (gC/m2/s) pft-level wood C loss */
   /* new variables for fire code */
   double pft_fire_closs;       /*  (gC/m2/s) total pft-level fire C loss  */

   double m_leafc_to_litr1c;        /*  leaf C mortality to litter 1 C (gC/m2/s)  */
   double m_leafc_to_litr2c;        /*  leaf C mortality to litter 2 C (gC/m2/s) */
   double m_leafc_to_litr3c;        /*  leaf C mortality to litter 3 C (gC/m2/s) */
   double m_frootc_to_litr1c;       /*  fine root C mortality to litter 1 C (gC/m2/s) */
   double m_frootc_to_litr2c;       /*  fine root C mortality to litter 2 C (gC/m2/s) */
   double m_frootc_to_litr3c;       /*  fine root C mortality to litter 3 C (gC/m2/s) */
   double m_livestemc_to_cwdc;      /*  live stem C mortality to coarse woody debris C (gC/m2/s) */
   double m_deadstemc_to_cwdc;      /*  dead stem C mortality to coarse woody debris C (gC/m2/s) */
   double m_livecrootc_to_cwdc;     /*  live coarse root C mortality to coarse woody debris C (gC/m2/s) */
   double m_deadcrootc_to_cwdc;     /*  dead coarse root C mortality to coarse woody debris C (gC/m2/s) */
   double m_leafc_storage_to_litr1c;        /*  leaf C storage mortality to litter 1 C (gC/m2/s) */
   double m_frootc_storage_to_litr1c;       /*  fine root C storage mortality to litter 1 C (gC/m2/s) */
   double m_livestemc_storage_to_litr1c;        /*  live stem C storage mortality to litter 1 C (gC/m2/s) */
   double m_deadstemc_storage_to_litr1c;        /*  dead stem C storage mortality to litter 1 C (gC/m2/s) */
   double m_livecrootc_storage_to_litr1c;       /*  live coarse root C storage mortality to litter 1 C (gC/m2/s) */
   double m_deadcrootc_storage_to_litr1c;       /*  dead coarse root C storage mortality to litter 1 C (gC/m2/s) */
   double m_gresp_storage_to_litr1c;        /*  growth respiration storage mortality to litter 1 C (gC/m2/s) */
   double m_leafc_xfer_to_litr1c;       /*  leaf C transfer mortality to litter 1 C (gC/m2/s) */
   double m_frootc_xfer_to_litr1c;      /*  fine root C transfer mortality to litter 1 C (gC/m2/s) */
   double m_livestemc_xfer_to_litr1c;       /*  live stem C transfer mortality to litter 1 C (gC/m2/s) */
   double m_deadstemc_xfer_to_litr1c;       /*  dead stem C transfer mortality to litter 1 C (gC/m2/s) */
   double m_livecrootc_xfer_to_litr1c;      /*  live coarse root C transfer mortality to litter 1 C (gC/m2/s) */
   double m_deadcrootc_xfer_to_litr1c;      /*  dead coarse root C transfer mortality to litter 1 C (gC/m2/s) */
   double m_gresp_xfer_to_litr1c;       /*  growth respiration transfer mortality to litter 1 C (gC/m2/s) */
   double hrv_leafc_to_litr1c;      /*  leaf C harvest mortality to litter 1 C (gC/m2/s)                          */
   double hrv_leafc_to_litr2c;      /*  leaf C harvest mortality to litter 2 C (gC/m2/s)                         */
   double hrv_leafc_to_litr3c;      /*  leaf C harvest mortality to litter 3 C (gC/m2/s)                         */
   double hrv_frootc_to_litr1c;     /*  fine root C harvest mortality to litter 1 C (gC/m2/s)                    */
   double hrv_frootc_to_litr2c;     /*  fine root C harvest mortality to litter 2 C (gC/m2/s)                    */
   double hrv_frootc_to_litr3c;     /*  fine root C harvest mortality to litter 3 C (gC/m2/s)                    */
   double hrv_livestemc_to_cwdc;        /*  live stem C harvest mortality to coarse woody debris C (gC/m2/s)         */
   double hrv_deadstemc_to_prod10c;     /*  dead stem C harvest mortality to 10-year product pool (gC/m2/s)         */
   double hrv_deadstemc_to_prod100c;        /*  dead stem C harvest mortality to 100-year product pool (gC/m2/s)         */
   double hrv_livecrootc_to_cwdc;       /*  live coarse root C harvest mortality to coarse woody debris C (gC/m2/s)  */
   double hrv_deadcrootc_to_cwdc;       /*  dead coarse root C harvest mortality to coarse woody debris C (gC/m2/s)  */
   double hrv_leafc_storage_to_litr1c;      /*  leaf C storage harvest mortality to litter 1 C (gC/m2/s)                 */
   double hrv_frootc_storage_to_litr1c;     /*  fine root C storage harvest mortality to litter 1 C (gC/m2/s)            */
   double hrv_livestemc_storage_to_litr1c;      /*  live stem C storage harvest mortality to litter 1 C (gC/m2/s)            */
   double hrv_deadstemc_storage_to_litr1c;      /*  dead stem C storage harvest mortality to litter 1 C (gC/m2/s)            */
   double hrv_livecrootc_storage_to_litr1c;     /*  live coarse root C storage harvest mortality to litter 1 C (gC/m2/s)     */
   double hrv_deadcrootc_storage_to_litr1c;     /*  dead coarse root C storage harvest mortality to litter 1 C (gC/m2/s)     */
   double hrv_gresp_storage_to_litr1c;      /*  growth respiration storage harvest mortality to litter 1 C (gC/m2/s)     */
   double hrv_leafc_xfer_to_litr1c;     /*  leaf C transfer harvest mortality to litter 1 C (gC/m2/s)                */
   double hrv_frootc_xfer_to_litr1c;        /*  fine root C transfer harvest mortality to litter 1 C (gC/m2/s)           */
   double hrv_livestemc_xfer_to_litr1c;     /*  live stem C transfer harvest mortality to litter 1 C (gC/m2/s)           */
   double hrv_deadstemc_xfer_to_litr1c;     /*  dead stem C transfer harvest mortality to litter 1 C (gC/m2/s)           */
   double hrv_livecrootc_xfer_to_litr1c;        /*  live coarse root C transfer harvest mortality to litter 1 C (gC/m2/s)    */
   double hrv_deadcrootc_xfer_to_litr1c;        /*  dead coarse root C transfer harvest mortality to litter 1 C (gC/m2/s)    */
   double hrv_gresp_xfer_to_litr1c;     /*  growth respiration transfer harvest mortality to litter 1 C (gC/m2/s)    */
   /* column-level fire fluxes */
   double m_deadstemc_to_cwdc_fire;     /*  dead stem C to coarse woody debris C by fire (gC/m2/s) */
   double m_deadcrootc_to_cwdc_fire;        /*  dead coarse root C to to woody debris C by fire (gC/m2/s) */
   double m_litr1c_to_fire;     /*  litter 1 C fire loss (gC/m2/s) */
   double m_litr2c_to_fire;     /*  litter 2 C fire loss (gC/m2/s) */
   double m_litr3c_to_fire;     /*  litter 3 C fire loss (gC/m2/s) */
   double m_cwdc_to_fire;       /*  coarse woody debris C fire loss (gC/m2/s) */
   /* litterfall fluxes */
   double livestemc_to_litr1c;      /*  livestem C litterfall to litter 1 C (gC/m2/s) */
   double livestemc_to_litr2c;      /*  livestem C litterfall to litter 2 C (gC/m2/s) */
   double livestemc_to_litr3c;      /*  livestem C litterfall to litter 3 C (gC/m2/s) */
   double leafc_to_litr1c;      /*  leaf C litterfall to litter 1 C (gC/m2/s) */
   double leafc_to_litr2c;      /*  leaf C litterfall to litter 2 C (gC/m2/s) */
   double leafc_to_litr3c;      /*  leaf C litterfall to litter 3 C (gC/m2/s) */
   double frootc_to_litr1c;     /*  fine root C litterfall to litter 1 C (gC/m2/s) */
   double frootc_to_litr2c;     /*  fine root C litterfall to litter 2 C (gC/m2/s) */
   double frootc_to_litr3c;     /*  fine root C litterfall to litter 3 C (gC/m2/s) */
   /* litterfall fluxes for prognostic crop model */
   double grainc_to_litr1c;     /*  grain C litterfall to litter 1 C (gC/m2/s) */
   double grainc_to_litr2c;     /*  grain C litterfall to litter 2 C (gC/m2/s) */
   double grainc_to_litr3c;     /*  grain C litterfall to litter 3 C (gC/m2/s) */
   /* decomposition fluxes */
   double cwdc_to_litr2c;       /*  decomp. of coarse woody debris C to litter 2 C (gC/m2/s) */
   double cwdc_to_litr3c;       /*  decomp. of coarse woody debris C to litter 3 C (gC/m2/s) */
   double litr1_hr;     /*  het. resp. from litter 1 C (gC/m2/s) */
   double litr1c_to_soil1c;     /*  decomp. of litter 1 C to SOM 1 C (gC/m2/s) */
   double litr2_hr;     /*  het. resp. from litter 2 C (gC/m2/s) */
   double litr2c_to_soil2c;     /*  decomp. of litter 2 C to SOM 2 C (gC/m2/s) */
   double litr3_hr;     /*  het. resp. from litter 3 C (gC/m2/s) */
   double litr3c_to_soil3c;     /*  decomp. of litter 3 C to SOM 3 C (gC/m2/s) */
   double soil1_hr;     /*  het. resp. from SOM 1 C (gC/m2/s) */
   double soil1c_to_soil2c;     /*  decomp. of SOM 1 C to SOM 2 C (gC/m2/s) */
   double soil2_hr;     /*  het. resp. from SOM 2 C (gC/m2/s) */
   double soil2c_to_soil3c;     /*  decomp. of SOM 2 C to SOM 3 C (gC/m2/s) */
   double soil3_hr;     /*  het. resp. from SOM 3 C (gC/m2/s) */
   double soil3c_to_soil4c;     /*  decomp. of SOM 3 C to SOM 4 C (gC/m2/s) */
   double soil4_hr;     /*  het. resp. from SOM 4 C (gC/m2/s) */
   /* dynamic landcover fluxes */
   double dwt_seedc_to_leaf;        /*  (gC/m2/s) seed source to PFT-level */
   double dwt_seedc_to_deadstem;        /*  (gC/m2/s) seed source to PFT-level */
   double dwt_conv_cflux;       /*  (gC/m2/s) conversion C flux (immediate loss to atm) */
   double dwt_prod10c_gain;     /*  (gC/m2/s) addition to 10-yr wood product pool */
   double dwt_prod100c_gain;        /*  (gC/m2/s) addition to 100-yr wood product pool */
   double dwt_frootc_to_litr1c;     /*  (gC/m2/s) fine root to litter due to landcover change */
   double dwt_frootc_to_litr2c;     /*  (gC/m2/s) fine root to litter due to landcover change */
   double dwt_frootc_to_litr3c;     /*  (gC/m2/s) fine root to litter due to landcover change */
   double dwt_livecrootc_to_cwdc;       /*  (gC/m2/s) live coarse root to CWD due to landcover change */
   double dwt_deadcrootc_to_cwdc;       /*  (gC/m2/s) dead coarse root to CWD due to landcover change */
   double dwt_closs;        /*  (gC/m2/s) total carbon loss from product pools and conversion */
   double landuseflux;      /*  (gC/m2/s) dwt_closs+product_closs */
   double landuptake;       /*  (gC/m2/s) nee-landuseflux */
   /* wood product pool loss fluxes */
   double prod10c_loss;     /*  (gC/m2/s) decomposition loss from 10-yr wood product pool */
   double prod100c_loss;        /*  (gC/m2/s) decomposition loss from 100-yr wood product pool */
   double product_closs;        /*  (gC/m2/s) total wood product carbon loss */
   /* summary (diagnostic) flux variables, not involved in mass balance */
   double lithr;        /*  (gC/m2/s) litter heterotrophic respiration  */
   double somhr;        /*  (gC/m2/s) soil organic matter heterotrophic respiration */
   double hr;       /*  (gC/m2/s) total heterotrophic respiration */
   double sr;       /*  (gC/m2/s) total soil respiration (HR + root resp) */
   double er;       /*  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic */
   double litfire;      /*  (gC/m2/s) litter fire losses */
   double somfire;      /*  (gC/m2/s) soil organic matter fire losses */
   double totfire;      /*  (gC/m2/s) total ecosystem fire losses */
   double nep;      /*  (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink */
   double nbp;      /*  (gC/m2/s) net biome production, includes fire, landuse, and harvest flux, positive for sink */
   double nee;      /*  (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source */
   double col_cinputs;      /*  (gC/m2/s) total column-level carbon inputs (for balance check) */
   double col_coutputs;     /*  (gC/m2/s) total column-level carbon outputs (for balance check)  */
   /* CLAMP summary (diagnostic) flux variables, not involved in mass balance */
   double cwdc_hr;      /*  (gC/m2/s) col-level coarse woody debris C heterotrophic respiration */
   double cwdc_loss;        /*  (gC/m2/s) col-level coarse woody debris C loss */
   double litterc_loss;     /*  (gC/m2/s) col-level litter C loss */
   /*  new variables for fire */
   double col_fire_closs;       /*  (gC/m2/s) total column-level fire C loss */
} cflux_type;


/*
 * pft nitrogen flux variables structure
 */
typedef struct nflux_type
{
    /* new variables for CN code */
    /* gap mortality fluxes */
   double m_leafn_to_litter;        /*  leaf N mortality (gN/m2/s) */
   double m_frootn_to_litter;       /*  fine root N mortality (gN/m2/s) */
   double m_leafn_storage_to_litter;        /*  leaf N storage mortality (gN/m2/s) */
   double m_frootn_storage_to_litter;       /*  fine root N storage mortality (gN/m2/s) */
   double m_livestemn_storage_to_litter;        /*  live stem N storage mortality (gN/m2/s) */
   double m_deadstemn_storage_to_litter;        /*  dead stem N storage mortality (gN/m2/s) */
   double m_livecrootn_storage_to_litter;       /*  live coarse root N storage mortality (gN/m2/s) */
   double m_deadcrootn_storage_to_litter;       /*  dead coarse root N storage mortality (gN/m2/s) */
   double m_leafn_xfer_to_litter;       /*  leaf N transfer mortality (gN/m2/s) */
   double m_frootn_xfer_to_litter;      /*  fine root N transfer mortality (gN/m2/s) */
   double m_livestemn_xfer_to_litter;       /*  live stem N transfer mortality (gN/m2/s) */
   double m_deadstemn_xfer_to_litter;       /*  dead stem N transfer mortality (gN/m2/s) */
   double m_livecrootn_xfer_to_litter;      /*  live coarse root N transfer mortality (gN/m2/s) */
   double m_deadcrootn_xfer_to_litter;      /*  dead coarse root N transfer mortality (gN/m2/s) */
   double m_livestemn_to_litter;        /*  live stem N mortality (gN/m2/s) */
   double m_deadstemn_to_litter;        /*  dead stem N mortality (gN/m2/s) */
   double m_livecrootn_to_litter;       /*  live coarse root N mortality (gN/m2/s) */
   double m_deadcrootn_to_litter;       /*  dead coarse root N mortality (gN/m2/s) */
   double m_retransn_to_litter;     /*  retranslocated N pool mortality (gN/m2/s) */
   /* harvest mortality fluxes */
   double hrv_leafn_to_litter;      /*  leaf N harvest mortality (gN/m2/s) */
   double hrv_frootn_to_litter;     /*  fine root N harvest mortality (gN/m2/s) */
   double hrv_leafn_storage_to_litter;      /*  leaf N storage harvest mortality (gN/m2/s) */
   double hrv_frootn_storage_to_litter;     /*  fine root N storage harvest mortality (gN/m2/s) */
   double hrv_livestemn_storage_to_litter;      /*  live stem N storage harvest mortality (gN/m2/s) */
   double hrv_deadstemn_storage_to_litter;      /*  dead stem N storage harvest mortality (gN/m2/s) */
   double hrv_livecrootn_storage_to_litter;     /*  live coarse root N storage harvest mortality (gN/m2/s) */
   double hrv_deadcrootn_storage_to_litter;     /*  dead coarse root N storage harvest mortality (gN/m2/s) */
   double hrv_leafn_xfer_to_litter;     /*  leaf N transfer harvest mortality (gN/m2/s) */
   double hrv_frootn_xfer_to_litter;        /*  fine root N transfer harvest mortality (gN/m2/s) */
   double hrv_livestemn_xfer_to_litter;     /*  live stem N transfer harvest mortality (gN/m2/s) */
   double hrv_deadstemn_xfer_to_litter;     /*  dead stem N transfer harvest mortality (gN/m2/s) */
   double hrv_livecrootn_xfer_to_litter;        /*  live coarse root N transfer harvest mortality (gN/m2/s) */
   double hrv_deadcrootn_xfer_to_litter;        /*  dead coarse root N transfer harvest mortality (gN/m2/s) */
   double hrv_livestemn_to_litter;      /*  live stem N harvest mortality (gN/m2/s) */
   double hrv_deadstemn_to_prod10n;     /*  dead stem N harvest to 10-year product pool (gN/m2/s) */
   double hrv_deadstemn_to_prod100n;        /*  dead stem N harvest to 100-year product pool (gN/m2/s) */
   double hrv_livecrootn_to_litter;     /*  live coarse root N harvest mortality (gN/m2/s) */
   double hrv_deadcrootn_to_litter;     /*  dead coarse root N harvest mortality (gN/m2/s) */
   double hrv_retransn_to_litter;       /*  retranslocated N pool harvest mortality (gN/m2/s) */
   /* fire mortality fluxes */
   double m_leafn_to_fire;      /*  leaf N fire loss (gN/m2/s) */
   double m_leafn_storage_to_fire;      /*  leaf N storage fire loss (gN/m2/s) */
   double m_leafn_xfer_to_fire;     /*  leaf N transfer fire loss (gN/m2/s) */
   double m_frootn_to_fire;     /*  fine root N fire loss (gN/m2/s) */
   double m_frootn_storage_to_fire;     /*  fine root N storage fire loss (gN/m2/s) */
   double m_frootn_xfer_to_fire;        /*  fine root N transfer fire loss (gN/m2/s) */
   double m_livestemn_to_fire;      /*  live stem N fire loss (gN/m2/s) */
   double m_livestemn_storage_to_fire;      /*  live stem N storage fire loss (gN/m2/s) */
   double m_livestemn_xfer_to_fire;     /*  live stem N transfer fire loss (gN/m2/s) */
   double m_deadstemn_to_fire;      /*  dead stem N fire loss (gN/m2/s) */
   double m_deadstemn_to_litter_fire;       /*  dead stem N fire mortality to litter (gN/m2/s) */
   double m_deadstemn_storage_to_fire;      /*  dead stem N storage fire loss (gN/m2/s) */
   double m_deadstemn_xfer_to_fire;     /*  dead stem N transfer fire loss (gN/m2/s) */
   double m_livecrootn_to_fire;     /*  live coarse root N fire loss (gN/m2/s) */
   double m_livecrootn_storage_to_fire;     /*  live coarse root N storage fire loss (gN/m2/s) */
   double m_livecrootn_xfer_to_fire;        /*  live coarse root N transfer fire loss (gN/m2/s) */
   double m_deadcrootn_to_fire;     /*  dead coarse root N fire loss (gN/m2/s) */
   double m_deadcrootn_to_litter_fire;      /*  dead coarse root N fire mortality to litter (gN/m2/s) */
   double m_deadcrootn_storage_to_fire;     /*  dead coarse root N storage fire loss (gN/m2/s) */
   double m_deadcrootn_xfer_to_fire;        /*  dead coarse root N transfer fire loss (gN/m2/s) */
   double m_retransn_to_fire;       /*  retranslocated N pool fire loss (gN/m2/s) */
   /* phenology fluxes from transfer pool */
   double grainn_xfer_to_grainn;        /*  grain N growth from storage for prognostic crop model (gN/m2/s) */
   double leafn_xfer_to_leafn;      /*  leaf N growth from storage (gN/m2/s) */
   double frootn_xfer_to_frootn;        /*  fine root N growth from storage (gN/m2/s) */
   double livestemn_xfer_to_livestemn;      /*  live stem N growth from storage (gN/m2/s) */
   double deadstemn_xfer_to_deadstemn;      /*  dead stem N growth from storage (gN/m2/s) */
   double livecrootn_xfer_to_livecrootn;        /*  live coarse root N growth from storage (gN/m2/s) */
   double deadcrootn_xfer_to_deadcrootn;        /*  dead coarse root N growth from storage (gN/m2/s) */
   /* litterfall fluxes */
   double livestemn_to_litter;      /*  livestem N to litter (gN/m2/s) */
   double grainn_to_food;       /*  grain N to food for prognostic crop (gN/m2/s) */
   double leafn_to_litter;      /*  leaf N litterfall (gN/m2/s) */
   double leafn_to_retransn;        /*  leaf N to retranslocated N pool (gN/m2/s) */
   double frootn_to_litter;     /*  fine root N litterfall (gN/m2/s) */
   /* allocation fluxes */
   double retransn_to_npool;        /*  deployment of retranslocated N (gN/m2/s)        */
   double sminn_to_npool;       /*  deployment of soil mineral N uptake (gN/m2/s) */
   double npool_to_grainn;      /*  allocation to grain N for prognostic crop (gN/m2/s) */
   double npool_to_grainn_storage;      /*  allocation to grain N storage for prognostic crop (gN/m2/s) */
   double npool_to_leafn;       /*  allocation to leaf N (gN/m2/s) */
   double npool_to_leafn_storage;       /*  allocation to leaf N storage (gN/m2/s) */
   double npool_to_frootn;      /*  allocation to fine root N (gN/m2/s) */
   double npool_to_frootn_storage;      /*  allocation to fine root N storage (gN/m2/s) */
   double npool_to_livestemn;       /*  allocation to live stem N (gN/m2/s) */
   double npool_to_livestemn_storage;       /*  allocation to live stem N storage (gN/m2/s) */
   double npool_to_deadstemn;       /*  allocation to dead stem N (gN/m2/s) */
   double npool_to_deadstemn_storage;       /*  allocation to dead stem N storage (gN/m2/s) */
   double npool_to_livecrootn;      /*  allocation to live coarse root N (gN/m2/s) */
   double npool_to_livecrootn_storage;      /*  allocation to live coarse root N storage (gN/m2/s) */
   double npool_to_deadcrootn;      /*  allocation to dead coarse root N (gN/m2/s) */
   double npool_to_deadcrootn_storage;      /*  allocation to dead coarse root N storage (gN/m2/s) */
   /* annual turnover of storage to transfer pools */
   double grainn_storage_to_xfer;       /*  grain N shift storage to transfer for prognostic crop (gN/m2/s) */
   double leafn_storage_to_xfer;        /*  leaf N shift storage to transfer (gN/m2/s) */
   double frootn_storage_to_xfer;       /*  fine root N shift storage to transfer (gN/m2/s) */
   double livestemn_storage_to_xfer;        /*  live stem N shift storage to transfer (gN/m2/s) */
   double deadstemn_storage_to_xfer;        /*  dead stem N shift storage to transfer (gN/m2/s) */
   double livecrootn_storage_to_xfer;       /*  live coarse root N shift storage to transfer (gN/m2/s) */
   double deadcrootn_storage_to_xfer;       /*  dead coarse root N shift storage to transfer (gN/m2/s) */
   /* turnover of livewood to deadwood, with retranslocation  */
   double livestemn_to_deadstemn;       /*  live stem N turnover (gN/m2/s) */
   double livestemn_to_retransn;        /*  live stem N to retranslocated N pool (gN/m2/s) */
   double livecrootn_to_deadcrootn;     /*  live coarse root N turnover (gN/m2/s) */
   double livecrootn_to_retransn;       /*  live coarse root N to retranslocated N pool (gN/m2/s) */
   /* summary (diagnostic) flux variables, not involved in mass balance */
   double ndeploy;      /*  total N deployed to growth and storage (gN/m2/s) */
   double pft_ninputs;      /*  total N inputs to pft-level (gN/m2/s) */
   double pft_noutputs;     /*  total N outputs from pft-level (gN/m2/s) */
   double wood_harvestn;        /*  total N losses to wood product pools (gN/m2/s) */
   /* new variables for fire code */
   double pft_fire_nloss;       /*  total pft-level fire N loss (gN/m2/s)  */

    /*  deposition fluxes */
   double ndep_to_sminn;        /*  atmospheric N deposition to soil mineral N (gN/m2/s) */
   double nfix_to_sminn;        /*  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)  */
   /* column-level gap mortality fluxes */
   double m_leafn_to_litr1n;        /*  leaf N mortality to litter 1 N (gC/m2/s) */
   double m_leafn_to_litr2n;        /*  leaf N mortality to litter 2 N (gC/m2/s) */
   double m_leafn_to_litr3n;        /*  leaf N mortality to litter 3 N (gC/m2/s) */
   double m_frootn_to_litr1n;       /*  fine root N mortality to litter 1 N (gN/m2/s) */
   double m_frootn_to_litr2n;       /*  fine root N mortality to litter 2 N (gN/m2/s) */
   double m_frootn_to_litr3n;       /*  fine root N mortality to litter 3 N (gN/m2/s) */
   double m_livestemn_to_cwdn;      /*  live stem N mortality to coarse woody debris N (gN/m2/s) */
   double m_deadstemn_to_cwdn;      /*  dead stem N mortality to coarse woody debris N (gN/m2/s) */
   double m_livecrootn_to_cwdn;     /*  live coarse root N mortality to coarse woody debris N (gN/m2/s) */
   double m_deadcrootn_to_cwdn;     /*  dead coarse root N mortality to coarse woody debris N (gN/m2/s) */
   double m_retransn_to_litr1n;     /*  retranslocated N pool mortality to litter 1 N (gN/m2/s) */
   double m_leafn_storage_to_litr1n;        /*  leaf N storage mortality to litter 1 N (gN/m2/s) */
   double m_frootn_storage_to_litr1n;       /*  fine root N storage mortality to litter 1 N (gN/m2/s) */
   double m_livestemn_storage_to_litr1n;        /*  live stem N storage mortality to litter 1 N (gN/m2/s) */
   double m_deadstemn_storage_to_litr1n;        /*  dead stem N storage mortality to litter 1 N (gN/m2/s) */
   double m_livecrootn_storage_to_litr1n;       /*  live coarse root N storage mortality to litter 1 N (gN/m2/s) */
   double m_deadcrootn_storage_to_litr1n;       /*  dead coarse root N storage mortality to litter 1 N (gN/m2/s) */
   double m_leafn_xfer_to_litr1n;       /*  leaf N transfer mortality to litter 1 N (gN/m2/s) */
   double m_frootn_xfer_to_litr1n;      /*  fine root N transfer mortality to litter 1 N (gN/m2/s) */
   double m_livestemn_xfer_to_litr1n;       /*  live stem N transfer mortality to litter 1 N (gN/m2/s) */
   double m_deadstemn_xfer_to_litr1n;       /*  dead stem N transfer mortality to litter 1 N (gN/m2/s) */
   double m_livecrootn_xfer_to_litr1n;      /*  live coarse root N transfer mortality to litter 1 N (gN/m2/s) */
   double m_deadcrootn_xfer_to_litr1n;      /*  dead coarse root N transfer mortality to litter 1 N (gN/m2/s) */
   /* column-level harvest fluxes */
   double hrv_leafn_to_litr1n;      /*  leaf N harvest mortality to litter 1 N (gC/m2/s) */
   double hrv_leafn_to_litr2n;      /*  leaf N harvest mortality to litter 2 N (gC/m2/s) */
   double hrv_leafn_to_litr3n;      /*  leaf N harvest mortality to litter 3 N (gC/m2/s) */
   double hrv_frootn_to_litr1n;     /*  fine root N harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_frootn_to_litr2n;     /*  fine root N harvest mortality to litter 2 N (gN/m2/s) */
   double hrv_frootn_to_litr3n;     /*  fine root N harvest mortality to litter 3 N (gN/m2/s) */
   double hrv_livestemn_to_cwdn;        /*  live stem N harvest mortality to coarse woody debris N (gN/m2/s) */
   double hrv_deadstemn_to_prod10n;     /*  dead stem N harvest mortality to 10-year product pool (gN/m2/s) */
   double hrv_deadstemn_to_prod100n;        /*  dead stem N harvest mortality to 100-year product pool (gN/m2/s) */
   double hrv_livecrootn_to_cwdn;       /*  live coarse root N harvest mortality to coarse woody debris N (gN/m2/s) */
   double hrv_deadcrootn_to_cwdn;       /*  dead coarse root N harvest mortality to coarse woody debris N (gN/m2/s) */
   double hrv_retransn_to_litr1n;       /*  retranslocated N pool harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_leafn_storage_to_litr1n;      /*  leaf N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_frootn_storage_to_litr1n;     /*  fine root N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_livestemn_storage_to_litr1n;      /*  live stem N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_deadstemn_storage_to_litr1n;      /*  dead stem N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_livecrootn_storage_to_litr1n;     /*  live coarse root N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_deadcrootn_storage_to_litr1n;     /*  dead coarse root N storage harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_leafn_xfer_to_litr1n;     /*  leaf N transfer harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_frootn_xfer_to_litr1n;        /*  fine root N transfer harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_livestemn_xfer_to_litr1n;     /*  live stem N transfer harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_deadstemn_xfer_to_litr1n;     /*  dead stem N transfer harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_livecrootn_xfer_to_litr1n;        /*  live coarse root N transfer harvest mortality to litter 1 N (gN/m2/s) */
   double hrv_deadcrootn_xfer_to_litr1n;        /*  dead coarse root N transfer harvest mortality to litter 1 N (gN/m2/s) */
   /* column-level fire fluxes */
   double m_deadstemn_to_cwdn_fire;     /*  dead stem N to coarse woody debris N by fire (gN/m2/s) */
   double m_deadcrootn_to_cwdn_fire;        /*  dead coarse root N to to woody debris N by fire (gN/m2/s) */
   double m_litr1n_to_fire;     /*  litter 1 N fire loss (gN/m2/s) */
   double m_litr2n_to_fire;     /*  litter 2 N fire loss (gN/m2/s) */
   double m_litr3n_to_fire;     /*  litter 3 N fire loss (gN/m2/s) */
   double m_cwdn_to_fire;       /*  coarse woody debris N fire loss (gN/m2/s) */
   /*  litterfall fluxes */
   double livestemn_to_litr1n;      /*  livestem N litterfall to litter 1 N (gN/m2/s) */
   double livestemn_to_litr2n;      /*  livestem N litterfall to litter 2 N (gN/m2/s) */
   double livestemn_to_litr3n;      /*  livestem N litterfall to litter 3 N (gN/m2/s) */
   double leafn_to_litr1n;      /*  leaf N litterfall to litter 1 N (gN/m2/s) */
   double leafn_to_litr2n;      /*  leaf N litterfall to litter 2 N (gN/m2/s) */
   double leafn_to_litr3n;      /*  leaf N litterfall to litter 3 N (gN/m2/s) */
   double frootn_to_litr1n;     /*  fine root N litterfall to litter 1 N (gN/m2/s) */
   double frootn_to_litr2n;     /*  fine root N litterfall to litter 2 N (gN/m2/s) */
   double frootn_to_litr3n;     /*  fine root N litterfall to litter 3 N (gN/m2/s) */
   /* litterfall fluxes for prognostic crop model */
   double grainn_to_litr1n;     /*  grain N litterfall to litter 1 N (gN/m2/s) */
   double grainn_to_litr2n;     /*  grain N litterfall to litter 2 N (gN/m2/s) */
   double grainn_to_litr3n;     /*  grain N litterfall to litter 3 N (gN/m2/s) */
   /* decomposition fluxes */
   double cwdn_to_litr2n;       /*  decomp. of coarse woody debris N to litter 2 N (gN/m2/s) */
   double cwdn_to_litr3n;       /*  decomp. of coarse woody debris N to litter 3 N (gN/m2/s) */
   double litr1n_to_soil1n;     /*  decomp. of litter 1 N to SOM 1 N (gN/m2/s) */
   double sminn_to_soil1n_l1;       /*  mineral N flux for decomp. of litter 1 to SOM 1 (gN/m2/s) */
   double litr2n_to_soil2n;     /*  decomp. of litter 2 N to SOM 2 N (gN/m2/s) */
   double sminn_to_soil2n_l2;       /*  mineral N flux for decomp. of litter 2 to SOM 2 (gN/m2/s) */
   double litr3n_to_soil3n;     /*  decomp. of litter 3 N to SOM 3 N (gN/m2/s) */
   double sminn_to_soil3n_l3;       /*  mineral N flux for decomp. of litter 3 to SOM 3 (gN/m2/s) */
   double soil1n_to_soil2n;     /*  decomp. of SOM 1 N to SOM 2 N (gN/m2/s) */
   double sminn_to_soil2n_s1;       /*  mineral N flux for decomp. of SOM 1 to SOM 2 (gN/m2/s) */
   double soil2n_to_soil3n;     /*  decomp. of SOM 2 N to SOM 3 N (gN/m2/s) */
   double sminn_to_soil3n_s2;       /*  mineral N flux for decomp. of SOM 2 to SOM 3 (gN/m2/s) */
   double soil3n_to_soil4n;     /*  decomp. of SOM 3 N to SOM 4 N (gN/m2/s) */
   double sminn_to_soil4n_s3;       /*  mineral N flux for decomp. of SOM 3 to SOM 4 (gN/m2/s) */
   double soil4n_to_sminn;      /*  N mineralization for decomp. of SOM 4 (gN/m2/s) */
   /* denitrification fluxes */
   double sminn_to_denit_l1s1;      /*  denitrification for decomp. of litter 1 to SOM 1 (gN/m2/s)  */
   double sminn_to_denit_l2s2;      /*  denitrification for decomp. of litter 2 to SOM 2 (gN/m2/s) */
   double sminn_to_denit_l3s3;      /*  denitrification for decomp. of litter 3 to SOM 3 (gN/m2/s) */
   double sminn_to_denit_s1s2;      /*  denitrification for decomp. of SOM 1 to SOM 2 (gN/m2/s) */
   double sminn_to_denit_s2s3;      /*  denitrification for decomp. of SOM 2 to SOM 3 (gN/m2/s) */
   double sminn_to_denit_s3s4;      /*  denitrification for decomp. of SOM 3 to SOM 4 (gN/m2/s) */
   double sminn_to_denit_s4;        /*  denitrification for decomp. of SOM 4 (gN/m2/s) */
   double sminn_to_denit_excess;        /*  denitrification from excess mineral N pool (gN/m2/s) */
   /* leaching fluxes */
   double sminn_leached;        /*  soil mineral N pool loss to leaching (gN/m2/s) */
   /* dynamic landcover fluxes */
   double dwt_seedn_to_leaf;        /*  (gN/m2/s) seed source to PFT-level */
   double dwt_seedn_to_deadstem;        /*  (gN/m2/s) seed source to PFT-level */
   double dwt_conv_nflux;       /*  (gN/m2/s) conversion N flux (immediate loss to atm) */
   double dwt_prod10n_gain;     /*  (gN/m2/s) addition to 10-yr wood product pool */
   double dwt_prod100n_gain;        /*  (gN/m2/s) addition to 100-yr wood product pool */
   double dwt_frootn_to_litr1n;     /*  (gN/m2/s) fine root to litter due to landcover change */
   double dwt_frootn_to_litr2n;     /*  (gN/m2/s) fine root to litter due to landcover change */
   double dwt_frootn_to_litr3n;     /*  (gN/m2/s) fine root to litter due to landcover change */
   double dwt_livecrootn_to_cwdn;       /*  (gN/m2/s) live coarse root to CWD due to landcover change */
   double dwt_deadcrootn_to_cwdn;       /*  (gN/m2/s) dead coarse root to CWD due to landcover change */
   double dwt_nloss;        /*  (gN/m2/s) total nitrogen loss from product pools and conversion */
   /* wood product pool loss fluxes */
   double prod10n_loss;     /*  (gN/m2/s) decomposition loss from 10-yr wood product pool */
   double prod100n_loss;        /*  (gN/m2/s) decomposition loss from 100-yr wood product pool */
   double product_nloss;        /*  (gN/m2/s) total wood product nitrogen loss */
   /* summary (diagnostic) flux variables, not involved in mass balance */
   double potential_immob;      /*  potential N immobilization (gN/m2/s) */
   double actual_immob;     /*  actual N immobilization (gN/m2/s) */
   double sminn_to_plant;       /*  plant uptake of soil mineral N (gN/m2/s) */
   double supplement_to_sminn;      /*  supplemental N supply (gN/m2/s) */
   double gross_nmin;       /*  gross rate of N mineralization (gN/m2/s) */
   double net_nmin;     /*  net rate of N mineralization (gN/m2/s) */
   double denit;        /*  total rate of denitrification (gN/m2/s) */
   double col_ninputs;      /*  column-level N inputs (gN/m2/s) */
   double col_noutputs;     /*  column-level N outputs (gN/m2/s) */
   /* new variables for fire */
   double col_fire_nloss;       /*  total column-level fire N loss (gN/m2/s) */
} nflux_type;

typedef struct estate_type
{
    double t_grnd;      /* ground temperature (Kelvin) */
   double t_grnd_u;     /* Urban ground temperature (Kelvin) */
   double t_grnd_r;     /* Rural ground temperature (Kelvin) */
   double dt_grnd;      /* change in t_grnd, last iteration (Kelvin) */
   double *t_soisno;     /* soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)  */
   double t_soi_10cm;       /* soil temperature in top 10cm of soil (Kelvin) */
   double *t_lake;       /* lake temperature (Kelvin)  (1:nlevlak)           */
   double *tssbef;       /* soil/snow temperature before update (-nlevsno+1:nlevgrnd)  */
   double thv;      /* virtual potential temperature (kelvin) */
   double hc_soi;       /* soil heat content (MJ/m2) */
   double hc_soisno;        /* soil plus snow heat content (MJ/m2) */
   double forc_t;       /* atm temperature, downscaled to column (Kelvin) */
   double forc_th;      /* atm potl temperature, downscaled to column (Kelvin) */
} estate_type;

typedef CNM_struct
{
    carbon_balance_type pcbal;      /* carbon balance structure */
    carbon_balance_type ccbal;      /* carbon balance structure */
    nitrogen_balance_type pnbal;    /* nitrogen balance structure */
    nitrogen_balance_type cnbal;    /* nitrogen balance structure */
    pstate_type ps;            /* physical state variables */
    epc_type pftcon;
    epv_type epv;              /* pft ecophysiological variables */
    cflux_type cf;             /* pft carbon flux */
    nflux_type nf;             /* pft nitrogen flux */
    cstate_type cs;            /* pft carbon state */
    nstate_type ns;             /* pft nitrogen state */
    estate_type es;             /* energy state */
} *CNM_struct;    
