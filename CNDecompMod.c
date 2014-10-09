void CNDecompAlloc ()
{
/* DESCRIPTION:
 * Module holding routines used in litter and soil decomposition model
 * for coupled carbon-nitrogen code.
 */
    int c, j;       /* indices */
    double dt;      /* decomp timestep (seconds) */
    double dtd;     /* decomp timestep (days) */
    double *fr;     /* column-level rooting fraction by soil depth */
    double frw;     /* rooting fraction weight */
    double t_scalar;    /* soil temperature scalar for decomp */
    double minpsi, maxpsi;  /* limits for soil water scalar for decomp */
    double psi;             /* temporary soilpsi for water scalar */
    double w_scalar;        /* soil water scalar for decomp */
    double rate_scalar;     /* combined rate scalar for decomp */
    double cn_l1;           /* C:N for litter 1 */
    double cn_l2(lbc:ubc);  /* C:N for litter 2 */
    double cn_l3(lbc:ubc)   /* C:N for litter 3 */
    double cn_s1;           /* C:N for SOM 1 */
    double cn_s2;           /* C:N for SOM 2 */
    double cn_s3;           /* C:N for SOM 3 */
    double cn_s4;           /* C:N for SOM 4 */
    double rf_l1s1;         /* respiration fraction litter 1 -> SOM 1 */
    double rf_l2s2;         /* respiration fraction litter 2 -> SOM 2 */
    double rf_l3s3;         /* respiration fraction litter 3 -> SOM 3 */
    double rf_s1s2;         /* respiration fraction SOM 1 -> SOM 2 */
    double rf_s2s3;         /* respiration fraction SOM 2 -> SOM 3 */
    double rf_s3s4;         /* respiration fraction SOM 3 -> SOM 4 */
    double k_l1;            /* decomposition rate constant litter 1 */
    double k_l2;            /* decomposition rate constant litter 2 */
    double k_l3;            /* decomposition rate constant litter 3 */
    double k_s1;            /* decomposition rate constant SOM 1 */
    double k_s2;            /* decomposition rate constant SOM 2 */
    double k_s3;            /* decomposition rate constant SOM 3 */
    double k_s4;            /* decomposition rate constant SOM 3 */
    double k_frag;          /* fragmentation rate constant CWD */
    double ck_l1;           /* corrected decomposition rate constant litter 1 */
    double ck_l2;           /* corrected decomposition rate constant litter 2 */
    double ck_l3;           /* corrected decomposition rate constant litter 3 */
    double ck_s1;           /* corrected decomposition rate constant SOM 1 */
    double ck_s2;           /* corrected decomposition rate constant SOM 2 */
    double ck_s3;           /* corrected decomposition rate constant SOM 3 */
    double ck_s4;           /* corrected decomposition rate constant SOM 3 */
    double ck_frag;         /* corrected fragmentation rate constant CWD */
    double cwd_fcel;            /* cellulose fraction of coarse woody debris */
    double cwd_flig;            /* lignin fraction of coarse woody debris */
    double cwdc_loss;           /* fragmentation rate for CWD carbon (gC/m2/s) */
    double cwdn_loss;           /* fragmentation rate for CWD nitrogen (gN/m2/s) */
    double plitr1c_loss(lbc:ubc);           /* potential C loss from litter 1 */
    double plitr2c_loss(lbc:ubc);           /* potential C loss from litter 2 */
    double plitr3c_loss(lbc:ubc);           /* potential C loss from litter 3 */
    double psoil1c_loss(lbc:ubc);           /* potential C loss from SOM 1 */
    double psoil2c_loss(lbc:ubc);           /* potential C loss from SOM 2 */
    double psoil3c_loss(lbc:ubc);           /* potential C loss from SOM 3 */
    double psoil4c_loss(lbc:ubc);           /* potential C loss from SOM 4 */
    double pmnf_l1s1(lbc:ubc);          /* potential mineral N flux, litter 1 -> SOM 1 */
    double pmnf_l2s2(lbc:ubc);          /* potential mineral N flux, litter 2 -> SOM 2 */
    double pmnf_l3s3(lbc:ubc);          /* potential mineral N flux, litter 3 -> SOM 3 */
    double pmnf_s1s2(lbc:ubc);          /* potential mineral N flux, SOM 1 -> SOM 2 */
    double pmnf_s2s3(lbc:ubc);          /* potential mineral N flux, SOM 2 -> SOM 3 */
    double pmnf_s3s4(lbc:ubc);          /* potential mineral N flux, SOM 3 -> SOM 4 */
    double pmnf_s4(lbc:ubc);            /* potential mineral N flux, SOM 4 */
    double immob(lbc:ubc);          /* potential N immobilization */
    double ratio;           /* temporary variable */
    double dnp;         /* denitrification proportion */
    int nlevdecomp;     /* bottom layer to consider for decomp controls */
    double spinup_scalar;           /* multiplier for AD_SPINUP algorithm */

/*
   t_soisno              => ces%t_soisno
   psisat                => cps%psisat
   soilpsi               => cps%soilpsi
   dz                    => cps%dz
   cwdc                  => ccs%cwdc
   litr1c                => ccs%litr1c
   litr2c                => ccs%litr2c
   litr3c                => ccs%litr3c
   soil1c                => ccs%soil1c
   soil2c                => ccs%soil2c
   soil3c                => ccs%soil3c
   soil4c                => ccs%soil4c
   cwdn                  => cns%cwdn
   litr1n                => cns%litr1n
   litr2n                => cns%litr2n
   litr3n                => cns%litr3n
   fpi                   => cps%fpi
   cwdc_to_litr2c        => ccf%cwdc_to_litr2c
   cwdc_to_litr3c        => ccf%cwdc_to_litr3c
   litr1_hr              => ccf%litr1_hr
   litr1c_to_soil1c      => ccf%litr1c_to_soil1c
   litr2_hr              => ccf%litr2_hr
   litr2c_to_soil2c      => ccf%litr2c_to_soil2c
   litr3_hr              => ccf%litr3_hr
   litr3c_to_soil3c      => ccf%litr3c_to_soil3c
   soil1_hr              => ccf%soil1_hr
   soil1c_to_soil2c      => ccf%soil1c_to_soil2c
   soil2_hr              => ccf%soil2_hr
   soil2c_to_soil3c      => ccf%soil2c_to_soil3c
   soil3_hr              => ccf%soil3_hr
   soil3c_to_soil4c      => ccf%soil3c_to_soil4c
   soil4_hr              => ccf%soil4_hr
   cwdn_to_litr2n        => cnf%cwdn_to_litr2n
   cwdn_to_litr3n        => cnf%cwdn_to_litr3n
   potential_immob       => cnf%potential_immob
   litr1n_to_soil1n      => cnf%litr1n_to_soil1n
   sminn_to_soil1n_l1    => cnf%sminn_to_soil1n_l1
   litr2n_to_soil2n      => cnf%litr2n_to_soil2n
   sminn_to_soil2n_l2    => cnf%sminn_to_soil2n_l2
   litr3n_to_soil3n      => cnf%litr3n_to_soil3n
   sminn_to_soil3n_l3    => cnf%sminn_to_soil3n_l3
   soil1n_to_soil2n      => cnf%soil1n_to_soil2n
   sminn_to_soil2n_s1    => cnf%sminn_to_soil2n_s1
   soil2n_to_soil3n      => cnf%soil2n_to_soil3n
   sminn_to_soil3n_s2    => cnf%sminn_to_soil3n_s2
   soil3n_to_soil4n      => cnf%soil3n_to_soil4n
   sminn_to_soil4n_s3    => cnf%sminn_to_soil4n_s3
   soil4n_to_sminn       => cnf%soil4n_to_sminn
   sminn_to_denit_l1s1   => cnf%sminn_to_denit_l1s1
   sminn_to_denit_l2s2   => cnf%sminn_to_denit_l2s2
   sminn_to_denit_l3s3   => cnf%sminn_to_denit_l3s3
   sminn_to_denit_s1s2   => cnf%sminn_to_denit_s1s2
   sminn_to_denit_s2s3   => cnf%sminn_to_denit_s2s3
   sminn_to_denit_s3s4   => cnf%sminn_to_denit_s3s4
   sminn_to_denit_s4     => cnf%sminn_to_denit_s4
   sminn_to_denit_excess => cnf%sminn_to_denit_excess
   gross_nmin            => cnf%gross_nmin
   net_nmin              => cnf%net_nmin
   rootfr                => pps%rootfr
   clandunit             => col%landunit
   itypelun              => lun%itype
*/

    /* set time steps */
    dt = (double)get_step_size();
    dtd = dt / secspday;

    /* set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0) */
    cn_s1 = 12.0;
    cn_s2 = 12.0;
    cn_s3 = 10.0;
    cn_s4 = 10.0;

    /* set respiration fractions for fluxes between compartments
     * (from Biome-BGC v4.2.0) */
    rf_l1s1 = 0.39;
    rf_l2s2 = 0.55;
    rf_l3s3 = 0.29;
    rf_s1s2 = 0.28;
    rf_s2s3 = 0.46;
    rf_s3s4 = 0.55;

    /* set the cellulose and lignin fractions for coarse woody debris */
    cwd_fcel = 0.76;
    cwd_flig = 0.24;

    /* set initial base rates for decomposition mass loss (1/day)
     * (from Biome-BGC v4.2.0, using three SOM pools)
     * Value inside log function is the discrete-time values for a
     * daily time step model, and the result of the log function is
     * the corresponding continuous-time decay rate (1/day), following
     * Olson, 1963. */
    k_l1 = -log(1.0 - 0.7);
    k_l2 = -log(1.0 - 0.07);
    k_l3 = -log(1.0 - 0.014);
    k_s1 = -log(1.0 - 0.07);
    k_s2 = -log(1.0-0.014);
    k_s3 = -log(1.0-0.0014);
    k_s4 = -log(1.0-0.0001);
    k_frag = -log(1.0-0.001);

    /* calculate the new discrete-time decay rate for model timestep */
   k_l1 = 1.0-exp(-k_l1*dtd);
   k_l2 = 1.0-exp(-k_l2*dtd);
   k_l3 = 1.0-exp(-k_l3*dtd);
   k_s1 = 1.0-exp(-k_s1*dtd);
   k_s2 = 1.0-exp(-k_s2*dtd);
   k_s3 = 1.0-exp(-k_s3*dtd);
   k_s4 = 1.0-exp(-k_s4*dtd);
   k_frag = 1.0-exp(-k_frag*dtd);
   
    /* The following code implements the acceleration part of the AD spinup
     * algorithm, by multiplying all of the SOM decomposition base rates by 10.0. */

    if (use_ad_spinup)
    {
        spinup_scalar = 20.;
        k_s1 = k_s1 * spinup_scalar;
        k_s2 = k_s2 * spinup_scalar;
        k_s3 = k_s3 * spinup_scalar;
        k_s4 = k_s4 * spinup_scalar;
    }

    /* calculate function to weight the temperature and water potential scalars
     * for decomposition control. */

    /* the following normalizes values in fr so that they
     * sum to 1.0 across top nlevdecomp levels on a column */
    frw = 0.;
    nlevdecomp = 5;
    fr = (double *) malloc (nlevdecomp * sizeof(double));
    for (j = 0; j < nlevdecomp; j++)
        frw = frw + ps->dz[j];
    for (j = 0; j < nlevdecomp; j++)
    {
        if (frw != 0.)
            fr[j] = ps->dz[j] / frw;
         else
            fr[j] = 0.;
    }

    /* calculate rate constant scalar for soil temperature
     * assuming that the base rate constants are assigned for non-moisture
     * limiting conditions at 25 C.
     * Peter Thornton: 3/13/09
     * Replaced the Lloyd and Taylor function with a Q10 formula
     * with Q10 = 1.5 as part of the modifications made to improve the
     * seasonal cycle of atmospheric CO2 concentration in global simulations.
     * This does not impact the base rates at 25 C, which are calibrated from
     * microcosm studies. */
    t_scalar = 0.;
    for (j = 0; j < nlevdecomp; j++)
        t_scalar = t_scalar + pow(1.5, (es->t_soisno[j] - (SHR_CONST_TKFRZ + 25.)) / 10.) * fr[j]

    /* calculate the rate constant scalar for soil water content.
     * Uses the log relationship with water potential given in
     * Andren, O., and K. Paustian, 1987. Barley straw decomposition in the
     * field: a comparison of models. Ecology, 68(5):1190-1200.
     * and supported by data in
     * Orchard, V.A., and F.J. Cook, 1983. Relationship between soil
     * respiration and soil moisture. Soil Biol. Biochem., 15(4):447-453. */

    minpsi = -10.0;
    w_scalar = 0.;
    for (j = 0; j < nlevdecomp; j++)
    {
        maxpsi = ps->psisat[j];
        psi = ps->soilpsi[j] < maxpsi ? ps->soilpsi[j] : maxpsi;
        /* decomp only if soilpsi is higher than minpsi */
        if (psi > minpsi)
            w_scalar = w_scalar + (log(minpsi / psi) / log(minpsi / maxpsi)) * fr[j];
    }

    /* set initial values for potential C and N fluxes */
    plitr1c_loss = 0.;
    plitr2c_loss = 0.;
    plitr3c_loss = 0.;
    psoil1c_loss = 0.;
    psoil2c_loss = 0.;
    psoil3c_loss = 0.;
    psoil4c_loss = 0.;
    pmnf_l1s1 = 0.;
    pmnf_l2s2 = 0.;
    pmnf_l3s3 = 0.;
    pmnf_s1s2 = 0.;
    pmnf_s2s3 = 0.;
    pmnf_s3s4 = 0.;
    pmnf_s4 = 0.;

    /* column loop to calculate potential decomp rates and total
     * immobilization demand. */

    /* calculate litter compartment C:N ratios */
    if (ns->litr1n > 0.)
        cn_l1 = cs->litr1c / ns->litr1n;
    if (ns->litr2n > 0.)
        cn_l2 = cs->litr2c / ns->litr2n;
    if (ns->litr3n > 0.)
        cn_l3 = cs->litr3c / ns->litr3n;

    /* calculate the final rate scalar as the product of temperature and water
     * rate scalars, and correct the base decomp rates */

    rate_scalar = t_scalar * w_scalar;
    ck_l1 = k_l1 * rate_scalar;
    ck_l2 = k_l2 * rate_scalar;
    ck_l3 = k_l3 * rate_scalar;
    ck_s1 = k_s1 * rate_scalar;
    ck_s2 = k_s2 * rate_scalar;
    ck_s3 = k_s3 * rate_scalar;
    ck_s4 = k_s4 * rate_scalar;
    ck_frag = k_frag * rate_scalar;

    /* calculate the non-nitrogen-limited fluxes
     * these fluxes include the  "/ dt" term to put them on a
     * per second basis, since the rate constants have been
     * calculated on a per timestep basis. */

    /* CWD fragmentation -> litter pools */
    cwdc_loss = cs->cwdc * ck_frag / dt;
    cf->cwdc_to_litr2c = cwdc_loss * cwd_fcel;
    cf->cwdc_to_litr3c = cwdc_loss * cwd_flig;
    cwdn_loss = ns->cwdn * ck_frag / dt;
    nf->cwdn_to_litr2n = cwdn_loss * cwd_fcel;
    nf->cwdn_to_litr3n = cwdn_loss * cwd_flig;

    /* litter 1 -> SOM 1 */
    if (cs->litr1c > 0.)
    {
        plitr1c_loss = cs->litr1c * ck_l1 / dt;
        ratio = 0.;
        if (ns->litr1n > 0.)
            ratio = cn_s1 / cn_l1;
        pmnf_l1s1 = (plitr1c_loss * (1.0 - rf_l1s1 - ratio)) / cn_s1;
    }

    /* litter 2 -> SOM 2 */
    if (cs->litr2c > 0.)
    {
        plitr2c_loss = cs->litr2c * ck_l2 / dt;
        ratio = 0.
        if (ns->litr2n > 0.)
            ratio = cn_s2 / cn_l2;
        pmnf_l2s2 = (plitr2c_loss * (1.0 - rf_l2s2 - ratio)) / cn_s2;
    }

    /* litter 3 -> SOM 3 */
    if (cs->litr3c > 0.)
    {
        plitr3c_loss = cs->litr3c * ck_l3 / dt;
        ratio = 0.;
        if (ns->litr3n > 0.)
            ratio = cn_s3 / cn_l3;
        pmnf_l3s3 = (plitr3c_loss * (1.0 - rf_l3s3 - ratio)) / cn_s3;
    }

    /* SOM 1 -> SOM 2 */
    if (cs->soil1c > 0.)
    {
        psoil1c_loss = cs->soil1c * ck_s1 / dt;
        pmnf_s1s2 = (psoil1c_loss * (1.0 - rf_s1s2 - (cn_s2 / cn_s1))) / cn_s2;
    }

    /* SOM 2 -> SOM 3 */
    if (cs->soil2c > 0.)
    {
        psoil2c_loss = cs->soil2c * ck_s2 / dt;
        pmnf_s2s3 = (psoil2c_loss * (1.0 - rf_s2s3 - (cn_s3 / cn_s2))) / cn_s3;
    }

    /* SOM 3 -> SOM 4 */
    if (cs->soil3c > 0.)
    {
        psoil3c_loss = cs->soil3c * ck_s3 / dt;
        pmnf_s3s4 = (psoil3c_loss * (1.0 - rf_s3s4 - (cn_s4 / cn_s3))) / cn_s4;
    }

    /* Loss from SOM 4 is entirely respiration (no downstream pool) */
    if (cs->soil4c > 0.)
    {
        psoil4c_loss = cs->soil4c * ck_s4 / dt;
        pmnf_s4 = -psoil4c_loss / cn_s4;
    }

    /* Sum up all the potential immobilization fluxes (positive pmnf flux)
     * and all the mineralization fluxes (negative pmnf flux) */

    immob = 0.;
    /* litter 1 -> SOM 1 */
    if (pmnf_l1s1 > 0.)
        immob = immob + pmnf_l1s1;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_l1s1;

    /* litter 2 -> SOM 2 */
    if (pmnf_l2s2 > 0.)
        immob = immob + pmnf_l2s2;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_l2s2;

    /* litter 3 -> SOM 3 */
    if (pmnf_l3s3 > 0.)
        immob = immob + pmnf_l3s3;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_l3s3;

    /* SOM 1 -> SOM 2 */
    if (pmnf_s1s2 > 0.)
        immob = immob + pmnf_s1s2;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_s1s2;

    /* SOM 2 -> SOM 3 */
    if (pmnf_s2s3 > 0.)
        immob = immob + pmnf_s2s3;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_s2s3;

    /* SOM 3 -> SOM 4 */
    if (pmnf_s3s4 > 0.)
        immob = immob + pmnf_s3s4;
    else
        nf->gross_nmin = nf->gross_nmin - pmnf_s3s4;
      end if

    /* SOM 4 */
    nf->gross_nmin = nf->gross_nmin - pmnf_s4;

    nf->potential_immob = immob;


    /* now that potential N immobilization is known, call allocation
     * to resolve the competition between plants and soil heterotrophs
     * for available soil mineral N resource. */
  
    CNAllocation(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp, filter_soilp, num_pcropp);

   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N

   dnp = 0.01

   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! upon return from CNAllocation, the fraction of potential immobilization
      ! has been set (cps%fpi). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.

      ! litter 1 fluxes (labile pool)
      if (cs->litr1c(c) > 0.) then
         if (pmnf_l1s1(c) > 0.) then
            plitr1c_loss(c) = plitr1c_loss(c) * ps->fpi(c)
            pmnf_l1s1(c) = pmnf_l1s1(c) * ps->fpi(c)
            sminn_to_denit_l1s1(c) = 0.
         else
            sminn_to_denit_l1s1(c) = -dnp * pmnf_l1s1(c)
         end if
         cf->litr1_hr(c) = rf_l1s1 * plitr1c_loss(c)
         cf->litr1c_to_soil1c(c) = (1. - rf_l1s1) * plitr1c_loss(c)
         if (ns->litr1n(c) > 0.) litr1n_to_soil1n(c) = plitr1c_loss(c) / cn_l1(c)
         sminn_to_soil1n_l1(c) = pmnf_l1s1(c)
         net_nmin(c) = net_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 fluxes (cellulose pool)
      if (cs->litr2c(c) > 0.) then
         if (pmnf_l2s2(c) > 0.) then
            plitr2c_loss(c) = plitr2c_loss(c) * ps->fpi(c)
            pmnf_l2s2(c) = pmnf_l2s2(c) * ps->fpi(c)
            sminn_to_denit_l2s2(c) = 0.
         else
            sminn_to_denit_l2s2(c) = -dnp * pmnf_l2s2(c)
         end if
         cf->litr2_hr(c) = rf_l2s2 * plitr2c_loss(c)
         cf->litr2c_to_soil2c(c) = (1. - rf_l2s2) * plitr2c_loss(c)
         if (ns->litr2n(c) > 0.) litr2n_to_soil2n(c) = plitr2c_loss(c) / cn_l2(c)
         sminn_to_soil2n_l2(c) = pmnf_l2s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 fluxes (lignin pool)
      if (cs->litr3c(c) > 0.) then
         if (pmnf_l3s3(c) > 0.) then
            plitr3c_loss(c) = plitr3c_loss(c) * ps->fpi(c)
            pmnf_l3s3(c) = pmnf_l3s3(c) * ps->fpi(c)
            sminn_to_denit_l3s3(c) = 0.
         else
            sminn_to_denit_l3s3(c) = -dnp * pmnf_l3s3(c)
         end if
         cf->litr3_hr(c) = rf_l3s3 * plitr3c_loss(c)
         cf->litr3c_to_soil3c(c) = (1. - rf_l3s3) * plitr3c_loss(c)
         if (ns->litr3n(c) > 0.) litr3n_to_soil3n(c) = plitr3c_loss(c) / cn_l3(c)
         sminn_to_soil3n_l3(c) = pmnf_l3s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 fluxes (fast rate soil organic matter pool)
      if (cs->soil1c(c) > 0.) then
         if (pmnf_s1s2(c) > 0.) then
            psoil1c_loss(c) = psoil1c_loss(c) * ps->fpi(c)
            pmnf_s1s2(c) = pmnf_s1s2(c) * ps->fpi(c)
            sminn_to_denit_s1s2(c) = 0.
         else
            sminn_to_denit_s1s2(c) = -dnp * pmnf_s1s2(c)
         end if
         cf->soil1_hr(c) = rf_s1s2 * psoil1c_loss(c)
         cf->soil1c_to_soil2c(c) = (1. - rf_s1s2) * psoil1c_loss(c)
         soil1n_to_soil2n(c) = psoil1c_loss(c) / cn_s1
         sminn_to_soil2n_s1(c) = pmnf_s1s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 fluxes (medium rate soil organic matter pool)
      if (cs->soil2c(c) > 0.) then
         if (pmnf_s2s3(c) > 0.) then
            psoil2c_loss(c) = psoil2c_loss(c) * ps->fpi(c)
            pmnf_s2s3(c) = pmnf_s2s3(c) * ps->fpi(c)
            sminn_to_denit_s2s3(c) = 0.
         else
            sminn_to_denit_s2s3(c) = -dnp * pmnf_s2s3(c)
         end if
         cf->soil2_hr(c) = rf_s2s3 * psoil2c_loss(c)
         cf->soil2c_to_soil3c(c) = (1. - rf_s2s3) * psoil2c_loss(c)
         soil2n_to_soil3n(c) = psoil2c_loss(c) / cn_s2
         sminn_to_soil3n_s2(c) = pmnf_s2s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 fluxes (slow rate soil organic matter pool)
      if (cs->soil3c(c) > 0.) then
         if (pmnf_s3s4(c) > 0.) then
            psoil3c_loss(c) = psoil3c_loss(c) * ps->fpi(c)
            pmnf_s3s4(c) = pmnf_s3s4(c) * ps->fpi(c)
            sminn_to_denit_s3s4(c) = 0.
         else
            sminn_to_denit_s3s4(c) = -dnp * pmnf_s3s4(c)
         end if
         cf->soil3_hr(c) = rf_s3s4 * psoil3c_loss(c)
         cf->soil3c_to_soil4c(c) = (1. - rf_s3s4) * psoil3c_loss(c)
         soil3n_to_soil4n(c) = psoil3c_loss(c) / cn_s3
         sminn_to_soil4n_s3(c) = pmnf_s3s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4 fluxes (slowest rate soil organic matter pool)
      if (cs->soil4c(c) > 0.) then
         cf->soil4_hr(c) = psoil4c_loss(c)
         soil4n_to_sminn(c) = psoil4c_loss(c) / cn_s4
         sminn_to_denit_s4(c) = -dnp * pmnf_s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s4(c)
      end if

   end do

   deallocate(fr)

end subroutine CNDecompAlloc
!-----------------------------------------------------------------------

end module CNDecompMod
