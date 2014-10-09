void CNNDeposition(nflux_type *nf, Model_Data DS, double t)
{
/*
 * DESCRIPTION:
 * On the radiation time step, update the nitrogen deposition rate
 * from atmospheric forcing. For now it is assumed that all the atmospheric
 * N deposition goes to the soil mineral N pool.
 * This could be updated later to divide the inputs between mineral N absorbed
 * directly into the canopy and mineral N entering the soil pool.
 */

    nf->ndep_to_sminn = Interpolation (&DS->Forcing[14][0], t); /* nitrogen deposition rate (gN/m2/s) */
}

void CNNFixation(pstate_type *ps, nflux_type *nf)
{
/*
 * DESCRIPTION:
 * On the radiation time step, update the nitrogen fixation rate
 * as a function of annual total NPP. This rate gets updated once per year.
 * All N fixation goes to the soil mineral N pool.
 */
   integer  :: c,fc                  ! indices
    double t;       /* temporary */
   real(r8) :: dayspyr               ! days per year

   cannsum_npp   => cps%cannsum_npp

   ! Assign local pointers to derived type arrays (out)
   nfix_to_sminn => cnf%nfix_to_sminn

   dayspyr = get_days_per_year()

    /* 
     * the value 0.001666 is set to give 100 TgN/yr when global
     * NPP = 60 PgC/yr.  (Cleveland et al., 1999)
     * Convert from gN/m2/yr -> gN/m2/s
     */

//  t = cannsum_npp(c) * 0.001666_r8 / (secspday * dayspyr)
    t = (1.8 * (1. - exp(-0.003 * ps->annsum_npp))) / (secspday * dayspyr);
    nf->nfix_to_sminn = t > 0. ? t : 0.;
    /* PET 2/14/05: commenting out the dependence on NPP, and
     * forcing Nfix to global constant = 0.4 gN/m2/yr
     * nfix_to_sminn(c) = 0.4 / (secspday*dayspyr) */
}
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNNLeaching
!
! !INTERFACE:
subroutine CNNLeaching(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, update the nitrogen leaching rate
! as a function of soluble mineral N and total soil water outflow.
!
! !USES:
   use clmtype
   use clm_varpar      , only : nlevsoi
   use clm_time_manager    , only : get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 6/9/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
   real(r8), pointer :: sminn(:)         ! (gN/m2) soil mineral N
!
! local pointers to implicit out scalars
!
   real(r8), pointer :: sminn_leached(:) ! rate of mineral N leaching (gN/m2/s)
!
! !OTHER LOCAL VARIABLES:
   integer  :: j,c,fc             ! indices
   real(r8) :: dt                 ! radiation time step (seconds)
   real(r8) :: tot_water(lbc:ubc) ! total column liquid water (kg water/m2)
   real(r8) :: sf                 ! soluble fraction of mineral N (unitless)
   real(r8) :: disn_conc          ! dissolved mineral N concentration
                                  ! (gN/kg water)

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   h2osoi_liq    => cws%h2osoi_liq
   qflx_drain    => cwf%qflx_drain
   sminn         => cns%sminn

   ! Assign local pointers to derived type arrays (out)
   sminn_leached => cnf%sminn_leached

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! Assume that 10% of the soil mineral N is in a soluble form
   sf = 0.1_r8

   ! calculate the total soil water
   tot_water(lbc:ubc) = 0._r8
   do j = 1,nlevsoi
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
      end do
   end do

   ! Loop through columns
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! calculate the dissolved mineral N concentration (gN/kg water)
      ! assumes that 10% of mineral nitrogen is soluble
      disn_conc = 0._r8
      if (tot_water(c) > 0._r8) then
         disn_conc = (sf * sminn(c))/tot_water(c)
      end if

      ! calculate the N leaching flux as a function of the dissolved
      ! concentration and the sub-surface drainage flux
      sminn_leached(c) = disn_conc * qflx_drain(c)

      ! limit the flux based on current sminn state
      ! only let at most the assumed soluble fraction
      ! of sminn be leached on any given timestep
      sminn_leached(c) = min(sminn_leached(c), (sf * sminn(c))/dt)
      
      ! limit the flux to a positive value
      sminn_leached(c) = max(sminn_leached(c), 0._r8)

   end do

end subroutine CNNLeaching

end module CNNDynamicsMod
