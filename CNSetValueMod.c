void CNZeroFluxes(pft_cflux_type *pcf, pft_nflux_type *pnf)
{
    /*
     * zero the pft-level C and N fluxes
     */
    CNSetPcf (0., pcf);
    CNSetPnf (0., pnf);
}

void CNZeroFluxes_dwt(column_cflux_type *ccf)
{
    int c, p;           /* indices */
//    type(column_type),   pointer :: cptr         ! pointer to column derived subtype

//    cptr => col
    /*
     * set column-level conversion and product pool fluxes
     * to 0 at the beginning of every timestep
     */

    /* C fluxes */
    ccf->dwt_seedc_to_leaf = 0.;
    ccf->dwt_seedc_to_deadstem = 0.;
    ccf->dwt_conv_cflux = 0.;
    ccf->dwt_prod10c_gain = 0.;
    ccf->dwt_prod100c_gain = 0.;
    ccf->dwt_frootc_to_litr1c = 0.;
    ccf->dwt_frootc_to_litr2c = 0.;
    ccf->dwt_frootc_to_litr3c = 0.;
    ccf->dwt_livecrootc_to_cwdc = 0.;
    ccf->dwt_deadcrootc_to_cwdc = 0.;
    /* N fluxes */
    cnf->dwt_seedn_to_leaf = 0.;
    cnf->dwt_seedn_to_deadstem = 0.;
    cnf->dwt_conv_nflux = 0.;
    cnf->dwt_prod10n_gain = 0.;
    cnf->dwt_prod100n_gain = 0.;
    cnf->dwt_frootn_to_litr1n = 0.;
    cnf->dwt_frootn_to_litr2n = 0.;
    cnf->dwt_frootn_to_litr3n = 0.;
    cnf->dwt_livecrootn_to_cwdn = 0.;
    cnf->dwt_deadcrootn_to_cwdn = 0.;

    pcs->dispvegc   = 0.;
    pcs->storvegc   = 0.;
    pcs->totpftc    = 0.;
    pns->dispvegn   = 0.;
    pns->storvegn   = 0.;
    pns->totvegn    = 0.;
    pns->totpftn    = 0.;
}    

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetPps
!
! !INTERFACE:
subroutine CNSetPps(num, filter, val, pps)
!
! !DESCRIPTION:
! Set pft physical state variables
! !USES:
    use clm_varpar  , only : numrad
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (pft_pstate_type), intent(inout) :: pps
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i,j     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pps%slasun(i) = val
      pps%slasha(i) = val
      pps%lncsun(i) = val
      pps%lncsha(i) = val
      pps%vcmxsun(i) = val
      pps%vcmxsha(i) = val
      pps%gdir(i) = val
   end do

   do j = 1,numrad
      do fi = 1,num
         i = filter(fi)
         pps%omega(i,j) = val
         pps%eff_kid(i,j) = val
         pps%eff_kii(i,j) = val
         pps%sun_faid(i,j) = val
         pps%sun_faii(i,j) = val
         pps%sha_faid(i,j) = val
         pps%sha_faii(i,j) = val
      end do
   end do

end subroutine CNSetPps
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetPepv
!
! !INTERFACE:
subroutine CNSetPepv (num, filter, val, pepv)
!
! !DESCRIPTION:
! Set pft ecophysiological variables
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (pft_epv_type), intent(inout) :: pepv
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pepv%dormant_flag(i) = val
      pepv%days_active(i) = val
      pepv%onset_flag(i) = val
      pepv%onset_counter(i) = val
      pepv%onset_gddflag(i) = val
      pepv%onset_fdd(i) = val
      pepv%onset_gdd(i) = val
      pepv%onset_swi(i) = val
      pepv%offset_flag(i) = val
      pepv%offset_counter(i) = val
      pepv%offset_fdd(i) = val
      pepv%offset_swi(i) = val
      pepv%lgsf(i) = val
      pepv%bglfr(i) = val
      pepv%bgtr(i) = val
      pepv%dayl(i) = val
      pepv%prev_dayl(i) = val
      pepv%annavg_t2m(i) = val
      pepv%tempavg_t2m(i) = val
      pepv%gpp(i) = val
      pepv%availc(i) = val
      pepv%xsmrpool_recover(i) = val
      if (use_c13) then
         pepv%xsmrpool_c13ratio(i) = val
      end if
      pepv%alloc_pnow(i) = val
      pepv%c_allometry(i) = val
      pepv%n_allometry(i) = val
      pepv%plant_ndemand(i) = val
      pepv%tempsum_potential_gpp(i) = val
      pepv%annsum_potential_gpp(i) = val
      pepv%tempmax_retransn(i) = val
      pepv%annmax_retransn(i) = val
      pepv%avail_retransn(i) = val
      pepv%plant_nalloc(i) = val
      pepv%plant_calloc(i) = val
      pepv%excess_cflux(i) = val
      pepv%downreg(i) = val
      pepv%prev_leafc_to_litter(i) = val
      pepv%prev_frootc_to_litter(i) = val
      pepv%tempsum_npp(i) = val
      pepv%annsum_npp(i) = val
      if (use_cndv) then
         pepv%tempsum_litfall(i) = val
         pepv%annsum_litfall(i) = val
      end if
   end do

end subroutine CNSetPepv
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetPcs
!
! !INTERFACE:
subroutine CNSetPcs (num, filter, val, pcs)
!
! !DESCRIPTION:
! Set pft carbon state variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (pft_cstate_type), intent(inout) :: pcs
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pcs%leafc(i) = val
      pcs%leafc_storage(i) = val
      pcs%leafc_xfer(i) = val
      pcs%frootc(i) = val
      pcs%frootc_storage(i) = val
      pcs%frootc_xfer(i) = val
      pcs%livestemc(i) = val
      pcs%livestemc_storage(i) = val
      pcs%livestemc_xfer(i) = val
      pcs%deadstemc(i) = val
      pcs%deadstemc_storage(i) = val
      pcs%deadstemc_xfer(i) = val
      pcs%livecrootc(i) = val
      pcs%livecrootc_storage(i) = val
      pcs%livecrootc_xfer(i) = val
      pcs%deadcrootc(i) = val
      pcs%deadcrootc_storage(i) = val
      pcs%deadcrootc_xfer(i) = val
      pcs%gresp_storage(i) = val
      pcs%gresp_xfer(i) = val
      pcs%cpool(i) = val
      pcs%xsmrpool(i) = val
      pcs%pft_ctrunc(i) = val
      pcs%dispvegc(i) = val
      pcs%storvegc(i) = val
      pcs%totvegc(i) = val
      pcs%totpftc(i) = val
      pcs%woodc(i) = val

      if ( crop_prog )then
         pcs%grainc(i)         = val
         pcs%grainc_storage(i) = val
         pcs%grainc_xfer(i)    = val
      end if
   end do

end subroutine CNSetPcs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetPns
!
! !INTERFACE:
subroutine CNSetPns(num, filter, val, pns)
!
! !DESCRIPTION:
! Set pft nitrogen state variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (pft_nstate_type), intent(inout) :: pns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pns%leafn(i) = val
      pns%leafn_storage(i) = val
      pns%leafn_xfer(i) = val
      pns%frootn(i) = val
      pns%frootn_storage(i) = val
      pns%frootn_xfer(i) = val
      pns%livestemn(i) = val
      pns%livestemn_storage(i) = val
      pns%livestemn_xfer(i) = val
      pns%deadstemn(i) = val
      pns%deadstemn_storage(i) = val
      pns%deadstemn_xfer(i) = val
      pns%livecrootn(i) = val
      pns%livecrootn_storage(i) = val
      pns%livecrootn_xfer(i) = val
      pns%deadcrootn(i) = val
      pns%deadcrootn_storage(i) = val
      pns%deadcrootn_xfer(i) = val
      pns%retransn(i) = val
      pns%npool(i) = val
      pns%pft_ntrunc(i) = val
      pns%dispvegn(i) = val
      pns%storvegn(i) = val
      pns%totvegn(i) = val
      pns%totpftn(i) = val
      if ( crop_prog )then
         pns%grainn(i)         = val
         pns%grainn_storage(i) = val
         pns%grainn_xfer(i)    = val
      end if
   end do

end subroutine CNSetPns
void CNSetPcf(double val, pft_cflux_type *pcf)
{
    pcf->m_leafc_to_litter = val;
    pcf->m_frootc_to_litter = val;
    pcf->m_leafc_storage_to_litter = val;
    pcf->m_frootc_storage_to_litter = val;
    pcf->m_livestemc_storage_to_litter = val;
    pcf->m_deadstemc_storage_to_litter = val;
    pcf->m_livecrootc_storage_to_litter = val;
    pcf->m_deadcrootc_storage_to_litter = val;
    pcf->m_leafc_xfer_to_litter = val;
    pcf->m_frootc_xfer_to_litter = val;
    pcf->m_livestemc_xfer_to_litter = val;
    pcf->m_deadstemc_xfer_to_litter = val;
    pcf->m_livecrootc_xfer_to_litter = val;
    pcf->m_deadcrootc_xfer_to_litter = val;
    pcf->m_livestemc_to_litter = val;
    pcf->m_deadstemc_to_litter = val;
    pcf->m_livecrootc_to_litter = val;
    pcf->m_deadcrootc_to_litter = val;
    pcf->m_gresp_storage_to_litter = val;
    pcf->m_gresp_xfer_to_litter = val;
    pcf->hrv_leafc_to_litter = val             ;
    pcf->hrv_leafc_storage_to_litter = val     ;
    pcf->hrv_leafc_xfer_to_litter = val        ;
    pcf->hrv_frootc_to_litter = val            ;
    pcf->hrv_frootc_storage_to_litter = val    ;
    pcf->hrv_frootc_xfer_to_litter = val       ;
    pcf->hrv_livestemc_to_litter = val         ;
    pcf->hrv_livestemc_storage_to_litter = val ;
    pcf->hrv_livestemc_xfer_to_litter = val    ;
    pcf->hrv_deadstemc_to_prod10c = val        ;
    pcf->hrv_deadstemc_to_prod100c = val       ;
    pcf->hrv_deadstemc_storage_to_litter = val ;
    pcf->hrv_deadstemc_xfer_to_litter = val    ;
    pcf->hrv_livecrootc_to_litter = val        ;
    pcf->hrv_livecrootc_storage_to_litter = val;
    pcf->hrv_livecrootc_xfer_to_litter = val   ;
    pcf->hrv_deadcrootc_to_litter = val        ;
    pcf->hrv_deadcrootc_storage_to_litter = val;
    pcf->hrv_deadcrootc_xfer_to_litter = val   ;
    pcf->hrv_gresp_storage_to_litter = val     ;
    pcf->hrv_gresp_xfer_to_litter = val        ;
    pcf->hrv_xsmrpool_to_atm = val                 ;
    pcf->m_leafc_to_fire = val;
    pcf->m_frootc_to_fire = val;
    pcf->m_leafc_storage_to_fire = val;
    pcf->m_frootc_storage_to_fire = val;
    pcf->m_livestemc_storage_to_fire = val;
    pcf->m_deadstemc_storage_to_fire = val;
    pcf->m_livecrootc_storage_to_fire = val;
    pcf->m_deadcrootc_storage_to_fire = val;
    pcf->m_leafc_xfer_to_fire = val;
    pcf->m_frootc_xfer_to_fire = val;
    pcf->m_livestemc_xfer_to_fire = val;
    pcf->m_deadstemc_xfer_to_fire = val;
    pcf->m_livecrootc_xfer_to_fire = val;
    pcf->m_deadcrootc_xfer_to_fire = val;
    pcf->m_livestemc_to_fire = val;
    pcf->m_deadstemc_to_fire = val;
    pcf->m_deadstemc_to_litter_fire = val;
    pcf->m_livecrootc_to_fire = val;
    pcf->m_deadcrootc_to_fire = val;
    pcf->m_deadcrootc_to_litter_fire = val;
    pcf->m_gresp_storage_to_fire = val;
    pcf->m_gresp_xfer_to_fire = val;
    pcf->leafc_xfer_to_leafc = val;
    pcf->frootc_xfer_to_frootc = val;
    pcf->livestemc_xfer_to_livestemc = val;
    pcf->deadstemc_xfer_to_deadstemc = val;
    pcf->livecrootc_xfer_to_livecrootc = val;
    pcf->deadcrootc_xfer_to_deadcrootc = val;
    pcf->leafc_to_litter = val;
    pcf->frootc_to_litter = val;
    pcf->leaf_mr = val;
    pcf->froot_mr = val;
    pcf->livestem_mr = val;
    pcf->livecroot_mr = val;
    pcf->leaf_curmr = val;
    pcf->froot_curmr = val;
    pcf->livestem_curmr = val;
    pcf->livecroot_curmr = val;
    pcf->leaf_xsmr = val;
    pcf->froot_xsmr = val;
    pcf->livestem_xsmr = val;
    pcf->livecroot_xsmr = val;
    pcf->psnsun_to_cpool = val;
    pcf->psnshade_to_cpool = val;
    pcf->cpool_to_xsmrpool = val;
    pcf->cpool_to_leafc = val;
    pcf->cpool_to_leafc_storage = val;
    pcf->cpool_to_frootc = val;
    pcf->cpool_to_frootc_storage = val;
    pcf->cpool_to_livestemc = val;
    pcf->cpool_to_livestemc_storage = val;
    pcf->cpool_to_deadstemc = val;
    pcf->cpool_to_deadstemc_storage = val;
    pcf->cpool_to_livecrootc = val;
    pcf->cpool_to_livecrootc_storage = val;
    pcf->cpool_to_deadcrootc = val;
    pcf->cpool_to_deadcrootc_storage = val;
    pcf->cpool_to_gresp_storage = val;
    pcf->cpool_leaf_gr = val;
    pcf->cpool_leaf_storage_gr = val;
    pcf->transfer_leaf_gr = val;
    pcf->cpool_froot_gr = val;
    pcf->cpool_froot_storage_gr = val;
    pcf->transfer_froot_gr = val;
    pcf->cpool_livestem_gr = val;
    pcf->cpool_livestem_storage_gr = val;
    pcf->transfer_livestem_gr = val;
    pcf->cpool_deadstem_gr = val;
    pcf->cpool_deadstem_storage_gr = val;
    pcf->transfer_deadstem_gr = val;
    pcf->cpool_livecroot_gr = val;
    pcf->cpool_livecroot_storage_gr = val;
    pcf->transfer_livecroot_gr = val;
    pcf->cpool_deadcroot_gr = val;
    pcf->cpool_deadcroot_storage_gr = val;
    pcf->transfer_deadcroot_gr = val;
    pcf->leafc_storage_to_xfer = val;
    pcf->frootc_storage_to_xfer = val;
    pcf->livestemc_storage_to_xfer = val;
    pcf->deadstemc_storage_to_xfer = val;
    pcf->livecrootc_storage_to_xfer = val;
    pcf->deadcrootc_storage_to_xfer = val;
    pcf->gresp_storage_to_xfer = val;
    pcf->livestemc_to_deadstemc = val;
    pcf->livecrootc_to_deadcrootc = val;
    pcf->gpp = val;
    pcf->mr = val;
    pcf->current_gr = val;
    pcf->transfer_gr = val;
    pcf->storage_gr = val;
    pcf->gr = val;
    pcf->ar = val;
    pcf->rr = val;
    pcf->npp = val;
    pcf->agnpp = val;
    pcf->bgnpp = val;
    pcf->litfall = val;
    pcf->vegfire = val;
    pcf->wood_harvestc = val;
    pcf->pft_cinputs = val;
    pcf->pft_coutputs = val;
    pcf->pft_fire_closs = val;
    pcf->frootc_alloc = val;
    pcf->frootc_loss = val;
    pcf->leafc_alloc = val;
    pcf->leafc_loss = val;
    pcf->woodc_alloc = val;
    pcf->woodc_loss = val;
}

void CNSetPnf(double val, pft_nflux_type *pnf)
{

    pnf->m_leafn_to_litter = val;
    pnf->m_frootn_to_litter = val;
    pnf->m_leafn_storage_to_litter = val;
    pnf->m_frootn_storage_to_litter = val;
    pnf->m_livestemn_storage_to_litter = val;
    pnf->m_deadstemn_storage_to_litter = val;
    pnf->m_livecrootn_storage_to_litter = val;
    pnf->m_deadcrootn_storage_to_litter = val;
    pnf->m_leafn_xfer_to_litter = val;
    pnf->m_frootn_xfer_to_litter = val;
    pnf->m_livestemn_xfer_to_litter = val;
    pnf->m_deadstemn_xfer_to_litter = val;
    pnf->m_livecrootn_xfer_to_litter = val;
    pnf->m_deadcrootn_xfer_to_litter = val;
    pnf->m_livestemn_to_litter = val;
    pnf->m_deadstemn_to_litter = val;
    pnf->m_livecrootn_to_litter = val;
    pnf->m_deadcrootn_to_litter = val;
    pnf->m_retransn_to_litter = val;
    pnf->hrv_leafn_to_litter = val             ;
    pnf->hrv_frootn_to_litter = val            ;
    pnf->hrv_leafn_storage_to_litter = val     ;
    pnf->hrv_frootn_storage_to_litter = val    ;
    pnf->hrv_livestemn_storage_to_litter = val ;
    pnf->hrv_deadstemn_storage_to_litter = val ;
    pnf->hrv_livecrootn_storage_to_litter = val;
    pnf->hrv_deadcrootn_storage_to_litter = val;
    pnf->hrv_leafn_xfer_to_litter = val        ;
    pnf->hrv_frootn_xfer_to_litter = val       ;
    pnf->hrv_livestemn_xfer_to_litter = val    ;
    pnf->hrv_deadstemn_xfer_to_litter = val    ;
    pnf->hrv_livecrootn_xfer_to_litter = val   ;
    pnf->hrv_deadcrootn_xfer_to_litter = val   ;
    pnf->hrv_livestemn_to_litter = val         ;
    pnf->hrv_deadstemn_to_prod10n = val        ;
    pnf->hrv_deadstemn_to_prod100n = val       ;
    pnf->hrv_livecrootn_to_litter = val        ;
    pnf->hrv_deadcrootn_to_litter = val        ;
    pnf->hrv_retransn_to_litter = val           ;
    pnf->m_leafn_to_fire = val;
    pnf->m_frootn_to_fire = val;
    pnf->m_leafn_storage_to_fire = val;
    pnf->m_frootn_storage_to_fire = val;
    pnf->m_livestemn_storage_to_fire = val;
    pnf->m_deadstemn_storage_to_fire = val;
    pnf->m_livecrootn_storage_to_fire = val;
    pnf->m_deadcrootn_storage_to_fire = val;
    pnf->m_leafn_xfer_to_fire = val;
    pnf->m_frootn_xfer_to_fire = val;
    pnf->m_livestemn_xfer_to_fire = val;
    pnf->m_deadstemn_xfer_to_fire = val;
    pnf->m_livecrootn_xfer_to_fire = val;
    pnf->m_deadcrootn_xfer_to_fire = val;
    pnf->m_livestemn_to_fire = val;
    pnf->m_deadstemn_to_fire = val;
    pnf->m_deadstemn_to_litter_fire = val;
    pnf->m_livecrootn_to_fire = val;
    pnf->m_deadcrootn_to_fire = val;
    pnf->m_deadcrootn_to_litter_fire = val;
    pnf->m_retransn_to_fire = val;
    pnf->leafn_xfer_to_leafn = val;
    pnf->frootn_xfer_to_frootn = val;
    pnf->livestemn_xfer_to_livestemn = val;
    pnf->deadstemn_xfer_to_deadstemn = val;
    pnf->livecrootn_xfer_to_livecrootn = val;
    pnf->deadcrootn_xfer_to_deadcrootn = val;
    pnf->leafn_to_litter = val;
    pnf->leafn_to_retransn = val;
    pnf->frootn_to_litter = val;
    pnf->retransn_to_npool = val;
    pnf->sminn_to_npool = val;
    pnf->npool_to_leafn = val;
    pnf->npool_to_leafn_storage = val;
    pnf->npool_to_frootn = val;
    pnf->npool_to_frootn_storage = val;
    pnf->npool_to_livestemn = val;
    pnf->npool_to_livestemn_storage = val;
    pnf->npool_to_deadstemn = val;
    pnf->npool_to_deadstemn_storage = val;
    pnf->npool_to_livecrootn = val;
    pnf->npool_to_livecrootn_storage = val;
    pnf->npool_to_deadcrootn = val;
    pnf->npool_to_deadcrootn_storage = val;
    pnf->leafn_storage_to_xfer = val;
    pnf->frootn_storage_to_xfer = val;
    pnf->livestemn_storage_to_xfer = val;
    pnf->deadstemn_storage_to_xfer = val;
    pnf->livecrootn_storage_to_xfer = val;
    pnf->deadcrootn_storage_to_xfer = val;
    pnf->livestemn_to_deadstemn = val;
    pnf->livestemn_to_retransn = val;
    pnf->livecrootn_to_deadcrootn = val;
    pnf->livecrootn_to_retransn = val;
    pnf->ndeploy = val;
    pnf->pft_ninputs = val;
    pnf->pft_noutputs = val;
    pnf->wood_harvestn = val;
    pnf->pft_fire_nloss = val;
}

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetCps
!
! !INTERFACE:
subroutine CNSetCps(num, filter, val, cps)
!
! !DESCRIPTION:
! Set column physical state variables
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (column_pstate_type), intent(inout) :: cps
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i,j     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      cps%decl(i) = val
      cps%coszen(i) = val
      cps%fpi(i) = val
      cps%fpg(i) = val
      cps%annsum_counter(i) = val
      cps%cannsum_npp(i) = val
      cps%cannavg_t2m(i) = val
      cps%wf(i) = val
      cps%me(i) = val
      cps%fire_prob(i) = val
      cps%mean_fire_prob(i) = val
      cps%fireseasonl(i) = val
      cps%farea_burned(i) = val
      cps%ann_farea_burned(i) = val
   end do

   do j = 1,nlevgrnd
      do fi = 1,num
         i = filter(fi)
         cps%bsw2(i,j) = val
         cps%psisat(i,j) = val
         cps%vwcsat(i,j) = val
         cps%soilpsi(i,j) = val
      end do
   end do

end subroutine CNSetCps
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetCcs
!
! !INTERFACE:
subroutine CNSetCcs(num, filter, val, ccs)
!
! !DESCRIPTION:
! Set column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (column_cstate_type), intent(inout) :: ccs
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      ccs%cwdc(i) = val
      ccs%litr1c(i) = val
      ccs%litr2c(i) = val
      ccs%litr3c(i) = val
      ccs%soil1c(i) = val
      ccs%soil2c(i) = val
      ccs%soil3c(i) = val
      ccs%soil4c(i) = val
      ccs%col_ctrunc(i) = val
      ccs%totlitc(i) = val
      ccs%totsomc(i) = val
      ccs%totecosysc(i) = val
      ccs%totcolc(i) = val

   end do

end subroutine CNSetCcs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetCns
!
! !INTERFACE:
subroutine CNSetCns(num, filter, val, cns)
!
! !DESCRIPTION:
! Set column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (column_nstate_type), intent(inout) :: cns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      cns%cwdn(i) = val
      cns%litr1n(i) = val
      cns%litr2n(i) = val
      cns%litr3n(i) = val
      cns%soil1n(i) = val
      cns%soil2n(i) = val
      cns%soil3n(i) = val
      cns%soil4n(i) = val
      cns%sminn(i) = val
      cns%col_ntrunc(i) = val
      cns%totlitn(i) = val
      cns%totsomn(i) = val
      cns%totecosysn(i) = val
      cns%totcoln(i) = val
   end do

end subroutine CNSetCns
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetCcf
!
! !INTERFACE:
subroutine CNSetCcf(num, filter, val, ccf)
!
! !DESCRIPTION:
! Set column carbon flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (column_cflux_type), intent(inout) :: ccf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      ccf%m_leafc_to_litr1c(i)                = val
      ccf%m_leafc_to_litr2c(i)                = val
      ccf%m_leafc_to_litr3c(i)                = val
      ccf%m_frootc_to_litr1c(i)               = val
      ccf%m_frootc_to_litr2c(i)               = val
      ccf%m_frootc_to_litr3c(i)               = val
      ccf%m_leafc_storage_to_litr1c(i)        = val
      ccf%m_frootc_storage_to_litr1c(i)       = val
      ccf%m_livestemc_storage_to_litr1c(i)    = val
      ccf%m_deadstemc_storage_to_litr1c(i)    = val
      ccf%m_livecrootc_storage_to_litr1c(i)   = val
      ccf%m_deadcrootc_storage_to_litr1c(i)   = val
      ccf%m_leafc_xfer_to_litr1c(i)           = val
      ccf%m_frootc_xfer_to_litr1c(i)          = val
      ccf%m_livestemc_xfer_to_litr1c(i)       = val
      ccf%m_deadstemc_xfer_to_litr1c(i)       = val
      ccf%m_livecrootc_xfer_to_litr1c(i)      = val
      ccf%m_deadcrootc_xfer_to_litr1c(i)      = val
      ccf%m_livestemc_to_cwdc(i)              = val
      ccf%m_deadstemc_to_cwdc(i)              = val
      ccf%m_livecrootc_to_cwdc(i)             = val
      ccf%m_deadcrootc_to_cwdc(i)             = val
      ccf%m_gresp_storage_to_litr1c(i)        = val
      ccf%m_gresp_xfer_to_litr1c(i)           = val
      ccf%hrv_leafc_to_litr1c(i)              = val             
      ccf%hrv_leafc_to_litr2c(i)              = val             
      ccf%hrv_leafc_to_litr3c(i)              = val             
      ccf%hrv_frootc_to_litr1c(i)             = val            
      ccf%hrv_frootc_to_litr2c(i)             = val            
      ccf%hrv_frootc_to_litr3c(i)             = val            
      ccf%hrv_livestemc_to_cwdc(i)            = val           
      ccf%hrv_deadstemc_to_prod10c(i)         = val        
      ccf%hrv_deadstemc_to_prod100c(i)        = val       
      ccf%hrv_livecrootc_to_cwdc(i)           = val          
      ccf%hrv_deadcrootc_to_cwdc(i)           = val          
      ccf%hrv_leafc_storage_to_litr1c(i)      = val     
      ccf%hrv_frootc_storage_to_litr1c(i)     = val    
      ccf%hrv_livestemc_storage_to_litr1c(i)  = val 
      ccf%hrv_deadstemc_storage_to_litr1c(i)  = val 
      ccf%hrv_livecrootc_storage_to_litr1c(i) = val
      ccf%hrv_deadcrootc_storage_to_litr1c(i) = val
      if ( crop_prog )then
         ccf%livestemc_to_litr1c(i) = val
         ccf%livestemc_to_litr2c(i) = val
         ccf%livestemc_to_litr3c(i) = val
         ccf%grainc_to_litr1c(i)    = val
         ccf%grainc_to_litr2c(i)    = val
         ccf%grainc_to_litr3c(i)    = val
      end if
      ccf%hrv_gresp_storage_to_litr1c(i)      = val
      ccf%hrv_leafc_xfer_to_litr1c(i)         = val
      ccf%hrv_frootc_xfer_to_litr1c(i)        = val
      ccf%hrv_livestemc_xfer_to_litr1c(i)     = val
      ccf%hrv_deadstemc_xfer_to_litr1c(i)     = val
      ccf%hrv_livecrootc_xfer_to_litr1c(i)    = val
      ccf%hrv_deadcrootc_xfer_to_litr1c(i)    = val
      ccf%hrv_gresp_xfer_to_litr1c(i)         = val
      ccf%m_deadstemc_to_cwdc_fire(i)         = val
      ccf%m_deadcrootc_to_cwdc_fire(i)        = val
      ccf%m_litr1c_to_fire(i)                 = val
      ccf%m_litr2c_to_fire(i)                 = val
      ccf%m_litr3c_to_fire(i)                 = val
      ccf%m_cwdc_to_fire(i)                   = val
      ccf%prod10c_loss(i)                     = val
      ccf%prod100c_loss(i)                    = val
      ccf%product_closs(i)                    = val
      ccf%leafc_to_litr1c(i)                  = val
      ccf%leafc_to_litr2c(i)                  = val
      ccf%leafc_to_litr3c(i)                  = val
      ccf%frootc_to_litr1c(i)                 = val
      ccf%frootc_to_litr2c(i)                 = val
      ccf%frootc_to_litr3c(i)                 = val
      ccf%cwdc_to_litr2c(i)                   = val
      ccf%cwdc_to_litr3c(i)                   = val
      ccf%litr1_hr(i)                         = val
      ccf%litr1c_to_soil1c(i)                 = val
      ccf%litr2_hr(i)                         = val
      ccf%litr2c_to_soil2c(i)                 = val
      ccf%litr3_hr(i)                         = val
      ccf%litr3c_to_soil3c(i)                 = val
      ccf%soil1_hr(i)                         = val
      ccf%soil1c_to_soil2c(i)                 = val
      ccf%soil2_hr(i)                         = val
      ccf%soil2c_to_soil3c(i)                 = val
      ccf%soil3_hr(i)                         = val
      ccf%soil3c_to_soil4c(i)                 = val
      ccf%soil4_hr(i)                         = val
      ccf%lithr(i)                            = val
      ccf%somhr(i)                            = val
      ccf%hr(i)                               = val
      ccf%sr(i)                               = val
      ccf%er(i)                               = val
      ccf%litfire(i)                          = val
      ccf%somfire(i)                          = val
      ccf%totfire(i)                          = val
      ccf%nep(i)                              = val
      ccf%nbp(i)                              = val
      ccf%nee(i)                              = val
      ccf%col_cinputs(i)                      = val
      ccf%col_coutputs(i)                     = val
      ccf%col_fire_closs(i)                   = val
      ccf%cwdc_hr(i)                          = val
      ccf%cwdc_loss(i)                        = val
      ccf%litterc_loss(i)                     = val

  end do

end subroutine CNSetCcf
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSetCnf
!
! !INTERFACE:
subroutine CNSetCnf(num, filter, val, cnf)
!
! !DESCRIPTION:
! Set column nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num
    integer , intent(in) :: filter(:)
    real(r8), intent(in) :: val
    type (column_nflux_type), intent(inout) :: cnf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in/out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: fi,i     ! loop index
!EOP
!------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      cnf%ndep_to_sminn(i) = val
      cnf%nfix_to_sminn(i) = val
      cnf%m_leafn_to_litr1n(i) = val
      cnf%m_leafn_to_litr2n(i) = val
      cnf%m_leafn_to_litr3n(i) = val
      cnf%m_frootn_to_litr1n(i) = val
      cnf%m_frootn_to_litr2n(i) = val
      cnf%m_frootn_to_litr3n(i) = val
      cnf%m_leafn_storage_to_litr1n(i) = val
      cnf%m_frootn_storage_to_litr1n(i) = val
      cnf%m_livestemn_storage_to_litr1n(i) = val
      cnf%m_deadstemn_storage_to_litr1n(i) = val
      cnf%m_livecrootn_storage_to_litr1n(i) = val
      cnf%m_deadcrootn_storage_to_litr1n(i) = val
      cnf%m_leafn_xfer_to_litr1n(i) = val
      cnf%m_frootn_xfer_to_litr1n(i) = val
      cnf%m_livestemn_xfer_to_litr1n(i) = val
      cnf%m_deadstemn_xfer_to_litr1n(i) = val
      cnf%m_livecrootn_xfer_to_litr1n(i) = val
      cnf%m_deadcrootn_xfer_to_litr1n(i) = val
      cnf%m_livestemn_to_cwdn(i) = val
      cnf%m_deadstemn_to_cwdn(i) = val
      cnf%m_livecrootn_to_cwdn(i) = val
      cnf%m_deadcrootn_to_cwdn(i) = val
      cnf%m_retransn_to_litr1n(i) = val
      cnf%hrv_leafn_to_litr1n(i) = val             
      cnf%hrv_leafn_to_litr2n(i) = val             
      cnf%hrv_leafn_to_litr3n(i) = val             
      cnf%hrv_frootn_to_litr1n(i) = val            
      cnf%hrv_frootn_to_litr2n(i) = val            
      cnf%hrv_frootn_to_litr3n(i) = val            
      cnf%hrv_livestemn_to_cwdn(i) = val           
      cnf%hrv_deadstemn_to_prod10n(i) = val        
      cnf%hrv_deadstemn_to_prod100n(i) = val       
      cnf%hrv_livecrootn_to_cwdn(i) = val          
      cnf%hrv_deadcrootn_to_cwdn(i) = val          
      cnf%hrv_retransn_to_litr1n(i) = val          
      cnf%hrv_leafn_storage_to_litr1n(i) = val     
      cnf%hrv_frootn_storage_to_litr1n(i) = val    
      cnf%hrv_livestemn_storage_to_litr1n(i) = val 
      cnf%hrv_deadstemn_storage_to_litr1n(i) = val 
      cnf%hrv_livecrootn_storage_to_litr1n(i) = val
      cnf%hrv_deadcrootn_storage_to_litr1n(i) = val
      cnf%hrv_leafn_xfer_to_litr1n(i) = val        
      cnf%hrv_frootn_xfer_to_litr1n(i) = val       
      cnf%hrv_livestemn_xfer_to_litr1n(i) = val    
      cnf%hrv_deadstemn_xfer_to_litr1n(i) = val    
      cnf%hrv_livecrootn_xfer_to_litr1n(i) = val   
      cnf%hrv_deadcrootn_xfer_to_litr1n(i) = val   
      cnf%m_deadstemn_to_cwdn_fire(i) = val
      cnf%m_deadcrootn_to_cwdn_fire(i) = val
      cnf%m_litr1n_to_fire(i) = val
      cnf%m_litr2n_to_fire(i) = val
      cnf%m_litr3n_to_fire(i) = val
      cnf%m_cwdn_to_fire(i) = val
      cnf%prod10n_loss(i) = val
      cnf%prod100n_loss(i) = val
      cnf%product_nloss(i) = val
      if ( crop_prog )then
         cnf%grainn_to_litr1n(i)    = val
         cnf%grainn_to_litr2n(i)    = val
         cnf%grainn_to_litr3n(i)    = val
         cnf%livestemn_to_litr1n(i) = val
         cnf%livestemn_to_litr2n(i) = val
         cnf%livestemn_to_litr3n(i) = val
      end if
      cnf%leafn_to_litr1n(i) = val
      cnf%leafn_to_litr2n(i) = val
      cnf%leafn_to_litr3n(i) = val
      cnf%frootn_to_litr1n(i) = val
      cnf%frootn_to_litr2n(i) = val
      cnf%frootn_to_litr3n(i) = val
      cnf%cwdn_to_litr2n(i) = val
      cnf%cwdn_to_litr3n(i) = val
      cnf%litr1n_to_soil1n(i) = val
      cnf%sminn_to_soil1n_l1(i) = val
      cnf%litr2n_to_soil2n(i) = val
      cnf%sminn_to_soil2n_l2(i) = val
      cnf%litr3n_to_soil3n(i) = val
      cnf%sminn_to_soil3n_l3(i) = val
      cnf%soil1n_to_soil2n(i) = val
      cnf%sminn_to_soil2n_s1(i) = val
      cnf%soil2n_to_soil3n(i) = val
      cnf%sminn_to_soil3n_s2(i) = val
      cnf%soil3n_to_soil4n(i) = val
      cnf%sminn_to_soil4n_s3(i) = val
      cnf%soil4n_to_sminn(i) = val
      cnf%sminn_to_denit_l1s1(i) = val
      cnf%sminn_to_denit_l2s2(i) = val
      cnf%sminn_to_denit_l3s3(i) = val
      cnf%sminn_to_denit_s1s2(i) = val
      cnf%sminn_to_denit_s2s3(i) = val
      cnf%sminn_to_denit_s3s4(i) = val
      cnf%sminn_to_denit_s4(i) = val
      cnf%sminn_to_denit_excess(i) = val
      cnf%sminn_leached(i) = val
      cnf%potential_immob(i) = val
      cnf%actual_immob(i) = val
      cnf%sminn_to_plant(i) = val
      cnf%supplement_to_sminn(i) = val
      cnf%gross_nmin(i) = val
      cnf%net_nmin(i) = val
      cnf%denit(i) = val
      cnf%col_ninputs(i) = val
      cnf%col_noutputs(i) = val
      cnf%col_fire_nloss(i) = val
   end do

end subroutine CNSetCnf
!-----------------------------------------------------------------------

end module CNSetValueMod
