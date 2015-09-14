module fv_phys_mod

use constants_mod,         only: grav, rdgas, rvgas, pi, cp_air, cp_vapor, hlv, radius, kappa
use time_manager_mod,      only: time_type, get_time
use lin_cld_microphys_mod, only: lin_cld_microphys_driver, sg_conv, qsmith
use hswf_mod,              only: Held_Suarez_Tend
use fv_sg_mod,             only: fv_dry_conv
use fv_update_phys_mod,    only: fv_update_phys
use fv_timing_mod,         only: timing_on, timing_off
use monin_obukhov_mod,     only: mon_obkv
use tracer_manager_mod,    only: get_tracer_index
use field_manager_mod,     only: MODEL_ATMOS
use fms_mod,               only: error_mesg, FATAL, file_exist, open_namelist_file,  &
                                 check_nml_error, mpp_pe, mpp_root_pe, close_file, &
                                 write_version_number, stdlog, mpp_error
use fv_mp_mod,             only: is_master, mp_reduce_max
use fv_diagnostics_mod,    only: prt_maxmin

use fv_arrays_mod,          only: fv_grid_type, fv_flags_type, fv_nest_type, fv_grid_bounds_type
use mpp_domains_mod,       only: domain2d
use diag_manager_mod,      only: register_diag_field, register_static_field, send_data

implicit none

  logical:: nudge_initialized = .false.
  logical:: sim_phys_initialized = .false.
  logical:: g0_sum_initialized = .false.
  real, allocatable, dimension(:,:) :: l_area
  logical:: master
  real, allocatable:: u0(:,:,:), v0(:,:,:), t0(:,:,:), dp(:,:,:), u_star(:,:)
  integer:: id_rain_k, id_vr_k
  real, allocatable:: prec_total(:,:)
  real:: global_area = -1.
  real:: missing_value = -1.e10

public :: fv_phys, fv_nudge

!---- version number -----
  character(len=128) :: version = '$Id: fv_phys.F90,v 20.0 2013/12/13 23:04:08 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  integer:: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, cld_amt
! For nudging the IC to a steady state:
  real :: tau_winds = 25.
  real :: tau_temp  = -1.
  real :: tau_press = -1.

  real, parameter:: c_liq = 4190.       ! heat capacity of water at 0C
  real, parameter:: c_ice = 2106.       ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
  real, parameter:: cp_vap = cp_vapor   ! 1846.
! For consistency, cv_vap derived FMS constants:
  real, parameter:: cv_vap = cp_vap - rvgas  ! 1384.5
  real, parameter:: cv_air = cp_air - rdgas
  real, parameter:: dc_vap = cv_vap - c_liq     ! = -2368.
  real, parameter:: hlv0 = 2.501e6   ! Emanual Appendix-2
  real, parameter:: Lv0 =  hlv0 - dc_vap*273.15   ! = 3.147782e6

  real::  zvir
  real:: tau_difz = 600.
  real:: area_ref = 6.25E8
  real:: solar_constant = 1367.
! Namelist for Sim_Phys
  logical:: do_k_momentum = .true.
  logical:: do_sedi_w = .true.
  logical:: diurnal_cycle    = .false.
  logical:: mixed_layer      = .false.
  logical:: gray_rad         = .false.
  logical:: strat_rad        = .false.
  logical:: do_abl = .false.
  logical:: do_mon_obkv = .true.
  logical:: do_lin_microphys = .false.
  logical:: do_k_warm_rain = .false.
  logical:: do_strat_forcing = .true.
  logical:: do_sg_conv       = .false.
  logical:: prog_cloud       = .true.
  logical:: zero_winds       = .false.  ! use this only for the doubly periodic domain
  logical:: do_reed_phys     = .false.
  logical:: reed_cond_only   = .false.
  real:: tau_zero  = 1.   ! time-scale to "zero" domain-averaged winds (days)
  real:: s_fac  =  0.1
  real:: heating_rate =  0.5    ! deg K per day, stratosphere, for gray-radiation
  real:: cooling_rate = -1.0    ! deg K per day
  real:: c0 = 1.E7              ! Ocean heat capabicity
  real:: low_c = 0.3            ! global mean *low* cloud fraction
  real:: sw_abs = 0.            ! fraction of the solar absorbed/reflected by the atm
  logical:: uniform_sst = .true.
  real:: sst0 = 300.    ! 302.15 for DCIMP case-51
  real:: shift_n = 0.
  logical:: do_t_strat = .false.
  real:: p_strat = 50.E2
  real:: t_strat = 200.
  real:: tau_strat = 10.     ! days
  real:: t_fac = 1.
  integer:: print_freq = 86400   ! sec

namelist /sim_phys_nml/mixed_layer, gray_rad, strat_rad, do_lin_microphys,   &
                       heating_rate, cooling_rate, uniform_sst, sst0, c0,    &
                       sw_abs, prog_cloud, low_c, do_sg_conv, diurnal_cycle, &
                       do_mon_obkv, do_t_strat, p_strat, t_strat, tau_strat, &
                       do_abl, t_fac, tau_difz, do_strat_forcing, s_fac,     &
                       shift_n, print_freq, zero_winds, tau_zero, tau_winds, &
                       tau_temp, tau_press, do_K_warm_rain, do_k_momentum,   &
                       do_sedi_w,  do_reed_phys, reed_cond_only

    integer :: id_rain, id_dqdt, id_dTdt, id_dudt, id_dvdt

contains

 subroutine fv_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                    u, v, w, pt, q, pe, delp, peln, pkz, pdt,    &
                    ua, va, phis, grid, ptop, ak, bk, ks, ps, pk,&
                    u_srf, v_srf, ts, delz, hydrostatic,         &
                    oro, rayf, p_ref, fv_sg_adj,                 &
                    do_Held_Suarez, gridstruct, flagstruct,      &
                    neststruct, nwat, bd, domain, Time, time_total)

    integer, INTENT(IN   ) :: npx, npy, npz
    integer, INTENT(IN   ) :: is, ie, js, je, ng, nq, nwat
    integer, INTENT(IN   ) :: fv_sg_adj
    real, INTENT(IN) :: p_ref, ptop
    real, INTENT(IN) :: oro(is:ie,js:je)

    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
    real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
    real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
    real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)
    real   , INTENT(INOUT) ::  pk (is   :ie     ,js   :je     ,npz+1)
    real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   ps(is-ng:ie+ng,js-ng:je+ng)

    real   , INTENT(IN   ) :: phis(is-ng:ie+ng,js-ng:je+ng)
    real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 2)
    real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
    integer, INTENT(IN   ) :: ks

    real   , INTENT(IN   )    :: pdt
    logical, INTENT(IN   )    :: rayf, do_Held_Suarez
    real, INTENT(inout):: u_srf(is:ie,js:je)
    real, INTENT(inout):: v_srf(is:ie,js:je)
    real, INTENT(inout)::    ts(is:ie,js:je)

    type(fv_grid_type) :: gridstruct
    type(fv_flags_type) :: flagstruct
    type(fv_nest_type) :: neststruct
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(domain2d), intent(INOUT) :: domain

    type (time_type), intent(in) :: Time
    real, INTENT(IN), optional:: time_total
    logical, intent(in) ::  hydrostatic
    real, intent(inout) ::  delz(is-ng:ie+ng,js-ng:je+ng,npz)
    real, allocatable:: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:), q_dt(:,:,:,:)
    real, allocatable :: rain(:,:)
    real, allocatable ::   pz(:,:), rdelp(:,:)
    real :: qdiag(1,1,1)
    logical moist_phys
    integer  isd, ied, jsd, jed
    integer  i, j, k, n

    real rdg
    logical :: first_step = .true.
    logical used

    master = is_master()
!   if (.not. sim_phys_initialized) call fv_phys_init(nwat, time)

    isd = is-ng;   ied = ie + ng
    jsd = js-ng;   jed = je + ng

     allocate ( u_dt(isd:ied,jsd:jed,npz) )
     allocate ( v_dt(isd:ied,jsd:jed,npz) )
     allocate ( t_dt(is:ie,js:je,npz) )
     allocate ( q_dt(is:ie,js:je,npz,nq) )

! Place the memory in the optimal shared mem space!
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,is,ie,js,je,nq,u_dt,v_dt,t_dt,q_dt) 
     do k=1, npz
        do j=jsd, jed
           do i=isd, ied
              u_dt(i,j,k) = 0.
              v_dt(i,j,k) = 0.
           enddo
        enddo

        do j=js, je
           do i=is, ie
              t_dt(i,j,k) = 0.
           enddo
        enddo
        do n=1,nq
        do j=js, je
           do i=is, ie
              q_dt(i,j,k,n) = 0.
           enddo
        enddo
        enddo
     enddo

     if ( fv_sg_adj > 0 ) then
          call fv_dry_conv(isd, ied, jsd, jed, is, ie, js, je, npz, nq, pdt,  &
                           fv_sg_adj, nwat, delp, pe, peln, pkz, pt, q, ua, va,  &
                           hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
     endif

    if( do_Held_Suarez ) then
       moist_phys = .false.
       call Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, grid,   &
                              delz, phis, hydrostatic, ak, bk, ks,    &
                              do_strat_forcing, rayf, master, Time, time_total)
    elseif (do_reed_phys) then

       moist_phys = .true.
       if (.not. hydrostatic) rdg = -rdgas / grav
       if (.not. allocated(pz)) allocate(pz(is:ie,npz))
       if (.not. allocated(rdelp)) allocate(rdelp(is:ie,npz))
       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
       zvir = rvgas/rdgas - 1.

       do j=js,je

          if (hydrostatic) then

             do k=1,npz
                do i=is,ie
                   !pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                   pz(i,k) = (pe(i,k+1,j)-pe(i,k,j))/(peln(i,k+1,j)-peln(i,k,j))
                enddo
             enddo

          else
             do k=1,npz
             do i=is,ie
				!Reed physics has no condensate
                pz(i,k) = rdg*delp(i,j,k)*pt(i,j,k)*    &
                     (1.+zvir*q(i,j,k,sphum))/delz(i,j,k)
!!$                pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!!$                     (1.+dp1(i,j,k))/delz(i,j,k)) )
             enddo
             enddo
                
          endif

          if (master .and. first_step) then
             do k=1,npz
                print*, pe(i,k,j), pz(is,k)
             enddo
             print*, pe(i,npz+1,j)
             first_step = .false.
          endif

          do k=1,npz
          do i=is,ie
             rdelp(i,k) = 1./delp(i,j,k)
          enddo
          enddo

          if (.not. allocated(rain)) allocate(rain(is:ie,js:je))
!!$          call reed_simple_physics_v5(ie-is+1, npz, pdt, gridstruct%agrid(is:ie,j,2), pt(is:ie,j,:), q(is:ie,j,:,sphum), &
!!$               ua(is:ie,j,:), va(is:ie,j,:), pz, pe(is:ie,:,j), delp(is:ie,j,:), &
!!$               rdelp, ps(is:ie,j), rain(is:ie,j), 0)
          call reed_simple_physics(ie-is+1, npz, pdt, pt(is:ie,j,:), q(is:ie,j,:,sphum), &
               ua(is:ie,j,:), va(is:ie,j,:), pz, pe(is:ie,:,j), delp(is:ie,j,:), &
               rdelp, pe(is:ie,npz+1,j), rain(is:ie,j), reed_cond_only, & 
               q_dt(is:ie,j,:,sphum), t_dt(is:ie,j,:), u_dt(is:ie,j,:), v_dt(is:ie,j,:))

       enddo

       if (id_rain > 0) then
          used=send_data(id_rain, rain, time)
       endif
    else
       moist_phys = .true.
                                             call timing_on('SIM_PHYS')
       call sim_phys(npx, npy, npz, is, ie, js, je, ng, nq, nwat, pk, pkz, &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va, pt, delz, q, &
                     pe, delp, peln, ts, oro, hydrostatic, pdt, grid, ak, bk, &
                     p_ref, fv_sg_adj, Time, time_total, gridstruct)
                                            call timing_off('SIM_PHYS')
    endif

    if (id_dudt > 0) then
       used=send_data(id_dudt, u_dt(is:ie,js:je,npz), time)
    endif
    if (id_dvdt > 0) then
       used=send_data(id_dvdt, v_dt(is:ie,js:je,npz), time)
    endif
    if (id_dtdt > 0) then
       used=send_data(id_dtdt, t_dt(:,:,:), time)
    endif
    if (id_dqdt > 0) then
       used=send_data(id_dqdt, q_dt(:,:,:,sphum), time)
    endif

                        call timing_on('UPDATE_PHYS')
    call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                         u, v, w, delp, pt, q, qdiag, ua, va, ps, pe, peln, pk, pkz,  &
                         ak, bk, phis, u_srf, v_srf, ts,  &
                         delz, hydrostatic, u_dt, v_dt, t_dt, q_dt, &
                         moist_phys, Time, .false., gridstruct, &
                         gridstruct%agrid(:,:,1), gridstruct%agrid(:,:,2), &
                         npx, npy, npz, flagstruct, neststruct, bd, domain, ptop)
                        call timing_off('UPDATE_PHYS')
    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )

 end subroutine fv_phys


 subroutine sim_phys(npx, npy, npz, is, ie, js, je, ng, nq, nwat, pk, pkz,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va,   &
                     pt, delz, q, pe, delp, peln, sst, oro, hydrostatic, &
                     pdt, agrid, ak, bk, p_ref, fv_sg_adj,    &
                     Time, time_total, gridstruct)
!-----------------------
! A simple moist physics
!-----------------------


 integer, INTENT(IN) :: npx, npy, npz
 integer, INTENT(IN) :: is, ie, js, je, ng, nq, nwat
 integer, INTENT(IN) :: fv_sg_adj
 real   , INTENT(IN) :: pdt
 real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
 real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
 real, INTENT(IN):: p_ref
 real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
 logical, INTENT(IN):: hydrostatic

 type(time_type), intent(in) :: Time
 real, INTENT(IN), optional:: time_total

 real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
 real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
 real, INTENT(INOUT)::  pk (is:ie,js:je,npz+1)
 real, INTENT(INOUT)::  pkz(is:ie,js:je,npz)
 real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
 real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
 real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )
 real, INTENT(INOUT):: delz(is-ng :ie+ng  ,js-ng  :je+ng  ,npz)
 real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: sst(is:ie,js:je)

 type(fv_grid_type) :: gridstruct

! Tendencies:
 real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)
 real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)
! Local
 real, dimension(is:ie,js:je):: flux_t, flux_q, flux_u, flux_v, delm
 logical:: phys_hydrostatic = .true.
 real, parameter:: tice = 273.16
 real, parameter:: f1 = 2./3.
 real, parameter:: f2 = 1./3.
 real, dimension(is:ie,js:je,npz):: u3, v3, t3, p3, dz, zfull
 real, dimension(is:ie,js:je,npz+1):: zhalf, phalf
 real, dimension(is:ie,js:je,npz,nq):: q3
 real, dimension(is:ie,js:je):: rain, snow, ice, graup, land, mu
 real, dimension(is:ie,js:je):: ps, qs, rho, clouds
 real, dimension(is:ie,js:je):: olr, lwu, lwd, sw_surf
 real, dimension(is:ie,npz):: den
 real, dimension(is:ie,npz):: t_dt_rad
 real:: sday, rrg, tvm, olrm, swab, sstm, clds, deg2rad
 real:: rdt, tmp, cooling, heating,utmp, vtmp
 real:: fac_sm, rate_u, rate_v, rate_t, rate_q
 real:: prec
 integer  i,j,k, km, iq, k_mp
 integer  isd, ied, jsd, jed
 integer  seconds, days
 logical print_diag, used

!  if (.not. sim_phys_initialized) call fv_phys_init(nwat)

   call get_time (time, seconds, days)

   if ( mod(seconds, print_freq)==0 ) then
        print_diag = .true.
   else
        print_diag = .false.
   endif

   km = npz
   isd = is-ng;   ied = ie + ng
   jsd = js-ng;   jed = je + ng

   zvir = rvgas/rdgas - 1.
! Factor for Small-Earth Approx.
   fac_sm = radius / 6371.0e3
   rrg  = rdgas / grav
   sday  = 24.*3600.*fac_sm
   rdt = 1. / pdt

   deg2rad = pi/180.

     sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

   if ( nq.ge.6 ) then
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
   endif


   do k=1,km+1
      do j=js,je
         do i=is,ie
            phalf(i,j,k) = pe(i,k,j)
         enddo
      enddo
   enddo

! Compute zfull, zhalf
   do j=js,je
      do i=is,ie
         zhalf(i,j,npz+1) = 0.
         ps(i,j) = pe(i,km+1,j)
      enddo
   enddo

  if ( hydrostatic ) then
   do k=km,1,-1
      do j=js,je
         do i=is,ie
                     tvm = rrg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))
               dz(i,j,k) = -tvm*(peln(i,k+1,j)-peln(i,k,j))
               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            zfull(i,j,k) = zhalf(i,j,k+1) + tvm*(1.-pe(i,k,j)/p3(i,j,k))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
         enddo
      enddo
   enddo
 else
   do k=1,km
      do j=js,je
         do i=is,ie
               dz(i,j,k) = delz(i,j,k)
               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
            zfull(i,j,k) = zhalf(i,j,k+1) - 0.5*dz(i,j,k)
         enddo
      enddo
   enddo
 endif

   do k=1, km
      do j=js,je
         do i=is,ie
            u3(i,j,k) = ua(i,j,k) + pdt*u_dt(i,j,k)
            v3(i,j,k) = va(i,j,k) + pdt*v_dt(i,j,k)
            t3(i,j,k) = pt(i,j,k) + pdt*t_dt(i,j,k)
         enddo 
      enddo 
   enddo

   do j=js,je
      do i=is,ie
         delm(i,j) =  delp(i,j,km)/grav
          rho(i,j) = -delm(i,j)/dz(i,j,km)
      enddo
   enddo
   
!----------
! Setup SST:
!----------

! Need to save sst in a restart file if mixed_layer = .T.
 if ( .not. mixed_layer ) then
   if ( uniform_sst ) then
        sst(:,:) = sst0
   else
      do j=js,je
         do i=is,ie
! Neale & Hoskins 2001
            if ( abs(agrid(i,j,2)-shift_n*deg2rad) <  pi/3. ) then
                 tmp = sin( 1.5*(agrid(i,j,2)-shift_n*deg2rad) ) ** 2
! control:
!                sst(i,j) = 273.16 + 27.*(1.-tmp)
! QOBS:
                 sst(i,j) = 273.16 + 27.*(1.-0.5*(tmp + tmp**2))
            else
                 sst(i,j) = 273.16
            endif
         enddo
      enddo
   endif
 endif


   do iq=1,nq
      do k=1,npz
         do j=js,je
            do i=is,ie
               q3(i,j,k,iq) = q(i,j,k,iq) + pdt*q_dt(i,j,k,iq)
            enddo
         enddo
      enddo
   enddo


!----------------------------
! Apply net radiative cooling
!----------------------------
 if ( gray_rad ) then
    if ( prog_cloud ) then
         call get_low_clouds( is,ie, js,je, km, q3(is,js,1,liq_wat), q3(is,js,1,ice_wat),   &
                                                q3(is,js,1,cld_amt), clouds )
    else
         clouds(:,:) = low_c
    endif

    do j=js,je
       do k=1, km
          do i=is,ie
             den(i,k) = -delp(i,j,k)/(grav*dz(i,j,k))
          enddo 
       enddo 
       call gray_radiation(seconds, is, ie, km, agrid(is,j,1), agrid(is,j,2), clouds(is,j), &
                           sst(is,j), pt(is:ie,j,1:km), ps(is,j), phalf(is:ie,j,1:km+1), &
                           dz(is:ie,j,1:km), den, t_dt_rad,              &
                           olr(is,j), lwu(is,j), lwd(is,j), sw_surf(is,j), ng)
       do k=1, km
          do i=is,ie
             t_dt(i,j,k) = t_dt(i,j,k) + t_dt_rad(i,k)
          enddo 
       enddo 
    enddo
    if ( print_diag ) then
         olrm = g0_sum(olr, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
         swab = g0_sum(sw_surf, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)

         if ( mixed_layer ) then
              sstm = g0_sum(sst, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
              call prt_maxmin('TS', sst, is, ie, js, je, 0,  1, 1.0)
         endif

         if ( prog_cloud ) &
         clds = g0_sum(clouds, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)

         if( master ) then
            write(*,*) 'Domain mean OLR =', olrm
            write(*,*) 'Domain mean SWA =', swab
            if ( mixed_layer ) &
            write(*,*) 'Domain mean SST (mixed-layer)=', sstm
            if ( prog_cloud ) &
            write(*,*) 'Domain mean low clouds fraction =', clds
         endif
    endif

 else

! Prescribed (non-interating) heating/cooling rate:
   rate_t = 1. / (tau_strat*86400. + pdt )
!  cooling = cooling_rate/sday


   do k=1, km
      do j=js, je
         do i=is, ie
            cooling = cooling_rate*(f1+f2*cos(agrid(i,j,2)))/sday
            if ( p3(i,j,k) >= 100.E2 ) then
               t_dt(i,j,k) = t_dt(i,j,k) + cooling*min(1.0, (p3(i,j,k)-100.E2)/200.E2)
            elseif ( do_t_strat ) then
               if ( p3(i,j,k) < p_strat) then
                    t_dt(i,j,k) = t_dt(i,j,k) + (t_strat - t3(i,j,k))*rate_t
               endif
            endif
         enddo 
      enddo 
   enddo      ! k-loop

 endif


 if ( strat_rad ) then
! Solar heating above 100 mb:
      heating = heating_rate / 86400.
      do k=1, km
         do j=js,je
            do i=is,ie
               if ( p3(i,j,k) < 100.E2 ) then
                  t_dt(i,j,k) = t_dt(i,j,k) + heating*(100.E2-p3(i,j,k))/100.E2
               endif
            enddo 
         enddo 
      enddo 
 endif


 do k=1, km
    do j=js,je
       do i=is,ie
          t3(i,j,k) = pt(i,j,k) + pdt*t_dt(i,j,k)
       enddo 
    enddo 
 enddo



!---------------
! Surface fluxes
!---------------

if( do_mon_obkv ) then

  call qsmith(ie-is+1, je-js+1, 1, sst, ps, q3(is:ie,js:je,km,sphum), qs)
  call qsmith(ie-is+1, je-js+1, 1, sst, ps, qs, qs) ! Iterate once

! Need to save ustar in a restart file (sim_phys)
! Because u_star is prognostic but not saved
!  simple physics will not reproduce across restarts.
  if ( .not. allocated(u_star) ) then
       allocate ( u_star(is:ie,js:je) )
!      u_star(:,:) = 1.E25  ! large enough to cause blowup if used
       u_star(:,:) = 1.E-3
  endif

                                                                             call timing_on('mon_obkv')
  call mon_obkv(zvir, ps, t3(is:ie,js:je, km), zfull(is:ie,js:je,km),     &
                rho, p3(is:ie,js:je,km), u3(is:ie,js:je,km), v3(is:ie,js:je,km), sst, &
                qs, q3(is:ie,js:je,km,sphum), flux_t, flux_q, flux_u, flux_v, u_star, &
                delm, pdt, mu, t_fac, master)
                                                                             call timing_off('mon_obkv')
!---------------------------------------------------
! delp/grav = delm = kg/m**2
! watts = J/s = N*m/s = kg * m**2 / s**3
! CP = J/kg/deg
! flux_t = w/m**2 = J / (s*m**2)
!---------------------------------------------------
   do j=js,je
      do i=is,ie
               rate_u = flux_u(i,j)/delm(i,j)
           u3(i,j,km) = u3(i,j,km) + pdt*rate_u
         u_dt(i,j,km) = u_dt(i,j,km) + rate_u
               rate_v = flux_v(i,j)/delm(i,j)
           v3(i,j,km) =   v3(i,j,km) + pdt*rate_v
         v_dt(i,j,km) = v_dt(i,j,km) + rate_v
               rate_t = flux_t(i,j)/(cp_air*delm(i,j))
         t_dt(i,j,km) = t_dt(i,j,km) + rate_t
           t3(i,j,km) =   t3(i,j,km) + rate_t*pdt
               rate_q = flux_q(i,j)/delm(i,j)
         q_dt(i,j,km,sphum) = q_dt(i,j,km,sphum) + rate_q
           q3(i,j,km,sphum) =   q3(i,j,km,sphum) + rate_q*pdt
      enddo
   enddo
endif

if ( zero_winds ) then

! Zero out (doubly periodic) domain averaged winds with time-scale tau_zero:
  rate_u = min(1., 1./(tau_zero*86400.))
!$OMP parallel do default(none) shared(is,ie,js,je,km,u3,v3,gridstruct,u_dt,v_dt,rate_u,master) &
!$OMP                          private(utmp, vtmp)
 do k=1, km
    utmp = g0_sum(u3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
    vtmp = g0_sum(v3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
! The following might have issue for multiple thread
    if ( master ) then
        if( sqrt(utmp**2+vtmp**2) > 1. ) then 
            write(*,*) k, 'Domain avg winds=', utmp, vtmp
        endif
    endif
    do j=js,je
       do i=is,ie
          u_dt(i,j,k) = u_dt(i,j,k) - utmp * rate_u
          v_dt(i,j,k) = v_dt(i,j,k) - vtmp * rate_u
       enddo 
    enddo 
 enddo

endif

 if ( do_abl ) then
! ABL: vertical diffusion:
! do j=js,je
!    do i=is,ie
!       mu(i,j) = -sqrt(area(i,j)/area_ref)*dz(i,j,km)/tau_difz
!    enddo
! enddo

  call pbl_diff(hydrostatic, pdt, is, ie, js, je, ng, km, nq, u3, v3, t3,      &
                w, q3, delp, p3, pe, sst, mu, dz, u_dt, v_dt, t_dt, q_dt, gridstruct%area, print_diag )
 endif

  if ( mixed_layer .and. gray_rad ) then
     do j=js, je
        do i=is, ie
           sst(i,j) = sst(i,j)+pdt*(sw_surf(i,j)+lwd(i,j)-flux_t(i,j)-hlv*flux_q(i,j)-lwu(i,j))/c0
        enddo 
     enddo
  endif

! if ( do_sg_conv .and. fv_sg_adj > 0 ) then
!       call sg_conv(is, ie, js, je, isd, ied, jsd, jed, km, nq, pdt, fv_sg_adj,   &
!                    delp, phalf, p3, zfull, zhalf, t3, q3, u3, v3, w, &
!                    u_dt, v_dt, t_dt, q_dt, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, &
!                    hydrostatic, phys_hydrostatic)
! endif

  if ( do_lin_microphys ) then
      land(:,:) = 0.
      k_mp = 1
      do k=2, km
         tmp = ak(k) + bk(k)*1000.E2
         if ( tmp > 20.E2 ) then
              k_mp = k
              exit
         endif
      enddo
      if(master .and. print_diag) write(*,*) 'k_mp=', k_mp, tmp*0.01
      k_mp = 1

                                                                                          call timing_on('lin_cld_mp')
      call lin_cld_microphys_driver(q3(is:ie,js:je,1:km,sphum), q3(is:ie,js:je,1:km,liq_wat), q3(is:ie,js:je,1:km,rainwat),  &
                                    q3(is:ie,js:je,1:km,ice_wat), q3(is:ie,js:je,1:km,snowwat), q3(is:ie,js:je,1:km,graupel),  &
                                    q3(is:ie,js:je,1:km,cld_amt), q3(is:ie,js:je,1:km,cld_amt),   &
                                  q_dt(is:ie,js:je,1:km,sphum), q_dt(is:ie,js:je,1:km,liq_wat), q_dt(is:ie,js:je,1:km,rainwat), &
                                  q_dt(is:ie,js:je,1:km,ice_wat), q_dt(is:ie,js:je,1:km,snowwat), q_dt(is:ie,js:je,1:km,graupel), &
                                  q_dt(is:ie,js:je,1:km,cld_amt), t_dt, t3, w(is:ie,js:je,1:km), u3, v3,     &
                                  u_dt(is:ie,js:je,1:km),v_dt(is:ie,js:je,1:km), dz, zhalf,     &
                                  delp(is:ie,js:je,1:km), gridstruct%area(is:ie,js:je),  &
                                  pdt, land, rain, snow, ice, graup, hydrostatic, phys_hydrostatic, &
                                  1,ie-is+1, 1,je-js+1, 1,km, k_mp,npz, Time )
                                                                                          call timing_off('lin_cld_mp')
  elseif ( do_K_warm_rain ) then

! This is a stand alone routine (modified from NCAR code). It updates q, delp, and pt directly.
! Do not use other components in this module with this code (because NO tendencies)

      call K_warm_rain(pdt, is, ie, js, je, ng, km, zvir, sphum, liq_wat, rainwat, u3, v3, w, u_dt, v_dt,   &
                       q, pt, q_dt, t_dt, delp, delz, rain, Time)

      if ( id_rain_k>0 ) then
           prec_total(:,:) = prec_total(:,:) + rain(:,:)
           used = send_data(id_rain_k, prec_total, time)
           prec = g_sum(prec_total, is, ie, js, je, gridstruct%area(is:ie,js:je), 1)
           if(master) write(*,*) 'HIWPP accumulated rain=', prec
           call prt_maxmin('HIWPP: W', w, is, ie, js, je, ng,  km, 1.)
      endif

  endif

 end subroutine sim_phys


 subroutine pbl_diff(hydrostatic, dt, is, ie, js, je, ng, km, nq, ua, va, &
                     ta, w, q, delp, pm,  pe, ts, mu, dz, udt, vdt, tdt, qdt, area, print_diag )
 logical, intent(in):: hydrostatic
 integer, intent(in):: is, ie, js, je, ng, km, nq
 real, intent(in):: dt
 real, intent(in), dimension(is:ie,js:je):: ts, mu
 real, intent(in), dimension(is:ie,js:je,km):: dz, pm
 real, intent(in)::   pe(is-1:ie+1 ,1:km+1,js-1:je+1)
 real, intent(inout), dimension(is:ie,js:je,km):: ua, va, ta, tdt
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, udt, vdt, delp
 real, intent(inout), dimension(is:ie,js:je,km,nq):: q, qdt
 logical, intent(in):: print_diag
 real, intent(in) :: area(is-ng:ie+ng,js-ng:je+ng)
! Local:
 real, dimension(is:ie,js:je):: pblh
 real, dimension(is:ie,km+1):: gh
 real, dimension(is:ie,km):: nu, a, f, r, q2, den
 real, dimension(is:ie,km):: gz, hd, te
 real, dimension(is:ie):: kz
 real, parameter:: mu_min = 1.E-5    ! ~ molecular viscosity
 real, parameter:: ustar2 = 1.E-4
 integer:: n, i, j, k
 real:: cv, rcv, rdt, tmp, tvm, tv_surf, surf_h, rin

 cv = cp_air - rdgas
 rcv = 1./cv
 rdt = 1./dt

! Put opnmp loop here
 do 1000 j=js, je

    do i=is, ie
       gh(i,km+1) = 0.
       pblh(i,j) = 0.
    enddo
    do k=km, 1, -1
       do i=is, ie
          gh(i,k) = gh(i,k+1) - dz(i,j,k)
       enddo
    enddo

! Locate PBL top (m)
    do 125 i=is, ie
       tv_surf = ts(i,j)*(1.+zvir*q(i,j,km,sphum))
       do k=km, 3, -1 
          tvm = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)-q(i,j,k,liq_wat)-    &
                q(i,j,k,ice_wat)-q(i,j,k,snowwat)-q(i,j,k,rainwat)-q(i,j,k,graupel))
!         tvm = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)-q(i,j,k,liq_wat)-q(i,j,k,ice_wat))
          tvm = tvm*(pe(i,km+1,j)/pm(i,j,k))**kappa
          rin = grav*(gh(i,k+1)-0.5*dz(i,j,k))*(tvm-tv_surf) / ( 0.5*(tv_surf+tvm)*  &
                (ua(i,j,k)**2+va(i,j,k)**2+ustar2) )
          if ( rin > 1. ) then
               pblh(i,j) = gh(i,k+1)
               goto 125
          endif
       enddo
125 continue

  do k=km, 1, -1
     do i=is, ie
        if ( gh(i,k)>6.E3 .or. (pblh(i,j) < -0.5*dz(i,j,km)) ) then
             nu(i,k) = mu_min
        else
            kz(i) = 0.5*mu(i,j)*gh(i,km)
            surf_h = s_fac*pblh(i,j)
            if ( gh(i,k) <= surf_h) then
                 kz(i) = mu(i,j)*gh(i,k)
                 nu(i,k) = kz(i)
            elseif (gh(i,k) <= pblh(i,j)) then
! Use Dargan's form:
                 nu(i,k) = kz(i)*gh(i,k)/surf_h*(1.-(gh(i,k)-surf_h)/(pblh(i,j)-surf_h))**2
            else
                 nu(i,k) = mu_min
            endif
        endif
     enddo
  enddo

 do k=1,km-1
    do i=is, ie
       a(i,k) = -nu(i,k) / ( 0.5*(dz(i,j,k)+dz(i,j,k+1)) )
    enddo
 enddo
 do i=is, ie
    a(i,km) = 0.
 enddo

 if ( .not. hydrostatic ) then
    do k=1, km
       do i=is,ie
          gz(i,k) = grav*(gh(i,k+1) - 0.5*dz(i,j,k))
             tmp  = gz(i,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)
             tvm  = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)) 
          hd(i,k) = cp_air*tvm + tmp
          te(i,k) =     cv*tvm + tmp
       enddo
    enddo
 endif

! u-diffusion
    do k=1,km
       do i=is, ie
          den(i,k) = delp(i,j,k)/dz(i,j,k)
           q2(i,k) = ua(i,j,k)*den(i,k)
            f(i,k) = -dz(i,j,k)*rdt
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
           q2(i,k) = q2(i,k) / den(i,k)
          udt(i,j,k) = udt(i,j,k) + rdt*(q2(i,k) - ua(i,j,k))
           ua(i,j,k) = q2(i,k)
       enddo
    enddo
!---
! v-diffusion
    do k=1,km
       do i=is, ie
           q2(i,k) = va(i,j,k)*den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
             q2(i,k) = q2(i,k) / den(i,k)
          vdt(i,j,k) = vdt(i,j,k) + rdt*(q2(i,k) - va(i,j,k))
           va(i,j,k) = q2(i,k)
       enddo
    enddo
!---
 if ( .not. hydrostatic ) then
! w-diffusion
    do k=1,km
       do i=is, ie
           q2(i,k) = w(i,j,k) * den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          w(i,j,k) = q2(i,k) / den(i,k)
       enddo
    enddo
 endif
!--- micro-physics
 do n=1, nq
! if ( (n.ne.rainwat) .and. (n.ne.snowwat) .and. (n.ne.graupel) .and. (n.ne.cld_amt) ) then
  if ( n.ne.cld_amt ) then
    do k=1,km
       do i=is, ie
          q2(i,k) = q(i,j,k,n)*den(i,k)
           r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          q2(i,k) = q2(i,k) / den(i,k)
          qdt(i,j,k,n) = qdt(i,j,k,n) + rdt*(q2(i,k) - q(i,j,k,n))
            q(i,j,k,n) = q2(i,k)
       enddo
    enddo
  endif
 enddo

! Diffusion of dry static energy
 if ( .not. hydrostatic ) then
    do k=1,km
       do i=is, ie
           q2(i,k) = hd(i,k)*den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          te(i,k) = te(i,k) + q2(i,k)/den(i,k) - hd(i,k)
          q2(i,k) = rcv*(te(i,k)-gz(i,k)-0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)) &
                  / (1.+zvir*q(i,j,k,sphum))
          tdt(i,j,k) = tdt(i,j,k) + rdt*(q2(i,k) - ta(i,j,k))
           ta(i,j,k) = q2(i,k)
       enddo
    enddo
 endif

1000 continue

 if ( print_diag ) then
     tmp = g0_sum(pblh, is, ie, js, je, 0, area(is:ie,js:je), 1)
     if (master) write(*,*) 'Mean PBL H (km)=', tmp*0.001
     call prt_maxmin('PBLH(km)', pblh, is, ie, js, je, 0,  1, 0.001)
!    call prt_maxmin('K_ABL (m^2/s)', mu, is, ie, js, je, 0,  1, 1.)
 endif

 end subroutine pbl_diff

 subroutine trid_dif2(a, v, r, q, i1, i2, km)
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 real, intent(in), dimension(i1:i2,km):: a, v, r
 real, intent(out),dimension(i1:i2,km):: q
! Local:
 real:: gam(i1:i2,km)
 real:: bet(i1:i2)
 integer:: i, k

! Zero diffusive fluxes at top and bottom:
! top: k=1
  do i=i1,i2
     bet(i) = -(v(i,1) + a(i,1))
     q(i,1) = r(i,1) / bet(i)   
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) =  a(i,k-1) / bet(i)
          bet(i) = -( a(i,k-1)+v(i,k)+a(i,k) + a(i,k-1)*gam(i,k))
! a(i,km) = 0
          q(i,k) = ( r(i,k) - a(i,k-1)*q(i,k-1) ) / bet(i)
     enddo   
  enddo       
               
  do k=km-1,1,-1
     do i=i1,i2  
        q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
     enddo                                                                                              
  enddo                                                                                                 
                                                                                                        
 end subroutine trid_dif2


 subroutine fv_nudge( npz, is, ie, js, je, ng, u, v, w, delz, delp, pt, dt, hydrostatic)
!
! Nudge the prognostic varaibles toward the IC
! This is only useful for generating balanced steady state IC

    real   , INTENT(IN   ) :: dt
    integer, INTENT(IN   ) :: npz, is, ie, js, je, ng
    logical, INTENT(IN   ) :: hydrostatic
    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delz(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real c_w, c_p, c_t
    integer  i, j, k

    if ( .not. nudge_initialized ) then

         if( tau_winds > 0. ) then
             allocate ( u0(is:ie,  js:je+1,npz) )
             allocate ( v0(is:ie+1,js:je  ,npz) )
             do k=1,npz
                do j=js,je+1
                   do i=is,ie
                      u0(i,j,k) = u(i,j,k) 
                   enddo
                enddo
                do j=js,je
                   do i=is,ie+1
                      v0(i,j,k) = v(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_press > 0. ) then
             allocate ( dp(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      dp(i,j,k) = delp(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_temp > 0. ) then
             allocate ( t0(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      t0(i,j,k) = pt(i,j,k)
                   enddo
                enddo
             enddo
         endif

         nudge_initialized = .true.
         if ( is_master() ) write(*,*) 'Nudging of IC initialized.'
         return
    endif

! Nudge winds to initial condition:

      do k=1,npz
         if( tau_winds > 0. ) then
            c_w = dt/tau_winds
            do j=js,je+1
               do i=is,ie
                  u(i,j,k) = (u(i,j,k)+c_w*u0(i,j,k)) / (1.+c_w)
               enddo
            enddo
            do j=js,je
               do i=is,ie+1
                  v(i,j,k) = (v(i,j,k)+c_w*v0(i,j,k)) / (1.+c_w)
               enddo
            enddo
            if ( .not.hydrostatic ) then
            do j=js,je
               do i=is,ie+1
                  w(i,j,k) = w(i,j,k) / (1.+c_w)
               enddo
            enddo
            endif
         endif
         if ( tau_temp > 0. ) then
             c_t = dt/tau_temp
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = (pt(i,j,k)+c_t*t0(i,j,k)) / (1.+c_t)
                enddo
             enddo
             if ( .not.hydrostatic ) then
! delz
             endif
         endif

         if ( tau_press > 0. ) then
             c_p = dt/tau_press
             do j=js,je
                do i=is,ie
                   delp(i,j,k) = (delp(i,j,k)+c_p*dp(i,j,k)) / (1.+c_p)
                enddo
             enddo
         endif
      enddo

 end subroutine fv_nudge

 subroutine gray_radiation(sec, is, ie, km, lon, lat, clouds, ts, temp, ps, phalf, &
                           delz, rho, t_dt, olr, lwu, lwd, sw_surf, ng)

! Gray-Radiation algorithms based on Frierson, Held, and Zurita-Gotor, 2006 JAS
! Note: delz is negative
! Coded by S.-J. Lin, June 20, 2012
      integer, intent(in):: sec
      integer, intent(in):: is, ie, km, ng
      real, dimension(is:ie):: ts
      real, intent(in), dimension(is:ie):: lon, lat
      real, intent(in), dimension(is-ng:ie+ng,km):: temp, phalf, delz, rho
      real, intent(in), dimension(is:ie):: ps, clouds
      real, intent(out), dimension(is:ie,km):: t_dt
      real, intent(out), dimension(is:ie):: olr, lwu, lwd, sw_surf
! local:
      real, dimension(is:ie,km+1):: ur, dr
      real, dimension(is:ie,km+1):: tau, lw
      real, dimension(is:ie,km):: delt, b
      real, dimension(is:ie):: tau0
      real, parameter:: t0e =  8.    ! Dargan value= 6.
      real, parameter:: t0p = 1.5
      real, parameter::   fl = 0.1
      real, parameter:: sbc = 5.6734E-8
      real, parameter::  cp = 1004.64
      real:: rs0 = 938.4            ! W m-2
      real:: dels  = 1.4
      real:: sig, sw_rad, solar_ang
      integer:: i, k

      do i=is, ie
         tau0(i) = t0e + (t0p-t0e) * sin(lat(i))**2
      enddo
! Annual & global mean solar_abs ~ 0.25 * 1367 * ( 1-0.3) ~ 240
! Earth cross section/total_area = 0.25; 0.3 = Net cloud reflection and atm abs

      if ( diurnal_cycle ) then
           sw_rad = solar_constant*(1. - sw_abs) 
           do i=is, ie
               solar_ang = 2*pi*real(sec)/86400. + lon(i)
              sw_surf(i) = sw_rad*(1.-clouds(i))*cos(lat(i))*max(0.,cos(solar_ang))
           enddo
      else
           sw_rad = (1./pi) * solar_constant*(1. - sw_abs)
           do i=is, ie
              sw_surf(i) = sw_rad*(1.-clouds(i))*cos(lat(i))
           enddo
      endif

      do k=1, km+1
         do i=is, ie
#ifndef STRAT_OFF
! Dargan version:
                 sig = phalf(i,k)/ps(i) 
            tau(i,k) = tau0(i)*( sig*fl + (1.-fl)*sig**4 )
#else
! SJL: less cooling for the stratosphere
            tau(i,k) = tau0(i) * (phalf(i,k)/ps(i))**4
#endif
         enddo
      enddo
      do k=1, km
         do i=is, ie
            delt(i,k) = tau(i,k+1) - tau(i,k)
               b(i,k) =  sbc*temp(i,k)**4
         enddo
      enddo

! top down integration:
      do i=is, ie
         dr(i,1) = 0.
      enddo
      do k=1, km
         do i=is, ie
            dr(i,k+1) = (dr(i,k)+delt(i,k)*(b(i,k)-0.5*dr(i,k)))/(1.+0.5*delt(i,k)) 
         enddo
      enddo
! Bottom up
      do i=is, ie
         ur(i,km+1) = sbc*ts(i)**4
             lwu(i) = ur(i,km+1)
      enddo
      do k=km, 1, -1
         do i=is, ie
            ur(i,k) = (ur(i,k+1)+delt(i,k)*(b(i,k)-0.5*ur(i,k+1)))/(1.+0.5*delt(i,k)) 
         enddo
      enddo

! Compute net long wave cooling rate:
      do k=1, km+1
         do i=is, ie
            lw(i,k) = ur(i,k) - dr(i,k)
         enddo
      enddo
      do k=1, km
         do i=is, ie
            t_dt(i,k) = (lw(i,k) - lw(i,k+1))/(cp*rho(i,k)*delz(i,k))
         enddo
      enddo

      do i=is, ie
         olr(i) = ur(i,1)
         lwd(i) = dr(i,km+1)
      enddo

 end subroutine gray_radiation



 subroutine get_low_clouds( is,ie, js,je, km, ql, qi, qa, clouds )
 integer, intent(in):: is,ie, js,je, km
 real, intent(in), dimension(is:ie,js:je,km):: ql, qi, qa 
 real, intent(out), dimension(is:ie,js:je):: clouds
 integer:: i, j, k

 do j=js, je
    do i=is, ie
       clouds(i,j) = 0.
       do k=2,km
          if ((ql(i,j,k)>1.E-6 .or. qi(i,j,k)>1.e-5) .and.  qa(i,j,k)>1.E-4 ) then
! Maximum overlap
               clouds(i,j) = max(clouds(i,j), qa(i,j,k))
          endif
       enddo
       clouds(i,j) = min( 1., clouds(i,j) )
    enddo
 enddo

 end subroutine get_low_clouds

 subroutine fv_phys_init(is, ie, js, je, nwat, time, axes)
 integer, intent(IN) :: is, ie, js, je
 integer, intent(IN) :: nwat
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: time
 integer :: unit, ierr, io
 real:: total_area

!   ----- read and write namelist -----
    if ( file_exist('input.nml')) then
         unit = open_namelist_file ('input.nml')
         read  (unit, nml=sim_phys_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'sim_phys_nml')
 10     call close_file (unit)
    endif

    if (nwat /= 6 .and. do_lin_microphys) then
       call mpp_error(FATAL, "Need nwat == 6 to run Lin Microphysics.")
    endif

    id_rain_k = register_diag_field ('SIM_PHYS', 'rain_k', axes(1:2), time,        &
                'accumulated rain_Kessler', 'mm', missing_value=missing_value )
    id_vr_k = register_diag_field ('SIM_PHYS', 'vr_rain_k', axes(1:3), time,        &
                'Terminal fall V_Kessler', 'm/s', missing_value=missing_value )

    allocate ( prec_total(is:ie,js:je) )
    prec_total(:,:) = 0.

    master = is_master()

    if ( master ) then
        total_area = 4.*pi*radius**2
        write(*,*) 'Total surface area =', total_area
    endif
    if (do_reed_phys) then
       id_rain = register_diag_field ( 'dynamics', 'rain', axes(1:2), time,        &
            'rainfall_rate', 'mm H20/d', missing_value=missing_value )
    endif
    id_dqdt    = register_diag_field ( 'dynamics', 'dqdt', axes(1:3), time, &
         'humidity tendency', 'kg/kg/s', missing_value=missing_value )
    id_dTdt    = register_diag_field ( 'dynamics', 'dtdt', axes(1:3), time, &
         'temperature tendency', 'K/s', missing_value=missing_value )
    id_dudt    = register_diag_field ( 'dynamics', 'dudt', axes(1:3), time, &
         'zonal wind tendency', 'm/s/s', missing_value=missing_value )
    id_dvdt    = register_diag_field ( 'dynamics', 'dvdt', axes(1:3), time, &
         'meridional wind tendency', 'm/s/s', missing_value=missing_value )
    !To do: surface fluxes, surface winds, SST

    sim_phys_initialized = .true.

 end subroutine fv_phys_init


 real function g0_sum(p, ifirst, ilast, jfirst, jlast, ngc, area, mode)
 use mpp_mod,           only: mpp_sum
      real, save :: global_area

! Fast version of globalsum
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      integer :: i,j
      real gsum

!-------------------------
! Quick local sum algorithm
!-------------------------
      if ( .not. g0_sum_initialized ) then
         allocate (l_area(ifirst:ilast,jfirst:jlast))
         global_area = 0.
         do j=jfirst,jlast
           do i=ifirst,ilast
             global_area = global_area + area(i,j)
             l_area(i,j) = area(i,j)
           enddo
         enddo
         call mpp_sum(global_area)
!        if ( mpp_pe().eq.mpp_root_pe() ) write(*,*) 'Global Area=',global_area
         g0_sum_initialized = .true.
      end if

      gsum = 0.
      do j=jfirst,jlast
        do i=ifirst,ilast
          gsum = gsum + p(i,j)*l_area(i,j)
        enddo
      enddo
      call mpp_sum(gsum)

      if ( mode==1 ) then
        g0_sum = gsum / global_area
      else
        g0_sum = gsum
      endif

! Make it reproducible by truncating to 32-bit precision:
      g0_sum = real(g0_sum, 4)

 end function g0_sum

 subroutine K_warm_rain(dt, is, ie, js, je, ng, km, zvir, sphum, liq_wat, rainwat, u, v, w, u_dt, v_dt, &
                        q, pt, q_dt, t_dt, dp, delz, rain, Time)
 real, intent(in):: dt ! time step
 real, intent(in):: zvir
 type (time_type), intent(in) :: Time
 integer, intent(in):: is, ie, js, je, km, ng
 integer, intent(in):: sphum, liq_wat, rainwat
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: dp, delz, pt, w, u_dt, v_dt
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km,*):: q
 real, intent(inout),dimension(is:ie,js:je, km,*):: q_dt
 real, intent(inout),dimension(is:ie,js:je, km):: u,v, t_dt
 real, intent(out), dimension(is:ie,js:je):: rain
! Local:
 real, allocatable:: vr_k(:,:,:)
 real, dimension(km):: t1, q0, q1, q2, q3, zm, rho, drym, dm, dz
 real, dimension(km):: vr, qa, qb, qc, m1, u1, v1, w1, dgz, cvm, cvn
 real:: ze(km+1)
 real:: rdt
 integer i,j,k
 logical used

 if ( id_vr_k>0 ) then
      allocate ( vr_k(is:ie,js:je,km) )
 else
      allocate ( vr_k(is:is,js:js,1) )
 endif

 rdt = 1. / dt

 rain(:,:) = 0.

!$OMP parallel do default(none) shared(is,ie,js,je,km,delz,dp,sphum,liq_wat,rainwat,pt, &
!$OMP                                  u,v,w,dt,id_vr_k,vr_k,rain,do_k_momentum,u_dt,   &
!$OMP                                  rdt,v_dt,do_sedi_w)                              &
!$OMP                          private(dz,vr,dgz,cvm,cvn,m1,u1,v1,w1,q0,qa,qb,qc,t1,ze, &
!$OMP                                  dm,q1,q2,q3,zm,rho,drym)
 do j=js, je

    do i=is, ie
       ze(km+1) = 0.
       do k=km, 1, -1
          ze(k) = ze(k+1) - delz(i,j,k)
       enddo
       do k=1,km
! Moist air mass:
          dm(k) = dp(i,j,k)
! Tracer mass:
          qa(k) = q(i,j,k,sphum)
          qb(k) = q(i,j,k,liq_wat)
          qc(k) = q(i,j,k,rainwat)
          q1(k) = qa(k) * dm(k)
          q2(k) = qb(k) * dm(k)
          q3(k) = qc(k) * dm(k)
! Dry air mass per unit area
          drym(k) = dm(k) - (q1(k)+q2(k)+q3(k))
          dz(k) = -delz(i,j,k)
! Dry air density
          rho(k) = drym(k)/(grav*dz(k))
! Convert to dry mixing ratios:
          q1(k) = q1(k) / drym(k)
          q2(k) = q2(k) / drym(k)
          q3(k) = q3(k) / drym(k)
          t1(k) = pt(i,j,k)
! For condensate transport
          q0(k) = q1(k) + q2(k) + q3(k) ! total water
          u1(k) = u(i,j,k)
          v1(k) = v(i,j,k)
          w1(k) = w(i,j,k)
       enddo

!      do k=1,km
!         zm(k)  = 0.5*(ze(k)+ze(k+1))
!      enddo
!      call KESSLER1( t1, q1, q2, q3, vr, rho, dt, zm, km )
! Time implicit upwind scheme
       call kessler_imp( t1, q1, q2, q3, vr, rho, dt, dz, km )
       if ( id_vr_k>0 ) then
          do k=1, km
             vr_k(i,j,k) = vr(k)
          enddo
       endif

! Retrive rain_flux from non-local changes in total water
       m1(1) = drym(1)*(q0(1)-(q1(1)+q2(1)+q3(1)))
       do k=2,km
          m1(k) = m1(k-1) + drym(k)*(q0(k)-(q1(k)+q2(k)+q3(k)))
       enddo
       rain(i,j) = max(0., m1(km)) / grav

    if ( do_k_momentum ) then
! Momentum transport
       do k=2,km
! Here dm is total moist mass at previous time step
          u1(k) = (dm(k)*u1(k) + m1(k-1)*u1(k-1))/(dm(k)+m1(k-1))
          v1(k) = (dm(k)*v1(k) + m1(k-1)*v1(k-1))/(dm(k)+m1(k-1))
          u_dt(i,j,k) = (u1(k)-u(i,j,k))*rdt
          v_dt(i,j,k) = (v1(k)-v(i,j,k))*rdt
       enddo
       if ( do_sedi_w ) then
! Momentum conservation: (note vr is downward)
! dm is the total moist mass before terminal fall
          w(i,j,1) = (dm(1)*w1(1)+m1(1)*vr(1))/(dm(1)-m1(1))
          do k=2, km
             w(i,j,k) = (dm(k)*w1(k)-m1(k-1)*vr(k-1)+m1(k)*vr(k))/(dm(k)+m1(k-1)-m1(k))
          enddo
       endif
    endif

#ifdef CONSV_HEAT
! Heat transport
       do k=1, km
          dgz(k) = -0.5*grav*delz(i,j,k)    ! > 0
          cvn(k) = cv_air + q1(k)*cv_vap + (q2(k)+q3(k))*c_liq
       enddo
       k = 1
       cvm(k) = cvn(k) + c_liq*m1(k)/drym(k)
       t1(k) = (drym(k)*cvm(k)*t1(k) + m1(k)*dgz(k) ) / (drym(k)*cvn(k) + m1(k)*c_liq)
       do k=2, km
          cvm(k) = cvn(k) - c_liq*(m1(k-1)-m1(k))/drym(k)
          t1(k) = ( drym(k)*cvm(k)*t1(k) + m1(k-1)*c_liq*t1(k-1) + dgz(k)*(m1(k-1)+m1(k)) )    &
                / ( drym(k)*cvn(k) + m1(k)*c_liq )
       enddo
#endif

       do k=1,km
!         t_dt(i,j,k) = (t1(k)-pt(i,j,k))*rdt
!         delz(i,j,k) = delz(i,j,k)*t1(k)/pt(i,j,k)
          pt(i,j,k) = t1(k)
! Convert back to moist mixing ratios
          q1(k) = q1(k) * drym(k)
          q2(k) = q2(k) * drym(k)
          q3(k) = q3(k) * drym(k)
! Update total air mass:
          dp(i,j,k) = drym(k) + q1(k)+q2(k)+q3(k)
! Update tracers:
          q(i,j,k,  sphum) = q1(k) / dp(i,j,k)
          q(i,j,k,liq_wat) = q2(k) / dp(i,j,k)
          q(i,j,k,rainwat) = q3(k) / dp(i,j,k)
! Update temperature
       enddo
    enddo
 enddo

 if ( id_vr_k>0 ) then
      used = send_data(id_vr_k, vr_k, time)
      call prt_maxmin('VR_K', vr_k, is, ie, js, je, 0,  km, 1.)
 endif
 deallocate (vr_k)

 end subroutine K_warm_rain

 SUBROUTINE kessler_imp( T, QV, QC, QR, vr, rho, dt, dz, NZ )
! T  - TEMPERATURE (K)
! QV - WATER VAPOR MIXING RATIO (GM/GM)
! QC - CLOUD WATER MIXING RATIO (GM/GM)
! QR - RAIN WATER MIXING RATIO (GM/GM)
! R  - DRY AIR DENSITY (GM/M^3)
! dt - TIME STEP (S)
! Z - HEIGHTS OF THERMODYNAMIC LEVELS IN THE GRID COLUMN (M)
! NZ - NUMBER OF THERMODYNAMIC LEVELS IN THE COLUMN
! VARIABLES IN THE GRID COLUMN ARE ORDERED FROM THE SURFACE TO THE TOP.
! k=1 is the top layer
! Dry mixing ratios?
! OUTPUT VARIABLES:
   integer, intent(in):: nz
   real, intent(in):: dt
   REAL, intent(in)   :: rho(NZ), dz(nz)
   REAL, intent(inout):: T(NZ), QV(NZ), QC(NZ), QR(NZ)
   real, intent(out):: vr(nz)   ! terminal fall speed of rain * dt
! Local:
!  real, parameter:: c_k = 1003.    ! heat capacity of the K scheme
   real, parameter:: c_k = cp_air
   real, parameter:: F2X = 17.27
   real, parameter:: F5 = 237.3*F2X*2.5E6/c_k
   real, parameter:: XK = .2875
   real, parameter:: cvd = 717.56 
   REAL R(NZ), pc(NZ), qt(nz)
   REAL ERN, QRPROD, PROD, QVS, cvm, cfl
   real:: rdt
   INTEGER K

   rdt = 1. / dt
   DO K=1,NZ
       R(K) = 0.001*rho(k)
      pc(K) = 3.8e2 / (rho(k)*rdgas*T(k))
! vr is the terminal velocity of rain water:
      vr(k) = 36.34 * SQRT(rho(nz)/rho(k)) * (QR(K)*R(K))**0.1364
   END DO

! Sedimentation term using implicit upstream differencing
   qt(1) =  dz(1)*qr(1) / (dz(1)+dt*vr(1))
   do k=2, nz
      qt(k) = (dz(k)*qr(k) + qt(k-1)*r(k-1)*dt*vr(k-1)/r(k))  &
            / (dz(k)       +                dt*vr(k))
   enddo

   DO K=1,NZ
#ifdef CONSV_HEAT
      cvm = cvd + qv(k)*cv_vap + (qc(k)+qr(k))*c_liq
#endif
! Autoconversion and accretion rates following K&W78 Eq. 2.13a,b
! Threshold = 1.e-3
      QRPROD = QC(K) - (QC(K)-dt*max(.001*(QC(K)-.001), 0.))    &
             / (1.+dt*2.2*QR(K)**.875)
      QC(K) = max(QC(K) - QRPROD, 0.)
      QR(K) = max(qt(k) + QRPROD, 0.)

! Saturation vapor mixing ratio (gm/gm) following K&W78 Eq. 2.11
      QVS = pc(K)*EXP(F2X*(T(K)-273.) /(T(K)- 36.))

      PROD = (QV(K)-QVS)/(1.+QVS*F5/(T(K)-36.)**2)
! Evaporation rate following K&W78 Eq3. 3.8-3.10
      ERN = min(dt*(((1.6+124.9*(R(K)*QR(K))**.2046)   &
          *(R(K)*QR(K))**.525)/(2.55E6*pc(K)             &
          /(3.8 *QVS)+5.4E5))*(DIM(QVS,QV(K))            &
          /(R(K)*QVS)),max(-PROD-QC(K),0.),QR(K))
#ifdef CONSV_HEAT
! NICAM way:
        T(k) = T(k) + (Lv0+dc_vap*T(k))*(max(PROD,-QC(k))-ERN) / cvm
#else
! Saturation adjustment following K&W78 Eq.2.14a,b
!       T(K) = T(K) + 2.5E6*(max(PROD,-QC(K))-ERN) / c_k
        T(K) = T(K) + 2.5E6*(max(PROD,-QC(K))-ERN) / cv_air
#endif
       QV(K) = max(QV(K)-max(PROD,-QC(K))+ERN,0.)
       QC(K) = QC(K)+max(PROD,-QC(K))
       QR(K) = QR(K)-ERN
   END DO

 END SUBROUTINE kessler_imp


 SUBROUTINE KESSLER1( T, QV, QC, QR, vr, rho, dt, Z, NZ )
! T  - TEMPERATURE (K)
! QV - WATER VAPOR MIXING RATIO (GM/GM)
! QC - CLOUD WATER MIXING RATIO (GM/GM)
! QR - RAIN WATER MIXING RATIO (GM/GM)
! R  - DRY AIR DENSITY (GM/M^3)
! dt - TIME STEP (S)
! Z - HEIGHTS OF THERMODYNAMIC LEVELS IN THE GRID COLUMN (M)
! NZ - NUMBER OF THERMODYNAMIC LEVELS IN THE COLUMN
! VARIABLES IN THE GRID COLUMN ARE ORDERED FROM THE SURFACE TO THE TOP.
! k=1 is the top layer
! Dry mixing ratios?
! OUTPUT VARIABLES:
   integer, intent(in):: nz
   real, intent(in):: dt
   REAL, intent(in)   :: rho(NZ), Z(NZ)
   REAL, intent(inout):: T(NZ), QV(NZ), QC(NZ), QR(NZ)
   REAL, intent(out):: vr(nz)
! Local:
!  real, parameter:: c_k = 1003.    ! heat capacity of the K scheme
   real, parameter:: c_k = cp_air
   real, parameter:: F2X = 17.27
   real, parameter:: F5 = 237.3*F2X*2.5E6/cp_air
   real, parameter:: XK = .2875
   real, parameter:: cvd = 717.56 
   REAL R(NZ), SED(NZ), pc(NZ)
   REAL ERN, QRPROD, PROD, QVS, cvm, cfl
   INTEGER K

   DO K=1,NZ
       R(K) = 0.001*rho(k)
      pc(K) = 3.8e2 / (rho(k)*rdgas*T(k))
      vr(K) = 36.34*(QR(K)*R(K))**0.1364* SQRT(rho(nz)/rho(k))
   END DO

! Sedimentation term using upstream differencing
!
   SED(1) = -dt*QR(1)*vr(1)/(.5*(Z(1)-Z(2)))
   DO K=2,NZ
#ifdef CHECK_CFL
      cfl = dt*velqr(k-1)/(z(k-1)-z(k))
      if (cfl > 1. ) then
          write(*,*) 'CFL error: K_MP scheme', z(k), velqr(k), cfl
          call mpp_error(FATAL, "CFL > 1 from K_MP")
      endif
#endif
      SED(K) = dt*(R(K-1)*QR(K-1)*vr(K-1)         &
                  -R(K  )*QR(K  )*vr(K  )) / (R(K)*(Z(K-1)-Z(K)))
   END DO

   DO K=1,NZ
      cvm = cvd + qv(k)*cv_vap + (qc(k)+qr(k))*c_liq
! Autoconversion and accretion rates following K&W78 Eq. 2.13a,b
! Threshold = 1.e-3
      QRPROD = QC(K) - (QC(K)-dt*max(.001*(QC(K)-.001), 0.))    &
             / (1.+dt*2.2*QR(K)**.875)
      QC(K) = max(QC(K)-QRPROD,0.)
      QR(K) = max(QR(K)+QRPROD+SED(K),0.)

! Saturation vapor mixing ratio (gm/gm) following K&W78 Eq. 2.11
      QVS = pc(K)*EXP(F2X*(T(K)-273.) /(T(K)- 36.))

      PROD = (QV(K)-QVS)/(1.+QVS*F5/(T(K)-36.)**2)
! Evaporation rate following K&W78 Eq3. 3.8-3.10
      ERN = min(dt*(((1.6+124.9*(R(K)*QR(K))**.2046)   &
          *(R(K)*QR(K))**.525)/(2.55E6*pc(K)             &
          /(3.8 *QVS)+5.4E5))*(DIM(QVS,QV(K))            &
          /(R(K)*QVS)),max(-PROD-QC(K),0.),QR(K))
#ifdef CONSV_HEAT
        T(k) = T(k) + (Lv0+dc_vap*T(k))*(max(PROD,-QC(k))-ERN) / cvm
#else
! Saturation adjustment following K&W78 Eq.2.14a,b
!       T(K) = T(K) + 2.5E6*(max(PROD,-QC(K))-ERN) / c_k
        T(K) = T(K) + 2.5E6*(max(PROD,-QC(K))-ERN) / cv_air
#endif
       QV(K) = max(QV(K)-max(PROD,-QC(K))+ERN,0.)
       QC(K) = QC(K)+max(PROD,-QC(K))
       QR(K) = QR(K)-ERN
   END DO

 END SUBROUTINE KESSLER1


 real function g_sum(p, ifirst, ilast, jfirst, jlast, area, mode)
!-------------------------
! Quick local sum algorithm
!-------------------------
 use mpp_mod,           only: mpp_sum
 integer, intent(IN) :: ifirst, ilast
 integer, intent(IN) :: jfirst, jlast
 integer, intent(IN) :: mode  ! if ==1 divided by area
 real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
 real, intent(IN) :: area(ifirst:ilast,jfirst:jlast)
 integer :: i,j
 real gsum

   if( global_area < 0. ) then
       global_area = 0.
       do j=jfirst,jlast
          do i=ifirst,ilast
             global_area = global_area + area(i,j)
          enddo
       enddo
       call mpp_sum(global_area)
   end if

   gsum = 0.
   do j=jfirst,jlast
      do i=ifirst,ilast
         gsum = gsum + p(i,j)*area(i,j)
      enddo
   enddo
   call mpp_sum(gsum)

   if ( mode==1 ) then
        g_sum = gsum / global_area
   else
        g_sum = gsum
   endif

 end function g_sum

end module fv_phys_mod
