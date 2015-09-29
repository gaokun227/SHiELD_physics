module atmosphere_mod
#include <fms_platform.h>

!-----------------------------------------------------------------------
!
! Interface for Cubed_Sphere fv dynamical core
!
!-----------------------------------------------------------------------

!-----------------
! FMS modules:
!-----------------
use block_control_mod,      only: block_control_type
use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use time_manager_mod,       only: time_type, get_time, set_time, operator(+) 
use fms_mod,                only: file_exist, open_namelist_file,    &
                                  close_file, error_mesg, FATAL,     &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  set_domain,   &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default, nullify_domain
use mpp_mod,                only: mpp_error, stdout, FATAL, NOTE, &
                                  input_nml_file, mpp_root_pe,    &
                                  mpp_npes, mpp_pe, mpp_chksum,   &
                                  mpp_get_current_pelist,         &
                                  mpp_set_current_pelist
use mpp_domains_mod,        only: domain2d
use xgrid_mod,              only: grid_box_type
use field_manager_mod,      only: MODEL_ATMOS
use tracer_manager_mod,     only: get_tracer_index, get_number_tracers, &
                                  NO_TRACER
use gfs_physics_driver_mod, only: state_fields_in, state_fields_out, kind_phys

!-----------------
! FV core modules:
!-----------------
use fv_arrays_mod,      only: fv_atmos_type
use fv_control_mod,     only: fv_init, fv_end, ngrids
use fv_eta_mod,         only: get_eta_level
use fv_io_mod,          only: fv_io_register_nudge_restart
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time, prt_maxmin
use fv_restart_mod,     only: fv_restart, fv_write_restart
use fv_timing_mod,      only: timing_on, timing_off
use fv_mp_mod,          only: switch_current_Atm 
use fv_sg_mod,          only: fv_dry_conv
use fv_update_phys_mod, only: fv_update_phys
#if defined (ATMOS_NUDGE)
use atmos_nudge_mod,    only: atmos_nudge_init, atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
use fv_climate_nudge_mod,only: fv_climate_nudge_init,fv_climate_nudge_end
#elif defined (ADA_NUDGE)
use fv_ada_nudge_mod,   only: fv_ada_nudge_init, fv_ada_nudge_end
#else
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init
use amip_interp_mod,    only: forecast_mode
#endif

use mpp_domains_mod, only:  mpp_get_data_domain, mpp_get_compute_domain
use boundary_mod, only: update_coarse_grid

implicit none
private

!--- driver routines 
public :: atmosphere_init, atmosphere_end, atmosphere_restart, &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
public :: atmosphere_resolution, atmosphere_boundary, &
          atmosphere_grid_center, atmosphere_domain, &
          atmosphere_control_data, atmosphere_pref, &
          get_atmosphere_axes, get_bottom_mass, &
          get_bottom_wind, get_stock_pe, &
          set_atmosphere_pelist, get_atmosphere_grid

!--- physics/radiation data exchange routines
public :: atmos_phys_driver_statein

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
character(len=7)   :: mod_name = 'atmos'

!---- private data ----
  type (time_type) :: Time_step_atmos
  public Atm

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: nq                       ! transported tracers
  integer :: sec, seconds, days
  integer :: id_dynam, id_fv_diag, id_dryconv
  logical :: cold_start = .false.       ! read in initial condition

  integer, dimension(:), allocatable :: id_tracerdt_dyn
  integer :: num_tracers = 0

  integer :: mytile = 1
  integer :: p_split = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  integer :: id_udt_dyn, id_vdt_dyn

  real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

!---dynamics tendencies for use in fv_dry_conv and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
  real, allocatable, dimension(:,:,:,:) :: q_dt
  real, allocatable :: pref(:,:), dum1d(:)

!---need to define do_adiabatic_init to satisfy a reference when nwp_nudge is not the default
#if defined(ATMOS_NUDGE) || defined(CLIMATE_NUDGE) || defined(ADA_NUDGE)
   logical :: do_adiabatic_init
#endif

contains



 subroutine atmosphere_init (Time_init, Time, Time_step, Grid_box)
   type (time_type),      intent(in)    :: Time_init, Time, Time_step
   type(grid_box_type),   intent(inout) :: Grid_box

!--- local variables ---
   integer :: i, n
   integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units
   real :: ps1, ps2

                    call timing_on('ATMOS_INIT')
   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

   call get_number_tracers(MODEL_ATMOS, num_prog= num_tracers)

   zvir = rvgas/rdgas - 1.

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----
   !NOTE do we still need the second file_exist call?
   cold_start = (.not.file_exist('INPUT/fv_core.res.nc') .and. .not.file_exist('INPUT/fv_core.res.tile1.nc'))

   call fv_init( Atm, dt_atmos, grids_on_this_pe, p_split )  ! allocates Atm components

   do n=1,ngrids
      if (grids_on_this_pe(n)) mytile = n
   enddo

!----- write version and namelist to log file -----
   call write_version_number ( version, tagname )

!-----------------------------------

   npx   = Atm(mytile)%npx
   npy   = Atm(mytile)%npy
   npz   = Atm(mytile)%npz
   ncnst = Atm(mytile)%ncnst
   pnats = Atm(mytile)%flagstruct%pnats

   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec

   isd = isc - Atm(mytile)%bd%ng
   ied = iec + Atm(mytile)%bd%ng
   jsd = jsc - Atm(mytile)%bd%ng
   jed = jec + Atm(mytile)%bd%ng

   nq = ncnst-pnats

   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   ! This data is only needed for the COARSEST grid.
   call switch_current_Atm(Atm(mytile))

   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   allocate(Grid_box%area  (   isc:iec  , jsc:jec  ))
   allocate(Grid_box%edge_w(              jsc:jec+1))
   allocate(Grid_box%edge_e(              jsc:jec+1))
   allocate(Grid_box%edge_s(   isc:iec+1           ))
   allocate(Grid_box%edge_n(   isc:iec+1           ))
   allocate(Grid_box%en1   (3, isc:iec  , jsc:jec+1))
   allocate(Grid_box%en2   (3, isc:iec+1, jsc:jec  ))
   allocate(Grid_box%vlon  (3, isc:iec  , jsc:jec  ))
   allocate(Grid_box%vlat  (3, isc:iec  , jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%dy    (   isc:iec+1, jsc:jec  )
   Grid_box%area  (   isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%area  (   isc:iec  , jsc:jec  )
   Grid_box%edge_w(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_w(              jsc:jec+1)
   Grid_box%edge_e(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_e(              jsc:jec+1)
   Grid_box%edge_s(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_s(   isc:iec+1)
   Grid_box%edge_n(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_n(   isc:iec+1)
   Grid_box%en1   (:, isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%en1   (:, isc:iec  , jsc:jec+1)
   Grid_box%en2   (:, isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%en2   (:, isc:iec+1, jsc:jec  )
   do i = 1,3
     Grid_box%vlon  (i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlon  (isc:iec ,  jsc:jec, i )
     Grid_box%vlat  (i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlat  (isc:iec ,  jsc:jec, i )
   enddo

!----- allocate and zero out the dynamics (and accumulated) tendencies
   allocate( u_dt(isd:ied,jsd:jed,npz), &
             v_dt(isd:ied,jsd:jed,npz), &
             t_dt(isc:iec,jsc:jec,npz), &
             q_dt(isc:iec,jsc:jec,npz,nq) )
!--- allocate pref
    allocate(pref(npz+1,2), dum1d(npz+1))

   call set_domain ( Atm(mytile)%domain )
   call fv_restart(Atm(mytile)%domain, Atm, dt_atmos, seconds, days, cold_start, Atm(mytile)%gridstruct%grid_type, grids_on_this_pe)

   fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
   call fv_diag_init(Atm(mytile:mytile), Atm(mytile)%atmos_axes, Time, npx, npy, npz, Atm(mytile)%flagstruct%p_ref)

!---------- reference profile -----------
    ps1 = 101325.
    ps2 =  81060.
    pref(npz+1,1) = ps1
    pref(npz+1,2) = ps2
    call get_eta_level ( npz, ps1, pref(1,1), dum1d, Atm(mytile)%ak, Atm(mytile)%bk )
    call get_eta_level ( npz, ps2, pref(1,2), dum1d, Atm(mytile)%ak, Atm(mytile)%bk )

!--- initialize nudging module ---
#if defined (ATMOS_NUDGE)
    call atmos_nudge_init ( Time, Atm(mytile)%atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(mytile)%flagstruct%nudge ) then
         call mpp_error(NOTE, 'Code compiled with atmospheric nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge) then
         call mpp_error(NOTE, 'Code compiled with and using atmospheric nudging')
    endif
    Atm(mytile)%flagstruct%nudge = do_atmos_nudge
#elif defined (CLIMATE_NUDGE)
    call fv_climate_nudge_init ( Time, Atm(mytile)%atmos_axes(1:3), flag=do_atmos_nudge )
    if ( do_atmos_nudge .and. Atm(1)%flagstruct%nudge ) then
         call mpp_error(NOTE, 'Code compiled with climate nudging, but fv_core_nml nudge is also set to .true.')
    elseif ( do_atmos_nudge ) then
         call mpp_error(NOTE, 'Code compiled with and using climate nudging')
    endif
    Atm(mytile)%flagstruct%nudge = do_atmos_nudge
#elif defined (ADA_NUDGE)
    if ( Atm(1)%flagstruct%nudge ) then
        call fv_ada_nudge_init( Time, Atm(mytile)%atmos_axes, npz, zvir, Atm(1)%ak, Atm(1)%bk, Atm(1)%ts, &
           Atm(1)%phis, Atm(1)%gridstruct, Atm(1)%ks, Atm(1)%npx, Atm(1)%neststruct, Atm(1)%bd, Atm(1)%domain)
        call mpp_error(NOTE, 'ADA nudging is active')
     endif
#else
   !Only do nudging on coarse grid for now
   if ( Atm(mytile)%flagstruct%nudge ) then
      call fv_nwp_nudge_init( Time, Atm(mytile)%atmos_axes, npz, zvir, Atm(1)%ak, Atm(1)%bk, Atm(1)%ts, &
           Atm(1)%phis, Atm(1)%gridstruct, Atm(1)%ks, Atm(1)%npx, Atm(1)%neststruct, Atm(1)%bd)
        call mpp_error(NOTE, 'NWP nudging is active')
   endif
#endif

! This call needs to be separate from the register nudging restarts after initialization
   call fv_io_register_nudge_restart ( Atm )

!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dy-core', flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_dryconv   = mpp_clock_id ('FV dry conv',flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',    flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

#ifdef TWOWAY_UPDATE_BEFORE_PHYSICS
   if (ngrids > 1 ) call mpp_error(NOTE, "PERFORMING TWO-WAY UPDATING BEFORE PHYSICS")
#endif
                    call timing_off('ATMOS_INIT')

   n = mytile
   call switch_current_Atm(Atm(n)) 
      
 end subroutine atmosphere_init



 subroutine atmosphere_dynamics ( Time )
   type(time_type),intent(in) :: Time
   integer :: itrac, n, psc
   integer :: k, w_diff, nt_dyn

!---- Call FV dynamics -----

   call mpp_clock_begin (id_dynam)

   n = mytile
   do psc=1,p_split
                    call timing_on('fv_dynamics')
!uc/vc only need be same on coarse grid? However BCs do need to be the same
     call fv_dynamics(npx, npy, npz, nq, Atm(n)%ng, dt_atmos/real(p_split),&
                      Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill,  &
                      Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
                      Atm(n)%ptop, Atm(n)%ks, nq,                          &
                      Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split,&
                      Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,           &
                      Atm(n)%flagstruct%hydrostatic,                       & 
                      Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,         &
                      Atm(n)%pe, Atm(n)%pk, Atm(n)%peln,                   &
                      Atm(n)%pkz, Atm(n)%phis, Atm(n)%q_con,               &
                      Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc,        &
                      Atm(n)%vc, Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx,         &
                      Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0,        &
                      Atm(n)%flagstruct%hybrid_z,                          &
                      Atm(n)%gridstruct, Atm(n)%flagstruct,                &
                      Atm(n)%neststruct, Atm(n)%idiag, Atm(n)%bd,          &
                      Atm(n)%parent_grid, Atm(n)%domain)

     call timing_off('fv_dynamics')

    if (ngrids > 1 .and. psc < p_split) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split
    call mpp_clock_end (id_dynam)

#ifdef TWOWAY_UPDATE_BEFORE_PHYSICS
    if (ngrids > 1) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)
       call timing_off('TWOWAY_UPDATE')
    endif   
   call nullify_domain()
#endif

!-----------------------------------------------------
!--- COMPUTE DRY CONVECTION 
!-----------------------------------------------------
!--- zero out tendencies 
    call mpp_clock_begin (id_dryconv)
    u_dt(:,:,:)   = 0 
    v_dt(:,:,:)   = 0 
    t_dt(:,:,:)   = 0 
    q_dt(:,:,:,:) = 0 

    w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
    if ( Atm(n)%flagstruct%fv_sg_adj > 0 ) then
      nt_dyn = nq
      if ( w_diff /= NO_TRACER ) then
        nt_dyn = nq - 1
      endif
      call fv_dry_conv( isd, ied, jsd, jed, isc, iec, jsc, jec, Atm(n)%npz, &
                        nt_dyn, dt_atmos, Atm(n)%flagstruct%fv_sg_adj,      &
                        Atm(n)%flagstruct%nwat, Atm(n)%delp, Atm(n)%pe,     &
                        Atm(n)%peln, Atm(n)%pkz, Atm(n)%pt, Atm(n)%q,       &
                        Atm(n)%ua, Atm(n)%va, Atm(n)%flagstruct%hydrostatic,&
                        Atm(n)%w, Atm(n)%delz, u_dt, v_dt, t_dt, q_dt )
    endif

    if ( .not. Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default(shared) private(k)
       do k=1, Atm(n)%npz
          Atm(n)%q(isc:iec,jsc:jec,k,w_diff) = Atm(n)%w(isc:iec,jsc:jec,k) + w0_big
          q_dt(:,:,k,w_diff) = 0.
        enddo
    endif

   call mpp_clock_end (id_dryconv)

 end subroutine atmosphere_dynamics


 subroutine atmosphere_end (Time, Grid_box )!rab, Radiation, Physics)
   type (time_type),      intent(in)    :: Time
   type(grid_box_type),   intent(inout) :: Grid_box
!rab   type (radiation_type), intent(inout) :: Radiation
!rab   type (physics_type),   intent(inout) :: Physics

  ! initialize domains for writing global physics data
   call set_domain ( Atm(mytile)%domain )


!--- end nudging module ---
#if defined (ATMOS_NUDGE)
   if ( Atm(mytile)%flagstruct%nudge ) call atmos_nudge_end
#elif defined (CLIMATE_NUDGE)
   if ( Atm(mytile)%flagstruct%nudge ) call fv_climate_nudge_end
#elif defined (ADA_NUDGE)
   if ( Atm(mytile)%flagstruct%nudge ) call fv_ada_nudge_end
#else
   if ( Atm(mytile)%flagstruct%nudge ) call fv_nwp_nudge_end
#endif

   call nullify_domain ( )
   call fv_end(Atm, grids_on_this_pe)
   deallocate (Atm)

   deallocate( u_dt, v_dt, t_dt, q_dt, pref, dum1d )

 end subroutine atmosphere_end



  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call fv_write_restart(Atm, grids_on_this_pe, timestamp)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>


 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if

 end subroutine atmosphere_resolution


 subroutine atmosphere_pref (p_ref)
   real, dimension(:,:), intent(inout) :: p_ref

   p_ref = pref

 end subroutine atmosphere_pref


 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro
   i1 = Atm(mytile)%bd%isc
   i2 = Atm(mytile)%bd%iec
   j1 = Atm(mytile)%bd%jsc
   j2 = Atm(mytile)%bd%jec
   kt = Atm(mytile)%npz

   if (present(p_hydro)) p_hydro = Atm(mytile)%flagstruct%phys_hydrostatic
   if (present(  hydro))   hydro = Atm(mytile)%flagstruct%hydrostatic

 end subroutine atmosphere_control_data


 subroutine atmosphere_grid_center (lon, lat)
!---------------------------------------------------------------
!    returns the longitude and latitude cell centers
!---------------------------------------------------------------
    real,    intent(out) :: lon(:,:), lat(:,:)   ! Unit: radian
! Local data:
    integer i,j

    do j=jsc,jec
       do i=isc,iec
          lon(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%agrid(i,j,1)
          lat(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%agrid(i,j,2)
       enddo
    end do

 end subroutine atmosphere_grid_center


 subroutine atmosphere_boundary (blon, blat, global)
!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------
    real,    intent(out) :: blon(:,:), blat(:,:)   ! Unit: radian
    logical, intent(in), optional :: global
! Local data:
    integer i,j

    if( PRESENT(global) ) then
      if (global) call mpp_error(FATAL, '==> global grid is no longer available &
                               & in the Cubed Sphere')
    endif

    do j=jsc,jec+1
       do i=isc,iec+1
          blon(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%grid(i,j,1)
          blat(i-isc+1,j-jsc+1) = Atm(mytile)%gridstruct%grid(i,j,2)
       enddo
    end do

 end subroutine atmosphere_boundary


 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mytile)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


 subroutine atmosphere_domain ( fv_domain )
   type(domain2d), intent(out) :: fv_domain
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mytile)%domain_for_coupler

 end subroutine atmosphere_domain


 subroutine get_atmosphere_grid (dxmax, dxmin)
   real, intent(out) :: dxmax, dxmin

   dxmax = Atm(1)%gridstruct%da_max
   dxmin = Atm(1)%gridstruct%da_min

 end subroutine get_atmosphere_grid


 subroutine get_atmosphere_axes ( axes )
   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----
   if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                               'get_atmosphere_axes in atmosphere_mod', &
                               'size of argument is incorrect', FATAL   )

   axes (1:size(axes(:))) = Atm(mytile)%atmos_axes (1:size(axes(:)))
 
 end subroutine get_atmosphere_axes



 subroutine get_bottom_mass ( t_bot, tr_bot, p_bot, z_bot, p_surf, slp )
!--------------------------------------------------------------
! returns temp, sphum, pres, height at the lowest model level
! and surface pressure
!--------------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: t_bot, p_bot, z_bot, p_surf
   real, intent(out), optional, dimension(isc:iec,jsc:jec):: slp
   real, intent(out), dimension(isc:iec,jsc:jec,nq):: tr_bot
   integer :: i, j, m, k, kr
   real    :: rrg, sigtop, sigbot
   real, dimension(isc:iec,jsc:jec) :: tref
   real, parameter :: tlaps = 6.5e-3

   rrg  = rdgas / grav

   do j=jsc,jec
      do i=isc,iec
         p_surf(i,j) = Atm(mytile)%ps(i,j)
         t_bot(i,j) = Atm(mytile)%pt(i,j,npz)
         p_bot(i,j) = Atm(mytile)%delp(i,j,npz)/(Atm(mytile)%peln(i,npz+1,j)-Atm(mytile)%peln(i,npz,j))
         z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*Atm(mytile)%q(i,j,npz,1)) *  &
                      (1. - Atm(mytile)%pe(i,npz,j)/p_bot(i,j))
      enddo
   enddo

   if ( present(slp) ) then
     ! determine 0.8 sigma reference level
     sigtop = Atm(mytile)%ak(1)/pstd_mks+Atm(mytile)%bk(1)
     do k = 1, npz 
        sigbot = Atm(mytile)%ak(k+1)/pstd_mks+Atm(mytile)%bk(k+1)
        if (sigbot+sigtop > 1.6) then
           kr = k  
           exit    
        endif   
        sigtop = sigbot
     enddo
     do j=jsc,jec
        do i=isc,iec
           ! sea level pressure
           tref(i,j) = Atm(mytile)%pt(i,j,kr) * (Atm(mytile)%delp(i,j,kr)/ &
                            ((Atm(mytile)%peln(i,kr+1,j)-Atm(mytile)%peln(i,kr,j))*Atm(mytile)%ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = Atm(mytile)%ps(i,j)*(1.+tlaps*Atm(mytile)%phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
   endif

! Copy tracers
   do m=1,nq
      do j=jsc,jec
         do i=isc,iec
            tr_bot(i,j,m) = Atm(mytile)%q(i,j,npz,m)
         enddo
      enddo
   enddo

 end subroutine get_bottom_mass


 subroutine get_bottom_wind ( u_bot, v_bot )
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------
   real, intent(out), dimension(isc:iec,jsc:jec):: u_bot, v_bot
   integer i, j

   do j=jsc,jec
      do i=isc,iec
         u_bot(i,j) = Atm(mytile)%u_srf(i,j)
         v_bot(i,j) = Atm(mytile)%v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind



 subroutine get_stock_pe(index, value)
   integer, intent(in) :: index
   real,   intent(out) :: value

#ifdef USE_STOCK
   include 'stock.inc' 
#endif

   real wm(isc:iec,jsc:jec)
   integer i,j,k
   real, pointer :: area(:,:)
 
   area => Atm(mytile)%gridstruct%area

   select case (index)

#ifdef USE_STOCK
   case (ISTOCK_WATER)
#else
   case (1)
#endif
     
!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j=jsc,jec
        do k=1,npz
           do i=isc,iec
! Warning: the following works only with AM2 physics: water vapor; cloud water, cloud ice.
              wm(i,j) = wm(i,j) + Atm(mytile)%delp(i,j,k) * ( Atm(mytile)%q(i,j,k,1) +    &
                                                         Atm(mytile)%q(i,j,k,2) +    &
                                                         Atm(mytile)%q(i,j,k,3) )
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j=jsc,jec
        do i=isc,iec
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

   case default
     value = 0.0
   end select

 end subroutine get_stock_pe 


 subroutine atmosphere_state_update (Time, Statein, Stateout, Atm_block)
   type(time_type),                      intent(in) :: Time
   type(state_fields_in),  dimension(:), intent(in) :: Statein
   type(state_fields_out), dimension(:), intent(in) :: Stateout
   type(block_control_type),             intent(in) :: Atm_block
   type(time_type) :: Time_prev, Time_next
!--- local variables ---
   integer :: i, j, ix, k, n, w_diff, nt_dyn
   integer :: nb, ibs, ibe, jbs, jbe
   real ::  rcp

   Time_prev = Time
   Time_next = Time + Time_step_atmos

   n = mytile

   call set_domain ( Atm(mytile)%domain )

!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do shared (u_dt, v_dt, t_dt, q_dt, Atm, Statein, Stateout) &
!$OMP            private (nb, ibs, ibe, jbs, jbe, i, j, ix)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     ix = 0 
     do j=jbs,jbe
      do i=ibs,ibe
       ix = ix + 1
       do k = 1, npz
         u_dt(i,j,npz+1-k)   = (Stateout(nb)%gu0(ix,k) - Statein(nb)%ugrs(ix,k))/dt_atmos
         v_dt(i,j,npz+1-k)   = (Stateout(nb)%gv0(ix,k) - Statein(nb)%vgrs(ix,k))/dt_atmos
         t_dt(i,j,npz+1-k)   = (Stateout(nb)%gt0(ix,k) - Statein(nb)%tgrs(ix,k))/dt_atmos
         q_dt(i,j,npz+1-k,:) = (Stateout(nb)%gq0(ix,k,1:nq) - Statein(nb)%qgrs(ix,k,1:nq))/dt_atmos
       enddo
      enddo
     enddo

!--- diagnostic tracers are being updated in-place
!--- tracer fields must be returned to the Atm structure
       Atm(mytile)%qdiag(i,j,npz:1:-1,:) = Stateout(nb)%gq0(ix,1:npz,nq+1:ncnst)

   enddo

   w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )
   nt_dyn = ncnst-pnats   !nothing more than nq
   if ( w_diff /= NO_TRACER ) then
      nt_dyn = nt_dyn - 1
   endif

!--- adjust w and heat tendency for non-hydrostatic case
    if ( .not.Atm(n)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
      rcp = 1. / cp_air
!$OMP parallel do default(shared) private(i, j, k)
       do k=1, Atm(n)%npz
         do j=jsc, jec
           do i=isc, iec
             Atm(n)%q(i,j,k,w_diff) = q_dt(i,j,k,w_diff) ! w tendency due to phys
! Heating due to loss of KE (vertical diffusion of w)
             t_dt(i,j,k) = t_dt(i,j,k) - q_dt(i,j,k,w_diff)*rcp*&
                                     (Atm(n)%w(i,j,k)+0.5*dt_atmos*q_dt(i,j,k,w_diff))
             Atm(n)%w(i,j,k) = Atm(n)%w(i,j,k) + dt_atmos*Atm(n)%q(i,j,k,w_diff)
           enddo
         enddo
       enddo
    endif

       call timing_on('FV_UPDATE_PHYS')
    call fv_update_phys( dt_atmos, isc, iec, jsc, jec, isd, ied, jsd, jed, Atm(n)%ng, nt_dyn, &
                         Atm(n)%u,  Atm(n)%v,   Atm(n)%w,  Atm(n)%delp, Atm(n)%pt,         &
                         Atm(n)%q,  Atm(n)%qdiag,                                          &
                         Atm(n)%ua, Atm(n)%va,  Atm(n)%ps, Atm(n)%pe,   Atm(n)%peln,       &
                         Atm(n)%pk, Atm(n)%pkz, Atm(n)%ak, Atm(n)%bk,   Atm(n)%phis,       &
                         Atm(n)%u_srf, Atm(n)%v_srf, Atm(n)%ts, Atm(n)%delz,               &
                         Atm(n)%flagstruct%hydrostatic, u_dt, v_dt, t_dt, q_dt,            &
                         .true., Time_next, Atm(n)%flagstruct%nudge, Atm(n)%gridstruct,    &
                         Atm(n)%gridstruct%agrid(:,:,1), Atm(n)%gridstruct%agrid(:,:,2),   &
                         Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct,            &
                         Atm(n)%neststruct, Atm(n)%bd, Atm(n)%domain, Atm(n)%ptop)
       call timing_off('FV_UPDATE_PHYS')

!--- nesting update after updating atmospheric variables with
!--- physics tendencies
#ifndef TWOWAY_UPDATE_BEFORE_PHYSICS
       call timing_on('TWOWAY_UPDATE')
    if (ngrids > 1) then
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)
    endif   
       call timing_off('TWOWAY_UPDATE')
   call nullify_domain()
#endif

  !---- diagnostics for FV dynamics -----
   call mpp_clock_begin(id_fv_diag)

   fv_time = Time_next
   call get_time (fv_time, seconds,  days)

#if !defined(ATMOS_NUDGE) && !defined(CLIMATE_NUDGE) && !defined(ADA_NUDGE)
   if ( .not.forecast_mode .and. Atm(mytile)%flagstruct%nudge .and. Atm(mytile)%flagstruct%na_init>0 ) then
        if(mod(seconds, 21600)==0)  call adiabatic_init_drv (Time_prev, Time_next)
   endif
#endif

   call nullify_domain()
   call timing_on('FV_DIAG')
   call fv_diag(Atm(mytile:mytile), zvir, fv_time, Atm(mytile)%flagstruct%print_freq)
   call timing_off('FV_DIAG')

   call mpp_clock_end(id_fv_diag)

 end subroutine atmosphere_state_update

 subroutine adiabatic_init_drv (Time_prev, Time_next)
   type (time_type),       intent(in) :: Time_prev, Time_next
   integer:: n
   real, allocatable, dimension(:,:,:):: u_dt, v_dt, t_dt
   real, allocatable:: q_dt(:,:,:,:)
   integer:: isd, ied, jsd, jed, ngc

   character(len=80) :: errstr

!---------------------------------------------------
! Call the adiabatic forward-backward initialization
!---------------------------------------------------
   write(errstr,'(A, I4, A)') 'Performing adiabatic nudging',  Atm(mytile)%flagstruct%na_init, ' times'
   call mpp_error(NOTE, errstr)

        ngc = Atm(mytile)%ng
        isd = isc - ngc
        ied = iec + ngc
        jsd = jsc - ngc
        jed = jec + ngc

     call timing_on('adiabatic_init')

        allocate ( u_dt(isd:ied,jsd:jed, npz) )
        allocate ( v_dt(isd:ied,jsd:jed, npz) )
        allocate ( t_dt(isc:iec,jsc:jec, npz) )
        allocate ( q_dt(isc:iec,jsc:jec, npz, nq) )

     do_adiabatic_init = .true.

     do n=1,Atm(mytile)%flagstruct%na_init
        call adiabatic_init(Atm, Time_next, -dt_atmos, u_dt, v_dt, t_dt, q_dt, .false.)  ! Backward in time one step
        fv_time = Time_prev
        call adiabatic_init(Atm, Time_next,  dt_atmos, u_dt, v_dt, t_dt, q_dt, .true. )  ! Forward to the original time
        fv_time = Time_next
     enddo

     do_adiabatic_init = .false.

        deallocate ( u_dt )
        deallocate ( v_dt )
        deallocate ( t_dt )
        deallocate ( q_dt )

     call timing_off('adiabatic_init')

 end subroutine adiabatic_init_drv

 subroutine adiabatic_init (Atm, Time, dt_init, u_dt, v_dt, t_dt, q_dt, do_nudge)
!
! Perform adiabatic forward-backward adiabatic initialization with nudging
!
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type (time_type),    intent(in) :: Time
    real, intent(in):: dt_init
    logical, intent(in):: do_nudge
    real, intent(inout), dimension(:,:,:):: u_dt, v_dt, t_dt
    real, intent(inout)::  q_dt(:,:,:,:)
!
    integer:: isd, ied, jsd, jed, ngc, n
    type(time_type) :: Time_next

    Time_next = Time + Time_step_atmos

    n = mytile
    ngc = Atm(mytile)%ng
    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

    call fv_dynamics(npx, npy, npz, nq, ngc,  dt_init, 0.,                            &
                     .false.,  Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(n)%ptop, Atm(n)%ks, ncnst, Atm(n)%flagstruct%n_split,        &
                     Atm(n)%flagstruct%q_split, Atm(n)%u, Atm(n)%v, Atm(n)%w,         &
                     Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,                      &
                     Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,                     &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz, Atm(n)%phis,      &
                     Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc,  &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy,                    &
                     Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z,    &
                     Atm(n)%gridstruct, Atm(n)%flagstruct,                            &
                     Atm(n)%neststruct, Atm(n)%idiag, Atm(n)%bd,                      &
                     Atm(n)%parent_grid, Atm(n)%domain)

! No large-scale nudging at "Time_prev"
    if ( do_nudge ) then
         u_dt = 0.  ;   v_dt = 0. ;  t_dt = 0. ;  q_dt = 0.
    call fv_update_phys( abs(dt_init), isc, iec, jsc, jec, isd, ied, jsd, jed, ngc, nq,  &
                         Atm(n)%u,  Atm(n)%v,   Atm(n)%w,  Atm(n)%delp, Atm(n)%pt,       &
                         Atm(n)%q,  Atm(n)%qdiag,                                        &
                         Atm(n)%ua, Atm(n)%va,  Atm(n)%ps, Atm(n)%pe,   Atm(n)%peln,     &
                         Atm(n)%pk, Atm(n)%pkz, Atm(n)%ak, Atm(n)%bk,   Atm(n)%phis,     &
                         Atm(n)%u_srf, Atm(n)%v_srf, Atm(n)%ts, Atm(n)%delz,             &
                         Atm(n)%flagstruct%hydrostatic, u_dt, v_dt, t_dt, q_dt,          &
                         .true., Time_next, Atm(n)%flagstruct%nudge, Atm(n)%gridstruct,  &
                         Atm(n)%gridstruct%agrid(:,:,1), Atm(n)%gridstruct%agrid(:,:,2), &
                         Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct,          &
                         Atm(n)%neststruct, Atm(n)%bd, Atm(n)%domain, Atm(n)%ptop)

    endif

 end subroutine adiabatic_init


!rab subroutine atmos_phys_driver_statein (Time, Statein, Atm_block)
!rab   type (time_type),                     intent(in)    :: Time
 subroutine atmos_phys_driver_statein (Statein, Atm_block)
   type (state_fields_in), dimension(:), intent(inout) :: Statein
   type (block_control_type),            intent(in)    :: Atm_block
!--- local variables
   integer :: nb, npz, ibs, ibe, jbs, jbe, i, j, k, ix, sphum
   real(kind=kind_phys) :: pk0inv

   pk0inv = 1.0_kind_phys/(100000._kind_phys**kappa)

   sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )

   npz = Atm_block%npz
!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do shared (Atm, Statein, npz, nq, ncnst, sphum) &
!$OMP            private (nb, ibs, ibe, jbs, jbe, i, j, ix)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)
     
     ix = 0 
     do j=jbs,jbe
      do i=ibs,ibe
       ix = ix + 1
       !  surface pressure
       Statein(nb)%pgr(ix) = Atm(mytile)%ps(i,j)
       !-- level geopotential height
       !--- GFS physics assumes geopotential in a column is always zero at
       !--- the surface and only the difference between layers is important.
       Statein(nb)%phii(ix,1) = 0.0
!
       do k = 1, npz
         !  level pressure
         Statein(nb)%prsi(ix,k)  = Atm(mytile)%pe(i,npz+2-k,j)
         !  exner function pressure
         Statein(nb)%prsik(ix,k) = Atm(mytile)%pk (i,j,npz+2-k)*pk0inv   ! level
         Statein(nb)%prslk(ix,k) = Atm(mytile)%pkz(i,j,npz+1-k)*pk0inv   ! layer
         !  layer temp, u, & v
         Statein(nb)%tgrs(ix,k)  = Atm(mytile)%pt(i,j,npz+1-k)
         Statein(nb)%ugrs(ix,k)  = Atm(mytile)%ua(i,j,npz+1-k)
         Statein(nb)%vgrs(ix,k)  = Atm(mytile)%va(i,j,npz+1-k)
         !  layer vertical pressure velocity
         Statein(nb)%vvl(ix,k)   = Atm(mytile)%omga(i,j,npz+1-k)
!SJL
!SJL IF WE ARE GOING TO USE TRACER LIMITING, IT SHOULD OCCUR BELOW
!SJL
         !  layer tracers:  sphum, advected, non-advected
         Statein(nb)%qgrs_rad(ix,k)    = Atm(mytile)%q(i,j,npz+1-k,sphum)
         Statein(nb)%tracer(ix,k,1:nq) = Atm(mytile)%q(i,j,npz+1-k,:)
         Statein(nb)%tracer(ix,k,nq+1:ncnst) = Atm(mytile)%qdiag(i,j,npz+1-k,:)
!SJL
!SJL IF WE ARE GOING TO USE TRACER LIMITING, IT SHOULD OCCUR ABOVE
!SJL
         !-- level geopotential height
         Statein(nb)%phii(ix,k+1) = Statein(nb)%phii(ix,k) - Atm(mytile)%delz(i,j,npz+1-k)*grav
       enddo

       !  level pressure
       Statein(nb)%prsi(ix,npz+1) = Atm(mytile)%pe(i,1,j)
       !  level exner function pressure
       Statein(nb)%prsik(ix,npz+1) = Atm(mytile)%pk(i,j,1)*pk0inv
      enddo
     enddo

     do k = 1, npz
       !  layer pressure
       Statein(nb)%prsl(:,k) = 0.5_kind_phys*(Statein(nb)%prsi(:,k) + Statein(nb)%prsi(:,k+1))
       !  layer geopotential height
       if (Atm(mytile)%flagstruct%gfs_phil) then
         Statein(nb)%phil(:,k) = 0.0
       else
         Statein(nb)%phil(:,k) = 0.5_kind_phys*(Statein(nb)%phii(:,k) + Statein(nb)%phii(:,k+1))
       endif
     enddo
     !  reset this parameter to 1 just to be safe
     Statein(nb)%adjtrc = 1.0_kind_phys
     !  set the physics version of qgrs which is the same as tracer if separate arrays are needed
!rab Statein(nb)%qgrs(:,:,:) = Statein(nb)%tracer(:,:,;)
   enddo

 end subroutine atmos_phys_driver_statein

end module atmosphere_mod
