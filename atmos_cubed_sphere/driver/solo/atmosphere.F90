!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY, radius
use fms_mod,       only: file_exist, open_namelist_file,   &
                         error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         close_file, set_domain, nullify_domain, mpp_pe, mpp_root_pe, &
                         mpp_error, FATAL, NOTE
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d
use mpp_io_mod,       only: mpp_close
!------------------
! FV specific codes:
!------------------
use fv_arrays_mod, only: fv_atmos_type
use fv_control_mod,only: fv_init, fv_end, ngrids
use fv_phys_mod,   only: fv_phys, fv_nudge, fv_phys_init
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_timing_mod,   only: timing_on, timing_off
use fv_restart_mod, only: fv_restart
use fv_dynamics_mod, only: fv_dynamics
use fv_nesting_mod, only: twoway_nesting
use lin_cld_microphys_mod, only: lin_cld_microphys_init, lin_cld_microphys_end
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end
use fv_mp_mod, only: switch_current_Atm
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition
integer :: mytile = 1
integer :: p_split = 1

type(fv_atmos_type), allocatable, target :: Atm(:)

logical, allocatable :: grids_on_this_pe(:)
!-----------------------------------------------------------------------

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer :: axes(4)
    integer :: ss, ds
    integer i,j, isc, iec, jsc, jec
    real:: zvir
    integer :: f_unit, ios, ierr, n

                                           call timing_on('ATMOS_INIT')
  !----- write version and namelist to log file -----

    call write_version_number ( version, tagname )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

  !----- initialize FV dynamical core -----
    cold_start = (.not.file_exist('INPUT/fv_core.res.nc') .and. .not.file_exist('INPUT/fv_core.res.tile1.nc'))

    call fv_init(Atm, dt_atmos, grids_on_this_pe, p_split)  ! allocates Atm components

    do n=1,ngrids
       if (grids_on_this_pe(n)) mytile = n
    enddo

    isc = Atm(1)%bd%isc
    iec = Atm(1)%bd%iec
    jsc = Atm(1)%bd%jsc
    jec = Atm(1)%bd%jec

                   call timing_on('fv_restart')
    call fv_restart(Atm(1)%domain, Atm, dt_atmos, seconds, days, cold_start, &
         Atm(1)%flagstruct%grid_type, grids_on_this_pe)
                   call timing_off('fv_restart')

     fv_time = time

     do n=1,ngrids

       call switch_current_Atm(Atm(n))
       if ( grids_on_this_pe(n)) then

          Atm(N)%flagstruct%moist_phys = .false. ! need this for fv_diag calendar
          call fv_diag_init(Atm(n:n), axes, Time, Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%flagstruct%p_ref)

       endif

    if ( Atm(n)%flagstruct%adiabatic .or. Atm(n)%flagstruct%do_Held_Suarez ) then
         zvir = 0.         ! no virtual effect
         Atm(n)%flagstruct%moist_phys = .false.
    else
         zvir = rvgas/rdgas - 1.
         if ( grids_on_this_pe(n)) call fv_phys_init(isc,iec,jsc,jec,Atm(n)%flagstruct%nwat, Time, axes)
         Atm(n)%flagstruct%moist_phys = .true.
         if ( Atm(n)%flagstruct%nwat==6 .and. grids_on_this_pe(n))    then
            call lin_cld_microphys_init(iec-isc+1, jec-jsc+1, Atm(n)%npz, axes, Time)
         endif
    endif


    if ( Atm(n)%flagstruct%nudge .and. grids_on_this_pe(n) )    &
         call fv_nwp_nudge_init( Time, axes, Atm(n)%npz, zvir, Atm(n)%ak, Atm(n)%bk, Atm(n)%ts, &
         Atm(n)%phis, Atm(n)%gridstruct, Atm(n)%ks, Atm(n)%npx, Atm(n)%neststruct, Atm(n)%bd)

!   if( nlev > 1 ) call hs_forcing_init ( axes, Time )

    enddo

!-----------------------------------------------------------------------

                                           call timing_off('ATMOS_INIT')
  end subroutine atmosphere_init


!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    integer :: n, sphum, p, nc
    integer :: psc ! p_split counter


                                           call timing_on('ATMOSPHERE')
    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds
    
    do psc=1,p_split

    do n=1,ngrids

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif

       call switch_current_Atm(Atm(n)) 

       call set_domain(Atm(n)%domain)  ! needed for diagnostic output done in fv_dynamics

       if ( Atm(n)%flagstruct%nudge_ic )     &
            call  fv_nudge(Atm(n)%npz, Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, Atm(n)%ng, &
            Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%delp, Atm(n)%pt, dt_atmos/real(p_split), Atm(n)%flagstruct%hydrostatic )

       !---- call fv dynamics -----
       if ( Atm(n)%flagstruct%adiabatic .or. Atm(n)%flagstruct%do_Held_Suarez ) then
          zvir = 0.         ! no virtual effect
       else
          zvir = rvgas/rdgas - 1.
       endif

                                              call timing_on('fv_dynamics')
       call fv_dynamics(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%ncnst, Atm(n)%ng,   & 
            dt_atmos/real(p_split), Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill, &
            Atm(n)%flagstruct%reproduce_sum, kappa,   &
            cp_air, zvir, Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, &
            Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split, &
            Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,       &
            Atm(n)%flagstruct%hydrostatic, Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps, &
            Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz,             &
            Atm(n)%phis, Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc,  &
            Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy,    &
            Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z, Atm(n)%gridstruct, Atm(n)%flagstruct, &
            Atm(n)%neststruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid, Atm(n)%domain, time_total)
                                              call timing_off('fv_dynamics')
    end do

    if (ngrids > 1 .and. psc < p_split) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split

    do n=1,ngrids

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif

       if(Atm(n)%npz /=1 .and. .not. Atm(n)%flagstruct%adiabatic)then

                                                         call timing_on('FV_PHYS')
           call fv_phys(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%bd%isc, Atm(n)%bd%iec, &
                    Atm(n)%bd%jsc, Atm(n)%bd%jec, Atm(n)%ng, Atm(n)%ncnst,                &
                    Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%pt, Atm(n)%q, Atm(n)%pe,   &
                    Atm(n)%delp, Atm(n)%peln, Atm(n)%pkz, dt_atmos,                 &
                    Atm(n)%ua, Atm(n)%va, Atm(n)%phis, Atm(n)%gridstruct%agrid,     &
                    Atm(n)%ptop, Atm(n)%ak, Atm(n)%bk, Atm(n)%ks, Atm(n)%ps, Atm(n)%pk, &
                    Atm(n)%u_srf, Atm(n)%v_srf,  Atm(n)%ts, Atm(n)%delz,            &
                    Atm(n)%flagstruct%hydrostatic, Atm(n)%oro, .false.,             &
                    Atm(n)%flagstruct%p_ref,                                        &
                    Atm(n)%flagstruct%fv_sg_adj, Atm(n)%flagstruct%do_Held_Suarez,  &
                    Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,        &
                    Atm(n)%flagstruct%nwat, Atm(n)%bd,                              &
                    Atm(n)%domain, fv_time, time_total)
                                                        call timing_off('FV_PHYS')
       endif

       call nullify_domain()
    end do

    if (ngrids > 1) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, kappa, cp_air, zvir, dt_atmos)
       call timing_off('TWOWAY_UPDATE')
    endif


  !---- diagnostics for FV dynamics -----

    do n=1,ngrids

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif

       call nullify_domain()
       call timing_on('FV_DIAG')
       call fv_diag(Atm(n:n), zvir, fv_time, Atm(n)%flagstruct%print_freq)
       
       call timing_off('FV_DIAG')
    end do

                                           call timing_off('ATMOSPHERE')
 end subroutine atmosphere


 subroutine atmosphere_end

   integer n

    call get_time (fv_time, seconds,  days)

    do n=1,ngrids
       if ( Atm(n)%flagstruct%moist_phys .and. Atm(n)%flagstruct%nwat==6 .and. grids_on_this_pe(N)) call lin_cld_microphys_end
    enddo

    call fv_end(Atm, grids_on_this_pe)
    deallocate(Atm)

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
        
   fv_domain = Atm(mytile)%domain
        
 end subroutine atmosphere_domain

end module atmosphere_mod
