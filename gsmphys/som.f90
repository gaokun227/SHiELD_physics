!>  \file som.f
!!  This file contains routines for slab ocean model
!
!!  Written by Baoqiang Xiang (baoqiang.xiang@noaa.gov)

!  ==========================================================  !!!!!
!            'module_som' description            !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up SST using a slab ocean model (SOM)            !
!                                                                      !
!                                                                      !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'som_init'   -- initialization SOM  data           !
!                                                                      !
!      'update_som'     -- set up SOM          !
!         inputs:                                                      !
!           ( im, dtf, xcosz, Grid, islmsk, netflxsfc, qflux_restore, qflux_adj, &
!           tsclim, tsfc)
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!========================================!
      module module_som     !
!........................................!
!
      use physparam,         only : kind_phys
      use GFS_typedefs,      only : GFS_control_type, GFS_grid_type
!
      implicit   none
      private

!  ---  version tag and last revision date
!
      public  som_init, update_som

      character(len=24)     :: mld_option         = 'const'   ! use 'const'
                                                              !     'obs'
      integer               :: nudge_method       = 1         ! climatology
                                                  ! 2         ! climatology plus initial anomaly with a decay time scale
      real(kind=kind_phys)  :: const_mld          = 50.       ! meter
      real(kind=kind_phys)  :: sst_restore_tscale = 1.        ! day
      real(kind=kind_phys)  :: start_lat          = -60.      ! latitude starting from?
      real(kind=kind_phys)  :: end_lat            = 60.       ! latitude ending with?

      namelist /SOM_nml/   &
       mld_option, nudge_method, const_mld, sst_restore_tscale, start_lat, end_lat

! =================
      contains
! =================


!> This subroutine is the initialization program for surface radiation
!! related quantities (albedo, emissivity, etc.)
!!\param me       print control flag
!>\section gen_sfc_init General Algorithm
!! @{
!-----------------------------------
      subroutine som_init                                               &
     &     ( Model, logunit )!  ---  inputs:
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for surface radiation     !
!  related quantities (albedo, emissivity, etc.)                        !
!                                                                       !
! usage:         call sfc_init                                          !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!      me                                         !
!                                                                       !
!  outputs: (none) to module variables only                             !
!                                                                       !
!  external module variables:                                           !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs:
      type (GFS_control_type),    intent(in)  :: Model
      integer,    intent(in)  :: logunit
!      integer, intent(in) :: me

!  ---  outputs: ( none )

!  ---  locals:
      integer    :: i, k
!     integer    :: ia, ja
      logical    :: file_exist
      character  :: cline*80
      integer    :: ios
      logical    :: exists
!
!===> ...  begin here
!
        !--- read in the namelist
      inquire (file=trim(Model%fn_nml), exist=exists)
      if (.not. exists) then
      write(6,*) 'GFS_namelist_read:: namelist file: ',trim(Model%fn_nml),' does not exist'
      stop
      else
      open (unit=Model%nlunit, file=Model%fn_nml, READONLY, status='OLD', iostat=ios)
      endif
      rewind(Model%nlunit)
      read (Model%nlunit, nml=som_nml)
      close (Model%nlunit)
!--- write namelist to log file ---
      if (Model%me == Model%master) then
       write(logunit, *) "============================================="
       write(logunit, *) "Slab Ocean Model"
       write(logunit, nml=som_nml)
      endif

      return
!...................................
      end subroutine som_init
!-----------------------------------
      subroutine update_som                                           &
          ( im, dtf, xcosz, Grid, islmsk,  netflxsfc, qflux_restore,   &
           qflux_adj, tsclim, ts_clim_iano, tsfc) 

!  ===================================================================  !
!                                                                       !
!  this program computes the updated SST based on a simple SOM model    !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in)                                :: im
      real,    intent(in)                                :: dtf
      type (GFS_grid_type),  intent (in)                 :: Grid
      real (kind=kind_phys), dimension(:), intent(in)    ::         &
           xcosz, netflxsfc, tsclim, ts_clim_iano
      integer,  dimension(:), intent(in)                 :: islmsk

!  ---  inoutputs
      real (kind=kind_phys), dimension(:), intent(inout) ::         &
           qflux_restore, qflux_adj, tsfc

!  ---  locals:
      real (kind=kind_phys) :: lat, mlcp, tau, alpha
      integer :: i

!
!===> ...  begin here
!
          mlcp = const_mld * 4.e6
          tau = sst_restore_tscale*86400.
          alpha = 1. + dtf/tau
          do i = 1, im
           lat = Grid%xlat(i) * 57.29578
           if (islmsk(i) == 0 ) then
            if (lat >= start_lat .and. lat<= end_lat) then
!
             if (nudge_method == 1) then 
!              qflux_restore(i) = (tsclim(i) - tsfc(i)) * mlcp / tau
              tsfc(i) = (tsfc(i) + netflxsfc(i)/mlcp*dtf + tsclim(i)/tau*dtf ) / alpha 
             elseif (nudge_method == 2) then
!              qflux_restore(i) = (ts_clim_iano(i) - tsfc(i)) * mlcp / tau
              tsfc(i) = (tsfc(i) + netflxsfc(i)/mlcp*dtf + ts_clim_iano(i)/tau*dtf ) / alpha 
             endif
!             tsfc(i) = tsfc(i) + (qflux_restore(i)+netflxsfc(i))/mlcp*dtf  ! explicit
           else
            if (nudge_method == 1) then
             tsfc(i) = tsclim(i)
            elseif (nudge_method == 2) then
             tsfc(i) = ts_clim_iano(i)
            endif
           endif
           endif
          enddo

      return
!...................................
      end subroutine update_som
!-----------------------------------

      end module module_som 
!=========================================
