!>  \file som.f
!!  This file contains routines for a Slab Ocean Model (SOM)
!
!!  Written by Baoqiang Xiang (baoqiang.xiang@noaa.gov)

!  ==========================================================  !!!!!
!                          'module_som' description            !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up SST using a slab ocean model (SOM)            !
!                                                                      !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'som_init'   -- initialization SOM  by setting some namelists   !
!                                                                      !
!      'update_som' -- update SST with the combined effect of net      !
!                      surface heat flux and the nudging term          !
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
!
      public  som_init, update_som

      character(len=24)     :: mld_option         = 'obs'     ! option to set ocean mixed layer depth (MLD)
                                                              ! using either 'obs' or 'const'
      real(kind=kind_phys)  :: mld_obs_ratio      = 1.        ! tunning parameter for observed MLD
      integer               :: restore_method     = 1         ! option 1: nudging toward observational climatology
                                                  ! 2         ! option 2: nudging toward observational climatology plus
                                                              !           initial anomaly with a decay time scale of FTSFS (90 days)
      real(kind=kind_phys)  :: const_mld          = 40.       ! constant ocean MLD (meter)
      real(kind=kind_phys)  :: sst_restore_tscale = 3.        ! restoring time scale (day)
      real(kind=kind_phys)  :: start_lat          = -60.      ! latitude starting from?
      real(kind=kind_phys)  :: end_lat            = 60.       ! latitude ending with?
                                                              ! beyond the latitude bands (start_lat:end_lat), using climatological SST or
                                                              ! climatological SST plus initial anomaly 

      namelist /SOM_nml/   &
       mld_option, mld_obs_ratio, restore_method, const_mld,  &
       sst_restore_tscale, start_lat, end_lat

! =================
      contains
! =================


!-----------------------------------
      subroutine som_init                                               &
     &     ( Model, logunit )!  ---  inputs:

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for SOM model             !
!                                                                       !
! usage:         call som_init                                          !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!
      implicit none

!  ---  inputs:
      type (GFS_control_type),    intent(in)  :: Model
      integer,    intent(in)  :: logunit

!  ---  outputs: ( none )

!  ---  locals:
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
!
      subroutine update_som                                            &
          ( im, dtf, Grid, islmsk,  netflxsfc, qflux_restore,          &
           qflux_adj, mldclim, tsclim, ts_clim_iano, tsfc) 

!  ===================================================================  !
!                                                                       !
!  this program computes the updated SST based on a simple SOM model    !
!  A q-flux method is developing now and will be implemented later      !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in)                                :: im
      real,    intent(in)                                :: dtf  ! model time step
      type (GFS_grid_type),  intent (in)                 :: Grid
      real (kind=kind_phys), dimension(:), intent(in)    ::                &
                                                            netflxsfc,     & ! net surface heat flux
                                                            mldclim,       & ! ocean MLD
                                                            tsclim,        & ! observed climatological SST
                                                            ts_clim_iano     ! observed climatological SST plus initial anomaly
      integer,  dimension(:), intent(in)                 :: islmsk

!  ---  inoutputs
      real (kind=kind_phys), dimension(:), intent(inout) ::                &
                                                            qflux_restore, & ! restoring flux for diagnosis purpose
                                                            qflux_adj,     & ! qflux to be added later
                                                            tsfc             ! model SST

!  ---  locals:
      real (kind=kind_phys) :: lat, mlcp, tau, alpha
      integer :: i

!
!===> ...  begin here
!
          tau = sst_restore_tscale*86400.
          alpha = 1. + dtf/tau
          do i = 1, im
           if (mld_option == 'const') then
            mlcp = const_mld * 1000.*4.e3    ! rho*Cp*mld
           elseif (mld_option == 'obs') then
            mlcp =  mld_obs_ratio* mldclim(i) * 1000.*4.e3   ! rho*Cp*mld
           else
            write(*,*) ' mld_option can only be const or obs now'
            call abort
           endif
           
           lat = Grid%xlat(i) * 57.29578
           if (islmsk(i) == 0 ) then
            if (lat >= start_lat .and. lat<= end_lat) then
!
             if (restore_method == 1) then 
              qflux_restore(i) = (tsclim(i) - tsfc(i)) * mlcp / tau  ! for diagnosis purpose only
              tsfc(i) = (tsfc(i) + netflxsfc(i)/mlcp*dtf + tsclim(i)/tau*dtf ) / alpha 
             elseif (restore_method == 2) then
              qflux_restore(i) = (ts_clim_iano(i) - tsfc(i)) * mlcp / tau 
              tsfc(i) = (tsfc(i) + netflxsfc(i)/mlcp*dtf + ts_clim_iano(i)/tau*dtf ) / alpha 
             endif
!             tsfc(i) = tsfc(i) + (qflux_restore(i)+netflxsfc(i))/mlcp*dtf  ! explicit
           else
            if (restore_method == 1) then
             tsfc(i) = tsclim(i)
            elseif (restore_method == 2) then
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
