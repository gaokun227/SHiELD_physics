module FV3GFS_io_mod

!-----------------------------------------------------------------------
!    gfs_physics_driver_mod defines the GFS physics routines used by
!    the GFDL FMS system to obtain tendencies and boundary fluxes due 
!    to the physical parameterizations and processes that drive 
!    atmospheric time tendencies for use by other components, namely
!    the atmospheric dynamical core.
!
!    NOTE: This module currently supports only the operational GFS
!          parameterizations as of September 2015.  Further development
!          is needed to support the full suite of physical 
!          parameterizations present in the GFS physics package.
!-----------------------------------------------------------------------
!
!--- FMS/GFDL modules
  use block_control_mod,  only: block_control_type
  use mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, &
                                mpp_chksum, NOTE, FATAL, mpp_get_current_pelist_name
  use fms_mod,            only: file_exist, stdout
  use fms_io_mod,         only: restart_file_type, free_restart_type, &
                                register_restart_field,               &
                                restore_state, save_restart
  use mpp_domains_mod,    only: domain2d
  use time_manager_mod,   only: time_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use fv_mp_mod,          only: is_master, mp_reduce_sum, mp_reduce_min, mp_reduce_max
!
!--- GFS physics modules
  use machine,            only: kind_phys
!--- variables needed for calculating 'sncovr'
  use namelist_soilveg,   only: salp_data, snupx

!
! --- variables needed for Noah MP init
!
  use noahmp_tables,      only: laim_table,saim_table,sla_table,      &
                                bexp_table,smcmax_table,smcwlt_table, &
                                dwsat_table,dksat_table,psisat_table, &
                                isurban_table,isbarren_table,         &
                                isice_table,iswater_table

!
!--- GFS_typedefs
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_diag_type, GFS_grid_type
  use ozne_def,           only: oz_coeff
!
!--- IPD typdefs
  use IPD_typedefs,       only: IPD_control_type, IPD_data_type, &
                                IPD_restart_type
!--- GFS physics constants
  use physcons,           only: pi => con_pi, RADIUS => con_rerth, rd => con_rd
!--- needed for dq3dt output
  use ozne_def,           only: oz_coeff
!--- needed for cold-start capability to initialize q2m
  use gfdl_cld_mp_mod,    only: wqs1, qsmith_init
!
!-----------------------------------------------------------------------
  implicit none
  private
 
  !--- public interfaces ---
  public  FV3GFS_restart_read, FV3GFS_restart_write
  public  FV3GFS_IPD_checksum
  public  gfdl_diag_register, gfdl_diag_output

  !--- GFDL filenames
  character(len=32)  :: fn_oro = 'oro_data.nc'
  character(len=32)  :: fn_srf = 'sfc_data.nc'
  character(len=32)  :: fn_phy = 'phy_data.nc'

  !--- GFDL FMS netcdf restart data types
  type(restart_file_type) :: Oro_restart
  type(restart_file_type) :: Sfc_restart
  type(restart_file_type) :: Phy_restart
 
  !--- GFDL FMS restart containers
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2, sfc_name3
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2, phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3, phy_var3
  !--- Noah MP restart containers
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3sn,sfc_var3eq,sfc_var3zn

!-RAB
  type data_subtype
    real(kind=kind_phys), dimension(:),   pointer :: var2 => NULL()
    real(kind=kind_phys), dimension(:,:), pointer :: var3 => NULL()
    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
  end type data_subtype
  !--- data type definition for use with GFDL FMS diagnostic manager until write component is working
  type gfdl_diag_type
    private
    integer :: id
    integer :: axes
    logical :: time_avg
    character(len=64)    :: mod_name
    character(len=64)    :: name
    character(len=128)   :: desc
    character(len=64)    :: unit
    real(kind=kind_phys) :: cnvfac
    type(data_subtype), dimension(:), allocatable :: data
!rab    real(kind=kind_phys), dimension(:),   pointer :: var2 => NULL()
!rab    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
   end type gfdl_diag_type
   real(kind=kind_phys) :: zhour
!
   integer :: tot_diag_idx = 0
   integer, parameter :: DIAG_SIZE = 250
   real(kind=kind_phys), parameter :: missing_value = 1.d30
   type(gfdl_diag_type), dimension(DIAG_SIZE) :: Diag
!-RAB

 
!--- miscellaneous other variables
  logical :: module_is_initialized = .FALSE.

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!--------------------
! FV3GFS_restart_read
!--------------------
  subroutine FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain)
    type(IPD_data_type),      intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),   intent(inout) :: IPD_Restart
    type(block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),   intent(inout) :: Model
    type(domain2d),           intent(in)    :: fv_domain
 
    !--- read in surface data from chgres 
    call sfc_prop_restart_read (IPD_Data%Sfcprop, Atm_block, Model, fv_domain)
 
    !--- read in 
    if (Model%sfc_override) call sfc_prop_override  (IPD_Data%Sfcprop, IPD_Data%Grid, Atm_block, Model, fv_domain)

    !--- read in physics restart data
    call phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain)

  end subroutine FV3GFS_restart_read

!---------------------
! FV3GFS_restart_write
!---------------------
  subroutine FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    type(IPD_data_type),         intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),      intent(inout) :: IPD_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(IPD_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: fv_domain
    character(len=32), optional, intent(in)    :: timestamp
 
    !--- write surface data from chgres 
    call sfc_prop_restart_write (IPD_Data%Sfcprop, Atm_block, Model, fv_domain, timestamp)
 
    !--- write physics restart data
    call phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)

  end subroutine FV3GFS_restart_write


!--------------------
! FV3GFS_IPD_checksum
!--------------------
 subroutine FV3GFS_IPD_checksum (Model, IPD_Data, Atm_block)
   !--- interface variables
   type(IPD_control_type),    intent(in) :: Model
   type(IPD_data_type),       intent(in) :: IPD_Data(:)
   type (block_control_type), intent(in) :: Atm_block
   !--- local variables
   integer :: outunit, j, i, ix, nb, isc, iec, jsc, jec, lev, ct, l, ntr
   integer :: nsfcprop2d, idx_opt
   real(kind=kind_phys), allocatable :: temp2d(:,:,:)
   real(kind=kind_phys), allocatable :: temp3d(:,:,:,:)
   character(len=32) :: name

   isc = Model%isc
   iec = Model%isc+Model%nx-1
   jsc = Model%jsc
   jec = Model%jsc+Model%ny-1
   lev = Model%levs

   ntr = size(IPD_Data(1)%Statein%qgrs,3)

   if(Model%lsm == Model%lsm_noahmp) then
     nsfcprop2d = 149  
   else
     nsfcprop2d = 100
   endif

   allocate (temp2d(isc:iec,jsc:jec,nsfcprop2d+Model%ntot3d+Model%nctp))
   allocate (temp3d(isc:iec,jsc:jec,1:lev,17+Model%ntot3d+2*ntr))

   temp2d = 0.
   temp3d = 0.

   do j=jsc,jec
     do i=isc,iec
       nb = Atm_block%blkno(i,j) 
       ix = Atm_block%ixp(i,j) 
       !--- statein pressure
       temp2d(i,j, 1) = IPD_Data(nb)%Statein%pgr(ix)
       temp2d(i,j, 2) = IPD_Data(nb)%Sfcprop%slmsk(ix)
       temp2d(i,j, 3) = IPD_Data(nb)%Sfcprop%tsfc(ix)
       temp2d(i,j, 4) = IPD_Data(nb)%Sfcprop%tisfc(ix)
       temp2d(i,j, 5) = IPD_Data(nb)%Sfcprop%snowd(ix)
       temp2d(i,j, 6) = IPD_Data(nb)%Sfcprop%zorl(ix)
       temp2d(i,j, 7) = IPD_Data(nb)%Sfcprop%fice(ix)
       temp2d(i,j, 8) = IPD_Data(nb)%Sfcprop%hprime(ix,1)
       temp2d(i,j, 9) = IPD_Data(nb)%Sfcprop%sncovr(ix)
       temp2d(i,j,10) = IPD_Data(nb)%Sfcprop%snoalb(ix)
       temp2d(i,j,11) = IPD_Data(nb)%Sfcprop%alvsf(ix)
       temp2d(i,j,12) = IPD_Data(nb)%Sfcprop%alnsf(ix)
       temp2d(i,j,13) = IPD_Data(nb)%Sfcprop%alvwf(ix)
       temp2d(i,j,14) = IPD_Data(nb)%Sfcprop%alnwf(ix)
       temp2d(i,j,15) = IPD_Data(nb)%Sfcprop%facsf(ix)
       temp2d(i,j,16) = IPD_Data(nb)%Sfcprop%facwf(ix)
       temp2d(i,j,17) = IPD_Data(nb)%Sfcprop%slope(ix)
       temp2d(i,j,18) = IPD_Data(nb)%Sfcprop%shdmin(ix)
       temp2d(i,j,19) = IPD_Data(nb)%Sfcprop%shdmax(ix)
       temp2d(i,j,20) = IPD_Data(nb)%Sfcprop%tg3(ix)
       temp2d(i,j,21) = IPD_Data(nb)%Sfcprop%vfrac(ix)
       temp2d(i,j,22) = IPD_Data(nb)%Sfcprop%vtype(ix)
       temp2d(i,j,23) = IPD_Data(nb)%Sfcprop%stype(ix)
       temp2d(i,j,24) = IPD_Data(nb)%Sfcprop%uustar(ix)
       temp2d(i,j,25) = IPD_Data(nb)%Sfcprop%oro(ix)
       temp2d(i,j,26) = IPD_Data(nb)%Sfcprop%oro_uf(ix)
       temp2d(i,j,27) = IPD_Data(nb)%Sfcprop%hice(ix)
       temp2d(i,j,28) = IPD_Data(nb)%Sfcprop%weasd(ix)
       temp2d(i,j,29) = IPD_Data(nb)%Sfcprop%canopy(ix)
       temp2d(i,j,30) = IPD_Data(nb)%Sfcprop%ffmm(ix)
       temp2d(i,j,31) = IPD_Data(nb)%Sfcprop%ffhh(ix)
       temp2d(i,j,32) = IPD_Data(nb)%Sfcprop%f10m(ix)
       temp2d(i,j,33) = IPD_Data(nb)%Sfcprop%tprcp(ix)
       temp2d(i,j,34) = IPD_Data(nb)%Sfcprop%srflag(ix)
       temp2d(i,j,35) = IPD_Data(nb)%Sfcprop%slc(ix,1)
       temp2d(i,j,36) = IPD_Data(nb)%Sfcprop%slc(ix,2)
       temp2d(i,j,37) = IPD_Data(nb)%Sfcprop%slc(ix,3)
       temp2d(i,j,38) = IPD_Data(nb)%Sfcprop%slc(ix,4)
       temp2d(i,j,39) = IPD_Data(nb)%Sfcprop%smc(ix,1)
       temp2d(i,j,40) = IPD_Data(nb)%Sfcprop%smc(ix,2)
       temp2d(i,j,41) = IPD_Data(nb)%Sfcprop%smc(ix,3)
       temp2d(i,j,42) = IPD_Data(nb)%Sfcprop%smc(ix,4)
       temp2d(i,j,43) = IPD_Data(nb)%Sfcprop%stc(ix,1)
       temp2d(i,j,44) = IPD_Data(nb)%Sfcprop%stc(ix,2)
       temp2d(i,j,45) = IPD_Data(nb)%Sfcprop%stc(ix,3)
       temp2d(i,j,46) = IPD_Data(nb)%Sfcprop%stc(ix,4)
       temp2d(i,j,47) = IPD_Data(nb)%Sfcprop%t2m(ix)
       temp2d(i,j,48) = IPD_Data(nb)%Sfcprop%q2m(ix)
       temp2d(i,j,49) = IPD_Data(nb)%Coupling%nirbmdi(ix)
       temp2d(i,j,50) = IPD_Data(nb)%Coupling%nirdfdi(ix)
       temp2d(i,j,51) = IPD_Data(nb)%Coupling%visbmdi(ix)
       temp2d(i,j,52) = IPD_Data(nb)%Coupling%visdfdi(ix)
       temp2d(i,j,53) = IPD_Data(nb)%Coupling%nirbmui(ix)
       temp2d(i,j,54) = IPD_Data(nb)%Coupling%nirdfui(ix)
       temp2d(i,j,55) = IPD_Data(nb)%Coupling%visbmui(ix)
       temp2d(i,j,56) = IPD_Data(nb)%Coupling%visdfui(ix)
       temp2d(i,j,57) = IPD_Data(nb)%Coupling%sfcdsw(ix)
       temp2d(i,j,58) = IPD_Data(nb)%Coupling%sfcnsw(ix)
       temp2d(i,j,59) = IPD_Data(nb)%Coupling%sfcdlw(ix)
       temp2d(i,j,60) = IPD_Data(nb)%Grid%xlon(ix)
       temp2d(i,j,61) = IPD_Data(nb)%Grid%xlat(ix)
       temp2d(i,j,62) = IPD_Data(nb)%Grid%xlat_d(ix)
       temp2d(i,j,63) = IPD_Data(nb)%Grid%sinlat(ix)
       temp2d(i,j,64) = IPD_Data(nb)%Grid%coslat(ix)
       temp2d(i,j,65) = IPD_Data(nb)%Grid%area(ix)
       temp2d(i,j,66) = IPD_Data(nb)%Grid%dx(ix)
       if (Model%ntoz > 0) then
         temp2d(i,j,67) = IPD_Data(nb)%Grid%ddy_o3(ix)
       endif
       if (Model%h2o_phys) then
         temp2d(i,j,68) = IPD_Data(nb)%Grid%ddy_h(ix)
       endif
       temp2d(i,j,69) = IPD_Data(nb)%Cldprop%cv(ix)
       temp2d(i,j,70) = IPD_Data(nb)%Cldprop%cvt(ix)
       temp2d(i,j,71) = IPD_Data(nb)%Cldprop%cvb(ix)
       temp2d(i,j,72) = IPD_Data(nb)%Radtend%sfalb(ix)
       temp2d(i,j,73) = IPD_Data(nb)%Radtend%coszen(ix)
       temp2d(i,j,74) = IPD_Data(nb)%Radtend%tsflw(ix)
       temp2d(i,j,75) = IPD_Data(nb)%Radtend%semis(ix)
       temp2d(i,j,76) = IPD_Data(nb)%Radtend%coszdg(ix)
       temp2d(i,j,77) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfxc
       temp2d(i,j,78) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfx0
       temp2d(i,j,79) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfxc
       temp2d(i,j,80) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
       temp2d(i,j,81) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfxc
       temp2d(i,j,82) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfx0
       temp2d(i,j,83) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfxc
       temp2d(i,j,84) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfx0

        idx_opt = 85 
       if (Model%lsm == Model%lsm_noahmp) then
        temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%snowxy(ix)
        temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%tvxy(ix)
        temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%tgxy(ix)
        temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%canicexy(ix)
        temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%canliqxy(ix)
        temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%eahxy(ix)
        temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%tahxy(ix)
        temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%cmxy(ix)
        temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%chxy(ix)
        temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%fwetxy(ix)
        temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%sneqvoxy(ix)
        temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%alboldxy(ix)
        temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%qsnowxy(ix)
        temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%wslakexy(ix)
        temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%zwtxy(ix)
        temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%waxy(ix)
        temp2d(i,j,idx_opt+16) = IPD_Data(nb)%Sfcprop%wtxy(ix)
        temp2d(i,j,idx_opt+17) = IPD_Data(nb)%Sfcprop%lfmassxy(ix)
        temp2d(i,j,idx_opt+18) = IPD_Data(nb)%Sfcprop%rtmassxy(ix)
        temp2d(i,j,idx_opt+19) = IPD_Data(nb)%Sfcprop%stmassxy(ix)
        temp2d(i,j,idx_opt+20) = IPD_Data(nb)%Sfcprop%woodxy(ix)
        temp2d(i,j,idx_opt+21) = IPD_Data(nb)%Sfcprop%stblcpxy(ix)
        temp2d(i,j,idx_opt+22) = IPD_Data(nb)%Sfcprop%fastcpxy(ix)
        temp2d(i,j,idx_opt+23) = IPD_Data(nb)%Sfcprop%xsaixy(ix)
        temp2d(i,j,idx_opt+24) = IPD_Data(nb)%Sfcprop%xlaixy(ix)
        temp2d(i,j,idx_opt+25) = IPD_Data(nb)%Sfcprop%taussxy(ix)
        temp2d(i,j,idx_opt+26) = IPD_Data(nb)%Sfcprop%smcwtdxy(ix)
        temp2d(i,j,idx_opt+27) = IPD_Data(nb)%Sfcprop%deeprechxy(ix)
        temp2d(i,j,idx_opt+28) = IPD_Data(nb)%Sfcprop%rechxy(ix)

        temp2d(i,j,idx_opt+29) = IPD_Data(nb)%Sfcprop%snicexy(ix,-2)
        temp2d(i,j,idx_opt+30) = IPD_Data(nb)%Sfcprop%snicexy(ix,-1)
        temp2d(i,j,idx_opt+31) = IPD_Data(nb)%Sfcprop%snicexy(ix,0)
        temp2d(i,j,idx_opt+32) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-2)
        temp2d(i,j,idx_opt+33) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-1)
        temp2d(i,j,idx_opt+34) = IPD_Data(nb)%Sfcprop%snliqxy(ix,0)
        temp2d(i,j,idx_opt+35) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-2)
        temp2d(i,j,idx_opt+36) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-1)
        temp2d(i,j,idx_opt+37) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,0)
        temp2d(i,j,idx_opt+38) = IPD_Data(nb)%Sfcprop%smoiseq(ix,1)
        temp2d(i,j,idx_opt+39) = IPD_Data(nb)%Sfcprop%smoiseq(ix,2)
        temp2d(i,j,idx_opt+40) = IPD_Data(nb)%Sfcprop%smoiseq(ix,3)
        temp2d(i,j,idx_opt+41) = IPD_Data(nb)%Sfcprop%smoiseq(ix,4)
        temp2d(i,j,idx_opt+42) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-2)
        temp2d(i,j,idx_opt+43) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-1)
        temp2d(i,j,idx_opt+44) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,0)
        temp2d(i,j,idx_opt+45) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,1)
        temp2d(i,j,idx_opt+46) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,2)
        temp2d(i,j,idx_opt+47) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,3)
        temp2d(i,j,idx_opt+48) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,4)
        idx_opt = 134
       endif

       if (Model%nstf_name(1) > 0) then
         temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%tref(ix)
         temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%z_c(ix)
         temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%c_0(ix)
         temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%c_d(ix)
         temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%w_0(ix)
         temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%w_d(ix)
         temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%xt(ix)
         temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%xs(ix)
         temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%xu(ix)
         temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%xz(ix)
         temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%zm(ix)
         temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%xtts(ix)
         temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%xzts(ix)
         temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%ifd(ix)
         temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%dt_cool(ix)
         temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%qrain(ix)
       endif

       do l = 1,Model%ntot2d
         temp2d(i,j,nsfcprop2d+l) = IPD_Data(nb)%Tbd%phy_f2d(ix,l)
       enddo

       do l = 1,Model%nctp
         temp2d(i,j,nsfcprop2d+Model%ntot2d+l) = IPD_Data(nb)%Tbd%phy_fctd(ix,l)
       enddo

       temp3d(i,j,:, 1) = IPD_Data(nb)%Statein%phii(ix,1:lev)
       temp3d(i,j,:, 2) = IPD_Data(nb)%Statein%prsi(ix,1:lev)
       temp3d(i,j,:, 3) = IPD_Data(nb)%Statein%prsik(ix,1:lev)
       temp3d(i,j,:, 4) = IPD_Data(nb)%Statein%phil(ix,:)
       temp3d(i,j,:, 5) = IPD_Data(nb)%Statein%prsl(ix,:)
       temp3d(i,j,:, 6) = IPD_Data(nb)%Statein%prslk(ix,:)
       temp3d(i,j,:, 7) = IPD_Data(nb)%Statein%ugrs(ix,:)
       temp3d(i,j,:, 8) = IPD_Data(nb)%Statein%vgrs(ix,:)
       temp3d(i,j,:, 9) = IPD_Data(nb)%Statein%vvl(ix,:)
       temp3d(i,j,:,10) = IPD_Data(nb)%Statein%tgrs(ix,:)
       temp3d(i,j,:,11) = IPD_Data(nb)%Stateout%gu0(ix,:)
       temp3d(i,j,:,12) = IPD_Data(nb)%Stateout%gv0(ix,:)
       temp3d(i,j,:,13) = IPD_Data(nb)%Stateout%gt0(ix,:)
       temp3d(i,j,:,14) = IPD_Data(nb)%Radtend%htrsw(ix,:)
       temp3d(i,j,:,15) = IPD_Data(nb)%Radtend%htrlw(ix,:)
       temp3d(i,j,:,16) = IPD_Data(nb)%Radtend%swhc(ix,:)
       temp3d(i,j,:,17) = IPD_Data(nb)%Radtend%lwhc(ix,:)
       do l = 1,Model%ntot3d
         temp3d(i,j,:,17+l) = IPD_Data(nb)%Tbd%phy_f3d(ix,:,l)
       enddo
       do l = 1,ntr
         temp3d(i,j,:,17+Model%ntot3d+l)     = IPD_Data(nb)%Statein%qgrs(ix,:,l)
         temp3d(i,j,:,17+Model%ntot3d+ntr+l) = IPD_Data(nb)%Stateout%gq0(ix,:,l)
       enddo
     enddo
   enddo

   outunit = stdout()
   do i = 1,nsfcprop2d+Model%ntot2d+Model%nctp
     write (name, '(i3.3,3x,4a)') i, ' 2d '
     write(outunit,100) name, mpp_chksum(temp2d(:,:,i:i))
   enddo
   do i = 1,17+Model%ntot3d+2*ntr
     write (name, '(i2.2,3x,4a)') i, ' 3d '
     write(outunit,100) name, mpp_chksum(temp3d(:,:,:,i:i))
   enddo
100 format("CHECKSUM::",A32," = ",Z20)

   end subroutine FV3GFS_IPD_checksum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!----------------------------------------------------------------------      
! sfc_prop_restart_read
!----------------------------------------------------------------------      
!    creates and populates a data type which is then used to "register"
!    restart variables with the GFDL FMS restart subsystem.
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  register_restart_field, restart_state, free_restart
!   
!    opens:  oro_data.tile?.nc, sfc_data.tile?.nc
!   
!----------------------------------------------------------------------      
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain)
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(inout) :: Model
    type (domain2d),           intent(in)    :: fv_domain
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar_o2, nvar_s2m, nvar_s2o, nvar_s3
    integer :: nvar_s2mp, nvar_s3mp,isnow
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()

    character(len=64) :: fname
    !--- local variables for sncovr calculation
    integer :: vegtyp
    logical :: mand
    real(kind=kind_phys) :: rsnow, tem
    !--- Noah MP
    integer              :: soiltyp,ns,imon,iter,imn
    real(kind=kind_phys) :: masslai, masssai,snd
    real(kind=kind_phys) :: ddz,expon,aa,bb,smc,func,dfunc,dx
    real(kind=kind_phys) :: bexp, smcmax, smcwlt,dwsat,dksat,psisat
        
    real(kind=kind_phys), dimension(-2:0) :: dzsno
    real(kind=kind_phys), dimension(-2:4) :: dzsnso

    real(kind=kind_phys), dimension(4), save :: zsoil,dzs
    data dzs   /0.1,0.3,0.6,1.0/
    data zsoil /-0.1,-0.4,-1.0,-2.0/

    
    if (Model%cplflx) then ! needs more variables
      nvar_s2m = 34
    else
    nvar_s2m = 32
    endif
    nvar_o2  = 19
    nvar_s2o = 18
    nvar_s3  = 3

    if (Model%lsm == Model%lsm_noahmp) then
      nvar_s2mp = 34       !mp 2D
      nvar_s3mp = 5        !mp 3D
    else
      nvar_s2mp = 0        !mp 2D
      nvar_s3mp = 0        !mp 3D
    endif

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)
 
    !--- OROGRAPHY FILE
    if (.not. allocated(oro_name2)) then
    !--- allocate the various containers needed for orography data
      allocate(oro_name2(nvar_o2))
      allocate(oro_var2(nx,ny,nvar_o2))
      oro_var2 = -9999._kind_phys

      oro_name2(1)  = 'stddev'     ! hprime(ix,1)
      oro_name2(2)  = 'convexity'  ! hprime(ix,2)
      oro_name2(3)  = 'oa1'        ! hprime(ix,3)
      oro_name2(4)  = 'oa2'        ! hprime(ix,4)
      oro_name2(5)  = 'oa3'        ! hprime(ix,5)
      oro_name2(6)  = 'oa4'        ! hprime(ix,6)
      oro_name2(7)  = 'ol1'        ! hprime(ix,7)
      oro_name2(8)  = 'ol2'        ! hprime(ix,8)
      oro_name2(9)  = 'ol3'        ! hprime(ix,9)
      oro_name2(10) = 'ol4'        ! hprime(ix,10)
      oro_name2(11) = 'theta'      ! hprime(ix,11)
      oro_name2(12) = 'gamma'      ! hprime(ix,12)
      oro_name2(13) = 'sigma'      ! hprime(ix,13)
      oro_name2(14) = 'elvmax'     ! hprime(ix,14)
      oro_name2(15) = 'orog_filt'  ! oro
      oro_name2(16) = 'orog_raw'   ! oro_uf
      oro_name2(17) = 'land_frac'  ! land fraction [0:1]
      !--- variables below here are optional
      oro_name2(18) = 'lake_frac'  ! lake fraction [0:1]
      oro_name2(19) = 'lake_depth' ! lake depth(m)
      !--- register the 2D fields
      do num = 1,nvar_o2
        var2_p => oro_var2(:,:,num)
        if (trim(oro_name2(num)) == 'lake_frac' .or. trim(oro_name2(num)) == 'lake_depth') then
          id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
        id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=fv_domain)
        endif
      enddo
      nullify(var2_p)
    endif

    fname = 'INPUT/'//trim(fn_oro)
    if (file_exist(fname)) then 
       
       !--- read the orography restart/data
       call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
       call restore_state(Oro_restart)

       Model%frac_grid = .false.
       !--- copy data into GFS containers
       do nb = 1, Atm_block%nblks
          !--- 2D variables
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- stddev
             Sfcprop(nb)%hprim(ix)     = oro_var2(i,j,1)
             !--- hprime(1:14)
             Sfcprop(nb)%hprime(ix,1)  = oro_var2(i,j,1)
             Sfcprop(nb)%hprime(ix,2)  = oro_var2(i,j,2)
             Sfcprop(nb)%hprime(ix,3)  = oro_var2(i,j,3)
             Sfcprop(nb)%hprime(ix,4)  = oro_var2(i,j,4)
             Sfcprop(nb)%hprime(ix,5)  = oro_var2(i,j,5)
             Sfcprop(nb)%hprime(ix,6)  = oro_var2(i,j,6)
             Sfcprop(nb)%hprime(ix,7)  = oro_var2(i,j,7)
             Sfcprop(nb)%hprime(ix,8)  = oro_var2(i,j,8)
             Sfcprop(nb)%hprime(ix,9)  = oro_var2(i,j,9)
             Sfcprop(nb)%hprime(ix,10) = oro_var2(i,j,10)
             Sfcprop(nb)%hprime(ix,11) = oro_var2(i,j,11)
             Sfcprop(nb)%hprime(ix,12) = oro_var2(i,j,12)
             Sfcprop(nb)%hprime(ix,13) = oro_var2(i,j,13)
             Sfcprop(nb)%hprime(ix,14) = oro_var2(i,j,14)
                  !--- oro
             Sfcprop(nb)%oro(ix)       = oro_var2(i,j,15)
                  !--- oro_uf
             Sfcprop(nb)%oro_uf(ix)    = oro_var2(i,j,16)
             Sfcprop(nb)%landfrac(ix)  = oro_var2(i,j,17) !land frac [0:1]
             Sfcprop(nb)%lakefrac(ix)  = oro_var2(i,j,18) !lake frac [0:1]
          enddo
       enddo

       if (nint(oro_var2(1,1,18)) == -9999._kind_phys) then ! lakefrac doesn't exist in the restart, need to create it
         if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - will computing lakefrac') 
         Model%frac_grid = .false.
       else
         Model%frac_grid = .true.
       endif
       
       if (Model%me == Model%master ) write(0,*)' resetting Model%frac_grid=',Model%frac_grid
       
       !--- deallocate containers and free restart container
       deallocate(oro_name2, oro_var2)
       call free_restart_type(Oro_restart)     
       

    else ! cold_start (no way yet to create orography on-the-fly)

       call mpp_error(NOTE,'No INPUT/oro_data.tile*.nc orographic data found; setting to 0')
       !--- copy data into GFS containers
       do nb = 1, Atm_block%nblks
          !--- 2D variables
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- stddev
             Sfcprop(nb)%hprim(ix)      = 0.0
             !--- hprime(1:14)
             Sfcprop(nb)%hprime(ix,1:14)  = 0.0
             !--- oro
             Sfcprop(nb)%oro(ix)        = 0.0
             !--- oro_uf
             Sfcprop(nb)%oro_uf(ix)     = 0.0
          enddo
       enddo

    endif
 
 
    !--- SURFACE FILE
    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_name3(nvar_s3+nvar_s3mp))

      allocate(sfc_var2(nx,ny,nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar_s3))
      sfc_var2 = -9999._kind_phys
      sfc_var3 = -9999._kind_phys
!
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))
        sfc_var3sn = -9999._kind_phys
        sfc_var3eq = -9999._kind_phys
        sfc_var3zn = -9999._kind_phys
      end if
 
      !--- names of the 2D variables to save
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsea'    !tsfc
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorl'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
      !--- variables below here are optional
      sfc_name2(32) = 'sncovr'
      if(Model%cplflx) then
        sfc_name2(33) = 'tsfcl'   !temp on land portion of a cell
        sfc_name2(34) = 'zorll'   !zorl on land portion of a cell
      end if

      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0) 
      sfc_name2(nvar_s2m+1)  = 'tref'
      sfc_name2(nvar_s2m+2)  = 'z_c'
      sfc_name2(nvar_s2m+3)  = 'c_0'
      sfc_name2(nvar_s2m+4)  = 'c_d'
      sfc_name2(nvar_s2m+5)  = 'w_0'
      sfc_name2(nvar_s2m+6)  = 'w_d'
      sfc_name2(nvar_s2m+7)  = 'xt'
      sfc_name2(nvar_s2m+8)  = 'xs'
      sfc_name2(nvar_s2m+9)  = 'xu'
      sfc_name2(nvar_s2m+10) = 'xv'
      sfc_name2(nvar_s2m+11) = 'xz'
      sfc_name2(nvar_s2m+12) = 'zm'
      sfc_name2(nvar_s2m+13) = 'xtts'
      sfc_name2(nvar_s2m+14) = 'xzts'
      sfc_name2(nvar_s2m+15) = 'd_conv'
      sfc_name2(nvar_s2m+16) = 'ifd'
      sfc_name2(nvar_s2m+17) = 'dt_cool'
      sfc_name2(nvar_s2m+18) = 'qrain'
!
! Only needed when Noah MP LSM is used - 29 2D
!
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name2(nvar_s2m+19) = 'snowxy'
        sfc_name2(nvar_s2m+20) = 'tvxy'
        sfc_name2(nvar_s2m+21) = 'tgxy'
        sfc_name2(nvar_s2m+22) = 'canicexy'
        sfc_name2(nvar_s2m+23) = 'canliqxy'
        sfc_name2(nvar_s2m+24) = 'eahxy'
        sfc_name2(nvar_s2m+25) = 'tahxy'
        sfc_name2(nvar_s2m+26) = 'cmxy'
        sfc_name2(nvar_s2m+27) = 'chxy'
        sfc_name2(nvar_s2m+28) = 'fwetxy'
        sfc_name2(nvar_s2m+29) = 'sneqvoxy'
        sfc_name2(nvar_s2m+30) = 'alboldxy'
        sfc_name2(nvar_s2m+31) = 'qsnowxy'
        sfc_name2(nvar_s2m+32) = 'wslakexy'
        sfc_name2(nvar_s2m+33) = 'zwtxy'
        sfc_name2(nvar_s2m+34) = 'waxy'
        sfc_name2(nvar_s2m+35) = 'wtxy'
        sfc_name2(nvar_s2m+36) = 'lfmassxy'
        sfc_name2(nvar_s2m+37) = 'rtmassxy'
        sfc_name2(nvar_s2m+38) = 'stmassxy'
        sfc_name2(nvar_s2m+39) = 'woodxy'
        sfc_name2(nvar_s2m+40) = 'stblcpxy'
        sfc_name2(nvar_s2m+41) = 'fastcpxy'
        sfc_name2(nvar_s2m+42) = 'xsaixy'
        sfc_name2(nvar_s2m+43) = 'xlaixy'
        sfc_name2(nvar_s2m+44) = 'taussxy'
        sfc_name2(nvar_s2m+45) = 'smcwtdxy'
        sfc_name2(nvar_s2m+46) = 'deeprechxy'
        sfc_name2(nvar_s2m+47) = 'rechxy'
        sfc_name2(nvar_s2m+48) = 'drainncprv'
        sfc_name2(nvar_s2m+49) = 'draincprv'
        sfc_name2(nvar_s2m+50) = 'dsnowprv'
        sfc_name2(nvar_s2m+51) = 'dgraupelprv'
        sfc_name2(nvar_s2m+52) = 'diceprv'

      endif
 
      !--- register the 2D fields
      do num = 1,nvar_s2m
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or.trim(sfc_name2(num)) == 'tsfcl'.or.trim(sfc_name2(num)) == 'zorll') then
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        endif
      enddo
      
      
      if (Model%nstf_name(1) > 0) then
        mand = .false.
        if (Model%nstf_name(2) == 0) mand = .true.
        do num = nvar_s2m+1,nvar_s2m+nvar_s2o
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
! Noah MP register only necessary only lsm = 2, not necessary has values
      if (nvar_s2mp > 0) then
        mand = .false.
        do num = nvar_s2m+nvar_s2o+1,nvar_s2m+nvar_s2o+nvar_s2mp
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif ! noahmp

      nullify(var2_p)
    
    
      !--- names of the 3D variables to save
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'
      !--- Noah MP
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name3(4) = 'snicexy'
        sfc_name3(5) = 'snliqxy'
        sfc_name3(6) = 'tsnoxy'
        sfc_name3(7) = 'smoiseq'
        sfc_name3(8) = 'zsnsoxy'
      endif

      !--- register the 3D fields
      do num = 1,nvar_s3
        var3_p => sfc_var3(:,:,:,num)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p, domain=fv_domain)
      enddo
      if (Model%lsm == Model%lsm_noahmp) then
        mand = .false.
        do num = nvar_s3+1,nvar_s3+3
          var3_p1 => sfc_var3sn(:,:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p1, domain=fv_domain,mandatory=mand)
        enddo
      
        var3_p2 => sfc_var3eq(:,:,:,7)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(7), var3_p2, domain=fv_domain,mandatory=mand)
      
        var3_p3 => sfc_var3zn(:,:,:,8)
        id_restart = register_restart_fIeld(Sfc_restart, fn_srf, sfc_name3(8), var3_p3, domain=fv_domain,mandatory=mand)
      
        nullify(var3_p1)
        nullify(var3_p2)
        nullify(var3_p3)
      
      endif   !mp

      nullify(var3_p)
      
    endif  ! if not allocated

!--- Noah MP define arbitrary value (number layers of snow) to indicate
!coldstart(sfcfile doesn't include noah mp fields) or not
    
 
    !--- read the surface restart/data
    fname = 'INPUT/'//trim(fn_srf)
    if (.not. file_exist(fname)) then ! cold start
      call mpp_error(NOTE,'No INPUT/sfc_data.tile*.nc surface data found; cold-starting land surface')
      !Need a namelist for options:
      ! 1. choice of sst (uniform, profiles) --- ML0 should relax to this
      ! 2. Choice of veg, soil type with certain soil T,q,ql
      ! How to fix day of year (for astronomy)?
      !--- place the data into the block GFS containers
      do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- 2D variables
             !--- slmsk
             Sfcprop(nb)%slmsk(ix)  = 0.
             !--- tsfc (tsea in sfc file)
             Sfcprop(nb)%tsfc(ix)   = 300. ! should specify some latitudinal profile
             !--- weasd (sheleg in sfc file)
             Sfcprop(nb)%weasd(ix)  = 0.0
             !--- tg3
             Sfcprop(nb)%tg3(ix)    = 290. !generic value, probably not good; real value latitude-dependent
             !--- zorl
             Sfcprop(nb)%zorl(ix)   = 0.1 ! typical ocean value; different values for different land surfaces (use a lookup table?)
             !--- alvsf
             Sfcprop(nb)%alvsf(ix)  = 0.06
             !--- alvwf
             Sfcprop(nb)%alvwf(ix)  = 0.06
             !--- alnsf
             Sfcprop(nb)%alnsf(ix)  = 0.06
             !--- alnwf
             Sfcprop(nb)%alnwf(ix)  = 0.06
             !--- facsf
             Sfcprop(nb)%facsf(ix)  = 0.0
             !--- facwf
             Sfcprop(nb)%facwf(ix)  = 0.0
             !--- vfrac
             Sfcprop(nb)%vfrac(ix)  = 0.0
             !--- canopy
             Sfcprop(nb)%canopy(ix) = 0.0
             !--- f10m
             Sfcprop(nb)%f10m(ix)   = 0.9
             !--- t2m
             Sfcprop(nb)%t2m(ix)    = Sfcprop(nb)%tsfc(ix)
             !--- q2m
             Sfcprop(nb)%q2m(ix)    = 0.0 ! initially dry atmosphere?
             !--- vtype
             Sfcprop(nb)%vtype(ix)  = 0.0 
             !--- stype
             Sfcprop(nb)%stype(ix)  = 0.0 
             !--- uustar
             Sfcprop(nb)%uustar(ix) = 0.5
             !--- ffmm
             Sfcprop(nb)%ffmm(ix)   = 10.
             !--- ffhh
             Sfcprop(nb)%ffhh(ix)   = 10.
             !--- hice
             Sfcprop(nb)%hice(ix)   = 0.0
             !--- fice
             Sfcprop(nb)%fice(ix)   = 0.0
             !--- tisfc
             Sfcprop(nb)%tisfc(ix)  = Sfcprop(nb)%tsfc(ix)
             !--- tprcp
             Sfcprop(nb)%tprcp(ix)  = 0.0
             !--- srflag
             Sfcprop(nb)%srflag(ix) = 0.0
             !--- snowd (snwdph in the file)
             Sfcprop(nb)%snowd(ix)  = 0.0
             !--- shdmin
             Sfcprop(nb)%shdmin(ix) = 0.0 !this and the next depend on the surface type
             !--- shdmax
             Sfcprop(nb)%shdmax(ix) = 0.0
             !--- slope
             Sfcprop(nb)%slope(ix)  = 0.0 ! also land-surface dependent
             !--- snoalb
             Sfcprop(nb)%snoalb(ix) = 0.0
             !--- sncovr
             Sfcprop(nb)%sncovr(ix) = 0.0
             !
             !--- NSSTM variables
             if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 1)) then
                !--- nsstm tref
                Sfcprop(nb)%tref(ix)    = Sfcprop(nb)%tsfc(ix)
                Sfcprop(nb)%xz(ix)      = 30.0d0
             endif
             if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 0)) then
                !return an error
                call mpp_error(FATAL, 'cold-starting does not support NSST.')
             endif

             !--- 3D variables
             ! these are all set to ocean values.
                !--- stc
                Sfcprop(nb)%stc(ix,:) = Sfcprop(nb)%tsfc(ix)
                !--- smc
                Sfcprop(nb)%smc(ix,:) = 1.0 
                !--- slc
                Sfcprop(nb)%slc(ix,:) = 1.0 
          enddo
       enddo
    else

    !--- read surface properties form sfc_data.tile*.nc

      call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
      call restore_state(Sfc_restart)       
 
      !--- place the data into the block GFS containers
      do nb = 1, Atm_block%nblks
         do ix = 1, Atm_block%blksz(nb)
           i = Atm_block%index(nb)%ii(ix) - isc + 1
           j = Atm_block%index(nb)%jj(ix) - jsc + 1

!--- 2D variables
!    ------------
           Sfcprop(nb)%slmsk(ix)  = sfc_var2(i,j,1)    !--- slmsk
           Sfcprop(nb)%tsfco(ix)  = sfc_var2(i,j,2)    !--- tsfc (tsea in sfc file)
           Sfcprop(nb)%weasd(ix)  = sfc_var2(i,j,3)    !--- weasd (sheleg in sfc file)
           Sfcprop(nb)%tg3(ix)    = sfc_var2(i,j,4)    !--- tg3
           Sfcprop(nb)%zorlo(ix)  = sfc_var2(i,j,5)    !--- zorl on ocean
           Sfcprop(nb)%alvsf(ix)  = sfc_var2(i,j,6)    !--- alvsf
           Sfcprop(nb)%alvwf(ix)  = sfc_var2(i,j,7)    !--- alvwf
           Sfcprop(nb)%alnsf(ix)  = sfc_var2(i,j,8)    !--- alnsf
           Sfcprop(nb)%alnwf(ix)  = sfc_var2(i,j,9)    !--- alnwf
           Sfcprop(nb)%facsf(ix)  = sfc_var2(i,j,10)   !--- facsf
           Sfcprop(nb)%facwf(ix)  = sfc_var2(i,j,11)   !--- facwf
           Sfcprop(nb)%vfrac(ix)  = sfc_var2(i,j,12)   !--- vfrac
           Sfcprop(nb)%canopy(ix) = sfc_var2(i,j,13)   !--- canopy
           Sfcprop(nb)%f10m(ix)   = sfc_var2(i,j,14)   !--- f10m
           Sfcprop(nb)%t2m(ix)    = sfc_var2(i,j,15)   !--- t2m
           Sfcprop(nb)%q2m(ix)    = sfc_var2(i,j,16)   !--- q2m
           Sfcprop(nb)%vtype(ix)  = sfc_var2(i,j,17)   !--- vtype
           Sfcprop(nb)%stype(ix)  = sfc_var2(i,j,18)   !--- stype
           Sfcprop(nb)%uustar(ix) = sfc_var2(i,j,19)   !--- uustar
           Sfcprop(nb)%ffmm(ix)   = sfc_var2(i,j,20)   !--- ffmm
           Sfcprop(nb)%ffhh(ix)   = sfc_var2(i,j,21)   !--- ffhh
           Sfcprop(nb)%hice(ix)   = sfc_var2(i,j,22)   !--- hice
           Sfcprop(nb)%fice(ix)   = sfc_var2(i,j,23)   !--- fice
           Sfcprop(nb)%tisfc(ix)  = sfc_var2(i,j,24)   !--- tisfc
           Sfcprop(nb)%tprcp(ix)  = sfc_var2(i,j,25)   !--- tprcp
           Sfcprop(nb)%srflag(ix) = sfc_var2(i,j,26)   !--- srflag
           Sfcprop(nb)%snowd(ix)  = sfc_var2(i,j,27)   !--- snowd (snwdph in the file)
           Sfcprop(nb)%shdmin(ix) = sfc_var2(i,j,28)   !--- shdmin
           Sfcprop(nb)%shdmax(ix) = sfc_var2(i,j,29)   !--- shdmax
           Sfcprop(nb)%slope(ix)  = sfc_var2(i,j,30)   !--- slope
           Sfcprop(nb)%snoalb(ix) = sfc_var2(i,j,31)   !--- snoalb
           Sfcprop(nb)%sncovr(ix) = sfc_var2(i,j,32)   !--- sncovr
           if(Model%cplflx) then
             Sfcprop(nb)%tsfcl(ix)  = sfc_var2(i,j,33) !--- sfcl  (temp on land portion of a cell)
             Sfcprop(nb)%zorll(ix)  = sfc_var2(i,j,34) !--- zorll (zorl on land portion of a cell)
           end if
           !
           !--- NSSTM variables
           if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 1)) then
             !--- nsstm tref
             Sfcprop(nb)%tref(ix)    = Sfcprop(nb)%tsfco(ix)
             Sfcprop(nb)%xz(ix)      = 30.0d0
           endif
           if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 0)) then
             Sfcprop(nb)%tref(ix)    = sfc_var2(i,j,nvar_s2m+1)  !--- nsstm tref
             Sfcprop(nb)%z_c(ix)     = sfc_var2(i,j,nvar_s2m+2)  !--- nsstm z_c
             Sfcprop(nb)%c_0(ix)     = sfc_var2(i,j,nvar_s2m+3)  !--- nsstm c_0
             Sfcprop(nb)%c_d(ix)     = sfc_var2(i,j,nvar_s2m+4)  !--- nsstm c_d
             Sfcprop(nb)%w_0(ix)     = sfc_var2(i,j,nvar_s2m+5)  !--- nsstm w_0
             Sfcprop(nb)%w_d(ix)     = sfc_var2(i,j,nvar_s2m+6)  !--- nsstm w_d
             Sfcprop(nb)%xt(ix)      = sfc_var2(i,j,nvar_s2m+7)  !--- nsstm xt
             Sfcprop(nb)%xs(ix)      = sfc_var2(i,j,nvar_s2m+8)  !--- nsstm xs
             Sfcprop(nb)%xu(ix)      = sfc_var2(i,j,nvar_s2m+9)  !--- nsstm xu
             Sfcprop(nb)%xv(ix)      = sfc_var2(i,j,nvar_s2m+10) !--- nsstm xv
             Sfcprop(nb)%xz(ix)      = sfc_var2(i,j,nvar_s2m+11) !--- nsstm xz
             Sfcprop(nb)%zm(ix)      = sfc_var2(i,j,nvar_s2m+12) !--- nsstm zm
             Sfcprop(nb)%xtts(ix)    = sfc_var2(i,j,nvar_s2m+13) !--- nsstm xtts
             Sfcprop(nb)%xzts(ix)    = sfc_var2(i,j,nvar_s2m+14) !--- nsstm xzts
             Sfcprop(nb)%d_conv(ix)  = sfc_var2(i,j,nvar_s2m+15) !--- nsstm d_conv
             Sfcprop(nb)%ifd(ix)     = sfc_var2(i,j,nvar_s2m+16) !--- nsstm ifd
             Sfcprop(nb)%dt_cool(ix) = sfc_var2(i,j,nvar_s2m+17) !--- nsstm dt_cool
             Sfcprop(nb)%qrain(ix)   = sfc_var2(i,j,nvar_s2m+18) !--- nsstm qrain
           endif

! Noah MP
! -------
           if (Model%lsm == Model%lsm_noahmp) then
             Sfcprop(nb)%snowxy(ix)     = sfc_var2(i,j,nvar_s2m+19)
             Sfcprop(nb)%tvxy(ix)       = sfc_var2(i,j,nvar_s2m+20)
             Sfcprop(nb)%tgxy(ix)       = sfc_var2(i,j,nvar_s2m+21)
             Sfcprop(nb)%canicexy(ix)   = sfc_var2(i,j,nvar_s2m+22)
             Sfcprop(nb)%canliqxy(ix)   = sfc_var2(i,j,nvar_s2m+23)
             Sfcprop(nb)%eahxy(ix)      = sfc_var2(i,j,nvar_s2m+24)
             Sfcprop(nb)%tahxy(ix)      = sfc_var2(i,j,nvar_s2m+25)
             Sfcprop(nb)%cmxy(ix)       = sfc_var2(i,j,nvar_s2m+26)
             Sfcprop(nb)%chxy(ix)       = sfc_var2(i,j,nvar_s2m+27)
             Sfcprop(nb)%fwetxy(ix)     = sfc_var2(i,j,nvar_s2m+28)
             Sfcprop(nb)%sneqvoxy(ix)   = sfc_var2(i,j,nvar_s2m+29)
             Sfcprop(nb)%alboldxy(ix)   = sfc_var2(i,j,nvar_s2m+30)
             Sfcprop(nb)%qsnowxy(ix)    = sfc_var2(i,j,nvar_s2m+31)
             Sfcprop(nb)%wslakexy(ix)   = sfc_var2(i,j,nvar_s2m+32)
             Sfcprop(nb)%zwtxy(ix)      = sfc_var2(i,j,nvar_s2m+33)
             Sfcprop(nb)%waxy(ix)       = sfc_var2(i,j,nvar_s2m+34)
             Sfcprop(nb)%wtxy(ix)       = sfc_var2(i,j,nvar_s2m+35)
             Sfcprop(nb)%lfmassxy(ix)   = sfc_var2(i,j,nvar_s2m+36)
             Sfcprop(nb)%rtmassxy(ix)   = sfc_var2(i,j,nvar_s2m+37)
             Sfcprop(nb)%stmassxy(ix)   = sfc_var2(i,j,nvar_s2m+38)
             Sfcprop(nb)%woodxy(ix)     = sfc_var2(i,j,nvar_s2m+39)
             Sfcprop(nb)%stblcpxy(ix)   = sfc_var2(i,j,nvar_s2m+40)
             Sfcprop(nb)%fastcpxy(ix)   = sfc_var2(i,j,nvar_s2m+41)
             Sfcprop(nb)%xsaixy(ix)     = sfc_var2(i,j,nvar_s2m+42)
             Sfcprop(nb)%xlaixy(ix)     = sfc_var2(i,j,nvar_s2m+43)
             Sfcprop(nb)%taussxy(ix)    = sfc_var2(i,j,nvar_s2m+44)
             Sfcprop(nb)%smcwtdxy(ix)   = sfc_var2(i,j,nvar_s2m+45)
             Sfcprop(nb)%deeprechxy(ix) = sfc_var2(i,j,nvar_s2m+46)
             Sfcprop(nb)%rechxy(ix)     = sfc_var2(i,j,nvar_s2m+47)
             Sfcprop(nb)%drainncprv(ix) = sfc_var2(i,j,nvar_s2m+48)
             Sfcprop(nb)%draincprv(ix)  = sfc_var2(i,j,nvar_s2m+49)
             Sfcprop(nb)%dsnowprv(ix)   = sfc_var2(i,j,nvar_s2m+50)
             Sfcprop(nb)%dgraupelprv(ix)= sfc_var2(i,j,nvar_s2m+51)
             Sfcprop(nb)%diceprv(ix)    = sfc_var2(i,j,nvar_s2m+52)
           endif


           !--- 3D variables
           do lsoil = 1,Model%lsoil
             Sfcprop(nb)%stc(ix,lsoil) = sfc_var3(i,j,lsoil,1)   !--- stc
             Sfcprop(nb)%smc(ix,lsoil) = sfc_var3(i,j,lsoil,2)   !--- smc
             Sfcprop(nb)%slc(ix,lsoil) = sfc_var3(i,j,lsoil,3)   !--- slc
           enddo
           
           if (Model%lsm == Model%lsm_noahmp) then
             do lsoil = -2, 0
               Sfcprop(nb)%snicexy(ix,lsoil) = sfc_var3sn(i,j,lsoil,4)
               Sfcprop(nb)%snliqxy(ix,lsoil) = sfc_var3sn(i,j,lsoil,5)
               Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc_var3sn(i,j,lsoil,6)
             enddo 
           
             do lsoil = 1, 4
               Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc_var3eq(i,j,lsoil,7)
             enddo 
           
             do lsoil = -2, 4
               Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc_var3zn(i,j,lsoil,8)
             enddo 
           endif

        enddo   !ix
      enddo    !nb    
      
      call mpp_error(NOTE, 'gfs_driver:: - after put to container ')
! so far: At cold start everything is 9999.0, warm start snowxy has values
!         but the 3D of snow fields are not available because not allocated yet.
!         ix,nb loops may be consolidate with the Noah MP isnowxy init
!         restore traditional vars first,we need some of them to init snow fields
!         snow depth to actual snow layers; so we can allocate and register
!         note zsnsoxy is from -2:4 - isnowxy is from 0:-2, but we need
!         exact snow layers to pass 3D fields correctly, snow layers are
!         different fro grid to grid, we have to init point by point/grid.
!         It has to be done after the weasd is available
!         sfc_var2(1,1,32) is the first; we need this to allocate snow related fields
        
      !--- if sncovr does not exist in the restart, need to create it
      if (nint(sfc_var2(1,1,32)) == -9999) then
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr') 
        !--- compute sncovr from existing variables
        !--- code taken directly from read_fix.f
        do nb = 1, Atm_block%nblks
           do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%sncovr(ix) = 0.0
              if (Sfcprop(nb)%slmsk(ix) > 0.001) then
                 vegtyp = Sfcprop(nb)%vtype(ix)
                 if (vegtyp == 0) vegtyp = 7
                 rsnow  = 0.001*Sfcprop(nb)%weasd(ix)/snupx(vegtyp)
                 if (0.001*Sfcprop(nb)%weasd(ix) < snupx(vegtyp)) then
                    Sfcprop(nb)%sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
                 else
                    Sfcprop(nb)%sncovr(ix) = 1.0
                 endif
              endif
           enddo
        enddo
      endif
      
      if(Model%frac_grid) then ! 3-way composite
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            tem = 1.0 - Sfcprop(nb)%landfrac(ix) - Sfcprop(nb)%fice(ix)
            Sfcprop(nb)%zorl(ix) = Sfcprop(nb)%zorll(ix) * Sfcprop(nb)%landfrac(ix) &
                                 + Sfcprop(nb)%zorll(ix) * Sfcprop(nb)%fice(ix)     & !zorl ice = zorl land
                                 + Sfcprop(nb)%zorlo(ix) * tem
            Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix) * Sfcprop(nb)%landfrac(ix) &
                                 + Sfcprop(nb)%tisfc(ix) * Sfcprop(nb)%fice(ix)     &
                                 + Sfcprop(nb)%tsfco(ix) * tem
          enddo
        enddo
      else     ! in this case ice fracion is fraction of water fraction
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
        !--- specify tsfcl/zorll from existing variable tsfco/zorlo
            Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix)
            Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorlo(ix)
            Sfcprop(nb)%zorl(ix)  = Sfcprop(nb)%zorlo(ix)
            Sfcprop(nb)%tsfc(ix)  = Sfcprop(nb)%tsfco(ix)
            if (Sfcprop(nb)%slmsk(ix) > 1.9) then
              Sfcprop(nb)%landfrac(ix) = 0.0
            else
              Sfcprop(nb)%landfrac(ix) = Sfcprop(nb)%slmsk(ix)
            endif
          enddo
        enddo
      endif ! if (Model%frac_grid)

      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%lakefrac(ix) > 0.0) then
            Sfcprop(nb)%oceanfrac(ix) = 0.0 ! lake & ocean don't coexist in a cell
          else
            Sfcprop(nb)%oceanfrac(ix) = 1.0 - Sfcprop(nb)%landfrac(ix)  !LHS:ocean frac [0:1]
          endif

        enddo
      enddo


      if (Model%lsm == Model%lsm_noahmp) then 
        if (sfc_var2(1,1,nvar_s2m+19) < -9990. ) then
        !--- initialize NOAH MP properties
          if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver:: - Cold start Noah MP ')

          do nb = 1, Atm_block%nblks
            do ix = 1, Atm_block%blksz(nb)

              Sfcprop(nb)%tvxy(ix)     = missing_value
              Sfcprop(nb)%tgxy(ix)     = missing_value
              Sfcprop(nb)%tahxy(ix)    = missing_value
              Sfcprop(nb)%canicexy(ix) = missing_value
              Sfcprop(nb)%canliqxy(ix) = missing_value
              Sfcprop(nb)%eahxy(ix)    = missing_value
              Sfcprop(nb)%cmxy(ix)     = missing_value
              Sfcprop(nb)%chxy(ix)     = missing_value
              Sfcprop(nb)%fwetxy(ix)   = missing_value
              Sfcprop(nb)%sneqvoxy(ix) = missing_value
              Sfcprop(nb)%alboldxy(ix) = missing_value
              Sfcprop(nb)%qsnowxy(ix)  = missing_value
              Sfcprop(nb)%wslakexy(ix) = missing_value
              Sfcprop(nb)%taussxy(ix)  = missing_value
              Sfcprop(nb)%waxy(ix)     = missing_value
              Sfcprop(nb)%wtxy(ix)     = missing_value
              Sfcprop(nb)%zwtxy(ix)    = missing_value
              Sfcprop(nb)%xlaixy(ix)   = missing_value
              Sfcprop(nb)%xsaixy(ix)   = missing_value

              Sfcprop(nb)%lfmassxy(ix) = missing_value
              Sfcprop(nb)%stmassxy(ix) = missing_value
              Sfcprop(nb)%rtmassxy(ix) = missing_value
              Sfcprop(nb)%woodxy(ix)   = missing_value
              Sfcprop(nb)%stblcpxy(ix) = missing_value
              Sfcprop(nb)%fastcpxy(ix) = missing_value
              Sfcprop(nb)%smcwtdxy(ix) = missing_value
              Sfcprop(nb)%deeprechxy(ix) = missing_value
              Sfcprop(nb)%rechxy(ix)     = missing_value
              Sfcprop(nb)%drainncprv(ix) = 0.
              Sfcprop(nb)%draincprv(ix)  = 0.
              Sfcprop(nb)%dsnowprv(ix)   = 0.
              Sfcprop(nb)%dgraupelprv(ix)= 0.
              Sfcprop(nb)%diceprv(ix)    = 0.

              Sfcprop(nb)%snowxy (ix)   = missing_value
              Sfcprop(nb)%snicexy(ix, -2:0) = missing_value
              Sfcprop(nb)%snliqxy(ix, -2:0) = missing_value
              Sfcprop(nb)%tsnoxy (ix, -2:0) = missing_value
              Sfcprop(nb)%smoiseq(ix,  1:4) = missing_value
              Sfcprop(nb)%zsnsoxy(ix, -2:4) = missing_value

              if (Sfcprop(nb)%slmsk(ix) > 0.01) then

                Sfcprop(nb)%tvxy(ix)     = Sfcprop(nb)%tsfcl(ix)
                Sfcprop(nb)%tgxy(ix)     = Sfcprop(nb)%tsfcl(ix)
                Sfcprop(nb)%tahxy(ix)    = Sfcprop(nb)%tsfcl(ix)

                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tvxy(ix)  = 273.15
                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tgxy(ix)  = 273.15
                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tahxy(ix) = 273.15
     
                Sfcprop(nb)%canicexy(ix) = 0.0
                Sfcprop(nb)%canliqxy(ix) = Sfcprop(nb)%canopy(ix)

                Sfcprop(nb)%eahxy(ix)    = 2000.0

    !     eahxy = psfc*qv/(0.622+qv); qv is mixing ratio, converted from sepcific
    !     humidity specific humidity /(1.0 - specific humidity)

                Sfcprop(nb)%cmxy(ix)     = 0.0
                Sfcprop(nb)%chxy(ix)     = 0.0
                Sfcprop(nb)%fwetxy(ix)   = 0.0
                Sfcprop(nb)%sneqvoxy(ix) = Sfcprop(nb)%weasd(ix)     ! mm
                Sfcprop(nb)%alboldxy(ix) = 0.65
                Sfcprop(nb)%qsnowxy(ix)  = 0.0

    !          if (Sfcprop(nb)%srflag(ix) > 0.001) Sfcprop(nb)%qsnowxy(ix) = Sfcprop(nb)%tprcp(ix)/Model%dtp
    !already set to 0.0
                Sfcprop(nb)%wslakexy(ix) = 0.0
                Sfcprop(nb)%taussxy(ix)  = 0.0


                Sfcprop(nb)%waxy(ix)     = 4900.0
                Sfcprop(nb)%wtxy(ix)     = Sfcprop(nb)%waxy(ix)
                Sfcprop(nb)%zwtxy(ix)    = (25.0 + 2.0) - Sfcprop(nb)%waxy(ix) / 1000.0 /0.2
    !
                vegtyp                   = Sfcprop(nb)%vtype(ix)
                if (vegtyp == 0) vegtyp = 7
                imn                      = Model%idate(2)

                if ((vegtyp == isbarren_table) .or. (vegtyp == isice_table) .or.  (vegtyp == isurban_table) .or. (vegtyp == iswater_table)) then

                  Sfcprop(nb)%xlaixy(ix)   = 0.0
                  Sfcprop(nb)%xsaixy(ix)   = 0.0

                  Sfcprop(nb)%lfmassxy(ix) = 0.0
                  Sfcprop(nb)%stmassxy(ix) = 0.0
                  Sfcprop(nb)%rtmassxy(ix) = 0.0

                  Sfcprop(nb)%woodxy   (ix) = 0.0       
                  Sfcprop(nb)%stblcpxy (ix) = 0.0      
                  Sfcprop(nb)%fastcpxy (ix) = 0.0     

                else

    !           print *, 'vegtyp', vegtyp
    !           print *, 'imn', imn
    !           print *, 'xlaixy', Sfcprop(nb)%xlaixy(ix) 

                  Sfcprop(nb)%xlaixy(ix)   = max(laim_table(vegtyp, imn),0.05)
    !           Sfcprop(nb)%xsaixy(ix)   = max(saim_table(vegtyp, imn),0.05)
                  Sfcprop(nb)%xsaixy(ix)   = max(Sfcprop(nb)%xlaixy(ix)*0.1,0.05)

                  masslai                  = 1000.0 / max(sla_table(vegtyp),1.0)
                  Sfcprop(nb)%lfmassxy(ix) = Sfcprop(nb)%xlaixy(ix)*masslai
                  masssai                  = 1000.0 / 3.0
                  Sfcprop(nb)%stmassxy(ix) = Sfcprop(nb)%xsaixy(ix)* masssai

                  Sfcprop(nb)%rtmassxy(ix) = 500.0      

                  Sfcprop(nb)%woodxy  (ix) = 500.0       
                  Sfcprop(nb)%stblcpxy(ix) = 1000.0      
                  Sfcprop(nb)%fastcpxy(ix) = 1000.0     

                endif  ! non urban ...

                if ( vegtyp == isice_table )  then
                  do lsoil = 1,Model%lsoil
                    Sfcprop(nb)%stc(ix,lsoil) = min(Sfcprop(nb)%stc(ix,lsoil),min(Sfcprop(nb)%tg3(ix),263.15))
                    Sfcprop(nb)%smc(ix,lsoil) = 1
                    Sfcprop(nb)%slc(ix,lsoil) = 0
                  enddo
                endif

                snd   = Sfcprop(nb)%snowd(ix)/1000.0  ! go to m from snwdph

                if (Sfcprop(nb)%weasd(ix) /= 0.0 .and. snd == 0.0 ) then
                  snd = Sfcprop(nb)%weasd(ix)/1000.0
                endif

                if (vegtyp == 15) then                      ! land ice in MODIS/IGBP
                  if ( Sfcprop(nb)%weasd(ix) < 0.1) then
                    Sfcprop(nb)%weasd(ix) = 0.1
                    snd                   = 0.01
                  endif
                endif

                if (snd < 0.025 ) then
                  Sfcprop(nb)%snowxy(ix)   = 0.0
                  dzsno(-2:0)              = 0.0
                elseif (snd >= 0.025 .and. snd <= 0.05 ) then
                  Sfcprop(nb)%snowxy(ix)   = -1.0
                  dzsno(0)                 = snd
                elseif (snd > 0.05 .and. snd <= 0.10 ) then
                  Sfcprop(nb)%snowxy(ix)   = -2.0
                  dzsno(-1)                = 0.5*snd
                  dzsno(0)                 = 0.5*snd
                elseif (snd > 0.10 .and. snd <= 0.25 ) then
                  Sfcprop(nb)%snowxy(ix)   = -2.0
                  dzsno(-1)                = 0.05
                  dzsno(0)                 = snd - 0.05
                elseif (snd > 0.25 .and. snd <= 0.45 ) then
                  Sfcprop(nb)%snowxy(ix)   = -3.0
                  dzsno(-2)                = 0.05
                  dzsno(-1)                = 0.5*(snd-0.05)
                  dzsno(0)                 = 0.5*(snd-0.05)
                elseif (snd > 0.45) then 
                  Sfcprop(nb)%snowxy(ix)   = -3.0
                  dzsno(-2)                = 0.05
                  dzsno(-1)                = 0.20
                  dzsno(0)                 = snd - 0.05 - 0.20
                else
                  call mpp_error(FATAL, 'problem with the logic assigning snow layers.') 
                endif

    !Now we have the snowxy field
    !snice + snliq + tsno allocation and compute them from what we have
                
    !
                Sfcprop(nb)%tsnoxy(ix,-2:0)  = 0.0
                Sfcprop(nb)%snicexy(ix,-2:0) = 0.0
                Sfcprop(nb)%snliqxy(ix,-2:0) = 0.0
                Sfcprop(nb)%zsnsoxy(ix,-2:4) = 0.0
             
                isnow = nint(Sfcprop(nb)%snowxy(ix))+1    ! snowxy <=0.0, dzsno >= 0.0
             
                do ns = isnow , 0
                  Sfcprop(nb)%tsnoxy(ix,ns)  = Sfcprop(nb)%tgxy(ix)
                  Sfcprop(nb)%snliqxy(ix,ns) = 0.0
                  Sfcprop(nb)%snicexy(ix,ns) = 1.00 * dzsno(ns) * Sfcprop(nb)%weasd(ix)/snd
                enddo
    !
    !zsnsoxy, all negative ?
    !
                do ns = isnow, 0
                  dzsnso(ns) = -dzsno(ns)
                enddo

                do ns = 1 , 4
                  dzsnso(ns) = -dzs(ns)
                enddo
    !
    !Assign to zsnsoxy
    !
                Sfcprop(nb)%zsnsoxy(ix,isnow) = dzsnso(isnow)
                do ns = isnow+1,4
                  Sfcprop(nb)%zsnsoxy(ix,ns) = Sfcprop(nb)%zsnsoxy(ix,ns-1) + dzsnso(ns)
                enddo
     
    !
    !smoiseq
    !Init water table related quantities here
    !
                soiltyp  = Sfcprop(nb)%stype(ix)


                if (soiltyp == 0) then
                  Sfcprop(nb)%stype(ix) = 16
                  soiltyp  = Sfcprop(nb)%stype(ix)
                endif

                bexp   = bexp_table(soiltyp)
                smcmax = smcmax_table(soiltyp)
                smcwlt = smcwlt_table(soiltyp)
                dwsat  = dwsat_table(soiltyp)
                dksat  = dksat_table(soiltyp)
                psisat = -psisat_table(soiltyp)

                if (vegtyp == isurban_table) then
                  smcmax = 0.45
                  smcwlt = 0.40
                endif
                 
                if ((bexp > 0.0) .and. (smcmax > 0.0) .and. (-psisat > 0.0 )) then
                  do ns = 1, Model%lsoil          
                    if ( ns == 1 )then
                      ddz = -zsoil(ns+1) * 0.5
                    elseif ( ns < Model%lsoil ) then
                      ddz = ( zsoil(ns-1) - zsoil(ns+1) ) * 0.5
                    else
                      ddz = zsoil(ns-1) - zsoil(ns)
                    endif
    !
    !Use newton-raphson method to find eq soil moisture
    !
                    expon = bexp + 1.
                    aa    = dwsat / ddz
                    bb    = dksat / smcmax ** expon
                 
                    smc = 0.5 * smcmax
                 
                    do iter = 1, 100
                      func  = (smc - smcmax) * aa +  bb * smc ** expon
                      dfunc = aa + bb * expon * smc ** bexp
                      dx    = func / dfunc
                      smc   = smc - dx
                      if ( abs (dx) < 1.e-6) exit
                    enddo                               ! iteration
                    Sfcprop(nb)%smoiseq(ix,ns) = min(max(smc,1.e-4),smcmax*0.99)
                  enddo                                 ! ddz soil layer
                else                                    ! bexp <= 0.0 
                  Sfcprop(nb)%smoiseq(ix,1:4) = smcmax
                endif                                   ! end the bexp condition
    !
                  Sfcprop(nb)%smcwtdxy(ix)   = smcmax
                  Sfcprop(nb)%deeprechxy(ix) = 0.0
                  Sfcprop(nb)%rechxy(ix)     = 0.0
     
              endif !end if slmsk>0.01 (land only)

            enddo ! ix
          enddo  ! nb
        endif    
      endif !if Noah MP cold start ends    

    endif


  end subroutine sfc_prop_restart_read

  subroutine sfc_prop_override(Sfcprop, Grid, Atm_block, Model, fv_domain)

    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type(GFS_grid_type),       intent(inout) :: Grid(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(in)    :: Model
    type (domain2d),           intent(in)    :: fv_domain
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb
    integer :: isc, iec, jsc, jec, npz, nx, ny, ios

    logical :: ideal_sst = .false.
    real(kind=kind_phys) :: sst_max = 300.
    real(kind=kind_phys) :: sst_min = 271.14 ! -2c --> sea ice
    integer :: sst_profile = 0

    logical :: ideal_land = .false.
    !Assuming modern veg/soil types
    ! sample Amazon settings; values for OKC in comments
    integer :: vegtype = 2 ! 12
    integer :: soiltype = 9 ! 8
    real(kind=kind_phys) :: vegfrac = 0.8 ! 0.25 -- 0.75 
    real(kind=kind_phys) :: zorl = 265 ! 15
    !uniform soil temperature and moisture for now
    real(kind=kind_phys) :: stc = 300. ! 310.
    real(kind=kind_phys) :: smc = 0.4 ! wet season vs. 0.08 dry ! 0.2 okc highly variable and patchy
    
    
    namelist /sfc_prop_override_nml/ &
         ideal_sst, sst_max, sst_profile, & !Aquaplanet SST options
         ideal_land, vegtype, soiltype, & ! idealized soil/veg options
         vegfrac, zorl, stc, smc
    
#ifdef INTERNAL_FILE_NML
    read(Model%input_nml_file, nml=sfc_prop_override_nml, iostat=ios)
#else
!       print *,' in sfcsub nlunit=',nlunit,' me=',me,' ialb=',ialb
    inquire (file=trim(Model%fn_nml), exist=exists)
    if (.not. exists) then
       write(6,*) 'sfc_prop_override:: namelist file: ',trim(Model%fn_nml),' does not exist'
       stop
    else
       open (unit=Model%nlunit, file=Model%fn_nml, READONLY, status='OLD', iostat=ios)
    endif
    rewind(Model%nlunit)
    read (Model%nlunit,sfc_prop_override_nml)
    close (Model%nlunit)
#endif

    call qsmith_init

    call mpp_error(NOTE, "Calling sfc_prop_override")

    if (ideal_sst) then
       do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- slmsk
             Sfcprop(nb)%slmsk(ix)  = 0.0
             !--- tsfc (tsea in sfc file)
             select case (sst_profile)
                case (0) 
                   Sfcprop(nb)%tsfc(ix)   = sst_max
                case (1) ! symmetric
                   Sfcprop(nb)%tsfc(ix)   = sst_min + (sst_max - sst_min)*Grid(nb)%coslat(ix)
                case default
                   call mpp_error(FATAL, "value of sst_profile not defined.")
             end select
             !--- zorl
             Sfcprop(nb)%zorl(ix)   = zorl
             !--- vfrac
             Sfcprop(nb)%vfrac(ix)  = 0.0
             if (Sfcprop(nb)%tsfc(ix) <= sst_min) then
                !--- hice
                Sfcprop(nb)%hice(ix)   = 1.0
                !--- fice
                Sfcprop(nb)%fice(ix)   = 1.0
                Sfcprop(nb)%tsfc(ix)   = sst_min
             else
                !--- hice
                Sfcprop(nb)%hice(ix)   = 0.0
                !--- fice
                Sfcprop(nb)%fice(ix)   = 0.0
             endif
             !--- tisfc
             Sfcprop(nb)%tisfc(ix)  = Sfcprop(nb)%tsfc(ix)             
             !--- t2m ! slt. unstable
             Sfcprop(nb)%t2m(ix)    = Sfcprop(nb)%t2m(ix) * 0.98
             !--- q2m ! use RH = 98% and assume ps = 1000 mb
             Sfcprop(nb)%q2m(ix)    = wqs1 (Sfcprop(nb)%t2m(ix), 1.e5/rd/Sfcprop(nb)%t2m(ix))
             !--- vtype
             Sfcprop(nb)%vtype(ix)  = 0
             !--- stype
             Sfcprop(nb)%stype(ix)  = 0
             !Override MLO properties also
             if (Model%do_ocean) then
                Sfcprop(nb)%ts_clim_iano(ix) = Sfcprop(nb)%tsfc(ix)
                Sfcprop(nb)%tsclim(ix) = Sfcprop(nb)%tsfc(ix)
                Sfcprop(nb)%ts_som(ix) = Sfcprop(nb)%tsfc(ix)
             endif
          enddo
       enddo


    elseif (ideal_land) then
       do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- slmsk
             Sfcprop(nb)%slmsk(ix)  = 1.0
             !--- tsfc (tsea in sfc file)
             Sfcprop(nb)%tsfc(ix)   = stc
             !--- weasd (sheleg in sfc file)
             Sfcprop(nb)%weasd(ix)  = 0.0 ! snow
             !--- tg3
             Sfcprop(nb)%tg3(ix)    = stc ! simple approach
             !--- zorl
             Sfcprop(nb)%zorl(ix)   = zorl
             !--- vfrac
             Sfcprop(nb)%vfrac(ix)  = vegfrac
             !--- canopy
             Sfcprop(nb)%canopy(ix) = 0.0 !this quantity is quite variable
             !--- t2m
             Sfcprop(nb)%t2m(ix)    = stc * 0.98 !slt unstable
             !--- q2m ! use RH = 98%
             Sfcprop(nb)%q2m(ix)    = wqs1 (Sfcprop(nb)%t2m(ix), 1.e5/rd/Sfcprop(nb)%t2m(ix))
             !--- vtype
             Sfcprop(nb)%vtype(ix)  = vegtype
             !--- stype
             Sfcprop(nb)%stype(ix)  = soiltype
             !--- hice
             Sfcprop(nb)%hice(ix)   = 0.0
             !--- fice
             Sfcprop(nb)%fice(ix)   = 0.0
             !--- tisfc
             Sfcprop(nb)%tisfc(ix)  = stc
             !--- snowd (snwdph in the file)
             Sfcprop(nb)%snowd(ix)  = 0.0
             !--- snoalb
             Sfcprop(nb)%snoalb(ix) = 0.5
             !--- sncovr
             Sfcprop(nb)%sncovr(ix) = 0.0
             !--- 3D variables
             do lsoil = 1,Model%lsoil
                !--- stc
                Sfcprop(nb)%stc(ix,lsoil) = stc
                !--- smc
                Sfcprop(nb)%smc(ix,lsoil) = smc
                !--- slc = smc
                Sfcprop(nb)%slc(ix,lsoil) = smc
             enddo
             
          enddo
       enddo

    endif


  end subroutine sfc_prop_override


!----------------------------------------------------------------------      
! sfc_prop_restart_write
!----------------------------------------------------------------------      
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate 
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------      
  subroutine sfc_prop_restart_write (Sfcprop, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, lsoil, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2m, nvar2o, nvar3
    integer :: nvar2mp, nvar3mp
    logical :: mand
    character(len=32) :: fn_srf = 'sfc_data.nc'
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()

    if (Model%cplflx) then ! needs more variables
      nvar2m = 34
    else
    nvar2m = 32
    endif
    nvar2o = 18
    nvar3  = 3
    nvar2mp = 0
    nvar3mp = 0
    if (Model%lsm == Model%lsm_noahmp) then
      nvar2mp = 34
      nvar3mp = 5
    endif

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)

    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar2m+nvar2o))
      allocate(sfc_name3(nvar3))
      allocate(sfc_var2(nx,ny,nvar2m+nvar2o))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
      allocate(sfc_name2(nvar2m+nvar2o+nvar2mp))
      allocate(sfc_name3(nvar3+nvar3mp))
      allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
      sfc_var2 = -9999._kind_phys
      sfc_var3 = -9999._kind_phys
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))

        sfc_var3sn = -9999._kind_phys
        sfc_var3eq = -9999._kind_phys
        sfc_var3zn = -9999._kind_phys
      endif


      !--- names of the 2D variables to save
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsea'    !tsfc
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorl'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
    !--- variables below here are optional
      sfc_name2(32) = 'sncovr'
      if (Model%cplflx) then
        sfc_name2(33) = 'tsfcl'   !temp on land portion of a cell
        sfc_name2(34) = 'zorll'   !zorl on land portion of a cell
      end if
      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      sfc_name2(nvar2m+1)  = 'tref'
      sfc_name2(nvar2m+2)  = 'z_c'
      sfc_name2(nvar2m+3)  = 'c_0'
      sfc_name2(nvar2m+4)  = 'c_d'
      sfc_name2(nvar2m+5)  = 'w_0'
      sfc_name2(nvar2m+6)  = 'w_d'
      sfc_name2(nvar2m+7)  = 'xt'
      sfc_name2(nvar2m+8)  = 'xs'
      sfc_name2(nvar2m+9)  = 'xu'
      sfc_name2(nvar2m+10) = 'xv'
      sfc_name2(nvar2m+11) = 'xz'
      sfc_name2(nvar2m+12) = 'zm'
      sfc_name2(nvar2m+13) = 'xtts'
      sfc_name2(nvar2m+14) = 'xzts'
      sfc_name2(nvar2m+15) = 'd_conv'
      sfc_name2(nvar2m+16) = 'ifd'
      sfc_name2(nvar2m+17) = 'dt_cool'
      sfc_name2(nvar2m+18) = 'qrain'

! Only needed when Noah MP LSM is used - 29 2D
      if(Model%lsm == Model%lsm_noahmp) then
        sfc_name2(nvar2m+19) = 'snowxy'
        sfc_name2(nvar2m+20) = 'tvxy'
        sfc_name2(nvar2m+21) = 'tgxy'
        sfc_name2(nvar2m+22) = 'canicexy'
        sfc_name2(nvar2m+23) = 'canliqxy'
        sfc_name2(nvar2m+24) = 'eahxy'
        sfc_name2(nvar2m+25) = 'tahxy'
        sfc_name2(nvar2m+26) = 'cmxy'
        sfc_name2(nvar2m+27) = 'chxy'
        sfc_name2(nvar2m+28) = 'fwetxy'
        sfc_name2(nvar2m+29) = 'sneqvoxy'
        sfc_name2(nvar2m+30) = 'alboldxy'
        sfc_name2(nvar2m+31) = 'qsnowxy'
        sfc_name2(nvar2m+32) = 'wslakexy'
        sfc_name2(nvar2m+33) = 'zwtxy'
        sfc_name2(nvar2m+34) = 'waxy'
        sfc_name2(nvar2m+35) = 'wtxy'
        sfc_name2(nvar2m+36) = 'lfmassxy'
        sfc_name2(nvar2m+37) = 'rtmassxy'
        sfc_name2(nvar2m+38) = 'stmassxy'
        sfc_name2(nvar2m+39) = 'woodxy'
        sfc_name2(nvar2m+40) = 'stblcpxy'
        sfc_name2(nvar2m+41) = 'fastcpxy'
        sfc_name2(nvar2m+42) = 'xsaixy'
        sfc_name2(nvar2m+43) = 'xlaixy'
        sfc_name2(nvar2m+44) = 'taussxy'
        sfc_name2(nvar2m+45) = 'smcwtdxy'
        sfc_name2(nvar2m+46) = 'deeprechxy'
        sfc_name2(nvar2m+47) = 'rechxy'
        sfc_name2(nvar2m+48) = 'drainncprv'
        sfc_name2(nvar2m+49) = 'draincprv'
        sfc_name2(nvar2m+50) = 'dsnowprv'
        sfc_name2(nvar2m+51) = 'dgraupelprv'
        sfc_name2(nvar2m+52) = 'diceprv'
      endif
 
      !--- register the 2D fields
      do num = 1,nvar2m
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or.trim(sfc_name2(num)) == 'tsfcl'.or.trim(sfc_name2(num)) == 'zorll') then
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        endif
      enddo
      if (Model%nstf_name(1) > 0) then
        mand = .false.
        if (Model%nstf_name(2) ==0) mand = .true.
        do num = nvar2m+1,nvar2m+nvar2o
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
      if (Model%lsm == Model%lsm_noahmp) then
        mand = .true.                  ! actually should be true since it is after cold start
        do num = nvar2m+nvar2o+1,nvar2m+nvar2o+nvar2mp
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
      nullify(var2_p)
 
      !--- names of the 3D variables to save
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name3(4) = 'snicexy'
        sfc_name3(5) = 'snliqxy'
        sfc_name3(6) = 'tsnoxy'
        sfc_name3(7) = 'smoiseq'
        sfc_name3(8) = 'zsnsoxy'
      endif
 
      !--- register the 3D fields
      do num = 1,nvar3
        var3_p => sfc_var3(:,:,:,num)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p, domain=fv_domain)
      enddo
      nullify(var3_p)

      if (Model%lsm == Model%lsm_noahmp) then
        mand = .true.
        do num = nvar3+1,nvar3+3
          var3_p1 => sfc_var3sn(:,:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p1, domain=fv_domain,mandatory=mand)
        enddo

        var3_p2 => sfc_var3eq(:,:,:,7)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(7), var3_p2, domain=fv_domain,mandatory=mand)

        var3_p3 => sfc_var3zn(:,:,:,8)
        id_restart = register_restart_fIeld(Sfc_restart, fn_srf, sfc_name3(8), var3_p3, domain=fv_domain,mandatory=mand)

        nullify(var3_p1)
        nullify(var3_p2)
        nullify(var3_p3)
      endif ! lsm = lsm_noahmp
    endif
   
   
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        !--- 2D variables
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        sfc_var2(i,j,1)  = Sfcprop(nb)%slmsk(ix) !--- slmsk
        sfc_var2(i,j,2)  = Sfcprop(nb)%tsfc(ix)  !--- tsfc (tsea in sfc file)
        sfc_var2(i,j,3)  = Sfcprop(nb)%weasd(ix) !--- weasd (sheleg in sfc file)
        sfc_var2(i,j,4)  = Sfcprop(nb)%tg3(ix)   !--- tg3
        sfc_var2(i,j,5)  = Sfcprop(nb)%zorl(ix)  !--- zorl
        sfc_var2(i,j,6)  = Sfcprop(nb)%alvsf(ix) !--- alvsf
        sfc_var2(i,j,7)  = Sfcprop(nb)%alvwf(ix) !--- alvwf
        sfc_var2(i,j,8)  = Sfcprop(nb)%alnsf(ix) !--- alnsf
        sfc_var2(i,j,9)  = Sfcprop(nb)%alnwf(ix) !--- alnwf
        sfc_var2(i,j,10) = Sfcprop(nb)%facsf(ix) !--- facsf
        sfc_var2(i,j,11) = Sfcprop(nb)%facwf(ix) !--- facwf
        sfc_var2(i,j,12) = Sfcprop(nb)%vfrac(ix) !--- vfrac
        sfc_var2(i,j,13) = Sfcprop(nb)%canopy(ix)!--- canopy
        sfc_var2(i,j,14) = Sfcprop(nb)%f10m(ix)  !--- f10m
        sfc_var2(i,j,15) = Sfcprop(nb)%t2m(ix)   !--- t2m
        sfc_var2(i,j,16) = Sfcprop(nb)%q2m(ix)   !--- q2m
        sfc_var2(i,j,17) = Sfcprop(nb)%vtype(ix) !--- vtype
        sfc_var2(i,j,18) = Sfcprop(nb)%stype(ix) !--- stype
        sfc_var2(i,j,19) = Sfcprop(nb)%uustar(ix)!--- uustar
        sfc_var2(i,j,20) = Sfcprop(nb)%ffmm(ix)  !--- ffmm
        sfc_var2(i,j,21) = Sfcprop(nb)%ffhh(ix)  !--- ffhh
        sfc_var2(i,j,22) = Sfcprop(nb)%hice(ix)  !--- hice
        sfc_var2(i,j,23) = Sfcprop(nb)%fice(ix)  !--- fice
        sfc_var2(i,j,24) = Sfcprop(nb)%tisfc(ix) !--- tisfc
        sfc_var2(i,j,25) = Sfcprop(nb)%tprcp(ix) !--- tprcp
        sfc_var2(i,j,26) = Sfcprop(nb)%srflag(ix)!--- srflag
        sfc_var2(i,j,27) = Sfcprop(nb)%snowd(ix) !--- snowd (snwdph in the file)
        sfc_var2(i,j,28) = Sfcprop(nb)%shdmin(ix)!--- shdmin
        sfc_var2(i,j,29) = Sfcprop(nb)%shdmax(ix)!--- shdmax
        sfc_var2(i,j,30) = Sfcprop(nb)%slope(ix) !--- slope
        sfc_var2(i,j,31) = Sfcprop(nb)%snoalb(ix)!--- snoalb
        sfc_var2(i,j,32) = Sfcprop(nb)%sncovr(ix)!--- sncovr
        if (Model%cplflx) then
          sfc_var2(i,j,33) = Sfcprop(nb)%tsfcl(ix) !--- tsfcl (temp on land)
          sfc_var2(i,j,34) = Sfcprop(nb)%zorll(ix) !--- zorll (zorl on land)
        end if
        !--- NSSTM variables
        if (Model%nstf_name(1) > 0) then
          sfc_var2(i,j,nvar2m+1) = Sfcprop(nb)%tref(ix)    !--- nsstm tref
          sfc_var2(i,j,nvar2m+2) = Sfcprop(nb)%z_c(ix)     !--- nsstm z_c
          sfc_var2(i,j,nvar2m+3) = Sfcprop(nb)%c_0(ix)     !--- nsstm c_0
          sfc_var2(i,j,nvar2m+4) = Sfcprop(nb)%c_d(ix)     !--- nsstm c_d
          sfc_var2(i,j,nvar2m+5) = Sfcprop(nb)%w_0(ix)     !--- nsstm w_0
          sfc_var2(i,j,nvar2m+6) = Sfcprop(nb)%w_d(ix)     !--- nsstm w_d
          sfc_var2(i,j,nvar2m+7) = Sfcprop(nb)%xt(ix)      !--- nsstm xt
          sfc_var2(i,j,nvar2m+8) = Sfcprop(nb)%xs(ix)      !--- nsstm xs
          sfc_var2(i,j,nvar2m+9) = Sfcprop(nb)%xu(ix)      !--- nsstm xu
          sfc_var2(i,j,nvar2m+10) = Sfcprop(nb)%xv(ix)     !--- nsstm xv
          sfc_var2(i,j,nvar2m+11) = Sfcprop(nb)%xz(ix)     !--- nsstm xz
          sfc_var2(i,j,nvar2m+12) = Sfcprop(nb)%zm(ix)     !--- nsstm zm
          sfc_var2(i,j,nvar2m+13) = Sfcprop(nb)%xtts(ix)   !--- nsstm xtts
          sfc_var2(i,j,nvar2m+14) = Sfcprop(nb)%xzts(ix)   !--- nsstm xzts
          sfc_var2(i,j,nvar2m+15) = Sfcprop(nb)%d_conv(ix) !--- nsstm d_conv
          sfc_var2(i,j,nvar2m+16) = Sfcprop(nb)%ifd(ix)    !--- nsstm ifd
          sfc_var2(i,j,nvar2m+17) = Sfcprop(nb)%dt_cool(ix)!--- nsstm dt_cool
          sfc_var2(i,j,nvar2m+18) = Sfcprop(nb)%qrain(ix)  !--- nsstm qrain
        endif
 
! Noah MP
        if (Model%lsm == Model%lsm_noahmp) then
          sfc_var2(i,j,nvar2m+19) = Sfcprop(nb)%snowxy(ix)
          sfc_var2(i,j,nvar2m+20) = Sfcprop(nb)%tvxy(ix)
          sfc_var2(i,j,nvar2m+21) = Sfcprop(nb)%tgxy(ix)
          sfc_var2(i,j,nvar2m+22) = Sfcprop(nb)%canicexy(ix)
          sfc_var2(i,j,nvar2m+23) = Sfcprop(nb)%canliqxy(ix)
          sfc_var2(i,j,nvar2m+24) = Sfcprop(nb)%eahxy(ix)
          sfc_var2(i,j,nvar2m+25) = Sfcprop(nb)%tahxy(ix)
          sfc_var2(i,j,nvar2m+26) = Sfcprop(nb)%cmxy(ix)
          sfc_var2(i,j,nvar2m+27) = Sfcprop(nb)%chxy(ix)
          sfc_var2(i,j,nvar2m+28) = Sfcprop(nb)%fwetxy(ix)
          sfc_var2(i,j,nvar2m+29) = Sfcprop(nb)%sneqvoxy(ix)
          sfc_var2(i,j,nvar2m+30) = Sfcprop(nb)%alboldxy(ix)
          sfc_var2(i,j,nvar2m+31) = Sfcprop(nb)%qsnowxy(ix)
          sfc_var2(i,j,nvar2m+32) = Sfcprop(nb)%wslakexy(ix)
          sfc_var2(i,j,nvar2m+33) = Sfcprop(nb)%zwtxy(ix)
          sfc_var2(i,j,nvar2m+34) = Sfcprop(nb)%waxy(ix)
          sfc_var2(i,j,nvar2m+35) = Sfcprop(nb)%wtxy(ix)
          sfc_var2(i,j,nvar2m+36) = Sfcprop(nb)%lfmassxy(ix)
          sfc_var2(i,j,nvar2m+37) = Sfcprop(nb)%rtmassxy(ix)
          sfc_var2(i,j,nvar2m+38) = Sfcprop(nb)%stmassxy(ix)
          sfc_var2(i,j,nvar2m+39) = Sfcprop(nb)%woodxy(ix)
          sfc_var2(i,j,nvar2m+40) = Sfcprop(nb)%stblcpxy(ix)
          sfc_var2(i,j,nvar2m+41) = Sfcprop(nb)%fastcpxy(ix)
          sfc_var2(i,j,nvar2m+42) = Sfcprop(nb)%xsaixy(ix)
          sfc_var2(i,j,nvar2m+43) = Sfcprop(nb)%xlaixy(ix)
          sfc_var2(i,j,nvar2m+44) = Sfcprop(nb)%taussxy(ix)
          sfc_var2(i,j,nvar2m+45) = Sfcprop(nb)%smcwtdxy(ix)
          sfc_var2(i,j,nvar2m+46) = Sfcprop(nb)%deeprechxy(ix)
          sfc_var2(i,j,nvar2m+47) = Sfcprop(nb)%rechxy(ix)
          sfc_var2(i,j,nvar2m+48) = Sfcprop(nb)%drainncprv(ix)
          sfc_var2(i,j,nvar2m+49) = Sfcprop(nb)%draincprv(ix)
          sfc_var2(i,j,nvar2m+50) = Sfcprop(nb)%dsnowprv(ix)
          sfc_var2(i,j,nvar2m+51) = Sfcprop(nb)%dgraupelprv(ix)
          sfc_var2(i,j,nvar2m+52) = Sfcprop(nb)%diceprv(ix)
        endif

        !--- 3D variables
        do lsoil = 1,Model%lsoil
          sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%stc(ix,lsoil) !--- stc
          sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smc(ix,lsoil) !--- smc
          sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%slc(ix,lsoil) !--- slc
        enddo
! 5 Noah MP 3D
        if (Model%lsm == Model%lsm_noahmp) then

          do lsoil = -2,0
            sfc_var3sn(i,j,lsoil,4) = Sfcprop(nb)%snicexy(ix,lsoil)
            sfc_var3sn(i,j,lsoil,5) = Sfcprop(nb)%snliqxy(ix,lsoil)
            sfc_var3sn(i,j,lsoil,6) = Sfcprop(nb)%tsnoxy(ix,lsoil)
          enddo

          do lsoil = 1,Model%lsoil
            sfc_var3eq(i,j,lsoil,7)  = Sfcprop(nb)%smoiseq(ix,lsoil)
          enddo

          do lsoil = -2,4
            sfc_var3zn(i,j,lsoil,8)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
          enddo

        endif  ! Noah MP
      enddo
    enddo
    call save_restart(Sfc_restart, timestamp)

  end subroutine sfc_prop_restart_write


!----------------------------------------------------------------------      
! phys_restart_read
!----------------------------------------------------------------------      
!    creates and populates a data type which is then used to "register"
!    restart variables with the GFDL FMS restart subsystem.
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  register_restart_field, restart_state, free_restart
!   
!    opens:  phys_data.tile?.nc
!   
!----------------------------------------------------------------------      
  subroutine phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d
    character(len=64) :: fname
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()


    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d
 
    !--- register the restart fields
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = 0.0_kind_phys
      phy_var3 = 0.0_kind_phys
      
      do num = 1,nvar2d
        var2_p => phy_var2(:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_Restart%name2d(num)), &
                                             var2_p, domain=fv_domain, mandatory=.false.)
      enddo
      do num = 1,nvar3d
        var3_p => phy_var3(:,:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_restart%name3d(num)), &
                                             var3_p, domain=fv_domain, mandatory=.false.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    fname = 'INPUT/'//trim(fn_phy)
    if (file_exist(fname)) then
      !--- read the surface restart/data
      call mpp_error(NOTE,'reading physics restart data from INPUT/phy_data.tile*.nc')
      call restore_state(Phy_restart)
    else
      call mpp_error(NOTE,'No physics restarts - cold starting physical parameterizations')
      return
    endif
 
    !--- place the data into the block GFS containers
    !--- phy_var* variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)            
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          IPD_Restart%data(nb,num)%var2p(ix) = phy_var2(i,j,num)
        enddo
      enddo
    enddo
    do num = 1,nvar3d
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)            
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            IPD_Restart%data(nb,num)%var3p(ix,k) = phy_var3(i,j,k,num)
          enddo
        enddo
      enddo
    enddo

  end subroutine phys_restart_read


!----------------------------------------------------------------------      
! phys_restart_write
!----------------------------------------------------------------------      
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate 
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------      
  subroutine phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()


    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d

    !--- register the restart fields 
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = 0.0_kind_phys
      phy_var3 = 0.0_kind_phys
      
      do num = 1,nvar2d
        var2_p => phy_var2(:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_Restart%name2d(num)), &
                                             var2_p, domain=fv_domain, mandatory=.false.)
      enddo
      do num = 1,nvar3d
        var3_p => phy_var3(:,:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_restart%name3d(num)), &
                                             var3_p, domain=fv_domain, mandatory=.false.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    !--- 2D variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)            
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          phy_var2(i,j,num) = IPD_Restart%data(nb,num)%var2p(ix)
        enddo
      enddo
    enddo
    !--- 3D variables
    do num = 1,nvar3d
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)            
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            phy_var3(i,j,k,num) = IPD_Restart%data(nb,num)%var3p(ix,k)
          enddo
        enddo
      enddo
    enddo

    call save_restart(Phy_restart, timestamp)

  end subroutine phys_restart_write

!-------------------------------------------------------------------------      
!--- gfdl_diag_register ---
!-------------------------------------------------------------------------      
!    creates and populates a data type which is then used to "register"
!    GFS physics diagnostic variables with the GFDL FMS diagnostic manager.
!    includes short & long names, units, conversion factors, etc.
!    there is no copying of data, but instead a clever use of pointers.
!    calls a GFDL FMS routine to register diagnositcs and compare against
!    the diag_table to determine what variables are to be output.
!
!    calls:  register_diag_field
!-------------------------------------------------------------------------      
!    Current sizes
!    13+NFXR - radiation
!    76+pl_coeff - physics
!-------------------------------------------------------------------------      
  subroutine gfdl_diag_register(Time, Sfcprop, Gfs_diag, Atm_block, axes, NFXR, ldiag3d, nkld)
    use physcons,  only: con_g
!--- subroutine interface variable definitions
    type(time_type),           intent(in) :: Time
    type(Gfs_sfcprop_type),    intent(in) :: Sfcprop(:)
    type(GFS_diag_type),       intent(in) :: Gfs_diag(:)
    type (block_control_type), intent(in) :: Atm_block
    integer, dimension(4),     intent(in) :: axes
    integer,                   intent(in) :: NFXR
    logical,                   intent(in) :: ldiag3d
    integer,                   intent(in) :: nkld
!--- local variables
    integer :: idx, num, nb, nblks, nx, ny, k
    integer, allocatable :: blksz(:)
    character(len=2) :: xtra
    real(kind=kind_phys), parameter :: cn_one = 1._kind_phys
    real(kind=kind_phys), parameter :: cn_100 = 100._kind_phys
    real(kind=kind_phys), parameter :: cn_th  = 1000._kind_phys
    real(kind=kind_phys), parameter :: cn_hr  = 3600._kind_phys

    nblks = Atm_block%nblks
    allocate (blksz(nblks))
    blksz(:) = Atm_block%blksz(:)

    Diag(:)%id = -99
    Diag(:)%axes = -99
    Diag(:)%cnvfac = 1.0_kind_phys
    Diag(:)%time_avg = .FALSE.

    idx = 0 

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ALBDOsfc'
    Diag(idx)%desc = 'surface albedo (%)'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%fluxr(:,3)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DLWRFsfc'
    Diag(idx)%desc = 'surface downward longwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,19)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRFsfc'
    Diag(idx)%desc = 'surface upward longwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,20)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRFsfc'
    Diag(idx)%desc = 'surface downward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRFsfc'
    Diag(idx)%desc = 'surface upward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRFtoa'
    Diag(idx)%desc = 'top of atmos downward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,23)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRFtoa'
    Diag(idx)%desc = 'top of atmos upward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRFtoa'
    Diag(idx)%desc = 'top of atmos upward longwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCclm'
    Diag(idx)%desc = 'atmos column total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,17)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDChcl'
    Diag(idx)%desc = 'high cloud level total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDClcl'
    Diag(idx)%desc = 'low cloud level total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCmcl'
    Diag(idx)%desc = 'mid cloud level total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,6)
    enddo

!--- accumulated diagnostics ---
    do num = 1,NFXR
      write (xtra,'(I2.2)') num 
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'fluxr_'//trim(xtra)
      Diag(idx)%desc = 'fluxr diagnostic '//trim(xtra)//' - GFS radiation'
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,num)
      enddo
    enddo

!--- averaged diagnostics ---
    do num = 1,nkld
      write (xtra,'(I2.2)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'cloud_'//trim(xtra)
      Diag(idx)%desc = 'cloud diagnostic '//trim(xtra)//' - GFS radiation'
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cloud(:,:,num)
      enddo
    enddo

!--- the next two appear to be appear to be coupling fields in gloopr
!--- each has four elements
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num 
!rab      idx = idx + 1
!rab      Diag(idx)%axes = 2
!rab      Diag(idx)%name = 'dswcmp_'//trim(xtra)
!rab      Diag(idx)%desc = 'dswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      Diag(idx)%unit = 'XXX'
!rab      Diag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswcmp(:,num)
!rab      enddo
!rab    enddo
!rab
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num 
!rab      idx = idx + 1
!rab      Diag(idx)%axes = 2
!rab      Diag(idx)%name = 'uswcmp_'//trim(xtra)
!rab      Diag(idx)%desc = 'uswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      Diag(idx)%unit = 'XXX'
!rab      Diag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswcmp(:,num)
!rab      enddo
!rab    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfxc'
    Diag(idx)%desc = 'total sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%upfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_dnfxc'
    Diag(idx)%desc = 'total sky downward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%dnfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfx0'
    Diag(idx)%desc = 'clear sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%upfx0
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lw_upfxc'
    Diag(idx)%desc = 'total sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topflw(:)%upfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lw_upfx0'
    Diag(idx)%desc = 'clear sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topflw(:)%upfx0
    enddo

!--- physics accumulated diagnostics ---
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'srunoff'
    Diag(idx)%desc = 'surface water runoff - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%srunoff(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evbsa'
    Diag(idx)%desc = 'evbsa - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%evbsa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evcwa'
    Diag(idx)%desc = 'evcwa - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%evcwa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snohfa'
    Diag(idx)%desc = 'snohfa - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snohfa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'transa'
    Diag(idx)%desc = 'transa - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%transa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sbsnoa'
    Diag(idx)%desc = 'sbsnoa - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%sbsnoa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snowca'
    Diag(idx)%desc = 'snowca - GFS lsm'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snowca(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'soilm'
    Diag(idx)%desc = 'total column soil moisture content [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%soilm(:)
      Diag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmin'
    Diag(idx)%desc = 'min temperature at 2m height'
    Diag(idx)%unit = 'k'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tmpmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmax'
    Diag(idx)%desc = 'max temperature at 2m height'
    Diag(idx)%unit = 'k'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tmpmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfc'
    Diag(idx)%desc = 'surface zonal momentum flux [N/m**2]'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dusfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfc'
    Diag(idx)%desc = 'surface meridional momentum flux [N/m**2]'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dtsfc'
    Diag(idx)%desc = 'surface sensible heat flux [W/m**2]'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dtsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dqsfc'
    Diag(idx)%desc = 'surface latent heat flux [W/m**2]'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dqsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totprcp'
    Diag(idx)%desc = 'surface precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gflux'
    Diag(idx)%desc = 'surface ground heat flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%gflux(:)
      Diag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dlwsfc'
    Diag(idx)%desc = 'time accumulated downward lw flux at surface- GFS physics'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ulwsfc'
    Diag(idx)%desc = 'time accumulated upward lw flux at surface- GFS physics'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'suntim'
    Diag(idx)%desc = 'sunshine duration time'
    Diag(idx)%unit = 's'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%suntim(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'runoff'
    Diag(idx)%desc = 'total water runoff'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%runoff(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ep'
    Diag(idx)%desc = 'potential evaporation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ep(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cldwrk'
    Diag(idx)%desc = 'cloud workfunction (valid only with sas)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cldwrk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dugwd'
    Diag(idx)%desc = 'surface zonal gravity wave stress [N/m**2]'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dugwd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvgwd'
    Diag(idx)%desc = 'surface meridional gravity wave stress [N/m**2]'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvgwd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psmean'
    Diag(idx)%desc = 'surface pressure'
    Diag(idx)%unit = 'kPa'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%psmean(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cnvprcp'
    Diag(idx)%desc = 'surface convective precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cnvprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmin'
    Diag(idx)%desc = 'minimum specific humidity'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%spfhmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmax'
    Diag(idx)%desc = 'maximum specific humidity'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%spfhmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u10mmax'
    Diag(idx)%desc = 'maximum (magnitude) u-wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u10mmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v10mmax'
    Diag(idx)%desc = 'maximum (magnitude) v-wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v10mmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wind10mmax'
    Diag(idx)%desc = 'maximum wind speed'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%wind10mmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rain'
    Diag(idx)%desc = 'total rain at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rain(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rainc'
    Diag(idx)%desc = 'convective rain at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rainc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ice'
    Diag(idx)%desc = 'ice fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snow'
    Diag(idx)%desc = 'snow fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snow(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'graupel'
    Diag(idx)%desc = 'graupel fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%graupel(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totice'
    Diag(idx)%desc = 'surface ice precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totsnw'
    Diag(idx)%desc = 'surface snow precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totsnw(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totgrp'
    Diag(idx)%desc = 'surface graupel precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totgrp(:)
    enddo

!--- physics instantaneous diagnostics ---
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u10m'
    Diag(idx)%desc = '10 meter u wind [m/s]'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v10m'
    Diag(idx)%desc = '10 meter v wind [m/s]'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dpt2m'
    Diag(idx)%desc = '2 meter dew point temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dpt2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'zlvl'
    Diag(idx)%desc = 'layer 1 height'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%zlvl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psurf'
    Diag(idx)%desc = 'surface pressure [Pa]'
    Diag(idx)%unit = 'Pa'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%psurf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hpbl'
    Diag(idx)%desc = 'surface planetary boundary layer height [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hpbl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hgamt'
    Diag(idx)%desc = 'ysu counter-gradient heat flux factor'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hgamt(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hfxpbl'
    Diag(idx)%desc = 'ysu entrainment heat flux factor'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hfxpbl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pwat'
    Diag(idx)%desc = 'atmos column precipitable water [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%pwat(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 't1'
    Diag(idx)%desc = 'layer 1 temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%t1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'q1'
    Diag(idx)%desc = 'layer 1 specific humidity'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%q1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u1'
    Diag(idx)%desc = 'layer 1 zonal wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v1'
    Diag(idx)%desc = 'layer 1 meridional wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'chh'
    Diag(idx)%desc = 'thermal exchange coefficient'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%chh(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cmm'
    Diag(idx)%desc = 'momentum exchange coefficient'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cmm(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dlwsfci'
    Diag(idx)%desc = 'instantaneous sfc downward lw flux'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ulwsfci'
    Diag(idx)%desc = 'instantaneous sfc upward lw flux'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dswsfci'
    Diag(idx)%desc = 'instantaneous sfc downward sw flux'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'uswsfci'
    Diag(idx)%desc = 'instantaneous sfc upward sw flux'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfci'
    Diag(idx)%desc = 'instantaneous u component of surface stress'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dusfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfci'
    Diag(idx)%desc = 'instantaneous v component of surface stress'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dtsfci'
    Diag(idx)%desc = 'instantaneous surface sensible heat flux'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dtsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dqsfci'
    Diag(idx)%desc = 'instantaneous surface latent heat flux'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dqsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gfluxi'
    Diag(idx)%desc = 'instantaneous surface ground heat flux'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%gfluxi(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'epi'
    Diag(idx)%desc = 'instantaneous surface potential evaporation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%epi(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'smcwlt2'
    Diag(idx)%desc = 'wiltimg point (volumetric)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%smcwlt2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'smcref2'
    Diag(idx)%desc = 'soil moisture threshold (volumetric)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%smcref2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wet1'
    Diag(idx)%desc = 'normalized soil wetness'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%wet1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sr'
    Diag(idx)%desc = 'ratio of snow to total precipitation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%sr(:)
    enddo

!--- three-dimensional variables that need to be handled special when writing 
    if (ldiag3d) then

    do num = 1,6
      write (xtra,'(I1)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dt3dt_'//trim(xtra)
      Diag(idx)%desc = 'temperature change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dt3dt(:,:,num)
      enddo
    enddo

    do num = 1,5+oz_coeff
      write (xtra,'(I1)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dq3dt_'//trim(xtra)
      Diag(idx)%desc = 'moisture change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dq3dt(:,:,num)
      enddo
    enddo

    do num = 1,4
      write (xtra,'(I1)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'du3dt_'//trim(xtra)
      Diag(idx)%desc = 'u momentum change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%du3dt(:,:,num)
      enddo
    enddo

    do num = 1,4
      write (xtra,'(I1)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dv3dt_'//trim(xtra)
      Diag(idx)%desc = 'v momentum change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
         Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dv3dt(:,:,num)
      enddo
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dkt_pbl'
    Diag(idx)%desc = 'instantaneous heat diffusion coefficient'
    Diag(idx)%unit = 'm**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dkt(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'flux_cg'
    Diag(idx)%desc = 'instantaneous counter-gradient heat flux in ysu'
    Diag(idx)%unit = 'K*m/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%flux_cg(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'flux_en'
    Diag(idx)%desc = 'instantaneous entrainment heat flux in ysu'
    Diag(idx)%unit = 'K*m/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%flux_en(:,:)
    enddo

!!$    idx = idx + 1
!!$    Diag(idx)%axes = 3
!!$    !Requires lgocart = .T.
!!$    Diag(idx)%name = 'dqdt_v'
!!$    Diag(idx)%desc = 'instantaneous total moisture tendency'
!!$    Diag(idx)%unit = 'XXX'
!!$    Diag(idx)%mod_name = 'gfs_phys'
!!$    allocate (Diag(idx)%data(nblks))
!!$    do nb = 1,nblks
!!$       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%Diag(nb)%dqdt_v(:,:,num)
!!$    enddo

!
!--- prognostic variable tendencies (T, u, v, sph, clwmr, o3)
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dtemp_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics temperature tendency'
!rab    Diag(idx)%unit = 'K/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'du_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics horizontal wind component tendency'
!rab    Diag(idx)%unit = 'm/s/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dv_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics meridional wind component tendency'
!rab    Diag(idx)%unit = 'm/s/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dsphum_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics specific humidity tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dclwmr_dt'
!rab    Diag(idx)%desc = 'GFS radiation/radiation cloud water mixing ratio tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'do3mr_dt'
!rab    Diag(idx)%desc = 'GFS radiation/radiation ozone mixing ratio tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'


    endif

!--- Surface diagnostics in gfs_sfc
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alnsf'
    Diag(idx)%desc = 'mean nir albedo with strong cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alnsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alnwf'
    Diag(idx)%desc = 'mean nir albedo with weak cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alnwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alvsf'
    Diag(idx)%desc = 'mean vis albedo with strong cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alvsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alvwf'
    Diag(idx)%desc = 'mean vis albedo with weak cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alvwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'canopy'
    Diag(idx)%desc = 'canopy water (cnwat in gfs data)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%canopy(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'f10m'
    Diag(idx)%desc = '10-meter wind speed divided by lowest model wind speed'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%f10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'facsf'
    Diag(idx)%desc = 'fractional coverage with strong cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%facsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'facwf'
    Diag(idx)%desc = 'fractional coverage with weak cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%facwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ffhh'
    Diag(idx)%desc = 'fh parameter from PBL scheme'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%ffhh(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ffmm'
    Diag(idx)%desc = 'fm parameter from PBL scheme'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%ffmm(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'uustar'
    Diag(idx)%desc = 'uustar surface frictional wind'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%uustar(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'slope'
    Diag(idx)%desc = 'surface slope type'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slope(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'fice'
    Diag(idx)%desc = 'surface ice concentration (ice=1; no ice=0) [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%fice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hice'
    Diag(idx)%desc = 'sea ice thickness (icetk in gfs_data)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%hice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snoalb'
    Diag(idx)%desc = 'maximum snow albedo in fraction (salbd?? in gfs data)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%snoalb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shdmax'
    Diag(idx)%desc = 'maximum fractional coverage of green vegetation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shdmin'
    Diag(idx)%desc = 'minimum fractional coverage of green vegetation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snowd'
    Diag(idx)%desc = 'surface snow depth [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_one/cn_th
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%snowd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'stype'
    Diag(idx)%desc = 'soil type in integer 1-9'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stype(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'q2m'
    Diag(idx)%desc = '2m specific humidity [kg/kg]'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%q2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 't2m'
    Diag(idx)%desc = '2m temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%t2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tsfc'
    Diag(idx)%desc = 'surface temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'qsfc'
    Diag(idx)%desc = 'surface specific humidity [kg/kg]'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%qsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tg3'
    Diag(idx)%desc = 'deep soil temperature'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tg3(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tisfc'
    Diag(idx)%desc = 'surface temperature over ice fraction'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tisfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tprcp'
    Diag(idx)%desc = 'total precipitation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'vtype'
    Diag(idx)%desc = 'vegetation type in integer 1-13'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%vtype(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'weasd'
    Diag(idx)%desc = 'surface snow water equivalent [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%weasd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'HGTsfc'
    Diag(idx)%desc = 'surface geopotential height [gpm]'
    Diag(idx)%unit = 'gpm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = con_g
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%oro(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SLMSKsfc'
    Diag(idx)%desc = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ZORLsfc'
    Diag(idx)%desc = 'surface roughness [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_one/cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%zorl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'VFRACsfc'
    Diag(idx)%desc = 'vegetation fraction'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%vfrac(:)
    enddo

    do num = 1,4
      write (xtra,'(I1)') num 
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'slc_'//trim(xtra)
      Diag(idx)%desc = 'liquid soil mositure at layer-'//trim(xtra)
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_sfc'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slc(:,num)
      enddo
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW1'
    Diag(idx)%desc = 'volumetric soil moisture 0-10cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW2'
    Diag(idx)%desc = 'volumetric soil moisture 10-40cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW3'
    Diag(idx)%desc = 'volumetric soil moisture 40-100cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW4'
    Diag(idx)%desc = 'volumetric soil moisture 100-200cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT1'
    Diag(idx)%desc = 'soil temperature 0-10cm [K]' 
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT2'
    Diag(idx)%desc = 'soil temperature 10-40cm [K]' 
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT3'
    Diag(idx)%desc = 'soil temperature 40-100cm [K]' 
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT4'
    Diag(idx)%desc = 'soil temperature 100-200cm [K]' 
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,4)
    enddo
!bqx+
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'netflxsfc'
    Diag(idx)%desc = 'net surface heat flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%netflxsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'qflux_restore'
    Diag(idx)%desc = 'restoring flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%qflux_restore(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tclim_iano'
    Diag(idx)%desc = 'climatological SST plus initial anomaly'
    Diag(idx)%unit = 'degree C'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tclim_iano(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'MLD'
    Diag(idx)%desc = 'ocean mixed layer depth'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%mld(:)
    enddo

    tot_diag_idx = idx

    if (idx > DIAG_SIZE) then
      call mpp_error(FATAL, 'gfs_driver::gfs_diag_register - need to increase DIAG_SIZE') 
    endif

    do idx = 1,tot_diag_idx
      if (diag(idx)%axes == -99) then
        call mpp_error(FATAL, 'gfs_driver::gfs_diag_register - attempt to register an undefined variable')
      endif
      Diag(idx)%id = register_diag_field (trim(Diag(idx)%mod_name), trim(Diag(idx)%name),  &
                                          axes(1:Diag(idx)%axes), Time, trim(Diag(idx)%desc), &
                                          trim(Diag(idx)%unit), missing_value=real(missing_value))
    enddo

  end subroutine gfdl_diag_register
!-------------------------------------------------------------------------      


!-------------------------------------------------------------------------      
!--- gfs_diag_output ---
!-------------------------------------------------------------------------      
!    routine to transfer the diagnostic data to the GFDL FMS diagnostic 
!    manager for eventual output to the history files.
!
!    calls:  send_data
!-------------------------------------------------------------------------      
!rab  subroutine gfdl_diag_output(Time, Gfs_diag, Statein, Stateout, Atm_block, &
!rab                             nx, ny, levs, ntcw, ntoz, dt, time_int)
  subroutine gfdl_diag_output(Time, Atm_block, IPD_Data, nx, ny, fprint, &
                             levs, ntcw, ntoz, dt, time_int, fhswr, fhlwr)
!--- subroutine interface variable definitions
    logical :: fprint
    type(time_type),           intent(in) :: Time
!rab    type(diagnostics),         intent(in) :: Gfs_diag
!rab    type(state_fields_in),     intent(in) :: Statein
!rab    type(state_fields_out),    intent(in) :: Stateout
    type (block_control_type), intent(in) :: Atm_block
    type(IPD_data_type),       intent(in) :: IPD_Data(:)
    integer,                   intent(in) :: nx, ny, levs, ntcw, ntoz
    real(kind=kind_phys),      intent(in) :: dt
    real(kind=kind_phys),      intent(in) :: time_int
    real(kind=kind_phys),      intent(in) :: fhswr, fhlwr
!--- local variables
    integer :: i, j, k, idx, nblks, nb, ix, ii, jj, kflip
    integer :: is_in, js_in, isc, jsc
    character(len=2) :: xtra
    real(kind=kind_phys), dimension(nx*ny) :: var2p
    real(kind=kind_phys), dimension(nx*ny,levs) :: var3p
    real(kind=kind_phys), dimension(nx,ny) :: var2, area, lat, lon, one, landmask, seamask
    real(kind=kind_phys), dimension(nx,ny,levs) :: var3
    real(kind=kind_phys) :: rdt, rtime_int, lcnvfac
    logical :: used

     nblks = Atm_block%nblks
     rdt = 1.0d0/dt
     rtime_int = 1.0d0/time_int

     isc = Atm_block%isc
     jsc = Atm_block%jsc
     is_in = Atm_block%isc
     js_in = Atm_block%jsc

     !Metrics
     do j = 1, ny
        jj = j + jsc -1
        do i = 1, nx
           ii = i + isc -1
           nb = Atm_block%blkno(ii,jj)
           ix = Atm_block%ixp(ii,jj)
           area(i,j) = IPD_Data(nb)%Grid%area(ix)
           lat(i,j)  = IPD_Data(nb)%Grid%xlat(ix)
           lon(i,j)  = IPD_Data(nb)%Grid%xlon(ix)
           one(i,j)  = 1.
           landmask(i,j) = IPD_Data(nb)%Sfcprop%slmsk(ix)
           seamask(i,j)  = 1. - landmask(i,j) 
        enddo
     enddo
     do idx = 1,tot_diag_idx
       if (Diag(idx)%id > 0) then
         lcnvfac = Diag(idx)%cnvfac
         if (trim(Diag(idx)%name) == 'DLWRFsfc' .or. trim(Diag(idx)%name) == 'ULWRFsfc' .or. &
             trim(Diag(idx)%name) == 'ULWRFtoa') then
           if (Diag(idx)%time_avg) lcnvfac = lcnvfac*min(rtime_int,1.0d0/fhlwr)
         elseif (trim(Diag(idx)%name) == 'DSWRFsfc' .or. trim(Diag(idx)%name) == 'USWRFsfc' .or. &
                 trim(Diag(idx)%name) == 'DSWRFtoa' .or. trim(Diag(idx)%name) == 'USWRFtoa') then
           if (Diag(idx)%time_avg) lcnvfac = lcnvfac*min(rtime_int,1.0d0/fhswr)
         elseif (trim(Diag(idx)%name) == 'TCDCclm' .or. trim(Diag(idx)%name) == 'TCDChcl' .or. &
                 trim(Diag(idx)%name) == 'TCDClcl' .or. trim(Diag(idx)%name) == 'TCDCmcl') then
           if (Diag(idx)%time_avg) lcnvfac = lcnvfac*min(rtime_int,max(1.0d0/fhswr,1.0d0/fhlwr))
         else
           if (Diag(idx)%time_avg) lcnvfac = lcnvfac*rtime_int
         endif
         if (Diag(idx)%axes == 2) then
           if (trim(Diag(idx)%name) == 'ALBDOsfc') then
             !--- albedos are actually a ratio of two radiation surface properties
             var2(1:nx,1:ny) = 0._kind_phys
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix) > 0._kind_phys) &
                   var2(i,j) = max(0._kind_phys,Diag(idx)%data(nb)%var2(ix)/Diag(idx)%data(nb)%var21(ix))*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%name) == 'gflux') then
             !--- need to "mask" gflux to output valid data over land/ice only
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix) /= 0) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%name) == 'soilm') then
             !--- need to "mask" soilm to have value only over land
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix) == 1) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           else
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           endif
!rab           used=send_data(Diag(idx)%id, var2, Time, is_in=is_in, js_in=js_in)
           used=send_data(Diag(idx)%id, var2, Time)
           !!!! Accumulated diagnostics --- lmh 19 sep 17
           if (fprint) then
           select case (trim(Diag(idx)%name))
           case('totprcp')
              call prt_gb_nh_sh_us('Total Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 86400.)
              call prt_gb_nh_sh_us('Land Precip  (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 86400.)
           case('totsnw')
              call prt_gb_nh_sh_us('Total Snowfall (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 777600.)
              call prt_gb_nh_sh_us('Land Snowfall  (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 777600.)
!           case('totgrp') ! Tiny??
!              call prt_gb_nh_sh_us('Total Icefall (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 172800.)
!              call prt_gb_nh_sh_us('Land Icefall  (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 172800.)
           case('dqsfc')
              call prt_gb_nh_sh_us('Total sfc LH flux  ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('dtsfc')
              call prt_gb_nh_sh_us('Total sfc SH flux  ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('DSWRFtoa')
              call prt_gb_nh_sh_us('TOA SW down ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('USWRFtoa')
              call prt_gb_nh_sh_us('TOA SW up ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('ULWRFtoa')
              call prt_gb_nh_sh_us('TOA LW up ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('t2m')
              call prt_gb_nh_sh_us('2-m T max ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MAX')
              call prt_gb_nh_sh_us('2-m T min ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MIN')
           case('tsfc')
              call prt_gb_nh_sh_us('sfc T max ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MAX')
              call prt_gb_nh_sh_us('sfc T min ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MIN')
              call prt_gb_nh_sh_us('SST max ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1., 'MAX')
              call prt_gb_nh_sh_us('SST min ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1., 'MIN')
           end select
           endif
         elseif (Diag(idx)%axes == 3) then
           !--- dt3dt variables ---- restored 16 feb 18 lmh
            do k=1,levs
            kflip=levs+1-k
            do j=1,ny
            jj = j + jsc -1
            do i=1,nx
               ii = i + isc - 1
               nb = Atm_block%blkno(ii,jj)
               ix = Atm_block%ixp(ii,jj)
               var3(i,j,k) = Diag(idx)%data(nb)%var3(ix,kflip)*lcnvfac
            enddo
            enddo
            enddo
            used=send_data(Diag(idx)%id, var3, Time)
!!$               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dt3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
!!$               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
#ifdef JUNK
           !--- dq3dt variables
           do num = 1,5+Mdl_parms%pl_coeff
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'dq3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dq3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
           enddo
           !--- du3dt and dv3dt variables
           do num = 1,4
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'du3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%du3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
             if (trim(Diag(idx)%name) == 'dv3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dv3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
           enddo
           if (trim(Diag(idx)%name) == 'dqdt_v') then
             var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dqdt_v(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- temperature tendency
           if (trim(Diag(idx)%name) == 'dtemp_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%tgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gt0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- horizontal wind component tendency
           if (trim(Diag(idx)%name) == 'du_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%ugrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gu0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- meridional wind component tendency
           if (trim(Diag(idx)%name) == 'dv_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%vgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gv0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- specific humidity tendency
           if (trim(Diag(idx)%name) == 'dsphum_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- cloud water mixing ration tendency
           if (trim(Diag(idx)%name) == 'dclwmr_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
           !--- ozone mixing ration tendency
           if (trim(Diag(idx)%name) == 'do3mr_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
          endif
#endif
        endif
     endif
    enddo


  end subroutine gfdl_diag_output
!-------------------------------------------------------------------------      
 subroutine prt_gb_nh_sh_us(qname, is,ie, js,je, a2, area, lon, lat, mask, fac, operation_in) !Prints averages/sums, or maxes/mins
  use physcons,    pi=>con_pi
  character(len=*), intent(in)::  qname
  integer, intent(in):: is, ie, js, je
  real(kind=kind_phys), intent(in), dimension(is:ie, js:je):: a2
  real(kind=kind_phys), intent(in), dimension(is:ie, js:je):: area, lon, lat, mask
  real, intent(in) :: fac
  character(len=*), intent(in), OPTIONAL :: operation_in
! Local:
  real(kind=kind_phys), parameter:: rad2deg = 180./pi
  real(kind=kind_phys) :: slat, slon
  real(kind=kind_phys):: t_eq, t_nh, t_sh, t_gb, t_us
  real(kind=kind_phys):: area_eq, area_nh, area_sh, area_gb, area_us
  integer:: i,j
  character(len=100) :: diagstr
  character(len=20)  :: diagstr1
  character(len=3)  :: operation

  if (present(operation_in)) then
     operation = operation_in(1:3)
  else
     operation = 'SUM'
  endif

     if (operation == "MAX") then
        t_eq =-1.e14   ;    t_nh =-1.e14;    t_sh =-1.e14;    t_gb =-1.e14;    t_us =-1.e14
        area_eq = 0.   ; area_nh = 0.   ; area_sh = 0.   ; area_gb = 0.   ; area_us = 0.
        do j=js,je
        do i=is,ie
           if (mask(i,j) <= 1.e-6) cycle

           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = 1.
           t_gb = max(t_gb,a2(i,j))
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = 1.
                t_eq = max(t_eq,a2(i,j))
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = 1.
                t_nh = max(t_nh,a2(i,j))
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = 1.
                t_sh = max(t_sh,a2(i,j))
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = 1.
              t_us = max(t_us,a2(i,j))
           endif
        enddo
        enddo

        call mp_reduce_max(   t_gb)
        call mp_reduce_max(   t_nh)
        call mp_reduce_max(   t_sh)
        call mp_reduce_max(   t_eq)
        call mp_reduce_max(   t_us)
     elseif (operation == "MIN") then
        t_eq = 1.e14   ;    t_nh = 1.e14;    t_sh = 1.e14;    t_gb = 1.e14;    t_us = 1.e14
        area_eq = 0.   ; area_nh = 0.   ; area_sh = 0.   ; area_gb = 0.   ; area_us = 0.
        do j=js,je
        do i=is,ie
           if (mask(i,j) <= 1.e-6) cycle

           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = 1.
           t_gb = min(t_gb,a2(i,j))
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = 1.
                t_eq = min(t_eq,a2(i,j))
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = 1.
                t_nh = min(t_nh,a2(i,j))
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = 1.
                t_sh = min(t_sh,a2(i,j))
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = 1.
              t_us = min(t_us,a2(i,j))
           endif
        enddo
        enddo

        call mp_reduce_min(   t_gb)
        call mp_reduce_min(   t_nh)
        call mp_reduce_min(   t_sh)
        call mp_reduce_min(   t_eq)
        call mp_reduce_min(   t_us)
     else
        t_eq = 0.   ;    t_nh = 0.;    t_sh = 0.;    t_gb = 0.;    t_us = 0.
        area_eq = 0.; area_nh = 0.; area_sh = 0.; area_gb = 0.; area_us = 0.
        operation = 'SUM'
        do j=js,je
        do i=is,ie
           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = area_gb + area(i,j)*mask(i,j)
           t_gb = t_gb + a2(i,j)*area(i,j)*mask(i,j)
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = area_eq + area(i,j)*mask(i,j)
                t_eq = t_eq + a2(i,j)*area(i,j)*mask(i,j)
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = area_nh + area(i,j)*mask(i,j)
                t_nh = t_nh + a2(i,j)*area(i,j)*mask(i,j)
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = area_sh + area(i,j)*mask(i,j)
                t_sh = t_sh + a2(i,j)*area(i,j)*mask(i,j)
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = area_us + area(i,j)*mask(i,j)
              t_us = t_us + a2(i,j)*area(i,j)*mask(i,j)
           endif
        enddo
        enddo

        call mp_reduce_sum(area_gb)
        call mp_reduce_sum(   t_gb)
        call mp_reduce_sum(area_nh)
        call mp_reduce_sum(   t_nh)
        call mp_reduce_sum(area_sh)
        call mp_reduce_sum(   t_sh)
        call mp_reduce_sum(area_eq)
        call mp_reduce_sum(   t_eq)
        call mp_reduce_sum(area_us)
        call mp_reduce_sum(   t_us)
     endif

     diagstr = trim(qname) // ' ' // trim(mpp_get_current_pelist_name()) // ' '
     !if (area_gb < 1.) then
     !   diagstr1 = ''
     !elseif( area_gb <= 4.*pi*RADIUS*RADIUS*.98) then
     !   write(diagstr1,101) 'Grid', t_gb/area_gb*fac
     !else
     write(diagstr1,101) 'GB', t_gb/area_gb*fac
     !endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_nh <= 1.e-6 .or. area_nh == area_gb) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'NH', t_nh/area_nh*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_sh <= 1.e-6 .or. area_sh == area_gb) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'SH', t_sh/area_sh*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_eq <= 1.e-6) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'EQ', t_eq/area_eq*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_us <= 1.e-6) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'US', t_us/area_us*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)

     if (is_master()) write(*,'(A)') trim(diagstr)

101  format(3x, A, ': ', F7.2)

   end subroutine prt_gb_nh_sh_us


end module FV3GFS_io_mod



