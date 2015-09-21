module gfs_physics_driver_mod

!-----------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!-----------------------------------------------------------------------

!--- FMS/GFDL modules ---
  use block_control_mod,  only: block_control_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use mpp_mod,            only: input_nml_file, mpp_pe, mpp_root_pe, &
                                mpp_error, mpp_chksum
  use field_manager_mod,  only: MODEL_ATMOS
  use fms_mod,            only: fms_init, stdout, stdlog, string,     &
                                mpp_clock_id, mpp_clock_begin,        &
                                mpp_clock_end, CLOCK_MODULE_DRIVER,   &
                                open_namelist_file, check_nml_error,  &
                                file_exist, open_file, close_file,    &
                                error_mesg, FATAL, WARNING, NOTE,     &
                                write_version_number
  use fms_io_mod,         only: restart_file_type, register_restart_field, &
                                restore_state, save_restart, &
                                get_mosaic_tile_file, read_data
  use time_manager_mod,   only: time_type, get_date, get_time, operator(-)
  use tracer_manager_mod, only: get_number_tracers

  use machine,            only: kind_phys
!--- NUOPC GFS Physics module routines ---
  use nuopc_physics,      only: nuopc_phys_init, nuopc_phys_run, &
                                nuopc_rad_run, nuopc_rad_update


!--- NUOPC GFS Physics module datatypes ---
  use nuopc_physics,      only: state_fields_in, state_fields_out,      &
                                model_parameters, dynamic_parameters,  &
                                sfc_properties, cloud_properties,       &
                                radiation_tendencies, interface_fields, &
                                diagnostics, tbd_ddt
!--- GFS Physics share module ---
  use module_CONSTANTS,   only: pi
  use physparam,          only: ipsd0
  use mersenne_twister,   only: random_setseed, random_index, random_stat
  use physcons,           only: dxmax, dxmin, dxinv


!-----------------------------------------------------------------------
  implicit none
  private

!--- public interfaces ---
  public  phys_rad_driver_init, phys_rad_setup_step, radiation_driver, &
          physics_driver, phys_rad_driver_end

!--- public NUOPC GFS datatypes and data typing ---
  public  state_fields_in, state_fields_out, kind_phys


!--- module private data ---
!--- NUOPC data types
  type(model_parameters) :: Mdl_parms
  type(tbd_ddt),              dimension(:), allocatable :: Tbd_data
  type(dynamic_parameters),   dimension(:), allocatable :: Dyn_parms
  type(diagnostics),          dimension(:), allocatable :: Gfs_diags
  type(sfc_properties),       dimension(:), allocatable :: Sfc_props
  type(cloud_properties),     dimension(:), allocatable :: Cld_props
  type(radiation_tendencies), dimension(:), allocatable :: Rad_tends
  type(interface_fields),     dimension(:), allocatable :: Intr_flds


!--- netcdf restart
  type(restart_file_type), pointer, save :: Phy_restart => NULL()
  type(restart_file_type), pointer, save :: Til_restart => NULL()


!--- mpp clocks ids
  integer :: gbphys_pre_clk, gbphys_clk, grrad_pre_clk, grrad_clk, diags_clk


!--- diagnostic field ids and var-names
  type gfdl_diag_type
    private
    integer :: id
    integer :: axes
    character(len=64)  :: name
    character(len=128) :: desc
    character(len=64)  :: unit
   end type gfdl_diag_type

   type(gfdl_diag_type), dimension(:), allocatable :: Diag

!--- miscellaneous other variables
  logical :: module_is_initialized = .FALSE.

!--- NAMELIST/CONFIGURATION parameters
!--- convenience variables
    integer :: lonr, latr
    integer :: nsswr, nslwr   ! trigger aids for radiation

!--- gfs initialization variables
    integer :: ipt = 1 
    integer :: latgfs = 1
    real(kind=kind_phys) :: xkzm_m = 1.0
    real(kind=kind_phys) :: xkzm_h = 1.0
    real(kind=kind_phys) :: xkzm_s = 1.0
    real(kind=kind_phys) :: evpco = 2.0e-5
    real(kind=kind_phys) :: psautco(2) = (/6.0e-4,3.0e-4/)
    real(kind=kind_phys) :: prautco(2) = (/.0e-4, 1.0e-4/)
    real(kind=kind_phys) :: wminco(2) = (/0.0e-5,1.0e-5/)
    real(kind=kind_phys) :: clstp
    real(kind=kind_phys) :: sup = 1.1
    real(kind=kind_phys) :: fhswr = 3600.
    real(kind=kind_phys) :: fhlwr = 3600.
    logical :: lprnt = .false.   ! control flag for diagnostic print out (rad)
    logical :: lssav = .true.    ! logical flag for store 3-d cloud field

!--- namelist parameters ---
    integer :: NFXR = 39
    integer :: ntoz = 2
    integer :: ntcw = 3
    integer :: ncld = 1
    integer :: levs = 64
    integer :: levr = 64
    integer :: me           ! set by call to mpp_pe
    integer :: lsoil = 4
    integer :: lsm = 1      ! NOAH LSM
    integer :: nmtvr = 14
    integer :: nrcm  = 2    ! when using ras, will be computed
    integer :: levozp = 80  ! read from global_o3prdlos.f77
    integer :: jcap  = 0    ! should not matter it is used by spherical 
    integer :: num_p3d = 3  ! Ferrier:3  Zhao:4 
    integer :: num_p2d = 1  ! Ferrier:1  Zhao:3
    integer :: npdf3d = 0   ! Zhao & pdfcld=.T.:3  -  else:0
    integer :: pl_coeff = 4
    integer :: ncw(2) = (/200,25/)
    real (kind=kind_phys) :: flgmin(2) = (/0.150,0.200/)
    real (kind=kind_phys) :: crtrh(3) = (/0.85,0.85,0.85/)
    real (kind=kind_phys) :: cdmbgwd(2) = (/1.0,1.0/)
    real (kind=kind_phys) :: ccwf(2) = (/1.0,1.0/)
    real (kind=kind_phys) :: dlqf(2) = (/0.5,0.5/)
    real (kind=kind_phys) :: ctei_rm(2) = (/10.0,10.0/)
    real (kind=kind_phys) :: cgwf(2) = (/0.1,0.1/)
    real (kind=kind_phys) :: prslrd0 = 200.
    logical :: ras = .false.
    logical :: pre_rad = .false.
    logical :: ldiag3d = .false.
    logical :: lgocart = .false. 
    logical :: cplflx = .false.  
    logical :: lssav_cpl = .false.
    logical :: flipv = .false.
    logical :: old_monin = .false. 
    logical :: cnvgwd = .false.
    logical :: shal_cnv = .true.
    logical :: sashal = .true.
    logical :: newsas = .true.
    logical :: cal_pre = .false.
    logical :: mom4ice = .false.
    logical :: mstrat = .false.
    logical :: trans_trac = .true.
    integer :: nst_fcst = 0
    logical :: moist_adj = .false.  ! Must be true to turn on moist convective
    integer :: thermodyn_id = 0     ! idvm/10
    integer :: sfcpress_id  = 1     ! idvm-(idvm/10)*10
    logical :: gen_coord_hybrid = .false. ! in scrpt, could be T or F
    logical :: lsidea = .false.  
    logical :: pdfcld = .false.
    logical :: shcnvcw = .false.
    logical :: redrag = .false.
    logical :: hybedmf = .false.
    logical :: dspheat = .false.

! Radiation option control parameters
    integer :: ictm = 0 
    integer :: isol = 0
    integer :: ico2 = 0 
    integer :: iaer = 0 
    integer :: ialb = 0 
    integer :: iems = 0 
    integer :: iovr_sw = 1
    integer :: iovr_lw = 1
    integer :: isubc_sw = 2
    integer :: isubc_lw = 2
    logical :: crick_proof = .false.
    logical :: ccnorm = .false.
    logical :: norad_precip = .false.  ! This is effective only for Ferrier/Moorthi

! rad_save
    integer :: iflip = 1         ! surface to toa

! interface props
    logical :: SW0 = .false.
    logical :: SWB = .false.
    logical :: LW0 = .false.
    logical :: LWB = .false.

!--- namelist ---
   namelist / gfs_physics_nml / norad_precip
!-----------------------------------------------------------------------

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
!--- phys_rad_driver_init ---
!    constructor for gfs_physics_driver_mod
!---------------------------------------------------------------------
  subroutine phys_rad_driver_init (Time, lon, lat, glon, glat, npz,       &
                                   axes, dx, dy, area, indxmin, indxmax,  &
                                   dt_phys, Atm_block, State_in, State_out)
    type(time_type),           intent(in) :: Time
    real,    dimension(:,:),   intent(in) :: lon, lat, dx, dy, area
    integer,                   intent(in) :: glon, glat, npz
    integer, dimension(4),     intent(in) :: axes
    real (kind=kind_phys), intent(in) :: indxmin, indxmax, dt_phys
    type (block_control_type), intent(in) :: Atm_block
    type (state_fields_in),    dimension(:), intent(inout) :: State_in
    type (state_fields_out),   dimension(:), intent(inout) :: State_out
!--- local variables
    integer ::  ierr, io, unit, logunit, outunit
    integer :: nb, ibs, ibe, jbs, jbe, ngptc
    integer :: i, j, ix, ntrac, ntp
    integer :: jdate(8) = (/1, 1, 1, 0, 0, 0, 0, 0/)
    integer :: idate(4) = (/0, 1, 1, 1/)
    integer :: kdt = 0
    integer :: nnp = 0
    real(kind=kind_phys) :: solhr = 0.0   
    real(kind=kind_phys) :: fhour = 0.
    real (kind=kind_phys) :: dxmaxin, dxminin, dxinvin
    logical :: sas_shal
    real (kind=kind_phys) :: si(npz+1)

!--- if routine has already been executed, return ---
    if (module_is_initialized) return

!--- verify that the modules used by this module that are called ---
!--- later in this subroutine have already been initialized ---
    call fms_init

    call get_date (Time, jdate(1), jdate(2), jdate(3),  &
                         jdate(5), jdate(6), jdate(7))
    idate(4) = jdate(1)
    idate(3) = jdate(2)
    idate(2) = jdate(3)
    idate(1) = jdate(5)

!--- read namelist  ---
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=gfs_physics_nml, iostat=io)
    ierr = check_nml_error(io,"gfs_physics_nml")
#else
    if ( file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=gfs_physics_nml, iostat=io, end=10)
      ierr = check_nml_error(io, 'gfs_physics_nml')
      enddo
 10   call close_file (unit)
    endif
#endif

!--- write version number and namelist to log file ---
    call write_version_number ('vers 1', 'gfs_physics_driver_mod')
    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe() ) write(logunit, nml=gfs_physics_nml)

!--- set some configurational parameters that are derived from others
    lonr = glon
    latr = glat
    nsswr = nint(fhswr/dt_phys)
    nslwr = nint(fhlwr/dt_phys)
    sas_shal = (sashal .and. (.not. ras))
 
!--- define the model dimensions on the local processor ---
    call get_number_tracers (MODEL_ATMOS, num_tracers=ntrac, num_prog=ntp)

!--- initialize the clocks ---
    grrad_pre_clk  = mpp_clock_id( '   GFS_Rad: pre        ', grain=CLOCK_MODULE_DRIVER )
    grrad_clk      = mpp_clock_id( '   GFS_Rad: grrad      ', grain=CLOCK_MODULE_DRIVER )
    gbphys_pre_clk = mpp_clock_id( '   GFS_Physics: pre    ', grain=CLOCK_MODULE_DRIVER )
    gbphys_clk     = mpp_clock_id( '   GFS_Physics: gbphys ', grain=CLOCK_MODULE_DRIVER )
    diags_clk      = mpp_clock_id( '   GFS: diagnostics    ', grain=CLOCK_MODULE_DRIVER )

!--- initialize physics ---
    unit = open_namelist_file ()

    me = mpp_pe()
    dxmaxin = indxmax
    dxminin = indxmin
    dxinvin = 1.0/(dxmaxin-dxminin)
    call nuopc_phys_init (Mdl_parms, ntcw, ncld, ntoz, ntrac, npz, me, lsoil, lsm, nmtvr, nrcm, levozp,  &
                          glon, glat, jcap, num_p3d, num_p2d, npdf3d, pl_coeff, ncw, crtrh, cdmbgwd,  &
                          ccwf, dlqf, ctei_rm, cgwf, prslrd0, ras, pre_rad, ldiag3d, lgocart,  &
                          lssav_cpl, flipv, old_monin, cnvgwd, shal_cnv, sashal, newsas, cal_pre, mom4ice,  &
                          mstrat, trans_trac, nst_fcst, moist_adj, thermodyn_id, sfcpress_id,  &
                          gen_coord_hybrid, npz, lsidea, pdfcld, shcnvcw, redrag, hybedmf, dspheat, &
                          dxmaxin, dxminin, dxinvin, &
                          ! For radiation
                          si, ictm, isol, ico2, iaer, ialb, iems,                    &
                          iovr_sw,iovr_lw,isubc_sw,isubc_lw,   &
                          sas_shal,crick_proof,ccnorm,norad_precip,jdate,iflip,dt_phys,unit)
    call close_file (unit)

!--- allocate and call the different storage items needed by GFS physics/radiation ---
    allocate( Tbd_data(Atm_block%nblks))
    allocate(Dyn_parms(Atm_block%nblks))
    allocate(Gfs_diags(Atm_block%nblks))
    allocate(Sfc_props(Atm_block%nblks))
    allocate(Cld_props(Atm_block%nblks))
    allocate(Rad_tends(Atm_block%nblks))
    allocate(Intr_flds(Atm_block%nblks))
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
      ngptc = (ibe - ibs + 1) * (jbe - jbs + 1) 

      call Tbd_data(nb)%set      (ngptc, Mdl_parms, xkzm_m, xkzm_h, xkzm_s, &
                                  evpco, psautco, prautco, wminco)
      call Dyn_parms(nb)%setrad  (ngptc, ngptc, kdt, jdate, solhr, fhlwr, fhswr, &
                                  lssav, ipt, lprnt, dt_phys)
      call Dyn_parms(nb)%setphys (ngptc, ngptc, solhr, kdt, lssav, latgfs, &
                                  dt_phys, dt_phys, clstp, nnp, fhour)
      call Gfs_diags(nb)%setrad  (ngptc, NFXR)
      call Gfs_diags(nb)%setphys (ngptc, Mdl_parms)
      call Sfc_props(nb)%setrad  (ngptc, Mdl_parms, .FALSE.)  ! last argument determines gsm vs atmos-only
      call Sfc_props(nb)%setphys (ngptc, Mdl_parms)
      call Cld_props(nb)%setrad  (ngptc, Mdl_parms, sup)
      call Cld_props(nb)%setphys (ngptc, Mdl_parms, sup)
      call Rad_tends(nb)%set     (ngptc, Mdl_parms)
      call Intr_flds(nb)%setrad  (ngptc, Mdl_parms, SW0, SWB, LW0, LWB)
      call Intr_flds(nb)%setphys (ngptc, Mdl_parms)
      call State_in(nb)%setrad   (ngptc, Mdl_parms)
      call State_in(nb)%setphys  (ngptc, Mdl_parms)
      call State_out(nb)%setphys (ngptc, Mdl_parms)

      !  populate static values that are grid dependent
      ix = 0 
      do j=1,jbe-jbs+1
       do i=1,ibe-ibs+1
        ix = ix + 1
        Dyn_parms(nb)%area(ix) = area(i,j)
        Dyn_parms(nb)%dx(ix)   = 0.5*(dx(i,j)+dx(i,j+1))
        Dyn_parms(nb)%dy(ix)   = 0.5*(dy(i,j)+dy(i+1,j))
        Dyn_parms(nb)%xlat(ix) = lat(i,j)*pi/180.0_kind_phys
        Dyn_parms(nb)%xlon(ix) = lon(i,j)*pi/180.0_kind_phys
        Dyn_parms(nb)%sinlat(ix) = sin(Dyn_parms(nb)%xlat(ix))
        Dyn_parms(nb)%coslat(ix) = sqrt(1.0_kind_phys - Dyn_parms(nb)%sinlat(ix)*Dyn_parms(nb)%sinlat(ix))
       enddo
      enddo
    enddo

!--- read in surface data from chgres ---
    call surface_props_input (Atm_block)

!--- initialize diagnostics ---
    call gfs_diag_register(Time, Atm_block, axes, NFXR)

!--- mark the module as initialized ---
      module_is_initialized = .true.

  end subroutine phys_rad_driver_init
!-----------------------------------------------------------------------



!-------------------------------------------------------------------------      
!--- phys_rad_setup_step ---
!-------------------------------------------------------------------------      
  subroutine phys_rad_setup_step (Time_init, Time, Time_next, Atm_block)
    type(time_type),            intent(in) :: Time_init, Time, Time_next
    type (block_control_type),  intent(in) :: Atm_block
!   local variables
    integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds
    integer :: i, j, k, nb, ix
    integer :: ibs, ibe, jbs, jbe, nx, ny
    integer :: sec, ipseed, jdate(8)
    integer :: kdt
    integer :: numrdm(lonr*latr*2)
    type (random_stat) :: stat
    real(kind=kind_phys) :: fhour
    real(kind=kind_phys) :: work1, work2

!--- set the date
    call get_date (Time, jdate(1), jdate(2), jdate(3),  &
                         jdate(5), jdate(6), jdate(7))

    call get_time(Time - Time_init, sec)
    fhour = real(sec)/3600.

!--- may need this to repopulate sfc properties for AMIP runs
!    call sfc_populate (Sfc_props)


!--- set up random seed index for the whole tile in a reproducible way
    if (isubc_lw==2 .or. isubc_sw==2) then
      ipseed = mod(nint(100.0*sqrt(real(sec))), ipsdlim) + 1 + ipsd0
      call random_setseed (ipseed, stat)
      call random_index (ipsdlim, numrdm, stat)
    endif

    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
      nx = ibe-ibs+1
      ny = jbe-jbs+1

!--- increment the time step number
      Dyn_parms(nb)%kdt      = Dyn_parms(nb)%kdt + 1
      Dyn_parms(nb)%nnp      = Dyn_parms(nb)%nnp + 1
!--- set the current forecast hour
      Dyn_parms(nb)%fhour    = fhour
      Dyn_parms(nb)%jdate(:) = jdate(:)
!--- set the solhr
      Dyn_parms(nb)%solhr    = real(jdate(5))
!--- radiation triggers
      Dyn_parms(nb)%lsswr = (mod(Dyn_parms(nb)%kdt, nsswr) == 1)
      Dyn_parms(nb)%lslwr = (mod(Dyn_parms(nb)%kdt, nslwr) == 1)
! **************  Ken Campana Stuff  ********************************
!...  set switch for saving convective clouds
      if(Dyn_parms(nb)%lsswr) then
        Dyn_parms(nb)%clstp = 1100                           !initialize,accumulate
      else
        Dyn_parms(nb)%clstp = 0100                           !accumulate
      endif
! **************  Ken Campana Stuff  ********************************

!--- set the random seeds for each column
      if (isubc_lw==2 .or. isubc_sw==2) then
        ix = 0 
        do j = 1,ny
          do i = 1,nx
            ix = ix + 1
            Dyn_parms(nb)%icsdsw(i) = numrdm(i+ibs-1 + (j+jbs-2)*lonr)
            Dyn_parms(nb)%icsdlw(i) = numrdm(i+ibs-1 + (j+jbs-2)*lonr + latr)

            if (Mdl_parms%num_p3d == 3) then
              work1 = (Dyn_parms(nb)%dx(ix) - dxmin) * dxinv
              work1 = max(0.0, min(1.0,work1))
              Cld_props(nb)%flgmin(ix) = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
            else
              Cld_props(nb)%flgmin(ix) = 0.0
            endif
          enddo
        enddo
      endif
    enddo

  end subroutine phys_rad_setup_step


!-------------------------------------------------------------------------      
!--- radiation_driver ---
!-------------------------------------------------------------------------      
  subroutine radiation_driver (Time, Time_next, Atm_block, Statein)
    type(time_type),                         intent(in) :: Time, Time_next
    type (block_control_type),               intent(in) :: Atm_block
    type(state_fields_in),     dimension(:), intent(in) :: Statein
!   local variables
    integer :: nb

!--- call the nuopc radiation loop---
    call mpp_clock_begin(grrad_clk)
    do nb = 1, Atm_block%nblks

!
!--- call the nuopc radiation routine for time-varying data ---
!        can this be run inside parallel-do or does it need to be outside
!
      call nuopc_rad_update (Mdl_parms, Dyn_parms(nb))
      call nuopc_rad_run (Statein(nb), Sfc_props(nb), Gfs_diags(nb), &
                          Intr_flds(nb), Cld_props(nb), Rad_tends(nb), &
                          Mdl_parms, Dyn_parms(nb))
    enddo
    call mpp_clock_end(grrad_clk)

  end subroutine radiation_driver
!-------------------------------------------------------------------------      
 


!-------------------------------------------------------------------------      
!--- physics_driver ---
!-------------------------------------------------------------------------      
  subroutine physics_driver (Time, Time_next, Atm_block, Statein, Stateout)
    type(time_type),                         intent(in)    :: Time, Time_next
    type (block_control_type),               intent(in)    :: Atm_block
    type(state_fields_in),     dimension(:), intent(in)    :: Statein
    type(state_fields_out),    dimension(:), intent(inout) :: Stateout
!   local variables
    integer :: nb

!--- call the nuopc physics loop---
    call mpp_clock_begin(gbphys_clk)
    do nb = 1, Atm_block%nblks

      Tbd_data(nb)%dpshc(:) = 0.3d0 * Statein(nb)%prsi(:,1)

      call nuopc_phys_run (Statein(nb), Stateout(nb), Sfc_props(nb), &
                           Gfs_diags(nb), Intr_flds(nb), Cld_props(nb), &
                           Rad_tends(nb), Mdl_parms, Tbd_data(nb), &
                           Dyn_parms(nb))
    enddo
    call mpp_clock_end(gbphys_clk)

    call mpp_clock_begin(diags_clk)
    call gfs_diag_output (Time, Atm_block)
    call mpp_clock_end(diags_clk)

  end subroutine physics_driver
!-------------------------------------------------------------------------      



!-------------------------------------------------------------------------      
!--- phys_rad_driver_end ---
!-------------------------------------------------------------------------      
  subroutine phys_rad_driver_end (Time, Atm_block)
    type(time_type),            intent(in) :: Time
    type (block_control_type),  intent(in) :: Atm_block

!--- need to figure this one out yet

  end subroutine phys_rad_driver_end
!-------------------------------------------------------------------------      



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------      
!--- surface_props_input ---
!-------------------------------------------------------------------------      
  subroutine surface_props_input (Atm_block, GSM)
    type (block_control_type), intent(in) :: Atm_block
    logical, intent(in), optional :: GSM
!   local variables
    integer :: i, j, ibs, ibe, jbs, jbe, nct
    integer :: nb, nx, ny, ngptc
    integer :: start(4), nread(4)
    character(len=32)  :: fn_srf = 'INPUT/sfc_data.nc'
    character(len=32)  :: fn_oro = 'INPUT/oro_data.nc'
    character(len=128) :: errmsg
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3 => NULL()
    logical :: exists

    call get_mosaic_tile_file (fn_srf, fn_srf, .FALSE.)
    call get_mosaic_tile_file (fn_oro, fn_oro, .FALSE.)

    inquire(file=trim(fn_srf), exist=exists)
    if (exists) then
      errmsg = 'opening file '//trim(fn_srf)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), NOTE)
    else
      errmsg = 'error opening file '//trim(fn_srf)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), FATAL)
    endif

    inquire(file=trim(fn_oro), exist=exists)
    if (exists) then
      errmsg = 'opening file '//trim(fn_oro)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), NOTE)
    else
      errmsg = 'error opening file '//trim(fn_oro)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), FATAL)
    endif

    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
!rab      ibs = Atm_block%ibs(nb)-Atm_block%isc+1
!rab      ibe = Atm_block%ibe(nb)-Atm_block%isc+1
!rab      jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
!rab      jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
      nx = (ibe - ibs + 1)
      ny = (jbe - jbs + 1) 
      ngptc = nx * ny
      start(1) = ibs
      start(2) = jbs
      start(3) = 1
      start(4) = 1
      nread(1) = nx
      nread(2) = ny
      nread(3) = 1
      nread(4) = 1

!--- OROGRAPHY FILE
!--- stddev
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,1)
      call read_data(fn_oro,'stddev',var2,start,nread)
      Sfc_props(nb)%hprim(1:ngptc) = Sfc_props(nb)%hprime2(1:ngptc,1)
!--- convexity
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,2)
      call read_data(fn_oro,'convexity',var2,start,nread)
!--- oa1
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,3)
      call read_data(fn_oro,'oa1',var2,start,nread)
!--- oa2
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,4)
      call read_data(fn_oro,'oa2',var2,start,nread)
!--- oa3
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,5)
      call read_data(fn_oro,'oa3',var2,start,nread)
!--- oa4
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,6)
      call read_data(fn_oro,'oa4',var2,start,nread)
!--- ol1
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,7)
      call read_data(fn_oro,'ol1',var2,start,nread)
!--- ol2
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,8)
      call read_data(fn_oro,'ol2',var2,start,nread)
!--- ol3
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,9)
      call read_data(fn_oro,'ol3',var2,start,nread)
!--- ol4
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,10)
      call read_data(fn_oro,'ol4',var2,start,nread)
!--- theta
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,11)
      call read_data(fn_oro,'theta',var2,start,nread)
!--- gamma
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,12)
      call read_data(fn_oro,'gamma',var2,start,nread)
!--- sigma
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,13)
      call read_data(fn_oro,'sigma',var2,start,nread)
!--- elvmax
      var2(1:nx,1:ny) => Sfc_props(nb)%hprime2(1:ngptc,14)
      call read_data(fn_oro,'elvmax',var2,start,nread)

!--- SURFACE FILE
!--- slmsk
      var2(1:nx,1:ny) => Sfc_props(nb)%slmsk(1:ngptc)
      call read_data(fn_srf,'slmsk',var2,start,nread)
!--- oro (orog_filt in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%oro(1:ngptc)
      call read_data(fn_srf,'orog_filt',var2,start,nread)
!--- oro_uf (orog_raw in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%oro_uf(1:ngptc)
      call read_data(fn_srf,'orog_raw',var2,start,nread)
!--- tsfc (tsea in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%tsfc(1:ngptc)
      call read_data(fn_srf,'tsea',var2,start,nread)
!--- weasd (sheleg in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%weasd(1:ngptc)
      call read_data(fn_srf,'sheleg',var2,start,nread)
!--- tg3
      var2(1:nx,1:ny) => Sfc_props(nb)%tg3(1:ngptc)
      call read_data(fn_srf,'tg3',var2,start,nread)
!--- zorl
      var2(1:nx,1:ny) => Sfc_props(nb)%zorl(1:ngptc)
      call read_data(fn_srf,'zorl',var2,start,nread)
!rab      if (present(GSM)) then
!rab        if (GSM) then
!---     alvsf
          var2(1:nx,1:ny) => Sfc_props(nb)%alvsf(1:ngptc)
          call read_data(fn_srf,'alvsf',var2,start,nread)
!---     alvwf
          var2(1:nx,1:ny) => Sfc_props(nb)%alvwf(1:ngptc)
          call read_data(fn_srf,'alvwf',var2,start,nread)
!---     alnsf
          var2(1:nx,1:ny) => Sfc_props(nb)%alnsf(1:ngptc)
          call read_data(fn_srf,'alnsf',var2,start,nread)
!---     alnwf
          var2(1:nx,1:ny) => Sfc_props(nb)%alnwf(1:ngptc)
          call read_data(fn_srf,'alnwf',var2,start,nread)
!---     facsf
          var2(1:nx,1:ny) => Sfc_props(nb)%facsf(1:ngptc)
          call read_data(fn_srf,'facsf',var2,start,nread)
!---     facwf
          var2(1:nx,1:ny) => Sfc_props(nb)%facwf(1:ngptc)
          call read_data(fn_srf,'facwf',var2,start,nread)
!rab        endif
!rab      endif
!--- vfrac
      var2(1:nx,1:ny) => Sfc_props(nb)%vfrac(1:ngptc)
      call read_data(fn_srf,'vfrac',var2,start,nread)
!--- canopy
      var2(1:nx,1:ny) => Sfc_props(nb)%canopy(1:ngptc)
      call read_data(fn_srf,'canopy',var2,start,nread)
!--- f10m
      var2(1:nx,1:ny) => Sfc_props(nb)%f10m(1:ngptc)
      call read_data(fn_srf,'f10m',var2,start,nread)
!--- t2m
      var2(1:nx,1:ny) => Sfc_props(nb)%t2m(1:ngptc)
      call read_data(fn_srf,'t2m',var2,start,nread)
!--- q2m
      var2(1:nx,1:ny) => Sfc_props(nb)%q2m(1:ngptc)
      call read_data(fn_srf,'q2m',var2,start,nread)
!--- vtype
      var2(1:nx,1:ny) => Sfc_props(nb)%vtype(1:ngptc)
      call read_data(fn_srf,'vtype',var2,start,nread)
!--- stype
      var2(1:nx,1:ny) => Sfc_props(nb)%stype(1:ngptc)
      call read_data(fn_srf,'stype',var2,start,nread)
!--- uustar
      var2(1:nx,1:ny) => Sfc_props(nb)%uustar(1:ngptc)
      call read_data(fn_srf,'uustar',var2,start,nread)
!--- ffmm
      var2(1:nx,1:ny) => Sfc_props(nb)%ffmm(1:ngptc)
      call read_data(fn_srf,'ffmm',var2,start,nread)
!--- ffhh
      var2(1:nx,1:ny) => Sfc_props(nb)%ffhh(1:ngptc)
      call read_data(fn_srf,'ffhh',var2,start,nread)
!--- hice
      var2(1:nx,1:ny) => Sfc_props(nb)%hice(1:ngptc)
      call read_data(fn_srf,'hice',var2,start,nread)
!--- fice
      var2(1:nx,1:ny) => Sfc_props(nb)%fice(1:ngptc)
      call read_data(fn_srf,'fice',var2,start,nread)
!--- tisfc
      var2(1:nx,1:ny) => Sfc_props(nb)%tisfc(1:ngptc)
      call read_data(fn_srf,'tisfc',var2,start,nread)
!--- tprcp
      var2(1:nx,1:ny) => Tbd_data(nb)%tprcp(1:ngptc)
      call read_data(fn_srf,'tprcp',var2,start,nread)
!--- srflag
      var2(1:nx,1:ny) => Tbd_data(nb)%srflag(1:ngptc)
      call read_data(fn_srf,'srflag',var2,start,nread)
!--- snowd (snwdph in the file)
      var2(1:nx,1:ny) => Sfc_props(nb)%snowd(1:ngptc)
      call read_data(fn_srf,'snwdph',var2,start,nread)
!--- shdmin
      var2(1:nx,1:ny) => Sfc_props(nb)%shdmin(1:ngptc)
      call read_data(fn_srf,'shdmin',var2,start,nread)
!--- shdmax
      var2(1:nx,1:ny) => Sfc_props(nb)%shdmax(1:ngptc)
      call read_data(fn_srf,'shdmax',var2,start,nread)
!--- slope
      var2(1:nx,1:ny) => Sfc_props(nb)%slope(1:ngptc)
      call read_data(fn_srf,'slope',var2,start,nread)
!--- snoalb
      var2(1:nx,1:ny) => Sfc_props(nb)%snoalb(1:ngptc)
      call read_data(fn_srf,'snoalb',var2,start,nread)
!
!--- 3D variables
      allocate(var3(1:nx,1:ny,1:Mdl_parms%lsoil))
      start(1) = ibs
      start(2) = jbs
      start(3) = 1
      start(4) = 1
      nread(1) = nx
      nread(2) = ny
      nread(3) = Mdl_parms%lsoil
      nread(4) = 1
!--- stc
      call read_data(fn_srf,'stc',var3,start,nread)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%stc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
!--- smc
      call read_data(fn_srf,'smc',var3,start,nread)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%smc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
!--- slc
      call read_data(fn_srf,'slc',var3,start,nread)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%slc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
      deallocate(var3)
    enddo

!--- nullify/deallocate any temporaries used
    nullify(var2)

  end subroutine surface_props_input


!-------------------------------------------------------------------------      
!--- gfs_diag_register ---
!    10+NFXR - radiation
!    76+pl_coeff - physics
!-------------------------------------------------------------------------      
  subroutine gfs_diag_register(Time, Atm_block, axes, NFXR)
    type(time_type),           intent(in) :: Time
    type (block_control_type), intent(in) :: Atm_block
    integer, dimension(4),     intent(in) :: axes
    integer,                   intent(in) :: NFXR
!--- local variables
    integer :: idx, num, axes_l
    character(len=2) :: xtra
    character(len=11) :: mod_name = 'gfs_physics'

    allocate(Diag(86+NFXR+Mdl_parms%pl_coeff))

    Diag(:)%id = 0

    idx = 0

!--- accumulated diagnostics ---
    do num = 1,NFXR
      write (xtra,'(I2)') num 
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'fluxr_'//trim(xtra)
      Diag(idx)%desc = 'fluxr diagnostic '//trim(xtra)//' - GFS radiation'
      Diag(idx)%unit = 'XXX'
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dswcmp'
    Diag(idx)%desc = 'dswcmp dagnostic - GFS radiation'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'uswcmp'
    Diag(idx)%desc = 'uswcmp dagnostic - GFS radiation'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfxc'
    Diag(idx)%desc = 'total sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_dnflx'
    Diag(idx)%desc = 'total sky downward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfx0'
    Diag(idx)%desc = 'clear sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lw_upfxc'
    Diag(idx)%desc = 'total sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfx0'
    Diag(idx)%desc = 'clear sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'srunoff'
    Diag(idx)%desc = 'surface water runoff - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evbsa'
    Diag(idx)%desc = 'evbsa - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evcwa'
    Diag(idx)%desc = 'evcwa - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snohfa'
    Diag(idx)%desc = 'snohfa - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'transa'
    Diag(idx)%desc = 'transa - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sbsnoa'
    Diag(idx)%desc = 'sbsnoa - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snowca'
    Diag(idx)%desc = 'snowca - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'soilm'
    Diag(idx)%desc = 'soil moisture - GFS lsm'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmin'
    Diag(idx)%desc = 'min temperature at 2m height - GFS physics'
    Diag(idx)%unit = 'k'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmax'
    Diag(idx)%desc = 'max temperature at 2m height - GFS physics'
    Diag(idx)%unit = 'k'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmin'
    Diag(idx)%desc = 'min temperature at 2m height - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfc'
    Diag(idx)%desc = 'u component of surface stress - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfc'
    Diag(idx)%desc = 'v component of surface stress - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dtsfc'
    Diag(idx)%desc = 'sensible heat flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dqsfc'
    Diag(idx)%desc = 'latent heat flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totprcp'
    Diag(idx)%desc = 'accumulated total precipitation - GFS physics'
    Diag(idx)%unit = 'kg/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gflux'
    Diag(idx)%desc = 'ground conductive heat flux - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dlwsfc'
    Diag(idx)%desc = 'time accumulated downward lw flux at surface- GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ulwsfc'
    Diag(idx)%desc = 'time accumulated upward lw flux at surface- GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'suntim'
    Diag(idx)%desc = 'sunshine duration time - GFS physics'
    Diag(idx)%unit = 's'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'runoff'
    Diag(idx)%desc = 'total water runoff - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ep'
    Diag(idx)%desc = 'potential evaporation - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cldwrk'
    Diag(idx)%desc = 'cloud workfunction (valaxes only with sas) - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dugwd'
    Diag(idx)%desc = 'vertically integrated u change by OGWD - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvgwd'
    Diag(idx)%desc = 'vertically integrated v change by OGWD - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psmean'
    Diag(idx)%desc = 'surface pressure - GFS physics'
    Diag(idx)%unit = 'kPa'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cnvprcp'
    Diag(idx)%desc = 'accumulated convective precipitation - GFS physics'
    Diag(idx)%unit = 'kg/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmin'
    Diag(idx)%desc = 'minimum specific humidity - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmax'
    Diag(idx)%desc = 'maximum specific humidity - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rain'
    Diag(idx)%desc = 'total rain at this time step - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rainc'
    Diag(idx)%desc = 'convective rain at this time step - GFS physics'
    Diag(idx)%unit = 'XXX'

    do num = 1,6
      write (xtra,'(I2)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dt3dt_'//trim(xtra)
      Diag(idx)%desc = 'temperature change due to physics '//trim(xtra)//' - GFS physics'
      Diag(idx)%unit = 'XXX'
    enddo

    do num = 1,5+Mdl_parms%pl_coeff
      write (xtra,'(I2)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dq3dt_'//trim(xtra)
      Diag(idx)%desc = 'moisture change due to physics '//trim(xtra)//' - GFS physics'
      Diag(idx)%unit = 'XXX'
    enddo

    do num = 1,4
      write (xtra,'(I2)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'du3dt_'//trim(xtra)
      Diag(idx)%desc = 'u momentum change due to physics '//trim(xtra)//' - GFS physics'
      Diag(idx)%unit = 'XXX'
    enddo

    do num = 1,4
      write (xtra,'(I2)') num 
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dv3dt_'//trim(xtra)
      Diag(idx)%desc = 'v momentum change due to physics '//trim(xtra)//' - GFS physics'
      Diag(idx)%unit = 'XXX'
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dqdt_v'
    Diag(idx)%desc = 'total moisture tendency - GFS physics'
    Diag(idx)%unit = 'XXX'

!--- instantaneous diagnostics ---
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u10m'
    Diag(idx)%desc = '10 meter u windspeed - GFS physics'
    Diag(idx)%unit = 'm/s'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v10m'
    Diag(idx)%desc = '10 meter v windspeed - GFS physics'
    Diag(idx)%unit = 'm/s'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'zlvl'
    Diag(idx)%desc = 'layer 1 height - GFS physics'
    Diag(idx)%unit = 'm'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psurf'
    Diag(idx)%desc = 'surface pressure - GFS physics'
    Diag(idx)%unit = 'Pa'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hpbl'
    Diag(idx)%desc = 'pbl height - GFS physics'
    Diag(idx)%unit = 'm'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pwat'
    Diag(idx)%desc = 'precipitable water - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 't1'
    Diag(idx)%desc = 'layer 1 temperature - GFS physics'
    Diag(idx)%unit = 'K'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'q1'
    Diag(idx)%desc = 'layer 1 specific humidity - GFS physics'
    Diag(idx)%unit = 'kg/kg'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u1'
    Diag(idx)%desc = 'layer 1 zonal wind - GFS physics'
    Diag(idx)%unit = 'm/s'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v1'
    Diag(idx)%desc = 'layer 1 meridional wind - GFS physics'
    Diag(idx)%unit = 'm/s'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'chh'
    Diag(idx)%desc = 'thermal exchange coefficient - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cmm'
    Diag(idx)%desc = 'momentum exchange coefficient - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dlwsfci'
    Diag(idx)%desc = 'instantaneous sfc downward lw flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ulwsfci'
    Diag(idx)%desc = 'instantaneous sfc upward lw flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dswsfci'
    Diag(idx)%desc = 'instantaneous sfc downward sw flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'uswsfci'
    Diag(idx)%desc = 'instantaneous sfc upward sw flux - GFS physics'
    Diag(idx)%unit = 'w/m**2'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfci'
    Diag(idx)%desc = 'instantaneous u component of surface stress - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfci'
    Diag(idx)%desc = 'instantaneous v component of surface stress - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dtsfci'
    Diag(idx)%desc = 'instantaneous surface sensible heat flux - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dqsfci'
    Diag(idx)%desc = 'instantaneous surface latent heat flux - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gfluxi'
    Diag(idx)%desc = 'instantaneous surface ground heat flux - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'epi'
    Diag(idx)%desc = 'instantaneous surface potential evaporation - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'smcwlt2'
    Diag(idx)%desc = 'wiltimg point (volumetric) - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'smcref2'
    Diag(idx)%desc = 'soil moisture threshold (volumetric) - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wet1'
    Diag(idx)%desc = 'normalized soil wetness - GFS physics'
    Diag(idx)%unit = 'XXX'

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sr'
    Diag(idx)%desc = 'ratio of snow to total precipitation - GFS physics'
    Diag(idx)%unit = 'XXX'


    do num = 1,size(Diag,1)
      axes_l = Diag(num)%axes
      Diag(num)%id = register_diag_field (mod_name, trim(Diag(num)%name),  &
                                           axes(1:axes_l), Time, trim(Diag(num)%desc), &
                                           trim(Diag(num)%unit), missing_value=1.0d-30)
    enddo

  end subroutine gfs_diag_register
!-------------------------------------------------------------------------      


!-------------------------------------------------------------------------      
!--- gfs_diag_output ---
!-------------------------------------------------------------------------      
  subroutine gfs_diag_output(Time, Atm_block)
    type(time_type),           intent(in) :: Time
    type (block_control_type), intent(in) :: Atm_block
!--- local variables
    integer :: nb, ibs, ibe, jbs, jbe, nx, ny, ngptc
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3 => NULL()
    logical :: used

!--- call the nuopc physics loop---
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)-Atm_block%isc+1
      ibe = Atm_block%ibe(nb)-Atm_block%isc+1
      jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
      jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
      nx = (ibe - ibs + 1)
      ny = (jbe - jbs + 1) 
      ngptc = nx*ny

      var2(1:nx,1:ny) => Gfs_diags(nb)%srunoff(1:ngptc)
      used=send_data(Diag(1)%id, var2, Time, is_in=ibs, js_in=jbs)
    enddo

  end subroutine gfs_diag_output
!-------------------------------------------------------------------------      

end module gfs_physics_driver_mod
