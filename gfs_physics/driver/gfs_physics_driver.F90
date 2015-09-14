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
  use time_manager_mod,   only: time_type, get_time, operator (-), &
                                time_manager_init, operator(*)
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
!-----------------------------------------------------------------------
  implicit none
  private

!--- public interfaces ---
  public  phys_rad_driver_init, radiation_driver, physics_driver, &
          phys_rad_driver_end

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
  subroutine phys_rad_driver_init (Time, lon, lat, nlev, axes, Atm_block, &
                                   State_in, State_out)
    type(time_type),           intent(in) :: Time
    real,    dimension(:,:),   intent(in) :: lon, lat
    integer,                   intent(in) :: nlev
    integer, dimension(4),     intent(in) :: axes
    type (block_control_type), intent(in) :: Atm_block
    type (state_fields_in),    dimension(:), intent(inout) :: State_in
    type (state_fields_out),   dimension(:), intent(inout) :: State_out
!--- local variables
    integer :: nb, ibs, ibe, jbs, jbe, ngptc
    integer :: i, j, ix
    integer :: ntrac, ntp
    integer          ::  ierr, io, unit, logunit, outunit
!--- gfs initialization variables
    integer :: kdt, jdate(8), ipt, latgfs, nnp, NFXR
    real(kind=kind_phys) :: xkzm_m, xkzm_h, xkzm_s, evpco
    real(kind=kind_phys) :: psautco(2), prautco(2), wminco(2)
    real(kind=kind_phys) :: slag, sdec, cdec, sup, clstp
    real(kind=kind_phys) :: solhr, solcon, dtp, dtf, dtlw, dtsw, deltim, fhour
    logical :: lsswr, lslwr, lssav, lprnt, SW0, SWB, LW0, LWB

!--- namelist parameters ---
    integer :: ntcw
    integer :: ncld
    integer :: ntoz
!rab    integer :: ntrac
    integer :: levr
    integer :: levs
    integer :: me
    integer :: lsoil
    integer :: lsm = 1
    integer :: nmtvr = 14
    integer :: nrcm
    integer :: levozp
    integer :: lonr = 768
    integer :: latr = 768
    integer :: jcap = 0
    integer :: num_p3d 
    integer :: num_p2d
    integer :: npdf3d
    integer :: pl_coeff = 5
    integer :: ncw(2)
    real (kind=kind_phys) :: si(nlev+1) ! (levr+1)
    real (kind=kind_phys) :: crtrh(3)
    real (kind=kind_phys) :: cdmbgwd(2)
    real (kind=kind_phys) :: ccwf(2)
    real (kind=kind_phys) :: dlqf(2)
    real (kind=kind_phys) :: ctei_rm(2)
    real (kind=kind_phys) :: cgwf(2)
    real (kind=kind_phys) :: prslrd0
    real (kind=kind_phys) :: dxmaxin, dxminin, dxinvin
    logical :: ras 
    logical :: pre_rad
    logical :: ldiag3d = .true.
    logical :: lgocart = .false.
    logical :: cplflx
    logical :: lssav_cpl = .false.
    logical :: flipv = .false.
    logical :: old_monin = .false.
    logical :: cnvgwd
    logical :: shal_cnv
    logical :: sashal
    logical :: newsas
    logical :: cal_pre
    logical :: mom4ice
    logical :: mstrat
    logical :: trans_trac
    integer :: nst_fcst = 0
    logical :: moist_adj = .false.
    integer :: thermodyn_id
    integer :: sfcpress_id
    logical :: gen_coord_hybrid
    logical :: lsidea
    logical :: pdfcld
    logical :: shcnvcw
    logical :: redrag
    logical :: hybedmf
    logical :: dspheat

! Radiation option control parameters
    integer :: ictm, isol, ico2, iaer, ialb, iems
    integer :: iovr_sw, iovr_lw, isubc_sw, isubc_lw
    logical :: sas_shal, crick_proof, ccnorm, norad_precip
! rad_savei
    integer :: idate(4) = (/1, 1, 0, 0/)
    integer :: iflip = 0         ! toa to surface

    levr = nlev

!--- namelist ---
   namelist / gfs_physics_nml / norad_precip


!--- if routine has already been executed, return ---
    if (module_is_initialized) return

!--- verify that the modules used by this module that are called ---
!--- later in this subroutine have already been initialized ---
    call fms_init
 
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
    call nuopc_phys_init (Mdl_parms, ntcw, ncld, ntoz, ntrac, levs, me, lsoil, lsm, nmtvr, nrcm, levozp,  &
                          lonr, latr, jcap, num_p3d, num_p2d, npdf3d, pl_coeff, ncw, crtrh, cdmbgwd,  &
                          ccwf, dlqf, ctei_rm, cgwf, prslrd0, ras, pre_rad, ldiag3d, lgocart,  &
                          lssav_cpl, flipv, old_monin, cnvgwd, shal_cnv, sashal, newsas, cal_pre, mom4ice,  &
                          mstrat, trans_trac, nst_fcst, moist_adj, thermodyn_id, sfcpress_id,  &
                          gen_coord_hybrid, levr, lsidea, pdfcld, shcnvcw, redrag, hybedmf, dspheat, &
                          dxmaxin, dxminin, dxinvin, &
                          ! For radiation
                          si, ictm, isol, ico2, iaer, ialb, iems,                    &
                          iovr_sw,iovr_lw,isubc_sw,isubc_lw,   &
                          sas_shal,crick_proof,ccnorm,norad_precip,idate,iflip, unit)
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
      call Dyn_parms(nb)%setrad  (ngptc, ngptc, kdt, jdate, solhr, solcon,  &
                                  dtlw, dtsw, lsswr, lslwr, lssav, ipt, lprnt, &
                                  deltim, slag, sdec, cdec)
      call Dyn_parms(nb)%setphys (ngptc, ngptc, solhr, kdt, lssav, latgfs, &
                                  dtp, dtf, clstp, nnp, fhour, slag, &
                                  sdec, cdec)
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
      do j=jbs,jbe
       do i=ibs,ibe
        ix = ix + 1
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
!--- radiation_driver ---
!-------------------------------------------------------------------------      
  subroutine radiation_driver (Time, Time_next, Atm_block, Statein)
    type(time_type),                         intent(in) :: Time, Time_next
    type (block_control_type),               intent(in) :: Atm_block
    type(state_fields_in),     dimension(:), intent(in) :: Statein
!   local variables
    integer :: nb

!need to populate anything needing help
!    call sfc_populate (Sfc_props)

!rab      if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
!rab!       ipseed = mod(nint(100.0*sqrt(fhour*3600)), ipsdlim) + 1 + ipsd0
!rab        ipseed = mod(nint(100.0*sqrt(phour*3600)), ipsdlim) + 1 + ipsd0
!rab
!rab        call random_setseed                                             &
!rab!  ---  inputs:
!rab     &     ( ipseed,                                                    &
!rab!  ---  outputs:
!rab     &       stat                                                       &
!rab     &      )
!rab        call random_index                                               &
!rab!  ---  inputs:
!rab     &     ( ipsdlim,                                                   &
!rab!  ---  outputs:
!rab     &       numrdm, stat                                               &
!rab     &     )
!rab
!rab        do k = 1, 2
!rab          do j = 1, lats_node_r
!rab            lat = global_lats_r(ipt_lats_node_r-1+j)
!rab
!rab            do i = 1, LONR
!rab              ixseed(i,j,k) = numrdm(i+(lat-1)*LONR+(k-1)*LATR)
!rab            enddo
!rab          enddo
!rab        enddo
!rab      endif

!--- call the nuopc radiation loop---
    call mpp_clock_begin(grrad_clk)
    do nb = 1, Atm_block%nblks


!rab          if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
!rab            do i = 1, njeff
!rab              Dyn_parms(nb)%icsdsw(i) = ixseed(lon+i-1,lan,1)
!rab              Dyn_parms(nb)%icsdlw(i) = ixseed(lon+i-1,lan,2)
!rab            enddo
!rab          endif


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

!rab   if (first) then
!--- set up the random number streams needed or each column
!rab        if (ras) call ras_init(levs, me)
!rab        if (.not. newsas .or. cal_pre) then  ! random number needed for RAS and old SAS
!rab          if (random_clds) then ! create random number tank
!rab!                                 -------------------------
!rab            seed0 = idate(1) + idate(2) + idate(3) + idate(4)
!rab
!rab            call random_setseed(seed0)
!rab            call random_number(wrk)
!rab            seed0 = seed0 + nint(wrk(1)*thousnd)
!rab            if (me == 0) print *,' seed0=',seed0,' idate=',idate,
!rab     &                           ' wrk=',wrk
!rab!
!rab            if (.not. allocated(rannum_tank))
!rab     &                allocate (rannum_tank(lonr,maxran,lats_node_r))
!rab            if (.not. allocated(rannum)) allocate (rannum(lonr*maxrs))
!rab            lonrbm = lonr / maxsub
!rab            if (me == 0) write(0,*)' maxran=',maxran,' maxrs=',maxrs,
!rab     &          'maxsub=',maxsub,' lonrbm=',lonrbm,
!rab     &          ' lats_node_r=',lats_node_r
!rab            do j=1,lats_node_r
!rab              iseedl = global_lats_r(ipt_lats_node_r-1+j) + seed0
!rab              call random_setseed(iseedl)
!rab              call random_number(rannum)
!rab
!rab!$omp parallel do  shared(j,lonr,lonrbm,rannum,rannum_tank)
!rab!$omp+private(nrc,nnp,i,ii,k,kk)
!rab              do nrc=1,maxrs
!rab             nnp = (nrc-1)*lonr
!rab                do k=1,maxsub
!rab                  kk = k - 1
!rab                  do i=1,lonr
!rab                    ii = kk*lonrbm + i
!rab                    if (ii > lonr) ii = ii - lonr
!rab                    rannum_tank(i,nrc+kk*maxrs,j) = rannum(ii+nnp)
!rab                  enddo
!rab                enddo
!rab              enddo
!rab            enddo
!rab            if (allocated(rannum)) deallocate (rannum)
!rab          endif
!rab        endif
!rab   endif


!--- call the nuopc physics loop---
    call mpp_clock_begin(gbphys_clk)
    do nb = 1, Atm_block%nblks

!rab      nnp = 1
!rab      if (.not. newsas .or. cal_pre) then
!rab        if (random_clds) then
!rab          nnr = (nnp-1)*nrcm
!rab          do j=1,nrcm
!rab            do i=1,njeff
!rab              Tbd_data(nb)%(rann(i,j)=rannum_tank(lon+i-1,indxr(nnr+j),lan)
!rab            enddo
!rab          enddo
!rab        else
!rab          do j=1,nrcm
!rab            do i=1,njeff
!rab              Tbd_data(nb)%rann(i,j) = 0.6    ! This is useful for debugging
!rab            enddo
!rab          enddo
!rab        endif
!rab      endif

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
  subroutine surface_props_input (Atm_block)
    type (block_control_type), intent(in) :: Atm_block
!   local variables
    integer :: i, j, ibs, ibe, jbs, jbe, nct
    integer :: funit, nb, nx, ny, ngptc
    integer :: start(4), nread(4)
    character(len=32)  :: fname = 'INPUT/gfs_sfc.nc'
    character(len=128) :: errmsg
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3 => NULL()
    logical :: exists

    call get_mosaic_tile_file (fname, fname, .FALSE.)

    inquire(file=trim(fname), exist=exists)
    if (exists) then
      errmsg = 'opening file '//trim(fname)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), NOTE)
    else
      errmsg = 'error opening file '//trim(fname)//' for input'
      call error_mesg ('phys_rad_driver_init', trim(errmsg), FATAL)
    endif

    funit = open_file(trim(fname))

    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)-Atm_block%isc+1
      ibe = Atm_block%ibe(nb)-Atm_block%isc+1
      jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
      jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
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
!--- slmsk
      var2(1:nx,1:ny) => Sfc_props(nb)%slmsk(1:ngptc)
      call read_data(fname,'slmsk',var2,start,nread)
!--- oro (orog in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%oro(1:ngptc)
      call read_data(fname,'orog',var2,start,nread)
!--- tsfc (tsea in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%tsfc(1:ngptc)
      call read_data(fname,'tsea',var2,start,nread)
!--- weasd (sheleg in sfc file)
      var2(1:nx,1:ny) => Sfc_props(nb)%weasd(1:ngptc)
      call read_data(fname,'sheleg',var2,start,nread)
!--- tg3
      var2(1:nx,1:ny) => Sfc_props(nb)%tg3(1:ngptc)
      call read_data(fname,'tg3',var2,start,nread)
!--- zorl
      var2(1:nx,1:ny) => Sfc_props(nb)%zorl(1:ngptc)
      call read_data(fname,'zorl',var2,start,nread)
!--- alvsf
      var2(1:nx,1:ny) => Sfc_props(nb)%alvsf(1:ngptc)
      call read_data(fname,'alvsf',var2,start,nread)
!--- alvwf
      var2(1:nx,1:ny) => Sfc_props(nb)%alvwf(1:ngptc)
      call read_data(fname,'alvwf',var2,start,nread)
!--- alnsf
      var2(1:nx,1:ny) => Sfc_props(nb)%alnsf(1:ngptc)
      call read_data(fname,'alnsf',var2,start,nread)
!--- alnwf
      var2(1:nx,1:ny) => Sfc_props(nb)%alnwf(1:ngptc)
      call read_data(fname,'alnwf',var2,start,nread)
!--- vfrac
      var2(1:nx,1:ny) => Sfc_props(nb)%vfrac(1:ngptc)
      call read_data(fname,'vfrac',var2,start,nread)
!--- canopy
      var2(1:nx,1:ny) => Sfc_props(nb)%canopy(1:ngptc)
      call read_data(fname,'canopy',var2,start,nread)
!--- f10m
      var2(1:nx,1:ny) => Sfc_props(nb)%f10m(1:ngptc)
      call read_data(fname,'f10m',var2,start,nread)
!--- t2m
      var2(1:nx,1:ny) => Sfc_props(nb)%t2m(1:ngptc)
      call read_data(fname,'t2m',var2,start,nread)
!--- q2m
      var2(1:nx,1:ny) => Sfc_props(nb)%q2m(1:ngptc)
      call read_data(fname,'q2m',var2,start,nread)
!--- vtype
      var2(1:nx,1:ny) => Sfc_props(nb)%vtype(1:ngptc)
      call read_data(fname,'vtype',var2,start,nread)
!--- stype
      var2(1:nx,1:ny) => Sfc_props(nb)%stype(1:ngptc)
      call read_data(fname,'stype',var2,start,nread)
!--- facsf
      var2(1:nx,1:ny) => Sfc_props(nb)%facsf(1:ngptc)
      call read_data(fname,'facsf',var2,start,nread)
!--- facwf
      var2(1:nx,1:ny) => Sfc_props(nb)%facwf(1:ngptc)
      call read_data(fname,'facwf',var2,start,nread)
!--- uustar
      var2(1:nx,1:ny) => Sfc_props(nb)%uustar(1:ngptc)
      call read_data(fname,'uustar',var2,start,nread)
!--- ffmm
      var2(1:nx,1:ny) => Sfc_props(nb)%ffmm(1:ngptc)
      call read_data(fname,'ffmm',var2,start,nread)
!--- ffhh
      var2(1:nx,1:ny) => Sfc_props(nb)%ffhh(1:ngptc)
      call read_data(fname,'ffhh',var2,start,nread)
!--- hice
      var2(1:nx,1:ny) => Sfc_props(nb)%hice(1:ngptc)
      call read_data(fname,'hice',var2,start,nread)
!--- fice
      var2(1:nx,1:ny) => Sfc_props(nb)%fice(1:ngptc)
      call read_data(fname,'fice',var2,start,nread)
!--- tisfc
      var2(1:nx,1:ny) => Sfc_props(nb)%tisfc(1:ngptc)
      call read_data(fname,'tisfc',var2,start,nread)
!--- tprcp
      var2(1:nx,1:ny) => Tbd_data(nb)%tprcp(1:ngptc)
      call read_data(fname,'tprcp',var2,start,nread)
!--- srflag
      var2(1:nx,1:ny) => Tbd_data(nb)%srflag(1:ngptc)
      call read_data(fname,'srflag',var2,start,nread)
!--- snowd (snwdph in the file)
      var2(1:nx,1:ny) => Sfc_props(nb)%snowd(1:ngptc)
      call read_data(fname,'snwdph',var2,start,nread)
!--- shdmin
      var2(1:nx,1:ny) => Sfc_props(nb)%shdmin(1:ngptc)
      call read_data(fname,'shdmin',var2,start,nread)
!--- shdmax
      var2(1:nx,1:ny) => Sfc_props(nb)%shdmax(1:ngptc)
      call read_data(fname,'shdmax',var2,start,nread)
!--- slope
      var2(1:nx,1:ny) => Sfc_props(nb)%slope(1:ngptc)
      call read_data(fname,'slope',var2,start,nread)
!--- snoalb
      var2(1:nx,1:ny) => Sfc_props(nb)%snoalb(1:ngptc)
      call read_data(fname,'snoalb',var2,start,nread)
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
      call read_data(fname,'stc',var3,start,nread)
      call read_data(fname,'stc',var3)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%stc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
!--- smc
      call read_data(fname,'smc',var3,start,nread)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%smc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
!--- slc
      call read_data(fname,'slc',var3,start,nread)
      do j = 1, ny
       do i = 1, nx
         nct = (j-1)*nx + i
         Tbd_data(nb)%slc(nct,1:Mdl_parms%lsoil) = var3(i,j,1:Mdl_parms%lsoil)
       enddo
      enddo
    enddo


!--- nullify/deallocate any temporaries used
      nullify(var2)
      deallocate(var3)

!--- close the file
      call close_file(funit)

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
