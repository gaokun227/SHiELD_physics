module IPD_typedefs
  use machine,  only: kind_phys

  use GFS_typedefs, only: IPD_control_type => GFS_control_type, &
                          IPD_init_type    => GFS_init_type,    &
                          GFS_statein_type,  GFS_stateout_type, &
                          GFS_sfcprop_type,  GFS_coupling_type, &
                          GFS_grid_type,     GFS_tbd_type,      &
                          GFS_cldprop_type,  GFS_radtend_type,  &
                          GFS_diag_type

!--------------------
!  IPD sub-containers
!--------------------
  type IPD_data_type
    type(GFS_statein_type)  :: Statein
    type(GFS_stateout_type) :: Stateout
    type(GFS_sfcprop_type)  :: Sfcprop
    type(GFS_coupling_type) :: Coupling
    type(GFS_grid_type)     :: Grid
    type(GFS_tbd_type)      :: Tbd
    type(GFS_cldprop_type)  :: Cldprop
    type(GFS_radtend_type)  :: Radtend
    type(GFS_diag_type)     :: Intdiag
  end type IPD_data_type


  type var_subtype
    real(kind=kind_phys), pointer :: var2p(:)   => null()  !< 2D data saved in packed format [dim(ix)]
    real(kind=kind_phys), pointer :: var3p(:,:) => null()  !< 3D data saved in packed format [dim(ix,levs)]
  end type var_subtype
    
!-------------------------------------------
! IPD_restart_type
!   data necessary for reproducible restarts
!-------------------------------------------
  type IPD_restart_type
    integer           :: num2d                    !< current number of registered 2D restart variables
    integer           :: num3d                    !< current number of registered 3D restart variables
    character(len=32), allocatable :: name2d(:)   !< variable name as it will appear in the restart file
    character(len=32), allocatable :: name3d(:)   !< variable name as it will appear in the restart file
    type(var_subtype), allocatable :: data(:,:)   !< holds pointers to data in packed format (allocated to (nblks,max(2d/3dfields))
  end type IPD_restart_type

!----------------------------------------
! IPD_diag_type
!   fields targetted as diagnostic output
!----------------------------------------
  type IPD_diag_type
    character(len=32)     :: name           !< variable name in source
    character(len=32)     :: output_name    !< output name for variable
    character(len=32)     :: mod_name       !< module name (e.g. physics, radiation, etc)
    character(len=32)     :: file_name      !< output file name for variable
    character(len=128)    :: desc           !< long description of field
    character(len=32)     :: unit           !< units associated with fields
    character(len=32)     :: type_stat_proc !< type of statistic processing:
                                            !< average, accumulation, maximal, minimal, etc.
    character(len=32)     :: level_type     !< vertical level of the field
    integer               :: level          !< vertical level(s)
    real(kind=kind_phys)  :: cnvfac         !< conversion factors to output in specified units
    real(kind=kind_phys)  :: zhour          !< forecast hour when bucket was last emptied for statistical processing
    real(kind=kind_phys)  :: fcst_hour      !< current forecast hour (same as fhour)
    type(var_subtype), allocatable :: data(:) !< holds pointers to data in packed format (allocated to nblks)
  end type IPD_diag_type

  public kind_phys
  public IPD_control_type
  public IPD_data_type
  public IPD_restart_type
  public IPD_diag_type
  public IPD_init_type

  CONTAINS
!*******************************************************************************************

end module IPD_typedefs
