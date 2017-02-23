module IPD_driver

  use IPD_typedefs, only: IPD_control_type,  IPD_data_type,    &
                          IPD_diag_type,     IPD_restart_type

  use GFS_driver,   only: IPD_init_type => GFS_init_type

  use GFS_driver,   only: GFS_initialize,     GFS_setup_step,  &
                          GFS_radiation_step, GFS_physics_step

       implicit none

!------------------------------------------------------!
!  IPD containers                                      !
!------------------------------------------------------!
!  type(GFS_control_type)              :: IPD_Control  !
!  type(IPD_data_type)     allocatable :: IPD_Data(:)  !
!  type(IPD_diag_type),                :: IPD_Diag(:)  !
!  type(IPD_restart_type),             :: IPD_Restart  !
!------------------------------------------------------!

!----------------
! Public Entities
!----------------
! functions
  public IPD_initialize, IPD_setup_step, IPD_radiation_step, IPD_physics_step

  CONTAINS
!*******************************************************************************************


!----------------
!  IPD Initialize 
!----------------
  subroutine IPD_initialize (IPD_control, IPD_Data, IPD_Diag, IPD_Restart, GFS_init_parm)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart
    type(IPD_init_type),    intent(in)    :: GFS_init_parm

    call GFS_initialize (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                         IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                         IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                         IPD_Data(:)%Intdiag, IPD_Diag, IPD_Restart, GFS_init_parm)

  end subroutine IPD_initialize


!---------------------------------------------
!  IPD setup step
!    surface data cycling, random streams, etc
!---------------------------------------------
  subroutine IPD_setup_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call GFS_setup_step (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                         IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                         IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                         IPD_Data(:)%Intdiag, IPD_Diag, IPD_Restart)

  end subroutine IPD_setup_step


!--------------------
!  IPD radiation step
!--------------------
  subroutine IPD_radiation_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call GFS_radiation_step (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                             IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                             IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                             IPD_Data%Intdiag, IPD_Diag, IPD_Restart)

  end subroutine IPD_radiation_step


!------------------
!  IPD physics step
!------------------
  subroutine IPD_physics_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart

    call GFS_physics_step (IPD_control, IPD_Data%Statein, IPD_Data%Stateout,   &
                           IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                           IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                           IPD_Data%Intdiag, IPD_Diag, IPD_Restart)

  end subroutine IPD_physics_step

end module IPD_driver
