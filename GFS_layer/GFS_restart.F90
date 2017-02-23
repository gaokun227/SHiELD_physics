module GFS_restarts

  use machine,       only: kind_phys
  use IPD_typedefs,  only: IPD_restart_type
  use GFS_typedefs,  only: GFS_tbd_type, Gfs_coupling_type, GFS_control_type

  public GFS_populate_IPD_restart

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine GFS_populate_IPD_restart (IPD_Restart, Tbd, Coupling, Model, blksz)
!----------------------------------------------------------------------------------------!
!   IPD_METADATA                                                                         !
!     IPD_Restart%num2d          [int*4  ]  number of 2D variables to output             !
!     IPD_Restart%num3d          [int*4  ]  number of 3D variables to output             !
!     IPD_Restart%name2d         [char=32]  variable name in restart file                !
!     IPD_Restart%name3d         [char=32]  variable name in restart file                !
!     IPD_Restart%fld2d(:,:,:)   [real*8 ]  pointer to 2D data (im,nblks,MAX_RSTRT)      !
!     IPD_Restart%fld3d(:,:,:,:) [real*8 ]  pointer to 3D data (im,levs,nblks,MAX_RSTRT) !
!----------------------------------------------------------------------------------------!
    type(IPD_restart_type),  intent(inout) :: IPD_Restart
    type(GFS_tbd_type),      intent(in)    :: Tbd(:)
    type(GFS_coupling_type), intent(in)    :: Coupling(:)
    type(GFS_control_type),  intent(in)    :: Model
    integer,                 intent(in)    :: blksz(:)
    !--- local variables
    integer :: nblks, num, nb, max_rstrt
    character(len=2) :: c2 = ''
    
    nblks = size(blksz)
    max_rstrt = size(IPD_Restart%name2d)

    IPD_Restart%num2d = Model%ntot2d + Model%nctp
    IPD_Restart%num3d = Model%ntot3d

    allocate (IPD_Restart%name2d(IPD_Restart%num2d))
    allocate (IPD_Restart%name3d(IPD_Restart%num3d))
    allocate (IPD_Restart%data(nblks,max(IPD_Restart%num2d,IPD_Restart%num3d)))

    IPD_Restart%name2d(:) = ' '
    IPD_Restart%name3d(:) = ' '

    !--- phy_f2d variables
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name2d(num) = 'phy_f2d_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num)%var2p => Tbd(nb)%phy_f2d(:,num)
      enddo
    enddo

    !--- phy_fctd variables
    do num = 1, Model%nctp
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name2d(num+Model%ntot2d) = 'phy_fctd_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num+Model%ntot2d)%var2p => Tbd(nb)%phy_fctd(:,num)
      enddo
    enddo

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      IPD_Restart%name3d(num) = 'phy_f3d_'//c2
      do nb = 1,nblks
        IPD_Restart%data(nb,num)%var3p => Tbd(nb)%phy_f3d(:,:,num)
      enddo
    enddo

  end subroutine GFS_populate_IPD_restart

end module GFS_restarts
