#ifdef OVERLOAD_R4
#define _GET_VAR1 get_var1_real 
#else
#define _GET_VAR1 get_var1_double
#endif

module external_ic_mod

   use external_sst_mod,   only: i_sst, j_sst, sst_ncep
   use fms_mod,            only: file_exist, read_data, field_exist, write_version_number
   use fms_mod,            only: open_namelist_file, check_nml_error, close_file
   use fms_mod,            only: get_mosaic_tile_file, read_data, error_mesg
   use fms_io_mod,         only: get_tile_string, field_size, free_restart_type
   use fms_io_mod,         only: restart_file_type, register_restart_field
   use fms_io_mod,         only: save_restart, restore_state, set_filename_appendix
   use mpp_mod,            only: mpp_error, FATAL, NOTE, mpp_pe
   use mpp_mod,            only: stdlog, input_nml_file
   use mpp_parameter_mod,  only: AGRID_PARAM=>AGRID
   use mpp_domains_mod,    only: mpp_get_tile_id, domain2d, mpp_update_domains
   use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
   use tracer_manager_mod, only: set_tracer_profile
   use field_manager_mod,  only: MODEL_ATMOS

   use constants_mod,     only: pi, omega, grav, kappa, rdgas, rvgas, cp_air
   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, R_GRID
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_grid_utils_mod, only: ptop_min, g_sum
   use fv_io_mod,         only: fv_io_read_tracers 
   use fv_mapz_mod,       only: mappm
   use fv_mp_mod,         only: ng, is_master, fill_corners, YDir
   use fv_surf_map_mod,   only: surfdrv
   use fv_timing_mod,     only: timing_on, timing_off
   use init_hydro_mod,    only: p_var
   use sim_nc_mod,        only: open_ncfile, close_ncfile, get_ncdim1, get_var1_double, get_var2_real,   &
                                get_var3_r4, get_var1_real
   use fv_nwp_nudge_mod,  only: T_is_Tv
! The "T" field in NCEP analysis is actually virtual temperature (Larry H. post processing)
! BEFORE 20051201

   use boundary_mod,      only: nested_grid_BC
   use mpp_domains_mod,       only: mpp_get_data_domain, mpp_get_global_domain, mpp_get_compute_domain

   implicit none
   private

   real, parameter:: zvir = rvgas/rdgas - 1.
   real :: deg2rad

   public get_external_ic, get_cubed_sphere_terrain

!---- version number -----
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

contains

   subroutine get_external_ic( Atm, fv_domain, cold_start )

      type(fv_atmos_type), intent(inout), target :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      logical, intent(IN) :: cold_start
      real:: alpha = 0.
      real rdg
      integer i,j,k,nq

      real, pointer, dimension(:,:,:) :: grid, agrid
      real, pointer, dimension(:,:) :: fC, f0

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie
      js  = Atm(1)%bd%js
      je  = Atm(1)%bd%je
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed

      grid  => Atm(1)%gridstruct%grid
      agrid => Atm(1)%gridstruct%agrid

      fC    => Atm(1)%gridstruct%fC
      f0    => Atm(1)%gridstruct%f0

! * Initialize coriolis param:
 
      do j=jsd,jed+1
         do i=isd,ied+1
            fc(i,j) = 2.*omega*( -1.*cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) + &
                                     sin(grid(i,j,2))*cos(alpha) )
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            f0(i,j) = 2.*omega*( -1.*cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) + &
                                     sin(agrid(i,j,2))*cos(alpha) )
         enddo
      enddo

      call mpp_update_domains( f0, fv_domain )
      if ( Atm(1)%gridstruct%cubed_sphere .and. .not. Atm(1)%neststruct%nested) call fill_corners(f0, Atm(1)%npx, Atm(1)%npy, YDir)
 
! Read in cubed_sphere terrain
      if ( Atm(1)%flagstruct%mountain ) then
           call get_cubed_sphere_terrain(Atm, fv_domain)
      else
           Atm(1)%phis = 0.
      endif
 
! Read in the specified external dataset and do all the needed transformation
      if ( Atm(1)%flagstruct%ncep_ic ) then
           nq = 1
                             call timing_on('NCEP_IC')
           call get_ncep_ic( Atm, fv_domain, nq )
                             call timing_off('NCEP_IC')
#ifndef NO_FV_TRACERS
           if (.not. cold_start) then
              call fv_io_read_tracers( fv_domain, Atm )
              if(is_master()) write(*,*) 'All tracers except sphum replaced by FV IC'
           endif
#endif
      elseif ( Atm(1)%flagstruct%nggps_ic ) then
                             call timing_on('NGGPS_IC')
           call get_nggps_ic( Atm, fv_domain )
                             call timing_off('NGGPS_IC')
      elseif ( Atm(1)%flagstruct%fv_diag_ic ) then
! Interpolate/remap diagnostic output from a FV model diagnostic output file on uniform lat-lon A grid:
               nq = size(Atm(1)%q,4)
! Needed variables: lon, lat, pfull(dim), zsurf, ps, ucomp, vcomp, w, temp, and all q
! delz not implemnetd yet; set make_nh = .true.
               call get_diag_ic( Atm, fv_domain, nq )
      else
! The following is to read in legacy lat-lon FV core restart file
!  is Atm%q defined in all cases?

           nq = size(Atm(1)%q,4)
           call get_fv_ic( Atm, fv_domain, nq )
      endif

      call prt_maxmin('PS', Atm(1)%ps, is, ie, js, je, ng, 1, 0.01)
      call prt_maxmin('T', Atm(1)%pt, is, ie, js, je, ng, Atm(1)%npz, 1.)
      call prt_maxmin('W', Atm(1)%w, is, ie, js, je, ng, Atm(1)%npz, 1.)
      call prt_maxmin('SPHUM', Atm(1)%q(:,:,:,1), is, ie, js, je, ng, Atm(1)%npz, 1.)
      call prt_maxmin('O3MR', Atm(1)%q(:,:,:,2), is, ie, js, je, ng, Atm(1)%npz, 1.)
      call prt_maxmin('CLWMR', Atm(1)%q(:,:,:,3), is, ie, js, je, ng, Atm(1)%npz, 1.)

      call p_var(Atm(1)%npz,  is, ie, js, je, Atm(1)%ak(1),  ptop_min,         &
                 Atm(1)%delp, Atm(1)%delz, Atm(1)%pt, Atm(1)%ps,               &
                 Atm(1)%pe,   Atm(1)%peln, Atm(1)%pk, Atm(1)%pkz,              &
                 kappa, Atm(1)%q, ng, Atm(1)%ncnst, Atm(1)%gridstruct%area_64, Atm(1)%flagstruct%dry_mass,           &
                 Atm(1)%flagstruct%adjust_dry_mass, Atm(1)%flagstruct%mountain, Atm(1)%flagstruct%moist_phys,   &
                 Atm(1)%flagstruct%hydrostatic, Atm(1)%flagstruct%nwat, Atm(1)%domain, Atm(1)%flagstruct%make_nh)

  end subroutine get_external_ic



  subroutine get_cubed_sphere_terrain( Atm, fv_domain )
    type(fv_atmos_type), intent(inout), target :: Atm(:)
    type(domain2d),      intent(inout) :: fv_domain
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)
    character(len=64)    :: fname
    character(len=3)  :: gn
    integer              ::  n
    integer              ::  jbeg, jend
    real ftop
    real, allocatable :: g_dat2(:,:,:)
    real, allocatable :: pt_coarse(:,:,:)
    integer isc_p, iec_p, jsc_p, jec_p, isg, ieg, jsg,jeg

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie
      js  = Atm(1)%bd%js
      je  = Atm(1)%bd%je
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed

    if (Atm(1)%grid_number > 1) then
       write(gn,'(A2, I1)') ".g", Atm(1)%grid_number
    else
       gn = ''
    end if

    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
                            ! always one at this point

    allocate( tile_id(ntileMe) )
    tile_id = mpp_get_tile_id( fv_domain )
    do n=1,ntileMe

       call get_tile_string(fname, 'INPUT/fv_core'//trim(gn)//'.res.tile', tile_id(n), '.nc' )

       
       if( file_exist(fname) ) then
          call read_data(fname, 'phis', Atm(n)%phis(is:ie,js:je),      &
                         domain=fv_domain, tile_count=n)
       else
          call surfdrv(  Atm(n)%npx, Atm(n)%npy, Atm(n)%gridstruct%grid_64, Atm(n)%gridstruct%agrid_64,   &
                         Atm(n)%gridstruct%area_64, Atm(n)%gridstruct%dx, Atm(n)%gridstruct%dy, &
                         Atm(n)%gridstruct%dxc, Atm(n)%gridstruct%dyc, Atm(n)%gridstruct%sin_sg, &
                         Atm(n)%phis, Atm(n)%flagstruct%stretch_fac, &
                         Atm(n)%neststruct%nested, Atm(n)%neststruct%npx_global, Atm(N)%domain, &
                         Atm(n)%flagstruct%grid_number, Atm(n)%bd )
          call mpp_error(NOTE,'terrain datasets generated using USGS data')
       endif

    end do
 
    call mpp_update_domains( Atm(1)%phis, Atm(1)%domain )
    if (Atm(1)%neststruct%nested) then
       call mpp_get_compute_domain( Atm(1)%parent_grid%domain, &
            isc_p,  iec_p,  jsc_p,  jec_p  )
!!$       call mpp_get_data_domain( Atm(1)%parent_grid%domain, &
!!$            isc_p,  iec_p,  jsc_p,  jec_p  )
       call mpp_get_global_domain( Atm(1)%parent_grid%domain, &
            isg, ieg, jsg, jeg)

       allocate(g_dat2( isg:ieg, jsg:jeg,1) )
       
       g_dat2 = 0.
       if (ANY(mpp_pe()==Atm(1)%parent_grid%pelist)) g_dat2(isc_p:iec_p,  jsc_p:jec_p  , 1) = Atm(1)%parent_grid%phis
       call nested_grid_BC(Atm(1)%phis, g_dat2(:,:,1), Atm(1)%neststruct%nest_domain, &
         Atm(1)%neststruct%ind_h, Atm(1)%neststruct%wt_h, 0, 0, &
         Atm(1)%npx, Atm(1)%npy, Atm(1)%bd, isg, ieg, jsg, jeg, proc_in=.true.)

       deallocate(g_dat2)
    end if
    ftop = g_sum(Atm(1)%domain, Atm(1)%phis(is:ie,js:je), is, ie, js, je, ng, Atm(1)%gridstruct%area_64, 1)
 
    call prt_maxmin('ZS', Atm(1)%phis,  is, ie, js, je, ng, 1, 1./grav)
    if(is_master()) write(*,*) 'mean terrain height (m)=', ftop/grav
 
    deallocate( tile_id )

  end subroutine get_cubed_sphere_terrain


  subroutine get_diag_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq
! local:
      character(len=128) :: fname, tracer_name
      real(kind=4), allocatable:: wk1(:), wk2(:,:), wk3(:,:,:)
      real, allocatable:: tp(:,:,:), qp(:,:,:,:)
      real, allocatable:: ua(:,:,:), va(:,:,:), wa(:,:,:)
      real, allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      real:: s2c(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je,4)
      integer, dimension(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je):: id1, id2, jdc
      real psc(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je)
      real gzc(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je)
      integer:: i, j, k, im, jm, km, npz, npt
      integer:: i1, i2, j1, ncid
      integer:: jbeg, jend
      integer tsize(3), tr_ind
      logical:: found

      integer  sphum, liq_wat, rainwat, ice_wat, snowwat, graupel

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie
      js  = Atm(1)%bd%js
      je  = Atm(1)%bd%je
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed

      deg2rad = pi/180.

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

      fname = Atm(1)%flagstruct%res_latlon_dynamics

      if( file_exist(fname) ) then
          call open_ncfile( fname, ncid )        ! open the file
          call get_ncdim1( ncid, 'lon',   tsize(1) )
          call get_ncdim1( ncid, 'lat',   tsize(2) )
          call get_ncdim1( ncid, 'pfull', tsize(3) )

          im = tsize(1); jm = tsize(2); km = tsize(3)

          if(is_master())  write(*,*) fname, ' FV_diag IC dimensions:', tsize

          allocate (  lon(im) )
          allocate (  lat(jm) )
 
          call _GET_VAR1 (ncid, 'lon', im, lon )
          call _GET_VAR1 (ncid, 'lat', jm, lat )

! Convert to radian
          do i=1,im
             lon(i) = lon(i) * deg2rad  ! lon(1) = 0.
          enddo
          do j=1,jm
             lat(j) = lat(j) * deg2rad
          enddo

          allocate ( ak0(km+1) )
          allocate ( bk0(km+1) )
! npz
          if ( npz /= km ) then
               call mpp_error(FATAL,'==>Error in get_diag_ic: vertical dim must be the same')
          else
               ak0(:) = Atm(1)%ak(:)
               bk0(:) = Atm(1)%bk(:)
          endif
      else
          call mpp_error(FATAL,'==> Error in get_diag_ic: Expected file '//trim(fname)//' for dynamics does not exist')
      endif

! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      call remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c , Atm(1)%gridstruct%agrid, Atm(1)%bd)

! Find bounding latitudes:
      jbeg = jm-1;         jend = 2
      do j=js,je
         do i=is,ie
              j1 = jdc(i,j)
            jbeg = min(jbeg, j1) 
            jend = max(jend, j1+1)
         enddo
      enddo

! remap surface pressure and height:
      allocate ( wk2(im,jbeg:jend) )
      call get_var3_r4( ncid, 'ps', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            psc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                       s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      call get_var3_r4( ncid, 'zsurf', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            gzc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                       s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo
      deallocate ( wk2 )

! Read in temperature:
      allocate ( wk3(1:im,jbeg:jend, 1:km) )
      call get_var3_r4( ncid, 'temp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate (  tp(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
        enddo
      enddo

! Read in all tracers:
    allocate ( qp(is:ie,js:je,km,nq) )
    qp = 1.E25
    do tr_ind=1, nq
       call get_tracer_names(MODEL_ATMOS, tr_ind, tracer_name)
       if (field_exist(fname,tracer_name)) then
           call get_var3_r4( ncid, tracer_name, 1,im, jbeg,jend, 1,km, wk3 )
           do k=1,km
              do j=js,je
                 do i=is,ie
                    i1 = id1(i,j)
                    i2 = id2(i,j)
                    j1 = jdc(i,j)
                    qp(i,j,k,tr_ind) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                 enddo
              enddo
           enddo
       endif
    enddo
    call remap_scalar(im, jm, km, npz, nq, nq, ak0, bk0, psc, gzc, tp, qp, Atm(1))
    deallocate ( tp )
    deallocate ( qp )

! Winds:
      call get_var3_r4( ncid, 'ucomp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate ( ua(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ua(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
          enddo
        enddo
      enddo

      call get_var3_r4( ncid, 'vcomp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate ( va(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            va(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
          enddo
        enddo
      enddo
      call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm(1))

      deallocate ( ua )
      deallocate ( va )

      if ( .not. Atm(1)%flagstruct%hydrostatic ) then
        if (field_exist(fname,'w')) then
           allocate ( wa(is:ie,js:je,km) )
           call get_var3_r4( ncid, 'w', 1,im, jbeg,jend, 1,km, wk3 )
           do k=1,km
              do j=js,je
                 do i=is,ie
                    i1 = id1(i,j)
                    i2 = id2(i,j)
                    j1 = jdc(i,j)
                    wa(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                                s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                 enddo
              enddo
           enddo
           call remap_wz(im, jm, km, npz, ng, ak0, bk0, psc, wa, Atm(1)%w, Atm(1))
           deallocate ( wa )
        else    ! use "w = - fac * omega" ?
           Atm(1)%w(:,:,:) = 0.
        endif
! delz:
        if (field_exist(fname,'delz')) then
           allocate ( wa(is:ie,js:je,km) )
           call get_var3_r4( ncid, 'delz', 1,im, jbeg,jend, 1,km, wk3 )
           do k=1,km
              do j=js,je
                 do i=is,ie
                    i1 = id1(i,j)
                    i2 = id2(i,j)
                    j1 = jdc(i,j)
                    wa(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                                s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                 enddo
              enddo
           enddo
           call remap_wz(im, jm, km, npz, 0,  ak0, bk0, psc, wa, Atm(1)%delz, Atm(1))
           deallocate ( wa )
        else    ! Force make = T
           Atm(1)%flagstruct%make_nh = .true.
        endif

      endif   ! hydrostatic test

      call close_ncfile ( ncid )
      deallocate ( wk3 )
      deallocate ( ak0 )
      deallocate ( bk0 )
      deallocate ( lat )
      deallocate ( lon )

  end subroutine get_diag_ic



  subroutine get_nggps_ic (Atm, fv_domain)
!    read in data after it has been preprocessed with 
!    NCEP/EMC orography maker and global_chgres
!    and has been horiztontally interpolated to the 
!    current cubed-sphere grid
!
!--- variables read in from 'gfs_ctrl.nc'
!       VCOORD  -  level information
!                   maps to 'ak & bk'
!--- variables read in from 'sfc_data.nc'
!       land_frac  -  land-sea-ice mask (L:0 / S:1)
!                     maps to 'oro'
!       TSEA       -  surface skin temperature (k)
!                     maps to 'ts'
!--- variables read in from 'gfs_data.nc'
!       ZS  -  surface height (m)
!              maps to 'ze0' for hybrid_z non-hydrostatic(?)
!              maps to phis = zs*grav for case where we use
!       PS  -  surface pressure (Pa)
!       DP  -  prognostic delta-p (Pa)
!              maps to 'delp'
!       P   -  prognostic pressure (Pa)
!              do I map this at all?
!       T   -  prognostic temperature (k)
!       U   -  prognostic zonal wind (m/s)
!              maps to 'ua'
!       V   -  prognostic meridional wind (m/s)
!              maps to 'va'
!       W   -  prognostic 'omega' (Pa/s)
!       Q   -  prognostic tracer fields (Specific Humidity, 
!                                        O3 mixing ratio,
!                                        Cloud mixing ratio)
!--- Namelist variables 
!       filtered_terrain  -  use orography maker filtered terrain mapping
!       ncep_terrain      -  use NCEP ICs interpolated terrain 
!       ncep_plevels      -  use NCEP pressure levels (implies no vertical remapping)


      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
! local:
      real, dimension(64):: ak_sj, bk_sj
      real, dimension(:), allocatable:: ak, bk
      real, dimension(:,:), allocatable:: wk2, zs, ps
      real, dimension(:,:,:), allocatable:: dp, t, ua, va, omga
      real, dimension(:,:,:,:), allocatable:: q
      real rdg
      integer:: n, npz, itoa, nt, ntprog, ntdiag, ntracers, ntrac
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed
      integer :: ios, ierr, unit, id_res
      type (restart_file_type) :: SFC_restart, GFS_restart
      character(len=6)  :: gn, stile_name
      character(len=64) :: tracer_name
      character(len=64) :: fn_gfs_ctl = 'gfs_ctrl.nc'
      character(len=64) :: fn_gfs_ics = 'gfs_data.nc'
      character(len=64) :: fn_sfc_ctl = 'sfc_ctrl.nc'
      character(len=64) :: fn_sfc_ics = 'sfc_data.nc'
      logical :: remap
      logical :: filtered_terrain = .true.
      logical :: ncep_terrain = .false.
      logical :: ncep_plevels = .false.
      integer :: levp = 64
      namelist /external_ic_nml/ filtered_terrain, ncep_terrain, ncep_plevels, levp
! This is activated by USE_GFSL63
! Thfollowing L63 setting is the same as NCEP GFS's L64 except the top
! 3 layers
      data ak_sj/25.00000,     100.00000,     200.00000,    &
                311.00000,     430.00000,     558.00000,    &
                700.00000,     863.05803,    1051.07995,    &
               1265.75194,    1510.71101,    1790.05098,    &
               2108.36604,    2470.78817,    2883.03811,    &
               3351.46002,    3883.05187,    4485.49315,    &
               5167.14603,    5937.04991,    6804.87379,    &
               7780.84698,    8875.64338,   10100.20534,    &
              11264.35673,   12190.64366,   12905.42546,    &
              13430.87867,   13785.88765,   13986.77987,    &
              14047.96335,   13982.46770,   13802.40331,    &
              13519.33841,   13144.59486,   12689.45608,    &
              12165.28766,   11583.57006,   10955.84778,    &
              10293.60402,    9608.08306,    8910.07678,    &
               8209.70131,    7516.18560,    6837.69250,    &
               6181.19473,    5552.39653,    4955.72632,    &
               4394.37629,    3870.38682,    3384.76586,    &
               2937.63489,    2528.37666,    2155.78385,    &
               1818.20722,    1513.68173,    1240.03585,    &
                994.99144,     776.23591,     581.48797,    &
                408.53400,     255.26520,     119.70243, 0. /

      data bk_sj/0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00000,       0.00000,       0.00000,    &
                 0.00201,       0.00792,       0.01755,    &
                 0.03079,       0.04751,       0.06761,    &
                 0.09097,       0.11746,       0.14690,    &
                 0.17911,       0.21382,       0.25076,    &
                 0.28960,       0.32994,       0.37140,    &
                 0.41353,       0.45589,       0.49806,    &
                 0.53961,       0.58015,       0.61935,    &
                 0.65692,       0.69261,       0.72625,    &
                 0.75773,       0.78698,       0.81398,    &
                 0.83876,       0.86138,       0.88192,    &
                 0.90050,       0.91722,       0.93223,    &
                 0.94565,       0.95762,       0.96827,    &
                 0.97771,       0.98608,       0.99347,  1./

      call mpp_error(NOTE,'Using external_IC::get_nggps_ic which is valid only for data which has been &
                          &horizontally interpolated to the current cubed-sphere grid')
#ifdef INTERNAL_FILE_NML
      read (input_nml_file,external_ic_nml,iostat=ios)
      ierr = check_nml_error(ios,'external_ic_nml')
#else
      unit=open_namelist_file()
      read (unit,external_ic_nml,iostat=ios)
      ierr = check_nml_error(ios,'external_ic_nml')
      call close_file(unit)
#endif

      unit = stdlog()
      call write_version_number ( 'NGGPS_release', 'get_nggps_ics' )
      write(unit, nml=external_ic_nml)

      remap = .true.
      if (ncep_plevels) then
        if (ncep_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use NCEP-defined terrain and pressure &
                              &levels (no vertical remapping)')
          remap = .false.
        else if (filtered_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use externally-generated, filtered terrain &
                              &and NCEP pressure levels (vertical remapping)')
        else if (.not. filtered_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use externally-generated, raw terrain &
                              &and NCEP pressure levels (vertical remapping)')
        endif
      else  ! (.not.ncep_plevels)
        if (ncep_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use NCEP-defined terrain and FV3 pressure &
                              &levels (vertical remapping)')
          remap = .false.
        else if (filtered_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use externally-generated, filtered terrain &
                              &and FV3 pressure levels (vertical remapping)')
        else if (.not. filtered_terrain) then
          call mpp_error(NOTE,'External_IC::get_nggps_ic -  use externally-generated, raw terrain &
                              &and FV3 pressure levels (vertical remapping)')
        endif
      endif

      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie
      js  = Atm(1)%bd%js
      je  = Atm(1)%bd%je
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed
      npz = Atm(1)%npz
      call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)
      ntdiag = ntracers-ntprog

!--- set the 'nestXX' appendix for all files using fms_io
      if (Atm(1)%grid_number > 1) then
         write(gn,'(A4, I2.2)') "nest", Atm(1)%grid_number
      else
         gn = ''
      end if
      call set_filename_appendix(gn)

!--- test for existence of the GFS control file
      if (.not. file_exist('INPUT/'//trim(fn_gfs_ctl), no_domain=.TRUE.)) then
        call mpp_error(FATAL,'==> Error in External_ic::get_nggps_ic: file '//trim(fn_gfs_ctl)//' for NGGPS IC does not exist')
      endif
      call mpp_error(NOTE,'==> External_ic::get_nggps_ic: using control file '//trim(fn_gfs_ctl)//' for NGGPS IC')

!--- read in the number of tracers in the NCEP NGGPS ICs
      call read_data ('INPUT/'//trim(fn_gfs_ctl), 'ntrac', ntrac, no_domain=.TRUE.)

!--- read in ak and bk from the gfs control file using fms_io read_data ---
      allocate (wk2(levp+1,2))
      allocate (ak(levp+1))
      allocate (bk(levp+1))
      call read_data('INPUT/'//trim(fn_gfs_ctl),'vcoord',wk2, no_domain=.TRUE.)
      ak(1:levp+1) = wk2(1:levp+1,1)
      bk(1:levp+1) = wk2(1:levp+1,2)
      deallocate (wk2)

      if (.not. file_exist('INPUT/'//trim(fn_sfc_ics), domain=Atm(1)%domain)) then
        call mpp_error(FATAL,'==> Error in External_ic::get_nggps_ic: tiled file '//trim(fn_sfc_ics)//' for NGGPS IC does not exist')
      endif
      call mpp_error(NOTE,'==> External_ic::get_nggps_ic: using tiled data file '//trim(fn_sfc_ics)//' for NGGPS IC')

      if (.not. file_exist('INPUT/'//trim(fn_gfs_ics), domain=Atm(1)%domain)) then
        call mpp_error(FATAL,'==> Error in External_ic::get_nggps_ic: tiled file '//trim(fn_gfs_ics)//' for NGGPS IC does not exist')
      endif
      call mpp_error(NOTE,'==> External_ic::get_nggps_ic: using tiled data file '//trim(fn_gfs_ics)//' for NGGPS IC')

      allocate (zs(is:ie,js:je))
      allocate (ps(is:ie,js:je))
      allocate (dp(is:ie,js:je,levp))
      allocate (t   (is:ie,js:je,levp))
      allocate (ua  (is:ie,js:je,levp))
      allocate (va  (is:ie,js:je,levp))
      allocate (omga(is:ie,js:je,levp))
      allocate (q (is:ie,js:je,levp,ntrac))
      do n = 1,size(Atm(:))
!--- read in surface temperature (k) and land-frac
        ! surface skin temperature
        id_res = register_restart_field (SFC_restart, fn_sfc_ics, 'tsea', Atm(n)%ts, domain=Atm(n)%domain)

        ! terrain surface height -- (needs to be transformed into phis = zs*grav)
        if (filtered_terrain) then
          id_res = register_restart_field (SFC_restart, fn_sfc_ics, 'orog_filt', Atm(n)%phis, domain=Atm(n)%domain)
        elseif (.not. filtered_terrain) then
          id_res = register_restart_field (SFC_restart, fn_sfc_ics, 'orog_raw', Atm(n)%phis, domain=Atm(n)%domain)
        endif

        if ( Atm(n)%flagstruct%fv_land ) then
          ! stddev
          id_res = register_restart_field (SFC_restart, fn_sfc_ics, 'sdtdev', Atm(n)%sgh, domain=Atm(n)%domain)
          ! land-frac
          id_res = register_restart_field (SFC_restart, fn_sfc_ics, 'land-frac', Atm(n)%oro, domain=Atm(n)%domain)
        endif

        ! NCEP IC surface height -- (needs to be transformed ino phis = zs*grav)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'zs', zs, domain=Atm(n)%domain)
     
        ! surface pressure (Pa)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'ps', ps, domain=Atm(n)%domain)

        ! prognostic delta-p (Pa)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'dp', dp, domain=Atm(n)%domain)

        ! prognostic temperature (k)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 't', t, domain=Atm(n)%domain)

        ! prognostic horizonal wind (m/s)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'u', ua, domain=Atm(n)%domain)

        ! prognostic meridional wind (m/s)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'v', va, domain=Atm(n)%domain)

        ! prognostic vertical velocity 'omga' (Pa/s)
        id_res = register_restart_field (GFS_restart, fn_gfs_ics, 'w', omga, domain=Atm(n)%domain)

        ! prognostic tracers
        do nt = 1, ntrac
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          id_res = register_restart_field (GFS_restart, fn_gfs_ics, trim(tracer_name), q(:,:,:,nt), &
                                           domain=Atm(n)%domain)
        enddo

        ! initialize all tracers to default values prior to being input
        do nt = 1, ntprog
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(:,:,:,n)  )
        enddo
        do nt = ntprog+1, ntracers
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%qdiag(:,:,:,nt)  )
        enddo

        ! read in the restart
        call restore_state (SFC_restart)
        call restore_state (GFS_restart)
        ! free the restart type to be re-used by the nest
!rab        call free_restart_type(SFC_restart)
!rab        call free_restart_type(GFS_restart)

        ! multiply NCEP ICs 'zs' and terrain 'phis' by gravity to be true geopotential
        zs = zs*grav
        Atm(n)%phis = max(Atm(n)%phis*grav, 0.0)
        
        if (ncep_terrain .and. ncep_plevels) then
          ! no vertical remapping necessary, store data into Atm(:)
          ! carve off top layer of atmosphere
          itoa = levp - npz + 1
          Atm(n)%ptop =ak(itoa)
          Atm(n)%ak(1:npz+1) = ak(itoa:levp+1)
          Atm(n)%bk(1:npz+1) = bk(itoa:levp+1)
          Atm(n)%ps  (is:ie,js:je) = ps(is:ie,js:je)
          Atm(n)%phis(is:ie,js:je) = max(zs(is:ie,js:je),0.0)
          Atm(n)%delp(is:ie,js:je,1:npz) = dp  (is:ie,js:je,itoa:levp)
          Atm(n)%pt  (is:ie,js:je,1:npz) = t   (is:ie,js:je,itoa:levp)
          Atm(n)%ua  (is:ie,js:je,1:npz) = ua  (is:ie,js:je,itoa:levp)
          Atm(n)%va  (is:ie,js:je,1:npz) = va  (is:ie,js:je,itoa:levp)
          Atm(n)%omga(is:ie,js:je,1:npz) = omga(is:ie,js:je,itoa:levp)
          Atm(n)%q   (is:ie,js:je,1:npz,1:ntrac) = q(is:ie,js:je,itoa:levp,1:ntrac)
 
          call get_w_from_omga(is, js, npz, ak(itoa:levp+1), bk(itoa:levp+1), ps(:,:), &
                               dp(:,:,itoa:levp), t(:,:,itoa:levp), omga(:,:,itoa:levp),Atm(n))

          ! populate the haloes of Atm(:)%phis
          call mpp_update_domains( Atm(n)%phis, Atm(n)%domain )
          ! map the A-grid winds onto the D-grid winds
          call cubed_a2d (Atm(n)%npx, Atm(n)%npy, npz, Atm(n)%ua, Atm(n)%va, Atm(n)%u, Atm(n)%v, &
                          Atm(n)%gridstruct, Atm(n)%domain, Atm(n)%bd )
        else
          ! set the pressure levels and ptop to be used
          if (ncep_plevels) then
            itoa = levp - npz + 1
            Atm(n)%ptop = ak(itoa)
            Atm(n)%ak(1:npz+1) = ak(itoa:levp+1)
            Atm(n)%bk(1:npz+1) = bk(itoa:levp+1)
          else
            Atm(n)%ptop = ak_sj(1)
            Atm(n)%ak(:) = ak_sj(:)
            Atm(n)%bk(:) = bk_sj(:)
          endif
          ! call vertical remapping algorithms
          call remap_scalar_nggps(is, js, levp, npz, ntprog, ntrac, ak(:), bk(:), ps, zs, &
                                  t(:,:,:), q(:,:,:,:),omga(:,:,:),Atm(n))
          call remap_winds (is, js, levp, npz, ak(:), bk(:), ps, ua(:,:,:), va(:,:,:), Atm(n))
        endif
      enddo

      Atm(1)%flagstruct%make_nh = .false.

      deallocate (ak)
      deallocate (bk)
      deallocate (zs)
      deallocate (ps)
      deallocate (dp)
      deallocate (t )
      deallocate (ua)
      deallocate (va)
      deallocate (q )

  end subroutine get_nggps_ic

  

  subroutine get_ncep_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq
! local:
      character(len=128) :: fname
      real(kind=4), allocatable:: wk1(:), wk2(:,:), wk3(:,:,:)
      real, allocatable:: tp(:,:,:), qp(:,:,:)
      real, allocatable:: ua(:,:,:), va(:,:,:)
      real, allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      real:: s2c(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je,4)
      integer, dimension(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je):: id1, id2, jdc
      real psc(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je)
      real gzc(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je)
      real tmean
      integer:: i, j, k, im, jm, km, npz, npt
      integer:: i1, i2, j1, ncid
      integer:: jbeg, jend
      integer tsize(3) 
      logical:: read_ts = .true.
      logical:: land_ts = .false.
      logical:: found
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = Atm(1)%bd%is
      ie  = Atm(1)%bd%ie
      js  = Atm(1)%bd%js
      je  = Atm(1)%bd%je
      isd = Atm(1)%bd%isd
      ied = Atm(1)%bd%ied
      jsd = Atm(1)%bd%jsd
      jed = Atm(1)%bd%jed

      deg2rad = pi/180.

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
! SJL: 20110716
!      Atm(1)%q = 0.

      fname = Atm(1)%flagstruct%res_latlon_dynamics

      if( file_exist(fname) ) then
          call open_ncfile( fname, ncid )        ! open the file
          call get_ncdim1( ncid, 'lon', tsize(1) )
          call get_ncdim1( ncid, 'lat', tsize(2) )
          call get_ncdim1( ncid, 'lev', tsize(3) )

          im = tsize(1); jm = tsize(2); km = tsize(3)

          if(is_master())  write(*,*) fname
          if(is_master())  write(*,*) ' NCEP IC dimensions:', tsize

          allocate (  lon(im) )
          allocate (  lat(jm) )
 
          call _GET_VAR1 (ncid, 'lon', im, lon )
          call _GET_VAR1 (ncid, 'lat', jm, lat )

! Convert to radian
          do i=1,im
             lon(i) = lon(i) * deg2rad  ! lon(1) = 0.
          enddo
          do j=1,jm
             lat(j) = lat(j) * deg2rad
          enddo

          allocate ( ak0(km+1) )
          allocate ( bk0(km+1) )
          call _GET_VAR1 (ncid, 'hyai', km+1, ak0, found )
          if ( .not. found )  ak0(:) = 0.

          call _GET_VAR1 (ncid, 'hybi', km+1, bk0 )

! Note: definition of NCEP hybrid is p(k) = a(k)*1.E5 + b(k)*ps
          ak0(:) = ak0(:) * 1.E5

! Limiter to prevent NAN at top during remapping
          if ( bk0(1) < 1.E-9 ) ak0(1) = max(1.e-9, ak0(1))
      else
          call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' for NCEP IC does not exist')
      endif

! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      call remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c , Atm(1)%gridstruct%agrid, Atm(1)%bd)

! Find bounding latitudes:
      jbeg = jm-1;         jend = 2
      do j=js,je
         do i=is,ie
              j1 = jdc(i,j)
            jbeg = min(jbeg, j1) 
            jend = max(jend, j1+1)
         enddo
      enddo

! remap surface pressure and height:

      allocate ( wk2(im,jbeg:jend) )
      call get_var3_r4( ncid, 'PS', 1,im, jbeg,jend, 1,1, wk2 )

      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            psc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                       s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      call get_var3_r4( ncid, 'PHIS', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            gzc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                       s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      deallocate ( wk2 )
      allocate ( wk2(im,jm) )

      if ( read_ts ) then       ! read skin temperature; could be used for SST

        call get_var2_real( ncid, 'TS', im, jm, wk2 )

        if ( .not. land_ts ) then
           allocate ( wk1(im) )

           do j=1,jm
! Read NCEP ORO (1; land; 0: ocean; 2: sea_ice)
              call get_var3_r4( ncid, 'ORO', 1,im, j,j, 1,1, wk1 )
              tmean = 0.
              npt = 0
              do i=1,im
                 if( abs(wk1(i)-1.) > 0.99 ) then   ! ocean or sea ice
                     tmean = tmean + wk2(i,j)
                     npt = npt + 1
                 endif
              enddo
!------------------------------------------------------
! Replace TS over interior land with zonal mean SST/Ice
!------------------------------------------------------
              if ( npt /= 0 ) then
                   tmean= tmean / real(npt)
                   do i=1,im
                      if( abs(wk1(i)-1.) <= 0.99 ) then  ! Land points
                          if ( i==1 ) then
                               i1 = im;     i2 = 2
                          elseif ( i==im ) then
                               i1 = im-1;   i2 = 1
                          else
                               i1 = i-1;    i2 = i+1
                          endif
                          if ( abs(wk1(i2)-1.)>0.99 ) then     ! east side has priority
                               wk2(i,j) = wk2(i2,j)
                          elseif ( abs(wk1(i1)-1.)>0.99 ) then ! west side
                               wk2(i,j) = wk2(i1,j)
                          else
                               wk2(i,j) = tmean
                          endif
                      endif
                   enddo
              endif
           enddo   ! j-loop
           deallocate ( wk1 )
        endif   !(.not.land_ts)

        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            Atm(1)%ts(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                             s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
          enddo
        enddo
        call prt_maxmin('SST_model', Atm(1)%ts, is, ie, js, je, 0, 1, 1.)

! Perform interp to FMS SST format/grid
#ifndef DYCORE_SOLO
        call ncep2fms(im, jm, lon, lat, wk2)
        if( is_master() ) then
          write(*,*) 'External_ic_mod: i_sst=', i_sst, ' j_sst=', j_sst
          call pmaxmin( 'SST_ncep_fms',  sst_ncep, i_sst, j_sst, 1.)
        endif
#endif
      endif  !(read_ts)

      deallocate ( wk2 )

! Read in temperature:
      allocate ( wk3(1:im,jbeg:jend, 1:km) )
      call get_var3_r4( ncid, 'T', 1,im, jbeg,jend, 1,km, wk3 )

      allocate (  tp(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
        enddo
      enddo

! Read in tracers: only sphum at this point
      call get_var3_r4( ncid, 'Q', 1,im, jbeg,jend, 1,km, wk3 )

      allocate ( qp(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            qp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
          enddo
        enddo
      enddo

      call remap_scalar(im, jm, km, npz, nq, nq, ak0, bk0, psc, gzc, tp, qp, Atm(1))
      deallocate ( tp )
      deallocate ( qp )

! Winds:
      call get_var3_r4( ncid, 'U', 1,im, jbeg,jend, 1,km, wk3 )

      allocate ( ua(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ua(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
          enddo
        enddo
      enddo

      call get_var3_r4( ncid, 'V', 1,im, jbeg,jend, 1,km, wk3 )
      call close_ncfile ( ncid )

      allocate ( va(is:ie,js:je,km) )
      do k=1,km
        do j=js,je
          do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            va(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                        s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
          enddo
        enddo
      enddo
      deallocate ( wk3 )

      call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm(1))

      deallocate ( ua )
      deallocate ( va )

      deallocate ( ak0 )
      deallocate ( bk0 )
      deallocate ( lat )
      deallocate ( lon )

  
  end subroutine get_ncep_ic



  subroutine get_fv_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq

      character(len=128) :: fname, tracer_name
      real, allocatable:: ps0(:,:), gz0(:,:), u0(:,:,:), v0(:,:,:), t0(:,:,:), dp0(:,:,:), q0(:,:,:,:)
      real, allocatable:: ua(:,:,:), va(:,:,:)
      real, allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      integer :: i, j, k, im, jm, km, npz, tr_ind
      integer tsize(3)
!     integer sphum, liq_wat, ice_wat, cld_amt       ! GFDL AM2 physics
      logical found

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read in lat-lon FV core restart file
      fname = Atm(1)%flagstruct%res_latlon_dynamics

      if( file_exist(fname) ) then
          call field_size(fname, 'T', tsize, field_found=found)
          if(is_master()) write(*,*) 'Using lat-lon FV restart:', fname 

          if ( found ) then
               im = tsize(1); jm = tsize(2); km = tsize(3)
               if(is_master())  write(*,*) 'External IC dimensions:', tsize
          else
               call mpp_error(FATAL,'==> Error in get_external_ic: field not found')
          endif

! Define the lat-lon coordinate:
          allocate (  lon(im) )
          allocate (  lat(jm) )

          do i=1,im
             lon(i) = (0.5 + real(i-1)) * 2.*pi/real(im)
          enddo

          do j=1,jm
             lat(j) = -0.5*pi + real(j-1)*pi/real(jm-1)   ! SP to NP 
          enddo
 
          allocate ( ak0(1:km+1) )
          allocate ( bk0(1:km+1) )
          allocate ( ps0(1:im,1:jm) )
          allocate ( gz0(1:im,1:jm) )
          allocate (  u0(1:im,1:jm,1:km) )
          allocate (  v0(1:im,1:jm,1:km) )
          allocate (  t0(1:im,1:jm,1:km) )
          allocate ( dp0(1:im,1:jm,1:km) )

          call read_data (fname, 'ak', ak0)
          call read_data (fname, 'bk', bk0)
          call read_data (fname, 'Surface_geopotential', gz0)
          call read_data (fname, 'U',     u0)
          call read_data (fname, 'V',     v0)
          call read_data (fname, 'T',     t0)
          call read_data (fname, 'DELP', dp0)

! Share the load
          if(is_master()) call pmaxmin( 'ZS_data', gz0, im,    jm, 1./grav)
          if(mpp_pe()==1) call pmaxmin( 'U_data',   u0, im*jm, km, 1.)
          if(mpp_pe()==1) call pmaxmin( 'V_data',   v0, im*jm, km, 1.)
          if(mpp_pe()==2) call pmaxmin( 'T_data',   t0, im*jm, km, 1.)
          if(mpp_pe()==3) call pmaxmin( 'DEL-P',   dp0, im*jm, km, 0.01)


      else
          call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' for dynamics does not exist')
      endif

! Read in tracers: only AM2 "physics tracers" at this point
      fname = Atm(1)%flagstruct%res_latlon_tracers

      if( file_exist(fname) ) then
          if(is_master()) write(*,*) 'Using lat-lon tracer restart:', fname 

          allocate ( q0(im,jm,km,Atm(1)%ncnst) )
          q0 = 0.

          do tr_ind = 1, nq
            call get_tracer_names(MODEL_ATMOS, tr_ind, tracer_name)
            if (field_exist(fname,tracer_name)) then
               call read_data(fname, tracer_name, q0(1:im,1:jm,1:km,tr_ind))
               call mpp_error(NOTE,'==>  Have read tracer '//trim(tracer_name)//' from '//trim(fname))
               cycle
            endif
          enddo
      else
          call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' for tracers does not exist')
      endif

! D to A transform on lat-lon grid:
      allocate (  ua(im,jm,km) )
      allocate (  va(im,jm,km) )

      call d2a3d(u0, v0,  ua,  va, im, jm, km, lon)

      deallocate ( u0 ) 
      deallocate ( v0 ) 

      if(mpp_pe()==4) call pmaxmin( 'UA', ua, im*jm, km, 1.)
      if(mpp_pe()==4) call pmaxmin( 'VA', va, im*jm, km, 1.)

      do j=1,jm
         do i=1,im
            ps0(i,j) = ak0(1)
         enddo
      enddo

      do k=1,km
         do j=1,jm
            do i=1,im
               ps0(i,j) = ps0(i,j) + dp0(i,j,k)
            enddo
         enddo
      enddo

  if (is_master()) call pmaxmin( 'PS_data (mb)', ps0, im, jm, 0.01)

! Horizontal interpolation to the cubed sphere grid center
! remap vertically with terrain adjustment

      call remap_xyz( im, 1, jm, jm, km, npz, nq, Atm(1)%ncnst, lon, lat, ak0, bk0,   &
                      ps0,  gz0, ua, va, t0, q0, Atm(1) )

      deallocate ( ak0 ) 
      deallocate ( bk0 ) 
      deallocate ( ps0 ) 
      deallocate ( gz0 ) 
      deallocate ( t0 ) 
      deallocate ( q0 ) 
      deallocate ( dp0 ) 
      deallocate ( ua ) 
      deallocate ( va ) 
      deallocate ( lat ) 
      deallocate ( lon ) 

  end subroutine get_fv_ic


#ifndef DYCORE_SOLO
 subroutine ncep2fms(im, jm, lon, lat, wk)

  integer, intent(in):: im, jm
  real,    intent(in):: lon(im), lat(jm)
  real(kind=4),    intent(in):: wk(im,jm)
! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  real:: delx, dely
  real:: xc, yc    ! "data" location
  real:: c1, c2, c3, c4
  integer i,j, i1, i2, jc, i0, j0, it, jt

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to "FMS" 1x1 SST data grid
! lon: 0.5, 1.5, ..., 359.5
! lat: -89.5, -88.5, ... , 88.5, 89.5

  delx = 360./real(i_sst)
  dely = 180./real(j_sst)

  jt = 1
  do 5000 j=1,j_sst

     yc = (-90. + dely * (0.5+real(j-1)))  * deg2rad
     if ( yc<lat(1) ) then
            jc = 1
            b1 = 0.
     elseif ( yc>lat(jm) ) then
            jc = jm-1
            b1 = 1.
     else
          do j0=jt,jm-1
          if ( yc>=lat(j0) .and. yc<=lat(j0+1) ) then
               jc = j0
               jt = j0
               b1 = (yc-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
     endif
222  continue
     it = 1

     do i=1,i_sst
        xc = delx * (0.5+real(i-1)) * deg2rad
       if ( xc>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (xc-lon(im)) * rdlon(im)
       elseif ( xc<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (xc+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=it,im-1
            if ( xc>=lon(i0) .and. xc<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               it = i0
               a1 = (xc-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
       endif

       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1
! Interpolated surface pressure
       sst_ncep(i,j) = c1*wk(i1,jc  ) + c2*wk(i2,jc  ) +    &
                       c3*wk(i2,jc+1) + c4*wk(i1,jc+1)
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine ncep2fms
#endif



 subroutine remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c, agrid, bd )

   type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: im, jm
  real,    intent(in):: lon(im), lat(jm)
  real,    intent(out):: s2c(bd%is:bd%ie,bd%js:bd%je,4)
  integer, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je):: id1, id2, jdc
  real,    intent(in):: agrid(bd%isd:bd%ied,bd%jsd:bd%jed,2)
! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
       endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef


 subroutine remap_scalar(im, jm, km, npz, nq, ncnst, ak0, bk0, psc, gzc, ta, qa, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc, gzc
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ta
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
! local:
  real, dimension(Atm%bd%is:Atm%bd%ie,km):: tp
  real, dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real, dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real pt0(km), gz(km+1), pk0(km+1)
  real qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real pst, p1, p2, alpha, rdg
  integer i,j,k, iq
  integer  sphum, o3mr, clwmr
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed


  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
  o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')
  clwmr   = get_tracer_index(MODEL_ATMOS, 'clwmr')

  if (mpp_pe()==1) then
    print *, 'sphum = ', sphum
    print *, 'clwmr = ', clwmr
    print *, ' o3mr = ', o3mr
    print *, 'ncnst = ', ncnst
  endif

  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif

!rab  if ( o3mr/=ncnst ) then
!rab       call mpp_error(FATAL,'number of q species is not match')
!rab  endif

  do 5000 j=js,je

     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo

    if ( T_is_Tv ) then
! The "T" field in NCEP analysis is actually virtual temperature (Larry H. post processing)
! BEFORE 20051201
       do k=1,km
          tp(i,k) = ta(i,j,k)
       enddo
    else
       do k=1,km
          tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
       enddo
    endif
! Tracers:

       do k=1,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
          pn0(i,k) = log(pe0(i,k))
            pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc(i,j)
       Atm%phis(i,j) = gzc(i,j)
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j)
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop


     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

!---------------
! map shpum, o3mr, clwmr tracers
!----------------
      do iq=1,ncnst
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         if (iq==sphum .and. Atm%flagstruct%ncep_ic ) then
             p1 = 200.E2
             p2 =  75.E2
! Blend model sphum with NCEP data
             do k=1,npz
                do i=is,ie
                   pst = 0.5*(pe1(i,k)+pe1(i,k+1))
                   if ( pst > p1 ) then
                        Atm%q(i,j,k,iq) = qn1(i,k)
                   elseif( pst > p2 ) then            ! p2 < pst < p1
                        alpha = (pst-p2)/(p1-p2)
                        Atm%q(i,j,k,1) = qn1(i,k)*alpha + Atm%q(i,j,k,1)*(1.-alpha)
                   endif
                enddo
             enddo
         else
             do k=1,npz
                do i=is,ie
                   Atm%q(i,j,k,iq) = qn1(i,k)
                enddo
             enddo
         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

      if ( .not. Atm%flagstruct%hydrostatic .and. Atm%flagstruct%ncep_ic ) then
! Replace delz with NCEP hydrostatic state
         rdg = -rdgas / grav
         do k=1,npz
            do i=is,ie
               atm%delz(i,j,k) = rdg*qn1(i,k)*(pn1(i,k+1)-pn1(i,k))
            enddo
         enddo
      endif

5000 continue

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01)

  if (is_master()) write(*,*) 'done remap_scalar'

 end subroutine remap_scalar


 subroutine remap_scalar_nggps(im, jm, km, npz, nq, ncnst, ak0, bk0, psc, gzc, ta, qa, omga, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc, gzc
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ta
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: omga
! local:
  real, dimension(Atm%bd%is:Atm%bd%ie,km):: tp, omgap
  real, dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real, dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real pt0(km), gz(km+1), pk0(km+1)
  real qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real pst, p1, p2, alpha, rdg, ak0p
  integer i,j,k, iq
  integer  sphum, o3mr, clwmr
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed


  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
  clwmr   = get_tracer_index(MODEL_ATMOS, 'clwmr')
  o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')

  if (mpp_pe()==1) then
    print *, 'sphum = ', sphum
    print *, 'clwmr = ', clwmr
    print *, ' o3mr = ', o3mr
    print *, 'ncnst = ', ncnst
  endif

  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif

!rab  if ( o3mwr/=ncnst ) then
!rab       call mpp_error(FATAL,'number of q species is not match')
!rab  endif

  do 5000 j=js,je

     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo

       do k=1,km
          omgap(i,k) = omga(i,j,k)
       enddo

       if ( T_is_Tv ) then
! The "T" field in NCEP analysis is actually virtual temperature (Larry H. post processing)
! BEFORE 20051201
         do k=1,km
            tp(i,k) = ta(i,j,k)
         enddo
      else
         do k=1,km
            tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
         enddo
      endif

! Tracers:

      if ( bk0(1) < 1.E-9 ) then 
         ak0p = max(1.e-9, ak0(1))
         pe0(i,1) = ak0p + bk0(1)*psc(i,j)
         pn0(i,1) = log(pe0(i,1))
           pk0(1) = pe0(i,1)**kappa
      else
         pe0(i,1) = ak0(1) + bk0(1)*psc(i,j)
         pn0(i,1) = log(pe0(i,1))
           pk0(1) = pe0(i,1)**kappa
      endif
      
      do k=2,km+1
         pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
         pn0(i,k) = log(pe0(i,k))
           pk0(k) = pe0(i,k)**kappa
      enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc(i,j)
       Atm%phis(i,j) = gzc(i,j)
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j)
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop


     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

!---------------
! map shpum, o3mr, clwmr tracers
!----------------
      do iq=1,ncnst
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         if (iq==sphum .and. Atm%flagstruct%nggps_ic ) then
             p1 = 200.E2
             p2 =  75.E2
! Blend model sphum with NCEP data
             do k=1,npz
                do i=is,ie
                   pst = 0.5*(pe1(i,k)+pe1(i,k+1))
                   if ( pst > p1 ) then
                        Atm%q(i,j,k,iq) = qn1(i,k)
                   elseif( pst > p2 ) then            ! p2 < pst < p1
                        alpha = (pst-p2)/(p1-p2)
                        Atm%q(i,j,k,1) = qn1(i,k)*alpha + Atm%q(i,j,k,1)*(1.-alpha)
                   endif
                enddo
             enddo
         else
             do k=1,npz
                do i=is,ie
                   Atm%q(i,j,k,iq) = qn1(i,k)
                enddo
             enddo
         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

      if ( .not. Atm%flagstruct%hydrostatic .and. Atm%flagstruct%nggps_ic ) then
! Replace delz with NCEP hydrostatic state
         rdg = -rdgas / grav
         do k=1,npz
            do i=is,ie
               atm%delz(i,j,k) = rdg*qn1(i,k)*(pn1(i,k+1)-pn1(i,k))
            enddo
         enddo
      endif

!-------------------------------------------------------------
! map omega
!-------------------------------------------------------------
      call mappm(km, pe0, omgap, npz, pe1, qn1, is,ie, 0, 11, Atm%ptop)
      do k=1,npz
         do i=is,ie
            atm%w(i,j,k) = qn1(i,k)/atm%delp(i,j,k)*atm%delz(i,j,k)
         enddo
      enddo

5000 continue

!  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01)
!  call prt_maxmin('W_model in remap_scalar_nggps', Atm%w, is, ie, js, je, ng, npz, 1. )

  if (is_master()) write(*,*) 'done remap_scalar_nggps'

 end subroutine remap_scalar_nggps


 subroutine get_w_from_omga(im, jm, npz, ak0, bk0, psc, dp, ta, omga, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, npz
  real,    intent(in):: ak0(npz+1), bk0(npz+1)
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,npz):: ta, dp, omga
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe0, pn0
  real rdg
  integer i,j,k
  integer :: is,  ie,  js,  je

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je


  do j = js, je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(pe0(i,1))
     enddo

     do k=2,npz+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
           pn0(i,k) = log(pe0(i,k))
        enddo
     enddo

     rdg = -rdgas / grav
     do k= 1, npz
        do i = is, ie
           atm%delz(i,j,k) = rdg*ta(i,j,k)*(pn0(i,k+1)-pn0(i,k))
        enddo
     enddo

     do k = 1, npz
        do i = is, ie
           atm%w(i,j,k) = omga(i,j,k)/dp(i,j,k)*atm%delz(i,j,k)
        enddo
     enddo

  enddo

!  call prt_maxmin('W_model', Atm%w, is, ie, js, je, ng, npz, 1.)

  if (is_master()) write(*,*) 'done get_w_from_omga'

 end subroutine get_w_from_omga



 subroutine remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in):: psc(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je)
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ua, va
! local:
  real, dimension(Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed,npz):: ut, vt   ! winds
  real, dimension(Atm%bd%is:Atm%bd%ie, km+1):: pe0
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1
  real, dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  integer i,j,k

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  do 5000 j=js,je

     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
        enddo
     enddo

     do k=1,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
       enddo
     enddo

!------
! map u
!------
      call mappm(km, pe0, ua(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm%ptop)
      do k=1,npz
         do i=is,ie
            ut(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, va(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm%ptop)
      do k=1,npz
         do i=is,ie
            vt(i,j,k) = qn1(i,k)
         enddo
      enddo

5000 continue

  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1.)
  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1.)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
  call cubed_a2d(Atm%npx, Atm%npy, npz, ut, vt, Atm%u, Atm%v, Atm%gridstruct, Atm%domain, Atm%bd )

  if (is_master()) write(*,*) 'done remap_winds'

 end subroutine remap_winds


 subroutine remap_wz(im, jm, km, npz, mg, ak0, bk0, psc, wa, wz, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz
  integer, intent(in):: mg     ! mg = 0 for delz; mg=3 for w
  real,    intent(in):: ak0(km+1), bk0(km+1)
  real,    intent(in):: psc(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je)
  real,    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: wa
  real,   intent(out):: wz(Atm%bd%is-mg:Atm%bd%ie+mg,Atm%bd%js-mg:Atm%bd%je+mg,npz)
! local:
  real, dimension(Atm%bd%is:Atm%bd%ie, km+1):: pe0
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1
  real, dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  integer i,j,k
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  do 5000 j=js,je

     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
        enddo
     enddo

     do k=1,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
       enddo
     enddo

!------
! map w
!------
      call mappm(km, pe0, wa(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm%ptop)
      do k=1,npz
         do i=is,ie
            wz(i,j,k) = qn1(i,k)
         enddo
      enddo

5000 continue

! call prt_maxmin('WZ', wz, is, ie, js, je, mg, npz, 1.)
! if (is_master()) write(*,*) 'done remap_wz'

 end subroutine remap_wz



  subroutine remap_xyz( im, jbeg, jend, jm, km, npz, nq, ncnst, lon, lat, ak0, bk0, ps0, gz0,   &
                        ua, va, ta, qa, Atm )

  type(fv_atmos_type), intent(inout), target :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  integer, intent(in):: jbeg, jend
  real,    intent(in):: lon(im), lat(jm), ak0(km+1), bk0(km+1)
  real,    intent(in):: gz0(im,jbeg:jend), ps0(im,jbeg:jend)
  real,    intent(in), dimension(im,jbeg:jend,km):: ua, va, ta
  real,    intent(in), dimension(im,jbeg:jend,km,ncnst):: qa

  real, pointer, dimension(:,:,:) :: agrid

! local:
  real, dimension(Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed,npz):: ut, vt   ! winds 
  real, dimension(Atm%bd%is:Atm%bd%ie,km):: up, vp, tp
  real, dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real pt0(km), gz(km+1), pk0(km+1)
  real qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real, dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real, dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1, c1, c2, c3, c4
  real:: gzc, psc, pst
  integer i,j,k, i1, i2, jc, i0, j0, iq
! integer  sphum, liq_wat, ice_wat, cld_amt
  integer  sphum
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  !!NOTE: Only Atm is used in this routine.
  agrid => Atm%gridstruct%agrid

  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
! liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
! ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
! cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

   if ( sphum/=1 ) then
        call mpp_error(FATAL,'SPHUM must be 1st tracer')
   endif

  pk0(1) = ak0(1)**kappa 

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(ak0(1))
     enddo

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif

111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

#ifndef DEBUG_REMAP
       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) i,j,a1, b1
       endif
#endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1

! Interpolated surface pressure
       psc = c1*ps0(i1,jc  ) + c2*ps0(i2,jc  ) +    &
             c3*ps0(i2,jc+1) + c4*ps0(i1,jc+1)

! Interpolated surface geopotential
       gzc = c1*gz0(i1,jc  ) + c2*gz0(i2,jc  ) +    &
             c3*gz0(i2,jc+1) + c4*gz0(i1,jc+1)

! 3D fields:
       do iq=1,ncnst
!          if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
          do k=1,km
             qp(i,k,iq) = c1*qa(i1,jc,  k,iq) + c2*qa(i2,jc,  k,iq) +  &
                          c3*qa(i2,jc+1,k,iq) + c4*qa(i1,jc+1,k,iq)
          enddo
!          endif
       enddo

       do k=1,km
          up(i,k) = c1*ua(i1,jc,  k) + c2*ua(i2,jc,  k) +  &
                    c3*ua(i2,jc+1,k) + c4*ua(i1,jc+1,k)
          vp(i,k) = c1*va(i1,jc,  k) + c2*va(i2,jc,  k) +  &
                    c3*va(i2,jc+1,k) + c4*va(i1,jc+1,k)
          tp(i,k) = c1*ta(i1,jc,  k) + c2*ta(i2,jc,  k) +  &
                    c3*ta(i2,jc+1,k) + c4*ta(i1,jc+1,k)
! Virtual effect:
          tp(i,k) = tp(i,k)*(1.+zvir*qp(i,k,sphum))
       enddo
! Tracers:

       do k=2,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc
          pn0(i,k) = log(pe0(i,k))
          pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc
       Atm%phis(i,j) = gzc
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc 
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k)) 
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop
 

! * Compute delp from ps
     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo
 
! Use kord=9 for winds; kord=11 for tracers
!------
! map u
!------
      call mappm(km, pe0, up, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            ut(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, vp, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            vt(i,j,k) = qn1(i,k)
         enddo
      enddo

!---------------
! map tracers
!----------------
      do iq=1,ncnst
! Note: AM2 physics tracers only
!         if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         do k=1,npz
            do i=is,ie
               Atm%q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
!         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

5000 continue

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01)
  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1.)
  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1.)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
  call cubed_a2d(Atm%npx, Atm%npy, npz, ut, vt, Atm%u, Atm%v, Atm%gridstruct, Atm%domain, Atm%bd )

  if (is_master()) write(*,*) 'done remap_xyz'

 end subroutine remap_xyz


 subroutine cubed_a2d( npx, npy, npz, ua, va, u, v, gridstruct, fv_domain, bd )

! Purpose; Transform wind on A grid to D grid

  use mpp_domains_mod,    only: mpp_update_domains

  type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: npx, npy, npz
  real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: ua, va
  real, intent(out):: u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
  real, intent(out):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: fv_domain
! local:
  real v3(3,bd%is-1:bd%ie+1,bd%js-1:bd%je+1)
  real ue(3,bd%is-1:bd%ie+1,bd%js:bd%je+1)    ! 3D winds at edges
  real ve(3,bd%is:bd%ie+1,bd%js-1:bd%je+1)    ! 3D winds at edges
  real, dimension(bd%is:bd%ie):: ut1, ut2, ut3
  real, dimension(bd%js:bd%je):: vt1, vt2, vt3
  integer i, j, k, im2, jm2

  real(kind=R_GRID), pointer, dimension(:,:,:)   :: vlon, vlat
  real(kind=R_GRID), pointer, dimension(:)       :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n
  real(kind=R_GRID), pointer, dimension(:,:,:,:) :: ew, es

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed

  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  edge_vect_w => gridstruct%edge_vect_w
  edge_vect_e => gridstruct%edge_vect_e
  edge_vect_s => gridstruct%edge_vect_s
  edge_vect_n => gridstruct%edge_vect_n
  
  ew => gridstruct%ew
  es => gridstruct%es

  call mpp_update_domains(ua, fv_domain, complete=.false.)
  call mpp_update_domains(va, fv_domain, complete=.true.)

    im2 = (npx-1)/2
    jm2 = (npy-1)/2

    do k=1, npz
! Compute 3D wind on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(1,i,j) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
             v3(2,i,j) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
             v3(3,i,j) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(1,i,j) = 0.5*(v3(1,i,j-1) + v3(1,i,j))
             ue(2,i,j) = 0.5*(v3(2,i,j-1) + v3(2,i,j))
             ue(3,i,j) = 0.5*(v3(3,i,j-1) + v3(3,i,j))
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(1,i,j) = 0.5*(v3(1,i-1,j) + v3(1,i,j))
             ve(2,i,j) = 0.5*(v3(2,i-1,j) + v3(2,i,j))
             ve(3,i,j) = 0.5*(v3(3,i-1,j) + v3(3,i,j))
          enddo
       enddo

! --- E_W edges (for v-wind):
     if (.not. gridstruct%nested) then
     if ( is==1) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(1,i,j-1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j-1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j-1)+(1.-edge_vect_w(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_w(j)*ve(1,i,j+1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j+1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j+1)+(1.-edge_vect_w(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(1,i,j-1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j-1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j-1)+(1.-edge_vect_e(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_e(j)*ve(1,i,j+1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j+1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j+1)+(1.-edge_vect_e(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(1,i-1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i-1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i-1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_s(i)*ue(1,i+1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i+1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i+1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(1,i-1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i-1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i-1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_n(i)*ue(1,i+1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i+1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i+1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     endif ! .not. nested

     do j=js,je+1
        do i=is,ie
           u(i,j,k) =  ue(1,i,j)*es(1,i,j,1) +  &
                       ue(2,i,j)*es(2,i,j,1) +  &
                       ue(3,i,j)*es(3,i,j,1)
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) = ve(1,i,j)*ew(1,i,j,2) +  &
                      ve(2,i,j)*ew(2,i,j,2) +  &
                      ve(3,i,j)*ew(3,i,j,2)
        enddo
     enddo
 
   enddo         ! k-loop

 end subroutine cubed_a2d



 subroutine d2a3d(u, v,  ua,   va,  im,  jm, km, lon)
      integer, intent(in):: im, jm, km           ! Dimensions
      real, intent(in ) :: lon(im)
      real, intent(in ), dimension(im,jm,km):: u, v
      real, intent(out), dimension(im,jm,km):: ua, va
! local
      real :: coslon(im),sinlon(im)    ! Sine and cosine in longitude
      integer i, j, k
      integer imh
      real un, vn, us, vs

      integer :: ks, ke

      imh = im/2

      do i=1,im
         sinlon(i) = sin(lon(i))
         coslon(i) = cos(lon(i))
      enddo

      do k=1,km
         do j=2,jm-1
            do i=1,im
               ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
            enddo
         enddo

         do j=2,jm-1
            do i=1,im-1
               va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(im,j,k) = 0.5*(v(im,j,k) + v(1,j,k))
         enddo

! Projection at SP
             us = 0.
             vs = 0.
             do i=1,imh
                us = us + (ua(i+imh,2,k)-ua(i,2,k))*sinlon(i)      &
                     + (va(i,2,k)-va(i+imh,2,k))*coslon(i)
                vs = vs + (ua(i+imh,2,k)-ua(i,2,k))*coslon(i)      &
                     + (va(i+imh,2,k)-va(i,2,k))*sinlon(i)
             enddo
             us = us/im
             vs = vs/im
             do i=1,imh
                ua(i,1,k)   = -us*sinlon(i) - vs*coslon(i)
                va(i,1,k)   =  us*coslon(i) - vs*sinlon(i)
                ua(i+imh,1,k)   = -ua(i,1,k)
                va(i+imh,1,k)   = -va(i,1,k)
             enddo

! Projection at NP
             un = 0.
             vn = 0.
             do i=1,imh
                un = un + (ua(i+imh,jm-1,k)-ua(i,jm-1,k))*sinlon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*coslon(i)
                vn = vn + (ua(i,jm-1,k)-ua(i+imh,jm-1,k))*coslon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*sinlon(i)
             enddo

             un = un/im
             vn = vn/im
             do i=1,imh
                ua(i,jm,k) = -un*sinlon(i) + vn*coslon(i)
                va(i,jm,k) = -un*coslon(i) - vn*sinlon(i)
                ua(i+imh,jm,k) = -ua(i,jm,k)
                va(i+imh,jm,k) = -va(i,jm,k)
             enddo
      enddo

  end subroutine d2a3d



  subroutine pmaxmin( qname, a, im, jm, fac )

      integer, intent(in):: im, jm
      character(len=*) :: qname
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jm
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jm
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin



 end module external_ic_mod

