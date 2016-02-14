      SUBROUTINE GCYCLE(ME,len,lsoil,IDATE,FHOUR,FHCYC,XLON,XLAT, &
     &                  Sfc_props,Cld_props,Tbd_data,ialb,use_ufo,nst_anl)
!
!RAB      USE MACHINE
      USE PHYSCONS, PI => con_PI
      use nuopc_physics, only: sfc_properties, cloud_properties, tbd_ddt
      use fms_io_mod,    only: open_namelist_file, close_file
      implicit none
!
      TYPE(sfc_properties),   intent(inout) :: Sfc_props
      TYPE(cloud_properties), intent(inout) :: Cld_props
      TYPE(Tbd_ddt),          intent(inout) :: Tbd_data
!
      INTEGER, intent(in) :: len, lsoil, IALB
      INTEGER, intent(in) :: ME, IDATE(4)
      logical, intent(in) :: use_ufo, nst_anl
      real(kind=kind_phys), intent(in) :: fhour, fhcyc
      real(kind=kind_phys), intent(in) :: XLON(len), XLAT(len)

      integer :: nlunit
!
!     Local variables
!     ---------------
      integer il, i, l
!
      real(kind=kind_phys) ::     SLMASK(len),       &
     &      RLA(len),           RLO(len),          &
     &      OROG(len),          OROG_UF(len),      &
     &      TSFFCS(len),        SNOFCS(len),       &
     &      ZORFCS(len),        ALBFCS(len,4),     &
     &      TG3FCS(len),        CNPFCS(len),       &
     &      SMCFCS(len,LSOIL),  STCFCS(len,LSOIL), &
     &      SLIFCS(len),        AISFCS(len),       &
     &      F10MFCS(len),       VEGFCS(len),       &
     &      VETFCS(len),        SOTFCS(len),       &
     &      ALFFCS(len,2),      CVFCS(len),        &
     &      CVBFCS(len),        CVTFCS(len),       &
     &      SMCFC1(len*LSOIL),  STCFC1(len*LSOIL), &
     &      ALBFC1(len*4),      ALFFC1(len*2),     &
!CluX add swdfcs, sihfcs, sicfcs
     &      SWDFCS(len),        SIHFCS(len),       &
     &      SICFCS(len),        SITFCS(len),       &
!CluX add vmnfcs, vmxfcs, slpfcs, absfcs, slcfc1, slcfcs
     &      VMNFCS(len),        VMXFCS(len),       &
     &      SLPFCS(len),        ABSFCS(len),       &
     &      SLCFC1(len*LSOIL),  SLCFCS(len,LSOIL)


      real(kind=kind_phys) :: sig1t, pifac
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     if (me .eq. 0) print *,' nlats=',nlats,' lonsinpe='
!    *,lonsinpe(0,1)

      sig1t = 0.0
!
      pifac = 180.0 / pi

      DO il = 1,len
!     print *,' calling gcycle for ilat',ilat,' me=',me,' nlats='
!    *,nlats,' lonsinpe=',lonsinpe(:,ilat)
!     if (ilat .eq. nlats) stop
!
          RLA(il)      = XLAT(il) * pifac
          RLO(il)      = XLON(il) * pifac
          OROG(il)     = Sfc_props%ORO(il)
          OROG_UF(il)  = Sfc_props%ORO_UF(il)
          TSFFCS(il)   = Sfc_props%TSFC(il)
          SNOFCS(il)   = Sfc_props%WEASD(il)
          ZORFCS(il)   = Sfc_props%ZORL(il)
          ALBFCS(il,1) = Sfc_props%ALVSF(il)
          ALBFCS(il,2) = Sfc_props%ALVWF(il)
          ALBFCS(il,3) = Sfc_props%ALNSF(il)
          ALBFCS(il,4) = Sfc_props%ALNWF(il)
          TG3FCS(il)   = Sfc_props%TG3(il)
          CNPFCS(il)   = Sfc_props%CANOPY(il)
          SMCFCS(il,:) = Tbd_data%SMC(il,:)
          STCFCS(il,:) = Tbd_data%STC(il,:)
          SLIFCS(il)   = Sfc_props%SLMSK(il)
          F10MFCS(il)  = Sfc_props%F10M(il)
          VEGFCS(il)   = Sfc_props%VFRAC(il)
          VETFCS(il)   = Sfc_props%VTYPE(il)
          SOTFCS(il)   = Sfc_props%STYPE(il)
          ALFFCS(il,1) = Sfc_props%FACSF(il)
          ALFFCS(il,2) = Sfc_props%FACWF(il)
          CVFCS(il)    = Cld_props%CV(il)
          CVBFCS(il)   = Cld_props%CVB(il)
          CVTFCS(il)   = Cld_props%CVT(il)
!CluX add swdfcs, sihfcs, sicfcs
          SWDFCS(il)   = Sfc_props%SNOWD(il)
          SIHFCS(il)   = Sfc_props%HICE(il)
          SICFCS(il)   = Sfc_props%FICE(il)
          SITFCS(il)   = Sfc_props%TISFC(il)
!CluX add slcfcs, vmnfcs, vmxfcs, slpfcs, absfcs
          SLCFCS(il,:) = Tbd_data%SLC(il,:)
          VMNFCS(il)   = Sfc_props%SHDMIN(il)
          VMXFCS(il)   = Sfc_props%SHDMAX(il)
          SLPFCS(il)   = Sfc_props%SLOPE(il)
          ABSFCS(il)   = Sfc_props%SNOALB(il)

!
          IF (SLIFCS(il) .LT. 0.1 .OR. SLIFCS(il) .GT. 1.5) THEN
             SLMASK(il) = 0
          ELSE
             SLMASK(il) = 1
          ENDIF

          IF (SLIFCS(il) .EQ. 2) THEN
            AISFCS(il) = 1.
          ELSE
            AISFCS(il) = 0.
          ENDIF

!     if (me .eq. 0)
!    &   print *,' len=',il,' rla=',rla(il),' rlo=',rlo(il)
      ENDDO       !-----END len LOOP-------------------------------
!
! check
!     print *,' total points = ',il
!
      do l=1,lsoil
        il = (l-1)*len
        do i=1,len
          SMCFC1(il+i) = SMCFCS(i,l)
          STCFC1(il+i) = STCFCS(i,l)
          SLCFC1(il+i) = SLCFCS(i,l)
        enddo
      enddo
      do l=1,4
        il = (l-1)*len
        do i=1,len
          ALBFC1(il+i) = ALBFCS(i,l)
        enddo
      enddo
      do l=1,2
        il = (l-1)*len
        do i=1,len
          ALFFC1(il+i) = ALFFCS(i,l)
        enddo
      enddo
! check
!     call mymaxmin(slifcs,len,len,1,'slifcs')
!     call mymaxmin(slmask,len,len,1,'slmsk')
!
      nlunit = open_namelist_file()
      CALL SFCCYCLE(1001,LEN,LSOIL,SIG1T,fhcyc,                         &
     &              idate(4), idate(2), idate(3), idate(1), fhour,     &
     &              RLA, RLO, SLMASK, OROG, OROG_UF, USE_UFO, nst_anl, &
     &              SIHFCS,   SICFCS, SITFCS,                          &
     &              SWDFCS,   SLCFC1,                                  &
     &              VMNFCS,   VMXFCS, SLPFCS, ABSFCS,                  &
     &              TSFFCS,   SNOFCS, ZORFCS, ALBFC1, TG3FCS,          &
     &              CNPFCS,   SMCFC1, STCFC1, SLIFCS, AISFCS, F10MFCS, &
     &              VEGFCS,   VETFCS, SOTFCS, ALFFC1,                  &
     &              CVFCS,    CVBFCS, CVTFCS, me, nlunit, ialb)
      call close_file(nlunit)
!
      do l=1,lsoil
        il = (l-1)*len
        do i=1,len
          SMCFCS(i,l) = SMCFC1(il+i)
          STCFCS(i,l) = STCFC1(il+i)
          SLCFCS(i,l) = SLCFC1(il+i)
        enddo
      enddo
      do l=1,4
        il = (l-1)*len
        do i=1,len
          ALBFCS(i,l) = ALBFC1(il+i)
        enddo
      enddo
      do l=1,2
        il = (l-1)*len
        do i=1,len
          ALFFCS(i,l) = ALFFC1(il+i)
        enddo
      enddo
!
      DO il=1,len
          Sfc_props%TSFC(il)   = TSFFCS(il)
          Sfc_props%WEASD(il)  = SNOFCS(il)
          Sfc_props%ZORL(il)   = ZORFCS(il)
          Sfc_props%ALVSF(il)  = ALBFCS(il,1)
          Sfc_props%ALVWF(il)  = ALBFCS(il,2)
          Sfc_props%ALNSF(il)  = ALBFCS(il,3)
          Sfc_props%ALNWF(il)  = ALBFCS(il,4)
          Sfc_props%TG3(il)    = TG3FCS(il)
          Sfc_props%CANOPY(il) = CNPFCS(il)
          Tbd_data%SMC(il,:)   = SMCFCS(il,:)
          Tbd_data%STC(il,:)   = STCFCS(il,:)
          Sfc_props%SLMSK(il)  = SLIFCS(il)
          Sfc_props%F10M(il)   = F10MFCS(il)
          Sfc_props%VFRAC(il)  = VEGFCS(il)
          Sfc_props%VTYPE(il)  = VETFCS(il)
          Sfc_props%STYPE(il)  = SOTFCS(il)
          Sfc_props%FACSF(il)  = ALFFCS(il,1)
          Sfc_props%FACWF(il)  = ALFFCS(il,2)
          Cld_props%CV(il)     = CVFCS(il)
          Cld_props%CVB(il)    = CVBFCS(il)
          Cld_props%CVT(il)    = CVTFCS(il)
          Sfc_props%SNOWD(il)  = SWDFCS(il)
          Sfc_props%HICE(il)   = SIHFCS(il)
          Sfc_props%FICE(il)   = SICFCS(il)
          Sfc_props%TISFC(il)  = SITFCS(il)
          Tbd_data%SLC(il,:)   = SLCFCS(il,:)
          Sfc_props%SHDMIN(il) = VMNFCS(il)
          Sfc_props%SHDMAX(il) = VMXFCS(il)
          Sfc_props%SLOPE(il)  = SLPFCS(il)
          Sfc_props%SNOALB(il) = ABSFCS(il)
      ENDDO       !-----END len LOOP-------------------------------
!
!     if (me .eq. 0) print*,'executed gcycle during hour=',fhour
      
      RETURN
      END

