      subroutine mfscuq(im,ix,km,kmscu,ntcw,ntrac1,delt,
     &   cnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,
     &   thlx,thvx,thlvx,gdx,thetae,
     &   krad,mrad,radmin,buo,
     &   use_shear,wush,
     &   use_tke_ent_det,tkemean,vez0fun,xmfd,
     &   tcdo,qcdo,ucdo,vcdo,xlamde,a1)
!
      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp
     &,             rv => con_rv, hvap => con_hvap
     &,             fv => con_fvirt
     &,             eps => con_eps, epsm1 => con_epsm1
!
      implicit none
!
      integer            im, ix,  km, kmscu, ntcw, ntrac1
!    &,                  me
      integer   krad(im), mrad(im)
!
      logical cnvflg(im)
      real(kind=kind_phys) delt
      real(kind=kind_phys) q1(ix,km,ntrac1),t1(ix,km),
     &                     u1(ix,km),      v1(ix,km),
     &                     plyr(im,km),    pix(im,km),
     &                     thlx(im,km),    wush(im,km),
     &                     thvx(im,km),    thlvx(im,km),
     &                     gdx(im),
     &                     zl(im,km),      zm(im,km),
     &                     thetae(im,km),  radmin(im),
     &                     tkemean(im), vez0fun(im),
     &                     buo(im,km), xmfd(im,km),
     &                     tcdo(im,km), qcdo(im,km,ntrac1),
     &                     ucdo(im,km), vcdo(im,km),
     &                     xlamde(im,km-1)
      logical use_tke_ent_det, use_shear
!
!  local variables and arrays
!
!
      integer   i,j,indx, k, n, kk, ndc
      integer   krad1(im)
!
      real(kind=kind_phys) dt2,     dz,      ce0,     cm,
     &                     gocp,    factor,  g,       tau,
     &                     b1,      f1,      bb1,     bb2,
     &                     a1,      a2,
     &                     cteit,   pgcon,
     &                     qmin,    qlmin,
     &                     xmmx,    tem,     tem1,    tem2,
     &                     ptem,    ptem1,   ptem2,
     &                     tkcrt,   cmxfac
!
      real(kind=kind_phys) elocp,   el2orc,  qs,      es,
     &                     tld,     gamma,   qld,     thdn,
     &                     thvd,    dq
!
      real(kind=kind_phys) wd2(im,km), thld(im,km),
     &                     qtx(im,km), qtd(im,km),
     &                     thlvd(im),  hrad(im),
     &                     xlamdem(im,km-1), ra1(im)
      real(kind=kind_phys) delz(im), xlamax(im), ce0t(im)
!
      real(kind=kind_phys) xlamavg(im),   sigma(im),
     &                     scaldfunc(im), sumx(im)
!
      logical totflg, flg(im)
!
      real(kind=kind_phys) actei, cldtime
!
c  physical parameters
      parameter(g=grav)
      parameter(gocp=g/cp)
      parameter(elocp=hvap/cp,el2orc=hvap*hvap/(rv*cp))
      parameter(ce0=0.4,cm=1.0,pgcon=0.55)
      parameter(tkcrt=2.,cmxfac=5.)
      parameter(qmin=1.e-8,qlmin=1.e-12)
      parameter(b1=0.45,f1=0.15)
      parameter(a2=0.5)
      parameter(cldtime=500.)
      parameter(actei = 0.7)
!     parameter(actei = 0.23)
!
!************************************************************************
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
      dt2 = delt
!
      do k = 1, km
        do i=1,im
          if(cnvflg(i)) then
            buo(i,k) = 0.
            wd2(i,k) = 0.
            qtx(i,k) = q1(i,k,1) + q1(i,k,ntcw)
          endif
        enddo
      enddo
!
      do i = 1, im
        if(cnvflg(i)) then
           hrad(i) = zm(i,krad(i))
           krad1(i) = krad(i)-1
        endif
      enddo
!
      do i = 1, im
        if(cnvflg(i)) then
          k    = krad(i)
          tem  = zm(i,k+1)-zm(i,k)
          tem1 = cldtime*radmin(i)/tem
          tem1 = max(tem1, -3.0)
          thld(i,k)= thlx(i,k) + tem1
          qtd(i,k) = qtx(i,k)
          thlvd(i) = thlvx(i,k) + tem1
          buo(i,k) = - g * tem1 / thvx(i,k)
        endif
      enddo
!
!  specify downdraft fraction
!
      do i=1,im
        if(cnvflg(i)) then
          ra1(i) = a1
        endif
      enddo
!
!  if the condition for cloud-top instability is met,
!    increase downdraft fraction
!
      do i = 1, im
        if(cnvflg(i)) then
           k = krad(i)
           tem = thetae(i,k) - thetae(i,k+1)
           tem1 = qtx(i,k) - qtx(i,k+1)
           if (tem > 0. .and. tem1 > 0.) then
             cteit= cp*tem/(hvap*tem1)
             if(cteit > actei) then
               ra1(i) = a2
             endif
           endif
        endif
      enddo
!
!   first-quess level of downdraft extension (mrad)
! 
      do i = 1, im
        flg(i) = cnvflg(i)
        mrad(i) = krad(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k < krad(i)) then
          if(thlvd(i) <= thlvx(i,k)) then
             mrad(i) = k
          else
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if (cnvflg(i)) then
          kk = krad(i)-mrad(i)
          if(kk < 1) cnvflg(i)=.false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!
!  compute entrainment rate
!
      do i = 1, im
        if ( cnvflg(i) ) then
          if ( use_tke_ent_det ) then
            ! kgao 12/08/2023: compute entrainment/detrainment rate based on pbl-mean tke
            ce0t(i) = ce0 * vez0fun(i)
            if ( tkemean(i) > tkcrt ) then
              tem = sqrt(tkemean(i)/tkcrt)
              tem1 = min(tem, cmxfac)
              tem2 = tem1 * ce0
              ce0t(i) = max(ce0t(i), tem2)
            endif
          else
            ce0t(i) = ce0
          endif
          k = mrad(i) + (krad(i)-mrad(i)) / 2
          k = max(k, mrad(i))
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0t(i) / delz(i)
        endif
      enddo
!
      do k = 1, kmscu
        do i=1,im
          if(cnvflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              if(mrad(i) == 1) then
                ptem = 1./(zm(i,k)+delz(i))
              else
                ptem = 1./(zm(i,k)-zm(i,mrad(i)-1)+delz(i))
              endif
              tem = max((hrad(i)-zm(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamde(i,k) = ce0t(i) * (ptem+ptem1)
            else
              xlamde(i,k) = xlamax(i)
            endif
!
            xlamdem(i,k) = cm * xlamde(i,k)
          endif
        enddo
      enddo
!
!  compute buoyancy for downdraft air parcel
!
      do k = kmscu,1,-1
        do i=1,im
          if(cnvflg(i) .and. k < krad(i)) then
            dz = zl(i,k+1) - zl(i,k)
            tem  = 0.5 * xlamde(i,k) * dz
            factor = 1. + tem
! 
            thld(i,k) = ((1.-tem)*thld(i,k+1)+tem*
     &                     (thlx(i,k)+thlx(i,k+1)))/factor
            qtd(i,k) = ((1.-tem)*qtd(i,k+1)+tem*
     &                     (qtx(i,k)+qtx(i,k+1)))/factor
!
            tld = thld(i,k) / pix(i,k)
            es = 0.01 * fpvs(tld)      ! fpvs in pa
            qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
            dq = qtd(i,k) - qs
!
            if (dq > 0.) then
              gamma = el2orc * qs / (tld**2)
              qld = dq / (1. + gamma)
              qtd(i,k) = qs + qld
              tem1 = 1. + fv * qs - qld
              thdn = thld(i,k) + pix(i,k) * elocp * qld
              thvd = thdn * tem1
            else
              tem1 = 1. + fv * qtd(i,k)
              thvd = thld(i,k) * tem1
            endif
            buo(i,k) = g * (1. - thvd / thvx(i,k))
!
          endif
        enddo
      enddo
!
!  compute downdraft velocity square(wd2)
!
!     tem = 1.-2.*f1
!     bb1 = 2. * b1 / tem
!     bb2 = 2. / tem
!  from Soares et al. (2004,QJRMS)
!     bb1 = 2.
!     bb2 = 4.
!
!  from Bretherton et al. (2004, MWR)
!     bb1 = 4.
!     bb2 = 2.
!
!  from our tuning
      bb1 = 2.0
      bb2 = 4.0
!
      do i = 1, im
        if(cnvflg(i)) then
          k = krad1(i)
          dz = zm(i,k+1) - zm(i,k)
!         tem = 0.25*bb1*(xlamde(i,k)+xlamde(i,k+1))*dz
          tem = 0.5*bb1*xlamde(i,k)*dz
          tem1 = bb2 * buo(i,k+1) * dz
          ptem1 = 1. + tem
          wd2(i,k) = tem1 / ptem1
        endif
      enddo
      do k = kmscu,1,-1
        do i = 1, im
          if(cnvflg(i) .and. k < krad1(i)) then
            dz    = zm(i,k+1) - zm(i,k)
            tem  = 0.25*bb1*(xlamde(i,k)+xlamde(i,k+1))*dz
            ! kgao 12/15/2023 - consider shear effect on wd diagnosis
            if (use_shear) then
              tem1 = max(wd2(i,k+1), 0.)
              tem1 = bb2*buo(i,k+1) - wush(i,k+1)*sqrt(tem1)
              tem2 = tem1 * dz
              ptem = (1. - tem) * wd2(i,k+1)
              ptem1 = 1. + tem
              wd2(i,k) = (ptem + tem2) / ptem1
            else
              tem1 = bb2 * buo(i,k+1) * dz
              ptem = (1. - tem) * wd2(i,k+1)
              ptem1 = 1. + tem
              wd2(i,k) = (ptem + tem1) / ptem1
            endif
          endif
        enddo
      enddo
c
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) mrad(i) = krad(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k < krad(i)) then
          if(wd2(i,k) > 0.) then
            mrad(i) = k
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
!
      do i=1,im
        if (cnvflg(i)) then
          kk = krad(i)-mrad(i)
          if(kk < 1) cnvflg(i)=.false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  update entrainment rate
!
      do i=1,im
        if(cnvflg(i)) then
          k = mrad(i) + (krad(i)-mrad(i)) / 2
          k = max(k, mrad(i))
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0t(i) / delz(i)
        endif
      enddo
!
      do k = 1, kmscu
        do i=1,im
          if(cnvflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              if(mrad(i) == 1) then
                ptem = 1./(zm(i,k)+delz(i))
              else
                ptem = 1./(zm(i,k)-zm(i,mrad(i)-1)+delz(i))
              endif
              tem = max((hrad(i)-zm(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamde(i,k) = ce0t(i) * (ptem+ptem1)
            else
              xlamde(i,k) = xlamax(i)
            endif
!
            xlamdem(i,k) = cm * xlamde(i,k)
          endif
        enddo
      enddo
!
!  compute entrainment rate averaged over the whole downdraft layers
!
      do i = 1, im
        xlamavg(i) = 0.
        sumx(i) = 0.
      enddo
      do k = kmscu, 1, -1
        do i = 1, im
          if(cnvflg(i) .and.
     &       (k >= mrad(i) .and. k < krad(i))) then
            dz = zl(i,k+1) - zl(i,k)
            xlamavg(i) = xlamavg(i) + xlamde(i,k) * dz
            sumx(i) = sumx(i) + dz
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
           xlamavg(i) = xlamavg(i) / sumx(i)
        endif
      enddo
!
!  compute downdraft mass flux
!
      do k = kmscu, 1, -1
        do i = 1, im
          if(cnvflg(i) .and.
     &      (k >= mrad(i) .and. k < krad(i))) then
              xmfd(i,k) = ra1(i) * sqrt(wd2(i,k))
          endif
        enddo
      enddo
!
!--- compute downdraft fraction as a function of mean entrainment rate
!        (Grell & Freitas, 2014)
!
      do i = 1, im
        if(cnvflg(i)) then
          tem = 0.2 / xlamavg(i)
          tem1 = 3.14 * tem * tem
          sigma(i) = tem1 / (gdx(i) * gdx(i))
          sigma(i) = max(sigma(i), 0.001)
          sigma(i) = min(sigma(i), 0.999)
        endif
      enddo
!
!--- compute scale-aware function based on Arakawa & Wu (2013)
!
      do i = 1, im
        if(cnvflg(i)) then
          if (sigma(i) > ra1(i)) then
            scaldfunc(i) = (1.-sigma(i)) * (1.-sigma(i))
            scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          else
            scaldfunc(i) = 1.0
          endif
        endif
      enddo
!
!  final scale-aware downdraft mass flux
!
      do k = kmscu, 1, -1
        do i = 1, im
          if(cnvflg(i) .and.
     &       (k >= mrad(i) .and. k < krad(i))) then
             xmfd(i,k) = scaldfunc(i) * xmfd(i,k)
             dz   = zl(i,k+1) - zl(i,k)
             xmmx = dz / dt2
             xmfd(i,k) = min(xmfd(i,k),xmmx)
          endif
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute downdraft property using updated entranment rate
!
      do i = 1, im
        if(cnvflg(i)) then
          k = krad(i)
          thld(i,k)= thlx(i,k)
        endif
      enddo
!
!     do i = 1, im
!       if(cnvflg(i)) then
!         k = krad(i)
!         ptem1 = max(qcdo(i,k,ntcw), 0.)
!         tld = thld(i,k) / pix(i,k)
!         tcdo(i,k) = tld +  elocp * ptem1
!         qcdo(i,k,1) = qcdo(i,k,1)+0.2*qcdo(i,k,1)
!         qcdo(i,k,ntcw) = qcdo(i,k,ntcw)+0.2*qcdo(i,k,ntcw)
!       endif
!     enddo
!
      do k = kmscu,1,-1
        do i=1,im
          if(cnvflg(i) .and. 
     &       (k >= mrad(i) .and. k < krad(i))) then
            dz = zl(i,k+1) - zl(i,k)
            tem  = 0.5 * xlamde(i,k) * dz
            factor = 1. + tem
!
            thld(i,k) = ((1.-tem)*thld(i,k+1)+tem*
     &                     (thlx(i,k)+thlx(i,k+1)))/factor
            qtd(i,k) = ((1.-tem)*qtd(i,k+1)+tem*
     &                     (qtx(i,k)+qtx(i,k+1)))/factor
!
            tld = thld(i,k) / pix(i,k)
            es = 0.01 * fpvs(tld)      ! fpvs in pa
            qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
            dq = qtd(i,k) - qs
!
            if (dq > 0.) then
              gamma = el2orc * qs / (tld**2)
              qld = dq / (1. + gamma)
              qtd(i,k) = qs + qld
              qcdo(i,k,1) = qs
              qcdo(i,k,ntcw) = qld
              tcdo(i,k) = tld + elocp * qld
            else
              qcdo(i,k,1) = qtd(i,k)
              qcdo(i,k,ntcw) = 0.
              tcdo(i,k) = tld
            endif
!
          endif
        enddo
      enddo
!
      do k = kmscu, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5 * xlamdem(i,k) * dz
              factor = 1. + tem
              ptem = tem - pgcon
              ptem1= tem + pgcon
!
              ucdo(i,k) = ((1.-tem)*ucdo(i,k+1)+ptem*u1(i,k+1)
     &                     +ptem1*u1(i,k))/factor
              vcdo(i,k) = ((1.-tem)*vcdo(i,k+1)+ptem*v1(i,k+1)
     &                     +ptem1*v1(i,k))/factor
            endif
          endif
        enddo
      enddo
!
      if(ntcw > 2) then
!
      do n = 2, ntcw-1
      do k = kmscu, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5 * xlamde(i,k) * dz
              factor = 1. + tem
! 
              qcdo(i,k,n) = ((1.-tem)*qcdo(i,k+1,n)+tem*
     &                       (q1(i,k,n)+q1(i,k+1,n)))/factor
            endif
          endif
        enddo
      enddo
      enddo
!
      endif
!
      ndc = ntrac1 - ntcw
!
      if(ndc > 0) then
!
      do n = ntcw+1, ntrac1
      do k = kmscu, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5 * xlamde(i,k) * dz
              factor = 1. + tem
! 
              qcdo(i,k,n) = ((1.-tem)*qcdo(i,k+1,n)+tem*
     &                       (q1(i,k,n)+q1(i,k+1,n)))/factor
            endif
          endif
        enddo
      enddo
      enddo
!
      endif
!
      return
      end
