      subroutine mfpbltq(im,ix,km,kmpbl,ntcw,ntrac1,delt,
     &   cnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,thlx,thvx,
     &   gdx,hpbl,kpbl,vpert,buo,use_shear,wush, 
     &   use_tke_ent_det,tkemean,vez0fun,xmf,
     &   tcko,qcko,ucko,vcko,xlamue,a1)
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
      integer              im, ix, km, kmpbl, ntcw, ntrac1
!    &,                    me
      integer              kpbl(im)
      logical              cnvflg(im)
      real(kind=kind_phys) delt
      real(kind=kind_phys) q1(ix,km,ntrac1),
     &                     t1(ix,km),  u1(ix,km), v1(ix,km),
     &                     plyr(im,km),pix(im,km),thlx(im,km),
     &                     thvx(im,km),zl(im,km), zm(im,km),
     &                     wush(im,km),
     &                     gdx(im),    hpbl(im),  vpert(im),
     &                     tkemean(im), vez0fun(im),
     &                     buo(im,km), xmf(im,km),
     &                     tcko(im,km),qcko(im,km,ntrac1),
     &                     ucko(im,km),vcko(im,km),
     &                     xlamue(im,km-1)
      logical use_tke_ent_det,use_shear
!
c  local variables and arrays
!
      integer   i, j, k, n, ndc
      integer   kpblx(im), kpbly(im)
!
      real(kind=kind_phys) dt2,     dz,      ce0,     cm,
     &                     factor,  gocp,
     &                     g,       b1,      f1,
     &                     bb1,     bb2,
     &                     alp,     vprtmax, a1,      pgcon,
     &                     qmin,    qlmin,   xmmx,    rbint,
     &                     tem,     tem1,    tem2,
     &                     ptem,    ptem1,   ptem2,
     &                     tkcrt,   cmxfac
!
      real(kind=kind_phys) elocp,   el2orc,  qs,      es,
     &                     tlu,     gamma,   qlu,
     &                     thup,    thvu,    dq
!
      real(kind=kind_phys) rbdn(im), rbup(im), hpblx(im),
     &                     xlamuem(im,km-1)
      real(kind=kind_phys) delz(im), xlamax(im), ce0t(im)
!
      real(kind=kind_phys) wu2(im,km), thlu(im,km),
     &                     qtx(im,km), qtu(im,km)
!
      real(kind=kind_phys) xlamavg(im),   sigma(im),
     &                     scaldfunc(im), sumx(im)
!
      logical totflg, flg(im)
!
!  physical parameters
      parameter(g=grav)
      parameter(gocp=g/cp)
      parameter(elocp=hvap/cp,el2orc=hvap*hvap/(rv*cp))
      parameter(ce0=0.4,cm=1.0)
      parameter(tkcrt=2.,cmxfac=5.)
      parameter(qmin=1.e-8,qlmin=1.e-12)
      parameter(alp=1.5,vprtmax=3.0,pgcon=0.55)
      parameter(b1=0.5,f1=0.15)
!
!************************************************************************
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
      dt2 = delt
!
      do k = 1, km
        do i=1,im
          if (cnvflg(i)) then
            buo(i,k) = 0.
            wu2(i,k) = 0.
            qtx(i,k) = q1(i,k,1) + q1(i,k,ntcw)
          endif
        enddo
      enddo
!
!  compute thermal excess
!
      do i=1,im
        if(cnvflg(i)) then
          ptem = alp * vpert(i)
          ptem = min(ptem, vprtmax)
          thlu(i,1)= thlx(i,1) + ptem
          qtu(i,1) = qtx(i,1)
          buo(i,1) = g * ptem / thvx(i,1)
        endif
      enddo
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
          k = kpbl(i) / 2
          k = max(k, 1)
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0t(i) / delz(i)
        endif
      enddo
!
      do k = 1, kmpbl
        do i=1,im
          if(cnvflg(i)) then
            if(k < kpbl(i)) then
              ptem = 1./(zm(i,k)+delz(i))
              tem = max((hpbl(i)-zm(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamue(i,k) = ce0t(i) * (ptem+ptem1)
            else
              xlamue(i,k) = xlamax(i)
            endif
!
            xlamuem(i,k) = cm * xlamue(i,k)
          endif
        enddo
      enddo
!
!  compute buoyancy for updraft air parcel
!
      do k = 2, kmpbl
        do i=1,im
          if(cnvflg(i)) then
            dz   = zl(i,k) - zl(i,k-1)
            tem  = 0.5 * xlamue(i,k-1) * dz
            factor = 1. + tem
!
            thlu(i,k) = ((1.-tem)*thlu(i,k-1)+tem*
     &                  (thlx(i,k-1)+thlx(i,k)))/factor
            qtu(i,k) = ((1.-tem)*qtu(i,k-1)+tem*
     &                  (qtx(i,k-1)+qtx(i,k)))/factor
!
            tlu = thlu(i,k) / pix(i,k)
            es = 0.01 * fpvs(tlu)      ! fpvs in pa
            qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
            dq = qtu(i,k) - qs
!
            if (dq > 0.) then
              gamma = el2orc * qs / (tlu**2)
              qlu = dq / (1. + gamma)
              qtu(i,k) = qs + qlu
              tem1 = 1. + fv * qs - qlu
              thup = thlu(i,k) + pix(i,k) * elocp * qlu
              thvu = thup * tem1
            else
              tem1 = 1. + fv * qtu(i,k)
              thvu = thlu(i,k) * tem1
            endif
            buo(i,k) = g * (thvu / thvx(i,k) - 1.)
!
          endif
        enddo
      enddo
!
!  compute updraft velocity square(wu2)
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
          dz   = zm(i,1)
          tem  = 0.5*bb1*xlamue(i,1)*dz
          tem1 = bb2 * buo(i,1) * dz
          ptem1 = 1. + tem
          wu2(i,1) = tem1 / ptem1
        endif
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if(cnvflg(i)) then
            dz    = zm(i,k) - zm(i,k-1)
            tem  = 0.25*bb1*(xlamue(i,k-1)+xlamue(i,k))*dz
            ! kgao 12/15/2023 - consider shear effect on wu diagnosis
            if (use_shear) then
              tem1 = max(wu2(i,k-1), 0.)
              tem1 = bb2 * buo(i,k) - wush(i,k) * sqrt(tem1)
              tem2 = tem1 * dz
              ptem = (1. - tem) * wu2(i,k-1)
              ptem1 = 1. + tem
              wu2(i,k) = (ptem + tem1) / ptem1
            else
              tem1 = bb2 * buo(i,k) * dz
              ptem = (1. - tem) * wu2(i,k-1)
              ptem1 = 1. + tem
              wu2(i,k) = (ptem + tem1) / ptem1
             endif
          endif
        enddo
      enddo
!
!  update pbl height as the height where updraft velocity vanishes
!
      do i=1,im
         flg(i)  = .true.
         kpblx(i) = 1
         kpbly(i) = kpbl(i)
         if(cnvflg(i)) then
           flg(i)  = .false.
           rbup(i) = wu2(i,1)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          rbup(i) = wu2(i,k)
          kpblx(i)= k
          flg(i)  = rbup(i).le.0.
        endif
      enddo
      enddo
      do i = 1,im
        if(cnvflg(i)) then
           k = kpblx(i)
           if(rbdn(i) <= 0.) then
              rbint = 0.
           elseif(rbup(i) >= 0.) then
              rbint = 1.
           else
              rbint = rbdn(i)/(rbdn(i)-rbup(i))
           endif
           hpblx(i) = zm(i,k-1) + rbint*(zm(i,k)-zm(i,k-1))
        endif
      enddo
!
      do i = 1,im
        if(cnvflg(i)) then
          if(kpblx(i) < kpbl(i)) then
            kpbl(i) = kpblx(i)
            hpbl(i) = hpblx(i)
          endif
          if(kpbl(i) <= 1) cnvflg(i)=.false.
        endif
      enddo
! 
!  update entrainment rate
!
      do i=1,im
        if(cnvflg(i)) then
          k = kpbl(i) / 2
          k = max(k, 1)
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0t(i) / delz(i)
        endif
      enddo
!
      do k = 1, kmpbl
        do i=1,im
          if(cnvflg(i) .and. kpblx(i) < kpbly(i)) then
!         if(cnvflg(i)) then
            if(k < kpbl(i)) then
              ptem = 1./(zm(i,k)+delz(i))
              tem = max((hpbl(i)-zm(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamue(i,k) = ce0t(i) * (ptem+ptem1)
            else 
              xlamue(i,k) = xlamax(i)
            endif
!
            xlamuem(i,k) = cm * xlamue(i,k)
          endif
        enddo
      enddo
!
!  compute entrainment rate averaged over the whole pbl
!
      do i = 1, im
        xlamavg(i) = 0.
        sumx(i) = 0.
      enddo
      do k = 1, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k < kpbl(i)) then
            dz = zl(i,k+1) - zl(i,k)
            xlamavg(i) = xlamavg(i) + xlamue(i,k) * dz
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
!  updraft mass flux as a function of updraft velocity profile
!
      do k = 1, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k < kpbl(i)) then
             xmf(i,k) = a1 * sqrt(wu2(i,k))
          endif
        enddo
      enddo
!
!--- compute updraft fraction as a function of mean entrainment rate
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
          if (sigma(i) > a1) then
            scaldfunc(i) = (1.-sigma(i)) * (1.-sigma(i))
            scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          else
            scaldfunc(i) = 1.0
          endif
        endif
      enddo
!
!  final scale-aware updraft mass flux
!
      do k = 1, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k < kpbl(i)) then
             xmf(i,k) = scaldfunc(i) * xmf(i,k)
             dz   = zl(i,k+1) - zl(i,k)
             xmmx = dz / dt2
             xmf(i,k) = min(xmf(i,k),xmmx)
          endif
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute updraft property using updated entranment rate
!
      do i=1,im
        if(cnvflg(i)) then
          thlu(i,1)= thlx(i,1)
        endif
      enddo
!
!     do i=1,im
!       if(cnvflg(i)) then
!         ptem1 = max(qcko(i,1,ntcw), 0.)
!         tlu = thlu(i,1) / pix(i,1)
!         tcko(i,1) = tlu +  elocp * ptem1
!       endif
!     enddo
!
      do k = 2, kmpbl
        do i=1,im
          if(cnvflg(i) .and. k <= kpbl(i)) then
            dz   = zl(i,k) - zl(i,k-1)
            tem  = 0.5 * xlamue(i,k-1) * dz
            factor = 1. + tem
!
            thlu(i,k) = ((1.-tem)*thlu(i,k-1)+tem*
     &                  (thlx(i,k-1)+thlx(i,k)))/factor
            qtu(i,k) = ((1.-tem)*qtu(i,k-1)+tem*
     &                  (qtx(i,k-1)+qtx(i,k)))/factor
!
            tlu = thlu(i,k) / pix(i,k)
            es = 0.01 * fpvs(tlu)      ! fpvs in pa
            qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
            dq = qtu(i,k) - qs
!
            if (dq > 0.) then
              gamma = el2orc * qs / (tlu**2)
              qlu = dq / (1. + gamma)
              qtu(i,k) = qs + qlu
              qcko(i,k,1) = qs
              qcko(i,k,ntcw) = qlu
              tcko(i,k) = tlu + elocp * qlu
            else
              qcko(i,k,1) = qtu(i,k)
              qcko(i,k,ntcw) = 0.
              tcko(i,k) = tlu
            endif
!
          endif
        enddo
      enddo
!
      do k = 2, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamuem(i,k-1) * dz
             factor = 1. + tem
             ptem = tem + pgcon
             ptem1= tem - pgcon
             ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*u1(i,k)
     &                    +ptem1*u1(i,k-1))/factor
             vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*v1(i,k)
     &                    +ptem1*v1(i,k-1))/factor
          endif
        enddo
      enddo
!
      if(ntcw > 2) then
!
      do n = 2, ntcw-1
      do k = 2, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem
! 
             qcko(i,k,n) = ((1.-tem)*qcko(i,k-1,n)+tem*
     &                    (q1(i,k,n)+q1(i,k-1,n)))/factor
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
      do k = 2, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem
! 
             qcko(i,k,n) = ((1.-tem)*qcko(i,k-1,n)+tem*
     &                    (q1(i,k,n)+q1(i,k-1,n)))/factor
          endif
        enddo
      enddo
      enddo
!
      endif
!
      return
      end
