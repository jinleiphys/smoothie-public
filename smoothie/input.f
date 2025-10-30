cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This file gives the informations of input and output
c           a + A -> b + x + A -> x + b + A -> b+anything
c           =====        =====        =====
c             a            x            b
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      module input
      implicit none
      logical,dimension(1:9999) ::  written
       real*8,private :: ecmbmin, ecmbmax
       integer,private :: bin
       real*8,private :: wbin
       real*8,private :: ecmbh
       logical :: prior
       logical :: printf
       logical :: icf
       logical :: HM,cdccebu,UT
       integer :: dwba,cdcc
       logical :: surfacebw, surfacesw
       logical :: ubx
       logical :: VF ! logical variable for the Vincent-Fortune method 
       real*8 :: rs
       logical :: ZR ! logical variable for the Zero-Range approximation
      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initialize()
c     parameter initialization
c       namelist /global/ hcm,rmatch,lmax,elab,thmin,thmax
c                        thinc,cutl,lmin,rmax,nrx,nrmatch,nx,icf
c
c       namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt,
c                         nameb,massb,zb,jb,namex,massx,zx,jx,
c						  qval,neb
c
c       namelist /outgoing/ kn,ecmb,bin,wbin
c
c       namelist /potential/ kp1,kp2,ptype,a1,a2,rc,uv,av,
c                           rv,uw,aw,rw,vsov,rsov,asov,uw1,aw1,rw1
c                           vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use precision
      use cleb_coefficient
      use gauss
      implicit none


       integer :: counter

       integer :: ptype
       character(len=1) :: kp1
       integer :: kp2
       real*8,dimension(0:10) :: uv ! potential depth of real part of W-S
       real*8 :: av,rv ! parameters of real part of W-S
       real*8 :: uw,aw,rw ! parameters of imaginary part of W-S
       real*8 :: uw1,aw1,rw1
       real*8 :: vsov,rsov,asov ! real part of spin-orbit potential for projectile
       real*8 :: vsow,rsow,asow ! imaginary part of spin-orbit potential for projectile
       real*8 ::  vd,avd,rvd  ! real part of surface potential
       real*8 ::  wd,awd,rwd  ! imaginary part of surface potential
       real*8 :: a1,a2    ! mass for radius conversion
       real*8 :: rc ! for coulomb potential
       real*8  :: ecut, dump
       real*8 ::  renv,renw,renvd,renwd,renvls,renwls ! AMoro
       logical :: simpson_int
      !  integer :: ios


c      character(len=5) :: decide1,decide2
      namelist /global/ hcm,lmax,elab,thmin,thmax,thinc,
     &                  lmin,nx,rmax,nr,printf,jtmax,jtmin,
     &                  dwba,lxmax,lxmin

      namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt,
     &                   nameb,massb,zb,jb,namex,massx,zx,jx,sbx,
     &                   lbx,nodes,be


      namelist /outgoing/ ecmbmin,ecmbmax,ecmbh

      namelist /potential/ kp1,ptype,a1,a2,rc,uv,av,
     &                      rv,uw,aw,rw,vsov,rsov,asov,
     &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

      ! namelist /vincent/ Pi_ext_ny,Pi_ext_ymax,VF

       written=.false.
       written(1)=.true.;
C       open (unit=5,file='test.in')
c-----------------------------------------------------------------------
c /global/
       hcm=0.05_dpreal;r1=15_dpreal;r2=60_dpreal;rmax=7000_dpreal
       lmax=30;elab=0.0_dpreal
       thmin=0.0_dpreal;thmax=180.0_dpreal;thinc=1.0_dpreal
       cutl=200;lmin=0;nx=30;prior=.true.
       nr=60;printf=.false.;icf=.false.
       jtmin=0.0_dpreal
       jtmax=-30.0_dpreal
       HM=.false.
       dwba=1
       cdcc=-99
       cdccebu=.false.
       gswf=""
       surfacebw=.false.
       surfacesw=.false.
       rs=-9999.0_dpreal
       lxmin=0
       lxmax=-99
       ubx=.false.
       adde=.false.
       UT=.false.
       simpson_int=.false.
       ZR=.false.
       read(5,nml=global)

       write(1,nml=global)

       if(ZR) dwba=1

       if(cutl>lmax) cutl=lmax
       if(lxmax>0) then
        cutl=lxmax
       else
        lxmax=max(lxmax,cutl)
       end if

      !  write(*,*) "lxmax=",lxmax

       if (jtmax<0) jtmax=2.0_dpreal*lmax

       rmatch=rmax
       irmatch=nint(rmatch/hcm)

C       if(cutl>lmax) cutl=lmax

       
      if(simpson_int) then 
          allocate(rr(1:irmatch),rrw(1:irmatch))
          call simpson(irmatch,0.0_dpreal,rmax,rr,rrw)
          nr=irmatch
      else 
          allocate(rr(1:nr),rrw(1:nr))
          call gauleg(nr,0.0_dpreal,rmax,rr,rrw)
      end if  

c-----------------------------------------------------------------------
c/system/
       namep='null';massp=0.0d0;zp=0.0d0;jp=0.0d0;namet='null'
       masst=0.0d0;zt=0.0d0;jt=0.0d0;nameb='null';massb=0.0d0
       zb=0.d0;jb=0.0d0;namex='null';massx=0.0d0;zx=0.0d0
       jx=0.0d0;qval=0.0d0;neb=1;lbx=0;nodes=1;be=0.0001
       ptyp=1;ptyt=1;ptyb=1;ptyx=1
       sbx=0.0_dpreal
       read(5,nml=system)

       write(1,nml=system)

       if(abs(jb+jx+jt)>0.0001) then
          write(*,*) "spin dependent calculation, SET DWBA=3"
          DWBA=3
       end if

       if (jb<0.1 .or. jx<0.1) sbx=max(jb,jx)

       be=-max(abs(be), abs(qval))
       qval=-max(abs(be), abs(qval))

       lbxmin=lbx
       lbxmax=lbx
c-----------------------------------------------------------------------
c/outgoing/
      ecmbmin=0.0d0;ecmbmax=0.0d0;ecmbh=0.0d0;bin=0;wbin=0.0_dpreal


      read(5,nml=outgoing)

      write(1,nml=outgoing)
      neb=nint((ecmbmax-ecmbmin)/ecmbh) + 1
      allocate(nkn(1:neb),necmb(1:neb),nbin(1:neb),nwbin(1:neb))



      do counter=1, neb
         nkn(counter)=counter
             necmb(counter)=ecmbmin+(counter-1)*ecmbh
             nbin(counter)=bin
             nwbin(counter)=wbin
      end do

      if(abs(necmb(neb)-ecmbmax)>0.00001) then
         write(*,10)
         stop
      end if
10    format('!!error!! please check the number of outgoing channels')

c-----------------------------------------------------------------------
c /potential/ parameter
c      if (ecmbh<0.000001) then
      
      counter=1
      nkp2=0
      nuv=0.0_dpreal
      nuw=0.0_dpreal
      nvsov=0.0_dpreal
      nvsow=0.0_dpreal
      nvd=0.0_dpreal
      nwd=0.0_dpreal
      
      do
       kp1=' ';kp2=1;ptype=0;a1=0.0d0;a2=0.0d0;uv=-99.0d0;av=0.1d0
       rv=0.0d0
       uw=0.0d0;aw=0.1d0;rw=0.0d0;rc=1.3d0;vsov=0.0d0;rsov=0.0d0
       asov=0.1d0; vsow=0.0d0;rsow=0.0d0;asow=0.1d0
       vd=0.0d0;avd=0.1d0;rvd=0.0d0;wd=0.0d0;awd=0.1d0;rwd=0.0d0
       uw1=0.0d0;aw1=0.1d0;rw1=0.0d0;ecut=-99.0_dpreal; dump=-99.0_dpreal
! AMoro
       renv=1.0; renw=1.0; renvd=1.0; renwd=1.0; renvls=1.0; renwls=1.0;
       efix=0.0
       read(5,nml=potential)
       if(ptype>40 .and. uw<0.00001 .and. uv(0)<0.0001) then
         uw=1.0_dpreal
         uv(0)=1.0_dpreal
       end if


       write(1,nml=potential)




       if(kp1==' ') exit
       nkp1(counter)=kp1;nptype(counter)=ptype;na1(counter)=a1
       nkp2(counter)=kp2
       na2(counter)=a2;nuv(:,counter)=uv;nav(counter)=av;nrv(counter)=rv
       nuw(counter)=uw;naw(counter)=aw;nrw(counter)=rw;nrc(counter)=rc
       nvsov(counter)=vsov;nrsov(counter)=rsov;nasov(counter)=asov
       nvsow(counter)=vsow;nrsow(counter)=rsow;nasow(counter)=asow
       nvd(counter)=vd;navd(counter)=avd;nrvd(counter)=rvd
       nwd(counter)=wd;nawd(counter)=awd;nrwd(counter)=rwd
       nuw1(counter)=uw1;naw1(counter)=aw1;nrw1(counter)=rw1
       necut(counter)=ecut; ndump(counter)=dump
       nrenv(counter)=renv; nrenw(counter)=renw;
       nrenvd(counter)=renvd; nrenwd(counter)=renwd
       nrenvls(counter)=renvls; nrenwls(counter)=renwls;

       nefix(counter)=efix
       counter=counter+1

      end do




c-----------------------------------------------------------------------
      call factorialgen(6*lmax+20)
      nthmax=nint((thmax-thmin)/thinc)

c-----------------------------------------------------------------------
c /vincent/ parameter
      Pi_ext_ny=nr
      Pi_ext_ymax=rmax
      VF=.false.
      ! read(5,nml=vincent,iostat=ios)
      ! write(1,nml=vincent)
      allocate(yy_mesh(1:Pi_ext_ny), yyw(1:Pi_ext_ny))
      call gauleg(Pi_ext_ny, 0.0_dpreal, Pi_ext_ymax, yy_mesh, yyw)
      if(prior) VF=.false.
      

      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check()
c     check the input file and give the local copy of input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use constants
      implicit none
      !  integer :: bin
      !  real*8 :: wbin
       integer :: counter

      !  integer :: ptype
      !  character(len=1) :: kp1
      !  integer :: kp2
      !  real*8,dimension(0:10) :: uv ! potential depth of real part of W-S
      !  real*8 :: av,rv ! parameters of real part of W-S
      !  real*8 :: uw,aw,rw ! parameters of imaginary part of W-S
      !  real*8 :: uw1,aw1,rw1 ! parameters of imaginary part of W-S
      !  real*8 :: vsov,rsov,asov ! real part of spin-orbit potential for projectile
      !  real*8 :: vsow,rsow,asow ! imaginary part of spin-orbit potential for projectile
      !  real*8 ::  vd,avd,rvd  ! real part of surface potential
      !  real*8 ::  wd,awd,rwd  ! imaginary part of surface potential
      !  real*8 :: a1,a2    ! mass for radius conversion
      !  real*8 :: rc ! for coulomb potential



!        namelist /global/ hcm,lmax,elab,thmin,thmax,thinc,icf,HM,
!      &                  cutl,lmin,nx,rmax,nr,prior,printf,jtmax,jtmin,printlx,
!      &                  dwba,cdcc,cdccebu,surfacebw,surfacesw,rs,lxmax,lxmin,ubx

!        namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt,
!      &                   nameb,massb,zb,jb,namex,massx,zx,jx,sbx,
!      &                   qval,neb,lbx,nodes,be, ptyp, ptyt, ptyb, ptyx,
!      &                  gswf ! added by AMoro

!        namelist /outgoing/ ecmbmin,ecmbmax,ecmbh,bin,wbin

!        namelist /potential/ kp1,kp2,ptype,a1,a2,rc,uv,av,
!      &                      rv,uw,aw,rw,vsov,rsov,asov,uw1,aw1,rw1,
!      &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd! AMoro
!        write(1,nml=global)

!        write(1,nml=system)

!        write(1,nml=outgoing)


! c /potential/ parameter-------------------------------------------------
!        do  counter=1,99
!         kp1=nkp1(counter);ptype=nptype(counter);a1=na1(counter)
!         kp2=nkp2(counter)
!         a2=na2(counter);uv=nuv(:,counter);av=nav(counter);rv=nrv(counter)
!         uw=nuw(counter);aw=naw(counter);rw=nrw(counter);rc=nrc(counter)
!         vsov=nvsov(counter);rsov=nrsov(counter);asov=nasov(counter)
!         vsow=nvsow(counter);rsow=nrsow(counter);asow=nasow(counter)
!         vd=nvd(counter);avd=navd(counter);rvd=nrvd(counter)
!         wd=nwd(counter);awd=nawd(counter);rwd=nrwd(counter)
!         uw1=nuw1(counter);aw1=naw1(counter);rw1=nrw1(counter)
!         if (kp2==0) exit
!         write(1,nml=potential)
!        end do




!        if(neb==0) then
!        write(*,24)
!        stop
!        end if
! 24     format('!!error!!, the number of neb can not be 0')




!        if (thmax==0. .OR. thinc==0.) then
!        write(*,21)
!        stop
!        end if
! 21     format('!!error!! Please check the value of thinc and thmax!')


c-----------------------------------------------------------------------
c *** Print physical constants

c *** Print reaction system
      write(*,40) namep,namet,nameb,namex,namet,nameb,namex,namet,nameb
40    format('Assuming the reaction has the form:',/,T5,A5,'+  ',A5,
     &         '-->  ',A5,'+  ',A5,"+  ",A5,'-->  ',A5,'+  (',A5,'+  ',A5,
     &         ')-->  ',A5,'+ anything' )
c ***print parameters

      write(*,70) nint(rmatch/hcm),hcm,rmatch
70    format('Centre-of-mass Range is ',I5,'*',f5.3,' fm.',
     &       ' Maximum at ',f8.3,' fm.' )

      write(*,80)0,lmax
80    format('Range of angular momentum l is ',I1,' <= l <=',I3)
      write(*,*)

c***print reaction systems
      write(*,90)
      write(*,91)
      write(*,90)
90    format('***************************************************************************')
91    format('*                          REACTION SYSTEMS                               *')

      write(*,100) namep,massp,zp,jp
100   format('Project=',A5,' MASS = ',F8.4,' Z = ',F5.1, ' Jp = ',f4.1)

      write(*,110)namet,masst,zt,jt
110   format('Target =',A5,' MASS = ',F8.4,' Z = ',F5.1, ' Jt = ',f4.1)

      write(*,120)
120   format('--------------------------------------------------------')

      write(*,130)nameb,massb,zb,jb
130   format('Targeb =',A5,' MASS = ',F8.4,' Z = ',F5.1, ' Jb = ',f4.1)

      write(*,140)namex,massx,zx,jx
140   format('Targex =',A5,' MASS = ',F8.4,' Z = ',F5.1, ' Jx = ',f4.1)

      write(*,90)
      write(*,*)

c***print potential
      write(*,90)
      write(*,160)
      write(*,90)
160   format('*                          DEFINED POTENTIALS                             *')
      write(*,170)
170   format(1X,'KP1#',2x,'KP2#',7X,'TYPE',12X,'SHAPE')

180   format(3X,A1,2x,I2,5x,'Coulomb',10X,'CHARGE',10x,'a1=',F8.3,
     &   '   a2=',
     &       F8.3,'   rc=',F6.3)

190   format(3x,A1,2x,I2,5x,'Volume',10x,I2,'=',A15,"v=",F7.3,'    av=',
     &       F6.3,'   rv=',F6.3,' w=',F7.3, ' aw=',F6.3, ' rw=',
     &        F6.3)

200   format(3x,A1,2x,I2,5x,'Projtl S.O.',5x,I2,'=',A15,'vsov=',F7.3,
     &       '  rsov=',F6.3,'  asov=',F6.3,'  vsow=',F6.3,
     &         '  rsow=',F6.3,'  asow=',F6.3)

210   format(3x,A1,2x,I2,5x,'Surface',9x,I2,'=',A15,'vd=',F7.3,
     &       '  rvd=',F6.3,'  avd=',F6.3,'  wd=',F7.3,
     &         '  rwd=',F6.3,'  awd=',F6.3)

      do counter=1,5
         write(*,180) nkp1(counter),nkp2(counter),na1(counter),
     &         na2(counter),nrc(counter)
         write(*,190) nkp1(counter),nkp2(counter),nptype(counter),
     & potype(nptype(counter)),nuv(0,counter),nav(counter),
     &      nrv(counter),nuw(counter),naw(counter),nrw(counter)
         if (nvsov(counter)/=0 .or. nvsow(counter)/=0 ) then
            write(*,200) nkp1(counter),nkp2(counter),nptype(counter),
     &      potype(nptype(counter)),
     &      nvsov(counter),
     &            nrsov(counter),nasov(counter),nvsow(counter),
     &      nrsow(counter),nasow(counter)
         end if
         if (nvd(counter)/=0 .or. nwd(counter)/=0 ) then
            write(*,210)nkp1(counter),nkp2(counter),nptype(counter),
     &      potype(nptype(counter)),
     &      nvd(counter),nrvd(counter),
     &            navd(counter),nwd(counter),nrwd(counter),nawd(counter)
         end if
         if(counter/=5) write(*,120)
         if(counter==5) write(*,90)
      end do
      write(*,*)


      end subroutine


c-----------------------------------------------------------------------
      function potype(a)
      integer :: a
      character(len=15) :: potype
      select case(a)
      case(1)
         potype='Woods-Saxon'
      case(2)
         potype='Gaussian'
      case(3)
         potype="kd02"
      case(4)
         potype="ch89"
      case(5)
         potype="bgpn"
      case(6)
         potype="yyq06"
      case(8)
         potype="yyh10"
      case(9)
         potype="bg69"
      case(10)
         potype="bg69"
      case(11)
         potype="AW-pot"
      case(12)
         potype="Angela-pot"
      case(13)
         potype="alpha-global"
      case(15)
         potype="Morillon's DOM"
      case(17)
         potype="Guo's triton potential"
      case(18)
         potype="alpha-global by Palumbo"
      case(19)
        potype="3h/3he GDP08"
      case(20)
        potype="alpha potential by Avrigeneau"
      case(22)
        potype="Fit to Palumbo alpha+118Sn with WS2"
      case(41)
         potype="read fort.41"
      case(42)
         potype="read fort.42"
      case(43)
         potype="read fort.43"
      case(44)
         potype="read fort.44"
      end select
      end function

c-----------------------------------------------------------------------
cWrite output files units and names
       subroutine fkind()
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       flkind(:) = ' '
       flkind(1)='local copy of input'
       written(1)=.TRUE.
       flkind(2)='phase-shifts for incoming channel'
       written(2)=.TRUE.
       flkind(3)='S-matrix for incoming channel'
       written(3)=.TRUE.
       flkind(4)='radial part of Wf of channel a'
       written(4)=.TRUE.
       flkind(7)='bound-state wf of channel p'
       written(7)=.TRUE.
       flkind(8)='channels coupling index'
       written(8)=.TRUE.

       flkind(101)='radial part of wf of channel b'
       written(101)=.TRUE.
       flkind(151)='radial part of wf of channel b in grids'
       written(151)=.TRUE.
       flkind(201)='regular part of Gx'
       written(201)=.TRUE.
       flkind(251)='irregular part of Gx'
       written(251)=.TRUE.

       flkind(11)='source term'
       if(printf) written(11)=.TRUE.
       flkind(12)='x-A wave function '
       if(printf) written(12)=.TRUE.

       flkind(13)='No term'
       if(printf) written(13)=.TRUE.
       
       flkind(14)='|\Psi_{xA}| in la'
       written(14)=.FALSE.


       flkind(19)='x-section J_parity distribution'
       written(19)=.FALSE.
       flkind(20)='x-section alpha_xA distribution'
       written(20)=.TRUE.
       flkind(21)='double x-section angular distribution'
       written(21)=.TRUE.
       flkind(22)='x-section energy distribution'
       written(22)=.TRUE.
       flkind(23)='x-section la distribution'
       written(23)=.TRUE.
       flkind(24)='x-section lb distribution'
       written(24)=.TRUE.


       flkind(30)='potential used in the calcuation'
       written(30)=.TRUE.
       flkind(31)='Wx(rx)'
       written(31)=.TRUE.
       flkind(32)='parameters of calculated potentials'
       written(32)=.TRUE.
       flkind(33)='fusion potential'
       
       if(VF) written(215)=.TRUE.
       flkind(215)='integrand in y-direction'


C
C       flkind(14)='x-section energy distribution in cm'
C       flkind(15)='partial differential x-sections'
C       flkind(16)='double differential x-section in cm'
C       flkind(17)='x-sections angular distribution in cm'
C       flkind(20)='NO x-section energy distribution in cm'
C       flkind(21)='partial differential NO x-sections'
C       flkind(22)='double differential NO x-section in cm'
C       flkind(23)='interface dsde in cm'
C       flkind(24)='interface dsdela'
C       flkind(30)='potential used in the calcuation'
C       flkind(31)='real part of kernel q'
C       flkind(32)='imag part of kernel q'
C       flkind(33)='modules of kernel q'
C       flkind(34)='real part of integral function'
C       flkind(35)='imag part of integral function'
C       flkind(36)='modules of integral function'
C       flkind(37)='real part of kernel qNO'
C       flkind(38)='imag part of kernel qNO'
C       flkind(39)='modules of kernel qNO'



C       flkind(43)='regular solution of Gx'
C       flkind(44)='irregular solution of Gx'
C       flkind(45)='radial part of b-bin wf of channel b'
C       flkind(46)='finite range correction Lambda'
C       flkind(51)='rholalblx(r)'
C       flkind(52)='rlalblx(r)'
C       flkind(62)='radial part of wf of channel b'
C       flkind(101)='transfer x-section for bound states'
C       flkind(201)='BF x-section energy distribution in cm'
C       flkind(202)='EB x-section energy distribution in cm'
C       flkind(203)='BF partial differential x-sections'
C       flkind(204)='EB partial differential x-sections'
C       flkind(205)='BF double differential x-section in cm'
C       flkind(206)='EB double differential x-section in cm'
C
C


       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine

C******************************************************************************
C FILE: omp_getEnvInfo.f
C DESCRIPTION:
C   OpenMP Example - Get Environment Information - Fortran Version
C   The master thread queries and prints selected environment information.
C AUTHOR: Blaise Barney  7/06
C LAST REVISED: 07/12/06
C******************************************************************************

      subroutine GETINFO()

      INTEGER NTHREADS, TID, OMP_GET_NUM_THREADS,
     +  OMP_GET_THREAD_NUM, OMP_GET_NUM_PROCS, OMP_GET_MAX_THREADS,
c     +  OMP_IN_PARALLEL, OMP_GET_DYNAMIC, OMP_GET_NESTED,
     +  PROCS, MAXT

C     These are for AIX compilations
C     INTEGER INPAR, DYNAMIC, NESTED
C     These are for non-AIX compilations
c      LOGICAL INPAR, DYNAMIC, NESTED

C     Start parallel region
!$OMP PARALLEL PRIVATE(NTHREADS, TID)

C     Obtain thread number
      TID = OMP_GET_THREAD_NUM()

C     Only master thread does this
      IF (TID .EQ. 0) THEN

        PRINT *, 'Thread',tid,'getting environment information'

C     Get environment information

        PROCS = OMP_GET_NUM_PROCS()
        NTHREADS = OMP_GET_NUM_THREADS()
        MAXT = OMP_GET_MAX_THREADS()
c        INPAR = OMP_IN_PARALLEL()
c        DYNAMIC = OMP_GET_DYNAMIC()
c        NESTED = OMP_GET_NESTED()
C       if (nnodes>NTHREADS  .or. nnodes==0)  nnodes=NTHREADS
C     Print environment information

        PRINT *, 'Number of processors = ', PROCS
        PRINT *, 'Number of threads = ', NTHREADS
        PRINT *, 'Max threads = ', MAXT
c        PRINT *, 'In parallel? = ', INPAR
c        PRINT *, 'Dynamic threads enabled? = ', DYNAMIC
c        PRINT *, 'Nested parallelism supported? = ', NESTED

      END IF

C     Done
!$OMP END PARALLEL

      END subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_cm2lab_input()
c     Generate cm2lab.in file for center-of-mass to lab transformation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      implicit none
      
      real*8 :: m1, m2, m3, m4
      real*8 :: ecmfmin, ecmfmax, ecmfh
      real*8 :: thcmmin, thcmmax, thcminc
      real*8 :: elabmin, elabmax, elabh
      real*8 :: thlabmin, thlabmax, thlabinc
      real*8 :: ptheta, pelab
      real*8 :: the_int_min, the_int_max
      integer :: nefac, lmax_local
      character*20 :: filecm
      
      write(*,*) 'Generating cm2lab.in file...'
      
      m1 = massp
      m2 = masst  
      m3 = massb
      m4 = masst + massx
      
      ! Note: elab is already available from systems module
      
      ecmfmin = necmb(1)
      ecmfmax = necmb(neb)
      if (neb > 1) then
         ecmfh = necmb(2) - necmb(1)
      else
         ecmfh = 1.0d0
      endif
      
      thcmmin = thmin
      thcmmax = thmax
      thcminc = thinc
      
      ! Set lmax for cm2lab (use 0 to avoid reading fort.20 file)
      lmax_local = min(lxmax, lmax)
      
      ! Set optional parameters
      the_int_min = -99.0d0
      the_int_max = -99.0d0
      nefac = 1
      filecm = ""  ! Empty means read from stdin
      
      elabmin = 1.0d0
      elabmax = 31.0d0
      elabh = 2.0d0
      
      thlabmin = 1.0d0
      thlabmax = 180.0d0
      thlabinc = 1.0d0
      
      ptheta = 2.0d0
      pelab = 2.0d0
      
      open(unit=99, file='cm2lab.in', status='replace')
      
      write(99,*) 'NAMELIST'
      write(99,*) '&cmsys  m1=',m1,' m2=',m2,' m3=',m3,' m4=',m4,
     &            ' elab=',elab,' ecmfmin=',ecmfmin,
     &            ' ecmfmax=',ecmfmax,' ecmfh=',ecmfh,
     &            ' thcmmin=',thcmmin,' thcmmax=',thcmmax,
     &            ' thcminc=',thcminc,' lmax=',lmax_local,' /'
      write(99,*) '&labsys      elabmin=',elabmin,' elabmax=',elabmax,
     &            ' elabh=',elabh,' thlabmin=',thlabmin,
     &            ' thlabmax=',thlabmax,' thlabinc=',thlabinc,
     &            ' ptheta=',ptheta,' pelab=',pelab,' /'
      
      close(99)
      
      write(*,*) 'cm2lab.in file generated successfully!'
      write(*,*) 'Parameters:'
      write(*,100) m1,m2,m3,m4
      write(*,110) elab
      write(*,120) ecmfmin,ecmfmax,ecmfh
      write(*,130) thcmmin,thcmmax,thcminc
      write(*,140) lmax_local
      
100   format('  Masses: m1=',F8.4,' m2=',F8.4,' m3=',F8.4,' m4=',F8.4)
110   format('  Lab energy:',F8.2,' MeV')
120   format('  CM energy range:',F8.2,'-',F8.2,' MeV (step:',F6.2,')')
130   format('  CM angle range:',F6.1,'-',F6.1,' deg (step:',F4.1,')')
140   format('  lmax:',I3)
      
      end subroutine




      end module
