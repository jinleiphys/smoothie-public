cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program cm2lab
c     this code is the used to transform the results of nonelastic
c     breakup cross section from center of mass frame to lab frame
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      call cmch2cmpar()
      call labsec()
      
      write(*,*) ''
      write(*,*) '=================================================='
      write(*,*) '               CALCULATION COMPLETED'
      write(*,*) '=================================================='
      write(*,*) ''
      write(*,*) 'Output Files Generated:'
      write(*,*) '--------------------------------------------------'
      write(*,*) ' fort.910  | Double cross sections in CM frame'
      write(*,*) ' fort.911  | Double cross sections in lab frame'
      write(*,*) ' fort.912  | Cross section energy distribution'
      write(*,*) ' fort.913  | Cross section angular distribution 1'
      write(*,*) ' fort.914  | Cross section angular distribution 2'
      write(*,*) ' fort.915  | Double cross sections for fixed angle',
     &           ' (ptheta)'
      write(*,*) ' fort.916  | Double cross sections for fixed energy',
     &           ' (pelab)'
      write(*,*) ' fort.917  | Double cross sections 3-D plot (lab)'
      write(*,*) ' fort.918  | Double cross sections 3-D plot (CM)'
      write(*,*) ' fort.919  | Cross sections lx distribution (CM)'
      write(*,*) ' fort.920  | Cross sections for selected angles (CM)'
      write(*,*) ' fort.921  | CM angular distribution'
      write(*,*) '--------------------------------------------------'
      write(*,*) ''
      write(*,*) 'Analysis complete. Check output files for results.'
      write(*,*) '=================================================='
      end program
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module systems
      real*8 :: m1,m2,m3,m4 ! mass for a A b B
c     IAV model: m1=massa, m2=masst, m3=massb, m4=massB*
      real*8 :: ecmi
      real*8 :: elabh,ecmfh,ecmh
      real*8,allocatable,dimension(:) :: ecmf ! energy for channel b
      real*8,allocatable,dimension(:) :: ecmb ! energy of particle b
      real*8,allocatable,dimension(:) :: elabb ! lab energy of particle b
      real*8,allocatable,dimension(:) :: thcm,thlab
      real*8 :: thlabmax,thlabmin,thlabinc
      integer :: necm,nelab,nthcm,nthlab ! maximum value used for allocating the matrix
      real*8 :: ptheta, pelab
      real*8:: amu=931.49432    ! MeV
      end module

      module xsec
      real*8,allocatable,dimension(:,:) :: xsecmf ! double differential x-section of channel energy in CM
      real*8,allocatable,dimension(:,:) :: xsecmb ! double differential x-section of b in CM
      real*8,allocatable,dimension(:,:) :: xyt ! for two dimension interpolation
      real*8,allocatable,dimension(:,:) :: xseclabb ! double differential x-section of b in Lab
      real*8,allocatable,dimension(:,:) :: xsec_laba_cme,dsdelx
      end module


      module r_grid
      real*8, allocatable, dimension(:) :: rr, rw
      end module
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmch2cmpar()
c     this subroutine is used to transform channel energy to particle
c     energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use xsec
      use r_grid
      implicit none
      real*8 :: elab
      real*8 :: ecmfmin,ecmfmax
      real*8 :: thcmmin,thcmmax,thcminc
      real*8 :: elabmin,elabmax
      real*8 :: pi
      real*8,dimension(1:100) :: quadw,quadx
      real*8,allocatable,dimension(:) :: sumsec
      real*8 :: xm,xl,fival,f0,integral
      integer :: i,j,l
      real*8 :: f1,f2,f3
      integer  :: r1,r2,r3
      integer :: lmax,ll,lx=0,lread
      character*80 :: line
      character*20 :: filecm
      real*8 :: iecm,jthcm, cmxsec! used for interpolation
      real*8 :: ex,xsec_int
      real*8 :: the_int_min, the_int_max,dsdw
      integer :: n_theta_int,nefac

c      ptheta=0.0d0
c      pelab=0.0d0
C     lmax=0
      namelist /cmsys/ m1,m2,m3,m4,elab,ecmfmin,ecmfmax,ecmfh,thcmmin,
     &                  thcmmax,thcminc,lmax, the_int_min, the_int_max,
     &                  nefac,filecm

      namelist /labsys/ elabmin,elabmax,elabh,thlabmax,thlabmin,
     &                  thlabinc,ptheta,pelab

      filecm=""
      ptheta=0.0d0
      pelab=0.0d0
      lmax=0
      the_int_min=-99
      the_int_max=-99

      write(*,*) ''
      write(*,*) '=================================================='
      write(*,*) '           CM2LAB v1.0'
      write(*,*) '   Center of Mass to Lab Frame Converter'
      write(*,*) '   Nonelastic Breakup Cross Section Analysis'
      write(*,*) '=================================================='
      write(*,*) ''
      write(*,*) 'IAV Model Mass Parameters:'
      write(*,*) '  m1 = massa  (projectile mass)'
      write(*,*) '  m2 = masst  (target mass)'
      write(*,*) '  m3 = massb  (breakup fragment mass)'
      write(*,*) '  m4 = massB* (excited target mass)'
      write(*,*) ''
      write(*,*) 'Angular Momentum Parameter:'
      write(*,*) '  lmax = lxmax (maximum angular momentum between',
     &           ' x and A)'
      write(*,*) ''
      write(*,*) 'Input File Structure:'
      write(*,*) '--------------------------------------------------'
      write(*,*) 'The input file must contain two namelists:'
      write(*,*) ''
      write(*,*) '&cmsys'
      write(*,*) '  m1, m2, m3, m4     = masses (amu)'
      write(*,*) '  elab               = laboratory energy (MeV)'
      write(*,*) '  ecmfmin, ecmfmax   = CM energy range (MeV)'
      write(*,*) '  ecmfh              = CM energy step (MeV)'
      write(*,*) '  thcmmin, thcmmax   = CM angle range (degrees)'
      write(*,*) '  thcminc            = CM angle step (degrees)'
      write(*,*) '  lmax               = max angular momentum in',
     &           ' x-A channel'
      write(*,*) '/'
      write(*,*) ''
      write(*,*) '&labsys'
      write(*,*) '  elabmin, elabmax   = lab energy range (MeV)'
      write(*,*) '  elabh              = lab energy step (MeV)'
      write(*,*) '  thlabmin, thlabmax = lab angle range (degrees)'
      write(*,*) '  thlabinc           = lab angle step (degrees)'
      write(*,*) '  ptheta             = angle for fixed-angle cross',
     &           ' section (degrees)'
      write(*,*) '  pelab              = energy for fixed-energy cross',
     &           ' section (MeV)'
      write(*,*) '/'
      write(*,*) '--------------------------------------------------'
      write(*,*) '=================================================='
      write(*,*) 'Reading input parameters from namelist...'

      read(5,nml=cmsys)
      read(5,nml=labsys)
      pi=acos(-1.0d0)
      ecmi=elab*m2/(m1+m2)
      necm=nint((ecmfmax-ecmfmin)/ecmfh)
      nthcm=nint((thcmmax-thcmmin)/thcminc)
      nelab=nint((elabmax-elabmin)/elabh)
      nthlab=nint((thlabmax-thlabmin)/thlabinc)

      write(*,*) "ecmi=",ecmi
      write(*,*) "necm=",necm
      write(*,*) "nthcm=",nthcm
      write(*,*) "nelab=",nelab
      write(*,*) "nthlab=",nthlab


      allocate(ecmf(0:necm),ecmb(0:necm),elabb(0:nelab))
      allocate(xyt(1:2,1:max(necm+1,nthcm+1)))
      allocate(thcm(0:nthcm),thlab(0:nthlab))
      allocate(xsecmf(0:necm,0:nthcm),xsecmb(0:necm,0:nthcm))
      allocate(xseclabb(0:nelab,0:nthlab))
      allocate(sumsec(0:nthcm))
      allocate(xsec_laba_cme(0:nthlab,0:necm))
      if(lmax/=0) allocate(dsdelx(0:lmax,0:necm))


c***
      if (filecm.ne."") then 
         write(*,*)'Reading double x-sections from file: ', filecm
         open(21,file=filecm,status='old')
      endif
      ecmh=ecmfh*m4/(m3+m4)
      do i=0,necm
      ecmf(i)=ecmfmin+i*ecmfh
      ecmb(i)=m4*ecmf(i)/(m3+m4)
C     write(*,*)"ecmf=",ecmf(i)
C     write(*,*)"ecmb=",ecmb(i)
         do j=0,nthcm
            read(21,*)thcm(j),xsecmf(i,j)
            xsecmb(i,j)=xsecmf(i,j)*(m3+m4)/m4
         end do
         read(21,*) line

!! read for lx distribution
       if(lmax/=0) then
         do lx=0,lmax
           read(20,*) lread,dsdelx(lx,i)
         end do
          read(20,*) line
       end if

10     continue
      end do

! write for lx distribution
       if (lmax/=0) then
         do lx=0, lmax
            do i=0,necm
               write(919,*)ecmfmin+i*ecmfh,dsdelx(lx,i)
            end do
            write(919,*) "&"
         end do

       end if




! AMM: Write double differential cross section in CM
      write(910,*) necm+1,nthcm+1
      do j=0,nthcm
      dsdw=0
      write(910,'(a,1f8.4)') '#thcm=',thcm(j)
      do i=0,necm
      write(910,'(1f10.4,1g15.7)') ecmb(i), xsecmb(i,j)
      dsdw=dsdw +  xsecmb(i,j)*ecmh
      enddo ! energy
        write(921,*) thcm(j),dsdw
      enddo ! angle



        xyt(1,1:necm+1) =  ecmb(0:necm)
        xyt(2,1:nthcm+1) = thcm(0:nthcm)


      do i=-20,necm+20
         iecm=ecmfmin+i*ecmfh
         do j=0,nthcm

           if (i>=0 .and. i<=necm ) then
              write(918,*)thcm(j),iecm,xsecmf(i,j)
           else
              write(918,*)thcm(j),iecm,xsecmf(i,j)
           end if
         end do
          write(918,*) ' '
      end do
  

!!!!! integrated the cm frame cross section for select angles
!!!!!  use simpson integral method


      if (the_int_max > 0) then
      n_theta_int= nint( (the_int_max-the_int_min)/ thcminc ) + 1
      allocate(rr(1:n_theta_int), rw(1:n_theta_int))

      call simpson(n_theta_int,the_int_min,the_int_max,rr,rw)

      do i=0, necm
      iecm=ecmfmin+i*ecmfh
      xsec_int=0.0d0
      do j =1, n_theta_int
         xsec_int=xsec_int+ 2.*pi * sin( rr(j)*pi/180. ) *
     +            (rw(j)*pi/180.) * xsecmf(i,j)
      end do

      write(920,*) iecm,xsec_int


      end do

      deallocate(rr,rw)


      end if


c***
      do i=0,nelab
      elabb(i)=elabmin+i*elabh
      end do
c***
      do i=0,nthlab
      thlab(i)=thlabmin+i*thlabinc
      end do



c      do i=0,necm
c         do j=0,nthcm
c            write(13,*)thcm(j),xsecmb(i,j)
c         end do
c      write(13,*)'&'
c      end do


      sumsec=0.0d0
c      do j=0,nthcm
c         do i=0,necm
c            sumsec(j)=sumsec(j)+xsecmf(i,j)
c         end do
c         sumsec(j)=sumsec(j)*ecmfh
c      end do



      do j=0,nthcm
         r1=0;r2=1;r3=2
         do
            if (r2>necm) exit
            f1=xsecmf(r1,j)
            f2=xsecmf(r2,j)
            f3=xsecmf(r3,j)
            sumsec(j)=sumsec(j)+(ecmfh/3.0)*(f1+4*f2+f3)
            r1=r1+2;r2=r2+2;r3=r3+2
         end do
      end do


      integral=0.0d0
      call gauleg(100,quadx,quadw)
      xm=(thcmmax+thcmmin)*0.5*pi/180.
      xl=(thcmmax-thcmmin)*0.5*pi/180.

      do i=1,100
      f0=sin((xm + xl*quadx(i)))*fival((xm + xl*quadx(i))*180./pi,thcm,
     &             sumsec,nthcm+1,1.0d0)*2*pi
      integral=integral+quadw(i)*f0
      end do
      integral=integral*xl
      write(*,*) 'inelastic breakup X-sec before transformation is',
     &              integral
      end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine labsec()
c     this subroutine is used to calculate the lab x-sections
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use xsec
      implicit none
      real*8 :: iecm,jthcm ! used for lab2cm
      real*8 :: ijcmsec ! interpolate results
      real*8 :: pi,sinthcm,sinthlab
      real*8,dimension(0:nthlab) :: sumsec,sumsec2
      real*8,dimension(1:100) :: quadw,quadx
      real*8 :: fival,f0,xm,xl,integral
      real*8,dimension(0:nelab) ::elabsec
      integer :: i,j
      real*8 :: f1,f2,f3
      integer  :: r1,r2,r3
      real*8 :: beta
      real*8,dimension(0:nthcm) :: xsec_lab,theta_lab
      real*8 :: fint2db,fint2db2,f2c,ffr4 ! 2-dimension interpolation
      real*8 :: gamma,jacobian ! for theta=180 handling
      pi=acos(-1.0d0)



         do i=0,nelab
           do j=0,nthlab
            call lab2cm(elabb(i),thlab(j),iecm,jthcm)
c            ijcmsec=fint2db(ecmb,thcm,xsecmb,iecm,jthcm,necm+1,nthcm+1,
c     &                      necm+1)



             ijcmsec=fint2db2(ecmb,thcm,xsecmb,iecm,jthcm,necm+1,
     +       nthcm+1,1.0d0)






            sinthcm=sin(jthcm*pi/180.)
            sinthlab=sin(thlab(j)*pi/180.)
C           write(*,*)"sinthlab=",sinthlab
            
            ! Special handling for theta_lab = 0 or 180 degrees
            if (abs(sinthlab) < 1.0d-8) then
               ! At theta_lab = 0 or 180, use limiting behavior
               if (abs(thlab(j)) < 1.0d-8) then
                  ! Forward scattering: theta_lab = 0
                  ! At theta_lab = 0, theta_cm = 0, so sin(theta_cm) = 0
                  ! Use L'Hopital's rule limit: lim(sin(theta_cm)/sin(theta_lab)) = 1
                  xseclabb(i,j) = ijcmsec
               else
                  ! Backward scattering: theta_lab = 180
                  ! For theta_lab = 180, we need the proper Jacobian from kinematics
                  ! At theta_lab = 180, the transformation involves gamma factor
                  gamma=(m1*m3*ecmi/m2/(m1+m2)/elabb(i))**0.5
                  ! The Jacobian at theta_lab = 180 is (1 + gamma)^2
                  ! This comes from the kinematic transformation
                  jacobian = (1.0d0 + gamma)**2
                  xseclabb(i,j) = ijcmsec * jacobian
               end if
            else
               ! Normal case
               xseclabb(i,j)=sinthcm*ijcmsec/sinthlab
            end if
C            write(*,*)"xseclabb(i,j)=",xseclabb(i,j)
C           stop
            write(911,*) thlab(j), xseclabb(i,j)
            if (abs(thlab(j)-ptheta)<0.0001) then                                                  ! fixed angle
               write(915,*)elabb(i),xseclabb(i,j)
            end if
            if (abs(elabb(i)-pelab)<0.0001) then                                                  ! fixed angle
               write(916,*)thlab(j),xseclabb(i,j)
            end if
            write(917,*)thlab(j),elabb(i),xseclabb(i,j)
         end do
          write(911,*)'&elab=',elabb(i)
         write(917,*) ' '
      end do
      sumsec=0.0d0
c      do j=0,nthlab
c         do i=0,nelab
c            sumsec(j)=sumsec(j)+xseclabb(i,j)
c         end do
c         sumsec(j)=sumsec(j)*elabh
c         write(12,*) thlab(j), sumsec(j)
c      end do

      do j=0,nthlab
         r1=0;r2=1;r3=2
         do
            if (r2>nelab) exit
            f1=xseclabb(r1,j)
            f2=xseclabb(r2,j)
            f3=xseclabb(r3,j)
            sumsec(j)=sumsec(j)+(elabh/3.0)*(f1+4*f2+f3)
            r1=r1+2;r2=r2+2;r3=r3+2
         end do
         write(913,*) thlab(j), sumsec(j)
      end do

c************* one dimension cm 2 lab
      do i=0,necm
      xsec_lab=0.0d0
      theta_lab=0.0d0

         do j=0,nthcm

            call cmth2labth(ecmb(i),thcm(j),theta_lab(j),beta)
            xsec_lab(j)=beta*xsecmb(i,j)

         end do


         do j=0,nthlab
             xsec_laba_cme(j,i)=fival(thlab(j),theta_lab,
     &             xsec_lab,nthcm+1,1.0d0)

         end do
C          stop

!         write(100,*)"&"
       end do


      sumsec2=0.0d0

      do j=0,nthlab
         r1=0;r2=1;r3=2
         f1=0.0d0
         f2=0.0d0
         f3=0.0d0
         do
            if (r2>necm) exit
            f1=xsec_laba_cme(j,r1)
            f2=xsec_laba_cme(j,r2)
            f3=xsec_laba_cme(j,r3)
            if(r3>necm) f3=0.0d0
            sumsec2(j)=sumsec2(j)+(ecmh/3.0)*(f1+4*f2+f3)
!            if(sumsec2(j)>1000) then
!               write(*,*) f1,f2,f3
!            end if
            r1=r1+2;r2=r2+2;r3=r3+2
         end do
         write(914,*) thlab(j), sumsec2(j)
      end do

c******************************************************************



      integral=0.0d0
      call gauleg(100,quadx,quadw)
      xm=(thlabmax+thlabmin)*0.5*pi/180.
      xl=(thlabmax-thlabmin)*0.5*pi/180.

      do i=1,100
      f0=sin((xm + xl*quadx(i)))*fival((xm + xl*quadx(i))*180./pi,thlab,
     &             sumsec,nthlab+1,1.0d0)*2*pi
      integral=integral+quadw(i)*f0
      end do
      integral=integral*xl
      write(*,*) 'inelastic breakup X-sec after transformation is',
     &              integral


      elabsec=0.0d0
c      f1=0.0d0
      do i=0,nelab

         do j=1,100
         f0=sin((xm + xl*quadx(j)))*fival((xm + xl*quadx(j))*180./pi,
     &             thlab,xseclabb(i,0:nthlab),nthlab+1,1.0d0)*2*pi
         elabsec(i)=elabsec(i)+quadw(j)*f0
         end do
         elabsec(i)=elabsec(i)*xl
         write(912,*)elabb(i),elabsec(i)
c         f1=f1+elabsec(i)
      end do
c       f1=elabh*f1
c       write(*,*) 'F1=',f1



      end subroutine
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lab2cm(ielab,jthlab,iecm,jthcm)
c     this subroutine is used to transform from lab to C.M.
c     Special handling for theta_lab=0 and theta_lab=180
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems,only :m1,m2,m3,m4,ecmi
      implicit none
      real*8 :: ielab,jthlab,iecm,jthcm,pi
      real*8 :: beta,gamma,sinlab,coslab,tanthcm,sincm
      real*8 :: eps
      
      pi=acos(-1.0d0)
      eps=1.0d-8  ! small number to avoid singularities
      
      ! Special handling for theta_lab = 0 degrees (forward scattering)
      if (abs(jthlab) < eps) then
         jthcm = 0.0d0
         beta = 1.0d0
         iecm = ielab
         return
      end if
      
      ! Special handling for theta_lab = 180 degrees (backscattering)  
      if (abs(jthlab - 180.0d0) < eps) then
         jthcm = 180.0d0
         gamma=(m1*m3*ecmi/m2/(m1+m2)/ielab)**0.5
         beta = (1.0d0 + gamma)**2
         iecm = ielab/beta/beta
         return
      end if
      
      ! Normal case - original algorithm
      sinlab=sin(jthlab*pi/180.)
      coslab=cos(jthlab*pi/180.)
      gamma=(m1*m3*ecmi/m2/(m1+m2)/ielab)**0.5
      
      ! Check for potential division by zero
      if (abs(coslab-gamma) < eps) then
         ! Use approximation near singular point
         jthcm = 90.0d0  ! perpendicular scattering
         beta = sqrt(1.0d0 + gamma**2)
      else
         tanthcm=sinlab/(coslab-gamma)
         jthcm=atan(tanthcm)
         if (jthcm < 0.0d0) jthcm=jthcm+pi
         jthcm=jthcm*180.0d0/pi
         
         sincm=sin(jthcm*pi/180.)
         if (abs(sinlab) < eps) then
            beta = 1.0d0
         else
            beta=sincm/sinlab
         end if
      end if
      
      if (abs(beta) < eps) beta = eps  ! avoid division by zero
      iecm=ielab/beta/beta

      end subroutine



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmth2labth(iecm,jthcm,jthlab,beta)
c     this subroutine is used to transform from C.M. to lab
c     Special handling for theta_cm=0 and theta_cm=180
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems,only :m1,m2,m3,m4,ecmi
      implicit none
      real*8 :: jthlab,iecm,jthcm,pi
      real*8 :: beta,gamma,sincm,coscm,tanthlab
      real*8 :: eps
      
      pi=acos(-1d0)
      eps=1.0d-8  ! small number to avoid singularities
      
      ! Special handling for theta_cm = 0 degrees (forward scattering)
      if (abs(jthcm) < eps) then
         jthlab = 0.0d0
         beta = 1.0d0
         return
      end if
      
      ! Special handling for theta_cm = 180 degrees (backscattering)
      if (abs(jthcm - 180.0d0) < eps) then
         jthlab = 180.0d0
         gamma=(m1*m3*ecmi/m2/(m1+m2)/iecm)**0.5
         beta = (1.0d0 + gamma)**2
         return
      end if
      
      ! Normal case - original algorithm
      sincm=sin(jthcm*pi/180.)
      coscm=cos(jthcm*pi/180.)
      gamma=(m1*m3*ecmi/m2/(m1+m2)/iecm)**0.5
      
      ! Check for potential division by zero
      if (abs(coscm+gamma) < eps) then
         ! Use approximation near singular point
         jthlab = 90.0d0  ! perpendicular scattering
         beta = sqrt(1.0d0 + gamma**2)
      else
         tanthlab=sincm/(coscm+gamma)
         jthlab=atan(tanthlab)
         if (jthlab < 0.0d0) jthlab=jthlab+pi
         jthlab=jthlab*180.0d0/pi
         
         ! Calculate beta with protection against division by zero
         if (abs(1+gamma*coscm) < eps) then
            beta = 1.0d0
         else
            beta=(1+gamma**2+2*gamma*coscm)**1.5/abs(1+gamma*coscm)
         end if
      end if
      
      if (abs(beta) < eps) beta = eps  ! avoid division by zero
      end subroutine


c *** ------------------------------------------------------
c ***  Simple 2-dim interpolation
c ***  ------------------------------------------------------
      function fint2db(xtab,ytab,fxytab,xbar,ybar,nnx,nny,mmx)
      implicit none
!      implicit real*8 (a-h,o-z)
      integer ix,iy,nnx,nny,mmx,ixref,iyref
      real*8 xtab(nnx),ytab(nny)
      real*8 fxytab(mmx,nny),xbar,ybar,fint2db
      real*8 p1,p2,p3,p4,u,t

!      write(95,'(100e12.4)') xtab(1:neset)
!      write(96,'(100e12.4)') fxytab(1:nnx,2)
!      return
      fint2db=0d0


      if(xbar>xtab(nnx) .or. xbar<xtab(1) .or. ybar>ytab(nny)
     +                  .or. ybar<ytab(1) ) return

      do iy=1,nny-1
         iyref=iy
         if((ybar.ge.ytab(iy)).and.(ybar.le.ytab(iy+1)))goto 120
      enddo


 120  do ix=1,nnx-1
         ixref=ix
         if ((xbar.ge.xtab(ix)).and.(xbar.le.xtab(ix+1)))goto 100
      enddo
!      write(*,*)'fint2db error'
      return
 100  p1=fxytab(ixref,iyref)
      p2=fxytab(ixref+1,iyref)
      p3=fxytab(ixref,iyref+1)
      p4=fxytab(ixref+1,iyref+1)



      t=(xbar-xtab(ixref))/(xtab(ixref+1)-xtab(ixref))
      u=(ybar-ytab(iyref))/(ytab(iyref+1)-ytab(iyref))

!      fint2db=(p1+p2+p3+p4)/4d0
c bilinear interpolation
       fint2db=(1-t)*(1-u)*p1+t*(1-u)*p2+t*u*p3+(1-t)*u*p4
!      write(96,'(10f12.4)') ytab(iyref),ybar,ytab(iyref+1)
!      write(96,'(3f12.4)') xtab(ixref),xbar,xtab(ixref+1)
!      write(96,'(2i4,5e12.4/)')ixref,iyref,p1,p2,p3,p4,fint2db
      end





************************************************************************
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      FUNCTION fival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END

c-----------------------------------------------------------------------

      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END




      FUNCTION FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      END




! two dimension interpolation function
! based on the method of fival function
! f(xbar,ybar) from f(x,y)
      function fint2db2(xtab,ytab,fxytab,xbar,ybar,nnx,nny,alpha)
        implicit none
        integer :: nnx,nny ! dimesnion of x and y
        real*8 ::  xtab(nnx),ytab(nny)  ! vector of x and y
        real*8 :: fxytab(nnx,nny) ! original function
        real*8 :: xbar,ybar,fint2db2 ! interpolted points and value
        real*8 :: alpha
        real*8 :: fival2d(nny)
        real*8 :: fival

        fint2db2=0d0


        if(xbar>xtab(nnx) .or. xbar<xtab(1) .or. ybar>ytab(nny)
     +                  .or. ybar<ytab(1) ) return


        call fival2(xbar,xtab,fxytab,nnx,nny,alpha,fival2d)

        fint2db2 = fival(ybar,ytab,fival2d,nny,alpha)


      end function



************************************************************************
*     subroutine modified for the two dimension interpolation
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      subroutine fival2(r,xv,fdis,ndm,ndm2,alpha,fival2d)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 fdis(ndm,ndm2),y1(ndm2),y2(ndm2),y3(ndm2),y4(ndm2)
      DIMENSION xv(ndm)
      real*8 fival2d(ndm2)  !ndm2 stands for the other dimension of 2-dimension interpolation
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0,:)
      y2=fdis(nst+1,:)
      y3=fdis(nst+2,:)
      y4=fdis(nst+3,:)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival2d=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival2d=fdis(ndm,:) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END subroutine
