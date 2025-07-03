      module bound
      real*8,dimension(:,:),allocatable :: phi_bx
      real*8,dimension(:,:),allocatable :: phi_xA
      logical :: precheck
      real*8 :: tau
      integer,dimension(10) :: a2bindex=-99,nn
      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chanbx()
c     informations of channel bx
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:hcm,irmatch,rr,rrw,nr
       use systems,only:zb,massb,zx,massx,be,nodes,gswf
       use constants,only:amu,hbarc,e2,zero
       use precision
       use pot
       use interpolation
       use channels
       use whittaker
       use input,only:written,surfacebw,rs
       implicit none
       integer :: ir,i,irmid
       real*8 :: mu,ecm
       real*8 :: const,r
       real*8 :: norm !  normalization parameter
       real*8 :: fl !  matching point wave function
       real*8,dimension(0:irmatch) :: fin
       real*8,dimension(1:irmatch) :: kl
       real*8 :: f0,psivpsi
       real*8 :: eta
       integer :: IE
       real*8,dimension(0:lbx) :: WK,WKD
       real*8 :: kbx
       real*8 :: fpin,fpout
       real*8 :: delta
       real*8 :: n,m
       real*8,dimension(0:irmatch) :: Vc,Vn,Vn1
       real*8 :: a1,a2,a13,rc,uv
       integer :: counter
       integer :: nco
       real*8 :: f1,f2
       real*8 :: ls,s,j
       real*8 :: fff,anc
       integer :: l,nch
       real*8 :: sqrtmean_r !from:Junzhe
       real*8 :: mean_r
      !  real*8 :: delta_r

      !  complex*16, dimension(1:irmatch) :: vvv

      printpot=.true.

       write(*,*)' ----------------------'
       write(*,*)' PROJECTILE BOUND STATE'
       write(*,*)' ----------------------'

        a1=0.0_dpreal;a2=0.0_dpreal;rc=0.0_dpreal;uv=0.0_dpreal
        do counter=1,99
          if (nkp1(counter)=='p') then
            if(nkp2(counter)==1) then
               a1=na1(counter);a2=na2(counter);rc=nrc(counter)
               uv=nuv(0,counter)
               exit
            end if
         end if
       end do

       if (allocated(phi_bx)) deallocate(phi_bx)
       if (allocated(vbx)) deallocate(vbx)
       allocate(phi_bx(0:irmatch,in2b%nchmax))
       allocate(vbx(0:irmatch,1:in2b%nchmax))
       phi_bx=0.0_dpreal
       vbx=0.0_dpreal


       a13=a1**(1./3.)+a2**(1./3.)

!! this part should be modified for mulitichannel case
       ecm=be
       mu=amu*(massx*massb)/(massb+massx)
       kbx=sqrt(-2*mu*ecm/(hbarc**2))
       eta=zx*zb*e2*mu/hbarc/hbarc/kbx

! Bound-state wf from external file (ONLY PRIOR FORM!)
       if (gswf.ne."") then
            call readwf(gswf,phi_bx,irmatch,hcm,in2b%nchmax)
            written(7)=.true.
       	    do ir=0,irmatch
       		write (7,*) hcm*ir, phi_bx(ir,1)
       	    end do
            write(7,*)"&"
            return
       endif


       do ir=0,irmatch ! for coulomb potential
         r=ir*hcm
         Vc(ir)=vcoul(r,zx*zb,rc*a13)
       end do

       do nch=1,in2b%nchmax

          j=in2b%j(nch)
          l=in2b%l(nch)
          s=in2b%s(nch)
          ls=-0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
          call potr('p',1,zx*zb,ls,0)
          vbx(0:irmatch,nch)=v(0:irmatch)
          Vn1=vbx(0:irmatch,nch)-Vc




      !!!!!! here just for test, one should remove later 
C       vvv=vbx(1:irmatch,1)
C       call thobasis(mu)
C       call tho_eigenvalue(mu,vvv)
      !!!!! test block 

          m=1.0d0
          irmid=nint(2.0_dpreal/hcm) ! arbitrary radius

          do  ! for nodes
             delta=0.0_dpreal
             n=m
             do !  for potential
                n=n*(1+delta)
                Vn(0:irmatch)=n*Vn1(0:irmatch)
                vbx(0:irmatch,nch)=Vc+Vn
c****
c      Numerov method to solve the differential radial equation
c      from zero
       phi_bx(0,nch)=0      ! boundary condition
       phi_bx(1,nch)=hcm**(l+1) ! arbitrary value

       ir=1; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*vbx(ir,nch)/hbarc**2
       phi_bx(2,nch)=2.*phi_bx(1,nch)-hcm**2*kl(1)*phi_bx(1,nch)

       ir=2; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*vbx(ir,nch)/hbarc**2
       const=hcm**2/12.

       do ir=2 ,irmid+1
        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-2.*mu*vbx(ir+1,nch)/hbarc**2
        phi_bx(ir+1,nch)=((2.-10.*const*kl(ir))*phi_bx(ir,nch)-(1.+const*kl(ir-1))*phi_bx(ir-1,nch))/(1.+const*kl(ir+1))
       end do !ir

       fl=phi_bx(irmid,nch)
       fpin=(-phi_bx(irmid+2,nch)+8.*phi_bx(irmid+1,nch)-8.*phi_bx(irmid-1,nch)+phi_bx(irmid-2,nch))/12./hcm
       fin(0:irmid)=phi_bx(0:irmid,nch)

c      from infinity
        IE=0

        call WHIT(eta,hcm*(irmatch),kbx,ecm,l,WK,WKD,IE)
        phi_bx(irmatch,nch)=WK(l)
        fff=WK(l)

        call WHIT(eta,hcm*(irmatch-1),kbx,ecm,l,WK,WKD,IE)
        phi_bx(irmatch-1,nch)=WK(l)

        ir=irmatch;r=ir*hcm
        kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*vbx(ir,nch)/hbarc**2

       	ir=irmatch-1;r=ir*hcm
       	kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*vbx(ir,nch)/hbarc**2

       	do ir=irmatch-1,irmid-2,-1
          kl(ir-1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir-1.)*hcm)**2-2.*mu*vbx(ir-1,nch)/hbarc**2
          phi_bx(ir-1,nch)=((2.-10.*const*kl(ir))*phi_bx(ir,nch)-(1.+const*kl(ir+1))*phi_bx(ir+1,nch))/(1.+const*kl(ir-1))
       	end do   ! ir

       	fpout=(-phi_bx(irmid+2,nch)+8.*phi_bx(irmid+1,nch)-8.*phi_bx(irmid-1,nch)+phi_bx(irmid-2,nch))/12./hcm


        norm=phi_bx(irmid,nch)/fl
        phi_bx(0:irmid,nch)=fin(0:irmid)*norm
        fpin=fpin*norm

        psivpsi=0.0d0
        do i=1,nr
         f0=FFR4(rr(i)/hcm,Vn,irmatch+1)*abs(FFR4(rr(i)/hcm,phi_bx(0:irmatch,nch),irmatch+1))**2
         psivpsi=psivpsi+f0*rrw(i)
        end do ! i
        delta=real(phi_bx(irmid,nch)*(fpout-fpin)/psivpsi)
        if (abs(delta)<1e-6) exit

      end do ! for potential
c control nodes
        nco=1
        do ir=0,irmatch-1
           f1=real(phi_bx(ir,nch))
           f2=real(phi_bx(ir+1,nch))
           if((f1*f2)<0.0d0) then
             nco=nco+1
           end if
        end do
!        write(*,*) "nco=",nco
        write(*,*) " Nb of nodes=",nco

c end control nodes
      if (nco==nodes) exit
      if (nco>nodes) m=m*0.8d0
      if (nco<nodes) m=m*1.2d0
      end do   ! for nodes

C        do ir=irmatch+1,irmatch
C           call WHIT(eta,hcm*(ir),kbx,ecm,lbx,WK,WKD,IE)
C           rwflp(ir)=WK(lbx)
C        end do


c      Normalization

      norm=0.0d0
      do i=1,nr
         norm=norm+abs(FFR4(rr(i)/hcm,phi_bx(0:irmatch,nch),irmatch+1))**2*rrw(i)
      end do
      norm=1.0d0/norm
      phi_bx=phi_bx*sqrt(norm)

      sqrtmean_r=0.0d0
      mean_r=0.0d0
      do i=1,nr
         sqrtmean_r=sqrtmean_r+abs(FFR4(rr(i)/hcm,phi_bx(0:irmatch,nch),irmatch+1))**2*rrw(i)*rr(i)**2
      end do

      do i=1,nr
            mean_r=mean_r+abs(FFR4(rr(i)/hcm,phi_bx(0:irmatch,nch),irmatch+1))**2*rrw(i)*rr(i)
      end do

      write(*,*)"Now print the mean root square of the radius"
      write(*,77)nch,sqrt(sqrtmean_r),sqrt(sqrtmean_r-mean_r**2)
77    format('@nch=',I3,F10.6,F10.6)

c****
      write(*,99)uv,n*uv,n
99    format(2x,"Adjust potential depth from ",F8.3," to ",F8.3,
     +                  " with scaling factor " F7.3)

       anc=phi_bx(irmatch,nch)/fff
       write(*,'(" ANC=",1f10.5)') anc

C       qval=ecm

       if(surfacebw) then
        do ir=1,irmatch
          call WHIT(eta,hcm*(ir),kbx,ecm,l,WK,WKD,IE)
          phi_bx(ir,nch)=WK(l)*anc
        end do
       end if



        if(rs>0.0_dpreal) then
         do ir=1, irmatch
           r=ir*hcm
           phi_bx(ir,nch)=phi_bx(ir,nch)*  exp(-0.693*( exp( (rs - r) /0.1) ) )
         end do
        end if
c-------------------
        written(7)=.true.
       	do ir=0,irmatch
       		write (7,*) hcm*ir, phi_bx(ir,nch)
       	end do
        write(7,*)"&"

       end do
      end subroutine chanbx
c-----------------------------------------------------------------------


! *** Read external WF and interpolates in required grid
!     Format:
!     1st line: custom header in free format
!     2nd line: nb of radial points to be read
!     Next lines: r wf  (two columns, free format)
!
!    Normalization: int |wf(r)|^2
! *** ------------------------------------------------------------------
      subroutine readwf(filename, wf,npt,h,nch)
      implicit none
      integer:: ir,nr,npt,ich,nch
      character*40 filename
      character*80 line,header
      real*8:: r,h, norm
      real*8,parameter:: alpha=0.0
      real*8:: y,fival,wf(0:npt,nch)
      real*8,allocatable:: rv(:), uwf(:)

      open(20,file=filename)
      read(20,'(a80)') header
      read(20,'(a80)') line
      read(line,*)  nr !, rstep,rstart

      if (nch.ne.1) then
        write(*,*)'readwf: nch=',nch,' but only nch=1 implemented'
        stop
      endif

      write(*,'(/," Reading external gs WF from file: ",a40)') filename
      write(*,'(3x,a80)') header
      allocate(rv(nr),uwf(nr))

      do ir=1,nr
         read(20,*) rv(ir),uwf(ir)
      enddo !ir

c interpolate in smoothie grid

!      allocate(wf(0:npt,nch))
      wf(0,:)=0

      do ich=1,nch
      norm=0
      do ir=1,npt
      r=h*ir
      if (r.lt.1e-4) r=1e-4
      if(r.gt.rv(nr)) cycle
      y=fival(r,rv,uwf,nr,alpha)
      norm=norm + y*y*h
      wf(ir,nch)=y
      write(99,*)r,y
      enddo !ich
      write(99,*)'&'
      enddo !ir
      write(*,'(3x,"->  Calculated norm:": 1f8.4,/)') norm

      deallocate(uwf,rv)
      close(20)


      end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pauli_state(counter)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      use gauss
      use precision
      implicit none
      integer :: counter,kin
      integer :: nch_pauli,i
      integer :: rc

      namelist /pauli/ precheck,a2bindex,nn,tau


      kin=5
      rewind(kin)
      precheck=.false.
      tau=1000.0_dpreal
      read(kin,nml=pauli,iostat=rc)
      if(rc/=0) return


      ! finding the number of pauli states
      nch_pauli=0
      do i=1,10
        if(a2bindex(i)>0) nch_pauli=nch_pauli+1
      end do
      write(*,*) "there are",nch_pauli , "need to project out"



      nr1=nr
      deallocate(rr,rrw)
      nr=nint(rmax/hcm)
      allocate(rr(1:nr),rrw(1:nr))
      if (.not. precheck) then
         allocate(phi_xA(1:nr,nch_pauli))
      end if
      call simpson(nr,0.0_dpreal,rmax,rr,rrw)
      call boundxA(counter)
      nr=nr1
      deallocate(rr,rrw)
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0_dpreal,rmax,rr,rrw)


      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundxA(counter)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
      use systems
      use precision
      use mesh
      use channels
      use green_function
      use pot
      implicit none
      real*8 :: E_1,E_2,E
      real*8 :: lambda,lambda1,lambda2
      integer :: l,counter
      real*8 :: s,j,ls
      integer :: alpha_xA,nex,ir
      logical :: nobound,ctrl
      real*8,dimension(1:nr) :: wf
      integer :: ich_pauli,i


      if (allocated(UxA)) deallocate(UxA)
      allocate(UxA(0:irmatch,1:out2b%nchmax))

      write(*,10)
      allocate_k=.false.
      ich_pauli=0
      do alpha_xA=1,out2b%nchmax

        if (.not. precheck) then
          ctrl=.true.
          do i=1,10
            if(alpha_xA == a2bindex(i)) ctrl=.false.
          end do
          if (ctrl) cycle
          if (.not. ctrl) ich_pauli=ich_pauli+1
        end if

        l=out2b%l(alpha_xA)
        s=out2b%s(alpha_xA)
        j=out2b%j(alpha_xA)
        ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
        call potr('x',counter,zx*zt,ls,0)
        UxA(0:irmatch,alpha_xA)=v(0:irmatch)

        write(*,110)l,s,j
        nex=1
        nobound=.false.
        E=0.0_dpreal
        do

          E_1=-0.1_dpreal
          call eigenvalue(E_1,alpha_xA,nex,lambda1,wf)
          E_2=E_1*1.1_dpreal
          call eigenvalue(E_2,alpha_xA,nex,lambda2,wf)

          do while(abs(E_1-E_2)>1.E-4)
            E=E_2-(E_1-E_2)*(lambda2-1.0d0)/(lambda1-lambda2)
            if (E>0 .or. abs(E)<1.e-3) nobound=.true.
            if (E>0 .or. abs(E)<1.e-3) exit
            call eigenvalue(E,alpha_xA,nex,lambda,wf)
            E_1=E_2
            lambda1=lambda2
            E_2=E
            lambda2=lambda
          end do
          if (.not. nobound) then
          do ir=1,nr
          write(109,*)rr(ir),wf(ir)
          end do
          write(109,*)"&"
          end if

          if (nobound) exit

          if(.not. ctrl) then
            if (nex==nn(ich_pauli)) then
               phi_xA(1:nr,ich_pauli)=wf
            end if
          end if

          write(*,220)nex,E
          nex=nex+1
        end do

      end do


10     format('SINGLE-PARTICLE eigenstates:')
110    format('o Configuration: (l=',I2,') (s=',f4.1,') (j=',f4.1,')' )
220    format('    #',I2,'   ebound=',f15.5)
      end subroutine
c-----------------------------------------------------------------------
      subroutine eigenvalue(E,alpha_xA,nex,lambda,wf)
      use mesh
      use channels
      implicit none
      real*8 :: E,lambda
      integer :: n,alpha_xA,nex
      integer :: info,lwork,loc(1),lx
      real*8, dimension(1:nr) :: wr,wi,wf
      real*8, dimension(1:4*nr) :: work
      real*8,dimension(1:nr,1:nr) :: vl,vr,A
      lx=out2b%l(alpha_xA)
      call Amat_bound(E,alpha_xA,A)
! prepare calling the LAPACK routine DGEEV
      lwork=4*nr !  the value is taken from the documentation of DGEEV
! the routine will destroy the original array A
! and return the real parts of the eigenvalues in WR and the imaginary parts in WI
! VR(:,i) contain the eigenstate (=wave function) of the right eigenvector of eigenvalue i
! VL(:,i) contain the eigenstate of the left eigenvector of eigenvalue i
      call DGEEV('N','V',nr,A,nr,WR,WI,VL,nr,VR,nr,WORK,LWORK,INFO)
      if (info/=0) then
         write(*,*) "problem happened when calling DGEEV"
         stop
      end if

! locate maximal eigenvector
! here assume that this eigenvector is real
! then return this eigenvalue in ETA
! and the eigenvector in the global variable WF (which should be of dimension nr)

      do n=1,nex-1
       loc=MAXLOC(wr) ! return value is type of array, must defined with array type
       wr(loc(1))=0
      end do
       loc=MAXLOC(wr)

       lambda=wr(loc(1))
       wf=vr(:,loc(1))
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this subroutine is used to calculate the G0 for negative energy
C     since we need to define G0 as real type
c     for the positive energy, G0 is defined as complex
      subroutine Amat_bound(ecmx,alpha_xA,A)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
      use systems
      use precision
      use mesh
      use pot
      use channels
      use interpolation
      use green_function
      implicit none
      real*8 :: ecmx
      real*8 :: mux
      integer :: ir,lx,irp,alpha_xA
      real*8,dimension(1:nr,1:nr) :: A,G0
      real*8 :: VxA

      mux=amu*(masst*massx)/(massx+masst)
      lx=out2b%l(alpha_xA)
      call G0_func_re_ir(ecmx,lx,mux,G0)
C      call G0_func_fourier(ecmx,lx,mux,G0)

      do ir=1,nr
         do irp=1,nr
           VxA=real(FFC(rr(irp)/hcm,UxA(0:irmatch,alpha_xA),irmatch+1))
           A(ir,irp)=rr(irp)**2 * real(G0(ir,irp)) * VxA * rrw(irp)
         end do
      end do


      end subroutine
c-----------------------------------------------------------------------


      end module



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
