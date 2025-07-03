      module green_function
      real*8,allocatable,dimension(:) :: kk,kw
      integer :: nk
      real*8,allocatable,dimension(:,:,:) :: fl
      logical :: allocate_k
      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this subroutine is used to calculate the G0 for negative energy
C     since we need to define G0 as real type
c     for the positive energy, G0 is defined as complex
      subroutine G0_func_re_ir(ecmx,lx,mux,G0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
      use systems
      use precision
      use mesh
      use pot
      use channels
      use interpolation
      implicit none
      real*8 :: ecmx
      complex*16 :: kx,rho,eta,zlmin,k2,NA
      real*8 :: mux
      integer :: ir,ifail,mode,NL,lx,irp
      complex*16,dimension(0:lmax) :: fc,gc,fcp,gcp,sig
      complex*16,dimension(1:nr,0:lmax) :: Re_func,Ir_func
      real*8,dimension(1:nr,1:nr) :: G0


C      mux=amu*(masst*massx)/(massx+masst)
      k2=2.*mux*ecmx/(hbarc**2)
      kx=sqrt(k2)

      eta=0.0
      zlmin=0.0_dpreal
      NL=lmax+1
      mode=2

      do ir=1,nr
        rho=rr(ir)*kx
        call COULCC(rho,eta,zlmin,NL,fc,gc,fcp,gcp,sig,mode,1,ifail)
        Re_func(ir,:)=fc
        Ir_func(ir,:)=fc+iu*gc
      end do

      NA=-2.0_dpreal*mux*iu*kx/hbarc/hbarc
      do ir=1,nr
         do irp=1,nr
           G0(ir,irp)=real(NA*Re_func(min(irp,ir),lx)*Ir_func(max(irp,ir),lx))
         end do
      end do


      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this subroutine is used to calculate the G0 for negative energy
c     since we need to define G0 as real type
c     for the positive energy, G0 is defined as complex
      subroutine G0_func_fourier(ecmx,lx,mu,G0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      use precision
      use channels
      use gauss
      use constants
      implicit none
      integer :: ir,irp,ik
      integer :: lx
      real*8 :: ecmx,integral,NA,G0k,mu
      real*8,dimension(1:nr,1:nr) :: G0
       NA=2.0_dpreal/pi

      if (.not. allocate_k) then
        nk=80
        allocate(kk(1:nk))
        allocate(kw(1:nk))
        allocate(fl(1:nr,1:nk,0:lmax))
        call gauleg(nk,0.0_dpreal,5.0_dpreal,kk,kw)
        call bassel()
        allocate_k=.true.
      end if

      G0=0.0_dpreal
      do ir=1,nr
        do irp=1,nr
          do ik=1,nk
            G0k=1.0_dpreal/(ecmx-hbarc**2 * kk(ik)**2 /2.0_dpreal/mu)
            integral=NA * kk(ik)**2 * fl(ir,ik,lx) * fl(irp,ik,lx) * G0k
            G0(ir,irp)=G0(ir,irp) + integral * kw(ik)
          end do
C          if(ir==20)  write(199,*)rr(irp),G0(ir,irp)
        end do
      end do
C      stop
      end subroutine
c-----------------------------------------------------------------------
      subroutine  bassel()
       use coulfunc
       use channels
       use precision
       use mesh
       implicit none
       real*8,dimension(0:lmax) :: GC,FCP,GCP,FC
       real*8 :: rho
       integer :: ir,ik,ifail
       ifail=0
       do ik=1,nk
         do ir=1,nr
            rho=kk(ik)*rr(ir)
            call coul90(rho,0.0_dpreal,0.0_dpreal,lmax,fc,gc,fcp,gcp,0,ifail)
            fl(ir,ik,0:lmax)=fc/rho
            if (ifail/=0) then
               write(*,*) 'coul90: ifail=',ifail; stop
            endif
         end do
! problem happened when lmax is larger than 70, when rho is less than 0.01
! the value of fc will be NAN
!         fc=>fl(1,1,0:lmax)
!         rho=kk(1)*rr(1)
!         call coul90(0.001d0,0.0d0,0.0d0,lmax,fc,gc,fcp,gcp,0,ifail)
!         write(*,*) "fl=", fl(1,1,0)
       end do

      end subroutine
c--------------------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  Green_func_nonlocalpot(ecmx,lx)
C     this subroutine is used to calculate the full green's function
C     with nonlocal potential, I would like used the matrix form to
C     perform the calcuation
C     G_x= G_c + G_c U_x G_x, where G_c is the Green function with
C     Coulomb potentia
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use constants
      use systems
      use precision
      use mesh
      use pot
      use channels
      use interpolation
      use coulfunc
      implicit none
      real*8 :: ecmx
      complex*16 ::NA
      real*8 :: kx,rho,eta,k2
      real*8 :: mux,z12,rx
      integer :: ir,ifail,lx,irp,INFO
      real*8,dimension(0:lmax) :: fc,gc,fcp,gcp
      real*8,dimension(1:nr,0:lmax) :: Re_func,Ir_func
      complex*16,dimension(1:nr,1:nr) :: G0,U,A,Gx_mat
      integer,dimension(1:nr) :: ipIV
      integer :: i


      mux=amu*(masst*massx)/(massx+masst)
      k2=2.*mux*ecmx/(hbarc**2)
      kx=sqrt(k2)
      z12=zt*zx
      eta=z12*e2*mux/hbarc/hbarc/kx

      do ir=1,nr
        rho=rr(ir)*kx
        call coul90(rho,eta,zero,lmax,fc,gc,fcp,gcp,0,ifail)
        Re_func(ir,:)=fc/rho
        Ir_func(ir,:)=fc/rho+iu*gc/rho
      end do

      NA=-2.0_dpreal*mux*iu*kx/hbarc/hbarc
      do ir=1,nr
         do irp=1,nr
           G0(ir,irp)=NA*Re_func(min(irp,ir),lx)*Ir_func(max(irp,ir),lx)
         end do
      end do

      U=0.0_dpreal

      do ir=1,nr
         rx=rr(ir)
         U(ir,ir)=FFC(rx/hcm,UxA(:,1),irmatch+1)
      end do


      call Amat(G0,U,A)


      Gx_mat=G0
      call ZGESV(nr,nr,A,nr,ipIV,Gx_mat,nr,INFO)
!          ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!  ZGESV computes the solution to a complex system of linear equations
!     A * X = B,
!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
! B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
! IPIV    (output) INTEGER array, dimension (N)
!              The pivot indices that define the permutation matrix
!              P; row i of the matrix was interchanged with row
!              IPIV(i).
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
      if (info/=0) then
         write(*,*) "problem happened when calling ZGESV"
         write(*,*) "info=",info
         stop
      end if

      i=(1+nr)/2
      write(*,*) "r=",rr(i)
      do ir=1,nr
        write(100,*)rr(ir),real(Gx_mat(ir,i)),imag(Gx_mat(ir,i))
      end do

      write(100,*) "&"






      end subroutine
c--------------------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Amat(G0,U,A)
C     this subroutine is used to calculate the A matrix used for
C     calculating the full green function with the integral equations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      use precision
      implicit none
      complex*16,dimension(1:nr,1:nr) :: G0, U, A
      complex*16 :: func
      integer ::  i, k, l
      A=0.0_dpreal

      do i=1,nr
        A(i,i)=1
        do l=1,nr
          func=0
          do k=1,nr
            func=func + rr(k)**2 * rr(l)**2 * rrw(k) * rrw(l) * G0(i,k) * U(k,l)
          end do
          A(i,l)=A(i,l)-func
        end do
      end do


      end subroutine
c--------------------------------------------------------------------------------

      end module
