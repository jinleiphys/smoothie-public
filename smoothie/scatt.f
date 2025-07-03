      module scatt
      implicit none
      complex*16,dimension(:,:),allocatable :: wf_a,wf_b
      real*8,dimension(:),allocatable :: nfc,ngc,nfcp,ngcp
      complex*16,dimension(:),allocatable :: smat_a, smat_b

      type greenf
         complex*16,allocatable,dimension(:,:) :: re,ir
      end type

      type(greenf) :: Gx
      private :: nfc,ngc,nfcp,ngcp,sch,matching,chanaout
      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chana()
c     informations of channel a
c the coupling coffcient for this channel should be | (l(jxjt)sxt) jxt, (lam 0) lam ; J
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       use channels
       use systems,only:zp,massp,zt,masst,elab
       use constants,only:amu,hbarc,e2,zero
       use precision
       use pot
       use coulfunc
       use lagrange_mesh_single_channel
       use input, only: surfacesw
       implicit none

       integer :: l,s,j,nch
       integer :: ifail ! for subroutine coul90
       real*8 :: mu,k,rho
       real*8 :: ecm
       real*8 :: eta !Sommerfeld parameter
       complex*16 :: sl,nl !  s-matrix and normalization parameter
       real*8,dimension(0:lmax) :: cph  !Coulomb phase-shift
       integer :: r0
       real*8 :: ls
       complex*16,dimension(1:5) :: wfmatch
       integer :: ir
       complex*16 :: hc,hc1  !H(+),H(-)

       if(allocated(nfc)) deallocate(nfc)
       if(allocated(ngc)) deallocate(ngc)
       if(allocated(ngcp)) deallocate(ngcp)
       if(allocated(nfcp)) deallocate(nfcp)
       if(allocated(wf_a)) deallocate(wf_a)
       if(allocated(UaA)) deallocate(UaA)
       if(allocated(smat_a)) deallocate(smat_a)


       allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
       allocate(wf_a(0:irmatch,1:inspect%nchmax))
       allocate(UaA(0:irmatch,1:inspect%nchmax))
       allocate(smat_a(1:inspect%nchmax))

       smat_a = 0.0_dpreal 
       ecm=elab*masst/(massp+masst)
       mu=amu*(masst*massp)/(massp+masst)
       k=sqrt(2.*mu*ecm/(hbarc**2))
       rho=(irmatch-2)*hcm*k
       eta=zp*zt*e2*mu/hbarc/hbarc/k


      !  write(*,*) "ka=",k
       if (k*hcm>0.2)  then
       write(*,*) 'warning!please decrease the value of hcm,',
     +     'it should be smaller than ', 0.2/k
c          stop
       end if


       call coul90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
       if (ifail/=0) then
       write(*,*) 'coul90: ifail=',ifail; stop
       endif

       call coulph(eta,cph,lmax)
       write(*,11)
       write(*,12)
       write(*,11)
11     format('***************************************************************************')
12     format('*                      ELASTIC S-MATRIX ELEMENTS                          *')

       ! insert test for lagrange_mesh_single_channel
c       call initial_lagrange_func(rmax-hcm*2.)
c       call T_and_Bloch(mu)
       ! can be removed

       write(3,13)
13     FORMAT(10x,'nch', 10x, 'real(sl)', 16x, 'aimag(sl)', 16x, 'abs(sl)')
       do nch=1,inspect%nchmax
          if (nch==1) printpot=.true.
          if (nch/=1) printpot=.false.
          l=inspect%lam(nch)
          s=jt
          j=inspect%j(nch)
          ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
          call potr('a',1,zp*zt,ls,0)
          UaA(0:irmatch,nch)=v

          r0=2*l
          call sch(r0,mu,ecm,UaA(0:irmatch,nch),l,wf_a(:,nch))
          wfmatch(1:5)=wf_a(irmatch-4:irmatch,nch)
          call matching(l,k,wfmatch,sl,nl)
          smat_a(nch) = sl 


          write(2,*) l, (atan2(aimag(sl), real(sl)) * 0.5_dpreal) * (180.0_dpreal / pi)

          wf_a(:,nch)=wf_a(:,nch)*nl

          if(surfacesw) then
            do ir=1, irmatch
              call coul90(ir*hcm*k,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
              hc=cmplx(ngc(l),nfc(l),kind=8)
              hc1=cmplx(ngc(l),-nfc(l),kind=8)
              wf_a(ir,nch)=0.5*iu*(hc1-sl*hc)
            end do
          end if


          wf_a(:,nch)=wf_a(:,nch)*exp(iu*cph(l))
          call chanaout(nch,wf_a(0:irmatch,nch),sl)

          ! insert for testing lagrange_mesh_single_channel
C         call R_matrix(l,mu,ecm,UaA(0:irmatch,nch),cph(l),ngc(l),ngcp(l),nfc(l),nfcp(l))
          ! can be removed


       end do


       call flush(3)
       write(*,11)
       write(*,*)

      end subroutine chana
c-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sch(r0,mu,ecm,vpot,l,rwfl)
c     mu! reduce mass
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh, only: irmatch,hcm
       use constants, only:hbarc,e2
       use precision
       implicit none
       integer,intent(in) :: r0
       real*8,intent(in) :: ecm ! energy
       complex*16,dimension(0:irmatch),intent(in) :: vpot
       integer,intent(in) :: l
       real*8,intent(in) :: mu  ! reduced mass
       integer :: ir
       real*8 :: r
       complex*16,dimension(0:irmatch),intent(out) :: rwfl   ! partial radial wave function
       complex*16,dimension(0:irmatch) :: kl
       complex*16,dimension(0:irmatch) :: Tx
       complex*16,dimension(0:irmatch) :: Wx


        kl=0.0d0
        rwfl=0.0d0
        Tx=0.0d0


c      Numerov method to solve the differential radial equation
       rwfl(r0)=0      ! boundary condition


       ir=r0+1; r=ir*hcm
       rwfl(ir)=hcm**(l+1)  ! arbitrary value
       if (l>100) rwfl(ir)=0.0000001_dpreal


       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)

       rwfl(r0+2)=2.*rwfl(r0+1)-hcm**2*kl(r0+1)*rwfl(r0+1)
       ir=r0+2; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)


       do ir=r0+2 ,irmatch-1

        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-2.*mu*Vpot(ir+1)/hbarc**2
        Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
        Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)-Wx(ir-1)
        rwfl(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))

       end do


      end subroutine sch

c----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matching(l,k,wf,sl,nl) !wf has dimension (5)
c     nl*wf=0.5*i*(H(-)-sl*H(+))
c     nl*wfp=0.5*i*k*(H'(-)-sl*H'(+))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use precision, only:pi,iu
       use mesh,only:hcm
       implicit none

       integer,intent(in) :: l
       real*8,intent(in) :: k
       complex*16,intent(in),dimension(1:5) :: wf ! wavefunction at rmatch
       complex*16 :: wfp !derivative of wf
       complex*16,intent(out) :: nl !Normalization parameter
       complex*16 :: hc,hc1  !H(+),H(-)
       complex*16 ::hcp,hcp1 ! derivatives of H(+),H(-)
       complex*16,intent(out) :: sl ! S-matrix

       hc=cmplx(ngc(l),nfc(l),kind=8)
       hc1=cmplx(ngc(l),-nfc(l),kind=8)
       hcp=cmplx(ngcp(l),nfcp(l),kind=8)
       hcp1=cmplx(ngcp(l),-nfcp(l),kind=8)

       wfp=(-wf(5)+8.*wf(4)-8.*wf(2)+wf(1))/12./hcm
       nl=(hc1*hcp*iu*k-hc*hcp1*iu*k)/(2.*(hcp*wf(3)*k-hc*wfp))
       sl=(hc1*wfp-hcp1*wf(3)*k)/(hc*wfp-hcp*wf(3)*k)


      end subroutine matching
c----------------------------------------------------------------------



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chanaout(nch,rwfl,sl)
c     output subroutine for channel a
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       implicit none
       integer,intent(in) :: nch
       integer :: ir
       complex*16,dimension(0:irmatch) :: rwfl
       real*8,dimension(0:irmatch) :: rwflr,rwfli ! real part and imaginary part
       complex*16,intent(in) :: sl

       rwflr=real(rwfl)
       rwfli=aimag(rwfl)
c-----------------------------------------------------------------------
c*** print elastic S-matrix elements.
C       write(3,*)"elastic S-matrix elements for channel a"
       write(3,*)nch,real(sl),aimag(sl),abs(sl)
C100    format(F13.9,2X,F13.9,2X,I3,2X)
c-----------------------------------------------------------------------

       write(*,101) nch,real(sl),aimag(sl)
101    format('nch=',I3,2X,'S-matrix = ('
     &           ,F13.9,',',F13.9,')')
c-----------------------------------------------------------------------
       write(4,102)nch
102    format('@nch=',I3)

       do ir=0,irmatch
          write (4,*) hcm*ir, rwflr(ir),rwfli(ir)
       end do
c----------------------------------------------------------------------
       write (4,*) "& "
      end subroutine chanaout




c-----------------------------------------------------------------------
c
c     a + A -> b + x + A ->     b + B -> b+anything
c     =====        =====        =====
c       a            x            b
c     informations of channel b
c     in this channel we assumed we ignore the spin-orbit interaction
      subroutine chanb(counter,ecm,chi)
c-----------------------------------------------------------------------
       use channels
       use mesh
       use systems,only:zb,massb,zt,masst,massx,zx
       use constants,only:amu,hbarc,e2,zero,convert
       use precision
       use pot
       use coulfunc
       use input, only:surfacesw
       implicit none

       integer :: l,nch
       real*8 :: ls,j,s
       integer :: ifail ! for sunroutine coul90
       real*8 :: mu,k,rho
       integer,intent(in) :: counter
       real*8,intent(in) :: ecm
       real*8 :: eta !Sommerfeld parameter
       complex*16 :: sl,nl !  s-matrix and normalization parameter
       integer :: r0
       integer :: ir
       complex*16,dimension(0:irmatch,1:outspect%nchmax) :: chi
       real*8,dimension(0:lmax) :: cph
       complex*16 :: hc,hc1  !H(+),H(-)

       if (allocated(nfc)) deallocate(nfc)
       if (allocated(ngc)) deallocate(ngc)
       if (allocated(nfcp)) deallocate(nfcp)
       if (allocated(ngcp)) deallocate(ngcp)
       if (allocated(UbB)) deallocate(UbB)
       if (allocated(smat_b)) deallocate(smat_b)

       allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
       allocate(UbB(0:irmatch,1:outspect%nchmax))
       allocate(smat_b(1:outspect%nchmax))

       smat_b=0.0_dpreal 
       mu=amu*((masst+massx)*massb)/(massb+masst+massx)
       k=sqrt(2.*mu*ecm/(hbarc**2.))
       rho=hcm*(irmatch-2)*k
       eta=zb*(zt+zx)*e2*mu/hbarc/hbarc/k

       call coulph(eta,cph,lmax)
       call coul90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
          if (ifail/=0) then
          write(*,*) 'coul90: ifail=',ifail; stop
       endif


       do nch=1,outspect%nchmax
         if (nch==1) printpot=.true.
         if (nch/=1) printpot=.false.
         l=outspect%lam(nch)
         s=jb
         j=outspect%j(nch)
         ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))

         ls=0.0_dpreal
         call potr('b',counter,zb*(zt+zx),ls,0)
         UbB(0:irmatch,nch)=V(0:irmatch)

         r0=2*l
         call sch(r0,mu,ecm,UbB(0:irmatch,nch),l,chi(0:irmatch,nch))
       	 call matching(l,k,chi(irmatch-4:irmatch,nch),sl,nl)
       	 
       	 smat_b(nch) = sl

       	 chi(0:irmatch,nch)=chi(0:irmatch,nch)*nl

       	 if(surfacesw) then
            do ir=1, irmatch
              call coul90(ir*hcm*k,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
              hc=cmplx(ngc(l),nfc(l),kind=8)
              hc1=cmplx(ngc(l),-nfc(l),kind=8)
              chi(ir,nch)=0.5*iu*(hc1-sl*hc)
            end do
          end if


       	 chi(0:irmatch,nch)=chi(0:irmatch,nch)*exp(iu*cph(l)) !  this phase is needed, details see the notes for coulomb phase



       	 if (counter==1) then
       	   write(100+counter,*)'& nch=',nch
       	   do ir=0,irmatch
       		  write (100+counter,*) hcm*ir, real(chi(ir,nch)),aimag(chi(ir,nch))
       	   end do
       	 end if
       end do
      end subroutine
c-----------------------------------------------------------------------
      subroutine screenchanb(counter,ecm)
        use channels
        use mesh
        use constants
        use precision
        use interpolation
        use systems,only:zb,massb,zt,masst,massx,zx
        implicit none

        integer :: nch
        integer :: counter
        real*8 :: ecm
        real*8 :: r
        integer :: ir
      !   real*8 :: Rs
        complex*16,dimension(0:irmatch,1:outspect%nchmax) :: chi

        if (allocated(wf_b)) deallocate(wf_b)

        allocate(wf_b(1:nr,1:outspect%nchmax))


        call chanb(counter,ecm,chi)


        do nch=1,outspect%nchmax

          do ir=1,nr
            r=rr(ir)
            wf_b(ir,nch)=FFC(r/hcm,chi(0:irmatch,nch),irmatch+1)
          end do



C          Rs=120.0_dpreal
C          n=4
C          do ir=1,nr
C            r=rr(ir)
C            wf_b(ir,nch)=wf_b(ir,nch)*exp(-(r/Rs)**n)
C          end do

          if (counter==1) then
             write(150+counter,*) "&nch=",nch
             do ir=1,nr
                write(150+counter,*) rr(ir),real(wf_b(ir,nch)), aimag(wf_b(ir,nch))
             end do
          end if

        end do
      end subroutine
c--------------------------------------------------------------------------------
      subroutine screenchanbrbx(counter,ecm)
        use channels
        use mesh
        use constants
        use precision
        use interpolation
        use systems,only:zb,massb,zt,masst,massx,zx
        implicit none

        integer :: nch
        integer :: counter
        real*8 :: ecm
      !   real*8 :: r
        integer :: ir
      !   real*8 :: Rs
        complex*16,dimension(0:irmatch,1:outspect%nchmax) :: chi

        if (allocated(wf_b)) deallocate(wf_b)

        allocate(wf_b(0:irmatch,1:outspect%nchmax))


        call chanb(counter,ecm,chi)

        wf_b=chi





          do nch=1,outspect%nchmax
          if (counter==1) then
             write(150+counter,*) "&nch=",nch
             do ir=0,irmatch
                write(150+counter,*) ir*hcm,real(wf_b(ir,nch)), aimag(wf_b(ir,nch))
             end do
          end if

          end do


      end subroutine




c-----------------------------------------------------------------------
c     a + A -> b + x + A -> x + b + A -> b+anything
c     =====        =====        =====
c       a            x            b
c     this subroutine is used to calculate the green function
c                        of x-channel
      subroutine greenfunc(counter,ecm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use channels
       use mesh
       use systems,only:zx,massx,zt,masst
       use constants,only:amu,hbarc,e2,zero,convert
       use pot
       use precision
       use coulfunc
       use interpolation
C       use bound
       implicit none
       integer :: l,nch
       real*8 :: s,j,ls
       integer :: ifail ! for subroutine coul90
       integer :: ir
       real*8 :: mu,k,rho,r
       real*8 :: eta !Sommerfeld parameter
       integer,intent(in) :: counter
       real*8,intent(in) :: ecm
       complex*16 :: sl,nl !  s-matrix and normalization parameter
       integer :: r0
       complex*16  :: hlp,flp
       complex*16 :: NW
       complex*16,dimension(0:irmatch) :: flx,hlx
       real*8,dimension(0:lmax) :: cph
      !  integer :: i,ich_pauli

       if (allocated(nfc)) deallocate(nfc)
       if (allocated(ngc)) deallocate(ngc)
       if (allocated(nfcp)) deallocate(nfcp)
       if (allocated(ngcp)) deallocate(ngcp)
       if (allocated(Gx%re)) deallocate(Gx%re)
       if (allocated(Gx%ir)) deallocate(Gx%ir)
       if (allocated(UxA)) deallocate(UxA)


       allocate(Gx%re(1:nr,1:out2b%nchmax))
       allocate(Gx%ir(1:nr,1:out2b%nchmax))
       allocate(UxA(0:irmatch,1:out2b%nchmax))
       allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
       UxA=0.0_dpreal
       flx=0.0_dpreal
       hlx=0.0_dpreal



       if(ecm>0) then

        mu=amu*(masst*massx)/(massx+masst)
        k=sqrt(2.*mu*ecm/(hbarc**2.))
        rho=hcm*(irmatch-2)*k
        eta=zx*zt*e2*mu/hbarc/hbarc/k

        call coulph(eta,cph,lmax)
        call coul90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
        if (ifail/=0) then
           write(*,*) 'coul90: ifail=',ifail; stop
        endif


       do nch=1,out2b%nchmax

         l=out2b%l(nch)
         s=out2b%s(nch)
         j=out2b%j(nch)
         ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))

          if (nch==1) printpot=.true.
          if (nch/=1) printpot=.false.
         call potr('x',counter,zx*zt,ls,0)
         UxA(0:irmatch,nch)=v(0:irmatch)
         
         if (adde) then 
         call potr('e',counter,0.0_dpreal,ls,0)
         UxA(0:irmatch,nch)=UxA(0:irmatch,nch)+v(0:irmatch)
         if (nch==1) then 
         do ir=1, irmatch 
             write(333,*) ir*hcm, real(UxA(ir,nch)),aimag(UxA(ir,nch))
         end do 
         write(333,*)"&"
         end if 
         end if 

         r0=2*l
       	 call sch(r0,mu,ecm,UxA(0:irmatch,nch),l,flx)
         call matching(l,k,flx(irmatch-4:irmatch),sl,nl)
       	 flx=flx*nl

         call hlxkr(counter,l,UxA(0:irmatch,nch),mu,ecm,cph(l),hlx)

         if (counter==1) then
           write(200+counter,*)'& nch=',nch
           do ir=0,irmatch
              write (200+counter,*) hcm*ir, real(flx(ir)), aimag(flx(ir))
           end do
         end if

         do ir=1,nr
           r=rr(ir)
           Gx%re(ir,nch)=FFC(r/hcm,flx,irmatch+1)
           Gx%ir(ir,nch)=FFC(r/hcm,hlx,irmatch+1)
         end do

       end do



c----------------------
       else
         mu=amu*(masst*massx)/(massx+masst)
         k=sqrt(2.*mu*abs(ecm)/(hbarc**2.))
       do nch=1,out2b%nchmax

          l=out2b%l(nch)
          l=out2b%l(nch)
          s=out2b%s(nch)
          j=out2b%j(nch)
          ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
          if (nch==1) printpot=.true.
          if (nch/=1) printpot=.false.
          call potr('x',counter,zx*zt,ls,0)
          UxA(0:irmatch,nch)=v(0:irmatch)

       !Pauli projection block
C        ich_pauli=-99
C        do i=1,10
C          if(nch == a2bindex(i)) ich_pauli=i
C        end do
C        if(ich_pauli>0) then
C          do ir=1,irmatch
C            UxA(ir,nch) = UxA(ir,nch) + tau*
!!!!!!!!!!!!!!!!!!need to finish with nonlocal potential !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!R-matrix method method is not ready for the negative energies !!!!!
C          end do
C        end if
       !finish Pauli projection block




          r0=2*l
          call sch(r0,mu,ecm,UxA(0:irmatch,nch),l,flx)
          call hlxkr(counter,l,UxA(0:irmatch,nch),mu,ecm,cph(l),hlx)

          hlp=(-hlx(irmatch)+8.*hlx(irmatch-1)-8.*hlx(irmatch-3)+hlx(irmatch-4))/12./hcm
          flp=(-flx(irmatch)+8.*flx(irmatch-1)-8.*flx(irmatch-3)+flx(irmatch-4))/12./hcm

          NW=flx(irmatch-2)*hlp-flp*hlx(irmatch-2)
          flx=-flx*k/NW

          if (counter==1) then
          write(200+counter,*)'& nch=',nch
       	  do ir=0,irmatch
             write (200+counter,*) hcm*ir, real(flx(ir)), aimag(flx(ir))
          end do
          end if

          do ir=1,nr
            r=rr(ir)
            Gx%re(ir,nch)=FFC(r/hcm,flx,irmatch+1)
            Gx%ir(ir,nch)=FFC(r/hcm,hlx,irmatch+1)
          end do

       end do

       end if



      end subroutine
c-----------------------------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hlxkr(counter,l,vpot,mu,ecm,cphase,hlx)
c     This subroutine is used to calculate the hlx(r)
c              for each energy point
c           hlx(l,irmatch)=ngc(l)+i*nfc(l)
c           use Numerov method from outside to inside
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       use channels,only:lmax
       use systems,only:zt,zx
       use constants, only:hbarc,e2,zero
       use precision
       use whittaker
       use coulfunc
       implicit none
       integer :: counter   ! outgoing channel index
       integer :: l
       real*8 :: ecm ! energy
       real*8 :: mu  ! reduced mass
       complex*16,dimension(0:irmatch) :: hlx
       complex*16,dimension(0:irmatch) :: Vpot
       real*8 :: cphase
       real*8 :: k         ! k=sqrt(2*mu*ecm/(hbarc**2))
       real*8 :: r,eta,rho
       real*8 :: z12 ! For coulomb potential and eta zt*zx
       real*8 :: const
       integer :: ir ,ifail
       complex*16,dimension(1:irmatch) :: kl
       real*8,dimension(0:lmax) :: WK,WKD
       real*8,dimension(0:lmax) :: nfc1,ngc1,nfcp1,ngcp1

       k=sqrt(2.*mu*abs(ecm)/(hbarc**2))
       z12=zt*zx
       eta=z12*e2*mu/hbarc/hbarc/k
       hlx=0.0d0

c----boundary condition
       if (ecm>0) then
          rho=hcm*(irmatch)*k
          call coul90(rho,eta,zero,lmax,nfc1,ngc1,nfcp1,ngcp1,0,ifail)
          hlx(irmatch)=cmplx(ngc1(l),nfc1(l),kind=8)

          rho=hcm*(irmatch-1)*k
          call coul90(rho,eta,zero,lmax,nfc1,ngc1,nfcp1,ngcp1,0,ifail)
          hlx(irmatch-1)=cmplx(ngc1(l),nfc1(l),kind=8)

       else

          ifail=0
          rho=hcm*irmatch
          call WHIT(eta,rho,k,ecm,lmax,WK,WKD,ifail)
          hlx(irmatch)=WK(l)

          rho=hcm*irmatch-hcm
          call WHIT(eta,rho,k,ecm,lmax,WK,WKD,ifail)
          hlx(irmatch-1)=WK(l)

       end if


c      Numerov method to solve the differential radial equation

        ir=irmatch;r=ir*hcm
       	kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2

       	ir=irmatch-1;r=ir*hcm
       	kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2

       	const=hcm**2/12.


       	do ir=irmatch-1,2,-1
           kl(ir-1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir-1.)*hcm)**2-2.*mu*Vpot(ir-1)/hbarc**2
           hlx(ir-1)=((2.-10.*const*kl(ir))*hlx(ir)-(1.+const*kl(ir+1))*hlx(ir+1))/(1.+const*kl(ir-1))
       	end do

        if (counter==1) then
        write(250+counter,*)'& l=',l
       	do ir=1,irmatch
       	  write (250+counter,*) hcm*ir, real(hlx(ir)),aimag(hlx(ir))
       	end do
       	end if


      end subroutine hlxkr
c-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bbin(counter)
c     bin method for b waves
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use constants
      use systems
      use mesh
      use interpolation
      implicit none
      integer,intent(in) :: counter
      integer :: ir,nch,i
      real*8 :: kb,kbmax,kbmin,k
      real*8 :: eta !Sommerfeld parameter
      real*8 :: ecmb,mub,ecm,deltak,ebinmax,ebinmin
      real*8 :: norm,r
      complex*16 :: nl
      complex*16,dimension(0:irmatch,1:outspect%nchmax) :: chi
c***

      if (allocated(wf_b)) deallocate(wf_b)
      allocate(wf_b(0:irmatch,1:outspect%nchmax))
      wf_b=0.0_dpreal
c***
c***

      ecmb=necmb(counter)
      mub=amu*((masst+massx)*massb)/(massb+masst+massx)
      kb=sqrt(2*mub*ecmb/(hbarc**2))
      eta=zb*(zt+zx)*e2*mub/hbarc/hbarc/kb
      kbmax=kb+nwbin(counter)/2.
      kbmin=kb-nwbin(counter)/2.
      ebinmax=(hbarc*kbmax)**2/2./mub
      ebinmin=(hbarc*kbmin)**2/2./mub
      if(kbmin<0) then
      write(*,10)
      stop
      end if
10    format('!!error!! kbmin<0,,please decrease the width of b-bin')


      write(*,*)'using the b-bin method'
      write(*,12)ebinmax,ebinmin,kbmax,kbmin,nbin(counter)
12    format('#bin',3X,'Ebmax=',f7.3,3X,'Ebmin=',f7.3,3x,'kbmax=',f7.3,
     &  3x,'kbmin=',f7.3,3x,'Number of bins=',I3)
      deltak=nwbin(counter)/nbin(counter)
      do i=1,nbin(counter)
            k=kbmin+deltak*i-deltak/2.
            ecm=(hbarc*k)**2/2/mub
            call chanb(counter,ecm,chi)
            wf_b=wf_b+deltak*chi
      end do
       wf_b=wf_b*sqrt(2/pi/nwbin(counter))
c*********normalization check
      write(*,13)
13    format('#bin',3x,'nch',6x,'Norm')

      do nch=1,outspect%nchmax
         norm=0.0d0
         do ir=1,nr
            r=rr(ir)
            norm=norm+abs(FFC(r/hcm,wf_b(0:irmatch,nch),irmatch+1))**2*rrw(ir)

         end do

         write(*,14) nch, norm
      end do
14    format(5x,I3,5x,f7.4)
c**********
c*****

       do nch=1,outspect%nchmax
        nl=sqrt(pi/2/nwbin(counter))
        wf_b(:,nch)= wf_b(:,nch)*nl

        if (counter==1) then
       	write(150+counter,*)'& nch=',nch

       	do ir=0, irmatch
       		write (150+counter,*) hcm*ir, real(wf_b(ir,nch)),
     &            aimag(wf_b(ir,nch))

       	end do

       	end if
       end do


      end subroutine


c-----------------------------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bbin_rb(counter)
c     bin method for b waves
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use constants
      use systems
      use mesh
      use interpolation
      implicit none
      integer,intent(in) :: counter
      integer :: ir,nch,i
      real*8 :: kb,kbmax,kbmin,k
      real*8 :: eta !Sommerfeld parameter
      real*8 :: ecmb,mub,ecm,deltak,ebinmax,ebinmin
      real*8 :: norm,r
      complex*16 :: nl
      complex*16,dimension(0:irmatch,1:outspect%nchmax) :: chi,chi_bin
c***

      if (allocated(wf_b)) deallocate(wf_b)
      allocate(wf_b(1:nr,1:outspect%nchmax))
      wf_b=0.0_dpreal
      chi_bin=0.0_dpreal
c***
c***

      ecmb=necmb(counter)
      mub=amu*((masst+massx)*massb)/(massb+masst+massx)
      kb=sqrt(2*mub*ecmb/(hbarc**2))
      eta=zb*(zt+zx)*e2*mub/hbarc/hbarc/kb
      kbmax=kb+nwbin(counter)/2.
      kbmin=kb-nwbin(counter)/2.
      ebinmax=(hbarc*kbmax)**2/2./mub
      ebinmin=(hbarc*kbmin)**2/2./mub
      if(kbmin<0) then
      write(*,10)
      stop
      end if
10    format('!!error!! kbmin<0,,please decrease the width of b-bin')


      write(*,*)'using the b-bin method'
      write(*,12)ebinmax,ebinmin,kbmax,kbmin,nbin(counter)
12    format('#bin',3X,'Ebmax=',f7.3,3X,'Ebmin=',f7.3,3x,'kbmax=',f7.3,
     &  3x,'kbmin=',f7.3,3x,'Number of bins=',I3)
      deltak=nwbin(counter)/nbin(counter)
      do i=1,nbin(counter)
            k=kbmin+deltak*i-deltak/2.
            ecm=(hbarc*k)**2/2/mub
            call chanb(counter,ecm,chi)
            chi_bin=chi_bin+deltak*chi
      end do
       chi_bin=chi_bin*sqrt(2/pi/nwbin(counter))
c*********normalization check
      write(*,13)
13    format('#bin',3x,'nch',6x,'Norm')

      do nch=1,outspect%nchmax
         norm=0.0d0
         do ir=1,nr
            r=rr(ir)
            norm=norm+abs(FFC(r/hcm,chi_bin(0:irmatch,nch),irmatch+1))**2*rrw(ir)

         end do

         write(*,14) nch, norm
      end do
14    format(5x,I3,5x,f7.4)
c**********
c*****

       do nch=1,outspect%nchmax
        nl=sqrt(pi/2/nwbin(counter))
        chi_bin(:,nch)= chi_bin(:,nch)*nl

      
       if (counter==1) 	write(150+counter,*)'& nch=',nch

       	do ir=1, nr
       	    r=rr(ir)
       	    wf_b(ir, nch) = FFC(r/hcm,chi_bin(0:irmatch,nch),irmatch+1)
       if (counter==1)	write (150+counter,*) r, real(wf_b(ir,nch)),
     &            aimag(wf_b(ir,nch))

       	end do


       end do


      end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chana_readin()
c     informations of channel a
c the coupling coffcient for this channel should be | (l(jxjt)sxt) jxt, (lam 0) lam ; J
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch
       use channels
       use pot
       use precision

       implicit none

       integer :: l,s,j,nch,alphain,ir
       real*8 :: ls
       complex*16 :: sl

       sl=0.0_dpreal


       if(allocated(wf_a)) deallocate(wf_a)
       if(allocated(UaA)) deallocate(UaA)


       allocate(wf_a(0:irmatch,1:inspect%nchmax))
       allocate(UaA(0:irmatch,1:inspect%nchmax))


        open (unit=553,file='cdccwf_chi.dat')
       do nch=1,inspect%nchmax

          l=inspect%lam(nch)
          s=jt
          j=inspect%j(nch)
          ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
          if (nch==1) printpot=.true.
          if (nch/=1) printpot=.false.
          call potr('a',1,zp*zt,ls,0)
          UaA(0:irmatch,nch)=v

          read(553,*) alphain
          if(nch/=alphain) write(*,*) "error!!!!!!!!!!"

          write(*,*) "read wave function, alphain=",nch
           do ir=0, irmatch

           read(553,*) wf_a(ir,nch)



           end do

          call chanaout(nch,wf_a(0:irmatch,nch),sl)

       end do

       write(*,150)
150   format('********************************************************')

      end subroutine
c-----------------------------------------------------------------------



      end module
