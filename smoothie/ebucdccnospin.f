      module ebucdccnospin
      use precision
      use angularmesh
      implicit none

        complex*16,allocatable,dimension(:,:,:) :: psicdcc

        type cdccwf_seprate
        complex*16,allocatable,dimension(:,:,:) :: phi,chi
        end type

        type(cdccwf_seprate) :: cdccwf

        real*8,dimension(:,:),target,allocatable :: Y_out_lx
        real*8,target,dimension(:,:,:,:),allocatable :: Y_in_la
        real*8,target,dimension(:,:,:,:),allocatable ::Y_out_lb
        real*8 :: mua,mub,mux,ka,kb,kx



        complex*16,allocatable,dimension(:,:,:,:) :: fun_in
        complex*16,allocatable,dimension(:,:,:,:) :: pot_in

        complex*16,dimension(:),allocatable :: smat_abar



      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine EBU_cdcc_Zero_Spin()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use cleb_coefficient
        use mesh
        use constants
        use systems
        use channels
        use bound
        use scatt
        use pot
        use gauss
        use spharm
        use input, only:prior,icf,written,HM
        use coulfunc
        implicit none
        real*8 :: ecm,ecmb,ecmx
        integer :: counter
        integer :: alpha_xA,nch
        real*8,dimension(1:out2b%nchmax)  :: xsec_alpha_xA
        integer :: alphabar,alpha_out
        integer :: nth,ix,ir
        real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kplb
        real*8 :: theta, cth,xsec_sum,xsec,t1,t2
        real*8 :: omp_get_wtime
        integer :: jp_index,abarmax
        real*8,dimension(0:nthmax) :: dsdedw_sum,dsdedw
        real*8 :: eta
        real*8,dimension(0:lmax) :: cph
        real*8,dimension(1:inspect%nchmax) :: xsecla,xsecla_lx
        integer :: lxmax1


        allocate(UbA(0:irmatch))



        mua=(massp * masst) * amu / (massp + masst)
        mub=(massb * (massx+masst)) * amu / (massb + massx+masst)
        mux=(massx * masst) * amu / (massx + masst)
        ecm=elab * masst / (masst + massp)
        ka=sqrt(2*mua*ecm/(hbarc**2))

        eta=zp*zt*e2*mua/hbarc/hbarc/ka
        call coulph(eta,cph,lmax)

        call ang_mesh()
        call readincdccwf_seprate()
c---- calculate the potential between x and b
        allocate(vbx(0:irmatch,1))
        vbx=0.0_dpreal
        call potr('p',1,zx*zb,0.0_dpreal,0)
        vbx(0:irmatch,1)=v(0:irmatch)
c*** need to fixed

c calculate the icf fusion potential
C        if (icf) then
C           allocate(Wfus(0:irmatch))
C           call icf_fus_pot(Wfus)
C        end if

c!!!


        write(*,90)elab*masst/(masst+massp)+qval
90     format("Threshold at Eb=",f9.4)

         lxmax1=min(lmax,lxmax)
        if(allocated(Y_in_la)) deallocate(Y_in_la)
        if(allocated(Y_out_lb)) deallocate(Y_out_lb)
        if(allocated(Y_out_lx)) deallocate(Y_out_lx)
        allocate(Y_in_la(1:nx,1:nr,1:nr,1:lmax_cdcc*lmax_cdcc+2*lmax_cdcc+1))
        allocate(Y_out_lb(1:nx,1:nr,1:nr,1:lmax*lmax+2*lmax+1))
        allocate(Y_out_lx(1:nx,1:lxmax1*lxmax1+2*lxmax1+1))
        write(*,91)
91      format("Now initial the value of Y function ")
        t1=omp_get_wtime()
!$OMP PARALLEL default(shared) private(ix)
!$OMP DO schedule(dynamic)
        do ix=1,nx
          call initial_Y_out_rbx(ix)
          call initial_Y_in_rbx(ix)
        end do
!$OMP END DO
!$OMP end PARALLEL
       t2=omp_get_wtime()
       write(*,101) t2-t1
       write(*,120)

C this part to initial kplb
!$OMP PARALLEL default(shared) private(nth,theta,cth)
!$OMP DO schedule(dynamic)
      do nth=0,nthmax
         theta=thmin+nth*thinc
         cth=cos(theta*pi/180.)
         call PLM(cth,lmax,lmax,lmax+1,kplb(0:lmax,0:lmax,nth))
      end do
!$OMP END DO
!$OMP end PARALLEL





        do counter=1,neb
          ecmb=necmb(counter)
          ecmx=ecm+qval-ecmb
          kb=sqrt(2*mub*ecmb/(hbarc**2))
          kx=sqrt(2*mux*abs(ecmx)/(hbarc**2))


          write(*,92) ecmb
92        format('Now calculating the x-sec at Eb=',F7.3)


C          call screenchanbrbx(counter,ecmb)
          call bbin(counter)
          call greenfunc(counter,ecmx)

          call potr('t',counter,zb*zt,0.0_dpreal,0)
          UbA=V

C          call boundxA()
C          call OMP_SET_NUM_THREADS(1)

          write(*,*) "HM=",HM

          if (HM) write(*,93)
          if (.not. HM) then
           write(*,94)
C          if (nx<100) then
C             write(*,95)
C             stop
C          end if

         end if
 93      format("Using CDCC in HM form")
 94      format("Using CDCC in IVA form")
C95      format("error! in the post form nx should be very large")


          call initial_func_pot_in_cdcc()
          write(*,96)
96        format( "calculating cross section alpha_xA dependence ",/)
          t1=omp_get_wtime()

          xsec_sum=0.0_dpreal
          xsec_alpha_xA=0.0_dpreal
          dsdedw_sum=0.0_dpreal
          xsecla=0.0_dpreal
          do alpha_xA=1, out2b%nchmax
            call sigma_J_pi_cdcc(alpha_xA,xsec_alpha_xA(alpha_xA),kplb,cph,dsdedw,xsecla_lx)
              xsec_sum=xsec_sum+xsec_alpha_xA(alpha_xA)
              dsdedw_sum=dsdedw_sum+dsdedw
              xsecla=xsecla+xsecla_lx
          end do  ! jp_index

         t2=omp_get_wtime()
         write(*,101) t2-t1
         write(*,120)

         write(22,*)necmb(counter), xsec_sum


          do alpha_xA=1,out2b%nchmax
            write(*,10) alpha_xA,xsec_alpha_xA(alpha_xA)
            write(20,*) alpha_xA,xsec_alpha_xA(alpha_xA)
          end do
            write(20,*)"&"


            do nch=1,inspect%nchmax
               write(23,*)inspect%lam(nch),xsecla(nch)
            end do
            write(23,*)"&"


          write(*,100)
          do nth=0,nthmax
             theta=thmin+nth*thinc
             write(*,20)theta,dsdedw_sum(nth)
             write(21,*)theta,dsdedw_sum(nth)
          end do

          write(21,*)"&"



10    format('alpha_xA = ',I3,5x,'X-sec = ',2F9.4)
20    format('theta =',F7.2,5x,'X-sec = ',2F9.4)


        end do
100   format("The double differential cross sections",/)
101   format(20x,'(CPU  time =',F12.2,2x,'seconds)')
120   format('--------------------------------------------------------')
       end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sigma_J_pi_cdcc(alpha_xA,xsec_sum,kplb,cph,dsdedw,xsec_la)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      implicit none

       integer :: abarmax,alphabar,alpha_xA
       integer :: alpha_out
       real*8,allocatable,dimension(:) :: xsec
       real*8 :: xsec_sum
       integer :: nth,nch,nchla
       real*8,dimension(0:nthmax) :: dsdedw
       real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kplb
       real*8,dimension(0:lmax) :: cph
       real*8,dimension(1:inspect%nchmax) :: xsec_la


        dsdedw=0.0_dpreal
        xsec_sum=0.0_dpreal
        call abar_index_cdcc(alpha_xA,abarmax) ! something need to do here
        if(abarmax==0) return






        allocate(xsec(1:abarmax))
        allocate(smat_abar(1:abarmax))

!$OMP PARALLEL default(shared) private(alphabar)
!$OMP DO schedule(dynamic)
        do alphabar=1,abarmax
C           write(*,*)"there are ",abar%cdccnchmax(alphabar), "cdcc channels"
          call dsde_cdcc(alphabar,xsec(alphabar))
        end do
!$OMP END DO
!$OMP end PARALLEL

           xsec_la=0.0_dpreal
          do alphabar=1,abarmax
             xsec_sum=xsec_sum+xsec(alphabar)
             nchla=in3b%alphaspect(abar%alpha_in(alphabar))
             do nch=1,inspect%nchmax
               if (nch==nchla)  xsec_la(nch)=xsec_la(nch)+xsec(alphabar)
             end do
          end do



          write(*,20) xsec_sum



 20   format("xsec=",f19.5,/)


!$OMP PARALLEL default(shared) private(nth)
!$OMP DO schedule(dynamic)
       do nth=0,nthmax

          call dsdedw_cdcc(dsdedw(nth),kplb(:,:,nth),cph,abarmax)

       end do
!$OMP END DO
!$OMP end PARALLEL


        write(331,*)"&alpha_xA=",alpha_xA
       do nth=0, nthmax
          write(331,*)thmin+nth*thinc, dsdedw(nth)
       end do

      deallocate(smat_abar)

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dsdedw_cdcc(xsec,kplb,cph,abarmax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use precision
      implicit none
      real*8 :: xsec,xsec_alpha_xA
      integer :: mxA,mbx,lx,alpha_xA,abarmax
      real*8,dimension(0:lmax,0:lmax) :: kplb
      real*8,dimension(0:lmax) :: cph


      xsec=0.0_dpreal
      do alpha_xA=1,out2b%nchmax
        lx=out2b%l(alpha_xA)
        do mxA=-lx,lx
          do mbx=-lbx,lbx
            call psiintegration_cdcc(alpha_xA,abarmax,mbx,mxA,kplb,cph,xsec_alpha_xA)
            xsec=xsec+xsec_alpha_xA
          end do
        end do
      end do

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psiintegration_cdcc(alpha_xA,abarmax,mbx,mxA,PL,cph,xsec_alpha_xA)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use precision
      use constants
      use spharm
      use interpolation
      use pot
      use cleb_coefficient
      use mesh
      use input,only:icf
      implicit none
      real*8 :: xsec_alpha_xA
      integer :: alphabar,alpha_in,alpha_out,alpha_xA
      integer :: lx,lb,LL,la
      integer :: mxA,mbx,mb,mJ,ma
      integer :: nchla,nchlb,nchlx,nchlbx,nchLL,abarmax
      complex*16  :: Psi
      real*8 :: Yla,Ylb,YY
      integer :: irx
      real*8 :: rx
      real*8 :: N
      real*8,dimension(0:lmax,0:lmax) :: PL
      real*8 :: CGin,CGout
      real*8,dimension(0:lmax) :: cph
      N=mua * mub * kb / 4 / pi**2 / hbarc**4 /ka /(2.0_dpreal*lbx+1.0_dpreal)
      N=N * kx * mux / 8. / pi**3 / hbarc**2
      N=N*10.0_dpreal ! 1fm2=10mb


      ma=0
      Psi=0.0_dpreal
      do alphabar=1,abarmax
         alpha_out=abar%alpha_out(alphabar)
         alpha_in=abar%alpha_in(alphabar)
         if(alpha_xA/=out3b%alpha2b(alpha_out)) cycle
         lx=out3b%l(alpha_out)
         lb=out3b%lam(alpha_out)
         la=in3b%lam(alpha_in)
         LL=nint(in3b%j(alpha_in))

         YY=0.0_dpreal
        do mJ=-LL,LL
          do mb=-lb,lb ! mb=0 since I choose \hat{r}b as z-direction
            if(mbx /= mJ) cycle
            if((mxA+mb)/=mJ) cycle
            nchla=la**2+la+ma+1
            nchlbx=lbx**2+lbx+mbx+1
            nchlx=lx**2+lx+mxA+1
            nchlb=lb**2+lb+mb+1
            nchLL=LL**2+LL+mJ+1

            CGin=cleb(lbx*2,mbx*2,la*2,ma*2,LL*2,mJ*2)
            CGout=cleb(lx*2,mxA*2,lb*2,mb*2,LL*2,mJ*2)

            yla=sqrt(1.0d0*(2.*la+1)/4./pi) !
            ylb=YLMC(lb,mb)*PL(lb,abs(mb))*exp(-iu*abs(mb)*pi)
            YY=yla*ylb*CGin*CGout+yy
          end do
        end do
        Psi=Psi+YY*smat_abar(alphabar)* iu**(la)  *exp(iu*cph(la)) !!! add coulomb phase ????????????????????????????????????????????????????????????/
      end do

      xsec_alpha_xA=0.0_dpreal

        xsec_alpha_xA=abs(Psi)**2

      xsec_alpha_xA=xsec_alpha_xA*N

      end subroutine
c-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine dsde_cdcc(alphabar,xsec)
c      this subroutine is used to calculate differential cross section
c      the formula is the same as the one in the note
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use channels
        use pot
        use interpolation
        use constants
        use mesh
        use input,only:prior,icf
        implicit none
        integer :: alpha_xA,alpha_out
        real*8 :: xsec
        real*8 :: N
        integer :: irx
        real*8 :: rx
        real*8 :: abs_R
        integer :: LL
        real*8 :: J_out
        integer :: alphabar
        complex*16 :: S_cdcc

        N=mua * mub * kb / 4 / pi**2 / hbarc**4 /ka
        N=N * kx * mux / 8. / pi**3 / hbarc**2
        N=N*10.0_dpreal ! 1fm2=10mb



        xsec=0.0_dpreal





        alpha_out=abar%alpha_out(alphabar)
        J_out=out3b%j(alpha_out)
        LL=nint(J_out)
        alpha_xA=out3b%alpha2b(alpha_out)



        call Smat(alphabar,S_cdcc)

        xsec=abs(S_cdcc)**2 *(2.0_dpreal*LL+1.0_dpreal)/4.0_dpreal/pi
        smat_abar(alphabar)=S_cdcc





        xsec=xsec*N/(2.*lbx+1.)



       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine Smat(alphabar,S_cdcc_sum)
c       this subroutine is used to calculate the R function
c       after performing the sum in alphacdcc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         use channels
         use mesh
         implicit none
         integer :: nchcdcc,alpha_cdcc,alphabar
         complex*16 :: S_cdcc_sum,s_cdcc

         S_cdcc_sum=0.0_dpreal
         do nchcdcc=1, abar%cdccnchmax(alphabar)
           alpha_cdcc=abar%alpha_cdcc(alphabar,nchcdcc)
           call smat_cdcc_loops(alphabar,alpha_cdcc,s_cdcc)
           S_cdcc_sum=S_cdcc_sum+s_cdcc
         end do



       end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine smat_cdcc_loops(alphabar,alpha_cdcc,s_cdcc)
c      this subroutine is used to calculate the R_func with
c      index ain aout acdcc in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use channels
        use mesh
        use scatt
        use systems
        use constants
        use input,only:prior,printf,HM
        implicit none
        integer :: alphabar
        integer :: alpha2b,alphaspect,alpha_out
        integer :: alpha_cdcc
        real*8 :: rx,rbx,rxp
        integer :: irx,irbx,irxp
        complex*16 :: lambda,lambdaNO
        complex*16 :: s_cdcc
        complex*16,dimension(1:nr) :: rho
        integer :: lb,lx
        complex*16 :: N
        real*8,dimension(1:nx,1:nr,1:nr) :: Gabar
        real*8,dimension(1:nx) :: G_alphabar
        complex*16,dimension(1:nx) :: func_in,vpot_in


        alpha_out=abar%alpha_out(alphabar)


        call G_alphaout_alphain_cdcc(alphabar,alpha_cdcc,Gabar)




        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)
        alpha2b=out3b%alpha2b(alpha_out)
        lx=out3b%l(alpha_out)

        rho=0.0_dpreal

        if (printf)   write(11,*) "& alpha_cdcc=",alpha_cdcc,"alpha_out=",alpha_out
        if (printf)   write(13,*) "& alpha_cdcc=",alpha_cdcc,"alpha_out=",alpha_out
        do irxp=1,nr
          rxp=rr(irxp)

          do irbx=1,nr
            rbx=rr(irbx)
            G_alphabar(1:nx)=Gabar(1:nx,irbx,irxp)
              call lambda_func_cdcc(alphabar,alpha_cdcc,G_alphabar,irbx,irxp,lambda)
              rho(irxp)=rho(irxp) + rbx**2 * lambda * rrw(irbx)
          end do

          if (printf) write(11,*) rxp,real(rho(irxp)),aimag(rho(irxp))


        end do


        N=2.**6 * pi**3 * iu**(-lx) /kx/ka/kb
        s_cdcc=0.0_dpreal
        do irx=1,nr
          rx=rr(irx)
          s_cdcc= s_cdcc + N*rx* Gx%re(irx,alpha2b) * rho(irx) *rrw(irx)

        end do




       end subroutine
c-----------------------------------------------------------------------








ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_func_pot_in_cdcc()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh
       use channels

       implicit none
       integer :: nchin




       if(allocated (fun_in)) deallocate(fun_in)
       if(allocated (pot_in)) deallocate(pot_in)
       allocate(fun_in(1:nx,1:nr,1:nr,1:incdcc%nchmax))
       allocate(pot_in(1:nx,1:nr,1:nr,1:incdcc%nchmax))

       pot_in=0.0_dpreal
       fun_in=0.0_dpreal


  !$OMP PARALLEL default(shared) private(nchin)
  !$OMP DO schedule(static)
          do nchin=1, incdcc%nchmax
            call initial_lambda_func_cdcc(nchin)
          end do
  !$OMP END DO
  !$OMP end PARALLEL



      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_lambda_func_cdcc(nchin)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      implicit none
      integer :: irbx,irx,nchin


        do irx=1,nr
          do irbx=1,nr
            call initial_lambda_func_nxloop_cdcc(irbx,irx,nchin)
          end do
        end do
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine initial_lambda_func_nxloop_cdcc(irbx,irx,acdcc)
c       this subroutine is used to calculate the lambda function that
c       appears in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use pot
        use scatt
        use bound
        use channels
        use interpolation
        use mesh
        use input,only:prior,HM
        implicit none
        integer :: irbx,irx
        integer :: acdcc
        integer :: ix
        real*8 :: x
        real*8 :: ra,rbx,rx,rb,p,q
        real*8 :: rbA
        complex*16 :: wfa
        real*8 :: wfbx
        integer ::nchbx,nchaA
        complex*16 :: func_in
        complex*16 :: vpot_in
        integer :: ainspect
        integer :: nrank


        p=massb/(massb+massx)
        q=masst/(masst+massx)
        rx=rr(irx)
        rbx=rr(irbx)

        if (HM) then

          do ix=1,nx

              x=angx(ix)
              ra=sqrt(rx**2 + p**2 * rbx**2 - 2.*p*rx*rbx*x)
              rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
              rbA=sqrt(rx**2 + rbx**2 - 2*rx*rbx*x)

              func_in=0.0_dpreal
              do nrank=1, incdcc%n(acdcc)
               wfa=FFC(ra/hcm,cdccwf%chi(0:irmatch,nrank,acdcc),irmatch+1) / ra
               wfbx=cdccwf%phi(irbx,nrank,acdcc)
               func_in=func_in + wfa*wfbx
              end do

              fun_in(ix,irbx,irx,acdcc)=func_in

          end do

        else


          do ix=1,nx
            x=angx(ix)
            ra=sqrt(rx**2 + p**2 * rbx**2 - 2.*p*rx*rbx*x)
            rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
            rbA=sqrt(rx**2 + rbx**2 - 2*rx*rbx*x)


            func_in=0.0_dpreal
            do nrank=1, incdcc%n(acdcc)
               wfa=FFC(ra/hcm,cdccwf%chi(0:irmatch,nrank,acdcc),irmatch+1) / ra
               wfbx=cdccwf%phi(irbx,nrank,acdcc)
               func_in=func_in + wfa*wfbx
            end do

            fun_in(ix,irbx,irx,acdcc)=func_in

            vpot_in=FFR4(rbx/hcm,vbx(0:irmatch,1),irmatch+1) + FFC(rbA/hcm,UbA,irmatch+1)
            pot_in(ix,irbx,irx,acdcc)=vpot_in
          end do

        end if

       end subroutine
c-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine lambda_func_cdcc(alphabar,alpha_cdcc,G_alphabar,irbx,irx,lambda)
c       this subroutine is used to calculate the DELTA function that
c       appears in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use pot
        use scatt
        use bound
        use channels
        use interpolation
        use mesh
        implicit none
        integer :: alphabar,nchin,alphaspect
        integer :: alpha_out,irbx,irx,alpha_cdcc
        complex*16 :: lambda
        integer :: ix
        complex*16 :: Vpost
        real*8 :: rx,rb,rbx,x,q
        integer :: nchbB ! for UbB potential
        integer :: lb
        complex*16 :: UUbB,wfb
        complex*16,dimension(1:nx) :: func_in,vpot_in
        real*8,dimension(1:nx) :: G_alphabar



        q=masst/(masst+massx)
        alpha_out=abar%alpha_out(alphabar)


        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)

        func_in(1:nx)=fun_in(1:nx,irbx,irx,alpha_cdcc)
        vpot_in(1:nx)=pot_in(1:nx,irbx,irx,alpha_cdcc)


        lambda=0.0_dpreal


        rx=rr(irx)
        rbx=rr(irbx)

        nchbB=out3b%alphaspect(alpha_out)



        do ix=1,nx

          x=angx(ix)
          rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
          wfb=FFC(rb/hcm,wf_b(0:irmatch,alphaspect),irmatch+1)*iu**(-lb)/ rb
          UUbB=FFC(rb/hcm,UbB(0:irmatch,nchbB),1+irmatch)
          vpost=vpot_in(ix) - UUbB


          lambda = lambda + G_alphabar(ix)  * vpost * angw(ix) * func_in(ix) * wfb


        end do

       end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine lambda_func_cdcc_HM(alphabar,alpha_cdcc,G_alphabar,irbx,irx,lambdaNO)
C       this subroutine is used to calculate the DELTA function that
C       appears in the notes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use pot
        use scatt
        use bound
        use channels
        use interpolation
        use mesh
        implicit none
        integer :: alphabar,alphaspect
        integer :: alpha_out,irbx,irx,alpha_cdcc
        complex*16 :: lambdaNO
        integer :: ix
        real*8 :: rx,rb,x,rbx,q
        integer :: nchxA ! for UxA potential
        integer :: LL,lb
        complex*16 :: UUx
        complex*16,dimension(1:nx) :: func_in,vpot_in
        complex*16 :: func,wfb
        real*8,dimension(1:nx) :: G_alphabar

        q=masst/(masst+massx)

        alpha_out=abar%alpha_out(alphabar)


        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)

        func_in(1:nx)=fun_in(1:nx,irbx,irx,alpha_cdcc)


        lambdaNO=0.0_dpreal

        rx=rr(irx)
        rbx=rr(irbx)

        do ix=1,nx

            x=angx(ix)
            rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
            wfb=FFC(rb/hcm,wf_b(0:irmatch,alphaspect),irmatch+1)* iu**(-lb) / rb
            func=G_alphabar(ix)  * angw(ix) * func_in(ix) * wfb
            lambdaNO = lambdaNO + func

        end do

       end subroutine
C-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine readincdccwf()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use channels
       use mesh
       implicit none

       integer :: ira, irbx, nch
       integer :: acdcc

       open (unit=551,file='cdccwf.dat')
       call alpha_cdcc_in()
       allocate(psicdcc(1:nr,0:irmatch,1:incdcc%nchmax))

       psicdcc=0.0_dpreal

         do nch=1, incdcc%nchmax
          read(551,*) acdcc
          write(*,*) acdcc

         do ira=0,irmatch
           do irbx=1, nr
            read(551,*) psicdcc(irbx,ira,nch)

           end do
         end do

        end do
       end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readincdccwf_seprate()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      implicit none
      integer :: ira, irbx, nch,l
      integer :: acdcc1,acdcc2,nmax,n


        open (unit=552,file='cdccwf_phi.dat')
        open (unit=553,file='cdccwf_chi.dat')
        open (unit=554,file='cdccwf_ranks.dat')
        call alpha_cdcc_in()

        do nch=1, incdcc%nchmax
           read(554,*) incdcc%n(nch)
        end do

        nmax=MAXVAL(incdcc%n(1:incdcc%nchmax))


        allocate(cdccwf%phi(1:nr,1:nmax,1:incdcc%nchmax))
        allocate(cdccwf%chi(0:irmatch,1:nmax,1:incdcc%nchmax))

        allocate(psicdcc(1:nr,0:irmatch,1:incdcc%nchmax))


        cdccwf%phi=0.0_dpreal
        cdccwf%chi=0.0_dpreal
        psicdcc=0.0_dpreal

        do nch=1, incdcc%nchmax
           read(552,*)acdcc1
           read(553,*)acdcc2
           write(*,*) "acdcc1=",acdcc1, "acdcc2=",acdcc2
C          write(999,* )nch
            l=incdcc%lam(nch)

          do n=1, incdcc%n(nch)

            do irbx=1, nr
              read(552,*) cdccwf%phi(irbx,n,nch)
               write(567,*) rr(irbx),real(cdccwf%phi(irbx,n,nch))*rr(irbx)
            end do
             write(567,*) "&"

              write(568,*) "&n=",nch
            do ira=0, irmatch
               read(553,*) cdccwf%chi(ira,n,nch)
               write(568,*) ira*hcm,real(cdccwf%chi(ira,n,nch)), imag(cdccwf%chi(ira,n,nch))
            end do



c----test
C           do irbx =1, nr
C              do ira=0, irmatch
C                 psicdcc(irbx,ira,nch) = psicdcc(irbx,ira,nch) + cdccwf%phi(irbx,n,nch)*cdccwf%chi(ira,n,nch)
C
C              end do
C           end do
c---- test

          end do

C---- test
C
C            do ira=0, irmatch
C            do irbx =1, nr
C                write(999,* )psicdcc(irbx,ira,nch)
C            end do
C         end do
C---- end test

        end do
      end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine G_alphaout_alphain_cdcc(alphabar,alpha_cdcc,Gabar)
C     this subroutine is used to calculate the G_alphabar(ix,irbx,irx)
c     coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      use cleb_coefficient
      implicit none
      integer :: irbx,irx,ix
      integer :: alphabar,alpha_cdcc, alpha_out
      integer :: ML,mx,mb
      integer :: nch_L,nch_la,nch_lbx,nch_lx,nch_lb
      integer :: LL,la,lb,lx,minl,l2b
      real*8 ::j_out
      real*8 :: CGin,CGout,Yout,CG_in_out
      real*8,dimension(1:nx) :: Yin
      real*8,dimension(1:nx,1:nr,1:nr) :: Gabar
      real*8,dimension(1:nx) :: yla,ylb,ylx
      real*8 :: func


        Gabar=0.0_dpreal
        alpha_out=abar%alpha_out(alphabar)

        j_out=out3b%j(alpha_out)
        LL=nint(j_out)

        la=incdcc%lam(alpha_cdcc)
        l2b=incdcc%l(alpha_cdcc)

        lx=out3b%l(alpha_out)
        lb=out3b%lam(alpha_out)
        minl=min(LL,la)


        do ML=-minl,minl

          nch_L=LL**2+LL+ML+1
          nch_lbx=l2b**2+l2b+0+1
          nch_la=la**2+la+ML+1
          CGin=cleb(l2b*2,0,la*2,ML*2,LL*2,ML*2)*sqrt( (2.*l2b+1.)/(4.0_dpreal*pi) )



          do mx=-lx,lx
            mb=ML-mx
            if(mb>lb .or. mb<-lb) cycle
              nch_lx=lx**2+lx+mx+1
              nch_lb=lb**2+lb+mb+1
              ylx=Y_out_lx(:,nch_lx)

               CGout=cleb(lx*2,mx*2,lb*2,mb*2,LL*2,ML*2)
               CG_in_out=CGin*CGout
              do irx=1,nr
                do irbx=1,nr
                  yla=Y_in_la(:,irbx,irx,nch_la)
                  ylb=Y_out_lb(:,irbx,irx,nch_lb)


                  do ix=1,nx


                    Yout=ylx(ix) *  ylb(ix)

                    func=CG_in_out*Yout*yla(ix)

                    Gabar(ix,irbx,irx)=Gabar(ix,irbx,irx)+func

                  end do

                end do
              end do

          end do

        end do



       Gabar= Gabar * 8 * pi**2 /(2.*LL+1.)


      end subroutine
c-----------------------------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_in_rbx(ix)
c
c              | 0 |             | rx*sqrt(1-x^2) |
c     vec{rbx}=| 0 |     vec{rx}=|        0       |
c              | rbx|            |       rx*x     |
c
c     vec{ra}=rx - p*rbx
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use spharm
      use precision
      use mesh
      use channels
      implicit none
      real*8 :: p
      integer :: ix
      integer :: irx,irbx,nch
      real*8 :: rbx,rbx_z,ra_z,ra_abs
      real*8 :: rx_z,rx,x
      integer ::  mla
      integer :: la
      real*8,dimension(0:lmax_cdcc,0:lmax_cdcc) :: PL

c-----define coefficients
      p=massb/(massb+massx)
C      q=masst/(masst+massx)


      x=angx(ix)

      do irx=1, nr
        rx=rr(irx)
        do irbx=1,nr

            rbx=rr(irbx)
            rbx_z=rbx
            rx_z=rx*x

            ra_z=rx_z-p*rbx_z
            ra_abs=sqrt(rx**2 + p**2 * rbx**2 - 2*p*rx*rbx*x)
            call PLM(ra_z/ra_abs,lmax_cdcc,lmax_cdcc,lmax_cdcc+1,PL)
            do la=0,lmax_cdcc
               do mla=-la,la

                 nch=la**2+la+mla+1
                 Y_in_la(ix,irbx,irx,nch)=YLMC(la,mla)*PL(la,abs(mla)) ! phi angle is 0 for ra since ra_x is positive

               end do
            end do



        end do ! irbx
      end do  ! irx

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_out_rbx(ix)
c              | 0 |             | rx*sqrt(1-x^2) |
c     vec{rbx}=| 0 |     vec{rx}=|        0       |
c              | rbx|            |       rx*x     |
c      rb=qrx-rbx
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use spharm
      use precision
      use mesh
      use channels
      implicit none
      integer :: ix
      integer :: l,m,irx,irbx,nch
      real*8 :: rbx,rbx_z,rb_z,rb_abs
      real*8 :: rx_z,rx,x
      integer ::  mlb
      integer :: lb
      real*8 :: q
      real*8,dimension(0:lmax,0:lmax) :: PL
c-----define coefficients
C      p=massb/(massb+massx)
      q=masst/(masst+massx)
      x=angx(ix)


c----Ylx
         call PLM(angx(ix),lmax,lmax,lmax+1,PL)

         do l=0,min(lmax,cutl)
           do m=-l,l

              nch=l**2+l+m+1
              Y_out_lx(ix,nch)=YLMC(l,-m)*PL(l,abs(m))*(-1.0_dpreal)**(m)
           end do
         end do

C------Ylb


       do irx=1, nr
         rx=rr(irx)
         do irbx=1,nr

             rbx=rr(irbx)
             rbx_z=rbx
             rx_z=rx*x

             rb_z=q*rx_z-rbx_z
             rb_abs=sqrt(q**2 * rx**2 +  rbx**2 - 2*q*rx*rbx*x)
             call PLM(rb_z/rb_abs,lmax,lmax,lmax+1,PL)
             do lb=0,lmax
                do mlb=-lb,lb

                  nch=lb**2+lb+mlb+1
                  Y_out_lb(ix,irbx,irx,nch)=YLMC(lb,mlb)*PL(lb,abs(mlb)) ! phi angle is 0 for rb since rb_x is positive

                end do
             end do



         end do ! irbx
       end do  ! irx




      end subroutine
c-----------------------------------------------------------------------




      end module
