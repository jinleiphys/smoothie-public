      module iavdwbarbx_lagrange
      use precision
      use angularmesh
      use mesh
      use fuspot
      implicit none

       real*8,target,dimension(:,:,:,:),allocatable :: Y_in_la
       real*8 :: mua,mub,mux,ka,kb,kx


       complex*16,allocatable,dimension(:,:,:,:) :: fun_in
       complex*16,allocatable,dimension(:,:,:,:) :: pot_in


       complex*16,dimension(:,:),allocatable :: Rabar
       real*8,allocatable,dimension(:) :: Wfus

       real*8,target,dimension(:,:,:,:),allocatable ::Y_out_lb
       real*8,target,dimension(:,:),allocatable ::Y_in_lbx
       public NEB_DWBA_Zero_Spin_rbx_lagrange

      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine NEB_DWBA_Zero_Spin_rbx_lagrange()
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
        use input, only:prior,icf
        use lagrange_mesh_source
        use coulfunc
        implicit none
        real*8 :: ecm,ecmb,ecmx
        integer :: counter
        integer :: alpha_xA,nch
        real*8,dimension(1:out2b%nchmax)  :: xsec_alpha_xA
        ! integer :: alphabar,alpha_out
        integer :: nth,ix
        real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kpla
        real*8 :: theta, cth,xsec_sum,t1,t2
        real*8 :: omp_get_wtime
        ! integer :: jp_index
        real*8,dimension(0:nthmax) :: dsdedw_sum,dsdedw
        real*8,dimension(1:inspect%nchmax) :: xsecla,xsecla_lx
        real*8,dimension(1:outspect%nchmax) :: xseclb, xseclb_lx
        real*8 :: eta, rho
        integer :: ifail
        real*8 :: store_size



        allocate(UbA(0:irmatch))



        mua=(massp * masst) * amu / (massp + masst)
        mub=(massb * (massx+masst)) * amu / (massb + massx+masst)
        mux=(massx * masst) * amu / (massx + masst)
        ecm=elab * masst / (masst + massp)
        ka=sqrt(2*mua*ecm/(hbarc**2))



        call ang_mesh()
        call chana()
C       call chana_readin()
        call chanbx()

        !initial Lagrange mesh
        call rmat_ini(nr,rmax)
        allocate(fc(0:lmax),gc(0:lmax),dfc(0:lmax),dgc(0:lmax))
        !

c calculate the icf fusion potential
        if (icf) then
           allocate(Wfus(0:irmatch))
           call icf_fus_pot(Wfus)
        end if

c!!!


        write(*,90)elab*masst/(masst+massp)+qval
90     format("Threshold at Eb=",f9.4)



        if(allocated(Y_in_la)) deallocate(Y_in_la)
        if(allocated(Y_in_lbx)) deallocate(Y_in_lbx)
        if(allocated(Y_out_lb)) deallocate(Y_out_lb)
        allocate(Y_in_la(1:nx,1:nr,1:nr,1:lmax*lmax+2*lmax+1))
        allocate(Y_out_lb(1:nx,1:nr,1:nr,1:lmax*lmax+2*lmax+1))
        allocate(Y_in_lbx(1:nx,1:lbx*lbx+2*lbx+1))

        write(*,91)
91      format("Now initial the value of Y function ")

        store_size=0.0_dpreal
        store_size=store_size+size(Y_in_la) * storage_size(Y_in_la) / 8.0 / 1048576.
        store_size=store_size+size(Y_out_lb) * storage_size(Y_out_lb) / 8.0 / 1048576.
        store_size=store_size+size(Y_in_lbx) * storage_size(Y_in_lbx) / 8.0 / 1048576.
        write(*,102) store_size
102   format("The storage size of Y functions are ",F10.2," MB")

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

C this part to initial kpla
!$OMP PARALLEL default(shared) private(nth,theta,cth)
!$OMP DO schedule(dynamic)
      do nth=0,nthmax
         theta=thmin+nth*thinc
         cth=cos(theta*pi/180.)
         call PLM(cth,lmax,lmax,lmax+1,kpla(0:lmax,0:lmax,nth))
      end do
!$OMP END DO
!$OMP end PARALLEL





        do counter=1,neb
          ecmb=necmb(counter)
          ecmx=ecm+qval-ecmb
          kb=sqrt(2*mub*ecmb/(hbarc**2))
          kx=sqrt(2*mux*abs(ecmx)/(hbarc**2))

          eta=zx*zt*e2*mux/hbarc/hbarc/kx
          rho=hcm*(irmatch)*kx
          call coul90(rho,eta,zero,lmax,fc,gc,dfc,dgc,0,ifail)   ! call for normalization of  inhomogeneous method

          write(*,92) ecmb
92        format('Now calculating the x-sec at Eb=',F7.3)


          if (nbin(counter)==0) then
             call screenchanbrbx(counter,ecmb)
          else
             call bbin(counter)
          end if



          call greenfunc(counter,ecmx) ! one could remove this, but add lines for the ubx potential,

          if(counter==1) then 
             printpot=.true.
          else 
             printpot=.false.
          end if
          call potr('t',counter,zb*zt,0.0_dpreal,0)
          UbA=V

C          call boundxA()
C          call OMP_SET_NUM_THREADS(1)


          if (prior) write(*,93)
          if (.not. prior) then
           write(*,94)
          end if
 93      format("Using finite range DWBA in IVA-prior form")
 94      format("Using finite range DWBA in IVA-post form")



          call initial_func_pot_in_rbx()
          write(*,96)
96        format( "calculating cross section alpha_xA dependence ",/)
          t1=omp_get_wtime()

          xsec_sum=0.0_dpreal
          xsec_alpha_xA=0.0_dpreal
          dsdedw_sum=0.0_dpreal
          xsecla=0.0_dpreal
          xseclb=0.0_dpreal
          printpot=.true.
          do alpha_xA=1, out2b%nchmax
              if(alpha_xA /= 1) printpot=.false.
              call sigma_J_pi_rbx(alpha_xA,xsec_alpha_xA(alpha_xA),kpla,dsdedw,xsecla_lx,xseclb_lx)
              xsec_sum=xsec_sum+xsec_alpha_xA(alpha_xA)
              dsdedw_sum=dsdedw_sum+dsdedw
              xsecla=xsecla+xsecla_lx
              xseclb=xseclb+xseclb_lx
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


          write(*,100)
          do nth=0,nthmax
             theta=thmin+nth*thinc
             write(*,20)theta,dsdedw_sum(nth)
             write(21,*)theta,dsdedw_sum(nth)
          end do

          write(21,*)"&"

          do nch=1,inspect%nchmax
             write(23,*)inspect%lam(nch),xsecla(nch)
          end do
            write(23,*)"&"
          do nch=1,outspect%nchmax
            write(24,*)outspect%lam(nch),xseclb(nch)
          end do
          write(24,*)"&"



10    format('alpha_xA = ',I3,5x,'X-sec = ',2F9.4)
20    format('theta =',F7.2,5x,'X-sec = ',2F9.4)


        end do
100   format("The double differential cross sections",/)
101   format(20x,'(CPU  time =',F12.2,2x,'seconds)')
120   format('--------------------------------------------------------')
       end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sigma_J_pi_rbx(alpha_xA,xsec_sum,kpla,dsdedw,xsec_la,xsec_lb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      implicit none

       integer :: abarmax,alphabar,alpha_xA
      !  integer :: alpha_out
       real*8,allocatable,dimension(:) :: xsec
       real*8 :: xsec_sum
       integer :: nth,nch,nchla,nchlb
       real*8,dimension(0:nthmax) :: dsdedw
       real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kpla
       real*8,dimension(1:inspect%nchmax) :: xsec_la
       real*8,dimension(1:outspect%nchmax) :: xsec_lb



        dsdedw=0.0_dpreal
        xsec_sum=0.0_dpreal
        call abar_index(alpha_xA,abarmax)
        if(abarmax==0) return






        allocate(xsec(1:abarmax))
        allocate(Rabar(1:nr,1:abarmax))

!$OMP PARALLEL default(shared) private(alphabar)
!$OMP DO schedule(dynamic)
        do alphabar=1,abarmax
          call dsde_rbx(alphabar,xsec(alphabar))
        end do
!$OMP END DO
!$OMP end PARALLEL



          xsec_la=0.0_dpreal
          xsec_lb=0.0_dpreal
         do alphabar=1,abarmax
            xsec_sum=xsec_sum+xsec(alphabar)
            nchla=in3b%alphaspect(abar%alpha_in(alphabar))
            nchlb=out3b%alphaspect(abar%alpha_out(alphabar))
            do nch=1,inspect%nchmax
              if (nch==nchla)  xsec_la(nch)=xsec_la(nch)+xsec(alphabar)
            end do
            do nch=1,outspect%nchmax
              if (nch==nchlb)  xsec_lb(nch)=xsec_lb(nch)+xsec(alphabar)
            end do
         end do



          write(*,20) xsec_sum



 20   format("xsec=",f19.5,/)


!$OMP PARALLEL default(shared) private(nth)
!$OMP DO schedule(dynamic)
       do nth=0,nthmax

          call dsdedwFR_rbx(dsdedw(nth),kpla(:,:,nth),abarmax,alpha_xA)

       end do
!$OMP END DO
!$OMP end PARALLEL

        write(331,*)"&alpha_xA=",alpha_xA
       do nth=0, nthmax
          write(331,*)thmin+nth*thinc, dsdedw(nth)
       end do
      deallocate(Rabar)


      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dsdedwFR_rbx(xsec,kpla,abarmax,alpha_xA)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use precision
      implicit none
      real*8 :: xsec,xsec_alpha_xA
      integer :: mxA,mbx,lx,alpha_xA,abarmax
      real*8,dimension(0:lmax,0:lmax) :: kpla


      xsec=0.0_dpreal
c      do alpha_xA=1,out2b%nchmax
        lx=out2b%l(alpha_xA)
        do mxA=-lx,lx
          do mbx=-lbx,lbx
            call psiintegration_rbx(alpha_xA,abarmax,mbx,mxA,kpla,xsec_alpha_xA)
            xsec=xsec+xsec_alpha_xA
          end do
        end do
c      end do

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psiintegration_rbx(alpha_xA,abarmax,mbx,mxA,PL,xsec_alpha_xA)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use precision
      use constants
      use spharm
      use interpolation
      use pot
      use cleb_coefficient
      use input,only:icf
      implicit none
      real*8 :: xsec_alpha_xA
      integer :: alphabar,alpha_xA,abarmax
      integer :: mxA,mbx
      complex*16,dimension(1:nr) :: Psi,Psi_parallel
      integer :: irx
      real*8 :: rx
      real*8,dimension(0:irmatch) :: Wx
      real*8 :: N,abs_R
      real*8,dimension(0:lmax,0:lmax) :: PL
      N=-mua * mub * kb / 4 / pi**3 / hbarc**4 /ka /(2.0_dpreal*lbx+1.0_dpreal)
      N=N*10.0_dpreal ! 1fm2=10mb


      Psi=0.0_dpreal
      do alphabar=1,abarmax
        CALL psiintegration_parallel(alpha_xA,alphabar,mbx,mxA,PL,Psi_parallel)
        Psi=Psi+Psi_parallel
      end do

      Wx=aimag(UxA(0:irmatch,1))
      if(icf) Wx=Wfus ! for incomplete fusion
      xsec_alpha_xA=0.0_dpreal
      do irx=1,nr
        rx=rr(irx)
        abs_R= abs(Psi(irx))**2
        xsec_alpha_xA=xsec_alpha_xA + rx**2 * FFR4(rx/hcm,Wx,irmatch+1) * abs_R  * rrw(irx)
      end do
      xsec_alpha_xA=xsec_alpha_xA*N

      end subroutine
c-----------------------------------------------------------------------

       subroutine psiintegration_parallel(alpha_xA,alphabar,mbx,mxA,PL,Psi)
         use channels
         use precision
         use constants
         use spharm
         use interpolation
         use pot
         use cleb_coefficient
         use input,only:icf
         implicit none
         integer :: alphabar,alpha_in,alpha_out,alpha_xA
         integer :: lx,lb,LL,la
         integer :: mxA,mbx,mb,mJ,ma
         integer :: nchla,nchlb,nchlx,nchlbx,nchLL
         complex*16,dimension(1:nr) :: Psi
         real*8 :: Yla,Ylb,YY
         real*8,dimension(0:lmax,0:lmax) :: PL
         real*8 :: CGin,CGout
         ma=0
         alpha_out=abar%alpha_out(alphabar)
         alpha_in=abar%alpha_in(alphabar)
         lx=out3b%l(alpha_out)
         lb=out3b%lam(alpha_out)
         la=in3b%lam(alpha_in)
         LL=nint(in3b%j(alpha_in))
         yla=sqrt(1.0d0*(2.*la+1)/4./pi)
         YY=0.0_dpreal
         do mJ=-LL,LL
          if(mbx /= mJ) cycle
          do mb=-lb,lb
            if((mxA+mb)/=mJ) cycle
            nchla=la**2+la+ma+1
            nchlbx=lbx**2+lbx+mbx+1
            nchlx=lx**2+lx+mxA+1
            nchlb=lb**2+lb+mb+1
            nchLL=LL**2+LL+mJ+1

            CGin=cleb(lbx*2,mbx*2,la*2,ma*2,LL*2,mJ*2)
            CGout=cleb(lx*2,mxA*2,lb*2,mb*2,LL*2,mJ*2)


            ylb=YLMC(lb,mb)*PL(lb,abs(mb))*exp(-iu*abs(mb)*pi)
            YY=yla*ylb*CGin*CGout+yy
          end do
         end do

         Psi=YY*Rabar(:,alphabar)




       end subroutine
c-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine dsde_rbx(alphabar,xsec)
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
        real*8,dimension(0:irmatch) :: Wx
        complex*16,dimension(1:nr) :: R_func
        real*8 :: abs_R
        integer :: LL
        real*8 :: J_out
        integer :: alphabar


        N=-mua * mub * kb / 4 / pi**3 / hbarc**4 /ka
        N=N*10.0_dpreal ! 1fm2=10mb



        xsec=0.0_dpreal





        alpha_out=abar%alpha_out(alphabar)
        J_out=out3b%j(alpha_out)
        LL=nint(J_out)
        alpha_xA=out3b%alpha2b(alpha_out)
        Wx=aimag(UxA(0:irmatch,alpha_xA))

        if(icf) Wx=Wfus ! for incomplete fusion

        call Rfunc_rbx(alphabar,R_func)
        Rabar(:,alphabar)=R_func
        do irx=1,nr
          rx=rr(irx)
          abs_R= abs(R_func(irx))**2 * (2.0_dpreal*LL+1.0_dpreal)/4.0_dpreal/pi
          xsec=xsec + rx**2 * FFR4(rx/hcm,Wx,irmatch+1) * abs_R  * rrw(irx)
        end do



      xsec=xsec*N/(2.*lbx+1.)



       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine Rfunc_rbx(alphabar,R_func)
c      this subroutine is used to calculate the R_func in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use channels
        use mesh
        use scatt
        use systems
        use constants
        use input,only:prior,printf
        use lagrange_mesh_source
        use pot, only: UxA
        use interpolation
        implicit none
        integer :: alphabar
        integer :: alpha2b,alphaspect,alpha_out,alpha_in
        real*8 :: rx,rbx,rxp
        integer :: irx,irbx,irxp
        complex*16 :: lambda,lambdaNO
        complex*16,dimension(1:nr) :: R_func,R_funcNO
        complex*16,dimension(1:nr) :: rho
        integer :: lb,lx,la
        ! real*8 :: N
        real*8,dimension(1:nx,1:nr,1:nr) :: Gabar
        real*8,dimension(1:nx) :: G_alphabar
        ! complex*16,dimension(1:nx) :: func_in,vpot_in
        complex*16,dimension(1:nr) :: cpot,csou
        complex*16 :: smat
        real*8 :: hm,ecmx,eta
        real*8 :: v_fact

        v_fact=sqrt(mua*kx/mux/ka) * kx / 4.0_dpreal / pi

        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)

        call G_alphaout_alphain_rbx(alphabar,Gabar)




        la=in3b%lam(alpha_in)
        lb=out3b%lam(alpha_out)
        lx=out3b%l(alpha_out)
        
C       write(*,*) "la=",la,"lb",lb,"lx=",lx
        alphaspect=out3b%alphaspect(alpha_out)
        alpha2b=out3b%alpha2b(alpha_out)

        rho=0.0_dpreal
        R_func=0.0_dpreal
        R_funcNO=0.0_dpreal

        if (printf)   write(11,*) "& alpha_in=",alpha_in,"alpha_out=",alpha_out
        if (printf)   write(13,*) "& alpha_in=",alpha_in,"alpha_out=",alpha_out
        do irxp=1,nr
          rxp=rr(irxp)

          do irbx=1,nr
            rbx=rr(irbx)
            G_alphabar(1:nx)=Gabar(1:nx,irbx,irxp)



            if(prior) then

               call lambda_func_prior_rbx(alphabar,G_alphabar,irbx,irxp,lambda,lambdaNO)
              rho(irxp)=rho(irxp) + rbx**2 * lambda * rrw(irbx)
              R_funcNO(irxp)=R_funcNO(irxp) + rbx**2 * lambdaNO * rrw(irbx)

            else
              call lambda_func_rbx(alphabar,G_alphabar,irbx,irxp,lambda)
              rho(irxp)=rho(irxp) + rbx**2 * lambda * rrw(irbx)
            end if
          end do



        
          if (printf) write(13,*) rxp,real(R_funcNO(irxp)),aimag(R_funcNO(irxp))

        end do

        do irx=1, nr
         rx=rr(irx)
         cpot(irx) = FFC(rx/hcm,UxA(0:irmatch,alpha2b),irmatch+1)
         csou(irx) = -rho(irx) * rx * 2.**4 * pi**2/ka/kb
         
          if (printf) write(11,*) rx,real(csou(irx)),aimag(csou(irx))
         
        end do
        
        


        R_func=0.0_dpreal
        eta=zx*zt*e2*mux/hbarc/hbarc/kx
        hm=hbarc**2/(2.0_dpreal*mux)
        ecmx= kx**2 * hbarc**2 / 2.0 / mux

        call rmat_inho(nr,rmax,cpot,csou,ecmx,eta,hm,lx,smat,R_func)

!         write(51,111)la,lb,lx, v_fact*real(smat),v_fact*aimag(smat)
!         write(52,112) la,lb,lx, v_fact*abs(smat)
! 111     format("la=",I3," lb=",I3," lx=", I3, " S-mat:","(",f12.4,",",f12.4,")")
! 112     format("la=",I3," lb=",I3," lx=", I3, " |S-mat|:",f12.4)


        do irx=1,nr
          rx=rr(irx)
          if(prior) R_func(irx)=R_func(irx)/rx+R_funcNO(irx)* 2.**4 * pi**2/ka/kb
          if (printf)   write(12,*) rx,real(R_func(irx)),aimag(R_func(irx))
        end do
        if (printf)   write(12,*)"&"
       end subroutine
c-----------------------------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_func_pot_in_rbx()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh
       use channels

       implicit none
       integer :: nchin
       real*8 :: store_size




       if(allocated (fun_in)) deallocate(fun_in)
       if(allocated (pot_in)) deallocate(pot_in)
       allocate(fun_in(1:nx,1:nr,1:nr,1:in3b%nchmax))
       allocate(pot_in(1:nx,1:nr,1:nr,1:in3b%nchmax))

       store_size=0.0_dpreal
       store_size=store_size+size(fun_in) * storage_size(fun_in) / 8.0 / 1048576.
       store_size=store_size+size(pot_in) * storage_size(pot_in) / 8.0 / 1048576.
       write(*,130) store_size
130   format("The storage size of function and potential are ",F10.2," MB")

       pot_in=0.0_dpreal
       fun_in=0.0_dpreal


  !$OMP PARALLEL default(shared) private(nchin)
  !$OMP DO schedule(static)
          do nchin=1, in3b%nchmax
            call initial_lambda_func_rbx(nchin)
          end do
  !$OMP END DO
  !$OMP end PARALLEL



      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_lambda_func_rbx(nchin)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      implicit none
      integer :: irbx,irx,nchin

        do irx=1,nr
          do irbx=1,nr
            call initial_lambda_func_nxloop_rbx(irbx,irx,nchin)
          end do
        end do
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine initial_lambda_func_nxloop_rbx(irbx,irx,alpha_in)
c       this subroutine is used to calculate the lambda function that
c       appears in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use pot
        use scatt
        use bound
        use channels
        use interpolation
        use mesh
        use input,only:prior
        implicit none
        integer :: irbx,irx
        integer :: alpha_in
        integer :: ix
        integer :: la
        real*8 :: x
        real*8 :: ra,rbx,rx,rb,p,q
        real*8 :: rbA
        complex*16 :: wfa
        real*8 :: wfbx
        integer ::nchbx,nchaA
        complex*16 :: func_in
        complex*16 :: vpot_in
        integer :: ainspect



        la=in3b%lam(alpha_in)
        ainspect=in3b%alphaspect(alpha_in)


        p=massb/(massb+massx)
        q=masst/(masst+massx)
        rx=rr(irx)
        rbx=rr(irbx)

        if (prior) then

          nchbx=in3b%alpha2b(alpha_in)
          nchaA=in3b%alphaspect(alpha_in)

          do ix=1,nx

              x=angx(ix)
              ra=sqrt(rx**2 + p**2 * rbx**2 - 2.*p*rx*rbx*x)
              rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
              rbA=sqrt(rx**2 + rbx**2 - 2*rx*rbx*x)

              wfa=FFC(ra/hcm,wf_a(0:irmatch,ainspect),irmatch+1)*iu**la / ra
              wfbx=FFR4(rbx/hcm,phi_bx(0:irmatch,nchbx),irmatch+1) / rbx
              func_in=wfa*wfbx
              fun_in(ix,irbx,irx,alpha_in)=func_in
              vpot_in=FFC(rbA/hcm,UbA,irmatch+1)-FFC(ra/hcm,UaA(0:irmatch,nchaA),irmatch+1)
              pot_in(ix,irbx,irx,alpha_in)=vpot_in

          end do

        else

          nchbx=in3b%alpha2b(alpha_in)
          do ix=1,nx
            x=angx(ix)
            ra=sqrt(rx**2 + p**2 * rbx**2 - 2.*p*rx*rbx*x)
            rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
            rbA=sqrt(rx**2 + rbx**2 - 2*rx*rbx*x)


            wfa=FFC(ra/hcm,wf_a(0:irmatch,ainspect),irmatch+1) * iu**la / ra
            wfbx=FFR4(rbx/hcm,phi_bx(0:irmatch,nchbx),irmatch+1) / rbx
            func_in=wfa*wfbx
            fun_in(ix,irbx,irx,alpha_in)=func_in

            vpot_in=FFR4(rbx/hcm,vbx(0:irmatch,nchbx),irmatch+1) + FFC(rbA/hcm,UbA,irmatch+1)
            pot_in(ix,irbx,irx,alpha_in)=vpot_in
          end do

        end if

       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine lambda_func_rbx(alphabar,G_alphabar,irbx,irx,lambda)
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
        integer :: alphabar,alphaspect
        integer :: alpha_out,irbx,irx,alpha_in
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
        alpha_in=abar%alpha_in(alphabar)

        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)

        func_in(1:nx)=fun_in(1:nx,irbx,irx,alpha_in)
        vpot_in(1:nx)=pot_in(1:nx,irbx,irx,alpha_in)


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
       subroutine lambda_func_prior_rbx(alphabar,G_alphabar,irbx,irx,lambda,lambdaNO)
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
        integer :: alpha_out,irbx,irx,alpha_in
        complex*16 :: lambda,lambdaNO
        integer :: ix
        real*8 :: rx,rb,x,rbx,q
        integer :: nchxA ! for UxA potential
        integer :: lb
        complex*16 :: UUx
        complex*16,dimension(1:nx) :: func_in,vpot_in
        complex*16 :: func,wfb
        real*8,dimension(1:nx) :: G_alphabar

        q=masst/(masst+massx)

        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)



        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)


        func_in(1:nx)=fun_in(1:nx,irbx,irx,alpha_in)
        vpot_in(1:nx)=pot_in(1:nx,irbx,irx,alpha_in)

        lambda=0.0_dpreal
        lambdaNO=0.0_dpreal





        rx=rr(irx)
        rbx=rr(irbx)

        nchxA=out3b%alpha2b(alpha_out)
        UUx=FFC(rx/hcm,UxA(0:irmatch,nchxA),irmatch+1)


        do ix=1,nx

            x=angx(ix)
            rb=sqrt(q**2 * rx**2 + rbx**2 - 2.*q*rx*rbx*x)
            wfb=FFC(rb/hcm,wf_b(0:irmatch,alphaspect),irmatch+1)* iu**(-lb) / rb

            func=G_alphabar(ix)  * angw(ix) * func_in(ix) * wfb
            lambda = lambda + func  * (vpot_in(ix) + UUx)
            lambdaNO = lambdaNO + func

        end do


       end subroutine
C-----------------------------------------------------------------------






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_in_rbx(ix)
c
c             | 0 |              | rbx*sqrt(1-x^2) |
c     vec{rx}=| 0 |     vec{rbx}=|        0       |
c             | rx|              |       rbx*x     |
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
      integer ::  mla,mbx
      integer :: la
      real*8,dimension(0:lmax,0:lmax) :: PL
      real*8,dimension(0:lbx,0:lbx)  :: plbx

c-----define coefficients
      p=massb/(massb+massx)
C      q=masst/(masst+massx)


      x=angx(ix)

C Ylbx

         call PLM(angx(ix),lbx,lbx,lbx+1,plbx)


           do mbx=-lbx,lbx

              nch=lbx**2+lbx+mbx+1
              Y_in_lbx(ix,nch)=YLMC(lbx,mbx)*PLbx(lbx,abs(mbx))
           end do




C Yla
      do irx=1, nr
        rx=rr(irx)
        do irbx=1,nr

            rbx=rr(irbx)
            rbx_z=rbx*x
            rx_z=rx

            ra_z=rx_z-p*rbx_z
            ra_abs=sqrt(rx**2 + p**2 * rbx**2 - 2*p*rx*rbx*x)
            call PLM(ra_z/ra_abs,lmax,lmax,lmax+1,PL)
            do la=0,lmax
               do mla=-la,la

                 nch=la**2+la+mla+1
                 Y_in_la(ix,irbx,irx,nch)=YLMC(la,mla)*PL(la,abs(mla)) *  (-1.0_dpreal)**(mla) ! phi angle is pi for ra since ra_x is negative

               end do
            end do



        end do ! irbx
      end do  ! irx

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_out_rbx(ix)
c             | 0 |              | rbx*sqrt(1-x^2) |
c     vec{rx}=| 0 |     vec{rbx}=|        0       |
c             | rx|             |       rbx*x     |
c      rb=qrx-rbx
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use spharm
      use precision
      use mesh
      use channels
      implicit none
      integer :: ix
      integer :: irx,irbx,nch
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
C        call PLM(angx(ix),min(lmax,cutl),min(lmax,cutl),min(lmax,cutl)+1,PLx)
C
C        do l=0,min(lmax,cutl)
C          do m=-l,l
C
C             nch=l**2+l+m+1
C             Y_out_lx(ix,nch)=YLMC(l,-m)*PLx(l,abs(m))*(-1.0_dpreal)**(m)
C          end do
C        end do
C

C------Ylb


       do irx=1, nr
         rx=rr(irx)
         rx_z=rx
         do irbx=1,nr

             rbx=rr(irbx)
             rbx_z=rbx*x


             rb_z=q*rx_z-rbx_z
             rb_abs=sqrt(q**2 * rx**2 +  rbx**2 - 2*q*rx*rbx*x)
             call PLM(rb_z/rb_abs,lmax,lmax,lmax+1,PL)
             do lb=0,lmax
                do mlb=-lb,lb

                  nch=lb**2+lb+mlb+1
                  Y_out_lb(ix,irbx,irx,nch)=YLMC(lb,-mlb)*PL(lb,abs(mlb)) !* (-1.0_dpreal)**(mlb) * (-1.0_dpreal)**(mlb) ! phi angle is pi for rb since rb_x is negative
                end do
             end do



         end do ! irbx
       end do  ! irx



      end subroutine
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine G_alphaout_alphain_rbx(alphabar,Gabar)
C     this subroutine is used to calculate the G_alphabar(ix,irbx,irx)
c     coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      use cleb_coefficient
      implicit none
      integer :: irbx,irx,ix
      integer :: alphabar,alpha_in, alpha_out
      integer :: ML,mbx,ma
      integer :: nch_L,nch_la,nch_lbx,nch_lb
      integer :: LL,la,lb,lx,minl
      real*8 ::j_out
      real*8 :: CGin,CGout,CG_in_out
      real*8,dimension(1:nx,1:nr,1:nr),target :: Gabar
      real*8,pointer,dimension(:) :: GG,yla,ylb,ylbx
      real*8 :: ylx


        Gabar=0.0_dpreal
        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)

        j_out=out3b%j(alpha_out)
        LL=nint(j_out)
        la=in3b%lam(alpha_in)
        lx=out3b%l(alpha_out)
        lb=out3b%lam(alpha_out)
        minl=min(LL,lb)
        ylx=sqrt( (2.*lx+1.)/(4.0_dpreal*pi) )


        do ML=-minl,minl
          nch_L=LL**2+LL+ML+1
          CGout=cleb(lx*2,0,lb*2,ML*2,LL*2,ML*2)*ylx
          do ma=-la,la
            do mbx=-lbx,lbx
            if ( (ma+mbx) /= mL) cycle
              nch_lbx=lbx**2+lbx+mbx+1
              nch_lb=lb**2+lb+mL+1
              nch_la=la**2+la+ma+1
              ylbx=>Y_in_lbx(:,nch_lbx)
              CGin=cleb(lbx*2,mbx*2,la*2,ma*2,LL*2,ML*2)
              CG_in_out=CGin*CGout
              do irx=1,nr
                do irbx=1,nr
                  yla=>Y_in_la(:,irbx,irx,nch_la)
                  ylb=>Y_out_lb(:,irbx,irx,nch_lb)
                  GG=>Gabar(:,irbx,irx)

                  do ix=1,nx


                    GG(ix)=GG(ix)+CG_in_out*ylb(ix)*yla(ix)*ylbx(ix)

                  end do

                end do
              end do
            end do
          end do

        end do



       Gabar= Gabar * 8 * pi**2 /(2.*LL+1.)


      end subroutine
c-----------------------------------------------------------------------




      end module
