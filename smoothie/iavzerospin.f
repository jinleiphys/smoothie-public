      module IAVzerospin
      use precision
      use angularmesh
      use mesh
      use fuspot
      implicit none

       real*8,dimension(:,:),target,allocatable :: Y_out_lx
       real*8,target,dimension(:,:,:,:),allocatable :: Y_in_la,Y_in_lbx
       real*8 :: mua,mub,mux,ka,kb,kx



       complex*16,allocatable,dimension(:,:,:,:) :: fun_in
       complex*16,allocatable,dimension(:,:,:,:) :: pot_in


       complex*16,dimension(:,:),allocatable :: Rabar
       complex*16,dimension(:,:),allocatable :: Rla
       real*8,allocatable,dimension(:) :: Wfus


      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine NEB_DWBA_Zero_Spin()
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
        use interpolation
        implicit none
        real*8 :: ecm,ecmb,ecmx
        integer :: counter
        integer :: alpha_xA
        real*8,dimension(1:out2b%nchmax)  :: xsec_alpha_xA
        ! integer :: alphabar,alpha_out
        integer :: nth,ix,ir
        real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kpla
        real*8 :: theta, cth,xsec_sum,t1,t2
        real*8 :: omp_get_wtime
        ! integer :: jp_index
        real*8,dimension(0:nthmax) :: dsdedw_sum,dsdedw
        integer :: lxmax1,nch, la
        real*8,dimension(1:inspect%nchmax) :: xsecla,xsecla_lx
        real*8,dimension(1:outspect%nchmax) :: xseclb, xseclb_lx
        real*8 :: store_size 

        allocate(UbA(0:irmatch))



        mua=(massp * masst) * amu / (massp + masst)
        mub=(massb * (massx+masst)) * amu / (massb + massx+masst)
        mux=(massx * masst) * amu / (massx + masst)
        ecm=elab * masst / (masst + massp)
        ka=sqrt(2*mua*ecm/(hbarc**2))
       
        if(prior .eqv. .false.) then 
           call ang_mesh_post()
        else
           call ang_mesh()
        end if 
        call chana()
        call chanbx()

c calculate the icf fusion potential
        if (icf) then
           allocate(Wfus(0:irmatch))
           call icf_fus_pot(Wfus)
        end if

c!!!


        write(*,90)elab*masst/(masst+massp)+qval
90     format("Threshold at Eb=",f9.4)

        lxmax1=min(lmax,lxmax)
        if(allocated(Y_in_la)) deallocate(Y_in_la)
        if(allocated(Y_in_lbx)) deallocate(Y_in_lbx)
        if(allocated(Y_out_lx)) deallocate(Y_out_lx)
        allocate(Y_in_la(1:nx,1:nr,1:nr,1:lmax*lmax+2*lmax+1))
        allocate(Y_in_lbx(1:nx,1:nr,1:nr,1:lbxmax*lbxmax+2*lbxmax+1))
        allocate(Y_out_lx(1:nx,1:lxmax1*lxmax1+2*lxmax1+1))  !  cutl ????
        write(*,91)
91      format("Now initial the value of Y functions ")
        store_size=0.0_dpreal
        store_size=store_size+size(Y_in_la) * storage_size(Y_in_la) / 8.0 / 1048576.
        store_size=store_size+size(Y_in_lbx) * storage_size(Y_in_lbx) / 8.0 / 1048576.
        store_size=store_size+size(Y_out_lx) * storage_size(Y_out_lx) / 8.0 / 1048576.
        write(*,102) store_size
102   format("The storage size of Y functions are ",F10.2," MB")
        t1=omp_get_wtime()
!$OMP PARALLEL default(shared) private(ix)
!$OMP DO schedule(dynamic)
        do ix=1,nx
          call initial_Y_out(ix)
          call initial_Y_in(ix)
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


        
        allocate(Rla(1:nr,0:lmax))


        do counter=1,neb
          ecmb=necmb(counter)
          ecmx=ecm+qval-ecmb
          kb=sqrt(2*mub*ecmb/(hbarc**2))
          kx=sqrt(2*mux*abs(ecmx)/(hbarc**2))

          Rla=0.0_dpreal 
          write(*,92) ecmb
92        format('Now calculating the x-sec at Eb=',F7.3)


          
          
          
          
          if (nbin(counter)==0) then
             call screenchanb(counter,ecmb)
          else
             call bbin_rb(counter)
          end if

          
          
          call greenfunc(counter,ecmx)

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
           if (nx<100) then
              write(*,95)
C              stop
           end if

         end if
 93      format("Using finite range DWBA in IVA-prior form")
 94      format("Using finite range DWBA in IVA-post form")
 95      format("error! in the post form nx should be very large")


          call initial_func_pot_in()
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
              call sigma_J_pi(alpha_xA,xsec_alpha_xA(alpha_xA),kpla,dsdedw,xsecla_lx,xseclb_lx)
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
            
            
         if (counter==1) then 
           do la=0, lmax
              write(14,*) "&la=",la
              write(15,*) "&la=",la
              do ir=1, nr 
                 write(14,*) rr(ir), abs(Rla(ir,la))
                 write(15,*) rr(ir), abs(Rla(ir,la))*FFR4(rr(ir)/hcm,aimag(UxA(0:irmatch,1)),irmatch+1)
              end do 
           
            end do 
         
         end if 



10    format('alpha_xA = ',I3,5x,'X-sec = ',2F9.4)
20    format('theta =',F7.2,5x,'X-sec = ',2F9.4)


        end do
100   format("The double differential cross sections",/)
101   format(20x,'(CPU  time =',F12.2,2x,'seconds)')
120   format('--------------------------------------------------------')
       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_func_pot_in()
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
130    format("The storage size of function and potential are ",F10.2," MB")

       pot_in=0.0_dpreal
       fun_in=0.0_dpreal


  !$OMP PARALLEL default(shared) private(nchin)
  !$OMP DO schedule(static)
          do nchin=1, in3b%nchmax
            call initial_lambda_func(nchin)
          end do
  !$OMP END DO
  !$OMP end PARALLEL



      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dsdedwFR(xsec,kpla,abarmax,alpha_xA)
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
            call psiintegration(alpha_xA,abarmax,mbx,mxA,kpla,xsec_alpha_xA)
            xsec=xsec+xsec_alpha_xA
          end do
        end do
c      end do

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psiintegration(alpha_xA,abarmax,mbx,mxA,PL,xsec_alpha_xA)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sigma_J_pi(alpha_xA,xsec_sum,kpla,dsdedw,xsec_la,xsec_lb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      implicit none

       integer :: abarmax,alphabar,alpha_xA
      !  integer :: alpha_out
       real*8,allocatable,dimension(:) :: xsec
       real*8 :: xsec_sum
       integer :: nth
       real*8,dimension(0:nthmax) :: dsdedw
       real*8,dimension(0:lmax,0:lmax,0:nthmax) :: kpla
       real*8,dimension(1:inspect%nchmax) :: xsec_la
       real*8,dimension(1:outspect%nchmax) :: xsec_lb
       integer :: nch,nchla,nchlb


        dsdedw=0.0_dpreal
        xsec_sum=0.0_dpreal
        call abar_index(alpha_xA,abarmax)
        if(abarmax==0) return






        allocate(xsec(1:abarmax))
        allocate(Rabar(1:nr,1:abarmax))

!$OMP PARALLEL default(shared) private(alphabar)
!$OMP DO schedule(dynamic)
        do alphabar=1,abarmax
          call dsde(alphabar,xsec(alphabar))
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

          call dsdedwFR(dsdedw(nth),kpla(:,:,nth),abarmax,alpha_xA)

       end do
!$OMP END DO
!$OMP end PARALLEL

      deallocate(Rabar)


      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine dsde(alphabar,xsec)
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
        integer :: alphabar, alpha_in 


        N=-mua * mub * kb / 4 / pi**3 / hbarc**4 /ka
        N=N*10.0_dpreal ! 1fm2=10mb



        xsec=0.0_dpreal





        alpha_out=abar%alpha_out(alphabar)
        J_out=out3b%j(alpha_out)
        LL=nint(J_out)
        alpha_xA=out3b%alpha2b(alpha_out)
        Wx=aimag(UxA(0:irmatch,alpha_xA))

        if(icf) Wx=Wfus ! for incomplete fusion

        call Rfunc(alphabar,R_func)
        Rabar(:,alphabar)=R_func
        do irx=1,nr
          rx=rr(irx)
          abs_R= abs(R_func(irx))**2 * (2.0_dpreal*LL+1.0_dpreal)/4.0_dpreal/pi
          xsec=xsec + rx**2 * FFR4(rx/hcm,Wx,irmatch+1) * abs_R  * rrw(irx)
        end do



      xsec=xsec*N/(2.*lbx+1.)
      
      
      
      alpha_in=abar%alpha_in(alphabar)
      
C     if(alpha_in == 1) write(*,*) "x-sec=",xsec



       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine Rfunc(alphabar,R_func)
c      this subroutine is used to calculate the R_func in the notes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use channels
        use mesh
        use scatt
        use systems
        use constants
        use input,only:prior,printf,ut
        implicit none
        integer :: alphabar
        integer :: alpha2b,alphaspect,alpha_out,alpha_in
        real*8 :: rx,rb,rxp
        integer :: irx,irb,irxp
        complex*16 :: lambda,lambdaNO
        complex*16,dimension(1:nr) :: R_func,R_funcNO
        complex*16,dimension(1:nr) :: rho
        integer :: lb, la 
        real*8 :: N
        real*8,dimension(1:nx,1:nr,1:nr) :: Gabar
        real*8,dimension(1:nx) :: G_alphabar
        ! complex*16,dimension(1:nx) :: func_in,vpot_in


        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)

        call G_alphaout_alphain(alphabar,Gabar)

      
        
        R_func=0.0_dpreal
        R_funcNO=0.0_dpreal

        la=in3b%lam(alpha_in)
        lb=out3b%lam(alpha_out)
        alphaspect=out3b%alphaspect(alpha_out)
        alpha2b=out3b%alpha2b(alpha_out)

        rho=0.0_dpreal
        do irxp=1,nr
          rxp=rr(irxp)

          do irb=1,nr
            rb=rr(irb)
            G_alphabar(1:nx)=Gabar(1:nx,irb,irxp)

            if(prior) then

               call lambda_func_prior(alphabar,G_alphabar,irb,irxp,lambda,lambdaNO)

              rho(irxp)=rho(irxp) + rb * iu**(-lb) * wf_b(irb,alphaspect) * lambda * rrw(irb)
              R_funcNO(irxp)=R_funcNO(irxp) + rb * iu**(-lb) * wf_b(irb,alphaspect) * lambdaNO * rrw(irb)
            else
              call lambda_func(alphabar,G_alphabar,irb,irxp,lambda)
              rho(irxp)=rho(irxp) + rb * iu**(-lb) * wf_b(irb,alphaspect) * lambda * rrw(irb)
            end if
          end do

          if (printf) write(11,*) rxp,real(rho(irxp)),aimag(rho(irxp))
          if (printf) write(13,*) rxp,8.0_dpreal*pi**2*real(R_funcNO(irxp)),aimag(R_funcNO(irxp))

        end do

        if (printf)   write(11,*) "&"
        if (printf)   write(13,*) "&"

        R_func=0.0_dpreal
        do irx=1,nr
          rx=rr(irx)

          do irxp=1,nr
            rxp=rr(irxp)
            N=-2.0_dpreal*mux*rxp* 2.0**7 * pi**4 /hbarc/hbarc/kx/rx/ka/kb
            R_func(irx)= R_func(irx) + N * Gx%re(min(irx,irxp),alpha2b) * Gx%ir(max(irx,irxp),alpha2b) * rho(irxp) *rrw(irxp)
          end do

          if(prior .and. (.not. UT)) R_func(irx)=R_func(irx)+R_funcNO(irx)* 2.0**7 * pi**4/ka/kb
          if (printf)   write(12,*) rx,real(R_func(irx)),aimag(R_func(irx))
        end do
        if (printf)   write(12,*)"&"
        Rla(:,la) =Rla(:,la)+ abs(R_func)**2 
       end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_lambda_func(nchin)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh
      implicit none
      integer :: irb,irx,nchin

        do irx=1,nr
          do irb=1,nr
            call initial_lambda_func_nxloop(irb,irx,nchin)
          end do
        end do
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine initial_lambda_func_nxloop(irb,irx,alpha_in)
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
        integer :: irb,irx
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
        rb=rr(irb)

        if (prior) then

          nchbx=in3b%alpha2b(alpha_in)
          nchaA=in3b%alphaspect(alpha_in)

          do ix=1,nx

              x=angx(ix)
              ra=sqrt((1.-p*q)**2 * rx**2 + p**2 * rb**2 + 2.*(1.-p*q)*p*rx*rb*x)
              rbx=sqrt(q**2 * rx**2 + rb**2 - 2.*q*rx*rb*x)
              rbA=sqrt((1.-q)**2 * rx**2 + rb**2 + 2*(1.-q)*rx*rb*x)

              wfa=FFC(ra/hcm,wf_a(0:irmatch,ainspect),irmatch+1)*iu**la / ra
              wfbx=FFR4(rbx/hcm,phi_bx(0:irmatch,nchbx),irmatch+1) / rbx
              func_in=wfa*wfbx
              fun_in(ix,irb,irx,alpha_in)=func_in
              vpot_in=FFC(rbA/hcm,UbA,irmatch+1)-FFC(ra/hcm,UaA(0:irmatch,nchaA),irmatch+1)
              pot_in(ix,irb,irx,alpha_in)=vpot_in

          end do

        else

          nchbx=in3b%alpha2b(alpha_in)
          do ix=1,nx
            x=angx(ix)
            ra=sqrt((1.-p*q)**2 * rx**2 + p**2 * rb**2 + 2.*(1.-p*q)*p*rx*rb*x)
            rbx=sqrt(q**2 * rx**2 + rb**2 - 2.*q*rx*rb*x)
            rbA=sqrt((1.-q)**2 * rx**2 + rb**2 + 2*(1.-q)*rx*rb*x)


            wfa=FFC(ra/hcm,wf_a(0:irmatch,ainspect),irmatch+1) * iu**la / ra
            wfbx=FFR4(rbx/hcm,phi_bx(0:irmatch,nchbx),irmatch+1) / rbx
            func_in=wfa*wfbx
            fun_in(ix,irb,irx,alpha_in)=func_in

            vpot_in=FFR4(rbx/hcm,vbx(0:irmatch,nchbx),irmatch+1) + FFC(rbA/hcm,UbA,irmatch+1)
            pot_in(ix,irb,irx,alpha_in)=vpot_in
          end do

        end if

       end subroutine
c-----------------------------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine lambda_func(alphabar,G_alphabar,irb,irx,lambda)
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
        integer :: alphabar
        integer :: alpha_out,irb,irx,alpha_in
        complex*16 :: lambda
        integer :: ix
        complex*16 :: Vpost
        real*8 :: rx,rb
        integer :: nchbB ! for UbB potential
        integer :: LL,lx
        complex*16 :: UUbB
        complex*16,dimension(1:nx) :: func_in,vpot_in
        real*8,dimension(1:nx) :: G_alphabar


        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)

        func_in(1:nx)=fun_in(1:nx,irb,irx,alpha_in)
        vpot_in(1:nx)=pot_in(1:nx,irb,irx,alpha_in)


        lambda=0.0_dpreal
        LL=nint(out3b%j(alpha_out))
        lx=out3b%l(alpha_out)

        rx=rr(irx)
        rb=rr(irb)

        nchbB=out3b%alphaspect(alpha_out)
        UUbB=FFC(rb/hcm,UbB(0:irmatch,nchbB),1+irmatch)


        do ix=1,nx


          vpost=vpot_in(ix) - UUbB


          lambda = lambda + G_alphabar(ix)  * vpost * angw(ix) * func_in(ix)


        end do
          lambda=lambda/(2.0_dpreal*LL+1.0_dpreal)
       end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine lambda_func_prior(alphabar,G_alphabar,irb,irx,lambda,lambdaNO)
C       this subroutine is used to calculate the DELTA function that
C       appears in the notes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use pot
        use scatt
        use bound
        use channels
        use interpolation
        use mesh
        use input, only:ubx
        implicit none
        integer :: alphabar
        integer :: alpha_out,irb,irx,alpha_in
        complex*16 :: lambda,lambdaNO
        integer :: ix
        real*8 :: rx,rb,rbx,q
        integer :: nchxA ! for UxA potential
        integer :: LL,lx
        complex*16 :: UUx
        complex*16,dimension(1:nx) :: func_in,vpot_in
        complex*16 :: func
        real*8,dimension(1:nx) :: G_alphabar
        complex*16 :: Wbx


        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)

        func_in(1:nx)=fun_in(1:nx,irb,irx,alpha_in)
        vpot_in(1:nx)=pot_in(1:nx,irb,irx,alpha_in)

        lambda=0.0_dpreal
        lambdaNO=0.0_dpreal
        LL=nint(out3b%j(alpha_out))
        lx=out3b%l(alpha_out)



        rx=rr(irx)
        rb=rr(irb)

        nchxA=out3b%alpha2b(alpha_out)

        UUx=FFC(rx/hcm,UxA(0:irmatch,nchxA),irmatch+1)

        q=masst/(masst+massx)

        do ix=1,nx
        
            if(ubx) then 
            rbx=sqrt(q**2 * rx**2 + rb**2 - 2.*q*rx*rb*angx(ix))
            Wbx=FFR4(rbx/hcm,vbx(0:irmatch,1),irmatch+1)*iu
            end if 

            func=G_alphabar(ix)  * angw(ix) * func_in(ix)
            if (ubx) then 
            lambda = lambda + func  * (vpot_in(ix)+wbx + UUx)
            else
            lambda = lambda + func  * (vpot_in(ix) + UUx)
            end if 
            lambdaNO = lambdaNO + func

        end do
        lambda=lambda/(2.0_dpreal*LL+1.0_dpreal)
        lambdaNO=lambdaNO/(2.0_dpreal*LL+1.0_dpreal)

       end subroutine
C-----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_in(ix)
c
c     vec{rbx}=q*rx-rb     vec{ra}=(1-pq)rx+prb
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use spharm
      use precision
      use mesh
      use channels
      implicit none
      real*8 :: p,q
      integer :: ix
      integer :: irx,irb,nch
      real*8 :: rbx_z,rbx_abs,ra_z,ra_abs
      real*8 :: rb_z,rx_z,rx,rb,x
      integer ::  mlbx,mla
      integer :: l2b,la
      real*8,dimension(0:lmax,0:lmax) :: PL
      real*8,dimension(0:lbx,0:lbx) :: plbx

c-----define coefficients
      p=massb/(massb+massx)
      q=masst/(masst+massx)



      do irx=1, nr
        rx=rr(irx)
        do irb=1,nr
          rb=rr(irb)
C          do ix=1,nx
            x=angx(ix)
            rb_z=rb
            rx_z=rx*x

            ra_z=(1-p*q)*rx_z+p*rb_z
            ra_abs=sqrt((1-p*q)**2 * rx**2 + p**2 * rb**2 + 2*(1-p*q)*p*rx*rb*x)
            call PLM(ra_z/ra_abs,lmax,lmax,lmax+1,PL)
            do la=0,lmax
               do mla=-la,la

                 nch=la**2+la+mla+1
                 Y_in_la(ix,irb,irx,nch)=YLMC(la,mla)*PL(la,abs(mla)) ! phi angle is 0 for ra since ra_x is positive

               end do
            end do

            rbx_z=q*rx_z-rb_z
            rbx_abs=sqrt(q**2 * rx**2 + rb**2 - 2*q*rx*rb*x)
            call PLM(rbx_z/rbx_abs,lbx,lbx,lbx+1,plbx)
C           do l2b=lbxmin,lbxmax
              l2b=lbx
              do mlbx=-l2b,l2b

                nch=l2b**2+l2b+mlbx+1
                Y_in_lbx(ix,irb,irx,nch)=YLMC(l2b,mlbx)*plbx(l2b,abs(mlbx)) ! phi angle is 0 since rbx_x is positive
              end do
C           end do

C          end do !ix

        end do ! irb
      end do  ! irx

      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_Y_out(ix)
c             | 0 |             | rx*sqrt(1-x^2) |
c     vec{rb}=| 0 |     vec{rx}=|        0       |
c             | rb|             |       rx*x     |
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use spharm
      use precision
      use mesh
      use channels
      implicit none
      integer :: ix
      integer :: l,m,nch
      real*8,dimension(0:min(lmax,cutl),0:min(lmax,cutl)) :: PL
!--


C      do ix=1,nx
         call PLM(angx(ix),min(lmax,cutl),min(lmax,cutl),min(lmax,cutl)+1,PL)

         do l=0,min(lmax,cutl)
           do m=-l,l

              nch=l**2+l+m+1
              Y_out_lx(ix,nch)=YLMC(l,-m)*PL(l,abs(m))*(-1.0_dpreal)**(m)
           end do
         end do

C      end do

      end subroutine
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine G_alphaout_alphain(alphabar,Gabar)
C     this subroutine is used to calculate the G_alphabar(ix,irb,irx)
c     coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      use cleb_coefficient
      implicit none
      integer :: irb,irx,ix
      integer :: alphabar,alpha_in, alpha_out
      integer :: ML,ma,mbx
      integer :: nch_L,nch_la,nch_lbx,nch_lx,nch_lb
      integer :: LL,la,l2b,lb,lx,minl
      real*8 ::j_out
      real*8 :: CGin,CGout,Yin
      real*8,dimension(1:nx) :: Yout
      real*8,dimension(1:nx,1:nr,1:nr),target :: Gabar
      real*8,pointer,dimension(:) :: yla,ylbx,ylx,GG
      real*8 :: func



        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)
        j_out=out3b%j(alpha_out)
        LL=nint(j_out)
        la=in3b%lam(alpha_in)
        l2b=in3b%l(alpha_in)
        lx=out3b%l(alpha_out)
        lb=out3b%lam(alpha_out)
        minl=min(LL,lx)
        Gabar=0.0_dpreal



        do ML=-minl,minl
          nch_lx=lx**2+lx+ML+1
          nch_lb=lb**2+lb+0+1
          nch_L=LL**2+LL+ML+1
          CGout=cleb(lx*2,ML*2,lb*2,0,LL*2,ML*2)
          ylx=>Y_out_lx(:,nch_lx)
          Yout(:)=ylx(:)*sqrt( (2.*lb+1.)/(4.0_dpreal*pi) )
          do ma=-la,la
            do mbx=-l2b,l2b
              if (ma+mbx/=ML) cycle
              nch_lbx=lbx**2+lbx+mbx+1
              nch_la=la**2+la+ma+1
              CGin=cleb(lbx*2,mbx*2,la*2,ma*2,LL*2,ML*2)

              do irx=1,nr
                do irb=1,nr
                  yla=>Y_in_la(:,irb,irx,nch_la)
                  ylbx=>Y_in_lbx(:,irb,irx,nch_lbx)
                  GG=>Gabar(:,irb,irx)
         
                  do ix=1,nx

                    Yin=yla(ix)*ylbx(ix)

                    func=CGout*CGin*Yout(ix)*Yin

                    GG(ix)=GG(ix)+func


                  end do

                end do
              end do
            end do
          end do

        end do

        ! Gabar = Gabar *  8.0_dpreal*pi**2 






      end subroutine
c-----------------------------------------------------------------------




      end module
