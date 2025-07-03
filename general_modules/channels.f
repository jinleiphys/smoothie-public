      module channels
      use precision
      use systems

       integer :: lmax !maximum value of l
       integer :: lmin !minimum value of l
       real*8 :: sbx ! coupling results of jb and jx
       integer :: cutl ! cut lx
       integer :: lxmax,lxmin  ! maximum value and minimum value of lx
       integer :: lbx,lbxmin,lbxmax
       real*8 :: jtmax,jtmin
       integer :: lmax_cdcc
       integer :: abarMLmax,alpha_in_MLmax
       integer,dimension(:,:),allocatable :: abarML,alpha_in_ML
! coupling coefficients of alpha_in and alpha_out
       Type alpha_bar
       integer,dimension(:),allocatable :: alpha_in,alpha_out
       integer,dimension(:,:),allocatable :: alpha_cdcc
       integer,dimension(:),allocatable :: cdccnchmax
       End type
       Type(alpha_bar) :: abar




!Jacobi set for 3b channels
       Type nch_3b
!  | (l s_2b)j_2b, (lam j? ) j_spect; J >
        integer :: nchmax
        integer,dimension(:),allocatable :: alpha2b
        integer,dimension(:),allocatable :: alphaspect
        integer,dimension(:),allocatable :: l,lam
        real*8,dimension(:),allocatable :: s2b,j2b
        real*8,dimension(:),allocatable  ::j,j_spect
       End type



! subsystem 2b channels
       Type nch_2b
!  |(ls)j>
         integer :: nchmax
         integer,dimension(:),allocatable :: l
         real*8,dimension(:),allocatable  :: s, j
       End type



!index for spectator
       Type nch_spect
! |(lam s) j_spect>
         integer :: nchmax
         integer,dimension(:),allocatable :: lam
         real*8,dimension(:),allocatable :: j
       end type

!Jacobi set for 3b channels
       Type nch_cdcc
!  | (l s_2b)j_2b, (lam j? ) j_spect; J >
        integer :: nchmax
        integer,dimension(:),allocatable :: l,lam,n
        real*8,dimension(:),allocatable :: s2b,j2b
        real*8,dimension(:),allocatable  ::j,j_spect
       End type




       type(nch_3b) :: in3b,out3b
       type(nch_2b) :: in2b,out2b
       type(nch_spect) :: inspect,outspect
       type(nch_cdcc) :: incdcc


       contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_2b_in()
c     index of subsystem of b and x
c     index of {(l_ij (ji jj)sij) Jij}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l
      real*8 :: s,J
      integer :: ns, nJ
      integer :: nch  ! channel index {(l (jb jx)s) J}

      in2b%nchmax=0
      do l=lbxmin,lbxmax
         do ns=nint(2.*abs(jb-jx)),nint(2.*(jb+jx)),2
           if(ns/=nint(2.*sbx)) cycle
            s=ns/2.0_dpreal
            do nJ=nint(2.*abs(l-s)),nint(2.*(l+s)),2
               if (nJ /= nint(2.*jp)) cycle
               J=nJ/2.0_dpreal
               in2b%nchmax=in2b%nchmax+1
            end do
         end do
      end do

       allocate(in2b%l(1:in2b%nchmax),in2b%s(1:in2b%nchmax))
       allocate(in2b%j(1:in2b%nchmax))
      !  write(*,*) "in2b%nchmax=",in2b%nchmax
       if(in2b%nchmax/=1) stop "error in incoming projectile LS couplings"
      nch=1
      do l=lbxmin,lbxmax
         do ns=nint(2.*abs(jb-jx)),nint(2.*(jb+jx)),2
            if(ns/=nint(2.*sbx)) cycle
            s=ns/2.0_dpreal
            do nJ=nint(2.*abs(l-s)),nint(2.*(l+s)),2
               if (nJ /= nint(2.*jp)) cycle
               J=nJ/2.0_dpreal
               in2b%l(nch)=l
               in2b%s(nch)=s
               in2b%j(nch)=J
               nch=nch+1
            end do
         end do
      end do
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_spect_in()
c     index of spectator
c     index of {(lam jt) J_spect}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: lam
      real*8 :: s,J
      integer :: nJ
      integer :: nch  ! channel index {(l jt) J}


      inspect%nchmax=0
      do lam=lmin,lmax
         s=jt
         do nJ=nint(2*abs(lam-s)),nint(2*(lam+s)) ,2
            J=nJ/2.0_dpreal
            inspect%nchmax=inspect%nchmax+1
         end do
      end do

       allocate(inspect%lam(1:inspect%nchmax))
       allocate(inspect%j(1:inspect%nchmax))

      nch=1
      do lam=lmin,lmax
        s=jt
        do nJ=nint(2*abs(lam-s)),nint(2*(lam+s)) ,2
           J=nJ/2.0_dpreal
           inspect%lam(nch)=lam
           inspect%j(nch)=J
           nch=nch+1
        end do
      end do

      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_3b_in()
c     index of  | \alpha >_{in}
c     index of {(alpha_2b (lam jt)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J
      integer :: nJ
      integer :: nch !index of {(l(jxjb)s2b)J2b (lam jt)J_spect ; J}
      integer :: nch2b,nchspect

      call alpha_2b_in()
      call alpha_spect_in()

      in3b%nchmax=0
      do nch2b=1,in2b%nchmax
         J2b=in2b%J(nch2b)
         l=in2b%l(nch2b)
         do nchspect=1,inspect%nchmax
            lambda=inspect%lam(nchspect)
            j_spect=inspect%j(nchspect)
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
C              j=nj/2.0_dpreal
      	      in3b%nchmax=in3b%nchmax+1
      	    end do
         end do
      end do
      write(8,10)in3b%nchmax
10    format('there are',I3,1X,'incoming channels')

       if (allocated(in3b%l)) deallocate(in3b%l)
       if (allocated(in3b%lam)) deallocate(in3b%lam)
       if (allocated(in3b%J)) deallocate(in3b%J)
       if (allocated(in3b%s2b)) deallocate(in3b%s2b)
       if (allocated(in3b%J2b)) deallocate(in3b%J2b)
       if (allocated(in3b%j_spect)) deallocate(in3b%j_spect)
       if (allocated(in3b%alpha2b)) deallocate(in3b%alpha2b)
       if (allocated(in3b%alphaspect)) deallocate(in3b%alphaspect)


       allocate(in3b%l(1:in3b%nchmax),in3b%lam(1:in3b%nchmax))
       allocate(in3b%J(1:in3b%nchmax),in3b%s2b(1:in3b%nchmax))
       allocate(in3b%J2b(1:in3b%nchmax),in3b%j_spect(1:in3b%nchmax))
       allocate(in3b%alpha2b(1:in3b%nchmax))
       allocate(in3b%alphaspect(1:in3b%nchmax))
      nch=1
      write(8,20)
      write(8,30)
      do nch2b=1,in2b%nchmax
         l=in2b%l(nch2b)
         s2b=in2b%s(nch2b)
         J2b=in2b%J(nch2b)
         do nchspect=1,inspect%nchmax
           lambda=inspect%lam(nchspect)
           j_spect=inspect%j(nchspect)
           jmin=max(jtmin,abs(J2b-J_spect))
           jmax=min(jtmax,J2b+J_spect)
           do nj=nint(2.*jmin),nint(2.*jmax),2
               j=nj/2.0_dpreal
      	       in3b%l(nch)=l
               in3b%lam(nch)=lambda
               in3b%J(nch)=J
               in3b%s2b(nch)=s2b
               in3b%J2b(nch)=J2b
               in3b%j_spect(nch)=j_spect
               in3b%alpha2b(nch)=nch2b
               in3b%alphaspect(nch)=nchspect
               write(8,40)nch,nch2b,l,jb,jx,s2b,J2b,
     +                       lambda,jt,j_spect,j
      	           nch=nch+1
      	    end do
         end do
      end do

20    format('---For incoming channels the coupling coefficients are')
30    format(' a3b','|', 'a2b','|','(',' l ','(',' jb ',' jx ',')',
     +      ' s2b ',')',' J2b ', '(',' lam ',' jt ',')',' J3 ',',',
     +       ' Jtot ')
40    format(I4,I4,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,
     +       f4.1,3x,I3,1x,f3.1,1x,f4.1,3x,f4.1)
      end subroutine
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c now is the set for outgoing channel index
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_2b_out()
c     index of subsystem of A and x
c     index of {(l_ij (ji jj)sij) Jij}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l
      real*8 :: s,J
      integer :: ns, nJ
      integer :: nch  ! channel index {(l (jt jx)s) J}

      out2b%nchmax=0
      do l=lxmin,min(lmax,lxmax)
         do ns=nint(2.*abs(jx-jt)),nint(2.*(jt+jx)),2
            s=ns/2.0_dpreal
            do nJ=nint(2.*abs(l-s)),nint(2.*(l+s)),2
               J=nJ/2.0_dpreal
               out2b%nchmax=out2b%nchmax+1
            end do
         end do
      end do

       allocate(out2b%l(1:out2b%nchmax),out2b%s(1:out2b%nchmax))
       allocate(out2b%j(1:out2b%nchmax))

      nch=1
      do l=lxmin,min(lmax,lxmax)
         do ns=nint(2.*abs(jx-jt)),nint(2.*(jt+jx)),2
            s=ns/2.0_dpreal
            do nJ=nint(2.*abs(l-s)),nint(2.*(l+s)),2
               J=nJ/2.0_dpreal
               out2b%l(nch)=l
               out2b%s(nch)=s
               out2b%j(nch)=J
               nch=nch+1
            end do
         end do
      end do
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_spect_out()
c     index of spectator
c     index of {(lam jb) J_spect}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: lam
      real*8 :: s,J
      integer ::  nJ
      integer :: nch  ! channel index {(l jt) J}


      outspect%nchmax=0
      do lam=lmin,lmax
         s=jb
         do nJ=nint(2*abs(lam-s)),nint(2*(lam+s)) ,2
            J=nJ/2.0_dpreal
            outspect%nchmax=outspect%nchmax+1
         end do
      end do

       allocate(outspect%lam(1:outspect%nchmax))
       allocate(outspect%j(1:outspect%nchmax))

      nch=1
      do lam=lmin,lmax
        s=jb
        do nJ=nint(2*abs(lam-s)),nint(2*(lam+s)) ,2
           J=nJ/2.0_dpreal
           outspect%lam(nch)=lam
           outspect%j(nch)=J
           nch=nch+1
        end do
      end do

      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_3b_out()
c     index of  | \alpha >_{out}
c     index of {(alpha_2b (lam jb)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J
      integer :: nJ
      integer :: nch !index of {(l(jxjt)s2b)J2b (lam jb)J_spect ; J}
      integer :: nch2b,nchspect

      call alpha_2b_out()
      call alpha_spect_out()

      out3b%nchmax=0
      do nch2b=1,out2b%nchmax
         J2b=out2b%J(nch2b)
         l=out2b%l(nch2b)
         do nchspect=1,outspect%nchmax
            lambda=outspect%lam(nchspect)
            j_spect=outspect%j(nchspect)
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2

C              j=nj/2.0_dpreal
      	      out3b%nchmax=out3b%nchmax+1
      	    end do
         end do
      end do
      write(8,10)out3b%nchmax
10    format('there are',I6,1X,'outgoing channels')



      if (allocated(out3b%l)) deallocate(out3b%l)
      if (allocated(out3b%lam)) deallocate(out3b%lam)
      if (allocated(out3b%J)) deallocate(out3b%J)
      if (allocated(out3b%s2b)) deallocate(out3b%s2b)
      if (allocated(out3b%J2b)) deallocate(out3b%J2b)
      if (allocated(out3b%j_spect)) deallocate(out3b%j_spect)
      if (allocated(out3b%alpha2b)) deallocate(out3b%alpha2b)
      if (allocated(out3b%alphaspect)) deallocate(out3b%alphaspect)


       allocate(out3b%l(1:out3b%nchmax),out3b%lam(1:out3b%nchmax))
       allocate(out3b%J(1:out3b%nchmax),out3b%s2b(1:out3b%nchmax))
       allocate(out3b%J2b(1:out3b%nchmax),out3b%j_spect(1:out3b%nchmax))
       allocate(out3b%alpha2b(1:out3b%nchmax))
       allocate(out3b%alphaspect(1:out3b%nchmax))
      nch=1
      write(8,20)
      write(8,30)
      do nch2b=1,out2b%nchmax
         l=out2b%l(nch2b)
         s2b=out2b%s(nch2b)
         J2b=out2b%J(nch2b)
         do nchspect=1,outspect%nchmax
           lambda=outspect%lam(nchspect)
           j_spect=outspect%j(nchspect)
           jmin=max(jtmin,abs(J2b-J_spect))
           jmax=min(jtmax,J2b+J_spect)
           do nj=nint(2.*jmin),nint(2.*jmax),2

               j=nj/2.0_dpreal
      	       out3b%l(nch)=l
               out3b%lam(nch)=lambda
               out3b%J(nch)=J
               out3b%s2b(nch)=s2b
               out3b%J2b(nch)=J2b
               out3b%j_spect(nch)=j_spect
               out3b%alpha2b(nch)=nch2b
               out3b%alphaspect(nch)=nchspect
               write(8,40)nch,nch2b,l,jt,jx,s2b,J2b,
     +                       lambda,jb,j_spect,j
      	           nch=nch+1
      	    end do
         end do
      end do

20    format('---For outgoing channels the coupling coefficients are')
30    format(' a3b ','|', ' a2b','|','(',' l ','(',' jt ',' jx ',')',
     +      ' s2b ',')',' J2b ', '(',' lam ',' jb ',')',' J3 ',',',
     +       ' Jtot ')
40    format(I5,I5,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,
     +       f4.1,3x,I3,1x,f3.1,1x,f4.1,3x,f4.1)
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine abar_index(alpha_xA,abarmax)
c     this file used to compute the abar_index
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: la,lb,lx,l2b
      real*8 :: J_in, J_out
      integer :: alpha_in, alpha_out,alpha_xA
      integer :: nch,abarmax


      abarmax=0
      do alpha_in=1, in3b%nchmax
        l2b=in3b%l(alpha_in)
        la=in3b%lam(alpha_in)
        J_in=in3b%j(alpha_in)
        do alpha_out=1,out3b%nchmax
           lx=out3b%l(alpha_out)
           lb=out3b%lam(alpha_out)
           J_out=out3b%j(alpha_out)
           if(out3b%alpha2b(alpha_out)/= alpha_xA) cycle
           if(nint(2.*J_in)/=nint(2.*J_out)) cycle
           if((-1)**(l2b+la) /= (-1)**(lb+lx)) cycle

           abarmax=abarmax+1
        end do
      end do
      if (abarmax/=0) then
      write(*,10) alpha_xA
      write(*,11)abarmax
      end if
10    format("For alpha_xA=", I4)
11    format("there are ",I8,"  channel couplings to be solved")
      if (allocated(abar%alpha_in)) deallocate(abar%alpha_in)
      if(allocated(abar%alpha_out)) deallocate(abar%alpha_out)
      allocate(abar%alpha_in(1:abarmax))
      allocate(abar%alpha_out(1:abarmax))


      nch=0
      do alpha_in=1, in3b%nchmax
        l2b=in3b%l(alpha_in)
        la=in3b%lam(alpha_in)
        J_in=in3b%j(alpha_in)

        do alpha_out=1,out3b%nchmax
           lx=out3b%l(alpha_out)
           lb=out3b%lam(alpha_out)
           J_out=out3b%j(alpha_out)
           if(out3b%alpha2b(alpha_out)/= alpha_xA) cycle
           if(nint(2.*J_in)/=nint(2.*J_out)) cycle
           if((-1)**(l2b+la) /= (-1)**(lb+lx)) cycle



           nch=nch+1
           abar%alpha_in(nch)=alpha_in
           abar%alpha_out(nch)=alpha_out
        end do

      end do



      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_cdcc_in()
c     index of  | \alpha >_{in}
c     index of {(alpha_2b (lam jt)j_spect ; J}
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l,lambda
      real*8 :: jmin,jmax
      real*8 :: s2b,j2b,j_spect,J,s
      integer :: nJ,nJ2b,nJ_spect
      integer :: nch !index of {(l(jxjb)s2b)J2b (lam jt)J_spect ; J}
      real*8 :: jtmin, jtmax
      integer,dimension(10) :: l2b
      integer :: nl
      namelist /cdccwf/ l2b,lmax_cdcc,jtmin, jtmax
      l2b=-99
      open (unit=551,file='iav.in')
      read(551,nml=cdccwf)
      s2b=0.0d0 ! by assuming all spin zero particle, need to improved
      s=0.0d0 ! spin of the target


      incdcc%nchmax=0
C     do l=0,l2bmax
      do nl=1,10
        l=l2b(nl)
        if(l<0) cycle
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
      	      incdcc%nchmax=incdcc%nchmax+1
      	    end do
           end do
         end do
        end do
      end do
      write(8,10)incdcc%nchmax
10    format('there are',I3,1X,'cdcc channels')

       if (allocated(incdcc%l)) deallocate(incdcc%l)
       if (allocated(incdcc%lam)) deallocate(incdcc%lam)
       if (allocated(incdcc%J)) deallocate(incdcc%J)
       if (allocated(incdcc%s2b)) deallocate(incdcc%s2b)
       if (allocated(incdcc%J2b)) deallocate(incdcc%J2b)
       if (allocated(incdcc%j_spect)) deallocate(incdcc%j_spect)
       if (allocated(incdcc%n)) deallocate(incdcc%n)



       allocate(incdcc%l(1:incdcc%nchmax),incdcc%lam(1:incdcc%nchmax))
       allocate(incdcc%J(1:incdcc%nchmax),incdcc%s2b(1:incdcc%nchmax))
       allocate(incdcc%J2b(1:incdcc%nchmax))
       allocate(incdcc%j_spect(1:incdcc%nchmax))
       allocate(incdcc%n(1:incdcc%nchmax))

       incdcc%n=0

      nch=1
      write(8,20)
      write(8,30)
C     do l=0,l2bmax
      do nl=1,10
        l=l2b(nl)
        if(l<0) cycle
        do nJ2b=nint(2.*(l+s2b)),nint(2.*abs(l-s2b)),2
           j2b=nJ2b/2.0d0
         do lambda=0,lmax_cdcc
           do nJ_spect=nint(2*abs(lambda-s)),nint(2*(lambda+s)),2
            J_spect=nJ_spect/2.0d0
            jmin=max(jtmin,abs(J2b-J_spect))
            jmax=min(jtmax,J2b+J_spect)
      	    do nj=nint(2.*jmin),nint(2.*jmax),2
               j=nj/2.0d0
      	       incdcc%l(nch)=l
               incdcc%lam(nch)=lambda
               incdcc%J(nch)=J
               incdcc%s2b(nch)=s2b
               incdcc%J2b(nch)=J2b
               incdcc%j_spect(nch)=j_spect

               write(8,40)nch,0,l,s,s,s2b,J2b,
     +                       lambda,s,j_spect,j
      	           nch=nch+1
      	    end do
         end do
      end do
      end do
      end do

20    format('---For CDCC channels the coupling coefficients are')
30    format(' a3b','|', 'a2b','|','(',' l ','(',' jb ',' jx ',')',
     +      ' s2b ',')',' J2b ', '(',' lam ',' jt ',')',' J3 ',',',
     +       ' Jtot ')
40    format(I4,I4,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,
     +       f4.1,3x,I3,1x,f3.1,1x,f4.1,3x,f4.1)
      end subroutine
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine abar_index_cdcc(alpha_xA,abarmax)
c     this file used to compute the abar_index
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: la,lb,lx,l2b,l2b_cdcc,la_cdcc
      real*8 :: J_in, J_out,J_CDCC
      integer :: alpha_in, alpha_out,alpha_xA,acdcc
      integer :: nch,abarmax,nchcdccmax,nchcdcc


      abarmax=0
      nchcdccmax=0
      do alpha_in=1, in3b%nchmax
        l2b=in3b%l(alpha_in)
        la=in3b%lam(alpha_in)
        J_in=in3b%j(alpha_in)
        do alpha_out=1,out3b%nchmax
           lx=out3b%l(alpha_out)
           lb=out3b%lam(alpha_out)
           J_out=out3b%j(alpha_out)
           if(out3b%alpha2b(alpha_out)/= alpha_xA) cycle
           if(nint(2.*J_in)/=nint(2.*J_out)) cycle
           if((-1)**(l2b+la) /= (-1)**(lb+lx)) cycle
           abarmax=abarmax+1

           nchcdcc=0
           do acdcc=1, incdcc%nchmax
             l2b_cdcc=incdcc%l(acdcc)
             la_cdcc=incdcc%lam(acdcc)
             J_cdcc=incdcc%j(acdcc)
            if(nint(2.*J_in)/=nint(2.*J_cdcc)) cycle
            if((-1)**(l2b+la) /= (-1)**(l2b_cdcc+la_cdcc)) cycle
            nchcdcc=nchcdcc+1
         end do
         nchcdccmax=max(nchcdcc,nchcdccmax)
        end do
      end do
      if (abarmax/=0) then
      write(*,10) alpha_xA
      write(*,11)abarmax
      end if
10    format("For alpha_xA=", I4)
11    format("there are ",I8,"  couplings channels need to solve")
      if(allocated(abar%alpha_in)) deallocate(abar%alpha_in)
      if(allocated(abar%alpha_out)) deallocate(abar%alpha_out)
      if(allocated(abar%alpha_cdcc)) deallocate(abar%alpha_cdcc)
      if(allocated(abar%cdccnchmax))  deallocate(abar%cdccnchmax)

      allocate(abar%alpha_in(1:abarmax))
      allocate(abar%alpha_out(1:abarmax))
      allocate(abar%alpha_cdcc(1:abarmax,1:nchcdccmax))
      allocate(abar%cdccnchmax(1:abarmax))



      nch=0
      do alpha_in=1, in3b%nchmax
        l2b=in3b%l(alpha_in)
        la=in3b%lam(alpha_in)
        J_in=in3b%j(alpha_in)

        do alpha_out=1,out3b%nchmax
           lx=out3b%l(alpha_out)
           lb=out3b%lam(alpha_out)
           J_out=out3b%j(alpha_out)
           if(out3b%alpha2b(alpha_out)/= alpha_xA) cycle
           if(nint(2.*J_in)/=nint(2.*J_out)) cycle
           if((-1)**(l2b+la) /= (-1)**(lb+lx)) cycle
           nch=nch+1
           abar%alpha_in(nch)=alpha_in
           abar%alpha_out(nch)=alpha_out

           nchcdcc=0
           do acdcc=1, incdcc%nchmax
             l2b_cdcc=incdcc%l(acdcc)
             la_cdcc=incdcc%lam(acdcc)
             J_cdcc=incdcc%j(acdcc)
             if(nint(2.*J_in)/=nint(2.*J_cdcc)) cycle
             if((-1)**(l2b+la) /= (-1)**(l2b_cdcc+la_cdcc)) cycle
             nchcdcc=nchcdcc+1
             abar%alpha_cdcc(nch,nchcdcc)=acdcc
        end do

        abar%cdccnchmax(nch)=nchcdcc
      end do
      end do



      end subroutine
c-----------------------------------------------------------------------


      end module channels
