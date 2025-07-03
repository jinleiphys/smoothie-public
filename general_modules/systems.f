ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module systems
c     system: identify the reaction systems
c           a + A -> b + x + A -> x + b + A -> b+anything
c           =====        =====        =====
c             a            x            b
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 :: zp, massp !define the charge and mass of a
       real*8 :: zt, masst !define the charge and mass of A
       real*8 :: zb, massb !define the charge and mass of b
       real*8 :: zx, massx !define the charge and mass of x
       real*8 :: qval,be !  for most case be=qval, but for 6He with di-neutron model ,be/=qval
       integer :: nodes
       character(len=5) :: namep,namet,nameb,namex
       character(len=40):: gswf ! added by AMoro
       real*8 :: jp,jt,jb,jx !spin of a,A,b,x
       integer :: neb
       integer :: ptyp, ptyt,ptyb,ptyx ! define the particle spin  
    

       real*8,allocatable,dimension(:) :: necmb   ! center of mass energy of outgoing b
       integer,allocatable,dimension(:) :: nkn
       integer,allocatable,dimension(:) :: nbin ! bin=0 single point, bin>0 b-bins method, bin is the number of Gauss points
       real*8,allocatable,dimension(:) :: nwbin  ! width of b-bin in the momentum space,kmax=k+wbin/2, kmin=k-wbin/2

       real*8 :: elab !energy of the reaction system

      end module systems
c----------------------------------------------------------------------
