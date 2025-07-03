      module mesh
      integer :: nr,nr1,nr2,nk
      real*8 :: r1,r2,rmax
c      wave functions calculated at intervals of HCM up to abs(RMATCH).
      real*8 :: rmatch
      real*8 :: hcm
      integer :: irmatch,irmax
      integer :: nx !  number of Gauss Points
      real*8  :: thmin,thmax,thinc
      integer :: nthmax  !nthmax=nint((thmax-thmin)/thinc)
      real*8,allocatable,dimension(:) :: rr, rrw
ccc use for bin method
      integer :: irbmatch
      real*8  :: rbmatch
ccc  mesh set for CDCC-NEB
      real*8,dimension(1:2) :: rmaxcdcc,rmatchcdcc
      integer,dimension(1:2) :: nrcdcc
      real*8,allocatable,dimension(:) :: rr_in, rrw_in
      real*8,allocatable,dimension(:) :: rr_out, rrw_out
ccc  mesh set for CDCC calculations 
      real*8 :: rxmax, rymax 
      integer :: nrx, nry 
ccc  mesh set for Vincent-Fortune calculations 
      integer :: Pi_ext_ny ! grid number for y direction
      real*8 :: Pi_ext_ymax! maximum value of y 
      real*8,allocatable :: yy_mesh(:),yyw(:) ! grids for y direction
      end module
