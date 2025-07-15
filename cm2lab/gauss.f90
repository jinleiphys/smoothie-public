      module gaussm3

                        ! Common Global variables within module !
!  	implicit none
       INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!   PRIVATE
!   REAL (dbp) :: newv
      REAL(dbp)  :: EPS, M_PI
      PARAMETER (EPS=3.0d-15)             !EPS is the relative precision
      PARAMETER (M_PI=3.141592654d0)      ! Pi value

!   PUBLIC :: newv, EPS, M_PI, n, xabsc, weig, dbp, qgss2d

      end module gaussm3
!   INTERFACE            

!  END INTERFACE

!   CONTAINS
!* This module has the following INTERNAL FUNCTIONS:
!* gauleg, qgauss, qgss3d, qgss2d, gsselm, identity_matrix
!* This module has the following INTERNAL SUBROUTINES:
!* linear_solver
!* They can call on each other without first specifying their type
!* NO INTERFACE nor EXTERNAL is required since they are INTERNAL functions
!Module contains:
!gauleg - Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature integration

!qgauss - N-point Gauss-Legendre single integration
!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
      SUBROUTINE  gauleg(ngp, xabsc, weig)
      use gaussm3
      implicit none
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


       m = (ngp + 1) / 2

!* Roots are symmetric in the interval - so only need to find half of them  */

       do i = 1, m            ! Loop over the desired roots */

            z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!     		write(*,*) 'pi=',M_PI
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */
            
        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                         ! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                 ! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop

      End subroutine gauleg



       ! subroutine to generate the grid points and width for Simpson integration
       subroutine simpson(nnu,xstart,xmax,rr,rw)
       implicit none 
       integer :: nxmx,nnu
       real*8 :: xmax,xstart
       real*8,dimension(1:nnu) :: rr,rw
       real*8,dimension(1:nnu+1) :: sxx,wxx
       real*8 :: dx, d43,d23
       integer :: nxmx1, nxmx2,i 
       
        nxmx=nnu+1
        dx=(xmax-xstart)/float(nxmx-1)
 
        d43=4.d0/3.d0
        d23=2.d0/3.d0
   
        wxx(1)=dx/3.d0
        wxx(nxmx)=dx/3.d0
        sxx(1)=xstart
        sxx(nxmx)=float(nxmx-1)*dx+xstart

 
        nxmx1=nxmx-1
        nxmx2=nxmx-2
 
        do 50 i=2,nxmx1,2
        sxx(i)=float(i-1)*dx+xstart
  50    wxx(i)=d43*dx
 
        do 55 i=3,nxmx2,2
        sxx(i)=float(i-1)*dx+xstart
  55    wxx(i)=d23*dx
  
  
        do i=1,nnu
          rr(i)=sxx(i+1)
          rw(i)=wxx(i+1)
        end do 
       end subroutine  

        


 
