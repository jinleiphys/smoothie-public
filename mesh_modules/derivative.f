      module derivative
! this module used to store the subroutine for n-order derivative
! based on the Central Finite difference coefficient
! https://en.wikipedia.org/wiki/Finite_difference_coefficient
      use precision
      implicit none

      contains


      subroutine first_derivative_r(f,df,n,h)
!     this is the subroutine used to compute the first derivative of
!     a function f with dimension of n
!     df is the first derivative with dimension of n
!     assuming the f function is real
!     h is the step size of function
!     this subroutine use 4th order accuracy
      implicit none
      integer :: n
      real*8 :: h
      real*8,dimension(1:n) :: df
      real*8,dimension(1:n) :: f
      integer :: i ! loop index

      ! for the first two point, one can use Forward finite difference
      !−25/12	4	−3	4/3	−1/4
      df(1)=(-25.*f(1)/12. + 4.*f(2) - 3.*f(3) + 4.*f(4)/3. -  f(5)/4.)/h
      df(2)=(-25.*f(2)/12. + 4.*f(3) - 3.*f(4) + 4.*f(5)/3. -  f(6)/4.)/h



      ! for the point from 3 to n-2, one can use Central finite difference
      !	1/12	−2/3	0	2/3	−1/12
      do i=3, n-2
        df(i)=  ( f(i-2)/12. - 2.*f(i-1)/3. + 2.*f(i+1)/3. - f(i+2)/12.)/h
      end do

      ! for the lst two point, one can use Backward finite difference, but with 3th order accuracy
      !−1/3	3/2	−3	11/6
      df(n)=(-1.*f(n-3)/3. + 3.*f(n-2)/2. - 3.*f(n-1) + 11.*f(n)/6. )/h
      df(n-1)=(-1.*f(n-4)/3. + 3.*f(n-3)/2. - 3.*f(n-2) + 11.*f(n-1)/6. )/h





      end subroutine






      subroutine second_derivative_r(f,d2f,n,h)
!     this is the subroutine used to compute the second derivative of
!     a function f with dimension of n
!     d2f is the second derivative with dimension of n
!     assuming the f function is real
!     h is the step size of function
!     this subroutine use 4th order accuracy
      implicit none
      integer :: n
      real*8 :: h
      real*8,dimension(1:n) :: d2f
      real*8,dimension(1:n) :: f
      integer :: i ! loop index

      ! for the first two point, one can use Forward finite difference
      !15/4	−77/6	107/6	−13	61/12	−5/6
      d2f(1)=(15.*f(1)/4. - 77.*f(2)/6. + 107.*f(3)/6. - 13.*f(4) +  61.*f(5)/12. - 5.*f(6)/6.)/h**2
      d2f(2)=(15.*f(2)/4. - 77.*f(3)/6. + 107.*f(4)/6. - 13.*f(5) +  61.*f(6)/12. - 5.*f(7)/6.)/h**2



      ! for the point from 3 to n-2, one can use Central finite difference
      !	−1/12	4/3	−5/2	4/3	−1/12
      do i=3, n-2
        d2f(i)=  (- f(i-2)/12. + 4.*f(i-1)/3. - 2.5*f(i) + 4.*f(i+1)/3. - f(i+2)/12.)/h**2
      end do

      ! for the lst two point, one can use Backward finite difference, but with 2th order accuracy
      !−1	4	−5	2
      d2f(n)= (-1.*f(n-3) + 4.*f(n-2) - 5.*f(n-1) + 2.*f(n) )/h**2
      d2f(n-1)= (-1.*f(n-4) + 4.*f(n-3) - 5.*f(n-2) + 2.*f(n-1) )/h**2

      end subroutine







      end module
