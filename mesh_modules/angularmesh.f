      module angularmesh
       use mesh
       use precision
       use gauss
       implicit none
       real*8,dimension(:),allocatable :: angx,angw

      contains

       subroutine ang_mesh()
         implicit none
        !  integer :: nx1,nx2
        !  integer :: ix
        !  real*8 :: theta

         allocate(angx(1:nx),angw(1:nx))


C         nx1=nint(nx/2.0_dpreal)
C         nx2=nx-nx1
C
C         call TRNS(nx1,nx2,nx,5.0/180.*pi,25.0/180.*pi,pi,angx,angw)
C
C         do ix=1, nx
C           angw(ix)=sin(angx(ix))*angw(ix)
C           angx(ix)=cos(angx(ix))
C         end do


        call gauleg(nx,-1.0_dpreal,1.0_dpreal,angx,angw)
C
C


  
C        if(prior==.false.) then 
C        call gauleg(nx,0.0_dpreal,1.0_dpreal,angx,angw)
C        do ix=1,nx
C
C          theta=0.25_dpreal*(3*angx(ix)**2+1)*angx(ix)*pi
C          angw(ix)= (2.25_dpreal*angx(ix)*angx(ix)*pi + 0.25_dpreal*pi) * sin(theta) * angw(ix)
C          angx(ix)=cos(theta)
C        end do
C        
C        end if 
C        



       end subroutine
       
       
       
       
       subroutine ang_mesh_post()
         implicit none
        !  integer :: nx1,nx2
         integer :: ix
         real*8 :: theta

         allocate(angx(1:nx),angw(1:nx))


C         nx1=nint(nx/2.0_dpreal)
C         nx2=nx-nx1
C
C         call TRNS(nx1,nx2,nx,5.0/180.*pi,25.0/180.*pi,pi,angx,angw)
C
C         do ix=1, nx
C           angw(ix)=sin(angx(ix))*angw(ix)
C           angx(ix)=cos(angx(ix))
C         end do


C       call gauleg(nx,-1.0_dpreal,1.0_dpreal,angx,angw)
C
C



         
         call gauleg(nx,0.0_dpreal,1.0_dpreal,angx,angw)
         do ix=1,nx
 
           theta=0.25_dpreal*(3*angx(ix)**2+1)*angx(ix)*pi
           angw(ix)= (2.25_dpreal*angx(ix)*angx(ix)*pi + 0.25_dpreal*pi) * sin(theta) * angw(ix)
           angx(ix)=cos(theta)
         end do
         




       end subroutine

      end module
