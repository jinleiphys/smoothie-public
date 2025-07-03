      module fuspot
      use systems
      use pot
      use precision
      use mesh
      use input,only:written
      contains
      subroutine icf_fus_pot(Wfus)
        implicit none
        real*8,dimension(:),allocatable :: Wfus
      !   real*8 :: cutr
      !   real*8 :: a1,a2,a13,rw  
        integer :: ir
        
C       do counter=1,10
C          if (nkp1(counter)=='f') then
C             a1=na1(counter)
C             a2=na2(counter)
C             rw=nrw(counter)
C             exit
C          end if
C       end do
C       a13=a1**(1./3.)+a2**(1./3.)
C       cutr=rw*a13
        
C             call potr('x',1,zx*zt,0.0_dpreal)
C
C             Wfus=aimag(v)
C             write(*,*) "rw=",rw
              call potr('f',1,zx*zt,0.0_dpreal,0)
              Wfus=aimag(v)
           written(33)=.TRUE.
           do ir=0,irmatch
C             if(hcm*ir > cutr) Wfus(ir)=0.0_dpreal
              write(33,*) hcm*ir, Wfus(ir)
           end do
           write(33,*) "&"   
      end subroutine 
      end module