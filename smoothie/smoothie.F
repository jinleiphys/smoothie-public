c-----------------------------------------------------------------------
c           a + A -> b + x + A -> x + b + A -> b+anything
c           =====        =====        =====
c             a            x            b
c-----------------------------------------------------------------------

      program smoothie
      use input
      use precision
      use scatt
      use channels
      use bound
      use IAVzerospin
      use IAVdwba
      use iavdwbarbx
      use iavdwbarbx_lagrange
      use iavdwbarb_lagrange, only: NEB_DWBA_Zero_Spin_rb_lagrange
      implicit none
      character*100 line,homedir
      real*8 :: t1, t2
      real*8 :: omp_get_wtime
      t1=omp_get_wtime()
      
      read(5,1005) line
 1005  format(a100)


      call print_logo()
! using the old style input
      call initialize()   ! Parameter initialization
      call check()        !check input and local copy
      call alpha_3b_in()
      call alpha_3b_out()


      select case(dwba)
        case(1)
          call NEB_DWBA_Zero_Spin_rbx()
        case(2)
          call NEB_DWBA_Zero_Spin()
        case(3)
          call NEB_DWBA_with_Spin()
        case(4)
          call NEB_DWBA_Zero_Spin_rbx_lagrange()
        case(5)
          call NEB_DWBA_Zero_Spin_rb_lagrange()       
        case default
          stop "please set dwba =1,2,3,4,5"
      end select


      call fkind()

      call write_cm2lab_input() 

      call init_eps
      call print_eps
! finish old style input
C         if (dwba==2 .and. (prior .eqv. .false.)) write(*,*) "ATTENTION!!!!!! nx > 100 are needed"



      call GETINFO()

      WRITE(*,*) 'smoothie version Date: ',VERDATE
      write(*,*) "smoothie git version: ", VERREV
      write(*,*) "smoothie compile Date:", COMPDATE
       
       
      t2=omp_get_wtime()
      write(*,101) t2-t1
101   format(20x,'(CPU  time =',F12.2,2x,'seconds)')


       

      end program smoothie
      
      subroutine print_logo()
      implicit none
  
      write(*,'(a)') "/*===========================================================================*\\"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   ███████╗███╗   ███╗ ██████╗  ██████╗ ████████╗██╗  ██╗██╗███████╗        |"
      write(*,'(a)') "|   ██╔════╝████╗ ████║██╔═══██╗██╔═══██╗╚══██╔══╝██║  ██║██║██╔════╝        |"
      write(*,'(a)') "|   ███████╗██╔████╔██║██║   ██║██║   ██║   ██║   ███████║██║█████╗          |"
      write(*,'(a)') "|   ╚════██║██║╚██╔╝██║██║   ██║██║   ██║   ██║   ██╔══██║██║██╔══╝          |"
      write(*,'(a)') "|   ███████║██║ ╚═╝ ██║╚██████╔╝╚██████╔╝   ██║   ██║  ██║██║███████╗        |"
      write(*,'(a)') "|   ╚══════╝╚═╝     ╚═╝ ╚═════╝  ╚═════╝    ╚═╝   ╚═╝  ╚═╝╚═╝╚══════╝        |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   **** Scattering Model of Optical Operator Theory for                     |"
      write(*,'(a)') "|   **** Ichimura-Austern-Vincent Equations                                  |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|*==========================================================================*|"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   **** IAV FORMALISM CALCULATIONS                                          |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   **** Non-elastic breakup calculations in inclusive breakup reactions     |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   **** Jin Lei and Antonio M. Moro                                         |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "|   **** a(=b+x) + A -> b + B, where B = (x+A)                               |"
      write(*,'(a)') "|                                                                            |"
      write(*,'(a)') "\\*==========================================================================*/"
      write(*,*)  

c     Display calculation status bar (decorative)
      write(*,'(a)') "╔════════════════════════════════════════════════════════════════════════╗"
      write(*,'(a)') "║ Initializing...  ▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▱▱▱▱▱▱▱▱▱▱  75%                        ║"
      write(*,'(a)') "╚════════════════════════════════════════════════════════════════════════╝"
      write(*,*)  


      end subroutine print_logo


