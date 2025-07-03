      module zerorangesource
      use precision
      use mesh 
      implicit none 
      real*8 ::  D0 ! Zero range factor 
      complex*16, allocatable :: Hp_a_pos(:,:)
      complex*16, allocatable :: Hp_b_pos(:,:), Hm_b_pos(:,:)
      complex*16, allocatable :: Hp_x_pos(:,:)
      complex*16, allocatable :: Hm_a_neg(:,:)
      complex*16, allocatable :: Hp_b_neg(:,:), Hm_b_neg(:,:)
      complex*16, allocatable :: Hp_x_neg(:,:)
      real*8,allocatable :: cph_a(:),cph_b(:),cph_x(:)  !Coulomb phase-shift
      real*8 :: ka_store, kb_store, kx_store ! wave number
c     Scaling factors for H+ and H- functions for each l value and channel
      complex*16, dimension(:), allocatable :: scale_plus_a, scale_minus_a
      complex*16, dimension(:), allocatable :: scale_plus_b, scale_minus_b
      complex*16, dimension(:), allocatable :: scale_plus_x, scale_minus_x
      
c     Flag to indicate if scaling factors have been initialized
      logical :: scaling_initialized = .false.  , ecmxpositive = .true.


      contains 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine init_Pi_ext(counter)
c     This subroutine initializes arrays for Pi_ext calculation
c     Should be called once before using Pi_ext
c     Uses COULCC directly for all Hankel functions except Hm_a_neg,
c     which requires special l-dependent handling when y values are large
      use systems
      use channels
      use constants, only: pi, iu, e2, hbarc, amu, dpreal
      use mesh
      use gauss
      use coulfunc
      implicit none
      
      integer, intent(in) :: counter
      integer :: ifail_total, ifail, iy, l, mode
      real*8 :: ka, kb, kx, c, y
      real*8 :: eta_a, eta_b, eta_x
      real*8 :: mua, mub, mux, ecm, ecmb, ecmx
      complex*16 :: z_pos, z_neg, zlmin
      
c     COULCC calculation arrays
      complex*16, dimension(:), allocatable :: fc, gc, fcp, gcp, sig_i
      
c     Variables for Hm_a_neg special handling
      real*8, dimension(:), allocatable :: y_threshold
      complex*16, dimension(:), allocatable :: scale_factors ! 
      complex*16 :: Hp_asymp, Hm_asymp, h_minus
      real*8, parameter :: min_value = 1.0d-50  ! Threshold to detect COULCC failure/precision loss
      
c     Error handling variables
      integer :: alloc_stat
      logical :: debug_output
      
c     Set debug level 
      debug_output = .false. 
      
c     Initialize error counter
      ifail_total = 0
      
c     Output basic mesh information if debug enabled
      if (debug_output) then
          write(*,*) 'DEBUG: Integration mesh set up with:'
          write(*,*) 'DEBUG: Pi_ext_ny=', Pi_ext_ny, 'Pi_ext_ymax=', 
     &               Pi_ext_ymax
          write(*,*) 'DEBUG: First point:', yy_mesh(1), 
     &               'Last point:', yy_mesh(Pi_ext_ny)
      endif
      
c     Calculate wavenumbers and related parameters
      mua = (massp * masst) * amu / (massp + masst)
      mub = (massb * (massx+masst)) * amu / (massb + massx+masst)
      mux = (massx * masst) * amu / (massx + masst)
      
c     Validate masses to prevent division by zero
      if (abs(massp+masst) < 1.0d-10 .or. 
     &    abs(massb+massx+masst) < 1.0d-10 .or. 
     &    abs(massx+masst) < 1.0d-10) then
          write(*,*) 'ERROR: Near-zero denominator in reduced mass calculation'
          return
      endif
      
      ecm = elab * masst / (masst + massp)
      ecmb = necmb(counter)
      ecmx = ecm + qval - ecmb
      
c     Check for negative arguments to sqrt
      if (ecm < 0.0d0 .or. ecmb < 0.0d0) then
          write(*,*) 'ERROR: Negative energy detected in wavenumber calculation'
          write(*,*) 'ecm=', ecm, 'ecmb=', ecmb
          return
      endif
      
      ka = sqrt(2*mua*ecm/(hbarc**2))
      kb = sqrt(2*mub*ecmb/(hbarc**2))
      
c     Special handling for ecmx potentially negative (bound state)
      if (ecmx >= 0.0d0) then
          kx = sqrt(2*mux*ecmx/(hbarc**2))
          ecmxpositive = .true.
      else
          kx = sqrt(2*mux*abs(ecmx)/(hbarc**2))
          ecmxpositive = .false.
          if (debug_output) then
             write(*,*) 'DEBUG: Negative ecmx=', ecmx, 
     &                  'treated as bound state'
          endif
      endif
      
c     Store wave numbers for later use
      ka_store = ka
      kb_store = kb
      kx_store = kx
      
c     Calculate Sommerfeld parameters with check for division by zero
      if (abs(ka) < 1.0d-10 .or. abs(kb) < 1.0d-10 .or. abs(kx) < 1.0d-10) then
          write(*,*) 'WARNING: Near-zero wavenumber detected, using small value'
          if (abs(ka) < 1.0d-10) ka = 1.0d-10
          if (abs(kb) < 1.0d-10) kb = 1.0d-10
          if (abs(kx) < 1.0d-10) kx = 1.0d-10
      endif
      
      eta_a = zp*zt*e2*mua/hbarc/hbarc/ka
      eta_b = zb*(zt+zx)*e2*mub/hbarc/hbarc/kb
      eta_x = zx*zt*e2*mux/hbarc/hbarc/kx
      
      if (debug_output) then
          write(*,*) 'DEBUG: Wavenumbers: ka=', ka, 'kb=', kb, 'kx=', kx
          write(*,*) 'DEBUG: Sommerfeld params: eta_a=', eta_a, 
     &               'eta_b=', eta_b, 'eta_x=', eta_x
      endif
      
c     Scaling factor for b coordinate
      c = masst/(masst+massx)
      
c     Allocate COULCC calculation and special handling arrays
      allocate(fc(0:max(lmax,lxmax)), 
     &         gc(0:max(lmax,lxmax)),
     &         fcp(0:max(lmax,lxmax)),
     &         gcp(0:max(lmax,lxmax)),
     &         sig_i(0:max(lmax,lxmax)),
     &         y_threshold(0:lmax),
     &         scale_factors(0:lmax), ! Changed dimension
     &         stat=alloc_stat)
      if (alloc_stat /= 0) then
          write(*,*) 'ERROR: Failed to allocate temporary arrays'
          return
      endif

c     Coulomb phase shifts
      if(allocated(cph_a)) deallocate(cph_a, cph_b, cph_x, stat=alloc_stat)
      if (alloc_stat /= 0) then
          write(*,*) 'WARNING: Error in deallocating phase shift arrays'
          if (allocated(fc)) deallocate(fc, gc, fcp, gcp, sig_i)
          if (allocated(y_threshold)) deallocate(y_threshold)
          if (allocated(scale_factors)) deallocate(scale_factors)
c         if (allocated(use_asymp)) deallocate(use_asymp) ! Removed
          return
      endif
      
      allocate(cph_a(0:lmax), cph_b(0:lmax), cph_x(0:lxmax), 
     &         stat=alloc_stat)
      if (alloc_stat /= 0) then
          write(*,*) 'ERROR: Failed to allocate phase shift arrays'
          if (allocated(fc)) deallocate(fc, gc, fcp, gcp, sig_i)
          if (allocated(y_threshold)) deallocate(y_threshold)
          if (allocated(scale_factors)) deallocate(scale_factors)
c         if (allocated(use_asymp)) deallocate(use_asymp) ! Removed
          return
      endif
      
      cph_a = 0.0_dpreal
      cph_b = 0.0_dpreal
      cph_x = 0.0_dpreal
      call coulph(eta_a, cph_a, lmax)
      call coulph(eta_b, cph_b, lmax)
      call coulph(eta_x, cph_x, lxmax)
      
c     Allocate Hankel function arrays with improved error handling
      if (allocated(Hp_a_pos)) then
          deallocate(Hp_a_pos, stat=alloc_stat)
          if (alloc_stat /= 0) 
     &        write(*,*) 'WARNING: Error deallocating Hp_a_pos'
      endif
      
      if (allocated(Hp_b_pos)) then
          deallocate(Hp_b_pos, Hm_b_pos, stat=alloc_stat)
          if (alloc_stat /= 0) 
     &        write(*,*) 'WARNING: Error deallocating Hp/m_b_pos'
      endif
      
      if (allocated(Hp_x_pos)) then
          deallocate(Hp_x_pos, stat=alloc_stat)
          if (alloc_stat /= 0) 
     &        write(*,*) 'WARNING: Error deallocating Hp_x_pos'
      endif
      
      if (allocated(Hm_a_neg)) then
          deallocate(Hm_a_neg, stat=alloc_stat)
          if (alloc_stat /= 0) 
     &        write(*,*) 'WARNING: Error deallocating Hm_a_neg'
      endif
      
      if (allocated(Hp_b_neg)) then
          deallocate(Hp_b_neg, Hm_b_neg, stat=alloc_stat)
          if (alloc_stat /= 0)
     &        write(*,*) 'WARNING: Error deallocating Hp/m_b_neg'
      endif
      
      if (allocated(Hp_x_neg)) then
          deallocate(Hp_x_neg, stat=alloc_stat)
          if (alloc_stat /= 0)
     &        write(*,*) 'WARNING: Error deallocating Hp_x_neg'
      endif
      
c     Allocate all arrays needed
      allocate(Hp_a_pos(Pi_ext_ny, 0:lmax), 
     &         Hp_b_pos(Pi_ext_ny, 0:lmax), 
     &         Hm_b_pos(Pi_ext_ny, 0:lmax),
     &         Hp_x_pos(Pi_ext_ny, 0:lxmax),
     &         Hm_a_neg(Pi_ext_ny, 0:lmax),
     &         Hp_b_neg(Pi_ext_ny, 0:lmax), 
     &         Hm_b_neg(Pi_ext_ny, 0:lmax),
     &         Hp_x_neg(Pi_ext_ny, 0:lxmax),
     &         stat=alloc_stat)
      
      if (alloc_stat /= 0) then
          write(*,*) 'ERROR: Failed to allocate Hankel function arrays'
          if (allocated(fc)) deallocate(fc, gc, fcp, gcp, sig_i)
          if (allocated(y_threshold)) deallocate(y_threshold)
          if (allocated(scale_factors)) deallocate(scale_factors)
          return
      endif

c     Pre-compute all Hankel functions along contours
      write(*,*) 'Pre-computing Hankel functions for Pi_ext...'

c-----------------------------------------------------------------------
c     First step: Determine thresholds and scaling factors for Hm_a_neg
c-----------------------------------------------------------------------
      write(*,*) 'Finding l-dependent thresholds for Hm_a_neg...'
      
c     Initialize thresholds to conservative values
      y_threshold = Pi_ext_ymax 
      scale_factors = (1.0d0, 0.0d0) ! Default scale factor
      
c     Mode for COULCC
      mode = 1
      zlmin = 0.0d0
      
c     For each l-value, find where COULCC starts to fail for Hm_a_neg
      do l = 0, lmax
  !        write(*,*) 'Testing COULCC limits for l =', l
         
         do iy = 1, Pi_ext_ny
            y = yy_mesh(iy)
            z_neg = cmplx(rmax, -y, kind=8)
            
c           Try COULCC
            call COULCC(ka * z_neg, cmplx(eta_a, 0.0d0, kind=8),
     &                  zlmin, l+1, fc, gc, fcp, gcp, sig_i, 
     &                  mode, 0, ifail)
c           Check if COULCC is returning unreliable values (fail flag or too small H-)
            if (ifail /= 0 .or. abs(gc(l) - iu * fc(l)) < min_value) then
c              We found where COULCC starts failing for this l
               y_threshold(l) = yy_mesh(iy-1)
               
c              If COULCC fails at the very first y-value, cannot compute scaling factor.
c              Use default factor = 1 and skip to next l.
               if (iy == 1) then
                  write(*,*) 'WARNING: COULCC fails at all y for l =', l
                  scale_factors(l) = (1.0d0, 0.0d0) ! Use default factor
                  exit ! Exit iy loop for this l
               endif
               
c              Scale factor is based on the last successful y value (iy-1)
               y = yy_mesh(iy-2) 
               z_neg = cmplx(rmax, -y, kind=8)
               
c              Get COULCC value at the last good point
               call COULCC(ka * z_neg, cmplx(eta_a, 0.0d0, kind=8),
     &                     zlmin, l+1, fc, gc, fcp, gcp, sig_i, 
     &                     mode, 0, ifail) 
               ! Ifail should be 0 here, otherwise logic is flawed. Add check?
               if (ifail /= 0) then
                  write(*,*) '*** LOGIC ERROR: COULCC failed at iy-1 for l=', l
                  scale_factors(l) = (1.0d0, 0.0d0) ! Fallback
               else
c                 COULCC H- value
                  h_minus = gc(l) - iu * fc(l)
                  
c                 Get asymptotic value at the same point
                  call coulomb_asymp_large(ka * z_neg, eta_a, l, Hp_asymp, Hm_asymp, cph_a(l))
                                          
                  
c                 Scale factor calculation
                  if (abs(Hm_asymp) > 1.0d-20) then
                     scale_factors(l) = h_minus / Hm_asymp ! Store single factor
                  else
                     scale_factors(l) = (1.0d0, 0.0d0)     ! Fallback if asymptotic is zero
                  endif
               endif

            !    write(*,*) 'For l =', l, ', threshold at y >=', y_threshold(l), " scale_factor = ", scale_factors(l)
               exit ! Exit iy loop for this l
            endif
            
c           If we reach the last y value without COULCC failing, we can use COULCC
c           for all y values for this l. Set threshold beyond mesh.
c           Still compute scaling factor at last point for consistency/potential use.
            if (iy == Pi_ext_ny) then
               write(*,*) 'COULCC works for all y values for l =', l
               y_threshold(l) = Pi_ext_ymax + 1.0d0  ! Set threshold beyond mesh
               
c              Compute scaling factor at the last point (iy)
               y = yy_mesh(iy)
               z_neg = cmplx(rmax, -y, kind=8)
               
c              Get COULCC value (already computed in this iteration)
               h_minus = gc(l) - iu * fc(l) 
               
c              Get asymptotic value at the same point
               call coulomb_asymp_large(ka * z_neg, eta_a, l, 
     &                                  Hp_asymp, Hm_asymp, cph_a(l))
               
c              Scale factor calculation
               if (abs(Hm_asymp) > 1.0d-20) then
                  scale_factors(l) = h_minus / Hm_asymp ! Store single factor
               else
                  scale_factors(l) = (1.0d0, 0.0d0)     ! Fallback
               endif
               ! No need to print scale factor here if COULCC is always used
            endif
         enddo  ! End loop over iy
      enddo  ! End loop over l

c-----------------------------------------------------------------------
c     Second step: Calculate all Hankel functions using appropriate method
c-----------------------------------------------------------------------
      write(*,*) 'Computing all Hankel functions...'
      
c     Loop over y values and calculate all functions
      do iy = 1, Pi_ext_ny
         y = yy_mesh(iy)
         
c        Upper plane (C2): z = Rmax + iy -----------------------------------------
         z_pos = cmplx(rmax, y, kind=8)
         
c        For channel a (H+)
         call COULCC(ka * z_pos, cmplx(eta_a, 0.0d0, kind=8),
     &               zlmin, lmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ifail_total = ifail_total + ifail
         
         if (ifail == 0) then
            do l = 0, lmax
               Hp_a_pos(iy, l) = gc(l) + iu * fc(l)   ! H^(+)
            enddo
         else
            write(*,*) 'WARNING: COULCC failed for Hp_a_pos at y =', y
            do l = 0, lmax
               call coulomb_asymp_large(ka * z_pos, eta_a, l, 
     &                                  Hp_asymp, Hm_asymp, cph_a(l))
               Hp_a_pos(iy, l) = Hp_asymp
            enddo
         endif
         
c        For channel b (H+ and H-)
         call COULCC(kb * c * z_pos, cmplx(eta_b, 0.0d0, kind=8),
     &               zlmin, lmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ifail_total = ifail_total + ifail
         
         if (ifail == 0) then
            do l = 0, lmax
               Hp_b_pos(iy, l) = gc(l) + iu * fc(l)   ! H^(+)
               Hm_b_pos(iy, l) = gc(l) - iu * fc(l)   ! H^(-)
            enddo
         else
            write(*,*) 'WARNING: COULCC failed for Hp/m_b_pos at y =', y
            do l = 0, lmax
               call coulomb_asymp_large(kb * c * z_pos, eta_b, l, 
     &                                  Hp_asymp, Hm_asymp, cph_b(l))
               Hp_b_pos(iy, l) = Hp_asymp
               Hm_b_pos(iy, l) = Hm_asymp
            enddo
         endif
         
c        For channel x (H+)
         call COULCC(kx * z_pos, cmplx(eta_x, 0.0d0, kind=8),
     &               zlmin, lxmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ifail_total = ifail_total + ifail
         
         if (ifail == 0) then
            do l = 0, lxmax
               Hp_x_pos(iy, l) = gc(l) + iu * fc(l)   ! H^(+)
            enddo
         else
            write(*,*) 'WARNING: COULCC failed for Hp_x_pos at y =', y
            do l = 0, lxmax
               call coulomb_asymp_large(kx * z_pos, eta_x, l, 
     &                                  Hp_asymp, Hm_asymp, cph_x(l))
               Hp_x_pos(iy, l) = Hp_asymp
            enddo
         endif
         
c        Lower plane (C3): z = Rmax - iy -----------------------------------------
         z_neg = cmplx(rmax, -y, kind=8)
         
c        For channel a (H-) - Using thresholding and scaling
         call COULCC(ka * z_neg, cmplx(eta_a, 0.0d0, kind=8),
     &               zlmin, lmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ! Note: ifail here might be non-zero, handled below
         
         do l = 0, lmax
c           Check if we are below the threshold AND COULCC call was successful for this z_neg
            if (y < y_threshold(l) .and. ifail == 0) then
c              Use COULCC result directly
               Hm_a_neg(iy, l) = gc(l) - iu * fc(l)   ! H^(-)
               
c              Verify result is not too small (potential COULCC precision issue)
               if (abs(Hm_a_neg(iy, l)) < min_value) then
c                 Use scaled asymptotic form even though below threshold due to precision loss
                  call coulomb_asymp_large(ka * z_neg, eta_a, l, 
     &                                     Hp_asymp, Hm_asymp, cph_a(l))
                  Hm_a_neg(iy, l) = Hm_asymp * scale_factors(l) ! Use stored factor
               endif
            else
c              Use scaled asymptotic form (either y>=threshold OR COULCC failed for this z_neg)
               if (ifail /= 0 .and. l==0) then 
                  ! Print warning only once per failed COULCC call
                  write(*,*) 'WARNING: COULCC failed for Hm_a_neg at y =', y,  ' (Using asymptotic for all l)'
               endif
               call coulomb_asymp_large(ka * z_neg, eta_a, l, 
     &                                  Hp_asymp, Hm_asymp, cph_a(l))
               Hm_a_neg(iy, l) = Hm_asymp * scale_factors(l) ! Use stored factor
            endif
         enddo ! end loop l for Hm_a_neg
         ifail_total = ifail_total + ifail ! Add ifail from this COULCC call
         
c        For channel b (H+ and H-)
         call COULCC(kb * c * z_neg, cmplx(eta_b, 0.0d0, kind=8),
     &               zlmin, lmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ifail_total = ifail_total + ifail
         
         if (ifail == 0) then
            do l = 0, lmax
               Hp_b_neg(iy, l) = gc(l) + iu * fc(l)   ! H^(+)
               Hm_b_neg(iy, l) = gc(l) - iu * fc(l)   ! H^(-)
            enddo
         else
            write(*,*) 'WARNING: COULCC failed for Hp/m_b_neg at y =', y
            do l = 0, lmax
               call coulomb_asymp_large(kb * c * z_neg, eta_b, l, 
     &                                  Hp_asymp, Hm_asymp, cph_b(l))
               Hp_b_neg(iy, l) = Hp_asymp
               Hm_b_neg(iy, l) = Hm_asymp
            enddo
         endif
         
c        For channel x (H+)
         call COULCC(kx * z_neg, cmplx(eta_x, 0.0d0, kind=8),
     &               zlmin, lxmax+1, fc, gc, fcp, gcp, sig_i, 
     &               mode, 0, ifail)
         ifail_total = ifail_total + ifail
         
         if (ifail == 0) then
            do l = 0, lxmax
               Hp_x_neg(iy, l) = gc(l) + iu * fc(l)   ! H^(+)
            enddo
         else
            write(*,*) 'WARNING: COULCC failed for Hp_x_neg at y =', y
            do l = 0, lxmax
               call coulomb_asymp_large(kx * z_neg, eta_x, l, 
     &                                  Hp_asymp, Hm_asymp, cph_x(l))
               Hp_x_neg(iy, l) = Hp_asymp
            enddo
         endif
      end do  ! End loop over iy

c     Perform consistency check on Hm_a_neg values (Optional - Keep commented unless debugging)
!      write(*,*) 'Performing consistency check on Hm_a_neg values...'
!      do iy = 1, Pi_ext_ny
!         do l = 0, lmax
!            if (abs(Hm_a_neg(iy, l)) < 1.0d-30) then ! Use a very small number
!               write(*,*) 'INFO: Near-zero Hm_a_neg detected post-calc at y =', 
!     &                    yy_mesh(iy), 'for l =', l, ' val=', Hm_a_neg(iy, l)
!               ! Could potentially force recalculation here if needed, but might indicate deeper issues
!            endif
!         enddo
!      enddo

      if (ifail_total > 0) then
         write(*,*) 'WARNING: Total COULCC failures encountered:', ifail_total, 
     &              ' (Asymptotic formulas used as fallback)'
      else
         write(*,*) 'All Hankel function calculations completed successfully.'
      endif
      
c     Clean up temporary arrays
      deallocate(fc, gc, fcp, gcp, sig_i)
      deallocate(y_threshold, scale_factors) ! Removed use_asymp

      write(*,*) 'Pi_ext initialization complete.'
      
      end subroutine init_Pi_ext
c-----------------------------------------------------------------------
      subroutine coulomb_asymp_large(z, eta, l, H_plus, H_minus, cph_l)
c     Asymptotic formula for Coulomb Hankel functions for large arguments
c     Parameters:
c       z - complex coordinate (rho = kr)
c       eta - Sommerfeld parameter
c       l - angular momentum
c       H_plus - output for H+ Hankel function
c       H_minus - output for H- Hankel function
c       cph_l - Coulomb phase shift for angular momentum l
c-----------------------------------------------------------------------
      use constants, only: pi, iu
      implicit none
      
      complex*16, intent(in) :: z
      real*8, intent(in) :: eta, cph_l
      integer, intent(in) :: l
      complex*16, intent(out) :: H_plus, H_minus
      
      complex*16 :: common_term, phase_plus, phase_minus
      
c     Common exponential term
      common_term =  exp(iu * (cph_l + 0.5d0*l*pi))
      
c     Phase terms for H+ and H-
      phase_plus = exp(iu * (z - eta * log(2.0d0*z) - 0.5d0*l*pi))
      phase_minus = exp(-iu * (z - eta * log(2.0d0*z) - 0.5d0*l*pi))
      
c     Final Hankel functions
      H_plus = common_term * phase_plus
      H_minus = common_term * phase_minus
      
      end subroutine coulomb_asymp_large


c-----------------------------------------------------------------------
      subroutine coulomb_asymp_small(z, eta, l, H_plus, H_minus, cph_l)
c     Asymptotic formula for Coulomb Hankel functions for small arguments
c     Only used in very specific cases or when COULCC fails
c     Parameters:
c       z - complex coordinate (rho = kr)
c       eta - Sommerfeld parameter
c       l - angular momentum
c       H_plus - output for H+ Hankel function
c       H_minus - output for H- Hankel function
c       cph_l - Coulomb phase shift for angular momentum l
c-----------------------------------------------------------------------
      use constants, only: pi, iu
      implicit none
      
      complex*16, intent(in) :: z
      real*8, intent(in) :: eta, cph_l
      integer, intent(in) :: l
      complex*16, intent(out) :: H_plus, H_minus
      
      real*8 :: factorial_l
      complex*16 :: z_factor, common_term
      integer :: i
      
c     Compute l factorial
      factorial_l = 1.0d0
      do i = 1, l
         factorial_l = factorial_l * i
      enddo
      
c     z^l factor
      z_factor = z**l
      
c     Common term
      common_term = (2.0d0**l) * factorial_l * exp(iu * cph_l) / 
     &              ((2.0d0*l + 1.0d0) * factorial_l)
      
c     Hankel functions approximations
      H_plus = common_term * z_factor
      H_minus = conjg(H_plus)
      
      end subroutine coulomb_asymp_small
c-----------------------------------------------------------------------
      subroutine Pi_ext(alphabar, Pi_ext_val)
c     This subroutine calculates the external part of the wavefunction
c     Pi_ext(r) = -f_(alpha_xA)(r)/4 * [∫_C2 dz U(z) + ∫_C3 dz L(z)]
c     where C2 and C3 are integration paths in upper and lower complex plane
c     Uses pre-computed Hankel functions for efficiency
      use systems
      use channels
      use interpolation
      use scatt
      use constants, only: pi, iu
      use mesh 
      implicit none
      integer, intent(in) :: alphabar
      complex*16, dimension(1:nr), intent(out) :: Pi_ext_val
      
      integer :: ir, alpha_xA, alpha_out, alpha_in
      integer :: lx, la, lb, LL
      real*8 :: y, c
      complex*16 :: Hz, U_val, L_val
      integer :: iy
      complex*16 :: z_pos, z_neg
      integer :: alphaspect_a, alphaspect_b
      complex*16 :: Sa, Sb
      
      ! Arrays to store integration results
      complex*16 :: C2_integral, C3_integral


      Pi_ext_val=0.0_dpreal
      if(.not. ecmxpositive) return
      
c     Get channel indices and angular momenta
      alpha_out = abar%alpha_out(alphabar)
      alpha_in = abar%alpha_in(alphabar)
      alpha_xA = out3b%alpha2b(alpha_out)
      
      lx = out3b%l(alpha_out)
      la = in3b%lam(alpha_in)
      lb = out3b%lam(alpha_out)
      LL = nint(out3b%j(alpha_out))
      
      alphaspect_a = in3b%alphaspect(alpha_in)
      alphaspect_b = out3b%alphaspect(alpha_out)
      
c     Get pre-computed S-matrices from scatt module
      Sa = smat_a(alphaspect_a)
      Sb = smat_b(alphaspect_b)
      
c     Scaling factor for b coordinate
      c = masst/(masst+massx)


      
c     Pre-compute the integrals which are independent of r
      C2_integral = 0.0_dpreal
      C3_integral = 0.0_dpreal
      
c     Integration over C2 (upper plane)
      do iy = 1, Pi_ext_ny
         y = yy_mesh(iy)
         z_pos = cmplx(rmax, y, kind=8)
         
         Hz = 1.0/z_pos
         
c        Calculate U(z) term using Hankel functions directly
         U_val = Hz * (-Sa * Hp_x_pos(iy, lx) * Hm_b_pos(iy, lb) * 
     &                 Hp_a_pos(iy, la) + 
     &                 Sa * Sb * Hp_x_pos(iy, lx) * Hp_b_pos(iy, lb) * 
     &                 Hp_a_pos(iy, la))
         
c        Add to integral with weight
         C2_integral = C2_integral + U_val  * yyw(iy)
      !    write(213,*) y,abs(U_val)
      end do
      ! write(213,*) "&"
      
c     Integration over C3 (lower plane)
      do iy = 1, Pi_ext_ny
         y = yy_mesh(iy)
         z_neg = cmplx(rmax, -y, kind=8)
         
         Hz = 1.0/z_neg
         
c        Calculate L(z) term using Hankel functions directly
         L_val = Hz * (Hp_x_neg(iy, lx) * Hm_b_neg(iy, lb) * 
     &                 Hm_a_neg(iy, la) - 
     &                 Sb * Hp_x_neg(iy, lx) * Hp_b_neg(iy, lb) * 
     &                 Hm_a_neg(iy, la))
         
c        Add to integral with weight (dz = -i*dy for lower contour)
         C3_integral = C3_integral + L_val * yyw(iy)
      !    write(214,*) y, abs(L_val)
         write(215,*)y, abs(Hp_x_neg(iy, lx)*Hp_b_neg(iy, lb)*Hm_a_neg(iy, la))
         write(216,*)y,abs(Hm_a_neg(iy, la))
      end do
      ! write(214,*) "&"
      write(215,*) "&alpha_in=",alpha_in,"alpha_out=",alpha_out
      write(216,*) "&alpha_in=",alpha_in,"alpha_out=",alpha_out
      
c     Now apply the results to each radial point
      do ir = 1, nr
         Pi_ext_val(ir) = -0.25_dpreal * Gx%re(ir,alpha_xA) * 
     &                    (iu*C2_integral - iu*C3_integral)
      end do
      
c     Output for debugging/visualization
      ! write(112,*) '&alpha_in=',alpha_in,'alpha_out=',alpha_out
      ! do ir = 1, nr
      !    write(112,*) rr(ir), real(Pi_ext_val(ir)), aimag(Pi_ext_val(ir))
      ! end do
      
      end subroutine Pi_ext
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine R_ext(alphabar, R_ext_val)
c     This subroutine calculates R^ext(r,alphabar) using the formula
c     R^ext = -32*pi^2*mu_x/(hbar^2*k_a*k_b*k_x) * (1/cr) * i^(la-lb) *
c             exp(i*sigma_la + i*sigma_lb) * G * D0 * Pi^ext(r)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use constants, only: pi, iu, hbarc, amu
      use mesh, only: rr, nr
      implicit none
      integer, intent(in) :: alphabar
      complex*16, dimension(1:nr), intent(out) :: R_ext_val
      
      integer :: ir, alpha_out, alpha_in
      integer :: la, lb
      real*8 :: r, c, mux
      real*8 :: ka, kb, kx
      real*8 :: Gabar
      complex*16, dimension(1:nr) :: Pi_ext_val
      complex*16 :: phase_factor
      real*8 :: prefactor
      
c     Get channel indices and angular momenta
      alpha_out = abar%alpha_out(alphabar)
      alpha_in = abar%alpha_in(alphabar)
      
      la = in3b%lam(alpha_in)
      lb = out3b%lam(alpha_out)
      
c     Calculate reduced masses
      mux = (massx * masst) * amu / (massx + masst)
      
c     Calculate energies and wave numbers
      ka = ka_store
      kb = kb_store
      kx = kx_store
      
c     Scaling factor for b coordinate
      c = masst/(masst+massx)
      
c     Calculate G factor
      call G_zerorange(alphabar, Gabar)
      
c     Calculate Pi_ext
      call Pi_ext(alphabar, Pi_ext_val)
      
c     Calculate the phase factor e^{i(sigma_la + sigma_lb)}
      phase_factor = exp(iu * (cph_a(la) + cph_b(lb)))
      
c     Calculate the prefactor
      prefactor = -32.0_dpreal * pi*pi * mux / (hbarc*hbarc * ka * kb * kx)
      
c     Calculate R_ext for all radial points
      do ir = 1, nr
        r = rr(ir)
        R_ext_val(ir) = prefactor * (1.0_dpreal/(c*r)) * (iu**(la-lb)) * 
     &                  phase_factor * Gabar * D0 * Pi_ext_val(ir)
      end do
c     Output for debugging/visualization
      ! write(213,*) '&alpha_in=',alpha_in,'alpha_out=',alpha_out
      ! do ir = 1, nr
      !   write(213,*) rr(ir), real(R_ext_val(ir)), aimag(R_ext_val(ir))
      ! end do
      
      end subroutine R_ext
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rho_func_zero(alphabar,rho_zero)
c     this subroutine is used to calculate the source term function in the zero
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use interpolation
      use scatt
      implicit none
      complex*16,dimension(1:nr) :: rho_zero
      integer :: ir, alphabar, ainspect,alphaspect
      real*8 :: r, c,  ra, rb
      real*8 :: Gabar
      integer :: alpha_in, alpha_out
      integer :: LL,la, lb
      complex*16 :: wfa,wfb

      rho_zero=0.0_dpreal

      alpha_out=abar%alpha_out(alphabar)
      alpha_in=abar%alpha_in(alphabar)
      LL=nint(out3b%j(alpha_out))
      la=in3b%lam(alpha_in)
      lb=out3b%lam(alpha_out)
      if (la/=LL) return 


      ainspect=in3b%alphaspect(alpha_in)
      alphaspect=out3b%alphaspect(alpha_out)

      c=masst/(masst+massx)
      ! call D0factor()
      call G_zerorange(alphabar,Gabar)



      do ir=1, nr 
         r=rr(ir)
         ra = r 
         rb = c * r 
         wfa=FFC(ra/hcm,wf_a(0:irmatch,ainspect),irmatch+1)*iu**la / ra
         wfb=FFC(rb/hcm,wf_b(0:irmatch,alphaspect),irmatch+1)*iu**(-lb)/ rb
         rho_zero(ir)= Gabar * wfa * wfb * D0
      end do 
      ! write(211,*) '&alpha_in=',alpha_in,'alpha_out=',alpha_out
      ! do ir=1, nr 
      !   r=rr(ir)
      !   write(211,*) r, real(rho_zero(ir)), aimag(rho_zero(ir))
      ! end do 





      end subroutine 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine D0factor()
c     this subroutine is used to calculate the D0 factor 
c     used in the zero range approximation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use pot
      use interpolation
      use bound 
      use systems
      implicit none 
      ! integer :: n,l
      real*8 :: r
      integer :: ir

      d0=0.d0
      do ir=1,nr
         r=rr(ir)
         d0=d0 + r * FFR4(r/hcm,vbx(0:irmatch,1),irmatch+1) * 
     &         FFR4(r/hcm,phi_bx(0:irmatch,1),irmatch+1) * rrw(ir)

      end do 

      d0 = sqrt(4*pi)*d0

      write(*,*)'The following SINGLE-PARTICLE FORM FACTORS are', 
     &                                            'constructed :'
     
      write(*,70)
      write(*,71) be,d0
70    format('k1',2x,'k2',5x,'be',7x,'D0')
71    format('p',3x,'1',3x,F7.4,1x,F10.4)
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine G_zerorange(alphabar,Gabar)
C     this subroutine is used to calculate the G factor in the zero 
c     range approximation 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use channels
      use mesh
      use cleb_coefficient
      implicit none
      integer :: alphabar,alpha_in, alpha_out
      integer :: ma
      integer :: LL,la,l2b,lb,lx
      real*8 ::j_out
      real*8 :: CG1, CG0
      real*8 :: Gabar
      real*8 :: func



        alpha_out=abar%alpha_out(alphabar)
        alpha_in=abar%alpha_in(alphabar)
        j_out=out3b%j(alpha_out)
        LL=nint(j_out)
        la=in3b%lam(alpha_in)
        l2b=in3b%l(alpha_in)
        lx=out3b%l(alpha_out)
        lb=out3b%lam(alpha_out)
        if (la/=LL) then
          write(*,*)'Error in G_zerorange: la/=LL'
          stop
        end if
        Gabar=0.0_dpreal

        CG0=cleb(lx*2,0,lb*2,0,la*2,0)
        func=hat(lx)*hat(lb)/ hat(la)**3 *hat(l2b)*CG0/sqrt(4.0_dpreal*pi)
        func=func * (-1)**(lx+lb-la)

        CG1=0.0_dpreal
        do ma=-la,la
          CG1= CG1 + cleb(l2b*2,0,la*2,ma*2,la*2,ma*2)
        end do

        Gabar=func*CG1







      end subroutine
c-----------------------------------------------------------------------






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function hat(n)
c     this function is used to calculate \hat{n}
c     hat=sqrt(convert*(2*n+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      use constants, only: convert
      implicit none 
      integer, intent(in) :: n
      if (n < 0) then
          write(*,*) 'Error: Negative argument in hat function'
         hat = 0.0_dpreal  ! Return safe value
      else
         hat = sqrt(convert*(2*n+1))
      endif
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      end module 