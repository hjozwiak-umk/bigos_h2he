module PROPAGATORS
   !! this modules contains the subroutines used by the propagator:
   !! 1. cenitrfugal_matrix - prepares the centrifugal term (Eq. (6.19))
   !! 2. pes_contribution   - prepares the interaction energy term (Eq. (6.21))
   !! 3. radtermvalue - returns the value of a radial coupling term at given R
   !! 4. calculate_log_der_matrix - calculates the log-derivative matrix (Eq. 6.29)
   !! 5. pes_diagonalization - diagonalizes the coupling matrix from rmin to rmax
   !! 6. numerov - renormalized Numerov's algorithm                                      
   !----------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use array_operations_mod, only: invert_symmetric_matrix, fill_symmetric_matrix
   use radial_coupling_terms_mod, only: get_radial_coupling_term_value
   implicit none
   contains
   !---------------------------------------------------------------------------!
      subroutine cenitrfugal_matrix(number_of_channels,jj,                     &
         channels_level_indices,channels_omega_values,centmatrix)    
         !! calculates the (R**2)*centrifugal matrix from Eq. (6.19)
         !! only called once at the beginning of the calculations
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: jj
            !! total angular momentum
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         real(dp), intent(out) :: centmatrix(number_of_channels, number_of_channels)
            !! (output) - (R**2)*centrifugal matrix
         !---------------------------------------------------------------------!
         integer(int32) :: omegatmp, omegaptmp, v1tmp, j1tmp, v1ptmp, j1ptmp,  &
            ii, ij
         real(dp) :: centtmp, delta1, delta2
         !---------------------------------------------------------------------!
         centmatrix  = 0

         do ii = 1, number_of_channels
            v1tmp = v1array(channels_level_indices(ii))
            j1tmp = j1array(channels_level_indices(ii))
            omegatmp = channels_omega_values(ii)
            delta1 = 0.d0
            if (omegatmp.eq.0) delta1 = 1.d0
            do ij = 1, ii
               v1ptmp = v1array(channels_level_indices(ij))
               j1ptmp = j1array(channels_level_indices(ij))
               omegaptmp = channels_omega_values(ij)
               delta2 = 0.d0
               if (omegaptmp.eq.0) delta2 = 1.d0
               if (v1tmp.ne.v1ptmp) cycle
               if (j1tmp.ne.j1ptmp) cycle
               if (abs(omegatmp-omegaptmp).gt.1) cycle
               if (omegatmp.eq.omegaptmp) then
                  centmatrix(ii, ij) = jj*(jj+1) +j1tmp*(j1tmp+1)-2*omegatmp**2
               else if (omegatmp-omegaptmp.eq.1) then
                  centtmp = -dsqrt(dfloat((jj*(jj+1)-&
                   omegatmp*(omegatmp-1))*(j1tmp*(j1tmp+1)&
                    -omegatmp*(omegatmp-1))))
                  centmatrix(ii, ij) = centtmp * dsqrt(1.d0+delta1)*dsqrt(1.d0+delta2)
               else if (omegatmp-omegaptmp.eq.-1) then
                  centtmp = -dsqrt(dfloat((jj*(jj+1)-omegatmp*(omegatmp+1))    &
                     *(j1tmp*(j1tmp+1)-omegatmp*(omegatmp+1))))
                  centmatrix(ii, ij) = centtmp*dsqrt(1.d0+delta1)*dsqrt(1.d0+delta2)
               endif               

            enddo
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(centmatrix, 'u')
         !---------------------------------------------------------------------!
      end subroutine
!------------------------------------------------------------------------------!
      subroutine pes_contribution(number_of_channels,jj,r,                     &
         channels_level_indices,channels_omega_values,                         &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients,                              &
         nonzero_terms_per_element,nonzero_legendre_indices,                   &
         nonzero_coupling_coefficients,vmatrix)
         !! calculates the contribution from the PES in (X) at given R
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: jj
            !! total angular momentum
         real(dp), intent(in) :: r
            !! intermolecular distance
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: number_of_nonzero_coupling_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element of the coupling matrix
         integer(int32), intent(in) :: number_of_nonzero_coupling_coefficients
            !! number of all non-zero algberaix coefficients in the whole coupling matrix
         integer(int32), intent(in) :: nonzero_terms_per_element(number_of_nonzero_coupling_matrix_elements)
            !! keeps the number of non-zero terms in the sum (Eq. (6.21)) for
            !! each non-zero element of W/V
         integer(int32), intent(in) :: nonzero_legendre_indices(number_of_nonzero_coupling_coefficients)
            !! holds the proper indices pointing to l1/l2/lltabs, which
            !! correspond to the non-vanishing elements of the sum  (Eq. (6.21))
            !! for each non-zero element of W/V
         real(dp), intent(in) :: nonzero_coupling_coefficients(number_of_nonzero_coupling_coefficients)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(out) :: vmatrix(number_of_channels, number_of_channels)
            !! (output) - the PES contribution to the coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: indvl, inonzero_coupling_matrix_elements, omegatmp, omegaptmp, l1,          &
            nonzerolam, v1tmp, j1tmp, v1ptmp, j1ptmp, ii, ij, il
         real(dp) :: erot, sumtemp, pscoeff, v
         !---------------------------------------------------------------------!
         vmatrix   = 0
         indvl     = 0
         inonzero_coupling_matrix_elements = 0
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            v1tmp    = v1array(channels_level_indices(ii))
            j1tmp    = j1array(channels_level_indices(ii))
            omegatmp = channels_omega_values(ii)
            erot     = elevel(channels_level_indices(ii))               
            do ij = 1, ii
               v1ptmp    = v1array(channels_level_indices(ij))
               j1ptmp    = j1array(channels_level_indices(ij))
               omegaptmp = channels_omega_values(ij)
               if (omegatmp.ne.omegaptmp) cycle
               inonzero_coupling_matrix_elements  = inonzero_coupling_matrix_elements + 1
               nonzerolam = nonzero_terms_per_element(inonzero_coupling_matrix_elements)
               sumtemp    = 0.d0
               do il = 1, nonzerolam
                  indvl   = indvl+1
                  l1      = l1tab(nonzero_legendre_indices(indvl))
                  pscoeff = nonzero_coupling_coefficients(indvl)
                  call get_radial_coupling_term_value(r,l1,v1tmp,j1tmp,v1ptmp,j1ptmp,v)
                  sumtemp = sumtemp+pscoeff*v
               enddo
               vmatrix(ii,ij) = -2*reducedmass*sumtemp
             enddo
             vmatrix(ii,ii) = vmatrix(ii,ii)-(2*reducedmass*erot) 
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(vmatrix,'u')
         !---------------------------------------------------------------------!
      end subroutine pes_contribution
!------------------------------------------------------------------------------!
     
!------------------------------------------------------------------------------!
      subroutine calculate_log_der_matrix(h,y_dim,tt_min,tt_n,tt_plus,r_n,r_plus,log_der_matrix)
         !! calculates the log-derivative matrix from Eq. (6.29)
         !! called by numerov and log_derivative at the end of the propagation
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: h
            !! step of the propagation
         integer(int32), intent(in) :: y_dim
            !! dimension of the log-derivative matrix
         real(dp), intent(in) :: tt_min(y_dim, y_dim)
            !! T-matrix at R_{max - 1}
         real(dp), intent(in) :: tt_n(y_dim, y_dim)
            !! T-matrix at R_{max}
         real(dp), intent(in) :: tt_plus(y_dim, y_dim)
            !! T-matrix at R_{max + 1} 
         real(dp), intent(in) :: r_n(y_dim, y_dim)
            !! R-matrix at R_{max} 
         real(dp), intent(in) :: r_plus(y_dim, y_dim)
            !! R-matrix at R_{max + 1}
         real(dp), intent(inout) :: log_der_matrix(y_dim, y_dim)
            !! log-derivative matrix
         !---------------------------------------------------------------------!
         integer(int32) :: i, j
         real(dp) :: a(y_dim,y_dim), b(y_dim,y_dim), c(y_dim,y_dim),           &
            d(y_dim,y_dim), e(y_dim,y_dim), f(y_dim,y_dim), g(y_dim,y_dim),    &
            ab(y_dim,y_dim), abc(y_dim,y_dim), de(y_dim,y_dim),                &
            def(y_dim,y_dim), bra(y_dim,y_dim)
         !---------------------------------------------------------------------!
         log_der_matrix = 0
         !---------------------------------------------------------------------!
         do i=1,y_dim
            do j=1,y_dim
               a(i,j) = - tt_plus(i,j)
               b(i,j) = - tt_plus(i,j)
               c(i,j) = r_plus(i,j)
               d(i,j) = - tt_min(i,j)
               e(i,j) = - tt_min (i,j)
               f(i,j) = r_n(i,j)
               g(i,j) = - tt_n(i,j)
            end do
            a(i,i) = 0.5 + a(i,i)
            b(i,i) = 1. + b(i,i)
            d(i,i) = 0.5 + d(i,i)
            e(i,i) = 1. + e(i,i)
            g(i,i) = 1. + g(i,i)
         end do
         !---------------------------------------------------------------------!
         call invert_symmetric_matrix(b)
         call fill_symmetric_matrix(b, "u")
         call invert_symmetric_matrix(e)
         call fill_symmetric_matrix(e, "u")
         call invert_symmetric_matrix(f)
         call fill_symmetric_matrix(f, "u")
         !---------------------------------------------------------------------!
         CALL DGEMM('N','N',y_dim,y_dim,y_dim,1.0d0,&
                    a,y_dim,b,y_dim,0.d0,ab,y_dim)
         CALL DGEMM('N','N',y_dim,y_dim,y_dim,1.0d0,&
                    ab,y_dim,c,y_dim,0.d0,abc,y_dim)

         CALL DGEMM('N','N',y_dim,y_dim,y_dim,1.0d0,&
                    d,y_dim,e,y_dim,0.d0,de,y_dim)
         CALL DGEMM('N','N',y_dim,y_dim,y_dim,1.0d0,&
                    de,y_dim,f,y_dim,0.d0,DEF,y_dim)
        !----------------------------------------------------------------------!
         bra = abc - def
         !---------------------------------------------------------------------!
         CALL DGEMM('N','N',y_dim,y_dim,y_dim,dble(1.d0/h),bra,&
                    y_dim,g,y_dim,0.d0,log_der_matrix,y_dim)
         !---------------------------------------------------------------------!
      end subroutine calculate_log_der_matrix
      !------------------------------------------------------------------------!
      subroutine numerov(channels_level_indices,channels_omega_values,         &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients,nonzero_terms_per_element,    &
         nonzero_legendre_indices,nonzero_coupling_coefficients,nsteps,        &
         number_of_channels,jj,log_der_matrix)
         !! renormalized Numerov propagator
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: number_of_nonzero_coupling_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element of the coupling matrix
         integer(int32), intent(in) :: number_of_nonzero_coupling_coefficients
            !! number of all non-zero algberaix coefficients in the whole coupling matrix
         integer(int32), intent(in) :: nonzero_terms_per_element(number_of_nonzero_coupling_matrix_elements)
            !! keeps the number of non-zero terms in the sum (Eq. (6.21)) for
            !! each non-zero element of W/V
         integer(int32), intent(in) :: nonzero_legendre_indices(number_of_nonzero_coupling_coefficients)
            !! holds the proper indices pointing to l1/l2/lltabs, which
            !! correspond to the non-vanishing elements of the sum  (Eq. (6.21))
            !! for each non-zero element of W/V
         real(dp), intent(in) :: nonzero_coupling_coefficients(number_of_nonzero_coupling_coefficients)
            !! holds the values of the non-zero algebraic coefficients
         integer(int32), intent(in) :: nsteps
            !! number of steps from rmin to rmax
         integer(int32), intent(in) :: jj
            !! total angular momentum
         real(dp), intent(inout) :: log_der_matrix(:,:)
            !! resulting log-derivative matrix at RMAX  
         !---------------------------------------------------------------------!
         integer(int32) :: i, ii, ij
         real(dp) :: start, finish, r, step_numerov
         real(dp), allocatable, dimension(:,:) :: umatrix, r_temp, rmatrix,    &
            t_minus, t_center, t_plus, r_center, r_plus, w_tmp, v_tmp
         real(dp), dimension(number_of_channels,number_of_channels) ::         &
            cent_mat, v_mat, w_mat, t_mat, u_left, u_right
         !---------------------------------------------------------------------!
         CALL CPU_TIME(start)
         !---------------------------------------------------------------------!
         ! Calculate the centrifugal term
         !---------------------------------------------------------------------!
         call cenitrfugal_matrix(number_of_channels,jj,channels_level_indices, &
            channels_omega_values,cent_mat)
         step_numerov = (rmax - rmin)/dble(nsteps - 1)
         call allocate_2d(R_temp,number_of_channels,number_of_channels)
         !---------------------------------------------------------------------!
         ! Calculate the PES contribution at rmin
         !---------------------------------------------------------------------!
         call pes_contribution(number_of_channels,jj,rmin,                     &
            channels_level_indices,channels_omega_values,                      &
            number_of_nonzero_coupling_matrix_elements,                        &
            number_of_nonzero_coupling_coefficients,                           &
            nonzero_terms_per_element,nonzero_legendre_indices,                &
            nonzero_coupling_coefficients,v_mat)
         !---------------------------------------------------------------------!
         ! Coupling matrix W at rmin
         !---------------------------------------------------------------------!
         w_mat = v_mat - (1./rmin**2.)*cent_mat
         !---------------------------------------------------------------------!
         ! T-matrix (Eq. 6.23)
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            do ij = 1, number_of_channels
               t_mat(ii, ij) = - ((step_numerov**2.)/12.) * (w_mat(ii, ij))
            enddo
            t_mat(ii, ii) = - ((step_numerov**2.)/12.)*&
                                  (2*reducedmass*ETOTAL())+t_mat(ii, ii)
         enddo
         !---------------------------------------------------------------------!
         ! U-matrix (Eq. 6.25)
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            do ij = 1, number_of_channels
               u_left(ii, ij) = - t_mat(ii, ij)
            end do
            u_left(ii, ii) = 1.d0 + u_left(ii, ii)
         end do

         call invert_symmetric_matrix(u_left)
         call fill_symmetric_matrix(u_left, 'u')
         !---------------------------------------------------------------------!
         ! R-matrix at rmin + 1 = U-matrix at rmin
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            do ij = 1, number_of_channels
               r_temp(ii, ij) = 12.d0*u_left(ii, ij)
            enddo
            r_temp(ii, ii) = r_temp(ii, ii) - 10.d0
         enddo
         !---------------------------------------------------------------------!
         ! Continue the propagation to rmax
         !---------------------------------------------------------------------!
         do i=2,nsteps
            !------------------------------------------------------------------!
            ! Coupling matrix W at R
            !------------------------------------------------------------------!
            R = rmin + (i-1)*step_numerov
            call allocate_2d(rmatrix,number_of_channels,number_of_channels)
            call allocate_2d(umatrix,number_of_channels,number_of_channels)
            call pes_contribution(number_of_channels,jj,R,                     &
               channels_level_indices,channels_omega_values,                   &
               number_of_nonzero_coupling_matrix_elements,                     &
               number_of_nonzero_coupling_coefficients,                        &
               nonzero_terms_per_element,nonzero_legendre_indices,             &
               nonzero_coupling_coefficients,v_mat)
            w_mat = v_mat - (1./R**2.)*cent_mat
            !------------------------------------------------------------------!
            ! T-matrix at R (6.18)
            !------------------------------------------------------------------!
            do ii = 1, number_of_channels
               do ij= 1, number_of_channels
                  t_mat(ii, ij) = - ((step_numerov**2.)/12.) * (w_mat(ii, ij))
               enddo
               t_mat(ii, ii) = - ((step_numerov**2.)/12.) *                    &
                  (2*reducedmass*ETOTAL()) + t_mat(ii, ii)
            enddo
            !------------------------------------------------------------------!
            ! U-matrix at R (6.20)
            !------------------------------------------------------------------!
            do ii = 1, number_of_channels
               do ij = 1, number_of_channels
                  u_left(ii, ij) = - t_mat(ii, ij)
               end do
               u_left(ii, ii) = 1.d0 + u_left(ii, ii)
            end do
            call invert_symmetric_matrix(u_left)
            call fill_symmetric_matrix(u_left, 'u')
            do ii = 1, number_of_channels
               do ij = 1, number_of_channels
                  umatrix(ii, ij) = 12.d0*u_left(ii, ij)
               enddo
               umatrix(ii, ii) = umatrix(ii, ii) - 10.d0
            enddo

            call invert_symmetric_matrix(r_temp)
            call fill_symmetric_matrix(r_temp, 'u')
            !------------------------------------------------------------------!
            ! Prepare T at Rmax - 1 and R at Rmax
            !------------------------------------------------------------------!
            if (i == nsteps - 1) then
               call allocate_2d(v_tmp,number_of_channels,number_of_channels)
               call allocate_2d(w_tmp,number_of_channels,number_of_channels)
               call allocate_2d(t_minus,number_of_channels,number_of_channels)

               call pes_contribution(number_of_channels,jj,r,                  &
                  channels_level_indices,channels_omega_values,                &
                  number_of_nonzero_coupling_matrix_elements,                  &
                  number_of_nonzero_coupling_coefficients,                     &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp - (1./R**2.)*cent_mat

               do ii = 1, number_of_channels
                  do ij= 1, number_of_channels
                     t_minus(ii, ij) = -((step_numerov**2.)/12.)*&
                                              (w_tmp(ii, ij))
                  enddo
                  t_minus(ii, ii) = -((step_numerov**2.)/12.) *                &
                     (2*reducedmass*ETOTAL()) + t_minus(ii, ii)
               enddo

               call allocate_2d(r_center,number_of_channels,number_of_channels)
               r_center = umatrix - r_temp
               !---------------------------------------------------------------!
               ! Prepare T at Rmax and R max + 1, and R at Rmax + 1
               !---------------------------------------------------------------!
            else if (i == nsteps) then
               call allocate_2d(v_tmp,number_of_channels,number_of_channels)
               call allocate_2d(w_tmp,number_of_channels,number_of_channels)
               call allocate_2d(t_center,number_of_channels,number_of_channels)

               call pes_contribution(number_of_channels,jj,r,                  &
                  channels_level_indices,channels_omega_values,                &
                  number_of_nonzero_coupling_matrix_elements,                  &
                  number_of_nonzero_coupling_coefficients,                     &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp - (1./R**2.)*cent_mat

               do ii = 1, number_of_channels
                  do ij= 1, number_of_channels
                     t_center(ii, ij) = -((step_numerov**2.)/12.)*             &
                                               (w_tmp(ii, ij))
                  enddo 
                  t_center(ii, ii) = -((step_numerov**2.)/12.) *               &
                     (2*reducedmass*ETOTAL()) + t_center(ii, ii)
               enddo

               call allocate_2d(r_plus,number_of_channels,number_of_channels)
               r_plus = umatrix - r_temp

               r = rmax + step_numerov
               call allocate_2d(v_tmp,number_of_channels,number_of_channels)
               call allocate_2d(t_plus,number_of_channels,number_of_channels)
               call allocate_2d(w_tmp,number_of_channels,number_of_channels)

               call pes_contribution(number_of_channels,jj,r,                  &
                  channels_level_indices,channels_omega_values,                &
                  number_of_nonzero_coupling_matrix_elements,                  &
                  number_of_nonzero_coupling_coefficients,                     &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp - (1./R**2.)*cent_mat

               do ii = 1, number_of_channels
                  do ij= 1, number_of_channels
                     t_plus(ii, ij) = -((step_numerov**2.)/12.) *              &
                        (w_tmp(ii, ij))
                  enddo
                  t_plus(ii, ii) = -((step_numerov**2.)/12.) *                 &
                     (2*reducedmass*ETOTAL()) + t_plus(ii, ii)
               enddo

            end if
            !------------------------------------------------------------------!
            ! R-matrix at R_{n+1} (Eq. 6.28)
            !------------------------------------------------------------------!
            rmatrix = umatrix - r_temp
            call allocate_2d(r_temp,number_of_channels,number_of_channels)
            !------------------------------------------------------------------!
            ! Move R_{n+1} to R_{n}
            !------------------------------------------------------------------!
            r_temp = rmatrix

         end do
         CALL CPU_TIME(finish)
         !---------------------------------------------------------------------!
         ! Eq. (6.29)
         !---------------------------------------------------------------------!
         call calculate_log_der_matrix(step_numerov,number_of_channels,        &
            t_minus,t_center,t_plus,r_center,r_plus,log_der_matrix)
         !---------------------------------------------------------------------!
         if (prntlvl.ge.2) then
            call write_message("Numerov propagator took " //                   &
               trim(adjustl(float_to_character(finish-start, "(E14.8)"))) //   &
               " seconds" )
         endif
         !---------------------------------------------------------------------!
      end subroutine numerov
!------------------------------------------------------------------------------!
end module PROPAGATORS
