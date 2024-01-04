module PROPAGATORS
   !! this modules contains the subroutines used by the propagator:
   !! 1. centrifugal_matrix - prepares the centrifugal term (Eq. (6.19))
   !! 2. pes_contribution   - prepares the interaction energy term (Eq. (6.21))
   !! 4. calculate_log_der_matrix - calculates the log-derivative matrix (Eq. 6.29)
   !! 6. numerov - renormalized Numerov's algorithm                                      
   !----------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use array_operations_mod, only: invert_symmetric_matrix, fill_symmetric_matrix
   use radial_coupling_terms_mod, only: get_radial_coupling_term_value
   implicit none
   contains
   !---------------------------------------------------------------------------!
   !                           Centrifugal matrix
   !---------------------------------------------------------------------------!
      subroutine centrifugal_matrix(total_angular_momentum_,&
         channel_indices_,channels_omega_values_,centrifugal_matrix_)    
         !! calculates the (R**2)*centrifugal matrix from the second term
         !! of Eq. 3 in "What are coupled equations?" section;
         !! Matrix elements are given in Eq. 4 and 6 of "Coupling Matrix" secion
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(:)
            !! holds all values of \bar{\Omega}
         real(dp), intent(out) :: centrifugal_matrix_(:,:)
            !! (output) - (R**2)*centrifugal matrix
         !---------------------------------------------------------------------!
         integer(int32) :: omega_, omega_prime_, v_, j_, v_prime_, j_prime_,  &
            channel_index_1_, channel_index_2_
         real(dp) :: centtmp, delta_1_, delta_2_
         !---------------------------------------------------------------------!
         centrifugal_matrix_  = 0

         do channel_index_1_ = 1, size(channel_indices_)
            v_ = v1array(channel_indices_(channel_index_1_))
            j_ = j1array(channel_indices_(channel_index_1_))
            omega_ = channels_omega_values_(channel_index_1_)
            delta_1_ = delta_for_zero_omega(omega_)
            do channel_index_2_ = 1, channel_index_1_
               v_prime_ = v1array(channel_indices_(channel_index_2_))
               j_prime_ = j1array(channel_indices_(channel_index_2_))
               omega_prime_ = channels_omega_values_(channel_index_2_)
               delta_2_ = delta_for_zero_omega(omega_prime_)
               !---------------------------------------------------------------!
               if (v_ /= v_prime_ .or. j_ /= j_prime_ .or.                     &
                  abs(omega_ - omega_prime_) > 1) then
                  cycle
               endif
               !---------------------------------------------------------------!
               if (omega_ == omega_prime_) then
                  !------------------------------------------------------------!
                  ! Eq. 4 in "Coupling Matrix" section
                  !------------------------------------------------------------!
                  centrifugal_matrix_(channel_index_1_, channel_index_2_)      &
                     = calculate_diagonal_centrifugal_element(                 &
                     total_angular_momentum_, j_, omega_)
               else
                  !------------------------------------------------------------!
                  ! Eq. 5 in "Coupling Matrix" section
                  !------------------------------------------------------------!
                  centrifugal_matrix_(channel_index_1_, channel_index_2_)      &
                     = calculate_offdiagonal_centrifugal_element(              &
                     total_angular_momentum_, j_, omega_, omega_prime_,        &
                     delta_1_, delta_2_)
               endif
               !---------------------------------------------------------------!
            enddo
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(centrifugal_matrix_, 'u')
         !---------------------------------------------------------------------!
      end subroutine
!------------------------------------------------------------------------------!
      function delta_for_zero_omega(omega_) result(delta_)
         !! Checks if the input value equals 0; used in the calculation
         !! of off-diagonal elements of the centrifugal matrix; see
         !! Eq. 5 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: omega_
            !! input value which is to be compared with 0
         real(dp) :: delta_
            !! (output) 1 if omega_ = 0, 0 otherwise
         !---------------------------------------------------------------------!
         if (omega_ == 0) then
           delta_ = 1.0_dp
         else
           delta_ = 0.0_dp
         endif
         !---------------------------------------------------------------------!
      end function delta_for_zero_omega
!------------------------------------------------------------------------------!
      function calculate_diagonal_centrifugal_element(total_angular_momentum_, &
         j_, omega_) result(diagonal_element_)
         !! calculates diagonal element of the centrifgual matrix, see
         !! Eq. 4 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         real(dp) :: diagonal_element_
            !! (output) value of the diagonal element of the coupling matrix
         !---------------------------------------------------------------------!
         diagonal_element_ = real(                                             &
            total_angular_momentum_ * (total_angular_momentum_ + 1)            &
            + j_ * (j_ + 1) - 2 * omega_ **2, dp)
         !---------------------------------------------------------------------!
      end function calculate_diagonal_centrifugal_element
!------------------------------------------------------------------------------!
      function calculate_offdiagonal_centrifugal_element(total_angular_momentum_, &
         j_, omega_, omega_prime_, delta_1_, delta_2_) result(offdiagonal_element_)
         !! calculates off-diagonal element of the centrifgual matrix, see
         !! Eq. 5 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         integer(int32), intent(in) :: omega_prime_
            !! \\(\bar{\Omega}'\\)
         real(dp), intent(in) :: delta_1_, delta_2_
            !! Kronecker delta functions determined earlier
         real(dp) :: offdiagonal_element_
            !! (output) value of the off-diagonal element of the coupling matrix
         !---------------------------------------------------------------------!
         offdiagonal_element_ = - sqrt(real(                                   &
            (total_angular_momentum_ * (total_angular_momentum_ + 1)           &
            - omega_ * omega_prime_) * (j_ * (j_ + 1) - omega_ * omega_prime_),&
            dp)) * sqrt( (1.0_dp + delta_1_) * (1.0_dp + delta_2_) )
         !---------------------------------------------------------------------!
      end function calculate_offdiagonal_centrifugal_element
   !---------------------------------------------------------------------------!
   !                      Interaction potential matrix
   !---------------------------------------------------------------------------!
      subroutine pes_contribution(total_angular_momentum_,                     &
         intermolecular_distance_, channel_indices_, channels_omega_values_,   &
         nonzero_terms_per_element, nonzero_legendre_indices,                  &
         nonzero_coupling_coefficients, vmatrix)
         !! calculates the contribution to the coupling matrix
         !! from the the interaction potential (PES);
         !! see Eq. 1 in "Coupling Matrix" section;
         !! diagonal contribution from wavevectors (see the last term in
         !! Eq. 3 of "What are coupled equations" section) is added
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(in) :: intermolecular_distance_
            !! intermolecular distance
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: nonzero_terms_per_element(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix
         integer(int32), intent(in) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(in) :: nonzero_coupling_coefficients(:)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(out) :: vmatrix(:,:)
            !! (output) - the interaction potential contribution to the coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_coupling_coefficients,                &
            count_nonzero_coupling_matrix_elements,                            &
            count_nonzero_legendre_terms, channel_index_1_, channel_index_2_,  &
            omega_, omega_prime_
         !---------------------------------------------------------------------!
         vmatrix   = 0
         count_nonzero_coupling_coefficients    = 0
         count_nonzero_coupling_matrix_elements = 0
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, size(channel_indices_)
            omega_ = channels_omega_values_(channel_index_1_)
            do channel_index_2_ = 1, channel_index_1_
               omega_prime_ = channels_omega_values_(channel_index_2_)
               !---------------------------------------------------------------!
               if (omega_ /= omega_prime_) cycle
               !---------------------------------------------------------------!
               count_nonzero_coupling_matrix_elements                          &
                  = count_nonzero_coupling_matrix_elements + 1
               !---------------------------------------------------------------!
               ! process a single matrix element:
               ! get number of  non-zero terms in Legendre expansion for this
               ! matrix element
               !---------------------------------------------------------------!
               count_nonzero_legendre_terms                                    &
                  = nonzero_terms_per_element(count_nonzero_coupling_matrix_elements)
               !---------------------------------------------------------------!
               ! implementation of Eq. 1 in "Coupling Matrix" section
               !---------------------------------------------------------------!
               vmatrix(channel_index_1_, channel_index_2_)                     &
                  = calculate_single_pes_matrix_element(                       &
                     intermolecular_distance_, channel_index_1_,               &
                     channel_index_2_, channel_indices_,                       &
                     count_nonzero_legendre_terms,                             &
                     count_nonzero_coupling_coefficients,                      &
                     nonzero_legendre_indices, nonzero_coupling_coefficients)
               !---------------------------------------------------------------!
             enddo
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(vmatrix,'u')
         !---------------------------------------------------------------------!
      end subroutine pes_contribution
!------------------------------------------------------------------------------!
      function calculate_single_pes_matrix_element(intermolecular_distance_,   &
         channel_index_1_, channel_index_2_, channel_indices_,                 &
         count_nonzero_legendre_terms, count_nonzero_coupling_coefficients,    &
         nonzero_legendre_indices, nonzero_coupling_coefficients) result(matrix_element_)
         !! Implementation of Eq. 1 in "Coupling Matrix" section;
         !! diagonal contribution from wavevectors (see the last term in
         !! Eq. 3 of "What are coupled equations" section) is added
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: intermolecular_distance_
            !! intermolecular distance
         integer(int32), intent(in) :: channel_index_1_, channel_index_2_
            !! indices identifying matrix element
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: count_nonzero_legendre_terms
            !! number of non-zero terms in Legendre expansion for a single
            !! matrix element of the interaction potential
         integer(int32), intent(in) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(in) :: nonzero_coupling_coefficients(:)
            !! holds the values of the non-zero algebraic coefficients
         integer(int32), intent(inout) :: count_nonzero_coupling_coefficients
            !! running index counting non-zero coupling coefficients,
            !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) in the whole matrix;
            !! incremented in this subroutine
         real(dp) :: matrix_element_
            !! matirx element of the interaction potential contribution
            !! to the coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: v_, j_, v_prime_, j_prime_, lambda_, lambda_index_
         real(dp) :: internal_energy_, sum_over_lambda_,                       &
            algebraic_coefficient_, radial_term_
         !---------------------------------------------------------------------!
         v_ = v1array(channel_indices_(channel_index_1_))
         j_ = j1array(channel_indices_(channel_index_1_))
         v_prime_ = v1array(channel_indices_(channel_index_2_))
         j_prime_ = j1array(channel_indices_(channel_index_2_))
         internal_energy_ = elevel(channel_indices_(channel_index_1_))
         !---------------------------------------------------------------------!
         sum_over_lambda_ = 0.0_dp
         do lambda_index_ = 1, count_nonzero_legendre_terms
            !------------------------------------------------------------------!
            count_nonzero_coupling_coefficients                                &
               = count_nonzero_coupling_coefficients + 1
            !------------------------------------------------------------------!
            lambda_ = l1tab(nonzero_legendre_indices(count_nonzero_coupling_coefficients))
            algebraic_coefficient_ = nonzero_coupling_coefficients(count_nonzero_coupling_coefficients)
            !------------------------------------------------------------------!
            call get_radial_coupling_term_value(intermolecular_distance_,      &
               lambda_, v_, j_, v_prime_, j_prime_, radial_term_)
            !------------------------------------------------------------------!
            sum_over_lambda_ = sum_over_lambda_ + algebraic_coefficient_ * radial_term_
            !------------------------------------------------------------------!
         enddo
         matrix_element_ =  2.0_dp * reducedmass * sum_over_lambda_
         !---------------------------------------------------------------------!
         if (channel_index_1_ == channel_index_2_) then
            matrix_element_ = matrix_element_                                  &
               - wavenumber_squared_from_energy(internal_energy_)
         endif
         !---------------------------------------------------------------------!
      end function calculate_single_pes_matrix_element
!------------------------------------------------------------------------------! ---------------> clean the following 2 functions
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
      subroutine numerov(channel_indices_,channels_omega_values_,         &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients,nonzero_terms_per_element,    &
         nonzero_legendre_indices,nonzero_coupling_coefficients,nsteps,        &
         number_of_channels_,total_angular_momentum_,log_der_matrix)
         !! renormalized Numerov propagator
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_
            !! size of the basis
         integer(int32), intent(in) :: channel_indices_(number_of_channels_)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(number_of_channels_)
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
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(inout) :: log_der_matrix(:,:)
            !! resulting log-derivative matrix at RMAX  
         !---------------------------------------------------------------------!
         integer(int32) :: i, channel_index_1_, channel_index_2_
         real(dp) :: start, finish, r, step_numerov
         real(dp), allocatable, dimension(:,:) :: umatrix, r_temp, rmatrix,    &
            t_minus, t_center, t_plus, r_center, r_plus, w_tmp, v_tmp
         real(dp), dimension(number_of_channels_,number_of_channels_) ::         &
            cent_mat, v_mat, w_mat, t_mat, u_left, u_right
         !---------------------------------------------------------------------!
         CALL CPU_TIME(start)
         !---------------------------------------------------------------------!
         ! Calculate the centrifugal term
         !---------------------------------------------------------------------!
         call centrifugal_matrix(total_angular_momentum_,channel_indices_, &
            channels_omega_values_,cent_mat)
         step_numerov = (rmax - rmin)/dble(nsteps - 1)
         call allocate_2d(R_temp,number_of_channels_,number_of_channels_)
         !---------------------------------------------------------------------!
         ! Calculate the PES contribution at rmin
         !---------------------------------------------------------------------!
         call pes_contribution(total_angular_momentum_,rmin,                     &
            channel_indices_,channels_omega_values_,                      &
            nonzero_terms_per_element,nonzero_legendre_indices,                &
            nonzero_coupling_coefficients,v_mat)
         !---------------------------------------------------------------------!
         ! Coupling matrix W at rmin
         !---------------------------------------------------------------------!
         w_mat = v_mat + (1./rmin**2.)*cent_mat
         !---------------------------------------------------------------------!
         ! T-matrix (Eq. 6.23)
         !---------------------------------------------------------------------!
         t_mat = (step_numerov**2.0_dp)/12.0_dp * w_mat
         !---------------------------------------------------------------------!
         ! U-matrix (Eq. 6.25)
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, number_of_channels_
            do channel_index_2_ = 1, number_of_channels_
               u_left(channel_index_1_, channel_index_2_) = - t_mat(channel_index_1_, channel_index_2_)
            end do
            u_left(channel_index_1_, channel_index_1_) = 1.d0 + u_left(channel_index_1_, channel_index_1_)
         end do

         call invert_symmetric_matrix(u_left)
         call fill_symmetric_matrix(u_left, 'u')
         !---------------------------------------------------------------------!
         ! R-matrix at rmin + 1 = U-matrix at rmin
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, number_of_channels_
            do channel_index_2_ = 1, number_of_channels_
               r_temp(channel_index_1_, channel_index_2_) = 12.d0*u_left(channel_index_1_, channel_index_2_)
            enddo
            r_temp(channel_index_1_, channel_index_1_) = r_temp(channel_index_1_, channel_index_1_) - 10.d0
         enddo

         !---------------------------------------------------------------------!
         ! Continue the propagation to rmax
         !---------------------------------------------------------------------!
         do i=2,nsteps
            !------------------------------------------------------------------!
            ! Coupling matrix W at R
            !------------------------------------------------------------------!
            R = rmin + (i-1)*step_numerov
            call allocate_2d(rmatrix,number_of_channels_,number_of_channels_)
            call allocate_2d(umatrix,number_of_channels_,number_of_channels_)
            call pes_contribution(total_angular_momentum_,R,                     &
               channel_indices_,channels_omega_values_,                   &
               nonzero_terms_per_element,nonzero_legendre_indices,             &
               nonzero_coupling_coefficients,v_mat)
            w_mat = v_mat + (1./R**2.)*cent_mat
            !------------------------------------------------------------------!
            ! T-matrix at R (6.18)
            !------------------------------------------------------------------!
            t_mat = (step_numerov**2.0_dp)/12.0_dp * w_mat
            !------------------------------------------------------------------!
            ! U-matrix at R (6.20)
            !------------------------------------------------------------------!
            do channel_index_1_ = 1, number_of_channels_
               do channel_index_2_ = 1, number_of_channels_
                  u_left(channel_index_1_, channel_index_2_) = - t_mat(channel_index_1_, channel_index_2_)
               end do
               u_left(channel_index_1_, channel_index_1_) = 1.d0 + u_left(channel_index_1_, channel_index_1_)
            end do
            call invert_symmetric_matrix(u_left)
            call fill_symmetric_matrix(u_left, 'u')
            do channel_index_1_ = 1, number_of_channels_
               do channel_index_2_ = 1, number_of_channels_
                  umatrix(channel_index_1_, channel_index_2_) = 12.d0*u_left(channel_index_1_, channel_index_2_)
               enddo
               umatrix(channel_index_1_, channel_index_1_) = umatrix(channel_index_1_, channel_index_1_) - 10.d0
            enddo

            call invert_symmetric_matrix(r_temp)
            call fill_symmetric_matrix(r_temp, 'u')
            !------------------------------------------------------------------!
            ! Prepare T at Rmax - 1 and R at Rmax
            !------------------------------------------------------------------!
            if (i == nsteps - 1) then
               call allocate_2d(v_tmp,number_of_channels_,number_of_channels_)
               call allocate_2d(w_tmp,number_of_channels_,number_of_channels_)
               call allocate_2d(t_minus,number_of_channels_,number_of_channels_)

               call pes_contribution(total_angular_momentum_,r,                  &
                  channel_indices_,channels_omega_values_,                &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp + (1./R**2.)*cent_mat
               t_minus = (step_numerov**2.0_dp)/12.0_dp * w_tmp
               
               call allocate_2d(r_center,number_of_channels_,number_of_channels_)
               r_center = umatrix - r_temp
               !---------------------------------------------------------------!
               ! Prepare T at Rmax and R max + 1, and R at Rmax + 1
               !---------------------------------------------------------------!
            else if (i == nsteps) then
               call allocate_2d(v_tmp,number_of_channels_,number_of_channels_)
               call allocate_2d(w_tmp,number_of_channels_,number_of_channels_)
               call allocate_2d(t_center,number_of_channels_,number_of_channels_)

               call pes_contribution(total_angular_momentum_,r,                  &
                  channel_indices_,channels_omega_values_,                &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp + (1./R**2.)*cent_mat
               t_center = (step_numerov**2.0_dp)/12.0_dp * w_tmp
               
               call allocate_2d(r_plus,number_of_channels_,number_of_channels_)
               r_plus = umatrix - r_temp

               r = rmax + step_numerov
               call allocate_2d(v_tmp,number_of_channels_,number_of_channels_)
               call allocate_2d(t_plus,number_of_channels_,number_of_channels_)
               call allocate_2d(w_tmp,number_of_channels_,number_of_channels_)

               call pes_contribution(total_angular_momentum_,r,                  &
                  channel_indices_,channels_omega_values_,                &
                  nonzero_terms_per_element,nonzero_legendre_indices,          &
                  nonzero_coupling_coefficients,v_tmp)
               w_tmp = v_tmp + (1./R**2.)*cent_mat
               t_plus = (step_numerov**2.0_dp)/12.0_dp * w_tmp

            end if
            !------------------------------------------------------------------!
            ! R-matrix at R_{n+1} (Eq. 6.28)
            !------------------------------------------------------------------!
            rmatrix = umatrix - r_temp
            call allocate_2d(r_temp,number_of_channels_,number_of_channels_)
            !------------------------------------------------------------------!
            ! Move R_{n+1} to R_{n}
            !------------------------------------------------------------------!
            r_temp = rmatrix

         end do
         CALL CPU_TIME(finish)
         !---------------------------------------------------------------------!
         ! Eq. (6.29)
         !---------------------------------------------------------------------!
         call calculate_log_der_matrix(step_numerov,number_of_channels_,        &
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
