module propagator_mod
   !! This modules contains the subroutines used by the renormalized
   !! Numerov propagator.                                   
   !----------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use data_mod
   use io_mod
   use array_operations_mod, only: invert_symmetric_matrix, fill_symmetric_matrix, add_scalar_to_diagonal
   use centrifugal_matrix_mod, only: calculate_centrifugal_matrix
   use pes_matrix_mod, only: calculate_pes_matrix
   use utility_functions_mod, only: time_count_summary
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: numerov
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
      subroutine numerov(number_of_channels_, channel_indices_,                &
         channels_omega_values_, nonzero_terms_per_element_,                   &
         nonzero_legendre_indices_, nonzero_algebraic_coefficients_,           &
         number_of_steps_, total_angular_momentum_, log_der_matrix_)
         !! renormalized Numerov propagator
         !! ...
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_
            !! size of the basis
         integer(int32), intent(in) :: channel_indices_(number_of_channels_)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(number_of_channels_)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: nonzero_terms_per_element_(:)
            !! keeps the number of non-zero terms in the sum (Eq. (6.21)) for
            !! each non-zero element of W/V
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds the proper indices pointing to l1/l2/lltabs, which
            !! correspond to the non-vanishing elements of the sum  (Eq. (6.21))
            !! for each non-zero element of W/V
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         integer(int32), intent(in) :: number_of_steps_
            !! number of steps from rmin to rmax
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(inout) :: log_der_matrix_(number_of_channels_,number_of_channels_)
            !! resulting log-derivative matrix at RMAX  
         !---------------------------------------------------------------------!
         integer(int32) :: i, channel_index_1_, channel_index_2_
         real(dp) :: start, finish, intermolecular_distance_, step_numerov_, calculation_time_
         real(dp), dimension(number_of_channels_,number_of_channels_) ::       &
            centrifugal_matrix_,  &
            t_matrix_minus_, t_matrix_, t_matrix_plus_, r_matrix_,             &
            r_matrix_rmax_, r_matrix_plus_
         !---------------------------------------------------------------------!
         CALL CPU_TIME(start)
         step_numerov_ = (rmax - rmin)/dble(number_of_steps_ - 1)
         intermolecular_distance_ = rmin
         !---------------------------------------------------------------------!
         ! Initial setup: calculate centrifugal matrix and R_matrix at Rmin + 1
         !---------------------------------------------------------------------!
         call initial_setup(number_of_channels_, step_numerov_, total_angular_momentum_,      &
            intermolecular_distance_, channel_indices_, channels_omega_values_,&
            nonzero_terms_per_element_, nonzero_legendre_indices_,             &
            nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_)
         !---------------------------------------------------------------------!
         ! Propagation loop
         !---------------------------------------------------------------------!
         do i=2, number_of_steps_ - 2
            !------------------------------------------------------------------!
            intermolecular_distance_ = rmin + (i-1)*step_numerov_
            !------------------------------------------------------------------!
            call general_propagation_step(number_of_channels_, step_numerov_, &
                  total_angular_momentum_, intermolecular_distance_, channel_indices_,    &
                  channels_omega_values_, nonzero_terms_per_element_, nonzero_legendre_indices_, &
                  nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_)
            !------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------!
         call handle_final_propagation_steps(number_of_channels_, step_numerov_, &
                  total_angular_momentum_, channel_indices_,    &
                  channels_omega_values_, nonzero_terms_per_element_, nonzero_legendre_indices_, &
                  nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_, &
                  t_matrix_minus_, t_matrix_, t_matrix_plus_, r_matrix_rmax_, r_matrix_plus_)
         !---------------------------------------------------------------------!
         CALL CPU_TIME(finish)
         !---------------------------------------------------------------------!
         ! Eq. (6.29)
         !---------------------------------------------------------------------!
         call calculate_log_der_matrix(step_numerov_,number_of_channels_,        &
            t_matrix_minus_,t_matrix_,t_matrix_plus_,r_matrix_rmax_,r_matrix_plus_,log_der_matrix_)
         !---------------------------------------------------------------------!
         call propagator_summary(rmin, rmax, number_of_steps_)
         !---------------------------------------------------------------------!
         if (prntlvl.ge.2) then
            call time_count_summary(start, finish, calculation_time_,          &
               "Propagation completed in ")
         endif
         !---------------------------------------------------------------------!
      end subroutine numerov
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine initial_setup(number_of_channels_, step_numerov_, total_angular_momentum_,   &
         intermolecular_distance_, channel_indices_, channels_omega_values_,   &
         nonzero_terms_per_element_, nonzero_legendre_indices_,                &
         nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_)
         !! ...
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_
            !! size of the basis
         real(dp), intent(in) :: step_numerov_
            !! step of the propagator
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(in) :: intermolecular_distance_
            !! intermolecular distance
         integer(int32), intent(in) :: channel_indices_(number_of_channels_)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(number_of_channels_)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: nonzero_terms_per_element_(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(inout) :: centrifugal_matrix_(number_of_channels_,number_of_channels_)
            !! (R**2)*centrifugal matrix - calculated once, will be used throughout the propagation
         real(dp), intent(inout) :: r_matrix_(number_of_channels_,number_of_channels_)
            !! R-matrix at Rmin
         !---------------------------------------------------------------------!   
         real(dp), dimension(number_of_channels_,number_of_channels_) ::       &
            pes_matrix_, coupling_matrix_, t_matrix_, u_matrix_
         !---------------------------------------------------------------------!   
         ! Calculate centrifugal matrix
         !---------------------------------------------------------------------!   
         call calculate_centrifugal_matrix(total_angular_momentum_,            &
            channel_indices_, channels_omega_values_, centrifugal_matrix_)
         !---------------------------------------------------------------------!   
         ! Calculate PES matrix at rmin
         !---------------------------------------------------------------------!   
         call calculate_pes_matrix(total_angular_momentum_,                    &
            intermolecular_distance_, channel_indices_,                        &
            channels_omega_values_, nonzero_terms_per_element_,                &
            nonzero_legendre_indices_, nonzero_algebraic_coefficients_, pes_matrix_) 
         !---------------------------------------------------------------------!   
         ! Merge centrifugal and PES matrix into Coupling matrix
         !---------------------------------------------------------------------! 
         call calculate_coupling_matrix(intermolecular_distance_, pes_matrix_, &
            centrifugal_matrix_, coupling_matrix_)
         !---------------------------------------------------------------------!
         ! Calculate initial T-matrix and U-matrix
         !---------------------------------------------------------------------!
         call calculate_t_matrix(step_numerov_, coupling_matrix_, t_matrix_)
         call calculate_u_matrix(t_matrix_, u_matrix_)
         !---------------------------------------------------------------------!
         ! Initialize R-matrix: R-matrix at rmin + 1 = U-matrix at rmin
         !---------------------------------------------------------------------!
         r_matrix_ = u_matrix_
         !---------------------------------------------------------------------!
      end subroutine initial_setup
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine general_propagation_step(number_of_channels_, step_numerov_, total_angular_momentum_,   &
         intermolecular_distance_, channel_indices_, channels_omega_values_,   &
         nonzero_terms_per_element_, nonzero_legendre_indices_,                &
         nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_, is_t_matrix_required_, t_matrix_returned_)
         !! ...
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_
            !! size of the basis
         real(dp), intent(in) :: step_numerov_
            !! step of the propagator
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(in) :: intermolecular_distance_
            !! intermolecular distance
         integer(int32), intent(in) :: channel_indices_(number_of_channels_)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(number_of_channels_)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: nonzero_terms_per_element_(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix;
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(in) :: centrifugal_matrix_(number_of_channels_,number_of_channels_)
            !! (R**2)*centrifugal matrix -
         real(dp), intent(inout) :: r_matrix_(number_of_channels_,number_of_channels_)
            !! on input: R-matrix at previous step
            !! on output: R-matrix at next step
         logical, intent(in), optional :: is_t_matrix_required_
            !! ...
         real(dp), intent(out), optional :: t_matrix_returned_(number_of_channels_,number_of_channels_)
            !! on input: R-matrix at previous step
            !! on output: R-matrix at next step
         !---------------------------------------------------------------------!   
         real(dp), dimension(number_of_channels_,number_of_channels_) ::       &
            pes_matrix_, coupling_matrix_, t_matrix_, u_matrix_, r_matrix_plus_
         !---------------------------------------------------------------------!   
         ! Calculate PES matrix at R
         !---------------------------------------------------------------------!   
         call calculate_pes_matrix(total_angular_momentum_,                    &
            intermolecular_distance_, channel_indices_,                        &
            channels_omega_values_, nonzero_terms_per_element_,                &
            nonzero_legendre_indices_, nonzero_algebraic_coefficients_, pes_matrix_) 
         !---------------------------------------------------------------------!   
         ! Merge centrifugal and PES matrix into Coupling matrix
         !---------------------------------------------------------------------! 
         call calculate_coupling_matrix(intermolecular_distance_, pes_matrix_, &
            centrifugal_matrix_, coupling_matrix_)
         !---------------------------------------------------------------------!
         ! Calculate T-matrix and U-matrix
         !---------------------------------------------------------------------!
         call calculate_t_matrix(step_numerov_, coupling_matrix_, t_matrix_)
         call calculate_u_matrix(t_matrix_, u_matrix_)
         !---------------------------------------------------------------------!
         ! Invert R matrix from previous step
         !---------------------------------------------------------------------!
         call invert_symmetric_matrix(r_matrix_)
         call fill_symmetric_matrix(r_matrix_, 'u')
         !------------------------------------------------------------------!
         ! R_{n+1} = U_{n} - R_{n}^{-1}
         !------------------------------------------------------------------!
         r_matrix_plus_ = u_matrix_ - r_matrix_
         !------------------------------------------------------------------!
         ! Move R_{n+1} to R_{n}
         !------------------------------------------------------------------!
         r_matrix_ = r_matrix_plus_
         !------------------------------------------------------------------!
         if (present(is_t_matrix_required_) .and. is_t_matrix_required_) then
            t_matrix_returned_ = t_matrix_
         endif
      !---------------------------------------------------------------------!
      end subroutine general_propagation_step
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine handle_final_propagation_steps(number_of_channels_,          &
         step_numerov_, total_angular_momentum_,   &
         channel_indices_, channels_omega_values_,   &
         nonzero_terms_per_element_, nonzero_legendre_indices_,                &
         nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_,      &
         t_matrix_minus_, t_matrix_, t_matrix_plus_, r_matrix_rmax_, r_matrix_plus_)
         !! Handles propagation at the last two grid points:
         !! R_{N-1} and R_{N}: provides T-matrix at N-1, N and N+1 points
         !! and the Ratio matrix at N and N+1 points
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_
            !! size of the basis
         real(dp), intent(in) :: step_numerov_
            !! step of the propagator
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: channel_indices_(number_of_channels_)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(number_of_channels_)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: nonzero_terms_per_element_(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix;
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(in) :: centrifugal_matrix_(number_of_channels_,number_of_channels_)
            !! \\(R^{2} \cdot\\) centrifugal matrix -
         real(dp), intent(inout) :: r_matrix_(number_of_channels_,number_of_channels_)
            !! Ratio matrix at N-1 step
         real(dp), intent(out) :: t_matrix_minus_(number_of_channels_,number_of_channels_)
            !! T-matrix at N-1 step
         real(dp), intent(out) :: t_matrix_(number_of_channels_,number_of_channels_)
            !! T-matrix at N step
         real(dp), intent(out) :: t_matrix_plus_(number_of_channels_,number_of_channels_)
            !! T-matrix at N+1 step
         real(dp), intent(out) :: r_matrix_rmax_(number_of_channels_,number_of_channels_)
            !! Ratio matrix at N step
         real(dp), intent(out) :: r_matrix_plus_(number_of_channels_,number_of_channels_)
            !! Ratio matrix at N+1 step
         !---------------------------------------------------------------------!
         logical :: is_t_matrix_required_
         real(dp) :: intermolecular_distance_
         real(dp), dimension(number_of_channels_,number_of_channels_) ::       &
            pes_matrix_, coupling_matrix_, u_matrix_
         !---------------------------------------------------------------------!
         is_t_matrix_required_ = .true.
         !---------------------------------------------------------------------!
         ! N - 1 step
         !---------------------------------------------------------------------!
         intermolecular_distance_ = rmax - step_numerov_
         call general_propagation_step(number_of_channels_, step_numerov_,     &
            total_angular_momentum_, intermolecular_distance_,                 &
            channel_indices_, channels_omega_values_,                          &
            nonzero_terms_per_element_, nonzero_legendre_indices_,             &
            nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_,   &
            is_t_matrix_required_, t_matrix_minus_)
         r_matrix_rmax_ = r_matrix_
         !---------------------------------------------------------------------!
         ! N  step
         !---------------------------------------------------------------------!
         intermolecular_distance_ = intermolecular_distance_ + step_numerov_
         call general_propagation_step(number_of_channels_, step_numerov_,     &
            total_angular_momentum_, intermolecular_distance_,                 &
            channel_indices_, channels_omega_values_,                          &
            nonzero_terms_per_element_, nonzero_legendre_indices_,             &
            nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_,   &
            is_t_matrix_required_, t_matrix_)
         r_matrix_plus_ = r_matrix_
         !---------------------------------------------------------------------!
         ! N + 1 step
         !---------------------------------------------------------------------!
         intermolecular_distance_ = intermolecular_distance_ + step_numerov_
         call general_propagation_step(number_of_channels_, step_numerov_,     &
            total_angular_momentum_, intermolecular_distance_,                 &
            channel_indices_, channels_omega_values_,                          &
            nonzero_terms_per_element_, nonzero_legendre_indices_,             &
            nonzero_algebraic_coefficients_, centrifugal_matrix_, r_matrix_,   &
            is_t_matrix_required_, t_matrix_plus_)
         !---------------------------------------------------------------------!
      end subroutine handle_final_propagation_steps
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_coupling_matrix(intermolecular_distance_,           &
         pes_matrix_, centrifugal_matrix_, coupling_matrix_)
         !! Combines the contribution from the interaction potential, total and
         !! and internal energy (pes_matrix_) with centrifugal matrix
         !! \\( W_{\mathrm{N}} = V_{\mathrm{N}} + 1/R^{2} L^{2} \\)
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: intermolecular_distance_
            !! intermolecular distance
         real(dp), intent(in) :: pes_matrix_(:,:)
            !! holds contribution from the interaction potential, total and
            !! and internal energy
         real(dp), intent(in) :: centrifugal_matrix_(:,:)
            !! \\(R^{2}\\) centrifugal matrix
         real(dp), intent(inout) :: coupling_matrix_(:,:)
            !! (output) Coupling (W) matrix
         !---------------------------------------------------------------------!
         coupling_matrix_ = pes_matrix_                                        &
            + (1.0_dp/intermolecular_distance_**2.0_dp) * centrifugal_matrix_
         !---------------------------------------------------------------------!
      end subroutine calculate_coupling_matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_t_matrix(step_, coupling_matrix_, t_matrix_)
         !! Calculates the T-matrix from the coupling matrix at grid point N:
         !! \\( T_{\mathrm{N}} = h^{2}/12 W_{\mathrm{N}} \\)
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: step_
            !! step of the propagator
         real(dp), intent(in) :: coupling_matrix_(:,:)
            !! Coupling (W) matrix v
         real(dp), intent(inout) :: t_matrix_(:,:)
            !! (output) T-matrix at grid point N
         !---------------------------------------------------------------------!
         t_matrix_ = (step_**2.0_dp)/12.0_dp * coupling_matrix_
         !---------------------------------------------------------------------!
      end subroutine calculate_t_matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_u_matrix(t_matrix_, u_matrix_)
         !! Calculates the U-matrix from T-matrix at grid point N:
         !! \\(U_{\mathrm{N}} = 12(\mathbf{I} - 10 T_{\mathrm{N}})^{-1} - 10 \mathbf{I}\\)
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: t_matrix_(:,:)
            !! T-matrix at grid point N
         real(dp), intent(inout) :: u_matrix_(:,:)
            !! (output) U-matrix at grid point N
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         !---------------------------------------------------------------------!
         u_matrix_ = - t_matrix_
         call add_scalar_to_diagonal(u_matrix_, 1.0_dp)
         call invert_symmetric_matrix(u_matrix_)
         call fill_symmetric_matrix(u_matrix_, 'u')
         u_matrix_ = 12.0_dp * u_matrix_
         call add_scalar_to_diagonal(u_matrix_, - 10.0_dp)
         !---------------------------------------------------------------------!
      end subroutine calculate_u_matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_log_der_matrix(step_, number_of_channels_,          &
         t_matrix_minus_, t_matrix_, t_matrix_plus_, r_matrix_, r_matrix_plus_,&
         log_der_matrix_)
         !! calculates the log-derivative matrix from
         !! \begin{equation}
         !! {Y}_{\rm N} = \frac{1}{h} \Biggl(\Bigl(\frac{1}{2}\mathbf{I}-{T}_{\rm{N}+1}\Bigr)
         !! \Bigl(\mathbf{I}-{T}_{\rm{N}+1}\Bigr)^{-1} {R}_{\rm{N}+1} -
         !! \Bigl(\frac{1}{2}\mathbf{I}-{T}_{\rm{N}-1}\Bigr)
         !! \Bigl(\mathbf{I}-\mathbf{T}_{\rm{N}-1}\Bigr)^{-1}\mathbf{R}_{\rm{N}}^{-1}
         !! \Biggr)\Bigl(\mathbf{I}-{T}_{\rm{N}}\Bigr) 
         !! \end{equation}
         !! called by numerov at the end of the propagation
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: step_
            !! propagator step
         integer(int32), intent(in) :: number_of_channels_
            !! number of scattering channels in the block
         real(dp), intent(in) :: t_matrix_minus_(number_of_channels_, number_of_channels_)
            !! T-matrix at R_{max - 1}
         real(dp), intent(in) :: t_matrix_(number_of_channels_, number_of_channels_)
            !! T-matrix at R_{max}
         real(dp), intent(in) :: t_matrix_plus_(number_of_channels_, number_of_channels_)
            !! T-matrix at R_{max + 1} 
         real(dp), intent(in) :: r_matrix_(number_of_channels_, number_of_channels_)
            !! R-matrix at R_{max} 
         real(dp), intent(in) :: r_matrix_plus_(number_of_channels_, number_of_channels_)
            !! R-matrix at R_{max + 1}
         real(dp), intent(inout) :: log_der_matrix_(number_of_channels_, number_of_channels_)
            !! log-derivative matrix
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_1_
         real(dp) :: matrix_a_(number_of_channels_,number_of_channels_),       &
            matrix_b_(number_of_channels_,number_of_channels_),                &
            matrix_c_(number_of_channels_,number_of_channels_),                &
            matrix_d_(number_of_channels_,number_of_channels_),                &
            matrix_e_(number_of_channels_,number_of_channels_),                &
            matrix_ab_(number_of_channels_,number_of_channels_),               &
            matrix_cd_(number_of_channels_,number_of_channels_),               &
            left_matrix_(number_of_channels_,number_of_channels_),             &
            right_matrix_(number_of_channels_,number_of_channels_),            &
            working_r_matrix_(number_of_channels_,number_of_channels_),        &
            matrix_difference_(number_of_channels_,number_of_channels_)
         !---------------------------------------------------------------------!
         log_der_matrix_ = 0
         !---------------------------------------------------------------------!
         ! First bracket: (1/2 I - T_{N+1})
         !---------------------------------------------------------------------!
         matrix_a_ = - t_matrix_plus_
         call add_scalar_to_diagonal(matrix_a_, 0.5_dp)
         !---------------------------------------------------------------------!
         ! Second bracket: (I - T_{N+1})^{-1}
         !---------------------------------------------------------------------!
         matrix_b_ = - t_matrix_plus_
         call add_scalar_to_diagonal(matrix_b_, 1.0_dp)
         call invert_symmetric_matrix(matrix_b_)
         call fill_symmetric_matrix(matrix_b_, "u")
         !---------------------------------------------------------------------!
         ! Third bracket: (1/2 I - T_{N-1})
         !---------------------------------------------------------------------!
         matrix_c_ = - t_matrix_minus_
         call add_scalar_to_diagonal(matrix_c_, 0.5_dp)
         !---------------------------------------------------------------------!
         ! Fourth bracket: (I - T_{N-1})^{-1}
         !---------------------------------------------------------------------!
         matrix_d_ = - t_matrix_minus_
         call add_scalar_to_diagonal(matrix_d_, 1.0_dp)
         call invert_symmetric_matrix(matrix_d_)
         call fill_symmetric_matrix(matrix_d_, "u")
         !---------------------------------------------------------------------!
         ! Last bracket: (I - T_{N})
         !---------------------------------------------------------------------!
         matrix_e_ = - t_matrix_
         call add_scalar_to_diagonal(matrix_e_, 1.0_dp)
         !---------------------------------------------------------------------!
         ! Copy R_{N} to another matrix (r_matrix_ is protected as intent(in))
         !---------------------------------------------------------------------!
         working_r_matrix_ = r_matrix_
         call invert_symmetric_matrix(working_r_matrix_)
         call fill_symmetric_matrix(working_r_matrix_, "u")
         !---------------------------------------------------------------------!
         ! The first term in the large bracket:
         ! (1/2 I - T_{N+1}) (I - T_{N+1})^{-1} R_{N+1}
         !---------------------------------------------------------------------!
         CALL DGEMM('N','N',number_of_channels_,number_of_channels_,           &
            number_of_channels_,1.0_dp,matrix_a_,number_of_channels_,matrix_b_,&
            number_of_channels_,0.0_dp,matrix_ab_,number_of_channels_)
         CALL DGEMM('N','N',number_of_channels_,number_of_channels_,           &
            number_of_channels_,1.0_dp,matrix_ab_,number_of_channels_,         &
            r_matrix_plus_,number_of_channels_,0.0_dp,left_matrix_,number_of_channels_)
         !---------------------------------------------------------------------!
         ! The second term in the large bracket:
         ! (1/2 I - T_{N-1}) (I - T_{N-1})^{-1} R_{N}^{-1}
         !---------------------------------------------------------------------!
         CALL DGEMM('N','N',number_of_channels_,number_of_channels_,           &
            number_of_channels_,1.0_dp,matrix_c_,number_of_channels_,matrix_d_,&
            number_of_channels_,0.0_dp,matrix_cd_,number_of_channels_)
         CALL DGEMM('N','N',number_of_channels_,number_of_channels_,           &
            number_of_channels_,1.0_dp,matrix_cd_,number_of_channels_,         &
            working_r_matrix_,number_of_channels_,0.0_dp,right_matrix_,number_of_channels_)
         !---------------------------------------------------------------------!
         ! Substract the two terms in the large bracket
         !---------------------------------------------------------------------!
         matrix_difference_ = left_matrix_ - right_matrix_
         !---------------------------------------------------------------------!
         ! Multiply the product by (I - T_{N})
         !---------------------------------------------------------------------!
         CALL DGEMM('N','N',number_of_channels_,number_of_channels_,           &
            number_of_channels_,1.0_dp/step_,matrix_difference_,               &
            number_of_channels_,matrix_e_,number_of_channels_,0.0_dp,          &
            log_der_matrix_,number_of_channels_)
         !---------------------------------------------------------------------!
      end subroutine calculate_log_der_matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine propagator_summary(r_min_, r_max_, number_of_steps_)
         !! Print a simple message after the propagation is finished
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: r_min_
            !! initial point on the \\(R\\) grid
         real(dp), intent(in) :: r_max_
            !! final point on the \\(R\\) grid
         integer(int32), intent(in) :: number_of_steps_
            !! number of steps on the \\(R\\) grid
         !---------------------------------------------------------------------!
         call write_message("Coupled equations were solved from " //           &
            trim(adjustl(float_to_character(r_min_, "(F10.4)")))// " a.u. to " &
            // trim(adjustl(float_to_character(r_max_, "(F10.4)")))//          &
            " a.u. in "// trim(adjustl(integer_to_character(number_of_steps_)))&
            // " steps ")
         call write_message("(constant dr = " //                               &
            trim(adjustl(float_to_character((rmax - rmin) /                    &
            real(number_of_steps_ - 1, dp), "(E14.8)"))) // " a.u.)")
         !---------------------------------------------------------------------!
      end subroutine propagator_summary
!------------------------------------------------------------------------------!
end module propagator_mod
