module pes_matrix_mod
   !! This module provides functions calculating the algebraic coefficients
   !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) entering the PES matrix,
   !! and the full PES matrix (see Eq. 1 in the "Coupling Matrix" section).
   !---------------------------------------------------------------------------!
   !! Subroutines preparing algebraic coefficients and determining the number
   !! of non-zero terms of the PES matrix are called only once,
   !! before Numerov propagator is initialized
   !---------------------------------------------------------------------------!
   !! Subroutines calculating the full PES matrix at desired R are called
   !! within the propagator loop
   !---------------------------------------------------------------------------!
   !! The data is organized as follows:
   !! - number of non-zero terms of the PES matrix due to
   !!   \\(\bar{\Omega} = \bar{\Omega}'\\) condition is saved as
   !!   "number_of_nonzero_pes_matrix_elements"
   !! - number of non-vanishing terms in the sum over \\(\lambda\\)
   !!   in Eq. 1 in the "Coupling Matrix" section is saved in a
   !!   "nonzero_terms_per_element" array which is of
   !!   "number_of_nonzero_pes_matrix_elements" size
   !! - a _total_ number of non-vanishing  \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\)
   !!   coefficients is saved in "number_of_nonzero_algebraic_coefficients"
   !! - _all_ non-vanishing  \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\)
   !!   coefficients are saved in "nonzero_algebraic_coefficients" array,
   !!   which is of "number_of_nonzero_algebraic_coefficients" size
   !! - corresponding \\(\lambda\\) value for each non-vanishing coefficient
   !!   is saved as an index to "l1tab" in the "nonzero_legendre_indices"
   !!   array (which is of "number_of_nonzero_algebraic_coefficients" size)
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_message,                &
      integer_to_character, time_count_summary
   use data_mod
   use io_mod
   use array_operations_mod, only: allocate_1d, fill_symmetric_matrix
   use math_functions_mod, only: percival_seaton_coefficient,                  &
      triangle_inequality_holds, is_sum_even, zero_projections_3j_condition
   use radial_coupling_terms_mod, only: get_radial_coupling_term_value
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: initialize_pes_matrix, calculate_pes_matrix
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
   ! Subroutines preparing algebraic coefficients and determining the number
   ! of non-zero terms of the PES matrix; these are called only once,
   ! before Numerov propagator is initialized
   !---------------------------------------------------------------------------!
      subroutine initialize_pes_matrix(channel_indices,                        &
         channels_omega_values, nonzero_terms_per_element,                     &
         nonzero_legendre_indices, nonzero_algebraic_coefficients)
         !! launches "check_nonzero_pes_matrix_elements" and
         !! "prepare_pes_matrix_elements" subroutines; called from the main program
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(inout), allocatable :: nonzero_terms_per_element(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix
         integer(int32), intent(inout), allocatable :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix;
         real(dp), intent(inout), allocatable :: nonzero_algebraic_coefficients(:)
            !! holds the values of the non-zero algebraic coefficients
         !---------------------------------------------------------------------!
         integer(int32) :: number_of_channels,                                 &
            number_of_nonzero_pes_matrix_elements,                             &
            number_of_nonzero_algebraic_coefficients
         real(dp) :: time_start_, time_stop_, total_time_
         !---------------------------------------------------------------------!
         call cpu_time(time_start_)
         !---------------------------------------------------------------------!
         call check_nonzero_pes_matrix_elements(channel_indices,               &
            channels_omega_values, number_of_nonzero_pes_matrix_elements,      &
            number_of_nonzero_algebraic_coefficients)
         !---------------------------------------------------------------------!
         call allocate_1d(nonzero_terms_per_element,number_of_nonzero_pes_matrix_elements)
         call allocate_1d(nonzero_algebraic_coefficients,number_of_nonzero_algebraic_coefficients)
         call allocate_1d(nonzero_legendre_indices,number_of_nonzero_algebraic_coefficients)
         !---------------------------------------------------------------------!
         call prepare_pes_matrix_elements(channel_indices,         &
            channels_omega_values, nonzero_terms_per_element,                  &
            nonzero_legendre_indices, nonzero_algebraic_coefficients)
         !---------------------------------------------------------------------!
         call cpu_time(time_stop_)
         !---------------------------------------------------------------------!
         if (prntlvl.ge.2) then
            !------------------------------------------------------------------!
            number_of_channels = size(channel_indices)
            call print_pes_matrix_elements_summary(                            &
               number_of_channels, number_of_nonzero_pes_matrix_elements,      &
               number_of_nonzero_algebraic_coefficients)
            !------------------------------------------------------------------!
            call time_count_summary(time_start_, time_stop_, total_time_,      &
               "-- PES matrix initialization completed in ")
            !------------------------------------------------------------------!
         endif
         !---------------------------------------------------------------------!
      end subroutine initialize_pes_matrix
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine check_nonzero_pes_matrix_elements(channel_indices,&
         channels_omega_values, number_of_nonzero_pes_matrix_elements,    &
         number_of_nonzero_algebraic_coefficients)
         !! checks the number of non-zero PES matrix elements due to
         !! the \bar{\Omega} = \bar{\Omega}' condition,
         !! "number_of_nonzero_pes_matrix_elements",
         !! and the total number of non-zero algebraic coefficients,
         !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\), in the whole matrix,
         !! "number_of_nonzero_algebraic_coefficients".
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(out) :: number_of_nonzero_pes_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element of the PES matrix
         integer(int32), intent(out) :: number_of_nonzero_algebraic_coefficients
            !! number of all non-zero algberaix coefficients in the whole PES matrix
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_pes_matrix_elements,             &
            count_nonzero_algebraic_coefficients, j_, j_prime_, omega_,         &
            omega_prime_, lambda_, channel_index_1_, channel_index_2_, legendre_term_index_
         !---------------------------------------------------------------------!
         count_nonzero_algebraic_coefficients = 0
         count_nonzero_pes_matrix_elements = 0
         do channel_index_1_ = 1, size(channel_indices)
            j_ = j1array(channel_indices(channel_index_1_))
            omega_ = channels_omega_values(channel_index_1_)
            do channel_index_2_ = 1, channel_index_1_
               j_prime_ = j1array(channel_indices(channel_index_2_))
               omega_prime_ = channels_omega_values(channel_index_2_)
               !---------------------------------------------------------------!
               if (omega_ /= omega_prime_) cycle
               !---------------------------------------------------------------!
               count_nonzero_pes_matrix_elements =                        &
                  count_nonzero_pes_matrix_elements + 1
               do legendre_term_index_ = 1, nterms
                  lambda_ = l1tab(legendre_term_index_)
                  if (.not. zero_projections_3j_condition(j_, j_prime_, lambda_)) cycle
                  count_nonzero_algebraic_coefficients =                        &
                     count_nonzero_algebraic_coefficients+1            
              enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
         number_of_nonzero_algebraic_coefficients = count_nonzero_algebraic_coefficients
         number_of_nonzero_pes_matrix_elements = count_nonzero_pes_matrix_elements
         !---------------------------------------------------------------------!
      end subroutine check_nonzero_pes_matrix_elements
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine prepare_pes_matrix_elements(channel_indices,      &
         channels_omega_values, nonzero_terms_per_element,                     &
         nonzero_legendre_indices, nonzero_algebraic_coefficients)
         !! prepares:
         !! -- nonzero_terms_per_element - number of non-vanishing terms in
         !!    the sum over \\(\lambda\\) in Eq. 1 in the "Coupling Matrix" section
         !! -- nonzero_legendre_indices - corresponding \\(\lambda\\) value for
         !!    each non-vanishing coefficient is saved as an index to "l1tab"
         !! -- nonzero_algebraic_coefficients --  holds _all_ non-vanishing
         !!    \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) coefficients
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(inout) :: nonzero_terms_per_element(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix
         integer(int32), intent(inout) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix;
         real(dp), intent(inout) :: nonzero_algebraic_coefficients(:)
            !! holds the values of the non-zero algebraic coefficients
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_pes_matrix_elements,             &
            count_nonzero_algebraic_coefficients, count_nonzero_legendre_terms, &
            j_, j_prime_, omega_, omega_prime_, lambda_, channel_index_1_,     &
            channel_index_2_, legendre_term_index_
         real(dp) :: pscoeff
         !---------------------------------------------------------------------!
         nonzero_terms_per_element        = 0
         nonzero_legendre_indices         = 0
         nonzero_algebraic_coefficients    = 0
         count_nonzero_algebraic_coefficients     = 0
         count_nonzero_pes_matrix_elements  = 0
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, size(channel_indices)
            j_     = j1array(channel_indices(channel_index_1_))
            omega_ = channels_omega_values(channel_index_1_)
            do channel_index_2_ = 1, channel_index_1_
               j_prime_     = j1array(channel_indices(channel_index_2_))
               omega_prime_ = channels_omega_values(channel_index_2_)
               if (omega_  /= omega_prime_) cycle
               !---------------------------------------------------------------!
               ! passed \bar{\Omega} = \bar{\Omega}' condition
               !---------------------------------------------------------------!
               count_nonzero_pes_matrix_elements =                        &
                  count_nonzero_pes_matrix_elements + 1
               !---------------------------------------------------------------!
               ! process a single matrix element:
               ! determine non-zero terms in the sum over legendre polynomials
               ! for this element; these are saved to ...
               !---------------------------------------------------------------!
               call process_single_matrix_element(j_, j_prime_, omega_,        &
                  count_nonzero_algebraic_coefficients,                         &
                  count_nonzero_legendre_terms, nonzero_legendre_indices,      &
                  nonzero_algebraic_coefficients)
               
               nonzero_terms_per_element(count_nonzero_pes_matrix_elements)&
                  = count_nonzero_legendre_terms
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine prepare_pes_matrix_elements
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine process_single_matrix_element(j_, j_prime_, omega_,           &
         count_nonzero_algebraic_coefficients, count_nonzero_legendre_terms,    &
         nonzero_legendre_indices, nonzero_algebraic_coefficients)
         !! calculates the non-zero algebraic coefficients
         !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) for a single matrix
         !! element - see Eq. (1) in the "Coupling matrix" section;
         !! algebraic coefficients are saved to nonzero_algebraic_coefficients
         !! array; corresponding indices to l1tab are saved to
         !! nonzero_legendre_indices array
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: j_
            !! pre-collisional rotational angular momentum
         integer(int32), intent(in) :: j_prime_
            !! post-collisional rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         integer(int32), intent(inout) :: count_nonzero_algebraic_coefficients
            !! running index counting non-zero coupling coefficients,
            !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) in the whole matrix;
            !! incremented in this subroutine
         integer(int32), intent(inout) :: count_nonzero_legendre_terms
            !! number of non-zero terms in Legendre expansion for a single
            !! matrix element of the interaction potential
         integer(int32), intent(inout) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the PES matrix;
         real(dp), intent(inout) :: nonzero_algebraic_coefficients(:)
            !! holds values of the non-zero algebraic coefficients
         !---------------------------------------------------------------------!
         integer(int32) :: legendre_term_index_, lambda_
         real(dp) :: pscoeff
         !---------------------------------------------------------------------!
         count_nonzero_legendre_terms = 0
         do legendre_term_index_ = 1, nterms
            lambda_ = l1tab(legendre_term_index_)
            !------------------------------------------------------------!
            ! check the condition on 3-j symbol with zero projections
            !------------------------------------------------------------!
            if (.not. zero_projections_3j_condition(j_, j_prime_, lambda_)) cycle
            !------------------------------------------------------------!
            count_nonzero_algebraic_coefficients =                        &
               count_nonzero_algebraic_coefficients + 1
            !------------------------------------------------------------!
            ! calculate the Percival-Seaton coefficient
            !------------------------------------------------------------!
            pscoeff = percival_seaton_coefficient(j_, j_prime_, lambda_, omega_)
            !------------------------------------------------------------!
            ! save the Percival-Seaton coefficient
            !------------------------------------------------------------!
            nonzero_algebraic_coefficients(                               &
               count_nonzero_algebraic_coefficients) = pscoeff
            !------------------------------------------------------------!
            ! save indices to l1tab corresponding to \\(\lambda\\)
            !------------------------------------------------------------!
            nonzero_legendre_indices(count_nonzero_algebraic_coefficients)&
               = legendre_term_index_
            !------------------------------------------------------------!
            count_nonzero_legendre_terms = count_nonzero_legendre_terms + 1
         enddo
         !---------------------------------------------------------------------!
      end subroutine process_single_matrix_element
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine print_pes_matrix_elements_summary(number_of_channels,    &
         number_of_nonzero_pes_matrix_elements,                           &
         number_of_nonzero_algebraic_coefficients)
         !! print a shor summary on the number of non-zero matrix elements
         !! of the PES matrix
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: number_of_nonzero_pes_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element
            !! of the PES matrix
         integer(int32), intent(in) :: number_of_nonzero_algebraic_coefficients
            !! number of all non-zero algberaix coefficients in the whole
            !! PES matrix
         !---------------------------------------------------------------------!
         call write_message(" - Size of the PES matrix: "//                    &
            integer_to_character(number_of_channels))
         call write_message(" - Number of non-zero elements " //               &
            "of the potential matrix: " // integer_to_character(               &
            number_of_nonzero_pes_matrix_elements))
         call write_message(" - Number of non-zero elements " //               &
            " of the PES matrix: " // integer_to_character(                    &
            number_of_nonzero_algebraic_coefficients))
         !---------------------------------------------------------------------!
      end subroutine print_pes_matrix_elements_summary
   !---------------------------------------------------------------------------!
   ! Subroutines calculating the full PES matrix at desired R;
   ! these are called within the propagator loop
   !---------------------------------------------------------------------------!
      subroutine calculate_pes_matrix(total_angular_momentum_,                 &
         intermolecular_distance_, channel_indices_, channels_omega_values_,   &
         nonzero_terms_per_element_, nonzero_legendre_indices_,                &
         nonzero_algebraic_coefficients_, vmatrix)
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
         integer(int32), intent(in) :: nonzero_terms_per_element_(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         real(dp), intent(out) :: vmatrix(:,:)
            !! (output) - the interaction potential contribution to the coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_algebraic_coefficients_,              &
            count_nonzero_coupling_matrix_elements,                            &
            count_nonzero_legendre_terms, channel_index_1_, channel_index_2_,  &
            omega_, omega_prime_
         !---------------------------------------------------------------------!
         vmatrix   = 0
         count_nonzero_algebraic_coefficients_    = 0
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
                  = nonzero_terms_per_element_(count_nonzero_coupling_matrix_elements)
               !---------------------------------------------------------------!
               ! implementation of Eq. 1 in "Coupling Matrix" section
               !---------------------------------------------------------------!
               vmatrix(channel_index_1_, channel_index_2_)                     &
                  = calculate_single_pes_matrix_element(                       &
                     intermolecular_distance_, channel_index_1_,               &
                     channel_index_2_, channel_indices_,                       &
                     count_nonzero_legendre_terms,                             &
                     count_nonzero_algebraic_coefficients_,                      &
                     nonzero_legendre_indices_, nonzero_algebraic_coefficients_)
               !---------------------------------------------------------------!
             enddo
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(vmatrix,'u')
         !---------------------------------------------------------------------!
      end subroutine calculate_pes_matrix
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function calculate_single_pes_matrix_element(intermolecular_distance_,   &
         channel_index_1_, channel_index_2_, channel_indices_,                 &
         count_nonzero_legendre_terms, count_nonzero_algebraic_coefficients_,    &
         nonzero_legendre_indices_, nonzero_algebraic_coefficients_) result(matrix_element_)
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
         integer(int32), intent(in) :: nonzero_legendre_indices_(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(in) :: nonzero_algebraic_coefficients_(:)
            !! holds the values of the non-zero algebraic coefficients
         integer(int32), intent(inout) :: count_nonzero_algebraic_coefficients_
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
            count_nonzero_algebraic_coefficients_                                &
               = count_nonzero_algebraic_coefficients_ + 1
            !------------------------------------------------------------------!
            lambda_ = l1tab(nonzero_legendre_indices_(count_nonzero_algebraic_coefficients_))
            algebraic_coefficient_ = nonzero_algebraic_coefficients_(count_nonzero_algebraic_coefficients_)
            !------------------------------------------------------------------!
            call get_radial_coupling_term_value(intermolecular_distance_,      &
               lambda_, v_, j_, v_prime_, j_prime_, radial_term_)
            !------------------------------------------------------------------!
            sum_over_lambda_ = sum_over_lambda_ + algebraic_coefficient_ * radial_term_
            !------------------------------------------------------------------!
         enddo
         matrix_element_ =  2.0_dp * reduced_mass * sum_over_lambda_
         !---------------------------------------------------------------------!
         if (channel_index_1_ == channel_index_2_) then
            matrix_element_ = matrix_element_                                  &
               - wavenumber_squared_from_energy(internal_energy_)
         endif
         !---------------------------------------------------------------------!
      end function calculate_single_pes_matrix_element
!------------------------------------------------------------------------------!
end module pes_matrix_mod
