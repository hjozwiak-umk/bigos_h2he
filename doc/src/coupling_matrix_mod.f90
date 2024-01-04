module coupling_matrix_mod
   !! this module provides functions calculating the algebraic coefficients
   !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) entering the coupling matrix:
   !! see Eq. 1 in the "Coupling Matrix" section. The data is organized as follows:
   !! - number of non-zero terms of the coupling matrix due to
   !!   \\(\bar{\Omega} = \bar{\Omega}'\\) condition is saved as
   !!   "number_of_nonzero_coupling_matrix_elements"
   !! - number of non-vanishing terms in the sum over \\(\lambda\\)
   !!   in Eq. 1 in the "Coupling Matrix" section is saved in a
   !!   "nonzero_terms_per_element" array which is of
   !!   "number_of_nonzero_coupling_matrix_elements" size
   !! - a _total_ number of non-vanishing  \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\)
   !!   coefficients is saved in "number_of_nonzero_coupling_coefficients"
   !! - _all_ non-vanishing  \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\)
   !!   coefficients are saved in "nonzero_coupling_coefficients" array,
   !!   which is of "number_of_nonzero_coupling_coefficients" size
   !! - corresponding \\(\lambda\\) value for each non-vanishing coefficient
   !!   is saved as an index to "l1tab" in the "nonzero_legendre_indices"
   !!   array (which is of "number_of_nonzero_coupling_coefficients" size)
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_message, integer_to_character
   use io_mod
   use math_functions_mod, only: percival_seaton_coefficient,                  &
      triangle_inequality_holds, is_sum_even, zero_projections_3j_condition
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
!------------------------------------------------------------------------------!
      subroutine check_nonzero_coupling_matrix_elements(channels_level_indices,&
         channels_omega_values, number_of_nonzero_coupling_matrix_elements,    &
         number_of_nonzero_coupling_coefficients)
         !! checks the number of non-zero coupling matrix elements due to
         !! the \bar{\Omega} = \bar{\Omega}' condition,
         !! "number_of_nonzero_coupling_matrix_elements",
         !! and the total number of non-zero algebraic coefficients,
         !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\), in the whole matrix,
         !! "number_of_nonzero_coupling_coefficients".
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(out) :: number_of_nonzero_coupling_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element of the coupling matrix
         integer(int32), intent(out) :: number_of_nonzero_coupling_coefficients
            !! number of all non-zero algberaix coefficients in the whole coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_coupling_matrix_elements,             &
            count_nonzero_coupling_coefficients, j_, j_prime_, omega_,         &
            omega_prime_, lambda_, channel_index_1_, channel_index_2_, legendre_term_index_
         !---------------------------------------------------------------------!
         count_nonzero_coupling_coefficients = 0
         count_nonzero_coupling_matrix_elements = 0
         do channel_index_1_ = 1, size(channels_level_indices)
            j_ = j1array(channels_level_indices(channel_index_1_))
            omega_ = channels_omega_values(channel_index_1_)
            do channel_index_2_ = 1, channel_index_1_
               j_prime_ = j1array(channels_level_indices(channel_index_2_))
               omega_prime_ = channels_omega_values(channel_index_2_)
               !---------------------------------------------------------------!
               if (omega_ /= omega_prime_) cycle
               !---------------------------------------------------------------!
               count_nonzero_coupling_matrix_elements =                        &
                  count_nonzero_coupling_matrix_elements + 1
               do legendre_term_index_ = 1, nterms
                  lambda_ = l1tab(legendre_term_index_)
                  if (.not. zero_projections_3j_condition(j_, j_prime_, lambda_)) cycle
                  count_nonzero_coupling_coefficients =                        &
                     count_nonzero_coupling_coefficients+1            
              enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
         number_of_nonzero_coupling_coefficients = count_nonzero_coupling_coefficients
         number_of_nonzero_coupling_matrix_elements = count_nonzero_coupling_matrix_elements
         !---------------------------------------------------------------------!
      end subroutine check_nonzero_coupling_matrix_elements
!------------------------------------------------------------------------------!
      subroutine prepare_coupling_matrix_elements(channels_level_indices,      &
         channels_omega_values, nonzero_terms_per_element,                     &
         nonzero_legendre_indices, nonzero_coupling_coefficients)
         !! prepares:
         !! -- nonzero_terms_per_element - number of non-vanishing terms in
         !!    the sum over \\(\lambda\\) in Eq. 1 in the "Coupling Matrix" section
         !! -- nonzero_legendre_indices - corresponding \\(\lambda\\) value for
         !!    each non-vanishing coefficient is saved as an index to "l1tab"
         !! -- nonzero_coupling_coefficients --  holds _all_ non-vanishing
         !!    \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) coefficients
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(inout) :: nonzero_terms_per_element(:)
            !! keeps the number of non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix
         integer(int32), intent(inout) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(inout) :: nonzero_coupling_coefficients(:)
            !! holds the values of the non-zero algebraic coefficients
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_coupling_matrix_elements,             &
            count_nonzero_coupling_coefficients, count_nonzero_legendre_terms, &
            j_, j_prime_, omega_, omega_prime_, lambda_, channel_index_1_,     &
            channel_index_2_, legendre_term_index_
         real(dp) :: pscoeff
         !---------------------------------------------------------------------!
         nonzero_terms_per_element        = 0
         nonzero_legendre_indices         = 0
         nonzero_coupling_coefficients    = 0
         count_nonzero_coupling_coefficients     = 0
         count_nonzero_coupling_matrix_elements  = 0
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, size(channels_level_indices)
            j_     = j1array(channels_level_indices(channel_index_1_))
            omega_ = channels_omega_values(channel_index_1_)
            do channel_index_2_ = 1, channel_index_1_
               j_prime_     = j1array(channels_level_indices(channel_index_2_))
               omega_prime_ = channels_omega_values(channel_index_2_)
               if (omega_  /= omega_prime_) cycle
               !---------------------------------------------------------------!
               ! passed \bar{\Omega} = \bar{\Omega}' condition
               !---------------------------------------------------------------!
               count_nonzero_coupling_matrix_elements =                        &
                  count_nonzero_coupling_matrix_elements + 1
               !---------------------------------------------------------------!
               ! process a single matrix element:
               ! determine non-zero terms in the sum over legendre polynomials
               ! for this element; these are saved to ...
               !---------------------------------------------------------------!
               call process_single_matrix_element(j_, j_prime_, omega_,        &
                  count_nonzero_coupling_coefficients,                         &
                  count_nonzero_legendre_terms, nonzero_legendre_indices,      &
                  nonzero_coupling_coefficients)
               
               nonzero_terms_per_element(count_nonzero_coupling_matrix_elements)&
                  = count_nonzero_legendre_terms
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine prepare_coupling_matrix_elements
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
      subroutine process_single_matrix_element(j_, j_prime_, omega_,           &
         count_nonzero_coupling_coefficients, count_nonzero_legendre_terms,    &
         nonzero_legendre_indices, nonzero_coupling_coefficients)
         !! calculates the non-zero algebraic coefficients
         !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) for a single matrix
         !! element - see Eq. (1) in the "Coupling matrix" section;
         !! algebraic coefficients are saved to nonzero_coupling_coefficients
         !! array; corresponding indices to l1tab are saved to
         !! nonzero_legendre_indices array
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: j_
            !! pre-collisional rotational angular momentum
         integer(int32), intent(in) :: j_prime_
            !! post-collisional rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         integer(int32), intent(inout) :: count_nonzero_coupling_coefficients
            !! running index counting non-zero coupling coefficients,
            !! \\( g\_{{\lambda},\gamma,\gamma'}^{Jp} \\) in the whole matrix;
            !! incremented in this subroutine
         integer(int32), intent(inout) :: count_nonzero_legendre_terms
            !! number of non-zero terms in Legendre expansion for a single
            !! matrix element of the interaction potential
         integer(int32), intent(inout) :: nonzero_legendre_indices(:)
            !! holds indices pointing to l1tab, which correspond to
            !! the non-vanishing elements of the sum over \\(\lambda\\)
            !! for each non-zero element of the coupling matrix;
         real(dp), intent(inout) :: nonzero_coupling_coefficients(:)
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
            count_nonzero_coupling_coefficients =                        &
               count_nonzero_coupling_coefficients + 1
            !------------------------------------------------------------!
            ! calculate the Percival-Seaton coefficient
            !------------------------------------------------------------!
            pscoeff = percival_seaton_coefficient(j_, j_prime_, lambda_, omega_)
            !------------------------------------------------------------!
            ! save the Percival-Seaton coefficient
            !------------------------------------------------------------!
            nonzero_coupling_coefficients(                               &
               count_nonzero_coupling_coefficients) = pscoeff
            !------------------------------------------------------------!
            ! save indices to l1tab corresponding to \\(\lambda\\)
            !------------------------------------------------------------!
            nonzero_legendre_indices(count_nonzero_coupling_coefficients)&
               = legendre_term_index_
            !------------------------------------------------------------!
            count_nonzero_legendre_terms = count_nonzero_legendre_terms + 1
         enddo
         !---------------------------------------------------------------------!
      end subroutine process_single_matrix_element
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
      subroutine print_coupling_matrix_elements_summary(number_of_channels,    &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients)
         !! print a shor summary on the number of non-zero matrix elements
         !! of the coupling matrix
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: number_of_nonzero_coupling_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element
            !! of the coupling matrix
         integer(int32), intent(in) :: number_of_nonzero_coupling_coefficients
            !! number of all non-zero algberaix coefficients in the whole
            !! coupling matrix
         !---------------------------------------------------------------------!
         call write_message(repeat("*", 90))
         call write_message(repeat(" ", 5) // "Size of the coupling matrix: "//&
            integer_to_character(number_of_channels))
         call write_message(repeat(" ", 5) // "Number of non-zero elements " //&
            "of the potential matrix: " // integer_to_character(               &
            number_of_nonzero_coupling_matrix_elements))
         call write_message(repeat(" ", 5) // "Number of non-zero elements " //&
            " of the coupling matrix: " // integer_to_character(               &
            number_of_nonzero_coupling_coefficients))
         call write_message(repeat('*', 90))
         !---------------------------------------------------------------------!
      end subroutine print_coupling_matrix_elements_summary
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
end module coupling_matrix_mod
