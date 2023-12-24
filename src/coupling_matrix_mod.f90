module PES_COUPLING_MATRIX
   !! this module provides functions calculating the algebraic coefficients
   !! entering the coupling matrix (Eq. ...)
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use supplementary_mod, only: write_error, write_message, integer_to_character
   use fwigxjpf, only: fwig3jj
   use io_mod
   use supplementary
   implicit none
   contains
!------------------------------------------------------------------------------!
      subroutine check_nonzero_coupling_matrix_elements(number_of_channels,    &
         channels_level_indices, channels_omega_values,                        &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients)
         !! checks the number of non-zero coupling matrix elements due to
         !! the \bar{\Omega} = \bar{\Omega}' condition
         !! (number_of_nonzero_coupling_matrix_elements variable),
         !! and the total number of non-zero algebraic coefficients that enter
         !! Eqs. (6.13)-(6.15) (number_of_nonzero_coupling_coefficients variable).
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(out) :: number_of_nonzero_coupling_matrix_elements
            !! number of non-zero terms in the sum () for each non-zero element of the coupling matrix
         integer(int32), intent(out) :: number_of_nonzero_coupling_coefficients
            !! number of all non-zero algberaix coefficients in the whole coupling matrix
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_coupling_matrix_elements,             &
            count_nonzero_coupling_coefficients, nil, j1tmp, j1ptmp, omegatmp, &
            omegaptmp, l1, ii, ij, il
         !---------------------------------------------------------------------!
         count_nonzero_coupling_coefficients = 0
         count_nonzero_coupling_matrix_elements = 0
         do ii = 1, number_of_channels
            j1tmp = j1array(channels_level_indices(ii))
            omegatmp = channels_omega_values(ii)
            do ij = 1, ii
               j1ptmp = j1array(channels_level_indices(ij))
               omegaptmp = channels_omega_values(ij)
               if (omegatmp.ne.omegaptmp) cycle
               count_nonzero_coupling_matrix_elements =                        &
                  count_nonzero_coupling_matrix_elements + 1
               do il = 1,nterms
                  l1 = l1tab(il)
                  if (triangle_inequality_holds(j1tmp,j1ptmp,l1).eq.0) cycle
                  if (is_sum_even(j1tmp,j1ptmp,l1).eq.0) cycle
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
      subroutine prepare_coupling_matrix_elements(number_of_channels,          &
         channels_level_indices, channels_omega_values,                        &
         number_of_nonzero_coupling_matrix_elements,                           &
         number_of_nonzero_coupling_coefficients,                              &
         nonzero_terms_per_element, nonzero_legendre_indices,                  &
         nonzero_coupling_coefficients)
         !! prepares:
         !! -- nonzero_terms_per_element - keeps the number of
         !!    non-zero terms in the sum (Eq. (6.21)) for each non-zero matrix
         !!    element of the coupling matrix
         !! -- nonzero_legendre_indices - holds the proper
         !!    indices in the range (0, nterms) pointing to l1/l2/lltabs, which
         !!    correspond to the non-vanishing elements of the sum (Eq. (6.21))
         !!    for each non-zero  matrix element of the coupling matrix
         !! -- nonzero_coupling_coefficients --  holds values of non-zero Percival-Seaton coefficients
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
         integer(int32), intent(inout) ::                                      &
            nonzero_terms_per_element(number_of_nonzero_coupling_matrix_elements)
            !! keeps the number of non-zero terms in the sum (Eq. (6.21)) for
            !! each non-zero element of the coupling matrix
         integer(int32), intent(inout) ::                                      &
            nonzero_legendre_indices(number_of_nonzero_coupling_coefficients)
            !! holds proper indices pointing to l1tab, which
            !! correspond to the non-vanishing elements of the sum  (Eq. (6.21))
            !! for each non-zero element of the coupling matrix
         real(dp), intent(inout) ::                                            &
            nonzero_coupling_coefficients(number_of_nonzero_coupling_coefficients)
            !! holds the values of the non-zero algebraic coefficients
         !---------------------------------------------------------------------!
         integer(int32) :: count_nonzero_coupling_matrix_elements,             &
            count_nonzero_coupling_coefficients, nonzero_legendre, j1tmp,      &
            j1ptmp, omegatmp, omegaptmp, l1, ii, ij, il
         real(dp) :: pscoeff
         !---------------------------------------------------------------------!
         nonzero_terms_per_element        = 0
         nonzero_legendre_indices         = 0
         nonzero_coupling_coefficients    = 0
         count_nonzero_coupling_coefficients     = 0
         count_nonzero_coupling_matrix_elements  = 0
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels      
            j1tmp = j1array(channels_level_indices(ii))
            omegatmp = channels_omega_values(ii)
            do ij = 1, ii
               j1ptmp = j1array(channels_level_indices(ij))
               omegaptmp = channels_omega_values(ij)
               if (omegatmp.ne.omegaptmp) cycle
               !---------------------------------------------------------------!
               ! passed \bar{\Omega} = \bar{\Omega}' condition
               !---------------------------------------------------------------!
               count_nonzero_coupling_matrix_elements =                        &
                  count_nonzero_coupling_matrix_elements + 1
               !---------------------------------------------------------------!
               ! count non-zero terms in the sum over legendre polynomials
               ! for this element
               !---------------------------------------------------------------!
               nonzero_legendre = 0
               do il = 1, nterms
                  l1 = l1tab(il)
                  if (triangle_inequality_holds(j1tmp,j1ptmp,l1).eq.0) cycle
                  if (is_sum_even(j1tmp,j1ptmp,l1).eq.0) cycle
                  !------------------------------------------------------------!
                  ! passed non-zero conditions
                  !------------------------------------------------------------!
                  count_nonzero_coupling_coefficients =                        &
                     count_nonzero_coupling_coefficients + 1
                  !------------------------------------------------------------!
                  ! calculate the Percival-Seaton coefficient
                  !------------------------------------------------------------!
                  pscoeff = dsqrt(real((2 * j1tmp + 1) * (2 * j1ptmp + 1), dp))&
                     * fwig3jj(2* j1tmp ,   2* j1ptmp  , 2* l1, 0, 0, 0)       &
                     * fwig3jj(2* j1tmp ,   2* j1ptmp  , 2* l1,                &
                            2 * omegatmp, -2 * omegatmp,     0)                &
                     * (-1.0_dp)**(omegatmp)
                  !------------------------------------------------------------!
                  nonzero_coupling_coefficients(                               &
                     count_nonzero_coupling_coefficients) = pscoeff
                  nonzero_legendre_indices(count_nonzero_coupling_coefficients)&
                     = il
                  nonzero_legendre = nonzero_legendre + 1
               enddo
               nonzero_terms_per_element(count_nonzero_coupling_matrix_elements)&
                  = nonzero_legendre
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine prepare_coupling_matrix_elements
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
end module PES_COUPLING_MATRIX
