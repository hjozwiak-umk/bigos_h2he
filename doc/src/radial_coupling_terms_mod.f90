module radial_coupling_terms_mod
   !! This module provides all functions that handle radial coupling terms
   !! of the PES. It covers:
   !---------------------------------------------------------------------------!
   !! (1) reading radial coupling terms from external file
   !!    ("read_radial_coupling_terms", "skip_header_lines",
   !!    "read_and_validate_lambda", "read_potential_data", "validate_r_range")
   !---------------------------------------------------------------------------!
   !! (2) reducing the number of read coupling terms to retain only necessary
   !!    couplings ("reduce_radial_coupling_terms", "print_pes_quantum_numbers",
   !!    "reduce_coupling_terms", "find_reduced_term")
   !---------------------------------------------------------------------------!
   !! (3) interpolation of radial coupling terms
   !!    ("interpolate_radial_coupling_terms")
   !---------------------------------------------------------------------------!
   !! (4) providing value of the interpolated radial coupling term
   !!    ("get_radial_coupling_term_value")
   !!--------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: file_io_status, write_error, write_message,&
      integer_to_character, float_to_character
   use math_utilities_mod, only: spline, ispline
   use global_variables_mod
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: read_radial_coupling_terms, reduce_radial_coupling_terms,         &
      interpolate_radial_coupling_terms, get_radial_coupling_term_value
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                            Reading procedures
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine read_radial_coupling_terms
         !! Reads the radial coupling terms from the external file.
         !! The file is assumed to be formatted as described in
         !! "Supplying radial terms" section.
         !! The read radial coupling terms are kept in
         !! "tabulated_coupling_terms"
         !---------------------------------------------------------------------!
         character(len = 200) :: err_message
         integer(int32) :: lambda_index_, io_status
         !---------------------------------------------------------------------!
         open (coupling_terms_file_unit, file=trim(coupling_terms_file_name),  &
            form='formatted', status='old', iostat = io_status,                &
            iomsg = err_message)
         call file_io_status(io_status, err_message, coupling_terms_file_unit, &
            'o')
         !---------------------------------------------------------------------!
         ! Skip the informative lines at the beginning                         
         !---------------------------------------------------------------------!
         call skip_header_lines
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, number_of_legendre_indices
            call read_and_validate_lambda(lambda_index_)
            call read_potential_data(lambda_index_)
         enddo
         !---------------------------------------------------------------------!
         close(coupling_terms_file_unit, iostat = io_status,                   &
            iomsg = err_message)
         call file_io_status(io_status, err_message, coupling_terms_file_unit, &
            'c')
         !---------------------------------------------------------------------!
         ! Check if supplied radial terms cover a sufficient range of R
         !---------------------------------------------------------------------!
         call validate_r_range
         !---------------------------------------------------------------------!
      end subroutine read_radial_coupling_terms
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine skip_header_lines
         !! Skips the first n_skip_lines (read on input) lines
         !! in the coupling_terms_file
         !---------------------------------------------------------------------!
         integer(int32) :: line_index_
         !---------------------------------------------------------------------!
         do line_index_ = 1, n_skip_lines
            read(coupling_terms_file_unit, *) 
         enddo
         !---------------------------------------------------------------------!
      end subroutine skip_header_lines
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine read_and_validate_lambda(lambda_index_)
         !! Reads the value of lambda and compares with expected value.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_index_
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_
         !---------------------------------------------------------------------!
         read (coupling_terms_file_unit, *) lambda_
         if (lambda_.ne.legendre_indices(lambda_index_)) then
            close(coupling_terms_file_unit)
            close(s_matrix_unit)
            if (print_partial_cross_sections) close(partial_file_unit)
            call write_error("read_radial_coupling_terms: lambda = " //        &
               trim(adjustl(integer_to_character(lambda_))) //                 &
               " differs from expected value in legendre_indices (" //         &
               trim(adjustl(integer_to_character(lambda_index_))) // ") = " // &
               trim(adjustl(integer_to_character(legendre_indices(lambda_index_)))))
         endif
         !---------------------------------------------------------------------!
      end subroutine read_and_validate_lambda
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine read_potential_data(lambda_index_)
         !! Reads the intermolecular distance and radial coupling terms formatted
         !! in columns by iterating over number of tabulated ]](R\\) points.
         !! Immediately converts \\(R\\) and radial coupling terms to a.u.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_index_
         !---------------------------------------------------------------------!
         integer(int32) :: r_index_, column_index_
         real(dp) :: read_line(total_number_of_coupling_terms + 1)
         !---------------------------------------------------------------------!
         do r_index_ = 1, number_of_r_points
            read(coupling_terms_file_unit, *) (read_line(column_index_),       &
               column_index_ = 1, total_number_of_coupling_terms + 1)
            r_grid(r_index_) = read_line(1) * radial_term_distance_converter
            do column_index_ = 1, total_number_of_coupling_terms
               tabulated_coupling_terms(r_index_, lambda_index_, column_index_)&
                  = read_line(column_index_ + 1) * radial_term_energy_converter
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine read_potential_data
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine validate_r_range
         !! Checks if read R values are consistent with r_min and r_max.
         !---------------------------------------------------------------------!
         if (r_min < r_grid(1)) then
            close(s_matrix_unit)
            if (print_partial_cross_sections) close(partial_file_unit)
            call write_error("r_min value provided by the user (" //           &
               trim(adjustl(float_to_character(r_min, "(F10.4)"))) //          &
               ") is smaller than r_min supplied in " //                       &
               trim(adjustl(coupling_terms_file_name)) // " ( "//              &
               trim(adjustl(float_to_character(r_grid(1), "(F10.4)"))) // ")")
         endif
         !---------------------------------------------------------------------!
         if (r_max > r_grid(number_of_r_points)) then
            close(s_matrix_unit)
            if (print_partial_cross_sections) close(partial_file_unit)
            call write_error("r_max value provided by the user (" //           &
               trim(adjustl(float_to_character(r_max, "(F10.4)"))) //          &
               ") is larger than r_max supplied in " //                        &
               trim(adjustl(coupling_terms_file_name)) // " ( "//  &
               trim(adjustl(float_to_character(r_grid(number_of_r_points),     &
               "(F10.4)"))) // ")")
         endif
         !---------------------------------------------------------------------!
      end subroutine validate_r_range
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                            Reducing procedures
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine reduce_radial_coupling_terms
         !! Reduces the tabulated_coupling_terms matrix to retain
         !! only the necessary coupling terms.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, radial_index_, coupling_index_
         !---------------------------------------------------------------------!
         if (total_number_of_coupling_terms                                    &
            /= minimal_number_of_coupling_terms) then
            call write_message("-- Reducing the number of the radial " //      &
               "coupling terms...")
            call print_pes_quantum_numbers("Original",                         &
               total_number_of_coupling_terms)
            call reduce_coupling_terms()
            call print_pes_quantum_numbers("Reduced",                          &
               minimal_number_of_coupling_terms)
            call write_message("-- Reduced "//                                 &
               trim(adjustl(integer_to_character(                              &
               total_number_of_coupling_terms))) // " radial terms to "//      &
               trim(adjustl(integer_to_character(                              &
               minimal_number_of_coupling_terms))))
         else
            !------------------------------------------------------------------!
            ! if there is nothing to be reduced,
            ! copy tabulated_coupling_terms to coupling_terms
            !------------------------------------------------------------------!
            coupling_terms = tabulated_coupling_terms
         endif
         !---------------------------------------------------------------------!
         deallocate(tabulated_coupling_terms)
         !---------------------------------------------------------------------!
      end subroutine reduce_radial_coupling_terms
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine reduce_coupling_terms
         !! Reduces the coupling terms based on the existence of couplings.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, r_index_, coupling_index_
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, number_of_legendre_indices
            do r_index_ = 1, number_of_r_points
               do coupling_index_ = 1, minimal_number_of_coupling_terms
                  coupling_terms(r_index_, lambda_index_, coupling_index_) =   &
                     find_reduced_term(r_index_, lambda_index_, coupling_index_)
               enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine reduce_coupling_terms
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function find_reduced_term(r_index_, lambda_index_, coupling_index_)     &
         result(reduced_term)
         !! Finds and returns the reduced term for the given indices.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: r_index_
            !! index interating over intermolecular grid         
         integer(int32), intent(in) :: lambda_index_
            !! index interating over Legendre expansion terms
         integer(int32), intent(in) :: coupling_index_
            !! index interating over necessary couplings
         real(dp) :: reduced_term
            !! (output) sought value of the coupling term
         !---------------------------------------------------------------------!
         integer(int32) :: column_index_
         !---------------------------------------------------------------------!
         reduced_term = 0.0_dp
         do column_index_ = 1, total_number_of_coupling_terms
         !---------------------------------------------------------------------!
         ! iterate over quantum numbers describing all couplings
         ! (v1/j1/v1p/rot_prime_couplings) until necessary couplings are found
         !---------------------------------------------------------------------!
            if ((reduced_rot_couplings(coupling_index_)                        &
                  == rot_couplings(column_index_)).and.                        &
               (reduced_rot_prime_couplings(coupling_index_)                   &
                  == rot_prime_couplings(column_index_)).and.                  &
               (reduced_vib_couplings(coupling_index_)                         &
                  == vib_couplings(column_index_)).and.                        &
               (reduced_vib_prime_couplings(coupling_index_)                   &
                  == vib_prime_couplings(column_index_))) then
               reduced_term = tabulated_coupling_terms(r_index_, lambda_index_,&
                  column_index_)
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_reduced_term
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_pes_quantum_numbers(set_type, col_count)
         !! Prints quantum numbers describing radial coupling terms of the PES
         !! based on the provided set type and column count.
         !---------------------------------------------------------------------!
         character(len=*), intent(in) :: set_type
            !! "Original" or "reduced" - describes the set of quantum numbers
         integer(int32), intent(in) :: col_count
            !! number of coupling terms
         !---------------------------------------------------------------------!
         integer(int32) :: column_index_
         !---------------------------------------------------------------------!
         if (print_level >= 3) then
            call write_message("*** " // trim(set_type) //                     &
               " number of quantum numbers describing radial coupling terms: " &
               // trim(adjustl(integer_to_character(col_count))) // " ***")
            call write_message("Set of quantum numbers:")
            call write_message("       v1  j1  v1`  j1`")
            select case(set_type)
               case("Original")
                  do column_index_ = 1, col_count
                     write(*,"(5X,2(2X,I2),2(2X,I2))")                         &
                        vib_couplings(column_index_),                          &
                        rot_couplings(column_index_),                          &
                        vib_prime_couplings(column_index_),                    &
                        rot_prime_couplings(column_index_)
                  enddo
               case("Reduced")
                  do column_index_ = 1, col_count
                     write(*,"(5X,2(2X,I2),2(2X,I2))")                         &
                        reduced_vib_couplings(column_index_),                  &
                        reduced_rot_couplings(column_index_),                  &
                        reduced_vib_prime_couplings(column_index_),            &
                        reduced_rot_prime_couplings(column_index_)
                  enddo
            end select
         endif
      end subroutine print_pes_quantum_numbers
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      !                        Interpolation procedure
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      subroutine interpolate_radial_coupling_terms
         !! Interpolates the radial coupling terms using cubic spline functions.
         !! The resulting spline coefficients for each coupling term
         !! are stored in coupling_terms_b_coeffs, coupling_terms_c_coeffs,
         !! and coupling_terms_d_coeffs matrices.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, coupling_index_
         real(dp) :: spline_coeff_b(number_of_r_points),                       &
            spline_coeff_c(number_of_r_points), spline_coeff_d(number_of_r_points)
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, number_of_legendre_indices
            do coupling_index_ = 1, minimal_number_of_coupling_terms
               !---------------------------------------------------------------!
               ! Compute spline coefficients for each coupling term
               !---------------------------------------------------------------!
               call SPLINE(number_of_r_points,r_grid,coupling_terms(           &
                  :,lambda_index_,coupling_index_), spline_coeff_b,            &
                  spline_coeff_c, spline_coeff_d)
               !---------------------------------------------------------------!
               ! Store coefficients in the respective matrices
               !---------------------------------------------------------------!
               coupling_terms_b_coeffs(:,lambda_index_,coupling_index_)        &
                  = spline_coeff_b
               coupling_terms_c_coeffs(:,lambda_index_,coupling_index_)        &
                  = spline_coeff_c
               coupling_terms_d_coeffs(:,lambda_index_,coupling_index_)        &
                  = spline_coeff_d
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine interpolate_radial_coupling_terms
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      !                     Radial coupling term value
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      subroutine get_radial_coupling_term_value(intermolecular_distance,       &
         lambda_, v_, j_, v_prime_, j_prime_, radial_term_value_)
         !! Returns the interpolated value of a specific radial coupling term
         !! at a given distance.
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: intermolecular_distance
            !! Intermolecular distance, \\(R\\)
         integer(int32), intent(in) :: lambda_
            !! Legendre expansion index
         integer(int32), intent(in) :: v_
            !! pre-collisional vibrational quantum number
         integer(int32), intent(in) :: j_
            !! pre-collisional rotational quantum number
         integer(int32), intent(in) :: v_prime_
            !! post-collisional vibrational quantum number
         integer(int32), intent(in) :: j_prime_
            !! post-collisional rotational quantum number
         real(dp), intent(out) :: radial_term_value_
            !! Value of the radial coupling coefficient
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index, coupling_index
         !---------------------------------------------------------------------!
         lambda_index = find_lambda_index(lambda_)
         !---------------------------------------------------------------------!
         if (lambda_index == 0) then
            call handle_lambda_index_error(lambda_)
            return
         endif
         !---------------------------------------------------------------------!
         coupling_index = find_coupling_index(v_, j_, v_prime_, j_prime_)
         !---------------------------------------------------------------------!
         if (coupling_index == 0) then
            call handle_coupling_index_error(v_, j_, v_prime_, j_prime_)
            return
         endif
         !---------------------------------------------------------------------!
         radial_term_value_ = ISPLINE(intermolecular_distance,                 &
            number_of_r_points, r_grid, coupling_terms(:, lambda_index,        &
            coupling_index), coupling_terms_b_coeffs(:, lambda_index,          &
            coupling_index), coupling_terms_c_coeffs(:, lambda_index,          &
            coupling_index), coupling_terms_d_coeffs(:, lambda_index,          &
            coupling_index))
         !---------------------------------------------------------------------!
      end subroutine get_radial_coupling_term_value
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      function find_lambda_index(lambda_) result(result_index_)
         !! Locates given \\(\lambda\\) value in legendre_indices.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_
            !! Legendre expansion index, \\(\lambda\\)
         integer(int32) :: result_index_
            !! Index pointing to  \\(\lambda\\) in legendre_indices
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_
         !---------------------------------------------------------------------!
         result_index_ = 0
         do lambda_index_ = 1, number_of_legendre_indices
            if (legendre_indices(lambda_index_) == lambda_) then
               result_index_ = lambda_index_
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_lambda_index
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      subroutine handle_lambda_index_error(lambda_)
         !! Handles error when \\(\lambda\\) is not found in legendre_indices.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_
            !! Legendre expansion index, \\(\lambda\\)
         !---------------------------------------------------------------------!
         call write_error("Radial coupling terms with lambda = " //            &
            trim(adjustl(integer_to_character(lambda_))) //                    &
            " not found in legendre_indices")
         !---------------------------------------------------------------------!
      end subroutine handle_lambda_index_error
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      function find_coupling_index(v_, j_, v_prime_, j_prime_) result(result_index_)
         !! Locates the correct quantum number that describes the v/j coupling.
         !! Note that coupling terms are symmetric with respect to the change
         !! of pre- and post-collisional quantum numbers.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: v_
            !! pre-collisional vibrational quantum number
         integer(int32), intent(in) :: j_
            !! pre-collisional rotational quantum number
         integer(int32), intent(in) :: v_prime_
            !! post-collisional vibrational quantum number
         integer(int32), intent(in) :: j_prime_
            !! post-collisional rotational quantum number
         integer(int32) :: result_index_
            !! Index pointing to  \\(v, j, v^{\prime}, j^{\prime}\\)
            !! in reduced_* arrays
         !---------------------------------------------------------------------!
         integer(int32) :: coupling_index_
         !---------------------------------------------------------------------!
         result_index_ = 0
         do coupling_index_ = 1, minimal_number_of_coupling_terms
            if ((((reduced_vib_couplings(coupling_index_) == v_).and.          &
               (reduced_rot_couplings(coupling_index_).eq.j_).and.             &
               (reduced_vib_prime_couplings(coupling_index_).eq.v_prime_).and. &
               (reduced_rot_prime_couplings(coupling_index_).eq.j_prime_))     &
               .or.                                                            &
               ((reduced_vib_couplings(coupling_index_).eq.v_prime_).and.      &
               (reduced_rot_couplings(coupling_index_).eq.j_prime_).and.       &
               (reduced_vib_prime_couplings(coupling_index_).eq.v_).and.       &
               (reduced_rot_prime_couplings(coupling_index_).eq.j_))) ) then
               result_index_ = coupling_index_
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_coupling_index
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!      
      subroutine handle_coupling_index_error(v_, j_, v_prime_, j_prime_)
         !! Handles error when the appropriate coupling term is not found.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: v_
            !! pre-collisional vibrational quantum number
         integer(int32), intent(in) :: j_
            !! pre-collisional rotational quantum number
         integer(int32), intent(in) :: v_prime_
            !! post-collisional vibrational quantum number
         integer(int32), intent(in) :: j_prime_
            !! post-collisional rotational quantum number
         !---------------------------------------------------------------------!
         call write_error("Coupling term with v = " //                         &
            trim(adjustl(integer_to_character(v_))) // ", j = " //             &
            trim(adjustl(integer_to_character(j_))) // ", v` = " //            &
            trim(adjustl(integer_to_character(v_prime_))) // ", j1` = " //     &
            trim(adjustl(integer_to_character(j_prime_))) // " not found")
         !---------------------------------------------------------------------!
      end subroutine handle_coupling_index_error
   !---------------------------------------------------------------------------!
end module radial_coupling_terms_mod
