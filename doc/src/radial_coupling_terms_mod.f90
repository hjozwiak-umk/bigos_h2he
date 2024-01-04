module radial_coupling_terms_mod
   !! This module provides all functions that handle radial coupling terms
   !! of the PES. It covers:
   !! 1. reading radial coupling terms from external file ("read_radial_coupling_terms",
   !!    "skip_header_lines", "read_and_validate_lambda", "read_potential_data",
   !!    "validate_r_range")
   !! 2. reducing the number of read coupling terms to retain only necessary
   !!    couplings ("reduce_radial_coupling_terms", "print_pes_quantum_numbers",
   !!    "reduce_coupling_terms", "find_reduced_term")
   !! 3. interpolation of radial coupling terms ("interpolate_radial_coupling_terms")
   !! 4. providing value of the interpolated radial coupling term ("get_radial_coupling_term_value")
   !!--------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use utility_functions_mod, only: file_io_status, write_error, write_message,&
      time_count_summary
   use math_functions_mod, only: spline, ispline
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: read_radial_coupling_terms, reduce_radial_coupling_terms,         &
      interpolate_radial_coupling_terms, get_radial_coupling_term_value
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
   !                            Reading procedures
   !---------------------------------------------------------------------------!
      subroutine read_radial_coupling_terms
         !! Reads the radial coupling terms from the external file.
         !! The file is assumed to be formatted as described in
         !! "Supplying radial terms" section.
         !! The read radial coupling terms are kept in vmat/read_vmat3D
         !---------------------------------------------------------------------!
         character(len = 200) :: err_message
         integer(int32) :: nrtmp, l1, iskip_, lambda_index_, ir, icol, io_status
         !---------------------------------------------------------------------!
         open (pes_file_unit, file=trim(potentialfile), form='formatted',      &
            status='old', iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, pes_file_unit, 'o')
         !---------------------------------------------------------------------!
         ! Skip the informative lines at the beginning                         
         !---------------------------------------------------------------------!
         call skip_header_lines
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, nterms
            call read_and_validate_lambda(lambda_index_)
            call read_potential_data(lambda_index_)
         enddo
         !---------------------------------------------------------------------!
         close(pes_file_unit, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, pes_file_unit, 'c')
         !---------------------------------------------------------------------!
         ! Check if supplied radial terms cover a sufficient range of R
         !---------------------------------------------------------------------!
         call validate_r_range
         !---------------------------------------------------------------------!
      end subroutine read_radial_coupling_terms
   !---------------------------------------------------------------------------!
      subroutine skip_header_lines
         !! Skips the first n_skip_lines (read on input) lines in the pes_file
         !---------------------------------------------------------------------!
         integer(int32) :: line_index_
         !---------------------------------------------------------------------!
         do line_index_ = 1, n_skip_lines
            read(pes_file_unit, *) 
         enddo
         !---------------------------------------------------------------------!
      end subroutine skip_header_lines
   !---------------------------------------------------------------------------!
      subroutine read_and_validate_lambda(lambda_index_)
         !! Reads the value of lambda and compares with expected value.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_index_
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_
         !---------------------------------------------------------------------!
         read (pes_file_unit, *) lambda_
         if (lambda_.ne.l1tab(lambda_index_)) then
            close(pes_file_unit)
            close(s_matrix_unit)
            if (ipart.eq.1) close(partial_file_unit)
            call write_error("read_radial_coupling_terms: lambda = " //                    &
               trim(adjustl(integer_to_character(lambda_))) //                 &
               " differs from expected value in l1tab (" //                    &
               trim(adjustl(integer_to_character(lambda_index_))) // ") = " // &
               trim(adjustl(integer_to_character(l1tab(lambda_index_)))))
         endif
         !---------------------------------------------------------------------!
      end subroutine read_and_validate_lambda
   !---------------------------------------------------------------------------!
      subroutine read_potential_data(lambda_index_)
         !! Reads the intermolecular distance and radial coupling terms formatted
         !! in columns by iterating over number of tabulated ]](R\\) points.
         !! Immediately converts \\(R\\) and radial coupling terms to a.u.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_index_
         !---------------------------------------------------------------------!
         integer(int32) :: r_index_, column_index_
         real(dp) :: read_line(totalcol + 1)
         !---------------------------------------------------------------------!
         do r_index_ = 1, nr
            read(pes_file_unit, *) (read_line(column_index_),                  &
               column_index_ = 1, totalcol + 1)
            rmat(r_index_) = read_line(1) * radial_term_distance_converter
            do column_index_ = 1, totalcol
               read_vmat3D(r_index_, lambda_index_, column_index_)             &
                  = read_line(column_index_ + 1) * radial_term_energy_converter
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine read_potential_data
   !---------------------------------------------------------------------------!
      subroutine validate_r_range
         !! Checks if read R values are consistent with rmin and rmax.
         !---------------------------------------------------------------------!
         if (rmin < rmat(1)) then
            close(s_matrix_unit)
            if (ipart.eq.1) close(partial_file_unit)
            call write_error("rmin value provided by the user (" //            &
               trim(adjustl(float_to_character(rmin, "(F10.4)"))) //           &
               ") is smaller than rmin supplied in " //                        &
               trim(adjustl(potentialfile)) // " ( "//                         &
               trim(adjustl(float_to_character(rmat(1), "(F10.4)"))) // ")")
         endif
         !---------------------------------------------------------------------!
         if (rmax > rmat(nr)) then
            close(s_matrix_unit)
            if (ipart.eq.1) close(partial_file_unit)
            call write_error("rmax value provided by the user (" //            &
               trim(adjustl(float_to_character(rmax, "(F10.4)"))) //           &
               ") is larger than rmax supplied in " //                         &
               trim(adjustl(potentialfile)) // " ( "//  &
               trim(adjustl(float_to_character(rmat(nr), "(F10.4)"))) // ")")
         endif
         !---------------------------------------------------------------------!
      end subroutine validate_r_range
   !---------------------------------------------------------------------------!
   !                            Reducing procedures
   !---------------------------------------------------------------------------!
      subroutine reduce_radial_coupling_terms
         !! Reduces the read_vmat3D matrix to retain only the necessary coupling terms.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, radial_index_, coupling_index_
         !---------------------------------------------------------------------!
         if (totalcol /= ncoupl) then
            call write_message("Reducing the number of the radial coupling terms...")
            call print_pes_quantum_numbers("Original", totalcol)
            call reduce_coupling_terms()
            call print_pes_quantum_numbers("Reduced", ncoupl)
            call write_message("Reduced "//                                    &
               trim(adjustl(integer_to_character(totalcol))) //                &
               " radial terms to "// trim(adjustl(integer_to_character(ncoupl))))
         else
            !------------------------------------------------------------------!
            ! if there is nothing to be reduced, copy read_vmat3d to vmat3d
            !------------------------------------------------------------------!
            vmat3D = read_vmat3D
         endif
         !---------------------------------------------------------------------!
         deallocate(read_vmat3D)
         !---------------------------------------------------------------------!
      end subroutine reduce_radial_coupling_terms
   !---------------------------------------------------------------------------!
      subroutine reduce_coupling_terms
         !! Reduces the coupling terms based on the existence of couplings.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, r_index_, coupling_index_
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, nterms
            do r_index_ = 1, nr
               do coupling_index_ = 1, ncoupl
                  vmat3D(r_index_, lambda_index_, coupling_index_) =           &
                     find_reduced_term(r_index_, lambda_index_, coupling_index_)
               enddo
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine reduce_coupling_terms
   !---------------------------------------------------------------------------!
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
         do column_index_ = 1, totalcol
         !---------------------------------------------------------------------!
         ! iterate over quantum numbers describing all couplings
         ! (v1/j1/v1p/j1ppes) until necessary couplings are found
         !---------------------------------------------------------------------!
            if ((reduced_j1pes(coupling_index_) == j1pes(column_index_)).and.  &
               (reduced_j1ppes(coupling_index_) == j1ppes(column_index_)).and. &
               (reduced_v1pes(coupling_index_) == v1pes(column_index_)).and.   &
               (reduced_v1ppes(coupling_index_) == v1ppes(column_index_))) then
               reduced_term = read_vmat3D(r_index_, lambda_index_, column_index_)
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_reduced_term
   !---------------------------------------------------------------------------!
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
         if (prntlvl >= 3) then
            call write_message("*** " // trim(set_type) //                     &
               " number of quantum numbers describing radial coupling terms: " &
               // trim(adjustl(integer_to_character(col_count))) // " ***")
            call write_message("Set of quantum numbers:")
            call write_message("       v1  j1  v1`  j1`")
            select case(set_type)
               case("Original")
                  do column_index_ = 1, col_count
                     write(*,"(5X,2(2X,I2),2(2X,I2))") v1pes(column_index_),   &
                        j1pes(column_index_), v1ppes(column_index_),           &
                        j1ppes(column_index_)
                  enddo
               case("Reduced")
                  do column_index_ = 1, col_count
                     write(*,"(5X,2(2X,I2),2(2X,I2))")                         &
                        reduced_v1pes(column_index_),                          &
                        reduced_j1pes(column_index_),                          &
                        reduced_v1ppes(column_index_),                         &
                        reduced_j1ppes(column_index_)
                  enddo
            end select
         endif
      end subroutine print_pes_quantum_numbers
   !---------------------------------------------------------------------------!
   !                        Interpolation procedure
   !---------------------------------------------------------------------------!
      subroutine interpolate_radial_coupling_terms
         !! Interpolates the radial coupling terms using cubic spline functions.
         !! The resulting spline coefficients for each coupling term
         !! are stored in bmat3D, cmat3D, and dmat3D matrices.
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_, coupling_index_
         !---------------------------------------------------------------------!
         do lambda_index_ = 1, nterms
            do coupling_index_ = 1, ncoupl
               !---------------------------------------------------------------!
               ! Compute spline coefficients for each coupling term
               !---------------------------------------------------------------!
               call SPLINE(nr,rmat,vmat3D(:,lambda_index_,coupling_index_),    &
                  spline_coeff_b, spline_coeff_c, spline_coeff_d)
               !---------------------------------------------------------------!
               ! Store coefficients in the respective matrices
               !---------------------------------------------------------------!
               bmat3D(:,lambda_index_,coupling_index_) = spline_coeff_b
               cmat3D(:,lambda_index_,coupling_index_) = spline_coeff_c
               dmat3D(:,lambda_index_,coupling_index_) = spline_coeff_d
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine interpolate_radial_coupling_terms
   !---------------------------------------------------------------------------!
   !                     Radial coupling term value
   !---------------------------------------------------------------------------!
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
         radial_term_value_ = ISPLINE(intermolecular_distance, nr, rmat,       &
            vmat3D(:, lambda_index, coupling_index),                           &
            bmat3D(:, lambda_index, coupling_index),                           &
            cmat3D(:, lambda_index, coupling_index),                           &
            dmat3D(:, lambda_index, coupling_index))
         !---------------------------------------------------------------------!
      end subroutine get_radial_coupling_term_value
   !---------------------------------------------------------------------------!
      function find_lambda_index(lambda_) result(result_index_)
         !! Locates given \\(\lambda\\) value in l1tab.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_
            !! Legendre expansion index, \\(\lambda\\)
         integer(int32) :: result_index_
            !! Index pointing to  \\(\lambda\\) in l1tab
         !---------------------------------------------------------------------!
         integer(int32) :: lambda_index_
         !---------------------------------------------------------------------!
         result_index_ = 0
         do lambda_index_ = 1, nterms
            if (l1tab(lambda_index_) == lambda_) then
               result_index_ = lambda_index_
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_lambda_index
   !---------------------------------------------------------------------------!
      subroutine handle_lambda_index_error(lambda_)
         !! Handles error when \\(\lambda\\) is not found in l1tab.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: lambda_
            !! Legendre expansion index, \\(\lambda\\)
         !---------------------------------------------------------------------!
         call write_error("Radial coupling terms with lambda = " //            &
            trim(adjustl(integer_to_character(lambda_))) // " not found in l1tab")
         !---------------------------------------------------------------------!
      end subroutine handle_lambda_index_error
   !---------------------------------------------------------------------------!
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
            !! Index pointing to  \\(v, j, v^{\prime}, j^{\prime}\\) in reduced_* arrays
         !---------------------------------------------------------------------!
         integer(int32) :: coupling_index_
         !---------------------------------------------------------------------!
         result_index_ = 0
         do coupling_index_ = 1, ncoupl
            if ((((reduced_v1pes(coupling_index_) == v_).and.                  &
               (reduced_j1pes(coupling_index_).eq.j_).and.                     &
               (reduced_v1ppes(coupling_index_).eq.v_prime_).and.              &
               (reduced_j1ppes(coupling_index_).eq.j_prime_))                  &
               .or.                                                            &
               ((reduced_v1pes(coupling_index_).eq.v_prime_).and.              &
               (reduced_j1pes(coupling_index_).eq.j_prime_).and.               &
               (reduced_v1ppes(coupling_index_).eq.v_).and.                    &
               (reduced_j1ppes(coupling_index_).eq.j_))) ) then
               result_index_ = coupling_index_
               exit
            endif
         enddo
         !---------------------------------------------------------------------!
      end function find_coupling_index
   !---------------------------------------------------------------------------!
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
!------------------------------------------------------------------------------!
end module radial_coupling_terms_mod
