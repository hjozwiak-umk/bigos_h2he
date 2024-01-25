module input_validation
   !! This module provides subroutines validating read variables' values.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_message,&
      incorrect_value, integer_to_character, to_lowercase
   use global_variables_mod
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine check_namelist_input
         !! Check variables read from namelist "input"
         !---------------------------------------------------------------------!
         logical :: coupling_terms_file_exists = .false.
         !---------------------------------------------------------------------!
         if (reduced_mass <= 0) then
            call incorrect_value("reduced_mass", reduced_mass, input_unit)
         endif

         if ((relative_energy_flag /= 0).and.(relative_energy_flag /= 1)) then
            call incorrect_value("relative_energy_flag", relative_energy_flag, &
               input_unit)
         endif

         if (energy < 0) then
            call incorrect_value("energy", energy, input_unit)
         endif

         if (r_min <= 0) then
            call incorrect_value("r_min", r_min, input_unit)
         endif

         if (r_max <= 0) then
            call incorrect_value("r_max", r_max, input_unit)
         endif

         if (r_max <= r_min) then
            call incorrect_value("r_max/r_min", r_max / r_min, input_unit)
         endif

         if (steps <= 0.0_dp) then
            call incorrect_value("steps", steps, input_unit)
         endif

         if (potential_depth < 0.0_dp) then
            call incorrect_value("potential_depth", potential_depth, input_unit)
         endif

         if (jtot_min < 0) then
            call incorrect_value("jtot_min", jtot_min, input_unit)
         endif

         if (jtot_max < 0) then

            if (consecutive_blocks_threshold < 0) then
               call write_message("jtot_max < 0:")
               call incorrect_value("consecutive_blocks_threshold",            &
                  consecutive_blocks_threshold, input_unit)
            endif

            if (elastic_xs_threshold < 0) then
               call write_message("jtot_max < 0:")
               call incorrect_value("elastic_xs_threshold",                    &
                  elastic_xs_threshold, input_unit)
            endif

            if (inelastic_xs_threshold < 0) then
               call write_message("jtot_max < 0:")
               call incorrect_value("inelastic_xs_threshold",                  &
                  inelastic_xs_threshold, input_unit)
            endif

         else
         
            if (jtot_max < jtot_min) then
               call write_message("jtot_max is smaller than jtot_min")
               call incorrect_value("jtot_max/jtot_min",                       &
                  real(jtot_max/jtot_min, dp), input_unit)
            endif

         endif

         if (number_of_basis_levels <= 0) then
            call incorrect_value("number_of_basis_levels",                     &
               number_of_basis_levels, input_unit)
         endif

         if (relative_energy_flag == 1) then

            if (initial_level <= 0) then
               call write_message("relative_energy_flag = 1:")
               call incorrect_value("initial_level", initial_level, input_unit)
            endif

            if (initial_level > number_of_basis_levels) then
               call write_message("relative_energy_flag = 1:")
               call write_message("number_of_basis_levels = " //               &
                  trim(adjustl(integer_to_character(number_of_basis_levels))))
               call incorrect_value("initial_level > number_of_basis_levels",  &
                  initial_level, input_unit)
            endif

         endif

         if (number_of_r_points <= 0) then
            call incorrect_value("number_of_r_points", number_of_r_points,     &
               input_unit)
         endif

         if (number_of_legendre_indices <= 0) then
            call incorrect_value("number_of_legendre_indices",                 &
               number_of_legendre_indices, input_unit)
         endif

         if (total_number_of_coupling_terms <= 0) then
            call incorrect_value("total_number_of_coupling_terms",             &
               total_number_of_coupling_terms, input_unit)
         endif

         if (n_skip_lines < 0) then
            call incorrect_value("n_skip_lines", n_skip_lines, input_unit)
         endif

         select case(to_lowercase(trim(adjustl(coupling_terms_r_unit))))
            case("bohr")
            case("angstrom")
            case default
               call incorrect_value("coupling_terms_r_unit",                   &
                  coupling_terms_r_unit, input_unit)
         end select

         inquire(file = coupling_terms_file_name,                              &
            exist = coupling_terms_file_exists)
         if (.not.(coupling_terms_file_exists)) then
            call write_error(trim(adjustl(coupling_terms_file_name)) //        &
               " does not exist")
         endif

         if (print_level < 0) then
            call incorrect_value("print_level", print_level, input_unit)
         endif
         !---------------------------------------------------------------------!
      end subroutine check_namelist_input
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine check_namelist_basis
         !! Check variables read from namelist "basis"
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_
         !---------------------------------------------------------------------!
         do level_index_ = 1, number_of_basis_levels
            if (vib_levels(level_index_) < 0) then
               call incorrect_value("vib_levels(" //                           &
                  integer_to_character(level_index_) // ")",                   &
                  vib_levels(level_index_), input_unit)
            endif

            if (rot_levels(level_index_) < 0) then
               call incorrect_value("rot_levels(" //                           &
                  integer_to_character(level_index_) // ")",                   &
                  rot_levels(level_index_), input_unit)
            endif

            if (internal_energies(level_index_) < 0.0_dp) then
               call incorrect_value("internal_energies(" //                    &
                  integer_to_character(level_index_) // ")",                   &
                  internal_energies(level_index_), input_unit)
            endif
         enddo

      end subroutine check_namelist_basis
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine check_namelist_potential
         !! Check variables read from namelist "potential"
         !---------------------------------------------------------------------!
         integer(int32) :: legendre_index_, column_index_
         !---------------------------------------------------------------------!
         do legendre_index_ = 1, number_of_legendre_indices
            if (legendre_indices(legendre_index_) < 0) then
               call incorrect_value("legendre_indices(" //                     &
                  integer_to_character(legendre_index_) // ")",                &
                  legendre_indices(legendre_index_), input_unit)
            endif
         enddo

         do column_index_ = 1, total_number_of_coupling_terms
            if (vib_couplings(column_index_) < 0) then
               call incorrect_value("vib_couplings(" //                        &
                  integer_to_character(column_index_) // ")",                  &
                  vib_couplings(column_index_), input_unit)
            endif

            if (rot_couplings(column_index_) < 0) then
               call incorrect_value("rot_couplings(" //                        &
                  integer_to_character(column_index_) // ")",                  &
                  rot_couplings(column_index_), input_unit)
            endif

            if (vib_prime_couplings(column_index_) < 0) then
               call incorrect_value("vp1pes(" //                               &
                  integer_to_character(column_index_) // ")",                  &
                  vib_prime_couplings(column_index_), input_unit)
            endif

            if (rot_prime_couplings(column_index_) < 0) then
               call incorrect_value("rot_prime_couplings(" //                  &
                  integer_to_character(column_index_) // ")",                  &
                  rot_prime_couplings(column_index_), input_unit)
            endif
         enddo

      end subroutine check_namelist_potential
   !---------------------------------------------------------------------------!
end module input_validation
