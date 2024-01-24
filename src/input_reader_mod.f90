module input_reader_mod
   !! This module provides following functions and subroutines:
   !! -- input_file - reads the input file prepared by the user
   !! -- input_summary - interprets and prints the input parameters
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: file_io_status, write_message,             &
      integer_to_character, float_to_character, to_lowercase
   use array_operations_mod, only: allocate_1d, allocate_3d
   use global_variables_mod
   use physics_utilities_mod, only: total_energy
   use input_validation
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: read_input_file
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine read_input_file
         !! reads the input file prepared by the user using NAMELIST feature
         !! the code uses 3 namelists: "input", "basis" and "potential"
         !---------------------------------------------------------------------!
         character(len = 200):: err_message
         integer(int32) :: level_index_1_, level_index_2_, coupling_index_,    &
            io_status
         !---------------------------------------------------------------------!
         namelist / INPUT / label, reduced_mass, relative_energy_flag, energy, &
            jtot_min, jtot_max, jtot_step, r_min, r_max, r_step, steps,        &
            potential_depth, consecutive_blocks_threshold,                     &
            elastic_xs_threshold, inelastic_xs_threshold,                      &
            number_of_basis_levels, initial_level, number_of_r_points,         &
            number_of_legendre_indices, total_number_of_coupling_terms,        &
            n_skip_lines, coupling_terms_r_unit, coupling_terms_file_name,     &
            s_matrix_file_name, print_partial_cross_sections,                  &
            partial_xs_file_name, print_level
         namelist / BASIS / vib_levels, rot_levels, internal_energies
         namelist / POTENTIAL / legendre_indices, vib_couplings, rot_couplings,&
            vib_prime_couplings, rot_prime_couplings
         !---------------------------------------------------------------------!
         open(unit=input_unit, action='read', form='formatted',                &
            access='sequential', status = 'old', iostat = io_status,           &
            iomsg = err_message)
         call file_io_status(io_status, err_message, input_unit, 'o')
         !---------------------------------------------------------------------!
         call write_message("-- Reading input file...")
         !---------------------------------------------------------------------!
         ! Read the "input" namelist:
         !---------------------------------------------------------------------!
         read(unit = input_unit, nml = INPUT, iostat = io_status,              &
            iomsg = err_message)
         call file_io_status(io_status, err_message, input_unit, 'r')
         !---------------------------------------------------------------------!
         ! Check if the variables from "input" namelist are supplied correctly:
         !---------------------------------------------------------------------!
         call check_namelist_input
         !---------------------------------------------------------------------!
         ! Set values/sizes of related variables
         !---------------------------------------------------------------------!
         if (jtot_max.eq.-1) jtot_max = 999999
         !---------------------------------------------------------------------!
         call allocate_1d(vib_levels,number_of_basis_levels)
         call allocate_1d(rot_levels,number_of_basis_levels)
         call allocate_1d(internal_energies,number_of_basis_levels)
         !---------------------------------------------------------------------!
         call allocate_1d(legendre_indices,number_of_legendre_indices)
         call allocate_1d(vib_couplings,total_number_of_coupling_terms)
         call allocate_1d(vib_prime_couplings,total_number_of_coupling_terms)
         call allocate_1d(rot_couplings,total_number_of_coupling_terms)
         call allocate_1d(rot_prime_couplings,total_number_of_coupling_terms)
         !---------------------------------------------------------------------!
         select case(to_lowercase(trim(adjustl(coupling_terms_r_unit))))
            case("bohr")
               radial_term_distance_converter = 1.0_dp
            case("angstrom")
               radial_term_distance_converter = bohr_to_angstrom
         end select
         !---------------------------------------------------------------------!
         radial_term_energy_converter = 1.0_dp / hartree_to_cm
         !---------------------------------------------------------------------!
         ! Read the "basis" namelist & check the values
         !---------------------------------------------------------------------!
         read(unit = input_unit, nml = BASIS, iostat = io_status,              &
            iomsg = err_message)
         call file_io_status(io_status, err_message, input_unit, 'r')
         call check_namelist_basis
         !---------------------------------------------------------------------!
         ! The code reads all the total_number_of_coupling_terms coupling terms,
         ! but some of them will not be used in the calculations. Here, the
         ! code prepares the arrays of minimal_number_of_coupling_terms size,
         ! that will hold only the necessary terms           !
         !---------------------------------------------------------------------!
         minimal_number_of_coupling_terms                                      &
            = number_of_basis_levels * (number_of_basis_levels + 1) / 2
         !---------------------------------------------------------------------!
         call allocate_1d(reduced_rot_couplings,                               &
            minimal_number_of_coupling_terms)
         call allocate_1d(reduced_rot_prime_couplings,                         &
            minimal_number_of_coupling_terms)
         call allocate_1d(reduced_vib_couplings,                               &
            minimal_number_of_coupling_terms)
         call allocate_1d(reduced_vib_prime_couplings,                         &
            minimal_number_of_coupling_terms)
         !---------------------------------------------------------------------!
         coupling_index_ = 0
         do level_index_1_ = 1, number_of_basis_levels
            do level_index_2_ = level_index_1_, number_of_basis_levels
               coupling_index_ = coupling_index_ + 1
               reduced_vib_couplings(coupling_index_)                          &
                  = vib_levels(level_index_1_)
               reduced_rot_couplings(coupling_index_)                          &
                  = rot_levels(level_index_1_)
               reduced_vib_prime_couplings(coupling_index_)                    &
                  = vib_levels(level_index_2_)
               reduced_rot_prime_couplings(coupling_index_)                    &
                  = rot_levels(level_index_2_)
            enddo
         enddo
         !---------------------------------------------------------------------!
         ! Read the "potential" namelist & check the values
         !---------------------------------------------------------------------!
         read(unit = input_unit, nml=POTENTIAL, iostat = io_status,            &
            iomsg = err_message)
         call file_io_status(io_status, err_message, input_unit, 'r')
         call check_namelist_potential
         !---------------------------------------------------------------------!
         close(input_unit)
         !---------------------------------------------------------------------!
         ! Set values/sizes of related variables
         !---------------------------------------------------------------------!
         call allocate_1d(r_grid,number_of_r_points)
         call allocate_3d(tabulated_coupling_terms,number_of_r_points,         &
            number_of_legendre_indices, total_number_of_coupling_terms)
         call allocate_3d(coupling_terms,number_of_r_points,                   &
            number_of_legendre_indices, minimal_number_of_coupling_terms)
         call allocate_3d(coupling_terms_b_coeffs,number_of_r_points,          &
            number_of_legendre_indices, minimal_number_of_coupling_terms)
         call allocate_3d(coupling_terms_c_coeffs,number_of_r_points,          &
            number_of_legendre_indices, minimal_number_of_coupling_terms)
         call allocate_3d(coupling_terms_d_coeffs,number_of_r_points,          &
            number_of_legendre_indices, minimal_number_of_coupling_terms)
         !---------------------------------------------------------------------!
         ! Summarize the input parameters
         !---------------------------------------------------------------------!
         call input_summary
         !---------------------------------------------------------------------!
      end subroutine read_input_file
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine input_summary
         !! summarize the input parameters for the current run
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_1_
         !---------------------------------------------------------------------!
         call write_message(" - User-supplied label: " // label)
         !---------------------------------------------------------------------!
         call write_message(" - Reduced mass: " //                             &
            trim(adjustl(float_to_character(reduced_mass, "(F10.4)"))) //      &
            " a.m.u.")
         !---------------------------------------------------------------------!
         call write_message(" - Energy levels in the basis set:")
         !---------------------------------------------------------------------!
         call write_message("   v       j            Energy (cm^{-1})")
         !---------------------------------------------------------------------!
         do level_index_1_ = 1,number_of_basis_levels
            write(*,"(I4,4X,I4,16X,F12.4)") vib_levels(level_index_1_),        &
               rot_levels(level_index_1_), internal_energies(level_index_1_)
         enddo
         !---------------------------------------------------------------------!
         if (jtot_max.ne.999999) then
            !------------------------------------------------------------------!
            call write_message(" - The equations will be solved " //           &
               "for total angular momentum J from " //                         &
               trim(adjustl(integer_to_character(jtot_min))) // " to "         &
               // trim(adjustl(integer_to_character(jtot_max))) //             &
               " with step "// trim(adjustl(integer_to_character(jtot_step))))
            !------------------------------------------------------------------!
         else
            !------------------------------------------------------------------!
            call write_message(" - The loop over JTOT will be performed from " &
               // trim(adjustl(integer_to_character(jtot_min))) //             &
               " with step " // trim(adjustl(                                  &
                  integer_to_character(jtot_step))) // " until " //            &
               trim(adjustl(integer_to_character(                              &
                  consecutive_blocks_threshold))) // " consecutive")
             call write_message("   total angular momentum blocks contribute " &
               //"less than")
             call write_message("   - " //                                     &
               trim(adjustl(float_to_character(elastic_xs_threshold,           &
               "(E10.4)"))) // " A^2 to the elastic XS")
             call write_message("   - " //                                     &
               trim(adjustl(float_to_character(inelastic_xs_threshold,         &
               "(E10.4)"))) // " A^2 to the inelastic XS")
            !------------------------------------------------------------------!
         endif
         !---------------------------------------------------------------------!
         if (relative_energy_flag.eq.0) then
            !------------------------------------------------------------------!
            call write_message(" - The calculations will be performed for the" &
               // " total energy equal to "                                    &
               // trim(adjustl(float_to_character(total_energy(), "(F10.4)"))) &
               //" cm-1")
            !------------------------------------------------------------------!
         else if(relative_energy_flag.eq.1) then
            !------------------------------------------------------------------!
            call write_message(" - Relative kinetic energy of the colliding"   &
               // " system: " // trim(adjustl(float_to_character(energy,       &
               "(F10.4)"))) // " cm-1")
            !------------------------------------------------------------------!
            call write_message(" - The kinetic energy is calculated with" //   &
               " respect to the v = " // trim(adjustl(integer_to_character(    &
               vib_levels(initial_level)))) // " j = " //                      &
               trim(adjustl(integer_to_character(rot_levels(initial_level))))  &
               // " level")
            call write_message("   with the internal energy of " //            &
               trim(adjustl(float_to_character(internal_energies(              &
               initial_level), "(F10.4)"))) // " cm-1.")
            !------------------------------------------------------------------!
            call write_message(" - This gives the total energy equal to " //   &
               trim(adjustl(float_to_character(total_energy(), "(F10.4)"))) // &
               " cm-1")
            !------------------------------------------------------------------!
         endif
         !---------------------------------------------------------------------!
         if (print_partial_cross_sections) then
            call write_message(" - Partial cross sections will be saved into " &
               // partial_xs_file_name )
         endif
         !---------------------------------------------------------------------!
         call write_message(" - S-matrix elements will be saved into " //      &
            s_matrix_file_name )
         !---------------------------------------------------------------------!
      end subroutine input_summary
   !---------------------------------------------------------------------------!
end module input_reader_mod
