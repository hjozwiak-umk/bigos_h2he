module input_reader_mod
   !! this module provides following functions and subroutines:
   !! 1. input_file - reads the input file prepared by the user
   !! 2. input_check - checks the variables supplied in the input file
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: file_io_status, write_error, write_message,&
      incorrect_value, integer_to_character, float_to_character
   use array_operations_mod, only: allocate_1d, allocate_2d, allocate_3d
   use data_mod
   use physics_utilities_mod, only: total_energy
   use input_validation
!------------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains

      subroutine read_input_file
         !! reads the input file prepared by the user using NAMELIST feature
         !! the code uses 3 namelists: input, basis and potential
         !------------------------------------------------------------------------!
         character(len = 200):: err_message
         integer(int32) :: level_index_1_, level_index_2_, icoupl, icol, il, io_status
         !------------------------------------------------------------------------!
         namelist / INPUT / label, reduced_mass, relative_energy_flag, energy, &
            jtotmin, jtotmax, jtotstep, rmin, rmax, dr, steps, vdepth,         &
            consecutive_blocks_threshold, elastic_xs_threshold,                &
            inelastic_xs_threshold, nlevel, initial, nr, nterms,               &
            total_number_of_coupling_terms, n_skip_lines, iunits,              &
            potentialfile, smatrixfile, print_partial_cross_sections,          &
            partialfile, prntlvl
         namelist / BASIS / v1array, j1array, elevel
         namelist / POTENTIAL / l1tab, v1pes, j1pes, v1ppes, j1ppes
!------------------------------------------------------------------------------!
         open(unit=5, action='read', form='formatted', access='sequential',    &
            status = 'old', iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'o')
         !---------------------------------------------------------------------!
         call write_message("-- Reading input file...")
         !---------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! Read the input namelist:                                                     !
!------------------------------------------------------------------------------!
         read(unit=5, nml=INPUT, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
!------------------------------------------------------------------------------!
! Check if the variables from input namelist are supplied correctly:           !
!------------------------------------------------------------------------------!
         call check_namelist_input

         if (jtotmax.eq.-1) jtotmax = 999999

         call allocate_1d(v1array,nlevel)
         call allocate_1d(j1array,nlevel)
         call allocate_1d(elevel,nlevel)

         call allocate_1d(l1tab,nterms)
         call allocate_1d(v1pes,total_number_of_coupling_terms)
         call allocate_1d(v1ppes,total_number_of_coupling_terms)
         call allocate_1d(j1pes,total_number_of_coupling_terms)
         call allocate_1d(j1ppes,total_number_of_coupling_terms)

         select case(iunits)
            case(0)
               radial_term_distance_converter = 1.0_dp
            case(1)
               radial_term_distance_converter = bohrtoangstrom
         end select

         radial_term_energy_converter = 1.0_dp / hartreetocm
!------------------------------------------------------------------------------!
! Read the basis namelist & check if the values were supplied correctly:       !
!------------------------------------------------------------------------------!
         read(unit=5, nml=BASIS, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
         call check_namelist_basis
!------------------------------------------------------------------------------!
! If itype = 2/4 the code reads all the total_number_of_coupling_terms coupling terms, but some of   !
! them will not be used in the calculations. Here, the code prepares           !
! the arrays of minimal_number_of_coupling_terms size, that will hold only the necessary terms           !
!------------------------------------------------------------------------------!
         minimal_number_of_coupling_terms = nlevel * (nlevel + 1) / 2

         call allocate_1d(reduced_j1pes,minimal_number_of_coupling_terms)
         call allocate_1d(reduced_j1ppes,minimal_number_of_coupling_terms)
         call allocate_1d(reduced_v1pes,minimal_number_of_coupling_terms)
         call allocate_1d(reduced_v1ppes,minimal_number_of_coupling_terms)

         icoupl = 0

         do level_index_1_ = 1, nlevel
            do level_index_2_ = level_index_1_, nlevel
               icoupl = icoupl + 1
               reduced_v1pes(icoupl)  = v1array(level_index_1_)
               reduced_j1pes(icoupl)  = j1array(level_index_1_)
               reduced_v1ppes(icoupl) = v1array(level_index_2_)
               reduced_j1ppes(icoupl) = j1array(level_index_2_)
            enddo
         enddo
!------------------------------------------------------------------------------!
! Read the potential namelist & check if the values were supplied correctly:   !
!------------------------------------------------------------------------------!
         read(unit=5, nml=POTENTIAL, iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, 5, 'r')
         call check_namelist_potential

         close(5)
!------------------------------------------------------------------------------!
! Prepare the arrays that are needed for interpolation of the coupling terms:  !
!------------------------------------------------------------------------------!
         call allocate_1d(rmat,nr)

         call allocate_3d(read_vmat3D,nr,nterms,total_number_of_coupling_terms)
         call allocate_3d(vmat3D,nr,nterms,minimal_number_of_coupling_terms)
         call allocate_3d(bmat3D,nr,nterms,minimal_number_of_coupling_terms)
         call allocate_3d(cmat3D,nr,nterms,minimal_number_of_coupling_terms)
         call allocate_3d(dmat3D,nr,nterms,minimal_number_of_coupling_terms)
!------------------------------------------------------------------------------!
! Summarize the input parameters:                                              !
!------------------------------------------------------------------------------!
         call input_summary

      end subroutine read_input_file
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
   subroutine input_summary
      !! summarize the input parameters for the current run
      !------------------------------------------------------------------------!
      integer(int32) :: level_index_1_
      !------------------------------------------------------------------------!
      call write_message(" - User-supplied label: " // label)
      !------------------------------------------------------------------------!
      call write_message(" - Reduced mass: " //                                &
         trim(adjustl(float_to_character(reduced_mass, "(F10.4)"))) // " a.m.u.")
      !------------------------------------------------------------------------!
      call write_message(" - Energy levels in the basis set:")
      !------------------------------------------------------------------------!
      call write_message("   v       j            Energy (cm^{-1})")
      !------------------------------------------------------------------------!
      do level_index_1_ = 1,nlevel
         write(*,"(I4,4X,I4,16X,F12.4)") v1array(level_index_1_),              &
            j1array(level_index_1_), elevel(level_index_1_)
      enddo
      !------------------------------------------------------------------------!
      if (jtotmax.ne.999999) then
         !---------------------------------------------------------------------!
         call write_message(" - The equations will be solved " //              &
            "for total angular momentum J from " //                            &
            trim(adjustl(integer_to_character(jtotmin))) // " to "             &
            // trim(adjustl(integer_to_character(jtotmax))) // " with step "// &
            trim(adjustl(integer_to_character(jtotstep))))
         !---------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------!
         call write_message(" - The loop over JTOT will be performed from " // &
            trim(adjustl(integer_to_character(jtotmin))) // " with step " //   &
            trim(adjustl(integer_to_character(jtotstep))) // " until " //      &
            trim(adjustl(integer_to_character(consecutive_blocks_threshold)))  &
            // " consecutive")
          call write_message("   total angular momentum blocks contribute less than")
          call write_message("   - " //                                        &
            trim(adjustl(float_to_character(elastic_xs_threshold, "(E10.4)"))) &
            // " A^2 to the elastic XS")
          call write_message("   - " //                                        &
            trim(adjustl(float_to_character(inelastic_xs_threshold, "(E10.4)")))&
            // " A^2 to the inelastic XS")
         !---------------------------------------------------------------------!
      endif
      !------------------------------------------------------------------------!
      if (relative_energy_flag.eq.0) then
         !---------------------------------------------------------------------!
         call write_message(" - The calculations will be performed for the total energy equal to "&
            // trim(adjustl(float_to_character(total_energy(), "(F10.4)")))//" cm-1")
         !---------------------------------------------------------------------!
      else if(relative_energy_flag.eq.1) then
         !---------------------------------------------------------------------!
         call write_message(" - Relative kinetic energy of the colliding system: " //&
            trim(adjustl(float_to_character(energy, "(F10.4)"))) // " cm-1")
         !---------------------------------------------------------------------!
         call write_message(" - The kinetic energy is calculated with respect to the" //&
            " v = " // trim(adjustl(integer_to_character(v1array(initial)))) //&
            " j = " // trim(adjustl(integer_to_character(j1array(initial)))) //&
            " level")
         call write_message("   with the internal energy of " //               &
            trim(adjustl(float_to_character(elevel(initial), "(F10.4)")))      &
            // " cm-1.")
         !---------------------------------------------------------------------!
         call write_message(" - This gives the total energy equal to " // &
            trim(adjustl(float_to_character(total_energy(), "(F10.4)"))) // " cm-1")
         !---------------------------------------------------------------------!
      endif
      !------------------------------------------------------------------------!
      if (print_partial_cross_sections) then
         call write_message(" - Partial cross sections will be saved into " // partialfile )
      endif
      !------------------------------------------------------------------------!
      call write_message(" - S-matrix elements will be saved into " // smatrixfile )
      !------------------------------------------------------------------------!
   end subroutine input_summary
!------------------------------------------------------------------------------!
end module input_reader_mod
