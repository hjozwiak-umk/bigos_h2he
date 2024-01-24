program SCATTERING
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_header,                              &
      write_message, float_to_character, integer_to_character,                 &
      time_count_summary, no_open_channels_message
   use array_operations_mod, only: append, allocate_1d, allocate_2d
   use global_variables_mod
   use physics_utilities_mod, only: count_open_basis_levels,                   &
      save_open_basis_levels, units_conversion
   use input_reader_mod, only: read_input_file
   use radial_coupling_terms_mod, only: read_radial_coupling_terms,            &
      reduce_radial_coupling_terms, interpolate_radial_coupling_terms
   use channels_mod, only: set_number_of_channels, set_body_fixed_channels,    &
      set_space_fixed_channels, count_open_channels_in_block,                  &
      calculate_largest_wavevector, prepare_wavevector_array,                  &
      print_short_block_summary, print_channels
   use pes_matrix_mod, only: initialize_pes_matrix
   use propagator_mod, only: numerov
   use boundary_conditions_mod, only: calculate_sf_matrix_from_bf_matrix,      &
      calculate_k_matrix, calculate_s_matrix
   use unitarity_check_mod, only: unitarity_check, print_final_unitarity_warning
   use state_to_state_cross_sections_mod, only:initialize_cross_section_arrays,&
      calculate_state_to_state_cross_section,  add_cross_sections,             &
      print_largest_partial_cross_sections, print_cross_sections_for_jtot,     &
      print_final_cross_sections,  determine_largest_cross_sections,           &
      check_cross_section_thresholds, save_partial_xs_file_header,             &
      save_partial_xs_single_block
   use save_s_matrix_mod, only: save_s_matrix_file_header, save_s_matrix_block_info
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   integer(int32) :: number_of_channels, size_even, size_odd,                  &
      number_of_open_basis_levels, jtot_, parity_exponent, nsteps,             &
      number_of_open_channels, max_elastic_index, max_inelastic_index_1,       &
      max_inelastic_index_2, iopen, iopen2
   !---------------------------------------------------------------------------!
   integer(int32) :: count_blocks = 0
   integer(int32) :: consecutive_blocks_thresholddiag = 0
   integer(int32) :: consecutive_blocks_thresholdoff = 0
   !---------------------------------------------------------------------------!
   real(dp) :: largest_wavevector, wavpotential_depth,                         &
      max_elastic_cross_section, max_inelastic_cross_section, time_total_start,&
      time_total_stop, time_total, time_init_stop, time_init, time_jtot_start, &
      time_jtot_stop, time_jtot, time_parity_start,time_parity_stop,           &
      time_parity
   logical :: unitarity_block_check
   logical :: terminate = .false.
   integer, allocatable :: channel_indices(:), channels_omega_values(:),&
      channel_l_values(:), open_basis_levels(:), nonzero_terms_per_element(:), &
      nonzero_legendre_indices(:), smatcheckarr(:)
   real(dp), allocatable :: basis_wavevectors(:), block_wavevectors(:),        &
      nonzero_algebraic_coefficients(:), accumulated_cross_sections(:),        &
      partial_cross_sections_block(:), partial_cross_sections_jtot(:)
   real(dp), allocatable :: BF_log_der_matrix(:,:), SF_log_der_matrix(:,:),    &
      k_matrix(:,:), s_matrix_real(:,:), s_matrix_imag(:,:)
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   ! Initizalization: start the time count
   !---------------------------------------------------------------------------!
   call cpu_time(time_total_start)
   !---------------------------------------------------------------------------!
   ! Initialize fwigxjpf library
   !---------------------------------------------------------------------------!
   call fwig_table_init(2*100, 9)
   call fwig_temp_init(2*100)
   !---------------------------------------------------------------------------!
   ! Print the header
   !---------------------------------------------------------------------------!
   call write_header("main")
   !---------------------------------------------------------------------------!
   ! Read the input file
   !---------------------------------------------------------------------------!
   call read_input_file
   !---------------------------------------------------------------------------!
   ! S-matrix file: write input parameters and basis levels
   !---------------------------------------------------------------------------!
   call save_s_matrix_file_header
   !---------------------------------------------------------------------------!
   ! Prepare the file with the partial XS
   !---------------------------------------------------------------------------!
   if (print_partial_cross_sections) then
      call save_partial_xs_file_header
   endif
   !---------------------------------------------------------------------------!
   ! Convert units: starting now, everything is in atomic units
   !---------------------------------------------------------------------------!
   call units_conversion
   !---------------------------------------------------------------------------!
   ! Read the radial terms of the potential from external file
   !---------------------------------------------------------------------------!
   call read_radial_coupling_terms
   !---------------------------------------------------------------------------!
   ! Reduce matrix elements that are not needed
   !---------------------------------------------------------------------------!
   call reduce_radial_coupling_terms
   !---------------------------------------------------------------------------!
   ! Interpolate radial terms                                            
   !---------------------------------------------------------------------------!
   call interpolate_radial_coupling_terms
   !---------------------------------------------------------------------------!
   ! Search for energetically accessible levels and prepare the arrays that are
   ! needed in the calculations of the state-to-state XS                       
   !---------------------------------------------------------------------------!
   number_of_open_basis_levels = count_open_basis_levels()
   call save_open_basis_levels(number_of_open_basis_levels, open_basis_levels, &
      basis_wavevectors)
   !---------------------------------------------------------------------------!
   ! Initialize arrays that save the state-to-state cross-sections
   !---------------------------------------------------------------------------!
   call initialize_cross_section_arrays(number_of_open_basis_levels,           &
      partial_cross_sections_block, partial_cross_sections_jtot,               &
      accumulated_cross_sections)
   !---------------------------------------------------------------------------!
   ! Initialization is finished                                                
   !---------------------------------------------------------------------------!
   call cpu_time(time_init_stop)
   if (print_level.ge.2) call time_count_summary(time_total_start,             &
      time_init_stop, time_init, "Initialization completed in ")
   !---------------------------------------------------------------------------!
   ! Loop over total angular momentum
   !---------------------------------------------------------------------------!
   call write_header("jtot_loop")
   !---------------------------------------------------------------------------!
   do jtot_ = jtot_min,jtot_max,jtot_step
      !------------------------------------------------------------------------!
      call write_header("block", opt_integer_ = jtot_)
      !------------------------------------------------------------------------!
      call cpu_time(time_jtot_start)
      !------------------------------------------------------------------------!
      partial_cross_sections_jtot = 0
      call set_number_of_channels(jtot_, size_even, size_odd)
      !------------------------------------------------------------------------!
      do parity_exponent = 0,1
         !---------------------------------------------------------------------!
         call cpu_time(time_parity_start)
         !---------------------------------------------------------------------!
         select case(parity_exponent)
            case(0)
               number_of_channels = size_even
            case(1)
               number_of_channels = size_odd
         end select
         !---------------------------------------------------------------------!
         if (number_of_channels == 0) cycle
         !---------------------------------------------------------------------!
         ! Summary of the current block
         !---------------------------------------------------------------------!
         count_blocks = count_blocks+1
         if (print_level.ge.1) then
            call print_short_block_summary(jtot_, parity_exponent,             &
               count_blocks, number_of_channels)
         endif
         !---------------------------------------------------------------------!
         ! Prepare of the basis for each J/p block:
         ! channels_omega_values holds all values of omega (BF_)
         ! channel_l_values holds all values of l (SF_)
         ! channel_indices holds the indices which refer to the basis arrays:
         !   --   v1level/j1level/internal_energies
         !---------------------------------------------------------------------!
         call allocate_1d(channels_omega_values,number_of_channels)
         call allocate_1d(channel_l_values,number_of_channels)
         call allocate_1d(channel_indices,number_of_channels)
         !---------------------------------------------------------------------!
         ! Prepare channels_omega_values, channel_indices and channel_l_values
         !---------------------------------------------------------------------!
         call set_body_fixed_channels(jtot_, parity_exponent, channel_indices, &
            channels_omega_values)
         call set_space_fixed_channels(jtot_, parity_exponent, channel_l_values)
         !---------------------------------------------------------------------!
         ! Print the BF quantum numbers on screen
         !---------------------------------------------------------------------!
         if (print_level.ge.1) call print_channels(parity_exponent,            &
            channel_indices, channels_omega_values)
         !---------------------------------------------------------------------!
         ! Determine the number of open (energetically accessible) channels
         !---------------------------------------------------------------------!
         number_of_open_channels = count_open_channels_in_block(channel_indices)
         !---------------------------------------------------------------------!
         ! If there are no open channels, skip this block
         !---------------------------------------------------------------------!
         if (number_of_open_channels == 0) then
            call no_open_channels_message(count_blocks)
            cycle
         endif
         !---------------------------------------------------------------------!
         ! Array of wavevectors in given block (to save in the S-matrix file)
         !---------------------------------------------------------------------!
         call allocate_1d(block_wavevectors, number_of_open_channels)
         call prepare_wavevector_array(channel_indices, block_wavevectors)
         !---------------------------------------------------------------------!
         ! Determine the largest wavevector in the block
         !---------------------------------------------------------------------!
         largest_wavevector = calculate_largest_wavevector(channel_indices)
         !---------------------------------------------------------------------!
         ! Determine the number of steps on the intermolecular (R) grid
         ! This is done either directly (if r_step> 0)
         ! or through the number of steps per half de Broglie wavelength
         !---------------------------------------------------------------------!
         wavpotential_depth = dsqrt(2*reduced_mass*potential_depth)
         if (r_step<= 0) then
            nsteps = nint((r_max-r_min)/PI*((largest_wavevector+wavpotential_depth)*steps))
         else
            nsteps = nint((r_max-r_min)/r_step)+1
         endif
         !---------------------------------------------------------------------!
         ! Prepare the PES matrix
         !---------------------------------------------------------------------!
         call initialize_pes_matrix(channel_indices, channels_omega_values,    &
            nonzero_terms_per_element, nonzero_legendre_indices,               &
            nonzero_algebraic_coefficients)
         !---------------------------------------------------------------------!
         ! Allocate asymptotic body-fixed (BF) log-derivative matrix
         ! and call the propagator
         !---------------------------------------------------------------------!
         call allocate_2d(BF_log_der_matrix, number_of_channels, number_of_channels)
         call numerov(number_of_channels, channel_indices,                     &
            channels_omega_values, nonzero_terms_per_element,                  &
            nonzero_legendre_indices, nonzero_algebraic_coefficients, nsteps,  &
            jtot_, BF_log_der_matrix)
         !---------------------------------------------------------------------!
         ! Allocate asymptotic space-fixed (SF) log-derivative matrix and
         ! transform the BF result to the SF frame
         !---------------------------------------------------------------------!
         call allocate_2d(SF_log_der_matrix, number_of_channels, number_of_channels)
         call calculate_sf_matrix_from_bf_matrix(number_of_channels, jtot_,    &
            channel_indices, channels_omega_values, channel_l_values,          &
            BF_log_der_matrix, SF_log_der_matrix)
         !---------------------------------------------------------------------!
         ! Get the K-matrix from log-derivative matrix
         ! (Eq. 4 in "Solution of the coupled equations")
         !---------------------------------------------------------------------!
         call allocate_2d(k_matrix, number_of_open_channels, number_of_open_channels)
         call calculate_k_matrix(number_of_channels, SF_log_der_matrix,        &
            number_of_open_channels, channel_indices, channel_l_values,        &
            r_max, k_matrix)
         !---------------------------------------------------------------------!
         ! Get the S-matrix from the K-matrix
         ! (Eq. 12 in "Solution of the coupled equations")
         !---------------------------------------------------------------------!
         call allocate_2d(s_matrix_real,number_of_open_channels,number_of_open_channels)
         call allocate_2d(s_matrix_imag,number_of_open_channels,number_of_open_channels)
         call calculate_s_matrix(number_of_open_channels,k_matrix,s_matrix_real,s_matrix_imag)
         !---------------------------------------------------------------------!
         ! S-matrix is written to the binary S-matrix file
         !---------------------------------------------------------------------!
         call save_s_matrix_block_info(jtot_, parity_exponent,                 &
            number_of_open_channels, channel_indices, channel_l_values,        &
            block_wavevectors, s_matrix_real, s_matrix_imag)
         !---------------------------------------------------------------------!
         ! Check if the S-matrices are unitary
         !---------------------------------------------------------------------!
         call unitarity_check(number_of_open_channels, s_matrix_real,          &
            s_matrix_imag, unitarity_block_check)
         !---------------------------------------------------------------------!
         ! If the unitary is not fulfilled, keep the information about this block
         !---------------------------------------------------------------------!
         if (.not.(unitarity_block_check)) then
            call append(smatcheckarr, jtot_)
         endif
         !---------------------------------------------------------------------!
         ! Calculate all available cross-sections
         !---------------------------------------------------------------------!
         call calculate_state_to_state_cross_section(jtot_, open_basis_levels, &
            basis_wavevectors,s_matrix_real,s_matrix_imag,channel_indices,     &
            channel_l_values,partial_cross_sections_block)
         !---------------------------------------------------------------------!
         ! Print the results from this parity block to the partial XS file
         ! and add the calculated partial XS to the partial_cross_sections_jtot array
         !---------------------------------------------------------------------!
         if (print_partial_cross_sections) then
            call save_partial_xs_single_block(jtot_, count_blocks,             &
               open_basis_levels, partial_cross_sections_block)
         endif
         call add_cross_sections(number_of_open_basis_levels,                  &
            partial_cross_sections_block, partial_cross_sections_jtot)
         !---------------------------------------------------------------------!
         ! Check the time after each parity block:
         !---------------------------------------------------------------------!
         call cpu_time(time_parity_stop)
         if (print_level.ge.2) call time_count_summary(time_parity_start,      &
            time_parity_stop, time_parity, "Parity block completed in ")
         !---------------------------------------------------------------------!
         ! ... end of the loop over parity                                      
         !---------------------------------------------------------------------!
         call write_message(repeat("-", 90))
      enddo
      !------------------------------------------------------------------------!
      ! Add the cross-sections from this Jtot block:
      !------------------------------------------------------------------------!
      call add_cross_sections(number_of_open_basis_levels,                     &
            partial_cross_sections_jtot, accumulated_cross_sections)
      !------------------------------------------------------------------------!
      ! Determine the largest partial elastic/inelastic XS in this Jtot block:
      !------------------------------------------------------------------------!
      call determine_largest_cross_sections(open_basis_levels,                 &
         partial_cross_sections_jtot, max_elastic_cross_section,               &
         max_inelastic_cross_section, max_elastic_index,                       &
         max_inelastic_index_1, max_inelastic_index_2)
      !------------------------------------------------------------------------!
      call print_largest_partial_cross_sections(jtot_,                         &
         max_elastic_cross_section, max_inelastic_cross_section,               &
         max_elastic_index, max_inelastic_index_1, max_inelastic_index_2,      &
         open_basis_levels)
      !------------------------------------------------------------------------!
      if (jtot_max == 999999) then
         call check_cross_section_thresholds(max_elastic_cross_section,        &
            max_inelastic_cross_section, consecutive_blocks_thresholddiag,     &
            consecutive_blocks_thresholdoff, terminate)
      endif
      !------------------------------------------------------------------------!
      ! Check the time after each JTOT block:                                  
      !------------------------------------------------------------------------!
      call cpu_time(time_jtot_stop)
      !------------------------------------------------------------------------!
      ! Print all the XS after current JTOT block                              
      !------------------------------------------------------------------------!
      if (print_level.ge.3) then
         call print_cross_sections_for_jtot(jtot_, open_basis_levels,          &
            accumulated_cross_sections)
      endif
      !------------------------------------------------------------------------!
      if (print_level.ge.2) call time_count_summary(time_jtot_start,           &
         time_jtot_stop, time_jtot, "Total angular momentum block completed in ")
      !------------------------------------------------------------------------!
      ! terminate the loop if elastic_xs_threshold/inelastic_xs_threshold
      ! condition is satisfied
      !------------------------------------------------------------------------!
      if (terminate) exit
   enddo
   !---------------------------------------------------------------------------!
   call write_header("loop_terminated")
   call write_header("summary")
   !---------------------------------------------------------------------------!
   ! if for some JTOTs the S-matrix did not fulfill the unitary check,
   ! these are listed here
   !---------------------------------------------------------------------------!
   if (allocated(smatcheckarr)) then
      call print_final_unitarity_warning(smatcheckarr)
   endif
   !---------------------------------------------------------------------------!
   ! Print all the calculated XS                                                
   !---------------------------------------------------------------------------!
   call print_final_cross_sections(open_basis_levels, accumulated_cross_sections)
   !---------------------------------------------------------------------------!
   call fwig_temp_free();
   call fwig_table_free();
   !---------------------------------------------------------------------------!
   ! Stop the time count
   !---------------------------------------------------------------------------!
   call cpu_time(time_total_stop)
   call time_count_summary(time_total_start, time_total_stop, time_total,      &
      "Total CPU time: ")
   !---------------------------------------------------------------------------!
   close(s_matrix_unit)
   !---------------------------------------------------------------------------!
   if (print_partial_cross_sections) then
      close(partial_file_unit)
   endif
   !---------------------------------------------------------------------------!
end program SCATTERING
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
