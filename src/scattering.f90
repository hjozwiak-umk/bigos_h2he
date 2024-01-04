program SCATTERING
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use radial_coupling_terms_mod, only: read_radial_coupling_terms,            &
      reduce_radial_coupling_terms, interpolate_radial_coupling_terms
   use channels_mod, only: set_number_of_channels, set_body_fixed_channels,    &
      set_space_fixed_channels, count_open_channels_in_block,                  &
      calculate_largest_wavenumber, print_channels
   use coupling_matrix_mod
   use PROPAGATORS
   use boundary_conditions_mod, only: calculate_sf_matrix_from_bf_matrix,      &
      calculate_k_matrix, calculate_s_matrix
   use unitarity_check_mod, only: unitarity_check
   use state_to_state_cross_sections_mod, only:                                &
      calculate_state_to_state_cross_section,                                  &
      print_largest_partial_cross_sections, check_cross_section_thresholds
   use utility_functions_mod, only: write_header, file_io_status,              &
      write_message, float_to_character, integer_to_character, time_count_summary
   use array_operations_mod, only: append
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   character(len=200) :: err_message, partial_line, xs_line
   integer(int32) :: number_of_nonzero_coupling_matrix_elements,               &
      number_of_nonzero_coupling_coefficients, number_of_channels, size_even,  &
      size_odd, number_of_open_basis_levels, iblock, jtot_, parity_exponent,      &
      parity_exponenttmp, nsteps, number_of_open_channels, ncacdiag, ncacoff,     &
      omegamax, lmin, lmax, ltmp, lmat_len, len_even, len_odd, jinddiag,       &
      jindoff1, jindoff2, ij, ilevel, iomega, iopen, iopen2, isize_, isize_2,  &
      icheck, icount, icount2, io_status
   real(dp) :: largest_wavevector, wavvdepth, maxXSdiag, maxXSoff, time_total_start,       &
      time_total_stop, time_total, time_init_stop, time_init, time_jtot_start, &
      time_jtot_stop, time_jtot, time_parity_start,time_parity_stop,           &
      time_parity, time_coupling_start, time_coupling_stop, time_coupling
   logical :: unitarity_block_check, terminate
   integer, allocatable :: channel_indices(:), channels_omega_values(:),&
      channel_l_values(:), open_basis_levels(:), nonzero_terms_per_element(:),&
      nonzero_legendre_indices(:), smatcheckarr(:)
   real(dp), allocatable :: wv(:), open_basis_wavevectors(:),                  &
      nonzero_coupling_coefficients(:), xs_total(:), xs_block(:), xs_jtot(:)
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
   open(11, file=trim(smatrixfile), form='unformatted', iostat = io_status,    &
      iomsg = err_message)
   call file_io_status(io_status, err_message, 11, "o")
   !---------------------------------------------------------------------------!
   write(11) label, 2, nlevel, reducedmass
   write(11) (v1array(ilevel), j1array(ilevel), ilevel = 1, nlevel)
   write(11) (elevel(ilevel), ilevel = 1, nlevel)
   write(11) initial, energy
   !---------------------------------------------------------------------------!
   ! Prepare the file with the partial XS
   !---------------------------------------------------------------------------!
   if (parity_exponent == 1) then
      open(12, file=trim(partialfile),form='formatted',status='unknown',       &
         iostat = io_status, iomsg = err_message)
      call file_io_status(io_status, err_message, 12, "o")
      !------------------------------------------------------------------------!
      call write_message( "  jtot  iblock  v1_f  j1_f  <-  v1_i  j1_i'" //     &
         repeat(" ", 14) // "K.E." // repeat(" ", 16) // "XS", unit_ = 12)
      !------------------------------------------------------------------------!
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
   call allocate_1d(open_basis_levels, number_of_open_basis_levels)
   call allocate_1d(open_basis_wavevectors, number_of_open_basis_levels)
   call save_open_basis_levels(number_of_open_basis_levels, open_basis_levels, &
      open_basis_wavevectors)
   !---------------------------------------------------------------------------!
   ! xs array summed over all blocks
   !---------------------------------------------------------------------------!
   call allocate_1d(xs_total,                                                  &
      number_of_open_basis_levels*number_of_open_basis_levels)
   !---------------------------------------------------------------------------!
   ! xs array a single JTOT value
   !---------------------------------------------------------------------------!
   call allocate_1d(xs_jtot,                                                   &
      number_of_open_basis_levels*number_of_open_basis_levels)
   !---------------------------------------------------------------------------!
   ! xs array a single parity block
   !---------------------------------------------------------------------------!
   call allocate_1d(xs_block,                                                  &
      number_of_open_basis_levels*number_of_open_basis_levels)
   !---------------------------------------------------------------------------!
   ! Initialization is finished                                                
   !---------------------------------------------------------------------------!
   call cpu_time(time_init_stop)
   if (prntlvl.ge.2) call time_count_summary(time_total_start, time_init_stop, &
      time_init, "Initialization completed in ")
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   ! Prepare J-blocks
   !---------------------------------------------------------------------------!
   ! If JTOTMAX=-1 is called, iterate until convergence is achieved:
   ! this is managed by ncacdiag and ncacoff
   !---------------------------------------------------------------------------!
   ncacdiag  = 0
   ncacoff   = 0
   iblock    = 0
   terminate = .false.
   !---------------------------------------------------------------------------!
   ! Loop over total angular momentum
   !---------------------------------------------------------------------------!
   call write_message(repeat("*", 90))
   call write_message(repeat(" ", 28) // "*** Loop over JTOT: ***")
   !---------------------------------------------------------------------------!
   do jtot_ = jtotmin,jtotmax,jtotstep
      !------------------------------------------------------------------------!
      call write_header("block", opt_integer_ = jtot_)
      !------------------------------------------------------------------------!
      call cpu_time(time_jtot_start)
      !------------------------------------------------------------------------!
      xs_jtot = 0
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
         iblock = iblock+1
         if (prntlvl.ge.1) then
            call write_message("Block number: " // integer_to_character(iblock))
            call write_message("jtot: " //                                     &
               trim(adjustl(integer_to_character(jtot_))) // " parity: " //    &
               trim(adjustl(integer_to_character((-1)**parity_exponent) )))
            call write_message("Number of scattering channels: " //            &
               integer_to_character(number_of_channels))
         endif
         !---------------------------------------------------------------------!
         ! Prepare of the basis for each J/p block:
         ! channels_omega_values holds all values of omega (BF_)
         ! channel_l_values holds all values of l (SF_)
         ! channel_indices holds the indices which refer to the basis arrays:
         !   --   v1level/j1level/elevel
         !---------------------------------------------------------------------!
         call allocate_1d(channels_omega_values,number_of_channels)
         call allocate_1d(channel_l_values,number_of_channels)
         call allocate_1d(channel_indices,number_of_channels)
         !---------------------------------------------------------------------!
         ! Prepare channels_omega_values, channel_indices and channel_l_values
         !---------------------------------------------------------------------!
         call set_body_fixed_channels(jtot_, parity_exponent, channel_indices,   &
            channels_omega_values)
         call set_space_fixed_channels(jtot_, parity_exponent, channel_l_values)
         !---------------------------------------------------------------------!
         ! Print the BF quantum numbers on screen
         !---------------------------------------------------------------------!
         if (prntlvl.ge.1) call print_channels(parity_exponent,                &
            channel_indices, channels_omega_values)
         !---------------------------------------------------------------------!
         ! Determine the number of open (energetically accessible) channels
         !---------------------------------------------------------------------!
         number_of_open_channels = count_open_channels_in_block(channel_indices)
         !---------------------------------------------------------------------!
         ! If there are no open channels, skip this block
         !---------------------------------------------------------------------!
         if (number_of_open_channels == 0) then
            call write_message(repeat('-', 90))
            call write_message("No open channels for block no." //             &
               integer_to_character(iblock) )
            call write_message(repeat('-', 90))
            cycle
         endif
         !---------------------------------------------------------------------!
         ! Determine the largest wavevector in the block
         !---------------------------------------------------------------------!
         largest_wavevector = calculate_largest_wavenumber(channel_indices)
         !---------------------------------------------------------------------!
         ! Determine the number of steps on the intermolecular (R) grid
         ! This is done either directly (if dr > 0)
         ! or through the number of steps per half de Broglie wavelength
         !---------------------------------------------------------------------!
         wavvdepth = dsqrt(2*reducedmass*vdepth)
         if (dr <= 0) then
            nsteps = nint((Rmax-Rmin)/PI*((largest_wavevector+wavvdepth)*steps))
         else
            nsteps = nint((Rmax-Rmin)/dr)+1
         endif
         !---------------------------------------------------------------------!
         ! Prepare the coupling matrix
         !---------------------------------------------------------------------!
         call cpu_time(time_coupling_start)
         call check_nonzero_coupling_matrix_elements(channel_indices,   &
            channels_omega_values, number_of_nonzero_coupling_matrix_elements, &
            number_of_nonzero_coupling_coefficients)
         call allocate_1d(nonzero_terms_per_element,number_of_nonzero_coupling_matrix_elements)
         call allocate_1d(nonzero_coupling_coefficients,number_of_nonzero_coupling_coefficients)
         call allocate_1d(nonzero_legendre_indices,number_of_nonzero_coupling_coefficients)
         call prepare_coupling_matrix_elements(channel_indices,         &
            channels_omega_values, nonzero_terms_per_element,                  &
            nonzero_legendre_indices, nonzero_coupling_coefficients)
         if (prntlvl.ge.2) call print_coupling_matrix_elements_summary(        &
            number_of_channels, number_of_nonzero_coupling_matrix_elements,    &
            number_of_nonzero_coupling_coefficients)
         call cpu_time(time_coupling_stop)
         if (prntlvl.ge.2) call write_message("Calculations of the coupling "//&
            "matrix took " //  trim(adjustl(float_to_character(                &
            time_coupling_stop-time_coupling_start,"(E14.8)"))) // " seconds")
         !---------------------------------------------------------------------!
         ! Prepare the log-derivative matrix (Eqs. 6.29 and 6.43)
         ! and the K-matrix (Eq. 6.53)
         !---------------------------------------------------------------------!
         call allocate_2d(BF_log_der_matrix, number_of_channels, number_of_channels)
         call allocate_2d(SF_log_der_matrix, number_of_channels, number_of_channels)
         call allocate_2d(k_matrix, number_of_open_channels, number_of_open_channels)
         !---------------------------------------------------------------------!
         ! Call the propagator:
         !---------------------------------------------------------------------!
         call numerov(channel_indices, channels_omega_values,           &
            number_of_nonzero_coupling_matrix_elements,                        &
            number_of_nonzero_coupling_coefficients, nonzero_terms_per_element,&
            nonzero_legendre_indices, nonzero_coupling_coefficients, nsteps,   &
            number_of_channels, jtot_, BF_log_der_matrix)
         call write_message("Coupled equations were solved from " //           &
            trim(adjustl(float_to_character(Rmin, "(F10.4)")))// " a.u. to "// &
            trim(adjustl(float_to_character(Rmax, "(F10.4)")))// " a.u. in "// &
            trim(adjustl(integer_to_character(nsteps)))// " steps (dr = " //   &
            trim(adjustl(float_to_character((rmax - rmin) / dble(nsteps - 1),  &
            "(E14.8)"))) // " a.u.)")
         !---------------------------------------------------------------------!
         ! Transform the log-derivative matrix to the SF frame
         !---------------------------------------------------------------------!
         call calculate_sf_matrix_from_bf_matrix(number_of_channels, jtot_,               &
            channel_indices, channels_omega_values, channel_l_values,  &
            BF_log_der_matrix, SF_log_der_matrix)
         !---------------------------------------------------------------------!
         ! Get the K-matrix from log-derivative matrix (Eq. 6.53)
         !---------------------------------------------------------------------!
         call calculate_k_matrix(number_of_channels, SF_log_der_matrix,        &
            number_of_open_channels, channel_indices, channel_l_values,&
            rmax, k_matrix)
         !---------------------------------------------------------------------!
         ! Get the S-matrix from the K-matrix (Eq. 6.57)
         !---------------------------------------------------------------------!
         call allocate_2d(s_matrix_real,number_of_open_channels,number_of_open_channels)
         call allocate_2d(s_matrix_imag,number_of_open_channels,number_of_open_channels)
         call calculate_s_matrix(number_of_open_channels,k_matrix,s_matrix_real,s_matrix_imag)
         !---------------------------------------------------------------------!
         ! Array of wavevectors (necessary for the XS calculations)
         !---------------------------------------------------------------------!
         call allocate_1d(wv, number_of_open_channels)
         do iopen = 1, number_of_open_channels
            wv(iopen) = dsqrt((2*reducedmass*&
               (ETOTAL()-elevel(channel_indices(iopen)))))/bohrtoangstrom
         enddo
         !---------------------------------------------------------------------!
         ! S-matrix is written to the binary S-matrix file
         !---------------------------------------------------------------------!
         write(11) jtot_, parity_exponent, number_of_open_channels
         write(11) (channel_indices(iopen), channel_l_values(iopen),&
                    wv(iopen), iopen = 1, number_of_open_channels)
         write(11) ((s_matrix_real(iopen,iopen2), iopen2 = 1, iopen),&
                    iopen=1, number_of_open_channels)
         write(11) ((s_matrix_imag(iopen,iopen2), iopen2 = 1, iopen),&
                    iopen=1, number_of_open_channels)
         !---------------------------------------------------------------------!
         ! Check if the S-matrices are unitary
         !---------------------------------------------------------------------!
         call unitarity_check(number_of_open_channels,s_matrix_real,s_matrix_imag,unitarity_block_check)
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
            open_basis_wavevectors,s_matrix_real,s_matrix_imag,channel_indices,&
            channel_l_values,xs_block)
         !---------------------------------------------------------------------!
         ! Print the results from this parity block to the partial XS file
         ! and add the calculated partial XS to the xs_jtot array
         !---------------------------------------------------------------------!
         do icount = 1, number_of_open_basis_levels
            do icount2 = 1, number_of_open_basis_levels
               if (parity_exponent == 1) then
                  write(partial_line,                                          &
                     "(I6,2X,I6,2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8,2X,E16.8)")   &
                     jtot_,iblock,v1array(open_basis_levels(icount2)),         &
                     j1array(open_basis_levels(icount2)),                      &
                     v1array(open_basis_levels(icount)),                       &
                     j1array(open_basis_levels(icount)),                       &
                     (etotal()-elevel(open_basis_levels(icount)))*hartreetocm, &
                     xs_block((icount-1)*number_of_open_basis_levels+icount2)
                  call write_message(partial_line, unit_ = 12)
               endif
               xs_jtot((icount-1)*number_of_open_basis_levels+icount2) =       &
                  xs_jtot((icount-1)*number_of_open_basis_levels+icount2)      &
                  + xs_block((icount-1)*number_of_open_basis_levels+icount2)
            enddo
         enddo
         !---------------------------------------------------------------------!
         ! Check the time after each parity block:
         !---------------------------------------------------------------------!
         call cpu_time(time_parity_stop)
         if (prntlvl.ge.2) call time_count_summary(time_parity_start,          &
            time_parity_stop, time_parity, "Parity block completed in ")
         !---------------------------------------------------------------------!
         ! ... end of the loop over parity                                      
         !---------------------------------------------------------------------!
         call write_message(repeat(" ", 43) // "***")
      enddo
      !------------------------------------------------------------------------!
      ! Add the cross-sections from this Jtot block:
      !------------------------------------------------------------------------!
      do icount = 1, number_of_open_basis_levels
         do icount2 = 1, number_of_open_basis_levels
            xs_total((icount-1)*number_of_open_basis_levels+icount2) =         &
               xs_total((icount-1)*number_of_open_basis_levels+icount2)        &
               + xs_jtot((icount-1)*number_of_open_basis_levels+icount2)
         enddo
      enddo
      !------------------------------------------------------------------------!
      ! Determine the largest partial elastic/inelastic XS in this Jtot block:
      !------------------------------------------------------------------------!
      jinddiag  = 0
      jindoff1  = 0
      jindoff2  = 0
      maxXSdiag = 0.0_dp
      maxXSoff  = 0.0_dp
      do icount = 1, number_of_open_basis_levels
         do icount2 = 1, number_of_open_basis_levels
            if (open_basis_levels(icount2) == open_basis_levels(icount)) then
               if ((xs_jtot((icount-1)*number_of_open_basis_levels+icount2)).gt.maxXSdiag) then
                  maxXSdiag = xs_jtot((icount-1)*number_of_open_basis_levels+icount2)
                  jinddiag = icount
               endif
            else
               if ((xs_jtot((icount-1)*number_of_open_basis_levels+icount2)).gt.maxXSoff) then
                  maxXSoff = xs_jtot((icount-1)*number_of_open_basis_levels+icount2)
                  jindoff1 = icount
                  jindoff2 = icount2
               endif
            endif
         enddo
      enddo
      !-----------------------------------------------------------------------!
      call print_largest_partial_cross_sections(jtot_, maxXSdiag, maxXSoff, jinddiag,     &
         jindoff1, jindoff2, open_basis_levels)
      !-----------------------------------------------------------------------!
      if (jtotmax == 999999) then
         call check_cross_section_thresholds(maxXSdiag, maxXSoff, ncacdiag, ncacoff, terminate)
      endif
      !------------------------------------------------------------------------!
      ! Check the time after each JTOT block:                                  
      !------------------------------------------------------------------------!
      call cpu_time(time_jtot_stop)
      !------------------------------------------------------------------------!
      ! Print all the XS after current JTOT block                              
      !------------------------------------------------------------------------!
      if (prntlvl.ge.3) then
         call write_message("Cross sections for J: "//                         &
            trim(adjustl(integer_to_character(jtot_))) // " and energy: " //   &
            trim(adjustl(float_to_character(ETOTAL()*hartreetocm, "(F10.4)"))) &
            // " cm-1")
         call write_message("  v1_f  j1_f  <-  v1_i  j1_i" // repeat(" ", 14)  &
            // "K.E." // repeat(" ", 16) // "XS")
         do icount = 1, number_of_open_basis_levels
            do icount2 = 1, number_of_open_basis_levels
               write(xs_line, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8,2X,E16.8)")   &
                  v1array(open_basis_levels(icount2)),                         &
                  j1array(open_basis_levels(icount2)),                         &
                  v1array(open_basis_levels(icount)),                          &
                  j1array(open_basis_levels(icount)),                          &
                  (etotal()-elevel(open_basis_levels(icount)))*hartreetocm,    &
                  xs_total((icount-1)*number_of_open_basis_levels+icount2)
               call write_message(xs_line)
            enddo
         enddo
      endif
      !------------------------------------------------------------------------!
      if (prntlvl.ge.2) call time_count_summary(time_jtot_start,               &
         time_jtot_stop, time_jtot, "JTOT block completed in ")
      !------------------------------------------------------------------------!
      ! terminate the loop if dtol/otol condition is satisfied
      !------------------------------------------------------------------------!
      if (terminate) exit
   enddo

   call write_message(repeat('*', 90))
   call write_message(repeat(" ", 31) // "Loop over JTOT finished")

   call write_message(repeat('*', 90))
   call write_message("*" // repeat(" ", 40) // "SUMMARY" // repeat(" ", 41)   &
      // "*")
   call write_message(repeat('*', 90))
   !---------------------------------------------------------------------------!
   ! if for some JTOTs the S-matrix did not fulfill the unitary check,
   ! these are listed here
   !---------------------------------------------------------------------------!
   if (allocated(smatcheckarr)) then
      print *
      call write_message(repeat("-", 90))
      call write_message(repeat(" ", 37) // "*** WARNING ***")
      call write_message(repeat("-", 90))
      call write_message("Check unitarity of the S-matrix in the following "// &
         "JTOT blocks:")
      do icheck=1, size(smatcheckarr)
         call write_message("JTOT:" // repeat(" ", 8) // integer_to_character( &
            smatcheckarr(icheck)))
      enddo
      call write_message(repeat("-", 90))
      print *
   endif
   !---------------------------------------------------------------------------!
   ! Print all the calculated XS                                                
   !---------------------------------------------------------------------------!
   call write_message("Final state-to-state XS")
   call write_message("  v1_f  j1_f  <-  v1_i  j1_i" // repeat(" ", 14) //     &
      "K.E." // repeat(" ", 16) // "XS")
   do icount = 1, number_of_open_basis_levels
      do icount2 = 1, number_of_open_basis_levels
         write(xs_line, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8,2X,E16.8)")         &
            v1array(open_basis_levels(icount2)),                               &
            j1array(open_basis_levels(icount2)),                               &
            v1array(open_basis_levels(icount)),                                &
            j1array(open_basis_levels(icount)),                                &
            (etotal()-elevel(open_basis_levels(icount)))*hartreetocm,          &
            xs_total((icount-1)*number_of_open_basis_levels+icount2)
         call write_message(xs_line)
      enddo
   enddo
   !---------------------------------------------------------------------------!
   call fwig_temp_free();
   call fwig_table_free();
   !---------------------------------------------------------------------------!
   ! Stop the time count
   !---------------------------------------------------------------------------!
   call cpu_time(time_total_stop)
   call time_count_summary(time_total_start, time_total_stop, time_total,      &
      "Total CPU time: ")
   close(11)
   close(12)
   !---------------------------------------------------------------------------!
end program SCATTERING
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
