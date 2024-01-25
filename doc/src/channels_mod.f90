module channels_mod
   !! This module provides subroutines that set the number of channels in the
   !! block, save quantum numbers for each channel (both in body- and
   !! space-fixed cases) and print quantum numbers on screen.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_message, write_warning, &
      integer_to_character, float_to_character
   use global_variables_mod
   use physics_utilities_mod, only: is_open, wavevector_squared_from_energy
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: set_number_of_channels, set_body_fixed_channels,                  &
      set_space_fixed_channels, count_open_channels_in_block,                  &
      prepare_wavevector_array, calculate_largest_wavevector,                  &
      calculate_number_of_steps, print_short_block_summary, print_channels
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine set_number_of_channels(total_angular_momentum_,               &
         number_of_channels_even_parity_block, number_of_channels_odd_parity_block)
         !! determine the number of scattering channels in each parity block 
         !! for given total angular momentum in both body-fixed and
         !! space-fixed frames
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum         
         integer(int32), intent(out) :: number_of_channels_even_parity_block
            !! number of channels in the p = 1 (even parity) block
         integer(int32), intent(out) :: number_of_channels_odd_parity_block
            !! number of channels in the p = -1 (odd parity) block
         !---------------------------------------------------------------------!
         integer(int32) :: number_of_channels_even_parity_block_sf,            &
            number_of_channels_odd_parity_block_sf
         !---------------------------------------------------------------------!
         ! body-fixed frame
         !---------------------------------------------------------------------!
         call calculate_number_of_channels_body_fixed(total_angular_momentum_, &
            number_of_channels_even_parity_block, number_of_channels_odd_parity_block)
         !---------------------------------------------------------------------!
         ! space-fixed frame
         !---------------------------------------------------------------------!
         call calculate_number_of_channels_space_fixed(total_angular_momentum_, &
            number_of_channels_even_parity_block_sf, number_of_channels_odd_parity_block_sf)
         !---------------------------------------------------------------------!
         ! Check if the number of channels is the same
         !---------------------------------------------------------------------!
         call check_number_of_channels(number_of_channels_even_parity_block,   &
            number_of_channels_even_parity_block_sf, "even")
         call check_number_of_channels(number_of_channels_odd_parity_block,    &
            number_of_channels_odd_parity_block_sf, "odd")
         !---------------------------------------------------------------------!
      end subroutine set_number_of_channels
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_number_of_channels_body_fixed(                      &
         total_angular_momentum_, number_of_channels_even_parity_block,        &
         number_of_channels_odd_parity_block)
         !! calculate number of channels in even and odd parity
         !! blocks in the body-fixed frame;
         !! in principle, \\(\bar{\Omega}\in \langle 0, \mathrm{min}(j, J)\)),
         !! but number of channels additionally depends on the
         !! sign of \\( p (-1)^{J} \\): channels with
         !! \\(\bar{\Omega}=0\\) values enter blocks with
         !! \\( p (-1)^{J} = + 1 \\) _only_;
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum      
         integer(int32), intent(inout) :: number_of_channels_even_parity_block
            !! number of channels in the p = 1 (even parity) block
         integer(int32), intent(inout) :: number_of_channels_odd_parity_block
            !! number of channels in the p = -1 (odd parity) block
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_, omega_max_
         !---------------------------------------------------------------------!
         number_of_channels_even_parity_block = 0
         number_of_channels_odd_parity_block  = 0
         do level_index_ = 1, number_of_basis_levels
            omega_max_ = min(rot_levels(level_index_), total_angular_momentum_)

            call update_channel_counts_body_fixed(omega_max_,                  &
               total_angular_momentum_, number_of_channels_even_parity_block,  &
               number_of_channels_odd_parity_block)

         end do
         !---------------------------------------------------------------------!
      end subroutine calculate_number_of_channels_body_fixed
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine update_channel_counts_body_fixed(omega_max_,                  &
         total_angular_momentum_, number_of_channels_even_parity_block,        &
         number_of_channels_odd_parity_block)
         !! updates number_of_channels_even and number_of_channels_odd
         !! in the body-fixed frame for given \\(J\\) and \\(\bar{Omega}_{max}\\)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: omega_max_
            !! largest value of \\(\bar{\Omega}\\), for given
            !! rotational and total angular momenta
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(inout) :: number_of_channels_even_parity_block
            !! number of channels for the p = 1 block
         integer(int32), intent(inout) :: number_of_channels_odd_parity_block
            !! number of channels for the p = -1 block
         !---------------------------------------------------------------------!
         if (mod(total_angular_momentum_, 2) == 0) then
            !------------------------------------------------------------------!
            ! Even J: channels with Omega = 0 only count in the even parity block
            !------------------------------------------------------------------!
            number_of_channels_even_parity_block                               &
               = number_of_channels_even_parity_block + omega_max_ + 1
            number_of_channels_odd_parity_block                                &
               = number_of_channels_odd_parity_block + omega_max_
         else
            !------------------------------------------------------------------!
            ! Odd J: channels with Omega = 0 only count in the odd parity block
            !------------------------------------------------------------------!
            number_of_channels_odd_parity_block                                &
               = number_of_channels_odd_parity_block + omega_max_ + 1
            number_of_channels_even_parity_block                               &
               = number_of_channels_even_parity_block + omega_max_
         endif
         !---------------------------------------------------------------------!
      end subroutine update_channel_counts_body_fixed
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_number_of_channels_space_fixed(                     &
         total_angular_momentum_, number_of_channels_even_parity_block,        &
         number_of_channels_odd_parity_block)
         !! calculate number of channels in even and odd parity
         !! blocks in the space-fixed frame based on available
         !! values of orbital angular momentum:
         !! \\( l \in \langle |j-J|, j+J \rangle \\);
         !! parity is defined as \\(p= (-1)^{j+l}\\)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum      
         integer(int32), intent(inout) :: number_of_channels_even_parity_block
            !! number of channels in the p = 1 (even parity) block
         integer(int32), intent(inout) :: number_of_channels_odd_parity_block
            !! number of channels in the p = -1 (odd parity) block
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_, l_min_, l_max_
         !---------------------------------------------------------------------!
         number_of_channels_even_parity_block = 0
         number_of_channels_odd_parity_block  = 0
         do level_index_ = 1, number_of_basis_levels
            l_min_ = abs(total_angular_momentum_ - rot_levels(level_index_))
            l_max_ = total_angular_momentum_ + rot_levels(level_index_)
            call update_channel_counts_space_fixed(l_min_, l_max_,             &
               level_index_, number_of_channels_even_parity_block,             &
               number_of_channels_odd_parity_block)
         end do
         !---------------------------------------------------------------------!
      end subroutine calculate_number_of_channels_space_fixed
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine update_channel_counts_space_fixed(l_min_, l_max_,             &
         level_index_, number_of_channels_even_parity_block,                   &
         number_of_channels_odd_parity_block)
         !! updates number_of_channels_even and number_of_channels_odd
         !! in the space-fixed frame for given
         !! range of orbital angular momentum
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: l_min_
            !! smallest value of \\( l = |j-J|\\)
         integer(int32), intent(in) :: l_max_
            !! largest value of \\( l = j+J\\)
         integer(int32), intent(in) :: level_index_
            !! index pointing to speceific \\(j\\) value in rot_levels
         integer(int32), intent(inout) :: number_of_channels_even_parity_block
            !! number of channels for the p = 1 block
         integer(int32), intent(inout) :: number_of_channels_odd_parity_block
            !! number of channels for the p = -1 block
         !---------------------------------------------------------------------!
         integer(int32) :: l_
         !---------------------------------------------------------------------!
         do l_ = l_min_, l_max_
            if (mod(rot_levels(level_index_) + l_, 2) == 0) then
               number_of_channels_even_parity_block =                          &
                  number_of_channels_even_parity_block + 1
            else
               number_of_channels_odd_parity_block =                           &
                  number_of_channels_odd_parity_block + 1
            endif
         end do
         !---------------------------------------------------------------------!
      end subroutine update_channel_counts_space_fixed
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine check_number_of_channels(number_of_channels_bf,               &
         number_of_channels_sf, parity_block)
         !! check if the number of channels is the same in body-fixed
         !! and space-fixed frames
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels_bf
            !! number of channels in the body-fixed frame
         integer(int32), intent(in) :: number_of_channels_sf
            !! number of channels in the space-fixed frame
         character(len=*), intent(in) :: parity_block
            !! "even" or "odd", for printing purposes
         !---------------------------------------------------------------------!
         if (number_of_channels_bf /= number_of_channels_sf) then
            call write_error("Different number of channels in "// parity_block &
            // " block (BF = " // trim(adjustl(integer_to_character(           &
            number_of_channels_bf))) // ", SF = " // trim(adjustl(             &
            integer_to_character(number_of_channels_sf))) // "); check "       &
            // "set_number_of_channels")
         endif
      end subroutine check_number_of_channels
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine set_body_fixed_channels(total_angular_momentum_,              &
         parity_exponent_, channel_indices, channels_omega_values)
         !! Prepares the channel_indices array which holds indices that refer
         !! to the basis arrays: v1level/j1level/internal_energies, and
         !! channels_omega_values which holds values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_, channel_index_, omega_max_,           &
            parity_term_exponent_
         !---------------------------------------------------------------------!
         ! due to construction of body-fixed basis states:
         ! |v j \bar{\Omega} J p > = N (|v j \bar{\Omega} J >
         !                         + p (-1)^{J} |v j -\bar{\Omega} J > )
         ! we are interested in the exponent of the "p (-1)^{J}" term
         !---------------------------------------------------------------------!
         parity_term_exponent_=mod(parity_exponent_ + total_angular_momentum_,2)
         !---------------------------------------------------------------------!
         channel_index_ = 0
         !---------------------------------------------------------------------!
         do level_index_ = 1, number_of_basis_levels
            omega_max_ = min(rot_levels(level_index_), total_angular_momentum_)
            call update_body_fixed_channels_info(omega_max_,                   &
               parity_term_exponent_, level_index_, channel_index_,            &
               channel_indices, channels_omega_values)
          enddo
         !---------------------------------------------------------------------!
      end subroutine set_body_fixed_channels
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine update_body_fixed_channels_info(omega_max_,                   &
         parity_term_exponent_, level_index_, channel_index_, channel_indices, &
         channels_omega_values)
         !! update channel_indices array which holds indices within the
         !! loop over level_index_ in set_body_fixed_channels
         !---------------------------------------------------------------------!
         integer(int32) :: omega_max_
            !! largest value of \\(\bar{\Omega}\\), for given
            !! rotational and total angular momenta
         integer(int32) :: parity_term_exponent_
            !! exponent of the \\(p (-1)^{J}\\) term
         integer(int32) :: level_index_
            !! indices pointing to the basis arrays
         integer(int32), intent(inout) :: channel_index_
            !! index pointing to the current value in channel_indices
            !! and channels_omega_values; incremented in this subroutine
         integer(int32), intent(inout) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: omega_, omega_start_
         !---------------------------------------------------------------------!
         ! if p (-1)^{J} = 1, \bar{\Omega} states enter the basis
         ! otherwise, \bar{\Omega} > 0; to avoid redundancy, we handle this
         ! with omega_start_ which is 0 if parity_term_exponent_ is 0,
         ! otherwise 1
         !---------------------------------------------------------------------!
         omega_start_ = parity_term_exponent_
         !---------------------------------------------------------------------!
         do omega_ = omega_start_, omega_max_
            channel_index_ = channel_index_ + 1
            if (channel_index_ > size(channel_indices)) then
               call write_error("channel_index_ out of bounds of " //          &
                  "channel_indices in set_body_fixed_channels.")
            end if
            channels_omega_values(channel_index_)  = omega_
            channel_indices(channel_index_) = level_index_
         enddo
         !---------------------------------------------------------------------!
      end subroutine update_body_fixed_channels_info
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine set_space_fixed_channels(total_angular_momentum_,             &
         parity_exponent_, channel_l_values)
         !! Prepares the channel_l_values array which holds values of orbital
         !! angular momentum, \\(l\\), a space-fixed-frame quantum number.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channel_l_values(:)
            !! holds all values of l
         !---------------------------------------------------------------------!
         integer :: level_index_, l_min_, l_max_, l_, channel_index_
         !---------------------------------------------------------------------!
         channel_index_ = 0

         do level_index_ = 1, number_of_basis_levels
            l_min_ = abs(total_angular_momentum_ - rot_levels(level_index_))
            l_max_ = total_angular_momentum_ + rot_levels(level_index_)

            do l_ = l_min_, l_max_
               if (mod(l_ + rot_levels(level_index_), 2)==parity_exponent_) then
                  channel_index_ = channel_index_ + 1
                  if (channel_index_ > size(channel_l_values)) then
                     call write_error("channel_index_ out of bounds of " //    &
                        "channel_l_values in set_space_fixed_channels.")
                  end if
                  channel_l_values(channel_index_) = l_
               endif
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine set_space_fixed_channels
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function count_open_channels_in_block(channel_indices)                   &
         result(number_of_open_channels_)
         !! counts the energetically accessible channels in the given block
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32) :: number_of_open_channels_
            !! (output) number of open channels
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         !---------------------------------------------------------------------!
         number_of_open_channels_ = 0
         do channel_index_ = 1, size(channel_indices)
            if (is_open(internal_energies(channel_indices(channel_index_)))) then
               number_of_open_channels_ = number_of_open_channels_ + 1
            endif
         enddo
         !---------------------------------------------------------------------!
      end function count_open_channels_in_block
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine prepare_wavevector_array(channel_indices_, wavevectors_)
         !! Prepare an array of wavevectors in a given block (in A^2)
         !! which are saved in the S-matrix file
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         real(dp), intent(inout) :: wavevectors_(:)
            !! array of wavevectors (in A^2)
         !---------------------------------------------------------------------!
         integer(int32) :: index_, wv_index_
         real(dp) :: internal_energy_
         !---------------------------------------------------------------------!
         wv_index_ = 0
         do index_ = 1, size(channel_indices_)
            internal_energy_ = internal_energies(channel_indices_(index_))
            if (is_open(internal_energy_)) then
               wv_index_ = wv_index_ + 1
               wavevectors_ (wv_index_) =                                      &
                  sqrt( wavevector_squared_from_energy(internal_energy_) )     &
                  / bohr_to_angstrom
            endif
         enddo
         !---------------------------------------------------------------------!
         if (wv_index_ /= size(wavevectors_)) then
            call write_error("prepare_wavevector_array: mismatch between " //  &
               "number of open channels and size of wavevectors array")
         endif
         !---------------------------------------------------------------------!
      end subroutine prepare_wavevector_array
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function calculate_largest_wavevector(channel_indices)                   &
         result(largest_wavevector_)
         !! Calculates the largest wavevector in the block;
         !! called only if there are any open channels
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         real(dp) :: largest_wavevector_
            !! (output) the largest wavevector in the block
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         real(dp) :: wavevector_
         !---------------------------------------------------------------------!
         wavevector_ = 0.0_dp
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(channel_indices)
            if (is_open(internal_energies(channel_indices(channel_index_)))) then
               wavevector_ = sqrt( wavevector_squared_from_energy(             &
                  internal_energies(channel_indices(channel_index_))) )
               largest_wavevector_ = max(largest_wavevector_, wavevector_)
            endif
         enddo
         !---------------------------------------------------------------------!
      end function calculate_largest_wavevector
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function calculate_number_of_steps(largest_wavevector_,                  &
         k_potential_depth_) result(number_of_steps_)
         !! Calculates the number of steps on the intermolecular (R) grid in
         !! current block. This is done either directly (if r_step > 0)
         !! or through the number of steps per half de Broglie wavelength
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: largest_wavevector_
            !! the largest wavevector in the block
         real(dp), intent(in) :: k_potential_depth_
            !! correction from the depth of the potential well converted
            !! to wavevector units
         integer(int32) :: number_of_steps_
            !! number of steps on the \\(R\\) grid
         !---------------------------------------------------------------------!
         if (r_step <= 0) then
            number_of_steps_ = nint( (r_max-r_min) / PI                        &
               * ((largest_wavevector_+k_potential_depth_)*steps))
         else
            number_of_steps_ = nint((r_max-r_min)/r_step)+1
         endif
         !---------------------------------------------------------------------!
      end function calculate_number_of_steps
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_short_block_summary(jtot_, parity_exponent_,            &
         count_blocks_, number_of_channels_)
         integer(int32), intent(in) :: jtot_
            !! total angular momentum
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(in) :: count_blocks_
            !! ...
         integer(int32), intent(in) :: number_of_channels_
            !! number of channels in the block
         !---------------------------------------------------------------------!
         call write_message(" - Block number: " //                             &
            integer_to_character(count_blocks_))
         call write_message(" - Total angular momentum: " //                   &
            trim(adjustl(integer_to_character(jtot_))))
         call write_message(" - Parity: " //                                   &
            trim(adjustl(integer_to_character((-1)**parity_exponent_) )))
         call write_message(" - Number of scattering channels: " //            &
            integer_to_character(number_of_channels_))
         !---------------------------------------------------------------------!
      end subroutine
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_channels(parity_exponent_, channel_indices,             &
         channels_omega_values)
         !! prints information about body-fixed channels on screen
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channel_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_, v_, j_, omega_, parity_
         real(dp) :: internal_energy_, wavevector_
         !---------------------------------------------------------------------!
         call write_message("   v1      j1     omega      p" // repeat(" ", 10)&
            // "E_vj" // repeat(" ", 16) // "wv")
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(channel_indices)
            v_               = vib_levels(channel_indices(channel_index_))
            j_               = rot_levels(channel_indices(channel_index_))
            omega_           = channels_omega_values(channel_index_)
            parity_          = (-1)**parity_exponent_
            internal_energy_ = internal_energies(channel_indices(channel_index_))
            !------------------------------------------------------------------!
            ! format for open channels:
            !------------------------------------------------------------------!
            if (is_open(internal_energy_)) then
               wavevector_ = sqrt( wavevector_squared_from_energy(             &
                  internal_energy_) )
               call write_channel_line(v_, j_, omega_, parity_,                &
                  internal_energy_, wavevector_)
            !------------------------------------------------------------------!
            ! format for closed channels:
            !------------------------------------------------------------------!
            else
               call write_channel_line(v_,j_,omega_,parity_,internal_energy_)
            endif
            !------------------------------------------------------------------!
         enddo
         !---------------------------------------------------------------------!
      end subroutine print_channels
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine write_channel_line(v_, j_, omega_, parity_, internal_energy_, & 
         wavevector_)
         ! Subroutine arguments
         integer(int32), intent(in) :: v_
            !! vibrational quantum number
         integer(int32), intent(in) :: j_
            !! rotational quantum number
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         integer(int32), intent(in) :: parity_
            !! parity of the block
         real(dp), intent(in) :: internal_energy_
            !! \\(E_{vj}\\)
         real(dp), intent(in), optional :: wavevector_
            !! (optional) if the channel is open, print information
            !! about the wavevector
         !---------------------------------------------------------------------!
         character(len=200) :: line_
         !---------------------------------------------------------------------!
         ! Check if wavevector is provided
         !---------------------------------------------------------------------!
         if (present(wavevector_)) then
            write(line_, "(X,I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,F14.8)")         &
               v_, j_, omega_, parity_, internal_energy_*hartree_to_cm,        &
               wavevector_/bohr_to_angstrom
         else
            write(line_, "(X,I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,'--------------')")&
               v_, j_, omega_, parity_, internal_energy_*hartree_to_cm
         endif
         !---------------------------------------------------------------------!
         ! Print the formatted line
         !---------------------------------------------------------------------!
         call write_message(line_)
         !---------------------------------------------------------------------!
      end subroutine write_channel_line
   !---------------------------------------------------------------------------!
end module channels_mod
