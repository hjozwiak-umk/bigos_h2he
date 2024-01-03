module channels_mod
   !! This module provides subroutines that set the number of channels in the
   !! block, save quantum numbers for each channel (both in body- and space-fixed
   !! cases) and print quantum numbers on screen
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use utility_functions_mod, only: write_error, write_message, write_warning, &
      integer_to_character, float_to_character
   !---------------------------------------------------------------------------!
   implicit none
   private
   public :: set_number_of_channels, set_body_fixed_channels,                  &
      set_space_fixed_channels, count_open_channels_in_block,                  &
      calculate_largest_wavenumber, print_channels
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      subroutine set_number_of_channels(total_angular_momentum_,               &
         number_of_channels_even_parity_block, number_of_channels_odd_parity_block)
         !! determine the number of scattering channels in each parity block 
         !! for given total angular momentum in both body-fixed and
         !! space-fixed frames
         !---------------------------------------------------------------------!
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
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
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
         do level_index_ = 1, nlevel
            omega_max_ = min(j1array(level_index_), total_angular_momentum_)

            call update_channel_counts_body_fixed(omega_max_, total_angular_momentum_,   &
               number_of_channels_even_parity_block, number_of_channels_odd_parity_block)

         end do
         !---------------------------------------------------------------------!
      end subroutine calculate_number_of_channels_body_fixed
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
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
            number_of_channels_even_parity_block = number_of_channels_even_parity_block + omega_max_ + 1
            number_of_channels_odd_parity_block = number_of_channels_odd_parity_block + omega_max_
         else
            !------------------------------------------------------------------!
            ! Odd J: channels with Omega = 0 only count in the odd parity block
            !------------------------------------------------------------------!
            number_of_channels_odd_parity_block = number_of_channels_odd_parity_block + omega_max_ + 1
            number_of_channels_even_parity_block = number_of_channels_even_parity_block + omega_max_
         endif
         !---------------------------------------------------------------------!
      end subroutine update_channel_counts_body_fixed
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
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
         do level_index_ = 1, nlevel
            l_min_ = abs(total_angular_momentum_ - j1array(level_index_))
            l_max_ = total_angular_momentum_ + j1array(level_index_)

            call update_channel_counts_space_fixed(l_min_, l_max_, level_index_, &
               number_of_channels_even_parity_block, number_of_channels_odd_parity_block)
               
         end do
         !---------------------------------------------------------------------!
      end subroutine calculate_number_of_channels_space_fixed
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine update_channel_counts_space_fixed(l_min_, l_max_, level_index_, &
         number_of_channels_even_parity_block, number_of_channels_odd_parity_block)
         !! updates number_of_channels_even and number_of_channels_odd
         !! in the space-fixed frame for given
         !! range of orbital angular momentum
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: l_min_
            !! smallest value of \\( l = |j-J|\\)
         integer(int32), intent(in) :: l_max_
            !! largest value of \\( l = j+J\\)
         integer(int32), intent(in) :: level_index_
            !! index pointing to speceific \\(j\\) value in j1array
         integer(int32), intent(inout) :: number_of_channels_even_parity_block
            !! number of channels for the p = 1 block
         integer(int32), intent(inout) :: number_of_channels_odd_parity_block
            !! number of channels for the p = -1 block
         !---------------------------------------------------------------------!
         integer(int32) :: l_
         !---------------------------------------------------------------------!
         do l_ = l_min_, l_max_
            if (mod(j1array(level_index_) + l_, 2) == 0) then
               number_of_channels_even_parity_block =                       &
                  number_of_channels_even_parity_block + 1
            else
               number_of_channels_odd_parity_block =                        &
                  number_of_channels_odd_parity_block + 1
            endif
         end do
         !---------------------------------------------------------------------!
      end subroutine update_channel_counts_space_fixed
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
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
            call write_error("Different number of channels in " // parity_block &
            // " block (BF = " // trim(adjustl(integer_to_character(number_of_channels_bf))) &
            // ", SF = " // trim(adjustl(integer_to_character(number_of_channels_sf))) &
            // "); check set_number_of_channels")
         endif
      end subroutine check_number_of_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine set_body_fixed_channels(total_angular_momentum_,              &
         parity_exponent_, channels_level_indices, channels_omega_values)
         !! Prepares the channels_level_indices array which holds indices that refer to the
         !! basis arrays: v1level/j1level/elevel, and channels_omega_values which holds values
         !! of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: level_index_, channel_index_, omega_max_, parity_term_exponent_
         !---------------------------------------------------------------------!
         ! due to construction of body-fixed basis states:
         ! |v j \bar{\Omega} J p > = N (|v j \bar{\Omega} J >
         !                         + p (-1)^{J} |v j -\bar{\Omega} J > )
         ! we are interested in the exponent of the "p (-1)^{J}" term
         !---------------------------------------------------------------------!
         parity_term_exponent_ = mod(parity_exponent_ + total_angular_momentum_, 2)
         !---------------------------------------------------------------------!
         channel_index_ = 0
         !---------------------------------------------------------------------!
         do level_index_ = 1, nlevel
            omega_max_ = min(j1array(level_index_), total_angular_momentum_)
            call update_body_fixed_channels_info(omega_max_, parity_term_exponent_,     &
               level_index_, channel_index_, channels_level_indices, channels_omega_values)
          enddo
         !---------------------------------------------------------------------!
      end subroutine set_body_fixed_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine update_body_fixed_channels_info(omega_max_, parity_term_exponent_, &
         level_index_, channel_index_, channels_level_indices, channels_omega_values)
         !! update channels_level_indices array which holds indices within the
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
            !! index pointing to the current value in channels_level_indices
            !! and channels_omega_values; incremented in this subroutine
         integer(int32), intent(inout) :: channels_level_indices(:)
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
            if (channel_index_ > size(channels_level_indices)) then
               call write_error("channel_index_ out of bounds of " //      &
                  "channels_level_indices in set_body_fixed_channels.")
            end if
            channels_omega_values(channel_index_)  = omega_
            channels_level_indices(channel_index_) = level_index_
         enddo
         !---------------------------------------------------------------------!
      end subroutine update_body_fixed_channels_info
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine set_space_fixed_channels(total_angular_momentum_,             &
         parity_exponent_, channels_l_values)
         !! Prepares the channels_l_values array which holds values of
         !! orbital angular momentum, \\(l\\), a space-fixed-frame quantum number.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_l_values(:)
            !! holds all values of l
         !---------------------------------------------------------------------!
         integer :: level_index_, l_min_, l_max_, l_, channel_index_
         !---------------------------------------------------------------------!
         channel_index_ = 0

         do level_index_ = 1, nlevel
            l_min_ = abs(total_angular_momentum_ - j1array(level_index_))
            l_max_ = total_angular_momentum_ + j1array(level_index_)

            do l_ = l_min_, l_max_
               if (mod(l_ + j1array(level_index_), 2) == parity_exponent_) then
                  channel_index_ = channel_index_ + 1
                  if (channel_index_ > size(channels_l_values)) then
                     call write_error("channel_index_ out of bounds of " //    &
                        "channels_l_values in set_space_fixed_channels.")
                  end if
                  channels_l_values(channel_index_) = l_
               endif
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine set_space_fixed_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      function count_open_channels_in_block(channels_level_indices)            &
         result(number_of_open_channels_)
         !! counts the energetically accessible channels in the given block
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32) :: number_of_open_channels_
            !! (output) number of open channels
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         !---------------------------------------------------------------------!
         number_of_open_channels_ = 0
         do channel_index_ = 1, size(channels_level_indices)
            if (is_open(elevel(channels_level_indices(channel_index_)))) then
               number_of_open_channels_ = number_of_open_channels_ + 1
            endif
         enddo
         !---------------------------------------------------------------------!
      end function count_open_channels_in_block
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function calculate_largest_wavenumber(channels_level_indices) result(largest_wavenumber_)
         !! Calculates the largest wave number in the block;
         !! called only if there are any open channels
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         real(dp) :: largest_wavenumber_
            !! (output) the largest wave number (wavmax) in the block
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         real(dp) :: wavenumber_
         !---------------------------------------------------------------------!
         wavenumber_ = 0.0_dp
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(channels_level_indices)
            if (is_open(elevel(channels_level_indices(channel_index_)))) then
               wavenumber_ = wavenumber_from_energy(elevel(channels_level_indices(channel_index_)))
               largest_wavenumber_ = max(largest_wavenumber_, wavenumber_)
            endif
         enddo
         !---------------------------------------------------------------------!
      end function calculate_largest_wavenumber
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine print_channels(parity_exponent_, channels_level_indices,      &
         channels_omega_values)
         !! prints information about body-fixed channels on screen
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: parity_exponent_
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_level_indices(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(:)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_, v_, j_, omega_, parity_
         real(dp) :: internal_energy_, wavenumber_
         !---------------------------------------------------------------------!
         call write_message("  v1      j1     omega      p" // repeat(" ", 10) &
            // "E_vj" // repeat(" ", 16) // "wv")
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(channels_level_indices)
            v_               = v1array(channels_level_indices(channel_index_))
            j_               = j1array(channels_level_indices(channel_index_))
            omega_           = channels_omega_values(channel_index_)
            parity_          = (-1)**parity_exponent_
            internal_energy_ = elevel(channels_level_indices(channel_index_))
            !------------------------------------------------------------------!
            ! format for open channels:
            !------------------------------------------------------------------!
            if (is_open(internal_energy_)) then
               wavenumber_ = wavenumber_from_energy(internal_energy_)
               call write_channel_line(v_, j_, omega_, parity_,                &
                  internal_energy_, wavenumber_)
            !------------------------------------------------------------------!
            ! format for closed channels:
            !------------------------------------------------------------------!
            else
               call write_channel_line(v_, j_, omega_, parity_, internal_energy_)
            endif
            !------------------------------------------------------------------!
         enddo
         !---------------------------------------------------------------------!
      end subroutine print_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine write_channel_line(v_, j_, omega_, parity_, internal_energy_, & 
         wavenumber_)
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
         real(dp), intent(in), optional :: wavenumber_
            !! (optional) if the channel is open, print information
            !! about the wavenumber
         !---------------------------------------------------------------------!
         character(len=200) :: line_
         !---------------------------------------------------------------------!
         ! Check if wavenumber is provided
         !---------------------------------------------------------------------!
         if (present(wavenumber_)) then
            write(line_, "(I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,F14.8)")           &
               v_, j_, omega_, parity_, internal_energy_*hartreetocm,          &
               wavenumber_/bohrtoangstrom
         else
            write(line_, "(I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,'--------------')")&
               v_, j_, omega_, parity_, internal_energy_*hartreetocm
         endif
         !---------------------------------------------------------------------!
         ! Print the formatted line
         !---------------------------------------------------------------------!
         call write_message(line_)
         !---------------------------------------------------------------------!
      end subroutine write_channel_line
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
end module channels_mod
