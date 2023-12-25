module channels_mod
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use utility_functions_mod, only: write_error, write_message, write_warning,     &
      integer_to_character, float_to_character
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      subroutine set_number_of_channels(jtot_, size_BF_even, size_BF_odd)
         !! determine the number of scattering channels in each parity block 
         !! for given JTOT; check both BF and SF frames
         !---------------------------------------------------------------------!
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: jtot_
            !! total angular momentum         
         integer(int32), intent(inout) :: size_BF_even
            !! number of channels in p = 1 block
         integer(int32), intent(inout) :: size_BF_odd
            !! number of channels in p = -1 block
         !---------------------------------------------------------------------!
         integer(int32) :: ii, j1tmp, omegamax, lmin, lmax, l_len, len_even,   &
            len_odd, size_SF_even, size_SF_odd
         !---------------------------------------------------------------------!
         ! BF frame
         !---------------------------------------------------------------------!
         size_BF_even = 0
         size_BF_odd = 0
         do ii = 1, nlevel
            j1tmp    = j1array(ii)
            omegamax = min(j1tmp, jtot_)
            size_BF_even = size_BF_even+omegamax+1
            size_BF_odd  = size_BF_odd+omegamax
         enddo
         !---------------------------------------------------------------------!
         ! SF frame
         !---------------------------------------------------------------------!
         size_SF_even = 0
         size_SF_odd = 0
         do ii = 1,nlevel
            j1tmp = j1array(ii)
            lmin = abs(jtot_-j1tmp)
            lmax = jtot_+j1tmp

            l_len = lmax-lmin+1
            len_even = int((l_len+1)/2)
            len_odd  = int((l_len-1)/2)
            size_SF_even = size_SF_even+len_even
            size_SF_odd  = size_SF_odd+len_odd
         enddo
         !---------------------------------------------------------------------!
         ! Check if the results are the same
         !---------------------------------------------------------------------!
         if (size_BF_even.ne.size_SF_even) then
            call write_error("Different number of channels in even block " //  &
               "(BF/SF); check set_number_of_channels")
         endif
         !---------------------------------------------------------------------!
         if (size_BF_odd.ne.size_SF_odd) then
            call write_error("Different number of channels in odd block " //   &
               "(BF/SF); check set_number_of_channels")
         endif
         !---------------------------------------------------------------------!
      end subroutine set_number_of_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine set_channels_level_indices_channels_omega_values(number_of_channels, jtot_, nparity, channels_level_indices,   &
         channels_omega_values)
         !! Prepares the channels_level_indices array which holds indices that refer to the
         !! basis arrays: v1level/j1level/elevel, and channels_omega_values which holds values
         !! of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! number of channels in the block
         integer(int32), intent(in) :: jtot_
            !! total angular momentum
         integer(int32), intent(in) :: nparity
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         integer(int32) :: isize_, omegatmp, omegamax, ij, iomega
         !---------------------------------------------------------------------!
         isize_ = 0
         !---------------------------------------------------------------------!
         do ij = 1,nlevel
            omegamax = min(j1array(ij), jtot_)
            if (nparity.eq.0) then
               do iomega = 1, omegamax + 1
                  if (isize_ >= number_of_channels) then
                     call write_error("Array index isize_ out of bounds of " //      &
                        "channels_l_values in set_channels_level_indices_channels_omega_values.")
                  end if
                  omegatmp = iomega - 1
                  channels_omega_values(isize_ + 1) = omegatmp
                  channels_level_indices(isize_ + 1) = ij
                  isize_ = isize_ + 1 
               enddo
            else
               do iomega = 1, omegamax
                  if (isize_ >= number_of_channels) then
                     call write_error("Array index isize_ out of bounds of " //      &
                        "channels_l_values in set_channels_level_indices_channels_omega_values.")
                  end if
                  omegatmp = iomega
                  channels_omega_values(isize_ + 1) = omegatmp
                  channels_level_indices(isize_ + 1) = ij
                  isize_ = isize_ + 1
               enddo
            endif
          enddo
         !---------------------------------------------------------------------!
      end subroutine set_channels_level_indices_channels_omega_values
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine set_channels_l_values(number_of_channels, jtot_, nparity, channels_l_values)
         !! Prepares the channels_l_values array which holds values of orbital angular momentum
         !! (l), an SF-frame quantum number.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! number of channels in the block
         integer(int32), intent(in) :: jtot_
            !! total angular momentum
         integer(int32), intent(in) :: nparity
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_l_values(number_of_channels)
            !! holds all values of l
         !---------------------------------------------------------------------!
         integer :: lmin, lmax, ltmp
         integer :: ij, il, isize_
         !---------------------------------------------------------------------!
         isize_ = 0

         do ij = 1, nlevel
            lmin = abs(jtot_ - j1array(ij))
            lmax = jtot_ + j1array(ij)

            do il = lmin, lmax
               ltmp = il
               
               if (mod(ltmp + j1array(ij) + jtot_, 2) == nparity) then
                  if (isize_ >= number_of_channels) then
                     call write_error("Array index isize_ out of bounds of " //      &
                        "channels_l_values in set_channels_l_values.")
                  end if
                  channels_l_values(isize_ + 1) = ltmp
                  isize_ = isize_ + 1
               endif
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine set_channels_l_values
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine set_number_of_open_channels_wavmax(number_of_channels,        &
         channels_level_indices, channels_omega_values,                        &
         number_of_open_channels, wavmax)
         !! Calculates the number of energetically open channels (number_of_open_channels) and
         !! the largest wave number (wavmax) in the block.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! number of channels in the block
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(inout) :: number_of_open_channels
            !! number of open channels
         real(dp), intent(inout) :: wavmax
            !! the largest wave number (wavmax) in the block
         !---------------------------------------------------------------------!
         integer(int32) :: isize_
         real(dp) :: wavenumber
         !---------------------------------------------------------------------!
         number_of_open_channels = 0
         wavmax = 0.0_dp
         !---------------------------------------------------------------------!
         do isize_ = 1, number_of_channels
            if ((ETOTAL() - elevel(channels_level_indices(isize_))) > 0.0_dp) then
               number_of_open_channels = number_of_open_channels + 1
               wavenumber = dsqrt(2 * reducedmass * (ETOTAL() - elevel(channels_level_indices(isize_))))
               wavmax = max(wavmax, wavenumber)
            endif
         enddo
         !---------------------------------------------------------------------!
      end subroutine set_number_of_open_channels_wavmax
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine print_bf_channels(number_of_channels, jtot_, nparity, channels_level_indices, channels_omega_values)
         !! prints BF quantum numbers on screen
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! number of channels in the block
         integer(int32), intent(in) :: jtot_
            !! total angular momentum
         integer(int32), intent(in) :: nparity
            !! parity exponent of the block (0 if p = +1, 1 if p = -1)
         integer(int32), intent(inout) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(inout) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         !---------------------------------------------------------------------!
         character(len = 200) :: line_
         integer(int32) :: isize_, itmp, v1tmp, j1tmp, omegatmp
         real(dp) :: erot, wavenumber
         !---------------------------------------------------------------------!
         call write_message("  v1      j1     omega      p" // repeat(" ", 11) //    &
            "E_j" // repeat(" ", 16) // "wv")
         !---------------------------------------------------------------------!
         do isize_ = 1, number_of_channels
            v1tmp    = v1array(channels_level_indices(isize_))
            j1tmp    = j1array(channels_level_indices(isize_))
            omegatmp = channels_omega_values(isize_)
            erot     = elevel(channels_level_indices(isize_))
            itmp     = (-1)**nparity
            if ((ETOTAL()-erot) <= 0.0_dp) then
               write(line_, "(I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,'--------------')") &
                  v1tmp, j1tmp, omegatmp, itmp*(-1)**jtot_, erot*hartreetocm
               call write_message(line_)
            else
               wavenumber = dsqrt(2 * reducedmass * (ETOTAL() - erot))
               write(line_, "(I4,4X,I4,6X,I4,5X,I2,2X,F12.4,4X,F14.8)")        &
                  v1tmp, j1tmp, omegatmp, itmp*(-1)**jtot_, erot*hartreetocm,     &
                  wavenumber/bohrtoangstrom
               call write_message(line_)
            endif
         enddo
         !---------------------------------------------------------------------!
      end subroutine print_bf_channels
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
end module channels_mod
