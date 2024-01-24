module save_s_matrix_mod
   !! This module provides procedures that save selective information
   !! to the S-matrix file
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use global_variables_mod
   use utility_functions_mod, only: file_io_status
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: save_s_matrix_file_header, save_s_matrix_block_info
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine save_s_matrix_file_header
         !! save "header" of the S-matrix file:
         !! -- label, "itype", number of levels in the basis, reduced mass
         !!    of the system
         !! -- vibrational and rotational quantum numbers
         !! -- rovibrational energies
         !! -- index pointing to the initial level and the kinetic/total energy
         !---------------------------------------------------------------------!
         character(len=200) :: err_message
         integer(int32) :: io_status, level_index_
         !---------------------------------------------------------------------!
         open(s_matrix_unit, file=trim(s_matrix_file_name), form='unformatted',&
            iostat = io_status, iomsg = err_message)
         call file_io_status(io_status, err_message, s_matrix_unit, "o")
         !---------------------------------------------------------------------!
         write(s_matrix_unit) label, 2, number_of_basis_levels, reduced_mass
         write(s_matrix_unit) (vib_levels(level_index_),                       &
            rot_levels(level_index_), level_index_ = 1, number_of_basis_levels)
         write(s_matrix_unit) (internal_energies(level_index_),                &
            level_index_ = 1, number_of_basis_levels)
         write(s_matrix_unit) initial_level, energy
         !---------------------------------------------------------------------!
      end subroutine save_s_matrix_file_header
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine save_s_matrix_block_info(total_angular_momentum,              &
         parity_exponent, number_of_open_channels, channel_indices,            &
         channel_l_values, block_wavevectors, s_matrix_real, s_matrix_imag)
         !! save information about current block
         !! -- total angular momentum, parity exponent, number of open channels
         !!    in the current block
         !! -- array of indices pointing to the basis arrays, array holding
         !!    \\(l\\) values, wavenumbers
         !! -- real part of the S-matrix
         !! -- imaginary part of the S-matrix
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum
            !! total angular momentum of the current block
         integer(int32), intent(in) :: parity_exponent
            !! parity exponent of the current block
         integer(int32), intent(in) :: number_of_open_channels
            !! number of open channels in the block
         integer(int32), intent(in) :: channel_indices(number_of_open_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channel_l_values(number_of_open_channels)
            !! holds all values of \\(l\\)
         real(dp), intent(in) :: block_wavevectors(number_of_open_channels)
            !! holds all values of wavevectors in the block
         real(dp), intent(in) :: s_matrix_real(number_of_open_channels,        &
            number_of_open_channels)
            !! real part of the S-matrix
         real(dp), intent(in) :: s_matrix_imag(number_of_open_channels,        &
            number_of_open_channels)
            !! imaginary part of the S-matrix
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_, channel_index_2_
         !---------------------------------------------------------------------!
         write(s_matrix_unit) total_angular_momentum, parity_exponent,         &
            number_of_open_channels
         write(s_matrix_unit) (channel_indices(channel_index_),                &
            channel_l_values(channel_index_),block_wavevectors(channel_index_),&
            channel_index_ = 1, number_of_open_channels)
         write(s_matrix_unit)((s_matrix_real(channel_index_, channel_index_2_),&
            channel_index_2_ = 1, channel_index_), channel_index_=1,           &
            number_of_open_channels)
         write(s_matrix_unit) ((s_matrix_imag(channel_index_,channel_index_2_),&
            channel_index_2_ = 1, channel_index_), channel_index_=1,           &
            number_of_open_channels)
         !---------------------------------------------------------------------!
      end subroutine save_s_matrix_block_info
   !---------------------------------------------------------------------------!
end module save_s_matrix_mod
