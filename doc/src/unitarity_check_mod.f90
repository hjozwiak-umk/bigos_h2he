module unitarity_check_mod
   !! This module contains subroutines that check the unitarity condition
   !! of the S-matrix (see Eq. (13) in "Solution of coupled equations")
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use data_mod, only: unitary_tolerance
   use utility_functions_mod, only: integer_to_character, float_to_character,  &
      write_warning, write_message
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      subroutine unitarity_check(number_of_open_channels, s_matrix_real,       &
         s_matrix_imag,is_unitary)
         !! checks the unitarity of the S-matrix
         !! (Eq. (13) in "Solution of coupled equations")
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_open_channels
            !! number of open channels
         real(dp), intent(in) :: s_matrix_real(number_of_open_channels,number_of_open_channels)
            !! real part of the S-matrix
         real(dp), intent(in) :: s_matrix_imag(number_of_open_channels,number_of_open_channels)
            !! imaginary part of the S-matrix
         logical, intent(inout) :: is_unitary
            !! (output) if .true. unitarity is fulfilled, .false. otherwise
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index
         real(dp) :: sum_of_squares(number_of_open_channels)
         !---------------------------------------------------------------------!
         is_unitary  = .true.
         !---------------------------------------------------------------------!
         call write_message("Check of the unitarity of the S-matrix...")
         !---------------------------------------------------------------------!
         ! Calculating sum of squares for each channel
         !---------------------------------------------------------------------!
         call calculate_sum_of_squares_for_each_channel(s_matrix_real,         &
            s_matrix_imag, sum_of_squares)
         !---------------------------------------------------------------------!
         ! Checking unitarity for each channel
         !---------------------------------------------------------------------!
         is_unitary = check_unitarity_for_each_channel(sum_of_squares)
         !---------------------------------------------------------------------!
         ! Handling the output message based on unitarity check
         !---------------------------------------------------------------------!
         call handle_unitarity_output_message(is_unitary, sum_of_squares)
         !---------------------------------------------------------------------!
      end subroutine unitarity_check
      !------------------------------------------------------------------------!
      subroutine calculate_sum_of_squares_for_each_channel(s_matrix_real,        &
         s_matrix_imag, sum_of_squares_)
         !! calculates the sum
         !! \\( \sum\_{\gamma'} \Bigl|{S}^{Jp}\_{\gamma, \gamma'}\Bigr|^{2} \\)
         !! for all \\(\gamma\\) channels
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: s_matrix_real(:,:)
            !! real part of the S-matrix
         real(dp), intent(in) :: s_matrix_imag(:,:)
            !! imaginary part of the S-matrix
         real(dp), intent(out) :: sum_of_squares_(:)
            !! (output) \\( \sum\_{\gamma'} \Bigl|{S}^{Jp}\_{\gamma, \gamma'}\Bigr|^{2} \\)
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(s_matrix_real, 1)
            sum_of_squares_(channel_index_) =                                  &
               sum(s_matrix_real(channel_index_, :)**2                         &
               + s_matrix_imag(channel_index_, :)**2)
         enddo
         !---------------------------------------------------------------------!
      end subroutine calculate_sum_of_squares_for_each_channel
      !------------------------------------------------------------------------!
      function check_unitarity_for_each_channel(sum_of_squares)                &
         result(is_unitary_)
         !! checks if the calculated sum of squares equals 1 for each channel
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: sum_of_squares(:)
            !! \\( \sum\_{\gamma'} \Bigl|{S}^{Jp}\_{\gamma, \gamma'}\Bigr|^{2} \\)
         logical :: is_unitary_
            !! (output) if .true. unitarity is fulfilled, .false. otherwise
         !---------------------------------------------------------------------!
         integer :: channel_index
         !---------------------------------------------------------------------!
         is_unitary_ = .true.
         do channel_index = 1, size(sum_of_squares)
            if (abs(sum_of_squares(channel_index) - 1.0_dp) > unitary_tolerance) then
               is_unitary_ = .false.
               exit
            endif
          end do
         !---------------------------------------------------------------------!
      end function check_unitarity_for_each_channel
      !------------------------------------------------------------------------!
      subroutine handle_unitarity_output_message(is_unitary, sum_of_squares)
         !! handle printing messages depending on the outcome of unitarity check
         !---------------------------------------------------------------------!
         logical, intent(in) :: is_unitary
            !! if .true. unitarity is fulfilled, .false. otherwise
         real(dp), intent(in) :: sum_of_squares(:)
            !! array holding
            !! \\(\sum\_{\gamma^{\prime}}|S\_{\gamma,\gamma^{\prime}}|^{2}\\)
            !! for each \\(\gamma\\)
         !---------------------------------------------------------------------!
         if (is_unitary) then
            call write_message("S-matrix unitary condition fulfilled")
         else
            call write_warning("Unitary condition is not fulfilled for one or more channels")
            call write_message("Consider increasing the STEPS parameter")
            call print_sum_of_squares(sum_of_squares)
         endif
         !---------------------------------------------------------------------!
      end subroutine handle_unitarity_output_message
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine print_sum_of_squares(sum_of_squares)
         !! print S-matrix on screen
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: sum_of_squares(:)
            !! array holding
            !! \\(\sum\_{\gamma^{\prime}}|S\_{\gamma,\gamma^{\prime}}|^{2}\\)
            !! for each \\(\gamma\\)
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_
         !---------------------------------------------------------------------!
         call write_message(repeat(" ", 3)// "row" // repeat(" ", 12)//        &
            "sum(S**2)")
         !---------------------------------------------------------------------!
         do channel_index_ = 1, size(sum_of_squares)
            call write_message(" " //                                          &
               trim(adjustl(integer_to_character(channel_index_, "(i5)"))) //  &
               repeat(" ", 8) //                                               &
               trim(adjustl(float_to_character(sum_of_squares(channel_index_), &
               "(E15.8)"))))
         enddo
         !---------------------------------------------------------------------!
      end subroutine print_sum_of_squares
   !---------------------------------------------------------------------------!
end module unitarity_check_mod
