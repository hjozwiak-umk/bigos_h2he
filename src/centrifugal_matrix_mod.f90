module centrifugal_matrix_mod
   !! This module calculates the centrifugal matrix - see the second term
   !! in Eq. 3 in "What are coupled equations?" section.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use array_operations_mod, only: invert_symmetric_matrix, fill_symmetric_matrix
   use global_variables_mod
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: calculate_centrifugal_matrix
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                           Centrifugal matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_centrifugal_matrix(total_angular_momentum_,         &
         channel_indices_, channels_omega_values_, centrifugal_matrix_)    
         !! calculates the (R**2)*centrifugal matrix from the second term
         !! of Eq. 3 in "What are coupled equations?" section;
         !! Matrix elements are given in Eq. 4 and 6 of "Coupling Matrix" secion
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values_(:)
            !! holds all values of \bar{\Omega}
         real(dp), intent(out) :: centrifugal_matrix_(:,:)
            !! (output) - (R**2)*centrifugal matrix
         !---------------------------------------------------------------------!
         integer(int32) :: omega_, omega_prime_, v_, j_, v_prime_, j_prime_,  &
            channel_index_1_, channel_index_2_
         real(dp) :: centtmp, delta_1_, delta_2_
         !---------------------------------------------------------------------!
         centrifugal_matrix_  = 0

         do channel_index_1_ = 1, size(channel_indices_)
            v_ = vib_levels(channel_indices_(channel_index_1_))
            j_ = rot_levels(channel_indices_(channel_index_1_))
            omega_ = channels_omega_values_(channel_index_1_)
            delta_1_ = delta_for_zero_omega(omega_)
            do channel_index_2_ = 1, channel_index_1_
               v_prime_ = vib_levels(channel_indices_(channel_index_2_))
               j_prime_ = rot_levels(channel_indices_(channel_index_2_))
               omega_prime_ = channels_omega_values_(channel_index_2_)
               delta_2_ = delta_for_zero_omega(omega_prime_)
               !---------------------------------------------------------------!
               if (v_ /= v_prime_ .or. j_ /= j_prime_ .or.                     &
                  abs(omega_ - omega_prime_) > 1) then
                  cycle
               endif
               !---------------------------------------------------------------!
               if (omega_ == omega_prime_) then
                  !------------------------------------------------------------!
                  ! Eq. 4 in "Coupling Matrix" section
                  !------------------------------------------------------------!
                  centrifugal_matrix_(channel_index_1_, channel_index_2_)      &
                     = calculate_diagonal_centrifugal_element(                 &
                     total_angular_momentum_, j_, omega_)
               else
                  !------------------------------------------------------------!
                  ! Eq. 5 in "Coupling Matrix" section
                  !------------------------------------------------------------!
                  centrifugal_matrix_(channel_index_1_, channel_index_2_)      &
                     = calculate_offdiagonal_centrifugal_element(              &
                     total_angular_momentum_, j_, omega_, omega_prime_,        &
                     delta_1_, delta_2_)
               endif
               !---------------------------------------------------------------!
            enddo
         enddo
         !---------------------------------------------------------------------!
         call fill_symmetric_matrix(centrifugal_matrix_, 'u')
         !---------------------------------------------------------------------!
      end subroutine calculate_centrifugal_matrix
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function delta_for_zero_omega(omega_) result(delta_)
         !! Checks if the input value equals 0; used in the calculation
         !! of off-diagonal elements of the centrifugal matrix; see
         !! Eq. 5 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: omega_
            !! input value which is to be compared with 0
         real(dp) :: delta_
            !! (output) 1 if omega_ = 0, 0 otherwise
         !---------------------------------------------------------------------!
         if (omega_ == 0) then
           delta_ = 1.0_dp
         else
           delta_ = 0.0_dp
         endif
         !---------------------------------------------------------------------!
      end function delta_for_zero_omega
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function calculate_diagonal_centrifugal_element(total_angular_momentum_, &
         j_, omega_) result(diagonal_element_)
         !! calculates diagonal element of the centrifgual matrix, see
         !! Eq. 4 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         real(dp) :: diagonal_element_
            !! (output) value of the diagonal element of the centrifgual matrix
         !---------------------------------------------------------------------!
         diagonal_element_ = real(                                             &
            total_angular_momentum_ * (total_angular_momentum_ + 1)            &
            + j_ * (j_ + 1) - 2 * omega_ **2, dp)
         !---------------------------------------------------------------------!
      end function calculate_diagonal_centrifugal_element
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function calculate_offdiagonal_centrifugal_element(                      &
         total_angular_momentum_, j_, omega_, omega_prime_, delta_1_, delta_2_)&
         result(offdiagonal_element_)
         !! calculates off-diagonal element of the centrifgual matrix, see
         !! Eq. 5 in "Coupling Matrix" section
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational angular momentum
         integer(int32), intent(in) :: omega_
            !! \\(\bar{\Omega}\\)
         integer(int32), intent(in) :: omega_prime_
            !! \\(\bar{\Omega}'\\)
         real(dp), intent(in) :: delta_1_, delta_2_
            !! Kronecker delta functions determined earlier
         real(dp) :: offdiagonal_element_
            !! (output) value of the off-diagonal element of
            !! the centrifgual matrix
         !---------------------------------------------------------------------!
         offdiagonal_element_ = - sqrt(real(                                   &
            (total_angular_momentum_ * (total_angular_momentum_ + 1)           &
            - omega_ * omega_prime_) * (j_ * (j_ + 1) - omega_ * omega_prime_),&
            dp)) * sqrt( (1.0_dp + delta_1_) * (1.0_dp + delta_2_) )
         !---------------------------------------------------------------------!
      end function calculate_offdiagonal_centrifugal_element
!------------------------------------------------------------------------------!
end module centrifugal_matrix_mod
