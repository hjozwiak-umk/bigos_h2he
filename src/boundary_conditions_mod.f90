module boundary_conditions_mod
   !! This module contains subroutines that transform the asymptotic
   !! log-derivative matrix into the scattering S-matrix
   !! (see "Solution of coupled equations" section)
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use fwigxjpf, only: fwig3jj
   use math_functions_mod, only: riccati_bessel_j, riccati_bessel_y,           &
      modified_bessel_k_ratio
   use utility_functions_mod, only: write_warning, write_header
   use array_operations_mod, only: invert_symmetric_matrix, fill_symmetric_matrix
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: calculate_sf_matrix_from_bf_matrix, calculate_k_matrix,           &
      calculate_s_matrix
   !---------------------------------------------------------------------------!
   contains
   !---------------------------------------------------------------------------!
      function p_coeff(total_angular_momentum_,j_,l_,omega_) result(p_coeff_)
         !! calculates the P coefficients from Eq. (3) in
         !! "Solution of coupled equations" 
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational quantum number
         integer(int32), intent(in) :: l_
            !! orbital angular momentum
         integer(int32), intent(in) :: omega_
            !! projection of j on the BF-Z axis
         real(dp) :: p_coeff_
            !! result - P function (Eq. (3) in "Solution of coupled equations")
         !---------------------------------------------------------------------!
         real(dp) :: delta_
         !---------------------------------------------------------------------!
         delta_ = 0.d0
         if (omega_ == 0) delta_ = 1.0_dp
         p_coeff_ = (-1.0_dp)**(total_angular_momentum_+omega_) * dsqrt(2.0_dp)&
            * dsqrt(real(2*l_+1, dp))                                          &
            * fwig3jj(2* j_, 2* total_angular_momentum_, 2* l_,                &
            2 * omega_, -2* omega_, 0) / dsqrt(1.0_dp + delta_)
         !---------------------------------------------------------------------!
      end function p_coeff
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine calculate_single_SF_element(number_of_channels,               &
         total_angular_momentum_, v_, j_, vp_, jp_, l_, lp_,                   &
         channel_indices, channels_omega_values, bf_matrix, sf_element)
         !! calculates a single space-fixed matrix element from Eq. (2) in
         !! "Solution of coupled equations" 
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: v_, j_, vp_, jp_
            !! vibrational and rotational quantum numbers
         integer(int32), intent(in) :: l_, lp_
            !! orbital angular momenta
         integer(int32), intent(in) :: channel_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         real(dp), intent(in) :: bf_matrix(number_of_channels, number_of_channels)
            !! matrix in the BF frame
         real(dp), intent(out) :: sf_element
            !! (output) matrix element in the SF frame
         !---------------------------------------------------------------------!
         integer(int32) :: channel_index_1_, channel_index_2_
         real(dp) :: p_coeff_outer, p_coeff_inner, sum_outer, sum_inner
         !---------------------------------------------------------------------!
         sum_outer = 0.0_dp
         do channel_index_1_ = 1, number_of_channels
            if (v1array(channel_indices(channel_index_1_)) /= v_ .or.                 &
               j1array(channel_indices(channel_index_1_)) /= j_) cycle
            p_coeff_outer = p_coeff(total_angular_momentum_, j_, l_,           &
               channels_omega_values(channel_index_1_))
            sum_inner = 0.0_dp
            do channel_index_2_ = 1, number_of_channels
               if (v1array(channel_indices(channel_index_2_)) /= vp_ .or.             &
                  j1array(channel_indices(channel_index_2_)) /= jp_) cycle
               p_coeff_inner = p_coeff(total_angular_momentum_, jp_, lp_,      &
                  channels_omega_values(channel_index_2_))
               sum_inner = sum_inner + p_coeff_inner * bf_matrix(channel_index_1_,channel_index_2_)
            end do
            sum_outer = sum_outer + p_coeff_outer * sum_inner
         end do
         sf_element = sum_outer
         !---------------------------------------------------------------------!
      end subroutine calculate_single_SF_element
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine calculate_sf_matrix_from_bf_matrix(number_of_channels,        &
         total_angular_momentum_, channel_indices,                      &
         channels_omega_values,channel_l_values, bf_matrix, sf_matrix)
         !! takes as an input matrix in the body-fixed frame and transforms it 
         !! to the spec-fixed frame; iterates over all matrix elements
         !! and calls calculate_single_SF_element
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: channel_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: channel_l_values(number_of_channels)
            !! holds all values of l
         real(dp), intent(in) :: bf_matrix(number_of_channels, number_of_channels)
            !! matrix in the BF frame
         real(dp), intent(inout) :: sf_matrix(number_of_channels, number_of_channels)
            !! (output) matrix in the SF frame
         !---------------------------------------------------------------------!
         integer(int32) :: l_, lp_, omega_, omegap_, v1_, j1_, v1p_, j1p_, channel_index_1_, channel_index_2_
         real(dp) :: single_sf_element
         !---------------------------------------------------------------------!
         do channel_index_1_ = 1, number_of_channels
            v1_ = v1array(channel_indices(channel_index_1_))
            j1_ = j1array(channel_indices(channel_index_1_))
            l_  = channel_l_values(channel_index_1_)
            do channel_index_2_ = 1, number_of_channels
               v1p_ = v1array(channel_indices(channel_index_2_))
               j1p_ = j1array(channel_indices(channel_index_2_))
               lp_  = channel_l_values(channel_index_2_)
               call calculate_single_SF_element(number_of_channels,            &
                  total_angular_momentum_, v1_, j1_, v1p_, j1p_, l_, lp_,      &
                  channel_indices, channels_omega_values, bf_matrix,    &
                  single_sf_element)
               sf_matrix(channel_index_1_, channel_index_2_) = single_sf_element
            enddo
         enddo
         !---------------------------------------------------------------------!         
      end subroutine calculate_sf_matrix_from_bf_matrix
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine calculate_k_matrix(number_of_channels, log_der_matrix,        &
         number_of_open_channels, channel_indices, channel_l_values,   &
         r_, k_matrix)
         !! calculates the K-matrix from log-derivative matrix using Eq. (4) in
         !! "Solution of coupled equations" 
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! y-matrix is of number_of_channels x number_of_channels size
         real(dp), intent(in) :: log_der_matrix(number_of_channels, number_of_channels)
            !! asymptotic log-derivative matrix
         integer(int32), intent(in) :: number_of_open_channels
            !! number of open channels
         integer(int32), intent(in) :: channel_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channel_l_values(number_of_channels)
            !! holds all values of l
         real(dp), intent(in) :: r_
            !! Rmax
         real(dp), intent(inout) :: k_matrix(number_of_open_channels, number_of_open_channels)
            !! K-matrix
         !---------------------------------------------------------------------!
         integer(int32) :: open_channel_index_, closed_channel_index_,         &
            channel_index_, status_, l_
         real(dp) :: wavenumber, x, j_element_, jp_element_, n_element_,       &
            np_element_, ratio
         integer(int32) :: open_channels_indices(number_of_open_channels)
         integer(int32) :: closed_channels_indices(number_of_channels-number_of_open_channels)
         real(dp) ::  diag_n_matrix(number_of_channels, number_of_channels),   &
            diag_np_matrix(number_of_channels, number_of_channels),            &
            diag_j_matrix(number_of_channels, number_of_open_channels),        &
            diag_jp_matrix(number_of_channels,number_of_open_channels)
         !---------------------------------------------------------------------!
         ! diag_j_matrix   -  diagonal J matrix (Eqs. 5, 7)
         ! diag_jp_matrix  -  diagonal J`matrix (derivative of J)
         ! diag_n_matrix   -  diagonal N matrix (Eqs. 6, 8)
         ! diag_np_matrix  -  diagonal N`matrix (derivative of N)
         !---------------------------------------------------------------------!
         diag_j_matrix  = 0
         diag_jp_matrix = 0
         diag_n_matrix  = 0
         diag_np_matrix = 0
         !---------------------------------------------------------------------!
         open_channel_index_   = 0
         closed_channel_index_ = 0
         !---------------------------------------------------------------------!
         ! save indices to open and closed channels
         ! this is because channels might not be sorted eneregetically
         !---------------------------------------------------------------------!
         do channel_index_ = 1, number_of_channels
            if (is_open(elevel(channel_indices(channel_index_)))) then
               open_channel_index_ = open_channel_index_+1
               open_channels_indices(open_channel_index_) = channel_index_
            else
               closed_channel_index_ = closed_channel_index_+1
               closed_channels_indices(closed_channel_index_) = channel_index_
            endif
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices (Eqs. 5-6)
         ! open channels:
         !---------------------------------------------------------------------!
         do open_channel_index_ = 1, number_of_open_channels
            wavenumber = sqrt( wavenumber_squared_from_energy(                 &
               elevel(channel_indices(open_channels_indices(open_channel_index_)))))
            x  = wavenumber*r_
            l_ = channel_l_values(open_channels_indices(open_channel_index_))
            call riccati_bessel_j(                                             &
               l_, x, j_element_, jp_element_)
            diag_j_matrix(open_channel_index_,open_channel_index_)             &
               = (wavenumber)**(-0.5d0)*j_element_
            diag_jp_matrix(open_channel_index_,open_channel_index_)            &
               = (wavenumber)**(0.5d0)*jp_element_

            call riccati_bessel_y(l_, x, n_element_, np_element_)
            diag_n_matrix(open_channel_index_,open_channel_index_)             &
               = (wavenumber)**(-0.5d0)*n_element_
            diag_np_matrix(open_channel_index_,open_channel_index_)            &
               = (wavenumber)**(0.5d0)*np_element_
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices (Eqs. 7-8)
         ! closed channels:
         !---------------------------------------------------------------------!
         do closed_channel_index_ = 1,number_of_channels-number_of_open_channels  
            wavenumber = sqrt( abs( wavenumber_squared_from_energy(            &
               elevel(channel_indices(closed_channels_indices(closed_channel_index_))))))
            x  = wavenumber*r_
            l_ = channel_l_values(closed_channels_indices(closed_channel_index_))
            call modified_bessel_k_ratio(l_,x,ratio)
            !------------------------------------------------------------------!
            ! substitution for closed channels, (Eqs. 10 - 11)             
            !------------------------------------------------------------------!
            diag_n_matrix(number_of_open_channels+closed_channel_index_,       &
               number_of_open_channels+closed_channel_index_) = 1.d0
            diag_np_matrix(number_of_open_channels+closed_channel_index_,      &
               number_of_open_channels+closed_channel_index_) = wavenumber*ratio
         enddo
         !---------------------------------------------------------------------! -----------------------> consider a separate function
         call DGEMM('N','N',number_of_channels,number_of_channels,             &
            number_of_channels,1.0d0,log_der_matrix,number_of_channels,diag_n_matrix,       &
            number_of_channels,-1.d0,diag_np_matrix,number_of_channels)
         call DGEMM('N','N',number_of_channels,number_of_open_channels,        &
            number_of_channels,-1.0d0,log_der_matrix,number_of_channels,diag_j_matrix,      &
            number_of_channels,1.d0,diag_jp_matrix,number_of_channels)
         !---------------------------------------------------------------------!
         call DGESV(number_of_channels,number_of_open_channels,diag_np_matrix, &
            number_of_channels,diag_j_matrix,diag_jp_matrix,number_of_channels,status_)
         !---------------------------------------------------------------------!
         k_matrix = diag_jp_matrix(1:number_of_open_channels, 1:number_of_open_channels)
         !---------------------------------------------------------------------!
      end subroutine calculate_k_matrix
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
      subroutine calculate_s_matrix(number_of_open_channels, k_matrix,         &
         s_matrix_real, s_matrix_imag)
         !! calculates S-matrix from open-open portion of the K-matrix using
         !! Eq. (12) in "Solution of coupled equations" 
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_open_channels
            !! number of open channels
         real(dp), intent(in) :: k_matrix(number_of_open_channels,number_of_open_channels)
            !! K-matrix
         real(dp), intent(inout) :: s_matrix_real(number_of_open_channels,number_of_open_channels)
            !! (output) real part of the S-matrix
         real(dp), intent(inout) :: s_matrix_imag(number_of_open_channels,number_of_open_channels)
            !! (output) imaginary part of the S-matrix
         !---------------------------------------------------------------------!
         integer(int32) :: open_channel_index_1_, open_channel_index_2_
         real(dp) :: s_tmp_matrix(number_of_open_channels,number_of_open_channels)
         !---------------------------------------------------------------------!
         s_matrix_real = 0
         s_matrix_imag = 0
         !---------------------------------------------------------------------!
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,0.5d0,k_matrix,number_of_open_channels,    &
            k_matrix,number_of_open_channels,0.d0,s_tmp_matrix,number_of_open_channels)
         !---------------------------------------------------------------------!
         do open_channel_index_1_ = 1, number_of_open_channels
            s_tmp_matrix(open_channel_index_1_, open_channel_index_1_) =       &
               s_tmp_matrix(open_channel_index_1_, open_channel_index_1_) + 0.5d0
         enddo
         !---------------------------------------------------------------------!
         call invert_symmetric_matrix(s_tmp_matrix)
         call fill_symmetric_matrix(s_tmp_matrix, 'u')
         !---------------------------------------------------------------------!
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,-1.0d0,s_tmp_matrix,number_of_open_channels,&
            k_matrix,number_of_open_channels,0.d0,s_matrix_imag,number_of_open_channels)
         !---------------------------------------------------------------------!
         do open_channel_index_1_ = 1, number_of_open_channels
            do open_channel_index_2_ = 1, number_of_open_channels
               s_matrix_real(open_channel_index_1_, open_channel_index_2_) =   &
                  s_tmp_matrix(open_channel_index_1_, open_channel_index_2_)
            enddo
            s_matrix_real(open_channel_index_1_, open_channel_index_1_) =      &
               s_matrix_real(open_channel_index_1_, open_channel_index_1_) - 1.d0
         enddo
         !---------------------------------------------------------------------!
      end subroutine calculate_s_matrix
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
end module boundary_conditions_mod 
