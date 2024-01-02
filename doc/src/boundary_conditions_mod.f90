module boundary_conditions_mod
   !! This module contains subroutines that transform the asymptotic
   !! log-derivative matrix into the scattering S-matrix:
   !! 1. pfunc      - calculates the P coefficients from Eq. (6.19)
   !! 2. transsum   - calculates the sum over \bar{\Omega} in Eq. (6.17-6.18)
   !! 3. bftosfmat  - transforms a given BF-matrix to the SF frame of rerence
   !! 4. calculate_k_matrix - transforms log-derivative matrix into reactance matrix
   !     (Eq.6.42)
   !! 5. kmattosmat - transfroms reactance matrix into scattering S-matrix
   !     (Eq.6.44)
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
      subroutine calculate_single_SF_element(number_of_channels,               &
         total_angular_momentum_, v_, j_, vp_, jp_, l_, lp_,                   &
         channels_level_indices, channels_omega_values, bf_matrix, sf_element)
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
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         real(dp), intent(in) :: bf_matrix(number_of_channels, number_of_channels)
            !! matrix in the BF frame
         real(dp), intent(out) :: sf_element
            !! (output) matrix element in the SF frame
         !---------------------------------------------------------------------!
         integer(int32) :: ii, ij, omtmp, omptmp, v1tmp, j1tmp, v1ptmp, j1ptmp
         real(dp) :: p_coeff_outer, p_coeff_inner, sum_outer, sum_inner
         !---------------------------------------------------------------------!
         sum_outer = 0.0_dp
         do ii = 1, number_of_channels
            if (v1array(channels_level_indices(ii)) /= v_ .or.                 &
               j1array(channels_level_indices(ii)) /= j_) cycle
            p_coeff_outer = p_coeff(total_angular_momentum_, j_, l_,           &
               channels_omega_values(ii))
            sum_inner = 0.0_dp
            do ij = 1, number_of_channels
               if (v1array(channels_level_indices(ij)) /= vp_ .or.             &
                  j1array(channels_level_indices(ij)) /= jp_) cycle
               p_coeff_inner = p_coeff(total_angular_momentum_, jp_, lp_,      &
                  channels_omega_values(ij))
               sum_inner = sum_inner + p_coeff_inner * bf_matrix(ii,ij)
            end do
            sum_outer = sum_outer + p_coeff_outer * sum_inner
         end do
         sf_element = sum_outer
         !---------------------------------------------------------------------!
      end subroutine calculate_single_SF_element
      !------------------------------------------------------------------------!
      subroutine calculate_sf_matrix_from_bf_matrix(number_of_channels,        &
         total_angular_momentum_, channels_level_indices,                      &
         channels_omega_values,channels_l_values, bf_matrix, sf_matrix)
         !! takes as an input matrix in the body-fixed frame and transforms it 
         !! to the spec-fixed frame; iterates over all matrix elements
         !! and calls calculate_single_SF_element
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_omega_values(number_of_channels)
            !! holds all values of \bar{\Omega}
         integer(int32), intent(in) :: channels_l_values(number_of_channels)
            !! holds all values of l
         real(dp), intent(in) :: bf_matrix(number_of_channels, number_of_channels)
            !! matrix in the BF frame
         real(dp), intent(inout) :: sf_matrix(number_of_channels, number_of_channels)
            !! (output) matrix in the SF frame
         !---------------------------------------------------------------------!
         integer(int32) :: l_, lp_, omega_, omegap_, v1_, j1_, v1p_, j1p_, ii, ij
         real(dp) :: single_sf_element
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            v1_ = v1array(channels_level_indices(ii))
            j1_ = j1array(channels_level_indices(ii))
            l_  = channels_l_values(ii)
            do ij = 1, number_of_channels
               v1p_ = v1array(channels_level_indices(ij))
               j1p_ = j1array(channels_level_indices(ij))
               lp_  = channels_l_values(ij)
               call calculate_single_SF_element(number_of_channels,            &
                  total_angular_momentum_, v1_, j1_, v1p_, j1p_, l_, lp_,      &
                  channels_level_indices, channels_omega_values, bf_matrix,    &
                  single_sf_element)
               sf_matrix(ii, ij) = single_sf_element
            enddo
         enddo
         !---------------------------------------------------------------------!         
      end subroutine calculate_sf_matrix_from_bf_matrix
      !------------------------------------------------------------------------!
      subroutine calculate_k_matrix(number_of_channels, log_der_matrix,        &
         number_of_open_channels, channels_level_indices, channels_l_values,   &
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
         integer(int32), intent(in) :: channels_level_indices(number_of_channels)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channels_l_values(number_of_channels)
            !! holds all values of l
         real(dp), intent(in) :: r_
            !! Rmax
         real(dp), intent(inout) :: k_matrix(number_of_open_channels, number_of_open_channels)
            !! K-matrix
         !---------------------------------------------------------------------!
         integer(int32) :: index_open, index_closed, ii, i_open, i_open2,      &
            i_closed, status_
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
         index_open   = 0
         index_closed = 0
         !---------------------------------------------------------------------!
         ! save indices to open and closed channels
         ! this is because channels might not be sorted eneregetically
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            if (is_open(elevel(channels_level_indices(ii)))) then
               index_open = index_open+1
               open_channels_indices(index_open) = ii
            else
               index_closed = index_closed+1
               closed_channels_indices(index_closed) = ii
            endif
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices (Eqs. 5-6)
         ! open channels:
         !---------------------------------------------------------------------!
         do i_open = 1, number_of_open_channels
            wavenumber = wavenumber_from_energy(                               &
               elevel(channels_level_indices(open_channels_indices(i_open))))
            x = wavenumber*r_
            call riccati_bessel_j(                                             &
               channels_l_values(open_channels_indices(i_open)), x, j_element_,&
               jp_element_)
            diag_j_matrix(i_open,i_open)  = (wavenumber)**(-0.5d0)*j_element_
            diag_jp_matrix(i_open,i_open) = (wavenumber)**(0.5d0)*jp_element_

            call riccati_bessel_y(                                             &
               channels_l_values(open_channels_indices(i_open)), x, n_element_,&
               np_element_)
            diag_n_matrix(i_open,i_open)  = (wavenumber)**(-0.5d0)*n_element_
            diag_np_matrix(i_open,i_open) = (wavenumber)**(0.5d0)*np_element_
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices (Eqs. 7-8)
         ! closed channels:
         !---------------------------------------------------------------------!
         do i_closed = 1,number_of_channels-number_of_open_channels  
            wavenumber = wavenumber_from_energy(                               &
               elevel(channels_level_indices(closed_channels_indices(i_closed))))
            x = wavenumber*r_
            call modified_bessel_k_ratio(channels_l_values(closed_channels_indices(i_closed)),x,ratio)
            !------------------------------------------------------------------!
            ! substitution for closed channels, (Eqs. 10 - 11)             
            !------------------------------------------------------------------!
            diag_n_matrix(number_of_open_channels+i_closed,number_of_open_channels+i_closed) = 1.d0
            diag_np_matrix(number_of_open_channels+i_closed,number_of_open_channels+i_closed) = (wavenumber)*ratio
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
      !------------------------------------------------------------------------!
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
         integer(int32) :: i_open, i_open2
         real(dp) :: s_tmp_matrix(number_of_open_channels,number_of_open_channels)
         !---------------------------------------------------------------------!
         s_matrix_real = 0
         s_matrix_imag = 0
         !---------------------------------------------------------------------!
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,0.5d0,k_matrix,number_of_open_channels,    &
            k_matrix,number_of_open_channels,0.d0,s_tmp_matrix,number_of_open_channels)
         !---------------------------------------------------------------------!
         do i_open = 1, number_of_open_channels
            s_tmp_matrix(i_open, i_open) = s_tmp_matrix(i_open, i_open) + 0.5d0
         enddo
         !---------------------------------------------------------------------!
         call invert_symmetric_matrix(s_tmp_matrix)
         call fill_symmetric_matrix(s_tmp_matrix, 'u')
         !---------------------------------------------------------------------!
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,-1.0d0,s_tmp_matrix,number_of_open_channels,&
            k_matrix,number_of_open_channels,0.d0,s_matrix_imag,number_of_open_channels)
         !---------------------------------------------------------------------!
         do i_open = 1, number_of_open_channels
            do i_open2 = 1, number_of_open_channels
               s_matrix_real(i_open, i_open2) = s_tmp_matrix(i_open, i_open2)
            enddo
            s_matrix_real(i_open, i_open) = s_matrix_real(i_open, i_open) - 1.d0
         enddo
         !---------------------------------------------------------------------!
      end subroutine calculate_s_matrix
      !------------------------------------------------------------------------!
      subroutine unitarity_check(number_of_open_channels,s_matrix_real,        &
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
         character(len = 200) :: line_
         integer(int32) :: channel_index
         real(dp) :: sum_of_squares(number_of_open_channels)
         !---------------------------------------------------------------------!
         is_unitary  = .true.
         !---------------------------------------------------------------------!
         call write_header("unitarity")
         !---------------------------------------------------------------------!
         do channel_index=1, number_of_open_channels
            sum_of_squares(channel_index) =                                    &
               sum((s_matrix_real(channel_index, :)**2)                        &
               + (s_matrix_imag(channel_index, :)**2))
         enddo
         !---------------------------------------------------------------------!
         do channel_index=1, number_of_open_channels
            if (dabs(sum_of_squares(channel_index)-1.d0).gt.unitary_tolerance) then
               call write_warning("Unitary condition is not fulfilled for channel no.")
               write(line_, "(1X,I5,8X,E15.8)") channel_index, sum_of_squares(channel_index)
               call write_message(trim(adjustl(line_)))
               call write_message("Consider increasing the STEPS parameter")
               is_unitary = .false.
            endif
         enddo
         !---------------------------------------------------------------------!
         if (is_unitary) then
            call write_message("S-matrix unitary condition fulfilled")
            call write_message(repeat(" ", 43) // "***")
         endif
         !---------------------------------------------------------------------!
      end subroutine unitarity_check
      !------------------------------------------------------------------------!
end module boundary_conditions_mod 
