module BOUNDARY_CONDITIONS
   !! This module contains subroutines that transform the asymptotic
   !! log-derivative matrix into the scattering S-matrix:
   !! 1. pfunc      - calculates the P coefficients from Eq. (6.19)
   !!  2. transsum   - calculates the sum over \bar{\Omega} in Eq. (6.17-6.18)
   !!  3. bftosfmat  - transforms a given BF-matrix to the SF frame of rerence
   !!  4. calculate_k_matrix - transforms log-derivative matrix into reactance matrix
   !     (Eq.6.42)
   !!  5. kmattosmat - transfroms reactance matrix into scattering S-matrix
   !     (Eq.6.44)
   !------------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use io_mod
   use fwigxjpf, only: fwig3jj
   use supplementary
   use supplementary_mod, only: write_warning
   use additional_mod, only: invert_symmetric_matrix, fill_symmetric_matrix
   !------------------------------------------------------------------------------!
   implicit none
   contains
   !------------------------------------------------------------------------------!
      function p_coeff(jtot_,j_,l_,omega_) result(p_coeff_)
         !! calculates the P coefficients from Eq. (6.19)
         !! It is called by numerov and log_derivative at the end of the propagation
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: jtot_
            !! total angular momentum
         integer(int32), intent(in) :: j_
            !! rotational quantum number
         integer(int32), intent(in) :: l_
            !! orbital angular momentum
         integer(int32), intent(in) :: omega_
            !! projection of j on the BF-Z axis
         real(dp) :: p_coeff_
            !! result (...)
         !---------------------------------------------------------------------!
         real(dp) :: delta_
         !---------------------------------------------------------------------!
         delta_ = 0.d0
         if (omega_ == 0) delta_ = 1.0_dp
         p_coeff_ = (-1.0_dp)**(jtot_+omega_) * dsqrt(2.0_dp)                  &
            * dsqrt(real(2*l_+1, dp)) * fwig3jj(2* j_, 2* jtot_, 2* l_,        &
            2 * omega_, -2* omega_, 0) / dsqrt(1.0_dp + delta_)
         !---------------------------------------------------------------------!
      end function p_coeff
   !---------------------------------------------------------------------------!
      subroutine transform_summation(number_of_channels, jtot_, v_, j_, vp_, jp_, l_, lp_,&
         channels_level_indices, channels_omega_values, bf_matrix, sf_element)
         !! performs the summation in Eqs. (6.17) and (6.18)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: jtot_
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
         real(dp) :: coeff_outer, coeff_inner, sum_outer, sum_inner
         !---------------------------------------------------------------------!
         sum_outer = 0.0_dp
         do ii = 1, number_of_channels
            if (v1array(channels_level_indices(ii)) /= v_ .or. j1array(channels_level_indices(ii)) /= j_) cycle
            coeff_outer = p_coeff(jtot_, j_, l_, channels_omega_values(ii))
            sum_inner = 0.0_dp
            do ij = 1, number_of_channels
               if (v1array(channels_level_indices(ij)) /= vp_ .or. j1array(channels_level_indices(ij)) /= jp_) cycle
               coeff_inner = p_coeff(jtot_, jp_, lp_, channels_omega_values(ij))
               sum_inner = sum_inner + coeff_inner * bf_matrix(ii,ij)
            end do
            sum_outer = sum_outer + coeff_outer * sum_inner
         end do
         sf_element = sum_outer
         !---------------------------------------------------------------------!
      end subroutine transform_summation
      !------------------------------------------------------------------------!
      subroutine bf_to_sf_transformation(number_of_channels,jtot_,             &
         channels_level_indices,channels_omega_values,channels_l_values,       &
         bf_matrix,sf_matrix)
         !! This subroutine takes an input matrix in the BF-frame and 
         !! transforms its elements to the SF-frame; calls transsum and pfunc
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_channels
            !! size of the basis
         integer(int32), intent(in) :: jtot_
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
         real(dp) :: coeff, sf_result
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels
            v1_ = v1array(channels_level_indices(ii))
            j1_ = j1array(channels_level_indices(ii))
            l_  = channels_l_values(ii)
            do ij = 1, number_of_channels
               v1p_ = v1array(channels_level_indices(ij))
               j1p_ = j1array(channels_level_indices(ij))
               lp_  = channels_l_values(ij)
               call transform_summation(number_of_channels, jtot_, v1_, j1_, v1p_, j1p_,  &
                  l_, lp_, channels_level_indices, channels_omega_values,      &
                  bf_matrix, sf_result)
               sf_matrix(ii, ij) = sf_result
            enddo
         enddo
         !---------------------------------------------------------------------!         
      end subroutine bf_to_sf_transformation
!------------------------------------------------------------------------------!
      subroutine calculate_k_matrix(number_of_channels, log_der_matrix,        &
         number_of_open_channels, channels_level_indices, channels_l_values,   &
         r_, k_matrix)
         !! implementation of Eq. 6.53
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
         integer(int32) :: index_open, index_closed, ii, iopen, iopen2,        &
            iclosed, status_
         real(dp) :: wavenumber, x, j_element_, jp_element_, n_element_, np_element_, ratio
         integer(int32) :: open_channels_indices(number_of_open_channels)
         integer(int32) :: closed_channels_indices(number_of_channels-number_of_open_channels)
         real(dp) ::  diag_n_matrix(number_of_channels, number_of_channels),   &
            diag_np_matrix(number_of_channels, number_of_channels),            &
            diag_j_matrix(number_of_channels, number_of_open_channels),        &
            diag_jp_matrix(number_of_channels,number_of_open_channels)
         !---------------------------------------------------------------------!
         ! diag_j_matrix   -  diagonal J matrix (Eqs. 6.37 & 6.39)
         ! diag_jp_matrix  -  diagonal J`matrix (derivative of J)
         ! diag_n_matrix   -  diagonal N matrix (Eqs. 6.30 & 6.40)
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
         !---------------------------------------------------------------------!
         do ii = 1, number_of_channels           
            if (elevel(channels_level_indices(ii)).le.ETOTAL()) then
               index_open = index_open+1
               open_channels_indices(index_open) = ii
            else
               index_closed = index_closed+1
               closed_channels_indices(index_closed) = ii
            endif
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices
         ! open channels:
         !---------------------------------------------------------------------!
         do iopen = 1, number_of_open_channels
            wavenumber = dsqrt(2*reducedmass*(ETOTAL()&
               -elevel(channels_level_indices(open_channels_indices(iopen)))))
            x = wavenumber*r_
            call BESSELJ(channels_l_values(open_channels_indices(iopen)),x,j_element_,jp_element_)
            diag_j_matrix(iopen,iopen) = (wavenumber)**(-0.5d0)*j_element_
            diag_jp_matrix(iopen,iopen) = (wavenumber)**(0.5d0)*jp_element_

            call BESSELY(channels_l_values(open_channels_indices(iopen)),x,n_element_,np_element_)
            diag_n_matrix(iopen,iopen) = (wavenumber)**(-0.5d0)*n_element_
            diag_np_matrix(iopen,iopen) = (wavenumber)**(0.5d0)*np_element_
         enddo
         !---------------------------------------------------------------------!
         ! Prepare J, J', N and N' matrices
         ! closed channels:
         !---------------------------------------------------------------------!
         do iclosed = 1,number_of_channels-number_of_open_channels
            wavenumber = dsqrt(dabs(2*reducedmass*(ETOTAL()&
               -elevel(channels_level_indices(closed_channels_indices(iclosed))))))
            x = wavenumber*r_
            call MODBESSELK(channels_l_values(closed_channels_indices(iclosed)),x,ratio)
            !------------------------------------------------------------------!
            ! substitution for closed channels, (Eqs. 6.44 - 6.45)             
            !------------------------------------------------------------------!
            diag_n_matrix(number_of_open_channels+iclosed,number_of_open_channels+iclosed) = 1.d0
            diag_np_matrix(number_of_open_channels+iclosed,number_of_open_channels+iclosed) = (wavenumber)*ratio
         enddo
         !---------------------------------------------------------------------!
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
          do iopen = 1,number_of_open_channels
             do iopen2 = 1,number_of_open_channels
                k_matrix(iopen,iopen2) = diag_jp_matrix(iopen,iopen2)
             enddo
          enddo
         !---------------------------------------------------------------------!
      end subroutine calculate_k_matrix
!------------------------------------------------------------------------------!
      subroutine calculate_s_matrix(number_of_open_channels, k_matrix,         &
         s_matrix_real, s_matrix_imag)
         !! implementation of Eq. 6.57
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
         integer(int32) :: iopen, iopen2
         real(dp) :: s_tmp_matrix(number_of_open_channels,number_of_open_channels)
         !---------------------------------------------------------------------!
         s_matrix_real = 0
         s_matrix_imag = 0
         !---------------------------------------------------------------------!
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,0.5d0,k_matrix,number_of_open_channels,    &
            k_matrix,number_of_open_channels,0.d0,s_tmp_matrix,number_of_open_channels)

         do iopen = 1, number_of_open_channels
            s_tmp_matrix(iopen, iopen) = s_tmp_matrix(iopen, iopen) + 0.5d0
         enddo

         !call INV_MAT(s_tmp_matrix, number_of_open_channels)
         call invert_symmetric_matrix(s_tmp_matrix)
         call fill_symmetric_matrix(s_tmp_matrix, 'u')
         
         call DGEMM('N','N',number_of_open_channels,number_of_open_channels,   &
            number_of_open_channels,-1.0d0,s_tmp_matrix,number_of_open_channels,&
            k_matrix,number_of_open_channels,0.d0,s_matrix_imag,number_of_open_channels)

         do iopen = 1, number_of_open_channels
            do iopen2 = 1, number_of_open_channels
               s_matrix_real(iopen, iopen2) = s_tmp_matrix(iopen, iopen2)
            enddo
            s_matrix_real(iopen, iopen) = s_matrix_real(iopen, iopen) - 1.d0
         enddo
!------------------------------------------------------------------------------!
      end subroutine calculate_s_matrix
!------------------------------------------------------------------------------!
      subroutine unitarity_check(number_of_open_channels,s_matrix_real,s_matrix_imag,totalcheck)
         !! checks the unitarity of the S-matrix (Eq. 6.15)
         !---------------------------------------------------------------------------!
         use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
         use io_mod
         implicit none
         !---------------------------------------------------------------------------!
         integer(int32), intent(in) :: number_of_open_channels
            !! number of open channels
         real(dp), intent(in) :: s_matrix_real(number_of_open_channels,number_of_open_channels)
            !! real part of the S-matrix
         real(dp), intent(in) :: s_matrix_imag(number_of_open_channels,number_of_open_channels)
            !! imaginary part of the S-matrix
         logical, intent(inout) :: totalcheck
            !! (output) if .true. unitary is fulfilled, .false. otherwise
         !---------------------------------------------------------------------------!
         character(len = 200) :: line_
         integer(int32) :: iopen, iopen2
         real(dp) :: unitarychecktmp, unitary_tol
         real(dp) :: unitarycheck(number_of_open_channels)
         !---------------------------------------------------------------------------!
         totalcheck  = .true.
         unitary_tol = 1e-6_dp 
         !---------------------------------------------------------------------------!
         if (prntlvl.ge.4) then
            call write_message("Check of the unitarity of the S-matrix:")
            call write_message(repeat(" ", 34) // "*** S-matrix elements: ***")
            call write_message("   ROW   COL" // repeat(" ", 15) // "S**2" //        &
               repeat(" ", 17) // "RE (S)" // repeat(" ", 17) // "IM (S)")
         endif
         !---------------------------------------------------------------------------!
         do iopen=1, number_of_open_channels
            unitarychecktmp = 0.0_dp
            do iopen2=1, number_of_open_channels
               unitarychecktmp = unitarychecktmp +  s_matrix_real(iopen,iopen2)**2.  &
                  + s_matrix_imag(iopen,iopen2)**2.
               !---------------------------------------------------------------------!
               if (prntlvl.ge.4) then
                  write(line_, "(1X,I5,1X,I5,8X,E15.8,8X,E15.8,8X,E15.8)") iopen,     &
                     iopen2, s_matrix_real(iopen,iopen2)**2.+s_matrix_imag(iopen,iopen2)**2.,  &
                     s_matrix_real(iopen,iopen2), s_matrix_imag(iopen,iopen2)
                  call write_message(line_)
               endif
               !---------------------------------------------------------------------!
            enddo
            unitarycheck(iopen) = unitarychecktmp
         enddo
         !---------------------------------------------------------------------------!
         do iopen=1, number_of_open_channels
            if (dabs(unitarycheck(iopen)-1.d0).gt.unitary_tol) then
               call write_warning("Unitary condition is not fulfilled for channel no.")
               write(line_, "(1X,I5,8X,E15.8)") iopen, unitarycheck(iopen)
               call write_message(trim(adjustl(line_)))
               call write_message("Consider increasing the STEPS parameter")
               totalcheck = .false.
            endif
         enddo

         if (prntlvl.ge.4) then
            call write_message(repeat(" ", 43) // "***")
            call write_message("Check of the unitarity of the S-matrix:")
            call write_message("   ROW" // repeat(" ", 12) // "sum(S**2)")
            do iopen=1, number_of_open_channels
               write(line_, "(1X,I5,8X,E15.8)") iopen, unitarycheck(iopen)
               call write_message(trim(adjustl(line_)))
            enddo
            call write_message(repeat(" ", 43) // "***")
         endif
         !---------------------------------------------------------------------------!
      end subroutine unitarity_check
!------------------------------------------------------------------------------!
end module BOUNDARY_CONDITIONS 
