module math_functions_mod
   !! this module holds 3 types of functions:
   !! -- geometric functions: triangle_inequality_holds and is_sum_even
   !! -- bessel functions: groups functions: riccati_bessel_j, bessely and modified_bessel_k_ratio
   !!    that call special functions from special_functions.f90 library
   !! -- interpolation procedures: spline and ispline functions for interpolating data
   !! -- additional functions: rctj, rcty, envj, msta1, msta2, ikv, gamma from
   !!    special_functions library
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_warning,                &
      integer_to_character, float_to_character, time_count_summary
   use special_functions_mod, only: rctj, rcty, envj, msta1, msta2, ikv, gamma
   !---------------------------------------------------------------------------!
   implicit none
	contains
   !---------------------------------------------------------------------------!
   !                           Geometric functions
   !---------------------------------------------------------------------------!
   function triangle_inequality_holds(x, y, z) result(triang)
      !! check if the triangle inequality for 3 variables hols
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: x, y, z
         !! variables to check the triangle inequality
      integer(int32) :: triang
         !! (out) result: 1 = true, 0 = false
      !------------------------------------------------------------------------!
      triang = 0
      if (x + y >= z .and. x + z >= y .and. y + z >= x) then
         triang = 1
      endif
      !------------------------------------------------------------------------!
   end function triangle_inequality_holds
   !---------------------------------------------------------------------------!
   function is_sum_even(x, y, z) result(sum_even)
      !! checks if the sum of 3 integers is an even integer
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: x, y, z
         !! variables to check if the sum is even
      integer(int32) :: sum_even
         !! (out) result: 1 = true, 0 = false
      !------------------------------------------------------------------------!
      sum_even = merge(1, 0, modulo(x + y + z, 2) == 0)
      !------------------------------------------------------------------------!
   end function is_sum_even
   !---------------------------------------------------------------------------!
   !                             Bessel functions
   ! these functions handle calling to specific subroutines from
   ! special_functions.f90 library
   !---------------------------------------------------------------------------!
   subroutine riccati_bessel_j(l_,x_,j_,jp_)
      !! calculates the Riccati-Bessel function of the first kind and its
      !! first derivative. Calls the rctj function from special_functions.f90
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: l_
         !! l - order of the Riccati-Bessel function of the first kind
      real(dp), intent(in) :: x_
         !! x - argument of the Riccati-Bessel function of the first kind
      real(dp), intent(inout) :: j_
         !! j_{l} (x) - Riccati-Bessel funciton of the first kind
      real(dp), intent(inout) :: jp_
         !! j_{l}' (x) - derivative of the Riccati-Bessel funciton of the first kind
      !------------------------------------------------------------------------!
      integer(int32) :: highest_order_
      real(dp), dimension(l_+1) :: rj_, dj_
      !------------------------------------------------------------------------!
      if(l_ == 0) then
         call rctj(l_+1, x_, highest_order_, rj_, dj_)
      else
         call rctj(l_, x_, highest_order_, rj_, dj_)
      endif

      if (highest_order_ < l_) then
         !---------------------------------------------------------------------!
         call write_warning("riccati_bessel_j: maximum order of Riccati-Bessel function:"&
         // trim(adjustl(integer_to_character(highest_order_))) // "is smaller than " //  &
         "requested order l = " // trim(adjustl(integer_to_character(l_))) )
         !---------------------------------------------------------------------!
         j_  = rj_(highest_order_)
         jp_ = dj_(highest_order_)
      else 
         j_  = rj_(l_+1)
         jp_ = dj_(l_+1)
      endif
      !------------------------------------------------------------------------!
   end subroutine riccati_bessel_j
   !---------------------------------------------------------------------------!
   subroutine riccati_bessel_y(l_, x_, y_, yp_)
      !! calculates the Riccati-Bessel function of the second kind and its
      !! first derivative. Calls the rcty function from special_functions.f90
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: l_
         !! l - order of the Riccati-Bessel function of the second kind
      real(dp), intent(in) :: x_
         !! x - argument of the Riccati-Bessel function of the second kind
      real(dp), intent(inout) :: y_
         !! y_{l} (x) - Riccati-Bessel funciton of the second kind
      real(dp), intent(inout) :: yp_
         !! y_{l}' (x) - derivative of the Riccati-Bessel funciton of the second kind
      !------------------------------------------------------------------------!
      integer(int32) :: highest_order_
      real(dp), dimension(l_+1) :: ry_, dy_
      !------------------------------------------------------------------------!
      if(l_ == 0) then
         call rcty(l_+1, x_, highest_order_, ry_, dy_)
      else
         call rcty(l_, x_, highest_order_, ry_, dy_)
      endif
      y_  = ry_(l_ + 1)
      yp_ = dy_(l_ + 1)
      !------------------------------------------------------------------------!
      if (highest_order_ < l_) then
         !---------------------------------------------------------------------!
         call write_warning("riccati_bessel_j: maximum order of Riccati-Bessel function:"&
         // trim(adjustl(integer_to_character(highest_order_))) // "is smaller than " //  &
         "requested order l = " // trim(adjustl(integer_to_character(l_))) )
         !---------------------------------------------------------------------!
         y_  = ry_(highest_order_)
         yp_ = dy_(highest_order_)
      else 
         y_  = ry_(l_ + 1)
         yp_ = dy_(l_ + 1)
      endif
      !------------------------------------------------------------------------!
   end subroutine riccati_bessel_y
   !---------------------------------------------------------------------------!
   subroutine modified_bessel_k_ratio(l_, x_, ratio_)
      !! calculates the ratio of the modified Bessel function of the second
      !! kind K_{l_ + 1/2}(x) and its first derivative (Eq. ...)
      !! Uses ikv function from special_functions.f90
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: l_
         !! l - order of the function (without the 1/2 factor!)
      real(dp), intent(in) :: x_
         !! x - argument of the function
      real(dp), intent(inout) :: ratio_
         !! ratio of the modified Bessel function of the second kind to its derivative
      !------------------------------------------------------------------------!
      integer(int32) :: highest_order_
      real(dp) :: order_, highest_order_real_
      real(dp), dimension(l_+1) :: bi_arr_, di_arr_, bk_arr_, dk_arr_
		!------------------------------------------------------------------------!
      order_ = real(l_, dp) + 0.5_dp
      call ikv(order_, x_, highest_order_real_, bi_arr_, di_arr_, bk_arr_, dk_arr_)
      highest_order_ = nint(highest_order_real_ - 0.5_dp)
      !------------------------------------------------------------------------!
      if (highest_order_ < l_) then
         !---------------------------------------------------------------------!
         call write_warning("modified_bessel_k_ratio: maximum order of modified Bessel function:"&
         // trim(adjustl(integer_to_character(highest_order_))) // "is smaller than " //  &
         "requested order l = " // trim(adjustl(integer_to_character(l_))) )
         !---------------------------------------------------------------------!
         ratio_ = dk_arr_(highest_order_) / bk_arr_(highest_order_)
      else 
         ratio_ = dk_arr_(l_ + 1) / bk_arr_(l_ + 1)
      endif
      !------------------------------------------------------------------------!
   end subroutine modified_bessel_k_ratio
   !---------------------------------------------------------------------------!
   !                          Interpolation procedures
   !---------------------------------------------------------------------------!
   subroutine spline (N_, x_, y_, b_, c_, d_)
      !! determines b, c and d coefficients of the cubic spline function
      !! y(x) = y_i + b_i * dx + c_i * dx^2 + d_i * dx^3,
      !! where dx = x - x_i, and x_i <= x < x_i+1.
      !! The algorithm is based on
      !! Gerald, C., and Wheatley, P., "Applied Numerical Analysis", Addison-Wesley, 1994.
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: N_
         !! number of grid points
      real(dp), intent(in) :: x_(N_)
         !! grid points (ascending order)
      real(dp), intent(in) :: y_(N_)
         !! tabulated values
      real(dp), intent(out) :: b_(N_), c_(N_), d_(N_)
         !! arrays with coefficients of the spline function
      !------------------------------------------------------------------------!
      integer(int32) :: i_, j_
      real(dp) :: w_
      real(dp) :: diff_x_(N_-1)
      !------------------------------------------------------------------------!
      ! check if the number of points is larger than 4
      !------------------------------------------------------------------------!
      if (N_ < 4) then
         call write_error("spline function called with " //                    &
            trim(adjustl(integer_to_character(N_))) // " points")
      endif
      !------------------------------------------------------------------------!
      ! check if x is sorted in ascending order
      !------------------------------------------------------------------------!
      do i_ = 2, N_
         if (x_(i_) <= x_(i_-1)) then
            call write_error("spline: x values are not in ascending order " // &
               " at index " // trim(adjustl(integer_to_character(i_))))
         endif
      end do
      !------------------------------------------------------------------------!
      diff_x_ = x_(2:N_) - x_(1:N_-1)
      !------------------------------------------------------------------------!
      b_(2:N_-1) = 2.0_dp * (diff_x_(1:N_-2) + diff_x_(2:N_-1))
      b_(1)  = - diff_x_(1)
      b_(N_) = - diff_x_(N_-1)

      c_(2:N_-1) = ( y_(3:N_) - y_(2:N_-1) ) / diff_x_(2:N_-1)                 &
         - ( y_(2:N_-1) - y_(1:N_-2) ) / diff_x_(1:N_-2)

      c_(1)  = c_(3)/(x_(4)-x_(2)) - c_(2)/(x_(3)-x_(1))
      c_(N_) = c_(N_-1)/(x_(N_)-x_(N_-2)) - c_(N_-2)/(x_(N_-1)-x_(N_-3))

      c_(1)  = c_(1)/(x_(4)-x_(1))*diff_x_(1)**2
      c_(N_) = -c_(N_)/(x_(N_)-x_(N_-3))*diff_x_(N_-1)**2

      do i_ = 2, N_
         w_ = diff_x_(i_-1)/b_(i_-1)
         b_(i_) = b_(i_) - w_*diff_x_(i_-1)
         c_(i_) = c_(i_) - w_*c_(i_-1)
      end do

      c_(N_) = c_(N_) / b_(N_)
      do j_ = 1, N_-1
         i_ = N_-j_
         c_(i_) = (c_(i_) - diff_x_(i_)*c_(i_+1)) / b_(i_)
      end do

      b_(1:N_-1) = ( y_(2:N_) - y_(1:N_-1) ) / diff_x_(1:N_-1)                 &
         - ( 2.0_dp * c_(1:N_-1) + c_(2:N_) ) * diff_x_(1:N_-1)

      d_(1:N_-1) = ( c_(2:N_) - c_(1:N_-1) ) / diff_x_(1:N_-1)

      c_ = c_ * 3.0_dp
        
   end subroutine spline
   !---------------------------------------------------------------------------!
   function ispline(u_, N_, x_, y_, b_, c_, d_) result(spl_result)
      !! returns interpolated value at guven u_ point
      !! number of points and ascending order of x is not checked since
      !! ispline is called after "spline" where these checks are done
      !------------------------------------------------------------------------!
      real(dp), intent(in) :: u_
         !! point at which the tabulated value is interpolated
      integer(int32), intent(in) :: N_
         !! number of grid points
      real(dp), intent(in) :: x_(N_)
         !! grid points
      real(dp), intent(in) :: y_(N_)
         !! tabulated values
      real(dp), intent(in) :: b_(N_), c_(N_), d_(N_)
         !! arrays with coefficients of the spline function
      real(dp) :: spl_result
         !! interpolated value at u_
      !------------------------------------------------------------------------!
      integer(int32) :: k_, l_, mid_
      real(dp) :: dx_
      !------------------------------------------------------------------------!
      if (u_ > x_(N_)) then
         call write_warning("ispline: point u_ = " //                          &
            trim(adjustl(float_to_character(u_))) // " exceeds the original "  &
            // "grid: x(N) = " // trim(adjustl(float_to_character(x_(N_)))))
         spl_result = y_(N_)
      else if (u_ < x_(1) ) then
         call write_warning("ispline: point u_ = " //                          &
            trim(adjustl(float_to_character(u_))) // " exceeds the original "  &
            // "grid: x(1) = " // trim(adjustl(float_to_character(x_(1)))))
         spl_result = y_(1)
      else
         !---------------------------------------------------------------------!
         l_ = 1
         k_ = N_+1
         do while (k_ > l_+1)

            mid_ = nint( (l_ + k_) / 2.0_dp )

            if (x_(mid_) > u_) then
                k_ = mid_
            else
                l_ = mid_
            endif

         end do

         dx_ = u_ - x_(l_)
         spl_result = y_(l_) + dx_ * (b_(l_) + dx_ * (c_(l_) + d_(l_) * dx_))
         !---------------------------------------------------------------------!
      endif
      !------------------------------------------------------------------------!
    	end function ispline
!------------------------------------------------------------------------------!
end module math_functions_mod
