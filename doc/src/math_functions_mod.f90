module math_functions_mod
   !! this module holds 4 types of functions:
   !! -- algebraic functions: percival_seaton_coefficient
   !! -- geometric functions: triangle_inequality_holds, is_sum_even, zero_projections_3j_condition
   !! -- bessel functions: groups functions: riccati_bessel_j, bessely and modified_bessel_k_ratio
   !!    that call special functions from special_functions.f90 library
   !! -- interpolation procedures: spline and ispline functions for interpolating data
   !! -- additional functions: rctj, rcty, envj, msta1, msta2, ikv, gamma from
   !!    special_functions library
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_error, write_warning,                &
      integer_to_character, float_to_character, time_count_summary
   use special_functions_mod, only: rctj, rcty
   !---------------------------------------------------------------------------!
   implicit none
	contains
   !---------------------------------------------------------------------------!
   !                           Algebraic functions
   !---------------------------------------------------------------------------!
   function percival_seaton_coefficient(j_, j_prime_, lambda_, omega_)         &
      result(percival_seaton_coefficient_)
      !! calculates Percival-Seaton coefficients (body-fixed variant)
      !! \begin{equation}
      !! \label{eq:algebraic_coeffs}
      !! g_{{\lambda},\gamma,\gamma'}^{Jp} = \delta_{\bar{\Omega},\bar{\Omega}'} (-1)^{\bar{\Omega}} \sqrt{(2j+1)(2j'+1)}
      !! \begin{pmatrix}
      !!   j & j' & \lambda \\ 0 & 0 & 0
      !! \end{pmatrix}
      !! \begin{pmatrix}
      !! j & j' & \lambda \\ \bar{\Omega} & -\bar{\Omega} & 0 \end{pmatrix}.
      !! \end{equation}
      !------------------------------------------------------------------------!
      use fwigxjpf, only: fwig3jj
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: j_
         !! pre-collisional rotational angular momentum
      integer(int32), intent(in) :: j_prime_
         !! post-collisional rotational angular momentum
      integer(int32), intent(in) :: omega_
         !! \\(\bar{\Omega}\\)
      integer(int32), intent(in) :: lambda_
         !! Legendre expansion coefficient \\( \lambda\\)
      real(dp) :: percival_seaton_coefficient_
         !! (out) result: percival seaton coefficient in the body-fixed frame
      !------------------------------------------------------------------------!
      percival_seaton_coefficient_ = (-1.0_dp)**(omega_) * sqrt(               &
         real((2 * j_ + 1)  * (2 * j_prime_ + 1), dp))                         &
         * fwig3jj(2* j_ ,   2* j_prime_  , 2* lambda_, 0, 0, 0)               &
         * fwig3jj(2* j_ ,   2* j_prime_  , 2* lambda_,                        &
            2 * omega_, -2 * omega_,     0)
      !------------------------------------------------------------------------!
   end function percival_seaton_coefficient
   !---------------------------------------------------------------------------!
   !                           Geometric functions
   !---------------------------------------------------------------------------!
   function triangle_inequality_holds(x, y, z) result(holds)
      !! check if the triangle inequality for 3 variables hols
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: x, y, z
         !! variables to check the triangle inequality
      logical :: holds
         !! (out) result: true/false
      !------------------------------------------------------------------------!
      holds = ( (x + y >= z) .and. (x + z >= y) .and. (y + z >= x) )
      !------------------------------------------------------------------------!
   end function triangle_inequality_holds
   !---------------------------------------------------------------------------!
   function is_sum_even(x, y, z) result(sum_even)
      !! checks if the sum of 3 integers is an even integer
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: x, y, z
         !! variables to check if the sum is even
      logical :: sum_even
         !! (out) result: true/false
      !------------------------------------------------------------------------!
      sum_even = (modulo(x + y + z, 2) == 0)
      !------------------------------------------------------------------------!
   end function is_sum_even
   !---------------------------------------------------------------------------!
   function zero_projections_3j_condition(x, y, z) result(is_valid)
      !! checks the condition for nonvanishing 3-j symbol with zero projections:
      !! triangle inequality on x,y,z and if the sum x+y+z is an even integer
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: x, y, z
         !! variables to check for 3-j symbol conditions
      logical :: is_valid
         !! (out) result: true/false if conditions are met
      !------------------------------------------------------------------------!
      is_valid = (triangle_inequality_holds(x, y, z) .and. is_sum_even(x, y, z))
      !------------------------------------------------------------------------!
   end function zero_projections_3j_condition
   !---------------------------------------------------------------------------!
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
      !! kind K_{l_ + 1/2}(x) and its first derivative (Eq. 8 in the
      !! "Solution of the coupled equations" section)
      !! Uses Temme's algorithm [N. M. Temme, J. Comput. Phys. 19 (1975) 324],
      !! implemented in "modified_bessel_temme_algorithm" subroutine;
      !! Unfortunately, the "ikv" function from special_functions
      !! library failed at large x_ values.
      !------------------------------------------------------------------------!
      integer(int32), intent(in) :: l_
         !! l - order of the function (without the 1/2 factor!)
      real(dp), intent(in) :: x_
         !! x - argument of the function
      real(dp), intent(inout) :: ratio_
         !! ratio of the modified Bessel function of the second kind to its derivative
      !------------------------------------------------------------------------!
      real(dp) :: order_, ck_, dk_, ek_
      !------------------------------------------------------------------------!
      order_ = real(l_, dp) + 0.5_dp
      call modified_bessel_temme_algorithm(order_, x_, ck_, dk_, ek_)
      ratio_ = dk_/ck_
      !------------------------------------------------------------------------!
   end subroutine modified_bessel_k_ratio
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   subroutine modified_bessel_temme_algorithm(v, x, ck, dk, ek)
      !! Implementation of the Temme's algorithm
      !! [N. M. Temme, J. Comput. Phys. 19 (1975) 324] to calculating
      !! modified Bessel functions of the second kind.
      !! This is a direct modernization of the "mbessk" subroutine
      !! in MOLSCAT:
      !! https://github.com/molscat/molscat/blob/master/source_code/rbessk.f
      !------------------------------------------------------------------------!
      real(dp), intent(in) :: v, x
      real(dp), intent(out) :: ck, dk, ek
      real(dp) :: a, b, c, d, e, f, g, h, p, q, s, sk, y, ak, ak1, ex, pi
      integer :: n, m, na, maxit
      real(dp), parameter :: eps = 1.0e-15_dp, xmin = 1.0_dp
      !------------------------------------------------------------------------!
      pi = acos(-1.0_dp)
      if (v < 0.0_dp .or. x <= 0.0_dp) then
         call write_error("modified_bessel_temme_algorithm: Invalid input values for v or x")
      endif

      na = int(v + 0.5_dp)
      a = v - na
      if (x < xmin) then
      ! Small x: Temme's series for small x
         b = x / 2.0_dp
         d = -log(b)
         e = a * d
         c = a * pi
         if (abs(c) < eps) then
          c = 1.0_dp
         else
          c = c / sin(c)
         endif
         if (abs(e) < eps) then
          s = 1.0_dp
         else
          s = sinh(e) / e
         endif
         e = exp(e)
         ! Compute the gamma function and its derivatives P and Q
         ! Replace RGAMMA(A,P,Q) with a modern equivalent if necessary
         g = e * rgamma(a, p, q)
         e = (e + 1.0_dp / e) / 2.0_dp
         f = c * (p * e + q * s * d)
         e = a * a
         p = 0.5_dp * g * c
         q = 0.5_dp / g
         c = 1.0_dp
         d = b * b
         ak = f
         ak1 = p
         do n = 1, maxit
          f = (f * n + p + q) / (n * n - e)
          c = c * d / n
          p = p / (n - a)
          q = q / (n + a)
          g = c * (p - n * f)
          h = c * f
          ak = ak + h
          ak1 = ak1 + g
          if (abs(h / ak) + abs(g / ak1) < eps) exit
         end do
         f = ak
         g = ak1 / b
         ex = 0.0_dp
      else
      ! Large x: Temme's PQ method for large x
         c = 0.25_dp - a * a
         g = 1.0_dp
         f = 0.0_dp
         e = x * cos(a * pi) / pi / eps
         do n = 1, maxit
          h = (2 * (n + x) * g - (n - 1 + c / n) * f) / (n + 1)
          f = g
          g = h
          if (h * n > e) exit
         end do

         p = f / g
         q = p
         b = x + x
         e = b - 2.0_dp
         do m = n, 1, -1
          p = (m - 1 + c / m) / (e + (m + 1) * (2.0_dp - p))
          q = p * (q + 1.0_dp)
         end do
         f = sqrt(pi / b) / (1.0_dp + q)
         g = f * (a + x + 0.5_dp - p) / x
         ex = x
      endif

      ! Upward recursion
      p = 0.0_dp
      if (na > 0) then
         y = 2.0_dp / x
         do n = 1, na
           h = y * (a + n) * g + f
           f = g
           g = h
           if (abs(f) > 4.0_dp) then
             p = p + 1.0_dp
             f = 0.0625_dp * f
             g = 0.0625_dp * g
           endif
         end do
      endif

      ck = f
      dk = (v / x) * f - g
      sk = sqrt(ck * ck + dk * dk)
      ck = ck / sk
      dk = dk / sk
      ek = log(sk) + p * log(16.0_dp) - ex
   end subroutine modified_bessel_temme_algorithm
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   function rgamma(x, odd, even) result(rgamma_val)
      !! Calculates 1/Gamma(1-X); modernized version of Molscat's
      !! rgamma function; see:
      !! https://github.com/molscat/molscat/blob/36fa8f93a92f851e9d84245dd6a972e2910541c5/source_code/rbesjy.f
      !!-----------------------------------------------------------------------!
      real(dp), intent(in) :: x
      real(dp), intent(out) :: odd, even
      real(dp) :: rgamma_val, x2, alfa, beta
      integer :: i
      real(dp), dimension(12), save :: b = [ &
      -0.283876542276024_dp, -0.076852840844786_dp, &
       0.001706305071096_dp,  0.001271927136655_dp, &
       0.000076309597586_dp, -0.000004971736704_dp, &
      -0.000000865920800_dp, -0.000000033126120_dp, &
       0.000000001745136_dp,  0.000000000242310_dp, &
       0.000000000009161_dp, -0.000000000000170_dp ]

      x2 = x * x * 8.0_dp
      alfa = -0.000000000000001_dp
      beta = 0.0_dp

      do i = 12, 2, -2
      beta = -(2 * alfa + beta)
      alfa = -beta * x2 - alfa + b(i)
      end do

      even = (beta / 2.0_dp + alfa) * x2 - alfa + 0.921870293650453_dp

      alfa = -0.000000000000034_dp
      beta = 0.0_dp

      do i = 11, 1, -2
      beta = -(2 * alfa + beta)
      alfa = -beta * x2 - alfa + b(i)
      end do

      odd = 2 * (alfa + beta)
      rgamma_val = odd * x + even
   end function rgamma
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
