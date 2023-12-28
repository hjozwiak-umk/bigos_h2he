module special_functions_mod
   !! This module holds rctj, rcty, envj, msta1, msta2, ikv, gamma from
   !! special_functions library, donwloaded from:
   !! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
   !!  Author:
   !!
   !!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
   !!    FORTRAN90 version by John Burkardt.
   !!
   !!  Reference:
   !!
   !!    Shanjie Zhang, Jianming Jin,
   !!    Computation of Special Functions,
   !!    Wiley, 1996,
   !!    ISBN: 0-471-11963-6,
   !!    LC: QA351.C45.
   !!--------------------------------------------------------------------------!
   implicit none
   !--------------------------------------------------------------------------!
   contains
   !--------------------------------------------------------------------------!
      subroutine rctj(n,x,nm,rj,dj)
      !------------------------------------------------------------------------! 
      !! computes Riccati-Bessel function of the first kind, and derivatives.
      !------------------------------------------------------------------------! 
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    18 July 2012
      !
      !  Author:
      !
      !    Shanjie Zhang,Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang,Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley,1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input,integer(kind = 4) N,the order of jn(x).
      !
      !    Input,real(kind = 8) X,the argument.
      !
      !    Output,integer(kind = 4) NM,the highest order computed.
      !
      !    Output,real(kind = 8) RJ(0:N),the values of x jn(x).
      !
      !    Output,real(kind = 8) DJ(0:N),the values of [x jn(x)]'.
      !------------------------------------------------------------------------! 
		  implicit none

		  integer(kind = 4) n

		  real(kind = 8) cs
		  real(kind = 8) dj(0:n)
		  real(kind = 8) f
		  real(kind = 8) f0
		  real(kind = 8) f1
		  integer(kind = 4) k
		  integer(kind = 4) m
        ! ----------------------------------------------------------------------
        ! commented because the functions are now in the same module
!		  integer(kind = 4) msta1 
!		  integer(kind = 4) msta2
        ! ----------------------------------------------------------------------
		  integer(kind = 4) nm
		  real(kind = 8) rj(0:n)
		  real(kind = 8) rj0
		  real(kind = 8) rj1
		  real(kind = 8) x

		  nm = n

		  if(abs(x) < 1.0D-100) then
			 do k = 0,n
				rj(k) = 0.0D+00
				dj(k) = 0.0D+00
			 end do
			 dj(0) = 1.0D+00
			 return
		  end if

		  rj(0) = sin(x)
		  rj(1) = rj(0) / x - cos(x)
		  rj0 = rj(0)
		  rj1 = rj(1)

		  if(2 <= n) then

			 m = msta1(x,200)

			 if(m < n) then
				nm = m
			 else
				m = msta2(x,n,15)
			 end if

			 f0 = 0.0D+00
			 f1 = 1.0D-100
			 do k = m,0,-1
				f =(2.0D+00 * k + 3.0D+00) * f1 / x - f0
				if(k <= nm) then
				  rj(k) = f
				end if
				f0 = f1
				f1 = f
			 end do

			 if(abs(rj1) < abs(rj0)) then
				cs = rj0 / f
			 else
				cs = rj1 / f0
			 end if

			 do k = 0,nm
				rj(k) = cs * rj(k)
			 end do

		  end if

		  dj(0) = cos(x)
		  do k = 1,nm
			 dj(k) = - k * rj(k) / x + rj(k-1)
		  end do

		  return
		end
      !------------------------------------------------------------------------!
		subroutine rcty(n,x,nm,ry,dy)
      !------------------------------------------------------------------------!
      !! computes Riccati-Bessel function of the second kind, and derivatives.
      !------------------------------------------------------------------------!
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    18 July 2012
      !
      !  Author:
      !
      !    Shanjie Zhang,Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang,Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley,1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input,integer(kind = 4) N,the order of yn(x).
      !
      !    Input,real(kind = 8) X,the argument.
      !
      !    Output,integer(kind = 4) NM,the highest order computed.
      !
      !    Output,real(kind = 8) RY(0:N),the values of x yn(x).
      !
      !    Output,real(kind = 8) DY(0:N),the values of [x yn(x)]'.
      !------------------------------------------------------------------------!
		  implicit none

		  integer(kind = 4) n

		  real(kind = 8) dy(0:n)
		  integer(kind = 4) k
		  integer(kind = 4) nm
		  real(kind = 8) rf0
		  real(kind = 8) rf1
		  real(kind = 8) rf2
		  real(kind = 8) ry(0:n)
		  real(kind = 8) x

		  nm = n

		  if(x < 1.0D-60) then
			 do k = 0,n
				ry(k) = -1.0D+300
				dy(k) = 1.0D+300
			 end do
			 ry(0) = -1.0D+00
			 dy(0) = 0.0D+00
			 return
		  end if

		  ry(0) = - cos(x)
		  ry(1) = ry(0) / x - sin(x)
		  rf0 = ry(0)
		  rf1 = ry(1)
		  do k = 2,n
			 rf2 =(2.0D+00 * k - 1.0D+00) * rf1 / x - rf0
			 if(1.0D+300 < abs(rf2)) then
				exit
			 end if
			 ry(k) = rf2
			 rf0 = rf1
			 rf1 = rf2
		  end do

		  nm = k - 1
		  dy(0) = sin(x)
		  do k = 1,nm
			 dy(k) = - k * ry(k) / x + ry(k-1)
		  end do

		  return
		end
      !------------------------------------------------------------------------!
      subroutine gamma ( x, ga )
      !------------------------------------------------------------------------!
      !! evaluates the Gamma function.
      !------------------------------------------------------------------------!
      !  Licensing:
      !
      !    The original FORTRAN77 version of this routine is copyrighted by 
      !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
      !    incorporate this routine into a user program that the copyright 
      !    is acknowledged.
      !
      !  Modified:
      !
      !    08 September 2007
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Shanjie Zhang, Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley, 1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) X, the argument.
      !    X must not be 0, or any negative integer.
      !
      !    Output, real ( kind = 8 ) GA, the value of the Gamma function.
      !------------------------------------------------------------------------!
        implicit none

        real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
          1.0D+00, &
          0.5772156649015329D+00, &
         -0.6558780715202538D+00, &
         -0.420026350340952D-01, &
          0.1665386113822915D+00, &
         -0.421977345555443D-01, &
         -0.96219715278770D-02, &
          0.72189432466630D-02, &
         -0.11651675918591D-02, &
         -0.2152416741149D-03, &
          0.1280502823882D-03, & 
         -0.201348547807D-04, &
         -0.12504934821D-05, &
          0.11330272320D-05, &
         -0.2056338417D-06, & 
          0.61160950D-08, &
          0.50020075D-08, &
         -0.11812746D-08, &
          0.1043427D-09, & 
          0.77823D-11, &
         -0.36968D-11, &
          0.51D-12, &
         -0.206D-13, &
         -0.54D-14, &
          0.14D-14, &
          0.1D-15 /)
        real ( kind = 8 ) ga
        real ( kind = 8 ) gr
        integer ( kind = 4 ) k
        integer ( kind = 4 ) m
        integer ( kind = 4 ) m1
        real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
        real ( kind = 8 ) r
        real ( kind = 8 ) x
        real ( kind = 8 ) z

        if ( x == aint ( x ) ) then

          if ( 0.0D+00 < x ) then
            ga = 1.0D+00
            m1 = int ( x ) - 1
            do k = 2, m1
              ga = ga * k
            end do
          else
            ga = 1.0D+300
          end if

        else

          if ( 1.0D+00 < abs ( x ) ) then
            z = abs ( x )
            m = int ( z )
            r = 1.0D+00
            do k = 1, m
              r = r * ( z - real ( k, kind = 8 ) )
            end do
            z = z - real ( m, kind = 8 )
          else
            z = x
          end if

          gr = g(26)
          do k = 25, 1, -1
            gr = gr * z + g(k)
          end do

          ga = 1.0D+00 / ( gr * z )

          if ( 1.0D+00 < abs ( x ) ) then
            ga = ga * r
            if ( x < 0.0D+00 ) then
              ga = - pi / ( x* ga * sin ( pi * x ) )
            end if
          end if

        end if

        return
      end
      !------------------------------------------------------------------------!
      subroutine ikv ( v, x, vm, bi, di, bk, dk )
      !------------------------------------------------------------------------!
      !! computes modified Bessel function Iv(x) and Kv(x) and their derivatives.
      !------------------------------------------------------------------------!
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    17 July 2012
      !
      !  Author:
      !
      !    Shanjie Zhang, Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang, Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley, 1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) V, the order of Iv(x) and Kv(x).
      !    V = N + V0.
      !
      !    Input, real ( kind = 8 ) X, the argument.
      !
      !    Output, real ( kind = 8 ) VM, the highest order computed.
      !
      !    Output, real ( kind = 8 ) BI(0:N), DI(0:N), BK(0:N), DK(0:N), the
      !    values of In+v0(x), In+v0'(x), Kn+v0(x), Kn+v0'(x).
      !------------------------------------------------------------------------!
        implicit none

        real ( kind = 8 ) a1
        real ( kind = 8 ) a2
        real ( kind = 8 ) bi(0:*)
        real ( kind = 8 ) bi0
        real ( kind = 8 ) bk(0:*)
        real ( kind = 8 ) bk0
        real ( kind = 8 ) bk1
        real ( kind = 8 ) bk2
        real ( kind = 8 ) ca
        real ( kind = 8 ) cb
        real ( kind = 8 ) cs
        real ( kind = 8 ) ct
        real ( kind = 8 ) di(0:*)
        real ( kind = 8 ) dk(0:*)
        real ( kind = 8 ) f
        real ( kind = 8 ) f1
        real ( kind = 8 ) f2
        real ( kind = 8 ) gan
        real ( kind = 8 ) gap
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k0
        integer ( kind = 4 ) m
        ! ----------------------------------------------------------------------
        ! commented because the functions are now in the same module
!		  integer(kind = 4) msta1 
!		  integer(kind = 4) msta2
        ! ----------------------------------------------------------------------
        integer ( kind = 4 ) n
        real ( kind = 8 ) pi
        real ( kind = 8 ) piv
        real ( kind = 8 ) r
        real ( kind = 8 ) r1
        real ( kind = 8 ) r2
        real ( kind = 8 ) sum
        real ( kind = 8 ) v
        real ( kind = 8 ) v0
        real ( kind = 8 ) v0n
        real ( kind = 8 ) v0p
        real ( kind = 8 ) vm
        real ( kind = 8 ) vt
        real ( kind = 8 ) w0
        real ( kind = 8 ) wa
        real ( kind = 8 ) ww
        real ( kind = 8 ) x
        real ( kind = 8 ) x2

        pi = 3.141592653589793D+00
        x2 = x * x
        n = int ( v )
        v0 = v - n
        if ( n == 0 ) then
          n = 1
        end if

        if ( x < 1.0D-100 ) then

          do k = 0, n
            bi(k) = 0.0D+00
            di(k) = 0.0D+00
            bk(k) = -1.0D+300
            dk(k) = 1.0D+300
          end do

          if ( v == 0.0D+00 ) then
            bi(0) = 1.0D+00
            di(1) = 0.5D+00
          end if

          vm = v
          return

        end if

        piv = pi * v0
        vt = 4.0D+00 * v0 * v0

        if ( v0 == 0.0D+00 ) then
          a1 = 1.0D+00
        else
          v0p = 1.0D+00 + v0
          call gamma ( v0p, gap )
          a1 = ( 0.5D+00 * x ) ** v0 / gap
        end if

        if ( x < 35.0D+00 ) then
          k0 = 14
        else if ( x < 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if
       
        if ( x <= 18.0D+00 ) then

          bi0 = 1.0D+00
          r = 1.0D+00
          do k = 1, 30
            r = 0.25D+00 * r * x2 / ( k * ( k + v0 ) )
            bi0 = bi0 + r
            if ( abs ( r / bi0 ) < 1.0D-15 ) then
              exit
            end if
          end do

          bi0 = bi0 * a1

        else

          ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
          sum = 1.0D+00
          r = 1.0D+00
          do k = 1, k0
            r = -0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
            sum = sum + r
          end do
          bi0 = ca * sum

        end if

        m = msta1 ( x, 200 )

        if ( m < n ) then
          n = m
        else
          m = msta2 ( x, n, 15 )
        end if

        f2 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 + f2
          if ( k <= n ) then
            bi(k) = f
          end if
          f2 = f1
          f1 = f
        end do

        cs = bi0 / f
        do k = 0, n
          bi(k) = cs * bi(k)
        end do

        di(0) = v0 / x * bi(0) + bi(1)
        do k = 1, n
          di(k) = - ( k + v0 ) / x * bi(k) + bi(k-1)
        end do

        if ( x <= 9.0D+00 ) then

          if ( v0 == 0.0D+00 ) then

            ct = - log ( 0.5D+00 * x ) - 0.5772156649015329D+00
            cs = 0.0D+00
            w0 = 0.0D+00
            r = 1.0D+00
            do k = 1, 50
              w0 = w0 + 1.0D+00 / k
              r = 0.25D+00 * r / ( k * k ) * x2
              cs = cs + r * ( w0 + ct )
              wa = abs ( cs )
              if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
                exit
              end if
              ww = wa
            end do

            bk0 = ct + cs

          else

            v0n = 1.0D+00 - v0
            call gamma ( v0n, gan )
            a2 = 1.0D+00 / ( gan * ( 0.5D+00 * x ) ** v0 )
            a1 = ( 0.5D+00 * x ) ** v0 / gap
            sum = a2 - a1
            r1 = 1.0D+00
            r2 = 1.0D+00
            do k = 1, 120
              r1 = 0.25D+00 * r1 * x2 / ( k * ( k - v0 ) )
              r2 = 0.25D+00 * r2 * x2 / ( k * ( k + v0 ) )
              sum = sum + a2 * r1 - a1 * r2
              wa = abs ( sum )
              if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
                exit
              end if
              ww = wa
            end do

            bk0 = 0.5D+00 * pi * sum / sin ( piv )

          end if

        else

          cb = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )
          sum = 1.0D+00
          r = 1.0D+00
          do k = 1, k0
            r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
            sum = sum + r
          end do
          bk0 = cb * sum

        end if

        bk1 = ( 1.0D+00 / x - bi(1) * bk0 ) / bi(0)
        bk(0) = bk0
        bk(1) = bk1
        do k = 2, n
          bk2 = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * bk1 + bk0
          bk(k) = bk2
          bk0 = bk1
          bk1 = bk2
        end do

        dk(0) = v0 / x * bk(0) - bk(1)
        do k = 1, n
          dk(k) = - ( k + v0 ) / x * bk(k) - bk(k-1)
        end do

        vm = n + v0

        return
      end
      !------------------------------------------------------------------------!
		function envj(n,x)
      !------------------------------------------------------------------------!
      !! utility function used by MSTA1 and MSTA2.
      !------------------------------------------------------------------------!
      !  Discussion:
      !
      !    ENVJ estimates -log(Jn(x)) from the estimate
      !    Jn(x) approx 1/sqrt(2*pi*n) *(e*x/(2*n))^n
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    14 January 2016
      !
      !  Author:
      !
      !    Shanjie Zhang,Jianming Jin
      !    Modifications suggested by Vincent Lafage,11 January 2016.
      !
      !  Reference:
      !
      !    Shanjie Zhang,Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley,1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input,integer(kind = 4) N,the order of the Bessel function.
      !
      !    Input,real(kind = 8) X,the absolute value of the argument.
      !
      !    Output,real(kind = 8) ENVJ,the value.
      !------------------------------------------------------------------------!
		  implicit none

		  real(kind = 8) envj
		  real(kind = 8) logten
		  integer(kind = 4) n
		  real(kind = 8) n_r8
		  real(kind = 8) r8_gamma_log
		  real(kind = 8) x
		!
		!  Original code
		!
		  if(.true.) then

			 envj = 0.5D+00 * log10(6.28D+00 * n) &
				- n * log10(1.36D+00 * x / n)
		!
		!  Modification suggested by Vincent Lafage.
		!
		  else

			 n_r8 = real(n,kind = 8)
			 logten = log(10.0D+00)
			 envj = r8_gamma_log(n_r8 + 1.0D+00) / logten - n_r8 * log10(x)

		  end if

		  return
		end
      !------------------------------------------------------------------------!
		function msta1(x,mp)
      !------------------------------------------------------------------------!
      !! determines a backward recurrence starting point for Jn(x).
      !------------------------------------------------------------------------!
      !  Discussion:
      !
      !    This procedure determines the starting point for backward  
      !    recurrence such that the magnitude of    
      !    Jn(x) at that point is about 10^(-MP).
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    08 July 2012
      !
      !  Author:
      !
      !    Shanjie Zhang,Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang,Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley,1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input,real(kind = 8) X,the argument.
      !
      !    Input,integer(kind = 4) MP,the negative logarithm of the 
      !    desired magnitude.
      !
      !    Output,integer(kind = 4) MSTA1,the starting point.
      !------------------------------------------------------------------------!
		  implicit none

		  real(kind = 8) a0
        !----------------------------------------------------------------------
        ! commented because the envj is in the same module
!		  real(kind = 8) envj
        !-----------------------------------------------------------------------
		  real(kind = 8) f
		  real(kind = 8) f0
		  real(kind = 8) f1
		  integer(kind = 4) it
		  integer(kind = 4) mp
		  integer(kind = 4) msta1
		  integer(kind = 4) n0
		  integer(kind = 4) n1
		  integer(kind = 4) nn
		  real(kind = 8) x

		  a0 = abs(x)
		  n0 = int(1.1D+00 * a0) + 1
		  f0 = envj(n0,a0) - mp
		  n1 = n0 + 5
		  f1 = envj(n1,a0) - mp
		  do it = 1,20       
			 nn = n1 -(n1 - n0) /(1.0D+00 - f0 / f1)                  
			 f = envj(nn,a0) - mp
			 if(abs(nn - n1) < 1) then
				exit
			 end if
			 n0 = n1
			 f0 = f1
			 n1 = nn
			 f1 = f
		  end do

		  msta1 = nn

		  return
		end
      !------------------------------------------------------------------------!
		function msta2(x,n,mp)
      !------------------------------------------------------------------------!
      !! determines a backward recurrence starting point for Jn(x).
      !------------------------------------------------------------------------!
      !  Discussion:
      !
      !    This procedure determines the starting point for a backward
      !    recurrence such that all Jn(x) has MP significant digits.
      !
      !    Jianming Jin supplied a modification to this code on 12 January 2016.
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program 
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    14 January 2016
      !
      !  Author:
      !
      !    Shanjie Zhang,Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang,Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley,1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input,real(kind = 8) X,the argument of Jn(x).
      !
      !    Input,integer(kind = 4) N,the order of Jn(x).
      !
      !    Input,integer(kind = 4) MP,the number of significant digits.
      !
      !    Output,integer(kind = 4) MSTA2,the starting point.
      !------------------------------------------------------------------------!
		  implicit none

		  real(kind = 8) a0
		  real(kind = 8) ejn
        !----------------------------------------------------------------------
        ! commented because the envj is in the same module
!		  real(kind = 8) envj
        !----------------------------------------------------------------------
		  real(kind = 8) f
		  real(kind = 8) f0
		  real(kind = 8) f1
		  real(kind = 8) hmp
		  integer(kind = 4) it
		  integer(kind = 4) mp
		  integer(kind = 4) msta2
		  integer(kind = 4) n
		  integer(kind = 4) n0
		  integer(kind = 4) n1
		  integer(kind = 4) nn
		  real(kind = 8) obj
		  real(kind = 8) x

		  a0 = abs(x)
		  hmp = 0.5D+00 * mp
		  ejn = envj(n,a0)

		  if(ejn <= hmp) then
			 obj = mp
		!
		!  Original code:
		!
		!   n0 = int(1.1D+00 * a0)
		!
		!  Updated code:
		!
			 n0 = int(1.1D+00 * a0) + 1
		  else
			 obj = hmp + ejn
			 n0 = n
		  end if

		  f0 = envj(n0,a0) - obj
		  n1 = n0 + 5
		  f1 = envj(n1,a0) - obj

		  do it = 1,20
			 nn = n1 -(n1 - n0) /(1.0D+00 - f0 / f1)
			 f = envj(nn,a0) - obj
			 if(abs(nn - n1) < 1) then
				exit
			 end if
			 n0 = n1
			 f0 = f1
			 n1 = nn
			 f1 = f
		  end do

		  msta2 = nn + 10

		  return
		end
      !------------------------------------------------------------------------!  
end module special_functions_mod
