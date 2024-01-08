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
