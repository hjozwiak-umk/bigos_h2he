module supplementary
	contains
   !------------------------------------------------------------------------------!
   function triangle_inequality_holds(x, y, z) result(triang)
      integer, intent(in) :: x, y, z
      integer :: triang

      triang = 0
      if (x + y >= z .and. x + z >= y .and. y + z >= x) then
         triang = 1
      endif
   end function triangle_inequality_holds
!------------------------------------------------------------------------------!
   function is_sum_even(x, y, z) result(sum_even)
      integer, intent(in) :: x, y, z
      integer :: sum_even

      sum_even = merge(1, 0, modulo(x + y + z, 2) == 0)
   end function is_sum_even
!------------------------------------------------------------------------------!
!********************************************************************
!********************************************************************
! Functions and subroutines below are taken from Molscat or from
! special_functions.f90 library
!********************************************************************
!********************************************************************
		subroutine BESSELJ(ll,X,BJ,BJP)
! This subroutines gives the Riccati-Bessel function of the first kind
! as well as its first derivative
! Based on the rctj function from special_functions.f90 library
! 
! Variables:
!
! Integer,in :: ll
! Double precision,in :: X
!
! Double precision,out :: BJ - Riccati-Bessel function \hat{j}_{ll}(X)
! Double precision,out :: BJP - First derivative of \hat{j}_{ll}(X)
			
			implicit none

			integer :: ll,NM
			double precision :: X,BJ,BJP
			double precision,dimension(ll+1) :: RJ,DJ

			if(ll.eq.0) then
				call rctj(ll+1,X,NM,RJ,DJ)
			else
				call rctj(ll,X,NM,RJ,DJ)
			endif

			if (NM.lt.ll) then
				BJ=RJ(NM)
				BJP=DJ(NM)
			else 
				BJ=RJ(ll+1)
				BJP=DJ(ll+1)
			endif

			return
			
		end subroutine BESSELJ
!********************************************************************
!********************************************************************
		subroutine BESSELY(ll,X,BY,BYP)
! This subroutines gives the Riccati-Bessel function of the second kind
! as well as its first derivative
! Based on the rcty function from special_functions.f90 library
! 
! Variables:
!
! Integer,in :: ll
! Double precision,in :: X
!
! Double precision,out :: BY - Riccati-Bessel function \hat{y}_{ll}(X)
! Double precision,out :: BYP - First derivative of \hat{y}_{ll}(X)
			
			implicit none

			integer :: ll,NM
			double precision :: X,BY,BYP
			double precision,dimension(ll+1) :: RY,DY

			if(ll.eq.0) then
				call rcty(ll+1,X,NM,RY,DY)
			else
				call rcty(ll,X,NM,RY,DY)
			endif
			BY=RY(ll+1)
			BYP=DY(ll+1)
			if (NM.lt.ll) then
				BY=RY(NM)
				BYP=DY(NM)
			else 
				BY=RY(ll+1)
				BYP=DY(ll+1)
			endif
			return
			
		end subroutine BESSELY
!********************************************************************
!********************************************************************
		subroutine MODBESSELK(ll,X,RATIO)
! This subroutines gives the Modified Bessel function/
! Bessel function of the third kind as well as its first derivative
! Based on the mbessk function from Molscat
! 
! Variables:
!
! Integer,in :: ll
! Double precision,in :: X
!
! Double precision CK
! Double precision DK
! Double precision EK
!
! Double precision,out :: BK - Mod. Bessel function
! K_{ll+1/2}(X) = CK * EXP(EK)
! Double precision,out :: BKP - First derivative of K_{ll+1/2}(X)
! K`_{ll+1/2}(X) = DK * EXP(EK)

			implicit none

			integer :: ll
			double precision :: V,RATIO,CK,DK,EK,X
			
			V=dfloat(ll)+0.5d0
			call mbessk(V,X,CK,DK,EK)
			RATIO=DK/CK
			return
			
		end subroutine MODBESSELK
!**************************************************************************************
    	SUBROUTINE SPLINE (x, y, b, c, d, N) 
!==============================================================================
! function looking for values of coefficients b, c and d 
! from cubic spline interpolation in from:
!           y(x) = y_i + b_i * dx + c_i * dx^2 + d_i * dx^3
!                   dx = x - x_i
! where i is defined by:
!                   x_i <= x < x_i+1
! alghoritm based on
!------------------------------------------------------------------------------
! input:
!   x        -  grid points (necessarily in raising order)  
!   y        -  value of function at grid points            
!   N        -  number of grid points (size of x)           
!
! output:
!   b, c, d  -  array with coefficients for                 
!               spline function
!
! comments:
!   alghoritm is written based on book:
!   Gerald, C., and Wheatley, P., "Applied Numerical Analysis", Addison-Wesley, 1994. 
!==============================================================================
    	implicit none

!      ---------------------------INPUT----------------------------------
        integer, intent(in) :: N
        double precision, dimension(N), intent(in) :: x, y

!      ---------------------------OUTPUT---------------------------------
        double precision, dimension(N), intent(out) :: b, c, d

!      ------------------SUPPLEMENTARY VARIABLES-------------------------
        double precision, dimension(N-1) :: diff_x
        double precision :: w
        integer :: i, j

!~~~~~~~~~~~~~~~~~~~~~~~~~~~CALCULATING S MATRIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        diff_x = x(2:N) - x(1:N-1)



        b(2:N-1) = 2.0_8*(diff_x(1:N-2) + diff_x(2:N-1))
        b(1) = -diff_x(1)
        b(N) = -diff_x(N-1)



        c(2:N-1) = ( y(3:N) - y(2:N-1) ) / diff_x(2:N-1) &
                - ( y(2:N-1) - y(1:N-2) ) / diff_x(1:N-2)

        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(N) = c(N-1)/(x(N)-x(N-2)) - c(N-2)/(x(N-1)-x(N-3))

        c(1) = c(1)/(x(4)-x(1))*diff_x(1)**2
        c(N) = -c(N)/(x(N)-x(N-3))*diff_x(N-1)**2

    
    
 !~~~~~~~~~~~~~~~~~~~~~~~~TRIDIAGONAL MATRIX ALGHORITM~~~~~~~~~~~~~~~~~~~~~~~~~~

        do i = 2,N
            w = diff_x(i-1)/b(i-1)
            b(i) = b(i) - w*diff_x(i-1)
            c(i) = c(i) - w*c(i-1)
        end do

        c(N) = c(N)/b(N)
        do j = 1,N-1
            i = N-j
            c(i) =(c(i) - diff_x(i)*c(i+1))/b(i)
        end do

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~SETTING COEFFICIENTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        b(1:N-1) = ( y(2:N) - y(1:N-1) ) / diff_x(1:N-1) &
                    - ( 2.0_8 * c(1:N-1) + c(2:N) ) * diff_x(1:N-1)


        d(1:N-1) = ( c(2:N) - c(1:N-1) ) / diff_x(1:N-1)


        c = c * 3.0_8
        
   	 END SUBROUTINE SPLINE
!**************************************************************************************
   	 DOUBLE PRECISION FUNCTION ISPLINE (u, x, y, b, c, d, N)
!==============================================================================
! function giving value of Spline in form:
!           y(x) = y_i + b_i * dx + c_i * dx^2 + d_i * dx^3
!                   dx = u - x_i
! where i is defined by:
!                   x_i <= x < x_i+1
!------------------------------------------------------------------------------
! input:
!   x        -  grid points (necessarily in raising order)  
!   y        -  value of function at grid points            
!   b, c, d  -  array with coefficients from                
!               SPLINE function (dim : N-1)                 
!   u        -  point, where we want to aproximate          
!               value of function                           
!
! output:
!   ISPLINE  -  value at u position
!==============================================================================
    	implicit none

!      ----------------------------INPUT---------------------------------
        integer :: N
        double precision, dimension(N) :: x, y
        double precision, dimension(N) :: b, c, d
        double precision :: u, dx

!      -----------------------SUPPLEMENTARY------------------------------
        integer :: l, k, mid
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~BINARY SEARCH OF i~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF (u .ge. x(N)) THEN
            ISPLINE = y(N)
        ELSE IF (u .lt. x(1)) THEN
            ISPLINE = y(1)
        END IF

        l = 1
        k = n+1
        do while (k .gt. l+1)

            mid = (l + k) / 2

            IF ( x(mid) .gt. u ) THEN
                k = mid
            ELSE
                l = mid
            END IF

        end do
        dx = u - x(l)
        ISPLINE = y(l) + dx*(b(l) + dx*(c(l) + d(l) * dx))

    	END FUNCTION ISPLINE
!--------------------------------------------------------------------------------
! From special_functions libraryL
!--------------------------------------------------------------------------------
		subroutine rctj(n,x,nm,rj,dj)

!*****************************************************************************80
!
!! RCTJ computes Riccati-Bessel function of the first kind,and derivatives.
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
!    Input,integer(kind = 4) N,the order of jn(x).
!
!    Input,real(kind = 8) X,the argument.
!
!    Output,integer(kind = 4) NM,the highest order computed.
!
!    Output,real(kind = 8) RJ(0:N),the values of x jn(x).
!
!    Output,real(kind = 8) DJ(0:N),the values of [x jn(x)]'.
!
		  implicit none

		  integer(kind = 4) n

		  real(kind = 8) cs
		  real(kind = 8) dj(0:n)
		  real(kind = 8) f
		  real(kind = 8) f0
		  real(kind = 8) f1
		  integer(kind = 4) k
		  integer(kind = 4) m
!		  integer(kind = 4) msta1
!		  integer(kind = 4) msta2
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
		!  write(*,*) n,size(rj)
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
!--------------------------------------------------------------------
		subroutine rcty(n,x,nm,ry,dy)

!*****************************************************************************80
!
!! RCTY computes Riccati-Bessel function of the second kind,and derivatives.
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
!
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
!-------------------------------------------------------------------------------
		function envj(n,x)

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
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
!
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
!-------------------------------------------------------------------
		function msta1(x,mp)
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
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
!
		  implicit none

		  real(kind = 8) a0
!		  real(kind = 8) envj
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
!------------------------------------------------------------------
		function msta2(x,n,mp)
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
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
!
		  implicit none

		  real(kind = 8) a0
		  real(kind = 8) ejn
!		  real(kind = 8) envj
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
!-------------------------------------------------------------
!---------------------------------------------------------------------------
! Modified Bessel Function taken from Molscat
!==============================================================================
		SUBROUTINE MBESSK(V,X,CK,DK,EK)
		IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
!  -----------------------------------------------------------------
!  THIS SUBROUTINE USES TEMME'S METHOD [ N.M.TEMME,J COMPUT PHYS
!  19(1975) 324-337 ] TO CALCULATE THE MODIFIED BESSEL FUNCTION
!
!  K(V,X) = CK * EXP(EK)
!
!  AND ITS FIRST DERIVATIVE WITH RESPECT TO X
!
!  D/DX K(V,X) = DK * EXP(EK)
!
!  FOR A GIVEN REAL ORDER V >= 0 AND REAL ARGUMENT X > 0.
!  NOTE THE EXPONENTIAL SCALING,WHICH IS USED TO AVOID
!  OVERFLOW OF K(V,X) FOR V >> X AND UNDERFLOW FOR V << X.
!  -----------------------------------------------------------------
!
		PARAMETER(EPS = 1.D-15)! CONSISTENT WITH RGAMMA
		PARAMETER(MAXIT = 1000)
!
		IF(V.LT.0.D0 .OR. X.LE.0.D0) STOP 'MBESSK 0'
		PI = ACOS(-1.D0)
		XMIN = 1.D0
!
!  BEGIN BY CALCULATING K(A,X) AND K(A+1,X) FOR |A| <= 1/2
!
		NA = INT(V+0.5D0)
		A = V-NA
		IF(X.LT.XMIN) THEN
!
!  USING TEMME'S SERIES FOR SMALL X
!
		   B = X/2.D0
		   D = -DLOG(B)
		   E = A*D
		   C = A*PI
		   IF(ABS(C).LT.EPS) THEN
		      C = 1.D0
		   ELSE
		      C = C/SIN(C)
		   ENDIF
		   IF(ABS(E).LT.EPS) THEN
		      S = 1.D0
		   ELSE
		      S = SINH(E)/E
		   ENDIF
		   E = EXP(E)
		   G = E*RGAMMA(A,P,Q)
		   E =(E+1.D0/E)/2.D0
		   F = C*(P*E+Q*S*D)
		   E = A*A
		   P = 0.5D0*G*C
		   Q = 0.5D0/G
		   C = 1.D0
		   D = B*B
		   AK = F
		   AK1 = P
		   DO N = 1,MAXIT
		      F =(F*N+P+Q)/(N*N-E)
		      C = C*D/N
		      P = P/(N-A)
		      Q = Q/(N+A)
		      G = C*(P-N*F)
		      H = C*F
		      AK = AK+H
		      AK1 = AK1+G
		      IF(H/AK+ABS(G)/AK1.LT.EPS) GOTO 1
		   ENDDO
		   STOP 'MBESSK 1'
	1     F = AK
		   G = AK1/B
		   EX = 0.D0
		ELSEIF(X.GE.XMIN) THEN
!
!  AND TEMME'S PQ METHOD FOR LARGE X
!
		   C = 0.25D0-A*A
		   G = 1.D0
		   F = 0.D0
		   E = X*COS(A*PI)/PI/EPS
		   DO N = 1,MAXIT
		      H =(2*(N+X)*G-(N-1+C/N)*F)/(N+1)
		      F = G
		      G = H
		      IF(H*N.GT.E) GOTO 2
		   ENDDO
		   STOP 'MBESSK 2'
	2     P = F/G
		   Q = P
		   B = X+X
		   E = B-2.D0
		   DO M = N,1,-1
		      P =(M-1+C/M)/(E+(M+1)*(2.D0-P))
		      Q = P*(Q+1.D0)
		   ENDDO
		   F = SQRT(PI/B)/(1.D0+Q)
		   G = F*(A+X+0.5D0-P)/X
		   EX = X
		ENDIF
!
!  NOW RECUR UPWARDS FROM K(A,X) TO K(V,X),
!  SCALING TO AVOID OVERFLOW ALONG THE WAY
!
		P = 0.D0
		IF(NA.GT.0) THEN
		   Y = 2.D0/X
		   DO N = 1,NA
		      H = Y*(A+N)*G+F
		      F = G
		      G = H
	3        IF(ABS(F).GT.4.D0) THEN
		         P = P+1.D0
		         F = 0.0625D0*F
		         G = 0.0625D0*G
		         GOTO 3
		      ENDIF
		   ENDDO
		ENDIF
		CK = F
		DK =(V/X)*F-G
		SK = SQRT(CK*CK+DK*DK)
		CK = CK/SK
		DK = DK/SK
		EK = DLOG(SK)+P*DLOG(16.D0)-EX
		RETURN
		END
!--------------------------------------------------
!==============================================================================
		FUNCTION RGAMMA(X,ODD,EVEN)
		IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
!  -----------------------------------------------------------------
!  DIRECT FORTRAN TRANSLATION OF TEMME'S ALGOL ROUTINE FOR COMPUTING
!  RGAMMA = 1/GAMMA(1-X),ALONG WITH ITS ODD AND EVEN PARTS,FOR
!  ABS(X) .LE. 0.5. [ N.M.TEMME,J COMPUT PHYS 19(1975) 324-337 ]
!  -----------------------------------------------------------------
!
		DIMENSION B(12)
		DATA B / -0.283876542276024D0,-0.076852840844786D0,&
		         +0.001706305071096D0,+0.001271927136655D0,&
		         +0.000076309597586D0,-0.000004971736704D0,&
		         -0.000000865920800D0,-0.000000033126120D0,&
		         +0.000000001745136D0,+0.000000000242310D0,&
		         +0.000000000009161D0,-0.000000000000170D0 /
		SAVE B
!
		X2 = X*X*8.D0
		ALFA = -0.000000000000001D0
		BETA = 0.D0
		DO I = 12,2,-2
		   BETA = -(2*ALFA+BETA)
		   ALFA = -BETA*X2-ALFA+B(I)
		ENDDO
		EVEN =(BETA/2.D0+ALFA)*X2-ALFA+0.921870293650453D0
		ALFA = -0.000000000000034D0
		BETA = 0.D0
		DO I = 11,1,-2
		   BETA = -(2*ALFA+BETA)
		   ALFA = -BETA*X2-ALFA+B(I)
		ENDDO
		ODD = 2*(ALFA+BETA)
		RGAMMA = ODD*X+EVEN
		RETURN
		END
end module supplementary
