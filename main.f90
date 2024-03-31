! ILLUSTRATING PROGRAM FOR SPLINE PROGRAMS
! AND SEVAL
!
   implicit none
   double precision x(10), y(10), b(10), c(10), d(10)
   double precision s, u, seval
   integer i, n

   n = 10
   do i = 1, n
      x(i) = i
      y(i) = x(i)**3
   end do

   call spline(n, x, y, b, c, d)

   u = 2.5
   s = seval(n, u, x, y, b, c, d)
   write (6, 2) u, s
2  format(2f10.5)
   stop
end

subroutine spline(n, x, y, b, c, d)
   implicit none
   integer n
   double precision x(n), y(n), b(n), c(n), d(n)
!
!  THE COEFFICIENTS B(I), C(I) AND D(I) ARE CALCULATED, 1=1,
!  2, ..., N, FOR CUBIC INTERPOLATION SPLINE
!
!  S(X) = Y(I)+B(I)*(X-X(I)) + C(I)*(X-X(I))**2 +
!  -fD(I)*(X - X(I))**3
!
!  FOR X(I) .LE. X .LE. X(I+1)
!
!  INPUT INFORMATION..
!
!  N = NUMBER OF SPECIFIED POINTS OR NODES (N .GE. 2)
!  X = ABSCISSUE OF NODES IN STRICTLY INCREASING ORDER
!  Y = ORDINATES OF NODES
!
!  OUTPUT...
!
!  B, C, D = ARRAYS OF SPLINE COEFFICIENTS DEFINITED ABOVE.
!
!  IF YOU DESIGNATE THE DIFFERENTIATION SYMBOL BY P, THEN
!
!  Y(I)= S(X(I))
!  B(I) = SP(X(I))
!  C(I) = SPP(X(I))/2
!  D(I) = SPPP(X(I))/6 (RIGHT HAND DERIVATIVE)
!
!  USING THE ACCOMPANYING SEVAL FUNCTION SUBROUTINE
!  YOU CAN CALCULATE SPLINE VALUES.
!
   integer nm1, ib, i
   real t

   nm1 = n - 1
   if (n .lt. 2) return
   if (n .lt. 3) go to 50
!
! BUILD A TRIDIAGONAL SYSTEM
! B = DIAGONAL, O = OVERDIAGONAL, C = RIGHT PARTS.
!
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i = 2, nm1
      d(i) = x(i + 1) - x(i)
      b(i) = 2.*(d(i - 1) + d(i))
      c(i + 1) = (y(i + 1) - y(i))/d(i)
      c(i) = c(i + 1) - c(i)
   end do
!
! BOUNDARY CONDITIONS. THIRD DERIVATIVES AT POINTS
! X(1) AND X(N) ARE CALCULATED USING DIVISIONED
! DIFFERENCES
!
   b(1) = -d(1)
   b(n) = -d(n - 1)
   c(1) = 0.
   c(n) = 0.
   if (n .eq. 3) go to 15
   c(1) = c(3)/(x(4) - x(2)) - c(2)/(x(3) - x(1))
   c(n) = c(n - 1)/(x(n) - x(n - 2)) - c(n - 2)/(x(n - 1) - x(n - 3))
   c(1) = c(1)*d(1)**2/(x(4) - x(1))
   c(n) = -c(n)*d(n - 1)**2/(x(n) - x(n - 3))
!
! STRAIGHT RUN
!
15 do i = 2, n
      t = d(i - 1)/b(i - 1)
      b(i) = b(i) - t*d(i - 1)
      c(i) = c(i) - t*c(i - 1)
   end do
!
! REVERSE SUBSTITUSTION
!
   c(n) = c(n)/b(n)
   do ib = 1, nm1
      i = n - ib
      c(i) = (c(i) - d(i)*c(i + 1))/b(i)
   end do
!
! C(I) NOW STORES THE VALUE OF SIGMA(I), DEFINED
! IN #4.4.
!
! CALCULATE COEFFICIENTS OF POLYNOMIALS
!
   b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
   do i = 1, nm1
      b(i) = (y(i + 1) - y(i))/d(i) - d(i)*(c(i + 1) + 2.*c(i))
      d(i) = (c(i + 1) - c(i))/d(i)
      c(i) = 3.*c(i)
   end do
   c(n) = 3.*c(n)
   d(n) = d(n - 1)
   return

50 b(1) = (y(2) - y(1))/(x(2) - x(1))
   c(1) = 0.
   d(1) = 0.
   b(2) = b(1)
   c(2) = 0.
   d(2) = 0.
   return
end subroutine spline

double precision function seval(n, u, x, y, b, c, d)
   implicit none
   integer n
   double precision u, x(n), y(n), b(n), c(n), d(n)
!
!THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
!
!SEVAL = Y(I)+B(I)*(U-X(I)) + C(I)*(U-X(I)))**2 + D(I)*(U-X(I))**3
!
!WHERE X(I) .LT. U .LT. X(I + 1). USING HORNER'S RULE
!
!IF U .LT. X(1), THEN THE VALUE 1 = 1 IS TAKEN.
!IF U .GE. X(N), THEN THE VALUE I = N IS TAKEN.
!
!INPUT..
!
!N = THE NUMBER OF DATA POINTS
!U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
!X, Y = THE ARRAYS OF DATA ABSCISSAS AND ORD1NATES
!B, C, D = ARRAYS OF SPLINE COEFFICIENTS, COMPUTED BY SPLINE SUBROUTINE
!
!IF U IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
!BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
!
   integer i, j, k
   real dx
   data i/1/
   if (i .ge. n) i = 1
   if (u .lt. x(i)) go to 10
   if (u .le. x(i + 1)) go to 30

!
! BINARY SEARCH
!
10 i = 1
   j = n + 1
20 k = (i + j)/2
   if (u .lt. x(k)) j = k
   if (u .ge. x(k)) i = k
   if (j .gt. i + 1) go to 20
!
! EVALUATE SPLINE
!
30 dx = u - x(i)
   seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
   return
end function seval
