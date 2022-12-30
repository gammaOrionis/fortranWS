program fractWS
 use, intrinsic :: iso_fortran_env, only: wp => real128

!
 use mod_math, only: ws, fdo, setNodesAndWeights, &
                     quadrature, setParmsWS, integrateWS, derive50
 implicit none

 integer :: i

 real(wp),  parameter       :: stepSize = 0.1_wp
 integer ,  parameter       :: imax = 150
 real(wp),  dimension(imax) :: positions
 
 integer ,  parameter       :: dmax = 75
 real(wp)                   :: res50(dmax)
  
 ! parameters of the fractional Woods-Saxon potential 
 real(wp)   , parameter :: alpha = -0.95_wp
 real(wp)   , parameter :: W = 7_wp
 real(wp)   , parameter :: a0 = 0.7_wp
 real(wp)   , parameter :: R0 = 6.8_wp
 ! additional parameters of the fractional damped oscillations 
 real(wp)   , parameter :: k = 1.0_wp
 real(wp)   , parameter :: m = 2.0_wp

 
 type (quadrature)      :: qGLAlpha

 real(wp)               :: result(imax)
 
 
! generate values
 do concurrent (i = 1 : imax)
  positions(i) = i*stepSize
 end do  
    

!  call fdo fractional damped oscillation
!y = fdo( alpha, W, a0, k, m, positions )

!open(1, file = 'testResults.txt', status = 'old')  


110 format (1x,1(1x,E44.33)) 
120 format (1x,2(1x,E44.33)) 
 

! qGLAlpha contains alpha, n, nodes, weights
qGLAlpha    =  setNodesAndWeights( alpha )

! first do the fractional integral  for Woods-Saxon
result = integrateWS( positions, qGLAlpha, setParmsWS(W, R0, a0) )
print 120, positions, result
 
! first ordinary derivative of result
! dmax < imax - 50
res50 = derive50( stepSize, dmax, result )
print 110,   res50

!close(1) 
 
end program fractWS