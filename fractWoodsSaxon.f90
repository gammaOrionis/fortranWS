program fractWS
 use, intrinsic :: iso_fortran_env, only: wp => real128

!
 use mod_math, only: ws, fdo, setNodesAndWeights, &
                     quadrature, integrateSin, testExactSin, &
                     setParmsWS, integrateWS, derive50
 implicit none

 integer :: i

 integer, parameter :: imax = 150
 integer, parameter :: dmax = 70
 
 real(wp), dimension(imax) :: positions

 real(wp)   , parameter :: step = 0.1_wp
 
 ! parameters of the fractional Woods-Saxon potential 
 real(wp)   , parameter :: alpha = -0.95_wp
 real(wp)   , parameter :: W = 7_wp
 real(wp)   , parameter :: a0 = 0.7_wp
 real(wp)   , parameter :: R0 = 6.8_wp
 real(wp)   , parameter :: k = 1.0_wp
 real(wp)   , parameter :: m = 2.0_wp

 integer    , parameter :: n = 50
 real(wp)   , parameter :: x = 4.9_wp
 
 type (quadrature) :: qGLAlpha

 real(wp) :: result(imax)
 real(wp) :: exact(imax)
 
 real(wp) :: res10(dmax)
 
! generate values
 do concurrent (i = 1 : imax)
  positions(i) = i*step
 end do  
    

!  call fdo fractional damped oscillation
!y = fdo( alpha, W, a0, k, m, positions )

!open(1, file = 'testResults.txt', status = 'old')  


110 format (1x,1(1x,E44.33)) 
120 format (1x,2(1x,E44.33)) 
130 format (1x,3(1x,E44.33)) 
 

! qGLAlpha contains alpha, n, nodes, weights
qGLAlpha    =  setNodesAndWeights( alpha )

! test tempered sine
!result = integrateSin( positions , qGLAlpha )
!exact  = testExactSin( positions , qGLAlpha)
!print 130, positions, result, result - exact
! write(1,130) positions(i), y(i)

! now do the job for Woods-Saxon
result = integrateWS( positions, qGLAlpha, setParmsWS(W, R0, a0) )
print 120, positions, result

! 
res10 = derive50( step, dmax, result )
print 110,   res10

!close(1) 
 
end program fractWS