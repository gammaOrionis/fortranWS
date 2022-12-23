program fractWS
 use, intrinsic :: iso_fortran_env, only: wp => real128

!
 use mod_math, only: ws, fdo, setNodesAndWeights, getInitialGuessNodes, &
                     dataFinDiffF, deriv, quadrature, integrateSin, testExactSin, &
                     setParmsWS, integrateWS
 implicit none

 integer :: i

 integer, parameter :: imax = 146

 
 real(wp), dimension(imax) :: position
 

 real(wp)   , parameter :: step = 0.1_wp
 real(wp)   , parameter :: dummy = 1.0E-11_wp
 
 
 real(wp)   , parameter :: alpha = -0.95_wp
 real(wp)   , parameter :: W = 7_wp
 real(wp)   , parameter :: a0 = 0.7_wp
 real(wp)   , parameter :: R0 = 6.8_wp
 real(wp)   , parameter :: k = 1.0_wp
 real(wp)   , parameter :: m = 2.0_wp
 

 integer    , parameter :: n = 50
 real(wp)   , parameter :: x = 4.9_wp
 
 type (quadrature) :: qGL

 real(wp) :: result(imax)
 real(wp) :: exact(imax)

 
 
! generate values
 do concurrent (i = 1 : imax)
  position(i) = i*step
 end do  
    

!  call fdo fractional damped oscillation
!y = fdo( alpha, W, a0, k, m, position )

!open(1, file = 'data1.txt', status = 'old')  


!print *, getInitialGuessNodes(alpha)

 120 format (1x,2(1x,E44.33)) 
 130 format (1x,3(1x,E44.33)) 
 

! qGL contains alpha, n, nodes, weights
qGL =  setNodesAndWeights( alpha )
result = integrateSin( position, qGL )
exact = testExactSin( position , qGL)
  
do i = 1, imax
    print 130, position(i), result(i), result(i)-exact(i)
    ! write(1,100) position(i), y(i)
end do  

! now do the job for Woods-Saxon
result = integrateWS( position, qGL, setParmsWS(W, R0, a0) )
do i = 1, imax
    print 120, position(i), result(i)
    ! write(1,100) position(i), y(i)
end do  


!close(1) 
 
!print *,  sum(dataFinDiffF)
!print *,  deriv(1.0_wp)  - exp(1.0_wp)
end program fractWS