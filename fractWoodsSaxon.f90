program fractWS
 use, intrinsic :: iso_fortran_env, only: wp => real128

!
 use mod_math, only: ws, fdo, setNodesAndWeights, getInitialGuessNodes, &
                     dataFinDiffF, deriv, quadrature, integrate
 implicit none

 integer :: i

 integer, parameter :: imax = 146

 
 real(wp), dimension(imax) :: position
 

 real(wp)   , parameter :: step = 0.1_wp
 real(wp)   , parameter :: dummy = 1.0E-11_wp
 
 
 real(wp)   , parameter :: alpha = -0.95_wp
 real(wp)   , parameter :: W = 7_wp
 real(wp)   , parameter :: a0 = 0.7_wp
 real(wp)   , parameter :: k = 1.0_wp
 real(wp)   , parameter :: m = 2.0_wp
 

 integer    , parameter :: n = 50
 real(wp)   , parameter :: x = 4.9_wp
 
 type (quadrature) :: qGL

 real(wp) :: result(2)

 
 
! generate values
 do concurrent (i = 1 : imax)
  position(i) = i*step
 end do  
    
 !resultn = ws( rho0n,  R0n,  a0n, position )
 !resultp = ws( rho0p,  R0p,  a0p, position )
   
! write results

!do i = 1, imax
!end do   
!1000 format (1x,f5.2, 6(1x,E14.8), 2x, 4(1x,E14.8))

!  call fdo fractional damped oscillation
!y = fdo( alpha, W, a0, k, m, position )

!open(1, file = 'data1.txt', status = 'old')  

do i = 1, imax
 ! print 100, position(i), y(i)
 ! write(1,100) position(i), y(i)
end do  
!100 format (1x,2(1x,E60.33)) 
!print *, getInitialGuessNodes(alpha)

! qGL contains alpha, n, nodes, weights
qGL =  setNodesAndWeights( alpha )
result = integrate([1.0_wp,1.0_wp], qGL)
  
print *, qGL
print *, sum(qGL%nodes* qGL%weights)
print *, result

!close(1) 
 
print *,  sum(dataFinDiffF)
print *,  deriv(1.0_wp)  - exp(1.0_wp)
end program fractWS