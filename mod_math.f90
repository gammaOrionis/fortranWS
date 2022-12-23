module mod_math

  use, intrinsic :: iso_fortran_env, only: wp => real128
  implicit none

  integer, parameter   :: nGL = 50        ! LaguerreL order n,  here fixed to 50
  integer, parameter   :: imax = 6        ! newton raphson iterations  4 -> error 10^-13 real64

  real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_wp

  type quadrature
    real(wp)                      :: alpha
    integer                       :: dimension
    real(wp), dimension(nGL)      :: nodes  
    real(wp), dimension(nGL)      :: weights
    real(wp), dimension(nGL)      :: expNodes
    
  end type 
  type paramsWS
    real(wp)                      :: rho0
    real(wp)                      :: R0
    real(wp)                      :: a0
  end type 

 
 
! Forward difference coefficients n = 50, precision up to 128
  real(wp), parameter ::  &
    dataFinDiffF(nGL) = (/    &
    -4.479205338329425057560471792964769091970601823967453829658902432_wp,&
    49._wp,&
    -588._wp,&
    6141.333333333333333333333333333333333333333333333333333333333333_wp,&
    -52969._wp,&
    381376.8_wp,&
    -2.330636e6_wp,&
    1.2271512e7_wp,&
    -5.637225825e7_wp,&
    2.282728482222222222222222222222222222222222222222222222222222222e8_wp,&
    -8.217822536e8_wp,&
    2.648719660363636363636363636363636363636363636363636363636363636e9_wp,&
    -7.688644569666666666666666666666666666666666666666666666666666667e9_wp,&
    2.019975259723076923076923076923076923076923076923076923076923077e10_wp,&
    -4.8232062324e10_wp,&
    1.050387135056e11_wp,&
    -2.092568120619375e11_wp,&
    3.823100234211176470588235294117647058823529411764705882352941176e11_wp,&
    -6.419032492008888888888888888888888888888888888888888888888888889e11_wp,&
    9.921939419781052631578947368421052631578947368421052631578947368e11_wp,&
    -1.4138763673188e12_wp,&
    1.859519938877333333333333333333333333333333333333333333333333333e12_wp,&
    -2.259086206735272727272727272727272727272727272727272727272727273e12_wp,&
    2.536667687714086956521739130434782608695652173913043478260869565e12_wp,&
    -2.6335543007865e12_wp,&
    2.52821212875504e12_wp,&
    -2.243975262208615384615384615384615384615384615384615384615384615e12_wp,&
    1.840736909191703703703703703703703703703703703703703703703703704e12_wp,&
    -1.394639954158e12_wp,&
    9.750871498750344827586206896551724137931034482758620689655172414e11_wp,&
    -6.283894965861333333333333333333333333333333333333333333333333333e11_wp,&
    3.727180156650322580645161290322580645161290322580645161290322581e11_wp,&
    -2.0310219994246875e11_wp,&
    1.014578482724545454545454545454545454545454545454545454545454545e11_wp,&
    -4.634060889952941176470588235294117647058823529411764705882352941e10_wp,&
    1.92928249296e10_wp,&
    -7.294355104555555555555555555555555555555555555555555555555555556e9_wp,&
    2.493614455027027027027027027027027027027027027027027027027027027e9_wp,&
    -7.667346385263157894736842105263157894736842105263157894736842105e8_wp,&
    2.107133983589743589743589743589743589743589743589743589743589744e8_wp,&
    -5.136139085e7_wp,&
    1.09994650243902439024390243902439024390243902439024390243902439e7_wp,&
    -2.045252e6_wp,&
    325205.0232558139534883720930232558139534883720930232558139534884_wp,&
    -43338.27272727272727272727272727272727272727272727272727272727273_wp,&
    4708.355555555555555555555555555555555555555555555555555555555556_wp,&
    -400.5217391304347826086956521739130434782608695652173913043478261_wp,&
    25.02127659574468085106382978723404255319148936170212765957446809_wp,&
    -1.020833333333333333333333333333333333333333333333333333333333333_wp,&
    0.02040816326530612244897959183673469387755102040816326530612244898_wp  /)



  
  
  
  
  
    ! first guess zeros of LaguerreL
  
    real(wp), parameter ::  &
    dataAlphaM1(nGL) = (/    &
    0.0, 0.0734188, 0.246193, 0.517944, 0.888921, 1.35949, 1.93012, &
    2.60140, 3.37399, 4.24872, 5.22648, 6.30832, 7.49541, 8.78904, &
    10.1907, 11.7019, 13.3245, 15.0603, 16.9117, 18.8808, 20.9703, &
    23.1830, 25.5219, 27.9906, 30.5928, 33.3326, 36.2145, 39.2436, &
    42.4256, 45.7666, 49.2736, 52.9545, 56.8180, 60.8742, 65.1346, &
    69.6121, 74.3220, 79.2821, 84.5134, 90.0411, 95.8958, 102.116, &
    108.748, 115.857, 123.525, 131.872, 141.077, 151.435, 163.512, 178.767 /)

  real(wp), parameter ::  &
    dataAlpha0(nGL) = (/    &
    0.0286305, 0.150883, 0.370949, 0.689091, 1.10563, 1.62096, 2.23561, &
    2.95018, 3.76540, 4.68209, 5.70120, 6.82379, 8.05106, 9.38435, &
    10.8251, 12.3750, 14.0358, 15.8094, 17.6981, 19.7041, 21.8302, &
    24.0792, 26.4541, 28.9584, 31.5959, 34.3707, 37.2875, 40.3513, &
    43.5677, 46.9430, 50.4843, 54.1992, 58.0968, 62.1871, 66.4814, &
    70.9929, 75.7370, 80.7314, 85.9972, 91.5597, 97.4496, 103.705, &
    110.374, 117.519, 125.225, 133.612, 142.858, 153.260, 165.386, 180.698  /)
  
  real(wp), parameter ::  &
    dataAlphaP1(nGL) = (/    &
    0.0719789, 0.241362, 0.507772, 0.871441, 1.33272, 1.89204, 2.54995, &
    3.30711, 4.16425, 5.12223, 6.18204, 7.34477, 8.61163, 9.98398, &
    11.4633, 13.0513, 14.7498, 16.5606, 18.4862, 20.5287, 22.6909, &
    24.9756, 27.3860, 29.9255, 32.5980, 35.4075, 38.3588, 41.4568, &
    44.7074, 48.1166, 51.6917, 55.4404, 59.3717, 63.4956, 67.8236, &
    72.3689, 77.1468, 82.1751, 87.4751, 93.0721, 98.9968, 105.287, &
    111.992, 119.174, 126.918, 135.344, 144.631, 155.077, 167.251, 182.620 /)
  
  contains
  
    pure function setNodesAndWeights(alpha) result(qGL)
    ! purpose : set nodes and weights 
      
    real(wp), intent(in) :: alpha           ! validity range -1 < alpha < +1   , no error checking 
    integer, parameter   :: n = nGL         ! LaguerreL order n here fixed to 50
     
    real(wp)             :: nodes(n)        ! nodes = zeroes of L[n, alpha, x_i]
    real(wp)             :: weights(n)      ! nodes = zeroes of L[n, alpha, x_i]
   
    integer              :: i               ! index 
   
    type ( quadrature ) :: qGL

    if(qGL%alpha .eq. alpha) return
   
    !! linear interpolation guess
    nodes = getInitialGuessNodes( alpha );
    
    !! newton raphson iteration to get node-accuracy
    do i = 1, imax
      nodes = iterateNodes(n, alpha, nodes)
    end do
    
    !! get weights
    weights = getWeights(n, alpha, nodes)

    !! finally get expNodes and finish
    
    qGL = quadrature( alpha, n, nodes, weights, exp(nodes) )
    
  end function setNodesAndWeights


  pure function getInitialGuessNodes( alpha ) result( x )
    ! purpose : approximate nodes for L[n, alpha, x_i] === 0
    
    real(wp), intent(in)   :: alpha                     ! validity range -1 < alpha < +1    
    real(wp)               :: x(size(dataAlpha0))       ! nodes = zeroes of L[nGL, alpha, x_i]
                                                        ! linear interpolation  
    if(alpha <= 0) then
      x = (1.0_wp + alpha) * dataAlpha0 - alpha * dataAlphaM1
    else
      x = (1.0_wp - alpha) * dataAlpha0 + alpha * dataAlphaP1
    end if

  end function getInitialGuessNodes



  pure  elemental  function iterateNodes( n, alpha, x ) result( x1 )
    
    integer,  intent(in) :: n
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: x

    integer  :: k
    real(wp) :: y0, y1, y2, x1

  
    ! start values
    y0 = 1.0_wp               ! L^alpha_0
    y1 = 1.0_wp + alpha - x   ! L^alpha_1

    do concurrent (k = 1:n-1)
      y2 = ((1 + k + k  + alpha -x)*y1 -(k + alpha)*y0)/(1 + k)
      y0 = y1
      y1 = y2
    end do 

    ! derivative D[L[n,alpha,x],x]
    y2 = -((n+alpha) * y0 - n*y1)/x

    ! next step in Newton
    x1 = x - y1/y2

  end function iterateNodes


  pure elemental function getWeights( n, alpha, x ) result( y1 )
    
    integer,  intent(in) :: n
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: x

    integer :: k
    real(wp) :: y0, y1, y2, fak


    ! start values
    
    y0 = 1.0_wp               ! L^alpha_0
    y1 = 1.0_wp + alpha - x   ! L^alpha_1

    ! loop up to L^alpha_{n+1}
    do concurrent (k = 1:n)
      y2 = ((1 + k + k + alpha -x)*y1 -(k + alpha)*y0)/(1 + k)
      y0 = y1
      y1 = y2
    end do 

    ! y2 = L[n+1,alpha,x] is known, now weights

      fak = 1.0_wp/(1 + n)
      fak = fak*fak*gamma(1 + n + alpha) / gamma(1.0_wp + n) 
      y1 = x*fak / (y2*y2)

  end function getWeights


  pure function deriv( x ) result( y )
    real(wp), intent(in) :: x
    
    real(wp) :: y
    integer  :: i
    real(wp), parameter  ::  h= 0.01_wp   ! stepsize 
    real(wp)  ::  res 

    res = 0.0_wp
    
    do i = 1, nGL
      res = res + dataFinDiffF(i)*exp( ( i-1 )*h + x)
    end do 
    res = res / h 
    
    y = res
    
  end function deriv

  ! 
  pure function setParmsWS( rho0, R0, a0 ) result(pWS)
    real(wp), intent(in) :: rho0 
    real(wp), intent(in) :: R0
    real(wp), intent(in) :: a0
    type ( paramsWS ) :: pWS
    
    pWS = paramsWS(rho0, R0, a0)
    
  end function setParmsWS

  pure elemental function ws( x, pWS ) result( y )
    real(wp), intent(in) :: x
    type ( paramsWS ) , intent(in) :: pWS
    
    real(wp) :: y
   
    y =  pWS%rho0/(1.0_wp + exp(( x - pWS%R0)/pWS%a0))
  end function ws

  pure function fdo( alpha, W, a0,  kreal, mreal,  x ) result( yreal )

    real(wp), intent(in) :: alpha, W, a0, kreal, mreal
    real(wp), intent(in) :: x(:)

    real(wp) :: yreal(size(x))

    complex(wp) :: k, m, factor
    complex(wp), parameter  :: im = (0.0_wp, 1.0_wp) 
    
    integer :: i
    
    i = size(x)
    k = kreal
    m = im*mreal
    factor = (k-m)**alpha

    yreal( 1:i ) =  a0**alpha * W * exp( -kreal * x(1:i) ) * real( factor * exp( m * x(1:i) ) )

  end function fdo

  pure elemental function integrateWS( x , qGL, pWS) result( y )
    type ( paramsWS )  , intent(in) :: pWS       ! parameters of the Woods-Saxon potential
    type ( quadrature ), intent(in) :: qGL       ! parameters of the Gauss-Lagurre quadrature
    real(wp), intent(in) :: x                    ! vector of input variable e.g. 0 .. 14.6
        
    real(wp) :: y                                ! vector of output
    
    
    ! y[a_, x_] =  NIntegrate[w[h] Exp[-h] Exp[+h] f[h+x],{h,0,Infinity}]
    ! with f[x_] = ws[x]
    
    y =  sum( qGL%weights *  ws( qGL%nodes + x , pWS) * qGL%expNodes)
  
  end function integrateWS



  pure function integrateSin( x , qGL) result( y )

    real(wp), intent(in) :: x(:)
    type ( quadrature ), intent(in) :: qGL
   
    real(wp) :: y(size(x))
  
    integer :: i

    ! y[a_, x_] =  NIntegrate[w[h] Exp[-h] f[h+x],{h,0,Infinity}]
    ! with f[x_] = Sin[x]
    
    do concurrent (i = 1 : size(x))   ! iterate all x-positions
      y(i) =  sum( qGL%weights * sin( qGL%nodes + x(i) ))
    end do   

  end function integrateSin

  pure function testExactSin( x , qGL) result( y )

  real(wp), intent(in) :: x(:)
  real(wp) :: factor
    type ( quadrature ), intent(in) :: qGL
  
    real(wp) :: y(size(x))
    
    integer :: i
    i = size(x)

    ! y[a_, x_] =  Integrate[w[h] Exp[-h] f[h+x],{h,0,Infinity}]
    factor = 2._wp**(0.5_wp*(-1.0_wp - qGL%alpha)) * gamma(1.0_wp + qGL%alpha)
  
    y(1:i) =   factor *  Cos(0.25_wp *(pi - qGL%alpha *  pi - 4* x(1:i))) 
  
  end function testExactSin
  
end module mod_math
