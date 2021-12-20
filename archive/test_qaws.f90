program test_qaws
implicit none

! double precision a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
! integer ier,integr,iwork,last,lenw,limit,neval
! dimension iwork(100),work(400)
! external f
! a = 0.0e0
! b = 1.0e0
! alfa = -0.5e0
! beta = -0.5e0
! integr = 1
! epsabs = 0.0e0
! epsrel = 1.0e-3
! limit = 100
! lenw = limit*4
! call dqaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,&
!             abserr,neval,ier,limit,lenw,last,iwork,work)
!  include write statements
end program test_qaws

double precision function f(x)
double precision x
f = sin(10.0e0*x)
return
end function f
