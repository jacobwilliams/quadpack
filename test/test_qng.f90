program test_qng
    use quadpack
implicit none

double precision a,abserr,b,f,epsabs,epsrel,result
integer ier,neval
external f
a = 0.0e0
b = 1.0e0
epsabs = 0.0e0
epsrel = 1.0e-3
call dqng(f,a,b,epsabs,epsrel,result,abserr,neval,ier)
!  include write statements
end program test_qng

double precision function f(x)
double precision x
f = exp(x)/(x*x+0.1e+01)
return
end function f
