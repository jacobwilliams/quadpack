program test_qng
    use quadpack
implicit none

double precision a,abserr,b,f,epsabs,epsrel,result
integer ier,neval
external f
a = 0.0d0
b = 1.0d0
epsabs = 0.0d0
epsrel = 1.0d-3
call dqng(f,a,b,epsabs,epsrel,result,abserr,neval,ier)

! result from maxima: quad_qags(exp(x)/(x*x+1), x, 0, 1);
write(*,'(1P,A25,1X,*(E13.6,1X))') 'dqng: result, error = ', result, result - 1.27072413983362d0
end program test_qng

double precision function f(x)
double precision x
f = exp(x)/(x*x+0.1d+01)
return
end function f
