program test_qawo
use quadpack
implicit none

double precision a,abserr,b,epsabs,epsrel,result,omega,work
integer ier,integr,iwork,last,leniw,lenw,limit,maxp1,neval
dimension iwork(200),work(925)

a = 0.0d0
b = 1.0d0
omega = 10.0d0
integr = 1
epsabs = 0.0d0
epsrel = 1.0d-9
limit = 100
leniw = limit*2
maxp1 = 21
lenw = limit*4+maxp1*25

call dqawo(f,a,b,omega,integr,epsabs,epsrel,result,abserr,&
           neval,ier,leniw,maxp1,lenw,last,iwork,work)

! result from maxima: quad_qags(exp(-x)*log(x)*cos(10*x), x, 0, 1);
write(*,'(1P,A25,1X,*(E13.6,1X))') 'dqawo: result, error = ', result, result - (-0.17763920651138d0)

contains

double precision function f(x)
double precision,intent(in) :: x
f = 0.0d0
if(x>0.0d0) f = exp(-x)*log(x)
end function f

end program test_qawo
