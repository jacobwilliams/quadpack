program test_qawc
use quadpack
implicit none

double precision a,abserr,b,c,epsabs,epsrel,result,work
integer ier,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)

a = -1.0d0
b = 1.0d0
c = 0.5d0
epsabs = 0.0d0
epsrel = 1.0d-9
limit = 100
lenw = limit*4
call dqawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,&
           ier,limit,lenw,last,iwork,work)

! maxima: quad_qawc((1/(x*x+1.0e-4)), x, 0.5, -1, 1);
write(*,'(1P,A25,1X,*(E13.6,1X))') 'dqawc: result, error = ', result, result - (-628.4617285065623d0)

contains

double precision function f(x)
implicit none
double precision,intent(in) :: x
f = 1.0d0/(x*x+1.0d-4)
end function f

end program test_qawc
