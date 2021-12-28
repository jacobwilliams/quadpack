program test_qagi
use quadpack
implicit none

double precision abserr,boun,epsabs,epsrel,result,work
integer ier,inf,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)

double precision,parameter :: pi = acos(-1.0d0)

boun = 0.0d0
inf = 1
epsabs = 0.0d0
epsrel = 1.0d-9
limit = 100
lenw = limit*4
call dqagi(f,boun,inf,epsabs,epsrel,result,abserr,neval,&
          ier,limit,lenw,last,iwork,work)

write(*,'(1P,A25,1X,2(E13.6,1X),I6)') 'dqagi: result, error = ', result, result - sqrt(2.0d0)*pi*log(2.0d0), neval

contains

double precision function f(x)
implicit none
double precision,intent(in) :: x
f = 0.0d0
if(x>0.0d0) f = sqrt(x)*log(x)/((x+1.0d0)*(x+2.0d0))
end function f

end program test_qagi