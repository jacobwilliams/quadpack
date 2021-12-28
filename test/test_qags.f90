program test_qags
use quadpack
implicit none

double precision a,abserr,b,epsabs,epsrel,result,work
integer ier,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)

a = 0.0d0
b = 1.0d0
epsabs = 0.0d0
epsrel = 1.0d-9
limit = 100
lenw = limit*4
call dqags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,&
          limit,lenw,last,iwork,work)

write(*,'(1P,A25,1X,2(E13.6,1X),I6)') 'dqags: result, error = ', result, result - 2.0d0, neval

contains

double precision function f(x)
implicit none
double precision,intent(in) :: x
f = 0.0d0
if(x>0.0d0) f = 1.0d0/sqrt(x)
end function f

end program test_qags
