program test_qagi
implicit none

double precision abserr,boun,epsabs,epsrel,f,result,work
integer ier,inf,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)
external f
boun = 0.0e0
inf = 1
epsabs = 0.0e0
epsrel = 1.0e-3
limit = 100
lenw = limit*4
call dqagi(f,boun,inf,epsabs,epsrel,result,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
!  include write statements
end program test_qagi

double precision function f(x)
double precision x
f = 0.0e0
if(x>0.0e0) f = sqrt(x)*log(x)/((x+1.0e0)*(x+2.0e0))
return
end function f
