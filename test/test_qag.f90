program test_qag
implicit none

double precision a,abserr,b,epsabs,epsrel,f,result,work
integer ier,iwork,key,last,lenw,limit,neval
dimension iwork(100),work(400)
external f
a = 0.0d0
b = 1.0d0
epsabs = 0.0d0
epsrel = 1.0d-3
key = 6
limit = 100
lenw = limit*4
call dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,&
         ier,limit,lenw,last,iwork,work)
!  include write statements
end program test_qag

double precision function f(x)
double precision x
f = 2.0d0/(2.0d0+sin(31.41592653589793d0*x))
return
end function f
