program test_qagp
implicit none

double precision a,abserr,b,epsabs,epsrel,f,points,result,work
integer ier,iwork,last,leniw,lenw,limit,neval,npts2
dimension iwork(204),points(4),work(404)
external f
a = 0.0e0
b = 1.0e0
npts2 = 4
points(1) = 1.0e0/7.0e0
points(2) = 2.0e0/3.0e0
limit = 100
leniw = limit*2+npts2
lenw = limit*4+npts2
call dqagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,&
          neval,ier,leniw,lenw,last,iwork,work)
!  include write statements
end program test_qagp

double precision function f(x)
double precision x
f = 0.0e+00
if(x/=1.0e0/7.0e0.and.x/=2.0e0/3.0e0) f = &
   abs(x-1.0e0/7.0e0)**(-0.25e0)* &
   abs(x-2.0e0/3.0e0)**(-0.55e0)
return
end function f
