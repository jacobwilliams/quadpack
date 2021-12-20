program test_qawc
    use quadpack
implicit none

double precision a,abserr,b,c,epsabs,epsrel,f,result,work
integer ier,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)
external f
a = -1.0e0
b = 1.0e0
c = 0.5e0
epsabs = 0.0e0
epsrel = 1.0e-3
limit = 100
lenw = limit*4
call dqawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
write(*,*) 'result = ', result
end program test_qawc

double precision function f(x)
implicit none
double precision x
f = 1.0e0/(x*x+1.0e-4)
return
end function f