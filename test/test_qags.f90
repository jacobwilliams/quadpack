program test_qags
    use quadpack
implicit none

double precision a,abserr,b,epsabs,epsrel,f,result,work
integer ier,iwork,last,lenw,limit,neval
dimension iwork(100),work(400)
external f
a = 0.0e0
b = 1.0e0
epsabs = 0.0e0
epsrel = 1.0e-3
limit = 100
lenw = limit*4
call dqags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,&
          limit,lenw,last,iwork,work)
write(*,*) 'result = ', result
end program test_qags

double precision function f(x)
implicit none
double precision x
f = 0.0e0
if(x>0.0e0) f = 1.0e0/sqrt(x)
return
end function f
