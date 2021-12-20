program test_qag
    use quadpack
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

write(*,*) 'result, error = ', result, result - 2.0d0/sqrt(3.0d0)

end program test_qag

double precision function f(x)
implicit none
double precision x
double precision,parameter :: pi = acos(-1.0d0)
f = 2.0d0/(2.0d0+sin(10.0d0 * pi * x))
return
end function f
