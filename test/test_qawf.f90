program test_qawf
    use quadpack
implicit none

double precision a,abserr,epsabs,f,result,omega,work
integer ier,integr,iwork,last,leniw,lenw,limit,limlst,&
        lst,maxp1,neval
dimension iwork(250),work(1025)
external f

a = 0.0d0
omega = 8.0d0
integr = 2
epsabs = 1.0d-3
limlst = 50
limit = 100
leniw = limit*2+limlst
maxp1 = 21
lenw = leniw*2+maxp1*25
call dqawf(f,a,omega,integr,epsabs,result,abserr,neval,&
          ier,limlst,lst,leniw,maxp1,lenw,iwork,work)
write(*,*) 'result = ', result
end program test_qawf

double precision function f(x)
implicit none
double precision x
if(x>0.0d0) then
    f = sin(50.0d0*x)/(x*sqrt(x))
else
    f = 0.0d0
end if
end function f
