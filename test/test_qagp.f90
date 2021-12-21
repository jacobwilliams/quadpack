program test_qagp
    use quadpack
    implicit none

    integer,parameter :: npts2 = 4
    integer,parameter :: limit = 100
    integer,parameter :: leniw = limit*2 + npts2
    integer,parameter :: lenw = limit*4 + npts2
    double precision,parameter :: epsabs = 0.0d0
    double precision,parameter :: epsrel = 1.0d-9

    double precision a,abserr,b,points,result,work
    integer ier,iwork,last,neval
    dimension iwork(leniw),points(npts2),work(lenw)

    a = 0.0d0
    b = 1.0d0
    points(1) = 1.0d0/7.0d0
    points(2) = 2.0d0/3.0d0
    call dqagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,&
               neval,ier,leniw,lenw,last,iwork,work)

    ! answer from maxima: quad_qags(abs(x-1/7)^(-0.25)*abs(x-2/3)^(-0.55), x, 0, 1);
    write(*,'(1P,A25,1X,*(E13.6,1X))') 'dqagp: result, error = ', result, result - 4.253687688108305D0

contains

    double precision function f(x)
    implicit none
    double precision,intent(in) :: x
    f = 0.0d+00
    if(x/=1.0d0/7.0d0.and.x/=2.0d0/3.0d0) f = &
        abs(x-1.0d0/7.0d0)**(-0.25d0)* &
        abs(x-2.0d0/3.0d0)**(-0.55d0)
    return
    end function f

end program test_qagp
