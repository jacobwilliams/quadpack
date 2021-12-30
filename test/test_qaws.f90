
    program test_qaws
    use quadpack
    implicit none

    real(wp), parameter :: a = 0.0_wp
    real(wp), parameter :: b = 1.0_wp
    real(wp), parameter :: alfa = -0.5_wp
    real(wp), parameter :: beta = -0.5_wp
    integer, parameter :: integr = 1
    real(wp), parameter :: epsabs = 0.0_wp
    real(wp), parameter :: epsrel = 1.0e-3_wp
    integer, parameter :: limit = 100
    integer, parameter :: lenw = limit*4
    real(wp), parameter :: ans = 0.535019056922365344126359_wp

    real(wp) :: result, abserr
    integer :: neval, ier, last
    integer, dimension(limit) :: iwork
    real(wp), dimension(lenw) :: work

    call dqaws(f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
               abserr, neval, ier, limit, lenw, last, iwork, work)

    write(*,'(1P,A25,1X,2(E13.6,1X),I6)') 'dqaws: result, error = ', result, result - ans, neval

    contains

    function f(x)
        implicit none
        real(wp) :: f
        real(wp),intent(in) :: x
        f = sin(10.0_wp*x)
    end function f

    end program test_qaws