#ifndef MOD_INCLUDE
    module quadpack_test_module
    use quadpack, wp => quadpack_RK
#endif

    implicit none

    public

    real(wp), parameter :: epsabs = 0.0_wp
    real(wp), parameter :: epsrel = 10**( log10(epsilon(1.0_wp)) / 2.0_wp )

contains

    subroutine check_result(routine, value, answer, neval)
        implicit none
        character(len=*),intent(in) :: routine
        real(wp),intent(in) :: value
        real(wp),intent(in) :: answer
        integer,intent(in) :: neval

        write (*, '(1P,A25,1X,2(E13.6,1X),I6)') routine, value, abs(value - answer), neval

        if (abs(value-answer)>epsrel) error stop 'TEST FAILED'

    end subroutine check_result

    subroutine test_qag
        implicit none

        real(wp) :: a, abserr, b, result, work
        integer :: ier, iwork, key, last, lenw, limit, neval
        dimension iwork(100), work(400)

        real(wp),parameter :: answer = 2.0_wp/sqrt(3.0_wp)

        a = 0.0_wp
        b = 1.0_wp
        key = 6
        limit = 100
        lenw = limit*4
        call dqag(f, a, b, epsabs, epsrel, key, result, abserr, neval, &
                  ier, limit, lenw, last, iwork, work)

        call check_result('dqag', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            real(wp), parameter :: pi = acos(-1.0_wp)
            f = 2.0_wp/(2.0_wp + sin(10.0_wp*pi*x))
        end function f

    end subroutine test_qag

    subroutine test_qagi
        implicit none

        real(wp) abserr, boun, result, work
        integer ier, inf, iwork, last, lenw, limit, neval
        dimension iwork(100), work(400)

        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp),parameter :: answer = sqrt(2.0_wp)*pi*log(2.0_wp)

        boun = 0.0_wp
        inf = 1
        limit = 100
        lenw = limit*4
        call dqagi(f, boun, inf, epsabs, epsrel, result, abserr, neval, &
                   ier, limit, lenw, last, iwork, work)

        call check_result('dqagi', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            f = 0.0_wp
            if (x > 0.0_wp) f = sqrt(x)*log(x)/((x + 1.0_wp)*(x + 2.0_wp))
        end function f

    end subroutine test_qagi

    subroutine test_qagp
        implicit none

        integer, parameter :: npts2 = 4
        integer, parameter :: limit = 100
        integer, parameter :: leniw = limit*2 + npts2
        integer, parameter :: lenw = limit*4 + npts2

        real(wp),parameter :: answer = 4.253687688108305_wp

        real(wp) a, abserr, b, points, result, work
        integer ier, iwork, last, neval
        dimension iwork(leniw), points(npts2), work(lenw)

        a = 0.0_wp
        b = 1.0_wp
        points(1) = 1.0_wp/7.0_wp
        points(2) = 2.0_wp/3.0_wp
        call dqagp(f, a, b, npts2, points, epsabs, epsrel, result, abserr, &
                   neval, ier, leniw, lenw, last, iwork, work)

        ! answer from maxima: quad_qags(abs(x-1/7)^(-0.25)*abs(x-2/3)^(-0.55), x, 0, 1);
        call check_result('dqagp', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            f = 0.0d+00
            if (x /= 1.0_wp/7.0_wp .and. x /= 2.0_wp/3.0_wp) f = &
                abs(x - 1.0_wp/7.0_wp)**(-0.25_wp)* &
                abs(x - 2.0_wp/3.0_wp)**(-0.55_wp)
            return
        end function f

    end subroutine test_qagp

    subroutine test_qags
        implicit none

        real(wp) a, abserr, b, result, work
        integer ier, iwork, last, lenw, limit, neval
        dimension iwork(100), work(400)

        real(wp),parameter :: answer = 2.0_wp

        a = 0.0_wp
        b = 1.0_wp
        limit = 100
        lenw = limit*4
        call dqags(f, a, b, epsabs, epsrel, result, abserr, neval, ier, &
                   limit, lenw, last, iwork, work)

        call check_result('dqags', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            f = 0.0_wp
            if (x > 0.0_wp) f = 1.0_wp/sqrt(x)
        end function f

    end subroutine test_qags

    subroutine test_qawc
        implicit none

        real(wp) a, abserr, b, c, result, work
        integer ier, iwork, last, lenw, limit, neval
        dimension iwork(100), work(400)

        real(wp),parameter :: answer = -628.4617285065623_wp

        a = -1.0_wp
        b = 1.0_wp
        c = 0.5_wp
        limit = 100
        lenw = limit*4
        call dqawc(f, a, b, c, epsabs, epsrel, result, abserr, neval, &
                   ier, limit, lenw, last, iwork, work)

        ! maxima: quad_qawc((1/(x*x+1.0e-4)), x, 0.5, -1, 1);
        call check_result('dqawc', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            f = 1.0_wp/(x*x + 1.0d-4)
        end function f

    end subroutine test_qawc

    subroutine test_qawf
        implicit none

        real(wp) a, abserr, result, omega, work
        integer ier, integr, iwork, last, leniw, lenw, limit, limlst, &
            lst, maxp1, neval
        dimension iwork(250), work(1025)

        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp),parameter :: answer = sqrt(29.0_wp*pi) - sqrt(21.0_wp*pi)

        a = 0.0_wp
        omega = 8.0_wp
        integr = 2
        limlst = 50
        limit = 100
        leniw = limit*2 + limlst
        maxp1 = 21
        lenw = leniw*2 + maxp1*25
        call dqawf(f, a, omega, integr, epsrel, result, abserr, neval, &
                   ier, limlst, lst, leniw, maxp1, lenw, iwork, work)

        ! wolfram alpha: integral of  sin(8*x)*sin(50*x)/(x*sqrt(x)) from 0 to infinity
        call check_result('dqawf', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            if (x > 0.0_wp) then
                f = sin(50.0_wp*x)/(x*sqrt(x))
            else
                f = 0.0_wp
            end if
        end function f

    end subroutine test_qawf

    subroutine test_qawo
        implicit none

        real(wp) a, abserr, b, result, omega, work
        integer ier, integr, iwork, last, leniw, lenw, limit, maxp1, neval
        dimension iwork(200), work(925)

        real(wp),parameter :: answer = -0.17763920651138_wp

        a = 0.0_wp
        b = 1.0_wp
        omega = 10.0_wp
        integr = 1
        limit = 100
        leniw = limit*2
        maxp1 = 21
        lenw = limit*4 + maxp1*25

        call dqawo(f, a, b, omega, integr, epsabs, epsrel, result, abserr, &
                   neval, ier, leniw, maxp1, lenw, last, iwork, work)

        ! result from maxima: quad_qags(exp(-x)*log(x)*cos(10*x), x, 0, 1);
        call check_result('dqawo', result, answer, neval)

    contains

        real(wp) function f(x)
            real(wp), intent(in) :: x
            f = 0.0_wp
            if (x > 0.0_wp) f = exp(-x)*log(x)
        end function f

    end subroutine test_qawo

    subroutine test_qaws
        implicit none

        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        real(wp), parameter :: alfa = -0.5_wp
        real(wp), parameter :: beta = -0.5_wp
        integer, parameter :: integr = 1
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: answer = 0.535019056922365344126359_wp

        real(wp) :: result, abserr
        integer :: neval, ier, last
        integer, dimension(limit) :: iwork
        real(wp), dimension(lenw) :: work

        call dqaws(f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
                   abserr, neval, ier, limit, lenw, last, iwork, work)

        call check_result('dqaws', result, answer, neval)

    contains

        function f(x)
            implicit none
            real(wp) :: f
            real(wp), intent(in) :: x
            f = sin(10.0_wp*x)
        end function f

    end subroutine test_qaws

    subroutine test_qng
        implicit none

        real(wp) a, abserr, b, result
        integer ier, neval

        real(wp), parameter :: answer = 1.27072413983362_wp

        a = 0.0_wp
        b = 1.0_wp
        call dqng(f, a, b, epsabs, epsrel, result, abserr, neval, ier)

        ! result from maxima: quad_qags(exp(x)/(x*x+1), x, 0, 1);
        call check_result('dqng', result, answer, neval)

    contains

        real(wp) function f(x)
            real(wp), intent(in) :: x
            f = exp(x)/(x*x + 1.0_wp)
            return
        end function f

        end subroutine test_qng

#ifndef MOD_INCLUDE
    end module quadpack_test_module
#endif
