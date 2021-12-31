#ifndef MOD_INCLUDE
module quadpack_test_module
    use quadpack, wp => quadpack_RK
#endif

    implicit none

    public

    real(wp), parameter :: epsabs = 0.0_wp
    real(wp), parameter :: epsrel = 10**(log10(epsilon(1.0_wp))/2.0_wp+1)

contains

    subroutine check_result(routine, value, answer, neval)
        implicit none

        character(len=*), intent(in) :: routine
        real(wp), intent(in) :: value
        real(wp), intent(in) :: answer
        integer, intent(in) :: neval

        write (*, '(1P,A25,1X,2(E13.6,1X),I6)') routine, value, abs(value - answer), neval

        if (abs(value - answer) > epsrel) then
            write(*,*) '  value  = ', value
            write(*,*) '  answer = ', answer
            write(*,*) '  epsrel = ', epsrel
            write(*,*) 'TEST FAILED'
            !error stop 'TEST FAILED'
        end if

    end subroutine check_result

    subroutine test_qag
        implicit none

        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        integer, parameter :: key = 6
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: answer = 2.0_wp/sqrt(3.0_wp)

        real(wp) :: abserr, result, work(400)
        integer :: ier, iwork(100), last, neval

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

        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp), parameter :: boun = 0.0_wp
        integer, parameter :: inf = 1
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: answer = sqrt(2.0_wp)*pi*log(2.0_wp)

        real(wp) :: abserr, result, work(400)
        integer :: ier, iwork(100), last, neval

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
        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        real(wp), parameter :: answer = 4.25368768812224946110743394858422_wp
        real(wp), dimension(npts2), parameter :: points = [1.0_wp/7.0_wp, &
                                                           2.0_wp/3.0_wp, &
                                                           0.0_wp, &
                                                           0.0_wp]

        real(wp) :: abserr, result, work(lenw)
        integer :: ier, iwork(leniw), last, neval

        call dqagp(f, a, b, npts2, points, epsabs, epsrel, result, abserr, &
                   neval, ier, leniw, lenw, last, iwork, work)

        ! answer from maxima: quad_qags(abs(x-1/7)^(-0.25)*abs(x-2/3)^(-0.55), x, 0, 1);
        call check_result('dqagp', result, answer, neval)

    contains

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            f = 0.0_wp
            if (x /= 1.0_wp/7.0_wp .and. x /= 2.0_wp/3.0_wp) f = &
                abs(x - 1.0_wp/7.0_wp)**(-0.25_wp)* &
                abs(x - 2.0_wp/3.0_wp)**(-0.55_wp)
            return
        end function f

    end subroutine test_qagp

    subroutine test_qags
        implicit none

        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: answer = 2.0_wp

        real(wp) :: abserr, result, work(400)
        integer :: ier, iwork(100), last, neval

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

        real(wp), parameter :: a = -1.0_wp
        real(wp), parameter :: b = 1.0_wp
        real(wp), parameter :: c = 0.5_wp
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: answer = -628.461728506562332312831199677428_wp

        real(wp) :: abserr, result, work(400)
        integer :: ier, iwork(100), last, neval

        call dqawc(f, a, b, c, epsabs, epsrel, result, abserr, &
                   neval, ier, limit, lenw, last, iwork, work)

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

        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: omega = 8.0_wp
        integer, parameter :: integr = 2
        integer, parameter :: limlst = 100
        integer, parameter :: limit = 500
        integer, parameter :: leniw = limit*2 + limlst
        integer, parameter :: maxp1 = 100
        integer, parameter :: lenw = leniw*2 + maxp1*25
        real(wp), parameter :: answer = sqrt(29.0_wp*pi) - sqrt(21.0_wp*pi)

        real(wp) :: abserr, result, work(lenw)
        integer :: ier, iwork(leniw), last, lst, neval

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

        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        real(wp), parameter :: omega = 10.0_wp
        integer, parameter :: integr = 1
        integer, parameter :: limit = 100
        integer, parameter :: leniw = limit*2
        integer, parameter :: maxp1 = 21
        integer, parameter :: lenw = limit*4 + maxp1*25
        real(wp), parameter :: answer = -0.177639206511388980501003222731069_wp

        real(wp) :: abserr, result, work(925)
        integer :: ier, iwork(200), last, neval

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

        real(wp), parameter :: a = 0.0_wp
        real(wp), parameter :: b = 1.0_wp
        real(wp), parameter :: answer = 1.27072413983362022013785374440150_wp

        real(wp) :: abserr, result
        integer :: ier, neval

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
