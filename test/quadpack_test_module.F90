#ifndef MOD_INCLUDE
module quadpack_test_module
    use quadpack, wp => quadpack_RK
#endif

    implicit none

    public

    real(wp), parameter :: epsabs = 0.0_wp
    real(wp), parameter :: epsrel = 10**(log10(epsilon(1.0_wp))/2.0_wp+1)

    abstract interface
        real(wp) function func(x)
            import
            implicit none
            real(wp), intent(in) :: x
        end function func
    end interface

contains

    pure real(wp) function Catalan()
        !! Catalan's constant
        integer :: k
        real(wp) :: term
        Catalan = 0.0_wp
        k = 0
        do
            term = (-1)**k / (2.0_wp*k + 1.0_wp)**2
            if (Catalan == Catalan + term) exit
            Catalan = Catalan + term
            k = k + 1
        end do
    end function Catalan

    subroutine check_result(routine, value, answer, neval)
        implicit none

        character(len=*), intent(in) :: routine
        real(wp), intent(in) :: value
        real(wp), intent(in) :: answer
        integer, intent(in) :: neval

        write (*, '(1P,A25,1X,2(E13.6,1X),I6)') &
            routine, value, abs(value - answer), neval

        if (abs(value - answer) > epsrel) then
            write(*,*) '  value  = ', value
            write(*,*) '  answer = ', answer
            write(*,*) '  epsrel = ', epsrel
            write(*,*) 'TEST FAILED'
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

    subroutine test_C
        !! Tests from:
        !!
        !! * David H. Bailey, Karthik Jeyabalan, and Xiaoye S. Li,
        !!   "A Comparison of Three High-Precision Quadrature Schemes",
        !!   Experimental Mathematics 14:3
        implicit none

        integer, parameter :: key = 6
        integer, parameter :: limit = 100
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: pi = acos(-1.0_wp)

        procedure(func), pointer :: test_func
        real(wp) :: a, b, abserr, result, work(lenw), answer, G
        integer :: ier, iwork(limit), last, neval, i
        character(len=:), allocatable :: casename

        G = Catalan()

        do i = 1, 8
            select case (i)
            case(1)
                casename = 'C(1)'
                a = 0.0_wp
                b = 1.0_wp
                test_func => f
                answer = pi/4.0_wp - pi*sqrt(2.0_wp)/2.0_wp + &
                         3.0_wp*sqrt(2.0_wp)*atan(sqrt(2.0_wp))/2.0_wp
            case(2)
                casename = 'f1'
                a = 0.0_wp
                b = 1.0_wp
                test_func => f1
                answer = 0.25_wp
            case(3)
                casename = 'f2'
                a = 0.0_wp
                b = 1.0_wp
                test_func => f2
                answer = (pi - 2.0_wp + 2.0_wp * log(2.0_wp)) / 12.0_wp
            case(4)
                casename = 'f3'
                a = 0.0_wp
                b = pi/2.0_wp
                test_func => f3
                answer = (exp(pi/2.0_wp) - 1.0_wp) / 2.0_wp
            case(5)
                casename = 'f4'
                a = 0.0_wp
                b = 1.0_wp
                test_func => f4
                answer = 5.0_wp * pi**2 / 96.0_wp
            case(6)
                casename = 'i1'
                a = 0.0_wp
                b = 1.0_wp
                test_func => i1
                answer = pi**2*(2.0_wp - sqrt(2.0_wp)) / 32.0_wp
            case(7)
                casename = 'i2'
                a = 0.0_wp
                b = pi / 4.0_wp
                test_func => i2
                answer = -pi**2 / 16.0_wp + pi * log(2.0_wp) / 4.0_wp + G
            case(8)
                casename = 'i3'
                a = 0.0_wp
                b = pi
                test_func => i3
                answer = pi**2 / 4.0_wp
            end select
            call test_case(casename, test_func, a, b, answer)
        end do

    contains

        subroutine test_case(casename, f, a, b, answer)
            implicit none
            character(len=*),intent(in) :: casename
            procedure(func) :: f
            real(wp),intent(in) :: a
            real(wp),intent(in) :: b
            real(wp),intent(in) :: answer

            call dqag(f, a, b, epsabs, epsrel, key, result, abserr, neval, &
                      ier, limit, lenw, last, iwork, work)
            call check_result('dqag '//casename, result, answer, neval)

            call dqng(f, a, b, epsabs, epsrel, result, abserr, neval, ier)
            call check_result('dqng '//casename, result, answer, neval)

        end subroutine test_case

        real(wp) function f(x)
            implicit none
            real(wp), intent(in) :: x
            real(wp), parameter :: aa = 1.0_wp
            f = atan(sqrt(x**2+aa**2))/(sqrt(x**2+aa**2)*(x**2+1.0_wp))
        end function f

        real(wp) function f1(x)
            implicit none
            real(wp), intent(in) :: x
            f1 = x * log(1.0_wp + x)
        end function f1

        real(wp) function f2(x)
            implicit none
            real(wp), intent(in) :: x
            f2 = x**2 * atan(x)
        end function f2

        real(wp) function f3(x)
            implicit none
            real(wp), intent(in) :: x
            f3 = exp(x) * cos(x)
        end function f3

        real(wp) function f4(x)
            implicit none
            real(wp), intent(in) :: x
            f4 = atan(sqrt(2.0_wp + x**2)) / ((1.0_wp + x**2)*sqrt(2.0_wp + x**2))
        end function f4

        real(wp) function i1(x)
            implicit none
            real(wp), intent(in) :: x
            i1 = (x**2 * log(x)) / ((x**2-1.0_wp)*(x**4+1.0_wp))
        end function i1

        real(wp) function i2(x)
            implicit none
            real(wp), intent(in) :: x
            i2 = x**2 / sin(x)**2
        end function i2

        real(wp) function i3(x)
            implicit none
            real(wp), intent(in) :: x
            i3 =  (x * sin(x)) / (1 + cos(x)**2)
        end function i3

    end subroutine test_C

    subroutine test_G
        !! Catalan's constant integrals
        !! See: https://en.wikipedia.org/wiki/Catalan%27s_constant
        implicit none

        integer, parameter :: key = 6
        integer, parameter :: limit = 200
        integer, parameter :: lenw = limit*4
        real(wp), parameter :: pi = acos(-1.0_wp)

        procedure(func), pointer :: test_func
        real(wp) :: a, b, abserr, result, work(lenw), G
        integer :: ier, iwork(limit), last, neval, i
        character(len=:), allocatable :: casename

        G = Catalan()

        do i = 2, 2
            select case (i)
            case(1)
                casename = 'g1'
                a = 0.0_wp
                b = 1.0_wp
                test_func => g1
            case(2)
                casename = 'g2'
                a = 0.0_wp
                b = pi / 4.0_wp
                test_func => g2
            case(3)
                casename = 'g3'
                a = - pi / 2.0_wp
                b = pi / 2.0_wp
                test_func => g3
            end select
            call test_case(casename, test_func, a, b, G)
        end do

    contains

        subroutine test_case(casename, f, a, b, answer)
            implicit none
            character(len=*),intent(in) :: casename
            procedure(func) :: f
            real(wp),intent(in) :: a
            real(wp),intent(in) :: b
            real(wp),intent(in) :: answer

            call dqag(f, a, b, epsabs, epsrel, key, result, abserr, neval, &
                      ier, limit, lenw, last, iwork, work)
            call check_result('dqag '//casename, result, answer, neval)

            call dqng(f, a, b, epsabs, epsrel, result, abserr, neval, ier)
            call check_result('dqng '//casename, result, answer, neval)

        end subroutine test_case

        real(wp) function g1(x)
            real(wp), intent(in) :: x
            g1 = - log(x) / (1.0_wp + x**2)
        end function g1

        real(wp) function g2(x)
            real(wp), intent(in) :: x
            g2 = x / (sin(x)*cos(x))
        end function g2

        real(wp) function g3(x)
            real(wp), intent(in) :: x
            g3 = x / sin(x) / 4.0_wp
        end function g3

    end subroutine test_G

#ifndef MOD_INCLUDE
end module quadpack_test_module
#endif
