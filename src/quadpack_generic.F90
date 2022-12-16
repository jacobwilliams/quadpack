
!********************************************************************************
!>
!  Modernized QUADPACK: a Fortran subroutine package for the numerical
!  computation of definite one-dimensional integrals
!
!### References
!  * Original version on [Netlib](http://www.netlib.org/quadpack/)
!
!### Authors
!  * Piessens, Robert. Applied Mathematics and Programming Division, K. U. Leuven
!  * de Doncker, Elise. Applied Mathematics and Programming Division, K. U. Leuven
!  * Kahaner, D. K., (NBS)
!  * Jacob Williams, Dec 2021. Modernized the Fortran 77 code from Netlib.

#ifndef MOD_INCLUDE
module quadpack_generic
    use iso_fortran_env, only: wp => real64 ! double precision by default
#endif

    implicit none

    private

    integer, parameter, public :: quadpack_RK = wp !! the real kind used in this module

    real(wp), dimension(5), parameter, private :: d1mach = [tiny(1.0_wp), &
                                                            huge(1.0_wp), &
                                                            real(radix(1.0_wp),kind(1.0_wp))**(-digits(1.0_wp)), &
                                                            epsilon(1.0_wp), &
                                                            log10(real(radix(1.0_wp), kind(1.0_wp)))] !! machine constants
    integer,parameter :: i1mach10 = radix(1.0_wp)
    integer,parameter :: i1mach14 = digits(1.0_wp)

    real(wp), parameter, private :: uflow = d1mach(1) !! the smallest positive magnitude.
    real(wp), parameter, private :: oflow = d1mach(2) !! the largest positive magnitude.
    real(wp), parameter, private :: epmach = d1mach(4) !! the largest relative spacing.
    real(wp), parameter, private :: pi = acos(-1.0_wp) !! pi

    integer, parameter, private :: limexp = 50 !! `limexp` is the maximum number of elements the epsilon
                                               !! table can contain. if this number is reached, the upper
                                               !! diagonal of the epsilon table is deleted.
                                               !! originally defined in [[dqelg]]. Was moved to be a module
                                               !! variable since various dimensions in other routines
                                               !! depend on the value

    abstract interface

        real(wp) function func(x)
        !! interface for user-supplied function.
            import :: wp
            implicit none
            real(wp), intent(in) :: x
        end function func

        real(wp) function weight_func(x, a, b, c, d, i)
         !! weight function interface for [[dqk15w]]
            import :: wp
            implicit none
            real(wp), intent(in) :: x
            real(wp), intent(in) :: a
            real(wp), intent(in) :: b
            real(wp), intent(in) :: c
            real(wp), intent(in) :: d
            integer, intent(in) :: i
        end function weight_func

    end interface

    ! by default, the double precision names are exported  (dqag, etc.)
    public :: dqag, dqage, dqagi, dqagie, dqagp, dqagpe, dqags, &
              dqagse, dqawc, dqawce, dqawf, dqawfe, dqawo, dqawoe, dqaws, &
              dqawse, dqc25c, dqc25f, dqc25s, dqcheb, dqk15, dqk15i, &
              dqk15w, dqk21, dqk31, dqk41, dqk51, dqk61, dqmomo, dqng
    public :: dquad
    public :: davint
    public :: dqnc79
    public :: dgauss8
    public :: dsimpson, dlobatto

    contains
!********************************************************************************

!********************************************************************************
!>
!  1D globally adaptive integrator using Gauss-Kronrod quadrature, oscillating integrand
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqag(f, a, b, Epsabs, Epsrel, Key, Result, Abserr, Neval, Ier, &
                    Limit, Lenw, Last, Iwork, Work)

        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if epsabs<=0
                                       !! and epsrel<max(50*rel.mach.acc.,0.5e-28),
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`
                                    !! `lenw` must be at least `limit*4`.
                                    !! if `lenw<limit*4`, the routine will end with
                                    !! ier = 6.
        integer, intent(in) :: Limit !! dimensioning parameter for `iwork`
                                     !! limit determines the maximum number of subintervals
                                     !! in the partition of the given integration interval
                                     !! (a,b), limit>=1.
                                     !! if limit<1, the routine will end with ier = 6.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`
                               !! on return
                               !! work(1), ..., work(last) contain the left end
                               !! points of the subintervals in the partition of
                               !!  (a,b),
                               !! `work(limit+1), ..., work(limit+last)` contain the
                               !!  right end points,
                               !! `work(limit*2+1), ..., work(limit*2+last)` contain
                               !!  the integral approximations over the subintervals,
                               !! work(limit*3+1), ..., work(limit*3+last) contain
                               !!  the error estimates.
        integer :: Iwork(Limit) !! vector of dimension at least `limit`, the first `k`
                                !! elements of which contain pointers to the error
                                !! estimates over the subintervals, such that
                                !! work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
                                !! form a decreasing sequence with k = last if
                                !! last<=(limit/2+2), and k = limit+1-last otherwise
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!         routine. it is assumed that the requested
                                    !!         accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!         the estimates for result and error are
                                    !!         less reliable. it is assumed that the
                                    !!         requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!         has been achieved. one can allow more
                                    !!         subdivisions by increasing the value of
                                    !!         limit (and taking the according dimension
                                    !!         adjustments into account). however, if
                                    !!         this yield no improvement it is advised
                                    !!         to analyze the integrand in order to
                                    !!         determine the integration difficulties.
                                    !!         if the position of a local difficulty can
                                    !!         be determined (i.e.singularity,
                                    !!         discontinuity within the interval) one
                                    !!         will probably gain from splitting up the
                                    !!         interval at this point and calling the
                                    !!         integrator on the subranges. if possible,
                                    !!         an appropriate special-purpose integrator
                                    !!         should be used which is designed for
                                    !!         handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!         detected, which prevents the requested
                                    !!         tolerance from being achieved.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!         at some points of the integration
                                    !!         interval.
                                    !! * ier = 6 the input is invalid, because
                                    !!         `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28_wp))`
                                    !!         or `limit<1` or `lenw<limit*4`.
                                    !!         `result`, `abserr`, `neval`, `last` are set
                                    !!         to zero.
                                    !!         except when lenw is invalid, `iwork(1)`,
                                    !!         `work(limit*2+1)` and `work(limit*3+1)` are
                                    !!         set to zero, `work(1)` is set to a and
                                    !!         `work(limit+1)` to `b`.
        integer, intent(in) :: Key !! key for choice of local integration rule.
                                   !! a gauss-kronrod pair is used with:
                                   !!
                                   !!  *  7 - 15 points if key<2,
                                   !!  * 10 - 21 points if key = 2,
                                   !!  * 15 - 31 points if key = 3,
                                   !!  * 20 - 41 points if key = 4,
                                   !!  * 25 - 51 points if key = 5,
                                   !!  * 30 - 61 points if key>5.
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, which
                                     !! determines the number of significant elements
                                     !! actually in the work arrays.
        integer, intent(out) :: Neval !! number of integrand evaluations

        integer :: lvl, l1, l2, l3

        ! check validity of lenw.
        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Limit >= 1 .and. Lenw >= Limit*4) then

            ! prepare call for dqage.

            l1 = Limit + 1
            l2 = Limit + l1
            l3 = Limit + l2

            call dqage(f, a, b, Epsabs, Epsrel, Key, Limit, Result, Abserr, Neval, &
                       Ier, Work(1), Work(l1), Work(l2), Work(l3), Iwork, Last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqag ', Ier, lvl)

    end subroutine dqag
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqag]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-reslt)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqage(f, a, b, Epsabs, Epsrel, Key, Limit, Result, Abserr, &
                     Neval, Ier, Alist, Blist, Rlist, Elist, Iord, Last)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and epsrel<max(50*rel.mach.acc.,0.5e-28),
                                       !! the routine will end with ier = 6.
        integer, intent(in) :: Key !! key for choice of local integration rule
                                   !!  a gauss-kronrod pair is used with
                                   !!
                                   !!  * 7 - 15 points if key<2,
                                   !!  * 10 - 21 points if key = 2,
                                   !!  * 15 - 31 points if key = 3,
                                   !!  * 20 - 41 points if key = 4,
                                   !!  * 25 - 51 points if key = 5,
                                   !!  * 30 - 61 points if key>5.
        integer, intent(in) :: Limit !! gives an upperbound on the number of subintervals
                                     !! in the partition of `(a,b)`, `limit>=1`.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !!  * ier = 0 normal and reliable termination of the
                                    !!    routine. it is assumed that the requested
                                    !!    accuracy has been achieved.
                                    !!  * ier>0 abnormal termination of the routine
                                    !!    the estimates for result and error are
                                    !!    less reliable. it is assumed that the
                                    !!    requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !!  * ier = 1 maximum number of subdivisions allowed
                                    !!          has been achieved. one can allow more
                                    !!          subdivisions by increasing the value
                                    !!          of limit.
                                    !!          however, if this yields no improvement it
                                    !!          is rather advised to analyze the integrand
                                    !!          in order to determine the integration
                                    !!          difficulties. if the position of a local
                                    !!          difficulty can be determined(e.g.
                                    !!          singularity, discontinuity within the
                                    !!          interval) one will probably gain from
                                    !!          splitting up the interval at this point
                                    !!          and calling the integrator on the
                                    !!          subranges. if possible, an appropriate
                                    !!          special-purpose integrator should be used
                                    !!          which is designed for handling the type of
                                    !!          difficulty involved.
                                    !!  * ier = 2 the occurrence of roundoff error is
                                    !!          detected, which prevents the requested
                                    !!          tolerance from being achieved.
                                    !!  * ier = 3 extremely bad integrand behaviour occurs
                                    !!          at some points of the integration
                                    !!          interval.
                                    !!  * ier = 6 the input is invalid, because
                                    !!          (epsabs<=0 and
                                    !!           epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
                                    !!          result, abserr, neval, last, rlist(1) ,
                                    !!          `elist(1)` and `iord(1)` are set to zero.
                                    !!          alist(1) and blist(1) are set to a and b
                                    !!          respectively.
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the moduli of the
                                              !! absolute error estimates on the subintervals
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the
                                              !! integral approximations on the subintervals
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! elements of which are pointers to the
                                            !! error estimates over the subintervals,
                                            !! such that `elist(iord(1))`, ...,
                                            !! `elist(iord(k))` form a decreasing sequence,
                                            !! with `k = last` if `last<=(limit/2+2)`, and
                                            !! `k = limit+1-last` otherwise
        integer, intent(out) :: Last !! number of subintervals actually produced in the
                                     !! subdivision process

        real(wp) :: area1, a1, b1, defab1, error1 !! variable for the left subinterval
        real(wp) :: area2, a2, b2, defab2, error2 !! variable for the right subinterval
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errmax !! `elist(maxerr)`
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        integer :: maxerr !! pointer to the interval with largest error estimate
        real(wp) :: resabs, defabs
        integer :: iroff1, iroff2, k, keyf, nrmax

        ! test on validity of parameters

        Ier = 0
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        Alist(1) = a
        Blist(1) = b
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        if (Epsabs <= 0.0_wp .and. Epsrel < max(50.0_wp*epmach, 0.5e-28_wp)) Ier = 6
        if (Ier /= 6) then

            ! first approximation to the integral

            keyf = Key
            if (Key <= 0) keyf = 1
            if (Key >= 7) keyf = 6
            Neval = 0
            select case (keyf)
            case (1); call dqk15(f, a, b, Result, Abserr, defabs, resabs)
            case (2); call dqk21(f, a, b, Result, Abserr, defabs, resabs)
            case (3); call dqk31(f, a, b, Result, Abserr, defabs, resabs)
            case (4); call dqk41(f, a, b, Result, Abserr, defabs, resabs)
            case (5); call dqk51(f, a, b, Result, Abserr, defabs, resabs)
            case (6); call dqk61(f, a, b, Result, Abserr, defabs, resabs)
            end select
            Last = 1
            Rlist(1) = Result
            Elist(1) = Abserr
            Iord(1) = 1

            ! test on accuracy.

            errbnd = max(Epsabs, Epsrel*abs(Result))
            if (Abserr <= 50.0_wp*epmach*defabs .and. Abserr > errbnd) Ier = 2
            if (Limit == 1) Ier = 1

            if (.not. (Ier /= 0 .or. (Abserr <= errbnd .and. Abserr /= resabs) &
                       .or. Abserr == 0.0_wp)) then

                ! initialization
                errmax = Abserr
                maxerr = 1
                area = Result
                errsum = Abserr
                nrmax = 1
                iroff1 = 0
                iroff2 = 0

                ! main do-loop

                do Last = 2, Limit

                    ! bisect the subinterval with the largest error estimate.

                    a1 = Alist(maxerr)
                    b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                    a2 = b1
                    b2 = Blist(maxerr)
                    select case (keyf)
                    case (1)
                        call dqk15(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk15(f, a2, b2, area2, error2, resabs, defab2)
                    case (2)
                        call dqk21(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk21(f, a2, b2, area2, error2, resabs, defab2)
                    case (3)
                        call dqk31(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk31(f, a2, b2, area2, error2, resabs, defab2)
                    case (4)
                        call dqk41(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk41(f, a2, b2, area2, error2, resabs, defab2)
                    case (5)
                        call dqk51(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk51(f, a2, b2, area2, error2, resabs, defab2)
                    case (6)
                        call dqk61(f, a1, b1, area1, error1, resabs, defab1)
                        call dqk61(f, a2, b2, area2, error2, resabs, defab2)
                    end select

                    ! improve previous approximations to integral
                    ! and error and test for accuracy.

                    Neval = Neval + 1
                    area12 = area1 + area2
                    erro12 = error1 + error2
                    errsum = errsum + erro12 - errmax
                    area = area + area12 - Rlist(maxerr)
                    if (defab1 /= error1 .and. defab2 /= error2) then
                        if (abs(Rlist(maxerr) - area12) <= 0.1e-4_wp*abs(area12) &
                            .and. erro12 >= 0.99_wp*errmax) iroff1 = iroff1 + 1
                        if (Last > 10 .and. erro12 > errmax) iroff2 = iroff2 + 1
                    end if
                    Rlist(maxerr) = area1
                    Rlist(Last) = area2
                    errbnd = max(Epsabs, Epsrel*abs(area))
                    if (errsum > errbnd) then

                        ! test for roundoff error and eventually set error flag.

                        if (iroff1 >= 6 .or. iroff2 >= 20) Ier = 2

                        ! set error flag in the case that the number of subintervals
                        ! equals limit.

                        if (Last == Limit) Ier = 1

                        ! set error flag in the case of bad integrand behaviour
                        ! at a point of the integration range.

                        if (max(abs(a1), abs(b2)) &
                            <= (1.0_wp + 100.0_wp*epmach) &
                            *(abs(a2) + 1000.0_wp*uflow)) Ier = 3
                    end if

                    ! append the newly-created intervals to the list.

                    if (error2 > error1) then
                        Alist(maxerr) = a2
                        Alist(Last) = a1
                        Blist(Last) = b1
                        Rlist(maxerr) = area2
                        Rlist(Last) = area1
                        Elist(maxerr) = error2
                        Elist(Last) = error1
                    else
                        Alist(Last) = a2
                        Blist(maxerr) = b1
                        Blist(Last) = b2
                        Elist(maxerr) = error1
                        Elist(Last) = error2
                    end if

                    ! call subroutine dqpsrt to maintain the descending ordering
                    ! in the list of error estimates and select the subinterval
                    ! with the largest error estimate (to be bisected next).

                    call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                    if (Ier /= 0 .or. errsum <= errbnd) exit  ! jump out of do-loop
                end do

                ! compute final result.

                Result = 0.0_wp
                do k = 1, Last
                    Result = Result + Rlist(k)
                end do
                Abserr = errsum
            end if
            if (keyf /= 1) Neval = (10*keyf + 1)*(2*Neval + 1)
            if (keyf == 1) Neval = 30*Neval + 15
        end if

    end subroutine dqage
!********************************************************************************

!********************************************************************************
!>
!  1D globally adaptive integrator, infinite intervals
!
!  the routine calculates an approximation result to a given
!  integral with one of the following forms:
!
!  * i = integral of `f` over `(bound, +infinity)`
!  * i = integral of `f` over `(-infinity, bound)`
!  * i = integral of `f` over `(-infinity, +infinity)`
!
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqagi(f, Bound, Inf, Epsabs, Epsrel, Result, Abserr, Neval, &
                     Ier, Limit, Lenw, Last, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(in) :: Bound !! finite bound of integration range
                                      !! (has no meaning if interval is doubly-infinite)
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if  `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`
                                    !! `lenw` must be at least `limit*4`.
                                    !! if `lenw<limit*4`, the routine will end
                                    !! with ier = 6.
        integer, intent(in) :: Limit !! dimensioning parameter for `iwork`
                                     !! limit determines the maximum number of subintervals
                                     !! in the partition of the given integration interval
                                     !! (a,b), `limit>=1`.
                                     !! if `limit<1`, the routine will end with ier = 6.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`
                               !! on return:
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !! * `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end points,
                               !! * `work(limit*2+1), ...,work(limit*2+last)` contain the
                               !!   integral approximations over the subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3)`
                               !!   contain the error estimates.
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!           routine. it is assumed that the requested
                                    !!           accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine. the
                                    !!           estimates for result and error are less
                                    !!           reliable. it is assumed that the requested
                                    !!           accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!         has been achieved. one can allow more
                                    !!         subdivisions by increasing the value of
                                    !!         limit (and taking the according dimension
                                    !!         adjustments into account). however, if
                                    !!         this yields no improvement it is advised
                                    !!         to analyze the integrand in order to
                                    !!         determine the integration difficulties. if
                                    !!         the position of a local difficulty can be
                                    !!         determined (e.g. singularity,
                                    !!         discontinuity within the interval) one
                                    !!         will probably gain from splitting up the
                                    !!         interval at this point and calling the
                                    !!         integrator on the subranges. if possible,
                                    !!         an appropriate special-purpose integrator
                                    !!         should be used, which is designed for
                                    !!         handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!         detected, which prevents the requested
                                    !!         tolerance from being achieved.
                                    !!         the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!         at some points of the integration
                                    !!         interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!         roundoff error is detected in the
                                    !!         extrapolation table.
                                    !!         it is assumed that the requested tolerance
                                    !!         cannot be achieved, and that the returned
                                    !!         result is the best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!         slowly convergent. it must be noted that
                                    !!         divergence can occur with any other value
                                    !!         of ier.
                                    !! * ier = 6 the input is invalid, because
                                    !!         `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!         or `limit<1` or `leniw<limit*4`.
                                    !!         `result`, `abserr`, `neval`, `last` are set to
                                    !!         zero. except when `limit` or `leniw` is
                                    !!         invalid, `iwork(1)`, `work(limit*2+1)` and
                                    !!         `work(limit*3+1)` are set to zero, `work(1)`
                                    !!         is set to `a` and `work(limit+1)` to `b`.
        integer, intent(in) :: Inf !! indicating the kind of integration range involved:
                                   !!
                                   !! * inf = 1 corresponds to `(bound,+infinity)`
                                   !! * inf = -1 corresponds to `(-infinity,bound)`
                                   !! * inf = 2 corresponds to `(-infinity,+infinity)`
        integer :: Iwork(Limit) !! vector of dimension at least `limit`, the first
                                 !! `k` elements of which contain pointers
                                 !! to the error estimates over the subintervals,
                                 !! such that `work(limit*3+iwork(1)),...,work(limit*3+iwork(k))`
                                 !! form a decreasing sequence, with `k = last`
                                 !! if `last<=(limit/2+2)`, and `k = limit+1-last` otherwise
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, which
                                     !! determines the number of significant elements
                                     !! actually in the work arrays.
        integer, intent(out) :: Neval !! number of integrand evaluations

        integer :: lvl, l1, l2, l3

        ! check validity of limit and lenw.
        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Limit >= 1 .and. Lenw >= Limit*4) then

            ! prepare call for dqagie.
            l1 = Limit + 1
            l2 = Limit + l1
            l3 = Limit + l2

            call dqagie(f, Bound, Inf, Epsabs, Epsrel, Limit, Result, Abserr, &
                        Neval, Ier, Work(1), Work(l1), Work(l2), Work(l3), Iwork, &
                        Last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqagi', Ier, lvl)

    end subroutine dqagi
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqagi]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  integral with one of the following forms:
!
!  * i = integral of `f` over `(bound, +infinity)`
!  * i = integral of `f` over `(-infinity, bound)`
!  * i = integral of `f` over `(-infinity, +infinity)`
!
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqagie(f, Bound, Inf, Epsabs, Epsrel, Limit, Result, Abserr, &
                      Neval, Ier, Alist, Blist, Rlist, Elist, Iord, Last)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        integer, intent(in) :: Limit !! gives an upper bound on the number of subintervals
                                     !! in the partition of `(a,b)`, `limit>=1`
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left
                                              !! end points of the subintervals in the partition
                                              !! of the transformed integration range (0,1).
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right
                                              !! end points of the subintervals in the partition
                                              !! of the transformed integration range (0,1).
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`,  the first
                                              !! `last` elements of which are the moduli of the
                                              !! absolute error estimates on the subintervals
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the integral
                                              !! approximations on the subintervals
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with `ier = 6`.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(in) :: Bound !! finite bound of integration range
                                      !! (has no meaning if interval is doubly-infinite)
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine. the
                                    !!   estimates for result and error are less
                                    !!   reliable. it is assumed that the requested
                                    !!   accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit (and taking the according dimension
                                    !!   adjustments into account). however,if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties.
                                    !!   if the position of a local difficulty can
                                    !!   be determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling the
                                    !!   integrator on the subranges. if possible,
                                    !!   an appropriate special-purpose integrator
                                    !!   should be used, which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table.
                                    !!   it is assumed that the requested tolerance
                                    !!   cannot be achieved, and that the returned
                                    !!   result is the best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                    !!   `result`, `abserr`, `neval`, `last`, `rlist(1)`,
                                    !!   `elist(1)` and `iord(1)` are set to zero.
                                    !!   `alist(1)` and `blist(1)` are set to 0
                                    !!   and 1 respectively.
        integer, intent(in) :: Inf !! indicating the kind of integration range involved
                                   !! * inf = 1  corresponds to `(bound,+infinity)`
                                   !! * inf = -1 corresponds to `(-infinity,bound)`
                                   !! * inf = 2  corresponds to `(-infinity,+infinity)`
        integer, intent(out) :: Iord(Limit) !! vector of dimension `limit`, the first `k`
                                            !! elements of which are pointers to the
                                            !! error estimates over the subintervals,
                                            !! such that `elist(iord(1)), ..., elist(iord(k))`
                                            !! form a decreasing sequence, with `k = last`
                                            !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                            !! otherwise
        integer, intent(out) :: Last !! number of subintervals actually produced
                                     !! in the subdivision process
        integer, intent(out) :: Neval !! number of integrand evaluations

        real(wp) :: area1, a1, b1, defab1, error1 !! variable for the left subinterval
        real(wp) :: area2, a2, b2, defab2, error2 !! variable for the right subinterval
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: errmax !! `elist(maxerr)`
        real(wp) :: erlast !! error on the interval currently subdivided
                           !! (before that subdivision has taken place)
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        real(wp) :: small !! length of the smallest interval considered up
                          !! to now, multiplied by 1.5
        real(wp) :: erlarg !! sum of the errors over the intervals larger
                           !! than the smallest interval considered up to now
        integer :: maxerr !! pointer to the interval with largest error estimate
        integer :: nres !! number of calls to the extrapolation routine
        integer :: numrl2 !! number of elements currently in rlist2. if an
                          !! appropriate approximation to the compounded
                          !! integral has been obtained, it is put in
                          !! rlist2(numrl2) after numrl2 has been increased
                          !! by one.
        logical :: extrap !! logical variable denoting that the routine
                          !! is attempting to perform extrapolation. i.e.
                          !! before subdividing the smallest interval we
                          !! try to decrease the value of erlarg.
        logical :: noext !! logical variable denoting that extrapolation
                         !! is no longer allowed (true-value)
        real(wp) :: rlist2(limexp + 2) !! array of dimension at least (`limexp+2`),
                                       !! containing the part of the epsilon table
                                       !! which is still needed for further computations.
        real(wp) :: abseps, boun, correc, defabs, dres, &
                    ertest, resabs, reseps, res3la(3)
        integer :: id, ierro, iroff1, iroff2, iroff3, &
                   jupbnd, k, ksgn, ktmin, nrmax

        ! test on validity of parameters

        Ier = 0
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        Alist(1) = 0.0_wp
        Blist(1) = 1.0_wp
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        if (Epsabs <= 0.0_wp .and. Epsrel < max(50.0_wp*epmach, 0.5e-28_wp)) &
            Ier = 6
        if (Ier == 6) return

        main : block

            ! first approximation to the integral

            ! determine the interval to be mapped onto (0,1).
            ! if inf = 2 the integral is computed as i = i1+i2, where
            ! i1 = integral of f over (-infinity,0),
            ! i2 = integral of f over (0,+infinity).

            boun = Bound
            if (Inf == 2) boun = 0.0_wp
            call dqk15i(f, boun, Inf, 0.0_wp, 1.0_wp, Result, Abserr, defabs, &
                        resabs)

            ! test on accuracy

            Last = 1
            Rlist(1) = Result
            Elist(1) = Abserr
            Iord(1) = 1
            dres = abs(Result)
            errbnd = max(Epsabs, Epsrel*dres)
            if (Abserr <= 100.0_wp*epmach*defabs .and. Abserr > errbnd) Ier = 2
            if (Limit == 1) Ier = 1
            if (Ier /= 0 .or. (Abserr <= errbnd .and. Abserr /= resabs) .or. &
                Abserr == 0.0_wp) exit main

            ! initialization

            rlist2(1) = Result
            errmax = Abserr
            maxerr = 1
            area = Result
            errsum = Abserr
            Abserr = oflow
            nrmax = 1
            nres = 0
            ktmin = 0
            numrl2 = 2
            extrap = .false.
            noext = .false.
            ierro = 0
            iroff1 = 0
            iroff2 = 0
            iroff3 = 0
            ksgn = -1
            if (dres >= (1.0_wp - 50.0_wp*epmach)*defabs) ksgn = 1

            ! main do-loop

            loop: do Last = 2, Limit

                ! bisect the subinterval with nrmax-th largest error estimate.

                a1 = Alist(maxerr)
                b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                a2 = b1
                b2 = Blist(maxerr)
                erlast = errmax
                call dqk15i(f, boun, Inf, a1, b1, area1, error1, resabs, defab1)
                call dqk15i(f, boun, Inf, a2, b2, area2, error2, resabs, defab2)

                ! improve previous approximations to integral
                ! and error and test for accuracy.

                area12 = area1 + area2
                erro12 = error1 + error2
                errsum = errsum + erro12 - errmax
                area = area + area12 - Rlist(maxerr)
                if (defab1 /= error1 .and. defab2 /= error2) then
                    if (abs(Rlist(maxerr) - area12) <= 0.1e-4_wp*abs(area12) .and. &
                        erro12 >= 0.99_wp*errmax) then
                        if (extrap) iroff2 = iroff2 + 1
                        if (.not. extrap) iroff1 = iroff1 + 1
                    end if
                    if (Last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
                end if
                Rlist(maxerr) = area1
                Rlist(Last) = area2
                errbnd = max(Epsabs, Epsrel*abs(area))

                ! test for roundoff error and eventually set error flag.

                if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) Ier = 2
                if (iroff2 >= 5) ierro = 3

                ! set error flag in the case that the number of
                ! subintervals equals limit.

                if (Last == Limit) Ier = 1

                ! set error flag in the case of bad integrand behaviour
                ! at some points of the integration range.

                if (max(abs(a1), abs(b2)) <= (1.0_wp + 100.0_wp*epmach) &
                    *(abs(a2) + 1000.0_wp*uflow)) Ier = 4

                ! append the newly-created intervals to the list.

                if (error2 > error1) then
                    Alist(maxerr) = a2
                    Alist(Last) = a1
                    Blist(Last) = b1
                    Rlist(maxerr) = area2
                    Rlist(Last) = area1
                    Elist(maxerr) = error2
                    Elist(Last) = error1
                else
                    Alist(Last) = a2
                    Blist(maxerr) = b1
                    Blist(Last) = b2
                    Elist(maxerr) = error1
                    Elist(Last) = error2
                end if

                ! call subroutine dqpsrt to maintain the descending ordering
                ! in the list of error estimates and select the subinterval
                ! with nrmax-th largest error estimate (to be bisected next).

                call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                if (errsum <= errbnd) then
                    ! compute global integral sum.
                    Result = sum(Rlist(1:Last))
                    Abserr = errsum
                    exit main
                end if
                if (Ier /= 0) exit
                if (Last == 2) then
                    small = 0.375_wp
                    erlarg = errsum
                    ertest = errbnd
                    rlist2(2) = area
                elseif (.not. (noext)) then
                    erlarg = erlarg - erlast
                    if (abs(b1 - a1) > small) erlarg = erlarg + erro12
                    if (.not. (extrap)) then
                        ! test whether the interval to be bisected next is the
                        ! smallest interval.
                        if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle loop
                        extrap = .true.
                        nrmax = 2
                    end if
                    if (ierro /= 3 .and. erlarg > ertest) then

                        ! the smallest interval has the largest error.
                        ! before bisecting decrease the sum of the errors over the
                        ! larger intervals (erlarg) and perform extrapolation.

                        id = nrmax
                        jupbnd = Last
                        if (Last > (2 + Limit/2)) jupbnd = Limit + 3 - Last
                        do k = id, jupbnd
                            maxerr = Iord(nrmax)
                            errmax = Elist(maxerr)
                            if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle loop
                            nrmax = nrmax + 1
                        end do
                    end if

                    ! perform extrapolation.

                    numrl2 = numrl2 + 1
                    rlist2(numrl2) = area
                    call dqelg(numrl2, rlist2, reseps, abseps, res3la, nres)
                    ktmin = ktmin + 1
                    if (ktmin > 5 .and. Abserr < 0.1e-02_wp*errsum) Ier = 5
                    if (abseps < Abserr) then
                        ktmin = 0
                        Abserr = abseps
                        Result = reseps
                        correc = erlarg
                        ertest = max(Epsabs, Epsrel*abs(reseps))
                        if (Abserr <= ertest) exit
                    end if

                    ! prepare bisection of the smallest interval.

                    if (numrl2 == 1) noext = .true.
                    if (Ier == 5) exit
                    maxerr = Iord(1)
                    errmax = Elist(maxerr)
                    nrmax = 1
                    extrap = .false.
                    small = small*0.5_wp
                    erlarg = errsum
                end if

            end do loop

            ! set final result and error estimate.

            if (Abserr /= oflow) then
                if ((Ier + ierro) /= 0) then
                    if (ierro == 3) Abserr = Abserr + correc
                    if (Ier == 0) Ier = 3
                    if (Result == 0.0_wp .or. area == 0.0_wp) then
                        if (Abserr > errsum) then
                            ! compute global integral sum.
                            Result = sum(Rlist(1:Last))
                            Abserr = errsum
                            exit main
                        end if
                        if (area == 0.0_wp) exit main
                    elseif (Abserr/abs(Result) > errsum/abs(area)) then
                        ! compute global integral sum.
                        Result = sum(Rlist(1:Last))
                        Abserr = errsum
                        exit main
                    end if
                end if

                ! test on divergence

                if (ksgn /= (-1) .or. max(abs(Result), abs(area)) > defabs*0.01_wp) then
                    if (0.01_wp > (Result/area) .or. &
                        (Result/area) > 100.0_wp .or. &
                        errsum > abs(area)) Ier = 6
                end if

            else
                ! compute global integral sum.
                Result = sum(Rlist(1:Last))
                Abserr = errsum
            end if

        end block main

        Neval = 30*Last - 15
        if (Inf == 2) Neval = 2*Neval
        if (Ier > 2) Ier = Ier - 1

    end subroutine dqagie
!********************************************************************************

!********************************************************************************
!>
!  1D globally adaptive integrator, singularities or discontinuities
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  break points of the integration interval, where local
!  difficulties of the integrand may occur (e.g.
!  singularities, discontinuities), are provided by the user.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqagp(f, a, b, Npts2, Points, Epsabs, Epsrel, Result, Abserr, &
                     Neval, Ier, Leniw, Lenw, Last, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        integer, intent(in) :: Npts2 !! number equal to two more than the number of
                                     !! user-supplied break points within the integration
                                     !! range, `npts>=2`.
                                     !! if `npts2<2`, the routine will end with ier = 6.
        real(wp), intent(in) :: Points(Npts2) !! vector of dimension npts2, the first `(npts2-2)`
                                              !! elements of which are the user provided break
                                              !! points. if these points do not constitute an
                                              !! ascending sequence there will be an automatic
                                              !! sorting.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine.
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties. if
                                    !!   the position of a local difficulty can be
                                    !!   determined (i.e. singularity,
                                    !!   discontinuity within the interval), it
                                    !!   should be supplied to the routine as an
                                    !!   element of the vector points. if necessary
                                    !!   an appropriate special-purpose integrator
                                    !!   must be used, which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table.
                                    !!   it is presumed that the requested
                                    !!   tolerance cannot be achieved, and that
                                    !!   the returned result is the best which
                                    !!   can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier>0.
                                    !! * ier = 6 the input is invalid because
                                    !!   `npts2<2` or
                                    !!   break points are specified outside
                                    !!   the integration range or
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!   `result`, `abserr`, `neval`, `last` are set to
                                    !!   zero. except when `leniw` or `lenw` or `npts2` is
                                    !!   invalid, `iwork(1)`, `iwork(limit+1)`,
                                    !!   `work(limit*2+1) and work(limit*3+1)`
                                    !!   are set to zero.
                                    !!   `work(1)` is set to a and `work(limit+1)`
                                    !!   to `b` (where `limit = (leniw-npts2)/2`).
        integer, intent(in) :: Leniw !! dimensioning parameter for `iwork`.
                                     !! `leniw` determines `limit = (leniw-npts2)/2`,
                                     !! which is the maximum number of subintervals in the
                                     !! partition of the given integration interval `(a,b)`,
                                     !! `leniw>=(3*npts2-2)`.
                                     !! if `leniw<(3*npts2-2)`, the routine will end with
                                     !! ier = 6.
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`.
                                    !! `lenw` must be at least `leniw*2-npts2`.
                                    !! if `lenw<leniw*2-npts2`, the routine will end
                                    !! with ier = 6.
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, which
                                     !! determines the number of significant elements
                                     !! actually in the work arrays.
        integer :: Iwork(Leniw) !! vector of dimension at least `leniw`. on return,
                                !! the first `k` elements of which contain
                                !! pointers to the error estimates over the
                                !! subintervals, such that
                                !! `work(limit*3+iwork(1)),...,work(limit*3+iwork(k))`
                                !! form a decreasing
                                !! sequence, with `k = last` if `last<=(limit/2+2)`, and
                                !! `k = limit+1-last` otherwise
                                !! `iwork(limit+1), ...,iwork(limit+last)` contain the
                                !! subdivision levels of the subintervals, i.e.
                                !! if `(aa,bb)` is a subinterval of `(p1,p2)`
                                !! where `p1` as well as `p2` is a user-provided
                                !! break point or integration limit, then `(aa,bb)` has
                                !! level `l` if `abs(bb-aa) = abs(p2-p1)*2**(-l)`,
                                !! `iwork(limit*2+1), ..., iwork(limit*2+npts2)` have
                                !! no significance for the user,
                                !! note that `limit = (leniw-npts2)/2`.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`.
                               !! on return:
                               !!
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !! * `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end points,
                               !! * `work(limit*2+1), ..., work(limit*2+last)` contain
                               !!   the integral approximations over the subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3+last)`
                               !!   contain the corresponding error estimates,
                               !! * `work(limit*4+1), ..., work(limit*4+npts2)`
                               !!   contain the integration limits and the
                               !!   break points sorted in an ascending sequence.
                               !!
                               !! note that `limit = (leniw-npts2)/2`.

        integer :: limit, lvl, l1, l2, l3, l4

        ! check validity of limit and lenw.
        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Leniw >= (3*Npts2 - 2) .and. Lenw >= (Leniw*2 - Npts2) .and. Npts2 >= 2) then

            ! prepare call for dqagpe.
            limit = (Leniw - Npts2)/2
            l1 = limit + 1
            l2 = limit + l1
            l3 = limit + l2
            l4 = limit + l3

            call dqagpe(f, a, b, Npts2, Points, Epsabs, Epsrel, limit, Result, &
                        Abserr, Neval, Ier, Work(1), Work(l1), Work(l2), Work(l3), &
                        Work(l4), Iwork(1), Iwork(l1), Iwork(l2), Last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqagp', Ier, lvl)

    end subroutine dqagp
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqagp]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`, hopefully
!  satisfying following claim for accuracy `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!  break points of the integration interval, where local difficulties
!  of the integrand may occur (e.g. singularities, discontinuities),provided by user.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqagpe(f, a, b, Npts2, Points, Epsabs, Epsrel, Limit, Result, &
                      Abserr, Neval, Ier, Alist, Blist, Rlist, Elist, Pts, &
                      Iord, Level, Ndin, Last)
        implicit none

        procedure(func) :: f
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left end points
                                              !! of the subintervals in the partition of the given
                                              !! integration range (a,b)
        real(wp), intent(out)  :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                               !! `last` elements of which are the right end points
                                               !! of the subintervals in the partition of the given
                                               !! integration range (a,b)
        real(wp), intent(out)  :: Elist(Limit) !! vector of dimension at least `limit`, the first
                                               !! `last` elements of which are the moduli of the
                                               !! absolute error estimates on the subintervals
        real(wp), intent(out)  :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                               !! `last` elements of which are the integral
                                               !! approximations on the subintervals
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if  `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(in) :: Points(Npts2) !! vector of dimension npts2, the first (npts2-2)
                                              !! elements of which are the user provided break
                                              !! points. if these points do not constitute an
                                              !! ascending sequence there will be an automatic
                                              !! sorting.
        real(wp), intent(out) :: Pts(Npts2) !! vector of dimension at least npts2, containing the
                                            !! integration limits and the break points of the
                                            !! interval in ascending sequence.
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine.
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties. if
                                    !!   the position of a local difficulty can be
                                    !!   determined (i.e. singularity,
                                    !!   discontinuity within the interval), it
                                    !!   should be supplied to the routine as an
                                    !!   element of the vector points. if necessary
                                    !!   an appropriate special-purpose integrator
                                    !!   must be used, which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table. it is presumed that
                                    !!   the requested tolerance cannot be
                                    !!   achieved, and that the returned result is
                                    !!   the best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier>0.
                                    !! * ier = 6 the input is invalid because
                                    !!   `npts2<2` or
                                    !!   break points are specified outside
                                    !!   the integration range or
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!   or `limit<npts2`.
                                    !!   `result`, `abserr`, `neval`, `last`, `rlist(1)`,
                                    !!   and elist(1) are set to zero. alist(1) and
                                    !!   blist(1) are set to `a` and `b` respectively.
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! elements of which are pointers to the
                                            !! error estimates over the subintervals,
                                            !! such that `elist(iord(1)), ..., elist(iord(k))`
                                            !! form a decreasing sequence, with `k = last`
                                            !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                            !! otherwise
        integer, intent(out) :: Last !! number of subintervals actually produced in the
                                     !! subdivisions process
        integer, intent(in) :: Limit !! gives an upper bound on the number of subintervals
                                     !! in the partition of `(a,b)`, `limit>=npts2`
                                     !! if `limit<npts2`, the routine will end with
                                     !! ier = 6.
        integer, intent(in) :: Npts2 !! number equal to two more than the number of
                                     !! user-supplied break points within the integration
                                     !! range, `npts2>=2`.
                                     !! if `npts2<2`, the routine will end with ier = 6.
        integer, intent(out) :: Ndin(Npts2) !! vector of dimension at least npts2, after first
                                            !! integration over the intervals `(pts(i)),pts(i+1)`,
                                            !! `i = 0,1, ..., npts2-2`, the error estimates over
                                            !! some of the intervals may have been increased
                                            !! artificially, in order to put their subdivision
                                            !! forward. if this happens for the subinterval
                                            !! numbered `k`, `ndin(k)` is put to 1, otherwise
                                            !! `ndin(k) = 0`.
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Level(Limit) !! vector of dimension at least `limit`, containing the
                                             !! subdivision levels of the subinterval, i.e. if
                                             !! `(aa,bb)` is a subinterval of `(p1,p2)` where `p1` as
                                             !! well as `p2` is a user-provided break point or
                                             !! integration limit, then `(aa,bb)` has level `l` if
                                             !! `abs(bb-aa) = abs(p2-p1)*2**(-l)`.

        real(wp) :: a, abseps, b, correc, defabs, &
                    dres, ertest, resa, reseps, Result, &
                    res3la(3), sign, temp, resabs
        integer :: i, id, ierro, ind1, ind2, ip1, iroff1, &
                   iroff2, iroff3, j, jlow, jupbnd, k, ksgn, ktmin, &
                   levcur, levmax, nint, nintp1, npts, nrmax
        real(wp) :: area1, a1, b1, defab1, error1 !! variable for the left subinterval
        real(wp) :: area2, a2, b2, defab2, error2 !! variable for the right subinterval
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: rlist2(limexp + 2) !! array of dimension at least `limexp+2`
                                       !! containing the part of the epsilon table which
                                       !! is still needed for further computations.
        real(wp) :: erlast !! error on the interval currently subdivided
                           !! (before that subdivision has taken place)
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: erlarg !! sum of the errors over the intervals larger
                           !! than the smallest interval considered up to now
        real(wp) :: errmax !! `elist(maxerr)`
        logical :: extrap !! logical variable denoting that the routine
                          !! is attempting to perform extrapolation. i.e.
                          !! before subdividing the smallest interval we
                          !! try to decrease the value of `erlarg`.
        logical :: noext !! logical variable denoting that extrapolation is
                         !! no longer allowed (true-value)
        integer :: maxerr !! pointer to the interval with largest error estimate
        integer :: nres !! number of calls to the extrapolation routine
        integer :: numrl2 !! number of elements in `rlist2`. if an appropriate
                          !! approximation to the compounded integral has
                          !! been obtained, it is put in `rlist2(numrl2)` after
                          !! `numrl2` has been increased by one.

        ! test on validity of parameters

        Ier = 0
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        Alist(1) = a
        Blist(1) = b
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        Level(1) = 0
        npts = Npts2 - 2
        if (Npts2 < 2 .or. Limit <= npts .or. &
            (Epsabs <= 0.0_wp .and. Epsrel < max(50.0_wp*epmach, 0.5e-28_wp))) &
            Ier = 6
        if (Ier == 6) return

        ! if any break points are provided, sort them into an
        ! ascending sequence.

        sign = 1.0_wp
        if (a > b) sign = -1.0_wp
        Pts(1) = min(a, b)
        if (npts /= 0) then
            do i = 1, npts
                Pts(i + 1) = Points(i)
            end do
        end if
        Pts(npts + 2) = max(a, b)
        nint = npts + 1
        a1 = Pts(1)
        if (npts /= 0) then
            nintp1 = nint + 1
            do i = 1, nint
                ip1 = i + 1
                do j = ip1, nintp1
                    if (Pts(i) > Pts(j)) then
                        temp = Pts(i)
                        Pts(i) = Pts(j)
                        Pts(j) = temp
                    end if
                end do
            end do
            if (Pts(1) /= min(a, b) .or. Pts(nintp1) /= max(a, b)) Ier = 6
            if (Ier == 6) return
        end if

        main : block

            ! compute first integral and error approximations.

            resabs = 0.0_wp
            do i = 1, nint
                b1 = Pts(i + 1)
                call dqk21(f, a1, b1, area1, error1, defabs, resa)
                Abserr = Abserr + error1
                Result = Result + area1
                Ndin(i) = 0
                if (error1 == resa .and. error1 /= 0.0_wp) Ndin(i) = 1
                resabs = resabs + defabs
                Level(i) = 0
                Elist(i) = error1
                Alist(i) = a1
                Blist(i) = b1
                Rlist(i) = area1
                Iord(i) = i
                a1 = b1
            end do
            errsum = 0.0_wp
            do i = 1, nint
                if (Ndin(i) == 1) Elist(i) = Abserr
                errsum = errsum + Elist(i)
            end do

            ! test on accuracy.

            Last = nint
            Neval = 21*nint
            dres = abs(Result)
            errbnd = max(Epsabs, Epsrel*dres)
            if (Abserr <= 100.0_wp*epmach*resabs .and. Abserr > errbnd) Ier = 2
            if (nint /= 1) then
            do i = 1, npts
                jlow = i + 1
                ind1 = Iord(i)
                do j = jlow, nint
                    ind2 = Iord(j)
                    if (Elist(ind1) <= Elist(ind2)) then
                        ind1 = ind2
                        k = j
                    end if
                end do
                if (ind1 /= Iord(i)) then
                    Iord(k) = Iord(i)
                    Iord(i) = ind1
                end if
            end do
            if (Limit < Npts2) Ier = 1
            end if
            if (Ier /= 0 .or. Abserr <= errbnd) exit main

            ! initialization

            rlist2(1) = Result
            maxerr = Iord(1)
            errmax = Elist(maxerr)
            area = Result
            nrmax = 1
            nres = 0
            numrl2 = 1
            ktmin = 0
            extrap = .false.
            noext = .false.
            erlarg = errsum
            ertest = errbnd
            levmax = 1
            iroff1 = 0
            iroff2 = 0
            iroff3 = 0
            ierro = 0
            Abserr = oflow
            ksgn = -1
            if (dres >= (1.0_wp - 50.0_wp*epmach)*resabs) ksgn = 1

            ! main do-loop

            loop: do Last = Npts2, Limit

                ! bisect the subinterval with the nrmax-th largest error
                ! estimate.

                levcur = Level(maxerr) + 1
                a1 = Alist(maxerr)
                b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                a2 = b1
                b2 = Blist(maxerr)
                erlast = errmax
                call dqk21(f, a1, b1, area1, error1, resa, defab1)
                call dqk21(f, a2, b2, area2, error2, resa, defab2)

                ! improve previous approximations to integral
                ! and error and test for accuracy.

                Neval = Neval + 42
                area12 = area1 + area2
                erro12 = error1 + error2
                errsum = errsum + erro12 - errmax
                area = area + area12 - Rlist(maxerr)
                if (defab1 /= error1 .and. defab2 /= error2) then
                    if (abs(Rlist(maxerr) - area12) <= 0.1e-4_wp*abs(area12) .and. &
                        erro12 >= 0.99_wp*errmax) then
                        if (extrap) iroff2 = iroff2 + 1
                        if (.not. extrap) iroff1 = iroff1 + 1
                    end if
                    if (Last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
                end if
                Level(maxerr) = levcur
                Level(Last) = levcur
                Rlist(maxerr) = area1
                Rlist(Last) = area2
                errbnd = max(Epsabs, Epsrel*abs(area))

                ! test for roundoff error and eventually set error flag.

                if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) Ier = 2
                if (iroff2 >= 5) ierro = 3

                ! set error flag in the case that the number of
                ! subintervals equals limit.

                if (Last == Limit) Ier = 1

                ! set error flag in the case of bad integrand behaviour
                ! at a point of the integration range

                if (max(abs(a1), abs(b2)) <= (1.0_wp + 100.0_wp*epmach) &
                    *(abs(a2) + 1000.0_wp*uflow)) Ier = 4

                ! append the newly-created intervals to the list.

                if (error2 > error1) then
                    Alist(maxerr) = a2
                    Alist(Last) = a1
                    Blist(Last) = b1
                    Rlist(maxerr) = area2
                    Rlist(Last) = area1
                    Elist(maxerr) = error2
                    Elist(Last) = error1
                else
                    Alist(Last) = a2
                    Blist(maxerr) = b1
                    Blist(Last) = b2
                    Elist(maxerr) = error1
                    Elist(Last) = error2
                end if

                ! call subroutine dqpsrt to maintain the descending ordering
                ! in the list of error estimates and select the subinterval
                ! with nrmax-th largest error estimate (to be bisected next).

                call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                ! ***jump out of do-loop
                if (errsum <= errbnd) then
                    ! compute global integral sum.
                    Result = sum(Rlist(1:Last))
                    Abserr = errsum
                    exit main
                end if
                ! ***jump out of do-loop
                if (Ier /= 0) exit loop
                if (.not. (noext)) then
                    erlarg = erlarg - erlast
                    if (levcur + 1 <= levmax) erlarg = erlarg + erro12
                    if (.not. (extrap)) then
                        ! test whether the interval to be bisected next is the
                        ! smallest interval.
                        if (Level(maxerr) + 1 <= levmax) cycle loop
                        extrap = .true.
                        nrmax = 2
                    end if
                    if (ierro /= 3 .and. erlarg > ertest) then
                        ! the smallest interval has the largest error.
                        ! before bisecting decrease the sum of the errors over
                        ! the larger intervals (erlarg) and perform extrapolation.
                        id = nrmax
                        jupbnd = Last
                        if (Last > (2 + Limit/2)) jupbnd = Limit + 3 - Last
                        do k = id, jupbnd
                            maxerr = Iord(nrmax)
                            errmax = Elist(maxerr)
                            ! ***jump out of do-loop
                            if (Level(maxerr) + 1 <= levmax) cycle loop
                            nrmax = nrmax + 1
                        end do
                    end if

                    ! perform extrapolation.

                    numrl2 = numrl2 + 1
                    rlist2(numrl2) = area
                    if (numrl2 > 2) then
                        call dqelg(numrl2, rlist2, reseps, abseps, res3la, nres)
                        ktmin = ktmin + 1
                        if (ktmin > 5 .and. Abserr < 0.1e-02_wp*errsum) Ier = 5
                        if (abseps < Abserr) then
                            ktmin = 0
                            Abserr = abseps
                            Result = reseps
                            correc = erlarg
                            ertest = max(Epsabs, Epsrel*abs(reseps))
                            ! ***jump out of do-loop
                            if (Abserr < ertest) exit loop
                        end if
                        ! prepare bisection of the smallest interval.
                        if (numrl2 == 1) noext = .true.
                        if (Ier >= 5) exit loop
                    end if
                    maxerr = Iord(1)
                    errmax = Elist(maxerr)
                    nrmax = 1
                    extrap = .false.
                    levmax = levmax + 1
                    erlarg = errsum
                end if

            end do loop

            ! set the final result.

            if (Abserr /= oflow) then
                if ((Ier + ierro) /= 0) then
                    if (ierro == 3) Abserr = Abserr + correc
                    if (Ier == 0) Ier = 3
                    if (Result == 0.0_wp .or. area == 0.0_wp) then
                        if (Abserr > errsum) then
                            ! compute global integral sum.
                            Result = sum(Rlist(1:Last))
                            Abserr = errsum
                            exit main
                        end if
                        if (area == 0.0_wp) exit main
                    elseif (Abserr/abs(Result) > errsum/abs(area)) then
                        ! compute global integral sum.
                        Result = sum(Rlist(1:Last))
                        Abserr = errsum
                        exit main
                    end if
                end if

                ! test on divergence.

                if (ksgn /= (-1) .or. max(abs(Result), abs(area)) &
                    > resabs*0.01_wp) then
                    if (0.01_wp > (Result/area) .or. (Result/area) > 100.0_wp .or. &
                        errsum > abs(area)) Ier = 6
                end if

            else
                ! compute global integral sum.
                Result = sum(Rlist(1:Last))
                Abserr = errsum
            end if

        end block main

        if (Ier > 2) Ier = Ier - 1
        Result = Result*sign

    end subroutine dqagpe
!********************************************************************************

!********************************************************************************
!>
!  1D globally adaptive integrator using interval subdivision and extrapolation
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqags(f, a, b, Epsabs, Epsrel, Result, Abserr, Neval, Ier, &
                     Limit, Lenw, Last, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand
                             !! function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more sub-
                                    !!   divisions by increasing the value of limit
                                    !!   (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties. if
                                    !!   the position of a local difficulty can be
                                    !!   determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling the
                                    !!   integrator on the subranges. if possible,
                                    !!   an appropriate special-purpose integrator
                                    !!   should be used, which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is detected,
                                    !!   which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour
                                    !!   occurs at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table. it is presumed that
                                    !!   the requested tolerance cannot be
                                    !!   achieved, and that the returned result is
                                    !!   the best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `(epsabs<=0` and
                                    !!   `epsrel<max(50*rel.mach.acc.,0.5e-28)`
                                    !!   or `limit<1` or `lenw<limit*4`.
                                    !!   `result`, `abserr`, `neval`, `last` are set to
                                    !!   zero. except when limit or lenw is invalid,
                                    !!   `iwork(1), work(limit*2+1)` and
                                    !!   `work(limit*3+1)` are set to zero, `work(1)`
                                    !!   is set to `a` and `work(limit+1)` to `b`.
        integer, intent(in) :: Limit !! dimensioning parameter for `iwork`.
                                     !! `limit` determines the maximum number of subintervals
                                     !! in the partition of the given integration interval
                                     !! `(a,b)`, `limit>=1`.
                                     !! if `limit<1`, the routine will end with ier = 6.
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`.
                                    !! `lenw` must be at least `limit*4`.
                                    !! if `lenw<limit*4`, the routine will end
                                    !! with ier = 6.
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, determines the
                                     !! number of significant elements actually in the `work`
                                     !! arrays.
        integer :: Iwork(Limit) !! vector of dimension at least `limit`, the first `k`
                                !! elements of which contain pointers
                                !! to the error estimates over the subintervals
                                !! such that `work(limit*3+iwork(1)),...,work(limit*3+iwork(k))`
                                !! form a decreasing sequence, with `k = last` if `last<=(limit/2+2)`,
                                !! and `k = limit+1-last` otherwise
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`.
                               !! on return:
                               !!
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end-points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !! * `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end-points,
                               !! * `work(limit*2+1), ..., work(limit*2+last)` contain
                               !!   the integral approximations over the subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3+last)`
                               !!   contain the error estimates.

        integer :: lvl, l1, l2, l3

        ! check validity of limit and lenw.

        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Limit >= 1 .and. Lenw >= Limit*4) then

            ! prepare call for dqagse.
            l1 = Limit + 1
            l2 = Limit + l1
            l3 = Limit + l2

            call dqagse(f, a, b, Epsabs, Epsrel, Limit, Result, Abserr, Neval, Ier, &
                        Work(1), Work(l1), Work(l2), Work(l3), Iwork, Last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqags', Ier, lvl)

    end subroutine dqags
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqags]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqagse(f, a, b, Epsabs, Epsrel, Limit, Result, Abserr, Neval, &
                      Ier, Alist, Blist, Rlist, Elist, Iord, Last)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand
                             !! function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        integer, intent(in) :: Limit !! gives an upperbound on the number of subintervals
                                     !! in the partition of `(a,b)`
        real(wp), intent(out) :: Result !! approximation to the integral
        integer, intent(out) :: Neval !! number of integrand evaluations
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left end points
                                              !! of the subintervals in the partition of the
                                              !! given integration range (a,b)
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right end points
                                              !! of the subintervals in the partition of the given
                                              !! integration range (a,b)
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the moduli of the
                                              !! absolute error estimates on the subintervals
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the integral
                                              !! approximations on the subintervals
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more sub-
                                    !!   divisions by increasing the value of limit
                                    !!   (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties. if
                                    !!   the position of a local difficulty can be
                                    !!   determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling the
                                    !!   integrator on the subranges. if possible,
                                    !!   an appropriate special-purpose integrator
                                    !!   should be used, which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour
                                    !!   occurs at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table.
                                    !!   it is presumed that the requested
                                    !!   tolerance cannot be achieved, and that the
                                    !!   returned result is the best which can be
                                    !!   obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `epsabs<=0` and
                                    !!   `epsrel<max(50*rel.mach.acc.,0.5e-28)`.
                                    !!   `result`, `abserr`, `neval`, `last`, `rlist(1)`,
                                    !!   `iord(1)` and `elist(1)` are set to zero.
                                    !!   `alist(1)` and `blist(1)` are set to a and b
                                    !!   respectively.
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! elements of which are pointers to the
                                            !! error estimates over the subintervals,
                                            !! such that `elist(iord(1)), ..., elist(iord(k))`
                                            !! form a decreasing sequence, with `k = last`
                                            !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                            !! otherwise
        integer, intent(out) :: Last !! number of subintervals actually produced in the
                                     !! subdivision process

        real(wp) :: abseps, correc, defabs, dres, &
                    ertest, resabs, reseps, res3la(3)
        integer :: id, ierro, iroff1, iroff2, iroff3, &
                   jupbnd, k, ksgn, ktmin, nrmax
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: area1, a1, b1, defab1, error1 !! variable for the left interval
        real(wp) :: area2, a2, b2, defab2, error2 !! variable for the right interval
        real(wp) :: rlist2(limexp + 2) !! array of dimension at least `limexp+2` containing
                                       !! the part of the epsilon table which is still
                                       !! needed for further computations.
        integer :: maxerr !! pointer to the interval with largest error estimate
        integer :: nres !! number of calls to the extrapolation routine
        integer :: numrl2 !! number of elements currently in `rlist2`. if an
                          !! appropriate approximation to the compounded
                          !! integral has been obtained it is put in
                          !! `rlist2(numrl2)` after `numrl2` has been increased
                          !! by one.
        real(wp) :: errmax !! elist(maxerr)
        real(wp) :: erlast !! error on the interval currently subdivided
                           !! (before that subdivision has taken place)
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        real(wp) :: small !! length of the smallest interval considered up
                          !! to now, multiplied by 1.5
        real(wp) :: erlarg !! sum of the errors over the intervals larger
                           !! than the smallest interval considered up to now
        logical :: extrap !! logical variable denoting that the routine is
                          !! attempting to perform extrapolation i.e. before
                          !! subdividing the smallest interval we try to
                          !! decrease the value of `erlarg`.
        logical :: noext !! logical variable denoting that extrapolation
                         !! is no longer allowed (true value)

        ! test on validity of parameters

        Ier = 0
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        Alist(1) = a
        Blist(1) = b
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        if (Epsabs <= 0.0_wp .and. Epsrel < max(50.0_wp*epmach, 0.5e-28_wp)) Ier = 6
        if (Ier == 6) return

        main : block

            ! first approximation to the integral

            ierro = 0
            call dqk21(f, a, b, Result, Abserr, defabs, resabs)

            ! test on accuracy.

            dres = abs(Result)
            errbnd = max(Epsabs, Epsrel*dres)
            Last = 1
            Rlist(1) = Result
            Elist(1) = Abserr
            Iord(1) = 1
            if (Abserr <= 100.0_wp*epmach*defabs .and. Abserr > errbnd) &
                Ier = 2
            if (Limit == 1) Ier = 1
            if (Ier /= 0 .or. (Abserr <= errbnd .and. Abserr /= resabs) .or. &
                Abserr == 0.0_wp) then
                Neval = 42*Last - 21
                return
            end if

            ! initialization

            rlist2(1) = Result
            errmax = Abserr
            maxerr = 1
            area = Result
            errsum = Abserr
            Abserr = oflow
            nrmax = 1
            nres = 0
            numrl2 = 2
            ktmin = 0
            extrap = .false.
            noext = .false.
            iroff1 = 0
            iroff2 = 0
            iroff3 = 0
            ksgn = -1
            if (dres >= (1.0_wp - 50.0_wp*epmach)*defabs) ksgn = 1

            ! main do-loop

            loop: do Last = 2, Limit

                ! bisect the subinterval with the nrmax-th largest error
                ! estimate.

                a1 = Alist(maxerr)
                b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                a2 = b1
                b2 = Blist(maxerr)
                erlast = errmax
                call dqk21(f, a1, b1, area1, error1, resabs, defab1)
                call dqk21(f, a2, b2, area2, error2, resabs, defab2)

                ! improve previous approximations to integral
                ! and error and test for accuracy.

                area12 = area1 + area2
                erro12 = error1 + error2
                errsum = errsum + erro12 - errmax
                area = area + area12 - Rlist(maxerr)
                if (defab1 /= error1 .and. defab2 /= error2) then
                    if (abs(Rlist(maxerr) - area12) <= 0.1e-4_wp*abs(area12) &
                        .and. erro12 >= 0.99_wp*errmax) then
                        if (extrap) iroff2 = iroff2 + 1
                        if (.not. extrap) iroff1 = iroff1 + 1
                    end if
                    if (Last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
                end if
                Rlist(maxerr) = area1
                Rlist(Last) = area2
                errbnd = max(Epsabs, Epsrel*abs(area))

                ! test for roundoff error and eventually set error flag.

                if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) Ier = 2
                if (iroff2 >= 5) ierro = 3

                ! set error flag in the case that the number of subintervals
                ! equals limit.

                if (Last == Limit) Ier = 1

                ! set error flag in the case of bad integrand behaviour
                ! at a point of the integration range.

                if (max(abs(a1), abs(b2)) <= (1.0_wp + 100.0_wp*epmach) &
                    *(abs(a2) + 1000.0_wp*uflow)) Ier = 4

                ! append the newly-created intervals to the list.

                if (error2 > error1) then
                    Alist(maxerr) = a2
                    Alist(Last) = a1
                    Blist(Last) = b1
                    Rlist(maxerr) = area2
                    Rlist(Last) = area1
                    Elist(maxerr) = error2
                    Elist(Last) = error1
                else
                    Alist(Last) = a2
                    Blist(maxerr) = b1
                    Blist(Last) = b2
                    Elist(maxerr) = error1
                    Elist(Last) = error2
                end if

                ! call subroutine dqpsrt to maintain the descending ordering
                ! in the list of error estimates and select the subinterval
                ! with nrmax-th largest error estimate (to be bisected next).

                call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                ! ***jump out of do-loop
                if (errsum <= errbnd) exit main
                ! ***jump out of do-loop
                if (Ier /= 0) exit loop
                if (Last == 2) then
                    small = abs(b - a)*0.375_wp
                    erlarg = errsum
                    ertest = errbnd
                    rlist2(2) = area
                elseif (.not. (noext)) then
                    erlarg = erlarg - erlast
                    if (abs(b1 - a1) > small) erlarg = erlarg + erro12
                    if (.not. (extrap)) then
                        ! test whether the interval to be bisected next is the
                        ! smallest interval.
                        if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle loop
                        extrap = .true.
                        nrmax = 2
                    end if
                    if (ierro /= 3 .and. erlarg > ertest) then
                        ! the smallest interval has the largest error.
                        ! before bisecting decrease the sum of the errors over the
                        ! larger intervals (erlarg) and perform extrapolation.
                        id = nrmax
                        jupbnd = Last
                        if (Last > (2 + Limit/2)) jupbnd = Limit + 3 - Last
                        do k = id, jupbnd
                            maxerr = Iord(nrmax)
                            errmax = Elist(maxerr)
                            ! ***jump out of do-loop
                            if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle loop
                            nrmax = nrmax + 1
                        end do
                    end if

                    ! perform extrapolation.

                    numrl2 = numrl2 + 1
                    rlist2(numrl2) = area
                    call dqelg(numrl2, rlist2, reseps, abseps, res3la, nres)
                    ktmin = ktmin + 1
                    if (ktmin > 5 .and. Abserr < 0.1e-02_wp*errsum) Ier = 5
                    if (abseps < Abserr) then
                        ktmin = 0
                        Abserr = abseps
                        Result = reseps
                        correc = erlarg
                        ertest = max(Epsabs, Epsrel*abs(reseps))
                        ! ***jump out of do-loop
                        if (Abserr <= ertest) exit loop
                    end if

                    ! prepare bisection of the smallest interval.

                    if (numrl2 == 1) noext = .true.
                    if (Ier == 5) exit loop
                    maxerr = Iord(1)
                    errmax = Elist(maxerr)
                    nrmax = 1
                    extrap = .false.
                    small = small*0.5_wp
                    erlarg = errsum
                end if
            end do loop

            ! set final result and error estimate.

            if (Abserr /= oflow) then
                if (Ier + ierro /= 0) then
                    if (ierro == 3) Abserr = Abserr + correc
                    if (Ier == 0) Ier = 3
                    if (Result == 0.0_wp .or. area == 0.0_wp) then
                        if (Abserr > errsum) exit main
                        if (area == 0.0_wp) then
                            if (Ier > 2) Ier = Ier - 1
                            Neval = 42*Last - 21
                            return
                        end if
                    elseif (Abserr/abs(Result) > errsum/abs(area)) then
                        exit main
                    end if
                end if

                ! test on divergence.

                if (ksgn /= (-1) .or. max(abs(Result), abs(area)) &
                    > defabs*0.01_wp) then
                    if (0.01_wp > (Result/area) .or. (Result/area) &
                        > 100.0_wp .or. errsum > abs(area)) Ier = 6
                end if
                if (Ier > 2) Ier = Ier - 1
                Neval = 42*Last - 21
                return
            end if

        end block main

        ! compute global integral sum.

        Result = sum(Rlist(1:Last))
        Abserr = errsum
        if (Ier > 2) Ier = Ier - 1
        Neval = 42*Last - 21

    end subroutine dqagse
!********************************************************************************

!********************************************************************************
!>
!  compute Cauchy principal value of `f(x)/(x-c)` over a finite interval
!
!  the routine calculates an approximation result to a
!  cauchy principal value i = integral of `f*w` over `(a,b)`
!  `(w(x) = 1/((x-c), c/=a, c/=b)`, hopefully satisfying
!  following claim for accuracy
!  `abs(i-result)<=max(epsabe,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawc(f, a, b, c, Epsabs, Epsrel, Result, Abserr, Neval, Ier, &
                     Limit, Lenw, Last, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! under limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: c !! parameter in the weight function, `c/=a`, `c/=b`.
                                  !! if `c = a` or `c = b`, the routine will end with
                                  !! ier = 6 .
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate or the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more sub-
                                    !!   divisions by increasing the value of limit
                                    !!   (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand in order to
                                    !!   determine the integration difficulties.
                                    !!   if the position of a local difficulty
                                    !!   can be determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling
                                    !!   appropriate integrators on the subranges.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `c = a` or `c = b` or
                                    !!   (`epsabs<=0` and `epsrel<max(50*rel.mach.acc.,0.5e-28)`)
                                    !!   or `limit<1` or `lenw<limit*4`.
                                    !!   `esult`, `abserr`, `neval`, `last` are set to
                                    !!   zero. except when `lenw` or `limit` is invalid,
                                    !!   `iwork(1)`, `work(limit*2+1)` and
                                    !!   `work(limit*3+1)` are set to zero, `work(1)`
                                    !!   is set to a and `work(limit+1)` to `b`.
        integer, intent(in) :: Limit !! dimensioning parameter for `iwork`.
                                     !! `limit` determines the maximum number of subintervals
                                     !! in the partition of the given integration interval
                                     !! `(a,b)`, `limit>=1`.
                                     !! if `limit<1`, the routine will end with ier = 6.
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`.
                                    !! `lenw` must be at least `limit*4`.
                                    !! if `lenw<limit*4`, the routine will end with
                                    !! ier = 6.
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, which
                                     !! determines the number of significant elements
                                     !! actually in the work arrays.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`.
                               !! on return:
                               !!
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !! * `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end points,
                               !! * `work(limit*2+1), ..., work(limit*2+last)` contain
                               !!   the integral approximations over the subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3+last)`
                               !!   contain the error estimates.
        integer :: Iwork(Limit) !! vector of dimension at least `limit`, the first `k`
                                !! elements of which contain pointers
                                !! to the error estimates over the subintervals,
                                !! such that `work(limit*3+iwork(1)),...,work(limit*3+iwork(k))`
                                !! form a decreasing sequence, with `k = last` if
                                !! `last<=(limit/2+2)`, and `k = limit+1-last` otherwise

        integer :: lvl, l1, l2, l3

        ! check validity of limit and lenw.
        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Limit >= 1 .and. Lenw >= Limit*4) then

            ! prepare call for dqawce.
            l1 = Limit + 1
            l2 = Limit + l1
            l3 = Limit + l2
            call dqawce(f, a, b, c, Epsabs, Epsrel, Limit, Result, Abserr, Neval, &
                        Ier, Work(1), Work(l1), Work(l2), Work(l3), Iwork, Last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqawc', Ier, lvl)

    end subroutine dqawc
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqawc]] but provides more information and control
!
!  the routine calculates an approximation result to a
!  cauchy principal value i = integral of `f*w` over `(a,b)`
!  `(w(x) = 1/(x-c), (c/=a, c/=b)`, hopefully satisfying
!  following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawce(f, a, b, c, Epsabs, Epsrel, Limit, Result, Abserr, Neval, &
                      Ier, Alist, Blist, Rlist, Elist, Iord, Last)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        integer, intent(in) :: Limit !! gives an upper bound on the number of subintervals
                                     !! in the partition of `(a,b)`, `limit>=1`
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more sub-
                                    !!   divisions by increasing the value of
                                    !!   limit. however, if this yields no
                                    !!   improvement it is advised to analyze the
                                    !!   the integrand, in order to determine the
                                    !!   the integration difficulties. if the
                                    !!   position of a local difficulty can be
                                    !!   determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling
                                    !!   appropriate integrators on the subranges.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !! * ier = 3 extremely bad integrand behaviour
                                    !!   occurs at some interior points of
                                    !!   the integration interval.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `c = a` or `c = b` or
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!   or `limit<1`.
                                    !!   `result`, `abserr`, `neval`, `rlist(1)`, `elist(1)`,
                                    !!   `iord(1)` and `last` are set to zero. `alist(1)`
                                    !!   and `blist(1)` are set to `a` and `b`
                                    !!   respectively.
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the integral
                                              !! approximations on the subintervals
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension `limit`, the first `last`
                                              !! elements of which are the moduli of the absolute
                                              !! error estimates on the subintervals
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! elements of which are pointers to the error
                                            !! estimates over the subintervals, so that
                                            !! `elist(iord(1)), ..., elist(iord(k))` with `k = last`
                                            !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                            !! otherwise, form a decreasing sequence
        integer, intent(out) :: Last !! number of subintervals actually produced in
                                     !! the subdivision process

        real(wp) :: aa, bb, c
        integer :: iroff1, iroff2, k, krule, nev, nrmax
        real(wp) :: area1, a1, b1, error1 !! variable for the left subinterval
        real(wp) :: area2, a2, b2, error2 !! variable for the right subinterval
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: errmax !! elist(maxerr)
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        integer :: maxerr !! pointer to the interval with largest error estimate

        ! test on validity of parameters

        Ier = 6
        Neval = 0
        Last = 0
        Alist(1) = a
        Blist(1) = b
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (.not. (c == a .or. c == b .or. (Epsabs <= 0.0_wp .and. Epsrel < max &
                                            (50.0_wp*epmach, 0.5e-28_wp)))) then

            ! first approximation to the integral

            aa = a
            bb = b
            if (a > b) then
                aa = b
                bb = a
            end if
            Ier = 0
            krule = 1
            call dqc25c(f, aa, bb, c, Result, Abserr, krule, Neval)
            Last = 1
            Rlist(1) = Result
            Elist(1) = Abserr
            Iord(1) = 1
            Alist(1) = a
            Blist(1) = b

            ! test on accuracy

            errbnd = max(Epsabs, Epsrel*abs(Result))
            if (Limit == 1) Ier = 1
            if (Abserr >= min(0.01_wp*abs(Result), errbnd) .and. Ier /= 1) then

                ! initialization

                Alist(1) = aa
                Blist(1) = bb
                Rlist(1) = Result
                errmax = Abserr
                maxerr = 1
                area = Result
                errsum = Abserr
                nrmax = 1
                iroff1 = 0
                iroff2 = 0

                ! main do-loop

                do Last = 2, Limit

                    ! bisect the subinterval with nrmax-th largest
                    ! error estimate.

                    a1 = Alist(maxerr)
                    b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                    b2 = Blist(maxerr)
                    if (c <= b1 .and. c > a1) b1 = 0.5_wp*(c + b2)
                    if (c > b1 .and. c < b2) b1 = 0.5_wp*(a1 + c)
                    a2 = b1
                    krule = 2
                    call dqc25c(f, a1, b1, c, area1, error1, krule, nev)
                    Neval = Neval + nev
                    call dqc25c(f, a2, b2, c, area2, error2, krule, nev)
                    Neval = Neval + nev

                    ! improve previous approximations to integral
                    ! and error and test for accuracy.

                    area12 = area1 + area2
                    erro12 = error1 + error2
                    errsum = errsum + erro12 - errmax
                    area = area + area12 - Rlist(maxerr)
                    if (abs(Rlist(maxerr) - area12) < 0.1e-4_wp*abs(area12) &
                        .and. erro12 >= 0.99_wp*errmax .and. krule == 0) &
                        iroff1 = iroff1 + 1
                    if (Last > 10 .and. erro12 > errmax .and. krule == 0) &
                        iroff2 = iroff2 + 1
                    Rlist(maxerr) = area1
                    Rlist(Last) = area2
                    errbnd = max(Epsabs, Epsrel*abs(area))
                    if (errsum > errbnd) then

                        ! test for roundoff error and eventually set error flag.

                        if (iroff1 >= 6 .and. iroff2 > 20) Ier = 2

                        ! set error flag in the case that number of interval
                        ! bisections exceeds limit.

                        if (Last == Limit) Ier = 1

                        ! set error flag in the case of bad integrand behaviour
                        ! at a point of the integration range.

                        if (max(abs(a1), abs(b2)) &
                            <= (1.0_wp + 100.0_wp*epmach) &
                            *(abs(a2) + 1000.0_wp*uflow)) Ier = 3
                    end if

                    ! append the newly-created intervals to the list.

                    if (error2 > error1) then
                        Alist(maxerr) = a2
                        Alist(Last) = a1
                        Blist(Last) = b1
                        Rlist(maxerr) = area2
                        Rlist(Last) = area1
                        Elist(maxerr) = error2
                        Elist(Last) = error1
                    else
                        Alist(Last) = a2
                        Blist(maxerr) = b1
                        Blist(Last) = b2
                        Elist(maxerr) = error1
                        Elist(Last) = error2
                    end if

                    ! call subroutine dqpsrt to maintain the descending ordering
                    ! in the list of error estimates and select the subinterval
                    ! with nrmax-th largest error estimate (to be bisected next).

                    call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                    ! ***jump out of do-loop
                    if (Ier /= 0 .or. errsum <= errbnd) exit
                end do

                ! compute final result.
                Result = 0.0_wp
                do k = 1, Last
                    Result = Result + Rlist(k)
                end do
                Abserr = errsum
            end if
            if (aa == b) Result = -Result
        end if

    end subroutine dqawce
!********************************************************************************

!********************************************************************************
!>
!  Fourier sine/cosine transform for user supplied interval `a` to `infinity`
!
!  the routine calculates an approximation result to a given
!  fourier integral i=integral of `f(x)*w(x)` over `(a,infinity)`
!  where `w(x) = cos(omega*x)` or `w(x) = sin(omega*x)`.
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=epsabs`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawf(f, a, Omega, Integr, Epsabs, Result, Abserr, Neval, Ier, &
                     Limlst, Lst, Leniw, Maxp1, Lenw, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: Omega !! parameter in the integrand weight function
        integer, intent(in) :: Integr !! indicates which of the weight functions is used:
                                      !!
                                      !! * integr = 1 `w(x) = cos(omega*x)`
                                      !! * integr = 2 `w(x) = sin(omega*x)`
                                      !!
                                      !! if `integr/=1 .and. integr/=2`, the routine
                                      !! will end with ier = 6.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested, `epsabs>0`.
                                       !! if `epsabs<=0`, the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out):: Ier !! * ier = 0 normal and reliable termination of the
                                   !!   routine. it is assumed that the requested
                                   !!   accuracy has been achieved.
                                   !! * ier>0 abnormal termination of the routine.
                                   !!   the estimates for integral and error are
                                   !!   less reliable. it is assumed that the
                                   !!   requested accuracy has not been achieved.
                                   !!
                                   !! error messages:
                                   !!
                                   !! `if omega/=0`:
                                   !!
                                   !! * ier = 1 maximum number of cycles allowed
                                   !!   has been achieved, i.e. of subintervals
                                   !!   `(a+(k-1)c,a+kc)` where
                                   !!   `c = (2*int(abs(omega))+1)*pi/abs(omega)`,
                                   !!   for `k = 1, 2, ..., lst`.
                                   !!   one can allow more cycles by increasing
                                   !!   the value of limlst (and taking the
                                   !!   according dimension adjustments into
                                   !!   account). examine the array iwork which
                                   !!   contains the error flags on the cycles, in
                                   !!   order to look for eventual local
                                   !!   integration difficulties.
                                   !!   if the position of a local difficulty
                                   !!   can be determined (e.g. singularity,
                                   !!   discontinuity within the interval) one
                                   !!   will probably gain from splitting up the
                                   !!   interval at this point and calling
                                   !!   appropriate integrators on the subranges.
                                   !! * ier = 4 the extrapolation table constructed for
                                   !!   convergence acceleration of the series
                                   !!   formed by the integral contributions over
                                   !!   the cycles, does not converge to within
                                   !!   the requested accuracy.
                                   !!   as in the case of ier = 1, it is advised
                                   !!   to examine the array iwork which contains
                                   !!   the error flags on the cycles.
                                   !! * ier = 6 the input is invalid because
                                   !!   `(integr/=1 and integr/=2)` or
                                   !!   `epsabs<=0` or `limlst<1` or
                                   !!   `leniw<(limlst+2)` or `maxp1<1` or
                                   !!   `lenw<(leniw*2+maxp1*25)`.
                                   !!   `result`, `abserr`, `neval`, `lst` are set to
                                   !!   zero.
                                   !! * ier = 7 bad integrand behaviour occurs within
                                   !!   one or more of the cycles. location and
                                   !!   type of the difficulty involved can be
                                   !!   determined from the first `lst` elements of
                                   !!   vector `iwork`.  here `lst` is the number of
                                   !!   cycles actually needed (see below):
                                   !!
                                   !!    * iwork(k) = 1 the maximum number of
                                   !!      subdivisions `(=(leniw-limlst)/2)` has
                                   !!      been achieved on the `k`th cycle.
                                   !!    * iwork(k) = 2 occurrence of roundoff error
                                   !!      is detected and prevents the
                                   !!      tolerance imposed on the `k`th
                                   !!      cycle, from being achieved
                                   !!      on this cycle.
                                   !!    * iwork(k) = 3 extremely bad integrand
                                   !!      behaviour occurs at some
                                   !!      points of the `k`th cycle.
                                   !!    * iwork(k) = 4 the integration procedure
                                   !!      over the `k`th cycle does
                                   !!      not converge (to within the
                                   !!      required accuracy) due to
                                   !!      roundoff in the extrapolation
                                   !!      procedure invoked on this
                                   !!      cycle. it is assumed that the
                                   !!      result on this interval is
                                   !!      the best which can be
                                   !!      obtained.
                                   !!    * iwork(k) = 5 the integral over the `k`th
                                   !!      cycle is probably divergent
                                   !!      or slowly convergent. it must
                                   !!      be noted that divergence can
                                   !!      occur with any other value of
                                   !!      `iwork(k)`.
                                   !!
                                   !! if `omega = 0` and `integr = 1`,
                                   !! the integral is calculated by means of [[dqagie]],
                                   !! and `ier = iwork(1)` (with meaning as described
                                   !! for `iwork(k),k = 1`).
        integer, intent(in) :: Limlst !! limlst gives an upper bound on the number of
                                      !! cycles, `limlst>=3`.
                                      !! if `limlst<3`, the routine will end with ier = 6.
        integer, intent(out) :: Lst !! on return, lst indicates the number of cycles
                                    !! actually needed for the integration.
                                    !! if `omega = 0`, then lst is set to 1.
        integer, intent(in) :: Leniw !! dimensioning parameter for `iwork`. on entry,
                                     !! `(leniw-limlst)/2` equals the maximum number of
                                     !! subintervals allowed in the partition of each
                                     !! cycle, `leniw>=(limlst+2)`.
                                     !! if `leniw<(limlst+2)`, the routine will end with
                                     !! ier = 6.
        integer, intent(in) :: Maxp1 !! maxp1 gives an upper bound on the number of
                                     !! chebyshev moments which can be stored, i.e. for
                                     !! the intervals of lengths `abs(b-a)*2**(-l)`,
                                     !! `l = 0,1, ..., maxp1-2, maxp1>=1`.
                                     !! if `maxp1<1`, the routine will end with ier = 6.
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`.
                                    !! `lenw` must be at least `leniw*2+maxp1*25`.
                                    !! if `lenw<(leniw*2+maxp1*25)`, the routine will
                                    !! end with ier = 6.
        integer :: Iwork(Leniw) !! vector of dimension at least `leniw`
                                !! on return, `iwork(k)` for `k = 1, 2, ..., lst`
                                !! contain the error flags on the cycles.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`
                               !! on return:
                               !!
                               !! * `work(1), ..., work(lst)` contain the integral
                               !!   approximations over the cycles,
                               !! * `work(limlst+1), ..., work(limlst+lst)` contain
                               !!   the error estimates over the cycles.
                               !!
                               !! further elements of work have no specific
                               !! meaning for the user.

        integer :: last, limit, ll2, lvl, l1, l2, l3, l4, l5, l6

        ! check validity of limlst, leniw, maxp1 and lenw.
        Ier = 6
        Neval = 0
        last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Limlst >= 3 .and. Leniw >= (Limlst + 2) .and. Maxp1 >= 1 .and. &
            Lenw >= (Leniw*2 + Maxp1*25)) then

            ! prepare call for dqawfe
            limit = (Leniw - Limlst)/2
            l1 = Limlst + 1
            l2 = Limlst + l1
            l3 = limit + l2
            l4 = limit + l3
            l5 = limit + l4
            l6 = limit + l5
            ll2 = limit + l1
            call dqawfe(f, a, Omega, Integr, Epsabs, Limlst, limit, Maxp1, Result, &
                        Abserr, Neval, Ier, Work(1), Work(l1), Iwork(1), Lst, &
                        Work(l2), Work(l3), Work(l4), Work(l5), Iwork(l1), &
                        Iwork(ll2), Work(l6))

            ! call error handler if necessary
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqawf', Ier, lvl)

    end subroutine dqawf
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqawf]] but provides more information and control
!
!  the routine calculates an approximation result to a
!  given fourier integral
!  i = integral of `f(x)*w(x)` over `(a,infinity)`
!  where `w(x)=cos(omega*x)` or `w(x)=sin(omega*x)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=epsabs`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawfe(f, a, Omega, Integr, Epsabs, Limlst, Limit, Maxp1, &
                      Result, Abserr, Neval, Ier, Rslst, Erlst, Ierlst, Lst, &
                      Alist, Blist, Rlist, Elist, Iord, Nnlog, Chebmo)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: Omega !! parameter in the weight function
        integer, intent(in) :: Integr !! indicates which weight function is used:
                                      !!
                                      !! * integr = 1  `w(x) = cos(omega*x)`
                                      !! * integr = 2  `w(x) = sin(omega*x)`
                                      !!
                                      !! if `integr/=1.and.integr/=2`, the routine will
                                      !! end with ier = 6.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested, `epsabs>0`
                                       !! if `epsabs<=0`, the routine will end with ier = 6.
        integer, intent(in) :: Limlst !! limlst gives an upper bound on the number of
                                      !! cycles, `limlst>=1`.
                                      !! if `limlst<3`, the routine will end with ier = 6.
        integer, intent(in) :: Limit !! gives an upper bound on the number of subintervals
                                     !! allowed in the partition of each cycle, `limit>=1`
                                     !! each cycle, `limit>=1`.
        integer, intent(in) :: Maxp1 !! gives an upper bound on the number of
                                     !! chebyshev moments which can be stored, i.e.
                                     !! for the intervals of lengths
                                     !! `abs(b-a)*2**(-l), `l=0,1, ..., maxp1-2, maxp1>=1``
        real(wp), intent(out) :: Result !! approximation to the integral `x`
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of
                                    !!   the routine. it is assumed that the
                                    !!   requested accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine. the
                                    !!   estimates for integral and error are less
                                    !!   reliable. it is assumed that the requested
                                    !!   accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! if `omega/=0`:
                                    !!
                                    !! * ier = 1 maximum number of `cycles` allowed
                                    !!   has been achieved., i.e. of subintervals
                                    !!   `(a+(k-1)c,a+kc)` where
                                    !!   `c = (2*int(abs(omega))+1)*pi/abs(omega)`,
                                    !!   for `k = 1, 2, ..., lst`.
                                    !!   one can allow more cycles by increasing
                                    !!   the value of limlst (and taking the
                                    !!   according dimension adjustments into
                                    !!   account).
                                    !!   examine the array `iwork` which contains
                                    !!   the error flags on the cycles, in order to
                                    !!   look for eventual local integration
                                    !!   difficulties. if the position of a local
                                    !!   difficulty can be determined (e.g.
                                    !!   singularity, discontinuity within the
                                    !!   interval) one will probably gain from
                                    !!   splitting up the interval at this point
                                    !!   and calling appropriate integrators on
                                    !!   the subranges.
                                    !! * ier = 4 the extrapolation table constructed for
                                    !!   convergence acceleration of the series
                                    !!   formed by the integral contributions over
                                    !!   the cycles, does not converge to within
                                    !!   the requested accuracy. as in the case of
                                    !!   ier = 1, it is advised to examine the
                                    !!   array `iwork` which contains the error
                                    !!   flags on the cycles.
                                    !! * ier = 6 the input is invalid because
                                    !!   (`integr/=1` and `integr/=2`) or
                                    !!   `epsabs<=0` or `limlst<3`.
                                    !!   `result`, `abserr`, `neval`, `lst` are set
                                    !!   to zero.
                                    !! * ier = 7 bad integrand behaviour occurs within one
                                    !!   or more of the cycles. location and type
                                    !!   of the difficulty involved can be
                                    !!   determined from the vector `ierlst`. here
                                    !!   `lst` is the number of cycles actually
                                    !!   needed (see below):
                                    !!
                                    !!    * ierlst(k) = 1 the maximum number of
                                    !!      subdivisions (= `limit`) has
                                    !!      been achieved on the `k`th
                                    !!      cycle.
                                    !!    * ierlst(k) = 2 occurrence of roundoff error
                                    !!      is detected and prevents the
                                    !!      tolerance imposed on the
                                    !!      `k`th cycle, from being
                                    !!      achieved.
                                    !!    * ierlst(k) = 3 extremely bad integrand
                                    !!      behaviour occurs at some
                                    !!      points of the `k`th cycle.
                                    !!    * ierlst(k) = 4 the integration procedure
                                    !!      over the `k`th cycle does
                                    !!      not converge (to within the
                                    !!      required accuracy) due to
                                    !!      roundoff in the
                                    !!      extrapolation procedure
                                    !!      invoked on this cycle. it
                                    !!      is assumed that the result
                                    !!      on this interval is the
                                    !!      best which can be obtained.
                                    !!    * ierlst(k) = 5 the integral over the `k`th
                                    !!      cycle is probably divergent
                                    !!      or slowly convergent. it
                                    !!      must be noted that
                                    !!      divergence can occur with
                                    !!      any other value of
                                    !!      `ierlst(k)`.
                                    !!
                                    !! if `omega = 0` and `integr = 1`,
                                    !! the integral is calculated by means of [[dqagie]]
                                    !! and `ier = ierlst(1)` (with meaning as described
                                    !! for `ierlst(k), k = 1`).
        real(wp), intent(out) :: Rslst(Limlst) !! vector of dimension at least `limlst`
                                               !! `rslst(k)` contains the integral contribution
                                               !! over the interval `(a+(k-1)c,a+kc)` where
                                               !! `c = (2*int(abs(omega))+1)*pi/abs(omega)`,
                                               !! `k = 1, 2, ..., lst`.
                                               !! note that, if `omega = 0`, `rslst(1)` contains
                                               !! the value of the integral over `(a,infinity)`.
        real(wp), intent(out) :: Erlst(Limlst) !! vector of dimension at least `limlst`
                                               !! `erlst(k)` contains the error estimate corresponding
                                               !! with `rslst(k)`.
        integer, intent(out) :: Ierlst(Limlst) !! vector of dimension at least `limlst`
                                               !! `ierlst(k)` contains the error flag corresponding
                                               !! with `rslst(k)`. for the meaning of the local error
                                               !! flags see description of output parameter `ier`.
        integer, intent(out) :: Lst !! number of subintervals needed for the integration
                                    !! if `omega = 0` then lst is set to 1.
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, providing
                                            !! space for the quantities needed in the subdivision
                                            !! process of each cycle
        integer, intent(out) :: Nnlog(Limit) !! vector of dimension at least `limit`, providing
                                             !! space for the quantities needed in the subdivision
                                             !! process of each cycle
        real(wp), intent(out) :: Chebmo(Maxp1, 25) !! array of dimension at least `(maxp1,25)`, providing
                                                   !! space for the chebyshev moments needed within the
                                                   !! cycles (see also routine [[dqc25f]])

        real(wp) :: abseps, correc, dl, dla, drl, ep, eps, fact, p1, reseps, res3la(3)
        integer :: ktmin, l, last, ll, momcom, nev, nres, numrl2
        real(wp) :: psum(limexp + 2) !! `psum` contains the part of the epsilon table
                                     !! which is still needed for further computations.
                                     !! each element of `psum` is a partial sum of the
                                     !! series which should sum to the value of the
                                     !! integral.
        real(wp) :: c1, c2 !! end points of subinterval (of length cycle)
        real(wp) :: cycle !! `(2*int(abs(omega))+1)*pi/abs(omega)`
        real(wp) :: errsum  !! sum of error estimates over the subintervals,
                            !! calculated cumulatively
        real(wp) :: epsa !! absolute tolerance requested over current
                         !! subinterval

        real(wp), parameter :: p = 0.9_wp

        ! test on validity of parameters

        Result = 0.0_wp
        Abserr = 0.0_wp
        Neval = 0
        Lst = 0
        Ier = 0
        if ((Integr /= 1 .and. Integr /= 2) .or. Epsabs <= 0.0_wp .or. &
            Limlst < 3) Ier = 6
        if (Ier == 6) return

        if (Omega == 0.0_wp) then
            ! integration by dqagie if omega is zero
            if (Integr == 1) call dqagie(f, 0.0_wp, 1, Epsabs, 0.0_wp, &
                                            Limit, Result, Abserr, Neval, &
                                            Ier, Alist, Blist, Rlist, Elist, &
                                            Iord, last)
            Rslst(1) = Result
            Erlst(1) = Abserr
            Ierlst(1) = Ier
            Lst = 1
            return
        end if

        main : block

            ! initializations

            l = abs(Omega)
            dl = 2*l + 1
            cycle = dl*pi/abs(Omega)
            Ier = 0
            ktmin = 0
            Neval = 0
            numrl2 = 0
            nres = 0
            c1 = a
            c2 = cycle + a
            p1 = 1.0_wp - p
            eps = Epsabs
            if (Epsabs > uflow/p1) eps = Epsabs*p1
            ep = eps
            fact = 1.0_wp
            correc = 0.0_wp
            Abserr = 0.0_wp
            errsum = 0.0_wp

            ! main do-loop

            do Lst = 1, Limlst

                ! integrate over current subinterval.

                dla = Lst
                epsa = eps*fact
                call dqawoe(f, c1, c2, Omega, Integr, epsa, 0.0_wp, Limit, Lst, &
                            Maxp1, Rslst(Lst), Erlst(Lst), nev, Ierlst(Lst), &
                            last, Alist, Blist, Rlist, Elist, Iord, Nnlog, &
                            momcom, Chebmo)
                Neval = Neval + nev
                fact = fact*p
                errsum = errsum + Erlst(Lst)
                drl = 50.0_wp*abs(Rslst(Lst))

                ! test on accuracy with partial sum

                if ((errsum + drl) <= Epsabs .and. Lst >= 6) exit main
                correc = max(correc, Erlst(Lst))
                if (Ierlst(Lst) /= 0) eps = max(ep, correc*p1)
                if (Ierlst(Lst) /= 0) Ier = 7
                if (Ier == 7 .and. (errsum + drl) <= correc*10.0_wp .and. &
                    Lst > 5) exit main
                numrl2 = numrl2 + 1
                if (Lst > 1) then
                    psum(numrl2) = psum(ll) + Rslst(Lst)
                    if (Lst /= 2) then

                        ! test on maximum number of subintervals

                        if (Lst == Limlst) Ier = 1

                        ! perform new extrapolation

                        call dqelg(numrl2, psum, reseps, abseps, res3la, nres)

                        ! test whether extrapolated result is influenced by roundoff

                        ktmin = ktmin + 1
                        if (ktmin >= 15 .and. Abserr <= 0.1e-02_wp*(errsum + drl)) &
                            Ier = 4
                        if (abseps <= Abserr .or. Lst == 3) then
                            Abserr = abseps
                            Result = reseps
                            ktmin = 0

                            ! if ier is not 0, check whether direct result (partial sum)
                            ! or extrapolated result yields the best integral
                            ! approximation

                            if ((Abserr + 10.0_wp*correc) <= Epsabs .or. &
                                (Abserr <= Epsabs .and. &
                                    10.0_wp*correc >= Epsabs)) exit
                        end if
                        if (Ier /= 0 .and. Ier /= 7) exit
                    end if
                else
                    psum(1) = Rslst(1)
                end if
                ll = numrl2
                c1 = c2
                c2 = c2 + cycle
            end do

            ! set final result and error estimate

            Abserr = Abserr + 10.0_wp*correc
            if (Ier == 0) return
            if (Result == 0.0_wp .or. psum(numrl2) == 0.0_wp) then
                if (Abserr > errsum) exit main
                if (psum(numrl2) == 0.0_wp) return
            end if
            if (Abserr/abs(Result) <= (errsum + drl)/abs(psum(numrl2))) &
                then
                if (Ier >= 1 .and. Ier /= 7) Abserr = Abserr + drl
                return
            end if

        end block main

        Result = psum(numrl2)
        Abserr = errsum + drl

    end subroutine dqawfe
!********************************************************************************

!********************************************************************************
!>
!  1D integration of `cos(omega*x)*f(x)` or `sin(omega*x)*f(x)`
!  over a finite interval, adaptive subdivision with extrapolation
!
!  the routine calculates an approximation result to a given
!  definite integral i=integral of `f(x)*w(x)` over `(a,b)`
!  where `w(x) = cos(omega*x)` or `w(x) = sin(omega*x)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawo(f, a, b, Omega, Integr, Epsabs, Epsrel, Result, Abserr, &
                     Neval, Ier, Leniw, Maxp1, Lenw, Last, Iwork, Work)
        implicit none

        procedure(func) :: f !! function subprogram defining the function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Omega !! parameter in the integrand weight function
        integer, intent(in) :: Integr !! indicates which of the weight functions is used
                                      !!
                                      !! * integr = 1  `w(x) = cos(omega*x)`
                                      !! * integr = 2  `w(x) = sin(omega*x)`
                                      !!
                                      !! if `integr/=1.and.integr/=2`, the routine will
                                      !! end with ier = 6.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if `epsabs<=0` and
                                       !! `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine.
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   `(= leniw/2)` has been achieved. one can
                                    !!   allow more subdivisions by increasing the
                                    !!   value of leniw (and taking the according
                                    !!   dimension adjustments into account).
                                    !!   however, if this yields no improvement it
                                    !!   is advised to analyze the integrand in
                                    !!   order to determine the integration
                                    !!   difficulties. if the position of a local
                                    !!   difficulty can be determined (e.g.
                                    !!   singularity, discontinuity within the
                                    !!   interval) one will probably gain from
                                    !!   splitting up the interval at this point
                                    !!   and calling the integrator on the
                                    !!   subranges. if possible, an appropriate
                                    !!   special-purpose integrator should be used
                                    !!   which is designed for handling the type of
                                    !!   difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some interior points of the
                                    !!   integration interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table. it is presumed that
                                    !!   the requested tolerance cannot be achieved
                                    !!   due to roundoff in the extrapolation
                                    !!   table, and that the returned result is
                                    !!   the best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of `ier`.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!   or `(integr/=1 and integr/=2)`,
                                    !!   or `leniw<2` or `maxp1<1` or
                                    !!   `lenw<leniw*2+maxp1*25`.
                                    !!   `result`, `abserr`, `neval`, `last` are set to
                                    !!   zero. except when `leniw`, `maxp1` or `lenw` are
                                    !!   invalid, `work(limit*2+1)`, `work(limit*3+1)`,
                                    !!   `iwork(1)`, `iwork(limit+1)` are set to zero,
                                    !!   `work(1)` is set to `a` and `work(limit+1)` to
                                    !!   `b`.
        integer, intent(in) :: Leniw !! dimensioning parameter for `iwork`.
                                     !! `leniw/2` equals the maximum number of subintervals
                                     !! allowed in the partition of the given integration
                                     !! interval `(a,b)`, `leniw>=2`.
                                     !! if `leniw<2`, the routine will end with ier = 6.
        integer, intent(in) :: Maxp1 !! gives an upper bound on the number of chebyshev
                                     !! moments which can be stored, i.e. for the
                                     !! intervals of lengths `abs(b-a)*2**(-l)`,
                                     !! `l=0,1, ..., maxp1-2, maxp1>=1`
                                     !! if `maxp1<1`, the routine will end with ier = 6.
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`
                                    !! `lenw` must be at least `leniw*2+maxp1*25`.
                                    !! if `lenw<(leniw*2+maxp1*25)`, the routine will
                                    !! end with ier = 6.
        integer, intent(out) :: Last !! on return, `last` equals the number of subintervals
                                     !! produced in the subdivision process, which
                                     !! determines the number of significant elements
                                     !! actually in the work arrays.
        integer :: Iwork(Leniw) !! vector of dimension at least leniw
                                !! on return, the first `k` elements of which contain
                                !! pointers to the error estimates over the
                                !! subintervals, such that
                                !! `work(limit*3+iwork(1)), .., work(limit*3+iwork(k))`
                                !! form a decreasing
                                !! sequence, with `limit = lenw/2` , and `k = last`
                                !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                !! otherwise.
                                !! furthermore, `iwork(limit+1), ..., iwork(limit+last)`
                                !! indicate the subdivision levels of the
                                !! subintervals, such that `iwork(limit+i) = l` means
                                !! that the subinterval numbered `i` is of length
                                !! `abs(b-a)*2**(1-l)`.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`.
                               !! on return:
                               !!
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !! * `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end points,
                               !! * `work(limit*2+1), ..., work(limit*2+last)` contain
                               !!   the integral approximations over the
                               !!   subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3+last)`
                               !!   contain the error estimates.
                               !! * `work(limit*4+1), ..., work(limit*4+maxp1*25)`
                               !!   provide space for storing the chebyshev moments.
                               !!
                               !! note that `limit = lenw/2`.

        integer :: limit, lvl, l1, l2, l3, l4, momcom

        ! check validity of leniw, maxp1 and lenw.
        Ier = 6
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (Leniw >= 2 .and. Maxp1 >= 1 .and. Lenw >= (Leniw*2 + Maxp1*25)) then
            ! prepare call for dqawoe
            limit = Leniw/2
            l1 = limit + 1
            l2 = limit + l1
            l3 = limit + l2
            l4 = limit + l3
            call dqawoe(f, a, b, Omega, Integr, Epsabs, Epsrel, limit, 1, Maxp1, &
                        Result, Abserr, Neval, Ier, Last, Work(1), Work(l1), &
                        Work(l2), Work(l3), Iwork(1), Iwork(l1), momcom, &
                        Work(l4))
            ! call error handler if necessary
            lvl = 0
        end if
        if (Ier == 6) lvl = 1
        if (Ier /= 0) call xerror('abnormal return from dqawo', Ier, lvl)

    end subroutine dqawo
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqawo]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  definite integral
!  i = integral of `f(x)*w(x)` over `(a,b)`
!  where `w(x) = cos(omega*x)` or `w(x)=sin(omega*x)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawoe(f, a, b, Omega, Integr, Epsabs, Epsrel, Limit, Icall, &
                      Maxp1, Result, Abserr, Neval, Ier, Last, Alist, Blist, &
                      Rlist, Elist, Iord, Nnlog, Momcom, Chebmo)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Omega !! parameter in the integrand weight function
        integer, intent(in) :: Integr !! indicates which of the weight functions is to be
                                      !! used:
                                      !!
                                      !! * integr = 1  `w(x) = cos(omega*x)`
                                      !! * integr = 2  `w(x) = sin(omega*x)`
                                      !!
                                      !! if `integr/=1` and `integr/=2`, the routine
                                      !! will end with ier = 6.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested.
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        integer, intent(in) :: Limit !! gives an upper bound on the number of subdivisions
                                     !! in the partition of `(a,b)`, `limit>=1`.
        integer, intent(in) :: Icall !! if dqawoe is to be used only once, icall must
                                     !! be set to 1.  assume that during this call, the
                                     !! chebyshev moments (for clenshaw-curtis integration
                                     !! of degree 24) have been computed for intervals of
                                     !! lengths `(abs(b-a))*2**(-l), l=0,1,2,...momcom-1`.
                                     !! if `icall>1` this means that dqawoe has been
                                     !! called twice or more on intervals of the same
                                     !! length `abs(b-a)`. the chebyshev moments already
                                     !! computed are then re-used in subsequent calls.
                                     !! if `icall<1`, the routine will end with ier = 6.
        integer, intent(in) :: Maxp1 !! gives an upper bound on the number of chebyshev
                                     !! moments which can be stored, i.e. for the
                                     !! intervals of lengths `abs(b-a)*2**(-l)`,
                                     !! `l=0,1, ..., maxp1-2, maxp1>=1`.
                                     !! if `maxp1<1`, the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the
                                    !!   requested accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine.
                                    !!   the estimates for integral and error are
                                    !!   less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit (and taking according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand, in order to
                                    !!   determine the integration difficulties.
                                    !!   if the position of a local difficulty can
                                    !!   be determined (e.g. singularity,
                                    !!   discontinuity within the interval) one
                                    !!   will probably gain from splitting up the
                                    !!   interval at this point and calling the
                                    !!   integrator on the subranges. if possible,
                                    !!   an appropriate special-purpose integrator
                                    !!   should be used which is designed for
                                    !!   handling the type of difficulty involved.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !!   the error may be under-estimated.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 4 the algorithm does not converge.
                                    !!   roundoff error is detected in the
                                    !!   extrapolation table.
                                    !!   it is presumed that the requested
                                    !!   tolerance cannot be achieved due to
                                    !!   roundoff in the extrapolation table,
                                    !!   and that the returned result is the
                                    !!   best which can be obtained.
                                    !! * ier = 5 the integral is probably divergent, or
                                    !!   slowly convergent. it must be noted that
                                    !!   divergence can occur with any other value
                                    !!   of ier>0.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28_wp))`
                                    !!   or `(integr/=1 and integr/=2)` or
                                    !!   `icall<1` or `maxp1<1`.
                                    !!   `result`, `abserr`, `neval`, `last`, `rlist(1)`,
                                    !!   `elist(1)`, `iord(1)` and `nnlog(1)` are set
                                    !!   to zero. `alist(1)` and `blist(1)` are set
                                    !!   to `a` and `b` respectively.
        integer, intent(out) :: Last !! on return, `last` equals the number of
                                     !! subintervals produces in the subdivision
                                     !! process, which determines the number of
                                     !! significant elements actually in the
                                     !! work arrays.
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the integral
                                              !! approximations on the subintervals
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the moduli of the
                                              !! absolute error estimates on the subintervals
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! elements of which are pointers to the error
                                            !! estimates over the subintervals,
                                            !! such that `elist(iord(1)), ..., elist(iord(k))`
                                            !! form a decreasing sequence, with
                                            !! `k = last if last<=(limit/2+2)`, and
                                            !! `k = limit+1-last` otherwise.
        integer, intent(out) :: Nnlog(Limit) !! vector of dimension at least `limit`, containing the
                                             !! subdivision levels of the subintervals, i.e.
                                             !! iwork(i) = l means that the subinterval
                                             !! numbered `i` is of length `abs(b-a)*2**(1-l)`
        integer, intent(inout) :: Momcom !! indicating that the chebyshev moments
                                         !! have been computed for intervals of lengths
                                         !! `(abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1`,
                                         !! `momcom<maxp1`
        real(wp), intent(inout) :: Chebmo(Maxp1, 25) !! array of dimension `(maxp1,25)`
                                                     !! containing the chebyshev moments

        real(wp) :: rlist2(limexp + 2) !! array of dimension at least `limexp+2`
                                       !! containing the part of the epsilon table
                                       !! which is still needed for further computations
        integer :: maxerr !! pointer to the interval with largest error estimate
        real(wp) :: errmax !! `elist(maxerr)`
        real(wp) :: erlast !! error on the interval currently subdivided
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        real(wp) :: a1, area1, b1, error1 !! variable for the left subinterval
        real(wp) :: a2, area2, b2, error2 !! variable for the right subinterval
        integer :: nres !! number of calls to the extrapolation routine
        integer :: numrl2 !! number of elements in `rlist2`. if an appropriate
                          !! approximation to the compounded integral has
                          !! been obtained it is put in `rlist2(numrl2)` after
                          !! `numrl2` has been increased by one
        real(wp) :: small !! length of the smallest interval considered
                          !! up to now, multiplied by 1.5
        real(wp) :: erlarg !! sum of the errors over the intervals larger
                           !! than the smallest interval considered up to now
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        logical :: extrap !! logical variable denoting that the routine is
                          !! attempting to perform extrapolation, i.e. before
                          !! subdividing the smallest interval we try to
                          !! decrease the value of erlarg
        logical :: noext !! logical variable denoting that extrapolation
                         !! is no longer allowed (true  value)

        real(wp) :: abseps, correc, defab1, defab2, &
                    defabs, domega, dres, ertest, resabs, &
                    reseps, res3la(3), width
        integer :: id, ierro, iroff1, iroff2, iroff3, &
                   jupbnd, k, ksgn, ktmin, nev, nrmax, nrmom
        logical :: extall, done, test

        ! test on validity of parameters

        Ier = 0
        Neval = 0
        Last = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        Alist(1) = a
        Blist(1) = b
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        Nnlog(1) = 0
        if ((Integr /= 1 .and. Integr /= 2) .or. &
            (Epsabs <= 0.0_wp .and. Epsrel < max(50.0_wp*epmach, 0.5e-28_wp)) &
            .or. Icall < 1 .or. Maxp1 < 1) Ier = 6
        if (Ier == 6) return
        done = .false.

        ! first approximation to the integral

        domega = abs(Omega)
        nrmom = 0
        if (Icall <= 1) Momcom = 0
        call dqc25f(f, a, b, domega, Integr, nrmom, Maxp1, 0, Result, Abserr, &
                    Neval, defabs, resabs, Momcom, Chebmo)

        ! test on accuracy.

        dres = abs(Result)
        errbnd = max(Epsabs, Epsrel*dres)
        Rlist(1) = Result
        Elist(1) = Abserr
        Iord(1) = 1
        if (Abserr <= 100.0_wp*epmach*defabs .and. Abserr > errbnd) &
            Ier = 2
        if (Limit == 1) Ier = 1
        if (Ier /= 0 .or. Abserr <= errbnd) then
            if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result
            return
        end if

        ! initializations

        errmax = Abserr
        maxerr = 1
        area = Result
        errsum = Abserr
        Abserr = oflow
        nrmax = 1
        extrap = .false.
        noext = .false.
        ierro = 0
        iroff1 = 0
        iroff2 = 0
        iroff3 = 0
        ktmin = 0
        small = abs(b - a)*0.75_wp
        nres = 0
        numrl2 = 0
        extall = .false.
        if (0.5_wp*abs(b - a)*domega <= 2.0_wp) then
            numrl2 = 1
            extall = .true.
            rlist2(1) = Result
        end if
        if (0.25_wp*abs(b - a)*domega <= 2.0_wp) extall = .true.
        ksgn = -1
        if (dres >= (1.0_wp - 50.0_wp*epmach)*defabs) ksgn = 1

        ! main do-loop

        do Last = 2, Limit

            ! bisect the subinterval with the nrmax-th largest
            ! error estimate.

            nrmom = Nnlog(maxerr) + 1
            a1 = Alist(maxerr)
            b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
            a2 = b1
            b2 = Blist(maxerr)
            erlast = errmax
            call dqc25f(f, a1, b1, domega, Integr, nrmom, Maxp1, 0, area1, &
                        error1, nev, resabs, defab1, Momcom, Chebmo)
            Neval = Neval + nev
            call dqc25f(f, a2, b2, domega, Integr, nrmom, Maxp1, 1, area2, &
                        error2, nev, resabs, defab2, Momcom, Chebmo)
            Neval = Neval + nev

            ! improve previous approximations to integral
            ! and error and test for accuracy.

            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area = area + area12 - Rlist(maxerr)
            if (defab1 /= error1 .and. defab2 /= error2) then
                if (abs(Rlist(maxerr) - area12) <= 0.1e-4_wp*abs(area12) &
                    .and. erro12 >= 0.99_wp*errmax) then
                    if (extrap) iroff2 = iroff2 + 1
                    if (.not. extrap) iroff1 = iroff1 + 1
                end if
                if (Last > 10 .and. erro12 > errmax) iroff3 = iroff3 + 1
            end if
            Rlist(maxerr) = area1
            Rlist(Last) = area2
            Nnlog(maxerr) = nrmom
            Nnlog(Last) = nrmom
            errbnd = max(Epsabs, Epsrel*abs(area))

            ! test for roundoff error and eventually set error flag.

            if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) Ier = 2
            if (iroff2 >= 5) ierro = 3

            ! set error flag in the case that the number of
            ! subintervals equals limit.

            if (Last == Limit) Ier = 1

            ! set error flag in the case of bad integrand behaviour
            ! at a point of the integration range.

            if (max(abs(a1), abs(b2)) <= (1.0_wp + 100.0_wp*epmach) &
                *(abs(a2) + 1000.0_wp*uflow)) Ier = 4

            ! append the newly-created intervals to the list.

            if (error2 > error1) then
                Alist(maxerr) = a2
                Alist(Last) = a1
                Blist(Last) = b1
                Rlist(maxerr) = area2
                Rlist(Last) = area1
                Elist(maxerr) = error2
                Elist(Last) = error1
            else
                Alist(Last) = a2
                Blist(maxerr) = b1
                Blist(Last) = b2
                Elist(maxerr) = error1
                Elist(Last) = error2
            end if

            ! call subroutine dqpsrt to maintain the descending ordering
            ! in the list of error estimates and select the subinterval
            ! with nrmax-th largest error estimate (to bisected next).

            call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
            if (errsum <= errbnd) then
                ! ***jump out of do-loop
                done = .true.
                exit
            end if
            if (Ier /= 0) exit
            if (Last == 2 .and. extall) then
                small = small*0.5_wp
                numrl2 = numrl2 + 1
                rlist2(numrl2) = area
            else
                if (noext) cycle
                test = .true.
                if (extall) then
                    erlarg = erlarg - erlast
                    if (abs(b1 - a1) > small) erlarg = erlarg + erro12
                    if (extrap) test = .false.
                end if

                if (test) then
                    ! test whether the interval to be bisected next is the
                    ! smallest interval.

                    width = abs(Blist(maxerr) - Alist(maxerr))
                    if (width > small) cycle
                    if (extall) then
                        extrap = .true.
                        nrmax = 2
                    else
                        ! test whether we can start with the extrapolation procedure
                        ! (we do this if we integrate over the next interval with
                        ! use of a gauss-kronrod rule - see subroutine dqc25f).
                        small = small*0.5_wp
                        if (0.25_wp*width*domega > 2.0_wp) cycle
                        extall = .true.
                        ertest = errbnd
                        erlarg = errsum
                        cycle
                    end if
                end if

                if (ierro /= 3 .and. erlarg > ertest) then

                    ! the smallest interval has the largest error.
                    ! before bisecting decrease the sum of the errors over
                    ! the larger intervals (erlarg) and perform extrapolation.

                    jupbnd = Last
                    if (Last > (Limit/2 + 2)) jupbnd = Limit + 3 - Last
                    id = nrmax
                    do k = id, jupbnd
                        maxerr = Iord(nrmax)
                        errmax = Elist(maxerr)
                        if (abs(Blist(maxerr) - Alist(maxerr)) > small) &
                            cycle
                        nrmax = nrmax + 1
                    end do
                end if

                ! perform extrapolation.

                numrl2 = numrl2 + 1
                rlist2(numrl2) = area
                if (numrl2 >= 3) then
                    call dqelg(numrl2, rlist2, reseps, abseps, res3la, nres)
                    ktmin = ktmin + 1
                    if (ktmin > 5 .and. Abserr < 0.1e-02_wp*errsum) Ier = 5
                    if (abseps < Abserr) then
                        ktmin = 0
                        Abserr = abseps
                        Result = reseps
                        correc = erlarg
                        ertest = max(Epsabs, Epsrel*abs(reseps))
                        ! ***jump out of do-loop
                        if (Abserr <= ertest) exit
                    end if

                    ! prepare bisection of the smallest interval.

                    if (numrl2 == 1) noext = .true.
                    if (Ier == 5) exit
                end if
                maxerr = Iord(1)
                errmax = Elist(maxerr)
                nrmax = 1
                extrap = .false.
                small = small*0.5_wp
                erlarg = errsum
                cycle
            end if
            ertest = errbnd
            erlarg = errsum
        end do

        final : block

            if (done) exit final

            ! set the final result.

            if (Abserr /= oflow .and. nres /= 0) then
                if (Ier + ierro /= 0) then
                    if (ierro == 3) Abserr = Abserr + correc
                    if (Ier == 0) Ier = 3
                    if (Result == 0.0_wp .or. area == 0.0_wp) then
                        if (Abserr > errsum) exit final
                        if (area == 0.0_wp) then
                            if (Ier > 2) Ier = Ier - 1
                            if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result
                            return
                        end if
                    elseif (Abserr/abs(Result) > errsum/abs(area)) then
                        exit final
                    end if
                end if

                ! test on divergence.

                if (ksgn /= (-1) .or. max(abs(Result), abs(area)) &
                    > defabs*0.01_wp) then
                    if (0.01_wp > (Result/area) .or. (Result/area) &
                        > 100.0_wp .or. errsum >= abs(area)) Ier = 6
                end if
                if (Ier > 2) Ier = Ier - 1
                if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result
                return
            end if

        end block final

        ! compute global integral sum.

        Result = 0.0_wp
        do k = 1, Last
            Result = Result + Rlist(k)
        end do
        Abserr = errsum
        if (Ier > 2) Ier = Ier - 1
        if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result

    end subroutine dqawoe
!********************************************************************************

!********************************************************************************
!>
!  1D integration of functions with powers and or logs over a finite interval
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of `f*w` over `(a,b)`,
!  (where `w` shows a singular behaviour at the end points
!  see parameter `integr`).
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqaws(f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
                     abserr, neval, ier, limit, lenw, last, iwork, work)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration, b>a
                                  !! if b<=a, the routine will end with ier = 6.
        real(wp), intent(in) :: alfa !! parameter in the integrand function, `alfa>(-1)`
                                     !! if `alfa<=(-1)`, the routine will end with
                                     !! ier = 6.
        real(wp), intent(in) :: beta !! parameter in the integrand function, `beta>(-1)`
                                     !! if `beta<=(-1)`, the routine will end with
                                     !! ier = 6.
        integer, intent(in) :: integr !! indicates which weight function is to be used:
                                      !!
                                      !! * = 1  `(x-a)**alfa*(b-x)**beta`
                                      !! * = 2  `(x-a)**alfa*(b-x)**beta*log(x-a)`
                                      !! * = 3  `(x-a)**alfa*(b-x)**beta*log(b-x)`
                                      !! * = 4  `(x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)`
                                      !!
                                      !! if `integr<1` or `integr>4`, the routine
                                      !! will end with ier = 6.
        real(wp), intent(in) :: epsabs !! absolute accuracy requested
        real(wp), intent(in) :: epsrel !! relative accuracy requested.
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: result !! approximation to the integral
        real(wp), intent(out) :: abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: neval !! number of integrand evaluations
        integer, intent(out) :: ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for the integral and error
                                    !!   are less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit (and taking the according dimension
                                    !!   adjustments into account). however, if
                                    !!   this yields no improvement it is advised
                                    !!   to analyze the integrand, in order to
                                    !!   determine the integration difficulties
                                    !!   which prevent the requested tolerance from
                                    !!   being achieved. in case of a jump
                                    !!   discontinuity or a local singularity
                                    !!   of algebraico-logarithmic type at one or
                                    !!   more interior points of the integration
                                    !!   range, one should proceed by splitting up
                                    !!   the interval at these points and calling
                                    !!   the integrator on the subranges.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `b<=a` or `alfa<=(-1)` or `beta<=(-1)` or
                                    !!   or `integr<1` or `integr>4` or
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5e-28))`
                                    !!   or `limit<2` or `lenw<limit*4`.
                                    !!   `result`, `abserr`, `neval`, `last` are set to
                                    !!   zero. except when `lenw` or `limit` is invalid
                                    !!   `iwork(1)`, `work(limit*2+1)` and
                                    !!   `work(limit*3+1)` are set to zero, `work(1)`
                                    !!   is set to `a` and `work(limit+1)` to `b`.
        integer, intent(in) :: limit !! dimensioning parameter for `iwork`
                                    !! limit determines the maximum number of
                                    !! subintervals in the partition of the given
                                    !! integration interval `(a,b)`, `limit>=2`.
                                    !! if `limit<2`, the routine will end with ier = 6.
        integer, intent(in) :: lenw !! dimensioning parameter for `work`
                                    !! `lenw` must be at least `limit*4`.
                                    !! if `lenw<limit*4`, the routine will end
                                    !! with ier = 6.
        integer, intent(out) :: last !! on return, `last` equals the number of
                                     !! subintervals produced in the subdivision process,
                                     !! which determines the significant number of
                                     !! elements actually in the work arrays.
        integer :: iwork(limit) !! vector of dimension limit, the first `k`
                                !! elements of which contain pointers
                                !! to the error estimates over the subintervals,
                                !! such that `work(limit*3+iwork(1))`, ...,
                                !! `work(limit*3+iwork(k))` form a decreasing
                                !! sequence with `k = last` if `last<=(limit/2+2)`,
                                !! and `k = limit+1-last` otherwise
        real(wp) :: work(lenw) !! on return:
                               !!
                               !! * `work(1), ..., work(last)` contain the left
                               !!   end points of the subintervals in the
                               !!   partition of `(a,b)`,
                               !!   `work(limit+1), ..., work(limit+last)` contain
                               !!   the right end points,
                               !! * `work(limit*2+1), ..., work(limit*2+last)`
                               !!   contain the integral approximations over
                               !!   the subintervals,
                               !! * `work(limit*3+1), ..., work(limit*3+last)`
                               !!   contain the error estimates.

        integer :: lvl, l1, l2, l3

        ! check validity of limit and lenw.
        ier = 6
        neval = 0
        last = 0
        result = 0.0_wp
        abserr = 0.0_wp
        if (limit >= 2 .and. lenw >= limit*4) then

            ! prepare call for dqawse.

            l1 = limit + 1
            l2 = limit + l1
            l3 = limit + l2

            call dqawse(f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, &
                        abserr, neval, ier, work(1), work(l1), work(l2), work(l3), iwork, last)

            ! call error handler if necessary.
            lvl = 0
        end if
        if (ier == 6) lvl = 1
        if (ier /= 0) call xerror('abnormal return from dqaws', ier, lvl)

    end subroutine dqaws
!********************************************************************************

!********************************************************************************
!>
!  same as [[dqaws]] but provides more information and control
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of f*w over `(a,b)`,
!  (where `w` shows a singular behaviour at the end points,
!  see parameter integr).
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd)

    subroutine dqawse(f, a, b, Alfa, Beta, Integr, Epsabs, Epsrel, Limit, &
                      Result, Abserr, Neval, Ier, Alist, Blist, Rlist, Elist, &
                      Iord, Last)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration, `b>a`.
                                  !! if `b<=a`, the routine will end with ier = 6.
        real(wp), intent(in) :: Alfa !! parameter in the weight function, `alfa>(-1)`
                                     !! if `alfa<=(-1)`, the routine will end with
                                     !! ier = 6.
        real(wp), intent(in) :: Beta !! parameter in the weight function, `beta>(-1)`
                                     !! if `beta<=(-1)`, the routine will end with
                                     !! ier = 6.
        integer, intent(in) :: Integr !! indicates which weight function is to be used:
                                      !!
                                      !! * = 1  `(x-a)**alfa*(b-x)**beta`
                                      !! * = 2  `(x-a)**alfa*(b-x)**beta*log(x-a)`
                                      !! * = 3  `(x-a)**alfa*(b-x)**beta*log(b-x)`
                                      !! * = 4  `(x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)`
                                      !!
                                      !! if `integr<1` or `integr>4`, the routine
                                      !! will end with ier = 6.
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested.
                                       !! if `epsabs<=0`
                                       !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                       !! the routine will end with ier = 6.
        integer, intent(in) :: Limit !! gives an upper bound on the number of subintervals
                                     !! in the partition of `(a,b)`, `limit>=2`
                                     !! if `limit<2`, the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine
                                    !!   the estimates for the integral and error
                                    !!   are less reliable. it is assumed that the
                                    !!   requested accuracy has not been achieved.
                                    !!   error messages
                                    !! * ier = 1 maximum number of subdivisions allowed
                                    !!   has been achieved. one can allow more
                                    !!   subdivisions by increasing the value of
                                    !!   limit. however, if this yields no
                                    !!   improvement, it is advised to analyze the
                                    !!   integrand in order to determine the
                                    !!   integration difficulties which prevent the
                                    !!   requested tolerance from being achieved.
                                    !!   in case of a jump discontinuity or a local
                                    !!   singularity of algebraico-logarithmic type
                                    !!   at one or more interior points of the
                                    !!   integration range, one should proceed by
                                    !!   splitting up the interval at these
                                    !!   points and calling the integrator on the
                                    !!   subranges.
                                    !! * ier = 2 the occurrence of roundoff error is
                                    !!   detected, which prevents the requested
                                    !!   tolerance from being achieved.
                                    !! * ier = 3 extremely bad integrand behaviour occurs
                                    !!   at some points of the integration
                                    !!   interval.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `b<=a` or `alfa<=(-1)` or `beta<=(-1)`, or
                                    !!   `integr<1` or `integr>4`, or
                                    !!   `epsabs<=0` and
                                    !!   `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                    !!   or `limit<2`.
                                    !!   `result`, `abserr`, `neval`, `rlist(1)`, `elist(1)`,
                                    !!   `iord(1)` and `last` are set to zero. `alist(1)`
                                    !!   and `blist(1)` are set to `a` and `b`
                                    !!   respectively.
        real(wp), intent(out) :: Alist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the left
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Blist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the right
                                              !! end points of the subintervals in the partition
                                              !! of the given integration range `(a,b)`
        real(wp), intent(out) :: Rlist(Limit) !! vector of dimension at least `limit`,the first
                                              !! `last` elements of which are the integral
                                              !! approximations on the subintervals
        real(wp), intent(out) :: Elist(Limit) !! vector of dimension at least `limit`, the first
                                              !! `last` elements of which are the moduli of the
                                              !! absolute error estimates on the subintervals
        integer, intent(out) :: Iord(Limit) !! vector of dimension at least `limit`, the first `k`
                                            !! of which are pointers to the error
                                            !! estimates over the subintervals, so that
                                            !! `elist(iord(1)), ..., elist(iord(k))` with `k = last`
                                            !! if `last<=(limit/2+2)`, and `k = limit+1-last`
                                            !! otherwise form a decreasing sequence
        integer, intent(out) :: Last !! number of subintervals actually produced in
                                     !! the subdivision process

        real(wp) :: a1, b1, area1, error1 !! variable for the left subinterval
        real(wp) :: a2, b2, area2, error2 !! variable for the right subinterval
        real(wp) :: area12 !! `area1 + area2`
        real(wp) :: erro12 !! `error1 + error2`
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: errbnd !! requested accuracy `max(epsabs,epsrel*abs(result))`
        real(wp) :: errmax !! `elist(maxerr)`
        real(wp) :: errsum !! sum of the errors over the subintervals
        integer :: maxerr !! pointer to the interval with largest error estimate
        real(wp) :: centre, resas1, resas2, rg(25), rh(25), ri(25), rj(25)
        integer :: iroff1, iroff2, k, nev, nrmax

        ! test on validity of parameters

        Ier = 6
        Neval = 0
        Last = 0
        Rlist(1) = 0.0_wp
        Elist(1) = 0.0_wp
        Iord(1) = 0
        Result = 0.0_wp
        Abserr = 0.0_wp
        if (.not. (b <= a .or. (Epsabs == 0.0_wp .and. Epsrel < max(50.0_wp* &
            epmach, 0.5e-28_wp)) .or. Alfa <= (-1.0_wp) .or. Beta <= (-1.0_wp) &
            .or. Integr < 1 .or. Integr > 4 .or. Limit < 2)) then
            Ier = 0

            ! compute the modified chebyshev moments.

            call dqmomo(Alfa, Beta, ri, rj, rg, rh, Integr)

            ! integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).

            centre = 0.5_wp*(b + a)
            call dqc25s(f, a, b, a, centre, Alfa, Beta, ri, rj, rg, rh, area1, error1, &
                        resas1, Integr, nev)
            Neval = nev
            call dqc25s(f, a, b, centre, b, Alfa, Beta, ri, rj, rg, rh, area2, error2, &
                        resas2, Integr, nev)
            Last = 2
            Neval = Neval + nev
            Result = area1 + area2
            Abserr = error1 + error2

            ! test on accuracy.

            errbnd = max(Epsabs, Epsrel*abs(Result))

            ! initialization

            if (error2 > error1) then
                Alist(1) = centre
                Alist(2) = a
                Blist(1) = b
                Blist(2) = centre
                Rlist(1) = area2
                Rlist(2) = area1
                Elist(1) = error2
                Elist(2) = error1
            else
                Alist(1) = a
                Alist(2) = centre
                Blist(1) = centre
                Blist(2) = b
                Rlist(1) = area1
                Rlist(2) = area2
                Elist(1) = error1
                Elist(2) = error2
            end if
            Iord(1) = 1
            Iord(2) = 2
            if (Limit == 2) Ier = 1
            if (Abserr > errbnd .and. Ier /= 1) then
                errmax = Elist(1)
                maxerr = 1
                nrmax = 1
                area = Result
                errsum = Abserr
                iroff1 = 0
                iroff2 = 0

                ! main do-loop

                do Last = 3, Limit

                    ! bisect the subinterval with largest error estimate.

                    a1 = Alist(maxerr)
                    b1 = 0.5_wp*(Alist(maxerr) + Blist(maxerr))
                    a2 = b1
                    b2 = Blist(maxerr)

                    call dqc25s(f, a, b, a1, b1, Alfa, Beta, ri, rj, rg, rh, area1, &
                                error1, resas1, Integr, nev)
                    Neval = Neval + nev
                    call dqc25s(f, a, b, a2, b2, Alfa, Beta, ri, rj, rg, rh, area2, &
                                error2, resas2, Integr, nev)
                    Neval = Neval + nev

                    ! improve previous approximations integral and error
                    ! and test for accuracy.

                    area12 = area1 + area2
                    erro12 = error1 + error2
                    errsum = errsum + erro12 - errmax
                    area = area + area12 - Rlist(maxerr)
                    if (a /= a1 .and. b /= b2) then
                        if (resas1 /= error1 .and. resas2 /= error2) then
                            ! test for roundoff error.
                            if (abs(Rlist(maxerr) - area12) &
                                < 0.1e-4_wp*abs(area12) .and. &
                                erro12 >= 0.99_wp*errmax) iroff1 = iroff1 + 1
                            if (Last > 10 .and. erro12 > errmax) &
                                iroff2 = iroff2 + 1
                        end if
                    end if
                    Rlist(maxerr) = area1
                    Rlist(Last) = area2

                    ! test on accuracy.

                    errbnd = max(Epsabs, Epsrel*abs(area))
                    if (errsum > errbnd) then

                        ! set error flag in the case that the number of interval
                        ! bisections exceeds limit.

                        if (Last == Limit) Ier = 1

                        ! set error flag in the case of roundoff error.

                        if (iroff1 >= 6 .or. iroff2 >= 20) Ier = 2

                        ! set error flag in the case of bad integrand behaviour
                        ! at interior points of integration range.

                        if (max(abs(a1), abs(b2)) &
                            <= (1.0_wp + 100.0_wp*epmach) &
                            *(abs(a2) + 1000.0_wp*uflow)) Ier = 3
                    end if

                    ! append the newly-created intervals to the list.

                    if (error2 > error1) then
                        Alist(maxerr) = a2
                        Alist(Last) = a1
                        Blist(Last) = b1
                        Rlist(maxerr) = area2
                        Rlist(Last) = area1
                        Elist(maxerr) = error2
                        Elist(Last) = error1
                    else
                        Alist(Last) = a2
                        Blist(maxerr) = b1
                        Blist(Last) = b2
                        Elist(maxerr) = error1
                        Elist(Last) = error2
                    end if

                    ! call subroutine dqpsrt to maintain the descending ordering
                    ! in the list of error estimates and select the subinterval
                    ! with largest error estimate (to be bisected next).

                    call dqpsrt(Limit, Last, maxerr, errmax, Elist, Iord, nrmax)
                    ! ***jump out of do-loop
                    if (Ier /= 0 .or. errsum <= errbnd) exit
                end do

                ! compute final result.
                Result = 0.0_wp
                do k = 1, Last
                    Result = Result + Rlist(k)
                end do
                Abserr = errsum
            end if
        end if

    end subroutine dqawse
!********************************************************************************

!********************************************************************************
!>
!  1D integral for Cauchy principal values using a 25 point quadrature rule
!
!  to compute i = integral of `f*w` over `(a,b)` with
!  error estimate, where `w(x) = 1/(x-c)`
!
!### History
!  * QUADPACK: date written 810101, revision date 830518 (yymmdd)

    subroutine dqc25c(f, a, b, c, Result, Abserr, Krul, Neval)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! left end point of the integration interval
        real(wp), intent(in) :: b !! right end point of the integration interval, `b>a`
        real(wp), intent(in) :: c !! parameter in the weight function
        real(wp), intent(out) :: Result !! approximation to the integral.
                                        !! `result` is computed by using a generalized
                                        !! clenshaw-curtis method if `c` lies within ten percent
                                        !! of the integration interval. in the other case the
                                        !! 15-point kronrod rule obtained by optimal addition
                                        !! of abscissae to the 7-point gauss rule, is applied.
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(inout) :: Krul !! key which is decreased by 1 if the 15-point
                                       !! gauss-kronrod scheme has been used
        integer, intent(out) :: Neval !! number of integrand evaluations

        real(wp) :: ak22, amom0, amom1, amom2, cc, &
                    p2, p3, p4, resabs, resasc, u
        integer :: i, isym, k, kp
        real(wp) :: fval(25) !! value of the function `f` at the points
                             !! `cos(k*pi/24)`, `k = 0, ..., 24`
        real(wp) :: cheb12(13) !! chebyshev series expansion coefficients,
                               !! for the function `f`, of degree 12
        real(wp) :: cheb24(25) !! chebyshev series expansion coefficients,
                               !! for the function `f`, of degree 24
        real(wp) :: res12 !! approximation to the integral corresponding
                          !! to the use of cheb12
        real(wp) :: res24 !! approximation to the integral corresponding
                          !! to the use of cheb24
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: centr !! mid point of the interval

        real(wp), dimension(11), parameter :: x = [(cos(k*pi/24.0_wp), k=1, 11)]
            !! the vector x contains the values `cos(k*pi/24)`,
            !! `k = 1, ..., 11`, to be used for the chebyshev series
            !! expansion of `f`

        ! check the position of c.

        cc = (2.0_wp*c - b - a)/(b - a)
        if (abs(cc) < 1.1_wp) then

            ! use the generalized clenshaw-curtis method.

            hlgth = 0.5_wp*(b - a)
            centr = 0.5_wp*(b + a)
            Neval = 25
            fval(1) = 0.5_wp*f(hlgth + centr)
            fval(13) = f(centr)
            fval(25) = 0.5_wp*f(centr - hlgth)
            do i = 2, 12
                u = hlgth*x(i - 1)
                isym = 26 - i
                fval(i) = f(u + centr)
                fval(isym) = f(centr - u)
            end do

            ! compute the chebyshev series expansion.

            call dqcheb(x, fval, cheb12, cheb24)

            ! the modified chebyshev moments are computed by forward
            ! recursion, using amom0 and amom1 as starting values.

            amom0 = log(abs((1.0_wp - cc)/(1.0_wp + cc)))
            amom1 = 2.0_wp + cc*amom0
            res12 = cheb12(1)*amom0 + cheb12(2)*amom1
            res24 = cheb24(1)*amom0 + cheb24(2)*amom1
            do k = 3, 13
                amom2 = 2.0_wp*cc*amom1 - amom0
                ak22 = (k - 2)*(k - 2)
                if ((k/2)*2 == k) amom2 = amom2 - 4.0_wp/(ak22 - 1.0_wp)
                res12 = res12 + cheb12(k)*amom2
                res24 = res24 + cheb24(k)*amom2
                amom0 = amom1
                amom1 = amom2
            end do
            do k = 14, 25
                amom2 = 2.0_wp*cc*amom1 - amom0
                ak22 = (k - 2)*(k - 2)
                if ((k/2)*2 == k) amom2 = amom2 - 4.0_wp/(ak22 - 1.0_wp)
                res24 = res24 + cheb24(k)*amom2
                amom0 = amom1
                amom1 = amom2
            end do
            Result = res24
            Abserr = abs(res24 - res12)
        else

            ! apply the 15-point gauss-kronrod scheme.

            ! dqwgtc - external function subprogram defining the weight function

            Krul = Krul - 1
            call dqk15w(f, dqwgtc, c, p2, p3, p4, kp, a, b, Result, Abserr, resabs, &
                        resasc)
            Neval = 15
            if (resasc == Abserr) Krul = Krul + 1
        end if

    end subroutine dqc25c
!********************************************************************************

!********************************************************************************
!>
!  1D integral for sin/cos integrand using a 25 point quadrature rule
!
!  to compute the integral i=integral of `f(x)` over `(a,b)`
!  where `w(x) = cos(omega*x)` or `w(x)=sin(omega*x)` and to
!  compute j = integral of `abs(f)` over `(a,b)`. for small value
!  of `omega` or small intervals `(a,b)` the 15-point gauss-kronrod
!  rule is used. otherwise a generalized clenshaw-curtis
!  method is used.
!
!### History
!  * QUADPACK: date written 810101, revision date 211011 (yymmdd)

    subroutine dqc25f(f, a, b, Omega, Integr, Nrmom, Maxp1, Ksave, Result, &
                      Abserr, Neval, Resabs, Resasc, Momcom, Chebmo)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand
                             !! function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Omega !! parameter in the weight function
        integer, intent(in) :: Integr !! indicates which weight function is to be used
                                      !!
                                      !! * integr = 1   `w(x) = cos(omega*x)`
                                      !! * integr = 2   `w(x) = sin(omega*x)`
        integer, intent(in) :: Nrmom !! the length of interval `(a,b)` is equal to the length
                                     !! of the original integration interval divided by
                                     !! `2**nrmom` (we suppose that the routine is used in an
                                     !! adaptive integration process, otherwise set
                                     !! nrmom = 0). `nrmom` must be zero at the first call.
        integer, intent(in) :: Maxp1 !! gives an upper bound on the number of chebyshev
                                     !! moments which can be stored, i.e. for the
                                     !! intervals of lengths `abs(bb-aa)*2**(-l)`,
                                     !! `l = 0,1,2, ..., maxp1-2`.
        integer, intent(in) :: Ksave !! key which is one when the moments for the
                                     !! current interval have been computed
        real(wp), intent(out) :: Result !! approximation to the integral i
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute
                                        !! error, which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`
        integer, intent(inout) :: Momcom !! for each interval length we need to compute the
                                         !! chebyshev moments. momcom counts the number of
                                         !! intervals for which these moments have already been
                                         !! computed. if `nrmom<momcom` or `ksave = 1`, the
                                         !! chebyshev moments for the interval `(a,b)` have
                                         !! already been computed and stored, otherwise we
                                         !! compute them and we increase momcom.
        real(wp), intent(inout) :: Chebmo(Maxp1, 25) !! array of dimension at least `(maxp1,25)` containing
                                                     !! the modified chebyshev moments for the first `momcom`
                                                     !! `momcom` interval lengths

        real(wp) :: ac, an, an2, as, asap, ass, conc, &
                    cons, cospar, d(25), d1(25), d2(25), &
                    estc, ests, parint, par2, par22, &
                    p2, p3, p4, sinpar, v(28)
        integer :: i, iers, isym, j, k, m, noequ, noeq1
        real(wp) :: centr !! mid point of the integration interval
        real(wp) :: hlgth !! half-length of the integration interval
        real(wp) :: fval(25) !! value of the function `f` at the points
                             !! `(b-a)*0.5*cos(k*pi/12) + (b+a)*0.5`, `k = 0, ..., 24`
        real(wp) :: cheb12(13) !! coefficients of the chebyshev series expansion
                               !! of degree 12, for the function `f`, in the
                               !! interval `(a,b)`
        real(wp) :: cheb24(25) !! coefficients of the chebyshev series expansion
                               !! of degree 24, for the function `f`, in the
                               !! interval `(a,b)`
        real(wp) :: resc12 !! approximation to the integral of
                           !! `cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))`
                           !! over `(-1,+1)`, using the chebyshev series
                           !! expansion of degree 12
        real(wp) :: resc24 !! approximation to the same integral, using the
                           !! chebyshev series expansion of degree 24
        real(wp) :: ress12 !! the analogue of `resc12` for the sine
        real(wp) :: ress24 !! the analogue of `resc24` for the sine

        real(wp), dimension(11), parameter :: x = [(cos(k*pi/24.0_wp), k=1, 11)]
            !! the vector x contains the values `cos(k*pi/24)`,
            !! `k = 1, ..., 11`, to be used for the chebyshev series
            !! expansion of `f`

        centr = 0.5_wp*(b + a)
        hlgth = 0.5_wp*(b - a)
        parint = Omega*hlgth

        ! compute the integral using the 15-point gauss-kronrod
        ! formula if the value of the parameter in the integrand
        ! is small.

        if (abs(parint) > 2.0_wp) then

            ! compute the integral using the generalized clenshaw-
            ! curtis method.

            conc = hlgth*cos(centr*Omega)
            cons = hlgth*sin(centr*Omega)
            Resasc = oflow
            Neval = 25

            ! check whether the chebyshev moments for this interval
            ! have already been computed.

            if (Nrmom >= Momcom .and. Ksave /= 1) then

                ! compute a new set of chebyshev moments.

                m = Momcom + 1
                par2 = parint*parint
                par22 = par2 + 2.0_wp
                sinpar = sin(parint)
                cospar = cos(parint)

                ! compute the chebyshev moments with respect to cosine.

                v(1) = 2.0_wp*sinpar/parint
                v(2) = (8.0_wp*cospar + (par2 + par2 - 8.0_wp)*sinpar/parint) &
                       /par2
                v(3) = (32.0_wp*(par2 - 12.0_wp)*cospar + (2.0_wp*((par2 - &
                       80.0_wp)*par2 + 192.0_wp)*sinpar)/parint)/(par2*par2)
                ac = 8.0_wp*cospar
                as = 24.0_wp*parint*sinpar
                if (abs(parint) > 24.0_wp) then

                    ! compute the chebyshev moments by means of forward
                    ! recursion.

                    an = 4.0_wp
                    do i = 4, 13
                        an2 = an*an
                        v(i) = ((an2 - 4.0_wp)*(2.0_wp*(par22 - an2 - an2)*v(i - 1) - &
                               ac) + as - par2*(an + 1.0_wp)*(an + 2.0_wp)*v(i - 2)) &
                               /(par2*(an - 1.0_wp)*(an - 2.0_wp))
                        an = an + 2.0_wp
                    end do
                else

                    ! compute the chebyshev moments as the solutions of a
                    ! boundary value problem with 1 initial value (v(3)) and 1
                    ! end value (computed using an asymptotic formula).

                    noequ = 25
                    noeq1 = noequ - 1
                    an = 6.0_wp
                    do k = 1, noeq1
                        an2 = an*an
                        d(k) = -2.0_wp*(an2 - 4.0_wp)*(par22 - an2 - an2)
                        d2(k) = (an - 1.0_wp)*(an - 2.0_wp)*par2
                        d1(k + 1) = (an + 3.0_wp)*(an + 4.0_wp)*par2
                        v(k + 3) = as - (an2 - 4.0_wp)*ac
                        an = an + 2.0_wp
                    end do
                    an2 = an*an
                    d(noequ) = -2.0_wp*(an2 - 4.0_wp)*(par22 - an2 - an2)
                    v(noequ + 3) = as - (an2 - 4.0_wp)*ac
                    v(4) = v(4) - 56.0_wp*par2*v(3)
                    ass = parint*sinpar
                    asap = (((((210.0_wp*par2 - 1.0_wp)*cospar - (105.0_wp* &
                            par2 - 63.0_wp)*ass)/an2 - (1.0_wp - 15.0_wp*par2) &
                            *cospar + 15.0_wp*ass)/an2 - cospar + 3.0_wp*ass) &
                            /an2 - cospar)/an2
                    v(noequ + 3) = v(noequ + 3) - 2.0_wp*asap*par2*(an - 1.0_wp) &
                                   *(an - 2.0_wp)

                    ! solve the tridiagonal system by means of gaussian
                    ! elimination with partial pivoting.

                    call dgtsl(noequ, d1, d, d2, v(4), iers)
                end if
                do j = 1, 13
                    Chebmo(m, 2*j - 1) = v(j)
                end do

                ! compute the chebyshev moments with respect to sine.

                v(1) = 2.0_wp*(sinpar - parint*cospar)/par2
                v(2) = (18.0_wp - 48.0_wp/par2)*sinpar/par2 + &
                       (-2.0_wp + 48.0_wp/par2)*cospar/parint
                ac = -24.0_wp*parint*cospar
                as = -8.0_wp*sinpar
                if (abs(parint) > 24.0_wp) then

                    ! compute the chebyshev moments by means of forward recursion.

                    an = 3.0_wp
                    do i = 3, 12
                        an2 = an*an
                        v(i) = ((an2 - 4.0_wp)*(2.0_wp*(par22 - an2 - an2)*v(i - 1) + &
                                as) + ac - par2*(an + 1.0_wp)*(an + 2.0_wp)*v(i - 2)) &
                                /(par2*(an - 1.0_wp)*(an - 2.0_wp))
                        an = an + 2.0_wp
                    end do
                else

                    ! compute the chebyshev moments as the solutions of a boundary
                    ! value problem with 1 initial value (v(2)) and 1 end value
                    ! (computed using an asymptotic formula).

                    an = 5.0_wp
                    do k = 1, noeq1
                        an2 = an*an
                        d(k) = -2.0_wp*(an2 - 4.0_wp)*(par22 - an2 - an2)
                        d2(k) = (an - 1.0_wp)*(an - 2.0_wp)*par2
                        d1(k + 1) = (an + 3.0_wp)*(an + 4.0_wp)*par2
                        v(k + 2) = ac + (an2 - 4.0_wp)*as
                        an = an + 2.0_wp
                    end do
                    an2 = an*an
                    d(noequ) = -2.0_wp*(an2 - 4.0_wp)*(par22 - an2 - an2)
                    v(noequ + 2) = ac + (an2 - 4.0_wp)*as
                    v(3) = v(3) - 42.0_wp*par2*v(2)
                    ass = parint*cospar
                    asap = (((((105.0_wp*par2 - 63.0_wp)*ass + (210.0_wp*par2 - &
                            1.0_wp)*sinpar)/an2 + (15.0_wp*par2 - 1.0_wp) &
                            *sinpar - 15.0_wp*ass)/an2 - 3.0_wp*ass - sinpar) &
                            /an2 - sinpar)/an2
                    v(noequ + 2) = v(noequ + 2) - 2.0_wp*asap*par2*(an - 1.0_wp) &
                                   *(an - 2.0_wp)

                    ! solve the tridiagonal system by means of gaussian
                    ! elimination with partial pivoting.

                    call dgtsl(noequ, d1, d, d2, v(3), iers)
                end if
                do j = 1, 12
                    Chebmo(m, 2*j) = v(j)
                end do
            end if
            if (Nrmom < Momcom) m = Nrmom + 1
            if (Momcom < (Maxp1 - 1) .and. Nrmom >= Momcom) Momcom = Momcom + 1

            ! compute the coefficients of the chebyshev expansions
            ! of degrees 12 and 24 of the function f.

            fval(1) = 0.5_wp*f(centr + hlgth)
            fval(13) = f(centr)
            fval(25) = 0.5_wp*f(centr - hlgth)
            do i = 2, 12
                isym = 26 - i
                fval(i) = f(hlgth*x(i - 1) + centr)
                fval(isym) = f(centr - hlgth*x(i - 1))
            end do
            call dqcheb(x, fval, cheb12, cheb24)

            ! compute the integral and error estimates.

            resc12 = cheb12(13)*Chebmo(m, 13)
            ress12 = 0.0_wp
            k = 11
            do j = 1, 6
                resc12 = resc12 + cheb12(k)*Chebmo(m, k)
                ress12 = ress12 + cheb12(k + 1)*Chebmo(m, k + 1)
                k = k - 2
            end do
            resc24 = cheb24(25)*Chebmo(m, 25)
            ress24 = 0.0_wp
            Resabs = abs(cheb24(25))
            k = 23
            do j = 1, 12
                resc24 = resc24 + cheb24(k)*Chebmo(m, k)
                ress24 = ress24 + cheb24(k + 1)*Chebmo(m, k + 1)
                Resabs = Resabs + abs(cheb24(k)) + abs(cheb24(k + 1))
                k = k - 2
            end do
            estc = abs(resc24 - resc12)
            ests = abs(ress24 - ress12)
            Resabs = Resabs*abs(hlgth)
            if (Integr == 2) then
                Result = conc*ress24 + cons*resc24
                Abserr = abs(conc*ests) + abs(cons*estc)
            else
                Result = conc*resc24 - cons*ress24
                Abserr = abs(conc*estc) + abs(cons*ests)
            end if
        else
            call dqk15w(f, dqwgtf, Omega, p2, p3, p4, Integr, a, b, Result, Abserr, &
                        Resabs, Resasc)
            Neval = 15
        end if

    end subroutine dqc25f
!********************************************************************************

!********************************************************************************
!>
!  25-point clenshaw-curtis integration
!
!  to compute i = integral of `f*w` over `(bl,br)`, with error
!  estimate, where the weight function `w` has a singular
!  behaviour of algebraico-logarithmic type at the points
!  `a` and/or `b`. `(bl,br)` is a part of `(a,b)`.
!
!### History
!  * QUADPACK: date written 810101, revision date 830518 (yymmdd)

    subroutine dqc25s(f, a, b, Bl, Br, Alfa, Beta, Ri, Rj, Rg, Rh, Result, Abserr, &
                      Resasc, Integr, Nev)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand f(x).
        real(wp), intent(in) :: a !! left end point of the original interval
        real(wp), intent(in) :: b !! right end point of the original interval, `b>a`
        real(wp), intent(in) :: Bl !! lower limit of integration, `bl>=a`
        real(wp), intent(in) :: Br !! upper limit of integration, `br<=b`
        real(wp), intent(in) :: Alfa !! parameter in the weight function
        real(wp), intent(in) :: Beta !! parameter in the weight function
        real(wp), intent(in) :: Ri(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine [[dqmomo]])
        real(wp), intent(in) :: Rj(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine [[dqmomo]])
        real(wp), intent(in) :: Rg(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine [[dqmomo]])
        real(wp), intent(in) :: Rh(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine [[dqmomo]])
        real(wp), intent(out) :: Result !! approximation to the integral
                                        !! `result` is computed by using a generalized
                                        !! clenshaw-curtis method if `b1 = a` or `br = b`.
                                        !! in all other cases the 15-point kronrod
                                        !! rule is applied, obtained by optimal addition of
                                        !! abscissae to the 7-point gauss rule.
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Resasc !! approximation to the integral of abs(f*w-i/(b-a))
        integer, intent(in) :: Integr !! which determines the weight function
                                      !! * = 1  `w(x) = (x-a)**alfa*(b-x)**beta`
                                      !! * = 2  `w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)`
                                      !! * = 3  `w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)`
                                      !! * = 4  `w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)`
        integer, intent(out) :: Nev !! number of integrand evaluations

        real(wp) :: cheb12(13) !! coefficients of the chebyshev series expansion
                               !! of degree 12, for the function `f`, in the
                               !! interval `(bl,br)`
        real(wp) :: cheb24(25) !! coefficients of the chebyshev series expansion
                               !! of degree 24, for the function `f`, in the
                               !! interval `(bl,br)`
        real(wp) :: fval(25) !! value of the function f at the points
                             !! `(br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5`
                             !! `k = 0, ..., 24`
        real(wp) :: res12 !! approximation to the integral obtained from `cheb12`
        real(wp) :: res24 !! approximation to the integral obtained from `cheb24`
        real(wp) :: hlgth !! half-length of the interval `(bl,br)`
        real(wp) :: centr !! mid point of the interval `(bl,br)`
        integer :: k !! counter for `x`
        real(wp) :: dc, factor, fix, resabs, u
        integer :: i, isym

        real(wp), dimension(11), parameter :: x = [(cos(k*pi/24.0_wp), k=1, 11)]
            !! the vector x contains the values `cos(k*pi/24)`,
            !! `k = 1, ..., 11`, to be used for the chebyshev series
            !! expansion of `f`

        Nev = 25
        if (Bl == a .and. (Alfa /= 0.0_wp .or. Integr == 2 .or. Integr == 4)) &
            then

            ! this part of the program is executed only if a = bl.

            ! compute the chebyshev series expansion of the
            ! following function
            ! f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
            !      *f(0.5*(br-a)*x+0.5*(br+a))

            hlgth = 0.5_wp*(Br - Bl)
            centr = 0.5_wp*(Br + Bl)
            fix = b - centr
            fval(1) = 0.5_wp*f(hlgth + centr)*(fix - hlgth)**Beta
            fval(13) = f(centr)*(fix**Beta)
            fval(25) = 0.5_wp*f(centr - hlgth)*(fix + hlgth)**Beta
            do i = 2, 12
                u = hlgth*x(i - 1)
                isym = 26 - i
                fval(i) = f(u + centr)*(fix - u)**Beta
                fval(isym) = f(centr - u)*(fix + u)**Beta
            end do
            factor = hlgth**(Alfa + 1.0_wp)
            Result = 0.0_wp
            Abserr = 0.0_wp
            res12 = 0.0_wp
            res24 = 0.0_wp
            if (Integr > 2) then

                ! compute the chebyshev series expansion of the
                ! following function
                ! f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)

                fval(1) = fval(1)*log(fix - hlgth)
                fval(13) = fval(13)*log(fix)
                fval(25) = fval(25)*log(fix + hlgth)
                do i = 2, 12
                    u = hlgth*x(i - 1)
                    isym = 26 - i
                    fval(i) = fval(i)*log(fix - u)
                    fval(isym) = fval(isym)*log(fix + u)
                end do
                call dqcheb(x, fval, cheb12, cheb24)

                ! integr = 3  (or 4)

                do i = 1, 13
                    res12 = res12 + cheb12(i)*Ri(i)
                    res24 = res24 + cheb24(i)*Ri(i)
                end do
                do i = 14, 25
                    res24 = res24 + cheb24(i)*Ri(i)
                end do
                if (Integr /= 3) then

                    ! integr = 4

                    dc = log(Br - Bl)
                    Result = res24*dc
                    Abserr = abs((res24 - res12)*dc)
                    res12 = 0.0_wp
                    res24 = 0.0_wp
                    do i = 1, 13
                        res12 = res12 + cheb12(i)*Rg(i)
                        res24 = res24 + cheb24(i)*Rg(i)
                    end do
                    do i = 14, 25
                        res24 = res24 + cheb24(i)*Rg(i)
                    end do
                end if
            else
                call dqcheb(x, fval, cheb12, cheb24)

                ! integr = 1  (or 2)

                do i = 1, 13
                    res12 = res12 + cheb12(i)*Ri(i)
                    res24 = res24 + cheb24(i)*Ri(i)
                end do
                do i = 14, 25
                    res24 = res24 + cheb24(i)*Ri(i)
                end do
                if (Integr /= 1) then

                    ! integr = 2

                    dc = log(Br - Bl)
                    Result = res24*dc
                    Abserr = abs((res24 - res12)*dc)
                    res12 = 0.0_wp
                    res24 = 0.0_wp
                    do i = 1, 13
                        res12 = res12 + cheb12(i)*Rg(i)
                        res24 = res12 + cheb24(i)*Rg(i)
                    end do
                    do i = 14, 25
                        res24 = res24 + cheb24(i)*Rg(i)
                    end do
                end if
            end if
            Result = (Result + res24)*factor
            Abserr = (Abserr + abs(res24 - res12))*factor
        elseif (Br == b .and. (Beta /= 0.0_wp .or. Integr == 3 .or. Integr == 4)) then

            ! this part of the program is executed only if b = br.

            ! compute the chebyshev series expansion of the
            ! following function
            ! f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
            !      *f(0.5*(b-bl)*x+0.5*(b+bl))

            hlgth = 0.5_wp*(Br - Bl)
            centr = 0.5_wp*(Br + Bl)
            fix = centr - a
            fval(1) = 0.5_wp*f(hlgth + centr)*(fix + hlgth)**Alfa
            fval(13) = f(centr)*(fix**Alfa)
            fval(25) = 0.5_wp*f(centr - hlgth)*(fix - hlgth)**Alfa
            do i = 2, 12
                u = hlgth*x(i - 1)
                isym = 26 - i
                fval(i) = f(u + centr)*(fix + u)**Alfa
                fval(isym) = f(centr - u)*(fix - u)**Alfa
            end do
            factor = hlgth**(Beta + 1.0_wp)
            Result = 0.0_wp
            Abserr = 0.0_wp
            res12 = 0.0_wp
            res24 = 0.0_wp
            if (Integr == 2 .or. Integr == 4) then

                ! compute the chebyshev series expansion of the
                ! following function
                ! f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))

                fval(1) = fval(1)*log(hlgth + fix)
                fval(13) = fval(13)*log(fix)
                fval(25) = fval(25)*log(fix - hlgth)
                do i = 2, 12
                    u = hlgth*x(i - 1)
                    isym = 26 - i
                    fval(i) = fval(i)*log(u + fix)
                    fval(isym) = fval(isym)*log(fix - u)
                end do
                call dqcheb(x, fval, cheb12, cheb24)

                ! integr = 2  (or 4)

                do i = 1, 13
                    res12 = res12 + cheb12(i)*Rj(i)
                    res24 = res24 + cheb24(i)*Rj(i)
                end do
                do i = 14, 25
                    res24 = res24 + cheb24(i)*Rj(i)
                end do
                if (Integr /= 2) then
                    dc = log(Br - Bl)
                    Result = res24*dc
                    Abserr = abs((res24 - res12)*dc)
                    res12 = 0.0_wp
                    res24 = 0.0_wp

                    ! integr = 4

                    do i = 1, 13
                        res12 = res12 + cheb12(i)*Rh(i)
                        res24 = res24 + cheb24(i)*Rh(i)
                    end do
                    do i = 14, 25
                        res24 = res24 + cheb24(i)*Rh(i)
                    end do
                end if
            else

                ! integr = 1  (or 3)

                call dqcheb(x, fval, cheb12, cheb24)
                do i = 1, 13
                    res12 = res12 + cheb12(i)*Rj(i)
                    res24 = res24 + cheb24(i)*Rj(i)
                end do
                do i = 14, 25
                    res24 = res24 + cheb24(i)*Rj(i)
                end do
                if (Integr /= 1) then

                    ! integr = 3

                    dc = log(Br - Bl)
                    Result = res24*dc
                    Abserr = abs((res24 - res12)*dc)
                    res12 = 0.0_wp
                    res24 = 0.0_wp
                    do i = 1, 13
                        res12 = res12 + cheb12(i)*Rh(i)
                        res24 = res24 + cheb24(i)*Rh(i)
                    end do
                    do i = 14, 25
                        res24 = res24 + cheb24(i)*Rh(i)
                    end do
                end if
            end if
            Result = (Result + res24)*factor
            Abserr = (Abserr + abs(res24 - res12))*factor
        else

            ! if a>bl and b<br, apply the 15-point gauss-kronrod
            ! scheme.

            ! dqwgts - external function subprogram defining
            ! the four possible weight functions

            call dqk15w(f, dqwgts, a, b, Alfa, Beta, Integr, Bl, Br, Result, Abserr, &
                        resabs, Resasc)
            Nev = 15
        end if

    end subroutine dqc25s
!********************************************************************************

!********************************************************************************
!>
!  chebyshev series expansion
!
!  this routine computes the chebyshev series expansion
!  of degrees 12 and 24 of a function using a
!  fast fourier transform method
!
!  * `f(x) = sum(k=1,..,13)` `(cheb12(k)*t(k-1,x))`
!  * `f(x) = sum(k=1,..,25)` `(cheb24(k)*t(k-1,x))`
!
!  where `t(k,x)` is the chebyshev polynomial of degree `k`.
!
!### See also
!  * [[dqc25c]], [[dqc25f]], [[dqc25s]]
!
!### History
!  * QUADPACK: revision date 830518 (yymmdd)

    subroutine dqcheb(x, Fval, Cheb12, Cheb24)
        implicit none

        real(wp), intent(in) :: x(11) !! vector of dimension 11 containing the
                                      !! values `cos(k*pi/24), k = 1, ..., 11`
        real(wp), intent(inout) :: Fval(25) !! vector of dimension 25 containing the
                                            !! function values at the points
                                            !! `(b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24`,
                                            !! where `(a,b)` is the approximation interval.
                                            !! `fval(1)` and `fval(25)` are divided by two
                                            !! (these values are destroyed at output).
        real(wp), intent(out) :: Cheb12(13) !! vector of dimension 13 containing the
                                            !! chebyshev coefficients for degree 12
        real(wp), intent(out) :: Cheb24(25) !! vector of dimension 25 containing the
                                            !! chebyshev coefficients for degree 24

        real(wp) :: alam, alam1, alam2, part1, part2, part3, v(12)
        integer :: i, j

        do i = 1, 12
            j = 26 - i
            v(i) = Fval(i) - Fval(j)
            Fval(i) = Fval(i) + Fval(j)
        end do
        alam1 = v(1) - v(9)
        alam2 = x(6)*(v(3) - v(7) - v(11))
        Cheb12(4) = alam1 + alam2
        Cheb12(10) = alam1 - alam2
        alam1 = v(2) - v(8) - v(10)
        alam2 = v(4) - v(6) - v(12)
        alam = x(3)*alam1 + x(9)*alam2
        Cheb24(4) = Cheb12(4) + alam
        Cheb24(22) = Cheb12(4) - alam
        alam = x(9)*alam1 - x(3)*alam2
        Cheb24(10) = Cheb12(10) + alam
        Cheb24(16) = Cheb12(10) - alam
        part1 = x(4)*v(5)
        part2 = x(8)*v(9)
        part3 = x(6)*v(7)
        alam1 = v(1) + part1 + part2
        alam2 = x(2)*v(3) + part3 + x(10)*v(11)
        Cheb12(2) = alam1 + alam2
        Cheb12(12) = alam1 - alam2
        alam = x(1)*v(2) + x(3)*v(4) + x(5)*v(6) + x(7)*v(8) + x(9)*v(10) &
               + x(11)*v(12)
        Cheb24(2) = Cheb12(2) + alam
        Cheb24(24) = Cheb12(2) - alam
        alam = x(11)*v(2) - x(9)*v(4) + x(7)*v(6) - x(5)*v(8) + x(3)*v(10) &
               - x(1)*v(12)
        Cheb24(12) = Cheb12(12) + alam
        Cheb24(14) = Cheb12(12) - alam
        alam1 = v(1) - part1 + part2
        alam2 = x(10)*v(3) - part3 + x(2)*v(11)
        Cheb12(6) = alam1 + alam2
        Cheb12(8) = alam1 - alam2
        alam = x(5)*v(2) - x(9)*v(4) - x(1)*v(6) - x(11)*v(8) + x(3)*v(10) &
               + x(7)*v(12)
        Cheb24(6) = Cheb12(6) + alam
        Cheb24(20) = Cheb12(6) - alam
        alam = x(7)*v(2) - x(3)*v(4) - x(11)*v(6) + x(1)*v(8) - x(9)*v(10) &
               - x(5)*v(12)
        Cheb24(8) = Cheb12(8) + alam
        Cheb24(18) = Cheb12(8) - alam
        do i = 1, 6
            j = 14 - i
            v(i) = Fval(i) - Fval(j)
            Fval(i) = Fval(i) + Fval(j)
        end do
        alam1 = v(1) + x(8)*v(5)
        alam2 = x(4)*v(3)
        Cheb12(3) = alam1 + alam2
        Cheb12(11) = alam1 - alam2
        Cheb12(7) = v(1) - v(5)
        alam = x(2)*v(2) + x(6)*v(4) + x(10)*v(6)
        Cheb24(3) = Cheb12(3) + alam
        Cheb24(23) = Cheb12(3) - alam
        alam = x(6)*(v(2) - v(4) - v(6))
        Cheb24(7) = Cheb12(7) + alam
        Cheb24(19) = Cheb12(7) - alam
        alam = x(10)*v(2) - x(6)*v(4) + x(2)*v(6)
        Cheb24(11) = Cheb12(11) + alam
        Cheb24(15) = Cheb12(11) - alam
        do i = 1, 3
            j = 8 - i
            v(i) = Fval(i) - Fval(j)
            Fval(i) = Fval(i) + Fval(j)
        end do
        Cheb12(5) = v(1) + x(8)*v(3)
        Cheb12(9) = Fval(1) - x(8)*Fval(3)
        alam = x(4)*v(2)
        Cheb24(5) = Cheb12(5) + alam
        Cheb24(21) = Cheb12(5) - alam
        alam = x(8)*Fval(2) - Fval(4)
        Cheb24(9) = Cheb12(9) + alam
        Cheb24(17) = Cheb12(9) - alam
        Cheb12(1) = Fval(1) + Fval(3)
        alam = Fval(2) + Fval(4)
        Cheb24(1) = Cheb12(1) + alam
        Cheb24(25) = Cheb12(1) - alam
        Cheb12(13) = v(1) - v(3)
        Cheb24(13) = Cheb12(13)
        alam = 1.0_wp/6.0_wp
        do i = 2, 12
            Cheb12(i) = Cheb12(i)*alam
        end do
        alam = 0.5_wp*alam
        Cheb12(1) = Cheb12(1)*alam
        Cheb12(13) = Cheb12(13)*alam
        do i = 2, 24
            Cheb24(i) = Cheb24(i)*alam
        end do
        Cheb24(1) = 0.5_wp*alam*Cheb24(1)
        Cheb24(25) = 0.5_wp*alam*Cheb24(25)

    end subroutine dqcheb
!********************************************************************************

!********************************************************************************
!>
!  the routine determines the limit of a given sequence of
!  approximations, by means of the epsilon algorithm of
!  p.wynn. an estimate of the absolute error is also given.
!  the condensed epsilon table is computed. only those
!  elements needed for the computation of the next diagonal
!  are preserved.
!
!### See also
!  *  [[dqagie]], [[dqagoe]], [[dqagpe]], [[dqagse]]
!
!### History
!  * QUADPACK: revision date 830518 (yymmdd).

    subroutine dqelg(n, Epstab, Result, Abserr, Res3la, Nres)
        implicit none

        integer, intent(inout) :: n !! epstab(n) contains the new element in the
                                    !! first column of the epsilon table.
        real(wp), intent(out) :: Abserr !! estimate of the absolute error computed from
                                        !! result and the 3 previous results
        real(wp), intent(inout) :: Epstab(limexp + 2) !! vector of dimension 52 containing the elements
                                                      !! of the two lower diagonals of the triangular
                                                      !! epsilon table. the elements are numbered
                                                      !! starting at the right-hand corner of the
                                                      !! triangle.
        real(wp), intent(out) :: Result !! resulting approximation to the integral
        real(wp), intent(inout) :: Res3la(3) !! vector of dimension 3 containing the last 3
                                             !! results
        integer, intent(inout) :: Nres !! number of calls to the routine
                                       !! (should be zero at first call)

        real(wp) :: delta1, delta2, delta3, epsinf, &
                    err1, err2, err3, e0, e1, e1abs, &
                    e2, e3, res, ss, tol1, tol2, tol3
        integer :: i, ib, ib2, ie, indx, k1, k2, k3, num

        integer :: newelm !! number of elements to be computed in the new diagonal
        real(wp) :: error !! `error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)`

        ! result is the element in the new diagonal with least value of error

        ! e0     - the 4 elements on which the computation of a new
        ! e1       element in the epsilon table is based
        ! e2
        ! e3                 e0
        !              e3    e1    new
        !                    e2

        Nres = Nres + 1
        Abserr = oflow
        Result = Epstab(n)
        if (n >= 3) then
            Epstab(n + 2) = Epstab(n)
            newelm = (n - 1)/2
            Epstab(n) = oflow
            num = n
            k1 = n
            do i = 1, newelm
                k2 = k1 - 1
                k3 = k1 - 2
                res = Epstab(k1 + 2)
                e0 = Epstab(k3)
                e1 = Epstab(k2)
                e2 = res
                e1abs = abs(e1)
                delta2 = e2 - e1
                err2 = abs(delta2)
                tol2 = max(abs(e2), e1abs)*epmach
                delta3 = e1 - e0
                err3 = abs(delta3)
                tol3 = max(e1abs, abs(e0))*epmach
                if (err2 > tol2 .or. err3 > tol3) then
                    e3 = Epstab(k1)
                    Epstab(k1) = e1
                    delta1 = e1 - e3
                    err1 = abs(delta1)
                    tol1 = max(e1abs, abs(e3))*epmach
                    ! if two elements are very close to each other, omit
                    ! a part of the table by adjusting the value of n
                    if (err1 > tol1 .and. err2 > tol2 .and. err3 > tol3) then
                        ss = 1.0_wp/delta1 + 1.0_wp/delta2 - 1.0_wp/delta3
                        epsinf = abs(ss*e1)
                        ! test to detect irregular behaviour in the table, and
                        ! eventually omit a part of the table adjusting the value
                        ! of n.
                        if (epsinf > 0.1e-03_wp) then
                            ! compute a new element and eventually adjust
                            ! the value of result.
                            res = e1 + 1.0_wp/ss
                            Epstab(k1) = res
                            k1 = k1 - 2
                            error = err2 + abs(res - e2) + err3
                            if (error <= Abserr) then
                                Abserr = error
                                Result = res
                            end if
                            cycle
                        end if
                    end if
                    n = i + i - 1
                    ! ***jump out of do-loop
                    exit
                else
                    ! if e0, e1 and e2 are equal to within machine
                    ! accuracy, convergence is assumed.
                    ! result = e2
                    ! abserr = abs(e1-e0)+abs(e2-e1)
                    Result = res
                    Abserr = err2 + err3
                    ! ***jump out of do-loop
                    Abserr = max(Abserr, 5.0_wp*epmach*abs(Result))
                    return
                end if
            end do

            ! shift the table.
            if (n == limexp) n = 2*(limexp/2) - 1
            ib = 1
            if ((num/2)*2 == num) ib = 2
            ie = newelm + 1
            do i = 1, ie
                ib2 = ib + 2
                Epstab(ib) = Epstab(ib2)
                ib = ib2
            end do
            if (num /= n) then
                indx = num - n + 1
                do i = 1, n
                    Epstab(i) = Epstab(indx)
                    indx = indx + 1
                end do
            end if
            if (Nres >= 4) then
                ! compute error estimate
                Abserr = abs(Result - Res3la(3)) + abs(Result - Res3la(2)) &
                         + abs(Result - Res3la(1))
                Res3la(1) = Res3la(2)
                Res3la(2) = Res3la(3)
                Res3la(3) = Result
            else
                Res3la(Nres) = Result
                Abserr = oflow
            end if
        end if

        Abserr = max(Abserr, 5.0_wp*epmach*abs(Result))

    end subroutine dqelg
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 15 point gauss-kronrod
!  rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)`, with error
!  estimate j = integral of `abs(f)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk15(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 15-point
                                        !! kronrod rule (resk) obtained by optimal addition
                                        !! of abscissae to the7-point gauss rule(resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should not exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))` over `(a,b)`

        real(wp) :: dhlgth, fc, fsum, fv1(7), fv2(7)
        integer :: j, jtw, jtwm1
        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 7-point gauss formula
        real(wp) :: resk !! result of the 15-point kronrod formula
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(4), parameter :: wg = [ &
                                             1.29484966168869693270611432679082018329e-1_wp, &
                                             2.79705391489276667901467771423779582487e-1_wp, &
                                             3.81830050505118944950369775488975133878e-1_wp, &
                                             4.17959183673469387755102040816326530612e-1_wp] !! weights of the 7-point gauss rule

        real(wp), dimension(8), parameter :: xgk = [ &
                                             9.91455371120812639206854697526328516642e-1_wp, &
                                             9.49107912342758524526189684047851262401e-1_wp, &
                                             8.64864423359769072789712788640926201211e-1_wp, &
                                             7.41531185599394439863864773280788407074e-1_wp, &
                                             5.86087235467691130294144838258729598437e-1_wp, &
                                             4.05845151377397166906606412076961463347e-1_wp, &
                                             2.07784955007898467600689403773244913480e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 15-point kronrod rule:
                                                                                            !!
                                                                                            !! * xgk(2), xgk(4), ... abscissae of the 7-point
                                                                                            !!   gauss rule
                                                                                            !! * xgk(1), xgk(3), ... abscissae which are optimally
                                                                                            !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             2.29353220105292249637320080589695919936e-2_wp, &
                                             6.30920926299785532907006631892042866651e-2_wp, &
                                             1.04790010322250183839876322541518017444e-1_wp, &
                                             1.40653259715525918745189590510237920400e-1_wp, &
                                             1.69004726639267902826583426598550284106e-1_wp, &
                                             1.90350578064785409913256402421013682826e-1_wp, &
                                             2.04432940075298892414161999234649084717e-1_wp, &
                                             2.09482141084727828012999174891714263698e-1_wp] !! weights of the 15-point kronrod rule

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 15-point kronrod approximation to
        ! the integral, and estimate the absolute error.

        fc = f(centr)
        resg = fc*wg(4)
        resk = fc*wgk(8)
        Resabs = abs(resk)
        do j = 1, 3
            jtw = j*2
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 4
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(8)*abs(fc - reskh)
        do j = 1, 7
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk15
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on (semi)infinite interval using a 15 point
!  gauss-kronrod quadrature rule, non-automatic
!
!  the original (infinite integration range is mapped
!  onto the interval (0,1) and (a,b) is a part of (0,1).
!  it is the purpose to compute:
!
!  * i = integral of transformed integrand over `(a,b)`,
!  * j = integral of abs(transformed integrand) over `(a,b)`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk15i(f, Boun, Inf, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: Boun !! finite bound of original integration
                                     !! range (set to zero if inf = +2)
        real(wp), intent(in) :: a !! lower limit for integration over subrange of `(0,1)`
        real(wp), intent(in) :: b !! upper limit for integration over subrange of `(0,1)`
        integer, intent(in) :: Inf !! * if inf = -1, the original interval is
                                   !!   `(-infinity,bound)`,
                                   !! * if inf = +1, the original interval is
                                   !!   `(bound,+infinity)`,
                                   !! * if inf = +2, the original interval is
                                   !!   `(-infinity,+infinity)` and
                                   !!
                                   !! the integral is computed as the sum of two
                                   !! integrals, one over `(-infinity,0)` and one over
                                   !! `(0,+infinity)`.
        real(wp), intent(out) :: Result !! approximation to the integral i.
                                        !! `result` is computed by applying the 15-point
                                        !! kronrod rule(resk) obtained by optimal addition
                                        !! of abscissae to the 7-point gauss rule(resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of
                                        !! `abs((transformed integrand)-i/(b-a))` over `(a,b)`

        real(wp) :: absc, dinf, fc, fsum, fv1(7), fv2(7)
        integer :: j
        real(wp) :: centr   !! mid point of the interval
        real(wp) :: hlgth   !! half-length of the interval
        real(wp) :: absc1   !! abscissa
        real(wp) :: absc2   !! abscissa
        real(wp) :: tabsc1  !! transformed abscissa
        real(wp) :: tabsc2  !! transformed abscissa
        real(wp) :: fval1   !! function value
        real(wp) :: fval2   !! function value
        real(wp) :: resg    !! result of the 7-point gauss formula
        real(wp) :: resk    !! result of the 15-point kronrod formula
        real(wp) :: reskh   !! approximation to the mean value of the transformed
                            !! integrand over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are supplied for the interval
        ! (-1,1).  because of symmetry only the positive abscissae and
        ! their corresponding weights are given.

        real(wp), dimension(8), parameter :: wg = [ &
                                             0.00000000000000000000000000000000000000e0_wp, &
                                             1.29484966168869693270611432679082018329e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp, &
                                             2.79705391489276667901467771423779582487e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp, &
                                             3.81830050505118944950369775488975133878e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp, &
                                             4.17959183673469387755102040816326530612e-1_wp] !! weights of the 7-point gauss rule, corresponding
                                                                                             !! to the abscissae `xgk(2), xgk(4), ...`.
                                                                                             !! `wg(1), wg(3), ...` are set to zero.

        real(wp), dimension(8), parameter :: xgk = [ &
                                             9.91455371120812639206854697526328516642e-1_wp, &
                                             9.49107912342758524526189684047851262401e-1_wp, &
                                             8.64864423359769072789712788640926201211e-1_wp, &
                                             7.41531185599394439863864773280788407074e-1_wp, &
                                             5.86087235467691130294144838258729598437e-1_wp, &
                                             4.05845151377397166906606412076961463347e-1_wp, &
                                             2.07784955007898467600689403773244913480e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 15-point kronrod rule:
                                                                                            !!
                                                                                            !! * xgk(2), xgk(4), ... abscissae of the 7-point
                                                                                            !!   gauss rule
                                                                                            !! * xgk(1), xgk(3), ... abscissae which are optimally
                                                                                            !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             2.29353220105292249637320080589695919936e-2_wp, &
                                             6.30920926299785532907006631892042866651e-2_wp, &
                                             1.04790010322250183839876322541518017444e-1_wp, &
                                             1.40653259715525918745189590510237920400e-1_wp, &
                                             1.69004726639267902826583426598550284106e-1_wp, &
                                             1.90350578064785409913256402421013682826e-1_wp, &
                                             2.04432940075298892414161999234649084717e-1_wp, &
                                             2.09482141084727828012999174891714263698e-1_wp] !! weights of the 15-point kronrod rule

        dinf = min(1, Inf)
        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        tabsc1 = Boun + dinf*(1.0_wp - centr)/centr
        fval1 = f(tabsc1)
        if (Inf == 2) fval1 = fval1 + f(-tabsc1)
        fc = (fval1/centr)/centr

        ! compute the 15-point kronrod approximation to
        ! the integral, and estimate the error.

        resg = wg(8)*fc
        resk = wgk(8)*fc
        Resabs = abs(resk)
        do j = 1, 7
            absc = hlgth*xgk(j)
            absc1 = centr - absc
            absc2 = centr + absc
            tabsc1 = Boun + dinf*(1.0_wp - absc1)/absc1
            tabsc2 = Boun + dinf*(1.0_wp - absc2)/absc2
            fval1 = f(tabsc1)
            fval2 = f(tabsc2)
            if (Inf == 2) then
                fval1 = fval1 + f(-tabsc1)
                fval2 = fval2 + f(-tabsc2)
            end if
            fval1 = (fval1/absc1)/absc1
            fval2 = (fval2/absc2)/absc2
            fv1(j) = fval1
            fv2(j) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(j)*fsum
            Resabs = Resabs + wgk(j)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(8)*abs(fc - reskh)
        do j = 1, 7
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resasc = Resasc*hlgth
        Resabs = Resabs*hlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk15i
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral with special singular weight functions using
!  a 15 point gauss-kronrod quadrature rule
!
!  to compute i = integral of `f*w` over `(a,b)`, with error
!  estimate j = integral of `abs(f*w)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 810101, revision date 830518 (yymmdd).

    subroutine dqk15w(f, w, p1, p2, p3, p4, Kp, a, b, Result, Abserr, Resabs, &
                      Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        procedure(weight_func) :: w !! function subprogram defining the integrand weight function `w(x)`.
        real(wp), intent(in) :: p1 !! parameter in the weight function
        real(wp), intent(in) :: p2 !! parameter in the weight function
        real(wp), intent(in) :: p3 !! parameter in the weight function
        real(wp), intent(in) :: p4 !! parameter in the weight function
        integer, intent(in) :: Kp !! key for indicating the type of weight function
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 15-point
                                        !! kronrod rule (resk) obtained by optimal addition
                                        !! of abscissae to the 7-point gauss rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral of `abs(f)`
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`

        real(wp) :: absc1, absc2, dhlgth, fc, fsum, fv1(7), fv2(7)
        integer :: j, jtw, jtwm1
        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 7-point gauss formula
        real(wp) :: resk !! result of the 15-point kronrod formula
        real(wp) :: reskh !! approximation to the mean value of f*w over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(8), parameter :: xgk = [ &
                                             9.91455371120812639206854697526328516642e-1_wp, &
                                             9.49107912342758524526189684047851262401e-1_wp, &
                                             8.64864423359769072789712788640926201211e-1_wp, &
                                             7.41531185599394439863864773280788407074e-1_wp, &
                                             5.86087235467691130294144838258729598437e-1_wp, &
                                             4.05845151377397166906606412076961463347e-1_wp, &
                                             2.07784955007898467600689403773244913480e-1_wp, &
                                             0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 15-point kronrod rule:
                                                                                            !!
                                                                                            !! * xgk(2), xgk(4), ... abscissae of the 7-point
                                                                                            !!   gauss rule
                                                                                            !! * xgk(1), xgk(3), ... abscissae which are optimally
                                                                                            !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             2.29353220105292249637320080589695919936e-2_wp, &
                                             6.30920926299785532907006631892042866651e-2_wp, &
                                             1.04790010322250183839876322541518017444e-1_wp, &
                                             1.40653259715525918745189590510237920400e-1_wp, &
                                             1.69004726639267902826583426598550284106e-1_wp, &
                                             1.90350578064785409913256402421013682826e-1_wp, &
                                             2.04432940075298892414161999234649084717e-1_wp, &
                                             2.09482141084727828012999174891714263698e-1_wp] !! weights of the 15-point kronrod rule

        real(wp), dimension(4), parameter :: wg = [ &
                                             1.29484966168869693270611432679082018329e-1_wp, &
                                             2.79705391489276667901467771423779582487e-1_wp, &
                                             3.81830050505118944950369775488975133878e-1_wp, &
                                             4.17959183673469387755102040816326530612e-1_wp] !! weights of the 7-point gauss rule

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 15-point kronrod approximation to the
        ! integral, and estimate the error.

        fc = f(centr)*w(centr, p1, p2, p3, p4, Kp)
        resg = wg(4)*fc
        resk = wgk(8)*fc
        Resabs = abs(resk)
        do j = 1, 3
            jtw = j*2
            absc = hlgth*xgk(jtw)
            absc1 = centr - absc
            absc2 = centr + absc
            fval1 = f(absc1)*w(absc1, p1, p2, p3, p4, Kp)
            fval2 = f(absc2)*w(absc2, p1, p2, p3, p4, Kp)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 4
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            absc1 = centr - absc
            absc2 = centr + absc
            fval1 = f(absc1)*w(absc1, p1, p2, p3, p4, Kp)
            fval2 = f(absc2)*w(absc2, p1, p2, p3, p4, Kp)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(8)*abs(fc - reskh)
        do j = 1, 7
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk15w
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 21 point
!  gauss-kronrod rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)`, with error
!  estimate j = integral of `abs(f)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk21(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 21-point
                                        !! kronrod rule (resk) obtained by optimal addition
                                        !! of abscissae to the 10-point gauss rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should not exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`
                                        !! over `(a,b)`

        real(wp) :: dhlgth, fc, fsum, fv1(10), fv2(10)
        integer :: j, jtw, jtwm1
        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 10-point gauss formula
        real(wp) :: resk !! result of the 21-point kronrod formula
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(5), parameter :: wg = [ &
                                             6.66713443086881375935688098933317928579e-2_wp, &
                                             1.49451349150580593145776339657697332403e-1_wp, &
                                             2.19086362515982043995534934228163192459e-1_wp, &
                                             2.69266719309996355091226921569469352860e-1_wp, &
                                             2.95524224714752870173892994651338329421e-1_wp] !! weights of the 10-point gauss rule

        real(wp), dimension(11), parameter :: xgk = [ &
                                              9.95657163025808080735527280689002847921e-1_wp, &
                                              9.73906528517171720077964012084452053428e-1_wp, &
                                              9.30157491355708226001207180059508346225e-1_wp, &
                                              8.65063366688984510732096688423493048528e-1_wp, &
                                              7.80817726586416897063717578345042377163e-1_wp, &
                                              6.79409568299024406234327365114873575769e-1_wp, &
                                              5.62757134668604683339000099272694140843e-1_wp, &
                                              4.33395394129247190799265943165784162200e-1_wp, &
                                              2.94392862701460198131126603103865566163e-1_wp, &
                                              1.48874338981631210884826001129719984618e-1_wp, &
                                              0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 21-point kronrod rule:
                                                                                             !!
                                                                                             !! * xgk(2), xgk(4), ...  abscissae of the 10-point
                                                                                             !!   gauss rule
                                                                                             !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                             !!   added to the 10-point gauss rule

        real(wp), dimension(11), parameter :: wgk = [ &
                                              1.16946388673718742780643960621920483962e-2_wp, &
                                              3.25581623079647274788189724593897606174e-2_wp, &
                                              5.47558965743519960313813002445801763737e-2_wp, &
                                              7.50396748109199527670431409161900093952e-2_wp, &
                                              9.31254545836976055350654650833663443900e-2_wp, &
                                              1.09387158802297641899210590325804960272e-1_wp, &
                                              1.23491976262065851077958109831074159512e-1_wp, &
                                              1.34709217311473325928054001771706832761e-1_wp, &
                                              1.42775938577060080797094273138717060886e-1_wp, &
                                              1.47739104901338491374841515972068045524e-1_wp, &
                                              1.49445554002916905664936468389821203745e-1_wp] !! weights of the 21-point kronrod rule

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 21-point kronrod approximation to
        ! the integral, and estimate the absolute error.

        resg = 0.0_wp
        fc = f(centr)
        resk = wgk(11)*fc
        Resabs = abs(resk)
        do j = 1, 5
            jtw = 2*j
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 5
            jtwm1 = 2*j - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(11)*abs(fc - reskh)
        do j = 1, 10
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk21
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 31 point
!  gauss-kronrod rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)` with error
!  estimate j = integral of `abs(f)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk31(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 31-point
                                        !! gauss-kronrod rule (resk), obtained by optimal
                                        !! addition of abscissae to the 15-point gauss
                                        !! rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the modulus,
                                        !! which should not exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`
                                        !! over `(a,b)`

        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 15-point gauss formula
        real(wp) :: resk !! result of the 31-point kronrod formula
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`
        real(wp) :: dhlgth, fc, fsum, fv1(15), fv2(15)
        integer :: j, jtw, jtwm1

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(8), parameter :: wg = [ &
                                             3.07532419961172683546283935772044177217e-2_wp, &
                                             7.03660474881081247092674164506673384667e-2_wp, &
                                             1.07159220467171935011869546685869303416e-1_wp, &
                                             1.39570677926154314447804794511028322521e-1_wp, &
                                             1.66269205816993933553200860481208811131e-1_wp, &
                                             1.86161000015562211026800561866422824506e-1_wp, &
                                             1.98431485327111576456118326443839324819e-1_wp, &
                                             2.02578241925561272880620199967519314839e-1_wp] !! weights of the 15-point gauss rule

        real(wp), dimension(16), parameter :: xgk = [ &
                                              9.98002298693397060285172840152271209073e-1_wp, &
                                              9.87992518020485428489565718586612581147e-1_wp, &
                                              9.67739075679139134257347978784337225283e-1_wp, &
                                              9.37273392400705904307758947710209471244e-1_wp, &
                                              8.97264532344081900882509656454495882832e-1_wp, &
                                              8.48206583410427216200648320774216851366e-1_wp, &
                                              7.90418501442465932967649294817947346862e-1_wp, &
                                              7.24417731360170047416186054613938009631e-1_wp, &
                                              6.50996741297416970533735895313274692547e-1_wp, &
                                              5.70972172608538847537226737253910641238e-1_wp, &
                                              4.85081863640239680693655740232350612866e-1_wp, &
                                              3.94151347077563369897207370981045468363e-1_wp, &
                                              2.99180007153168812166780024266388962662e-1_wp, &
                                              2.01194093997434522300628303394596207813e-1_wp, &
                                              1.01142066918717499027074231447392338787e-1_wp, &
                                              0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 31-point kronrod rule:
                                                                                             !!
                                                                                             !! * xgk(2), xgk(4), ...  abscissae of the 15-point
                                                                                             !!   gauss rule
                                                                                             !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                             !!   added to the 15-point gauss rule

        real(wp), dimension(16), parameter :: wgk = [ &
                                              5.37747987292334898779205143012764981831e-3_wp, &
                                              1.50079473293161225383747630758072680946e-2_wp, &
                                              2.54608473267153201868740010196533593973e-2_wp, &
                                              3.53463607913758462220379484783600481226e-2_wp, &
                                              4.45897513247648766082272993732796902233e-2_wp, &
                                              5.34815246909280872653431472394302967716e-2_wp, &
                                              6.20095678006706402851392309608029321904e-2_wp, &
                                              6.98541213187282587095200770991474757860e-2_wp, &
                                              7.68496807577203788944327774826590067221e-2_wp, &
                                              8.30805028231330210382892472861037896016e-2_wp, &
                                              8.85644430562117706472754436937743032123e-2_wp, &
                                              9.31265981708253212254868727473457185619e-2_wp, &
                                              9.66427269836236785051799076275893351367e-2_wp, &
                                              9.91735987217919593323931734846031310596e-2_wp, &
                                              1.00769845523875595044946662617569721916e-1_wp, &
                                              1.01330007014791549017374792767492546771e-1_wp] !! weights of the 31-point kronrod rule

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 31-point kronrod approximation to
        ! the integral, and estimate the absolute error.

        fc = f(centr)
        resg = wg(8)*fc
        resk = wgk(16)*fc
        Resabs = abs(resk)
        do j = 1, 7
            jtw = j*2
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 8
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(16)*abs(fc - reskh)
        do j = 1, 15
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk31
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 41 point
!  gauss-kronrod rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)`, with error
!  estimate j = integral of `abs(f)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk41(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 41-point
                                        !! gauss-kronrod rule (resk) obtained by optimal
                                        !! addition of abscissae to the 20-point gauss
                                        !! rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should not exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of abs(f-i/(b-a))
                                        !! over `(a,b)`

        real(wp) :: dhlgth, fc, fsum, fv1(20), fv2(20)
        integer :: j, jtw, jtwm1
        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 20-point gauss formula
        real(wp) :: resk !! result of the 41-point kronrod formula
        real(wp) :: reskh !! approximation to mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(10), parameter :: wg = [ &
                                              1.76140071391521183118619623518528163621e-2_wp, &
                                              4.06014298003869413310399522749321098791e-2_wp, &
                                              6.26720483341090635695065351870416063516e-2_wp, &
                                              8.32767415767047487247581432220462061002e-2_wp, &
                                              1.01930119817240435036750135480349876167e-1_wp, &
                                              1.18194531961518417312377377711382287005e-1_wp, &
                                              1.31688638449176626898494499748163134916e-1_wp, &
                                              1.42096109318382051329298325067164933035e-1_wp, &
                                              1.49172986472603746787828737001969436693e-1_wp, &
                                              1.52753387130725850698084331955097593492e-1_wp] !! weights of the 20-point gauss rule

        real(wp), dimension(21), parameter :: xgk = [ &
                                              9.98859031588277663838315576545863010000e-1_wp, &
                                              9.93128599185094924786122388471320278223e-1_wp, &
                                              9.81507877450250259193342994720216944567e-1_wp, &
                                              9.63971927277913791267666131197277221912e-1_wp, &
                                              9.40822633831754753519982722212443380274e-1_wp, &
                                              9.12234428251325905867752441203298113049e-1_wp, &
                                              8.78276811252281976077442995113078466711e-1_wp, &
                                              8.39116971822218823394529061701520685330e-1_wp, &
                                              7.95041428837551198350638833272787942959e-1_wp, &
                                              7.46331906460150792614305070355641590311e-1_wp, &
                                              6.93237656334751384805490711845931533386e-1_wp, &
                                              6.36053680726515025452836696226285936743e-1_wp, &
                                              5.75140446819710315342946036586425132814e-1_wp, &
                                              5.10867001950827098004364050955250998425e-1_wp, &
                                              4.43593175238725103199992213492640107840e-1_wp, &
                                              3.73706088715419560672548177024927237396e-1_wp, &
                                              3.01627868114913004320555356858592260615e-1_wp, &
                                              2.27785851141645078080496195368574624743e-1_wp, &
                                              1.52605465240922675505220241022677527912e-1_wp, &
                                              7.65265211334973337546404093988382110048e-2_wp, &
                                              0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 41-point gauss-kronrod rule:
                                                                                             !!
                                                                                             !! * xgk(2), xgk(4), ...  abscissae of the 20-point
                                                                                             !!   gauss rule
                                                                                             !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                             !!   added to the 20-point gauss rule

        real(wp), dimension(21), parameter :: wgk = [ &
                                              3.07358371852053150121829324603098748803e-3_wp, &
                                              8.60026985564294219866178795010234725213e-3_wp, &
                                              1.46261692569712529837879603088683561639e-2_wp, &
                                              2.03883734612665235980102314327547051228e-2_wp, &
                                              2.58821336049511588345050670961531429995e-2_wp, &
                                              3.12873067770327989585431193238007378878e-2_wp, &
                                              3.66001697582007980305572407072110084875e-2_wp, &
                                              4.16688733279736862637883059368947380440e-2_wp, &
                                              4.64348218674976747202318809261075168421e-2_wp, &
                                              5.09445739237286919327076700503449486648e-2_wp, &
                                              5.51951053482859947448323724197773291948e-2_wp, &
                                              5.91114008806395723749672206485942171364e-2_wp, &
                                              6.26532375547811680258701221742549805858e-2_wp, &
                                              6.58345971336184221115635569693979431472e-2_wp, &
                                              6.86486729285216193456234118853678017155e-2_wp, &
                                              7.10544235534440683057903617232101674129e-2_wp, &
                                              7.30306903327866674951894176589131127606e-2_wp, &
                                              7.45828754004991889865814183624875286161e-2_wp, &
                                              7.57044976845566746595427753766165582634e-2_wp, &
                                              7.63778676720807367055028350380610018008e-2_wp, &
                                              7.66007119179996564450499015301017408279e-2_wp] !! weights of the 41-point gauss-kronrod rule

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 41-point gauss-kronrod approximation to
        ! the integral, and estimate the absolute error.

        resg = 0.0_wp
        fc = f(centr)
        resk = wgk(21)*fc
        Resabs = abs(resk)
        do j = 1, 10
            jtw = j*2
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 10
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(21)*abs(fc - reskh)
        do j = 1, 20
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk41
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 51 point
!  gauss-kronrod rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)` with error
!  estimate j = integral of `abs(f)` over `(a,b)`
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk51(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subroutine defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i.
                                        !! `result` is computed by applying the 51-point
                                        !! kronrod rule (resk) obtained by optimal addition
                                        !! of abscissae to the 25-point gauss rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should not exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`
                                        !! over `(a,b)`

        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 25-point gauss formula
        real(wp) :: resk !! result of the 51-point kronrod formula
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        real(wp) :: dhlgth, fc, fsum, fv1(25), fv2(25)
        integer :: j, jtw, jtwm1

        ! the abscissae and weights are given for the interval (-1,1).
        ! because of symmetry only the positive abscissae and their
        ! corresponding weights are given.

        real(wp), dimension(13), parameter :: wg = [ &
                                              1.13937985010262879479029641132347736033e-2_wp, &
                                              2.63549866150321372619018152952991449360e-2_wp, &
                                              4.09391567013063126556234877116459536608e-2_wp, &
                                              5.49046959758351919259368915404733241601e-2_wp, &
                                              6.80383338123569172071871856567079685547e-2_wp, &
                                              8.01407003350010180132349596691113022902e-2_wp, &
                                              9.10282619829636498114972207028916533810e-2_wp, &
                                              1.00535949067050644202206890392685826988e-1_wp, &
                                              1.08519624474263653116093957050116619340e-1_wp, &
                                              1.14858259145711648339325545869555808641e-1_wp, &
                                              1.19455763535784772228178126512901047390e-1_wp, &
                                              1.22242442990310041688959518945851505835e-1_wp, &
                                              1.23176053726715451203902873079050142438e-1_wp] !! weights of the 25-point gauss rule

        real(wp), dimension(26), parameter :: xgk = [ &
                                              9.99262104992609834193457486540340593705e-1_wp, &
                                              9.95556969790498097908784946893901617258e-1_wp, &
                                              9.88035794534077247637331014577406227072e-1_wp, &
                                              9.76663921459517511498315386479594067745e-1_wp, &
                                              9.61614986425842512418130033660167241692e-1_wp, &
                                              9.42974571228974339414011169658470531905e-1_wp, &
                                              9.20747115281701561746346084546330631575e-1_wp, &
                                              8.94991997878275368851042006782804954175e-1_wp, &
                                              8.65847065293275595448996969588340088203e-1_wp, &
                                              8.33442628760834001421021108693569569461e-1_wp, &
                                              7.97873797998500059410410904994306569409e-1_wp, &
                                              7.59259263037357630577282865204360976388e-1_wp, &
                                              7.17766406813084388186654079773297780598e-1_wp, &
                                              6.73566368473468364485120633247622175883e-1_wp, &
                                              6.26810099010317412788122681624517881020e-1_wp, &
                                              5.77662930241222967723689841612654067396e-1_wp, &
                                              5.26325284334719182599623778158010178037e-1_wp, &
                                              4.73002731445714960522182115009192041332e-1_wp, &
                                              4.17885382193037748851814394594572487093e-1_wp, &
                                              3.61172305809387837735821730127640667422e-1_wp, &
                                              3.03089538931107830167478909980339329200e-1_wp, &
                                              2.43866883720988432045190362797451586406e-1_wp, &
                                              1.83718939421048892015969888759528415785e-1_wp, &
                                              1.22864692610710396387359818808036805532e-1_wp, &
                                              6.15444830056850788865463923667966312817e-2_wp, &
                                              0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 51-point kronrod rule
                                                                                             !!
                                                                                             !! * xgk(2), xgk(4), ...  abscissae of the 25-point
                                                                                             !!   gauss rule
                                                                                             !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                             !!   added to the 25-point gauss rule

        real(wp), dimension(26), parameter :: wgk = [ &
                                              1.98738389233031592650785188284340988943e-3_wp, &
                                              5.56193213535671375804023690106552207018e-3_wp, &
                                              9.47397338617415160720771052365532387165e-3_wp, &
                                              1.32362291955716748136564058469762380776e-2_wp, &
                                              1.68478177091282982315166675363363158404e-2_wp, &
                                              2.04353711458828354565682922359389736788e-2_wp, &
                                              2.40099456069532162200924891648810813929e-2_wp, &
                                              2.74753175878517378029484555178110786148e-2_wp, &
                                              3.07923001673874888911090202152285856009e-2_wp, &
                                              3.40021302743293378367487952295512032257e-2_wp, &
                                              3.71162714834155435603306253676198759960e-2_wp, &
                                              4.00838255040323820748392844670756464014e-2_wp, &
                                              4.28728450201700494768957924394951611020e-2_wp, &
                                              4.55029130499217889098705847526603930437e-2_wp, &
                                              4.79825371388367139063922557569147549836e-2_wp, &
                                              5.02776790807156719633252594334400844406e-2_wp, &
                                              5.23628858064074758643667121378727148874e-2_wp, &
                                              5.42511298885454901445433704598756068261e-2_wp, &
                                              5.59508112204123173082406863827473468203e-2_wp, &
                                              5.74371163615678328535826939395064719948e-2_wp, &
                                              5.86896800223942079619741758567877641398e-2_wp, &
                                              5.97203403241740599790992919325618538354e-2_wp, &
                                              6.05394553760458629453602675175654271623e-2_wp, &
                                              6.11285097170530483058590304162927119227e-2_wp, &
                                              6.14711898714253166615441319652641775865e-2_wp, &
                                              6.15808180678329350787598242400645531904e-2_wp] !! weights of the 51-point kronrod rule.

        centr = 0.5_wp*(a + b)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 51-point kronrod approximation to
        ! the integral, and estimate the absolute error.

        fc = f(centr)
        resg = wg(13)*fc
        resk = wgk(26)*fc
        Resabs = abs(resk)
        do j = 1, 12
            jtw = j*2
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 13
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(26)*abs(fc - reskh)
        do j = 1, 25
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk51
!********************************************************************************

!********************************************************************************
!>
!  estimate 1D integral on finite interval using a 61 point
!  gauss-kronrod rule and give error estimate, non-automatic
!
!  to compute i = integral of `f` over `(a,b)` with error
!  estimate j = integral of `abs(f)` over `(a,b)`.
!
!### History
!  * QUADPACK: date written 800101, revision date 830518 (yymmdd).

    subroutine dqk61(f, a, b, Result, Abserr, Resabs, Resasc)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! `result` is computed by applying the 61-point
                                        !! kronrod rule (resk) obtained by optimal addition of
                                        !! abscissae to the 30-point gauss rule (resg).
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(out) :: Resabs !! approximation to the integral j
        real(wp), intent(out) :: Resasc !! approximation to the integral of `abs(f-i/(b-a))`

        real(wp) :: dhlgth, fc, fsum, fv1(30), fv2(30)
        integer :: j, jtw, jtwm1
        real(wp) :: centr !! mid point of the interval
        real(wp) :: hlgth !! half-length of the interval
        real(wp) :: absc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 30-point gauss rule
        real(wp) :: resk !! result of the 61-point kronrod rule
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the
        ! interval (-1,1). because of symmetry only the positive
        ! abscissae and their corresponding weights are given.

        real(wp), dimension(15), parameter :: wg = [ &
                                              7.96819249616660561546588347467362245048e-3_wp, &
                                              1.84664683110909591423021319120472690962e-2_wp, &
                                              2.87847078833233693497191796112920436396e-2_wp, &
                                              3.87991925696270495968019364463476920332e-2_wp, &
                                              4.84026728305940529029381404228075178153e-2_wp, &
                                              5.74931562176190664817216894020561287971e-2_wp, &
                                              6.59742298821804951281285151159623612374e-2_wp, &
                                              7.37559747377052062682438500221907341538e-2_wp, &
                                              8.07558952294202153546949384605297308759e-2_wp, &
                                              8.68997872010829798023875307151257025768e-2_wp, &
                                              9.21225222377861287176327070876187671969e-2_wp, &
                                              9.63687371746442596394686263518098650964e-2_wp, &
                                              9.95934205867952670627802821035694765299e-2_wp, &
                                              1.01762389748405504596428952168554044633e-1_wp, &
                                              1.02852652893558840341285636705415043868e-1_wp] !! weights of the 30-point gauss rule

        real(wp), dimension(31), parameter :: xgk = [ &
                                              9.99484410050490637571325895705810819469e-1_wp, &
                                              9.96893484074649540271630050918695283341e-1_wp, &
                                              9.91630996870404594858628366109485724851e-1_wp, &
                                              9.83668123279747209970032581605662801940e-1_wp, &
                                              9.73116322501126268374693868423706884888e-1_wp, &
                                              9.60021864968307512216871025581797662930e-1_wp, &
                                              9.44374444748559979415831324037439121586e-1_wp, &
                                              9.26200047429274325879324277080474004086e-1_wp, &
                                              9.05573307699907798546522558925958319569e-1_wp, &
                                              8.82560535792052681543116462530225590057e-1_wp, &
                                              8.57205233546061098958658510658943856821e-1_wp, &
                                              8.29565762382768397442898119732501916439e-1_wp, &
                                              7.99727835821839083013668942322683240736e-1_wp, &
                                              7.67777432104826194917977340974503131695e-1_wp, &
                                              7.33790062453226804726171131369527645669e-1_wp, &
                                              6.97850494793315796932292388026640068382e-1_wp, &
                                              6.60061064126626961370053668149270753038e-1_wp, &
                                              6.20526182989242861140477556431189299207e-1_wp, &
                                              5.79345235826361691756024932172540495907e-1_wp, &
                                              5.36624148142019899264169793311072794164e-1_wp, &
                                              4.92480467861778574993693061207708795644e-1_wp, &
                                              4.47033769538089176780609900322854000162e-1_wp, &
                                              4.00401254830394392535476211542660633611e-1_wp, &
                                              3.52704725530878113471037207089373860654e-1_wp, &
                                              3.04073202273625077372677107199256553531e-1_wp, &
                                              2.54636926167889846439805129817805107883e-1_wp, &
                                              2.04525116682309891438957671002024709524e-1_wp, &
                                              1.53869913608583546963794672743255920419e-1_wp, &
                                              1.02806937966737030147096751318000592472e-1_wp, &
                                              5.14718425553176958330252131667225737491e-2_wp, &
                                              0.00000000000000000000000000000000000000e0_wp] !! abscissae of the 61-point kronrod rule:
                                                                                             !!
                                                                                             !! * `xgk(2), xgk(4)`  ... abscissae of the 30-point
                                                                                             !!   gauss rule
                                                                                             !! * `xgk(1), xgk(3)`  ... optimally added abscissae
                                                                                             !!   to the 30-point gauss rule

        real(wp), dimension(31), parameter :: wgk = [ &
                                              1.38901369867700762455159122675969968105e-3, &
                                              3.89046112709988405126720184451550327852e-3, &
                                              6.63070391593129217331982636975016813363e-3, &
                                              9.27327965951776342844114689202436042127e-3, &
                                              1.18230152534963417422328988532505928963e-2, &
                                              1.43697295070458048124514324435800101958e-2, &
                                              1.69208891890532726275722894203220923686e-2, &
                                              1.94141411939423811734089510501284558514e-2, &
                                              2.18280358216091922971674857383389934015e-2, &
                                              2.41911620780806013656863707252320267604e-2, &
                                              2.65099548823331016106017093350754143665e-2, &
                                              2.87540487650412928439787853543342111447e-2, &
                                              3.09072575623877624728842529430922726353e-2, &
                                              3.29814470574837260318141910168539275106e-2, &
                                              3.49793380280600241374996707314678750972e-2, &
                                              3.68823646518212292239110656171359677370e-2, &
                                              3.86789456247275929503486515322810502509e-2, &
                                              4.03745389515359591119952797524681142161e-2, &
                                              4.19698102151642461471475412859697577901e-2, &
                                              4.34525397013560693168317281170732580746e-2, &
                                              4.48148001331626631923555516167232437574e-2, &
                                              4.60592382710069881162717355593735805947e-2, &
                                              4.71855465692991539452614781810994864829e-2, &
                                              4.81858617570871291407794922983045926058e-2, &
                                              4.90554345550297788875281653672381736059e-2, &
                                              4.97956834270742063578115693799423285392e-2, &
                                              5.04059214027823468408930856535850289022e-2, &
                                              5.08817958987496064922974730498046918534e-2, &
                                              5.12215478492587721706562826049442082511e-2, &
                                              5.14261285374590259338628792157812598296e-2, &
                                              5.14947294294515675583404336470993075327e-2] !! weights of the 61-point kronrod rule

        centr = 0.5_wp*(b + a)
        hlgth = 0.5_wp*(b - a)
        dhlgth = abs(hlgth)

        ! compute the 61-point kronrod approximation to the
        ! integral, and estimate the absolute error.

        resg = 0.0_wp
        fc = f(centr)
        resk = wgk(31)*fc
        Resabs = abs(resk)
        do j = 1, 15
            jtw = j*2
            absc = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 15
            jtwm1 = j*2 - 1
            absc = hlgth*xgk(jtwm1)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum = fval1 + fval2
            resk = resk + wgk(jtwm1)*fsum
            Resabs = Resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do
        reskh = resk*0.5_wp
        Resasc = wgk(31)*abs(fc - reskh)
        do j = 1, 30
            Resasc = Resasc + wgk(j) &
                     *(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do
        Result = resk*hlgth
        Resabs = Resabs*dhlgth
        Resasc = Resasc*dhlgth
        Abserr = abs((resk - resg)*hlgth)
        if (Resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk61
!********************************************************************************

!********************************************************************************
!>
!  1D integration of `k`-th degree Chebyshev polynomial times a function with singularities
!
!  this routine computes modified chebsyshev moments. the `k`-th
!  modified chebyshev moment is defined as the integral over
!  `(-1,1)` of `w(x)*t(k,x)`, where `t(k,x)` is the chebyshev
!  polynomial of degree `k`.
!
!### History
!  * QUADPACK: date written 820101, revision date 830518 (yymmdd).

    subroutine dqmomo(Alfa, Beta, Ri, Rj, Rg, Rh, Integr)
        implicit none

        real(wp), intent(in) :: Alfa !! parameter in the weight function `w(x)`, `alfa>(-1)`
        real(wp), intent(in) :: Beta !! parameter in the weight function `w(x)`, `beta>(-1)`
        real(wp), intent(out) :: Ri(25) !! `i(k)` is the integral over (-1,1) of
                                        !! `(1+x)**alfa*t(k-1,x), k = 1, ..., 25`.
        real(wp), intent(out) :: Rj(25) !! `rj(k)` is the integral over (-1,1) of
                                        !! `(1-x)**beta*t(k-1,x), k = 1, ..., 25`.
        real(wp), intent(out) :: Rg(25) !! `rg(k)` is the integral over (-1,1) of
                                        !! `(1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25`.
        real(wp), intent(out) :: Rh(25) !! `rh(k)` is the integral over (-1,1) of
                                        !! `(1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25`.
        integer, intent(in) :: Integr !! input parameter indicating the modified
                                      !! moments to be computed:
                                      !!
                                      !! * integr = 1 compute `ri`, `rj`
                                      !! * integr = 2 compute `ri`, `rj`, `rg`
                                      !! * integr = 3 compute `ri`, `rj`, `rh`
                                      !! * integr = 4 compute `ri`, `rj`, `rg`, `rh`

        real(wp) :: alfp1, alfp2, an, anm1, betp1, betp2, ralf, rbet
        integer :: i, im1

        main : block

            alfp1 = Alfa + 1.0_wp
            betp1 = Beta + 1.0_wp
            alfp2 = Alfa + 2.0_wp
            betp2 = Beta + 2.0_wp
            ralf = 2.0_wp**alfp1
            rbet = 2.0_wp**betp1

            ! compute ri, rj using a forward recurrence relation.

            Ri(1) = ralf/alfp1
            Rj(1) = rbet/betp1
            Ri(2) = Ri(1)*Alfa/alfp2
            Rj(2) = Rj(1)*Beta/betp2
            an = 2.0_wp
            anm1 = 1.0_wp
            do i = 3, 25
                Ri(i) = -(ralf + an*(an - alfp2)*Ri(i - 1))/(anm1*(an + alfp1))
                Rj(i) = -(rbet + an*(an - betp2)*Rj(i - 1))/(anm1*(an + betp1))
                anm1 = an
                an = an + 1.0_wp
            end do
            if (Integr /= 1) then
                if (Integr /= 3) then

                    ! compute rg using a forward recurrence relation.

                    Rg(1) = -Ri(1)/alfp1
                    Rg(2) = -(ralf + ralf)/(alfp2*alfp2) - Rg(1)
                    an = 2.0_wp
                    anm1 = 1.0_wp
                    im1 = 2
                    do i = 3, 25
                        Rg(i) = -(an*(an - alfp2)*Rg(im1) - an*Ri(im1) + anm1*Ri(i)) &
                                /(anm1*(an + alfp1))
                        anm1 = an
                        an = an + 1.0_wp
                        im1 = i
                    end do
                    if (Integr == 2) exit main
                end if

                ! compute rh using a forward recurrence relation.

                Rh(1) = -Rj(1)/betp1
                Rh(2) = -(rbet + rbet)/(betp2*betp2) - Rh(1)
                an = 2.0_wp
                anm1 = 1.0_wp
                im1 = 2
                do i = 3, 25
                    Rh(i) = -(an*(an - betp2)*Rh(im1) - an*Rj(im1) + anm1*Rj(i)) &
                            /(anm1*(an + betp1))
                    anm1 = an
                    an = an + 1.0_wp
                    im1 = i
                end do
                do i = 2, 25, 2
                    Rh(i) = -Rh(i)
                end do
            end if

        end block main

        do i = 2, 25, 2
            Rj(i) = -Rj(i)
        end do

    end subroutine dqmomo
!********************************************************************************

!********************************************************************************
!>
!  1D non-adaptive automatic integrator
!
!  the routine calculates an approximation result to a
!  given definite integral i = integral of `f` over `(a,b)`,
!  hopefully satisfying following claim for accuracy
!  `abs(i-result)<=max(epsabs,epsrel*abs(i))`.
!
!### History
!  * QUADPACK: date written 800101, revision date 810101 (yymmdd),
!    kahaner,david,nbs - modified (2/82)

    subroutine dqng(f, a, b, Epsabs, Epsrel, Result, Abserr, Neval, Ier)
        implicit none

        procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accuracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                      !! if `epsabs<=0`
                                      !! and `epsrel<max(50*rel.mach.acc.,0.5e-28)`,
                                      !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral i
                                        !! result is obtained by applying the 21-point
                                        !! gauss-kronrod rule (res21) obtained by optimal
                                        !! addition of abscissae to the 10-point gauss rule
                                        !! (res10), or by applying the 43-point rule (res43)
                                        !! obtained by optimal addition of abscissae to the
                                        !! 21-point gauss-kronrod rule, or by applying the
                                        !! 87-point rule (res87) obtained by optimal addition
                                        !! of abscissae to the 43-point rule.
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        integer, intent(out) :: Neval !! number of integrand evaluations
        integer, intent(out) :: Ier !! * ier = 0 normal and reliable termination of the
                                    !!   routine. it is assumed that the requested
                                    !!   accuracy has been achieved.
                                    !! * ier>0 abnormal termination of the routine. it is
                                    !!   assumed that the requested accuracy has
                                    !!   not been achieved.
                                    !!
                                    !! error messages:
                                    !!
                                    !! * ier = 1 the maximum number of steps has been
                                    !!   executed. the integral is probably too
                                    !!   difficult to be calculated by dqng.
                                    !! * ier = 6 the input is invalid, because
                                    !!   `epsabs<=0` and
                                    !!   `epsrel<max(50*rel.mach.acc.,0.5e-28)`.
                                    !!   `result`, `abserr` and `neval` are set to zero.

        real(wp) :: dhlgth, fval1, fval2, fv1(5), fv2(5), fv3(5), fv4(5), reskh
        integer :: ipx, k, l
        real(wp) :: centr !! mid point of the integration interval
        real(wp) :: hlgth !! half-length of the integration interval
        real(wp) :: fcentr !! function value at mid point
        real(wp) :: absc !! abscissa
        real(wp) :: fval !! function value
        real(wp) :: savfun(21) !! array of function values which have already been computed
        real(wp) :: res10 !! 10-point gauss result
        real(wp) :: res21 !! 21-point kronrod result
        real(wp) :: res43 !! 43-point result
        real(wp) :: res87 !! 87-point result
        real(wp) :: resabs !! approximation to the integral of `abs(f)`
        real(wp) :: resasc !! approximation to the integral of `abs(f-i/(b-a))`

        ! the following data statements contain the
        ! abscissae and weights of the integration rules used.

        real(wp), dimension(5), parameter :: x1 = [ &
                                             9.73906528517171720077964012084452053428e-1_wp, &
                                             8.65063366688984510732096688423493048528e-1_wp, &
                                             6.79409568299024406234327365114873575769e-1_wp, &
                                             4.33395394129247190799265943165784162200e-1_wp, &
                                             1.48874338981631210884826001129719984618e-1_wp] !! abscissae common to the 10-, 21-, 43- and 87-point rule

        real(wp), dimension(5), parameter :: w10 = [ &
                                             6.66713443086881375935688098933317928579e-2_wp, &
                                             1.49451349150580593145776339657697332403e-1_wp, &
                                             2.19086362515982043995534934228163192459e-1_wp, &
                                             2.69266719309996355091226921569469352860e-1_wp, &
                                             2.95524224714752870173892994651338329421e-1_wp] !! weights of the 10-point formula

        real(wp), dimension(5), parameter :: x2 = [ &
                                             9.95657163025808080735527280689002847921e-1_wp, &
                                             9.30157491355708226001207180059508346225e-1_wp, &
                                             7.80817726586416897063717578345042377163e-1_wp, &
                                             5.62757134668604683339000099272694140843e-1_wp, &
                                             2.94392862701460198131126603103865566163e-1_wp] !! abscissae common to the 21-, 43- and 87-point rule

        real(wp), dimension(5), parameter :: w21a = [ &
                                             3.25581623079647274788189724593897606174e-2_wp, &
                                             7.50396748109199527670431409161900093952e-2_wp, &
                                             1.09387158802297641899210590325804960272e-1_wp, &
                                             1.34709217311473325928054001771706832761e-1_wp, &
                                             1.47739104901338491374841515972068045524e-1_wp] !! weights of the 21-point formula for abscissae x1

        real(wp), dimension(6), parameter :: w21b = [ &
                                             1.16946388673718742780643960621920483962e-2_wp, &
                                             5.47558965743519960313813002445801763737e-2_wp, &
                                             9.31254545836976055350654650833663443900e-2_wp, &
                                             1.23491976262065851077958109831074159512e-1_wp, &
                                             1.42775938577060080797094273138717060886e-1_wp, &
                                             1.49445554002916905664936468389821203745e-1_wp] !! weights of the 21-point formula for abscissae x2

        ! 43 and 87 coefficients are computed via the algorithm in the quadpack
        ! manual, section 2.2.2.
        !TODO: They need to be regenerated with the same precision as the others.
        real(wp), dimension(11), parameter :: x3 = [ &
                                              0.999333360901932081394099323919911_wp, &
                                              0.987433402908088869795961478381209_wp, &
                                              0.954807934814266299257919200290473_wp, &
                                              0.900148695748328293625099494069092_wp, &
                                              0.825198314983114150847066732588520_wp, &
                                              0.732148388989304982612354848755461_wp, &
                                              0.622847970537725238641159120344323_wp, &
                                              0.499479574071056499952214885499755_wp, &
                                              0.364901661346580768043989548502644_wp, &
                                              0.222254919776601296498260928066212_wp, &
                                              0.074650617461383322043914435796506_wp] !! abscissae common to the 43- and 87-point rule

        real(wp), dimension(10), parameter :: w43a = [ &
                                              0.016296734289666564924281974617663_wp, &
                                              0.037522876120869501461613795898115_wp, &
                                              0.054694902058255442147212685465005_wp, &
                                              0.067355414609478086075553166302174_wp, &
                                              0.073870199632393953432140695251367_wp, &
                                              0.005768556059769796184184327908655_wp, &
                                              0.027371890593248842081276069289151_wp, &
                                              0.046560826910428830743339154433824_wp, &
                                              0.061744995201442564496240336030883_wp, &
                                              0.071387267268693397768559114425516_wp] !! weights of the 43-point formula for abscissae x1, x3

        real(wp), dimension(12), parameter :: w43b = [ &
                                              0.001844477640212414100389106552965_wp, &
                                              0.010798689585891651740465406741293_wp, &
                                              0.021895363867795428102523123075149_wp, &
                                              0.032597463975345689443882222526137_wp, &
                                              0.042163137935191811847627924327955_wp, &
                                              0.050741939600184577780189020092084_wp, &
                                              0.058379395542619248375475369330206_wp, &
                                              0.064746404951445885544689259517511_wp, &
                                              0.069566197912356484528633315038405_wp, &
                                              0.072824441471833208150939535192842_wp, &
                                              0.074507751014175118273571813842889_wp, &
                                              0.074722147517403005594425168280423_wp] !! weights of the 43-point formula for abscissae x3

        real(wp), dimension(22), parameter :: x4 = [ &
                                              0.999902977262729234490529830591582_wp, &
                                              0.997989895986678745427496322365960_wp, &
                                              0.992175497860687222808523352251425_wp, &
                                              0.981358163572712773571916941623894_wp, &
                                              0.965057623858384619128284110607926_wp, &
                                              0.943167613133670596816416634507426_wp, &
                                              0.915806414685507209591826430720050_wp, &
                                              0.883221657771316501372117548744163_wp, &
                                              0.845710748462415666605902011504855_wp, &
                                              0.803557658035230982788739474980964_wp, &
                                              0.757005730685495558328942793432020_wp, &
                                              0.706273209787321819824094274740840_wp, &
                                              0.651589466501177922534422205016736_wp, &
                                              0.593223374057961088875273770349144_wp, &
                                              0.531493605970831932285268948562671_wp, &
                                              0.466763623042022844871966781659270_wp, &
                                              0.399424847859218804732101665817923_wp, &
                                              0.329874877106188288265053371824597_wp, &
                                              0.258503559202161551802280975429025_wp, &
                                              0.185695396568346652015917141167606_wp, &
                                              0.111842213179907468172398359241362_wp, &
                                              0.037352123394619870814998165437704_wp] !! abscissae of the 87-point rule

        real(wp), dimension(21), parameter :: w87a = [ &
                                              0.008148377384149172900002878448190_wp, &
                                              0.018761438201562822243935059003794_wp, &
                                              0.027347451050052286161582829741283_wp, &
                                              0.033677707311637930046581056957588_wp, &
                                              0.036935099820427907614589586742499_wp, &
                                              0.002884872430211530501334156248695_wp, &
                                              0.013685946022712701888950035273128_wp, &
                                              0.023280413502888311123409291030404_wp, &
                                              0.030872497611713358675466394126442_wp, &
                                              0.035693633639418770719351355457044_wp, &
                                              0.000915283345202241360843392549948_wp, &
                                              0.005399280219300471367738743391053_wp, &
                                              0.010947679601118931134327826856808_wp, &
                                              0.016298731696787335262665703223280_wp, &
                                              0.021081568889203835112433060188190_wp, &
                                              0.025370969769253827243467999831710_wp, &
                                              0.029189697756475752501446154084920_wp, &
                                              0.032373202467202789685788194889595_wp, &
                                              0.034783098950365142750781997949596_wp, &
                                              0.036412220731351787562801163687577_wp, &
                                              0.037253875503047708539592001191226_wp] !! weights of the 87-point formula for abscissae x1, x2, x3

        real(wp), dimension(23), parameter :: w87b = [ &
                                              0.000274145563762072350016527092881_wp, &
                                              0.001807124155057942948341311753254_wp, &
                                              0.004096869282759164864458070683480_wp, &
                                              0.006758290051847378699816577897424_wp, &
                                              0.009549957672201646536053581325377_wp, &
                                              0.012329447652244853694626639963780_wp, &
                                              0.015010447346388952376697286041943_wp, &
                                              0.017548967986243191099665352925900_wp, &
                                              0.019938037786440888202278192730714_wp, &
                                              0.022194935961012286796332102959499_wp, &
                                              0.024339147126000805470360647041454_wp, &
                                              0.026374505414839207241503786552615_wp, &
                                              0.028286910788771200659968002987960_wp, &
                                              0.030052581128092695322521110347341_wp, &
                                              0.031646751371439929404586051078883_wp, &
                                              0.033050413419978503290785944862689_wp, &
                                              0.034255099704226061787082821046821_wp, &
                                              0.035262412660156681033782717998428_wp, &
                                              0.036076989622888701185500318003895_wp, &
                                              0.036698604498456094498018047441094_wp, &
                                              0.037120549269832576114119958413599_wp, &
                                              0.037334228751935040321235449094698_wp, &
                                              0.037361073762679023410321241766599_wp] !! weights of the 87-point formula for abscissae x4

        ! test on validity of parameters

        Result = 0.0_wp
        Abserr = 0.0_wp
        Neval = 0
        Ier = 6
        if (Epsabs > 0.0_wp .or. Epsrel >= max(50.0_wp*epmach, 0.5e-28_wp)) &
            then
            hlgth = 0.5_wp*(b - a)
            dhlgth = abs(hlgth)
            centr = 0.5_wp*(b + a)
            fcentr = f(centr)
            Neval = 21
            Ier = 1

            ! compute the integral using the 10- and 21-point formula.

            do l = 1, 3
                select case (l)
                case (2)

                    ! compute the integral using the 43-point formula.

                    res43 = w43b(12)*fcentr
                    Neval = 43
                    do k = 1, 10
                        res43 = res43 + savfun(k)*w43a(k)
                    end do
                    do k = 1, 11
                        ipx = ipx + 1
                        absc = hlgth*x3(k)
                        fval = f(absc + centr) + f(centr - absc)
                        res43 = res43 + fval*w43b(k)
                        savfun(ipx) = fval
                    end do

                    ! test for convergence.

                    Result = res43*hlgth
                    Abserr = abs((res43 - res21)*hlgth)
                case (3)

                    ! compute the integral using the 87-point formula.

                    res87 = w87b(23)*fcentr
                    Neval = 87
                    do k = 1, 21
                        res87 = res87 + savfun(k)*w87a(k)
                    end do
                    do k = 1, 22
                        absc = hlgth*x4(k)
                        res87 = res87 + w87b(k)*(f(absc + centr) + f(centr - absc))
                    end do
                    Result = res87*hlgth
                    Abserr = abs((res87 - res43)*hlgth)
                case default
                    res10 = 0.0_wp
                    res21 = w21b(6)*fcentr
                    resabs = w21b(6)*abs(fcentr)
                    do k = 1, 5
                        absc = hlgth*x1(k)
                        fval1 = f(centr + absc)
                        fval2 = f(centr - absc)
                        fval = fval1 + fval2
                        res10 = res10 + w10(k)*fval
                        res21 = res21 + w21a(k)*fval
                        resabs = resabs + w21a(k)*(abs(fval1) + abs(fval2))
                        savfun(k) = fval
                        fv1(k) = fval1
                        fv2(k) = fval2
                    end do
                    ipx = 5
                    do k = 1, 5
                        ipx = ipx + 1
                        absc = hlgth*x2(k)
                        fval1 = f(centr + absc)
                        fval2 = f(centr - absc)
                        fval = fval1 + fval2
                        res21 = res21 + w21b(k)*fval
                        resabs = resabs + w21b(k)*(abs(fval1) + abs(fval2))
                        savfun(ipx) = fval
                        fv3(k) = fval1
                        fv4(k) = fval2
                    end do

                    ! test for convergence.

                    Result = res21*hlgth
                    resabs = resabs*dhlgth
                    reskh = 0.5_wp*res21
                    resasc = w21b(6)*abs(fcentr - reskh)
                    do k = 1, 5
                        resasc = resasc + w21a(k) &
                                 *(abs(fv1(k) - reskh) + abs(fv2(k) - reskh)) &
                                 + w21b(k) &
                                 *(abs(fv3(k) - reskh) + abs(fv4(k) - reskh))
                    end do
                    Abserr = abs((res21 - res10)*hlgth)
                    resasc = resasc*dhlgth
                end select
                if (resasc /= 0.0_wp .and. Abserr /= 0.0_wp) &
                    Abserr = resasc*min(1.0_wp, (200.0_wp*Abserr/resasc)**1.5_wp)
                if (resabs > uflow/(50.0_wp*epmach)) &
                    Abserr = max((epmach*50.0_wp)*resabs, Abserr)
                if (Abserr <= max(Epsabs, Epsrel*abs(Result))) Ier = 0
                ! ***jump out of do-loop
                if (Ier == 0) return
            end do
        end if
        call xerror('abnormal return from dqng ', Ier, 0)

    end subroutine dqng
!********************************************************************************

!********************************************************************************
!>
!  this routine maintains the descending ordering in the
!  list of the local error estimated resulting from the
!  interval subdivision process. at each call two error
!  estimates are inserted using the sequential search
!  method, top-down for the largest error estimate and
!  bottom-up for the smallest error estimate.
!
!### See also
!  *  [[dqage]], [[dqagie]], [[dqagpe]], [[dqawse]]
!
!### History
!  * QUADPACK: revision date 810101 (yymmdd)

    subroutine dqpsrt(Limit, Last, Maxerr, Ermax, Elist, Iord, Nrmax)
        implicit none

        integer, intent(in) :: Limit !! maximum number of error estimates the list can contain
        integer, intent(in) :: Last !! number of error estimates currently in the list
        integer, intent(inout) :: Maxerr !! `maxerr` points to the `nrmax`-th largest error
                                         !! estimate currently in the list
        real(wp), intent(out) :: Ermax !! `nrmax`-th largest error estimate
                                       !! `ermax = elist(maxerr)`
        real(wp), intent(in) :: Elist(Last) !! vector of dimension `last` containing
                                            !! the error estimates
        integer, intent(inout) :: Iord(Last) !! vector of dimension `last`, the first `k` elements
                                             !! of which contain pointers to the error
                                             !! estimates, such that
                                             !! `elist(iord(1)),...,  elist(iord(k))`
                                             !! form a decreasing sequence, with
                                             !! `k = last` if `last<=(limit/2+2)`, and
                                             !! `k = limit+1-last` otherwise
        integer, intent(inout) :: Nrmax !! `maxerr = iord(nrmax)`

        real(wp) :: errmax, errmin
        integer :: i, ibeg, ido, isucc, j, jbnd, jupbn, k

        main : block

            ! check whether the list contains more than
            ! two error estimates.

            if (Last > 2) then

                ! this part of the routine is only executed if, due to a
                ! difficult integrand, subdivision increased the error
                ! estimate. in the normal case the insert procedure should
                ! start after the nrmax-th largest error estimate.

                errmax = Elist(Maxerr)
                if (Nrmax /= 1) then
                    ido = Nrmax - 1
                    do i = 1, ido
                        isucc = Iord(Nrmax - 1)
                        ! ***jump out of do-loop
                        if (errmax <= Elist(isucc)) exit
                        Iord(Nrmax) = isucc
                        Nrmax = Nrmax - 1
                    end do
                end if

                ! compute the number of elements in the list to be maintained
                ! in descending order. this number depends on the number of
                ! subdivisions still allowed.

                jupbn = Last
                if (Last > (Limit/2 + 2)) jupbn = Limit + 3 - Last
                errmin = Elist(Last)

                ! insert errmax by traversing the list top-down,
                ! starting comparison from the element elist(iord(nrmax+1)).

                jbnd = jupbn - 1
                ibeg = Nrmax + 1
                if (ibeg <= jbnd) then
                    do i = ibeg, jbnd
                        isucc = Iord(i)
                        ! ***jump out of do-loop
                        if (errmax >= Elist(isucc)) then
                            ! insert errmin by traversing the list bottom-up.
                            Iord(i - 1) = Maxerr
                            k = jbnd
                            do j = i, jbnd
                                isucc = Iord(k)
                                ! ***jump out of do-loop
                                if (errmin < Elist(isucc)) then
                                    Iord(k + 1) = Last
                                    exit main
                                end if
                                Iord(k + 1) = isucc
                                k = k - 1
                            end do
                            Iord(i) = Last
                            exit main
                        end if
                        Iord(i - 1) = isucc
                    end do
                end if
                Iord(jbnd) = Maxerr
                Iord(jupbn) = Last
            else
                Iord(1) = 1
                Iord(2) = 2
            end if

        end block main

        ! set maxerr and ermax.
        Maxerr = Iord(Nrmax)
        Ermax = Elist(Maxerr)

    end subroutine dqpsrt
!********************************************************************************

!********************************************************************************
!>
!  this function subprogram is used together with the
!  routine [[qawc]] and defines the weight function.
!
!### See also
!  * [[dqk15w]]
!
!### History
!  * QUADPACK: revision date 810101 (yymmdd)
!
!### Keywords
!  * weight function, cauchy principal value

    real(wp) function dqwgtc(x, c, p2, p3, p4, Kp)
        implicit none

        real(wp), intent(in) :: c
        real(wp), intent(in) :: p2
        real(wp), intent(in) :: p3
        real(wp), intent(in) :: p4
        real(wp), intent(in) :: x
        integer, intent(in) :: Kp

        dqwgtc = 1.0_wp/(x - c)

    end function dqwgtc
!********************************************************************************

!********************************************************************************
!>
!  cos or sin in weight function
!
!### See also
!  * [[dqk15w]]
!
!### History
!  * QUADPACK: revision date 810101 (yymmdd)

    real(wp) function dqwgtf(x, Omega, p2, p3, p4, Integr)
        implicit none

        real(wp), intent(in) :: x
        real(wp), intent(in) :: Omega
        real(wp), intent(in) :: p2
        real(wp), intent(in) :: p3
        real(wp), intent(in) :: p4
        integer, intent(in) :: Integr

        if (Integr == 2) then
            dqwgtf = sin(Omega*x)
        else
            dqwgtf = cos(Omega*x)
        end if

    end function dqwgtf
!********************************************************************************

!********************************************************************************
!>
!  this function subprogram is used together with the
!  routine [[dqaws]] and defines the weight function.
!
!### See also
!  * [[dqk15w]]
!
!### History
!  * QUADPACK: revision date 810101 (yymmdd)

    real(wp) function dqwgts(x, a, b, Alfa, Beta, Integr)
        implicit none

        real(wp), intent(in) :: x
        real(wp), intent(in) :: a
        real(wp), intent(in) :: b
        real(wp), intent(in) :: Alfa
        real(wp), intent(in) :: Beta
        integer, intent(in) :: Integr

        real(wp) :: bmx, xma

        xma = x - a
        bmx = b - x

        dqwgts = xma**Alfa*bmx**Beta
        select case (Integr)
        case (1)
        case (3)
            dqwgts = dqwgts*log(bmx)
        case (4)
            dqwgts = dqwgts*log(xma)*log(bmx)
        case default
            dqwgts = dqwgts*log(xma)
        end select

    end function dqwgts
!********************************************************************************

!********************************************************************************
!>
!  dgtsl given a general tridiagonal matrix and a right hand
!  side will find the solution.
!
!### History
!  * linpack. this version dated 08/14/78.
!    jack dongarra, argonne national laboratory.

    subroutine dgtsl(n, c, d, e, b, info)
        implicit none

        integer, intent(in) :: n !! the order of the tridiagonal matrix.
        integer, intent(out) :: info !! * = 0 normal value.
                                     !! * = `k` if the `k`-th element of the diagonal becomes
                                     !!   exactly zero.  the subroutine returns when
                                     !!   this is detected.
        real(wp), intent(inout) :: c(n) !! the subdiagonal of the tridiagonal matrix.
                                        !! `c(2)` through `c(n) `should contain the subdiagonal.
                                        !! on output `c` is destroyed.
        real(wp), intent(inout) :: d(n) !! the diagonal of the tridiagonal matrix.
                                        !! on output `d` is destroyed.
        real(wp), intent(inout) :: e(n) !! the superdiagonal of the tridiagonal matrix.
                                        !! `e(1)` through `e(n-1)` should contain the superdiagonal.
                                        !! on output `e` is destroyed.
        real(wp), intent(inout) :: b(n) !! input: is the right hand side vector..
                                        !! output: the solution vector.

        integer :: k, kb, kp1, nm1, nm2
        real(wp) :: t

        info = 0
        c(1) = d(1)
        nm1 = n - 1

        if (nm1 >= 1) then
            d(1) = e(1)
            e(1) = 0.0_wp
            e(n) = 0.0_wp

            do k = 1, nm1
                kp1 = k + 1

                ! find the largest of the two rows

                if (abs(c(kp1)) >= abs(c(k))) then
                    ! interchange row
                    t = c(kp1)
                    c(kp1) = c(k)
                    c(k) = t
                    t = d(kp1)
                    d(kp1) = d(k)
                    d(k) = t
                    t = e(kp1)
                    e(kp1) = e(k)
                    e(k) = t
                    t = b(kp1)
                    b(kp1) = b(k)
                    b(k) = t
                end if

                ! zero elements
                if (c(k) == 0.0_wp) then
                    info = k
                    return
                end if

                t = -c(kp1)/c(k)
                c(kp1) = d(kp1) + t*d(k)
                d(kp1) = e(kp1) + t*e(k)
                e(kp1) = 0.0_wp
                b(kp1) = b(kp1) + t*b(k)
            end do

        end if

        if (c(n) == 0.0_wp) then
            info = n
        else
            ! back solve
            nm2 = n - 2
            b(n) = b(n)/c(n)
            if (n /= 1) then
                b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
                if (nm2 >= 1) then
                    do kb = 1, nm2
                        k = nm2 - kb + 1
                        b(k) = (b(k) - d(k)*b(k + 1) - e(k)*b(k + 2))/c(k)
                    end do
                end if
            end if
        end if

    end subroutine dgtsl
!********************************************************************************

!********************************************************************************
!>
!  This subroutine attempts to calculate the integral of `f(x)`
!  over the interval `a` to `b` with relative error not
!  exceeding `epsil`.
!
!  The result is obtained using a sequence of 1,3,7,15,31,63,
!  127, and 255 point interlacing formulae (no integrand
!  evaluations are wasted) of respective degree 1,5,11,23,
!  47,95,191 and 383. the formulae are based on the optimal
!  extension of the 3-point gauss formula.
!
!### See also
!  * Details of the formulae are given in "The optimum addition of points
!    to quadrature formulae" by t.n.l. patterson, maths. comp.
!    vol 22,847-856,1968.
!  * QUAD From [NSWC Mathematical Library](https://github.com/jacobwilliams/nswc)

subroutine dquad(f, a, b, result, epsil, npts, icheck)
    implicit none

    procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
    real(wp),intent(in) :: a !! lower limit of integration.
    real(wp),intent(in) :: b !! upper limit of integration.
    real(wp),intent(out) :: result !! the value of the integral to the
                                   !! specified relative accuracy.
    real(wp),intent(in) :: epsil !! relative accuracy required. when the relative
                                 !! difference of two successive formulae does not
                                 !! exceed `epsil` the last formula computed is taken
                                 !! as the result.
    integer,intent(out) :: npts !! number integrand evaluations.
    integer,intent(out) :: icheck !! on exit normally `icheck=0`. however if convergence
                                  !! to the accuracy requested is not achieved `icheck=1`
                                  !! on exit.

    real(wp) :: acum, diff, funct(127), fzero, sum, x
    integer :: i, inew, iold, j
    real(wp),dimension(8) :: results !! this array holds the results obtained by
                                     !! the 1,3,7, etc., point formulae. the number of
                                     !! formulae computed depends on `epsil`.
    integer :: k !! `results(k)` holds the value of the integral to the
                 !! specified relative accuracy.

    !>
    ! abscissae and weights of quadrature rules are stacked in
    ! array `p` in the order in which they are needed.
    real(wp),dimension(381),parameter :: p = [  7.74596669241483377035853079956479922167e-1_wp, &
                                                5.55555555555555555555555555555555555556e-1_wp, &
                                                8.88888888888888888888888888888888888889e-1_wp, &
                                                2.68488089868333440728569280666709624761e-1_wp, &
                                                9.60491268708020283423507092629079962670e-1_wp, &
                                                1.04656226026467265193823857192073038242e-1_wp, &
                                                4.34243749346802558002071502844627817283e-1_wp, &
                                                4.01397414775962222905051818618431878727e-1_wp, &
                                                4.50916538658474142345110087045570916539e-1_wp, &
                                                1.34415255243784220359968764802491520513e-1_wp, &
                                                5.16032829970797396969201205678609837136e-2_wp, &
                                                2.00628529376989021033931873331359306159e-1_wp, &
                                                9.93831963212755022208512841307951444370e-1_wp, &
                                                1.70017196299402603390274174026535252385e-2_wp, &
                                                8.88459232872256998890420167258502892651e-1_wp, &
                                                9.29271953151245376858942226541688263538e-2_wp, &
                                                6.21102946737226402940687443816594795012e-1_wp, &
                                                1.71511909136391380787353165019717217859e-1_wp, &
                                                2.23386686428966881628203986843998040091e-1_wp, &
                                                2.19156858401587496403693161643773747710e-1_wp, &
                                                2.25510499798206687386422549155949744906e-1_wp, &
                                                6.72077542959907035404010635813430091802e-2_wp, &
                                                2.58075980961766535646461187652328497046e-2_wp, &
                                                1.00314278611795578771293642695006079161e-1_wp, &
                                                8.43456573932110624631492964416019854788e-3_wp, &
                                                4.64628932617579865414046429639417161231e-2_wp, &
                                                8.57559200499903511541865204367976552400e-2_wp, &
                                                1.09578421055924638236688360572517068437e-1_wp, &
                                                9.99098124967667597662226062412998227686e-1_wp, &
                                                2.54478079156187441540278232983103810087e-3_wp, &
                                                9.81531149553740106867361888547025995016e-1_wp, &
                                                1.64460498543878109337883880689799875528e-2_wp, &
                                                9.29654857429740056670125725933373526769e-1_wp, &
                                                3.59571033071293220967778262209699862374e-2_wp, &
                                                8.36725938168868735502753818110221989775e-1_wp, &
                                                5.69795094941233574121973665457200316724e-2_wp, &
                                                7.02496206491527078609800156008001394343e-1_wp, &
                                                7.68796204990035310427051900809456411508e-2_wp, &
                                                5.31319743644375623972103438052468706781e-1_wp, &
                                                9.36271099812644736166587803392598658389e-2_wp, &
                                                3.31135393257976833092640782248746539410e-1_wp, &
                                                1.05669893580234809743815890442168534725e-1_wp, &
                                                1.12488943133186625745843327560318993879e-1_wp, &
                                                1.11956873020953456880143562321223860344e-1_wp, &
                                                1.12755256720768691607149869983804955967e-1_wp, &
                                                3.36038771482077305417339884731735403814e-2_wp, &
                                                1.29038001003512656259766532186329120125e-2_wp, &
                                                5.01571393058995374136795474239510758613e-2_wp, &
                                                4.21763044155885483908422682357386192911e-3_wp, &
                                                2.32314466399102694432564889365852548106e-2_wp, &
                                                4.28779600250077344929123037819815802239e-2_wp, &
                                                5.47892105279628650322175309941558213286e-2_wp, &
                                                1.26515655623006801137260909998182196593e-3_wp, &
                                                8.22300795723592966925778441546773952923e-3_wp, &
                                                1.79785515681282703328960466708609587502e-2_wp, &
                                                2.84897547458335486125060947723978716475e-2_wp, &
                                                3.84398102494555320386403467778787096784e-2_wp, &
                                                4.68135549906280124026480823343486642930e-2_wp, &
                                                5.28349467901165198620766563965308399269e-2_wp, &
                                                5.59784365104763194075533785872269074002e-2_wp, &
                                                9.99872888120357611937956782213944071260e-1_wp, &
                                                3.63221481845530659693580600240556307992e-4_wp, &
                                                9.97206259372221959076452532976228304987e-1_wp, &
                                                2.57904979468568827242779555856155526923e-3_wp, &
                                                9.88684757547429479938528919613635431554e-1_wp, &
                                                6.11550682211724633967828383326055155253e-3_wp, &
                                                9.72182874748581796578058835234688013989e-1_wp, &
                                                1.04982469096213218982728445836355320904e-2_wp, &
                                                9.46342858373402905148496208230196252152e-1_wp, &
                                                1.54067504665594978021308263315475287125e-2_wp, &
                                                9.10371156957004292497790670606627802042e-1_wp, &
                                                2.05942339159127111491885619503196295807e-2_wp, &
                                                8.63907938193690477146415857372833975090e-1_wp, &
                                                2.58696793272147469107582662448480815698e-2_wp, &
                                                8.06940531950217611856307980888497524441e-1_wp, &
                                                3.10735511116879648798843878245423584976e-2_wp, &
                                                7.39756044352694758677217797247847849281e-1_wp, &
                                                3.60644327807825726401071605896068916356e-2_wp, &
                                                6.62909660024780595461015255689389143141e-1_wp, &
                                                4.07155101169443189338940956005120803688e-2_wp, &
                                                5.77195710052045814843690955654189188852e-1_wp, &
                                                4.49145316536321974142542482618307358856e-2_wp, &
                                                4.83618026945841027562153280531749528761e-1_wp, &
                                                4.85643304066731987159471181667515286036e-2_wp, &
                                                3.83359324198730346916485193850312924770e-1_wp, &
                                                5.15832539520484587768091008575259100889e-2_wp, &
                                                2.77749822021824315065356412191446337302e-1_wp, &
                                                5.39054993352660639268769548863627639088e-2_wp, &
                                                1.68235251552207464982313275440102194714e-1_wp, &
                                                5.54814043565593639878384079955474248395e-2_wp, &
                                                5.63443130465927899719678607894467994099e-2_wp, &
                                                5.62776998312543012725953494255420385181e-2_wp, &
                                                5.63776283603847173876625571652345456628e-2_wp, &
                                                1.68019385741038652708694177373376419512e-2_wp, &
                                                6.45190005017573692280509776823864801062e-3_wp, &
                                                2.50785696529497687068397738442843404553e-2_wp, &
                                                2.10881524572663287933255325908005307552e-3_wp, &
                                                1.16157233199551347269849538868063638578e-2_wp, &
                                                2.14389800125038672464561593340624586806e-2_wp, &
                                                2.73946052639814325161087655093506901318e-2_wp, &
                                                6.32607319362633544219014096675880699298e-4_wp, &
                                                4.11150397865469304717026799389472424747e-3_wp, &
                                                8.98927578406413572328060374118804325340e-3_wp, &
                                                1.42448773729167743063415662436440605523e-2_wp, &
                                                1.92199051247277660193202803314218350072e-2_wp, &
                                                2.34067774953140062013240419700257395196e-2_wp, &
                                                2.64174733950582599310383282311985688836e-2_wp, &
                                                2.79892182552381597037766893004181239916e-2_wp, &
                                                1.80739564445388357820333919514772193888e-4_wp, &
                                                1.28952408261041739209850869778722441219e-3_wp, &
                                                3.05775341017553113613138395354134040323e-3_wp, &
                                                5.24912345480885912513384612635322646208e-3_wp, &
                                                7.70337523327974184816597819689326816907e-3_wp, &
                                                1.02971169579563555236864641070254134718e-2_wp, &
                                                1.29348396636073734547339558742365283615e-2_wp, &
                                                1.55367755558439824399284170162975429371e-2_wp, &
                                                1.80322163903912863200530999857265918070e-2_wp, &
                                                2.03577550584721594669470211177738968197e-2_wp, &
                                                2.24572658268160987071271218144441916129e-2_wp, &
                                                2.42821652033365993579735587740315274638e-2_wp, &
                                                2.57916269760242293884045503660307978007e-2_wp, &
                                                2.69527496676330319634384774240575382488e-2_wp, &
                                                2.77407021782796819939192039890754553228e-2_wp, &
                                                2.81388499156271506362976747068974890301e-2_wp, &
                                                9.99982430354891598580012135905109717915e-1_wp, &
                                                5.05360952078625176246656006337139648434e-5_wp, &
                                                9.99598799671910683251967529211801629987e-1_wp, &
                                                3.77746646326984660274364525157659292846e-4_wp, &
                                                9.98316635318407392530634580111074984770e-1_wp, &
                                                9.38369848542381500794044394681832138117e-4_wp, &
                                                9.95724104698407188509439459018460213288e-1_wp, &
                                                1.68114286542146990631373023491466618281e-3_wp, &
                                                9.91495721178106132398500079082519841189e-1_wp, &
                                                2.56876494379402037312771598563833315664e-3_wp, &
                                                9.85371499598520371113758241326513834962e-1_wp, &
                                                3.57289278351729964938448769864570199506e-3_wp, &
                                                9.77141514639705714156395810916629371363e-1_wp, &
                                                4.67105037211432174740543340826718946450e-3_wp, &
                                                9.66637851558416567092279836370846960853e-1_wp, &
                                                5.84344987583563950755951196450566504689e-3_wp, &
                                                9.53730006425761136414748643963112198908e-1_wp, &
                                                7.07248999543355546804631626841303341137e-3_wp, &
                                                9.38320397779592883654822310657872070243e-1_wp, &
                                                8.34283875396815770558412424167922936020e-3_wp, &
                                                9.20340025470012420729821382965612468142e-1_wp, &
                                                9.64117772970253669529830300284767390288e-3_wp, &
                                                8.99744899776940036638633212194468142956e-1_wp, &
                                                1.09557333878379016480327257363071595543e-2_wp, &
                                                8.76513414484705269741626645388423610417e-1_wp, &
                                                1.22758305600827700869663307413667617882e-2_wp, &
                                                8.50644494768350279757827407542049433990e-1_wp, &
                                                1.35915710097655467895729161814962317789e-2_wp, &
                                                8.22156254364980407372527142399375938309e-1_wp, &
                                                1.48936416648151820348103959267637767075e-2_wp, &
                                                7.91084933799848361434638057884175040395e-1_wp, &
                                                1.61732187295777199419479627980342182818e-2_wp, &
                                                7.57483966380513637926269606413039215349e-1_wp, &
                                                1.74219301594641737471522631397278549267e-2_wp, &
                                                7.21423085370098915484976184424530392547e-1_wp, &
                                                1.86318482561387901863140395332782911045e-2_wp, &
                                                6.82987431091079228087077605443637571318e-1_wp, &
                                                1.97954950480974994880277229389153128227e-2_wp, &
                                                6.42276642509759513774113624213729383798e-1_wp, &
                                                2.09058514458120238522218505878770859167e-2_wp, &
                                                5.99403930242242892974251049643553400441e-1_wp, &
                                                2.19563663053178249392605004207807929855e-2_wp, &
                                                5.54495132631932548866381362001869387185e-1_wp, &
                                                2.29409642293877487608005319195974357365e-2_wp, &
                                                5.07687757533716602154783137518047824630e-1_wp, &
                                                2.38540521060385400804460326687470805434e-2_wp, &
                                                4.59130011989832332873501971840246609692e-1_wp, &
                                                2.46905247444876769090608353528487841618e-2_wp, &
                                                4.08979821229888672409031653482169654497e-1_wp, &
                                                2.54457699654647658125743963445742965154e-2_wp, &
                                                3.57403837831532152376214925551056574778e-1_wp, &
                                                2.61156733767060976804988093771272602809e-2_wp, &
                                                3.04576441556714043335324049984830586514e-1_wp, &
                                                2.66966229274503599061546992881962515319e-2_wp, &
                                                2.50678730303483176612957105310757374530e-1_wp, &
                                                2.71855132296247918192086027320328453777e-2_wp, &
                                                1.95897502711100153915460230694341454649e-1_wp, &
                                                2.75797495664818730348687126189110696657e-2_wp, &
                                                1.40424233152560174593819634863430055039e-1_wp, &
                                                2.78772514766137016085237966902996263720e-2_wp, &
                                                8.44540400837108837101821672793851125821e-2_wp, &
                                                2.80764557938172466068478485336831566215e-2_wp, &
                                                2.81846489497456943393973278703614550567e-2_wp, &
                                                2.81763190330166021306535805326311346689e-2_wp, &
                                                2.81888141801923586938312785882097958145e-2_wp, &
                                                8.40096928705193263543470886866882097559e-3_wp, &
                                                3.22595002508786846140254888664674399963e-3_wp, &
                                                1.25392848264748843534198869221421702276e-2_wp, &
                                                1.05440762286331677224956681256723093434e-3_wp, &
                                                5.80786165997756736349247694340318193339e-3_wp, &
                                                1.07194900062519336232280796670312293403e-2_wp, &
                                                1.36973026319907162580543827546753450659e-2_wp, &
                                                3.16303660822264476886001542319765673695e-4_wp, &
                                                2.05575198932734652358557179891967892437e-3_wp, &
                                                4.49463789203206786164030187059407895106e-3_wp, &
                                                7.12243868645838715317078312182203027615e-3_wp, &
                                                9.60995256236388300966014016571091750361e-3_wp, &
                                                1.17033887476570031006620209850128697598e-2_wp, &
                                                1.32087366975291299655191641155992844418e-2_wp, &
                                                1.39946091276190798518883446502090619958e-2_wp, &
                                                9.03727346587511492612048292799447801127e-5_wp, &
                                                6.44762041305724779327197260132661244643e-4_wp, &
                                                1.52887670508776556838105789798193451205e-3_wp, &
                                                2.62456172740442956256692394303736650452e-3_wp, &
                                                3.85168761663987092408298909845685878251e-3_wp, &
                                                5.14855847897817776184323205351270717418e-3_wp, &
                                                6.46741983180368672736697793711826418082e-3_wp, &
                                                7.76838777792199121996420850814877146855e-3_wp, &
                                                9.01610819519564316002654999286329590351e-3_wp, &
                                                1.01788775292360797334735105588869484098e-2_wp, &
                                                1.12286329134080493535635609072220958065e-2_wp, &
                                                1.21410826016682996789867793870157637319e-2_wp, &
                                                1.28958134880121146942022751830153989003e-2_wp, &
                                                1.34763748338165159817192387120287691244e-2_wp, &
                                                1.38703510891398409969596019945377276614e-2_wp, &
                                                1.40694249578135753181488373534487445151e-2_wp, &
                                                2.51578703842806614886029901874368269190e-5_wp, &
                                                1.88873264506504913660930569062668820773e-4_wp, &
                                                4.69184924247850409754566477203398287419e-4_wp, &
                                                8.40571432710722463646844648204542489678e-4_wp, &
                                                1.28438247189701017680511226368885244509e-3_wp, &
                                                1.78644639175864982468103287043436779759e-3_wp, &
                                                2.33552518605716087370269795035052675936e-3_wp, &
                                                2.92172493791781975377975593711547903293e-3_wp, &
                                                3.53624499771677773402315813405234465284e-3_wp, &
                                                4.17141937698407885279206212083887894115e-3_wp, &
                                                4.82058886485126834764915150142383212497e-3_wp, &
                                                5.47786669391895082401636286815357973430e-3_wp, &
                                                6.13791528004138504348316537068338089359e-3_wp, &
                                                6.79578550488277339478645809074811588946e-3_wp, &
                                                7.44682083240759101740519796338188835376e-3_wp, &
                                                8.08660936478885997097398139901710914092e-3_wp, &
                                                8.71096507973208687357613156986392746334e-3_wp, &
                                                9.31592412806939509315701976663914555223e-3_wp, &
                                                9.89774752404874974401386146945765641137e-3_wp, &
                                                1.04529257229060119261109252939385429584e-2_wp, &
                                                1.09781831526589124696302502103903964927e-2_wp, &
                                                1.14704821146938743804002659597987178682e-2_wp, &
                                                1.19270260530192700402230163343735402717e-2_wp, &
                                                1.23452623722438384545304176764243920809e-2_wp, &
                                                1.27228849827323829062871981722871482577e-2_wp, &
                                                1.30578366883530488402494046885636301405e-2_wp, &
                                                1.33483114637251799530773496440981257659e-2_wp, &
                                                1.35927566148123959096043013660164226889e-2_wp, &
                                                1.37898747832409365174343563094555348329e-2_wp, &
                                                1.39386257383068508042618983451498131860e-2_wp, &
                                                1.40382278969086233034239242668415783107e-2_wp, &
                                                1.40881595165083010653267902663155673344e-2_wp, &
                                                9.99997596379748464620231592559093837611e-1_wp, &
                                                6.93793643241082671695382297169979368601e-6_wp, &
                                                9.99943996207054375763853646470050626596e-1_wp, &
                                                5.32752936697806131253524393895881823770e-5_wp, &
                                                9.99760490924432047330447933438138365417e-1_wp, &
                                                1.35754910949228719729842895656339874910e-4_wp, &
                                                9.99380338025023581928079338774322759519e-1_wp, &
                                                2.49212400482997294024537662868023009356e-4_wp, &
                                                9.98745614468095114703528542397791959986e-1_wp, &
                                                3.89745284473282293215563879845838727539e-4_wp, &
                                                9.97805354495957274561833338685736105778e-1_wp, &
                                                5.54295314930374714917732120266906130439e-4_wp, &
                                                9.96514145914890273848684083613153803279e-1_wp, &
                                                7.40282804244503330463160177700222594979e-4_wp, &
                                                9.94831502800621000519130529785414200225e-1_wp, &
                                                9.45361516858525382463015198607451979300e-4_wp, &
                                                9.92721344282788615328202203758497351413e-1_wp, &
                                                1.16748411742995940769333157872940045783e-3_wp, &
                                                9.90151370400770159180535140748087193102e-1_wp, &
                                                1.40490799565514464271521123296916900291e-3_wp, &
                                                9.87092527954034067189898792468859039993e-1_wp, &
                                                1.65611272815445260521682786451135109534e-3_wp, &
                                                9.83518657578632728761664630770795617152e-1_wp, &
                                                1.91971297101387241252271734466970358673e-3_wp, &
                                                9.79406281670862683806133521363753397925e-1_wp, &
                                                2.19440692536383883880291840868628867052e-3_wp, &
                                                9.74734459752402667760726712997609707570e-1_wp, &
                                                2.47895822665756793067821535745476374906e-3_wp, &
                                                9.69484659502459231770908123207442170150e-1_wp, &
                                                2.77219576459345099399521424961083418592e-3_wp, &
                                                9.63640621569812132520974048832142316972e-1_wp, &
                                                3.07301843470257832340783765226605973620e-3_wp, &
                                                9.57188216109860962736208621751374728884e-1_wp, &
                                                3.38039799108692038234993039038885672945e-3_wp, &
                                                9.50115297521294876557842262038304179472e-1_wp, &
                                                3.69337791702565081825729998764452535617e-3_wp, &
                                                9.42411565191083059812560025758972247897e-1_wp, &
                                                4.01106872407502339888993614903965571565e-3_wp, &
                                                9.34068436157725787999477771530264179420e-1_wp, &
                                                4.33264096809298285453769983324695296414e-3_wp, &
                                                9.25078932907075652364132996222672693491e-1_wp, &
                                                4.65731729975685477727794484849624969667e-3_wp, &
                                                9.15437587155765040643953616154536973514e-1_wp, &
                                                4.98436456476553860120001022162080486896e-3_wp, &
                                                9.05140358813261595189303779754262290451e-1_wp, &
                                                5.31308660518705656628804340372923963811e-3_wp, &
                                                8.94184568335559022859352159222674193953e-1_wp, &
                                                5.64281810138444415845460587311671071412e-3_wp, &
                                                8.82568840247341906841695404228946666934e-1_wp, &
                                                5.97291956550816580494729856935913899149e-3_wp, &
                                                8.70293055548113905851151444154923420039e-1_wp, &
                                                6.30277344908575871716398763418949052534e-3_wp, &
                                                8.57358310886232156525126596087163923324e-1_wp, &
                                                6.63178124290188789412200734180398266358e-3_wp, &
                                                8.43766882672708601038314138625718101532e-1_wp, &
                                                6.95936140939042293944507544479114448976e-3_wp, &
                                                8.29522194637401400178105088351227616660e-1_wp, &
                                                7.28494798055380706387981147534993110085e-3_wp, &
                                                8.14628787655137413435816577891367083540e-1_wp, &
                                                7.60798966571905658321739694223386579593e-3_wp, &
                                                7.99092290960841401799803164024282388556e-1_wp, &
                                                7.92794933429484911025254235115728574858e-3_wp, &
                                                7.82919394118283016385180478369806362244e-1_wp, &
                                                8.24430376303286803055059706535356438929e-3_wp, &
                                                7.66117819303760090716674093891474570508e-1_wp, &
                                                8.55654356130768961917293275004918273728e-3_wp, &
                                                7.48696293616936602822828737479369222926e-1_wp, &
                                                8.86417320948249426411429453091759055196e-3_wp, &
                                                7.30664521242181261329306715350070027793e-1_wp, &
                                                9.16671116356078840670519648472888628456e-3_wp, &
                                                7.12033155362252034586679081013994469857e-1_wp, &
                                                9.46368999383006529427243113943215866506e-3_wp, &
                                                6.92813769779114702894651485928486730921e-1_wp, &
                                                9.75465653631741146108293452735497379607e-3_wp, &
                                                6.73018830230418479198879472689545414663e-1_wp, &
                                                1.00391720440568407981810290438378080094e-2_wp, &
                                                6.52661665410017496100770934689234627423e-1_wp, &
                                                1.03168123309476216819207000244181912440e-2_wp, &
                                                6.31756437711194230413584623172536712454e-1_wp, &
                                                1.05871679048851979309428189932402399185e-2_wp, &
                                                6.10318113715186400155578672320162394224e-1_wp, &
                                                1.08498440893373140990245263318076192187e-2_wp, &
                                                5.88362434447662541434367386275547111879e-1_wp, &
                                                1.11044611340069265369994188454572096386e-2_wp, &
                                                5.65905885423654422622970392231343950219e-1_wp, &
                                                1.13506543159805966017344840804968802477e-2_wp, &
                                                5.42965666498311490492303133422203430532e-1_wp, &
                                                1.15880740330439525684239776012385794172e-2_wp, &
                                                5.19559661537457021992914143047305013398e-1_wp, &
                                                1.18163858908302357632247900084966241627e-2_wp, &
                                                4.95706407918761460170111534008667847416e-1_wp, &
                                                1.20352707852795626304498694306103606393e-2_wp, &
                                                4.71425065871658876934088018252224136473e-1_wp, &
                                                1.22444249816119858986292063324627371480e-2_wp, &
                                                4.46735387662028473742222281592907967623e-1_wp, &
                                                1.24435601907140352631495031087115129475e-2_wp, &
                                                4.21657686626163300056304726883310969563e-1_wp, &
                                                1.26324036435420787645405441085200317588e-2_wp, &
                                                3.96212806057615939182521394284924513267e-1_wp, &
                                                1.28106981638773619668417039218064387909e-2_wp, &
                                                3.70422087950078230137537383958155880174e-1_wp, &
                                                1.29782022395373992858421803348245496762e-2_wp, &
                                                3.44307341599438022776622416041385263462e-1_wp, &
                                                1.31346900919601528363813260381779658443e-2_wp, &
                                                3.17890812068476683181739338725980798218e-1_wp, &
                                                1.32799517439305306503775089710281336690e-2_wp, &
                                                2.91195148518246681963691099017626573079e-1_wp, &
                                                1.34137930851100985129663776085717215632e-2_wp, &
                                                2.64243372410926761944948292977628978728e-1_wp, &
                                                1.35360359349562136136653091890522717067e-2_wp, &
                                                2.37058845589829727212668030348623871778e-1_wp, &
                                                1.36465181025712914283998912158692590540e-2_wp, &
                                                2.09665238243181194766342717964439602895e-1_wp, &
                                                1.37450934430018966322520540025550273779e-2_wp, &
                                                1.82086496759252198246399488588060039322e-1_wp, &
                                                1.38316319095064286764959688535114143323e-2_wp, &
                                                1.54346811481378108692446779987579230421e-1_wp, &
                                                1.39060196013254612635312215253609885781e-2_wp, &
                                                1.26470584372301966850663538758563345841e-1_wp, &
                                                1.39681588065169385157277797674326721757e-2_wp, &
                                                9.84823965981192020902757578971386695319e-2_wp, &
                                                1.40179680394566088098722249688496041850e-2_wp, &
                                                7.04069760428551790632968760555968372924e-2_wp, &
                                                1.40553820726499642771679253311023986914e-2_wp, &
                                                4.22691647653636032124048988444769492564e-2_wp, &
                                                1.40803519625536613248458411104536513059e-2_wp, &
                                                1.40938864107824626141884882355263630430e-2_wp, &
                                                1.40928450691604083549592735386756230351e-2_wp, &
                                                1.40944070900961793469156392941048979072e-2_wp ]

    icheck = 0

    ! check for trivial case.
    if (a == b) then
        ! trivial case
        result = 0.0_wp
        npts = 0
        return
    else
        ! scale factors.
        sum = (b + a)/2.0_wp
        diff = (b - a)/2.0_wp
        ! 1-point gauss
        fzero = f(sum)
        results(1) = 2.0_wp*fzero*diff
        i = 0
        iold = 0
        inew = 1
        k = 2
        acum = 0.0_wp
        do
            ! contribution from new function values.
            iold = iold + inew
            do j = inew, iold
                i = i + 1
                x = p(i)*diff
                funct(j) = f(sum + x) + f(sum - x)
                i = i + 1
                acum = acum + p(i)*funct(j)
            end do
            inew = iold + 1
            i = i + 1
            results(k) = (acum + p(i)*fzero)*diff
            ! check for convergence.
            if (abs(results(k) - results(k - 1)) <= epsil*abs(results(k))) exit
            if (k == 8) then
                ! convergence not achieved.
                icheck = 1
                exit
            else
                k = k + 1
                acum = 0.0_wp
                ! contribution from function values already computed.
                do j = 1, iold
                    i = i + 1
                    acum = acum + p(i)*funct(j)
                end do
            end if
        end do
        result = results(k)
    end if

    ! normal termination.
    npts = inew + iold

    end subroutine dquad
!********************************************************************************

!********************************************************************************
!>
!  Integrate a function tabulated at arbitrarily spaced
!  abscissas using overlapping parabolas.
!
!  DAVINT integrates a function tabulated at arbitrarily spaced
!  abscissas.  The limits of integration need not coincide
!  with the tabulated abscissas.
!
!  A method of overlapping parabolas fitted to the data is used
!  provided that there are at least 3 abscissas between the
!  limits of integration.  DAVINT also handles two special cases.
!  If the limits of integration are equal, DAVINT returns a
!  result of zero regardless of the number of tabulated values.
!  If there are only two function values, DAVINT uses the
!  trapezoid rule.
!
!### References
!  * R. E. Jones, "Approximate integrator of functions
!    tabulated at arbitrarily spaced abscissas",
!    Report SC-M-69-335, Sandia Laboratories, 1969.
!  * Original program from *Numerical Integration* by Davis & Rabinowitz
!    Adaptation and modifications by Rondall E Jones.
!
!### Author
!  * Jones, R. E., (SNLA)
!
!### Revision history
!  * 690901  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, Jan 2022 : modernized this procedure.

    subroutine davint(x,y,n,xlo,xup,ans,ierr)

    implicit none

    real(wp),dimension(:),intent(in) :: x !! array of abscissas, which must be in increasing order.
    real(wp),dimension(:),intent(in) :: y !! array of function values. i.e., `y(i)=func(x(i))`
    integer,intent(in) :: n !! The integer number of function values supplied.
                            !! `N >= 2` unless `XLO = XUP`.
    real(wp),intent(in) :: xlo !! lower limit of integration
    real(wp),intent(in) :: xup !! upper limit of integration.  Must have `XLO <= XUP`
    real(wp),intent(out) :: ans !! computed approximate value of integral
    integer,intent(out) :: ierr !! A status code:
                                !!
                                !! * Normal Code
                                !!    * =1 Means the requested integration was performed.
                                !! * Abnormal Codes
                                !!    * =2 Means `XUP` was less than `XLO`.
                                !!    * =3 Means the number of `X(I)` between `XLO` and `XUP`
                                !!      (inclusive) was less than 3 and neither of the two
                                !!      special cases described in the abstract occurred.
                                !!      No integration was performed.
                                !!    * =4 Means the restriction `X(I+1)>X(I)` was violated.
                                !!    * =5 Means the number `N` of function values was < 2.
                                !!
                                !! ANS is set to zero if `IERR` = 2, 3, 4, or 5.

    integer :: i , inlft , inrt , istart , istop
    real(wp) :: a , b , c , ca , cb , cc , fl , fr , r3 , &
                rp5 , slope , sum , syl , syl2 , syl3 , syu , &
                syu2 , syu3 , term1 , term2 , term3 , x1 , &
                x12 , x13 , x2 , x23 , x3

    ierr = 1
    ans = 0.0_wp

    ! error checks and trivial cases:
    if (xlo == xup) return
    if (xlo > xup) then
        ierr = 2
        call xerror('the upper limit of integration was not greater '//&
                    'than the lower limit.',4,1)
        return
    end if
    if (n < 2) then
        ierr = 5
        call xerror('less than two function values were supplied.', &
                     4,1)
        return
    end if
    do i = 2 , n
        if ( x(i)<=x(i-1) ) then
            ierr = 4
            call xerror('the abscissas were not strictly increasing.  must have '&
                        //'x(i-1) < x(i) for all i.',4,1)
            return
        end if
        if ( x(i)>xup ) exit
    enddo

    if ( n<3 ) then

        ! special n=2 case
        slope = (y(2)-y(1))/(x(2)-x(1))
        fl = y(1) + slope*(xlo-x(1))
        fr = y(2) + slope*(xup-x(2))
        ans = 0.5_wp*(fl+fr)*(xup-xlo)

    elseif ( x(n-2)<xlo ) then

        ierr = 3
        call xerror('there were less than three function values '&
                    //'between the limits of integration.',4,1)

    elseif ( x(3)<=xup ) then

        i = 1
        do
            if ( x(i)>=xlo ) then
                inlft = i
                i = n
                do
                    if ( x(i)<=xup ) then
                        inrt = i
                        if ( (inrt-inlft)>=2 ) then
                            istart = inlft
                            if ( inlft==1 ) istart = 2
                            istop = inrt
                            if ( inrt==n ) istop = n - 1
                            r3 = 3.0_wp
                            rp5 = 0.5_wp
                            sum = 0.0_wp
                            syl = xlo
                            syl2 = syl*syl
                            syl3 = syl2*syl
                            do i = istart , istop
                                x1 = x(i-1)
                                x2 = x(i)
                                x3 = x(i+1)
                                x12 = x1 - x2
                                x13 = x1 - x3
                                x23 = x2 - x3
                                term1 = y(i-1)/(x12*x13)
                                term2 = -y(i)/(x12*x23)
                                term3 = y(i+1)/(x13*x23)
                                a = term1 + term2 + term3
                                b = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3
                                c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
                                if ( i>istart ) then
                                    ca = 0.5_wp*(a+ca)
                                    cb = 0.5_wp*(b+cb)
                                    cc = 0.5_wp*(c+cc)
                                else
                                    ca = a
                                    cb = b
                                    cc = c
                                endif
                                syu = x2
                                syu2 = syu*syu
                                syu3 = syu2*syu
                                sum = sum + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2) + cc*(syu-syl)
                                ca = a
                                cb = b
                                cc = c
                                syl = syu
                                syl2 = syu2
                                syl3 = syu3
                            enddo
                            syu = xup
                            ans = sum + ca*(syu**3-syl3)/r3 + cb*rp5*(syu**2-syl2) + cc*(syu-syl)
                        else
                            ierr = 3
                            call xerror('there were less than three function values '&
                                        //'between the limits of integration.',4,1)
                        endif
                        return
                    endif
                    i = i - 1
                end do
            endif
            i = i + 1
        end do

    else
        ierr = 3
        call xerror('there were less than three function values '&
                    //'between the limits of integration.',4,1)
    endif

    end subroutine davint
!********************************************************************************

!********************************************************************************
!>
!  Integrate a function using a 7-point adaptive Newton-Cotes
!  quadrature rule.
!
!  DQNC79 is a general purpose program for evaluation of
!  one dimensional integrals of user defined functions.
!  DQNC79 will pick its own points for evaluation of the
!  integrand and these will vary from problem to problem.
!  Thus, DQNC79 is not designed to integrate over data sets.
!  Moderately smooth integrands will be integrated efficiently
!  and reliably.  For problems with strong singularities,
!  oscillations etc., the user may wish to use more sophis-
!  ticated routines such as those in QUADPACK.  One measure
!  of the reliability of DQNC79 is the output parameter `K`,
!  giving the number of integrand evaluations that were needed.
!
!### Author
!  * Kahaner, D. K., (NBS)
!  * Jones, R. E., (SNLA)
!
!### Revision history  (YYMMDD)
!  * 790601  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890911  Removed unnecessary intrinsics.  (WRB)
!  * 890911  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920218  Code redone to parallel QNC79.  (WRB)
!  * 930120  Increase array size 80->99, and KMX 2000->5000 for SUN -r8 wordlength.  (RWC)
!  * Jacob Williams, Jan 2022 : modernized the SLATEC procedure. added quad-precision coefficients.
!
!@note This one has a lot of failures in the test cases.

    subroutine dqnc79(fun,a,b,err,ans,ierr,k)

    implicit none

    procedure(func) :: fun !! function subprogram defining the integrand function `f(x)`.
    real(wp),intent(in) :: a !! lower limit of integration
    real(wp),intent(in) :: b !! upper limit of integration (may be less than `A`)
    real(wp),intent(in) :: err !! a requested error tolerance.  Normally, pick a value
                               !! `0 < ERR < 1.0e-8`.
    real(wp),intent(out) :: ans !! computed value of the integral.  Hopefully, `ANS` is
                                !! accurate to within `ERR *` integral of `ABS(FUN(X))`.
    integer,intent(out) :: ierr !! a status code:
                                !!
                                !!  * Normal codes
                                !!    * **1** `ANS` most likely meets requested error tolerance.
                                !!    * **-1** `A` and `B` are too nearly equal to
                                !!      allow normal integration. `ANS` is set to zero.
                                !!  * Abnormal code
                                !!    * **2**  `ANS` probably does not meet requested error tolerance.
    integer,intent(out) :: k !! the number of function evaluations actually used to do
                             !! the integration.  A value of `K > 1000` indicates a
                             !! difficult problem; other programs may be more efficient.
                             !! `DQNC79` will gracefully give up if `K` exceeds 5000.

    real(wp),parameter :: w1 = 41.0_wp/140.0_wp
    real(wp),parameter :: w2 = 216.0_wp/140.0_wp
    real(wp),parameter :: w3 = 27.0_wp/140.0_wp
    real(wp),parameter :: w4 = 272.0_wp/140.0_wp
    real(wp),parameter :: sq2 = sqrt(2.0_wp)
    real(wp),parameter :: ln2 = log(2.0_wp)
    integer,parameter :: nbits = int(d1mach(5)*i1mach14/0.30102000_wp)  !! is 0.30102000 supposed to be log10(2.0_wp) ???
    integer,parameter :: nlmx = min(99, (nbits*4)/5)
    integer,parameter :: nlmn = 2
    integer,parameter :: kml = 7
    integer,parameter :: kmx = 5000      !! JW : is this the max function evals? should be an input
    integer,parameter :: array_size = 99 !! JW : what is this magic number 99 array size ??
                                         !! does it depend on the number of function evals ?
                                         !! (see comment in revision history)

    real(wp) :: ae,area,bank,blocal,c,ce,ee,ef,eps,q13,q7,q7l,test,tol,vr
    integer :: i,l,lmn,lmx,nib
    real(wp),dimension(13) :: f
    real(wp),dimension(array_size) :: aa,f1,f2,f3,f4,f5,f6,f7,hh,q7r,vl
    integer,dimension(array_size) :: lr

    ans = 0.0_wp
    ierr = 1
    if ( a==b ) return ! JW : this was an error return in the original code

    ce = 0.0_wp
    lmx = nlmx
    lmn = nlmn
    if ( b/=0.0_wp ) then
        if ( sign(1.0_wp,b)*a>0.0_wp ) then
            c = abs(1.0_wp-a/b)
            if ( c<=0.1_wp ) then
                if ( c<=0.0_wp ) then
                    ierr = -1
                    call xerror('a and b are too nearly equal to allow normal integration. '&
                                //'ans is set to zero and ierr to -1.',-1,-1)
                    return
                end if
                nib = 0.5_wp - log(c)/ln2
                lmx = min(nlmx,nbits-nib-4)
                if ( lmx<2 ) then
                    call xerror('a and b are too nearly equal to allow normal integration. '&
                                //'ans is set to zero and ierr to -1.',-1,-1)
                    return
                end if
                lmn = min(lmn,lmx)
            endif
        endif
    endif
    tol = max(abs(err),2.0_wp**(5-nbits))
    if ( err==0.0_wp ) tol = sqrt(epmach)
    eps = tol
    hh(1) = (b-a)/12.0_wp
    aa(1) = a
    lr(1) = 1
    do i = 1 , 11 , 2
        f(i) = fun(a+(i-1)*hh(1))
    enddo
    blocal = b
    f(13) = fun(blocal)
    k = 7
    l = 1
    area = 0.0_wp
    q7 = 0.0_wp
    ef = 256.0_wp/255.0_wp
    bank = 0.0_wp

    loop : do

        ! compute refined estimates, estimate the error, etc.
        do i = 2 , 12 , 2
            f(i) = fun(aa(l)+(i-1)*hh(l))
        enddo
        k = k + 6

        ! compute left and right half estimates
        q7l = hh(l)*((w1*(f(1)+f(7))+w2*(f(2)+f(6)))+(w3*(f(3)+f(5))+w4*f(4)))
        q7r(l) = hh(l)*((w1*(f(7)+f(13))+w2*(f(8)+f(12)))+(w3*(f(9)+f(11))+w4*f(10)))

        ! update estimate of integral of absolute value
        area = area + (abs(q7l)+abs(q7r(l))-abs(q7))

        ! do not bother to test convergence before minimum refinement level
        if ( l>=lmn ) then

            ! estimate the error in new value for whole interval, q13
            q13 = q7l + q7r(l)
            ee = abs(q7-q13)*ef

            ! compute nominal allowed error
            ae = eps*area

            ! borrow from bank account, but not too much
            test = min(ae+0.8_wp*bank,10.0_wp*ae)

            ! don't ask for excessive accuracy
            test = max(test,tol*abs(q13),0.00003_wp*tol*area)   ! jw : should change ?

            ! now, did this interval pass or not?
            if ( ee<=test ) then
                ! on good intervals accumulate the theoretical estimate
                ce = ce + (q7-q13)/255.0_wp
            else
                ! consider the left half of next deeper level
                if ( k>kmx ) lmx = min(kml,lmx)
                if ( l<lmx ) then
                    call f200()
                    cycle loop
                end if
                ! have hit maximum refinement level -- penalize the cumulative error
                ce = ce + (q7-q13)
            endif

            ! update the bank account.  don't go into debt.
            bank = bank + (ae-ee)
            if ( bank<0.0_wp ) bank = 0.0_wp

            ! did we just finish a left half or a right half?
            if ( lr(l)<=0 ) then
                ! proceed to right half at this level
                vl(l) = q13
                call f300()
                cycle loop
            else
                ! left and right halves are done, so go back up a level
                vr = q13
                do
                    if ( l<=1 ) then
                        !   exit
                        ans = vr
                        if ( abs(ce)>2.0_wp*tol*area ) then
                            ierr = 2
                            call xerror('ans is probably insufficiently accurate.',2,1)
                        endif
                        return
                    else
                        if ( l<=17 ) ef = ef*sq2
                        eps = eps*2.0_wp
                        l = l - 1
                        if ( lr(l)<=0 ) then
                            vl(l) = vl(l+1) + vr
                            call f300()
                            cycle loop
                        else
                            vr = vl(l+1) + vr
                        endif
                    endif
                end do
            endif
        endif

        call f200()

    end do loop

    contains

        subroutine f200()
            l = l + 1
            eps = eps*0.5_wp
            if ( l<=17 ) ef = ef/sq2
            hh(l) = hh(l-1)*0.5_wp
            lr(l) = -1
            aa(l) = aa(l-1)
            q7 = q7l
            f1(l) = f(7)
            f2(l) = f(8)
            f3(l) = f(9)
            f4(l) = f(10)
            f5(l) = f(11)
            f6(l) = f(12)
            f7(l) = f(13)
            f(13) = f(7)
            f(11) = f(6)
            f(9) = f(5)
            f(7) = f(4)
            f(5) = f(3)
            f(3) = f(2)
        end subroutine f200

        subroutine f300()
            q7 = q7r(l-1)
            lr(l) = 1
            aa(l) = aa(l) + 12.0_wp*hh(l)
            f(1) = f1(l)
            f(3) = f2(l)
            f(5) = f3(l)
            f(7) = f4(l)
            f(9) = f5(l)
            f(11) = f6(l)
            f(13) = f7(l)
        end subroutine f300

    end subroutine dqnc79
!********************************************************************************

!********************************************************************************
!>
!  Integrate a real function of one variable over a finite
!  interval using an adaptive 8-point Legendre-Gauss
!  algorithm.
!
!  Intended primarily for high accuracy
!  integration or integration of smooth functions.
!
!### See also
!  * Original SLATEC sourcecode from: http://www.netlib.org/slatec/src/dgaus8.f
!
!### History
!  * Author: Jones, R. E., (SNLA)
!  * 810223  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890911  Removed unnecessary intrinsics.  (WRB)
!  * 890911  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * Jacob Williams : Jan 2022 : refactored SLATEC routine to modern Fortran.

    subroutine dgauss8( f, a, b, error_tol, ans, ierr, err)

    implicit none

    procedure(func) :: f !! function subprogram defining the integrand function `f(x)`.
    real(wp),intent(in)   :: a          !! lower bound of the integration
    real(wp),intent(in)   :: b          !! upper bound of the integration
    real(wp),intent(in)   :: error_tol  !! is a requested pseudorelative error tolerance.  normally
                                        !! pick a value of abs(error_tol) so that
                                        !! `dtol < abs(error_tol) <= 1.0e-3` where dtol is the larger
                                        !! of `1.0e-18 `and the real unit roundoff `d1mach(4)`.
                                        !! `ans` will normally have no more error than `abs(error_tol)`
                                        !! times the integral of the absolute value of `f(x)`.  usually,
                                        !! smaller values of error_tol yield more accuracy and require
                                        !! more function evaluations.
    real(wp),intent(out)  :: ans        !! computed value of integral
    integer,intent(out)   :: ierr       !! status code:
                                        !!
                                        !!  * normal codes:
                                        !!    * 1 : `ans` most likely meets requested error tolerance,
                                        !!      or `a=b`.
                                        !!    * -1 : `a` and `b` are too nearly equal to allow normal
                                        !!      integration. `ans` is set to zero.
                                        !!  * abnormal code:
                                        !!    * 2 : `ans` probably does not meet requested error tolerance.
    real(wp),intent(out)  :: err        !! an estimate of the absolute error in `ans`.
                                        !! the estimated error is solely for information to the user and
                                        !! should not be used as a correction to the computed integral.

    ! note: see also dqnc79 for some clues about the purpose of these numbers...
    real(wp),parameter  :: sq2    = sqrt(2.0_wp)
    real(wp),parameter  :: ln2    = log(2.0_wp)
    integer,parameter   :: kmx    = 5000
    integer,parameter   :: kml    = 6
    real(wp),parameter  :: magic  = 0.30102000_wp   !! is 0.30102000 supposed to be log10(2.0_wp) ???
    integer,parameter   :: iwork  = 60              !! size of the work arrays. ?? Why 60 ??
    integer,parameter   :: nbits  = int(d1mach(5)*i1mach14/magic)
    integer,parameter   :: nlmn   = 1
    integer,parameter   :: nlmx   = min(60,(nbits*5)/8)

    integer :: k !! number of function evaluations
    integer                   :: l,lmn,lmx,mxl,nib
    real(wp)                  :: ae,area,c,ee,ef,eps,est,gl,glr,tol
    real(wp),dimension(iwork) :: aa,hh,vl,gr
    integer,dimension(iwork)  :: lr

    ans = 0.0_wp
    ierr = 1
    err = 0.0_wp
    if (a == b) return

    aa = 0.0_wp
    hh = 0.0_wp
    vl = 0.0_wp
    gr = 0.0_wp
    lr = 0
    lmx = nlmx
    lmn = nlmn
    if (b /= 0.0_wp) then
        if (sign(1.0_wp,b)*a > 0.0_wp) then
            c = abs(1.0_wp-a/b)
            if (c <= 0.1_wp) then
                if (c <= 0.0_wp) return
                nib = int(0.5_wp - log(c)/ln2)
                lmx = min(nlmx,nbits-nib-7)
                if (lmx < 1) then
                    ! a and b are too nearly equal to allow
                    ! normal integration [ans is set to zero]
                    ierr = -1
                    return
                end if
                lmn = min(lmn,lmx)
            end if
        end if
    end if
    if (error_tol == 0.0_wp) then
        tol = sqrt(epmach)
    else
        tol = max(abs(error_tol),2.0_wp**(5-nbits))/2.0_wp
    end if
    eps = tol
    hh(1) = (b-a)/4.0_wp
    aa(1) = a
    lr(1) = 1
    l = 1
    est = g(aa(l)+2.0_wp*hh(l),2.0_wp*hh(l))
    k = 8
    area = abs(est)
    ef = 0.5_wp
    mxl = 0

    !compute refined estimates, estimate the error, etc.
    main : do

        gl = g(aa(l)+hh(l),hh(l))
        gr(l) = g(aa(l)+3.0_wp*hh(l),hh(l))
        k = k + 16
        area = area + (abs(gl)+abs(gr(l))-abs(est))
        glr = gl + gr(l)
        ee = abs(est-glr)*ef
        ae = max(eps*area,tol*abs(glr))
        if (ee-ae > 0.0_wp) then
            !consider the left half of this level
            if (k > kmx) lmx = kml
            if (l >= lmx) then
                mxl = 1
            else
                l = l + 1
                eps = eps*0.5_wp
                ef = ef/sq2
                hh(l) = hh(l-1)*0.5_wp
                lr(l) = -1
                aa(l) = aa(l-1)
                est = gl
                cycle main
            end if
        end if

        err = err + (est-glr)
        if (lr(l) > 0) then
            !return one level
            ans = glr
            do
                if (l <= 1) exit main ! finished
                l = l - 1
                eps = eps*2.0_wp
                ef = ef*sq2
                if (lr(l) <= 0) then
                    vl(l) = vl(l+1) + ans
                    est = gr(l-1)
                    lr(l) = 1
                    aa(l) = aa(l) + 4.0_wp*hh(l)
                    cycle main
                end if
                ans = vl(l+1) + ans
            end do
        else
            !proceed to right half at this level
            vl(l) = glr
            est = gr(l-1)
            lr(l) = 1
            aa(l) = aa(l) + 4.0_wp*hh(l)
            cycle main
        end if

    end do main

    if ((mxl/=0) .and. (abs(err)>2.0_wp*tol*area)) ierr = 2 ! ans is probably insufficiently accurate

    contains

    !************************************************************************************
    !>
    !  This is the 8-point formula from the original SLATEC routine
    !  [DGAUS8](http://www.netlib.org/slatec/src/dgaus8.f).
    !
    !@note replaced coefficients with high-precision ones from:
    !      http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php

        function g(x, h)

        implicit none

        real(wp),intent(in) :: x
        real(wp),intent(in) :: h
        real(wp)            :: g

        !> abscissae:
        real(wp),parameter :: x1 = 0.18343464249564980493947614236018398066675781291297378231718847_wp
        real(wp),parameter :: x2 = 0.52553240991632898581773904918924634904196424312039285775085709_wp
        real(wp),parameter :: x3 = 0.79666647741362673959155393647583043683717173161596483207017029_wp
        real(wp),parameter :: x4 = 0.96028985649753623168356086856947299042823523430145203827163977_wp

        !> weights:
        real(wp),parameter :: w1 = 0.36268378337836198296515044927719561219414603989433054052482306_wp
        real(wp),parameter :: w2 = 0.31370664587788728733796220198660131326032899900273493769026394_wp
        real(wp),parameter :: w3 = 0.22238103445337447054435599442624088443013087005124956472590928_wp
        real(wp),parameter :: w4 = 0.10122853629037625915253135430996219011539409105168495705900369_wp

        g = h * ( w1*( f(x-x1*h) + f(x+x1*h) ) + &
                  w2*( f(x-x2*h) + f(x+x2*h) ) + &
                  w3*( f(x-x3*h) + f(x+x3*h) ) + &
                  w4*( f(x-x4*h) + f(x+x4*h) ) )

        end function g
    !************************************************************************************

    end subroutine dgauss8
!********************************************************************************

!********************************************************************************
!>
!  Numerically evaluate integral using adaptive Simpson rule.
!
!### See also
!  * W. Gander and W. Gautschi, "Adaptive Quadrature - Revisited",
!    BIT Vol. 40, No. 1, March 2000, pp. 84--101.

    recursive subroutine dsimpson (f, a, b, error_tol, ans, ierr)

    implicit none

    procedure(func)      :: f !! function subprogram defining the integrand function `f(x)`.
    real(wp),intent(in)  :: a !! lower bound of the integration
    real(wp),intent(in)  :: b !! upper bound of the integration
    real(wp),intent(in)  :: error_tol !! relative error tolerance
    real(wp),intent(out) :: ans !! computed value of integral
    integer,intent(out)  :: ierr !! status code:
                                 !!
                                 !! * 1 = success
                                 !! * 2 = requested accuracy may not be satisfied

    real(wp) :: bma,is,tol,fa,fm,fb
    real(wp),dimension(5) :: yy
    integer :: k !! number of calls to the recursive function

    real(wp),parameter :: eps = epsilon(1.0_wp)
    real(wp),dimension(5),parameter :: c = [.9501_wp, .2311_wp, .6068_wp, .4860_wp, .8913_wp]
    integer,parameter :: kmax = 10000 !! maximum number of calls to the recursive function (probably should be an input)

    k = 0
    ierr = 1
    bma = b-a
    tol = max(eps, error_tol)

    fa    = f(a)
    fm    = f((a+b)/2.0_wp)
    fb    = f(b)
    yy(1) = f(a+c(1)*bma )
    yy(2) = f(a+c(2)*bma )
    yy(3) = f(a+c(3)*bma )
    yy(4) = f(a+c(4)*bma )
    yy(5) = f(a+c(5)*bma )

    is = bma/8.0_wp * (fa+fm+fb+sum(yy))
    if (is==0.0_wp) is = bma
    is = is*tol/eps

    call adaptive_simpson_step(a,b,fa,fm,fb,is,ans)

    contains

    recursive subroutine adaptive_simpson_step (a,b,fa,fm,fb,is,ans)
        !!  Recursive function used by adaptive_simpson.
        !!  Tries to approximate the integral of f(x) from a to b
        !!  to an appropriate relative error.

        implicit none

        real(wp),intent(in)   :: a
        real(wp),intent(in)   :: b
        real(wp),intent(in)   :: fa
        real(wp),intent(in)   :: fm
        real(wp),intent(in)   :: fb
        real(wp),intent(in)   :: is
        real(wp),intent(out)  :: ans

        real(wp) :: m,h,fml,fmr,i1,i2,q1,q2

        k = k + 1
        if (k>kmax) then
            ierr = 2
            ans = 0.0_wp
            return
        end if
        m = (a + b)/2.0_wp
        h = (b - a)/4.0_wp
        fml = f(a + h)
        fmr = f(b - h)
        i1 = h/1.5_wp * (fa + 4.0_wp*fm + fb)
        i2 = h/3.0_wp * (fa + 4.0_wp*(fml + fmr) + 2.0_wp*fm + fb)
        i1 = (16.0_wp*i2 - i1)/15.0_wp

        if ( (is + (i1-i2) == is) .or. (m <= a) .or. (b <= m) ) then

            if ( ((m <= a) .or. (b<=m)) .and. (ierr==1) ) ierr = 2
            ans = i1

        else

            if (ierr==1) call adaptive_simpson_step (a,m,fa,fml,fm,is,q1)
            if (ierr==1) call adaptive_simpson_step (m,b,fm,fmr,fb,is,q2)

            if (ierr==1) then
                ans = q1 + q2
            else
                ans = i1
            end if

        end if

        end subroutine adaptive_simpson_step
    !**************************************************************

    end subroutine dsimpson
!********************************************************************************

!********************************************************************************
!>
!  Numerically evaluate integral using adaptive Lobatto rule
!
!### See also
!  * W. Gander and W. Gautschi, "Adaptive Quadrature - Revisited",
!    BIT Vol. 40, No. 1, March 2000, pp. 84--101.

    recursive subroutine dlobatto (f, a, b, error_tol, ans, ierr)

    procedure(func)      :: f !! function subprogram defining the integrand function `f(x)`.
    real(wp),intent(in)  :: a !! lower bound of the integration
    real(wp),intent(in)  :: b !! upper bound of the integration
    real(wp),intent(in)  :: error_tol !! relative error tolerance
    real(wp),intent(out) :: ans !! computed value of integral
    integer,intent(out)  :: ierr !! status code:
                                 !!
                                 !! * 1 = success
                                 !! * 2 = requested accuracy may not be satisfied

    real(wp) :: m,h,s,erri1,erri2,is,tol,fa,fb,i1,i2,r
    real(wp),dimension(13) :: x,y
    integer :: i
    integer :: k !! number of calls to the recursive function

    integer,parameter :: kmax = 10000 !! maximum number of calls to the recursive function (probably should be an input)
    real(wp),parameter :: eps  = epsilon(1.0_wp)
    real(wp),parameter :: alpha = sqrt(2.0_wp/3.0_wp)
    real(wp),parameter :: beta  = 1.0_wp/sqrt(5.0_wp)
    real(wp),parameter :: x1   = .94288241569547971905635175843185720232_wp
    real(wp),parameter :: x2   = .64185334234578130578123554132903188354_wp
    real(wp),parameter :: x3   = .23638319966214988028222377349205292599_wp
    real(wp),dimension(7) :: c = [ .015827191973480183087169986733305510591_wp, &
                                   .094273840218850045531282505077108171960_wp, &
                                   .15507198733658539625363597980210298680_wp, &
                                   .18882157396018245442000533937297167125_wp, &
                                   .19977340522685852679206802206648840246_wp, &
                                   .22492646533333952701601768799639508076_wp, &
                                   .24261107190140773379964095790325635233_wp ]
    k = 0
    ierr = 1
    tol = max(eps, error_tol)
    m = (a+b)/2.0_wp
    h = (b-a)/2.0_wp

    x = [a, m-x1*h, m-alpha*h, m-x2*h, m-beta*h, m-x3*h, m, m+x3*h, m+beta*h, m+x2*h, m+alpha*h, m+x1*h, b]
    do i=1,13
        y(i) = f(x(i))
    end do

    fa=y(1)
    fb=y(13)
    i2=(h/6.0_wp)*(y(1)+y(13)+5.0_wp*(y(5)+y(9)))
    i1=(h/1470.0_wp)*(77.0_wp*(y(1)+y(13))+432.0_wp*(y(3)+y(11))+625.0_wp*(y(5)+y(9))+672.0_wp*y(7))

    is = h*(c(1)*(y(1)+y(13)) + &
            c(2)*(y(2)+y(12)) + &
            c(3)*(y(3)+y(11)) + &
            c(4)*(y(4)+y(10)) + &
            c(5)*(y(5)+y(9))  + &
            c(6)*(y(6)+y(8))  + &
            c(7)*y(7))

    s = sign(1.0_wp,is)
    if (s==0.0_wp) s = 1.0_wp
    erri1 = abs(i1-is)
    erri2 = abs(i2-is)
    r = 1.0_wp
    if (erri2/=0.0_wp) r=erri1/erri2
    if (r>0.0_wp .and. r<1.0_wp) tol=tol/r
    is=s*abs(is)*tol/eps
    if (is==0.0_wp) is=b-a

    call adaptive_lobatto_step(a,b,fa,fb,is,ans)

    contains

    recursive subroutine adaptive_lobatto_step(a,b,fa,fb,is,ans)

        !!  Recursive function used by adaptive_lobatto.
        !!  Tries to approximate the integral of f(x) from a to b
        !!  to an appropriate relative error.

        implicit none

        real(wp),intent(in)   :: a
        real(wp),intent(in)   :: b
        real(wp),intent(in)   :: fa
        real(wp),intent(in)   :: fb
        real(wp),intent(in)   :: is
        real(wp),intent(out)  :: ans

        real(wp) :: h,m,mll,ml,mr,mrr,fmll,fml,fm,fmr,fmrr,i2,i1
        real(wp),dimension(6) :: q

        k = k + 1
        if (k>kmax) then
            ierr = 2
            ans = 0.0_wp
            return
        end if
        h   = (b-a)/2.0_wp
        m   = (a+b)/2.0_wp
        mll = m-alpha*h
        ml  = m-beta*h
        mr  = m+beta*h
        mrr = m+alpha*h

        fmll = f(mll)
        fml  = f(ml)
        fm   = f(m)
        fmr  = f(mr)
        fmrr = f(mrr)

        i2 = (h/6.0_wp)*(fa+fb+5.0_wp*(fml+fmr))
        i1 = (h/1470.0_wp)*(77.0_wp*(fa+fb)+432.0_wp*(fmll+fmrr)+625.0_wp*(fml+fmr)+672.0_wp*fm)

        if ( (is+(i1-i2)==is) .or. (mll<=a) .or. (b<=mrr) ) then

            if (((m <= a) .or. (b<=m)) .and. (ierr==1)) ierr = 2
            ans = i1

        else

            if (ierr==1) call adaptive_lobatto_step(a,mll,fa,fmll,    is,q(1))
            if (ierr==1) call adaptive_lobatto_step(mll,ml,fmll,fml,  is,q(2))
            if (ierr==1) call adaptive_lobatto_step(ml,m,fml,fm,      is,q(3))
            if (ierr==1) call adaptive_lobatto_step(m,mr,fm,fmr,      is,q(4))
            if (ierr==1) call adaptive_lobatto_step(mr,mrr,fmr,fmrr,  is,q(5))
            if (ierr==1) call adaptive_lobatto_step(mrr,b,fmrr,fb,    is,q(6))

            if (ierr==1) then
                ans = sum(q)
            else
                ans = i1
            end if

        end if

    end subroutine adaptive_lobatto_step
!**************************************************************

    end subroutine dlobatto
!********************************************************************************

!********************************************************************************
!>
!  XERROR processes a diagnostic message, in a manner
!  determined by the value of LEVEL and the current value
!  of the library error control flag, KONTRL.
!  (See subroutine XSETF for details.)
!
!     Examples
!```fortran
!  call xerror('smooth -- num was zero.',1,2)
!  call xerror('integ  -- less than full accuracy achieved.',2,1)
!  call xerror('rooter -- actual zero of f found before interval fully collapsed.',3,0)
!  call xerror('exp    -- underflows being set to zero.',1,-1)
!```
!
!### History
!  * Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!  * Latest SLATEC revision ---  19 MAR 1980
!  * Jacob Williams, Dec 2021 : rewrite simple version for new quadpack
!
!### References
!  * Jones R.E., Kahaner D.K., "Xerror, the slatec error-handling package",
!    sand82-0800, sandia laboratories, 1982.

    subroutine xerror(messg, nerr, level)
        use,intrinsic :: iso_fortran_env, only: error_unit
        implicit none

        character(len=*), intent(in) :: messg !! message to be processed
        integer, intent(in) :: nerr  !! the error number associated with this message.
                                     !! NERR must not be zero.
        integer, intent(in) :: level !! error category:
                                     !!  * =2 means this is an unconditionally fatal error.
                                     !!  * =1 means this is a recoverable error.  (I.e., it is
                                     !!    non-fatal if XSETF has been appropriately called.)
                                     !!  * =0 means this is a warning message only.
                                     !!  * =-1 means this is a warning message which is to be
                                     !!    printed at most once, regardless of how many
                                     !!    times this call is executed.

        write (error_unit, '(I5,1X,A)') nerr, messg
        if (level == 2) error stop

    end subroutine xerror
!********************************************************************************

#ifndef MOD_INCLUDE
!********************************************************************************
end module quadpack_generic
!********************************************************************************
#endif