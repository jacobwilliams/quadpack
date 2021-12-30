
!********************************************************************************
!>
!  Modernized QUADPACK
!
!### References
!  * Original version on NETLIB: http://www.netlib.org/quadpack/
!
!### Authors
!  * Piessens, Robert. Applied Mathematics and Programming Division, K. U. Leuven
!  * de Doncker, Elise. Applied Mathematics and Programming Division, K. U. Leuven
!  * Kahaner, D. K., (NBS)
!  * Jacob Williams, Dec 2021. Modernized the Fortran 77 code from Netlib.

module quadpack

    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp), dimension(5), parameter :: d1mach = [tiny(1.0_wp), &
                                                   huge(1.0_wp), &
                                                   real(radix(1.0_wp), &
                                                   kind(1.0_wp))**(-digits(1.0_wp)), &
                                                   epsilon(1.0_wp), &
                                                   log10(real(radix(1.0_wp), &
                                                   kind(1.0_wp)))] !! machine constants

    real(wp), parameter, private :: uflow = d1mach(1) !! the smallest positive magnitude.
    real(wp), parameter, private :: oflow = d1mach(2) !! the largest positive magnitude.
    real(wp), parameter, private :: epmach = d1mach(4) !! the largest relative spacing.
    real(wp), parameter, private :: pi = acos(-1.0_wp) !! pi

    integer, parameter, private :: limexp = 50 !! limexp is the maximum number of elements the epsilon
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

contains
!********************************************************************************

!********************************************************************************
!>
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

        procedure(func) :: f !! function subprogam defining the integrand function `f(x)`.
        real(wp), intent(in) :: a !! lower limit of integration
        real(wp), intent(out) :: Abserr !! estimate of the modulus of the absolute error,
                                        !! which should equal or exceed `abs(i-result)`
        real(wp), intent(in) :: b !! upper limit of integration
        real(wp), intent(in) :: Epsabs !! absolute accoracy requested
        real(wp), intent(in) :: Epsrel !! relative accuracy requested
                                       !! if epsabs<=0
                                       !! and epsrel<max(50*rel.mach.acc.,0.5e-28),
                                       !! the routine will end with ier = 6.
        real(wp), intent(out) :: Result !! approximation to the integral
        integer, intent(in) :: Lenw !! dimensioning parameter for `work`
                                    !! lenw must be at least limit*4.
                                    !! if lenw<limit*4, the routine will end with
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
                                    !!         determine the integration difficulaties.
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
                                    !!         (epsabs<=0 and
                                    !!          epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
                                    !!         or limit<1 or lenw<limit*4.
                                    !!         result, abserr, neval, last are set
                                    !!         to zero.
                                    !!         except when lenw is invalid, iwork(1),
                                    !!         work(limit*2+1) and work(limit*3+1) are
                                    !!         set to zero, work(1) is set to a and
                                    !!         work(limit+1) to b.
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
                                     !! produced in the subdiviosion process, which
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
        if (Ier /= 0) call xerror('abnormal return from dqag ', 26, Ier, lvl)

    end subroutine dqag
!********************************************************************************

!********************************************************************************
!>
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
        real(wp), intent(in) :: b !! uppwer limit of integration
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
                                            !! such that elist(iord(1)), ...,
                                            !! elist(iord(k)) form a decreasing sequence,
                                            !! with k = last if last<=(limit/2+2), and
                                            !! k = limit+1-last otherwise
        integer, intent(out) :: Last !! number of subintervals actually produced in the
                                     !! subdivision process

        real(wp) :: area1, a1, b1, defab1, error1 !! variable for the left subinterval
        real(wp) :: area2, a2, b2, defab2, error2 !! variable for the right subinterval
        real(wp) :: area !! sum of the integrals over the subintervals
        real(wp) :: area12 !! area1 + area2
        real(wp) :: erro12 !! error1 + error2
        real(wp) :: errsum !! sum of the errors over the subintervals
        real(wp) :: errmax !! elist(maxerr)
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
                                    !!         `(epsabs<=0 and
                                    !!         epsrel<max(50*rel.mach.acc.,0.5e-28))
                                    !!         or limit<1 or leniw<limit*4`.
                                    !!         `result`, `abserr`, `neval`, `last` are set to
                                    !!         zero. exept when limit or leniw is
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
        if (Ier /= 0) call xerror('abnormal return from dqagi', 26, Ier, lvl)

    end subroutine dqagi
!********************************************************************************

!********************************************************************************
!>
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
            Abserr == 0.0_wp) goto 400

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

        main: do Last = 2, Limit

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
            if (errsum <= errbnd) goto 300
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
                    if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle main
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
                        if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle main
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
        end do main

        ! set final result and error estimate.

        if (Abserr /= oflow) then
            if ((Ier + ierro) /= 0) then
                if (ierro == 3) Abserr = Abserr + correc
                if (Ier == 0) Ier = 3
                if (Result == 0.0_wp .or. area == 0.0_wp) then
                    if (Abserr > errsum) goto 300
                    if (area == 0.0_wp) goto 400
                elseif (Abserr/abs(Result) > errsum/abs(area)) then
                    goto 300
                end if
            end if

            ! test on divergence

            if (ksgn /= (-1) .or. max(abs(Result), abs(area)) > defabs*0.01_wp) then
                if (0.01_wp > (Result/area) .or. &
                    (Result/area) > 100.0_wp .or. &
                    errsum > abs(area)) Ier = 6
            end if
            goto 400
        end if

        ! compute global integral sum.

300     Result = sum(Rlist(1:Last))
        Abserr = errsum
400     Neval = 30*Last - 15
        if (Inf == 2) Neval = 2*Neval
        if (Ier > 2) Ier = Ier - 1

    end subroutine dqagie
!********************************************************************************

!********************************************************************************
!>
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
                                    !!   zero. exept when `leniw` or `lenw` or `npts2` is
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
        if (Ier /= 0) call xerror('abnormal return from dqagp', 26, Ier, lvl)

    end subroutine dqagp
!********************************************************************************

!********************************************************************************
!>
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
                    dres, ertest, resa, esabs, reseps, Result, &
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
        if (Ier /= 0 .or. Abserr <= errbnd) goto 400

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

        main: do Last = Npts2, Limit

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
            if (errsum <= errbnd) goto 300
            ! ***jump out of do-loop
            if (Ier /= 0) exit main
            if (.not. (noext)) then
                erlarg = erlarg - erlast
                if (levcur + 1 <= levmax) erlarg = erlarg + erro12
                if (.not. (extrap)) then
                    ! test whether the interval to be bisected next is the
                    ! smallest interval.
                    if (Level(maxerr) + 1 <= levmax) cycle main
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
                        if (Level(maxerr) + 1 <= levmax) cycle main
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
                        if (Abserr < ertest) exit main
                    end if
                    ! prepare bisection of the smallest interval.
                    if (numrl2 == 1) noext = .true.
                    if (Ier >= 5) exit main
                end if
                maxerr = Iord(1)
                errmax = Elist(maxerr)
                nrmax = 1
                extrap = .false.
                levmax = levmax + 1
                erlarg = errsum
            end if
        end do main

        ! set the final result.

        if (Abserr /= oflow) then
            if ((Ier + ierro) /= 0) then
                if (ierro == 3) Abserr = Abserr + correc
                if (Ier == 0) Ier = 3
                if (Result == 0.0_wp .or. area == 0.0_wp) then
                    if (Abserr > errsum) goto 300
                    if (area == 0.0_wp) goto 400
                elseif (Abserr/abs(Result) > errsum/abs(area)) then
                    goto 300
                end if
            end if

            ! test on divergence.

            if (ksgn /= (-1) .or. max(abs(Result), abs(area)) &
                > resabs*0.01_wp) then
                if (0.01_wp > (Result/area) .or. (Result/area) > 100.0_wp .or. &
                    errsum > abs(area)) Ier = 6
            end if
            goto 400
        end if

        ! compute global integral sum.

300     Result = sum(Rlist(1:Last))
        Abserr = errsum
400     if (Ier > 2) Ier = Ier - 1
        Result = Result*sign

    end subroutine dqagpe
!********************************************************************************

!********************************************************************************
!>
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
                                     !! produced in the subdivision process, detemines the
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
        if (Ier /= 0) call xerror('abnormal return from dqags', 26, Ier, lvl)

    end subroutine dqags
!********************************************************************************

!********************************************************************************
!>
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
        if (Ier /= 6) then

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
            else

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

                main: do Last = 2, Limit

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
                    if (errsum <= errbnd) goto 50
                    ! ***jump out of do-loop
                    if (Ier /= 0) exit main
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
                            if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle main
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
                                if (abs(Blist(maxerr) - Alist(maxerr)) > small) cycle main
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
                            if (Abserr <= ertest) exit main
                        end if

                        ! prepare bisection of the smallest interval.

                        if (numrl2 == 1) noext = .true.
                        if (Ier == 5) exit main
                        maxerr = Iord(1)
                        errmax = Elist(maxerr)
                        nrmax = 1
                        extrap = .false.
                        small = small*0.5_wp
                        erlarg = errsum
                    end if
                end do main

                ! set final result and error estimate.

                if (Abserr /= oflow) then
                    if (Ier + ierro /= 0) then
                        if (ierro == 3) Abserr = Abserr + correc
                        if (Ier == 0) Ier = 3
                        if (Result == 0.0_wp .or. area == 0.0_wp) then
                            if (Abserr > errsum) goto 50
                            if (area == 0.0_wp) then
                                if (Ier > 2) Ier = Ier - 1
                                Neval = 42*Last - 21
                                return
                            end if
                        elseif (Abserr/abs(Result) > errsum/abs(area)) then
                            goto 50
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
            end if

            ! compute global integral sum.

50          Result = sum(Rlist(1:Last))
            Abserr = errsum
            if (Ier > 2) Ier = Ier - 1
            Neval = 42*Last - 21
        end if

    end subroutine dqagse
!********************************************************************************

!********************************************************************************
!>
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
                                    !!   zero. exept when `lenw` or `limit` is invalid,
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
        if (Ier /= 0) call xerror('abnormal return from dqawc', 26, Ier, lvl)

    end subroutine dqawc
!********************************************************************************

!********************************************************************************
!>
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
                                   !!   convergence accelaration of the series
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
                                   !!   determined from the first lst elements of
                                   !!   vector iwork.  here `lst` is the number of
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
                                    !! lenw must be at least `leniw*2+maxp1*25`.
                                    !! if `lenw<(leniw*2+maxp1*25)`, the routine will
                                    !! end with ier = 6.
        integer :: Iwork(Leniw) !! vector of dimension at least `leniw`
                                !! on return, `iwork(k)` for `k = 1, 2, ..., lst`
                                !! contain the error flags on the cycles.
        real(wp) :: Work(Lenw) !! vector of dimension at least `lenw`
                               !! on return:
                               !!
                               !! * work(1), ..., work(lst) contain the integral
                               !!   approximations over the cycles,
                               !! * work(limlst+1), ..., work(limlst+lst) contain
                               !!   the error extimates over the cycles.
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
        if (Ier /= 0) call xerror('abnormal return from dqawf', 26, Ier, lvl)

    end subroutine dqawf
!********************************************************************************

!********************************************************************************
!>
!  the routine calculates an approximation result to a
!  given fourier integal
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
                                    !! * ier = 1 maximum number of  cycles  allowed
                                    !!   has been achieved., i.e. of subintervals
                                    !!   (a+(k-1)c,a+kc) where
                                    !!   c = (2*int(abs(omega))+1)*pi/abs(omega),
                                    !!   for k = 1, 2, ..., lst.
                                    !!   one can allow more cycles by increasing
                                    !!   the value of limlst (and taking the
                                    !!   according dimension adjustments into
                                    !!   account).
                                    !!   examine the array iwork which contains
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
                                    !!   array iwork which contains the error
                                    !!   flags on the cycles.
                                    !! * ier = 6 the input is invalid because
                                    !!   (integr/=1 and integr/=2) or
                                    !!   epsabs<=0 or limlst<3.
                                    !!   result, abserr, neval, lst are set
                                    !!   to zero.
                                    !! * ier = 7 bad integrand behaviour occurs within one
                                    !!   or more of the cycles. location and type
                                    !!   of the difficulty involved can be
                                    !!   determined from the vector ierlst. here
                                    !!   lst is the number of cycles actually
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
        real(wp), intent(out) :: Rslst(Limlst) !! vector of dimension at least limlst
                                               !! rslst(k) contains the integral contribution
                                               !! over the interval (a+(k-1)c,a+kc) where
                                               !! c = (2*int(abs(omega))+1)*pi/abs(omega),
                                               !! k = 1, 2, ..., lst.
                                               !! note that, if omega = 0, rslst(1) contains
                                               !! the value of the integral over (a,infinity).
        real(wp), intent(out) :: Erlst(Limlst) !! vector of dimension at least limlst
                                               !! erlst(k) contains the error estimate corresponding
                                               !! with rslst(k).
        integer, intent(out) :: Ierlst(Limlst) !! vector of dimension at least limlst
                                               !! ierlst(k) contains the error flag corresponding
                                               !! with rslst(k). for the meaning of the local error
                                               !! flags see description of output parameter ier.
        integer, intent(out) :: Lst !! number of subintervals needed for the integration
                                    !! if omega = 0 then lst is set to 1.
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
        if (Ier /= 6) then
            if (Omega /= 0.0_wp) then

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

                    if ((errsum + drl) <= Epsabs .and. Lst >= 6) goto 50
                    correc = max(correc, Erlst(Lst))
                    if (Ierlst(Lst) /= 0) eps = max(ep, correc*p1)
                    if (Ierlst(Lst) /= 0) Ier = 7
                    if (Ier == 7 .and. (errsum + drl) <= correc*10.0_wp .and. &
                        Lst > 5) goto 50
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
                    if (Abserr > errsum) goto 50
                    if (psum(numrl2) == 0.0_wp) return
                end if
                if (Abserr/abs(Result) <= (errsum + drl)/abs(psum(numrl2))) &
                    then
                    if (Ier >= 1 .and. Ier /= 7) Abserr = Abserr + drl
                    return
                end if
            else

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
50          Result = psum(numrl2)
            Abserr = errsum + drl
        end if

    end subroutine dqawfe
!********************************************************************************

!********************************************************************************
!>
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

        if (Ier == 6) lvl = 0
        if (Ier /= 0) call xerror('abnormal return from dqawo', 26, Ier, lvl)

    end subroutine dqawo
!********************************************************************************

!********************************************************************************
!>
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
        logical :: extall

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
        if (Ier /= 6) then

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
            else

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
                    ! ***jump out of do-loop
                    if (errsum <= errbnd) goto 50
                    if (Ier /= 0) goto 40
                    if (Last == 2 .and. extall) then
                        small = small*0.5_wp
                        numrl2 = numrl2 + 1
                        rlist2(numrl2) = area
                    else
                        if (noext) goto 20
                        if (extall) then
                            erlarg = erlarg - erlast
                            if (abs(b1 - a1) > small) erlarg = erlarg + erro12
                            if (extrap) goto 5
                        end if

                        ! test whether the interval to be bisected next is the
                        ! smallest interval.

                        width = abs(Blist(maxerr) - Alist(maxerr))
                        if (width > small) goto 20
                        if (extall) then
                            extrap = .true.
                            nrmax = 2
                        else

                            ! test whether we can start with the extrapolation procedure
                            ! (we do this if we integrate over the next interval with
                            ! use of a gauss-kronrod rule - see subroutine dqc25f).

                            small = small*0.5_wp
                            if (0.25_wp*width*domega > 2.0_wp) goto 20
                            extall = .true.
                            goto 10
                        end if
5                       if (ierro /= 3 .and. erlarg > ertest) then

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
                                    goto 20
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
                                if (Abserr <= ertest) goto 40
                            end if

                            ! prepare bisection of the smallest interval.

                            if (numrl2 == 1) noext = .true.
                            if (Ier == 5) goto 40
                        end if
                        maxerr = Iord(1)
                        errmax = Elist(maxerr)
                        nrmax = 1
                        extrap = .false.
                        small = small*0.5_wp
                        erlarg = errsum
                        goto 20
                    end if
10                  ertest = errbnd
                    erlarg = errsum
20              end do

                ! set the final result.

40              if (Abserr /= oflow .and. nres /= 0) then
                    if (Ier + ierro /= 0) then
                        if (ierro == 3) Abserr = Abserr + correc
                        if (Ier == 0) Ier = 3
                        if (Result == 0.0_wp .or. area == 0.0_wp) then
                            if (Abserr > errsum) goto 50
                            if (area == 0.0_wp) then
                                if (Ier > 2) Ier = Ier - 1
                                if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result
                                return
                            end if
                        elseif (Abserr/abs(Result) > errsum/abs(area)) then
                            goto 50
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
            end if

            ! compute global integral sum.

50          Result = 0.0_wp
            do k = 1, Last
                Result = Result + Rlist(k)
            end do
            Abserr = errsum
            if (Ier > 2) Ier = Ier - 1
            if (Integr == 2 .and. Omega < 0.0_wp) Result = -Result
        end if

    end subroutine dqawoe
!********************************************************************************

!********************************************************************************
!>
!  integration of functions having algebraico-logarithmic
!  end point singularities.
!
!  the routine calculates an approximation result to a given
!  definite integral i = integral of f*w over `(a,b)`,
!  (where w shows a singular behaviour at the end points
!  see parameter integr).
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
                                    !!   `(epsabs<=0 and epsrel<max(50*rel.mach.acc.,0.5d-28))`
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
        if (ier /= 0) call xerror('abnormal return from dqaws', 26, ier, lvl)

    end subroutine dqaws
!********************************************************************************

!********************************************************************************
!>
!  integration of functions having algebraico-logarithmic
!  end point singularities
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
!  integration rules for the computation of cauchy
!  principal value integrals
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
        real(wp), intent(out) :: Result !! approximation to the integral
                                        !! result is computed by using a generalized
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
                               !! for the function f, of degree 12
        real(wp) :: cheb24(25) !! chebyshev series expansion coefficients,
                               !! for the function f, of degree 24
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
                               !! of degree 24, for the function f, in the
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
                    v(3) = v(3) - 0.42d+02*par2*v(2)
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
                                       !! method (computed in subroutine dqmomo)
        real(wp), intent(in) :: Rj(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine dqmomo)
        real(wp), intent(in) :: Rg(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine dqmomo)
        real(wp), intent(in) :: Rh(25) !! modified chebyshev moments for the application
                                       !! of the generalized clenshaw-curtis
                                       !! method (computed in subroutine dqmomo)
        real(wp), intent(out) :: Result !! approximation to the integral
                                        !! result is computed by using a generalized
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
                               !! of degree 12, for the function f, in the
                               !! interval (bl,br)
        real(wp) :: cheb24(25) !! coefficients of the chebyshev series expansion
                               !! of degree 24, for the function f, in the
                               !! interval (bl,br)
        real(wp) :: fval(25) !! value of the function f at the points
                             !! (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
                             !! k = 0, ..., 24
        real(wp) :: res12 !! approximation to the integral obtained from cheb12
        real(wp) :: res24 !! approximation to the integral obtained from cheb24
        real(wp) :: hlgth !! half-length of the interval (bl,br)
        real(wp) :: centr !! mid point of the interval (bl,br)
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
                        if (epsinf > 0.1d-03) then
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
                    !goto 200
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
!  15-point gauss-kronrod rules
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
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(4), parameter :: wg = [ &
                                             0.129484966168869693270611432679082_wp, &
                                             0.279705391489276667901467771423780_wp, &
                                             0.381830050505118944950369775488975_wp, &
                                             0.417959183673469387755102040816327_wp] !! weights of the 7-point gauss rule

        real(wp), dimension(8), parameter :: xgk = [ &
                                             0.991455371120812639206854697526329_wp, &
                                             0.949107912342758524526189684047851_wp, &
                                             0.864864423359769072789712788640926_wp, &
                                             0.741531185599394439863864773280788_wp, &
                                             0.586087235467691130294144838258730_wp, &
                                             0.405845151377397166906606412076961_wp, &
                                             0.207784955007898467600689403773245_wp, &
                                             0.000000000000000000000000000000000_wp] !! abscissae of the 15-point kronrod rule:
                                                                                     !!
                                                                                     !! * xgk(2), xgk(4), ...  abscissae of the 7-point
                                                                                     !!   gauss rule
                                                                                     !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                     !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             0.022935322010529224963732008058970_wp, &
                                             0.063092092629978553290700663189204_wp, &
                                             0.104790010322250183839876322541518_wp, &
                                             0.140653259715525918745189590510238_wp, &
                                             0.169004726639267902826583426598550_wp, &
                                             0.190350578064785409913256402421014_wp, &
                                             0.204432940075298892414161999234649_wp, &
                                             0.209482141084727828012999174891714_wp] !! weights of the 15-point kronrod rule

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
!  15-point transformed gauss-kronrod rules
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

        procedure(func) :: f !! fuction subprogram defining the integrand function `f(x)`.
        real(wp), intent(in) :: Boun !! finite bound of original integration
                                     !! range (set to zero if inf = +2)
        real(wp), intent(in) :: a !! lower limit for integration over subrange of (0,1)
        real(wp), intent(in) :: b !! upper limit for integration over subrange of (0,1)
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
                                             0.0_wp, &
                                             0.129484966168869693270611432679082_wp, &
                                             0.0_wp, &
                                             0.279705391489276667901467771423780_wp, &
                                             0.0_wp, &
                                             0.381830050505118944950369775488975_wp, &
                                             0.0_wp, &
                                             0.417959183673469387755102040816327_wp] !! weights of the 7-point gauss rule, corresponding
                                                                                     !! to the abscissae `xgk(2), xgk(4), ...`.
                                                                                     !! `wg(1), wg(3), ...` are set to zero.

        real(wp), dimension(8), parameter :: xgk = [ &
                                             0.991455371120812639206854697526329_wp, &
                                             0.949107912342758524526189684047851_wp, &
                                             0.864864423359769072789712788640926_wp, &
                                             0.741531185599394439863864773280788_wp, &
                                             0.586087235467691130294144838258730_wp, &
                                             0.405845151377397166906606412076961_wp, &
                                             0.207784955007898467600689403773245_wp, &
                                             0.000000000000000000000000000000000_wp] !! abscissae of the 15-point kronrod rule:
                                                                                     !!
                                                                                     !! * xgk(2), xgk(4), ... abscissae of the 7-point
                                                                                     !!   gauss rule
                                                                                     !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                     !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             0.022935322010529224963732008058970_wp, &
                                             0.063092092629978553290700663189204_wp, &
                                             0.104790010322250183839876322541518_wp, &
                                             0.140653259715525918745189590510238_wp, &
                                             0.169004726639267902826583426598550_wp, &
                                             0.190350578064785409913256402421014_wp, &
                                             0.204432940075298892414161999234649_wp, &
                                             0.209482141084727828012999174891714_wp] !! weights of the 15-point kronrod rule

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
!  15-point gauss-kronrod rules
!
!  to compute i = integral of f*w over `(a,b)`, with error
!  estimate j = integral of abs(f*w) over `(a,b)`
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
                                             0.9914553711208126_wp, 0.9491079123427585_wp, &
                                             0.8648644233597691_wp, 0.7415311855993944_wp, &
                                             0.5860872354676911_wp, 0.4058451513773972_wp, &
                                             0.2077849550078985_wp, 0.0000000000000000_wp] !! abscissae of the 15-point gauss-kronrod rule:
                                                                                           !!
                                                                                           !! * xgk(2), xgk(4), ... abscissae of the 7-point
                                                                                           !!   gauss rule
                                                                                           !! * xgk(1), xgk(3), ... abscissae which are optimally
                                                                                           !!   added to the 7-point gauss rule

        real(wp), dimension(8), parameter :: wgk = [ &
                                             0.2293532201052922d-01, 0.6309209262997855d-01, &
                                             0.1047900103222502_wp, 0.1406532597155259_wp, &
                                             0.1690047266392679_wp, 0.1903505780647854_wp, &
                                             0.2044329400752989_wp, 0.2094821410847278_wp] !! weights of the 15-point gauss-kronrod rule

        real(wp), dimension(4), parameter :: wg = [ &
                                             0.1294849661688697_wp, 0.2797053914892767_wp, &
                                             0.3818300505051189_wp, 0.4179591836734694_wp] !! weights of the 7-point gauss rule

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
!  21-point gauss-kronrod rules
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
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(5), parameter :: wg = [ &
                                             0.066671344308688137593568809893332_wp, &
                                             0.149451349150580593145776339657697_wp, &
                                             0.219086362515982043995534934228163_wp, &
                                             0.269266719309996355091226921569469_wp, &
                                             0.295524224714752870173892994651338_wp] !! weights of the 10-point gauss rule

        real(wp), dimension(11), parameter :: xgk = [ &
                                              0.995657163025808080735527280689003_wp, &
                                              0.973906528517171720077964012084452_wp, &
                                              0.930157491355708226001207180059508_wp, &
                                              0.865063366688984510732096688423493_wp, &
                                              0.780817726586416897063717578345042_wp, &
                                              0.679409568299024406234327365114874_wp, &
                                              0.562757134668604683339000099272694_wp, &
                                              0.433395394129247190799265943165784_wp, &
                                              0.294392862701460198131126603103866_wp, &
                                              0.148874338981631210884826001129720_wp, &
                                              0.000000000000000000000000000000000_wp] !! abscissae of the 21-point kronrod rule:
                                                                                      !!
                                                                                      !! * xgk(2), xgk(4), ...  abscissae of the 10-point
                                                                                      !!   gauss rule
                                                                                      !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                      !!   added to the 10-point gauss rule

        real(wp), dimension(11), parameter :: wgk = [ &
                                              0.011694638867371874278064396062192_wp, &
                                              0.032558162307964727478818972459390_wp, &
                                              0.054755896574351996031381300244580_wp, &
                                              0.075039674810919952767043140916190_wp, &
                                              0.093125454583697605535065465083366_wp, &
                                              0.109387158802297641899210590325805_wp, &
                                              0.123491976262065851077958109831074_wp, &
                                              0.134709217311473325928054001771707_wp, &
                                              0.142775938577060080797094273138717_wp, &
                                              0.147739104901338491374841515972068_wp, &
                                              0.149445554002916905664936468389821_wp] !! weights of the 21-point kronrod rule

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
!  31-point gauss-kronrod rules
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
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(8), parameter :: wg = [ &
                                             0.030753241996117268354628393577204_wp, &
                                             0.070366047488108124709267416450667_wp, &
                                             0.107159220467171935011869546685869_wp, &
                                             0.139570677926154314447804794511028_wp, &
                                             0.166269205816993933553200860481209_wp, &
                                             0.186161000015562211026800561866423_wp, &
                                             0.198431485327111576456118326443839_wp, &
                                             0.202578241925561272880620199967519_wp] !! weights of the 15-point gauss rule

        real(wp), dimension(16), parameter :: xgk = [ &
                                              0.998002298693397060285172840152271_wp, &
                                              0.987992518020485428489565718586613_wp, &
                                              0.967739075679139134257347978784337_wp, &
                                              0.937273392400705904307758947710209_wp, &
                                              0.897264532344081900882509656454496_wp, &
                                              0.848206583410427216200648320774217_wp, &
                                              0.790418501442465932967649294817947_wp, &
                                              0.724417731360170047416186054613938_wp, &
                                              0.650996741297416970533735895313275_wp, &
                                              0.570972172608538847537226737253911_wp, &
                                              0.485081863640239680693655740232351_wp, &
                                              0.394151347077563369897207370981045_wp, &
                                              0.299180007153168812166780024266389_wp, &
                                              0.201194093997434522300628303394596_wp, &
                                              0.101142066918717499027074231447392_wp, &
                                              0.000000000000000000000000000000000_wp] !! abscissae of the 31-point kronrod rule:
                                                                                      !!
                                                                                      !! * xgk(2), xgk(4), ...  abscissae of the 15-point
                                                                                      !!   gauss rule
                                                                                      !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                      !!   added to the 15-point gauss rule

        real(wp), dimension(16), parameter :: wgk = [ &
                                              0.005377479872923348987792051430128_wp, &
                                              0.015007947329316122538374763075807_wp, &
                                              0.025460847326715320186874001019653_wp, &
                                              0.035346360791375846222037948478360_wp, &
                                              0.044589751324764876608227299373280_wp, &
                                              0.053481524690928087265343147239430_wp, &
                                              0.062009567800670640285139230960803_wp, &
                                              0.069854121318728258709520077099147_wp, &
                                              0.076849680757720378894432777482659_wp, &
                                              0.083080502823133021038289247286104_wp, &
                                              0.088564443056211770647275443693774_wp, &
                                              0.093126598170825321225486872747346_wp, &
                                              0.096642726983623678505179907627589_wp, &
                                              0.099173598721791959332393173484603_wp, &
                                              0.100769845523875595044946662617570_wp, &
                                              0.101330007014791549017374792767493_wp] !! weights of the 31-point kronrod rule

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
!  41-point gauss-kronrod rules
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
        real(wp), intent(out) :: Resasc !! approximation to the integal of abs(f-i/(b-a))
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
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(10), parameter :: wg = [ &
                                              0.017614007139152118311861962351853_wp, &
                                              0.040601429800386941331039952274932_wp, &
                                              0.062672048334109063569506535187042_wp, &
                                              0.083276741576704748724758143222046_wp, &
                                              0.101930119817240435036750135480350_wp, &
                                              0.118194531961518417312377377711382_wp, &
                                              0.131688638449176626898494499748163_wp, &
                                              0.142096109318382051329298325067165_wp, &
                                              0.149172986472603746787828737001969_wp, &
                                              0.152753387130725850698084331955098_wp] !! weights of the 20-point gauss rule

        real(wp), dimension(21), parameter :: xgk = [ &
                                              0.998859031588277663838315576545863_wp, &
                                              0.993128599185094924786122388471320_wp, &
                                              0.981507877450250259193342994720217_wp, &
                                              0.963971927277913791267666131197277_wp, &
                                              0.940822633831754753519982722212443_wp, &
                                              0.912234428251325905867752441203298_wp, &
                                              0.878276811252281976077442995113078_wp, &
                                              0.839116971822218823394529061701521_wp, &
                                              0.795041428837551198350638833272788_wp, &
                                              0.746331906460150792614305070355642_wp, &
                                              0.693237656334751384805490711845932_wp, &
                                              0.636053680726515025452836696226286_wp, &
                                              0.575140446819710315342946036586425_wp, &
                                              0.510867001950827098004364050955251_wp, &
                                              0.443593175238725103199992213492640_wp, &
                                              0.373706088715419560672548177024927_wp, &
                                              0.301627868114913004320555356858592_wp, &
                                              0.227785851141645078080496195368575_wp, &
                                              0.152605465240922675505220241022678_wp, &
                                              0.076526521133497333754640409398838_wp, &
                                              0.000000000000000000000000000000000_wp] !! abscissae of the 41-point gauss-kronrod rule:
                                                                                      !!
                                                                                      !! * xgk(2), xgk(4), ...  abscissae of the 20-point
                                                                                      !!   gauss rule
                                                                                      !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                      !!   added to the 20-point gauss rule

        real(wp), dimension(21), parameter :: wgk = [ &
                                              0.003073583718520531501218293246031_wp, &
                                              0.008600269855642942198661787950102_wp, &
                                              0.014626169256971252983787960308868_wp, &
                                              0.020388373461266523598010231432755_wp, &
                                              0.025882133604951158834505067096153_wp, &
                                              0.031287306777032798958543119323801_wp, &
                                              0.036600169758200798030557240707211_wp, &
                                              0.041668873327973686263788305936895_wp, &
                                              0.046434821867497674720231880926108_wp, &
                                              0.050944573923728691932707670050345_wp, &
                                              0.055195105348285994744832372419777_wp, &
                                              0.059111400880639572374967220648594_wp, &
                                              0.062653237554781168025870122174255_wp, &
                                              0.065834597133618422111563556969398_wp, &
                                              0.068648672928521619345623411885368_wp, &
                                              0.071054423553444068305790361723210_wp, &
                                              0.073030690332786667495189417658913_wp, &
                                              0.074582875400499188986581418362488_wp, &
                                              0.075704497684556674659542775376617_wp, &
                                              0.076377867672080736705502835038061_wp, &
                                              0.076600711917999656445049901530102_wp] !! weights of the 41-point gauss-kronrod rule

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
        if (Resasc /= 0.0_wp .and. Abserr /= 0._wp) &
            Abserr = Resasc*min(1.0_wp, (200.0_wp*Abserr/Resasc)**1.5_wp)
        if (Resabs > uflow/(50.0_wp*epmach)) &
            Abserr = max((epmach*50.0_wp)*Resabs, Abserr)

    end subroutine dqk41
!********************************************************************************

!********************************************************************************
!>
!  51-point gauss-kronrod rules
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
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(13), parameter :: wg = [ &
                                              0.011393798501026287947902964113235_wp, &
                                              0.026354986615032137261901815295299_wp, &
                                              0.040939156701306312655623487711646_wp, &
                                              0.054904695975835191925936891540473_wp, &
                                              0.068038333812356917207187185656708_wp, &
                                              0.080140700335001018013234959669111_wp, &
                                              0.091028261982963649811497220702892_wp, &
                                              0.100535949067050644202206890392686_wp, &
                                              0.108519624474263653116093957050117_wp, &
                                              0.114858259145711648339325545869556_wp, &
                                              0.119455763535784772228178126512901_wp, &
                                              0.122242442990310041688959518945852_wp, &
                                              0.123176053726715451203902873079050_wp] !! weights of the 25-point gauss rule

        real(wp), dimension(26), parameter :: xgk = [ &
                                              0.999262104992609834193457486540341_wp, &
                                              0.995556969790498097908784946893902_wp, &
                                              0.988035794534077247637331014577406_wp, &
                                              0.976663921459517511498315386479594_wp, &
                                              0.961614986425842512418130033660167_wp, &
                                              0.942974571228974339414011169658471_wp, &
                                              0.920747115281701561746346084546331_wp, &
                                              0.894991997878275368851042006782805_wp, &
                                              0.865847065293275595448996969588340_wp, &
                                              0.833442628760834001421021108693570_wp, &
                                              0.797873797998500059410410904994307_wp, &
                                              0.759259263037357630577282865204361_wp, &
                                              0.717766406813084388186654079773298_wp, &
                                              0.673566368473468364485120633247622_wp, &
                                              0.626810099010317412788122681624518_wp, &
                                              0.577662930241222967723689841612654_wp, &
                                              0.526325284334719182599623778158010_wp, &
                                              0.473002731445714960522182115009192_wp, &
                                              0.417885382193037748851814394594572_wp, &
                                              0.361172305809387837735821730127641_wp, &
                                              0.303089538931107830167478909980339_wp, &
                                              0.243866883720988432045190362797452_wp, &
                                              0.183718939421048892015969888759528_wp, &
                                              0.122864692610710396387359818808037_wp, &
                                              0.061544483005685078886546392366797_wp, &
                                              0.000000000000000000000000000000000_wp] !! abscissae of the 51-point kronrod rule"
                                                                                      !!
                                                                                      !! * xgk(2), xgk(4), ...  abscissae of the 25-point
                                                                                      !!   gauss rule
                                                                                      !! * xgk(1), xgk(3), ...  abscissae which are optimally
                                                                                      !!   added to the 25-point gauss rule

        real(wp), dimension(26), parameter :: wgk = [ &
                                              0.001987383892330315926507851882843_wp, &
                                              0.005561932135356713758040236901066_wp, &
                                              0.009473973386174151607207710523655_wp, &
                                              0.013236229195571674813656405846976_wp, &
                                              0.016847817709128298231516667536336_wp, &
                                              0.020435371145882835456568292235939_wp, &
                                              0.024009945606953216220092489164881_wp, &
                                              0.027475317587851737802948455517811_wp, &
                                              0.030792300167387488891109020215229_wp, &
                                              0.034002130274329337836748795229551_wp, &
                                              0.037116271483415543560330625367620_wp, &
                                              0.040083825504032382074839284467076_wp, &
                                              0.042872845020170049476895792439495_wp, &
                                              0.045502913049921788909870584752660_wp, &
                                              0.047982537138836713906392255756915_wp, &
                                              0.050277679080715671963325259433440_wp, &
                                              0.052362885806407475864366712137873_wp, &
                                              0.054251129888545490144543370459876_wp, &
                                              0.055950811220412317308240686382747_wp, &
                                              0.057437116361567832853582693939506_wp, &
                                              0.058689680022394207961974175856788_wp, &
                                              0.059720340324174059979099291932562_wp, &
                                              0.060539455376045862945360267517565_wp, &
                                              0.061128509717053048305859030416293_wp, &
                                              0.061471189871425316661544131965264_wp, &
                                              0.061580818067832935078759824240066_wp] !! weights of the 51-point kronrod rule.
                                                                                      !! note: `wgk(26)` was calculated from
                                                                                      !! the values of `wgk(1..25)`

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
!  61-point gauss-kronrod rules
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
        real(wp) :: dabsc !! abscissa
        real(wp) :: fval1 !! function value
        real(wp) :: fval2 !! function value
        real(wp) :: resg !! result of the 30-point gauss rule
        real(wp) :: resk !! result of the 61-point kronrod rule
        real(wp) :: reskh !! approximation to the mean value of `f` over `(a,b)`, i.e. to `i/(b-a)`

        ! the abscissae and weights are given for the
        ! interval (-1,1). because of symmetry only the positive
        ! abscissae and their corresponding weights are given.
        !
        ! gauss quadrature weights and kronrod quadrature abscissae and weights
        ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
        ! bell labs, nov. 1981.

        real(wp), dimension(15), parameter :: wg = [ &
                                              0.007968192496166605615465883474674_wp, &
                                              0.018466468311090959142302131912047_wp, &
                                              0.028784707883323369349719179611292_wp, &
                                              0.038799192569627049596801936446348_wp, &
                                              0.048402672830594052902938140422808_wp, &
                                              0.057493156217619066481721689402056_wp, &
                                              0.065974229882180495128128515115962_wp, &
                                              0.073755974737705206268243850022191_wp, &
                                              0.080755895229420215354694938460530_wp, &
                                              0.086899787201082979802387530715126_wp, &
                                              0.092122522237786128717632707087619_wp, &
                                              0.096368737174644259639468626351810_wp, &
                                              0.099593420586795267062780282103569_wp, &
                                              0.101762389748405504596428952168554_wp, &
                                              0.102852652893558840341285636705415_wp] !! weigths of the 30-point gauss rule

        real(wp), dimension(31), parameter :: xgk = [ &
                                              0.999484410050490637571325895705811_wp, &
                                              0.996893484074649540271630050918695_wp, &
                                              0.991630996870404594858628366109486_wp, &
                                              0.983668123279747209970032581605663_wp, &
                                              0.973116322501126268374693868423707_wp, &
                                              0.960021864968307512216871025581798_wp, &
                                              0.944374444748559979415831324037439_wp, &
                                              0.926200047429274325879324277080474_wp, &
                                              0.905573307699907798546522558925958_wp, &
                                              0.882560535792052681543116462530226_wp, &
                                              0.857205233546061098958658510658944_wp, &
                                              0.829565762382768397442898119732502_wp, &
                                              0.799727835821839083013668942322683_wp, &
                                              0.767777432104826194917977340974503_wp, &
                                              0.733790062453226804726171131369528_wp, &
                                              0.697850494793315796932292388026640_wp, &
                                              0.660061064126626961370053668149271_wp, &
                                              0.620526182989242861140477556431189_wp, &
                                              0.579345235826361691756024932172540_wp, &
                                              0.536624148142019899264169793311073_wp, &
                                              0.492480467861778574993693061207709_wp, &
                                              0.447033769538089176780609900322854_wp, &
                                              0.400401254830394392535476211542661_wp, &
                                              0.352704725530878113471037207089374_wp, &
                                              0.304073202273625077372677107199257_wp, &
                                              0.254636926167889846439805129817805_wp, &
                                              0.204525116682309891438957671002025_wp, &
                                              0.153869913608583546963794672743256_wp, &
                                              0.102806937966737030147096751318001_wp, &
                                              0.051471842555317695833025213166723_wp, &
                                              0.000000000000000000000000000000000_wp] !! abscissae of the 61-point kronrod rule:
                                                                                      !!
                                                                                      !! * `xgk(2), xgk(4)`  ... abscissae of the 30-point
                                                                                      !!   gauss rule
                                                                                      !! * `xgk(1), xgk(3)`  ... optimally added abscissae
                                                                                      !!   to the 30-point gauss rule

        real(wp), dimension(31), parameter :: wgk = [ &
                                              0.001389013698677007624551591226760_wp, &
                                              0.003890461127099884051267201844516_wp, &
                                              0.006630703915931292173319826369750_wp, &
                                              0.009273279659517763428441146892024_wp, &
                                              0.011823015253496341742232898853251_wp, &
                                              0.014369729507045804812451432443580_wp, &
                                              0.016920889189053272627572289420322_wp, &
                                              0.019414141193942381173408951050128_wp, &
                                              0.021828035821609192297167485738339_wp, &
                                              0.024191162078080601365686370725232_wp, &
                                              0.026509954882333101610601709335075_wp, &
                                              0.028754048765041292843978785354334_wp, &
                                              0.030907257562387762472884252943092_wp, &
                                              0.032981447057483726031814191016854_wp, &
                                              0.034979338028060024137499670731468_wp, &
                                              0.036882364651821229223911065617136_wp, &
                                              0.038678945624727592950348651532281_wp, &
                                              0.040374538951535959111995279752468_wp, &
                                              0.041969810215164246147147541285970_wp, &
                                              0.043452539701356069316831728117073_wp, &
                                              0.044814800133162663192355551616723_wp, &
                                              0.046059238271006988116271735559374_wp, &
                                              0.047185546569299153945261478181099_wp, &
                                              0.048185861757087129140779492298305_wp, &
                                              0.049055434555029778887528165367238_wp, &
                                              0.049795683427074206357811569379942_wp, &
                                              0.050405921402782346840893085653585_wp, &
                                              0.050881795898749606492297473049805_wp, &
                                              0.051221547849258772170656282604944_wp, &
                                              0.051426128537459025933862879215781_wp, &
                                              0.051494729429451567558340433647099_wp] !! weights of the 61-point kronrod rule

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
            dabsc = hlgth*xgk(jtw)
            fval1 = f(centr - dabsc)
            fval2 = f(centr + dabsc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum = fval1 + fval2
            resg = resg + wg(j)*fsum
            resk = resk + wgk(jtw)*fsum
            Resabs = Resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        do j = 1, 15
            jtwm1 = j*2 - 1
            dabsc = hlgth*xgk(jtwm1)
            fval1 = f(centr - dabsc)
            fval2 = f(centr + dabsc)
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
!  this routine computes modified chebsyshev moments. the k-th
!  modified chebyshev moment is defined as the integral over
!  (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
!  polynomial of degree k.
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
                if (Integr == 2) goto 100
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

100     do i = 2, 25, 2
            Rj(i) = -Rj(i)
        end do

    end subroutine dqmomo
!********************************************************************************

!********************************************************************************
!>
!  non-adaptive integration
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
        !
        ! gauss-kronrod-patterson quadrature coefficients for use in
        ! quadpack routine qng.  these coefficients were calculated with
        ! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.

        real(wp), dimension(5), parameter :: x1 = [ &
                                             0.973906528517171720077964012084452_wp, &
                                             0.865063366688984510732096688423493_wp, &
                                             0.679409568299024406234327365114874_wp, &
                                             0.433395394129247190799265943165784_wp, &
                                             0.148874338981631210884826001129720_wp] !! abscissae common to the 10-, 21-, 43- and 87-point rule

        real(wp), dimension(5), parameter :: w10 = [ &
                                             0.066671344308688137593568809893332_wp, &
                                             0.149451349150580593145776339657697_wp, &
                                             0.219086362515982043995534934228163_wp, &
                                             0.269266719309996355091226921569469_wp, &
                                             0.295524224714752870173892994651338_wp] !! weights of the 10-point formula

        real(wp), dimension(5), parameter :: x2 = [ &
                                             0.995657163025808080735527280689003_wp, &
                                             0.930157491355708226001207180059508_wp, &
                                             0.780817726586416897063717578345042_wp, &
                                             0.562757134668604683339000099272694_wp, &
                                             0.294392862701460198131126603103866_wp] !! abscissae common to the 21-, 43- and 87-point rule

        real(wp), dimension(5), parameter :: w21a = [ &
                                             0.032558162307964727478818972459390_wp, &
                                             0.075039674810919952767043140916190_wp, &
                                             0.109387158802297641899210590325805_wp, &
                                             0.134709217311473325928054001771707_wp, &
                                             0.147739104901338491374841515972068_wp] !! weights of the 21-point formula for abscissae x1

        real(wp), dimension(6), parameter :: w21b = [ &
                                             0.011694638867371874278064396062192_wp, &
                                             0.054755896574351996031381300244580_wp, &
                                             0.093125454583697605535065465083366_wp, &
                                             0.123491976262065851077958109831074_wp, &
                                             0.142775938577060080797094273138717_wp, &
                                             0.149445554002916905664936468389821_wp] !! weights of the 21-point formula for abscissae x2

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
        call xerror('abnormal return from dqng ', 26, Ier, 0)

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
                    if (errmax >= Elist(isucc)) goto 100
                    Iord(i - 1) = isucc
                end do
            end if
            Iord(jbnd) = Maxerr
            Iord(jupbn) = Last
        else
            Iord(1) = 1
            Iord(2) = 2
        end if
        goto 300

        ! insert errmin by traversing the list bottom-up.

100     Iord(i - 1) = Maxerr
        k = jbnd
        do j = i, jbnd
            isucc = Iord(k)
            ! ***jump out of do-loop
            if (errmin < Elist(isucc)) goto 200
            Iord(k + 1) = isucc
            k = k - 1
        end do
        Iord(i) = Last
        goto 300

200     Iord(k + 1) = Last

        ! set maxerr and ermax.
300     Maxerr = Iord(Nrmax)
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
!  XERROR processes a diagnostic message, in a manner
!  determined by the value of LEVEL and the current value
!  of the library error control flag, KONTRL.
!  (See subroutine XSETF for details.)
!
!     Examples
!```fortran
!  call xerror('smooth -- num was zero.',23,1,2)
!  call xerror('integ  -- less than full accuracy achieved.',43,2,1)
!  call xerror('rooter -- actual zero of f found before interval fully collapsed.',65,3,0)
!  call xerror('exp    -- underflows being set to zero.',39,1,-1)
!```
!
!### History
!  * Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!  * Latest SLATEC revision ---  19 MAR 1980
!  * Jacob Williams, Dec 2021 : rewrite simple version for new quadpack
!
!### References
!  * Jones R.E., Kahaner D.K., "Xerror, the slatec error-
!    handling package", sand82-0800, sandia laboratories,
!    1982.

    subroutine xerror(messg, nmessg, nerr, level)
        implicit none

        character(len=*), intent(in) :: messg !! message to be processed
        integer, intent(in) :: nmessg !! the actual number of characters in MESSG
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

        !call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)

        write (*, *) nerr, messg(1:nmessg)
        if (level == 2) error stop

    end subroutine xerror
!********************************************************************************

!********************************************************************************
end module quadpack
!********************************************************************************