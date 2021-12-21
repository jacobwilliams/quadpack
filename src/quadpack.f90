module quadpack

    use iso_fortran_env, only: wp => real64

    implicit none

    abstract interface
        real(wp) function func(x)
        !! interface for user-supplied function.
        import :: wp
        implicit none
        real(wp),intent(in) :: x
        end function func
    end interface

    contains
!********************************************************************************

!********************************************************************************

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             integrand examinator, globally adaptive,
!             gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)le.max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!            f      - real(wp)
!                     function subprogam defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            epsabs - real(wp)
!                     absolute accoracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                       7 - 15 points if key<2,
!                      10 - 21 points if key = 2,
!                      15 - 31 points if key = 3,
!                      20 - 41 points if key = 4,
!                      25 - 51 points if key = 5,
!                      30 - 61 points if key>5.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                      error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yield no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulaties.
!                             if the position of a local difficulty can
!                             be determined (i.e.singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or limit<1 or lenw<limit*4.
!                             result, abserr, neval, last are set
!                             to zero.
!                             except when lenw is invalid, iwork(1),
!                             work(limit*2+1) and work(limit*3+1) are
!                             set to zero, work(1) is set to a and
!                             work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit>=1.
!                    if limit<1, the routine will end with ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end with
!                    ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdiviosion process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first k
!                    elements of which contain pointers to the error
!                    estimates over the subintervals, such that
!                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
!                    form a decreasing sequence with k = last if
!                    last<=(limit/2+2), and k = limit+1-last otherwise
!
!            work  - real(wp)
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left end
!                    points of the subintervals in the partition of
!                     (a,b),
!                    work(limit+1), ..., work(limit+last) contain the
!                     right end points,
!                    work(limit*2+1), ..., work(limit*2+last) contain
!                     the integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3+last) contain
!                     the error estimates.

      SUBROUTINE DQAG(F,A,B,Epsabs,Epsrel,Key,Result,Abserr,Neval,Ier, &
                      Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

      real(wp) A , Abserr , B , Epsabs , Epsrel , Result , &
                       Work
      INTEGER Ier , Iwork , Key , Last , Lenw , Limit , lvl , l1 , l2 , &
              l3 , Neval

      DIMENSION Iwork(Limit) , Work(Lenw)

      procedure(func) :: f

    ! check validity of lenw.

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN

        ! prepare call for dqage.

         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2

         CALL DQAGE(F,A,B,Epsabs,Epsrel,Key,Limit,Result,Abserr,Neval, &
                    Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)

        ! call error handler if necessary.

         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqag ',26,Ier,lvl)

      END SUBROUTINE DQAG
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGE(F,A,B,Epsabs,Epsrel,Key,Limit,Result,Abserr, &
                       Neval,Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             integrand examinator, globally adaptive,
!             gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral   i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-reslt)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                          7 - 15 points if key<2,
!                         10 - 21 points if key = 2,
!                         15 - 31 points if key = 3,
!                         20 - 41 points if key = 4,
!                         25 - 51 points if key = 5,
!                         30 - 61 points if key>5.
!
!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b), limit>=1.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value
!                             of limit.
!                             however, if this yields no improvement it
!                             is rather advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined(e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                             result, abserr, neval, last, rlist(1) ,
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!            alist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            blist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            rlist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the
!                      integral approximations on the subintervals
!
!            elist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the moduli of the
!                      absolute error estimates on the subintervals
!
!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the
!                      error estimates over the subintervals,
!                      such that elist(iord(1)), ...,
!                      elist(iord(k)) form a decreasing sequence,
!                      with k = last if last<=(limit/2+2), and
!                      k = limit+1-last otherwise
!
!            last    - integer
!                      number of subintervals actually produced in the
!                      subdivision process
!

      real(wp) A , Abserr , Alist , area , area1 , area12 , &
                       area2 , a1 , a2 , B , Blist , b1 , b2 , &
                       defabs , defab1 , defab2 , &
                       Elist , epmach , Epsabs , Epsrel , errbnd , &
                       errmax , error1 , error2 , erro12 , errsum , &
                       resabs , Result , Rlist , uflow
      INTEGER Ier , Iord , iroff1 , iroff2 , k , Key , keyf , Last , &
              Limit , maxerr , Neval , nrmax
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
                , Rlist(Limit)
!
      procedure(func) :: f
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                      (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
!
!           machine dependent constants
!           ---------------------------
!
!           epmach  is the largest relative spacing.
!           uflow  is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      Iord(1) = 0
      IF ( Epsabs<=0.0_wp .AND. Epsrel<max(50.0_wp*epmach,0.5e-28_wp) ) &
           Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         keyf = Key
         IF ( Key<=0 ) keyf = 1
         IF ( Key>=7 ) keyf = 6
         Neval = 0
         IF ( keyf==1 ) CALL DQK15(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==2 ) CALL DQK21(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==3 ) CALL DQK31(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==4 ) CALL DQK41(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==5 ) CALL DQK51(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==6 ) CALL DQK61(F,A,B,Result,Abserr,defabs,resabs)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
!
!           test on accuracy.
!
         errbnd = max(Epsabs,Epsrel*abs(Result))
         IF ( Abserr<=50.0_wp*epmach*defabs .AND. Abserr>errbnd ) &
              Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( .NOT.(Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) &
              .OR. Abserr==0.0_wp) ) THEN
!
!           initialization
!           --------------
!
!
            errmax = Abserr
            maxerr = 1
            area = Result
            errsum = Abserr
            nrmax = 1
            iroff1 = 0
            iroff2 = 0
!
!           main do-loop
!           ------------
!
            DO Last = 2 , Limit
!
!           bisect the subinterval with the largest error estimate.
!
               a1 = Alist(maxerr)
               b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               IF ( keyf==1 ) CALL DQK15(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==2 ) CALL DQK21(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==3 ) CALL DQK31(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==4 ) CALL DQK41(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==5 ) CALL DQK51(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==6 ) CALL DQK61(F,a1,b1,area1,error1,resabs, &
                    defab1)
               IF ( keyf==1 ) CALL DQK15(F,a2,b2,area2,error2,resabs, &
                    defab2)
               IF ( keyf==2 ) CALL DQK21(F,a2,b2,area2,error2,resabs, &
                    defab2)
               IF ( keyf==3 ) CALL DQK31(F,a2,b2,area2,error2,resabs, &
                    defab2)
               IF ( keyf==4 ) CALL DQK41(F,a2,b2,area2,error2,resabs, &
                    defab2)
               IF ( keyf==5 ) CALL DQK51(F,a2,b2,area2,error2,resabs, &
                    defab2)
               IF ( keyf==6 ) CALL DQK61(F,a2,b2,area2,error2,resabs, &
                    defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               Neval = Neval + 1
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
                  IF ( abs(Rlist(maxerr)-area12)<=0.1e-4_wp*abs(area12) &
                       .AND. erro12>=0.99_wp*errmax ) iroff1 = iroff1 +&
                       1
                  IF ( Last>10 .AND. erro12>errmax ) iroff2 = iroff2 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = max(Epsabs,Epsrel*abs(area))
               IF ( errsum>errbnd ) THEN
!
!           test for roundoff error and eventually set error flag.
!
                  IF ( iroff1>=6 .OR. iroff2>=20 ) Ier = 2
!
!           set error flag in the case that the number of subintervals
!           equals limit.
!
                  IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range.
!
                  IF ( max(abs(a1),abs(b2)) &
                       <=(1.0_wp+100.0_wp*epmach) &
                       *(abs(a2)+1000.0_wp*uflow) ) Ier = 3
               ENDIF
!
!           append the newly-created intervals to the list.
!
               IF ( error2>error1 ) THEN
                  Alist(maxerr) = a2
                  Alist(Last) = a1
                  Blist(Last) = b1
                  Rlist(maxerr) = area2
                  Rlist(Last) = area1
                  Elist(maxerr) = error2
                  Elist(Last) = error1
               ELSE
                  Alist(Last) = a2
                  Blist(maxerr) = b1
                  Blist(Last) = b2
                  Elist(maxerr) = error1
                  Elist(Last) = error2
               ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with the largest error estimate (to be bisected next).
!
               CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0_wp
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
         IF ( keyf/=1 ) Neval = (10*keyf+1)*(2*Neval+1)
         IF ( keyf==1 ) Neval = 30*Neval + 15
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGI(F,Bound,Inf,Epsabs,Epsrel,Result,Abserr,Neval, &
                       Ier,Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, infinite intervals,
!             general-purpose, transformation, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. -k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!            or i = integral of f over (-infinity,bound)
!            or i = integral of f over (-infinity,+infinity)
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        integration over infinite intervals
!        standard fortran subroutine
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            bound  - real(wp)
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - integer
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier>0 abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is assumed that the requested tolerance
!                             cannot be achieved, and that the returned
!                             result is the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                              or limit<1 or leniw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero. exept when limit or leniw is
!                             invalid, iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit>=1.
!                    if limit<1, the routine will end with ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end
!                    with ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first
!                    k elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that work(limit*3+iwork(1)),... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2), and
!                    k = limit+1-last otherwise
!
!            work  - real(wp)
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                     end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ...,work(limit*2+last) contain the
!                     integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3)
!                     contain the error estimates.

!
      real(wp) Abserr , Bound , Epsabs , Epsrel , Result , &
                       Work
      INTEGER Ier , Inf , Iwork , Last , Lenw , Limit , lvl , l1 , l2 , &
              l3 , Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of limit and lenw.
!

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for dqagie.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL DQAGIE(F,Bound,Inf,Epsabs,Epsrel,Limit,Result,Abserr, &
                     Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,&
                     Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqagi',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGIE(F,Bound,Inf,Epsabs,Epsrel,Limit,Result,Abserr, &
                        Neval,Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, infinite intervals,
!             general-purpose, transformation, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math & progr. div - k.u.leuven
!           de doncker,elise,appl. math & progr. div - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!            or i = integral of f over (-infinity,bound)
!            or i = integral of f over (-infinity,+infinity),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i))
!***description
!
! integration over infinite intervals
! standard fortran subroutine
!
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            bound  - real(wp)
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - real(wp)
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=1
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier>0 abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however,if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is assumed that the requested tolerance
!                             cannot be achieved, and that the returned
!                             result is the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                             result, abserr, neval, last, rlist(1),
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to 0
!                             and 1 respectively.
!
!            alist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            blist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            rlist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real(wp)
!                     vector of dimension at least limit,  the first
!                     last elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            iord   - integer
!                     vector of dimension limit, the first k
!                     elements of which are pointers to the
!                     error estimates over the subintervals,
!                     such that elist(iord(1)), ..., elist(iord(k))
!                     form a decreasing sequence, with k = last
!                     if last<=(limit/2+2), and k = limit+1-last
!                     otherwise
!
!            last   - integer
!                     number of subintervals actually produced
!                     in the subdivision process
!

      real(wp) abseps , Abserr , Alist , area , area1 , area12 ,&
                       area2 , a1 , a2 , Blist , boun , Bound , b1 , &
                       b2 , correc , defabs , defab1 , defab2 , &
                       dres , Elist , epmach , Epsabs ,&
                       Epsrel , erlarg , erlast , errbnd , errmax , &
                       error1 , error2 , erro12 , errsum , ertest , &
                       oflow , resabs , reseps , Result , res3la , &
                       Rlist , rlist2 , small , uflow
      INTEGER id , Ier , ierro , Inf , Iord , iroff1 , iroff2 , iroff3 ,&
              jupbnd , k , ksgn , ktmin , Last , Limit , maxerr , &
              Neval , nres , nrmax , numrl2
      LOGICAL extrap , noext
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
                , res3la(3) , Rlist(Limit) , rlist2(52)
!
      procedure(func) :: f
!
!            the dimension of rlist2 is determined by the value of
!            limexp in subroutine dqelg.
!
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least (limexp+2),
!                       containing the part of the epsilon table
!                       wich is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained, it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered up
!                       to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true-value)
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!           oflow is the largest positive magnitude.
!

      epmach = D1MACH(4)
!
!           test on validity of parameters
!           -----------------------------
!
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
      IF ( Epsabs<=0.0_wp .AND. Epsrel<max(50.0_wp*epmach,0.5e-28_wp) ) &
           Ier = 6
      IF ( Ier==6 ) return
!
!
!           first approximation to the integral
!           -----------------------------------
!
!           determine the interval to be mapped onto (0,1).
!           if inf = 2 the integral is computed as i = i1+i2, where
!           i1 = integral of f over (-infinity,0),
!           i2 = integral of f over (0,+infinity).
!
      boun = Bound
      IF ( Inf==2 ) boun = 0.0_wp
      CALL DQK15I(F,boun,Inf,0.0_wp,1.0_wp,Result,Abserr,defabs, &
                  resabs)
!
!           test on accuracy
!
      Last = 1
      Rlist(1) = Result
      Elist(1) = Abserr
      Iord(1) = 1
      dres = abs(Result)
      errbnd = max(Epsabs,Epsrel*dres)
      IF ( Abserr<=100.0_wp*epmach*defabs .AND. Abserr>errbnd ) Ier = 2
      IF ( Limit==1 ) Ier = 1
      IF ( Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) .OR. &
           Abserr==0.0_wp ) GOTO 400
!
!           initialization
!           --------------
!
      uflow = D1MACH(1)
      oflow = D1MACH(2)
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
      extrap = .FALSE.
      noext = .FALSE.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      IF ( dres>=(1.0_wp-50.0_wp*epmach)*defabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
      DO Last = 2 , Limit
!
!           bisect the subinterval with nrmax-th largest error estimate.
!
         a1 = Alist(maxerr)
         b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
         a2 = b1
         b2 = Blist(maxerr)
         erlast = errmax
         CALL DQK15I(F,boun,Inf,a1,b1,area1,error1,resabs,defab1)
         CALL DQK15I(F,boun,Inf,a2,b2,area2,error2,resabs,defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
         area12 = area1 + area2
         erro12 = error1 + error2
         errsum = errsum + erro12 - errmax
         area = area + area12 - Rlist(maxerr)
         IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
            IF ( abs(Rlist(maxerr)-area12)<=0.1e-4_wp*abs(area12) .AND. &
                 erro12>=0.99_wp*errmax ) THEN
               IF ( extrap ) iroff2 = iroff2 + 1
               IF ( .NOT.extrap ) iroff1 = iroff1 + 1
            ENDIF
            IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
         ENDIF
         Rlist(maxerr) = area1
         Rlist(Last) = area2
         errbnd = max(Epsabs,Epsrel*abs(area))
!
!           test for roundoff error and eventually set error flag.
!
         IF ( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
         IF ( iroff2>=5 ) ierro = 3
!
!           set error flag in the case that the number of
!           subintervals equals limit.
!
         IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at some points of the integration range.
!
         IF ( max(abs(a1),abs(b2))<=(1.0_wp+100.0_wp*epmach) &
              *(abs(a2)+1000.0_wp*uflow) ) Ier = 4
!
!           append the newly-created intervals to the list.
!
         IF ( error2>error1 ) THEN
            Alist(maxerr) = a2
            Alist(Last) = a1
            Blist(Last) = b1
            Rlist(maxerr) = area2
            Rlist(Last) = area1
            Elist(maxerr) = error2
            Elist(Last) = error1
         ELSE
            Alist(Last) = a2
            Blist(maxerr) = b1
            Blist(Last) = b2
            Elist(maxerr) = error1
            Elist(Last) = error2
         ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with nrmax-th largest error estimate (to be bisected next).
!
         CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
         IF ( errsum<=errbnd ) GOTO 300
         IF ( Ier/=0 ) GOTO 200
         IF ( Last==2 ) THEN
            small = 0.375_wp
            erlarg = errsum
            ertest = errbnd
            rlist2(2) = area
         ELSEIF ( .NOT.(noext) ) THEN
            erlarg = erlarg - erlast
            IF ( abs(b1-a1)>small ) erlarg = erlarg + erro12
            IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
               IF ( abs(Blist(maxerr)-Alist(maxerr))>small ) GOTO 100
               extrap = .TRUE.
               nrmax = 2
            ENDIF
            IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors over the
!           larger intervals (erlarg) and perform extrapolation.
!
               id = nrmax
               jupbnd = Last
               IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
               DO k = id , jupbnd
                  maxerr = Iord(nrmax)
                  errmax = Elist(maxerr)
                  IF ( abs(Blist(maxerr)-Alist(maxerr))>small ) &
                       GOTO 100
                  nrmax = nrmax + 1
               ENDDO
            ENDIF
!
!           perform extrapolation.
!
            numrl2 = numrl2 + 1
            rlist2(numrl2) = area
            CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
            ktmin = ktmin + 1
            IF ( ktmin>5 .AND. Abserr<0.1e-02_wp*errsum ) Ier = 5
            IF ( abseps<Abserr ) THEN
               ktmin = 0
               Abserr = abseps
               Result = reseps
               correc = erlarg
               ertest = max(Epsabs,Epsrel*abs(reseps))
               IF ( Abserr<=ertest ) GOTO 200
            ENDIF
!
!            prepare bisection of the smallest interval.
!
            IF ( numrl2==1 ) noext = .TRUE.
            IF ( Ier==5 ) GOTO 200
            maxerr = Iord(1)
            errmax = Elist(maxerr)
            nrmax = 1
            extrap = .FALSE.
            small = small*0.5_wp
            erlarg = errsum
         ENDIF
 100  ENDDO
!
!           set final result and error estimate.
!           ------------------------------------
!
 200  IF ( Abserr/=oflow ) THEN
         IF ( (Ier+ierro)/=0 ) THEN
            IF ( ierro==3 ) Abserr = Abserr + correc
            IF ( Ier==0 ) Ier = 3
            IF ( Result==0.0_wp .OR. area==0.0_wp ) THEN
               IF ( Abserr>errsum ) GOTO 300
               IF ( area==0.0_wp ) GOTO 400
            ELSEIF ( Abserr/abs(Result)>errsum/abs(area) ) THEN
               GOTO 300
            ENDIF
         ENDIF
!
!           test on divergence
!
         IF ( ksgn/=(-1) .OR. max(abs(Result),abs(area)) &
              >defabs*0.01_wp ) THEN
            IF ( 0.01_wp>(Result/area) .OR. (Result/area)>100.0_wp .OR. &
                 errsum>abs(area) ) Ier = 6
         ENDIF
         GOTO 400
      ENDIF
!
!           compute global integral sum.
!
 300  Result = 0.0_wp
      DO k = 1 , Last
         Result = Result + Rlist(k)
      ENDDO
      Abserr = errsum
 400  Neval = 30*Last - 15
      IF ( Inf==2 ) Neval = 2*Neval
      IF ( Ier>2 ) Ier = Ier - 1
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGP(F,A,B,Npts2,Points,Epsabs,Epsrel,Result,Abserr, &
                       Neval,Ier,Leniw,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             singularities at user specified points,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            break points of the integration interval, where local
!            difficulties of the integrand may occur (e.g.
!            singularities, discontinuities), are provided by the user.
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            npts2  - integer
!                     number equal to two more than the number of
!                     user-supplied break points within the integration
!                     range, npts>=2.
!                     if npts2<2, the routine will end with ier = 6.
!
!            points - real(wp)
!                     vector of dimension npts2, the first (npts2-2)
!                     elements of which are the user provided break
!                     points. if these points do not constitute an
!                     ascending sequence there will be an automatic
!                     sorting.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (i.e. singularity,
!                             discontinuity within the interval), it
!                             should be supplied to the routine as an
!                             element of the vector points. if necessary
!                             an appropriate special-purpose integrator
!                             must be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is presumed that the requested
!                             tolerance cannot be achieved, and that
!                             the returned result is the best which
!                             can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier>0.
!                         = 6 the input is invalid because
!                             npts2<2 or
!                             break points are specified outside
!                             the integration range or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             result, abserr, neval, last are set to
!                             zero. exept when leniw or lenw or npts2 is
!                             invalid, iwork(1), iwork(limit+1),
!                             work(limit*2+1) and work(limit*3+1)
!                             are set to zero.
!                             work(1) is set to a and work(limit+1)
!                             to b (where limit = (leniw-npts2)/2).
!
!         dimensioning parameters
!            leniw - integer
!                    dimensioning parameter for iwork
!                    leniw determines limit = (leniw-npts2)/2,
!                    which is the maximum number of subintervals in the
!                    partition of the given integration interval (a,b),
!                    leniw>=(3*npts2-2).
!                    if leniw<(3*npts2-2), the routine will end with
!                    ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least leniw*2-npts2.
!                    if lenw<leniw*2-npts2, the routine will end
!                    with ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least leniw. on return,
!                    the first k elements of which contain
!                    pointers to the error estimates over the
!                    subintervals, such that work(limit*3+iwork(1)),...,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2), and
!                    k = limit+1-last otherwise
!                    iwork(limit+1), ...,iwork(limit+last) contain the
!                     subdivision levels of the subintervals, i.e.
!                     if (aa,bb) is a subinterval of (p1,p2)
!                     where p1 as well as p2 is a user-provided
!                     break point or integration limit, then (aa,bb) has
!                     level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
!                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have
!                     no significance for the user,
!                    note that limit = (leniw-npts2)/2.
!
!            work  - real(wp)
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                     end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ..., work(limit*2+last) contain
!                     the integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3+last)
!                     contain the corresponding error estimates,
!                    work(limit*4+1), ..., work(limit*4+npts2)
!                     contain the integration limits and the
!                     break points sorted in an ascending sequence.
!                    note that limit = (leniw-npts2)/2.
!

!
      real(wp) A , Abserr , B , Epsabs , Epsrel , Points , &
                       Result , Work
      INTEGER Ier , Iwork , Last , Leniw , Lenw , limit , lvl , l1 , &
              l2 , l3 , l4 , Neval , Npts2
!
      DIMENSION Iwork(Leniw) , Points(Npts2) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of limit and lenw.
!

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Leniw>=(3*Npts2-2) .AND. Lenw>=(Leniw*2-Npts2) .AND. &
           Npts2>=2 ) THEN
!
!         prepare call for dqagpe.
!
         limit = (Leniw-Npts2)/2
         l1 = limit + 1
         l2 = limit + l1
         l3 = limit + l2
         l4 = limit + l3
!
         CALL DQAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,limit,Result, &
                     Abserr,Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3)&
                     ,Work(l4),Iwork(1),Iwork(l1),Iwork(l2),Last)
!
!         call error handler if necessary.
!
         lvl = 0
      else
         write(*,*) Leniw, (3*Npts2-2)
         write(*,*) Lenw, (Leniw*2-Npts2)
         write(*,*) Npts2,  2
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqagp',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,Limit,Result, &
                        Abserr,Neval,Ier,Alist,Blist,Rlist,Elist,Pts, &
                        Iord,Level,Ndin,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             singularities at user specified points,
!             extrapolation, globally adaptive.
!***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b), hopefully
!            satisfying following claim for accuracy abs(i-result)<=
!            max(epsabs,epsrel*abs(i)). break points of the integration
!            interval, where local difficulties of the integrand may
!            occur(e.g. singularities,discontinuities),provided by user.
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            npts2  - integer
!                     number equal to two more than the number of
!                     user-supplied break points within the integration
!                     range, npts2>=2.
!                     if npts2<2, the routine will end with ier = 6.
!
!            points - real(wp)
!                     vector of dimension npts2, the first (npts2-2)
!                     elements of which are the user provided break
!                     points. if these points do not constitute an
!                     ascending sequence there will be an automatic
!                     sorting.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=npts2
!                     if limit<npts2, the routine will end with
!                     ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (i.e. singularity,
!                             discontinuity within the interval), it
!                             should be supplied to the routine as an
!                             element of the vector points. if necessary
!                             an appropriate special-purpose integrator
!                             must be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table. it is presumed that
!                             the requested tolerance cannot be
!                             achieved, and that the returned result is
!                             the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier>0.
!                         = 6 the input is invalid because
!                             npts2<2 or
!                             break points are specified outside
!                             the integration range or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or limit<npts2.
!                             result, abserr, neval, last, rlist(1),
!                             and elist(1) are set to zero. alist(1) and
!                             blist(1) are set to a and b respectively.
!
!            alist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            blist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            rlist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            pts    - real(wp)
!                     vector of dimension at least npts2, containing the
!                     integration limits and the break points of the
!                     interval in ascending sequence.
!
!            level  - integer
!                     vector of dimension at least limit, containing the
!                     subdivision levels of the subinterval, i.e. if
!                     (aa,bb) is a subinterval of (p1,p2) where p1 as
!                     well as p2 is a user-provided break point or
!                     integration limit, then (aa,bb) has level l if
!                     abs(bb-aa) = abs(p2-p1)*2**(-l).
!
!            ndin   - integer
!                     vector of dimension at least npts2, after first
!                     integration over the intervals (pts(i)),pts(i+1),
!                     i = 0,1, ..., npts2-2, the error estimates over
!                     some of the intervals may have been increased
!                     artificially, in order to put their subdivision
!                     forward. if this happens for the subinterval
!                     numbered k, ndin(k) is put to 1, otherwise
!                     ndin(k) = 0.
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the
!                     error estimates over the subintervals,
!                     such that elist(iord(1)), ..., elist(iord(k))
!                     form a decreasing sequence, with k = last
!                     if last<=(limit/2+2), and k = limit+1-last
!                     otherwise
!
!            last   - integer
!                     number of subintervals actually produced in the
!                     subdivisions process
!

      real(wp) A , abseps , Abserr , Alist , area , area1 , &
                       area12 , area2 , a1 , a2 , B , Blist , b1 , b2 , &
                       correc , defabs , defab1 , defab2 , &
                       min , dres , Elist , epmach , &
                       Epsabs , Epsrel , erlarg , erlast , errbnd , &
                       errmax , error1 , erro12 , error2 , errsum , &
                       ertest , oflow , Points , Pts , resa , &
                       resabs , reseps , Result , res3la , Rlist , &
                       rlist2 , sign , temp , uflow
      INTEGER i , id , Ier , ierro , ind1 , ind2 , Iord , ip1 , iroff1 ,&
              iroff2 , iroff3 , j , jlow , jupbnd , k , ksgn , ktmin , &
              Last , levcur , Level , levmax , Limit , maxerr , Ndin , &
              Neval , nint , nintp1 , npts , Npts2 , nres , nrmax , &
              numrl2
      LOGICAL extrap , noext
!
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
                , Level(Limit) , Ndin(Npts2) , Points(Npts2) , &
                Pts(Npts2) , res3la(3) , Rlist(Limit) , rlist2(52)
!
      procedure(func) :: f
!
!            the dimension of rlist2 is determined by the value of
!            limexp in subroutine epsalg (rlist2 should be of dimension
!            (limexp+2) at least).
!
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2
!                       containing the part of the epsilon table which
!                       is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements in rlist2. if an appropriate
!                       approximation to the compounded integral has
!                       been obtained, it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one.
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation is
!                       no longer allowed (true-value)
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!           oflow is the largest positive magnitude.
!

      epmach = D1MACH(4)
!
!            test on validity of parameters
!            -----------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      Iord(1) = 0
      Level(1) = 0
      npts = Npts2 - 2
      IF ( Npts2<2 .OR. Limit<=npts .OR.  &
           (Epsabs<=0.0_wp .AND. Epsrel<max(50.0_wp*epmach,0.5e-28_wp)) &
           ) then
            write(*,*) 'Npts2=',Npts2
            write(*,*) 'Epsabs=',Epsabs
            write(*,*) 'Epsrel=',Epsabs
            write(*,*) 'max(50.0_wp*epmach,0.5e-28_wp)=',max(50.0_wp*epmach,0.5e-28_wp)
            write(*,*) 'oops'
            Ier = 6
      end if
      IF ( Ier==6 ) return
!
!            if any break points are provided, sort them into an
!            ascending sequence.
!
      sign = 1.0_wp
      IF ( A>B ) sign = -1.0_wp
      Pts(1) = min(A,B)
      IF ( npts/=0 ) THEN
         DO i = 1 , npts
            Pts(i+1) = Points(i)
         ENDDO
      ENDIF
      Pts(npts+2) = max(A,B)
      nint = npts + 1
      a1 = Pts(1)
      IF ( npts/=0 ) THEN
         nintp1 = nint + 1
         DO i = 1 , nint
            ip1 = i + 1
            DO j = ip1 , nintp1
               IF ( Pts(i)>Pts(j) ) THEN
                  temp = Pts(i)
                  Pts(i) = Pts(j)
                  Pts(j) = temp
               ENDIF
            ENDDO
         ENDDO
         IF ( Pts(1)/=min(A,B) .OR. Pts(nintp1)/=max(A,B) ) Ier = 6
         IF ( Ier==6 ) return
      ENDIF
!
!            compute first integral and error approximations.
!            ------------------------------------------------
!
      resabs = 0.0_wp
      DO i = 1 , nint
         b1 = Pts(i+1)
         CALL DQK21(F,a1,b1,area1,error1,defabs,resa)
         Abserr = Abserr + error1
         Result = Result + area1
         Ndin(i) = 0
         IF ( error1==resa .AND. error1/=0.0_wp ) Ndin(i) = 1
         resabs = resabs + defabs
         Level(i) = 0
         Elist(i) = error1
         Alist(i) = a1
         Blist(i) = b1
         Rlist(i) = area1
         Iord(i) = i
         a1 = b1
      ENDDO
      errsum = 0.0_wp
      DO i = 1 , nint
         IF ( Ndin(i)==1 ) Elist(i) = Abserr
         errsum = errsum + Elist(i)
      ENDDO
!
!           test on accuracy.
!
      Last = nint
      Neval = 21*nint
      dres = abs(Result)
      errbnd = max(Epsabs,Epsrel*dres)
      IF ( Abserr<=100.0_wp*epmach*resabs .AND. Abserr>errbnd ) Ier = 2
      IF ( nint/=1 ) THEN
         DO i = 1 , npts
            jlow = i + 1
            ind1 = Iord(i)
            DO j = jlow , nint
               ind2 = Iord(j)
               IF ( Elist(ind1)<=Elist(ind2) ) THEN
                  ind1 = ind2
                  k = j
               ENDIF
            ENDDO
            IF ( ind1/=Iord(i) ) THEN
               Iord(k) = Iord(i)
               Iord(i) = ind1
            ENDIF
         ENDDO
         IF ( Limit<Npts2 ) Ier = 1
      ENDIF
      IF ( Ier/=0 .OR. Abserr<=errbnd ) GOTO 400
!
!           initialization
!           --------------
!
      rlist2(1) = Result
      maxerr = Iord(1)
      errmax = Elist(maxerr)
      area = Result
      nrmax = 1
      nres = 0
      numrl2 = 1
      ktmin = 0
      extrap = .FALSE.
      noext = .FALSE.
      erlarg = errsum
      ertest = errbnd
      levmax = 1
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ierro = 0
      uflow = D1MACH(1)
      oflow = D1MACH(2)
      Abserr = oflow
      ksgn = -1
      IF ( dres>=(1.0_wp-50.0_wp*epmach)*resabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
      DO Last = Npts2 , Limit
!
!           bisect the subinterval with the nrmax-th largest error
!           estimate.
!
         levcur = Level(maxerr) + 1
         a1 = Alist(maxerr)
         b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
         a2 = b1
         b2 = Blist(maxerr)
         erlast = errmax
         CALL DQK21(F,a1,b1,area1,error1,resa,defab1)
         CALL DQK21(F,a2,b2,area2,error2,resa,defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
         Neval = Neval + 42
         area12 = area1 + area2
         erro12 = error1 + error2
         errsum = errsum + erro12 - errmax
         area = area + area12 - Rlist(maxerr)
         IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
            IF ( abs(Rlist(maxerr)-area12)<=0.1e-4_wp*abs(area12) .AND. &
                 erro12>=0.99_wp*errmax ) THEN
               IF ( extrap ) iroff2 = iroff2 + 1
               IF ( .NOT.extrap ) iroff1 = iroff1 + 1
            ENDIF
            IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
         ENDIF
         Level(maxerr) = levcur
         Level(Last) = levcur
         Rlist(maxerr) = area1
         Rlist(Last) = area2
         errbnd = max(Epsabs,Epsrel*abs(area))
!
!           test for roundoff error and eventually set error flag.
!
         IF ( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
         IF ( iroff2>=5 ) ierro = 3
!
!           set error flag in the case that the number of
!           subintervals equals limit.
!
         IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range
!
         IF ( max(abs(a1),abs(b2))<=(1.0_wp+100.0_wp*epmach) &
              *(abs(a2)+1000.0_wp*uflow) ) Ier = 4
!
!           append the newly-created intervals to the list.
!
         IF ( error2>error1 ) THEN
            Alist(maxerr) = a2
            Alist(Last) = a1
            Blist(Last) = b1
            Rlist(maxerr) = area2
            Rlist(Last) = area1
            Elist(maxerr) = error2
            Elist(Last) = error1
         ELSE
            Alist(Last) = a2
            Blist(maxerr) = b1
            Blist(Last) = b2
            Elist(maxerr) = error1
            Elist(Last) = error2
         ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with nrmax-th largest error estimate (to be bisected next).
!
         CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
         IF ( errsum<=errbnd ) GOTO 300
! ***jump out of do-loop
         IF ( Ier/=0 ) GOTO 200
         IF ( .NOT.(noext) ) THEN
            erlarg = erlarg - erlast
            IF ( levcur+1<=levmax ) erlarg = erlarg + erro12
            IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
               IF ( Level(maxerr)+1<=levmax ) GOTO 100
               extrap = .TRUE.
               nrmax = 2
            ENDIF
            IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors over
!           the larger intervals (erlarg) and perform extrapolation.
!
               id = nrmax
               jupbnd = Last
               IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
               DO k = id , jupbnd
                  maxerr = Iord(nrmax)
                  errmax = Elist(maxerr)
! ***jump out of do-loop
                  IF ( Level(maxerr)+1<=levmax ) GOTO 100
                  nrmax = nrmax + 1
               ENDDO
            ENDIF
!
!           perform extrapolation.
!
            numrl2 = numrl2 + 1
            rlist2(numrl2) = area
            IF ( numrl2>2 ) THEN
               CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
               ktmin = ktmin + 1
               IF ( ktmin>5 .AND. Abserr<0.1e-02_wp*errsum ) Ier = 5
               IF ( abseps<Abserr ) THEN
                  ktmin = 0
                  Abserr = abseps
                  Result = reseps
                  correc = erlarg
                  ertest = max(Epsabs,Epsrel*abs(reseps))
! ***jump out of do-loop
                  IF ( Abserr<ertest ) GOTO 200
               ENDIF
!
!           prepare bisection of the smallest interval.
!
               IF ( numrl2==1 ) noext = .TRUE.
               IF ( Ier>=5 ) GOTO 200
            ENDIF
            maxerr = Iord(1)
            errmax = Elist(maxerr)
            nrmax = 1
            extrap = .FALSE.
            levmax = levmax + 1
            erlarg = errsum
         ENDIF
 100  ENDDO
!
!           set the final result.
!           ---------------------
!
!
 200  IF ( Abserr/=oflow ) THEN
         IF ( (Ier+ierro)/=0 ) THEN
            IF ( ierro==3 ) Abserr = Abserr + correc
            IF ( Ier==0 ) Ier = 3
            IF ( Result==0.0_wp .OR. area==0.0_wp ) THEN
               IF ( Abserr>errsum ) GOTO 300
               IF ( area==0.0_wp ) GOTO 400
            ELSEIF ( Abserr/abs(Result)>errsum/abs(area) ) THEN
               GOTO 300
            ENDIF
         ENDIF
!
!           test on divergence.
!
         IF ( ksgn/=(-1) .OR. max(abs(Result),abs(area)) &
              >resabs*0.01_wp ) THEN
            IF ( 0.01_wp>(Result/area) .OR. (Result/area)>100.0_wp .OR. &
                 errsum>abs(area) ) Ier = 6
         ENDIF
         GOTO 400
      ENDIF
!
!           compute global integral sum.
!
 300  Result = 0.0_wp
      DO k = 1 , Last
         Result = Result + Rlist(k)
      ENDDO
      Abserr = errsum
 400  IF ( Ier>2 ) Ier = Ier - 1
      Result = Result*sign
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGS(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier, &
                       Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             (end-point) singularities, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & prog. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral  i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of limit
!                             (and taking the according dimension
!                             adjustments into account. however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour
!                             occurs at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table. it is presumed that
!                             the requested tolerance cannot be
!                             achieved, and that the returned result is
!                             the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp)
!                             or limit<1 or lenw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero.except when limit or lenw is invalid,
!                             iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit>=1.
!                    if limit<1, the routine will end with ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end
!                    with ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, detemines the
!                    number of significant elements actually in the work
!                    arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first k
!                    elements of which contain pointers
!                    to the error estimates over the subintervals
!                    such that work(limit*3+iwork(1)),... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2),
!                    and k = limit+1-last otherwise
!
!            work  - real(wp)
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                     end-points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end-points,
!                    work(limit*2+1), ..., work(limit*2+last) contain
!                     the integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3+last)
!                     contain the error estimates.
!

!
!
      real(wp) A , Abserr , B , Epsabs , Epsrel , Result , &
                       Work
      INTEGER Ier , Iwork , Last , Lenw , Limit , lvl , l1 , l2 , l3 , &
              Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of limit and lenw.
!

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for dqagse.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL DQAGSE(F,A,B,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier, &
                     Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqags',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAGSE(F,A,B,Epsabs,Epsrel,Limit,Result,Abserr,Neval, &
                        Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, general-purpose,
!             (end point) singularities, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b)
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of limit
!                             (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour
!                             occurs at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs<=0 and
!                             epsrel<max(50*rel.mach.acc.,0.5e-28_wp).
!                             result, abserr, neval, last, rlist(1),
!                             iord(1) and elist(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!            alist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left end points
!                     of the subintervals in the partition of the
!                     given integration range (a,b)
!
!            blist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            rlist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the
!                     error estimates over the subintervals,
!                     such that elist(iord(1)), ..., elist(iord(k))
!                     form a decreasing sequence, with k = last
!                     if last<=(limit/2+2), and k = limit+1-last
!                     otherwise
!
!            last   - integer
!                     number of subintervals actually produced in the
!                     subdivision process
!

!
      real(wp) A , abseps , Abserr , Alist , area , area1 , &
                       area12 , area2 , a1 , a2 , B , Blist , b1 , b2 , &
                       correc , defabs , defab1 , defab2 , &
                       dres , Elist , epmach , Epsabs ,&
                       Epsrel , erlarg , erlast , errbnd , errmax , &
                       error1 , error2 , erro12 , errsum , ertest , &
                       oflow , resabs , reseps , Result , res3la , &
                       Rlist , rlist2 , small , uflow
      INTEGER id , Ier , ierro , Iord , iroff1 , iroff2 , iroff3 , &
              jupbnd , k , ksgn , ktmin , Last , Limit , maxerr , &
              Neval , nres , nrmax , numrl2
      LOGICAL extrap , noext
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
                , res3la(3) , Rlist(Limit) , rlist2(52)
!
      procedure(func) :: f
!
!            the dimension of rlist2 is determined by the value of
!            limexp in subroutine dqelg (rlist2 should be of dimension
!            (limexp+2) at least).
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered up
!                       to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!           oflow is the largest positive magnitude.
!

      epmach = D1MACH(4)
!
!            test on validity of parameters
!            ------------------------------
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      IF ( Epsabs<=0.0_wp .AND. Epsrel<max(50.0_wp*epmach,0.5e-28_wp) ) &
           Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         uflow = D1MACH(1)
         oflow = D1MACH(2)
         ierro = 0
         CALL DQK21(F,A,B,Result,Abserr,defabs,resabs)
!
!           test on accuracy.
!
         dres = abs(Result)
         errbnd = max(Epsabs,Epsrel*dres)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         IF ( Abserr<=100.0_wp*epmach*defabs .AND. Abserr>errbnd ) &
              Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) .OR. &
              Abserr==0.0_wp ) THEN
            Neval = 42*Last - 21
            return
         ELSE
!
!           initialization
!           --------------
!
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
            extrap = .FALSE.
            noext = .FALSE.
            iroff1 = 0
            iroff2 = 0
            iroff3 = 0
            ksgn = -1
            IF ( dres>=(1.0_wp-50.0_wp*epmach)*defabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
            DO Last = 2 , Limit
!
!           bisect the subinterval with the nrmax-th largest error
!           estimate.
!
               a1 = Alist(maxerr)
               b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               erlast = errmax
               CALL DQK21(F,a1,b1,area1,error1,resabs,defab1)
               CALL DQK21(F,a2,b2,area2,error2,resabs,defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
                  IF ( abs(Rlist(maxerr)-area12)<=0.1e-4_wp*abs(area12) &
                       .AND. erro12>=0.99_wp*errmax ) THEN
                     IF ( extrap ) iroff2 = iroff2 + 1
                     IF ( .NOT.extrap ) iroff1 = iroff1 + 1
                  ENDIF
                  IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = max(Epsabs,Epsrel*abs(area))
!
!           test for roundoff error and eventually set error flag.
!
               IF ( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
               IF ( iroff2>=5 ) ierro = 3
!
!           set error flag in the case that the number of subintervals
!           equals limit.
!
               IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range.
!
               IF ( max(abs(a1),abs(b2))<=(1.0_wp+100.0_wp*epmach) &
                    *(abs(a2)+1000.0_wp*uflow) ) Ier = 4
!
!           append the newly-created intervals to the list.
!
               IF ( error2>error1 ) THEN
                  Alist(maxerr) = a2
                  Alist(Last) = a1
                  Blist(Last) = b1
                  Rlist(maxerr) = area2
                  Rlist(Last) = area1
                  Elist(maxerr) = error2
                  Elist(Last) = error1
               ELSE
                  Alist(Last) = a2
                  Blist(maxerr) = b1
                  Blist(Last) = b2
                  Elist(maxerr) = error1
                  Elist(Last) = error2
               ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with nrmax-th largest error estimate (to be bisected next).
!
               CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( errsum<=errbnd ) GOTO 50
! ***jump out of do-loop
               IF ( Ier/=0 ) GOTO 40
               IF ( Last==2 ) THEN
                  small = abs(B-A)*0.375_wp
                  erlarg = errsum
                  ertest = errbnd
                  rlist2(2) = area
               ELSEIF ( .NOT.(noext) ) THEN
                  erlarg = erlarg - erlast
                  IF ( abs(b1-a1)>small ) erlarg = erlarg + erro12
                  IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
                     IF ( abs(Blist(maxerr)-Alist(maxerr))>small ) GOTO 20
                     extrap = .TRUE.
                     nrmax = 2
                  ENDIF
                  IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors over the
!           larger intervals (erlarg) and perform extrapolation.
!
                     id = nrmax
                     jupbnd = Last
                     IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
                     DO k = id , jupbnd
                        maxerr = Iord(nrmax)
                        errmax = Elist(maxerr)
! ***jump out of do-loop
                        IF ( abs(Blist(maxerr)-Alist(maxerr))>small ) GOTO 20
                        nrmax = nrmax + 1
                     ENDDO
                  ENDIF
!
!           perform extrapolation.
!
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
                  CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
                  ktmin = ktmin + 1
                  IF ( ktmin>5 .AND. Abserr<0.1e-02_wp*errsum ) Ier = 5
                  IF ( abseps<Abserr ) THEN
                     ktmin = 0
                     Abserr = abseps
                     Result = reseps
                     correc = erlarg
                     ertest = max(Epsabs,Epsrel*abs(reseps))
! ***jump out of do-loop
                     IF ( Abserr<=ertest ) GOTO 40
                  ENDIF
!
!           prepare bisection of the smallest interval.
!
                  IF ( numrl2==1 ) noext = .TRUE.
                  IF ( Ier==5 ) GOTO 40
                  maxerr = Iord(1)
                  errmax = Elist(maxerr)
                  nrmax = 1
                  extrap = .FALSE.
                  small = small*0.5_wp
                  erlarg = errsum
               ENDIF
 20         ENDDO
!
!           set final result and error estimate.
!           ------------------------------------
!
 40         IF ( Abserr/=oflow ) THEN
               IF ( Ier+ierro/=0 ) THEN
                  IF ( ierro==3 ) Abserr = Abserr + correc
                  IF ( Ier==0 ) Ier = 3
                  IF ( Result==0.0_wp .OR. area==0.0_wp ) THEN
                     IF ( Abserr>errsum ) GOTO 50
                     IF ( area==0.0_wp ) THEN
                        IF ( Ier>2 ) Ier = Ier - 1
                        Neval = 42*Last - 21
                        return
                     ENDIF
                  ELSEIF ( Abserr/abs(Result)>errsum/abs(area) ) THEN
                     GOTO 50
                  ENDIF
               ENDIF
!
!           test on divergence.
!
               IF ( ksgn/=(-1) .OR. max(abs(Result),abs(area)) &
                    >defabs*0.01_wp ) THEN
                  IF ( 0.01_wp>(Result/area) .OR. (Result/area) &
                       >100.0_wp .OR. errsum>abs(area) ) Ier = 6
               ENDIF
               IF ( Ier>2 ) Ier = Ier - 1
               Neval = 42*Last - 21
               return
            ENDIF
         ENDIF
!
!           compute global integral sum.
!
 50      Result = 0.0_wp
         DO k = 1 , Last
            Result = Result + Rlist(k)
         ENDDO
         Abserr = errsum
         IF ( Ier>2 ) Ier = Ier - 1
         Neval = 42*Last - 21
      ENDIF
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWC(F,A,B,C,Epsabs,Epsrel,Result,Abserr,Neval,Ier, &
                       Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             cauchy principal value,
!             clenshaw-curtis, globally adaptive
!***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a
!            cauchy principal value i = integral of f*w over (a,b)
!            (w(x) = 1/((x-c), c/=a, c/=b), hopefully satisfying
!            following claim for accuracy
!            abs(i-result)<=max(epsabe,epsrel*abs(i)).
!***description
!
!        computation of a cauchy principal value
!        standard fortran subroutine
!        real(wp) version
!
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     under limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            c      - parameter in the weight function, c/=a, c/=b.
!                     if c = a or c = b, the routine will end with
!                     ier = 6 .
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate or the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of limit
!                             (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty
!                             can be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or limit<1 or lenw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero. exept when lenw or limit is invalid,
!                             iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit>=1.
!                    if limit<1, the routine will end with ier = 6.
!
!           lenw   - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end with
!                    ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first k
!                    elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that work(limit*3+iwork(1)), ... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2),
!                    and k = limit+1-last otherwise
!
!            work  - real(wp)
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                     end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ..., work(limit*2+last) contain
!                     the integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3+last)
!                     contain the error estimates.
!

!
      real(wp) A , Abserr , B , C , Epsabs , Epsrel , &
                       Result , Work
      INTEGER Ier , Iwork , Last , Lenw , Limit , lvl , l1 , l2 , l3 , &
              Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of limit and lenw.
!

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for dqawce.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
         CALL DQAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval, &
                     Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqawc',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval,&
                        Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             cauchy principal value, clenshaw-curtis method
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***  purpose  the routine calculates an approximation result to a
!              cauchy principal value i = integral of f*w over (a,b)
!              (w(x) = 1/(x-c), (c/=a, c/=b), hopefully satisfying
!              following claim for accuracy
!              abs(i-result)<=max(epsabs,epsrel*abs(i))
!***description
!
!        computation of a cauchy principal value
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            c      - real(wp)
!                     parameter in the weight function, c/=a, c/=b
!                     if c = a or c = b, the routine will end with
!                     ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=1
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement it is advised to analyze the
!                             the integrand, in order to determine the
!                             the integration difficulties. if the
!                             position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour
!                             occurs at some interior points of
!                             the integration interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or limit<1.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero. alist(1)
!                             and blist(1) are set to a and b
!                             respectively.
!
!            alist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            blist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            rlist   - real(wp)
!                      vector of dimension at least limit, the first
!                       last  elements of which are the integral
!                      approximations on the subintervals
!
!            elist   - real(wp)
!                      vector of dimension limit, the first  last
!                      elements of which are the moduli of the absolute
!                      error estimates on the subintervals
!
!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the error
!                      estimates over the subintervals, so that
!                      elist(iord(1)), ..., elist(iord(k)) with k = last
!                      if last<=(limit/2+2), and k = limit+1-last
!                      otherwise, form a decreasing sequence
!
!            last    - integer
!                      number of subintervals actually produced in
!                      the subdivision process
!

!
      real(wp) A , aa , Abserr , Alist , area , area1 , area12 ,&
                       area2 , a1 , a2 , B , bb , Blist , b1 , b2 , C , &
                       abs , Elist , epmach , Epsabs ,&
                       Epsrel , errbnd , errmax , error1 , erro12 , &
                       error2 , errsum , Result , Rlist , uflow
      INTEGER Ier , Iord , iroff1 , iroff2 , k , krule , Last , Limit , &
              maxerr , nev , Neval , nrmax
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) , &
                Elist(Limit) , Iord(Limit)
!
      procedure(func) :: f
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
!
!           test on validity of parameters
!           ------------------------------
!
      Ier = 6
      Neval = 0
      Last = 0
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      Iord(1) = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( .NOT.(C==A .OR. C==B .OR. (Epsabs<=0.0_wp .AND. Epsrel<max &
           (50.0_wp*epmach,0.5e-28_wp))) ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         aa = A
         bb = B
         IF ( A>B ) THEN
            aa = B
            bb = A
         ENDIF
         Ier = 0
         krule = 1
         CALL DQC25C(F,aa,bb,C,Result,Abserr,krule,Neval)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         Alist(1) = A
         Blist(1) = B
!
!           test on accuracy
!
         errbnd = max(Epsabs,Epsrel*abs(Result))
         IF ( Limit==1 ) Ier = 1
         IF ( Abserr>=min(0.01_wp*abs(Result),errbnd) .AND. Ier/=1 ) THEN
!
!           initialization
!           --------------
!
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
!
!           main do-loop
!           ------------
!
            DO Last = 2 , Limit
!
!           bisect the subinterval with nrmax-th largest
!           error estimate.
!
               a1 = Alist(maxerr)
               b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
               b2 = Blist(maxerr)
               IF ( C<=b1 .AND. C>a1 ) b1 = 0.5_wp*(C+b2)
               IF ( C>b1 .AND. C<b2 ) b1 = 0.5_wp*(a1+C)
               a2 = b1
               krule = 2
               CALL DQC25C(F,a1,b1,C,area1,error1,krule,nev)
               Neval = Neval + nev
               CALL DQC25C(F,a2,b2,C,area2,error2,krule,nev)
               Neval = Neval + nev
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( abs(Rlist(maxerr)-area12)<0.1e-4_wp*abs(area12) &
                    .AND. erro12>=0.99_wp*errmax .AND. krule==0 ) &
                    iroff1 = iroff1 + 1
               IF ( Last>10 .AND. erro12>errmax .AND. krule==0 ) &
                    iroff2 = iroff2 + 1
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = max(Epsabs,Epsrel*abs(area))
               IF ( errsum>errbnd ) THEN
!
!           test for roundoff error and eventually set error flag.
!
                  IF ( iroff1>=6 .AND. iroff2>20 ) Ier = 2
!
!           set error flag in the case that number of interval
!           bisections exceeds limit.
!
                  IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range.
!
                  IF ( max(abs(a1),abs(b2))    &
                       <=(1.0_wp+100.0_wp*epmach)  &
                       *(abs(a2)+1000.0_wp*uflow) ) Ier = 3
               ENDIF
!
!           append the newly-created intervals to the list.
!
               IF ( error2>error1 ) THEN
                  Alist(maxerr) = a2
                  Alist(Last) = a1
                  Blist(Last) = b1
                  Rlist(maxerr) = area2
                  Rlist(Last) = area1
                  Elist(maxerr) = error2
                  Elist(Last) = error1
               ELSE
                  Alist(Last) = a2
                  Blist(maxerr) = b1
                  Blist(Last) = b2
                  Elist(maxerr) = error1
                  Elist(Last) = error2
               ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with nrmax-th largest error estimate (to be bisected next).
!
               CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0_wp
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
         IF ( aa==B ) Result = -Result
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWF(F,A,Omega,Integr,Epsabs,Result,Abserr,Neval,Ier, &
                       Limlst,Lst,Leniw,Maxp1,Lenw,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,fourier
!             integral, integration between zeros with dqawoe,
!             convergence acceleration with dqelg
!***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            fourier integral i=integral of f(x)*w(x) over (a,infinity)
!            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=epsabs.
!***description
!
!        computation of fourier integrals
!        standard fortran subroutine
!        real(wp) version
!
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            omega  - real(wp)
!                     parameter in the integrand weight function
!
!            integr - integer
!                     indicates which of the weight functions is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine
!                     will end with ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested, epsabs>0.
!                     if epsabs<=0, the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                    if omega/=0
!                     ier = 1 maximum number of cycles allowed
!                             has been achieved, i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ..., lst.
!                             one can allow more cycles by increasing
!                             the value of limlst (and taking the
!                             according dimension adjustments into
!                             account). examine the array iwork which
!                             contains the error flags on the cycles, in
!                             order to look for eventual local
!                             integration difficulties.
!                             if the position of a local difficulty
!                             can be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 4 the extrapolation table constructed for
!                             convergence accelaration of the series
!                             formed by the integral contributions over
!                             the cycles, does not converge to within
!                             the requested accuracy.
!                             as in the case of ier = 1, it is advised
!                             to examine the array iwork which contains
!                             the error flags on the cycles.
!                         = 6 the input is invalid because
!                             (integr/=1 and integr/=2) or
!                              epsabs<=0 or limlst<1 or
!                              leniw<(limlst+2) or maxp1<1 or
!                              lenw<(leniw*2+maxp1*25).
!                              result, abserr, neval, lst are set to
!                              zero.
!                         = 7 bad integrand behaviour occurs within
!                             one or more of the cycles. location and
!                             type of the difficulty involved can be
!                             determined from the first lst elements of
!                             vector iwork.  here lst is the number of
!                             cycles actually needed (see below).
!                             iwork(k) = 1 the maximum number of
!                                          subdivisions (=(leniw-limlst)
!                                          /2) has been achieved on the
!                                          k th cycle.
!                                      = 2 occurrence of roundoff error
!                                          is detected and prevents the
!                                          tolerance imposed on the k th
!                                          cycle, from being achieved
!                                          on this cycle.
!                                      = 3 extremely bad integrand
!                                          behaviour occurs at some
!                                          points of the k th cycle.
!                                      = 4 the integration procedure
!                                          over the k th cycle does
!                                          not converge (to within the
!                                          required accuracy) due to
!                                          roundoff in the extrapolation
!                                          procedure invoked on this
!                                          cycle. it is assumed that the
!                                          result on this interval is
!                                          the best which can be
!                                          obtained.
!                                      = 5 the integral over the k th
!                                          cycle is probably divergent
!                                          or slowly convergent. it must
!                                          be noted that divergence can
!                                          occur with any other value of
!                                          iwork(k).
!                    if omega = 0 and integr = 1,
!                    the integral is calculated by means of dqagie,
!                    and ier = iwork(1) (with meaning as described
!                    for iwork(k),k = 1).
!
!         dimensioning parameters
!            limlst - integer
!                     limlst gives an upper bound on the number of
!                     cycles, limlst>=3.
!                     if limlst<3, the routine will end with ier = 6.
!
!            lst    - integer
!                     on return, lst indicates the number of cycles
!                     actually needed for the integration.
!                     if omega = 0, then lst is set to 1.
!
!            leniw  - integer
!                     dimensioning parameter for iwork. on entry,
!                     (leniw-limlst)/2 equals the maximum number of
!                     subintervals allowed in the partition of each
!                     cycle, leniw>=(limlst+2).
!                     if leniw<(limlst+2), the routine will end with
!                     ier = 6.
!
!            maxp1  - integer
!                     maxp1 gives an upper bound on the number of
!                     chebyshev moments which can be stored, i.e. for
!                     the intervals of lengths abs(b-a)*2**(-l),
!                     l = 0,1, ..., maxp1-2, maxp1>=1.
!                     if maxp1<1, the routine will end with ier = 6.
!            lenw   - integer
!                     dimensioning parameter for work
!                     lenw must be at least leniw*2+maxp1*25.
!                     if lenw<(leniw*2+maxp1*25), the routine will
!                     end with ier = 6.
!
!         work arrays
!            iwork  - integer
!                     vector of dimension at least leniw
!                     on return, iwork(k) for k = 1, 2, ..., lst
!                     contain the error flags on the cycles.
!
!            work   - real(wp)
!                     vector of dimension at least
!                     on return,
!                     work(1), ..., work(lst) contain the integral
!                      approximations over the cycles,
!                     work(limlst+1), ..., work(limlst+lst) contain
!                      the error extimates over the cycles.
!                     further elements of work have no specific
!                     meaning for the user.
!

!
      real(wp) A , Abserr , Epsabs , Omega , Result , Work
      INTEGER Ier , Integr , Iwork , last , Leniw , Lenw , limit , &
              Limlst , ll2 , lvl , Lst , l1 , l2 , l3 , l4 , l5 , l6 , &
              Maxp1 , Neval
!
      DIMENSION Iwork(Leniw) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of limlst, leniw, maxp1 and lenw.
!

      Ier = 6
      Neval = 0
      last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Limlst>=3 .AND. Leniw>=(Limlst+2) .AND. Maxp1>=1 .AND. &
           Lenw>=(Leniw*2+Maxp1*25) ) THEN
!
!         prepare call for dqawfe
!
         limit = (Leniw-Limlst)/2
         l1 = Limlst + 1
         l2 = Limlst + l1
         l3 = limit + l2
         l4 = limit + l3
         l5 = limit + l4
         l6 = limit + l5
         ll2 = limit + l1
         CALL DQAWFE(F,A,Omega,Integr,Epsabs,Limlst,limit,Maxp1,Result, &
                     Abserr,Neval,Ier,Work(1),Work(l1),Iwork(1),Lst, &
                     Work(l2),Work(l3),Work(l4),Work(l5),Iwork(l1), &
                     Iwork(ll2),Work(l6))
!
!         call error handler if necessary
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqawf',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWFE(F,A,Omega,Integr,Epsabs,Limlst,Limit,Maxp1, &
                        Result,Abserr,Neval,Ier,Rslst,Erlst,Ierlst,Lst, &
                        Alist,Blist,Rlist,Elist,Iord,Nnlog,Chebmo)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             fourier integrals,
!             integration between zeros with dqawoe,
!             convergence acceleration with dqelg
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a
!            given fourier integal
!            i = integral of f(x)*w(x) over (a,infinity)
!            where w(x)=cos(omega*x) or w(x)=sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=epsabs.
!***description
!
!        computation of fourier integrals
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to
!                     be declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            omega  - real(wp)
!                     parameter in the weight function
!
!            integr - integer
!                     indicates which weight function is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine will
!                     end with ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested, epsabs>0
!                     if epsabs<=0, the routine will end with ier = 6.
!
!            limlst - integer
!                     limlst gives an upper bound on the number of
!                     cycles, limlst>=1.
!                     if limlst<3, the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     allowed in the partition of each cycle, limit>=1
!                     each cycle, limit>=1.
!
!            maxp1  - integer
!                     gives an upper bound on the number of
!                     chebyshev moments which can be stored, i.e.
!                     for the intervals of lengths abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1>=1
!
!         on return
!            result - real(wp)
!                     approximation to the integral x
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - ier = 0 normal and reliable termination of
!                             the routine. it is assumed that the
!                             requested accuracy has been achieved.
!                     ier>0 abnormal termination of the routine. the
!                             estimates for integral and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                    if omega/=0
!                     ier = 1 maximum number of  cycles  allowed
!                             has been achieved., i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ..., lst.
!                             one can allow more cycles by increasing
!                             the value of limlst (and taking the
!                             according dimension adjustments into
!                             account).
!                             examine the array iwork which contains
!                             the error flags on the cycles, in order to
!                             look for eventual local integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling appropriate integrators on
!                             the subranges.
!                         = 4 the extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions over
!                             the cycles, does not converge to within
!                             the requested accuracy. as in the case of
!                             ier = 1, it is advised to examine the
!                             array iwork which contains the error
!                             flags on the cycles.
!                         = 6 the input is invalid because
!                             (integr/=1 and integr/=2) or
!                              epsabs<=0 or limlst<3.
!                              result, abserr, neval, lst are set
!                              to zero.
!                         = 7 bad integrand behaviour occurs within one
!                             or more of the cycles. location and type
!                             of the difficulty involved can be
!                             determined from the vector ierlst. here
!                             lst is the number of cycles actually
!                             needed (see below).
!                             ierlst(k) = 1 the maximum number of
!                                           subdivisions (= limit) has
!                                           been achieved on the k th
!                                           cycle.
!                                       = 2 occurrence of roundoff error
!                                           is detected and prevents the
!                                           tolerance imposed on the
!                                           k th cycle, from being
!                                           achieved.
!                                       = 3 extremely bad integrand
!                                           behaviour occurs at some
!                                           points of the k th cycle.
!                                       = 4 the integration procedure
!                                           over the k th cycle does
!                                           not converge (to within the
!                                           required accuracy) due to
!                                           roundoff in the
!                                           extrapolation procedure
!                                           invoked on this cycle. it
!                                           is assumed that the result
!                                           on this interval is the
!                                           best which can be obtained.
!                                       = 5 the integral over the k th
!                                           cycle is probably divergent
!                                           or slowly convergent. it
!                                           must be noted that
!                                           divergence can occur with
!                                           any other value of
!                                           ierlst(k).
!                    if omega = 0 and integr = 1,
!                    the integral is calculated by means of dqagie
!                    and ier = ierlst(1) (with meaning as described
!                    for ierlst(k), k = 1).
!
!            rslst  - real(wp)
!                     vector of dimension at least limlst
!                     rslst(k) contains the integral contribution
!                     over the interval (a+(k-1)c,a+kc) where
!                     c = (2*int(abs(omega))+1)*pi/abs(omega),
!                     k = 1, 2, ..., lst.
!                     note that, if omega = 0, rslst(1) contains
!                     the value of the integral over (a,infinity).
!
!            erlst  - real(wp)
!                     vector of dimension at least limlst
!                     erlst(k) contains the error estimate corresponding
!                     with rslst(k).
!
!            ierlst - integer
!                     vector of dimension at least limlst
!                     ierlst(k) contains the error flag corresponding
!                     with rslst(k). for the meaning of the local error
!                     flags see description of output parameter ier.
!
!            lst    - integer
!                     number of subintervals needed for the integration
!                     if omega = 0 then lst is set to 1.
!
!            alist, blist, rlist, elist - real(wp)
!                     vector of dimension at least limit,
!
!            iord, nnlog - integer
!                     vector of dimension at least limit, providing
!                     space for the quantities needed in the subdivision
!                     process of each cycle
!
!            chebmo - real(wp)
!                     array of dimension at least (maxp1,25), providing
!                     space for the chebyshev moments needed within the
!                     cycles
!

!
      real(wp) A , abseps , Abserr , Alist , Blist , Chebmo , &
                       correc , cycle , c1 , c2 , dl , dla , &
                       drl , Elist , Erlst , ep , eps ,&
                       epsa , Epsabs , errsum , fact , Omega , &
                       p1 , psum , reseps , Result , res3la , &
                       Rlist , Rslst , uflow
      INTEGER Ier , Ierlst , Integr , Iord , ktmin , l , last , Lst , &
              Limit , Limlst , ll , Maxp1 , momcom , nev , Neval , &
              Nnlog , nres , numrl2
!
      DIMENSION Alist(Limit) , Blist(Limit) , Chebmo(Maxp1,25) , &
                Elist(Limit) , Erlst(Limlst) , Ierlst(Limlst) , &
                Iord(Limit) , Nnlog(Limit) , psum(52) , res3la(3) , &
                Rlist(Limit) , Rslst(Limlst)
!
      procedure(func) :: f
!
!
!            the dimension of  psum  is determined by the value of
!            limexp in subroutine dqelg (psum must be of dimension
!            (limexp+2) at least).
!
!           list of major variables
!           -----------------------
!
!           c1, c2    - end points of subinterval (of length cycle)
!           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
!           psum      - vector of dimension at least (limexp+2)
!                       (see routine dqelg)
!                       psum contains the part of the epsilon table
!                       which is still needed for further computations.
!                       each element of psum is a partial sum of the
!                       series which should sum to the value of the
!                       integral.
!           errsum    - sum of error estimates over the subintervals,
!                       calculated cumulatively
!           epsa      - absolute tolerance requested over current
!                       subinterval
!           chebmo    - array containing the modified chebyshev
!                       moments (see also routine dqc25f)

    real(wp),parameter :: p = 0.9_wp
    real(wp),parameter :: pi = acos(-1.0_wp)

!
!           test on validity of parameters
!           ------------------------------
!

      Result = 0.0_wp
      Abserr = 0.0_wp
      Neval = 0
      Lst = 0
      Ier = 0
      IF ( (Integr/=1 .AND. Integr/=2) .OR. Epsabs<=0.0_wp .OR. &
           Limlst<3 ) Ier = 6
      IF ( Ier/=6 ) THEN
         IF ( Omega/=0.0_wp ) THEN
!
!           initializations
!           ---------------
!
            l = abs(Omega)
            dl = 2*l + 1
            cycle = dl*pi/abs(Omega)
            Ier = 0
            ktmin = 0
            Neval = 0
            numrl2 = 0
            nres = 0
            c1 = A
            c2 = cycle + A
            p1 = 1.0_wp - p
            uflow = D1MACH(1)
            eps = Epsabs
            IF ( Epsabs>uflow/p1 ) eps = Epsabs*p1
            ep = eps
            fact = 1.0_wp
            correc = 0.0_wp
            Abserr = 0.0_wp
            errsum = 0.0_wp
!
!           main do-loop
!           ------------
!
            DO Lst = 1 , Limlst
!
!           integrate over current subinterval.
!
               dla = Lst
               epsa = eps*fact
               CALL DQAWOE(F,c1,c2,Omega,Integr,epsa,0.0_wp,Limit,Lst, &
                           Maxp1,Rslst(Lst),Erlst(Lst),nev,Ierlst(Lst), &
                           last,Alist,Blist,Rlist,Elist,Iord,Nnlog, &
                           momcom,Chebmo)
               Neval = Neval + nev
               fact = fact*p
               errsum = errsum + Erlst(Lst)
               drl = 50.0_wp*abs(Rslst(Lst))
!
!           test on accuracy with partial sum
!
               IF ( (errsum+drl)<=Epsabs .AND. Lst>=6 ) GOTO 50
               correc = max(correc,Erlst(Lst))
               IF ( Ierlst(Lst)/=0 ) eps = max(ep,correc*p1)
               IF ( Ierlst(Lst)/=0 ) Ier = 7
               IF ( Ier==7 .AND. (errsum+drl)<=correc*10.0_wp .AND. &
                    Lst>5 ) GOTO 50
               numrl2 = numrl2 + 1
               IF ( Lst>1 ) THEN
                  psum(numrl2) = psum(ll) + Rslst(Lst)
                  IF ( Lst/=2 ) THEN
!
!           test on maximum number of subintervals
!
                     IF ( Lst==Limlst ) Ier = 1
!
!           perform new extrapolation
!
                     CALL DQELG(numrl2,psum,reseps,abseps,res3la,nres)
!
!           test whether extrapolated result is influenced by roundoff
!
                     ktmin = ktmin + 1
                     IF ( ktmin>=15 .AND. Abserr<=0.1e-02_wp*(errsum+drl) ) &
                          Ier = 4
                     IF ( abseps<=Abserr .OR. Lst==3 ) THEN
                        Abserr = abseps
                        Result = reseps
                        ktmin = 0
!
!           if ier is not 0, check whether direct result (partial sum)
!           or extrapolated result yields the best integral
!           approximation
!
                        IF ( (Abserr+10.0_wp*correc)<=Epsabs .OR. &
                             (Abserr<=Epsabs .AND. &
                             10.0_wp*correc>=Epsabs) ) GOTO 20
                     ENDIF
                     IF ( Ier/=0 .AND. Ier/=7 ) GOTO 20
                  ENDIF
               ELSE
                  psum(1) = Rslst(1)
               ENDIF
               ll = numrl2
               c1 = c2
               c2 = c2 + cycle
            ENDDO
!
!         set final result and error estimate
!         -----------------------------------
!
 20         Abserr = Abserr + 10.0_wp*correc
            IF ( Ier==0 ) return
            IF ( Result==0.0_wp .OR. psum(numrl2)==0.0_wp ) THEN
               IF ( Abserr>errsum ) GOTO 50
               IF ( psum(numrl2)==0.0_wp ) return
            ENDIF
            IF ( Abserr/abs(Result)<=(errsum+drl)/abs(psum(numrl2)) ) &
                 THEN
               IF ( Ier>=1 .AND. Ier/=7 ) Abserr = Abserr + drl
               return
            ENDIF
         ELSE
!
!           integration by dqagie if omega is zero
!           --------------------------------------
!
            IF ( Integr==1 ) CALL DQAGIE(F,0.0_wp,1,Epsabs,0.0_wp, &
                 Limit,Result,Abserr,Neval,Ier,Alist,Blist,Rlist,Elist, &
                 Iord,last)
            Rslst(1) = Result
            Erlst(1) = Abserr
            Ierlst(1) = Ier
            Lst = 1
            return
         ENDIF
 50      Result = psum(numrl2)
         Abserr = errsum + drl
      ENDIF
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWO(F,A,B,Omega,Integr,Epsabs,Epsrel,Result,Abserr, &
                       Neval,Ier,Leniw,Maxp1,Lenw,Last,Iwork,Work)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             integrand with oscillatory cos or sin factor,
!             clenshaw-curtis method, (end point) singularities,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i=integral of f(x)*w(x) over (a,b)
!            where w(x) = cos(omega*x)
!            or w(x) = sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of oscillatory integrals
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the function
!                     f(x).  the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            omega  - real(wp)
!                     parameter in the integrand weight function
!
!            integr - integer
!                     indicates which of the weight functions is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine will
!                     end with ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if epsabs<=0 and
!                     epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of  integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier>0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             (= leniw/2) has been achieved. one can
!                             allow more subdivisions by increasing the
!                             value of leniw (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand in
!                             order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some interior points of the
!                             integration interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table. it is presumed that
!                             the requested tolerance cannot be achieved
!                             due to roundoff in the extrapolation
!                             table, and that the returned result is
!                             the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or (integr/=1 and integr/=2),
!                             or leniw<2 or maxp1<1 or
!                             lenw<leniw*2+maxp1*25.
!                             result, abserr, neval, last are set to
!                             zero. except when leniw, maxp1 or lenw are
!                             invalid, work(limit*2+1), work(limit*3+1),
!                             iwork(1), iwork(limit+1) are set to zero,
!                             work(1) is set to a and work(limit+1) to
!                             b.
!
!         dimensioning parameters
!            leniw  - integer
!                     dimensioning parameter for iwork.
!                     leniw/2 equals the maximum number of subintervals
!                     allowed in the partition of the given integration
!                     interval (a,b), leniw>=2.
!                     if leniw<2, the routine will end with ier = 6.
!
!            maxp1  - integer
!                     gives an upper bound on the number of chebyshev
!                     moments which can be stored, i.e. for the
!                     intervals of lengths abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1>=1
!                     if maxp1<1, the routine will end with ier = 6.
!
!            lenw   - integer
!                     dimensioning parameter for work
!                     lenw must be at least leniw*2+maxp1*25.
!                     if lenw<(leniw*2+maxp1*25), the routine will
!                     end with ier = 6.
!
!            last   - integer
!                     on return, last equals the number of subintervals
!                     produced in the subdivision process, which
!                     determines the number of significant elements
!                     actually in the work arrays.
!
!         work arrays
!            iwork  - integer
!                     vector of dimension at least leniw
!                     on return, the first k elements of which contain
!                     pointers to the error estimates over the
!                     subintervals, such that work(limit*3+iwork(1)), ..
!                     work(limit*3+iwork(k)) form a decreasing
!                     sequence, with limit = lenw/2 , and k = last
!                     if last<=(limit/2+2), and k = limit+1-last
!                     otherwise.
!                     furthermore, iwork(limit+1), ..., iwork(limit+
!                     last) indicate the subdivision levels of the
!                     subintervals, such that iwork(limit+i) = l means
!                     that the subinterval numbered i is of length
!                     abs(b-a)*2**(1-l).
!
!            work   - real(wp)
!                     vector of dimension at least lenw
!                     on return
!                     work(1), ..., work(last) contain the left
!                      end points of the subintervals in the
!                      partition of (a,b),
!                     work(limit+1), ..., work(limit+last) contain
!                      the right end points,
!                     work(limit*2+1), ..., work(limit*2+last) contain
!                      the integral approximations over the
!                      subintervals,
!                     work(limit*3+1), ..., work(limit*3+last)
!                      contain the error estimates.
!                     work(limit*4+1), ..., work(limit*4+maxp1*25)
!                      provide space for storing the chebyshev moments.
!                     note that limit = lenw/2.
!

!
      real(wp) A , Abserr , B , Epsabs , Epsrel , Omega , &
                       Result , Work
      INTEGER Ier , Integr , Iwork , Last , limit , Lenw , Leniw , lvl ,&
              l1 , l2 , l3 , l4 , Maxp1 , momcom , Neval
!
      DIMENSION Iwork(Leniw) , Work(Lenw)
!
      procedure(func) :: f
!
!         check validity of leniw, maxp1 and lenw.
!

      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( Leniw>=2 .AND. Maxp1>=1 .AND. Lenw>=(Leniw*2+Maxp1*25) ) THEN
!
!         prepare call for dqawoe
!
         limit = Leniw/2
         l1 = limit + 1
         l2 = limit + l1
         l3 = limit + l2
         l4 = limit + l3
         CALL DQAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,limit,1,Maxp1, &
                     Result,Abserr,Neval,Ier,Last,Work(1),Work(l1), &
                     Work(l2),Work(l3),Iwork(1),Iwork(l1),momcom, &
                     Work(l4))
!
!         call error handler if necessary
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 0
      IF ( Ier/=0 ) CALL XERROR('abnormal return from dqawo',26,Ier,lvl)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,Limit,Icall, &
                        Maxp1,Result,Abserr,Neval,Ier,Last,Alist,Blist, &
                        Rlist,Elist,Iord,Nnlog,Momcom,Chebmo)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             integrand with oscillatory cos or sin factor,
!             clenshaw-curtis method, (end point) singularities,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral
!            i = integral of f(x)*w(x) over (a,b)
!            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of oscillatory integrals
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration
!
!            omega  - real(wp)
!                     parameter in the integrand weight function
!
!            integr - integer
!                     indicates which of the weight functions is to be
!                     used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1 and integr/=2, the routine
!                     will end with ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subdivisions
!                     in the partition of (a,b), limit>=1.
!
!            icall  - integer
!                     if dqawoe is to be used only once, icall must
!                     be set to 1.  assume that during this call, the
!                     chebyshev moments (for clenshaw-curtis integration
!                     of degree 24) have been computed for intervals of
!                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
!                     if icall>1 this means that dqawoe has been
!                     called twice or more on intervals of the same
!                     length abs(b-a). the chebyshev moments already
!                     computed are then re-used in subsequent calls.
!                     if icall<1, the routine will end with ier = 6.
!
!            maxp1  - integer
!                     gives an upper bound on the number of chebyshev
!                     moments which can be stored, i.e. for the
!                     intervals of lenghts abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1>=1.
!                     if maxp1<1, the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                   - ier>0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand, in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is presumed that the requested
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table,
!                             and that the returned result is the
!                             best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier>0.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp))
!                             or (integr/=1 and integr/=2) or
!                             icall<1 or maxp1<1.
!                             result, abserr, neval, last, rlist(1),
!                             elist(1), iord(1) and nnlog(1) are set
!                             to zero. alist(1) and blist(1) are set
!                             to a and b respectively.
!
!            last  -  integer
!                     on return, last equals the number of
!                     subintervals produces in the subdivision
!                     process, which determines the number of
!                     significant elements actually in the
!                     work arrays.
!            alist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            blist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            rlist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the error
!                     estimates over the subintervals,
!                     such that elist(iord(1)), ...,
!                     elist(iord(k)) form a decreasing sequence, with
!                     k = last if last<=(limit/2+2), and
!                     k = limit+1-last otherwise.
!
!            nnlog  - integer
!                     vector of dimension at least limit, containing the
!                     subdivision levels of the subintervals, i.e.
!                     iwork(i) = l means that the subinterval
!                     numbered i is of length abs(b-a)*2**(1-l)
!
!         on entry and return
!            momcom - integer
!                     indicating that the chebyshev moments
!                     have been computed for intervals of lengths
!                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
!                     momcom<maxp1
!
!            chebmo - real(wp)
!                     array of dimension (maxp1,25) containing the
!                     chebyshev moments
!

!
      real(wp) A , abseps , Abserr , Alist , area , area1 , &
                       area12 , area2 , a1 , a2 , B , Blist , b1 , b2 , &
                       Chebmo , correc , defab1 , defab2 , &
                       defabs , domega , dres , Elist ,&
                       epmach , Epsabs , Epsrel , erlarg , erlast , &
                       errbnd , errmax , error1 , erro12 , error2 , &
                       errsum , ertest , oflow , Omega , resabs , &
                       reseps , Result , res3la , Rlist , rlist2 , &
                       small , uflow , width
      INTEGER Icall , id , Ier , ierro , Integr , Iord , iroff1 , &
              iroff2 , iroff3 , jupbnd , k , ksgn , ktmin , Last , &
              Limit , maxerr , Maxp1 , Momcom , nev , Neval , Nnlog , &
              nres , nrmax , nrmom , numrl2
      LOGICAL extrap , noext , extall
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) , &
                Elist(Limit) , Iord(Limit) , rlist2(52) , res3la(3) , &
                Chebmo(Maxp1,25) , Nnlog(Limit)
!
      procedure(func) :: f
!
!            the dimension of rlist2 is determined by  the value of
!            limexp in subroutine dqelg (rlist2 should be of
!            dimension (limexp+2) at least).
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2
!                       containing the part of the epsilon table
!                       which is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements in rlist2. if an appropriate
!                       approximation to the compounded integral has
!                       been obtained it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation, i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true  value)
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!           oflow is the largest positive magnitude.
!

      epmach = D1MACH(4)
!
!         test on validity of parameters
!         ------------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      Iord(1) = 0
      Nnlog(1) = 0
      IF ( (Integr/=1 .AND. Integr/=2) .OR. &
           (Epsabs<=0.0_wp .AND. Epsrel<max(50.0_wp*epmach,0.5e-28_wp)) &
           .OR. Icall<1 .OR. Maxp1<1 ) Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         domega = abs(Omega)
         nrmom = 0
         IF ( Icall<=1 ) Momcom = 0
         CALL DQC25F(F,A,B,domega,Integr,nrmom,Maxp1,0,Result,Abserr, &
                     Neval,defabs,resabs,Momcom,Chebmo)
!
!           test on accuracy.
!
         dres = abs(Result)
         errbnd = max(Epsabs,Epsrel*dres)
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         IF ( Abserr<=100.0_wp*epmach*defabs .AND. Abserr>errbnd ) &
              Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( Ier/=0 .OR. Abserr<=errbnd ) THEN
            IF ( Integr==2 .AND. Omega<0.0_wp ) Result = -Result
            return
         ELSE
!
!           initializations
!           ---------------
!
            uflow = D1MACH(1)
            oflow = D1MACH(2)
            errmax = Abserr
            maxerr = 1
            area = Result
            errsum = Abserr
            Abserr = oflow
            nrmax = 1
            extrap = .FALSE.
            noext = .FALSE.
            ierro = 0
            iroff1 = 0
            iroff2 = 0
            iroff3 = 0
            ktmin = 0
            small = abs(B-A)*0.75_wp
            nres = 0
            numrl2 = 0
            extall = .FALSE.
            IF ( 0.5_wp*abs(B-A)*domega<=2.0_wp ) THEN
               numrl2 = 1
               extall = .TRUE.
               rlist2(1) = Result
            ENDIF
            IF ( 0.25_wp*abs(B-A)*domega<=2.0_wp ) extall = .TRUE.
            ksgn = -1
            IF ( dres>=(1.0_wp-50.0_wp*epmach)*defabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
            DO Last = 2 , Limit
!
!           bisect the subinterval with the nrmax-th largest
!           error estimate.
!
               nrmom = Nnlog(maxerr) + 1
               a1 = Alist(maxerr)
               b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               erlast = errmax
               CALL DQC25F(F,a1,b1,domega,Integr,nrmom,Maxp1,0,area1, &
                           error1,nev,resabs,defab1,Momcom,Chebmo)
               Neval = Neval + nev
               CALL DQC25F(F,a2,b2,domega,Integr,nrmom,Maxp1,1,area2, &
                           error2,nev,resabs,defab2,Momcom,Chebmo)
               Neval = Neval + nev
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
                  IF ( abs(Rlist(maxerr)-area12)<=0.1e-4_wp*abs(area12) &
                       .AND. erro12>=0.99_wp*errmax ) THEN
                     IF ( extrap ) iroff2 = iroff2 + 1
                     IF ( .NOT.extrap ) iroff1 = iroff1 + 1
                  ENDIF
                  IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               Nnlog(maxerr) = nrmom
               Nnlog(Last) = nrmom
               errbnd = max(Epsabs,Epsrel*abs(area))
!
!           test for roundoff error and eventually set error flag.
!
               IF ( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
               IF ( iroff2>=5 ) ierro = 3
!
!           set error flag in the case that the number of
!           subintervals equals limit.
!
               IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range.
!
               IF ( max(abs(a1),abs(b2))<=(1.0_wp+100.0_wp*epmach) &
                    *(abs(a2)+1000.0_wp*uflow) ) Ier = 4
!
!           append the newly-created intervals to the list.
!
               IF ( error2>error1 ) THEN
                  Alist(maxerr) = a2
                  Alist(Last) = a1
                  Blist(Last) = b1
                  Rlist(maxerr) = area2
                  Rlist(Last) = area1
                  Elist(maxerr) = error2
                  Elist(Last) = error1
               ELSE
                  Alist(Last) = a2
                  Blist(maxerr) = b1
                  Blist(Last) = b2
                  Elist(maxerr) = error1
                  Elist(Last) = error2
               ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with nrmax-th largest error estimate (to bisected next).
!
               CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( errsum<=errbnd ) GOTO 50
               IF ( Ier/=0 ) GOTO 40
               IF ( Last==2 .AND. extall ) THEN
                  small = small*0.5_wp
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
               ELSE
                  IF ( noext ) GOTO 20
                  IF ( extall ) THEN
                     erlarg = erlarg - erlast
                     IF ( abs(b1-a1)>small ) erlarg = erlarg + erro12
                     IF ( extrap ) GOTO 5
                  ENDIF
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
                  width = abs(Blist(maxerr)-Alist(maxerr))
                  IF ( width>small ) GOTO 20
                  IF ( extall ) THEN
                     extrap = .TRUE.
                     nrmax = 2
                  ELSE
!
!           test whether we can start with the extrapolation procedure
!           (we do this if we integrate over the next interval with
!           use of a gauss-kronrod rule - see subroutine dqc25f).
!
                     small = small*0.5_wp
                     IF ( 0.25_wp*width*domega>2.0_wp ) GOTO 20
                     extall = .TRUE.
                     GOTO 10
                  ENDIF
 5                IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors over
!           the larger intervals (erlarg) and perform extrapolation.
!
                     jupbnd = Last
                     IF ( Last>(Limit/2+2) ) jupbnd = Limit + 3 - Last
                     id = nrmax
                     DO k = id , jupbnd
                        maxerr = Iord(nrmax)
                        errmax = Elist(maxerr)
                        IF ( abs(Blist(maxerr)-Alist(maxerr))>small ) &
                             GOTO 20
                        nrmax = nrmax + 1
                     ENDDO
                  ENDIF
!
!           perform extrapolation.
!
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
                  IF ( numrl2>=3 ) THEN
                     CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
                     ktmin = ktmin + 1
                     IF ( ktmin>5 .AND. Abserr<0.1e-02_wp*errsum ) Ier = 5
                     IF ( abseps<Abserr ) THEN
                        ktmin = 0
                        Abserr = abseps
                        Result = reseps
                        correc = erlarg
                        ertest = max(Epsabs,Epsrel*abs(reseps))
! ***jump out of do-loop
                        IF ( Abserr<=ertest ) GOTO 40
                     ENDIF
!
!           prepare bisection of the smallest interval.
!
                     IF ( numrl2==1 ) noext = .TRUE.
                     IF ( Ier==5 ) GOTO 40
                  ENDIF
                  maxerr = Iord(1)
                  errmax = Elist(maxerr)
                  nrmax = 1
                  extrap = .FALSE.
                  small = small*0.5_wp
                  erlarg = errsum
                  GOTO 20
               ENDIF
 10            ertest = errbnd
               erlarg = errsum
 20         ENDDO
!
!           set the final result.
!           ---------------------
!
 40         IF ( Abserr/=oflow .AND. nres/=0 ) THEN
               IF ( Ier+ierro/=0 ) THEN
                  IF ( ierro==3 ) Abserr = Abserr + correc
                  IF ( Ier==0 ) Ier = 3
                  IF ( Result==0.0_wp .OR. area==0.0_wp ) THEN
                     IF ( Abserr>errsum ) GOTO 50
                     IF ( area==0.0_wp ) THEN
                        IF ( Ier>2 ) Ier = Ier - 1
                        IF ( Integr==2 .AND. Omega<0.0_wp ) Result = -Result
                        return
                     ENDIF
                  ELSEIF ( Abserr/abs(Result)>errsum/abs(area) ) THEN
                     GOTO 50
                  ENDIF
               ENDIF
!
!           test on divergence.
!
               IF ( ksgn/=(-1) .OR. max(abs(Result),abs(area)) &
                    >defabs*0.01_wp ) THEN
                  IF ( 0.01_wp>(Result/area) .OR. (Result/area) &
                       >100.0_wp .OR. errsum>=abs(area) ) Ier = 6
               ENDIF
               IF ( Ier>2 ) Ier = Ier - 1
               IF ( Integr==2 .AND. Omega<0.0_wp ) Result = -Result
               return
            ENDIF
         ENDIF
!
!           compute global integral sum.
!
 50      Result = 0.0_wp
         DO k = 1 , Last
            Result = Result + Rlist(k)
         ENDDO
         Abserr = errsum
         IF ( Ier>2 ) Ier = Ier - 1
         IF ( Integr==2 .AND. Omega<0.0_wp ) Result = -Result
      ENDIF
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit, &
                        Result,Abserr,Neval,Ier,Alist,Blist,Rlist,Elist,&
                        Iord,Last)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  automatic integrator, special-purpose,
!             algebraico-logarithmic end point singularities,
!             clenshaw-curtis method
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f*w over (a,b),
!            (where w shows a singular behaviour at the end points,
!            see parameter integr).
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        integration of functions having algebraico-logarithmic
!        end point singularities
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!            f      - real(wp)
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared external in the driver program.
!
!            a      - real(wp)
!                     lower limit of integration
!
!            b      - real(wp)
!                     upper limit of integration, b>a
!                     if b<=a, the routine will end with ier = 6.
!
!            alfa   - real(wp)
!                     parameter in the weight function, alfa>(-1)
!                     if alfa<=(-1), the routine will end with
!                     ier = 6.
!
!            beta   - real(wp)
!                     parameter in the weight function, beta>(-1)
!                     if beta<=(-1), the routine will end with
!                     ier = 6.
!
!            integr - integer
!                     indicates which weight function is to be used
!                     = 1  (x-a)**alfa*(b-x)**beta
!                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!                     if integr<1 or integr>4, the routine
!                     will end with ier = 6.
!
!            epsabs - real(wp)
!                     absolute accuracy requested
!            epsrel - real(wp)
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=2
!                     if limit<2, the routine will end with ier = 6.
!
!         on return
!            result - real(wp)
!                     approximation to the integral
!
!            abserr - real(wp)
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier>0 abnormal termination of the routine
!                             the estimates for the integral and error
!                             are less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement, it is advised to analyze the
!                             integrand in order to determine the
!                             integration difficulties which prevent the
!                             requested tolerance from being achieved.
!                             in case of a jump discontinuity or a local
!                             singularity of algebraico-logarithmic type
!                             at one or more interior points of the
!                             integration range, one should proceed by
!                             splitting up the interval at these
!                             points and calling the integrator on the
!                             subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b<=a or alfa<=(-1) or beta<=(-1), or
!                             integr<1 or integr>4, or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                             or limit<2.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero. alist(1)
!                             and blist(1) are set to a and b
!                             respectively.
!
!            alist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            blist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            rlist  - real(wp)
!                     vector of dimension at least limit,the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real(wp)
!                     vector of dimension at least limit, the first
!                      last  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            iord   - integer
!                     vector of dimension at least limit, the first k
!                     of which are pointers to the error
!                     estimates over the subintervals, so that
!                     elist(iord(1)), ..., elist(iord(k)) with k = last
!                     if last<=(limit/2+2), and k = limit+1-last
!                     otherwise form a decreasing sequence
!
!            last   - integer
!                     number of subintervals actually produced in
!                     the subdivision process
!

!
      real(wp) A , Abserr , Alfa , Alist , area , area1 , &
                       area12 , area2 , a1 , a2 , B , Beta , Blist , &
                       b1 , b2 , centre , &
                       Elist , epmach , Epsabs , Epsrel , errbnd , &
                       errmax , error1 , erro12 , error2 , errsum , &
                       resas1 , resas2 , Result , rg , rh , ri , rj , &
                       Rlist , uflow
      INTEGER Ier , Integr , Iord , iroff1 , iroff2 , k , Last , Limit ,&
              maxerr , nev , Neval , nrmax
!
      procedure(func) :: f
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) , &
                Elist(Limit) , Iord(Limit) , ri(25) , rj(25) , rh(25) , &
                rg(25)
!
!            list of major variables
!            -----------------------
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Ier = 6
      Neval = 0
      Last = 0
      Rlist(1) = 0.0_wp
      Elist(1) = 0.0_wp
      Iord(1) = 0
      Result = 0.0_wp
      Abserr = 0.0_wp
      IF ( .NOT.(B<=A .OR. (Epsabs==0.0_wp .AND. Epsrel<max(50.0_wp* &
           epmach,0.5e-28_wp)) .OR. Alfa<=(-1.0_wp) .OR. Beta<=(-1.0_wp) &
           .OR. Integr<1 .OR. Integr>4 .OR. Limit<2) ) THEN
         Ier = 0
!
!           compute the modified chebyshev moments.
!
         CALL DQMOMO(Alfa,Beta,ri,rj,rg,rh,Integr)
!
!           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
!
         centre = 0.5_wp*(B+A)
         CALL DQC25S(F,A,B,A,centre,Alfa,Beta,ri,rj,rg,rh,area1,error1, &
                     resas1,Integr,nev)
         Neval = nev
         CALL DQC25S(F,A,B,centre,B,Alfa,Beta,ri,rj,rg,rh,area2,error2, &
                     resas2,Integr,nev)
         Last = 2
         Neval = Neval + nev
         Result = area1 + area2
         Abserr = error1 + error2
!
!           test on accuracy.
!
         errbnd = max(Epsabs,Epsrel*abs(Result))
!
!           initialization
!           --------------
!
         IF ( error2>error1 ) THEN
            Alist(1) = centre
            Alist(2) = A
            Blist(1) = B
            Blist(2) = centre
            Rlist(1) = area2
            Rlist(2) = area1
            Elist(1) = error2
            Elist(2) = error1
         ELSE
            Alist(1) = A
            Alist(2) = centre
            Blist(1) = centre
            Blist(2) = B
            Rlist(1) = area1
            Rlist(2) = area2
            Elist(1) = error1
            Elist(2) = error2
         ENDIF
         Iord(1) = 1
         Iord(2) = 2
         IF ( Limit==2 ) Ier = 1
         IF ( Abserr>errbnd .AND. Ier/=1 ) THEN
            errmax = Elist(1)
            maxerr = 1
            nrmax = 1
            area = Result
            errsum = Abserr
            iroff1 = 0
            iroff2 = 0
!
!            main do-loop
!            ------------
!
            DO Last = 3 , Limit
!
!           bisect the subinterval with largest error estimate.
!
               a1 = Alist(maxerr)
               b1 = 0.5_wp*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
!
               CALL DQC25S(F,A,B,a1,b1,Alfa,Beta,ri,rj,rg,rh,area1, &
                           error1,resas1,Integr,nev)
               Neval = Neval + nev
               CALL DQC25S(F,A,B,a2,b2,Alfa,Beta,ri,rj,rg,rh,area2, &
                           error2,resas2,Integr,nev)
               Neval = Neval + nev
!
!           improve previous approximations integral and error
!           and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( A/=a1 .AND. B/=b2 ) THEN
                  IF ( resas1/=error1 .AND. resas2/=error2 ) THEN
!
!           test for roundoff error.
!
                     IF ( abs(Rlist(maxerr)-area12) &
                          <0.1e-4_wp*abs(area12) .AND. &
                          erro12>=0.99_wp*errmax ) iroff1 = iroff1 + 1
                     IF ( Last>10 .AND. erro12>errmax ) &
                          iroff2 = iroff2 + 1
                  ENDIF
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
!
!           test on accuracy.
!
               errbnd = max(Epsabs,Epsrel*abs(area))
               IF ( errsum>errbnd ) THEN
!
!           set error flag in the case that the number of interval
!           bisections exceeds limit.
!
                  IF ( Last==Limit ) Ier = 1
!
!
!           set error flag in the case of roundoff error.
!
                  IF ( iroff1>=6 .OR. iroff2>=20 ) Ier = 2
!
!           set error flag in the case of bad integrand behaviour
!           at interior points of integration range.
!
                  IF ( max(abs(a1),abs(b2))    &
                       <=(1.0_wp+100.0_wp*epmach)  &
                       *(abs(a2)+1000.0_wp*uflow) ) Ier = 3
               ENDIF
!
!           append the newly-created intervals to the list.
!
               IF ( error2>error1 ) THEN
                  Alist(maxerr) = a2
                  Alist(Last) = a1
                  Blist(Last) = b1
                  Rlist(maxerr) = area2
                  Rlist(Last) = area1
                  Elist(maxerr) = error2
                  Elist(Last) = error1
               ELSE
                  Alist(Last) = a2
                  Blist(maxerr) = b1
                  Blist(Last) = b2
                  Elist(maxerr) = error1
                  Elist(Last) = error2
               ENDIF
!
!           call subroutine dqpsrt to maintain the descending ordering
!           in the list of error estimates and select the subinterval
!           with largest error estimate (to be bisected next).
!
               CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0_wp
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQC25C(F,A,B,C,Result,Abserr,Krul,Neval)
      IMPLICIT NONE

!***date written   810101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  25-point clenshaw-curtis integration
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f*w over (a,b) with
!            error estimate, where w(x) = 1/(x-c)
!***description
!
!        integration rules for the computation of cauchy
!        principal value integrals
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!           f      - real(wp)
!                    function subprogram defining the integrand function
!                    f(x). the actual name for f needs to be declared
!                    external  in the driver program.
!
!           a      - real(wp)
!                    left end point of the integration interval
!
!           b      - real(wp)
!                    right end point of the integration interval, b>a
!
!           c      - real(wp)
!                    parameter in the weight function
!
!           result - real(wp)
!                    approximation to the integral
!                    result is computed by using a generalized
!                    clenshaw-curtis method if c lies within ten percent
!                    of the integration interval. in the other case the
!                    15-point kronrod rule obtained by optimal addition
!                    of abscissae to the 7-point gauss rule, is applied.
!
!           abserr - real(wp)
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed abs(i-result)
!
!           krul   - integer
!                    key which is decreased by 1 if the 15-point
!                    gauss-kronrod scheme has been used
!
!           neval  - integer
!                    number of integrand evaluations
!
!.......................................................................

!
      real(wp) A , Abserr , ak22 , amom0 , amom1 , amom2 , B , &
                       C , cc , centr , cheb12 , cheb24 , &
                       fval , hlgth , p2 , p3 , p4 , &
                       resabs , resasc , Result , res12 , res24 , u , x
      INTEGER i , isym , k , kp , Krul , Neval
!
      DIMENSION x(11) , fval(25) , cheb12(13) , cheb24(25)
!
      procedure(func) :: f
!
!           the vector x contains the values cos(k*pi/24),
!           k = 1, ..., 11, to be used for the chebyshev series
!           expansion of f
!
      DATA x(1)/0.991444861373810411144557526928563_wp/
      DATA x(2)/0.965925826289068286749743199728897_wp/
      DATA x(3)/0.923879532511286756128183189396788_wp/
      DATA x(4)/0.866025403784438646763723170752936_wp/
      DATA x(5)/0.793353340291235164579776961501299_wp/
      DATA x(6)/0.707106781186547524400844362104849_wp/
      DATA x(7)/0.608761429008720639416097542898164_wp/
      DATA x(8)/0.500000000000000000000000000000000_wp/
      DATA x(9)/0.382683432365089771728459984030399_wp/
      DATA x(10)/0.258819045102520762348898837624048_wp/
      DATA x(11)/0.130526192220051591548406227895489_wp/
!
!           list of major variables
!           ----------------------
!           fval   - value of the function f at the points
!                    cos(k*pi/24),  k = 0, ..., 24
!           cheb12 - chebyshev series expansion coefficients,
!                    for the function f, of degree 12
!           cheb24 - chebyshev series expansion coefficients,
!                    for the function f, of degree 24
!           res12  - approximation to the integral corresponding
!                    to the use of cheb12
!           res24  - approximation to the integral corresponding
!                    to the use of cheb24
!           dqwgtc - external function subprogram defining
!                    the weight function
!           hlgth  - half-length of the interval
!           centr  - mid point of the interval
!
!
!           check the position of c.
!

      cc = (2.0_wp*C-B-A)/(B-A)
      IF ( abs(cc)<1.1_wp ) THEN
!
!           use the generalized clenshaw-curtis method.
!
         hlgth = 0.5_wp*(B-A)
         centr = 0.5_wp*(B+A)
         Neval = 25
         fval(1) = 0.5_wp*F(hlgth+centr)
         fval(13) = F(centr)
         fval(25) = 0.5_wp*F(centr-hlgth)
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)
            fval(isym) = F(centr-u)
         ENDDO
!
!           compute the chebyshev series expansion.
!
         CALL DQCHEB(x,fval,cheb12,cheb24)
!
!           the modified chebyshev moments are computed by forward
!           recursion, using amom0 and amom1 as starting values.
!
         amom0 = log(abs((1.0_wp-cc)/(1.0_wp+cc)))
         amom1 = 2.0_wp + cc*amom0
         res12 = cheb12(1)*amom0 + cheb12(2)*amom1
         res24 = cheb24(1)*amom0 + cheb24(2)*amom1
         DO k = 3 , 13
            amom2 = 2.0_wp*cc*amom1 - amom0
            ak22 = (k-2)*(k-2)
            IF ( (k/2)*2==k ) amom2 = amom2 - 4.0_wp/(ak22-1.0_wp)
            res12 = res12 + cheb12(k)*amom2
            res24 = res24 + cheb24(k)*amom2
            amom0 = amom1
            amom1 = amom2
         ENDDO
         DO k = 14 , 25
            amom2 = 2.0_wp*cc*amom1 - amom0
            ak22 = (k-2)*(k-2)
            IF ( (k/2)*2==k ) amom2 = amom2 - 4.0_wp/(ak22-1.0_wp)
            res24 = res24 + cheb24(k)*amom2
            amom0 = amom1
            amom1 = amom2
         ENDDO
         Result = res24
         Abserr = abs(res24-res12)
      ELSE
!
!           apply the 15-point gauss-kronrod scheme.
!
         Krul = Krul - 1
         CALL DQK15W(F,DQWGTC,C,p2,p3,p4,kp,A,B,Result,Abserr,resabs, &
                     resasc)
         Neval = 15
         IF ( resasc==Abserr ) Krul = Krul + 1
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQC25F(F,A,B,Omega,Integr,Nrmom,Maxp1,Ksave,Result, &
                        Abserr,Neval,Resabs,Resasc,Momcom,Chebmo)
      IMPLICIT NONE

!***date written   810101   (yymmdd)
!***revision date  211011   (yymmdd)
!***keywords  integration rules for functions with cos or sin
!             factor, clenshaw-curtis, gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute the integral i=integral of f(x) over (a,b)
!            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
!            compute j = integral of abs(f) over (a,b). for small value
!            of omega or small intervals (a,b) the 15-point gauss-kronro
!            rule is used. otherwise a generalized clenshaw-curtis
!            method is used.
!***description
!
!        integration rules for functions with cos or sin factor
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!         on entry
!           f      - real(wp)
!                    function subprogram defining the integrand
!                    function f(x). the actual name for f needs to
!                    be declared external in the calling program.
!
!           a      - real(wp)
!                    lower limit of integration
!
!           b      - real(wp)
!                    upper limit of integration
!
!           omega  - real(wp)
!                    parameter in the weight function
!
!           integr - integer
!                    indicates which weight function is to be used
!                       integr = 1   w(x) = cos(omega*x)
!                       integr = 2   w(x) = sin(omega*x)
!
!           nrmom  - integer
!                    the length of interval (a,b) is equal to the length
!                    of the original integration interval divided by
!                    2**nrmom (we suppose that the routine is used in an
!                    adaptive integration process, otherwise set
!                    nrmom = 0). nrmom must be zero at the first call.
!
!           maxp1  - integer
!                    gives an upper bound on the number of chebyshev
!                    moments which can be stored, i.e. for the
!                    intervals of lengths abs(bb-aa)*2**(-l),
!                    l = 0,1,2, ..., maxp1-2.
!
!           ksave  - integer
!                    key which is one when the moments for the
!                    current interval have been computed
!
!         on return
!           result - real(wp)
!                    approximation to the integral i
!
!           abserr - real(wp)
!                    estimate of the modulus of the absolute
!                    error, which should equal or exceed abs(i-result)
!
!           neval  - integer
!                    number of integrand evaluations
!
!           resabs - real(wp)
!                    approximation to the integral j
!
!           resasc - real(wp)
!                    approximation to the integral of abs(f-i/(b-a))
!
!         on entry and return
!           momcom - integer
!                    for each interval length we need to compute the
!                    chebyshev moments. momcom counts the number of
!                    intervals for which these moments have already been
!                    computed. if nrmom<momcom or ksave = 1, the
!                    chebyshev moments for the interval (a,b) have
!                    already been computed and stored, otherwise we
!                    compute them and we increase momcom.
!
!           chebmo - real(wp)
!                    array of dimension at least (maxp1,25) containing
!                    the modified chebyshev moments for the first momcom
!                    momcom interval lengths
!
! ......................................................................

!
      real(wp) A , Abserr , ac , an , an2 , as , asap , ass , &
                       B , centr , Chebmo , cheb12 , cheb24 , conc , &
                       cons , cospar , d ,&
                       d1 , d2 , estc , ests , fval , &
                       hlgth , oflow , Omega , parint , par2 , par22 , &
                       p2 , p3 , p4 , Resabs , Resasc , resc12 , &
                       resc24 , ress12 , ress24 , Result , sinpar , v , &
                       x
      INTEGER i , iers , Integr , isym , j , k , Ksave , m , Momcom , &
              Neval , Maxp1 , noequ , noeq1 , Nrmom
!
      DIMENSION Chebmo(Maxp1,25) , cheb12(13) , cheb24(25) , d(25) , &
                d1(25) , d2(25) , fval(25) , v(28) , x(11)
!
      procedure(func) :: f
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ...,11, to be used for the chebyshev expansion of f
!
      DATA x(1)/0.991444861373810411144557526928563_wp/
      DATA x(2)/0.965925826289068286749743199728897_wp/
      DATA x(3)/0.923879532511286756128183189396788_wp/
      DATA x(4)/0.866025403784438646763723170752936_wp/
      DATA x(5)/0.793353340291235164579776961501299_wp/
      DATA x(6)/0.707106781186547524400844362104849_wp/
      DATA x(7)/0.608761429008720639416097542898164_wp/
      DATA x(8)/0.500000000000000000000000000000000_wp/
      DATA x(9)/0.382683432365089771728459984030399_wp/
      DATA x(10)/0.258819045102520762348898837624048_wp/
      DATA x(11)/0.130526192220051591548406227895489_wp/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fval   - value of the function f at the points
!                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
!           cheb12 - coefficients of the chebyshev series expansion
!                    of degree 12, for the function f, in the
!                    interval (a,b)
!           cheb24 - coefficients of the chebyshev series expansion
!                    of degree 24, for the function f, in the
!                    interval (a,b)
!           resc12 - approximation to the integral of
!                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
!                    over (-1,+1), using the chebyshev series
!                    expansion of degree 12
!           resc24 - approximation to the same integral, using the
!                    chebyshev series expansion of degree 24
!           ress12 - the analogue of resc12 for the sine
!           ress24 - the analogue of resc24 for the sine
!
!
!           machine dependent constant
!           --------------------------
!
!           oflow is the largest positive magnitude.
!

      oflow = D1MACH(2)
!
      centr = 0.5_wp*(B+A)
      hlgth = 0.5_wp*(B-A)
      parint = Omega*hlgth
!
!           compute the integral using the 15-point gauss-kronrod
!           formula if the value of the parameter in the integrand
!           is small.
!
      IF ( abs(parint)>2.0_wp ) THEN
!
!           compute the integral using the generalized clenshaw-
!           curtis method.
!
         conc = hlgth*cos(centr*Omega)
         cons = hlgth*sin(centr*Omega)
         Resasc = oflow
         Neval = 25
!
!           check whether the chebyshev moments for this interval
!           have already been computed.
!
         IF ( Nrmom>=Momcom .AND. Ksave/=1 ) THEN
!
!           compute a new set of chebyshev moments.
!
            m = Momcom + 1
            par2 = parint*parint
            par22 = par2 + 2.0_wp
            sinpar = sin(parint)
            cospar = cos(parint)
!
!           compute the chebyshev moments with respect to cosine.
!
            v(1) = 2.0_wp*sinpar/parint
            v(2) = (8.0_wp*cospar+(par2+par2-8.0_wp)*sinpar/parint) &
                   /par2
            v(3) = (32.0_wp*(par2-12.0_wp)*cospar+(2.0_wp*((par2- &
                   80.0_wp)*par2+192.0_wp)*sinpar)/parint)/(par2*par2)
            ac = 8.0_wp*cospar
            as = 24.0_wp*parint*sinpar
            IF ( abs(parint)>24.0_wp ) THEN
!
!           compute the chebyshev moments by means of forward
!           recursion.
!
               an = 4.0_wp
               DO i = 4 , 13
                  an2 = an*an
                  v(i) = ((an2-4.0_wp)*(2.0_wp*(par22-an2-an2)*v(i-1)-&
                         ac)+as-par2*(an+1.0_wp)*(an+2.0_wp)*v(i-2)) &
                         /(par2*(an-1.0_wp)*(an-2.0_wp))
                  an = an + 2.0_wp
               ENDDO
            ELSE
!
!           compute the chebyshev moments as the solutions of a
!           boundary value problem with 1 initial value (v(3)) and 1
!           end value (computed using an asymptotic formula).
!
               noequ = 25
               noeq1 = noequ - 1
               an = 6.0_wp
               DO k = 1 , noeq1
                  an2 = an*an
                  d(k) = -2.0_wp*(an2-4.0_wp)*(par22-an2-an2)
                  d2(k) = (an-1.0_wp)*(an-2.0_wp)*par2
                  d1(k+1) = (an+3.0_wp)*(an+4.0_wp)*par2
                  v(k+3) = as - (an2-4.0_wp)*ac
                  an = an + 2.0_wp
               ENDDO
               an2 = an*an
               d(noequ) = -2.0_wp*(an2-4.0_wp)*(par22-an2-an2)
               v(noequ+3) = as - (an2-4.0_wp)*ac
               v(4) = v(4) - 56.0_wp*par2*v(3)
               ass = parint*sinpar
               asap = (((((210.0_wp*par2-1.0_wp)*cospar-(105.0_wp* &
                      par2-63.0_wp)*ass)/an2-(1.0_wp-15.0_wp*par2) &
                      *cospar+15.0_wp*ass)/an2-cospar+3.0_wp*ass) &
                      /an2-cospar)/an2
               v(noequ+3) = v(noequ+3) - 2.0_wp*asap*par2*(an-1.0_wp) &
                            *(an-2.0_wp)
!
!           solve the tridiagonal system by means of gaussian
!           elimination with partial pivoting.
!
!***        call to dgtsl must be replaced by call to
!***        real(wp) version of linpack routine sgtsl
!
               CALL DGTSL(noequ,d1,d,d2,v(4),iers)
            ENDIF
            DO j = 1 , 13
               Chebmo(m,2*j-1) = v(j)
            ENDDO
!
!           compute the chebyshev moments with respect to sine.
!
            v(1) = 2.0_wp*(sinpar-parint*cospar)/par2
            v(2) = (18.0_wp-48.0_wp/par2)*sinpar/par2 + &
                   (-2.0_wp+48.0_wp/par2)*cospar/parint
            ac = -24.0_wp*parint*cospar
            as = -8.0_wp*sinpar
            IF ( abs(parint)>24.0_wp ) THEN
!
!           compute the chebyshev moments by means of forward recursion.
!
               an = 3.0_wp
               DO i = 3 , 12
                  an2 = an*an
                  v(i) = ((an2-4.0_wp)*(2.0_wp*(par22-an2-an2)*v(i-1)+&
                         as)+ac-par2*(an+1.0_wp)*(an+2.0_wp)*v(i-2)) &
                         /(par2*(an-1.0_wp)*(an-2.0_wp))
                  an = an + 2.0_wp
               ENDDO
            ELSE
!
!           compute the chebyshev moments as the solutions of a boundary
!           value problem with 1 initial value (v(2)) and 1 end value
!           (computed using an asymptotic formula).
!
               an = 0.5D+01
               DO k = 1 , noeq1
                  an2 = an*an
                  d(k) = -2.0_wp*(an2-4.0_wp)*(par22-an2-an2)
                  d2(k) = (an-1.0_wp)*(an-2.0_wp)*par2
                  d1(k+1) = (an+3.0_wp)*(an+4.0_wp)*par2
                  v(k+2) = ac + (an2-4.0_wp)*as
                  an = an + 2.0_wp
               ENDDO
               an2 = an*an
               d(noequ) = -2.0_wp*(an2-4.0_wp)*(par22-an2-an2)
               v(noequ+2) = ac + (an2-4.0_wp)*as
               v(3) = v(3) - 0.42D+02*par2*v(2)
               ass = parint*cospar
               asap = (((((105.0_wp*par2-63.0_wp)*ass+(210.0_wp*par2-&
                      1.0_wp)*sinpar)/an2+(15.0_wp*par2-1.0_wp) &
                      *sinpar-15.0_wp*ass)/an2-3.0_wp*ass-sinpar) &
                      /an2-sinpar)/an2
               v(noequ+2) = v(noequ+2) - 2.0_wp*asap*par2*(an-1.0_wp) &
                            *(an-2.0_wp)
!
!           solve the tridiagonal system by means of gaussian
!           elimination with partial pivoting.
!
!***        call to dgtsl must be replaced by call to
!***        real(wp) version of linpack routine sgtsl
!
               CALL DGTSL(noequ,d1,d,d2,v(3),iers)
            ENDIF
            DO j = 1 , 12
               Chebmo(m,2*j) = v(j)
            ENDDO
         ENDIF
         IF ( Nrmom<Momcom ) m = Nrmom + 1
         IF ( Momcom<(Maxp1-1) .AND. Nrmom>=Momcom ) Momcom = Momcom + 1
!
!           compute the coefficients of the chebyshev expansions
!           of degrees 12 and 24 of the function f.
!
         fval(1) = 0.5_wp*F(centr+hlgth)
         fval(13) = F(centr)
         fval(25) = 0.5_wp*F(centr-hlgth)
         DO i = 2 , 12
            isym = 26 - i
            fval(i) = F(hlgth*x(i-1)+centr)
            fval(isym) = F(centr-hlgth*x(i-1))
         ENDDO
         CALL DQCHEB(x,fval,cheb12,cheb24)
!
!           compute the integral and error estimates.
!
         resc12 = cheb12(13)*Chebmo(m,13)
         ress12 = 0.0_wp
         k = 11
         DO j = 1 , 6
            resc12 = resc12 + cheb12(k)*Chebmo(m,k)
            ress12 = ress12 + cheb12(k+1)*Chebmo(m,k+1)
            k = k - 2
         ENDDO
         resc24 = cheb24(25)*Chebmo(m,25)
         ress24 = 0.0_wp
         Resabs = abs(cheb24(25))
         k = 23
         DO j = 1 , 12
            resc24 = resc24 + cheb24(k)*Chebmo(m,k)
            ress24 = ress24 + cheb24(k+1)*Chebmo(m,k+1)
            Resabs = Resabs + abs(cheb24(k)) + abs(cheb24(k+1))
            k = k - 2
         ENDDO
         estc = abs(resc24-resc12)
         ests = abs(ress24-ress12)
         Resabs = Resabs*abs(hlgth)
         IF ( Integr==2 ) THEN
            Result = conc*ress24 + cons*resc24
            Abserr = abs(conc*ests) + abs(cons*estc)
         ELSE
            Result = conc*resc24 - cons*ress24
            Abserr = abs(conc*estc) + abs(cons*ests)
         ENDIF
      ELSE
         CALL DQK15W(F,DQWGTF,Omega,p2,p3,p4,Integr,A,B,Result,Abserr, &
                     Resabs,Resasc)
         Neval = 15
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQC25S(F,A,B,Bl,Br,Alfa,Beta,Ri,Rj,Rg,Rh,Result,Abserr,&
                        Resasc,Integr,Nev)
      IMPLICIT NONE

!***date written   810101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  25-point clenshaw-curtis integration
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f*w over (bl,br), with error
!            estimate, where the weight function w has a singular
!            behaviour of algebraico-logarithmic type at the points
!            a and/or b. (bl,br) is a part of (a,b).
!***description
!
!        integration rules for integrands having algebraico-logarithmic
!        end point singularities
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!           f      - real(wp)
!                    function subprogram defining the integrand
!                    f(x). the actual name for f needs to be declared
!                    external  in the driver program.
!
!           a      - real(wp)
!                    left end point of the original interval
!
!           b      - real(wp)
!                    right end point of the original interval, b>a
!
!           bl     - real(wp)
!                    lower limit of integration, bl>=a
!
!           br     - real(wp)
!                    upper limit of integration, br<=b
!
!           alfa   - real(wp)
!                    parameter in the weight function
!
!           beta   - real(wp)
!                    parameter in the weight function
!
!           ri,rj,rg,rh - real(wp)
!                    modified chebyshev moments for the application
!                    of the generalized clenshaw-curtis
!                    method (computed in subroutine dqmomo)
!
!           result - real(wp)
!                    approximation to the integral
!                    result is computed by using a generalized
!                    clenshaw-curtis method if b1 = a or br = b.
!                    in all other cases the 15-point kronrod
!                    rule is applied, obtained by optimal addition of
!                    abscissae to the 7-point gauss rule.
!
!           abserr - real(wp)
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed abs(i-result)
!
!           resasc - real(wp)
!                    approximation to the integral of abs(f*w-i/(b-a))
!
!           integr - integer
!                    which determines the weight function
!                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
!                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
!                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
!                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
!                                 log(b-x)
!
!           nev    - integer
!                    number of integrand evaluations

!
      real(wp) A , Abserr , Alfa , B , Beta , Bl , Br , centr , &
                       cheb12 , cheb24 , dc , factor ,&
                       fix , fval , hlgth , resabs , Resasc , Result , &
                       res12 , res24 , Rg , Rh , Ri , Rj , u , &
                       x
      INTEGER i , Integr , isym , Nev
!
      DIMENSION cheb12(13) , cheb24(25) , fval(25) , Rg(25) , Rh(25) , &
                Ri(25) , Rj(25) , x(11)
!
      procedure(func) :: f
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ..., 11, to be used for the computation of the
!           chebyshev series expansion of f.
!
      DATA x(1)/0.991444861373810411144557526928563_wp/
      DATA x(2)/0.965925826289068286749743199728897_wp/
      DATA x(3)/0.923879532511286756128183189396788_wp/
      DATA x(4)/0.866025403784438646763723170752936_wp/
      DATA x(5)/0.793353340291235164579776961501299_wp/
      DATA x(6)/0.707106781186547524400844362104849_wp/
      DATA x(7)/0.608761429008720639416097542898164_wp/
      DATA x(8)/0.500000000000000000000000000000000_wp/
      DATA x(9)/0.382683432365089771728459984030399_wp/
      DATA x(10)/0.258819045102520762348898837624048_wp/
      DATA x(11)/0.130526192220051591548406227895489_wp/
!
!           list of major variables
!           -----------------------
!
!           fval   - value of the function f at the points
!                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
!                    k = 0, ..., 24
!           cheb12 - coefficients of the chebyshev series expansion
!                    of degree 12, for the function f, in the
!                    interval (bl,br)
!           cheb24 - coefficients of the chebyshev series expansion
!                    of degree 24, for the function f, in the
!                    interval (bl,br)
!           res12  - approximation to the integral obtained from cheb12
!           res24  - approximation to the integral obtained from cheb24
!           dqwgts - external function subprogram defining
!                    the four possible weight functions
!           hlgth  - half-length of the interval (bl,br)
!           centr  - mid point of the interval (bl,br)
!

      Nev = 25
      IF ( Bl==A .AND. (Alfa/=0.0_wp .OR. Integr==2 .OR. Integr==4) ) &
           THEN
!
!           this part of the program is executed only if a = bl.
!           ----------------------------------------------------
!
!           compute the chebyshev series expansion of the
!           following function
!           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
!                  *f(0.5*(br-a)*x+0.5*(br+a))
!
         hlgth = 0.5_wp*(Br-Bl)
         centr = 0.5_wp*(Br+Bl)
         fix = B - centr
         fval(1) = 0.5_wp*F(hlgth+centr)*(fix-hlgth)**Beta
         fval(13) = F(centr)*(fix**Beta)
         fval(25) = 0.5_wp*F(centr-hlgth)*(fix+hlgth)**Beta
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)*(fix-u)**Beta
            fval(isym) = F(centr-u)*(fix+u)**Beta
         ENDDO
         factor = hlgth**(Alfa+1.0_wp)
         Result = 0.0_wp
         Abserr = 0.0_wp
         res12 = 0.0_wp
         res24 = 0.0_wp
         IF ( Integr>2 ) THEN
!
!           compute the chebyshev series expansion of the
!           following function
!           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
!
            fval(1) = fval(1)*log(fix-hlgth)
            fval(13) = fval(13)*log(fix)
            fval(25) = fval(25)*log(fix+hlgth)
            DO i = 2 , 12
               u = hlgth*x(i-1)
               isym = 26 - i
               fval(i) = fval(i)*log(fix-u)
               fval(isym) = fval(isym)*log(fix+u)
            ENDDO
            CALL DQCHEB(x,fval,cheb12,cheb24)
!
!           integr = 3  (or 4)
!
            DO i = 1 , 13
               res12 = res12 + cheb12(i)*Ri(i)
               res24 = res24 + cheb24(i)*Ri(i)
            ENDDO
            DO i = 14 , 25
               res24 = res24 + cheb24(i)*Ri(i)
            ENDDO
            IF ( Integr/=3 ) THEN
!
!           integr = 4
!
               dc = log(Br-Bl)
               Result = res24*dc
               Abserr = abs((res24-res12)*dc)
               res12 = 0.0_wp
               res24 = 0.0_wp
               DO i = 1 , 13
                  res12 = res12 + cheb12(i)*Rg(i)
                  res24 = res24 + cheb24(i)*Rg(i)
               ENDDO
               DO i = 14 , 25
                  res24 = res24 + cheb24(i)*Rg(i)
               ENDDO
            ENDIF
         ELSE
            CALL DQCHEB(x,fval,cheb12,cheb24)
!
!           integr = 1  (or 2)
!
            DO i = 1 , 13
               res12 = res12 + cheb12(i)*Ri(i)
               res24 = res24 + cheb24(i)*Ri(i)
            ENDDO
            DO i = 14 , 25
               res24 = res24 + cheb24(i)*Ri(i)
            ENDDO
            IF ( Integr/=1 ) THEN
!
!           integr = 2
!
               dc = log(Br-Bl)
               Result = res24*dc
               Abserr = abs((res24-res12)*dc)
               res12 = 0.0_wp
               res24 = 0.0_wp
               DO i = 1 , 13
                  res12 = res12 + cheb12(i)*Rg(i)
                  res24 = res12 + cheb24(i)*Rg(i)
               ENDDO
               DO i = 14 , 25
                  res24 = res24 + cheb24(i)*Rg(i)
               ENDDO
            ENDIF
         ENDIF
         Result = (Result+res24)*factor
         Abserr = (Abserr+abs(res24-res12))*factor
      ELSEIF ( Br==B .AND. (Beta/=0.0_wp .OR. Integr==3 .OR. Integr==4) ) THEN
!
!           this part of the program is executed only if b = br.
!           ----------------------------------------------------
!
!           compute the chebyshev series expansion of the
!           following function
!           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
!                *f(0.5*(b-bl)*x+0.5*(b+bl))
!
         hlgth = 0.5_wp*(Br-Bl)
         centr = 0.5_wp*(Br+Bl)
         fix = centr - A
         fval(1) = 0.5_wp*F(hlgth+centr)*(fix+hlgth)**Alfa
         fval(13) = F(centr)*(fix**Alfa)
         fval(25) = 0.5_wp*F(centr-hlgth)*(fix-hlgth)**Alfa
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)*(fix+u)**Alfa
            fval(isym) = F(centr-u)*(fix-u)**Alfa
         ENDDO
         factor = hlgth**(Beta+1.0_wp)
         Result = 0.0_wp
         Abserr = 0.0_wp
         res12 = 0.0_wp
         res24 = 0.0_wp
         IF ( Integr==2 .OR. Integr==4 ) THEN
!
!           compute the chebyshev series expansion of the
!           following function
!           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
!
            fval(1) = fval(1)*log(hlgth+fix)
            fval(13) = fval(13)*log(fix)
            fval(25) = fval(25)*log(fix-hlgth)
            DO i = 2 , 12
               u = hlgth*x(i-1)
               isym = 26 - i
               fval(i) = fval(i)*log(u+fix)
               fval(isym) = fval(isym)*log(fix-u)
            ENDDO
            CALL DQCHEB(x,fval,cheb12,cheb24)
!
!           integr = 2  (or 4)
!
            DO i = 1 , 13
               res12 = res12 + cheb12(i)*Rj(i)
               res24 = res24 + cheb24(i)*Rj(i)
            ENDDO
            DO i = 14 , 25
               res24 = res24 + cheb24(i)*Rj(i)
            ENDDO
            IF ( Integr/=2 ) THEN
               dc = log(Br-Bl)
               Result = res24*dc
               Abserr = abs((res24-res12)*dc)
               res12 = 0.0_wp
               res24 = 0.0_wp
!
!           integr = 4
!
               DO i = 1 , 13
                  res12 = res12 + cheb12(i)*Rh(i)
                  res24 = res24 + cheb24(i)*Rh(i)
               ENDDO
               DO i = 14 , 25
                  res24 = res24 + cheb24(i)*Rh(i)
               ENDDO
            ENDIF
         ELSE
!
!           integr = 1  (or 3)
!
            CALL DQCHEB(x,fval,cheb12,cheb24)
            DO i = 1 , 13
               res12 = res12 + cheb12(i)*Rj(i)
               res24 = res24 + cheb24(i)*Rj(i)
            ENDDO
            DO i = 14 , 25
               res24 = res24 + cheb24(i)*Rj(i)
            ENDDO
            IF ( Integr/=1 ) THEN
!
!           integr = 3
!
               dc = log(Br-Bl)
               Result = res24*dc
               Abserr = abs((res24-res12)*dc)
               res12 = 0.0_wp
               res24 = 0.0_wp
               DO i = 1 , 13
                  res12 = res12 + cheb12(i)*Rh(i)
                  res24 = res24 + cheb24(i)*Rh(i)
               ENDDO
               DO i = 14 , 25
                  res24 = res24 + cheb24(i)*Rh(i)
               ENDDO
            ENDIF
         ENDIF
         Result = (Result+res24)*factor
         Abserr = (Abserr+abs(res24-res12))*factor
      ELSE
!
!           if a>bl and b<br, apply the 15-point gauss-kronrod
!           scheme.
!
!
         CALL DQK15W(F,DQWGTS,A,B,Alfa,Beta,Integr,Bl,Br,Result,Abserr, &
                     resabs,Resasc)
         Nev = 15
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQCHEB(X,Fval,Cheb12,Cheb24)
      IMPLICIT NONE

!### See also
!  *  dqc25c,dqc25f,dqc25s
!***revision date  830518   (yymmdd)
!***keywords  chebyshev series expansion, fast fourier transform
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this routine computes the chebyshev series expansion
!            of degrees 12 and 24 of a function using a
!            fast fourier transform method
!            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
!            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
!            where t(k,x) is the chebyshev polynomial of degree k.
!***description
!
!        chebyshev series expansion
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!          on entry
!           x      - real(wp)
!                    vector of dimension 11 containing the
!                    values cos(k*pi/24), k = 1, ..., 11
!
!           fval   - real(wp)
!                    vector of dimension 25 containing the
!                    function values at the points
!                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
!                    where (a,b) is the approximation interval.
!                    fval(1) and fval(25) are divided by two
!                    (these values are destroyed at output).
!
!          on return
!           cheb12 - real(wp)
!                    vector of dimension 13 containing the
!                    chebyshev coefficients for degree 12
!
!           cheb24 - real(wp)
!                    vector of dimension 25 containing the
!                    chebyshev coefficients for degree 24
!
!
      real(wp) alam , alam1 , alam2 , Cheb12 , Cheb24 , Fval , &
                       part1 , part2 , part3 , v , X
      INTEGER i , j
!
      DIMENSION Cheb12(13) , Cheb24(25) , Fval(25) , v(12) , X(11)
!

      DO i = 1 , 12
         j = 26 - i
         v(i) = Fval(i) - Fval(j)
         Fval(i) = Fval(i) + Fval(j)
      ENDDO
      alam1 = v(1) - v(9)
      alam2 = X(6)*(v(3)-v(7)-v(11))
      Cheb12(4) = alam1 + alam2
      Cheb12(10) = alam1 - alam2
      alam1 = v(2) - v(8) - v(10)
      alam2 = v(4) - v(6) - v(12)
      alam = X(3)*alam1 + X(9)*alam2
      Cheb24(4) = Cheb12(4) + alam
      Cheb24(22) = Cheb12(4) - alam
      alam = X(9)*alam1 - X(3)*alam2
      Cheb24(10) = Cheb12(10) + alam
      Cheb24(16) = Cheb12(10) - alam
      part1 = X(4)*v(5)
      part2 = X(8)*v(9)
      part3 = X(6)*v(7)
      alam1 = v(1) + part1 + part2
      alam2 = X(2)*v(3) + part3 + X(10)*v(11)
      Cheb12(2) = alam1 + alam2
      Cheb12(12) = alam1 - alam2
      alam = X(1)*v(2) + X(3)*v(4) + X(5)*v(6) + X(7)*v(8) + X(9)*v(10) &
             + X(11)*v(12)
      Cheb24(2) = Cheb12(2) + alam
      Cheb24(24) = Cheb12(2) - alam
      alam = X(11)*v(2) - X(9)*v(4) + X(7)*v(6) - X(5)*v(8) + X(3)*v(10)&
             - X(1)*v(12)
      Cheb24(12) = Cheb12(12) + alam
      Cheb24(14) = Cheb12(12) - alam
      alam1 = v(1) - part1 + part2
      alam2 = X(10)*v(3) - part3 + X(2)*v(11)
      Cheb12(6) = alam1 + alam2
      Cheb12(8) = alam1 - alam2
      alam = X(5)*v(2) - X(9)*v(4) - X(1)*v(6) - X(11)*v(8) + X(3)*v(10)&
             + X(7)*v(12)
      Cheb24(6) = Cheb12(6) + alam
      Cheb24(20) = Cheb12(6) - alam
      alam = X(7)*v(2) - X(3)*v(4) - X(11)*v(6) + X(1)*v(8) - X(9)*v(10)&
             - X(5)*v(12)
      Cheb24(8) = Cheb12(8) + alam
      Cheb24(18) = Cheb12(8) - alam
      DO i = 1 , 6
         j = 14 - i
         v(i) = Fval(i) - Fval(j)
         Fval(i) = Fval(i) + Fval(j)
      ENDDO
      alam1 = v(1) + X(8)*v(5)
      alam2 = X(4)*v(3)
      Cheb12(3) = alam1 + alam2
      Cheb12(11) = alam1 - alam2
      Cheb12(7) = v(1) - v(5)
      alam = X(2)*v(2) + X(6)*v(4) + X(10)*v(6)
      Cheb24(3) = Cheb12(3) + alam
      Cheb24(23) = Cheb12(3) - alam
      alam = X(6)*(v(2)-v(4)-v(6))
      Cheb24(7) = Cheb12(7) + alam
      Cheb24(19) = Cheb12(7) - alam
      alam = X(10)*v(2) - X(6)*v(4) + X(2)*v(6)
      Cheb24(11) = Cheb12(11) + alam
      Cheb24(15) = Cheb12(11) - alam
      DO i = 1 , 3
         j = 8 - i
         v(i) = Fval(i) - Fval(j)
         Fval(i) = Fval(i) + Fval(j)
      ENDDO
      Cheb12(5) = v(1) + X(8)*v(3)
      Cheb12(9) = Fval(1) - X(8)*Fval(3)
      alam = X(4)*v(2)
      Cheb24(5) = Cheb12(5) + alam
      Cheb24(21) = Cheb12(5) - alam
      alam = X(8)*Fval(2) - Fval(4)
      Cheb24(9) = Cheb12(9) + alam
      Cheb24(17) = Cheb12(9) - alam
      Cheb12(1) = Fval(1) + Fval(3)
      alam = Fval(2) + Fval(4)
      Cheb24(1) = Cheb12(1) + alam
      Cheb24(25) = Cheb12(1) - alam
      Cheb12(13) = v(1) - v(3)
      Cheb24(13) = Cheb12(13)
      alam = 1.0_wp/6.0_wp
      DO i = 2 , 12
         Cheb12(i) = Cheb12(i)*alam
      ENDDO
      alam = 0.5_wp*alam
      Cheb12(1) = Cheb12(1)*alam
      Cheb12(13) = Cheb12(13)*alam
      DO i = 2 , 24
         Cheb24(i) = Cheb24(i)*alam
      ENDDO
      Cheb24(1) = 0.5_wp*alam*Cheb24(1)
      Cheb24(25) = 0.5_wp*alam*Cheb24(25)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQELG(N,Epstab,Result,Abserr,Res3la,Nres)
      IMPLICIT NONE

!### See also
!  *  dqagie,dqagoe,dqagpe,dqagse
!***revision date  830518   (yymmdd)
!***keywords  epsilon algorithm, convergence acceleration,
!             extrapolation
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  the routine determines the limit of a given sequence of
!            approximations, by means of the epsilon algorithm of
!            p.wynn. an estimate of the absolute error is also given.
!            the condensed epsilon table is computed. only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***description
!
!           epsilon algorithm
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!              n      - integer
!                       epstab(n) contains the new element in the
!                       first column of the epsilon table.
!
!              epstab - real(wp)
!                       vector of dimension 52 containing the elements
!                       of the two lower diagonals of the triangular
!                       epsilon table. the elements are numbered
!                       starting at the right-hand corner of the
!                       triangle.
!
!              result - real(wp)
!                       resulting approximation to the integral
!
!              abserr - real(wp)
!                       estimate of the absolute error computed from
!                       result and the 3 previous results
!
!              res3la - real(wp)
!                       vector of dimension 3 containing the last 3
!                       results
!
!              nres   - integer
!                       number of calls to the routine
!                       (should be zero at first call)
!
!
      real(wp) Abserr , delta1 , delta2 , delta3 , &
                       epmach , epsinf , Epstab , &
                       error , err1 , err2 , err3 , e0 , e1 , e1abs , &
                       e2 , e3 , oflow , res , Result , Res3la , ss , &
                       tol1 , tol2 , tol3
      INTEGER i , ib , ib2 , ie , indx , k1 , k2 , k3 , limexp , N , &
              newelm , Nres , num
      DIMENSION Epstab(52) , Res3la(3)
!
!           list of major variables
!           -----------------------
!
!           e0     - the 4 elements on which the computation of a new
!           e1       element in the epsilon table is based
!           e2
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           oflow is the largest positive magnitude.
!           limexp is the maximum number of elements the epsilon
!           table can contain. if this number is reached, the upper
!           diagonal of the epsilon table is deleted.
!

      epmach = D1MACH(4)
      oflow = D1MACH(2)
      Nres = Nres + 1
      Abserr = oflow
      Result = Epstab(N)
      IF ( N>=3 ) THEN
         limexp = 50
         Epstab(N+2) = Epstab(N)
         newelm = (N-1)/2
         Epstab(N) = oflow
         num = N
         k1 = N
         DO i = 1 , newelm
            k2 = k1 - 1
            k3 = k1 - 2
            res = Epstab(k1+2)
            e0 = Epstab(k3)
            e1 = Epstab(k2)
            e2 = res
            e1abs = abs(e1)
            delta2 = e2 - e1
            err2 = abs(delta2)
            tol2 = max(abs(e2),e1abs)*epmach
            delta3 = e1 - e0
            err3 = abs(delta3)
            tol3 = max(e1abs,abs(e0))*epmach
            IF ( err2>tol2 .OR. err3>tol3 ) THEN
               e3 = Epstab(k1)
               Epstab(k1) = e1
               delta1 = e1 - e3
               err1 = abs(delta1)
               tol1 = max(e1abs,abs(e3))*epmach
!
!           if two elements are very close to each other, omit
!           a part of the table by adjusting the value of n
!
               IF ( err1>tol1 .AND. err2>tol2 .AND. err3>tol3 ) THEN
                  ss = 1.0_wp/delta1 + 1.0_wp/delta2 - 1.0_wp/delta3
                  epsinf = abs(ss*e1)
!
!           test to detect irregular behaviour in the table, and
!           eventually omit a part of the table adjusting the value
!           of n.
!
                  IF ( epsinf>0.1D-03 ) THEN
!
!           compute a new element and eventually adjust
!           the value of result.
!
                     res = e1 + 1.0_wp/ss
                     Epstab(k1) = res
                     k1 = k1 - 2
                     error = err2 + abs(res-e2) + err3
                     IF ( error<=Abserr ) THEN
                        Abserr = error
                        Result = res
                     ENDIF
                     GOTO 50
                  ENDIF
               ENDIF
               N = i + i - 1
! ***jump out of do-loop
               GOTO 100
            ELSE
!
!           if e0, e1 and e2 are equal to within machine
!           accuracy, convergence is assumed.
!           result = e2
!           abserr = abs(e1-e0)+abs(e2-e1)
!
               Result = res
               Abserr = err2 + err3
! ***jump out of do-loop
               GOTO 200
            ENDIF
 50      ENDDO
!
!           shift the table.
!
 100     IF ( N==limexp ) N = 2*(limexp/2) - 1
         ib = 1
         IF ( (num/2)*2==num ) ib = 2
         ie = newelm + 1
         DO i = 1 , ie
            ib2 = ib + 2
            Epstab(ib) = Epstab(ib2)
            ib = ib2
         ENDDO
         IF ( num/=N ) THEN
            indx = num - N + 1
            DO i = 1 , N
               Epstab(i) = Epstab(indx)
               indx = indx + 1
            ENDDO
         ENDIF
         IF ( Nres>=4 ) THEN
!
!           compute error estimate
!
            Abserr = abs(Result-Res3la(3)) + abs(Result-Res3la(2)) &
                     + abs(Result-Res3la(1))
            Res3la(1) = Res3la(2)
            Res3la(2) = Res3la(3)
            Res3la(3) = Result
         ELSE
            Res3la(Nres) = Result
            Abserr = oflow
         ENDIF
      ENDIF
 200  Abserr = max(Abserr,0.5D+01*epmach*abs(Result))
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK15(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the calling program.
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!

!
      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(7) , fv2(7) , wg(4) , wgk(8) , xgk(8)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.129484966168869693270611432679082_wp/
      DATA wg(2)/0.279705391489276667901467771423780_wp/
      DATA wg(3)/0.381830050505118944950369775488975_wp/
      DATA wg(4)/0.417959183673469387755102040816327_wp/
!
      DATA xgk(1)/0.991455371120812639206854697526329_wp/
      DATA xgk(2)/0.949107912342758524526189684047851_wp/
      DATA xgk(3)/0.864864423359769072789712788640926_wp/
      DATA xgk(4)/0.741531185599394439863864773280788_wp/
      DATA xgk(5)/0.586087235467691130294144838258730_wp/
      DATA xgk(6)/0.405845151377397166906606412076961_wp/
      DATA xgk(7)/0.207784955007898467600689403773245_wp/
      DATA xgk(8)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.022935322010529224963732008058970_wp/
      DATA wgk(2)/0.063092092629978553290700663189204_wp/
      DATA wgk(3)/0.104790010322250183839876322541518_wp/
      DATA wgk(4)/0.140653259715525918745189590510238_wp/
      DATA wgk(5)/0.169004726639267902826583426598550_wp/
      DATA wgk(6)/0.190350578064785409913256402421014_wp/
      DATA wgk(7)/0.204432940075298892414161999234649_wp/
      DATA wgk(8)/0.209482141084727828012999174891714_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      Resabs = abs(resk)
      DO j = 1 , 3
         jtw = j*2
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 4
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(8)*abs(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j) &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK15I(F,Boun,Inf,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  15-point transformed gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the original (infinite integration range is mapped
!            onto the interval (0,1) and (a,b) is a part of (0,1).
!            it is the purpose to compute
!            i = integral of transformed integrand over (a,b),
!            j = integral of abs(transformed integrand) over (a,b).
!***description
!
!           integration rule
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       fuction subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the calling program.
!
!              boun   - real(wp)
!                       finite bound of original integration
!                       range (set to zero if inf = +2)
!
!              inf    - integer
!                       if inf = -1, the original interval is
!                                   (-infinity,bound),
!                       if inf = +1, the original interval is
!                                   (bound,+infinity),
!                       if inf = +2, the original interval is
!                                   (-infinity,+infinity) and
!                       the integral is computed as the sum of two
!                       integrals, one over (-infinity,0) and one over
!                       (0,+infinity).
!
!              a      - real(wp)
!                       lower limit for integration over subrange
!                       of (0,1)
!
!              b      - real(wp)
!                       upper limit for integration over subrange
!                       of (0,1)
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule(resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule(resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integral of
!                       abs((transformed integrand)-i/(b-a)) over (a,b)
!

!
      real(wp) A , absc , absc1 , absc2 , Abserr , B , Boun , &
                       centr , dinf , min , &
                       epmach , fc , fsum , fval1 , fval2 , fv1 , &
                       fv2 , hlgth , Resabs , Resasc , resg , resk , &
                       reskh , Result , tabsc1 , tabsc2 , uflow , wg , &
                       wgk , xgk
      INTEGER Inf , j
      procedure(func) :: f
!
      DIMENSION fv1(7) , fv2(7) , xgk(8) , wgk(8) , wg(8)
!
!           the abscissae and weights are supplied for the interval
!           (-1,1).  because of symmetry only the positive abscissae and
!           their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
      DATA wg(1)/0.0_wp/
      DATA wg(2)/0.129484966168869693270611432679082_wp/
      DATA wg(3)/0.0_wp/
      DATA wg(4)/0.279705391489276667901467771423780_wp/
      DATA wg(5)/0.0_wp/
      DATA wg(6)/0.381830050505118944950369775488975_wp/
      DATA wg(7)/0.0_wp/
      DATA wg(8)/0.417959183673469387755102040816327_wp/
!
      DATA xgk(1)/0.991455371120812639206854697526329_wp/
      DATA xgk(2)/0.949107912342758524526189684047851_wp/
      DATA xgk(3)/0.864864423359769072789712788640926_wp/
      DATA xgk(4)/0.741531185599394439863864773280788_wp/
      DATA xgk(5)/0.586087235467691130294144838258730_wp/
      DATA xgk(6)/0.405845151377397166906606412076961_wp/
      DATA xgk(7)/0.207784955007898467600689403773245_wp/
      DATA xgk(8)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.022935322010529224963732008058970_wp/
      DATA wgk(2)/0.063092092629978553290700663189204_wp/
      DATA wgk(3)/0.104790010322250183839876322541518_wp/
      DATA wgk(4)/0.140653259715525918745189590510238_wp/
      DATA wgk(5)/0.169004726639267902826583426598550_wp/
      DATA wgk(6)/0.190350578064785409913256402421014_wp/
      DATA wgk(7)/0.204432940075298892414161999234649_wp/
      DATA wgk(8)/0.209482141084727828012999174891714_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
      dinf = MIN0(1,Inf)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      tabsc1 = Boun + dinf*(1.0_wp-centr)/centr
      fval1 = F(tabsc1)
      IF ( Inf==2 ) fval1 = fval1 + F(-tabsc1)
      fc = (fval1/centr)/centr
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the error.
!
      resg = wg(8)*fc
      resk = wgk(8)*fc
      Resabs = abs(resk)
      DO j = 1 , 7
         absc = hlgth*xgk(j)
         absc1 = centr - absc
         absc2 = centr + absc
         tabsc1 = Boun + dinf*(1.0_wp-absc1)/absc1
         tabsc2 = Boun + dinf*(1.0_wp-absc2)/absc2
         fval1 = F(tabsc1)
         fval2 = F(tabsc2)
         IF ( Inf==2 ) fval1 = fval1 + F(-tabsc1)
         IF ( Inf==2 ) fval2 = fval2 + F(-tabsc2)
         fval1 = (fval1/absc1)/absc1
         fval2 = (fval2/absc2)/absc2
         fv1(j) = fval1
         fv2(j) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(j)*fsum
         Resabs = Resabs + wgk(j)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(8)*abs(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resasc = Resasc*hlgth
      Resabs = Resabs*hlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp )    &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK15W(F,W,P1,P2,P3,P4,Kp,A,B,Result,Abserr,Resabs, &
                        Resasc)
      IMPLICIT NONE

!***date written   810101   (yymmdd)
!***revision date  830518   (mmddyy)
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f*w over (a,b), with error
!                           estimate
!                       j = integral of abs(f*w) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!             on entry
!              f      - real(wp)
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the driver program.
!
!              w      - real(wp)
!                       function subprogram defining the integrand
!                       weight function w(x). the actual name for w
!                       needs to be declared external in the
!                       calling program.
!
!              p1, p2, p3, p4 - real(wp)
!                       parameters in the weight function
!
!              kp     - integer
!                       key for indicating the type of weight function
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule (resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral of abs(f)
!
!              resasc - real(wp)
!                       approximation to the integral of abs(f-i/(b-a))
!
!

!
      real(wp) A , absc , absc1 , absc2 , Abserr , B , centr , &
                       abs , dhlgth , min , epmach ,&
                       fc , fsum , fval1 , fval2 , fv1 , fv2 , &
                       hlgth , P1 , P2 , P3 , P4 , Resabs , Resasc , &
                       resg , resk , reskh , Result , uflow , W , wg , &
                       wgk , xgk
      INTEGER j , jtw , jtwm1 , Kp
      procedure(func) :: f
      external :: W
!
      DIMENSION fv1(7) , fv2(7) , xgk(8) , wgk(8) , wg(4)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point gauss-kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ... abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point gauss-kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
           , xgk(8)/0.9914553711208126_wp , 0.9491079123427585_wp , &
           0.8648644233597691_wp , 0.7415311855993944_wp , &
           0.5860872354676911_wp , 0.4058451513773972_wp , &
           0.2077849550078985_wp , 0.0000000000000000_wp/
!
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
           , wgk(8)/0.2293532201052922D-01 , 0.6309209262997855D-01 , &
           0.1047900103222502_wp , 0.1406532597155259_wp , &
           0.1690047266392679_wp , 0.1903505780647854_wp , &
           0.2044329400752989_wp , 0.2094821410847278_wp/
!
      DATA wg(1) , wg(2) , wg(3) , wg(4)/0.1294849661688697_wp , &
           0.2797053914892767_wp , 0.3818300505051189_wp , &
           0.4179591836734694_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f*w over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 15-point kronrod approximation to the
!           integral, and estimate the error.
!
      fc = F(centr)*W(centr,P1,P2,P3,P4,Kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      Resabs = abs(resk)
      DO j = 1 , 3
         jtw = j*2
         absc = hlgth*xgk(jtw)
         absc1 = centr - absc
         absc2 = centr + absc
         fval1 = F(absc1)*W(absc1,P1,P2,P3,P4,Kp)
         fval2 = F(absc2)*W(absc2,P1,P2,P3,P4,Kp)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 4
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         absc1 = centr - absc
         absc2 = centr + absc
         fval1 = F(absc1)*W(absc1,P1,P2,P3,P4,Kp)
         fval2 = F(absc2)*W(absc2,P1,P2,P3,P4,Kp)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(8)*abs(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK21(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  21-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the driver program.
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 21-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 10-point gauss rule (resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!

!
      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(10) , fv2(10) , wg(5) , wgk(11) , xgk(11)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point gauss rule
!
!           wgk    - weights of the 21-point kronrod rule
!
!           wg     - weights of the 10-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.066671344308688137593568809893332_wp/
      DATA wg(2)/0.149451349150580593145776339657697_wp/
      DATA wg(3)/0.219086362515982043995534934228163_wp/
      DATA wg(4)/0.269266719309996355091226921569469_wp/
      DATA wg(5)/0.295524224714752870173892994651338_wp/
!
      DATA xgk(1)/0.995657163025808080735527280689003_wp/
      DATA xgk(2)/0.973906528517171720077964012084452_wp/
      DATA xgk(3)/0.930157491355708226001207180059508_wp/
      DATA xgk(4)/0.865063366688984510732096688423493_wp/
      DATA xgk(5)/0.780817726586416897063717578345042_wp/
      DATA xgk(6)/0.679409568299024406234327365114874_wp/
      DATA xgk(7)/0.562757134668604683339000099272694_wp/
      DATA xgk(8)/0.433395394129247190799265943165784_wp/
      DATA xgk(9)/0.294392862701460198131126603103866_wp/
      DATA xgk(10)/0.148874338981631210884826001129720_wp/
      DATA xgk(11)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.011694638867371874278064396062192_wp/
      DATA wgk(2)/0.032558162307964727478818972459390_wp/
      DATA wgk(3)/0.054755896574351996031381300244580_wp/
      DATA wgk(4)/0.075039674810919952767043140916190_wp/
      DATA wgk(5)/0.093125454583697605535065465083366_wp/
      DATA wgk(6)/0.109387158802297641899210590325805_wp/
      DATA wgk(7)/0.123491976262065851077958109831074_wp/
      DATA wgk(8)/0.134709217311473325928054001771707_wp/
      DATA wgk(9)/0.142775938577060080797094273138717_wp/
      DATA wgk(10)/0.147739104901338491374841515972068_wp/
      DATA wgk(11)/0.149445554002916905664936468389821_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point gauss formula
!           resk   - result of the 21-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 21-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      resg = 0.0_wp
      fc = F(centr)
      resk = wgk(11)*fc
      Resabs = abs(resk)
      DO j = 1 , 5
         jtw = 2*j
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 5
         jtwm1 = 2*j - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(11)*abs(fc-reskh)
      DO j = 1 , 10
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK31(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  31-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the calling program.
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 31-point
!                       gauss-kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point gauss
!                       rule (resg).
!
!              abserr - double precison
!                       estimate of the modulus of the modulus,
!                       which should not exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!

      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(15) , fv2(15) , xgk(16) , wgk(16) , wg(8)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 31-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 15-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 15-point gauss rule
!
!           wgk    - weights of the 31-point kronrod rule
!
!           wg     - weights of the 15-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.030753241996117268354628393577204_wp/
      DATA wg(2)/0.070366047488108124709267416450667_wp/
      DATA wg(3)/0.107159220467171935011869546685869_wp/
      DATA wg(4)/0.139570677926154314447804794511028_wp/
      DATA wg(5)/0.166269205816993933553200860481209_wp/
      DATA wg(6)/0.186161000015562211026800561866423_wp/
      DATA wg(7)/0.198431485327111576456118326443839_wp/
      DATA wg(8)/0.202578241925561272880620199967519_wp/
!
      DATA xgk(1)/0.998002298693397060285172840152271_wp/
      DATA xgk(2)/0.987992518020485428489565718586613_wp/
      DATA xgk(3)/0.967739075679139134257347978784337_wp/
      DATA xgk(4)/0.937273392400705904307758947710209_wp/
      DATA xgk(5)/0.897264532344081900882509656454496_wp/
      DATA xgk(6)/0.848206583410427216200648320774217_wp/
      DATA xgk(7)/0.790418501442465932967649294817947_wp/
      DATA xgk(8)/0.724417731360170047416186054613938_wp/
      DATA xgk(9)/0.650996741297416970533735895313275_wp/
      DATA xgk(10)/0.570972172608538847537226737253911_wp/
      DATA xgk(11)/0.485081863640239680693655740232351_wp/
      DATA xgk(12)/0.394151347077563369897207370981045_wp/
      DATA xgk(13)/0.299180007153168812166780024266389_wp/
      DATA xgk(14)/0.201194093997434522300628303394596_wp/
      DATA xgk(15)/0.101142066918717499027074231447392_wp/
      DATA xgk(16)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.005377479872923348987792051430128_wp/
      DATA wgk(2)/0.015007947329316122538374763075807_wp/
      DATA wgk(3)/0.025460847326715320186874001019653_wp/
      DATA wgk(4)/0.035346360791375846222037948478360_wp/
      DATA wgk(5)/0.044589751324764876608227299373280_wp/
      DATA wgk(6)/0.053481524690928087265343147239430_wp/
      DATA wgk(7)/0.062009567800670640285139230960803_wp/
      DATA wgk(8)/0.069854121318728258709520077099147_wp/
      DATA wgk(9)/0.076849680757720378894432777482659_wp/
      DATA wgk(10)/0.083080502823133021038289247286104_wp/
      DATA wgk(11)/0.088564443056211770647275443693774_wp/
      DATA wgk(12)/0.093126598170825321225486872747346_wp/
      DATA wgk(13)/0.096642726983623678505179907627589_wp/
      DATA wgk(14)/0.099173598721791959332393173484603_wp/
      DATA wgk(15)/0.100769845523875595044946662617570_wp/
      DATA wgk(16)/0.101330007014791549017374792767493_wp/
!
!
!           list of major variables
!           -----------------------
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 15-point gauss formula
!           resk   - result of the 31-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 31-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      Resabs = abs(resk)
      DO j = 1 , 7
         jtw = j*2
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 8
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(16)*abs(fc-reskh)
      DO j = 1 , 15
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK41(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  41-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the calling program.
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 41-point
!                       gauss-kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point gauss
!                       rule (resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integal of abs(f-i/(b-a))
!                       over (a,b)
!

!
      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(20) , fv2(20) , xgk(21) , wgk(21) , wg(10)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 41-point gauss-kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 20-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 20-point gauss rule
!
!           wgk    - weights of the 41-point gauss-kronrod rule
!
!           wg     - weights of the 20-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.017614007139152118311861962351853_wp/
      DATA wg(2)/0.040601429800386941331039952274932_wp/
      DATA wg(3)/0.062672048334109063569506535187042_wp/
      DATA wg(4)/0.083276741576704748724758143222046_wp/
      DATA wg(5)/0.101930119817240435036750135480350_wp/
      DATA wg(6)/0.118194531961518417312377377711382_wp/
      DATA wg(7)/0.131688638449176626898494499748163_wp/
      DATA wg(8)/0.142096109318382051329298325067165_wp/
      DATA wg(9)/0.149172986472603746787828737001969_wp/
      DATA wg(10)/0.152753387130725850698084331955098_wp/
!
      DATA xgk(1)/0.998859031588277663838315576545863_wp/
      DATA xgk(2)/0.993128599185094924786122388471320_wp/
      DATA xgk(3)/0.981507877450250259193342994720217_wp/
      DATA xgk(4)/0.963971927277913791267666131197277_wp/
      DATA xgk(5)/0.940822633831754753519982722212443_wp/
      DATA xgk(6)/0.912234428251325905867752441203298_wp/
      DATA xgk(7)/0.878276811252281976077442995113078_wp/
      DATA xgk(8)/0.839116971822218823394529061701521_wp/
      DATA xgk(9)/0.795041428837551198350638833272788_wp/
      DATA xgk(10)/0.746331906460150792614305070355642_wp/
      DATA xgk(11)/0.693237656334751384805490711845932_wp/
      DATA xgk(12)/0.636053680726515025452836696226286_wp/
      DATA xgk(13)/0.575140446819710315342946036586425_wp/
      DATA xgk(14)/0.510867001950827098004364050955251_wp/
      DATA xgk(15)/0.443593175238725103199992213492640_wp/
      DATA xgk(16)/0.373706088715419560672548177024927_wp/
      DATA xgk(17)/0.301627868114913004320555356858592_wp/
      DATA xgk(18)/0.227785851141645078080496195368575_wp/
      DATA xgk(19)/0.152605465240922675505220241022678_wp/
      DATA xgk(20)/0.076526521133497333754640409398838_wp/
      DATA xgk(21)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.003073583718520531501218293246031_wp/
      DATA wgk(2)/0.008600269855642942198661787950102_wp/
      DATA wgk(3)/0.014626169256971252983787960308868_wp/
      DATA wgk(4)/0.020388373461266523598010231432755_wp/
      DATA wgk(5)/0.025882133604951158834505067096153_wp/
      DATA wgk(6)/0.031287306777032798958543119323801_wp/
      DATA wgk(7)/0.036600169758200798030557240707211_wp/
      DATA wgk(8)/0.041668873327973686263788305936895_wp/
      DATA wgk(9)/0.046434821867497674720231880926108_wp/
      DATA wgk(10)/0.050944573923728691932707670050345_wp/
      DATA wgk(11)/0.055195105348285994744832372419777_wp/
      DATA wgk(12)/0.059111400880639572374967220648594_wp/
      DATA wgk(13)/0.062653237554781168025870122174255_wp/
      DATA wgk(14)/0.065834597133618422111563556969398_wp/
      DATA wgk(15)/0.068648672928521619345623411885368_wp/
      DATA wgk(16)/0.071054423553444068305790361723210_wp/
      DATA wgk(17)/0.073030690332786667495189417658913_wp/
      DATA wgk(18)/0.074582875400499188986581418362488_wp/
      DATA wgk(19)/0.075704497684556674659542775376617_wp/
      DATA wgk(20)/0.076377867672080736705502835038061_wp/
      DATA wgk(21)/0.076600711917999656445049901530102_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 20-point gauss formula
!           resk   - result of the 41-point kronrod formula
!           reskh  - approximation to mean value of f over (a,b), i.e.
!                    to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 41-point gauss-kronrod approximation to
!           the integral, and estimate the absolute error.
!
      resg = 0.0_wp
      fc = F(centr)
      resk = wgk(21)*fc
      Resabs = abs(resk)
      DO j = 1 , 10
         jtw = j*2
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 10
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(21)*abs(fc-reskh)
      DO j = 1 , 20
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0._wp )  &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK51(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  51-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real(wp) version
!
!           parameters
!            on entry
!              f      - real(wp)
!                       function subroutine defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared external in the calling program.
!
!              a      - real(wp)
!                       lower limit of integration
!
!              b      - real(wp)
!                       upper limit of integration
!
!            on return
!              result - real(wp)
!                       approximation to the integral i
!                       result is computed by applying the 51-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point gauss rule (resg).
!
!              abserr - real(wp)
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real(wp)
!                       approximation to the integral j
!
!              resasc - real(wp)
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!

!
      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(25) , fv2(25) , xgk(26) , wgk(26) , wg(13)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 51-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 25-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 25-point gauss rule
!
!           wgk    - weights of the 51-point kronrod rule
!
!           wg     - weights of the 25-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.011393798501026287947902964113235_wp/
      DATA wg(2)/0.026354986615032137261901815295299_wp/
      DATA wg(3)/0.040939156701306312655623487711646_wp/
      DATA wg(4)/0.054904695975835191925936891540473_wp/
      DATA wg(5)/0.068038333812356917207187185656708_wp/
      DATA wg(6)/0.080140700335001018013234959669111_wp/
      DATA wg(7)/0.091028261982963649811497220702892_wp/
      DATA wg(8)/0.100535949067050644202206890392686_wp/
      DATA wg(9)/0.108519624474263653116093957050117_wp/
      DATA wg(10)/0.114858259145711648339325545869556_wp/
      DATA wg(11)/0.119455763535784772228178126512901_wp/
      DATA wg(12)/0.122242442990310041688959518945852_wp/
      DATA wg(13)/0.123176053726715451203902873079050_wp/
!
      DATA xgk(1)/0.999262104992609834193457486540341_wp/
      DATA xgk(2)/0.995556969790498097908784946893902_wp/
      DATA xgk(3)/0.988035794534077247637331014577406_wp/
      DATA xgk(4)/0.976663921459517511498315386479594_wp/
      DATA xgk(5)/0.961614986425842512418130033660167_wp/
      DATA xgk(6)/0.942974571228974339414011169658471_wp/
      DATA xgk(7)/0.920747115281701561746346084546331_wp/
      DATA xgk(8)/0.894991997878275368851042006782805_wp/
      DATA xgk(9)/0.865847065293275595448996969588340_wp/
      DATA xgk(10)/0.833442628760834001421021108693570_wp/
      DATA xgk(11)/0.797873797998500059410410904994307_wp/
      DATA xgk(12)/0.759259263037357630577282865204361_wp/
      DATA xgk(13)/0.717766406813084388186654079773298_wp/
      DATA xgk(14)/0.673566368473468364485120633247622_wp/
      DATA xgk(15)/0.626810099010317412788122681624518_wp/
      DATA xgk(16)/0.577662930241222967723689841612654_wp/
      DATA xgk(17)/0.526325284334719182599623778158010_wp/
      DATA xgk(18)/0.473002731445714960522182115009192_wp/
      DATA xgk(19)/0.417885382193037748851814394594572_wp/
      DATA xgk(20)/0.361172305809387837735821730127641_wp/
      DATA xgk(21)/0.303089538931107830167478909980339_wp/
      DATA xgk(22)/0.243866883720988432045190362797452_wp/
      DATA xgk(23)/0.183718939421048892015969888759528_wp/
      DATA xgk(24)/0.122864692610710396387359818808037_wp/
      DATA xgk(25)/0.061544483005685078886546392366797_wp/
      DATA xgk(26)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.001987383892330315926507851882843_wp/
      DATA wgk(2)/0.005561932135356713758040236901066_wp/
      DATA wgk(3)/0.009473973386174151607207710523655_wp/
      DATA wgk(4)/0.013236229195571674813656405846976_wp/
      DATA wgk(5)/0.016847817709128298231516667536336_wp/
      DATA wgk(6)/0.020435371145882835456568292235939_wp/
      DATA wgk(7)/0.024009945606953216220092489164881_wp/
      DATA wgk(8)/0.027475317587851737802948455517811_wp/
      DATA wgk(9)/0.030792300167387488891109020215229_wp/
      DATA wgk(10)/0.034002130274329337836748795229551_wp/
      DATA wgk(11)/0.037116271483415543560330625367620_wp/
      DATA wgk(12)/0.040083825504032382074839284467076_wp/
      DATA wgk(13)/0.042872845020170049476895792439495_wp/
      DATA wgk(14)/0.045502913049921788909870584752660_wp/
      DATA wgk(15)/0.047982537138836713906392255756915_wp/
      DATA wgk(16)/0.050277679080715671963325259433440_wp/
      DATA wgk(17)/0.052362885806407475864366712137873_wp/
      DATA wgk(18)/0.054251129888545490144543370459876_wp/
      DATA wgk(19)/0.055950811220412317308240686382747_wp/
      DATA wgk(20)/0.057437116361567832853582693939506_wp/
      DATA wgk(21)/0.058689680022394207961974175856788_wp/
      DATA wgk(22)/0.059720340324174059979099291932562_wp/
      DATA wgk(23)/0.060539455376045862945360267517565_wp/
      DATA wgk(24)/0.061128509717053048305859030416293_wp/
      DATA wgk(25)/0.061471189871425316661544131965264_wp/
!       note: wgk (26) was calculated from the values of wgk(1..25)
      DATA wgk(26)/0.061580818067832935078759824240066_wp/
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 25-point gauss formula
!           resk   - result of the 51-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(A+B)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 51-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      Resabs = abs(resk)
      DO j = 1 , 12
         jtw = j*2
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 13
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(26)*abs(fc-reskh)
      DO j = 1 , 25
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQK61(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  61-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!        integration rule
!        standard fortran subroutine
!        real(wp) version
!
!
!        parameters
!         on entry
!           f      - real(wp)
!                    function subprogram defining the integrand
!                    function f(x). the actual name for f needs to be
!                    declared external in the calling program.
!
!           a      - real(wp)
!                    lower limit of integration
!
!           b      - real(wp)
!                    upper limit of integration
!
!         on return
!           result - real(wp)
!                    approximation to the integral i
!                    result is computed by applying the 61-point
!                    kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point gauss rule (resg).
!
!           abserr - real(wp)
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed abs(i-result)
!
!           resabs - real(wp)
!                    approximation to the integral j
!
!           resasc - real(wp)
!                    approximation to the integral of abs(f-i/(b-a))
!
!

!
      real(wp) A , dabsc , Abserr , B , centr , dhlgth , &
                       min , epmach , fc , fsum ,&
                       fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , &
                       Resasc , resg , resk , reskh , Result , uflow , &
                       wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      procedure(func) :: f
!
      DIMENSION fv1(30) , fv2(30) , xgk(31) , wgk(31) , wg(15)
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point gauss rule
!
!           wgk   - weights of the 61-point kronrod rule
!
!           wg    - weigths of the 30-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
      DATA wg(1)/0.007968192496166605615465883474674_wp/
      DATA wg(2)/0.018466468311090959142302131912047_wp/
      DATA wg(3)/0.028784707883323369349719179611292_wp/
      DATA wg(4)/0.038799192569627049596801936446348_wp/
      DATA wg(5)/0.048402672830594052902938140422808_wp/
      DATA wg(6)/0.057493156217619066481721689402056_wp/
      DATA wg(7)/0.065974229882180495128128515115962_wp/
      DATA wg(8)/0.073755974737705206268243850022191_wp/
      DATA wg(9)/0.080755895229420215354694938460530_wp/
      DATA wg(10)/0.086899787201082979802387530715126_wp/
      DATA wg(11)/0.092122522237786128717632707087619_wp/
      DATA wg(12)/0.096368737174644259639468626351810_wp/
      DATA wg(13)/0.099593420586795267062780282103569_wp/
      DATA wg(14)/0.101762389748405504596428952168554_wp/
      DATA wg(15)/0.102852652893558840341285636705415_wp/
!
      DATA xgk(1)/0.999484410050490637571325895705811_wp/
      DATA xgk(2)/0.996893484074649540271630050918695_wp/
      DATA xgk(3)/0.991630996870404594858628366109486_wp/
      DATA xgk(4)/0.983668123279747209970032581605663_wp/
      DATA xgk(5)/0.973116322501126268374693868423707_wp/
      DATA xgk(6)/0.960021864968307512216871025581798_wp/
      DATA xgk(7)/0.944374444748559979415831324037439_wp/
      DATA xgk(8)/0.926200047429274325879324277080474_wp/
      DATA xgk(9)/0.905573307699907798546522558925958_wp/
      DATA xgk(10)/0.882560535792052681543116462530226_wp/
      DATA xgk(11)/0.857205233546061098958658510658944_wp/
      DATA xgk(12)/0.829565762382768397442898119732502_wp/
      DATA xgk(13)/0.799727835821839083013668942322683_wp/
      DATA xgk(14)/0.767777432104826194917977340974503_wp/
      DATA xgk(15)/0.733790062453226804726171131369528_wp/
      DATA xgk(16)/0.697850494793315796932292388026640_wp/
      DATA xgk(17)/0.660061064126626961370053668149271_wp/
      DATA xgk(18)/0.620526182989242861140477556431189_wp/
      DATA xgk(19)/0.579345235826361691756024932172540_wp/
      DATA xgk(20)/0.536624148142019899264169793311073_wp/
      DATA xgk(21)/0.492480467861778574993693061207709_wp/
      DATA xgk(22)/0.447033769538089176780609900322854_wp/
      DATA xgk(23)/0.400401254830394392535476211542661_wp/
      DATA xgk(24)/0.352704725530878113471037207089374_wp/
      DATA xgk(25)/0.304073202273625077372677107199257_wp/
      DATA xgk(26)/0.254636926167889846439805129817805_wp/
      DATA xgk(27)/0.204525116682309891438957671002025_wp/
      DATA xgk(28)/0.153869913608583546963794672743256_wp/
      DATA xgk(29)/0.102806937966737030147096751318001_wp/
      DATA xgk(30)/0.051471842555317695833025213166723_wp/
      DATA xgk(31)/0.000000000000000000000000000000000_wp/
!
      DATA wgk(1)/0.001389013698677007624551591226760_wp/
      DATA wgk(2)/0.003890461127099884051267201844516_wp/
      DATA wgk(3)/0.006630703915931292173319826369750_wp/
      DATA wgk(4)/0.009273279659517763428441146892024_wp/
      DATA wgk(5)/0.011823015253496341742232898853251_wp/
      DATA wgk(6)/0.014369729507045804812451432443580_wp/
      DATA wgk(7)/0.016920889189053272627572289420322_wp/
      DATA wgk(8)/0.019414141193942381173408951050128_wp/
      DATA wgk(9)/0.021828035821609192297167485738339_wp/
      DATA wgk(10)/0.024191162078080601365686370725232_wp/
      DATA wgk(11)/0.026509954882333101610601709335075_wp/
      DATA wgk(12)/0.028754048765041292843978785354334_wp/
      DATA wgk(13)/0.030907257562387762472884252943092_wp/
      DATA wgk(14)/0.032981447057483726031814191016854_wp/
      DATA wgk(15)/0.034979338028060024137499670731468_wp/
      DATA wgk(16)/0.036882364651821229223911065617136_wp/
      DATA wgk(17)/0.038678945624727592950348651532281_wp/
      DATA wgk(18)/0.040374538951535959111995279752468_wp/
      DATA wgk(19)/0.041969810215164246147147541285970_wp/
      DATA wgk(20)/0.043452539701356069316831728117073_wp/
      DATA wgk(21)/0.044814800133162663192355551616723_wp/
      DATA wgk(22)/0.046059238271006988116271735559374_wp/
      DATA wgk(23)/0.047185546569299153945261478181099_wp/
      DATA wgk(24)/0.048185861757087129140779492298305_wp/
      DATA wgk(25)/0.049055434555029778887528165367238_wp/
      DATA wgk(26)/0.049795683427074206357811569379942_wp/
      DATA wgk(27)/0.050405921402782346840893085653585_wp/
      DATA wgk(28)/0.050881795898749606492297473049805_wp/
      DATA wgk(29)/0.051221547849258772170656282604944_wp/
      DATA wgk(30)/0.051426128537459025933862879215781_wp/
      DATA wgk(31)/0.051494729429451567558340433647099_wp/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           dabsc  - abscissa
!           fval*  - function value
!           resg   - result of the 30-point gauss rule
!           resk   - result of the 61-point kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!
      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
      centr = 0.5_wp*(B+A)
      hlgth = 0.5_wp*(B-A)
      dhlgth = abs(hlgth)
!
!           compute the 61-point kronrod approximation to the
!           integral, and estimate the absolute error.
!

      resg = 0.0_wp
      fc = F(centr)
      resk = wgk(31)*fc
      Resabs = abs(resk)
      DO j = 1 , 15
         jtw = j*2
         dabsc = hlgth*xgk(jtw)
         fval1 = F(centr-dabsc)
         fval2 = F(centr+dabsc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(abs(fval1)+abs(fval2))
      ENDDO
      DO j = 1 , 15
         jtwm1 = j*2 - 1
         dabsc = hlgth*xgk(jtwm1)
         fval1 = F(centr-dabsc)
         fval2 = F(centr+dabsc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(abs(fval1)+abs(fval2))
      ENDDO
      reskh = resk*0.5_wp
      Resasc = wgk(31)*abs(fc-reskh)
      DO j = 1 , 30
         Resasc = Resasc + wgk(j)           &
                  *(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = abs((resk-resg)*hlgth)
      IF ( Resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
           Abserr = Resasc*min(1.0_wp,(200.0_wp*Abserr/Resasc) &
           **1.5_wp)
      IF ( Resabs>uflow/(50.0_wp*epmach) )  &
           Abserr = max((epmach*50.0_wp)*Resabs,Abserr)
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQMOMO(Alfa,Beta,Ri,Rj,Rg,Rh,Integr)
      IMPLICIT NONE

!***date written   820101   (yymmdd)
!***revision date  830518   (yymmdd)
!***keywords  modified chebyshev moments
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this routine computes modified chebsyshev moments. the k-th
!            modified chebyshev moment is defined as the integral over
!            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
!            polynomial of degree k.
!***description
!
!        modified chebyshev moments
!        standard fortran subroutine
!        real(wp) version
!
!        parameters
!           alfa   - real(wp)
!                    parameter in the weight function w(x), alfa>(-1)
!
!           beta   - real(wp)
!                    parameter in the weight function w(x), beta>(-1)
!
!           ri     - real(wp)
!                    vector of dimension 25
!                    ri(k) is the integral over (-1,1) of
!                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
!
!           rj     - real(wp)
!                    vector of dimension 25
!                    rj(k) is the integral over (-1,1) of
!                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
!
!           rg     - real(wp)
!                    vector of dimension 25
!                    rg(k) is the integral over (-1,1) of
!                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
!
!           rh     - real(wp)
!                    vector of dimension 25
!                    rh(k) is the integral over (-1,1) of
!                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
!
!           integr - integer
!                    input parameter indicating the modified
!                    moments to be computed
!                    integr = 1 compute ri, rj
!                           = 2 compute ri, rj, rg
!                           = 3 compute ri, rj, rh
!                           = 4 compute ri, rj, rg, rh
!

!
      real(wp) Alfa , alfp1 , alfp2 , an , anm1 , Beta , betp1 ,&
                       betp2 , ralf , rbet , Rg , Rh , Ri , Rj
      INTEGER i , im1 , Integr
!
      DIMENSION Rg(25) , Rh(25) , Ri(25) , Rj(25)
!
!

      alfp1 = Alfa + 1.0_wp
      betp1 = Beta + 1.0_wp
      alfp2 = Alfa + 2.0_wp
      betp2 = Beta + 2.0_wp
      ralf = 2.0_wp**alfp1
      rbet = 2.0_wp**betp1
!
!           compute ri, rj using a forward recurrence relation.
!
      Ri(1) = ralf/alfp1
      Rj(1) = rbet/betp1
      Ri(2) = Ri(1)*Alfa/alfp2
      Rj(2) = Rj(1)*Beta/betp2
      an = 2.0_wp
      anm1 = 1.0_wp
      DO i = 3 , 25
         Ri(i) = -(ralf+an*(an-alfp2)*Ri(i-1))/(anm1*(an+alfp1))
         Rj(i) = -(rbet+an*(an-betp2)*Rj(i-1))/(anm1*(an+betp1))
         anm1 = an
         an = an + 1.0_wp
      ENDDO
      IF ( Integr/=1 ) THEN
         IF ( Integr/=3 ) THEN
!
!           compute rg using a forward recurrence relation.
!
            Rg(1) = -Ri(1)/alfp1
            Rg(2) = -(ralf+ralf)/(alfp2*alfp2) - Rg(1)
            an = 2.0_wp
            anm1 = 1.0_wp
            im1 = 2
            DO i = 3 , 25
               Rg(i) = -(an*(an-alfp2)*Rg(im1)-an*Ri(im1)+anm1*Ri(i)) &
                       /(anm1*(an+alfp1))
               anm1 = an
               an = an + 1.0_wp
               im1 = i
            ENDDO
            IF ( Integr==2 ) GOTO 100
         ENDIF
!
!           compute rh using a forward recurrence relation.
!
         Rh(1) = -Rj(1)/betp1
         Rh(2) = -(rbet+rbet)/(betp2*betp2) - Rh(1)
         an = 2.0_wp
         anm1 = 1.0_wp
         im1 = 2
         DO i = 3 , 25
            Rh(i) = -(an*(an-betp2)*Rh(im1)-an*Rj(im1)+anm1*Rj(i)) &
                    /(anm1*(an+betp1))
            anm1 = an
            an = an + 1.0_wp
            im1 = i
         ENDDO
         DO i = 2 , 25 , 2
            Rh(i) = -Rh(i)
         ENDDO
      ENDIF
 100  DO i = 2 , 25 , 2
         Rj(i) = -Rj(i)
      ENDDO
      END
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQNG(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier)
      IMPLICIT NONE

!***date written   800101   (yymmdd)
!***revision date  810101   (yymmdd)
!***keywords  automatic integrator, smooth integrand,
!             non-adaptive, gauss-kronrod(patterson)
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl math & progr. div. - k.u.leuven
!           kahaner,david,nbs - modified (2/82)
!***purpose  the routine calculates an approximation result to a
!            given definite integral i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
! non-adaptive integration
! standard fortran subroutine
! real(wp) version
!
!           f      - real(wp)
!                    function subprogram defining the integrand function
!                    f(x). the actual name for f needs to be declared
!                    external in the driver program.
!
!           a      - real(wp)
!                    lower limit of integration
!
!           b      - real(wp)
!                    upper limit of integration
!
!           epsabs - real(wp)
!                    absolute accuracy requested
!           epsrel - real(wp)
!                    relative accuracy requested
!                    if  epsabs<=0
!                    and epsrel<max(50*rel.mach.acc.,0.5e-28_wp),
!                    the routine will end with ier = 6.
!
!         on return
!           result - real(wp)
!                    approximation to the integral i
!                    result is obtained by applying the 21-point
!                    gauss-kronrod rule (res21) obtained by optimal
!                    addition of abscissae to the 10-point gauss rule
!                    (res10), or by applying the 43-point rule (res43)
!                    obtained by optimal addition of abscissae to the
!                    21-point gauss-kronrod rule, or by applying the
!                    87-point rule (res87) obtained by optimal addition
!                    of abscissae to the 43-point rule.
!
!           abserr - real(wp)
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed abs(i-result)
!
!           neval  - integer
!                    number of integrand evaluations
!
!           ier    - ier = 0 normal and reliable termination of the
!                            routine. it is assumed that the requested
!                            accuracy has been achieved.
!                    ier>0 abnormal termination of the routine. it is
!                            assumed that the requested accuracy has
!                            not been achieved.
!           error messages
!                    ier = 1 the maximum number of steps has been
!                            executed. the integral is probably too
!                            difficult to be calculated by dqng.
!                        = 6 the input is invalid, because
!                            epsabs<=0 and
!                            epsrel<max(50*rel.mach.acc.,0.5e-28_wp).
!                            result, abserr and neval are set to zero.
!

!
      real(wp) A , absc , Abserr , B , centr , dhlgth , &
                       min , epmach , Epsabs , &
                       Epsrel , fcentr , fval , fval1 , fval2 , &
                       fv1 , fv2 , fv3 , fv4 , hlgth , Result , res10 , &
                       res21 , res43 , res87 , resabs , resasc , reskh ,&
                       savfun , uflow , w10 , w21a , w21b , w43a , &
                       w43b , w87a , w87b , x1 , x2 , x3 , x4
      INTEGER Ier , ipx , k , l , Neval
      procedure(func) :: f
!
      DIMENSION fv1(5) , fv2(5) , fv3(5) , fv4(5) , x1(5) , x2(5) , &
                x3(11) , x4(22) , w10(5) , w21a(5) , w21b(6) , w43a(10) &
                , w43b(12) , w87a(21) , w87b(23) , savfun(21)
!
!           the following data statements contain the
!           abscissae and weights of the integration rules used.
!
!           x1      abscissae common to the 10-, 21-, 43- and 87-
!                   point rule
!           x2      abscissae common to the 21-, 43- and 87-point rule
!           x3      abscissae common to the 43- and 87-point rule
!           x4      abscissae of the 87-point rule
!           w10     weights of the 10-point formula
!           w21a    weights of the 21-point formula for abscissae x1
!           w21b    weights of the 21-point formula for abscissae x2
!           w43a    weights of the 43-point formula for abscissae x1, x3
!           w43b    weights of the 43-point formula for abscissae x3
!           w87a    weights of the 87-point formula for abscissae x1,
!                   x2, x3
!           w87b    weights of the 87-point formula for abscissae x4
!
!
! gauss-kronrod-patterson quadrature coefficients for use in
! quadpack routine qng.  these coefficients were calculated with
! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
!
      DATA x1(1)/0.973906528517171720077964012084452_wp/
      DATA x1(2)/0.865063366688984510732096688423493_wp/
      DATA x1(3)/0.679409568299024406234327365114874_wp/
      DATA x1(4)/0.433395394129247190799265943165784_wp/
      DATA x1(5)/0.148874338981631210884826001129720_wp/
      DATA w10(1)/0.066671344308688137593568809893332_wp/
      DATA w10(2)/0.149451349150580593145776339657697_wp/
      DATA w10(3)/0.219086362515982043995534934228163_wp/
      DATA w10(4)/0.269266719309996355091226921569469_wp/
      DATA w10(5)/0.295524224714752870173892994651338_wp/
!
      DATA x2(1)/0.995657163025808080735527280689003_wp/
      DATA x2(2)/0.930157491355708226001207180059508_wp/
      DATA x2(3)/0.780817726586416897063717578345042_wp/
      DATA x2(4)/0.562757134668604683339000099272694_wp/
      DATA x2(5)/0.294392862701460198131126603103866_wp/
      DATA w21a(1)/0.032558162307964727478818972459390_wp/
      DATA w21a(2)/0.075039674810919952767043140916190_wp/
      DATA w21a(3)/0.109387158802297641899210590325805_wp/
      DATA w21a(4)/0.134709217311473325928054001771707_wp/
      DATA w21a(5)/0.147739104901338491374841515972068_wp/
      DATA w21b(1)/0.011694638867371874278064396062192_wp/
      DATA w21b(2)/0.054755896574351996031381300244580_wp/
      DATA w21b(3)/0.093125454583697605535065465083366_wp/
      DATA w21b(4)/0.123491976262065851077958109831074_wp/
      DATA w21b(5)/0.142775938577060080797094273138717_wp/
      DATA w21b(6)/0.149445554002916905664936468389821_wp/
!
      DATA x3(1)/0.999333360901932081394099323919911_wp/
      DATA x3(2)/0.987433402908088869795961478381209_wp/
      DATA x3(3)/0.954807934814266299257919200290473_wp/
      DATA x3(4)/0.900148695748328293625099494069092_wp/
      DATA x3(5)/0.825198314983114150847066732588520_wp/
      DATA x3(6)/0.732148388989304982612354848755461_wp/
      DATA x3(7)/0.622847970537725238641159120344323_wp/
      DATA x3(8)/0.499479574071056499952214885499755_wp/
      DATA x3(9)/0.364901661346580768043989548502644_wp/
      DATA x3(10)/0.222254919776601296498260928066212_wp/
      DATA x3(11)/0.074650617461383322043914435796506_wp/
      DATA w43a(1)/0.016296734289666564924281974617663_wp/
      DATA w43a(2)/0.037522876120869501461613795898115_wp/
      DATA w43a(3)/0.054694902058255442147212685465005_wp/
      DATA w43a(4)/0.067355414609478086075553166302174_wp/
      DATA w43a(5)/0.073870199632393953432140695251367_wp/
      DATA w43a(6)/0.005768556059769796184184327908655_wp/
      DATA w43a(7)/0.027371890593248842081276069289151_wp/
      DATA w43a(8)/0.046560826910428830743339154433824_wp/
      DATA w43a(9)/0.061744995201442564496240336030883_wp/
      DATA w43a(10)/0.071387267268693397768559114425516_wp/
      DATA w43b(1)/0.001844477640212414100389106552965_wp/
      DATA w43b(2)/0.010798689585891651740465406741293_wp/
      DATA w43b(3)/0.021895363867795428102523123075149_wp/
      DATA w43b(4)/0.032597463975345689443882222526137_wp/
      DATA w43b(5)/0.042163137935191811847627924327955_wp/
      DATA w43b(6)/0.050741939600184577780189020092084_wp/
      DATA w43b(7)/0.058379395542619248375475369330206_wp/
      DATA w43b(8)/0.064746404951445885544689259517511_wp/
      DATA w43b(9)/0.069566197912356484528633315038405_wp/
      DATA w43b(10)/0.072824441471833208150939535192842_wp/
      DATA w43b(11)/0.074507751014175118273571813842889_wp/
      DATA w43b(12)/0.074722147517403005594425168280423_wp/
!
      DATA x4(1)/0.999902977262729234490529830591582_wp/
      DATA x4(2)/0.997989895986678745427496322365960_wp/
      DATA x4(3)/0.992175497860687222808523352251425_wp/
      DATA x4(4)/0.981358163572712773571916941623894_wp/
      DATA x4(5)/0.965057623858384619128284110607926_wp/
      DATA x4(6)/0.943167613133670596816416634507426_wp/
      DATA x4(7)/0.915806414685507209591826430720050_wp/
      DATA x4(8)/0.883221657771316501372117548744163_wp/
      DATA x4(9)/0.845710748462415666605902011504855_wp/
      DATA x4(10)/0.803557658035230982788739474980964_wp/
      DATA x4(11)/0.757005730685495558328942793432020_wp/
      DATA x4(12)/0.706273209787321819824094274740840_wp/
      DATA x4(13)/0.651589466501177922534422205016736_wp/
      DATA x4(14)/0.593223374057961088875273770349144_wp/
      DATA x4(15)/0.531493605970831932285268948562671_wp/
      DATA x4(16)/0.466763623042022844871966781659270_wp/
      DATA x4(17)/0.399424847859218804732101665817923_wp/
      DATA x4(18)/0.329874877106188288265053371824597_wp/
      DATA x4(19)/0.258503559202161551802280975429025_wp/
      DATA x4(20)/0.185695396568346652015917141167606_wp/
      DATA x4(21)/0.111842213179907468172398359241362_wp/
      DATA x4(22)/0.037352123394619870814998165437704_wp/
      DATA w87a(1)/0.008148377384149172900002878448190_wp/
      DATA w87a(2)/0.018761438201562822243935059003794_wp/
      DATA w87a(3)/0.027347451050052286161582829741283_wp/
      DATA w87a(4)/0.033677707311637930046581056957588_wp/
      DATA w87a(5)/0.036935099820427907614589586742499_wp/
      DATA w87a(6)/0.002884872430211530501334156248695_wp/
      DATA w87a(7)/0.013685946022712701888950035273128_wp/
      DATA w87a(8)/0.023280413502888311123409291030404_wp/
      DATA w87a(9)/0.030872497611713358675466394126442_wp/
      DATA w87a(10)/0.035693633639418770719351355457044_wp/
      DATA w87a(11)/0.000915283345202241360843392549948_wp/
      DATA w87a(12)/0.005399280219300471367738743391053_wp/
      DATA w87a(13)/0.010947679601118931134327826856808_wp/
      DATA w87a(14)/0.016298731696787335262665703223280_wp/
      DATA w87a(15)/0.021081568889203835112433060188190_wp/
      DATA w87a(16)/0.025370969769253827243467999831710_wp/
      DATA w87a(17)/0.029189697756475752501446154084920_wp/
      DATA w87a(18)/0.032373202467202789685788194889595_wp/
      DATA w87a(19)/0.034783098950365142750781997949596_wp/
      DATA w87a(20)/0.036412220731351787562801163687577_wp/
      DATA w87a(21)/0.037253875503047708539592001191226_wp/
      DATA w87b(1)/0.000274145563762072350016527092881_wp/
      DATA w87b(2)/0.001807124155057942948341311753254_wp/
      DATA w87b(3)/0.004096869282759164864458070683480_wp/
      DATA w87b(4)/0.006758290051847378699816577897424_wp/
      DATA w87b(5)/0.009549957672201646536053581325377_wp/
      DATA w87b(6)/0.012329447652244853694626639963780_wp/
      DATA w87b(7)/0.015010447346388952376697286041943_wp/
      DATA w87b(8)/0.017548967986243191099665352925900_wp/
      DATA w87b(9)/0.019938037786440888202278192730714_wp/
      DATA w87b(10)/0.022194935961012286796332102959499_wp/
      DATA w87b(11)/0.024339147126000805470360647041454_wp/
      DATA w87b(12)/0.026374505414839207241503786552615_wp/
      DATA w87b(13)/0.028286910788771200659968002987960_wp/
      DATA w87b(14)/0.030052581128092695322521110347341_wp/
      DATA w87b(15)/0.031646751371439929404586051078883_wp/
      DATA w87b(16)/0.033050413419978503290785944862689_wp/
      DATA w87b(17)/0.034255099704226061787082821046821_wp/
      DATA w87b(18)/0.035262412660156681033782717998428_wp/
      DATA w87b(19)/0.036076989622888701185500318003895_wp/
      DATA w87b(20)/0.036698604498456094498018047441094_wp/
      DATA w87b(21)/0.037120549269832576114119958413599_wp/
      DATA w87b(22)/0.037334228751935040321235449094698_wp/
      DATA w87b(23)/0.037361073762679023410321241766599_wp/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fcentr - function value at mid point
!           absc   - abscissa
!           fval   - function value
!           savfun - array of function values which have already been
!                    computed
!           res10  - 10-point gauss result
!           res21  - 21-point kronrod result
!           res43  - 43-point result
!           res87  - 87-point result
!           resabs - approximation to the integral of abs(f)
!           resasc - approximation to the integral of abs(f-i/(b-a))
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!

      epmach = D1MACH(4)
      uflow = D1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Result = 0.0_wp
      Abserr = 0.0_wp
      Neval = 0
      Ier = 6
      IF ( Epsabs>0.0_wp .OR. Epsrel>=max(50.0_wp*epmach,0.5e-28_wp) ) &
           THEN
         hlgth = 0.5_wp*(B-A)
         dhlgth = abs(hlgth)
         centr = 0.5_wp*(B+A)
         fcentr = F(centr)
         Neval = 21
         Ier = 1
!
!          compute the integral using the 10- and 21-point formula.
!
         DO l = 1 , 3
            SELECT CASE (l)
            CASE (2)
!
!          compute the integral using the 43-point formula.
!
               res43 = w43b(12)*fcentr
               Neval = 43
               DO k = 1 , 10
                  res43 = res43 + savfun(k)*w43a(k)
               ENDDO
               DO k = 1 , 11
                  ipx = ipx + 1
                  absc = hlgth*x3(k)
                  fval = F(absc+centr) + F(centr-absc)
                  res43 = res43 + fval*w43b(k)
                  savfun(ipx) = fval
               ENDDO
!
!          test for convergence.
!
               Result = res43*hlgth
               Abserr = abs((res43-res21)*hlgth)
            CASE (3)
!
!          compute the integral using the 87-point formula.
!
               res87 = w87b(23)*fcentr
               Neval = 87
               DO k = 1 , 21
                  res87 = res87 + savfun(k)*w87a(k)
               ENDDO
               DO k = 1 , 22
                  absc = hlgth*x4(k)
                  res87 = res87 + w87b(k)*(F(absc+centr)+F(centr-absc))
               ENDDO
               Result = res87*hlgth
               Abserr = abs((res87-res43)*hlgth)
            CASE DEFAULT
               res10 = 0.0_wp
               res21 = w21b(6)*fcentr
               resabs = w21b(6)*abs(fcentr)
               DO k = 1 , 5
                  absc = hlgth*x1(k)
                  fval1 = F(centr+absc)
                  fval2 = F(centr-absc)
                  fval = fval1 + fval2
                  res10 = res10 + w10(k)*fval
                  res21 = res21 + w21a(k)*fval
                  resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
                  savfun(k) = fval
                  fv1(k) = fval1
                  fv2(k) = fval2
               ENDDO
               ipx = 5
               DO k = 1 , 5
                  ipx = ipx + 1
                  absc = hlgth*x2(k)
                  fval1 = F(centr+absc)
                  fval2 = F(centr-absc)
                  fval = fval1 + fval2
                  res21 = res21 + w21b(k)*fval
                  resabs = resabs + w21b(k)*(abs(fval1)+abs(fval2))
                  savfun(ipx) = fval
                  fv3(k) = fval1
                  fv4(k) = fval2
               ENDDO
!
!          test for convergence.
!
               Result = res21*hlgth
               resabs = resabs*dhlgth
               reskh = 0.5_wp*res21
               resasc = w21b(6)*abs(fcentr-reskh)
               DO k = 1 , 5
                  resasc = resasc + w21a(k) &
                           *(abs(fv1(k)-reskh)+abs(fv2(k)-reskh)) &
                           + w21b(k)        &
                           *(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
               ENDDO
               Abserr = abs((res21-res10)*hlgth)
               resasc = resasc*dhlgth
            END SELECT
            IF ( resasc/=0.0_wp .AND. Abserr/=0.0_wp ) &
                 Abserr = resasc*min(1.0_wp,(200.0_wp*Abserr/resasc) &
                 **1.5_wp)
            IF ( resabs>uflow/(50.0_wp*epmach) )   &
                 Abserr = max((epmach*50.0_wp)*resabs,Abserr)
            IF ( Abserr<=max(Epsabs,Epsrel*abs(Result)) ) Ier = 0
! ***jump out of do-loop
            IF ( Ier==0 ) return
         ENDDO
      ENDIF
      CALL XERROR('abnormal return from dqng ',26,Ier,0)
    end
!********************************************************************************

!********************************************************************************
      SUBROUTINE DQPSRT(Limit,Last,Maxerr,Ermax,Elist,Iord,Nrmax)
      IMPLICIT NONE

!### See also
!  *  dqage,dqagie,dqagpe,dqawse
!***revision date  810101   (yymmdd)
!***keywords  sequential sorting
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this routine maintains the descending ordering in the
!            list of the local error estimated resulting from the
!            interval subdivision process. at each call two error
!            estimates are inserted using the sequential search
!            method, top-down for the largest error estimate and
!            bottom-up for the smallest error estimate.
!***description
!
!           ordering routine
!           standard fortran subroutine
!           real(wp) version
!
!           parameters (meaning at output)
!              limit  - integer
!                       maximum number of error estimates the list
!                       can contain
!
!              last   - integer
!                       number of error estimates currently in the list
!
!              maxerr - integer
!                       maxerr points to the nrmax-th largest error
!                       estimate currently in the list
!
!              ermax  - real(wp)
!                       nrmax-th largest error estimate
!                       ermax = elist(maxerr)
!
!              elist  - real(wp)
!                       vector of dimension last containing
!                       the error estimates
!
!              iord   - integer
!                       vector of dimension last, the first k elements
!                       of which contain pointers to the error
!                       estimates, such that
!                       elist(iord(1)),...,  elist(iord(k))
!                       form a decreasing sequence, with
!                       k = last if last<=(limit/2+2), and
!                       k = limit+1-last otherwise
!
!              nrmax  - integer
!                       maxerr = iord(nrmax)
!
!
      real(wp) Elist , Ermax , errmax , errmin
      INTEGER i , ibeg , ido , Iord , isucc , j , jbnd , jupbn , k , &
              Last , Limit , Maxerr , Nrmax
      DIMENSION Elist(Last) , Iord(Last)
!
!           check whether the list contains more than
!           two error estimates.
!

      IF ( Last>2 ) THEN
!
!           this part of the routine is only executed if, due to a
!           difficult integrand, subdivision increased the error
!           estimate. in the normal case the insert procedure should
!           start after the nrmax-th largest error estimate.
!
         errmax = Elist(Maxerr)
         IF ( Nrmax/=1 ) THEN
            ido = Nrmax - 1
            DO i = 1 , ido
               isucc = Iord(Nrmax-1)
! ***jump out of do-loop
               IF ( errmax<=Elist(isucc) ) GOTO 50
               Iord(Nrmax) = isucc
               Nrmax = Nrmax - 1
            ENDDO
         ENDIF
!
!           compute the number of elements in the list to be maintained
!           in descending order. this number depends on the number of
!           subdivisions still allowed.
!
 50      jupbn = Last
         IF ( Last>(Limit/2+2) ) jupbn = Limit + 3 - Last
         errmin = Elist(Last)
!
!           insert errmax by traversing the list top-down,
!           starting comparison from the element elist(iord(nrmax+1)).
!
         jbnd = jupbn - 1
         ibeg = Nrmax + 1
         IF ( ibeg<=jbnd ) THEN
            DO i = ibeg , jbnd
               isucc = Iord(i)
! ***jump out of do-loop
               IF ( errmax>=Elist(isucc) ) GOTO 100
               Iord(i-1) = isucc
            ENDDO
         ENDIF
         Iord(jbnd) = Maxerr
         Iord(jupbn) = Last
      ELSE
         Iord(1) = 1
         Iord(2) = 2
      ENDIF
      GOTO 300
!
!           insert errmin by traversing the list bottom-up.
!
 100  Iord(i-1) = Maxerr
      k = jbnd
      DO j = i , jbnd
         isucc = Iord(k)
! ***jump out of do-loop
         IF ( errmin<Elist(isucc) ) GOTO 200
         Iord(k+1) = isucc
         k = k - 1
      ENDDO
      Iord(i) = Last
      GOTO 300
 200  Iord(k+1) = Last
!
!           set maxerr and ermax.
!
 300  Maxerr = Iord(Nrmax)
      Ermax = Elist(Maxerr)
      END
!********************************************************************************

!********************************************************************************
      real(wp) FUNCTION DQWGTC(X,C,P2,P3,P4,Kp)
      IMPLICIT NONE

!### See also
!  * dqk15w
!***revision date  810101   (yymmdd)
!***keywords  weight function, cauchy principal value
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this function subprogram is used together with the
!            routine qawc and defines the weight function.
!
      real(wp) C , P2 , P3 , P4 , X
      INTEGER Kp

      DQWGTC = 1.0_wp/(X-C)
      END
!********************************************************************************

!********************************************************************************
      real(wp) FUNCTION DQWGTF(X,Omega,P2,P3,P4,Integr)
      IMPLICIT NONE

!### See also
!  *   dqk15w
!***revision date 810101   (yymmdd)
!***keywords  cos or sin in weight function
!***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. * progr. div. - k.u.leuven
!
      real(wp) Omega , omx , P2 , P3 , P4 , X
      INTEGER Integr

      omx = Omega*X
      IF ( Integr==2 ) THEN
         DQWGTF = sin(omx)
      ELSE
         DQWGTF = cos(omx)
      ENDIF
      END
!********************************************************************************

!********************************************************************************
      real(wp) FUNCTION DQWGTS(X,A,B,Alfa,Beta,Integr)
      IMPLICIT NONE

!### See also
!  * dqk15w
!***revision date  810101   (yymmdd)
!***keywords  weight function, algebraico-logarithmic
!             end-point singularities
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this function subprogram is used together with the
!            routine dqaws and defines the weight function.
!
      real(wp) A , Alfa , B , Beta , bmx , X , xma
      INTEGER Integr

      xma = X - A
      bmx = B - X
      DQWGTS = xma**Alfa*bmx**Beta
      SELECT CASE (Integr)
      CASE (1)
      CASE (3)
         DQWGTS = DQWGTS*log(bmx)
      CASE (4)
         DQWGTS = DQWGTS*log(xma)*log(bmx)
      CASE DEFAULT
         DQWGTS = DQWGTS*log(xma)
      END SELECT
      END

   !===================================================================
   ! additional routines
   !===================================================================

   !********************************************************************************
   !>
   !     dgtsl given a general tridiagonal matrix and a right hand
   !     side will find the solution.
   !
   !     on entry
   !
   !        n       integer
   !                is the order of the tridiagonal matrix.
   !
   !        c       real(wp)(n)
   !                is the subdiagonal of the tridiagonal matrix.
   !                c(2) through c(n) should contain the subdiagonal.
   !                on output c is destroyed.
   !
   !        d       real(wp)(n)
   !                is the diagonal of the tridiagonal matrix.
   !                on output d is destroyed.
   !
   !        e       real(wp)(n)
   !                is the superdiagonal of the tridiagonal matrix.
   !                e(1) through e(n-1) should contain the superdiagonal.
   !                on output e is destroyed.
   !
   !        b       real(wp)(n)
   !                is the right hand side vector.
   !
   !     on return
   !
   !        b       is the solution vector.
   !
   !        info    integer
   !                = 0 normal value.
   !                = k if the k-th element of the diagonal becomes
   !                    exactly zero.  the subroutine returns when
   !                    this is detected.
   !
   !     linpack. this version dated 08/14/78 .
   !     jack dongarra, argonne national laboratory.

    subroutine dgtsl(n,c,d,e,b,info)
    implicit none
    integer n,info
    real(wp) c(*),d(*),e(*),b(*)

    integer k,kb,kp1,nm1,nm2
    real(wp) t

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
        return
    end if

    ! back solve

    nm2 = n - 2
    b(n) = b(n)/c(n)
    if (n /= 1) then
        b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
        if (nm2 >= 1) then
            do kb = 1, nm2
                k = nm2 - kb + 1
                b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
            end do
        end if
    end if

    end subroutine dgtsl
!********************************************************************************

!********************************************************************************
!>
!
! This function is intended to replace the old D1MACH by using F90
! intrinsic functions.
!
! The traditional D1MACH constants are:
!
!  * D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  * D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  * D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  * D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  * D1MACH( 5) = LOG10(B)

    function d1mach(i)
    implicit none

    integer,intent(in) :: i
    real(wp) :: d1mach

    real(wp), dimension(5), parameter :: d1mach_values = [ tiny(1.0d0), &
                                                           huge(1.0d0), &
                                                           real(radix(1.0d0),&
                                                           kind(1.0d0))**(-digits(1.0d0)), &
                                                           epsilon(1.0d0), &
                                                           log10(real(radix(1.0d0),&
                                                           kind(1.0d0))) ]

    if (i<1 .or. i>5) then
        write (*,'(1x,''d1mach(i) - i out of bounds, i ='',i10)') i
        error stop ' d1mach(i) - i out of bounds'
    end if
    d1mach = d1mach_values(i)

    end function d1mach
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

    subroutine xerror(messg,nmessg,nerr,level)
    implicit none
    character(len=*),intent(in) :: messg !! message to be processed
    integer,intent(in) :: nmessg !! the actual number of characters in MESSG
    integer,intent(in) :: nerr  !! the error number associated with this message.
                                !! NERR must not be zero.
    integer,intent(in) :: level !! error category:
                                !!  * =2 means this is an unconditionally fatal error.
                                !!  * =1 means this is a recoverable error.  (I.e., it is
                                !!    non-fatal if XSETF has been appropriately called.)
                                !!  * =0 means this is a warning message only.
                                !!  * =-1 means this is a warning message which is to be
                                !!    printed at most once, regardless of how many
                                !!    times this call is executed.

    !call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)

    write(*,*) nerr, messg(1:nmessg)
    if (level==2) error stop

    end subroutine xerror
!********************************************************************************

end module quadpack
