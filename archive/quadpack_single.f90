!*==QAG.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAG(F,A,B,Epsabs,Epsrel,Key,Result,Abserr,Neval,Ier,   &
                   & Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAG8593
!*** Start of declarations inserted by SPAG
      INTEGER Last
!*** End of declarations inserted by SPAG
!***begin prologue  qag
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
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
!        real version
!
!            f      - real
!                     function subprogam defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
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
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!                    produced in the subdivision process, which
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
!            work  - real
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
!
!***references  (none)
!***routines called  qage,xerror
!***end prologue  qag
!
      REAL A , Abserr , B , Epsabs , Epsrel , F , Result , Work
      INTEGER Ier , Iwork , Key , Lenw , Limit , lvl , l1 , l2 , l3 ,   &
            & Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of lenw.
!
!***first executable statement  qag
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for qage.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL QAGE(F,A,B,Epsabs,Epsrel,Key,Limit,Result,Abserr,Neval,   &
                 & Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qag ',26,Ier,lvl)
      END
!*==QAGE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGE(F,A,B,Epsabs,Epsrel,Key,Limit,Result,Abserr,Neval,&
                    & Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE
!*--QAGE8783
!***begin prologue  qage
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
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
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
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
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28),
!                             result, abserr, neval, last, rlist(1) ,
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!            alist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            blist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            rlist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the
!                      integral approximations on the subintervals
!
!            elist   - real
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
!***references  (none)
!***routines called  qk15,qk21,qk31,qk41,qk51,qk61,qpsrt,r1mach
!***end prologue  qage
!
      REAL A , Abserr , Alist , area , area1 , area12 , area2 , a1 ,    &
         & a2 , B , Blist , b1 , b2 , defabs , defab1 , defab2 ,        &
         & R1MACH , Elist , epmach , Epsabs , Epsrel , errbnd , errmax ,&
         & error1 , error2 , erro12 , errsum , F , resabs , Result ,    &
         & Rlist , uflow
      INTEGER Ier , Iord , iroff1 , iroff2 , k , Key , keyf , Last ,    &
            & Limit , maxerr , Neval , nrmax
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
              & , Rlist(Limit)
!
      EXTERNAL F
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
!***first executable statement  qage
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      IF ( Epsabs<=0.0E+00 .AND. Epsrel<AMAX1(0.5E+02*epmach,0.5E-14) ) &
         & Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         keyf = Key
         IF ( Key<=0 ) keyf = 1
         IF ( Key>=7 ) keyf = 6
         Neval = 0
         IF ( keyf==1 ) CALL QK15(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==2 ) CALL QK21(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==3 ) CALL QK31(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==4 ) CALL QK41(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==5 ) CALL QK51(F,A,B,Result,Abserr,defabs,resabs)
         IF ( keyf==6 ) CALL QK61(F,A,B,Result,Abserr,defabs,resabs)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
!
!           test on accuracy.
!
         errbnd = AMAX1(Epsabs,Epsrel*ABS(Result))
         IF ( Abserr<=0.5E+02*epmach*defabs .AND. Abserr>errbnd )       &
            & Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( .NOT.(Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs)   &
            & .OR. Abserr==0.0E+00) ) THEN
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
               b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               IF ( keyf==1 ) CALL QK15(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==2 ) CALL QK21(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==3 ) CALL QK31(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==4 ) CALL QK41(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==5 ) CALL QK51(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==6 ) CALL QK61(F,a1,b1,area1,error1,resabs,    &
                                      & defab1)
               IF ( keyf==1 ) CALL QK15(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
               IF ( keyf==2 ) CALL QK21(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
               IF ( keyf==3 ) CALL QK31(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
               IF ( keyf==4 ) CALL QK41(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
               IF ( keyf==5 ) CALL QK51(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
               IF ( keyf==6 ) CALL QK61(F,a2,b2,area2,error2,resabs,    &
                                      & defab2)
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
                  IF ( ABS(Rlist(maxerr)-area12)<=0.1E-04*ABS(area12)   &
                     & .AND. erro12>=0.99E+00*errmax ) iroff1 = iroff1 +&
                     & 1
                  IF ( Last>10 .AND. erro12>errmax ) iroff2 = iroff2 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
               IF ( errsum>errbnd ) THEN
!
!           test for roundoff error and eventually
!           set error flag.
!
                  IF ( iroff1>=6 .OR. iroff2>=20 ) Ier = 2
!
!           set error flag in the case that the number of
!           subintervals equals limit.
!
                  IF ( Last==Limit ) Ier = 1
!
!           set error flag in the case of bad integrand behaviour
!           at a point of the integration range.
!
                  IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach) &
                     & *(ABS(a2)+0.1E+04*uflow) ) Ier = 3
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with the largest error estimate (to be
!           bisected next).
!
               CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0E+00
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
         IF ( keyf/=1 ) Neval = (10*keyf+1)*(2*Neval+1)
         IF ( keyf==1 ) Neval = 30*Neval + 15
      ENDIF
      END
!*==QAGI.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGI(F,Bound,Inf,Epsabs,Epsrel,Result,Abserr,Neval,Ier,&
                    & Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAGI9146
!*** Start of declarations inserted by SPAG
      REAL Bound
      INTEGER Inf , Last
!*** End of declarations inserted by SPAG
!***begin prologue  qagi
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a1,h2a4a1
!***keywords  automatic integrator, infinite intervals,
!             general-purpose, transformation, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. -k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!                    or i = integral of f over (-infinity,bound)
!                    or i = integral of f over (-infinity,+infinity)
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        integration over infinite intervals
!        standard fortran subroutine
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            bound  - real
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - integer
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!            work  - real
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
!***references  (none)
!***routines called  qagie,xerror
!***end prologue  qagi
!
      REAL Abserr , Epsabs , Epsrel , F , Result , Work
      INTEGER Ier , Iwork , Lenw , Limit , lvl , l1 , l2 , l3 , Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limit and lenw.
!
!***first executable statement  qagi
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for qagie.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL QAGIE(F,Bound,Inf,Epsabs,Epsrel,Limit,Result,Abserr,Neval,&
                  & Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qagi',26,Ier,lvl)
      END
!*==QAGIE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGIE(F,Bound,Inf,Epsabs,Epsrel,Limit,Result,Abserr,   &
                     & Neval,Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE
!*--QAGIE9344
!***begin prologue  qagie
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a1,h2a4a1
!***keywords  automatic integrator, infinite intervals,
!             general-purpose, transformation, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math & progr. div - k.u.leuven
!           de doncker,elise,appl. math & progr. div - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!                    or i = integral of f over (-infinity,bound)
!                    or i = integral of f over (-infinity,+infinity),
!                    hopefully satisfying following claim for accuracy
!                    abs(i-result)<=max(epsabs,epsrel*abs(i))
!***description
!
! integration over infinite intervals
! standard fortran subroutine
!
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            bound  - real
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - real
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=1
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28),
!                             result, abserr, neval, last, rlist(1),
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to 0
!                             and 1 respectively.
!
!            alist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            blist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            rlist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real
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
!***references  (none)
!***routines called  qelg,qk15i,qpsrt,r1mach
!***end prologue  qagie
!
      REAL abseps , Abserr , Alist , area , area1 , area12 , area2 ,    &
         & a1 , a2 , Blist , boun , Bound , b1 , b2 , correc , defabs , &
         & defab1 , defab2 , dres , R1MACH , Elist , epmach , Epsabs ,  &
         & Epsrel , erlarg , erlast , errbnd , errmax , error1 ,        &
         & error2 , erro12 , errsum , ertest , F , oflow , resabs ,     &
         & reseps , Result , res3la , Rlist , rlist2 , small , uflow
      INTEGER id , Ier , ierro , Inf , Iord , iroff1 , iroff2 , iroff3 ,&
            & jupbnd , k , ksgn , ktmin , Last , Limit , maxerr ,       &
            & Neval , nres , nrmax , numrl2
      LOGICAL extrap , noext
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
              & , res3la(3) , Rlist(Limit) , rlist2(52)
!
      EXTERNAL F
!
!            the dimension of rlist2 is determined by the value of
!            limexp in subroutine qelg.
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
      epmach = R1MACH(4)
!
!           test on validity of parameters
!           -----------------------------
!
!***first executable statement  qagie
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      Alist(1) = 0.0E+00
      Blist(1) = 0.1E+01
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      IF ( Epsabs<=0.0E+00 .AND. Epsrel<AMAX1(0.5E+02*epmach,0.5E-14) ) &
         & Ier = 6
      IF ( Ier==6 ) GOTO 99999
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
      IF ( Inf==2 ) boun = 0.0E+00
      CALL QK15I(F,boun,Inf,0.0E+00,0.1E+01,Result,Abserr,defabs,resabs)
!
!           test on accuracy
!
      Last = 1
      Rlist(1) = Result
      Elist(1) = Abserr
      Iord(1) = 1
      dres = ABS(Result)
      errbnd = AMAX1(Epsabs,Epsrel*dres)
      IF ( Abserr<=1.0E+02*epmach*defabs .AND. Abserr>errbnd ) Ier = 2
      IF ( Limit==1 ) Ier = 1
      IF ( Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) .OR.       &
         & Abserr==0.0E+00 ) GOTO 400
!
!           initialization
!           --------------
!
      uflow = R1MACH(1)
      oflow = R1MACH(2)
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
      IF ( dres>=(0.1E+01-0.5E+02*epmach)*defabs ) ksgn = 1
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
         b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
         a2 = b1
         b2 = Blist(maxerr)
         erlast = errmax
         CALL QK15I(F,boun,Inf,a1,b1,area1,error1,resabs,defab1)
         CALL QK15I(F,boun,Inf,a2,b2,area2,error2,resabs,defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
         area12 = area1 + area2
         erro12 = error1 + error2
         errsum = errsum + erro12 - errmax
         area = area + area12 - Rlist(maxerr)
         IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
            IF ( ABS(Rlist(maxerr)-area12)<=0.1E-04*ABS(area12) .AND.   &
               & erro12>=0.99E+00*errmax ) THEN
               IF ( extrap ) iroff2 = iroff2 + 1
               IF ( .NOT.extrap ) iroff1 = iroff1 + 1
            ENDIF
            IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
         ENDIF
         Rlist(maxerr) = area1
         Rlist(Last) = area2
         errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
!
!           test for roundoff error and eventually
!           set error flag.
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
         IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach)          &
            & *(ABS(a2)+0.1E+04*uflow) ) Ier = 4
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
         CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
         IF ( errsum<=errbnd ) GOTO 300
         IF ( Ier/=0 ) GOTO 200
         IF ( Last==2 ) THEN
            small = 0.375E+00
            erlarg = errsum
            ertest = errbnd
            rlist2(2) = area
         ELSEIF ( .NOT.(noext) ) THEN
            erlarg = erlarg - erlast
            IF ( ABS(b1-a1)>small ) erlarg = erlarg + erro12
            IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
               IF ( ABS(Blist(maxerr)-Alist(maxerr))>small ) GOTO 100
               extrap = .TRUE.
               nrmax = 2
            ENDIF
            IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors
!           over the larger intervals (erlarg) and perform
!           extrapolation.
!
               id = nrmax
               jupbnd = Last
               IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
               DO k = id , jupbnd
                  maxerr = Iord(nrmax)
                  errmax = Elist(maxerr)
                  IF ( ABS(Blist(maxerr)-Alist(maxerr))>small ) GOTO 100
                  nrmax = nrmax + 1
               ENDDO
            ENDIF
!
!           perform extrapolation.
!
            numrl2 = numrl2 + 1
            rlist2(numrl2) = area
            CALL QELG(numrl2,rlist2,reseps,abseps,res3la,nres)
            ktmin = ktmin + 1
            IF ( ktmin>5 .AND. Abserr<0.1E-02*errsum ) Ier = 5
            IF ( abseps<Abserr ) THEN
               ktmin = 0
               Abserr = abseps
               Result = reseps
               correc = erlarg
               ertest = AMAX1(Epsabs,Epsrel*ABS(reseps))
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
            small = small*0.5E+00
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
            IF ( Result==0.0E+00 .OR. area==0.0E+00 ) THEN
               IF ( Abserr>errsum ) GOTO 300
               IF ( area==0.0E+00 ) GOTO 400
            ELSEIF ( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
               GOTO 300
            ENDIF
         ENDIF
!
!           test on divergence
!
         IF ( ksgn/=(-1) .OR. AMAX1(ABS(Result),ABS(area))              &
            & >defabs*0.1E-01 ) THEN
            IF ( 0.1E-01>(Result/area) .OR. (Result/area)>0.1E+03 .OR.  &
               & errsum>ABS(area) ) Ier = 6
         ENDIF
         GOTO 400
      ENDIF
!
!           compute global integral sum.
!
 300  Result = 0.0E+00
      DO k = 1 , Last
         Result = Result + Rlist(k)
      ENDDO
      Abserr = errsum
 400  Neval = 30*Last - 15
      IF ( Inf==2 ) Neval = 2*Neval
      IF ( Ier>2 ) Ier = Ier - 1
99999 END
!*==QAGP.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGP(F,A,B,Npts2,Points,Epsabs,Epsrel,Result,Abserr,   &
                    & Neval,Ier,Leniw,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAGP9814
!*** Start of declarations inserted by SPAG
      INTEGER l4 , Last
!*** End of declarations inserted by SPAG
!***begin prologue  qagp
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
!***keywords  automatic integrator, general-purpose,
!             singularities at user specified points,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            break points of the integration interval, where local
!            difficulties of the integrand may occur(e.g. singularities,
!            discontinuities), are provided by the user.
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            npts2  - integer
!                     number equal to two more than the number of
!                     user-supplied break points within the integration
!                     range, npts>=2.
!                     if npts2<2, the routine will end with ier = 6.
!
!            points - real
!                     vector of dimension npts2, the first (npts2-2)
!                     elements of which are the user provided break
!                     points. if these points do not constitute an
!                     ascending sequence there will be an automatic
!                     sorting.
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!            work  - real
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
!***references  (none)
!***routines called  qagpe,xerror
!***end prologue  qagp
!
      REAL A , Abserr , B , Epsabs , Epsrel , F , Points , Result , Work
      INTEGER Ier , Iwork , Leniw , Lenw , limit , lvl , l1 , l2 , l3 , &
            & Neval , Npts2
!
      DIMENSION Iwork(Leniw) , Points(Npts2) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limit and lenw.
!
!***first executable statement  qagp
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Leniw>=(3*Npts2-2) .AND. Lenw>=(Leniw*2-Npts2) .AND.         &
         & Npts2>=2 ) THEN
!
!         prepare call for qagpe.
!
         limit = (Leniw-Npts2)/2
         l1 = limit + 1
         l2 = limit + l1
         l3 = limit + l2
         l4 = limit + l3
!
         CALL QAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,limit,Result,      &
                  & Abserr,Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3),&
                  & Work(l4),Iwork(1),Iwork(l1),Iwork(l2),Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qagp',26,Ier,lvl)
      END
!*==QAGPE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,Limit,Result,   &
                     & Abserr,Neval,Ier,Alist,Blist,Rlist,Elist,Pts,    &
                     & Iord,Level,Ndin,Last)
      IMPLICIT NONE
!*--QAGPE10046
!***begin prologue  qagpe
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
!***keywords  automatic integrator, general-purpose,
!             singularities at user specified points,
!             extrapolation, globally adaptive.
!***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f over (a,b),hopefully
!            satisfying following claim for accuracy abs(i-result)<=
!            max(epsabs,epsrel*abs(i)). break points of the integration
!            interval, where local difficulties of the integrand may
!            occur(e.g. singularities,discontinuities),provided by user.
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            npts2  - integer
!                     number equal to two more than the number of
!                     user-supplied break points within the integration
!                     range, npts2>=2.
!                     if npts2<2, the routine will end with ier = 6.
!
!            points - real
!                     vector of dimension npts2, the first (npts2-2)
!                     elements of which are the user provided break
!                     points. if these points do not constitute an
!                     ascending sequence there will be an automatic
!                     sorting.
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=npts2
!                     if limit<npts2, the routine will end with
!                     ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
!                             or limit<npts2.
!                             result, abserr, neval, last, rlist(1),
!                             and elist(1) are set to zero. alist(1) and
!                             blist(1) are set to a and b respectively.
!
!            alist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            blist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            rlist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            pts    - real
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
!***references  (none)
!***routines called  qelg,qk21,qpsrt,r1mach
!***end prologue  qagpe
      REAL A , abseps , Abserr , Alist , area , area1 , area12 , area2 ,&
         & a1 , a2 , B , Blist , b1 , b2 , correc , defabs , defab1 ,   &
         & defab2 , dres , R1MACH , Elist , epmach , Epsabs , Epsrel ,  &
         & erlarg , erlast , errbnd , errmax , error1 , erro12 ,        &
         & error2 , errsum , ertest , F , oflow , Points , Pts , resa , &
         & resabs , reseps , Result , res3la , Rlist , rlist2 , sign ,  &
         & temp , uflow
      INTEGER i , id , Ier , ierro , ind1 , ind2 , Iord , ip1 , iroff1 ,&
            & iroff2 , iroff3 , j , jlow , jupbnd , k , ksgn , ktmin ,  &
            & Last , levcur , Level , levmax , Limit , maxerr , Ndin ,  &
            & Neval , nint , nintp1 , npts , Npts2 , nres , nrmax ,     &
            & numrl2
      LOGICAL extrap , noext
!
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
              & , Level(Limit) , Ndin(Npts2) , Points(Npts2) ,          &
              & Pts(Npts2) , res3la(3) , Rlist(Limit) , rlist2(52)
!
      EXTERNAL F
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
!           numrl2    - number of elements in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained, it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
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
!***first executable statement  qagpe
      epmach = R1MACH(4)
!
!            test on validity of parameters
!            -----------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      Level(1) = 0
      npts = Npts2 - 2
      IF ( Npts2<2 .OR. Limit<=npts .OR.                                &
         & (Epsabs<=0.0E+00 .AND. Epsrel<AMAX1(0.5E+02*epmach,0.5E-14)) &
         & ) Ier = 6
      IF ( Ier/=6 ) THEN
!
!            if any break points are provided, sort them into an
!            ascending sequence.
!
         sign = 1.0E+00
         IF ( A>B ) sign = -1.0E+00
         Pts(1) = AMIN1(A,B)
         IF ( npts/=0 ) THEN
            DO i = 1 , npts
               Pts(i+1) = Points(i)
            ENDDO
         ENDIF
         Pts(npts+2) = AMAX1(A,B)
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
            IF ( Pts(1)/=AMIN1(A,B) .OR. Pts(nintp1)/=AMAX1(A,B) )      &
               & Ier = 6
            IF ( Ier==6 ) GOTO 99999
         ENDIF
!
!            compute first integral and error approximations.
!            ------------------------------------------------
!
         resabs = 0.0E+00
         DO i = 1 , nint
            b1 = Pts(i+1)
            CALL QK21(F,a1,b1,area1,error1,defabs,resa)
            Abserr = Abserr + error1
            Result = Result + area1
            Ndin(i) = 0
            IF ( error1==resa .AND. error1/=0.0E+00 ) Ndin(i) = 1
            resabs = resabs + defabs
            Level(i) = 0
            Elist(i) = error1
            Alist(i) = a1
            Blist(i) = b1
            Rlist(i) = area1
            Iord(i) = i
            a1 = b1
         ENDDO
         errsum = 0.0E+00
         DO i = 1 , nint
            IF ( Ndin(i)==1 ) Elist(i) = Abserr
            errsum = errsum + Elist(i)
         ENDDO
!
!           test on accuracy.
!
         Last = nint
         Neval = 21*nint
         dres = ABS(Result)
         errbnd = AMAX1(Epsabs,Epsrel*dres)
         IF ( Abserr<=0.1E+03*epmach*resabs .AND. Abserr>errbnd )       &
            & Ier = 2
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
         IF ( Ier/=0 .OR. Abserr<=errbnd ) GOTO 99999
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
         uflow = R1MACH(1)
         oflow = R1MACH(2)
         Abserr = oflow
         ksgn = -1
         IF ( dres>=(0.1E+01-0.5E+02*epmach)*resabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
         DO Last = Npts2 , Limit
!
!           bisect the subinterval with the nrmax-th largest
!           error estimate.
!
            levcur = Level(maxerr) + 1
            a1 = Alist(maxerr)
            b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
            a2 = b1
            b2 = Blist(maxerr)
            erlast = errmax
            CALL QK21(F,a1,b1,area1,error1,resa,defab1)
            CALL QK21(F,a2,b2,area2,error2,resa,defab2)
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
               IF ( ABS(Rlist(maxerr)-area12)<=0.1E-04*ABS(area12) .AND.&
                  & erro12>=0.99E+00*errmax ) THEN
                  IF ( extrap ) iroff2 = iroff2 + 1
                  IF ( .NOT.extrap ) iroff1 = iroff1 + 1
               ENDIF
               IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
            ENDIF
            Level(maxerr) = levcur
            Level(Last) = levcur
            Rlist(maxerr) = area1
            Rlist(Last) = area2
            errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
!
!           test for roundoff error and eventually
!           set error flag.
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
            IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach)       &
               & *(ABS(a2)+0.1E+04*uflow) ) Ier = 4
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
            CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
            IF ( errsum<=errbnd ) GOTO 150
! ***jump out of do-loop
            IF ( Ier/=0 ) GOTO 100
            IF ( .NOT.(noext) ) THEN
               erlarg = erlarg - erlast
               IF ( levcur+1<=levmax ) erlarg = erlarg + erro12
               IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
                  IF ( Level(maxerr)+1<=levmax ) GOTO 50
                  extrap = .TRUE.
                  nrmax = 2
               ENDIF
               IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors
!           over the larger intervals (erlarg) and perform
!           extrapolation.
!
                  id = nrmax
                  jupbnd = Last
                  IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
                  DO k = id , jupbnd
                     maxerr = Iord(nrmax)
                     errmax = Elist(maxerr)
! ***jump out of do-loop
                     IF ( Level(maxerr)+1<=levmax ) GOTO 50
                     nrmax = nrmax + 1
                  ENDDO
               ENDIF
!
!           perform extrapolation.
!
               numrl2 = numrl2 + 1
               rlist2(numrl2) = area
               IF ( numrl2>2 ) THEN
                  CALL QELG(numrl2,rlist2,reseps,abseps,res3la,nres)
                  ktmin = ktmin + 1
                  IF ( ktmin>5 .AND. Abserr<0.1E-02*errsum ) Ier = 5
                  IF ( abseps<Abserr ) THEN
                     ktmin = 0
                     Abserr = abseps
                     Result = reseps
                     correc = erlarg
                     ertest = AMAX1(Epsabs,Epsrel*ABS(reseps))
! ***jump out of do-loop
                     IF ( Abserr<ertest ) GOTO 100
                  ENDIF
!
!           prepare bisection of the smallest interval.
!
                  IF ( numrl2==1 ) noext = .TRUE.
                  IF ( Ier>=5 ) GOTO 100
               ENDIF
               maxerr = Iord(1)
               errmax = Elist(maxerr)
               nrmax = 1
               extrap = .FALSE.
               levmax = levmax + 1
               erlarg = errsum
            ENDIF
 50      ENDDO
!
!           set the final result.
!           ---------------------
!
!
 100     IF ( Abserr/=oflow ) THEN
            IF ( (Ier+ierro)/=0 ) THEN
               IF ( ierro==3 ) Abserr = Abserr + correc
               IF ( Ier==0 ) Ier = 3
               IF ( Result==0.0E+00 .OR. area==0.0E+00 ) THEN
                  IF ( Abserr>errsum ) GOTO 150
                  IF ( area==0.0E+00 ) GOTO 200
               ELSEIF ( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
                  GOTO 150
               ENDIF
            ENDIF
!
!           test on divergence.
!
            IF ( ksgn/=(-1) .OR. AMAX1(ABS(Result),ABS(area))           &
               & >resabs*0.1E-01 ) THEN
               IF ( 0.1E-01>(Result/area) .OR. (Result/area)            &
                  & >0.1E+03 .OR. errsum>ABS(area) ) Ier = 6
            ENDIF
            GOTO 200
         ENDIF
!
!           compute global integral sum.
!
 150     Result = 0.0E+00
         DO k = 1 , Last
            Result = Result + Rlist(k)
         ENDDO
         Abserr = errsum
      ENDIF
 200  IF ( Ier>2 ) Ier = Ier - 1
      Result = Result*sign
99999 END
!*==QAGS.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGS(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier,Limit,&
                    & Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAGS10628
!*** Start of declarations inserted by SPAG
      INTEGER Last
!*** End of declarations inserted by SPAG
!***begin prologue  qags
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
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
!        real version
!
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28)
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
!            work  - real
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
!***references  (none)
!***routines called  qagse,xerror
!***end prologue  qags
!
!
      REAL A , Abserr , B , Epsabs , Epsrel , F , Result , Work
      INTEGER Ier , Iwork , Lenw , Limit , lvl , l1 , l2 , l3 , Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limit and lenw.
!
!***first executable statement  qags
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for qagse.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL QAGSE(F,A,B,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,  &
                  & Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qags',26,Ier,lvl)
      END
!*==QAGSE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAGSE(F,A,B,Epsabs,Epsrel,Limit,Result,Abserr,Neval,   &
                     & Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE
!*--QAGSE10824
!***begin prologue  qagse
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
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
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b)
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                             epsrel<max(50*rel.mach.acc.,0.5d-28).
!                             result, abserr, neval, last, rlist(1),
!                             iord(1) and elist(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!            alist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left end points
!                     of the subintervals in the partition of the
!                     given integration range (a,b)
!
!            blist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right end points
!                     of the subintervals in the partition of the given
!                     integration range (a,b)
!
!            rlist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real
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
!***references  (none)
!***routines called  qelg,qk21,qpsrt,r1mach
!***end prologue  qagse
!
      REAL A , abseps , Abserr , Alist , area , area1 , area12 , area2 ,&
         & a1 , a2 , B , Blist , b1 , b2 , correc , defabs , defab1 ,   &
         & defab2 , R1MACH , dres , Elist , epmach , Epsabs , Epsrel ,  &
         & erlarg , erlast , errbnd , errmax , error1 , error2 ,        &
         & erro12 , errsum , ertest , F , oflow , resabs , reseps ,     &
         & Result , res3la , Rlist , rlist2 , small , uflow
      INTEGER id , Ier , ierro , Iord , iroff1 , iroff2 , iroff3 ,      &
            & jupbnd , k , ksgn , ktmin , Last , Limit , maxerr ,       &
            & Neval , nres , nrmax , numrl2
      LOGICAL extrap , noext
!
      DIMENSION Alist(Limit) , Blist(Limit) , Elist(Limit) , Iord(Limit)&
              & , res3la(3) , Rlist(Limit) , rlist2(52)
!
      EXTERNAL F
!
!            the dimension of rlist2 is determined by the value of
!            limexp in subroutine qelg (rlist2 should be of dimension
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
!           rlist2    - array of dimension at least limexp+2
!                       containing the part of the epsilon table
!                       which is still needed for further computations
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
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation
!                       i.e. before subdividing the smallest interval
!                       we try to decrease the value of erlarg.
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
!***first executable statement  qagse
      epmach = R1MACH(4)
!
!            test on validity of parameters
!            ------------------------------
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      IF ( Epsabs<=0.0E+00 .AND. Epsrel<AMAX1(0.5E+02*epmach,0.5E-14) ) &
         & Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         uflow = R1MACH(1)
         oflow = R1MACH(2)
         ierro = 0
         CALL QK21(F,A,B,Result,Abserr,defabs,resabs)
!
!           test on accuracy.
!
         dres = ABS(Result)
         errbnd = AMAX1(Epsabs,Epsrel*dres)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         IF ( Abserr<=1.0E+02*epmach*defabs .AND. Abserr>errbnd )       &
            & Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) .OR.    &
            & Abserr==0.0E+00 ) THEN
            Neval = 42*Last - 21
            GOTO 99999
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
            IF ( dres>=(0.1E+01-0.5E+02*epmach)*defabs ) ksgn = 1
!
!           main do-loop
!           ------------
!
            DO Last = 2 , Limit
!
!           bisect the subinterval with the nrmax-th largest
!           error estimate.
!
               a1 = Alist(maxerr)
               b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               erlast = errmax
               CALL QK21(F,a1,b1,area1,error1,resabs,defab1)
               CALL QK21(F,a2,b2,area2,error2,resabs,defab2)
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( defab1/=error1 .AND. defab2/=error2 ) THEN
                  IF ( ABS(Rlist(maxerr)-area12)<=0.1E-04*ABS(area12)   &
                     & .AND. erro12>=0.99E+00*errmax ) THEN
                     IF ( extrap ) iroff2 = iroff2 + 1
                     IF ( .NOT.extrap ) iroff1 = iroff1 + 1
                  ENDIF
                  IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
!
!           test for roundoff error and eventually
!           set error flag.
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
               IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach)    &
                  & *(ABS(a2)+0.1E+04*uflow) ) Ier = 4
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
               CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( errsum<=errbnd ) GOTO 50
! ***jump out of do-loop
               IF ( Ier/=0 ) GOTO 40
               IF ( Last==2 ) THEN
                  small = ABS(B-A)*0.375E+00
                  erlarg = errsum
                  ertest = errbnd
                  rlist2(2) = area
               ELSEIF ( .NOT.(noext) ) THEN
                  erlarg = erlarg - erlast
                  IF ( ABS(b1-a1)>small ) erlarg = erlarg + erro12
                  IF ( .NOT.(extrap) ) THEN
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
                     IF ( ABS(Blist(maxerr)-Alist(maxerr))>small )      &
                        & GOTO 20
                     extrap = .TRUE.
                     nrmax = 2
                  ENDIF
                  IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors
!           over the larger intervals (erlarg) and perform
!           extrapolation.
!
                     id = nrmax
                     jupbnd = Last
                     IF ( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
                     DO k = id , jupbnd
                        maxerr = Iord(nrmax)
                        errmax = Elist(maxerr)
! ***jump out of do-loop
                        IF ( ABS(Blist(maxerr)-Alist(maxerr))>small )   &
                           & GOTO 20
                        nrmax = nrmax + 1
                     ENDDO
                  ENDIF
!
!           perform extrapolation.
!
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
                  CALL QELG(numrl2,rlist2,reseps,abseps,res3la,nres)
                  ktmin = ktmin + 1
                  IF ( ktmin>5 .AND. Abserr<0.1E-02*errsum ) Ier = 5
                  IF ( abseps<Abserr ) THEN
                     ktmin = 0
                     Abserr = abseps
                     Result = reseps
                     correc = erlarg
                     ertest = AMAX1(Epsabs,Epsrel*ABS(reseps))
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
                  small = small*0.5E+00
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
                  IF ( Result==0.0E+00 .OR. area==0.0E+00 ) THEN
                     IF ( Abserr>errsum ) GOTO 50
                     IF ( area==0.0E+00 ) THEN
                        IF ( Ier>2 ) Ier = Ier - 1
                        Neval = 42*Last - 21
                        GOTO 99999
                     ENDIF
                  ELSEIF ( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
                     GOTO 50
                  ENDIF
               ENDIF
!
!           test on divergence.
!
               IF ( ksgn/=(-1) .OR. AMAX1(ABS(Result),ABS(area))        &
                  & >defabs*0.1E-01 ) THEN
                  IF ( 0.1E-01>(Result/area) .OR. (Result/area)         &
                     & >0.1E+03 .OR. errsum>ABS(area) ) Ier = 6
               ENDIF
               IF ( Ier>2 ) Ier = Ier - 1
               Neval = 42*Last - 21
               GOTO 99999
            ENDIF
         ENDIF
!
!           compute global integral sum.
!
 50      Result = 0.0E+00
         DO k = 1 , Last
            Result = Result + Rlist(k)
         ENDDO
         Abserr = errsum
         IF ( Ier>2 ) Ier = Ier - 1
         Neval = 42*Last - 21
      ENDIF
99999 END
!*==QAWC.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWC(F,A,B,C,Epsabs,Epsrel,Result,Abserr,Neval,Ier,    &
                    & Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAWC11299
!*** Start of declarations inserted by SPAG
      INTEGER Last
!*** End of declarations inserted by SPAG
!***begin prologue  qawc
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1,j4
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
!        real version
!
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     under limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            c      - parameter in the weight function, c/=a, c/=b.
!                     if c = a or c = b, the routine will end with
!                     ier = 6 .
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!            work  - real
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
!***references  (none)
!***routines called  qawce,xerror
!***end prologue  qawc
!
      REAL A , Abserr , B , C , Epsabs , Epsrel , F , Result , Work
      INTEGER Ier , Iwork , Lenw , Limit , lvl , l1 , l2 , l3 , Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limit and lenw.
!
!***first executable statement  qawc
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for qawce.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
         CALL QAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,&
                  & Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qawc',26,Ier,lvl)
      END
!*==QAWCE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval, &
                     & Ier,Alist,Blist,Rlist,Elist,Iord,Last)
      IMPLICIT NONE
!*--QAWCE11483
!*** Start of declarations inserted by SPAG
      REAL erro12
!*** End of declarations inserted by SPAG
!***begin prologue  qawce
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1,j4
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
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            c      - real
!                     parameter in the weight function, c/=a, c/=b
!                     if c = a or c = b, the routine will end with
!                     ier = 6.
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=1
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
!                             or limit<1.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero. alist(1)
!                             and blist(1) are set to a and b
!                             respectively.
!
!            alist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            blist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            rlist   - real
!                      vector of dimension at least limit, the first
!                       last  elements of which are the integral
!                      approximations on the subintervals
!
!            elist   - real
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
!***references  (none)
!***routines called  qc25c,qpsrt,r1mach
!***end prologue  qawce
!
      REAL A , aa , Abserr , Alist , area , area1 , area12 , area2 ,    &
         & a1 , a2 , B , bb , Blist , b1 , b2 , C , R1MACH , Elist ,    &
         & epmach , Epsabs , Epsrel , errbnd , errmax , error1 ,        &
         & error2 , errsum , F , Result , Rlist , uflow
      INTEGER Ier , Iord , iroff1 , iroff2 , k , krule , Last , Limit , &
            & maxerr , nev , Neval , nrmax
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) ,            &
              & Elist(Limit) , Iord(Limit)
!
      EXTERNAL F
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
!***first executable statement  qawce
      epmach = R1MACH(4)
      uflow = R1MACH(1)
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
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( .NOT.(C==A .OR. C==B .OR. (Epsabs<=0.0E+00 .AND. Epsrel<AMAX1&
         & (0.5E+02*epmach,0.5E-14))) ) THEN
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
         CALL QC25C(F,aa,bb,C,Result,Abserr,krule,Neval)
         Last = 1
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         Alist(1) = A
         Blist(1) = B
!
!           test on accuracy
!
         errbnd = AMAX1(Epsabs,Epsrel*ABS(Result))
         IF ( Limit==1 ) Ier = 1
         IF ( Abserr>=AMIN1(0.1E-01*ABS(Result),errbnd) .AND. Ier/=1 )  &
            & THEN
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
               b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
               b2 = Blist(maxerr)
               IF ( C<=b1 .AND. C>a1 ) b1 = 0.5E+00*(C+b2)
               IF ( C>b1 .AND. C<b2 ) b1 = 0.5E+00*(a1+C)
               a2 = b1
               krule = 2
               CALL QC25C(F,a1,b1,C,area1,error1,krule,nev)
               Neval = Neval + nev
               CALL QC25C(F,a2,b2,C,area2,error2,krule,nev)
               Neval = Neval + nev
!
!           improve previous approximations to integral
!           and error and test for accuracy.
!
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - Rlist(maxerr)
               IF ( ABS(Rlist(maxerr)-area12)<0.1E-04*ABS(area12) .AND. &
                  & erro12>=0.99E+00*errmax .AND. krule==0 )            &
                  & iroff1 = iroff1 + 1
               IF ( Last>10 .AND. erro12>errmax .AND. krule==0 )        &
                  & iroff2 = iroff2 + 1
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
               IF ( errsum>errbnd ) THEN
!
!           test for roundoff error and eventually
!           set error flag.
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
                  IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach) &
                     & *(ABS(a2)+0.1E+04*uflow) ) Ier = 3
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
               CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0E+00
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
         IF ( aa==B ) Result = -Result
      ENDIF
      END
!*==QAWF.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWF(F,A,Omega,Integr,Epsabs,Result,Abserr,Neval,Ier,  &
                    & Limlst,Lst,Leniw,Maxp1,Lenw,Iwork,Work)
      IMPLICIT NONE
!*--QAWF11822
!*** Start of declarations inserted by SPAG
      INTEGER Iwork , last , Lenw , ll2
!*** End of declarations inserted by SPAG
!***begin prologue  qawf
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a1
!***keywords  automatic integrator, special-purpose,fourier
!             integral, integration between zeros with dqawoe,
!             convergence acceleration with dqext
!***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            fourier integral
!            i = integral of f(x)*w(x) over (a,infinity)
!            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=epsabs.
!***description
!
!        computation of fourier integrals
!        standard fortran subroutine
!        real version
!
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            omega  - real
!                     parameter in the integrand weight function
!
!            integr - integer
!                     indicates which of the weight functions is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine
!                     will end with ier = 6.
!
!            epsabs - real
!                     absolute accuracy requested, epsabs>0.
!                     if epsabs<=0, the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!            work   - real
!                     vector of dimension at least
!                     on return,
!                     work(1), ..., work(lst) contain the integral
!                      approximations over the cycles,
!                     work(limlst+1), ..., work(limlst+lst) contain
!                      the error extimates over the cycles.
!                     further elements of work have no specific
!                     meaning for the user.
!
!***references  (none)
!***routines called  qawfe,xerror
!***end prologue  qawf
!
      REAL A , Abserr , Epsabs , F , Omega , Result , Work
      INTEGER Ier , Integr , Leniw , limit , Limlst , lvl , Lst , l1 ,  &
            & l2 , l3 , l4 , l5 , l6 , Maxp1 , Neval
!
      DIMENSION Iwork(Leniw) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limlst, leniw, maxp1 and lenw.
!
!***first executable statement  qawf
      Ier = 6
      Neval = 0
      last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limlst>=3 .AND. Leniw>=(Limlst+2) .AND. Maxp1>=1 .AND.       &
         & Lenw>=(Leniw*2+Maxp1*25) ) THEN
!
!         prepare call for qawfe
!
         limit = (Leniw-Limlst)/2
         l1 = Limlst + 1
         l2 = Limlst + l1
         l3 = limit + l2
         l4 = limit + l3
         l5 = limit + l4
         l6 = limit + l5
         ll2 = limit + l1
         CALL QAWFE(F,A,Omega,Integr,Epsabs,Limlst,limit,Maxp1,Result,  &
                  & Abserr,Neval,Ier,Work(1),Work(l1),Iwork(1),Lst,     &
                  & Work(l2),Work(l3),Work(l4),Work(l5),Iwork(l1),      &
                  & Iwork(ll2),Work(l6))
!
!         call error handler if necessary
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qawf',26,Ier,lvl)
      END
!*==QAWFE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWFE(F,A,Omega,Integr,Epsabs,Limlst,Limit,Maxp1,      &
                     & Result,Abserr,Neval,Ier,Rslst,Erlst,Ierlst,Lst,  &
                     & Alist,Blist,Rlist,Elist,Iord,Nnlog,Chebmo)
      IMPLICIT NONE
!*--QAWFE12062
!*** Start of declarations inserted by SPAG
      REAL F
      INTEGER last , Limlst , momcom
!*** End of declarations inserted by SPAG
!***begin prologue  qawfe
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a1
!***keywords  automatic integrator, special-purpose,
!             fourier integrals,
!             integration between zeros with dqawoe,
!             convergence acceleration with dqelg
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a
!            given fourier integal
!            i = integral of f(x)*w(x) over (a,infinity)
!             where w(x) = cos(omega*x) or w(x) = sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=epsabs.
!***description
!
!        computation of fourier integrals
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to
!                     be declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            omega  - real
!                     parameter in the weight function
!
!            integr - integer
!                     indicates which weight function is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine will
!                     end with ier = 6.
!
!            epsabs - real
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
!            result - real
!                     approximation to the integral x
!
!            abserr - real
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
!            rslst  - real
!                     vector of dimension at least limlst
!                     rslst(k) contains the integral contribution
!                     over the interval (a+(k-1)c,a+kc) where
!                     c = (2*int(abs(omega))+1)*pi/abs(omega),
!                     k = 1, 2, ..., lst.
!                     note that, if omega = 0, rslst(1) contains
!                     the value of the integral over (a,infinity).
!
!            erlst  - real
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
!            alist, blist, rlist, elist - real
!                     vector of dimension at least limit,
!
!            iord, nnlog - integer
!                     vector of dimension at least limit, providing
!                     space for the quantities needed in the subdivision
!                     process of each cycle
!
!            chebmo - real
!                     array of dimension at least (maxp1,25), providing
!                     space for the chebyshev moments needed within the
!                     cycles
!
!***references  (none)
!***routines called  qagie,qawoe,qelg,r1mach
!***end prologue  qawfe
!
      REAL A , abseps , Abserr , Alist , Blist , Chebmo , correc ,      &
         & cycle , c1 , c2 , dl , dla , drl , Elist , ep , eps , epsa , &
         & Epsabs , Erlst , errsum , fact , Omega , p , pi , p1 , psum ,&
         & reseps , Result , res3la , Rlist , Rslst , R1MACH , uflow
      INTEGER Ier , Ierlst , Integr , Iord , ktmin , l , Lst , Limit ,  &
            & ll , Maxp1 , nev , Neval , Nnlog , nres , numrl2
!
      DIMENSION Alist(Limit) , Blist(Limit) , Chebmo(Maxp1,25) ,        &
              & Elist(Limit) , Erlst(Limlst) , Ierlst(Limlst) ,         &
              & Iord(Limit) , Nnlog(Limit) , psum(52) , res3la(3) ,     &
              & Rlist(Limit) , Rslst(Limlst)
!
      EXTERNAL F
!
!
!            the dimension of  psum  is determined by the value of
!            limexp in subroutine qelg (psum must be
!            of dimension (limexp+2) at least).
!
!           list of major variables
!           -----------------------
!
!           c1, c2    - end points of subinterval (of length
!                       cycle)
!           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
!           psum      - vector of dimension at least (limexp+2)
!                       (see routine qelg)
!                       psum contains the part of the epsilon
!                       table which is still needed for further
!                       computations.
!                       each element of psum is a partial sum of
!                       the series which should sum to the value of
!                       the integral.
!           errsum    - sum of error estimates over the
!                       subintervals, calculated cumulatively
!           epsa      - absolute tolerance requested over current
!                       subinterval
!           chebmo    - array containing the modified chebyshev
!                       moments (see also routine qc25f)
!
      DATA p/0.9E+00/ , pi/0.31415926535897932E+01/
!
!           test on validity of parameters
!           ------------------------------
!
!***first executable statement  qawfe
      Result = 0.0E+00
      Abserr = 0.0E+00
      Neval = 0
      Lst = 0
      Ier = 0
      IF ( (Integr/=1 .AND. Integr/=2) .OR. Epsabs<=0.0E+00 .OR.        &
         & Limlst<3 ) Ier = 6
      IF ( Ier/=6 ) THEN
         IF ( Omega/=0.0E+00 ) THEN
!
!           initializations
!           ---------------
!
            l = ABS(Omega)
            dl = 2*l + 1
            cycle = dl*pi/ABS(Omega)
            Ier = 0
            ktmin = 0
            Neval = 0
            numrl2 = 0
            nres = 0
            c1 = A
            c2 = cycle + A
            p1 = 0.1E+01 - p
            eps = Epsabs
            uflow = R1MACH(1)
            IF ( Epsabs>uflow/p1 ) eps = Epsabs*p1
            ep = eps
            fact = 0.1E+01
            correc = 0.0E+00
            Abserr = 0.0E+00
            errsum = 0.0E+00
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
               CALL QAWOE(F,c1,c2,Omega,Integr,epsa,0.0E+00,Limit,Lst,  &
                        & Maxp1,Rslst(Lst),Erlst(Lst),nev,Ierlst(Lst),  &
                        & last,Alist,Blist,Rlist,Elist,Iord,Nnlog,      &
                        & momcom,Chebmo)
               Neval = Neval + nev
               fact = fact*p
               errsum = errsum + Erlst(Lst)
               drl = 0.5E+02*ABS(Rslst(Lst))
!
!           test on accuracy with partial sum
!
               IF ( errsum+drl<=Epsabs .AND. Lst>=6 ) GOTO 50
               correc = AMAX1(correc,Erlst(Lst))
               IF ( Ierlst(Lst)/=0 ) eps = AMAX1(ep,correc*p1)
               IF ( Ierlst(Lst)/=0 ) Ier = 7
               IF ( Ier==7 .AND. (errsum+drl)<=correc*0.1E+02 .AND.     &
                  & Lst>5 ) GOTO 50
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
                     CALL QELG(numrl2,psum,reseps,abseps,res3la,nres)
!
!           test whether extrapolated result is influenced by
!           roundoff
!
                     ktmin = ktmin + 1
                     IF ( ktmin>=15 .AND. Abserr<=0.1E-02*(errsum+drl) )&
                        & Ier = 4
                     IF ( abseps<=Abserr .OR. Lst==3 ) THEN
                        Abserr = abseps
                        Result = reseps
                        ktmin = 0
!
!           if ier is not 0, check whether direct result (partial
!           sum) or extrapolated result yields the best integral
!           approximation
!
                        IF ( (Abserr+0.1E+02*correc)<=Epsabs .OR.       &
                           & (Abserr<=Epsabs .AND.                      &
                           & 0.1E+02*correc>=Epsabs) ) GOTO 20
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
 20         Abserr = Abserr + 0.1E+02*correc
            IF ( Ier==0 ) GOTO 99999
            IF ( Result==0.0E+00 .OR. psum(numrl2)==0.0E+00 ) THEN
               IF ( Abserr>errsum ) GOTO 50
               IF ( psum(numrl2)==0.0E+00 ) GOTO 99999
            ENDIF
            IF ( Abserr/ABS(Result)<=(errsum+drl)/ABS(psum(numrl2)) )   &
               & THEN
               IF ( Ier>=1 .AND. Ier/=7 ) Abserr = Abserr + drl
               GOTO 99999
            ENDIF
         ELSE
!
!           integration by qagie if omega is zero
!           --------------------------------------
!
            IF ( Integr==1 ) CALL QAGIE(F,0.0E+00,1,Epsabs,0.0E+00,     &
                                      & Limit,Result,Abserr,Neval,Ier,  &
                                      & Alist,Blist,Rlist,Elist,Iord,   &
                                      & last)
            Rslst(1) = Result
            Erlst(1) = Abserr
            Ierlst(1) = Ier
            Lst = 1
            GOTO 99999
         ENDIF
 50      Result = psum(numrl2)
         Abserr = errsum + drl
      ENDIF
99999 END
!*==QAWO.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWO(F,A,B,Omega,Integr,Epsabs,Epsrel,Result,Abserr,   &
                    & Neval,Ier,Leniw,Maxp1,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAWO12446
!*** Start of declarations inserted by SPAG
      INTEGER Iwork , Last , Lenw , limit
      REAL Work
!*** End of declarations inserted by SPAG
!***begin prologue  qawo
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
!***keywords  automatic integrator, special-purpose,
!             integrand with oscillatory cos or sin factor,
!             clenshaw-curtis method, (end point) singularities,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral
!            i = integral of f(x)*w(x) over (a,b)
!            where w(x) = cos(omega*x) or w(x) = sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of oscillatory integrals
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the function
!                     f(x).  the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            omega  - real
!                     parameter in the integrand weight function
!
!            integr - integer
!                     indicates which of the weight functions is used
!                     integr = 1      w(x) = cos(omega*x)
!                     integr = 2      w(x) = sin(omega*x)
!                     if integr/=1.and.integr/=2, the routine will
!                     end with ier = 6.
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if epsabs<=0 and
!                     epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!            work   - real
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
!***references  (none)
!***routines called  qawoe,xerror
!***end prologue  qawo
!
      REAL A , Abserr , B , Epsabs , Epsrel , F , Omega , Result
      INTEGER Ier , Integr , Leniw , lvl , l1 , l2 , l3 , l4 , Maxp1 ,  &
            & momcom , Neval
!
      DIMENSION Iwork(Leniw) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of leniw, maxp1 and lenw.
!
!***first executable statement  qawo
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Leniw>=2 .AND. Maxp1>=1 .AND. Lenw>=(Leniw*2+Maxp1*25) ) THEN
!
!         prepare call for qawoe
!
         limit = Leniw/2
         l1 = limit + 1
         l2 = limit + l1
         l3 = limit + l2
         l4 = limit + l3
         CALL QAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,limit,1,Maxp1,     &
                  & Result,Abserr,Neval,Ier,Last,Work(1),Work(l1),      &
                  & Work(l2),Work(l3),Iwork(1),Iwork(l1),momcom,Work(l4)&
                  & )
!
!         call error handler if necessary
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qawo',26,Ier,lvl)
      END
!*==QAWOE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,Limit,Icall,    &
                     & Maxp1,Result,Abserr,Neval,Ier,Last,Alist,Blist,  &
                     & Rlist,Elist,Iord,Nnlog,Momcom,Chebmo)
      IMPLICIT NONE
!*--QAWOE12679
!***begin prologue  qawoe
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
!***keywords  automatic integrator, special-purpose,
!             integrand with oscillatory cos or sin factor,
!             clenshaw-curtis method, (end point) singularities,
!             extrapolation, globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral
!            i = integral of f(x)*w(x) over (a,b)
!            where w(x) = cos(omega*x) or w(x) = sin(omega*x),
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of oscillatory integrals
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration
!
!            omega  - real
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
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
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
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
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
!            alist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            blist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            rlist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real
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
!            chebmo - real
!                     array of dimension (maxp1,25) containing the
!                     chebyshev moments
!
!***references  (none)
!***routines called  qc25f,qelg,qpsrt,r1mach
!***end prologue  qawoe
!
      REAL A , abseps , Abserr , Alist , area , area1 , area12 , area2 ,&
         & a1 , a2 , B , Blist , b1 , b2 , Chebmo , correc , defab1 ,   &
         & defab2 , defabs , domega , R1MACH , dres , Elist , epmach ,  &
         & Epsabs , Epsrel , erlarg , erlast , errbnd , errmax ,        &
         & error1 , erro12 , error2 , errsum , ertest , F , oflow ,     &
         & Omega , resabs , reseps , Result , res3la , Rlist , rlist2 , &
         & small , uflow , width
      INTEGER Icall , id , Ier , ierro , Integr , Iord , iroff1 ,       &
            & iroff2 , iroff3 , jupbnd , k , ksgn , ktmin , Last ,      &
            & Limit , maxerr , Maxp1 , Momcom , nev , Neval , Nnlog ,   &
            & nres , nrmax , nrmom , numrl2
      LOGICAL extrap , noext , extall
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) ,            &
              & Elist(Limit) , Iord(Limit) , rlist2(52) , res3la(3) ,   &
              & Chebmo(Maxp1,25) , Nnlog(Limit)
!
      EXTERNAL F
!
!            the dimension of rlist2 is determined by  the value of
!            limexp in subroutine qelg (rlist2 should be of
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
!                       is no longer allowed (true value)
!
!            machine dependent constants
!            ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!           oflow is the largest positive magnitude.
!
!***first executable statement  qawoe
      epmach = R1MACH(4)
!
!         test on validity of parameters
!         ------------------------------
!
      Ier = 0
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      Alist(1) = A
      Blist(1) = B
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      Nnlog(1) = 0
      IF ( (Integr/=1 .AND. Integr/=2) .OR.                             &
         & (Epsabs<=0.0E+00 .AND. Epsrel<AMAX1(0.5E+02*epmach,0.5E-14)) &
         & .OR. Icall<1 .OR. Maxp1<1 ) Ier = 6
      IF ( Ier/=6 ) THEN
!
!           first approximation to the integral
!           -----------------------------------
!
         domega = ABS(Omega)
         nrmom = 0
         IF ( Icall<=1 ) Momcom = 0
         CALL QC25F(F,A,B,domega,Integr,nrmom,Maxp1,0,Result,Abserr,    &
                  & Neval,defabs,resabs,Momcom,Chebmo)
!
!           test on accuracy.
!
         dres = ABS(Result)
         errbnd = AMAX1(Epsabs,Epsrel*dres)
         Rlist(1) = Result
         Elist(1) = Abserr
         Iord(1) = 1
         IF ( Abserr<=0.1E+03*epmach*defabs .AND. Abserr>errbnd )       &
            & Ier = 2
         IF ( Limit==1 ) Ier = 1
         IF ( Ier/=0 .OR. Abserr<=errbnd ) THEN
            IF ( Integr==2 .AND. Omega<0.0E+00 ) Result = -Result
            GOTO 99999
         ELSE
!
!           initializations
!           ---------------
!
            uflow = R1MACH(1)
            oflow = R1MACH(2)
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
            small = ABS(B-A)*0.75E+00
            nres = 0
            numrl2 = 0
            extall = .FALSE.
            IF ( 0.5E+00*ABS(B-A)*domega<=0.2E+01 ) THEN
               numrl2 = 1
               extall = .TRUE.
               rlist2(1) = Result
            ENDIF
            IF ( 0.25E+00*ABS(B-A)*domega<=0.2E+01 ) extall = .TRUE.
            ksgn = -1
            IF ( dres>=(0.1E+01-0.5E+02*epmach)*defabs ) ksgn = 1
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
               b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
               erlast = errmax
               CALL QC25F(F,a1,b1,domega,Integr,nrmom,Maxp1,0,area1,    &
                        & error1,nev,resabs,defab1,Momcom,Chebmo)
               Neval = Neval + nev
               CALL QC25F(F,a2,b2,domega,Integr,nrmom,Maxp1,1,area2,    &
                        & error2,nev,resabs,defab2,Momcom,Chebmo)
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
                  IF ( ABS(Rlist(maxerr)-area12)<=0.1E-04*ABS(area12)   &
                     & .AND. erro12>=0.99E+00*errmax ) THEN
                     IF ( extrap ) iroff2 = iroff2 + 1
                     IF ( .NOT.extrap ) iroff1 = iroff1 + 1
                  ENDIF
                  IF ( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
               Nnlog(maxerr) = nrmom
               Nnlog(Last) = nrmom
               errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
!
!           test for roundoff error and eventually
!           set error flag
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
               IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach)    &
                  & *(ABS(a2)+0.1E+04*uflow) ) Ier = 4
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
               CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( errsum<=errbnd ) GOTO 50
               IF ( Ier/=0 ) GOTO 40
               IF ( Last==2 .AND. extall ) THEN
                  small = small*0.5E+00
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
               ELSE
                  IF ( noext ) GOTO 20
                  IF ( extall ) THEN
                     erlarg = erlarg - erlast
                     IF ( ABS(b1-a1)>small ) erlarg = erlarg + erro12
                     IF ( extrap ) GOTO 5
                  ENDIF
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
                  width = ABS(Blist(maxerr)-Alist(maxerr))
                  IF ( width>small ) GOTO 20
                  IF ( extall ) THEN
                     extrap = .TRUE.
                     nrmax = 2
                  ELSE
!
!           test whether we can start with the extrapolation
!           procedure (we do this if we integrate over the
!           next interval with use of a gauss-kronrod rule - see
!           subroutine qc25f).
!
                     small = small*0.5E+00
                     IF ( 0.25E+00*width*domega>0.2E+01 ) GOTO 20
                     extall = .TRUE.
                     GOTO 10
                  ENDIF
 5                IF ( ierro/=3 .AND. erlarg>ertest ) THEN
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors
!           over the larger intervals (erlarg) and perform
!           extrapolation.
!
                     jupbnd = Last
                     IF ( Last>(Limit/2+2) ) jupbnd = Limit + 3 - Last
                     id = nrmax
                     DO k = id , jupbnd
                        maxerr = Iord(nrmax)
                        errmax = Elist(maxerr)
                        IF ( ABS(Blist(maxerr)-Alist(maxerr))>small )   &
                           & GOTO 20
                        nrmax = nrmax + 1
                     ENDDO
                  ENDIF
!
!           perform extrapolation.
!
                  numrl2 = numrl2 + 1
                  rlist2(numrl2) = area
                  IF ( numrl2>=3 ) THEN
                     CALL QELG(numrl2,rlist2,reseps,abseps,res3la,nres)
                     ktmin = ktmin + 1
                     IF ( ktmin>5 .AND. Abserr<0.1E-02*errsum ) Ier = 5
                     IF ( abseps<Abserr ) THEN
                        ktmin = 0
                        Abserr = abseps
                        Result = reseps
                        correc = erlarg
                        ertest = AMAX1(Epsabs,Epsrel*ABS(reseps))
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
                  small = small*0.5E+00
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
                  IF ( Result==0.0E+00 .OR. area==0.0E+00 ) THEN
                     IF ( Abserr>errsum ) GOTO 50
                     IF ( area==0.0E+00 ) THEN
                        IF ( Ier>2 ) Ier = Ier - 1
                        IF ( Integr==2 .AND. Omega<0.0E+00 ) Result = - &
                           & Result
                        GOTO 99999
                     ENDIF
                  ELSEIF ( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
                     GOTO 50
                  ENDIF
               ENDIF
!
!           test on divergence.
!
               IF ( ksgn/=(-1) .OR. AMAX1(ABS(Result),ABS(area))        &
                  & >defabs*0.1E-01 ) THEN
                  IF ( 0.1E-01>(Result/area) .OR. (Result/area)         &
                     & >0.1E+03 .OR. errsum>=ABS(area) ) Ier = 6
               ENDIF
               IF ( Ier>2 ) Ier = Ier - 1
               IF ( Integr==2 .AND. Omega<0.0E+00 ) Result = -Result
               GOTO 99999
            ENDIF
         ENDIF
!
!           compute global integral sum.
!
 50      Result = 0.0E+00
         DO k = 1 , Last
            Result = Result + Rlist(k)
         ENDDO
         Abserr = errsum
         IF ( Ier>2 ) Ier = Ier - 1
         IF ( Integr==2 .AND. Omega<0.0E+00 ) Result = -Result
      ENDIF
99999 END
!*==QAWS.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWS(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Result,      &
                    & Abserr,Neval,Ier,Limit,Lenw,Last,Iwork,Work)
      IMPLICIT NONE
!*--QAWS13248
!*** Start of declarations inserted by SPAG
      INTEGER Last
!*** End of declarations inserted by SPAG
!***begin prologue  qaws
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
!***keywords  automatic integrator, special-purpose,
!             algebraico-logarithmic end-point singularities,
!             clenshaw-curtis, globally adaptive
!***author  piessens,robert,appl. math. & progr. div. -k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral i = integral of f*w over (a,b),
!            (where w shows a singular behaviour at the end points
!            see parameter integr).
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!***description
!
!        integration of functions having algebraico-logarithmic
!        end point singularities
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration, b>a
!                     if b<=a, the routine will end with ier = 6.
!
!            alfa   - real
!                     parameter in the integrand function, alfa>(-1)
!                     if alfa<=(-1), the routine will end with
!                     ier = 6.
!
!            beta   - real
!                     parameter in the integrand function, beta>(-1)
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
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand, in order to
!                             determine the integration difficulties
!                             which prevent the requested tolerance from
!                             being achieved. in case of a jump
!                             discontinuity or a local singularity
!                             of algebraico-logarithmic type at one or
!                             more interior points of the integration
!                             range, one should proceed by splitting up
!                             the interval at these points and calling
!                             the integrator on the subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b<=a or alfa<=(-1) or beta<=(-1) or
!                             or integr<1 or integr>4 or
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
!                             or limit<2 or lenw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero. except when lenw or limit is invalid
!                             iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit  - integer
!                     dimensioning parameter for iwork
!                     limit determines the maximum number of
!                     subintervals in the partition of the given
!                     integration interval (a,b), limit>=2.
!                     if limit<2, the routine will end with ier = 6.
!
!            lenw   - integer
!                     dimensioning parameter for work
!                     lenw must be at least limit*4.
!                     if lenw<limit*4, the routine will end
!                     with ier = 6.
!
!            last   - integer
!                     on return, last equals the number of
!                     subintervals produced in the subdivision process,
!                     which determines the significant number of
!                     elements actually in the work arrays.
!
!         work arrays
!            iwork  - integer
!                     vector of dimension limit, the first k
!                     elements of which contain pointers
!                     to the error estimates over the subintervals,
!                     such that work(limit*3+iwork(1)), ...,
!                     work(limit*3+iwork(k)) form a decreasing
!                     sequence with k = last if last<=(limit/2+2),
!                     and k = limit+1-last otherwise
!
!            work   - real
!                     vector of dimension lenw
!                     on return
!                     work(1), ..., work(last) contain the left
!                      end points of the subintervals in the
!                      partition of (a,b),
!                     work(limit+1), ..., work(limit+last) contain
!                      the right end points,
!                     work(limit*2+1), ..., work(limit*2+last)
!                      contain the integral approximations over
!                      the subintervals,
!                     work(limit*3+1), ..., work(limit*3+last)
!                      contain the error estimates.
!
!***references  (none)
!***routines called  qawse,xerror
!***end prologue  qaws
!
      REAL A , Abserr , Alfa , B , Beta , Epsabs , Epsrel , F , Result ,&
         & Work
      INTEGER Ier , Integr , Iwork , Lenw , Limit , lvl , l1 , l2 , l3 ,&
            & Neval
!
      DIMENSION Iwork(Limit) , Work(Lenw)
!
      EXTERNAL F
!
!         check validity of limit and lenw.
!
!***first executable statement  qaws
      Ier = 6
      Neval = 0
      Last = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( Limit>=2 .AND. Lenw>=Limit*4 ) THEN
!
!         prepare call for qawse.
!
         l1 = Limit + 1
         l2 = Limit + l1
         l3 = Limit + l2
!
         CALL QAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit,Result,  &
                  & Abserr,Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3),&
                  & Iwork,Last)
!
!         call error handler if necessary.
!
         lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERROR('abnormal return from  qaws',26,Ier,lvl)
      END
!*==QAWSE.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit,      &
                     & Result,Abserr,Neval,Ier,Alist,Blist,Rlist,Elist, &
                     & Iord,Last)
      IMPLICIT NONE
!*--QAWSE13458
!***begin prologue  qawse
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1
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
!        real version
!
!        parameters
!         on entry
!            f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real
!                     lower limit of integration
!
!            b      - real
!                     upper limit of integration, b>a
!                     if b<=a, the routine will end with ier = 6.
!
!            alfa   - real
!                     parameter in the weight function, alfa>(-1)
!                     if alfa<=(-1), the routine will end with
!                     ier = 6.
!
!            beta   - real
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
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit>=2
!                     if limit<2, the routine will end with ier = 6.
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
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
!                              epsrel<max(50*rel.mach.acc.,0.5d-28),
!                             or limit<2.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero. alist(1)
!                             and blist(1) are set to a and b
!                             respectively.
!
!            alist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            blist  - real
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (a,b)
!
!            rlist  - real
!                     vector of dimension at least limit,the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real
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
!***references  (none)
!***routines called  qc25s,qmomo,qpsrt,r1mach
!***end prologue  qawse
!
      REAL A , Abserr , Alfa , Alist , area , area1 , area12 , area2 ,  &
         & a1 , a2 , B , Beta , Blist , b1 , b2 , centre , R1MACH ,     &
         & Elist , epmach , Epsabs , Epsrel , errbnd , errmax , error1 ,&
         & erro12 , error2 , errsum , F , resas1 , resas2 , Result ,    &
         & rg , rh , ri , rj , Rlist , uflow
      INTEGER Ier , Integr , Iord , iroff1 , iroff2 , k , Last , Limit ,&
            & maxerr , nev , Neval , nrmax
!
      EXTERNAL F
!
      DIMENSION Alist(Limit) , Blist(Limit) , Rlist(Limit) ,            &
              & Elist(Limit) , Iord(Limit) , ri(25) , rj(25) , rh(25) , &
              & rg(25)
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
!***first executable statement  qawse
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Ier = 6
      Neval = 0
      Last = 0
      Rlist(1) = 0.0E+00
      Elist(1) = 0.0E+00
      Iord(1) = 0
      Result = 0.0E+00
      Abserr = 0.0E+00
      IF ( .NOT.(B<=A .OR. (Epsabs==0.0E+00 .AND. Epsrel<AMAX1(0.5E+02* &
         & epmach,0.5E-14)) .OR. Alfa<=(-0.1E+01) .OR. Beta<=(-0.1E+01) &
         & .OR. Integr<1 .OR. Integr>4 .OR. Limit<2) ) THEN
         Ier = 0
!
!           compute the modified chebyshev moments.
!
         CALL QMOMO(Alfa,Beta,ri,rj,rg,rh,Integr)
!
!           integrate over the intervals (a,(a+b)/2)
!           and ((a+b)/2,b).
!
         centre = 0.5E+00*(B+A)
         CALL QC25S(F,A,B,A,centre,Alfa,Beta,ri,rj,rg,rh,area1,error1,  &
                  & resas1,Integr,nev)
         Neval = nev
         CALL QC25S(F,A,B,centre,B,Alfa,Beta,ri,rj,rg,rh,area2,error2,  &
                  & resas2,Integr,nev)
         Last = 2
         Neval = Neval + nev
         Result = area1 + area2
         Abserr = error1 + error2
!
!           test on accuracy.
!
         errbnd = AMAX1(Epsabs,Epsrel*ABS(Result))
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
               b1 = 0.5E+00*(Alist(maxerr)+Blist(maxerr))
               a2 = b1
               b2 = Blist(maxerr)
!
               CALL QC25S(F,A,B,a1,b1,Alfa,Beta,ri,rj,rg,rh,area1,      &
                        & error1,resas1,Integr,nev)
               Neval = Neval + nev
               CALL QC25S(F,A,B,a2,b2,Alfa,Beta,ri,rj,rg,rh,area2,      &
                        & error2,resas2,Integr,nev)
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
                     IF ( ABS(Rlist(maxerr)-area12)<0.1E-04*ABS(area12) &
                        & .AND. erro12>=0.99E+00*errmax )               &
                        & iroff1 = iroff1 + 1
                     IF ( Last>10 .AND. erro12>errmax )                 &
                        & iroff2 = iroff2 + 1
                  ENDIF
               ENDIF
               Rlist(maxerr) = area1
               Rlist(Last) = area2
!
!           test on accuracy.
!
               errbnd = AMAX1(Epsabs,Epsrel*ABS(area))
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
                  IF ( AMAX1(ABS(a1),ABS(b2))<=(0.1E+01+0.1E+03*epmach) &
                     & *(ABS(a2)+0.1E+04*uflow) ) Ier = 3
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
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with largest error estimate (to be
!           bisected next).
!
               CALL QPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
! ***jump out of do-loop
               IF ( Ier/=0 .OR. errsum<=errbnd ) GOTO 20
            ENDDO
!
!           compute final result.
!           ---------------------
!
 20         Result = 0.0E+00
            DO k = 1 , Last
               Result = Result + Rlist(k)
            ENDDO
            Abserr = errsum
         ENDIF
      ENDIF
      END
!*==QC25C.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QC25C(F,A,B,C,Result,Abserr,Krul,Neval)
      IMPLICIT NONE
!*--QC25C13840
!***begin prologue  qc25c
!***date written   810101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a2,j4
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
!        real version
!
!        parameters
!           f      - real
!                    function subprogram defining the integrand function
!                    f(x). the actual name for f needs to be declared
!                    e x t e r n a l  in the driver program.
!
!           a      - real
!                    left end point of the integration interval
!
!           b      - real
!                    right end point of the integration interval, b>a
!
!           c      - real
!                    parameter in the weight function
!
!           result - real
!                    approximation to the integral
!                    result is computed by using a generalized
!                    clenshaw-curtis method if c lies within ten percent
!                    of the integration interval. in the other case the
!                    15-point kronrod rule obtained by optimal addition
!                    of abscissae to the 7-point gauss rule, is applied.
!
!           abserr - real
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
!***references  (none)
!***routines called  qcheb,qk15w,qwgtc
!***end prologue  qc25c
!
      REAL A , Abserr , ak22 , amom0 , amom1 , amom2 , B , C , cc ,     &
         & centr , cheb12 , cheb24 , QWGTC , F , fval , hlgth , p2 ,    &
         & p3 , p4 , resabs , resasc , Result , res12 , res24 , u , x
      INTEGER i , isym , k , kp , Krul , Neval
!
      DIMENSION x(11) , fval(25) , cheb12(13) , cheb24(25)
!
      EXTERNAL F , QWGTC
!
!           the vector x contains the values cos(k*pi/24),
!           k = 1, ..., 11, to be used for the chebyshev series
!           expansion of f
!
      DATA x(1) , x(2) , x(3) , x(4) , x(5) , x(6) , x(7) , x(8) ,      &
         & x(9) , x(10) , x(11)/0.9914448613738104E+00 ,                &
         & 0.9659258262890683E+00 , 0.9238795325112868E+00 ,            &
         & 0.8660254037844386E+00 , 0.7933533402912352E+00 ,            &
         & 0.7071067811865475E+00 , 0.6087614290087206E+00 ,            &
         & 0.5000000000000000E+00 , 0.3826834323650898E+00 ,            &
         & 0.2588190451025208E+00 , 0.1305261922200516E+00/
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
!           qwgtc - external function subprogram defining
!                    the weight function
!           hlgth  - half-length of the interval
!           centr  - mid point of the interval
!
!
!           check the position of c.
!
!***first executable statement  qc25c
      cc = (0.2E+01*C-B-A)/(B-A)
      IF ( ABS(cc)<0.11E+01 ) THEN
!
!           use the generalized clenshaw-curtis method.
!
         hlgth = 0.5E+00*(B-A)
         centr = 0.5E+00*(B+A)
         Neval = 25
         fval(1) = 0.5E+00*F(hlgth+centr)
         fval(13) = F(centr)
         fval(25) = 0.5E+00*F(centr-hlgth)
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)
            fval(isym) = F(centr-u)
         ENDDO
!
!           compute the chebyshev series expansion.
!
         CALL QCHEB(x,fval,cheb12,cheb24)
!
!           the modified chebyshev moments are computed
!           by forward recursion, using amom0 and amom1
!           as starting values.
!
         amom0 = ALOG(ABS((0.1E+01-cc)/(0.1E+01+cc)))
         amom1 = 0.2E+01 + cc*amom0
         res12 = cheb12(1)*amom0 + cheb12(2)*amom1
         res24 = cheb24(1)*amom0 + cheb24(2)*amom1
         DO k = 3 , 13
            amom2 = 0.2E+01*cc*amom1 - amom0
            ak22 = (k-2)*(k-2)
            IF ( (k/2)*2==k ) amom2 = amom2 - 0.4E+01/(ak22-0.1E+01)
            res12 = res12 + cheb12(k)*amom2
            res24 = res24 + cheb24(k)*amom2
            amom0 = amom1
            amom1 = amom2
         ENDDO
         DO k = 14 , 25
            amom2 = 0.2E+01*cc*amom1 - amom0
            ak22 = (k-2)*(k-2)
            IF ( (k/2)*2==k ) amom2 = amom2 - 0.4E+01/(ak22-0.1E+01)
            res24 = res24 + cheb24(k)*amom2
            amom0 = amom1
            amom1 = amom2
         ENDDO
         Result = res24
         Abserr = ABS(res24-res12)
      ELSE
!
!           apply the 15-point gauss-kronrod scheme.
!
         Krul = Krul - 1
         CALL QK15W(F,QWGTC,C,p2,p3,p4,kp,A,B,Result,Abserr,resabs,     &
                  & resasc)
         Neval = 15
         IF ( resasc==Abserr ) Krul = Krul + 1
      ENDIF
      END
!*==QC25F.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QC25F(F,A,B,Omega,Integr,Nrmom,Maxp1,Ksave,Result,     &
                     & Abserr,Neval,Resabs,Resasc,Momcom,Chebmo)
      IMPLICIT NONE
!*--QC25F14001
!***begin prologue  qc25f
!***date written   810101   (yymmdd)
!***revision date  211011   (yymmdd)
!***category no.  h2a2a2
!***keywords  integration rules for functions with cos or sin
!             factor, clenshaw-curtis, gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute the integral i=integral of f(x) over (a,b)
!            where w(x) = cos(omega*x) or (wx)=sin(omega*x)
!            and to compute j=integral of abs(f) over (a,b). for small
!            value of omega or small intervals (a,b) 15-point gauss-
!            kronrod rule used. otherwise generalized clenshaw-curtis us
!***description
!
!        integration rules for functions with cos or sin factor
!        standard fortran subroutine
!        real version
!
!        parameters
!         on entry
!           f      - real
!                    function subprogram defining the integrand
!                    function f(x). the actual name for f needs to
!                    be declared e x t e r n a l in the calling program.
!
!           a      - real
!                    lower limit of integration
!
!           b      - real
!                    upper limit of integration
!
!           omega  - real
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
!           result - real
!                    approximation to the integral i
!
!           abserr - real
!                    estimate of the modulus of the absolute
!                    error, which should equal or exceed abs(i-result)
!
!           neval  - integer
!                    number of integrand evaluations
!
!           resabs - real
!                    approximation to the integral j
!
!           resasc - real
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
!           chebmo - real
!                    array of dimension at least (maxp1,25) containing
!                    the modified chebyshev moments for the first momcom
!                    momcom interval lengths
!
!***references  (none)
!***routines called  qcheb,qk15w,qwgtf,r1mach,sgtsl
!***end prologue  qc25f
!
      REAL A , Abserr , ac , an , an2 , as , asap , ass , B , centr ,   &
         & Chebmo , cheb12 , cheb24 , conc , cons , cospar , d , QWGTF ,&
         & d1 , R1MACH , d2 , estc , ests , F , fval , hlgth , oflow ,  &
         & Omega , parint , par2 , par22 , p2 , p3 , p4 , Resabs ,      &
         & Resasc , resc12 , resc24 , ress12 , ress24 , Result ,        &
         & sinpar , v , x
      INTEGER i , iers , Integr , isym , j , k , Ksave , m , Maxp1 ,    &
            & Momcom , Neval , noequ , noeq1 , Nrmom
!
      DIMENSION Chebmo(Maxp1,25) , cheb12(13) , cheb24(25) , d(25) ,    &
              & d1(25) , d2(25) , fval(25) , v(28) , x(11)
!
      EXTERNAL F , QWGTF
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ...,11, to be used for the chebyshev expansion of f
!
      DATA x(1) , x(2) , x(3) , x(4) , x(5) , x(6) , x(7) , x(8) ,      &
         & x(9) , x(10) , x(11)/0.9914448613738104E+00 ,                &
         & 0.9659258262890683E+00 , 0.9238795325112868E+00 ,            &
         & 0.8660254037844386E+00 , 0.7933533402912352E+00 ,            &
         & 0.7071067811865475E+00 , 0.6087614290087206E+00 ,            &
         & 0.5000000000000000E+00 , 0.3826834323650898E+00 ,            &
         & 0.2588190451025208E+00 , 0.1305261922200516E+00/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fval   - value of the function f at the points
!                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5,
!                    k = 0, ..., 24
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
!***first executable statement  qc25f
      oflow = R1MACH(2)
!
      centr = 0.5E+00*(B+A)
      hlgth = 0.5E+00*(B-A)
      parint = Omega*hlgth
!
!           compute the integral using the 15-point gauss-kronrod
!           formula if the value of the parameter in the integrand
!           is small.
!
      IF ( ABS(parint)>0.2E+01 ) THEN
!
!           compute the integral using the generalized clenshaw-
!           curtis method.
!
         conc = hlgth*COS(centr*Omega)
         cons = hlgth*SIN(centr*Omega)
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
            par22 = par2 + 0.2E+01
            sinpar = SIN(parint)
            cospar = COS(parint)
!
!           compute the chebyshev moments with respect to cosine.
!
            v(1) = 0.2E+01*sinpar/parint
            v(2) = (0.8E+01*cospar+(par2+par2-0.8E+01)*sinpar/parint)   &
                 & /par2
            v(3) = (0.32E+02*(par2-0.12E+02)*cospar+(0.2E+01*((par2-    &
                 & 0.80E+02)*par2+0.192E+03)*sinpar)/parint)/(par2*par2)
            ac = 0.8E+01*cospar
            as = 0.24E+02*parint*sinpar
            IF ( ABS(parint)>0.24E+02 ) THEN
!
!           compute the chebyshev moments by means of forward
!           recursion.
!
               an = 0.4E+01
               DO i = 4 , 13
                  an2 = an*an
                  v(i) = ((an2-0.4E+01)*(0.2E+01*(par22-an2-an2)*v(i-1)-&
                       & ac)+as-par2*(an+0.1E+01)*(an+0.2E+01)*v(i-2))  &
                       & /(par2*(an-0.1E+01)*(an-0.2E+01))
                  an = an + 0.2E+01
               ENDDO
            ELSE
!
!           compute the chebyshev moments as the
!           solutions of a boundary value problem with 1
!           initial value (v(3)) and 1 end value (computed
!           using an asymptotic formula).
!
               noequ = 25
               noeq1 = noequ - 1
               an = 0.6E+01
               DO k = 1 , noeq1
                  an2 = an*an
                  d(k) = -0.2E+01*(an2-0.4E+01)*(par22-an2-an2)
                  d2(k) = (an-0.1E+01)*(an-0.2E+01)*par2
                  d1(k+1) = (an+0.3E+01)*(an+0.4E+01)*par2
                  v(k+3) = as - (an2-0.4E+01)*ac
                  an = an + 0.2E+01
               ENDDO
               an2 = an*an
               d(noequ) = -0.2E+01*(an2-0.4E+01)*(par22-an2-an2)
               v(noequ+3) = as - (an2-0.4E+01)*ac
               v(4) = v(4) - 0.56E+02*par2*v(3)
               ass = parint*sinpar
               asap = (((((0.210E+03*par2-0.1E+01)*cospar-(0.105E+03*   &
                    & par2-0.63E+02)*ass)/an2-(0.1E+01-0.15E+02*par2)   &
                    & *cospar+0.15E+02*ass)/an2-cospar+0.3E+01*ass)     &
                    & /an2-cospar)/an2
               v(noequ+3) = v(noequ+3) - 0.2E+01*asap*par2*(an-0.1E+01) &
                          & *(an-0.2E+01)
!
!           solve the tridiagonal system by means of gaussian
!           elimination with partial pivoting.
!
               CALL SGTSL(noequ,d1,d,d2,v(4),iers)
            ENDIF
            DO j = 1 , 13
               Chebmo(m,2*j-1) = v(j)
            ENDDO
!
!           compute the chebyshev moments with respect to sine.
!
            v(1) = 0.2E+01*(sinpar-parint*cospar)/par2
            v(2) = (0.18E+02-0.48E+02/par2)*sinpar/par2 +               &
                 & (-0.2E+01+0.48E+02/par2)*cospar/parint
            ac = -0.24E+02*parint*cospar
            as = -0.8E+01*sinpar
            IF ( ABS(parint)>0.24E+02 ) THEN
!
!           compute the chebyshev moments by means of
!           forward recursion.
!
               an = 0.3E+01
               DO i = 3 , 12
                  an2 = an*an
                  v(i) = ((an2-0.4E+01)*(0.2E+01*(par22-an2-an2)*v(i-1)+&
                       & as)+ac-par2*(an+0.1E+01)*(an+0.2E+01)*v(i-2))  &
                       & /(par2*(an-0.1E+01)*(an-0.2E+01))
                  an = an + 0.2E+01
               ENDDO
            ELSE
!
!           compute the chebyshev moments as the
!           solutions of a boundary value problem with 1
!           initial value (v(2)) and 1 end value (computed
!           using an asymptotic formula).
!
               an = 0.5E+01
               DO k = 1 , noeq1
                  an2 = an*an
                  d(k) = -0.2E+01*(an2-0.4E+01)*(par22-an2-an2)
                  d2(k) = (an-0.1E+01)*(an-0.2E+01)*par2
                  d1(k+1) = (an+0.3E+01)*(an+0.4E+01)*par2
                  v(k+2) = ac + (an2-0.4E+01)*as
                  an = an + 0.2E+01
               ENDDO
               an2 = an*an
               d(noequ) = -0.2E+01*(an2-0.4E+01)*(par22-an2-an2)
               v(noequ+2) = ac + (an2-0.4E+01)*as
               v(3) = v(3) - 0.42E+02*par2*v(2)
               ass = parint*cospar
               asap = (((((0.105E+03*par2-0.63E+02)*ass+(0.210E+03*par2-&
                    & 0.1E+01)*sinpar)/an2+(0.15E+02*par2-0.1E+01)      &
                    & *sinpar-0.15E+02*ass)/an2-0.3E+01*ass-sinpar)     &
                    & /an2-sinpar)/an2
               v(noequ+2) = v(noequ+2) - 0.2E+01*asap*par2*(an-0.1E+01) &
                          & *(an-0.2E+01)
!
!           solve the tridiagonal system by means of gaussian
!           elimination with partial pivoting.
!
               CALL SGTSL(noequ,d1,d,d2,v(3),iers)
            ENDIF
            DO j = 1 , 12
               Chebmo(m,2*j) = v(j)
            ENDDO
         ENDIF
         IF ( Nrmom<Momcom ) m = Nrmom + 1
         IF ( Momcom<Maxp1-1 .AND. Nrmom>=Momcom ) Momcom = Momcom + 1
!
!           compute the coefficients of the chebyshev expansions
!           of degrees 12 and 24 of the function f.
!
         fval(1) = 0.5E+00*F(centr+hlgth)
         fval(13) = F(centr)
         fval(25) = 0.5E+00*F(centr-hlgth)
         DO i = 2 , 12
            isym = 26 - i
            fval(i) = F(hlgth*x(i-1)+centr)
            fval(isym) = F(centr-hlgth*x(i-1))
         ENDDO
         CALL QCHEB(x,fval,cheb12,cheb24)
!
!           compute the integral and error estimates.
!
         resc12 = cheb12(13)*Chebmo(m,13)
         ress12 = 0.0E+00
         k = 11
         DO j = 1 , 6
            resc12 = resc12 + cheb12(k)*Chebmo(m,k)
            ress12 = ress12 + cheb12(k+1)*Chebmo(m,k+1)
            k = k - 2
         ENDDO
         resc24 = cheb24(25)*Chebmo(m,25)
         ress24 = 0.0E+00
         Resabs = ABS(cheb24(25))
         k = 23
         DO j = 1 , 12
            resc24 = resc24 + cheb24(k)*Chebmo(m,k)
            ress24 = ress24 + cheb24(k+1)*Chebmo(m,k+1)
            Resabs = Resabs + ABS(cheb24(k)) + ABS(cheb24(k+1))
            k = k - 2
         ENDDO
         estc = ABS(resc24-resc12)
         ests = ABS(ress24-ress12)
         Resabs = Resabs*ABS(hlgth)
         IF ( Integr==2 ) THEN
            Result = conc*ress24 + cons*resc24
            Abserr = ABS(conc*ests) + ABS(cons*estc)
         ELSE
            Result = conc*resc24 - cons*ress24
            Abserr = ABS(conc*estc) + ABS(cons*ests)
         ENDIF
      ELSE
         CALL QK15W(F,QWGTF,Omega,p2,p3,p4,Integr,A,B,Result,Abserr,    &
                  & Resabs,Resasc)
         Neval = 15
      ENDIF
      END
!*==QC25S.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QC25S(F,A,B,Bl,Br,Alfa,Beta,Ri,Rj,Rg,Rh,Result,Abserr, &
                     & Resasc,Integr,Nev)
      IMPLICIT NONE
!*--QC25S14357
!***begin prologue  qc25s
!***date written   810101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a2
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
!        real version
!
!        parameters
!           f      - real
!                    function subprogram defining the integrand
!                    f(x). the actual name for f needs to be declared
!                    e x t e r n a l  in the driver program.
!
!           a      - real
!                    left end point of the original interval
!
!           b      - real
!                    right end point of the original interval, b>a
!
!           bl     - real
!                    lower limit of integration, bl>=a
!
!           br     - real
!                    upper limit of integration, br<=b
!
!           alfa   - real
!                    parameter in the weight function
!
!           beta   - real
!                    parameter in the weight function
!
!           ri,rj,rg,rh - real
!                    modified chebyshev moments for the application
!                    of the generalized clenshaw-curtis
!                    method (computed in subroutine dqmomo)
!
!           result - real
!                    approximation to the integral
!                    result is computed by using a generalized
!                    clenshaw-curtis method if b1 = a or br = b.
!                    in all other cases the 15-point kronrod
!                    rule is applied, obtained by optimal addition of
!                    abscissae to the 7-point gauss rule.
!
!           abserr - real
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed abs(i-result)
!
!           resasc - real
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
!***references  (none)
!***routines called  qcheb,qk15w
!***end prologue  qc25s
!
      REAL A , Abserr , Alfa , B , Beta , Bl , Br , centr , cheb12 ,    &
         & cheb24 , dc , F , factor , fix , fval , hlgth , resabs ,     &
         & Resasc , Result , res12 , res24 , Rg , Rh , Ri , Rj , u ,    &
         & QWGTS , x
      INTEGER i , Integr , isym , Nev
!
      DIMENSION cheb12(13) , cheb24(25) , fval(25) , Rg(25) , Rh(25) ,  &
              & Ri(25) , Rj(25) , x(11)
!
      EXTERNAL F , QWGTS
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ..., 11, to be used for the computation of the
!           chebyshev series expansion of f.
!
      DATA x(1) , x(2) , x(3) , x(4) , x(5) , x(6) , x(7) , x(8) ,      &
         & x(9) , x(10) , x(11)/0.9914448613738104E+00 ,                &
         & 0.9659258262890683E+00 , 0.9238795325112868E+00 ,            &
         & 0.8660254037844386E+00 , 0.7933533402912352E+00 ,            &
         & 0.7071067811865475E+00 , 0.6087614290087206E+00 ,            &
         & 0.5000000000000000E+00 , 0.3826834323650898E+00 ,            &
         & 0.2588190451025208E+00 , 0.1305261922200516E+00/
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
!           qwgts - external function subprogram defining
!                    the four possible weight functions
!           hlgth  - half-length of the interval (bl,br)
!           centr  - mid point of the interval (bl,br)
!
!***first executable statement  qc25s
      Nev = 25
      IF ( Bl==A .AND. (Alfa/=0.0E+00 .OR. Integr==2 .OR. Integr==4) )  &
         & THEN
!
!           this part of the program is executed only if a = bl.
!           ----------------------------------------------------
!
!           compute the chebyshev series expansion of the
!           following function
!           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
!                  *f(0.5*(br-a)*x+0.5*(br+a))
!
         hlgth = 0.5E+00*(Br-Bl)
         centr = 0.5E+00*(Br+Bl)
         fix = B - centr
         fval(1) = 0.5E+00*F(hlgth+centr)*(fix-hlgth)**Beta
         fval(13) = F(centr)*(fix**Beta)
         fval(25) = 0.5E+00*F(centr-hlgth)*(fix+hlgth)**Beta
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)*(fix-u)**Beta
            fval(isym) = F(centr-u)*(fix+u)**Beta
         ENDDO
         factor = hlgth**(Alfa+0.1E+01)
         Result = 0.0E+00
         Abserr = 0.0E+00
         res12 = 0.0E+00
         res24 = 0.0E+00
         IF ( Integr>2 ) THEN
!
!           compute the chebyshev series expansion of the
!           following function
!           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
!
            fval(1) = fval(1)*ALOG(fix-hlgth)
            fval(13) = fval(13)*ALOG(fix)
            fval(25) = fval(25)*ALOG(fix+hlgth)
            DO i = 2 , 12
               u = hlgth*x(i-1)
               isym = 26 - i
               fval(i) = fval(i)*ALOG(fix-u)
               fval(isym) = fval(isym)*ALOG(fix+u)
            ENDDO
            CALL QCHEB(x,fval,cheb12,cheb24)
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
               dc = ALOG(Br-Bl)
               Result = res24*dc
               Abserr = ABS((res24-res12)*dc)
               res12 = 0.0E+00
               res24 = 0.0E+00
               DO i = 1 , 13
                  res12 = res12 + cheb12(i)*Rg(i)
                  res24 = res24 + cheb24(i)*Rg(i)
               ENDDO
               DO i = 14 , 25
                  res24 = res24 + cheb24(i)*Rg(i)
               ENDDO
            ENDIF
         ELSE
            CALL QCHEB(x,fval,cheb12,cheb24)
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
               dc = ALOG(Br-Bl)
               Result = res24*dc
               Abserr = ABS((res24-res12)*dc)
               res12 = 0.0E+00
               res24 = 0.0E+00
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
         Abserr = (Abserr+ABS(res24-res12))*factor
      ELSEIF ( Br==B .AND. (Beta/=0.0E+00 .OR. Integr==3 .OR. Integr==4)&
             & ) THEN
!
!           this part of the program is executed only if b = br.
!           ----------------------------------------------------
!
!           compute the chebyshev series expansion of the
!           following function
!           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
!                *f(0.5*(b-bl)*x+0.5*(b+bl))
!
         hlgth = 0.5E+00*(Br-Bl)
         centr = 0.5E+00*(Br+Bl)
         fix = centr - A
         fval(1) = 0.5E+00*F(hlgth+centr)*(fix+hlgth)**Alfa
         fval(13) = F(centr)*(fix**Alfa)
         fval(25) = 0.5E+00*F(centr-hlgth)*(fix-hlgth)**Alfa
         DO i = 2 , 12
            u = hlgth*x(i-1)
            isym = 26 - i
            fval(i) = F(u+centr)*(fix+u)**Alfa
            fval(isym) = F(centr-u)*(fix-u)**Alfa
         ENDDO
         factor = hlgth**(Beta+0.1E+01)
         Result = 0.0E+00
         Abserr = 0.0E+00
         res12 = 0.0E+00
         res24 = 0.0E+00
         IF ( Integr==2 .OR. Integr==4 ) THEN
!
!           compute the chebyshev series expansion of the
!           following function
!           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
!
            fval(1) = fval(1)*ALOG(hlgth+fix)
            fval(13) = fval(13)*ALOG(fix)
            fval(25) = fval(25)*ALOG(fix-hlgth)
            DO i = 2 , 12
               u = hlgth*x(i-1)
               isym = 26 - i
               fval(i) = fval(i)*ALOG(u+fix)
               fval(isym) = fval(isym)*ALOG(fix-u)
            ENDDO
            CALL QCHEB(x,fval,cheb12,cheb24)
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
               dc = ALOG(Br-Bl)
               Result = res24*dc
               Abserr = ABS((res24-res12)*dc)
               res12 = 0.0E+00
               res24 = 0.0E+00
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
            CALL QCHEB(x,fval,cheb12,cheb24)
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
               dc = ALOG(Br-Bl)
               Result = res24*dc
               Abserr = ABS((res24-res12)*dc)
               res12 = 0.0E+00
               res24 = 0.0E+00
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
         Abserr = (Abserr+ABS(res24-res12))*factor
      ELSE
!
!           if a>bl and b<br, apply the 15-point gauss-kronrod
!           scheme.
!
!
         CALL QK15W(F,QWGTS,A,B,Alfa,Beta,Integr,Bl,Br,Result,Abserr,   &
                  & resabs,Resasc)
         Nev = 15
      ENDIF
      END
!*==QCHEB.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QCHEB(X,Fval,Cheb12,Cheb24)
      IMPLICIT NONE
!*--QCHEB14699
!***begin prologue  qcheb
!***refer to  qc25c,qc25f,qc25s
!***routines called  (none)
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
!        real version
!
!        parameters
!          on entry
!           x      - real
!                    vector of dimension 11 containing the
!                    values cos(k*pi/24), k = 1, ..., 11
!
!           fval   - real
!                    vector of dimension 25 containing the
!                    function values at the points
!                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
!                    where (a,b) is the approximation interval.
!                    fval(1) and fval(25) are divided by two
!                    (these values are destroyed at output).
!
!          on return
!           cheb12 - real
!                    vector of dimension 13 containing the
!                    chebyshev coefficients for degree 12
!
!           cheb24 - real
!                    vector of dimension 25 containing the
!                    chebyshev coefficients for degree 24
!
!***end prologue  qcheb
!
      REAL alam , alam1 , alam2 , Cheb12 , Cheb24 , Fval , part1 ,      &
         & part2 , part3 , v , X
      INTEGER i , j
!
      DIMENSION Cheb12(13) , Cheb24(25) , Fval(25) , v(12) , X(11)
!
!***first executable statement  qcheb
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
           & + X(11)*v(12)
      Cheb24(2) = Cheb12(2) + alam
      Cheb24(24) = Cheb12(2) - alam
      alam = X(11)*v(2) - X(9)*v(4) + X(7)*v(6) - X(5)*v(8) + X(3)*v(10)&
           & - X(1)*v(12)
      Cheb24(12) = Cheb12(12) + alam
      Cheb24(14) = Cheb12(12) - alam
      alam1 = v(1) - part1 + part2
      alam2 = X(10)*v(3) - part3 + X(2)*v(11)
      Cheb12(6) = alam1 + alam2
      Cheb12(8) = alam1 - alam2
      alam = X(5)*v(2) - X(9)*v(4) - X(1)*v(6) - X(11)*v(8) + X(3)*v(10)&
           & + X(7)*v(12)
      Cheb24(6) = Cheb12(6) + alam
      Cheb24(20) = Cheb12(6) - alam
      alam = X(7)*v(2) - X(3)*v(4) - X(11)*v(6) + X(1)*v(8) - X(9)*v(10)&
           & - X(5)*v(12)
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
      alam = 0.1E+01/0.6E+01
      DO i = 2 , 12
         Cheb12(i) = Cheb12(i)*alam
      ENDDO
      alam = 0.5E+00*alam
      Cheb12(1) = Cheb12(1)*alam
      Cheb12(13) = Cheb12(13)*alam
      DO i = 2 , 24
         Cheb24(i) = Cheb24(i)*alam
      ENDDO
      Cheb24(1) = 0.5E+00*alam*Cheb24(1)
      Cheb24(25) = 0.5E+00*alam*Cheb24(25)
      END
!*==QELG.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QELG(N,Epstab,Result,Abserr,Res3la,Nres)
      IMPLICIT NONE
!*--QELG14849
!***begin prologue  qelg
!***refer to  qagie,qagoe,qagpe,qagse
!***routines called  r1mach
!***revision date  830518   (yymmdd)
!***keywords  epsilon algorithm, convergence acceleration,
!             extrapolation
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  the routine determines the limit of a given sequence of
!            approximations, by means of the epsilon algorithm of
!            p. wynn. an estimate of the absolute error is also given.
!            the condensed epsilon table is computed. only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***description
!
!           epsilon algorithm
!           standard fortran subroutine
!           real version
!
!           parameters
!              n      - integer
!                       epstab(n) contains the new element in the
!                       first column of the epsilon table.
!
!              epstab - real
!                       vector of dimension 52 containing the elements
!                       of the two lower diagonals of the triangular
!                       epsilon table. the elements are numbered
!                       starting at the right-hand corner of the
!                       triangle.
!
!              result - real
!                       resulting approximation to the integral
!
!              abserr - real
!                       estimate of the absolute error computed from
!                       result and the 3 previous results
!
!              res3la - real
!                       vector of dimension 3 containing the last 3
!                       results
!
!              nres   - integer
!                       number of calls to the routine
!                       (should be zero at first call)
!
!***end prologue  qelg
!
      REAL Abserr , delta1 , delta2 , delta3 , R1MACH , epmach ,        &
         & epsinf , Epstab , error , err1 , err2 , err3 , e0 , e1 ,     &
         & e1abs , e2 , e3 , oflow , res , Result , Res3la , ss , tol1 ,&
         & tol2 , tol3
      INTEGER i , ib , ib2 , ie , indx , k1 , k2 , k3 , limexp , N ,    &
            & newelm , Nres , num
      DIMENSION Epstab(52) , Res3la(3)
!
!           list of major variables
!           -----------------------
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
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
!***first executable statement  qelg
      epmach = R1MACH(4)
      oflow = R1MACH(2)
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
            e1abs = ABS(e1)
            delta2 = e2 - e1
            err2 = ABS(delta2)
            tol2 = AMAX1(ABS(e2),e1abs)*epmach
            delta3 = e1 - e0
            err3 = ABS(delta3)
            tol3 = AMAX1(e1abs,ABS(e0))*epmach
            IF ( err2>tol2 .OR. err3>tol3 ) THEN
               e3 = Epstab(k1)
               Epstab(k1) = e1
               delta1 = e1 - e3
               err1 = ABS(delta1)
               tol1 = AMAX1(e1abs,ABS(e3))*epmach
!
!           if two elements are very close to each other, omit
!           a part of the table by adjusting the value of n
!
               IF ( err1>tol1 .AND. err2>tol2 .AND. err3>tol3 ) THEN
                  ss = 0.1E+01/delta1 + 0.1E+01/delta2 - 0.1E+01/delta3
                  epsinf = ABS(ss*e1)
!
!           test to detect irregular behaviour in the table, and
!           eventually omit a part of the table adjusting the value
!           of n.
!
                  IF ( epsinf>0.1E-03 ) THEN
!
!           compute a new element and eventually adjust
!           the value of result.
!
                     res = e1 + 0.1E+01/ss
                     Epstab(k1) = res
                     k1 = k1 - 2
                     error = err2 + ABS(res-e2) + err3
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
            Abserr = ABS(Result-Res3la(3)) + ABS(Result-Res3la(2))      &
                   & + ABS(Result-Res3la(1))
            Res3la(1) = Res3la(2)
            Res3la(2) = Res3la(3)
            Res3la(3) = Result
         ELSE
            Res3la(Nres) = Result
            Abserr = oflow
         ENDIF
      ENDIF
 200  Abserr = AMAX1(Abserr,0.5E+01*epmach*ABS(Result))
      END
!*==QK15.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK15(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK1515046
!***begin prologue  qk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk15
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , &
         & resg , resk , reskh , Result , R1MACH , uflow , wg , wgk ,   &
         & xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8)/0.9914553711208126E+00 , 0.9491079123427585E+00 ,   &
         & 0.8648644233597691E+00 , 0.7415311855993944E+00 ,            &
         & 0.5860872354676911E+00 , 0.4058451513773972E+00 ,            &
         & 0.2077849550078985E+00 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8)/0.2293532201052922E-01 , 0.6309209262997855E-01 ,   &
         & 0.1047900103222502E+00 , 0.1406532597155259E+00 ,            &
         & 0.1690047266392679E+00 , 0.1903505780647854E+00 ,            &
         & 0.2044329400752989E+00 , 0.2094821410847278E+00/
      DATA wg(1) , wg(2) , wg(3) , wg(4)/0.1294849661688697E+00 ,       &
         & 0.2797053914892767E+00 , 0.3818300505051189E+00 ,            &
         & 0.4179591836734694E+00/
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
!***first executable statement  qk15
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(8)*ABS(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK15I.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK15I(F,Boun,Inf,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK15I15210
!***begin prologue  qk15i
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a2,h2a4a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       fuction subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              boun   - real
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
!              a      - real
!                       lower limit for integration over subrange
!                       of (0,1)
!
!              b      - real
!                       upper limit for integration over subrange
!                       of (0,1)
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule(resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule(resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integral of
!                       abs((transformed integrand)-i/(b-a)) over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk15i
!
      REAL A , absc , absc1 , absc2 , Abserr , B , Boun , centr , dinf ,&
         & R1MACH , epmach , F , fc , fsum , fval1 , fval2 , fv1 , fv2 ,&
         & hlgth , Resabs , Resasc , resg , resk , reskh , Result ,     &
         & tabsc1 , tabsc2 , uflow , wg , wgk , xgk
      INTEGER Inf , j , MIN0
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8)/0.9914553711208126E+00 , 0.9491079123427585E+00 ,   &
         & 0.8648644233597691E+00 , 0.7415311855993944E+00 ,            &
         & 0.5860872354676911E+00 , 0.4058451513773972E+00 ,            &
         & 0.2077849550078985E+00 , 0.0000000000000000E+00/
!
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8)/0.2293532201052922E-01 , 0.6309209262997855E-01 ,   &
         & 0.1047900103222502E+00 , 0.1406532597155259E+00 ,            &
         & 0.1690047266392679E+00 , 0.1903505780647854E+00 ,            &
         & 0.2044329400752989E+00 , 0.2094821410847278E+00/
!
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) ,      &
         & wg(8)/0.0000000000000000E+00 , 0.1294849661688697E+00 ,      &
         & 0.0000000000000000E+00 , 0.2797053914892767E+00 ,            &
         & 0.0000000000000000E+00 , 0.3818300505051189E+00 ,            &
         & 0.0000000000000000E+00 , 0.4179591836734694E+00/
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
!***first executable statement  qk15i
      epmach = R1MACH(4)
      uflow = R1MACH(1)
      dinf = MIN0(1,Inf)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      tabsc1 = Boun + dinf*(0.1E+01-centr)/centr
      fval1 = F(tabsc1)
      IF ( Inf==2 ) fval1 = fval1 + F(-tabsc1)
      fc = (fval1/centr)/centr
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the error.
!
      resg = wg(8)*fc
      resk = wgk(8)*fc
      Resabs = ABS(resk)
      DO j = 1 , 7
         absc = hlgth*xgk(j)
         absc1 = centr - absc
         absc2 = centr + absc
         tabsc1 = Boun + dinf*(0.1E+01-absc1)/absc1
         tabsc2 = Boun + dinf*(0.1E+01-absc2)/absc2
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
         Resabs = Resabs + wgk(j)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(8)*ABS(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resasc = Resasc*hlgth
      Resabs = Resabs*hlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.E0 )                         &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK15W.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK15W(F,W,P1,P2,P3,P4,Kp,A,B,Result,Abserr,Resabs,     &
                     & Resasc)
      IMPLICIT NONE
!*--QK15W15400
!***begin prologue  qk15w
!***date written   810101   (yymmdd)
!***revision date  830518   (mmddyy)
!***category no.  h2a2a2
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
!           real version
!
!           parameters
!             on entry
!              f      - real
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the driver program.
!
!              w      - real
!                       function subprogram defining the integrand
!                       weight function w(x). the actual name for w
!                       needs to be declared e x t e r n a l in the
!                       calling program.
!
!              p1, p2, p3, p4 - real
!                       parameters in the weight function
!
!              kp     - integer
!                       key for indicating the type of weight function
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule (resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral of abs(f)
!
!              resasc - real
!                       approximation to the integral of abs(f-i/(b-a))
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk15w
!
      REAL A , absc , absc1 , absc2 , Abserr , B , centr , dhlgth ,     &
         & R1MACH , epmach , F , fc , fsum , fval1 , fval2 , fv1 , fv2 ,&
         & hlgth , P1 , P2 , P3 , P4 , Resabs , Resasc , resg , resk ,  &
         & reskh , Result , uflow , W , wg , wgk , xgk
      INTEGER j , jtw , jtwm1 , Kp
      EXTERNAL F , W
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
         & , xgk(8)/0.9914553711208126E+00 , 0.9491079123427585E+00 ,   &
         & 0.8648644233597691E+00 , 0.7415311855993944E+00 ,            &
         & 0.5860872354676911E+00 , 0.4058451513773972E+00 ,            &
         & 0.2077849550078985E+00 , 0.0000000000000000E+00/
!
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8)/0.2293532201052922E-01 , 0.6309209262997855E-01 ,   &
         & 0.1047900103222502E+00 , 0.1406532597155259E+00 ,            &
         & 0.1690047266392679E+00 , 0.1903505780647854E+00 ,            &
         & 0.2044329400752989E+00 , 0.2094821410847278E+00/
!
      DATA wg(1) , wg(2) , wg(3) , wg(4)/0.1294849661688697E+00 ,       &
         & 0.2797053914892767E+00 , 0.3818300505051189E+00 ,            &
         & 0.4179591836734694E+00/
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
!***first executable statement  qk15w
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 15-point kronrod approximation to the
!           integral, and estimate the error.
!
      fc = F(centr)*W(centr,P1,P2,P3,P4,Kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(8)*ABS(fc-reskh)
      DO j = 1 , 7
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK21.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK21(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK2115581
!*** Start of declarations inserted by SPAG
      REAL Resasc
!*** End of declarations inserted by SPAG
!***begin prologue  qk21
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the driver program.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 21-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 10-point gauss rule (resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk21
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , resg ,   &
         & resk , reskh , Result , R1MACH , uflow , wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8) , xgk(9) , xgk(10) , xgk(11)                        &
         & /0.9956571630258081E+00 , 0.9739065285171717E+00 ,           &
         & 0.9301574913557082E+00 , 0.8650633666889845E+00 ,            &
         & 0.7808177265864169E+00 , 0.6794095682990244E+00 ,            &
         & 0.5627571346686047E+00 , 0.4333953941292472E+00 ,            &
         & 0.2943928627014602E+00 , 0.1488743389816312E+00 ,            &
         & 0.0000000000000000E+00/
!
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8) , wgk(9) , wgk(10) , wgk(11)                        &
         & /0.1169463886737187E-01 , 0.3255816230796473E-01 ,           &
         & 0.5475589657435200E-01 , 0.7503967481091995E-01 ,            &
         & 0.9312545458369761E-01 , 0.1093871588022976E+00 ,            &
         & 0.1234919762620659E+00 , 0.1347092173114733E+00 ,            &
         & 0.1427759385770601E+00 , 0.1477391049013385E+00 ,            &
         & 0.1494455540029169E+00/
!
      DATA wg(1) , wg(2) , wg(3) , wg(4) ,                              &
         & wg(5)/0.6667134430868814E-01 , 0.1494513491505806E+00 ,      &
         & 0.2190863625159820E+00 , 0.2692667193099964E+00 ,            &
         & 0.2955242247147529E+00/
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
!***first executable statement  qk21
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 21-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      resg = 0.0E+00
      fc = F(centr)
      resk = wgk(11)*fc
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(11)*ABS(fc-reskh)
      DO j = 1 , 10
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK31.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK31(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK3115757
!***begin prologue  qk31
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 31-point
!                       gauss-kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point gauss
!                       rule (resg).
!
!              abserr - real
!                       estimate of the modulus of the modulus,
!                       which should not exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk31
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , &
         & resg , resk , reskh , Result , R1MACH , uflow , wg , wgk ,   &
         & xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8) , xgk(9) , xgk(10) , xgk(11) , xgk(12) , xgk(13) ,  &
         & xgk(14) , xgk(15) , xgk(16)/0.9980022986933971E+00 ,         &
         & 0.9879925180204854E+00 , 0.9677390756791391E+00 ,            &
         & 0.9372733924007059E+00 , 0.8972645323440819E+00 ,            &
         & 0.8482065834104272E+00 , 0.7904185014424659E+00 ,            &
         & 0.7244177313601700E+00 , 0.6509967412974170E+00 ,            &
         & 0.5709721726085388E+00 , 0.4850818636402397E+00 ,            &
         & 0.3941513470775634E+00 , 0.2991800071531688E+00 ,            &
         & 0.2011940939974345E+00 , 0.1011420669187175E+00 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8) , wgk(9) , wgk(10) , wgk(11) , wgk(12) , wgk(13) ,  &
         & wgk(14) , wgk(15) , wgk(16)/0.5377479872923349E-02 ,         &
         & 0.1500794732931612E-01 , 0.2546084732671532E-01 ,            &
         & 0.3534636079137585E-01 , 0.4458975132476488E-01 ,            &
         & 0.5348152469092809E-01 , 0.6200956780067064E-01 ,            &
         & 0.6985412131872826E-01 , 0.7684968075772038E-01 ,            &
         & 0.8308050282313302E-01 , 0.8856444305621177E-01 ,            &
         & 0.9312659817082532E-01 , 0.9664272698362368E-01 ,            &
         & 0.9917359872179196E-01 , 0.1007698455238756E+00 ,            &
         & 0.1013300070147915E+00/
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) ,      &
         & wg(8)/0.3075324199611727E-01 , 0.7036604748810812E-01 ,      &
         & 0.1071592204671719E+00 , 0.1395706779261543E+00 ,            &
         & 0.1662692058169939E+00 , 0.1861610000155622E+00 ,            &
         & 0.1984314853271116E+00 , 0.2025782419255613E+00/
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
!
!***first executable statement  qk31
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 31-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(16)*ABS(fc-reskh)
      DO j = 1 , 15
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK41.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK41(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK4115932
!***begin prologue  qk41
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 41-point
!                       gauss-kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point gauss
!                       rule (resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integal of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk41
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , &
         & resg , resk , reskh , Result , R1MACH , uflow , wg , wgk ,   &
         & xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8) , xgk(9) , xgk(10) , xgk(11) , xgk(12) , xgk(13) ,  &
         & xgk(14) , xgk(15) , xgk(16) , xgk(17) , xgk(18) , xgk(19) ,  &
         & xgk(20) , xgk(21)/0.9988590315882777E+00 ,                   &
         & 0.9931285991850949E+00 , 0.9815078774502503E+00 ,            &
         & 0.9639719272779138E+00 , 0.9408226338317548E+00 ,            &
         & 0.9122344282513259E+00 , 0.8782768112522820E+00 ,            &
         & 0.8391169718222188E+00 , 0.7950414288375512E+00 ,            &
         & 0.7463319064601508E+00 , 0.6932376563347514E+00 ,            &
         & 0.6360536807265150E+00 , 0.5751404468197103E+00 ,            &
         & 0.5108670019508271E+00 , 0.4435931752387251E+00 ,            &
         & 0.3737060887154196E+00 , 0.3016278681149130E+00 ,            &
         & 0.2277858511416451E+00 , 0.1526054652409227E+00 ,            &
         & 0.7652652113349733E-01 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8) , wgk(9) , wgk(10) , wgk(11) , wgk(12) , wgk(13) ,  &
         & wgk(14) , wgk(15) , wgk(16) , wgk(17) , wgk(18) , wgk(19) ,  &
         & wgk(20) , wgk(21)/0.3073583718520532E-02 ,                   &
         & 0.8600269855642942E-02 , 0.1462616925697125E-01 ,            &
         & 0.2038837346126652E-01 , 0.2588213360495116E-01 ,            &
         & 0.3128730677703280E-01 , 0.3660016975820080E-01 ,            &
         & 0.4166887332797369E-01 , 0.4643482186749767E-01 ,            &
         & 0.5094457392372869E-01 , 0.5519510534828599E-01 ,            &
         & 0.5911140088063957E-01 , 0.6265323755478117E-01 ,            &
         & 0.6583459713361842E-01 , 0.6864867292852162E-01 ,            &
         & 0.7105442355344407E-01 , 0.7303069033278667E-01 ,            &
         & 0.7458287540049919E-01 , 0.7570449768455667E-01 ,            &
         & 0.7637786767208074E-01 , 0.7660071191799966E-01/
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) ,      &
         & wg(8) , wg(9) , wg(10)/0.1761400713915212E-01 ,              &
         & 0.4060142980038694E-01 , 0.6267204833410906E-01 ,            &
         & 0.8327674157670475E-01 , 0.1019301198172404E+00 ,            &
         & 0.1181945319615184E+00 , 0.1316886384491766E+00 ,            &
         & 0.1420961093183821E+00 , 0.1491729864726037E+00 ,            &
         & 0.1527533871307259E+00/
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
!***first executable statement  qk41
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 41-point gauss-kronrod approximation to
!           the integral, and estimate the absolute error.
!
      resg = 0.0E+00
      fc = F(centr)
      resk = wgk(21)*fc
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(21)*ABS(fc-reskh)
      DO j = 1 , 20
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.E+00 )                       &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK51.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK51(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK5116119
!***begin prologue  qk51
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
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
!           real version
!
!           parameters
!            on entry
!              f      - real
!                       function subroutine defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!            on return
!              result - real
!                       approximation to the integral i
!                       result is computed by applying the 51-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point gauss rule (resg).
!
!              abserr - real
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real
!                       approximation to the integral j
!
!              resasc - real
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk51
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , &
         & resg , resk , reskh , Result , R1MACH , uflow , wg , wgk ,   &
         & xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8) , xgk(9) , xgk(10) , xgk(11) , xgk(12) , xgk(13) ,  &
         & xgk(14)/0.9992621049926098E+00 , 0.9955569697904981E+00 ,    &
         & 0.9880357945340772E+00 , 0.9766639214595175E+00 ,            &
         & 0.9616149864258425E+00 , 0.9429745712289743E+00 ,            &
         & 0.9207471152817016E+00 , 0.8949919978782754E+00 ,            &
         & 0.8658470652932756E+00 , 0.8334426287608340E+00 ,            &
         & 0.7978737979985001E+00 , 0.7592592630373576E+00 ,            &
         & 0.7177664068130844E+00 , 0.6735663684734684E+00/
      DATA xgk(15) , xgk(16) , xgk(17) , xgk(18) , xgk(19) , xgk(20) ,  &
         & xgk(21) , xgk(22) , xgk(23) , xgk(24) , xgk(25) , xgk(26)    &
         & /0.6268100990103174E+00 , 0.5776629302412230E+00 ,           &
         & 0.5263252843347192E+00 , 0.4730027314457150E+00 ,            &
         & 0.4178853821930377E+00 , 0.3611723058093878E+00 ,            &
         & 0.3030895389311078E+00 , 0.2438668837209884E+00 ,            &
         & 0.1837189394210489E+00 , 0.1228646926107104E+00 ,            &
         & 0.6154448300568508E-01 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8) , wgk(9) , wgk(10) , wgk(11) , wgk(12) , wgk(13) ,  &
         & wgk(14)/0.1987383892330316E-02 , 0.5561932135356714E-02 ,    &
         & 0.9473973386174152E-02 , 0.1323622919557167E-01 ,            &
         & 0.1684781770912830E-01 , 0.2043537114588284E-01 ,            &
         & 0.2400994560695322E-01 , 0.2747531758785174E-01 ,            &
         & 0.3079230016738749E-01 , 0.3400213027432934E-01 ,            &
         & 0.3711627148341554E-01 , 0.4008382550403238E-01 ,            &
         & 0.4287284502017005E-01 , 0.4550291304992179E-01/
      DATA wgk(15) , wgk(16) , wgk(17) , wgk(18) , wgk(19) , wgk(20) ,  &
         & wgk(21) , wgk(22) , wgk(23) , wgk(24) , wgk(25) , wgk(26)    &
         & /0.4798253713883671E-01 , 0.5027767908071567E-01 ,           &
         & 0.5236288580640748E-01 , 0.5425112988854549E-01 ,            &
         & 0.5595081122041232E-01 , 0.5743711636156783E-01 ,            &
         & 0.5868968002239421E-01 , 0.5972034032417406E-01 ,            &
         & 0.6053945537604586E-01 , 0.6112850971705305E-01 ,            &
         & 0.6147118987142532E-01 , 0.6158081806783294E-01/
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) ,      &
         & wg(8) , wg(9) , wg(10) , wg(11) , wg(12) , wg(13)            &
         & /0.1139379850102629E-01 , 0.2635498661503214E-01 ,           &
         & 0.4093915670130631E-01 , 0.5490469597583519E-01 ,            &
         & 0.6803833381235692E-01 , 0.8014070033500102E-01 ,            &
         & 0.9102826198296365E-01 , 0.1005359490670506E+00 ,            &
         & 0.1085196244742637E+00 , 0.1148582591457116E+00 ,            &
         & 0.1194557635357848E+00 , 0.1222424429903100E+00 ,            &
         & 0.1231760537267155E+00/
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
!***first executable statement  qk51
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(A+B)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 51-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
      fc = F(centr)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      Resabs = ABS(resk)
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
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
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
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(26)*ABS(fc-reskh)
      DO j = 1 , 25
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QK61.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QK61(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK6116313
!***begin prologue  qk61
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  61-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of dabs(f) over (a,b)
!***description
!
!        integration rule
!        standard fortran subroutine
!        real version
!
!
!        parameters
!         on entry
!           f      - real
!                    function subprogram defining the integrand
!                    function f(x). the actual name for f needs to be
!                    declared e x t e r n a l in the calling program.
!
!           a      - real
!                    lower limit of integration
!
!           b      - real
!                    upper limit of integration
!
!         on return
!           result - real
!                    approximation to the integral i
!                    result is computed by applying the 61-point
!                    kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point gauss rule (resg).
!
!           abserr - real
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed dabs(i-result)
!
!           resabs - real
!                    approximation to the integral j
!
!           resasc - real
!                    approximation to the integral of dabs(f-i/(b-a))
!
!
!***references  (none)
!***routines called  r1mach
!***end prologue  qk61
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc ,   &
         & fsum , fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , &
         & resg , resk , reskh , Result , R1MACH , uflow , wg , wgk ,   &
         & xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
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
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) &
         & , xgk(8) , xgk(9) , xgk(10)/0.9994844100504906E+00 ,         &
         & 0.9968934840746495E+00 , 0.9916309968704046E+00 ,            &
         & 0.9836681232797472E+00 , 0.9731163225011263E+00 ,            &
         & 0.9600218649683075E+00 , 0.9443744447485600E+00 ,            &
         & 0.9262000474292743E+00 , 0.9055733076999078E+00 ,            &
         & 0.8825605357920527E+00/
      DATA xgk(11) , xgk(12) , xgk(13) , xgk(14) , xgk(15) , xgk(16) ,  &
         & xgk(17) , xgk(18) , xgk(19) , xgk(20)                        &
         & /0.8572052335460611E+00 , 0.8295657623827684E+00 ,           &
         & 0.7997278358218391E+00 , 0.7677774321048262E+00 ,            &
         & 0.7337900624532268E+00 , 0.6978504947933158E+00 ,            &
         & 0.6600610641266270E+00 , 0.6205261829892429E+00 ,            &
         & 0.5793452358263617E+00 , 0.5366241481420199E+00/
      DATA xgk(21) , xgk(22) , xgk(23) , xgk(24) , xgk(25) , xgk(26) ,  &
         & xgk(27) , xgk(28) , xgk(29) , xgk(30) , xgk(31)              &
         & /0.4924804678617786E+00 , 0.4470337695380892E+00 ,           &
         & 0.4004012548303944E+00 , 0.3527047255308781E+00 ,            &
         & 0.3040732022736251E+00 , 0.2546369261678898E+00 ,            &
         & 0.2045251166823099E+00 , 0.1538699136085835E+00 ,            &
         & 0.1028069379667370E+00 , 0.5147184255531770E-01 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) &
         & , wgk(8) , wgk(9) , wgk(10)/0.1389013698677008E-02 ,         &
         & 0.3890461127099884E-02 , 0.6630703915931292E-02 ,            &
         & 0.9273279659517763E-02 , 0.1182301525349634E-01 ,            &
         & 0.1436972950704580E-01 , 0.1692088918905327E-01 ,            &
         & 0.1941414119394238E-01 , 0.2182803582160919E-01 ,            &
         & 0.2419116207808060E-01/
      DATA wgk(11) , wgk(12) , wgk(13) , wgk(14) , wgk(15) , wgk(16) ,  &
         & wgk(17) , wgk(18) , wgk(19) , wgk(20)                        &
         & /0.2650995488233310E-01 , 0.2875404876504129E-01 ,           &
         & 0.3090725756238776E-01 , 0.3298144705748373E-01 ,            &
         & 0.3497933802806002E-01 , 0.3688236465182123E-01 ,            &
         & 0.3867894562472759E-01 , 0.4037453895153596E-01 ,            &
         & 0.4196981021516425E-01 , 0.4345253970135607E-01/
      DATA wgk(21) , wgk(22) , wgk(23) , wgk(24) , wgk(25) , wgk(26) ,  &
         & wgk(27) , wgk(28) , wgk(29) , wgk(30) , wgk(31)              &
         & /0.4481480013316266E-01 , 0.4605923827100699E-01 ,           &
         & 0.4718554656929915E-01 , 0.4818586175708713E-01 ,            &
         & 0.4905543455502978E-01 , 0.4979568342707421E-01 ,            &
         & 0.5040592140278235E-01 , 0.5088179589874961E-01 ,            &
         & 0.5122154784925877E-01 , 0.5142612853745903E-01 ,            &
         & 0.5149472942945157E-01/
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) ,      &
         & wg(8)/0.7968192496166606E-02 , 0.1846646831109096E-01 ,      &
         & 0.2878470788332337E-01 , 0.3879919256962705E-01 ,            &
         & 0.4840267283059405E-01 , 0.5749315621761907E-01 ,            &
         & 0.6597422988218050E-01 , 0.7375597473770521E-01/
      DATA wg(9) , wg(10) , wg(11) , wg(12) , wg(13) , wg(14) , wg(15)  &
         & /0.8075589522942022E-01 , 0.8689978720108298E-01 ,           &
         & 0.9212252223778613E-01 , 0.9636873717464426E-01 ,            &
         & 0.9959342058679527E-01 , 0.1017623897484055E+00 ,            &
         & 0.1028526528935588E+00/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
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
!***first executable statement  qk61
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(B+A)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           compute the 61-point kronrod approximation to the
!           integral, and estimate the absolute error.
!
      resg = 0.0E+00
      fc = F(centr)
      resk = wgk(31)*fc
      Resabs = ABS(resk)
      DO j = 1 , 15
         jtw = j*2
         absc = hlgth*xgk(jtw)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtw) = fval1
         fv2(jtw) = fval2
         fsum = fval1 + fval2
         resg = resg + wg(j)*fsum
         resk = resk + wgk(jtw)*fsum
         Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
      ENDDO
      DO j = 1 , 15
         jtwm1 = j*2 - 1
         absc = hlgth*xgk(jtwm1)
         fval1 = F(centr-absc)
         fval2 = F(centr+absc)
         fv1(jtwm1) = fval1
         fv2(jtwm1) = fval2
         fsum = fval1 + fval2
         resk = resk + wgk(jtwm1)*fsum
         Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(31)*ABS(fc-reskh)
      DO j = 1 , 30
         Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                      &
         & Abserr = Resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/Resasc)        &
         & **1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) )                              &
         & Abserr = AMAX1((epmach*0.5E+02)*Resabs,Abserr)
      END
!*==QMOMO.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QMOMO(Alfa,Beta,Ri,Rj,Rg,Rh,Integr)
      IMPLICIT NONE
!*--QMOMO16517
!***begin prologue  qmomo
!***date written   810101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a2a1,c3a2
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
!        real version
!
!        parameters
!           alfa   - real
!                    parameter in the weight function w(x), alfa>(-1)
!
!           beta   - real
!                    parameter in the weight function w(x), beta>(-1)
!
!           ri     - real
!                    vector of dimension 25
!                    ri(k) is the integral over (-1,1) of
!                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
!
!           rj     - real
!                    vector of dimension 25
!                    rj(k) is the integral over (-1,1) of
!                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
!
!           rg     - real
!                    vector of dimension 25
!                    rg(k) is the integral over (-1,1) of
!                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
!
!           rh     - real
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
!***references  (none)
!***routines called  (none)
!***end prologue  qmomo
!
      REAL Alfa , alfp1 , alfp2 , an , anm1 , Beta , betp1 , betp2 ,    &
         & ralf , rbet , Rg , Rh , Ri , Rj
      INTEGER i , im1 , Integr
!
      DIMENSION Rg(25) , Rh(25) , Ri(25) , Rj(25)
!
!
!***first executable statement  qmomo
      alfp1 = Alfa + 0.1E+01
      betp1 = Beta + 0.1E+01
      alfp2 = Alfa + 0.2E+01
      betp2 = Beta + 0.2E+01
      ralf = 0.2E+01**alfp1
      rbet = 0.2E+01**betp1
!
!           compute ri, rj using a forward recurrence relation.
!
      Ri(1) = ralf/alfp1
      Rj(1) = rbet/betp1
      Ri(2) = Ri(1)*Alfa/alfp2
      Rj(2) = Rj(1)*Beta/betp2
      an = 0.2E+01
      anm1 = 0.1E+01
      DO i = 3 , 25
         Ri(i) = -(ralf+an*(an-alfp2)*Ri(i-1))/(anm1*(an+alfp1))
         Rj(i) = -(rbet+an*(an-betp2)*Rj(i-1))/(anm1*(an+betp1))
         anm1 = an
         an = an + 0.1E+01
      ENDDO
      IF ( Integr/=1 ) THEN
         IF ( Integr/=3 ) THEN
!
!           compute rg using a forward recurrence relation.
!
            Rg(1) = -Ri(1)/alfp1
            Rg(2) = -(ralf+ralf)/(alfp2*alfp2) - Rg(1)
            an = 0.2E+01
            anm1 = 0.1E+01
            im1 = 2
            DO i = 3 , 25
               Rg(i) = -(an*(an-alfp2)*Rg(im1)-an*Ri(im1)+anm1*Ri(i))   &
                     & /(anm1*(an+alfp1))
               anm1 = an
               an = an + 0.1E+01
               im1 = i
            ENDDO
            IF ( Integr==2 ) GOTO 100
         ENDIF
!
!           compute rh using a forward recurrence relation.
!
         Rh(1) = -Rj(1)/betp1
         Rh(2) = -(rbet+rbet)/(betp2*betp2) - Rh(1)
         an = 0.2E+01
         anm1 = 0.1E+01
         im1 = 2
         DO i = 3 , 25
            Rh(i) = -(an*(an-betp2)*Rh(im1)-an*Rj(im1)+anm1*Rj(i))      &
                  & /(anm1*(an+betp1))
            anm1 = an
            an = an + 0.1E+01
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
!*==QNG.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QNG(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier)
      IMPLICIT NONE
!*--QNG16647
!*** Start of declarations inserted by SPAG
      REAL w21b
!*** End of declarations inserted by SPAG
!***begin prologue  qng
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
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
! real version
!
!           f      - real version
!                    function subprogram defining the integrand function
!                    f(x). the actual name for f needs to be declared
!                    e x t e r n a l in the driver program.
!
!           a      - real version
!                    lower limit of integration
!
!           b      - real version
!                    upper limit of integration
!
!           epsabs - real
!                    absolute accuracy requested
!           epsrel - real
!                    relative accuracy requested
!                    if  epsabs<=0
!                    and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                    the routine will end with ier = 6.
!
!         on return
!           result - real
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
!           abserr - real
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
!                            epsrel<max(50*rel.mach.acc.,0.5d-28).
!                            result, abserr and neval are set to zero.
!
!***references  (none)
!***routines called  r1mach,xerror
!***end prologue  qng
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , Epsabs ,   &
         & Epsrel , F , fcentr , fval , fval1 , fval2 , fv1 , fv2 ,     &
         & fv3 , fv4 , hlgth , Result , res10 , res21 , res43 , res87 , &
         & resabs , resasc , reskh , R1MACH , savfun , uflow , w10 ,    &
         & w21a , w43a , w43b , w87a , w87b , x1 , x2 , x3 , x4
      INTEGER Ier , ipx , k , l , Neval
      EXTERNAL F
!
      DIMENSION fv1(5) , fv2(5) , fv3(5) , fv4(5) , x1(5) , x2(5) ,     &
              & x3(11) , x4(22) , w10(5) , w21a(5) , w21b(6) , w43a(10) &
              & , w43b(12) , w87a(21) , w87b(23) , savfun(21)
!
!           the following data statements contain the
!           abscissae and weights of the integration rules used.
!
!           x1      abscissae common to the 10-, 21-, 43-
!                   and 87-point rule
!           x2      abscissae common to the 21-, 43- and
!                   87-point rule
!           x3      abscissae common to the 43- and 87-point
!                   rule
!           x4      abscissae of the 87-point rule
!           w10     weights of the 10-point formula
!           w21a    weights of the 21-point formula for
!                   abscissae x1
!           w21b    weights of the 21-point formula for
!                   abscissae x2
!           w43a    weights of the 43-point formula for
!                   abscissae x1, x3
!           w43b    weights of the 43-point formula for
!                   abscissae x3
!           w87a    weights of the 87-point formula for
!                   abscissae x1, x2, x3
!           w87b    weights of the 87-point formula for
!                   abscissae x4
!
      DATA x1(1) , x1(2) , x1(3) , x1(4) ,                              &
         & x1(5)/0.9739065285171717E+00 , 0.8650633666889845E+00 ,      &
         & 0.6794095682990244E+00 , 0.4333953941292472E+00 ,            &
         & 0.1488743389816312E+00/
      DATA x2(1) , x2(2) , x2(3) , x2(4) ,                              &
         & x2(5)/0.9956571630258081E+00 , 0.9301574913557082E+00 ,      &
         & 0.7808177265864169E+00 , 0.5627571346686047E+00 ,            &
         & 0.2943928627014602E+00/
      DATA x3(1) , x3(2) , x3(3) , x3(4) , x3(5) , x3(6) , x3(7) ,      &
         & x3(8) , x3(9) , x3(10) , x3(11)/0.9993333609019321E+00 ,     &
         & 0.9874334029080889E+00 , 0.9548079348142663E+00 ,            &
         & 0.9001486957483283E+00 , 0.8251983149831142E+00 ,            &
         & 0.7321483889893050E+00 , 0.6228479705377252E+00 ,            &
         & 0.4994795740710565E+00 , 0.3649016613465808E+00 ,            &
         & 0.2222549197766013E+00 , 0.7465061746138332E-01/
      DATA x4(1) , x4(2) , x4(3) , x4(4) , x4(5) , x4(6) , x4(7) ,      &
         & x4(8) , x4(9) , x4(10) , x4(11) , x4(12) , x4(13) , x4(14) , &
         & x4(15) , x4(16) , x4(17) , x4(18) , x4(19) , x4(20) , x4(21) &
         & , x4(22)/0.9999029772627292E+00 , 0.9979898959866787E+00 ,   &
         & 0.9921754978606872E+00 , 0.9813581635727128E+00 ,            &
         & 0.9650576238583846E+00 , 0.9431676131336706E+00 ,            &
         & 0.9158064146855072E+00 , 0.8832216577713165E+00 ,            &
         & 0.8457107484624157E+00 , 0.8035576580352310E+00 ,            &
         & 0.7570057306854956E+00 , 0.7062732097873218E+00 ,            &
         & 0.6515894665011779E+00 , 0.5932233740579611E+00 ,            &
         & 0.5314936059708319E+00 , 0.4667636230420228E+00 ,            &
         & 0.3994248478592188E+00 , 0.3298748771061883E+00 ,            &
         & 0.2585035592021616E+00 , 0.1856953965683467E+00 ,            &
         & 0.1118422131799075E+00 , 0.3735212339461987E-01/
      DATA w10(1) , w10(2) , w10(3) , w10(4) , w10(5)                   &
         & /0.6667134430868814E-01 , 0.1494513491505806E+00 ,           &
         & 0.2190863625159820E+00 , 0.2692667193099964E+00 ,            &
         & 0.2955242247147529E+00/
      DATA w21a(1) , w21a(2) , w21a(3) , w21a(4) , w21a(5)              &
         & /0.3255816230796473E-01 , 0.7503967481091995E-01 ,           &
         & 0.1093871588022976E+00 , 0.1347092173114733E+00 ,            &
         & 0.1477391049013385E+00/
      DATA w21b(1) , w21b(2) , w21b(3) , w21b(4) , w21b(5) , w21b(6)    &
         & /0.1169463886737187E-01 , 0.5475589657435200E-01 ,           &
         & 0.9312545458369761E-01 , 0.1234919762620659E+00 ,            &
         & 0.1427759385770601E+00 , 0.1494455540029169E+00/
      DATA w43a(1) , w43a(2) , w43a(3) , w43a(4) , w43a(5) , w43a(6) ,  &
         & w43a(7) , w43a(8) , w43a(9) , w43a(10)                       &
         & /0.1629673428966656E-01 , 0.3752287612086950E-01 ,           &
         & 0.5469490205825544E-01 , 0.6735541460947809E-01 ,            &
         & 0.7387019963239395E-01 , 0.5768556059769796E-02 ,            &
         & 0.2737189059324884E-01 , 0.4656082691042883E-01 ,            &
         & 0.6174499520144256E-01 , 0.7138726726869340E-01/
      DATA w43b(1) , w43b(2) , w43b(3) , w43b(4) , w43b(5) , w43b(6) ,  &
         & w43b(7) , w43b(8) , w43b(9) , w43b(10) , w43b(11) , w43b(12) &
         & /0.1844477640212414E-02 , 0.1079868958589165E-01 ,           &
         & 0.2189536386779543E-01 , 0.3259746397534569E-01 ,            &
         & 0.4216313793519181E-01 , 0.5074193960018458E-01 ,            &
         & 0.5837939554261925E-01 , 0.6474640495144589E-01 ,            &
         & 0.6956619791235648E-01 , 0.7282444147183321E-01 ,            &
         & 0.7450775101417512E-01 , 0.7472214751740301E-01/
      DATA w87a(1) , w87a(2) , w87a(3) , w87a(4) , w87a(5) , w87a(6) ,  &
         & w87a(7) , w87a(8) , w87a(9) , w87a(10) , w87a(11) , w87a(12) &
         & , w87a(13) , w87a(14) , w87a(15) , w87a(16) , w87a(17) ,     &
         & w87a(18) , w87a(19) , w87a(20) , w87a(21)                    &
         & /0.8148377384149173E-02 , 0.1876143820156282E-01 ,           &
         & 0.2734745105005229E-01 , 0.3367770731163793E-01 ,            &
         & 0.3693509982042791E-01 , 0.2884872430211531E-02 ,            &
         & 0.1368594602271270E-01 , 0.2328041350288831E-01 ,            &
         & 0.3087249761171336E-01 , 0.3569363363941877E-01 ,            &
         & 0.9152833452022414E-03 , 0.5399280219300471E-02 ,            &
         & 0.1094767960111893E-01 , 0.1629873169678734E-01 ,            &
         & 0.2108156888920384E-01 , 0.2537096976925383E-01 ,            &
         & 0.2918969775647575E-01 , 0.3237320246720279E-01 ,            &
         & 0.3478309895036514E-01 , 0.3641222073135179E-01 ,            &
         & 0.3725387550304771E-01/
      DATA w87b(1) , w87b(2) , w87b(3) , w87b(4) , w87b(5) , w87b(6) ,  &
         & w87b(7) , w87b(8) , w87b(9) , w87b(10) , w87b(11) , w87b(12) &
         & , w87b(13) , w87b(14) , w87b(15) , w87b(16) , w87b(17) ,     &
         & w87b(18) , w87b(19) , w87b(20) , w87b(21) , w87b(22) ,       &
         & w87b(23)/0.2741455637620724E-03 , 0.1807124155057943E-02 ,   &
         & 0.4096869282759165E-02 , 0.6758290051847379E-02 ,            &
         & 0.9549957672201647E-02 , 0.1232944765224485E-01 ,            &
         & 0.1501044734638895E-01 , 0.1754896798624319E-01 ,            &
         & 0.1993803778644089E-01 , 0.2219493596101229E-01 ,            &
         & 0.2433914712600081E-01 , 0.2637450541483921E-01 ,            &
         & 0.2828691078877120E-01 , 0.3005258112809270E-01 ,            &
         & 0.3164675137143993E-01 , 0.3305041341997850E-01 ,            &
         & 0.3425509970422606E-01 , 0.3526241266015668E-01 ,            &
         & 0.3607698962288870E-01 , 0.3669860449845609E-01 ,            &
         & 0.3712054926983258E-01 , 0.3733422875193504E-01 ,            &
         & 0.3736107376267902E-01/
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fcentr - function value at mid point
!           absc   - abscissa
!           fval   - function value
!           savfun - array of function values which
!                    have already been computed
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
!***first executable statement  qng
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
!           test on validity of parameters
!           ------------------------------
!
      Result = 0.0E+00
      Abserr = 0.0E+00
      Neval = 0
      Ier = 6
      IF ( Epsabs>0.0E+00 .OR. Epsrel>=AMAX1(0.5E-14,0.5E+02*epmach) )  &
         & THEN
         hlgth = 0.5E+00*(B-A)
         dhlgth = ABS(hlgth)
         centr = 0.5E+00*(B+A)
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
               Abserr = ABS((res43-res21)*hlgth)
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
               Abserr = ABS((res87-res43)*hlgth)
            CASE DEFAULT
               res10 = 0.0E+00
               res21 = w21b(6)*fcentr
               resabs = w21b(6)*ABS(fcentr)
               DO k = 1 , 5
                  absc = hlgth*x1(k)
                  fval1 = F(centr+absc)
                  fval2 = F(centr-absc)
                  fval = fval1 + fval2
                  res10 = res10 + w10(k)*fval
                  res21 = res21 + w21a(k)*fval
                  resabs = resabs + w21a(k)*(ABS(fval1)+ABS(fval2))
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
                  resabs = resabs + w21b(k)*(ABS(fval1)+ABS(fval2))
                  savfun(ipx) = fval
                  fv3(k) = fval1
                  fv4(k) = fval2
               ENDDO
!
!          test for convergence.
!
               Result = res21*hlgth
               resabs = resabs*dhlgth
               reskh = 0.5E+00*res21
               resasc = w21b(6)*ABS(fcentr-reskh)
               DO k = 1 , 5
                  resasc = resasc + w21a(k)                             &
                         & *(ABS(fv1(k)-reskh)+ABS(fv2(k)-reskh))       &
                         & + w21b(k)                                    &
                         & *(ABS(fv3(k)-reskh)+ABS(fv4(k)-reskh))
               ENDDO
               Abserr = ABS((res21-res10)*hlgth)
               resasc = resasc*dhlgth
            END SELECT
            IF ( resasc/=0.0E+00 .AND. Abserr/=0.0E+00 )                &
               & Abserr = resasc*AMIN1(0.1E+01,(0.2E+03*Abserr/resasc)  &
               & **1.5E+00)
            IF ( resabs>uflow/(0.5E+02*epmach) )                        &
               & Abserr = AMAX1((epmach*0.5E+02)*resabs,Abserr)
            IF ( Abserr<=AMAX1(Epsabs,Epsrel*ABS(Result)) ) Ier = 0
! ***jump out of do-loop
            IF ( Ier==0 ) GOTO 99999
         ENDDO
      ENDIF
      CALL XERROR('abnormal return from  qng ',26,Ier,0)
99999 END
!*==QPSRT.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      SUBROUTINE QPSRT(Limit,Last,Maxerr,Ermax,Elist,Iord,Nrmax)
      IMPLICIT NONE
!*--QPSRT16994
!***begin prologue  qpsrt
!***refer to  qage,qagie,qagpe,qagse,qawce,qawse,qawoe
!***routines called  (none)
!***keywords  sequential sorting
!***description
!
! 1.        qpsrt
!           ordering routine
!              standard fortran subroutine
!              real version
!
! 2.        purpose
!              this routine maintains the descending ordering
!              in the list of the local error estimates resulting from
!              the interval subdivision process. at each call two error
!              estimates are inserted using the sequential search
!              method, top-down for the largest error estimate
!              and bottom-up for the smallest error estimate.
!
! 3.        calling sequence
!              call qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!
!           parameters (meaning at output)
!              limit  - integer
!                       maximum number of error estimates the list
!                       can contain
!
!              last   - integer
!                       number of error estimates currently
!                       in the list
!
!              maxerr - integer
!                       maxerr points to the nrmax-th largest error
!                       estimate currently in the list
!
!              ermax  - real
!                       nrmax-th largest error estimate
!                       ermax = elist(maxerr)
!
!              elist  - real
!                       vector of dimension last containing
!                       the error estimates
!
!              iord   - integer
!                       vector of dimension last, the first k
!                       elements of which contain pointers
!                       to the error estimates, such that
!                       elist(iord(1)),... , elist(iord(k))
!                       form a decreasing sequence, with
!                       k = last if last<=(limit/2+2), and
!                       k = limit+1-last otherwise
!
!              nrmax  - integer
!                       maxerr = iord(nrmax)
!
! 4.        no subroutines or functions needed
!***end prologue  qpsrt
!
      REAL Elist , Ermax , errmax , errmin
      INTEGER i , ibeg , ido , Iord , isucc , j , jbnd , jupbn , k ,    &
            & Last , Limit , Maxerr , Nrmax
      DIMENSION Elist(Last) , Iord(Last)
!
!           check whether the list contains more than
!           two error estimates.
!
!***first executable statement  qpsrt
      IF ( Last>2 ) THEN
!
!           this part of the routine is only executed
!           if, due to a difficult integrand, subdivision
!           increased the error estimate. in the normal case
!           the insert procedure should start after the
!           nrmax-th largest error estimate.
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
!           compute the number of elements in the list to
!           be maintained in descending order. this number
!           depends on the number of subdivisions still
!           allowed.
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
!*==QWGTC.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      REAL FUNCTION QWGTC(X,C,P2,P3,P4,Kp)
      IMPLICIT NONE
!*--QWGTC17135
!***begin prologue  qwgtc
!***refer to qk15w
!***routines called  (none)
!***revision date  810101   (yymmdd)
!***keywords  weight function, cauchy principal value
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this function subprogram is used together with the
!            routine qawc and defines the weight function.
!***end prologue  qwgtc
!
      REAL C , P2 , P3 , P4 , X
      INTEGER Kp
!***first executable statement
      QWGTC = 0.1E+01/(X-C)
      END
!*==QWGTF.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      REAL FUNCTION QWGTF(X,Omega,P2,P3,P4,Integr)
      IMPLICIT NONE
!*--QWGTF17155
!***begin prologue  qwgtf
!***refer to   qk15w
!***routines called  (none)
!***revision date 810101   (yymmdd)
!***keywords  cos or sin in weight function
!***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. * progr. div. - k.u.leuven
!***end prologue  qwgtf
!
      REAL Omega , omx , P2 , P3 , P4 , X
      INTEGER Integr
!***first executable statement
      omx = Omega*X
      IF ( Integr==2 ) THEN
         QWGTF = SIN(omx)
      ELSE
         QWGTF = COS(omx)
      ENDIF
      END
!*==QWGTS.spg  processed by SPAG 6.72Dc at 04:31 on  7 Dec 2021
      REAL FUNCTION QWGTS(X,A,B,Alfa,Beta,Integr)
      IMPLICIT NONE
!*--QWGTS17178
!***begin prologue  qwgts
!***refer to qk15w
!***routines called  (none)
!***revision date  810101   (yymmdd)
!***keywords  weight function, algebraico-logarithmic
!             end-point singularities
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this function subprogram is used together with the
!            routine qaws and defines the weight function.
!***end prologue  qwgts
!
      REAL A , Alfa , B , Beta , bmx , X , xma
      INTEGER Integr
!***first executable statement
      xma = X - A
      bmx = B - X
      QWGTS = xma**Alfa*bmx**Beta
      SELECT CASE (Integr)
      CASE (1)
      CASE (3)
         QWGTS = QWGTS*ALOG(bmx)
      CASE (4)
         QWGTS = QWGTS*ALOG(xma)*ALOG(bmx)
      CASE DEFAULT
         QWGTS = QWGTS*ALOG(xma)
      END SELECT
      END
