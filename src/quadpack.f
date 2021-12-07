      subroutine dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     *    limit,lenw,last,iwork,work)
c***begin prologue  dqag
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result)le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c            f      - double precision
c                     function subprogam defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accoracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                       7 - 15 points if key.lt.2,
c                      10 - 21 points if key = 2,
c                      15 - 31 points if key = 3,
c                      20 - 41 points if key = 4,
c                      25 - 51 points if key = 5,
c                      30 - 61 points if key.gt.5.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c                      error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yield no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulaties.
c                             if the position of a local difficulty can
c                             be determined (i.e.singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set
c                             to zero.
c                             except when lenw is invalid, iwork(1),
c                             work(limit*2+1) and work(limit*3+1) are
c                             set to zero, work(1) is set to a and
c                             work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdiviosion process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers to the error
c                    estimates over the subintervals, such that
c                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
c                    form a decreasing sequence with k = last if
c                    last.le.(limit/2+2), and k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left end
c                    points of the subintervals in the partition of
c                     (a,b),
c                    work(limit+1), ..., work(limit+last) contain the
c                     right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last) contain
c                     the error estimates.
c
c***references  (none)
c***routines called  dqage,xerror
c***end prologue  dqag
      double precision a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,key,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of lenw.
c
c***first executable statement  dqag
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqage.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqag ,26,ier,lvl)
      return
      end
      subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqage
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral   i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                          7 - 15 points if key.lt.2,
c                         10 - 21 points if key = 2,
c                         15 - 31 points if key = 3,
c                         20 - 41 points if key = 4,
c                         25 - 51 points if key = 5,
c                         30 - 61 points if key.gt.5.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value
c                             of limit.
c                             however, if this yields no improvement it
c                             is rather advised to analyze the integrand
c                             in order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined(e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1) ,
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the
c                      integral approximations on the subintervals
c
c            elist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the moduli of the
c                      absolute error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the
c                      error estimates over the subintervals,
c                      such that elist(iord(1)), ...,
c                      elist(iord(k)) form a decreasing sequence,
c                      with k = last if last.le.(limit/2+2), and
c                      k = limit+1-last otherwise
c
c            last    - integer
c                      number of subintervals actually produced in the
c                      subdivision process
c
c***references  (none)
c***routines called  d1mach,dqk15,dqk21,dqk31,
c                    dqk41,dqk51,dqk61,dqpsrt
c***end prologue  dqage
c
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                      (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach  is the largest relative spacing.
c           uflow  is the smallest positive magnitude.
c
c***first executable statement  dqage
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
c
c           test on accuracy.
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
      if(abserr.le.0.5d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)
     *  .or.abserr.eq.0.0d+00) go to 60
c
c           initialization
c           --------------
c
c
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 30 last = 2,limit
c
c           bisect the subinterval with the largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(dabs(rlist(maxerr)-area12).le.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 8
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*
     *  epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with the largest error estimate (to be bisected next).
c
   20   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
c
c           compute final result.
c           ---------------------
c
   40 result = 0.0d+00
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
      subroutine dqagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,
     *   ier,limit,lenw,last,iwork,work)
c***begin prologue  dqagi
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. -k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c            or i = integral of f over (-infinity,bound)
c            or i = integral of f over (-infinity,+infinity)
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration over infinite intervals
c        standard fortran subroutine
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - double precision
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - integer
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                              or limit.lt.1 or leniw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. exept when limit or leniw is
c                             invalid, iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first
c                    k elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ...,work(limit*2+last) contain the
c                     integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3)
c                     contain the error estimates.
c***references  (none)
c***routines called  dqagie,xerror
c***end prologue  dqagi
c
      double precision abserr,bound,epsabs,epsrel,f,result,work
      integer ier,inf,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqagi
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqagie.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
       lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqagi,26,ier,lvl)
      return
      end
      subroutine dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagie
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math & progr. div - k.u.leuven
c           de doncker,elise,appl. math & progr. div - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c            or i = integral of f over (-infinity,bound)
c            or i = integral of f over (-infinity,+infinity),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c integration over infinite intervals
c standard fortran subroutine
c
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - double precision
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - double precision
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however,if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1),
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to 0
c                             and 1 respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit,  the first
c                     last elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced
c                     in the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk15i,dqpsrt
c***end prologue  dqagie
      double precision abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,blist,boun,bound,b1,b2,correc,dabs,defabs,defab1,defab2,
     *  dmax1,dres,d1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs,
     *  reseps,result,res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,inf,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg.
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least (limexp+2),
c                       containing the part of the epsilon table
c                       wich is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained, it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagie
       epmach = d1mach(4)
c
c           test on validity of parameters
c           -----------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = 0.0d+00
      blist(1) = 0.1d+01
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     *  ier = 6
       if(ier.eq.6) go to 999
c
c
c           first approximation to the integral
c           -----------------------------------
c
c           determine the interval to be mapped onto (0,1).
c           if inf = 2 the integral is computed as i = i1+i2, where
c           i1 = integral of f over (-infinity,0),
c           i2 = integral of f over (0,+infinity).
c
      boun = bound
      if(inf.eq.2) boun = 0.0d+00
      call dqk15i(f,boun,inf,0.0d+00,0.1d+01,result,abserr,
     *  defabs,resabs)
c
c           test on accuracy
c
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0d+00) go to 130
c
c           initialization
c           --------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
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
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
        call dqk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2)go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at some points of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
        if(errsum.le.errbnd) go to 115
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
        if(abserr.le.ertest) go to 100
c
c            prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = 0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if((ier+ierro).eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00)go to 105
      if(abserr.gt.errsum)go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area))go to 115
c
c           test on divergence
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03.
     *or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
        result = result+rlist(k)
  120 continue
      abserr = errsum
  130 neval = 30*last-15
      if(inf.eq.2) neval = 2*neval
      if(ier.gt.2) ier=ier-1
  999 return
      end
      subroutine dqagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,
     *   neval,ier,leniw,lenw,last,iwork,work)
c***begin prologue  dqagp
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            break points of the integration interval, where local
c            difficulties of the integrand may occur (e.g.
c            singularities, discontinuities), are provided by the user.
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            npts2  - integer
c                     number equal to two more than the number of
c                     user-supplied break points within the integration
c                     range, npts.ge.2.
c                     if npts2.lt.2, the routine will end with ier = 6.
c
c            points - double precision
c                     vector of dimension npts2, the first (npts2-2)
c                     elements of which are the user provided break
c                     points. if these points do not constitute an
c                     ascending sequence there will be an automatic
c                     sorting.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (i.e. singularity,
c                             discontinuity within the interval), it
c                             should be supplied to the routine as an
c                             element of the vector points. if necessary
c                             an appropriate special-purpose integrator
c                             must be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that
c                             the returned result is the best which
c                             can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid because
c                             npts2.lt.2 or
c                             break points are specified outside
c                             the integration range or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             result, abserr, neval, last are set to
c                             zero. exept when leniw or lenw or npts2 is
c                             invalid, iwork(1), iwork(limit+1),
c                             work(limit*2+1) and work(limit*3+1)
c                             are set to zero.
c                             work(1) is set to a and work(limit+1)
c                             to b (where limit = (leniw-npts2)/2).
c
c         dimensioning parameters
c            leniw - integer
c                    dimensioning parameter for iwork
c                    leniw determines limit = (leniw-npts2)/2,
c                    which is the maximum number of subintervals in the
c                    partition of the given integration interval (a,b),
c                    leniw.ge.(3*npts2-2).
c                    if leniw.lt.(3*npts2-2), the routine will end with
c                    ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least leniw*2-npts2.
c                    if lenw.lt.leniw*2-npts2, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least leniw. on return,
c                    the first k elements of which contain
c                    pointers to the error estimates over the
c                    subintervals, such that work(limit*3+iwork(1)),...,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c                    iwork(limit+1), ...,iwork(limit+last) contain the
c                     subdivision levels of the subintervals, i.e.
c                     if (aa,bb) is a subinterval of (p1,p2)
c                     where p1 as well as p2 is a user-provided
c                     break point or integration limit, then (aa,bb) has
c                     level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
c                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have
c                     no significance for the user,
c                    note that limit = (leniw-npts2)/2.
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the corresponding error estimates,
c                    work(limit*4+1), ..., work(limit*4+npts2)
c                     contain the integration limits and the
c                     break points sorted in an ascending sequence.
c                    note that limit = (leniw-npts2)/2.
c
c***references  (none)
c***routines called  dqagpe,xerror
c***end prologue  dqagp
c
      double precision a,abserr,b,epsabs,epsrel,f,points,result,work
      integer ier,iwork,last,leniw,lenw,limit,lvl,l1,l2,l3,l4,neval,
     *  npts2
c
      dimension iwork(leniw),points(npts2),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqagp
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(leniw.lt.(3*npts2-2).or.lenw.lt.(leniw*2-npts2).or.npts2.lt.2)
     *  go to 10
c
c         prepare call for dqagpe.
c
      limit = (leniw-npts2)/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
c
      call dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),work(l4),
     *  iwork(1),iwork(l1),iwork(l2),last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqagp,26,ier,lvl)
      return
      end
      subroutine dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,
     *   abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin,
     *   last)
c***begin prologue  dqagpe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive.
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b), hopefully
c            satisfying following claim for accuracy abs(i-result).le.
c            max(epsabs,epsrel*abs(i)). break points of the integration
c            interval, where local difficulties of the integrand may
c            occur(e.g. singularities,discontinuities),provided by user.
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            npts2  - integer
c                     number equal to two more than the number of
c                     user-supplied break points within the integration
c                     range, npts2.ge.2.
c                     if npts2.lt.2, the routine will end with ier = 6.
c
c            points - double precision
c                     vector of dimension npts2, the first (npts2-2)
c                     elements of which are the user provided break
c                     points. if these points do not constitute an
c                     ascending sequence there will be an automatic
c                     sorting.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.npts2
c                     if limit.lt.npts2, the routine will end with
c                     ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (i.e. singularity,
c                             discontinuity within the interval), it
c                             should be supplied to the routine as an
c                             element of the vector points. if necessary
c                             an appropriate special-purpose integrator
c                             must be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid because
c                             npts2.lt.2 or
c                             break points are specified outside
c                             the integration range or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.npts2.
c                             result, abserr, neval, last, rlist(1),
c                             and elist(1) are set to zero. alist(1) and
c                             blist(1) are set to a and b respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            pts    - double precision
c                     vector of dimension at least npts2, containing the
c                     integration limits and the break points of the
c                     interval in ascending sequence.
c
c            level  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subinterval, i.e. if
c                     (aa,bb) is a subinterval of (p1,p2) where p1 as
c                     well as p2 is a user-provided break point or
c                     integration limit, then (aa,bb) has level l if
c                     abs(bb-aa) = abs(p2-p1)*2**(-l).
c
c            ndin   - integer
c                     vector of dimension at least npts2, after first
c                     integration over the intervals (pts(i)),pts(i+1),
c                     i = 0,1, ..., npts2-2, the error estimates over
c                     some of the intervals may have been increased
c                     artificially, in order to put their subdivision
c                     forward. if this happens for the subinterval
c                     numbered k, ndin(k) is put to 1, otherwise
c                     ndin(k) = 0.
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivisions process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagpe
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,dmax1,dmin1,
     *  dres,d1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,
     *  errmax,error1,erro12,error2,errsum,ertest,f,oflow,points,pts,
     *  resa,resabs,reseps,result,res3la,rlist,rlist2,sign,temp,uflow
      integer i,id,ier,ierro,ind1,ind2,iord,ip1,iroff1,iroff2,iroff3,j,
     *  jlow,jupbnd,k,ksgn,ktmin,last,levcur,level,levmax,limit,maxerr,
     *  ndin,neval,nint,nintp1,npts,npts2,nres,nrmax,numrl2
      logical extrap,noext
c
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  level(limit),ndin(npts2),points(npts2),pts(npts2),res3la(3),
     *  rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine epsalg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table which
c                       is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an appropriate
c                       approximation to the compounded integral has
c                       been obtained, it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one.
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation is
c                       no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagpe
      epmach = d1mach(4)
c
c            test on validity of parameters
c            -----------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      level(1) = 0
      npts = npts2-2
      if(npts2.lt.2.or.limit.le.npts.or.(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))) ier = 6
      if(ier.eq.6) go to 999
c
c            if any break points are provided, sort them into an
c            ascending sequence.
c
      sign = 1.0d+00
      if(a.gt.b) sign = -1.0d+00
      pts(1) = dmin1(a,b)
      if(npts.eq.0) go to 15
      do 10 i = 1,npts
        pts(i+1) = points(i)
   10 continue
   15 pts(npts+2) = dmax1(a,b)
      nint = npts+1
      a1 = pts(1)
      if(npts.eq.0) go to 40
      nintp1 = nint+1
      do 20 i = 1,nint
        ip1 = i+1
        do 20 j = ip1,nintp1
          if(pts(i).le.pts(j)) go to 20
          temp = pts(i)
          pts(i) = pts(j)
          pts(j) = temp
   20 continue
      if(pts(1).ne.dmin1(a,b).or.pts(nintp1).ne.dmax1(a,b)) ier = 6
      if(ier.eq.6) go to 999
c
c            compute first integral and error approximations.
c            ------------------------------------------------
c
   40 resabs = 0.0d+00
      do 50 i = 1,nint
        b1 = pts(i+1)
        call dqk21(f,a1,b1,area1,error1,defabs,resa)
        abserr = abserr+error1
        result = result+area1
        ndin(i) = 0
        if(error1.eq.resa.and.error1.ne.0.0d+00) ndin(i) = 1
        resabs = resabs+defabs
        level(i) = 0
        elist(i) = error1
        alist(i) = a1
        blist(i) = b1
        rlist(i) = area1
        iord(i) = i
        a1 = b1
   50 continue
      errsum = 0.0d+00
      do 55 i = 1,nint
        if(ndin(i).eq.1) elist(i) = abserr
        errsum = errsum+elist(i)
   55 continue
c
c           test on accuracy.
c
      last = nint
      neval = 21*nint
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      if(abserr.le.0.1d+03*epmach*resabs.and.abserr.gt.errbnd) ier = 2
      if(nint.eq.1) go to 80
      do 70 i = 1,npts
        jlow = i+1
        ind1 = iord(i)
        do 60 j = jlow,nint
          ind2 = iord(j)
          if(elist(ind1).gt.elist(ind2)) go to 60
          ind1 = ind2
          k = j
   60   continue
        if(ind1.eq.iord(i)) go to 70
        iord(k) = iord(i)
        iord(i) = ind1
   70 continue
      if(limit.lt.npts2) ier = 1
   80 if(ier.ne.0.or.abserr.le.errbnd) go to 210
c
c           initialization
c           --------------
c
      rlist2(1) = result
      maxerr = iord(1)
      errmax = elist(maxerr)
      area = result
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
      uflow = d1mach(1)
      oflow = d1mach(2)
      abserr = oflow
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*resabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 160 last = npts2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        levcur = level(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk21(f,a1,b1,area1,error1,resa,defab1)
        call dqk21(f,a2,b2,area2,error2,resa,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+42
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 95
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 90
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   90   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   95   level(maxerr) = levcur
        level(last) = levcur
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 100
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 110
  100   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
  110   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 190
c ***jump out of do-loop
        if(ier.ne.0) go to 170
        if(noext) go to 160
        erlarg = erlarg-erlast
        if(levcur+1.le.levmax) erlarg = erlarg+erro12
        if(extrap) go to 120
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(level(maxerr)+1.le.levmax) go to 160
        extrap = .true.
        nrmax = 2
  120   if(ierro.eq.3.or.erlarg.le.ertest) go to 140
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over
c           the larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 130 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(level(maxerr)+1.le.levmax) go to 160
          nrmax = nrmax+1
  130   continue
c
c           perform extrapolation.
c
  140   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.le.2) go to 155
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 150
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.lt.ertest) go to 170
c
c           prepare bisection of the smallest interval.
c
  150   if(numrl2.eq.1) noext = .true.
        if(ier.ge.5) go to 170
  155   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        levmax = levmax+1
        erlarg = errsum
  160 continue
c
c           set the final result.
c           ---------------------
c
c
  170 if(abserr.eq.oflow) go to 190
      if((ier+ierro).eq.0) go to 180
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00)go to 175
      if(abserr.gt.errsum)go to 190
      if(area.eq.0.0d+00) go to 210
      go to 180
  175 if(abserr/dabs(result).gt.errsum/dabs(area))go to 190
c
c           test on divergence.
c
  180 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     *  resabs*0.1d-01) go to 210
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03.or.
     *  errsum.gt.dabs(area)) ier = 6
      go to 210
c
c           compute global integral sum.
c
  190 result = 0.0d+00
      do 200 k = 1,last
        result = result+rlist(k)
  200 continue
      abserr = errsum
  210 if(ier.gt.2) ier = ier-1
      result = result*sign
  999 return
      end
      subroutine dqags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)
c***begin prologue  dqags
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end-point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & prog. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral  i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account. however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero.except when limit or lenw is invalid,
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, detemines the
c                    number of significant elements actually in the work
c                    arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end-points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end-points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c***references  (none)
c***routines called  dqagse,xerror
c***end prologue  dqags
c
c
      double precision a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqags
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqagse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqags,26,ier,lvl)
      return
      end
      subroutine dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagse
c
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,d1mach,dmax1,
     *  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax,
     *  error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result,
     *  res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     * res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2 containing
c                       the part of the epsilon table which is still
c                       needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagse
      epmach = d1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     *   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      ierro = 0
      call dqk21(f,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0d+00) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
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
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk21(f,a1,b1,area1,error1,resabs,defab1)
        call dqk21(f,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = dabs(b-a)*0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     * .or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
      subroutine dqawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)
c***begin prologue  dqawc
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,j4
c***keywords  automatic integrator, special-purpose,
c             cauchy principal value,
c             clenshaw-curtis, globally adaptive
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a
c            cauchy principal value i = integral of f*w over (a,b)
c            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
c            following claim for accuracy
c            abs(i-result).le.max(epsabe,epsrel*abs(i)).
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     under limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            c      - parameter in the weight function, c.ne.a, c.ne.b.
c                     if c = a or c = b, the routine will end with
c                     ier = 6 .
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate or the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. exept when lenw or limit is invalid,
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c           lenw   - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)), ... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - double precision
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c***references  (none)
c***routines called  dqawce,xerror
c***end prologue  dqawc
c
      double precision a,abserr,b,c,epsabs,epsrel,f,result,work
      integer ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqawc
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqawce.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      call dqawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,
     *  work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqawc,26,ier,lvl)
      return
      end
      subroutine dqawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,
     *   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqawce
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,j4
c***keywords  automatic integrator, special-purpose,
c             cauchy principal value, clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***  purpose  the routine calculates an approximation result to a
c              cauchy principal value i = integral of f*w over (a,b)
c              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
c              following claim for accuracy
c              abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            c      - double precision
c                     parameter in the weight function, c.ne.a, c.ne.b
c                     if c = a or c = b, the routine will end with
c                     ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement it is advised to analyze the
c                             the integrand, in order to determine the
c                             the integration difficulties. if the
c                             position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some interior points of
c                             the integration interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - double precision
c                      vector of dimension at least limit, the first
c                       last  elements of which are the integral
c                      approximations on the subintervals
c
c            elist   - double precision
c                      vector of dimension limit, the first  last
c                      elements of which are the moduli of the absolute
c                      error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the error
c                      estimates over the subintervals, so that
c                      elist(iord(1)), ..., elist(iord(k)) with k = last
c                      if last.le.(limit/2+2), and k = limit+1-last
c                      otherwise, form a decreasing sequence
c
c            last    - integer
c                      number of subintervals actually produced in
c                      the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqc25c,dqpsrt
c***end prologue  dqawce
c
      double precision a,aa,abserr,alist,area,area1,area12,area2,a1,a2,
     *  b,bb,blist,b1,b2,c,dabs,dmax1,d1mach,elist,epmach,epsabs,epsrel,
     *  errbnd,errmax,error1,erro12,error2,errsum,f,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,krule,last,limit,maxerr,nev,
     *  neval,nrmax
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqawce
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(c.eq.a.or.c.eq.b.or.(epsabs.le.0.0d+00.and
     *  .epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      aa=a
      bb=b
      if (a.le.b) go to 10
      aa=b
      bb=a
10    ier=0
      krule = 1
      call dqc25c(f,aa,bb,c,result,abserr,krule,neval)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      alist(1) = a
      blist(1) = b
c
c           test on accuracy
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
      if(limit.eq.1) ier = 1
      if(abserr.lt.dmin1(0.1d-01*dabs(result),errbnd)
     *  .or.ier.eq.1) go to 70
c
c           initialization
c           --------------
c
      alist(1) = aa
      blist(1) = bb
      rlist(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 40 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        b2 = blist(maxerr)
        if(c.le.b1.and.c.gt.a1) b1 = 0.5d+00*(c+b2)
        if(c.gt.b1.and.c.lt.b2) b1 = 0.5d+00*(a1+c)
        a2 = b1
        krule = 2
        call dqc25c(f,a1,b1,c,area1,error1,krule,nev)
        neval = neval+nev
        call dqc25c(f,a2,b2,c,area2,error2,krule,nev)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(dabs(rlist(maxerr)-area12).lt.0.1d-04*dabs(area12)
     *    .and.erro12.ge.0.99d+00*errmax.and.krule.eq.0)
     *    iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax.and.krule.eq.0)
     *    iroff2 = iroff2+1
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 15
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1.ge.6.and.iroff2.gt.20) ier = 2
c
c           set error flag in the case that number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)
     *    *(dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   15   if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30    call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 50
   40 continue
c
c           compute final result.
c           ---------------------
c
   50 result = 0.0d+00
      do 60 k=1,last
        result = result+rlist(k)
   60 continue
      abserr = errsum
   70 if (aa.eq.b) result=-result
  999 return
      end
      subroutine dqawf(f,a,omega,integr,epsabs,result,abserr,neval,ier,
     *   limlst,lst,leniw,maxp1,lenw,iwork,work)
c***begin prologue  dqawf
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,fourier
c             integral, integration between zeros with dqawoe,
c             convergence acceleration with dqelg
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            fourier integral i=integral of f(x)*w(x) over (a,infinity)
c            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested, epsabs.gt.0.
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of cycles allowed
c                             has been achieved, i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account). examine the array iwork which
c                             contains the error flags on the cycles, in
c                             order to look for eventual local
c                             integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence accelaration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy.
c                             as in the case of ier = 1, it is advised
c                             to examine the array iwork which contains
c                             the error flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.1 or
c                              leniw.lt.(limlst+2) or maxp1.lt.1 or
c                              lenw.lt.(leniw*2+maxp1*25).
c                              result, abserr, neval, lst are set to
c                              zero.
c                         = 7 bad integrand behaviour occurs within
c                             one or more of the cycles. location and
c                             type of the difficulty involved can be
c                             determined from the first lst elements of
c                             vector iwork.  here lst is the number of
c                             cycles actually needed (see below).
c                             iwork(k) = 1 the maximum number of
c                                          subdivisions (=(leniw-limlst)
c                                          /2) has been achieved on the
c                                          k th cycle.
c                                      = 2 occurrence of roundoff error
c                                          is detected and prevents the
c                                          tolerance imposed on the k th
c                                          cycle, from being achieved
c                                          on this cycle.
c                                      = 3 extremely bad integrand
c                                          behaviour occurs at some
c                                          points of the k th cycle.
c                                      = 4 the integration procedure
c                                          over the k th cycle does
c                                          not converge (to within the
c                                          required accuracy) due to
c                                          roundoff in the extrapolation
c                                          procedure invoked on this
c                                          cycle. it is assumed that the
c                                          result on this interval is
c                                          the best which can be
c                                          obtained.
c                                      = 5 the integral over the k th
c                                          cycle is probably divergent
c                                          or slowly convergent. it must
c                                          be noted that divergence can
c                                          occur with any other value of
c                                          iwork(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie,
c                    and ier = iwork(1) (with meaning as described
c                    for iwork(k),k = 1).
c
c         dimensioning parameters
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.3.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            lst    - integer
c                     on return, lst indicates the number of cycles
c                     actually needed for the integration.
c                     if omega = 0, then lst is set to 1.
c
c            leniw  - integer
c                     dimensioning parameter for iwork. on entry,
c                     (leniw-limlst)/2 equals the maximum number of
c                     subintervals allowed in the partition of each
c                     cycle, leniw.ge.(limlst+2).
c                     if leniw.lt.(limlst+2), the routine will end with
c                     ier = 6.
c
c            maxp1  - integer
c                     maxp1 gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e. for
c                     the intervals of lengths abs(b-a)*2**(-l),
c                     l = 0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, iwork(k) for k = 1, 2, ..., lst
c                     contain the error flags on the cycles.
c
c            work   - double precision
c                     vector of dimension at least
c                     on return,
c                     work(1), ..., work(lst) contain the integral
c                      approximations over the cycles,
c                     work(limlst+1), ..., work(limlst+lst) contain
c                      the error extimates over the cycles.
c                     further elements of work have no specific
c                     meaning for the user.
c
c***references  (none)
c***routines called  dqawfe,xerror
c***end prologue  dqawf
c
       double precision a,abserr,epsabs,f,omega,result,work
       integer ier,integr,iwork,last,leniw,lenw,limit,limlst,ll2,lvl,
     *  lst,l1,l2,l3,l4,l5,l6,maxp1,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of limlst, leniw, maxp1 and lenw.
c
c***first executable statement  dqawf
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limlst.lt.3.or.leniw.lt.(limlst+2).or.maxp1.lt.1.or.lenw.lt.
     *   (leniw*2+maxp1*25)) go to 10
c
c         prepare call for dqawfe
c
      limit = (leniw-limlst)/2
      l1 = limlst+1
      l2 = limlst+l1
      l3 = limit+l2
      l4 = limit+l3
      l5 = limit+l4
      l6 = limit+l5
      ll2 = limit+l1
      call dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,result,
     *  abserr,neval,ier,work(1),work(l1),iwork(1),lst,work(l2),
     *  work(l3),work(l4),work(l5),iwork(l1),iwork(ll2),work(l6))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqawf,26,ier,lvl)
      return
      end
      subroutine dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,
     *   result,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)
c***begin prologue  dqawfe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,
c             fourier integrals,
c             integration between zeros with dqawoe,
c             convergence acceleration with dqelg
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a
c            given fourier integal
c            i = integral of f(x)*w(x) over (a,infinity)
c            where w(x)=cos(omega*x) or w(x)=sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to
c                     be declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            omega  - double precision
c                     parameter in the weight function
c
c            integr - integer
c                     indicates which weight function is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested, epsabs.gt.0
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.1.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     allowed in the partition of each cycle, limit.ge.1
c                     each cycle, limit.ge.1.
c
c            maxp1  - integer
c                     gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e.
c                     for the intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c
c         on return
c            result - double precision
c                     approximation to the integral x
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - ier = 0 normal and reliable termination of
c                             the routine. it is assumed that the
c                             requested accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine. the
c                             estimates for integral and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of  cycles  allowed
c                             has been achieved., i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account).
c                             examine the array iwork which contains
c                             the error flags on the cycles, in order to
c                             look for eventual local integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling appropriate integrators on
c                             the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence acceleration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy. as in the case of
c                             ier = 1, it is advised to examine the
c                             array iwork which contains the error
c                             flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.3.
c                              result, abserr, neval, lst are set
c                              to zero.
c                         = 7 bad integrand behaviour occurs within one
c                             or more of the cycles. location and type
c                             of the difficulty involved can be
c                             determined from the vector ierlst. here
c                             lst is the number of cycles actually
c                             needed (see below).
c                             ierlst(k) = 1 the maximum number of
c                                           subdivisions (= limit) has
c                                           been achieved on the k th
c                                           cycle.
c                                       = 2 occurrence of roundoff error
c                                           is detected and prevents the
c                                           tolerance imposed on the
c                                           k th cycle, from being
c                                           achieved.
c                                       = 3 extremely bad integrand
c                                           behaviour occurs at some
c                                           points of the k th cycle.
c                                       = 4 the integration procedure
c                                           over the k th cycle does
c                                           not converge (to within the
c                                           required accuracy) due to
c                                           roundoff in the
c                                           extrapolation procedure
c                                           invoked on this cycle. it
c                                           is assumed that the result
c                                           on this interval is the
c                                           best which can be obtained.
c                                       = 5 the integral over the k th
c                                           cycle is probably divergent
c                                           or slowly convergent. it
c                                           must be noted that
c                                           divergence can occur with
c                                           any other value of
c                                           ierlst(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie
c                    and ier = ierlst(1) (with meaning as described
c                    for ierlst(k), k = 1).
c
c            rslst  - double precision
c                     vector of dimension at least limlst
c                     rslst(k) contains the integral contribution
c                     over the interval (a+(k-1)c,a+kc) where
c                     c = (2*int(abs(omega))+1)*pi/abs(omega),
c                     k = 1, 2, ..., lst.
c                     note that, if omega = 0, rslst(1) contains
c                     the value of the integral over (a,infinity).
c
c            erlst  - double precision
c                     vector of dimension at least limlst
c                     erlst(k) contains the error estimate corresponding
c                     with rslst(k).
c
c            ierlst - integer
c                     vector of dimension at least limlst
c                     ierlst(k) contains the error flag corresponding
c                     with rslst(k). for the meaning of the local error
c                     flags see description of output parameter ier.
c
c            lst    - integer
c                     number of subintervals needed for the integration
c                     if omega = 0 then lst is set to 1.
c
c            alist, blist, rlist, elist - double precision
c                     vector of dimension at least limit,
c
c            iord, nnlog - integer
c                     vector of dimension at least limit, providing
c                     space for the quantities needed in the subdivision
c                     process of each cycle
c
c            chebmo - double precision
c                     array of dimension at least (maxp1,25), providing
c                     space for the chebyshev moments needed within the
c                     cycles
c
c***references  (none)
c***routines called  d1mach,dqagie,dqawoe,dqelg
c***end prologue  dqawfe
c
      double precision a,abseps,abserr,alist,blist,chebmo,correc,cycle,
     *  c1,c2,dabs,dl,dla,dmax1,drl,d1mach,elist,erlst,ep,eps,epsa,
     *  epsabs,errsum,f,fact,omega,p,pi,p1,psum,reseps,result,res3la,
     *  rlist,rslst,uflow
      integer ier,ierlst,integr,iord,ktmin,l,last,lst,limit,limlst,ll,
     *    maxp1,momcom,nev,neval,nnlog,nres,numrl2
c
      dimension alist(limit),blist(limit),chebmo(maxp1,25),elist(limit),
     *  erlst(limlst),ierlst(limlst),iord(limit),nnlog(limit),psum(52),
     *  res3la(3),rlist(limit),rslst(limlst)
c
      external f
c
c
c            the dimension of  psum  is determined by the value of
c            limexp in subroutine dqelg (psum must be of dimension
c            (limexp+2) at least).
c
c           list of major variables
c           -----------------------
c
c           c1, c2    - end points of subinterval (of length cycle)
c           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
c           psum      - vector of dimension at least (limexp+2)
c                       (see routine dqelg)
c                       psum contains the part of the epsilon table
c                       which is still needed for further computations.
c                       each element of psum is a partial sum of the
c                       series which should sum to the value of the
c                       integral.
c           errsum    - sum of error estimates over the subintervals,
c                       calculated cumulatively
c           epsa      - absolute tolerance requested over current
c                       subinterval
c           chebmo    - array containing the modified chebyshev
c                       moments (see also routine dqc25f)
c
      data p/0.9d+00/
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c
c           test on validity of parameters
c           ------------------------------
c
c***first executable statement  dqawfe
      result = 0.0d+00
      abserr = 0.0d+00
      neval = 0
      lst = 0
      ier = 0
      if((integr.ne.1.and.integr.ne.2).or.epsabs.le.0.0d+00.or.
     *  limlst.lt.3) ier = 6
      if(ier.eq.6) go to 999
      if(omega.ne.0.0d+00) go to 10
c
c           integration by dqagie if omega is zero
c           --------------------------------------
c
      if(integr.eq.1) call dqagie(f,0.0d+00,1,epsabs,0.0d+00,limit,
     *  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      rslst(1) = result
      erlst(1) = abserr
      ierlst(1) = ier
      lst = 1
      go to 999
c
c           initializations
c           ---------------
c
   10 l = dabs(omega)
      dl = 2*l+1
      cycle = dl*pi/dabs(omega)
      ier = 0
      ktmin = 0
      neval = 0
      numrl2 = 0
      nres = 0
      c1 = a
      c2 = cycle+a
      p1 = 0.1d+01-p
      uflow = d1mach(1)
      eps = epsabs
      if(epsabs.gt.uflow/p1) eps = epsabs*p1
      ep = eps
      fact = 0.1d+01
      correc = 0.0d+00
      abserr = 0.0d+00
      errsum = 0.0d+00
c
c           main do-loop
c           ------------
c
      do 50 lst = 1,limlst
c
c           integrate over current subinterval.
c
        dla = lst
        epsa = eps*fact
        call dqawoe(f,c1,c2,omega,integr,epsa,0.0d+00,limit,lst,maxp1,
     *  rslst(lst),erlst(lst),nev,ierlst(lst),last,alist,blist,rlist,
     *  elist,iord,nnlog,momcom,chebmo)
        neval = neval+nev
        fact = fact*p
        errsum = errsum+erlst(lst)
        drl = 0.5d+02*dabs(rslst(lst))
c
c           test on accuracy with partial sum
c
        if((errsum+drl).le.epsabs.and.lst.ge.6) go to 80
        correc = dmax1(correc,erlst(lst))
        if(ierlst(lst).ne.0) eps = dmax1(ep,correc*p1)
        if(ierlst(lst).ne.0) ier = 7
        if(ier.eq.7.and.(errsum+drl).le.correc*0.1d+02.and.
     *  lst.gt.5) go to 80
        numrl2 = numrl2+1
        if(lst.gt.1) go to 20
        psum(1) = rslst(1)
        go to 40
   20   psum(numrl2) = psum(ll)+rslst(lst)
        if(lst.eq.2) go to 40
c
c           test on maximum number of subintervals
c
        if(lst.eq.limlst) ier = 1
c
c           perform new extrapolation
c
        call dqelg(numrl2,psum,reseps,abseps,res3la,nres)
c
c           test whether extrapolated result is influenced by roundoff
c
        ktmin = ktmin+1
        if(ktmin.ge.15.and.abserr.le.0.1d-02*(errsum+drl)) ier = 4
        if(abseps.gt.abserr.and.lst.ne.3) go to 30
        abserr = abseps
        result = reseps
        ktmin = 0
c
c           if ier is not 0, check whether direct result (partial sum)
c           or extrapolated result yields the best integral
c           approximation
c
        if((abserr+0.1d+02*correc).le.epsabs.or.
     *  (abserr.le.epsabs.and.0.1d+02*correc.ge.epsabs)) go to 60
   30   if(ier.ne.0.and.ier.ne.7) go to 60
   40   ll = numrl2
        c1 = c2
        c2 = c2+cycle
   50 continue
c
c         set final result and error estimate
c         -----------------------------------
c
   60 abserr = abserr+0.1d+02*correc
      if(ier.eq.0) go to 999
      if(result.ne.0.0d+00.and.psum(numrl2).ne.0.0d+00) go to 70
      if(abserr.gt.errsum) go to 80
      if(psum(numrl2).eq.0.0d+00) go to 999
   70 if(abserr/dabs(result).gt.(errsum+drl)/dabs(psum(numrl2)))
     *  go to 80
      if(ier.ge.1.and.ier.ne.7) abserr = abserr+drl
      go to 999
   80 result = psum(numrl2)
      abserr = errsum+drl
  999 return
      end
      subroutine dqawo(f,a,b,omega,integr,epsabs,epsrel,result,abserr,
     *   neval,ier,leniw,maxp1,lenw,last,iwork,work)
c***begin prologue  dqawo
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i=integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x)
c            or w(x) = sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the function
c                     f(x).  the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if epsabs.le.0 and
c                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of  integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             (= leniw/2) has been achieved. one can
c                             allow more subdivisions by increasing the
c                             value of leniw (and taking the according
c                             dimension adjustments into account).
c                             however, if this yields no improvement it
c                             is advised to analyze the integrand in
c                             order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some interior points of the
c                             integration interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be achieved
c                             due to roundoff in the extrapolation
c                             table, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2),
c                             or leniw.lt.2 or maxp1.lt.1 or
c                             lenw.lt.leniw*2+maxp1*25.
c                             result, abserr, neval, last are set to
c                             zero. except when leniw, maxp1 or lenw are
c                             invalid, work(limit*2+1), work(limit*3+1),
c                             iwork(1), iwork(limit+1) are set to zero,
c                             work(1) is set to a and work(limit+1) to
c                             b.
c
c         dimensioning parameters
c            leniw  - integer
c                     dimensioning parameter for iwork.
c                     leniw/2 equals the maximum number of subintervals
c                     allowed in the partition of the given integration
c                     interval (a,b), leniw.ge.2.
c                     if leniw.lt.2, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of subintervals
c                     produced in the subdivision process, which
c                     determines the number of significant elements
c                     actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, the first k elements of which contain
c                     pointers to the error estimates over the
c                     subintervals, such that work(limit*3+iwork(1)), ..
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence, with limit = lenw/2 , and k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise.
c                     furthermore, iwork(limit+1), ..., iwork(limit+
c                     last) indicate the subdivision levels of the
c                     subintervals, such that iwork(limit+i) = l means
c                     that the subinterval numbered i is of length
c                     abs(b-a)*2**(1-l).
c
c            work   - double precision
c                     vector of dimension at least lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last) contain
c                      the integral approximations over the
c                      subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c                     work(limit*4+1), ..., work(limit*4+maxp1*25)
c                      provide space for storing the chebyshev moments.
c                     note that limit = lenw/2.
c
c***references  (none)
c***routines called  dqawoe,xerror
c***end prologue  dqawo
c
       double precision a,abserr,b,epsabs,epsrel,f,omega,result,work
       integer ier,integr,iwork,last,limit,lenw,leniw,lvl,l1,l2,l3,l4,
     *  maxp1,momcom,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of leniw, maxp1 and lenw.
c
c***first executable statement  dqawo
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(leniw.lt.2.or.maxp1.lt.1.or.lenw.lt.(leniw*2+maxp1*25))
     *  go to 10
c
c         prepare call for dqawoe
c
      limit = leniw/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
      call dqawoe(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result,
     *   abserr,neval,ier,last,work(1),work(l1),work(l2),work(l3),
     *   iwork(1),iwork(l1),momcom,work(l4))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 0
      if(ier.ne.0) call xerror(26habnormal return from dqawo,26,ier,lvl)
      return
      end
      subroutine dqawoe (f,a,b,omega,integr,epsabs,epsrel,limit,icall,
     *  maxp1,result,abserr,neval,ier,last,alist,blist,rlist,elist,iord,
     *   nnlog,momcom,chebmo)
c***begin prologue  dqawoe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral
c            i = integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration
c
c            omega  - double precision
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is to be
c                     used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1 and integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subdivisions
c                     in the partition of (a,b), limit.ge.1.
c
c            icall  - integer
c                     if dqawoe is to be used only once, icall must
c                     be set to 1.  assume that during this call, the
c                     chebyshev moments (for clenshaw-curtis integration
c                     of degree 24) have been computed for intervals of
c                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
c                     if icall.gt.1 this means that dqawoe has been
c                     called twice or more on intervals of the same
c                     length abs(b-a). the chebyshev moments already
c                     computed are then re-used in subsequent calls.
c                     if icall.lt.1, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lenghts abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the
c                             requested accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved due to
c                             roundoff in the extrapolation table,
c                             and that the returned result is the
c                             best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2) or
c                             icall.lt.1 or maxp1.lt.1.
c                             result, abserr, neval, last, rlist(1),
c                             elist(1), iord(1) and nnlog(1) are set
c                             to zero. alist(1) and blist(1) are set
c                             to a and b respectively.
c
c            last  -  integer
c                     on return, last equals the number of
c                     subintervals produces in the subdivision
c                     process, which determines the number of
c                     significant elements actually in the
c                     work arrays.
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the error
c                     estimates over the subintervals,
c                     such that elist(iord(1)), ...,
c                     elist(iord(k)) form a decreasing sequence, with
c                     k = last if last.le.(limit/2+2), and
c                     k = limit+1-last otherwise.
c
c            nnlog  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subintervals, i.e.
c                     iwork(i) = l means that the subinterval
c                     numbered i is of length abs(b-a)*2**(1-l)
c
c         on entry and return
c            momcom - integer
c                     indicating that the chebyshev moments
c                     have been computed for intervals of lengths
c                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
c                     momcom.lt.maxp1
c
c            chebmo - double precision
c                     array of dimension (maxp1,25) containing the
c                     chebyshev moments
c
c***references  (none)
c***routines called  d1mach,dqc25f,dqelg,dqpsrt
c***end prologue  dqawoe
c
      double precision a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,chebmo,correc,dabs,defab1,defab2,defabs,dmax1,
     *  domega,d1mach,dres,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,erro12,error2,errsum,ertest,f,oflow,
     *  omega,resabs,reseps,result,res3la,rlist,rlist2,small,uflow,width
      integer icall,id,ier,ierro,integr,iord,iroff1,iroff2,iroff3,
     *  jupbnd,k,ksgn,ktmin,last,limit,maxerr,maxp1,momcom,nev,neval,
     *  nnlog,nres,nrmax,nrmom,numrl2
      logical extrap,noext,extall
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),rlist2(52),res3la(3),chebmo(maxp1,25),nnlog(limit)
c
      external f
c
c            the dimension of rlist2 is determined by  the value of
c            limexp in subroutine dqelg (rlist2 should be of
c            dimension (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table
c                       which is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an appropriate
c                       approximation to the compounded integral has
c                       been obtained it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one
c           small     - length of the smallest interval considered
c                       up to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation, i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true  value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqawoe
      epmach = d1mach(4)
c
c         test on validity of parameters
c         ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      nnlog(1) = 0
      if((integr.ne.1.and.integr.ne.2).or.(epsabs.le.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.icall.lt.1.or.
     *  maxp1.lt.1) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      domega = dabs(omega)
      nrmom = 0
      if (icall.gt.1) go to 5
      momcom = 0
    5 call dqc25f(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr,
     *  neval,defabs,resabs,momcom,chebmo)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.0.1d+03*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.abserr.le.errbnd) go to 200
c
c           initializations
c           ---------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ktmin = 0
      small = dabs(b-a)*0.75d+00
      nres = 0
      numrl2 = 0
      extall = .false.
      if(0.5d+00*dabs(b-a)*domega.gt.0.2d+01) go to 10
      numrl2 = 1
      extall = .true.
      rlist2(1) = result
   10 if(0.25d+00*dabs(b-a)*domega.le.0.2d+01) extall = .true.
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 140 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        nrmom = nnlog(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqc25f(f,a1,b1,domega,integr,nrmom,maxp1,0,
     *  area1,error1,nev,resabs,defab1,momcom,chebmo)
        neval = neval+nev
        call dqc25f(f,a2,b2,domega,integr,nrmom,maxp1,1,
     *  area2,error2,nev,resabs,defab2,momcom,chebmo)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 25
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 20
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   20   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   25   rlist(maxerr) = area1
        rlist(last) = area2
        nnlog(maxerr) = nrmom
        nnlog(last) = nrmom
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)
     *  *(dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 30
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 40
   30   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to bisected next).
c
   40   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
      if(errsum.le.errbnd) go to 170
      if(ier.ne.0) go to 150
        if(last.eq.2.and.extall) go to 120
        if(noext) go to 140
        if(.not.extall) go to 50
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 70
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
   50   width = dabs(blist(maxerr)-alist(maxerr))
        if(width.gt.small) go to 140
        if(extall) go to 60
c
c           test whether we can start with the extrapolation procedure
c           (we do this if we integrate over the next interval with
c           use of a gauss-kronrod rule - see subroutine dqc25f).
c
        small = small*0.5d+00
        if(0.25d+00*width*domega.gt.0.2d+01) go to 140
        extall = .true.
        go to 130
   60   extrap = .true.
        nrmax = 2
   70   if(ierro.eq.3.or.erlarg.le.ertest) go to 90
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over
c           the larger intervals (erlarg) and perform extrapolation.
c
        jupbnd = last
        if (last.gt.(limit/2+2)) jupbnd = limit+3-last
        id = nrmax
        do 80 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 140
          nrmax = nrmax+1
   80   continue
c
c           perform extrapolation.
c
   90   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.lt.3) go to 110
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 100
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 150
c
c           prepare bisection of the smallest interval.
c
  100   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 150
  110   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 140
  120   small = small*0.5d+00
        numrl2 = numrl2+1
        rlist2(numrl2) = area
  130   ertest = errbnd
        erlarg = errsum
  140 continue
c
c           set the final result.
c           ---------------------
c
  150 if(abserr.eq.oflow.or.nres.eq.0) go to 170
      if(ier+ierro.eq.0) go to 165
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 160
      if(abserr.gt.errsum) go to 170
      if(area.eq.0.0d+00) go to 190
      go to 165
  160 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 170
c
c           test on divergence.
c
  165 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     * defabs*0.1d-01) go to 190
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     * .or.errsum.ge.dabs(area)) ier = 6
      go to 190
c
c           compute global integral sum.
c
  170 result = 0.0d+00
      do 180 k=1,last
        result = result+rlist(k)
  180 continue
      abserr = errsum
  190 if (ier.gt.2) ier=ier-1
  200 if (integr.eq.2.and.omega.lt.0.0d+00) result=-result
  999 return
      end
      subroutine dqaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,
     *   abserr,neval,ier,limit,lenw,last,iwork,work)
c***begin prologue  dqaws
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end-point singularities,
c             clenshaw-curtis, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. -k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - double precision
c                     parameter in the integrand function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - double precision
c                     parameter in the integrand function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal tepination of the routine
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties
c                             which prevent the requested tolerance from
c                             being achieved. in case of a jump
c                             discontinuity or a local singularity
c                             of algebraico-logarithmic type at one or
c                             more interior points of the integration
c                             range, one should proceed by splitting up
c                             the interval at these points and calling
c                             the integrator on the subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1) or
c                             or integr.lt.1 or integr.gt.4 or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.2 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. except when lenw or limit is invalid
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit  - integer
c                     dimensioning parameter for iwork
c                     limit determines the maximum number of
c                     subintervals in the partition of the given
c                     integration interval (a,b), limit.ge.2.
c                     if limit.lt.2, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least limit*4.
c                     if lenw.lt.limit*4, the routine will end
c                     with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of
c                     subintervals produced in the subdivision process,
c                     which determines the significant number of
c                     elements actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension limit, the first k
c                     elements of which contain pointers
c                     to the error estimates over the subintervals,
c                     such that work(limit*3+iwork(1)), ...,
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence with k = last if last.le.(limit/2+2),
c                     and k = limit+1-last otherwise
c
c            work   - double precision
c                     vector of dimension lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last)
c                      contain the integral approximations over
c                      the subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c
c***references  (none)
c***routines called  dqawse,xerror
c***end prologue  dqaws
c
      double precision a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
      integer ier,integr,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqaws
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.2.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqawse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result,
     *  abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from dqaws,26,ier,lvl)
      return
      end
      subroutine dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,
     *   result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqawse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end point singularities,
c             clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points,
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - double precision
c                     parameter in the weight function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - double precision
c                     parameter in the weight function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.2
c                     if limit.lt.2, the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement, it is advised to analyze the
c                             integrand in order to determine the
c                             integration difficulties which prevent the
c                             requested tolerance from being achieved.
c                             in case of a jump discontinuity or a local
c                             singularity of algebraico-logarithmic type
c                             at one or more interior points of the
c                             integration range, one should proceed by
c                             splitting up the interval at these
c                             points and calling the integrator on the
c                             subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1), or
c                             integr.lt.1 or integr.gt.4, or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             or limit.lt.2.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit,the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     of which are pointers to the error
c                     estimates over the subintervals, so that
c                     elist(iord(1)), ..., elist(iord(k)) with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise form a decreasing sequence
c
c            last   - integer
c                     number of subintervals actually produced in
c                     the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqc25s,dqmomo,dqpsrt
c***end prologue  dqawse
c
      double precision a,abserr,alfa,alist,area,area1,area12,area2,a1,
     *  a2,b,beta,blist,b1,b2,centre,dabs,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,erro12,error2,errsum,f,
     *  resas1,resas2,result,rg,rh,ri,rj,rlist,uflow
      integer ier,integr,iord,iroff1,iroff2,k,last,limit,maxerr,nev,
     *  neval,nrmax
c
      external f
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),ri(25),rj(25),rh(25),rg(25)
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqawse
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(b.le.a.or.(epsabs.eq.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.alfa.le.(-0.1d+01).
     *  or.beta.le.(-0.1d+01).or.integr.lt.1.or.integr.gt.4.or.
     *  limit.lt.2) go to 999
      ier = 0
c
c           compute the modified chebyshev moments.
c
      call dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c
c           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
c
      centre = 0.5d+00*(b+a)
      call dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
      neval = nev
      call dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
      last = 2
      neval = neval+nev
      result = area1+area2
      abserr = error1+error2
c
c           test on accuracy.
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
c
c           initialization
c           --------------
c
      if(error2.gt.error1) go to 10
      alist(1) = a
      alist(2) = centre
      blist(1) = centre
      blist(2) = b
      rlist(1) = area1
      rlist(2) = area2
      elist(1) = error1
      elist(2) = error2
      go to 20
   10 alist(1) = centre
      alist(2) = a
      blist(1) = b
      blist(2) = centre
      rlist(1) = area2
      rlist(2) = area1
      elist(1) = error2
      elist(2) = error1
   20 iord(1) = 1
      iord(2) = 2
      if(limit.eq.2) ier = 1
      if(abserr.le.errbnd.or.ier.eq.1) go to 999
      errmax = elist(1)
      maxerr = 1
      nrmax = 1
      area = result
      errsum = abserr
      iroff1 = 0
      iroff2 = 0
c
c            main do-loop
c            ------------
c
      do 60 last = 3,limit
c
c           bisect the subinterval with largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
c
        call dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
        neval = neval+nev
        call dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
        neval = neval+nev
c
c           improve previous approximations integral and error
c           and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(a.eq.a1.or.b.eq.b2) go to 30
        if(resas1.eq.error1.or.resas2.eq.error2) go to 30
c
c           test for roundoff error.
c
        if(dabs(rlist(maxerr)-area12).lt.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
   30   rlist(maxerr) = area1
        rlist(last) = area2
c
c           test on accuracy.
c
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 35
c
c           set error flag in the case that the number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c
c           set error flag in the case of roundoff error.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case of bad integrand behaviour
c           at interior points of integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   35   if(error2.gt.error1) go to 40
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 50
   40   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with largest error estimate (to be bisected next).
c
   50   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if (ier.ne.0.or.errsum.le.errbnd) go to 70
   60 continue
c
c           compute final result.
c           ---------------------
c
   70 result = 0.0d+00
      do 80 k=1,last
        result = result+rlist(k)
   80 continue
      abserr = errsum
  999 return
      end
      subroutine dqc25c(f,a,b,c,result,abserr,krul,neval)
c***begin prologue  dqc25c
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2,j4
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b) with
c            error estimate, where w(x) = 1/(x-c)
c***description
c
c        integration rules for the computation of cauchy
c        principal value integrals
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           f      - double precision
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - double precision
c                    left end point of the integration interval
c
c           b      - double precision
c                    right end point of the integration interval, b.gt.a
c
c           c      - double precision
c                    parameter in the weight function
c
c           result - double precision
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if c lies within ten percent
c                    of the integration interval. in the other case the
c                    15-point kronrod rule obtained by optimal addition
c                    of abscissae to the 7-point gauss rule, is applied.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           krul   - integer
c                    key which is decreased by 1 if the 15-point
c                    gauss-kronrod scheme has been used
c
c           neval  - integer
c                    number of integrand evaluations
c
c.......................................................................
c***references  (none)
c***routines called  dqcheb,dqk15w,dqwgtc
c***end prologue  dqc25c
c
      double precision a,abserr,ak22,amom0,amom1,amom2,b,c,cc,centr,
     *  cheb12,cheb24,dabs,dlog,dqwgtc,f,fval,hlgth,p2,p3,p4,resabs,
     *  resasc,result,res12,res24,u,x
      integer i,isym,k,kp,krul,neval
c
      dimension x(11),fval(25),cheb12(13),cheb24(25)
c
      external f,dqwgtc
c
c           the vector x contains the values cos(k*pi/24),
c           k = 1, ..., 11, to be used for the chebyshev series
c           expansion of f
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           ----------------------
c           fval   - value of the function f at the points
c                    cos(k*pi/24),  k = 0, ..., 24
c           cheb12 - chebyshev series expansion coefficients,
c                    for the function f, of degree 12
c           cheb24 - chebyshev series expansion coefficients,
c                    for the function f, of degree 24
c           res12  - approximation to the integral corresponding
c                    to the use of cheb12
c           res24  - approximation to the integral corresponding
c                    to the use of cheb24
c           dqwgtc - external function subprogram defining
c                    the weight function
c           hlgth  - half-length of the interval
c           centr  - mid point of the interval
c
c
c           check the position of c.
c
c***first executable statement  dqc25c
      cc = (0.2d+01*c-b-a)/(b-a)
      if(dabs(cc).lt.0.11d+01) go to 10
c
c           apply the 15-point gauss-kronrod scheme.
c
      krul = krul-1
      call dqk15w(f,dqwgtc,c,p2,p3,p4,kp,a,b,result,abserr,
     *  resabs,resasc)
      neval = 15
      if (resasc.eq.abserr) krul = krul+1
      go to 50
c
c           use the generalized clenshaw-curtis method.
c
   10 hlgth = 0.5d+00*(b-a)
      centr = 0.5d+00*(b+a)
      neval = 25
      fval(1) = 0.5d+00*f(hlgth+centr)
      fval(13) = f(centr)
      fval(25) = 0.5d+00*f(centr-hlgth)
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)
        fval(isym) = f(centr-u)
   20 continue
c
c           compute the chebyshev series expansion.
c
      call dqcheb(x,fval,cheb12,cheb24)
c
c           the modified chebyshev moments are computed by forward
c           recursion, using amom0 and amom1 as starting values.
c
      amom0 = dlog(dabs((0.1d+01-cc)/(0.1d+01+cc)))
      amom1 = 0.2d+01+cc*amom0
      res12 = cheb12(1)*amom0+cheb12(2)*amom1
      res24 = cheb24(1)*amom0+cheb24(2)*amom1
      do 30 k=3,13
        amom2 = 0.2d+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4d+01/(ak22-0.1d+01)
        res12 = res12+cheb12(k)*amom2
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   30 continue
      do 40 k=14,25
        amom2 = 0.2d+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4d+01/(ak22-0.1d+01)
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   40 continue
      result = res24
      abserr = dabs(res24-res12)
   50 return
      end
      subroutine dqc25f(f,a,b,omega,integr,nrmom,maxp1,ksave,result,
     *   abserr,neval,resabs,resasc,momcom,chebmo)
c***begin prologue  dqc25f
c***date written   810101   (yymmdd)
c***revision date  211011   (yymmdd)
c***category no.  h2a2a2
c***keywords  integration rules for functions with cos or sin
c             factor, clenshaw-curtis, gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute the integral i=integral of f(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
c            compute j = integral of abs(f) over (a,b). for small value
c            of omega or small intervals (a,b) the 15-point gauss-kronro
c            rule is used. otherwise a generalized clenshaw-curtis
c            method is used.
c***description
c
c        integration rules for functions with cos or sin factor
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to
c                    be declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c           omega  - double precision
c                    parameter in the weight function
c
c           integr - integer
c                    indicates which weight function is to be used
c                       integr = 1   w(x) = cos(omega*x)
c                       integr = 2   w(x) = sin(omega*x)
c
c           nrmom  - integer
c                    the length of interval (a,b) is equal to the length
c                    of the original integration interval divided by
c                    2**nrmom (we suppose that the routine is used in an
c                    adaptive integration process, otherwise set
c                    nrmom = 0). nrmom must be zero at the first call.
c
c           maxp1  - integer
c                    gives an upper bound on the number of chebyshev
c                    moments which can be stored, i.e. for the
c                    intervals of lengths abs(bb-aa)*2**(-l),
c                    l = 0,1,2, ..., maxp1-2.
c
c           ksave  - integer
c                    key which is one when the moments for the
c                    current interval have been computed
c
c         on return
c           result - double precision
c                    approximation to the integral i
c
c           abserr - double precision
c                    estimate of the modulus of the absolute
c                    error, which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of abs(f-i/(b-a))
c
c         on entry and return
c           momcom - integer
c                    for each interval length we need to compute the
c                    chebyshev moments. momcom counts the number of
c                    intervals for which these moments have already been
c                    computed. if nrmom.lt.momcom or ksave = 1, the
c                    chebyshev moments for the interval (a,b) have
c                    already been computed and stored, otherwise we
c                    compute them and we increase momcom.
c
c           chebmo - double precision
c                    array of dimension at least (maxp1,25) containing
c                    the modified chebyshev moments for the first momcom
c                    momcom interval lengths
c
c ......................................................................
c***references  (none)
c***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf
c***end prologue  dqc25f
c
      double precision a,abserr,ac,an,an2,as,asap,ass,b,centr,chebmo,
     *  cheb12,cheb24,conc,cons,cospar,d,dabs,dcos,dsin,dqwgtf,d1,
     *  d1mach,d2,estc,ests,f,fval,hlgth,oflow,omega,parint,par2,par22,
     *  p2,p3,p4,resabs,resasc,resc12,resc24,ress12,ress24,result,
     *  sinpar,v,x
      integer i,iers,integr,isym,j,k,ksave,m,momcom,neval,maxp1,
     *  noequ,noeq1,nrmom
c
      dimension chebmo(maxp1,25),cheb12(13),cheb24(25),d(25),d1(25),
     *  d2(25),fval(25),v(28),x(11)
c
      external f,dqwgtf
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ...,11, to be used for the chebyshev expansion of f
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fval   - value of the function f at the points
c                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (a,b)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (a,b)
c           resc12 - approximation to the integral of
c                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
c                    over (-1,+1), using the chebyshev series
c                    expansion of degree 12
c           resc24 - approximation to the same integral, using the
c                    chebyshev series expansion of degree 24
c           ress12 - the analogue of resc12 for the sine
c           ress24 - the analogue of resc24 for the sine
c
c
c           machine dependent constant
c           --------------------------
c
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqc25f
      oflow = d1mach(2)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      parint = omega*hlgth
c
c           compute the integral using the 15-point gauss-kronrod
c           formula if the value of the parameter in the integrand
c           is small.
c
      if(dabs(parint).gt.0.2d+01) go to 10
      call dqk15w(f,dqwgtf,omega,p2,p3,p4,integr,a,b,result,
     *  abserr,resabs,resasc)
      neval = 15
      go to 170
c
c           compute the integral using the generalized clenshaw-
c           curtis method.
c
   10 conc = hlgth*dcos(centr*omega)
      cons = hlgth*dsin(centr*omega)
      resasc = oflow
      neval = 25
c
c           check whether the chebyshev moments for this interval
c           have already been computed.
c
      if(nrmom.lt.momcom.or.ksave.eq.1) go to 120
c
c           compute a new set of chebyshev moments.
c
      m = momcom+1
      par2 = parint*parint
      par22 = par2+0.2d+01
      sinpar = dsin(parint)
      cospar = dcos(parint)
c
c           compute the chebyshev moments with respect to cosine.
c
      v(1) = 0.2d+01*sinpar/parint
      v(2) = (0.8d+01*cospar+(par2+par2-0.8d+01)*sinpar/parint)/par2
      v(3) = (0.32d+02*(par2-0.12d+02)*cospar+(0.2d+01*
     *  ((par2-0.80d+02)*par2+0.192d+03)*sinpar)/parint)/(par2*par2)
      ac = 0.8d+01*cospar
      as = 0.24d+02*parint*sinpar
      if(dabs(parint).gt.0.24d+02) go to 30
c
c           compute the chebyshev moments as the solutions of a
c           boundary value problem with 1 initial value (v(3)) and 1
c           end value (computed using an asymptotic formula).
c
      noequ = 25
      noeq1 = noequ-1
      an = 0.6d+01
      do 20 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+3) = as-(an2-0.4d+01)*ac
        an = an+0.2d+01
   20 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+3) = as-(an2-0.4d+01)*ac
      v(4) = v(4)-0.56d+02*par2*v(3)
      ass = parint*sinpar
      asap = (((((0.210d+03*par2-0.1d+01)*cospar-(0.105d+03*par2
     *  -0.63d+02)*ass)/an2-(0.1d+01-0.15d+02*par2)*cospar
     *  +0.15d+02*ass)/an2-cospar+0.3d+01*ass)/an2-cospar)/an2
      v(noequ+3) = v(noequ+3)-0.2d+01*asap*par2*(an-0.1d+01)*
     *   (an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl must be replaced by call to
c***        double precision version of linpack routine sgtsl
c
      call dgtsl(noequ,d1,d,d2,v(4),iers)
      go to 50
c
c           compute the chebyshev moments by means of forward
c           recursion.
c
   30 an = 0.4d+01
      do 40 i = 4,13
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)-ac)
     *  +as-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))/
     *  (par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   40 continue
   50 do 60 j = 1,13
        chebmo(m,2*j-1) = v(j)
   60 continue
c
c           compute the chebyshev moments with respect to sine.
c
      v(1) = 0.2d+01*(sinpar-parint*cospar)/par2
      v(2) = (0.18d+02-0.48d+02/par2)*sinpar/par2
     *  +(-0.2d+01+0.48d+02/par2)*cospar/parint
      ac = -0.24d+02*parint*cospar
      as = -0.8d+01*sinpar
      if(dabs(parint).gt.0.24d+02) go to 80
c
c           compute the chebyshev moments as the solutions of a boundary
c           value problem with 1 initial value (v(2)) and 1 end value
c           (computed using an asymptotic formula).
c
      an = 0.5d+01
      do 70 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
        d2(k) = (an-0.1d+01)*(an-0.2d+01)*par2
        d1(k+1) = (an+0.3d+01)*(an+0.4d+01)*par2
        v(k+2) = ac+(an2-0.4d+01)*as
        an = an+0.2d+01
   70 continue
      an2 = an*an
      d(noequ) = -0.2d+01*(an2-0.4d+01)*(par22-an2-an2)
      v(noequ+2) = ac+(an2-0.4d+01)*as
      v(3) = v(3)-0.42d+02*par2*v(2)
      ass = parint*cospar
      asap = (((((0.105d+03*par2-0.63d+02)*ass+(0.210d+03*par2
     *  -0.1d+01)*sinpar)/an2+(0.15d+02*par2-0.1d+01)*sinpar-
     *  0.15d+02*ass)/an2-0.3d+01*ass-sinpar)/an2-sinpar)/an2
      v(noequ+2) = v(noequ+2)-0.2d+01*asap*par2*(an-0.1d+01)
     *  *(an-0.2d+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
c***        call to dgtsl must be replaced by call to
c***        double precision version of linpack routine sgtsl
c
      call dgtsl(noequ,d1,d,d2,v(3),iers)
      go to 100
c
c           compute the chebyshev moments by means of forward recursion.
c
   80 an = 0.3d+01
      do 90 i = 3,12
        an2 = an*an
        v(i) = ((an2-0.4d+01)*(0.2d+01*(par22-an2-an2)*v(i-1)+as)
     *  +ac-par2*(an+0.1d+01)*(an+0.2d+01)*v(i-2))
     *  /(par2*(an-0.1d+01)*(an-0.2d+01))
        an = an+0.2d+01
   90 continue
  100 do 110 j = 1,12
        chebmo(m,2*j) = v(j)
  110 continue
  120 if (nrmom.lt.momcom) m = nrmom+1
       if (momcom.lt.(maxp1-1).and.nrmom.ge.momcom) momcom = momcom+1
c
c           compute the coefficients of the chebyshev expansions
c           of degrees 12 and 24 of the function f.
c
      fval(1) = 0.5d+00*f(centr+hlgth)
      fval(13) = f(centr)
      fval(25) = 0.5d+00*f(centr-hlgth)
      do 130 i = 2,12
        isym = 26-i
        fval(i) = f(hlgth*x(i-1)+centr)
        fval(isym) = f(centr-hlgth*x(i-1))
  130 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           compute the integral and error estimates.
c
      resc12 = cheb12(13)*chebmo(m,13)
      ress12 = 0.0d+00
      k = 11
      do 140 j = 1,6
        resc12 = resc12+cheb12(k)*chebmo(m,k)
        ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
        k = k-2
  140 continue
      resc24 = cheb24(25)*chebmo(m,25)
      ress24 = 0.0d+00
      resabs = dabs(cheb24(25))
      k = 23
      do 150 j = 1,12
        resc24 = resc24+cheb24(k)*chebmo(m,k)
        ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
        resabs = resabs+dabs(cheb24(k))+dabs(cheb24(k+1))
        k = k-2
  150 continue
      estc = dabs(resc24-resc12)
      ests = dabs(ress24-ress12)
      resabs = resabs*dabs(hlgth)
      if(integr.eq.2) go to 160
      result = conc*resc24-cons*ress24
      abserr = dabs(conc*estc)+dabs(cons*ests)
      go to 170
  160 result = conc*ress24+cons*resc24
      abserr = dabs(conc*ests)+dabs(cons*estc)
  170 return
      end
      subroutine dqc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result,
     *   abserr,resasc,integr,nev)
c***begin prologue  dqc25s
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (bl,br), with error
c            estimate, where the weight function w has a singular
c            behaviour of algebraico-logarithmic type at the points
c            a and/or b. (bl,br) is a part of (a,b).
c***description
c
c        integration rules for integrands having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           f      - double precision
c                    function subprogram defining the integrand
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - double precision
c                    left end point of the original interval
c
c           b      - double precision
c                    right end point of the original interval, b.gt.a
c
c           bl     - double precision
c                    lower limit of integration, bl.ge.a
c
c           br     - double precision
c                    upper limit of integration, br.le.b
c
c           alfa   - double precision
c                    parameter in the weight function
c
c           beta   - double precision
c                    parameter in the weight function
c
c           ri,rj,rg,rh - double precision
c                    modified chebyshev moments for the application
c                    of the generalized clenshaw-curtis
c                    method (computed in subroutine dqmomo)
c
c           result - double precision
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if b1 = a or br = b.
c                    in all other cases the 15-point kronrod
c                    rule is applied, obtained by optimal addition of
c                    abscissae to the 7-point gauss rule.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           resasc - double precision
c                    approximation to the integral of abs(f*w-i/(b-a))
c
c           integr - integer
c                    which determines the weight function
c                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
c                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
c                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
c                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
c                                 log(b-x)
c
c           nev    - integer
c                    number of integrand evaluations
c***references  (none)
c***routines called  dqcheb,dqk15w
c***end prologue  dqc25s
c
      double precision a,abserr,alfa,b,beta,bl,br,centr,cheb12,cheb24,
     *  dabs,dc,dlog,f,factor,fix,fval,hlgth,resabs,resasc,result,res12,
     *  res24,rg,rh,ri,rj,u,dqwgts,x
      integer i,integr,isym,nev
c
      dimension cheb12(13),cheb24(25),fval(25),rg(25),rh(25),ri(25),
     *  rj(25),x(11)
c
      external f,dqwgts
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ..., 11, to be used for the computation of the
c           chebyshev series expansion of f.
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           -----------------------
c
c           fval   - value of the function f at the points
c                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
c                    k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (bl,br)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (bl,br)
c           res12  - approximation to the integral obtained from cheb12
c           res24  - approximation to the integral obtained from cheb24
c           dqwgts - external function subprogram defining
c                    the four possible weight functions
c           hlgth  - half-length of the interval (bl,br)
c           centr  - mid point of the interval (bl,br)
c
c***first executable statement  dqc25s
      nev = 25
      if(bl.eq.a.and.(alfa.ne.0.0d+00.or.integr.eq.2.or.integr.eq.4))
     * go to 10
      if(br.eq.b.and.(beta.ne.0.0d+00.or.integr.eq.3.or.integr.eq.4))
     * go to 140
c
c           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod
c           scheme.
c
c
      call dqk15w(f,dqwgts,a,b,alfa,beta,integr,bl,br,
     *    result,abserr,resabs,resasc)
      nev = 15
      go to 270
c
c           this part of the program is executed only if a = bl.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
c                  *f(0.5*(br-a)*x+0.5*(br+a))
c
   10 hlgth = 0.5d+00*(br-bl)
      centr = 0.5d+00*(br+bl)
      fix = b-centr
      fval(1) = 0.5d+00*f(hlgth+centr)*(fix-hlgth)**beta
      fval(13) = f(centr)*(fix**beta)
      fval(25) = 0.5d+00*f(centr-hlgth)*(fix+hlgth)**beta
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix-u)**beta
        fval(isym) = f(centr-u)*(fix+u)**beta
   20 continue
      factor = hlgth**(alfa+0.1d+01)
      result = 0.0d+00
      abserr = 0.0d+00
      res12 = 0.0d+00
      res24 = 0.0d+00
      if(integr.gt.2) go to 70
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 1  (or 2)
c
      do 30 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   30 continue
      do 40 i=14,25
        res24 = res24+cheb24(i)*ri(i)
   40 continue
      if(integr.eq.1) go to 130
c
c           integr = 2
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 50 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res12+cheb24(i)*rg(i)
   50 continue
      do 60 i=14,25
        res24 = res24+cheb24(i)*rg(i)
   60 continue
      go to 130
c
c           compute the chebyshev series expansion of the
c           following function
c           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
c
   70 fval(1) = fval(1)*dlog(fix-hlgth)
      fval(13) = fval(13)*dlog(fix)
      fval(25) = fval(25)*dlog(fix+hlgth)
      do 80 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*dlog(fix-u)
        fval(isym) = fval(isym)*dlog(fix+u)
   80 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 3  (or 4)
c
      do 90 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   90 continue
      do 100 i=14,25
        res24 = res24+cheb24(i)*ri(i)
  100 continue
      if(integr.eq.3) go to 130
c
c           integr = 4
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 110 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
  110 continue
      do 120 i=14,25
        res24 = res24+cheb24(i)*rg(i)
  120 continue
  130 result = (result+res24)*factor
      abserr = (abserr+dabs(res24-res12))*factor
      go to 270
c
c           this part of the program is executed only if b = br.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
c                *f(0.5*(b-bl)*x+0.5*(b+bl))
c
  140 hlgth = 0.5d+00*(br-bl)
      centr = 0.5d+00*(br+bl)
      fix = centr-a
      fval(1) = 0.5d+00*f(hlgth+centr)*(fix+hlgth)**alfa
      fval(13) = f(centr)*(fix**alfa)
      fval(25) = 0.5d+00*f(centr-hlgth)*(fix-hlgth)**alfa
      do 150 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix+u)**alfa
        fval(isym) = f(centr-u)*(fix-u)**alfa
  150 continue
      factor = hlgth**(beta+0.1d+01)
      result = 0.0d+00
      abserr = 0.0d+00
      res12 = 0.0d+00
      res24 = 0.0d+00
      if(integr.eq.2.or.integr.eq.4) go to 200
c
c           integr = 1  (or 3)
c
      call dqcheb(x,fval,cheb12,cheb24)
      do 160 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  160 continue
      do 170 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  170 continue
      if(integr.eq.1) go to 260
c
c           integr = 3
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 180 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  180 continue
      do 190 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  190 continue
      go to 260
c
c           compute the chebyshev series expansion of the
c           following function
c           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
c
  200 fval(1) = fval(1)*dlog(hlgth+fix)
      fval(13) = fval(13)*dlog(fix)
      fval(25) = fval(25)*dlog(fix-hlgth)
      do 210 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*dlog(u+fix)
        fval(isym) = fval(isym)*dlog(fix-u)
  210 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 2  (or 4)
c
      do 220 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  220 continue
      do 230 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  230 continue
      if(integr.eq.2) go to 260
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
c
c           integr = 4
c
      do 240 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  240 continue
      do 250 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  250 continue
  260 result = (result+res24)*factor
      abserr = (abserr+dabs(res24-res12))*factor
  270 return
      end
      subroutine dqcheb(x,fval,cheb12,cheb24)
c***begin prologue  dqcheb
c***refer to  dqc25c,dqc25f,dqc25s
c***routines called  (none)
c***revision date  830518   (yymmdd)
c***keywords  chebyshev series expansion, fast fourier transform
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes the chebyshev series expansion
c            of degrees 12 and 24 of a function using a
c            fast fourier transform method
c            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
c            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
c            where t(k,x) is the chebyshev polynomial of degree k.
c***description
c
c        chebyshev series expansion
c        standard fortran subroutine
c        double precision version
c
c        parameters
c          on entry
c           x      - double precision
c                    vector of dimension 11 containing the
c                    values cos(k*pi/24), k = 1, ..., 11
c
c           fval   - double precision
c                    vector of dimension 25 containing the
c                    function values at the points
c                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
c                    where (a,b) is the approximation interval.
c                    fval(1) and fval(25) are divided by two
c                    (these values are destroyed at output).
c
c          on return
c           cheb12 - double precision
c                    vector of dimension 13 containing the
c                    chebyshev coefficients for degree 12
c
c           cheb24 - double precision
c                    vector of dimension 25 containing the
c                    chebyshev coefficients for degree 24
c
c***end prologue  dqcheb
c
      double precision alam,alam1,alam2,cheb12,cheb24,fval,part1,part2,
     *  part3,v,x
      integer i,j
c
      dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)
c
c***first executable statement  dqcheb
      do 10 i=1,12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   10 continue
      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)
     *  +x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)
     *  +x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)
     *  -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)
     *  -x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam
      do 20 i=1,6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   20 continue
      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam
      do 30 i=1,3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   30 continue
      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 0.1d+01/0.6d+01
      do 40 i=2,12
        cheb12(i) = cheb12(i)*alam
   40 continue
      alam = 0.5d+00*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam
      do 50 i=2,24
        cheb24(i) = cheb24(i)*alam
   50 continue
      cheb24(1) = 0.5d+00*alam*cheb24(1)
      cheb24(25) = 0.5d+00*alam*cheb24(25)
      return
      end
      subroutine dqelg(n,epstab,result,abserr,res3la,nres)
c***begin prologue  dqelg
c***refer to  dqagie,dqagoe,dqagpe,dqagse
c***routines called  d1mach
c***revision date  830518   (yymmdd)
c***keywords  epsilon algorithm, convergence acceleration,
c             extrapolation
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine determines the limit of a given sequence of
c            approximations, by means of the epsilon algorithm of
c            p.wynn. an estimate of the absolute error is also given.
c            the condensed epsilon table is computed. only those
c            elements needed for the computation of the next diagonal
c            are preserved.
c***description
c
c           epsilon algorithm
c           standard fortran subroutine
c           double precision version
c
c           parameters
c              n      - integer
c                       epstab(n) contains the new element in the
c                       first column of the epsilon table.
c
c              epstab - double precision
c                       vector of dimension 52 containing the elements
c                       of the two lower diagonals of the triangular
c                       epsilon table. the elements are numbered
c                       starting at the right-hand corner of the
c                       triangle.
c
c              result - double precision
c                       resulting approximation to the integral
c
c              abserr - double precision
c                       estimate of the absolute error computed from
c                       result and the 3 previous results
c
c              res3la - double precision
c                       vector of dimension 3 containing the last 3
c                       results
c
c              nres   - integer
c                       number of calls to the routine
c                       (should be zero at first call)
c
c***end prologue  dqelg
c
      double precision abserr,dabs,delta1,delta2,delta3,dmax1,d1mach,
     *  epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
     *  oflow,res,result,res3la,ss,tol1,tol2,tol3
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
c
c           list of major variables
c           -----------------------
c
c           e0     - the 4 elements on which the computation of a new
c           e1       element in the epsilon table is based
c           e2
c           e3                 e0
c                        e3    e1    new
c                              e2
c           newelm - number of elements to be computed in the new
c                    diagonal
c           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
c           result - the element in the new diagonal with least value
c                    of error
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           oflow is the largest positive magnitude.
c           limexp is the maximum number of elements the epsilon
c           table can contain. if this number is reached, the upper
c           diagonal of the epsilon table is deleted.
c
c***first executable statement  dqelg
      epmach = d1mach(4)
      oflow = d1mach(2)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = dabs(e1)
        delta2 = e2-e1
        err2 = dabs(delta2)
        tol2 = dmax1(dabs(e2),e1abs)*epmach
        delta3 = e1-e0
        err3 = dabs(delta3)
        tol3 = dmax1(e1abs,dabs(e0))*epmach
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10
c
c           if e0, e1 and e2 are equal to within machine
c           accuracy, convergence is assumed.
c           result = e2
c           abserr = abs(e1-e0)+abs(e2-e1)
c
        result = res
        abserr = err2+err3
c ***jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = dabs(delta1)
        tol1 = dmax1(e1abs,dabs(e3))*epmach
c
c           if two elements are very close to each other, omit
c           a part of the table by adjusting the value of n
c
        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 0.1d+01/delta1+0.1d+01/delta2-0.1d+01/delta3
        epsinf = dabs(ss*e1)
c
c           test to detect irregular behaviour in the table, and
c           eventually omit a part of the table adjusting the value
c           of n.
c
        if(epsinf.gt.0.1d-03) go to 30
   20   n = i+i-1
c ***jump out of do-loop
        go to 50
c
c           compute a new element and eventually adjust
c           the value of result.
c
   30   res = e1+0.1d+01/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+dabs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        result = res
   40 continue
c
c           shift the table.
c
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = result
      abserr = oflow
      go to 100
c
c           compute error estimate
c
   90 abserr = dabs(result-res3la(3))+dabs(result-res3la(2))
     *  +dabs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = dmax1(abserr,0.5d+01*epmach*dabs(result))
      return
      end
      subroutine dqk15(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk15
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the7-point gauss rule(resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.1294849661 6886969327 0611432679 082 d0 /
      data wg  (  2) / 0.2797053914 8927666790 1467771423 780 d0 /
      data wg  (  3) / 0.3818300505 0511894495 0369775488 975 d0 /
      data wg  (  4) / 0.4179591836 7346938775 5102040816 327 d0 /
c
      data xgk (  1) / 0.9914553711 2081263920 6854697526 329 d0 /
      data xgk (  2) / 0.9491079123 4275852452 6189684047 851 d0 /
      data xgk (  3) / 0.8648644233 5976907278 9712788640 926 d0 /
      data xgk (  4) / 0.7415311855 9939443986 3864773280 788 d0 /
      data xgk (  5) / 0.5860872354 6769113029 4144838258 730 d0 /
      data xgk (  6) / 0.4058451513 7739716690 6606412076 961 d0 /
      data xgk (  7) / 0.2077849550 0789846760 0689403773 245 d0 /
      data xgk (  8) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0229353220 1052922496 3732008058 970 d0 /
      data wgk (  2) / 0.0630920926 2997855329 0700663189 204 d0 /
      data wgk (  3) / 0.1047900103 2225018383 9876322541 518 d0 /
      data wgk (  4) / 0.1406532597 1552591874 5189590510 238 d0 /
      data wgk (  5) / 0.1690047266 3926790282 6583426598 550 d0 /
      data wgk (  6) / 0.1903505780 6478540991 3256402421 014 d0 /
      data wgk (  7) / 0.2044329400 7529889241 4161999234 649 d0 /
      data wgk (  8) / 0.2094821410 8472782801 2999174891 714 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk15i
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a2,h2a4a2
c***keywords  15-point transformed gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the original (infinite integration range is mapped
c            onto the interval (0,1) and (a,b) is a part of (0,1).
c            it is the purpose to compute
c            i = integral of transformed integrand over (a,b),
c            j = integral of abs(transformed integrand) over (a,b).
c***description
c
c           integration rule
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       fuction subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              boun   - double precision
c                       finite bound of original integration
c                       range (set to zero if inf = +2)
c
c              inf    - integer
c                       if inf = -1, the original interval is
c                                   (-infinity,bound),
c                       if inf = +1, the original interval is
c                                   (bound,+infinity),
c                       if inf = +2, the original interval is
c                                   (-infinity,+infinity) and
c                       the integral is computed as the sum of two
c                       integrals, one over (-infinity,0) and one over
c                       (0,+infinity).
c
c              a      - double precision
c                       lower limit for integration over subrange
c                       of (0,1)
c
c              b      - double precision
c                       upper limit for integration over subrange
c                       of (0,1)
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule(resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule(resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of
c                       abs((transformed integrand)-i/(b-a)) over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15i
c
      double precision a,absc,absc1,absc2,abserr,b,boun,centr,dabs,dinf,
     *  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk,
     *  xgk
      integer inf,j
      external f
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
c
c           the abscissae and weights are supplied for the interval
c           (-1,1).  because of symmetry only the positive abscissae and
c           their corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule, corresponding
c                    to the abscissae xgk(2), xgk(4), ...
c                    wg(1), wg(3), ... are set to zero.
c
      data wg(1) / 0.0d0 /
      data wg(2) / 0.1294849661 6886969327 0611432679 082d0 /
      data wg(3) / 0.0d0 /
      data wg(4) / 0.2797053914 8927666790 1467771423 780d0 /
      data wg(5) / 0.0d0 /
      data wg(6) / 0.3818300505 0511894495 0369775488 975d0 /
      data wg(7) / 0.0d0 /
      data wg(8) / 0.4179591836 7346938775 5102040816 327d0 /
c
      data xgk(1) / 0.9914553711 2081263920 6854697526 329d0 /
      data xgk(2) / 0.9491079123 4275852452 6189684047 851d0 /
      data xgk(3) / 0.8648644233 5976907278 9712788640 926d0 /
      data xgk(4) / 0.7415311855 9939443986 3864773280 788d0 /
      data xgk(5) / 0.5860872354 6769113029 4144838258 730d0 /
      data xgk(6) / 0.4058451513 7739716690 6606412076 961d0 /
      data xgk(7) / 0.2077849550 0789846760 0689403773 245d0 /
      data xgk(8) / 0.0000000000 0000000000 0000000000 000d0 /
c
      data wgk(1) / 0.0229353220 1052922496 3732008058 970d0 /
      data wgk(2) / 0.0630920926 2997855329 0700663189 204d0 /
      data wgk(3) / 0.1047900103 2225018383 9876322541 518d0 /
      data wgk(4) / 0.1406532597 1552591874 5189590510 238d0 /
      data wgk(5) / 0.1690047266 3926790282 6583426598 550d0 /
      data wgk(6) / 0.1903505780 6478540991 3256402421 014d0 /
      data wgk(7) / 0.2044329400 7529889241 4161999234 649d0 /
      data wgk(8) / 0.2094821410 8472782801 2999174891 714d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           tabsc* - transformed abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of the transformed
c                    integrand over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15i
      epmach = d1mach(4)
      uflow = d1mach(1)
      dinf = min0(1,inf)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      tabsc1 = boun+dinf*(0.1d+01-centr)/centr
      fval1 = f(tabsc1)
      if(inf.eq.2) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the error.
c
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
        tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
        fval1 = f(tabsc1)
        fval2 = f(tabsc2)
        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
        if(inf.eq.2) fval2 = fval2+f(-tabsc2)
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(dabs(fval1)+dabs(fval2))
   10 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d0) abserr = resasc*
     * dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     * ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr,
     *   resabs,resasc)
c***begin prologue  dqk15w
c***date written   810101   (yymmdd)
c***revision date  830518   (mmddyy)
c***category no.  h2a2a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b), with error
c                           estimate
c                       j = integral of abs(f*w) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c             on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              w      - double precision
c                       function subprogram defining the integrand
c                       weight function w(x). the actual name for w
c                       needs to be declared e x t e r n a l in the
c                       calling program.
c
c              p1, p2, p3, p4 - double precision
c                       parameters in the weight function
c
c              kp     - integer
c                       key for indicating the type of weight function
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral of abs(f)
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15w
c
      double precision a,absc,absc1,absc2,abserr,b,centr,dabs,dhlgth,
     *  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,w,wg,wgk,
     *  xgk
      integer j,jtw,jtwm1,kp
      external f,w
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point gauss-kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ... abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point gauss-kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/
     *     0.9914553711208126d+00,     0.9491079123427585d+00,
     *     0.8648644233597691d+00,     0.7415311855993944d+00,
     *     0.5860872354676911d+00,     0.4058451513773972d+00,
     *     0.2077849550078985d+00,     0.0000000000000000d+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/
     *     0.2293532201052922d-01,     0.6309209262997855d-01,
     *     0.1047900103222502d+00,     0.1406532597155259d+00,
     *     0.1690047266392679d+00,     0.1903505780647854d+00,
     *     0.2044329400752989d+00,     0.2094821410847278d+00/
c
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697d+00,    0.2797053914892767d+00,
     *     0.3818300505051189d+00,    0.4179591836734694d+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f*w over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15w
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 15-point kronrod approximation to the
c           integral, and estimate the error.
c
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*
     *  0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk21(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data wg  (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data wg  (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data wg  (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data wg  (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data xgk (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data xgk (  2) / 0.9739065285 1717172007 7964012084 452 d0 /
      data xgk (  3) / 0.9301574913 5570822600 1207180059 508 d0 /
      data xgk (  4) / 0.8650633666 8898451073 2096688423 493 d0 /
      data xgk (  5) / 0.7808177265 8641689706 3717578345 042 d0 /
      data xgk (  6) / 0.6794095682 9902440623 4327365114 874 d0 /
      data xgk (  7) / 0.5627571346 6860468333 9000099272 694 d0 /
      data xgk (  8) / 0.4333953941 2924719079 9265943165 784 d0 /
      data xgk (  9) / 0.2943928627 0146019813 1126603103 866 d0 /
      data xgk ( 10) / 0.1488743389 8163121088 4826001129 720 d0 /
      data xgk ( 11) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data wgk (  2) / 0.0325581623 0796472747 8818972459 390 d0 /
      data wgk (  3) / 0.0547558965 7435199603 1381300244 580 d0 /
      data wgk (  4) / 0.0750396748 1091995276 7043140916 190 d0 /
      data wgk (  5) / 0.0931254545 8369760553 5065465083 366 d0 /
      data wgk (  6) / 0.1093871588 0229764189 9210590325 805 d0 /
      data wgk (  7) / 0.1234919762 6206585107 7958109831 074 d0 /
      data wgk (  8) / 0.1347092173 1147332592 8054001771 707 d0 /
      data wgk (  9) / 0.1427759385 7706008079 7094273138 717 d0 /
      data wgk ( 10) / 0.1477391049 0133849137 4841515972 068 d0 /
      data wgk ( 11) / 0.1494455540 0291690566 4936468389 821 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk31(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk31
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  31-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 31-point
c                       gauss-kronrod rule (resk), obtained by optimal
c                       addition of abscissae to the 15-point gauss
c                       rule (resg).
c
c              abserr - double precison
c                       estimate of the modulus of the modulus,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk31
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 31-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 15-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 15-point gauss rule
c
c           wgk    - weights of the 31-point kronrod rule
c
c           wg     - weights of the 15-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0307532419 9611726835 4628393577 204 d0 /
      data wg  (  2) / 0.0703660474 8810812470 9267416450 667 d0 /
      data wg  (  3) / 0.1071592204 6717193501 1869546685 869 d0 /
      data wg  (  4) / 0.1395706779 2615431444 7804794511 028 d0 /
      data wg  (  5) / 0.1662692058 1699393355 3200860481 209 d0 /
      data wg  (  6) / 0.1861610000 1556221102 6800561866 423 d0 /
      data wg  (  7) / 0.1984314853 2711157645 6118326443 839 d0 /
      data wg  (  8) / 0.2025782419 2556127288 0620199967 519 d0 /
c
      data xgk (  1) / 0.9980022986 9339706028 5172840152 271 d0 /
      data xgk (  2) / 0.9879925180 2048542848 9565718586 613 d0 /
      data xgk (  3) / 0.9677390756 7913913425 7347978784 337 d0 /
      data xgk (  4) / 0.9372733924 0070590430 7758947710 209 d0 /
      data xgk (  5) / 0.8972645323 4408190088 2509656454 496 d0 /
      data xgk (  6) / 0.8482065834 1042721620 0648320774 217 d0 /
      data xgk (  7) / 0.7904185014 4246593296 7649294817 947 d0 /
      data xgk (  8) / 0.7244177313 6017004741 6186054613 938 d0 /
      data xgk (  9) / 0.6509967412 9741697053 3735895313 275 d0 /
      data xgk ( 10) / 0.5709721726 0853884753 7226737253 911 d0 /
      data xgk ( 11) / 0.4850818636 4023968069 3655740232 351 d0 /
      data xgk ( 12) / 0.3941513470 7756336989 7207370981 045 d0 /
      data xgk ( 13) / 0.2991800071 5316881216 6780024266 389 d0 /
      data xgk ( 14) / 0.2011940939 9743452230 0628303394 596 d0 /
      data xgk ( 15) / 0.1011420669 1871749902 7074231447 392 d0 /
      data xgk ( 16) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0053774798 7292334898 7792051430 128 d0 /
      data wgk (  2) / 0.0150079473 2931612253 8374763075 807 d0 /
      data wgk (  3) / 0.0254608473 2671532018 6874001019 653 d0 /
      data wgk (  4) / 0.0353463607 9137584622 2037948478 360 d0 /
      data wgk (  5) / 0.0445897513 2476487660 8227299373 280 d0 /
      data wgk (  6) / 0.0534815246 9092808726 5343147239 430 d0 /
      data wgk (  7) / 0.0620095678 0067064028 5139230960 803 d0 /
      data wgk (  8) / 0.0698541213 1872825870 9520077099 147 d0 /
      data wgk (  9) / 0.0768496807 5772037889 4432777482 659 d0 /
      data wgk ( 10) / 0.0830805028 2313302103 8289247286 104 d0 /
      data wgk ( 11) / 0.0885644430 5621177064 7275443693 774 d0 /
      data wgk ( 12) / 0.0931265981 7082532122 5486872747 346 d0 /
      data wgk ( 13) / 0.0966427269 8362367850 5179907627 589 d0 /
      data wgk ( 14) / 0.0991735987 2179195933 2393173484 603 d0 /
      data wgk ( 15) / 0.1007698455 2387559504 4946662617 570 d0 /
      data wgk ( 16) / 0.1013300070 1479154901 7374792767 493 d0 /
c
c
c           list of major variables
c           -----------------------
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 15-point gauss formula
c           resk   - result of the 31-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c***first executable statement  dqk31
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 31-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(16)*dabs(fc-reskh)
      do 20 j=1,15
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk41(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk41
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  41-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 41-point
c                       gauss-kronrod rule (resk) obtained by optimal
c                       addition of abscissae to the 20-point gauss
c                       rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integal of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk41
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 41-point gauss-kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 20-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 20-point gauss rule
c
c           wgk    - weights of the 41-point gauss-kronrod rule
c
c           wg     - weights of the 20-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0176140071 3915211831 1861962351 853 d0 /
      data wg  (  2) / 0.0406014298 0038694133 1039952274 932 d0 /
      data wg  (  3) / 0.0626720483 3410906356 9506535187 042 d0 /
      data wg  (  4) / 0.0832767415 7670474872 4758143222 046 d0 /
      data wg  (  5) / 0.1019301198 1724043503 6750135480 350 d0 /
      data wg  (  6) / 0.1181945319 6151841731 2377377711 382 d0 /
      data wg  (  7) / 0.1316886384 4917662689 8494499748 163 d0 /
      data wg  (  8) / 0.1420961093 1838205132 9298325067 165 d0 /
      data wg  (  9) / 0.1491729864 7260374678 7828737001 969 d0 /
      data wg  ( 10) / 0.1527533871 3072585069 8084331955 098 d0 /
c
      data xgk (  1) / 0.9988590315 8827766383 8315576545 863 d0 /
      data xgk (  2) / 0.9931285991 8509492478 6122388471 320 d0 /
      data xgk (  3) / 0.9815078774 5025025919 3342994720 217 d0 /
      data xgk (  4) / 0.9639719272 7791379126 7666131197 277 d0 /
      data xgk (  5) / 0.9408226338 3175475351 9982722212 443 d0 /
      data xgk (  6) / 0.9122344282 5132590586 7752441203 298 d0 /
      data xgk (  7) / 0.8782768112 5228197607 7442995113 078 d0 /
      data xgk (  8) / 0.8391169718 2221882339 4529061701 521 d0 /
      data xgk (  9) / 0.7950414288 3755119835 0638833272 788 d0 /
      data xgk ( 10) / 0.7463319064 6015079261 4305070355 642 d0 /
      data xgk ( 11) / 0.6932376563 3475138480 5490711845 932 d0 /
      data xgk ( 12) / 0.6360536807 2651502545 2836696226 286 d0 /
      data xgk ( 13) / 0.5751404468 1971031534 2946036586 425 d0 /
      data xgk ( 14) / 0.5108670019 5082709800 4364050955 251 d0 /
      data xgk ( 15) / 0.4435931752 3872510319 9992213492 640 d0 /
      data xgk ( 16) / 0.3737060887 1541956067 2548177024 927 d0 /
      data xgk ( 17) / 0.3016278681 1491300432 0555356858 592 d0 /
      data xgk ( 18) / 0.2277858511 4164507808 0496195368 575 d0 /
      data xgk ( 19) / 0.1526054652 4092267550 5220241022 678 d0 /
      data xgk ( 20) / 0.0765265211 3349733375 4640409398 838 d0 /
      data xgk ( 21) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0030735837 1852053150 1218293246 031 d0 /
      data wgk (  2) / 0.0086002698 5564294219 8661787950 102 d0 /
      data wgk (  3) / 0.0146261692 5697125298 3787960308 868 d0 /
      data wgk (  4) / 0.0203883734 6126652359 8010231432 755 d0 /
      data wgk (  5) / 0.0258821336 0495115883 4505067096 153 d0 /
      data wgk (  6) / 0.0312873067 7703279895 8543119323 801 d0 /
      data wgk (  7) / 0.0366001697 5820079803 0557240707 211 d0 /
      data wgk (  8) / 0.0416688733 2797368626 3788305936 895 d0 /
      data wgk (  9) / 0.0464348218 6749767472 0231880926 108 d0 /
      data wgk ( 10) / 0.0509445739 2372869193 2707670050 345 d0 /
      data wgk ( 11) / 0.0551951053 4828599474 4832372419 777 d0 /
      data wgk ( 12) / 0.0591114008 8063957237 4967220648 594 d0 /
      data wgk ( 13) / 0.0626532375 5478116802 5870122174 255 d0 /
      data wgk ( 14) / 0.0658345971 3361842211 1563556969 398 d0 /
      data wgk ( 15) / 0.0686486729 2852161934 5623411885 368 d0 /
      data wgk ( 16) / 0.0710544235 5344406830 5790361723 210 d0 /
      data wgk ( 17) / 0.0730306903 3278666749 5189417658 913 d0 /
      data wgk ( 18) / 0.0745828754 0049918898 6581418362 488 d0 /
      data wgk ( 19) / 0.0757044976 8455667465 9542775376 617 d0 /
      data wgk ( 20) / 0.0763778676 7208073670 5502835038 061 d0 /
      data wgk ( 21) / 0.0766007119 1799965644 5049901530 102 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 20-point gauss formula
c           resk   - result of the 41-point kronrod formula
c           reskh  - approximation to mean value of f over (a,b), i.e.
c                    to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk41
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 41-point gauss-kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr)
      resk = wgk(21)*fc
      resabs = dabs(resk)
      do 10 j=1,10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(21)*dabs(fc-reskh)
      do 20 j=1,20
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk51(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk51
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  51-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
c              f      - double precision
c                       function subroutine defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 51-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 25-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk51
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 51-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 25-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 25-point gauss rule
c
c           wgk    - weights of the 51-point kronrod rule
c
c           wg     - weights of the 25-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0113937985 0102628794 7902964113 235 d0 /
      data wg  (  2) / 0.0263549866 1503213726 1901815295 299 d0 /
      data wg  (  3) / 0.0409391567 0130631265 5623487711 646 d0 /
      data wg  (  4) / 0.0549046959 7583519192 5936891540 473 d0 /
      data wg  (  5) / 0.0680383338 1235691720 7187185656 708 d0 /
      data wg  (  6) / 0.0801407003 3500101801 3234959669 111 d0 /
      data wg  (  7) / 0.0910282619 8296364981 1497220702 892 d0 /
      data wg  (  8) / 0.1005359490 6705064420 2206890392 686 d0 /
      data wg  (  9) / 0.1085196244 7426365311 6093957050 117 d0 /
      data wg  ( 10) / 0.1148582591 4571164833 9325545869 556 d0 /
      data wg  ( 11) / 0.1194557635 3578477222 8178126512 901 d0 /
      data wg  ( 12) / 0.1222424429 9031004168 8959518945 852 d0 /
      data wg  ( 13) / 0.1231760537 2671545120 3902873079 050 d0 /
c
      data xgk (  1) / 0.9992621049 9260983419 3457486540 341 d0 /
      data xgk (  2) / 0.9955569697 9049809790 8784946893 902 d0 /
      data xgk (  3) / 0.9880357945 3407724763 7331014577 406 d0 /
      data xgk (  4) / 0.9766639214 5951751149 8315386479 594 d0 /
      data xgk (  5) / 0.9616149864 2584251241 8130033660 167 d0 /
      data xgk (  6) / 0.9429745712 2897433941 4011169658 471 d0 /
      data xgk (  7) / 0.9207471152 8170156174 6346084546 331 d0 /
      data xgk (  8) / 0.8949919978 7827536885 1042006782 805 d0 /
      data xgk (  9) / 0.8658470652 9327559544 8996969588 340 d0 /
      data xgk ( 10) / 0.8334426287 6083400142 1021108693 570 d0 /
      data xgk ( 11) / 0.7978737979 9850005941 0410904994 307 d0 /
      data xgk ( 12) / 0.7592592630 3735763057 7282865204 361 d0 /
      data xgk ( 13) / 0.7177664068 1308438818 6654079773 298 d0 /
      data xgk ( 14) / 0.6735663684 7346836448 5120633247 622 d0 /
      data xgk ( 15) / 0.6268100990 1031741278 8122681624 518 d0 /
      data xgk ( 16) / 0.5776629302 4122296772 3689841612 654 d0 /
      data xgk ( 17) / 0.5263252843 3471918259 9623778158 010 d0 /
      data xgk ( 18) / 0.4730027314 4571496052 2182115009 192 d0 /
      data xgk ( 19) / 0.4178853821 9303774885 1814394594 572 d0 /
      data xgk ( 20) / 0.3611723058 0938783773 5821730127 641 d0 /
      data xgk ( 21) / 0.3030895389 3110783016 7478909980 339 d0 /
      data xgk ( 22) / 0.2438668837 2098843204 5190362797 452 d0 /
      data xgk ( 23) / 0.1837189394 2104889201 5969888759 528 d0 /
      data xgk ( 24) / 0.1228646926 1071039638 7359818808 037 d0 /
      data xgk ( 25) / 0.0615444830 0568507888 6546392366 797 d0 /
      data xgk ( 26) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0019873838 9233031592 6507851882 843 d0 /
      data wgk (  2) / 0.0055619321 3535671375 8040236901 066 d0 /
      data wgk (  3) / 0.0094739733 8617415160 7207710523 655 d0 /
      data wgk (  4) / 0.0132362291 9557167481 3656405846 976 d0 /
      data wgk (  5) / 0.0168478177 0912829823 1516667536 336 d0 /
      data wgk (  6) / 0.0204353711 4588283545 6568292235 939 d0 /
      data wgk (  7) / 0.0240099456 0695321622 0092489164 881 d0 /
      data wgk (  8) / 0.0274753175 8785173780 2948455517 811 d0 /
      data wgk (  9) / 0.0307923001 6738748889 1109020215 229 d0 /
      data wgk ( 10) / 0.0340021302 7432933783 6748795229 551 d0 /
      data wgk ( 11) / 0.0371162714 8341554356 0330625367 620 d0 /
      data wgk ( 12) / 0.0400838255 0403238207 4839284467 076 d0 /
      data wgk ( 13) / 0.0428728450 2017004947 6895792439 495 d0 /
      data wgk ( 14) / 0.0455029130 4992178890 9870584752 660 d0 /
      data wgk ( 15) / 0.0479825371 3883671390 6392255756 915 d0 /
      data wgk ( 16) / 0.0502776790 8071567196 3325259433 440 d0 /
      data wgk ( 17) / 0.0523628858 0640747586 4366712137 873 d0 /
      data wgk ( 18) / 0.0542511298 8854549014 4543370459 876 d0 /
      data wgk ( 19) / 0.0559508112 2041231730 8240686382 747 d0 /
      data wgk ( 20) / 0.0574371163 6156783285 3582693939 506 d0 /
      data wgk ( 21) / 0.0586896800 2239420796 1974175856 788 d0 /
      data wgk ( 22) / 0.0597203403 2417405997 9099291932 562 d0 /
      data wgk ( 23) / 0.0605394553 7604586294 5360267517 565 d0 /
      data wgk ( 24) / 0.0611285097 1705304830 5859030416 293 d0 /
      data wgk ( 25) / 0.0614711898 7142531666 1544131965 264 d0 /
c       note: wgk (26) was calculated from the values of wgk(1..25)
      data wgk ( 26) / 0.0615808180 6783293507 8759824240 066 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 25-point gauss formula
c           resk   - result of the 51-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk51
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 51-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = dabs(resk)
      do 10 j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(26)*dabs(fc-reskh)
      do 20 j=1,25
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqk61(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk61
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  61-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of dabs(f) over (a,b)
c***description
c
c        integration rule
c        standard fortran subroutine
c        double precision version
c
c
c        parameters
c         on entry
c           f      - double precision
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to be
c                    declared e x t e r n a l in the calling program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c         on return
c           result - double precision
c                    approximation to the integral i
c                    result is computed by applying the 61-point
c                    kronrod rule (resk) obtained by optimal addition of
c                    abscissae to the 30-point gauss rule (resg).
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed dabs(i-result)
c
c           resabs - double precision
c                    approximation to the integral j
c
c           resasc - double precision
c                    approximation to the integral of dabs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk61
c
      double precision a,dabsc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     *  resg,resk,reskh,result,uflow,wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
c
c           the abscissae and weights are given for the
c           interval (-1,1). because of symmetry only the positive
c           abscissae and their corresponding weights are given.
c
c           xgk   - abscissae of the 61-point kronrod rule
c                   xgk(2), xgk(4)  ... abscissae of the 30-point
c                   gauss rule
c                   xgk(1), xgk(3)  ... optimally added abscissae
c                   to the 30-point gauss rule
c
c           wgk   - weights of the 61-point kronrod rule
c
c           wg    - weigths of the 30-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0079681924 9616660561 5465883474 674 d0 /
      data wg  (  2) / 0.0184664683 1109095914 2302131912 047 d0 /
      data wg  (  3) / 0.0287847078 8332336934 9719179611 292 d0 /
      data wg  (  4) / 0.0387991925 6962704959 6801936446 348 d0 /
      data wg  (  5) / 0.0484026728 3059405290 2938140422 808 d0 /
      data wg  (  6) / 0.0574931562 1761906648 1721689402 056 d0 /
      data wg  (  7) / 0.0659742298 8218049512 8128515115 962 d0 /
      data wg  (  8) / 0.0737559747 3770520626 8243850022 191 d0 /
      data wg  (  9) / 0.0807558952 2942021535 4694938460 530 d0 /
      data wg  ( 10) / 0.0868997872 0108297980 2387530715 126 d0 /
      data wg  ( 11) / 0.0921225222 3778612871 7632707087 619 d0 /
      data wg  ( 12) / 0.0963687371 7464425963 9468626351 810 d0 /
      data wg  ( 13) / 0.0995934205 8679526706 2780282103 569 d0 /
      data wg  ( 14) / 0.1017623897 4840550459 6428952168 554 d0 /
      data wg  ( 15) / 0.1028526528 9355884034 1285636705 415 d0 /
c
      data xgk (  1) / 0.9994844100 5049063757 1325895705 811 d0 /
      data xgk (  2) / 0.9968934840 7464954027 1630050918 695 d0 /
      data xgk (  3) / 0.9916309968 7040459485 8628366109 486 d0 /
      data xgk (  4) / 0.9836681232 7974720997 0032581605 663 d0 /
      data xgk (  5) / 0.9731163225 0112626837 4693868423 707 d0 /
      data xgk (  6) / 0.9600218649 6830751221 6871025581 798 d0 /
      data xgk (  7) / 0.9443744447 4855997941 5831324037 439 d0 /
      data xgk (  8) / 0.9262000474 2927432587 9324277080 474 d0 /
      data xgk (  9) / 0.9055733076 9990779854 6522558925 958 d0 /
      data xgk ( 10) / 0.8825605357 9205268154 3116462530 226 d0 /
      data xgk ( 11) / 0.8572052335 4606109895 8658510658 944 d0 /
      data xgk ( 12) / 0.8295657623 8276839744 2898119732 502 d0 /
      data xgk ( 13) / 0.7997278358 2183908301 3668942322 683 d0 /
      data xgk ( 14) / 0.7677774321 0482619491 7977340974 503 d0 /
      data xgk ( 15) / 0.7337900624 5322680472 6171131369 528 d0 /
      data xgk ( 16) / 0.6978504947 9331579693 2292388026 640 d0 /
      data xgk ( 17) / 0.6600610641 2662696137 0053668149 271 d0 /
      data xgk ( 18) / 0.6205261829 8924286114 0477556431 189 d0 /
      data xgk ( 19) / 0.5793452358 2636169175 6024932172 540 d0 /
      data xgk ( 20) / 0.5366241481 4201989926 4169793311 073 d0 /
      data xgk ( 21) / 0.4924804678 6177857499 3693061207 709 d0 /
      data xgk ( 22) / 0.4470337695 3808917678 0609900322 854 d0 /
      data xgk ( 23) / 0.4004012548 3039439253 5476211542 661 d0 /
      data xgk ( 24) / 0.3527047255 3087811347 1037207089 374 d0 /
      data xgk ( 25) / 0.3040732022 7362507737 2677107199 257 d0 /
      data xgk ( 26) / 0.2546369261 6788984643 9805129817 805 d0 /
      data xgk ( 27) / 0.2045251166 8230989143 8957671002 025 d0 /
      data xgk ( 28) / 0.1538699136 0858354696 3794672743 256 d0 /
      data xgk ( 29) / 0.1028069379 6673703014 7096751318 001 d0 /
      data xgk ( 30) / 0.0514718425 5531769583 3025213166 723 d0 /
      data xgk ( 31) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0013890136 9867700762 4551591226 760 d0 /
      data wgk (  2) / 0.0038904611 2709988405 1267201844 516 d0 /
      data wgk (  3) / 0.0066307039 1593129217 3319826369 750 d0 /
      data wgk (  4) / 0.0092732796 5951776342 8441146892 024 d0 /
      data wgk (  5) / 0.0118230152 5349634174 2232898853 251 d0 /
      data wgk (  6) / 0.0143697295 0704580481 2451432443 580 d0 /
      data wgk (  7) / 0.0169208891 8905327262 7572289420 322 d0 /
      data wgk (  8) / 0.0194141411 9394238117 3408951050 128 d0 /
      data wgk (  9) / 0.0218280358 2160919229 7167485738 339 d0 /
      data wgk ( 10) / 0.0241911620 7808060136 5686370725 232 d0 /
      data wgk ( 11) / 0.0265099548 8233310161 0601709335 075 d0 /
      data wgk ( 12) / 0.0287540487 6504129284 3978785354 334 d0 /
      data wgk ( 13) / 0.0309072575 6238776247 2884252943 092 d0 /
      data wgk ( 14) / 0.0329814470 5748372603 1814191016 854 d0 /
      data wgk ( 15) / 0.0349793380 2806002413 7499670731 468 d0 /
      data wgk ( 16) / 0.0368823646 5182122922 3911065617 136 d0 /
      data wgk ( 17) / 0.0386789456 2472759295 0348651532 281 d0 /
      data wgk ( 18) / 0.0403745389 5153595911 1995279752 468 d0 /
      data wgk ( 19) / 0.0419698102 1516424614 7147541285 970 d0 /
      data wgk ( 20) / 0.0434525397 0135606931 6831728117 073 d0 /
      data wgk ( 21) / 0.0448148001 3316266319 2355551616 723 d0 /
      data wgk ( 22) / 0.0460592382 7100698811 6271735559 374 d0 /
      data wgk ( 23) / 0.0471855465 6929915394 5261478181 099 d0 /
      data wgk ( 24) / 0.0481858617 5708712914 0779492298 305 d0 /
      data wgk ( 25) / 0.0490554345 5502977888 7528165367 238 d0 /
      data wgk ( 26) / 0.0497956834 2707420635 7811569379 942 d0 /
      data wgk ( 27) / 0.0504059214 0278234684 0893085653 585 d0 /
      data wgk ( 28) / 0.0508817958 9874960649 2297473049 805 d0 /
      data wgk ( 29) / 0.0512215478 4925877217 0656282604 944 d0 /
      data wgk ( 30) / 0.0514261285 3745902593 3862879215 781 d0 /
      data wgk ( 31) / 0.0514947294 2945156755 8340433647 099 d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           dabsc  - abscissa
c           fval*  - function value
c           resg   - result of the 30-point gauss rule
c           resk   - result of the 61-point kronrod rule
c           reskh  - approximation to the mean value of f
c                    over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(b+a)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 61-point kronrod approximation to the
c           integral, and estimate the absolute error.
c
c***first executable statement  dqk61
      resg = 0.0d+00
      fc = f(centr)
      resk = wgk(31)*fc
      resabs = dabs(resk)
      do 10 j=1,15
        jtw = j*2
        dabsc = hlgth*xgk(jtw)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
  15    continue
      reskh = resk*0.5d+00
      resasc = wgk(31)*dabs(fc-reskh)
      do 20 j=1,30
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c***begin prologue  dqmomo
c***date written   820101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,c3a2
c***keywords  modified chebyshev moments
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes modified chebsyshev moments. the k-th
c            modified chebyshev moment is defined as the integral over
c            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
c            polynomial of degree k.
c***description
c
c        modified chebyshev moments
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           alfa   - double precision
c                    parameter in the weight function w(x), alfa.gt.(-1)
c
c           beta   - double precision
c                    parameter in the weight function w(x), beta.gt.(-1)
c
c           ri     - double precision
c                    vector of dimension 25
c                    ri(k) is the integral over (-1,1) of
c                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
c
c           rj     - double precision
c                    vector of dimension 25
c                    rj(k) is the integral over (-1,1) of
c                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
c
c           rg     - double precision
c                    vector of dimension 25
c                    rg(k) is the integral over (-1,1) of
c                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           rh     - double precision
c                    vector of dimension 25
c                    rh(k) is the integral over (-1,1) of
c                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           integr - integer
c                    input parameter indicating the modified
c                    moments to be computed
c                    integr = 1 compute ri, rj
c                           = 2 compute ri, rj, rg
c                           = 3 compute ri, rj, rh
c                           = 4 compute ri, rj, rg, rh
c
c***references  (none)
c***routines called  (none)
c***end prologue  dqmomo
c
      double precision alfa,alfp1,alfp2,an,anm1,beta,betp1,betp2,ralf,
     *  rbet,rg,rh,ri,rj
      integer i,im1,integr
c
      dimension rg(25),rh(25),ri(25),rj(25)
c
c
c***first executable statement  dqmomo
      alfp1 = alfa+0.1d+01
      betp1 = beta+0.1d+01
      alfp2 = alfa+0.2d+01
      betp2 = beta+0.2d+01
      ralf = 0.2d+01**alfp1
      rbet = 0.2d+01**betp1
c
c           compute ri, rj using a forward recurrence relation.
c
      ri(1) = ralf/alfp1
      rj(1) = rbet/betp1
      ri(2) = ri(1)*alfa/alfp2
      rj(2) = rj(1)*beta/betp2
      an = 0.2d+01
      anm1 = 0.1d+01
      do 20 i=3,25
        ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
        rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
   20 continue
      if(integr.eq.1) go to 70
      if(integr.eq.3) go to 40
c
c           compute rg using a forward recurrence relation.
c
      rg(1) = -ri(1)/alfp1
      rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 30 i=3,25
        rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/
     *  (anm1*(an+alfp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   30 continue
      if(integr.eq.2) go to 70
c
c           compute rh using a forward recurrence relation.
c
   40 rh(1) = -rj(1)/betp1
      rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 50 i=3,25
        rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+
     *  anm1*rj(i))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   50 continue
      do 60 i=2,25,2
        rh(i) = -rh(i)
   60 continue
   70 do 80 i=2,25,2
        rj(i) = -rj(i)
   80 continue
   90 return
      end
      subroutine dqng(f,a,b,epsabs,epsrel,result,abserr,neval,ier)
c***begin prologue  dqng
c***date written   800101   (yymmdd)
c***revision date  810101   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, smooth integrand,
c             non-adaptive, gauss-kronrod(patterson)
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl math & progr. div. - k.u.leuven
c           kahaner,david,nbs - modified (2/82)
c***purpose  the routine calculates an approximation result to a
c            given definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c non-adaptive integration
c standard fortran subroutine
c double precision version
c
c           f      - double precision
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l in the driver program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c           epsabs - double precision
c                    absolute accuracy requested
c           epsrel - double precision
c                    relative accuracy requested
c                    if  epsabs.le.0
c                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                    the routine will end with ier = 6.
c
c         on return
c           result - double precision
c                    approximation to the integral i
c                    result is obtained by applying the 21-point
c                    gauss-kronrod rule (res21) obtained by optimal
c                    addition of abscissae to the 10-point gauss rule
c                    (res10), or by applying the 43-point rule (res43)
c                    obtained by optimal addition of abscissae to the
c                    21-point gauss-kronrod rule, or by applying the
c                    87-point rule (res87) obtained by optimal addition
c                    of abscissae to the 43-point rule.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           ier    - ier = 0 normal and reliable termination of the
c                            routine. it is assumed that the requested
c                            accuracy has been achieved.
c                    ier.gt.0 abnormal termination of the routine. it is
c                            assumed that the requested accuracy has
c                            not been achieved.
c           error messages
c                    ier = 1 the maximum number of steps has been
c                            executed. the integral is probably too
c                            difficult to be calculated by dqng.
c                        = 6 the input is invalid, because
c                            epsabs.le.0 and
c                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                            result, abserr and neval are set to zero.
c
c***references  (none)
c***routines called  d1mach,xerror
c***end prologue  dqng
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,epsabs,epsrel,f,fcentr,fval,fval1,fval2,fv1,fv2,
     *  fv3,fv4,hlgth,result,res10,res21,res43,res87,resabs,resasc,
     *  reskh,savfun,uflow,w10,w21a,w21b,w43a,w43b,w87a,w87b,x1,x2,x3,x4
      integer ier,ipx,k,l,neval
      external f
c
      dimension fv1(5),fv2(5),fv3(5),fv4(5),x1(5),x2(5),x3(11),x4(22),
     *  w10(5),w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23),
     *  savfun(21)
c
c           the following data statements contain the
c           abscissae and weights of the integration rules used.
c
c           x1      abscissae common to the 10-, 21-, 43- and 87-
c                   point rule
c           x2      abscissae common to the 21-, 43- and 87-point rule
c           x3      abscissae common to the 43- and 87-point rule
c           x4      abscissae of the 87-point rule
c           w10     weights of the 10-point formula
c           w21a    weights of the 21-point formula for abscissae x1
c           w21b    weights of the 21-point formula for abscissae x2
c           w43a    weights of the 43-point formula for abscissae x1, x3
c           w43b    weights of the 43-point formula for abscissae x3
c           w87a    weights of the 87-point formula for abscissae x1,
c                   x2, x3
c           w87b    weights of the 87-point formula for abscissae x4
c
c
c gauss-kronrod-patterson quadrature coefficients for use in
c quadpack routine qng.  these coefficients were calculated with
c 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
c
      data x1    (  1) / 0.9739065285 1717172007 7964012084 452 d0 /
      data x1    (  2) / 0.8650633666 8898451073 2096688423 493 d0 /
      data x1    (  3) / 0.6794095682 9902440623 4327365114 874 d0 /
      data x1    (  4) / 0.4333953941 2924719079 9265943165 784 d0 /
      data x1    (  5) / 0.1488743389 8163121088 4826001129 720 d0 /
      data w10   (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data w10   (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data w10   (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data w10   (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data w10   (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data x2    (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data x2    (  2) / 0.9301574913 5570822600 1207180059 508 d0 /
      data x2    (  3) / 0.7808177265 8641689706 3717578345 042 d0 /
      data x2    (  4) / 0.5627571346 6860468333 9000099272 694 d0 /
      data x2    (  5) / 0.2943928627 0146019813 1126603103 866 d0 /
      data w21a  (  1) / 0.0325581623 0796472747 8818972459 390 d0 /
      data w21a  (  2) / 0.0750396748 1091995276 7043140916 190 d0 /
      data w21a  (  3) / 0.1093871588 0229764189 9210590325 805 d0 /
      data w21a  (  4) / 0.1347092173 1147332592 8054001771 707 d0 /
      data w21a  (  5) / 0.1477391049 0133849137 4841515972 068 d0 /
      data w21b  (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data w21b  (  2) / 0.0547558965 7435199603 1381300244 580 d0 /
      data w21b  (  3) / 0.0931254545 8369760553 5065465083 366 d0 /
      data w21b  (  4) / 0.1234919762 6206585107 7958109831 074 d0 /
      data w21b  (  5) / 0.1427759385 7706008079 7094273138 717 d0 /
      data w21b  (  6) / 0.1494455540 0291690566 4936468389 821 d0 /
c
      data x3    (  1) / 0.9993333609 0193208139 4099323919 911 d0 /
      data x3    (  2) / 0.9874334029 0808886979 5961478381 209 d0 /
      data x3    (  3) / 0.9548079348 1426629925 7919200290 473 d0 /
      data x3    (  4) / 0.9001486957 4832829362 5099494069 092 d0 /
      data x3    (  5) / 0.8251983149 8311415084 7066732588 520 d0 /
      data x3    (  6) / 0.7321483889 8930498261 2354848755 461 d0 /
      data x3    (  7) / 0.6228479705 3772523864 1159120344 323 d0 /
      data x3    (  8) / 0.4994795740 7105649995 2214885499 755 d0 /
      data x3    (  9) / 0.3649016613 4658076804 3989548502 644 d0 /
      data x3    ( 10) / 0.2222549197 7660129649 8260928066 212 d0 /
      data x3    ( 11) / 0.0746506174 6138332204 3914435796 506 d0 /
      data w43a  (  1) / 0.0162967342 8966656492 4281974617 663 d0 /
      data w43a  (  2) / 0.0375228761 2086950146 1613795898 115 d0 /
      data w43a  (  3) / 0.0546949020 5825544214 7212685465 005 d0 /
      data w43a  (  4) / 0.0673554146 0947808607 5553166302 174 d0 /
      data w43a  (  5) / 0.0738701996 3239395343 2140695251 367 d0 /
      data w43a  (  6) / 0.0057685560 5976979618 4184327908 655 d0 /
      data w43a  (  7) / 0.0273718905 9324884208 1276069289 151 d0 /
      data w43a  (  8) / 0.0465608269 1042883074 3339154433 824 d0 /
      data w43a  (  9) / 0.0617449952 0144256449 6240336030 883 d0 /
      data w43a  ( 10) / 0.0713872672 6869339776 8559114425 516 d0 /
      data w43b  (  1) / 0.0018444776 4021241410 0389106552 965 d0 /
      data w43b  (  2) / 0.0107986895 8589165174 0465406741 293 d0 /
      data w43b  (  3) / 0.0218953638 6779542810 2523123075 149 d0 /
      data w43b  (  4) / 0.0325974639 7534568944 3882222526 137 d0 /
      data w43b  (  5) / 0.0421631379 3519181184 7627924327 955 d0 /
      data w43b  (  6) / 0.0507419396 0018457778 0189020092 084 d0 /
      data w43b  (  7) / 0.0583793955 4261924837 5475369330 206 d0 /
      data w43b  (  8) / 0.0647464049 5144588554 4689259517 511 d0 /
      data w43b  (  9) / 0.0695661979 1235648452 8633315038 405 d0 /
      data w43b  ( 10) / 0.0728244414 7183320815 0939535192 842 d0 /
      data w43b  ( 11) / 0.0745077510 1417511827 3571813842 889 d0 /
      data w43b  ( 12) / 0.0747221475 1740300559 4425168280 423 d0 /
c
      data x4    (  1) / 0.9999029772 6272923449 0529830591 582 d0 /
      data x4    (  2) / 0.9979898959 8667874542 7496322365 960 d0 /
      data x4    (  3) / 0.9921754978 6068722280 8523352251 425 d0 /
      data x4    (  4) / 0.9813581635 7271277357 1916941623 894 d0 /
      data x4    (  5) / 0.9650576238 5838461912 8284110607 926 d0 /
      data x4    (  6) / 0.9431676131 3367059681 6416634507 426 d0 /
      data x4    (  7) / 0.9158064146 8550720959 1826430720 050 d0 /
      data x4    (  8) / 0.8832216577 7131650137 2117548744 163 d0 /
      data x4    (  9) / 0.8457107484 6241566660 5902011504 855 d0 /
      data x4    ( 10) / 0.8035576580 3523098278 8739474980 964 d0 /
      data x4    ( 11) / 0.7570057306 8549555832 8942793432 020 d0 /
      data x4    ( 12) / 0.7062732097 8732181982 4094274740 840 d0 /
      data x4    ( 13) / 0.6515894665 0117792253 4422205016 736 d0 /
      data x4    ( 14) / 0.5932233740 5796108887 5273770349 144 d0 /
      data x4    ( 15) / 0.5314936059 7083193228 5268948562 671 d0 /
      data x4    ( 16) / 0.4667636230 4202284487 1966781659 270 d0 /
      data x4    ( 17) / 0.3994248478 5921880473 2101665817 923 d0 /
      data x4    ( 18) / 0.3298748771 0618828826 5053371824 597 d0 /
      data x4    ( 19) / 0.2585035592 0216155180 2280975429 025 d0 /
      data x4    ( 20) / 0.1856953965 6834665201 5917141167 606 d0 /
      data x4    ( 21) / 0.1118422131 7990746817 2398359241 362 d0 /
      data x4    ( 22) / 0.0373521233 9461987081 4998165437 704 d0 /
      data w87a  (  1) / 0.0081483773 8414917290 0002878448 190 d0 /
      data w87a  (  2) / 0.0187614382 0156282224 3935059003 794 d0 /
      data w87a  (  3) / 0.0273474510 5005228616 1582829741 283 d0 /
      data w87a  (  4) / 0.0336777073 1163793004 6581056957 588 d0 /
      data w87a  (  5) / 0.0369350998 2042790761 4589586742 499 d0 /
      data w87a  (  6) / 0.0028848724 3021153050 1334156248 695 d0 /
      data w87a  (  7) / 0.0136859460 2271270188 8950035273 128 d0 /
      data w87a  (  8) / 0.0232804135 0288831112 3409291030 404 d0 /
      data w87a  (  9) / 0.0308724976 1171335867 5466394126 442 d0 /
      data w87a  ( 10) / 0.0356936336 3941877071 9351355457 044 d0 /
      data w87a  ( 11) / 0.0009152833 4520224136 0843392549 948 d0 /
      data w87a  ( 12) / 0.0053992802 1930047136 7738743391 053 d0 /
      data w87a  ( 13) / 0.0109476796 0111893113 4327826856 808 d0 /
      data w87a  ( 14) / 0.0162987316 9678733526 2665703223 280 d0 /
      data w87a  ( 15) / 0.0210815688 8920383511 2433060188 190 d0 /
      data w87a  ( 16) / 0.0253709697 6925382724 3467999831 710 d0 /
      data w87a  ( 17) / 0.0291896977 5647575250 1446154084 920 d0 /
      data w87a  ( 18) / 0.0323732024 6720278968 5788194889 595 d0 /
      data w87a  ( 19) / 0.0347830989 5036514275 0781997949 596 d0 /
      data w87a  ( 20) / 0.0364122207 3135178756 2801163687 577 d0 /
      data w87a  ( 21) / 0.0372538755 0304770853 9592001191 226 d0 /
      data w87b  (  1) / 0.0002741455 6376207235 0016527092 881 d0 /
      data w87b  (  2) / 0.0018071241 5505794294 8341311753 254 d0 /
      data w87b  (  3) / 0.0040968692 8275916486 4458070683 480 d0 /
      data w87b  (  4) / 0.0067582900 5184737869 9816577897 424 d0 /
      data w87b  (  5) / 0.0095499576 7220164653 6053581325 377 d0 /
      data w87b  (  6) / 0.0123294476 5224485369 4626639963 780 d0 /
      data w87b  (  7) / 0.0150104473 4638895237 6697286041 943 d0 /
      data w87b  (  8) / 0.0175489679 8624319109 9665352925 900 d0 /
      data w87b  (  9) / 0.0199380377 8644088820 2278192730 714 d0 /
      data w87b  ( 10) / 0.0221949359 6101228679 6332102959 499 d0 /
      data w87b  ( 11) / 0.0243391471 2600080547 0360647041 454 d0 /
      data w87b  ( 12) / 0.0263745054 1483920724 1503786552 615 d0 /
      data w87b  ( 13) / 0.0282869107 8877120065 9968002987 960 d0 /
      data w87b  ( 14) / 0.0300525811 2809269532 2521110347 341 d0 /
      data w87b  ( 15) / 0.0316467513 7143992940 4586051078 883 d0 /
      data w87b  ( 16) / 0.0330504134 1997850329 0785944862 689 d0 /
      data w87b  ( 17) / 0.0342550997 0422606178 7082821046 821 d0 /
      data w87b  ( 18) / 0.0352624126 6015668103 3782717998 428 d0 /
      data w87b  ( 19) / 0.0360769896 2288870118 5500318003 895 d0 /
      data w87b  ( 20) / 0.0366986044 9845609449 8018047441 094 d0 /
      data w87b  ( 21) / 0.0371205492 6983257611 4119958413 599 d0 /
      data w87b  ( 22) / 0.0373342287 5193504032 1235449094 698 d0 /
      data w87b  ( 23) / 0.0373610737 6267902341 0321241766 599 d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fcentr - function value at mid point
c           absc   - abscissa
c           fval   - function value
c           savfun - array of function values which have already been
c                    computed
c           res10  - 10-point gauss result
c           res21  - 21-point kronrod result
c           res43  - 43-point result
c           res87  - 87-point result
c           resabs - approximation to the integral of abs(f)
c           resasc - approximation to the integral of abs(f-i/(b-a))
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqng
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      result = 0.0d+00
      abserr = 0.0d+00
      neval = 0
      ier = 6
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     *  go to 80
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      centr = 0.5d+00*(b+a)
      fcentr = f(centr)
      neval = 21
      ier = 1
c
c          compute the integral using the 10- and 21-point formula.
c
      do 70 l = 1,3
      go to (5,25,45),l
    5 res10 = 0.0d+00
      res21 = w21b(6)*fcentr
      resabs = w21b(6)*dabs(fcentr)
      do 10 k=1,5
        absc = hlgth*x1(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1+fval2
        res10 = res10+w10(k)*fval
        res21 = res21+w21a(k)*fval
        resabs = resabs+w21a(k)*(dabs(fval1)+dabs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
   10 continue
      ipx = 5
      do 15 k=1,5
        ipx = ipx+1
        absc = hlgth*x2(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1+fval2
        res21 = res21+w21b(k)*fval
        resabs = resabs+w21b(k)*(dabs(fval1)+dabs(fval2))
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
   15 continue
c
c          test for convergence.
c
      result = res21*hlgth
      resabs = resabs*dhlgth
      reskh = 0.5d+00*res21
      resasc = w21b(6)*dabs(fcentr-reskh)
      do 20 k = 1,5
        resasc = resasc+w21a(k)*(dabs(fv1(k)-reskh)+dabs(fv2(k)-reskh))
     *                  +w21b(k)*(dabs(fv3(k)-reskh)+dabs(fv4(k)-reskh))
   20 continue
      abserr = dabs((res21-res10)*hlgth)
      resasc = resasc*dhlgth
      go to 65
c
c          compute the integral using the 43-point formula.
c
   25 res43 = w43b(12)*fcentr
      neval = 43
      do 30 k=1,10
        res43 = res43+savfun(k)*w43a(k)
   30 continue
      do 40 k=1,11
        ipx = ipx+1
        absc = hlgth*x3(k)
        fval = f(absc+centr)+f(centr-absc)
        res43 = res43+fval*w43b(k)
        savfun(ipx) = fval
   40 continue
c
c          test for convergence.
c
      result = res43*hlgth
      abserr = dabs((res43-res21)*hlgth)
      go to 65
c
c          compute the integral using the 87-point formula.
c
   45 res87 = w87b(23)*fcentr
      neval = 87
      do 50 k=1,21
        res87 = res87+savfun(k)*w87a(k)
   50 continue
      do 60 k=1,22
        absc = hlgth*x4(k)
        res87 = res87+w87b(k)*(f(absc+centr)+f(centr-absc))
   60 continue
      result = res87*hlgth
      abserr = dabs((res87-res43)*hlgth)
   65 if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if (resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      if (abserr.le.dmax1(epsabs,epsrel*dabs(result))) ier = 0
c ***jump out of do-loop
      if (ier.eq.0) go to 999
   70 continue
   80 call xerror(26habnormal return from dqng ,26,ier,0)
  999 return
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  dqpsrt
c***refer to  dqage,dqagie,dqagpe,dqawse
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  sequential sorting
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine maintains the descending ordering in the
c            list of the local error estimated resulting from the
c            interval subdivision process. at each call two error
c            estimates are inserted using the sequential search
c            method, top-down for the largest error estimate and
c            bottom-up for the smallest error estimate.
c***description
c
c           ordering routine
c           standard fortran subroutine
c           double precision version
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - double precision
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - double precision
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k elements
c                       of which contain pointers to the error
c                       estimates, such that
c                       elist(iord(1)),...,  elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c***end prologue  dqpsrt
c
      double precision elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed if, due to a
c           difficult integrand, subdivision increased the error
c           estimate. in the normal case the insert procedure should
c           start after the nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to be maintained
c           in descending order. this number depends on the number of
c           subdivisions still allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
      double precision function dqwgtc(x,c,p2,p3,p4,kp)
c***begin prologue  dqwgtc
c***refer to dqk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, cauchy principal value
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine qawc and defines the weight function.
c***end prologue  dqwgtc
c
      double precision c,p2,p3,p4,x
      integer kp
c***first executable statement  dqwgtc
      dqwgtc = 0.1d+01/(x-c)
      return
      end
      double precision function dqwgtf(x,omega,p2,p3,p4,integr)
c***begin prologue  dqwgtf
c***refer to   dqk15w
c***routines called  (none)
c***revision date 810101   (yymmdd)
c***keywords  cos or sin in weight function
c***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. * progr. div. - k.u.leuven
c***end prologue  dqwgtf
c
      double precision dcos,dsin,omega,omx,p2,p3,p4,x
      integer integr
c***first executable statement  dqwgtf
      omx = omega*x
      go to(10,20),integr
   10 dqwgtf = dcos(omx)
      go to 30
   20 dqwgtf = dsin(omx)
   30 return
      end
      double precision function dqwgts(x,a,b,alfa,beta,integr)
c***begin prologue  dqwgts
c***refer to dqk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, algebraico-logarithmic
c             end-point singularities
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine dqaws and defines the weight function.
c***end prologue  dqwgts
c
      double precision a,alfa,b,beta,bmx,dlog,x,xma
      integer integr
c***first executable statement  dqwgts
      xma = x-a
      bmx = b-x
      dqwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 dqwgts = dqwgts*dlog(xma)
      go to 40
   20 dqwgts = dqwgts*dlog(bmx)
      go to 40
   30 dqwgts = dqwgts*dlog(xma)*dlog(bmx)
   40 return
      end
      subroutine qag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     *    limit,lenw,last,iwork,work)
c***begin prologue  qag
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result)le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c            f      - real
c                     function subprogam defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                       7 - 15 points if key.lt.2,
c                      10 - 21 points if key = 2,
c                      15 - 31 points if key = 3,
c                      20 - 41 points if key = 4,
c                      25 - 51 points if key = 5,
c                      30 - 61 points if key.gt.5.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c                      error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yield no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulaties.
c                             if the position of a local difficulty can
c                             be determined (i.e.singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set
c                             to zero.
c                             except when lenw is invalid, iwork(1),
c                             work(limit*2+1) and work(limit*3+1) are
c                             set to zero, work(1) is set to a and
c                             work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers to the error
c                    estimates over the subintervals, such that
c                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
c                    form a decreasing sequence with k = last if
c                    last.le.(limit/2+2), and k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left end
c                    points of the subintervals in the partition of
c                     (a,b),
c                    work(limit+1), ..., work(limit+last) contain the
c                     right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last) contain
c                     the error estimates.
c
c***references  (none)
c***routines called  qage,xerror
c***end prologue  qag
c
      real a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,key,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of lenw.
c
c***first executable statement  qag
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qage.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call qage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qag ,
     *  26,ier,lvl)
      return
      end
      subroutine qage(f,a,b,epsabs,epsrel,key,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qage
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             integrand examinator, globally adaptive,
c             gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral   i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            key    - integer
c                     key for choice of local integration rule
c                     a gauss-kronrod pair is used with
c                          7 - 15 points if key.lt.2,
c                         10 - 21 points if key = 2,
c                         15 - 31 points if key = 3,
c                         20 - 41 points if key = 4,
c                         25 - 51 points if key = 5,
c                         30 - 61 points if key.gt.5.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for result and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value
c                             of limit.
c                             however, if this yields no improvement it
c                             is rather advised to analyze the integrand
c                             in order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined(e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1) ,
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the
c                      integral approximations on the subintervals
c
c            elist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the moduli of the
c                      absolute error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the
c                      error estimates over the subintervals,
c                      such that elist(iord(1)), ...,
c                      elist(iord(k)) form a decreasing sequence,
c                      with k = last if last.le.(limit/2+2), and
c                      k = limit+1-last otherwise
c
c            last    - integer
c                      number of subintervals actually produced in the
c                      subdivision process
c
c***references  (none)
c***routines called  qk15,qk21,qk31,qk41,qk51,qk61,qpsrt,r1mach
c***end prologue  qage
c
      real a,abserr,alist,area,area1,area12,area2,a1,a2,b,blist,
     *  b1,b2,defabs,defab1,defab2,r1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,
     *  limit,maxerr,neval,nrmax
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                      (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach  is the largest relative spacing.
c           uflow  is the smallest positive magnitude.
c
c***first executable statement  qage
      epmach = r1mach(4)
      uflow = r1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      if(epsabs.le.0.0e+00.and.
     *  epsrel.lt.amax1(0.5e+02*epmach,0.5e-14)) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call qk15(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.2) call qk21(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.3) call qk31(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.4) call qk41(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.5) call qk51(f,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.6) call qk61(f,a,b,result,abserr,defabs,resabs)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
c
c           test on accuracy.
c
      errbnd = amax1(epsabs,epsrel*abs(result))
      if(abserr.le.0.5e+02*epmach*defabs.and.abserr.gt.
     *  errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)
     *  .or.abserr.eq.0.0e+00) go to 60
c
c           initialization
c           --------------
c
c
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 30 last = 2,limit
c
c           bisect the subinterval with the largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call qk15(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.2) call qk21(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.3) call qk31(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.4) call qk41(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.5) call qk51(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.6) call qk61(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.1) call qk15(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.2) call qk21(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.3) call qk31(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.4) call qk41(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.5) call qk51(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.6) call qk61(f,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(abs(rlist(maxerr)-area12).le.0.1e-04*abs(area12)
     *  .and.erro12.ge.0.99e+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 8
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*
     *  epmach)*(abs(a2)+0.1e+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with the largest error estimate (to be
c           bisected next).
c
   20   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
c
c           compute final result.
c           ---------------------
c
   40 result = 0.0e+00
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
      subroutine qagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,
     *   ier,limit,lenw,last,iwork,work)
c***begin prologue  qagi
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. -k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c                    or i = integral of f over (-infinity,bound)
c                    or i = integral of f over (-infinity,+infinity)
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration over infinite intervals
c        standard fortran subroutine
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - real
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - integer
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                              or limit.lt.1 or leniw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. exept when limit or leniw is
c                             invalid, iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first
c                    k elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ...,work(limit*2+last) contain the
c                     integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3)
c                     contain the error estimates.
c***references  (none)
c***routines called  qagie,xerror
c***end prologue  qagi
c
      real   abserr,  epsabs,epsrel,f,result,work
      integer ier,iwork,    lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qagi
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qagie.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call qagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qagi,
     *  26,ier,lvl)
      return
      end
      subroutine qagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qagie
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math & progr. div - k.u.leuven
c           de doncker,elise,appl. math & progr. div - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c                    or i = integral of f over (-infinity,bound)
c                    or i = integral of f over (-infinity,+infinity),
c                    hopefully satisfying following claim for accuracy
c                    abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c integration over infinite intervals
c standard fortran subroutine
c
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - real
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - real
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however,if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1),
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to 0
c                             and 1 respectively.
c
c            alist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            blist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            rlist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real
c                     vector of dimension at least limit,  the first
c                     last elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced
c                     in the subdivision process
c
c***references  (none)
c***routines called  qelg,qk15i,qpsrt,r1mach
c***end prologue  qagie
c
      real abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,blist,boun,bound,b1,b2,correc,defabs,defab1,defab2,
     *  dres,r1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs,
     *  reseps,result,res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,inf,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine qelg.
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least (limexp+2),
c                       containing the part of the epsilon table
c                       wich is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained, it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
       epmach = r1mach(4)
c
c           test on validity of parameters
c           -----------------------------
c
c***first executable statement  qagie
      ier = 0
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      alist(1) = 0.0e+00
      blist(1) = 0.1e+01
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      if(epsabs.le.0.0e+00.and.epsrel.lt.amax1(0.5e+02*epmach,0.5e-14))
     *  ier = 6
      if(ier.eq.6) go to 999
c
c
c           first approximation to the integral
c           -----------------------------------
c
c           determine the interval to be mapped onto (0,1).
c           if inf = 2 the integral is computed as i = i1+i2, where
c           i1 = integral of f over (-infinity,0),
c           i2 = integral of f over (0,+infinity).
c
      boun = bound
      if(inf.eq.2) boun = 0.0e+00
      call qk15i(f,boun,inf,0.0e+00,0.1e+01,result,abserr,
     *  defabs,resabs)
c
c           test on accuracy
c
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = abs(result)
      errbnd = amax1(epsabs,epsrel*dres)
      if(abserr.le.1.0e+02*epmach*defabs.and.abserr.gt.
     *  errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0e+00) go to 130
c
c           initialization
c           --------------
c
      uflow = r1mach(1)
      oflow = r1mach(2)
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
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
      if(dres.ge.(0.1e+01-0.5e+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
        call qk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2)go to 15
        if(abs(rlist(maxerr)-area12).gt.0.1e-04*abs(area12)
     *  .or.erro12.lt.0.99e+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at some points of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)*
     *  (abs(a2)+0.1e+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   30   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
        if(errsum.le.errbnd) go to 115
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors
c           over the larger intervals (erlarg) and perform
c           extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1e-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = amax1(epsabs,epsrel*abs(reseps))
        if(abserr.le.ertest) go to 100
c
c            prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5e+00
        erlarg = errsum
        go to 90
   80   small = 0.375e+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if((ier+ierro).eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0e+00.and.area.ne.0.0e+00)go to 105
      if(abserr.gt.errsum)go to 115
      if(area.eq.0.0e+00) go to 130
      go to 110
  105 if(abserr/abs(result).gt.errsum/abs(area))go to 115
c
c           test on divergence
c
  110 if(ksgn.eq.(-1).and.amax1(abs(result),abs(area)).le.
     * defabs*0.1e-01) go to 130
      if(0.1e-01.gt.(result/area).or.(result/area).gt.0.1e+03.
     *or.errsum.gt.abs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0e+00
      do 120 k = 1,last
        result = result+rlist(k)
  120 continue
      abserr = errsum
  130 neval = 30*last-15
      if(inf.eq.2) neval = 2*neval
      if(ier.gt.2) ier=ier-1
  999 return
      end
      subroutine qagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,
     *   neval,ier,leniw,lenw,last,iwork,work)
c***begin prologue  qagp
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            break points of the integration interval, where local
c            difficulties of the integrand may occur(e.g. singularities,
c            discontinuities), are provided by the user.
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            npts2  - integer
c                     number equal to two more than the number of
c                     user-supplied break points within the integration
c                     range, npts.ge.2.
c                     if npts2.lt.2, the routine will end with ier = 6.
c
c            points - real
c                     vector of dimension npts2, the first (npts2-2)
c                     elements of which are the user provided break
c                     points. if these points do not constitute an
c                     ascending sequence there will be an automatic
c                     sorting.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (i.e. singularity,
c                             discontinuity within the interval), it
c                             should be supplied to the routine as an
c                             element of the vector points. if necessary
c                             an appropriate special-purpose integrator
c                             must be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that
c                             the returned result is the best which
c                             can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid because
c                             npts2.lt.2 or
c                             break points are specified outside
c                             the integration range or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             result, abserr, neval, last are set to
c                             zero. exept when leniw or lenw or npts2 is
c                             invalid, iwork(1), iwork(limit+1),
c                             work(limit*2+1) and work(limit*3+1)
c                             are set to zero.
c                             work(1) is set to a and work(limit+1)
c                             to b (where limit = (leniw-npts2)/2).
c
c         dimensioning parameters
c            leniw - integer
c                    dimensioning parameter for iwork
c                    leniw determines limit = (leniw-npts2)/2,
c                    which is the maximum number of subintervals in the
c                    partition of the given integration interval (a,b),
c                    leniw.ge.(3*npts2-2).
c                    if leniw.lt.(3*npts2-2), the routine will end with
c                    ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least leniw*2-npts2.
c                    if lenw.lt.leniw*2-npts2, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least leniw. on return,
c                    the first k elements of which contain
c                    pointers to the error estimates over the
c                    subintervals, such that work(limit*3+iwork(1)),...,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c                    iwork(limit+1), ...,iwork(limit+last) contain the
c                     subdivision levels of the subintervals, i.e.
c                     if (aa,bb) is a subinterval of (p1,p2)
c                     where p1 as well as p2 is a user-provided
c                     break point or integration limit, then (aa,bb) has
c                     level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
c                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have
c                     no significance for the user,
c                    note that limit = (leniw-npts2)/2.
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the corresponding error estimates,
c                    work(limit*4+1), ..., work(limit*4+npts2)
c                     contain the integration limits and the
c                     break points sorted in an ascending sequence.
c                    note that limit = (leniw-npts2)/2.
c
c***references  (none)
c***routines called  qagpe,xerror
c***end prologue  qagp
c
      real a,abserr,b,epsabs,epsrel,f,points,result,work
      integer ier,iwork,leniw,lenw,limit,lvl,l1,l2,l3,neval,npts2
c
      dimension iwork(leniw),points(npts2),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qagp
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(leniw.lt.(3*npts2-2).or.lenw.lt.(leniw*2-npts2).or.npts2.lt.2)
     *  go to 10
c
c         prepare call for qagpe.
c
      limit = (leniw-npts2)/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
c
      call qagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),work(l4),
     *  iwork(1),iwork(l1),iwork(l2),last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qagp,
     *  26,ier,lvl)
      return
      end
      subroutine qagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,
     *   abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin,
     *   last)
c***begin prologue  qagpe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, general-purpose,
c             singularities at user specified points,
c             extrapolation, globally adaptive.
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),hopefully
c            satisfying following claim for accuracy abs(i-result).le.
c            max(epsabs,epsrel*abs(i)). break points of the integration
c            interval, where local difficulties of the integrand may
c            occur(e.g. singularities,discontinuities),provided by user.
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            npts2  - integer
c                     number equal to two more than the number of
c                     user-supplied break points within the integration
c                     range, npts2.ge.2.
c                     if npts2.lt.2, the routine will end with ier = 6.
c
c            points - real
c                     vector of dimension npts2, the first (npts2-2)
c                     elements of which are the user provided break
c                     points. if these points do not constitute an
c                     ascending sequence there will be an automatic
c                     sorting.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.npts2
c                     if limit.lt.npts2, the routine will end with
c                     ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (i.e. singularity,
c                             discontinuity within the interval), it
c                             should be supplied to the routine as an
c                             element of the vector points. if necessary
c                             an appropriate special-purpose integrator
c                             must be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid because
c                             npts2.lt.2 or
c                             break points are specified outside
c                             the integration range or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.npts2.
c                             result, abserr, neval, last, rlist(1),
c                             and elist(1) are set to zero. alist(1) and
c                             blist(1) are set to a and b respectively.
c
c            alist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            blist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            pts    - real
c                     vector of dimension at least npts2, containing the
c                     integration limits and the break points of the
c                     interval in ascending sequence.
c
c            level  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subinterval, i.e. if
c                     (aa,bb) is a subinterval of (p1,p2) where p1 as
c                     well as p2 is a user-provided break point or
c                     integration limit, then (aa,bb) has level l if
c                     abs(bb-aa) = abs(p2-p1)*2**(-l).
c
c            ndin   - integer
c                     vector of dimension at least npts2, after first
c                     integration over the intervals (pts(i)),pts(i+1),
c                     i = 0,1, ..., npts2-2, the error estimates over
c                     some of the intervals may have been increased
c                     artificially, in order to put their subdivision
c                     forward. if this happens for the subinterval
c                     numbered k, ndin(k) is put to 1, otherwise
c                     ndin(k) = 0.
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivisions process
c
c***references  (none)
c***routines called  qelg,qk21,qpsrt,r1mach
c***end prologue  qagpe
      real a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,defabs,defab1,defab2,
     *  dres,r1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,
     *  errmax,error1,erro12,error2,errsum,ertest,f,oflow,points,pts,
     *  resa,resabs,reseps,result,res3la,rlist,rlist2,sign,temp,
     *  uflow
      integer i,id,ier,ierro,ind1,ind2,iord,ip1,iroff1,iroff2,
     *  iroff3,j,jlow,jupbnd,k,ksgn,ktmin,last,levcur,level,levmax,
     *  limit,maxerr,ndin,neval,nint,nintp1,npts,npts2,nres,
     *  nrmax,numrl2
      logical extrap,noext
c
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  level(limit),ndin(npts2),points(npts2),pts(npts2),res3la(3),
     *  rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine epsalg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table which
c                       is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained, it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation is
c                       no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  qagpe
      epmach = r1mach(4)
c
c            test on validity of parameters
c            -----------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      level(1) = 0
      npts = npts2-2
      if(npts2.lt.2.or.limit.le.npts.or.(epsabs.le.0.0e+00.and.
     *  epsrel.lt.amax1(0.5e+02*epmach,0.5e-14))) ier = 6
      if(ier.eq.6) go to 210
c
c            if any break points are provided, sort them into an
c            ascending sequence.
c
      sign = 1.0e+00
      if(a.gt.b) sign = -1.0e+00
      pts(1) = amin1(a,b)
      if(npts.eq.0) go to 15
      do 10 i = 1,npts
        pts(i+1) = points(i)
   10 continue
   15 pts(npts+2) = amax1(a,b)
      nint = npts+1
      a1 = pts(1)
      if(npts.eq.0) go to 40
      nintp1 = nint+1
      do 20 i = 1,nint
        ip1 = i+1
        do 20 j = ip1,nintp1
          if(pts(i).le.pts(j)) go to 20
          temp = pts(i)
          pts(i) = pts(j)
          pts(j) = temp
   20 continue
      if(pts(1).ne.amin1(a,b).or.pts(nintp1).ne.
     *  amax1(a,b)) ier = 6
      if(ier.eq.6) go to 999
c
c            compute first integral and error approximations.
c            ------------------------------------------------
c
   40 resabs = 0.0e+00
      do 50 i = 1,nint
        b1 = pts(i+1)
        call qk21(f,a1,b1,area1,error1,defabs,resa)
        abserr = abserr+error1
        result = result+area1
        ndin(i) = 0
        if(error1.eq.resa.and.error1.ne.0.0e+00) ndin(i) = 1
        resabs = resabs+defabs
        level(i) = 0
        elist(i) = error1
        alist(i) = a1
        blist(i) = b1
        rlist(i) = area1
        iord(i) = i
        a1 = b1
   50 continue
      errsum = 0.0e+00
      do 55 i = 1,nint
        if(ndin(i).eq.1) elist(i) = abserr
        errsum = errsum+elist(i)
   55 continue
c
c           test on accuracy.
c
      last = nint
      neval = 21*nint
      dres = abs(result)
      errbnd = amax1(epsabs,epsrel*dres)
      if(abserr.le.0.1e+03*epmach*resabs.and.abserr.gt.
     *  errbnd) ier = 2
      if(nint.eq.1) go to 80
      do 70 i = 1,npts
        jlow = i+1
        ind1 = iord(i)
        do 60 j = jlow,nint
          ind2 = iord(j)
          if(elist(ind1).gt.elist(ind2)) go to 60
          ind1 = ind2
          k = j
   60   continue
        if(ind1.eq.iord(i)) go to 70
        iord(k) = iord(i)
        iord(i) = ind1
   70 continue
      if(limit.lt.npts2) ier = 1
   80 if(ier.ne.0.or.abserr.le.errbnd) go to 999
c
c           initialization
c           --------------
c
      rlist2(1) = result
      maxerr = iord(1)
      errmax = elist(maxerr)
      area = result
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
      uflow = r1mach(1)
      oflow = r1mach(2)
      abserr = oflow
      ksgn = -1
      if(dres.ge.(0.1e+01-0.5e+02*epmach)*resabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 160 last = npts2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        levcur = level(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk21(f,a1,b1,area1,error1,resa,defab1)
        call qk21(f,a2,b2,area2,error2,resa,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        neval = neval+42
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 95
        if(abs(rlist(maxerr)-area12).gt.0.1e-04*abs(area12)
     *  .or.erro12.lt.0.99e+00*errmax) go to 90
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   90   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   95   level(maxerr) = levcur
        level(last) = levcur
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)*
     *  (abs(a2)+0.1e+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 100
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 110
  100   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
  110   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 190
c ***jump out of do-loop
        if(ier.ne.0) go to 170
        if(noext) go to 160
        erlarg = erlarg-erlast
        if(levcur+1.le.levmax) erlarg = erlarg+erro12
        if(extrap) go to 120
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(level(maxerr)+1.le.levmax) go to 160
        extrap = .true.
        nrmax = 2
  120   if(ierro.eq.3.or.erlarg.le.ertest) go to 140
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors
c           over the larger intervals (erlarg) and perform
c           extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 130 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(level(maxerr)+1.le.levmax) go to 160
          nrmax = nrmax+1
  130   continue
c
c           perform extrapolation.
c
  140   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.le.2) go to 155
        call qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1e-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 150
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = amax1(epsabs,epsrel*abs(reseps))
c ***jump out of do-loop
        if(abserr.lt.ertest) go to 170
c
c           prepare bisection of the smallest interval.
c
  150   if(numrl2.eq.1) noext = .true.
        if(ier.ge.5) go to 170
  155   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        levmax = levmax+1
        erlarg = errsum
  160 continue
c
c           set the final result.
c           ---------------------
c
c
  170 if(abserr.eq.oflow) go to 190
      if((ier+ierro).eq.0) go to 180
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0e+00.and.area.ne.0.0e+00)go to 175
      if(abserr.gt.errsum)go to 190
      if(area.eq.0.0e+00) go to 210
      go to 180
  175 if(abserr/abs(result).gt.errsum/abs(area))go to 190
c
c           test on divergence.
c
  180 if(ksgn.eq.(-1).and.amax1(abs(result),abs(area)).le.
     *  resabs*0.1e-01) go to 210
      if(0.1e-01.gt.(result/area).or.(result/area).gt.0.1e+03.or.
     *  errsum.gt.abs(area)) ier = 6
      go to 210
c
c           compute global integral sum.
c
  190 result = 0.0e+00
      do 200 k = 1,last
        result = result+rlist(k)
  200 continue
      abserr = errsum
  210 if(ier.gt.2) ier = ier - 1
      result = result*sign
 999  return
      end
      subroutine qags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)
c***begin prologue  qags
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end-point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & prog. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral  i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account. however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be
c                             achieved, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero.except when limit or lenw is invalid,
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, detemines the
c                    number of significant elements actually in the work
c                    arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end-points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end-points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c
c
c***references  (none)
c***routines called  qagse,xerror
c***end prologue  qags
c
c
      real a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qags
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qagse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call qagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *  ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qags,
     *  26,ier,lvl)
      return
      end
      subroutine qagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,
     *   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, general-purpose,
c             (end point) singularities, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c
c***references  (none)
c***routines called  qelg,qk21,qpsrt,r1mach
c***end prologue  qagse
c
      real a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,correc,defabs,defab1,defab2,r1mach,
     *  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,
     *  errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs,
     *  reseps,result,res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     * res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine qelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table
c                       which is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered
c                       up to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation
c                       i.e. before subdividing the smallest interval
c                       we try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  qagse
      epmach = r1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      if(epsabs.le.0.0e+00.and.epsrel.lt.amax1(0.5e+02*epmach,0.5e-14))
     *   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = r1mach(1)
      oflow = r1mach(2)
      ierro = 0
      call qk21(f,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = abs(result)
      errbnd = amax1(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0e+02*epmach*defabs.and.abserr.gt.
     *  errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0e+00) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
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
      if(dres.ge.(0.1e+01-0.5e+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk21(f,a1,b1,area1,error1,resabs,defab1)
        call qk21(f,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(abs(rlist(maxerr)-area12).gt.0.1e-04*abs(area12)
     *  .or.erro12.lt.0.99e+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)*
     *  (abs(a2)+0.1e+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   30   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors
c           over the larger intervals (erlarg) and perform
c           extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1e-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = amax1(epsabs,epsrel*abs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5e+00
        erlarg = errsum
        go to 90
   80   small = abs(b-a)*0.375e+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0e+00.and.area.ne.0.0e+00) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0e+00) go to 130
      go to 110
  105 if(abserr/abs(result).gt.errsum/abs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.amax1(abs(result),abs(area)).le.
     * defabs*0.1e-01) go to 130
      if(0.1e-01.gt.(result/area).or.(result/area).gt.0.1e+03
     * .or.errsum.gt.abs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0e+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
      subroutine qawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)
c***begin prologue  qawc
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,j4
c***keywords  automatic integrator, special-purpose,
c             cauchy principal value,
c             clenshaw-curtis, globally adaptive
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a
c            cauchy principal value i = integral of f*w over (a,b)
c            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
c            following claim for accuracy
c            abs(i-result).le.max(epsabe,epsrel*abs(i)).
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        real version
c
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     under limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            c      - parameter in the weight function, c.ne.a, c.ne.b.
c                     if c = a or c = b, the routine will end with
c                     ier = 6 .
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate or the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. exept when lenw or limit is invalid,
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c           lenw   - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)), ... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c***references  (none)
c***routines called  qawce,xerror
c***end prologue  qawc
c
      real a,abserr,b,c,epsabs,epsrel,f,result,work
      integer ier,iwork,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qawc
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qawce.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      call qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,
     *  work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qawc,
     *  26,ier,lvl)
      return
      end
      subroutine qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,
     *   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qawce
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,j4
c***keywords  automatic integrator, special-purpose,
c             cauchy principal value, clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***  purpose  the routine calculates an approximation result to a
c              cauchy principal value i = integral of f*w over (a,b)
c              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
c              following claim for accuracy
c              abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            c      - real
c                     parameter in the weight function, c.ne.a, c.ne.b
c                     if c = a or c = b, the routine will end with
c                     ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement it is advised to analyze the
c                             the integrand, in order to determine the
c                             the integration difficulties. if the
c                             position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some interior points of
c                             the integration interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the integral
c                      approximations on the subintervals
c
c            elist   - real
c                      vector of dimension limit, the first  last
c                      elements of which are the moduli of the absolute
c                      error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the error
c                      estimates over the subintervals, so that
c                      elist(iord(1)), ..., elist(iord(k)) with k = last
c                      if last.le.(limit/2+2), and k = limit+1-last
c                      otherwise, form a decreasing sequence
c
c            last    - integer
c                      number of subintervals actually produced in
c                      the subdivision process
c
c***references  (none)
c***routines called  qc25c,qpsrt,r1mach
c***end prologue  qawce
c
      real a,aa,abserr,alist,area,area1,area12,area2,a1,a2,b,bb,blist,
     *  b1,b2,c,r1mach,elist,epmach,epsabs,epsrel,errbnd,errmax,error1,
     *  error2,errsum,f,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,krule,last,limit,maxerr,nev,
     *  neval,nrmax
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qawce
      epmach = r1mach(4)
      uflow = r1mach(1)
c
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(c.eq.a.or.c.eq.b.or.(epsabs.le.0.0e+00.and
     *  .epsrel.lt.amax1(0.5e+02*epmach,0.5e-14))) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      aa=a
      bb=b
      if (a.le.b) go to 10
      aa=b
      bb=a
10    ier=0
      krule = 1
      call qc25c(f,aa,bb,c,result,abserr,krule,neval)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      alist(1) = a
      blist(1) = b
c
c           test on accuracy
c
      errbnd = amax1(epsabs,epsrel*abs(result))
      if(limit.eq.1) ier = 1
      if(abserr.lt.amin1(0.1e-01*abs(result),errbnd)
     *  .or.ier.eq.1) go to 70
c
c           initialization
c           --------------
c
      alist(1) = aa
      blist(1) = bb
      rlist(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 40 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        b2 = blist(maxerr)
        if(c.le.b1.and.c.gt.a1) b1 = 0.5e+00*(c+b2)
        if(c.gt.b1.and.c.lt.b2) b1 = 0.5e+00*(a1+c)
        a2 = b1
        krule = 2
        call qc25c(f,a1,b1,c,area1,error1,krule,nev)
        neval = neval+nev
        call qc25c(f,a2,b2,c,area2,error2,krule,nev)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(abs(rlist(maxerr)-area12).lt.0.1e-04*abs(area12)
     *    .and.erro12.ge.0.99e+00*errmax.and.krule.eq.0)
     *    iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax.and.krule.eq.0)
     *    iroff2 = iroff2+1
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 15
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1.ge.6.and.iroff2.gt.20) ier = 2
c
c           set error flag in the case that number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)
     *    *(abs(a2)+0.1e+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   15   if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   30    call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 50
   40 continue
c
c           compute final result.
c           ---------------------
c
   50 result = 0.0e+00
      do 60 k=1,last
        result = result+rlist(k)
   60 continue
      abserr = errsum
   70 if (aa.eq.b) result=-result
  999 return
      end
      subroutine qawf(f,a,omega,integr,epsabs,result,abserr,neval,ier,
     *   limlst,lst,leniw,maxp1,lenw,iwork,work)
c***begin prologue  qawf
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,fourier
c             integral, integration between zeros with dqawoe,
c             convergence acceleration with dqext
c***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            fourier integral
c            i = integral of f(x)*w(x) over (a,infinity)
c            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        real version
c
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            omega  - real
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested, epsabs.gt.0.
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of cycles allowed
c                             has been achieved, i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account). examine the array iwork which
c                             contains the error flags on the cycles, in
c                             order to look for eventual local
c                             integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence accelaration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy.
c                             as in the case of ier = 1, it is advised
c                             to examine the array iwork which contains
c                             the error flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.1 or
c                              leniw.lt.(limlst+2) or maxp1.lt.1 or
c                              lenw.lt.(leniw*2+maxp1*25).
c                              result, abserr, neval, lst are set to
c                              zero.
c                         = 7 bad integrand behaviour occurs within
c                             one or more of the cycles. location and
c                             type of the difficulty involved can be
c                             determined from the first lst elements of
c                             vector iwork.  here lst is the number of
c                             cycles actually needed (see below).
c                             iwork(k) = 1 the maximum number of
c                                          subdivisions (=(leniw-limlst)
c                                          /2) has been achieved on the
c                                          k th cycle.
c                                      = 2 occurrence of roundoff error
c                                          is detected and prevents the
c                                          tolerance imposed on the k th
c                                          cycle, from being achieved
c                                          on this cycle.
c                                      = 3 extremely bad integrand
c                                          behaviour occurs at some
c                                          points of the k th cycle.
c                                      = 4 the integration procedure
c                                          over the k th cycle does
c                                          not converge (to within the
c                                          required accuracy) due to
c                                          roundoff in the extrapolation
c                                          procedure invoked on this
c                                          cycle. it is assumed that the
c                                          result on this interval is
c                                          the best which can be
c                                          obtained.
c                                      = 5 the integral over the k th
c                                          cycle is probably divergent
c                                          or slowly convergent. it must
c                                          be noted that divergence can
c                                          occur with any other value of
c                                          iwork(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie,
c                    and ier = iwork(1) (with meaning as described
c                    for iwork(k),k = 1).
c
c         dimensioning parameters
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.3.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            lst    - integer
c                     on return, lst indicates the number of cycles
c                     actually needed for the integration.
c                     if omega = 0, then lst is set to 1.
c
c            leniw  - integer
c                     dimensioning parameter for iwork. on entry,
c                     (leniw-limlst)/2 equals the maximum number of
c                     subintervals allowed in the partition of each
c                     cycle, leniw.ge.(limlst+2).
c                     if leniw.lt.(limlst+2), the routine will end with
c                     ier = 6.
c
c            maxp1  - integer
c                     maxp1 gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e. for
c                     the intervals of lengths abs(b-a)*2**(-l),
c                     l = 0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, iwork(k) for k = 1, 2, ..., lst
c                     contain the error flags on the cycles.
c
c            work   - real
c                     vector of dimension at least
c                     on return,
c                     work(1), ..., work(lst) contain the integral
c                      approximations over the cycles,
c                     work(limlst+1), ..., work(limlst+lst) contain
c                      the error extimates over the cycles.
c                     further elements of work have no specific
c                     meaning for the user.
c
c***references  (none)
c***routines called  qawfe,xerror
c***end prologue  qawf
c
       real a,abserr,epsabs,f,omega,result,work
       integer ier,integr,leniw,limit,limlst,lvl,lst,l1,l2,l3,l4,l5,l6,
     *  maxp1,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of limlst, leniw, maxp1 and lenw.
c
c***first executable statement  qawf
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limlst.lt.3.or.leniw.lt.(limlst+2).or.maxp1.lt.1.or.lenw.lt.
     *   (leniw*2+maxp1*25)) go to 10
c
c         prepare call for qawfe
c
      limit = (leniw-limlst)/2
      l1 = limlst+1
      l2 = limlst+l1
      l3 = limit+l2
      l4 = limit+l3
      l5 = limit+l4
      l6 = limit+l5
      ll2 = limit+l1
      call qawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,result,
     *  abserr,neval,ier,work(1),work(l1),iwork(1),lst,work(l2),
     *  work(l3),work(l4),work(l5),iwork(l1),iwork(ll2),work(l6))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qawf,
     *  26,ier,lvl)
      return
      end
      subroutine qawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,
     *   result,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)
c***begin prologue  qawfe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1
c***keywords  automatic integrator, special-purpose,
c             fourier integrals,
c             integration between zeros with dqawoe,
c             convergence acceleration with dqelg
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a
c            given fourier integal
c            i = integral of f(x)*w(x) over (a,infinity)
c             where w(x) = cos(omega*x) or w(x) = sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to
c                     be declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            omega  - real
c                     parameter in the weight function
c
c            integr - integer
c                     indicates which weight function is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested, epsabs.gt.0
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.1.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     allowed in the partition of each cycle, limit.ge.1
c                     each cycle, limit.ge.1.
c
c            maxp1  - integer
c                     gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e.
c                     for the intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c
c         on return
c            result - real
c                     approximation to the integral x
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - ier = 0 normal and reliable termination of
c                             the routine. it is assumed that the
c                             requested accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine. the
c                             estimates for integral and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of  cycles  allowed
c                             has been achieved., i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account).
c                             examine the array iwork which contains
c                             the error flags on the cycles, in order to
c                             look for eventual local integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling appropriate integrators on
c                             the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence acceleration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy. as in the case of
c                             ier = 1, it is advised to examine the
c                             array iwork which contains the error
c                             flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.3.
c                              result, abserr, neval, lst are set
c                              to zero.
c                         = 7 bad integrand behaviour occurs within one
c                             or more of the cycles. location and type
c                             of the difficulty involved can be
c                             determined from the vector ierlst. here
c                             lst is the number of cycles actually
c                             needed (see below).
c                             ierlst(k) = 1 the maximum number of
c                                           subdivisions (= limit) has
c                                           been achieved on the k th
c                                           cycle.
c                                       = 2 occurrence of roundoff error
c                                           is detected and prevents the
c                                           tolerance imposed on the
c                                           k th cycle, from being
c                                           achieved.
c                                       = 3 extremely bad integrand
c                                           behaviour occurs at some
c                                           points of the k th cycle.
c                                       = 4 the integration procedure
c                                           over the k th cycle does
c                                           not converge (to within the
c                                           required accuracy) due to
c                                           roundoff in the
c                                           extrapolation procedure
c                                           invoked on this cycle. it
c                                           is assumed that the result
c                                           on this interval is the
c                                           best which can be obtained.
c                                       = 5 the integral over the k th
c                                           cycle is probably divergent
c                                           or slowly convergent. it
c                                           must be noted that
c                                           divergence can occur with
c                                           any other value of
c                                           ierlst(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie
c                    and ier = ierlst(1) (with meaning as described
c                    for ierlst(k), k = 1).
c
c            rslst  - real
c                     vector of dimension at least limlst
c                     rslst(k) contains the integral contribution
c                     over the interval (a+(k-1)c,a+kc) where
c                     c = (2*int(abs(omega))+1)*pi/abs(omega),
c                     k = 1, 2, ..., lst.
c                     note that, if omega = 0, rslst(1) contains
c                     the value of the integral over (a,infinity).
c
c            erlst  - real
c                     vector of dimension at least limlst
c                     erlst(k) contains the error estimate corresponding
c                     with rslst(k).
c
c            ierlst - integer
c                     vector of dimension at least limlst
c                     ierlst(k) contains the error flag corresponding
c                     with rslst(k). for the meaning of the local error
c                     flags see description of output parameter ier.
c
c            lst    - integer
c                     number of subintervals needed for the integration
c                     if omega = 0 then lst is set to 1.
c
c            alist, blist, rlist, elist - real
c                     vector of dimension at least limit,
c
c            iord, nnlog - integer
c                     vector of dimension at least limit, providing
c                     space for the quantities needed in the subdivision
c                     process of each cycle
c
c            chebmo - real
c                     array of dimension at least (maxp1,25), providing
c                     space for the chebyshev moments needed within the
c                     cycles
c
c***references  (none)
c***routines called  qagie,qawoe,qelg,r1mach
c***end prologue  qawfe
c
      real a,abseps,abserr,alist,blist,chebmo,correc,cycle,
     *  c1,c2,dl,dla,drl,elist,ep,eps,epsa,epsabs,erlst,
     *  errsum,fact,omega,p,pi,p1,psum,reseps,result,res3la,rlist,rslst
     *  ,r1mach,uflow
      integer ier,ierlst,integr,iord,ktmin,l,lst,limit,ll,maxp1,
     *    nev,neval,nnlog,nres,numrl2
c
      dimension alist(limit),blist(limit),chebmo(maxp1,25),elist(limit),
     *  erlst(limlst),ierlst(limlst),iord(limit),nnlog(limit),psum(52),
     *  res3la(3),rlist(limit),rslst(limlst)
c
      external f
c
c
c            the dimension of  psum  is determined by the value of
c            limexp in subroutine qelg (psum must be
c            of dimension (limexp+2) at least).
c
c           list of major variables
c           -----------------------
c
c           c1, c2    - end points of subinterval (of length
c                       cycle)
c           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
c           psum      - vector of dimension at least (limexp+2)
c                       (see routine qelg)
c                       psum contains the part of the epsilon
c                       table which is still needed for further
c                       computations.
c                       each element of psum is a partial sum of
c                       the series which should sum to the value of
c                       the integral.
c           errsum    - sum of error estimates over the
c                       subintervals, calculated cumulatively
c           epsa      - absolute tolerance requested over current
c                       subinterval
c           chebmo    - array containing the modified chebyshev
c                       moments (see also routine qc25f)
c
      data p/0.9e+00/,pi/0.31415926535897932e+01/
c
c           test on validity of parameters
c           ------------------------------
c
c***first executable statement  qawfe
      result = 0.0e+00
      abserr = 0.0e+00
      neval = 0
      lst = 0
      ier = 0
      if((integr.ne.1.and.integr.ne.2).or.epsabs.le.0.0e+00.or.
     *  limlst.lt.3) ier = 6
      if(ier.eq.6) go to 999
      if(omega.ne.0.0e+00) go to 10
c
c           integration by qagie if omega is zero
c           --------------------------------------
c
      if(integr.eq.1) call qagie(f,0.0e+00,1,epsabs,0.0e+00,limit,
     *  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      rslst(1) = result
      erlst(1) = abserr
      ierlst(1) = ier
      lst = 1
      go to 999
c
c           initializations
c           ---------------
c
   10 l = abs(omega)
      dl = 2*l+1
      cycle = dl*pi/abs(omega)
      ier = 0
      ktmin = 0
      neval = 0
      numrl2 = 0
      nres = 0
      c1 = a
      c2 = cycle+a
      p1 = 0.1e+01-p
      eps = epsabs
      uflow = r1mach(1)
      if(epsabs.gt.uflow/p1) eps = epsabs*p1
      ep = eps
      fact = 0.1e+01
      correc = 0.0e+00
      abserr = 0.0e+00
      errsum = 0.0e+00
c
c           main do-loop
c           ------------
c
      do 50 lst = 1,limlst
c
c           integrate over current subinterval.
c
        dla = lst
        epsa = eps*fact
        call qawoe(f,c1,c2,omega,integr,epsa,0.0e+00,limit,lst,maxp1,
     *  rslst(lst),erlst(lst),nev,ierlst(lst),last,alist,blist,rlist,
     *  elist,iord,nnlog,momcom,chebmo)
        neval = neval+nev
        fact = fact*p
        errsum = errsum+erlst(lst)
        drl = 0.5e+02*abs(rslst(lst))
c
c           test on accuracy with partial sum
c
        if(errsum+drl.le.epsabs.and.lst.ge.6) go to 80
        correc = amax1(correc,erlst(lst))
        if(ierlst(lst).ne.0) eps = amax1(ep,correc*p1)
        if(ierlst(lst).ne.0) ier = 7
        if(ier.eq.7.and.(errsum+drl).le.correc*0.1e+02.and.
     *  lst.gt.5) go to 80
        numrl2 = numrl2+1
        if(lst.gt.1) go to 20
        psum(1) = rslst(1)
        go to 40
   20   psum(numrl2) = psum(ll)+rslst(lst)
        if(lst.eq.2) go to 40
c
c           test on maximum number of subintervals
c
        if(lst.eq.limlst) ier = 1
c
c           perform new extrapolation
c
        call qelg(numrl2,psum,reseps,abseps,res3la,nres)
c
c           test whether extrapolated result is influenced by
c           roundoff
c
        ktmin = ktmin+1
        if(ktmin.ge.15.and.abserr.le.0.1e-02*(errsum+drl)) ier = 4
        if(abseps.gt.abserr.and.lst.ne.3) go to 30
        abserr = abseps
        result = reseps
        ktmin = 0
c
c           if ier is not 0, check whether direct result (partial
c           sum) or extrapolated result yields the best integral
c           approximation
c
        if((abserr+0.1e+02*correc).le.epsabs.or.
     *  (abserr.le.epsabs.and.0.1e+02*correc.ge.epsabs)) go to 60
   30   if(ier.ne.0.and.ier.ne.7) go to 60
   40   ll = numrl2
        c1 = c2
        c2 = c2+cycle
   50 continue
c
c         set final result and error estimate
c         -----------------------------------
c
   60 abserr = abserr+0.1e+02*correc
      if(ier.eq.0) go to 999
      if(result.ne.0.0e+00.and.psum(numrl2).ne.0.0e+00) go to 70
      if(abserr.gt.errsum) go to 80
      if(psum(numrl2).eq.0.0e+00) go to 999
   70 if(abserr/abs(result).gt.(errsum+drl)/abs(psum(numrl2)))
     *  go to 80
      if(ier.ge.1.and.ier.ne.7) abserr = abserr+drl
      go to 999
   80 result = psum(numrl2)
      abserr = errsum+drl
  999 return
      end
      subroutine qawo(f,a,b,omega,integr,epsabs,epsrel,result,abserr,
     *   neval,ier,leniw,maxp1,lenw,last,iwork,work)
c***begin prologue  qawo
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral
c            i = integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x) = sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the function
c                     f(x).  the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            omega  - real
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if epsabs.le.0 and
c                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of  integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             (= leniw/2) has been achieved. one can
c                             allow more subdivisions by increasing the
c                             value of leniw (and taking the according
c                             dimension adjustments into account).
c                             however, if this yields no improvement it
c                             is advised to analyze the integrand in
c                             order to determine the integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling the integrator on the
c                             subranges. if possible, an appropriate
c                             special-purpose integrator should be used
c                             which is designed for handling the type of
c                             difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some interior points of the
c                             integration interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table. it is presumed that
c                             the requested tolerance cannot be achieved
c                             due to roundoff in the extrapolation
c                             table, and that the returned result is
c                             the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2),
c                             or leniw.lt.2 or maxp1.lt.1 or
c                             lenw.lt.leniw*2+maxp1*25.
c                             result, abserr, neval, last are set to
c                             zero. except when leniw, maxp1 or lenw are
c                             invalid, work(limit*2+1), work(limit*3+1),
c                             iwork(1), iwork(limit+1) are set to zero,
c                             work(1) is set to a and work(limit+1) to
c                             b.
c
c         dimensioning parameters
c            leniw  - integer
c                     dimensioning parameter for iwork.
c                     leniw/2 equals the maximum number of subintervals
c                     allowed in the partition of the given integration
c                     interval (a,b), leniw.ge.2.
c                     if leniw.lt.2, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least leniw*2+maxp1*25.
c                     if lenw.lt.(leniw*2+maxp1*25), the routine will
c                     end with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of subintervals
c                     produced in the subdivision process, which
c                     determines the number of significant elements
c                     actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension at least leniw
c                     on return, the first k elements of which contain
c                     pointers to the error estimates over the
c                     subintervals, such that work(limit*3+iwork(1)), ..
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence, with limit = lenw/2 , and k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise.
c                     furthermore, iwork(limit+1), ..., iwork(limit+
c                     last) indicate the subdivision levels of the
c                     subintervals, such that iwork(limit+i) = l means
c                     that the subinterval numbered i is of length
c                     abs(b-a)*2**(1-l).
c
c            work   - real
c                     vector of dimension at least lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last) contain
c                      the integral approximations over the
c                      subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c                     work(limit*4+1), ..., work(limit*4+maxp1*25)
c                      provide space for storing the chebyshev moments.
c                     note that limit = lenw/2.
c
c***references  (none)
c***routines called  qawoe,xerror
c***end prologue  qawo
c
       real a,abserr,b,epsabs,epsrel,f,omega,result
       integer ier,integr,leniw,lvl,l1,l2,l3,l4,maxp1,momcom,neval
c
       dimension iwork(leniw),work(lenw)
c
       external f
c
c         check validity of leniw, maxp1 and lenw.
c
c***first executable statement  qawo
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(leniw.lt.2.or.maxp1.lt.1.or.lenw.lt.(leniw*2+maxp1*25))
     *  go to 10
c
c         prepare call for qawoe
c
      limit = leniw/2
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      l4 = limit+l3
      call qawoe(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result,
     *  abserr,neval,ier,last,work(1),work(l1),work(l2),work(l3),
     *  iwork(1),iwork(l1),momcom,work(l4))
c
c         call error handler if necessary
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qawo,
     *  26,ier,lvl)
      return
      end
      subroutine qawoe (f,a,b,omega,integr,epsabs,epsrel,limit,icall,
     *  maxp1,result,abserr,neval,ier,last,alist,blist,rlist,elist,iord,
     *   nnlog,momcom,chebmo)
c***begin prologue  qawoe
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             integrand with oscillatory cos or sin factor,
c             clenshaw-curtis method, (end point) singularities,
c             extrapolation, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral
c            i = integral of f(x)*w(x) over (a,b)
c            where w(x) = cos(omega*x) or w(x) = sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of oscillatory integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            omega  - real
c                     parameter in the integrand weight function
c
c            integr - integer
c                     indicates which of the weight functions is to be
c                     used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1 and integr.ne.2, the routine
c                     will end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subdivisions
c                     in the partition of (a,b), limit.ge.1.
c
c            icall  - integer
c                     if dqawoe is to be used only once, icall must
c                     be set to 1.  assume that during this call, the
c                     chebyshev moments (for clenshaw-curtis integration
c                     of degree 24) have been computed for intervals of
c                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
c                     if icall.gt.1 this means that dqawoe has been
c                     called twice or more on intervals of the same
c                     length abs(b-a). the chebyshev moments already
c                     computed are then re-used in subsequent calls.
c                     if icall.lt.1, the routine will end with ier = 6.
c
c            maxp1  - integer
c                     gives an upper bound on the number of chebyshev
c                     moments which can be stored, i.e. for the
c                     intervals of lenghts abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1.
c                     if maxp1.lt.1, the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the
c                             requested accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine.
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved due to
c                             roundoff in the extrapolation table,
c                             and that the returned result is the
c                             best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.gt.0.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or (integr.ne.1 and integr.ne.2) or
c                             icall.lt.1 or maxp1.lt.1.
c                             result, abserr, neval, last, rlist(1),
c                             elist(1), iord(1) and nnlog(1) are set
c                             to zero. alist(1) and blist(1) are set
c                             to a and b respectively.
c
c            last  -  integer
c                     on return, last equals the number of
c                     subintervals produces in the subdivision
c                     process, which determines the number of
c                     significant elements actually in the
c                     work arrays.
c            alist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the error
c                     estimates over the subintervals,
c                     such that elist(iord(1)), ...,
c                     elist(iord(k)) form a decreasing sequence, with
c                     k = last if last.le.(limit/2+2), and
c                     k = limit+1-last otherwise.
c
c            nnlog  - integer
c                     vector of dimension at least limit, containing the
c                     subdivision levels of the subintervals, i.e.
c                     iwork(i) = l means that the subinterval
c                     numbered i is of length abs(b-a)*2**(1-l)
c
c         on entry and return
c            momcom - integer
c                     indicating that the chebyshev moments
c                     have been computed for intervals of lengths
c                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
c                     momcom.lt.maxp1
c
c            chebmo - real
c                     array of dimension (maxp1,25) containing the
c                     chebyshev moments
c
c***references  (none)
c***routines called  qc25f,qelg,qpsrt,r1mach
c***end prologue  qawoe
c
      real a,abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,b,blist,b1,b2,chebmo,correc,defab1,defab2,defabs,
     *  domega,r1mach,dres,elist,epmach,epsabs,epsrel,erlarg,
     *  erlast,errbnd,errmax,error1,erro12,error2,errsum,ertest,f,oflow,
     *  omega,resabs,reseps,result,res3la,rlist,rlist2,small,uflow,width
      integer icall,id,ier,ierro,integr,iord,iroff1,iroff2,iroff3,
     *  jupbnd,k,ksgn,ktmin,last,limit,maxerr,maxp1,momcom,nev,
     *  neval,nnlog,nres,nrmax,nrmom,numrl2
      logical extrap,noext,extall
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),rlist2(52),res3la(3),chebmo(maxp1,25),nnlog(limit)
c
      external f
c
c            the dimension of rlist2 is determined by  the value of
c            limexp in subroutine qelg (rlist2 should be of
c            dimension (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2
c                       containing the part of the epsilon table
c                       which is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements in rlist2. if an appropriate
c                       approximation to the compounded integral has
c                       been obtained it is put in rlist2(numrl2) after
c                       numrl2 has been increased by one
c           small     - length of the smallest interval considered
c                       up to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation, i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  qawoe
      epmach = r1mach(4)
c
c         test on validity of parameters
c         ------------------------------
c
      ier = 0
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      nnlog(1) = 0
      if((integr.ne.1.and.integr.ne.2).or.(epsabs.le.0.0e+00.and.
     *  epsrel.lt.amax1(0.5e+02*epmach,0.5e-14)).or.icall.lt.1.or.
     *  maxp1.lt.1) ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      domega = abs(omega)
      nrmom = 0
      if (icall.gt.1) go to 5
      momcom = 0
    5 call qc25f(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr,
     *  neval,defabs,resabs,momcom,chebmo)
c
c           test on accuracy.
c
      dres = abs(result)
      errbnd = amax1(epsabs,epsrel*dres)
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.0.1e+03*epmach*defabs.and.abserr.gt.
     *  errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.abserr.le.errbnd) go to 200
c
c           initializations
c           ---------------
c
      uflow = r1mach(1)
      oflow = r1mach(2)
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ktmin = 0
      small = abs(b-a)*0.75e+00
      nres = 0
      numrl2 = 0
      extall = .false.
      if(0.5e+00*abs(b-a)*domega.gt.0.2e+01) go to 10
      numrl2 = 1
      extall = .true.
      rlist2(1) = result
   10 if(0.25e+00*abs(b-a)*domega.le.0.2e+01) extall = .true.
      ksgn = -1
      if(dres.ge.(0.1e+01-0.5e+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 140 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest
c           error estimate.
c
        nrmom = nnlog(maxerr)+1
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qc25f(f,a1,b1,domega,integr,nrmom,maxp1,0,
     *  area1,error1,nev,resabs,defab1,momcom,chebmo)
        neval = neval+nev
        call qc25f(f,a2,b2,domega,integr,nrmom,maxp1,1,
     *  area2,error2,nev,resabs,defab2,momcom,chebmo)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 25
        if(abs(rlist(maxerr)-area12).gt.0.1e-04*abs(area12)
     *  .or.erro12.lt.0.99e+00*errmax) go to 20
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   20   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   25   rlist(maxerr) = area1
        rlist(last) = area2
        nnlog(maxerr) = nrmom
        nnlog(last) = nrmom
        errbnd = amax1(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually
c           set error flag
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)
     *  *(abs(a2)+0.1e+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 30
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 40
   30   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   40   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
      if(errsum.le.errbnd) go to 170
      if(ier.ne.0) go to 150
        if(last.eq.2.and.extall) go to 120
        if(noext) go to 140
        if(.not.extall) go to 50
        erlarg = erlarg-erlast
        if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 70
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
   50   width = abs(blist(maxerr)-alist(maxerr))
        if(width.gt.small) go to 140
        if(extall) go to 60
c
c           test whether we can start with the extrapolation
c           procedure (we do this if we integrate over the
c           next interval with use of a gauss-kronrod rule - see
c           subroutine qc25f).
c
        small = small*0.5e+00
        if(0.25e+00*width*domega.gt.0.2e+01) go to 140
        extall = .true.
        go to 130
   60   extrap = .true.
        nrmax = 2
   70   if(ierro.eq.3.or.erlarg.le.ertest) go to 90
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors
c           over the larger intervals (erlarg) and perform
c           extrapolation.
c
        jupbnd = last
        if (last.gt.(limit/2+2)) jupbnd = limit+3-last
        id = nrmax
        do 80 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 140
          nrmax = nrmax+1
   80   continue
c
c           perform extrapolation.
c
   90   numrl2 = numrl2+1
        rlist2(numrl2) = area
        if(numrl2.lt.3) go to 110
        call qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1e-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 100
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = amax1(epsabs,epsrel*abs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 150
c
c           prepare bisection of the smallest interval.
c
  100   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 150
  110   maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5e+00
        erlarg = errsum
        go to 140
  120   small = small*0.5e+00
        numrl2 = numrl2+1
        rlist2(numrl2) = area
  130   ertest = errbnd
        erlarg = errsum
  140 continue
c
c           set the final result.
c           ---------------------
c
  150 if(abserr.eq.oflow.or.nres.eq.0) go to 170
      if(ier+ierro.eq.0) go to 165
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0e+00.and.area.ne.0.0e+00) go to 160
      if(abserr.gt.errsum) go to 170
      if(area.eq.0.0e+00) go to 190
      go to 165
  160 if(abserr/abs(result).gt.errsum/abs(area)) go to 170
c
c           test on divergence.
c
  165 if(ksgn.eq.(-1).and.amax1(abs(result),abs(area)).le.
     * defabs*0.1e-01) go to 190
      if(0.1e-01.gt.(result/area).or.(result/area).gt.0.1e+03
     * .or.errsum.ge.abs(area)) ier = 6
      go to 190
c
c           compute global integral sum.
c
  170 result = 0.0e+00
      do 180 k=1,last
        result = result+rlist(k)
  180 continue
      abserr = errsum
  190 if (ier.gt.2) ier=ier-1
  200 if (integr.eq.2.and.omega.lt.0.0e+00) result=-result
  999 return
      end
      subroutine qaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,
     *   abserr,neval,ier,limit,lenw,last,iwork,work)
c***begin prologue  qaws
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end-point singularities,
c             clenshaw-curtis, globally adaptive
c***author  piessens,robert,appl. math. & progr. div. -k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - real
c                     parameter in the integrand function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - real
c                     parameter in the integrand function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
c                     will end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand, in order to
c                             determine the integration difficulties
c                             which prevent the requested tolerance from
c                             being achieved. in case of a jump
c                             discontinuity or a local singularity
c                             of algebraico-logarithmic type at one or
c                             more interior points of the integration
c                             range, one should proceed by splitting up
c                             the interval at these points and calling
c                             the integrator on the subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1) or
c                             or integr.lt.1 or integr.gt.4 or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.2 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. except when lenw or limit is invalid
c                             iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit  - integer
c                     dimensioning parameter for iwork
c                     limit determines the maximum number of
c                     subintervals in the partition of the given
c                     integration interval (a,b), limit.ge.2.
c                     if limit.lt.2, the routine will end with ier = 6.
c
c            lenw   - integer
c                     dimensioning parameter for work
c                     lenw must be at least limit*4.
c                     if lenw.lt.limit*4, the routine will end
c                     with ier = 6.
c
c            last   - integer
c                     on return, last equals the number of
c                     subintervals produced in the subdivision process,
c                     which determines the significant number of
c                     elements actually in the work arrays.
c
c         work arrays
c            iwork  - integer
c                     vector of dimension limit, the first k
c                     elements of which contain pointers
c                     to the error estimates over the subintervals,
c                     such that work(limit*3+iwork(1)), ...,
c                     work(limit*3+iwork(k)) form a decreasing
c                     sequence with k = last if last.le.(limit/2+2),
c                     and k = limit+1-last otherwise
c
c            work   - real
c                     vector of dimension lenw
c                     on return
c                     work(1), ..., work(last) contain the left
c                      end points of the subintervals in the
c                      partition of (a,b),
c                     work(limit+1), ..., work(limit+last) contain
c                      the right end points,
c                     work(limit*2+1), ..., work(limit*2+last)
c                      contain the integral approximations over
c                      the subintervals,
c                     work(limit*3+1), ..., work(limit*3+last)
c                      contain the error estimates.
c
c***references  (none)
c***routines called  qawse,xerror
c***end prologue  qaws
c
      real a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
      integer ier,integr,iwork,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qaws
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.2.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qawse.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call qawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result,
     *  abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qaws,
     *  26,ier,lvl)
      return
      end
      subroutine qawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,
     *   result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qawse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end point singularities,
c             clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points,
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - real
c                     parameter in the weight function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - real
c                     parameter in the weight function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
c                     will end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.2
c                     if limit.lt.2, the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement, it is advised to analyze the
c                             integrand in order to determine the
c                             integration difficulties which prevent the
c                             requested tolerance from being achieved.
c                             in case of a jump discontinuity or a local
c                             singularity of algebraico-logarithmic type
c                             at one or more interior points of the
c                             integration range, one should proceed by
c                             splitting up the interval at these
c                             points and calling the integrator on the
c                             subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1), or
c                             integr.lt.1 or integr.gt.4, or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             or limit.lt.2.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - real
c                     vector of dimension at least limit,the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     of which are pointers to the error
c                     estimates over the subintervals, so that
c                     elist(iord(1)), ..., elist(iord(k)) with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise form a decreasing sequence
c
c            last   - integer
c                     number of subintervals actually produced in
c                     the subdivision process
c
c***references  (none)
c***routines called  qc25s,qmomo,qpsrt,r1mach
c***end prologue  qawse
c
      real a,abserr,alfa,alist,area,area1,area12,
     *  area2,a1,a2,b,beta,blist,b1,b2,centre,
     *  r1mach,elist,epmach,epsabs,epsrel,errbnd,errmax,
     *  error1,erro12,error2,errsum,f,resas1,resas2,result,rg,rh,ri,rj,
     *  rlist,uflow
      integer ier,int,integr,iord,iroff1,iroff2,k,last,
     *  limit,maxerr,nev,neval,nrmax
c
      external f
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),ri(25),rj(25),rh(25),rg(25)
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qawse
      epmach = r1mach(4)
      uflow = r1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(b.le.a.or.(epsabs.eq.0.0e+00.and.
     *  epsrel.lt.amax1(0.5e+02*epmach,0.5e-14)).or.alfa.le.(-0.1e+01).
     *  or.beta.le.(-0.1e+01).or.integr.lt.1.or.integr.gt.4.or.
     *  limit.lt.2) go to 999
      ier = 0
c
c           compute the modified chebyshev moments.
c
      call qmomo(alfa,beta,ri,rj,rg,rh,integr)
c
c           integrate over the intervals (a,(a+b)/2)
c           and ((a+b)/2,b).
c
      centre = 0.5e+00*(b+a)
      call qc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
      neval = nev
      call qc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
      last = 2
      neval = neval+nev
      result = area1+area2
      abserr = error1+error2
c
c           test on accuracy.
c
      errbnd = amax1(epsabs,epsrel*abs(result))
c
c           initialization
c           --------------
c
      if(error2.gt.error1) go to 10
      alist(1) = a
      alist(2) = centre
      blist(1) = centre
      blist(2) = b
      rlist(1) = area1
      rlist(2) = area2
      elist(1) = error1
      elist(2) = error2
      go to 20
   10 alist(1) = centre
      alist(2) = a
      blist(1) = b
      blist(2) = centre
      rlist(1) = area2
      rlist(2) = area1
      elist(1) = error2
      elist(2) = error1
   20 iord(1) = 1
      iord(2) = 2
      if(limit.eq.2) ier = 1
      if(abserr.le.errbnd.or.ier.eq.1) go to 999
      errmax = elist(1)
      maxerr = 1
      nrmax = 1
      area = result
      errsum = abserr
      iroff1 = 0
      iroff2 = 0
c
c            main do-loop
c            ------------
c
      do 60 last = 3,limit
c
c           bisect the subinterval with largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
c
        call qc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
        neval = neval+nev
        call qc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
        neval = neval+nev
c
c           improve previous approximations integral and error
c           and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(a.eq.a1.or.b.eq.b2) go to 30
        if(resas1.eq.error1.or.resas2.eq.error2) go to 30
c
c           test for roundoff error.
c
        if(abs(rlist(maxerr)-area12).lt.0.1e-04*abs(area12)
     *  .and.erro12.ge.0.99e+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
   30   rlist(maxerr) = area1
        rlist(last) = area2
c
c           test on accuracy.
c
        errbnd = amax1(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 35
c
c           set error flag in the case that the number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c
c           set error flag in the case of roundoff error.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case of bad integrand behaviour
c           at interior points of integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)*
     *  (abs(a2)+0.1e+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   35   if(error2.gt.error1) go to 40
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 50
   40   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with largest error estimate (to be
c           bisected next).
c
   50   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if (ier.ne.0.or.errsum.le.errbnd) go to 70
   60 continue
c
c           compute final result.
c           ---------------------
c
   70 result = 0.0e+00
      do 80 k=1,last
        result = result+rlist(k)
   80 continue
      abserr = errsum
  999 return
      end
      subroutine qc25c(f,a,b,c,result,abserr,krul,neval)
c***begin prologue  qc25c
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2,j4
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b) with
c            error estimate, where w(x) = 1/(x-c)
c***description
c
c        integration rules for the computation of cauchy
c        principal value integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c           f      - real
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - real
c                    left end point of the integration interval
c
c           b      - real
c                    right end point of the integration interval, b.gt.a
c
c           c      - real
c                    parameter in the weight function
c
c           result - real
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if c lies within ten percent
c                    of the integration interval. in the other case the
c                    15-point kronrod rule obtained by optimal addition
c                    of abscissae to the 7-point gauss rule, is applied.
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           krul   - integer
c                    key which is decreased by 1 if the 15-point
c                    gauss-kronrod scheme has been used
c
c           neval  - integer
c                    number of integrand evaluations
c
c***references  (none)
c***routines called  qcheb,qk15w,qwgtc
c***end prologue  qc25c
c
      real a,abserr,ak22,amom0,amom1,amom2,b,c,cc,
     *  centr,cheb12,cheb24,qwgtc,f,fval,hlgth,p2,p3,p4,
     *  resabs,resasc,result,res12,res24,u,x
      integer i,isym,k,kp,krul,neval
c
      dimension x(11),fval(25),cheb12(13),cheb24(25)
c
      external f,qwgtc
c
c           the vector x contains the values cos(k*pi/24),
c           k = 1, ..., 11, to be used for the chebyshev series
c           expansion of f
c
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),
     *  x(11)/
     *     0.9914448613738104e+00,     0.9659258262890683e+00,
     *     0.9238795325112868e+00,     0.8660254037844386e+00,
     *     0.7933533402912352e+00,     0.7071067811865475e+00,
     *     0.6087614290087206e+00,     0.5000000000000000e+00,
     *     0.3826834323650898e+00,     0.2588190451025208e+00,
     *     0.1305261922200516e+00/
c
c           list of major variables
c           ----------------------
c           fval   - value of the function f at the points
c                    cos(k*pi/24),  k = 0, ..., 24
c           cheb12 - chebyshev series expansion coefficients,
c                    for the function f, of degree 12
c           cheb24 - chebyshev series expansion coefficients,
c                    for the function f, of degree 24
c           res12  - approximation to the integral corresponding
c                    to the use of cheb12
c           res24  - approximation to the integral corresponding
c                    to the use of cheb24
c           qwgtc - external function subprogram defining
c                    the weight function
c           hlgth  - half-length of the interval
c           centr  - mid point of the interval
c
c
c           check the position of c.
c
c***first executable statement  qc25c
      cc = (0.2e+01*c-b-a)/(b-a)
      if(abs(cc).lt.0.11e+01) go to 10
c
c           apply the 15-point gauss-kronrod scheme.
c
      krul = krul-1
      call qk15w(f,qwgtc,c,p2,p3,p4,kp,a,b,result,abserr,
     *  resabs,resasc)
      neval = 15
      if (resasc.eq.abserr) krul = krul+1
      go to 50
c
c           use the generalized clenshaw-curtis method.
c
   10 hlgth = 0.5e+00*(b-a)
      centr = 0.5e+00*(b+a)
      neval = 25
      fval(1) = 0.5e+00*f(hlgth+centr)
      fval(13) = f(centr)
      fval(25) = 0.5e+00*f(centr-hlgth)
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)
        fval(isym) = f(centr-u)
   20 continue
c
c           compute the chebyshev series expansion.
c
      call qcheb(x,fval,cheb12,cheb24)
c
c           the modified chebyshev moments are computed
c           by forward recursion, using amom0 and amom1
c           as starting values.
c
      amom0 = alog(abs((0.1e+01-cc)/(0.1e+01+cc)))
      amom1 = 0.2e+01+cc*amom0
      res12 = cheb12(1)*amom0+cheb12(2)*amom1
      res24 = cheb24(1)*amom0+cheb24(2)*amom1
      do 30 k=3,13
        amom2 = 0.2e+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4e+01/(ak22-0.1e+01)
        res12 = res12+cheb12(k)*amom2
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   30 continue
      do 40 k=14,25
        amom2 = 0.2e+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4e+01/
     *  (ak22-0.1e+01)
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   40 continue
      result = res24
      abserr = abs(res24-res12)
   50 return
      end
      subroutine qc25f(f,a,b,omega,integr,nrmom,maxp1,ksave,result,
     *   abserr,neval,resabs,resasc,momcom,chebmo)
c***begin prologue  qc25f
c***date written   810101   (yymmdd)
c***revision date  211011   (yymmdd)
c***category no.  h2a2a2
c***keywords  integration rules for functions with cos or sin
c             factor, clenshaw-curtis, gauss-kronrod
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute the integral i=integral of f(x) over (a,b)
c            where w(x) = cos(omega*x) or (wx)=sin(omega*x)
c            and to compute j=integral of abs(f) over (a,b). for small
c            value of omega or small intervals (a,b) 15-point gauss-
c            kronrod rule used. otherwise generalized clenshaw-curtis us
c***description
c
c        integration rules for functions with cos or sin factor
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c           f      - real
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to
c                    be declared e x t e r n a l in the calling program.
c
c           a      - real
c                    lower limit of integration
c
c           b      - real
c                    upper limit of integration
c
c           omega  - real
c                    parameter in the weight function
c
c           integr - integer
c                    indicates which weight function is to be used
c                       integr = 1   w(x) = cos(omega*x)
c                       integr = 2   w(x) = sin(omega*x)
c
c           nrmom  - integer
c                    the length of interval (a,b) is equal to the length
c                    of the original integration interval divided by
c                    2**nrmom (we suppose that the routine is used in an
c                    adaptive integration process, otherwise set
c                    nrmom = 0). nrmom must be zero at the first call.
c
c           maxp1  - integer
c                    gives an upper bound on the number of chebyshev
c                    moments which can be stored, i.e. for the
c                    intervals of lengths abs(bb-aa)*2**(-l),
c                    l = 0,1,2, ..., maxp1-2.
c
c           ksave  - integer
c                    key which is one when the moments for the
c                    current interval have been computed
c
c         on return
c           result - real
c                    approximation to the integral i
c
c           abserr - real
c                    estimate of the modulus of the absolute
c                    error, which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           resabs - real
c                    approximation to the integral j
c
c           resasc - real
c                    approximation to the integral of abs(f-i/(b-a))
c
c         on entry and return
c           momcom - integer
c                    for each interval length we need to compute the
c                    chebyshev moments. momcom counts the number of
c                    intervals for which these moments have already been
c                    computed. if nrmom.lt.momcom or ksave = 1, the
c                    chebyshev moments for the interval (a,b) have
c                    already been computed and stored, otherwise we
c                    compute them and we increase momcom.
c
c           chebmo - real
c                    array of dimension at least (maxp1,25) containing
c                    the modified chebyshev moments for the first momcom
c                    momcom interval lengths
c
c***references  (none)
c***routines called  qcheb,qk15w,qwgtf,r1mach,sgtsl
c***end prologue  qc25f
c
      real a,abserr,ac,an,an2,as,asap,ass,b,centr,chebmo,
     *  cheb12,cheb24,conc,cons,cospar,d,qwgtf,
     *  d1,r1mach,d2,estc,ests,f,fval,hlgth,oflow,omega,parint,par2,
     *  par22,p2,p3,p4,resabs,resasc,resc12,resc24,ress12,ress24,
     *  result,sinpar,v,x
      integer i,iers,integr,isym,j,k,ksave,m,maxp1,momcom,neval,
     *  noequ,noeq1,nrmom
c
      dimension chebmo(maxp1,25),cheb12(13),cheb24(25),d(25),d1(25),
     *  d2(25),fval(25),v(28),x(11)
c
      external f,qwgtf
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ...,11, to be used for the chebyshev expansion of f
c
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),
     *  x(10),x(11)/
     *     0.9914448613738104e+00,     0.9659258262890683e+00,
     *     0.9238795325112868e+00,     0.8660254037844386e+00,
     *     0.7933533402912352e+00,     0.7071067811865475e+00,
     *     0.6087614290087206e+00,     0.5000000000000000e+00,
     *     0.3826834323650898e+00,     0.2588190451025208e+00,
     *     0.1305261922200516e+00/
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fval   - value of the function f at the points
c                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5,
c                    k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (a,b)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (a,b)
c           resc12 - approximation to the integral of
c                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
c                    over (-1,+1), using the chebyshev series
c                    expansion of degree 12
c           resc24 - approximation to the same integral, using the
c                    chebyshev series expansion of degree 24
c           ress12 - the analogue of resc12 for the sine
c           ress24 - the analogue of resc24 for the sine
c
c
c           machine dependent constant
c           --------------------------
c
c           oflow is the largest positive magnitude.
c
c***first executable statement  qc25f
      oflow = r1mach(2)
c
      centr = 0.5e+00*(b+a)
      hlgth = 0.5e+00*(b-a)
      parint = omega*hlgth
c
c           compute the integral using the 15-point gauss-kronrod
c           formula if the value of the parameter in the integrand
c           is small.
c
      if(abs(parint).gt.0.2e+01) go to 10
      call qk15w(f,qwgtf,omega,p2,p3,p4,integr,a,b,result,
     *  abserr,resabs,resasc)
      neval = 15
      go to 170
c
c           compute the integral using the generalized clenshaw-
c           curtis method.
c
   10 conc = hlgth*cos(centr*omega)
      cons = hlgth*sin(centr*omega)
      resasc = oflow
      neval = 25
c
c           check whether the chebyshev moments for this interval
c           have already been computed.
c
      if(nrmom.lt.momcom.or.ksave.eq.1) go to 120
c
c           compute a new set of chebyshev moments.
c
      m = momcom+1
      par2 = parint*parint
      par22 = par2+0.2e+01
      sinpar = sin(parint)
      cospar = cos(parint)
c
c           compute the chebyshev moments with respect to cosine.
c
      v(1) = 0.2e+01*sinpar/parint
      v(2) = (0.8e+01*cospar+(par2+par2-0.8e+01)*sinpar/
     *  parint)/par2
      v(3) = (0.32e+02*(par2-0.12e+02)*cospar+(0.2e+01*
     *  ((par2-0.80e+02)*par2+0.192e+03)*sinpar)/
     *  parint)/(par2*par2)
      ac = 0.8e+01*cospar
      as = 0.24e+02*parint*sinpar
      if(abs(parint).gt.0.24e+02) go to 30
c
c           compute the chebyshev moments as the
c           solutions of a boundary value problem with 1
c           initial value (v(3)) and 1 end value (computed
c           using an asymptotic formula).
c
      noequ = 25
      noeq1 = noequ-1
      an = 0.6e+01
      do 20 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2e+01*(an2-0.4e+01)*(par22-an2-an2)
        d2(k) = (an-0.1e+01)*(an-0.2e+01)*par2
        d1(k+1) = (an+0.3e+01)*(an+0.4e+01)*par2
        v(k+3) = as-(an2-0.4e+01)*ac
        an = an+0.2e+01
   20 continue
      an2 = an*an
      d(noequ) = -0.2e+01*(an2-0.4e+01)*(par22-an2-an2)
      v(noequ+3) = as-(an2-0.4e+01)*ac
      v(4) = v(4)-0.56e+02*par2*v(3)
      ass = parint*sinpar
      asap = (((((0.210e+03*par2-0.1e+01)*cospar-(0.105e+03*par2
     *  -0.63e+02)*ass)/an2-(0.1e+01-0.15e+02*par2)*cospar
     *  +0.15e+02*ass)/an2-cospar+0.3e+01*ass)/an2-cospar)/an2
      v(noequ+3) = v(noequ+3)-0.2e+01*asap*par2*(an-0.1e+01)*
     *   (an-0.2e+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
      call sgtsl(noequ,d1,d,d2,v(4),iers)
      go to 50
c
c           compute the chebyshev moments by means of forward
c           recursion.
c
   30 an = 0.4e+01
      do 40 i = 4,13
        an2 = an*an
        v(i) = ((an2-0.4e+01)*(0.2e+01*(par22-an2-an2)*v(i-1)-ac)
     *  +as-par2*(an+0.1e+01)*(an+0.2e+01)*v(i-2))/
     *  (par2*(an-0.1e+01)*(an-0.2e+01))
        an = an+0.2e+01
   40 continue
   50 do 60 j = 1,13
        chebmo(m,2*j-1) = v(j)
   60 continue
c
c           compute the chebyshev moments with respect to sine.
c
      v(1) = 0.2e+01*(sinpar-parint*cospar)/par2
      v(2) = (0.18e+02-0.48e+02/par2)*sinpar/par2
     *  +(-0.2e+01+0.48e+02/par2)*cospar/parint
      ac = -0.24e+02*parint*cospar
      as = -0.8e+01*sinpar
      if(abs(parint).gt.0.24e+02) go to 80
c
c           compute the chebyshev moments as the
c           solutions of a boundary value problem with 1
c           initial value (v(2)) and 1 end value (computed
c           using an asymptotic formula).
c
      an = 0.5e+01
      do 70 k = 1,noeq1
        an2 = an*an
        d(k) = -0.2e+01*(an2-0.4e+01)*(par22-an2-an2)
        d2(k) = (an-0.1e+01)*(an-0.2e+01)*par2
        d1(k+1) = (an+0.3e+01)*(an+0.4e+01)*par2
        v(k+2) = ac+(an2-0.4e+01)*as
        an = an+0.2e+01
   70 continue
      an2 = an*an
      d(noequ) = -0.2e+01*(an2-0.4e+01)*(par22-an2-an2)
      v(noequ+2) = ac+(an2-0.4e+01)*as
      v(3) = v(3)-0.42e+02*par2*v(2)
      ass = parint*cospar
      asap = (((((0.105e+03*par2-0.63e+02)*ass+(0.210e+03*par2
     *  -0.1e+01)*sinpar)/an2+(0.15e+02*par2-0.1e+01)*sinpar-
     *  0.15e+02*ass)/an2-0.3e+01*ass-sinpar)/an2-sinpar)/an2
      v(noequ+2) = v(noequ+2)-0.2e+01*asap*par2*(an-0.1e+01)
     *  *(an-0.2e+01)
c
c           solve the tridiagonal system by means of gaussian
c           elimination with partial pivoting.
c
      call sgtsl(noequ,d1,d,d2,v(3),iers)
      go to 100
c
c           compute the chebyshev moments by means of
c           forward recursion.
c
   80 an = 0.3e+01
      do 90 i = 3,12
        an2 = an*an
        v(i) = ((an2-0.4e+01)*(0.2e+01*(par22-an2-an2)*v(i-1)+as)
     *  +ac-par2*(an+0.1e+01)*(an+0.2e+01)*v(i-2))
     *  /(par2*(an-0.1e+01)*(an-0.2e+01))
        an = an+0.2e+01
   90 continue
  100 do 110 j = 1,12
        chebmo(m,2*j) = v(j)
  110 continue
  120 if (nrmom.lt.momcom) m = nrmom+1
       if (momcom.lt.maxp1-1.and.nrmom.ge.momcom) momcom = momcom+1
c
c           compute the coefficients of the chebyshev expansions
c           of degrees 12 and 24 of the function f.
c
      fval(1) = 0.5e+00*f(centr+hlgth)
      fval(13) = f(centr)
      fval(25) = 0.5e+00*f(centr-hlgth)
      do 130 i = 2,12
        isym = 26-i
        fval(i) = f(hlgth*x(i-1)+centr)
        fval(isym) = f(centr-hlgth*x(i-1))
  130 continue
      call qcheb(x,fval,cheb12,cheb24)
c
c           compute the integral and error estimates.
c
      resc12 = cheb12(13)*chebmo(m,13)
      ress12 = 0.0e+00
      k = 11
      do 140 j = 1,6
        resc12 = resc12+cheb12(k)*chebmo(m,k)
        ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
        k = k-2
  140 continue
      resc24 = cheb24(25)*chebmo(m,25)
      ress24 = 0.0e+00
      resabs = abs(cheb24(25))
      k = 23
      do 150 j = 1,12
        resc24 = resc24+cheb24(k)*chebmo(m,k)
        ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
        resabs = resabs+abs(cheb24(k))+abs(cheb24(k+1))
        k = k-2
  150 continue
      estc = abs(resc24-resc12)
      ests = abs(ress24-ress12)
      resabs = resabs*abs(hlgth)
      if(integr.eq.2) go to 160
      result = conc*resc24-cons*ress24
      abserr = abs(conc*estc)+abs(cons*ests)
      go to 170
  160 result = conc*ress24+cons*resc24
      abserr = abs(conc*ests)+abs(cons*estc)
  170 return
      end
      subroutine qc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result,
     *   abserr,resasc,integr,nev)
c***begin prologue  qc25s
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (bl,br), with error
c            estimate, where the weight function w has a singular
c            behaviour of algebraico-logarithmic type at the points
c            a and/or b. (bl,br) is a part of (a,b).
c***description
c
c        integration rules for integrands having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        real version
c
c        parameters
c           f      - real
c                    function subprogram defining the integrand
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - real
c                    left end point of the original interval
c
c           b      - real
c                    right end point of the original interval, b.gt.a
c
c           bl     - real
c                    lower limit of integration, bl.ge.a
c
c           br     - real
c                    upper limit of integration, br.le.b
c
c           alfa   - real
c                    parameter in the weight function
c
c           beta   - real
c                    parameter in the weight function
c
c           ri,rj,rg,rh - real
c                    modified chebyshev moments for the application
c                    of the generalized clenshaw-curtis
c                    method (computed in subroutine dqmomo)
c
c           result - real
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if b1 = a or br = b.
c                    in all other cases the 15-point kronrod
c                    rule is applied, obtained by optimal addition of
c                    abscissae to the 7-point gauss rule.
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           resasc - real
c                    approximation to the integral of abs(f*w-i/(b-a))
c
c           integr - integer
c                    which determines the weight function
c                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
c                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
c                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
c                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
c                                 log(b-x)
c
c           nev    - integer
c                    number of integrand evaluations
c
c***references  (none)
c***routines called  qcheb,qk15w
c***end prologue  qc25s
c
      real a,abserr,alfa,b,beta,bl,br,centr,cheb12,cheb24,
     *  dc,f,factor,fix,fval,hlgth,resabs,resasc,
     *  result,res12,res24,rg,rh,ri,rj,u,qwgts,x
      integer i,integr,isym,nev
c
      dimension cheb12(13),cheb24(25),fval(25),rg(25),rh(25),ri(25),
     *  rj(25),x(11)
c
      external f,qwgts
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ..., 11, to be used for the computation of the
c           chebyshev series expansion of f.
c
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),
     *  x(11)/
     *     0.9914448613738104e+00,     0.9659258262890683e+00,
     *     0.9238795325112868e+00,     0.8660254037844386e+00,
     *     0.7933533402912352e+00,     0.7071067811865475e+00,
     *     0.6087614290087206e+00,     0.5000000000000000e+00,
     *     0.3826834323650898e+00,     0.2588190451025208e+00,
     *     0.1305261922200516e+00/
c
c           list of major variables
c           -----------------------
c
c           fval   - value of the function f at the points
c                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
c                    k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (bl,br)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (bl,br)
c           res12  - approximation to the integral obtained from cheb12
c           res24  - approximation to the integral obtained from cheb24
c           qwgts - external function subprogram defining
c                    the four possible weight functions
c           hlgth  - half-length of the interval (bl,br)
c           centr  - mid point of the interval (bl,br)
c
c***first executable statement  qc25s
      nev = 25
      if(bl.eq.a.and.(alfa.ne.0.0e+00.or.integr.eq.2.or.integr.eq.4))
     * go to 10
      if(br.eq.b.and.(beta.ne.0.0e+00.or.integr.eq.3.or.integr.eq.4))
     * go to 140
c
c           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod
c           scheme.
c
c
      call qk15w(f,qwgts,a,b,alfa,beta,integr,bl,br,
     *    result,abserr,resabs,resasc)
      nev = 15
      go to 270
c
c           this part of the program is executed only if a = bl.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
c                  *f(0.5*(br-a)*x+0.5*(br+a))
c
   10 hlgth = 0.5e+00*(br-bl)
      centr = 0.5e+00*(br+bl)
      fix = b-centr
      fval(1) = 0.5e+00*f(hlgth+centr)*(fix-hlgth)**beta
      fval(13) = f(centr)*(fix**beta)
      fval(25) = 0.5e+00*f(centr-hlgth)*(fix+hlgth)**beta
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix-u)**beta
        fval(isym) = f(centr-u)*(fix+u)**beta
   20 continue
      factor = hlgth**(alfa+0.1e+01)
      result = 0.0e+00
      abserr = 0.0e+00
      res12 = 0.0e+00
      res24 = 0.0e+00
      if(integr.gt.2) go to 70
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 1  (or 2)
c
      do 30 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   30 continue
      do 40 i=14,25
        res24 = res24+cheb24(i)*ri(i)
   40 continue
      if(integr.eq.1) go to 130
c
c           integr = 2
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 50 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res12+cheb24(i)*rg(i)
   50 continue
      do 60 i=14,25
        res24 = res24+cheb24(i)*rg(i)
   60 continue
      go to 130
c
c           compute the chebyshev series expansion of the
c           following function
c           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
c
   70 fval(1) = fval(1)*alog(fix-hlgth)
      fval(13) = fval(13)*alog(fix)
      fval(25) = fval(25)*alog(fix+hlgth)
      do 80 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*alog(fix-u)
        fval(isym) = fval(isym)*alog(fix+u)
   80 continue
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 3  (or 4)
c
      do 90 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   90 continue
      do 100 i=14,25
        res24 = res24+cheb24(i)*ri(i)
  100 continue
      if(integr.eq.3) go to 130
c
c           integr = 4
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 110 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
  110 continue
      do 120 i=14,25
        res24 = res24+cheb24(i)*rg(i)
  120 continue
  130 result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor
      go to 270
c
c           this part of the program is executed only if b = br.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
c                *f(0.5*(b-bl)*x+0.5*(b+bl))
c
  140 hlgth = 0.5e+00*(br-bl)
      centr = 0.5e+00*(br+bl)
      fix = centr-a
      fval(1) = 0.5e+00*f(hlgth+centr)*(fix+hlgth)**alfa
      fval(13) = f(centr)*(fix**alfa)
      fval(25) = 0.5e+00*f(centr-hlgth)*(fix-hlgth)**alfa
      do 150 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix+u)**alfa
        fval(isym) = f(centr-u)*(fix-u)**alfa
  150 continue
      factor = hlgth**(beta+0.1e+01)
      result = 0.0e+00
      abserr = 0.0e+00
      res12 = 0.0e+00
      res24 = 0.0e+00
      if(integr.eq.2.or.integr.eq.4) go to 200
c
c           integr = 1  (or 3)
c
      call qcheb(x,fval,cheb12,cheb24)
      do 160 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  160 continue
      do 170 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  170 continue
      if(integr.eq.1) go to 260
c
c           integr = 3
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 180 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  180 continue
      do 190 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  190 continue
      go to 260
c
c           compute the chebyshev series expansion of the
c           following function
c           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
c
  200 fval(1) = fval(1)*alog(hlgth+fix)
      fval(13) = fval(13)*alog(fix)
      fval(25) = fval(25)*alog(fix-hlgth)
      do 210 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*alog(u+fix)
        fval(isym) = fval(isym)*alog(fix-u)
  210 continue
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 2  (or 4)
c
      do 220 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  220 continue
      do 230 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  230 continue
      if(integr.eq.2) go to 260
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
c
c           integr = 4
c
      do 240 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  240 continue
      do 250 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  250 continue
  260 result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor
  270 return
      end
      subroutine qcheb(x,fval,cheb12,cheb24)
c***begin prologue  qcheb
c***refer to  qc25c,qc25f,qc25s
c***routines called  (none)
c***revision date  830518   (yymmdd)
c***keywords  chebyshev series expansion, fast fourier transform
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes the chebyshev series expansion
c            of degrees 12 and 24 of a function using a
c            fast fourier transform method
c            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
c            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
c            where t(k,x) is the chebyshev polynomial of degree k.
c***description
c
c        chebyshev series expansion
c        standard fortran subroutine
c        real version
c
c        parameters
c          on entry
c           x      - real
c                    vector of dimension 11 containing the
c                    values cos(k*pi/24), k = 1, ..., 11
c
c           fval   - real
c                    vector of dimension 25 containing the
c                    function values at the points
c                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
c                    where (a,b) is the approximation interval.
c                    fval(1) and fval(25) are divided by two
c                    (these values are destroyed at output).
c
c          on return
c           cheb12 - real
c                    vector of dimension 13 containing the
c                    chebyshev coefficients for degree 12
c
c           cheb24 - real
c                    vector of dimension 25 containing the
c                    chebyshev coefficients for degree 24
c
c***end prologue  qcheb
c
      real alam,alam1,alam2,cheb12,cheb24,
     *  fval,part1,part2,part3,v,x
      integer i,j
c
      dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)
c
c***first executable statement  qcheb
      do 10 i=1,12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   10 continue
      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)
     *  +x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)
     *  +x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)
     *  -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)
     *  -x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam
      do 20 i=1,6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   20 continue
      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam
      do 30 i=1,3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   30 continue
      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 0.1e+01/0.6e+01
      do 40 i=2,12
        cheb12(i) = cheb12(i)*alam
   40 continue
      alam = 0.5e+00*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam
      do 50 i=2,24
        cheb24(i) = cheb24(i)*alam
   50 continue
      cheb24(1) = 0.5e+00*alam*cheb24(1)
      cheb24(25) = 0.5e+00*alam*cheb24(25)
      return
      end
      subroutine qelg(n,epstab,result,abserr,res3la,nres)
c***begin prologue  qelg
c***refer to  qagie,qagoe,qagpe,qagse
c***routines called  r1mach
c***revision date  830518   (yymmdd)
c***keywords  epsilon algorithm, convergence acceleration,
c             extrapolation
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine determines the limit of a given sequence of
c            approximations, by means of the epsilon algorithm of
c            p. wynn. an estimate of the absolute error is also given.
c            the condensed epsilon table is computed. only those
c            elements needed for the computation of the next diagonal
c            are preserved.
c***description
c
c           epsilon algorithm
c           standard fortran subroutine
c           real version
c
c           parameters
c              n      - integer
c                       epstab(n) contains the new element in the
c                       first column of the epsilon table.
c
c              epstab - real
c                       vector of dimension 52 containing the elements
c                       of the two lower diagonals of the triangular
c                       epsilon table. the elements are numbered
c                       starting at the right-hand corner of the
c                       triangle.
c
c              result - real
c                       resulting approximation to the integral
c
c              abserr - real
c                       estimate of the absolute error computed from
c                       result and the 3 previous results
c
c              res3la - real
c                       vector of dimension 3 containing the last 3
c                       results
c
c              nres   - integer
c                       number of calls to the routine
c                       (should be zero at first call)
c
c***end prologue  qelg
c
      real abserr,delta1,delta2,delta3,r1mach,
     *  epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
     *  oflow,res,result,res3la,ss,tol1,tol2,tol3
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
c
c           list of major variables
c           -----------------------
c
c           e0     - the 4 elements on which the
c           e1       computation of a new element in
c           e2       the epsilon table is based
c           e3                 e0
c                        e3    e1    new
c                              e2
c           newelm - number of elements to be computed in the new
c                    diagonal
c           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
c           result - the element in the new diagonal with least value
c                    of error
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           oflow is the largest positive magnitude.
c           limexp is the maximum number of elements the epsilon
c           table can contain. if this number is reached, the upper
c           diagonal of the epsilon table is deleted.
c
c***first executable statement  qelg
      epmach = r1mach(4)
      oflow = r1mach(2)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = abs(e1)
        delta2 = e2-e1
        err2 = abs(delta2)
        tol2 = amax1(abs(e2),e1abs)*epmach
        delta3 = e1-e0
        err3 = abs(delta3)
        tol3 = amax1(e1abs,abs(e0))*epmach
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10
c
c           if e0, e1 and e2 are equal to within machine
c           accuracy, convergence is assumed.
c           result = e2
c           abserr = abs(e1-e0)+abs(e2-e1)
c
        result = res
        abserr = err2+err3
c ***jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = abs(delta1)
        tol1 = amax1(e1abs,abs(e3))*epmach
c
c           if two elements are very close to each other, omit
c           a part of the table by adjusting the value of n
c
        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 0.1e+01/delta1+0.1e+01/delta2-0.1e+01/delta3
        epsinf = abs(ss*e1)
c
c           test to detect irregular behaviour in the table, and
c           eventually omit a part of the table adjusting the value
c           of n.
c
        if(epsinf.gt.0.1e-03) go to 30
   20   n = i+i-1
c ***jump out of do-loop
        go to 50
c
c           compute a new element and eventually adjust
c           the value of result.
c
   30   res = e1+0.1e+01/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+abs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        result = res
   40 continue
c
c           shift the table.
c
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = result
      abserr = oflow
      go to 100
c
c           compute error estimate
c
   90 abserr = abs(result-res3la(3))+abs(result-res3la(2))
     *  +abs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = amax1(abserr,0.5e+01*epmach*abs(result))
      return
      end
      subroutine qk15(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk15
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the7-point gauss rule(resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk15
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/
     *     0.9914553711208126e+00,   0.9491079123427585e+00,
     *     0.8648644233597691e+00,   0.7415311855993944e+00,
     *     0.5860872354676911e+00,   0.4058451513773972e+00,
     *     0.2077849550078985e+00,   0.0e+00              /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/
     *     0.2293532201052922e-01,   0.6309209262997855e-01,
     *     0.1047900103222502e+00,   0.1406532597155259e+00,
     *     0.1690047266392679e+00,   0.1903505780647854e+00,
     *     0.2044329400752989e+00,   0.2094821410847278e+00/
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697e+00,   0.2797053914892767e+00,
     *     0.3818300505051189e+00,   0.4179591836734694e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk15
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = abs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk15i
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a2,h2a4a2
c***keywords  15-point transformed gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the original (infinite integration range is mapped
c            onto the interval (0,1) and (a,b) is a part of (0,1).
c            it is the purpose to compute
c            i = integral of transformed integrand over (a,b),
c            j = integral of abs(transformed integrand) over (a,b).
c***description
c
c           integration rule
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       fuction subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              boun   - real
c                       finite bound of original integration
c                       range (set to zero if inf = +2)
c
c              inf    - integer
c                       if inf = -1, the original interval is
c                                   (-infinity,bound),
c                       if inf = +1, the original interval is
c                                   (bound,+infinity),
c                       if inf = +2, the original interval is
c                                   (-infinity,+infinity) and
c                       the integral is computed as the sum of two
c                       integrals, one over (-infinity,0) and one over
c                       (0,+infinity).
c
c              a      - real
c                       lower limit for integration over subrange
c                       of (0,1)
c
c              b      - real
c                       upper limit for integration over subrange
c                       of (0,1)
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule(resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule(resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integral of
c                       abs((transformed integrand)-i/(b-a)) over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk15i
c
      real a,absc,absc1,absc2,abserr,b,boun,centr,
     *  dinf,r1mach,epmach,f,fc,fsum,fval1,fval2,fv1,
     *  fv2,hlgth,resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,
     *  uflow,wg,wgk,xgk
      integer inf,j,min0
      external f
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
c
c           the abscissae and weights are supplied for the interval
c           (-1,1).  because of symmetry only the positive abscissae and
c           their corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule, corresponding
c                    to the abscissae xgk(2), xgk(4), ...
c                    wg(1), wg(3), ... are set to zero.
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),
     *  xgk(8)/
     *     0.9914553711208126e+00,     0.9491079123427585e+00,
     *     0.8648644233597691e+00,     0.7415311855993944e+00,
     *     0.5860872354676911e+00,     0.4058451513773972e+00,
     *     0.2077849550078985e+00,     0.0000000000000000e+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),
     *  wgk(8)/
     *     0.2293532201052922e-01,     0.6309209262997855e-01,
     *     0.1047900103222502e+00,     0.1406532597155259e+00,
     *     0.1690047266392679e+00,     0.1903505780647854e+00,
     *     0.2044329400752989e+00,     0.2094821410847278e+00/
c
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/
     *     0.0000000000000000e+00,     0.1294849661688697e+00,
     *     0.0000000000000000e+00,     0.2797053914892767e+00,
     *     0.0000000000000000e+00,     0.3818300505051189e+00,
     *     0.0000000000000000e+00,     0.4179591836734694e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           tabsc* - transformed abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of the transformed
c                    integrand over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk15i
      epmach = r1mach(4)
      uflow = r1mach(1)
      dinf = min0(1,inf)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      tabsc1 = boun+dinf*(0.1e+01-centr)/centr
      fval1 = f(tabsc1)
      if(inf.eq.2) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the error.
c
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)
      do 10 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(0.1e+01-absc1)/absc1
        tabsc2 = boun+dinf*(0.1e+01-absc2)/absc2
        fval1 = f(tabsc1)
        fval2 = f(tabsc2)
        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
        if(inf.eq.2) fval2 = fval2+f(-tabsc2)
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
   10 continue
      reskh = resk*0.5e+00
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.e0) abserr = resasc*
     * amin1(0.1e+01,(0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     * ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr,
     *   resabs,resasc)
c***begin prologue  qk15w
c***date written   810101   (yymmdd)
c***revision date  830518   (mmddyy)
c***category no.  h2a2a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b), with error
c                           estimate
c                       j = integral of abs(f*w) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c             on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              w      - real
c                       function subprogram defining the integrand
c                       weight function w(x). the actual name for w
c                       needs to be declared e x t e r n a l in the
c                       calling program.
c
c              p1, p2, p3, p4 - real
c                       parameters in the weight function
c
c              kp     - integer
c                       key for indicating the type of weight function
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral of abs(f)
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk15w
c
      real a,absc,absc1,absc2,abserr,b,centr,dhlgth,
     *  r1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,
     *  hlgth,p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,
     *  w,wg,wgk,xgk
      integer j,jtw,jtwm1,kp
      external f,w
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point gauss-kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ... abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point gauss-kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),
     *  xgk(8)/
     *     0.9914553711208126e+00,     0.9491079123427585e+00,
     *     0.8648644233597691e+00,     0.7415311855993944e+00,
     *     0.5860872354676911e+00,     0.4058451513773972e+00,
     *     0.2077849550078985e+00,     0.0000000000000000e+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),
     *  wgk(8)/
     *     0.2293532201052922e-01,     0.6309209262997855e-01,
     *     0.1047900103222502e+00,     0.1406532597155259e+00,
     *     0.1690047266392679e+00,     0.1903505780647854e+00,
     *     0.2044329400752989e+00,     0.2094821410847278e+00/
c
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697e+00,    0.2797053914892767e+00,
     *     0.3818300505051189e+00,    0.4179591836734694e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f*w over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk15w
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 15-point kronrod approximation to the
c           integral, and estimate the error.
c
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1((epmach*
     *  0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk21(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk21
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resg,resk,reskh,result,r1mach,uflow,wg,wgk,
     *  xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),
     *  xgk(8),xgk(9),xgk(10),xgk(11)/
     *         0.9956571630258081e+00,     0.9739065285171717e+00,
     *     0.9301574913557082e+00,     0.8650633666889845e+00,
     *     0.7808177265864169e+00,     0.6794095682990244e+00,
     *     0.5627571346686047e+00,     0.4333953941292472e+00,
     *     0.2943928627014602e+00,     0.1488743389816312e+00,
     *     0.0000000000000000e+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),
     *  wgk(8),wgk(9),wgk(10),wgk(11)/
     *     0.1169463886737187e-01,     0.3255816230796473e-01,
     *     0.5475589657435200e-01,     0.7503967481091995e-01,
     *     0.9312545458369761e-01,     0.1093871588022976e+00,
     *     0.1234919762620659e+00,     0.1347092173114733e+00,
     *     0.1427759385770601e+00,     0.1477391049013385e+00,
     *     0.1494455540029169e+00/
c
      data wg(1),wg(2),wg(3),wg(4),wg(5)/
     *     0.6667134430868814e-01,     0.1494513491505806e+00,
     *     0.2190863625159820e+00,     0.2692667193099964e+00,
     *     0.2955242247147529e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk21
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0e+00
      fc = f(centr)
      resk = wgk(11)*fc
      resabs = abs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(11)*abs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk31(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk31
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  31-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 31-point
c                       gauss-kronrod rule (resk), obtained by optimal
c                       addition of abscissae to the 15-point gauss
c                       rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the modulus,
c                       which should not exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk31
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 31-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 15-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 15-point gauss rule
c
c           wgk    - weights of the 31-point kronrod rule
c
c           wg     - weights of the 15-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),
     *  xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),
     *  xgk(16)/
     *     0.9980022986933971e+00,   0.9879925180204854e+00,
     *     0.9677390756791391e+00,   0.9372733924007059e+00,
     *     0.8972645323440819e+00,   0.8482065834104272e+00,
     *     0.7904185014424659e+00,   0.7244177313601700e+00,
     *     0.6509967412974170e+00,   0.5709721726085388e+00,
     *     0.4850818636402397e+00,   0.3941513470775634e+00,
     *     0.2991800071531688e+00,   0.2011940939974345e+00,
     *     0.1011420669187175e+00,   0.0e+00               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),
     *  wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),
     *  wgk(16)/
     *     0.5377479872923349e-02,   0.1500794732931612e-01,
     *     0.2546084732671532e-01,   0.3534636079137585e-01,
     *     0.4458975132476488e-01,   0.5348152469092809e-01,
     *     0.6200956780067064e-01,   0.6985412131872826e-01,
     *     0.7684968075772038e-01,   0.8308050282313302e-01,
     *     0.8856444305621177e-01,   0.9312659817082532e-01,
     *     0.9664272698362368e-01,   0.9917359872179196e-01,
     *     0.1007698455238756e+00,   0.1013300070147915e+00/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/
     *     0.3075324199611727e-01,   0.7036604748810812e-01,
     *     0.1071592204671719e+00,   0.1395706779261543e+00,
     *     0.1662692058169939e+00,   0.1861610000155622e+00,
     *     0.1984314853271116e+00,   0.2025782419255613e+00/
c
c
c           list of major variables
c           -----------------------
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 15-point gauss formula
c           resk   - result of the 31-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk31
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 31-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = abs(resk)
      do 10 j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(16)*abs(fc-reskh)
      do 20 j=1,15
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk41(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk41
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  41-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 41-point
c                       gauss-kronrod rule (resk) obtained by optimal
c                       addition of abscissae to the 20-point gauss
c                       rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integal of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk41
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,
     *  resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 41-point gauss-kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 20-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 20-point gauss rule
c
c           wgk    - weights of the 41-point gauss-kronrod rule
c
c           wg     - weights of the 20-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),
     *  xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),
     *  xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/
     *     0.9988590315882777e+00,   0.9931285991850949e+00,
     *     0.9815078774502503e+00,   0.9639719272779138e+00,
     *     0.9408226338317548e+00,   0.9122344282513259e+00,
     *     0.8782768112522820e+00,   0.8391169718222188e+00,
     *     0.7950414288375512e+00,   0.7463319064601508e+00,
     *     0.6932376563347514e+00,   0.6360536807265150e+00,
     *     0.5751404468197103e+00,   0.5108670019508271e+00,
     *     0.4435931752387251e+00,   0.3737060887154196e+00,
     *     0.3016278681149130e+00,   0.2277858511416451e+00,
     *     0.1526054652409227e+00,   0.7652652113349733e-01,
     *     0.0e+00               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),
     *  wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),
     *  wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/
     *     0.3073583718520532e-02,   0.8600269855642942e-02,
     *     0.1462616925697125e-01,   0.2038837346126652e-01,
     *     0.2588213360495116e-01,   0.3128730677703280e-01,
     *     0.3660016975820080e-01,   0.4166887332797369e-01,
     *     0.4643482186749767e-01,   0.5094457392372869e-01,
     *     0.5519510534828599e-01,   0.5911140088063957e-01,
     *     0.6265323755478117e-01,   0.6583459713361842e-01,
     *     0.6864867292852162e-01,   0.7105442355344407e-01,
     *     0.7303069033278667e-01,   0.7458287540049919e-01,
     *     0.7570449768455667e-01,   0.7637786767208074e-01,
     *     0.7660071191799966e-01/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/
     *     0.1761400713915212e-01,    0.4060142980038694e-01,
     *     0.6267204833410906e-01,    0.8327674157670475e-01,
     *     0.1019301198172404e+00,    0.1181945319615184e+00,
     *     0.1316886384491766e+00,    0.1420961093183821e+00,
     *     0.1491729864726037e+00,    0.1527533871307259e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 20-point gauss formula
c           resk   - result of the 41-point kronrod formula
c           reskh  - approximation to mean value of f over (a,b), i.e.
c                    to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk41
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 41-point gauss-kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0e+00
      fc = f(centr)
      resk = wgk(21)*fc
      resabs = abs(resk)
      do 10 j=1,10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(21)*abs(fc-reskh)
      do 20 j=1,20
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk51(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk51
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  51-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c            on entry
c              f      - real
c                       function subroutine defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 51-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 25-point gauss rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral j
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk51
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 51-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 25-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 25-point gauss rule
c
c           wgk    - weights of the 51-point kronrod rule
c
c           wg     - weights of the 25-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),
     *  xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/
     *     0.9992621049926098e+00,   0.9955569697904981e+00,
     *     0.9880357945340772e+00,   0.9766639214595175e+00,
     *     0.9616149864258425e+00,   0.9429745712289743e+00,
     *     0.9207471152817016e+00,   0.8949919978782754e+00,
     *     0.8658470652932756e+00,   0.8334426287608340e+00,
     *     0.7978737979985001e+00,   0.7592592630373576e+00,
     *     0.7177664068130844e+00,   0.6735663684734684e+00/
       data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21),
     *  xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/
     *     0.6268100990103174e+00,   0.5776629302412230e+00,
     *     0.5263252843347192e+00,   0.4730027314457150e+00,
     *     0.4178853821930377e+00,   0.3611723058093878e+00,
     *     0.3030895389311078e+00,   0.2438668837209884e+00,
     *     0.1837189394210489e+00,   0.1228646926107104e+00,
     *     0.6154448300568508e-01,   0.0e+00               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),
     *  wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/
     *     0.1987383892330316e-02,   0.5561932135356714e-02,
     *     0.9473973386174152e-02,   0.1323622919557167e-01,
     *     0.1684781770912830e-01,   0.2043537114588284e-01,
     *     0.2400994560695322e-01,   0.2747531758785174e-01,
     *     0.3079230016738749e-01,   0.3400213027432934e-01,
     *     0.3711627148341554e-01,   0.4008382550403238e-01,
     *     0.4287284502017005e-01,   0.4550291304992179e-01/
       data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)
     *  ,wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/
     *     0.4798253713883671e-01,   0.5027767908071567e-01,
     *     0.5236288580640748e-01,   0.5425112988854549e-01,
     *     0.5595081122041232e-01,   0.5743711636156783e-01,
     *     0.5868968002239421e-01,   0.5972034032417406e-01,
     *     0.6053945537604586e-01,   0.6112850971705305e-01,
     *     0.6147118987142532e-01,   0.6158081806783294e-01/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),
     *  wg(10),wg(11),wg(12),wg(13)/
     *     0.1139379850102629e-01,    0.2635498661503214e-01,
     *     0.4093915670130631e-01,    0.5490469597583519e-01,
     *     0.6803833381235692e-01,    0.8014070033500102e-01,
     *     0.9102826198296365e-01,    0.1005359490670506e+00,
     *     0.1085196244742637e+00,    0.1148582591457116e+00,
     *     0.1194557635357848e+00,    0.1222424429903100e+00,
     *     0.1231760537267155e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 25-point gauss formula
c           resk   - result of the 51-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk51
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 51-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      fc = f(centr)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = abs(resk)
      do 10 j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(26)*abs(fc-reskh)
      do 20 j=1,25
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qk61(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk61
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  61-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of dabs(f) over (a,b)
c***description
c
c        integration rule
c        standard fortran subroutine
c        real version
c
c
c        parameters
c         on entry
c           f      - real
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to be
c                    declared e x t e r n a l in the calling program.
c
c           a      - real
c                    lower limit of integration
c
c           b      - real
c                    upper limit of integration
c
c         on return
c           result - real
c                    approximation to the integral i
c                    result is computed by applying the 61-point
c                    kronrod rule (resk) obtained by optimal addition of
c                    abscissae to the 30-point gauss rule (resg).
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed dabs(i-result)
c
c           resabs - real
c                    approximation to the integral j
c
c           resasc - real
c                    approximation to the integral of dabs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk61
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
c
c           the abscissae and weights are given for the
c           interval (-1,1). because of symmetry only the positive
c           abscissae and their corresponding weights are given.
c
c           xgk   - abscissae of the 61-point kronrod rule
c                   xgk(2), xgk(4)  ... abscissae of the 30-point
c                   gauss rule
c                   xgk(1), xgk(3)  ... optimally added abscissae
c                   to the 30-point gauss rule
c
c           wgk   - weights of the 61-point kronrod rule
c
c           wg    - weigths of the 30-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),
     *   xgk(9),xgk(10)/
     *     0.9994844100504906e+00,     0.9968934840746495e+00,
     *     0.9916309968704046e+00,     0.9836681232797472e+00,
     *     0.9731163225011263e+00,     0.9600218649683075e+00,
     *     0.9443744447485600e+00,     0.9262000474292743e+00,
     *     0.9055733076999078e+00,     0.8825605357920527e+00/
      data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),
     *  xgk(17),xgk(18),xgk(19),xgk(20)/
     *     0.8572052335460611e+00,     0.8295657623827684e+00,
     *     0.7997278358218391e+00,     0.7677774321048262e+00,
     *     0.7337900624532268e+00,     0.6978504947933158e+00,
     *     0.6600610641266270e+00,     0.6205261829892429e+00,
     *     0.5793452358263617e+00,     0.5366241481420199e+00/
      data xgk(21),xgk(22),xgk(23),xgk(24),
     *  xgk(25),xgk(26),xgk(27),xgk(28),xgk(29),xgk(30),xgk(31)/
     *     0.4924804678617786e+00,     0.4470337695380892e+00,
     *     0.4004012548303944e+00,     0.3527047255308781e+00,
     *     0.3040732022736251e+00,     0.2546369261678898e+00,
     *     0.2045251166823099e+00,     0.1538699136085835e+00,
     *     0.1028069379667370e+00,     0.5147184255531770e-01,
     *     0.0e+00                   /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),
     *  wgk(9),wgk(10)/
     *     0.1389013698677008e-02,     0.3890461127099884e-02,
     *     0.6630703915931292e-02,     0.9273279659517763e-02,
     *     0.1182301525349634e-01,     0.1436972950704580e-01,
     *     0.1692088918905327e-01,     0.1941414119394238e-01,
     *     0.2182803582160919e-01,     0.2419116207808060e-01/
      data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),
     *  wgk(17),wgk(18),wgk(19),wgk(20)/
     *     0.2650995488233310e-01,     0.2875404876504129e-01,
     *     0.3090725756238776e-01,     0.3298144705748373e-01,
     *     0.3497933802806002e-01,     0.3688236465182123e-01,
     *     0.3867894562472759e-01,     0.4037453895153596e-01,
     *     0.4196981021516425e-01,     0.4345253970135607e-01/
      data wgk(21),wgk(22),wgk(23),wgk(24),
     *  wgk(25),wgk(26),wgk(27),wgk(28),wgk(29),wgk(30),wgk(31)/
     *     0.4481480013316266e-01,     0.4605923827100699e-01,
     *     0.4718554656929915e-01,     0.4818586175708713e-01,
     *     0.4905543455502978e-01,     0.4979568342707421e-01,
     *     0.5040592140278235e-01,     0.5088179589874961e-01,
     *     0.5122154784925877e-01,     0.5142612853745903e-01,
     *     0.5149472942945157e-01/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/
     *     0.7968192496166606e-02,     0.1846646831109096e-01,
     *     0.2878470788332337e-01,     0.3879919256962705e-01,
     *     0.4840267283059405e-01,     0.5749315621761907e-01,
     *     0.6597422988218050e-01,     0.7375597473770521e-01/
      data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/
     *     0.8075589522942022e-01,     0.8689978720108298e-01,
     *     0.9212252223778613e-01,     0.9636873717464426e-01,
     *     0.9959342058679527e-01,     0.1017623897484055e+00,
     *     0.1028526528935588e+00/
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 30-point gauss rule
c           resk   - result of the 61-point kronrod rule
c           reskh  - approximation to the mean value of f
c                    over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk61
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(b+a)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 61-point kronrod approximation to the
c           integral, and estimate the absolute error.
c
      resg = 0.0e+00
      fc = f(centr)
      resk = wgk(31)*fc
      resabs = abs(resk)
      do 10 j=1,15
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  15    continue
      reskh = resk*0.5e+00
      resasc = wgk(31)*abs(fc-reskh)
      do 20 j=1,30
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
      subroutine qmomo(alfa,beta,ri,rj,rg,rh,integr)
c***begin prologue  qmomo
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,c3a2
c***keywords  modified chebyshev moments
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes modified chebsyshev moments. the k-th
c            modified chebyshev moment is defined as the integral over
c            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
c            polynomial of degree k.
c***description
c
c        modified chebyshev moments
c        standard fortran subroutine
c        real version
c
c        parameters
c           alfa   - real
c                    parameter in the weight function w(x), alfa.gt.(-1)
c
c           beta   - real
c                    parameter in the weight function w(x), beta.gt.(-1)
c
c           ri     - real
c                    vector of dimension 25
c                    ri(k) is the integral over (-1,1) of
c                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
c
c           rj     - real
c                    vector of dimension 25
c                    rj(k) is the integral over (-1,1) of
c                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
c
c           rg     - real
c                    vector of dimension 25
c                    rg(k) is the integral over (-1,1) of
c                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           rh     - real
c                    vector of dimension 25
c                    rh(k) is the integral over (-1,1) of
c                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           integr - integer
c                    input parameter indicating the modified
c                    moments to be computed
c                    integr = 1 compute ri, rj
c                           = 2 compute ri, rj, rg
c                           = 3 compute ri, rj, rh
c                           = 4 compute ri, rj, rg, rh
c***references  (none)
c***routines called  (none)
c***end prologue  qmomo
c
      real alfa,alfp1,alfp2,an,anm1,beta,betp1,
     *  betp2,ralf,rbet,rg,rh,ri,rj
      integer i,im1,integr
c
      dimension rg(25),rh(25),ri(25),rj(25)
c
c
c***first executable statement  qmomo
      alfp1 = alfa+0.1e+01
      betp1 = beta+0.1e+01
      alfp2 = alfa+0.2e+01
      betp2 = beta+0.2e+01
      ralf = 0.2e+01**alfp1
      rbet = 0.2e+01**betp1
c
c           compute ri, rj using a forward recurrence relation.
c
      ri(1) = ralf/alfp1
      rj(1) = rbet/betp1
      ri(2) = ri(1)*alfa/alfp2
      rj(2) = rj(1)*beta/betp2
      an = 0.2e+01
      anm1 = 0.1e+01
      do 20 i=3,25
        ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/
     *  (anm1*(an+alfp1))
        rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/
     *  (anm1*(an+betp1))
        anm1 = an
        an = an+0.1e+01
   20 continue
      if(integr.eq.1) go to 70
      if(integr.eq.3) go to 40
c
c           compute rg using a forward recurrence relation.
c
      rg(1) = -ri(1)/alfp1
      rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
      an = 0.2e+01
      anm1 = 0.1e+01
      im1 = 2
      do 30 i=3,25
        rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/
     *  (anm1*(an+alfp1))
        anm1 = an
        an = an+0.1e+01
        im1 = i
   30 continue
      if(integr.eq.2) go to 70
c
c           compute rh using a forward recurrence relation.
c
   40 rh(1) = -rj(1)/betp1
      rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
      an = 0.2e+01
      anm1 = 0.1e+01
      im1 = 2
      do 50 i=3,25
        rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+
     *  anm1*rj(i))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1e+01
        im1 = i
   50 continue
      do 60 i=2,25,2
        rh(i) = -rh(i)
   60 continue
   70 do 80 i=2,25,2
        rj(i) = -rj(i)
   80 continue
   90 return
      end
      subroutine qng(f,a,b,epsabs,epsrel,result,abserr,neval,ier)
c***begin prologue  qng
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, smooth integrand,
c             non-adaptive, gauss-kronrod(patterson)
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl math & progr. div. - k.u.leuven
c           kahaner,david,nbs - modified (2/82)
c***purpose  the routine calculates an approximation result to a
c            given definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c non-adaptive integration
c standard fortran subroutine
c real version
c
c           f      - real version
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l in the driver program.
c
c           a      - real version
c                    lower limit of integration
c
c           b      - real version
c                    upper limit of integration
c
c           epsabs - real
c                    absolute accuracy requested
c           epsrel - real
c                    relative accuracy requested
c                    if  epsabs.le.0
c                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                    the routine will end with ier = 6.
c
c         on return
c           result - real
c                    approximation to the integral i
c                    result is obtained by applying the 21-point
c                    gauss-kronrod rule (res21) obtained by optimal
c                    addition of abscissae to the 10-point gauss rule
c                    (res10), or by applying the 43-point rule (res43)
c                    obtained by optimal addition of abscissae to the
c                    21-point gauss-kronrod rule, or by applying the
c                    87-point rule (res87) obtained by optimal addition
c                    of abscissae to the 43-point rule.
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           ier    - ier = 0 normal and reliable termination of the
c                            routine. it is assumed that the requested
c                            accuracy has been achieved.
c                    ier.gt.0 abnormal termination of the routine. it is
c                            assumed that the requested accuracy has
c                            not been achieved.
c           error messages
c                    ier = 1 the maximum number of steps has been
c                            executed. the integral is probably too
c                            difficult to be calculated by dqng.
c                        = 6 the input is invalid, because
c                            epsabs.le.0 and
c                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                            result, abserr and neval are set to zero.
c
c***references  (none)
c***routines called  r1mach,xerror
c***end prologue  qng
c
      real a,absc,abserr,b,centr,dhlgth,epmach,epsabs,epsrel,f,fcentr,
     *  fval,fval1,fval2,fv1,fv2,fv3,fv4,hlgth,result,res10,res21,res43,
     *  res87,resabs,resasc,reskh,r1mach,savfun,uflow,w10,w21a,w43a,
     *  w43b,w87a,w87b,x1,x2,x3,x4
      integer ier,ipx,k,l,neval
      external f
c
      dimension fv1(5),fv2(5),fv3(5),fv4(5),x1(5),x2(5),x3(11),x4(22),
     *  w10(5),w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23),
     *  savfun(21)
c
c           the following data statements contain the
c           abscissae and weights of the integration rules used.
c
c           x1      abscissae common to the 10-, 21-, 43-
c                   and 87-point rule
c           x2      abscissae common to the 21-, 43- and
c                   87-point rule
c           x3      abscissae common to the 43- and 87-point
c                   rule
c           x4      abscissae of the 87-point rule
c           w10     weights of the 10-point formula
c           w21a    weights of the 21-point formula for
c                   abscissae x1
c           w21b    weights of the 21-point formula for
c                   abscissae x2
c           w43a    weights of the 43-point formula for
c                   abscissae x1, x3
c           w43b    weights of the 43-point formula for
c                   abscissae x3
c           w87a    weights of the 87-point formula for
c                   abscissae x1, x2, x3
c           w87b    weights of the 87-point formula for
c                   abscissae x4
c
      data x1(1),x1(2),x1(3),x1(4),x1(5)/
     *     0.9739065285171717e+00,     0.8650633666889845e+00,
     *     0.6794095682990244e+00,     0.4333953941292472e+00,
     *     0.1488743389816312e+00/
      data x2(1),x2(2),x2(3),x2(4),x2(5)/
     *     0.9956571630258081e+00,     0.9301574913557082e+00,
     *     0.7808177265864169e+00,     0.5627571346686047e+00,
     *     0.2943928627014602e+00/
      data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),
     *  x3(9),x3(10),x3(11)/
     *     0.9993333609019321e+00,     0.9874334029080889e+00,
     *     0.9548079348142663e+00,     0.9001486957483283e+00,
     *     0.8251983149831142e+00,     0.7321483889893050e+00,
     *     0.6228479705377252e+00,     0.4994795740710565e+00,
     *     0.3649016613465808e+00,     0.2222549197766013e+00,
     *     0.7465061746138332e-01/
      data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),
     *  x4(10),x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),
     *  x4(19),x4(20),x4(21),x4(22)/   0.9999029772627292e+00,
     *     0.9979898959866787e+00,     0.9921754978606872e+00,
     *     0.9813581635727128e+00,     0.9650576238583846e+00,
     *     0.9431676131336706e+00,     0.9158064146855072e+00,
     *     0.8832216577713165e+00,     0.8457107484624157e+00,
     *     0.8035576580352310e+00,     0.7570057306854956e+00,
     *     0.7062732097873218e+00,     0.6515894665011779e+00,
     *     0.5932233740579611e+00,     0.5314936059708319e+00,
     *     0.4667636230420228e+00,     0.3994248478592188e+00,
     *     0.3298748771061883e+00,     0.2585035592021616e+00,
     *     0.1856953965683467e+00,     0.1118422131799075e+00,
     *     0.3735212339461987e-01/
      data w10(1),w10(2),w10(3),w10(4),w10(5)/
     *     0.6667134430868814e-01,     0.1494513491505806e+00,
     *     0.2190863625159820e+00,     0.2692667193099964e+00,
     *     0.2955242247147529e+00/
      data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/
     *     0.3255816230796473e-01,     0.7503967481091995e-01,
     *     0.1093871588022976e+00,     0.1347092173114733e+00,
     *     0.1477391049013385e+00/
      data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/
     *     0.1169463886737187e-01,     0.5475589657435200e-01,
     *     0.9312545458369761e-01,     0.1234919762620659e+00,
     *     0.1427759385770601e+00,     0.1494455540029169e+00/
      data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7),
     *  w43a(8),w43a(9),w43a(10)/      0.1629673428966656e-01,
     *     0.3752287612086950e-01,     0.5469490205825544e-01,
     *     0.6735541460947809e-01,     0.7387019963239395e-01,
     *     0.5768556059769796e-02,     0.2737189059324884e-01,
     *     0.4656082691042883e-01,     0.6174499520144256e-01,
     *     0.7138726726869340e-01/
      data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),
     *  w43b(7),w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/
     *     0.1844477640212414e-02,     0.1079868958589165e-01,
     *     0.2189536386779543e-01,     0.3259746397534569e-01,
     *     0.4216313793519181e-01,     0.5074193960018458e-01,
     *     0.5837939554261925e-01,     0.6474640495144589e-01,
     *     0.6956619791235648e-01,     0.7282444147183321e-01,
     *     0.7450775101417512e-01,     0.7472214751740301e-01/
      data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),
     *  w87a(7),w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),
     *  w87a(13),w87a(14),w87a(15),w87a(16),w87a(17),w87a(18),
     *  w87a(19),w87a(20),w87a(21)/
     *     0.8148377384149173e-02,     0.1876143820156282e-01,
     *     0.2734745105005229e-01,     0.3367770731163793e-01,
     *     0.3693509982042791e-01,     0.2884872430211531e-02,
     *     0.1368594602271270e-01,     0.2328041350288831e-01,
     *     0.3087249761171336e-01,     0.3569363363941877e-01,
     *     0.9152833452022414e-03,     0.5399280219300471e-02,
     *     0.1094767960111893e-01,     0.1629873169678734e-01,
     *     0.2108156888920384e-01,     0.2537096976925383e-01,
     *     0.2918969775647575e-01,     0.3237320246720279e-01,
     *     0.3478309895036514e-01,     0.3641222073135179e-01,
     *     0.3725387550304771e-01/
      data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7),
     *  w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14),
     *  w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),
     *  w87b(21),w87b(22),w87b(23)/    0.2741455637620724e-03,
     *     0.1807124155057943e-02,     0.4096869282759165e-02,
     *     0.6758290051847379e-02,     0.9549957672201647e-02,
     *     0.1232944765224485e-01,     0.1501044734638895e-01,
     *     0.1754896798624319e-01,     0.1993803778644089e-01,
     *     0.2219493596101229e-01,     0.2433914712600081e-01,
     *     0.2637450541483921e-01,     0.2828691078877120e-01,
     *     0.3005258112809270e-01,     0.3164675137143993e-01,
     *     0.3305041341997850e-01,     0.3425509970422606e-01,
     *     0.3526241266015668e-01,     0.3607698962288870e-01,
     *     0.3669860449845609e-01,     0.3712054926983258e-01,
     *     0.3733422875193504e-01,     0.3736107376267902e-01/
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fcentr - function value at mid point
c           absc   - abscissa
c           fval   - function value
c           savfun - array of function values which
c                    have already been computed
c           res10  - 10-point gauss result
c           res21  - 21-point kronrod result
c           res43  - 43-point result
c           res87  - 87-point result
c           resabs - approximation to the integral of abs(f)
c           resasc - approximation to the integral of abs(f-i/(b-a))
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qng
      epmach = r1mach(4)
      uflow = r1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      result = 0.0e+00
      abserr = 0.0e+00
      neval = 0
      ier = 6
      if(epsabs.le.0.0e+00.and.epsrel.lt.amax1(0.5e-14,0.5e+02*epmach))
     *  go to 80
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
      centr = 0.5e+00*(b+a)
      fcentr = f(centr)
      neval = 21
      ier = 1
c
c          compute the integral using the 10- and 21-point formula.
c
      do 70 l = 1,3
      go to (5,25,45),l
    5 res10 = 0.0e+00
      res21 = w21b(6)*fcentr
      resabs = w21b(6)*abs(fcentr)
      do 10 k=1,5
        absc = hlgth*x1(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1+fval2
        res10 = res10+w10(k)*fval
        res21 = res21+w21a(k)*fval
        resabs = resabs+w21a(k)*(abs(fval1)+abs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
   10 continue
      ipx = 5
      do 15 k=1,5
        ipx = ipx+1
        absc = hlgth*x2(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1+fval2
        res21 = res21+w21b(k)*fval
        resabs = resabs+w21b(k)*(abs(fval1)+abs(fval2))
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
   15 continue
c
c          test for convergence.
c
      result = res21*hlgth
      resabs = resabs*dhlgth
      reskh = 0.5e+00*res21
      resasc = w21b(6)*abs(fcentr-reskh)
      do 20 k = 1,5
        resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh))
     *                  +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
   20 continue
      abserr = abs((res21-res10)*hlgth)
      resasc = resasc*dhlgth
      go to 65
c
c          compute the integral using the 43-point formula.
c
   25 res43 = w43b(12)*fcentr
      neval = 43
      do 30 k=1,10
        res43 = res43+savfun(k)*w43a(k)
   30 continue
      do 40 k=1,11
        ipx = ipx+1
        absc = hlgth*x3(k)
        fval = f(absc+centr)+f(centr-absc)
        res43 = res43+fval*w43b(k)
        savfun(ipx) = fval
   40 continue
c
c          test for convergence.
c
      result = res43*hlgth
      abserr = abs((res43-res21)*hlgth)
      go to 65
c
c          compute the integral using the 87-point formula.
c
   45 res87 = w87b(23)*fcentr
      neval = 87
      do 50 k=1,21
        res87 = res87+savfun(k)*w87a(k)
   50 continue
      do 60 k=1,22
        absc = hlgth*x4(k)
        res87 = res87+w87b(k)*(f(absc+centr)+f(centr-absc))
   60 continue
      result = res87*hlgth
      abserr = abs((res87-res43)*hlgth)
   65 if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if (resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      if (abserr.le.amax1(epsabs,epsrel*abs(result))) ier = 0
c ***jump out of do-loop
      if (ier.eq.0) go to 999
   70 continue
   80 call xerror(26habnormal return from  qng ,26,ier,0)
  999 return
      end
      subroutine qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  qpsrt
c***refer to  qage,qagie,qagpe,qagse,qawce,qawse,qawoe
c***routines called  (none)
c***keywords  sequential sorting
c***description
c
c 1.        qpsrt
c           ordering routine
c              standard fortran subroutine
c              real version
c
c 2.        purpose
c              this routine maintains the descending ordering
c              in the list of the local error estimates resulting from
c              the interval subdivision process. at each call two error
c              estimates are inserted using the sequential search
c              method, top-down for the largest error estimate
c              and bottom-up for the smallest error estimate.
c
c 3.        calling sequence
c              call qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently
c                       in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - real
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - real
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k
c                       elements of which contain pointers
c                       to the error estimates, such that
c                       elist(iord(1)),... , elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c 4.        no subroutines or functions needed
c***end prologue  qpsrt
c
      real elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  qpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed
c           if, due to a difficult integrand, subdivision
c           increased the error estimate. in the normal case
c           the insert procedure should start after the
c           nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to
c           be maintained in descending order. this number
c           depends on the number of subdivisions still
c           allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
      real function qwgtc(x,c,p2,p3,p4,kp)
c***begin prologue  qwgtc
c***refer to qk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, cauchy principal value
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine qawc and defines the weight function.
c***end prologue  qwgtc
c
      real c,p2,p3,p4,x
      integer kp
c***first executable statement
      qwgtc = 0.1e+01/(x-c)
      return
      end
      real function qwgtf(x,omega,p2,p3,p4,integr)
c***begin prologue  qwgtf
c***refer to   qk15w
c***routines called  (none)
c***revision date 810101   (yymmdd)
c***keywords  cos or sin in weight function
c***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. * progr. div. - k.u.leuven
c***end prologue  qwgtf
c
      real omega,omx,p2,p3,p4,x
      integer integr
c***first executable statement
      omx = omega*x
      go to(10,20),integr
   10 qwgtf = cos(omx)
      go to 30
   20 qwgtf = sin(omx)
   30 return
      end
      real function qwgts(x,a,b,alfa,beta,integr)
c***begin prologue  qwgts
c***refer to qk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, algebraico-logarithmic
c             end-point singularities
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine qaws and defines the weight function.
c***end prologue  qwgts
c
      real a,alfa,b,beta,bmx,x,xma
      integer integr
c***first executable statement
      xma = x-a
      bmx = b-x
      qwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 qwgts = qwgts*alog(xma)
      go to 40
   20 qwgts = qwgts*alog(bmx)
      go to 40
   30 qwgts = qwgts*alog(xma)*alog(bmx)
   40 return
      end
