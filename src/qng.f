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
