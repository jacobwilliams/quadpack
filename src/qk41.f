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
