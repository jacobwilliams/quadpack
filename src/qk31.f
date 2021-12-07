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
