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
