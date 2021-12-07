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
