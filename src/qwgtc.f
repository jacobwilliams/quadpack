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
