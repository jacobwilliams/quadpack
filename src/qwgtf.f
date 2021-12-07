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
