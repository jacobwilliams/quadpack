![quadpack2](media/logo.png)
============

![Build Status](https://github.com/jacobwilliams/quadpack/actions/workflows/CI.yml/badge.svg)

### Description

QUADPACK is a Fortran subroutine package for the numerical
computation of definite one-dimensional integrals (numerical quadrature). Development of this library, which had ceased in the 1980s, has been restarted. The original code is being modernized, and new methods are being added.  The goal is a comprehensive and modern Fortran library that includes both classic and state-of-the-art methods for numerical integration.

### Overview

The original QUADPACK code (written in the early 1980s) has been extensively refactored:

* It has been converted from FORTRAN 77 fixed form to modern free form syntax. This includes elimination of most GOTOs and other obsolescent language features.
* It is now a single stand-alone module, and has no dependencies on any other code from SLATEC or LINPACK.
* The SLATEC docstrings have been converted to [Ford](https://github.com/Fortran-FOSS-Programmers/ford) style, which allows for auto-generation of the [API docs](https://jacobwilliams.github.io/quadpack/).
* Some typos have been corrected in the comments.
* General code cleanup and formatting.
* Added automated unit testing in GitHub CI.
* The separate routines for single and double precision versions have been eliminated. The library now exports a single (`real32`), double (`real64`) and quadruple (`real128`) precision interface using the same code by employing a preprocessor scheme.
* New procedures not present in the original QUADPACK have been added.
* The coefficients have been regenerated with full quadruple precision. *(Note: this has not yet been done for all the coefficients in **DQNG**)*
* Some bugs have been fixed in the original code. Note that this version includes the recent (Oct 2021) updates (see [here](https://github.com/scipy/scipy/issues/14807) and [here](https://github.com/scipy/scipy/pull/14836)) reported by the Scipy project.

### To do list

 - [ ] Remove the last remaining GOTOs.
 - [ ] Additional docstring cleanups.
 - [ ] Add more unit tests.
 - [ ] In the unit tests, the "truth" values for the cases without analytical solutions need to be regenerated with some more precision so we have the exact results for the quad precision test.
 - [ ] Reconsider the `error stop` in `xerror` for fatal errors. Probably the library shouldn't be stopping the calling problem.

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and test cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `quadpack` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
quadpack = { git="https://github.com/jacobwilliams/quadpack.git" }
```

### Example

A simple example is given here (see the `test` folder for more examples):

```fortran
subroutine test_qag
  use quadpack, only: dqag
  use iso_fortran_env, only: wp => real64 ! double precision
  implicit none

  real(wp), parameter :: a = 0.0_wp
  real(wp), parameter :: b = 1.0_wp
  integer, parameter :: key = 6
  integer, parameter :: limit = 100
  integer, parameter :: lenw = limit*4
  real(wp), parameter :: answer = 2.0_wp/sqrt(3.0_wp)

  real(wp) :: abserr, result, work(lenw)
  integer :: ier, iwork(limit), last, neval

  call dqag(f, a, b, epsabs, epsrel, key, result, &
            abserr, neval, ier, limit, lenw, last, &
            iwork, work)

  write(*,'(1P,A,1X,*(E13.6,1X))') &
        'result, error = ', result, abs(result-answer)

contains

  real(wp) function f(x)
    implicit none
    real(wp), intent(in) :: x
    real(wp), parameter :: pi = acos(-1.0_wp)
    f = 2.0_wp/(2.0_wp + sin(10.0_wp*pi*x))
  end function f

end subroutine test_qag
```

Which outputs:

```text
 result, error = 1.154701E+00 2.220446E-16
```

### Survey of procedures

The following list gives an overview of the QUADPACK integrators.
The routine names for the double precision versions are preceded
by the letter `D`, and the quadruple precision versions are preceded by `Q`.

- **QNG**  : Is a simple non-adaptive automatic integrator, based on
         a sequence of rules with increasing degree of algebraic
         precision ([Patterson, 1968](https://www.ams.org/journals/mcom/1968-22-104/S0025-5718-68-99866-9/S0025-5718-68-99866-9.pdf)).

- **QAG**  : Is a simple globally adaptive integrator using the
         strategy of Aind (Piessens, 1973). It is possible to
         choose between 6 pairs of Gauss-Kronrod quadrature
         formulae for the rule evaluation component. The pairs
         of high degree of precision are suitable for handling
         integration difficulties due to a strongly oscillating
         integrand.

- **QAGS** : Is an integrator based on globally adaptive interval
         subdivision in connection with extrapolation ([de Doncker,
         1978](https://dl.acm.org/doi/10.1145/1053402.1053403)) by the Epsilon algorithm ([Wynn, 1956](https://www.jstor.org/stable/2002183)).

- **QAGP** : Serves the same purposes as QAGS, but also allows
         for eventual user-supplied information, i.e. the
         abscissae of internal singularities, discontinuities
         and other difficulties of the integrand function.
         The algorithm is a modification of that in QAGS.

- **QAGI** : Handles integration over infinite intervals. The
         infinite range is mapped onto a finite interval and
         then the same strategy as in QAGS is applied.

- **QAWO** : Is a routine for the integration of `COS(OMEGA*X)*F(X)`
         or `SIN(OMEGA*X)*F(X)` over a finite interval `(A,B)`.
         `OMEGA` is is specified by the user
         The rule evaluation component is based on the
         modified Clenshaw-Curtis technique.
         An adaptive subdivision scheme is used connected with
         an extrapolation procedure, which is a modification
         of that in QAGS and provides the possibility to deal
         even with singularities in F.

- **QAWF** : Calculates the Fourier cosine or Fourier sine
         transform of `F(X)`, for user-supplied interval `(A,INFINITY)`, `OMEGA`, and `F`. The procedure of QAWO is
         used on successive finite intervals, and convergence
         acceleration by means of the Epsilon algorithm ([Wynn,
         1956](https://www.jstor.org/stable/2002183)) is applied to the series of the integral
         contributions.

- **QAWS** : Integrates `W(X)*F(X)` over `(A,B)` with `A<B` finite,
         and   `W(X) = ((X-A)**ALFA)*((B-X)**BETA)*V(X)`
         where `V(X) = 1 or LOG(X-A) or LOG(B-X)`
                        or `LOG(X-A)*LOG(B-X)`
         and   `ALFA>(-1), BETA>(-1)`.
         The user specifies `A`, `B`, `ALFA`, `BETA` and the type of
         the function `V`.
         A globally adaptive subdivision strategy is applied,
         with modified Clenshaw-Curtis integration on the
         subintervals which contain `A` or `B`.

- **QAWC** : Computes the Cauchy Principal Value of `F(X)/(X-C)`
         over a finite interval `(A,B)` and for
         user-determined `C`.
         The strategy is globally adaptive, and modified
         Clenshaw-Curtis integration is used on the subranges
         which contain the point `X = C`.

   Each of the routines above also has a "more detailed" version
with a name ending in E, as QAGE.  These provide more
information and control than the easier versions.

   The preceding routines are all automatic.  That is, the user
inputs his problem and an error tolerance.  The routine
attempts to perform the integration to within the requested
absolute or relative error.
   There are, in addition, a number of non-automatic integrators.
These are most useful when the problem is such that the
user knows that a fixed rule will provide the accuracy
required.  Typically they return an error estimate but make
no attempt to satisfy any particular input error request.

  * **QK15**, **QK21**, **QK31**, **QK41**, **QK51**, **QK61**:
       Estimate the integral on [a,b] using 15, 21,..., 61
       point rule and return an error estimate.
  * **QK15I**: 15 point rule for (semi)infinite interval.
  * **QK15W**: 15 point rule for special singular weight functions.
  * **QC25C**: 25 point rule for Cauchy Principal Values
  * **QC25F**: 25 point rule for sin/cos integrand.
  * **QMOMO**: Integrates k-th degree Chebyshev polynomial times
        function with various explicit singularities.

### Other procedures

The following procedures were not in the original QUADPACK, but are included in the new library:

- **QUAD** : The result is obtained using a sequence of 1, 3, 7, 15, 31, 63, 127, and 255 point interlacing formulae. The formulae are based on the optimal extension of the 3-point gauss formula. See: [Patterson, 1968](https://www.ams.org/journals/mcom/1968-22-104/S0025-5718-68-99866-9/S0025-5718-68-99866-9.pdf). See also **QNG**. This code is based on QUAD from [NSWC Mathematical Library](https://github.com/jacobwilliams/nswc), with the addition of full quadruple-precision coefficients.

- **AVINT** : Integrates a function tabulated at arbitrarily spaced abscissas using overlapping parabolas. This procedure was originally from [SLATEC](http://www.netlib.org/slatec/src/davint.f).

- **QNC79** :  Integrate a function over a finite interval using a 7-point adaptive Newton-Cotes quadrature rule. This procedure was originally from [SLATEC](http://www.netlib.org/slatec/src/dqnc79.f).

- **GAUSS8** : Integrate a function over a finite interval using an adaptive 8-point Legendre-Gauss algorithm. This procedure was originally from [SLATEC](http://www.netlib.org/slatec/src/dgauss8.f).

- **SIMPSON** : Integrate a function over a finite interval using an adaptive Simpson rule. See: [Gander & Gautschi, 2000](https://www.researchgate.net/publication/226706221_Adaptive_Quadrature-Revisited).

- **LOBATTO** : Integrate a function over a finite interval using an adaptive Lobatto rule. See: [Gander & Gautschi, 2000](https://www.researchgate.net/publication/226706221_Adaptive_Quadrature-Revisited).

### Guidelines for the use of QUADPACK

Here it is not our purpose to investigate the question when
automatic quadrature should be used. We shall rather attempt
to help the user who already made the decision to use QUADPACK,
with selecting an appropriate routine or a combination of
several routines for handling his problem.

For both quadrature over finite and over infinite intervals,
one of the first questions to be answered by the user is
related to the amount of computer time he wants to spend,
versus his -own- time which would be needed, for example, for
manual subdivision of the interval or other analytic
manipulations.

1.  The user may not care about computer time, or not be
willing to do any analysis of the problem. especially when
only one or a few integrals must be calculated, this attitude
can be perfectly reasonable. In this case it is clear that
either the most sophisticated of the routines for finite
intervals, QAGS, must be used, or its analogue for infinite
intervals, GAGI. These routines are able to cope with
rather difficult, even with improper integrals.
This way of proceeding may be expensive. But the integrator
is supposed to give you an answer in return, with additional
information in the case of a failure, through its error
estimate and flag. Yet it must be stressed that the programs
cannot be totally reliable.

2. The user may want to examine the integrand function.
If bad local difficulties occur, such as a discontinuity, a
singularity, derivative singularity or high peak at one or
more points within the interval, the first advice is to
split up the interval at these points. The integrand must
then be examined over each of the subintervals separately,
so that a suitable integrator can be selected for each of
them. If this yields problems involving relative accuracies
to be imposed on -finite- subintervals, one can make use of
QAGP, which must be provided with the positions of the local
difficulties. However, if strong singularities are present
and a high accuracy is requested, application of QAGS on the
subintervals may yield a better result.

For quadrature over finite intervals we thus dispose of QAGS
and
- QNG for well-behaved integrands,
- QAG for functions with an oscillating behaviour of a non
    specific type,
- QAWO for functions, eventually singular, containing a
    factor `COS(OMEGA*X)` or `SIN(OMEGA*X)` where OMEGA is known,
- QAWS for integrands with Algebraico-Logarithmic end point
    singularities of known type,
- QAWC for Cauchy Principal Values.

### Remark

On return, the work arrays in the argument lists of the
adaptive integrators contain information about the interval
subdivision process and hence about the integrand behaviour:
the end points of the subintervals, the local integral
contributions and error estimates, and eventually other
characteristics. For this reason, and because of its simple
globally adaptive nature, the routine QAG in particular is
well-suited for integrand examination. Difficult spots can
be located by investigating the error estimates on the
subintervals.

For infinite intervals we provide only one general-purpose
routine, QAGI. It is based on the QAGS algorithm applied
after a transformation of the original interval into (0,1).
Yet it may eventuate that another type of transformation is
more appropriate, or one might prefer to break up the
original interval and use QAGI only on the infinite part
and so on. These kinds of actions suggest a combined use of
different QUADPACK integrators. Note that, when the only
difficulty is an integrand singularity at the finite
integration limit, it will in general not be necessary to
break up the interval, as QAGI deals with several types of
singularity at the boundary point of the integration range.
It also handles slowly convergent improper integrals, on
the condition that the integrand does not oscillate over
the entire infinite interval. If it does we would advise
to sum succeeding positive and negative contributions to
the integral -e.g. integrate between the zeros- with one
or more of the finite-range integrators, and apply
convergence acceleration eventually by means of QUADPACK
subroutine QELG which implements the Epsilon algorithm.
Such quadrature problems include the Fourier transform as
a special case. Yet for the latter we have an automatic
integrator available, QAWF.

### Documentation

The API documentation for the current `master` branch can be found [here](https://jacobwilliams.github.io/quadpack/).  This is generated by processing the source files with [FORD](https://github.com/Fortran-FOSS-Programmers/ford). Note that the procedures listed in the API documentation are the double precision version (`DQNG`, etc.)

### License

The original Quadpack was a public domain work of the United States government. The modifications are released under a permissive (BSD-3) license.

### References

  * R. Piessens, E. deDoncker-Kapenga, C. Uberhuber, D. Kahaner
    [Quadpack: a Subroutine Package for Automatic Integration](https://link.springer.com/book/10.1007/978-3-642-61786-7)
    Springer Verlag, 1983. Series in Computational Mathematics v.1
    515.43/Q1S 100394Z
  * Paola Favati, Grazia Lotti, Francesco Romani, [Algorithm 691: Improving QUADPACK automatic integration routines](https://dl.acm.org/doi/abs/10.1145/108556.108580), ACM Transactions on Mathematical Software, Volume 17, Issue 2, June 1991, pp 218-232.
  * Original SLATEC code from [Netlib](http://www.netlib.org/quadpack/). Last modified 11 Oct 2021.
  * W. Gander and W. Gautschi, "[Adaptive Quadrature - Revisited](https://www.researchgate.net/publication/226706221_Adaptive_Quadrature-Revisited)", BIT Vol. 40, No. 1, March 2000, pp. 84--101.

### Other versions

There are other versions of Quadpack out there. There are at least two projects to provide module interface to the unmodified Fortran 77 code (see [nshaffer/modern_quadpack](https://github.com/nshaffer/modern_quadpack) and [ivan-pi/quadpack](https://github.com/ivan-pi/quadpack)). The license for these are not specified.  Another fixed to free conversion can be found at [John Burkardt's site](https://people.math.sc.edu/Burkardt/f_src/quadpack_double/quadpack_double.f90) (this is not an aggressive modernization though and also has an LGPL license). Also note that the Quadpack code in [SLATEC](http://www.netlib.org/slatec/src/) is slightly modified from the stand-alone one at [Netlib](http://www.netlib.org/quadpack/). It is not known if these modifications were anything significant.

### See also

 * [quadrature-fortran](https://github.com/jacobwilliams/quadrature-fortran)

### Keywords
  * survey of integrators, guidelines for selection,
    quadpack, automatic integrator, general-purpose,
    integrand examinator, globally adaptive,
    gauss-kronrod, infinite intervals, transformation,
    extrapolation, singularities at user specified points,
    (end-point) singularities, cauchy principal value,
    clenshaw-curtis method, special-purpose, fourier integral,
    integration between zeros, convergence acceleration,
    integrand with oscillatory cos or sin factor,
    (end point) singularities, 25-point clenshaw-curtis integration,
    smooth integrand, non-adaptive, gauss-kronrod (patterson),
    epsilon algorithm, algebraico-logarithmic end point singularities,
    chebyshev series expansion, fast fourier transform
