module quadpack_quad
    !!
    !!@note For this module, `wp` is `real128` (double precision).
    !!
    use iso_fortran_env, only: wp => real128
#if !defined(NOQUAD)
#define MOD_INCLUDE=1
#include "quadpack_generic.F90"
#endif
end module quadpack_quad