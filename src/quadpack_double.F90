module quadpack_double
    !!
    !!@note For this module, `wp` is `real64` (double precision).
    !!
use iso_fortran_env, only: wp => real64
#define MOD_INCLUDE=1
#include "quadpack.F90"
end module quadpack_double