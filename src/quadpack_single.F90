module quadpack_single
    !!
    !!@note For this module, `wp` is `real32` (single precision).
    !!
    use iso_fortran_env, only: wp => real32
#define MOD_INCLUDE=1
#include "quadpack_generic.F90"
end module quadpack_single