module quadpack_test_module_quad
#if !defined(NOQUAD)
    use quadpack_quad, wp => quadpack_RK
#define MOD_INCLUDE=1
#include "quadpack_test_module.F90"
#endif
end module quadpack_test_module_quad
