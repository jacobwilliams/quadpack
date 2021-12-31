module quadpack_test_module_quad
    use quadpack_quad, wp => quadpack_RK
#define MOD_INCLUDE=1
#include "quadpack_test_module.F90"
    end module quadpack_test_module_quad
