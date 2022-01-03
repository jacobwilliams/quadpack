module quadpack_test_module_double
    use quadpack_double, wp => quadpack_RK
#define MOD_INCLUDE=1
#include "quadpack_test_module.F90"
end module quadpack_test_module_double
