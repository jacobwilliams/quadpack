program quadpack_tests_quad
#if !defined(NOQUAD)
    use quadpack_quad, wp => quadpack_RK
    use quadpack_test_module_quad
    write(*,*) ''
    write(*,*) ' quadpack_tests : quad'
    write(*,*) ''
#define MOD_INCLUDE=1
#include "quadpack_tests.F90"
#endif
end program quadpack_tests_quad
