program quadpack_tests_single
    use quadpack_single, wp => quadpack_RK
    use quadpack_test_module_single
    write(*,*) ''
    write(*,*) ' quadpack_tests : Single'
    write(*,*) ''
#define MOD_INCLUDE=1
#include "quadpack_tests.F90"
end program quadpack_tests_single
