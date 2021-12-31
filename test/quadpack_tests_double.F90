program quadpack_tests_double
    use quadpack_double, wp => quadpack_RK
    use quadpack_test_module_double
    write(*,*) ''
    write(*,*) ' quadpack_tests : Double'
    write(*,*) ''
#define MOD_INCLUDE=1
#include "quadpack_tests.F90"
    end program quadpack_tests_double
