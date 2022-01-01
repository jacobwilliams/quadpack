#ifndef MOD_INCLUDE
    program quadpack_tests
    use quadpack_test_module
    implicit none
    write(*,*) ''
    write(*,*) ' quadpack_tests : Default'
    write(*,*) ''
#endif

    call test_qag()
    call test_qagi()
    call test_qagp()
    call test_qags()
    call test_qawc()
    call test_qawf()
    call test_qawo()
    call test_qaws()
    call test_qng()

#ifndef MOD_INCLUDE
    end program quadpack_tests
#endif