#ifndef MOD_INCLUDE
    program quadpack_tests
    use quadpack_test_module
    implicit none
    write(*,*) ''
    write(*,*) ' quadpack_tests : Default (kind=',wp,')'
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
    call test_C()
    !call test_G()
    call test_davint()

#ifndef MOD_INCLUDE
    end program quadpack_tests
#endif