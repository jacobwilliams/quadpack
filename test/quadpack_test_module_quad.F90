module quadpack_test_module_quad
#if !defined(NOQUAD)
    use quadpack_quad, wp => quadpack_RK
#if !defined(__GFORTRAN__)
use quadpack_quad, dqag => qqag, &
                   dqage => qqage, &
                   dqagi => qqagi, &
                   dqagie => qqagie, &
                   dqagp => qqagp, &
                   dqagpe => qqagpe, &
                   dqags => qqags, &
                   dqagse => qqagse, &
                   dqawc => qqawc, &
                   dqawce => qqawce, &
                   dqawf => qqawf, &
                   dqawfe => qqawfe, &
                   dqawo => qqawo, &
                   dqawoe => qqawoe, &
                   dqaws => qqaws, &
                   dqawse => qqawse, &
                   dqc25c => qqc25c, &
                   dqc25f => qqc25f, &
                   dqc25s => qqc25s, &
                   dqcheb => qqcheb, &
                   dqk15 => qqk15, &
                   dqk15i => qqk15i, &
                   dqk15w => qqk15w, &
                   dqk21 => qqk21, &
                   dqk31 => qqk31, &
                   dqk41 => qqk41, &
                   dqk51 => qqk51, &
                   dqk61 => qqk61, &
                   dqmomo => qqmomo, &
                   dqng => qqng
#endif
#define MOD_INCLUDE=1
#include "quadpack_test_module.F90"
#endif
    end module quadpack_test_module_quad
