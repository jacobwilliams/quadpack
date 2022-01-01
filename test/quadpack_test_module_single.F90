module quadpack_test_module_single
    use quadpack_single, wp => quadpack_RK, &
                         dqag => qag, &
                         dqage => qage, &
                         dqagi => qagi, &
                         dqagie => qagie, &
                         dqagp => qagp, &
                         dqagpe => qagpe, &
                         dqags => qags, &
                         dqagse => qagse, &
                         dqawc => qawc, &
                         dqawce => qawce, &
                         dqawf => qawf, &
                         dqawfe => qawfe, &
                         dqawo => qawo, &
                         dqawoe => qawoe, &
                         dqaws => qaws, &
                         dqawse => qawse, &
                         dqc25c => qc25c, &
                         dqc25f => qc25f, &
                         dqc25s => qc25s, &
                         dqcheb => qcheb, &
                         dqk15 => qk15, &
                         dqk15i => qk15i, &
                         dqk15w => qk15w, &
                         dqk21 => qk21, &
                         dqk31 => qk31, &
                         dqk41 => qk41, &
                         dqk51 => qk51, &
                         dqk61 => qk61, &
                         dqmomo => qmomo, &
                         dqng => qng
#define MOD_INCLUDE=1
#include "quadpack_test_module.F90"
    end module quadpack_test_module_single
