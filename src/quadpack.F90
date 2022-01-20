!********************************************************************************
!>
!  Main QUADPACK module.
!
!  This module simply imports the specific modules for
!  single, double, and quad precision.

    module quadpack

    use quadpack_single, only: qag => dqag, qage => dqage, qagi => dqagi, qagie => dqagie, qagp => dqagp, &
                               qagpe => dqagpe, qags => dqags, qagse => dqagse, qawc => dqawc, &
                               qawce => dqawce, qawf => dqawf, qawfe => dqawfe, qawo => dqawo, &
                               qawoe => dqawoe, qaws => dqaws, qawse => dqawse, qc25c => dqc25c, &
                               qc25f => dqc25f, qc25s => dqc25s, qcheb => dqcheb, qk15 => dqk15, &
                               qk15i => dqk15i, qk15w => dqk15w, qk21 => dqk21, qk31 => dqk31, &
                               qk41 => dqk41, qk51 => dqk51, qk61 => dqk61, qmomo => dqmomo, qng => dqng, &
                               quad => dquad, &
                               avint => davint, &
                               qnc79 => dqnc79, &
                               gauss8 => dgauss8

    use quadpack_double

#if !defined(NOQUAD)
    use quadpack_quad, only: qqag => dqag, qqage => dqage, qqagi => dqagi, qqagie => dqagie, qqagp => dqagp, &
                             qqagpe => dqagpe, qqags => dqags, qqagse => dqagse, qqawc => dqawc, &
                             qqawce => dqawce, qqawf => dqawf, qqawfe => dqawfe, qqawo => dqawo, &
                             qqawoe => dqawoe, qqaws => dqaws, qqawse => dqawse, qqc25c => dqc25c, &
                             qqc25f => dqc25f, qqc25s => dqc25s, qqcheb => dqcheb, qqk15 => dqk15, &
                             qqk15i => dqk15i, qqk15w => dqk15w, qqk21 => dqk21, qqk31 => dqk31, &
                             qqk41 => dqk41, qqk51 => dqk51, qqk61 => dqk61, qqmomo => dqmomo, qqng => dqng, &
                             qquad => dquad, &
                             qavint => davint, &
                             qqnc79 => dqnc79, &
                             qgauss8 => dgauss8
#endif

!********************************************************************************
    end module quadpack
!********************************************************************************