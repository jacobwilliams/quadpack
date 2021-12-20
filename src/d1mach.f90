FUNCTION d1mach(i)
!!!-------------------------------------------------------------------
!!!
!!! This function is  intended to replace  the old D1MACH by using F90
!!! intrinsic functions.
!!!
!!!
!!! The traditional D1MACH constants are ...
!!!
!!!
!!! -- DOUBLE-PRECISION MACHINE CONSTANTS --
!!!
!!! D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!!!
!!! D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!!!
!!! D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!!!
!!! D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!!!
!!! D1MACH( 5) = LOG10(B)
!!!
!!!-------------------------------------------------------------------
! .. Implicit None Statement ..
  IMPLICIT NONE
! ..
! .. Function Return Value ..
  DOUBLE PRECISION :: d1mach
! ..
! .. Scalar Arguments ..
  INTEGER :: i
! ..
! .. Local Scalars ..
  LOGICAL, SAVE :: qfirst_call = .TRUE.
! ..
! .. Local Arrays ..
  DOUBLE PRECISION, SAVE :: d1mach_values(5)
! ..
! .. Intrinsic Functions ..
  INTRINSIC digits, epsilon, huge, kind, log10, radix, real, tiny
! ..
! .. Executable Statements ..
  IF (i<1 .OR. i>5) THEN
    WRITE (*,'(1x,''D1MACH(I) - I out of bounds, I ='',i10)') i
    STOP ' D1MACH(I) - I out of bounds'
  END IF
  IF (qfirst_call) THEN
    d1mach_values = (/ tiny(1.0D0), huge(1.0D0), &
      real(radix(1.0D0),kind(1.0D0))**(-digits(1.0D0)), &
      epsilon(1.0D0), log10(real(radix(1.0D0),kind(1.0D0))) /)
    qfirst_call = .FALSE.
  END IF
  d1mach = d1mach_values(i)
  RETURN
END FUNCTION d1mach
