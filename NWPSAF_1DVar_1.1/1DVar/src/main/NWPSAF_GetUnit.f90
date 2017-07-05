!+ Supplies an unused fortran I/O unit number.

SUBROUTINE NWPSAF_GetUnit( &
  unit,              & ! out
  ErrorStatus )        ! Out  

! Description:
!   Supplies an unused fortran I/O unit number.
!
!   If a unit number cant be supplied then NWPSAF_Report is called
!   with a warning.
!
! Method:
!   The module GenMod_FortranIO keeps track of which unit numbers have been 
!   assigned.
!
! Owner:John Bray
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       21/07/97 Creation from VarMod_FortranIO 1.8 JB
! 1.2       04/09/97 ErrorStatus INTENT(OUT)
! 1.5       16/02/98 ErrorStatus OK if no problem. DPM
!
! Code Description:
!   Language:       Fortran 90.
!   Software Standards: GTDP 8
!
! Parent module: GenMod_FortranIO
!
! Declarations:

! Modules used:

USE NWPSAFMod_Params, ONLY :       &
    LB, &
    UB, &
    UnitsInUse,  &
    StatusOK,                           & ! Good ErrorStatus value.
    StatusWarning                         ! Warning error ErrorStatus value.

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'

! Subroutine arguments

INTEGER, INTENT(OUT) :: unit        ! Unit number supplied
INTEGER, INTENT(OUT) :: ErrorStatus ! Error status flag

! Local constants

CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_GetUnit"

!- End of header

unit = LB

DO WHILE (UnitsInUse(unit) .AND. (unit < UB) )
  unit = unit + 1
END DO

IF (unit == UB) THEN
  ! No unassigned unit numbers found. Report a warning message.
  unit        = 0
  ErrorStatus = StatusWarning

  CALL NWPSAF_Report(                                        &
    RoutineName,                                             &
    (/ 'No unused unit numbers available for fortran IO' /), &
    ErrorStatus=ErrorStatus)
ELSE
  UnitsInUse(unit) = .TRUE.
  ErrorStatus = StatusOK
END IF

END SUBROUTINE NWPSAF_GetUnit
