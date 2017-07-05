!+ Routine to handle error reports

SUBROUTINE NWPSAF_Report( &
  NameOfRoutine,          & ! in
  Message,                & ! in
  ErrorStatus)              ! in

! Description:
!   Writes out fatal, warning and info error messages. Execution is terminated
!   in the event of a fatal error.
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0       07/11/16 Simplified from Gen code. Fiona Smith
!
! Code Description:
!   Language:       Fortran 90.
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------
!
! Modules used:

USE NWPSAFMod_Params, ONLY : &
    StatusFatal,        &
    StatusWarning

IMPLICIT NONE

!* Subroutine arguments
CHARACTER(LEN=*), INTENT(IN)  :: NameOfRoutine ! Calling this one
CHARACTER(LEN=*), INTENT(IN)  :: Message(:)    ! Message to output
INTEGER, INTENT(IN), OPTIONAL :: ErrorStatus   ! Input error code

! Local variables
INTEGER                       :: lines         ! number of message lines
INTEGER                       :: line
INTEGER                       :: ErrorReport  
!------------------------------------------------------------------------------

lines = SIZE(Message, 1)

IF( present(ErrorStatus)) THEN
  ErrorReport=ErrorStatus
ELSE
  ErrorReport = 0
ENDIF

IF (ErrorReport == StatusWarning) THEN
  WRITE(*,'(A)') "WARNING: " // TRIM(NameOfRoutine)
  DO line = 1, lines
    IF (Message(line) /= "") WRITE(*,'(A)') TRIM(Message(line))
  END DO
ELSE IF (ErrorReport  == StatusFatal ) THEN
  WRITE(*,'(A)') "FATAL ERROR: " // TRIM(NameOfRoutine)
  DO line = 1, lines
    IF (Message(line) /= "") WRITE(*,'(A)') TRIM(Message(line))
  END DO
  STOP
ELSE ! Anything else is OK, just print message
  WRITE(*,'(A)') TRIM(NameOfRoutine)
  DO line = 1, lines
    IF (Message(line) /= "") WRITE(*,'(A)') TRIM(Message(line))
  END DO
END IF

END SUBROUTINE NWPSAF_Report
