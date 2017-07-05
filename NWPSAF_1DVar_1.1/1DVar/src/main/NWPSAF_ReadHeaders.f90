SUBROUTINE NWPSAF_ReadHeaders ( &
    NumVars,       & ! in
    String,        & ! inout
    Vars,          & ! out
    Status)          ! inout

!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2004, EUMETSAT, All Rights Reserved.
!
!------------------------------------------------------------------------------
! Description: Read in colon-separated fields in observation and
!              background file headers.
!
!
! Method:
!     Obvious.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 2.2      30/04/02  Original Version
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

IMPLICIT NONE

! Subroutine Arguments:

INTEGER          , INTENT(IN)    :: NumVars
CHARACTER (LEN=*), INTENT(INOUT) :: String
REAL             , INTENT(OUT)   :: Vars(NumVars)
INTEGER          , INTENT(OUT)   :: Status

! Internal Variables:

INTEGER :: I
INTEGER :: II

! Program:

Status = 0
Vars(:) = 0.

DO I=1,NumVars
  II=Index(String,':')
  IF (II == 0) Then
     Status = -1
     EXIT
  END IF
  String=String(II+1:)
  READ(String,*,IOSTAT=Status) Vars(I)
  IF (Status /= 0) EXIT
END DO


END Subroutine NWPSAF_ReadHeaders


