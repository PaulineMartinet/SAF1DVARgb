SUBROUTINE NWPSAF_CHOLESKY(U, &
     V, &
     N, &
     q  )

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
!
! Solves the Linear equation Uq=V for q where U is a symmetric
! positive definite matrix and v and q are vectors of length N.
!  
! The method follows that in Golub and Van Loan although this is
! pretty standard. 
! 
! If U is not positive definite this will be detected by the program
! and flagged as an error.  U is assumed to be symmetric as only the
! upper triangle is in fact used.
!
! History:
!   1.0   24Nov99  Program written.  Andrew Collard.
! 

USE NWPSAFMod_Params, ONLY : &
    StatusFatal

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'

INTEGER, INTENT(IN)  :: N  
REAL, INTENT(IN)     :: U(N,N)    
REAL, INTENT(IN)     :: v(N)
REAL, INTENT(OUT)    :: q(N)      ! This is the result

CHARACTER(LEN=*), PARAMETER :: RoutineName = "NWPSAF_Cholesky"
CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

INTEGER :: J, K

REAL :: G(N,N)   ! The Cholesky Triangle Matrix
REAL :: x(N)     ! Temporary array used in calculating G
                 ! and in forward and backward substituting.
!-----------------------------------------------------------------------------

! Determine the Cholesky triangle matrix.    

DO J = 1,N
  x(J:N) = U(J:N,J)
  IF (J /= 1) THEN
    DO K = 1,J-1
      x(J:N) = x(J:N) - G(J,K)*G(J:N,K)
    END DO
  END IF
  IF (x(J) <= 0) THEN
    Errormessage(1) = 'U matrix is not positive definite'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusFatal )
  ENDIF
  G(J:N,J) = x(J:N)/SQRT(x(J))
END DO

! Solve Gx=v for x by forward substitution 

x=v
x(1)=x(1)/G(1,1)
DO J = 2,N
  x(J) = (x(J) - DOT_PRODUCT(G(J,1:J-1),x(1:J-1)))/G(J,J)
END DO

! Solve G^T.q=x for q by backward substitution 

q=x
q(N)=q(N)/G(N,N)
DO J = N-1, 1, -1 
  q(J) = (q(J) - DOT_PRODUCT(G(J+1:N,J),q(J+1:N)))/G(J,J)
END DO

END SUBROUTINE NWPSAF_CHOLESKY
