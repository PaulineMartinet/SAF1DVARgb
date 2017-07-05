!+ Invert a band diagonal matrix.

Subroutine NWPSAF_BandInverse ( &
  n,           & ! in
  p,           & ! in
  pout,        & ! in
  A,           & ! in
  AI           ) ! out

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
! Description: Invert a matrix.
!
!   Variables with intent in:
!
!     n:               Size of the matrix being inverted
!
!   Variables with intent inout:
!
!     A:               Real matrix (assumed square and symmetrical) 
!                      overwritten by its inverse
!
!   Variables with intent out:
!
!     Status:          0: ok, 1: A not positive definite.
!
! Method:
!   Uses Band Cholesky decomposition - a method particularly suitable for 
!   real symmetric band diagonal matrices.  See Golub and van Loan pp.155-6.
!
! Version Date     Comment
! ------- -------- -------
! 1.0     10/01/00 Original version using Cholesky decomposition.  
!                                                     Andrew Collard.
!                                                     Met Office
! 2.0     31/01/00 Band diagonal version.  Andrew Collard.
! 2.1     31/01/00 Code rearranged so that G is no longer calculated
!                  n times!  ADC.  
! 3.0.5   14/04/04 Remove unnecessary STOP.           Andrew Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header ---------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY : &
    StatusFatal

Implicit None

INCLUDE 'NWPSAF_Report.interface'

! Subroutine arguments
Integer, Intent(IN)  :: n            ! order of A
Integer, Intent(IN)  :: p            ! Number of bands in A
Integer, Intent(IN)  :: pout         ! Number of bands in A^-1
Real, Intent(IN)     :: A(0:p,n)     ! Input matrix
Real, Intent(OUT)    :: AI(0:pout,n) ! Inverse of A


! Local constants
CHARACTER(LEN=*), Parameter :: RoutineName = "NWPSAF_BandInverse"
CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

! local variables

INTEGER :: I, J, K, LAMBDA

REAL ::    G(0:P,N) ! The Cholesky Triangle Matrix
REAL :: Gtmp(P)
REAL ::    Q(N)
REAL ::    V(N)
REAL ::    X(0:P)   ! Temporary array used in calculating G
REAL ::   XX(N)     ! Temporary array used in forward and 
                    ! backward substituting.

!-------------------------------------------------------------------------------

! Determine the Cholesky triangle matrix.    

G(:,:)=0.
DO J = 1,N
  Lambda=MIN(P,N-J)
  x(0:Lambda) = A(0:Lambda,J)
  DO K = MAX(1,J-P),J-1
    Lambda=MIN(P,N-J,P+K-J)
    x(0:Lambda) = x(0:Lambda) - G(J-K,K)*G(J-K:Lambda+J-K,K)
  END DO
  IF (x(0) <= 0) THEN
    Errormessage(1) = 'U matrix is not positive definite'
    CALL NWPSAF_Report( &
      RoutineName,  & ! in
      ErrorMessage, & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
  Lambda=MIN(P,N-J)
  G(0:Lambda,J) = x(0:Lambda)/SQRT(x(0))
END DO 

AI(:,:)=0.

DO I=1,N

  V(:) = 0.
  V(I) = 1.

! Solve G.(xx)=v for xx by forward substitution 
  xx=v
  xx(1)=xx(1)/G(0,1)
  DO J = 2,N
    Lambda = MAX(1,J-P)
    DO K=Lambda,J-1 
      Gtmp(K-Lambda+1)=G(J-K,K)
    END DO
    xx(J) = (xx(J) - DOT_PRODUCT(Gtmp(1:J-Lambda),xx(Lambda:J-1)))/G(0,J)
  END DO

! Solve G^T.q=xx for q by backward substitution 
  q=xx
  q(N)=q(N)/G(0,N)
  DO J = N-1, 1, -1 
    Lambda = MIN(P,N-J)
    q(J) = (q(J) - DOT_PRODUCT(G(1:Lambda,J),q(J+1:J+Lambda)))/G(0,J)
  END DO

  Lambda = MIN(POUT,N-I)
  AI(0:Lambda,I) = Q(I:I+Lambda)

END DO

END Subroutine NWPSAF_BandInverse
