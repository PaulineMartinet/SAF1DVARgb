!+ Invert a matrix.

Subroutine NWPSAF_SatMatInv ( &
  N,           & ! in
  M,           & ! in
  A,           & ! inout
  Matrix)        ! optional inout     

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
!     N:  Size of the matrix being inverted
!     M:  If MATRIX is not present this is the same as N, else this is
!         the other dimension of MATRIX.
!
!   Variables with intent inout:
!
!     A:               Real matrix (assumed square and symmetrical) 
!                      overwritten by its inverse if MATRIX is not 
!                      present.
!
!   Variables with intent out:
!
!     Status:          0: ok, 1: A is not positive definite.
!
!   Optional variables with intent inout:
!
!     Matrix:          If present this input matrix is replaced by
!                      (Matrix).A^-1 on exit (leaving A unchanged). 
!
!
! Method:
!   Uses Cholesky decomposition - a method particularly suitable for real    
!   symmetric matrices.  
!
!   Cholesky decomposition solves the Linear equation UQ=V for Q 
!   where U is a symmetric positive definite matrix and U and Q 
!   are vectors of length N.
!  
!   The method follows that in Golub and Van Loan although this is
!   pretty standard. 
! 
!   If U is not positive definite this will be detected by the program
!   and flagged as an error.  U is assumed to be symmetric as only the
!   upper triangle is in fact used.
!
!   If the the optional parameter Matrix is present, it is replaced by 
!   (Matrix).A^-1 on exit and A is left unchanged.
!
! Version Date     Comment
! ------- -------- -------
! 1.0     10/01/00 Original version using Cholesky decomposition.  
!                                                       Andrew Collard
!                                                           Met Office
! 1.1     31/01/00 Code rearranged so that G is no longer calulated
!                  n times!  ADC.
! 1.2     22/02/00 Added option to calculate (Matrix)^T.A^-1.  ADC.
! 1.3     02/02/01 Reversed indices on MATRIX.  A. Collard.
! 3.0.5   23/03/04 Fixed MATRIX indices order.  A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY : &
    StatusFatal

Implicit None

INCLUDE 'NWPSAF_Report.interface'

! Subroutine arguments
Integer, Intent(IN)    :: n          ! order of A
Integer, Intent(IN)    :: m          ! Order of optional matrix, if required
Real, Intent(INOUT)    :: A(n,n)     ! square mx, overwritten by its inverse
REAL, OPTIONAL, INTENT(INOUT) :: Matrix(N,M) 
    ! If Matrix is present on entry, it is replaced by (Matrix).A^-1 on exit 


! Local constants
CHARACTER(LEN=*), Parameter :: RoutineName = "NWPSAF_SatMatInv"
CHARACTER(LEN=80)  :: ErrorMessage(2) ! Message for NWPSAF_Report

! local variables

INTEGER :: I, J, K, MM 

REAL :: G(N,N)   ! The Cholesky Triangle Matrix
REAL :: Q(N)
REAL :: TMP(N,M)
REAL :: V(N)
REAL :: x(N)     ! Temporary array used in calculating G
                 ! and in forward and backward substituting.

!-----------------------------------------------------------------------------

ErrorMessage = ' '

! Determine the Cholesky triangle matrix.    

DO J = 1,N
  x(J:N) = A(J:N,J)
  IF (J /= 1) THEN
    DO K = 1,J-1
      x(J:N) = x(J:N) - G(J,K)*G(J:N,K)
    END DO
  END IF
  IF (x(J) <= 0) THEN
    Errormessage(1) = 'Matrix is not positive definite'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusFatal )   
  END IF
  G(J:N,J) = x(J:N)/SQRT(x(J))
END DO 
 
! Now solve the equation G.G^T.q=v for the set of
! vectors, v, with one element = 1 and the rest zero.
! The solutions q are brought together at the end to form
! the inverse of G.G^T (i.e., the inverse of A).

IF (PRESENT(Matrix)) THEN
  MM=M
ELSE
  ! Make sure that the dimensions of TMP were correctly specified
  IF (M /= N) THEN
    Errormessage(1) = '2nd and 3rd Arguments of routine must be'
    Errormessage(2) = 'identical if the MATRIX option is not present'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusFatal )   
  END IF
  MM=N
END IF

Main_Loop : DO I=1,MM 
  IF (.NOT.(PRESENT(Matrix))) THEN
    V(:) = 0.
    V(I) = 1.
  ELSE
    V(:) = Matrix(:,I)
  END IF

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

   TMP(:,I) = Q(:)
END DO Main_Loop

IF (.NOT.(PRESENT(Matrix))) THEN
  A(:,:)=TMP(:,:)
ELSE
  Matrix(:,:)=TMP(:,:)
END IF
   
END Subroutine NWPSAF_SatMatInv
