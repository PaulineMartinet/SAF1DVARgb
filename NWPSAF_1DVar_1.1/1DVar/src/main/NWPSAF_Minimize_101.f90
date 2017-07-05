Subroutine NWPSAF_Minimize_101( &
     Bmatrix,          & ! in
     DeltaBT,          & ! in  
     Num_Chans,        & ! in
     HTX_Dimension,    & ! in
     H_Matrix,         & ! in
     H_Matrix_T,       & ! in
     R_SubMatrix,      & ! in
     Delta_Profile,    & ! inout
     Status)             ! out

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
! Description:
!   
!   Updates the profile vector Delta_Profile according to Rodgers (1976), 
!   Eqn. 101:
!
!   x_(n+1) = xb + B.Hn'.Q
!   Q = U^-1.V
!   where: x is an atmospheric state vector, subscripted b=background,n=nth
!           iteration
!          U = (Hn.B.Hn'+R)
!          V = (ym-y(xn) - H'.(xb-xn))
!          B is the background error covariance matrix
!          R is the combined forward model and ob error covariance matrix
!
!   Q = U^-1.V is solved by Cholesky decomposition.
!
!   This routine should be used when:
!     1) The length of the observation vector is less than the length 
!        of the state vector, 
!     2) Where no additional cost function terms are provided and 
!     3) where Newtonian minimisation is desired.
!
! References: 
!
!   Rodgers, Retrieval of atmospheric temperature and composition from
!   remote measurements of thermal radiation, Rev. Geophys.Sp.Phys. 14, 1976.
!
!   Rodgers, Inverse Methods for Atmospheres: Theory and Practice.  World 
!            Scientific Publishing, 2000.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     22/02/00 Code for calculation of U and V split off from main
!                  minimisation routine.     A. Collard.  Met. Office.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type, &
    R_Full_Matrix, &
    R_Band_Diagonal, & 
    R_Eigenvectors

USE NWPSAFMod_RTModel, ONLY : &
    Num_ProfElementsUsed

USE NWPSAFMod_Params, ONLY : &
    StatusOK,                &
    StatusFatal

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_Cholesky.interface'

! Subroutine arguments:
INTEGER, INTENT(IN) :: Num_Chans         ! No. of channels to be used 
INTEGER, INTENT(IN) :: HTX_Dimension     ! Dimension of HTX Array
REAL, INTENT(IN) :: Bmatrix(:,:)         ! B-matrix
REAL, INTENT(IN) :: DeltaBT(:)           ! y-y(x)
REAL, INTENT(INOUT) :: Delta_Profile(:)  !profile increments 
REAL, INTENT(IN) :: H_Matrix(:,:)        ! Jacobian
REAL, INTENT(IN) :: H_Matrix_T(:,:)      ! (Jacobian)^T
TYPE(R_Matrix_Type), INTENT(IN) :: R_SubMatrix  ! R Matrix for selected chans.
INTEGER, INTENT(OUT) :: Status

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Minimize_101"

! Local variables:
CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

INTEGER :: Element, I

REAL :: HB(Num_Chans, Num_ProfElementsUsed) ! Scratch vector
REAL :: HTX(Num_Chans, HTX_Dimension) ! Scratch vector (eigenvectors case only)
REAL :: Q(Num_Chans)                  ! Q = U^-1.V
REAL :: U (Num_Chans, Num_Chans)      ! U = H.B.H^T + R
REAL :: V (Num_Chans)                 ! V = (y-y(x_n))-H^T(xb-x_n)

!-----------------------------------------------------------------------------

Status=StatusOK

!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.  
!---------------------------------------------------------------------------

HB = MATMUL(H_matrix,Bmatrix)
U = MATMUL(HB,H_matrix_T)

V = DeltaBT + MATMUL(H_matrix,Delta_Profile)

!---------------------------------------------------------------------------
! 1.1. Add the R matrix into the U matrix.
!---------------------------------------------------------------------------

SELECT CASE (R_SubMatrix % RType)

  CASE (R_Full_Matrix)
    U = U + R_SubMatrix % Matrix(:,:)

  CASE (R_Band_Diagonal)
    DO Element = 0,R_SubMatrix % Num_Elements
      DO I=1,Num_Chans-Element
        U(I,I+Element) = U(I,I+Element) + R_SubMatrix % Matrix(Element,I)
        U(I+Element,I) = U(I,I+Element) 
      END DO
    END DO

  CASE (R_Eigenvectors)
    ! The automatic array HTX should have a second dimension 
    ! R_Matrix % Num_Elements.  Check that this is the case.
    IF (HTX_Dimension /= R_SubMatrix % Num_Elements) THEN
      Errormessage(1) = 'Incorrect value for HTX_Dimension'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in    
    END IF
    DO I=1,R_SubMatrix % Num_Elements
      HTX(:,I)=R_SubMatrix % EigenValues(I) * R_SubMatrix % Matrix(:,I)
    END DO
    U = U + MATMUL(TRANSPOSE(R_SubMatrix % Matrix(:,:)),HTX)

  CASE DEFAULT
    Errormessage(1) = 'Unknown Format for Observation Errors'
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 

END SELECT

! Calculate Q=(U^-1).V
!------
CALL NWPSAF_CHOLESKY( &
  U,         & !in
  V,         & !in
  Num_Chans, & !in
  Q          ) !out


! Delta profile is (HB)^T.Q
!------
Delta_Profile = MATMUL(TRANSPOSE(HB), Q)

End Subroutine NWPSAF_Minimize_101
