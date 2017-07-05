SUBROUTINE NWPSAF_Calculate_Cost_Function( &
     BMatrixInverse,   &  ! in
     R_SubMatrix,      &  ! inout
     Delta_Profile,    &  ! in
     DeltaObs,         &  ! in
     Work_Dimension,   &  ! in
     Eqn_101,          &  ! in
     Guess_Profile,    &  ! in
     JCost,            &  ! out
     Status,           &  ! out 
     H_Matrix_T,       &  ! optional in 
     JCost_Gradient)      ! optional out

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
!   Calculates the cost function.  Contained in NWPSAF_Minimise
!
! History:
!   1.0   15Dec99  Stripped out of NWPSAF_Minimise.  Andrew Collard.
!   1.1   28Feb00  Rewritten as a stand alone subroutine rather than being
!                  simply CONTAINed in NWPSAF_Minimise.
!   1.2   04Feb00  When Rogers's Eqn. 101 is being used, the inverse of 
!                  the observational error covariance matrix is only 
!                  needed in calculating the cost function which in turn is 
!                  only needed for a convergence test.  Rather than do a
!                  largely unneccessary inversion, here we simply assume
!                  the R-matrix is diagonal if the logical argument Eqn_101
!                  is .TRUE. 
!   1.3   22Nov00  Fixed a bug in the evaluation of Scratch_VectorO in the 
!                  band diagonal (not diagonal) case.  A. Collard.
!   2.2   14Mar02  Fixed calculation of JCost_Gradient. 
!                  Fixed calculation of Jo in Eqn_101, full matrix case.
!   2.3   23May02  Changed H_Matrix_T to optional argument.
!   3.2   04Feb04  Add in Additional Cost Function.  A. Collard.
!

USE NWPSAFMod_Params, ONLY :    &
    StatusFatal,              &
    StatusWarning,            &
    Additional_Cost_Function, &
    No_Additional_Cost_Function, &
    CloudBoundaries

USE NWPSAFMod_RTModel, ONLY :    &
    Num_ProfElementsUsed

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type, &
    R_Full_Matrix, &
    R_Band_Diagonal, & 
    R_Eigenvectors, &
    Max_Elements

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_SatMatInv.interface'
INCLUDE 'NWPSAF_BandInverse.interface'
INCLUDE 'NWPSAF_BandMultiply.interface'
INCLUDE 'NWPSAF_AdditionalCost_Cloud.interface'

! Subroutine arguments:
REAL, INTENT(IN)    :: BmatrixInverse(:,:)  ! inverse of B-matrix
TYPE(R_Matrix_Type), INTENT(INOUT) :: R_SubMatrix ! chan-selected R matrix
REAL, INTENT(IN)    :: Delta_Profile(:)     ! Profile increments 
                                            ! (on minimisation levels)
REAL, INTENT(IN)    :: DeltaObs(:)           ! Obs-Ret BT difference
INTEGER, INTENT(IN) :: Work_Dimension
LOGICAL, INTENT(IN) :: Eqn_101       ! Rogers's Eqn. 101 is being used.
REAL, INTENT(IN)    :: Guess_Profile(:) ! Guess prof., used in Additional Cost
REAL, INTENT(OUT)   :: JCost         ! Cost function value
INTEGER, INTENT(OUT):: Status        ! Status of calculation (0=OK)
REAL, OPTIONAL, INTENT(IN)  :: H_Matrix_T(:,:)   ! Transpose of Jacobian
REAL, OPTIONAL, INTENT(OUT) :: JCost_Gradient(:) ! Gradient of cost function


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Calculate_Cost_Function"

! Local variables:

CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

INTEGER :: I

REAL, ALLOCATABLE :: AI(:,:)
REAL :: Jb, Jo, J_Additional ! Background, Observation and Additional 
                             ! components of the cost function.
REAL :: J_Additional_Gradient(Num_ProfElementsUsed) ! Gradient of J_Additional
REAL :: Scratch_vectorB (Num_ProfElementsUsed)
REAL :: Scratch_vectorO (Work_Dimension)

!-----------------------------------------------------------------------------

Status = 0
 
!-----------------------------------------------------------------------------
! 1. Calculate the background part of the cost function
!-----------------------------------------------------------------------------

Scratch_VectorB(:) = MATMUL(BmatrixInverse(:,:), Delta_Profile(:) )
Jb = DOT_PRODUCT(Delta_Profile(:) , Scratch_VectorB(:) )

!-----------------------------------------------------------------------------
! 2. Calculate the observational part of the cost function.
!    There are three ways of representing the observational error
!    covariance matrix which are treated separately.  In each case the
!    matrix is inverted if required. 
!-----------------------------------------------------------------------------

Jo = 0.

SELECT CASE (R_SubMatrix % RType)

!-----------------------------------------------------------------------------
! 2.1. Full matrix representation.
!-----------------------------------------------------------------------------

  CASE (R_Full_Matrix)
    IF (Eqn_101 .AND. .NOT. R_SubMatrix % Inverse) THEN
      ! This is the "poor man's cost funtion" calculation
      DO I = 1, R_SubMatrix % Num_Chans
        Jo = Jo + DeltaObs(I) * DeltaObs(I) / R_SubMatrix % Matrix(I,I)
      END DO
    ELSE    
      ! This is the full calculation 
      Scratch_VectorO(:) = DeltaObs(:)
      IF (.NOT. R_SubMatrix % Inverse) THEN
        ALLOCATE(AI(R_SubMatrix % Num_Chans,R_SubMatrix % Num_Chans))
        AI = R_SubMatrix % Matrix(:,:)
        CALL NWPSAF_SatMatInv(     &
          R_SubMatrix % Num_Chans, &
          1,                       &
          AI(:,:),                 &
          Matrix = Scratch_VectorO(:) )
        DEALLOCATE(AI)
      ELSE 
        Scratch_VectorO(:) = MATMUL(R_SubMatrix % Matrix(:,:), DeltaObs(:) )
      END IF
      Jo = DOT_PRODUCT(DeltaObs(:), Scratch_VectorO(:))
    END IF
!-----------------------------------------------------------------------------
! 2.2. Band diagonal matrix representation.
!-----------------------------------------------------------------------------
  CASE (R_Band_Diagonal)
     !  Diagonal Case 
    IF (R_SubMatrix % Num_Elements == 0 .OR. Eqn_101) THEN
      IF (R_SubMatrix % Inverse) THEN
        Scratch_VectorO(:) = DeltaObs(:) * R_SubMatrix % Matrix(0,:)
      ELSE
        Scratch_VectorO(:) = DeltaObs(:) / R_SubMatrix % Matrix(0,:)
      END IF
      Jo = DOT_PRODUCT(DeltaObs(:),Scratch_VectorO(:))
    ELSE  ! Non-Diagonal Case
      IF (.NOT.(R_SubMatrix % Inverse)) THEN
        ALLOCATE(AI(0:Max_Elements, R_SubMatrix % Num_Chans))
        CALL NWPSAF_BandInverse(     &
          R_SubMatrix % Num_Chans,   & ! in
          R_SubMatrix % Num_Elements,& ! in
          Max_Elements,              & ! in
          R_SubMatrix % Matrix(:,:), & ! in
          AI(:,:)                    ) ! out
        R_SubMatrix % Num_Elements = Max_Elements
        DEALLOCATE(R_SubMatrix % Matrix)
        ALLOCATE(R_SubMatrix % Matrix(0:Max_Elements,R_SubMatrix % Num_Chans))
        R_SubMatrix % Matrix(:,:) = AI(:,:)
        DEALLOCATE(AI)
        R_SubMatrix % Inverse = .TRUE.
      END IF
      CALL NWPSAF_BandMultiply (    &
        R_SubMatrix % Matrix(:,:),  & ! in
        DeltaObs(:),                 & ! in
        R_SubMatrix % Num_Chans,    & ! in
        R_SubMatrix % Num_Elements, & ! in
        1,                          & ! in
        Scratch_VectorO(:) )          ! out
      Jo = DOT_PRODUCT(DeltaObs(:), Scratch_VectorO(:))
    END IF
!-----------------------------------------------------------------------------
! 2.2. Eigenvector representation.
!-----------------------------------------------------------------------------
  CASE (R_Eigenvectors)
    Scratch_VectorO(:) = MATMUL(R_SubMatrix%Matrix(:,:), DeltaObs(:))
    IF (R_SubMatrix % Inverse) THEN
      Jo = DOT_PRODUCT(Scratch_VectorO(:), &
           Scratch_VectorO(:) * R_SubMatrix % EigenValues(:))
      IF (PRESENT(JCost_Gradient)) THEN
        IF (PRESENT(H_Matrix_T)) THEN
          JCost_Gradient(:) = MATMUL(                                     &
            MATMUL(H_Matrix_T(:,:),TRANSPOSE(R_SubMatrix % Matrix(:,:))), &
            Scratch_VectorO(:) * R_SubMatrix % EigenValues(:)             )
        ELSE
          Errormessage(1) = &
            'H_Matrix_T must be passed in to calculate JCost_Gradient'
          CALL NWPSAF_Report( &
            RoutineName,            & ! in
            ErrorMessage,           & ! in
            ErrorStatus=StatusFatal ) ! in
        END IF
      END IF
    ELSE
      Jo = DOT_PRODUCT(Scratch_VectorO(:), &
                       Scratch_VectorO(:) / R_SubMatrix % EigenValues(:) )
      IF (PRESENT(JCost_Gradient)) THEN
        IF (PRESENT(H_Matrix_T)) THEN 
          JCost_Gradient(:) = MATMUL(                                     &
            MATMUL(H_Matrix_T(:,:),TRANSPOSE(R_SubMatrix % Matrix(:,:))), &
                   Scratch_VectorO(:) / R_SubMatrix % EigenValues(:)      )
        ELSE
          Errormessage(1) = &
            'H_Matrix_T must be passed in to calculate JCost_Gradient'
          CALL NWPSAF_Report( &
            RoutineName,            & ! in
            ErrorMessage,           & ! in
            ErrorStatus=StatusFatal ) ! in
        END IF
      END IF
    END IF
  CASE DEFAULT
    Errormessage(1) = 'Unknown Format for Observation Errors'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
END SELECT

!-----------------------------------------------------------------------------
! 3. Add on any additional cost function terms.
!-----------------------------------------------------------------------------
IF (Additional_Cost_Function == No_Additional_Cost_Function) THEN
  J_Additional = 0.
  J_Additional_Gradient(:) = 0.
ELSE IF (Additional_Cost_Function == CloudBoundaries) THEN
  CALL NWPSAF_AdditionalCost_Cloud( &
    Guess_Profile,             &  ! in
    J_Additional,              &  ! out
    J_Additional_Gradient)        ! out
ELSE
  Errormessage(1) = 'No Additional Cost Functions are defined'
  CALL NWPSAF_Report( &
    RoutineName,              & ! in
    ErrorMessage,             & ! in
    ErrorStatus=StatusWarning ) ! in
   J_Additional = 0.
END IF

!-----------------------------------------------------------------------------
! 4. Calculate the complete cost function.
!-----------------------------------------------------------------------------

Jcost = (Jo + Jb)/2. + J_Additional

IF (PRESENT(JCost_Gradient)) THEN
  IF (R_SubMatrix % RType == R_Eigenvectors) THEN
    JCost_Gradient(:) =  Scratch_VectorB(:) - JCost_Gradient(:)
  ELSE
    IF (PRESENT(H_Matrix_T)) THEN 
      JCost_Gradient(:) = &
        Scratch_VectorB(:) - MATMUL(H_Matrix_T(:,:),Scratch_VectorO(:))
    ELSE
      Errormessage(1) = &
        'H_Matrix_T must be passed in to calculate JCost_Gradient'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in    
      END IF
   ENDIF
   JCost_Gradient(:) = JCost_Gradient(:) + J_Additional_Gradient(:)
ENDIF

End Subroutine NWPSAF_Calculate_Cost_Function
