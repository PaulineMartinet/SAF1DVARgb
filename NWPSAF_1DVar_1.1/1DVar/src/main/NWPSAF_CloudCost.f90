Subroutine NWPSAF_CloudCost( &
     MeasurementObs, &   ! in
     Bmatrix,    &       ! in
     R_Matrix,   &       ! in
     UsedChans,  &       ! in
     RT_Params,  &       ! inout
     Cost,       &       ! out
     RTErrorCode      )  ! out

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

! The difference between the observed and background radiances over the
! given channels are used to calculate a cost.  
!
! Cloud Cost, Jc, is calculated from : 
!
!   Jc   =   ( BT differences )   U   (BT differences)transpose
!
! with BT differences a vector, size of the number of channels chosen, of 
!                     ( observed BT - background BT)
!
!      U is the inverse of ( H.B.Ht + E + F )   
!         H is the partial derivatives of background BT's w.r.t the background
!         Ht is the transpose of H
!         B is the error covariance matrix
!         E the instrumental errors
!         F the forward model errors
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY : &
    StatusFatal, &
    CalcRadiance

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type, &
     Num_ProfElementsUsed, &
     FastmodelMode_Gradient, &
     BackGrProf

USE NWPSAFMod_CovarianceMatrices, ONLY : &
     R_Matrix_Type, &
     R_Full_Matrix, &
     R_Band_Diagonal, & 
     R_Eigenvectors

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_Cholesky.interface'
INCLUDE 'NWPSAF_Fastmodel_Interface.interface'
INCLUDE 'NWPSAF_RMatrix_ChanSelect.interface'

! Subroutine arguments:
REAL, INTENT(IN)     :: MeasurementObs(:) ! Observed BTs
REAL, INTENT(IN)     :: BMatrix(:,:)  ! Background Error Covariance Matrix
TYPE(R_Matrix_Type)        , INTENT(IN) :: R_Matrix  ! R-Matrix
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans ! channels to use
TYPE(RTParams_Type)     , INTENT(INOUT) :: RT_Params ! RT model input data
REAL, INTENT(OUT)    :: Cost                   ! the calculated cloud cost 
INTEGER, INTENT(OUT) :: RTErrorCode            ! Error code from RT model

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_CloudCost"

! Local variables:
CHARACTER (LEN=80):: ErrorMessage(1)
INTEGER :: element
INTEGER :: Fastmodel_Mode       !   Mode in which fastmodel is called
INTEGER :: I
REAL :: DeltaObs (UsedChans % NumChans)
REAL :: HB (UsedChans % NumChans,Num_ProfElementsUsed)
REAL :: U (UsedChans % NumChans, UsedChans % NumChans)
REAL :: Scratch_Matrix (UsedChans % NumChans)

REAL, ALLOCATABLE :: HTX (:,:)

TYPE(R_Matrix_Type) :: R_SubMatrix  ! Chan selected R Matrix

INTEGER :: WhichProf    

!- End of header ------------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Get the H matrix.
!-----------------------------------------------------------------------

WhichProf = BackGrProf
Fastmodel_Mode = FastmodelMode_Gradient
ALLOCATE(RT_Params % H_matrix(UsedChans % NumChans,Num_ProfElementsUsed))
ALLOCATE(RT_Params % H_matrix_T(Num_ProfElementsUsed,UsedChans % NumChans))
CALL NWPSAF_Fastmodel_interface(     &
  Fastmodel_Mode,          & ! in
  RT_Params,               & ! in
  WhichProf,               & ! in
  UsedChans,               & ! in
  RTerrorCode)               ! out

!-----------------------------------------------------------------------
! 2. Calculate H.B.H'
!-----------------------------------------------------------------------

HB = MATMUL(RT_Params % H_matrix,Bmatrix)
U = MATMUL(HB,RT_Params % H_matrix_T)

!-----------------------------------------------------------------------
! 3. Add the R matrix into the U matrix.
!-----------------------------------------------------------------------

CALL NWPSAF_RMatrix_ChanSelect( &
  UsedChans, & ! in
  R_Matrix,  & ! in
  R_SubMatrix) ! out

SELECT CASE (R_SubMatrix % RType)
   
CASE (R_Full_Matrix)
  U = U + R_SubMatrix % Matrix(:,:)
   
CASE (R_Band_Diagonal)
  DO Element = 0,R_SubMatrix % Num_Elements
    DO I=1, UsedChans % NumChans - Element
      U(I,I+Element) = U(I,I+Element) + R_SubMatrix % Matrix(Element,I)
      U(I+Element,I) = U(I,I+Element) 
    END DO
  END DO
  
CASE (R_Eigenvectors)
  ALLOCATE(HTX(Num_ProfElementsUsed, R_SubMatrix % Num_Elements))
  DO I=1,R_SubMatrix % Num_Elements
    HTX(:,I)=R_SubMatrix % EigenValues(I) * &
         R_SubMatrix % Matrix(:,I)
  END DO
  U = U + MATMUL(TRANSPOSE(R_SubMatrix % Matrix(:,:)),HTX)
   
CASE DEFAULT
  Errormessage(1) = 'Unknown Format for Observation Errors'
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=StatusFatal )
   
END SELECT

!-----------------------------------------------------------------------
! 4. Calculate the difference between the observed radiances and
!    those expected from the background.
!-----------------------------------------------------------------------

If ( CalcRadiance ) Then
  DeltaObs(:) = &
    MeasurementObs(UsedChans % Channels(1 : UsedChans % NumChans)) - &
    RT_Params % TotalRadiances(1 : UsedChans % NumChans)
Else
  DeltaObs(:) = &
    MeasurementObs(UsedChans % Channels(1 : UsedChans % NumChans)) - &
    RT_Params % TotalBTs(1 : UsedChans % NumChans)
End If
!-----------------------------------------------------------------------
! 5. Calculate U^-1.DeltaObs using Cholesky decomposition and 
!    then premultiply by DeltaObs^T to get the cloud cost.
!-----------------------------------------------------------------------

CALL NWPSAF_CHOLESKY( &
  U,                    & !in
  DeltaObs,         & !in
  UsedChans % NumChans, & !in
  Scratch_Matrix)              !out

Cost = 0.5*DOT_PRODUCT(DeltaObs,Scratch_matrix)

!-----------------------------------------------------------------------
! 6. Clean up.
!-----------------------------------------------------------------------

IF (R_SubMatrix % RType == R_Eigenvectors) DEALLOCATE(HTX)
DEALLOCATE(R_SubMatrix % Matrix)
DEALLOCATE(RT_Params % H_matrix)
DEALLOCATE(RT_Params % H_matrix_T)
  
End Subroutine NWPSAF_CloudCost
