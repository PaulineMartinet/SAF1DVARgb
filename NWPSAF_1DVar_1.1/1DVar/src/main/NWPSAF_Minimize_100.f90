!+ Find a solution to the satellite sounding inverse problem

Subroutine NWPSAF_Minimize_100( &
     BmatrixInverse,   & ! in
     DeltaBT,          & ! in  
     Num_Chans,        & ! in
     Work_Dimension,   & ! in
     H_Matrix,         & ! in
     H_Matrix_T,       & ! in
     R_SubMatrix,      & ! inout
     Delta_Profile     ) ! inout

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
!   Eqn. 100, extended to allow for additional cost function terms.
!
!   x_(n+1) = xb + U^-1.V
!   where U=(B^-1 + H^T R^-1 H + J2)
!         V=H^T R^-1 [(y-y(x_n))+H(x_n-xb)] - J1
!   and   J_extra=J0+J1.(x-xb)+(x-xb)^T.J2.(x-xb) is the additional cost 
!         function
!
!   When J_extra is zero this is simply Rogers (1976), Eqn. 100.
!
!   U^-1.V is solved using Cholesky decomposition.
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
!                  minimisation routine.  Andrew Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type,             &     
    R_Full_Matrix,             &
    R_Band_Diagonal,           & 
    R_Eigenvectors,            &
    Max_Elements,              &
    Analysis_Error_Covariance, &
    Prop_Measurement_Noise

USE NWPSAFMod_Params, ONLY :    &
    StatusFatal,              &
    Additional_Cost_Function, &
    No_Additional_Cost_Function

USE NWPSAFMod_RTModel, ONLY :    &
    Num_ProfElementsUsed 

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_SatMatInv.interface'
INCLUDE 'NWPSAF_BandInverse.interface'
INCLUDE 'NWPSAF_BandMultiply.interface'
INCLUDE 'NWPSAF_Cholesky.interface'

! Subroutine arguments:
INTEGER, INTENT(IN) :: Num_Chans           ! No. of channels to be used 
INTEGER, INTENT(IN) :: Work_Dimension   
REAL, INTENT(IN)    :: BmatrixInverse(:,:) ! inverse of B-matrix
REAL, INTENT(IN)    :: DeltaBT(:)          ! y-y(x)
REAL, INTENT(INOUT) :: Delta_Profile(:)    ! profile increments 
REAL, INTENT(IN)    :: H_Matrix(:,:)       ! Jacobian
REAL, INTENT(IN)    :: H_Matrix_T(:,:)     ! (Jacobian)^T
TYPE(R_Matrix_Type), INTENT(INOUT) :: R_SubMatrix  ! R Matrix for selected 
                                                   ! channels

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Minimize_100"

! Local variables:
CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

INTEGER :: I
REAL, ALLOCATABLE :: AI(:,:)    ! Temporary matrix for inversion

!xxx REAL, ALLOCATABLE :: J1(:)   ! 1st derivative of additional cost function 
!xxx                              ! wrt profile
!xxx REAL, ALLOCATABLE :: J2(:,:) ! 2nd derivative of additional cost function 
!xxx                              ! wrt profile

REAL :: HTR(Num_ProfElementsUsed, Work_Dimension)      ! Scratch vector
REAL :: U (Num_ProfElementsUsed, Num_ProfElementsUsed) ! U = H.B.H^T + R
REAL :: V (Num_ProfElementsUsed)                 ! V = (y-y(x_n))-H^T(xb-x_n)

!-----------------------------------------------------------------------------

!---------------------------------------------------------------------------
! 0.1. For now do not allow additional cost function terms.  This section 
!      should be replaced by a function call when the appropriate code is 
!      written.
!---------------------------------------------------------------------------

IF (Additional_Cost_Function > No_Additional_Cost_Function) &
    Additional_Cost_Function = No_Additional_Cost_Function


!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.  In all cases, the matrix is tested to determine whether it is 
!    stored as an inverse and inverted if not.
!---------------------------------------------------------------------------

SELECT CASE (R_SubMatrix % RType)

  !-------------------------------------------------------------------------
  ! 1.1.1 Full matrix case.  The inverse, if required, is determined using 
  !     Cholesky decomposition.
  !-------------------------------------------------------------------------

  CASE (R_Full_Matrix)
    IF (.NOT. R_SubMatrix % Inverse) THEN
      CALL NWPSAF_SatMatInv(      &
        Num_Chans,                &
        Num_Chans,                &
        R_SubMatrix % Matrix(:,:) )
      R_SubMatrix % Inverse = .TRUE.
    END IF
     HTR = MATMUL(H_matrix_T, R_SubMatrix % Matrix(:,:))
     U = MATMUL(HTR, H_matrix)
     V = MATMUL(HTR,DeltaBT)     
     V = V + MATMUL(U,Delta_Profile)

     !-----------------------------------------------------------------------
     ! 1.1.2 Band-diagonal matrix case.  The inverse, if required, is 
     !       determined using Cholesky decomposition.
     !-----------------------------------------------------------------------
  CASE (R_Band_Diagonal)
    ! Diagonal Case
    IF (R_SubMatrix % Num_Elements == 0) THEN
      IF (R_SubMatrix % Inverse) THEN
        DO I=1,R_SubMatrix % Num_Chans
          HTR(:,I) = H_Matrix_T(:,I) * R_SubMatrix % Matrix(0,I)
        END DO
      ELSE
        DO I=1,R_SubMatrix % Num_Chans
          HTR(:,I) = H_Matrix_T(:,I) / R_SubMatrix % Matrix(0,I)
        END DO
      END IF
    ELSE
    ! Non-Diagonal Case
      IF (.NOT. R_SubMatrix % Inverse) THEN
        ALLOCATE(AI(0:Max_Elements, Num_Chans))
        CALL NWPSAF_BandInverse (         &
          R_SubMatrix % Num_Chans,   & ! in
          R_SubMatrix % Num_Elements,& ! in
          Max_Elements,              & ! in
          R_SubMatrix % Matrix(:,:), & ! in
          AI                         ) ! out
        R_SubMatrix % Num_Elements = Max_Elements
        DEALLOCATE(R_SubMatrix % Matrix)
        ALLOCATE(R_SubMatrix % Matrix(0:Max_Elements, Num_Chans))
        R_SubMatrix % Matrix(:,:) = AI
        DEALLOCATE(AI)
        R_SubMatrix % Inverse = .TRUE.
      END IF
      CALL NWPSAF_BandMultiply( &
        R_SubMatrix % Matrix(:,:), &
        H_Matrix_T,                &
        R_SubMatrix % Num_Chans,   &
        R_SubMatrix % Num_Elements,&
        Num_ProfElementsUsed,      &
        HTR)     
    END IF
    U = MATMUL(HTR, H_matrix)
    V = MATMUL(HTR,DeltaBT)
    V = V + MATMUL(U,Delta_Profile)

    !------------------------------------------------------------------------
    ! 1.1.3 R-matrix represented by eigenvectors case.
    !------------------------------------------------------------------------
  CASE (R_Eigenvectors)
    ! The automatic array HTR should have a second dimension 
    ! R_Matrix % Num_Elements.  Check that this is the case.
    IF (Work_Dimension /= R_SubMatrix % Num_Elements) THEN
      Errormessage(1) = 'Incorrect value for HTR_Dimension'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
    HTR = MATMUL(H_matrix_T, R_SubMatrix % Matrix(:,:))
    DO I=1,R_SubMatrix % Num_Elements
      HTR(:,I)=HTR(:,I)/R_SubMatrix % EigenValues(I)
    END DO
    U = MATMUL(HTR, TRANSPOSE(HTR))
    V = MATMUL(TRANSPOSE(R_SubMatrix % Matrix(:,:)),DeltaBT)
    V = MATMUL(HTR, V)
    V = V + MATMUL(U,Delta_Profile)

    !------------------------------------------------------------------------
    ! 1.1.4 No other R-matrix representations allowed.
    !------------------------------------------------------------------------

  CASE DEFAULT
    Errormessage(1) = 'Unknown Format for Observation Errors'
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 

END SELECT

! Add on additional terms to U and V due to additional cost function
! terms.

!---------------------------------------------------------------------------
! 1.2 Add on inverse of B-matrix to U.
!---------------------------------------------------------------------------

Prop_Measurement_Noise = U  ! This needs to be pre and post multiplied 
                            ! by A^-1 still (done in calling routine)

U = U + BmatrixInverse

Analysis_Error_Covariance = U  ! This is still the inverse of A

!---------------------------------------------------------------------------
! 1.3 Add on additional cost function terms to U and V if desired.
!---------------------------------------------------------------------------

!xxx IF (Additional_Cost_Function > No_Additional_Cost_Function) THEN
!xxx    CALL NWPSAF_Addition_Cost_Function( &
!xxx         J1, &   ! out
!xxx         J2)     ! out
!xxx    U = U + J2
!xxx    V = V - J1
!xxx END IF

!---------------------------------------------------------------------------
! 1.4 Calculate new profile increment.
!---------------------------------------------------------------------------

CALL NWPSAF_CHOLESKY(   &
  U,                    & !in              
  V,                    & !in
  Num_ProfElementsUsed, & !in
  Delta_Profile         ) !out

End Subroutine NWPSAF_Minimize_100
