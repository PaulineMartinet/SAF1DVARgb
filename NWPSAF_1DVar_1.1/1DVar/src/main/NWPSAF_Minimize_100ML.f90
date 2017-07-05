!+ Find a solution to the satellite sounding inverse problem

Subroutine NWPSAF_Minimize_100ML( &
     Observation,       & ! in 
     RT_Params,        & ! in
     BmatrixInverse,   & ! in
     Back,             & ! in
     UsedChans,        & ! in
     Work_Dimension,   & ! in
     R_SubMatrix,      & ! inout
     DeltaObs,         & ! inout  
     Delta_Profile,    & ! inout
     Gamma,            & ! inout
     JOld,             & ! inout
     Status,           & ! out
     Out_of_Range     )  ! out

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
!   Eqn. 100, extended to allow for additional cost function terms and
!   Marquardt-Levenberg descent.
!
!   x_(n+1) = xb + U^-1.V
!   where U=(B^-1 + H^T R^-1 H + J2 + gamma I)
!         V=H^T R^-1 [(y-y(x_n))+H(x_n-xb)] + gamma (x_n-xb) - J1 + J2(x_n-x_0)
!   and   J1 and J2 are the first and second derivatives of the additional
!         cost function respectively.
!         gamma = is the multiplier used in the Marquart-Levenberg routine
!                                    (set to zero for simple inverse Hessian).
!
!   When gamma and J_extra are zero this is simply Rogers (1976), Eqn. 100.
!
!   When gamma -> infinity, x_(n+1) -> x_n, i.e, the step size is zero.
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
! 1.2     25/02/00 Marquardt-Levenberg version.  
!                  Brought more code in from NWPSAF_Minimisation so that 
!                  the inner loop of the Marquardt-Levenberg algorithm 
!                  entirely in this subroutine.
!                                               Andrew Collard
!                                               Met Office
! 2.3     23/05/02 Add code to prevent the main loop running away.
!                  Modified call to NWPSAF_Calculate_Cost_Function.
!                  Extra comment line.          A. Collard.
! 3.0.1   15/07/03 Add M. Szyndel's test for maximum gamma.  A. Collard.
!                  Fix errant reseting of Gamma in descent loop.
!                                                            A. Collard.
! 3.0.2   04/02/04 Included additional cost function. 
!                  Pass RT_Params%Guess to NWPSAF_Calculate_Cost_Function
!                  to be used in additional cost function calc. ADC.
! 3.0.4   04/03/04 Rename NWPSAF_MaxIterations, MaxIterations.  Pass presures
!                  to NWPSAF_CheckIteration.                      A. Collard.
! 3.0.5   08/04/04 Remove tolerance from JCost test and allow minimisation
!                  to continue if the maximum number of ML iterations is
!                  performed.
! 3.0.6   21/06/04 Remove pressure argument to NWPSAF_CheckIteration.
!                                                          A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
! 32        14/08/13 Added call to NWPSAF_LWP_to_Layers for update of clw profile
!                                                        A. Andersson (CM SAF)
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type, &
    R_Full_Matrix, &
    R_Band_Diagonal, & 
    R_Eigenvectors, &
    Max_Elements,   &
    Analysis_Error_Covariance, &
    Prop_Measurement_Noise

USE NWPSAFMod_Params, ONLY :    &
     StatusFatal, &
     GeneralMode, &
     VerboseMode, &
     Additional_Cost_Function, &
     CloudBoundaries, &
     Retrieved_Elements, &
     Max_ML_Iterations, &
     Gamma_Factor, &
     GammaMax, &
     MwClwRetrieval, &
     Retrieve_LWP, &
     UsePCs, &
     CalcRadiance

USE NWPSAFMod_RTmodel, ONLY : &
    RTParams_Type,        &
    Num_ProfElementsUsed, &
    FastmodelMode_Forward, &
    GuessProf, &
    Num_RTlevels, &
    Prof_LWP, &
    Prof_FirstCLW, &
    Prof_lastCLW

Use NWPSAFMod_LiquidWater, Only : &
    LWP_to_Layers

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_AdditionalCost_Cloud.interface'
INCLUDE 'NWPSAF_Calculate_Cost_Function.interface'
INCLUDE 'NWPSAF_CheckIteration.interface'
INCLUDE 'NWPSAF_Fastmodel_Interface.interface'
INCLUDE 'NWPSAF_SatMatInv.interface'
INCLUDE 'NWPSAF_BandInverse.interface'
INCLUDE 'NWPSAF_BandMultiply.interface'
INCLUDE 'NWPSAF_Cholesky.interface'
INCLUDE 'NWPSAF_CloudStructure.interface'

! Subroutine arguments:
REAL, INTENT(IN)     :: Observation(:)                   ! Observation BT's
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params         ! RT Model Data
REAL, INTENT(IN)     :: BmatrixInverse(:,:)             ! inverse of B-matrix
REAL, INTENT(IN)     :: Back(:)                         ! Background profile
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans    ! Chans to be used
INTEGER, INTENT(IN) :: Work_Dimension              ! Dimension of work array
TYPE(R_Matrix_Type), INTENT(INOUT) :: R_SubMatrix  ! R Matrix for selected 
                                                   ! channels
REAL, INTENT(INOUT)  :: DeltaObs(:)                ! y-y(x)
REAL, INTENT(INOUT)  :: Delta_Profile(:)           ! profile increments 
REAL, INTENT(INOUT)  :: Gamma          ! Used in Marquart-Levenberg minimisation
REAL, INTENT(INOUT)  :: JOld           ! Cost function value
INTEGER, INTENT(OUT) :: Status
LOGICAL, INTENT(OUT) :: Out_of_Range   ! Gross limts exceeded in minimisation

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Minimize_100ML"

! Local variables:
CHARACTER(LEN=80)  :: ErrorMessage(1)  ! Message for NWPSAF_Report

INTEGER :: I
INTEGER :: Fastmodel_Mode        ! Mode in which fastmodel is called
INTEGER :: Iter_ML               ! Number of M-L iterations

REAL, ALLOCATABLE :: AI(:,:)     ! Temporary matrix for inversion
REAL              :: JCost       ! Cost function value
REAL              :: J_Add       ! Additional cost function value
REAL :: J1(Num_ProfElementsUsed)       ! Gradient of J_Add
REAL :: J2(Num_ProfElementsUsed,Num_ProfElementsUsed) ! 2nd Deriv of J_Add

REAL :: New_Delta_Profile (Num_ProfElementsUsed)   
REAL :: V     (Num_ProfElementsUsed)             ! V = (y-y(x_n))-H^T(xb-x_n)
REAL :: VSave (Num_ProfElementsUsed)             

REAL :: HTR   (Num_ProfElementsUsed, Work_Dimension)       ! Scratch matrix
REAL :: U     (Num_ProfElementsUsed, Num_ProfElementsUsed) ! U = H.B.H^T + R
REAL :: USave (Num_ProfElementsUsed, Num_ProfElementsUsed) 

LOGICAL :: Profile_Variables_Reset

INTEGER :: WhichProf
REAL :: cloud_structure(Num_RTlevels)
!-----------------------------------------------------------------------------
!---------------------------------------------------------------------------
! 1. Calculate the U and V vectors for the three forms of R matrix allowed
!    for.  In all cases, the matrix is tested to determine whether it is 
!    stored as an inverse and inverted if not.
!---------------------------------------------------------------------------

SELECT CASE (R_SubMatrix % RType)

   !------------------------------------------------------------------------
   ! 1.1.1 Full matrix case. 
   !------------------------------------------------------------------------

  CASE (R_Full_Matrix)
    IF (.NOT. R_SubMatrix % Inverse) THEN
      CALL NWPSAF_SatMatInv(      &
        UsedChans % NumChans,     &
        UsedChans % NumChans,     &
        R_SubMatrix % Matrix(:,:) )
      R_SubMatrix % Inverse = .TRUE.
    END IF
    HTR = MATMUL(RT_Params % H_matrix_T(:,:), R_SubMatrix % Matrix(:,:))
    U = MATMUL(HTR, RT_Params % H_matrix(:,:) )
    V = MATMUL(HTR,DeltaObs)
    V = V + MATMUL(U,Delta_Profile)

  !------------------------------------------------------------------------
  ! 1.1.2 Band-diagonal matrix case.  
  !------------------------------------------------------------------------
  CASE (R_Band_Diagonal)
    IF (R_SubMatrix % Num_Elements == 0) THEN
      ! Diagonal Case
      IF (R_SubMatrix % Inverse) THEN
        DO I=1,R_SubMatrix % Num_Chans
          HTR(:,I) = RT_Params % H_Matrix_T(:,I) * R_SubMatrix % Matrix(0,I)
        END DO
      ELSE
        DO I=1,R_SubMatrix % Num_Chans
          HTR(:,I) = RT_Params % H_Matrix_T(:,I) / R_SubMatrix % Matrix(0,I)
        END DO
      END IF

    ELSE
      ! Non-Diagonal Case
      IF (.NOT. R_SubMatrix % Inverse) THEN
        ALLOCATE(AI(0:Max_Elements, UsedChans % NumChans))
        CALL NWPSAF_BandInverse (         &
          R_SubMatrix % Num_Chans,   & ! in
          R_SubMatrix % Num_Elements,& ! in
          Max_Elements,              & ! in
          R_SubMatrix % Matrix(:,:), & ! in
          AI                        ) ! out
        R_SubMatrix % Num_Elements = Max_Elements
        DEALLOCATE(R_SubMatrix % Matrix)
        ALLOCATE(R_SubMatrix % Matrix(0:Max_Elements, UsedChans % NumChans))
        R_SubMatrix % Matrix(:,:) = AI
        DEALLOCATE(AI)
        R_SubMatrix % Inverse = .TRUE.
      END IF
      CALL NWPSAF_BandMultiply( &
        R_SubMatrix % Matrix(:,:),   &
        RT_Params % H_Matrix_T(:,:), &
        R_SubMatrix % Num_Chans,     &
        R_SubMatrix % Num_Elements,  &
        Num_ProfElementsUsed,        &
        HTR)     
    END IF
    U = MATMUL(HTR, RT_Params % H_matrix(:,:) )
    V = MATMUL(HTR,DeltaObs)
    V = V + MATMUL(U,Delta_Profile)

  !------------------------------------------------------------------------
  ! 1.1.3 R-matrix represented by eigenvectors case.
  !------------------------------------------------------------------------
  CASE (R_Eigenvectors)
    ! The automatic array HTR should have a second dimension 
    ! R_Matrix % Num_Elements.  Check that this is the case.
    IF (Work_Dimension /= R_SubMatrix % Num_Elements) THEN
      Errormessage(1) = 'Incorrect value for Work_Dimension'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
   END IF
   HTR = MATMUL(RT_Params % H_matrix_T(:,:), R_SubMatrix % Matrix(:,:))
   DO I=1,R_SubMatrix % Num_Elements
     HTR(:,I)=HTR(:,I)/R_SubMatrix % EigenValues(I)
   END DO
   U = MATMUL(HTR, TRANSPOSE(HTR))
   V = MATMUL(TRANSPOSE(R_SubMatrix % Matrix(:,:)),DeltaObs)
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
! terms and the use of the Marquardt-Levenberg method.

!---------------------------------------------------------------------------
! 1.2 Add on inverse of B-matrix to U.
!---------------------------------------------------------------------------

Prop_Measurement_Noise = U  ! This needs to be pre and post multiplied 
                            ! by A^-1 still.
U = U + BmatrixInverse

Analysis_Error_Covariance = U  ! This is still the inverse of A

!---------------------------------------------------------------------------
! 2. Start Marquardt-Levenberg loop
!---------------------------------------------------------------------------

JCost = 1.e10
Gamma = Gamma/Gamma_Factor
USave = U
VSave = V
Iter_ML = 0

DescentLoop : DO WHILE ( JCost > JOld .AND. &
                         Iter_ML <= Max_ML_Iterations .AND. Gamma <= GammaMax )
  Iter_ML = Iter_ML + 1
  IF (GeneralMode  >= VerboseMode) WRITE(*,*) 'Iter_ML=',Iter_ML

  !------------------------------------------------------------------------
  ! 2.1 Add on Marquardt-Levenberg terms to U and V.
  !------------------------------------------------------------------------
  Gamma = Gamma*Gamma_Factor

  DO I=1,Num_ProfElementsUsed
    U(I,I) = USave(I,I) + Gamma
  END DO
  V = VSave + Gamma*Delta_Profile  

  !---------------------------------------------------------------------------
  ! 1.3 Add on additional cost function terms to U and V if required.
  !---------------------------------------------------------------------------
   
  IF (Additional_Cost_Function == CloudBoundaries) THEN
    CALL NWPSAF_AdditionalCost_Cloud( &
      RT_Params % RTGuess,       &  ! in 
      J_Add,                     &  ! out
      J1,                        &  ! out
      J2                         )  ! optional out 
    U = U + J2
    V = V - J1 + MATMUL(J2,Delta_Profile)
  END IF

  !------------------------------------------------------------------------
  ! 2.2 Calculate new profile increment.
  !------------------------------------------------------------------------

  CALL NWPSAF_CHOLESKY( &
    U,                    & !in
    V,                    & !in
    Num_ProfElementsUsed, & !in
    New_Delta_Profile)      !out

  !------------------------------------------------------------------------
  ! 2.3. Calculate the radiances for the new profile.
  !------------------------------------------------------------------------
  RT_Params % RTGuess(Retrieved_Elements) = Back + New_Delta_Profile
  WhichProf = GuessProf
  Fastmodel_Mode = FastmodelMode_Forward
   
  CALL NWPSAF_CheckIteration( &
    RT_Params % RTGuess,      & ! inout
    New_Delta_Profile,        & ! inout
    Profile_Variables_Reset,  & ! out
    Out_of_Range )              ! out

   
  IF (Profile_Variables_Reset) THEN
    IF ( GeneralMode >= VerboseMode ) WRITE(*,*) 'Profile Variables Reset'
  END IF
   
  IF(MwClwRetrieval .AND. Retrieve_LWP) THEN
    !Spread LWP to layers to take account of modified LWP.
    CALL NWPSAF_CloudStructure(RT_params,cloud_structure)
    RT_Params%RTguess(Prof_FirstCLW:Prof_LastCLW) = &
      LWP_to_Layers( RT_Params % RTguess(Prof_LWP), cloud_structure )
  END IF

  CALL NWPSAF_Fastmodel_Interface(     &
    Fastmodel_Mode,          & ! in
    RT_Params,               & ! in
    WhichProf,               & ! in
    UsedChans,               & ! in
    Status                   ) ! in

   !------------------------------------------------------------------------
   ! 2.3. Calculate the new cost function.
   !------------------------------------------------------------------------
 
  If (UsePCs)  Then
    DeltaObs = Observation(UsedChans % Channels(1 : UsedChans % NumChans)) - &
              RT_Params % PCScores(1 : UsedChans % NumChans)
  Else If (CalcRadiance) Then
    DeltaObs = Observation(UsedChans % Channels(1 : UsedChans % NumChans)) - &
              RT_Params % TotalRadiances(1 : UsedChans % NumChans)
  Else
    DeltaObs = Observation(UsedChans % Channels(1 : UsedChans % NumChans)) - &
              RT_Params % TotalBTs(1 : UsedChans % NumChans)
  End If


  CALL NWPSAF_Calculate_Cost_Function( &
    BMatrixInverse,              & ! in 
    R_SubMatrix,                 & ! inout
    New_Delta_Profile,           & ! in 
    DeltaObs,                    & ! in           
    Work_Dimension,              & ! in 
    .FALSE.,                     & ! in 
    RT_Params % RTGuess,         & ! in
    JCost,                       & ! out
    Status)                        ! out

END DO DescentLoop

Delta_Profile = New_Delta_Profile
Gamma = Gamma / Gamma_Factor
JOld = JCost

End Subroutine NWPSAF_Minimize_100ML
