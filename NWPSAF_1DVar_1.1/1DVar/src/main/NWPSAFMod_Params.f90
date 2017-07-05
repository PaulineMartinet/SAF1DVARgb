!+ Parameters for global access in NWPSAF_1DVar

Module NWPSAFMod_Params

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
! Description: contains parameters for global access in NWPSAF_1DVar
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     25/05/99 Original version modified from ATOVSMod_Params.  
!                                         A.D. Collard  Met. Office
! 1.1     01/03/01 Added ATOVS Roger Saunders
! 3.0.1   15/07/03 Add cloud retrieval types (from M. Szyndel) A. Collard.
! 3.0.2   04/02/04 Add Ret_CTP and Ret_CloudCover and additional cost
!                  function.   A. Collard.
! 3.0.3   26/02/04 Add RTModelToUse                            A. Collard.
! 3.0.4   02/03/04 Changes to make instrument setup generic.   A. Collard.
! 3.0.5   29/03/04 Remove all references to orbit height.     
!                  Move some minimisation thresholds here.     A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date      Comment
! ------  --------  -------
! 19      30/03/10  Add MaxLevs. Fiona Hilton
! 25      20/02/12  Add emissivity related variables. P. Weston
! 28      22/02/12  Add MwClwRetrieval, Lqtotal and Ret_CLW. TR Sreerekha
! 25      29/03/12  Add Ozone_Present. P. Weston
! 32      07/12/12  Added WindspeedSD, LwpSD, Ret_UWind, Ret_VWind
!                   A. Andersson (CM SAF)
! 32      15/11/13  Add Output_Dir. P. Weston
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------


USE NWPSAFMod_Channellist, ONLY: &
     MaxChanUsed

USE NWPSAFMod_Constants, ONLY: &
     MissData_R, &
     MissData_I

 
!--------------------------------------
! 0. Reporting and Fortran_IO Options
!--------------------------------------

INTEGER, PARAMETER  :: StatusOK      = 0  ! Good ErrorStatus value.
INTEGER, PARAMETER  :: StatusWarning = -1 ! Non fatal error ErrorStatus value.
INTEGER, PARAMETER  :: StatusFatal   = 1  ! Fatal error ErrorStatus value.

INTEGER, PARAMETER :: LB = 10               ! Lower bound
INTEGER, PARAMETER :: UB = 64               ! Upper bound 
                                            ! (max value on HP is 63)
 
INTEGER, PRIVATE   :: i                     ! used in declaration of 
                                            ! UnitsInUse, below.
 
! Global Arrays:

LOGICAL, SAVE  :: UnitsInUse (LB:UB)  & ! If T then unit is in use. 
                       = (/ (.FALSE., i = LB, UB) /) 

!----------------------
!1. Control Options
!------------------
!----------------------
! 1.1. General Control
!----------------------

! This controls the amount of diagnostic output:
INTEGER, PARAMETER   :: OperationalMode = 0
INTEGER, PARAMETER   :: ProductionMode  = 10 ! for non-operational production
INTEGER, PARAMETER   :: DiagnosticMode  = 20 ! for normal program development
INTEGER, PARAMETER   :: DebugMode       = 30 ! for testing and debugging
INTEGER, PARAMETER   :: QuickNDirtyMode = 35 ! obsolete
INTEGER, PARAMETER   :: VerboseMode     = 40 ! for detailed program tracing
INTEGER , SAVE :: GeneralMode = ProductionMode ! mode of run

!---------------------------------
! 1.2. General Processing Options
!---------------------------------

INTEGER, PARAMETER :: &
  No_Minimisation     = 0, &
  Newtonian           = 1, &
  Marquardt_Levenberg = 2

INTEGER, PARAMETER :: &
  No_Additional_Cost_Function    = 0, &
  CloudBoundaries                = 1

LOGICAL, SAVE :: Perform1DVar  = .TRUE.
LOGICAL, SAVE :: DetectCloud   = .TRUE.      ! flag for cloud retrieval
INTEGER, SAVE :: MaxIterations =  10         ! No. iterations for IASI inversion
INTEGER, SAVE :: Max_ML_Iterations = 10      ! Max. no of inner loop iterations
                                             ! for Marquardt-Levenberg 
                                             ! minimisation      
REAL, SAVE    :: Gamma_Factor = 10.          ! Change in gamma between ML iters
REAL, SAVE    :: GammaMax = 1.0E25           ! Maximum value of Gamma for 
                                             ! Marquardt-Levenberg
REAL, SAVE    :: DeltaJ = 0.01               ! Max change in cost for 
                                             ! convergence
REAL, SAVE    :: SmallJCost_Gradient = 1.    ! Convergence criterion for the
                                             ! cost function gradient (ML only)
INTEGER, SAVE :: Minimisation_Method = Newtonian
                                           ! Newtonian or Marquardt-Levenberg
                                           ! (Initialised to missing)
INTEGER, SAVE :: Additional_Cost_Function = No_Additional_Cost_Function

LOGICAL, SAVE :: DoTextrapolation   = .TRUE. ! Extrapolate temperature above 
                                             ! model top
LOGICAL, SAVE :: Ozone_Present = .FALSE.   ! True if ozone appears in the 
                                           ! background profile
LOGICAL, SAVE :: CloudyRetrieval = .FALSE. ! True if cloud parameters are 
                                           ! to be retrieved
LOGICAL, SAVE :: MwClwRetrieval  = .FALSE. ! True if cloud parameters are 
                                           ! to be retrieved
LOGICAL, SAVE :: Retrieve_qtotal = .FALSE.
LOGICAL, SAVE :: Retrieve_LWP = .TRUE.

LOGICAL, SAVE :: Read_CLW_Background = .FALSE.

LOGICAL, SAVE :: Allow_Eqn_101 = .TRUE.   
LOGICAL, SAVE :: Force_Eqn_101 = .FALSE.   
LOGICAL, SAVE :: EnhancedDiagnostics = .TRUE. !Output Jacobians and Averaging Kernels
                                              !Along with QC information about 
                                              !Whether ob minimised OK.   


INTEGER, SAVE :: FirstOb  = 0  ! First Observation to be processed
INTEGER, SAVE :: LastOb   = 0  ! Last Observation to be processed

CHARACTER(LEN=200), SAVE :: Output_Dir='.' ! Output directory

! Set up for PC-RTTOV: AIRS and IASI only
LOGICAL, SAVE :: UsePCs        = .FALSE.      ! flag simulation of PCs
INTEGER, SAVE :: ipcreg        = 3            ! code for how many predictors
                                              ! (see RTTOV manual for details)
INTEGER, SAVE :: NPCScores     = 400          ! how many pcs

LOGICAL, SAVE :: UseRRs        = .FALSE.      !Setting this will set UsePCs to true

LOGICAL, SAVE :: CalcRadiance  = .FALSE.      !BTs by default (works also for 
                                              !reconstructed radiances/BTs)
INTEGER, SAVE :: Gas_Units     = 2            !ppmv moist air, 0=dry air
                                              !See RTTOV documentation
LOGICAL, SAVE :: Legacy_Setting = .FALSE.
!----------------------------------------
! 1.4. Profile Information
!----------------------------------------

INTEGER, PARAMETER :: MaxLevs=200 !Make sure this is big enough to fit whole
                                  !profile

INTEGER, SAVE :: Ret_FirstQ = 0  
INTEGER, SAVE :: Ret_LastQ = 0    
INTEGER, SAVE :: Ret_q2 = 0
INTEGER, SAVE :: Ret_CTP = 0
INTEGER, SAVE :: Ret_CloudCover = 0
INTEGER, SAVE :: Ret_CLW = 0
INTEGER, SAVE :: Ret_UWind = 0
INTEGER, SAVE :: Ret_VWind = 0

INTEGER, ALLOCATABLE, SAVE :: Retrieved_Elements(:) 

!-------------------------------
! 1.4.1. B-matrix information
!-------------------------------
INTEGER, ALLOCATABLE, SAVE :: B_ElementsUsed(:) 
REAL,SAVE :: WindspeedSD = 1.4 ! Standard deviation for wind speed retrieval 
                               ! (value for each component u,v)
REAL,SAVE :: LwpSD = 0.2       ! Standard deviation for LWP retrieval; 
                               ! default 0.2

!----------------------
! 1.4.2. Humidity Units
!----------------------
INTEGER, PARAMETER :: Humidity_PPMV    = 1
INTEGER, PARAMETER :: Humidity_MassMix = 2
INTEGER, PARAMETER :: Humidity_RH      = 3
INTEGER, SAVE :: Humidity_Units

!----------------------------------------
! 1.5. Auxiliary datasets
!----------------------------------------

INTEGER, PARAMETER :: NumFiles = 3
CHARACTER (LEN=200), SAVE :: Coeffs_Dir
CHARACTER (LEN=17), PARAMETER :: filenames(NumFiles) = &
  (/'Bmatrix          ', 'Rmatrix          ','ChannelChoice.dat'  /)
INTEGER, PARAMETER :: FileType_Bmatrix       = 1
INTEGER, PARAMETER :: FileType_Rmatrix       = 2
INTEGER, PARAMETER :: FileType_ChannelChoice = 3

!--------------------------------------
! 1.6. Cloud Detection and retrieval.
!-------------------------------------

REAL, SAVE :: CostThresh_Land             = MissData_R
REAL, SAVE :: CostThresh_Sea              = MissData_R
REAL, SAVE :: CostThresh_IRWindow_Land    = MissData_R
REAL, SAVE :: CostThresh_IRWindow_Sea     = MissData_R
REAL, SAVE :: CloudAbsThresh_IRWindow     = MissData_R
REAL, SAVE :: HighCloudAbsThresh_IRWindow = MissData_R

REAL, SAVE :: Cloud_Min_Pressure  =    100.0 ! ( hPa )
REAL, SAVE :: CloudFree           =    1.e-6 ! Minimum cloud amount
REAL, SAVE :: ThickCloud          =     1.00 ! Maximum cloud amount

!------------------------------------
! 1.7. Emissivity
!------------------------------------

LOGICAL, SAVE :: Use_EmisAtlas = .FALSE.
CHARACTER (LEN=80), SAVE :: Atlas_Dir = 'EmisAtlas'
INTEGER, SAVE :: Atlas_Ver = 100 ! 100 is default, 200 is CNRM atlas

!------------------------------------
!2. Processing and surface type codes
!------------------------------------

!2.1) For RT Models:
!----
INTEGER, PARAMETER :: &
  RTland  = 0,        & !   RTTOV code for land
  RTsea   = 1,        & !   RTTOV code for sea
  RTice   = 2,        & !   RTTOV code for ice
  RTmixed = 3           !   RTTOV code for land and sea ice

CHARACTER(LEN=5), PARAMETER :: &
  RTsurftype_text(0:3) = (/ "land ", "sea  ", "ice  ", "mixed" /)

!2.2) For 1DVar:
!----


INTEGER, PARAMETER :: &
  sea           = 1,  &
  seaice        = 2,  &
  land          = 3,  &
  surf_highland = 4,  &
  surf_mismatch = 5

CHARACTER(LEN=8), PARAMETER :: &
  SurfType_text(5) = (/ "sea     ", "seaice  ", "land    ", &
                        "highland", "mismatch" /)


!3. Processing options/classes/types
!-----------------------------------
! These constants define codes for each of the recognised conditions for
! satellite sounding, and are also used to access the array containing
! the channels for each retrieval type.
INTEGER, PARAMETER :: &
  CloudType_Clear     = 1, &
  CloudType_IrCloudy  = 2, &
  CloudType_MwCloudy  = 3, &
  CloudType_Rain      = 4, &
  CloudType_HighCloud = 5

INTEGER, PARAMETER :: &
  QC_BadRawBT       = 1, &
  QC_BackgroundProf = 2

!Output QC information
INTEGER, PARAMETER :: &
  ObAccepted        = 0, &
  ObNotProcessed    = 2, &
  ObNotConverged    = 1

!4. Retrieval Output File Units
!----------------------------------- 

INTEGER, SAVE :: FileUnit_Retrieved_Profiles
INTEGER, SAVE :: FileUnit_Retrieved_BTs
INTEGER, SAVE :: FileUnit_MinimisationLog
INTEGER, SAVE :: FileUnit_BTMinimisationLog
INTEGER, SAVE :: FileUnit_AMatrix
INTEGER, SAVE :: FileUnit_AmMatrix
INTEGER, SAVE :: FileUnit_ProfileQC
INTEGER, SAVE :: FileUnit_AveragingKernel
INTEGER, SAVE :: FileUnit_Jacobian
INTEGER, SAVE :: FileUnit_RetJacobian


! 5. Namelist
!------------------------------------

NAMELIST / Control / &
     FirstOb,                     &
     LastOb,                      &
     GeneralMode,                 &
     Perform1DVar,                &
     DetectCloud,                 &
     MaxIterations,               &
     DeltaJ,                      &
     SmallJCost_Gradient,         &
     Max_ML_Iterations,           & 
     Gamma_Factor,                &
     GammaMax,                    &
     DoTextrapolation,            & 
     Minimisation_Method,         &
     Allow_Eqn_101,               &
     Force_Eqn_101,               &
     MaxChanUsed,                 &
     CostThresh_Land,             &
     CostThresh_Sea,              &
     CostThresh_IRWindow_Land,    &
     CostThresh_IRWindow_Sea,     &
     CloudAbsThresh_IRWindow,     &
     HighCloudAbsThresh_IRWindow, &
     Additional_Cost_Function,    &
     Cloud_Min_Pressure,          &
     Use_EmisAtlas,               &
     Atlas_Dir,                   &
     Atlas_Ver,                   &
     Retrieve_qtotal,             &
     Read_CLW_Background,         &
     UseRRs,                      &
     EnhancedDiagnostics,         &
     Gas_Units,                   &
     Legacy_Setting


END MODULE NWPSAFMod_Params
