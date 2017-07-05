!+ Loop over the given data records,
!  processing one spot at a time.

Subroutine NWPSAF_ProcessData( &
     Obs,           & ! inout
     Background,    & ! inout
     FirstOb,       & ! in
     LastOb,        & ! in
     RT_Params,     & ! inout
     Total_channels ) ! in

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
! Contains a loop over the given data records (spots), processing them one spot
! at a time. 
!
! Method:
!
! For each record:
!
! (1) Translate the data into a suitable form for processing. i.e. variables
!     inserted into the RT model's profile vector in the correct format.
!     This is carried out in NWPSAF_TranslateDataIn.
!
! (2) Call the radiative transfer model to calculate background brightness
!     temperatures for all appropriate channels.
!
! (3) Cloud parameters are estimated from the background and radiances.
!
! (4) In NWPSAF_CloudyOrNot, test the measured radiances for the presence
!     of cloud. This is done in the Cloud_Cost routine. A threshold test
!     on the cloud cost is used to decide if the spot is clear or cloudy.
!     This is skipped over if the cloud retrieval is in place.
!
! (5) Set up for 1DVAR processing based on the result of NWPSAF_CloudyOrNot.
!     If 1DVar processing is not required, then this and sections 5, 6 
!     are skipped over.
!
! (6) The next stage is the iterative loop, in NWPSAF_Minimize. This is
!     exited on convergence or on reaching the maximum number of iterations.
!
! (7) This profile is converted to the required format for output in 
!     NWPSAF_TranslateDataOut.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     20/05/99 Original code based on ATOVS_ProcessData, stripped
!                  down but with options to add OPS code back in when 
!                  needed.  Andrew D. Collard
! 1.2     25/05/00 Many undocumented changes in the last year.  ADC
! 1.3     25/05/00 RTGuess and RTBack are ALLOCATEd here.  (Moved from 
!                  NWPSAF_SetUpBackground).  They are also DEALLOCATEd
!                  at the end. ADC
! 1.4     09/11/01 Copied NWPSAF_Read_Observations here instead of reading
!                  in all obs at the start which can cause memory 
!                  problems.
! 2.2     13/03/02 Removed some unneccessary arguments in call to 
!                  NWPSAF_Minimize.  
!                  Check for zero channels                   A. Collard.
! 2.3     22/05/02 Removed SatID_Convert (obsolete).         A. Collard.
! 3.0.1   15/07/03 Add M. Szyndel's cloud processing.        A. Collard.
! 3.0.3   01/03/04 Add channel constants to RT_Params.       A. Collard.
! 3.0.4   02/03/04 Rename ProcType, CloudType.               A. Collard.
! 3.0.5   29/03/04 Remove RT_Params % SatViewAngle.          A. Collard.
! 3.0.6   18/06/04 Add RTParams % RT1stGuess.                A. Collard.
! 
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    B_Reference, B_Reference_Inverse, &
    Bmatrix_sea,Bmatrix_land, &
    R_Matrix_Type, &
    R_reference

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type, &
    ModelOb_Type

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type, &
    BackChans,             &
    DetectCloudChans,      &
    ChannelChoice

USE NWPSAFMod_Constants, ONLY : &
    MissData_R

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,         &
    DebugMode,           &
    VerboseMode,         &
    StatusFatal,         &
    StatusWarning,       &
    DetectCloud,         &
    Perform1DVar,        &
    CloudyRetrieval,     &
    EnhancedDiagnostics, &
    ObAccepted,          &
    ObNotProcessed,      & 
    ObNotConverged,      &
    MaxIterations,       &
    UsePCs,              &
    CalcRadiance,        &
  !Surface types:
    RTsea,RTland,RTice,RTsurftype_text,      &
    sea,seaice,land,surf_highland,surf_mismatch, &
  !Cloud processing codes:
    CloudType_Clear,     &
    CloudType_IrCloudy,  &
    CloudType_MwCloudy,  &
    CloudType_Rain,      &
    CloudType_HighCloud, &
  !Various constants
    QC_BackgroundProf, &
    FileUnit_Retrieved_Profiles, &
    FileUnit_Retrieved_BTs,      &
    FileUnit_AMatrix,            &
    FileUnit_AmMatrix,           &
    FileUnit_ProfileQC,          &
    FileUnit_Jacobian,           &
    FileUnit_RetJacobian,        &
    FileUnit_AveragingKernel,    &
    FileUnit_MinimisationLog,    &
    FileUnit_BTMinimisationLog

USE NWPSAFMod_RTmodel, ONLY : &
    RTParams_Type, &
    Num_RTLevels, &
    ProfSize, &
    FastmodelMode_CleanUp, &
    FastmodelMode_Forward, &
    BackGrProf, &
    GuessProf, &
    NeitherProfile

USE rttov_const, ONLY : &
    inst_id_iasi,        &
    inst_id_airs

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_CloudyOrNot.interface'
INCLUDE 'NWPSAF_Fastmodel_Interface.interface'
INCLUDE 'NWPSAF_CO2Slice.interface'
INCLUDE 'NWPSAF_Minimize.interface'
INCLUDE 'NWPSAF_Read_Observations.interface'
INCLUDE 'NWPSAF_SetUpBackground.interface'
INCLUDE 'NWPSAF_TranslateDataIn.interface'
INCLUDE 'NWPSAF_TranslateDataOut.interface'

! Subroutine arguments:

TYPE(Ob_type),       INTENT(INOUT) :: Obs            ! Observed/Retrieval data 
TYPE(ModelOb_type),  INTENT(INOUT) :: Background     ! Background data
INTEGER,             INTENT(IN)    :: FirstOb        ! First ob to process
INTEGER,             INTENT(IN)    :: LastOb         ! Last ob to process
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params      ! Info for RT Model
INTEGER,             INTENT(IN)    :: Total_channels ! Number of channels

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = "NWPSAF_ProcessData"
CHARACTER(LEN=*), PARAMETER :: cutoff = "------------------------------&
                                        &------------------------------&
                                        &--------------------"
CHARACTER(LEN=*), PARAMETER :: Reject_format = &
  "(' Ob',I6,' (',F7.2,',',F7.2,';',A6,'):'"

! Local variables:
                             ! Switches and data flags:
LOGICAL :: Contaminated
LOGICAL :: highcloud
LOGICAL :: MwCloudy
LOGICAL :: IrCloudy
LOGICAL :: Clear
LOGICAL :: DoCloudDetection
LOGICAL :: SurfaceMismatch     !   Flag for NESDIS /= Model surface type
LOGICAL :: HighLand            !   Flag for land surface > specified limit
LOGICAL :: Valid_data
LOGICAL :: DoOneDVar
                                ! Miscellaneous:
REAL    :: BTDifference         !   Brightness temperature difference (K)
INTEGER :: Fastmodel_Mode       !   Mode in which fastmodel is called
INTEGER :: channel              !   Loop variable

!xxx INTEGER :: Band, Bsurf
INTEGER :: CloudType             !index for accessing ChannelChoice
INTEGER :: i
INTEGER :: Bmatrix_surface ! index for B_ErrCovar (depends on surface type)
INTEGER :: obnumber
INTEGER :: surftype
INTEGER :: ErrorCode
INTEGER :: ProfileQC
INTEGER :: RTerrorcode
REAL, POINTER :: MeasurementCalc(:)        !Forward Modelled BTs/Rads/PCs
REAL    :: MeasurementObs(Total_channels)  !Observed (bias corrected) BTs/Rads/PCs
LOGICAL :: GoodChannels(Total_channels)
LOGICAL :: GoodBackChannels(Total_channels)   !Channel data flags
CHARACTER(LEN=160) :: Output_format
CHARACTER(LEN=80)  :: Message(2)
CHARACTER(LEN=80)  :: Output_text

REAL, POINTER :: B_Matrix(:,:)
REAL, POINTER :: B_Matrix_Inverse(:,:)
TYPE(R_Matrix_Type), POINTER :: R_Matrix
TYPE(ChannelSelection_Type) :: UsedChans
  
INTEGER :: WhichProf

!------------------------------------------------------------------------------

!-------------------
! 1. Initialize
!-------------------

!1.1) Initialize variables and allocate array space
!----
Nullify(RT_Params % PCscores)
Nullify(RT_Params % TotalBTs)
Nullify(RT_Params % TotalRadiances)
IF ( UsePCs ) THEN
  ALLOCATE(RT_Params % PCScores(Total_Channels))
  MeasurementCalc=> RT_Params % PCscores
ELSE IF ( CalcRadiance ) THEN
  ALLOCATE(RT_Params % TotalRadiances(Total_Channels))
  MeasurementCalc=> RT_Params % TotalRadiances
ELSE
  ALLOCATE(RT_Params % TotalBTs(Total_Channels))
  MeasurementCalc=> RT_Params % TotalBTs
ENDIF

IF (CloudyRetrieval) THEN
  ALLOCATE(RT_Params % TotalRadiances(Total_Channels))
  ALLOCATE(RT_Params % CloudyRadiances(Total_Channels,2*Num_RTLevels+2))
  ALLOCATE(RT_Params % tc1(Total_Channels))
  ALLOCATE(RT_Params % tc2(Total_Channels))
  ALLOCATE(RT_Params % bcon1(Total_Channels))
  ALLOCATE(RT_Params % bcon2(Total_Channels))
END IF
NULLIFY(R_Matrix)

Output_text = ''
Message(:) = ''

! Allocate RT_Params % RTBack and RT_Params % RTGuess for the fastmodel 
! interface and RT_Params % RT1stGuess so we can keep track of the
! 1st guess where it differs from the background (e.g., clody retrievals)

ALLOCATE(RT_Params % RTBack(ProfSize))
ALLOCATE(RT_Params % RT1stGuess(ProfSize))
ALLOCATE(RT_Params % RTGuess(ProfSize))

! Other data
RT_Params % SatZenithAngle    = 0.0
RT_Params % SatAzimAngle      = 0.0
RT_Params % SolarZenAngle     = 0.0
RT_Params % SatSolarAzimAngle = 0.0

!--------------------
! 2. Process profiles
!--------------------

Loop: DO obnumber = FirstOb, LastOb

  WRITE(*,*) 'Processing Ob number ',obnumber

  !Initialize data flags
  SurfaceMismatch = .FALSE.
  HighLand        = .FALSE.
  Contaminated    = .FALSE.
  MwCloudy        = .FALSE.
  highcloud       = .FALSE.
  IrCloudy        = .FALSE.
  Clear           = .FALSE.
  Valid_data      = .TRUE.

  !Initialize switches
  DoOneDVar        = Perform1DVar
  DoCloudDetection = DetectCloud

  IF (UsePCs .and. DoCloudDetection) THEN
    DoCloudDetection = .false.
    Message(1)='Cloud detection only works for Brightness Temps / Radiances for now'
    Message(2)='Switching off Cloud Detection'
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      Message,                & ! in
      ErrorStatus=StatusWarning ) ! in 
    Message(:)=''
  END IF

  !Initialize QC information
  ProfileQC       = ObAccepted
  !--------------------------------
  ! 2.1 Read in next Ob
  !--------------------------------

  CALL NWPSAF_Read_Observations ( &
    Obnumber,     & ! in
    Background,   & ! inout
    Obs,          & ! inout
    RT_Params)      ! inout

  !2.1.1) Check setup if PC assimilation requested
  IF ( UsePCs .and.  &
       .not. size(RT_Params % SubTypeChoice(:)) == 1        .and. &
       .not. ( RT_Params % SubTypeChoice(1) == inst_id_iasi .or.  &
               RT_Params % SubTypeChoice(1) == inst_id_airs )     ) THEN
    Message(1)='PC Assimilation not possible for chosen instrument(s)'
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      Message,                & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  ! --------------------------------
  ! 2.2 Process satellite identifier
  ! --------------------------------

  DO I=1,RT_Params % Num_SatIDs
    IF (RT_Params % SatID(I) % WMO == Obs % SatID) THEN
      RT_Params % SatIndex = I
      EXIT
    END IF
    IF (I == RT_Params % Num_SatIDs) THEN
      Message(1) = 'Required satellite cannot be processed'
      WRITE(Message(2),FMT='(A,I4,A,I6)') &
          'Required WMO Satellite number is ',Obs % SatID,&
          ', Obnumber is ',obnumber
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        Message,                & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO

  ! -------------------------
  ! 2.3 Identify surface type
  ! -------------------------
  
  ! Surface type is taken from the observation file.  
  ! Set surftype and the appropriate RTSurfaceType in each case.

  surftype = Obs% Surface
  SELECT CASE( surftype )
  CASE(Sea)
    RT_Params % RTSurfaceType = RTsea 
  CASE(SeaIce)
    RT_Params % RTSurfaceType = RTIce 
  CASE(Land)
    RT_Params % RTSurfaceType = RTLand
  CASE(Surf_HighLand)
    RT_Params % RTSurfaceType = RTLand 
  CASE(Surf_Mismatch)
    RT_Params % RTSurfaceType = RTLand 
  CASE DEFAULT
    RT_Params % RTSurfaceType = RTLand 
  END SELECT

  ! ------------------------
  ! 2.4 Transform Input Data
  ! ------------------------

  !Initialise
  MeasurementCalc  = MissData_R

  !Set up background / observation data.

  CALL NWPSAF_TranslateDataIn( &
    Obs,               & !in
    RT_Params,         & !inout
    MeasurementObs(:), & !out
    GoodChannels(:),   & !out
    ErrorCode          ) !out

  !Set BTs to missing data and skip to next observation if data is bad
  IF ( ErrorCode /= 0 ) THEN
    write(*,*) 'Bad observation data, ErrorCode =',ErrorCode
    Valid_Data = .FALSE.
    MeasurementCalc(:) = MissData_R
    ProfileQC       = ObNotProcessed
    Write(FileUnit_ProfileQC, '(2I8)') obnumber, ProfileQC
    CYCLE Loop
  END IF

  !--
  ! 2.5) Set up background profile information
  !--
  
  CALL NWPSAF_SetupBackground( &
    BackGround,      & ! in
    obnumber,        & ! in
    RT_Params,       & ! in
    Valid_Data )       ! out

  ! This has been moved from NWPSAF_Minimise as we want to use the guess
  ! values for CTP and CloudFrac in the RT model calculations rather than
  ! background.
  RT_Params % RTguess = RT_Params % RTBack(:)
  RT_Params % RT1stguess = RT_Params % RTGuess(:)

  IF ( .NOT. Valid_Data ) ErrorCode = QC_BackgroundProf

  ! --------------------------------------
  ! 2.6 Background brightness temperatures
  ! --------------------------------------

  !2.6.1) Set channels for forward model
  !------

  UsedChans % NumChans = BackChans % numchans
  UsedChans % Channels => BackChans % channels(1 : BackChans % numchans)

  !2.6.2) Calculate background brightness temperatures
  !------

  WhichProf = GuessProf ! This is the background profile unless 
                        ! CloudyRetrieval.
  Fastmodel_Mode = FastmodelMode_Forward

  CALL NWPSAF_Fastmodel_interface( &
    Fastmodel_Mode,                & ! in
    RT_Params,                     & ! inout
    WhichProf,                     & ! in
    UsedChans,                     & ! in
    RTErrorCode)                     ! out

  !2.6.3) Set error flags
  !------
  IF ( RTErrorCode > 0 ) THEN
    IF ( GeneralMode >= DebugMode ) THEN
      WRITE( Output_Text, FMT='(A,I3)' ) &
          ' error returned from RT Model =',RTErrorCode
      Output_format = Reject_format // ",'" // Output_text // "')"
      WRITE(*,Output_format) Obs% Id, &
          Obs % latitude % value, &
          Obs % longitude % value, &
          RTsurftype_text(RT_Params % RTsurfacetype)
      ProfileQC       = ObNotProcessed
      Write(FileUnit_ProfileQC, '(2I8)') obnumber, ProfileQC
    END IF
    IF ( RTErrorCode >= 20 ) THEN
      MeasurementCalc = MissData_R
      CYCLE
    ELSE
      DoCloudDetection = .FALSE.
      DoOneDVar        = .FALSE.
    END IF
  END IF

   ! ------------
   ! 2.7) Set up the R-matrix
   ! ------------
   IF (.NOT.(RT_Params % SatID(RT_Params % SatIndex) % R_Matrix_Present)) THEN
     DoCloudDetection = .FALSE.
     DoOneDVar = .FALSE.
     WRITE(*,*) 'Missing R-matrix for obnumber',obnumber
     WRITE(*,*) RT_Params % SatIndex
     CYCLE
   ELSE
     !Select correct R_matrix
     R_matrix => R_reference(RT_Params % SatIndex)
   END IF


  !------
  !2.8) Obs minus background check
  !------

  GoodBackChannels(:) = .FALSE.
  DO i = 1, UsedChans % NumChans
    channel = Usedchans % Channels(i)
    BTdifference = MeasurementCalc(i) - MeasurementObs(channel)
    IF (CloudyRetrieval) THEN
      ! For cloudy retrievals, the channel brightness difference test is 
      ! relaxed as the background at this point is cloud free, whilst the 
      ! observation may well be cloudy.
      IF ( ABS( BTdifference ) <= 100.0 ) GoodBackChannels(channel) = .TRUE.
    ELSE IF (UsePCs) THEN
      ! For PC Retrieval, if any given PC is too far from the threshold, switch off
      ! 1D-Var. Who knows what would happen!
      IF ( i == 1 .and. ABS(BTdifference) > R_matrix % Diagonal(channel)*50.0 ) THEN
      ! ************EXPERIMENTAL SETTING************
      ! Alternatively, could try something like
      ! IF ( i == 1 .and. ABS(BTdifference) > 2000.0 ) THEN
        DoOneDVar = .FALSE.
        IF ( GeneralMode >= DebugMode ) THEN
          WRITE( UNIT=Message(1), FMT='(A,F10.3,A,I5,A,F10.3)' ) &
                ' Switching off 1DVar, O-B too large: ',BTdifference
          WRITE( UNIT=Message(2), FMT='(A,I5,A,F10.3)' ) &
              ' for PC number ', channel, ' Threshold: ', &
                R_matrix % Diagonal(channel)*50.0
          WRITE(*,*) Message(1)
          WRITE(*,*) Message(2)
          ProfileQC       = ObNotProcessed
          Write(FileUnit_ProfileQC, '(2I8)') obnumber, ProfileQC
        END IF
        EXIT
      ELSE
        GoodBackChannels(channel) = .TRUE.
      END IF
      GoodBackChannels(channel+1:UsedChans%NumChans) = .TRUE.
    ELSE
      ! Normal case
      IF ( ABS( BTdifference ) <= 20.0 ) GoodBackChannels(channel) = .TRUE.
    END IF
  END DO
  
  Background% BriTemp(Usedchans % Channels) = &
                                         MeasurementCalc(1:UsedChans % NumChans)
 
   ! ------------
   ! 2.9 Set up the B-Matrix
   ! ------------
  
  IF ( DoCloudDetection .OR. DoOneDVar .OR. CloudyRetrieval) THEN
    SELECT CASE( RT_Params % RTSurfaceType )
    CASE(RTsea)
      Bmatrix_surface = Bmatrix_sea
    CASE(RTice,RTland)
      Bmatrix_surface = Bmatrix_land
    CASE DEFAULT
      SurfaceMismatch = .TRUE.
      RT_Params % RTSurfaceType = RTland
      Bmatrix_surface = Bmatrix_land
    END SELECT
    B_Matrix         => B_Reference(:,:,Bmatrix_surface)
    B_Matrix_Inverse => B_Reference_Inverse(:,:,Bmatrix_surface)
  END IF


  IF (CloudyRetrieval) THEN
    
    !2.9.1) For cloudy retrievals, determine first guess and background of 
    !       cloud top pressure and fraction using CO2 slicing and relative 
    !       humidity thresholding respectively.
    !       Brightness temperatures are also recalculated from the overcast 
    !       radiances to correspond to the first guess (i.e. the CO2 slice 
    !       cloud parameters)
    !------ 
  
    CALL NWPSAF_CO2Slice(  &
      MeasurementObs(:),   &  ! in
      RT_Params,           &  ! inout
      R_Matrix,            &  ! in
      UsedChans,           &  ! in
      Total_Channels       )     ! in
  
    !2.9.2) Transfer recalculated values into 1DVar variables
    !------
  
    GoodBackChannels(:) = .FALSE.
    DO i = 1, UsedChans % NumChans
      channel = Usedchans % Channels(i)
      BTdifference = MeasurementCalc(i) - MeasurementObs(channel)
      IF ( ABS( BTdifference ) <= 20.0 ) GoodBackChannels(channel) = .TRUE.
    END DO

    Background% BriTemp(Usedchans % Channels) = &
                                         MeasurementCalc(1:UsedChans % NumChans)

  END IF

  ! --------------------------
  ! 2.10 Select processing type
  ! --------------------------


  !2.10.1) Calculate a cloud cost function
  !------
  DetectCloudBlock: IF ( DoCloudDetection ) THEN
    UsedChans % NumChans    =  DetectCloudChans % numchans
    UsedChans % Channels => &
          DetectCloudChans % channels(1 : UsedChans % NumChans)
    WhichProf = BackGrProf
    CALL NWPSAF_CloudyOrNot( &
      Obs,                & ! in
      MeasurementObs(:),  & ! in
      B_Matrix,           & ! in
      R_Matrix,           & ! in
      UsedChans,          & ! in
      RT_Params,          & ! in        
      Valid_data,         & ! inout
      IrCloudy,           & ! inout
      HighCloud,          & ! inout
      RTerrorcode         ) ! out
    IF ( RTerrorcode > 0 ) THEN
      IF ( GeneralMode >= DebugMode ) THEN
        WRITE( UNIT=Output_text, FMT='(A,I3)' ) &
          ' error from RT model (evaluating cloud cost) = ',RTerrorcode
        Output_format = Reject_format // ",'" // Output_text // "')"
      END IF
    ELSEIF ( IrCloudy ) THEN
      Clear = .FALSE.  
    ELSE
      Clear = .TRUE.
    END IF
  ELSE
    IrCloudy = .TRUE.
  END IF DetectCloudBlock


  !2.10.2) Select processing type
  !------
  ! Processing options are chosen with an order of preference that
  ! roughly correlates with the degree of data 'contamination' (in the
  ! event that some spots will belong to more than one category).
  ! This generally just affects channel selection
  
  IF ( SurfaceMismatch ) THEN
    surftype = surf_mismatch
  ELSEIF ( HighLand ) THEN
    surftype = surf_highland
  END IF
  IF ( Contaminated ) THEN
    CloudType = CloudType_Rain
  ELSEIF ( MwCloudy ) THEN
    CloudType = CloudType_MwCloudy
  ELSEIF ( highcloud ) THEN
    CloudType = CloudType_HighCloud
  ELSEIF ( IrCloudy ) THEN
    CloudType = CloudType_IrCloudy
  ELSEIF ( Clear ) THEN
    CloudType = CloudType_Clear
  ELSE
    CloudType = CloudType_IrCloudy
  END IF
  
  !2.10.3) Set choice of channels for iterations
  !------

  IF ( ChannelChoice(surftype,CloudType) % numchans > 0 ) THEN
    UsedChans % NumChans = ChannelChoice(surftype,CloudType) % numchans
    UsedChans % Channels => &
      ChannelChoice(surftype,CloudType) % channels(1:UsedChans % NumChans)
    IF ( GeneralMode >= VerboseMode ) THEN
      WRITE(*,*) 'Choice of processing channels: ', &
        UsedChans % Channels(1 : UsedChans % NumChans)
      WRITE(*,*)
    END IF
  ELSE
    CYCLE
  END IF

  !reject observation if process channels flagged
  Check_channels: DO channel = 1, UsedChans % NumChans
    IF ( .NOT. GoodBackChannels(UsedChans % Channels(channel)) .OR. &
      ChannelChoice(surftype,CloudType) % numchans == 0) THEN
      Valid_data = .FALSE.
      IF ( GeneralMode >= DebugMode ) THEN
        WRITE( UNIT=Output_text,FMT='(A,I5,A)' ) &
          ' channel',UsedChans % Channels(channel), ' failed o-b check'
        WRITE(*,*) Output_Text
        ProfileQC       = ObNotProcessed
        Write(FileUnit_ProfileQC, '(2I8)') obnumber, ProfileQC
      END IF
      EXIT Check_channels
    END IF
  END DO Check_channels
  
  ! ---------------------------
  ! 2.11 Perform 1DVar retrieval
  ! ---------------------------

  Perform1DVarBlock: IF ( DoOneDVar .AND. Valid_data ) THEN
    
    RTErrorCode=0

    CALL NWPSAF_Minimize( &
      Obs,                & ! inout
      MeasurementObs(:),  & ! in
      B_matrix,           & ! in
      B_matrix_Inverse,   & ! in
      R_Matrix,           & ! in
      obnumber,           & ! in
      UsedChans,          & ! in
      RT_Params,          & ! inout
      RTerrorcode         ) ! out

   
    ! Reset channel choices in case they were reset in NWPSAF_Minimize
    !  (not currently possible)

    UsedChans % NumChans = ChannelChoice(surftype,CloudType) % numchans
    UsedChans % Channels => &
      ChannelChoice(surftype,CloudType) % channels(1:UsedChans % NumChans)
    
    !2.10.1) Set error flags
    !------
    
    IF ( RTerrorcode /= 0 ) THEN
      Valid_data = .FALSE.
      IF ( GeneralMode >= DebugMode ) THEN
        WRITE( UNIT=Output_text, FMT='(A,I3)' ) &
              ' error from RTTOVK = ',RTerrorcode
        WRITE(*,*) Output_Text
      END IF
    END IF
    
    ! --------------------------------------
    ! 2.11 QC and tidy up results for output
    ! --------------------------------------
    IF ( Valid_data ) THEN
        
      !2.11.1) Output retrieval data 
      !-------
      
      CALL NWPSAF_TranslateDataOut( &
        obnumber,       &   ! in
        MeasurementObs, & ! in
        Background,     & ! in
        RT_Params,      & ! in
        UsedChans,      & ! in
        Obs)              ! inout

      IF ( Obs % NIter > MaxIterations ) ProfileQC=ObNotConverged

    ELSE

      ProfileQC=ObNotProcessed

    END IF

    Write(FileUnit_ProfileQC, '(2I8)') obnumber, ProfileQC
    
  END IF Perform1DVarBlock

  IF ( GeneralMode >= VerboseMode ) THEN
    WRITE(*,'(A)') 'Brightness temperatures:'
    WRITE(*,*)
    WRITE(*,'(A)') 'Channel      Obs     Back      Ret      O-R      O-B'
    DO i = 1, UsedChans % NumChans
      channel = UsedChans % Channels(i)
      WRITE(*,'(I7,5F9.2)') channel, &
        MeasurementObs(channel),Background%BriTemp(channel), &
        MeasurementCalc(i), &
        MeasurementObs(channel)-MeasurementCalc(i), &
        MeasurementObs(channel)-Background%BriTemp(channel)
    END DO
  END IF
   
END DO Loop

! ----------
! 3. Tidy up
! ----------

!3.1) Deallocate arrays allocated in ProcessData
!----
DEALLOCATE(RT_Params % RTBack)
DEALLOCATE(RT_Params % RTGuess)
DEALLOCATE(RT_Params % RT1stGuess)
IF (ASSOCIATED(RT_Params % TotalBTs) ) DEALLOCATE(RT_Params % TotalBTs)
IF (ASSOCIATED(RT_Params % TotalRadiances) ) DEALLOCATE(RT_Params % TotalRadiances)
IF (ASSOCIATED(RT_Params % PCScores) ) DEALLOCATE(RT_Params % PCScores)
IF (CloudyRetrieval) THEN
  DEALLOCATE (RT_Params % CloudyRadiances)
  DEALLOCATE(RT_Params % tc1)
  DEALLOCATE(RT_Params % tc2)
  DEALLOCATE(RT_Params % bcon1)
  DEALLOCATE(RT_Params % bcon2)
END IF
NULLIFY(UsedChans % Channels)
NULLIFY(R_Matrix)
NULLIFY(B_Matrix)
NULLIFY(B_Matrix_Inverse)

!3.2) Deallocate arrays in Fastmodel
!----

Fastmodel_Mode = FastmodelMode_CleanUp
WhichProf = NeitherProfile
UsedChans % NumChans = 0
CALL NWPSAF_Fastmodel_interface( &
  Fastmodel_Mode,          & ! in
  RT_Params,               & ! in
  WhichProf,               & ! in
  UsedChans,               & ! in
  RTErrorCode)               ! out

!3.3) Close Output Files
!----

CLOSE(FileUnit_ProfileQC)
CLOSE(FileUnit_Retrieved_Profiles)
CLOSE(FileUnit_Retrieved_BTs)

IF (GeneralMode >= DebugMode) THEN
  CLOSE(FileUnit_MinimisationLog)
  CLOSE(FileUnit_BTMinimisationLog)
  IF (EnhancedDiagnostics) THEN
    CLOSE(FileUnit_AMatrix)
    CLOSE(FileUnit_AmMatrix)
    CLOSE(FileUnit_AveragingKernel)
    CLOSE(FileUnit_Jacobian)
    CLOSE(FileUnit_RetJacobian)
  END IF
END IF
   
!3.3) Exit message
!----
IF ( GeneralMode >= DebugMode ) THEN
   Message(1) = 'Finished processing'
   IF ( GeneralMode >= DebugMode ) THEN
      Message(2) = cutoff
   ELSE
      Message(2) = ' '
   END IF
   CALL NWPSAF_Report( RoutineName,Message(1:2) )
END IF
   
END SUBROUTINE NWPSAF_ProcessData
