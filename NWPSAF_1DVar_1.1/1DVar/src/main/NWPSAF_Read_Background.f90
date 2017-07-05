SUBROUTINE NWPSAF_Read_Background ( &
     NumObs,       & ! in
     NumLevels,    & ! out
     Background)     ! inout

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
!------------------------------------------------------------------------------
! Description: Read in background (a priori) data
!
! Method:
!    Reads in a priori data from file and assigns to Background.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.1      15/07/99  Original Version
! 1.2      02/03/01  Added ozone profile R. Saunders
! 2.2      30/04/02  Removed the need for explicit formats in input files.
!                    (Calls new subroutine NWPSAF_ReadHeaders).  ADC.
! 2.3      27/05/02  Remove status component of ModelOB structure.  ADC.
! 3.0.1    17/07/03  Allow reading in of cloud parameters.          ADC.
! 3.0.2    30/09/03  Allow cloud parameters to be optional.         ADC.
! 3.0.4    04/03/04  Read number of levels and pressures here.   A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date      Comment
! ------  --------  -------
! 19      30/03/10  Fix problems with swapping array order when compiling
!                   with GFortran. Fiona Hilton
! 28      22/02/12  Include changes to read in cloud liquid water profile
!                   TR Sreerekha
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Constants, ONLY : &
     MissData_R

USE NWPSAFMod_ObsInfo, ONLY : &
     ElementHeader_type,    &
     ModelOb_type

USE NWPSAFMod_Params, ONLY: &
     StatusFatal,         &
     StatusWarning,       &
     Humidity_Units,      &
     CloudyRetrieval,     &
     GeneralMode,         &
     VerboseMode,         &
     MaxLevs,             &
     Read_CLW_Background, &
     Ozone_Present

USE NWPSAFMod_RTmodel, ONLY : &
     ProfSize,     &
     Prof_FirstT,  &
     Prof_LastT,   &
     Prof_FirstQ,  &
     Prof_LastQ,   &
     Prof_FirstO3, &
     Prof_LastO3,  &
     Prof_FirstCLW,&
     Prof_LastCLw, &
     Prof_LWP,     &
     Prof_T2,      &
     Prof_q2,      &
     Prof_Tstar,   &
     Prof_pstar,   &
     Prof_uwind,   &
     Prof_vwind,   &
     Prof_CTP,     &
     Prof_CloudCover

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_OpenFile.interface'
INCLUDE 'NWPSAF_ReadHeaders.interface'

! Subroutine Arguments 

INTEGER, INTENT(IN)  :: NumObs         ! # of soundings in Observation structure
INTEGER, INTENT(OUT) :: NumLevels      ! # of vertical levels
TYPE(ModelOB_type), INTENT(INOUT) :: Background    ! BackGround at Ob

! Declarations

INTEGER :: I, J

INTEGER :: fileunit              ! I/O unit number
INTEGER :: NumProf               ! # of profile in Background file
INTEGER :: NumHeaderLines = 10   ! # of lines of informational header at 
                                 ! start of background file
INTEGER :: NumSubHeaderLines = 3 ! # of lines of informational header in 
                                 ! background file for each profile
INTEGER :: ReadStatus            ! IOSTAT keyword value

CHARACTER(LEN=70) :: ErrorMessage(2)  ! Message for NWPSAF_Report
CHARACTER(LEN=*),  PARAMETER :: RoutineName = 'NWPSAF_Read_Background'
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER(LEN=80) :: Header                    ! Header for ObsFileName
CHARACTER(LEN=80) :: ObsFileName = "Background.dat"
CHARACTER(LEN=20) :: LevNumChar  ! Character Representation of Level Number
CHARACTER(LEN=20) :: ProfNumChar ! Character Representation of Profile Number
CHARACTER(LEN=3)  :: Status = "OLD"
CHARACTER(LEN=80) :: String = " "              ! Header for ObsFileName

TYPE(ElementHeader_type) :: Dummy_header  ! holds initialization data

REAL :: Vars(3)

REAL :: temparray(maxlevs) 

LOGICAL :: convert = .false.

!-----------------------------------------------------------------------------

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

ErrorMessage(:) = ''

Dummy_header % NumLev = 0

!-----------------------------------
!1. Initialize BackGround structure
!-----------------------------------


!-------------------------------
!1.1 Set up some variable headers
!-------------------------------

Background % header % NumObsLocal = NumObs

!--------------------------------------------------------
!1.2 Open the Observations File and Read in Header Data
!--------------------------------------------------------

CALL NWPSAF_OpenFile( &
  Trim(ObsFileName),   & ! in
  Access,              & ! in
  Action,              & ! in
  Status,              & ! in
  Form,                & ! inout
  fileunit )             ! out

! Allow NumHeaderLines Lines for an informational header

ReadHeader : DO i = 1, NumHeaderLines
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (ReadStatus /= 0) THEN
    Errormessage(1) = ' Error Reading Background File Header '
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
END DO ReadHeader

! Read in number of Background Profiles Present
! If there is only one, use it to process all observations
! If there are NumObs observations, each observation has its own background.
! Else report a fatal error


READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (ReadStatus == 0) &
  CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
NumProf = NINT(Vars(1))
IF (ReadStatus /= 0 .OR. (NumProf /= 1 .AND. NumProf /= NumObs)) THEN
  Errormessage(1) = ' Error Reading Background File: No. of Profiles '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 

END IF

! Read in number of levels per profile
! Note: Number of levels should be the same for each observation and
! should match the number of levels in the RT model.

READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (ReadStatus == 0) &
     CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
NumLevels = NINT(Vars(1))
IF (ReadStatus /= 0 .OR. NumLevels > 1.E6 .OR. NumLevels <= 0.) THEN
  Errormessage(1) = ' Error Reading Background File: No. of Levels '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 

END IF
IF ( NumLevels > maxlevs) THEN
  Errormessage(1) = ' Variable temparray is not big enough for full profile '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

! Read in the units used for humidity.  The units used during minimisation 
! are log(kg/kg).
READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (ReadStatus == 0) &
     CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
Humidity_Units = NINT(Vars(1))
IF (ReadStatus /= 0 .OR. Humidity_Units<=0 .OR. Humidity_Units>3 ) THEN
   Errormessage(1) = &
        ' Error Reading Background File: Unsupported Humidity Units '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

! Allocate relevent parts of the BackGround structure and initialise 
! associated header elements.

ALLOCATE( Background % Pstar ( NumObs ))
ALLOCATE( Background % rh2 ( NumObs ))
ALLOCATE( Background % t2 ( NumObs ))
ALLOCATE( Background % Tskin ( NumObs ))
ALLOCATE( Background % u10 ( NumObs ))
ALLOCATE( Background % v10 ( NumObs ))
ALLOCATE( Background % CTP ( NumObs ))
ALLOCATE( Background % CldFrac ( NumObs ))
ALLOCATE( Background % LWP ( NumObs ))

Background % Pstar(:) = MissData_R
Background % rh2(:)   = MissData_R
Background % t2(:)    = MissData_R
Background % Tskin(:) = MissData_R

Background % header % Pstar = Dummy_header
Background % header % rh2   = Dummy_header
Background % header % t2    = Dummy_header
Background % header % Tskin = Dummy_header

ALLOCATE( Background % p ( NumObs, NumLevels ))
ALLOCATE( Background % rh ( NumObs, NumLevels ))
ALLOCATE( Background % t ( NumObs, NumLevels ))
ALLOCATE( Background % ozone ( NumObs, NumLevels ))
ALLOCATE( Background % clw ( NumObs, NumLevels ))

Background % header % rh           = Dummy_header
Background % header % rh % NumLev  = NumLevels
Background % header % t            = Dummy_header
Background % header % t %  NumLev = NumLevels
Background % header % ozone        = Dummy_header
Background % header % ozone % NumLev = NumLevels
Background % header % clw           = Dummy_header
Background % header % clw % NumLev = NumLevels
!------------------------------------------------------------
!1.3 Loop through observations and read in elements
!------------------------------------------------------------

!
! Read first sub-header line here (this is to allow the last line in the
! profile entry to be tested before being used).
!
READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
IF (ReadStatus /= 0) THEN
  WRITE (ProfNumChar, FMT='(I20)') I 
  Errormessage(1) = ' Error Reading Background File SubHeader #1 '
  CALL NWPSAF_Report( &
  RoutineName,            & ! in
  ErrorMessage,           & ! in
  ErrorStatus=StatusFatal ) ! in 
END IF

ReadData : DO I = 1, NumProf

  ReadSubHeader: DO J = 1,NumSubHeaderLines-1
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
    IF (ReadStatus /= 0) THEN
      WRITE (ProfNumChar, FMT='(I20)') I 
      Errormessage(1) = ' Error Reading Background File SubHeader #'&
        //Trim(ProfNumChar)
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO ReadSubHeader
  
  ReadProfile: DO J = 1, NumLevels
    IF (Read_CLW_Background) THEN
      READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) &
            Background % p(I,J), & 
            Background % t(I,J), & 
            Background % rh(I,J),&
            Background % ozone(I,J),&
            Background % clw(I,J)
    ELSE
      READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) &
            Background % p(I,J), & 
            Background % t(I,J), & 
            Background % rh(I,J),&
            Background % ozone(I,J)
    ENDIF
    IF (ReadStatus /= 0 .OR. Background % t(I,J) < 0. .OR. &
        Background % rh(I,J) < 0.) THEN
      WRITE (LevNumChar, FMT='(I20)') J 
      WRITE (ProfNumChar, FMT='(I20)') I 
      Errormessage(1) = ' Error Reading Profile #'//Trim(ProfNumChar)//&
            ' at level '//Trim(LevNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
    ENDIF
  END DO ReadProfile
  
  ! Check that the level order is correct (smallest pressure first) and
  ! reverse if need be
  ! Convert pressure levels from hPa to Pa if largest pressure is < 2000.
  IF ( MAXVAL(Background % P(I,:) ) < 2000.0 ) THEN
    Background % P(I,:)=Background % P(I,:)*100.0
    convert=.true.
  END IF

  IF (Background % P(I,1) > Background % P(I,2)) THEN
    temparray(1:NumLevels)=Background % p(I,1:NumLevels)
    Background % p(I,1:NumLevels) = temparray(NumLevels:1:-1)
    temparray(1:NumLevels)=Background % t(I,1:NumLevels)
    Background % t(I,1:NumLevels) = temparray(NumLevels:1:-1)
    temparray(1:NumLevels)=Background % rh(I,1:NumLevels)
    Background % rh(I,1:NumLevels) = temparray(NumLevels:1:-1)
    temparray(1:NumLevels)=Background % ozone(I,1:NumLevels)
    Background % ozone(I,1:NumLevels) = temparray(NumLevels:1:-1)
    IF (Read_CLW_Background) THEN
      temparray(1:NumLevels)=Background % clw(I,1:NumLevels)
      Background % clw(I,1:NumLevels) = temparray(NumLevels:1:-1)
    ENDIF
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  Background % t2(I) = Vars(1)
      
  IF (ReadStatus /= 0) THEN
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading Surface Temperature for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  Background % rh2(I) = Vars(1)
  IF (ReadStatus /= 0) THEN
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading Surface Humidity for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  Background % TSkin(I) = Vars(1)
  IF (ReadStatus /= 0) THEN
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading Skin Temperature for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
   CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  !Convert from hPa to Pa if necessary
  Background % PStar(I) = Vars(1)
  If ( Convert ) Background % PStar(I) = Background % PStar(I) *100.0
  IF (ReadStatus /= 0) THEN 
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading Surface Pressure for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  Background % u10(I) = Vars(1)
  IF (ReadStatus /= 0) THEN
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading 10m U-wind for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
       CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  Background % v10(I) = Vars(1)
  IF (ReadStatus /= 0) THEN
    WRITE (ProfNumChar, FMT='(I20)') I 
    Errormessage(1) = ' Error Reading 10m V-wind for Profile #'&
         //Trim(ProfNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) & 
   CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  IF (ReadStatus /= 0) THEN
    IF (GeneralMode >= VerboseMode) THEN
      WRITE (ProfNumChar, FMT='(I20)') I 
      Errormessage(1) = ' Error Reading CTP for Profile #'&
            //Trim(ProfNumChar)
      Errormessage(2) = ' Cloud Fraction Set to zero and '// &
            'CTP set to surface pressure'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusWarning ) ! in 
    END IF
    Background % CTP(I) = Background % PStar(I)
    Background % CldFrac(I) = 0.0
  ELSE
    ! Read in as hPa and convert to Pa
    Background % CTP(I) = Vars(1)
    IF ( Convert ) Background % CTP(I) = Background % CTP(I) *100.0 
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
    IF (ReadStatus == 0) &
      CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
    IF (ReadStatus /= 0) THEN
      IF (GeneralMode >= VerboseMode) THEN
        WRITE (ProfNumChar, FMT='(I20)') I 
        Errormessage(1) = ' Error Reading Cloud Fraction for Profile #'&
            //Trim(ProfNumChar)
        Errormessage(2) = ' Cloud Fraction Set to zero'
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          ErrorMessage,           & ! in
          ErrorStatus=StatusWarning ) ! in 
      END IF
      Background % CldFrac(I) = 0.0
    END IF
    Background % CldFrac(I) = Vars(1)
    IF (I < NumProf) THEN
      READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
      IF (ReadStatus /= 0) THEN
        WRITE (ProfNumChar, FMT='(I20)') I 
        Errormessage(1) = ' Error Reading Background File SubHeader #'&
            //Trim(ProfNumChar)
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          ErrorMessage,           & ! in
          ErrorStatus=StatusFatal ) ! in 
      END IF
    END IF
  END IF


  ! For CO2Slice to work background must be cloud free
  IF (CloudyRetrieval) THEN
    Background % CldFrac(I) = 0.0
    Background % CTP(I) = Background % PStar(I)
  END IF

  IF (NumProf == 1) THEN
    DO J=1,NumLevels
      Background % P(:,J)  = Background % P(1,J)
      Background % t(:,J)  = Background % t(1,J)
      Background % rh(:,J) = Background % rh(1,J)
      Background % ozone(:,J) = Background % ozone(1,J)
      IF (Read_CLW_Background) THEN
        Background % clw(:,J)   = Background % clw(1,J)
      ELSE
        Background % clw(:,J)   = 0
      ENDIF
    END DO
    Background % t2(:)      = Background % t2(1)
    Background % rh2(:)     = Background % rh2(1)
    Background % Tskin(:)   = Background % Tskin(1)
    Background % PStar(:)   = Background % PStar(1)
    Background % u10(:)     = Background % u10(1)
    Background % v10(:)     = Background % v10(1)
    Background % CTP(:)     = Background % CTP(1)
    Background % CldFrac(:) = Background % CldFrac(1)
  END IF

END DO ReadData

! Set Ozone_present if data in ozone profile is greater than 0

IF (ANY(Background % ozone(:,:) > 0)) Ozone_Present=.TRUE.

!-----------------------------------------------------------------
! 2. Set up positions in RT_Params % RTBack/RTGuess arrays.  These
!    are set up here as they are dependent on the number of levels
!    and are needed in the NWPSAF_SetUpRetrievals subroutine.
!-----------------------------------------------------------------

!Temperature
Prof_FirstT     = 1
Prof_LastT      = NumLevels
! Humidity
Prof_FirstQ     = NumLevels + 1
Prof_LastQ      = 2*NumLevels
! Ozone
Prof_FirstO3    = 2*NumLevels + 1
Prof_LastO3     = 3*NumLevels
! Cloud Liquid Water (clw) profile
Prof_FirstCLW   = 3*NumLevels + 1
Prof_LastCLW    = 4*NumLevels

! ******** (Add additional profiles here!) **********
! Cloud liquid water
Prof_LWP        = 4*NumLevels + 1
! 1.5m temperature and humidity 
Prof_T2         = 4*NumLevels + 2
Prof_q2         = 4*NumLevels + 3
! Surface pressure
Prof_pstar      = 4*NumLevels + 4
! 10m wind
Prof_uwind      = 4*NumLevels + 5
Prof_vwind      = 4*NumLevels + 6
! Skin temperature
Prof_Tstar      = 4*NumLevels + 7
! Cloud top pressure and effective cloud fraction
Prof_CTP        = 4*NumLevels + 8
Prof_CloudCover = 4*NumLevels + 9
! Total size of vector
ProfSize = 4*NumLevels + 9


END SUBROUTINE NWPSAF_Read_Background
