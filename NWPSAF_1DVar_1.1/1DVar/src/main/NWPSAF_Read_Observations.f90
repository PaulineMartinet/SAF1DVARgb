SUBROUTINE NWPSAF_Read_Observations ( &
    ObNumber, &      ! in
    BackGrModelOb, & ! inout
    Observations, &  ! inout
    RT_Params)       ! inout
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
! Description: Read in Observations for 1DVar
!
! Outputs:
!   Observations - The Observation structure field for use in the 
!                  NWPSAF 1DVar
!
! Method:
!    Provides observation information to the standalone NWPSAF 1DVar
!    program by initialising all required variables in the Observations
!    structure and reading in the observations from the file
!    NWPSAF_Observations.dat
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.1      15/07/99  Original Version
! 2.2      30/04/02  Removed the need for explicit formats in input files.
!                    (Calls new subroutine NWPSAF_ReadHeaders).  ADC.
! 2.3      22/05/02  Changes to allow more flexibility on setting up 
!                    multi-instruements such as ATOVS.         
!                    Remove status part of Ob structure.       ADC.
! 3.0.4    04/03/04  Set RT_Params % RTCoeffs to 0.            A. Collard.
! 3.0.5    29/03/04  Remove all references to orbit height and replace SatView
!                    with SatZenith.                           A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
! 25        20/02/12 Added step to read in date to ob structure. P.Weston.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Constants, ONLY : &
     MissData_R, &
     MissData_I

USE NWPSAFMod_CovarianceMatrices, ONLY : &
     Default_SatID_Text

USE NWPSAFMod_ObsInfo, ONLY : &
     Coord_type,            &
     ElementHeader_type,    &
     ModelOb_type,          &
     OB_type

USE NWPSAFMod_Params, ONLY :  &
     StatusFatal,             &
     UsePCs,                  &
     NPCScores,               &
     CalcRadiance,            &
     NumChannels => MaxChanUsed

USE NWPSAFMod_RTModel, ONLY : &
     RTParams_Type

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_OpenFile.interface'
INCLUDE 'NWPSAF_ReadHeaders.interface'

! Subroutine Arguments 

INTEGER, INTENT(IN) :: ObNumber               ! Observation Number to Read In
TYPE(ModelOB_type), INTENT(INOUT) :: BackGrModelOB !BackGround at Ob
TYPE(OB_type), INTENT(INOUT) :: Observations     ! Observation structure
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params  ! RTModel Parameters Structure

! Declarations

INTEGER :: IComp = 0
INTEGER :: I, II, Instr
INTEGER :: First_Channel, Last_Channel

INTEGER, SAVE :: fileunit      ! I/O unit number
INTEGER, SAVE :: NumObs        ! # of soundings in Observation structure
INTEGER, SAVE :: NumObRead = 0 ! # of observations read in
INTEGER :: NumHeaderLines = 10 ! # of lines of informational header in 
                               ! ObsFileName 
INTEGER :: NumHeaderLines2 = 3 ! # of lines of 2nd header in ObsFileName 

INTEGER :: ReadStatus          ! IOSTAT keyword value
INTEGER :: WMO, WMO_Old

LOGICAL, SAVE :: FirstCall = .TRUE.

CHARACTER(LEN=50) :: ErrorMessage(2)  ! Message for NWPSAF_Report
CHARACTER(LEN=*),  PARAMETER :: RoutineName = 'NWPSAF_Read_Observations'
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER(LEN=80) :: Header                    ! Header for ObsFileName
CHARACTER(LEN=80) :: ObsFileName = "ObsFile.dat"
CHARACTER(LEN=8)  :: ObsNumChar     ! Character Representation of 
                                    ! Observation Number
CHARACTER(LEN=3)  :: Status = "OLD"
CHARACTER(LEN=80) :: String = " "              ! Header for ObsFileName

TYPE(ElementHeader_type) :: Dummy_header  ! holds initialization data
TYPE(Coord_type)         :: Dummy_Coord   ! holds initialization data

REAL :: Vars(3)
CHARACTER(LEN=10)  :: VarChar
!-----------------------------------------------------------------------------

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

ErrorMessage(:) = ''

Dummy_header% NumLev = 0

Dummy_Coord% Value    = MissData_R

Vars(:) = 0.0

!-----------------------------------
!1. Initialize observation structure
!-----------------------------------

 
!--------------------------------------------------------
!1.1 Open the Observations File and Read in Header Data
!--------------------------------------------------------

IF (FirstCall) THEN

  CALL NWPSAF_OpenFile( Trim(ObsFileName),   & ! in
    Access,              & ! in
    Action,              & ! in
    Status,              & ! in
    Form,                & ! inout
    fileunit )             ! out

  ! Allow NumHeaderLines Lines for an informational header

  ReadHeader : DO i = 1, NumHeaderLines
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
    IF (ReadStatus /= 0) THEN
        Errormessage(1) = ' Error Reading Observations File Header '
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO ReadHeader
  
  ! Read in number of observations
  
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  NumObs = NINT(Vars(1))
  IF (ReadStatus /= 0 .OR. NumObs <= 0) THEN
    Errormessage(1) = ' Error Reading Obs File: No. of Obs '
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  Observations% header% NumObsLocal = NumObs
  
  ! Read in number of channels per observation
  ! Note: Number of channels is assumed to be the same for each observation
  
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  NumChannels = NINT(Vars(1))
  IF (ReadStatus /= 0 .OR. NumChannels <= 0) THEN
    Errormessage(1) = ' Error Reading Obs File: No. of Channels '
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
  
  ! Read in number of separate instruments that make up the obs type
  ! Note: This is assumed to be the same for each observation
  
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
  RT_Params % Num_Instruments = NINT(Vars(1))
  IF (ReadStatus /= 0 .OR. RT_Params % Num_Instruments <= 0) THEN
    Errormessage(1) = ' Error Reading Obs File: No. of Instruments '
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  ! ALLOCATE the arrays to store the files containing RT coefficients,
  ! and sub-instrument specifications.

  ALLOCATE(RT_Params % RTCoeffs(RT_Params % Num_Instruments))
  ! Set to zero for now to let RT model's hard-coded path be used.
  RT_Params % RTCoeffs(:) = 0
  ALLOCATE(RT_Params % SeriesChoice(RT_Params % Num_Instruments)) 
  ALLOCATE(RT_Params % PlatformChoice(RT_Params % Num_Instruments))
  ALLOCATE(RT_Params % SubTypeChoice(RT_Params % Num_Instruments))

  ! Three more header lines
  
  ReadHeader2 : DO i = 1, NumHeaderLines2
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
    IF (ReadStatus /= 0) THEN
      Errormessage(1) = ' Error Reading Observations File Header '
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO ReadHeader2

  ! Now read in instrument definitions

  ! Start by looking to see if there is any information given on the
  ! number of composite instruments (this would be in the last record
  ! read in).
  
  String = ADJUSTL(Header)
  IF (String(1:8) == 'Composit') THEN
    IF (ReadStatus == 0) & 
      CALL NWPSAF_ReadHeaders (1,String,Vars(1:1),ReadStatus)
    IF (ReadStatus == 0) & 
      RT_Params % Num_SatIDs = NINT(Vars(1))
    IF (ReadStatus /= 0 .OR. RT_Params % Num_SatIDs <= 0) THEN
      Errormessage(1) = ' Error Reading Number of Composite Instruments '
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
    ALLOCATE(RT_Params % SatID(RT_Params % Num_SatIDs))
    ! Now read in the names of the composite instruments
    DO I = 1, RT_Params % Num_SatIDs
      READ(UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
      IF (ReadStatus == 0) THEN
        RT_Params % SatID(I) % SatID_Text = TRIM(ADJUSTL(String))
      ELSE
        Errormessage(1) = ' Error Reading Composite Instrument Names '
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          ErrorMessage,           & ! in
          ErrorStatus=StatusFatal ) ! in 
      ENDIF
    END DO
    ! This is a separator
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  ELSE
    RT_Params % Num_SatIDs = 1
    ALLOCATE(RT_Params % SatID(RT_Params % Num_SatIDs))
    RT_Params % SatID(1) % SatID_Text = Default_SatID_Text
  END IF

  IF (String(1:5) == 'Units' ) THEN
    IF (ReadStatus == 0) THEN
      II=Index(String,':')
      IF (II == 0) Then
        ReadStatus = -1
      ELSE
        VarChar=TRIM(ADJUSTL(String(II+1:)))
      END IF
    END IF
    IF (ReadStatus == 0) THEN
      IF ( VarChar(1:8) == 'PC Score' ) UsePCs=.true.
      IF ( VarChar(1:8) == 'Radiance' ) CalcRadiance=.true.
    ELSE
        Errormessage(1) = ' Error Reading Observation Units '
        CALL NWPSAF_Report(         &
          RoutineName,              & ! in
          ErrorMessage,             & ! in
          ErrorStatus=StatusFatal ) ! in 
    END IF
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  END IF

  ! Now read in data on composite instruments 
  ALLOCATE(RT_Params % NumChans(RT_Params % Num_Instruments))
  ALLOCATE(RT_Params % First_Channel_for_Instrument(RT_Params % Num_Instruments))
  ALLOCATE(RT_Params % Instrument_Number(RT_Params % Num_SatIDs, NumChannels))
  ALLOCATE(RT_Params % Absolute_Channel_Number(NumChannels))
  RT_Params % NumChans(:) = 0
  RT_Params % First_Channel_for_Instrument(:) = 0
  RT_Params % Instrument_Number(:,:) = 0
  RT_Params % Absolute_Channel_Number(:) = 0  ! This will result in the
                                              ! default RTTOV option of
                                              ! all channels being used.

  IComp = 0
  WMO_Old = 0
  DO I = 1, RT_Params % Num_Instruments
    READ (UNIT=fileunit, FMT=* , IOSTAT=ReadStatus) &
          RT_Params % SeriesChoice(I),   &
          RT_Params % PlatformChoice(I), &
          RT_Params % SubTypeChoice(I),  &
          First_Channel,                 &
          Last_Channel,                  &
          WMO

    IF (ReadStatus /= 0 .OR. RT_Params % SeriesChoice(I) < 0 .OR. &
        RT_Params % PlatformChoice(I) < 0 .OR. &
        RT_Params % SubTypeChoice(I) < 0 .OR. &
        First_Channel <= 0 .OR. Last_Channel > NumChannels .OR. &
        Last_Channel < First_Channel) THEN
      Errormessage(1) = ' Error Reading Instrument Associations '
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF

    IF (WMO /= WMO_Old) THEN
      IComp = IComp + 1
      IF (IComp > RT_Params % Num_SatIDs) THEN
        Errormessage(1) = 'Trying to read too many composite instruments'
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          ErrorMessage,           & ! in
          ErrorStatus=StatusFatal ) ! in 
      END IF
      RT_Params % SatID(IComp) % WMO = WMO
      RT_Params % SatID(IComp) % First_Instr = I
      RT_Params % SatID(IComp) % Num_Channels = NumChannels
      WMO_Old = WMO
    END IF
    IF (IComp > 0) RT_Params % SatID(IComp) % Last_Instr = I

    RT_Params % First_Channel_for_Instrument(I) = First_Channel
    RT_Params % Instrument_Number( IComp,First_Channel:Last_Channel) = I
    RT_Params % NumChans(I) = Last_Channel - First_Channel + 1

  END DO

  ! Check that all channels are assigned an instrument
  
  DO Instr = 1, RT_Params % Num_SatIDs
    IF (.NOT.(ALL(RT_Params % Instrument_Number(Instr,:) == 0))) THEN
      IF ( ANY((RT_Params % Instrument_Number(Instr,:) <= 0)  .OR. &
           ( RT_Params % Instrument_Number(Instr,:) >              &
             RT_Params % Num_Instruments)                  ) ) THEN
        Errormessage(1) = 'Not all channels have been assigned an instrument '
        WRITE(Errormessage(2),FMT='(A,I3)') 'for instrument ',Instr
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          ErrorMessage,           & ! in
          ErrorStatus=StatusFatal ) ! in 
      END IF
    END IF
  END DO

  ! Check to see if we have any information on the true instrument
  ! channel numbers.  This is for cases where we don't want to have
  ! all the instrument's channels in the observation file.
  ! This choice applies to all compoiste instruments.
  ! * Please note * This information is only used in the initialisation
  ! of RTTOV and at all other times the channel numbers should be treated
  ! as consecutive integers (this includes ChannelChoice and R-Matrix).
  ! 
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  String = ADJUSTL(String)
  IF (String(1:8) == 'Channels' .or. String(1:3) == 'PCs' ) THEN
    IF (ReadStatus == 0) READ(UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) &
      RT_Params % Absolute_Channel_Number(:) 
    ! This line reads a seperator
    IF (ReadStatus == 0) &
      READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
    IF (ReadStatus /= 0 .OR. &
        ANY(RT_Params % Absolute_Channel_Number(:) <= 0)) THEN
        Errormessage(1) = ' Error Reading Absolute Instrument Channels '
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END IF


  ! Allocate relevent parts of the Observations structure and initialise 
  ! associated header elements.

  Observations% header% BriTemp         = Dummy_header
  Observations% header% Radiance        = Dummy_header
  Observations% header% PCScore         = Dummy_header
  Nullify(Observations % BriTemp)
  Nullify(Observations % Radiance)
  Nullify(Observations % PCScore)

  IF (UsePCs) THEN
    Observations% header% PCScore % NumLev = NumChannels
    ALLOCATE( Observations% PCScore( NumChannels ))
    NPCScores=NumChannels
  ELSE IF (CalcRadiance) THEN
    Observations% header% Radiance% NumLev = NumChannels
    ALLOCATE( Observations% Radiance( NumChannels ))
  ELSE
    Observations% header% BriTemp% NumLev = NumChannels
    ALLOCATE( Observations% BriTemp( NumChannels ))
  END IF

  ! Also allocate the BackGrModelOb brightness temperaure 

  ALLOCATE( BackGrModelOb% BriTemp( NumChannels ))

  FirstCall = .FALSE.

END IF

IF (UsePCs) THEN
  Observations% PCScore(:)  = MissData_R
ELSE IF (CalcRadiance) THEN
  Observations% Radiance(:) = MissData_R
ELSE
  Observations% BriTemp(:)  = MissData_R
END IF
Observations% Date(:)     = MissData_I
Observations% Elevation   = MissData_R
Observations% ID          = MissData_I
Observations% Latitude    = Dummy_Coord
Observations% Longitude   = Dummy_Coord
Observations% ObsType     = MissData_I
Observations% SatId       = MissData_I
Observations% SatZenith   = MissData_R
Observations% SolarZenith = MissData_R
Observations% Surface     = MissData_I

Observations% header% Date         = Dummy_Header
Observations% header% Elevation    = Dummy_Header
Observations% header% Latitude     = Dummy_Header
Observations% header% Longitude    = Dummy_Header
Observations% header% ObsType      = Dummy_Header
Observations% header% SatId        = Dummy_Header
Observations% header% SatZenith    = Dummy_Header
Observations% header% SolarZenith  = Dummy_Header
Observations% header% Surface      = Dummy_Header

!------------------------------------------------------------
!1.2 Loop through observations and read in elements
!------------------------------------------------------------

ReadData: DO WHILE (NumObRead < ObNumber)

  NumObRead = NumObRead + 1

  ! Observation identifiers:
  !-------
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observations% ID      = NINT(Vars(1))
  Observations% ObsType = NINT(Vars(2))
  Observations% SatID   = NINT(Vars(3))
  IF (ReadStatus /= 0) THEN
    WRITE (ObsNumChar, FMT='(I8)') ObNumber
    Errormessage(1) = ' Error Reading Observation Header #'//Trim(ObsNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
  IF (Observations%ObsType <= 0) Observations % ObsType = MissData_I
  IF (Observations%SatID <= 0) Observations % SatID = MissData_I

  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (String(1:4) == 'Year') THEN
    ! Temporal information (required if reading in emissivity atlas)
    !-------
    IF (ReadStatus == 0) &
      CALL NWPSAF_ReadHeaders (3,String,Vars(1:3),ReadStatus)
    Observations% Date(1) = NINT(Vars(1))
    Observations% Date(2) = NINT(Vars(2))
    Observations% Date(3) = NINT(Vars(3))
    IF (ReadStatus /= 0) THEN
      WRITE (ObsNumChar, FMT='(I8)') ObNumber
      Errormessage(1) = ' Error Reading Observation Header #'// &
        Trim(ObsNumChar)
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
    IF (Observations%Date(1) <= 0) Observations % Date(1) = MissData_I
    IF (Observations%Date(2) <= 0 .OR. Observations%Date(2) > 12) &
      Observations % Date(2) = MissData_I
    IF (Observations%Date(3) <= 0 .OR. Observations%Date(3) > 31) &
      Observations % Date(3) = MissData_I
    READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  END IF

  ! Navigation information:
  !-------
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observations% Latitude% Value  = Vars(1)
  Observations% Longitude% Value = Vars(2)
  Observations% Elevation        = Vars(3)
  IF (ReadStatus /= 0) THEN
    WRITE (ObsNumChar, FMT='(I8)') ObNumber
    Errormessage(1) = ' Error Reading Observation Header #'//Trim(ObsNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
  IF (ABS(Observations % Latitude % Value) > 90.) &
    Observations % Latitude % Value = MissData_R
  IF (ABS(Observations% Longitude% Value) > 180.) &
    Observations % Longitude % Value = MissData_R
  IF (Observations% Elevation < -1000.) Observations % Elevation = MissData_R

  ! More navigation information:
  !-------
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus == 0) &
    CALL NWPSAF_ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observations% Surface     = NINT(Vars(1))
  Observations% SatZenith   = Vars(2) 
  Observations% SolarZenith = Vars(3)
  IF (ReadStatus /= 0) THEN
    WRITE (ObsNumChar, FMT='(I8)') I 
    Errormessage(1) = ' Error Reading Observation Header #'//Trim(ObsNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
  IF (Observations % Surface < 0) Observations % Surface = MissData_I
  IF (Observations % SatZenith > 90. .OR. Observations % SatZenith < 0.) &
    Observations % SatZenith = MissData_R
  IF (Observations%SolarZenith > 180. .OR. & 
    Observations%SolarZenith < 0.) Observations % Surface = MissData_R
  
  ! Header message:
  !-------
  READ (UNIT=fileunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (ReadStatus /= 0) THEN
    Errormessage(1) = ' Error Reading Observations File Header '
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  ! Measurement Itself: 
  !-------

  IF (UsePCs) THEN
    READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) Observations% PCScore(:)
  ELSE IF (CalcRadiance) THEN
    READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) Observations% Radiance(:)
  ELSE
    READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) Observations% BriTemp(:)
  END IF
  IF (ReadStatus /= 0) THEN
    WRITE (ObsNumChar, FMT='(I8)') ObNumber
    Errormessage(1) = ' Error Reading Observation #'//Trim(ObsNumChar)
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

END DO ReadData

! 
!  If the last observation has been read in, close the file. 
!

IF (ObNumber == NumObs) CLOSE(FileUnit)

END SUBROUTINE NWPSAF_Read_Observations
