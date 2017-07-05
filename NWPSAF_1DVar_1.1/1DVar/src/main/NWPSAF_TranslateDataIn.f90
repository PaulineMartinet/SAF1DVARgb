!+ Transform input data into a form that can be processed by RTIASI.  
!  IASI Only.

Subroutine NWPSAF_TranslateDataIn( &
     Obs,            & ! in
     RT_Params,      & ! inout
     Observation,     & ! out
     GoodChannels,   & ! out
     ErrorCode)        ! out

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
! Description: Transform the input data into a form that can be better 
!              processed by the Fastmodel_Interface and the minimisation 
!              routines.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     02/01/01 Original code based on ATOVS_TranslateDataIn.
!                                                             A.D. Collard
!                                                             Met. Office
! 3.0.5   29/03/04 Remove references to SatView and dispense with 
!                  orbit heights.
!
! Ticket  Date     Comment
! ------- -------- -------
! 13      10/02/09 Added latitude and elevation for RTTOV9. E. Pavelin.
! 25      20/02/12 Added longitude for RTTOV10.             P. Weston.

!
! Bugs:
!
!   Bias correction, surface emissivity and some aspects of quality 
!   control are not currently implemented in this code.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type

USE NWPSAFMod_Channellist, ONLY :  &
    BackChans

USE NWPSAFMod_Constants, ONLY : &
    MaxTemperature, &
    MinTemperature, & 
    MinRadiance, &
    MaxRadiance

USE NWPSAFMod_Params, ONLY : &
    GeneralMode, &
    VerboseMode, &
    QC_BadRawBT, &
    UsePCs, &
    CalcRadiance

USE NWPSAFMod_RTmodel, ONLY : &
    RTParams_Type, &
    Num_Profs, &
    NumPredChansPC

IMPLICIT NONE

! Subroutine arguments:
TYPE(Ob_type),       INTENT(INOUT) :: Obs            ! Observed/Retrieval data 
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params      ! RT Model data
INTEGER,             INTENT(OUT)   :: ErrorCode
LOGICAL,             INTENT(OUT)   :: GoodChannels(:)
REAL,                INTENT(OUT)   :: Observation(:)  ! BT/Rad observations

! Local constants:
REAL, PARAMETER    :: EarthRadius = 6371229.  ! Radius of the Earth 
                                              ! (in meters).
REAL, PARAMETER    :: Pi = 3.14159265358979323846  ! Pi  

REAL :: MaxToCheck
REAL :: MinToCheck

! Local variables:
INTEGER :: i
INTEGER :: channel                       ! Loop counter
REAL, POINTER :: RawObs(:)               ! Raw (uncorrected) BT's/Radiances
! LOGICAL :: ValidData
LOGICAL :: DoBias

!-----------------------------------------------------------------------------


!1) Initialise
!--

!Data flag
ErrorCode = 0
DoBias = .TRUE.

! Temporal information
RT_Params % Date              = Obs% Date

! Angles, latitude and elevation
RT_Params % SatZenithAngle    = Obs% SatZenith
RT_Params % SolarZenAngle     = Obs% SolarZenith
! These currently not provided in ObsFile:
!RT_Params % SatAzimAngle      = Obs% SatAzimth
!RT_Params % SatSolarAzimAngle = Obs% SolarAzimth
RT_Params % SatAzimAngle      = 0.0
RT_Params % SatSolarAzimAngle = 0.0
RT_Params % Latitude          = Obs% Latitude % Value
RT_Params % Longitude         = Obs% Longitude % Value
RT_Params % Elevation         = Obs% Elevation


!Observations
IF ( UsePCs ) THEN
  RawObs => Obs % PCScore(:)
ELSE IF ( CalcRadiance ) THEN
  RawObs => Obs % Radiance(:)
  MaxToCheck=MaxRadiance
  MinToCheck=MinRadiance
ELSE
  RawObs => Obs % BriTemp(:)
  MaxToCheck=MaxTemperature
  MinToCheck=MinTemperature
END IF

!2) Quality control raw brightness temperature data (use gross limit check)
!--
IF (.not. UsePCs) THEN !not able to do this check for PC scores
  GoodChannels(:) = .FALSE.
  DO i = 1, BackChans% numchans
    channel = BackChans% channels(i)
    IF ( RawObs(channel) <= MaxToCheck .AND. &
          RawObs(channel) >= MinToCheck      ) THEN
        GoodChannels(channel) = .TRUE.
    END IF
  END DO

  IF (ANY(.NOT.GoodChannels(BackChans% channels(1:BackChans% numchans)))) THEN
    ErrorCode = QC_BadRawBT
    IF ( GeneralMode >= VerboseMode ) THEN
        WRITE(*,*) 'INVALID DATA:'
        DO i = 1, BackChans% numchans
          channel = BackChans% channels(i)
          IF ( .NOT. GoodChannels(channel) ) &
                WRITE(*,*) 'Bad BT for channel ',channel,' = ',RawObs(channel) 
                WRITE(*,*) 'Limits: ', MinToCheck, MaxToCheck
        END DO
        WRITE(*,*)
    END IF
  END IF
END IF

!3) Perform bias correction
!--

IF ( DoBias ) THEN
  !! ********* No Bias Correction for now **************
  Observation(:) = RawObs(:)
ELSE
  Observation(:) = RawObs(:)
END IF


!4) Emissivities
!--
IF (UsePCs) THEN
  ALLOCATE(RT_Params % RTEmissivity(NumPredChansPC * Num_Profs))
ELSE
  ALLOCATE(RT_Params % RTEmissivity(Backchans % numchans * Num_Profs))
END IF
! Not sure if this makes sense for addradrec though...


END SUBROUTINE NWPSAF_TranslateDataIn
