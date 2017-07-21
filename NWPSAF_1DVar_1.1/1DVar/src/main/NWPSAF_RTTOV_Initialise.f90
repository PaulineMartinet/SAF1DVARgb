Subroutine NWPSAF_RTTOV_Initialise (&
  RT_Params,   & !inout
  RT_Coefs,    & !out
  RT_Opts,     & !out
  Soft_Limits  ) !out
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
! Description: Interface between IASI 1DVar code and RTTOV11 for initialisation
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
!         27/09/16 New. Fiona Smith
!
! Code Description:
!   Language:   Fortran 90
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     GeneralMode,        &
     DebugMode,          &
     Read_CLW_Background,&
     UsePCs,             &
     CalcRadiance,       &
     ipcreg,             &
     Ozone_Present,      &
     Cloud_Min_Pressure, &
     Legacy_Setting,     &
    !PM
     retrieval_in_log,   &
     MwClwRetrieval
     !PM

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type,      &
     Limit_Type,         &
     ProfSize,           &
     Prof_FirstT,        &
     Prof_FirstQ,        &
     Prof_FirstO3,       &
     Prof_LastT,         &
     Prof_LastQ,         &
     Prof_LastO3,        &
     Prof_T2,            &
     Prof_q2,            &
     Prof_Tstar,         &
     Prof_pstar,         &
     Prof_uwind,         &
     Prof_vwind,         &
     Prof_CTP,           &
     Prof_CloudCover,    &
     Num_RTlevels,       &
     UseModelLevels,     &
     maxchans,           &
     maxsensors,         &
     NumPredChansPC

USE NWPSAFMod_Constants, ONLY: &
     q_mixratio_to_ppmv
   
USE rttov_const, ONLY :   &
     gas_id_watervapour,  &
     gas_id_ozone

USE rttov_types, ONLY : &
    rttov_coefs,        &
    rttov_options

USE parkind1, ONLY : &
    jpim, &
    jprb

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'rttov_read_coefs.interface'
INCLUDE 'rttov_user_options_checkinput.interface'

! Subroutine arguments:
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params
TYPE(rttov_coefs), POINTER              :: RT_coefs(:)
TYPE(rttov_options), POINTER            :: RT_opts(:)
TYPE(Limit_Type), INTENT(OUT)           :: Soft_Limits

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV_Initialise"

! Local variables:
INTEGER(kind=jpim) :: ErrorStatus     ! RTTOV Error code
INTEGER            :: I
INTEGER            :: Actual_Num_Instr
INTEGER, POINTER   :: First_Instr_Chan
INTEGER, POINTER   :: Num_Instr_Chan
INTEGER            :: Last_Instr_Chan
INTEGER(kind=jpim) :: Valid_Channels(maxchans,maxsensors)
INTEGER(kind=jpim) :: RTInstrument(3,RT_Params%Num_Instruments)
INTEGER            :: wv_pos, o3_pos

CHARACTER(LEN=80)  :: ErrorMessage(4)  ! Message for NWPSAF_Report
LOGICAL            :: Verbosity


!------------------------------------------------------------------------------

!----- Translate verbosity mode into RTTOV Verbosity flag
Verbosity=.False.
IF ( GeneralMode >= DebugMode ) THEN
  ErrorMessage(1) = 'Initialising RTTOV'
  CALL NWPSAF_Report( RoutineName,ErrorMessage(1:1) )
  Verbosity=.True.
END IF

!------------------------------------------------------------------------------
!1) Perform checks and set up channel/instrument lists
!------------------------------------------------------------------------------

!1.1) Basic checks
!----

! Check that there aren't too many instruments
IF (RT_Params % Num_Instruments > maxsensors) THEN
  ErrorMessage(1)='RT_Params % Num_Instruments > maxsensors'
  WRITE( UNIT=ErrorMessage(2),FMT='(A,I5,A,I5)' ) &
    'Currently they are ',RT_Params % Num_Instruments,' and ',maxsensors
  ErrorMessage(3)='Either reduce the number of instruments required in'// &
    ' the observations file '
  ErrorMessage(4)='or increase maxsensors in NWPSAFMod_RTModel.f90'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

!1.2) Set up channel arrays for reading in coefficients
!----

IF (UsePCs) THEN
  !We've already checked we have only one instrument, in ProcessData 
  !We want to read in all the channels as it is not known in advance which
  !are required for the PC calculation
  Valid_Channels(:,:) = 0
  Actual_Num_Instr=1
ELSE
  Valid_Channels(:,:)=0
  ! Note: Channel numbers will be zero if not defined in obsfile.
  ! RTTOV default is to use all channels in this case.
  DO I = 1, RT_Params % Num_Instruments
    First_Instr_Chan => RT_Params % First_Channel_for_Instrument(I)
    Num_Instr_Chan   => RT_Params % NumChans(I)
    IF (Num_Instr_Chan > 0) Actual_Num_Instr = I
    Last_Instr_Chan  =  First_Instr_Chan + Num_Instr_Chan - 1
    Valid_Channels(1:Num_Instr_Chan,I) = RT_Params % &
      Absolute_Channel_Number(First_Instr_Chan:Last_Instr_Chan)
    IF (maxchans > Num_Instr_Chan) &
      Valid_Channels(Num_Instr_Chan+1:maxchans,I) = 0
  END DO
END IF

!1.3) Set up RTTOV11 Instrument definition
!----

! Reset number of instruments to prevent some unnecessary initialisations
! (this isn't the perfect solution yet)
RT_Params % Num_Instruments = Actual_Num_Instr

RTInstrument(1,:) = RT_Params % SeriesChoice(1:Actual_Num_Instr)
RTInstrument(2,:) = RT_Params % PlatformChoice(1:Actual_Num_Instr)
RTInstrument(3,:) = RT_Params % SubTypeChoice(1:Actual_Num_instr)

ALLOCATE(RT_coefs(RT_Params % Num_Instruments))

ALLOCATE(RT_opts(RT_Params % Num_Instruments))

!1.4) Set up RTTOV11_Opts structure
!----

! Initialise options to default values 

!The following options are version dependent
#ifdef _RTTOV12
  RT_opts(:) % rt_all % plane_parallel    = .false.
  RT_opts(:) % rt_ir % ir_sea_emis_model  = 2 
  RT_Opts(:) % rt_ir % ir_scatt_model     = 2
  RT_Opts(:) % rt_ir % vis_scatt_model    = 1
  RT_Opts(:) % rt_ir % dom_nstreams       = 8
  RT_Opts(:) % rt_ir % dom_accuracy       = 0.
  RT_Opts(:) % rt_ir % dom_opdep_threshold = 0.
#endif

RT_opts(:) % config % do_checkinput     = .true.
RT_opts(:) % config % apply_reg_limits  = .true.
RT_opts(:) % config % verbose           = Verbosity

IF (CalcRadiance) THEN
  RT_opts(:) % rt_all % switchrad        = .false.
ELSE
  RT_opts(:) % rt_all % switchrad        = .true.
END IF
RT_opts(:) % rt_all % addrefrac          = .false.
RT_opts(:) % rt_all % use_q2m            = .true. !****** legacy_setting=.false.
!PM
!RT_opts(:) % rt_all % do_lambertian      = .false.
!PM
RT_opts(:) % rt_mw % fastem_version       = 6  
!PM
!RT_opts(:) % rt_mw % supply_foam_fraction = .false.
!PM

IF (Read_CLW_Background) THEN
  RT_opts(:) % rt_mw % clw_data         = .true.
  !PM
ELSEIF (MwClwRetrieval) THEN
  RT_opts(:) % rt_mw % clw_data         = .true.
  !PM
ELSE
  RT_opts(:) % rt_mw % clw_data         = .false.
END IF

RT_opts(:) % rt_ir % addsolar           = .false.
RT_opts(:) % rt_ir % do_nlte_correction = .false.
RT_opts(:) % rt_ir % addaerosl          = .false.
RT_opts(:) % rt_ir % addclouds          = .false.
RT_opts(:) % rt_ir % user_aer_opt_param = .false.
RT_opts(:) % rt_ir % user_cld_opt_param = .false.
RT_opts(:) % rt_ir % cldstr_threshold   = 0.001_jprb
RT_opts(:) % rt_ir % ozone_data         = .false.
RT_opts(:) % rt_ir % co2_data           = .false.
RT_opts(:) % rt_ir % n2o_data           = .false.
RT_opts(:) % rt_ir % co_data            = .false.
RT_opts(:) % rt_ir % ch4_data           = .false.

RT_opts(:) % interpolation % addinterp  = .false.
RT_opts(:) % interpolation % interp_mode = 5          !****** legacy_setting=1
RT_opts(:) % interpolation % reg_limit_extrap = .true.
RT_opts(:) % interpolation % lgradp     = .false.
RT_opts(:) % interpolation % spacetop   = .true.

!Set up PC forward modelling
RT_opts(:) % rt_ir % pc % addpc         = .false.
RT_opts(:) % rt_ir % pc % addradrec     = .false.
RT_opts(:) % rt_ir % pc % ipcbnd        = 1_jpim
RT_opts(:) % rt_ir % pc % ipcreg        = 1_jpim

!Remember we are assuming only one instrument for PC assimilation
IF ( UsePCs ) THEN
  RT_opts(1) % rt_ir % pc % addpc  = .true.
  RT_opts(1) % rt_ir % pc % ipcreg = ipcreg
  RT_opts(1) % rt_ir % pc % ipcbnd = 1_jpim
  RT_opts(1) % rt_all % addrefrac  = .true.
  IF (.not. RT_opts(1) % rt_ir % pc % addradrec) THEN
    !Just to do this to avoid an error being generated in
    !rttov_user_options_checkinput
    RT_opts(:) % rt_all % switchrad         = .false.
  ENDIF
END IF

!Legacy settings that have a large effect on the brightness temperatures
!Other things you might want to look at are ozone_data, fastem_version and 
!gas_units. The latter is controlled via namelist also
IF (legacy_setting) THEN
  RT_opts(:) % rt_all % use_q2m            = .false.
  RT_opts(:) % interpolation % interp_mode = 1 
END IF
!------------------------------------------------------------------------------
!2) Call RTTOV_SETUP to read in required coefficient files for each instrument
!------------------------------------------------------------------------------

!2.1) Read in coefficients
!----

DO I=1,RT_Params%Num_Instruments

  IF (ANY(Valid_Channels(:,I) > 0) .and. .not. UsePCs ) THEN
    CALL RTTOV_READ_COEFS(                                        &
          ErrorStatus,                                            & ! out
          RT_coefs(I),                                            & ! out
          RT_opts(I),                                             & ! in
          channels = Valid_Channels(1:RT_Params % NumChans(I),I), & ! in
          instrument = RTInstrument(:,I)                          ) ! in
  ELSE
    CALL RTTOV_READ_COEFS(                &
          ErrorStatus,                    & ! out
          RT_coefs(I),                    & ! out
          RT_opts(I),                     & ! in
          instrument = RTInstrument(:,I)  ) ! in
  END IF

  IF (ErrorStatus /= 0) THEN
    ErrorMessage(1)='Error in RTTOV_READ_COEFS'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
      'Error is Code ',ErrorStatus
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

!2.2) Check requested RT_opts against available options in coefs file
!----

  Call RTTOV_USER_OPTIONS_CHECKINPUT( &
    ErrorStatus, & ! out
    RT_opts(I),  & ! in
    RT_coefs(I)  ) ! in

  IF (ErrorStatus /= 0) THEN 
    ErrorMessage(1)='RTTOV Incompatibility between options and coefs file'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorStatus
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

END DO

!2.3) Adjust ozone settings depending on coefficients and background profiles
!----

!Note: this only checks the coefs of the first instrument.
IF ( Ozone_Present .AND. &
     ( ANY(RT_coefs(1) % coef % fmv_gas_id(:) == gas_id_ozone))) THEN
  RT_opts(:) % rt_ir % ozone_data = .true.
END IF

!2.4) Get number of predictor channels for PCRTTOV. 
!----

IF (UsePCs) THEN
  NumPredChansPC = SIZE( RT_coefs(1) % coef_pccomp % pcreg( &
                           RT_opts(1) % rt_ir  % pc % ipcbnd,  &
                           RT_opts(1) % rt_ir  % pc % ipcreg ) % predictindex )
END IF

!-----------------------------------------------------------------------------
!3) Set up limits to use in CheckIteration so profile lies within RTTOV limits
!-----------------------------------------------------------------------------

!NOTE: There are no limits set for CLW

! Determine whether pressure level interpolation is needed.
! Note this may not work well if different coeff levels are used
! for different instruments..
IF ( ANY(RT_coefs(:)%coef%nlevels /= Num_RTLevels) ) THEN
  UseModelLevels=.TRUE.
  DO i=1,RT_Params%Num_Instruments
     RT_opts(i) % interpolation % addinterp = .TRUE.
     RT_opts(i) % config % apply_reg_limits = .TRUE.
     RT_opts(i) % interpolation % lgradp = .TRUE.
    IF(GeneralMode >= DebugMode) Then
       Write(*,*)'N.B.: Using profile interpolation for instrument ',i
    END IF
  END DO

ELSE

  ! Set up soft limits
  ALLOCATE( Soft_Limits % Minimum(ProfSize))
  ALLOCATE( Soft_Limits % Maximum(ProfSize))

  ! Default values:
  Soft_Limits % Minimum(:) = 0.
  Soft_Limits % Maximum(:) = 1.e10

  wv_pos = RT_coefs(1)%coef % fmv_gas_pos( gas_id_watervapour )
  IF(RT_opts(1) % rt_ir % ozone_data) o3_pos = &
    RT_coefs(1) % coef % fmv_gas_pos( gas_id_ozone )


  !3.1) Profile minima (include allowable tolerances)
  !----

  ! Mininum temperature
  Soft_Limits % Minimum(Prof_FirstT:Prof_LastT) = &
    RT_coefs(1) % coef % lim_prfl_tmin(:)-0.5
  ! Minimum water vapour
  IF (retrieval_in_log) THEN
  Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = &
    LOG(RT_coefs(1) % coef % lim_prfl_gmin(:,wv_pos)/ &
    q_mixratio_to_ppmv)-0.5
  ELSE
  Soft_Limits % Minimum(Prof_FirstQ:Prof_LastQ) = RT_coefs(1) % coef % lim_prfl_gmin(:,wv_pos)/ &
    q_mixratio_to_ppmv 
  ENDIF
  ! Minimum ozone (if present)
  IF(RT_opts(1) % rt_ir % ozone_data) THEN
    Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = &
      RT_coefs(1) % coef % lim_prfl_gmin(:,o3_pos)*0.8
  ELSE
    Soft_Limits % Minimum(Prof_FirstO3:Prof_LastO3) = -9999.99
  END IF
  ! Minimum surface variables
  Soft_Limits % Minimum(Prof_T2) = &
    RT_coefs(1) % coef % lim_prfl_tmin(Num_RTLevels)-0.5
    !PM
    IF (retrieval_in_log) THEN
  Soft_Limits % Minimum(Prof_Q2) = &
    LOG(RT_coefs(1) % coef % lim_prfl_gmin(Num_RTLevels,wv_pos)/ &
    q_mixratio_to_ppmv)-0.5
    ELSE
   Soft_Limits % Minimum(Prof_Q2) = &
    RT_coefs(1) % coef % lim_prfl_gmin(Num_RTLevels,wv_pos)/ &
    q_mixratio_to_ppmv   
    ENDIF
    
    
    !PM
  Soft_Limits % Minimum(Prof_PStar) = 300.
  Soft_Limits % Minimum(Prof_UWind) = -100.
  Soft_Limits % Minimum(Prof_VWind) = -100.
  Soft_Limits % Minimum(Prof_TStar) = &
    RT_coefs(1) % coef % lim_prfl_tmin(Num_RTLevels)-0.5
  Soft_Limits % Minimum(Prof_CTP) = Cloud_Min_Pressure
  Soft_Limits % Minimum(Prof_CloudCover) = 0.

  !3.2) Profile maxima
  !----

  ! Maximum temperature
  Soft_Limits % Maximum(Prof_FirstT:Prof_LastT) = &
    RT_coefs(1) % coef % lim_prfl_tmax(:)+0.5
  ! Maximum water vapour
  IF (retrieval_in_log) THEN
  Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
    LOG(RT_coefs(1) % coef % lim_prfl_gmax(:,wv_pos)/ &
    q_mixratio_to_ppmv)+0.5
  ELSE
  Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
    RT_coefs(1) % coef % lim_prfl_gmax(:,wv_pos)/ &
    q_mixratio_to_ppmv
  ENDIF
  ! Maximum ozone (if present)
  IF (RT_opts(1) % rt_ir % ozone_data) Then
    Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
      RT_coefs(1) % coef % lim_prfl_gmax(:,o3_pos)*1.2
  ELSE
    Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = 99999.99
  END IF
  ! Maximum surface variables
  Soft_Limits % Maximum(Prof_T2) = &
    RT_coefs(1) % coef % lim_prfl_tmax(Num_RTLevels)+0.5
 IF (retrieval_in_log) THEN
  Soft_Limits % Maximum(Prof_Q2) = &
    LOG(RT_coefs(1) % coef % lim_prfl_gmax(Num_RTLevels,wv_pos)/ &
    q_mixratio_to_ppmv)+0.5
 ELSE
   Soft_Limits % Maximum(Prof_Q2) = &
    RT_coefs(1) % coef % lim_prfl_gmax(Num_RTLevels,wv_pos)/ &
    q_mixratio_to_ppmv
 ENDIF
  ! (This is the RTTOV hard limit:)
  Soft_Limits % Maximum(Prof_PStar) = 1200.
  Soft_Limits % Maximum(Prof_UWind) = 100.
  Soft_Limits % Maximum(Prof_VWind) = 100.
  Soft_Limits % Maximum(Prof_TStar) = &
    RT_coefs(1) % coef % lim_prfl_tmax(Num_RTLevels)+0.5
  Soft_Limits % Maximum(Prof_CTP) = 1200.
  Soft_Limits % Maximum(Prof_CloudCover) = 1.0

  ! These are the RTTOV hard limits:
  !PM
  IF (retrieval_in_log) THEN
  Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
    MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),-2.99578)
  ELSE
  Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ) = &
    MIN(Soft_Limits % Maximum(Prof_FirstQ:Prof_LastQ),exp(-2.99578)/q_mixratio_to_ppmv)
  ENDIF
!PM
  Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3) = &
    MIN(Soft_Limits % Maximum(Prof_FirstO3:Prof_LastO3),20.)

    IF (retrieval_in_log) THEN
  Soft_Limits % Maximum(Prof_Q2) = &
    MIN(Soft_Limits % Maximum(Prof_Q2),-2.99578)
    ELSE
   Soft_Limits % Maximum(Prof_Q2) = &
    MIN(Soft_Limits % Maximum(Prof_Q2),exp(-2.99578)/q_mixratio_to_ppmv)  
    ENDIF

  !3.3) Check for profile mismatches
  !----

  DO I=1,ProfSize
    IF (Soft_Limits % Minimum(I) >= Soft_Limits % Maximum(I)) THEN
      WRITE( UNIT=ErrorMessage(1),FMT='(A)' ) &
        'Mismatch in RTTOV soft limits'
      WRITE( UNIT=ErrorMessage(2),FMT=* )  I, &
        Soft_Limits % Minimum(I), &
        Soft_Limits % Maximum(I)
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO

END IF ! (UseModelLevels)


END SUBROUTINE NWPSAF_RTTOV_Initialise
