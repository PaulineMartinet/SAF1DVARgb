Subroutine NWPSAF_RTTOV12_Interface ( &
     Fastmodel_Mode,          & ! in
     RT_Params,               & ! inout
     WhichProf,               & ! in
     UsedChans,               & ! in
     ErrorCode)                 ! out
       
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
! Description: Interface between NWPSAF 1DVar code and RTTOV12
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 31      18/01/13 Original version, based on RTTOV10_interface P. Weston.
! 32      28/01/14 Changes to allow for HDF5 coefficients and channel 
!                  selections to be used. P. Weston
!
! Code Description:
!   Language:		Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_Params, ONLY : &
     StatusFatal,        &
     StatusWarning,      & 
     GeneralMode,        &
     DebugMode,          &
     CloudyRetrieval,    &
     RTsea,              &
     Use_EmisAtlas,      &
     Atlas_Dir,          &
     UsePCs,             &
     NPCScores,          &
     CalcRadiance,       &
     Gas_Units

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type,            &
     FastmodelMode_CleanUp, FastmodelMode_Initialise, &
     FastmodelMode_Forward, FastmodelMode_Gradient,   &
     GuessProf, BackGrProf, &
     Prof_FirstT,  Prof_LastT, &
     Prof_FirstQ,  Prof_LastQ, &
     Prof_FirstO3, Prof_LastO3, &
     Prof_FirstCLW,Prof_LastCLW, &
     Prof_T2, Prof_q2, &
     Prof_Tstar, Prof_pstar, &
     Prof_uwind, Prof_vwind, &
     Prof_CTP, Prof_CloudCover,  &
     Num_RTlevels,       &
     Num_Profs,          &
     Soft_Limits,        &
     UseModelLevels,     &
     RT_coefs,        &
     RT_opts,         &
     RT_Chanprof
 
USE NWPSAFMod_Constants, ONLY : &
    q_mixratio_to_ppmv

USE rttov_const, ONLY :   &
     errorstatus_fatal,   &
     sensor_id_mw,        &
     sensor_id_po

USE rttov_types, ONLY:   &
     rttov_profile,      &
     rttov_radiance,     &
     rttov_transmission, &
     rttov_emissivity,   &
     rttov_pccomp

USE mod_rttov_emis_atlas, ONLY : &
    atlas_type_ir, &
    atlas_type_mw,       &
    rttov_emis_atlas_data

Use parkind1, Only : jpim,jprb

IMPLICIT NONE

#ifdef _CompileEmissAtlas
INCLUDE 'rttov_setup_emis_atlas.interface'
INCLUDE 'rttov_deallocate_emis_atlas.interface'
INCLUDE 'rttov_get_emis.interface'
#endif
INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'rttov_direct.interface'
INCLUDE 'rttov_k.interface'
INCLUDE 'rttov_dealloc_coefs.interface'
INCLUDE 'rttov_alloc_prof.interface'
INCLUDE 'NWPSAF_RTTOV_Initialise.interface'
INCLUDE 'NWPSAF_RTTOV12_Allocate.interface'
INCLUDE 'NWPSAF_RTTOV12_GetHMatrix.interface'
include "rttov_print_opts.interface"
include "rttov_print_profile.interface"
include "rttov_print_info.interface"

! Subroutine arguments:

INTEGER, INTENT(IN)                     :: Fastmodel_Mode  !Forward/Gradient etc
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params   ! Info for RT Model
INTEGER, INTENT(IN)                     :: WhichProf   ! Which profile to use 
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans
INTEGER, INTENT(OUT)                    :: ErrorCode


! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV12_Interface"

! Local variables:

INTEGER :: I
INTEGER :: IChan1, IChan2
INTEGER :: DeallocateError
INTEGER :: First_Instr, Last_Instr
INTEGER :: Instrument
INTEGER :: NumInstChans
INTEGER :: SatIndex
INTEGER :: First_Chan_Pos !) Position of first and last channels for  
INTEGER :: Last_Chan_Pos  !) current instrument in RT_Params % TotalBTs
INTEGER, Allocatable :: Inst_Chans (:)

CHARACTER(LEN=80) :: Message(4)  ! Message for NWPSAF_Report

REAL, POINTER :: RTProf(:)


! Arrays used in RTTOV interface 
!--------

TYPE(rttov_profile)                 :: Profiles(Num_Profs)
TYPE(rttov_profile), ALLOCATABLE    :: Profiles_K(:)
TYPE(rttov_profile)                 :: Profiles_K_PC(NPCScores)
TYPE(rttov_emissivity), ALLOCATABLE :: Surf_Emiss(:) 
TYPE(rttov_emissivity), ALLOCATABLE :: Surf_Emiss_K(:) 
LOGICAL, ALLOCATABLE                :: calcemiss(:)
TYPE(rttov_radiance)                :: Radiance
TYPE(rttov_radiance)                :: Radiance_K
Type(rttov_transmission)            :: Transmission
Type(rttov_transmission)            :: transmission_K
TYPE(rttov_pccomp)                  :: PCcomp
TYPE(rttov_pccomp)                  :: PCcomp_K
TYPE(rttov_emis_atlas_data)         :: emis_atlas               ! Data structure for emissivity atlas

INTEGER, PARAMETER :: ASW_ALLOCATE = 1
INTEGER, PARAMETER :: ASW_DEALLOCATE = 0

!---------------------------------------------------------------

IF ( GeneralMode >= DebugMode ) THEN
  WRITE(UNIT=Message(1),FMT='("RTTOV12 Called in MODE = ",i3)') &
  Fastmodel_Mode
  CALL NWPSAF_Report( RoutineName,Message(1:1) )
END IF

ErrorCode = 0
Message(:)=' '

!------------------------------------------------------------------------------
! ****** SETUP MODE ******
!------------------------------------------------------------------------------
RTTOV_FastmodelMode : IF ( Fastmodel_Mode == FastmodelMode_Initialise ) THEN
  CALL NWPSAF_RTTOV_Initialise( &
    RT_Params,   &
    RT_coefs, &
    RT_opts,  &
    Soft_Limits  )
!------------------------------------------------------------------------------
! ****** CLEAN UP MODE ******
!------------------------------------------------------------------------------
ELSE IF (FastModel_Mode == FastModelMode_CleanUp) THEN RTTOV_FastmodelMode

  DO I=1,RT_Params % Num_Instruments

    Call RTTOV_DEALLOC_COEFS( ErrorCode,    & ! out
                              RT_coefs(I))  ! in
    IF (ErrorCode /= 0) THEN 
      Message(1)='Error in RTTOV_DEALLOC_COEFS'
      WRITE( UNIT=Message(2),FMT='(A,I3)' ) 'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        Message,                & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END DO

  IF (ASSOCIATED(RT_coefs)) Deallocate(RT_coefs)
  IF (ASSOCIATED(RT_opts)) Deallocate(RT_opts)
  IF (.not. UseModelLevels) THEN
     Deallocate( Soft_Limits % Maximum )
     Deallocate( Soft_Limits % Minimum )
  END IF


!------------------------------------------------------------------------------
! ****** CALL FORWARD OR GRADIENT CODE ******
!------------------------------------------------------------------------------
ELSE RTTOV_FastmodelMode

  SatIndex = RT_Params % SatIndex
  First_Instr = RT_Params % SatID(SatIndex) % First_Instr
  Last_Instr = RT_Params % SatID(SatIndex) % Last_Instr

  !---------------------------------------------------------------------
  !  1.) Put profile into the RTModel vector 
  !---------------------------------------------------------------------
  ! Note that the profile structure contains a lot of optional components
  ! only those that requested in the options structure can be set.
  !---------------------------------------------------------------------

  IF (WhichProf == BackGrProf) THEN
    RTProf => RT_Params % RTBack
  ELSE IF (WhichProf == GuessProf) THEN
    RTProf => RT_Params % RTGuess
  ELSE
    WRITE( UNIT=Message(1),FMT='(A,I2,A)' )  &
           'Incorrect profile (',WhichProf,') specified'
    CALL NWPSAF_Report( &
      RoutineName,              & ! in
      Message,                  & ! in
      ErrorStatus=StatusWarning ) ! in 
    Message(:)=''
   END IF   

  !1.1) Allocate profile and nullify/initialise to zero
  !----

  !Note: assumes consistent options required for composite instruments, and only
  !looks at the first. e.g. can't use CLW profile just for MHS as it will look
  !at HIRS first.
  CALL RTTOV_ALLOC_PROF(&
    ErrorCode,          &
    1,                  &
    Profiles(1),        &
    Num_RTLevels,       & 
    RT_Opts(1),         &
    asw = ASW_ALLOCATE, &
    coefs = RT_coefs(1),&
    init = .true. )


  IF (ErrorCode /= 0) THEN 
    Message(1)='Error in RTTOV_ALLOC_PROF'
    WRITE( UNIT=Message(2),FMT='(A,I3)' ) 'Error is Code ',ErrorCode
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      Message,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  !1.2) Set basic information
  !----
  profiles(1) % nlevels = Num_RTLevels
  profiles(1) % nlayers = Num_RTLevels - 1
  profiles(1) % gas_Units = Gas_Units ! ppmv 

  ! 1.3) Fill in the atmospheric profile variables
  !----
  profiles(1) % p(:) = RT_Params%Pressure_Pa(:)/100.0
  profiles(1) % t(:) = RTProf(Prof_FirstT : Prof_LastT)
  profiles(1) % q(:) = EXP(RTProf(Prof_FirstQ:Prof_LastQ))*q_mixratio_to_ppmv
  IF ( ANY(RT_opts(:) % rt_ir % ozone_data) ) THEN
    profiles(1) % o3(:) = RTProf(Prof_FirstO3 : Prof_LastO3)
  END IF
  IF( ANY(RT_opts(:) % rt_mw % clw_data) ) THEN
    profiles(1) % clw(:) = RTProf(Prof_FirstCLW : Prof_LastCLW)
  END IF

  ! 1.4)  Fill in the Surface Variables
  !----
  profiles(1) % s2m % t = RTProf(Prof_T2)
  profiles(1) % s2m % q = EXP(RTProf(Prof_Q2)) * q_mixratio_to_ppmv
  profiles(1) % s2m % p = RTProf(Prof_PStar)
  profiles(1) % s2m % u = RTProf(Prof_UWind)
  profiles(1) % s2m % v = RTProf(Prof_VWind)

  !1.5) Fill in the Skin Variables
  !----
  profiles(1)%skin%t = RTProf(Prof_TStar)
  profiles(1)%skin%fastem(:) = 0.0
  profiles(1)%skin%surftype = RT_Params % RTSurfaceType
  Profiles(1)%skin%watertype = 1 ! Default to ocean water
   
  !1.6) Fill in the Cloud Variables
  !----
  ! This is cloud top pressure (make sure is strictly within the 
  ! atmosphere as otherwise RTTOV may crash!)
  ! Should not be required with opts % apply_reg_limits, but maintained for now
  profiles(1)%ctp = RTProf(Prof_CTP)
  profiles(1)%ctp = MIN(profiles(1)%ctp,RTProf(Prof_PStar)*1._JPRB)
  profiles(1)%ctp = MAX(profiles(1)%ctp, &
                        MINVAL(RT_Params % Pressure_Pa(:) / 100.*1._JPRB))
  ! This is cloud fraction (make sure is strictly within the 
  ! range 0-1 as otherwise RTTOV may crash!)
  ! Should not be required with opts % apply_reg_limits, but maintained for now
  profiles(1)%cfraction = RTProf(Prof_CloudCover)
  profiles(1)%cfraction = MIN(profiles(1)%cfraction,1._JPRB)
  profiles(1)%cfraction = MAX(profiles(1)%cfraction,0._JPRB)

  !1.7) Set geometry etc
  !----
  Profiles(1)%zenangle = RT_Params % SatZenithAngle
  Profiles(1)%azangle = RT_Params % SatAzimAngle
  Profiles(1)%sunzenangle = RT_Params % SolarZenAngle
  Profiles(1)%sunazangle = RT_Params % SatSolarAzimAngle
  Profiles(1)%date = RT_Params % Date
  Profiles(1)%elevation = RT_Params % Elevation / 1000._JPRB
  Profiles(1)%latitude = RT_Params % Latitude
  Profiles(1)%longitude = RT_Params % Longitude

  !1.7) Components of profile currently not in use
  !----
  !do not use arrays unless you have changed the code to ask for that option
  !as they will not have been allocated
  !chemical constituents
  !profiles(1) % co2(:) = 0.0
  !profiles(1) % n2o(:) = 0.0
  !profiles(1) % co(:)  = 0.0
  !profiles(1) % ch4(:) = 0.0
  !profiles(1) % so2(:) = 0.0

  !ice particles
  profiles(1) % ice_scheme = 1 ! Ice scheme 
  profiles(1) % idg = 1 ! Ice effective diameter scheme
  !profiles(1) % icede(:) = 0.0 ! Ice effective diameter scheme

  !aerosols
  profiles(1) % mmr_cldaer = .false. !units of aerosol input
  !profiles(1) % aerosols(:,:) = 0.0

  !clouds
  !profiles(1) % cloud(:,:) = 0.0
  !profiles(1) % cfrac(:) = 0.0

  !surface variables
  profiles(1) % s2m % wfetc = 100000.0 ! Wind fetch
  profiles(1) % s2m % o = 0.0 ! Not even used inside rttov!

  !skin variables
  profiles(1) % skin % salinity = 0 !should we set this?
  profiles(1) % skin % foam_fraction = 0.0
  profiles(1) % skin % snow_fraction = 0.0

  !Magnetic filed info for Zeeman splitting
  profiles(1) % Be = 0.0
  profiles(1) % cosbk = 0.0

  !----------------------------------------------------------------------------
  ! 2) Loop over instruments: do remaining setup and call RTTOV forward/K code
  !----------------------------------------------------------------------------
  Instrument_Loop : DO Instrument = First_Instr, Last_Instr

    !2.1) Initialise arrays which need to be sized depending on input
    !----
    IF ( UsePCs ) THEN
      NumInstChans = &
        SIZE(RT_coefs(First_Instr) % coef_pccomp % pcreg( &
          RT_opts(First_Instr) % rt_ir  % pc % ipcbnd,  &
          RT_opts(First_Instr) % rt_ir  % pc % ipcreg ) % predictindex )
      ALLOCATE(Inst_Chans(NumInstChans))
      Inst_Chans =  &
        RT_coefs(First_Instr) % coef_pccomp % pcreg( &
          RT_opts(First_Instr) % rt_ir % pc % ipcbnd, &
          RT_opts(First_Instr) % rt_ir % pc % ipcreg )  % predictindex
      First_Chan_Pos = 1
      Last_Chan_Pos = NPCScores
    ELSE
      NumInstChans = COUNT(RT_Params % Instrument_Number( &
        SatIndex,UsedChans % Channels) == Instrument)
      ALLOCATE(Inst_Chans(NumInstChans))
      IChan2 = 1
      First_Chan_Pos = 0
      Last_Chan_Pos = 0
      DO IChan1 = 1, UsedChans % NumChans
        IF (RT_Params % Instrument_Number(SatIndex, &
            UsedChans % Channels(IChan1)) == Instrument) THEN
          Inst_Chans(IChan2) = UsedChans % Channels(IChan1) - &
            RT_Params % First_Channel_for_Instrument(Instrument) + 1
          IF (First_Chan_Pos == 0) First_Chan_Pos = IChan1
          Last_Chan_Pos = IChan1
          IChan2 = IChan2 + 1
        END IF
      END DO
    END IF 

    IF (First_Chan_Pos == 0) CYCLE

    Allocate(RT_Chanprof(NumInstChans))
    RT_Chanprof(1:NumInstChans)%chan = Inst_Chans(1:NumInstChans)
    RT_Chanprof(1:NumInstChans)%prof = 1_JPIM

    !The inclusion of FastModel_Mode in the following call determines whether
    !the K variables are alloacted
    CALL NWPSAF_RTTOV12_Allocate ( &
      FastModel_Mode, & !in
      ASW_ALLOCATE,   & !in
      Instrument,     & !in
      NumInstChans,   & !in
      NPCScores,      & !in
      UsePCs,         & !in
      Radiance,       & !inout
      Radiance_K,     & !inout
      Transmission,   & !inout
      Transmission_K, & !inout
      Profiles_K,     & !inout
      Profiles_K_PC,  & !inout
      PCcomp,         & !inout
      PCcomp_K ,      & !inout
      Surf_Emiss,     & !inout
      Surf_Emiss_K,   & !inout
      CalcEmiss       ) !inout

    ! 2.2) Initialise emissivity arrays
    !---

    IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN

      !It would be nice not to have to read the file for every single ob, but
      !the problem is that we don't know in advance what months we would need to
      !read in. We only know the month on-the-fly when we read in the data.
      !We could read in all 12 months, but this may cause memory problems.
      !Perhaps we can think of a better way to do this in the future.
#ifdef _CompileEmissAtlas
      IF ( RT_Coefs(Instrument)%coef%id_sensor == sensor_id_mw  .OR. &
           RT_Coefs(Instrument)%coef%id_sensor == sensor_id_po  )  THEN ! MW
        CALL RTTOV_SETUP_EMIS_ATLAS(   &
          ErrorCode,                 & ! out
          RT_opts(Instrument),       & ! in
          Profiles(1) % Date(2),     & ! in
          atlas_type_mw,             & ! in
          emis_atlas,                & ! inout
          path=Atlas_Dir,            & ! in
          coefs=RT_coefs(Instrument) ) ! in
      ELSE ! IR
        CALL RTTOV_SETUP_EMIS_ATLAS(   &
          ErrorCode,                 & ! out
          RT_opts(Instrument),       & ! in
          Profiles(1) % Date(2),     & ! in
          atlas_type_ir,             & ! in
          emis_atlas,                & ! inout
          path=Atlas_Dir,            & ! in
          coefs=RT_coefs(Instrument) ) ! in
      END IF

      IF (ErrorCode >= errorstatus_fatal) THEN 
        Message(1)='Error in rttov_atlas_setup'
        WRITE( UNIT=Message(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorCode
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          Message,                & ! in
          ErrorStatus=StatusFatal ) ! in 
      END IF                

      Surf_Emiss(1:NumInstChans)%emis_in = 0.
      Surf_Emiss(1:NumInstChans)%emis_out = 0.
      Calcemiss(1:NumInstChans) = .FALSE.

      CALL RTTOV_GET_EMIS(      &
        ErrorCode,            & ! out
        RT_opts(Instrument),  & ! in
        RT_Chanprof,          & ! in
        Profiles,             & ! in
        RT_Coefs(Instrument), & ! in
        emis_atlas,           & ! in
        Surf_Emiss%emis_in    ) ! out

      IF (ErrorCode >= errorstatus_fatal) THEN 
        Message(1)='Error in rttov_get_emis'
        WRITE( UNIT=Message(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        Message,                & ! in
        ErrorStatus=StatusFatal ) ! in 
      END IF
#endif
    ELSE

      Surf_Emiss(1:NumInstChans)%emis_in = 0.
      Surf_Emiss(1:NumInstChans)%emis_out = 0.
      Calcemiss(1:NumInstChans) = .TRUE.

    END IF

    !Retain information for radiance to BT conversion
    IF (CloudyRetrieval) THEN
      RT_Params % tc1(First_Chan_Pos:Last_Chan_Pos) = &
            RT_Coefs(Instrument)%coef % ff_bco(Inst_Chans(1:NumInstChans))
      RT_Params % tc2(First_Chan_Pos:Last_Chan_Pos) = &
            RT_Coefs(Instrument)%coef % ff_bcs(Inst_Chans(1:NumInstChans))
      RT_Params % bcon1(First_Chan_Pos:Last_Chan_Pos) = &
            RT_Coefs(Instrument)%coef % fc_planck_c1 * &
            RT_Coefs(Instrument)%coef % ff_cwn(UsedChans%Channels)**3
      RT_Params % bcon2(First_Chan_Pos:Last_Chan_Pos) = &
            RT_Coefs(Instrument)%coef % fc_planck_c2 * &
            RT_Coefs(Instrument)%coef % ff_cwn(UsedChans%Channels)
    END IF

    !--------------------------------------------------------------------------
    !4) CALL RTTOV
    !--------------------------------------------------------------------------

    ForwardOrGradient: SELECT CASE ( Fastmodel_Mode )

    !4.1) RTTOV Direct call
    !----
    CASE( FastmodelMode_Forward ) ForwardOrGradient

!      CALL rttov_print_profile(profiles(1), lu=6)
!      CALL rttov_print_opts(RT_opts(Instrument), lu=6)
!      CALL rttov_print_info(RT_Coefs(Instrument), lu=6)

      IF ( UsePCs ) THEN
        CALL RTTOV_DIRECT( &
          ErrorCode,                            & ! out
          RT_Chanprof,                          & ! in
          RT_opts(Instrument),                  & ! in
          Profiles,                             & ! in
          RT_Coefs(Instrument),                 & ! in
          transmission,                         & ! inout
          Radiance,                             & ! inout
          calcemis=calcemiss(1:NumInstChans),   & ! in
          emissivity=Surf_Emiss(1:NumInstChans),& ! inout
          pccomp = PCcomp                       ) ! PC scores
      ELSE
        CALL RTTOV_DIRECT( &
          ErrorCode,                            & ! out
          RT_Chanprof,                          & ! in
          RT_opts(Instrument),                  & ! in
          Profiles,                             & ! in
          RT_Coefs(Instrument),                 & ! in
          transmission,                         & ! inout
          Radiance,                             & ! inout
          calcemis=calcemiss(1:NumInstChans),   & ! in
          emissivity=Surf_Emiss(1:NumInstChans) ) ! inout
      ENDIF

      IF (ErrorCode >= errorstatus_fatal) THEN 
        Message(1)='Error in RTTOV_DIRECT'
        WRITE( UNIT=Message(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorCode
        CALL NWPSAF_Report( &
          RoutineName,            & ! in
          Message,                & ! in
          ErrorStatus=StatusFatal ) ! in 
      END IF

      !NOTE: We deal with the radiance output from the direct call outside the
      !case statement because it is the same for forward and gradient calls

    !4.2) RTTOV K call
    !----
     CASE ( FastmodelMode_Gradient ) ForwardOrGradient

      !---- Initialise Emissivity_K arrays
      Surf_Emiss_K(:)%emis_in = 0.
      Surf_Emiss_K(:)%emis_out = 0.

      !CALL rttov_print_profile(profiles(1), lu=6)
      !CALL rttov_print_opts(RT_opts(Instrument), lu=6)

      IF ( UsePCs ) THEN
        ! Initialise pccomp_k%pcscores to 1
        pccomp_k%pcscores(:)=1.0
        ! For addradrec
        !pccomp_k%total_pccomp(:)=1.0 !switchrad == .false.
        !pccomp_k%bt_pccomp(:)=1.0 !switchrad == .true.
        CALL RTTOV_K( &
          ErrorCode   ,                             & ! out
          RT_Chanprof,                              & ! in
          RT_opts(Instrument),                      & ! in
          Profiles,                                 & ! in
          Profiles_K,                               & ! inout
          RT_Coefs(Instrument),                     & ! in
          transmission,                             & ! inout
          transmission_K,                           & ! inout
          Radiance,                                 & ! inout
          Radiance_K,                               & ! inout
          calcemis=calcemiss(1:NumInstChans),       & ! in
          emissivity=Surf_Emiss(1:NumInstChans),    & ! inout
          emissivity_k=Surf_Emiss_K(1:NumInstChans),& ! inout
          PCcomp=PCcomp,                            & ! inout
          PCcomp_K=PCcomp_K,                        & ! inout
          Profiles_K_PC=Profiles_K_PC               ) ! inout 
!      CALL rttov_print_profile(profiles_k_pc(1), lu=6)
      ELSE
        ! Initialise radiance_k % bt to 1K.
        ! Used if RTTOV_Opts%switchrad == .true.
        Radiance_k % bt = 1.0
        ! Initialise radiance_k % total to 1 radiance unit.
        ! Used if RTTOV_Opts%switchrad == .false.
        Radiance_k % total = 1.0
        CALL RTTOV_K( &
          ErrorCode   ,                             & ! out
          RT_Chanprof,                              & ! in
          RT_opts(Instrument),                      & ! in
          Profiles,                                 & ! in
          Profiles_K,                               & ! inout
          RT_Coefs(Instrument),                     & ! in
          transmission,                             & ! inout
          transmission_K,                           & ! inout
          Radiance,                                 & ! inout
          Radiance_K,                               & ! inout
          calcemis=calcemiss(1:NumInstChans),       & ! in
          emissivity=Surf_Emiss(1:NumInstChans),    & ! inout
          emissivity_k=Surf_Emiss_K(1:NumInstChans) ) ! inout
      ENDIF

      IF (ErrorCode >= errorstatus_fatal) THEN 
        Message(1)='Error in RTTOV_K'
        WRITE( UNIT=Message(2),FMT='(A,I3)' ) &
          'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        Message,                & ! in
        ErrorStatus=StatusFatal ) ! in 
      END IF



      !----
      CALL NWPSAF_RTTOV12_GetHMatrix (  &
             Profiles,        & !in
             Profiles_K,      & !in
             Profiles_K_PC,   & !in
             Instrument,      & !in
             NumInstChans,    & !in
             First_Chan_Pos,  & !in
             Last_Chan_Pos,   & !in
             RT_Params        ) !inout

    CASE DEFAULT ForwardOrGradient

      WRITE( UNIT=Message(1),FMT='(A)' ) &
        'Incorrect Fastmodel Mode in Fastmodel Interface'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        Message,                & ! in
        ErrorStatus=StatusFatal ) ! in 

    END SELECT ForwardOrGradient

    !4.3) Unpack Radiance/BT from RTTOV/RTTOVK into RT_Params structure
    !----

    
    IF ( UsePCs ) THEN
      RT_Params % PCScores(1:NPCScores) = PCcomp % PCscores (1:NPCScores)
    ELSE IF ( CalcRadiance ) THEN
      RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
        Radiance%total(1:NumInstChans) 
    ELSE
      RT_Params % TotalBTs(First_Chan_Pos:Last_Chan_Pos) = &
        Radiance%bt(1:NumInstChans)
      IF (CloudyRetrieval) THEN
        RT_Params % TotalRadiances(First_Chan_Pos:Last_Chan_Pos) = &
        Radiance%total(1:NumInstChans) 
        RT_Params % CloudyRadiances(First_Chan_Pos:Last_Chan_Pos, &
          1:profiles(1) % nlayers) = &
          TRANSPOSE(Radiance%overcast(:,1:NumInstChans))
      END IF
    END IF

    RT_Params % RTEmissivity(1:NumInstChans) = Surf_Emiss(1:NumInstChans)%emis_out

    !4.4) Tidy up and deallocate
    !-------
    DEALLOCATE ( Inst_Chans )
    DEALLOCATE ( RT_Chanprof )

    !The ASW_DEALLOCATE switch tells the routine to deallocate anything that was
    !previously allocated depending on FastModel_Mode.
    CALL NWPSAF_RTTOV12_Allocate ( &
      FastModel_Mode, & !in
      ASW_DEALLOCATE, & !in
      Instrument,     & !in
      NumInstChans,   & !in
      NPCScores,      & !in
      UsePCs,         & !in
      Radiance,       & !inout
      Radiance_K,     & !inout
      Transmission,   & !inout
      Transmission_K, & !inout
      Profiles_K,     & !inout
      Profiles_K_PC,  & !inout
      PCcomp,         & !inout
      PCcomp_K,       & !inout
      Surf_Emiss,     & !inout
      Surf_Emiss_K,   & !inout
      CalcEmiss       ) !inout

    IF (Use_EmisAtlas .AND. Profiles(1)%skin%surftype /= RTsea) THEN
#ifdef _CompileEmissAtlas
      CALL rttov_deallocate_emis_atlas( emis_atlas )
#endif
    END IF

  END DO Instrument_Loop


  ! Deallocate profile arrays

  CALL RTTOV_ALLOC_PROF(  &
    DeallocateError,      &
    1,                    &
    Profiles(1),          &
    Num_RTLevels,         & 
    RT_Opts(1),           &
    asw = ASW_DEALLOCATE, &
    coefs = RT_coefs(1)   )

  IF (DeallocateError >= errorstatus_fatal) THEN 
    Message(1)='Error in RTTOV_ALLOC_PROF (Deallocate)'
    WRITE( UNIT=Message(2),FMT='(A,I3)' ) &
      'Error is Code ',DeallocateError
  CALL NWPSAF_Report(       &
    RoutineName,            & ! in
    Message,                & ! in
    ErrorStatus=StatusFatal ) ! in 
  END IF


ENDIF RTTOV_FastmodelMode

End Subroutine NWPSAF_RTTOV12_Interface
