Subroutine NWPSAF_RTTOV11_Allocate ( &
  FastModel_Mode, & !in
  ASW,            & !in
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
! Description: Interface between IASI 1DVar code and RTTOV
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
!         11/04/12 New. Code taken from IASI_RTTOV_Allocate. Fiona Smith
!
! Code Description:
!   Language:   Fortran 95
!   Software Standards:
!
! End of header -------------------------------------------------------

USE NWPSAFMod_Params, ONLY : &
  StatusFatal,               &
  MwClwRetrieval

USE NWPSAFMod_RTmodel, ONLY : &
  FastmodelMode_Gradient,  &
  Num_RTlevels,            &
  RT_coefs,             &
  RT_opts

USE rttov_types, ONLY:  &
  profile_type,      &
  radiance_type,     &
  transmission_type, &
  rttov_pccomp, &
  rttov_emissivity

Use parkind1, ONLY : &
  jpim

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'rttov_alloc_pccomp.interface'
INCLUDE 'rttov_alloc_rad.interface'
INCLUDE 'rttov_alloc_transmission.interface'
INCLUDE 'rttov_alloc_prof.interface'

! Subroutine arguments:
INTEGER, INTENT(IN)             :: Fastmodel_Mode  ! Forward/Gradient
INTEGER, INTENT(IN)             :: ASW             ! Allocate/Dealloc
INTEGER, INTENT(IN)             :: Instrument
INTEGER, INTENT(IN)             :: NumInstChans
INTEGER, INTENT(IN)             :: NPCScores
LOGICAL, INTENT(IN)             :: UsePCs
TYPE(profile_type), ALLOCATABLE, INTENT(INOUT)     :: Profiles_K(:)
TYPE(profile_type), INTENT(INOUT)                  :: Profiles_K_PC(:)
TYPE(radiance_type), INTENT(INOUT)                 :: Radiance
TYPE(radiance_type), INTENT(INOUT)                 :: Radiance_K
TYPE(transmission_type), INTENT(INOUT)             :: Transmission
TYPE(transmission_type), INTENT(INOUT)             :: Transmission_K
TYPE(rttov_pccomp), INTENT(INOUT)                  :: PCcomp
TYPE(rttov_pccomp), INTENT(INOUT)                  :: PCcomp_K
TYPE(rttov_emissivity), ALLOCATABLE, INTENT(INOUT) :: Surf_Emiss(:) 
TYPE(rttov_emissivity), ALLOCATABLE, INTENT(INOUT) :: Surf_Emiss_K(:)
LOGICAL, ALLOCATABLE, INTENT(INOUT)                :: calcemiss(:)
 
! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV11_Allocate"

! Local variables:
INTEGER(kind=jpim) :: ErrorCode
INTEGER            :: Num_RTLayers

CHARACTER(LEN=80) :: ErrorMessage(2)  ! Message for NWPSAF_Report

INTEGER, PARAMETER :: ASW_ALLOCATE = 1
INTEGER, PARAMETER :: ASW_DEALLOCATE = 0
INTEGER :: I
!---------------------------------------------------------------

ErrorCode = 0
ErrorMessage(:)=' '
Num_RTLayers = Num_RTLevels-1

!---------------------------------------------
!1) These are required for forward and K calls
!---------------------------------------------

!----- Allocate/Deallocate Emissivity arrays
IF ( ASW == ASW_ALLOCATE ) THEN
  ALLOCATE(Surf_Emiss(NumInstChans))
  ALLOCATE(CalcEmiss(NumInstChans))
ELSE IF ( ASW == ASW_DEALLOCATE ) THEN
  DEALLOCATE(Surf_Emiss)
  DEALLOCATE(CalcEmiss)
END IF


!----- Allocate/Deallocate Radiance array and nullify/initialise to zero
CALL RTTOV_ALLOC_RAD( &
  ErrorCode,          & ! out
  NumInstChans,       & ! in
  Radiance,           & ! in/out
  Num_RTLayers,       & ! in
  ASW,                & ! in
  init=.TRUE.         ) ! in

IF (ErrorCode /= 0) THEN
  ErrorMessage(1)='Error in RTTOV_ALLOC_RAD'
  WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
    'Error is Code ',ErrorCode
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF
!----- Allocate/Deallocate transmission arrays
CALL RTTOV_ALLOC_TRANSMISSION( &
  ErrorCode,    &
  Transmission, &
  Num_RTLayers, &
  NumInstChans, &
  ASW,          &
  init = .TRUE. )

IF (ErrorCode /= 0) THEN
  ErrorMessage(1)='Error in RTTOV_ALLOC_TRANSMISSION'
  WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
    'Error is Code ',ErrorCode
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

!---- Allocate/Deallocate PC component arrays
IF ( UsePCs ) THEN
  CALL rttov_alloc_pccomp( &
    ErrorCode,   &
    PCcomp,      &
    NPCScores,   &
    ASW,         &
    init = .TRUE.)

  IF (ErrorCode /= 0) THEN
    ErrorMessage(1)='Error in RTTOV_ALLOC_PCCOMP'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
      'Error is Code ',ErrorCode
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF
END IF


!--------------------------------------
!2) These are required for K calls only
!--------------------------------------

IF ( FastModel_Mode == FastModelMode_Gradient ) THEN

  IF (ASW == ASW_ALLOCATE) THEN
    ALLOCATE(Surf_Emiss_K(NumInstChans))
    ALLOCATE ( Profiles_K(NumInstChans))
  END IF 

  !---- Allocate/Dealloacate radiance_K
  CALL RTTOV_ALLOC_RAD( &
    ErrorCode,    & ! out
    NumInstChans, & ! in
    Radiance_K,   & ! in/out
    Num_RTLayers, & ! in
    ASW,          & ! in
    init=.TRUE.   ) ! in
  IF (ErrorCode /= 0) THEN
    ErrorMessage(1)='Error in RTTOV_ALLOC_RAD'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
    'Error is Code ',ErrorCode
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  !---- Allocate/Deallocate transmission_K
  CALL rttov_alloc_transmission( &
    ErrorCode,     &
    transmission_K, &
    Num_RTLayers,   &
    NumInstChans,   &
    ASW,            &
    init=.TRUE.     )

  IF (ErrorCode /= 0) THEN
    ErrorMessage(1)='Error in RTTOV_ALLOC_TRANSMISSION'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
    'Error is Code ',ErrorCode
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
  END IF

  IF ( UsePCs ) THEN
    !---- Allocate/deallocate pccomp_k
    CALL rttov_alloc_pccomp( &
      ErrorCode,   &
      PCcomp_K,     &
      NPCScores,    &
      ASW,          &
      init = .TRUE. )


    IF (ErrorCode /= 0) THEN
      ErrorMessage(1)='Error in RTTOV_ALLOC_PCCOMP'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
        'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END IF


  CALL rttov_alloc_prof ( &
    ErrorCode,     &
    NumInstChans,  & !nprofiles
    profiles_K,    &
    Num_RTLevels,  &
    RT_opts(Instrument),  &
    ASW,           &
    RT_coefs(Instrument), &
    init=.true.    )
    
!PM Case where we retrieve LWP values but start from a clear background
   IF( ALL(RT_opts(:) % rt_mw % clw_data .eqv. .FALSE.) .and. MwClwRetrieval ) THEN
      Do I=1,NumInstChans
         Allocate( Profiles_K(I) % clw(Num_RTLevels) )
         Profiles_K(I) % clw(:) = 0.
      END DO
    END IF    
!PM    

  IF (ErrorCode /= 0) THEN
    ErrorMessage(1)='Error in RTTOV_ALLOC_PROF'
    WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
      'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
  END IF

  IF (UsePCs) THEN
    CALL rttov_alloc_prof ( &
      ErrorCode,     &
      NPCScores,     & !nprofiles
      profiles_K_PC, &
      Num_RTLevels,  &
      RT_opts(Instrument),  &
      ASW,           &
      RT_coefs(Instrument), &
      init=.true.    )
    IF (ErrorCode /= 0) THEN
      ErrorMessage(1)='Error in RTTOV_ALLOC_PROF'
      WRITE( UNIT=ErrorMessage(2),FMT='(A,I3)' ) &
        'Error is Code ',ErrorCode
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF
  END IF


  IF (ASW == ASW_DEALLOCATE) THEN
    DEALLOCATE(Surf_Emiss_K)
    DEALLOCATE (Profiles_K)
  END IF


END IF

END SUBROUTINE NWPSAF_RTTOV11_Allocate

