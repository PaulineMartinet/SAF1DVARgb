PROGRAM pc_to_spec

!------------------------------------------------------------------------------
! Description: Read in Observations for 1DVar, Convert PC scores to Radiances
!              or BTs and write them out again.
!              Assumes the use of PC-RTTOV PC coefficients.
!
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

Use rttov_types, Only:   &
    rttov_options, &
    rttov_coefs

IMPLICIT NONE

include 'rttov_read_coefs.interface'

TYPE observations_type
  REAL, ALLOCATABLE :: BriTemp(:)
  REAL, ALLOCATABLE :: Radiances(:)
  REAL, ALLOCATABLE :: PC_scores(:)
  INTEGER :: Date(3)
  REAL    :: Elevation
  INTEGER :: ID
  REAL    :: Latitude
  REAL    :: Longitude
  INTEGER :: ObsType
  INTEGER :: SatId
  REAL    :: SatZenith
  REAL    :: SolarZenith
  INTEGER :: Surface
END TYPE

! Declarations

!General
INTEGER :: i,j,ii

!Input/Output file info and storage variables
CHARACTER(LEN=*), PARAMETER  :: infile = &
  '/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/output_simrad_rt11/70L_200_PCFromRR_truth.dat'
CHARACTER(LEN=*), PARAMETER  :: outfile = &
  '/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/output_simrad_rt11/70L_8461_RRFromPCFromRR_truth.dat'

INTEGER, PARAMETER :: inunit=11
INTEGER, PARAMETER :: outunit=12      !File units for I/O

INTEGER :: ReadStatus          ! IOSTAT keyword value
INTEGER :: NumHeaderLines = 10 ! # of lines of informational header in
                               ! ObsFileName
INTEGER :: NumHeaderLines2 = 3 ! # of lines of 2nd header in ObsFileName
INTEGER :: obnumber            ! # of observations read in
REAL :: Vars(3)                ! For splitting up header lines
CHARACTER(LEN=80) :: String    ! Observation information
CHARACTER(LEN=80) :: Header    ! Header line text

!Observation information
TYPE(observations_type) :: Observation
INTEGER :: First_Channel
INTEGER :: Last_Channel
INTEGER :: NumObs               ! # of soundings in Observation structure
INTEGER, ALLOCATABLE :: Channel_Number(:)
INTEGER :: Series
INTEGER :: Platform
INTEGER :: WMO
INTEGER :: Num_Instruments
INTEGER :: Subtype
INTEGER, PARAMETER :: NumChannels = 8461
INTEGER, PARAMETER  :: ipcreg=3
INTEGER, PARAMETER  :: ipcbnd=1

!RTTOV Variables
TYPE(rttov_coefs) :: coef
TYPE(rttov_options) :: opts
INTEGER :: Err
INTEGER :: NPCScores
CHARACTER(LEN=*), PARAMETER  :: coeffile = &
  '/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/test_simspec/rtcoef_metop_2_iasi.H5'
CHARACTER(LEN=*), PARAMETER  :: pcfile = &
  '/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/test_simspec/pccoef_metop_2_iasi.dat'
REAL, POINTER :: evec(:,:)

LOGICAL :: convert_to_bt=.false.
CHARACTER(LEN=10)  :: VarChar

!-----------------------------------------------------------------------------

!-------------------------------------------------
!0. Initialise constants and variables, open files
!-------------------------------------------------

Vars(:) = 0.0
obnumber = 0

WRITE(*,*) 'Opening Files'

OPEN( UNIT   = inunit, &
      FILE   = infile,&
      ACTION = "READ", &
      STATUS = "OLD", &
      FORM   = "FORMATTED" )

OPEN( UNIT   = outunit, &
      FILE   = outfile,&
      ACTION = "WRITE", &
      STATUS = "REPLACE", &
      FORM   = "FORMATTED" )


!----- Write obs file header
!----- Start with 10 lines for comments
Write(outunit,'(a)') 'This is a simulated clear-air IASI observation dataset.'
Write(outunit,'(a)') 'Generated by sim_rad_iasi_rttov10_70L.F90, using RTTOV10 on 101 levels.'
Write(outunit,'(a)') 'Atmospheric state profiles were taken from the Met Office'
Write(outunit,'(a)') '70-level dataset Model_20101123_0000_Sea.dat interpolated inside RTTOV10.'
Write(outunit,'(a)') 'No observation errors have been added'
Write(outunit,'(a)') 'The observations have been converted to PC space using PC-RTTOV coeffs'
Write(outunit,'(a)') 'using the program rad_to_pc.f90, then back to RRs using pc_to_rad.f90'
Write(outunit,'(a)') ''
Write(outunit,'(a)') ''
Write(outunit,'(a)') ''


!---------------------------------
! 1. Read in RTTOV PC coefficients
!---------------------------------

opts % rt_ir % pc % ipcreg = ipcreg
opts % rt_ir % pc % addpc  = .true.
opts % rt_ir % pc % ipcbnd = ipcbnd
WRITE(*,*) 'Opening RTTOV Coefficient files'

Call rttov_read_coefs( err,                &
                       coef,               &
                       opts,               &
                       file_coef=coeffile, &
                       file_pccoef=pcfile  )
IF (Err /= 0) CALL error_trap(' Error Reading RTTOV coeffs ')


evec => coef % coef_pccomp % eigen(1) % eigenvectors(:,:)

!-----------------------------------
!1. Initialize observation structure
!-----------------------------------


!--------------------------------------------------------
!1.1 Open the Observations File and Read in Header info
!--------------------------------------------------------

WRITE(*,*) 'Reading input file header'


! Allow NumHeaderLines Lines for an informational header
ReadHeader : DO i = 1, NumHeaderLines
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (ReadStatus /= 0) CALL error_trap(' Error With Header1 ')
END DO ReadHeader

! Read in number of observations
READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (Readstatus /= 0) CALL error_trap(' Problem with numobs string ')

Vars(1) = 0.
CALL ReadHeaders (1,String,Vars(1:1),ReadStatus)
IF (Readstatus /= 0) CALL error_trap(' Problem with numobs string (1)')
NumObs = NINT(Vars(1))

! Read in number of channels per observation
! Note: Number of channels is assumed to be the same for each observation
READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (Readstatus /= 0) CALL error_trap(' Problem with numchans string ')

Vars(1) = 0.
CALL ReadHeaders (1,String,Vars(1:1),ReadStatus)
IF (Readstatus /= 0) CALL error_trap(' Problem with numobs string (1)')
NPCScores = NINT(Vars(1))

! Read in number of separate instruments that make up the obs type
! Note: This is assumed to be the same for each observation
READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
IF (Readstatus /= 0) CALL error_trap(' Problem with numinst string (0)')

Vars(1) = 0.
CALL ReadHeaders (1,String,Vars(1:1),ReadStatus)
IF (Readstatus /= 0) CALL error_trap(' Problem with numobs string (0)')
Num_Instruments = NINT(Vars(1))


! Three more header lines
ReadHeader2 : DO i = 1, NumHeaderLines2
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (Readstatus /= 0) CALL error_trap(' Problem with header2')
END DO ReadHeader2


String = ADJUSTL(Header)
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
    IF ( VarChar(1:8) /= 'PC Score' ) THEN
      Write(*,*) 'Observations are already Radiance/BT!'
      STOP
    END IF
  ELSE
      CALL error_trap('Error Reading Observation Units')
  END IF
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
END IF


ALLOCATE(Channel_Number(NPCScores))
Channel_Number(:) = 0
READ (UNIT=inunit, FMT=* , IOSTAT=ReadStatus) &
  Series,   &
  Platform, &
  SubType,  &
  First_Channel,  &
  Last_Channel,   &
  WMO

IF (ReadStatus /= 0 .OR. Series < 0 .OR. &
    Platform < 0 .OR. SubType < 0 .OR. &
    First_Channel <= 0 .OR. Last_Channel > NumChannels .OR. &
    Last_Channel < First_Channel) &
  CALL error_trap(' Error reading satellite information')


READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
String = ADJUSTL(String)
IF (String(1:8) == 'Channels' .or. String(1:3) == 'PCs' ) THEN
  READ(UNIT=inunit, FMT=*, IOSTAT=ReadStatus) Channel_Number(:)

  ! This line reads a seperator
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (ReadStatus /= 0 .OR. ANY(Channel_Number(:) <= 0)) &
    CALL error_trap('Error reading channel numbers')
END IF

!Write this information back out to the output file
Write(outunit,'(a,i6)') 'Number of Observations in File:',NumObs
Write(outunit,'(a,i8)') 'No. of Chans per Observation:',NumChannels
Write(outunit,'(a)') 'Number of instruments making up observations : 1'
Write(outunit,'(a)') '*** In the following Series, Platform and Instrument are defined  ***'
Write(outunit,'(a)') '*** according to the relevant RT Model definitions (if required): ***'
If (convert_to_bt ) THEN
  Write(outunit,'(a)') 'Units: BT'
Else
  Write(outunit,'(a)') 'Units: Radiance'
End If
Write(outunit,'(a)') 'Sat. Series   Platform   Instrument First_Channel   Last_Channel  Sat ID'
Write(outunit,'(i5,i13,i12,i12,i16,i12)') Series,Platform,Subtype,1,NumChannels,WMO
Write(outunit,'(a)') 'Channels'
Write(outunit,'(16i5)') (/ (i,i=1,NumChannels) /)
Write(outunit,'(a)') '----------------------------------------------------------------------'


! Allocate relevent parts of the Observations structure and initialise
! associated header elements.

ALLOCATE( Observation% BriTemp  ( NumChannels ) )
ALLOCATE( Observation% Radiances( NumChannels ) )
ALLOCATE( Observation% PC_scores( NPCScores   ) )
Observation% BriTemp(:)  = 0.0
Observation% Radiances(:)= 0.0
Observation% PC_scores(:)= 0.0
Observation% Date(:)     = 0
Observation% Elevation   = 0.0
Observation% ID          = 0
Observation% Latitude    = 0.0
Observation% Longitude   = 0.0
Observation% ObsType     = 0
Observation% SatId       = 0
Observation% SatZenith   = 0.0
Observation% SolarZenith = 0.0
Observation% Surface     = 0

!--------------------------------------------------------
!1.2 Read the observations in
!--------------------------------------------------------
WRITE(*,*) 'Reading observations'


ReadData: DO WHILE (obnumber < NumObs)

  obnumber = obnumber + 1
  Observation%PC_scores(:)=0.0

  ! Observation identifiers:
  !-------
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus /= 0) CALL error_trap('Error reading first ob string',num=obnumber)
  CALL ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observation% ID      = NINT(Vars(1))
  Observation% ObsType = NINT(Vars(2))
  Observation% SatID   = NINT(Vars(3))

  IF (ReadStatus /= 0) CALL error_trap('Error reading first ob string',num=obnumber)

  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
  ! Temporal information (required if reading in emissivity atlas)
  !-------
  IF (String(1:4) == 'Year') THEN
    IF (ReadStatus /= 0) CALL error_trap('Error reading date',num=obnumber)
      CALL ReadHeaders (3,String,Vars(1:3),ReadStatus)
      Observation% Date(1) = NINT(Vars(1))
      Observation% Date(2) = NINT(Vars(2))
      Observation% Date(3) = NINT(Vars(3))
      IF (ReadStatus /= 0) CALL error_trap('Error reading date info',num=obnumber)
      READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
  END IF

  ! Navigation information:
  !-------
  IF (ReadStatus /= 0) CALL error_trap('Error reading lat/lon',num=obnumber)
  CALL ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observation% Latitude  = Vars(1)
  Observation% Longitude = Vars(2)
  Observation% Elevation = Vars(3)
  IF (ReadStatus /= 0) CALL error_trap('Error reading lat/lon',num=obnumber)

  ! More navigation information:
  !-------
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) String
  IF (ReadStatus /= 0) CALL error_trap('Error reading surface info',num=obnumber)
  CALL ReadHeaders (3,String,Vars(1:3),ReadStatus)
  Observation% Surface     = NINT(Vars(1))
  Observation% SatZenith   = Vars(2)
  Observation% SolarZenith = Vars(3)
  IF (ReadStatus /= 0) CALL error_trap('Error reading surface info',num=obnumber)

  ! Header message:
  !-------
  READ (UNIT=inunit, FMT='(A80)', IOSTAT=ReadStatus) Header
  IF (ReadStatus /= 0) CALL error_trap('Error reading header line',num=obnumber)

  ! Brightness temperatures:
  !-------
  READ (UNIT=inunit, FMT=*, IOSTAT=ReadStatus) Observation% PC_Scores(:)
  IF (ReadStatus /= 0) CALL error_trap('Error reading BTs',num=obnumber)

!--------------------------------------------------------
!2) Convert data to PCs
!--------------------------------------------------------

!--------------------------------------------------------
!2.1 BT to Radiance
!--------------------------------------------------------
  Observation%Radiances(:)=0.0
  WRITE(*,*) 'Converting to Radiance'
  DO i=1, NumChannels
    DO j=1, NPCScores
      Observation%Radiances(i)= Observation%Radiances(i) + &
        evec(i,j) * Observation%PC_scores(j)
    END DO
    Observation%Radiances(i)= Observation%Radiances(i) * &
                              coef%coef_pccomp%noise_in(i)
  END DO

!--------------------------------------------------------
!2.2 Radiance to PCs
!--------------------------------------------------------

  IF ( convert_to_bt ) THEN
    WRITE(*,*) 'Converting to BT'
    Observation%BriTemp(:) = coef%coef%planck2(:) /  &
      LOG(1 + coef%coef%planck1(:) / Observation%Radiances(:)) - &
      coef%coef%ff_bco(:) / coef%coef%ff_bcs(:)
  END IF


!--------------------------------------------------------
!3) Write the data back out
!--------------------------------------------------------


  write(*,'(A,I6)')'Ob number:', obnumber
  write(outunit,'(a,i15,a,i11,a,i6)') &
    'Obs ID:',Observation%ID,' Obs Type:',Observation%ObsType, &
    ' Satellite ID:', Observation%SatID
  !HERE WE WOULD WRITE OUT THE DATE INFO
  write(outunit,'(a,f9.3,a,f9.3,a,f7.1)') &
    'Latitude:',Observation%Latitude,' Longitude:',Observation%Longitude, &
    ' Elevation:', Observation%Elevation
  write(outunit,'(a,i4,a,f9.3,a,f9.3)') &
    'Surface Type:',Observation%surface, &
    ' Sat Zen Angle:',Observation%satzenith, &
    ' Solar Zen. Ang.:',Observation%solarzenith
  IF ( convert_to_bt ) THEN
    write(outunit,'(a)') 'Brightness Temperatures:'
    write(outunit,'(6f13.3)') Observation%BriTemp(:)
  ELSE
    write(outunit,'(a)') 'Radiances:'
    write(outunit,'(6f13.3)') Observation%Radiances(:)
  END IF


END DO ReadData


!--------------------------------------------------------
!4) Tidy Up
!--------------------------------------------------------

WRITE(*,*) 'All done, tidying up now'

DEALLOCATE(Observation%BriTemp)
DEALLOCATE(Observation%Radiances)
DEALLOCATE(Observation%PC_Scores)


CLOSE(inunit)
CLOSE(outunit)

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE ReadHeaders ( &
    NumVars,       & ! in
    String,        & ! inout
    Vars,          & ! out
    Status)          ! inout

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
! Description: Read in colon-separated fields in observation and
!              background file headers.
!
!
! Method:
!     Obvious.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 2.2      30/04/02  Original Version
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

IMPLICIT NONE

! Subroutine Arguments:

INTEGER          , INTENT(IN)    :: NumVars
CHARACTER (LEN=*), INTENT(INOUT) :: String
REAL             , INTENT(OUT)   :: Vars(NumVars)
INTEGER          , INTENT(OUT)   :: Status

! Internal Variables:

INTEGER :: I
INTEGER :: II

! Program:

Status = 0
Vars(:) = 0.

DO I=1,NumVars
  II=Index(String,':')
  IF (II == 0) Then
     Status = -1
     EXIT
  END IF
  String=String(II+1:)
  READ(String,*,IOSTAT=Status) Vars(I)
  IF (Status /= 0) EXIT
END DO


END Subroutine ReadHeaders

!-----------------------------------------------------------------------------
SUBROUTINE error_trap (message, num)

IMPLICIT NONE

CHARACTER (LEN=*), INTENT(IN) :: message
INTEGER, OPTIONAL, INTENT(IN) :: num

CHARACTER (LEN=100) :: ErrMsg

IF (PRESENT(num)) THEN
    WRITE (ErrMsg, '(A,I8)') , message, ObNumber
ELSE
  ErrMsg=message
ENDIF

Write(*,'(A)') ErrMsg
STOP

END SUBROUTINE error_trap
!-----------------------------------------------------------------------------


END PROGRAM pc_to_spec


