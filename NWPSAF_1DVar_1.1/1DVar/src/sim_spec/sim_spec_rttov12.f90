PROGRAM sim_spec

! Description:
! Reads Met Office 70L model profiles and generates
! simulated IASI brightness temperatures based on these profiles using RTTOV9.
! Model profiles are interpolated onto 43 RTTOV levels inside RTTOV.
!
! The input file is in the same format as the background file - in fact it is the "truth"
! background file. There are two, one for sea points and one for land
!   /data/local/frfh/sim_spec_bg/Model_20101123_0000_Sea.dat
!   /data/local/frfh/sim_spec_bg/Model_20101123_0000_Land.dat
!
! Produces output file suitable for use in Met Office 1D-Var.
!
! Clear only for now, as nobody understands RTTOV cloud processing in RTTOV9!
!
!----- Imported Type Definitions
Use rttov_types, Only:   &
    rttov_profile,        &
    rttov_transmission,   &
    rttov_radiance,       &
    rttov_pccomp,        &
    rttov_coefs,        &
    rttov_options,      &
    rttov_chanprof,    &
    rttov_emissivity

Use parkind1, Only : &
    jpim,            &
    jprb

Use rttov_const, Only : &
    gas_id_watervapour, &
    gas_id_ozone

IMPLICIT NONE

include "rttov_direct.interface"
include "rttov_alloc_transmission.interface"
include "rttov_alloc_rad.interface"
include "rttov_alloc_pccomp.interface"
include "rttov_alloc_prof.interface"
include "rttov_dealloc_coefs.interface"
include "rttov_read_coefs.interface"
include "rttov_user_options_checkinput.interface"
include "rttov_print_profile.interface"
include "rttov_print_opts.interface"

!----- Local parameters:
integer(kind=jpim), parameter :: sat_id=4
integer(kind=jpim), parameter :: max_iasi_channels=8461   ! Total number of IASI channels
!-----TEST----------
integer(kind=jpim), parameter   :: idim=1                ! Max. no. of profiles to read - sea
!---------NEW FILE MAR 2011
!integer(kind=jpim), parameter   :: idim=4348              ! Max. no. of profiles to read - sea
! integer(kind=jpim), parameter   :: idim=2738             ! Max. no. of profiles to read - land
!-------------------------------
integer(kind=jpim), parameter :: nwplevels=70             ! No. of model levels
integer(kind=jpim), parameter :: nwplayers=69             ! No. of model levels - 1
integer(kind=jpim), parameter :: nslay = 4                ! No. of model surface layers

INTEGER(kind=jpim), PARAMETER :: ipcreg=3               !Number of PC predictors
INTEGER(kind=jpim), PARAMETER :: npcscores=200            !number of PC scores

real(kind=jprb), parameter    :: H2O_MassMixToPPMV = 1.6078e6

integer, parameter :: METOP_ID = 10
integer, parameter :: IASI_ID = 16
integer, parameter :: METOP_SAT_ID = 2

LOGICAL, parameter :: add_error=.false.
LOGICAL, parameter :: add_bias=.false.
LOGICAL, parameter :: usepcs=.true.
LOGICAL, parameter :: produceradiances=.true.
LOGICAL, parameter :: reconstruct=.false.
Integer, parameter :: surftype=1   ! 0=land    1=sea    2=snow/ice


CHARACTER(len=200) :: biasfile='/home/h03/frfh/meto1dvar/idl/bias.txt'
!CHARACTER(len=200) :: noisefile='/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/IASI_COEFFS_DIR_RAD/Rmatrix_8461_instnoise_PCRTTOV.out'
CHARACTER(len=200) :: noisefile='/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/IASI_COEFFS_DIR_BT/Rmatrix_8461_instnoise_PCRTTOV.out'

CHARACTER(len=200) :: profilefile='/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/Sample_Background/truth_70L.dat'
CHARACTER(len=200) :: outfile='/data/local/frfh/NWPSAF_1DVar/r1508_HTFRTC_LHR/1DVar/output_simrad_rt12/70L_200_PC_truth.dat'

!----- Miscellaneous variables
integer(kind=jpim)            :: jlev                              ! Level index variable
integer(kind=jpim)            :: iasi_channels(max_iasi_channels)  ! List of IASI channels
integer(kind=jpim)            :: i,j                               ! Loop index variable
integer                       :: surftype_1dvar
Integer                       :: Rtype
Real, Allocatable             :: Rdiag(:)
Real, Allocatable             :: Rread(:)
Real, Allocatable             :: Rchannels(:)
Real                          :: rndg
Real                          :: lat, lon
Integer                       :: ii,jj
Character(len=100)            :: text, string1, string2
real(kind=jprb)               :: qppmv(nwplevels)
real(kind=jprb)               :: pprof(nwplevels)
real(kind=jprb)               :: tprof(nwplevels)
real(kind=jprb)               :: qprof(nwplevels)
real(kind=jprb)               :: oprof(nwplevels)
real(kind=jprb)               :: s2qppmv
integer(kind=jpim)            :: n_chans             ! No. of channels to use
integer(kind=jpim)            :: numinstchans        ! No. of channels to use
integer(kind=jpim)            :: nchans_out             !For post-processing
integer                       :: numRchannels
integer                       :: numRelements
integer                       :: Inverse
integer                       :: biased_chans(183)
real                          :: bias(183)
INTEGER :: check
INTEGER :: chan

!----- RTTOV_SETUP interface
integer(kind=jpim), parameter     :: err_unit = 0        ! Logical unit for error messages
integer(kind=jpim), parameter     :: verbosity_level = 3 ! Maximum verbosity
integer(kind=jpim)                :: ninst               ! Number of instruments
integer(kind=jpim)                :: nprofiles           ! Number of profiles per call
integer(kind=jpim), allocatable   :: instrument(:,:)     ! Platform/sat/instr id
integer(kind=jpim), allocatable   :: channels(:,:)       ! Channels to use
integer(kind=jpim), allocatable   :: errorstatus(:)      ! Return error codes
Type(rttov_coefs), allocatable    :: coef(:)             ! RTTOV coefficients
Type(rttov_options), allocatable  :: options(:)        ! RTTOV options

!----- RTTOV inputs
Type(rttov_chanprof), allocatable  :: chanprof(:)       ! Channels per profile
Type(rttov_profile), allocatable   :: profiles(:)         ! Atmospheric profiles
Type(rttov_radiance)               :: radiance            ! radiance structure
Type(rttov_pccomp)                 :: pccomp              ! PC scores
Type(rttov_transmission)           :: transmission        ! transmission structure
TYPE(rttov_emissivity)             :: Surf_Emiss(max_iasi_channels) 
logical                            :: calcemis(max_iasi_channels) ! Emissivity calc. switch
integer(kind=jpim)                 :: qunits              ! (1=ppmv, 2=kg/kg, 3=RH)

INTEGER, PARAMETER :: ASW_ALLOCATE = 1
INTEGER, PARAMETER :: ASW_DEALLOCATE = 0



!-------------------------------------------------------------------------------------------

!-----------------
!1) General Set Up
!-----------------

!----- Generate IASI channel selection

!      Currently selects all IASI channels
n_chans=max_iasi_channels
Do i=1,max_iasi_channels
  iasi_channels(i) = i
End Do

!----- Acquire observation errors
IF (add_error .and. .not. (UsePCs .and. .not. Reconstruct) ) THEN
  Open(2,file=TRIM(noisefile))
  READ (2,'(A20)' ) text
  READ (2,*) Rtype, numRchannels, numRelements, Inverse 
  IF (numRchannels < n_chans) THEN
    write(*,*) 'Not enough elements in R matrix'
    STOP
  END IF
  Allocate(Rchannels(numRchannels))
  READ(2,*) Rchannels
  check=0
  Do i=1, n_chans
    chan=iasi_channels(i)
    If (ANY(iasi_channels == chan) ) check=check+1
  End Do
  IF (check /= n_chans) THEN
    write(*,*) 'Wrong elements in R matrix'
    STOP
  END IF
  Allocate(Rread(numRchannels))
  Read(2,*) Rread
  Close(2)
print*, Rread(1), Rread(8461)
  Allocate(Rdiag(n_chans))
  Rdiag=SQRT(Rread(iasi_channels))
  Deallocate(Rread)
END IF


!----- Get bias
IF (add_bias) THEN
  Open(2,file=biasfile)
  do i=1,183
    Read(2,*) biased_chans(i), bias(i)
  end do
  Close(2)
ENDIF

!------------------------------------------------------------------------------
!1.1) Set up RTTOV inputs
!------------------------------------------------------------------------------

ninst = 1
Allocate( instrument(3,ninst) )
instrument(1,1) = METOP_ID
instrument(2,1) = METOP_SAT_ID    ! MetOp-A
instrument(3,1) = IASI_ID
Allocate( coef(ninst) )
Allocate( errorstatus(ninst) )
nprofiles = 1
Allocate (options(ninst))
!1.4) Set up RTTOV11_Opts structure
!----

! Initialise options to default values 

options(:) % rt_all % plane_parallel    = .false.
options(:) % rt_ir % ir_sea_emis_model  = 2 
options(:) % rt_ir % ir_scatt_model     = 2
options(:) % rt_ir % vis_scatt_model    = 1
options(:) % rt_ir % dom_nstreams       = 8
options(:) % rt_ir % dom_accuracy       = 0.
options(:) % rt_ir % dom_opdep_threshold = 0.

options(:) % config % do_checkinput     = .true.
options(:) % config % apply_reg_limits  = .true.
options(:) % config % verbose           = Verbosity_level

IF (produceradiances) THEN
  options(:) % rt_all % switchrad       = .false.
ELSE
  options(:) % rt_all % switchrad       = .true.
END IF
options(:) % rt_all % addrefrac         = .false.
options(:) % rt_all % use_q2m           = .true.
options(:) % rt_all % do_lambertian     = .false.

options(:) % rt_mw % fastem_version     = 6 
options(:) % rt_mw % supply_foam_fraction = .false.
options(:) % rt_mw % clw_data         = .false.

options(:) % rt_ir % addsolar           = .false.
options(:) % rt_ir % do_nlte_correction = .false.
options(:) % rt_ir % addaerosl          = .false.
options(:) % rt_ir % addclouds          = .false.
options(:) % rt_ir % user_aer_opt_param = .false.
options(:) % rt_ir % user_cld_opt_param = .false.
options(:) % rt_ir % cldstr_threshold   = 0.001_jprb
options(:) % rt_ir % ozone_data         = .true.
options(:) % rt_ir % co2_data           = .false.
options(:) % rt_ir % n2o_data           = .false.
options(:) % rt_ir % co_data            = .false.
options(:) % rt_ir % ch4_data           = .false.

options(:) % interpolation % addinterp  = .false.
options(:) % interpolation % interp_mode = 5 
options(:) % interpolation % reg_limit_extrap = .true.
options(:) % interpolation % lgradp     = .false.
options(:) % interpolation % spacetop   = .true.

!Set up PC forward modelling
options(:) % rt_ir % pc % addpc         = .false.
options(:) % rt_ir % pc % addradrec     = reconstruct
options(:) % rt_ir % pc % ipcbnd        = 1_jpim
options(:) % rt_ir % pc % ipcreg        = 1_jpim

!Remember we are assuming only one instrument for PC assimilation
IF ( UsePCs ) THEN
  options(1) % rt_ir % pc % addpc  = .true.
  options(1) % rt_ir % pc % ipcreg = ipcreg
  options(1) % rt_ir % pc % ipcbnd = 1_jpim
  options(1) % rt_all % addrefrac  = .true.
  IF (.not. options(1) % rt_ir % pc % addradrec) THEN
    !Just to do this to avoid an error being generated in
    !rttov_user_options_checkinput
    options(:) % rt_all % switchrad         = .false.
  ENDIF
END IF


!---- Channel selection

IF ( UsePCs ) THEN
  NumInstChans =max_iasi_channels
  Allocate( channels(numinstchans,ninst) )
  channels(:,1) = 0
ELSE
  NumInstChans = n_chans
  Allocate( channels(n_chans,ninst) )
  channels(1:numinstchans,1)=iasi_channels(1:n_chans)
END IF



!3.2) Allocate RTTOV structures and initialise (mostly to zero)
!----


IF ( coef(1)%coef%nlevels /= NWPLevels) THEN
  options(:) % interpolation % addinterp = .TRUE.
  options(:) % config % apply_reg_limits = .TRUE.
  options(:) % interpolation % lgradp = .TRUE.
END IF


IF (UsePCs) THEN
  CALL RTTOV_READ_COEFS(         &
    ErrorStatus(1),                 & ! out
    coef(1),                     & ! out
    options(1),                  & ! in
    instrument = instrument(:,1) ) ! in
ELSE  
  CALL RTTOV_READ_COEFS(         &
    ErrorStatus(1),                 & ! out
    coef(1),                     & ! out
    options(1),                  & ! in
    channels = channels(:,1),    & ! in
    instrument = instrument(:,1) ) ! in
END IF
If(errorstatus(1) /= 0) Then
  Write(*,*) 'Error detected in RTTOV Read Coeffs!'
  Stop
Endif

Call RTTOV_USER_OPTIONS_CHECKINPUT( &
  ErrorStatus(1), & ! out
  options(1),  & ! in
  coef(1)  ) ! in

If(errorstatus(1) /= 0) Then
  Write(*,*) 'Error in User Options Checkinput!'
  Stop
Endif

!  Allocate(Rread(8461))
!  Do i=1,8461 
!    Rread(i)=coef(1) % coef_pccomp % noise_in(i)*coef(1) % coef_pccomp % noise_in(i)
!  End DO
!  Open(3,file='/home/h03/frfh/iasi/channels/instnoise_PCRTTOV_Rad_variance_rt12.txt',status='REPLACE')
!  Write(3,*) 'IASI on Metop1'
!  Write(3,*) '2  8461     1     0'
!  Write(3, '(10I8)') (/ (i,i=1,8461) /)
!  write(3,'(6F13.6)') Rread
!  CLOSE(3)
!  Deallocate(Rread)


!Reset channel selection for PCs
IF ( UsePCs ) THEN
  NumInstChans = SIZE(coef(1) % coef_pccomp % pcreg( &
                        options(1) % rt_ir % pc % ipcbnd, &
                        options(1) % rt_ir % pc % ipcreg ) % predictindex )
  Deallocate( channels )
  Allocate( channels(numinstchans,ninst) )
  channels(1:numinstchans,1) = coef(1) % coef_pccomp % pcreg( &
                        options(1) % rt_ir % pc % ipcbnd, &
                        options(1) % rt_ir % pc % ipcreg ) % predictindex(:)
END IF

ALLOCATE(chanprof(NumInstChans))
chanprof(1:NumInstChans)%chan =channels(1:NumInstChans,1)
chanprof(1:NumInstChans)%prof = 1_JPIM




!----- Open 70-level dataset

Open(1,file=profilefile,status='OLD')

!----- Open output files
Open(3,file=TRIM(outfile),status='REPLACE')


!----- Write obs file header
!----- Start with 10 lines for comments
Write(3,'(a)') 'This is a simulated clear-air IASI observation dataset.'
Write(3,'(a)') 'Generated by sim_rad_iasi_rttov12_70L.F90, using RTTOV on 101 levels.'
Write(3,'(a)') 'Truth file Model_20110315_1200_Sea.dat on 70L'
Write(3,'(a)') ''
!!!!!!!!!!!!!!!!!!!!!!!
Write(3,'(a)') 'Instrument noise added from instnoise_PCRTTOV.txt'
!!!!!!!!!!!!!!!!!!!!!!!!
Write(3,'(a)') ''
Write(3,'(a)') ''
Write(3,'(a)') ''
Write(3,'(a)') ''
Write(3,'(a)') ''
Write(3,'(a,i6)') 'Number of Observations in File:',idim
IF ( UsePCs .and. .not. reconstruct ) THEN
  Write(3,'(a,i8)') 'No. of Chans per Observation:',NPCScores
ELSE IF (reconstruct) THEN
  Write(3,'(a,i8)') 'No. of Chans per Observation:',max_iasi_channels
ELSE
  Write(3,'(a,i8)') 'No. of Chans per Observation:',numinstchans
ENDIF
Write(3,'(a)') 'Number of instruments making up observations : 1'
IF ( UsePCs .and. .not. reconstruct ) THEN
  Write(3,'(a)') '*** In the following Series, Platform and Instrument are defined  ***'
  Write(3,'(a)') '*** according to the relevant RT Model definitions (if required): ***'
  Write(3,'(a)') 'Units: PC Score'
  Write(3,'(a)') 'Sat. Series   Platform   Instrument First_Channel   Last_Channel  Sat ID'
  Write(3,'(i5,i13,i12,i12,i16,i12)') METOP_ID,METOP_SAT_ID,IASI_ID,1,NPCScores,sat_id
  Write(3,'(a)') 'Channels'
  Write(3,'(16i5)') (/ (i,i=1,NPCscores) /)
ELSE IF (reconstruct) THEN
  Write(3,'(a)') '*** In the following Series, Platform and Instrument are defined  ***'
  Write(3,'(a)') '*** according to the relevant RT Model definitions (if required): ***'
  IF (produceradiances) THEN
    Write(3,'(a)') 'Units: Radiance'
  ELSE
    Write(3,'(a)') 'Units: BT'
  END IF
  Write(3,'(a)') 'Sat. Series   Platform   Instrument First_Channel   Last_Channel  Sat ID'
  Write(3,'(i5,i13,i12,i12,i16,i12)') METOP_ID,METOP_SAT_ID,IASI_ID,1,max_iasi_channels,sat_id
  Write(3,'(a)') 'Channels'
  Write(3,'(16i5)') iasi_channels
ELSE
  Write(3,'(a)') '*** In the following Series, Platform and Instrument are defined  ***'
  Write(3,'(a)') '*** according to the relevant RT Model definitions (if required): ***'
  IF (produceradiances) THEN
    Write(3,'(a)') 'Units: Radiance'
  ELSE
    Write(3,'(a)') 'Units: BT'
  END IF
  Write(3,'(a)') 'Sat. Series   Platform   Instrument First_Channel   Last_Channel  Sat ID'
  Write(3,'(i5,i13,i12,i12,i16,i12)') METOP_ID,METOP_SAT_ID,IASI_ID,1,numinstchans,sat_id
  Write(3,'(a)') 'Channels'
  Write(3,'(16i5)') channels(1:numinstchans,1)
END IF
Write(3,'(a)') '----------------------------------------------------------------------'


!----- Read input file header and get units for humidity
do i=1,13
  read(1,'(A80)') text
end do

ii=Index(text,':')
jj=Index(text,'(')
string1=text(ii+1:jj-1)
read(string1,*) qunits

!----- Allocate profile structure
Allocate( profiles(nprofiles) )
CALL rttov_alloc_prof ( &
     Errorstatus(1),     &
     nprofiles,     &
     profiles,      &
     NWPLevels,  &
     options(1),    &
     asw=ASW_ALLOCATE,  &
     coefs=coef(1),      &
     init=.true.    )

!-----------------------------------------------
!2) Loop over observations and do the processing
!-----------------------------------------------

do i=1,idim

  !2.1) Read single profile
  !----

  !----- Read profile

  read(1,'(A80)') text
  read(1,'(A80)') text

  ii=Index(text,'Lat')
  jj=Index(text,'Lon')
  string1=text(ii+3:jj-1)
  string2=text(jj+3:)
  read(string1,*) lat
  read(string2,*) lon

  !   write(*,'(A,I6,A,F8.3,A,F8.3)')'Profile no:',i,' - Lon:',lon,' Lat:',lat

  read(1,'(A80)') text

  do j=1,nwplevels
    read(1,*)  pprof(j), tprof(j),qprof(j),oprof(j)
  end do

  IF (pprof(1) > pprof(2) ) THEN
    do j=1,nwplevels
      profiles(1)%p(j)=pprof(nwplevels-j+1)
      profiles(1)%t(j)=tprof(nwplevels-j+1)
      profiles(1)%q(j)=qprof(nwplevels-j+1)
      profiles(1)%o3(j)=oprof(nwplevels-j+1)
    end do
  ELSE
    profiles(1)%p(:)=pprof(:)
    profiles(1)%t(:)=tprof(:)
    profiles(1)%q(:)=qprof(:)
    profiles(1)%o3(:)=oprof(:)
  END IF
  IF ( ANY (profiles(1)%o3(:) <= 0.0) ) THEN 
    options(:) % rt_ir % ozone_data         = .false.
  END IF

  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%s2m%t
  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%s2m%q
  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%skin%t
  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%s2m%p
  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%s2m%u
  read(1,'(A80)') text
  ii=Index(text,':')
  text=text(ii+1:)
  read(text,*) profiles(1)%s2m%v

  !convert pressure units if required
  If (profiles(1)%s2m%p > 2000.0 ) THEN 
      profiles(1)%s2m%p=profiles(1)%s2m%p/100.0
      profiles(1)%p(:)=profiles(1)%p(:)/100.0
  END IF

  ! 2.2) Convert humidity units
  !-----

  IF (qunits==2) THEN
    qppmv(:)=0.0

    do jlev = 1,nwplevels
      qppmv(jlev)=profiles(1)%q(jlev)* H2O_MassMixToPPMV
    enddo
    profiles(1)%q(:)=qppmv(:)
    s2qppmv=profiles(1)%s2m%q * H2O_MassMixToPPMV
    profiles(1)%s2m%q          = s2qppmv
  END IF

  ! 2.3) Assign surface type and emissivity
  !-----

  profiles(1)%skin%surftype = surftype


  Surf_Emiss(1:NumInstChans)%emis_in = 0.
  Surf_Emiss(1:NumInstChans)%emis_out = 0.
  Calcemis(1:NumInstChans) = .TRUE.


  !2.4) Fill profile
  !----
  profiles(1)%nlevels        = nwplevels
  profiles(1)%nlayers        = NWPLayers
  profiles(1)%idg            = 1 ! Ice scheme - not used
  profiles(1)%ice_scheme     = 1 ! Crystal shape - not used
  !profiles(1)%clw(:)         = sclw(:) ! CLW only needed for microwave
  profiles(1)%skin%fastem(:) = 0.0 ! Initialise FASTEM coefs to zero
  Profiles(1)%skin%watertype = 1 ! Default to ocean water for now
  profiles(1) % skin % salinity = 0 !should we set this?
  profiles(1) % skin % foam_fraction = 0.0
  profiles(1) % skin % snow_fraction = 0.0
  profiles(1)%s2m%o          = 0.0
  profiles(1)%s2m%wfetc      = 100000.0 ! Wind fetch
  profiles(1)%zenangle       = 0.
  profiles(1)%azangle        = 0.
  Profiles(1)%sunzenangle    = 0.
  Profiles(1)%sunazangle     = 0.
  profiles(1)%ctp            = 500.0
  profiles(1)%cfraction      = 0.
  Profiles(1)% Be            = 0    ! zeeman - must initialise
  Profiles(1)% cosbk         = 0    ! zeeman - must initialise
  profiles(1) % gas_Units    = 2 ! ppmv !*****REMEMBER TO SET TO 2!!!!! *****
  profiles(1)%skin%fastem(:) = 0.0
  Profiles(1)%date           = 0
  Profiles(1)%elevation      = 0.
  Profiles(1)%latitude       = lat
  Profiles(1)%longitude      = lon
  profiles(1) % mmr_cldaer = .false. !units of aerosol input
  profiles(1) % s2m % wfetc = 100000.0 ! Wind fetch
  profiles(1) % s2m % o = 0.0 ! Not even used inside rttov!


!2.5) Call RTTOV
!----

!CALL RTTOV_PRINT_PROFILE(profiles(1))
!CALL RTTOV_PRINT_OPTS(options(1))

  !----- Allocate radiance
  Call RTTOV_ALLOC_RAD( &
      errorstatus(1),  & ! out
      n_Chans,         & ! in
      Radiance,        & ! in/out
      nwplevels,       & ! in
      ASW_ALLOCATE,    &
      init=.true.      ) ! in

  !----- Allocate transmission
  CALL rttov_alloc_transmission( &
      errorstatus(1),           &
      transmission,             &
      nwplevels,                &
      n_chans,                  &
      ASW_ALLOCATE,             &
      init=.true.)


  !----- Allocate PC Component structure
  IF ( UsePCs ) THEN

    IF (reconstruct) THEN
      CALL rttov_alloc_pccomp( Errorstatus(1), &
                              PCcomp, &
                              NPCScores, &
                              ASW_ALLOCATE, &
                              init = .TRUE., &
                              nchannels_rec=max_iasi_channels)
    ELSE
      CALL rttov_alloc_pccomp( Errorstatus(1), &
                              PCcomp, &
                              NPCScores, &
                              ASW_ALLOCATE, &
                              init = .TRUE.)
    ENDIF
  END IF

    IF ( UsePCs ) THEN
      IF ( reconstruct ) THEN
        CALL RTTOV_DIRECT(Errorstatus(1),                   & ! out
                          chanprof,                      & ! in
                          options(1),                    & ! in
                          Profiles,                      & ! in
                          coef(1),                       & ! in
                          transmission,                  & ! inout
                          Radiance,                      & ! inout
                          calcemis=calcemis,                     & ! in
                          emissivity=Surf_Emiss(1:NumInstChans),  & ! in
                          pccomp = PCcomp,               & ! PC scores
                          channels_rec=iasi_channels     )
      ELSE
        CALL RTTOV_DIRECT(Errorstatus(1),                   & ! out
                          chanprof,                      & ! in
                          options(1),                    & ! in
                          Profiles,                      & ! in
                          coef(1),                       & ! in
                          transmission,                  & ! inout
                          Radiance,                      & ! inout
                          calcemis=calcemis,                     & ! in
                          emissivity=Surf_Emiss(1:NumInstChans),  & ! in
                          pccomp = PCcomp                ) ! PC scores

      ENDIF
    ELSE
      CALL RTTOV_DIRECT(Errorstatus(1),                   & ! out
                          chanprof,                      & ! in
                          options(1),                    & ! in
                          Profiles,                      & ! in
                          coef(1),                       & ! in
                          transmission,                  & ! inout
                          Radiance,                      & ! inout
                          calcemis=calcemis,                     & ! in
                          emissivity=Surf_Emiss(1:NumInstChans))! in
    END IF

    If (any(errorstatus(:) /= 0)) Then
      Write(*,*) 'Error detected in RTTOV!'
    Endif


  !2.6) Add simulated errors to BTs
  !----

  nchans_out=Numinstchans
  !Transfer reconstructed radiances to radiance array
  IF (UsePCs .and. reconstruct) THEN
    IF (produceradiances) THEN
      Deallocate(radiance%total)
      ALLOCATE(radiance%total(max_iasi_channels))
      radiance%total=pccomp%total_pccomp
    ELSE
      Deallocate(radiance%bt)
      ALLOCATE(radiance%bt(max_iasi_channels ))
      radiance%bt=pccomp%bt_pccomp
    ENDIF
    nchans_out=max_iasi_channels
    deallocate(channels)
    allocate(channels(max_iasi_channels,1))
    channels(:,1)=iasi_channels
  ENDIF


  IF (add_error) THEN
    IF ( UsePCs .and. .not. reconstruct ) THEN
      do j=1,npcscores
        Call Random_Gaussian(rndg)
        pccomp%PCscores(j)=pccomp%PCscores(j) + rndg
      end do
    ELSE IF ( produceradiances ) THEN
      do j=1,nchans_out
        Call Random_Gaussian(rndg)
        radiance%total(j) = radiance%total(j) + rndg*Rdiag(channels(j,1))
      end do
!         IF ( ANY( radiance%total(:) < 0.0 ) ) &
!           write(*,*) 'negative radiance for observation', i
    ELSE
      do j=1,nchans_out
        Call Random_Gaussian(rndg)
        IF ( j == 8 ) print*, rndg*Rdiag(channels(j,1))
        radiance%bt(j) = radiance%bt(j) + rndg*Rdiag(channels(j,1))
      end do
    ENDIF
  END IF

  IF (add_bias) THEN
    radiance%bt(biased_chans) = radiance%bt(biased_chans)-bias
  END IF

!2.7) Write profile information to file
!----
  if(surftype.eq.0) surftype_1dvar=3
  if(surftype.eq.1) surftype_1dvar=1
  if(surftype.eq.2) surftype_1dvar=2
  write(3,'(a,i15,a,i11,a,i6)') &
    'Obs ID:',i,' Obs Type:',3,' Satellite ID:',sat_id
  write(3,'(a,f9.3,a,f9.3,a,f7.1)') &
    'Latitude:',lat,' Longitude:',lon,' Elevation:', 0.0
  write(3,'(a,i4,a,f9.3,a,f9.3)') &
    'Surface Type:',surftype_1dvar, &
    ' Sat Zen Angle:',profiles(1)%zenangle,' Solar Zen. Ang.:',0.0
  IF ( UsePCs .and. .not. reconstruct ) THEN
    write(3,'(a)') 'PC Scores:'
    write(3,'(6e13.5)') pccomp%PCscores(1:NPCscores)
  ELSEIF ( ProduceRadiances ) THEN
    write(3,'(a)') 'Radiances:'
    write(3,'(6f13.3)') radiance%total(1:NChans_out)
  ELSE
    write(3,'(a)') 'Brightness Temperatures:'
    write(3,'(6f13.3)') radiance%bt(1:NChans_out)
  END IF

  Call RTTOV_ALLOC_RAD(&
      ErrorStatus(1), & ! out
      N_Chans,        & ! in
      Radiance,       & ! in/out
      NwpLevels,      & ! in
      ASW_DEALLOCATE  ) ! in


  CALL rttov_alloc_transmission( &
      ErrorStatus(1),&
      transmission,  &
      NwpLevels,     &
      N_Chans,       &
      ASW_DEALLOCATE )



  !----- Allocate PC Component structure
  IF ( UsePCs ) THEN
    CALL rttov_alloc_pccomp( Errorstatus(1), &
                              PCcomp, &
                              NPCScores, &
                              ASW_DEALLOCATE)
  END IF

end do

!2.8) Close files
!----

Close(1)
!2 already closed
Close(3)
Close(4)


!2.9) Tidy up
!----

Call RTTOV_DEALLOC_COEFS( &
     ErrorStatus(1),     & ! out
     Coef(1)             ) ! out

IF ( UsePCs ) THEN
  CALL rttov_alloc_pccomp( ErrorStatus(1), &
                           PCcomp, &
                           NPCScores, &
                           ASW_DEALLOCATE )
END IF


!Deallocate profile structure
!----- Allocate profile structure
CALL rttov_alloc_prof ( &
     ErrorStatus(1),    &
     nprofiles,     &
     profiles,      &
     NWPLevels,  &
     options(1),    &
     ASW_DEALLOCATE,  &
     coef(1)      )
Deallocate(profiles)

!Deallocate other stuff
Deallocate(channels)
Deallocate(chanprof)
Deallocate(instrument)
Deallocate(coef)
Deallocate(errorstatus)

CONTAINS

SUBROUTINE Random_Gaussian(RandG)
! Generates random Gaussian number with std dev=1
  Real, Intent(Out)  :: RandG
  Integer, parameter :: nrand=12
  Real               :: harvest(nrand)
  Call Random_Number(harvest)
  harvest(:) = harvest(:)
  RandG = sqrt(12.0/Real(nrand))*Sum(harvest)-sqrt(3.0*Real(nrand))
END SUBROUTINE Random_Gaussian

END PROGRAM sim_spec


