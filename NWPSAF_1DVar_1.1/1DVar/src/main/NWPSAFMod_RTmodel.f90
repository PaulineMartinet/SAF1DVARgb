!+ routines and constants for RTModel

MODULE NWPSAFMod_RTmodel

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
! Description: contains information for RT Models
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     21/10/99 Original module.      A.D. Collard.  Met Office.
! 1.1     01/03/01 Added RTTOV           R.W. Saunders
! 2.3     22/05/02 Changes to allow more flexibility on setting up 
!                  multi-instruments such as ATOVS.          ADC.
! 3.0.1   15/07/03 Add cloudy stuff (M. Szyndel)             ADC.
! 3.0.3   26/02/04 Add Gastropod stuff (V. Sherlock) and RT Model
!                  names.                                  A. Collard.
! 3.0.4   02/03/04 Changes in order to generalise adding new instruments.
!                                                          A. Collard.
! 3.0.5   29/03/04 Remove all references to orbit height and SatViewAngle.  
!                                                          A. Collard.
! 3.0.6   18/06/04 Add RT1stGuess to RTParams_Type.        A. Collard.
! 3.0.7   05/04/05 Added fastmodel identity for RTTOV8.    E. Pavelin.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date     Comment
! 13      29/01/09 Added fastmodel identity for RTTOV9.    E. Pavelin.
! 25      20/02/12 Added fastmodel identity for RTTOV10 with
!                  longitude, date and UseModelLevels      P. Weston.
! 28      22/02/12 Added Prof_FirstCLW and Prof_LastCLW    TR Sreerekha
! 31      18/01/13 Added fastmodel identity for RTTOV11    P. Weston.
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------
USE rttov_types, ONLY : &
    rttov_coefs,        &
    rttov_options,      &
    rttov_chanprof

#ifdef _RTTOV12
USE mod_rttov_emis_atlas, ONLY : &
    rttov_emis_atlas_data
#endif

IMPLICIT NONE

!--------------------------------------------
! 1. RTTOV structures required throughout code
!---------------------------------------------

! Coefficient structure for RTTOV
TYPE(rttov_coefs), POINTER :: RT_Coefs(:)

! Options list for RTTOV
TYPE(rttov_options), POINTER :: RT_Opts(:)

! Channels and profiles structure for RTTOV
TYPE(rttov_chanprof), POINTER :: RT_Chanprof(:)

! Emissivity atlas structure for RTTOV12 only
#ifdef _RTTOV12
TYPE(rttov_emis_atlas_data), SAVE :: emis_atlas
#endif

!-----------------------------------------------
!1. Setup fastmodel control parameters
!-----------------------------------------------

LOGICAL, SAVE :: UseModelLevels = .FALSE.
Integer, PARAMETER :: MaxChans   = 10000    ! Max. no. of channels
Integer, PARAMETER :: MaxSensors = 30       ! Max no sensors to be used

!-----------------------------------------------
!2. Types of fastmodel calls required
!-----------------------------------------------

INTEGER, PARAMETER :: &
     FastmodelMode_Initialise = 0, &
     FastmodelMode_Forward    = 1, &
     FastmodelMode_Gradient   = 2, &
     FastmodelMode_CleanUp    = 3

INTEGER, PARAMETER :: &
     NeitherProfile = 0, &
     BackGrProf     = 1, &
     GuessProf      = 2

!----------------------------------------------------------------
! 3. Positions of various parameters in the RT model vector.
!----------------------------------------------------------------

INTEGER, SAVE :: ProfSize

INTEGER, SAVE :: Prof_FirstT
INTEGER, SAVE :: Prof_FirstQ
INTEGER, SAVE :: Prof_FirstO3 
INTEGER, SAVE :: Prof_FirstCLW
INTEGER, SAVE :: Prof_LastT
INTEGER, SAVE :: Prof_LastQ
INTEGER, SAVE :: Prof_LastO3 
INTEGER, SAVE :: Prof_LastCLW
INTEGER, SAVE :: Prof_T2 
INTEGER, SAVE :: Prof_q2 
INTEGER, SAVE :: Prof_Tstar
INTEGER, SAVE :: Prof_pstar
INTEGER, SAVE :: Prof_uwind
INTEGER, SAVE :: Prof_vwind
INTEGER, SAVE :: Prof_CTP
INTEGER, SAVE :: Prof_CloudCover 
INTEGER, SAVE :: Prof_LWP

INTEGER, SAVE :: Num_ProfElementsUsed ! Number of elements to be used in
                                      ! retrievals
INTEGER, SAVE :: Num_WetLevels  ! Number of humidity levels - used for CLW
                                ! retrievals

INTEGER, SAVE :: NumPredChansPC ! Number of Predictor Channels for PC_RTTOV

!-----------------------------------------------
!4. Setup fastmodel parameters for RTTOV
!-----------------------------------------------

! Set up some RT parameters

! Maximum number of satellite series

INTEGER, PARAMETER :: Max_Series = 13

! Maximum number of platforms per series

INTEGER, PARAMETER :: Max_Platforms = 20 

! Maximum number of instruments per platform

INTEGER, PARAMETER :: Max_SubTypes = 26

! Input profile information
INTEGER, SAVE :: Num_RTLevels      
INTEGER, SAVE :: Num_Profs = 1    


TYPE Limit_Type
   REAL, POINTER :: Minimum(:)
   REAL, POINTER :: Maximum(:)
ENDTYPE Limit_Type
TYPE(Limit_Type), SAVE :: Soft_Limits


!  The following structure is used to communicate radiative transfer data
!  to the RT model.  They contain variables that may possibly change between
!  observations (including a possible change of satellite).  This is currently
!  set up for RTTOV7 and needs slight modifications to work for RTTOV6 and 
!  earlier.

!
TYPE SatIDInfo_Type
   CHARACTER (LEN=20) :: SatID_Text
   INTEGER :: WMO         ! WMO Number 
   INTEGER :: First_Instr
   INTEGER :: Last_Instr
   INTEGER :: Num_Channels
   LOGICAL :: R_Matrix_Present
ENDTYPE SatIDInfo_Type

TYPE RTParams_Type
   
   ! Choices of satellites, instruments, etc. 
   INTEGER :: SatIndex     
   INTEGER, POINTER :: SeriesChoice(:)       
   INTEGER, POINTER :: PlatformChoice(:)        
   INTEGER, POINTER :: SubTypeChoice(:)      
   INTEGER, POINTER :: NumChans(:)      
   INTEGER, POINTER :: RTCoeffs(:)
   !PM
   REAL, POINTER :: view_angle(:)
   !PM
   
   ! In the observation vector, which channel belongs to which instrument?
   INTEGER :: Num_Instruments
   INTEGER, POINTER :: Instrument_Number(:,:)
   INTEGER, POINTER :: Absolute_Channel_Number(:)
   INTEGER, POINTER :: First_Channel_for_Instrument(:)

   ! SatID stuff.  This indicates component instruments of combined 
   ! instruments such as ATOVS.
   INTEGER :: Num_SatIDs
   TYPE (SatIDInfo_Type), POINTER :: SatID(:)

   ! Input profile information
   REAL, POINTER :: RTBack(:)
   REAL, POINTER :: RT1stGuess(:)
   REAL, POINTER :: RTGuess(:)
   REAL, POINTER :: RTEmissivity(:)
   INTEGER :: RTSurfaceType     

   ! Temporal information
   INTEGER :: Date(3)  ! Year, Month, Day

   ! Geometrical information
   REAL :: SatZenithAngle
   REAL :: SolarZenAngle !solar zenith angle
   REAL :: SatAzimAngle
   REAL :: SatSolarAzimAngle
   REAL :: Latitude
   REAL :: Longitude
   REAL :: Elevation

   ! Output BT/Radiance/Gradient information
   REAL, POINTER :: TotalBTs(:)
   REAL, POINTER :: TotalRadiances(:)
   REAL, POINTER :: PCScores(:)
   REAL, POINTER :: CloudyRadiances(:,:)
   REAL, POINTER :: H_Matrix(:,:)
   REAL, POINTER :: H_Matrix_T(:,:)

   ! Channel constants (used in NWPSAF_CO2Slice)
  
   REAL, POINTER :: tc1(:)
   REAL, POINTER :: tc2(:)
   REAL, POINTER :: bcon1(:)
   REAL, POINTER :: bcon2(:)

   ! Pressure levels for current observation in hPa (for QSAT)
   REAL, POINTER :: Pressure_Pa(:)

END TYPE RTParams_Type




END MODULE NWPSAFMod_RTmodel



