!+ Data types and information 

MODULE NWPSAFMod_ObsInfo

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
!  This module currently contains:
!
!  a) Derived data types for Observations
!
!  b) Derived data types for Model equivalent of Observations
!
! History:
! 
! Version   Date     Comment
! -------   ----     -------
! 1.0       04/01/01 Pared down version of original OPSMod_ObsInfo to 
!                    be used purely for standalone NWPSAF_1DVar.
!                                            A. Collard  Met Office
! 1.1       02/03/01 Modified ozone R. Saunders
! 2.3       27/05/02 Removed Status from ElementHeaderType.  Andrew Collard.
! 3.0.1     15/07/03 Added Cloud parameters (from M. Szyndel)  A. Collard.
! 3.0.5     29/03/04 SatView becomes SatZenith               A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
! 25        20/02/12 Added date to ob structure  P.Weston.
! 28        22/02/12 Added LWP and clw to the ob structure and 
!           22/02/12 model_ob structures  TR.Sreerekha
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: GTDP8
!
! End of header -------------------------------------------------------

IMPLICIT NONE

! Global TYPE Definitions:

!-------------------------------------------------------------------   
! Element Header data type used in OBHeader, CXHeader and ModelOBHeader 
!-------------------------------------------------------------------   
 TYPE ElementHeader_type 

   INTEGER :: NumLev   ! Number of levels of data

 End type ElementHeader_type 

!-------------------------------------------------------------------   
!  Element data type used in OB_type
!-----------------------------------
  TYPE Element_type
    REAL :: Value  
  END TYPE Element_type

!-------------------------------------------------------------------   
!  Coord data type used in OB_type
!----------------------------------
  TYPE coord_type
    REAL :: Value  
 ! The current code does not make any allowance for possible errors in the
 ! observation position, so no space is provided to characterise such errors.
 ! A derived type is used so that a future version can add OBErr etc if desired.
  END TYPE coord_type

!-------------------------------------------------------------------   
! ModelOB Header data types used in ModelOB_type
!------------------------------------------------
! There is a one to one correspondence between the items 
! in ModelOBheader and ModelOB
! There is a one to one correspondence between the items 
! in ModelOBheader and OBheader except for those variables in ObHeader
! described as MetaData
!
  TYPE ModelOBheader_type
   INTEGER                   :: NumObsTotal  ! Number of Obs in total
   INTEGER                   :: NumObsLocal  ! Number of Obs locally
! Surface variables
   TYPE (ElementHeader_type) :: pstar        ! DPI for pstar
   TYPE (ElementHeader_type) :: t2           ! DPI for surface temperature
   TYPE (ElementHeader_type) :: rh2          ! DPI for surface humidity
   TYPE (ElementHeader_type) :: uv10         ! DPI for surface wind
   TYPE (ElementHeader_type) :: Tskin        ! DPI for Surface radiative temp
   TYPE (ElementHeader_type) :: CTP          ! DPI for Cloud Top Pressure
   TYPE (ElementHeader_type) :: Cldfrac      ! DPI for Cloud Fraction
   TYPE (ElementHeader_type) :: LWP          ! DPI for MW Liquid water path
   ! upper-level variables,             (NumObs,NumLev)
   ! single level uppper level variables are catered for here with NumLev=1
   TYPE (ElementHeader_type) :: p            ! DPI for pressure                
   TYPE (ElementHeader_type) :: t            ! DPI for temperature
   TYPE (ElementHeader_type) :: rh           ! DPI for rh
   TYPE (ElementHeader_type) :: ozone        ! DPI for ozone
   TYPE (ElementHeader_type) :: clw          ! DPI for clw profile
   TYPE (ElementHeader_type) :: BriTemp      ! DPI for brightness temps (K) 
   TYPE (ElementHeader_type) :: Radiance     ! DPI for radiances
   TYPE (ElementHeader_type) :: PCScore      ! DPI for PC scores
  END TYPE ModelOBheader_type

!-------------------------------------------------------------------   
! ModelOB data types
!--------------------
! There is a one to one correspondence between the items 
! in ModelOBheader and ModelOB
! There is a one to one correspondence between the items 
! in ModelOB and OB except for those variables in Ob
! described as MetaData
!
  TYPE ModelOB_type
    TYPE (ModelOBheader_type) :: header  
! Surface variables,                   (NumObs)
    REAL, POINTER             :: pstar     (:) ! pressure at model surface Pa
    REAL, POINTER             :: t2        (:) ! 2m temperature            K
    REAL, POINTER             :: rh2       (:) ! 2m relative humidity      %
    REAL, POINTER             :: u10       (:) ! 10m westerly wind         m/s
    REAL, POINTER             :: v10       (:) ! 10m southerly wind        m/s
    REAL, POINTER             :: Tskin     (:) ! Surface radiative temp. K
    REAL, POINTER             :: CTP       (:) ! Cloud Top Pressure
    REAL, POINTER             :: Cldfrac   (:) ! Cloud Fraction
    REAL, POINTER             :: LWP       (:) ! MW liquid water path
    ! upper-level variables,   (NumObs,NumLev) except for BriTemp which is (NumLev)
! single level upper level variables are catered for here with NumLev=1
    REAL, POINTER             :: p       (:,:) ! pressure                Pa 
    REAL, POINTER             :: t       (:,:) ! temperature              K 
    REAL, POINTER             :: rh      (:,:) ! relative humidity        % 
    REAL, POINTER             :: ozone   (:,:) ! ozone                 kg/Kg
    REAL, POINTER             :: clw     (:,:) ! clw profile           ppmv
    REAL, POINTER             :: BriTemp (:)   ! brightness temps         K
    REAL, POINTER             :: Radiance (:)  ! mW m-2 sr-1 cm
    REAL, POINTER             :: PCScore (:)   
  END TYPE ModelOB_type

!-------------------------------------------------------------------   
! OB header types used in OB_type
!-----------------------------------
! There is a one to one correspondence between the items 
! in OBheader and OB
!
  TYPE OBheader_type
   INTEGER                 :: NumObsTotal  ! Number of Obs on this PE
   INTEGER                 :: NumObsLocal  ! Number of Obs overall

! Start of MetaData
! The following three have been re-instated as header elements.  ADC 21/7/99
   TYPE (ElementHeader_type) :: Latitude     ! DPI for Latitude
   TYPE (ElementHeader_type) :: Longitude    ! DPI for Longitude
   TYPE (ElementHeader_type) :: ObsTYPE      ! DPI for ObsTYPE
   TYPE (ElementHeader_type) :: SatId        ! DPI for Satellite id
   TYPE (ElementHeader_type) :: Date         ! DPI for Date

! 1DVAR specific MetaData information (satellite sounding)
   TYPE (ElementHeader_type) :: ScanLine       ! DPI for scanline number
   TYPE (ElementHeader_type) :: ScanPosition   ! DPI for scan position
   TYPE (ElementHeader_type) :: surface        ! DPI for surface
   TYPE (ElementHeader_type) :: elevation      ! DPI for elevation
   TYPE (ElementHeader_type) :: LocalZenith    ! DPI for local zenith angle
   TYPE (ElementHeader_type) :: LocalAzimuth   ! DPI for local azimuth angle
   TYPE (ElementHeader_type) :: SatZenith      ! DPI for SatZenith 
   TYPE (ElementHeader_type) :: SatAzimth      ! DPI for SatAzimth
   TYPE (ElementHeader_type) :: SolarZenith    ! DPI for SolarZenith
   TYPE (ElementHeader_type) :: SolarAzimth    ! DPI for SolarAzimth  
   TYPE (ElementHeader_type) :: NIter          ! DPI for NIter
   TYPE (ElementHeader_type) :: CldCost        ! DPI for CldCost
   TYPE (ElementHeader_type) :: J              ! DPI for J (1dvar cost function)
   TYPE (ElementHeader_type) :: coloz          ! DPI for total col ozone   
   TYPE (ElementHeader_type) :: IREmiss        ! DPI for IR surf emissivities
   TYPE (ElementHeader_type) :: JCost          ! DPI for JCost
   TYPE (ElementHeader_type) :: JCost_Gradient ! DPI for JCost


! End of MetaData

! Surface variables
   TYPE (ElementHeader_type) :: pstar        ! DPI for pstar
   TYPE (ElementHeader_type) :: t2           ! DPI for t2
   TYPE (ElementHeader_type) :: rh2          ! DPI for rh2
   TYPE (ElementHeader_type) :: uv10         ! DPI for u & v 10 
   TYPE (ElementHeader_type) :: Tskin        ! DPI for Surface radiative temp
   TYPE (ElementHeader_type) :: CTP          ! DPI for Cloud Top Pressure
   TYPE (ElementHeader_type) :: Cldfrac      ! DPI for Cloud Fraction
   TYPE (ElementHeader_type) :: LWP          ! DPI for MW liquid water path
! upper-level variables,             (NumObs,NumLev)
! single level uppper level variables are catered for here with NumLev=1
   TYPE (ElementHeader_type) :: t            ! DPI for temperature
   TYPE (ElementHeader_type) :: rh           ! DPI for rh
   TYPE (ElementHeader_type) :: ozone        ! DPI for ozone
   TYPE (ElementHeader_type) :: clw          ! DPI for cloud liquid water
   TYPE (ElementHeader_type) :: BriTemp      ! DPI for brightness temps (K) 
   TYPE (ElementHeader_type) :: Radiance     ! DPI for radiances
   TYPE (ElementHeader_type) :: PCScore      ! DPI for PC scores
  END TYPE OBheader_type

!-------------------------------------------------------------------   
! OB type
!-----------------------------------
! There is a one to one correspondence between the items 
! in OBheader and OB
!

  TYPE OB_type
    TYPE (OBheader_type)     :: header       

! Start of MetaData
    TYPE (coord_type) :: Latitude  ! Degrees (-90 to 90) 
    TYPE (coord_type) :: Longitude ! Degrees (-180 to 180)
    INTEGER           :: ObsTYPE    
    INTEGER           :: Id        ! Unique id for each ob
    INTEGER           :: SatId     ! Satellite id 
    INTEGER           :: Date(3)   ! Year, Month, Day

! 1DVAR specific MetaData information (satellite sounding)

   INTEGER :: ScanLine       ! scanline number
   INTEGER :: ScanPosition   ! scan position
   INTEGER :: surface        ! type from the report
   REAL    :: elevation      ! elevation (m) from report
   REAL    :: LocalZenith    ! local zenith angle
   REAL    :: LocalAzimuth   ! local azimuth angle
   REAL    :: SatZenith      ! satellite zenith angle
   REAL    :: SatAzimth      ! satellite azimuth angle
   REAL    :: SolarZenith    ! solar zenith angle
   REAL    :: SolarAzimth    ! solar azimuth angle
   INTEGER :: NIter          ! number of iterations
   REAL    :: CldCost        ! 1DVAR cloud detection
   REAL    :: J              ! 1DVAR cost function
   REAL    :: coloz          ! total column ozone (Du)  
   REAL    :: IREmiss        ! infra-red surface emissivity
   REAL    :: JCost          ! JCost
   REAL    :: JCost_Gradient ! JCost_Gradient


! End of MetaData

! Surface variables,                     
    TYPE (Element_type) :: pstar   ! pressure at model surf  Pa
    TYPE (Element_type) :: t2      ! 2m temperature          K
    TYPE (Element_type) :: rh2     ! 2m relative humidity    %
    TYPE (Element_type) :: u10     ! 10m westerly wind       m/s
    TYPE (Element_type) :: v10     ! 10m southerly wind      m/s
    TYPE (Element_type) :: Tskin   ! Surface radiative temp. K
    TYPE (Element_type) :: CTP     ! Cloud Top Pressure
    TYPE (Element_type) :: Cldfrac ! Cloud Fraction
    TYPE (Element_type) :: LWP     ! MW liquid water path
! upper-level variables,             (NumLev)
! single level uppper level variables are catered for here with NumLev=1
    TYPE (Element_type), POINTER :: t     (:) ! temperature         K
    TYPE (Element_type), POINTER :: rh    (:) ! relative humidity   %
    TYPE (Element_type), POINTER :: ozone (:) ! ozone mmr           kg/Kg 
    TYPE (Element_type), POINTER :: clw   (:) ! cloud liquid water  kg/Kg 
    REAL, POINTER             :: BriTemp (:)   ! brightness temps         K
    REAL, POINTER             :: Radiance (:)  ! mW m-2 sr-1 cm
    REAL, POINTER             :: PCScore (:)        
  END TYPE OB_type

END MODULE NWPSAFMod_ObsInfo



