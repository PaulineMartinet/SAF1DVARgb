!+ routines and constants for things to do with LWP Retrieval

MODULE NWPSAFMod_LiquidWater

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
! Ticket  Date     Comment
!         07/10/16 New, from varions Layer to LWP Conversion subroutines. FSmith
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

IMPLICIT NONE


!----------------------------------------------------------------------
! Saturation Vapour Pressure Calculations
!----------------------------------------------------------------------

! Teten's formula parameters that are needed to compute
!  saturation water vapor pressure.:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.
REAL,    PARAMETER  ::  awater = 17.269
REAL,    PARAMETER  ::  bwater = 35.86
REAL,    PARAMETER  ::  aice   = 21.875
REAL,    PARAMETER  ::  bice   = 7.66
REAL,    PARAMETER  ::  Tmelt  = 273.16
REAL,    PARAMETER  ::  Cteten = 610.78
LOGICAL, PARAMETER  ::  bothphases=.TRUE.


CONTAINS

FUNCTION SVP (T)
 
  IMPLICIT NONE
  REAL, INTENT(In)  :: T    ! Temperature in Kelvin
  REAL              :: SVP  ! Pa       

  REAL            :: aconst
  REAL            :: bconst

!  Use Tetens' formula to calculate Saturated Vapour Pressure:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.
! Define constants to be used depending on phase (water or ice)

  aconst = awater
  bconst = bwater

  IF (bothphases.AND. T < Tmelt) THEN
        aconst = aice  
        bconst = bice  
  ENDIF
  SVP = Cteten * exp( aconst*(T-Tmelt)/(T-bconst) ) 

END FUNCTION SVP

FUNCTION SVP_Deriv (T)

IMPLICIT NONE

!Arguments
REAL, INTENT(IN) :: T         ! Temperature in Kelvin
REAL             :: SVP_deriv !Dln(es)_DT

REAL            :: aconst
REAL            :: bconst
REAL            :: aw           ! T-bconst

!  Use Tetens' formula:
!       Murray, F.W., 1966: On the computation of saturation
!       vapor pressure, J. Applied Meteorol., 6, 203-204.
!  Define constants to be used depending on phase (water or ice)

aconst = awater
bconst = bwater

IF (T < Tmelt.AND.bothphases) THEN
   aconst = aice
   bconst = bice
ENDIF

! Compute 

aw = T-bconst
svp_deriv = aconst*( Tmelt-bconst )/( aw*aw )

END FUNCTION SVP_Deriv

!----------------------------------------------------------------------------
! LWP to Layers and back
!----------------------------------------------------------------------

FUNCTION LWP_to_Layers(LWP_new,cloud_structure)

  IMPLICIT NONE

  REAL, INTENT(IN)  :: cloud_structure(:)
  REAL, INTENT(IN)  :: LWP_new

  REAL :: LWP_to_Layers(size(cloud_structure))

  LWP_to_Layers(:) = LWP_new * cloud_structure(:)

END FUNCTION LWP_to_Layers



FUNCTION Layers_to_LWP(Pressure,CLW)

  USE rttov_const, ONLY : &
    gravity

  IMPLICIT NONE

  REAL, INTENT(IN)  :: Pressure(:) ! Pressure in Pa on levels
  REAL, INTENT(IN)  :: clw(:)      ! cloud liquid water on the same levels

  INTEGER  :: NumLevs
  INTEGER  :: i
  REAL     :: Layers_to_LWP

  NumLevs=Size(clw(:))
  ! Not sure why we integrate from level 2?
  ! This seems to be at odds with the LWP_to_Layers code...
  Layers_to_LWP=0.0
  DO i=2,NumLevs
     Layers_to_LWP = Layers_to_LWP + (Pressure(i)-Pressure(i-1))*(clw(i)+clw(i-1))
  ENDDO
  
  Layers_to_LWP=(0.5*Layers_to_LWP)/gravity

END FUNCTION Layers_to_LWP


FUNCTION LayerK_to_LWPK (clw_k, cloud_structure)

  REAL, INTENT(IN)  :: clw_k(:)           ! cloud liquid water jacobian 
  REAL, INTENT(IN)  :: cloud_structure(:) ! on the same levels as above

  INTEGER  :: NumLevs
  INTEGER  :: i
  REAL     :: LayerK_to_LWPK

  NumLevs=Size(clw_k(:))
  LayerK_to_LWPK=0.0
  ! Not sure why we integrate from level 2?
  ! This seems to be at odds with the LWP_to_Layers code...
  DO i=2, NumLevs    
     LayerK_to_LWPK = LayerK_to_LWPK +  0.5*(clw_k(i)+clw_k(i-1))*Cloud_Structure(i)
  ENDDO

END FUNCTION LayerK_to_LWPK

END MODULE NWPSAFMod_LiquidWater



