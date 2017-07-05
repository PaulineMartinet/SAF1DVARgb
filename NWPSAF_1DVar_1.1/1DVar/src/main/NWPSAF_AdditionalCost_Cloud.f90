SUBROUTINE NWPSAF_AdditionalCost_Cloud( &
     Guess_Profile,    &  ! in
     JAdd,             &  ! out
     JAdd_Gradient,    &  ! out
     JAdd_Gradient2     ) ! optional out

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
!
!   Calculates an additional cost function and gradient if the retrieved
!   cloud top pressure falls outside of the range 0-PStar or the Cloud
!   fraction falls outside of the range 0-1.
!
!   The form of the additional cost function is dJ=ScaleFactor*dx^3 where
!   dx is how far the value falls ouside of the bounds and the scale factor
!   is chosen so that dJ becomes large when dx is greater than a suitably 
!   small value.
!
!
! History:
!   3.2   04/02/2004  New routine to stop CTP and Cloud Fraction straying
!                     out of bounds.
!


USE NWPSAFMod_Constants, ONLY : &
    ScaleFactorP, &
    ScaleFactorF

USE NWPSAFMod_Params, ONLY :    &
    Ret_CTP,                  &
    Ret_CloudCover,           &
    Cloud_Min_Pressure,       &
    CloudFree,                &
    ThickCloud,               &
    GeneralMode,              &
    VerboseMode

USE NWPSAFMod_RTmodel, ONLY : &
    Prof_CloudCover, &
    Prof_CTP, &
    Prof_PStar

IMPLICIT NONE

! Subroutine arguments:
REAL, INTENT(IN)    :: Guess_Profile(:)     ! Guess Profile (on RT levels)
REAL, INTENT(OUT)   :: JAdd                 ! Cost function value
REAL, INTENT(OUT)   :: JAdd_Gradient(:)     ! Cost function gradient
REAL, INTENT(OUT), OPTIONAL :: JAdd_Gradient2(:,:)  ! Cost function 2nd Deriv.

! Local variables:

REAL :: CloudCover
REAL :: CTP
REAL :: PStar

!---------------------------------------------------------------------
! 0. Initialise
!---------------------------------------------------------------------

CloudCover = Guess_Profile(Prof_CloudCover)
CTP = Guess_Profile(Prof_CTP)
PStar = Guess_Profile(Prof_PStar)

JAdd             = 0.
JAdd_Gradient(:) = 0.
IF (PRESENT(JAdd_Gradient2)) JAdd_Gradient2(:,:) = 0.

!---------------------------------------------------------------------
! 1. Cloud Top Pressure Additional Cost Function
!---------------------------------------------------------------------

IF ( CTP > PStar ) THEN
  JAdd = JAdd + ScaleFactorP * (CTP-PStar)**3
  JAdd_Gradient(Ret_CTP) = 3.*ScaleFactorP * (CTP-PStar)**2
  IF (PRESENT(JAdd_Gradient2)) JAdd_Gradient2(Ret_CTP,Ret_CTP) = &
       6.*ScaleFactorP * (CTP-PStar)
  If (GeneralMode >= VerboseMode) Then
     Write(*,*) ' Cloud Top Pressure ',CTP,' too high: '
     Write(*,*) ' Additional Cost Function increased to ',JAdd
     Write(*,*) ' Gradient increased by ',JAdd_Gradient(Ret_CTP)
     Write(*,*)
  END IF
ELSE IF ( CTP < Cloud_Min_Pressure ) THEN
  JAdd = JAdd + ScaleFactorP * (Cloud_Min_Pressure - CTP)**3
  JAdd_Gradient(Ret_CTP) = -3.*ScaleFactorP * (Cloud_Min_Pressure - CTP)**2
  IF (PRESENT(JAdd_Gradient2)) JAdd_Gradient2(Ret_CTP,Ret_CTP) = &
       6.*ScaleFactorP * (Cloud_Min_Pressure - CTP)
  If (GeneralMode >= VerboseMode) Then
     Write(*,*) ' Cloud Top Pressure ',CTP,' too low: '
     Write(*,*) ' Additional Cost Function increased to ',JAdd
     Write(*,*) ' Gradient increased by ',JAdd_Gradient(Ret_CTP)
     Write(*,*)
  END IF
END IF

!---------------------------------------------------------------------
! 2. Cloud Cover Additional Cost Function
!---------------------------------------------------------------------

IF ( CloudCover > ThickCloud ) THEN
  write(*,*) JAdd, ScaleFactorF, CloudCover, ThickCloud
  JAdd = JAdd + ScaleFactorF * (CloudCover-ThickCloud)**3
  JAdd_Gradient(Ret_CloudCover) = 3.*ScaleFactorF * (CloudCover-ThickCloud)**2
  IF (PRESENT(JAdd_Gradient2)) &
       JAdd_Gradient2(Ret_CloudCover,Ret_CloudCover) = &
       6.*ScaleFactorF * (CloudCover-ThickCloud)
  If (GeneralMode >= VerboseMode) Then
     Write(*,*) ' Cloud Cover ',CloudCover,' too high: '
     Write(*,*) ' Additional Cost Function increased to ',JAdd
     Write(*,*) ' Gradient increased by ',JAdd_Gradient(Ret_CloudCover)
     Write(*,*)
  END IF
ELSE IF ( CloudCover < CloudFree ) THEN
  JAdd = JAdd + ScaleFactorF * (CloudFree - CloudCover)**3
  JAdd_Gradient(Ret_CloudCover) = &
       -3.*ScaleFactorF * (CloudFree - CloudCover)**2
  IF (PRESENT(JAdd_Gradient2)) &
       JAdd_Gradient2(Ret_CloudCover,Ret_CloudCover) = &
       6.*ScaleFactorF * (CloudFree - CloudCover)
  If (GeneralMode >= VerboseMode) Then
     Write(*,*) ' Cloud Cover ',CloudCover,' too low: '
     Write(*,*) ' Additional Cost Function increased to ',JAdd
     Write(*,*) ' Gradient increased by ',JAdd_Gradient(Ret_CloudCover)
     Write(*,*)
  END IF
END IF


End Subroutine NWPSAF_AdditionalCost_Cloud
