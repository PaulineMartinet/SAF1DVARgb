!+  Make a cloud structure
Subroutine NWPSAF_CloudStructure(&
     RT_Params,     & !inout
     cloud_structure) !out
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
!      Make a cloud structure                                           
!      Create Normalized background cloud liquid water profile 
!      The cloud structure MUST NOT CHANGE during the minimization process
!      LWP_back MUST be computed before calling SSMI_CloudStructure
!
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
! 28     21/02/12 Original code based on SSMI_Cloud_Structure.  TR Sreerekha.
! 31     24/07/13 Bug fix: Test that background not guess CLW 
!                 is greater than 0.                            P Weston
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------
!

! Modules used:

USE NWPSAFMod_Constants, ONLY :  & 
    epsilon

USE NWPSAFMod_RTmodel, ONLY : &
    Num_RTlevels, &
    RTParams_Type, &
    Prof_FirstT, &
    Prof_LastT, &
    Prof_FirstCLW, &
    Prof_LastCLW, &
    Prof_FirstQ, &
    Prof_LastQ, &
    Prof_LWP

USE NWPSAFMod_LiquidWater, ONLY : &
    SVP

USE rttov_const, ONLY : &
    gravity

IMPLICIT NONE

TYPE(RTParams_Type)     , INTENT(IN)  :: RT_Params ! RT Model Data
REAL,                     INTENT(OUT) :: cloud_structure(:)

!local
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_CloudStructure"
REAL , PARAMETER             :: ZTOPCL = 253  
REAL , PARAMETER             :: ZRHCL = 0.8

REAL                         :: Zv(Num_RTlevels)                                
REAL                         :: ZSUMP
REAL                         :: ESAT,QSAT
REAL                         :: RH
REAL                         :: deltaP
INTEGER                      :: i
INTEGER                      :: KCT ! cloud top index                  
INTEGER                      :: KCB ! cloud bottom index                      
LOGICAL                      :: LDCL !
LOGICAL                      :: Lrh !Useful for debugging
REAL                         :: Temperature_Profile(Num_RTlevels)
REAL                         :: Specific_Q_Profile(Num_RTlevels)


!----------------End of header -----------------------------------------

cloud_structure(:) = 0.0
Lrh=.FALSE.
KCB=-99
KCT=-99

LDCL = .FALSE.
IF(RT_Params % RTBack(Prof_LWP) >= 1e-20) THEN 
  cloud_structure(:) = &
        RT_Params % RTBack(Prof_FirstCLW:Prof_LastCLW)/RT_Params % RTBack(Prof_LWP)
  LDCL = .TRUE.
ENDIF

IF (.NOT.LDCL) THEN
!   Based on L. Phalippou approach in SSMI 1DVAR --1996 code
!      developed at ECMWF

!   Make a cloud structure based on relative humidity
!   and temperature

    Lrh=.FALSE.
    Zv(:)=0
    Temperature_Profile(:)=RT_Params % RTguess(Prof_FirstT:Prof_LastT)
    Specific_Q_Profile(:) = RT_Params % RTguess(Prof_FirstQ:Prof_LastQ)
    DO i=1,Num_RTlevels
      ESAT=SVP(Temperature_Profile(i))
      IF ( ESAT  > RT_Params % Pressure_Pa(i)) THEN
        QSAT    = 1.0
      ELSE
        QSAT  = epsilon / ( RT_Params % Pressure_Pa(i)/ESAT - (1.-epsilon)) ! spec.hum (kg/kg)
      ENDIF  
      RH = Specific_Q_Profile(i)  /(QSAT*1000.0)
      IF ( Temperature_Profile(i) > ZTOPCL .AND. RH >= ZRHCL )  THEN
        Zv(I)=1.
        Lrh=.TRUE.
      ENDIF
    ENDDO

  !Cloud top and bottom

  KCT=0
  KCB=0
  DO i=1,Num_RTlevels
    IF (Zv(i) > 0.0) THEN
      KCT=I
      EXIT       
    ENDIF
  ENDDO
  DO i=Num_RTlevels,1,-1
    IF (Zv(i) > 0.0) THEN
      KCB=I
      EXIT            
    ENDIF
  ENDDO
  !If only 1 cloudy level forces a second cloudy levels
  IF ( KCT == KCB .AND. KCT >  2 ) THEN
    KCT=KCT-1
    Zv(KCT)=1.
  ENDIF

  !If still no cloud, force a structure on the last levels
  IF (KCT == 0) THEN
    KCT=Num_RTlevels-5
    KCB=Num_RTlevels-3
    Zv(KCT:KCB) = 1.0 
  ENDIF

  !Normalize the structure
  ZSUMP=0.
  DO i=1,Num_RTlevels-1
    deltaP=RT_Params % Pressure_Pa(i+1)-RT_Params % Pressure_Pa(i)
    ZSUMP=ZSUMP+deltaP*(Zv(i)+Zv(i+1))
  ENDDO
  ZSUMP=(0.5*ZSUMP)*100./gravity
  Zv(:) = Zv(:)/ZSUMP
  cloud_structure(:)= Zv(:)

ENDIF


End Subroutine NWPSAF_CloudStructure  
