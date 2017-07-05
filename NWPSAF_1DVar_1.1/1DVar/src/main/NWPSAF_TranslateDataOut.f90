!+ Converts data back into a form suitable for output.

SUBROUTINE NWPSAF_TranslateDataOut ( &
     obnumber,   & ! in
     MeasurementObs, & ! in
     Background, & ! in
     RT_Params,  & ! inout
     UsedChans,  & ! in
     Obs)          ! inout

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
! Converts data back into a form suitable for output.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! NWPSAF_TranslateDataOut:
! 1.0     06/01/00 Orginal code based on ATOVS_TranslateDataOut.  
!                                                            A. Collard.
!                                                            Met. Office
! 1.1     01/03/01 Add ozone Roger Saunders
! 3.0.1   15/07/03 Output cloud retievals (M. Szyndel).      A. Collard.
! 3.0.4   05/03/04 Replace Plevels_RTmodel_Pa with Background%p.  A. Collard.
! 3.0.5   01/04/04 Tidied up output.                         A. Collard.
! 3.1.1   06/10/04 Obs% t2 % value set to T2 rather than TStar!  A. Collard. 
! 3.1.2   21/04/05 Corrected units for 1st guess CTP output.    E. Pavelin
!
!  Ticket Date     Comment
!--------------------------
!  28     22/02/12  Modified to output LWP if cloud liquid 
!                   retrieval is switched on TR Sreerekha
!  32     07/12/12  Added output for wind speed  A. Andersson (CM SAF)
!  32     18/11/13  Changes to output pressure units and formatting. P. Weston
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type, &
    ModelOb_Type

USE NWPSAFMod_Params, ONLY : &
    CloudyRetrieval,             &
    FileUnit_Retrieved_Profiles, &
    FileUnit_Retrieved_BTs,      &
    Humidity_Units,              & 
    Humidity_PPMV,               &   
    Humidity_MassMix,            & 
    Humidity_RH,                 &
    MwClwRetrieval,              &
    Retrieve_LWP,                &
    Ret_UWind,                   &
    Ret_VWind,                   &
    UsePCs,                      &
    CalcRadiance

USE NWPSAFMod_RTmodel, ONLY : &
    RTParams_Type,   &
    Num_RTlevels,    &      
    Prof_FirstQ,     &        
    Prof_FirstT,     &       
    Prof_FirstO3,    &       
    Prof_LastT,      &       
    Prof_LastO3,     &       
    Prof_T2,         &           
    Prof_Q2,         &           
    Prof_Tstar,      &        
    Prof_pstar,      &
    Prof_CTP,        &
    Prof_CloudCover, &
    Prof_LWP,        &
    Prof_Uwind,      &
    Prof_Vwind

USE NWPSAFMod_Constants, ONLY : &
    q_mixratio_to_ppmv

IMPLICIT NONE

INCLUDE 'NWPSAF_QSAT.interface'

! Subroutine arguments:
INTEGER, INTENT(IN)                     :: obnumber   ! Observation Number
REAL,    INTENT(IN)                     :: MeasurementObs(:) ! Observation 
TYPE(ModelOb_type), INTENT(IN)          :: Background ! Background data
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params  ! RT Model Input Data
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans  ! Chans to be used
TYPE(Ob_type), INTENT(INOUT)            :: Obs        ! Observed/Retrieval data 


! Local variables:
INTEGER :: level
INTEGER :: element
INTEGER :: IChan
REAL :: Pstar(1)
REAL :: IASIlevel_rh_out(Num_RTlevels)
REAL :: qsaturated(Num_RTlevels)
REAL :: SurfaceT(1)
REAL, POINTER :: MeasurementRet(:)

!-----------------------------------------------------------------------------

!----------------------------------------------
! 1. Convert temperature and humidity and ozone
!----------------------------------------------

SELECT CASE (Humidity_Units)
  CASE(Humidity_RH) 

     !Convert q to rh 
     !--------------------------------------------------------------------
     !Note: After conversion, it is necessary to constrain rh to < 1.0
     !      due to rounding errors.

     CALL NWPSAF_QSAT(                                         &
          qsaturated(1:Num_RTlevels),                   &
          RT_Params % RTguess(Prof_FirstT:Prof_LastT),  &
          Background % p(obnumber,1:Num_RTlevels),      & 
          Num_RTlevels ) 
     DO level = 1, Num_RTlevels
        IASIlevel_rh_out(level) = &
             EXP(RT_Params % RTguess(Prof_FirstQ+level-1)) / &
             qsaturated(level)
        IASIlevel_rh_out(level) = MIN( 1.0, IASIlevel_rh_out(level) )
     END DO
     Pstar(1) = RT_Params % RTguess(Prof_Pstar) * 100.
     SurfaceT(1)=RT_Params % RTguess(Prof_T2)
     CALL NWPSAF_QSAT( qsaturated(1),                        &
          SurfaceT,                                   &
          Pstar(1),                                   &
          1 )
     RT_Params % RTguess(Prof_q2) = &
          1.0 * EXP(RT_Params % RTguess(Prof_q2)) / qsaturated(1)
     RT_Params % RTguess(Prof_q2) = &
          MIN( 1.0, RT_Params % RTguess(Prof_q2) )

  CASE (Humidity_PPMV)

     !Convert q to ppmv 
     DO level = 1, Num_RTlevels
        IASIlevel_rh_out(level) = &
             EXP(RT_Params % RTguess(Prof_FirstQ+level-1)) * q_mixratio_to_ppmv
     END DO
     RT_Params % RTguess(Prof_q2) = &
          EXP(RT_Params % RTguess(Prof_q2)) * q_mixratio_to_ppmv
     Pstar(1) = RT_Params % RTguess(Prof_Pstar) * 100.0
  CASE (Humidity_MassMix)
     
     !Convert q to ppmv 
     DO level = 1, Num_RTlevels
        IASIlevel_rh_out(level) = &
             EXP(RT_Params % RTguess(Prof_FirstQ+level-1))
     END DO
     RT_Params % RTguess(Prof_q2) = EXP(RT_Params % RTguess(Prof_q2)) 
     Pstar(1) = RT_Params % RTguess(Prof_Pstar) * 100.0
END SELECT

!------------------------
! 2. Assign Obs variables
!------------------------
!2.1) Surface data
!----
Obs% pstar % value = RT_Params % RTguess(Prof_Pstar)
Obs% t2 % value    = RT_Params % RTguess(Prof_T2)
Obs% rh2 % value   = RT_Params % RTguess(Prof_q2)
Obs% Tskin % value = RT_Params % RTguess(Prof_Tstar)
Obs% CTP% value = RT_Params % RTguess(Prof_CTP)
Obs% CldFrac% value = RT_Params % RTguess(Prof_CloudCover)
Obs% LWP% value = RT_Params % RTguess(Prof_LWP)

!2.2) Multi-level data
!----

Obs% t(:)% value = &
     RT_Params % RTguess(Prof_FirstT:Prof_LastT)
Obs% rh(:)% value = IASIlevel_rh_out(1:Num_RTlevels)
Obs% ozone(:)% value =&
     RT_Params % RTguess(Prof_FirstO3:Prof_LastO3)

!-----------------
! 4. Output to files
! ------------------

! ---------------------------
! 4.1. Output Profiles First
! ---------------------------

WRITE(FileUnit_Retrieved_Profiles, *) 'Observation = ',obnumber
WRITE(FileUnit_Retrieved_Profiles, *) &
'                         Retrieval                        Background' 
SELECT CASE (Humidity_Units)
  CASE(Humidity_RH)
     WRITE(FileUnit_Retrieved_Profiles, *) &
          ' Pressure (hPa) T (K)    Rel.Hum.     Ozone   '// &
          '     T (K)   Rel.Hum.      Ozone'
     DO Element=1,Num_RTLevels
        WRITE(FileUnit_Retrieved_Profiles,FMT='(F12.3,2(F10.3,E12.4,E12.4))') &
             Background % p(Obnumber,Element)/100.0, &
             Obs % t(Element) % Value,               & 
             Obs % rh(Element) % Value,              &
             Obs % ozone(Element) % Value,           &
             Background % t(obnumber, Element ),     & 
             Background % rh(obnumber, Element),     & 
             Background % ozone(obnumber, Element) 
     END DO
  CASE(Humidity_PPMV)
     WRITE(FileUnit_Retrieved_Profiles, *) &
          ' Pressure (hPa) T (K)    q (ppmv)     Ozone   '// &
          '     T (K)   q (ppmv)      Ozone'
     !xxx     DO Element=1,Num_RTLevels
     DO Element=Num_RTLevels,1,-1
        WRITE(FileUnit_Retrieved_Profiles,FMT='(F12.3,2(F10.3,E12.4,E12.4))') &
             Background % p(Obnumber,Element)/100.0,       &
             Obs % t(Element) % Value,                     & 
             Obs % rh(Element) % Value,                    &
             Obs % ozone(Element) % Value,                 &
             RT_Params % RTBack(Prof_FirstT+Element-1),    & 
             EXP(RT_Params % RTBack(Prof_FirstQ+Element-1))*q_mixratio_to_ppmv,&
             RT_Params % RTBack(Prof_FirstO3+Element-1)
     END DO
  CASE(Humidity_MassMix)
     WRITE(FileUnit_Retrieved_Profiles, *) &
          ' Pressure (hPa) T (K)    q (kg/kg)    Ozone   '// &
          '     T (K)   q (kg/kg)     Ozone'
     DO Element=1,Num_RTLevels
        WRITE(FileUnit_Retrieved_Profiles,FMT='(F12.3,2(F10.3,E12.4,E12.4))') &
             Background % p(Obnumber,Element)/100.0,                    &
             Obs % t(Element) % Value,                                  & 
             Obs % rh(Element) % Value,                                 &
             Obs % ozone(Element) % Value,                              &
             RT_Params % RTBack(Prof_FirstT+Element-1),                 & 
             EXP(RT_Params % RTBack(Prof_FirstQ+Element-1)),            &
             RT_Params % RTBack(Prof_FirstO3+Element-1)
     END DO
END SELECT
  
WRITE(FileUnit_Retrieved_Profiles, &
     FMT='(''Surface Temperature (K):      '', 2F10.3)') &
     Obs% t2% Value, &
     Background% t2(obnumber)
SELECT CASE (Humidity_Units)
CASE(Humidity_RH)
   WRITE(FileUnit_Retrieved_Profiles, &
        FMT='(''Surface Relative Humidity:    '', 2F10.3)') &
        Obs% rh2% Value, &
        Background% rh2(obnumber)
CASE(Humidity_PPMV)
   WRITE(FileUnit_Retrieved_Profiles, &
        FMT='(''Surface Humidity (ppmv):      '', 2F10.3)') &
        Obs% rh2% Value,                                    &
        EXP(RT_Params % RTBack(Prof_q2))* q_mixratio_to_ppmv
CASE(Humidity_MassMix)
   WRITE(FileUnit_Retrieved_Profiles, &
        FMT='(''Surface Humidity (kg/kg):     '', 2F10.3)') &
        Obs% rh2% Value,                          &
        EXP(RT_Params % RTBack(Prof_q2))
END SELECT

IF (MwClwRetrieval .AND. Retrieve_LWP) THEN
   ! N.B. 1st guess CTP and Cloud fraction reported, not background.
   WRITE(FileUnit_Retrieved_Profiles, &
        FMT='(''LWP    (kg/m2):              '', 2F10.3)') &
        Obs% LWP % Value, &
        RT_Params % RT1stguess(Prof_LWP)
   !WRITE(FileUnit_Retrieved_Profiles,FMT=*) 
END IF

WRITE(FileUnit_Retrieved_Profiles, &
     FMT='(''Skin Temperature (K):         '', 2F10.3)') &
     Obs% Tskin% Value, &
     Background% Tskin(obnumber)
  
WRITE(FileUnit_Retrieved_Profiles, &
     FMT='(''Surface Pressure (hPa):        '', 2F9.3)') &
     Pstar(1)/100.0, Background% Pstar(obnumber)/100.0

!wind speed
IF (Ret_UWind /= 0 .AND. Ret_VWind /=0 ) THEN
   WRITE(FileUnit_Retrieved_Profiles, &
   FMT='(''Wind Speed: (u,v,U10,U10Back): '', 4F8.3)') &
        RT_Params % RTguess(Prof_Uwind),RT_Params % RTguess(Prof_Vwind), &
        SQRT(RT_Params % RTguess(Prof_Uwind)**2+RT_Params % RTguess(Prof_Vwind)**2), &
        SQRT(Background % u10(obnumber)**2+Background % v10(obnumber)**2)
END IF

WRITE(FileUnit_Retrieved_Profiles,FMT=*) 

IF (CloudyRetrieval) THEN
   ! N.B. 1st guess CTP and Cloud fraction reported, not background.

  WRITE(FileUnit_Retrieved_Profiles, &
       FMT='(''Cloud Top Pressure (hPa):      '', 2F9.0)') &
       Obs% CTP% Value, &
       RT_Params % RT1stguess(Prof_CTP)
  WRITE(FileUnit_Retrieved_Profiles,FMT=*) 

  WRITE(FileUnit_Retrieved_Profiles, &
       FMT='(''Cloud Fraction :               '', 2F9.5)') &
       Obs% CldFrac% Value, &
       RT_Params % RT1stguess(Prof_CloudCover)
  WRITE(FileUnit_Retrieved_Profiles,FMT=*) 
END IF



WRITE(FileUnit_Retrieved_Profiles, &
     FMT='(''No. of Iterations:             '',I3)') &
     Obs% NIter
IF (ABS(Obs % Jcost) >= 1.e4) Obs % Jcost = -999.999
IF (Obs % Jcost_Gradient >= 1.e4) Obs % Jcost_Gradient = -999.999
WRITE(FileUnit_Retrieved_Profiles, &
     FMT='(''Normalised Cost Function:    '','// &
     'F8.3,'' Normalised Gradient:    '', F8.3)') &
     Obs% Jcost, Obs% Jcost_Gradient
WRITE(FileUnit_Retrieved_Profiles,FMT=*) &
     '--------------------------------------------------------------'

! ----------------------------------------
! 5.2. Now Output Brightness Temperatures
! ----------------------------------------

IF ( UsePCs ) THEN
  MeasurementRet=> RT_Params % PCscores
ELSE IF ( CalcRadiance ) THEN
  MeasurementRet=> RT_Params % TotalRadiances
ELSE
  MeasurementRet=> RT_Params % TotalBTs
ENDIF

WRITE(FileUnit_Retrieved_BTs, *) 'Observation = ',obnumber
WRITE(FileUnit_Retrieved_BTs, *) 'Number of Channels Used = ',&
     UsedChans % NumChans
WRITE(FileUnit_Retrieved_BTs, *) &
'Channel  Background  Observed   Retrieved          ' 
DO IChan=1,UsedChans % NumChans
   WRITE(FileUnit_Retrieved_BTs,FMT='(I7,3F11.3)') &
        UsedChans % Channels(IChan), &
        BackGround % BriTemp(UsedChans % Channels(IChan)), &
        MeasurementObs(UsedChans % Channels(IChan)), &
        MeasurementRet(IChan)
END DO

END SUBROUTINE NWPSAF_TranslateDataOut
