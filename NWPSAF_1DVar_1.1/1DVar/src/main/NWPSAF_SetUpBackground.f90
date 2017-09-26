!+ Set up the background profile vectors

SUBROUTINE NWPSAF_SetUpBackground ( &
     BackGround,      & ! in
     obnumber,        & ! in
     RT_Params,       & ! inout
     Valid_data )       ! out

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
! Description: Insert background data into the 1DVar input vector, performing
!              any conversions and extra calculations, plus quality control.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     22/02/01 Original version based on ATOVS_SetUpBackground.
!                                                     A.D. Collard Met Office
! 1.1     02/03/01 Rearranged to include ozone R. Saunders
! 3.0.1   15/07/03 Include cloud data (M.D.E.Szyndel).  A.D. Collard.
! 3.0.4   05/03/04 Replace PLevels_RTModel_Pa with Background%P. 
!                  Remove QC on background temperature profile. A. Collard.
! 3.0.5   01/08/07 Fix problem with newer NagWARE f95 compilers. E. Pavelin.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
!  28       22/02/12 Added cloud liquid water to teh background TR Sreerekha
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------
! Modules used:

USE NWPSAFMod_ObsInfo, ONLY : &
    ModelOb_Type

USE NWPSAFMod_Constants, ONLY : &
    MinSurfaceP, &
    MaxSurfaceP, &
    Tolerance, &
    MissData_R, &
    Min_q,        &
    q_mixratio_to_ppmv 

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,        &
    VerboseMode,        &
    DoTextrapolation,   &
    RTsea,              &
    Humidity_Units,     & 
    Humidity_PPMV,      &   
    Humidity_MassMix,   & 
    Humidity_RH  , &
    !PM
    retrieval_in_log
    !PM

USE NWPSAFMod_RTmodel, ONLY : &
    RTParams_Type, &
    Num_RTlevels, &      
    Prof_FirstT, &        
    Prof_Firstq, &        
    Prof_FirstO3, & 
    Prof_FirstCLW, & 
    Prof_LastT, &        
    Prof_Lastq, &        
    Prof_LastO3, &
    Prof_LastCLW, & 
    Prof_LWP, &
    Prof_CTP, &           
    Prof_CloudCover, &    
    Prof_T2, &            
    Prof_Tstar, &         
    Prof_q2, &            
    Prof_Pstar, &         
    Prof_Uwind, &         
    Prof_Vwind

USE NWPSAFMod_LiquidWater, ONLY : &
    Layers_to_LWP

IMPLICIT NONE

INCLUDE 'NWPSAF_Check_Temperatures.interface'
INCLUDE 'NWPSAF_StratosExtrap.interface'
INCLUDE 'NWPSAF_QSAT.interface'

! Subroutine arguments:
TYPE(ModelOb_type),INTENT(IN) :: Background ! Background data
INTEGER, INTENT(IN) :: obnumber        ! Profile number
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params ! Data for RT Model
LOGICAL, INTENT(OUT) :: Valid_Data      ! Error flag


! Local variables:
INTEGER :: level           ! loop variable
INTEGER :: Prof_T950hPa    ! Position of 950hPa temperature in RT profile
INTEGER :: Prof_T1000hPa   ! Position of 1000hPa temperature in RT profile
INTEGER :: Tmp(1)          ! Temporary array used with MINLOC intrinsic.
REAL    :: rhumidity       ! temporary storage for relative humidity
REAL    :: Tstar           ! surface temperature
REAL    :: NewTemp         ! temperature storage
REAL    :: Pstar           ! surface pressure
REAL    :: SurfaceP(1)     ! surface pressure
REAL    :: SurfaceT(1)     ! surface temperature
REAL, ALLOCATABLE :: qsaturated(:) ! specific humidity at saturation

!-----------------------------------------------------------------------------

!-------------
!0. Initialize
!-------------

Valid_data = .TRUE.
RT_Params % RTBack(:)    = MissData_R


ALLOCATE(qsaturated(Num_RTLevels+1))

!-----------
!1. Pressure
!-----------
!Note that I/P pressures are in Pa, while RT model pressures are in hPa

Pstar = Background % pstar(obnumber)*0.01
RT_Params % RTBack(Prof_Pstar) = Pstar
IF ( (Pstar > MaxSurfaceP) .OR. (Pstar < MinSurfaceP) ) THEN
   Valid_Data = .FALSE.
   IF ( GeneralMode >= VerboseMode ) THEN
      WRITE(*,*) 'Bad surface pressure: ',Pstar
   END IF
END IF

RT_Params % Pressure_Pa => Background % p(obnumber,:)

!---------------
!2. Temperatures
!---------------

Tmp = MINLOC(ABS(RT_Params % Pressure_Pa -  95000.))
Prof_T950hPa = Tmp(1) + Prof_FirstT - 1

Tmp = MINLOC(ABS(RT_Params % Pressure_Pa - 100000.))
Prof_T1000hPa = Tmp(1) + Prof_FirstT - 1

!2.1) Assign model values
!----

RT_Params % RTBack(Prof_FirstT:Prof_LastT) = &
     Background % t(obnumber,:)

!surface air temperature
RT_Params % RTBack(Prof_T2) = Background % t2(obnumber)

!surface skin temperature
RT_Params % RTBack(Prof_Tstar) = Background % Tskin(obnumber)

!2.2) Check model values
!----
!Replace any low levels that have missing data with 2m value. Also,
!determine the top extent of the interpolated levels.
!Note: any levels not covered in the OPS background interpolation routine
!     (i.e. out of range for interpolation) will be missing data here

DO level = Prof_LastT, Prof_FirstT, -1
   IF ( ABS(RT_Params % RTBack(level) - MissData_R) <= Tolerance ) THEN
      RT_Params % RTBack(level) = RT_Params % RTBack(Prof_T2)
   ELSE !finish when all missing levels above the surface have been covered
      EXIT
   END IF
END DO

Valid_data = Valid_Data .AND. NWPSAF_Check_Temperatures(  &
     RT_Params % RTBack(Prof_FirstT:Prof_LastT),        &
     Num_RTlevels,                                      &
     RT_Params % RTBack(Prof_T2),                       &
     RT_Params % RTBack(Prof_Tstar) )

!2.3) Extrapolate temperature profile to 0.01mb.
!----

IF ( Valid_Data .AND. DoTextrapolation ) THEN
   CALL NWPSAF_StratosExtrap(         &
        RT_Params % Pressure_Pa(:), &                 ! in
        RT_Params % RTBack(Prof_FirstT:Prof_LastT)  ) ! inout
END IF

!2.4) Reset low level temperatures over seaice and cold, low land
!----
Tstar = RT_Params % RTBack(Prof_Tstar)
IF ( RT_Params % RTSurfaceType /= RTSea .AND. &
     Tstar < 271.4 .AND. Pstar > 950.0 ) THEN
   NewTemp = RT_Params % RTBack(Prof_T950hPa)
   IF ( Pstar > 1000.0 ) THEN
      NewTemp = MAX(NewTemp,RT_Params % RTBack(Prof_T1000hPa))
   END IF
   NewTemp = MIN(NewTemp,271.4)
   IF ( RT_Params % RTBack(Prof_T1000hPa) < NewTemp ) &
        RT_Params % RTBack(Prof_T1000hPa) = NewTemp
   IF ( RT_Params % RTBack(Prof_T2) < NewTemp ) &
        RT_Params % RTBack(Prof_T2) = NewTemp
   IF ( RT_Params % RTBack(Prof_Tstar) < NewTemp ) &
        RT_Params % RTBack(Prof_Tstar) = NewTemp
END IF

!----------------------
!3. Specific humidities
!----------------------
!Note: The minimisation is done for Ln(mass mixing ratio), so 
!  RT_Params % RTBack elements for humidity and ozone are changed accordingly.
!
!Initialize
!PM
IF (retrieval_in_log) THEN
    RT_Params % RTBack(Prof_FirstQ:Prof_LastQ) = LOG(min_q)
ELSE
RT_Params % RTBack(Prof_FirstQ:Prof_LastQ) = min_q
ENDIF
!PM
Humidities: IF ( Valid_Data ) THEN

   !3.1) Convert surface humidity to Ln(q)
   !----
   SurfaceT(1) = RT_Params % RTBack(Prof_T2)
   SurfaceP(1) = Background % pstar(obnumber)
   CALL NWPSAF_QSAT( qsaturated(Num_RTLevels+1), &
        SurfaceT,                         &
        SurfaceP,                         &
        1 )
   rhumidity = Background % rh2(obnumber)
   SELECT CASE (Humidity_Units)
     CASE (Humidity_PPMV)
     !PM
     IF (retrieval_in_log) THEN
     !PM
        RT_Params % RTBack(Prof_q2) = &
             LOG(rhumidity / q_mixratio_to_ppmv)
     ELSE
        RT_Params % RTBack(Prof_q2) = &
             rhumidity / q_mixratio_to_ppmv
    ENDIF
     !PM
     CASE(Humidity_MassMix)
     IF (retrieval_in_log) THEN
        RT_Params % RTBack(Prof_q2) = LOG(rhumidity)
     ELSE
        RT_Params % RTBack(Prof_q2) = rhumidity
     ENDIF
     CASE(Humidity_RH)
     IF (retrieval_in_log) THEN
        RT_Params % RTBack(Prof_q2) = &
             LOG(qsaturated(Num_RTLevels+1) * rhumidity)
    ELSE
        RT_Params % RTBack(Prof_q2) = &
             qsaturated(Num_RTLevels+1) * rhumidity
    ENDIF
    !PM
   END SELECT
          
   !3.2) Convert humidity from given units to Ln(q) 
   !----
   
   CALL NWPSAF_QSAT( &
        qsaturated(1:Num_RTLevels),                 &
        RT_Params % RTBack(Prof_FirstT:Prof_LastT), &
        Background % P (Obnumber,1:Num_RTlevels),   &
        Num_RTLevels)        

   ! Convert QSAT (kg/kg) to input units
   IF (Humidity_Units == Humidity_PPMV) &
        qsaturated(1:Num_RTLevels) = &
        qsaturated(1:Num_RTLevels) * q_mixratio_to_ppmv
   
   DO level = 1, Num_RTLevels
      rhumidity = Background % &
           rh(obnumber,level)
      SELECT CASE (Humidity_Units)
        CASE (Humidity_PPMV)
        !PM
        IF (retrieval_in_log) THEN
           RT_Params % RTBack(Prof_FirstQ+level-1) = &
                LOG(rhumidity / q_mixratio_to_ppmv)
        ELSE
          RT_Params % RTBack(Prof_FirstQ+level-1) = &
                rhumidity / q_mixratio_to_ppmv
        ENDIF
        CASE(Humidity_MassMix)
       IF (retrieval_in_log) THEN        
           RT_Params % RTBack(Prof_FirstQ+level-1) = LOG(rhumidity)
       ELSE
           RT_Params % RTBack(Prof_FirstQ+level-1) = rhumidity
       ENDIF
        CASE(Humidity_RH)
        IF (retrieval_in_log) THEN
           RT_Params % RTBack(Prof_FirstQ+level-1) = &
                LOG(qsaturated(level) * rhumidity)
        ELSE
            RT_Params % RTBack(Prof_FirstQ+level-1) = &
                qsaturated(level) * rhumidity
        ENDIF
        !PM
      END SELECT
   END DO

END IF Humidities

!--------
!4. Ozone
!--------
!4.1) Assign model values (change units here if necessary)
!----
RT_Params % RTBack(Prof_FirstO3:Prof_LastO3) = &
     Background % ozone(obnumber,:)

!---------------------
!5. Cloud Liquid Water
!---------------------
RT_Params % RTBack(Prof_FirstCLW:Prof_LastCLw) = Background % clw(obnumber,:)
! Calculate background LWP from the background cloud liquid water
RT_Params % RTBack(Prof_LWP) = &
  Layers_to_LWP(RT_Params % RTBack(Prof_FirstCLW:Prof_LastCLW), &
                RT_Params % Pressure_Pa                         )

!------------------
!6. QC profile data
!------------------

!6.1) Check temperatures are within limits
!----
Valid_data = Valid_Data .AND. NWPSAF_Check_Temperatures(       &
     RT_Params % RTBack(Prof_FirstT:Prof_LastT),             &  
     Num_RTlevels,                                           &
     RT_Params % RTBack(Prof_T2),                            &
     RT_Params % RTBack(Prof_Tstar) )


!6.2) Check all specific humidities for super-saturation
!----

Check_humidities: IF ( Valid_Data ) THEN
   
   !Levels
   !------
   DO level = 1, Num_RTlevels
!PM
      !Constrain to less than qsaturated
     IF (retrieval_in_log) THEN 
      RT_Params % RTBack(level+Prof_FirstQ-1) = &
           MIN(RT_Params % RTBack(level+Prof_FirstQ-1), &
           LOG(qsaturated(level)))
      ELSE
        RT_Params % RTBack(level+Prof_FirstQ-1) = &
           MIN(RT_Params % RTBack(level+Prof_FirstQ-1), &
           qsaturated(level))
     ENDIF
      !Constrain to a minimum value
      IF (retrieval_in_log) THEN
      RT_Params % RTBack(level+Prof_FirstQ-1) = &
           MAX( RT_Params % RTBack(level+Prof_FirstQ-1),LOG(min_q))
      ELSE
        RT_Params % RTBack(level+Prof_FirstQ-1) = &
           MAX( RT_Params % RTBack(level+Prof_FirstQ-1),min_q)   
      ENDIF
!PM
      
   END DO
   
   !Surface
   !-------
   !Constrain to a maximum of 110% of qsaturated
 !PM  
   IF (retrieval_in_log) THEN
   RT_Params % RTBack(Prof_q2) = MIN(RT_Params % RTBack(Prof_q2), &
        LOG(1.1*qsaturated(Num_RTLevels+1)))
   ELSE
   RT_Params % RTBack(Prof_q2) = MIN(RT_Params % RTBack(Prof_q2), &
        1.1*qsaturated(Num_RTLevels+1))   
   ENDIF
   !PM
   
   !Constrain to a minimum value
!PM
   IF (retrieval_in_log) THEN
   RT_Params % RTBack(Prof_q2) = MAX( RT_Params % RTBack(Prof_q2), LOG(Min_q) )
   ELSE
   RT_Params % RTBack(Prof_q2) = MAX(RT_Params % RTBack(Prof_q2), Min_q )
   ENDIF
!PM
   
   
END IF Check_humidities

DEALLOCATE(qsaturated)


IF (.NOT. Valid_Data) THEN
   IF ( GeneralMode >= VerboseMode ) THEN
      WRITE(*,*) 'INVALID DATA: Bad extrapolated profile'
      write(*,*) 'Obnumber=',Obnumber
      WRITE(*,*) 't: ',(RT_Params % RTBack(level+Prof_FirstT-1), &
           level = 1,Num_RTlevels)
      WRITE(*,*) 'q: ',(RT_Params % RTBack(level+Prof_FirstQ-1), &
           level = 1,Num_RTlevels)
   END IF
END IF

!-------
!7. Wind
!-------

RT_Params % RTBack(Prof_uwind) = Background % u10(obnumber)
RT_Params % RTBack(Prof_vwind) = Background % v10(obnumber)


!--------
!8. Cloud
!--------

RT_Params % RTBack(Prof_cloudcover) = Background % CldFrac(obnumber)
RT_Params % RTBack(Prof_CTP)        = Background % CTP(obnumber)*0.01

END SUBROUTINE NWPSAF_SetUpBackground
