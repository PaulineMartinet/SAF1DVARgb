FUNCTION NWPSAF_check_temperatures( &
     Level_temperatures,        &  ! in
     Number_of_levels,          &  ! in
     Surface_air_temperature,   &  ! in
     Surface_skin_temperature ) &  ! in
     RESULT (Valid_data)
  
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
! Description: Checks whether input temperatures are inside gross limits
  
! Modules used:

  USE NWPSAFMod_Params, ONLY : &
       GeneralMode, &
       VerboseMode

  USE NWPSAFMod_Constants, ONLY : &
       MaxTemperature, &
       MinTemperature

  IMPLICIT NONE

  ! Subroutine arguments:
  INTEGER, INTENT(IN) :: Number_of_levels
  REAL,    INTENT(IN) :: Surface_air_temperature
  REAL,    INTENT(IN) :: Surface_skin_temperature
  REAL,    INTENT(IN) :: Level_temperatures(:)

  ! Local variables:
  INTEGER :: level        ! Loop variable
  REAL    :: Tdata        ! Temporary storage for temperature data
  LOGICAL :: Valid_data   ! Result of function = .FALSE. for bad data

!----------------------------------------------------------------------------

  !---------------------------
  !1. Level temperatures check
  !---------------------------

  Valid_Data = .TRUE.
  
  DO level = 1,Number_of_levels
    Tdata = Level_temperatures(level)
    IF ( (Tdata > MaxTemperature) .OR. (Tdata < MinTemperature) ) THEN
       Valid_Data = .FALSE.
       IF ( GeneralMode >= VerboseMode ) &
            WRITE(*,*) 'INVALID DATA: Bad temperature = ',Tdata, &
            ' level ',level
    END IF
 END DO


 !--------------------------------
 !2. Surface air temperature check
 !--------------------------------

 Tdata = Surface_air_temperature
 IF ( (Tdata > MaxTemperature) .OR. (Tdata < MinTemperature) ) THEN
    Valid_Data = .FALSE.
    IF ( GeneralMode >= VerboseMode ) &
         WRITE(*,*) 'INVALID DATA: Bad surface air temperature = ',Tdata
 END IF
 
 !---------------------------------
 !3. Surface skin temperature check
 !---------------------------------
 
 Tdata = Surface_skin_temperature
 IF ( (Tdata > MaxTemperature) .OR. (Tdata < MinTemperature) ) THEN
    Valid_Data = .FALSE.
    IF ( GeneralMode >= VerboseMode ) &
         WRITE(*,*) 'INVALID DATA: Bad surface skin temperature = ',Tdata
 END IF
 
END FUNCTION NWPSAF_check_temperatures
