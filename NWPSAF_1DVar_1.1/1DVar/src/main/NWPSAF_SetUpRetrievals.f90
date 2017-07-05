    SUBROUTINE NWPSAF_SetUpRetrievals()

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
! Description: Reads in retrieval options from input file and sets up 
!        the relevant arrays accordingly.
!
! As, in general, the profile vector to be retrieved will not 
! contain all the levels given in the B-matrix which in turn will
! not contain all the levels used in the radiative transfer code,
! here we set up two vectors (Retrieved_Elements and B_ElementsUsed)
! of indices used to indicate which levels are being used.  These 
! vectors are both Num_ProfElementsUsed long.   
!
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.0      08/11/01  Original Version.  A.D. Collard.  Met Office.
! 3.0.1    17/07/03  Add flag for cloudy retrievals.   A. Collard.
! 3.0.2    04/02/04  Add Ret_CTP and Ret_CloudCover.   A. Collard.
! 3.0.3    02/03/04  Add check that RTModel is suitable for cloudy
!                    retrievals.                       A. Collard. 
! 3.0.4    08/03/04  Rewrite to use a (simpler) namelist.  A. Collard.
! 3.0.5    14/04/04  Remove automatic use of cloud additional cost function
!                    in doing CloudyRetrieval.         A. Collard.
! 3.0.6    05/04/05  Added support for RTTOV8.         E. Pavelin.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! Ticket  Date      Comment
! ------- --------  -------
! 28      22/02/12  Added options for cloud liquid water retrieval. TR Sreerekha
! 32      07/12/12  Added options for wind speed retrieval with internally
!                   set B-matrix entry                  A. Andersson (CM SAF)
! 31      13/06/13  Added support for RTTOV11.          P. Weston
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY:      &
     StatusFatal,              &
     DetectCloud,              &
     Retrieved_Elements,       &
     Ret_FirstQ,               &
     Ret_LastQ,                &
     Ret_Q2,                   &      
     Ret_CloudCover,           &      
     Ret_CTP,                  & 
     Ret_CLW,                  &
     Ret_UWind,                &
     Ret_VWind,                &
     B_ElementsUsed,           & 
     CloudyRetrieval,          &
     MwClwRetrieval,           &
     Retrieve_LWP,             &
     CalcRadiance

USE NWPSAFMod_RTModel, ONLY: &
     ProfSize,             &
     Num_ProfElementsUsed, &
     Prof_FirstT, &
     Prof_FirstQ, &
     Prof_FirstO3, &
     Prof_T2, &
     Prof_q2, &
     Prof_Tstar, &
     Prof_pstar, &
     Prof_uwind, &
     Prof_vwind, &
     Prof_LWP,   &
     Prof_CTP,   &
     Prof_CloudCover, &
     Num_WetLevels

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_OpenFile.interface'

! Subroutine arguments
!---- None ----

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "NWPSAF_SetUpRetrievals"
CHARACTER(LEN=*), PARAMETER :: Filename = 'Retrieval.NL'
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER (LEN=5) :: Status = "OLD"

! Local variables:
CHARACTER(LEN=80)  :: ErrorMessage(3)=' '  ! Message for NWPSAF_Report

INTEGER :: I 
INTEGER :: fileunit             ! I/O unit number
INTEGER :: ReadStatus

INTEGER :: B_ElementsUsed_Tmp(1000)     = -999
INTEGER :: Retrieved_Elements_Tmp(1000) = -999

! Namelist
!-----------

! For profile variables, the array elements here are
! First (lowest pressure) level, number of levels and 
! the B-matrix index corresponding to the first level
! for this variable.
INTEGER :: Temperature(3) = (/ 0, 0, 0 /)        
INTEGER :: Humidity(3)    = (/ 0, 0, 0 /)   
INTEGER :: Ozone(3)       = (/ 0, 0, 0 /)
! For single level variables, the array elements here are
! Do we retrieve it (0=No,1=Yes) and the B-matrix index 
! corresponding to this variable.
INTEGER :: Surface_Temperature(2) = (/ 0, 0 /)
INTEGER :: Surface_Humidity(2)    = (/ 0, 0 /)
INTEGER :: Surface_Pressure(2)    = (/ 0, 0 /)
INTEGER :: Skin_Temperature(2)    = (/ 0, 0 /)
INTEGER :: Cloud_Liquid_Water(2)  = (/ 0, 0 /)          
INTEGER :: Cloud_Top_Pressure(2)  = (/ 0, 0 /)
INTEGER :: Cloud_Fraction(2)      = (/ 0, 0 /)
INTEGER :: Surface_UWind(2)       = (/ 0, 0 /)
INTEGER :: Surface_VWind(2)       = (/ 0, 0 /)

NAMELIST / Retrieval / &
     Temperature,         &
     Humidity,            &
     Ozone,               &
     Surface_Temperature, &
     Surface_Humidity,    &
     Surface_Pressure,    &
     Skin_Temperature,    &
     Cloud_Liquid_Water,  &               
     Cloud_Top_Pressure,  &
     Cloud_Fraction,      &
     Surface_UWind,       &
     Surface_VWind
     

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

ErrorMessage(:) = ''
Num_ProfElementsUsed = 0

!--------------------------------------------------------
!1. Open the Retrieval File and Read in Header Data
!--------------------------------------------------------

CALL NWPSAF_OpenFile( Trim(FileName),      & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     fileunit )             ! out

! Read Namelist:

IF (FileUnit <= 0) THEN
  Errormessage(1) = ' Error opening Retrieval namelist file'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

Read( fileunit, nml=Retrieval, IOSTAT=ReadStatus) 

IF (ReadStatus /= 0) THEN
  Errormessage(1) = ' Error Reading Retrieval Namelist'
  Errormessage(2) = ' You may need to remove comments if compiling with F90'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

CLOSE(fileunit)

! Temperature profile

IF (ALL(Temperature(1:3) /= 0)) THEN
  ! Temperature(1) is the first temperature element to be retrieved (counting
  ! from the top of atmosphere).
  ! Temperature(2) is the number of elements to retrieve
  ! Temperature(3) is the number of the B-matrix element corresponding to
  ! Temperature(1)
  DO I=1, Temperature(2)
    Retrieved_Elements_Tmp(Num_ProfElementsUsed+I) = (Prof_FirstT-1) + &
          (Temperature(1)-1) + I 
    B_ElementsUsed_Tmp(Num_ProfElementsUsed+I) = Temperature(3) + I - 1
  END DO
  Num_ProfElementsUsed = Num_ProfElementsUsed + Temperature(2)
END IF

! Humidity Profile

IF (ALL(Humidity(1:3) /= 0)) THEN
  ! Humidity(1) is the first humidity element to be retrieved (counting
  ! from the top of atmosphere).
  ! Humidity(2) is the number of elements to retrieve
  ! Humidity(3) is the number of the B-matrix element corresponding to
  ! Humidity(1)
  DO I=1, Humidity(2)
    Retrieved_Elements_Tmp(Num_ProfElementsUsed+I) = (Prof_FirstQ-1) + &
          (Humidity(1)-1) + I 
    B_ElementsUsed_Tmp(Num_ProfElementsUsed+I) = Humidity(3) + I - 1
  END DO
  Ret_FirstQ = Num_ProfElementsUsed+1
  Num_ProfElementsUsed = Num_ProfElementsUsed + Humidity(2)
  Ret_LastQ = Num_ProfElementsUsed
  Num_WetLevels = Humidity(2)
END IF
 
! Ozone profile

IF (ALL(Ozone(1:3) /= 0)) THEN
  ! Ozone(1) is the first ozone element to be retrieved (counting from the 
  ! top of atmosphere).
  ! Ozone(2) is the number of elements to retrieve
  ! Ozone(3) is the number of the B-matrix element corresponding to
  ! Ozone(1)
  DO I=1, Ozone(2)
    Retrieved_Elements_Tmp(Num_ProfElementsUsed+I) = (Prof_FirstO3-1) + &
          (Ozone(1)-1) + I 
    B_ElementsUsed_Tmp(Num_ProfElementsUsed+I) = Ozone(3) + I - 1
  END DO
  Num_ProfElementsUsed = Num_ProfElementsUsed + Ozone(2)
END IF

! Surface Temperature
IF (ALL(Surface_Temperature(1:2) /= 0)) THEN
  ! Surface_Temperature(1) is non-zero if this element is to be retrieved.
  ! Surface_Temperature(2) is the number of the B-matrix element 
  ! corresponding to Surface_Temperature.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_T2 
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Surface_Temperature(2) 
END IF

! Surface Humidity
IF (ALL(Surface_Humidity(1:2) /= 0)) THEN
  ! Surface_Humidity(1) is non-zero if this element is to be retrieved.
  ! Surface_Humidity(2) is the number of the B-matrix element 
  ! corresponding to Surface_Humidity.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_q2 
  Ret_Q2 = Num_ProfElementsUsed
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Surface_Humidity(2) 
END IF

! Surface Pressure
IF (ALL(Surface_Pressure(1:2) /= 0)) THEN
  ! Surface_Pressure(1) is non-zero if this element is to be retrieved.
  ! Surface_Pressure(2) is the number of the B-matrix element 
  ! corresponding to Surface_Pressure.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_PStar 
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Surface_Pressure(2) 
END IF

! Skin Temperature
IF (ALL(Skin_Temperature(1:2) /= 0)) THEN
  ! Skin_Temperature(1) is non-zero if this element is to be retrieved.
  ! Skin_Temperature(2) is the number of the B-matrix element 
  ! corresponding to Skin_Temperature.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_TStar
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Skin_Temperature(2) 
END IF

! Cloud Liquid Water
IF (Cloud_Liquid_Water(1) /= 0) THEN
  ! Cloud_Liquid_Water(1) is non-zero if this element is to be retrieved.
  ! Cloud_Liquid_Water(2) is the number of the B-matrix element 
  ! corresponding to Cloud_Liquid_Water.
  MwClwRetrieval = .TRUE.
  IF (Retrieve_LWP) THEN
    Num_ProfElementsUsed = Num_ProfElementsUsed + 1
    Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_LWP 
    B_ElementsUsed_Tmp(Num_ProfElementsUsed) =  Cloud_Liquid_Water(2)
    Ret_CLW = Num_ProfElementsUsed
  END IF
 END IF

! Cloud_Top_Pressure (doesn't need a B element)
IF (Cloud_Top_Pressure(1) /= 0) THEN  
  ! Cloud_Top_Pressure(1) is non-zero if this element is to be retrieved.
  ! Cloud_Top_Pressure(2) is the number of the B-matrix element 
  ! corresponding to Cloud_Top_Pressure.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_CTP
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Cloud_Top_Pressure(2) 
  CloudyRetrieval = .TRUE.
  DetectCloud = .FALSE.
  Ret_CTP = Num_ProfElementsUsed
END IF

! Cloud Fraction (doesn't need a B element)
IF (Cloud_Fraction(1) /= 0) THEN
  ! Cloud_Fraction(1) is non-zero if this element is to be retrieved.
  ! Cloud_Fraction(2) is the number of the B-matrix element 
  ! corresponding to Cloud_Fraction.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_CloudCover
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Cloud_Fraction(2) 
  CloudyRetrieval = .TRUE.
  DetectCloud = .FALSE.
  Ret_CloudCover = Num_ProfElementsUsed
END IF

! Surface U-Wind
IF (ALL(Surface_UWind(1:2) /= 0)) THEN
  ! Surface_UWind(1) is non-zero if this element is to be retrieved.
  ! Surface_UWind(2) is the number of the B-matrix element 
  ! corresponding to Surface_UWind.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_uwind
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Surface_UWind(2) 
END IF

! Surface V-Wind
IF (ALL(Surface_VWind(1:2) /= 0)) THEN
  ! Surface_VWind(1) is non-zero if this element is to be retrieved.
  ! Surface_VWind(2) is the number of the B-matrix element 
  ! corresponding to Surface_VWind.
  Num_ProfElementsUsed = Num_ProfElementsUsed + 1
  Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_vwind 
  B_ElementsUsed_Tmp(Num_ProfElementsUsed) = Surface_VWind(2) 
END IF

! Add entry for Surface U-Wind if B-matrix element has not been added from
! B-matrix file; i.e. Surface_UWind(2) == 0
IF (Surface_UWind(1) /= 0 .AND. Surface_UWind(2)==0) THEN
  ! Surface_UWind(1) is non-zero if this element is to be retrieved.
  ! Surface_UWind(2) is the number of the B-matrix element 
  ! corresponding to UWind
    Num_ProfElementsUsed = Num_ProfElementsUsed + 1
    Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_Uwind
    B_ElementsUsed_Tmp(Num_ProfElementsUsed) =  Surface_UWind(2)
    Ret_UWind = Num_ProfElementsUsed
 END IF

! Add entry for Surface V-Wind if B-matrix element has not been added from
! B-matrix file; i.e. Surface_VWind(2) == 0 
IF (Surface_VWind(1) /= 0 .AND. Surface_VWind(2)==0) THEN
  ! Surface_VWind(1) is non-zero if this element is to be retrieved.
  ! Surface_VWind(2) is the number of the B-matrix element 
  ! corresponding to UWind
    Num_ProfElementsUsed = Num_ProfElementsUsed + 1
    Retrieved_Elements_Tmp(Num_ProfElementsUsed) = Prof_Vwind
    B_ElementsUsed_Tmp(Num_ProfElementsUsed) =  Surface_VWind(2)
    Ret_VWind = Num_ProfElementsUsed
END IF

! Allocate arrays and transfer contents of _Tmp variables.

! The array Retrieved_Elements indicates the  elements of the full 
! profile vector (contained in RTBack or RTGuess) that are to be 
! used in the retrieval step.  

ALLOCATE(Retrieved_Elements(Num_ProfElementsUsed))
Retrieved_Elements = Retrieved_Elements_Tmp(1:Num_ProfElementsUsed)

!  The B_ElementsUsed array defines the elements in the B matrix
!  that are actually retrieved.

ALLOCATE(B_ElementsUsed(Num_ProfElementsUsed))
B_ElementsUsed = B_ElementsUsed_Tmp(1:Num_ProfElementsUsed)

!----------------------------------------------------------------
! Check a few things.
!---------------------------------------------------------------

! Check that all elements of Retrieved_Elements and B_ElementsUsed 
! are valid

IF (ANY(Retrieved_Elements(:) <= 0 .OR. Retrieved_Elements(:) > ProfSize)) THEN
  WRITE(*,*) 'Retrieved_Elements(:) = ',Retrieved_Elements(:)
  Errormessage(1) = ' Error setting up the Retrieved_Elements array '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
ENDIF

IF (ANY(B_ElementsUsed(:) < 0)) THEN
  WRITE(*,*) 'B_ElementsUsed(:)=',B_ElementsUsed(:)
  Errormessage(1) = ' Error setting up the B_ElementsUsed array '
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
ENDIF


! Ensure that if one cloud parameter is retreived, so is the other

IF (ANY(Retrieved_Elements == Prof_CTP) .NEQV. &
    ANY(Retrieved_Elements == Prof_CloudCover)) THEN
  Errormessage(1) = ' Only one cloud parameter has been requested! '
  Errormessage(2) = ' - please choose neither or both'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
ENDIF

IF ( CloudyRetrieval .and. CalcRadiance ) THEN
  Errormessage(1) = ' For now, CloudyRetrieval mode requires observations '
  Errormessage(2) = ' in brightness temperature units'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

END SUBROUTINE NWPSAF_SetUpRetrievals
 
