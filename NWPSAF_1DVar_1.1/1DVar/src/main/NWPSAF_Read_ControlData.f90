SUBROUTINE NWPSAF_Read_ControlData ()

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
! Description: Reads file containing main control data for NWPSAF 1DVar
!
! Method:
!   Reads control data and makes sure the correct sections of the structures
!   defined in NWPSAFMod_Params and NWPSAFMod_RTModel are assigned.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.0      22/10/99  Original Version 
! 1.1      09/11/00  Added FirstOb and LastOb.  ADC 
! 1.2      10/01/01  Revised for standalone NWPSAF 1DVar.
! 1.3      01/03/01  Added ATOVS Roger Saunders
! 1.4      16/04/02  Fixed small bug where MaxCloudChans was mistyped 
!                    MaxRetChans.   ADC
! 2.3      28/05/02  Changed to Namelist version.   A. Collard.
! 3.0.1    15/07/03  Add cloud retrievals (from M. Szyndel). A. Collard.
! 3.0.4    02/03/04  Changes for more generic instrument treatment. A. Collard.
! 3.0.5    29/03/04  Remove SoundingType_Text.  Add in check for forcing
!                    use of Eqn_101.                         A. Collard.
! 3.1.0    05/04/05  Added support for RTTOV8.               E. Pavelin.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date     Comment
! 13      30/01/09 Added RTTOV9. E. Pavelin.
! 25      20/02/12 Added RTTOV10. P. Weston.
! 28      22/02/12 Added whether Lqtotal is 0 or 1 for cloud liquid water
!                  retrievals.
! 31      18/01/13 Added RTTOV11. P. Weston.
!
! Code Description:
!   Language:           Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Params, ONLY:       &
     StatusFatal,                 &
     StatusWarning,               &
     FirstOb, LastOb,             &
     Minimisation_Method,         &
     Newtonian,                   &
     Additional_Cost_Function,    &
     No_Additional_Cost_Function, &
     Perform1DVar,                &
     Force_Eqn_101,               &
     Control,                     &
     Retrieve_qtotal,             &
     Retrieve_LWP,                &
     Output_Dir,                  &
     Coeffs_Dir

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_OpenFile.interface'


! Declarations

CHARACTER(LEN=100) :: ErrorMessage(2)  ! Message for NWPSAF_Report
CHARACTER(LEN=*), PARAMETER :: Filename = 'ControlData.NL'
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NWPSAF_Read_ControlData'
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER (LEN=5) :: Status = "OLD"

INTEGER :: fileunit              ! I/O unit number
INTEGER :: ReadStatus
!INTEGER :: OpenStatus

!-----------------------------------------------------------------------------

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

ErrorMessage(:) = ''

!--------------------------------------------------------
!1. Set directory locations by capturing environment variables
!--------------------------------------------------------

CALL get_environment_variable("COEFFS_DIR", Coeffs_Dir)
CALL get_environment_variable("OUTPUT_DIR", Output_Dir)

!--------------------------------------------------------
!2. Open the Control File and Read in Header Data
!--------------------------------------------------------
 
CALL NWPSAF_OpenFile( Trim(FileName),      & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     fileunit )             ! out

! Read Namelist:

IF (FileUnit <= 0) THEN
  Errormessage(1) = 'Error opening Control Namelist'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 

END IF

Read( fileunit, nml=Control, IOSTAT=ReadStatus) 

IF (ReadStatus /= 0) THEN
  Errormessage(1) = ' Error Reading Control Namelist'
  Errormessage(2) = ' You may need to remove comments. Also make sure there is a blank line at the end of the file'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

CLOSE(fileunit)

!------------------------------------------------------
! 3. Check some of the variables from the namelist
!------------------------------------------------------

! First observation to be processed:
!------
IF (FirstOb < 0 ) THEN
  Errormessage(1) = ' Error Reading FirstOb'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 

END IF

! Last observation to be processed:
!------
IF (LastOb < FirstOb ) THEN
  Errormessage(1) = ' Error Reading LastOb'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
END IF

!------------------------------------------------------
! 3. Check some more of the variables from the namelist
!------------------------------------------------------

! Minimisation method:
!------
IF (Minimisation_Method < 0 .OR. Minimisation_Method > 2) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' Minimisation Method',Minimisation_Method,' is not allowed'
  CALL NWPSAF_Report( &
    RoutineName,            & ! in
    ErrorMessage,           & ! in
    ErrorStatus=StatusFatal ) ! in 
ELSE IF (Minimisation_Method == 0) THEN
  Errormessage(1) = &
       ' Minimisation Method 0 chosen - turning off Minimisation'
  Perform1DVar = .FALSE.
  CALL NWPSAF_Report( &
    RoutineName,              & ! in
    ErrorMessage,             & ! in
    ErrorStatus=StatusWarning ) ! in 
END IF

IF ( Force_Eqn_101 .AND. Minimisation_Method /= Newtonian ) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' If Force_Eqn_101 is .TRUE., Minimisation_Method must be 1'
  WRITE(Errormessage(1),FMT=*)  &
       ' Resetting Minimisation_Method'
  CALL NWPSAF_Report( &
    RoutineName,              & ! in
    ErrorMessage,             & ! in
    ErrorStatus=StatusWarning ) ! in 
  Minimisation_Method = Newtonian
END IF

IF ( Force_Eqn_101 .AND. &
     Additional_Cost_Function /= No_Additional_Cost_Function ) THEN
  WRITE(Errormessage(1),FMT=*)  &
       ' If Force_Eqn_101 is .TRUE., Additional_Cost_Function must be 1'
  WRITE(Errormessage(1),FMT=*)  &
       ' Resetting Additional_Cost_Function'
  CALL NWPSAF_Report( &
    RoutineName,              & ! in
    ErrorMessage,             & ! in
    ErrorStatus=StatusWarning ) ! in 

  Additional_Cost_Function = No_Additional_Cost_Function
END IF

! Check the type of Mw cloud liquid water retrieval
!------

Retrieve_LWP = .not. Retrieve_qtotal 

END SUBROUTINE NWPSAF_Read_ControlData
