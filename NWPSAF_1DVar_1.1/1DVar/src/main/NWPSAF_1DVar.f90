!+ 1DVar satellite sounding for IASI, main routine.

SUBROUTINE NWPSAF_1DVar ( &
     Obs,         & !inout
     Background,  & !inout
     RT_Params)     !inout

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
! Description: main routine for 1DVar processing of satellite sounding
!              data. Currently set up to process IASI.
!
! History:
!
! Version Date     Comment
! ------- -------  -------
! 1.1     26/5/99  Initial version adapted from ATOVS_1DVar.  A.D. Collard
! 2.0     25/5/00  Many undocumented changes in the last year.  ADC
! 2.1     25/5/00  Added in code to DEALLOCATE and NULLIFY BackChans and
!                  DetectCloudChans.   ADC
! 2.2     9/11/00  FirstOb and LastOb now passed through module (after
!                  being read in from ControlData.in) rather than being 
!                  defined here.  ADC.
! 2.3     1/03/01  Added ATOVS/ozone  Roger Saunders
! 2.3     28/5/02  Remove call to NWPSAF_CheckInput.    Andrew Collard.
! 3.0.1   15/07/03 Add cloud retrieval sounding type (M Szyndel) A. Collard.
! 3.0.4   02/03/04 Modifications to make instrument data more generic and to
!                  fix LastOb specification.                   A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type, &
    ModelOb_Type

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,    &
    DiagnosticMode, &
    StatusFatal,    &
    StatusWarning,  &
    FirstOb,        &
    LastOb,         &
    UsePCs,         &
    CalcRadiance

USE NWPSAFMod_RTmodel, ONLY :    &
    RTParams_Type

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_Initialise.interface'
INCLUDE 'NWPSAF_ProcessData.interface'

! Subroutine arguments:
TYPE(Ob_type),     INTENT(INOUT)  :: Obs        ! Observed/Retrieval data 
TYPE(ModelOb_type),INTENT(INOUT)  :: Background ! Background data
TYPE(RTParams_Type),INTENT(INOUT) :: RT_Params

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_1DVar"

! Local variables:
INTEGER :: Total_Channels
CHARACTER(LEN=80) :: Message(2) ! Message for NWPSAF_Report

!- End of header for NWPSAF_1Dvar_routine ---------------------------------------


!-----------------------
!1. Initialise variables
!-----------------------

Message(:)      = ' '

!--------------------------------------
! 2. Set data type and initialise data
!--------------------------------------

FirstOb = MAX(1,FirstOb)
! LastOb=0 zero is its default (use all) value.
IF (LastOb == 0) LastOb = Obs% header% NumObsLocal

IF (LastOb > Obs% header% NumObsLocal) THEN
   Message(1) = 'Resetting LastOb to be less than the number'
   Message(2) = 'of obs in the Obsfile'
   CALL NWPSAF_Report( RoutineName, & ! in
     Message(:),                    & ! in
     ErrorStatus= StatusWarning     ) ! in
   LastOb  = Obs% header% NumObsLocal
   Message(2) = ''
END IF

IF (LastOb < FirstOb) THEN
   Message(1) = 'LastOb is less than FirstOb'
   Message(2) = ''
   CALL NWPSAF_Report( RoutineName, & ! in
     Message(:),                    & ! in
     ErrorStatus= StatusFatal       ) ! in
END IF

!2.1) Read coefficient data, initialise global variables
!----

CALL NWPSAF_Initialise ( &
  RT_Params ) ! in

!---------------
!3. Process data
!---------------
!The bulk of all processing is done here. 

IF ( UsePCs ) THEN
  Total_Channels = Obs% header% PCScore% NumLev
ELSE IF ( CalcRadiance ) THEN
  Total_Channels = Obs% header% Radiance% NumLev
ELSE
  Total_Channels = Obs% header% BriTemp% NumLev
END IF
  
CALL NWPSAF_ProcessData( &
  Obs,            & ! inout
  Background,     & ! in
  FirstOb,        & ! in
  LastOb,         & ! in
  RT_Params,      & ! inout
  Total_channels  ) ! in

!----------------------------------------
!4. Final Output Message
!----------------------------------------

IF ( GeneralMode >= DiagnosticMode ) THEN
  Message(1) = '**** Finished 1DVar ****'
  CALL NWPSAF_Report( &
    RoutineName, &
    Message(:)   )
END IF

END SUBROUTINE NWPSAF_1DVar
