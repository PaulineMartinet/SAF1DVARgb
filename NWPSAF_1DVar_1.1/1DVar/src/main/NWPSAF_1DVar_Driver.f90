PROGRAM NWPSAF_1DVar_Driver

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
! Description: Sets up inputs for the NWPSAF_1DVar subroutine.
!
! Method:
!
!  Reads in control data, background data and observations and allocates 
!  various arrays to be used in the NWPSAF_1DVar code.
!
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.0      10/01/01  Original version.   Based on ATOVS and IASI Ops code.
!                                                  A. Collard.  Met Office.
! 1.1      01/03/01  Add ATOVS/ozone Roger Saunders
! 2.3      27/05/02  Remove unwanted structure elements in Observation and
!                    background structures.              A. Collard.
! 3.0.1    15/07/03  Add MSG data and initialise cloud.  A. Collard.
! 3.0.4    02/03/04  Modifications to make instrument data more generic.
!                                                        A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Constants, ONLY : &
     MissData_R, &
     MissData_I

USE NWPSAFMod_RTModel, ONLY : &
     RTParams_Type, &
     Num_RTLevels

USE NWPSAFMod_ObsInfo, ONLY : &
     ElementHeader_type,    &
     Element_type,          &
     ModelOB_type,          &
     OB_type

IMPLICIT NONE

INCLUDE 'NWPSAF_1DVar.interface'
INCLUDE 'NWPSAF_Read_ControlData.interface'
INCLUDE 'NWPSAF_Read_Observations.interface'
INCLUDE 'NWPSAF_Read_Background.interface'
INCLUDE 'NWPSAF_SetUpRetrievals.interface'
INCLUDE 'NWPSAF_DeAllocate.interface'

! Subroutine arguments:

! Declarations

TYPE(ElementHeader_type) :: Dummy_header  ! holds initialization data
TYPE(Element_type)       :: Dummy_element ! holds initialization data

TYPE(OB_type), TARGET :: Observations     ! Observation structure
TYPE(ModelOB_type)    :: BackGrModelOB    ! BackGround at Ob
TYPE(RTParams_type)   :: RT_Params        ! Info for the RT Model

!-----------------------------------------------------------------------------

!-------------------------------------
!0. Initialise constants and variables
!-------------------------------------

CALL NWPSAF_Read_ControlData () 

! Read in Observational and Background Data

CALL NWPSAF_Read_Observations ( &
     0,             & ! in
     BackGrModelOb, & ! out
     Observations,  & ! inout
     RT_Params)       ! inout

CALL NWPSAF_Read_Background ( &
     Observations% header% NumObsLocal, & !in
     Num_RTLevels,                      & !out
     BackGrModelOB )                      !inout

! Set up the parameters that are to be retrieved.

CALL NWPSAF_SetupRetrievals () ! no arguments

!-----------------------------------
!1. Initialize observation structure
!-----------------------------------

!--------------------------------------------
!1.1 Allocate space for surface and meta data
!--------------------------------------------

Dummy_header % NumLev = 0

Observations % header % NIter          = Dummy_header
Observations % header % CldCost        = Dummy_header
Observations % header % t2             = Dummy_header
Observations % header % rh2            = Dummy_header
Observations % header % TSkin          = Dummy_header
Observations % header % PStar          = Dummy_header
Observations % header % Jcost          = Dummy_header
Observations % header % Jcost_Gradient = Dummy_header

!Surface data
!------------

Dummy_header % NumLev = 1
Observations % header % TSkin   = Dummy_header
Observations % header % CTP     = Dummy_header
Observations % header % CldFrac = Dummy_header

!-------------------------------------------------------
!1.2 Allocate space for multi-level data and set headers
!-------------------------------------------------------

!T
!-
Dummy_header% NumLev = Num_RTLevels
ALLOCATE( Observations% t(Dummy_header% NumLev))
Observations% header% t = Dummy_header

!rh and CLW
!----------

Dummy_header% NumLev = Num_RTLevels
ALLOCATE( Observations% rh(Dummy_header% NumLev))
Observations% header% rh = Dummy_header

!Ozone
!-
Dummy_header% NumLev = Num_RTLevels
ALLOCATE( Observations% ozone(Dummy_header% NumLev))
Observations% header% ozone = Dummy_header

!-------------------
!1.3 Initialize data
!-------------------
!All fields are given initial values of missing data unless they are
!flags or PGE values. Most flags are initialized to zero.
!Note that in 1DVar, multi-level element flags are set according to a
!fixed pattern that depends on the channel selection for each
!observation. This is to ensure that some of these elements do not pass
!quality control. Here they are initialized to be rejected.

Observations % NIter           = MissData_I
Observations % CldCost         = MissData_R
Observations % coloz           = MissData_R

Dummy_element % Value    = MissData_R

Observations% t2      = Dummy_element
Observations% rh2     = Dummy_element
Observations% TSkin   = Dummy_element
Observations% CTP     = Dummy_element
Observations% CldFrac = Dummy_element

Observations% Jcost          = MissData_R
Observations% Jcost_Gradient = MissData_R

Observations% t(:)   = Dummy_element
Observations% rh(:)  = Dummy_element
Observations% ozone(:)  = Dummy_element

!----------------
!2. Perform 1Dvar
!----------------

CALL NWPSAF_1DVar(Observations,     & ! inout
                  BackGrModelOB,    & ! inout
                  RT_Params)          ! inout

!---------------------
!3. DEALLOCATE Arrays
!---------------------

CALL NWPSAF_DeAllocate ( &
     Observations,     &  ! inout
     BackGrModelOb,     & ! inout
     RT_Params)           ! inout
 
END PROGRAM NWPSAF_1DVar_Driver
