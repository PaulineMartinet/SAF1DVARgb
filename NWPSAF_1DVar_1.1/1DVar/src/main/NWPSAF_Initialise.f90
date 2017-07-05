!+ Sets up the necessary coefficient arrays and global variables for
! NWPSAF_1DVar

Subroutine NWPSAF_Initialise ( &
     RT_Params ) ! inout

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
! Description: Initialisation for NWPSAF_1DVar. Handles IASI only for now.
!              Space is allocated for global arrays, and subroutines
!              are called to read in the data.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     02/06/99 Original code started from ATOVS_Initialise. A. Collard.
! 2.2     13/03/02 Removed (already commented out) references to 
!                  NWPSAF_GetBmatrix.           
!                  Set number of channels to MaxChanUsed (allows a subset 
!                  to be initialised).                          A. Collard.
! 2.3     27/05/02 Updated call to NWPSAF_InitBmatrix.            A. Collard.
! 3.0.4   02/03/04 Tidy up file naming.                         A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date      Comment
! ------  --------  -------
! 32      15/11/13  Add Output_Dir. P. Weston
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -----------------------------------------------------
!--

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
     ChannelSelection_Type

USE NWPSAFMod_Params, ONLY : &
     GeneralMode, &
     DebugMode,   &
     StatusFatal, &
     StatusOK,    &
     MaxChanUsed, &    
     Coeffs_Dir, &
     numfiles, &
     filenames, &
     FileType_Rmatrix, &
     FileType_Bmatrix, &
     FileType_ChannelChoice, &
     FileUnit_Retrieved_Profiles, &
     FileUnit_Retrieved_BTs, &
     FileUnit_MinimisationLog, &
     FileUnit_BTMinimisationLog, &
     FileUnit_AMatrix, &
     FileUnit_AmMatrix, &
     FileUnit_ProfileQC, &
     FileUnit_Jacobian, &
     FileUnit_RetJacobian, &
     FileUnit_AveragingKernel, &
     EnhancedDiagnostics, &
     Output_Dir

USE NWPSAFMod_RTmodel, ONLY :    &
     FastmodelMode_Initialise, &
     RTParams_Type, &  
     NeitherProfile

IMPLICIT NONE

INCLUDE 'NWPSAF_FreeUnit.interface'
INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_Channellist.interface'
INCLUDE 'NWPSAF_Fastmodel_Interface.interface'
INCLUDE 'NWPSAF_InitBmatrix.interface'
INCLUDE 'NWPSAF_InitRmatrix.interface'
INCLUDE 'NWPSAF_OpenFile.interface'

! Subroutine Arguments ::

TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Initialise"

! Local variables:
INTEGER :: unit_number      ! loop counter
INTEGER :: file             ! loop counter
INTEGER :: ErrorCode
INTEGER :: ErrStatRep       ! error status for NWPSAF_Report
INTEGER :: Fastmodel_Mode
INTEGER :: fileunit_Bmatrix        = 0 ! unit number for B-matrix file
INTEGER :: fileunit_Rmatrix        = 0 ! unit number for R-matrix file
INTEGER :: fileunit_ChannelChoice  = 0
CHARACTER(LEN=240):: file_pathname !complete pathname of file
CHARACTER(LEN=10) :: Access = "SEQUENTIAL"
CHARACTER(LEN=8)  :: Action = "READ"
CHARACTER(LEN=11) :: Form   = "FORMATTED"
CHARACTER(LEN=7)  :: Status = "OLD"
INTEGER :: fileunits(numfiles)            ! unit numbers given to the files 
CHARACTER(LEN=80) :: ErrorMessage(1)      ! Message for NWPSAF_Report
CHARACTER(LEN=80) :: Message(numfiles+1)  ! Message for NWPSAF_Report
TYPE(ChannelSelection_Type) :: UsedChans  ! Used for passing to
                                          ! Fastmodel_Interface only
INTEGER :: WhichProf

!-----------------------------------------------------------------------------

! -----------------------
! 1. Initialise variables
! -----------------------

ErrorMessage(:) = ' '
ErrStatRep      = StatusOK

! ------------------------------------------------------
! 2. Initialise control information and coefficient data
! ------------------------------------------------------

!2.1) Open coefficient files
!----
!All files, filetypes etc are declared in NWPSAFMod_Params. To use different
!files then minor modifications of this code and that module are needed.

!2.1.1) First do all files except for the RT coefficient files
!----

Message(:) = ' '

IF ( GeneralMode >= DebugMode ) THEN
  Message(1) = 'Opening coefficient files:'
  CALL NWPSAF_Report( RoutineName,Message(1:1) )
END IF

!Loop over files: construct pathname and call the file opening routine
get_files : DO file = 1, numfiles
  file_pathname = ' '
  WRITE(file_pathname,FMT=*) Trim(Coeffs_Dir),'/',Trim(filenames(file))
  CALL NWPSAF_OpenFile( Trim(AdjustL(file_pathname)), & ! in
      Access,              & ! in
      Action,              & ! in
      Status,              & ! in
      Form,                & ! inout
      fileunits(file) )      ! out
END DO get_files

!Assign unit variables with the correct units

fileunit_Bmatrix        = fileunits(FileType_Bmatrix) 
fileunit_Rmatrix        = fileunits(FileType_Rmatrix)
fileunit_ChannelChoice  = fileunits(FileType_ChannelChoice) 


!2.2) Open files for retrieval output
!----

file_pathname=Trim(Output_Dir)//'/Retrieved_Profiles.dat'
Action = "WRITE"
STATUS="UNKNOWN"
CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     FileUnit_Retrieved_Profiles )      ! out

file_pathname=Trim(Output_Dir)//'/Retrieved_BTs.dat'
Action = "WRITE"
STATUS="UNKNOWN"
CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     FileUnit_Retrieved_BTs )      ! out

!List of whether profiles resulted in a retrieval or not
file_pathname = Trim(Output_Dir)//'/ProfileQC.dat'
Action = "WRITE"
STATUS="UNKNOWN"
CALL NWPSAF_OpenFile( Trim(AdjustL(file_pathname)), & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     FileUnit_ProfileQC )      ! out

IF (GeneralMode >= DebugMode) THEN
   file_pathname=Trim(Output_Dir)//'/Minimisation.log'
   Action = "WRITE"
   STATUS="UNKNOWN"
   CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     FileUnit_MinimisationLog )      ! out

   file_pathname=Trim(Output_Dir)//'/Minimisation_BT.log'
   Action = "WRITE"
   STATUS="UNKNOWN"
   CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                     Access,              & ! in
                     Action,              & ! in
                     Status,              & ! in
                     Form,                & ! inout
                     FileUnit_BTMinimisationLog )      ! out

  IF ( EnhancedDiagnostics ) THEN

    ! This is where the analysis error covariance matrix
    ! is output
    file_pathname=Trim(Output_Dir)//'/A-Matrix.out'
    Action = "WRITE"
    STATUS="UNKNOWN"
    CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                      Access,              & ! in
                      Action,              & ! in
                      Status,              & ! in
                      Form,                & ! inout
                      FileUnit_AMatrix     ) ! out

    ! And this is where the measurement noise matrix
    ! is output
    file_pathname=Trim(Output_Dir)//'/Am-Matrix.out'
    Action = "WRITE"
    STATUS="UNKNOWN"
    CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                      Access,              & ! in
                      Action,              & ! in
                      Status,              & ! in
                      Form,                & ! inout
                      FileUnit_AmMatrix    ) ! out

    ! Jacobians of background profile
    file_pathname=Trim(Output_Dir)//'/BgJacobian.out'
    Action = "WRITE"
    STATUS="UNKNOWN"
    CALL NWPSAF_OpenFile( Trim(file_pathname), & ! in
                      Access,              & ! in
                      Action,              & ! in
                      Status,              & ! in
                      Form,                & ! inout
                      FileUnit_Jacobian    ) ! out

    ! Jacobians of retrieved profile
    file_pathname = Trim(Output_Dir)//'/RetJacobian.out'
      Action = "WRITE"
      STATUS="UNKNOWN"
      CALL NWPSAF_OpenFile( Trim(AdjustL(file_pathname)), & ! in
                        Access,              & ! in
                        Action,              & ! in
                        Status,              & ! in
                        Form,                & ! inout
                        FileUnit_RetJacobian    ) ! out

    ! Averaging Kernels from retrieval
    file_pathname = Trim(Output_Dir)//'/AveragingKernel.out'
      Action = "WRITE"
      STATUS="UNKNOWN"
      CALL NWPSAF_OpenFile( Trim(AdjustL(file_pathname)), & ! in
                        Access,                  & ! in
                        Action,                  & ! in
                        Status,                  & ! in
                        Form,                    & ! inout
                        FileUnit_AveragingKernel ) ! out
  END IF

END IF

!2.3) Construct channel arrays
!----

CALL NWPSAF_Channellist(FileUnit_ChannelChoice)

!2.4) Initialise Fastmodel
!----

IF ( GeneralMode >= DebugMode ) THEN
   Message(1) = 'Initialising RT Model'
   CALL NWPSAF_Report( RoutineName,Message(1:1) )
END IF
 
Fastmodel_Mode = FastModelMode_Initialise
WhichProf = NeitherProfile
UsedChans % NumChans = 0

CALL NWPSAF_Fastmodel_Interface( &
  Fastmodel_Mode,          & ! in
  RT_Params,               & ! inout
  WhichProf,               & ! in
  UsedChans,               & ! in
  ErrorCode)                 ! out

IF ( GeneralMode >= DebugMode ) THEN
   Message(1) = 'RT Model Initialised'
   CALL NWPSAF_Report( RoutineName,Message(1:1) )
END IF

!2.5) Read in biases
!----

! To be defined by user

!2.6) Set up R-matrix
!----
!
! Each satellite initialised for the RT model should have an R-Matrix 
!
!The current R-matrix file may have a greater number of elements than
!that required, but this does not affect the program.

IF ( fileunit_Rmatrix /= 0 ) THEN
   CALL NWPSAF_InitRmatrix( &
        fileunit_Rmatrix,                & ! in
        MaxChanUsed,                     & ! in
        RT_Params  )                       ! inout
ELSE
   ErrStatRep = StatusFatal
   ErrorMessage(1) = 'R-matrix not read in'
   CALL NWPSAF_Report( RoutineName,  &
        ErrorMessage, &
        ErrStatRep    )
END IF

!2.8) Set up B-matrix
!----
!The B-matrix is passed in through the NWPSAFMod_Params module
CALL NWPSAF_InitBmatrix( &
     fileunit_Bmatrix)        ! in

!----------
!3. Tidy up
!----------

close_files : DO file = numfiles, 1, -1
   unit_number = fileunits(file)
   IF (unit_number > 0) THEN
      CLOSE(unit_number)
      IF ( GeneralMode >= DebugMode ) THEN
         WRITE( UNIT=Message(file),FMT='(A,I4)') &
              'Closed file with unit number', unit_number
      END IF
      CALL NWPSAF_FreeUnit(unit_number,ErrStatRep)
   END IF
END DO close_files
IF ( GeneralMode >= DebugMode ) THEN
   CALL NWPSAF_Report( RoutineName,Message(1:numfiles) )
END IF


End Subroutine NWPSAF_Initialise
