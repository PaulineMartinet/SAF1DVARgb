!+ Sets up the observation (E) and forward model (F) combined 
!+ error matrix (E+F)

Subroutine NWPSAF_InitRmatrix( &
  file_unit,         & ! in
  total_channels,    & ! in
  RT_Params)           ! inout


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
! Description: Sets up the observation (E) and forward model (F) combined 
! error covariance
!              matrix R = E + F for each satellite
!
! The observational error covarince matrix can be defined as either 
!     1) A full matrix
!     2) A band diagonal matrix
!     3) A set of eigenvectors
! The structure element R_Reference(R_index)%Rtype defines which of these
! is being used for each satellite datatype.
!
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     23/09/99 Original code based on ATOVS_InitRMatrix.
!                  Changed R_matrices to a structure
! 1.1     28/02/00 For band diagonal case, the second index in the structure
!                  starts at zero, with the zeroth index referring to the
!                  diagonal.    A. Collard.  Met Office.
! 2.3     22/05/02 Changes to allow for better use of composite instruments.
!                                                       A. Collard.
! 3.0.5.  29/03/04 Change to allow default instrument.  A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    Default_SatID_Text, &
    R_Full_Matrix,      &
    R_Band_Diagonal,    &
    R_Eigenvectors,     &
    R_Reference

USE NWPSAFMod_RTModel, ONLY : &
    RTParams_Type

USE NWPSAFMod_Params, ONLY : &
    GeneralMode, &
    DebugMode,   &
    VerboseMode, &
    StatusFatal, &
    StatusOK,    &
    StatusWarning, &
    UsePCs

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'

! Subroutine arguments:
INTEGER, INTENT(IN)    :: file_unit  !I/O unit number
INTEGER, INTENT(IN)    :: total_channels
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params  ! Data used in 
                                           ! RT model interface

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_InitRmatrix"

! Local variables:
INTEGER :: numRchannels  ! no. of channels represented in the matrix
INTEGER :: numRelements  ! no. of channels, bands or eigenvalues
INTEGER :: Inverse     ! Variable set to unity if I/P matrix is inverted
INTEGER :: band_i      ! loop counter
INTEGER :: channel_i   ! loop counter
INTEGER :: channel_j   ! loop counter
INTEGER :: nchan       ! number of channels to put into R_reference
INTEGER :: readstatus  ! input error flag
INTEGER :: ErrStatRep  ! error status flag
INTEGER :: RType       ! Form of error data
INTEGER :: MatricesFound
INTEGER :: R_index
LOGICAL :: End_of_File = .FALSE.
LOGICAL :: GotMatrix(RT_Params % Num_SatIDs)
INTEGER, ALLOCATABLE :: Channel_number(:) ! array to hold data file channel #'s 
REAL,    ALLOCATABLE :: matrix(:,:)
REAL,    ALLOCATABLE :: eigenvalues(:)
CHARACTER(LEN=20) :: Instrument_Name
CHARACTER(LEN=80) :: Message(3)       ! message for NWPSAF_Report

!------------------------------------------------------------------------------

readstatus = 0
ErrStatRep = StatusOK
Message(:)      = ' '

MatricesFound = 0
GotMatrix(:) = .FALSE.

ALLOCATE(R_Reference(RT_Params % Num_SatIDs))
R_Reference(:) % Num_Chans = 0

Read_matrices: DO  ! read data until end of file reached 

  READ (file_unit,FMT='(A20)', IOSTAT=readstatus ) Instrument_Name
  IF (readstatus < 0) THEN
    End_of_File = .TRUE.
  ELSE IF ( readstatus > 0 ) THEN
    Message(1) = 'Error Reading Instrument Name'
    ErrStatRep = StatusFatal
    EXIT Read_matrices
  END IF

  IF (.NOT.(End_Of_File)) THEN
    READ (file_unit,FMT=*, IOSTAT=readstatus ) &
          Rtype, numRchannels, numRelements, Inverse 
    IF (readstatus < 0) THEN
      End_of_File = .TRUE.
    ELSE IF ( readstatus /= 0 ) THEN
      Message(1)='Error Reading Num_Channels etc.'
      ErrStatRep = StatusFatal
      EXIT Read_matrices
    END IF
  END IF

  IF ( Rtype == R_Full_Matrix .AND. numRchannels /= numRelements ) THEN
    Message(1) = 'Full Matrix is not Square'
    ErrStatRep = StatusFatal
    EXIT Read_matrices
  END IF
      
  IF (numRchannels * numRelements > 10000000) THEN
    Message(1) = 'R matrix is too big'
    WRITE( UNIT=Message(2),FMT='(A29,2I7)' ) &
           'numRchannels, numRelements =',numRchannels, numRelements
    ErrStatRep = StatusFatal
    EXIT Read_matrices
  END IF
      
  IF ( .NOT.(End_of_File) ) THEN
    ALLOCATE( Channel_number(numRchannels) )
    ! These are absolute instrument channel numbers and are only used
    ! to ensure the same channels are used in the observations and the
    ! R-matrix.
    READ (file_unit,*,IOSTAT=readstatus)  &
      ( Channel_number(channel_i),channel_i = 1,numRchannels )
    IF ( readstatus /= 0 ) THEN
      Message(1)= 'R matrix I/O error (second read)'
      ErrStatRep = StatusFatal
      EXIT Read_matrices
    ELSE
      ALLOCATE( matrix(numRchannels,numRelements) )
      Read_data : DO channel_i = 1, numRelements
        READ (file_unit,*,IOSTAT=readstatus) & 
          ( matrix(channel_j,channel_i),channel_j = 1,numRchannels )
        IF ( readstatus /= 0 ) THEN
          Message(1)='R matrix I/O error (third read) at element'
          WRITE( UNIT=Message(2),FMT='(A1,I4,A1,I4,A1)' ) &
            '(',channel_i,',', channel_j,')'
          ErrStatRep = StatusFatal
          EXIT Read_matrices
        ELSE IF ( Rtype == R_Eigenvectors ) THEN
          ALLOCATE( eigenvalues(numRelements) )
          READ (file_unit,'(10F8.4)',IOSTAT=readstatus) & 
            ( eigenvalues(band_i),band_i = 1,numRelements )
          IF ( readstatus /= 0 ) THEN
            Message(1) = 'R matrix I/O error (fourth read)'
            ErrStatRep = StatusFatal
            EXIT Read_matrices
          END IF
        END IF
      END DO Read_data
    END IF
      
    Check_satellite: DO R_index = 1, RT_Params % Num_SatIDs

      IF ( Trim(Instrument_Name) == &
           Trim(RT_Params % SatID(R_index) % SatID_Text) .OR. &
           Trim(RT_Params % SatID(R_index) % SatID_Text) == &
           Trim(Default_SatID_Text)) THEN
        IF (GeneralMode >= VerboseMode) &
          WRITE(*,*) 'Reading in R-matrix for ',Trim(Instrument_Name)
        IF (numRchannels /= RT_Params % SatID(R_index) % Num_Channels) THEN
          IF ( GeneralMode >= DebugMode ) THEN
            Message(1)='Incorrect number of elements for R-matrix in data file:'
            WRITE( UNIT=Message(2),FMT='(A,I5,A,I5)' ) &
              'Elements=',numRchannels,', Required=', &
              RT_Params % SatID(R_index) % Num_Channels
            CALL NWPSAF_Report(         &
              RoutineName,              &
              Message,                  &
              ErrorStatus=StatusWarning )
            Message(:) = ' '
          END IF
          EXIT Check_satellite   ! Don't use this satellite
        ELSEIF ( numRchannels < total_channels ) THEN
          Message(1) = 'Encountered R-matrix with too few elements:'
          WRITE( UNIT=Message(2),FMT='(A,I5,A,I5)' ) &
            'Elements=',numRchannels,', Required=',total_channels
          CALL NWPSAF_Report(         &
            RoutineName,              &
            Message,                  &
            ErrorStatus=StatusWarning )
          Message(:) = ' '
          EXIT Check_satellite   ! Don't use this satellite
        ELSEIF ( ANY(channel_number(:) &
                   /= RT_Params % Absolute_Channel_Number(:)) .AND. &
                 ANY(RT_Params % Absolute_Channel_Number(:) /= 0)   ) THEN
          Message(1) = 'Absolute channel numbers in R-matrix'
          Message(2) = 'are inconsistent with observations'
          CALL NWPSAF_Report(         &
            RoutineName,              &
             Message,                 &
            ErrorStatus=StatusWarning )
          Message(:) = ' '
          IF (GeneralMode >= VerboseMode) THEN
            WRITE(*,*) 'Absolute channel numbers in R-Matrix:'
            WRITE(*,*) channel_number(:)
            WRITE(*,*) 'Absolute channel numbers in ObsFile.dat:'
            WRITE(*,*) RT_Params % Absolute_Channel_Number(:)
          END IF
          EXIT Check_satellite   ! Don't use this satellite
        END IF
            
        IF ( .NOT. GotMatrix(R_index) ) THEN
          RT_Params % SatID(R_index) % R_Matrix_Present = .TRUE.
          MatricesFound = MatricesFound + 1
          GotMatrix(R_index) = .TRUE.
          R_Reference(R_index) % Rtype = Rtype
          nchan = Min ( numRchannels, total_channels )
          R_Reference(R_index) % Num_Chans = nchan
          ALLOCATE(R_Reference(R_index) % Diagonal(nchan))
          IF (Inverse == 1) THEN
            R_Reference(R_index) % Inverse = .TRUE.
          ELSE
            R_Reference(R_index) % Inverse = .FALSE.
          END IF
          SELECT CASE (RTYPE)
          CASE (R_Full_Matrix)      
            ALLOCATE(R_Reference(R_index) % Matrix(nchan,nchan))
            R_Reference(R_index) % Matrix(1:nchan,1:nchan) = &
                matrix(1:nchan,1:nchan)
            IF (Inverse == 1 .and. (UsePCs) ) THEN
              Message(1) = 'Need matrix diagonal for PC QC processing'
              Message(2) = 'Setting to 1000.0'
              CALL NWPSAF_Report(         &
                RoutineName,              &
                Message,                  &
                ErrorStatus=StatusWarning )
              Message(:)=''
              DO channel_i=1,nchan
                R_Reference(R_index) % Diagonal(channel_i) = 1000.0
              END DO
            ELSEIF (Inverse == 0) THEN 
              DO channel_i=1,nchan
                R_Reference(R_index) % Diagonal(channel_i) = &
                  R_Reference(R_index) % Matrix(channel_i,channel_i)
              END DO
            ENDIF
          CASE (R_Band_Diagonal)  
              
            ! numRelements is the number of bands in the banded matrix
            ! e.g., 1=diagonal, 2=tridiagonal etc.  In the 
            ! R_Reference structure 0=diagonal, 1=tridiagonal etc. 
            ! and Num_Elements is reduced accordingly
            R_Reference(R_index) % Num_Elements = numRelements - 1
            ALLOCATE(R_Reference(R_index) % Matrix(nchan,0:numRelements-1))
            R_Reference(R_index) % Matrix(1:nchan,0:numRelements-1) = &
               matrix(1:nchan,1:numRelements)
            R_Reference(R_index) % Diagonal(1:nchan) = &
               R_Reference(R_index) % Matrix(1:nchan,0)
          CASE (R_Eigenvectors) 
            ! numRelements is the number of eigenvectors
            R_Reference(R_index) % Num_Elements = numRelements
            ALLOCATE(R_Reference(R_index) % Matrix(nchan,numRelements))
            ALLOCATE(R_Reference(R_index) % Eigenvalues(numRelements))
            R_Reference(R_index) % Matrix(1:nchan,1:numRelements) = &
              matrix(1:nchan,1:numRelements)
            R_Reference(R_index) % Eigenvalues(:) = EigenValues(:)
            DEALLOCATE( eigenvalues )
            IF ( UsePCs ) THEN
              message(1) = 'Need matrix diagonal for PC QC processing'
              message(2) = 'Setting to 1000.0'
              CALL NWPSAF_Report(         &
                RoutineName,              &
                Message,                  &
                ErrorStatus=StatusWarning )
              Message(:)=''
              R_Reference(R_index) % Diagonal(1:nchan) = 1000.0
              END IF
          CASE DEFAULT
          END SELECT
          
          IF ( GeneralMode >= DebugMode ) THEN
            WRITE(UNIT=Message(1),FMT='(A,A,A)') &
                  'R matrix for ',Trim(Instrument_Name),' read in'               
          END IF
        END IF
            
      END IF
        
    END DO Check_satellite
      
    DEALLOCATE( Channel_number )
    DEALLOCATE( matrix )
      
  ELSE

    IF ( MatricesFound < RT_Params % Num_SatIDs ) THEN 
      ! didn't get all the data wanted
      Message(1) = 'EOF before all required data read in'
      Message(2) = 'satellites with no R-matrix will have reduced processing'
      ErrStatRep = StatusWarning
      !Reject satellites with no matrix for 1DVar processing
      WHERE(.NOT. GotMatrix(:) )
          RT_Params % SatID(:) % R_Matrix_Present = .FALSE.
      ENDWHERE
    END IF
    EXIT Read_matrices
      
  END IF
  
END DO Read_matrices


IF ( ErrStatRep /= StatusOK ) THEN
  CALL NWPSAF_Report(         &
    RoutineName,              &
    Message,                  &
    ErrorStatus=ErrStatRep )
ELSE
  IF ( GeneralMode >= DebugMode ) THEN
    WRITE(*,*) 'R-matrices initialised ok'
  END IF
END IF

End Subroutine NWPSAF_InitRmatrix
