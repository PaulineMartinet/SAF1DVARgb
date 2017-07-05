!+ Cut down the Rmatrix to the required channels.

Subroutine NWPSAF_RMatrix_ChanSelect( &
     UsedChans, &  ! in
     R_Matrix,  &  ! in
     R_SubMatrix)  ! out

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
! ALLOCATE and setup R_SubMatrix - the observational + forward 
! model error covariance matrix after channel selection.
!------------------------------------------------------------------------
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     08/06/00 Original code taken from NWPSAF_Minimize.
!                                                 A. Collard.
! 1.1     21/11/00 Fixed bug in off-diagonal part of band-diagonal code.
!                                                 A. Collard.
! 2.2     03/05/02 Remove references to R_SubMatrix % Channel_Number.  
!                                                 A. Collard
! 3.0.6   23/06/04 Add check for input matrix being inverted. A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type,   &
    R_Full_Matrix,   &
    R_Band_Diagonal, & 
    R_Eigenvectors

USE NWPSAFMod_Params, ONLY : &
    StatusFatal,     &
    StatusWarning,   &
    UsePCs

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'

! Subroutine arguments:
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans ! channels to use
TYPE(R_Matrix_Type), INTENT(IN)  :: R_Matrix     ! R matrix
TYPE(R_Matrix_Type), INTENT(OUT) :: R_SubMatrix  ! Channel selected R matrix

! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RMatrix_ChanSelect"

! Local variables:
CHARACTER(LEN=80)  :: ErrorMessage(2)  ! Message for NWPSAF_Report

INTEGER I, J, Offset

!------------------------------------------------------------------------
! 0.1. Initialise variables.
!-----------------------------------------------------------------------

ErrorMessage(:) = ' '

!------------------------------------------------------------------------
! 1.0. Set up R_SubMatrix
!-----------------------------------------------------------------------

SELECT CASE (R_Matrix % RType)

  !---------------------------------------------------------------------
  ! 1.1.1. Channel selection from a full matrix is straightforward.
  !---------------------------------------------------------------------
  CASE (R_Full_Matrix)

    ! Check for cases where the matrix is input as an inverse and
    ! channel selection is not appropriate.
    IF (R_Matrix%Inverse .AND. &
        UsedChans % NumChans /= R_SubMatrix % Num_Chans) THEN
      Errormessage(1) = &
        'If channel selection is being performed the input R-matrix '
      Errormessage(2) = &
        'cannot be given as an inverse for the full matrix case'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF

    ! Initialise R_SubMatrix 
    ALLOCATE(R_SubMatrix % Matrix(UsedChans % NumChans,UsedChans % NumChans))
    ALLOCATE(R_SubMatrix % Diagonal(UsedChans % NumChans))
    R_SubMatrix % Rtype          = R_Matrix % RType
    R_SubMatrix % Num_Chans      = UsedChans % NumChans
    R_SubMatrix % Num_Elements   = UsedChans % NumChans
    R_SubMatrix % Inverse        = .FALSE.
    NULLIFY(R_SubMatrix % EigenValues)
    
    ! Do channel selection
    DO I=1, UsedChans % NumChans
      R_SubMatrix % Matrix(I,:) = &
        R_Matrix%Matrix (UsedChans % Channels(I), &
                         UsedChans % Channels(1 : UsedChans % NumChans) )
      R_SubMatrix % Diagonal(I) = R_Matrix%Diagonal(UsedChans % Channels(I))
    END DO

    !---------------------------------------------------------------------
    ! 1.1.2. Channel selection from a band diagonal matrix is complicated.
    !    Note that the MATRIX array is transposed in this step - this is 
    !    because some of the later routines expect the band-diagonal 
    !    representation to be this way round.
    !---------------------------------------------------------------------

  CASE (R_Band_Diagonal) 

    ! Check for cases where the matrix is input as an inverse and
    ! channel selection is not appropriate.
    IF (R_Matrix % Inverse .AND. R_Matrix % Num_Elements /= 0 .AND. &
        UsedChans % NumChans /= R_SubMatrix % Num_Chans) THEN
      Errormessage(1) = &
        'If channel selection is being performed the input R-matrix '
      Errormessage(2) = &
        'cannot be given as an inverse for the band-diagonal matrix case'
      CALL NWPSAF_Report( &
        RoutineName,            & ! in
        ErrorMessage,           & ! in
        ErrorStatus=StatusFatal ) ! in 
    END IF

    ! Initialise R_SubMatrix 
    ALLOCATE(R_SubMatrix % Matrix(0:R_Matrix % Num_Elements,UsedChans % NumChans))
    ALLOCATE(R_SubMatrix % Diagonal(UsedChans % NumChans))
    R_SubMatrix % Matrix(:,:) = 0.
    R_SubMatrix % Rtype          = R_Matrix % RType
    R_SubMatrix % Num_Chans      = UsedChans % NumChans
    R_SubMatrix % Num_Elements   = R_Matrix % Num_Elements
    R_SubMatrix % Inverse        = R_Matrix % Inverse
    NULLIFY(R_SubMatrix % EigenValues)

    ! Do channel selection ...
    ! ... Start by getting the diagonal elements for each channel
    R_SubMatrix % Matrix(0,:) = &
      R_Matrix%Matrix(UsedChans % Channels(1 : UsedChans % NumChans),0)
    R_SubMatrix % Diagonal(:) = R_Matrix%Diagonal(UsedChans % Channels(:))
    
    ! Now do the off-diagonal elements 
    IF (R_Matrix%Num_Elements >= 1) THEN
      DO I = 1, UsedChans % NumChans - 1
        Offset = 1
        J = UsedChans % Channels(I+Offset) - UsedChans % Channels(I)
        DO WHILE (J <= UsedChans % NumChans - 1 .AND. &
                  J <= R_Matrix % Num_Elements   .AND. & 
                  I+Offset < UsedChans % NumChans )
          R_SubMatrix % Matrix(Offset,I) = &
             R_Matrix % Matrix(UsedChans % Channels(I),J)
          Offset = Offset + 1
          J = UsedChans % Channels(I+Offset) - UsedChans % Channels(I)
        END DO
      END DO
    END IF

    !----------------------------------------------------------------------
    ! 1.1.3. Channel selection for the eigenvector case (although this is 
    !        not advised).
    !----------------------------------------------------------------------
  CASE (R_EigenVectors) 

    ! Initialise R_SubMatrix 
    ALLOCATE(R_SubMatrix % Matrix(UsedChans % NumChans,R_Matrix % Num_Elements))
    ALLOCATE(R_SubMatrix % Diagonal(UsedChans % NumChans))
    R_SubMatrix % Rtype          = R_Matrix % RType
    R_SubMatrix % Num_Chans      = UsedChans % NumChans
    R_SubMatrix % Num_Elements   = UsedChans % NumChans
    R_SubMatrix % Inverse        = R_Matrix % Inverse
    R_SubMatrix % EigenValues => R_Matrix % EigenValues

    ! DO channel selection, but write a warning message
    Errormessage(1) = 'Channel selection is inadvisable when the'
    Errormessage(2) = 'R-matrix is specified as eigenvectors'
    CALL NWPSAF_Report( &
      RoutineName,              & ! in
      ErrorMessage,             & ! in
      ErrorStatus=StatusWarning ) ! in 
    R_SubMatrix % Matrix(:,:) = &
      R_Matrix%Matrix(UsedChans % Channels(1 : UsedChans % NumChans),:)
    IF ( UsePCs ) THEN
      Errormessage(1) = 'Need matrix diagonal for PC QC processing'
      Errormessage(2) = 'Setting to 1000.0'
    CALL NWPSAF_Report( &
      RoutineName,              & ! in
      ErrorMessage,             & ! in
      ErrorStatus=StatusWarning ) ! in 
      R_Submatrix % Diagonal(:) = 1000.0
    END IF
    !---------------------------------------------------------------------
    ! 1.1.4. There are no other formats for observation errors.
    !---------------------------------------------------------------------

  CASE DEFAULT 
    Errormessage(1) = 'Unknown Format for Observation Errors'
    Errormessage(2) = ''
    CALL NWPSAF_Report( &
      RoutineName,            & ! in
      ErrorMessage,           & ! in
      ErrorStatus=StatusFatal ) ! in 
END SELECT
   
END Subroutine NWPSAF_RMatrix_ChanSelect
