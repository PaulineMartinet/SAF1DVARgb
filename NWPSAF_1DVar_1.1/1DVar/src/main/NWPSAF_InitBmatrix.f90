!+ Sets up the background error covariance matrix (B matrix)

Subroutine NWPSAF_InitBmatrix( &
  file_unit)         ! in

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
! Description: Sets up the background error covariance matrix (B matrix).
!  
! Method: Allocate space.
!         Read in each version of the matrix from the given file.
!         Remove unwanted elements
!         Copy and make up the cloud detection versions if required.
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     11/11/99 Original version from ATOVS_InitBmatrix.      A.Collard.
! 1.1     22/02/00 The gaps left in B_Reference and B_Reference_Inverse where 
!                  elements are not to be used are now removed entirely.
!                                                  A. Collard.  Met Office
! 2.3     27/05/02 Remove some unwanted code.      A. collard.  Met Office
! 3.1     21/07/03 Changes for cloudy retrievals.  A. Collard.
! 3.0.6   18/06/04 Set Binv to 1.e-10 not 0 for cloud.         A. Collard.
! 3.1.1   05/10/04 Allocated BFromFile size changed to Num_Elements (right)
!                  from Num_ProfElementsUsed (wrong).          A. Collard. 
!
! Ticket  Date     Comment
! ------- -------- -------
! 28      20/02/12 Added Bmatrix diagonal elements for LWP retrieval 
!                  TR Sreerekha
! 32      07/12/12 Added options to set diagonal elements for wind speed 
!                  retrieval internally if not defined from Bmatrix file.
!                  Get LwpSD and WindspeedSD from NWPSAFMod_Params
!                  A. Andersson (CM SAF)
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    B_Reference, &
    B_Reference_Inverse, &
    BmatrixTypes

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,        &
    DebugMode,          &
    StatusFatal,        &
    StatusWarning,      &
    StatusOK,           &
    B_ElementsUsed,     & 
    Retrieved_Elements, &
    CloudyRetrieval,    &
    WindspeedSD,        &
    LwpSD

USE NWPSAFMod_RTmodel, ONLY : &
    Num_ProfElementsUsed, &
    Prof_CTP,             &
    Prof_CloudCover,      &
    Prof_LWP,             &
    Prof_Uwind,           &
    Prof_Vwind

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_SatMatInv.interface'

! Subroutine arguments:

INTEGER, INTENT(IN) :: file_unit        ! I/O unit number
 
! Local constants:
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_InitBmatrix"

! Local variables:
INTEGER :: i,j            ! loop counters
INTEGER :: type           ! Code for surface type
INTEGER :: Num_elements   ! number of B matrix elements
INTEGER :: ReadStatus     ! IOSTAT keyword value
INTEGER :: ErrStatRep     ! error status

CHARACTER (LEN=80) :: BmatrixTitle
CHARACTER (LEN=80) :: Message(1)
CHARACTER (LEN=80) :: ErrorMessage(2)

REAL, ALLOCATABLE  :: Bfromfile(:,:)
REAL, ALLOCATABLE  :: B_Comp(:,:)


!-----------------------------------------------------------------------------


!1) Allocate and initialize data
!--

ErrStatRep = StatusOK
ErrorMessage(:) = ' '
Message(:) = ' '

ALLOCATE( B_Reference(Num_ProfElementsUsed,Num_ProfElementsUsed,BmatrixTypes) )
B_Reference(:,:,:)         = 0.0

ReadAllB : DO type = 1, BmatrixTypes  ! currently type 1=sea,2=land

  !2) Read in B-matrix data from file
  !--
  READ (UNIT=file_unit, FMT='(A80)', IOSTAT=ReadStatus) BmatrixTitle
  ! This is an extra title line compared to current ATOVS processing:
  READ (UNIT=file_unit, FMT='(A80)', IOSTAT=ReadStatus) BmatrixTitle
  IF ( ReadStatus < 0 ) THEN
    ErrStatRep = StatusWarning
    ErrorMessage(1) = 'Insufficient data to make all required B-matrices'
    EXIT
  END IF
  READ (UNIT=file_unit,FMT=*, IOSTAT=ReadStatus) Num_elements
  ALLOCATE( Bfromfile(Num_Elements,Num_Elements) )
  ReadData : DO i = 1, Num_elements
    READ (UNIT=file_unit, FMT='(5E16.8)', IOSTAT=ReadStatus ) & 
      (Bfromfile(j,i), j = 1,Num_elements)
    IF (ReadStatus /= 0) THEN
      Errormessage(1) = 'B matrix I/O error reading in data'
      ErrStatRep = StatusFatal
      EXIT ReadAllB
    END IF
  END DO ReadData
    
     
  !3) Transfer matrix into 1DVar arrays
  !--

  ! First check that all the required retrieval elements are in the B-matrix
  ! If this is not the case, trigger a fatal error and advise that BMatrix and
  ! Retrieval.in are incompatible.
   
  IF (ANY(B_ElementsUsed(:) > Num_Elements)) THEN
    ErrorMessage(1) = 'B Matrix too small - check that Retrieval.NL and'
    ErrorMessage(2) = 'Bmatrix are compatible!'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusFatal )
  END IF

  !
  ! Now load the required elements from Bfromfile into B_Reference
  !
  DO i = 1, Num_ProfElementsUsed
    ! Cloudy a priori errors are set to moderate values here so the inverse
    ! that follows is numerically stable.  They get reset to large values
    ! after the matrix inversion.
    IF ( CloudyRetrieval .AND. B_ElementsUsed(i) == 0 .AND. &
        ( Retrieved_Elements(i) == Prof_CTP .OR. &
          Retrieved_Elements(i) == Prof_CloudCover ) ) THEN
      B_Reference(:,i,type) = 0.0
      B_Reference(i,:,type) = 0.0
      B_Reference(i,i,type) = 1.0
    ELSE IF ( B_ElementsUsed(i) == 0. .AND. &
              Retrieved_Elements(i) == Prof_LWP) THEN
      B_Reference(:,i,type) = 0.0
      B_Reference(i,:,type) = 0.0
      B_Reference(i,i,type) = LwpSD * LwpSD
    ELSE IF ( B_ElementsUsed(i) == 0. .AND. &
              (Retrieved_Elements(i) == Prof_Uwind .OR. &
               Retrieved_Elements(i) == Prof_Vwind)) THEN
      B_Reference(:,i,type) = 0.0
      B_Reference(i,:,type) = 0.0
      B_Reference(i,i,type) = WindspeedSD * WindspeedSD
    ELSE
      WHERE(B_ElementsUsed(:) > 0) &
        B_Reference(i,:,type) = Bfromfile(B_ElementsUsed(i),B_ElementsUsed(:))
    END IF
  END DO



  DEALLOCATE(Bfromfile)
   
END DO ReadAllB

!5) Invert B matrix 
!--
IF ( ErrStatRep == StatusOK ) THEN
  ALLOCATE( B_Reference_Inverse(Num_ProfElementsUsed, &
                                Num_ProfElementsUsed, BmatrixTypes) )
  B_Reference_Inverse(:,:,:) = 0.
  DO type = 1, BmatrixTypes
    ALLOCATE( B_Comp(Num_ProfElementsUsed,Num_ProfElementsUsed) )
    B_Comp(:,:) = B_Reference(:,:,type)
    ! invert
    CALL NWPSAF_SatMatInv(   &
      Num_ProfElementsUsed,  & ! in
      Num_ProfElementsUsed , & ! in
      B_Comp)                 ! inout
      
    ! store
    B_Reference_Inverse(:,:,type) = B_Comp(:,:)
    DEALLOCATE( B_Comp )

    ! Change the cloudy retrieval errors to large (small in the inverse)
    ! values.  This has an effect of the cost function gradient test in  
    ! NWPSAF_Minimize.

    IF (CloudyRetrieval) THEN
      DO i = 1, Num_ProfElementsUsed
        IF ( B_ElementsUsed(i) == 0 .AND. &
             ( Retrieved_Elements(i) == Prof_CTP .OR. &
               Retrieved_Elements(i) == Prof_CloudCover ) ) THEN
          B_Reference(i,i,type) = 1.0e10
          B_Reference_Inverse(I,I,type) = 1.0e-10
        END IF
      END DO
    END IF
      
  END DO
END IF

!6) Finished
!--

! We can now deallocate the B_ElementsUsed arrays
DEALLOCATE(B_ElementsUsed)

IF (ErrStatRep /= StatusOK) THEN        ! will stop program if status fatal
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=ErrStatRep )
ELSEIF ( GeneralMode >= DebugMode ) THEN
  Message(1) = 'B matrices initialised ok'
  CALL NWPSAF_Report( &
    RoutineName, &
    Message(1:1) )
END IF

End Subroutine NWPSAF_InitBmatrix
 
