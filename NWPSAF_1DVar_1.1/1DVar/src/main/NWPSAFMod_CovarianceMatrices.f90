!+ Background and Observation Error Covariance Matrices.

Module NWPSAFMod_CovarianceMatrices

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
! Description: contains background and observation error covariance matrices.
!
! Owner: Manager of operational observation processing and assimilation
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     07/02/01 Original version, split out from now-defunct 
!                  AIRSMod_Variables.       A. Collard  
! 1.3     15/08/01 Added Bmatrix with type Bmatrix_type (from ATOVS code)
!                                                        Andrew Collard.
! 2.2     14/03/02 Removed unused Bmatrix structures.  
!                  Removed Channel_Number from R_Matrix_Type.  A. Collard.
! 2.3     27/05/02 Removed unwanted Bmatrix variables.  A. Collard.
! 3.0.5.  31/03/04 Changed BmatrixTypes from 1 to 2 to be consistent with
!                  NWPSAF_ProcessData.                    A. Collard.
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
!
! End of header -------------------------------------------------------

!--------
!B-matrix (background error covariance matrix)
!--------

INTEGER, PARAMETER :: BmatrixTypes = 2

!Indices for B-matrix type
INTEGER, PARAMETER :: Bmatrix_sea = 1         ! Code for sea
INTEGER, PARAMETER :: Bmatrix_land = 2        ! Code for land and sea ice

!B
REAL, ALLOCATABLE, TARGET, SAVE :: B_Reference(:,:,:)
REAL, ALLOCATABLE, TARGET, SAVE :: B_Reference_Inverse(:,:,:)


!--------
! R-matrix (observation error covariance matrix)
!--------

INTEGER, PARAMETER :: Max_Elements = 20 ! Number of off-diagonals to store
                                        ! when taking a band-diagonal matrix
                                        ! inverse.

TYPE R_Matrix_Type
   INTEGER          :: Rtype
   REAL, POINTER    :: Matrix(:,:)
   REAL, POINTER    :: Diagonal(:)
   REAL, POINTER    :: EigenValues(:)
   INTEGER          :: Num_Chans
   INTEGER          :: Num_Elements
   LOGICAL          :: Inverse
End TYPE R_Matrix_Type

REAL, ALLOCATABLE, SAVE         :: R_original(:,:,:)
TYPE (R_Matrix_Type), ALLOCATABLE, TARGET, SAVE :: R_reference(:)

! Definition of R-Matrix Types
INTEGER, PARAMETER ::  R_Full_Matrix   = 1
INTEGER, PARAMETER ::  R_Band_Diagonal = 2 
INTEGER, PARAMETER ::  R_Eigenvectors  = 3 


! Matrices resulting from Retrieval Inversion
!---------------------------------------------

REAL, ALLOCATABLE, SAVE :: Analysis_Error_Covariance(:,:)
REAL, ALLOCATABLE, SAVE :: Prop_Measurement_Noise(:,:)

! Satellite Identity Default Value (used in R-matrix)
!------------------------------------------------------

CHARACTER (LEN=20):: Default_SatID_Text = 'Generic Instrument'


END MODULE NWPSAFMod_CovarianceMatrices
