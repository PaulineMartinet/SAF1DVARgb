Subroutine NWPSAF_BandMultiply ( &
     A_band,                & ! in
     Matrix,                & ! in
     Num_Elements_A,        & ! in
     Num_Bands,             & ! in
     Num_Elements_M,        & ! in
     Matrix_Out)              ! out

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
!
! Description:
!
!  Multiply the matrix, Matrix, by the band diagonal matrix, A_band.
!
! References: 
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     23/02/00 Original Version.  Andrew Collard.  Met Office.
! 1.1     02/03/00 Finally got it to give the right answer!  ADC.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments:

INTEGER, INTENT(IN) :: Num_Elements_A      
INTEGER, INTENT(IN) :: Num_Bands           ! Number of off-diagonals in A_band
INTEGER, INTENT(IN) :: Num_Elements_M      
REAL, INTENT(IN)  :: &
     A_band(0:Num_Bands,Num_Elements_A)             ! Band diagonal matrix
REAL, INTENT(IN)  :: &
     Matrix(Num_Elements_M, Num_Elements_A)         ! Other matrix
REAL, INTENT(OUT) :: &
     Matrix_Out(Num_Elements_M, Num_Elements_A)    ! = Matrix#A_Band

! Local constants:

! Local variables:

INTEGER :: I, J, K, Max_Bound

!-----------------------------------------------------------------------------

DO I=1,Num_Elements_M
   DO J=1,Num_Elements_A
      Max_Bound = MIN(Num_Bands, Num_Elements_A - J)
      Matrix_Out(I,J) = DOT_PRODUCT(Matrix(I,J:J+Max_Bound), &
           A_Band(0:Max_Bound,J))
      DO K=MAX(1,J-Num_Elements_A),MIN(J-1,Num_Bands)
         Matrix_Out(I,J) = Matrix_Out(I,J) + &
              Matrix(I,J-K)*A_band(K,J-K)
      END DO
   END DO
END DO

End Subroutine NWPSAF_BandMultiply
