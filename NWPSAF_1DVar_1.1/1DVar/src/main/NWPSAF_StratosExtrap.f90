!+ Extrapolate the background profile beyond the model top.

SUBROUTINE NWPSAF_StratosExtrap( &
     Pressures, & ! in
     Tprofile)    !inout

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
! If needed extrapolate the temperature profile into the stratosphere.  
! In this case assume a simple isothermal extrapolation (not very good!)
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     15/01/11 Original code Andrew D. Collard.
! 3.0.4   04/03/02 Change to allow variable pressures. A.D. Collard
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Constants, ONLY : &
    MissData_R, &
    Tolerance,  &
    Upper_Transition

USE NWPSAFMod_RTmodel, ONLY : &
    Num_RTlevels

IMPLICIT NONE

! Subroutine arguments:
REAL,    INTENT(IN)    :: Pressures(:)  ! in Pa
REAL,    INTENT(INOUT) :: Tprofile(:)

! Local variables:
INTEGER :: level       ! Loop variable
REAL    :: Pressure    ! Temporary storage for pressure level
REAL    :: temp        ! temperature at one level

!- End of header --------------------------------------------------------------

DO level = Num_RTlevels,1,-1
   Pressure = Pressures(level) / 100.   ! convert to hPa
   IF ( Pressure > Upper_Transition ) THEN
      temp = Tprofile(level)
   ELSE IF (ABS(Tprofile(level) - MissData_R) <= Tolerance) THEN
      Tprofile(level) = temp
   END IF
END DO


END SUBROUTINE NWPSAF_StratosExtrap
