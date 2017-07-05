!+ Check that profile elements are within limits.

Subroutine NWPSAF_CheckIteration( &
  Guess_Prof,              & ! inout
  Profile_incs,            & ! inout
  Profile_Variables_Reset, & ! out
  Out_of_Range )             ! out

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
!   Reset the profile elements if they are outside of the limits that allow 
!   the forward model to produce reasonable results.  
!   Also if specific humidities are super-saturated, surface emissivity is 
!   outside 0.0:1.0, cloud top pressure is outside of 200 hPa (?) to MSLP, 
!   cloud amount is outside 0.0:1.0. then those elements in the profile are 
!   set to be within the limits.
!   In all the above cases  Profile_Variables_Reset is set to .TRUE.
!   
!   If the temperatures are outside of the gross limits allowed, 
!        Out_of_range is set to .TRUE. 
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     12/05/00 Original code developed from ATOVS_CheckIteration
!                  Andrew Collard.
! 1.1     7/12/00  Relax the requirements vs. soft limits.  
!                  (Add 0.5 tolerance in Ln(q) and allow supersaturation)
! 3.0.1   15/07/03 Add Cloud stuff (M. Szyndel).   Andrew Collard.
! 3.0.2   11/02/04 Relax cloud constraints.        Andrew Collard.
! 3.0.4   05/03/04 Replace PLevels_RTmodel_Pa with Pressure.  A. Collard.
! 3.0.6   16/06/04 Remove soft-limit tolerances to fastmodel interfaces
!                                                             A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket  Date     Comment
! ------- -------- -------
! 25      20/02/12 Don't check profile against softlimits if UseModelLevels 
!                  is true (i.e. RTTOV is interpolating the profile)  P.Weston.
!
! Modules used:

USE NWPSAFMod_Constants, ONLY : &
    MaxTemperature, &
    MinTemperature

USE NWPSAFMod_RTmodel, ONLY : &
    ProfSize, &
    Num_ProfElementsUsed, &
    Num_RTlevels, &
    Soft_Limits, &
    Prof_FirstT, &
    Prof_T2, &
    ! Prof_q2, &
    Prof_Tstar, &
    Prof_Pstar, &
    Prof_CTP, &
    Prof_CloudCover, &
    UseModelLevels

USE NWPSAFMod_Params, ONLY : &
    GeneralMode,        &
    VerboseMode,        &
    Retrieved_Elements, &
    Cloud_Min_Pressure, &
    Cloudfree,          &
    ThickCloud

IMPLICIT NONE

! Subroutine arguments:
REAL,    INTENT(INOUT) :: Guess_Prof(:)       ! Input (guess) profile
REAL,    INTENT(INOUT) :: Profile_incs(:)     ! This is the increment vector
                                              ! used in the retrieval
LOGICAL, INTENT(OUT)   :: Profile_Variables_Reset
LOGICAL, INTENT(OUT)   :: Out_of_range

! Local variables:
REAL    :: Profile_Data
INTEGER :: element, element1
INTEGER :: level
REAL :: Original_Prof(ProfSize)

!----------------------------------------------------------------------------

! ------------------
! 1. Initialize data
! ------------------

Profile_Variables_Reset = .FALSE.
Out_of_range = .FALSE.

!Save incoming profile
DO element1 = 1,Num_ProfElementsUsed
   element = Retrieved_Elements(element1)
   Original_prof(element) = Guess_Prof(element)
END DO

! ---------------------
! 2. Check Temperatures
! ---------------------

!Check temperature profile within limits, else fail
DO level = 1,Num_RTlevels
   IF ( (Guess_Prof(Prof_FirstT+level-1) > MaxTemperature) .OR. &
        (Guess_Prof(Prof_FirstT+level-1) < MinTemperature)      ) THEN
      Out_of_range = .TRUE.
      IF ( GeneralMode >= VerboseMode ) THEN
         WRITE(*,*) 'INVALID DATA: Temperature outside limits at level ',&
            level, Guess_Prof(Prof_FirstT+level-1)
      ELSE
         EXIT
      END IF
   END IF
END DO

!Check surface temps within limits, else fail
IF ( Guess_Prof(Prof_T2)  > MaxTemperature  .OR. &
     Guess_Prof(Prof_T2)  < MinTemperature  .OR. &
     Guess_Prof(Prof_Tstar) > MaxTemperature  .OR. &
     Guess_Prof(Prof_Tstar) < MinTemperature ) THEN
   Out_of_range = .TRUE.
   IF ( GeneralMode >= VerboseMode ) THEN
      WRITE(*,*) 'INVALID DATA: Temperature outside limits at the surface'
   END IF
END IF


! ----------------------------------------------------
! 3. Constrain profile to within range of validity for 
!    the fastmodel and perform further constraints on
!    humidities, emissivity, cloud variables
! ----------------------------------------------------

Constrain: IF ( .NOT. Out_of_range ) THEN
   
   ! Constrain profiles to within bounds allowed by the fastmodel
   ! Can only do this if not interpolating to RTTOV levels, but if you are
   ! interpolating it shouldn't be necessary as profiles % apply_reg_limits
   ! is set to true if using model levels inside NWPSAF_RTTOV10_Interface

   IF (.NOT. UseModelLevels) THEN   
     DO element1 = 1,Num_ProfElementsUsed
        element = Retrieved_Elements(element1)
        Guess_Prof(element) = &
          MIN(Guess_Prof(element),Soft_Limits % Maximum(element))
        Guess_Prof(element) = &
          MAX(Guess_Prof(element),Soft_Limits % Minimum(element))
     END DO
   END IF
      
   !Constrain cloud pressure and cloud amount (allow some tolerance).
   Profile_Data = MIN(Guess_Prof(Prof_CTP), Guess_Prof(Prof_Pstar)*1.1) ! >= P*
   Guess_Prof(Prof_CTP) = MAX(Profile_Data,Cloud_Min_Pressure*0.9) ! <= 200hPa?
   Profile_Data = MIN(Guess_Prof(Prof_CloudCover),ThickCloud * 1.1) ! <= 1
   Guess_Prof(Prof_CloudCover) = MAX(Profile_Data,Cloudfree - 0.1)  ! >= 0
   
   ! ---------------------------------------------------------------
   ! 4. Update profile increments and indicate whether any profile
   !    elements have been reset.
   ! ---------------------------------------------------------------
   
   !Change profile increments to match changes to profile
   DO element1 = 1,Num_ProfElementsUsed
      element = Retrieved_Elements(element1)
      IF (Original_prof(element) /= Guess_Prof(element)) THEN
         Profile_Variables_Reset = .TRUE.
        IF ( GeneralMode >= VerboseMode ) THEN
          WRITE(*,*) 'Check Iteration: Elements Reset: ', &
              element, Original_prof(element),Guess_Prof(element)
        END IF
         Profile_incs(element1) = Profile_incs(element1) + &
              Guess_Prof(element)     - &
              Original_Prof(element)
      END IF
   END DO
   
END IF Constrain


End Subroutine NWPSAF_CheckIteration
