! Uses minimum residual method to calculate cloud top height
! and cloud fraction

Subroutine NWPSAF_CO2Slice(   &
           BTObserved,      & !in
           RT_Params,       & !inout
           R_Matrix,        & !in
           UsedChans,       & !in
           Total_Channels)    !in

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
!   This routine is designed for the evaluation of cloud top
!   height and cloud fraction.  This requires RTTOV to be the
!   RT model.  Also UsedChans should be the same as that used 
!   for the previous call to the RT model.
!
!
! History:
!
!  Version Date       Comment
! ----------------------------------- 
!  1.0     4/04/01    Original code Matthew D E Szyndel.
!  2.0    16/04/02    Rewrite for version 2.0 of 1DVAR code
!                                       Matthew D E Szyndel.
!  3.0.1  21/07/03    Added into new NWPSAF_1DVar with some small changes.  
!                                                             Andrew Collard.
!  3.0.3  01/03/04    Move channel constants into RT_Params.  Andrew Collard.
!  3.0.4  05/03/04    Replace RTIASI_Refernce_Pressures (obsolete) with
!                     RT_Params%Pressure_Pa.                  Andrew Collard.
!  3.0.6. 18/06/04    Set RTParams % RT1stGuess.              Andrew Collard.
!  3.0.7. 06/12/05    Corrected indexing of R_matrix and
!                     fixed factor of 2 in Press_Inc          Ed Pavelin
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: GTDP 8
!
! Declarations:
! Modules used:
! Imported Array Variables with intent (in):

USE NWPSAFMod_RTmodel, ONLY:&
  RTParams_Type,          &
  Num_RTLevels,           &
  Prof_CloudCover,        &
  Prof_CTP

USE NWPSAFMod_Channellist, ONLY : &
  DetectCloudChans,    &
  ChannelSelection_Type

USE NWPSAFMod_CovarianceMatrices, ONLY : &
  R_Matrix_Type

USE NWPSAFMod_Params, ONLY : &
  GeneralMode,       &
  VerboseMode,       &
  StatusFatal

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'

! Subroutine Arguments

REAL, INTENT(IN)                  :: BTObserved(:) ! Observed BTs
TYPE(RTParams_Type), INTENT(INOUT):: RT_Params     ! Info for RT Model
TYPE(R_Matrix_Type), INTENT(IN)   :: R_Matrix      ! R matrix
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans  ! Chans to be used
INTEGER, INTENT(IN)               :: Total_Channels 

! Local Constants

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NWPSAF_CO2Slice'

! Local Variables

INTEGER :: IChan                             ! Channel index (sequential)
INTEGER :: JChan                             ! Channel index (sequential)
INTEGER :: Channel                           ! Channel index (non-sequential)
INTEGER :: IPress_Lev                        ! Pressure level index
INTEGER :: MinLevel(1)                       ! Pressure level of cloud top
INTEGER :: TropLevel                         ! Highest possible level for cloud

REAL    :: zteff                             ! Effective BT
REAL    :: zteff2                            ! bcon2/Effective BT
REAL    :: sigma                             ! Standard deviation of modelled 
                                             ! radiance
REAL    :: sigma2(Total_Channels)            ! sigma^2
REAL    :: ObservedRadiance(Total_channels)  ! Observed Radiances
REAL    :: MixedRadiance(Total_Channels,Num_RTLevels-1)
                                             ! Superpostion radiances
REAL    :: Cld_Frac_Vec_Num(Num_RTLevels-1)  ! Intermediate values in 
REAL    :: Cld_Frac_Vec_Den(Num_RTLevels-1)  ! calculation of cloud
REAL    :: Cld_Frac_Vec_Num_Inc              ! fraction
REAL    :: Cld_Frac_Vec_Den_Inc              !
REAL    :: Cld_Frac_Vec(Num_RTLevels-1)      ! Cloud fraction vector
REAL    :: Cost_Func(Num_RTLevels-1)         ! Cost function
REAL    :: Cost_Func_Inc                     ! Used in calculation
REAL    :: Cost_Func_Inc_Root                ! of cost function
REAL    :: Pressures(Num_RTLevels-1)         ! Pressure levels in hPa
REAL    :: Press_Inc_Num                     ! used to interpolate
REAL    :: Press_Inc_Den                     ! cloud top from
REAL    :: Press_Inc                         ! pressure levels
REAL    :: Max_Press                         ! CTP above which cloud free 
                                             ! conditions are assumed
CHARACTER(LEN=80)   :: ErrorMessage(2)   ! Message for NWPSAF_Report

! End of header ------------------------------------------------

ErrorMessage(:) = ''
TropLevel       = 22
Cld_Frac_Vec(:) = 0.00
Max_Press = 950.0
Pressures(:) = RT_Params % Pressure_Pa(1:Num_RTLevels-1) / 100.  ! Pa => hPa

! ------ CO2 Slicing Routine                     ---------------
! ------This converts observed BTs into radiances---------------

DO IChan=1,UsedChans % NumChans

  Channel = UsedChans % Channels(IChan)

  zteff = RT_Params % tc2(Ichan) * BTObserved(Channel) + RT_Params % tc1(IChan)

  ObservedRadiance(IChan) = &
    RT_Params % bcon1(IChan) / (exp(RT_Params % bcon2(IChan)/zteff) - 1.0)

  zteff2 = &
    LOG(1.0 + (RT_Params % bcon1(IChan) / RT_Params % TotalRadiances(IChan)))

  sigma = RT_Params % TotalRadiances(IChan) * &
         ( RT_Params % bcon1(IChan) + RT_Params % TotalRadiances(IChan) ) * &
          RT_Params % tc2(IChan) / &
         ( RT_Params % bcon1(IChan) * RT_Params % bcon2(IChan) ) * zteff2*zteff2

  sigma2(IChan)=sigma*sigma*R_Matrix % Matrix(Channel,0)
END DO

! ------This finds the optimum cloud fraction at each level-------------
! ------A specific channel choice list for the minimum residual---------
! ------retrieval has now been implemented.                 ------------

DO IPress_Lev = TropLevel,Num_RTLevels-1
  Cld_Frac_Vec_Num(IPress_Lev)= 0.0
  Cld_Frac_Vec_Den(IPress_Lev)= 0.0

  JCHAN=1
  DO IChan = 1,UsedChans % NumChans
    IF (UsedChans % Channels(IChan) < DetectCloudChans % Channels(JChan)) THEN
      CYCLE
    ELSE IF (UsedChans % Channels(IChan) == &
         DetectCloudChans % Channels(JChan)) THEN
      Cld_Frac_Vec_Num_Inc =                                          & 
                ( RT_Params % CloudyRadiances(IChan,IPress_Lev) -     &
                  RT_Params % TotalRadiances(IChan)               ) * &
                ( ObservedRadiance(IChan) -                           & 
                  RT_Params % TotalRadiances(IChan)               ) / &
                sigma2(IChan)

      Cld_Frac_Vec_Den_Inc =                                          &
                ( RT_Params % CloudyRadiances(IChan,IPress_Lev) -     &
                  RT_Params % TotalRadiances(IChan)               ) * &
                ( RT_Params % CloudyRadiances(IChan,IPress_Lev) -     &
                  RT_Params % TotalRadiances(IChan)               ) / &
                sigma2(IChan)

      Cld_Frac_Vec_Num(IPress_Lev) = Cld_Frac_Vec_Num(IPress_Lev) + &
                                     Cld_Frac_Vec_Num_Inc

      Cld_Frac_Vec_Den(IPress_Lev) = Cld_Frac_Vec_Den(IPress_Lev) + &
                                     Cld_Frac_Vec_Den_Inc

      JChan = JChan + 1
    ELSE
      ErrorMessage(1) = 'Missing Minimum Residual Channel.'
      ErrorMessage(2) = &
         'Somewhat surprising as BackChans is set to include DetectCloudChans!'
      CALL NWPSAF_Report( &
        RoutineName,            & !in
        ErrorMessage,           & !in
        ErrorStatus=StatusFatal ) !in
    END IF
    IF (JChan > DetectCloudChans % NumChans) EXIT
  END DO
  Cld_Frac_Vec(IPress_Lev) = Cld_Frac_Vec_Num(IPress_Lev)&
                           / Cld_Frac_Vec_Den(IPress_Lev)
  Cld_Frac_Vec(IPress_Lev) = MIN(Cld_Frac_Vec(IPress_Lev),1.00)
  Cld_Frac_Vec(IPress_Lev) = MAX(Cld_Frac_Vec(IPress_Lev),0.001)
END DO

! ------This finds the cost function at each level--------------

DO IPress_Lev = TropLevel,Num_RTLevels-1
  Cost_Func(IPress_Lev) = 0.0
  JChan=1
  DO IChan = 1, UsedChans % NumChans
    IF (UsedChans % Channels(IChan) < DetectCloudChans % Channels(JChan)) THEN
      CYCLE
    ELSE IF (UsedChans % Channels(IChan) == &
             DetectCloudChans % Channels(JChan)) THEN
      MixedRadiance(IChan,IPress_Lev)= &
        RT_Params % TotalRadiances(IChan) * (1.0 - Cld_Frac_Vec(IPress_Lev)) + &
        RT_Params % CloudyRadiances(IChan,IPress_Lev) * Cld_Frac_Vec(IPress_Lev)
      Cost_Func_Inc_Root = &
        ObservedRadiance(IChan) - MixedRadiance(IChan,IPress_Lev)
      Cost_Func_Inc = Cost_Func_Inc_Root * Cost_Func_Inc_Root / sigma2(IChan)
      Cost_Func(IPress_Lev) = Cost_Func(IPress_Lev) + Cost_Func_Inc
      JChan = JChan + 1
    ELSE
      ErrorMessage(1) = 'Missing Cloud Detection Channel.'
      ErrorMessage(2) = 'Very surprising as this has already been checked!'
      CALL NWPSAF_Report( &
        RoutineName,           & !in
        ErrorMessage,          & !in
        ErrorStatus=StatusFatal) !in
    END IF
    IF (JChan > DetectCloudChans % NumChans) EXIT
  END DO
END DO

! ------This finds the minimum cost function pressure --------------
!       and interpolates cloud top pressure off of RTTOV
!       pressure levels with dP = -J'(P)/J''(P) from 2nd
! ------order Taylor expansion about minimum level    --------------

MinLevel = MINLOC(Cost_Func(TropLevel:Num_RTLevels-1)) + TropLevel -1

IF (MinLevel(1) >= (TropLevel + 1) .AND. MinLevel(1) < Num_RTLevels-1) THEN

  Press_Inc_Num = (Cost_Func(MinLevel(1)+1)-Cost_Func(MinLevel(1)-1))&
                * (Pressures(MinLevel(1)+1)-Pressures(MinLevel(1))) &
                * (Pressures(MinLevel(1))-Pressures(MinLevel(1)-1))

  Press_Inc_Den = Cost_Func(MinLevel(1)+1) * &
                (Pressures(MinLevel(1))-Pressures(MinLevel(1)-1))&  
                - Cost_Func(MinLevel(1)) * &
                (Pressures(MinLevel(1)+1)-Pressures(MinLevel(1)-1))&
                + Cost_Func(MinLevel(1)-1) * &
                (Pressures(MinLevel(1)+1)-Pressures(MinLevel(1)))

  Press_Inc = Press_Inc_Num / (2.0*Press_Inc_Den)
ELSE IF (MinLevel(1) == TropLevel .OR. MinLevel(1) == Num_RTLevels-1) THEN
  Press_Inc = 0.0 
ELSE
  ErrorMessage(1) = &
       'Error in CO2 Slicing routine: Best pressure outside troposphere!'
  CALL NWPSAF_Report( &
    RoutineName,           & !in
    ErrorMessage,          & !in
    ErrorStatus=StatusFatal)!in
END IF

!------------------- Now we set the first guess to the CO2 slicing profile -----

RT_Params % RTguess(Prof_CTP) = &
     (Pressures(MinLevel(1)) + Press_Inc)
IF (RT_Params % RTguess(Prof_CTP) < Max_Press) THEN
   RT_Params % RTguess(Prof_CloudCover) = Cld_Frac_Vec(MinLevel(1))
ELSE
   RT_Params % RTguess(Prof_CloudCover) = 0.001
END IF

! Set the background to median values
RT_Params % RTback(Prof_CTP) =        500.
RT_Params % RTback(Prof_CloudCover) = 0.5 

! And set the 1st guess.
RT_Params % RT1stGuess(Prof_CTP) = RT_Params % RTguess(Prof_CTP)     
RT_Params % RT1stGuess(Prof_CloudCover) = RT_Params % RTguess(Prof_CloudCover)


IF (GeneralMode >= VerboseMode) THEN
  WRITE(*,*) 'First Guess Cloud Cover, CTP = ', &
       RT_Params % RTguess(Prof_CloudCover), RT_Params % RTguess(Prof_CTP)
END IF

!------------------- And finally estimate the first guess BTs ------------------
!------------------- from the last call of RTTOV              ------------------

DO IChan = 1,UsedChans % NumChans
  Channel = UsedChans % Channels(IChan)
  MixedRadiance(IChan,MinLevel(1))= RT_Params % TotalRadiances(IChan)&
                  * (1.0 - Cld_Frac_Vec(MinLevel(1)))&
                  + RT_Params % CloudyRadiances(IChan,MinLevel(1))&
                  * Cld_Frac_Vec(MinLevel(1))
  RT_Params % TotalRadiances(IChan) = MixedRadiance(IChan,MinLevel(1))
  zteff = RT_Params % bcon2(IChan) / ALOG( 1.0 + RT_Params % bcon1(IChan) / &
          RT_Params % TotalRadiances(IChan))
  RT_Params % TotalBTs(IChan) = &
    (zteff - RT_Params % tc1(IChan)) / RT_Params % tc2(IChan)
END DO

END SUBROUTINE NWPSAF_CO2Slice
