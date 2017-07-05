!+ Applies the "HIRS" cloud detection test.

Subroutine NWPSAF_CloudyOrNot( &
          Obs,              & ! inout
          MeasurementObs,   & ! in
          BMatrix,          & ! in
          R_Matrix,         & ! in
          UsedChans,        & ! in
          RT_Params,        & ! inout
          Valid_Data,       & ! inout
          Cloudy,           & ! inout
          HighCloud,        & ! inout
          RTerrorcode     ) ! out

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
!   This routine is designed for the detection of cloudy spots present in
!   IR sounding data, and uses a method similar that in the variational
!   assimilation of radiances for NWP. The result of calling this routine
!   is a flag (clear or cloudy) that is then passed on to the retrieval step
!   in the main routine.
!
!   There is also a HighCloud flag that checks whether there is high, thick 
!   cloud and hence the situation is hopeless! 
!
! Method:
!
! 1) Calculate cloud cost function.
! 2) A decision is made as to whether the spot is significantly contaminated
!    with cloud, by comparing the cost with a predetermined threshold value
!    which effectively marks the boundary between a declared clear/cloudy
!    sounding.
!
! Reference :
!   English, S.J., J.R. Eyre and J.A. Smith (1999). A cloud-detection scheme 
!     for use with satellite sounding radiances in the context of data 
!     assimilation for numerical weather prediction. 
!     Q.J.R. Meteorol. Soc., 125, 2359--2378. 
! 
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     08/06/00 Original code developed from ATOVS_CloudyOrNot. A. Collard
! 1.1     18/01/01 Fixed bug where U matrix wasn't being inverted! A. Collard.
! 1.2     18/01/01 Tidied up.  A. Collard.   Met Office.
! 2.2     08/05/02 Check for absence of window channel.   A. Collard.
! 3.0.5   31/03/04 Fix to stop absence 
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------
! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    ChannelSelection_Type, &
    Window_Channel

USE NWPSAFMod_Constants, ONLY : &
    MissData_R

USE NWPSAFMod_CovarianceMatrices, ONLY : &
    R_Matrix_Type

USE NWPSAFMod_ObsInfo, ONLY : &
    Ob_Type

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type

USE NWPSAFMod_Params, ONLY : &
    GeneralMode, &
    VerboseMode, &
    CalcRadiance, &
    RTsea, &
    CostThresh_land, &
    CostThresh_sea, &
    CostThresh_IRWindow_land, &
    CostThresh_IRWindow_sea, &
    CloudAbsThresh_IRWindow, &
    HighCloudAbsThresh_IRWindow  

IMPLICIT NONE

INCLUDE 'NWPSAF_CloudCost.interface'

! Subroutine arguments:

TYPE(Ob_type), INTENT(INOUT) :: Obs           ! Observed/Retrieval data 
REAL,    INTENT(IN)          :: MeasurementObs(:) ! Observed BTs
REAL, INTENT(IN)             :: BMatrix(:,:)  ! Background Error Covariance
TYPE(R_Matrix_Type), INTENT(IN)         :: R_Matrix   ! R matrix
TYPE(ChannelSelection_Type), INTENT(IN) :: UsedChans  ! Chans to be used
TYPE(RTParams_Type), INTENT(INOUT)      :: RT_Params  ! Info for RT Model
LOGICAL, INTENT(INOUT) :: Valid_Data                  ! Data flag
LOGICAL, INTENT(INOUT) :: Cloudy                      ! Data flag
LOGICAL, INTENT(INOUT) :: HighCloud                   ! Data flag
INTEGER, INTENT(OUT)   :: RTerrorcode                 ! Error from RT Model 

! Local variables:
REAL :: DeltaObs
REAL :: Cost                    ! Value from the cloud cost routine

!-----------------------------------------------------------------------------

!
!1) Calculate the cloudy cost function
!--
  
CALL NWPSAF_CloudCost( &
     MeasurementObs,       & ! in
     Bmatrix,          & ! in
     R_Matrix,         & ! in
     UsedChans,        & ! in
     RT_Params,        & ! inout
     Cost,             & ! out
     RTErrorCode)        ! out

! Scale cost by number of channels   
Cost = Cost / UsedChans % NumChans

CalcCost: IF ( RTerrorcode == 0 ) THEN
   Obs% CldCost = Cost
   
!2) Flag as cloudy if threshold breached
!--
! The IRWindow threshold applies to both land and sea, but we may want to
! upgrade this with separate values for each.

   IF (Window_Channel > 0) THEN
      IF ( CalcRadiance ) THEN
        DeltaObs = MeasurementObs(UsedChans % Channels(Window_Channel)) - &
           RT_Params % TotalRadiances(Window_Channel)
      ELSE
        DeltaObs = MeasurementObs(UsedChans % Channels(Window_Channel)) - &
          RT_Params % TotalBTs(Window_Channel)
      END IF
   ELSE
      DeltaObs = -1.e10  ! Set to large -ve value so as not to
                             ! result in an automatic clear case.
   END IF

   IF ( DeltaObs < CloudAbsThresh_IRWindow .AND. Window_Channel > 0) THEN
      Cloudy = .TRUE.
      IF ( DeltaObs < HighCloudAbsThresh_IRWindow ) &
           HighCloud = .TRUE.
   ELSE
      Setflag: IF ( RT_Params % RTSurfaceType /= RTsea ) THEN
         IF ( Cost >= CostThresh_land .AND. &
              DeltaObs < CostThresh_IRWindow_land ) THEN
            Cloudy = .TRUE.
         ELSE
            Cloudy = .FALSE.
         END IF
      ELSE
         IF ( Cost >= CostThresh_sea .AND. &
              DeltaObs < CostThresh_IRWindow_sea ) THEN
            Cloudy = .TRUE.
         ELSE
            Cloudy = .FALSE.
         END IF
      END IF Setflag
   END IF
   
ELSE
   Obs% CldCost = MissData_R
   Valid_data = .FALSE.
END IF CalcCost

IF ( GeneralMode >= VerboseMode ) THEN
   WRITE(*,*) 'Cost, CostThresh_Sea, CostThresh_Land:', &
        Cost, CostThresh_Sea, CostThresh_Land
   WRITE(*,*) 'CostThresh_IRWindow_land, CostThresh_IRWindow_Sea:', &
        CostThresh_IRWindow_land, CostThresh_IRWindow_Sea
   write(*,*) 'Number of Channels:',UsedChans % NumChans
   write(*,*) 'DeltaObs:',DeltaObs
   IF (Cloudy) THEN
      WRITE(*,*) 'Test Result is Cloudy'
   ELSE
      WRITE(*,*) 'Test Result is Clear'
   END IF
END IF

End Subroutine NWPSAF_CloudyOrNot
