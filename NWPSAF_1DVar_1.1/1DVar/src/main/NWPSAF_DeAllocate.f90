SUBROUTINE NWPSAF_DeAllocate ( &
     Observations, &  ! inout
     BackGrModelOb, & ! inout
     RT_Params)       ! inout

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
!------------------------------------------------------------------------------
! Description: DEALLOCATE arrays allocated elsewhere in the NWPSAF_1DVar code.
!
! Method:
!  
! History:
!
! Version  Date      Comment
! -------  ----      -------
! 1.0      16/01/01  Original Version
! 2.2      03/05/02  Add deallocation of RT_Params % Absolute_Channel_Number
!                    Remove deallocation of R_Reference % Channel_Number
! 2.3      22/05/02  Remove references to SatID_Convert.  A. Collard.
! 3.0.4    03/03/04  Changes for generalisation of RT models. A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!

! Modules used:

USE NWPSAFMod_ObsInfo, ONLY : &
     Ob_Type, &
     ModelOb_Type

USE NWPSAFMod_Channellist, ONLY: &
     BackChans, &
     DetectCloudChans, &
     ChannelChoice, &
     MaxSurfaces, &
     MaxClasses

USE NWPSAFMod_Params, ONLY: &
     Retrieved_Elements

USE NWPSAFMod_RTModel, ONLY : &
     RTParams_Type

USE NWPSAFMod_CovarianceMatrices, ONLY : &
     R_Eigenvectors, &
     R_Reference

IMPLICIT NONE 

! Subroutine Arguments

TYPE(Ob_type),      INTENT(INOUT) :: Observations  ! Observed/Retrieval data 
TYPE(ModelOb_type) ,INTENT(INOUT) :: BackGrModelOb ! Background data
TYPE(RTParams_type),INTENT(INOUT) :: RT_Params     ! RT model data

! Local variables

INTEGER I, J

!-------------------------------------------------------------------------
! Start DeAllocating
!-------------------------------------------------------------------------

! DEALLOCATE SoundingType dependent arrays

DEALLOCATE(RT_Params % SeriesChoice) 
DEALLOCATE(RT_Params % PlatformChoice)
DEALLOCATE(RT_Params % SubTypeChoice)
DEALLOCATE(RT_Params % Instrument_Number)
DEALLOCATE(RT_Params % Absolute_Channel_Number)
DEALLOCATE(RT_Params % First_Channel_for_Instrument)
DEALLOCATE(RT_Params % RTEmissivity)
!PM
DEALLOCATE(RT_Params % view_angle)
!PM

DEALLOCATE(RT_Params % SatID)
NULLIFY(RT_Params % Pressure_Pa)

! DEALLOCATE other SoundingType dependent arrays
DEALLOCATE(Retrieved_Elements)

! DEALLOCATE Observation Elements

IF (ASSOCIATED(Observations % BriTemp) ) DEALLOCATE(Observations % BriTemp)
IF (ASSOCIATED(Observations % Radiance) ) DEALLOCATE(Observations % Radiance)
IF (ASSOCIATED(Observations % PCScore) ) DEALLOCATE(Observations % PCScore)

! DEALLOCATE BackGrModelOb Elements

DEALLOCATE( BackGrModelOb % Pstar)
DEALLOCATE( BackGrModelOb % rh2)
DEALLOCATE( BackGrModelOb % t2)
DEALLOCATE( BackGrModelOb % Tskin)
DEALLOCATE( BackGrModelOb % u10)
DEALLOCATE( BackGrModelOb % v10)

DEALLOCATE( BackGrModelOb % rh)
DEALLOCATE( BackGrModelOb % t)
DEALLOCATE( BackGrModelOb % p)

IF (ASSOCIATED(BackGrModelOb % BriTemp) ) DEALLOCATE( BackGrModelOb % BriTemp)

!  Free up background and cloud detection channel pointers

IF (ASSOCIATED(BackChans % Channels)) DEALLOCATE(BackChans % Channels)
IF (ASSOCIATED(DetectCloudChans % Channels)) &
     DEALLOCATE(DetectCloudChans % Channels)

! Free up retrieval channel pointers

RetrievalChannelLoop: DO i = 1, MaxSurfaces
   DO j = 1, MaxClasses
      IF ( ChannelChoice(i,j) % numchans > 0 ) THEN
         DEALLOCATE( ChannelChoice(i,j) % channels )
         NULLIFY( ChannelChoice(i,j) % channels )
      END IF
   END DO
END DO RetrievalChannelLoop
   
! Deallocate R-Matrix structure elements
!----

RMatrixLoop: DO i = 1, RT_Params % Num_SatIDs
   IF (R_Reference(i) % Num_Chans /= 0) THEN
      DEALLOCATE(R_Reference(i) % Matrix)
      IF (R_Reference(i) % Rtype == R_Eigenvectors) &
           DEALLOCATE(R_Reference(i) % Eigenvalues)
   END IF
END DO RMatrixLoop
DEALLOCATE(R_Reference)

END SUBROUTINE NWPSAF_DeAllocate
