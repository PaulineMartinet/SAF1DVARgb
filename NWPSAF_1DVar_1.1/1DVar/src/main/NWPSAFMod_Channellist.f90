Module NWPSAFMod_Channellist

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
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     27/10/99 Original version.  Andrew Collard.   Met Office.
! 3.0.4   03/03/04 Add MaxChanUsed    Andrew Collard.
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------


TYPE ChannelSelection_Type
  INTEGER :: NumChans
  INTEGER, POINTER :: Channels(:)
END TYPE

INTEGER, PARAMETER :: MaxSurfaces = 5
INTEGER, PARAMETER :: MaxClasses = 5

TYPE(ChannelSelection_Type), SAVE :: BackChans
TYPE(ChannelSelection_Type), SAVE :: DetectCloudChans
TYPE(ChannelSelection_Type), SAVE :: ChannelChoice(MaxSurfaces,MaxClasses)

INTEGER, SAVE :: Window_Channel
INTEGER, SAVE :: MaxChanUsed = 10000  ! Max number of channels allowed. 

End Module NWPSAFMod_Channellist
