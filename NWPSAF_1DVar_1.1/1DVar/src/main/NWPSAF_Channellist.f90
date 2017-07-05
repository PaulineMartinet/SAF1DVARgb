!+ Read in a channellist file and construct channel arrays

Subroutine NWPSAF_Channellist(fileunit)

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
! Description: Reads in channels to be used in various retrieval situations.
!              The use of each channel that is used at all is specified in 
!              the ChannelChoice.dat file through a binary code (described 
!              in the header to the file and in this subroutine).
!              Also reads in which channels are to be used in monitoring 
!              (stored in the BackChans structure) and in cloud detection
!              (stored in DetectCloudChans).
!
! Method:  Integer codes for each channel used in retrieval are stored in the 
!          ChannelChoice.dat file.  Each bit in the integer corresponds
!          to a sounding option as follows:
!
!          Bit Number        Sounding Option
!             1                  Sea
!             2                  SeaIce
!             3                  Land
!             4                  Highland
!             5                  Surface Mismatch
!             6                  Clear
!             7                  IRCloudy
!             8                  CloudFlag2
!             9                  CloudFlag3
!            10                  High Cloud
!
!   So if the code is 1023 = 1111111111 in binary, then the channel may be
!   used in all cases.  If the code is 33 = 100001 then the channel is only
!   used for clear skies above the sea.
! 
!   Extra codes may be added as required.  This may be especially useful if a 
!   more complicated situation were required (e.g., the channel can be used 
!   over land or with high cloud but not both) or for defining it's use in 
!   certain lat-long boxes.
!
!   A second column in the ChannelChoice.dat file stores information on 
!   the which channels are to be used in monitoring and cloud detection.
!   The method is similar to the one above with the first bit indicating 
!   whether a channel is used for monitoring and the second whether it
!   is used in cloud detection.  If this value is negative, it is also used
!   as a window channel during cloud detection.
!
! Bugs/Features:
!   
!   The functionality from ATOVS_Channellist of being able to read in the 
!   top and base retrieval levels is temporarily lost.  These are simply
!   set to default values here.  This can be reinstated later if deemed 
!   neccessary
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     15/10/99  Original version.  A. Collard
! 1.1     14/01/00  Added code for background and cloud detection
!                   channel selection here.  A. Collard.
! 1.2     07/02/00  File unit number now passewd from calling routine.
! 1.3     11/04/00  When observational errors are correlated, the channels
!                   must be sorted into order.  NWPSAF_IntegerSort has
!                   been added to do this.      A. Collard.
! 2.2     13/03/02  No longer replace the last channel with a window channel
!                   if the maximum number of channels has already been 
!                   assigned.                   
!                   Check if number of background and/or cloud detection
!                   channels are zero.  
!                   Issue warning if no window channel.         A. Collard.
! 2.3     07/06/02  Window channel specification is now via the 
!                   ChannelChoice.dat file (specified through a negative
!                   value in the cloud detection column) 
! 3.0.4   03/03/04  Removed MaxRetChans, MaxBackChans, MaxCloudChans. Now
!                   MaxChanUsed is the upper limit. Introduce BTEST.  
!                                                              A. Collard.
! Hereafter: Changes made under FCM
!
! Ticket  Date      Comment
! ------  --------  -------
! 19      30/03/10  Fix problems with swapping array order when compiling
!                   with GFortran. Fiona Hilton
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_Channellist, ONLY : &
    MaxSurfaces, &
    MaxClasses, &
    ChannelChoice, &
    Window_Channel, &
    BackChans, &
    DetectCloudChans

USE NWPSAFMod_Params, ONLY : &
     GeneralMode,   &
     VerboseMode,   &
     StatusFatal,   &
     StatusWarning, &
     Perform1DVar, &
     DetectCloud, &
     MaxChanUsed

IMPLICIT NONE

INCLUDE 'NWPSAF_Report.interface'
INCLUDE 'NWPSAF_IntegerSort.interface'

! Subroutine Arguments 

INTEGER, INTENT(IN) :: FileUnit

! Local constants:

CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_Channellist"

! Local Variables:

CHARACTER(LEN=80) :: ErrorMessage(5)  ! Message for NWPSAF_Report

INTEGER ::  I, IChannel, II, J
INTEGER ::  BG_and_Cloud_Mask_In
INTEGER ::  Channel_Number_In
INTEGER ::  Channel_Mask_In
INTEGER, ALLOCATABLE :: BG_and_Cloud_Mask(:)
INTEGER, ALLOCATABLE :: Channel_Number(:)
INTEGER, ALLOCATABLE :: Channel_Mask(:)
INTEGER, POINTER :: Index(:)  ! Used in ordering channels
INTEGER ::  Max_Channel_Chosen
INTEGER ::  Num_Channels
INTEGER ::  Num_Channels_to_Use
INTEGER ::  temparray(MaxChanUsed)
INTEGER :: ReadStatus            ! IOSTAT keyword value
LOGICAL, ALLOCATABLE::  Use_Channel(:)


INTEGER, ALLOCATABLE :: BackChans_Tmp(:)
INTEGER, ALLOCATABLE :: DetectCloudChans_Tmp(:)

!------------------------------------------------------------------------------

ErrorMessage(:)=' '

!1) Initialise variables
!-----------------------
!

Window_Channel = 0

! Set up channel Choice array 

DO i = 1, MaxSurfaces
  DO j = 1, MaxClasses
    ChannelChoice(i,j) % numchans = 0
    NULLIFY(ChannelChoice(i,j) % channels)
  ENDDO
ENDDO

!
! Assume minimal data in the header for now, i.e., the 
! the number of channels entered.
! Later versions will include a lot more!
! 

! Find out the number of channels in the channel selection file 
!--------
READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) Num_Channels
IF (ReadStatus /= 0 .OR. Num_Channels < 0 .OR. &
    Num_Channels > MaxChanUsed) THEN
  Errormessage(1) = ' Error Reading Channel Selection File'
  WRITE(Errormessage(2),'(a,i5,a,i5)') ' Num_Channels = ',Num_Channels, &
    ' Max Chans = ',MaxChanUsed
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=StatusFatal )
ENDIF

! Allocate various arrays
!-------
ALLOCATE (BG_and_Cloud_Mask(Num_Channels))
ALLOCATE (Channel_Number(Num_Channels))
ALLOCATE (Channel_Mask(Num_Channels))
ALLOCATE (Use_Channel(Num_Channels))
ALLOCATE (BackChans_Tmp(Num_Channels))
ALLOCATE (DetectCloudChans_Tmp(Num_Channels))
BackChans % numchans = 0
DetectCloudChans % numchans = 0

I=0

! Read in channel selection information.
!-------
DO WHILE (I < Num_Channels)
  READ (UNIT=fileunit, FMT=*, IOSTAT=ReadStatus) &
    Channel_Number_In, Channel_Mask_In, BG_and_Cloud_Mask_In
  IF (ReadStatus == -1) THEN
    Num_Channels = I
    Errormessage(1) = ' End of Channel Selection File Read'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusWarning )
  ELSE IF (ReadStatus /= 0) THEN
    Errormessage(1) = ' Error Reading Channel Selection File'
    CALL NWPSAF_Report(       &
      RoutineName,            &
      ErrorMessage,           &
      ErrorStatus=StatusFatal )
  ENDIF

! If the channel will never be used, go on to the next channel.  Otherwise
! assign it a channel number.

   IF (Channel_Number_In <= 0 .OR. &
       Channel_Number_In > MaxChanUsed .OR. &
       (Channel_Mask_In <= 0 .AND. BG_and_Cloud_Mask_In <= 0)) THEN
     Num_Channels = Num_Channels - 1
     CYCLE
   ELSE
     I = I+1
   ENDIF
   Channel_Number(I)    = Channel_Number_In
   Channel_Mask(I)      = Channel_Mask_In
   BG_and_Cloud_Mask(I) = BG_and_Cloud_Mask_In

ENDDO 

CLOSE (fileunit)


! Compare each surface type and class with the channel specifiers

Max_Channel_Chosen = 0

SurfaceType_Loop : DO i = 1, MaxSurfaces
  Class_loop :  DO j = 1, MaxClasses
    Num_Channels_to_Use = 0
    DO ichannel = 1, Num_Channels
      Use_Channel(ichannel) = &
        BTEST(Channel_Mask(ichannel),i-1) .AND. BTEST(Channel_Mask(ichannel),j+4)
      IF (Use_Channel(ichannel)) Num_Channels_to_Use = Num_Channels_to_Use + 1
      IF (ichannel > Max_Channel_Chosen) Max_Channel_Chosen = ichannel
      IF (Num_Channels_to_Use > MaxChanUsed) THEN
        Errormessage(1) = ' Too many retrieval channels requested'
        CALL NWPSAF_Report(       &
          RoutineName,            &
          ErrorMessage,           &
          ErrorStatus=StatusFatal )
      END IF
    ENDDO

! The channels that may be used in this situation are now flagged in 
! Use_Channel.
! 
! Now loop through Use_Channel and pick out the channels to be used
!
    ChannelChoice(i,j) % numchans = Num_Channels_to_Use
    ALLOCATE (ChannelChoice(i,j) % Channels(Num_Channels_to_Use))
    II = 1
    DO ichannel = 1, Max_Channel_Chosen
      IF (Use_Channel(ichannel)) THEN
        ChannelChoice(i,j) % Channels(II) = Channel_Number(ichannel)
        II = II + 1
      ENDIF
    ENDDO
! Sort the channels.  This is required for those cases where there are
! interchannel correlated errors
    ALLOCATE(Index(Num_Channels_to_Use))
    CALL NWPSAF_IntegerSort( &
      ChannelChoice(i,j) % Channels(:), & ! in
      Index) ! inout
    temparray(1:Num_Channels_to_Use)=ChannelChoice(i,j) % Channels(Index)
    ChannelChoice(i,j) % Channels(:) = temparray(1:Num_Channels_to_Use)
    DEALLOCATE(Index)
  ENDDO Class_Loop
ENDDO SurfaceType_Loop


! Check whether the channel will be used for Background Monitoring
! or Cloud Detection.  
!-------
DO ichannel = 1, Num_Channels
  IF (ABS(BG_and_Cloud_Mask(ichannel)) == 1 .OR.  &
      ABS(BG_and_Cloud_Mask(ichannel)) == 3 .AND. &
      BackChans % numchans < MaxChanUsed          ) THEN
    BackChans % numchans = BackChans % numchans + 1
    BackChans_Tmp(BackChans % numchans) = Channel_Number(ichannel)

  ! If retrievals are to be done, make sure that all of the chosen 
  ! retrieval channels are in the array of Background channels too.
  !-------
  ELSE IF (Perform1DVar .AND. ichannel <= Max_Channel_Chosen .AND. &
           Channel_Mask(ichannel) > 0) THEN
    BackChans % numchans = BackChans % numchans + 1
    BackChans_Tmp(BackChans % numchans) = Channel_Number(ichannel)
  ENDIF
  IF (ABS(BG_and_Cloud_Mask(ichannel)) == 2 .OR.  &
      ABS(BG_and_Cloud_Mask(ichannel)) == 3 .AND. &
      DetectCloudChans % numchans < MaxChanUsed   ) THEN
    DetectCloudChans % numchans = DetectCloudChans % numchans + 1
    DetectCloudChans_Tmp(DetectCloudChans % numchans) = Channel_Number(ichannel)
    IF (Window_Channel == 0 .AND. BG_and_Cloud_Mask(ichannel) < 0) &
      Window_Channel = DetectCloudChans % numchans
  ENDIF
ENDDO

!  Report errors in setting up channels.
!----------
IF (Window_Channel == 0) THEN
  Errormessage(1) = ' No Window Channel Set'
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=StatusWarning )
ELSE IF (GeneralMode >= VerboseMode) THEN
   WRITE(*,*) 'Window Channel is ',DetectCloudChans_Tmp(Window_Channel)
END IF

IF (BackChans % numchans == 0) THEN
  Errormessage(1) = ' No Background Channels have been set up'
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=StatusFatal )
ENDIF

IF (DetectCloud .AND. DetectCloudChans % numchans == 0) THEN
   Errormessage(1) = ' No Cloud Detection Channels have been set up'
  CALL NWPSAF_Report(       &
    RoutineName,            &
    ErrorMessage,           &
    ErrorStatus=StatusFatal )
ENDIF

ALLOCATE(BackChans % Channels(BackChans % numchans))
ALLOCATE(DetectCloudChans % Channels(DetectCloudChans % numchans))

! Sort the channels in the BackChans and DetectCloudChans arrays.  
! This is required for those cases where there are
! interchannel correlated errors
!-------

ALLOCATE(Index(BackChans % numchans))
CALL NWPSAF_IntegerSort( &
  BackChans_Tmp(1:BackChans % numchans), & ! in
  Index) ! inout
BackChans % Channels(:) = BackChans_Tmp(Index)
DEALLOCATE(Index)
ALLOCATE(Index(DetectCloudChans % numchans))
CALL NWPSAF_IntegerSort( &
  DetectCloudChans_Tmp(1:DetectCloudChans % numchans), & ! in
  Index) ! inout
DetectCloudChans % Channels(:) = DetectCloudChans_Tmp(Index)
IF (Window_Channel /= 0) Window_Channel = Index(Window_Channel)
DEALLOCATE(Index)

DEALLOCATE (BG_and_Cloud_Mask)
DEALLOCATE (Channel_Number)
DEALLOCATE (Channel_Mask)
DEALLOCATE (Use_Channel)
DEALLOCATE(BackChans_Tmp)
DEALLOCATE(DetectCloudChans_Tmp)
NULLIFY (Index)

End Subroutine NWPSAF_Channellist
