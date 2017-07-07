!+ Physical and other constants for global access in NWPSAF_1DVar

Module NWPSAFMod_Constants

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
!---------------------------------------------------------------------------
! Description:
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.1     30/09/98 Original code developed from GLOSSMod_Constants.
! 1.1              Andrew J. Smith
! 1.8     12/11/98 Updated thresholds for cloud test. AJS
! 1.9     18/11/98 Threshold for scattering index increased from 10 to 25. RJR
! 1.10    19/11/98 Removed ozone profile. AJSmith
! 1.16    06/01/99 New AMSU scattering index threshold value. AJSmith
! 1.25    08/03/99 Added reference temperature profile. AJSmith
! NWPSAFMod_Constants:
! 1.0     12/06/00 Original version based on ATOVSMod_Constants.  A Collard 
! 1.1     12/06/00 Moved cloud cost thresholds to NWPSAFMod_Params.  A Collard 
! 1.2     21/11/00 Decreased min_q from 3.e-06 to 1.e-08 (was above largest
!                  value that RTIASI could cope with!).  A. Collard.
!                                                        Met Office.
! 1.3     01/03/01 Corrected orbital heights added NOAA-16 R. Saunders
! 1.4     01/03/01 Added Pi                               M.D.E. Szyndel
! 3.0.1   15/07/03 Add Aqua and Geostationary orbit heights.   A.D. Collard.
! 3.0.4   04/03/04 Remove surplus parameters. Add Gastropod limits. A. Collard.
!
! Hereafter: Changes made under FCM
!
! Ticket    Date     Comment
! ------    ----     -------
! 28        22/02/12 Added parameters for cloud liquid water 
!                    retrieval.   TR Sreerekha
!
! Code Description:
!   Language:       Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:
USE rttov_const, ONLY : &
    mair, &
    mh2o

IMPLICIT NONE
SAVE

!1) Missing data
!---------------
!Assigned to OPS/UM/VAR values


INTEGER, PARAMETER :: &
         MissData_I = -9999     ! integer missing data indicator

REAL,    PARAMETER :: &  
         MissData_R = -9999., &  ! real missing data indicator
         Tolerance  = 0.01       ! Margin of error when testing for missing data


!2) Miscellaneous constants
!--------------------------

REAL, PARAMETER :: &
  MaxSurfaceP         =   1200.0,    &  ! ( hPa )
  MinSurfaceP         =    300.0,    &  ! ( hPa )
  Min_q               =      1.0E-8, &  ! ( kg / kg )
  MaxTemperature      =    340.0,    &  ! ( K ) 
  !PM
  MinTemperature      =     3.0,    &  ! ( K )
  !PM
  MaxRadiance         =    150.0,    &  ! (RTTOV radiance units mW sr-1 m-2 cm)
  MinRadiance         =      0.0,    & 
  Upper_Transition    =     10.0,    &  ! Define transition region for
  Lower_Transition    =     30.0,    &  ! stratospheric extrapolation (hPa)
!  Pi                  =      3.14159265358979323846,&   ! Pi
  ScaleFactorP        =     1.0,   &    ! Additional cost function coefficient
                                        ! for Cloud Top Pressure
  ScaleFactorF        =     1.e12       ! Additional cost function coefficient
                                        ! for Cloud Cover
REAL, PARAMETER              :: RH_1=0.95 ! below RH_1 no Cloud liquid water
REAL, PARAMETER              :: RH_2=1.05 ! above RH_2, water vapor content remains fixed
REAL, PARAMETER              :: SplitQtotal=0.50 ! split factor for water vapor and cloud liquid water


REAL, PARAMETER :: InverseTol = TINY( 0.0 )

! 3. Conversion constants 
REAL            :: epsilon   = mh2o/mair ! = Mw / Ma, the
                                         !   ratio of the molecular weight of
                                         !   water to that of air 

REAL, PARAMETER :: q_mixratio_to_ppmv  = 1.60771704e+6

END MODULE NWPSAFMod_Constants
