!+ Split total water content (qtotal) into q and ql
Subroutine NWPSAF_Qtot_to_q_ql(&
     qtotal,        & !In
     t,             & !In
     wpress,        & !In
     QtotOption,    & !In
     q,             & !Out
     ql             ) !Out
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
! IF QtotOption=1 : Split total water content (qtotal) into 
!   water vapor content (q) and (ql) cloud liquid water content.
! If QtotOption=0 : Compute derivatives: (q) =dq/dqtotal and (ql) = dql/dqtotal
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
!   28    24/02/12 Original code based on SSMI_Qtot_to_q_ql.f90  TR Sreerekha.
!-----------------------------------------------------------------------
!
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:

USE NWPSAFMod_RTmodel, ONLY : &
    Num_WetLevels

USE NWPSAFMod_Constants, ONLY :  & 
    epsilon,       &
    RH_1,                   &
    RH_2,                   &
    SplitQtotal

USE NWPSAFMod_LiquidWater, ONLY : &
    SVP

IMPLICIT NONE

!Subroutine Arguments
REAL, INTENT(IN)     :: qtotal(:) ! Num_WetLevels = Ret_LastQ - Ret_FirstQ + 1
REAL, INTENT(IN)     :: t(:)
REAL, INTENT(IN)     :: wpress(:)
INTEGER, INTENT(IN)  :: QtotOption
REAL, INTENT(OUT)    :: q(:)
REAL, INTENT(OUT)    :: ql(:)


!local
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_qtot_to_q_ql" 
INTEGER             ::  i ! loop index
REAL                :: qsat(Num_WetLevels)
REAL                :: RH_qtotal(Num_WetLevels)
REAL                :: qvfix(Num_WetLevels)
REAL                :: qexcess(Num_WetLevels)
REAL                :: esat(Num_WetLevels)
REAL                :: SESAT


!----------------End of header -----------------------------------------

q(:)=0.0
ql(:)=0.0

 DO i=1,Num_WetLevels
    SESAT=SVP(t(i)) ! SESAT is in Pa.
    esat(i)=SESAT
 ENDDO

 WHERE (esat(:) > wPress(:))
    qsat(:) = 1.0
 ELSEWHERE
    qsat(:)= epsilon / ( wpress(:)/esat(:) - (1.-epsilon))
 ENDWHERE

 RH_qtotal = qtotal / qsat

 IF (QtotOption == 1) THEN

    ! Split qtotal into q and ql 

    qvfix = qsat * (RH_1 + SplitQtotal * (RH_2 - RH_1))      
    qexcess = qtotal - RH_1 * qsat

    WHERE ( RH_qtotal < RH_1 ) 
       q  = qtotal
    ENDWHERE

    WHERE ( RH_qtotal >= RH_1 .AND. RH_qtotal < RH_2 )
        q= RH_1 *qsat + SplitQtotal* qexcess
    ENDWHERE

    WHERE ( RH_qtotal >= RH_2 ) 
       q=qvfix
    ENDWHERE

    ql = qtotal - q


 ELSE IF (QtotOption == 0) THEN
  
    ! Compute derivates
    ! q = dq/dqtotal and ql = dql/dqtotal

    WHERE ( RH_qtotal < RH_1 )
        q  = 1.0       
        ql = 0.0
    ENDWHERE
    WHERE ( RH_qtotal >= RH_1 .AND. RH_qtotal < RH_2 ) 
        q = SplitQtotal
        ql = (1. - SplitQtotal)
    ENDWHERE
    WHERE ( RH_qtotal >= RH_2 )
        q  = 0.0       
        ql = 1.0
    ENDWHERE

 ELSE IF (QtotOption == 2) THEN

    ! Compute derivatives of q and ql with respect to qsat
    !   q = dq/dqsat  and ql = dql/dqsat

    WHERE ( RH_qtotal(:) < RH_1 )
        q(:)  = 0.0
    ENDWHERE
    WHERE ( RH_qtotal(:) >= RH_1 .AND. RH_qtotal(:) < RH_2 )
        q(:) = (1. - SplitQtotal) *  RH_1
    ENDWHERE
    WHERE ( RH_qtotal(:) >= RH_2 )
        q(:)  = (RH_1 + SplitQtotal * (RH_2 - RH_1))
    ENDWHERE

    ql(:) = -q(:)

 ELSE

    WRITE(6,*) ' QtotOption not well defined'
    WRITE(6,*) ' QtotOption =', QtotOption
    STOP

 ENDIF

End Subroutine NWPSAF_Qtot_to_q_ql 
