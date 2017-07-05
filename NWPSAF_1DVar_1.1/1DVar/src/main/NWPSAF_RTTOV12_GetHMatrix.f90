SUBROUTINE NWPSAF_RTTOV12_GetHMatrix( &
  Profiles,        & !in
  Profiles_K,      & !in
  Profiles_K_PC,   & !in
  Instrument,      & !in
  NumInstChans,    & !in
  First_Chan_Pos,  & !in
  Last_Chan_Pos,   & !in
  RT_Params        ) !inout
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
! Description: Interface between IASI 1DVar code and RTTOV10
!
! History:
!
! Ticket  Date     Comment
! ------- -------- -------
!         11/04/12 New. Code taken from IASI_RTTOV10_Allocate. Fiona Smith
!
! Code Description:
!   Language:   Fortran 90
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Modules used:


USE NWPSAFMod_Params, ONLY : &
     Retrieved_Elements, &
     Ret_FirstQ,         &
     Ret_LastQ,          &
     Ret_Q2,             &
     Retrieve_QTotal,    &
     Retrieve_LWP,       &
     MwClwRetrieval,     &
     UsePCs,             &
     NPCScores

USE NWPSAFMod_RTmodel, ONLY : &
     RTParams_Type,         &
     ProfSize,              &
     Prof_FirstT,           &
     Prof_LastT,            &
     Prof_FirstQ,           &
     Prof_LastQ,            &
     Prof_FirstO3,          &
     Prof_LastO3,           &
     Prof_LastCLW,          &
     Prof_T2,               &
     Prof_Q2,               &
     Prof_VWind,            &
     Prof_UWind,            &
     Prof_PStar,            &
     Prof_CloudCover,       &
     Prof_LWP,              &
     Prof_CTP,              &
     Prof_Tstar,            &
     Num_RTlevels,          &
     Num_WetLevels,         &
     Num_ProfElementsUsed,  &
     RT_opts

USE NWPSAFMod_Constants, ONLY :  &
    epsilon

USE rttov_types, ONLY: &
    rttov_profile

USE NWPSAFMod_LiquidWater, ONLY: &
    LayerK_to_LWPK, &
    SVP, &
    SVP_Deriv

IMPLICIT NONE

INCLUDE 'NWPSAF_CloudStructure.interface'
INCLUDE 'NWPSAF_Qtot_to_q_ql.interface'

!Subroutine Arguments
TYPE(rttov_profile),  INTENT(IN)    :: Profiles(:)
TYPE(rttov_profile),  INTENT(IN)    :: Profiles_K(:)
TYPE(rttov_profile),  INTENT(IN)    :: Profiles_K_PC(:)
INTEGER,             INTENT(IN)    :: Instrument
INTEGER,             INTENT(IN)    :: NumInstChans
INTEGER,             INTENT(IN)    :: First_Chan_Pos
INTEGER,             INTENT(IN)    :: Last_Chan_Pos
TYPE(RTParams_Type), INTENT(INOUT) :: RT_Params

! Local constants
INTEGER, PARAMETER :: err_unit = 0
CHARACTER (LEN=*), PARAMETER :: RoutineName = "NWPSAF_RTTOV11_GetHMatrix"


! Local variables
INTEGER            :: i,j
REAL, ALLOCATABLE  :: RTModel_Jacobian (:,:)

! Variables used for cloud liquid water retrieval
INTEGER         :: QtotalOption
REAL            :: cloud_jac(Num_RTlevels)
REAL            :: cloud_structure(Num_RTlevels)
REAL            :: wt(Num_WetLevels)
REAL            :: wqtotal(Num_WetLevels)
REAL            :: wq(Num_WetLevels)
REAL            :: wql(Num_WetLevels)
REAL            :: wpress(Num_WetLevels)
REAL            :: esat(Num_WetLevels)
REAL            :: SESAT
REAL            :: SDlnes_DT
REAL            :: Dlnes_DT(Num_WetLevels)
REAL            :: qsat(Num_WetLevels)
REAL            :: Dqsat_dT(Num_WetLevels)
REAL            :: delta_eps
REAL            :: RTModel_Jacobian_Add(ProfSize, NumInstChans)
!--------------------------------------------------------------------------------


!----------------------------------------------------------------------------
! 1) Allocate Jacobian array
!----------------------------------------------------------------------------
IF ( UsePCs ) THEN
  ALLOCATE(RTModel_Jacobian(ProfSize, NPCScores))
ELSE
  ALLOCATE(RTModel_Jacobian(ProfSize, NumInstChans))
END IF

!2) Unpack the output K arrays from RTTOV into RTModel_Jacobian
!------

!2.1) PC Retrieval
!--------
PCRetrievalSwitch: IF (UsePCs) THEN

  ! First do the atmospheric profile variables:
  !----
  DO i=1,NPCSCores
    RTModel_Jacobian(Prof_FirstT:Prof_LastT, i) = &
          profiles_k_PC(i)%t(:)
    RTModel_Jacobian(Prof_FirstQ:Prof_LastQ, i) = &
          profiles_k_PC(i)%q(:)
    IF(RT_opts(Instrument)% rt_ir  % ozone_data) THEN
      RTModel_Jacobian(Prof_FirstO3:Prof_LastO3, i) = &
          profiles_k_PC(i)%o3(:)
    END IF
  END DO

  ! Now do the Surface Variables
  !-----
  RTModel_Jacobian(Prof_T2,1:NPCSCores) = &
        profiles_k_PC(1:NPCScores)%s2m%t
  RTModel_Jacobian(Prof_Q2,1:NPCSCores) = &
        profiles_k_PC(1:NPCScores)%s2m%q
  RTModel_Jacobian(Prof_PStar,1:NPCSCores) = &
        profiles_k_PC(1:NPCScores)%s2m%p
  RTModel_Jacobian(Prof_UWind,1:NPCSCores) = &
        profiles_k_PC(1:NPCSCores)%s2m%u
  RTModel_Jacobian(Prof_VWind,1:NPCSCores) = &
        profiles_k_PC(1:NPCSCores)%s2m%v

  ! Now do the Skin Variables
  !----
  RTModel_Jacobian(Prof_TStar,1:NPCSCores) = &
        profiles_k_PC(1:NPCScores)%skin%t

!2.2) Normal Retrieval
!--------
ELSE PCRetrievalSwitch

  ! First do the atmospheric profile variables:
  !----
  DO i=1,NumInstChans
    RTModel_Jacobian(Prof_FirstT:Prof_LastT, i) = &
          profiles_k(i)%t(:)
    RTModel_Jacobian(Prof_FirstQ:Prof_LastQ, i) = &
          profiles_k(i)%q(:)

    IF(RT_opts(Instrument) % rt_ir % ozone_data) THEN
      RTModel_Jacobian(Prof_FirstO3:Prof_LastO3, i) = &
            profiles_k(i)%o3(:)
    END IF
  END DO

  ! Cloud liquid water retrieval is actually handled as
  ! LWP or qtotal and is only for microwave instruments.
  ! Gradients transfer into CLW for a LWP retrieval but for qtotal
  ! The gradients move into the q and temperature parts of the
  ! structure

  IF (MwClwRetrieval .AND. Retrieve_LWP) THEN ! LWP retrieval

    !Generate cloud structure, and integrate up the clw profile jacobian
    CALL NWPSAF_CloudStructure( RT_params,         &
                                cloud_structure(:) )
    DO i=1,NumInstChans
      !this rather unnecesary step avoids a kind mismatch compile error with pgf90
      cloud_jac(:)=profiles_k(i)%clw(:)
      RTModel_Jacobian(Prof_LWP, i) = &
        LayerK_to_LWPK( cloud_jac(:), cloud_structure(:) )
    END DO

  ELSE IF (MwClwRetrieval .AND. Retrieve_Qtotal) THEN ! qtotal retrieval

    !First compute dBT/dQt derivative
    !-------------------------------
    wt(1:Num_WetLevels )= &
      RT_Params % RTguess(Prof_LastT - Num_WetLevels + 1:Prof_LastT)
    wpress(1:Num_WetLevels )= &
      RT_Params % Pressure_Pa(Prof_LastT - Num_WetLevels + 1:Prof_LastT)
    wqtotal(1:Num_WetLevels)= &
      EXP(RT_Params % RTguess( &
                        Prof_LastQ-Num_WetLevels + 1:Prof_LastQ)) + &
          RT_Params % RTguess( &
                        Prof_LastCLW- Num_WetLevels + 1:Prof_LastCLW)

    ! Compute derivatives  q = dq/dqtotal and ql = dql/dqtotal
    QtotalOption=0
    CALL NWPSAF_Qtot_to_q_ql( wqtotal(1:Num_WetLevels),&
                              wt(1:Num_WetLevels),     &
                              wpress(1:Num_WetLevels), &
                              QtotalOption,            &
                              wq(1:Num_WetLevels),     &
                              wql(1:Num_WetLevels)     )

    DO j=1,Num_WetLevels
      ! Kmat= dTB/dln(qtotal)=qtotal dTB/dqtotal
      DO i=1,NumInstChans
        RTModel_Jacobian( &
                    Prof_LastQ- Num_WetLevels + 1:Prof_LastQ,i) =   &
          profiles_k(i-First_Chan_Pos+1)% q(Prof_LastT -            &
                            Num_WetLevels+1:Prof_LastT)   * wq(j) + &
          profiles_k(i-First_Chan_Pos+1)% clw(Prof_LastT -          &
                            Num_WetLevels + 1:Prof_LastT) * wql(j)
      END DO
    END DO

    !Now compute effect on Temperature derivative
    !--------------------------------------------

    !Compute derivatives WRT qsat q = dq/dqsat  and ql = dql/dqsat
    QtotalOption=2

    CALL NWPSAF_Qtot_to_q_ql( wqtotal(1:Num_WetLevels),&
                              wt(1:Num_WetLevels),     &
                              wpress(1:Num_WetLevels), &
                              QtotalOption,            &
                              wq(1:Num_WetLevels),     &
                              wql(1:Num_WetLevels)     )


    delta_eps=1./epsilon -1.
    DO i=1,Num_WetLevels
      Sesat = svp(wt(i)) ! SESAT is in Pa.
      esat(i) = SESAT
      SDlnes_DT = svp_deriv(wt(i))
      Dlnes_DT(i) = SDlnes_DT
    ENDDO

    WHERE (esat(:) > wPress(:))
      qsat(:)     = 1.0
      Dqsat_DT(:) = 0.0
    ELSEWHERE
      qsat(:)     = epsilon / ( wpress(:)/esat(:) - (1.-epsilon))
      Dqsat_dT(:) = qsat(:) * Dlnes_DT(:) *(1.+delta_eps*qsat(:))
    ENDWHERE

    RTModel_Jacobian_Add(:, :) = 0.0

    DO j=1,Num_WetLevels
      DO i=1,NumInstChans
        RTModel_Jacobian_Add( &
                      Prof_LastQ-Num_WetLevels+1:Prof_LastQ,i) =    &
          profiles_k(i)%q(Prof_LastT -             &
                Num_WetLevels+1:Prof_LastT) *wq(j)  * dqsat_dT(j) + &
          profiles_k(i)%clw(Prof_LastT -           &
                Num_WetLevels+1:Prof_LastT) *wql(j) * dqsat_dT(j)
      ENDDO
    ENDDO

    ! Add this to the temperature jacobian
    DO i=1,NumInstChans
      RTModel_Jacobian( Prof_LastT- Num_WetLevels + 1:Prof_LastT, i) = &
        RTModel_Jacobian( Prof_LastT-Num_WetLevels+1:Prof_LastT, i) +    &
        RTModel_Jacobian_Add(Prof_LastQ-Num_WetLevels+1:Prof_LastQ, i)
    ENDDO
  ENDIF !qtotal/LWP retrieval


  ! Now do the Surface Variables
  !----
  RTModel_Jacobian(Prof_T2,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%s2m%t
  RTModel_Jacobian(Prof_Q2,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%s2m%q
  RTModel_Jacobian(Prof_PStar,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%s2m%p
  RTModel_Jacobian(Prof_UWind,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%s2m%u
  RTModel_Jacobian(Prof_VWind,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%s2m%v

  ! Now do the Skin Variables
  !----
  RTModel_Jacobian(Prof_TStar,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%skin%t

  ! Now do the Cloud Variables
  !----
  RTModel_Jacobian(Prof_CTP,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%ctp
  RTModel_Jacobian(Prof_CloudCover,1:NumInstChans) = &
        profiles_k(1:NumInstChans)%cfraction

END IF PCRetrievalSwitch


! 3) Pick out the Jacobian elements to be used by the minimisation
!-------
RT_Params % H_matrix_T(1:Num_ProfElementsUsed,First_Chan_Pos:Last_Chan_Pos) = &
   RTModel_Jacobian(Retrieved_Elements,:)


! Water Vapour Jacobians must be converted from
! ppmv to kg/kg to log(kg/kg) - the unit conversion cancels, then:
! dy/d(ln q) = dy/dq * q/dy   (dy = 1K)
!-----------------------------------------------------------------
IF (Ret_FirstQ > 0) THEN
  DO I = Ret_FirstQ, Ret_LastQ
    IF (Retrieved_Elements(I) > 0) &
      RT_Params % H_matrix_T(I,First_Chan_Pos:Last_Chan_Pos) = RT_Params % H_matrix_T(I,First_Chan_Pos:Last_Chan_Pos) * &
        profiles(1)%q(Retrieved_Elements(I)-Prof_FirstQ+1)
  END DO
END IF

! This is the surface humidity Jacobian
IF (Ret_q2 > 0) THEN
  IF (Retrieved_Elements(Ret_q2) > 0) &
      RT_Params % H_matrix_T(Ret_q2,First_Chan_Pos:Last_Chan_Pos) = &
        RT_Params % H_matrix_T(Ret_q2,First_Chan_Pos:Last_Chan_Pos) * profiles(1)%s2m%q
END IF
RT_Params % H_matrix = TRANSPOSE(RT_Params % H_matrix_T)

! 4) Deallocate Jacobian array
!-----
DEALLOCATE(RTModel_Jacobian)


END SUBROUTINE NWPSAF_RTTOV12_GetHMatrix