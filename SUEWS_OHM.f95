MODULE OHM_module
  ! USE allocateArray
  ! USE data_in
  ! USE defaultNotUsed
  ! USE gis_data
  ! USE sues_data
  ! USE time

  IMPLICIT NONE
CONTAINS
  !========================================================================================
  SUBROUTINE OHM(qn1,qn1_store,qn1_av_store,&
       qn1_S,qn1_S_store,qn1_S_av_store,&
       nsh,&
       sfr,nsurf,&
       HDDday,&
       OHM_coef,&
       OHM_threshSW,OHM_threshWD,&
       soilmoist,soilstoreCap,state,&
       BldgSurf,WaterSurf,&
       SnowUse,SnowFrac,&
       DiagQS,&
       qs,deltaQi)
    ! Made by HCW Jan 2015 to replace OHMnew (no longer needed).
    ! Calculates net storage heat flux (QS) from Eq 4, Grimmond et al. 1991, Atm Env.
    ! Accounts for variable timesteps in dQ*/dt term.
    ! BareSoilSurfFraction removed so bare soil is now handled like other surfaces.
    ! Snow part changed from summer wet to winter wet coefficients.
    ! Changed -333 checks to -999 checks and added error handling
    ! Gradient now calculated for t-1 (was previously calculated for t-2).
    ! Modified by TS 07 Aug 2017
    !  1. interface changed to account for explict passing
    !  2. calculation refactorization.
    ! Modified by HCW 25 Feb 2015
    !  Adapted q1,q2,q3 & r1,r2,r3 for multiple grids
    ! Modified by HCW 14 Dec 2016
    !  Thresholds for Summer/Winter and Wet/Dry now provided in input files
    !  Calculation of dqndt now uses hourly running mean rather than instantaneous values
    ! To Do:
    !   - No canyons implemented at the moment [OHM_coef(nsurf+1,,)]
    !========================================================================================


    IMPLICIT NONE

    REAL(KIND(1d0)),INTENT(in)::&
         qn1,& ! net all-wave radiation
         qn1_S,& !  net all-wave radiation over snow
         sfr(nsurf),& ! surface fractions
         SnowFrac(nsurf),& ! snow fractions of each surface
         HDDday,& ! HDDday=HDD(id-1,4) HDD at the begining of today (id-1)
         OHM_coef(9,4,3),& ! OHM coefficients
         OHM_threshSW(9),OHM_threshWD(9),& ! OHM thresholds
         soilmoist(nsurf),& ! soil moisture
         soilstoreCap(nsurf),&! capacity of soil store
         state(nsurf) ! wetness status
    INTEGER,INTENT(in)::&
         nsurf,& ! number of surfaces
         nsh,& ! number of timesteps in one hour
         BldgSurf,WaterSurf,& ! code for specific surfaces
         SnowUse,& ! option for snow related calculations
         DiagQS ! diagnostic option
    REAL(KIND(1d0)),INTENT(inout)::&
         qn1_store(nsh),qn1_av_store(2*nsh+1),&
         qn1_S_store(nsh),qn1_S_av_store(2*nsh+1)
    REAL(KIND(1d0)),INTENT(out):: &
         qs,& ! storage heat flux
         deltaQi(nsurf+2) ! storage heat flux of snow surfaces

    REAL(KIND(1d0)):: a1,a2,a3 ! OHM coefficients of grid

    ! REAL(KIND(1d0)):: nsh_nna ! number of timesteps per hour with non -999 values (used for spinup)

    REAL(KIND(1d0)):: dqndt    !Rate of change of net radiation [W m-2 h-1] at t-1
    ! REAL(KIND(1d0)):: surfrac  !Surface fraction accounting for SnowFrac if appropriate

    ! REAL(KIND(1d0)):: qn1_av, qn1_S_av    !Average net radiation over previous hour [W m-2]
    REAL(KIND(1d0)):: deltaQi0 ! temporarily store

    ! REAL(KIND(1d0)):: qn1_store0(nsh), qn1_av_store0(2*nsh+1) ! temporarily store

    !These are now provided in SiteInfo (OHMthresh for Summer/Winter and Wet/Dry)
    !!real(kind(1d0)):: OHM_TForSummer = 5  !Use summer coefficients if 5-day Tair >= 5 degC
    !real(kind(1d0)):: OHM_TForSummer = 10  !Use summer coefficients if 5-day Tair >= 10 degC - modified for UK HCW 14 Dec 2015
    !real(kind(1d0)):: OHM_SMForWet = 0.9  !Use wet coefficients if SM close to soil capacity

    CALL OHM_coef_cal(sfr,nsurf,&
         HDDday,OHM_coef,OHM_threshSW,OHM_threshWD,&
         soilmoist,soilstoreCap,state,&
         BldgSurf,WaterSurf,&
         SnowUse,SnowFrac,&
         a1,a2,a3)
    ! WRITE(*,*) '----- OHM coeffs new-----'
    ! WRITE(*,*) a1,a2,a3


    ! WRITE(*,*) '----- OHM coeffs -----'
    ! WRITE(*,*) a1,a2,a3

    ! Old OHM calculations (up to v2016a)
    !! Calculate radiation part ------------------------------------------------------------
    !qs=NAN              !qs  = Net storage heat flux  [W m-2]
    !if(qn1>-999) then   !qn1 = Net all-wave radiation [W m-2]
    !   !if(q1>-999.and.q3>-999) then
    !      !dqndt = 0.5*(q3-q1)*nsh_real                !gradient at t-2
    !      dqndt = 0.5*(qn1-q2_grids(Gridiv))*nsh_real   !gradient at t-1
    !
    !      !Calculate net storage heat flux
    !      qs = qn1*a1 + dqndt*a2 + a3   !Eq 4, Grimmond et al. 1991
    !   !endif
    !   !q1=q2  !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
    !   !q2=q3  !q2 = net radiation at t-1
    !   !q3=qn1  !q3 = net radiation at t   (at t-1 when q3 used in next timestep)
    !   q1_grids(Gridiv) = q2_grids(Gridiv) !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
    !   q2_grids(Gridiv) = q3_grids(Gridiv) !q2 = net radiation at t-1
    !   q3_grids(Gridiv) = qn1              !q3 = net radiation at t (at t-1 when q3 used in next timestep)
    !else
    !   call ErrorHint(21,'Bad value for qn1 found during OHM calculation',qn1,NotUsed,notUsedI)
    !endif


    ! New OHM calculations (v2017a onwards) using running mean (HCW Dec 2016)
    ! Calculate radiation part ------------------------------------------------------------
    qs=-999              !qs  = Net storage heat flux  [W m-2]
    IF(qn1>-999) THEN   !qn1 = Net all-wave radiation [W m-2]
       ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
       CALL OHM_dqndt_cal(nsh,qn1,qn1_store,qn1_av_store,dqndt)

       ! Calculate net storage heat flux
       CALL OHM_QS_cal(qn1,dqndt,a1,a2,a3,qs)
       IF(DiagQS==1) WRITE(*,*) 'qs: ',qs,'qn1:',qn1,'dqndt: ',dqndt

    ELSE
       CALL ErrorHint(21,'In SUEWS_OHM.f95: bad value for qn found during qs calculation.',qn1,-55.55,-55)
    ENDIF

    !write(*,*) qs
    !write(*,*) '--------------------'

    ! Do snow calculations separately -----
    ! Added by LJ in August 2013
    IF(snowUse==1) THEN
       deltaQi=-999
       IF(qn1_S>-999) THEN
          ! Old OHM calculations (commented out HCW Dec 2016)
          !!if(r1>-999.and.r3>-999) then
          !   !dqndt = 0.5*(r3-r1)*nsh_real    !gradient at t-2
          !   dqndt = 0.5*(qn1_S-r2_grids(Gridiv))*nsh_real     !gradient at t-1
          !   ! Calculate net storage heat flux for snow surface (winter wet conditions HCW 15/01/2015)
          !   deltaQi = qn1_S*OHM_coef(nsurf+2,3,1) + dqndt*OHM_coef(nsurf+2,3,2) + OHM_coef(nsurf+2,3,3)
          !!endif
          !r1_grids(Gridiv)=r2_grids(Gridiv)
          !r2_grids(Gridiv)=r3_grids(Gridiv)
          !r3_grids(Gridiv)=qn1_S
          ! New OHM calculations
          ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
          CALL OHM_dqndt_cal(nsh,qn1_S,qn1_S_store,qn1_S_av_store,dqndt)

          ! Calculate net storage heat flux for snow surface (winter wet conditions)
          CALL OHM_QS_cal(qn1_S,dqndt,&
               OHM_coef(nsurf+2,3,1),OHM_coef(nsurf+2,3,2),OHM_coef(nsurf+2,3,3),&
               deltaQi0)
          deltaQi=deltaQi0


       ELSE
          CALL ErrorHint(21,'In SUEWS_OHM.f95: bad value for qn(snow) found during qs calculation.',qn1_S,-55.55,-55)
       ENDIF

    ENDIF

    RETURN
  ENDSUBROUTINE OHM
  !========================================================================================

  SUBROUTINE OHM_coef_cal(sfr,nsurf,&
       HDDday,OHM_coef,OHM_threshSW,OHM_threshWD,&
       soilmoist,soilstoreCap,state,&
       BldgSurf,WaterSurf,&
       SnowUse,SnowFrac,&
       a1,a2,a3)
    IMPLICIT NONE
    INTEGER , INTENT(in) :: &
         nsurf,& ! number of surfaces
         SnowUse,& ! option for snow related calculations
         BldgSurf,WaterSurf ! code for specific surfaces
    REAL(KIND(1d0)), INTENT(in) :: &
         sfr(nsurf),& ! surface cover fractions
         SnowFrac(nsurf),& ! snow fractions of each surface
         HDDday,& ! HDDday=HDD(id-1,4) HDD at the begining of today (id-1)
         OHM_coef(9,4,3),&
         OHM_threshSW(9),OHM_threshWD(9),& ! OHM thresholds
         soilmoist(nsurf),& ! soil moisture
         soilstoreCap(nsurf),&! capacity of soil store
         state(nsurf) ! wetness status
    REAL(KIND(1d0)), INTENT(out):: a1,a2,a3

    REAL(KIND(1d0)) :: surfrac
    INTEGER :: i,ii,is

    ! OHM coefficients --------
    ! Set to zero initially
    a1=0   ![-]
    a2=0   ![h]
    a3=0   ![W m-2]
    ! -------------------------

    ! Loop through surface types ----------------------------------------------------------
    DO is=1,nsurf
       surfrac=sfr(is)

       ! Use 5-day running mean Tair to decide whether it is summer or winter ----------------
       IF(HDDday >= OHM_threshSW(is)) THEN !Summer
          ii=0
       ELSE          !Winter
          ii=2
       ENDIF

       IF(state(is) > 0) THEN     !Wet surface
          i=ii+1
       ELSE                    !Dry surface
          i=ii+2
          ! If the surface is dry but SM is close to capacity, use coefficients for wet surfaces
          IF(is>BldgSurf.AND.is/=WaterSurf)THEN    !Wet soil (i.e. EveTr, DecTr, Grass, BSoil surfaces)
             IF(soilmoist(is)/soilstoreCap(is) > OHM_threshWD(is) ) THEN
                i=ii+1
             ENDIF
          ENDIF
       ENDIF

       ! If snow, adjust surface fractions accordingly
       IF(SnowUse==1.AND.is/=BldgSurf.AND.is/=WaterSurf) THEN   ! ?? Why is BldgSurf excluded here?
          surfrac=surfrac*(1-SnowFrac(is))
       ENDIF

       ! Calculate the areally-weighted OHM coefficients
       a1 = a1+surfrac*OHM_coef(is,i,1)
       a2 = a2+surfrac*OHM_coef(is,i,2)
       a3 = a3+surfrac*OHM_coef(is,i,3)

    ENDDO  !end of loop over surface types ------------------------------------------------
  END SUBROUTINE OHM_coef_cal

  SUBROUTINE OHM_dqndt_cal(nsh,qn1,qn1_store,qn1_av_store,dqndt)
    IMPLICIT NONE
    INTEGER, INTENT(in)            :: nsh ! number of timesteps in one hour
    REAL(KIND(1d0)), INTENT(in) ::qn1
    REAL(KIND(1d0)), INTENT(inout) :: &
         qn1_store(nsh),& ! instantaneous qn1 values for previous hour
         qn1_av_store(2*nsh+1) ! average qn1 values for previous hour
    REAL(KIND(1d0)), INTENT(out)   :: dqndt !dQ* per dt for 60 min

    REAL(KIND(1d0)) :: qn1_av
    INTEGER :: nsh_nna

    ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
    IF(nsh > 1) THEN
       qn1_store=CSHIFT(qn1_store,1) ! shift to left with one place
       qn1_store(nsh)=qn1
       nsh_nna = COUNT(qn1_store/=-999, dim=1) !Find how many are not -999s  !bug fixed HCW 08 Feb 2017
       qn1_av = SUM(qn1_store, mask=qn1_store/= -999)/nsh_nna
    ELSEIF(nsh==1) THEN
       qn1_store(:) = qn1
       qn1_av = qn1
    ENDIF
    ! Store hourly average values (calculated every timestep) for previous 2 hours
    IF(nsh > 1) THEN
       qn1_av_store=CSHIFT(qn1_av_store,1)
       qn1_av_store(2*nsh+1) = qn1_av
    ELSEIF(nsh==1) THEN
       qn1_av_store(:) = qn1_av
    ENDIF
    ! Calculate dQ* per dt for 60 min (using running mean Q* at t hours and (t-2) hours)
    IF(ANY(qn1_av_store == -999)) THEN
       dqndt=0  ! Set dqndt term to zero for spinup
    ELSE
       dqndt=0.5*(qn1_av_store((2*nsh+1))-qn1_av_store(1))
    ENDIF

  END SUBROUTINE OHM_dqndt_cal

  SUBROUTINE OHM_QS_cal(qn1,dqndt,a1,a2,a3,qs)
    IMPLICIT NONE
    REAL(KIND(1d0)), INTENT(in) :: qn1,dqndt,a1,a2,a3
    REAL(KIND(1d0)), INTENT(out):: qs
    qs = qn1*a1 + dqndt*a2 + a3   !Eq 4, Grimmond et al. 1991

  END SUBROUTINE OHM_QS_cal


END MODULE OHM_module
