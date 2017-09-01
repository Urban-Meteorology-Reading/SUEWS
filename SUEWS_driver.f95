!========================================================================================
SUBROUTINE SUEWS_cal_Qn(&
     nsurf,NetRadiationMethod,snowUse,ldown_option,id,&!input
     DecidSurf,ConifSurf,GrassSurf,Diagnose,&
     snow_obs,ldown_obs,fcld_obs,&
     dectime,ZENITH_deg,avKdn,Temp_C,avRH,Press_hPa,qn1_obs,&
     SnowAlb,&
     AlbedoChoice,DiagQN,&
     alb,albDecTr,DecidCap,albEveTr,albGrass,surf,&!inout
     snowFrac,ldown,fcld,&!output
     qn1,qn1_SF,qn1_S,kclear,kup,lup,TSURF)



  IMPLICIT NONE

  INTEGER,INTENT(in)::nsurf
  INTEGER,INTENT(in)::NetRadiationMethod
  INTEGER,INTENT(in)::snowUse
  INTEGER,INTENT(in)::ldown_option
  INTEGER,INTENT(in)::id
  INTEGER,INTENT(in)::DecidSurf
  INTEGER,INTENT(in)::ConifSurf
  INTEGER,INTENT(in)::GrassSurf
  INTEGER,INTENT(in)::Diagnose
  REAL(KIND(1d0)),INTENT(in)::snow_obs
  REAL(KIND(1d0)),INTENT(in)::ldown_obs
  REAL(KIND(1d0)),INTENT(in)::fcld_obs
  REAL(KIND(1d0)),INTENT(in)::dectime
  REAL(KIND(1d0)),INTENT(in)::ZENITH_deg
  REAL(KIND(1d0)),INTENT(in)::avKdn
  REAL(KIND(1d0)),INTENT(in)::Temp_C
  REAL(KIND(1d0)),INTENT(in)::avRH
  REAL(KIND(1d0)),INTENT(in)::Press_hPa
  REAL(KIND(1d0)),INTENT(in)::qn1_obs
  REAL(KIND(1d0)),INTENT(in)::SnowAlb
  REAL(KIND(1d0)),INTENT(in)::AlbedoChoice
  REAL(KIND(1d0)),INTENT(in)::DiagQN

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout):: alb
  REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)::albDecTr
  REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)::DecidCap
  REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)::albEveTr
  REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)::albGrass
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(inout):: surf

  REAL(KIND(1d0)),INTENT(out)::snowFrac
  REAL(KIND(1d0)),INTENT(out)::ldown
  REAL(KIND(1d0)),INTENT(out)::fcld
  REAL(KIND(1d0)),INTENT(out)::qn1
  REAL(KIND(1d0)),INTENT(out)::qn1_SF
  REAL(KIND(1d0)),INTENT(out)::qn1_S
  REAL(KIND(1d0)),INTENT(out)::kclear
  REAL(KIND(1d0)),INTENT(out)::kup
  REAL(KIND(1d0)),INTENT(out)::lup
  REAL(KIND(1d0)),INTENT(out)::TSURF


  IF(NetRadiationMethod>0)THEN

     IF (snowUse==0) snowFrac=snow_obs

     IF(ldown_option==1) THEN !Observed ldown provided as forcing
        ldown=ldown_obs
     ELSE
        ldown=-9              !to be filled in NARP
     ENDIF

     IF(ldown_option==2) THEN !observed cloud fraction provided as forcing
        fcld=fcld_obs
     ENDIF

     !write(*,*) DecidCap(id), id, it, imin, 'Calc - near start'

     ! Update variables that change daily and represent seasonal variability
     alb(DecidSurf)    = albDecTr(id) !Change deciduous albedo
     surf(6,DecidSurf) = DecidCap(id) !Change current storage capacity of deciduous trees
     ! Change EveTr and Grass albedo too
     alb(ConifSurf) = albEveTr(id)
     alb(GrassSurf) = albGrass(id)

     IF(Diagnose==1) WRITE(*,*) 'Calling NARP...'
     CALL NARP(&
                                ! input:
          dectime,ZENITH_deg,avKdn,Temp_C,avRH,Press_hPa,qn1_obs,&
          SnowAlb,&
          AlbedoChoice,ldown_option,NetRadiationMethod,DiagQN,&
                                ! output:
          qn1,qn1_SF,qn1_S,kclear,kup,LDown,lup,fcld,TSURF)
     !Temp_C,kclear,fcld,dectime,avkdn,avRH,qn1,kup,ldown,lup,tsurf,&
     !AlbedoChoice,ldown_option,Press_hPa,Ea_hPa,qn1_obs,&
     !zenith_deg,NetRadiationMethod,
  ELSE
     snowFrac = snow_obs
     qn1      = qn1_obs
     qn1_sf   = qn1_obs
     qn1_s    = qn1_obs
  ENDIF

END SUBROUTINE SUEWS_cal_Qn
!========================================================================================

!========================================================================================
SUBROUTINE SUEWS_cal_Qs(&
     nsurf,&
     StorageHeatMethod,&
     OHMIncQF,&
     Gridiv,&
     id,&
     Diagnose,&
     sfr,&
     OHM_coef,&
     OHM_threshSW,OHM_threshWD,&
     soilmoist,&
     soilstoreCap,&
     state,&
     nsh,&
     BldgSurf,&
     WaterSurf,&
     SnowUse,&
     DiagQS,&
     HDDday,&
     MetForcingData_grid,&
     qn1,qn1_bup,&
     alb,&
     emis,&
     cp,&
     kk,&
     ch,&
     AnthropHeatMethod,&
     qn1_store,&
     qn1_S_store,&
     qn1_av_store,&
     qn1_S_av_store,&
     surf,&
     qn1_S,&
     snowFrac,&
     qs,&
     deltaQi,&
     a1,&
     a2,&
     a3&
     )


  IMPLICIT NONE

  INTEGER,INTENT(in)::nsurf
  INTEGER,INTENT(in)::StorageHeatMethod
  INTEGER,INTENT(in)::OHMIncQF
  INTEGER,INTENT(in)::Gridiv
  INTEGER,INTENT(in)::id
  INTEGER,INTENT(in)::Diagnose
  INTEGER,INTENT(in)::nsh       ! number of timesteps in one hour
  INTEGER,INTENT(in)::BldgSurf  ! code for specific surfaces
  INTEGER,INTENT(in)::WaterSurf ! code for specific surfaces
  INTEGER,INTENT(in)::SnowUse   ! option for snow related calculations
  INTEGER,INTENT(in)::DiagQS    ! diagnostic option
  INTEGER,INTENT(in):: AnthropHeatMethod !< AnthropHeat option [-]


  REAL(KIND(1d0)),INTENT(in)::OHM_coef(9,4,3)                 ! OHM coefficients
  REAL(KIND(1d0)),INTENT(in)::OHM_threshSW(9),OHM_threshWD(9) ! OHM thresholds
  REAL(KIND(1d0)),INTENT(in)::soilmoist(nsurf)                ! soil moisture
  REAL(KIND(1d0)),INTENT(in)::soilstoreCap(nsurf)             ! capacity of soil store
  REAL(KIND(1d0)),INTENT(in)::state(nsurf) ! wetness status


  REAL(KIND(1d0)),INTENT(in)::HDDday  ! HDDday=HDD(id-1,4) HDD at the begining of today (id-1)
  REAL(KIND(1d0)),INTENT(in)::qn1
  REAL(KIND(1d0)),INTENT(in)::qn1_bup

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
  REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::MetForcingData_grid !< met forcing array of grid
  REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::alb  !< albedo [-]
  REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::emis !< emissivity [-]
  REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::cp   !< heat capacity [J m-3 K-1]
  REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::kk   !< thermal conductivity [W m-1 K-1]
  REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::ch   !< bulk transfer coef [J m-3 K-1]


  REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout) ::qn1_store
  REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout) ::qn1_S_store !< stored qn1 [W m-2]

  REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_av_store
  REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_S_av_store !< average net radiation over previous hour [W m-2]
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(inout)::surf
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(out)::deltaQi(nsurf+2) ! storage heat flux of snow surfaces

  REAL(KIND(1d0)),INTENT(out)::qn1_S
  REAL(KIND(1d0)),INTENT(out)::snowFrac
  REAL(KIND(1d0)),INTENT(out):: qs ! storage heat flux

  REAL(KIND(1d0)),INTENT(out):: a1 !< AnOHM coefficients of grid [-]
  REAL(KIND(1d0)),INTENT(out):: a2 !< AnOHM coefficients of grid [h]
  REAL(KIND(1d0)),INTENT(out):: a3 !< AnOHM coefficients of grid [W m-2]


  IF(StorageHeatMethod==1) THEN           !Use OHM to calculate QS
     IF(OHMIncQF == 1) THEN      !Calculate QS using QSTAR+QF
        IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
        CALL OHM(qn1,qn1_store,qn1_av_store,&
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
     ELSEIF(OHMIncQF == 0) THEN  !Calculate QS using QSTAR
        ! qn1=qn1_bup
        IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
        CALL OHM(qn1_bup,qn1_store,qn1_av_store,&
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
     ENDIF
  ENDIF

  ! use AnOHM to calculate QS, TS 14 Mar 2016
  IF (StorageHeatMethod==3) THEN
     IF ( OHMIncQF == 1 ) THEN    !Calculate QS using QSTAR+QF
        IF(Diagnose==1) WRITE(*,*) 'Calling AnOHM...'
        CALL AnOHM(qn1,qn1_store,qn1_av_store,&
             MetForcingData_grid,state/surf(6,:),&
             alb, emis, cp, kk, ch,&
             sfr,nsurf,nsh,AnthropHeatMethod,id,Gridiv,&
             a1,a2,a3,qs)
     ELSEIF(OHMIncQF == 0) THEN   !Calculate QS using QSTAR
        ! qn1=qn1_bup
        IF(Diagnose==1) WRITE(*,*) 'Calling AnOHM...'
        CALL AnOHM(qn1_bup,qn1_store,qn1_av_store,&
             MetForcingData_grid,state/surf(6,:),&
             alb, emis, cp, kk, ch,&
             sfr,nsurf,nsh,AnthropHeatMethod,id,Gridiv,&
             a1,a2,a3,qs)
     END IF
  END IF

  ! !Calculate QS using ESTM
  ! IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
  !    !CALL ESTM_v2016(QSestm,iMB)
  !    IF(Diagnose==1) WRITE(*,*) 'Calling ESTM...'
  !    CALL ESTM_v2016(QSestm,Gridiv,ir)  ! iMB corrected to Gridiv, TS 09 Jun 2016
  !    QS=QSestm   ! Use ESTM qs
  ! ENDIF

END SUBROUTINE SUEWS_cal_Qs
