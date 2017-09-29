!========================================================================================
! SUEWS driver subroutines
! TS 31 Aug 2017: initial version
MODULE SUEWS_Driver


  IMPLICIT NONE

CONTAINS
  !=============net all-wave radiation=====================================
  SUBROUTINE SUEWS_cal_Qn(&
       NetRadiationMethod,snowUse,ldown_option,id,&!input
       Diagnose,snow_obs,ldown_obs,fcld_obs,&
       dectime,ZENITH_deg,avKdn,Temp_C,avRH,Press_hPa,qn1_obs,&
       SnowAlb,AlbedoChoice,DiagQN,&
       NARP_G,NARP_TRANS_SITE,NARP_EMIS_SNOW,IceFrac,sfr,emis,&
       alb,albDecTr,DecidCap,albEveTr,albGrass,surf,&!inout
       snowFrac,ldown,fcld,&!output
       qn1,qn1_SF,qn1_S,kclear,kup,lup,tsurf)
    USE NARP_MODULE, ONLY: NARP

    IMPLICIT NONE
    INTEGER,PARAMETER ::nsurf     = 7 ! number of surface types
    INTEGER,PARAMETER ::ConifSurf = 3 !New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER ::DecidSurf = 4 !New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER ::GrassSurf = 5

    INTEGER,INTENT(in)::NetRadiationMethod
    INTEGER,INTENT(in)::snowUse
    INTEGER,INTENT(in)::ldown_option
    INTEGER,INTENT(in)::id
    ! INTEGER,INTENT(in)::DecidSurf
    ! INTEGER,INTENT(in)::ConifSurf
    ! INTEGER,INTENT(in)::GrassSurf
    INTEGER,INTENT(in)::Diagnose
    INTEGER,INTENT(in)::AlbedoChoice
    INTEGER,INTENT(in)::DiagQN

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
    REAL(KIND(1d0)),INTENT(in)::NARP_EMIS_SNOW
    REAL(KIND(1d0)),INTENT(in)::NARP_TRANS_SITE

    REAL(KIND(1d0)),DIMENSION(365),INTENT(in)::NARP_G
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in):: IceFrac
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in):: sfr
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in):: emis

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)  ::alb
    REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)  ::albDecTr
    REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)  ::DecidCap
    REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)  ::albEveTr
    REAL(KIND(1d0)),DIMENSION(0:366),INTENT(inout)  ::albGrass
    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(inout)::surf

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowFrac

    REAL(KIND(1d0)),INTENT(out)::ldown
    REAL(KIND(1d0)),INTENT(out)::fcld
    REAL(KIND(1d0)),INTENT(out)::qn1
    REAL(KIND(1d0)),INTENT(out)::qn1_SF
    REAL(KIND(1d0)),INTENT(out)::qn1_S
    REAL(KIND(1d0)),INTENT(out)::kclear
    REAL(KIND(1d0)),INTENT(out)::kup
    REAL(KIND(1d0)),INTENT(out)::lup
    REAL(KIND(1d0)),INTENT(out)::tsurf


    REAL(KIND(1d0)),DIMENSION(nsurf):: lup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf):: kup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf):: tsurf_ind
    REAL(KIND(1d0)),DIMENSION(nsurf):: qn1_ind

    REAL(KIND(1d0)),PARAMETER::NAN=-999


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
            nsurf,sfr,snowFrac,alb,emis,IceFrac,&! input:
            NARP_G,NARP_TRANS_SITE,NARP_EMIS_SNOW,&
            dectime,ZENITH_deg,avKdn,Temp_C,avRH,Press_hPa,qn1_obs,&
            SnowAlb,&
            AlbedoChoice,ldown_option,NetRadiationMethod,DiagQN,&
            qn1,qn1_SF,qn1_S,kclear,kup,LDown,lup,fcld,tsurf)! output:
       !Temp_C,kclear,fcld,dectime,avkdn,avRH,qn1,kup,ldown,lup,tsurf,&
       !AlbedoChoice,ldown_option,Press_hPa,Ea_hPa,qn1_obs,&
       !zenith_deg,NetRadiationMethod,
    ELSE ! NetRadiationMethod==0
       snowFrac  = snow_obs
       qn1       = qn1_obs
       qn1_sf    = qn1_obs
       qn1_s     = qn1_obs
       ldown     = NAN
       lup       = NAN
       kup       = NAN
       tsurf     = NAN
       lup_ind   = NAN
       kup_ind   = NAN
       tsurf_ind = NAN
       qn1_ind   = NAN
       Fcld      = NAN
    ENDIF

    IF(ldown_option==1) THEN
       Fcld = NAN
    ENDIF

  END SUBROUTINE SUEWS_cal_Qn
  !========================================================================

  !=============storage heat flux=========================================
  SUBROUTINE SUEWS_cal_Qs(&
       StorageHeatMethod,OHMIncQF,Gridiv,id,Diagnose,sfr,&!input
       OHM_coef,OHM_threshSW,OHM_threshWD,&
       soilmoist,soilstoreCap,state,nsh,SnowUse,DiagQS,&
       HDDday,MetForcingData_grid,qf,qn1_bup,alb,emis,cp,kk,ch,&
       AnthropHeatMethod,&
       qn1_store,qn1_S_store,qn1_av_store,qn1_S_av_store,&!inout
       surf,qn1_S,snowFrac,qs,&
       deltaQi,a1,a2,a3)!output
    USE AnOHM_module, ONLY: anohm

    IMPLICIT NONE

    INTEGER,PARAMETER :: nsurf    = 7      ! number of surface types
    INTEGER,PARAMETER ::BldgSurf  = 2      !New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER ::WaterSurf = 7

    INTEGER,INTENT(in)::StorageHeatMethod
    INTEGER,INTENT(in)::OHMIncQF
    INTEGER,INTENT(in)::Gridiv
    INTEGER,INTENT(in)::id
    INTEGER,INTENT(in)::Diagnose
    INTEGER,INTENT(in)::nsh                ! number of timesteps in one hour
    INTEGER,INTENT(in)::SnowUse            ! option for snow related calculations
    INTEGER,INTENT(in)::DiagQS             ! diagnostic option
    INTEGER,INTENT(in):: AnthropHeatMethod !< AnthropHeat option [-]


    REAL(KIND(1d0)),INTENT(in)::OHM_coef(9,4,3)                 ! OHM coefficients
    REAL(KIND(1d0)),INTENT(in)::OHM_threshSW(9) ! OHM thresholds
    REAL(KIND(1d0)),INTENT(in)::OHM_threshWD(9) ! OHM thresholds
    REAL(KIND(1d0)),INTENT(in)::soilmoist(nsurf)                ! soil moisture
    REAL(KIND(1d0)),INTENT(in)::soilstoreCap(nsurf)             ! capacity of soil store
    REAL(KIND(1d0)),INTENT(in)::state(nsurf) ! wetness status


    REAL(KIND(1d0)),INTENT(in)::HDDday  ! HDDday=HDD(id-1,4) HDD at the begining of today (id-1)
    REAL(KIND(1d0)),INTENT(in)::qf
    REAL(KIND(1d0)),INTENT(in)::qn1_bup

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::alb  !< albedo [-]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::emis !< emissivity [-]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::cp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::kk   !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::ch   !< bulk transfer coef [J m-3 K-1]

    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::MetForcingData_grid !< met forcing array of grid

    REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout) ::qn1_store
    REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout) ::qn1_S_store !< stored qn1 [W m-2]

    REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_av_store
    REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_S_av_store !< average net radiation over previous hour [W m-2]
    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(inout)::surf

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::deltaQi ! storage heat flux of snow surfaces
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowFrac

    REAL(KIND(1d0)),INTENT(out)::qn1_S
    REAL(KIND(1d0)),INTENT(out):: qs ! storage heat flux
    REAL(KIND(1d0)),INTENT(out):: a1 !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out):: a2 !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out):: a3 !< AnOHM coefficients of grid [W m-2]



    IF(StorageHeatMethod==1) THEN           !Use OHM to calculate QS
       IF(OHMIncQF == 1) THEN      !Calculate QS using QSTAR+QF
          IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
          CALL OHM(qf+qn1_bup,qn1_store,qn1_av_store,&
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
          CALL AnOHM(qf+qn1_bup,qn1_store,qn1_av_store,&
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
  !=======================================================================

  !==========================water balance================================
  SUBROUTINE SUEWS_cal_Water(&
       Diagnose,&!input
       snowUse,NonWaterFraction,addPipes,addImpervious,addVeg,addWaterBody,&
       state,soilmoist,sfr,surf,WaterDist,nsh_real,&
       drain_per_tstep,&  !output
       drain,AddWaterRunoff,&
       AdditionalWater,runoffPipes,runoff_per_interval,&
       addWater,stateOld,soilmoistOld)

    IMPLICIT NONE
    INTEGER,PARAMETER :: nsurf=7! number of surface types
    INTEGER,PARAMETER::WaterSurf = 7
    INTEGER,INTENT(in) ::Diagnose
    INTEGER,INTENT(in) ::snowUse

    REAL(KIND(1d0)),INTENT(in)::NonWaterFraction
    REAL(KIND(1d0)),INTENT(in)::addPipes
    REAL(KIND(1d0)),INTENT(in)::addImpervious
    REAL(KIND(1d0)),INTENT(in)::addVeg
    REAL(KIND(1d0)),INTENT(in)::addWaterBody
    REAL(KIND(1d0)),INTENT(in)::nsh_real !nsh cast as a real for use in calculations

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)          ::state
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)          ::soilmoist
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)          ::sfr
    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in)        ::surf
    REAL(KIND(1d0)),DIMENSION(nsurf+1,nsurf-1),INTENT(in)::WaterDist

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: drain         !Drainage of surface type "is" [mm]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: AddWaterRunoff!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: addWater      !water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: stateOld
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: soilmoistOld

    REAL(KIND(1d0)),INTENT(out)::drain_per_tstep
    REAL(KIND(1d0)),INTENT(out)::AdditionalWater
    REAL(KIND(1d0)),INTENT(out)::runoffPipes
    REAL(KIND(1d0)),INTENT(out)::runoff_per_interval
    INTEGER:: is

    ! Retain previous surface state and soil moisture state
    stateOld     = state     !State of each surface [mm] for the previous timestep
    soilmoistOld = soilmoist !Soil moisture of each surface [mm] for the previous timestep


    !============= Grid-to-grid runoff =============
    ! Calculate additional water coming from other grids
    ! i.e. the variables addImpervious, addVeg, addWaterBody, addPipes
    !call RunoffFromGrid(GridFromFrac)  !!Need to code between-grid water transfer

    ! Sum water coming from other grids (these are expressed as depths over the whole surface)
    AdditionalWater = addPipes+addImpervious+addVeg+addWaterBody  ![mm]

    ! Initialise runoff in pipes
    runoffPipes         = addPipes !Water flowing in pipes from other grids. No need for scaling??
    !! CHECK p_i
    runoff_per_interval = addPipes !pipe plor added to total runoff.


    !================== Drainage ===================
    ! Calculate drainage for each soil subsurface (excluding water body)
    IF(Diagnose==1) WRITE(*,*) 'Calling Drainage...'

    IF (NonWaterFraction/=0) THEN !Soil states only calculated if soil exists. LJ June 2017
       DO is=1,nsurf-1

          CALL drainage(&
               is,&! input:
               state(is),&
               surf(6,is),&
               surf(2,is),&
               surf(3,is),&
               surf(4,is),&
               nsh_real,&
               drain(is))! output

          ! !HCW added and changed to surf(6,is) here 20 Feb 2015
          ! drain_per_tstep=drain_per_tstep+(drain(is)*sfr(is)/NonWaterFraction)   !No water body included
       ENDDO
       drain_per_tstep=DOT_PRODUCT(drain(1:nsurf-1),sfr(1:nsurf-1))/NonWaterFraction !No water body included
    ELSE
       drain(1:nsurf-1)=0
    ENDIF

    drain(WaterSurf) = 0  ! Set drainage from water body to zero

    ! Distribute water within grid, according to WithinGridWaterDist matrix (Cols 1-7)
    IF(Diagnose==1) WRITE(*,*) 'Calling ReDistributeWater...'
    ! CALL ReDistributeWater
    !Calculates AddWater(is)
    CALL ReDistributeWater(&
                                ! input:
         nsurf,& ! surface type number
         WaterSurf,&
         snowUse,&
         WaterDist,  &
         sfr,   &!
         Drain,&
                                ! output:
         AddWaterRunoff,&
         addWater)

  END SUBROUTINE SUEWS_cal_Water
  !=======================================================================

  !===============initialize sensible heat flux============================
  SUBROUTINE SUEWS_init_QH(&
       qh_obs,avdens,avcp,h_mod,qn1,dectime,&!input
       H_init)!output

    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(in)::qh_obs
    REAL(KIND(1d0)),INTENT(in)::avdens
    REAL(KIND(1d0)),INTENT(in)::avcp
    REAL(KIND(1d0)),INTENT(in)::h_mod
    REAL(KIND(1d0)),INTENT(in)::qn1
    REAL(KIND(1d0)),INTENT(in)::dectime
    REAL(KIND(1d0)),INTENT(out)::H_init


    REAL(KIND(1d0)),PARAMETER::NAN=-999
    INTEGER,PARAMETER::notUsedI=-999

    ! Calculate kinematic heat flux (w'T') from sensible heat flux [W m-2] from observed data (if available) or LUMPS
    IF(qh_obs/=NAN) THEN   !if(qh_obs/=NAN) qh=qh_obs   !Commented out by HCW 04 Mar 2015
       H_init=qh_obs/(avdens*avcp)  !Use observed value
    ELSE
       IF(h_mod/=NAN) THEN
          H_init = h_mod/(avdens*avcp)   !Use LUMPS value
       ELSE
          H_init=(qn1*0.2)/(avdens*avcp)   !If LUMPS has had a problem, we still need a value
          CALL ErrorHint(38,'LUMPS unable to calculate realistic value for H_mod.',h_mod, dectime, notUsedI)
       ENDIF
    ENDIF

  END SUBROUTINE SUEWS_init_QH
  !========================================================================

  !================latent heat flux and surface wetness===================
  ! TODO: optimise the structure of this function
  SUBROUTINE SUEWS_cal_QE(&
       Diagnose,&!input
       id,tstep,imin,it,ity,snowfractionchoice,snowCalcSwitch,DayofWeek,CRWmin,CRWmax,&
       nsh_real,lvS_J_kg,lv_j_kg,avdens,waterdens,avRh,Press_hPa,Temp_C,&
       RAsnow,psyc_hPa,avcp,sIce_hPa,&
       PervFraction,vegfraction,addimpervious,qn1_SF,qf,qs,vpd_hPa,s_hPa,&
       ResistSurf,ra,rb,tstep_real,snowdensmin,precip,PipeCapacity,RunoffToWater,&
       NonWaterFraction,wu_EveTr,wu_DecTr,wu_Grass,addVeg,addWaterBody,SnowLimPaved,SnowLimBuild,&
       SurfaceArea,drain,WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
       freezmelt,freezstate,freezstatevol,Qm_Melt,Qm_rain,Tsurf_ind,sfr,StateLimit,surf,&
       runoff_per_interval,& ! inout:
       state,soilmoist,SnowPack,snowFrac,MeltWaterStore,&
       SnowDepth,iceFrac,addwater,addwaterrunoff,SnowDens,SurplusEvap,&
       snowProf,& ! output:
       runoffSnow,runoff,runoffSoil,chang,changSnow,SnowToSurf,snowD,ev_snow,SnowRemoval,&
       evap,rss_nsurf,p_mm,rss,qe,state_per_tstep,NWstate_per_tstep,qeOut,&
       swe,ev,chSnow_per_interval,ev_per_tstep,qe_per_tstep,runoff_per_tstep,&
       surf_chang_per_tstep,runoffPipes,mwstore,runoffwaterbody,FlowChange,&
       runoffAGimpervious_m3,runoffAGveg_m3,runoffWaterBody_m3,runoffPipes_m3)

    IMPLICIT NONE

    INTEGER,PARAMETER::nsurf=7! number of surface types
    ! INTEGER,PARAMETER::PavSurf   = 1  !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER::BldgSurf  = 2  !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER::ConifSurf = 3  !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER::DecidSurf = 4  !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER::GrassSurf = 5
    ! INTEGER,PARAMETER::BSoilSurf = 6!New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER::WaterSurf = 7

    INTEGER,INTENT(in) ::Diagnose
    INTEGER,INTENT(in) ::id
    INTEGER,INTENT(in) ::tstep
    INTEGER,INTENT(in) ::imin
    INTEGER,INTENT(in) ::it
    INTEGER,INTENT(in) ::ity !Evaporation calculated according to Rutter (1) or Shuttleworth (2)
    INTEGER,INTENT(in) ::snowfractionchoice

    INTEGER,DIMENSION(nsurf),INTENT(in)::snowCalcSwitch
    INTEGER,DIMENSION(366,2),INTENT(in)::DayofWeek

    REAL(KIND(1d0)),INTENT(in)::CRWmin
    REAL(KIND(1d0)),INTENT(in)::CRWmax
    REAL(KIND(1d0)),INTENT(in)::nsh_real
    REAL(KIND(1d0)),INTENT(in)::lvS_J_kg
    REAL(KIND(1d0)),INTENT(in)::lv_j_kg
    REAL(KIND(1d0)),INTENT(in)::avdens
    REAL(KIND(1d0)),INTENT(in)::waterdens
    REAL(KIND(1d0)),INTENT(in)::avRh
    REAL(KIND(1d0)),INTENT(in)::Press_hPa
    REAL(KIND(1d0)),INTENT(in)::Temp_C
    REAL(KIND(1d0)),INTENT(in)::RAsnow
    REAL(KIND(1d0)),INTENT(in)::psyc_hPa
    REAL(KIND(1d0)),INTENT(in)::avcp
    REAL(KIND(1d0)),INTENT(in)::sIce_hPa
    REAL(KIND(1d0)),INTENT(in)::PervFraction
    REAL(KIND(1d0)),INTENT(in)::vegfraction
    REAL(KIND(1d0)),INTENT(in)::addimpervious
    REAL(KIND(1d0)),INTENT(in)::qn1_SF
    REAL(KIND(1d0)),INTENT(in)::qf
    REAL(KIND(1d0)),INTENT(in)::qs
    REAL(KIND(1d0)),INTENT(in)::vpd_hPa
    REAL(KIND(1d0)),INTENT(in)::s_hPa
    REAL(KIND(1d0)),INTENT(in)::ResistSurf
    REAL(KIND(1d0)),INTENT(in)::ra
    REAL(KIND(1d0)),INTENT(in)::rb
    REAL(KIND(1d0)),INTENT(in)::tstep_real
    REAL(KIND(1d0)),INTENT(in)::snowdensmin
    REAL(KIND(1d0)),INTENT(in)::precip
    REAL(KIND(1d0)),INTENT(in)::PipeCapacity
    REAL(KIND(1d0)),INTENT(in)::RunoffToWater
    REAL(KIND(1d0)),INTENT(in)::NonWaterFraction
    REAL(KIND(1d0)),INTENT(in)::wu_EveTr!Water use for evergreen trees/shrubs [mm]
    REAL(KIND(1d0)),INTENT(in)::wu_DecTr!Water use for deciduous trees/shrubs [mm]
    REAL(KIND(1d0)),INTENT(in)::wu_Grass!Water use for grass [mm]
    REAL(KIND(1d0)),INTENT(in)::addVeg!Water from vegetated surfaces of other grids [mm] for whole surface area
    REAL(KIND(1d0)),INTENT(in)::addWaterBody!Water from water surface of other grids [mm] for whole surface area
    REAL(KIND(1d0)),INTENT(in)::SnowLimPaved
    REAL(KIND(1d0)),INTENT(in)::SnowLimBuild
    REAL(KIND(1d0)),INTENT(in)::SurfaceArea

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::drain
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::WetThresh
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::stateOld
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::mw_ind
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::soilstorecap
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::rainonsnow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::freezmelt
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::freezstate
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::freezstatevol
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::Qm_Melt
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::Qm_rain
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::Tsurf_ind
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::StateLimit !Limit for state of each surface type [mm] (specified in input files)

    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in)::surf

    !Updated status: input and output
    REAL(KIND(1d0)),INTENT(inout)::runoff_per_interval! Total water transported to each grid for grid-to-grid connectivity

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::state
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::soilmoist
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowPack
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::snowFrac
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::MeltWaterStore
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowDepth
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::iceFrac
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::addwater
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::addwaterrunoff
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowDens
    REAL(KIND(1d0)),DIMENSION(2),INTENT(inout)    ::SurplusEvap        !Surplus for evaporation in 5 min timestep

    ! output:
    REAL(KIND(1d0)), DIMENSION(0:23,2),INTENT(out):: snowProf

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoffSnow !Initialize for runoff caused by snowmelting
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoff
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoffSoil
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::chang
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::changSnow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::SnowToSurf
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowD
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::ev_snow
    REAL(KIND(1d0)),DIMENSION(2),INTENT(out)::SnowRemoval
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::evap
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::rss_nsurf

    REAL(KIND(1d0)),INTENT(out)::p_mm!Inputs to surface water balance
    REAL(KIND(1d0)),INTENT(out)::rss
    REAL(KIND(1d0)),INTENT(out)::qe ! latent heat flux [W m-2]
    REAL(KIND(1d0)),INTENT(out)::state_per_tstep
    REAL(KIND(1d0)),INTENT(out)::NWstate_per_tstep
    REAL(KIND(1d0)),INTENT(out)::qeOut
    REAL(KIND(1d0)),INTENT(out)::swe
    REAL(KIND(1d0)),INTENT(out)::ev
    REAL(KIND(1d0)),INTENT(out)::chSnow_per_interval
    REAL(KIND(1d0)),INTENT(out)::ev_per_tstep
    REAL(KIND(1d0)),INTENT(out)::qe_per_tstep
    REAL(KIND(1d0)),INTENT(out)::runoff_per_tstep
    REAL(KIND(1d0)),INTENT(out)::surf_chang_per_tstep
    REAL(KIND(1d0)),INTENT(out)::runoffPipes
    REAL(KIND(1d0)),INTENT(out)::mwstore
    REAL(KIND(1d0)),INTENT(out)::runoffwaterbody
    REAL(KIND(1d0)),INTENT(out)::FlowChange
    REAL(KIND(1d0)),INTENT(out)::runoffAGimpervious_m3
    REAL(KIND(1d0)),INTENT(out)::runoffAGveg_m3
    REAL(KIND(1d0)),INTENT(out)::runoffWaterBody_m3
    REAL(KIND(1d0)),INTENT(out)::runoffPipes_m3


    ! local:
    INTEGER:: is
    REAL(KIND(1d0))::runoffAGveg
    REAL(KIND(1d0))::runoffAGimpervious
    REAL(KIND(1d0))::surplusWaterBody
    REAL(KIND(1d0))::pin!Rain per time interval
    REAL(KIND(1d0))::sae
    REAL(KIND(1d0))::vdrc
    REAL(KIND(1d0))::sp
    REAL(KIND(1d0))::numPM
    REAL(KIND(1d0))::tlv


    tlv=lv_J_kg/tstep_real !Latent heat of vapourisation per timestep

    pin=MAX(0.,Precip)!Initiate rain data [mm]
    ! IF(Precip>0) THEN   !Initiate rain data [mm]
    !    pin=Precip
    ! ELSE
    !    pin=0
    ! ENDIF

    ! Initialize the output variables
    qe=0
    ev=0
    ev_snow=0
    ev_per_tstep=0
    surf_chang_per_tstep=0
    runoff_per_tstep=0
    state_per_tstep=0
    NWstate_per_tstep=0
    qeOut=0
    runoffwaterbody=0
    chSnow_per_interval=0

    runoffAGveg        = 0
    runoffAGimpervious = 0
    surplusWaterBody   = 0
    runoffSoil         = 0
    runoff             = 0
    chang              = 0
    SurplusEvap        = 0

    !========= these need to be wrapped================================
    sae   = s_hPa*(qn1_SF+qf-qs)    !s_haPa - slope of svp vs t curve. qn1 changed to qn1_SF, lj in May 2013
    vdrc  = vpd_hPa*avdens*avcp
    sp    = s_hPa/psyc_hPa
    numPM = sae+vdrc/ra
    !write(*,*) numPM, sae, vdrc/ra, s_hPA+psyc_hPa, NumPM/(s_hPA+psyc_hPa)
    !========= these need to be wrapped end================================

    IF(Diagnose==1) WRITE(*,*) 'Calling evap_SUEWS and SoilStore...'
    DO is=1,nsurf   !For each surface in turn
       IF (snowCalcSwitch(is)==1) THEN
          IF (sfr(is)/=0) THEN
             IF(Diagnose==1) WRITE(*,*) 'Calling SnowCalc...'
             CALL SnowCalc(&
                  id,& !input
                  tstep,imin,it,is,&
                  snowfractionchoice,ity,CRWmin,CRWmax,nsh_real,lvS_J_kg,lv_j_kg,avdens,waterdens,&
                  avRh,Press_hPa,Temp_C,RAsnow,psyc_hPa,avcp,sIce_hPa,&
                  PervFraction,vegfraction,addimpervious,&
                  numPM,s_hPa,ResistSurf,sp,ra,rb,tlv,snowdensmin,precip,&
                  PipeCapacity,RunoffToWater,runoffAGimpervious,runoffAGveg,&
                  addVeg,surplusWaterBody,SnowLimPaved,SnowLimBuild,drain,&
                  WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
                  freezmelt,freezstate,freezstatevol,&
                  Qm_Melt,Qm_rain,Tsurf_ind,sfr,DayofWeek,surf,&
                  SnowPack,&!inout
                  snowFrac,MeltWaterStore,SnowDepth,iceFrac,addwater,addwaterrunoff,SnowDens,&
                  runoffSnow,& ! output
                  runoff,runoffSoil,chang,changSnow,SnowToSurf,state,snowD,ev_snow,soilmoist,&
                  SnowRemoval,snowProf,swe,ev,chSnow_per_interval,ev_per_tstep,qe_per_tstep,&
                  runoff_per_tstep,surf_chang_per_tstep,runoffPipes,mwstore,runoffwaterbody,&
                  FlowChange)
          ELSE
             snowFrac(is) = 0
             SnowDens(is) = 0
             SnowPack(is) = 0
          ENDIF
       ELSE

          !Calculates ev [mm]
          CALL Evap_SUEWS(&
               ity,&! input: !Evaporation calculated according to Rutter (1) or Shuttleworth (2)
               state(is),& ! wetness status
               WetThresh(is),&!When State > WetThresh, rs=0 limit in SUEWS_evap [mm] (specified in input files)
               surf(6,is),& ! = surf(6,is), current storage capacity [mm]
               numPM,&!numerator of P-M eqn
               s_hPa,&!Vapour pressure versus temperature slope in hPa
               psyc_hPa,&!Psychometric constant in hPa
               ResistSurf,&!Surface resistance
               sp,&!Term in calculation of E
               ra,&!Aerodynamic resistance
               rb,&!Boundary layer resistance
               tlv,&!Latent heat of vaporization per timestep [J kg-1 s-1], (tlv=lv_J_kg/tstep_real)
               rss,&! output:
               ev,&
               qe) ! latent heat flux [W m-2]


          rss_nsurf(is) = rss !Store rss for each surface
          ! CALL soilstore    !Surface water balance and soil store updates (can modify ev, updates state)
          !Surface water balance and soil store updates (can modify ev, updates state)
          CALL soilstore(&
               is,& ! input: ! surface type
               sfr,&! surface fractions
               PipeCapacity,&!Capacity of pipes to transfer water
               RunoffToWater,&!Fraction of surface runoff going to water body
               pin,&!Rain per time interval
               wu_EveTr,&!Water use for evergreen trees/shrubs [mm]
               wu_DecTr,&!Water use for deciduous trees/shrubs [mm]
               wu_Grass,&!Water use for grass [mm]
               AddWater,&!Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
               addImpervious,&!Water from impervious surfaces of other grids [mm] for whole surface area
               nsh_real,&!nsh cast as a real for use in calculations
               stateOld,&!Wetness status of each surface type from previous timestep [mm]
               AddWaterRunoff,&!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
               PervFraction,&! sum of surface cover fractions for impervious surfaces
               addVeg,&!Water from vegetated surfaces of other grids [mm] for whole surface area
               soilstoreCap,&!Capacity of soil store for each surface [mm]
               addWaterBody,&!Water from water surface of other grids [mm] for whole surface area
               FlowChange,&!Difference between the input and output flow in the water body
               StateLimit,&!Limit for state of each surface type [mm] (specified in input files)
               runoffAGimpervious,&!  inout:!Above ground runoff from impervious surface [mm] for whole surface area
               surplusWaterBody,&!Extra runoff that goes to water body [mm] as specified by RunoffToWater
               runoffAGveg,&!Above ground runoff from vegetated surfaces [mm] for whole surface area
               runoffPipes,&!Runoff in pipes [mm] for whole surface area
               ev,&!Evaporation
               soilmoist,&!Soil moisture of each surface type [mm]
               SurplusEvap,&!Surplus for evaporation in 5 min timestep
               runoffWaterBody,&!Above ground runoff from water surface [mm] for whole surface area
               runoff_per_interval,&! Total water transported to each grid for grid-to-grid connectivity
               p_mm,&!output: !Inputs to surface water balance
               chang,&!Change in state [mm]
               runoff,&!Runoff from each surface type [mm]
               drain,&!Drainage of each surface type [mm]
               state&!Wetness status of each surface type [mm]
               )

          evap(is)     = ev !Store ev for each surface

          ! Sum evaporation from different surfaces to find total evaporation [mm]
          ev_per_tstep = ev_per_tstep+evap(is)*sfr(is)

          ! Sum change from different surfaces to find total change to surface state
          surf_chang_per_tstep = surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)
          ! Sum runoff from different surfaces to find total runoff
          runoff_per_tstep     = runoff_per_tstep+runoff(is)*sfr(is)
          ! Calculate total state (including water body)
          state_per_tstep      = state_per_tstep+(state(is)*sfr(is))
          ! Calculate total state (excluding water body)

          IF (NonWaterFraction/=0 .AND. is.NE.WaterSurf) THEN
             NWstate_per_tstep=NWstate_per_tstep+(state(is)*sfr(is)/NonWaterFraction)
          ENDIF

          ChangSnow(is)  = 0
          runoffSnow(is) = 0

       ENDIF
    ENDDO  !end loop over surfaces


    ! Convert evaporation to latent heat flux [W m-2]
    qe_per_tstep = ev_per_tstep*tlv
    qeOut        = qe_per_tstep

    ! Calculate volume of water that will move between grids
    ! Volume [m3] = Depth relative to whole area [mm] / 1000 [mm m-1] * SurfaceArea [m2]
    ! Need to use these volumes when converting back to addImpervious, AddVeg and AddWater
    runoffAGimpervious_m3 = runoffAGimpervious/1000 *SurfaceArea
    runoffAGveg_m3        = runoffAGveg/1000 *SurfaceArea
    runoffWaterBody_m3    = runoffWaterBody/1000 *SurfaceArea
    runoffPipes_m3        = runoffPipes/1000 *SurfaceArea

  END SUBROUTINE SUEWS_cal_QE
  !========================================================================

  !===============sensible heat flux======================================
  SUBROUTINE SUEWS_cal_QH(&
       opt_QH,&
       qn1,&
       qf,&
       QmRain,&
       qeOut,&
       qs,&
       QmFreez,&
       qm,&
       avdens,&
       avcp,&
       tsurf,&
       Temp_C,&
       ra,&
       qh)
    IMPLICIT NONE

    INTEGER,INTENT(in) :: opt_QH ! option for QH calculation: 1, residual; 2, resistance-based

    REAL(KIND(1d0)),INTENT(in)::qn1
    REAL(KIND(1d0)),INTENT(in)::qf
    REAL(KIND(1d0)),INTENT(in)::QmRain
    REAL(KIND(1d0)),INTENT(in)::qeOut
    REAL(KIND(1d0)),INTENT(in)::qs
    REAL(KIND(1d0)),INTENT(in)::QmFreez
    REAL(KIND(1d0)),INTENT(in)::qm
    REAL(KIND(1d0)),INTENT(in)::avdens
    REAL(KIND(1d0)),INTENT(in)::avcp
    REAL(KIND(1d0)),INTENT(in)::tsurf
    REAL(KIND(1d0)),INTENT(in)::Temp_C
    REAL(KIND(1d0)),INTENT(in)::ra


    REAL(KIND(1d0)),INTENT(out)::qh

    REAL(KIND(1d0)),PARAMETER::NAN=-999

    ! ! Calculate QH using resistance method (for testing HCW 06 Jul 2016)

    SELECT CASE (opt_QH)
    CASE (1)
       ! Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
       qh=(qn1+qf+QmRain)-(qeOut+qs+Qm+QmFreez)     !qh=(qn1+qf+QmRain+QmFreez)-(qeOut+qs+Qm)

    CASE (2)
       IF(ra/=0) THEN
          qh = avdens*avcp*(tsurf-Temp_C)/ra
       ELSE
          qh=NAN
       ENDIF
    END SELECT

  END SUBROUTINE SUEWS_cal_QH
  !========================================================================

  !===============Resistance Calculations=======================
  SUBROUTINE SUEWS_cal_Resistance(&
       StabilityMethod,&!input:
       Diagnose,AerodynamicResistanceMethod,RoughLenHeatMethod,snowUse,&
       id,it,gsModel,SMDMethod,&
       qh_obs,avdens,avcp,h_mod,qn1_bup,dectime,zzd,z0M,zdm,&
       avU1,Temp_C,L_mod,UStar,VegFraction,&
       avkdn,Kmax,G1,G2,G3,G4,G5,G6,S1,S2,TH,TL,dq,&
       xsmd,vsmd,MaxConductance,LAIMax,LAI_id,snowFrac,sfr,&
       Tstar,&!output
       psim,gsc,ResistSurf,RA,RAsnow,rb)

    IMPLICIT NONE
    INTEGER,PARAMETER :: nsurf=7! number of surface types


    INTEGER,PARAMETER::ConifSurf = 3  !New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER::DecidSurf = 4  !New surface classes: Grass = 5th/7 surfaces
    INTEGER,PARAMETER::GrassSurf = 5

    INTEGER,PARAMETER::WaterSurf = 7



    INTEGER,INTENT(in)::StabilityMethod
    INTEGER,INTENT(in)::Diagnose
    INTEGER,INTENT(in)::AerodynamicResistanceMethod
    INTEGER,INTENT(in)::RoughLenHeatMethod
    INTEGER,INTENT(in)::snowUse
    INTEGER,INTENT(in)::id
    INTEGER,INTENT(in)::it       !time: day of year and hour
    INTEGER,INTENT(in)::gsModel  !Choice of gs parameterisation (1 = Ja11, 2 = Wa16)
    INTEGER,INTENT(in)::SMDMethod!Method of measured soil moisture

    REAL(KIND(1d0)),INTENT(in)::qh_obs
    REAL(KIND(1d0)),INTENT(in)::avdens
    REAL(KIND(1d0)),INTENT(in)::avcp
    REAL(KIND(1d0)),INTENT(in)::h_mod
    REAL(KIND(1d0)),INTENT(in)::qn1_bup
    REAL(KIND(1d0)),INTENT(in)::dectime    !Decimal time
    REAL(KIND(1d0)),INTENT(in)::zzd        !Active measurement height (meas. height-displac. height)
    REAL(KIND(1d0)),INTENT(in)::z0M        !Aerodynamic roughness length
    REAL(KIND(1d0)),INTENT(in)::zdm        !Displacement height
    REAL(KIND(1d0)),INTENT(in)::avU1       !Average wind speed
    REAL(KIND(1d0)),INTENT(in)::Temp_C     !Air temperature
    REAL(KIND(1d0)),INTENT(in)::VegFraction!Fraction of vegetation
    REAL(KIND(1d0)),INTENT(in)::avkdn      !Average downwelling shortwave radiation
    REAL(KIND(1d0)),INTENT(in)::Kmax       !Annual maximum hourly solar radiation
    REAL(KIND(1d0)),INTENT(in)::G1         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::G2         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::G3         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::G4         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::G5         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::G6         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::S1         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::S2         !Fitted parameters related to surface res. calculations
    REAL(KIND(1d0)),INTENT(in)::TH         !Maximum temperature limit
    REAL(KIND(1d0)),INTENT(in)::TL         !Minimum temperature limit
    REAL(KIND(1d0)),INTENT(in)::dq         !Specific humidity deficit
    REAL(KIND(1d0)),INTENT(in)::xsmd       !Measured soil moisture deficit
    REAL(KIND(1d0)),INTENT(in)::vsmd       !Soil moisture deficit for vegetated surfaces only (what about BSoil?)

    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::MaxConductance!Max conductance [mm s-1]
    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::LAIMax        !Max LAI [m2 m-2]
    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::LAI_id        !=LAI_id(id-1,:), LAI for each veg surface [m2 m-2]

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::snowFrac      !Surface fraction of snow cover
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr           !Surface fractions [-]

    REAL(KIND(1d0)),INTENT(out)::Tstar     !T*
    REAL(KIND(1d0)),INTENT(out)::UStar     !Friction velocity
    REAL(KIND(1d0)),INTENT(out)::psim      !Stability function of momentum
    REAL(KIND(1d0)),INTENT(out)::gsc       !Surface Layer Conductance
    REAL(KIND(1d0)),INTENT(out)::ResistSurf!Surface resistance
    REAL(KIND(1d0)),INTENT(out)::RA        !Aerodynamic resistance [s m^-1]
    REAL(KIND(1d0)),INTENT(out)::RAsnow    !Aerodynamic resistance for snow [s m^-1]
    REAL(KIND(1d0)),INTENT(out)::rb        !boundary layer resistance shuttleworth

    REAL(KIND(1d0))::H_init !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
    REAL(KIND(1d0))::L_mod  !Obukhov length

    ! Get first estimate of sensible heat flux. Modified by HCW 26 Feb 2015
    CALL SUEWS_init_QH(&
         qh_obs,avdens,avcp,h_mod,qn1_bup,dectime,&
         H_init)

    IF(Diagnose==1) WRITE(*,*) 'Calling STAB_lumps...'
    !u* and Obukhov length out
    CALL STAB_lumps(&
         StabilityMethod,&  ! input
         dectime,& !Decimal time
         zzd,&     !Active measurement height (meas. height-displac. height)
         z0M,&     !Aerodynamic roughness length
         zdm,&     !Displacement height
         avU1,&    !Average wind speed
         Temp_C,&  !Air temperature
         H_init,& !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
         L_mod,&! output: !Obukhov length
         Tstar,& !T*
         UStar,& !Friction velocity
         psim)!Stability function of momentum

    IF(Diagnose==1) WRITE(*,*) 'Calling AerodynamicResistance...'
    CALL AerodynamicResistance(&
         ZZD,&! input:
         z0m,&
         AVU1,&
         L_mod,&
         UStar,&
         VegFraction,&
         AerodynamicResistanceMethod,&
         StabilityMethod,&
         RoughLenHeatMethod,&
         RA) ! output:

    IF (snowUse==1) THEN
       IF(Diagnose==1) WRITE(*,*) 'Calling AerodynamicResistance for snow...'
       CALL AerodynamicResistance(&
            ZZD,&! input:
            z0m,&
            AVU1,&
            L_mod,&
            UStar,&
            VegFraction,&
            AerodynamicResistanceMethod,&
            StabilityMethod,&
            3,&
            RAsnow)     ! output:
    ENDIF

    IF(Diagnose==1) WRITE(*,*) 'Calling SurfaceResistance...'
    ! CALL SurfaceResistance(id,it)   !qsc and surface resistance out
    CALL  SurfaceResistance(&
         id,it,&! input:
         SMDMethod,&
         ConifSurf,&
         DecidSurf,&
         GrassSurf,&
         WaterSurf,&
         snowFrac,&
         sfr,&
         nsurf,&
         avkdn,&
         Temp_C,&
         dq,&
         xsmd,&
         vsmd,&
         MaxConductance,&
         LAIMax,&
         LAI_id,&
         gsModel,&
         Kmax,&
         G1,&
         G2,&
         G3,&
         G4,&
         G5,&
         G6,&
         TH,&
         TL,&
         S1,&
         S2,&
         gsc,&! output:
         ResistSurf)

    IF(Diagnose==1) WRITE(*,*) 'Calling BoundaryLayerResistance...'
    CALL BoundaryLayerResistance(&
         zzd,&! input:     !Active measurement height (meas. height-displac. height)
         z0M,&     !Aerodynamic roughness length
         avU1,&    !Average wind speed
         UStar,&  ! input/output:
         rb)  ! output:

  END SUBROUTINE SUEWS_cal_Resistance
  !========================================================================

  !==============Update output arrays=========================
  SUBROUTINE SUEWS_update_output(&
       ReadLinesMetdata,ncolumnsDataOut,ncolumnsDataOutSnow,NumberOfGrids,&
       Gridiv,&
       iy,&
       iy_prev_t,&
       id,&
       id_prev_t,&
       it,&
       imin,&
       SNOWuse,&
       ir,&
       AdditionalWater,&
       alb,&
       avkdn,&
       avU10_ms,&
       azimuth,&
       chSnow_per_interval,&
       dectime,&
       drain_per_tstep,&
       E_mod,&
       ev_per_tstep,&
       ext_wu,&
       Fc,&
       Fc_build,&
       Fc_metab,&
       Fc_photo,&
       Fc_respi,&
       Fc_traff,&
       fcld,&
       FlowChange,&
       freezMelt,&
       h_mod,&
       int_wu,&
       kup,&
       kup_ind_snow,&
       l_mod,&
       LAI,&
       ldown,&
       lup,&
       MeltWaterStore,&
       mw_ind,&
       mwh,&
       MwStore,&
       nsh_real,&
       NWstate_per_tstep,&
       Precip,&
       q2_gkg,&
       qeOut,&
       qf,&
       qh,&
       QH_r,&
       Qm,&
       Qm_freezState,&
       Qm_melt,&
       Qm_rain,&
       QmFreez,&
       QmRain,&
       qn1_bup,&
       qn1_ind_snow,&
       qn1_S,&
       qn1_SF,&
       qs,&
       RA,&
       rainOnSnow,&
       resistsurf,&
       runoff_per_tstep,&
       runoffAGimpervious,&
       runoffAGveg,&
       runoffPipes,&
       runoffSoil_per_tstep,&
       runoffWaterBody,&
       sfr,&
       smd,&
       smd_nsurf,&
       SnowAlb,&
       SnowDens,&
       snowDepth,&
       SnowRemoval,&
       SoilState,&
       state,&
       state_per_tstep,&
       surf_chang_per_tstep,&
       swe,&
       t2_C,&
       tot_chang_per_tstep,&
       tsurf,&
       Tsurf_ind_snow,&
       UStar,&
       wu_DecTr,&
       wu_EveTr,&
       wu_Grass,&
       z0m,&
       zdm,&
       zenith_deg,&
       SnowFrac,&
       SnowPack,&
       dataOut,dataOutSnow)
    IMPLICIT NONE

    INTEGER,PARAMETER::ndays    = 366
    INTEGER,PARAMETER::nvegsurf = 3
    INTEGER,PARAMETER::nsurf    = 7
    REAL(KIND(1d0)),PARAMETER :: NAN=-999

    INTEGER,INTENT(in) ::ReadLinesMetdata
    INTEGER,INTENT(in) ::ncolumnsDataOut
    INTEGER,INTENT(in) ::ncolumnsDataOutSnow
    INTEGER,INTENT(in) ::NumberOfGrids
    INTEGER,INTENT(in) :: Gridiv
    INTEGER,INTENT(in) :: iy
    INTEGER,INTENT(in) :: iy_prev_t
    INTEGER,INTENT(in) :: id
    INTEGER,INTENT(in) :: id_prev_t
    INTEGER,INTENT(in) :: it
    INTEGER,INTENT(in) :: imin
    INTEGER,INTENT(in) :: SNOWuse
    INTEGER,INTENT(in) :: ir

    REAL(KIND(1d0)),INTENT(in) :: AdditionalWater
    REAL(KIND(1d0)),INTENT(in) :: alb(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: avkdn
    REAL(KIND(1d0)),INTENT(in) :: avU10_ms
    REAL(KIND(1d0)),INTENT(in) :: azimuth
    REAL(KIND(1d0)),INTENT(in) :: chSnow_per_interval
    REAL(KIND(1d0)),INTENT(in) :: dectime
    REAL(KIND(1d0)),INTENT(in) :: drain_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: E_mod
    REAL(KIND(1d0)),INTENT(in) :: ev_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: ext_wu
    REAL(KIND(1d0)),INTENT(in) :: Fc
    REAL(KIND(1d0)),INTENT(in) :: Fc_build
    REAL(KIND(1d0)),INTENT(in) :: Fc_metab
    REAL(KIND(1d0)),INTENT(in) :: Fc_photo
    REAL(KIND(1d0)),INTENT(in) :: Fc_respi
    REAL(KIND(1d0)),INTENT(in) :: Fc_traff
    REAL(KIND(1d0)),INTENT(in) :: fcld
    REAL(KIND(1d0)),INTENT(in) :: FlowChange
    REAL(KIND(1d0)),INTENT(in) :: freezMelt(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: h_mod
    REAL(KIND(1d0)),INTENT(in) :: int_wu
    REAL(KIND(1d0)),INTENT(in) :: kup
    REAL(KIND(1d0)),INTENT(in) :: kup_ind_snow(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: l_mod
    REAL(KIND(1d0)),INTENT(in) :: LAI(-4:ndays, nvegsurf)
    REAL(KIND(1d0)),INTENT(in) :: ldown
    REAL(KIND(1d0)),INTENT(in) :: lup
    REAL(KIND(1d0)),INTENT(in) :: MeltWaterStore(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: mw_ind(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: mwh
    REAL(KIND(1d0)),INTENT(in) :: MwStore
    REAL(KIND(1d0)),INTENT(in) :: nsh_real
    REAL(KIND(1d0)),INTENT(in) :: NWstate_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: Precip
    REAL(KIND(1d0)),INTENT(in) :: q2_gkg
    REAL(KIND(1d0)),INTENT(in) :: qeOut
    REAL(KIND(1d0)),INTENT(in) :: qf
    REAL(KIND(1d0)),INTENT(in) :: qh
    REAL(KIND(1d0)),INTENT(in) :: QH_r
    REAL(KIND(1d0)),INTENT(in) :: Qm
    REAL(KIND(1d0)),INTENT(in) :: Qm_freezState(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: Qm_melt(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: Qm_rain(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: QmFreez
    REAL(KIND(1d0)),INTENT(in) :: QmRain
    REAL(KIND(1d0)),INTENT(in) :: qn1_bup
    REAL(KIND(1d0)),INTENT(in) :: qn1_ind_snow(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: qn1_S
    REAL(KIND(1d0)),INTENT(in) :: qn1_SF
    REAL(KIND(1d0)),INTENT(in) :: qs
    REAL(KIND(1d0)),INTENT(in) :: RA
    REAL(KIND(1d0)),INTENT(in) :: rainOnSnow(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: resistsurf
    REAL(KIND(1d0)),INTENT(in) :: runoff_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: runoffAGimpervious
    REAL(KIND(1d0)),INTENT(in) :: runoffAGveg
    REAL(KIND(1d0)),INTENT(in) :: runoffPipes
    REAL(KIND(1d0)),INTENT(in) :: runoffSoil_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: runoffWaterBody
    REAL(KIND(1d0)),INTENT(in) :: sfr(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: smd
    REAL(KIND(1d0)),INTENT(in) :: smd_nsurf(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: SnowAlb
    REAL(KIND(1d0)),INTENT(in) :: SnowDens(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: snowDepth(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: SnowRemoval(2)
    REAL(KIND(1d0)),INTENT(in) :: SoilState
    REAL(KIND(1d0)),INTENT(in) :: state(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: state_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: surf_chang_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: swe
    REAL(KIND(1d0)),INTENT(in) :: t2_C
    REAL(KIND(1d0)),INTENT(in) :: tot_chang_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: tsurf
    REAL(KIND(1d0)),INTENT(in) :: Tsurf_ind_snow(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: UStar
    REAL(KIND(1d0)),INTENT(in) :: wu_DecTr
    REAL(KIND(1d0)),INTENT(in) :: wu_EveTr
    REAL(KIND(1d0)),INTENT(in) :: wu_Grass
    REAL(KIND(1d0)),INTENT(in) :: z0m
    REAL(KIND(1d0)),INTENT(in) :: zdm
    REAL(KIND(1d0)),INTENT(in) :: zenith_deg
    REAL(KIND(1d0)),INTENT(in) :: SnowFrac(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: SnowPack(nsurf)


    REAL(KIND(1d0)),INTENT(inout) :: dataOut(ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids)
    REAL(KIND(1d0)),INTENT(inout) :: dataOutSnow(ReadLinesMetdata,ncolumnsDataOutSnow,NumberOfGrids)

    ! INTEGER:: is
    REAL(KIND(1d0)):: LAI_wt

    ! the variables below with '_x' endings stand for 'exported' values
    REAL(KIND(1d0))::ResistSurf_x
    REAL(KIND(1d0))::l_mod_x
    REAL(KIND(1d0))::bulkalbedo
    REAL(KIND(1d0))::qh_x
    REAL(KIND(1d0))::qh_r_x
    REAL(KIND(1d0))::qeOut_x
    REAL(KIND(1d0))::qs_x
    REAL(KIND(1d0))::surf_chang_per_tstep_x
    REAL(KIND(1d0))::tot_chang_per_tstep_x
    REAL(KIND(1d0))::SoilState_x
    REAL(KIND(1d0))::smd_x
    REAL(KIND(1d0))::smd_nsurf_x(nsurf)
    REAL(KIND(1d0))::state_x(nsurf)


    !=====================================================================
    !====================== Prepare data for output ======================

    ! Check surface composition (HCW 21 Jul 2016)
    ! if totally water surface, set output to -999 for columns that do not exist
    !IF(sfr(WaterSurf)==1) THEN
    !   NWstate_per_tstep = NAN    !no non-water surface
    !   smd = NAN                  !no soil store beneath water surface
    !   SoilState = NAN
    !   runoffSoil_per_tstep = NAN
    !   drain_per_tstep = NAN      !no drainage from water surf
    !ENDIF
    !These removed as now all these get a value of 0. LJ 15/06/2017

    ! Remove non-existing surface type from surface and soil outputs   ! Added back in with NANs by HCW 24 Aug 2016
    ! DO is=1,nsurf
    !    IF (sfr(is)<0.00001) THEN
    !       state_x(is)= NAN
    !       smd_nsurf_x(is)= NAN
    !       !runoffOut(is)= NAN
    !       !runoffSoilOut(is)= NAN
    !    ELSE
    !       state_x(is)=state(is)
    !       smd_nsurf_x(is)=smd_nsurf(is)
    !       !runoffOut(is)=runoff(is)
    !       !runoffSoilOut(is)=runoffSoil(is)
    !    ENDIF
    ! ENDDO
    state_x=UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr<0.00001), field=state)
    smd_nsurf_x=UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr<0.00001), field=smd_nsurf)


    ! Remove negative state   !ErrorHint added, then commented out by HCW 16 Feb 2015 as should never occur
    !if(st_per_interval<0) then
    !   call ErrorHint(63,'SUEWS_Calculations: st_per_interval < 0',st_per_interval,NotUsed,NotUsedI)
    !   !st_per_interval=0
    !endif

    !Set limits on output data to avoid formatting issues ---------------------------------
    ! Set limits as +/-9999, depending on sign of value
    ! errorHints -> warnings commented out as writing these slow down the code considerably when there are many instances
    ! IF(ResistSurf > 9999) THEN
    !    !CALL errorHint(6,'rs set to 9999 s m-1 in output; calculated value > 9999 s m-1',ResistSurf,notUsed,notUsedI)
    !    ResistSurf_x=MIN(9999.,ResistSurf)
    ! ENDIF
    ResistSurf_x=MIN(9999.,ResistSurf)

    ! IF(l_mod > 9999) THEN
    !    !CALL errorHint(6,'Lob set to 9999 m in output; calculated value > 9999 m',L_mod,notUsed,notUsedI)
    !    l_mod_x=MIN(9999.,l_mod)
    ! ELSEIF(l_mod < -9999) THEN
    !    !CALL errorHint(6,'Lob set to -9999 m in output; calculated value < -9999 m',L_mod,notUsed,notUsedI)
    !    l_mod_x=MAX(-9999.,l_mod)
    ! ENDIF
    l_mod_x=MAX(MIN(9999.,l_mod), -9999.)

    ! Set NA values   !!Why only these variables??  !!ErrorHints here too - error hints can be very slow here
    ! IF(ABS(qh)>pNAN) qh=NAN
    ! IF(ABS(qh_r)>pNAN) qh_r=NAN
    ! IF(ABS(qeOut)>pNAN) qeOut=NAN
    ! IF(ABS(qs)>pNAN) qs=NAN
    ! ! IF(ABS(ch_per_interval)>pNAN) ch_per_interval=NAN
    ! IF(ABS(surf_chang_per_tstep)>pNAN) surf_chang_per_tstep=NAN
    ! IF(ABS(tot_chang_per_tstep)>pNAN) tot_chang_per_tstep=NAN
    ! IF(ABS(SoilState)>pNAN) SoilState=NAN
    ! IF(ABS(smd)>pNAN) smd=NAN
    ! casting invalid values to NANs
    qh_x                   = set_nan(qh)
    qh_r_x                 = set_nan(qh_r)
    qeOut_x                = set_nan(qeOut)
    qs_x                   = set_nan(qs)
    surf_chang_per_tstep_x = set_nan(surf_chang_per_tstep)
    tot_chang_per_tstep_x  = set_nan(tot_chang_per_tstep)
    SoilState_x            = set_nan(SoilState)
    smd_x                  = set_nan(smd)

    ! this part is now handled in `SUEWS_cal_SoilMoist`
    ! ! If measured smd is used, set components to -999 and smd output to measured one
    ! IF (SMDMethod>0) THEN
    !    !  smd_nsurf=NAN
    !    smd_nsurf_x=NAN
    !    smd_x=xsmd
    ! ENDIF

    ! Calculate areally-weighted LAI
    IF(iy == (iy_prev_t+1) .AND. (id-1) == 0) THEN   !Check for start of next year and avoid using LAI(id-1) as this is at the start of the year
       !  LAI_wt=0
       !  DO is=1,nvegsurf
       !     LAI_wt=LAI_wt+LAI(id_prev_t,is)*sfr(is+2)
       !  ENDDO
       LAI_wt=DOT_PRODUCT(LAI(id_prev_t,:),sfr(1+2:nvegsurf+2))
    ELSE
       !  LAI_wt=0
       !  DO is=1,nvegsurf
       !     LAI_wt=LAI_wt+LAI(id-1,is)*sfr(is+2)
       !  ENDDO
       LAI_wt=DOT_PRODUCT(LAI(id-1,:),sfr(1+2:nvegsurf+2))
    ENDIF

    ! Calculate areally-weighted albedo
    bulkalbedo=DOT_PRODUCT(alb,sfr)
    ! bulkalbedo = 0
    ! DO is=1,nsurf
    !    bulkalbedo = bulkalbedo + alb(is)*sfr(is)
    ! ENDDO

    ! NB: this part needs to be reconsidered for calculation logic:
    ! where to put it? seems to be better placed somewhere near TS, 20170927
    ! ! Save qh and qe for CBL in next iteration
    ! IF(Qh_choice==1) THEN   !use QH and QE from SUEWS
    !    qhforCBL(Gridiv) = qh
    !    qeforCBL(Gridiv) = qeOut
    ! ELSEIF(Qh_choice==2)THEN   !use QH and QE from LUMPS
    !    qhforCBL(Gridiv) = h_mod
    !    qeforCBL(Gridiv) = e_mod
    ! ELSEIF(qh_choice==3)THEN  !use QH and QE from OBS
    !    qhforCBL(Gridiv) = qh_obs
    !    qeforCBL(Gridiv) = qe_obs
    !    IF(qh_obs<-900.OR.qe_obs<-900)THEN  ! observed data has a problem
    !       CALL ErrorHint(22,'Unrealistic observed qh or qe_value.',qh_obs,qe_obs,qh_choice)
    !    ENDIF
    ! ENDIF



    !====================== update output arrays ==============================
    !Define the overall output matrix to be printed out step by step
    dataOut(ir,1:ncolumnsDataOut,Gridiv)=[&
         REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,&   !5
         avkdn,kup,ldown,lup,tsurf,&
         qn1_bup,qf,qs_x,qh_x,qeOut_x,&
         h_mod,e_mod,qh_r_x,&
         precip,ext_wu,ev_per_tstep,runoff_per_tstep,tot_chang_per_tstep_x,&
         surf_chang_per_tstep_x,state_per_tstep,NWstate_per_tstep,drain_per_tstep,smd_x,&
         FlowChange/nsh_real,AdditionalWater,&
         runoffSoil_per_tstep,runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,&
         int_wu,wu_EveTr,wu_DecTr,wu_Grass,&
         smd_nsurf_x(1:nsurf-1),&
         state_x(1:nsurf),&
         zenith_deg,azimuth,bulkalbedo,Fcld,&
         LAI_wt,z0m,zdm,&
         UStar,l_mod_x,ra,ResistSurf_x,&
         Fc,&
         Fc_photo,Fc_respi,Fc_metab,Fc_traff,Fc_build,&
         qn1_SF,qn1_S,SnowAlb,&
         Qm,QmFreez,QmRain,swe,mwh,MwStore,chSnow_per_interval,&
         SnowRemoval(1:2),&
         t2_C,q2_gkg,avU10_ms& ! surface-level diagonostics
         ]

    IF (snowUse==1) THEN
       dataOutSnow(ir,1:ncolumnsDataOutSnow,Gridiv)=[&
            REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,& !5
            SnowPack(1:nsurf),mw_ind(1:nsurf),Qm_melt(1:nsurf),            & !26
            Qm_rain(1:nsurf),Qm_freezState(1:nsurf),snowFrac(1:(nsurf-1)), & !46
            rainOnSnow(1:nsurf),                                           & !53
            qn1_ind_snow(1:nsurf),kup_ind_snow(1:nsurf),freezMelt(1:nsurf),& !74
            MeltWaterStore(1:nsurf),SnowDens(1:nsurf),                     & !88
            snowDepth(1:nsurf),Tsurf_ind_snow(1:nsurf)]
    END IF
    !====================update output arrays end==============================

  END SUBROUTINE SUEWS_update_output
  !========================================================================

  !===============set variable of invalid value to NAN====================================
  ELEMENTAL FUNCTION set_nan(x) RESULT(xx)
    IMPLICIT NONE
    REAL(KIND(1d0)),PARAMETER::pNAN=999
    REAL(KIND(1d0)),PARAMETER::NAN=-999
    REAL(KIND(1d0)),INTENT(in)::x
    REAL(KIND(1d0))::xx

    IF(ABS(x)>pNAN) THEN
       xx=NAN
    ELSE
       xx=x
    ENDIF

  END FUNCTION set_nan
  !========================================================================

END MODULE SUEWS_Driver
