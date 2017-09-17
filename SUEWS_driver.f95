!========================================================================================
! SUEWS driver subroutines
! TS 31 Aug 2017: initial version



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
     a3)


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

  REAL(KIND(1d0)),DIMENSION(6,nsurf+2),INTENT(out)::deltaQi ! storage heat flux of snow surfaces

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


!==========================water balance================================
SUBROUTINE SUEWS_cal_water(&
     Diagnose,&!input
     nsurf,&
     snowUse,&
     NonWaterFraction,addPipes,addImpervious,addVeg,addWaterBody,&
     state,&
     sfr,&
     surf,&
     WaterDist,&
     nsh_real,&
     drain_per_tstep,&  !inout
     drain,&         !output
     AddWaterRunoff,&
     AdditionalWater,runoffPipes,runoff_per_interval,&
     addWater)

  IMPLICIT NONE
  INTEGER,INTENT(in) ::Diagnose
  INTEGER,INTENT(in) ::nsurf
  INTEGER,INTENT(in) ::snowUse

  REAL(KIND(1d0)),INTENT(in) ::NonWaterFraction,addPipes,addImpervious,addVeg,addWaterBody
  REAL(KIND(1d0)),INTENT(in)::nsh_real !nsh cast as a real for use in calculations

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)          ::state
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)          ::sfr
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in)        ::surf
  REAL(KIND(1d0)),DIMENSION(nsurf+1,nsurf-1),INTENT(in)::WaterDist

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: drain         !Drainage of surface type "is" [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: AddWaterRunoff!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: addWater      !water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]

  REAL(KIND(1d0)),INTENT(out)::drain_per_tstep
  REAL(KIND(1d0)),INTENT(out)::AdditionalWater
  REAL(KIND(1d0)),INTENT(out)::runoffPipes
  REAL(KIND(1d0)),INTENT(out)::runoff_per_interval
  INTEGER,PARAMETER ::      WaterSurf = 7
  INTEGER :: is

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

END SUBROUTINE SUEWS_cal_water


!==========================estimate QH================================
SUBROUTINE SUEWS_cal_Hinit(&
     qh_obs,avdens,avcp,h_mod,qn1,dectime,&
     H_init)

  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(in)::qh_obs,avdens,avcp,h_mod,qn1,dectime
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

END SUBROUTINE SUEWS_cal_Hinit


!================latent heat flux and surface wetness===================
! TODO: optimise the structure of this function
SUBROUTINE SUEWS_cal_QE(&
!
  !in
     Diagnose,&
     id,&
     nsurf,&
     tstep,&
     imin,&
     it,&
     ity,&
     snowfractionchoice,&
     ConifSurf,&
     BSoilSurf,&
     BldgSurf,&
     PavSurf,&
     WaterSurf,&
     DecidSurf,&
     GrassSurf,&
     snowCalcSwitch,&
     CRWmin,&
     CRWmax,&
     nsh_real,&
     lvS_J_kg,&
     lv_j_kg,&
     avdens,&
     waterdens,&
     avRh,&
     Press_hPa,&
     Temp_C,&
     RAsnow,&
     psyc_hPa,&
     avcp,&
     sIce_hPa,&
     PervFraction,&
     vegfraction,&
     addimpervious,&
     numPM,&
     s_hPa,&
     ResistSurf,&
     sp,&
     ra,&
     rb,&
     tlv,&
     snowdensmin,&
     precip,&
     PipeCapacity,&
     RunoffToWater,&
    !  runoffAGimpervious,&
    !  runoffAGveg,&
     surpluswaterbody,&
     pin,&
     NonWaterFraction,&
     wu_EveTr,&
     wu_DecTr,&
     wu_Grass,&
     addVeg,&
     addWaterBody,&
     SnowLimPaved,&
     SnowLimBuild,&
     SurfaceArea,&
     drain,&
     WetThresh,&
     stateold,&
     mw_ind,&
     soilstorecap,&
     rainonsnow,&
     freezmelt,&
     freezstate,&
     freezstatevol,&
     Qm_Melt,&
     Qm_rain,&
     Tsurf_ind,&
     sfr,&
     StateLimit,&
     DayofWeek,&
     surf,&
                                !inout:
     runoff_per_interval,&
     snowPack,&
     snowFrac,&
     MeltWaterStore,&
     SnowDepth,&
     iceFrac,&
     addwater,&
     addwaterrunoff,&
     SnowDens,&
     SurplusEvap,&
                                !out:
     snowProf,&
     runoffSnow,&
     runoff,&
     runoffSoil,&
     chang,&
     changSnow,&
     SnowToSurf,&
     state,&
     snowD,&
     ev_snow,&
     soilmoist,&
     SnowRemoval,&
     evap,&
     rss_nsurf,&
     p_mm,&
     rss,&
     qe,&
     state_per_tstep,&
     NWstate_per_tstep,&
     qeOut,&
     swe,&
     ev,&
     chSnow_per_interval,&
     ev_per_tstep,&
     qe_per_tstep,&
     runoff_per_tstep,&
     surf_chang_per_tstep,&
     runoffPipes,&
     mwstore,&
     runoffwaterbody,&
     FlowChange,&
     runoffAGimpervious_m3,&
     runoffAGveg_m3,&
     runoffWaterBody_m3,&
     runoffPipes_m3)

  IMPLICIT NONE
  INTEGER,INTENT(in) ::Diagnose
  INTEGER,INTENT(in) ::id
  INTEGER,INTENT(in) ::nsurf
  INTEGER,INTENT(in) ::tstep
  INTEGER,INTENT(in) ::imin
  INTEGER,INTENT(in) ::it
  INTEGER,INTENT(in) ::ity !Evaporation calculated according to Rutter (1) or Shuttleworth (2)
  INTEGER,INTENT(in) ::snowfractionchoice
  INTEGER,INTENT(in) ::ConifSurf
  INTEGER,INTENT(in) ::BSoilSurf
  INTEGER,INTENT(in) ::BldgSurf
  INTEGER,INTENT(in) ::PavSurf
  INTEGER,INTENT(in) ::WaterSurf
  INTEGER,INTENT(in) ::DecidSurf! surface type code
  INTEGER,INTENT(in) ::GrassSurf! surface type code

  INTEGER,INTENT(in),DIMENSION(nsurf)::snowCalcSwitch

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
  REAL(KIND(1d0)),INTENT(in)::numPM
  REAL(KIND(1d0)),INTENT(in)::s_hPa
  REAL(KIND(1d0)),INTENT(in)::ResistSurf
  REAL(KIND(1d0)),INTENT(in)::sp
  REAL(KIND(1d0)),INTENT(in)::ra
  REAL(KIND(1d0)),INTENT(in)::rb
  REAL(KIND(1d0)),INTENT(in)::tlv
  REAL(KIND(1d0)),INTENT(in)::snowdensmin
  REAL(KIND(1d0)),INTENT(in)::precip
  REAL(KIND(1d0)),INTENT(in)::PipeCapacity
  REAL(KIND(1d0)),INTENT(in)::RunoffToWater
  ! REAL(KIND(1d0)),INTENT(in)::runoffAGimpervious
  ! REAL(KIND(1d0)),INTENT(in)::runoffAGveg
  REAL(KIND(1d0)),INTENT(in)::surpluswaterbody
  REAL(KIND(1d0)),INTENT(in)::pin!Rain per time interval
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
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::stateold
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

  REAL(KIND(1d0)),DIMENSION(366,2),INTENT(in)  ::DayofWeek
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in)::surf

  !Updated status: input and output
  REAL(KIND(1d0)),INTENT(inout)::runoff_per_interval! Total water transported to each grid for grid-to-grid connectivity

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::snowPack
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
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::state
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowD
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::ev_snow
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::soilmoist
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::SnowRemoval
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
  INTEGER :: is
  REAL(KIND(1d0))::runoffAGveg,runoffAGimpervious

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

runoffAGveg=0
runoffAGimpervious=0



  IF(Diagnose==1) WRITE(*,*) 'Calling evap_SUEWS and SoilStore...'
  DO is=1,nsurf   !For each surface in turn
     IF (snowCalcSwitch(is)==1) THEN
        IF (sfr(is)/=0) THEN
           IF(Diagnose==1) WRITE(*,*) 'Calling SnowCalc...'
           CALL SnowCalc(&
                id,& !input
                nsurf,&
                tstep,&
                imin,&
                it,&
                is,&
                snowfractionchoice,&
                ConifSurf,&
                BSoilSurf,&
                BldgSurf,&
                PavSurf,&
                WaterSurf,&
                ity,&
                CRWmin,&
                CRWmax,&
                nsh_real,&
                lvS_J_kg,&
                lv_j_kg,&
                avdens,&
                waterdens,&
                avRh,&
                Press_hPa,&
                Temp_C,&
                RAsnow,&
                psyc_hPa,&
                avcp,&
                sIce_hPa,&
                PervFraction,&
                vegfraction,&
                addimpervious,&
                numPM,&
                s_hPa,&
                ResistSurf,&
                sp,&
                ra,&
                rb,&
                tlv,&
                snowdensmin,&
                precip,&
                PipeCapacity,&
                RunoffToWater,&
                runoffAGimpervious,&
                runoffAGveg,&
                surpluswaterbody,&
                SnowLimPaved,&
                SnowLimBuild,&
                drain,&
                WetThresh,&
                stateold,&
                mw_ind,&
                soilstorecap,&
                rainonsnow,&
                freezmelt,&
                freezstate,&
                freezstatevol,&
                Qm_Melt,&
                Qm_rain,&
                Tsurf_ind,&
                sfr,&
                DayofWeek,&
                surf,&
                snowPack,&!inout
                snowFrac,&
                MeltWaterStore,&
                SnowDepth,&
                iceFrac,&
                addwater,&
                addwaterrunoff,&
                SnowDens,&
                runoffSnow,& ! output
                runoff,&
                runoffSoil,&
                chang,&
                changSnow,&
                SnowToSurf,&
                state,&
                snowD,&
                ev_snow,&
                soilmoist,&
                SnowRemoval,&
                snowProf,&
                swe,&
                ev,&
                chSnow_per_interval,&
                ev_per_tstep,&
                qe_per_tstep,&
                runoff_per_tstep,&
                surf_chang_per_tstep,&
                runoffPipes,&
                mwstore,&
                runoffwaterbody,&
                FlowChange)
        ELSE
           snowFrac(is) = 0
           SnowDens(is) = 0
           SnowPack(is) = 0
        ENDIF
     ELSE

        !Calculates ev [mm]
        CALL Evap_SUEWS(&

                                ! input:
             ity,&!Evaporation calculated according to Rutter (1) or Shuttleworth (2)
             state(is),& ! wetness status
             WetThresh(is),&!When State > WetThresh, rs=0 limit in SUEWS_evap [mm] (specified in input files)
             surf(6,is),& ! = surf(is,6), current storage capacity [mm]
             numPM,&!numerator of P-M eqn
             s_hPa,&!Vapour pressure versus temperature slope in hPa
             psyc_hPa,&!Psychometric constant in hPa
             ResistSurf,&!Surface resistance
             sp,&!Term in calculation of E
             ra,&!Aerodynamic resistance
             rb,&!Boundary layer resistance
             tlv,&!Latent heat of vaporization per timestep [J kg-1 s-1], (tlv=lv_J_kg/tstep_real)

                                ! output:
             rss,&
             ev,&
             qe) ! latent heat flux [W m-2]


        rss_nsurf(is) = rss !Store rss for each surface
        ! CALL soilstore    !Surface water balance and soil store updates (can modify ev, updates state)
        !Surface water balance and soil store updates (can modify ev, updates state)
        CALL soilstore(&
                                ! input:
             nsurf,& ! number of surface types
             is,& ! surface type
             PavSurf,&! surface type code
             BldgSurf,&! surface type code
             WaterSurf,&! surface type code
             ConifSurf,&! surface type code
             BSoilSurf,&! surface type code
             DecidSurf,&! surface type code
             GrassSurf,&! surface type code
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
                                !  inout:
             runoffAGimpervious,&!Above ground runoff from impervious surface [mm] for whole surface area
             surplusWaterBody,&!Extra runoff that goes to water body [mm] as specified by RunoffToWater
             runoffAGveg,&!Above ground runoff from vegetated surfaces [mm] for whole surface area
             runoffPipes,&!Runoff in pipes [mm] for whole surface area
             ev,&!Evaporation
             soilmoist,&!Soil moisture of each surface type [mm]
             SurplusEvap,&!Surplus for evaporation in 5 min timestep
             runoffWaterBody,&!Above ground runoff from water surface [mm] for whole surface area
             runoff_per_interval,&! Total water transported to each grid for grid-to-grid connectivity
                                !  output:
             p_mm,&!Inputs to surface water balance
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

        IF (NonWaterFraction/=0) THEN
           IF (is.NE.WaterSurf) NWstate_per_tstep=NWstate_per_tstep+(state(is)*sfr(is)/NonWaterFraction)
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

SUBROUTINE SUEWS_cal_QH(&
     opt,&
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

  INTEGER,INTENT(in) :: opt ! option for QH calculation: 1, residual; 2, resistance-based

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

  SELECT CASE (opt)
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
