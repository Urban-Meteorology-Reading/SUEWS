!========================================================================================
! SUEWS driver subroutines
! TS 31 Aug 2017: initial version
! TS 02 Oct 2017: added `SUEWS_cal_Main` as the generic wrapper
! TS 03 Oct 2017: added `SUEWS_cal_AnthropogenicEmission`
MODULE SUEWS_Driver
  USE meteo,ONLY:qsatf
  USE AtmMoistStab_module,ONLY:LUMPS_cal_AtmMoist,STAB_lumps,stab_fn_heat,stab_fn_mom
  USE NARP_MODULE,ONLY:NARP_cal_SunPosition
  USE AnOHM_module,ONLY:AnOHM
  USE ESTM_module,ONLY:ESTM
  USE Snow_module,ONLY:SnowCalc,Snow_cal_MeltHeat
  USE DailyState_module,ONLY:SUEWS_cal_DailyState,update_DailyState
  USE WaterDist_module,ONLY:drainage,soilstore,&
       SUEWS_cal_SoilMoist,SUEWS_update_SoilMoist,&
       ReDistributeWater,SUEWS_cal_HorizontalSoilWater,&
       SUEWS_cal_WaterUse
  USE ctrl_output,ONLY:varList
  USE allocateArray,ONLY:&
       ndays,nsurf,nvegsurf,&
       PavSurf,BldgSurf,ConifSurf,DecidSurf,GrassSurf,BSoilSurf,WaterSurf,&
       ivConif,ivDecid,ivGrass,&
       ncolumnsDataOutSUEWS,ncolumnsDataOutSnow,&
       ncolumnsDataOutESTM,ncolumnsDataOutDailyState


  IMPLICIT NONE

CONTAINS
  ! ===================MAIN CALCULATION WRAPPER FOR ENERGY AND WATER FLUX===========
  SUBROUTINE SUEWS_cal_Main(&
       AerodynamicResistanceMethod,AH_MIN,AHProf_tstep,AH_SLOPE_Cooling,& ! input&inout in alphabetical order
       AH_SLOPE_Heating,alb,albDecTr,albEveTr,albGrass,AlbMax_DecTr,&
       AlbMax_EveTr,AlbMax_Grass,AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
       alpha_bioCO2,alpha_enh_bioCO2,alt,avkdn,avRh,avU1,BaseT,BaseTe,&
       BaseTHDD,beta_bioCO2,beta_enh_bioCO2,bldgH,CapMax_dec,CapMin_dec,&
       chAnOHM,cpAnOHM,CRWmax,CRWmin,DayWat,DayWatPer,&
       DecidCap,dectime,DecTreeH,Diagnose,DiagQN,DiagQS,DRAINRT,&
       dt_since_start,dqndt,qn1_av,dqnsdt,qn1_s_av,&
       EF_umolCO2perJ,emis,EmissionsMethod,EnEF_v_Jkm,endDLS,EveTreeH,FAIBldg,&
       FAIDecTree,FAIEveTree,Faut,FcEF_v_kgkm,fcld_obs,FlowChange,&
       FrFossilFuel_Heat,FrFossilFuel_NonHeat,G1,G2,G3,G4,G5,G6,GDD,&
       GDDFull,Gridiv,gsModel,HDD,HumActivity_tstep,&
       IceFrac,id,id_prev_t,Ie_a,Ie_end,Ie_m,Ie_start,imin,&
       InternalWaterUse_h,IrrFracConif,IrrFracDecid,IrrFracGrass,it,ity,&
       iy,iy_prev_t,kkAnOHM,Kmax,LAI,LAICalcYes,LAIMax,LAIMin,LAI_obs,&
       LAIPower,LAIType,lat,ldown_obs,lng,MaxConductance,MaxQFMetab,&
       MeltWaterStore,MetForcingData_grid,MinQFMetab,min_res_bioCO2,&
       NARP_EMIS_SNOW,NARP_TRANS_SITE,NetRadiationMethod,&
       NumCapita,OHM_coef,OHMIncQF,OHM_threshSW,&
       OHM_threshWD,PipeCapacity,PopDensDaytime,&
       PopDensNighttime,PopProf_tstep,PorMax_dec,PorMin_dec,porosity,&
       Precip,PrecipLimit,PrecipLimitAlb,Press_hPa,QF0_BEU,Qf_A,Qf_B,&
       Qf_C,qh_obs,qn1_obs,&
       RadMeltFact,RAINCOVER,RainMaxRes,resp_a,resp_b,&
       RoughLenHeatMethod,RoughLenMomMethod,RunoffToWater,S1,S2,&
       SatHydraulicConduct,SDDFull,sfr,SMDMethod,SnowAlb,SnowAlbMax,&
       SnowAlbMin,snowD,SnowDens,SnowDensMax,SnowDensMin,SnowfallCum,snowFrac,&
       SnowLimBuild,SnowLimPaved,snow_obs,SnowPack,SnowProf,snowUse,SoilDepth,&
       soilmoist,soilstoreCap,StabilityMethod,startDLS,state,StateLimit,&
       StorageHeatMethod,surf,SurfaceArea,Tair24HR,tau_a,tau_f,tau_r,&
       T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,TempMeltFact,TH,&
       theta_bioCO2,timezone,TL,TrafficRate,TrafficUnits,&
       TraffProf_tstep,Ts5mindata_ir,tstep,veg_type,&
       WaterDist,WaterUseMethod,WetThresh,WU_Day,WUProfA_tstep,&
       WUProfM_tstep,xsmd,Z,z0m_in,zdm_in,&
       datetimeLine,dataOutLineSUEWS,dataOutLineSnow,dataOutLineESTM,&!output
       DailyStateLine)!output

    IMPLICIT NONE

    ! input variables
    INTEGER,INTENT(IN)::AerodynamicResistanceMethod
    INTEGER,INTENT(IN)::Diagnose
    INTEGER,INTENT(IN)::DiagQN
    INTEGER,INTENT(IN)::DiagQS
    INTEGER,INTENT(IN)::startDLS
    INTEGER,INTENT(IN)::endDLS
    INTEGER,INTENT(IN)::EmissionsMethod
    INTEGER,INTENT(IN)::Gridiv
    INTEGER,INTENT(IN)::gsModel
    INTEGER,INTENT(IN)::id
    INTEGER,INTENT(IN)::id_prev_t
    INTEGER,INTENT(IN)::Ie_end
    INTEGER,INTENT(IN)::Ie_start
    INTEGER,INTENT(IN)::imin
    INTEGER,INTENT(IN)::it
    INTEGER,INTENT(IN)::ity
    INTEGER,INTENT(IN)::iy
    INTEGER,INTENT(IN)::iy_prev_t
    INTEGER,INTENT(IN)::LAICalcYes
    INTEGER,INTENT(IN)::NetRadiationMethod
    INTEGER,INTENT(IN)::OHMIncQF
    INTEGER,INTENT(IN)::RoughLenHeatMethod
    INTEGER,INTENT(IN)::RoughLenMomMethod
    INTEGER,INTENT(IN)::SMDMethod
    INTEGER,INTENT(IN)::snowUse
    INTEGER,INTENT(IN)::StabilityMethod
    INTEGER,INTENT(IN)::StorageHeatMethod
    INTEGER,INTENT(IN)::tstep
    INTEGER,INTENT(in)::dt_since_start ! time since simulation starts [s]
    INTEGER,INTENT(IN)::veg_type
    INTEGER,INTENT(IN)::WaterUseMethod

    REAL(KIND(1D0)),INTENT(IN)::AlbMax_DecTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMax_EveTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMax_Grass
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_DecTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_EveTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_Grass
    REAL(KIND(1D0)),INTENT(IN)::alt
    REAL(KIND(1D0)),INTENT(IN)::avkdn
    REAL(KIND(1D0)),INTENT(IN)::avRh
    REAL(KIND(1D0)),INTENT(IN)::avU1
    REAL(KIND(1D0)),INTENT(IN)::BaseTHDD
    REAL(KIND(1D0)),INTENT(IN)::bldgH
    REAL(KIND(1D0)),INTENT(IN)::CapMax_dec
    REAL(KIND(1D0)),INTENT(IN)::CapMin_dec
    REAL(KIND(1D0)),INTENT(IN)::CRWmax
    REAL(KIND(1D0)),INTENT(IN)::CRWmin
    REAL(KIND(1D0)),INTENT(IN)::dectime
    REAL(KIND(1D0)),INTENT(IN)::DecTreeH
    REAL(KIND(1D0)),INTENT(IN)::DRAINRT
    REAL(KIND(1D0)),INTENT(IN)::EF_umolCO2perJ
    REAL(KIND(1D0)),INTENT(IN)::EnEF_v_Jkm
    REAL(KIND(1D0)),INTENT(IN)::EveTreeH
    REAL(KIND(1D0)),INTENT(IN)::FAIBldg
    REAL(KIND(1D0)),INTENT(IN)::FAIDecTree
    REAL(KIND(1D0)),INTENT(IN)::FAIEveTree
    REAL(KIND(1D0)),INTENT(IN)::Faut
    REAL(KIND(1D0)),INTENT(IN)::FcEF_v_kgkm
    REAL(KIND(1D0)),INTENT(IN)::fcld_obs
    REAL(KIND(1D0)),INTENT(IN)::FlowChange
    REAL(KIND(1D0)),INTENT(IN)::FrFossilFuel_Heat
    REAL(KIND(1D0)),INTENT(IN)::FrFossilFuel_NonHeat
    REAL(KIND(1D0)),INTENT(IN)::G1
    REAL(KIND(1D0)),INTENT(IN)::G2
    REAL(KIND(1D0)),INTENT(IN)::G3
    REAL(KIND(1D0)),INTENT(IN)::G4
    REAL(KIND(1D0)),INTENT(IN)::G5
    REAL(KIND(1D0)),INTENT(IN)::G6
    REAL(KIND(1D0)),INTENT(IN)::InternalWaterUse_h
    REAL(KIND(1D0)),INTENT(IN)::IrrFracConif
    REAL(KIND(1D0)),INTENT(IN)::IrrFracDecid
    REAL(KIND(1D0)),INTENT(IN)::IrrFracGrass
    REAL(KIND(1D0)),INTENT(IN)::Kmax
    REAL(KIND(1D0)),INTENT(IN)::LAI_obs
    REAL(KIND(1D0)),INTENT(IN)::lat
    REAL(KIND(1D0)),INTENT(IN)::ldown_obs
    REAL(KIND(1D0)),INTENT(IN)::lng
    REAL(KIND(1D0)),INTENT(IN)::MaxQFMetab
    REAL(KIND(1D0)),INTENT(IN)::MinQFMetab
    REAL(KIND(1D0)),INTENT(IN)::NARP_EMIS_SNOW
    REAL(KIND(1D0)),INTENT(IN)::NARP_TRANS_SITE
    REAL(KIND(1D0)),INTENT(IN)::NumCapita
    REAL(KIND(1D0)),INTENT(IN)::PipeCapacity
    REAL(KIND(1D0)),INTENT(IN)::PopDensDaytime
    REAL(KIND(1D0)),INTENT(IN)::PopDensNighttime
    REAL(KIND(1D0)),INTENT(IN)::PorMax_dec
    REAL(KIND(1D0)),INTENT(IN)::PorMin_dec
    REAL(KIND(1D0)),INTENT(IN)::Precip
    REAL(KIND(1D0)),INTENT(IN)::PrecipLimit
    REAL(KIND(1D0)),INTENT(IN)::PrecipLimitAlb
    REAL(KIND(1D0)),INTENT(IN)::Press_hPa
    REAL(KIND(1D0)),INTENT(IN)::qh_obs
    REAL(KIND(1D0)),INTENT(IN)::qn1_obs
    REAL(KIND(1D0)),INTENT(IN)::RadMeltFact
    REAL(KIND(1D0)),INTENT(IN)::RAINCOVER
    REAL(KIND(1D0)),INTENT(IN)::RainMaxRes
    REAL(KIND(1D0)),INTENT(IN)::RunoffToWater
    REAL(KIND(1D0)),INTENT(IN)::S1
    REAL(KIND(1D0)),INTENT(IN)::S2
    REAL(KIND(1D0)),INTENT(IN)::SnowAlbMax
    REAL(KIND(1D0)),INTENT(IN)::SnowAlbMin
    REAL(KIND(1D0)),INTENT(IN)::SnowDensMax
    REAL(KIND(1D0)),INTENT(IN)::SnowDensMin
    REAL(KIND(1D0)),INTENT(IN)::SnowLimBuild
    REAL(KIND(1D0)),INTENT(IN)::SnowLimPaved
    REAL(KIND(1D0)),INTENT(IN)::snow_obs
    REAL(KIND(1D0)),INTENT(IN)::SurfaceArea
    REAL(KIND(1D0)),INTENT(IN)::tau_a
    REAL(KIND(1D0)),INTENT(IN)::tau_f
    REAL(KIND(1D0)),INTENT(IN)::tau_r
    REAL(KIND(1D0)),INTENT(IN)::Temp_C
    REAL(KIND(1D0)),INTENT(IN)::TempMeltFact
    REAL(KIND(1D0)),INTENT(IN)::TH
    REAL(KIND(1D0)),INTENT(IN)::timezone
    REAL(KIND(1D0)),INTENT(IN)::TL
    REAL(KIND(1D0)),INTENT(IN)::TrafficUnits
    REAL(KIND(1D0)),INTENT(IN)::xsmd
    REAL(KIND(1D0)),INTENT(IN)::Z
    REAL(KIND(1D0)),INTENT(IN)::z0m_in
    REAL(KIND(1D0)),INTENT(IN)::zdm_in

    INTEGER,DIMENSION(NVEGSURF),INTENT(IN)::LAIType

    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::AH_MIN
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::AH_SLOPE_Cooling
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::AH_SLOPE_Heating
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::QF0_BEU
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::Qf_A
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::Qf_B
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::Qf_C
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::T_CRITIC_Cooling
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::T_CRITIC_Heating
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)               ::TrafficRate
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)               ::Ie_a
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)               ::Ie_m
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)               ::MaxConductance
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN)               ::DayWat
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN)               ::DayWatPer
    REAL(KIND(1D0)),DIMENSION(nsurf+1),INTENT(IN)         ::OHM_threshSW
    REAL(KIND(1D0)),DIMENSION(nsurf+1),INTENT(IN)         ::OHM_threshWD
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::chAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::cpAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::emis
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::kkAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::SatHydraulicConduct
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::sfr
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::snowD
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::SoilDepth
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::soilstoreCap
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::StateLimit
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)           ::WetThresh
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::alpha_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::alpha_enh_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::BaseT
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::BaseTe
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::beta_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::beta_enh_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::GDDFull
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::LAIMax
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::LAIMin
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::min_res_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::resp_a
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::resp_b
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::SDDFull
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(IN)          ::snowProf
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)        ::theta_bioCO2
    REAL(KIND(1D0)),DIMENSION(4,NVEGSURF),INTENT(IN)      ::LAIPower
    REAL(KIND(1D0)),DIMENSION(nsurf+1,4,3),INTENT(IN)     ::OHM_coef
    REAL(KIND(1D0)),DIMENSION(NSURF+1,NSURF-1),INTENT(IN) ::WaterDist
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)               ::Ts5mindata_ir
    REAL(KIND(1D0)),DIMENSION(:,:),INTENT(IN)             ::MetForcingData_grid
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::AHProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::HumActivity_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::PopProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::TraffProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::WUProfA_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN) ::WUProfM_tstep

    ! inout variables
    REAL(KIND(1D0)),INTENT(INOUT)                             ::SnowfallCum
    REAL(KIND(1D0)),INTENT(INOUT)                             ::SnowAlb
    REAL(KIND(1d0)),INTENT(INOUT)                             ::qn1_av
    REAL(KIND(1d0)),INTENT(INOUT)                             ::dqndt
    REAL(KIND(1d0)),INTENT(INOUT)                             ::qn1_s_av
    REAL(KIND(1d0)),INTENT(INOUT)                             ::dqnsdt
    REAL(KIND(1d0)),DIMENSION(24*3600/tstep),INTENT(INOUT)    ::Tair24HR
    ! REAL(KIND(1D0)),DIMENSION(2*3600/tstep+1),INTENT(INOUT)   ::qn1_av_store_grid
    ! REAL(KIND(1D0)),DIMENSION(3600/tstep),INTENT(INOUT)       ::qn1_store_grid
    ! REAL(KIND(1D0)),DIMENSION(2*3600/tstep+1),INTENT(INOUT)   ::qn1_S_av_store_grid
    ! REAL(KIND(1D0)),DIMENSION(3600/tstep),INTENT(INOUT)       ::qn1_S_store_grid
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albDecTr
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albEveTr
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albGrass
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::DecidCap
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::porosity
    REAL(KIND(1D0)),DIMENSION(0:NDAYS,5),INTENT(INOUT)        ::GDD
    REAL(KIND(1D0)),DIMENSION(0:NDAYS,9),INTENT(INOUT)        ::WU_Day
    REAL(KIND(1D0)),DIMENSION(6,NSURF),INTENT(INOUT)          ::surf
    REAL(KIND(1D0)),DIMENSION(-4:NDAYS,6),INTENT(INOUT)       ::HDD
    REAL(KIND(1D0)),DIMENSION(-4:NDAYS,NVEGSURF),INTENT(INOUT)::LAI
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::alb
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::IceFrac
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::MeltWaterStore
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::SnowDens
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::snowFrac
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::SnowPack
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::soilmoist
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::state

    ! output variables
    REAL(KIND(1D0)),DIMENSION(5),INTENT(OUT)                           ::datetimeLine
    REAL(KIND(1D0)),DIMENSION(ncolumnsDataOutSUEWS-5),INTENT(OUT)      ::dataOutLineSUEWS
    REAL(KIND(1D0)),DIMENSION(ncolumnsDataOutSnow-5),INTENT(OUT)       ::dataOutLineSnow
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutESTM-5),INTENT(OUT)       ::dataOutLineESTM
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutDailyState-5),INTENT(OUT) ::DailyStateLine

    ! local variables
    REAL(KIND(1D0))::a1
    REAL(KIND(1D0))::a2
    REAL(KIND(1D0))::a3
    REAL(KIND(1D0))::AdditionalWater
    REAL(KIND(1D0))::avU10_ms
    REAL(KIND(1D0))::azimuth
    REAL(KIND(1D0))::chSnow_per_interval

    REAL(KIND(1D0))::dens_dry
    REAL(KIND(1d0))::deltaLAI
    REAL(KIND(1D0))::drain_per_tstep
    REAL(KIND(1D0))::Ea_hPa
    REAL(KIND(1D0))::E_mod
    REAL(KIND(1D0))::es_hPa
    REAL(KIND(1D0))::ev
    REAL(KIND(1D0))::ev_per_tstep
    REAL(KIND(1D0))::ext_wu
    REAL(KIND(1D0))::Fc
    REAL(KIND(1D0))::Fc_anthro
    REAL(KIND(1D0))::Fc_biogen
    REAL(KIND(1D0))::Fc_build
    REAL(KIND(1D0))::fcld
    REAL(KIND(1D0))::Fc_metab
    REAL(KIND(1D0))::Fc_photo
    REAL(KIND(1D0))::Fc_respi
    REAL(KIND(1D0))::Fc_traff
    REAL(KIND(1D0))::fwh
    REAL(KIND(1D0))::gsc
    REAL(KIND(1D0))::H_mod
    REAL(KIND(1D0))::int_wu
    REAL(KIND(1D0))::kclear
    REAL(KIND(1D0))::kup
    REAL(KIND(1D0))::ldown
    REAL(KIND(1D0))::lup
    REAL(KIND(1D0))::L_mod
    REAL(KIND(1D0))::mwh
    REAL(KIND(1D0))::mwstore
    REAL(KIND(1D0))::NWstate_per_tstep
    REAL(KIND(1D0))::planF
    REAL(KIND(1D0))::p_mm
    REAL(KIND(1D0))::zL
    REAL(KIND(1D0))::q2_gkg
    REAL(KIND(1D0))::qeOut
    REAL(KIND(1D0))::qe_per_tstep
    REAL(KIND(1D0))::qf
    REAL(KIND(1D0))::QF_SAHP
    REAL(KIND(1D0))::qh
    REAL(KIND(1D0))::qh_residual
    REAL(KIND(1D0))::qh_resist
    REAL(KIND(1D0))::Qm
    REAL(KIND(1D0))::QmFreez
    REAL(KIND(1D0))::QmRain
    REAL(KIND(1D0))::qn1
    REAL(KIND(1D0))::qn1_S
    REAL(KIND(1D0))::qn1_SF
    REAL(KIND(1D0))::qs
    REAL(KIND(1D0))::RA
    REAL(KIND(1D0))::ResistSurf
    REAL(KIND(1D0))::rss
    REAL(KIND(1d0))::runoffAGveg
    REAL(KIND(1d0))::runoffAGimpervious
    REAL(KIND(1D0))::runoff_per_tstep
    REAL(KIND(1D0))::runoffPipes
    REAL(KIND(1D0))::runoffPipes_m3
    REAL(KIND(1D0))::runoffSoil_per_tstep
    REAL(KIND(1D0))::runoffwaterbody
    REAL(KIND(1D0))::runoffWaterBody_m3
    REAL(KIND(1D0))::smd
    REAL(KIND(1D0))::SoilState
    REAL(KIND(1D0))::state_per_tstep
    REAL(KIND(1D0))::surf_chang_per_tstep
    REAL(KIND(1D0))::swe
    REAL(KIND(1D0))::t2_C
    REAL(KIND(1D0))::TempVeg
    REAL(KIND(1D0))::tot_chang_per_tstep
    REAL(KIND(1D0))::TStar
    REAL(KIND(1D0))::tsurf
    REAL(KIND(1D0))::UStar
    REAL(KIND(1D0))::VPD_Pa
    REAL(KIND(1D0))::WUAreaDecTr_m2
    REAL(KIND(1D0))::WUAreaEveTr_m2
    REAL(KIND(1D0))::WUAreaGrass_m2
    REAL(KIND(1D0))::WUAreaTotal_m2
    REAL(KIND(1D0))::wu_DecTr
    REAL(KIND(1D0))::wu_EveTr
    REAL(KIND(1D0))::wu_Grass
    REAL(KIND(1D0))::wu_m3
    REAL(KIND(1D0))::z0m
    REAL(KIND(1D0))::zdm
    REAL(KIND(1D0))::ZENITH_deg
    REAL(KIND(1D0))::Zh


    REAL(KIND(1D0)),DIMENSION(2)::SnowRemoval
    REAL(KIND(1D0)),DIMENSION(NSURF)::chang
    REAL(KIND(1D0)),DIMENSION(NSURF)::changSnow
    REAL(KIND(1D0)),DIMENSION(NSURF)::evap
    REAL(KIND(1D0)),DIMENSION(NSURF)::ev_snow
    REAL(KIND(1D0)),DIMENSION(NSURF)::FreezMelt
    REAL(KIND(1d0)),DIMENSION(nsurf)::kup_ind_snow
    REAL(KIND(1D0)),DIMENSION(NSURF)::mw_ind
    REAL(KIND(1D0)),DIMENSION(NSURF)::Qm_freezState
    REAL(KIND(1D0)),DIMENSION(NSURF)::Qm_melt
    REAL(KIND(1D0)),DIMENSION(NSURF)::Qm_rain
    REAL(KIND(1D0)),DIMENSION(NSURF)::qn1_ind_snow
    REAL(KIND(1D0)),DIMENSION(NSURF)::rainOnSnow
    REAL(KIND(1D0)),DIMENSION(NSURF)::rss_nsurf
    REAL(KIND(1D0)),DIMENSION(NSURF)::runoff
    REAL(KIND(1D0)),DIMENSION(NSURF)::runoffSnow
    REAL(KIND(1D0)),DIMENSION(NSURF)::runoffSoil
    REAL(KIND(1D0)),DIMENSION(NSURF)::smd_nsurf
    REAL(KIND(1D0)),DIMENSION(NSURF)::SnowToSurf
    REAL(KIND(1D0)),DIMENSION(NSURF)::snowDepth

    REAL(KIND(1d0)),DIMENSION(nsurf)::Tsurf_ind_snow

    INTEGER,DIMENSION(NSURF)::snowCalcSwitch
    INTEGER,DIMENSION(3)    ::dayofWeek_id
    INTEGER::DLS

    REAL(KIND(1D0))::avcp
    REAL(KIND(1D0))::avdens
    REAL(KIND(1D0))::dq
    REAL(KIND(1D0))::lv_J_kg
    REAL(KIND(1D0))::lvS_J_kg
    REAL(KIND(1D0))::psyc_hPa
    REAL(KIND(1D0))::qe
    REAL(KIND(1D0))::RAsnow
    REAL(KIND(1D0))::rb
    REAL(KIND(1D0))::runoff_per_interval
    REAL(KIND(1D0))::s_hPa
    REAL(KIND(1D0))::sIce_hpa
    REAL(KIND(1D0))::SoilMoistCap
    REAL(KIND(1D0))::veg_fr
    REAL(KIND(1D0))::VegPhenLumps
    REAL(KIND(1D0))::VPd_hpa
    REAL(KIND(1D0))::vsmd
    REAL(KIND(1D0))::ZZD

    REAL(KIND(1D0)),DIMENSION(NSURF)::deltaQi
    REAL(KIND(1D0)),DIMENSION(NSURF)::drain
    REAL(KIND(1D0)),DIMENSION(NSURF)::FreezState
    REAL(KIND(1D0)),DIMENSION(NSURF)::FreezStateVol
    REAL(KIND(1D0)),DIMENSION(NSURF)::soilmoistOld
    REAL(KIND(1D0)),DIMENSION(NSURF)::stateOld
    REAL(KIND(1D0)),DIMENSION(NSURF)::tsurf_ind

    ! TODO: TS 25 Oct 2017
    ! the `add-*` variables are not used currently as grid-to-grid connection is NOT set up.
    ! set these variables as zero.
    REAL(KIND(1D0))::addImpervious=0
    REAL(KIND(1D0))::addPipes=0
    REAL(KIND(1D0))::addVeg=0
    REAL(KIND(1D0))::addWaterBody=0
    REAL(KIND(1D0)),DIMENSION(NSURF)::AddWater=0
    REAL(KIND(1D0)),DIMENSION(NSURF)::AddWaterRunoff=0

    ! values that are derived from tstep
    INTEGER::nsh ! number of timesteps per hour
    REAL(KIND(1D0))::nsh_real ! nsh in type real
    REAL(KIND(1D0))::tstep_real ! tstep in type real

    ! values that are derived from sfr (surface fractions)
    REAL(KIND(1D0))::VegFraction
    REAL(KIND(1D0))::ImpervFraction
    REAL(KIND(1D0))::PervFraction
    REAL(KIND(1D0))::NonWaterFraction

    ! calculate tstep related VARIABLES
    CALL SUEWS_cal_tstep(&
         tstep,& ! input
         nsh, nsh_real, tstep_real) ! output

    ! calculate surface fraction related VARIABLES
    CALL SUEWS_cal_surf(&
         sfr,& !input
         VegFraction,ImpervFraction,PervFraction,NonWaterFraction) ! output

    ! calculate dayofweek information
    CALL SUEWS_cal_weekday(&
         iy,id,lat,& !input
         dayofWeek_id) !output

    ! calculate dayofweek information
    CALL SUEWS_cal_DLS(&
         id,startDLS,endDLS,& !input
         DLS) !output


    !==============main calculation start=======================
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_RoughnessParameters...'
    IF(Diagnose==1) PRINT*, 'z0m_in =',z0m_in
    CALL SUEWS_cal_RoughnessParameters(&
         RoughLenMomMethod,sfr,&!input
         bldgH,EveTreeH,DecTreeH,&
         porosity(id),FAIBldg,FAIEveTree,FAIDecTree,&
         z0m_in,zdm_in,Z,&
         planF,&!output
         Zh,z0m,zdm,ZZD)

    ! Calculate sun position
    IF(Diagnose==1) WRITE(*,*) 'Calling NARP_cal_SunPosition...'
    CALL NARP_cal_SunPosition(&
         REAL(iy,KIND(1d0)),&!input:
         dectime-tstep/2,&! sun position at middle of timestep before
         timezone,lat,lng,alt,&
         azimuth,zenith_deg)!output:


    !Call the SUEWS_cal_DailyState routine to get surface characteristics ready
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_DailyState...'
    CALL SUEWS_cal_DailyState(&
         iy,id,it,imin,tstep,DayofWeek_id,&!input
         WaterUseMethod,snowUse,Ie_start,Ie_end,&
         LAICalcYes,LAIType,&
         nsh_real,avkdn,Temp_C,Precip,BaseTHDD,&
         lat,Faut,LAI_obs,tau_a,tau_f,tau_r,&
         SnowDensMax,SnowDensMin,SnowAlbMin,&
         AlbMax_DecTr,AlbMax_EveTr,AlbMax_Grass,&
         AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
         CapMax_dec,CapMin_dec,PorMax_dec,PorMin_dec,&
         Ie_a,Ie_m,DayWatPer,DayWat,SnowPack,&
         BaseT,BaseTe,GDDFull,SDDFull,LAIMin,LAIMax,LAIPower,&
         SnowAlb,DecidCap,albDecTr,albEveTr,albGrass,&!inout
         porosity,GDD,HDD,SnowDens,LAI,WU_Day,&
         deltaLAI)!output


    !Calculation of density and other water related parameters
    IF(Diagnose==1) WRITE(*,*) 'Calling LUMPS_cal_AtmMoist...'
    CALL LUMPS_cal_AtmMoist(&
         Temp_C,Press_hPa,avRh,dectime,&! input:
         lv_J_kg,lvS_J_kg,&! output:
         es_hPa,Ea_hPa,VPd_hpa,VPD_Pa,dq,dens_dry,avcp,avdens)


    !======== Calculate soil moisture =========
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_update_SoilMoist...'
    CALL SUEWS_update_SoilMoist(&
         NonWaterFraction,&!input
         soilstoreCap,sfr,soilmoist,&
         SoilMoistCap,SoilState,&!output
         vsmd,smd)


    ! ===================NET ALLWAVE RADIATION================================
    CALL SUEWS_cal_Qn(&
         NetRadiationMethod,snowUse,id,&!input
         Diagnose,snow_obs,ldown_obs,fcld_obs,&
         dectime,ZENITH_deg,avKdn,Temp_C,avRH,Ea_hPa,qn1_obs,&
         SnowAlb,DiagQN,&
         NARP_TRANS_SITE,NARP_EMIS_SNOW,IceFrac,sfr,emis,&
         alb,albDecTr,DecidCap,albEveTr,albGrass,surf,&!inout
         snowFrac,ldown,fcld,&!output
         qn1,qn1_SF,qn1_S,kclear,kup,lup,tsurf,&
         qn1_ind_snow,kup_ind_snow,Tsurf_ind_snow,Tsurf_ind)


    ! ===================ANTHROPOGENIC HEAT FLUX================================
    CALL SUEWS_cal_AnthropogenicEmission(&
         AH_MIN,AHProf_tstep,AH_SLOPE_Cooling,AH_SLOPE_Heating,alpha_bioCO2,&
         alpha_enh_bioCO2,avkdn,beta_bioCO2,beta_enh_bioCO2,DayofWeek_id,&
         Diagnose,DLS,EF_umolCO2perJ,EmissionsMethod,EnEF_v_Jkm,Fc,Fc_anthro,Fc_biogen,&
         Fc_build,FcEF_v_kgkm,Fc_metab,Fc_photo,Fc_respi,Fc_traff,FrFossilFuel_Heat,&
         FrFossilFuel_NonHeat,HDD,HumActivity_tstep,id,imin,it,LAI, LaiMax,LaiMin,&
         MaxQFMetab,MinQFMetab,min_res_bioCO2,nsh,NumCapita,&
         PopDensDaytime,PopDensNighttime,PopProf_tstep,QF,QF0_BEU,Qf_A,Qf_B,Qf_C,QF_SAHP,&
         resp_a,resp_b,sfr,snowFrac,T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,&
         theta_bioCO2,TrafficRate,TrafficUnits,TraffProf_tstep)


    ! =================STORAGE HEAT FLUX=======================================
    CALL SUEWS_cal_Qs(&
         StorageHeatMethod,OHMIncQF,Gridiv,&!input
         id,tstep,dt_since_start,Diagnose,sfr,&
         OHM_coef,OHM_threshSW,OHM_threshWD,&
         soilmoist,soilstoreCap,state,nsh,SnowUse,DiagQS,&
         HDD,MetForcingData_grid,Ts5mindata_ir,qf,qn1,&
         avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown,&
         bldgh,alb,emis,cpAnOHM,kkAnOHM,chAnOHM,EmissionsMethod,&
         Tair24HR,qn1_av,dqndt,qn1_s_av,dqnsdt,&!inout
surf,&
         qn1_S,snowFrac,dataOutLineESTM,qs,&!output
         deltaQi,a1,a2,a3)


    !==================Energy related to snow melting/freezing processes=======
    IF(Diagnose==1) WRITE(*,*) 'Calling MeltHeat'
    CALL Snow_cal_MeltHeat(&
         snowUse,&!input
         lvS_J_kg,lv_J_kg,tstep_real,RadMeltFact,TempMeltFact,SnowAlbMax,&
         SnowDensMin,Temp_C,Precip,PrecipLimit,PrecipLimitAlb,&
         nsh_real,sfr,Tsurf_ind,Tsurf_ind_snow,state,qn1_ind_snow,&
         kup_ind_snow,Meltwaterstore,deltaQi,&
         SnowPack,snowFrac,SnowAlb,SnowDens,SnowfallCum,&!inout
         mwh,fwh,Qm,QmFreez,QmRain,&! output
         veg_fr,snowCalcSwitch,Qm_melt,Qm_freezState,Qm_rain,FreezMelt,&
         FreezState,FreezStateVol,rainOnSnow,SnowDepth,mw_ind,&
         dataOutLineSnow)!output



    !==========================Turbulent Fluxes================================
    IF(Diagnose==1) WRITE(*,*) 'Calling LUMPS_cal_QHQE...'
    !Calculate QH and QE from LUMPS
    CALL LUMPS_cal_QHQE(&
         veg_type,& !input
         snowUse,id,qn1,qf,qs,Qm,Temp_C,Veg_Fr,avcp,Press_hPa,lv_J_kg,&
         tstep_real,DRAINRT,nsh_real,&
         Precip,RainMaxRes,RAINCOVER,sfr,LAI,LAImax,LAImin,&
         H_mod,E_mod,psyc_hPa,s_hPa,sIce_hpa,TempVeg,VegPhenLumps)!output


    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_WaterUse...'
    !Gives the external and internal water uses per timestep
    CALL SUEWS_cal_WaterUse(&
         nsh_real,& ! input:
         SurfaceArea,sfr,&
         IrrFracConif,IrrFracDecid,IrrFracGrass,&
         dayofWeek_id,WUProfA_tstep,WUProfM_tstep,&
         InternalWaterUse_h,HDD(id-1,:),WU_Day(id-1,:),&
         WaterUseMethod,NSH,it,imin,DLS,&
         WUAreaEveTr_m2,WUAreaDecTr_m2,& ! output:
         WUAreaGrass_m2,WUAreaTotal_m2,&
         wu_EveTr,wu_DecTr,wu_Grass,wu_m3,int_wu,ext_wu)


    !===============Resistance Calculations=======================
    CALL SUEWS_cal_Resistance(&
         StabilityMethod,&!input:
         Diagnose,AerodynamicResistanceMethod,RoughLenHeatMethod,snowUse,&
         id,it,gsModel,SMDMethod,&
         qh_obs,avdens,avcp,h_mod,qn1,dectime,zzd,z0m,zdm,&
         avU1,Temp_C,VegFraction,avkdn,&
         Kmax,&
         g1,g2,g3,g4,&
         g5,g6,s1,s2,&
         th,tl,&
         dq,xsmd,vsmd,MaxConductance,LAIMax,LAI(id-1,:),snowFrac,sfr,&
         UStar,TStar,L_mod,&!output
         zL,gsc,ResistSurf,RA,RAsnow,rb)


    !============= calculate water balance =============
    CALL SUEWS_cal_Water(&
         Diagnose,&!input
         snowUse,NonWaterFraction,addPipes,addImpervious,addVeg,addWaterBody,&
         state,soilmoist,sfr,surf,WaterDist,nsh_real,&
         drain_per_tstep,&  !output
         drain,AddWaterRunoff,&
         AdditionalWater,runoffPipes,runoff_per_interval,&
         AddWater,stateOld,soilmoistOld)
    !============= calculate water balance end =============

    !======== Evaporation and surface state ========
    CALL SUEWS_cal_QE(&
         Diagnose,&!input
         id,tstep,imin,it,ity,snowCalcSwitch,DayofWeek_id,CRWmin,CRWmax,&
         nsh_real,dectime,lvS_J_kg,lv_j_kg,avdens,avRh,Press_hPa,Temp_C,&
         RAsnow,psyc_hPa,avcp,sIce_hPa,&
         PervFraction,vegfraction,addimpervious,qn1_SF,qf,qs,vpd_hPa,s_hPa,&
         ResistSurf,RA,rb,tstep_real,snowdensmin,precip,PipeCapacity,RunoffToWater,&
         NonWaterFraction,wu_EveTr,wu_DecTr,wu_Grass,addVeg,addWaterBody,SnowLimPaved,SnowLimBuild,&
         SurfaceArea,FlowChange,drain,WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
         freezmelt,freezstate,freezstatevol,Qm_Melt,Qm_rain,Tsurf_ind,sfr,&
         StateLimit,AddWater,addwaterrunoff,surf,snowD,&
         runoff_per_interval,state,soilmoist,SnowPack,snowFrac,MeltWaterStore,&! inout:
         iceFrac,SnowDens,&
         snowProf,& ! output:
         runoffSnow,runoff,runoffSoil,chang,changSnow,&
         snowDepth,SnowToSurf,ev_snow,SnowRemoval,&
         evap,rss_nsurf,p_mm,rss,qe,state_per_tstep,NWstate_per_tstep,qeOut,&
         swe,ev,chSnow_per_interval,ev_per_tstep,qe_per_tstep,runoff_per_tstep,&
         surf_chang_per_tstep,runoffPipes,mwstore,runoffwaterbody,&
         runoffAGveg,runoffAGimpervious,runoffWaterBody_m3,runoffPipes_m3)
    !======== Evaporation and surface state end========

    !============ Sensible heat flux ===============
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_QH...'
    CALL SUEWS_cal_QH(&
         1,&
         qn1,qf,QmRain,qeOut,qs,QmFreez,qm,avdens,avcp,tsurf,Temp_C,RA,&
         qh,qh_residual,qh_resist)!output
    !============ Sensible heat flux end===============

    !=== Horizontal movement between soil stores ===
    ! Now water is allowed to move horizontally between the soil stores
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_HorizontalSoilWater...'
    CALL SUEWS_cal_HorizontalSoilWater(&
         sfr,&! input: ! surface fractions
         SoilStoreCap,&!Capacity of soil store for each surface [mm]
         SoilDepth,&!Depth of sub-surface soil store for each surface [mm]
         SatHydraulicConduct,&!Saturated hydraulic conductivity for each soil subsurface [mm s-1]
         SurfaceArea,&!Surface area of the study area [m2]
         NonWaterFraction,&! sum of surface cover fractions for all except water surfaces
         tstep_real,& !tstep cast as a real for use in calculations
         SoilMoist,&! inout:!Soil moisture of each surface type [mm]
         runoffSoil,&!Soil runoff from each soil sub-surface [mm]
         runoffSoil_per_tstep&!  output:!Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)
         )

    !========== Calculate soil moisture ============
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_SoilMoist...'
    CALL SUEWS_cal_SoilMoist(&
         SMDMethod,xsmd,NonWaterFraction,SoilMoistCap,&!input
         SoilStoreCap,surf_chang_per_tstep,&
         soilmoist,soilmoistOld,sfr,&
         smd,smd_nsurf,tot_chang_per_tstep,SoilState)!output


    !============ surface-level diagonostics ===============
    IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_Diagnostics...'
    CALL SUEWS_cal_Diagnostics(&
         dectime,&!input
         avU1,Temp_C,&
                                ! NB: resistance-based QH is used to calculate diagnostics
                                ! as tsurf is assumed to be equal to Tair during night
                                ! implying a constant near-surface air temperature profile.
         tsurf,qh_resist,&
         Press_hPa,qe,&
         veg_fr,z0m,avdens,avcp,lv_J_kg,tstep_real,&
         RoughLenHeatMethod,StabilityMethod,&
         avU10_ms,t2_C,q2_gkg,L_MOD)!output
    !============ surface-level diagonostics end ===============


    !==============main calculation end=======================

    !==============translation of  output variables into output array===========
    CALL SUEWS_update_outputLine(&
         AdditionalWater,alb,avkdn,avU10_ms,azimuth,&!input
         chSnow_per_interval,dectime,&
         drain_per_tstep,E_mod,ev_per_tstep,ext_wu,Fc,Fc_build,fcld,&
         Fc_metab,Fc_photo,Fc_respi,Fc_traff,FlowChange,&
         h_mod,id,id_prev_t,imin,int_wu,it,iy,iy_prev_t,&
         kup,LAI,ldown,l_mod,lup,mwh,MwStore,&
         nsh_real,NWstate_per_tstep,Precip,q2_gkg,&
         qeOut,qf,qh,qh_resist,Qm,QmFreez,&
         QmRain,qn1,qn1_S,qn1_SF,qs,RA,&
         resistsurf,runoffAGimpervious,runoffAGveg,&
         runoff_per_tstep,runoffPipes,runoffSoil_per_tstep,&
         runoffWaterBody,sfr,smd,smd_nsurf,SnowAlb,SnowRemoval,&
         state,state_per_tstep,surf_chang_per_tstep,swe,t2_C,&
         tot_chang_per_tstep,tsurf,UStar,wu_DecTr,&
         wu_EveTr,wu_Grass,z0m,zdm,zenith_deg,&
         datetimeLine,dataOutLineSUEWS)!output

    ! model state:

    ! daily state:
    CALL update_DailyState(&
         iy,id,it,imin,nsh_real,&!input
         GDD,HDD,LAI,&
         DecidCap,albDecTr,albEveTr,albGrass,porosity,&
         WU_Day,&
         deltaLAI,VegPhenLumps,&
         SnowAlb,SnowDens,&
         a1,a2,a3,&
         DailyStateLine)!out

    !==============translation end ================

  END SUBROUTINE SUEWS_cal_Main
  ! ================================================================================

  ! ===================ANTHROPOGENIC HEAT + CO2 FLUX================================
  SUBROUTINE SUEWS_cal_AnthropogenicEmission(&
       AH_MIN,AHProf_tstep,AH_SLOPE_Cooling,AH_SLOPE_Heating,alpha_bioCO2,&
       alpha_enh_bioCO2,avkdn,beta_bioCO2,beta_enh_bioCO2,dayofWeek_id,&
       Diagnose,DLS,EF_umolCO2perJ,EmissionsMethod,EnEF_v_Jkm,Fc,Fc_anthro,Fc_biogen,&
       Fc_build,FcEF_v_kgkm,Fc_metab,Fc_photo,Fc_respi,Fc_traff,FrFossilFuel_Heat,&
       FrFossilFuel_NonHeat,HDD,HumActivity_tstep,id,imin,it,LAI, LaiMax,LaiMin,&
       MaxQFMetab,MinQFMetab,min_res_bioCO2,nsh,NumCapita,&
       PopDensDaytime,PopDensNighttime,PopProf_tstep,QF,QF0_BEU,Qf_A,Qf_B,Qf_C,QF_SAHP,&
       resp_a,resp_b,sfr,snowFrac,T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,&
       theta_bioCO2,TrafficRate,TrafficUnits,TraffProf_tstep)

    IMPLICIT NONE

    INTEGER,INTENT(in)::Diagnose
    INTEGER,INTENT(in)::EmissionsMethod
    INTEGER,INTENT(in)::id
    INTEGER,INTENT(in)::it
    INTEGER,INTENT(in)::imin
    INTEGER,INTENT(in)::DLS
    INTEGER,INTENT(in)::nsh
    ! INTEGER,INTENT(in)::notUsedI
    INTEGER,DIMENSION(3),INTENT(in)::dayofWeek_id
    REAL(KIND(1d0)),DIMENSION(-4:ndays, 6),INTENT(in)::HDD
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::Qf_A
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::Qf_B
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::Qf_C
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::AH_MIN
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::AH_SLOPE_Heating
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::AH_SLOPE_Cooling
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::T_CRITIC_Heating
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::T_CRITIC_Cooling
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::TrafficRate
    REAL(KIND(1d0)),DIMENSION(2),INTENT(in)::QF0_BEU
    REAL(KIND(1d0)),DIMENSION(24*nsh,2),INTENT(in)::AHProf_tstep
    REAL(KIND(1d0)),DIMENSION(24*nsh,2),INTENT(in)::HumActivity_tstep
    REAL(KIND(1d0)),DIMENSION(24*nsh,2),INTENT(in)::TraffProf_tstep
    REAL(KIND(1d0)),DIMENSION(24*nsh,2),INTENT(in)::PopProf_tstep
    REAL(KIND(1D0)),INTENT(in)::EF_umolCO2perJ
    REAL(KIND(1D0)),INTENT(in)::FcEF_v_kgkm
    REAL(KIND(1D0)),INTENT(in)::EnEF_v_Jkm
    REAL(KIND(1D0)),INTENT(in)::TrafficUnits
    REAL(KIND(1D0)),INTENT(in)::FrFossilFuel_Heat
    REAL(KIND(1D0)),INTENT(in)::FrFossilFuel_NonHeat
    REAL(KIND(1D0)),INTENT(in)::MinQFMetab
    REAL(KIND(1D0)),INTENT(in)::MaxQFMetab
    REAL(KIND(1D0)),INTENT(in)::NumCapita
    REAL(KIND(1D0)),INTENT(in)::PopDensDaytime
    REAL(KIND(1D0)),INTENT(in)::PopDensNighttime
    REAL(KIND(1D0)),INTENT(in)::Temp_C
    REAL(KIND(1D0)),INTENT(out)::QF
    REAL(KIND(1D0)),INTENT(out)::QF_SAHP
    REAL(KIND(1D0)),INTENT(out)::Fc_anthro
    REAL(KIND(1D0)),INTENT(out)::Fc_metab
    REAL(KIND(1D0)),INTENT(out)::Fc_traff
    REAL(KIND(1D0)),INTENT(out)::Fc_build
    REAL(KIND(1d0)),INTENT(in)::avkdn
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::snowFrac
    REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf),INTENT(in)::LAI
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::LaiMin
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in):: LaiMax
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::alpha_bioCO2
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::beta_bioCO2
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::theta_bioCO2
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::alpha_enh_bioCO2
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::beta_enh_bioCO2
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::resp_a
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::resp_b
    REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::min_res_bioCO2
    REAL(KIND(1D0)),INTENT(out)::Fc_biogen
    REAL(KIND(1D0)),INTENT(out)::Fc_respi
    REAL(KIND(1D0)),INTENT(out)::Fc_photo
    REAL(KIND(1D0)),INTENT(out)::Fc

    INTEGER,PARAMETER :: notUsedI=-999
    REAL(KIND(1D0)),PARAMETER::notUsed=-999


    !ih=it-DLS           !Moved to subroutine AnthropogenicEmissions MH 29 June 2017
    !IF(ih<0) ih=23

    IF(EmissionsMethod>0 .AND. EmissionsMethod<=6)THEN
       IF(Diagnose==1) WRITE(*,*) 'Calling AnthropogenicEmissions...'
       CALL AnthropogenicEmissions(&
            EmissionsMethod,&
            id,it,imin,DLS,nsh,dayofWeek_id,ndays,&
            EF_umolCO2perJ,FcEF_v_kgkm,EnEF_v_Jkm,TrafficUnits,&
            FrFossilFuel_Heat,FrFossilFuel_NonHeat,&
            MinQFMetab,MaxQFMetab,&
            NumCapita,PopDensDaytime,PopDensNighttime,&
            Temp_C,HDD,Qf_A,Qf_B,Qf_C,&
            AH_MIN,AH_SLOPE_Heating,AH_SLOPE_Cooling,&
            T_CRITIC_Heating,T_CRITIC_Cooling,&
            TrafficRate,&
            QF0_BEU,QF_SAHP,&
            Fc_anthro,Fc_metab,Fc_traff,Fc_build,&
            AHProf_tstep,HumActivity_tstep,TraffProf_tstep,PopProf_tstep,&
            notUsed,notUsedI)

       !  qn1_bup=qn1
       Fc_anthro=0
       Fc_metab=0
       Fc_traff=0
       Fc_build=0
       Fc_biogen=0
       Fc_respi=0
       Fc_photo=0
    ELSEIF(EmissionsMethod>=11)THEN
       IF(Diagnose==1) WRITE(*,*) 'Calling AnthropogenicEmissions...'
       CALL AnthropogenicEmissions(&
            EmissionsMethod,&
            id,it,imin,DLS,nsh,dayofWeek_id,ndays,&
            EF_umolCO2perJ,FcEF_v_kgkm,EnEF_v_Jkm,TrafficUnits,&
            FrFossilFuel_Heat,FrFossilFuel_NonHeat,&
            MinQFMetab,MaxQFMetab,&
            NumCapita,PopDensDaytime,PopDensNighttime,&
            Temp_C,HDD,Qf_A,Qf_B,Qf_C,&
            AH_MIN,AH_SLOPE_Heating,AH_SLOPE_Cooling,&
            T_CRITIC_Heating,T_CRITIC_Cooling,&
            TrafficRate,&
            QF0_BEU,QF_SAHP,&
            Fc_anthro,Fc_metab,Fc_traff,Fc_build,&
            AHProf_tstep,HumActivity_tstep,TraffProf_tstep,PopProf_tstep,&
            notUsed,notUsedI)

       !  qn1_bup=qn1
    ELSEIF(EmissionsMethod==0)THEN
       IF(Diagnose==1) WRITE(*,*) 'Calling AnthropogenicEmissions...'
       CALL AnthropogenicEmissions(&
            EmissionsMethod,&
            id,it,imin,DLS,nsh,dayofWeek_id,ndays,&
            EF_umolCO2perJ,FcEF_v_kgkm,EnEF_v_Jkm,TrafficUnits,&
            FrFossilFuel_Heat,FrFossilFuel_NonHeat,&
            MinQFMetab,MaxQFMetab,&
            NumCapita,PopDensDaytime,PopDensNighttime,&
            Temp_C,HDD,Qf_A,Qf_B,Qf_C,&
            AH_MIN,AH_SLOPE_Heating,AH_SLOPE_Cooling,&
            T_CRITIC_Heating,T_CRITIC_Cooling,&
            TrafficRate,&
            QF0_BEU,QF_SAHP,&
            Fc_anthro,Fc_metab,Fc_traff,Fc_build,&
            AHProf_tstep,HumActivity_tstep,TraffProf_tstep,PopProf_tstep,&
            notUsed,notUsedI)

       !  qn1_bup=qn1
       !  qn1=qn1+qf
    ELSE
       CALL ErrorHint(73,'RunControl.nml:EmissionsMethod unusable',notUsed,notUsed,EmissionsMethod)
    ENDIF
    ! -- qn1 is now QSTAR+QF (net all-wave radiation + anthropogenic heat flux)
    ! -- qn1_bup is QSTAR only
    IF(EmissionsMethod>=1) qf = QF_SAHP

    IF(EmissionsMethod>=11) THEN
       ! Calculate CO2 fluxes from biogenic components
       IF(Diagnose==1) WRITE(*,*) 'Calling CO2_biogen...'
       CALL CO2_biogen(EmissionsMethod,id,ndays,ivConif,ivDecid,ivGrass,ConifSurf,DecidSurf,GrassSurf,BSoilSurf,&
            snowFrac,nsurf,NVegSurf,avkdn,Temp_C,sfr,LAI,LaiMin,LaiMax,&
            alpha_bioCO2,beta_bioCO2,theta_bioCO2,alpha_enh_bioCO2,beta_enh_bioCO2,&
            resp_a,resp_b,min_res_bioCO2,Fc_biogen,Fc_respi,Fc_photo,&
            notUsed,notUsedI)
    ENDIF
    ! Sum anthropogenic and biogenic CO2 flux components to find overall CO2 flux
    Fc = Fc_anthro + Fc_biogen

    ! =================STORAGE HEAT FLUX=======================================

  END SUBROUTINE SUEWS_cal_AnthropogenicEmission
  ! ================================================================================

  !=============net all-wave radiation=====================================
  SUBROUTINE SUEWS_cal_Qn(&
       NetRadiationMethod,snowUse,id,&!input
       Diagnose,snow_obs,ldown_obs,fcld_obs,&
       dectime,ZENITH_deg,avKdn,Temp_C,avRH,ea_hPa,qn1_obs,&
       SnowAlb,DiagQN,&
       NARP_TRANS_SITE,NARP_EMIS_SNOW,IceFrac,sfr,emis,&
       alb,albDecTr,DecidCap,albEveTr,albGrass,surf,&!inout
       snowFrac,ldown,fcld,&!output
       qn1,qn1_SF,qn1_S,kclear,kup,lup,tsurf,&
       qn1_ind_snow,kup_ind_snow,Tsurf_ind_snow,Tsurf_ind)
    USE NARP_MODULE, ONLY: RadMethod,NARP

    IMPLICIT NONE
    ! INTEGER,PARAMETER ::nsurf     = 7 ! number of surface types
    ! INTEGER,PARAMETER ::ConifSurf = 3 !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER ::DecidSurf = 4 !New surface classes: Grass = 5th/7 surfaces
    ! INTEGER,PARAMETER ::GrassSurf = 5

    INTEGER,INTENT(in)::NetRadiationMethod
    INTEGER,INTENT(in)::snowUse
    INTEGER,INTENT(in)::id
    INTEGER,INTENT(in)::Diagnose
    INTEGER,INTENT(in)::DiagQN

    REAL(KIND(1d0)),INTENT(in)::snow_obs
    REAL(KIND(1d0)),INTENT(in)::ldown_obs
    REAL(KIND(1d0)),INTENT(in)::fcld_obs
    REAL(KIND(1d0)),INTENT(in)::dectime
    REAL(KIND(1d0)),INTENT(in)::ZENITH_deg
    REAL(KIND(1d0)),INTENT(in)::avKdn
    REAL(KIND(1d0)),INTENT(in)::Temp_C
    REAL(KIND(1d0)),INTENT(in)::avRH
    REAL(KIND(1d0)),INTENT(in)::ea_hPa
    REAL(KIND(1d0)),INTENT(in)::qn1_obs
    REAL(KIND(1d0)),INTENT(in)::SnowAlb
    REAL(KIND(1d0)),INTENT(in)::NARP_EMIS_SNOW
    REAL(KIND(1d0)),INTENT(in)::NARP_TRANS_SITE


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

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out) ::qn1_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out) ::kup_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out) ::Tsurf_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: tsurf_ind


    REAL(KIND(1d0)),DIMENSION(nsurf):: lup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf):: kup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf):: qn1_ind

    REAL(KIND(1d0)),PARAMETER::NAN=-999
    INTEGER :: NetRadiationMethodX
    INTEGER::AlbedoChoice,ldown_option


    CALL RadMethod(&
         NetRadiationMethod,&!inout
         snowUse,&!input
         NetRadiationMethodX,AlbedoChoice,ldown_option)!output

    IF(NetRadiationMethodX>0)THEN

       ! IF (snowUse==0) snowFrac=snow_obs
       IF (snowUse==0) snowFrac=0

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
            NARP_TRANS_SITE,NARP_EMIS_SNOW,&
            dectime,ZENITH_deg,avKdn,Temp_C,avRH,ea_hPa,qn1_obs,&
            SnowAlb,&
            AlbedoChoice,ldown_option,NetRadiationMethodX,DiagQN,&
            qn1,qn1_SF,qn1_S,kclear,kup,LDown,lup,fcld,tsurf,&! output:
            qn1_ind_snow,kup_ind_snow,Tsurf_ind_snow,Tsurf_ind)

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
       StorageHeatMethod,OHMIncQF,Gridiv,&!input
       id,tstep,dt_since_start,Diagnose,sfr,&
       OHM_coef,OHM_threshSW,OHM_threshWD,&
       soilmoist,soilstoreCap,state,nsh,SnowUse,DiagQS,&
       HDD,MetForcingData_grid,Ts5mindata_ir,qf,qn1,&
       avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown,&
       bldgh,alb,emis,cpAnOHM,kkAnOHM,chAnOHM,EmissionsMethod,&
       Tair24HR,qn1_av,dqndt,qn1_s_av,dqnsdt,&!inout
       surf,&
       qn1_S,snowFrac,dataOutLineESTM,qs,&!output
       deltaQi,a1,a2,a3)

    IMPLICIT NONE

    INTEGER,INTENT(in)  ::StorageHeatMethod
    INTEGER,INTENT(in)  ::OHMIncQF
    INTEGER,INTENT(in)  ::Gridiv
    INTEGER,INTENT(in)  ::id
    INTEGER,INTENT(in)  ::tstep ! time step [s]
    INTEGER, INTENT(in) ::dt_since_start  ! time since simulation starts [s]
    INTEGER,INTENT(in)  ::Diagnose
    INTEGER,INTENT(in)  ::nsh              ! number of timesteps in one hour
    INTEGER,INTENT(in)  ::SnowUse          ! option for snow related calculations
    INTEGER,INTENT(in)  ::DiagQS           ! diagnostic option
    INTEGER,INTENT(in)  :: EmissionsMethod !< AnthropHeat option [-]


    REAL(KIND(1d0)),INTENT(in)::OHM_coef(nsurf+1,4,3)                 ! OHM coefficients
    REAL(KIND(1d0)),INTENT(in)::OHM_threshSW(nsurf+1) ! OHM thresholds
    REAL(KIND(1d0)),INTENT(in)::OHM_threshWD(nsurf+1) ! OHM thresholds
    REAL(KIND(1d0)),INTENT(in)::soilmoist(nsurf)                ! soil moisture
    REAL(KIND(1d0)),INTENT(in)::soilstoreCap(nsurf)             ! capacity of soil store
    REAL(KIND(1d0)),INTENT(in)::state(nsurf) ! wetness status


    REAL(KIND(1d0)),DIMENSION(-4:ndays, 6),INTENT(in)::HDD
    REAL(KIND(1d0)),INTENT(in)::qf
    REAL(KIND(1d0)),INTENT(in)::qn1
    REAL(KIND(1d0)),INTENT(in)::avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown
    REAL(KIND(1d0)),INTENT(in)::bldgh

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::alb  !< albedo [-]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::emis !< emissivity [-]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::cpAnOHM   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::kkAnOHM   !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::chAnOHM   !< bulk transfer coef [J m-3 K-1]

    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::MetForcingData_grid !< met forcing array of grid

    REAL(KIND(1d0)),DIMENSION(:),INTENT(in)::Ts5mindata_ir

    REAL(KIND(1d0)),DIMENSION(24*nsh),INTENT(inout)::Tair24HR
    REAL(KIND(1d0)),INTENT(inout)                  ::qn1_av
    REAL(KIND(1d0)),INTENT(inout)                  ::dqndt!Rate of change of net radiation [W m-2 h-1] at t-1
    REAL(KIND(1d0)),INTENT(inout)                  ::qn1_s_av
    REAL(KIND(1d0)),INTENT(inout)                  ::dqnsdt !Rate of change of net radiation [W m-2 h-1] at t-1
    ! REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout)   ::qn1_store_grid
    ! REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout)   ::qn1_S_store_grid !< stored qn1 [W m-2]

    ! REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_av_store_grid
    ! REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_S_av_store_grid !< average net radiation over previous hour [W m-2]
    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(inout)::surf


    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::deltaQi ! storage heat flux of snow surfaces
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowFrac

    REAL(KIND(1d0)),DIMENSION(27),INTENT(out):: dataOutLineESTM
    REAL(KIND(1d0)),INTENT(out)::qn1_S
    REAL(KIND(1d0)),INTENT(out):: qs ! storage heat flux
    REAL(KIND(1d0)),INTENT(out):: a1 !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out):: a2 !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out):: a3 !< AnOHM coefficients of grid [W m-2]


    REAL(KIND(1d0))::HDDday ! HDDday=HDD(id-1,4) HDD at the begining of today (id-1)
    REAL(KIND(1d0))::qn1_use ! qn used in OHM calculations

    ! initialise output variables
    deltaQi=0
    snowFrac=0
    qn1_S=0
    dataOutLineESTM=-999
    qs=-999
    a1=-999
    a2=-999
    a3=-999


    ! calculate qn if qf should be included
    IF(OHMIncQF == 1) THEN
       qn1_use= qf+qn1
    ELSEIF(OHMIncQF == 0) THEN
       qn1_use= qn1
    ENDIF

    IF(StorageHeatMethod==1) THEN           !Use OHM to calculate QS
       HDDday=HDD(id-1,4)
       IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
       CALL OHM(qn1,qn1_av,dqndt,&
            qn1_S,qn1_s_av,dqnsdt,&
            tstep,dt_since_start,&
            sfr,nsurf,&
            HDDday,&
            OHM_coef,&
            OHM_threshSW,OHM_threshWD,&
            soilmoist,soilstoreCap,state,&
            BldgSurf,WaterSurf,&
            SnowUse,SnowFrac,&
            DiagQS,&
            a1,a2,a3,qs,deltaQi)

    ENDIF

    ! use AnOHM to calculate QS, TS 14 Mar 2016
    IF (StorageHeatMethod==3) THEN
       IF(Diagnose==1) WRITE(*,*) 'Calling AnOHM...'
       ! CALL AnOHM(qn1_use,qn1_store_grid,qn1_av_store_grid,qf,&
       !      MetForcingData_grid,state/surf(6,:),&
       !      alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&
       !      sfr,nsurf,nsh,EmissionsMethod,id,Gridiv,&
       !      a1,a2,a3,qs,deltaQi)
       CALL AnOHM(&
            tstep,dt_since_start,&
            qn1_use,qn1_av,dqndt,qf,&
            MetForcingData_grid,state/surf(6,:),&
            alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
            sfr,nsurf,EmissionsMethod,id,Gridiv,&
            a1,a2,a3,qs,deltaQi)! output

    END IF


    ! !Calculate QS using ESTM
    IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
       !    !CALL ESTM(QSestm,iMB)
       IF(Diagnose==1) WRITE(*,*) 'Calling ESTM...'
       CALL ESTM(&
            Gridiv,&!input
            nsh,tstep,&
            avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown,&
            bldgh,Ts5mindata_ir,&
            Tair24HR,&!inout
            dataOutLineESTM,QS)!output
       !    CALL ESTM(QSestm,Gridiv,ir)  ! iMB corrected to Gridiv, TS 09 Jun 2016
       !    QS=QSestm   ! Use ESTM qs
    ENDIF

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
       AddWater,stateOld,soilmoistOld)

    IMPLICIT NONE
    ! INTEGER,PARAMETER :: nsurf=7! number of surface types
    ! INTEGER,PARAMETER ::WaterSurf = 7
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
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out):: AddWater      !water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
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
    runoffPipes         = addPipes !Water flowing in pipes from other grids. QUESTION: No need for scaling?
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
       drain_per_tstep=0
    ENDIF

    drain(WaterSurf) = 0  ! Set drainage from water body to zero

    ! Distribute water within grid, according to WithinGridWaterDist matrix (Cols 1-7)
    IF(Diagnose==1) WRITE(*,*) 'Calling ReDistributeWater...'
    ! CALL ReDistributeWater
    !Calculates AddWater(is)
    CALL ReDistributeWater(&
         nsurf,& ! input:
         WaterSurf, snowUse, WaterDist,  sfr,   Drain,&
         AddWaterRunoff,&  ! output:
         AddWater)

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
       id,tstep,imin,it,ity,snowCalcSwitch,dayofWeek_id,CRWmin,CRWmax,&
       nsh_real,dectime,lvS_J_kg,lv_j_kg,avdens,avRh,Press_hPa,Temp_C,&
       RAsnow,psyc_hPa,avcp,sIce_hPa,&
       PervFraction,vegfraction,addimpervious,qn1_SF,qf,qs,vpd_hPa,s_hPa,&
       ResistSurf,RA,rb,tstep_real,snowdensmin,precip,PipeCapacity,RunoffToWater,&
       NonWaterFraction,wu_EveTr,wu_DecTr,wu_Grass,addVeg,addWaterBody,SnowLimPaved,SnowLimBuild,&
       SurfaceArea,FlowChange,drain,WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
       freezmelt,freezstate,freezstatevol,Qm_Melt,Qm_rain,Tsurf_ind,sfr,&
       StateLimit,AddWater,addwaterrunoff,surf,snowD,&
       runoff_per_interval,state,soilmoist,SnowPack,snowFrac,MeltWaterStore,&! inout:
       iceFrac,SnowDens,&
       snowProf,& ! output:
       runoffSnow,runoff,runoffSoil,chang,changSnow,&
       snowDepth,SnowToSurf,ev_snow,SnowRemoval,&
       evap,rss_nsurf,p_mm,rss,qe,state_per_tstep,NWstate_per_tstep,qeOut,&
       swe,ev,chSnow_per_interval,ev_per_tstep,qe_per_tstep,runoff_per_tstep,&
       surf_chang_per_tstep,runoffPipes,mwstore,runoffwaterbody,&
       runoffAGveg,runoffAGimpervious,runoffWaterBody_m3,runoffPipes_m3)

    IMPLICIT NONE

    INTEGER,INTENT(in) ::Diagnose
    INTEGER,INTENT(in) ::id
    INTEGER,INTENT(in) ::tstep
    INTEGER,INTENT(in) ::imin
    INTEGER,INTENT(in) ::it
    INTEGER,INTENT(in) ::ity !Evaporation calculated according to Rutter (1) or Shuttleworth (2)
    ! INTEGER,INTENT(in) ::snowfractionchoice

    INTEGER,DIMENSION(nsurf),INTENT(in)::snowCalcSwitch
    INTEGER,DIMENSION(3),INTENT(in)::dayofWeek_id

    REAL(KIND(1d0)),INTENT(in)::CRWmin
    REAL(KIND(1d0)),INTENT(in)::CRWmax
    REAL(KIND(1d0)),INTENT(in)::nsh_real
    REAL(KIND(1d0)),INTENT(in)::dectime
    REAL(KIND(1d0)),INTENT(in)::lvS_J_kg
    REAL(KIND(1d0)),INTENT(in)::lv_j_kg
    REAL(KIND(1d0)),INTENT(in)::avdens
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
    REAL(KIND(1d0)),INTENT(in)::RA
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
    REAL(KIND(1d0)),INTENT(in)::FlowChange

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
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::snowD
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::StateLimit !Limit for state of each surface type [mm] (specified in input files)
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::AddWater
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::addwaterrunoff
    REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in)::surf
    REAL(KIND(1d0)), DIMENSION(0:23,2),INTENT(in):: snowProf

    !Updated status: input and output
    REAL(KIND(1d0)),INTENT(inout)::runoff_per_interval! Total water transported to each grid for grid-to-grid connectivity

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::state
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::soilmoist
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowPack
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::snowFrac
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::MeltWaterStore

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::iceFrac
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowDens
    REAL(KIND(1d0)),DIMENSION(2)    ::SurplusEvap        !Surplus for evaporation in 5 min timestep


    ! output:
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoffSnow !Initialize for runoff caused by snowmelting
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoff
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoffSoil
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::chang
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::changSnow
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::snowDepth
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::SnowToSurf
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
    REAL(KIND(1d0)),INTENT(out)::runoffWaterBody_m3
    REAL(KIND(1d0)),INTENT(out)::runoffPipes_m3
    REAL(KIND(1d0)),INTENT(out)::runoffAGveg
    REAL(KIND(1d0)),INTENT(out)::runoffAGimpervious


    ! local:
    INTEGER:: is

    REAL(KIND(1d0))::surplusWaterBody
    REAL(KIND(1d0))::pin!Rain per time interval
    REAL(KIND(1d0))::sae
    REAL(KIND(1d0))::vdrc
    REAL(KIND(1d0))::sp
    REAL(KIND(1d0))::numPM
    REAL(KIND(1d0))::tlv
    REAL(KIND(1d0))::runoffAGimpervious_m3
    REAL(KIND(1d0))::runoffAGveg_m3


    tlv=lv_J_kg/tstep_real !Latent heat of vapourisation per timestep

    pin=MAX(0.,Precip)!Initiate rain data [mm]


    ! Initialize the output variables
    qe                   = 0
    ev                   = 0
    swe                  = 0
    ev_snow              = 0
    ev_per_tstep         = 0
    surf_chang_per_tstep = 0
    runoff_per_tstep     = 0
    state_per_tstep      = 0
    NWstate_per_tstep    = 0
    qeOut                = 0
    runoffwaterbody      = 0
    chSnow_per_interval  = 0
    mwstore              = 0
    runoffAGveg          = 0
    runoffAGimpervious   = 0
    surplusWaterBody     = 0
    runoffSoil           = 0
    runoff               = 0
    chang                = 0
    SurplusEvap          = 0
    SnowRemoval          = 0

    !========= these need to be wrapped================================
    sae   = s_hPa*(qn1_SF+qf-qs)    !s_haPa - slope of svp vs t curve. qn1 changed to qn1_SF, lj in May 2013
    vdrc  = vpd_hPa*avdens*avcp
    sp    = s_hPa/psyc_hPa
    numPM = sae+vdrc/RA
    !write(*,*) numPM, sae, vdrc/RA, s_hPA+psyc_hPa, NumPM/(s_hPA+psyc_hPa)
    !========= these need to be wrapped end================================

    IF(Diagnose==1) WRITE(*,*) 'Calling evap_SUEWS and SoilStore...'
    DO is=1,nsurf   !For each surface in turn
       IF (snowCalcSwitch(is)==1) THEN
          IF (sfr(is)/=0) THEN
             IF(Diagnose==1) WRITE(*,*) 'Calling SnowCalc...'
             CALL SnowCalc(&
                  id,& !input
                  tstep,imin,it,dectime,is,&
                  ity,CRWmin,CRWmax,nsh_real,lvS_J_kg,lv_j_kg,avdens,&
                  avRh,Press_hPa,Temp_C,RAsnow,psyc_hPa,avcp,sIce_hPa,&
                  PervFraction,vegfraction,addimpervious,&
                  numPM,s_hPa,ResistSurf,sp,RA,rb,tlv,snowdensmin,SnowProf,precip,&
                  PipeCapacity,RunoffToWater,runoffAGimpervious,runoffAGveg,&
                  addVeg,surplusWaterBody,SnowLimPaved,SnowLimBuild,FlowChange,drain,&
                  WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
                  freezmelt,freezstate,freezstatevol,&
                  Qm_Melt,Qm_rain,Tsurf_ind,sfr,dayofWeek_id,surf,snowD,&
                  AddWater,addwaterrunoff,&
                  SnowPack,SurplusEvap,&!inout
                  snowFrac,MeltWaterStore,iceFrac,SnowDens,&
                  runoffSnow,& ! output
                  runoff,runoffSoil,chang,changSnow,SnowToSurf,state,ev_snow,soilmoist,&
                  SnowDepth,SnowRemoval,swe,ev,chSnow_per_interval,&
                  ev_per_tstep,qe_per_tstep,runoff_per_tstep,surf_chang_per_tstep,&
                  runoffPipes,mwstore,runoffwaterbody)
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
               WetThresh(is),&!When State > WetThresh, RS=0 limit in SUEWS_evap [mm] (specified in input files)
               surf(6,is),& ! = surf(6,is), current storage capacity [mm]
               numPM,&!numerator of P-M eqn
               s_hPa,&!Vapour pressure versus temperature slope in hPa
               psyc_hPa,&!Psychometric constant in hPa
               ResistSurf,&!Surface resistance
               sp,&!Term in calculation of E
               RA,&!Aerodynamic resistance
               rb,&!Boundary layer resistance
               tlv,&!Latent heat of vaporization per timestep [J kg-1 s-1], (tlv=lv_J_kg/tstep_real)
               rss,&! output:
               ev,& ! evapotranspiration [mm]
               qe) ! latent heat flux [W m-2]


          rss_nsurf(is) = rss !Store rss for each surface

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
               drain,&!Drainage of each surface type [mm]
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
               state&!Wetness status of each surface type [mm]
               )

          evap(is) = ev !Store ev for each surface

          ! Sum evaporation from different surfaces to find total evaporation [mm]
          ev_per_tstep = ev_per_tstep+evap(is)*sfr(is)
          ! Sum change from different surfaces to find total change to surface state
          surf_chang_per_tstep = surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)
          ! Sum runoff from different surfaces to find total runoff
          runoff_per_tstep = runoff_per_tstep+runoff(is)*sfr(is)
          ! Calculate total state (including water body)
          state_per_tstep = state_per_tstep+state(is)*sfr(is)
          ! Calculate total state (excluding water body)

          IF (NonWaterFraction/=0 .AND. is/=WaterSurf) THEN
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
       QHMethod,&!input
       qn1,qf,QmRain,qeOut,qs,QmFreez,qm,avdens,avcp,tsurf,Temp_C,RA,&
       qh,qh_residual,qh_resist)!output
    IMPLICIT NONE

    INTEGER,INTENT(in) :: QHMethod ! option for QH calculation: 1, residual; 2, resistance-based

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
    REAL(KIND(1d0)),INTENT(in)::RA


    REAL(KIND(1d0)),INTENT(out)::qh
    REAL(KIND(1d0)),INTENT(out)::qh_resist
    REAL(KIND(1d0)),INTENT(out)::qh_residual

    REAL(KIND(1d0)),PARAMETER::NAN=-999

    ! Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
    qh_residual=(qn1+qf+QmRain)-(qeOut+qs+Qm+QmFreez)     !qh=(qn1+qf+QmRain+QmFreez)-(qeOut+qs+Qm)

    ! ! Calculate QH using resistance method (for testing HCW 06 Jul 2016)
    ! Aerodynamic-Resistance-based method
    IF(RA/=0) THEN
       qh_resist = avdens*avcp*(tsurf-Temp_C)/RA
    ELSE
       qh_resist=NAN
    ENDIF

    ! choose output QH
    SELECT CASE (QHMethod)
    CASE (1)
       qh= qh_residual
    CASE (2)
       qh= qh_resist
    END SELECT


  END SUBROUTINE SUEWS_cal_QH
  !========================================================================

  !===============Resistance Calculations=======================
  SUBROUTINE SUEWS_cal_Resistance(&
       StabilityMethod,&!input:
       Diagnose,AerodynamicResistanceMethod,RoughLenHeatMethod,snowUse,&
       id,it,gsModel,SMDMethod,&
       qh_obs,avdens,avcp,h_mod,qn1,dectime,zzd,z0m,zdm,&
       avU1,Temp_C,VegFraction,&
       avkdn,Kmax,G1,G2,G3,G4,G5,G6,S1,S2,TH,TL,dq,&
       xsmd,vsmd,MaxConductance,LAIMax,LAI_id,snowFrac,sfr,&
       UStar,TStar,L_mod,&!output
       zL,gsc,ResistSurf,RA,RAsnow,rb)

    IMPLICIT NONE

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
    REAL(KIND(1d0)),INTENT(in)::qn1
    REAL(KIND(1d0)),INTENT(in)::dectime    !Decimal time
    REAL(KIND(1d0)),INTENT(in)::zzd        !Active measurement height (meas. height-displac. height)
    REAL(KIND(1d0)),INTENT(in)::z0m        !Aerodynamic roughness length
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
    REAL(KIND(1d0)),INTENT(in)::vsmd       !Soil moisture deficit for vegetated surfaces only (QUESTION: what about BSoil?)

    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::MaxConductance!Max conductance [mm s-1]
    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::LAIMax        !Max LAI [m2 m-2]
    REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::LAI_id        !=LAI_id(id-1,:), LAI for each veg surface [m2 m-2]

    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::snowFrac      !Surface fraction of snow cover
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr           !Surface fractions [-]

    REAL(KIND(1d0)),INTENT(out)::TStar     !T*
    REAL(KIND(1d0)),INTENT(out)  ::UStar     !Friction velocity
    ! REAL(KIND(1d0)),INTENT(out)::psim      !Stability function of momentum
    REAL(KIND(1d0)),INTENT(out)  ::zL        !
    REAL(KIND(1d0)),INTENT(out)  ::gsc       !Surface Layer Conductance
    REAL(KIND(1d0)),INTENT(out)  ::ResistSurf!Surface resistance
    REAL(KIND(1d0)),INTENT(out)  ::RA        !Aerodynamic resistance [s m^-1]
    REAL(KIND(1d0)),INTENT(out)  ::RAsnow    !Aerodynamic resistance for snow [s m^-1]
    REAL(KIND(1d0)),INTENT(out)  ::rb        !boundary layer resistance shuttleworth
    REAL(KIND(1d0)),INTENT(out)  ::L_mod     !Obukhov length
    REAL(KIND(1d0))              ::H_init    !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity


    ! Get first estimate of sensible heat flux. Modified by HCW 26 Feb 2015
    CALL SUEWS_init_QH(&
         qh_obs,avdens,avcp,h_mod,qn1,dectime,&
         H_init)

    IF(Diagnose==1) WRITE(*,*) 'Calling STAB_lumps...'
    !u* and Obukhov length out
    CALL STAB_lumps(&
         StabilityMethod,&  ! input
         dectime,& !Decimal time
         zzd,&     !Active measurement height (meas. height-displac. height)
         z0m,&     !Aerodynamic roughness length
         zdm,&     !Displacement height
         avU1,&    !Average wind speed
         Temp_C,&  !Air temperature
         H_init,& !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
         L_mod,&! output: !Obukhov length
         TStar,& !T*
         UStar,& !Friction velocity
         zL)!Stability scale

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
    CALL SurfaceResistance(&
         id,it,&! input:
         SMDMethod,snowFrac,sfr,avkdn,Temp_C,dq,xsmd,vsmd,MaxConductance,&
         LAIMax,LAI_id,gsModel,Kmax,&
         G1,G2,G3,G4,G5,G6,TH,TL,S1,S2,&
         gsc,ResistSurf)! output:

    IF(Diagnose==1) WRITE(*,*) 'Calling BoundaryLayerResistance...'
    CALL BoundaryLayerResistance(&
         zzd,&! input:     !Active measurement height (meas. height-displac. height)
         z0m,&     !Aerodynamic roughness length
         avU1,&    !Average wind speed
         UStar,&  ! input/output:
         rb)  ! output:

  END SUBROUTINE SUEWS_cal_Resistance
  !========================================================================

  !==============Update output arrays=========================
  SUBROUTINE SUEWS_update_outputLine(&
       AdditionalWater,alb,avkdn,avU10_ms,azimuth,&!input
       chSnow_per_interval,dectime,&
       drain_per_tstep,E_mod,ev_per_tstep,ext_wu,Fc,Fc_build,fcld,&
       Fc_metab,Fc_photo,Fc_respi,Fc_traff,FlowChange,&
       h_mod,id,id_prev_t,imin,int_wu,it,iy,iy_prev_t,&
       kup,LAI,ldown,l_mod,lup,mwh,&
       MwStore,&
       nsh_real,NWstate_per_tstep,Precip,q2_gkg,&
       qeOut,qf,qh,qh_resist,Qm,QmFreez,&
       QmRain,qn1,qn1_S,qn1_SF,qs,RA,&
       resistsurf,runoffAGimpervious,runoffAGveg,&
       runoff_per_tstep,runoffPipes,runoffSoil_per_tstep,&
       runoffWaterBody,sfr,smd,smd_nsurf,SnowAlb,SnowRemoval,&
       state,state_per_tstep,surf_chang_per_tstep,swe,t2_C,&
       tot_chang_per_tstep,tsurf,UStar,wu_DecTr,&
       wu_EveTr,wu_Grass,z0m,zdm,zenith_deg,&
       datetimeLine,dataOutLineSUEWS)!output
    IMPLICIT NONE

    REAL(KIND(1d0)),PARAMETER :: NAN=-999
    INTEGER,INTENT(in) :: iy
    INTEGER,INTENT(in) :: iy_prev_t
    INTEGER,INTENT(in) :: id
    INTEGER,INTENT(in) :: id_prev_t
    INTEGER,INTENT(in) :: it
    INTEGER,INTENT(in) :: imin

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
    REAL(KIND(1d0)),INTENT(in) :: h_mod
    REAL(KIND(1d0)),INTENT(in) :: int_wu
    REAL(KIND(1d0)),INTENT(in) :: kup
    REAL(KIND(1d0)),INTENT(in) :: l_mod
    REAL(KIND(1d0)),INTENT(in) :: LAI(-4:ndays, nvegsurf)
    REAL(KIND(1d0)),INTENT(in) :: ldown
    REAL(KIND(1d0)),INTENT(in) :: lup
    REAL(KIND(1d0)),INTENT(in) :: mwh
    REAL(KIND(1d0)),INTENT(in) :: MwStore
    REAL(KIND(1d0)),INTENT(in) :: nsh_real
    REAL(KIND(1d0)),INTENT(in) :: NWstate_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: Precip
    REAL(KIND(1d0)),INTENT(in) :: q2_gkg
    REAL(KIND(1d0)),INTENT(in) :: qeOut
    REAL(KIND(1d0)),INTENT(in) :: qf
    REAL(KIND(1d0)),INTENT(in) :: qh
    REAL(KIND(1d0)),INTENT(in) :: qh_resist
    REAL(KIND(1d0)),INTENT(in) :: Qm
    REAL(KIND(1d0)),INTENT(in) :: QmFreez
    REAL(KIND(1d0)),INTENT(in) :: QmRain
    REAL(KIND(1d0)),INTENT(in) :: qn1
    REAL(KIND(1d0)),INTENT(in) :: qn1_S
    REAL(KIND(1d0)),INTENT(in) :: qn1_SF
    REAL(KIND(1d0)),INTENT(in) :: qs
    REAL(KIND(1d0)),INTENT(in) :: RA
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
    REAL(KIND(1d0)),INTENT(in) :: SnowRemoval(2)
    REAL(KIND(1d0)),INTENT(in) :: state(nsurf)
    REAL(KIND(1d0)),INTENT(in) :: state_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: surf_chang_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: swe
    REAL(KIND(1d0)),INTENT(in) :: t2_C
    REAL(KIND(1d0)),INTENT(in) :: tot_chang_per_tstep
    REAL(KIND(1d0)),INTENT(in) :: tsurf
    REAL(KIND(1d0)),INTENT(in) :: UStar
    REAL(KIND(1d0)),INTENT(in) :: wu_DecTr
    REAL(KIND(1d0)),INTENT(in) :: wu_EveTr
    REAL(KIND(1d0)),INTENT(in) :: wu_Grass
    REAL(KIND(1d0)),INTENT(in) :: z0m
    REAL(KIND(1d0)),INTENT(in) :: zdm
    REAL(KIND(1d0)),INTENT(in) :: zenith_deg


    REAL(KIND(1D0)),DIMENSION(5),INTENT(OUT)::datetimeLine
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutSUEWS-5),INTENT(out) :: dataOutLineSUEWS
    ! REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutSnow-5),INTENT(out) :: dataOutLineSnow
    ! REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutESTM-5),INTENT(out) :: dataOutLineESTM

    ! INTEGER:: is
    REAL(KIND(1d0)):: LAI_wt

    ! the variables below with '_x' endings stand for 'exported' values
    REAL(KIND(1d0))::ResistSurf_x
    REAL(KIND(1d0))::l_mod_x
    REAL(KIND(1d0))::bulkalbedo
    REAL(KIND(1d0))::smd_nsurf_x(nsurf)
    REAL(KIND(1d0))::state_x(nsurf)

    !=====================================================================
    !====================== Prepare data for output ======================
    ! values outside of reasonable range are set as NAN-like numbers. TS 10 Jun 2018

    ! Remove non-existing surface type from surface and soil outputs   ! Added back in with NANs by HCW 24 Aug 2016
    state_x=UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr<0.00001), field=state)
    smd_nsurf_x=UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr<0.00001), field=smd_nsurf)

    ResistSurf_x=MIN(9999.,ResistSurf)

    l_mod_x=MAX(MIN(9999.,l_mod), -9999.)

    ! Calculate areally-weighted LAI
    IF(iy == (iy_prev_t+1) .AND. (id-1) == 0) THEN   !Check for start of next year and avoid using LAI(id-1) as this is at the start of the year
       LAI_wt=DOT_PRODUCT(LAI(id_prev_t,:),sfr(1+2:nvegsurf+2))
    ELSE
       LAI_wt=DOT_PRODUCT(LAI(id-1,:),sfr(1+2:nvegsurf+2))
    ENDIF

    ! Calculate areally-weighted albedo
    bulkalbedo=DOT_PRODUCT(alb,sfr)

    ! NB: this part needs to be reconsidered for calculation logic. TS, 27 Sep 2018
    ! TODO: this part should be reconnected to an improved CBL interface. TS 10 Jun 2018
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



    !====================== update output line ==============================
    ! date & time:
    datetimeLine=[&
         REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),&
         REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime]
    !Define the overall output matrix to be printed out step by step
    dataOutLineSUEWS=[&
         avkdn,kup,ldown,lup,tsurf,&
         qn1,qf,qs,qh,qeOut,&
         h_mod,e_mod,qh_resist,&
         precip,ext_wu,ev_per_tstep,runoff_per_tstep,tot_chang_per_tstep,&
         surf_chang_per_tstep,state_per_tstep,NWstate_per_tstep,drain_per_tstep,smd,&
         FlowChange/nsh_real,AdditionalWater,&
         runoffSoil_per_tstep,runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,&
         int_wu,wu_EveTr,wu_DecTr,wu_Grass,&
         smd_nsurf_x(1:nsurf-1),&
         state_x(1:nsurf),&
         zenith_deg,azimuth,bulkalbedo,Fcld,&
         LAI_wt,z0m,zdm,&
         UStar,l_mod,RA,ResistSurf,&
         Fc,&
         Fc_photo,Fc_respi,Fc_metab,Fc_traff,Fc_build,&
         qn1_SF,qn1_S,SnowAlb,&
         Qm,QmFreez,QmRain,swe,mwh,MwStore,chSnow_per_interval,&
         SnowRemoval(1:2),&
         t2_C,q2_gkg,avU10_ms& ! surface-level diagonostics
         ]
    ! set invalid values to NAN
    ! dataOutLineSUEWS=set_nan(dataOutLineSUEWS)


    !====================update output line end==============================

  END SUBROUTINE SUEWS_update_outputLine
  !========================================================================

  !==============Update output arrays=========================
  SUBROUTINE SUEWS_update_output(&
       SnowUse,storageheatmethod,&!input
       ReadLinesMetdata,NumberOfGrids,&
       ir,gridiv,datetimeLine,dataOutLineSUEWS,dataOutLineSnow,dataOutLineESTM,&!input
       dataOutSUEWS,dataOutSnow,dataOutESTM)!inout
    IMPLICIT NONE

    INTEGER,INTENT(in) ::ReadLinesMetdata
    INTEGER,INTENT(in) ::NumberOfGrids
    INTEGER,INTENT(in) ::Gridiv
    INTEGER,INTENT(in) ::SnowUse
    INTEGER,INTENT(in) ::storageheatmethod
    INTEGER,INTENT(in) ::ir

    REAL(KIND(1d0)),DIMENSION(5),INTENT(in) :: datetimeLine
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutSUEWS-5),INTENT(in) :: dataOutLineSUEWS
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutESTM-5),INTENT(in) :: dataOutLineESTM
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutSnow-5),INTENT(in) :: dataOutLineSnow


    REAL(KIND(1d0)),INTENT(inout) :: dataOutSUEWS(ReadLinesMetdata,ncolumnsDataOutSUEWS,NumberOfGrids)
    REAL(KIND(1d0)),INTENT(inout) :: dataOutSnow(ReadLinesMetdata,ncolumnsDataOutSnow,NumberOfGrids)
    REAL(KIND(1d0)),INTENT(inout) :: dataOutESTM(ReadLinesMetdata,ncolumnsDataOutESTM,NumberOfGrids)


    !====================== update output arrays ==============================
    !Define the overall output matrix to be printed out step by step
    dataOutSUEWS(ir,1:ncolumnsDataOutSUEWS,Gridiv)=[datetimeLine,set_nan(dataOutLineSUEWS)]
    ! ! set invalid values to NAN
    ! dataOutSUEWS(ir,6:ncolumnsDataOutSUEWS,Gridiv)=set_nan(dataOutSUEWS(ir,6:ncolumnsDataOutSUEWS,Gridiv))

    IF (snowUse==1) THEN
       dataOutSnow(ir,1:ncolumnsDataOutSnow,Gridiv)=[datetimeLine,set_nan(dataOutLineSnow)]
    END IF

    IF (storageheatmethod==4) THEN
       dataOutESTM(ir,1:ncolumnsDataOutESTM,Gridiv)=[datetimeLine,set_nan(dataOutLineESTM)]
    END IF

    !====================update output arrays end==============================

  END SUBROUTINE SUEWS_update_output
  !========================================================================


  SUBROUTINE SUEWS_cal_Diagnostics(&
       dectime,&!input
       avU1,Temp_C,&
       tsurf,qh,&
       Press_hPa,qe,&
       veg_fr,z0m,avdens,avcp,lv_J_kg,tstep_real,&
       RoughLenHeatMethod,StabilityMethod,&
       avU10_ms,t2_C,q2_gkg,L_MOD)!output
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(in) ::dectime
    REAL(KIND(1d0)),INTENT(in) ::avU1,Temp_C
    REAL(KIND(1d0)),INTENT(in) ::tsurf,qh
    REAL(KIND(1d0)),INTENT(in) ::Press_hPa,qe
    REAL(KIND(1d0)),INTENT(in) :: veg_fr,z0m,avdens,avcp,lv_J_kg,tstep_real

    ! INTEGER,INTENT(in)         :: opt ! 0 for momentum, 1 for temperature, 2 for humidity
    INTEGER,INTENT(in)         :: RoughLenHeatMethod,StabilityMethod

    REAL(KIND(1d0)),INTENT(out):: avU10_ms,t2_C,q2_gkg,L_MOD
    REAL(KIND(1d0))::tlv,z2zd,zdm,H_init,TStar,zL,UStar
    REAL(KIND(1d0)),PARAMETER::k=0.4

    tlv=lv_J_kg/tstep_real !Latent heat of vapourisation per timestep
    z2zd=2 ! height at 2m assuming Displacement height is ZERO
    zdm=0 ! assuming Displacement height is ZERO

    ! get !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
    CALL SUEWS_init_QH(&
         qh,avdens,avcp,qh,0d0,dectime,& ! use qh as qh_obs to initialise H_init
         H_init)

    ! redo the calculation for stability correction
    CALL STAB_lumps(&
                                ! input
         StabilityMethod,&
         dectime,& !Decimal time
         z2zd,&     !Active measurement height (meas. height-displac. height)
         z0m,&     !Aerodynamic roughness length
         zdm,&     !Displacement height
         avU1,&    !Average wind speed
         Temp_C,&  !Air temperature
         h_init,    & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
                                ! output:
         L_MOD,& !Obukhov length
         TStar,& !T*
         UStar,& !Friction velocity
         zL)!Stability scale


    ! wind speed:
    CALL diagSfc(0d0,0d0,UStar,veg_fr,z0m,L_mod,k,avdens,avcp,tlv,avU10_ms,0,RoughLenHeatMethod,StabilityMethod)
    ! temperature:
    CALL diagSfc(tsurf,qh,UStar,veg_fr,z0m,L_mod,k,avdens,avcp,tlv,t2_C,1,RoughLenHeatMethod,StabilityMethod)
    ! humidity:
    CALL diagSfc(qsatf(tsurf,Press_hPa)*1000,& ! Saturation specific humidity at surface in g/kg
         qe,UStar,veg_fr,z0m,L_mod,k,avdens,avcp,tlv,q2_gkg,2,RoughLenHeatMethod,StabilityMethod)

  END SUBROUTINE SUEWS_cal_Diagnostics


  ! Calculate tstep-derived variables
  SUBROUTINE SUEWS_cal_tstep(&
       tstep,& ! input
       nsh, nsh_real, tstep_real) ! output
    IMPLICIT NONE
    INTEGER,INTENT(in)::tstep ! number of timesteps per hour
    ! values that are derived from tstep
    INTEGER,INTENT(out)::nsh ! number of timesteps per hour
    REAL(KIND(1D0)),INTENT(out)::nsh_real ! nsh in type real
    REAL(KIND(1D0)),INTENT(out)::tstep_real ! tstep in type real
    nsh=3600/tstep
    nsh_real=nsh*1.0
    tstep_real=tstep*1.0

  END SUBROUTINE SUEWS_cal_tstep

  SUBROUTINE SUEWS_cal_surf(&
       sfr,& !input
       vegfraction,ImpervFraction,PervFraction,NonWaterFraction) ! output
    IMPLICIT NONE

    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)::sfr
    REAL(KIND(1D0)),INTENT(OUT)::VegFraction
    REAL(KIND(1D0)),INTENT(OUT)::ImpervFraction
    REAL(KIND(1D0)),INTENT(OUT)::PervFraction
    REAL(KIND(1D0)),INTENT(OUT)::NonWaterFraction


    VegFraction=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)
    ImpervFraction=sfr(PavSurf)+sfr(BldgSurf)
    PervFraction=1-ImpervFraction
    NonWaterFraction=1 - sfr(WaterSurf)

  END SUBROUTINE SUEWS_cal_surf


  SUBROUTINE SUEWS_cal_weekday(&
       iy,id,lat,& !input
       dayofWeek_id) !output
    IMPLICIT NONE

    INTEGER,INTENT(in) :: iy  ! year
    INTEGER,INTENT(in) :: id  ! day of year
    REAL(KIND(1d0)),INTENT(in):: lat

    INTEGER,DIMENSION(3),INTENT(OUT) ::dayofWeek_id

    INTEGER::wd
    INTEGER::mb
    INTEGER::date
    INTEGER::seas



    CALL day2month(id,mb,date,seas,iy,lat) !Calculate real date from doy
    CALL Day_of_Week(date,mb,iy,wd)        !Calculate weekday (1=Sun, ..., 7=Sat)

    dayofWeek_id(1)=wd      !Day of week
    dayofWeek_id(2)=mb      !Month
    dayofweek_id(3)=seas    !Season

  END SUBROUTINE SUEWS_cal_weekday


  SUBROUTINE SUEWS_cal_DLS(&
       id,startDLS,endDLS,& !input
       DLS) !output
    IMPLICIT NONE

    INTEGER, INTENT(in) :: id,startDLS,endDLS
    INTEGER, INTENT(out) :: DLS

    DLS=0
    IF ( id>startDLS .AND. id<endDLS ) dls=1

  END SUBROUTINE SUEWS_cal_DLS

  SUBROUTINE diagSfc(&
       xSurf,xFlux,us,VegFraction,z0m,L_mod,k,avdens,avcp,tlv,&
       xDiag,opt,RoughLenHeatMethod,StabilityMethod)
    ! TS 05 Sep 2017: improved interface
    ! TS 20 May 2017: calculate surface-level diagonostics


    IMPLICIT NONE

    REAL(KIND(1d0)),INTENT(in) :: xSurf,xFlux,us,VegFraction,z0m,L_mod,k,avdens,avcp,tlv
    REAL(KIND(1d0)),INTENT(out):: xDiag
    INTEGER,INTENT(in)         :: opt ! 0 for momentum, 1 for temperature, 2 for humidity
    INTEGER,INTENT(in)         :: RoughLenHeatMethod,StabilityMethod

    REAL(KIND(1d0))            :: &
         psymz2,psymz10,psymz0,psyhz2,psyhz0,& ! stability correction functions
         z0h,& ! Roughness length for heat
         z2zd,z10zd!stability correction functions
    REAL(KIND(1d0)),PARAMETER :: muu=1.46e-5 !molecular viscosity
    REAL(KIND(1d0)),PARAMETER :: nan=-999



    !***************************************************************
    ! log-law based stability corrections:
    ! Roughness length for heat
    IF (RoughLenHeatMethod==1) THEN !Brutasert (1982) z0h=z0/10(see Grimmond & Oke, 1986)
       z0h=z0m/10
    ELSEIF (RoughLenHeatMethod==2) THEN ! Kawai et al. (2007)
       !z0h=z0m*exp(2-(1.2-0.9*veg_fr**0.29)*(us*z0m/muu)**0.25)
       ! Changed by HCW 05 Nov 2015 (veg_fr includes water; VegFraction = veg + bare soil)
       z0h=z0m*EXP(2-(1.2-0.9*VegFraction**0.29)*(us*z0m/muu)**0.25)
    ELSEIF (RoughLenHeatMethod==3) THEN
       z0h=z0m*EXP(-20.) ! Voogt and Grimmond, JAM, 2000
    ELSEIF (RoughLenHeatMethod==4) THEN
       z0h=z0m*EXP(2-1.29*(us*z0m/muu)**0.25) !See !Kanda and Moriwaki (2007),Loridan et al. (2010)
    ENDIF

    ! z0h=z0m/5

    ! zX-z0
    z2zd=2!+z0h   ! set lower limit as z0h to prevent arithmetic error, zd=0
    z10zd=10!+z0m ! set lower limit as z0m to prevent arithmetic error, zd=0

    ! stability correction functions
    ! momentum:
    psymz10=stab_fn_mom(StabilityMethod,z10zd/L_mod,z10zd/L_mod)
    psymz2=stab_fn_mom(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
    psymz0=stab_fn_mom(StabilityMethod,z0m/L_mod,z0m/L_mod)

    ! heat and vapor: assuming both are the same
    psyhz2=stab_fn_heat(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
    psyhz0=stab_fn_heat(StabilityMethod,z0h/L_mod,z0h/L_mod)
    !***************************************************************
    IF ( xSurf==nan ) THEN
       ! xSurf can be nan e.g. when TSurf is not calculated
       ! if so xDiag is set as nan as well
       xDiag=nan
    ELSE
       SELECT CASE (opt)
       CASE (0) ! wind (momentum) at 10 m
          xDiag=us/k*(LOG(z10zd/z0m)-psymz10+psymz0) ! Brutsaert (2005), p51, eq.2.54

       CASE (1) ! temperature at 2 m
          xDiag=xSurf-xFlux/(k*us*avdens*avcp)*(LOG(z2zd/z0h)-psyhz2+psyhz0) ! Brutsaert (2005), p51, eq.2.55
          !  IF ( ABS((LOG(z2zd/z0h)-psyhz2+psyhz0))>10 ) THEN
          !     PRINT*, '#####################################'
          !     PRINT*, 'xSurf',xSurf
          !     PRINT*, 'xFlux',xFlux
          !     PRINT*, 'k*us*avdens*avcp',k*us*avdens*avcp
          !     PRINT*, 'k',k
          !     PRINT*, 'us',us
          !     PRINT*, 'avdens',avdens
          !     PRINT*, 'avcp',avcp
          !     PRINT*, 'xFlux/X',xFlux/(k*us*avdens*avcp)
          !     PRINT*, 'stab',(LOG(z2zd/z0h)-psyhz2+psyhz0)
          !     PRINT*, 'LOG(z2zd/z0h)',LOG(z2zd/z0h)
          !     PRINT*, 'z2zd',z2zd,'L_mod',L_mod,'z0h',z0h
          !     PRINT*, 'z2zd/L_mod',z2zd/L_mod
          !     PRINT*, 'psyhz2',psyhz2
          !     PRINT*, 'psyhz0',psyhz0
          !     PRINT*, 'psyhz2-psyhz0',psyhz2-psyhz0
          !     PRINT*, 'xDiag',xDiag
          !     PRINT*, '*************************************'
          !  END IF


       CASE (2) ! humidity at 2 m
          xDiag=xSurf-xFlux/(k*us*avdens*tlv)*(LOG(z2zd/z0h)-psyhz2+psyhz0) ! Brutsaert (2005), p51, eq.2.56

       END SELECT


    END IF

  END SUBROUTINE diagSfc



  !===============set variable of invalid value to NAN=====================
  ELEMENTAL FUNCTION set_nan(x) RESULT(xx)
    IMPLICIT NONE
    REAL(KIND(1d0)),PARAMETER::pNAN=30000 ! 30000 to prevent water_state being filtered out as it can be large
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

  !===============the functions below are only for test in f2py conversion===
  FUNCTION square(x) RESULT(xx)
    IMPLICIT NONE
    REAL(KIND(1d0)),PARAMETER::pNAN=9999
    REAL(KIND(1d0)),PARAMETER::NAN=-999
    REAL(KIND(1d0)),INTENT(in)::x
    REAL(KIND(1d0))::xx

    xx=x**2+nan/pNAN
    xx=x**2

  END FUNCTION square

  FUNCTION square_real(x) RESULT(xx)
    IMPLICIT NONE
    REAL,PARAMETER::pNAN=9999
    REAL,PARAMETER::NAN=-999
    REAL,INTENT(in)::x
    REAL::xx

    xx=x**2+nan/pNAN
    xx=x**2

  END FUNCTION square_real

  SUBROUTINE output_name_n(i,name,group,aggreg)
    ! used by f2py module `SuPy` to handle output names
    IMPLICIT NONE
    ! the dimension is potentially incorrect,
    ! which should be consistent with that in output module
    INTEGER,INTENT(in) :: i
    CHARACTER(len = 15),INTENT(out) :: name,group,aggreg

    INTEGER :: n
    n=SIZE(varList, dim=1)
    IF ( i<n .AND.i>0  ) THEN
       name   = TRIM(varList(i)%header)
       group  = TRIM(varList(i)%group)
       aggreg = TRIM(varList(i)%aggreg)
    ELSE
       name   = ''
       group  = ''
       aggreg = ''
    END IF


    ! DO i = 1, SIZE(varList, dim=1), 1
    !    names(i)=TRIM(varList(i)%header)
    !    ! names(i,1)=trim(varList(i)%header)
    !    ! names(i,2)=trim(varList(i)%unit)
    !    ! ! names(i,3)=trim(varList(i)%fmt)
    !    ! names(i,3)=trim(varList(i)%longNm)
    !    ! names(i,4)=trim(varList(i)%aggreg)
    !    ! names(i,5)=trim(varList(i)%group)
    !    ! print*, names(i,:)
    ! END DO
    ! print*, varList

  END SUBROUTINE output_name_n


  SUBROUTINE output_size(n)
    ! used by f2py module `SuPy` to get size of the output list
    IMPLICIT NONE
    ! the dimension is potentially incorrect,
    ! which should be consistent with that in output module
    INTEGER,INTENT(out) :: n


    n=SIZE(varList, dim=1)

  END SUBROUTINE output_size

  ! SUBROUTINE output_names(names)
  !   ! used by f2py module `SuPy` to handle output names
  !   IMPLICIT NONE
  !   ! the dimension is potentially incorrect,
  !   ! which should be consistent with that in output module
  !   ! CHARACTER(len = 50),DIMENSION(300,5),INTENT(out) :: names
  !   ! CHARACTER(len = *),DIMENSION(300),INTENT(out) :: names
  !   CHARACTER(len = 12),INTENT(out) :: names(300)
  !
  !   INTEGER :: i,n,err,stat
  !   ! name=TRIM(varList(1)%header)
  !   n=SIZE(varList, dim=1)
  !   ! ALLOCATE(names(10,n), stat=err)
  !   ! IF (stat /= 0) PRINT *, "names: Allocation request denied"
  !
  !
  !   DO i = 1, n
  !      names(i)=TRIM(varList(i)%header)
  !      ! names(i,1)=trim(varList(i)%header)
  !      ! names(i,2)=trim(varList(i)%unit)
  !      ! ! names(i,3)=trim(varList(i)%fmt)
  !      ! names(i,3)=trim(varList(i)%longNm)
  !      ! names(i,4)=trim(varList(i)%aggreg)
  !      ! names(i,5)=trim(varList(i)%group)
  !      ! print*, names(i,:)
  !   END DO
  !   ! print*, varList
  !   ! IF (ALLOCATED(names)) DEALLOCATE(names)
  !   ! IF (stat /= 0) PRINT *, "names: Deallocation request denied"
  !
  ! END SUBROUTINE output_names

  SUBROUTINE SUEWS_write_model_state(&
       AerodynamicResistanceMethod,AH_MIN,AHProf_tstep,AH_SLOPE_Cooling,& ! input&inout in alphabetical order
       AH_SLOPE_Heating,alb,albDecTr,albEveTr,albGrass,AlbMax_DecTr,&
       AlbMax_EveTr,AlbMax_Grass,AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
       alpha_bioCO2,alpha_enh_bioCO2,alt,avkdn,avRh,avU1,BaseT,BaseTe,&
       BaseTHDD,beta_bioCO2,beta_enh_bioCO2,bldgH,CapMax_dec,CapMin_dec,&
       chAnOHM,cpAnOHM,CRWmax,CRWmin,DayWat,DayWatPer,&
       DecidCap,dectime,DecTreeH,Diagnose,DiagQN,DiagQS,DRAINRT,&
       EF_umolCO2perJ,emis,EmissionsMethod,EnEF_v_Jkm,endDLS,EveTreeH,FAIBldg,&
       FAIDecTree,FAIEveTree,Faut,FcEF_v_kgkm,fcld_obs,FlowChange,&
       FrFossilFuel_Heat,FrFossilFuel_NonHeat,G1,G2,G3,G4,G5,G6,GDD,&
       GDDFull,Gridiv,gsModel,HDD,HumActivity_tstep,&
       IceFrac,id,id_prev_t,Ie_a,Ie_end,Ie_m,Ie_start,imin,&
       InternalWaterUse_h,IrrFracConif,IrrFracDecid,IrrFracGrass,it,ity,&
       iy,iy_prev_t,kkAnOHM,Kmax,LAI,LAICalcYes,LAIMax,LAIMin,LAI_obs,&
       LAIPower,LAIType,lat,ldown_obs,lng,MaxConductance,MaxQFMetab,&
       MeltWaterStore,MetForcingData_grid,MinQFMetab,min_res_bioCO2,&
       NARP_EMIS_SNOW,NARP_TRANS_SITE,NetRadiationMethod,&
       NumCapita,OHM_coef,OHMIncQF,OHM_threshSW,&
       OHM_threshWD,PipeCapacity,PopDensDaytime,&
       PopDensNighttime,PopProf_tstep,PorMax_dec,PorMin_dec,porosity,&
       Precip,PrecipLimit,PrecipLimitAlb,Press_hPa,QF0_BEU,Qf_A,Qf_B,&
       Qf_C,qh_obs,qn1_av_store_grid,qn1_obs,qn1_S_av_store_grid,qn1_S_store_grid,&
       qn1_store_grid,RadMeltFact,RAINCOVER,RainMaxRes,resp_a,resp_b,&
       RoughLenHeatMethod,RoughLenMomMethod,RunoffToWater,S1,S2,&
       SatHydraulicConduct,SDDFull,sfr,SMDMethod,SnowAlb,SnowAlbMax,&
       SnowAlbMin,snowD,SnowDens,SnowDensMax,SnowDensMin,SnowfallCum,snowFrac,&
       SnowLimBuild,SnowLimPaved,snow_obs,SnowPack,SnowProf,snowUse,SoilDepth,&
       soilmoist,soilstoreCap,StabilityMethod,startDLS,state,StateLimit,&
       StorageHeatMethod,surf,SurfaceArea,Tair24HR,tau_a,tau_f,tau_r,&
       T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,TempMeltFact,TH,&
       theta_bioCO2,timezone,TL,TrafficRate,TrafficUnits,&
       TraffProf_tstep,Ts5mindata_ir,tstep,veg_type,&
       WaterDist,WaterUseMethod,WetThresh,WU_Day,WUProfA_tstep,&
       WUProfM_tstep,xsmd,Z)

    USE data_in,ONLY:FileOutputPath

    IMPLICIT NONE

    INTEGER,INTENT(IN)::AerodynamicResistanceMethod
    INTEGER,INTENT(IN)::Diagnose
    INTEGER,INTENT(IN)::DiagQN
    INTEGER,INTENT(IN)::DiagQS
    INTEGER,INTENT(IN)::startDLS
    INTEGER,INTENT(IN)::endDLS
    INTEGER,INTENT(IN)::EmissionsMethod
    INTEGER,INTENT(IN)::Gridiv
    INTEGER,INTENT(IN)::gsModel
    INTEGER,INTENT(IN)::id
    INTEGER,INTENT(IN)::id_prev_t
    INTEGER,INTENT(IN)::Ie_end
    INTEGER,INTENT(IN)::Ie_start
    INTEGER,INTENT(IN)::imin
    INTEGER,INTENT(IN)::it
    INTEGER,INTENT(IN)::ity
    INTEGER,INTENT(IN)::iy
    INTEGER,INTENT(IN)::iy_prev_t
    INTEGER,INTENT(IN)::LAICalcYes
    INTEGER,INTENT(IN)::NetRadiationMethod
    INTEGER,INTENT(IN)::OHMIncQF
    INTEGER,INTENT(IN)::RoughLenHeatMethod
    INTEGER,INTENT(IN)::RoughLenMomMethod
    INTEGER,INTENT(IN)::SMDMethod
    INTEGER,INTENT(IN)::snowUse
    INTEGER,INTENT(IN)::StabilityMethod
    INTEGER,INTENT(IN)::StorageHeatMethod
    INTEGER,INTENT(IN)::tstep
    INTEGER,INTENT(IN)::veg_type
    INTEGER,INTENT(IN)::WaterUseMethod

    REAL(KIND(1D0)),INTENT(IN)::AlbMax_DecTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMax_EveTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMax_Grass
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_DecTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_EveTr
    REAL(KIND(1D0)),INTENT(IN)::AlbMin_Grass
    REAL(KIND(1D0)),INTENT(IN)::alt
    REAL(KIND(1D0)),INTENT(IN)::avkdn
    REAL(KIND(1D0)),INTENT(IN)::avRh
    REAL(KIND(1D0)),INTENT(IN)::avU1
    REAL(KIND(1D0)),INTENT(IN)::BaseTHDD
    REAL(KIND(1D0)),INTENT(IN)::bldgH
    REAL(KIND(1D0)),INTENT(IN)::CapMax_dec
    REAL(KIND(1D0)),INTENT(IN)::CapMin_dec
    REAL(KIND(1D0)),INTENT(IN)::CRWmax
    REAL(KIND(1D0)),INTENT(IN)::CRWmin
    REAL(KIND(1D0)),INTENT(IN)::dectime
    REAL(KIND(1D0)),INTENT(IN)::DecTreeH
    REAL(KIND(1D0)),INTENT(IN)::DRAINRT
    REAL(KIND(1D0)),INTENT(IN)::EF_umolCO2perJ
    REAL(KIND(1D0)),INTENT(IN)::EnEF_v_Jkm
    REAL(KIND(1D0)),INTENT(IN)::EveTreeH
    REAL(KIND(1D0)),INTENT(IN)::FAIBldg
    REAL(KIND(1D0)),INTENT(IN)::FAIDecTree
    REAL(KIND(1D0)),INTENT(IN)::FAIEveTree
    REAL(KIND(1D0)),INTENT(IN)::Faut
    REAL(KIND(1D0)),INTENT(IN)::FcEF_v_kgkm
    REAL(KIND(1D0)),INTENT(IN)::fcld_obs
    REAL(KIND(1D0)),INTENT(IN)::FlowChange
    REAL(KIND(1D0)),INTENT(IN)::FrFossilFuel_Heat
    REAL(KIND(1D0)),INTENT(IN)::FrFossilFuel_NonHeat
    REAL(KIND(1D0)),INTENT(IN)::G1
    REAL(KIND(1D0)),INTENT(IN)::G2
    REAL(KIND(1D0)),INTENT(IN)::G3
    REAL(KIND(1D0)),INTENT(IN)::G4
    REAL(KIND(1D0)),INTENT(IN)::G5
    REAL(KIND(1D0)),INTENT(IN)::G6
    REAL(KIND(1D0)),INTENT(IN)::InternalWaterUse_h
    REAL(KIND(1D0)),INTENT(IN)::IrrFracConif
    REAL(KIND(1D0)),INTENT(IN)::IrrFracDecid
    REAL(KIND(1D0)),INTENT(IN)::IrrFracGrass
    REAL(KIND(1D0)),INTENT(IN)::Kmax
    REAL(KIND(1D0)),INTENT(IN)::LAI_obs
    REAL(KIND(1D0)),INTENT(IN)::lat
    REAL(KIND(1D0)),INTENT(IN)::ldown_obs
    REAL(KIND(1D0)),INTENT(IN)::lng
    REAL(KIND(1D0)),INTENT(IN)::MaxQFMetab
    REAL(KIND(1D0)),INTENT(IN)::MinQFMetab
    REAL(KIND(1D0)),INTENT(IN)::NARP_EMIS_SNOW
    REAL(KIND(1D0)),INTENT(IN)::NARP_TRANS_SITE
    REAL(KIND(1D0)),INTENT(IN)::NumCapita
    REAL(KIND(1D0)),INTENT(IN)::PipeCapacity
    REAL(KIND(1D0)),INTENT(IN)::PopDensDaytime
    REAL(KIND(1D0)),INTENT(IN)::PopDensNighttime
    REAL(KIND(1D0)),INTENT(IN)::PorMax_dec
    REAL(KIND(1D0)),INTENT(IN)::PorMin_dec
    REAL(KIND(1D0)),INTENT(IN)::Precip
    REAL(KIND(1D0)),INTENT(IN)::PrecipLimit
    REAL(KIND(1D0)),INTENT(IN)::PrecipLimitAlb
    REAL(KIND(1D0)),INTENT(IN)::Press_hPa
    REAL(KIND(1D0)),INTENT(IN)::qh_obs
    REAL(KIND(1D0)),INTENT(IN)::qn1_obs
    REAL(KIND(1D0)),INTENT(IN)::RadMeltFact
    REAL(KIND(1D0)),INTENT(IN)::RAINCOVER
    REAL(KIND(1D0)),INTENT(IN)::RainMaxRes
    REAL(KIND(1D0)),INTENT(IN)::RunoffToWater
    REAL(KIND(1D0)),INTENT(IN)::S1
    REAL(KIND(1D0)),INTENT(IN)::S2
    REAL(KIND(1D0)),INTENT(IN)::SnowAlbMax
    REAL(KIND(1D0)),INTENT(IN)::SnowAlbMin
    REAL(KIND(1D0)),INTENT(IN)::SnowDensMax
    REAL(KIND(1D0)),INTENT(IN)::SnowDensMin
    REAL(KIND(1D0)),INTENT(IN)::SnowLimBuild
    REAL(KIND(1D0)),INTENT(IN)::SnowLimPaved
    REAL(KIND(1D0)),INTENT(IN)::snow_obs
    REAL(KIND(1D0)),INTENT(IN)::SurfaceArea
    REAL(KIND(1D0)),INTENT(IN)::tau_a
    REAL(KIND(1D0)),INTENT(IN)::tau_f
    REAL(KIND(1D0)),INTENT(IN)::tau_r
    REAL(KIND(1D0)),INTENT(IN)::Temp_C
    REAL(KIND(1D0)),INTENT(IN)::TempMeltFact
    REAL(KIND(1D0)),INTENT(IN)::TH
    REAL(KIND(1D0)),INTENT(IN)::timezone
    REAL(KIND(1D0)),INTENT(IN)::TL
    REAL(KIND(1D0)),INTENT(IN)::TrafficUnits
    REAL(KIND(1D0)),INTENT(IN)::xsmd
    REAL(KIND(1D0)),INTENT(IN)::Z

    INTEGER,DIMENSION(NVEGSURF),INTENT(IN)::LAIType

    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::AH_MIN
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::AH_SLOPE_Cooling
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::AH_SLOPE_Heating
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::QF0_BEU
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::Qf_A
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::Qf_B
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::Qf_C
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::T_CRITIC_Cooling
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::T_CRITIC_Heating
    REAL(KIND(1D0)),DIMENSION(2),INTENT(IN)                   ::TrafficRate
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)                   ::Ie_a
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)                   ::Ie_m
    REAL(KIND(1D0)),DIMENSION(3),INTENT(IN)                   ::MaxConductance
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::AHProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::HumActivity_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::PopProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::TraffProf_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::WUProfA_tstep
    REAL(KIND(1D0)),DIMENSION(24*3600/tstep,2),INTENT(IN)     ::WUProfM_tstep
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN)                   ::DayWat
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN)                   ::DayWatPer
    REAL(KIND(1D0)),DIMENSION(nsurf+1),INTENT(IN)             ::OHM_threshSW
    REAL(KIND(1D0)),DIMENSION(nsurf+1),INTENT(IN)             ::OHM_threshWD
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::chAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::cpAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::emis
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::kkAnOHM
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::SatHydraulicConduct
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::sfr
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::snowD
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::SoilDepth
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::soilstoreCap
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::StateLimit
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(IN)               ::WetThresh
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::alpha_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::alpha_enh_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::BaseT
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::BaseTe
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::beta_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::beta_enh_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::GDDFull
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::LAIMax
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::LAIMin
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::min_res_bioCO2
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::resp_a
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::resp_b
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::SDDFull
    REAL(KIND(1D0)),DIMENSION(NVEGSURF),INTENT(IN)            ::theta_bioCO2
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in)                   ::Ts5mindata_ir
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in)              ::snowProf
    REAL(KIND(1D0)),DIMENSION(NSURF+1,NSURF-1),INTENT(IN)     ::WaterDist
    REAL(KIND(1D0)),DIMENSION(nsurf+1,4,3),INTENT(IN)         ::OHM_coef
    REAL(KIND(1D0)),DIMENSION(4,NVEGSURF),INTENT(IN)          ::LAIPower
    REAL(KIND(1D0)),DIMENSION(:,:),INTENT(IN)                 ::MetForcingData_grid

    REAL(KIND(1D0)),INTENT(INOUT)                             ::SnowfallCum
    REAL(KIND(1D0)),INTENT(INOUT)                             ::SnowAlb
    REAL(KIND(1d0)),DIMENSION(24*3600/tstep),INTENT(inout)    ::Tair24HR
    REAL(KIND(1D0)),DIMENSION(2*3600/tstep+1),INTENT(INOUT)   ::qn1_av_store_grid
    REAL(KIND(1D0)),DIMENSION(2*3600/tstep+1),INTENT(INOUT)   ::qn1_S_av_store_grid
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albDecTr
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albEveTr
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::albGrass
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::DecidCap
    REAL(KIND(1D0)),DIMENSION(0:NDAYS),INTENT(INOUT)          ::porosity
    REAL(KIND(1D0)),DIMENSION(0:NDAYS,5),INTENT(INOUT)        ::GDD
    REAL(KIND(1D0)),DIMENSION(0:NDAYS,9),INTENT(INOUT)        ::WU_Day
    REAL(KIND(1D0)),DIMENSION(6,NSURF),INTENT(INOUT)          ::surf
    REAL(KIND(1D0)),DIMENSION(-4:NDAYS,6),INTENT(INOUT)       ::HDD
    REAL(KIND(1D0)),DIMENSION(-4:NDAYS,NVEGSURF),INTENT(INOUT)::LAI
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::alb
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::IceFrac
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::MeltWaterStore
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::SnowDens
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::snowFrac
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::SnowPack
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::soilmoist
    REAL(KIND(1D0)),DIMENSION(NSURF),INTENT(INOUT)            ::state
    REAL(KIND(1D0)),DIMENSION(3600/tstep),INTENT(INOUT)       ::qn1_S_store_grid
    REAL(KIND(1D0)),DIMENSION(3600/tstep),INTENT(INOUT)       ::qn1_store_grid


    INTEGER :: fn,xx=0
    CHARACTER(len=1000) :: FileOut
    CHARACTER(len=20)::str_grid

    WRITE(str_grid,'(i4)') Gridiv
    FileOut=TRIM(FileOutputPath)//&
         'init_cond_'//&
         TRIM(ADJUSTL(str_grid))//&
         '.txt'
    ! write out data
    IF ( xx==0 ) THEN
       fn=50
       OPEN(fn,file=TRIM(FileOut),status='unknown')!,err=112)

       WRITE(fn,*)'AerodynamicResistanceMethod',AerodynamicResistanceMethod
       WRITE(fn,*)'Diagnose',Diagnose
       WRITE(fn,*)'DiagQN',DiagQN
       WRITE(fn,*)'DiagQS',DiagQS
       WRITE(fn,*)'startDLS',startDLS
       WRITE(fn,*)'endDLS',endDLS
       WRITE(fn,*)'EmissionsMethod',EmissionsMethod
       WRITE(fn,*)'Gridiv',Gridiv
       WRITE(fn,*)'gsModel',gsModel
       WRITE(fn,*)'id',id
       WRITE(fn,*)'id_prev_t',id_prev_t
       WRITE(fn,*)'Ie_end',Ie_end
       WRITE(fn,*)'Ie_start',Ie_start
       WRITE(fn,*)'imin',imin
       WRITE(fn,*)'it',it
       WRITE(fn,*)'ity',ity
       WRITE(fn,*)'iy',iy
       WRITE(fn,*)'iy_prev_t',iy_prev_t
       WRITE(fn,*)'LAICalcYes',LAICalcYes
       WRITE(fn,*)'NetRadiationMethod',NetRadiationMethod
       WRITE(fn,*)'OHMIncQF',OHMIncQF
       WRITE(fn,*)'RoughLenHeatMethod',RoughLenHeatMethod
       WRITE(fn,*)'RoughLenMomMethod',RoughLenMomMethod
       WRITE(fn,*)'SMDMethod',SMDMethod
       WRITE(fn,*)'snowUse',snowUse
       WRITE(fn,*)'StabilityMethod',StabilityMethod
       WRITE(fn,*)'StorageHeatMethod',StorageHeatMethod
       WRITE(fn,*)'tstep',tstep
       WRITE(fn,*)'veg_type',veg_type
       WRITE(fn,*)'WaterUseMethod',WaterUseMethod
       WRITE(fn,*)'AlbMax_DecTr',AlbMax_DecTr
       WRITE(fn,*)'AlbMax_EveTr',AlbMax_EveTr
       WRITE(fn,*)'AlbMax_Grass',AlbMax_Grass
       WRITE(fn,*)'AlbMin_DecTr',AlbMin_DecTr
       WRITE(fn,*)'AlbMin_EveTr',AlbMin_EveTr
       WRITE(fn,*)'AlbMin_Grass',AlbMin_Grass
       WRITE(fn,*)'alt',alt
       WRITE(fn,*)'avkdn',avkdn
       WRITE(fn,*)'avRh',avRh
       WRITE(fn,*)'avU1',avU1
       WRITE(fn,*)'BaseTHDD',BaseTHDD
       WRITE(fn,*)'bldgH',bldgH
       WRITE(fn,*)'CapMax_dec',CapMax_dec
       WRITE(fn,*)'CapMin_dec',CapMin_dec
       WRITE(fn,*)'CRWmax',CRWmax
       WRITE(fn,*)'CRWmin',CRWmin
       WRITE(fn,*)'dectime',dectime
       WRITE(fn,*)'DecTreeH',DecTreeH
       WRITE(fn,*)'DRAINRT',DRAINRT
       WRITE(fn,*)'EF_umolCO2perJ',EF_umolCO2perJ
       WRITE(fn,*)'EnEF_v_Jkm',EnEF_v_Jkm
       WRITE(fn,*)'EveTreeH',EveTreeH
       WRITE(fn,*)'FAIBldg',FAIBldg
       WRITE(fn,*)'FAIDecTree',FAIDecTree
       WRITE(fn,*)'FAIEveTree',FAIEveTree
       WRITE(fn,*)'Faut',Faut
       WRITE(fn,*)'FcEF_v_kgkm',FcEF_v_kgkm
       WRITE(fn,*)'fcld_obs',fcld_obs
       WRITE(fn,*)'FlowChange',FlowChange
       WRITE(fn,*)'FrFossilFuel_Heat',FrFossilFuel_Heat
       WRITE(fn,*)'FrFossilFuel_NonHeat',FrFossilFuel_NonHeat
       WRITE(fn,*)'G1',G1
       WRITE(fn,*)'G2',G2
       WRITE(fn,*)'G3',G3
       WRITE(fn,*)'G4',G4
       WRITE(fn,*)'G5',G5
       WRITE(fn,*)'G6',G6
       WRITE(fn,*)'InternalWaterUse_h',InternalWaterUse_h
       WRITE(fn,*)'IrrFracConif',IrrFracConif
       WRITE(fn,*)'IrrFracDecid',IrrFracDecid
       WRITE(fn,*)'IrrFracGrass',IrrFracGrass
       WRITE(fn,*)'Kmax',Kmax
       WRITE(fn,*)'LAI_obs',LAI_obs
       WRITE(fn,*)'lat',lat
       WRITE(fn,*)'ldown_obs',ldown_obs
       WRITE(fn,*)'lng',lng
       WRITE(fn,*)'MaxQFMetab',MaxQFMetab
       WRITE(fn,*)'MinQFMetab',MinQFMetab
       WRITE(fn,*)'NARP_EMIS_SNOW',NARP_EMIS_SNOW
       WRITE(fn,*)'NARP_TRANS_SITE',NARP_TRANS_SITE
       WRITE(fn,*)'NumCapita',NumCapita
       WRITE(fn,*)'PipeCapacity',PipeCapacity
       WRITE(fn,*)'PopDensDaytime',PopDensDaytime
       WRITE(fn,*)'PopDensNighttime',PopDensNighttime
       WRITE(fn,*)'PorMax_dec',PorMax_dec
       WRITE(fn,*)'PorMin_dec',PorMin_dec
       WRITE(fn,*)'Precip',Precip
       WRITE(fn,*)'PrecipLimit',PrecipLimit
       WRITE(fn,*)'PrecipLimitAlb',PrecipLimitAlb
       WRITE(fn,*)'Press_hPa',Press_hPa
       WRITE(fn,*)'qh_obs',qh_obs
       WRITE(fn,*)'qn1_obs',qn1_obs
       WRITE(fn,*)'RadMeltFact',RadMeltFact
       WRITE(fn,*)'RAINCOVER',RAINCOVER
       WRITE(fn,*)'RainMaxRes',RainMaxRes
       WRITE(fn,*)'RunoffToWater',RunoffToWater
       WRITE(fn,*)'S1',S1
       WRITE(fn,*)'S2',S2
       WRITE(fn,*)'SnowAlbMax',SnowAlbMax
       WRITE(fn,*)'SnowAlbMin',SnowAlbMin
       WRITE(fn,*)'SnowDensMax',SnowDensMax
       WRITE(fn,*)'SnowDensMin',SnowDensMin
       WRITE(fn,*)'SnowLimBuild',SnowLimBuild
       WRITE(fn,*)'SnowLimPaved',SnowLimPaved
       WRITE(fn,*)'snow_obs',snow_obs
       WRITE(fn,*)'SurfaceArea',SurfaceArea
       WRITE(fn,*)'tau_a',tau_a
       WRITE(fn,*)'tau_f',tau_f
       WRITE(fn,*)'tau_r',tau_r
       WRITE(fn,*)'Temp_C',Temp_C
       WRITE(fn,*)'TempMeltFact',TempMeltFact
       WRITE(fn,*)'TH',TH
       WRITE(fn,*)'timezone',timezone
       WRITE(fn,*)'TL',TL
       WRITE(fn,*)'TrafficUnits',TrafficUnits
       WRITE(fn,*)'xsmd',xsmd
       WRITE(fn,*)'Z',Z
       WRITE(fn,*)'LAIType',LAIType
       WRITE(fn,*)'AH_MIN',AH_MIN
       WRITE(fn,*)'AH_SLOPE_Cooling',AH_SLOPE_Cooling
       WRITE(fn,*)'AH_SLOPE_Heating',AH_SLOPE_Heating
       WRITE(fn,*)'QF0_BEU',QF0_BEU
       WRITE(fn,*)'Qf_A',Qf_A
       WRITE(fn,*)'Qf_B',Qf_B
       WRITE(fn,*)'Qf_C',Qf_C
       WRITE(fn,*)'T_CRITIC_Cooling',T_CRITIC_Cooling
       WRITE(fn,*)'T_CRITIC_Heating',T_CRITIC_Heating
       WRITE(fn,*)'TrafficRate',TrafficRate
       WRITE(fn,*)'Ie_a',Ie_a
       WRITE(fn,*)'Ie_m',Ie_m
       WRITE(fn,*)'MaxConductance',MaxConductance
       WRITE(fn,*)'AHProf_tstep',AHProf_tstep
       WRITE(fn,*)'HumActivity_tstep',HumActivity_tstep
       WRITE(fn,*)'PopProf_tstep',PopProf_tstep
       WRITE(fn,*)'TraffProf_tstep',TraffProf_tstep
       WRITE(fn,*)'WUProfA_tstep',WUProfA_tstep
       WRITE(fn,*)'WUProfM_tstep',WUProfM_tstep
       WRITE(fn,*)'DayWat',DayWat
       WRITE(fn,*)'DayWatPer',DayWatPer
       WRITE(fn,*)'OHM_threshSW',OHM_threshSW
       WRITE(fn,*)'OHM_threshWD',OHM_threshWD
       WRITE(fn,*)'chAnOHM',chAnOHM
       WRITE(fn,*)'cpAnOHM',cpAnOHM
       WRITE(fn,*)'emis',emis
       WRITE(fn,*)'kkAnOHM',kkAnOHM
       WRITE(fn,*)'SatHydraulicConduct',SatHydraulicConduct
       WRITE(fn,*)'sfr',sfr
       WRITE(fn,*)'snowD',snowD
       WRITE(fn,*)'SoilDepth',SoilDepth
       WRITE(fn,*)'soilstoreCap',soilstoreCap
       WRITE(fn,*)'StateLimit',StateLimit
       WRITE(fn,*)'WetThresh',WetThresh
       WRITE(fn,*)'alpha_bioCO2',alpha_bioCO2
       WRITE(fn,*)'alpha_enh_bioCO2',alpha_enh_bioCO2
       WRITE(fn,*)'BaseT',BaseT
       WRITE(fn,*)'BaseTe',BaseTe
       WRITE(fn,*)'beta_bioCO2',beta_bioCO2
       WRITE(fn,*)'beta_enh_bioCO2',beta_enh_bioCO2
       WRITE(fn,*)'GDDFull',GDDFull
       WRITE(fn,*)'LAIMax',LAIMax
       WRITE(fn,*)'LAIMin',LAIMin
       WRITE(fn,*)'min_res_bioCO2',min_res_bioCO2
       WRITE(fn,*)'resp_a',resp_a
       WRITE(fn,*)'resp_b',resp_b
       WRITE(fn,*)'SDDFull',SDDFull
       WRITE(fn,*)'theta_bioCO2',theta_bioCO2
       WRITE(fn,*)'Ts5mindata_ir',Ts5mindata_ir
       WRITE(fn,*)'snowProf',snowProf
       WRITE(fn,*)'WaterDist',WaterDist
       WRITE(fn,*)'OHM_coef',OHM_coef
       WRITE(fn,*)'LAIPower',LAIPower
       WRITE(fn,*)'MetForcingData_grid',MetForcingData_grid
       WRITE(fn,*)'SnowfallCum',SnowfallCum
       WRITE(fn,*)'SnowAlb',SnowAlb
       WRITE(fn,*)'Tair24HR',Tair24HR
       WRITE(fn,*)'qn1_av_store_grid',qn1_av_store_grid
       WRITE(fn,*)'qn1_S_av_store_grid',qn1_S_av_store_grid
       WRITE(fn,*)'albDecTr',albDecTr
       WRITE(fn,*)'albEveTr',albEveTr
       WRITE(fn,*)'albGrass',albGrass
       WRITE(fn,*)'DecidCap',DecidCap
       WRITE(fn,*)'porosity',porosity
       WRITE(fn,*)'GDD',GDD
       WRITE(fn,*)'WU_Day',WU_Day
       WRITE(fn,*)'surf',surf
       WRITE(fn,*)'HDD',HDD
       WRITE(fn,*)'LAI',LAI
       WRITE(fn,*)'alb',alb
       WRITE(fn,*)'IceFrac',IceFrac
       WRITE(fn,*)'MeltWaterStore',MeltWaterStore
       WRITE(fn,*)'SnowDens',SnowDens
       WRITE(fn,*)'snowFrac',snowFrac
       WRITE(fn,*)'SnowPack',SnowPack
       WRITE(fn,*)'soilmoist',soilmoist
       WRITE(fn,*)'state',state
       WRITE(fn,*)'qn1_S_store_grid',qn1_S_store_grid
       WRITE(fn,*)'qn1_store_grid',qn1_store_grid

       CLOSE (fn)
       PRINT*, 'initial condition has been written out for ', Gridiv
       xx=1 ! change flag for completion
    ELSE
       PRINT*, 'initial condition output has been done before for ', Gridiv
       PRINT*, 'see:', FileOut

    END IF

  END SUBROUTINE SUEWS_write_model_state


END MODULE SUEWS_Driver
