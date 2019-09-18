!========================================================================================
! SUEWS driver subroutines
! TS 31 Aug 2017: initial version
! TS 02 Oct 2017: added `SUEWS_cal_Main` as the generic wrapper
! TS 03 Oct 2017: added `SUEWS_cal_AnthropogenicEmission`
MODULE SUEWS_Driver
   USE meteo, ONLY: qsatf, RH2qa, qa2RH
   USE AtmMoistStab_module, ONLY: LUMPS_cal_AtmMoist, STAB_lumps, stab_psi_heat, stab_psi_mom
   USE NARP_MODULE, ONLY: NARP_cal_SunPosition
   USE AnOHM_module, ONLY: AnOHM
   USE resist_module, ONLY: AerodynamicResistance, BoundaryLayerResistance, SurfaceResistance, &
                            cal_z0V, SUEWS_cal_RoughnessParameters
   USE ESTM_module, ONLY: ESTM
   USE Snow_module, ONLY: SnowCalc, Snow_cal_MeltHeat, SnowUpdate
   USE DailyState_module, ONLY: SUEWS_cal_DailyState, update_DailyState
   USE WaterDist_module, ONLY: drainage, soilstore, &
                               SUEWS_cal_SoilState, SUEWS_update_SoilMoist, &
                               ReDistributeWater, SUEWS_cal_HorizontalSoilWater, &
                               SUEWS_cal_WaterUse
   USE ctrl_output, ONLY: varListAll
   USE DailyState_module, ONLY: SUEWS_update_DailyState
   use lumps_module, only: LUMPS_cal_QHQE
   use evap_module, only: evap_SUEWS
   use rsl_module, only: RSLProfile
   use anemsn_module, only: AnthropogenicEmissions
   use CO2_module, only: CO2_biogen
   use evap_module, only: evap_SUEWS
   USE allocateArray, ONLY: &
      nsurf, nvegsurf, &
      PavSurf, BldgSurf, ConifSurf, DecidSurf, GrassSurf, BSoilSurf, WaterSurf, &
      ivConif, ivDecid, ivGrass, &
      ncolumnsDataOutSUEWS, ncolumnsDataOutSnow, &
      ncolumnsDataOutESTM, ncolumnsDataOutDailyState, &
      ncolumnsDataOutRSL
   use moist, only: avcp, avdens, lv_J_kg

   IMPLICIT NONE

CONTAINS
   ! ===================MAIN CALCULATION WRAPPER FOR ENERGY AND WATER FLUX===========
   SUBROUTINE SUEWS_cal_Main( &
      AerodynamicResistanceMethod, AH_MIN, AHProf_24hr, AH_SLOPE_Cooling, & ! input&inout in alphabetical order
      AH_SLOPE_Heating, &
      alb, AlbMax_DecTr, AlbMax_EveTr, AlbMax_Grass, &
      AlbMin_DecTr, AlbMin_EveTr, AlbMin_Grass, &
      alpha_bioCO2, alpha_enh_bioCO2, alt, avkdn, avRh, avU1, BaseT, BaseTe, &
      BaseTHDD, beta_bioCO2, beta_enh_bioCO2, bldgH, CapMax_dec, CapMin_dec, &
      chAnOHM, CO2PointSource, cpAnOHM, CRWmax, CRWmin, DayWat, DayWatPer, &
      DecTreeH, Diagnose, DiagQN, DiagQS, DRAINRT, &
      dt_since_start, dqndt, qn1_av, dqnsdt, qn1_s_av, &
      EF_umolCO2perJ, emis, EmissionsMethod, EnEF_v_Jkm, endDLS, EveTreeH, FAIBldg, &
      FAIDecTree, FAIEveTree, Faut, FcEF_v_kgkm, fcld_obs, FlowChange, &
      FrFossilFuel_Heat, FrFossilFuel_NonHeat, G1, G2, G3, G4, G5, G6, GDD_id, &
      GDDFull, Gridiv, gsModel, HDD_id, HumActivity_24hr, &
      IceFrac, id, Ie_a, Ie_end, Ie_m, Ie_start, imin, &
      InternalWaterUse_h, IrrFracConif, IrrFracDecid, IrrFracGrass, isec, it, EvapMethod, &
      iy, kkAnOHM, Kmax, LAI_id, LAICalcYes, LAIMax, LAIMin, LAI_obs, &
      LAIPower, LAIType, lat, ldown_obs, lng, MaxConductance, MaxFCMetab, MaxQFMetab, &
      SnowWater, MetForcingData_grid, MinFCMetab, MinQFMetab, min_res_bioCO2, &
      NARP_EMIS_SNOW, NARP_TRANS_SITE, NetRadiationMethod, &
      NumCapita, OHM_coef, OHMIncQF, OHM_threshSW, &
      OHM_threshWD, PipeCapacity, PopDensDaytime, &
      PopDensNighttime, PopProf_24hr, PorMax_dec, PorMin_dec, &
      Precip, PrecipLimit, PrecipLimitAlb, Press_hPa, &
      QF0_BEU, Qf_A, Qf_B, Qf_C, &
      qn1_obs, qh_obs, qs_obs, qf_obs, &
      RadMeltFact, RAINCOVER, RainMaxRes, resp_a, resp_b, &
      RoughLenHeatMethod, RoughLenMomMethod, RunoffToWater, S1, S2, &
      SatHydraulicConduct, SDDFull, sfr, SMDMethod, SnowAlb, SnowAlbMax, &
      SnowAlbMin, SnowPackLimit, SnowDens, SnowDensMax, SnowDensMin, SnowfallCum, snowFrac, &
      SnowLimBldg, SnowLimPaved, snowFrac_obs, SnowPack, SnowProf_24hr, snowUse, SoilDepth, &
      soilstore_id, SoilStoreCap, StabilityMethod, startDLS, state_id, StateLimit, &
      StorageHeatMethod, StoreDrainPrm, SurfaceArea, Tair_av, tau_a, tau_f, tau_r, &
      T_CRITIC_Cooling, T_CRITIC_Heating, Temp_C, TempMeltFact, TH, &
      theta_bioCO2, timezone, TL, TrafficRate, TrafficUnits, &
      TraffProf_24hr, Ts5mindata_ir, tstep, tstep_prev, veg_type, &
      WaterDist, WaterUseMethod, WetThresh, wu_m3, &
      WUDay_id, DecidCap_id, albDecTr_id, albEveTr_id, albGrass_id, porosity_id, &
      WUProfA_24hr, WUProfM_24hr, xsmd, Z, z0m_in, zdm_in, &
      datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, dataoutLineRSL, &!output
      DailyStateLine)!output

      IMPLICIT NONE

      ! ########################################################################################
      ! input variables
      INTEGER, INTENT(IN)::AerodynamicResistanceMethod
      INTEGER, INTENT(IN)::Diagnose
      INTEGER, INTENT(IN)::DiagQN
      INTEGER, INTENT(IN)::DiagQS
      INTEGER, INTENT(IN)::startDLS
      INTEGER, INTENT(IN)::endDLS
      INTEGER, INTENT(IN)::EmissionsMethod
      INTEGER, INTENT(IN)::Gridiv
      INTEGER, INTENT(IN)::gsModel
      INTEGER, INTENT(IN)::id
      INTEGER, INTENT(IN)::Ie_end
      INTEGER, INTENT(IN)::Ie_start
      INTEGER, INTENT(IN)::isec
      INTEGER, INTENT(IN)::imin
      INTEGER, INTENT(IN)::it
      INTEGER, INTENT(IN)::EvapMethod
      INTEGER, INTENT(IN)::iy
      INTEGER, INTENT(IN)::LAICalcYes
      INTEGER, INTENT(IN)::NetRadiationMethod
      INTEGER, INTENT(IN)::OHMIncQF
      INTEGER, INTENT(IN)::RoughLenHeatMethod
      INTEGER, INTENT(IN)::RoughLenMomMethod
      INTEGER, INTENT(IN)::SMDMethod
      INTEGER, INTENT(IN)::snowUse
      INTEGER, INTENT(IN)::StabilityMethod
      INTEGER, INTENT(IN)::StorageHeatMethod
      INTEGER, INTENT(IN)::tstep
      INTEGER, INTENT(IN)::tstep_prev ! tstep size of the previous step
      INTEGER, INTENT(in)::dt_since_start ! time since simulation starts [s]
      INTEGER, INTENT(IN)::veg_type
      INTEGER, INTENT(IN)::WaterUseMethod

      REAL(KIND(1D0)), INTENT(IN)::AlbMax_DecTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMax_EveTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMax_Grass
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_DecTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_EveTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_Grass
      REAL(KIND(1D0)), INTENT(IN)::alt
      REAL(KIND(1D0)), INTENT(IN)::avkdn
      REAL(KIND(1D0)), INTENT(IN)::avRh
      REAL(KIND(1D0)), INTENT(IN)::avU1
      REAL(KIND(1D0)), INTENT(IN)::BaseTHDD
      REAL(KIND(1D0)), INTENT(IN)::bldgH
      REAL(KIND(1D0)), INTENT(IN)::CapMax_dec
      REAL(KIND(1D0)), INTENT(IN)::CapMin_dec
      REAL(KIND(1D0)), INTENT(IN)::CO2PointSource
      REAL(KIND(1D0)), INTENT(IN)::CRWmax
      REAL(KIND(1D0)), INTENT(IN)::CRWmin
      REAL(KIND(1D0)), INTENT(IN)::DecTreeH
      REAL(KIND(1D0)), INTENT(IN)::DRAINRT
      REAL(KIND(1D0)), INTENT(IN)::EF_umolCO2perJ
      REAL(KIND(1D0)), INTENT(IN)::EnEF_v_Jkm
      REAL(KIND(1D0)), INTENT(IN)::EveTreeH
      REAL(KIND(1D0)), INTENT(IN)::FAIBldg
      REAL(KIND(1D0)), INTENT(IN)::FAIDecTree
      REAL(KIND(1D0)), INTENT(IN)::FAIEveTree
      REAL(KIND(1D0)), INTENT(IN)::Faut
      REAL(KIND(1D0)), INTENT(IN)::fcld_obs
      REAL(KIND(1D0)), INTENT(IN)::FlowChange
      REAL(KIND(1D0)), INTENT(IN)::FrFossilFuel_Heat
      REAL(KIND(1D0)), INTENT(IN)::FrFossilFuel_NonHeat
      REAL(KIND(1D0)), INTENT(IN)::G1
      REAL(KIND(1D0)), INTENT(IN)::G2
      REAL(KIND(1D0)), INTENT(IN)::G3
      REAL(KIND(1D0)), INTENT(IN)::G4
      REAL(KIND(1D0)), INTENT(IN)::G5
      REAL(KIND(1D0)), INTENT(IN)::G6
      REAL(KIND(1D0)), INTENT(IN)::InternalWaterUse_h
      REAL(KIND(1D0)), INTENT(IN)::IrrFracConif
      REAL(KIND(1D0)), INTENT(IN)::IrrFracDecid
      REAL(KIND(1D0)), INTENT(IN)::IrrFracGrass
      REAL(KIND(1D0)), INTENT(IN)::Kmax
      REAL(KIND(1D0)), INTENT(IN)::LAI_obs
      REAL(KIND(1D0)), INTENT(IN)::lat
      REAL(KIND(1D0)), INTENT(IN)::ldown_obs
      REAL(KIND(1D0)), INTENT(IN)::lng
      REAL(KIND(1D0)), INTENT(IN)::MaxFCMetab
      REAL(KIND(1D0)), INTENT(IN)::MaxQFMetab
      REAL(KIND(1D0)), INTENT(IN)::MinFCMetab
      REAL(KIND(1D0)), INTENT(IN)::MinQFMetab
      REAL(KIND(1D0)), INTENT(IN)::NARP_EMIS_SNOW
      REAL(KIND(1D0)), INTENT(IN)::NARP_TRANS_SITE
      REAL(KIND(1D0)), INTENT(IN)::PipeCapacity
      REAL(KIND(1D0)), INTENT(IN)::PopDensNighttime
      REAL(KIND(1D0)), INTENT(IN)::PorMax_dec
      REAL(KIND(1D0)), INTENT(IN)::PorMin_dec
      REAL(KIND(1D0)), INTENT(IN)::Precip
      REAL(KIND(1D0)), INTENT(IN)::PrecipLimit
      REAL(KIND(1D0)), INTENT(IN)::PrecipLimitAlb
      REAL(KIND(1D0)), INTENT(IN)::Press_hPa
      REAL(KIND(1D0)), INTENT(IN)::qh_obs
      REAL(KIND(1D0)), INTENT(IN)::qn1_obs
      REAL(KIND(1D0)), INTENT(IN)::qs_obs
      REAL(KIND(1D0)), INTENT(IN)::qf_obs
      REAL(KIND(1D0)), INTENT(IN)::RadMeltFact
      REAL(KIND(1D0)), INTENT(IN)::RAINCOVER
      REAL(KIND(1D0)), INTENT(IN)::RainMaxRes
      REAL(KIND(1D0)), INTENT(IN)::RunoffToWater
      REAL(KIND(1D0)), INTENT(IN)::S1
      REAL(KIND(1D0)), INTENT(IN)::S2
      REAL(KIND(1D0)), INTENT(IN)::SnowAlbMax
      REAL(KIND(1D0)), INTENT(IN)::SnowAlbMin
      REAL(KIND(1D0)), INTENT(IN)::SnowDensMax
      REAL(KIND(1D0)), INTENT(IN)::SnowDensMin
      REAL(KIND(1D0)), INTENT(IN)::SnowLimBldg
      REAL(KIND(1D0)), INTENT(IN)::SnowLimPaved
      REAL(KIND(1D0)), INTENT(IN)::snowFrac_obs
      REAL(KIND(1D0)), INTENT(IN)::SurfaceArea
      REAL(KIND(1D0)), INTENT(IN)::tau_a
      REAL(KIND(1D0)), INTENT(IN)::tau_f
      REAL(KIND(1D0)), INTENT(IN)::tau_r
      REAL(KIND(1D0)), INTENT(IN)::Temp_C
      REAL(KIND(1D0)), INTENT(IN)::TempMeltFact
      REAL(KIND(1D0)), INTENT(IN)::TH
      REAL(KIND(1D0)), INTENT(IN)::timezone
      REAL(KIND(1D0)), INTENT(IN)::TL
      REAL(KIND(1D0)), INTENT(IN)::TrafficUnits
      REAL(KIND(1D0)), INTENT(IN)::wu_m3
      REAL(KIND(1D0)), INTENT(IN)::xsmd
      REAL(KIND(1D0)), INTENT(IN)::Z
      REAL(KIND(1D0)), INTENT(IN)::z0m_in
      REAL(KIND(1D0)), INTENT(IN)::zdm_in

      INTEGER, DIMENSION(NVEGSURF), INTENT(IN)::LAIType

      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_MIN
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_SLOPE_Cooling
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_SLOPE_Heating
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::FcEF_v_kgkm
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::QF0_BEU
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_A
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_B
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_C
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::NumCapita
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::PopDensDaytime
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::T_CRITIC_Cooling
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::T_CRITIC_Heating
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::TrafficRate
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::Ie_a
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::Ie_m
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::MaxConductance
      REAL(KIND(1D0)), DIMENSION(7), INTENT(IN)               ::DayWat
      REAL(KIND(1D0)), DIMENSION(7), INTENT(IN)               ::DayWatPer
      REAL(KIND(1D0)), DIMENSION(nsurf + 1), INTENT(IN)         ::OHM_threshSW
      REAL(KIND(1D0)), DIMENSION(nsurf + 1), INTENT(IN)         ::OHM_threshWD
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::chAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::cpAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::emis
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::kkAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SatHydraulicConduct
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::sfr
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SnowPackLimit
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SoilDepth
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SoilStoreCap
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::StateLimit
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::WetThresh
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::alpha_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::alpha_enh_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::BaseT
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::BaseTe
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::beta_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::beta_enh_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::GDDFull
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::LAIMax
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::LAIMin
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::min_res_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::resp_a
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::resp_b
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::SDDFull
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN)          ::SnowProf_24hr
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::theta_bioCO2
      REAL(KIND(1D0)), DIMENSION(4, NVEGSURF), INTENT(IN)      ::LAIPower
      REAL(KIND(1D0)), DIMENSION(nsurf + 1, 4, 3), INTENT(IN)     ::OHM_coef
      REAL(KIND(1D0)), DIMENSION(NSURF + 1, NSURF - 1), INTENT(IN) ::WaterDist
      REAL(KIND(1d0)), DIMENSION(:), INTENT(IN)               ::Ts5mindata_ir
      REAL(KIND(1D0)), DIMENSION(:, :), INTENT(IN)             ::MetForcingData_grid

      ! diurnal profile values for 24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::AHProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::HumActivity_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::PopProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::TraffProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::WUProfA_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::WUProfM_24hr
      ! ########################################################################################

      ! ########################################################################################
      ! inout variables
      ! OHM related:
      REAL(KIND(1d0)), INTENT(INOUT) ::qn1_av
      REAL(KIND(1d0)), INTENT(INOUT) ::dqndt
      REAL(KIND(1d0)), INTENT(INOUT) ::qn1_s_av
      REAL(KIND(1d0)), INTENT(INOUT) ::dqnsdt

      ! snow related:
      REAL(KIND(1D0)), INTENT(INOUT)                  ::SnowfallCum
      REAL(KIND(1D0)), INTENT(INOUT)                  ::SnowAlb
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::IceFrac
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowWater
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowDens
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::snowFrac
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowPack

      ! water balance related:
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::soilstore_id
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::state_id
      REAL(KIND(1D0)), DIMENSION(6, NSURF), INTENT(INOUT) ::StoreDrainPrm

      ! phenology related:
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::alb
      REAL(KIND(1d0)), DIMENSION(5), INTENT(INOUT)       ::GDD_id !Growing Degree Days (see SUEWS_DailyState.f95)
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(INOUT)::LAI_id !LAI for each veg surface [m2 m-2]
      REAL(KIND(1d0)), INTENT(INOUT)                    :: DecidCap_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albDecTr_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albEveTr_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albGrass_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: porosity_id

      ! anthropogenic heat related:
      REAL(KIND(1d0)), DIMENSION(12), INTENT(INOUT) ::HDD_id !Heating Degree Days (see SUEWS_DailyState.f95)

      ! water use related:
      REAL(KIND(1d0)), DIMENSION(9), INTENT(INOUT)  ::WUDay_id

      ! ESTM related:
      REAL(KIND(1d0)), INTENT(INOUT)    ::Tair_av
      ! ########################################################################################

      ! ########################################################################################
      ! output variables
      REAL(KIND(1D0)), DIMENSION(5), INTENT(OUT)                           ::datetimeLine
      REAL(KIND(1D0)), DIMENSION(ncolumnsDataOutSUEWS - 5), INTENT(OUT)      ::dataOutLineSUEWS
      REAL(KIND(1D0)), DIMENSION(ncolumnsDataOutSnow - 5), INTENT(OUT)       ::dataOutLineSnow
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutESTM - 5), INTENT(OUT)       ::dataOutLineESTM
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutRSL - 5), INTENT(OUT):: dataoutLineRSL ! Wind speed array
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutDailyState - 5), INTENT(OUT) ::DailyStateLine
      ! ########################################################################################

      ! ########################################################################################
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
      REAL(KIND(1D0))::Fc_point
      REAL(KIND(1D0))::Fc_respi
      REAL(KIND(1D0))::Fc_traff
      REAL(KIND(1D0))::fwh
      REAL(KIND(1D0))::gfunc
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
      REAL(KIND(1D0))::qn1_snowfree
      REAL(KIND(1D0))::qs
      REAL(KIND(1D0))::RA
      REAL(KIND(1D0))::ResistSurf
      REAL(KIND(1D0))::RH2
      ! REAL(KIND(1D0))::rss
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
      REAL(KIND(1D0))::tskin_C
      REAL(KIND(1D0))::TempVeg
      REAL(KIND(1D0))::tot_chang_per_tstep
      REAL(KIND(1D0))::TStar
      REAL(KIND(1D0))::tsurf
      REAL(KIND(1D0))::UStar
      REAL(KIND(1D0))::VPD_Pa
      ! REAL(KIND(1D0))::WUAreaDecTr_m2
      ! REAL(KIND(1D0))::WUAreaEveTr_m2
      ! REAL(KIND(1D0))::WUAreaGrass_m2
      ! REAL(KIND(1D0))::WUAreaTotal_m2
      REAL(KIND(1D0))::wu_DecTr
      REAL(KIND(1D0))::wu_EveTr
      REAL(KIND(1D0))::wu_Grass
      REAL(KIND(1D0))::z0m
      REAL(KIND(1D0))::zdm
      REAL(KIND(1D0))::ZENITH_deg
      REAL(KIND(1D0))::Zh

      ! REAL(KIND(1d0)), DIMENSION(nvegsurf):: LAI_id_prev !LAI for each veg surface [m2 m-2]

      REAL(KIND(1D0)), DIMENSION(2)::SnowRemoval
      REAL(KIND(1D0)), DIMENSION(NSURF)::chang
      REAL(KIND(1D0)), DIMENSION(NSURF)::changSnow
      REAL(KIND(1D0)), DIMENSION(NSURF)::evap
      REAL(KIND(1D0)), DIMENSION(NSURF)::ev_snow
      REAL(KIND(1D0)), DIMENSION(NSURF)::FreezMelt
      REAL(KIND(1d0)), DIMENSION(nsurf)::kup_ind_snow
      REAL(KIND(1D0)), DIMENSION(NSURF)::mw_ind
      REAL(KIND(1D0)), DIMENSION(NSURF)::Qm_freezState
      REAL(KIND(1D0)), DIMENSION(NSURF)::Qm_melt
      REAL(KIND(1D0)), DIMENSION(NSURF)::Qm_rain
      REAL(KIND(1D0)), DIMENSION(NSURF)::qn1_ind_snow
      REAL(KIND(1D0)), DIMENSION(NSURF)::rainOnSnow
      REAL(KIND(1D0)), DIMENSION(NSURF)::rss_nsurf
      REAL(KIND(1D0)), DIMENSION(NSURF)::runoff
      REAL(KIND(1D0)), DIMENSION(NSURF)::runoffSnow
      REAL(KIND(1D0)), DIMENSION(NSURF)::runoffSoil
      REAL(KIND(1D0)), DIMENSION(NSURF)::smd_nsurf
      REAL(KIND(1D0)), DIMENSION(NSURF)::SnowToSurf
      REAL(KIND(1D0)), DIMENSION(NSURF)::snowDepth

      REAL(KIND(1d0)), DIMENSION(nsurf)::Tsurf_ind_snow

      INTEGER, DIMENSION(NSURF)::snowCalcSwitch
      INTEGER, DIMENSION(3)    ::dayofWeek_id
      INTEGER::DLS

      ! REAL(KIND(1D0))::avcp
      ! REAL(KIND(1D0))::avdens
      REAL(KIND(1D0))::dq
      ! REAL(KIND(1D0))::lv_J_kg
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

      REAL(KIND(1D0)), DIMENSION(NSURF)::deltaQi
      REAL(KIND(1D0)), DIMENSION(NSURF)::drain
      REAL(KIND(1D0)), DIMENSION(NSURF)::FreezState
      REAL(KIND(1D0)), DIMENSION(NSURF)::FreezStateVol
      REAL(KIND(1D0)), DIMENSION(NSURF)::soilstoreOld
      REAL(KIND(1D0)), DIMENSION(NSURF)::stateOld
      REAL(KIND(1D0)), DIMENSION(NSURF)::tsurf_ind

      ! TODO: TS 25 Oct 2017
      ! the `add-*` variables are not used currently as grid-to-grid connection is NOT set up.
      ! set these variables as zero.
      REAL(KIND(1D0))::addImpervious = 0
      REAL(KIND(1D0))::addPipes = 0
      REAL(KIND(1D0))::addVeg = 0
      REAL(KIND(1D0))::addWaterBody = 0
      REAL(KIND(1D0)), DIMENSION(NSURF)::AddWater = 0
      REAL(KIND(1D0)), DIMENSION(NSURF)::AddWaterRunoff = 0

      ! values that are derived from tstep
      INTEGER::nsh ! number of timesteps per hour
      REAL(KIND(1D0))::nsh_real ! nsh in type real
      REAL(KIND(1D0))::tstep_real ! tstep in type real
      REAL(KIND(1D0))::dectime

      ! values that are derived from sfr (surface fractions)
      REAL(KIND(1D0))::VegFraction
      REAL(KIND(1D0))::ImpervFraction
      REAL(KIND(1D0))::PervFraction
      REAL(KIND(1D0))::NonWaterFraction

      ! snow related temporary values
      ! REAL(KIND(1D0))::alb0
      REAL(KIND(1D0))::alb1

      ! ########################################################################################
      ! temporary variables to save values for inout varialbes
      ! OHM related:
      REAL(KIND(1D0))::qn1_av_prev, qn1_av_next
      REAL(KIND(1D0))::dqndt_prev, dqndt_next
      REAL(KIND(1D0))::qn1_s_av_prev, qn1_s_av_next
      REAL(KIND(1D0))::dqnsdt_prev, dqnsdt_next

      ! snow related:
      REAL(KIND(1D0))::SnowfallCum_prev, SnowfallCum_next
      REAL(KIND(1D0))::SnowAlb_prev, SnowAlb_next

      REAL(KIND(1D0)), DIMENSION(NSURF)::IceFrac_prev, IceFrac_next
      REAL(KIND(1D0)), DIMENSION(NSURF)::SnowWater_prev, SnowWater_next
      REAL(KIND(1D0)), DIMENSION(NSURF)::SnowDens_prev, SnowDens_next
      REAL(KIND(1D0)), DIMENSION(NSURF)::SnowFrac_prev, SnowFrac_next
      REAL(KIND(1D0)), DIMENSION(NSURF)::SnowPack_prev, SnowPack_next

      ! water balance related:
      REAL(KIND(1D0)), DIMENSION(NSURF)   ::soilstore_id_prev, soilstore_id_next
      REAL(KIND(1D0)), DIMENSION(NSURF)   ::state_id_prev, state_id_next
      REAL(KIND(1D0)), DIMENSION(6, NSURF)::StoreDrainPrm_prev, StoreDrainPrm_next

      ! phenology related:
      REAL(KIND(1D0)), DIMENSION(NSURF)   ::alb_prev, alb_next
      REAL(KIND(1d0)), DIMENSION(5)       ::GDD_id_prev, GDD_id_next
      REAL(KIND(1d0)), DIMENSION(nvegsurf)::LAI_id_prev, LAI_id_next

      REAL(KIND(1D0))::DecidCap_id_prev, DecidCap_id_next
      REAL(KIND(1D0))::albDecTr_id_prev, albDecTr_id_next
      REAL(KIND(1D0))::albEveTr_id_prev, albEveTr_id_next
      REAL(KIND(1D0))::albGrass_id_prev, albGrass_id_next
      REAL(KIND(1D0))::porosity_id_prev, porosity_id_next

      REAL(KIND(1D0))::Tair_av_prev, Tair_av_next
      ! ########################################################################################

      ! Related to RSL wind profiles
      INTEGER, PARAMETER :: nz = 90   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy
      ! REAL(KIND(1d0)), DIMENSION(nz):: zarrays ! Height array

      ! ########################################################################################
      ! save initial values of inout variables
      qn1_av_prev = qn1_av
      dqndt_prev = dqndt
      qn1_s_av_prev = qn1_s_av
      dqnsdt_prev = dqnsdt
      SnowfallCum_prev = SnowfallCum
      SnowAlb_prev = SnowAlb
      IceFrac_prev = IceFrac
      SnowWater_prev = SnowWater
      SnowDens_prev = SnowDens
      SnowFrac_prev = SnowFrac
      SnowPack_prev = SnowPack
      soilstore_id_prev = soilstore_id
      state_id_prev = state_id
      StoreDrainPrm_prev = StoreDrainPrm
      alb_prev = alb
      GDD_id_prev = GDD_id
      LAI_id_prev = LAI_id
      DecidCap_id_prev = DecidCap_id
      albDecTr_id_prev = albDecTr_id
      albEveTr_id_prev = albEveTr_id
      albGrass_id_prev = albGrass_id
      porosity_id_prev = porosity_id
      Tair_av_prev = Tair_av

      ! ########################################################################################

      ! calculate dectime
      CALL SUEWS_cal_dectime( &
         id, it, imin, isec, & ! input
         dectime) ! output

      ! calculate tstep related VARIABLES
      CALL SUEWS_cal_tstep( &
         tstep, & ! input
         nsh, nsh_real, tstep_real) ! output

      ! calculate surface fraction related VARIABLES
      CALL SUEWS_cal_surf( &
         sfr, & !input
         VegFraction, ImpervFraction, PervFraction, NonWaterFraction) ! output

      ! calculate dayofweek information
      CALL SUEWS_cal_weekday( &
         iy, id, lat, & !input
         dayofWeek_id) !output

      ! calculate dayofweek information
      CALL SUEWS_cal_DLS( &
         id, startDLS, endDLS, & !input
         DLS) !output

      ! calculate mean air temperature of past 24 hours
      Tair_av_next = cal_tair_av(Tair_av_prev, dt_since_start, tstep, temp_c)

      !==============main calculation start=======================

      !==============surface roughness calculation=======================
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_RoughnessParameters...'
      IF (Diagnose == 1) PRINT *, 'z0m_in =', z0m_in
      CALL SUEWS_cal_RoughnessParameters( &
         RoughLenMomMethod, sfr, &!input
         bldgH, EveTreeH, DecTreeH, &
         porosity_id, FAIBldg, FAIEveTree, FAIDecTree, &
         z0m_in, zdm_in, Z, &
         planF, &!output
         Zh, z0m, zdm, ZZD)

      !=================Calculate sun position=================
      IF (Diagnose == 1) WRITE (*, *) 'Calling NARP_cal_SunPosition...'
      CALL NARP_cal_SunPosition( &
         REAL(iy, KIND(1d0)), &!input:
         dectime - tstep/2, &! sun position at middle of timestep before
         timezone, lat, lng, alt, &
         azimuth, zenith_deg)!output:

      !=================Call the SUEWS_cal_DailyState routine to get surface characteristics ready=================
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_DailyState...'
      CALL SUEWS_cal_DailyState( &
         iy, id, it, imin, isec, tstep, tstep_prev, dt_since_start, DayofWeek_id, &!input
         WaterUseMethod, Ie_start, Ie_end, &
         LAICalcYes, LAIType, &
         nsh_real, avkdn, Temp_C, Precip, BaseTHDD, &
         lat, Faut, LAI_obs, &
         AlbMax_DecTr, AlbMax_EveTr, AlbMax_Grass, &
         AlbMin_DecTr, AlbMin_EveTr, AlbMin_Grass, &
         CapMax_dec, CapMin_dec, PorMax_dec, PorMin_dec, &
         Ie_a, Ie_m, DayWatPer, DayWat, &
         BaseT, BaseTe, GDDFull, SDDFull, LAIMin, LAIMax, LAIPower, &
         GDD_id, HDD_id, LAI_id, LAI_id_prev, WUDay_id, &!inout
         DecidCap_id, albDecTr_id, albEveTr_id, albGrass_id, porosity_id, &
         deltaLAI)!output

      !=================Calculation of density and other water related parameters=================
      IF (Diagnose == 1) WRITE (*, *) 'Calling LUMPS_cal_AtmMoist...'
      CALL LUMPS_cal_AtmMoist( &
         Temp_C, Press_hPa, avRh, dectime, &! input:
         lv_J_kg, lvS_J_kg, &! output:
         es_hPa, Ea_hPa, VPd_hpa, VPD_Pa, dq, dens_dry, avcp, avdens)

      !======== Calculate soil moisture =========
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_update_SoilMoist...'
      CALL SUEWS_update_SoilMoist( &
         NonWaterFraction, &!input
         SoilStoreCap, sfr, soilstore_id, &
         SoilMoistCap, SoilState, &!output
         vsmd, smd)

      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_WaterUse...'
      !=================Gives the external and internal water uses per timestep=================
      CALL SUEWS_cal_WaterUse( &
         nsh_real, & ! input:
         wu_m3, SurfaceArea, sfr, &
         IrrFracConif, IrrFracDecid, IrrFracGrass, &
         DayofWeek_id, WUProfA_24hr, WUProfM_24hr, &
         InternalWaterUse_h, HDD_id, WUDay_id, &
         WaterUseMethod, NSH, it, imin, DLS, &
         wu_EveTr, wu_DecTr, wu_Grass, int_wu, ext_wu)! output:

      ! ===================ANTHROPOGENIC HEAT AND CO2 FLUX======================
      CALL SUEWS_cal_AnthropogenicEmission( &
         AH_MIN, AHProf_24hr, AH_SLOPE_Cooling, AH_SLOPE_Heating, CO2PointSource, &! input:
         dayofWeek_id, Diagnose, DLS, EF_umolCO2perJ, EmissionsMethod, EnEF_v_Jkm, &
         FcEF_v_kgkm, FrFossilFuel_Heat, FrFossilFuel_NonHeat, HDD_id, HumActivity_24hr, &
         id, imin, it, MaxFCMetab, MaxQFMetab, MinFCMetab, MinQFMetab, nsh, NumCapita, &
         PopDensDaytime, PopDensNighttime, PopProf_24hr, QF, QF0_BEU, Qf_A, Qf_B, Qf_C, &
         QF_obs, QF_SAHP, sfr, snowFrac, SurfaceArea, T_CRITIC_Cooling, T_CRITIC_Heating, &
         Temp_C, TrafficRate, TrafficUnits, TraffProf_24hr, &
         Fc_anthro, Fc_build, Fc_metab, Fc_point, Fc_traff)! output:

      IF (snowUse == 1) CALL SnowUpdate( &
         nsurf, tstep, &!input
         Temp_C, &
         tau_a, &
         tau_f, &
         tau_r, &
         SnowDensMax, &
         SnowDensMin, &
         SnowAlbMax, &
         SnowAlbMin, &
         SnowPack_prev, SnowAlb_prev, SnowDens_prev, &
         SnowAlb_next, SnowDens_next & ! output
         )

      ! ===================NET ALLWAVE RADIATION================================
      CALL SUEWS_cal_Qn( &
         NetRadiationMethod, snowUse, &!input
         Diagnose, snowFrac_obs, ldown_obs, fcld_obs, &
         dectime, ZENITH_deg, avKdn, Temp_C, avRH, ea_hPa, qn1_obs, &
         SnowAlb_next, snowFrac_prev, DiagQN, &
         NARP_TRANS_SITE, NARP_EMIS_SNOW, IceFrac, sfr, emis, &
         alb, albDecTr_id, DecidCap_id, albEveTr_id, albGrass_id, &!inout
         StoreDrainPrm, &!inout
         ldown, fcld, &!output
         qn1, qn1_snowfree, qn1_S, kclear, kup, lup, tsurf, &
         qn1_ind_snow, kup_ind_snow, Tsurf_ind_snow, Tsurf_ind, &
         alb1, snowFrac_next)

      ! ========================================================================
      ! N.B.: the following parts involves snow-related calculations.
      SnowfallCum = SnowfallCum_prev
      IceFrac = IceFrac_prev
      SnowWater = SnowWater_prev
      SnowPack = SnowPack_prev
      SnowDens = SnowDens_next
      SnowFrac = SnowFrac_next
      SnowAlb = SnowAlb_next

      ! =================STORAGE HEAT FLUX=======================================
      CALL SUEWS_cal_Qs( &
         StorageHeatMethod, qs_obs, OHMIncQF, Gridiv, &!input
         id, tstep, dt_since_start, Diagnose, sfr, &
         OHM_coef, OHM_threshSW, OHM_threshWD, &
         soilstore_id, SoilStoreCap, state_id, SnowUse, DiagQS, &
         HDD_id, MetForcingData_grid, Ts5mindata_ir, qf, qn1, &
         avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown, &
         bldgh, alb, emis, cpAnOHM, kkAnOHM, chAnOHM, EmissionsMethod, &
         Tair_av_next, qn1_av_prev, dqndt_prev, qn1_s_av_prev, dqnsdt_prev, &
         StoreDrainPrm, &!inout
         qn1_S, snowFrac, dataOutLineESTM, qs, &!output
         qn1_av_next, dqndt_next, qn1_s_av_next, dqnsdt_next, &
         deltaQi, a1, a2, a3)

      !==================Energy related to snow melting/freezing processes=======
      IF (Diagnose == 1) WRITE (*, *) 'Calling MeltHeat'
      CALL Snow_cal_MeltHeat( &
         snowUse, &!input
         lvS_J_kg, lv_J_kg, tstep_real, RadMeltFact, TempMeltFact, SnowAlbMax, &
         SnowDensMin, Temp_C, Precip, PrecipLimit, PrecipLimitAlb, &
         nsh_real, sfr, Tsurf_ind, Tsurf_ind_snow, state_id, qn1_ind_snow, &
         kup_ind_snow, SnowWater, deltaQi, alb1, &
         SnowPack, snowFrac, SnowAlb, SnowDens, SnowfallCum, &!inout
         mwh, Qm, QmFreez, QmRain, &! output
         veg_fr, snowCalcSwitch, Qm_melt, Qm_freezState, Qm_rain, FreezMelt, &
         FreezState, FreezStateVol, rainOnSnow, SnowDepth, mw_ind, &
         dataOutLineSnow)!output

      !==========================Turbulent Fluxes================================
      IF (Diagnose == 1) WRITE (*, *) 'Calling LUMPS_cal_QHQE...'
      !Calculate QH and QE from LUMPS
      CALL LUMPS_cal_QHQE( &
         veg_type, & !input
         snowUse, qn1, qf, qs, Qm, Temp_C, Veg_Fr, avcp, Press_hPa, lv_J_kg, &
         tstep_real, DRAINRT, nsh_real, &
         Precip, RainMaxRes, RAINCOVER, sfr, LAI_id, LAImax, LAImin, &
         H_mod, & !output
         E_mod, psyc_hPa, s_hPa, sIce_hpa, TempVeg, VegPhenLumps)

      !============= calculate water balance =============
      CALL SUEWS_cal_Water( &
         Diagnose, &!input
         snowUse, NonWaterFraction, addPipes, addImpervious, addVeg, addWaterBody, &
         state_id, soilstore_id, sfr, StoreDrainPrm, WaterDist, nsh_real, &
         drain_per_tstep, &  !output
         drain, AddWaterRunoff, &
         AdditionalWater, runoffPipes, runoff_per_interval, &
         AddWater, stateOld, soilstoreOld)
      !============= calculate water balance end =============

      !===============Resistance Calculations=======================
      CALL SUEWS_cal_Resistance( &
         StabilityMethod, &!input:
         Diagnose, AerodynamicResistanceMethod, RoughLenHeatMethod, snowUse, &
         id, it, gsModel, SMDMethod, &
         qh_obs, avdens, avcp, h_mod, qn1, dectime, zzd, z0m, zdm, &
         avU1, Temp_C, VegFraction, avkdn, &
         Kmax, &
         g1, g2, g3, g4, &
         g5, g6, s1, s2, &
         th, tl, &
         dq, xsmd, vsmd, MaxConductance, LAIMax, LAI_id, snowFrac, sfr, &
         UStar, TStar, L_mod, &!output
         zL, gsc, ResistSurf, RA, RAsnow, rb)

      !======== Evaporation and surface state_id ========
      CALL SUEWS_cal_QE( &
         Diagnose, snowuse, &!input
         tstep, imin, it, EvapMethod, snowCalcSwitch, dayofWeek_id, CRWmin, CRWmax, &
         dectime, lvS_J_kg, lv_j_kg, avdens, avRh, Press_hPa, Temp_C, &
         RAsnow, psyc_hPa, avcp, sIce_hPa, &
         PervFraction, vegfraction, addimpervious, qn1_snowfree, qf, qs, vpd_hPa, s_hPa, &
         ResistSurf, RA, rb, snowdensmin, precip, PipeCapacity, RunoffToWater, &
         NonWaterFraction, wu_EveTr, wu_DecTr, wu_Grass, addVeg, addWaterBody, SnowLimPaved, SnowLimBldg, &
         SurfaceArea, FlowChange, drain, WetThresh, stateOld, mw_ind, SoilStoreCap, rainonsnow, &
         freezmelt, freezstate, freezstatevol, Qm_Melt, Qm_rain, Tsurf_ind, sfr, &
         StateLimit, AddWater, addwaterrunoff, StoreDrainPrm, SnowPackLimit, SnowProf_24hr, &
         runoff_per_interval, state_id, & !inout:
         soilstore_id, SnowPack, snowFrac, SnowWater, iceFrac, SnowDens, &
         runoffSnow, runoff, runoffSoil, chang, changSnow, &! output:
         snowDepth, SnowToSurf, ev_snow, SnowRemoval, &
         evap, rss_nsurf, p_mm, qe, state_per_tstep, NWstate_per_tstep, qeOut, &
         swe, ev, chSnow_per_interval, ev_per_tstep, qe_per_tstep, runoff_per_tstep, &
         surf_chang_per_tstep, runoffPipes, mwstore, runoffwaterbody, &
         runoffAGveg, runoffAGimpervious, runoffWaterBody_m3, runoffPipes_m3)
      !======== Evaporation and surface state_id end========

      !============ Sensible heat flux ===============
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_QH...'
      CALL SUEWS_cal_QH( &
         1, &
         qn1, qf, QmRain, qeOut, qs, QmFreez, qm, avdens, avcp, tsurf, Temp_C, RA, &
         qh, qh_residual, qh_resist)!output
      !============ Sensible heat flux end===============

      SnowfallCum_next = SnowfallCum
      SnowAlb_next = SnowAlb
      IceFrac_next = IceFrac
      SnowWater_next = SnowWater
      SnowDens_next = SnowDens
      SnowFrac_next = SnowFrac
      SnowPack_next = SnowPack

      ! N.B.: snow-related calculations end here.
      !===================================================

      !=== Horizontal movement between soil stores ===
      ! Now water is allowed to move horizontally between the soil stores
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_HorizontalSoilWater...'
      CALL SUEWS_cal_HorizontalSoilWater( &
         sfr, &! input: ! surface fractions
         SoilStoreCap, &!Capacity of soil store for each surface [mm]
         SoilDepth, &!Depth of sub-surface soil store for each surface [mm]
         SatHydraulicConduct, &!Saturated hydraulic conductivity for each soil subsurface [mm s-1]
         SurfaceArea, &!Surface area of the study area [m2]
         NonWaterFraction, &! sum of surface cover fractions for all except water surfaces
         tstep_real, & !tstep cast as a real for use in calculations
         soilstore_id, &! inout:!Soil moisture of each surface type [mm]
         runoffSoil, &!Soil runoff from each soil sub-surface [mm]
         runoffSoil_per_tstep &!  output:!Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)
         )

      !========== Calculate soil moisture ============
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_SoilState...'
      CALL SUEWS_cal_SoilState( &
         SMDMethod, xsmd, NonWaterFraction, SoilMoistCap, &!input
         SoilStoreCap, surf_chang_per_tstep, &
         soilstore_id, soilstoreOld, sfr, &
         smd, smd_nsurf, tot_chang_per_tstep, SoilState)!output

      !============ surface-level diagonostics ===============
      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_Diagnostics...'
      CALL SUEWS_cal_Diagnostics( &
         dectime, &!input
         avU1, Temp_C, avRH, Press_hPa, &
         qh, qe, &
         VegFraction, z, z0m, zdm, RA, avdens, avcp, lv_J_kg, tstep_real, &
         RoughLenHeatMethod, StabilityMethod, &
         avU10_ms, t2_C, q2_gkg, tskin_C, RH2)!output

      IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_cal_Diagnostics...'
      CALL RSLProfile( &
         UStar, &!input
         L_mod, sfr, Zh, planF, &
         StabilityMethod, Temp_C, avRh, Press_hPa, z, TStar, qeOut, &
         dataoutLineRSL)!output

      !============ surface-level diagonostics end ===============

      ! ============ BIOGENIC CO2 FLUX =======================
      CALL SUEWS_cal_BiogenCO2( &
         alpha_bioCO2, alpha_enh_bioCO2, avkdn, avRh, beta_bioCO2, beta_enh_bioCO2, BSoilSurf, &! input:
         ConifSurf, DecidSurf, dectime, Diagnose, EmissionsMethod, Fc_anthro, G1, G2, G3, G4, &
         G5, G6, gfunc, GrassSurf, gsmodel, id, it, ivConif, ivDecid, ivGrass, Kmax, LAI_id, LAIMin, &
         LAIMax, MaxConductance, min_res_bioCO2, nsurf, NVegSurf, Press_hPa, resp_a, &
         resp_b, S1, S2, sfr, SMDMethod, snowFrac, t2_C, Temp_C, theta_bioCO2, TH, TL, vsmd, xsmd, &
         Fc, Fc_biogen, Fc_photo, Fc_respi)! output:

      !==============main calculation end=======================

      !==============================================================
      ! update inout variables with new values
      qn1_av = qn1_av_next
      dqndt = dqndt_next
      qn1_s_av = qn1_s_av_next
      dqnsdt = dqnsdt_next
      SnowfallCum = SnowfallCum_next
      SnowAlb = SnowAlb_next
      IceFrac = IceFrac_next
      SnowWater = SnowWater_next
      SnowDens = SnowDens_next
      SnowFrac = SnowFrac_next
      SnowPack = SnowPack_next

      ! soilstore_id = soilstore_id_next
      ! state_id = state_id_next
      ! StoreDrainPrm = StoreDrainPrm_next
      ! alb = alb_next
      ! GDD_id = GDD_id_next
      ! LAI_id = LAI_id_next
      ! DecidCap_id = DecidCap_id_next
      ! albDecTr_id = albDecTr_id_next
      ! albEveTr_id = albEveTr_id_next
      ! albGrass_id = albGrass_id_next
      ! porosity_id = porosity_id_next
      Tair_av = Tair_av_next

      !==============translation of  output variables into output array===========
      CALL SUEWS_update_outputLine( &
         AdditionalWater, alb, avkdn, avU10_ms, azimuth, &!input
         chSnow_per_interval, dectime, &
         drain_per_tstep, E_mod, ev_per_tstep, ext_wu, Fc, Fc_build, fcld, &
         Fc_metab, Fc_photo, Fc_respi, Fc_point, Fc_traff, FlowChange, &
         h_mod, id, imin, int_wu, it, iy, &
         kup, LAI_id, ldown, l_mod, lup, mwh, &
         MwStore, &
         nsh_real, NWstate_per_tstep, Precip, q2_gkg, &
         qeOut, qf, qh, qh_resist, Qm, QmFreez, &
         QmRain, qn1, qn1_S, qn1_snowfree, qs, RA, &
         resistsurf, RH2, runoffAGimpervious, runoffAGveg, &
         runoff_per_tstep, runoffPipes, runoffSoil_per_tstep, &
         runoffWaterBody, sfr, smd, smd_nsurf, SnowAlb, SnowRemoval, &
         state_id, state_per_tstep, surf_chang_per_tstep, swe, t2_C, tskin_C, &
         tot_chang_per_tstep, tsurf, UStar, wu_DecTr, &
         wu_EveTr, wu_Grass, z0m, zdm, zenith_deg, &
         datetimeLine, dataOutLineSUEWS)!output

      ! model state_id:

      ! daily state_id:
      CALL update_DailyState( &
         it, imin, nsh_real, &!input
         GDD_id, HDD_id, LAI_id, &
         DecidCap_id, &
         albDecTr_id, &
         albEveTr_id, &
         albGrass_id, &
         porosity_id, &
         WUDay_id, &
         deltaLAI, VegPhenLumps, &
         SnowAlb, SnowDens, &
         a1, a2, a3, &
         DailyStateLine)!out

      !==============translation end ================

   END SUBROUTINE SUEWS_cal_Main
   ! ================================================================================

   ! ===================ANTHROPOGENIC HEAT + CO2 FLUX================================
   SUBROUTINE SUEWS_cal_AnthropogenicEmission( &
      AH_MIN, AHProf_24hr, AH_SLOPE_Cooling, AH_SLOPE_Heating, CO2PointSource, &! input:
      dayofWeek_id, Diagnose, DLS, EF_umolCO2perJ, EmissionsMethod, EnEF_v_Jkm, &
      FcEF_v_kgkm, FrFossilFuel_Heat, FrFossilFuel_NonHeat, HDD_id, HumActivity_24hr, &
      id, imin, it, MaxFCMetab, MaxQFMetab, MinFCMetab, MinQFMetab, nsh, NumCapita, &
      PopDensDaytime, PopDensNighttime, PopProf_24hr, QF, QF0_BEU, Qf_A, Qf_B, Qf_C, &
      QF_obs, QF_SAHP, sfr, snowFrac, SurfaceArea, T_CRITIC_Cooling, T_CRITIC_Heating, &
      Temp_C, TrafficRate, TrafficUnits, TraffProf_24hr, &
      Fc_anthro, Fc_build, Fc_metab, Fc_point, Fc_traff)! output:

      IMPLICIT NONE

      INTEGER, INTENT(in)::Diagnose
      INTEGER, INTENT(in)::DLS
      INTEGER, INTENT(in)::EmissionsMethod
      INTEGER, INTENT(in)::id
      INTEGER, INTENT(in)::it
      INTEGER, INTENT(in)::imin
      INTEGER, INTENT(in)::nsh
      INTEGER, DIMENSION(3), INTENT(in)::dayofWeek_id

      REAL(KIND(1d0)), DIMENSION(6, 2), INTENT(in)::HDD_id

      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::AH_MIN
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::AH_SLOPE_Heating
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::AH_SLOPE_Cooling
      REAL(KIND(1D0)), DIMENSION(2), INTENT(in)::FcEF_v_kgkm
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::NumCapita
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::PopDensDaytime
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::QF0_BEU
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::Qf_A
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::Qf_B
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::Qf_C
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::T_CRITIC_Heating
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::T_CRITIC_Cooling
      REAL(KIND(1d0)), DIMENSION(2), INTENT(in)::TrafficRate

      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::AHProf_24hr
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::HumActivity_24hr
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::TraffProf_24hr
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::PopProf_24hr

      REAL(KIND(1D0)), INTENT(in)::CO2PointSource
      REAL(KIND(1D0)), INTENT(in)::EF_umolCO2perJ
      REAL(KIND(1D0)), INTENT(in)::EnEF_v_Jkm
      REAL(KIND(1D0)), INTENT(in)::FrFossilFuel_Heat
      REAL(KIND(1D0)), INTENT(in)::FrFossilFuel_NonHeat
      REAL(KIND(1D0)), INTENT(in)::MaxFCMetab
      REAL(KIND(1D0)), INTENT(in)::MaxQFMetab
      REAL(KIND(1D0)), INTENT(in)::MinFCMetab
      REAL(KIND(1D0)), INTENT(in)::MinQFMetab
      REAL(KIND(1D0)), INTENT(in)::PopDensNighttime
      REAL(KIND(1D0)), INTENT(in)::QF_obs
      REAL(KIND(1D0)), INTENT(in)::Temp_C
      REAL(KIND(1D0)), INTENT(in)::TrafficUnits

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::snowFrac
      REAL(KIND(1D0)), INTENT(IN)::SurfaceArea

      REAL(KIND(1D0)), INTENT(out)::Fc_anthro
      REAL(KIND(1D0)), INTENT(out)::Fc_build
      REAL(KIND(1D0)), INTENT(out)::Fc_metab
      REAL(KIND(1D0)), INTENT(out)::Fc_point
      REAL(KIND(1D0)), INTENT(out)::Fc_traff
      REAL(KIND(1D0)), INTENT(out)::QF
      REAL(KIND(1D0)), INTENT(out)::QF_SAHP

      INTEGER, PARAMETER :: notUsedI = -999
      REAL(KIND(1D0)), PARAMETER::notUsed = -999

      IF (EmissionsMethod == 0) THEN ! use observed qf
         qf = QF_obs
      ELSEIF ((EmissionsMethod > 0 .AND. EmissionsMethod <= 6) .OR. EmissionsMethod >= 11) THEN
         CALL AnthropogenicEmissions( &
            CO2PointSource, EmissionsMethod, &
            id, it, imin, DLS, nsh, DayofWeek_id, &
            EF_umolCO2perJ, FcEF_v_kgkm, EnEF_v_Jkm, TrafficUnits, &
            FrFossilFuel_Heat, FrFossilFuel_NonHeat, &
            MinFCMetab, MaxFCMetab, MinQFMetab, MaxQFMetab, &
            NumCapita, PopDensDaytime, PopDensNighttime, &
            Temp_C, HDD_id, Qf_A, Qf_B, Qf_C, &
            AH_MIN, AH_SLOPE_Heating, AH_SLOPE_Cooling, &
            T_CRITIC_Heating, T_CRITIC_Cooling, &
            TrafficRate, &
            QF0_BEU, QF_SAHP, &
            Fc_anthro, Fc_metab, Fc_traff, Fc_build, Fc_point, &
            AHProf_24hr, HumActivity_24hr, TraffProf_24hr, PopProf_24hr, SurfaceArea)

      ELSE
         CALL ErrorHint(73, 'RunControl.nml:EmissionsMethod unusable', notUsed, notUsed, EmissionsMethod)
      ENDIF

      IF (EmissionsMethod >= 1) qf = QF_SAHP

      IF (EmissionsMethod >= 0 .AND. EmissionsMethod <= 6) THEN
         Fc_anthro = 0
         Fc_metab = 0
         Fc_traff = 0
         Fc_build = 0
         Fc_point = 0
      ENDIF

   END SUBROUTINE SUEWS_cal_AnthropogenicEmission
   ! ================================================================================

   !==============BIOGENIC CO2 flux==================================================
   SUBROUTINE SUEWS_cal_BiogenCO2( &
      alpha_bioCO2, alpha_enh_bioCO2, avkdn, avRh, beta_bioCO2, beta_enh_bioCO2, BSoilSurf, &! input:
      ConifSurf, DecidSurf, dectime, Diagnose, EmissionsMethod, Fc_anthro, G1, G2, G3, G4, &
      G5, G6, gfunc, GrassSurf, gsmodel, id, it, ivConif, ivDecid, ivGrass, Kmax, LAI_id, LAIMin, &
      LAIMax, MaxConductance, min_res_bioCO2, nsurf, NVegSurf, Press_hPa, resp_a, &
      resp_b, S1, S2, sfr, SMDMethod, snowFrac, t2_C, Temp_C, theta_bioCO2, TH, TL, vsmd, xsmd, &
      Fc, Fc_biogen, Fc_photo, Fc_respi)! output:

      IMPLICIT NONE

      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::alpha_bioCO2
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::alpha_enh_bioCO2
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::beta_bioCO2
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::beta_enh_bioCO2
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::LAI_id
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::LAIMin
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::LAIMax
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::min_res_bioCO2
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::resp_a
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::resp_b
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(in)::theta_bioCO2

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::snowFrac

      REAL(KIND(1d0)), DIMENSION(3), INTENT(in)::MaxConductance

      INTEGER, INTENT(in)::BSoilSurf
      INTEGER, INTENT(in)::ConifSurf
      INTEGER, INTENT(in)::DecidSurf
      INTEGER, INTENT(in)::Diagnose
      INTEGER, INTENT(in)::EmissionsMethod
      INTEGER, INTENT(in)::GrassSurf
      INTEGER, INTENT(in)::gsmodel
      INTEGER, INTENT(in)::id
      INTEGER, INTENT(in)::it
      INTEGER, INTENT(in)::ivConif
      INTEGER, INTENT(in)::ivDecid
      INTEGER, INTENT(in)::ivGrass
      INTEGER, INTENT(in)::nsurf
      INTEGER, INTENT(in)::NVegSurf
      INTEGER, INTENT(in)::SMDMethod

      REAL(KIND(1d0)), INTENT(in)::avkdn
      REAL(KIND(1d0)), INTENT(in)::avRh
      REAL(KIND(1d0)), INTENT(in)::dectime
      REAL(KIND(1d0)), INTENT(in)::Fc_anthro
      REAL(KIND(1d0)), INTENT(in)::G1
      REAL(KIND(1d0)), INTENT(in)::G2
      REAL(KIND(1d0)), INTENT(in)::G3
      REAL(KIND(1d0)), INTENT(in)::G4
      REAL(KIND(1d0)), INTENT(in)::G5
      REAL(KIND(1d0)), INTENT(in)::G6
      REAL(KIND(1d0)), INTENT(in)::gfunc
      REAL(KIND(1d0)), INTENT(in)::Kmax
      REAL(KIND(1d0)), INTENT(in)::Press_hPa
      REAL(KIND(1d0)), INTENT(in)::S1
      REAL(KIND(1d0)), INTENT(in)::S2
      REAL(KIND(1d0)), INTENT(in)::t2_C
      REAL(KIND(1d0)), INTENT(in)::Temp_C
      REAL(KIND(1d0)), INTENT(in)::TH
      REAL(KIND(1d0)), INTENT(in)::TL
      REAL(KIND(1d0)), INTENT(in)::vsmd
      REAL(KIND(1d0)), INTENT(in)::xsmd

      REAL(KIND(1d0)), INTENT(out)::Fc_biogen
      REAL(KIND(1d0)), INTENT(out)::Fc_photo
      REAL(KIND(1d0)), INTENT(out)::Fc_respi
      REAL(KIND(1d0)), INTENT(out)::Fc

      REAL(KIND(1d0))::gfunc2
      REAL(KIND(1d0))::dq
      REAL(KIND(1d0))::t2
      REAL(KIND(1d0))::dummy1
      REAL(KIND(1d0))::dummy2
      REAL(KIND(1d0))::dummy3
      REAL(KIND(1d0))::dummy4
      REAL(KIND(1d0))::dummy5
      REAL(KIND(1d0))::dummy6
      REAL(KIND(1d0))::dummy7
      REAL(KIND(1d0))::dummy8
      REAL(KIND(1d0))::dummy9
      REAL(KIND(1d0))::dummy10
      REAL(KIND(1d0))::dummy11

      IF (EmissionsMethod >= 11) THEN

         IF (gsmodel == 3 .OR. gsmodel == 4) THEN ! With modelled 2 meter temperature
            ! Call LUMPS_cal_AtmMoist for dq and SurfaceResistance for gfunc with 2 meter temperature
            ! If modelled 2 meter temperature is too different from measured air temperature then
            ! use temp_c
            IF (ABS(Temp_C - t2_C) > 5) THEN
               t2 = Temp_C
            ELSE
               t2 = t2_C
            ENDIF

            CALL LUMPS_cal_AtmMoist( &
               t2, Press_hPa, avRh, dectime, &! input:
               dummy1, dummy2, &! output:
               dummy3, dummy4, dummy5, dummy6, dq, dummy7, dummy8, dummy9)

            CALL SurfaceResistance( &
               id, it, &! input:
               SMDMethod, snowFrac, sfr, avkdn, t2, dq, xsmd, vsmd, MaxConductance, &
               LAIMax, LAI_id, gsModel, Kmax, &
               G1, G2, G3, G4, G5, G6, TH, TL, S1, S2, &
               gfunc2, dummy10, dummy11)! output:
         ENDIF

         ! Calculate CO2 fluxes from biogenic components
         IF (Diagnose == 1) WRITE (*, *) 'Calling CO2_biogen...'
         CALL CO2_biogen( &
            alpha_bioCO2, alpha_enh_bioCO2, avkdn, beta_bioCO2, beta_enh_bioCO2, BSoilSurf, &! input:
            ConifSurf, DecidSurf, dectime, EmissionsMethod, gfunc, gfunc2, GrassSurf, gsmodel, &
            id, it, ivConif, ivDecid, ivGrass, LAI_id, LAIMin, LAIMax, min_res_bioCO2, nsurf, &
            NVegSurf, resp_a, resp_b, sfr, snowFrac, t2, Temp_C, theta_bioCO2, &
            Fc_biogen, Fc_photo, Fc_respi)! output:
      ENDIF

      IF (EmissionsMethod >= 0 .AND. EmissionsMethod <= 6) THEN
         Fc_biogen = 0
         Fc_photo = 0
         Fc_respi = 0
      ENDIF

      Fc = Fc_anthro + Fc_biogen

   END SUBROUTINE SUEWS_cal_BiogenCO2
   !========================================================================

   !=============net all-wave radiation=====================================
   SUBROUTINE SUEWS_cal_Qn( &
      NetRadiationMethod, snowUse, &!input
      Diagnose, snowFrac_obs, ldown_obs, fcld_obs, &
      dectime, ZENITH_deg, avKdn, Temp_C, avRH, ea_hPa, qn1_obs, &
      SnowAlb, snowFrac_prev, DiagQN, &
      NARP_TRANS_SITE, NARP_EMIS_SNOW, IceFrac, sfr, emis, &
      alb, albDecTr_id, DecidCap_id, albEveTr_id, albGrass_id, &!inout
      StoreDrainPrm, &!inout
      ldown, fcld, &!output
      qn1, qn1_snowfree, qn1_S, kclear, kup, lup, tsurf, &
      qn1_ind_snow, kup_ind_snow, Tsurf_ind_snow, Tsurf_ind, &
      alb1, snowFrac_next)
      USE NARP_MODULE, ONLY: RadMethod, NARP

      IMPLICIT NONE
      ! INTEGER,PARAMETER ::nsurf     = 7 ! number of surface types
      ! INTEGER,PARAMETER ::ConifSurf = 3 !New surface classes: Grass = 5th/7 surfaces
      ! INTEGER,PARAMETER ::DecidSurf = 4 !New surface classes: Grass = 5th/7 surfaces
      ! INTEGER,PARAMETER ::GrassSurf = 5

      INTEGER, INTENT(in)::NetRadiationMethod
      INTEGER, INTENT(in)::snowUse
      INTEGER, INTENT(in)::Diagnose
      INTEGER, INTENT(in)::DiagQN

      REAL(KIND(1d0)), INTENT(in)::snowFrac_obs
      REAL(KIND(1d0)), INTENT(in)::ldown_obs
      REAL(KIND(1d0)), INTENT(in)::fcld_obs
      REAL(KIND(1d0)), INTENT(in)::dectime
      REAL(KIND(1d0)), INTENT(in)::ZENITH_deg
      REAL(KIND(1d0)), INTENT(in)::avKdn
      REAL(KIND(1d0)), INTENT(in)::Temp_C
      REAL(KIND(1d0)), INTENT(in)::avRH
      REAL(KIND(1d0)), INTENT(in)::ea_hPa
      REAL(KIND(1d0)), INTENT(in)::qn1_obs
      REAL(KIND(1d0)), INTENT(in)::SnowAlb
      REAL(KIND(1d0)), INTENT(in)::NARP_EMIS_SNOW
      REAL(KIND(1d0)), INTENT(in)::NARP_TRANS_SITE

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in):: IceFrac
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in):: sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in):: emis

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)  ::alb
      REAL(KIND(1d0)), INTENT(in)  ::albDecTr_id
      REAL(KIND(1d0)), INTENT(in)  ::DecidCap_id
      REAL(KIND(1d0)), INTENT(in)  ::albEveTr_id
      REAL(KIND(1d0)), INTENT(in)  ::albGrass_id

      REAL(KIND(1d0)), DIMENSION(6, nsurf), INTENT(inout)::StoreDrainPrm

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::snowFrac_prev
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::snowFrac_next
      REAL(KIND(1d0)), DIMENSION(nsurf)::snowFrac

      REAL(KIND(1d0)), INTENT(out)::ldown
      REAL(KIND(1d0)), INTENT(out)::fcld
      REAL(KIND(1d0)), INTENT(out)::qn1
      REAL(KIND(1d0)), INTENT(out)::qn1_snowfree
      REAL(KIND(1d0)), INTENT(out)::qn1_S
      REAL(KIND(1d0)), INTENT(out)::kclear
      REAL(KIND(1d0)), INTENT(out)::kup
      REAL(KIND(1d0)), INTENT(out)::lup
      REAL(KIND(1d0)), INTENT(out)::tsurf
      REAL(KIND(1d0)), INTENT(out)::alb1
      REAL(KIND(1d0))::alb0

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out) ::qn1_ind_snow
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out) ::kup_ind_snow
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out) ::Tsurf_ind_snow
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: tsurf_ind

      REAL(KIND(1d0)), DIMENSION(nsurf):: lup_ind
      REAL(KIND(1d0)), DIMENSION(nsurf):: kup_ind
      REAL(KIND(1d0)), DIMENSION(nsurf):: qn1_ind

      REAL(KIND(1d0)), PARAMETER::NAN = -999
      INTEGER :: NetRadiationMethodX
      INTEGER::AlbedoChoice, ldown_option

      CALL RadMethod( &
         NetRadiationMethod, &!inout
         snowUse, &!input
         NetRadiationMethodX, AlbedoChoice, ldown_option)!output

      snowFrac = snowFrac_prev
      IF (NetRadiationMethodX > 0) THEN

         ! IF (snowUse==0) snowFrac=snowFrac_obs
         IF (snowUse == 0) snowFrac = 0

         IF (ldown_option == 1) THEN !Observed ldown provided as forcing
            ldown = ldown_obs
         ELSE
            ldown = -9              !to be filled in NARP
         ENDIF

         IF (ldown_option == 2) THEN !observed cloud fraction provided as forcing
            fcld = fcld_obs
         ENDIF

         !write(*,*) DecidCap(id), id, it, imin, 'Calc - near start'

         ! Update variables that change daily and represent seasonal variability
         alb(DecidSurf) = albDecTr_id !Change deciduous albedo
         StoreDrainPrm(6, DecidSurf) = DecidCap_id !Change current storage capacity of deciduous trees
         ! Change EveTr and Grass albedo too
         alb(ConifSurf) = albEveTr_id
         alb(GrassSurf) = albGrass_id

         IF (Diagnose == 1) WRITE (*, *) 'Calling NARP...'
         IF (Diagqn == 1) WRITE (*, *) 'NetRadiationMethodX:', NetRadiationMethodX
         IF (Diagqn == 1) WRITE (*, *) 'AlbedoChoice:', AlbedoChoice

         CALL NARP( &
            nsurf, sfr, snowFrac, alb, emis, IceFrac, &! input:
            NARP_TRANS_SITE, NARP_EMIS_SNOW, &
            dectime, ZENITH_deg, avKdn, Temp_C, avRH, ea_hPa, qn1_obs, &
            SnowAlb, &
            AlbedoChoice, ldown_option, NetRadiationMethodX, DiagQN, &
            qn1, qn1_snowfree, qn1_S, kclear, kup, LDown, lup, fcld, tsurf, &! output:
            qn1_ind_snow, kup_ind_snow, Tsurf_ind_snow, Tsurf_ind, alb0, alb1)

      ELSE ! NetRadiationMethod==0
         snowFrac = snowFrac_obs
         qn1 = qn1_obs
         qn1_snowfree = qn1_obs
         qn1_s = qn1_obs
         ldown = NAN
         lup = NAN
         kup = NAN
         tsurf = NAN
         lup_ind = NAN
         kup_ind = NAN
         tsurf_ind = NAN
         qn1_ind = NAN
         Fcld = NAN
      ENDIF
      snowFrac_next = snowFrac

      IF (ldown_option == 1) THEN
         Fcld = NAN
      ENDIF

   END SUBROUTINE SUEWS_cal_Qn
   !========================================================================

   !=============storage heat flux=========================================
   SUBROUTINE SUEWS_cal_Qs( &
      StorageHeatMethod, qs_obs, OHMIncQF, Gridiv, &!input
      id, tstep, dt_since_start, Diagnose, sfr, &
      OHM_coef, OHM_threshSW, OHM_threshWD, &
      soilstore_id, SoilStoreCap, state_id, SnowUse, DiagQS, &
      HDD_id, MetForcingData_grid, Ts5mindata_ir, qf, qn1, &
      avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown, &
      bldgh, alb, emis, cpAnOHM, kkAnOHM, chAnOHM, EmissionsMethod, &
      Tair_av, qn1_av_prev, dqndt_prev, qn1_s_av_prev, dqnsdt_prev, &
      StoreDrainPrm, &!inout
      qn1_S, snowFrac, dataOutLineESTM, qs, &!output
      qn1_av_next, dqndt_next, qn1_s_av_next, dqnsdt_next, &
      deltaQi, a1, a2, a3)

      IMPLICIT NONE

      INTEGER, INTENT(in)  ::StorageHeatMethod
      INTEGER, INTENT(in)  ::OHMIncQF
      INTEGER, INTENT(in)  ::Gridiv
      INTEGER, INTENT(in)  ::id
      INTEGER, INTENT(in)  ::tstep ! time step [s]
      INTEGER, INTENT(in) ::dt_since_start  ! time since simulation starts [s]
      INTEGER, INTENT(in)  ::Diagnose
      ! INTEGER, INTENT(in)  ::nsh              ! number of timesteps in one hour
      INTEGER, INTENT(in)  ::SnowUse          ! option for snow related calculations
      INTEGER, INTENT(in)  ::DiagQS           ! diagnostic option
      INTEGER, INTENT(in)  :: EmissionsMethod !< AnthropHeat option [-]

      REAL(KIND(1d0)), INTENT(in)::OHM_coef(nsurf + 1, 4, 3)                 ! OHM coefficients
      REAL(KIND(1d0)), INTENT(in)::OHM_threshSW(nsurf + 1) ! OHM thresholds
      REAL(KIND(1d0)), INTENT(in)::OHM_threshWD(nsurf + 1) ! OHM thresholds
      REAL(KIND(1d0)), INTENT(in)::soilstore_id(nsurf)                ! soil moisture
      REAL(KIND(1d0)), INTENT(in)::SoilStoreCap(nsurf)             ! capacity of soil store
      REAL(KIND(1d0)), INTENT(in)::state_id(nsurf) ! wetness status

      REAL(KIND(1d0)), DIMENSION(12), INTENT(in)::HDD_id
      REAL(KIND(1d0)), INTENT(in)::qf
      REAL(KIND(1d0)), INTENT(in)::qn1
      REAL(KIND(1d0)), INTENT(in)::qs_obs
      REAL(KIND(1d0)), INTENT(in)::avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown
      REAL(KIND(1d0)), INTENT(in)::bldgh

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::alb  !< albedo [-]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::emis !< emissivity [-]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::cpAnOHM   !< heat capacity [J m-3 K-1]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::kkAnOHM   !< thermal conductivity [W m-1 K-1]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::chAnOHM   !< bulk transfer coef [J m-3 K-1]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::snowFrac

      REAL(KIND(1d0)), DIMENSION(:, :), INTENT(in)::MetForcingData_grid !< met forcing array of grid

      REAL(KIND(1d0)), DIMENSION(:), INTENT(in)::Ts5mindata_ir

      REAL(KIND(1d0)), INTENT(in)::Tair_av
      REAL(KIND(1d0)), INTENT(in)::qn1_av_prev
      REAL(KIND(1d0)), INTENT(out)::qn1_av_next
      REAL(KIND(1d0)), INTENT(in)::dqndt_prev  ! Rate of change of net radiation [W m-2 h-1] at t-1
      REAL(KIND(1d0)), INTENT(out)::dqndt_next  ! Rate of change of net radiation [W m-2 h-1] at t-1
      REAL(KIND(1d0)), INTENT(in)::qn1_s_av_prev
      REAL(KIND(1d0)), INTENT(out)::qn1_s_av_next
      REAL(KIND(1d0)), INTENT(in)::dqnsdt_prev ! Rate of change of net radiation [W m-2 h-1] at t-1
      REAL(KIND(1d0)), INTENT(out)::dqnsdt_next ! Rate of change of net radiation [W m-2 h-1] at t-1
      ! REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout)   ::qn1_store_grid
      ! REAL(KIND(1d0)),DIMENSION(nsh),INTENT(inout)   ::qn1_S_store_grid !< stored qn1 [W m-2]

      ! REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_av_store_grid
      ! REAL(KIND(1d0)),DIMENSION(2*nsh+1),INTENT(inout)::qn1_S_av_store_grid !< average net radiation over previous hour [W m-2]
      REAL(KIND(1d0)), DIMENSION(6, nsurf), INTENT(inout)::StoreDrainPrm

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::deltaQi ! storage heat flux of snow surfaces

      REAL(KIND(1d0)), DIMENSION(27), INTENT(out):: dataOutLineESTM
      REAL(KIND(1d0)), INTENT(out)::qn1_S
      REAL(KIND(1d0)), INTENT(out):: qs ! storage heat flux
      REAL(KIND(1d0)), INTENT(out):: a1 !< AnOHM coefficients of grid [-]
      REAL(KIND(1d0)), INTENT(out):: a2 !< AnOHM coefficients of grid [h]
      REAL(KIND(1d0)), INTENT(out):: a3 !< AnOHM coefficients of grid [W m-2]

      REAL(KIND(1d0))::Tair_mav_5d ! Tair_mav_5d=HDD(id-1,4) HDD at the begining of today (id-1)
      REAL(KIND(1d0))::qn1_use ! qn used in OHM calculations

      REAL(KIND(1d0)):: moist_surf(nsurf) !< non-dimensional surface wetness status (0-1) [-]

      ! initialise output variables
      !deltaQi = 0
      !snowFrac = 0
      !qn1_S = 0
      dataOutLineESTM = -999
      qs = -999
      a1 = -999
      a2 = -999
      a3 = -999

      ! calculate qn if qf should be included
      IF (OHMIncQF == 1) THEN
         qn1_use = qf + qn1
      ELSEIF (OHMIncQF == 0) THEN
         qn1_use = qn1
      ENDIF

      IF (StorageHeatMethod == 0) THEN           !Use observed QS
         qs = qs_obs

      ELSEIF (StorageHeatMethod == 1) THEN           !Use OHM to calculate QS
         Tair_mav_5d = HDD_id(10)
         IF (Diagnose == 1) WRITE (*, *) 'Calling OHM...'
         CALL OHM(qn1_use, qn1_av_prev, dqndt_prev, qn1_av_next, dqndt_next, &
                  qn1_S, qn1_s_av_prev, dqnsdt_prev, qn1_s_av_next, dqnsdt_next, &
                  tstep, dt_since_start, &
                  sfr, nsurf, &
                  Tair_mav_5d, &
                  OHM_coef, &
                  OHM_threshSW, OHM_threshWD, &
                  soilstore_id, SoilStoreCap, state_id, &
                  BldgSurf, WaterSurf, &
                  SnowUse, SnowFrac, &
                  DiagQS, &
                  a1, a2, a3, qs, deltaQi)

         ! use AnOHM to calculate QS, TS 14 Mar 2016
      ELSEIF (StorageHeatMethod == 3) THEN
         IF (Diagnose == 1) WRITE (*, *) 'Calling AnOHM...'
         ! CALL AnOHM(qn1_use,qn1_store_grid,qn1_av_store_grid,qf,&
         !      MetForcingData_grid,state_id/StoreDrainPrm(6,:),&
         !      alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&
         !      sfr,nsurf,nsh,EmissionsMethod,id,Gridiv,&
         !      a1,a2,a3,qs,deltaQi)
         moist_surf = state_id/StoreDrainPrm(6, :)
         CALL AnOHM( &
            tstep, dt_since_start, &
            qn1_use, qn1_av_prev, dqndt_prev, qf, &
            MetForcingData_grid, moist_surf, &
            alb, emis, cpAnOHM, kkAnOHM, chAnOHM, &! input
            sfr, nsurf, EmissionsMethod, id, Gridiv, &
            qn1_av_next, dqndt_next, &
            a1, a2, a3, qs, deltaQi)! output

         ! !Calculate QS using ESTM
      ELSEIF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
         !    !CALL ESTM(QSestm,iMB)
         IF (Diagnose == 1) WRITE (*, *) 'Calling ESTM...'
         CALL ESTM( &
            Gridiv, &!input
            tstep, &
            avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown, &
            bldgh, Ts5mindata_ir, &
            Tair_av, &
            dataOutLineESTM, QS)!output
         !    CALL ESTM(QSestm,Gridiv,ir)  ! iMB corrected to Gridiv, TS 09 Jun 2016
         !    QS=QSestm   ! Use ESTM qs
      ENDIF

   END SUBROUTINE SUEWS_cal_Qs
   !=======================================================================

   !==========================water balance================================
   SUBROUTINE SUEWS_cal_Water( &
      Diagnose, &!input
      snowUse, NonWaterFraction, addPipes, addImpervious, addVeg, addWaterBody, &
      state_id, soilstore_id, sfr, StoreDrainPrm, WaterDist, nsh_real, &
      drain_per_tstep, &  !output
      drain, AddWaterRunoff, &
      AdditionalWater, runoffPipes, runoff_per_interval, &
      AddWater, stateOld, soilstoreOld)

      IMPLICIT NONE
      ! INTEGER,PARAMETER :: nsurf=7! number of surface types
      ! INTEGER,PARAMETER ::WaterSurf = 7
      INTEGER, INTENT(in) ::Diagnose
      INTEGER, INTENT(in) ::snowUse

      REAL(KIND(1d0)), INTENT(in)::NonWaterFraction
      REAL(KIND(1d0)), INTENT(in)::addPipes
      REAL(KIND(1d0)), INTENT(in)::addImpervious
      REAL(KIND(1d0)), INTENT(in)::addVeg
      REAL(KIND(1d0)), INTENT(in)::addWaterBody
      REAL(KIND(1d0)), INTENT(in)::nsh_real !nsh cast as a real for use in calculations

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)          ::state_id
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)          ::soilstore_id
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)          ::sfr
      REAL(KIND(1d0)), DIMENSION(6, nsurf), INTENT(in)        ::StoreDrainPrm
      REAL(KIND(1d0)), DIMENSION(nsurf + 1, nsurf - 1), INTENT(in)::WaterDist

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: drain         !Drainage of surface type "is" [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: AddWaterRunoff!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: AddWater      !water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: stateOld
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out):: soilstoreOld

      REAL(KIND(1d0)), INTENT(out)::drain_per_tstep
      REAL(KIND(1d0)), INTENT(out)::AdditionalWater
      REAL(KIND(1d0)), INTENT(out)::runoffPipes
      REAL(KIND(1d0)), INTENT(out)::runoff_per_interval
      INTEGER:: is

      ! Retain previous surface state_id and soil moisture state_id
      stateOld = state_id     !state_id of each surface [mm] for the previous timestep
      soilstoreOld = soilstore_id !Soil moisture of each surface [mm] for the previous timestep

      !============= Grid-to-grid runoff =============
      ! Calculate additional water coming from other grids
      ! i.e. the variables addImpervious, addVeg, addWaterBody, addPipes
      !call RunoffFromGrid(GridFromFrac)  !!Need to code between-grid water transfer

      ! Sum water coming from other grids (these are expressed as depths over the whole surface)
      AdditionalWater = addPipes + addImpervious + addVeg + addWaterBody  ![mm]

      ! Initialise runoff in pipes
      runoffPipes = addPipes !Water flowing in pipes from other grids. QUESTION: No need for scaling?
      !! CHECK p_i
      runoff_per_interval = addPipes !pipe plor added to total runoff.

      !================== Drainage ===================
      ! Calculate drainage for each soil subsurface (excluding water body)
      IF (Diagnose == 1) WRITE (*, *) 'Calling Drainage...'

      IF (NonWaterFraction /= 0) THEN !Soil states only calculated if soil exists. LJ June 2017
         DO is = 1, nsurf - 1

            CALL drainage( &
               is, &! input:
               state_id(is), &
               StoreDrainPrm(6, is), &
               StoreDrainPrm(2, is), &
               StoreDrainPrm(3, is), &
               StoreDrainPrm(4, is), &
               nsh_real, &
               drain(is))! output

            ! !HCW added and changed to StoreDrainPrm(6,is) here 20 Feb 2015
            ! drain_per_tstep=drain_per_tstep+(drain(is)*sfr(is)/NonWaterFraction)   !No water body included
         ENDDO
         drain_per_tstep = DOT_PRODUCT(drain(1:nsurf - 1), sfr(1:nsurf - 1))/NonWaterFraction !No water body included
      ELSE
         drain(1:nsurf - 1) = 0
         drain_per_tstep = 0
      ENDIF

      drain(WaterSurf) = 0  ! Set drainage from water body to zero

      ! Distribute water within grid, according to WithinGridWaterDist matrix (Cols 1-7)
      IF (Diagnose == 1) WRITE (*, *) 'Calling ReDistributeWater...'
      ! CALL ReDistributeWater
      !Calculates AddWater(is)
      CALL ReDistributeWater( &
         snowUse, WaterDist, sfr, Drain, &! input:
         AddWaterRunoff, AddWater)! output

   END SUBROUTINE SUEWS_cal_Water
   !=======================================================================

   !===============initialize sensible heat flux============================
   SUBROUTINE SUEWS_init_QH( &
      qh_obs, avdens, avcp, h_mod, qn1, dectime, &!input
      H_init)!output

      IMPLICIT NONE
      REAL(KIND(1d0)), INTENT(in)::qh_obs
      REAL(KIND(1d0)), INTENT(in)::avdens
      REAL(KIND(1d0)), INTENT(in)::avcp
      REAL(KIND(1d0)), INTENT(in)::h_mod
      REAL(KIND(1d0)), INTENT(in)::qn1
      REAL(KIND(1d0)), INTENT(in)::dectime
      REAL(KIND(1d0)), INTENT(out)::H_init

      REAL(KIND(1d0)), PARAMETER::NAN = -999
      INTEGER, PARAMETER::notUsedI = -999

      ! Calculate kinematic heat flux (w'T') from sensible heat flux [W m-2] from observed data (if available) or LUMPS
      IF (qh_obs /= NAN) THEN   !if(qh_obs/=NAN) qh=qh_obs   !Commented out by HCW 04 Mar 2015
         H_init = qh_obs/(avdens*avcp)  !Use observed value
      ELSE
         IF (h_mod /= NAN) THEN
            H_init = h_mod/(avdens*avcp)   !Use LUMPS value
         ELSE
            H_init = (qn1*0.2)/(avdens*avcp)   !If LUMPS has had a problem, we still need a value
            CALL ErrorHint(38, 'LUMPS unable to calculate realistic value for H_mod.', h_mod, dectime, notUsedI)
         ENDIF
      ENDIF

   END SUBROUTINE SUEWS_init_QH
   !========================================================================

   !================latent heat flux and surface wetness===================
   ! TODO: optimise the structure of this function
   SUBROUTINE SUEWS_cal_QE( &
      Diagnose, &!input
      snowuse, &
      tstep, imin, it, EvapMethod, snowCalcSwitch, dayofWeek_id, CRWmin, CRWmax, &
      dectime, lvS_J_kg, lv_j_kg, avdens, avRh, Press_hPa, Temp_C, &
      RAsnow, psyc_hPa, avcp, sIce_hPa, &
      PervFraction, vegfraction, addimpervious, qn1_snowfree, qf, qs, vpd_hPa, s_hPa, &
      ResistSurf, RA, rb, snowdensmin, precip, PipeCapacity, RunoffToWater, &
      NonWaterFraction, wu_EveTr, wu_DecTr, wu_Grass, addVeg, addWaterBody, SnowLimPaved, SnowLimBldg, &
      SurfaceArea, FlowChange, drain, WetThresh, stateOld, mw_ind, SoilStoreCap, rainonsnow, &
      freezmelt, freezstate, freezstatevol, Qm_Melt, Qm_rain, Tsurf_ind, sfr, &
      StateLimit, AddWater, addwaterrunoff, StoreDrainPrm, SnowPackLimit, SnowProf_24hr, &
      runoff_per_interval, state_id, soilstore_id, SnowPack, snowFrac, SnowWater, &! inout:
      iceFrac, SnowDens, &
      runoffSnow, runoff, runoffSoil, chang, changSnow, &! output:
      snowDepth, SnowToSurf, ev_snow, SnowRemoval, &
      evap, rss_nsurf, p_mm, qe, state_per_tstep, NWstate_per_tstep, qeOut, &
      swe, ev, chSnow_per_interval, ev_per_tstep, qe_per_tstep, runoff_per_tstep, &
      surf_chang_per_tstep, runoffPipes, mwstore, runoffwaterbody, &
      runoffAGveg, runoffAGimpervious, runoffWaterBody_m3, runoffPipes_m3)

      IMPLICIT NONE

      INTEGER, INTENT(in) ::Diagnose
      INTEGER, INTENT(in) ::snowuse
      INTEGER, INTENT(in) ::tstep
      INTEGER, INTENT(in) ::imin
      INTEGER, INTENT(in) ::it
      INTEGER, INTENT(in) ::EvapMethod !Evaporation calculated according to Rutter (1) or Shuttleworth (2)

      INTEGER, DIMENSION(nsurf), INTENT(in)::snowCalcSwitch
      INTEGER, DIMENSION(3), INTENT(in)::dayofWeek_id

      REAL(KIND(1d0)), INTENT(in)::CRWmin
      REAL(KIND(1d0)), INTENT(in)::CRWmax
      REAL(KIND(1d0)), INTENT(in)::dectime
      REAL(KIND(1d0)), INTENT(in)::lvS_J_kg
      REAL(KIND(1d0)), INTENT(in)::lv_j_kg
      REAL(KIND(1d0)), INTENT(in)::avdens
      REAL(KIND(1d0)), INTENT(in)::avRh
      REAL(KIND(1d0)), INTENT(in)::Press_hPa
      REAL(KIND(1d0)), INTENT(in)::Temp_C
      REAL(KIND(1d0)), INTENT(in)::RAsnow
      REAL(KIND(1d0)), INTENT(in)::psyc_hPa
      REAL(KIND(1d0)), INTENT(in)::avcp
      REAL(KIND(1d0)), INTENT(in)::sIce_hPa
      REAL(KIND(1d0)), INTENT(in)::PervFraction
      REAL(KIND(1d0)), INTENT(in)::vegfraction
      REAL(KIND(1d0)), INTENT(in)::addimpervious
      REAL(KIND(1d0)), INTENT(in)::qn1_snowfree
      REAL(KIND(1d0)), INTENT(in)::qf
      REAL(KIND(1d0)), INTENT(in)::qs
      REAL(KIND(1d0)), INTENT(in)::vpd_hPa
      REAL(KIND(1d0)), INTENT(in)::s_hPa
      REAL(KIND(1d0)), INTENT(in)::ResistSurf
      REAL(KIND(1d0)), INTENT(in)::RA
      REAL(KIND(1d0)), INTENT(in)::rb
      REAL(KIND(1d0)), INTENT(in)::snowdensmin
      REAL(KIND(1d0)), INTENT(in)::precip
      REAL(KIND(1d0)), INTENT(in)::PipeCapacity
      REAL(KIND(1d0)), INTENT(in)::RunoffToWater
      REAL(KIND(1d0)), INTENT(in)::NonWaterFraction
      REAL(KIND(1d0)), INTENT(in)::wu_EveTr!Water use for evergreen trees/shrubs [mm]
      REAL(KIND(1d0)), INTENT(in)::wu_DecTr!Water use for deciduous trees/shrubs [mm]
      REAL(KIND(1d0)), INTENT(in)::wu_Grass!Water use for grass [mm]
      REAL(KIND(1d0)), INTENT(in)::addVeg!Water from vegetated surfaces of other grids [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(in)::addWaterBody!Water from water surface of other grids [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(in)::SnowLimPaved
      REAL(KIND(1d0)), INTENT(in)::SnowLimBldg
      REAL(KIND(1d0)), INTENT(in)::SurfaceArea
      REAL(KIND(1d0)), INTENT(in)::FlowChange

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::drain
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::WetThresh
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::stateOld
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::mw_ind
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::SoilStoreCap
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::rainonsnow
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::freezmelt
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::freezstate
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::freezstatevol
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::Qm_Melt
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::Qm_rain
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::Tsurf_ind
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::SnowPackLimit
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::StateLimit !Limit for state_id of each surface type [mm] (specified in input files)
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::AddWater
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::addwaterrunoff
      REAL(KIND(1d0)), DIMENSION(6, nsurf), INTENT(in)::StoreDrainPrm
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in):: SnowProf_24hr

      !Updated status: input and output
      REAL(KIND(1d0)), INTENT(inout)::runoff_per_interval! Total water transported to each grid for grid-to-grid connectivity

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::state_id
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::soilstore_id
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::SnowPack
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::snowFrac
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::SnowWater
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::iceFrac
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::SnowDens

      ! output:
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::runoffSnow !Initialize for runoff caused by snowmelting
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::runoff
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::runoffSoil
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::chang
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::changSnow
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::snowDepth
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::SnowToSurf
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::ev_snow
      REAL(KIND(1d0)), DIMENSION(2), INTENT(out)::SnowRemoval
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::evap
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::rss_nsurf

      REAL(KIND(1d0)), INTENT(out)::p_mm!Inputs to surface water balance
      ! REAL(KIND(1d0)),INTENT(out)::rss
      REAL(KIND(1d0)), INTENT(out)::qe ! latent heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(out)::state_per_tstep
      REAL(KIND(1d0)), INTENT(out)::NWstate_per_tstep
      REAL(KIND(1d0)), INTENT(out)::qeOut
      REAL(KIND(1d0)), INTENT(out)::swe
      REAL(KIND(1d0)), INTENT(out)::ev
      REAL(KIND(1d0)), INTENT(out)::chSnow_per_interval
      REAL(KIND(1d0)), INTENT(out)::ev_per_tstep
      REAL(KIND(1d0)), INTENT(out)::qe_per_tstep
      REAL(KIND(1d0)), INTENT(out)::runoff_per_tstep
      REAL(KIND(1d0)), INTENT(out)::surf_chang_per_tstep
      REAL(KIND(1d0)), INTENT(out)::runoffPipes
      REAL(KIND(1d0)), INTENT(out)::mwstore
      REAL(KIND(1d0)), INTENT(out)::runoffwaterbody
      REAL(KIND(1d0)), INTENT(out)::runoffWaterBody_m3
      REAL(KIND(1d0)), INTENT(out)::runoffPipes_m3
      REAL(KIND(1d0)), INTENT(out)::runoffAGveg
      REAL(KIND(1d0)), INTENT(out)::runoffAGimpervious

      ! local:
      INTEGER:: is

      REAL(KIND(1d0)), DIMENSION(2)    ::SurplusEvap        !Surplus for evaporation in 5 min timestep
      REAL(KIND(1d0))::surplusWaterBody
      REAL(KIND(1d0))::pin!Rain per time interval
      ! REAL(KIND(1d0))::sae
      ! REAL(KIND(1d0))::vdrc
      ! REAL(KIND(1d0))::sp
      ! REAL(KIND(1d0))::numPM
      REAL(KIND(1d0))::qn_e
      REAL(KIND(1d0))::tlv
      REAL(KIND(1d0))::runoffAGimpervious_m3
      REAL(KIND(1d0))::runoffAGveg_m3
      REAL(KIND(1d0))::nsh_real
      REAL(KIND(1d0))::tstep_real
      REAL(KIND(1d0))::ev_tot
      REAL(KIND(1d0))::qe_tot
      REAL(KIND(1d0))::surf_chang_tot
      REAL(KIND(1d0))::runoff_tot
      REAL(KIND(1d0))::chSnow_tot

      REAL(KIND(1d0)), DIMENSION(7)::capStore ! current storage capacity [mm]

      tstep_real = tstep*1.d0
      nsh_real = 3600/tstep_real

      capStore = 0 !initialise capStore

      tlv = lv_J_kg/tstep_real !Latent heat of vapourisation per timestep

      pin = MAX(0., Precip)!Initiate rain data [mm]

      ! Initialize the output variables
      qe = 0
      ev = 0
      qe_tot = 0
      ev_tot = 0
      swe = 0
      ev_snow = 0
      ev_per_tstep = 0
      qe_per_tstep = 0
      surf_chang_per_tstep = 0
      runoff_per_tstep = 0
      state_per_tstep = 0
      NWstate_per_tstep = 0
      qeOut = 0
      runoffwaterbody = 0
      chSnow_per_interval = 0
      mwstore = 0
      runoffAGveg = 0
      runoffAGimpervious = 0
      surplusWaterBody = 0
      runoffSoil = 0
      runoff = 0
      chang = 0
      SurplusEvap = 0
      SnowRemoval = 0

      ! net available energy for evaporation
      qn_e = qn1_snowfree + qf - qs ! qn1 changed to qn1_snowfree, lj in May 2013

      IF (Diagnose == 1) WRITE (*, *) 'Calling evap_SUEWS and SoilStore...'
      DO is = 1, nsurf   !For each surface in turn
         IF (snowCalcSwitch(is) == 1) THEN ! snow calculation
            IF (sfr(is) /= 0) THEN
               ! IF (Diagnose == 1) WRITE (*, *) 'Calling SnowCalc...'
               CALL SnowCalc( &
                  tstep, imin, it, dectime, is, &!input
                  EvapMethod, CRWmin, CRWmax, nsh_real, lvS_J_kg, lv_j_kg, avdens, &
                  avRh, Press_hPa, Temp_C, RAsnow, psyc_hPa, avcp, sIce_hPa, &
                  PervFraction, vegfraction, addimpervious, &
                  vpd_hPa, qn_e, s_hPa, ResistSurf, RA, rb, tlv, snowdensmin, SnowProf_24hr, precip, &
                  PipeCapacity, RunoffToWater, &
                  addVeg, SnowLimPaved, SnowLimBldg, FlowChange, drain, &
                  WetThresh, stateOld, mw_ind, SoilStoreCap, rainonsnow, &
                  freezmelt, freezstate, freezstatevol, &
                  Qm_Melt, Qm_rain, Tsurf_ind, sfr, dayofWeek_id, StoreDrainPrm, SnowPackLimit, &
                  AddWater, addwaterrunoff, &
                  soilstore_id, SnowPack, SurplusEvap, &!inout
                  snowFrac, SnowWater, iceFrac, SnowDens, &
                  runoffAGimpervious, runoffAGveg, surplusWaterBody, &
                  rss_nsurf, runoffSnow, & ! output
                  runoff, runoffSoil, chang, changSnow, SnowToSurf, state_id, ev_snow, &
                  SnowDepth, SnowRemoval, swe, ev, chSnow_tot, &
                  ev_tot, qe_tot, runoff_tot, surf_chang_tot, &
                  runoffPipes, mwstore, runoffwaterbody)

               !Actual updates here as xx_tstep variables not taken as input to snowcalc
               ev_per_tstep = ev_per_tstep + ev_tot
               qe_per_tstep = qe_per_tstep + qe_tot
               runoff_per_tstep = runoff_per_tstep + runoff_tot
               surf_chang_per_tstep = surf_chang_per_tstep + surf_chang_tot
               chSnow_per_interval = chSnow_per_interval + chSnow_tot
            ELSE
               snowFrac(is) = 0
               SnowDens(is) = 0
               SnowPack(is) = 0
            ENDIF
         ELSE ! snow-free calculation

            capStore(is) = StoreDrainPrm(6, is)
            !Calculates ev [mm]
            CALL Evap_SUEWS( &
               EvapMethod, state_id(is), WetThresh(is), capStore(is), &!input
               vpd_hPa, avdens, avcp, qn_e, s_hPa, psyc_hPa, ResistSurf, RA, rb, tlv, &
               rss_nsurf(is), ev, qe) !output

            !Surface water balance and soil store updates (can modify ev, updates state_id)
            CALL soilstore( &
               is, sfr, PipeCapacity, RunoffToWater, pin, & ! input:
               wu_EveTr, wu_DecTr, wu_Grass, drain, AddWater, addImpervious, nsh_real, stateOld, AddWaterRunoff, &
               PervFraction, addVeg, SoilStoreCap, addWaterBody, FlowChange, StateLimit, runoffAGimpervious, surplusWaterBody, &
               runoffAGveg, runoffPipes, ev, soilstore_id, SurplusEvap, runoffWaterBody, &
               p_mm, chang, runoff, state_id)!output:

            !Store ev for each surface
            evap(is) = ev

            ! Sum evaporation from different surfaces to find total evaporation [mm]
            ev_per_tstep = ev_per_tstep + evap(is)*sfr(is)

            ! Sum latent heat flux from different surfaces to find total latent heat flux
            qe_per_tstep = qe_per_tstep + qe*sfr(is)

            ! Sum change from different surfaces to find total change to surface state_id
            surf_chang_per_tstep = surf_chang_per_tstep + (state_id(is) - stateOld(is))*sfr(is)

            ! Sum runoff from different surfaces to find total runoff
            runoff_per_tstep = runoff_per_tstep + runoff(is)*sfr(is)

            ! Calculate total state_id (including water body)
            state_per_tstep = state_per_tstep + state_id(is)*sfr(is)

            ! sum The total runoff from the area !!Check (HCW)
            runoff_per_interval = runoff_per_interval + (runoff(is)*sfr(is))

            IF (NonWaterFraction /= 0 .AND. is /= WaterSurf) THEN
               NWstate_per_tstep = NWstate_per_tstep + (state_id(is)*sfr(is)/NonWaterFraction)
            ENDIF

            ChangSnow(is) = 0
            runoffSnow(is) = 0

         ENDIF
      ENDDO  !end loop over surfaces

      qeOut = qe_per_tstep

      ! Calculate volume of water that will move between grids
      ! Volume [m3] = Depth relative to whole area [mm] / 1000 [mm m-1] * SurfaceArea [m2]
      ! Need to use these volumes when converting back to addImpervious, AddVeg and AddWater
      runoffAGimpervious_m3 = runoffAGimpervious/1000*SurfaceArea
      runoffAGveg_m3 = runoffAGveg/1000*SurfaceArea
      runoffWaterBody_m3 = runoffWaterBody/1000*SurfaceArea
      runoffPipes_m3 = runoffPipes/1000*SurfaceArea

   END SUBROUTINE SUEWS_cal_QE
   !========================================================================

   !===============sensible heat flux======================================
   SUBROUTINE SUEWS_cal_QH( &
      QHMethod, &!input
      qn1, qf, QmRain, qeOut, qs, QmFreez, qm, avdens, avcp, tsurf, Temp_C, RA, &
      qh, qh_residual, qh_resist)!output
      IMPLICIT NONE

      INTEGER, INTENT(in) :: QHMethod ! option for QH calculation: 1, residual; 2, resistance-based

      REAL(KIND(1d0)), INTENT(in)::qn1
      REAL(KIND(1d0)), INTENT(in)::qf
      REAL(KIND(1d0)), INTENT(in)::QmRain
      REAL(KIND(1d0)), INTENT(in)::qeOut
      REAL(KIND(1d0)), INTENT(in)::qs
      REAL(KIND(1d0)), INTENT(in)::QmFreez
      REAL(KIND(1d0)), INTENT(in)::qm
      REAL(KIND(1d0)), INTENT(in)::avdens
      REAL(KIND(1d0)), INTENT(in)::avcp
      REAL(KIND(1d0)), INTENT(in)::tsurf
      REAL(KIND(1d0)), INTENT(in)::Temp_C
      REAL(KIND(1d0)), INTENT(in)::RA

      REAL(KIND(1d0)), INTENT(out)::qh
      REAL(KIND(1d0)), INTENT(out)::qh_resist
      REAL(KIND(1d0)), INTENT(out)::qh_residual

      REAL(KIND(1d0)), PARAMETER::NAN = -999

      ! Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
      qh_residual = (qn1 + qf + QmRain) - (qeOut + qs + Qm + QmFreez)     !qh=(qn1+qf+QmRain+QmFreez)-(qeOut+qs+Qm)

      ! ! Calculate QH using resistance method (for testing HCW 06 Jul 2016)
      ! Aerodynamic-Resistance-based method
      IF (RA /= 0) THEN
         qh_resist = avdens*avcp*(tsurf - Temp_C)/RA
      ELSE
         qh_resist = NAN
      ENDIF

      ! choose output QH
      SELECT CASE (QHMethod)
      CASE (1)
         qh = qh_residual
      CASE (2)
         qh = qh_resist
      END SELECT

   END SUBROUTINE SUEWS_cal_QH
   !========================================================================

   !===============Resistance Calculations=======================
   SUBROUTINE SUEWS_cal_Resistance( &
      StabilityMethod, &!input:
      Diagnose, AerodynamicResistanceMethod, RoughLenHeatMethod, snowUse, &
      id, it, gsModel, SMDMethod, &
      qh_obs, avdens, avcp, h_mod, qn1, dectime, zzd, z0m, zdm, &
      avU1, Temp_C, VegFraction, &
      avkdn, Kmax, G1, G2, G3, G4, G5, G6, S1, S2, TH, TL, dq, &
      xsmd, vsmd, MaxConductance, LAIMax, LAI_id, snowFrac, sfr, &
      UStar, TStar, L_mod, &!output
      zL, gsc, ResistSurf, RA, RAsnow, rb)

      IMPLICIT NONE

      INTEGER, INTENT(in)::StabilityMethod
      INTEGER, INTENT(in)::Diagnose
      INTEGER, INTENT(in)::AerodynamicResistanceMethod
      INTEGER, INTENT(in)::RoughLenHeatMethod
      INTEGER, INTENT(in)::snowUse
      INTEGER, INTENT(in)::id
      INTEGER, INTENT(in)::it       !time: day of year and hour
      INTEGER, INTENT(in)::gsModel  !Choice of gs parameterisation (1 = Ja11, 2 = Wa16)
      INTEGER, INTENT(in)::SMDMethod!Method of measured soil moisture

      REAL(KIND(1d0)), INTENT(in)::qh_obs
      REAL(KIND(1d0)), INTENT(in)::avdens
      REAL(KIND(1d0)), INTENT(in)::avcp
      REAL(KIND(1d0)), INTENT(in)::h_mod
      REAL(KIND(1d0)), INTENT(in)::qn1
      REAL(KIND(1d0)), INTENT(in)::dectime    !Decimal time
      REAL(KIND(1d0)), INTENT(in)::zzd        !Active measurement height (meas. height-displac. height)
      REAL(KIND(1d0)), INTENT(in)::z0m        !Aerodynamic roughness length
      REAL(KIND(1d0)), INTENT(in)::zdm        !Displacement height
      REAL(KIND(1d0)), INTENT(in)::avU1       !Average wind speed
      REAL(KIND(1d0)), INTENT(in)::Temp_C     !Air temperature
      REAL(KIND(1d0)), INTENT(in)::VegFraction!Fraction of vegetation
      REAL(KIND(1d0)), INTENT(in)::avkdn      !Average downwelling shortwave radiation
      REAL(KIND(1d0)), INTENT(in)::Kmax       !Annual maximum hourly solar radiation
      REAL(KIND(1d0)), INTENT(in)::G1         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::G2         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::G3         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::G4         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::G5         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::G6         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::S1         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::S2         !Fitted parameters related to surface res. calculations
      REAL(KIND(1d0)), INTENT(in)::TH         !Maximum temperature limit
      REAL(KIND(1d0)), INTENT(in)::TL         !Minimum temperature limit
      REAL(KIND(1d0)), INTENT(in)::dq         !Specific humidity deficit
      REAL(KIND(1d0)), INTENT(in)::xsmd       !Measured soil moisture deficit
      REAL(KIND(1d0)), INTENT(in)::vsmd       !Soil moisture deficit for vegetated surfaces only (QUESTION: what about BSoil?)

      REAL(KIND(1d0)), DIMENSION(3), INTENT(in) ::MaxConductance!Max conductance [mm s-1]
      REAL(KIND(1d0)), DIMENSION(3), INTENT(in) ::LAIMax        !Max LAI [m2 m-2]
      REAL(KIND(1d0)), DIMENSION(3), INTENT(in) ::LAI_id        !=LAI_id(id-1,:), LAI for each veg surface [m2 m-2]

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::snowFrac      !Surface fraction of snow cover
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr           !Surface fractions [-]

      REAL(KIND(1d0)), INTENT(out)::TStar     !T*
      REAL(KIND(1d0)), INTENT(out)  ::UStar     !Friction velocity
      ! REAL(KIND(1d0)),INTENT(out)::psim      !Stability function of momentum
      REAL(KIND(1d0)), INTENT(out)  ::zL        !
      REAL(KIND(1d0)), INTENT(out)  ::gsc       !Surface Layer Conductance
      REAL(KIND(1d0)), INTENT(out)  ::ResistSurf!Surface resistance
      REAL(KIND(1d0)), INTENT(out)  ::RA        !Aerodynamic resistance [s m^-1]
      REAL(KIND(1d0)), INTENT(out)  ::RAsnow    !Aerodynamic resistance for snow [s m^-1]
      REAL(KIND(1d0)), INTENT(out)  ::rb        !boundary layer resistance shuttleworth
      REAL(KIND(1d0)), INTENT(out)  ::L_mod     !Obukhov length
      REAL(KIND(1d0))              ::gfunc     !gdq*gtemp*gs*gq for photosynthesis calculations
      REAL(KIND(1d0))              ::H_init    !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity

      ! Get first estimate of sensible heat flux. Modified by HCW 26 Feb 2015
      CALL SUEWS_init_QH( &
         qh_obs, avdens, avcp, h_mod, qn1, dectime, &
         H_init)

      IF (Diagnose == 1) WRITE (*, *) 'Calling STAB_lumps...'
      !u* and Obukhov length out
      CALL STAB_lumps( &
         StabilityMethod, &  ! input
         dectime, & !Decimal time
         zzd, &     !Active measurement height (meas. height-displac. height)
         z0m, &     !Aerodynamic roughness length
         zdm, &     !Displacement height
         avU1, &    !Average wind speed
         Temp_C, &  !Air temperature
         H_init, & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
         L_mod, &! output: !Obukhov length
         TStar, & !T*
         UStar, & !Friction velocity
         zL)!Stability scale

      IF (Diagnose == 1) WRITE (*, *) 'Calling AerodynamicResistance...'
      CALL AerodynamicResistance( &
         ZZD, &! input:
         z0m, &
         AVU1, &
         L_mod, &
         UStar, &
         VegFraction, &
         AerodynamicResistanceMethod, &
         StabilityMethod, &
         RoughLenHeatMethod, &
         RA) ! output:

      IF (snowUse == 1) THEN
         IF (Diagnose == 1) WRITE (*, *) 'Calling AerodynamicResistance for snow...'
         CALL AerodynamicResistance( &
            ZZD, &! input:
            z0m, &
            AVU1, &
            L_mod, &
            UStar, &
            VegFraction, &
            AerodynamicResistanceMethod, &
            StabilityMethod, &
            3, &
            RAsnow)     ! output:
      ENDIF

      IF (Diagnose == 1) WRITE (*, *) 'Calling SurfaceResistance...'
      ! CALL SurfaceResistance(id,it)   !qsc and surface resistance out
      CALL SurfaceResistance( &
         id, it, &! input:
         SMDMethod, snowFrac, sfr, avkdn, Temp_C, dq, xsmd, vsmd, MaxConductance, &
         LAIMax, LAI_id, gsModel, Kmax, &
         G1, G2, G3, G4, G5, G6, TH, TL, S1, S2, &
         gfunc, gsc, ResistSurf)! output:

      IF (Diagnose == 1) WRITE (*, *) 'Calling BoundaryLayerResistance...'
      CALL BoundaryLayerResistance( &
         zzd, &! input:     !Active measurement height (meas. height-displac. height)
         z0m, &     !Aerodynamic roughness length
         avU1, &    !Average wind speed
         UStar, &  ! input/output:
         rb)  ! output:

   END SUBROUTINE SUEWS_cal_Resistance
   !========================================================================

   !==============Update output arrays=========================
   SUBROUTINE SUEWS_update_outputLine( &
      AdditionalWater, alb, avkdn, avU10_ms, azimuth, &!input
      chSnow_per_interval, dectime, &
      drain_per_tstep, E_mod, ev_per_tstep, ext_wu, Fc, Fc_build, fcld, &
      Fc_metab, Fc_photo, Fc_respi, Fc_point, Fc_traff, FlowChange, &
      h_mod, id, imin, int_wu, it, iy, &
      kup, LAI_id, ldown, l_mod, lup, mwh, &
      MwStore, &
      nsh_real, NWstate_per_tstep, Precip, q2_gkg, &
      qeOut, qf, qh, qh_resist, Qm, QmFreez, &
      QmRain, qn1, qn1_S, qn1_snowfree, qs, RA, &
      resistsurf, RH2, runoffAGimpervious, runoffAGveg, &
      runoff_per_tstep, runoffPipes, runoffSoil_per_tstep, &
      runoffWaterBody, sfr, smd, smd_nsurf, SnowAlb, SnowRemoval, &
      state_id, state_per_tstep, surf_chang_per_tstep, swe, t2_C, tskin_C, &
      tot_chang_per_tstep, tsurf, UStar, wu_DecTr, &
      wu_EveTr, wu_Grass, z0m, zdm, zenith_deg, &
      datetimeLine, dataOutLineSUEWS)!output
      IMPLICIT NONE

      REAL(KIND(1d0)), PARAMETER :: NAN = -999
      INTEGER, INTENT(in) :: iy
      ! INTEGER,INTENT(in) :: iy_prev_t
      INTEGER, INTENT(in) :: id
      ! INTEGER,INTENT(in) :: id_prev_t
      INTEGER, INTENT(in) :: it
      INTEGER, INTENT(in) :: imin
      !  INTEGER, INTENT(in) :: Gridiv
      REAL(KIND(1d0)), INTENT(in) :: AdditionalWater
      REAL(KIND(1d0)), INTENT(in) :: alb(nsurf)
      REAL(KIND(1d0)), INTENT(in) :: avkdn
      REAL(KIND(1d0)), INTENT(in) :: avU10_ms
      REAL(KIND(1d0)), INTENT(in) :: azimuth
      REAL(KIND(1d0)), INTENT(in) :: chSnow_per_interval
      REAL(KIND(1d0)), INTENT(in) :: dectime
      REAL(KIND(1d0)), INTENT(in) :: drain_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: E_mod
      REAL(KIND(1d0)), INTENT(in) :: ev_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: ext_wu
      REAL(KIND(1d0)), INTENT(in) :: Fc
      REAL(KIND(1d0)), INTENT(in) :: Fc_build
      REAL(KIND(1d0)), INTENT(in) :: Fc_metab
      REAL(KIND(1d0)), INTENT(in) :: Fc_photo
      REAL(KIND(1d0)), INTENT(in) :: Fc_respi
      REAL(KIND(1d0)), INTENT(in) :: Fc_point
      REAL(KIND(1d0)), INTENT(in) :: Fc_traff
      REAL(KIND(1d0)), INTENT(in) :: fcld
      REAL(KIND(1d0)), INTENT(in) :: FlowChange
      REAL(KIND(1d0)), INTENT(in) :: h_mod
      REAL(KIND(1d0)), INTENT(in) :: int_wu
      REAL(KIND(1d0)), INTENT(in) :: kup
      REAL(KIND(1d0)), INTENT(in) :: l_mod
      REAL(KIND(1d0)), INTENT(in) :: LAI_id(nvegsurf)
      REAL(KIND(1d0)), INTENT(in) :: ldown
      REAL(KIND(1d0)), INTENT(in) :: lup
      REAL(KIND(1d0)), INTENT(in) :: mwh
      REAL(KIND(1d0)), INTENT(in) :: MwStore
      REAL(KIND(1d0)), INTENT(in) :: nsh_real
      REAL(KIND(1d0)), INTENT(in) :: NWstate_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: Precip
      REAL(KIND(1d0)), INTENT(in) :: q2_gkg
      REAL(KIND(1d0)), INTENT(in) :: qeOut
      REAL(KIND(1d0)), INTENT(in) :: qf
      REAL(KIND(1d0)), INTENT(in) :: qh
      REAL(KIND(1d0)), INTENT(in) :: qh_resist
      REAL(KIND(1d0)), INTENT(in) :: Qm
      REAL(KIND(1d0)), INTENT(in) :: QmFreez
      REAL(KIND(1d0)), INTENT(in) :: QmRain
      REAL(KIND(1d0)), INTENT(in) :: qn1
      REAL(KIND(1d0)), INTENT(in) :: qn1_S
      REAL(KIND(1d0)), INTENT(in) :: qn1_snowfree
      REAL(KIND(1d0)), INTENT(in) :: qs
      REAL(KIND(1d0)), INTENT(in) :: RA
      REAL(KIND(1d0)), INTENT(in) :: resistsurf
      REAL(KIND(1d0)), INTENT(in) :: RH2
      REAL(KIND(1d0)), INTENT(in) :: runoff_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: runoffAGimpervious
      REAL(KIND(1d0)), INTENT(in) :: runoffAGveg
      REAL(KIND(1d0)), INTENT(in) :: runoffPipes
      REAL(KIND(1d0)), INTENT(in) :: runoffSoil_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: runoffWaterBody
      REAL(KIND(1d0)), INTENT(in) :: sfr(nsurf)
      REAL(KIND(1d0)), INTENT(in) :: smd
      REAL(KIND(1d0)), INTENT(in) :: smd_nsurf(nsurf)
      REAL(KIND(1d0)), INTENT(in) :: SnowAlb
      REAL(KIND(1d0)), INTENT(in) :: SnowRemoval(2)
      REAL(KIND(1d0)), INTENT(in) :: state_id(nsurf)
      REAL(KIND(1d0)), INTENT(in) :: state_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: surf_chang_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: swe
      REAL(KIND(1d0)), INTENT(in) :: t2_C
      REAL(KIND(1d0)), INTENT(in) :: tskin_C
      REAL(KIND(1d0)), INTENT(in) :: tot_chang_per_tstep
      REAL(KIND(1d0)), INTENT(in) :: tsurf
      REAL(KIND(1d0)), INTENT(in) :: UStar
      REAL(KIND(1d0)), INTENT(in) :: wu_DecTr
      REAL(KIND(1d0)), INTENT(in) :: wu_EveTr
      REAL(KIND(1d0)), INTENT(in) :: wu_Grass
      REAL(KIND(1d0)), INTENT(in) :: z0m
      REAL(KIND(1d0)), INTENT(in) :: zdm
      REAL(KIND(1d0)), INTENT(in) :: zenith_deg

      REAL(KIND(1D0)), DIMENSION(5), INTENT(OUT)::datetimeLine
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutSUEWS - 5), INTENT(out) :: dataOutLineSUEWS
      ! REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutSnow-5),INTENT(out) :: dataOutLineSnow
      ! REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutESTM-5),INTENT(out) :: dataOutLineESTM
      ! INTEGER:: is
      REAL(KIND(1d0)):: LAI_wt
      REAL(KIND(1d0)):: RH2_pct ! RH2 in percentage

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
      state_x = UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr < 0.00001), field=state_id)
      smd_nsurf_x = UNPACK(SPREAD(NAN, dim=1, ncopies=SIZE(sfr)), mask=(sfr < 0.00001), field=smd_nsurf)

      ResistSurf_x = MIN(9999., ResistSurf)

      l_mod_x = MAX(MIN(9999., l_mod), -9999.)

      ! Calculate areally-weighted LAI
      ! IF(iy == (iy_prev_t  +1) .AND. (id-1) == 0) THEN   !Check for start of next year and avoid using LAI(id-1) as this is at the start of the year
      !    LAI_wt=DOT_PRODUCT(LAI(id_prev_t,:),sfr(1+2:nvegsurf+2))
      ! ELSE
      !    LAI_wt=DOT_PRODUCT(LAI(id-1,:),sfr(1+2:nvegsurf+2))
      ! ENDIF

      LAI_wt = DOT_PRODUCT(LAI_id(:), sfr(1 + 2:nvegsurf + 2))

      ! Calculate areally-weighted albedo
      bulkalbedo = DOT_PRODUCT(alb, sfr)

      ! convert RH2 to a percentage form
      RH2_pct = RH2*100.0

      ! NB: this part needs to be reconsidered for calculation logic. TS, 27 Sep 2018
      ! TODO: this part should be reconnected to an improved CBL interface. TS 10 Jun 2018
      ! uncomment the CBL.HWJ，21 Aug 2019
      ! Save qh and qe for CBL in next iteration
      ! IF(Qh_choice==1) THEN   !use QH and QE from SUEWS
      !     qhforCBL(Gridiv) = qh
      !     qeforCBL(Gridiv) = qeOut
      !  ELSEIF(Qh_choice==2)THEN   !use QH and QE from LUMPS
      !     qhforCBL(Gridiv) = h_mod
      !     qeforCBL(Gridiv) = e_mod
      !  ELSEIF(qh_choice==3)THEN  !use QH and QE from OBS
      !     qhforCBL(Gridiv) = qh_obs
      !     qeforCBL(Gridiv) = qe_obs
      !     IF(qh_obs<-900.OR.qe_obs<-900)THEN  ! observed data has a problem
      !        CALL ErrorHint(22,'Unrealistic observed qh or qe_value.',qh_obs,qe_obs,qh_choice)
      !     ENDIF
      ! ENDIF

      !====================== update output line ==============================
      ! date & time:
      datetimeLine = [ &
                     REAL(iy, KIND(1D0)), REAL(id, KIND(1D0)), &
                     REAL(it, KIND(1D0)), REAL(imin, KIND(1D0)), dectime]
      !Define the overall output matrix to be printed out step by step
      dataOutLineSUEWS = [ &
                         avkdn, kup, ldown, lup, tsurf, &
                         qn1, qf, qs, qh, qeOut, &
                         h_mod, e_mod, qh_resist, &
                         precip, ext_wu, ev_per_tstep, runoff_per_tstep, tot_chang_per_tstep, &
                         surf_chang_per_tstep, state_per_tstep, NWstate_per_tstep, drain_per_tstep, smd, &
                         FlowChange/nsh_real, AdditionalWater, &
                         runoffSoil_per_tstep, runoffPipes, runoffAGimpervious, runoffAGveg, runoffWaterBody, &
                         int_wu, wu_EveTr, wu_DecTr, wu_Grass, &
                         smd_nsurf_x(1:nsurf - 1), &
                         state_x(1:nsurf), &
                         zenith_deg, azimuth, bulkalbedo, Fcld, &
                         LAI_wt, z0m, zdm, &
                         UStar, l_mod, RA, ResistSurf, &
                         Fc, &
                         Fc_photo, Fc_respi, Fc_metab, Fc_traff, Fc_build, Fc_point, &
                         qn1_snowfree, qn1_S, SnowAlb, &
                         Qm, QmFreez, QmRain, swe, mwh, MwStore, chSnow_per_interval, &
                         SnowRemoval(1:2), &
                         tskin_C, t2_C, q2_gkg, avU10_ms, RH2_pct & ! surface-level diagonostics
                         ]
      ! set invalid values to NAN
      ! dataOutLineSUEWS = set_nan(dataOutLineSUEWS)

      !====================update output line end==============================

   END SUBROUTINE SUEWS_update_outputLine
   !========================================================================

   !==============Update output arrays=========================
   SUBROUTINE SUEWS_update_output( &
      SnowUse, storageheatmethod, &!input
      ReadLinesMetdata, NumberOfGrids, &
      ir, gridiv, &
      datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, dataoutLineURSL, &!input
      dataOutSUEWS, dataOutSnow, dataOutESTM, dataOutRSL)!inout
      IMPLICIT NONE

      INTEGER, INTENT(in) ::ReadLinesMetdata
      INTEGER, INTENT(in) ::NumberOfGrids
      INTEGER, INTENT(in) ::Gridiv
      INTEGER, INTENT(in) ::SnowUse
      INTEGER, INTENT(in) ::storageheatmethod
      INTEGER, INTENT(in) ::ir

      REAL(KIND(1d0)), DIMENSION(5), INTENT(in) :: datetimeLine
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutSUEWS - 5), INTENT(in) :: dataOutLineSUEWS
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutESTM - 5), INTENT(in) :: dataOutLineESTM
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutSnow - 5), INTENT(in) :: dataOutLineSnow
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutRSL - 5), INTENT(in) :: dataoutLineURSL

      REAL(KIND(1d0)), INTENT(inout) :: dataOutSUEWS(ReadLinesMetdata, ncolumnsDataOutSUEWS, NumberOfGrids)
      REAL(KIND(1d0)), INTENT(inout) :: dataOutSnow(ReadLinesMetdata, ncolumnsDataOutSnow, NumberOfGrids)
      REAL(KIND(1d0)), INTENT(inout) :: dataOutESTM(ReadLinesMetdata, ncolumnsDataOutESTM, NumberOfGrids)
      REAL(KIND(1d0)), INTENT(inout) :: dataOutRSL(ReadLinesMetdata, ncolumnsDataOutRSL, NumberOfGrids)

      !====================== update output arrays ==============================
      !Define the overall output matrix to be printed out step by step
      dataOutSUEWS(ir, 1:ncolumnsDataOutSUEWS, Gridiv) = [datetimeLine, set_nan(dataOutLineSUEWS)]
      dataOutRSL(ir, 1:ncolumnsDataOutRSL, Gridiv) = [datetimeLine, set_nan(dataoutLineURSL)]
      ! ! set invalid values to NAN
      ! dataOutSUEWS(ir,6:ncolumnsDataOutSUEWS,Gridiv)=set_nan(dataOutSUEWS(ir,6:ncolumnsDataOutSUEWS,Gridiv))

      IF (snowUse == 1) THEN
         dataOutSnow(ir, 1:ncolumnsDataOutSnow, Gridiv) = [datetimeLine, set_nan(dataOutLineSnow)]
      END IF

      IF (storageheatmethod == 4) THEN
         dataOutESTM(ir, 1:ncolumnsDataOutESTM, Gridiv) = [datetimeLine, set_nan(dataOutLineESTM)]
      END IF

      !====================update output arrays end==============================

   END SUBROUTINE SUEWS_update_output

   !========================================================================
   SUBROUTINE SUEWS_cal_Diagnostics( &
      dectime, &!input
      avU1, Temp_C, avRH, Press_hPa, &
      qh, qe, &
      VegFraction, zMeas, z0m, zdm, RA, avdens, avcp, lv_J_kg, tstep_real, &
      RoughLenHeatMethod, StabilityMethod, &
      avU10_ms, t2_C, q2_gkg, tskin_C, RH2)!output
      ! TS 03 Aug 2018: added limit on q2 by restricting RH2_max to 100%
      ! TS 31 Jul 2018: removed dependence on surface variables (Tsurf, qsat)
      ! TS 26 Jul 2018: improved the calculation logic
      ! TS 05 Sep 2017: improved interface
      ! TS 20 May 2017: calculate surface-level diagonostics
      IMPLICIT NONE
      REAL(KIND(1d0)), INTENT(in) ::dectime
      REAL(KIND(1d0)), INTENT(in) ::avU1, Temp_C, avRH
      REAL(KIND(1d0)), INTENT(in) ::qh
      REAL(KIND(1d0)), INTENT(in) ::Press_hPa, qe
      REAL(KIND(1d0)), INTENT(in) :: VegFraction, z0m, RA, avdens, avcp, lv_J_kg, tstep_real
      REAL(KIND(1d0)), INTENT(in) :: zMeas! height for measurement
      REAL(KIND(1d0)), INTENT(in) :: zdm ! displacement height

      ! INTEGER,INTENT(in)         :: opt ! 0 for momentum, 1 for temperature, 2 for humidity
      INTEGER, INTENT(in)         :: RoughLenHeatMethod, StabilityMethod

      REAL(KIND(1d0)), INTENT(out):: avU10_ms, t2_C, q2_gkg, tskin_C, RH2
      REAL(KIND(1d0))::qa_gkg
      REAL(KIND(1d0)), PARAMETER::k = 0.4

      ! wind speed:
      CALL diagSfc( &
         0, &
         zMeas, avU1, 0d0, 10d0, avU10_ms, &
         VegFraction, &
         z0m, zdm, avdens, avcp, lv_J_kg, &
         avU1, Temp_C, qh, &
         RoughLenHeatMethod, StabilityMethod, tstep_real, dectime)

      ! temperature at 2 m agl:
      CALL diagSfc( &
         1, &
         zMeas, Temp_C, qh, 2d0, t2_C, &
         VegFraction, &
         z0m, zdm, avdens, avcp, lv_J_kg, &
         avU1, Temp_C, qh, &
         RoughLenHeatMethod, StabilityMethod, tstep_real, dectime)

      ! skin temperature:
      tskin_C = qh/(avdens*avcp)*RA + temp_C

      ! humidity:
      qa_gkg = RH2qa(avRH/100, Press_hPa, Temp_c)
      CALL diagSfc( &
         2, &
         zMeas, qa_gkg, qe, 2d0, q2_gkg, &
         VegFraction, &
         z0m, zdm, avdens, avcp, lv_J_kg, &
         avU1, Temp_C, qh, &
         RoughLenHeatMethod, StabilityMethod, tstep_real, dectime)
      ! re-examine if the diagnostic RH2 > 100% ?
      RH2 = qa2RH(q2_gkg, Press_hPa, Temp_c)
      IF (RH2 > 1) THEN
         ! if so, limit RH2 to 100%
         RH2 = 1d0
         ! and adjust the diagnostic q2_gkg
         q2_gkg = RH2qa(RH2, Press_hPa, Temp_c)
      END IF

   END SUBROUTINE SUEWS_cal_Diagnostics

   ! calculate several surface fraction related parameters
   SUBROUTINE SUEWS_cal_surf( &
      sfr, & !input
      vegfraction, ImpervFraction, PervFraction, NonWaterFraction) ! output
      IMPLICIT NONE

      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)::sfr
      REAL(KIND(1D0)), INTENT(OUT)::VegFraction
      REAL(KIND(1D0)), INTENT(OUT)::ImpervFraction
      REAL(KIND(1D0)), INTENT(OUT)::PervFraction
      REAL(KIND(1D0)), INTENT(OUT)::NonWaterFraction

      VegFraction = sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf)
      ImpervFraction = sfr(PavSurf) + sfr(BldgSurf)
      PervFraction = 1 - ImpervFraction
      NonWaterFraction = 1 - sfr(WaterSurf)

   END SUBROUTINE SUEWS_cal_surf

   SUBROUTINE diagSfc( &
      opt, &
      zMeas, xMeas, xFlux, zDiag, xDiag, &
      VegFraction, &
      z0m, zd, avdens, avcp, lv_J_kg, &
      avU1, Temp_C, qh, &
      RoughLenHeatMethod, StabilityMethod, tstep_real, dectime)
      ! TS 31 Jul 2018: removed dependence on surface variables (Tsurf, qsat)
      ! TS 26 Jul 2018: improved the calculation logic
      ! TS 05 Sep 2017: improved interface
      ! TS 20 May 2017: calculate surface-level diagonostics

      IMPLICIT NONE
      REAL(KIND(1d0)), INTENT(in) :: dectime
      REAL(KIND(1d0)), INTENT(in) :: qh ! sensible heat flux
      REAL(KIND(1d0)), INTENT(in) :: z0m, avdens, avcp, lv_J_kg, tstep_real
      REAL(KIND(1d0)), INTENT(in) :: avU1, Temp_C ! atmospheric level variables
      REAL(KIND(1d0)), INTENT(in) :: zDiag ! height for diagonostics
      REAL(KIND(1d0)), INTENT(in) :: zMeas! height for measurement
      REAL(KIND(1d0)), INTENT(in) :: zd ! displacement height
      REAL(KIND(1d0)), INTENT(in) :: xMeas ! measurement at height
      REAL(KIND(1d0)), INTENT(in) :: xFlux!
      REAL(KIND(1d0)), INTENT(in) :: VegFraction ! vegetation fraction

      INTEGER, INTENT(in)         :: opt ! 0 for momentum, 1 for temperature, 2 for humidity
      INTEGER, INTENT(in)         :: RoughLenHeatMethod, StabilityMethod

      REAL(KIND(1d0)), INTENT(out):: xDiag

      REAL(KIND(1d0)) :: L_mod
      REAL(KIND(1d0)) :: psimz0, psihzDiag, psihzMeas, psihz0, psimzDiag ! stability correction functions
      REAL(KIND(1d0)) :: z0h ! Roughness length for heat
      REAL(KIND(1d0)) :: zDiagzd! height for diagnositcs
      REAL(KIND(1d0)) :: zMeaszd
      REAL(KIND(1d0)) :: tlv, H_kms, TStar, zL, UStar
      REAL(KIND(1d0)), PARAMETER :: muu = 1.46e-5 !molecular viscosity
      REAL(KIND(1d0)), PARAMETER :: nan = -999
      REAL(KIND(1d0)), PARAMETER :: zdm = 0 ! assuming Displacement height is ZERO
      REAL(KIND(1d0)), PARAMETER::k = 0.4

      tlv = lv_J_kg/tstep_real !Latent heat of vapourisation per timestep
      zDiagzd = zDiag + z0m ! height at hgtX assuming Displacement height is ZERO; set lower limit as z0 to prevent arithmetic error, zd=0

      ! get !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
      CALL SUEWS_init_QH( &
         qh, avdens, avcp, qh, 0d0, dectime, & ! use qh as qh_obs to initialise H_init
         H_kms)!output

      ! redo the calculation for stability correction
      CALL STAB_lumps( &
         ! input
         StabilityMethod, &
         dectime, & !Decimal time
         zDiagzd, &     !Active measurement height (meas. height-displac. height)
         z0m, &     !Aerodynamic roughness length
         zdm, &     !Displacement height
         avU1, &    !Average wind speed
         Temp_C, &  !Air temperature
         H_kms, & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
         ! output:
         L_MOD, & !Obukhov length
         TStar, & !T*
         UStar, & !Friction velocity
         zL)!Stability scale

      !***************************************************************
      ! log-law based stability corrections:
      ! Roughness length for heat
      z0h = cal_z0V(RoughLenHeatMethod, z0m, VegFraction, UStar)

      ! stability correction functions
      ! momentum:
      psimzDiag = stab_psi_mom(StabilityMethod, zDiagzd/L_mod, zDiagzd/L_mod)
      ! psimz2=stab_fn_mom(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
      psimz0 = stab_psi_mom(StabilityMethod, z0m/L_mod, z0m/L_mod)

      ! heat and vapor: assuming both are the same
      ! psihz2=stab_fn_heat(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
      psihz0 = stab_psi_heat(StabilityMethod, z0h/L_mod, z0h/L_mod)

      !***************************************************************
      SELECT CASE (opt)
      CASE (0) ! wind (momentum) at hgtX=10 m
         zDiagzd = zDiag + z0m! set lower limit as z0h to prevent arithmetic error, zd=0

         ! stability correction functions
         ! momentum:
         psimzDiag = stab_psi_mom(StabilityMethod, zDiagzd/L_mod, zDiagzd/L_mod)
         psimz0 = stab_psi_mom(StabilityMethod, z0m/L_mod, z0m/L_mod)
         xDiag = UStar/k*(LOG(zDiagzd/z0m) - psimzDiag + psimz0) ! Brutsaert (2005), p51, eq.2.54

      CASE (1) ! temperature at hgtX=2 m
         zMeaszd = zMeas - zd
         zDiagzd = zDiag + z0h! set lower limit as z0h to prevent arithmetic error, zd=0

         ! heat and vapor: assuming both are the same
         psihzMeas = stab_psi_heat(StabilityMethod, zMeaszd/L_mod, zMeaszd/L_mod)
         psihzDiag = stab_psi_heat(StabilityMethod, zDiagzd/L_mod, zDiagzd/L_mod)
         ! psihz0=stab_fn_heat(StabilityMethod,z0h/L_mod,z0h/L_mod)
         xDiag = xMeas + xFlux/(k*UStar*avdens*avcp)*(LOG(zMeaszd/zDiagzd) - (psihzMeas - psihzDiag)) ! Brutsaert (2005), p51, eq.2.55
         !  IF ( ABS((LOG(z2zd/z0h)-psihz2+psihz0))>10 ) THEN
         !     PRINT*, '#####################################'
         !     PRINT*, 'xSurf',xSurf
         !     PRINT*, 'xFlux',xFlux
         !     PRINT*, 'k*us*avdens*avcp',k*us*avdens*avcp
         !     PRINT*, 'k',k
         !     PRINT*, 'us',us
         !     PRINT*, 'avdens',avdens
         !     PRINT*, 'avcp',avcp
         !     PRINT*, 'xFlux/X',xFlux/(k*us*avdens*avcp)
         !     PRINT*, 'stab',(LOG(z2zd/z0h)-psihz2+psihz0)
         !     PRINT*, 'LOG(z2zd/z0h)',LOG(z2zd/z0h)
         !     PRINT*, 'z2zd',z2zd,'L_mod',L_mod,'z0h',z0h
         !     PRINT*, 'z2zd/L_mod',z2zd/L_mod
         !     PRINT*, 'psihz2',psihz2
         !     PRINT*, 'psihz0',psihz0
         !     PRINT*, 'psihz2-psihz0',psihz2-psihz0
         !     PRINT*, 'xDiag',xDiag
         !     PRINT*, '*************************************'
         !  END IF

      CASE (2) ! humidity at hgtX=2 m
         zMeaszd = zMeas - zd
         zDiagzd = zDiag + z0h! set lower limit as z0h to prevent arithmetic error, zd=0

         ! heat and vapor: assuming both are the same
         psihzMeas = stab_psi_heat(StabilityMethod, zMeaszd/L_mod, zMeaszd/L_mod)
         psihzDiag = stab_psi_heat(StabilityMethod, zDiagzd/L_mod, zDiagzd/L_mod)
         ! psihz0=stab_fn_heat(StabilityMethod,z0h/L_mod,z0h/L_mod)

         xDiag = xMeas + xFlux/(k*UStar*avdens*tlv)*(LOG(zMeaszd/zDiagzd) - (psihzMeas - psihzDiag)) ! Brutsaert (2005), p51, eq.2.56

      END SELECT

   END SUBROUTINE diagSfc

   !===============set variable of invalid value to NAN=====================
   ELEMENTAL FUNCTION set_nan(x) RESULT(xx)
      IMPLICIT NONE
      REAL(KIND(1d0)), PARAMETER::pNAN = 30000 ! 30000 to prevent water_state being filtered out as it can be large
      REAL(KIND(1d0)), PARAMETER::NAN = -999
      REAL(KIND(1d0)), INTENT(in)::x
      REAL(KIND(1d0))::xx

      IF (ABS(x) > pNAN) THEN
         xx = NAN
      ELSE
         xx = x
      ENDIF

   END FUNCTION set_nan
   !========================================================================

   !===============the functions below are only for test in f2py conversion===
   FUNCTION square(x) RESULT(xx)
      IMPLICIT NONE
      REAL(KIND(1d0)), PARAMETER::pNAN = 9999
      REAL(KIND(1d0)), PARAMETER::NAN = -999
      REAL(KIND(1d0)), INTENT(in)::x
      REAL(KIND(1d0))::xx

      xx = x**2 + nan/pNAN
      xx = x**2

   END FUNCTION square

   FUNCTION square_real(x) RESULT(xx)
      IMPLICIT NONE
      REAL, PARAMETER::pNAN = 9999
      REAL, PARAMETER::NAN = -999
      REAL, INTENT(in)::x
      REAL::xx

      xx = x**2 + nan/pNAN
      xx = x**2

   END FUNCTION square_real

   SUBROUTINE output_name_n(i, name, group, aggreg)
      ! used by f2py module `SuPy` to handle output names
      IMPLICIT NONE
      ! the dimension is potentially incorrect,
      ! which should be consistent with that in output module
      INTEGER, INTENT(in) :: i
      CHARACTER(len=15), INTENT(out) :: name, group, aggreg

      INTEGER :: nVar
      nVar = SIZE(varListAll, dim=1)
      IF (i < nVar .AND. i > 0) THEN
         name = TRIM(varListAll(i)%header)
         group = TRIM(varListAll(i)%group)
         aggreg = TRIM(varListAll(i)%aggreg)
      ELSE
         name = ''
         group = ''
         aggreg = ''
      END IF

   END SUBROUTINE output_name_n

   SUBROUTINE output_size(nVar)
      ! used by f2py module `SuPy` to get size of the output list
      IMPLICIT NONE
      ! the dimension is potentially incorrect,
      ! which should be consistent with that in output module
      INTEGER, INTENT(out) :: nVar

      nVar = SIZE(varListAll, dim=1)

   END SUBROUTINE output_size

   SUBROUTINE SUEWS_cal_multitsteps( &
      MetForcingBlock, len_sim, &
      AerodynamicResistanceMethod, AH_MIN, AHProf_24hr, AH_SLOPE_Cooling, & ! input&inout in alphabetical order
      AH_SLOPE_Heating, &
      alb, AlbMax_DecTr, AlbMax_EveTr, AlbMax_Grass, &
      AlbMin_DecTr, AlbMin_EveTr, AlbMin_Grass, &
      alpha_bioCO2, alpha_enh_bioCO2, alt, BaseT, BaseTe, &
      BaseTHDD, beta_bioCO2, beta_enh_bioCO2, bldgH, CapMax_dec, CapMin_dec, &
      chAnOHM, CO2PointSource, cpAnOHM, CRWmax, CRWmin, DayWat, DayWatPer, &
      DecTreeH, Diagnose, DiagQN, DiagQS, DRAINRT, &
      dt_since_start, dqndt, qn1_av, dqnsdt, qn1_s_av, &
      EF_umolCO2perJ, emis, EmissionsMethod, EnEF_v_Jkm, endDLS, EveTreeH, FAIBldg, &
      FAIDecTree, FAIEveTree, Faut, FcEF_v_kgkm, FlowChange, &
      FrFossilFuel_Heat, FrFossilFuel_NonHeat, G1, G2, G3, G4, G5, G6, GDD_id, &
      GDDFull, Gridiv, gsModel, HDD_id, HumActivity_24hr, &
      IceFrac, Ie_a, Ie_end, Ie_m, Ie_start, &
      InternalWaterUse_h, IrrFracConif, IrrFracDecid, IrrFracGrass, EvapMethod, &
      kkAnOHM, Kmax, LAI_id, LAICalcYes, LAIMax, LAIMin, &
      LAIPower, LAIType, lat, lng, MaxConductance, MaxFCMetab, MaxQFMetab, &
      SnowWater, MinFCMetab, MinQFMetab, min_res_bioCO2, &
      NARP_EMIS_SNOW, NARP_TRANS_SITE, NetRadiationMethod, &
      NumCapita, OHM_coef, OHMIncQF, OHM_threshSW, &
      OHM_threshWD, PipeCapacity, PopDensDaytime, &
      PopDensNighttime, PopProf_24hr, PorMax_dec, PorMin_dec, &
      PrecipLimit, PrecipLimitAlb, &
      QF0_BEU, Qf_A, Qf_B, Qf_C, &
      RadMeltFact, RAINCOVER, RainMaxRes, resp_a, resp_b, &
      RoughLenHeatMethod, RoughLenMomMethod, RunoffToWater, S1, S2, &
      SatHydraulicConduct, SDDFull, sfr, SMDMethod, SnowAlb, SnowAlbMax, &
      SnowAlbMin, SnowPackLimit, SnowDens, SnowDensMax, SnowDensMin, SnowfallCum, snowFrac, &
      SnowLimBldg, SnowLimPaved, SnowPack, SnowProf_24hr, snowUse, SoilDepth, &
      soilstore_id, SoilStoreCap, StabilityMethod, startDLS, state_id, StateLimit, &
      StorageHeatMethod, StoreDrainPrm, SurfaceArea, Tair_av, tau_a, tau_f, tau_r, &
      T_CRITIC_Cooling, T_CRITIC_Heating, TempMeltFact, TH, &
      theta_bioCO2, timezone, TL, TrafficRate, TrafficUnits, &
      TraffProf_24hr, Ts5mindata_ir, tstep, tstep_prev, veg_type, &
      WaterDist, WaterUseMethod, WetThresh, &
      WUDay_id, DecidCap_id, albDecTr_id, albEveTr_id, albGrass_id, porosity_id, &
      WUProfA_24hr, WUProfM_24hr, Z, z0m_in, zdm_in, &
      dataOutBlockSUEWS, dataOutBlockSnow, dataOutBlockESTM, dataOutBlockRSL, &!output
      DailyStateBlock)

      IMPLICIT NONE
      ! input:
      ! met forcing block
      REAL(KIND(1D0)), DIMENSION(len_sim, 24), INTENT(IN) ::MetForcingBlock
      INTEGER, INTENT(IN) :: len_sim
      ! input variables
      INTEGER, INTENT(IN)::AerodynamicResistanceMethod
      INTEGER, INTENT(IN)::Diagnose
      INTEGER, INTENT(IN)::DiagQN
      INTEGER, INTENT(IN)::DiagQS
      INTEGER, INTENT(IN)::startDLS
      INTEGER, INTENT(IN)::endDLS
      INTEGER, INTENT(IN)::EmissionsMethod
      INTEGER, INTENT(IN)::Gridiv
      INTEGER, INTENT(IN)::gsModel
      INTEGER, INTENT(IN)::Ie_end
      INTEGER, INTENT(IN)::Ie_start
      INTEGER, INTENT(IN)::EvapMethod
      INTEGER, INTENT(IN)::LAICalcYes
      INTEGER, INTENT(IN)::NetRadiationMethod
      INTEGER, INTENT(IN)::OHMIncQF
      INTEGER, INTENT(IN)::RoughLenHeatMethod
      INTEGER, INTENT(IN)::RoughLenMomMethod
      INTEGER, INTENT(IN)::SMDMethod
      INTEGER, INTENT(IN)::snowUse
      INTEGER, INTENT(IN)::StabilityMethod
      INTEGER, INTENT(IN)::StorageHeatMethod
      INTEGER, INTENT(IN)::tstep
      INTEGER, INTENT(IN)::tstep_prev ! tstep size of the previous step
      INTEGER, INTENT(in)::dt_since_start ! time since simulation starts [s]
      INTEGER, INTENT(IN)::veg_type
      INTEGER, INTENT(IN)::WaterUseMethod

      INTEGER, DIMENSION(NVEGSURF), INTENT(IN)::LAIType

      REAL(KIND(1D0)), INTENT(IN)::AlbMax_DecTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMax_EveTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMax_Grass
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_DecTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_EveTr
      REAL(KIND(1D0)), INTENT(IN)::AlbMin_Grass
      REAL(KIND(1D0)), INTENT(IN)::alt
      ! REAL(KIND(1D0)),INTENT(IN)::avkdn
      ! REAL(KIND(1D0)),INTENT(IN)::avRh
      ! REAL(KIND(1D0)),INTENT(IN)::avU1
      REAL(KIND(1D0)), INTENT(IN)::BaseTHDD
      REAL(KIND(1D0)), INTENT(IN)::bldgH
      REAL(KIND(1D0)), INTENT(IN)::CapMax_dec
      REAL(KIND(1D0)), INTENT(IN)::CapMin_dec
      REAL(KIND(1D0)), INTENT(IN)::CO2PointSource
      REAL(KIND(1D0)), INTENT(IN)::CRWmax
      REAL(KIND(1D0)), INTENT(IN)::CRWmin
      REAL(KIND(1D0)), INTENT(IN)::DecTreeH
      REAL(KIND(1D0)), INTENT(IN)::DRAINRT
      REAL(KIND(1D0)), INTENT(IN)::EF_umolCO2perJ
      REAL(KIND(1D0)), INTENT(IN)::EnEF_v_Jkm
      REAL(KIND(1D0)), INTENT(IN)::EveTreeH
      REAL(KIND(1D0)), INTENT(IN)::FAIBldg
      REAL(KIND(1D0)), INTENT(IN)::FAIDecTree
      REAL(KIND(1D0)), INTENT(IN)::FAIEveTree
      REAL(KIND(1D0)), INTENT(IN)::Faut
      ! REAL(KIND(1D0)),INTENT(IN)::fcld_obs
      REAL(KIND(1D0)), INTENT(IN)::FlowChange
      REAL(KIND(1D0)), INTENT(IN)::FrFossilFuel_Heat
      REAL(KIND(1D0)), INTENT(IN)::FrFossilFuel_NonHeat
      REAL(KIND(1D0)), INTENT(IN)::G1
      REAL(KIND(1D0)), INTENT(IN)::G2
      REAL(KIND(1D0)), INTENT(IN)::G3
      REAL(KIND(1D0)), INTENT(IN)::G4
      REAL(KIND(1D0)), INTENT(IN)::G5
      REAL(KIND(1D0)), INTENT(IN)::G6
      REAL(KIND(1D0)), INTENT(IN)::InternalWaterUse_h
      REAL(KIND(1D0)), INTENT(IN)::IrrFracConif
      REAL(KIND(1D0)), INTENT(IN)::IrrFracDecid
      REAL(KIND(1D0)), INTENT(IN)::IrrFracGrass
      REAL(KIND(1D0)), INTENT(IN)::Kmax
      ! REAL(KIND(1D0)),INTENT(IN)::LAI_obs
      REAL(KIND(1D0)), INTENT(IN)::lat
      ! REAL(KIND(1D0)),INTENT(IN)::ldown_obs
      REAL(KIND(1D0)), INTENT(IN)::lng
      REAL(KIND(1D0)), INTENT(IN)::MaxFCMetab
      REAL(KIND(1D0)), INTENT(IN)::MaxQFMetab
      REAL(KIND(1D0)), INTENT(IN)::MinFCMetab
      REAL(KIND(1D0)), INTENT(IN)::MinQFMetab
      REAL(KIND(1D0)), INTENT(IN)::NARP_EMIS_SNOW
      REAL(KIND(1D0)), INTENT(IN)::NARP_TRANS_SITE
      REAL(KIND(1D0)), INTENT(IN)::PipeCapacity
      REAL(KIND(1D0)), INTENT(IN)::PopDensNighttime
      REAL(KIND(1D0)), INTENT(IN)::PorMax_dec
      REAL(KIND(1D0)), INTENT(IN)::PorMin_dec
      ! REAL(KIND(1D0)),INTENT(IN)::Precip
      REAL(KIND(1D0)), INTENT(IN)::PrecipLimit
      REAL(KIND(1D0)), INTENT(IN)::PrecipLimitAlb
      ! REAL(KIND(1D0)),INTENT(IN)::Press_hPa
      ! REAL(KIND(1D0)),INTENT(IN)::qh_obs
      ! REAL(KIND(1D0)),INTENT(IN)::qn1_obs
      ! REAL(KIND(1D0)),INTENT(IN)::qs_obs
      ! REAL(KIND(1D0)),INTENT(IN)::qf_obs
      REAL(KIND(1D0)), INTENT(IN)::RadMeltFact
      REAL(KIND(1D0)), INTENT(IN)::RAINCOVER
      REAL(KIND(1D0)), INTENT(IN)::RainMaxRes
      REAL(KIND(1D0)), INTENT(IN)::RunoffToWater
      REAL(KIND(1D0)), INTENT(IN)::S1
      REAL(KIND(1D0)), INTENT(IN)::S2
      REAL(KIND(1D0)), INTENT(IN)::SnowAlbMax
      REAL(KIND(1D0)), INTENT(IN)::SnowAlbMin
      REAL(KIND(1D0)), INTENT(IN)::SnowDensMax
      REAL(KIND(1D0)), INTENT(IN)::SnowDensMin
      REAL(KIND(1D0)), INTENT(IN)::SnowLimBldg
      REAL(KIND(1D0)), INTENT(IN)::SnowLimPaved
      ! REAL(KIND(1D0)),INTENT(IN)::snowFrac_obs
      REAL(KIND(1D0)), INTENT(IN)::SurfaceArea
      REAL(KIND(1D0)), INTENT(IN)::tau_a
      REAL(KIND(1D0)), INTENT(IN)::tau_f
      REAL(KIND(1D0)), INTENT(IN)::tau_r
      ! REAL(KIND(1D0)),INTENT(IN)::Temp_C
      REAL(KIND(1D0)), INTENT(IN)::TempMeltFact
      REAL(KIND(1D0)), INTENT(IN)::TH
      REAL(KIND(1D0)), INTENT(IN)::timezone
      REAL(KIND(1D0)), INTENT(IN)::TL
      REAL(KIND(1D0)), INTENT(IN)::TrafficUnits
      ! REAL(KIND(1D0)),INTENT(IN)::xsmd
      REAL(KIND(1D0)), INTENT(IN)::Z
      REAL(KIND(1D0)), INTENT(IN)::z0m_in
      REAL(KIND(1D0)), INTENT(IN)::zdm_in

      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_MIN
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_SLOPE_Cooling
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::AH_SLOPE_Heating
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::FcEF_v_kgkm
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::QF0_BEU
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_A
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_B
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Qf_C
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::Numcapita
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::PopDensDaytime
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::T_CRITIC_Cooling
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::T_CRITIC_Heating
      REAL(KIND(1D0)), DIMENSION(2), INTENT(IN)               ::TrafficRate
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::Ie_a
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::Ie_m
      REAL(KIND(1D0)), DIMENSION(3), INTENT(IN)               ::MaxConductance
      REAL(KIND(1D0)), DIMENSION(7), INTENT(IN)               ::DayWat
      REAL(KIND(1D0)), DIMENSION(7), INTENT(IN)               ::DayWatPer
      REAL(KIND(1D0)), DIMENSION(nsurf + 1), INTENT(IN)         ::OHM_threshSW
      REAL(KIND(1D0)), DIMENSION(nsurf + 1), INTENT(IN)         ::OHM_threshWD
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::chAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::cpAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::emis
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::kkAnOHM
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SatHydraulicConduct
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::sfr
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SnowPackLimit
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SoilDepth
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::SoilStoreCap
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::StateLimit
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(IN)           ::WetThresh
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::alpha_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::alpha_enh_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::BaseT
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::BaseTe
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::beta_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::beta_enh_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::GDDFull
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::LAIMax
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::LAIMin
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::min_res_bioCO2
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::resp_a
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::resp_b
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::SDDFull
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN)          ::SnowProf_24hr
      REAL(KIND(1D0)), DIMENSION(NVEGSURF), INTENT(IN)        ::theta_bioCO2
      REAL(KIND(1D0)), DIMENSION(4, NVEGSURF), INTENT(IN)      ::LAIPower
      REAL(KIND(1D0)), DIMENSION(nsurf + 1, 4, 3), INTENT(IN)     ::OHM_coef
      REAL(KIND(1D0)), DIMENSION(NSURF + 1, NSURF - 1), INTENT(IN) ::WaterDist
      REAL(KIND(1d0)), DIMENSION(:), INTENT(IN)               ::Ts5mindata_ir

      ! diurnal profile values for 24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::AHProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::HumActivity_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::PopProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::TraffProf_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::WUProfA_24hr
      REAL(KIND(1D0)), DIMENSION(0:23, 2), INTENT(IN) ::WUProfM_24hr
      ! ########################################################################################

      ! ########################################################################################
      ! inout variables
      ! OHM related:
      REAL(KIND(1d0)), INTENT(INOUT) ::qn1_av
      REAL(KIND(1d0)), INTENT(INOUT) ::dqndt
      REAL(KIND(1d0)), INTENT(INOUT) ::qn1_s_av
      REAL(KIND(1d0)), INTENT(INOUT) ::dqnsdt

      ! snow related:
      REAL(KIND(1D0)), INTENT(INOUT)                  ::SnowfallCum
      REAL(KIND(1D0)), INTENT(INOUT)                  ::SnowAlb
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::IceFrac
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowWater
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowDens
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::snowFrac
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT) ::SnowPack

      ! water balance related:
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::soilstore_id
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::state_id
      REAL(KIND(1D0)), DIMENSION(6, NSURF), INTENT(INOUT) ::StoreDrainPrm

      ! phenology related:
      REAL(KIND(1D0)), DIMENSION(NSURF), INTENT(INOUT)   ::alb
      REAL(KIND(1d0)), DIMENSION(5), INTENT(INOUT)       ::GDD_id !Growing Degree Days (see SUEWS_DailyState.f95)
      REAL(KIND(1d0)), DIMENSION(nvegsurf), INTENT(INOUT)::LAI_id !LAI for each veg surface [m2 m-2]
      REAL(KIND(1d0)), INTENT(INOUT)                    :: DecidCap_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albDecTr_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albEveTr_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: albGrass_id
      REAL(KIND(1d0)), INTENT(INOUT)                    :: porosity_id

      ! anthropogenic heat related:
      REAL(KIND(1d0)), DIMENSION(12), INTENT(INOUT) ::HDD_id !Heating Degree Days (see SUEWS_DailyState.f95)

      ! water use related:
      REAL(KIND(1d0)), DIMENSION(9), INTENT(INOUT)  ::WUDay_id

      ! ESTM related:
      REAL(KIND(1d0)), INTENT(INOUT)    ::Tair_av
      ! ########################################################################################

      ! ########################################################################################
      ! output variables
      ! REAL(KIND(1D0)),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT) ::datetimeBlock
      REAL(KIND(1D0)), DIMENSION(len_sim, ncolumnsDataOutSUEWS), INTENT(OUT) ::dataOutBlockSUEWS
      REAL(KIND(1D0)), DIMENSION(len_sim, ncolumnsDataOutSnow), INTENT(OUT) ::dataOutBlockSnow
      REAL(KIND(1d0)), DIMENSION(len_sim, ncolumnsDataOutESTM), INTENT(OUT) ::dataOutBlockESTM
      REAL(KIND(1d0)), DIMENSION(len_sim, ncolumnsDataOutRSL), INTENT(OUT) ::dataOutBlockRSL
      REAL(KIND(1d0)), DIMENSION(len_sim, ncolumnsDataOutDailyState), INTENT(OUT) ::DailyStateBlock
      ! ########################################################################################

      ! internal temporal iteration related variables
      INTEGER::dt_since_start_x ! time since simulation starts [s]

      ! model output blocks of the same size as met forcing block

      ! local variables
      ! length of met forcing block
      INTEGER :: ir
      ! met forcing variables
      INTEGER :: iy
      INTEGER :: id
      INTEGER :: it
      INTEGER :: imin
      INTEGER :: isec
      INTEGER, PARAMETER :: gridiv_x = 1 ! a dummy gridiv as this routine is only one grid
      REAL(KIND(1D0))::qn1_obs
      REAL(KIND(1D0))::qh_obs
      REAL(KIND(1D0))::qe_obs
      REAL(KIND(1D0))::qs_obs
      REAL(KIND(1D0))::qf_obs
      REAL(KIND(1D0))::avu1
      REAL(KIND(1D0))::avrh
      REAL(KIND(1D0))::Temp_C
      REAL(KIND(1D0))::Press_hPa
      REAL(KIND(1D0))::Precip
      REAL(KIND(1D0))::avkdn
      REAL(KIND(1D0))::snowFrac_obs
      REAL(KIND(1D0))::ldown_obs
      REAL(KIND(1D0))::fcld_obs
      REAL(KIND(1D0))::wu_m3
      REAL(KIND(1D0))::xsmd
      REAL(KIND(1D0))::LAI_obs
      REAL(KIND(1D0))::kdiff
      REAL(KIND(1D0))::kdir
      REAL(KIND(1D0))::wdir

      REAL(KIND(1D0)), DIMENSION(5)::datetimeLine
      REAL(KIND(1D0)), DIMENSION(ncolumnsDataOutSUEWS - 5)::dataOutLineSUEWS
      REAL(KIND(1D0)), DIMENSION(ncolumnsDataOutSnow - 5)::dataOutLineSnow
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutESTM - 5)::dataOutLineESTM
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutRSL - 5)::dataOutLineRSL
      REAL(KIND(1d0)), DIMENSION(ncolumnsDataOutDailyState - 5)::DailyStateLine

      REAL(KIND(1D0)), DIMENSION(len_sim, ncolumnsDataOutSUEWS, 1) ::dataOutBlockSUEWS_X
      REAL(KIND(1D0)), DIMENSION(len_sim, ncolumnsDataOutSnow, 1) ::dataOutBlockSnow_X
      REAL(KIND(1d0)), DIMENSION(len_sim, ncolumnsDataOutESTM, 1) ::dataOutBlockESTM_X
      REAL(KIND(1d0)), DIMENSION(len_sim, ncolumnsDataOutRSL, 1) ::dataOutBlockRSL_X
      ! REAL(KIND(1d0)),DIMENSION(len_sim,ncolumnsDataOutDailyState,1) ::DailyStateBlock_X

      ! REAL(KIND(1D0)),DIMENSION(:,:)          ::MetForcingData_grid

      ! get initial dt_since_start_x from dt_since_start, dt_since_start_x is used for Qn averaging. TS 28 Nov 2018
      dt_since_start_x = dt_since_start

      DO ir = 1, len_sim, 1
         ! =============================================================================
         ! === Translate met data from MetForcingBlock to variable names used in model ==
         ! =============================================================================
         iy = INT(MetForcingBlock(ir, 1)) !Integer variables
         id = INT(MetForcingBlock(ir, 2))
         it = INT(MetForcingBlock(ir, 3))
         imin = INT(MetForcingBlock(ir, 4))
         isec = 0 ! NOT used by SUEWS but by WRF-SUEWS via the cal_main interface
         qn1_obs = MetForcingBlock(ir, 5)      !Real values (kind(1d0))
         qh_obs = MetForcingBlock(ir, 6)
         qe_obs = MetForcingBlock(ir, 7)
         qs_obs = MetForcingBlock(ir, 8)
         qf_obs = MetForcingBlock(ir, 9)
         avu1 = MetForcingBlock(ir, 10)
         avrh = MetForcingBlock(ir, 11)
         Temp_C = MetForcingBlock(ir, 12)
         Press_hPa = MetForcingBlock(ir, 13)
         Precip = MetForcingBlock(ir, 14)
         avkdn = MetForcingBlock(ir, 15)
         snowFrac_obs = MetForcingBlock(ir, 16)
         ldown_obs = MetForcingBlock(ir, 17)
         fcld_obs = MetForcingBlock(ir, 18)
         wu_m3 = MetForcingBlock(ir, 19)
         xsmd = MetForcingBlock(ir, 20)
         LAI_obs = MetForcingBlock(ir, 21)
         kdiff = MetForcingBlock(ir, 22)
         kdir = MetForcingBlock(ir, 23)
         wdir = MetForcingBlock(ir, 24)

         CALL SUEWS_cal_Main( &
            AerodynamicResistanceMethod, AH_MIN, AHProf_24hr, AH_SLOPE_Cooling, & ! input&inout in alphabetical order
            AH_SLOPE_Heating, &
            alb, AlbMax_DecTr, AlbMax_EveTr, AlbMax_Grass, &
            AlbMin_DecTr, AlbMin_EveTr, AlbMin_Grass, &
            alpha_bioCO2, alpha_enh_bioCO2, alt, avkdn, avRh, avU1, BaseT, BaseTe, &
            BaseTHDD, beta_bioCO2, beta_enh_bioCO2, bldgH, CapMax_dec, CapMin_dec, &
            chAnOHM, CO2PointSource, cpAnOHM, CRWmax, CRWmin, DayWat, DayWatPer, &
            DecTreeH, Diagnose, DiagQN, DiagQS, DRAINRT, &
            dt_since_start_x, dqndt, qn1_av, dqnsdt, qn1_s_av, &
            EF_umolCO2perJ, emis, EmissionsMethod, EnEF_v_Jkm, endDLS, EveTreeH, FAIBldg, &
            FAIDecTree, FAIEveTree, Faut, FcEF_v_kgkm, fcld_obs, FlowChange, &
            FrFossilFuel_Heat, FrFossilFuel_NonHeat, G1, G2, G3, G4, G5, G6, GDD_id, &
            GDDFull, gridiv_x, gsModel, HDD_id, HumActivity_24hr, &
            IceFrac, id, Ie_a, Ie_end, Ie_m, Ie_start, imin, &
            InternalWaterUse_h, IrrFracConif, IrrFracDecid, IrrFracGrass, isec, it, EvapMethod, &
            iy, kkAnOHM, Kmax, LAI_id, LAICalcYes, LAIMax, LAIMin, LAI_obs, &
            LAIPower, LAIType, lat, ldown_obs, lng, MaxConductance, MaxFCMetab, MaxQFMetab, &
            SnowWater, MetForcingBlock, MinFCMetab, MinQFMetab, min_res_bioCO2, &
            NARP_EMIS_SNOW, NARP_TRANS_SITE, NetRadiationMethod, &
            NumCapita, OHM_coef, OHMIncQF, OHM_threshSW, &
            OHM_threshWD, PipeCapacity, PopDensDaytime, &
            PopDensNighttime, PopProf_24hr, PorMax_dec, PorMin_dec, &
            Precip, PrecipLimit, PrecipLimitAlb, Press_hPa, &
            QF0_BEU, Qf_A, Qf_B, Qf_C, &
            qn1_obs, qh_obs, qs_obs, qf_obs, &
            RadMeltFact, RAINCOVER, RainMaxRes, resp_a, resp_b, &
            RoughLenHeatMethod, RoughLenMomMethod, RunoffToWater, S1, S2, &
            SatHydraulicConduct, SDDFull, sfr, SMDMethod, SnowAlb, SnowAlbMax, &
            SnowAlbMin, SnowPackLimit, SnowDens, SnowDensMax, SnowDensMin, SnowfallCum, snowFrac, &
            SnowLimBldg, SnowLimPaved, snowFrac_obs, SnowPack, SnowProf_24hr, snowUse, SoilDepth, &
            soilstore_id, SoilStoreCap, StabilityMethod, startDLS, state_id, StateLimit, &
            StorageHeatMethod, StoreDrainPrm, SurfaceArea, Tair_av, tau_a, tau_f, tau_r, &
            T_CRITIC_Cooling, T_CRITIC_Heating, Temp_C, TempMeltFact, TH, &
            theta_bioCO2, timezone, TL, TrafficRate, TrafficUnits, &
            TraffProf_24hr, Ts5mindata_ir, tstep, tstep_prev, veg_type, &
            WaterDist, WaterUseMethod, WetThresh, wu_m3, &
            WUDay_id, DecidCap_id, albDecTr_id, albEveTr_id, albGrass_id, porosity_id, &
            WUProfA_24hr, WUProfM_24hr, xsmd, Z, z0m_in, zdm_in, &
            datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, dataOutLineRSL, & !output
            DailyStateLine)!output

         ! update dt_since_start_x for next iteration, dt_since_start_x is used for Qn averaging. TS 28 Nov 2018
         dt_since_start_x = dt_since_start_x + tstep

         !============ update and write out SUEWS_cal_DailyState ===============
         ! only works at the last timestep of a day
         ! CALL SUEWS_update_DailyState(&
         !      id, datetimeLine, &!input
         !      Gridiv, 1, &
         !      DailyStateLine, &
         !      DailyStateBlock_X)!inout
         DailyStateBlock(ir, :) = [datetimeLine, DailyStateLine]

         !============ write out results ===============
         ! works at each timestep
         CALL SUEWS_update_output( &
            SnowUse, storageheatmethod, &!input
            len_sim, 1, &
            ir, gridiv_x, datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, dataoutLineRSL, &!input
            dataOutBlockSUEWS_X, dataOutBlockSnow_X, dataOutBlockESTM_X, dataOutBlockRSL_X)!inout

      END DO

      dataOutBlockSUEWS = dataOutBlockSUEWS_X(:, :, 1)
      dataOutBlockSnow = dataOutBlockSnow_X(:, :, 1)
      dataOutBlockESTM = dataOutBlockESTM_X(:, :, 1)
      dataOutBlockRSL = dataOutBlockRSL_X(:, :, 1)
      ! DailyStateBlock=DailyStateBlock_X(:,:,1)

   END SUBROUTINE SUEWS_cal_multitsteps

   ! a wrapper of NARP_cal_SunPosition used by supy
   SUBROUTINE SUEWS_cal_sunposition( &
      year, idectime, UTC, locationlatitude, locationlongitude, locationaltitude, & !input
      sunazimuth, sunzenith) !output
      IMPLICIT NONE

      REAL(KIND(1D0)), INTENT(in) :: year, idectime, UTC, &
                                     locationlatitude, locationlongitude, locationaltitude
      REAL(KIND(1D0)), INTENT(out) ::sunazimuth, sunzenith

      CALL NARP_cal_SunPosition( &
         year, idectime, UTC, locationlatitude, locationlongitude, locationaltitude, &
         sunazimuth, sunzenith)

   END SUBROUTINE SUEWS_cal_sunposition

   ! function func(arg) result(retval)
   !    implicit none
   !    type :: arg
   !    type :: retval

   ! end function func

   function cal_tair_av(tair_av_prev, dt_since_start, tstep, temp_c) result(tair_av_next)
      ! calculate mean air temperature of past 24 hours
      ! TS, 17 Sep 2019
      implicit none
      real(KIND(1D0)), intent(in) :: tair_av_prev
      real(KIND(1D0)), intent(in) ::  temp_c
      integer, intent(in) ::  dt_since_start
      integer, intent(in) ::  tstep

      real(KIND(1D0)) ::tair_av_next

      real(KIND(1D0)), parameter:: len_day_s = 24*3600 ! day length in seconds
      real(KIND(1D0)):: len_cal_s ! length of average period in seconds
      real(KIND(1D0)):: temp_k ! temp in K

      ! determine the average period
      if (dt_since_start > len_day_s) then
         ! if simulation has been running over one day
         len_cal_s = len_day_s
      else
         ! if simulation has been running less than one day
         len_cal_s = dt_since_start + tstep
      end if
      temp_k = temp_c + 273.15
      tair_av_next = tair_av_prev*(len_cal_s - tstep*1.)/len_cal_s + temp_k*tstep/len_cal_s

   end function cal_tair_av

END MODULE SUEWS_Driver
