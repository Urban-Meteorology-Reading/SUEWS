!This subroutine does the actual calculations of the SUEWS code (mostly old SUEWS_Temporal).
!Made by LJ and HW Oct 2014
!Gives in the grid ID (Gridiv) and number of line in the met forcing data to be analyzed (ir)
!Last modification
!
!Last modification:
! TS  02 May 2018 - added explict interfaces
! LJ  15 Jun 2017 - Add parts where NonWaterFraction=0 is used as a divider so that the code won't crash with 100% water
! TS  20 May 2017 - Add surface-level diagonostics: T2, Q2, U10
! HCW 09 Dec 2016 - Add zenith and azimuth to output file
! HCw 24 Aug 2016 - smd and state for each surface set to NAN in output file if surface does not exist.
!                 - Added Fc to output file
! HCW 21 Jul 2016 - Set soil variables to -999 in output when grid is 100% water surface.
! HCW 29 Jun 2016 - Commented out StateDay and SoilMoistDay as creates jumps and should not be needed.
!                   Would not work unless each met block consists of a whole day for each grid.
! TS 09 Mar 2016  - Added AnOHM subroutine to calculate heat storage
! HCW 10 Mar 2016 - Calculation of soil moisture deficit of vegetated surfaces added (vsmd)
! LJ 2 Jan 2016   - Calculation of snow fraction moved from SUEWS_Calculations to SUEWS_Snow
! HCW 12 Nov 2015 - Added z0m and zdm to output file
! HCW 05 Nov 2015 - Changed Kawai et al. (2007) z0v calculation so VegFraction(veg+soil) rather than veg_fr(veg+soil+water) is used.
! HCW 29 Jun 2015 - Added albEveTr and albGrass
! HCW 25 Jun 2015 - Fixed bug in LAI calculation at year change.
!                 - Changed AlbDec to use id value (not (id-1) value) to avoid issues at year change.
! HCW 27 Apr 2015 - Correction to tot_chang_per_tstep calculation (water balance should now close)
! HCW 16 Feb 2015 - Updated water balance calculations
!                 - Corrected area-averaged calculations (soil moisture, drain, two versions of state with/out water)
!                 - Replaced soilmoist_state variable with SoilState (as seems to be duplicate)
! HCW 15 Jan 2015 - Added switch OHMIncQF to calculate QS with (1, default) or without (0) QF added to QSTAR
! To do
!      - add iy and imin to output files (may impact LUMPS_RunoffFromGrid)
!      - phase out _per_interval water balance variables
!      - renormalise by NonWaterFraction where necessary
!      - Update Snow subroutines similarly in terms of water balance
!==================================================================

SUBROUTINE SUEWS_Calculations(Gridiv, ir, iMB, irMax)
  USE data_in, ONLY: diagnose, ah_min, ah_slope_cooling, ah_slope_heating, &
       alt, avkdn, avrh, avu1, basethdd, diagqn, diagqs, drainrt, ef_umolco2perj, &
       emissionsmethod, enef_v_jkm, enddls, fcef_v_kgkm, fcld_obs, &
       frfossilfuel_heat, frfossilfuel_nonheat, ity, &
       laicalcyes, lai_obs, lat, ldown_obs, lng, maxqfmetab, minqfmetab, &
       netradiationmethod, numcapita, ohmincqf, popdensdaytime, popdensnighttime, &
       precip, press_hpa, qf0_beu, qf_a, qf_b, qf_c, qh_obs, qn1_obs, &
       raincover, rainmaxres, &
       roughlenmommethod, smdmethod, snow_obs, snowuse, startdls, &
       storageheatmethod, t_critic_cooling, t_critic_heating, temp_c, &
       timezone, trafficrate, trafficunits, waterusemethod, xsmd
  USE time, ONLY : iy, iy_prev_t, id, id_prev_t, it, imin, dectime
  USE allocateArray, ONLY: &
       ahprof_tstep, alb, albdectr, albevetr, albgrass, &
       AlbMax_DecTr, AlbMax_EveTr, AlbMax_grass, &
       AlbMin_dectr, AlbMin_evetr, AlbMin_grass, &
       alpha_bioco2, alpha_enh_bioco2, baset, basete, &
       beta_bioco2, beta_enh_bioco2, capmax_dec, capmin_dec, &
       chanohm, cpanohm, decidcap, emis, gdd, gddfull, hdd, &
       humactivity_tstep, icefrac, kkanohm, &
       lai, laimax, laimin, laipower, laitype, maxconductance, &
       meltwaterstore, metforcingdata_grid, min_res_bioco2, &
       narp_emis_snow, narp_trans_site, &
       ohm_coef, ohm_threshsw, ohm_threshwd, popprof_tstep, &
       pormax_dec, pormin_dec, porosity, &
       qn1_av_store_grid, qn1_s_av_store_grid, qn1_s_store_grid, qn1_store_grid, &
       resp_a, resp_b, sathydraulicconduct, sddfull, &
       sfr, snowd, snowdens, snowfrac, snowpack, &
       soildepth, soilmoist, soilstorecap, state, statelimit, &
       surf, tair24hr, theta_bioco2, traffprof_tstep, ts5mindata_ir, &
       waterdist, wetthresh, wu_day, wuprofa_tstep, wuprofm_tstep, &
       datetimeline, dataoutlinesuews, dataoutlinesnow, &
       dataoutlineestm, dailystateline, dataoutdailystate, &
       dataoutsuews, dataoutsnow, dataoutestm
  USE sues_data, ONLY: &
       aerodynamicresistancemethod, daywat, daywatper, faut, flowchange, &
       ie_a, ie_end, ie_m, ie_start, internalwateruse_h, &
       irrfracconif, irrfracdecid, irrfracgrass,&
       pipecapacity, roughlenheatmethod, runofftowater, stabilitymethod,&
       surfacearea, tstep
  USE snowMod, ONLY: &
       crwmax, crwmin, preciplimit, preciplimitalb, radmeltfact,&
       snowalb, snowAlbMax, snowAlbMin, &
       snowdensmax, snowdensmin, snowfallcum, snowlimbuild, &
       snowlimpaved, snowprof,&
       tau_a, tau_f, tau_r, tempmeltfact
  USE gis_data, ONLY: &
       bldgh, dectreeh, evetreeh, faibldg, faidectree, faievetree, veg_type
  USE initial, ONLY : NumberOfGrids, ReadLinesMetdata
  USE mod_z, ONLY : z, z0m_in, zdm_in
  USE SUEWS_Driver, ONLY: suews_update_output
  USE resist, ONLY: g1, g2, g3, g4, g5, g6, gsmodel, kmax, s1, s2, th, tl
  USE DailyState_module, ONLY : SUEWS_update_DailyState
  USE SUEWS_Driver,ONLY: SUEWS_cal_Main

  IMPLICIT NONE

  INTEGER :: Gridiv
  INTEGER :: ir
  INTEGER :: iMB
  INTEGER :: irMax

  !==================================================================


  !Translate all data to the variables used in the model calculations
  IF(Diagnose==1) WRITE(*, *) 'Calling SUEWS_Translate...'
  CALL SUEWS_Translate(Gridiv, ir, iMB)


  IF(Diagnose==1) print*, 'Calling SUEWS_cal_Main...'
  IF(Diagnose==1) print*, 'z0m_in =',z0m_in

  CALL SUEWS_cal_Main(&
       AerodynamicResistanceMethod, AH_MIN, AHProf_tstep, AH_SLOPE_Cooling, & ! input&inout in alphabetical order
       AH_SLOPE_Heating, alb, albDecTr, albEveTr, albGrass, AlbMax_DecTr, &
       AlbMax_EveTr, AlbMax_Grass, AlbMin_DecTr, AlbMin_EveTr, AlbMin_Grass, &
       alpha_bioCO2, alpha_enh_bioCO2, alt, avkdn, avRh, avU1, BaseT, BaseTe, &
       BaseTHDD, beta_bioCO2, beta_enh_bioCO2, bldgH, CapMax_dec, CapMin_dec, &
       chAnOHM, cpAnOHM, CRWmax, CRWmin, DayWat, DayWatPer, &
       DecidCap, dectime, DecTreeH, Diagnose, DiagQN, DiagQS, DRAINRT, &
       EF_umolCO2perJ, emis, EmissionsMethod, EnEF_v_Jkm, endDLS, EveTreeH, FAIBldg, &
       FAIDecTree, FAIEveTree, Faut, FcEF_v_kgkm, fcld_obs, FlowChange, &
       FrFossilFuel_Heat, FrFossilFuel_NonHeat, G1, G2, G3, G4, G5, G6, GDD, &
       GDDFull, Gridiv, gsModel, HDD, HumActivity_tstep, &
       IceFrac, id, id_prev_t, Ie_a, Ie_end, Ie_m, Ie_start, imin, &
       InternalWaterUse_h, IrrFracConif, IrrFracDecid, IrrFracGrass, it, ity, &
       iy, iy_prev_t, kkAnOHM, Kmax, LAI, LAICalcYes, LAIMax, LAIMin, LAI_obs, &
       LAIPower, LAIType, lat, ldown_obs, lng, MaxConductance, MaxQFMetab, &
       MeltWaterStore, MetForcingData_grid, MinQFMetab, min_res_bioCO2, &
       NARP_EMIS_SNOW, NARP_TRANS_SITE, NetRadiationMethod, &
       NumCapita, OHM_coef, OHMIncQF, OHM_threshSW, &
       OHM_threshWD, PipeCapacity, PopDensDaytime, &
       PopDensNighttime, PopProf_tstep, PorMax_dec, PorMin_dec, porosity, &
       Precip, PrecipLimit, PrecipLimitAlb, Press_hPa, QF0_BEU, Qf_A, Qf_B, &
       Qf_C, qh_obs, qn1_av_store_grid, qn1_obs, qn1_S_av_store_grid, qn1_S_store_grid, &
       qn1_store_grid, RadMeltFact, RAINCOVER, RainMaxRes, resp_a, resp_b, &
       RoughLenHeatMethod, RoughLenMomMethod, RunoffToWater, S1, S2, &
       SatHydraulicConduct, SDDFull, sfr, SMDMethod, SnowAlb, SnowAlbMax, &
       SnowAlbMin, snowD, SnowDens, SnowDensMax, SnowDensMin, SnowfallCum, snowFrac, &
       SnowLimBuild, SnowLimPaved, snow_obs, SnowPack, SnowProf, snowUse, SoilDepth, &
       soilmoist, soilstoreCap, StabilityMethod, startDLS, state, StateLimit, &
       StorageHeatMethod, surf, SurfaceArea, Tair24HR, tau_a, tau_f, tau_r, &
       T_CRITIC_Cooling, T_CRITIC_Heating, Temp_C, TempMeltFact, TH, &
       theta_bioCO2, timezone, TL, TrafficRate, TrafficUnits, &
       TraffProf_tstep, Ts5mindata_ir, tstep, veg_type, &
       WaterDist, WaterUseMethod, WetThresh, WU_Day, WUProfA_tstep, &
       WUProfM_tstep, xsmd, z, z0m_in, zdm_in, &
       datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, &!output
       DailyStateLine)!output


  !============ update and write out SUEWS_cal_DailyState ===============
  ! only works at the last timestep of a day
  CALL SUEWS_update_DailyState(&
       iy, id, it, imin, dectime, &!input
       Gridiv, NumberOfGrids, &
       DailyStateLine, &
       dataOutDailyState)!inout

  !============ write out results ===============
  ! works at each timestep
  CALL SUEWS_update_output(&
       SnowUse, storageheatmethod, &!input
       ReadLinesMetdata, NumberOfGrids, &
       ir, gridiv, datetimeLine, dataOutLineSUEWS, dataOutLineSnow, dataOutLineESTM, &!input
       dataOutSUEWS, dataOutSnow, dataOutESTM)!inout


  ! NB: CBL disabled for the moment for interface improvement
  ! can be decoupled from SUEWS? TS 28 Sep 2017
  ! IF(CBLuse>=1)THEN ! If CBL is used, calculated Temp_C and RH are replaced with the obs.
  !    IF(Diagnose==1) WRITE(*,*) 'Calling CBL...'
  !    CALL CBL(ir,iMB,Gridiv)   !ir=1 indicates first row of each met data block
  ! ENDIF


  ! NB: SOLWEIG can be treated as a separate part: can be decoupled from SUEWS? TS 28 Sep 2017
  ! ===================SOLWEIG OUTPUT ========================================
  ! IF (SOLWEIGuse==1) THEN
  !    IF (OutInterval==imin) THEN
  !       IF (RunForGrid==-999) THEN
  !          IF(Diagnose==1) WRITE(*,*) 'Calling SOLWEIG_2014a_core...'
  !          CALL SOLWEIG_2014a_core(iMB)
  !          SolweigCount=SolweigCount+1
  !       ELSE
  !          IF (Gridiv == RunForGrid) THEN
  !             IF(Diagnose==1) WRITE(*,*) 'Calling SOLWEIG_2014a_core...'
  !             CALL SOLWEIG_2014a_core(iMB)
  !             SolweigCount=SolweigCount+1
  !          ENDIF
  !       ENDIF
  !    ENDIF
  ! ELSE
  ! SOLWEIGpoi_out=0 ! NB: turn off SOLWEIG for the moment
  ! ENDIF
  ! ===================SOLWEIG END================================

  IF(Diagnose==1) WRITE(*, *) 'Calling SUEWS_TranslateBack...'
  CALL SUEWS_TranslateBack(Gridiv, ir, irMax)

END SUBROUTINE SUEWS_Calculations
