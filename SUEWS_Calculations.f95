!This subroutine does the actual calculations of the SUEWS code (mostly old SUEWS_Temporal).
!Made by LJ and HW Oct 2014
!Gives in the grid ID (Gridiv) and number of line in the met forcing data to be analyzed (ir)
!Last modification
!
!Last modification:
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

SUBROUTINE SUEWS_Calculations(Gridiv,ir,iMB,irMax)

  USE data_in
  USE time
  USE NARP_MODULE
  USE defaultNotUsed
  USE allocateArray
  USE sues_data
  USE snowMod
  USE gis_data
  USE initial
  USE moist
  USE mod_z
  USE mod_k
  USE solweig_module
  USE WhereWhen
  USE SUEWS_Driver
  USE VegPhenogy
  USE resist


  IMPLICIT NONE

  INTEGER:: Gridiv
  INTEGER::ir
  ! INTEGER::ih
  INTEGER::iMB
  ! LOGICAL        :: debug=.FALSE.
  ! REAL(KIND(1d0)):: idectime
  !real(kind(1d0)):: SnowDepletionCurve  !for SUEWS_Snow - not needed here (HCW 24 May 2016)
  ! REAL(KIND(1d0)):: LAI_wt
  ! REAL(KIND(1d0))::xBo



  INTEGER:: irMax

  !==================================================================
  !==================================================================


  !Translate all data to the variables used in the model calculations
  IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_Translate...'
  CALL SUEWS_Translate(Gridiv,ir,iMB)

  ! !==============main calculation start=======================
  ! IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_RoughnessParameters...'
  ! ! CALL SUEWS_cal_RoughnessParameters(Gridiv) ! Added by HCW 11 Nov 2014
  ! CALL SUEWS_cal_RoughnessParameters(&
  !      RoughLenMomMethod,sfr,areaZh,&!input
  !      bldgH,EveTreeH,DecTreeH,&
  !      porosity(id),FAIBldg,FAIEveTree,FAIDecTree,Z,&
  !      planF,&!output
  !      Zh,Z0m,Zdm,ZZD)
  !
  ! ! Calculate sun position
  ! IF(Diagnose==1) WRITE(*,*) 'Calling NARP_cal_SunPosition...'
  ! CALL NARP_cal_SunPosition(&
  !      year,&!input:
  !      dectime-halftimestep,&! sun position at middle of timestep before
  !      timezone,lat,lng,alt,&
  !      azimuth,zenith_deg)!output:
  ! !write(*,*) DateTime, timezone,lat,lng,alt,azimuth,zenith_deg
  !
  !
  ! !Call the SUEWS_cal_DailyState routine to get surface characteristics ready
  ! IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_DailyState...'
  ! ! CALL SUEWS_cal_DailyState(Gridiv)
  ! CALL SUEWS_cal_DailyState(&
  !      iy,id,it,imin,Gridiv,tstep,&!input
  !      WaterUseMethod,snowUse,Ie_start,Ie_end,&
  !      ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids,LAICalcYes,LAIType,&
  !      nsh_real,avkdn,Temp_C,Precip,BaseTHDD,&
  !      lat,Faut,LAI_obs,tau_a,tau_f,tau_r,&
  !      SnowDensMax,SnowDensMin,SnowAlbMin,&
  !      alBMax_DecTr,alBMax_EveTr,alBMax_Grass,&
  !      AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
  !      CapMax_dec,CapMin_dec,PorMax_dec,PorMin_dec,&
  !      Ie_a,Ie_m,DayWatPer,DayWat,SnowPack,&
  !      BaseT,BaseTe,GDDFull,SDDFull,LAIMin,LAIMax,LAIPower,dataOut,&
  !      a1,a2,a3,tstepcount,SnowAlb,DecidCap,albDecTr,albEveTr,albGrass,&!inout
  !      porosity,GDD,HDD,SnowDens,LAI,&
  !      DayofWeek,WU_Day,xBo)!output
  !
  !
  ! !Calculation of density and other water related parameters
  ! IF(Diagnose==1) WRITE(*,*) 'Calling LUMPS_cal_AtmMoist...'
  ! CALL LUMPS_cal_AtmMoist(&
  !      Temp_C,Press_hPa,avRh,dectime,&! input:
  !      lv_J_kg,lvS_J_kg,&! output:
  !      es_hPa,Ea_hPa,VPd_hpa,VPD_Pa,dq,dens_dry,avcp,avdens)
  !
  !
  ! !======== Calculate soil moisture =========
  ! CALL SUEWS_update_SoilMoist(&
  !      nsurf,ConifSurf,DecidSurf,GrassSurf,&!input
  !      NonWaterFraction,&
  !      soilstoreCap,sfr,soilmoist,&
  !      SoilMoistCap,SoilState,&!output
  !      vsmd,smd)
  !
  !
  ! ! ===================NET ALLWAVE RADIATION================================
  ! CALL SUEWS_cal_Qn(&
  !      NetRadiationMethod,snowUse,ldown_option,id,&!input
  !      Diagnose,snow_obs,ldown_obs,fcld_obs,&
  !      dectime,ZENITH_deg,avKdn,Temp_C,avRH,Press_hPa,qn1_obs,&
  !      SnowAlb,AlbedoChoice,DiagQN,&
  !      NARP_G,NARP_TRANS_SITE,NARP_EMIS_SNOW,IceFrac,sfr,emis,&
  !      alb,albDecTr,DecidCap,albEveTr,albGrass,surf,&!inout
  !      snowFrac,ldown,fcld,&!output
  !      qn1,qn1_SF,qn1_S,kclear,kup,lup,tsurf)
  !
  !
  ! ! ===================ANTHROPOGENIC HEAT FLUX================================
  ! ! ih=it-DLS
  ! ! IF(ih<0) ih=23
  ! !
  ! ! IF(AnthropHeatMethod==1) THEN
  ! !    IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_1...'
  ! !    CALL SAHP_1(qf_sahp,QF_SAHP_base,QF_SAHP_heat,id,ih,imin)
  ! !    qn1_bup=qn1
  ! !    qn1=qn1+QF_SAHP
  ! ! ELSEIF(AnthropHeatMethod==2) THEN
  ! !    IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_2...'
  ! !    CALL SAHP_2(qf_sahp,QF_SAHP_base,QF_SAHP_heat,id,ih,imin)
  ! !    qn1_bup=qn1
  ! !    qn1=qn1+QF_SAHP
  ! ! ELSEIF(AnthropHeatMethod==3) THEN
  ! !    IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_3...'
  ! !    CALL SAHP_3(qf_sahp,id,ih,imin)
  ! !    qn1_bup=qn1
  ! !    qn1=qn1+QF_SAHP
  ! ! ELSE
  ! !    qn1_bup=qn1
  ! !    qn1=qn1+qf
  ! ! ENDIF
  ! !
  ! ! ! -- qn1 is now QSTAR+QF (net all-wave radiation + anthropogenic heat flux)
  ! ! ! -- qn1_bup is QSTAR only
  ! ! IF(AnthropHeatMethod>=1) THEN
  ! !    qf=QF_SAHP
  ! ! ENDIF
  ! !
  ! ! ! Calculate CO2 fluxes from anthropogenic components
  ! ! IF(Diagnose==1) WRITE(*,*) 'Calling CO2_anthro...'
  ! ! CALL CO2_anthro(id,ih,imin)
  ! ! ! For the purpose of turbulent fluxes, remove QF from the net all-wave radiation
  ! ! qn1=qn1_bup  !Remove QF from QSTAR
  ! ! ! -- qn1 is now QSTAR only
  ! qf=0 ! AH disabled for the moment, TS 28 Sep 2017
  ! ! ! ===================ANTHROPOGENIC HEAT FLUX END================================
  !
  !
  ! ! =================STORAGE HEAT FLUX=======================================
  ! CALL SUEWS_cal_Qs(&
  !      StorageHeatMethod,OHMIncQF,Gridiv,id,Diagnose,sfr,&!input
  !      OHM_coef,OHM_threshSW,OHM_threshWD,&
  !      soilmoist,soilstoreCap,state,nsh,SnowUse,DiagQS,&
  !      HDD,MetForcingData,qf,qn1_bup,&
  !      alb,emis,cpAnOHM,kkAnOHM,chAnOHM,&
  !      AnthropHeatMethod,&
  !      qn1_store,qn1_S_store,qn1_av_store,qn1_S_av_store,&!inout
  !      surf,qn1_S,snowFrac,qs,&
  !      deltaQi,a1,a2,a3)!output
  !
  !
  !
  ! !==================Energy related to snow melting/freezing processes=======
  ! IF(Diagnose==1) WRITE(*,*) 'Calling MeltHeat'
  ! CALL Snow_cal_MeltHeat(&
  !      snowUse,&!input
  !      bldgsurf,nsurf,PavSurf,WaterSurf,&
  !      lvS_J_kg,lv_J_kg,tstep_real,&
  !      RadMeltFact,TempMeltFact,SnowAlbMax,SnowDensMin,&
  !      Temp_C,Precip,PrecipLimit,PrecipLimitAlb,&
  !      nsh_real,waterdens,&
  !      sfr,Tsurf_ind,state,qn1_ind_snow,Meltwaterstore,deltaQi,&
  !      SnowPack,snowFrac,SnowAlb,SnowDens,& ! inout
  !      mwh,fwh,Qm,QmFreez,QmRain,CumSnowfall,snowCalcSwitch,&!output
  !      Qm_melt,Qm_freezState,Qm_rain,FreezMelt,FreezState,FreezStateVol,&
  !      rainOnSnow,SnowDepth,mw_ind,veg_fr)
  !
  !
  ! !==========================Turbulent Fluxes================================
  ! IF(Diagnose==1) WRITE(*,*) 'Calling LUMPS_cal_QHQE...'
  ! !Calculate QH and QE from LUMPS
  ! CALL LUMPS_cal_QHQE(&
  !      veg_type,& !input
  !      snowUse,&
  !      qn1_bup,&
  !      qf,&
  !      qs,&
  !      Qm,&
  !      Temp_C,&
  !      Veg_Fr,&
  !      avcp,&
  !      Press_hPa,&
  !      lv_J_kg,&
  !      tstep_real,&
  !      DRAINRT,&
  !      nsh_real,&
  !      Precip,&
  !      RainMaxRes,&
  !      RAINCOVER,&
  !      sfr(ivConif+2:ivGrass+2),&
  !      LAI(id-1,Gridiv),&
  !      LAImax,&
  !      LAImin,&
  !      H_mod,& !output
  !      E_mod,&
  !      psyc_hPa,&
  !      s_hPa,&
  !      sIce_hpa,&
  !      TempVeg,&
  !      VegPhenLumps)
  !
  !
  ! IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_WaterUse...'
  ! !Gives the external and internal water uses per timestep
  ! CALL SUEWS_cal_WaterUse(&
  !      nsh_real,&! input:
  !      SurfaceArea,&
  !      sfr,&
  !      IrrFracConif,&
  !      IrrFracDecid,&
  !      IrrFracGrass,&
  !      DayofWeek(id,:),&
  !      WUProfA_tstep,&
  !      WUProfM_tstep,&
  !      InternalWaterUse_h,&
  !      HDD(id-1,:),&
  !      WU_Day(id-1,:),&
  !      WaterUseMethod,&
  !      NSH,it,imin,DLS,&
  !      OverUse,&
  !      WUAreaEveTr_m2,&!  output:
  !      WUAreaDecTr_m2,&
  !      WUAreaGrass_m2,&
  !      WUAreaTotal_m2,&
  !      wu_EveTr,&
  !      wu_DecTr,&
  !      wu_Grass,&
  !      wu_m3,&
  !      int_wu,&
  !      ext_wu)
  !
  !
  ! !===============Resistance Calculations=======================
  ! CALL SUEWS_cal_Resistance(&
  !      StabilityMethod,&!input:
  !      Diagnose,AerodynamicResistanceMethod,RoughLenHeatMethod,snowUse,&
  !      id,it,gsModel,SMDMethod,&
  !      qh_obs,avdens,avcp,h_mod,qn1_bup,dectime,zzd,z0M,zdm,&
  !      avU1,Temp_C,L_mod,UStar,VegFraction,avkdn,&
  !      Kmax,&
  !      g1,g2,g3,g4,&
  !      g5,g6,s1,s2,&
  !      th,tl,&
  !     !  SurfaceChar(Gridiv,c_GsKmax),&
  !     !  SurfaceChar(Gridiv,c_GsG1),SurfaceChar(Gridiv,c_GsG2),&
  !     !  SurfaceChar(Gridiv,c_GsG3),SurfaceChar(Gridiv,c_GsG4),&
  !     !  SurfaceChar(Gridiv,c_GsG5),SurfaceChar(Gridiv,c_GsG6),&
  !     !  SurfaceChar(Gridiv,c_GsS1),SurfaceChar(Gridiv,c_GsS2),&
  !     !  SurfaceChar(Gridiv,c_GsTH),SurfaceChar(Gridiv,c_GsTL),&
  !      dq,xsmd,vsmd,MaxConductance,LAIMax,LAI(id-1,:),snowFrac,sfr,&
  !      Tstar,&!output
  !      psim,gsc,ResistSurf,RA,RAsnow,rb)
  !
  !
  ! !========= CO2-related calculations ================================
  ! ! Calculate CO2 fluxes from biogenic components
  ! ! IF(Diagnose==1) WRITE(*,*) 'Calling CO2_biogen...'
  ! ! CALL CO2_biogen
  ! ! ! Sum anthropogenic and biogenic CO2 flux components to find overall CO2 flux
  ! ! Fc = Fc_anthro + Fc_biogen
  ! Fc=0 ! disbled for the moment. TS 28 Sep 2017
  ! !========= CO2-related calculations end================================
  !
  !
  ! !============= calculate water balance =============
  ! CALL SUEWS_cal_Water(&
  !      Diagnose,&!input
  !      snowUse,NonWaterFraction,addPipes,addImpervious,addVeg,addWaterBody,&
  !      state,soilmoist,sfr,surf,WaterDist,nsh_real,&
  !      drain_per_tstep,&  !output
  !      drain,AddWaterRunoff,&
  !      AdditionalWater,runoffPipes,runoff_per_interval,&
  !      addWater,stateOld,soilmoistOld)
  ! !============= calculate water balance end =============
  !
  ! !======== Evaporation and surface state ========
  ! CALL SUEWS_cal_QE(&
  !      Diagnose,&!input
  !      id,tstep,imin,it,ity,snowfractionchoice,snowCalcSwitch,DayofWeek,CRWmin,CRWmax,&
  !      nsh_real,lvS_J_kg,lv_j_kg,avdens,waterdens,avRh,Press_hPa,Temp_C,&
  !      RAsnow,psyc_hPa,avcp,sIce_hPa,&
  !      PervFraction,vegfraction,addimpervious,qn1_SF,qf,qs,vpd_hPa,s_hPa,&
  !      ResistSurf,ra,rb,tstep_real,snowdensmin,precip,PipeCapacity,RunoffToWater,&
  !      NonWaterFraction,wu_EveTr,wu_DecTr,wu_Grass,addVeg,addWaterBody,SnowLimPaved,SnowLimBuild,&
  !      SurfaceArea,drain,WetThresh,stateOld,mw_ind,soilstorecap,rainonsnow,&
  !      freezmelt,freezstate,freezstatevol,Qm_Melt,Qm_rain,Tsurf_ind,sfr,StateLimit,surf,&
  !      runoff_per_interval,& ! inout:
  !      state,soilmoist,SnowPack,snowFrac,MeltWaterStore,&
  !      SnowDepth,iceFrac,addwater,addwaterrunoff,SnowDens,SurplusEvap,&
  !      snowProf,& ! output:
  !      runoffSnow,runoff,runoffSoil,chang,changSnow,SnowToSurf,snowD,ev_snow,SnowRemoval,&
  !      evap,rss_nsurf,p_mm,rss,qe,state_per_tstep,NWstate_per_tstep,qeOut,&
  !      swe,ev,chSnow_per_interval,ev_per_tstep,qe_per_tstep,runoff_per_tstep,&
  !      surf_chang_per_tstep,runoffPipes,mwstore,runoffwaterbody,FlowChange,&
  !      runoffAGimpervious_m3,runoffAGveg_m3,runoffWaterBody_m3,runoffPipes_m3)
  !
  ! !======== Evaporation and surface state end========
  !
  ! !============ Sensible heat flux ===============
  ! CALL SUEWS_cal_QH(&
  !      1,&
  !      qn1_bup,&
  !      qf,&
  !      QmRain,&
  !      qeOut,&
  !      qs,&
  !      QmFreez,&
  !      qm,&
  !      avdens,&
  !      avcp,&
  !      tsurf,&
  !      Temp_C,&
  !      ra,&
  !      qh)
  ! !============ Sensible heat flux end===============
  !
  ! !=== Horizontal movement between soil stores ===
  ! ! Now water is allowed to move horizontally between the soil stores
  ! IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_cal_HorizontalSoilWater...'
  ! ! CALL SUEWS_cal_HorizontalSoilWater
  ! CALL SUEWS_cal_HorizontalSoilWater(&
  !      sfr,&! input: ! surface fractions
  !      SoilStoreCap,&!Capacity of soil store for each surface [mm]
  !      SoilDepth,&!Depth of sub-surface soil store for each surface [mm]
  !      SatHydraulicConduct,&!Saturated hydraulic conductivity for each soil subsurface [mm s-1]
  !      SurfaceArea,&!Surface area of the study area [m2]
  !      NonWaterFraction,&! sum of surface cover fractions for all except water surfaces
  !      tstep_real,& !tstep cast as a real for use in calculations
  !      SoilMoist,&! inout:!Soil moisture of each surface type [mm]
  !      runoffSoil,&!Soil runoff from each soil sub-surface [mm]
  !      runoffSoil_per_tstep&!  output:!Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)
  !      )
  !
  ! !========== Calculate soil moisture ============
  ! CALL SUEWS_cal_SoilMoist(&
  !      nsurf,SMDMethod,&
  !      xsmd,NonWaterFraction,SoilMoistCap,SoilStoreCap,surf_chang_per_tstep,&
  !      soilmoist,soilmoistOld,sfr,smd,smd_nsurf,tot_chang_per_tstep,&
  !      SoilState)
  !
  !
  ! !============ surface-level diagonostics ===============
  ! CALL SUEWS_cal_Diagnostics(&
  !      tsurf,qh,&
  !      Press_hPa,qeOut,&
  !      UStar,veg_fr,z0m,L_mod,k,avdens,avcp,lv_J_kg,tstep_real,&
  !      RoughLenHeatMethod,StabilityMethod,&
  !      avU10_ms,t2_C,q2_gkg)!output
  ! !============ surface-level diagonostics end ===============
  !
  ! !============ update and write out SUEWS_cal_DailyState ===============
  ! ! only works at the last timestep of a day
  ! CALL SUEWS_update_DailyState(&
  !      iy,id,it,imin,&!input
  !      GDD,HDD,LAI,&
  !      DecidCap,albDecTr,albEveTr,albGrass,porosity,&
  !      WU_Day,&
  !      nsh_real,deltaLAI,VegPhenLumps,&
  !      SnowAlb,SnowDens,&
  !      xBo,a1,a2,a3,&
  !      Gridiv,GridIDmatrix,&!input
  !      FileCode,FileOutputPath,&
  !      DailyStateFirstOpen)
  !
  !
  ! !==============main calculation end=======================
  CALL SUEWS_cal_Main(&
       a1,a2,a3,addImpervious,AdditionalWater,addPipes,addVeg,addWater,addWaterBody,&
       AddWaterRunoff,AerodynamicResistanceMethod,AH_MIN,AHProf_tstep,&
       AH_SLOPE_Cooling,AH_SLOPE_Heating,alb,albDecTr,AlbedoChoice,albEveTr,albGrass,&
       alBMax_DecTr,alBMax_EveTr,alBMax_Grass,AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
       alpha_bioCO2,alpha_enh_bioCO2,alt,areaZh,avdens,avkdn,avRh,&
       avU1,avU10_ms,azimuth,BaseT,BaseTe,BaseTHDD,beta_bioCO2,beta_enh_bioCO2,&
       BiogenCO2Code,bldgH,CapMax_dec,CapMin_dec,chang,changSnow,chAnOHM,&
       chSnow_per_interval,cpAnOHM,CRWmax,CRWmin,CumSnowfall,&
       dataOut,dataOutESTM,DayofWeek,DayWat,DayWatPer,DecidCap,dectime,DecTreeH,dens_dry,&
       Diagnose,DiagQN,DiagQS,DLS,drain_per_tstep,DRAINRT,Ea_hPa,EF_umolCO2perJ,emis,&
       EmissionsMethod,E_mod,EnEF_v_Jkm,es_hPa,ev,evap,EveTreeH,ev_per_tstep,ev_snow,&
       ext_wu,FAIBldg,FAIDecTree,FAIEveTree,Faut,Fc,Fc_anthro,Fc_biogen,Fc_build,&
       FcEF_v_kgkm,fcld,fcld_obs,Fc_metab,Fc_photo,Fc_respi,Fc_traff,&
       FlowChange,FreezMelt,FrFossilFuel_Heat,FrFossilFuel_NonHeat,fwh,&
       G1,G2,G3,G4,G5,G6,GDD,GDDFull,Gridiv,gsc,gsModel,halftimestep,HDD,&
       H_mod,HumActivity_tstep,IceFrac,id,Ie_a,Ie_end,Ie_m,Ie_start,imin,&
       InternalWaterUse_h,int_wu,IrrFracConif,IrrFracDecid,IrrFracGrass,it,ity,iy,k,&
       kclear,kkAnOHM,Kmax,kup,kup_ind_snow,LAI,LAICalcYes,LAIMax,LAIMin,LAI_obs,LAIPower,LAIType,&
       lat,ldown,ldown_obs,ldown_option,L_mod,lng,lup,MaxConductance,MaxQFMetab,&
       Meltwaterstore,MetForcingData,MinQFMetab,min_res_bioCO2,mwh,mw_ind,mwstore,&
       NARP_EMIS_SNOW,NARP_G,NARP_TRANS_SITE,ncolumnsDataOut,NetRadiationMethod,&
       NonWaterFraction,nsh,nsh_real,NumberOfGrids,NumCapita,&
       NWstate_per_tstep,OHM_coef,OHMIncQF,OHM_threshSW,OHM_threshWD,OverUse,&
       PervFraction,PipeCapacity,planF,p_mm,PopDensDaytime,PopDensNighttime,&
       PopProf_tstep,PorMax_dec,PorMin_dec,porosity,Precip,PrecipLimit,PrecipLimitAlb,&
       Press_hPa,psim,q2_gkg,qeOut,qe_per_tstep,qf,QF0_BEU,Qf_A,Qf_B,Qf_C,QF_SAHP,qh,&
       qh_obs,qh_r,Qm,QmFreez,Qm_freezState,Qm_melt,Qm_rain,QmRain,qn1,qn1_av_store,&
       qn1_ind_snow,qn1_obs,qn1_S,qn1_S_av_store,qn1_SF,qn1_S_store,qn1_store,qs,RA,&
       RadMeltFact,RAINCOVER,RainMaxRes,rainOnSnow,ReadLinesMetdata,ResistSurf,resp_a,&
       resp_b,RoughLenHeatMethod,RoughLenMomMethod,rss,rss_nsurf,runoff,&
       runoffAGimpervious,runoffAGveg,runoff_per_interval,runoff_per_tstep,&
       runoffPipes,runoffPipes_m3,runoffSnow,runoffSoil,runoffSoil_per_tstep,&
       RunoffToWater,runoffwaterbody,runoffWaterBody_m3,S1,S2,SatHydraulicConduct,&
       SDDFull,sfr,smd,SMDMethod,smd_nsurf,SnowAlb,SnowAlbMax,SnowAlbMin,&
       snowD,SnowDens,SnowDensMax,SnowDensMin,SnowDepth,snowFrac,&
       snowfractionchoice,SnowLimBuild,SnowLimPaved,snow_obs,SnowPack,snowProf,&
       SnowRemoval,SnowToSurf,snowUse,SoilDepth,soilmoist,SoilState,soilstoreCap,&
       StabilityMethod,state,StateLimit,state_per_tstep,StorageHeatMethod,surf,&
       SurfaceArea,surf_chang_per_tstep,SurplusEvap,swe,t2_C,Tair24HR,tau_a,tau_f,tau_r,&
       T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,TempMeltFact,TempVeg,TH,theta_bioCO2,&
       timezone,TL,tot_chang_per_tstep,TrafficRate,TrafficUnits,TraffProf_tstep,Ts5mindata_ir,&
       Tstar,tstep,tstepcount,tstep_real,tsurf,tsurf_ind,Tsurf_ind_snow,UStar,VegFraction,veg_type,VPD_Pa,&
       waterdens,WaterDist,WaterUseMethod,WetThresh,WUAreaDecTr_m2,WUAreaEveTr_m2,&
       WUAreaGrass_m2,WUAreaTotal_m2,WU_Day,wu_DecTr,wu_EveTr,wu_Grass,wu_m3,&
       WUProfA_tstep,WUProfM_tstep,xBo,xsmd,year,Z,Z0m,Zdm,zenith_deg,Zh)



  !============ update and write out SUEWS_cal_DailyState ===============
  ! only works at the last timestep of a day
  CALL SUEWS_update_DailyState(&
       iy,id,it,imin,&!input
       GDD,HDD,LAI,&
       DecidCap,albDecTr,albEveTr,albGrass,porosity,&
       WU_Day,&
       nsh_real,deltaLAI,VegPhenLumps,&
       SnowAlb,SnowDens,&
       xBo,a1,a2,a3,&
       Gridiv,GridIDmatrix,&
       FileCode,FileOutputPath,&
       DailyStateFirstOpen)!inout

  !============ write out results ===============
  ! works at each timestep
  CALL SUEWS_update_output(&
       AdditionalWater,alb,avkdn,avU10_ms,azimuth,&!input
       chSnow_per_interval,dectime,&
       drain_per_tstep,E_mod,ev_per_tstep,ext_wu,Fc,Fc_build,fcld,&
       Fc_metab,Fc_photo,Fc_respi,Fc_traff,FlowChange,freezMelt,&
       Gridiv,h_mod,id,id_prev_t,imin,int_wu,ir,it,iy,iy_prev_t,&
       kup,kup_ind_snow,LAI,ldown,l_mod,lup,MeltWaterStore,mwh,&
       mw_ind,MwStore,ncolumnsDataOut,ncolumnsDataOutSnow,&
       nsh_real,NumberOfGrids,NWstate_per_tstep,Precip,q2_gkg,&
       qeOut,qf,qh,QH_r,Qm,QmFreez,Qm_freezState,Qm_melt,Qm_rain,&
       QmRain,qn1,qn1_ind_snow,qn1_S,qn1_SF,qs,RA,rainOnSnow,&
       ReadLinesMetdata,resistsurf,runoffAGimpervious,runoffAGveg,&
       runoff_per_tstep,runoffPipes,runoffSoil_per_tstep,&
       runoffWaterBody,sfr,smd,smd_nsurf,SnowAlb,SnowDens,&
       snowDepth,SnowFrac,SnowPack,SnowRemoval,SNOWuse,SoilState,&
       state,state_per_tstep,surf_chang_per_tstep,swe,t2_C,&
       tot_chang_per_tstep,tsurf,Tsurf_ind_snow,UStar,wu_DecTr,&
       wu_EveTr,wu_Grass,z0m,zdm,zenith_deg,&
       dataOut,dataOutSnow)!inout



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

  !write(*,*) DecidCap(id), id, it, imin, 'Calc - before translate back'
  !write(*,*) iy, id, it, imin, 'Calc - before translate back'
  !if(Gridiv==1)  write(*,*) iy, id, it, imin, HDD(id-1,5), HDD(id,5), HDD(id-1,6), HDD(id,6)
  !if(id==12) pause
  !write(*,*) ' '

  IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_TranslateBack...'
  CALL SUEWS_TranslateBack(Gridiv,ir,irMax)


!!!if((id <=3 .or. id > 364).and. it == 0 .and. imin == 0) pause
!!!if((id <=3 .or. id > 364).and. it == 23 .and. imin == 55) pause
!!!if((id >=100 .or. id < 103).and. it == 0 .and. imin == 0) pause

  ! write(*,*) '------------'

END SUBROUTINE SUEWS_Calculations
