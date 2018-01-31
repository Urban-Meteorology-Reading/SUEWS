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
  USE DailyState_module,ONLY:SUEWS_update_DailyState


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


  CALL SUEWS_cal_Main(&
       AerodynamicResistanceMethod,AH_MIN,AHProf_tstep,AH_SLOPE_Cooling,& ! input&inout in alphabetical order
       AH_SLOPE_Heating,alb,albDecTr,albEveTr,albGrass,alBMax_DecTr,&
       alBMax_EveTr,alBMax_Grass,AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
       alpha_bioCO2,alpha_enh_bioCO2,alt,avkdn,avRh,avU1,BaseT,BaseTe,&
       BaseTHDD,beta_bioCO2,beta_enh_bioCO2,bldgH,CapMax_dec,CapMin_dec,&
       chAnOHM,cpAnOHM,CRWmax,CRWmin,CumSnowfall,DayofWeek,DayWat,DayWatPer,&
       DecidCap,dectime,DecTreeH,Diagnose,DiagQN,DiagQS,DLS,DRAINRT,&
       EF_umolCO2perJ,emis,EmissionsMethod,EnEF_v_Jkm,EveTreeH,FAIBldg,&
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
       Qf_C,qh_obs,qn1_av_store,qn1_obs,qn1_S_av_store,qn1_S_store,&
       qn1_store,RadMeltFact,RAINCOVER,RainMaxRes,resp_a,resp_b,&
       RoughLenHeatMethod,RoughLenMomMethod,RunoffToWater,S1,S2,&
       SatHydraulicConduct,SDDFull,sfr,SMDMethod,SnowAlb,SnowAlbMax,&
       SnowAlbMin,snowD,SnowDens,SnowDensMax,SnowDensMin,snowFrac,&
       SnowLimBuild,SnowLimPaved,snow_obs,SnowPack,SnowProf,snowUse,SoilDepth,&
       soilmoist,soilstoreCap,StabilityMethod,state,StateLimit,&
       StorageHeatMethod,surf,SurfaceArea,Tair24HR,tau_a,tau_f,tau_r,&
       T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,TempMeltFact,TH,&
       theta_bioCO2,timezone,TL,TrafficRate,TrafficUnits,&
       TraffProf_tstep,Ts5mindata_ir,tstep,veg_type,&
       WaterDist,WaterUseMethod,WetThresh,WU_Day,WUProfA_tstep,&
       WUProfM_tstep,xsmd,year,Z,&
       datetimeLine,dataOutLineSUEWS,dataOutLineSnow,dataOutLineESTM,&!output
       DailyStateLine)!output


  !============ update and write out SUEWS_cal_DailyState ===============
  ! only works at the last timestep of a day
  CALL SUEWS_update_DailyState(&
       iy,id,it,imin,dectime,&!input
       Gridiv,NumberOfGrids,nsh_real,&
       DailyStateLine,&
       dataOutDailyState)!inout

  !============ write out results ===============
  ! works at each timestep
  CALL SUEWS_update_output(&
       SnowUse,storageheatmethod,&!input
       ReadLinesMetdata,NumberOfGrids,&
       ir,gridiv,datetimeLine,dataOutLineSUEWS,dataOutLineSnow,dataOutLineESTM,&!input
       dataOutSUEWS,dataOutSnow,dataOutESTM)!inout


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
