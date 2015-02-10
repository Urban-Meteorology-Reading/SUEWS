!Translates new input arrays to actual model variables
! Runs at the start of SUEWS_Calculations to transfer correct info
!for each grid to the actual model variables
!Made by HW&LJ Oct 2014
!Last modified: HCW 28 Nov 2014
! To Do:
!	- Adjust model to calculate LAI per surface
!	- Adjust model for SM per surface (measured characteristics)
!       - Code BareSoil properly
!       - Code separate heights and FAI for evergreen and deciduous trees
!===================================================================================
subroutine SUEWS_Translate(Gridiv,ir,iMB)

  use Initial
  use allocateArray   ! defines: SiteInfo, sfr, (PavSurf, BldgSurf, etc) 
  use ColNamesInputFiles
  use ColNamesModelDailyState
  use data_in	      ! defines: lat, lng, PopDaytime, PopNighttime, DayLightSavingDay, QF variables
  use defaultnotUsed
  use gis_data        ! defines: areaZh, VegFraction, veg_fr, veg_type, BldgH, TreeH, FAIBldg, FAITree, Alt
  use mod_z	      ! defines: z0m, zdm	
  use ohm_calc        ! defines: OHM_coef
  use resist	      ! defines: G1-G6, TH, TL, S1, S2, Kmax
  use run_info	      ! for file_qs (still needed??)	 
  use snowMod	      ! defines: alb_snow, etc		
  use sues_data       ! defines: SurfaceArea, IrrFractionTrees, IrrFracConif, IrrFracDecid, IrrFracGrass, Irrigation variables
  use time
  
  IMPLICIT NONE

  integer::Gridiv,&   ! Index of the analysed grid (Gridcounter)
           ir,&       ! Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
           iMB,&      !	Chunk of met data
           id_prev

  integer::iv, j
  integer::i, ii, iii
  
  real (Kind(1d0)):: skip
  real (Kind(1d0)):: FCskip = -9  !NULL value used for output to FileChoices	


  !write(*,*) '---- SUEWS_Translate ----' 
  !write(*,*) 'Year:', SurfaceChar(Gridiv,c_Year)
  !write(*,*) 'Grid:', SurfaceChar(Gridiv,c_Grid)
  !write(*,*) 'Gridiv:', Gridiv
  !write(*,*) 'Met block (iMB or iv):',iMB
  !write(*,*) 'Met line (ir):',ir
  !write(*,*) '----' 
  
  ! =================================================================================
  ! ======= Translate inputs from SurfaceChar to variable names used in model =======
  ! =================================================================================
  ! ---- Latitude and longitude
  lat = SurfaceChar(Gridiv,c_lat)
  lng = SurfaceChar(Gridiv,c_lng)
  ! ---- Altitude [m] 
  Alt = SurfaceChar(Gridiv,c_Alt)
  ! ---- Surface area [ha]
  SurfaceArea_ha = SurfaceChar(Gridiv,c_Area)
  ! Change from ha to m2 (was in RunControlByGridByYear)
  SurfaceArea = SurfaceArea_ha*10000 !Change surface area from ha to m^2
 
  ! ---- Surface fractions (previously in LUMPS_gis_read)
  sfr(PavSurf)   = SurfaceChar(Gridiv,c_FrPaved)  ! Paved
  sfr(BldgSurf)  = SurfaceChar(Gridiv,c_FrBuilt)  ! Built
  sfr(ConifSurf) = SurfaceChar(Gridiv,c_FrEveTr)  ! Everg
  sfr(DecidSurf) = SurfaceChar(Gridiv,c_FrDecTr)  ! Decid
  sfr(GrassSurf) = SurfaceChar(Gridiv,c_FrGrass)  ! Grass
  sfr(BSoilSurf) = SurfaceChar(Gridiv,c_FrBSoil)  ! BSoil
  sfr(WaterSurf) = SurfaceChar(Gridiv,c_FrWater)  ! Water

  !write(*,*) '----------'
  !write(*,*) Gridiv
  !write(*,*) sum(sfr)
  !write(*,*) sfr
  !write(*,*) '----------'

  ! Check the surface fractions add up to 1 (or close to 1)
  if(sum(sfr)>1.001.or.sum(sfr)<0.999) call ErrorHint(10,'SiteSelect.txt - check surface fractions',sum(sfr),notUsed,notUsedI)
  
  ! ---- Irrigated fractions
  IrrFracConif = SurfaceChar(Gridiv,c_IrrEveTrFrac)  ! Everg
  IrrFracDecid = SurfaceChar(Gridiv,c_IrrDecTrFrac)  ! Decid
  IrrFracGrass = SurfaceChar(Gridiv,c_IrrGrassFrac)  ! Grass
  
  IrrFractionTrees = IrrFracConif+IrrFracDecid	!! change this after testing (delete IrrFractionTrees)
  

  ! ---------------------------------------------------------------------------------  
  ! --------- Surface cover calculations (previously in LUMPS_gis_read) -------------
   
  ! ---- Buildings and trees fraction ----
  areaZh = (sfr(BldgSurf) + sfr(ConifSurf) + sfr(DecidSurf))
  
  ! ---- Vegetated fraction ----
  VegFraction = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf))
  !VegFraction = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))
  
  ! ---- Vegetated fraction (for LUMPS) ----
  ! For LUMPS, vegetated fraction includes Water and Bare soil surfaces
  IF(veg_type==1) THEN          ! area vegetated
  	veg_fr = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf) + sfr(WaterSurf))
  ELSEIF(veg_type==2) THEN      ! area irrigated
        veg_fr = (IrrFracConif*sfr(ConifSurf) + IrrFracDecid*sfr(DecidSurf) + IrrFracGrass*sfr(GrassSurf))
  ENDIF
  
  ! ---------------------------------------------------------------------------------

  ! ---- Heights & frontal areas
  BldgH = SurfaceChar(Gridiv,c_HBuilt)	  ! Building height [m]
  TreeH = (SurfaceChar(Gridiv,c_HEveTr) + SurfaceChar(Gridiv,c_HDecTr))/2	! Tree height [m]
  !! Code separate heights for evergreen and deciduous trees later !!
  FAIBldg = SurfaceChar(Gridiv,c_FAIBuilt)  ! Frontal area index for buildings
  FAITree = (SurfaceChar(Gridiv,c_FAIEveTr ) + SurfaceChar(Gridiv,c_FAIDecTr ))/2 ! Frontal area index for trees
  !! Code separate heights for evergreen and deciduous trees later !!
  z0m = SurfaceChar(Gridiv,c_z0m) 	  ! Roughness length [m]
  zdm = SurfaceChar(Gridiv,c_zdm) 	  ! Displacement height [m]
  
  ! ---- Population
  PopDensDaytime   = SurfaceChar(Gridiv,c_PopDensDay)	! Daytime population density
  PopDensNighttime = SurfaceChar(Gridiv,c_PopDensNight)  ! Night-time population density
  
  NumCapita = PopDensDaytime      ! Pop density [ha-1] !! Use Daytime pop density for NumCapita for testing! 
  
  ! ---- Albedo
  alb(1:nsurf) = SurfaceChar(Gridiv,c_Alb)
  alb_snow     = SurfaceChar(Gridiv,c_SnowAlb)
  
  ! ---- Emissivity
  emis(1:nsurf) = SurfaceChar(Gridiv,c_Emis)
  emis_snow = SurfaceChar(Gridiv,c_SnowEmis)
  
  ! ---- Storage capacities
  Surf(1,1:nsurf) = SurfaceChar(Gridiv,c_StorMin)   ! Minimum	
  Surf(5,1:nsurf) = SurfaceChar(Gridiv,c_StorMax)   ! Maximum	
  
  ! ---- Drainage
  Surf(2,1:nsurf) = SurfaceChar(Gridiv,c_DrEq)      ! Drainage equation
  Surf(3,1:nsurf) = SurfaceChar(Gridiv,c_DrCoef1)   ! Drainage coef 1
  Surf(4,1:nsurf) = SurfaceChar(Gridiv,c_DrCoef2)   ! Drainage coef 2
  
  ! ---- Soil store capacity
  SoilStoreCap(1:nsurf) = SurfaceChar(Gridiv,c_SoilStCap)      
  
  ! ---- Limit of SWE (each surface except Water)
  snowD(1:(nsurf-1)) = SurfaceChar(Gridiv,c_SnowLimPat(1:(nsurf-1)))
  
  ! ---- Snow limit for removal (only impervious surfaces)
  SnowLimPaved = SurfaceChar(Gridiv,c_SnowLimRem(PavSurf)) 
  SnowLimBuild = SurfaceChar(Gridiv,c_SnowLimRem(BldgSurf))
        
  ! ---- Soil characteristics (each surface except Water)
  VolSoilMoistCap    (1:(nsurf-1)) = SurfaceChar(Gridiv,c_VolSMCap(1:(nsurf-1))) ! Volumetric soil moisture capacity [m3 m-3]
  SatHydraulicConduct(1:(nsurf-1)) = SurfaceChar(Gridiv,c_KSat(1:(nsurf-1))) ! Hydraulic conductivity of saturated soil [mm s-1] 
  !SoilDensity          (1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilDens(1:(nsurf-1))) ! Soil density [units??]
  ! Not yet implemented in model
  !InfiltrationRate  (1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilInfRate(1:(nsurf-1))) ! Infiltration rate [mm h-1]
  
  !! Observed soil characteristics
  !SoilDensity  (1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilDens(1:(nsurf-1))) ! Soil density [units??]
  !SoilDepthMeas(1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSMDepth(1:(nsurf-1)))
  !SmCap        (1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSMMax(1:(nsurf-1)))
  !SoilRocks    (1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSNRFrac(1:(nsurf-1)))
  !!Obs soil characteristics now in SUEWS_Soil, i.e. per surface; single value was given previously in FunctionalTypes
  !!Take first row here for testing !! Need to alter model later...
  SoilDensity    = SurfaceChar(Gridiv,c_SoilDens(1))
  SoilDepthMeas = SurfaceChar(Gridiv,c_ObsSMDepth(1))
  SmCap         = SurfaceChar(Gridiv,c_ObsSMMax(1))
  SoilRocks     = SurfaceChar(Gridiv,c_ObsSNRFrac(1))
  
  ! ---- Vegetation characteristics (pervious surfaces)
  BaseT  (1:nvegsurf) = SurfaceChar(Gridiv,c_BaseT)   	    
  BaseTe (1:nvegsurf) = SurfaceChar(Gridiv,c_BaseTe)   
  GDDFull(1:nvegsurf) = SurfaceChar(Gridiv,c_GDDFull)  
  SDDFull(1:nvegsurf) = SurfaceChar(Gridiv,c_SDDFull)  
  LAIMin (1:nvegsurf) = SurfaceChar(Gridiv,c_LAIMin)  
  LAIMax (1:nvegsurf) = SurfaceChar(Gridiv,c_LAIMax)  
  MaxConductance(1:nvegsurf) = SurfaceChar(Gridiv,c_GsMax)  
  
  ! ---- LAI characteristics  !! Model currently uses one single Eq & set of powers; adjust later for each veg type
  !write(*,*) SurfaceChar(Gridiv,c_LAIEq)  
  !write(*,*) SurfaceChar(Gridiv,c_LeafGP1)  
  !write(*,*) SurfaceChar(Gridiv,c_LeafGP2)  
  !write(*,*) SurfaceChar(Gridiv,c_LeafOP1)  
  !write(*,*) SurfaceChar(Gridiv,c_LeafOP2)  
  !! Set all rows the same for now in input file & use Decid row in model
  LAItype = SurfaceChar(Gridiv,c_LAIEq(ivDecid))
  LAIPower(1) = SurfaceChar(Gridiv,c_LeafGP1(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(2) = SurfaceChar(Gridiv,c_LeafGP2(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(3) = SurfaceChar(Gridiv,c_LeafOP1(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(4) = SurfaceChar(Gridiv,c_LeafOP2(ivDecid)) ! 4 powers (not 4 veg types!)  
  
  ! ---- LUMPS-related parameters
  DRAINRT    = SurfaceChar(Gridiv,c_LUMPSDr)  	 ! LUMPS Drainage rate [mm h-1]
  RAINCOVER  = SurfaceChar(Gridiv,c_LUMPSCover)   ! LUMPS Limit when surface totally wet [mm]
  RAINMAXRES = SurfaceChar(Gridiv,c_LUMPSMaxRes)  ! LUMPS Maximum water bucket reservoir
  
  ! ---- NARP-related parameters
  TRANS_SITE = SurfaceChar(Gridiv,c_NARPTrans)    ! NARP atmospheric transmissivity
   
  ! ---- Snow-related characteristics
  RadMeltFact    = SurfaceChar(Gridiv,c_SnowRMFactor)  
  TempMeltFact   = SurfaceChar(Gridiv,c_SnowTMFactor)  
  albSnowMin     = SurfaceChar(Gridiv,c_SnowAlbMin)  
  albSnowMax	 = SurfaceChar(Gridiv,c_SnowAlbMax) 
  tau_a		 = SurfaceChar(Gridiv,c_Snowtau_a) 
  tau_f		 = SurfaceChar(Gridiv,c_Snowtau_f) 
  PrecipLimitAlb = SurfaceChar(Gridiv,c_SnowPlimAlb) 
  DensSnowMin	 = SurfaceChar(Gridiv,c_SnowSDMin) 
  DensSnowMax	 = SurfaceChar(Gridiv,c_SnowSDMax)
  tau_r   	 = SurfaceChar(Gridiv,c_Snowtau_r)
  CRWMin 	 = SurfaceChar(Gridiv,c_SnowCRWMin)
  CRWMax	 = SurfaceChar(Gridiv,c_SnowCRWMax)
  PrecipLimit 	 = SurfaceChar(Gridiv,c_SnowPLimSnow)

  ! ---- Conductance parameters
  G1   = SurfaceChar(Gridiv,c_GsG1)  
  G2   = SurfaceChar(Gridiv,c_GsG2)  
  G3   = SurfaceChar(Gridiv,c_GsG3)  
  G4   = SurfaceChar(Gridiv,c_GsG4)  
  G5   = SurfaceChar(Gridiv,c_GsG5)  
  G6   = SurfaceChar(Gridiv,c_GsG6)
  TH   = SurfaceChar(Gridiv,c_GsTH)
  TL   = SurfaceChar(Gridiv,c_GsTL)  
  S1   = SurfaceChar(Gridiv,c_GsS1)  
  S2   = SurfaceChar(Gridiv,c_GsS2)  
  Kmax = SurfaceChar(Gridiv,c_GsKmax)  
  
  
  ! ---- Pipe capacity (was from SiteSpecificParam.txt)
  PipeCapacity = SurfaceChar(Gridiv,c_PipeCapacity)
  
  ! ---- Water flows (was from SiteSpecificParam.txt)
  FlowChange = SurfaceChar(Gridiv,c_FlowChange)
  RunoffToWater = SurfaceChar(Gridiv,c_RunoffToWater)                  
  
  ! ---- Daylight saving (was from ModelledYears.txt)
  DayLightSavingDay(1) = SurfaceChar(Gridiv,c_StartDLS)
  DayLightSavingDay(2) = SurfaceChar(Gridiv,c_EndDLS)
                  
  ! ---- OHM coeffs (was in SUEWS_OHMnew.f95, subroutine OHMinitialize)                 
  OHM_coef = 0 ! Initialise OHM_coef 
  ! Surface types in OHM_coef: Paved, Roof, Conif, Decid, Grass, BareSoil, Water, CANYON, Snow
  ! No canyon in SurfaceChar, so 
  !  transfer coeffs for surface types 1-7,
  !  then skip row in OHM_coef (canyon),
  !  then transfer coeffs for snow surface (8th surface in SurfaceChar; 9th surface in OHM_Coefs)
  ! Summer wet
  OHM_coef(1:nsurf,1,1) = SurfaceChar(Gridiv,c_a1_SWet(1:nsurf)) !1:nsurf a1 Summer wet
  OHM_coef(nsurf+2,1,1) = SurfaceChar(Gridiv,c_a1_SWet(nsurf+1)) !Snow    a1 Summer wet 
  OHM_coef(1:nsurf,1,2) = SurfaceChar(Gridiv,c_a2_SWet(1:nsurf)) !1:nsurf a2 Summer wet
  OHM_coef(nsurf+2,1,2) = SurfaceChar(Gridiv,c_a2_SWet(nsurf+1)) !Snow    a2 Summer wet 
  OHM_coef(1:nsurf,1,3) = SurfaceChar(Gridiv,c_a3_SWet(1:nsurf)) !1:nsurf a3 Summer wet
  OHM_coef(nsurf+2,1,3) = SurfaceChar(Gridiv,c_a3_SWet(nsurf+1)) !Snow    a3 Summer wet 
  ! Summer dry
  OHM_coef(1:nsurf,2,1) = SurfaceChar(Gridiv,c_a1_SDry(1:nsurf)) !1:nsurf a1 Summer dry
  OHM_coef(nsurf+2,2,1) = SurfaceChar(Gridiv,c_a1_SDry(nsurf+1)) !Snow    a1 Summer dry 
  OHM_coef(1:nsurf,2,2) = SurfaceChar(Gridiv,c_a2_SDry(1:nsurf)) !1:nsurf a2 Summer dry
  OHM_coef(nsurf+2,2,2) = SurfaceChar(Gridiv,c_a2_SDry(nsurf+1)) !Snow    a2 Summer dry 
  OHM_coef(1:nsurf,2,3) = SurfaceChar(Gridiv,c_a3_SDry(1:nsurf)) !1:nsurf a3 Summer dry
  OHM_coef(nsurf+2,2,3) = SurfaceChar(Gridiv,c_a3_SDry(nsurf+1)) !Snow    a3 Summer dry 
  ! Winter wet
  OHM_coef(1:nsurf,3,1) = SurfaceChar(Gridiv,c_a1_WWet(1:nsurf)) !1:nsurf a1 Winter wet
  OHM_coef(nsurf+2,3,1) = SurfaceChar(Gridiv,c_a1_WWet(nsurf+1)) !Snow    a1 Winter wet 
  OHM_coef(1:nsurf,3,2) = SurfaceChar(Gridiv,c_a2_WWet(1:nsurf)) !1:nsurf a2 Winter wet
  OHM_coef(nsurf+2,3,2) = SurfaceChar(Gridiv,c_a2_WWet(nsurf+1)) !Snow    a2 Winter wet 
  OHM_coef(1:nsurf,3,3) = SurfaceChar(Gridiv,c_a3_WWet(1:nsurf)) !1:nsurf a3 Winter wet
  OHM_coef(nsurf+2,3,3) = SurfaceChar(Gridiv,c_a3_WWet(nsurf+1)) !Snow    a3 Winter wet 
  ! Winter dry
  OHM_coef(1:nsurf,4,1) = SurfaceChar(Gridiv,c_a1_WDry(1:nsurf)) !1:nsurf a1 Winter dry
  OHM_coef(nsurf+2,4,1) = SurfaceChar(Gridiv,c_a1_WDry(nsurf+1)) !Snow    a1 Winter dry 
  OHM_coef(1:nsurf,4,2) = SurfaceChar(Gridiv,c_a2_WDry(1:nsurf)) !1:nsurf a2 Winter dry
  OHM_coef(nsurf+2,4,2) = SurfaceChar(Gridiv,c_a2_WDry(nsurf+1)) !Snow    a2 Winter dry 
  OHM_coef(1:nsurf,4,3) = SurfaceChar(Gridiv,c_a3_WDry(1:nsurf)) !1:nsurf a3 Winter dry
  OHM_coef(nsurf+2,4,3) = SurfaceChar(Gridiv,c_a3_WDry(nsurf+1)) !Snow    a3 Winter dry
   
  ! ---- QF coeffs (was in SUEWS_SAHP.f95, subroutine SAHP_Coefs)                 
  BaseTHDD = -999 ! Initialise  QF coeffs 
  QF_A=0
  QF_B=0
  QF_C=0
  AH_min=0
  AH_slope=0
  T_Critic=0
 
  BaseTHDD = SurfaceChar(Gridiv,c_BaseTHDD)
  QF_A = SurfaceChar(Gridiv,(/c_QF_A1,c_QF_A2/))
  QF_B = SurfaceChar(Gridiv,(/c_QF_B1,c_QF_B2/))
  QF_C = SurfaceChar(Gridiv,(/c_QF_C1,c_QF_C2/))
  AH_min = SurfaceChar(Gridiv,c_AHMin) 
  AH_slope = SurfaceChar(Gridiv,c_AHSlope)       
  T_Critic = SurfaceChar(Gridiv,c_TCritic)
  
  ! ---- Irrigation
  Ie_start        	 = SurfaceChar(Gridiv,c_IeStart)
  Ie_end          	 = SurfaceChar(Gridiv,c_IeEnd)
  InternalWaterUse_h  = SurfaceChar(Gridiv,c_IntWU)
  Faut	            	 = SurfaceChar(Gridiv,c_Faut)
  Ie_a		    	 = SurfaceChar(Gridiv,c_Ie_a)
  Ie_m		    	 = SurfaceChar(Gridiv,c_Ie_m)
  DayWat	    	 = SurfaceChar(Gridiv,c_DayWat)
  DayWatPer        	 = SurfaceChar(Gridiv,c_DayWatPer)
    
  ! ---- Hourly profiles
  AHProf(0:23,1) = SurfaceChar(Gridiv,c_HrProfEnUseWD) ! Anthropogenic heat, weekdays
  AHProf(0:23,2) = SurfaceChar(Gridiv,c_HrProfEnUseWE) ! Anthropogenic heat, weekends
  WUProfM(0:23,1)  = SurfaceChar(Gridiv,c_HrProfWUManuWD)  ! Water use, manual, weekdays
  WUProfM(0:23,2)  = SurfaceChar(Gridiv,c_HrProfWUManuWE)  ! Water use, manual, weekends
  WUProfA(0:23,1)  = SurfaceChar(Gridiv,c_HrProfWUAutoWD)  ! Water use, automatic, weekdays
  WUProfA(0:23,2)  = SurfaceChar(Gridiv,c_HrProfWUAutoWE)  ! Water use, automatic, weekends
  SnowProf(0:23,1)  = SurfaceChar(Gridiv,c_HrProfSnowCWD)  ! Snow clearing, weekdays
  SnowProf(0:23,2)  = SurfaceChar(Gridiv,c_HrProfSnowCWE)  ! Snow clearing, weekends
    
  ! ---- Profiles at the resolution of model timestep
  AHProf_tstep(:,1)  = TstepProfiles(Gridiv,cTP_EnUseWD,:) ! Anthropogenic heat, weekdays
  AHProf_tstep(:,2)  = TstepProfiles(Gridiv,cTP_EnUseWE,:) ! Anthropogenic heat, weekends
  WUProfM_tstep(:,1) = TstepProfiles(Gridiv,cTP_WUManuWD,:) ! Water use, manual, weekdays
  WUProfM_tstep(:,2) = TstepProfiles(Gridiv,cTP_WUManuWE,:) ! Water use, manual, weekends
  WUProfA_tstep(:,1) = TstepProfiles(Gridiv,cTP_WUAutoWD,:) ! Water use, automatic, weekdays
  WUProfA_tstep(:,2) = TstepProfiles(Gridiv,cTP_WUAutoWE,:) ! Water use, automatic, weekends
  
  ! ---- Within-grid water distribution
  ! N.B. Rows and columns of WaterDist are the other way round to the input info
  !! Model currently does not include above-ground flow from the Water surface
  !! - Probably should adjust WaterDist to have nsurf columns so that Water can behave like the other surfaces.
  !! Model returns an error if both ToRunoff and ToSoilStore are non-zero. (Now done in CodeMatchDist.) Why?? 
  WaterDist(PavSurf,  1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToPaved(1:(nsurf-1)))
  WaterDist(BldgSurf, 1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToBuilt(1:(nsurf-1)))
  WaterDist(ConifSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToEveTr(1:(nsurf-1)))
  WaterDist(DecidSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToDecTr(1:(nsurf-1)))
  WaterDist(GrassSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToGrass(1:(nsurf-1)))
  WaterDist(BSoilSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToBSoil(1:(nsurf-1)))
  WaterDist(WaterSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToWater(1:(nsurf-1)))
  ! Runoff or SoilStore row !! Change later to allow both Runoff and SoilStore??
  !! N.B. model may not have have been allowing drainage to soilstore - all counted as runoff??
  do iv = 1,(nsurf-1)
     if(SurfaceChar(Gridiv,c_WGToRunoff(iv)) /= 0) then
        WaterDist((nsurf+1),iv) = SurfaceChar(Gridiv,c_WGToRunoff(iv))
     else   
        WaterDist((nsurf+1),iv) = SurfaceChar(Gridiv,c_WGToSoilStore(iv))
     endif   
  enddo
  
  !! ---- Between-grid water distribution
  !!! Need to make these larger than MaxGrids (and recode), as each grid can have 8 connections
  !!GridConnections(1,) = SurfaceChar(Gridiv,c_Grid)
  !!GridConnectionsFrac() = SurfaceChar(Gridiv,55)
  !!GridConnections(2,) = SurfaceChar(Gridiv,54)
  !
  !! Fraction of water from each grid
  !!! N.B. will need to check input files are correctly set up
  !GridToFrac(1:NConns) = (SurfaceChar(Gridiv,55:69:2))
  !! Grid where water goes to
  !GridTo(1:NConns)     = (SurfaceChar(Gridiv,54:68:2))
  !! Come back to this later
  
  
  !write(*,*) 'QFProf'
  !write(*,*) AHProf
  !write(*,*) 'WaterProf'
  !write(*,*) WUProfM
  !write(*,*) WUProfA
  !write(*,*) 'SnowProf'
  !write(*,*) SnowProf
  
  !write(*,*) '------ IRRIGATION ------'
  !write(*,*) Ie_start, Ie_end    
  !write(*,*) InternalWaterUse_h
  !write(*,*) Faut
  !write(*,*) Ie_a
  !write(*,*) Ie_m
  !write(*,*) DayWat
  !write(*,*) DayWatPer
  !write(*,*) WUAreaEveTr_m2,WUAreaDecTr_m2,WUAreaGrass_m2
  !write(*,*) WUAreaTrees_m2
  !write(*,*) '------ ---------- ------'      
  
  !write(*,*) '------ QF COEFFS ------'
  !write(*,*) BaseTHDD
  !write(*,*) QF_A
  !write(*,*) QF_B
  !write(*,*) QF_C
  !write(*,*) AH_min
  !write(*,*) AH_slope
  !write(*,*) T_Critic 
  !write(*,*) '------ --------- ------'        
  
  ! ================================================================================= 
   
  ! Make other assignments using this input information  !! NEEDS CHECKING!
  ! Max albedo for DecTr    
  AlbMax_dec = alb(DecidSurf)
  ! Min & max Storage capacities for DecTr
  CapMin_dec = surf(1,DecidSurf)
  CapMax_dec = surf(5,DecidSurf)
  ! Max GDD & SDD
  GDDMax=0
  SDDMax=0
  do iv=1,nvegsurf
     GDDMax=max(GDDFull(iv),GDDMax)	!! HCW - what is this doing??
     SDDMax=min(SDDFull(iv),SDDMax)
  enddo
  file_qs=.false.  !! HCW - not sure whether this is still needed - what is it for??
       

 !-----------------------------------------------------
  !-----------------------------------------------------
  !NARP_CONFIGURATION if net radiation is to be modelled

  IF(NetRadiationChoice>0)THEN
    NARP_LAT = SurfaceChar(Gridiv,c_lat)
    NARP_LONG = SurfaceChar(Gridiv,c_lng)    ! New sun_position_v2 use degrees FL
    NARP_YEAR = int(SurfaceChar(Gridiv,c_Year))
    NARP_TZ = TIMEZONE                           !not every 5-min
    NARP_EMIS_SNOW = SurfaceChar(Gridiv,c_SnowEmis)
           
    NARP_TRANS_SITE = TRANS_SITE

    !INTERVAL IS ONLY RELEVANT TO LUPCORR
    !ALL OTHER CALCULATIONS ARE INTERVAL INDEPENDENT
    !NB FOR INTERVALS LONGER THAN 15 MINUTES ERRORS IN KCLEAR WILL BE GREATER

    NARP_NPERHOUR=MAX(3600/INTERVAL,1) !!Check this
    IF(ALLOCATED(NARP_KDOWN_HR)) DEALLOCATE(NARP_KDOWN_HR)
    ALLOCATE(NARP_KDOWN_HR(NARP_NPERHOUR))
    NARP_KDOWN_HR=0.

    IF (ldown_option==4.or.ldown_option==5) then !Added by LJ
      !INIITIALIZE SMITH DAY OF YEAR GRID G
     ! NARP_G=SMITHLAMBDA(NINT(LAT))
    ENDIF
  ENDIF


  !When SUEWS_Translate is called from InitialState (ir=0), translate inputs 
!  write(*,*) 'Grid, Year:', SurfaceChar(Gridiv,c_Grid), SurfaceChar(Gridiv,c_Year)
!  write(*,*) 'met forcing file index:',ir
  if(ir == 0) then
     write(*,*) 'This should be seen only when called from InitialState and ir is 0. ir:',ir
     write(*,*) ''
     
     ! =============================================================================
     ! === Translate inputs from ModelDailyState to variable names used in model ===
     ! =============================================================================  
     
     ! Get id_prev from ModelDailyState
     id_prev = ModelDailyState(Gridiv,cMDS_id_prev)
     !write(*,*) 'id_prev',id_prev
     
     porosity = ModelDailyState(Gridiv,cMDS_porosity)
     albDec   = ModelDailyState(Gridiv,cMDS_albDec)
     DecidCap = ModelDailyState(Gridiv,cMDS_DecidCap)
     CumSnowfall = ModelDailyState(Gridiv,cMDS_CumSnowfall)
     
     ! ---- LAI
     lai=0
     lai(id_prev,ivConif)  = ModelDailyState(Gridiv,cMDS_LAIInitialEveTr) 
     lai(id_prev,ivDecid)  = ModelDailyState(Gridiv,cMDS_LAIInitialDecTr)
     lai(id_prev,ivGrass) = ModelDailyState(Gridiv,cMDS_LAIInitialGrass)
               
     ! ---- Growing degree days, GDD
     GDD = 0
     GDD(:,1)=0
     GDD(:,2)=0
     GDD(:,3) = ModelDailyState(Gridiv,cMDS_GDDMin)
     GDD(:,4) = ModelDailyState(Gridiv,cMDS_GDDMax)
     GDD(:,5)=0
     GDD(id_prev,1) = ModelDailyState(Gridiv,cMDS_GDD1_0)
     GDD(id_prev,2) = ModelDailyState(Gridiv,cMDS_GDD2_0)
  
     ! ---- Heating degree days, HDD
     HDD = 0       
     HDD(id_prev,1) = ModelDailyState(Gridiv,cMDS_HDD1)		! 1 = Heating
     HDD(id_prev,2) = ModelDailyState(Gridiv,cMDS_HDD2)		! 2 = Cooling
     HDD(id_prev-3,3) = ModelDailyState(Gridiv,cMDS_TempCOld3)	! 3 will become average 
     HDD(id_prev-2,3) = ModelDailyState(Gridiv,cMDS_TempCOld2)
     HDD(id_prev-1,3) = ModelDailyState(Gridiv,cMDS_TempCOld1)
     HDD(id_prev,3)   = ModelDailyState(Gridiv,cMDS_TempC)
								! 4 = 5 day running mean  
        							! 5 = daily precip total     
     HDD(id_prev,6) = ModelDailyState(Gridiv,cMDS_DaysSinceRain)  ! 6 = days since rain
    
    
     ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
     HDD_grids(:,:,Gridiv) = HDD(:,:) 
     GDD_grids(:,:,Gridiv) = GDD(:,:) 
     LAI_grids(:,:,Gridiv) = LAI(:,:) 
     WU_Day_grids(:,:,Gridiv) = WU_day(:,:)
     AlbDec_grids(:,Gridiv) = AlbDec(:)   
     DecidCap_grids(:,Gridiv) = DecidCap(:) 
     Porosity_grids(:,Gridiv) = Porosity(:)     
    
     ! ---- Snow density of each surface
     DensSnow(1:nsurf) = ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf))
         
     ! =============================================================================
     ! === Translate inputs from ModelOutputData to variable names used in model ===
     ! =============================================================================     
    
     ! ---- Above-ground state
     State(1:nsurf) = ModelOutputData(0,cMOD_State(1:nsurf),Gridiv)    
     ! ---- Below-ground state
     SoilMoist(1:nsurf) = ModelOutputData(0,cMOD_SoilState(1:nsurf),Gridiv)    
     ! ---- Snow fraction
     SnowFrac(1:nsurf)  = ModelOutputData(0,cMOD_SnowFrac(1:nsurf), Gridiv)    
     ! ---- Snow water equivalent in snowpack
     SnowPack(1:nsurf)  = ModelOutputData(0,cMOD_SnowPack(1:nsurf), Gridiv)    
     ! ---- Liquid (melted) water in snowpack
     MeltWaterStore(1:nsurf)  = ModelOutputData(0,cMOD_SnowWaterState(1:nsurf), Gridiv)   
              
  endif  !ir = 0

     
  
  !Write FileChoices.txt (was in SUEWS_Initial.f95) once per grid per year
  
  if (ir==1.and.iMB==1) then
  
  write(*,*) 'Writing to FileChoices for first chunk of met data per year per grid'
  
  FileChoices=trim(FileOutputPath)//trim(FileCode)//'_FileChoices.txt'
  open(12,file=FileChoices,position='append')
  
  write(12,*) ''
  write(12,*) '===================================================================='
  write(12,'(g6.2,i6,g6.2,i6)') 'YEAR:',int(SurfaceChar(Gridiv,c_Year)),'GRID:',int(SurfaceChar(Gridiv,c_Grid))
  
  write(12,*)'--------SUEWS_FunctionalTypes.txt---------------------------- '
  write(12,*)'!Paved  Bldgs   EveTr  DecTr  Grass   BSoil  Water  Snow        -9 not applicable  '
  write(12,120) (alb(iv),iv=1,nsurf), alb_snow,' albedo -1'             ! 1
  write(12,120) (emis(iv),iv=1,nsurf), emis_snow, 'emis -2'               ! 2
  write(12,120) FCskip, FCskip, (baseT(iv),iv=1,nVegsurf),FCskip,FCskip ,'BaseT'  ! 3
  write(12,120) FCskip, FCskip, (baseTe(iv),iv=1,nVegsurf),FCskip,FCskip, 'BaseTe'   ! 4
  write(12,120) (Surf(1,iv),iv=1,nsurf), FCskip ,'storage capacity minimum/default' !5
  write(12,120) (Surf(5,iv),iv=1,nsurf), FCskip ,'storage capacity maximum'       ! 
  write(12,120) (Surf(2,iv),iv=1,nsurf), FCskip ,'drain equation' 
  write(12,120) (Surf(3,iv),iv=1,nsurf), FCskip ,'dr coef1'     ! 7    ! 
  write(12,120) (Surf(4,iv),iv=1,nsurf), FCskip ,'dr coef2'      ! 8    !  
  write(12,'(7f8.0, 2g10.4)') FCskip, FCskip, (GDDFull(iv),iv=1,nVegsurf),FCskip,FCskip,FCskip, 'GDDFull '  ! 10
  write(12,'(7f8.0, 2g10.4)') FCskip, FCskip,(SDDFull(iv),iv=1,nVegsurf),FCskip,FCskip,FCskip,'SDDFull'  ! 11
  write(12,120) FCskip, FCskip, (LAImin(iv),iv=1,nVegsurf),FCskip,FCskip,FCskip,'LAI min'! 12  !
  write(12,120) FCskip, FCskip, (LAImax(iv),iv=1,nVegsurf),FCskip,FCskip,FCskip,'LAI max'  ! 13 ! 
  write(12,120) FCskip, FCskip, (MaxConductance(iv),iv=1,nVegsurf),FCskip,FCskip,FCskip,'MaxCond'  ! 14
  write(12,'(8f8.0, 2g10.4)') (soilstoreCap(iv),iv=1,nsurf), FCskip, 'soilstoreCap'     ! 15
  write(12,120) (VolSoilMoistCap(iv),iv=1,nsurf),FCskip,'VolSoilMoistCap'  ! 16
  write(12,120) (SatHydraulicConduct(iv),iv=1,nsurf),FCskip,'SatHydraulicConduct'  ! 17
  write(12,'(5g10.4, g10.5, g20.0)') G1,G2,G3,G4,G5,G6,' conductance parameters'   ! 17 
  write(12,'(4g10.2,g10.5,g20.0)') TH,TL,S1,S2, Kmax,' conductance parameters'   ! 18 
  write(12,*)  ! FCskip header Soil related
  write(12,'(g12.6,g10.2,f6.0,3g8.3,g8.1)')SoilDensity,SoilDepthMeas, SoilRocks,SmCap,'Soil'  !24
  write(12,*)  ! FCskip header LUMPS related
  write(12,'(10g8.2)')DRAINRT,RAINCOVER,RAINMAXRES,'LUMPS (1)drainage rate,adjust alpha/beta wet surface(3)Max water bucket'  ! 26
  write(12,*)! FCskip header snow related  
  write(12,'(10g8.2)')RadMeltFact,TempMeltFact,albSnowMin,albSnowMax,tau_a,tau_f,PrecipLimitAlb
  write(12,'(10g8.2)')densSnowMin,densSnowMax,tau_r,CRWmin,CRWmax,PrecipLimit    
  write(12,*)! FCskip header NARP related  
  write(12,'(10g8.2)') TRANS_SITE, 'tran_site'
  write(12,*)! FCskip header LAI related  
  write(12,'(10g8.2)') LAItype, (laiPower(iv),iv=1,2)
  
  120	 format (10g10.2)  
  
 
  !! Check what to do about OHM canyons !!
  
  write(12,*) '---Number of rows in OHM_Coefficients.txt =',nlinesOHMCoefficients
  
  write(12,*)'-------','Select OHM','----------------------'
  write(12,'(8g12.4)')  SurfaceChar(Gridiv,c_OHMCode_SWet)
  write(12,'(8g12.4)')  SurfaceChar(Gridiv,c_OHMCode_SDry)
  write(12,'(8g12.4)')  SurfaceChar(Gridiv,c_OHMCode_WWet)
  write(12,'(8g12.4)')  SurfaceChar(Gridiv,c_OHMCode_WDry)
  
  write(12,*)' OHM coefficients----------------------'
  do i=1,4
     do ii=1, nsurf+1
     	write(12,'(2i4,3g10.3)') ii,i, (OHM_coef(ii,i,iii),iii=1,3)
     enddo
  enddo
  
  write(12,*)'----------','Within-grid water distribution','----------'
  !write(12,*)'To !Paved Built  EveTr  DecTr Grass BSoil  Water Runoff/SoilStore'
  do iv=1,(nsurf-1)
    write(12,'(i4, 8f6.2)')iv,(WaterDist(j,iv),j=1,nsurf+1)
  enddo  
  
  write(12,*)'----------','Hourly Profiles','----------'
  write(12,'(24f6.2, a20)') AHProf(0:23,1), 'Hourly anthropogenic heat profile WD'
  write(12,'(24f6.2, a20)') AHProf(0:23,2), 'Hourly anthropogenic heat profile WE'
  write(12,'(24f6.2, a20)') WUProfM(0:23,1),'Hourly manual water use profile WD'  
  write(12,'(24f6.2, a20)') WUProfM(0:23,2),'Hourly manual water use profile WE'  
  write(12,'(24f6.2, a20)') WUProfA(0:23,1),'Hourly automatic water use profile WD'  
  write(12,'(24f6.2, a20)') WUProfA(0:23,2),'Hourly automatic water use profile WE'  
  write(12,'(24f6.2, a20)') SnowProf(0:23,1), 'Hourly snow clearing profile WD'  
  write(12,'(24f6.2, a20)') SnowProf(0:23,2), 'Hourly snow clearing profile WE'  
  
  write(12,*)'----------','Anthropogenic Heat','----------'
  write(12,'(2g12.3)') 'NumCapita = ',NumCapita
  write(12,'(2g12.3)') 'BaseTHDD = ',BaseTHDD
  write(12,'(3g12.3)') 'QF_A = ',QF_A
  write(12,'(3g12.3)') 'QF_B = ',QF_B
  write(12,'(3g12.3)') 'QF_C = ',QF_C  
  write(12,'(2g12.3)') 'AH_Min = ',AH_Min
  write(12,'(2g12.3)') 'AH_Slope = ',AH_Slope
  write(12,'(2g12.3)') 'T_critic = ',T_critic
  
  write(12,*)'----------','Site specific','----------'
  write(12,130) 'Lat = ',lat
  write(12,130) 'Lon = ',lng
  write(12,130) 'SurfaceArea = ',SurfaceArea
  write(12,130) 'RunoffToWater = ',RunoffToWater
  write(12,130) 'WaterUseAreaGrass = ',WUAreaGrass_m2  
  !write(12,130) 'WaterUseAreaTrees = ',WUAreaTrees_m2  
  write(12,130) 'FlowChange = ',FlowChange
  write(12,130) 'PipeCapacity = ',PipeCapacity  
  write(12,130) 'Faut = ',Faut  
  write(12,130) 'Ie_start = ',Ie_start
  write(12,130) 'Ie_end = ',Ie_end
  write(12,130) 'Ie_a = ',Ie_a
  write(12,130) 'Ie_m = ',Ie_m
  write(12,130) 'DayWat = ',DayWat
  write(12,130) 'DayWatPer = ',DayWatPer
  write(12,130) 'InternalWaterUse = ',InternalWaterUse_h
  write(12,130) 'SnowLimBuild = ',SnowLimBuild
  write(12,130) 'SnowLimPaved = ',SnowLimPaved
  
  130 	 format (g20.1,7g10.3) 
  
  write(12,*) '===================================================================='
  
  close(12)

  endif


  !For each row of the met forcing file (ir), translate correct info for each grid 
  ! into model variables
  if (ir>0) then
     ! =============================================================================
     ! === Translate met data from MetForcingData to variable names used in model ==
     ! ============================================================================= 
  
     iy =    int(MetForcingData(ir,1,Gridiv))  !Integer variables
     id =    int(MetForcingData(ir,2,Gridiv))
     it =    int(MetForcingData(ir,3,Gridiv))
     imin =  int(MetForcingData(ir,4,Gridiv))
     qn1_obs =   MetForcingData(ir,5,Gridiv)   !Real values (kind(1d0))
     qh_obs =    MetForcingData(ir,6,Gridiv)
     qe_obs =    MetForcingData(ir,7,Gridiv)
     qs =        MetForcingData(ir,8,Gridiv)
     qf =        MetForcingData(ir,9,Gridiv)
     avu1 =      MetForcingData(ir,10,Gridiv)
     avrh =      MetForcingData(ir,11,Gridiv)
     Temp_C =    MetForcingData(ir,12,Gridiv)
     Press_hPa = MetForcingData(ir,13,Gridiv)
     Precip = MetForcingData(ir,14,Gridiv)
     avkdn =     MetForcingData(ir,15,Gridiv)
     snow_obs =  MetForcingData(ir,16,Gridiv)
     ldown_obs = MetForcingData(ir,17,Gridiv)
     fcld_obs =  MetForcingData(ir,18,Gridiv)
     wu_m3 =     MetForcingData(ir,19,Gridiv)
     xsmd =      MetForcingData(ir,20,Gridiv)
     lai_obs =    MetForcingData(ir,21,Gridiv)
     kdiff =     MetForcingData(ir,22,Gridiv)
     kdir =      MetForcingData(ir,23,Gridiv)
     wdir =      MetForcingData(ir,24,Gridiv)

     ! Calculate dectime
     dectime = real(id,kind(1d0))+real(it,kind(1d0))/24+real(imin,kind(1d0))/(60*24)

     !write(*,*) 'In SUEWS_translate'
     !write(*,*) 'imin',imin             
     !write(*,*) 'it',it             
     !write(*,*) 'id',id
     !write(*,*) 'iy',iy

     ! =============================================================================
     ! === Translate values from ModelDailyState to variable names used in model ===
     ! =============================================================================

     porosity(id) = ModelDailyState(Gridiv,cMDS_porosity)
     albDec(id)   = ModelDailyState(Gridiv,cMDS_albDec)
     DecidCap(id) = ModelDailyState(Gridiv,cMDS_DecidCap)
     CumSnowfall = ModelDailyState(Gridiv,cMDS_CumSnowfall)
     
     ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
     HDD(:,:) = HDD_grids(:,:,Gridiv)
     GDD(:,:) = GDD_grids(:,:,Gridiv)
     LAI(:,:) = LAI_grids(:,:,Gridiv)
     WU_day(:,:) = WU_Day_grids(:,:,Gridiv)
     AlbDec(:)   = AlbDec_grids(:,Gridiv)
     DecidCap(:) = DecidCap_grids(:,Gridiv) 
     Porosity(:) = Porosity_grids(:,Gridiv)
           
     ! ---- Snow density of each surface
     DensSnow(1:nsurf) = ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf))
         
     ! =============================================================================
     ! === Translate values from ModelOutputData to variable names used in model ===
     ! =============================================================================     
    
     ! ---- Above-ground state
     State(1:nsurf) = ModelOutputData(ir-1,cMOD_State(1:nsurf),Gridiv)    
     ! ---- Below-ground state
     SoilMoist(1:nsurf) = ModelOutputData(ir-1,cMOD_SoilState(1:nsurf),Gridiv)    
     ! ---- Snow fraction
     SnowFrac(1:nsurf)  = ModelOutputData(ir-1,cMOD_SnowFrac(1:nsurf), Gridiv)    
     ! ---- Snow water equivalent in snowpack
     SnowPack(1:nsurf)  = ModelOutputData(ir-1,cMOD_SnowPack(1:nsurf), Gridiv)    
     ! ---- Liquid (melted) water in snowpack
     MeltWaterStore(1:nsurf)  = ModelOutputData(ir-1,cMOD_SnowWaterState(1:nsurf), Gridiv)    
   
  endif !ir>0

  return
 end subroutine SUEWS_Translate
 

!Translates model variables to arrays for each grid
! Runs at the end of SUEWS_Calculations to transfer correct info
!for each grid to arrays
!Made by HW Nov 2014
!Last modified: HCW 28 Nov 2014
! To Do:
!===================================================================== 
subroutine SUEWS_TranslateBack(Gridiv,ir,iMB,irMax)

  use Initial
  use allocateArray   ! defines: SiteInfo, sfr, (PavSurf, BldgSurf, etc) 
  use data_in	      ! defines: lat, lng, PopDaytime, PopNighttime, DayLightSavingDay, QF variables
  use sues_data       ! defines: SurfaceArea, IrrFractionTrees, IrrFracConif, IrrFracDecid, IrrFracGrass, Irrigation variables
  use gis_data        ! defines: areaZh, VegFraction, veg_fr, veg_type, BldgH, TreeH, FAIBldg, FAITree, Alt
  use mod_z	      ! defines: z0m, zdm	
  use snowMod	      ! defines: alb_snow, etc		
  use resist	      ! defines: G1-G6, TH, TL, S1, S2, Kmax
  use ohm_calc        ! defines: OHM_coef
  use defaultnotUsed
  use time
  use ColNamesInputFiles
  use ColNamesModelDailyState
  use run_info	      ! for file_qs (still needed??)	 

  IMPLICIT NONE

  integer::Gridiv,&   ! Index of the analysed grid (Gridcounter)
           ir,&       ! Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
           iMB,&      ! Chunk of met data
           irMax      ! Last row in current chunk of met data
           
           
  integer::iv, j
  integer::i, ii, iii 
 
  ! =============================================================================
  ! === Translate values from variable names used in model to ModelDailyState ===
  ! =============================================================================

  ModelDailyState(Gridiv,cMDS_porosity) = porosity(id) 
  ModelDailyState(Gridiv,cMDS_albDec)   = albDec(id)   
  ModelDailyState(Gridiv,cMDS_DecidCap) = DecidCap(id) 
  ModelDailyState(Gridiv,cMDS_CumSnowfall) = CumSnowfall
     
  ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
  HDD_grids(:,:,Gridiv) = HDD(:,:) 
  GDD_grids(:,:,Gridiv) = GDD(:,:) 
  LAI_grids(:,:,Gridiv) = LAI(:,:) 
  WU_Day_grids(:,:,Gridiv) = WU_day(:,:)
  AlbDec_grids(:,Gridiv) = AlbDec(:)   
  DecidCap_grids(:,Gridiv) = DecidCap(:) 
  Porosity_grids(:,Gridiv) = Porosity(:)    

  ! ---- Snow density of each surface
  ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf)) = DensSnow(1:nsurf) 
         
  ! =============================================================================
  ! === Translate values from variable names used in model to ModelOutputData ===
  ! =============================================================================     
    
  ModelOutputData(ir,cMOD_State(1:nsurf),Gridiv) = State(1:nsurf) 
  ModelOutputData(ir,cMOD_SoilState(1:nsurf),Gridiv) = SoilMoist(1:nsurf) 
  ModelOutputData(ir,cMOD_SnowFrac(1:nsurf), Gridiv) = SnowFrac(1:nsurf)  
  ModelOutputData(ir,cMOD_SnowPack(1:nsurf), Gridiv) = SnowPack(1:nsurf)  
  ModelOutputData(ir,cMOD_SnowWaterState(1:nsurf), Gridiv) = MeltWaterStore(1:nsurf)  
  
  if(ir==irMax) then   !Store variables ready for next chunk of met data
     ModelOutputData(0,cMOD_State(1:nsurf),GridCounter) = State(1:nsurf) 
     ModelOutputData(0,cMOD_SoilState(1:nsurf),GridCounter) = SoilMoist(1:nsurf) 
     ModelOutputData(0,cMOD_SnowFrac(1:nsurf),GridCounter) = SnowFrac(1:nsurf)  
     ModelOutputData(0,cMOD_SnowPack(1:nsurf),GridCounter) = SnowPack(1:nsurf)  
     ModelOutputData(0,cMOD_SnowWaterState(1:nsurf),GridCounter) = MeltWaterStore(1:nsurf)  
  endif   
           
  return
endsubroutine SUEWS_TranslateBack
 