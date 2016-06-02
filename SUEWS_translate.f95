!SUEWS_Translate
!Translates - new input arrays (v2014b) to existing model variables
!           - between arrays for different grids and the model variables
!Made by HW&LJ Oct 2014
!-----------------------------------------------------------------------------------
!Last modified: LJ 06 Jul 2015
! Changed to read snowAlb from ModelDailyState instead of SurfaceChar. Location also moved.
!Last modified: HCW 03 Jul 2015
! Use PopDensNighttime by default (not PopDensDaytime)
!Last modified: HCW 26 Jun 2015
! Translation of DailyState variables from the corresponding '_grids' arrays moved
! earlier in code in order to fix bug in DecidCap, AlbDec, Porosity.
!Last modified: HCW 28 Nov 2014
!
! To Do:
!       - Check observed soil moisture works correctly!!
!       - Adjust model to allow water to runoff and sub-surface soil store for each surface type
!	- Adjust model to calculate LAI per surface
!	- Adjust model for SM per surface (measured characteristics)
!===================================================================================
 subroutine SUEWS_Translate(Gridiv,ir,iMB)

  use allocateArray
  use ColNamesInputFiles
  use ColNamesModelDailyState
  use data_in        !defines: lat, lng, PopDaytime, PopNighttime, DayLightSavingDay, QF variables
  use defaultnotUsed
  use gis_data       !defines: areaZh, VegFraction, veg_fr, veg_type, BldgH, TreeH, FAIBldg, FAITree, Alt
  use initial
  use mod_z          !defines: z0m, zdm
  use resist         !defines: G1-G6, TH, TL, S1, S2, Kmax
  use snowMod        !defines: SnowAlb, etc
  use sues_data      !defines: SurfaceArea, IrrFracConif, IrrFracDecid, IrrFracGrass, Irrigation variables
  use time
  
  IMPLICIT NONE

  integer::Gridiv,&   !Index of the analysed grid (Gridcounter)
           ir,&       !Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
           iMB,&      !Chunk of met data
           id_prev

  integer::iv, j
    
  !real (Kind(1d0)):: FCskip = -9   !NULL value used for output to FileChoices	
  real (Kind(1d0)):: FCskip = -999  !NULL value used for output to FileChoices	(changed by HCW 24 May 2016)
  
  real(kind(1d0)):: z0m_in, zdm_in  !Values of z0m and zdm provided in SiteSelect input file (do not get updated unlike z0d and z0m)
  
  character(len=20):: grid_txt
  character(len=4):: year_txt
  character(len=12)::SsG_YYYY !Site, grid, year string 
  
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
  sfr(BldgSurf)  = SurfaceChar(Gridiv,c_FrBldgs)  ! Bldgs
  sfr(ConifSurf) = SurfaceChar(Gridiv,c_FrEveTr)  ! Everg
  sfr(DecidSurf) = SurfaceChar(Gridiv,c_FrDecTr)  ! Decid
  sfr(GrassSurf) = SurfaceChar(Gridiv,c_FrGrass)  ! Grass
  sfr(BSoilSurf) = SurfaceChar(Gridiv,c_FrBSoil)  ! BSoil
  sfr(WaterSurf) = SurfaceChar(Gridiv,c_FrWater)  ! Water

  ! Check the surface fractions add up to 1 (or close to 1)
  if(sum(sfr)>1.001.or.sum(sfr)<0.999) call ErrorHint(10,'SiteSelect.txt - check surface fractions',sum(sfr),notUsed,notUsedI)
  
  ! ---- Irrigated fractions
  IrrFracConif = SurfaceChar(Gridiv,c_IrrEveTrFrac)  ! Everg
  IrrFracDecid = SurfaceChar(Gridiv,c_IrrDecTrFrac)  ! Decid
  IrrFracGrass = SurfaceChar(Gridiv,c_IrrGrassFrac)  ! Grass
  
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
  
  ImpervFraction =  (sfr(PavSurf) + sfr(BldgSurf))
  PervFraction =  (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf) + sfr(WaterSurf)) 
  NonWaterFraction = (sfr(PavSurf) + sfr(BldgSurf) + sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf)) 
  ! ---------------------------------------------------------------------------------

  ! ---- Heights & frontal areas
  BldgH = SurfaceChar(Gridiv,c_HBldgs)        ! Building height [m]
  EveTreeH = SurfaceChar(Gridiv,c_HEveTr)     ! Evergreen tree height [m]
  DecTreeH = SurfaceChar(Gridiv,c_HDecTr)     ! Deciduous tree height [m]
  TreeH = (EveTreeH*sfr(ConifSurf) + DecTreeH*sfr(DecidSurf))/(sfr(ConifSurf)+sfr(DecidSurf)) ! Average tree height [m]
  FAIBldg = SurfaceChar(Gridiv,c_FAIBldgs)    ! Frontal area index for buildings
  FAIEveTree = SurfaceChar(Gridiv,c_FAIEveTr) ! Frontal area index for evergreen trees
  FAIDecTree = SurfaceChar(Gridiv,c_FAIDecTr) ! Frontal area index for deciduous trees
  FAITree = (FAIEveTree*sfr(ConifSurf) + FAIDecTree*sfr(DecidSurf))/(sfr(ConifSurf)+sfr(DecidSurf)) ! Frontal area index for trees
  
  z0m = SurfaceChar(Gridiv,c_z0m)     ! Roughness length [m]
  zdm = SurfaceChar(Gridiv,c_zdm)     ! Displacement height [m]
  ! z0m and zdm can vary in time depending on z0method selected. Save the input values here
  z0m_in = z0m
  zdm_in = zdm
  
  ! ---- Population
  PopDensDaytime   = SurfaceChar(Gridiv,c_PopDensDay)    ! Daytime population density [ha-1]
  PopDensNighttime = SurfaceChar(Gridiv,c_PopDensNight)  ! Night-time population density [ha-1]
  
  NumCapita = PopDensNighttime      ! Pop density [ha-1]     !!Use Night-time pop density for NumCapita for testing! 
  
  ! ---- Albedo [-]
  alb(1:nsurf) = SurfaceChar(Gridiv,c_AlbMax)   !Use maximum albedos as default value (albmin for veg surfaces handled below)     

  ! ---- Set min & max albedo for vegetated surfaces (min albedo not currently used for NonVeg or Water surfaces)
  AlbMin_EveTr = SurfaceChar(Gridiv,c_AlbMin(ConifSurf))
  AlbMax_EveTr = SurfaceChar(Gridiv,c_AlbMax(ConifSurf))
  AlbMin_DecTr = SurfaceChar(Gridiv,c_AlbMin(DecidSurf))
  AlbMax_DecTr = SurfaceChar(Gridiv,c_AlbMax(DecidSurf))
  AlbMin_Grass = SurfaceChar(Gridiv,c_AlbMin(GrassSurf))
  AlbMax_Grass = SurfaceChar(Gridiv,c_AlbMax(GrassSurf))
  
  ! ---- Emissivity [-] 
  emis(1:nsurf) = SurfaceChar(Gridiv,c_Emis)
  emis_snow = SurfaceChar(Gridiv,c_SnowEmis)
  
  ! ---- Storage capacities [mm]
  surf(1,1:nsurf) = SurfaceChar(Gridiv,c_StorMin)   ! Minimum	
  surf(5,1:nsurf) = SurfaceChar(Gridiv,c_StorMax)   ! Maximum	
  surf(6,1:nsurf) = surf(1,1:nsurf)  !Set storage capacities for all surface to minimum (DecTr changes with time in Calculations).

  ! ---- Set min & max storage capacities for DecTr
  CapMin_dec = surf(1,DecidSurf)
  CapMax_dec = surf(5,DecidSurf)
  ! ---- Set min & max porosity for DecTr
  PorMin_dec = 0.2
  PorMax_dec = 0.6
  
  ! ---- Threshold for wet evaporation [mm]
  WetThresh(1:nsurf) = SurfaceChar(Gridiv,c_WetThresh)
  
  ! ---- Limit for state [mm]
  StateLimit(1:nsurf) = SurfaceChar(Gridiv,c_StateLimit)
  
  ! ---- Drainage
  surf(2,1:nsurf) = SurfaceChar(Gridiv,c_DrEq)      ! Drainage equation
  surf(3,1:nsurf) = SurfaceChar(Gridiv,c_DrCoef1)   ! Drainage coef 1
  surf(4,1:nsurf) = SurfaceChar(Gridiv,c_DrCoef2)   ! Drainage coef 2
   
  ! ---- Limit of SWE (each surface except Water)
  snowD(1:(nsurf-1)) = SurfaceChar(Gridiv,c_SnowLimPat(1:(nsurf-1)))
  
  ! ---- Snow limit for removal (only impervious surfaces)
  SnowLimPaved = SurfaceChar(Gridiv,c_SnowLimRem(PavSurf)) 
  SnowLimBuild = SurfaceChar(Gridiv,c_SnowLimRem(BldgSurf))
  !SnowLimBSoil = SurfaceChar(Gridiv,c_SnowLimRem(BSoilSurf))   !Snow clearing not applicable to bare soil surface
  
  ! ---- Soil characteristics (each surface except Water)
  SoilDepth(1:(nsurf-1))    = SurfaceChar(Gridiv,c_SoilDepth(1:(nsurf-1))) ! Depth of sub-surface soil store [mm]
  SoilStoreCap(1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilStCap(1:(nsurf-1))) ! Soil store capacity [mm]
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
  SoilDensity    = SurfaceChar(Gridiv,c_SoilDens(1))   !!Not sure this works correctly - need to check
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
  !! Set all rows the same for now in input file & use Decid row in model
  LAItype = int(SurfaceChar(Gridiv,c_LAIEq(ivDecid)))
  LAIPower(1) = SurfaceChar(Gridiv,c_LeafGP1(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(2) = SurfaceChar(Gridiv,c_LeafGP2(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(3) = SurfaceChar(Gridiv,c_LeafOP1(ivDecid)) ! 4 powers (not 4 veg types!)
  LAIPower(4) = SurfaceChar(Gridiv,c_LeafOP2(ivDecid)) ! 4 powers (not 4 veg types!)  
  
  ! ---- LUMPS-related parameters
  DRAINRT    = SurfaceChar(Gridiv,c_LUMPSDr)      ! LUMPS Drainage rate [mm h-1]
  RAINCOVER  = SurfaceChar(Gridiv,c_LUMPSCover)   ! LUMPS Limit when surface totally wet [mm]
  RAINMAXRES = SurfaceChar(Gridiv,c_LUMPSMaxRes)  ! LUMPS Maximum water bucket reservoir [mm]
  
  ! ---- NARP-related parameters
  TRANS_SITE = SurfaceChar(Gridiv,c_NARPTrans)    ! NARP atmospheric transmissivity
   
  ! ---- Snow-related characteristics
  RadMeltFact    = SurfaceChar(Gridiv,c_SnowRMFactor)  
  TempMeltFact   = SurfaceChar(Gridiv,c_SnowTMFactor)  
  SnowAlbMin     = SurfaceChar(Gridiv,c_SnowAlbMin)
  SnowAlbMax     = SurfaceChar(Gridiv,c_SnowAlbMax)
  tau_a          = SurfaceChar(Gridiv,c_Snowtau_a)
  tau_f          = SurfaceChar(Gridiv,c_Snowtau_f)
  PrecipLimitAlb = SurfaceChar(Gridiv,c_SnowPlimAlb) 
  SnowDensMin    = SurfaceChar(Gridiv,c_SnowSDMin)
  SnowDensMax    = SurfaceChar(Gridiv,c_SnowSDMax)
  tau_r          = SurfaceChar(Gridiv,c_Snowtau_r)
  CRWMin         = SurfaceChar(Gridiv,c_SnowCRWMin)
  CRWMax         = SurfaceChar(Gridiv,c_SnowCRWMax)
  PrecipLimit    = SurfaceChar(Gridiv,c_SnowPLimSnow)

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
  DayLightSavingDay(1) = int(SurfaceChar(Gridiv,c_StartDLS))
  DayLightSavingDay(2) = int(SurfaceChar(Gridiv,c_EndDLS))
                  
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
  Ie_start           = int(SurfaceChar(Gridiv,c_IeStart))
  Ie_end             = int(SurfaceChar(Gridiv,c_IeEnd))
  InternalWaterUse_h = SurfaceChar(Gridiv,c_IntWU)
  Faut               = SurfaceChar(Gridiv,c_Faut)
  Ie_a = SurfaceChar(Gridiv,c_Ie_a)   !Automatic irrigation model coefficients [mm d-1]; [mm d-1 degC-1]; [mm d-2]
  Ie_m = SurfaceChar(Gridiv,c_Ie_m)   !Manual irrigation model coefficients [mm d-1]; [mm d-1 degC-1]; [mm d-2]
  DayWat             = SurfaceChar(Gridiv,c_DayWat)
  DayWatPer          = SurfaceChar(Gridiv,c_DayWatPer)
    
  ! ---- Hourly profiles
  AHProf(0:23,1)   = SurfaceChar(Gridiv,c_HrProfEnUseWD)   ! Anthropogenic heat, weekdays
  AHProf(0:23,2)   = SurfaceChar(Gridiv,c_HrProfEnUseWE)   ! Anthropogenic heat, weekends
  WUProfM(0:23,1)  = SurfaceChar(Gridiv,c_HrProfWUManuWD)  ! Water use, manual, weekdays
  WUProfM(0:23,2)  = SurfaceChar(Gridiv,c_HrProfWUManuWE)  ! Water use, manual, weekends
  WUProfA(0:23,1)  = SurfaceChar(Gridiv,c_HrProfWUAutoWD)  ! Water use, automatic, weekdays
  WUProfA(0:23,2)  = SurfaceChar(Gridiv,c_HrProfWUAutoWE)  ! Water use, automatic, weekends
  SnowProf(0:23,1) = SurfaceChar(Gridiv,c_HrProfSnowCWD)   ! Snow clearing, weekdays
  SnowProf(0:23,2) = SurfaceChar(Gridiv,c_HrProfSnowCWE)   ! Snow clearing, weekends
    
  ! ---- Profiles at the resolution of model time step
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
  ! Model returns an error if both ToRunoff and ToSoilStore are non-zero (in CodeMatchDist) 
  ! For impervious surfaces, water goes to runoff; for pervious surfaces, water goes to soilstore
  WaterDist(PavSurf,  1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToPaved(1:(nsurf-1)))
  WaterDist(BldgSurf, 1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToBldgs(1:(nsurf-1)))
  WaterDist(ConifSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToEveTr(1:(nsurf-1)))
  WaterDist(DecidSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToDecTr(1:(nsurf-1)))
  WaterDist(GrassSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToGrass(1:(nsurf-1)))
  WaterDist(BSoilSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToBSoil(1:(nsurf-1)))
  WaterDist(WaterSurf,1:(nsurf-1)) = SurfaceChar(Gridiv,c_WGToWater(1:(nsurf-1)))
  ! Runoff or SoilStore row   !!Change later to allow both Runoff and SoilStore
  do iv = 1,(nsurf-1)
     if(SurfaceChar(Gridiv,c_WGToRunoff(iv)) /= 0) then
        WaterDist((nsurf+1),iv) = SurfaceChar(Gridiv,c_WGToRunoff(iv))
     else   
        WaterDist((nsurf+1),iv) = SurfaceChar(Gridiv,c_WGToSoilStore(iv))
     endif   
  enddo

  ! Access required DailyState variables for the current grid (moved HCW 26 Jun 2015)
  HDD(:,:)    = HDD_grids(:,:,Gridiv)
  GDD(:,:)    = GDD_grids(:,:,Gridiv)
  LAI(:,:)    = LAI_grids(:,:,Gridiv)
  WU_day(:,:) = WU_Day_grids(:,:,Gridiv)
  AlbDecTr(:) = AlbDecTr_grids(:,Gridiv)
  DecidCap(:) = DecidCap_grids(:,Gridiv) 
  Porosity(:) = Porosity_grids(:,Gridiv)
  AlbEveTr(:) = AlbEveTr_grids(:,Gridiv)
  AlbGrass(:) = AlbGrass_grids(:,Gridiv)
  SnowAlb = ModelDailyState(Gridiv,cMDS_SnowAlb)
  
  !! ---- Between-grid water distribution
  !!! Need to make these larger than MaxNumberOfGrids (and recode), as each grid can have 8 connections
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
   
  ! ================================================================================= 
   
        
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

    ! Commented out HCW 04 Mar 2015
    !NARP_NPERHOUR=MAX(3600/t_INTERVAL,1) !!Check this
    !IF(ALLOCATED(NARP_KDOWN_HR)) DEALLOCATE(NARP_KDOWN_HR)
    !ALLOCATE(NARP_KDOWN_HR(NARP_NPERHOUR))
    !NARP_KDOWN_HR=0.

    !IF (ldown_option==4.or.ldown_option==5) then !Added by LJ
    !  INIITIALIZE SMITH DAY OF YEAR GRID G
    !  NARP_G=SMITHLAMBDA(NINT(LAT))
    !ENDIF
  ENDIF


  !=================================================================================
  ! When SUEWS_Translate is called from InitialState (ir=0), inputs need translating
  if(ir == 0) then
     !write(*,*) 'This should be seen only when called from InitialState and ir is 0. ir:',ir
          
     ! =============================================================================
     ! === Translate inputs from ModelDailyState to variable names used in model ===
     ! =============================================================================  
     
     ! Get id_prev from ModelDailyState
     id_prev = int(ModelDailyState(Gridiv,cMDS_id_prev))
         
     porosity = ModelDailyState(Gridiv,cMDS_porosity)
     albDecTr   = ModelDailyState(Gridiv,cMDS_albDecTr)
     albEveTr   = ModelDailyState(Gridiv,cMDS_albEveTr)
     albGrass   = ModelDailyState(Gridiv,cMDS_albGrass)
     DecidCap = ModelDailyState(Gridiv,cMDS_DecidCap)
     CumSnowfall = ModelDailyState(Gridiv,cMDS_CumSnowfall)
     SnowAlb    = ModelDailyState(Gridiv,cMDS_SnowAlb)

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
     HDD(id_prev,1) = ModelDailyState(Gridiv,cMDS_HDD1)         ! 1 = Heating
     HDD(id_prev,2) = ModelDailyState(Gridiv,cMDS_HDD2)         ! 2 = Cooling
     HDD(id_prev-3,3) = ModelDailyState(Gridiv,cMDS_TempCOld3)  ! 3 will become average
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
     AlbDecTr_grids(:,Gridiv) = AlbDecTr(:)   
     AlbEveTr_grids(:,Gridiv) = AlbEveTr(:)   
     AlbGrass_grids(:,Gridiv) = AlbGrass(:)   
     DecidCap_grids(:,Gridiv) = DecidCap(:) 
     Porosity_grids(:,Gridiv) = Porosity(:)

     ! ---- Snow density of each surface
     SnowDens(1:nsurf) = ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf))
         
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
  !=================================================================================
  
  
  !=========================== Write FileChoices.txt ===============================
  !=================================================================================
  ! Do once per grid per year (was in SUEWS_Initial.f95)
  if (ir==1.and.iMB==1) then   !For first row of first block only
     !write(*,*) 'Writing to FileChoices for first chunk of met data per year per grid'
     FileChoices=trim(FileOutputPath)//trim(FileCode)//'_FileChoices.txt'
     open(12,file=FileChoices,position='append')
      
     write(grid_txt,'(I5)') int(SurfaceChar(Gridiv,c_Grid))
     write(year_txt,'(I4)') int(SurfaceChar(Gridiv,c_Year)) 
     write(SsG_YYYY,'(A12)') trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(adjustl(year_txt))
     
     !write(12,*) '--------------------------------------------------------------------------------'
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Surface characteristics'//' -----'
     ! Characteristics that apply to some or all surface types 
     write(12,'(8a10,a16)') 'Paved','Bldgs','EveTr','DecTr','Grass','BSoil','Water','Snow', ' SurfType' 
     write(12,120) (sfr(iv),iv=1,nsurf),FCskip, ' SurfFr'
     write(12,120) FCskip,FCskip,IrrFracConif,IrrFracDecid,IrrFracGrass,FCskip,FCskip,FCskip, ' IrrFr'
     write(12,120) FCskip,FCskip,WUAreaEveTr_m2,WUAreaDecTr_m2,WUAreaGrass_m2,FCskip,FCskip,FCskip, ' WaterUseArea'
     write(12,120) FCskip,BldgH,EveTreeH,DecTreeH,FCskip,FCskip,FCskip,FCskip, ' H'
     write(12,120) FCskip,FAIBldg,FAIEveTree,FAIDecTree,FCskip,FCskip,FCskip,FCskip, ' FAI'
     write(12,120) FCskip,FCskip,AlbMin_EveTr,AlbMin_DecTr,AlbMin_Grass,FCskip,FCskip,SnowAlbMin, ' AlbedoMin'
     write(12,120) FCskip,FCskip,AlbMax_EveTr,AlbMax_DecTr,AlbMax_Grass,FCskip,FCskip,SnowAlbMax, ' AlbedoMax'
     !write(12,120) (alb(iv),iv=1,nsurf),SnowAlb, ' Albedo'   ! This is instantaneous value (not provided as input)
     write(12,120) (emis(iv),iv=1,nsurf),emis_snow, ' Emissivity'
     write(12,120) FCskip, FCskip,(baseT(iv),iv=1,nvegsurf), FCskip, FCskip, FCskip, ' BaseT' 
     write(12,120) FCskip, FCskip,(baseTe(iv),iv=1,nvegsurf),FCskip, FCskip, FCskip, ' BaseTe'    
     write(12,120) (Surf(1,iv),iv=1,nsurf), FCskip ,' StorageMin' 
     write(12,120) (Surf(5,iv),iv=1,nsurf), FCskip ,' StorageMax'      
     write(12,120) (WetThresh(iv),iv=1,nsurf), FCskip, ' WetThreshold' 
     write(12,120) (StateLimit(iv),iv=1,nsurf), FCskip,' StateLimit' 
     write(12,120) (Surf(2,iv),iv=1,nsurf), FCskip, ' DrainageEq'    !real
     write(12,120) (Surf(3,iv),iv=1,nsurf), FCskip, ' DrainageCoef1'    
     write(12,120) (Surf(4,iv),iv=1,nsurf), FCskip, ' DrainageCoef2'      
     write(12,120) FCskip,FCskip,(GDDFull(iv),iv=1,nvegsurf),FCskip,FCskip,FCskip, ' GDDFull'  
     write(12,120) FCskip,FCskip,(SDDFull(iv),iv=1,nvegsurf),FCskip,FCskip,FCskip, ' SDDFull'  
     write(12,120) FCskip,FCskip,(LAImin(iv),iv=1,nvegsurf),FCskip,FCskip,FCskip, ' LAIMin'    
     write(12,120) FCskip,FCskip,(LAImax(iv),iv=1,nvegsurf),FCskip,FCskip,FCskip, ' LAIMax'    
     write(12,'(3f10.3,i10,  4f10.3,a16)') FCskip,FCskip,FCskip,LAItype,FCskip,FCskip,FCskip,FCskip, ' LAIEq'   !integer
     write(12,'(3f10.3,f10.5,4f10.3,a16)') FCskip,FCskip,FCskip,LAIPower(1),FCskip,FCskip,FCskip,FCskip, ' LAI_LeafGP1'
     write(12,'(3f10.3,f10.5,4f10.3,a16)') FCskip,FCskip,FCskip,LAIPower(2),FCskip,FCskip,FCskip,FCskip, ' LAI_LeafGP2'
     write(12,'(3f10.3,f10.5,4f10.3,a16)') FCskip,FCskip,FCskip,LAIPower(3),FCskip,FCskip,FCskip,FCskip, ' LAI_LeafOP1'
     write(12,'(3f10.3,f10.5,4f10.3,a16)') FCskip,FCskip,FCskip,LAIPower(4),FCskip,FCskip,FCskip,FCskip, ' LAI_LeafOP2'
     write(12,120) FCskip,FCskip,(MaxConductance(iv),iv=1,nvegsurf),FCskip,FCskip,FCskip, ' MaxCond' 
     write(12,120) (SoilDepth(iv),iv=1,(nsurf-1)),FCskip,FCskip, ' SoilDepth' 
     write(12,120) (soilstoreCap(iv),iv=1,(nsurf-1)), FCskip, FCskip, ' SoilStoreCap'
     write(12,'(6f10.5,2f10.3,a16)') (SatHydraulicConduct(iv),iv=1,(nsurf-1)),FCskip,FCskip, ' SatHydraulicConduct'
     ! Not currently coded, but add these later: SoilDensity, InfiltrationRate, OBS_SMDept, OBS_SMCap, OBS_SoilNotRocks
     write(12,120) (snowD(iv),iv=1,(nsurf-1)),FCskip,FCskip, ' SnowLimPatch'    
     write(12,120) SnowLimPaved,SnowLimBuild,FCskip,FCskip,FCskip,FCskip,FCskip,FCskip, ' SnowLimRemove'  
     write(12,120) (OHM_coef(1:nsurf,1,1)),OHM_coef(nsurf+2,1,1), ' OHM_a1_Sum_Wet'
     write(12,120) (OHM_coef(1:nsurf,2,1)),OHM_coef(nsurf+2,2,1), ' OHM_a1_Sum_Dry'    
     write(12,120) (OHM_coef(1:nsurf,3,1)),OHM_coef(nsurf+2,3,1), ' OHM_a1_Win_Wet'    
     write(12,120) (OHM_coef(1:nsurf,4,1)),OHM_coef(nsurf+2,4,1), ' OHM_a1_Win_Dry'    
     write(12,120) (OHM_coef(1:nsurf,1,2)),OHM_coef(nsurf+2,1,2), ' OHM_a2_Sum_Wet'
     write(12,120) (OHM_coef(1:nsurf,2,2)),OHM_coef(nsurf+2,2,2), ' OHM_a2_Sum_Dry'
     write(12,120) (OHM_coef(1:nsurf,3,2)),OHM_coef(nsurf+2,3,2), ' OHM_a2_Win_Wet'
     write(12,120) (OHM_coef(1:nsurf,4,2)),OHM_coef(nsurf+2,4,2), ' OHM_a2_Win_Dry'
     write(12,120) (OHM_coef(1:nsurf,1,3)),OHM_coef(nsurf+2,1,3), ' OHM_a3_Sum_Wet'        
     write(12,120) (OHM_coef(1:nsurf,2,3)),OHM_coef(nsurf+2,2,3), ' OHM_a3_Sum_Dry'       
     write(12,120) (OHM_coef(1:nsurf,3,3)),OHM_coef(nsurf+2,3,3), ' OHM_a3_Win_Wet'  
     write(12,120) (OHM_coef(1:nsurf,4,3)),OHM_coef(nsurf+2,4,3), ' OHM_a3_Win_Dry'
         
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Snow parameters'//' -----'
     write(12,'(a12,11a10)') 'Grid','RadMeltF','TempMeltF','tau_a','tau_f','PLimAlb','SDensMin','SDensMax', &
                         'tau_r','CRWMin','CRWMax','PLimSnow'
     write(12,'(a12,11f10.4)') SsG_YYYY,RadMeltFact,TempMeltFact,tau_a,tau_f,PrecipLimitAlb,SnowDensMin,SnowDensMax, &
                           tau_r,CRWmin,CRWmax,PrecipLimit    

     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Conductance parameters'//' -----'
     write(12,'(a12,11a10)') 'Grid','G1','G2','G3','G4','G5','G6','TH','TL','S1','S2','Kmax'
     write(12,'(a12,11f10.3)') SsG_YYYY,G1,G2,G3,G4,G5,G6,TH,TL,S1,S2,Kmax
     
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Energy-use parameters'//' -----'
     write(12,'(a12,11a10)') 'Grid','NumCapita','BaseTHDD','QF_A_WD','QF_A_WE','QF_B_WD','QF_B_WE','QF_C_WD','QF_C_WE', & 
                         'AH_Min','AH_Slope','T_critic'
     write(12,'(a12,11f10.3)') SsG_YYYY,NumCapita,BaseTHDD,QF_A(1:2),QF_B(1:2),QF_C(1:2), &
                           AH_Min,AH_Slope,T_critic

     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Water-use parameters'//' -----'
     write(12,'(a12,10a10)') 'Grid','IeStart','IeEnd','IntWatUse','Faut', &
                         'Ie_a1','Ie_a2','Ie_a3','Ie_m1','Ie_m2','Ie_m3'
     write(12,'(a12,2i10,8f10.3)') SsG_YYYY,Ie_start,Ie_end,InternalWaterUse_h,Faut, &
                               Ie_a(1:3),Ie_m(1:3)
     
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Weekly profiles'//' -----'                          
     write(12,'(a12,7a10,  a16)') 'Grid','1_Sun','2_Mon','3_Tue','4_Wed','5_Thu','6_Fri','7_Sat', ' DayOfWeek'
     write(12,'(a12,7f10.3,a16)') SsG_YYYY,DayWat(1:7), ' Irr allowed'
     write(12,'(a12,7f10.3,a16)') SsG_YYYY,DayWatPer(1:7), ' Frac properties'
     
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Hourly profiles'//' -----'
     write(12,'(a12,24i10,a20)') 'Grid',(iv,iv=0,23),'HourOfDay'
     write(12,121) SsG_YYYY,AHProf(0:23,1), ' Anthropogenic heat WD'
     write(12,121) SsG_YYYY,AHProf(0:23,2), ' Anthropogenic heat  WE'
     write(12,121) SsG_YYYY,WUProfM(0:23,1),' Manual water use WD'  
     write(12,121) SsG_YYYY,WUProfM(0:23,2),' Manual water use WE'  
     write(12,121) SsG_YYYY,WUProfA(0:23,1),' Automatic water use WD'  
     write(12,121) SsG_YYYY,WUProfA(0:23,2),' Automatic water use WE'  
     write(12,121) SsG_YYYY,SnowProf(0:23,1), ' Snow clearing WD'  
     write(12,121) SsG_YYYY,SnowProf(0:23,2), ' Snow clearing WE'  
         
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Within-grid water distribution'//' -----'
     write(12,'(9a10)') 'ToPaved','ToBldgs','ToEveTr','ToDecTr','ToGrass','ToBSoil','ToWater','ToROorSS'
     do iv=1,(nsurf-1)
        write(12,'(8f10.4)')(WaterDist(j,iv),j=1,nsurf+1)
     enddo              
                 
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Other parameters'//' -----'
     write(12,'(a12,7a10)') 'Grid','FlowChange','ROToWater','PipeCap', &   ! Water-related
                            'DrRate','Cover','MaxRes',&   ! LUMPS-related
                            'Trans'   ! NARP-related
     write(12,'(a12,7f10.3)') SsG_YYYY,FlowChange,RunoffToWater,PipeCapacity, &
                              DRAINRT,RAINCOVER,RAINMAXRES, &
                              Trans_Site
    
     write(12,*) '----- '//trim(adjustl(SsG_YYYY))//' Site parameters'//' -----'    
     write(12,'(a12,9a10)') 'Grid','lat','lon','alt','SurfA_ha','NumCapita','z0_input','zd_input','StartDLS','EndDLS'     
     write(12,'(a12,3f10.4,f10.2,3f10.4,2i10)') SsG_YYYY,lat,lng,alt,SurfaceArea_ha,NumCapita,z0m_in,zdm_in,DayLightSavingDay(1:2)     
     
     write(12,*) ''
     
     120  format (8f10.3, a16)  !format (10g10.2)
     121  format (a12,24f10.4, a20)
     
     close(12)
      
     
     !==============================================================================
     ! Check input values are reasonable ===========================================
    
     ! Coefficients for anthropogenic heat models ----------------------------------
     if(AnthropHeatChoice==1) then   !Loridan et al. (2011) calculation
        if(AH_min==0.and.Ah_slope==0.and.T_Critic==0) then
           call ErrorHint(53,'Check QF calculation coefficients.',notUsed,notUsed,AnthropHeatChoice)
        endif
  
     elseif(AnthropHeatChoice==2) then   !Jarvi et al. (2011) calculation 
        if(sum(QF_A)==0.and.sum(QF_B)==0.and.sum(QF_C)==0) then
           call ErrorHint(54,'Check QF calculation coefficients.',notUsed,notUsed,AnthropHeatChoice)
        endif
     endif 

     ! Morphometric parameters -----------------------------------------------------
     if(z0_method==1) then   !z0, zd values provided in input file
        if(z0m<0.00001) call ErrorHint(5,'Check value of z0 in input data.',z0m,notUsed,notUsedI)
        if(zdm<0.00001) call ErrorHint(6,'Check value of zd in input data.',zdm,notUsed,notUsedI)
        zzd=z-zdm
     elseif(z0_method==3) then   !z0, zd calculated using FAI in input file, check FAIs reasonable
        if(FAIBLdg<0) call ErrorHint(1,'Check value of FAI for Buildings in input data',FAIBldg,notUsed,notUsedI)
        if(FAITree<0) call ErrorHint(2,'Check value of FAI for Trees in input data (weighted EveTr+DecTr)',FAITree,notUsed,notUsedI)
     endif 
         
  endif   !End for first row of first block only ===================================
     
  
  !=================================================================================
  !For each row of the met forcing file (ir), translate correct info for each grid 
  ! into model variables
  if (ir>0) then
     ! =============================================================================
     ! === Translate met data from MetForcingData to variable names used in model ==
     ! ============================================================================= 
     iy =    int(MetForcingData(ir, 1,Gridiv))  !Integer variables
     id =    int(MetForcingData(ir, 2,Gridiv))
     it =    int(MetForcingData(ir, 3,Gridiv))
     imin =  int(MetForcingData(ir, 4,Gridiv))
     qn1_obs =   MetForcingData(ir, 5,Gridiv)   !Real values (kind(1d0))
     qh_obs =    MetForcingData(ir, 6,Gridiv)
     qe_obs =    MetForcingData(ir, 7,Gridiv)
     qs =        MetForcingData(ir, 8,Gridiv)
     qf =        MetForcingData(ir, 9,Gridiv)
     avu1 =      MetForcingData(ir,10,Gridiv)
     avrh =      MetForcingData(ir,11,Gridiv)
     Temp_C =    MetForcingData(ir,12,Gridiv)
     Press_hPa = MetForcingData(ir,13,Gridiv)
     Precip =    MetForcingData(ir,14,Gridiv)
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

     ! =============================================================================
     ! === Translate values from ModelDailyState to variable names used in model ===
     ! =============================================================================
     porosity(id) = ModelDailyState(Gridiv,cMDS_porosity)
     albDecTr(id) = ModelDailyState(Gridiv,cMDS_albDecTr)
     albEveTr(id) = ModelDailyState(Gridiv,cMDS_albEveTr)
     albGrass(id) = ModelDailyState(Gridiv,cMDS_albGrass)
     DecidCap(id) = ModelDailyState(Gridiv,cMDS_DecidCap)
     CumSnowfall  = ModelDailyState(Gridiv,cMDS_CumSnowfall)
     ! ---- Snow density of each surface
     SnowDens(1:nsurf) = ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf))
         
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
    
  endif !ir>0   !===================================================================
  
  ! --------------------------------------------------------------------------------
  ! Check Initial Conditions are reasonable ----------------------------------------
  if (ir==1.and.iMB==1) then   !For first row of first block only
     call CheckInitial
  endif   
  ! --------------------------------------------------------------------------------
    
  return
 end subroutine SUEWS_Translate
!===================================================================================
 
 !SUEWS_TranslateBack
!Translates model variables to arrays for each grid
!Runs at the end of SUEWS_Calculations to store correct info for each grid
!Made by HW Nov 2014
!-----------------------------------------------------------------------------------
!Last modified:LJ 14 Sep 2015
!              HCW 28 Nov 2014
!===================================================================================
 subroutine SUEWS_TranslateBack(Gridiv,ir,irMax)

  use allocateArray   
  use ColNamesInputFiles
  use ColNamesModelDailyState
  use data_in
  use defaultnotUsed
  use gis_data
  use Initial
  use mod_z
  use resist
  use snowMod
  use sues_data
  use time
  
  IMPLICIT NONE

  integer::Gridiv,&   ! Index of the analysed grid (Gridcounter)
           ir,&       ! Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
           irMax      ! Last row in current chunk of met data
           
  ! =============================================================================
  ! === Translate values from variable names used in model to ModelDailyState ===
  ! =============================================================================

  ModelDailyState(Gridiv,cMDS_porosity) = porosity(id) 
  ModelDailyState(Gridiv,cMDS_albDecTr) = albDecTr(id)   
  ModelDailyState(Gridiv,cMDS_albEveTr) = albEveTr(id)   
  ModelDailyState(Gridiv,cMDS_albGrass) = albGrass(id)   
  ModelDailyState(Gridiv,cMDS_DecidCap) = DecidCap(id) 
  ModelDailyState(Gridiv,cMDS_CumSnowfall) = CumSnowfall
       
  ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
  HDD_grids(:,:,Gridiv) = HDD(:,:) 
  GDD_grids(:,:,Gridiv) = GDD(:,:) 
  LAI_grids(:,:,Gridiv) = LAI(:,:) 
  WU_Day_grids(:,:,Gridiv) = WU_day(:,:)
  AlbDecTr_grids(:,Gridiv) = AlbDecTr(:)
  AlbEveTr_grids(:,Gridiv) = AlbEveTr(:)
  AlbGrass_grids(:,Gridiv) = AlbGrass(:)
  DecidCap_grids(:,Gridiv) = DecidCap(:) 
  Porosity_grids(:,Gridiv) = Porosity(:)    

  ! ---- Snow density of each surface
  ModelDailyState(Gridiv,cMDS_SnowDens(1:nsurf)) = SnowDens(1:nsurf)
  ModelDailyState(Gridiv,cMDS_SnowAlb) = SnowAlb

         
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
!=================================================================================== 