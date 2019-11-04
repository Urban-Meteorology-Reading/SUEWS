!SUEWS_Translate
!Translates - new input arrays (v2014b) to existing model variables
!           - between arrays for different grids and the model variables
!Made by HW&LJ Oct 2014
!-----------------------------------------------------------------------------------
! MH 21 Jun 2017 : Added anthropogenic CO2 charasteristic
! MH 16 Jun 2017 : Added biogenic CO2 characteristic
! HCW 13 Dec 2016 : LAIPower and LAIType for all vegetation types now used (previously only DecTr were used)
! HCW 12 Dec 2016 : Switched sign of lng so that input should be -ve for W, +ve for E, as is conventional
!Last modified HCW 26 Aug 2016
! NumCapita now uses average of day and night pop density, unless only one is specified
!Last modified HCW 06 Jul 2016
! Checks on ESTM fractions
!  - default setting to first ESTM Class code if surface not present and ESTM fractions do not sum to 1.
!Last modified HCW 29 Jun 2016
! Removed SoilMoistDay and StateDay
!Last modified: HCW 16 Jun 2016
! ESTM development for 7 surface types + snow, allowing 3x Paved classes and 5x Bldgs classes
! Currently surface characteristics are averaged here; probably want to average QS instead.
!Last modified: TS 13 Apr 2016
! Added AnOHM required variables.
!Last modified: LJ 06 Jul 2015
! Changed to read SnowAlb from ModelDailyState instead of SurfaceChar. Location also moved.
!Last modified: HCW 03 Jul 2015
! Use PopDensNighttime by default (not PopDensDaytime)
!Last modified: HCW 26 Jun 2015
! Translation of DailyState variables from the corresponding '_grids' arrays moved
! earlier in code in order to fix bug in DecidCap, AlbDec, Porosity.
!Last modified: HCW 28 Nov 2014
!
! To Do:
!       - Add AnOHM and ESTM info to FileChoices
!       - Check observed soil moisture works correctly!!
!       - Adjust model to allow water to runoff and sub-surface soil store for each surface type
!  - Adjust model to calculate LAI per surface
!  - Adjust model for SM per surface (measured characteristics)
!===================================================================================
SUBROUTINE SUEWS_Translate(Gridiv, ir, iMB)
   USE allocateArray
   USE ColNamesInputFiles
   USE ColNamesModelDailyState
   USE data_in
   USE defaultnotUsed, ONLY: notused, notusedi
   USE gis_data, ONLY: &
      areazh, vegfraction, veg_type, veg_fr, impervfraction, pervfraction, &
      nonwaterfraction, bldgh, evetreeh, dectreeh, treeh, faibldg, faievetree, faidectree, faitree
   USE mod_z, ONLY: z, z0m, z0m_in, zdm, zdm_in, zzd
   USE resist, ONLY: g1, g2, g3, g4, g5, g6, th, tl, s1, s2, kmax, gsmodel
   USE snowMod, ONLY: &
      snowlimpaved, SnowLimBldg, radmeltfact, tempmeltfact, &
      snowAlbMin, snowAlbMax, tau_a, tau_f, preciplimitalb, &
      snowdensmin, snowdensmax, tau_r, crwmin, crwmax, &
      preciplimit, SnowProf_24hr, snowalb, snowfallcum
   USE sues_data, ONLY: &
      surfacearea_ha, surfacearea, irrfracconif, irrfracdecid, irrfracgrass, &
      soildensity, soildepthmeas, smcap, soilrocks, pipecapacity, flowchange, &
      runofftowater, ie_start, ie_end, internalwateruse_h, faut, &
      ie_a, ie_m, daywat, daywatper, &
      wuareaevetr_m2, wuareadectr_m2, wuareagrass_m2, &
      AerodynamicResistanceMethod, RoughLenHeatMethod, StabilityMethod, tstep, tstep_prev

   USE time, ONLY: iy, id, it, imin, isec, dectime, dt_since_start
   USE ESTM_data
   USE WhereWhen, ONLY: gridid, gridid_text
   USE ESTM_module, ONLY: ESTM_translate

   IMPLICIT NONE

   INTEGER::Gridiv, &   !Index of the analysed grid (Gridcounter)
             ir, &       !Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
             iMB, &      !Chunk of met data
             id_prev

   INTEGER::iv, j, i
   !real (Kind(1d0)):: FCskip = -9   !NULL value used for output to FileChoices
   REAL(KIND(1d0)):: FCskip = -999  !NULL value used for output to FileChoices        (changed by HCW 24 May 2016)

   ! REAL(KIND(1d0)):: z0m_in, zdm_in  !Values of z0m and zdm provided in SiteSelect input file (do not get updated unlike z0d and z0m)

   CHARACTER(len=20):: grid_txt
   CHARACTER(len=4):: year_txt
   CHARACTER(len=12)::SsG_YYYY !Site, grid, year string

   CHARACTER(len=4):: iy_text
   CHARACTER(len=3):: id_text
   CHARACTER(len=2):: it_text, imin_text

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
   ! GridID = GridIDmatrix(Gridiv) ! also in SUEWS_Program so deleted here. TS 10 Jun 2018
   ! ---- Latitude and longitude
   lat = SurfaceChar(Gridiv, c_lat)
   lng = SurfaceChar(Gridiv, c_lng)
   ! ---- Timezone
   TIMEZONE = SurfaceChar(Gridiv, c_tz)
   ! ---- Altitude [m]
   Alt = SurfaceChar(Gridiv, c_Alt)
   ! ---- Measurement height [m]
   z = SurfaceChar(Gridiv, c_z)
   ! ---- Surface area [ha]
   SurfaceArea_ha = SurfaceChar(Gridiv, c_Area)
   ! Change from ha to m2 (was in RunControlByGridByYear)
   SurfaceArea = SurfaceArea_ha*10000 !Change surface area from ha to m^2

   ! ---- Surface fractions (previously in LUMPS_gis_read)
   sfr(PavSurf) = SurfaceChar(Gridiv, c_FrPaved)  ! Paved
   sfr(BldgSurf) = SurfaceChar(Gridiv, c_FrBldgs)  ! Bldgs
   sfr(ConifSurf) = SurfaceChar(Gridiv, c_FrEveTr)  ! Everg
   sfr(DecidSurf) = SurfaceChar(Gridiv, c_FrDecTr)  ! Decid
   sfr(GrassSurf) = SurfaceChar(Gridiv, c_FrGrass)  ! Grass
   sfr(BSoilSurf) = SurfaceChar(Gridiv, c_FrBSoil)  ! BSoil
   sfr(WaterSurf) = SurfaceChar(Gridiv, c_FrWater)  ! Water

   ! Check the surface fractions add up to 1 (or close to 1)
   IF (SUM(sfr) > 1.001 .OR. SUM(sfr) < 0.999) &
      CALL ErrorHint(10, 'Surface fractions (Fr_) should add up to 1.', SUM(sfr), notUsed, notUsedI)

   ! ---- Irrigated fractions
   IrrFracConif = SurfaceChar(Gridiv, c_IrrEveTrFrac)  ! Everg
   IrrFracDecid = SurfaceChar(Gridiv, c_IrrDecTrFrac)  ! Decid
   IrrFracGrass = SurfaceChar(Gridiv, c_IrrGrassFrac)  ! Grass

   ! ---------------------------------------------------------------------------------
   ! --------- Surface cover calculations (previously in LUMPS_gis_read) -------------

   ! ---- Buildings and trees fraction ----
   areaZh = (sfr(BldgSurf) + sfr(ConifSurf) + sfr(DecidSurf))

   ! ---- Vegetated fraction ----
   VegFraction = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf))
   !VegFraction = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))

   ! ---- Vegetated fraction (for LUMPS) ----
   ! For LUMPS, vegetated fraction includes Water and Bare soil surfaces
   IF (veg_type == 1) THEN          ! area vegetated
      veg_fr = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf) + sfr(WaterSurf))
   ELSEIF (veg_type == 2) THEN      ! area irrigated
      veg_fr = (IrrFracConif*sfr(ConifSurf) + IrrFracDecid*sfr(DecidSurf) + IrrFracGrass*sfr(GrassSurf))
   ENDIF

   ImpervFraction = (sfr(PavSurf) + sfr(BldgSurf))
   PervFraction = (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf) + sfr(WaterSurf))
   NonWaterFraction = (sfr(PavSurf) + sfr(BldgSurf) + sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BSoilSurf))
   ! ---------------------------------------------------------------------------------

   ! ---- Heights & frontal areas
   BldgH = SurfaceChar(Gridiv, c_HBldgs)                                                            ! Building height [m]
   EveTreeH = SurfaceChar(Gridiv, c_HEveTr)                                                            ! Evergreen tree height [m]
   DecTreeH = SurfaceChar(Gridiv, c_HDecTr)                                                            ! Deciduous tree height [m]
   IF (sfr(ConifSurf) + sfr(DecidSurf) > 0.) THEN ! avoid arithmetic error
      TreeH = (EveTreeH*sfr(ConifSurf) + DecTreeH*sfr(DecidSurf))/(sfr(ConifSurf) + sfr(DecidSurf))     ! Average tree height [m]
   ELSE
      TreeH = 1.
   END IF

   FAIBldg = SurfaceChar(Gridiv, c_FAIBldgs)                                                          ! Frontal area index for buildings
   FAIEveTree = SurfaceChar(Gridiv, c_FAIEveTr)                                                          ! Frontal area index for evergreen trees
   FAIDecTree = SurfaceChar(Gridiv, c_FAIDecTr)                                                          ! Frontal area index for deciduous trees
   IF (sfr(ConifSurf) + sfr(DecidSurf) > 0.) THEN ! avoid arithmetic error
      FAITree = (FAIEveTree*sfr(ConifSurf) + FAIDecTree*sfr(DecidSurf))/(sfr(ConifSurf) + sfr(DecidSurf)) ! Frontal area index for trees
   ELSE
      FAITree = 1.
   END IF

   z0m = SurfaceChar(Gridiv, c_z0m)                                                               ! Roughness length [m]
   zdm = SurfaceChar(Gridiv, c_zdm)                                                               ! Displacement height [m]
   ! z0m and zdm can vary in time depending on z0method selected. Save the input values here
   z0m_in = z0m
   zdm_in = zdm

   ! ---- Population density [ha-1]
   ! Weekend fraction added to daytime population density
   PopDensDaytime = SurfaceChar(Gridiv, (/c_PopDensDay, c_PopDensDay/))   ! Daytime population density [ha-1]
   PopDensNighttime = SurfaceChar(Gridiv, c_PopDensNight) ! Night-time population density [ha-1]
   IF (PopDensDaytime(1) >= 0 .AND. PopDensNighttime < 0) PopDensNighttime = PopDensDaytime(1)  !If only daytime data provided, use them
   IF (PopDensDaytime(1) < 0 .AND. PopDensNighttime >= 0) PopDensDaytime(1) = PopDensNighttime  !If only night-time data provided, use them
   PopDensDaytime(2) = PopDensNighttime + (PopDensDaytime(1) - PopDensNighttime)*SurfaceChar(Gridiv, c_FrPDDwe) !Use weekend fraction to daytime population
   ! the following part has been moved into `SUEWS_cal_Main` as `NumCapita` can be derived there
   ! IF (PopDensDaytime(1) >= 0 .AND. PopDensNighttime >= 0) NumCapita(1) = (PopDensDaytime(1) + PopDensNighttime)/2  !If both, use average
   ! IF (PopDensDaytime(2) >= 0 .AND. PopDensNighttime >= 0) NumCapita(2) = (PopDensDaytime(2) + PopDensNighttime)/2  !If both, use average

   ! ! IF (PopDensDaytime >= 0 .AND. PopDensNighttime >= 0) NumCapita = (PopDensDaytime + PopDensNighttime)/2  !If both, use average ! moved to `AnthropogenicEmissions`, TS 27 Dec 2018

   ! ---- Traffic rate
   TrafficRate = SurfaceChar(Gridiv, (/c_TrafficRate_WD, c_TrafficRate_WE/)) ! Mean traffic rate within modelled area
   ! ---- Building energy use
   QF0_BEU = SurfaceChar(Gridiv, (/c_QF0_BEU_WD, c_QF0_BEU_WE/))   ! Building energy use within modelled area

   ! ---- Albedo [-]
   alb(1:nsurf) = SurfaceChar(Gridiv, c_AlbMax)   !Use maximum albedos as default value (AlbMin for veg surfaces handled below)

   ! ---- Set min & max albedo for vegetated surfaces (min albedo not currently used for NonVeg or Water surfaces)
   AlbMin_EveTr = SurfaceChar(Gridiv, c_AlbMin(ConifSurf))
   AlbMax_EveTr = SurfaceChar(Gridiv, c_AlbMax(ConifSurf))
   AlbMin_DecTr = SurfaceChar(Gridiv, c_AlbMin(DecidSurf))
   AlbMax_DecTr = SurfaceChar(Gridiv, c_AlbMax(DecidSurf))
   AlbMin_Grass = SurfaceChar(Gridiv, c_AlbMin(GrassSurf))
   AlbMax_Grass = SurfaceChar(Gridiv, c_AlbMax(GrassSurf))

   ! ---- Emissivity [-]
   emis(1:nsurf) = SurfaceChar(Gridiv, c_Emis)
   emis_snow = SurfaceChar(Gridiv, c_SnowEmis)

   ! ---- Storage capacities [mm]
   StoreDrainPrm(1, 1:nsurf) = SurfaceChar(Gridiv, c_StorMin)   ! Minimum
   StoreDrainPrm(5, 1:nsurf) = SurfaceChar(Gridiv, c_StorMax)   ! Maximum
   StoreDrainPrm(6, 1:nsurf) = StoreDrainPrm(1, 1:nsurf)  !Set storage capacities for all surface to minimum (DecTr changes with time in Calculations).

   ! ---- Set min & max storage capacities for DecTr
   CapMin_dec = StoreDrainPrm(1, DecidSurf)
   CapMax_dec = StoreDrainPrm(5, DecidSurf)
   ! ---- Set min & max porosity for DecTr
   PorMin_dec = SurfaceChar(Gridiv, c_PorosityMin(ivDecid))   ! Minimum
   PorMax_dec = SurfaceChar(Gridiv, c_PorosityMax(ivDecid))   ! Minimum

   ! ---- Threshold for wet evaporation [mm]
   WetThresh(1:nsurf) = SurfaceChar(Gridiv, c_WetThresh)

   ! ---- Limit for state [mm]
   StateLimit(1:nsurf) = SurfaceChar(Gridiv, c_StateLimit)

   ! ---- Water depth [mm]
   WaterDepth = SurfaceChar(Gridiv, c_WaterDepth)

   ! ---- Drainage
   StoreDrainPrm(2, 1:nsurf) = SurfaceChar(Gridiv, c_DrEq)      ! Drainage equation
   StoreDrainPrm(3, 1:nsurf) = SurfaceChar(Gridiv, c_DrCoef1)   ! Drainage coef 1
   StoreDrainPrm(4, 1:nsurf) = SurfaceChar(Gridiv, c_DrCoef2)   ! Drainage coef 2

   ! ---- Limit of SWE (each surface except Water)
   SnowPackLimit(1:(nsurf - 1)) = SurfaceChar(Gridiv, c_SnowLimPat(1:(nsurf - 1)))

   ! ---- Snow limit for removal (only impervious surfaces)
   SnowLimPaved = SurfaceChar(Gridiv, c_SnowLimRem(PavSurf))
   SnowLimBldg = SurfaceChar(Gridiv, c_SnowLimRem(BldgSurf))
   !SnowLimBSoil = SurfaceChar(Gridiv,c_SnowLimRem(BSoilSurf))   !Snow clearing not applicable to bare soil surface

   ! ---- Soil characteristics (each surface except Water)
   SoilDepth(1:(nsurf - 1)) = SurfaceChar(Gridiv, c_SoilDepth(1:(nsurf - 1))) ! Depth of sub-surface soil store [mm]
   SoilStoreCap(1:(nsurf - 1)) = SurfaceChar(Gridiv, c_SoilStCap(1:(nsurf - 1))) ! Soil store capacity [mm]
   SatHydraulicConduct(1:(nsurf - 1)) = SurfaceChar(Gridiv, c_KSat(1:(nsurf - 1)))      ! Hydraulic conductivity of saturated soil [mm s-1]
   !SoilDensity(1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilDens(1:(nsurf-1))) ! Soil density [kg m-3]
   ! Not yet implemented in model
   !InfiltrationRate  (1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilInfRate(1:(nsurf-1))) ! Infiltration rate [mm h-1]

   !! Observed soil characteristics
   !SoilDensity  (1:(nsurf-1)) = SurfaceChar(Gridiv,c_SoilDens(1:(nsurf-1))) ! Soil density [kg m-3]
   !SoilDepthMeas(1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSMDepth(1:(nsurf-1)))
   !SmCap        (1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSMMax(1:(nsurf-1)))
   !SoilRocks    (1:(nsurf-1)) = SurfaceChar(Gridiv,c_ObsSNRFrac(1:(nsurf-1)))
   !!Obs soil characteristics now in SUEWS_Soil, i.e. per surface; single value was given previously in FunctionalTypes
   !!Take first row here for testing !! Need to alter model later...
   SoilDensity = SurfaceChar(Gridiv, c_SoilDens(1)) !!Not sure this works correctly - need to check
   SoilDepthMeas = SurfaceChar(Gridiv, c_ObsSMDepth(1))
   SmCap = SurfaceChar(Gridiv, c_ObsSMMax(1))
   SoilRocks = SurfaceChar(Gridiv, c_ObsSNRFrac(1))

   ! ---- Vegetation characteristics (pervious surfaces)
   BaseT(1:nvegsurf) = SurfaceChar(Gridiv, c_BaseT)
   BaseTe(1:nvegsurf) = SurfaceChar(Gridiv, c_BaseTe)
   GDDFull(1:nvegsurf) = SurfaceChar(Gridiv, c_GDDFull)
   SDDFull(1:nvegsurf) = SurfaceChar(Gridiv, c_SDDFull)
   LAIMin(1:nvegsurf) = SurfaceChar(Gridiv, c_LAIMin)
   LAIMax(1:nvegsurf) = SurfaceChar(Gridiv, c_LAIMax)
   MaxConductance(1:nvegsurf) = SurfaceChar(Gridiv, c_GsMax)

   alpha_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_alpha_bioCO2)
   beta_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_beta_bioCO2)
   theta_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_theta_bioCO2)
   alpha_enh_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_alpha_enh_bioCO2)
   beta_enh_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_beta_enh_bioCO2)
   resp_a(1:nvegsurf) = SurfaceChar(Gridiv, c_resp_a)
   resp_b(1:nvegsurf) = SurfaceChar(Gridiv, c_resp_b)
   min_res_bioCO2(1:nvegsurf) = SurfaceChar(Gridiv, c_min_res_bioCO2)

   ! ---- LAI characteristics (updated HCW 13 Dec 2016)
   LAItype(1:nvegsurf) = INT(SurfaceChar(Gridiv, c_LAIEq(1:nvegsurf)))
   LAIPower(1, 1:nvegsurf) = SurfaceChar(Gridiv, c_LeafGP1(1:nvegsurf))
   LAIPower(2, 1:nvegsurf) = SurfaceChar(Gridiv, c_LeafGP2(1:nvegsurf))
   LAIPower(3, 1:nvegsurf) = SurfaceChar(Gridiv, c_LeafOP1(1:nvegsurf))
   LAIPower(4, 1:nvegsurf) = SurfaceChar(Gridiv, c_LeafOP2(1:nvegsurf))

   ! ---- LUMPS-related parameters
   DRAINRT = SurfaceChar(Gridiv, c_LUMPSDr)      ! LUMPS Drainage rate [mm h-1]
   RAINCOVER = SurfaceChar(Gridiv, c_LUMPSCover)   ! LUMPS Limit when surface totally wet [mm]
   RAINMAXRES = SurfaceChar(Gridiv, c_LUMPSMaxRes)  ! LUMPS Maximum water bucket reservoir [mm]

   ! ---- NARP-related parameters
   TRANS_SITE = SurfaceChar(Gridiv, c_NARPTrans)    ! NARP atmospheric transmissivity

   ! ---- Snow-related characteristics
   RadMeltFact = SurfaceChar(Gridiv, c_SnowRMFactor)
   TempMeltFact = SurfaceChar(Gridiv, c_SnowTMFactor)
   SnowAlbMin = SurfaceChar(Gridiv, c_SnowAlbMin)
   SnowAlbMax = SurfaceChar(Gridiv, c_SnowAlbMax)
   tau_a = SurfaceChar(Gridiv, c_Snowtau_a)
   tau_f = SurfaceChar(Gridiv, c_Snowtau_f)
   PrecipLimitAlb = SurfaceChar(Gridiv, c_SnowPlimAlb)
   SnowDensMin = SurfaceChar(Gridiv, c_SnowSDMin)
   SnowDensMax = SurfaceChar(Gridiv, c_SnowSDMax)
   tau_r = SurfaceChar(Gridiv, c_Snowtau_r)
   CRWMin = SurfaceChar(Gridiv, c_SnowCRWMin)
   CRWMax = SurfaceChar(Gridiv, c_SnowCRWMax)
   PrecipLimit = SurfaceChar(Gridiv, c_SnowPLimSnow)

   ! ---- Conductance parameters
   G1 = SurfaceChar(Gridiv, c_GsG1)
   G2 = SurfaceChar(Gridiv, c_GsG2)
   G3 = SurfaceChar(Gridiv, c_GsG3)
   G4 = SurfaceChar(Gridiv, c_GsG4)
   G5 = SurfaceChar(Gridiv, c_GsG5)
   G6 = SurfaceChar(Gridiv, c_GsG6)
   TH = SurfaceChar(Gridiv, c_GsTH)
   TL = SurfaceChar(Gridiv, c_GsTL)
   S1 = SurfaceChar(Gridiv, c_GsS1)
   S2 = SurfaceChar(Gridiv, c_GsS2)
   Kmax = SurfaceChar(Gridiv, c_GsKmax)
   gsModel = INT(SurfaceChar(Gridiv, c_gsModel))

   ! ---- Pipe capacity (was from SiteSpecificParam.txt)
   PipeCapacity = SurfaceChar(Gridiv, c_PipeCapacity)

   ! ---- Water flows (was from SiteSpecificParam.txt)
   FlowChange = SurfaceChar(Gridiv, c_FlowChange)
   RunoffToWater = SurfaceChar(Gridiv, c_RunoffToWater)

   ! ---- Daylight saving (was from ModelledYears.txt)
   startDLS = INT(SurfaceChar(Gridiv, c_StartDLS))
   endDLS = INT(SurfaceChar(Gridiv, c_EndDLS))

   ! ---- OHM coeffs (was in SUEWS_OHMnew.f95, subroutine OHMinitialize)
   OHM_coef = 0 ! Initialise OHM_coef
   ! Surface types in OHM_coef: Paved, Roof, Conif, Decid, Grass, BareSoil, Water, CANYON, Snow
   ! No canyon in SurfaceChar, so
   !  transfer coeffs for surface types 1-7,
   !  then skip row in OHM_coef (canyon),
   !  then transfer coeffs for snow surface (8th surface in SurfaceChar; 9th surface in OHM_Coefs)
   ! Summer wet
   OHM_coef(1:nsurf, 1, 1) = SurfaceChar(Gridiv, c_a1_SWet(1:nsurf)) !1:nsurf a1 Summer wet
   OHM_coef(nsurf + 1, 1, 1) = SurfaceChar(Gridiv, c_a1_SWet(nsurf + 1)) !Snow    a1 Summer wet
   OHM_coef(1:nsurf, 1, 2) = SurfaceChar(Gridiv, c_a2_SWet(1:nsurf)) !1:nsurf a2 Summer wet
   OHM_coef(nsurf + 1, 1, 2) = SurfaceChar(Gridiv, c_a2_SWet(nsurf + 1)) !Snow    a2 Summer wet
   OHM_coef(1:nsurf, 1, 3) = SurfaceChar(Gridiv, c_a3_SWet(1:nsurf)) !1:nsurf a3 Summer wet
   OHM_coef(nsurf + 1, 1, 3) = SurfaceChar(Gridiv, c_a3_SWet(nsurf + 1)) !Snow    a3 Summer wet
   ! Summer dry
   OHM_coef(1:nsurf, 2, 1) = SurfaceChar(Gridiv, c_a1_SDry(1:nsurf)) !1:nsurf a1 Summer dry
   OHM_coef(nsurf + 1, 2, 1) = SurfaceChar(Gridiv, c_a1_SDry(nsurf + 1)) !Snow    a1 Summer dry
   OHM_coef(1:nsurf, 2, 2) = SurfaceChar(Gridiv, c_a2_SDry(1:nsurf)) !1:nsurf a2 Summer dry
   OHM_coef(nsurf + 1, 2, 2) = SurfaceChar(Gridiv, c_a2_SDry(nsurf + 1)) !Snow    a2 Summer dry
   OHM_coef(1:nsurf, 2, 3) = SurfaceChar(Gridiv, c_a3_SDry(1:nsurf)) !1:nsurf a3 Summer dry
   OHM_coef(nsurf + 1, 2, 3) = SurfaceChar(Gridiv, c_a3_SDry(nsurf + 1)) !Snow    a3 Summer dry
   ! Winter wet
   OHM_coef(1:nsurf, 3, 1) = SurfaceChar(Gridiv, c_a1_WWet(1:nsurf)) !1:nsurf a1 Winter wet
   OHM_coef(nsurf + 1, 3, 1) = SurfaceChar(Gridiv, c_a1_WWet(nsurf + 1)) !Snow    a1 Winter wet
   OHM_coef(1:nsurf, 3, 2) = SurfaceChar(Gridiv, c_a2_WWet(1:nsurf)) !1:nsurf a2 Winter wet
   OHM_coef(nsurf + 1, 3, 2) = SurfaceChar(Gridiv, c_a2_WWet(nsurf + 1)) !Snow    a2 Winter wet
   OHM_coef(1:nsurf, 3, 3) = SurfaceChar(Gridiv, c_a3_WWet(1:nsurf)) !1:nsurf a3 Winter wet
   OHM_coef(nsurf + 1, 3, 3) = SurfaceChar(Gridiv, c_a3_WWet(nsurf + 1)) !Snow    a3 Winter wet
   ! Winter dry
   OHM_coef(1:nsurf, 4, 1) = SurfaceChar(Gridiv, c_a1_WDry(1:nsurf)) !1:nsurf a1 Winter dry
   OHM_coef(nsurf + 1, 4, 1) = SurfaceChar(Gridiv, c_a1_WDry(nsurf + 1)) !Snow    a1 Winter dry
   OHM_coef(1:nsurf, 4, 2) = SurfaceChar(Gridiv, c_a2_WDry(1:nsurf)) !1:nsurf a2 Winter dry
   OHM_coef(nsurf + 1, 4, 2) = SurfaceChar(Gridiv, c_a2_WDry(nsurf + 1)) !Snow    a2 Winter dry
   OHM_coef(1:nsurf, 4, 3) = SurfaceChar(Gridiv, c_a3_WDry(1:nsurf)) !1:nsurf a3 Winter dry
   OHM_coef(nsurf + 1, 4, 3) = SurfaceChar(Gridiv, c_a3_WDry(nsurf + 1)) !Snow    a3 Winter dry
   ! OHM thresholds
   OHM_threshSW(1:nsurf) = SurfaceChar(Gridiv, c_OHMThresh_SW(1:nsurf)) !1:nsurf
   OHM_threshSW(nsurf + 1) = SurfaceChar(Gridiv, c_OHMThresh_SW(nsurf + 1)) !Snow
   OHM_threshWD(1:nsurf) = SurfaceChar(Gridiv, c_OHMThresh_WD(1:nsurf)) !1:nsurf
   OHM_threshWD(nsurf + 1) = SurfaceChar(Gridiv, c_OHMThresh_WD(nsurf + 1)) !Snow

   ! ---- ESTM characteristics -------------------------
   ! HCW 16 Jun 2016
   ! Wall fraction for ESTM (in SiteSelect.txt)
   IF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
      AreaWall = SurfaceChar(Gridiv, c_AreaWall)
      fwall = AreaWall/SurfaceArea

      ! Get surface fractions for ESTM classes for Bldgs and Paved surfaces
      ESTMsfr_Paved = SurfaceChar(Gridiv, c_Fr_ESTMClass_Paved)   !Dim 3
      ESTMsfr_Bldgs = SurfaceChar(Gridiv, c_Fr_ESTMClass_Bldgs)   !Dim 5
      !Check these sum to 1 and are consistent with sfr of Paved and Bldgs surface types
      IF (sfr(PavSurf) > 0) THEN  !If surface exists, ESTM fractions must be correct
         IF (SUM(ESTMsfr_Paved) > 1.001 .OR. SUM(ESTMsfr_Paved) < 0.999) THEN
            CALL ErrorHint(10, 'Surface fractions (Fr_ESTMClass_Paved) should sum to 1.', SUM(ESTMsfr_Paved), notUsed, notUsedI)
         ENDIF
      ELSEIF (sfr(PavSurf) == 0) THEN !If surface does not exist, ESTM fraction does not matter
         IF (SUM(ESTMsfr_Paved) > 1.001 .OR. SUM(ESTMsfr_Paved) < 0.999) THEN   !If ESTM fractions do not sum to 1, set here
            ESTMsfr_Paved(1) = 1.000
            ESTMsfr_Paved(2:3) = 0.000
            CALL ErrorHint(67, 'ESTM Paved classes do not sum to 1 (but no Paved surface present).', &
                           SUM(ESTMsfr_Paved), notUsed, notUsedI)
         ENDIF
      ENDIF
      IF (sfr(BldgSurf) > 0) THEN
         IF (SUM(ESTMsfr_Bldgs) > 1.001 .OR. SUM(ESTMsfr_Bldgs) < 0.999) THEN
            CALL ErrorHint(10, 'Surface fractions (Fr_ESTMClass_Bldgs) should sum to 1.', SUM(ESTMsfr_Bldgs), notUsed, notUsedI)
         ENDIF
      ELSEIF (sfr(BldgSurf) == 0) THEN !If surface does not exist, ESTM fraction does not matter
         IF (SUM(ESTMsfr_Bldgs) > 1.001 .OR. SUM(ESTMsfr_Bldgs) < 0.999) THEN   !If ESTM fractions do not sum to 1, set here
            ESTMsfr_Bldgs(1) = 1.000
            ESTMsfr_Bldgs(2:5) = 0.000
            CALL ErrorHint(67, 'ESTM Bldgs classes do not sum to 1 (but no Bldgs surface present).', &
                           SUM(ESTMsfr_Bldgs), notUsed, notUsedI)
         ENDIF
      ENDIF

      ! ===== PAVED =====
      ! First combine characteristics of the 3x Paved classes
      IF (SurfaceChar(Gridiv, c_ESTMCode(PavSurf)) == 0) THEN   ! If Code = 0, use multiple classes
         ! Get characteristics of each Paved class
         DO i = 1, 3
            zSurf_Paved(:, i) = SurfaceChar(Gridiv, (/c_Surf_thick1_Paved(i), c_Surf_thick2_Paved(i), c_Surf_thick3_Paved(i), &
                                                      c_Surf_thick4_Paved(i), c_Surf_thick5_Paved(i)/))
            kSurf_Paved(:, i) = SurfaceChar(Gridiv, (/c_Surf_k1_Paved(i), c_Surf_k2_Paved(i), c_Surf_k3_Paved(i), &
                                                      c_Surf_k4_Paved(i), c_Surf_k5_Paved(i)/))
            rSurf_Paved(:, i) = SurfaceChar(Gridiv, (/c_Surf_rhoCp1_Paved(i), c_Surf_rhoCp2_Paved(i), c_Surf_rhoCp3_Paved(i), &
                                                      c_Surf_rhoCp4_Paved(i), c_Surf_rhoCp5_Paved(i)/))
         ENDDO
         ! Average characteristics of each Paved class according to surface fractions (these sum to 1)
         zSurf_SUEWSsurfs(:, PavSurf) = zSurf_Paved(:, 1)*ESTMsfr_Paved(1) &
                                        + zSurf_Paved(:, 2)*ESTMsfr_Paved(2) &
                                        + zSurf_Paved(:, 3)*ESTMsfr_Paved(3)
         kSurf_SUEWSsurfs(:, PavSurf) = kSurf_Paved(:, 1)*ESTMsfr_Paved(1) &
                                        + kSurf_Paved(:, 2)*ESTMsfr_Paved(2) &
                                        + kSurf_Paved(:, 3)*ESTMsfr_Paved(3)
         rSurf_SUEWSsurfs(:, PavSurf) = rSurf_Paved(:, 1)*ESTMsfr_Paved(1) &
                                        + rSurf_Paved(:, 2)*ESTMsfr_Paved(2) &
                                        + rSurf_Paved(:, 3)*ESTMsfr_Paved(3)
      ELSEIF (SurfaceChar(Gridiv, c_ESTMCode(PavSurf)) /= 0) THEN   !Otherwise use single values
         zSurf_SUEWSsurfs(:, PavSurf) = SurfaceChar(Gridiv, &
                                                    (/c_Surf_thick1(PavSurf), c_Surf_thick2(PavSurf), c_Surf_thick3(PavSurf), &
                                                      c_Surf_thick4(PavSurf), c_Surf_thick5(PavSurf)/))
         kSurf_SUEWSsurfs(:, PavSurf) = SurfaceChar(Gridiv, &
                                                    (/c_Surf_k1(PavSurf), c_Surf_k2(PavSurf), c_Surf_k3(PavSurf), &
                                                      c_Surf_k4(PavSurf), c_Surf_k5(PavSurf)/))
         rSurf_SUEWSsurfs(:, PavSurf) = SurfaceChar(Gridiv, &
                                                    (/c_Surf_rhoCp1(PavSurf), c_Surf_rhoCp2(PavSurf), c_Surf_rhoCp3(PavSurf), &
                                                      c_Surf_rhoCp4(PavSurf), c_Surf_rhoCp5(PavSurf)/))
      ENDIF

      ! ===== BLDGS =====
      ! Combine characteristics of 5x Bldgs classes into one
      IF (SurfaceChar(Gridiv, c_ESTMCode(BldgSurf)) == 0) THEN   ! If Code = 0, use multiple classes
         ! Get characteristics of each Bldgs class
         DO i = 1, 5
            zSurf_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Surf_thick1_Bldgs(i), c_Surf_thick2_Bldgs(i), c_Surf_thick3_Bldgs(i), &
                                                      c_Surf_thick4_Bldgs(i), c_Surf_thick5_Bldgs(i)/))
            kSurf_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Surf_k1_Bldgs(i), c_Surf_k2_Bldgs(i), c_Surf_k3_Bldgs(i), &
                                                      c_Surf_k4_Bldgs(i), c_Surf_k5_Bldgs(i)/))
            rSurf_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Surf_rhoCp1_Bldgs(i), c_Surf_rhoCp2_Bldgs(i), c_Surf_rhoCp3_Bldgs(i), &
                                                      c_Surf_rhoCp4_Bldgs(i), c_Surf_rhoCp5_Bldgs(i)/))
            zwall_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Wall_thick1_Bldgs(i), c_Wall_thick2_Bldgs(i), c_Wall_thick3_Bldgs(i), &
                                                      c_Wall_thick4_Bldgs(i), c_Wall_thick5_Bldgs(i)/))
            kwall_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Wall_k1_Bldgs(i), c_Wall_k2_Bldgs(i), c_Wall_k3_Bldgs(i), &
                                                      c_Wall_k4_Bldgs(i), c_Wall_k5_Bldgs(i)/))
            rwall_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Wall_rhoCp1_Bldgs(i), c_Wall_rhoCp2_Bldgs(i), c_Wall_rhoCp3_Bldgs(i), &
                                                      c_Wall_rhoCp4_Bldgs(i), c_Wall_rhoCp5_Bldgs(i)/))
            zibld_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Internal_thick1_Bldgs(i), c_Internal_thick2_Bldgs(i), &
                                                      c_Internal_thick3_Bldgs(i), &
                                                      c_Internal_thick4_Bldgs(i), c_Internal_thick5_Bldgs(i)/))
            kibld_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Internal_k1_Bldgs(i), c_Internal_k2_Bldgs(i), c_Internal_k3_Bldgs(i), &
                                                      c_Internal_k4_Bldgs(i), c_Internal_k5_Bldgs(i)/))
            ribld_Bldgs(:, i) = SurfaceChar(Gridiv, (/c_Internal_rhoCp1_Bldgs(i), c_Internal_rhoCp2_Bldgs(i), &
                                                      c_Internal_rhoCp3_Bldgs(i), &
                                                      c_Internal_rhoCp4_Bldgs(i), c_Internal_rhoCp5_Bldgs(i)/))
            nroom_Bldgs(i) = SurfaceChar(Gridiv, c_nroom_Bldgs(i))
            alb_ibld_Bldgs(i) = SurfaceChar(Gridiv, c_alb_ibld_Bldgs(i))
            em_ibld_Bldgs(i) = SurfaceChar(Gridiv, c_em_ibld_Bldgs(i))
            CH_iwall_Bldgs(i) = SurfaceChar(Gridiv, c_CH_iwall_Bldgs(i))
            CH_iroof_Bldgs(i) = SurfaceChar(Gridiv, c_CH_iroof_Bldgs(i))
            CH_ibld_Bldgs(i) = SurfaceChar(Gridiv, c_CH_ibld_Bldgs(i))
         ENDDO
         ! Average characteristics of each Bldgs class according to surface fractions (these sum to 1)
         zSurf_SUEWSsurfs(:, BldgSurf) = zSurf_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                                         + zSurf_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                                         + zSurf_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                                         + zSurf_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                                         + zSurf_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         kSurf_SUEWSsurfs(:, BldgSurf) = kSurf_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                                         + kSurf_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                                         + kSurf_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                                         + kSurf_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                                         + kSurf_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         rSurf_SUEWSsurfs(:, BldgSurf) = rSurf_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                                         + rSurf_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                                         + rSurf_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                                         + rSurf_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                                         + rSurf_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         !Wall
         zwall = zwall_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + zwall_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + zwall_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + zwall_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + zwall_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         kwall = kwall_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + kwall_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + kwall_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + kwall_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + kwall_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         rwall = rwall_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + rwall_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + rwall_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + rwall_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + rwall_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         !Internal
         zibld = zibld_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + zibld_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + zibld_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + zibld_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + zibld_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         kibld = kibld_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + kibld_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + kibld_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + kibld_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + kibld_Bldgs(:, 5)*ESTMsfr_Bldgs(5)
         ribld = ribld_Bldgs(:, 1)*ESTMsfr_Bldgs(1) &
                 + ribld_Bldgs(:, 2)*ESTMsfr_Bldgs(2) &
                 + ribld_Bldgs(:, 3)*ESTMsfr_Bldgs(3) &
                 + ribld_Bldgs(:, 4)*ESTMsfr_Bldgs(4) &
                 + ribld_Bldgs(:, 5)*ESTMsfr_Bldgs(5)

         nroom = nroom_Bldgs(1)*ESTMsfr_Bldgs(1) &
                 + nroom_Bldgs(2)*ESTMsfr_Bldgs(2) &
                 + nroom_Bldgs(3)*ESTMsfr_Bldgs(3) &
                 + nroom_Bldgs(4)*ESTMsfr_Bldgs(4) &
                 + nroom_Bldgs(5)*ESTMsfr_Bldgs(5)
         alb_ibld = alb_ibld_Bldgs(1)*ESTMsfr_Bldgs(1) &
                    + alb_ibld_Bldgs(2)*ESTMsfr_Bldgs(2) &
                    + alb_ibld_Bldgs(3)*ESTMsfr_Bldgs(3) &
                    + alb_ibld_Bldgs(4)*ESTMsfr_Bldgs(4) &
                    + alb_ibld_Bldgs(5)*ESTMsfr_Bldgs(5)
         em_ibld = em_ibld_Bldgs(1)*ESTMsfr_Bldgs(1) &
                   + em_ibld_Bldgs(2)*ESTMsfr_Bldgs(2) &
                   + em_ibld_Bldgs(3)*ESTMsfr_Bldgs(3) &
                   + em_ibld_Bldgs(4)*ESTMsfr_Bldgs(4) &
                   + em_ibld_Bldgs(5)*ESTMsfr_Bldgs(5)
         CH_iwall = CH_iwall_Bldgs(1)*ESTMsfr_Bldgs(1) &
                    + CH_iwall_Bldgs(2)*ESTMsfr_Bldgs(2) &
                    + CH_iwall_Bldgs(3)*ESTMsfr_Bldgs(3) &
                    + CH_iwall_Bldgs(4)*ESTMsfr_Bldgs(4) &
                    + CH_iwall_Bldgs(5)*ESTMsfr_Bldgs(5)
         CH_iroof = CH_iroof_Bldgs(1)*ESTMsfr_Bldgs(1) &
                    + CH_iroof_Bldgs(2)*ESTMsfr_Bldgs(2) &
                    + CH_iroof_Bldgs(3)*ESTMsfr_Bldgs(3) &
                    + CH_iroof_Bldgs(4)*ESTMsfr_Bldgs(4) &
                    + CH_iroof_Bldgs(5)*ESTMsfr_Bldgs(5)
         CH_ibld = CH_ibld_Bldgs(1)*ESTMsfr_Bldgs(1) &
                   + CH_ibld_Bldgs(2)*ESTMsfr_Bldgs(2) &
                   + CH_ibld_Bldgs(3)*ESTMsfr_Bldgs(3) &
                   + CH_ibld_Bldgs(4)*ESTMsfr_Bldgs(4) &
                   + CH_ibld_Bldgs(5)*ESTMsfr_Bldgs(5)

      ELSEIF (SurfaceChar(Gridiv, c_ESTMCode(BldgSurf)) /= 0) THEN   !Otherwise use single values
         zSurf_SUEWSsurfs(:, BldgSurf) = SurfaceChar(Gridiv, (/c_Surf_thick1(BldgSurf), c_Surf_thick2(BldgSurf), &
                                                               c_Surf_thick3(BldgSurf), &
                                                               c_Surf_thick4(BldgSurf), c_Surf_thick5(BldgSurf)/))
         kSurf_SUEWSsurfs(:, BldgSurf) = SurfaceChar(Gridiv, (/c_Surf_k1(BldgSurf), c_Surf_k2(BldgSurf), c_Surf_k3(BldgSurf), &
                                                               c_Surf_k4(BldgSurf), c_Surf_k5(BldgSurf)/))
         rSurf_SUEWSsurfs(:, BldgSurf) = SurfaceChar(Gridiv, (/c_Surf_rhoCp1(BldgSurf), c_Surf_rhoCp2(BldgSurf), &
                                                               c_Surf_rhoCp3(BldgSurf), &
                                                               c_Surf_rhoCp4(BldgSurf), c_Surf_rhoCp5(BldgSurf)/))
         zwall = SurfaceChar(Gridiv, (/c_Wall_thick1, c_Wall_thick2, c_Wall_thick3, c_Wall_thick4, c_Wall_thick5/))
         kwall = SurfaceChar(Gridiv, (/c_Wall_k1, c_Wall_k2, c_Wall_k3, c_Wall_k4, c_Wall_k5/))
         rwall = SurfaceChar(Gridiv, (/c_Wall_rhoCp1, c_Wall_rhoCp2, c_Wall_rhoCp3, c_Wall_rhoCp4, c_Wall_rhoCp5/))
         zibld = SurfaceChar(Gridiv, &
                             (/c_Internal_thick1, c_Internal_thick2, c_Internal_thick3, c_Internal_thick4, c_Internal_thick5/))
         kibld = SurfaceChar(Gridiv, (/c_Internal_k1, c_Internal_k2, c_Internal_k3, c_Internal_k4, c_Internal_k5/))
         ribld = SurfaceChar(Gridiv, &
                             (/c_Internal_rhoCp1, c_Internal_rhoCp2, c_Internal_rhoCp3, c_Internal_rhoCp4, c_Internal_rhoCp5/))

         nroom = SurfaceChar(Gridiv, c_nroom)
         alb_ibld = SurfaceChar(Gridiv, c_alb_ibld)
         em_ibld = SurfaceChar(Gridiv, c_em_ibld)
         CH_iwall = SurfaceChar(Gridiv, c_CH_iwall)
         CH_iroof = SurfaceChar(Gridiv, c_CH_iroof)
         CH_ibld = SurfaceChar(Gridiv, c_CH_ibld)
      ENDIF

      !For other surfaces, only one ESTM class
      DO iv = ConifSurf, nsurfIncSnow
         zSurf_SUEWSsurfs(:, iv) = SurfaceChar(Gridiv, (/c_Surf_thick1(iv), c_Surf_thick2(iv), c_Surf_thick3(iv), &
                                                         c_Surf_thick4(iv), c_Surf_thick5(iv)/))
         kSurf_SUEWSsurfs(:, iv) = SurfaceChar(Gridiv, (/c_Surf_k1(iv), c_Surf_k2(iv), c_Surf_k3(iv), &
                                                         c_Surf_k4(iv), c_Surf_k5(iv)/))
         rSurf_SUEWSsurfs(:, iv) = SurfaceChar(Gridiv, (/c_Surf_rhoCp1(iv), c_Surf_rhoCp2(iv), c_Surf_rhoCp3(iv), &
                                                         c_Surf_rhoCp4(iv), c_Surf_rhoCp5(iv)/))
      ENDDO

      ! Now combine SUEWS surfaces into ESTM facets
      !Surface fractions for ESTM facets (moved from SUEWS_ESTM_initials HCW 16 Jun 2016)
      !roof = Bldgs
      froof = sfr(BldgSurf)
      !ground = all except Bldgs
      fground = sfr(PavSurf) + sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) + sfr(BsoilSurf) + sfr(WaterSurf)
      !veg = EveTr, DecTr, Grass
      fveg = sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf)

      ! Ground = all except buildings (exclude snow at the moment)
      zground = 0
      kground = 0
      rground = 0
      DO iv = 1, nsurf
         IF (iv /= BldgSurf .AND. fground /= 0) THEN   !Bldgs surface excluded from ground facet
            zground = zground + zSurf_SUEWSsurfs(:, iv)*sfr(iv)/fground   !Normalised by ground fraction
            kground = kground + kSurf_SUEWSsurfs(:, iv)*sfr(iv)/fground   !Normalised by ground fraction
            rground = rground + rSurf_SUEWSsurfs(:, iv)*sfr(iv)/fground   !Normalised by ground fraction
         ELSEIF (fground == 0.) THEN !check fground==0 (or HW==0) scenario to avoid division-by-zero error, TS 21 Jul 2016
            zground = zground + 0.01
            kground = kground + 0.01
            rground = rground + 0.01
            ! PRINT*, zground
            ! PRINT*, kground
            ! PRINT*, rground
         ENDIF
      ENDDO
      ! Roof = buildings
      zroof = zSurf_SUEWSsurfs(:, BldgSurf)
      kroof = kSurf_SUEWSsurfs(:, BldgSurf)
      rroof = rSurf_SUEWSsurfs(:, BldgSurf)

      ! the following initialisation is problematic: TS 01 Mar 2019
      ! what would happen if zground(5)>0? Nground is initialised NOWHERE!
      ! initialise these variables as 5 so if z_sfc(5)>0 happens, the numbers are still correct, TS 06 Aug 2019
      Nground = 5
      Nroof = 5
      Nwall = 5
      Nibld = 5
      DO i = 1, 5
         IF (zground(i) <= 0) THEN
            Nground = i - 1
            EXIT
         ENDIF
      ENDDO
      DO i = 1, 5
         IF (zroof(i) <= 0) THEN
            Nroof = i - 1
            EXIT
         ENDIF
      ENDDO
      DO i = 1, 5
         IF (zwall(i) <= 0) THEN
            Nwall = i - 1
            EXIT
         ENDIF
      ENDDO
      DO i = 1, 5
         IF (zibld(i) <= 0) THEN
            Nibld = i - 1
            EXIT
         ENDIF
      ENDDO
   ENDIF ! ESTM related translation finished here.

   ! ---- AnOHM related ------------------------------
   IF (StorageHeatMethod == 3) THEN
      cpAnOHM(1:nsurf) = SurfaceChar(Gridiv, c_cpAnOHM) ! AnOHM TS
      kkAnOHM(1:nsurf) = SurfaceChar(Gridiv, c_kkAnOHM) ! AnOHM TS
      chAnOHM(1:nsurf) = SurfaceChar(Gridiv, c_chAnOHM) ! AnOHM TS

      ! cp and k are estimated from ESTM coefficients:
      ! cpAnOHM(1:nsurf)=rSurf_SUEWSsurfs(1,1:nsurf)
      ! kkAnOHM(1:nsurf)=kSurf_SUEWSsurfs(1,1:nsurf)
      ! IF ( ir ==1 .AND. iMb ==1) THEN
      !    PRINT*, 'StoreDrainPrm',PavSurf,':'
      !    PRINT'(a10,x,5f10.2)', 'Depth',zSurf_SUEWSsurfs(:,i)
      !    PRINT'(a10,x,5es10.2)', 'RhoCp',rSurf_SUEWSsurfs(:,i)
      !    PRINT'(a10,x,5es10.2)', 'avg_RhoCp',cpAnOHM(i)
      !    PRINT'(a10,x,5es10.2)', 'k',kSurf_SUEWSsurfs(:,i)
      !    PRINT'(a10,x,5es10.2)', 'avg_k',kkAnOHM(i)
      !    PRINT'(a10,x,5es10.2)', 'avg_Ch',chAnOHM(i)
      !
      ! END IF
      ! DO i = 1, nsurf, 1
      !    ! filter out invalid z values
      !    WHERE (  zSurf_SUEWSsurfs(:,i) == -999. ) zSurf_SUEWSsurfs(:,i)=0
      !
      !    ! cp: weight-averaged by depth
      !    cpAnOHM(i)=DOT_PRODUCT(rSurf_SUEWSsurfs(:,i),zSurf_SUEWSsurfs(:,i))/SUM(zSurf_SUEWSsurfs(:,i))
      !   !  IF ( i==PavSurf .AND. ir ==1 .AND. iMb ==1) THEN
      !   !     PRINT*, 'StoreDrainPrm',i,':'
      !   !     PRINT'(a10,x,5f10.2)', 'Depth',zSurf_SUEWSsurfs(:,i)
      !   !     PRINT'(a10,x,5es10.2)', 'RhoCp',rSurf_SUEWSsurfs(:,i)
      !   !     PRINT'(a10,x,5es10.2)', 'avg_RhoCp',cpAnOHM(i)
      !    !
      !   !  END IF
      !
      !    ! 1/k: weight-averaged by depth
      !    kkAnOHM(i)=DOT_PRODUCT(1/kSurf_SUEWSsurfs(:,i),zSurf_SUEWSsurfs(:,i))/SUM(zSurf_SUEWSsurfs(:,i))
      !    kkAnOHM(i)=1/kkAnOHM(i)
      !   !  IF ( i==PavSurf .AND. ir ==1 .AND. iMb ==1) THEN
      !   !     PRINT'(a10,x,5es10.2)', 'k',kSurf_SUEWSsurfs(:,i)
      !   !     PRINT'(a10,x,5es10.2)', 'avg_k',kkAnOHM(i)
      !   !     PRINT'(a10,x,5es10.2)', 'avg_Ch',chAnOHM(i)
      !    !
      !   !     PRINT'(a10,x,7f10.2)', 'fractions:',sfr
      !    !
      !   !  END IF
      !
      !
      !    ! restore invalid z values
      !    WHERE (  zSurf_SUEWSsurfs(:,i) == 0 ) zSurf_SUEWSsurfs(:,i)=nan
      !   !  IF ( i==PavSurf .AND. ir ==1 .AND. iMb ==1) THEN
      !   !     PRINT'(a10,x,5f10.2)', 'Depth',zSurf_SUEWSsurfs(:,i)
      !   !     PRINT*, '*****************'
      !   !  END IF
      !
      ! END DO
   END IF

   ! ---- QF coeffs (was in SUEWS_SAHP.f95, subroutine SAHP_Coefs)
   BaseTHDD = -999 ! Initialise QF coeffs
   QF_A = 0
   QF_B = 0
   QF_C = 0
   AH_min = 0
   T_CRITIC_Heating = 0
   T_CRITIC_Cooling = 0
   AH_slope_Heating = 0
   AH_slope_Cooling = 0

   BaseTHDD = SurfaceChar(Gridiv, c_BaseTHDD)
   QF_A = SurfaceChar(Gridiv, (/c_QF_A1, c_QF_A2/))
   QF_B = SurfaceChar(Gridiv, (/c_QF_B1, c_QF_B2/))
   QF_C = SurfaceChar(Gridiv, (/c_QF_C1, c_QF_C2/))
   AH_min = SurfaceChar(Gridiv, (/c_AHMin_WD, c_AHMin_WE/))
   AH_slope_Heating = SurfaceChar(Gridiv, (/c_AHSlopeHeating_WD, c_AHSlopeHeating_WE/))
   AH_slope_Cooling = SurfaceChar(Gridiv, (/c_AHSlopeCooling_WD, c_AHSlopeCooling_WE/))
   T_Critic_Heating = SurfaceChar(Gridiv, (/c_TCriticHeating_WD, c_TCriticHeating_WE/))
   T_Critic_Cooling = SurfaceChar(Gridiv, (/c_TCriticCooling_WD, c_TCriticCooling_WE/))
   EnProfWD = SurfaceChar(Gridiv, c_EnProfWD)
   EnProfWE = SurfaceChar(Gridiv, c_EnProfWE)
   CO2mWD = SurfaceChar(Gridiv, c_CO2mWD)
   CO2mWE = SurfaceChar(Gridiv, c_CO2mWE)
   TraffProfWD = SurfaceChar(Gridiv, c_TraffProfWD)
   TraffProfWE = SurfaceChar(Gridiv, c_TraffProfWE)
   PopProfWD = SurfaceChar(Gridiv, c_PopProfWD)
   PopProfWE = SurfaceChar(Gridiv, c_PopProfWE)
   MinQFMetab = SurfaceChar(Gridiv, c_MinQFMetab)
   MaxQFMetab = SurfaceChar(Gridiv, c_MaxQFMetab)
   MinFCMetab = SurfaceChar(Gridiv, c_MinFCMetab)
   MaxFCMetab = SurfaceChar(Gridiv, c_MaxFCMetab)
   FrFossilFuel_heat = SurfaceChar(Gridiv, c_FrFossilFuel_heat)
   FrFossilFuel_NonHeat = SurfaceChar(Gridiv, c_FrFossilFuel_NonHeat)
   EF_umolCO2perJ = SurfaceChar(Gridiv, c_EF_umolCO2perJ)
   EnEF_v_Jkm = SurfaceChar(Gridiv, c_EnEF_v_Jkm)
   FcEF_v_kgkm = SurfaceChar(Gridiv, (/c_FcEF_v_kgkmWD, c_FcEF_v_kgkmWE/))
   CO2PointSource = SurfaceChar(Gridiv, c_CO2PointSource)
   TrafficUnits = SurfaceChar(Gridiv, c_TrafficUnits)

   ! ---- Irrigation
   Ie_start = INT(SurfaceChar(Gridiv, c_IeStart))
   Ie_end = INT(SurfaceChar(Gridiv, c_IeEnd))
   InternalWaterUse_h = SurfaceChar(Gridiv, c_IntWU)
   Faut = SurfaceChar(Gridiv, c_Faut)
   Ie_a = SurfaceChar(Gridiv, c_Ie_a)   !Automatic irrigation model coefficients [mm d-1]; [mm d-1 degC-1]; [mm d-2]
   Ie_m = SurfaceChar(Gridiv, c_Ie_m)   !Manual irrigation model coefficients [mm d-1]; [mm d-1 degC-1]; [mm d-2]
   DayWat = SurfaceChar(Gridiv, c_DayWat)
   DayWatPer = SurfaceChar(Gridiv, c_DayWatPer)

   ! ---- Hourly profiles
   AHProf_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfEnUseWD)       ! Anthropogenic heat, weekdays
   AHProf_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HrProfEnUseWE)       ! Anthropogenic heat, weekends
   WUProfM_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfWUManuWD)      ! Water use, manual, weekdays
   WUProfM_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HrProfWUManuWE)      ! Water use, manual, weekends
   WUProfA_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfWUAutoWD)      ! Water use, automatic, weekdays
   WUProfA_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HrProfWUAutoWE)      ! Water use, automatic, weekends
   SnowProf_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfSnowCWD)       ! Snow clearing, weekdays
   SnowProf_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HrProfSnowCWE)       ! Snow clearing, weekends
   HumActivity_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfHumActivityWD) ! Human activity, weekdays
   HumActivity_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HrProfHumActivityWE) ! Human activity, weekends
   TraffProf_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HrProfTraffWD)       ! Traffic, weekdays
   TraffProf_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HRProfTraffWE)       ! Traffic, weekends
   PopProf_24hr(0:23, 1) = SurfaceChar(Gridiv, c_HRProfPopWD)         ! Population, weekdays
   PopProf_24hr(0:23, 2) = SurfaceChar(Gridiv, c_HRProfPopWE)         ! Population, weekends

   ! ---- Profiles at the resolution of model time step
   ! AHProf_tstep(:,1)      = TstepProfiles(Gridiv,cTP_EnUseWD,:)       ! Anthropogenic heat, weekdays
   ! AHProf_tstep(:,2)      = TstepProfiles(Gridiv,cTP_EnUseWE,:)       ! Anthropogenic heat, weekends
   ! WUProfM_tstep(:,1)     = TstepProfiles(Gridiv,cTP_WUManuWD,:)      ! Water use, manual, weekdays
   ! WUProfM_tstep(:,2)     = TstepProfiles(Gridiv,cTP_WUManuWE,:)      ! Water use, manual, weekends
   ! WUProfA_tstep(:,1)     = TstepProfiles(Gridiv,cTP_WUAutoWD,:)      ! Water use, automatic, weekdays
   ! WUProfA_tstep(:,2)     = TstepProfiles(Gridiv,cTP_WUAutoWE,:)      ! Water use, automatic, weekends
   ! HumActivity_tstep(:,1) = TstepProfiles(Gridiv,cTP_HumActivityWD,:) ! Human activity, weekdays
   ! HumActivity_tstep(:,2) = TstepProfiles(Gridiv,cTP_HumActivityWE,:) ! Human activity, weekends
   ! TraffProf_tstep(:,1)   = TstepProfiles(Gridiv,cTP_TraffProfWD,:)   !Traffic, weekdays
   ! TraffProf_tstep(:,2)   = TstepProfiles(Gridiv,cTP_TraffProfWE,:)   !Traffic, weekends
   ! PopProf_tstep(:,1)     = TstepProfiles(Gridiv,cTP_PopProfWD,:)     !Population, weekdays
   ! PopProf_tstep(:,2)     = TstepProfiles(Gridiv,cTP_PopProfWE,:)     !Population, weekends

   ! ---- Within-grid water distribution
   ! N.B. Rows and columns of WaterDist are the other way round to the input info
   !! Model currently does not include above-ground flow from the Water surface
   !! - Probably should adjust WaterDist to have nsurf columns so that Water can behave like the other surfaces.
   ! Model returns an error if both ToRunoff and ToSoilStore are non-zero (in CodeMatchDist)
   ! For impervious surfaces, water goes to runoff; for pervious surfaces, water goes to soilstore
   WaterDist(PavSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToPaved(1:(nsurf - 1)))
   WaterDist(BldgSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToBldgs(1:(nsurf - 1)))
   WaterDist(ConifSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToEveTr(1:(nsurf - 1)))
   WaterDist(DecidSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToDecTr(1:(nsurf - 1)))
   WaterDist(GrassSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToGrass(1:(nsurf - 1)))
   WaterDist(BSoilSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToBSoil(1:(nsurf - 1)))
   WaterDist(WaterSurf, 1:(nsurf - 1)) = SurfaceChar(Gridiv, c_WGToWater(1:(nsurf - 1)))
   ! Runoff or SoilStore row   !!Change later to allow both Runoff and SoilStore
   DO iv = 1, (nsurf - 1)
      IF (SurfaceChar(Gridiv, c_WGToRunoff(iv)) /= 0) THEN
         WaterDist((nsurf + 1), iv) = SurfaceChar(Gridiv, c_WGToRunoff(iv))
      ELSE
         WaterDist((nsurf + 1), iv) = SurfaceChar(Gridiv, c_WGToSoilStore(iv))
      ENDIF
   ENDDO

   ! Access required DailyState variables for the current grid (moved HCW 26 Jun 2015)
   ! HDD(:,:)    = HDD_grids(:,:,Gridiv)
   ! GDD(:,:)    = GDD_grids(:,:,Gridiv)
   ! LAI(:,:)    = LAI_grids(:,:,Gridiv)
   ! WUDay(:,:) = WUDay_grids(:,:,Gridiv)
   ! AlbDecTr(:) = AlbDecTr_grids(:,Gridiv)
   ! DecidCap(:) = DecidCap_grids(:,Gridiv)
   ! Porosity(:) = Porosity_grids(:,Gridiv)
   ! AlbEveTr(:) = AlbEveTr_grids(:,Gridiv)
   ! AlbGrass(:) = AlbGrass_grids(:,Gridiv)
   ! SnowAlb = ModelDailyState(Gridiv, cMDS_SnowAlb)

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
   ! load snow related properties for NARP
   if (snowuse == 1) NARP_EMIS_SNOW = SurfaceChar(Gridiv, c_SnowEmis)
   !NARP_CONFIGURATION if net radiation is to be modelled
   IF (NetRadiationMethod > 0) THEN
      NARP_LAT = SurfaceChar(Gridiv, c_lat)
      NARP_LONG = SurfaceChar(Gridiv, c_lng)    ! New sun_position_v2 use degrees FL
      NARP_YEAR = INT(SurfaceChar(Gridiv, c_Year))
      NARP_TZ = TIMEZONE                           !not every 5-min
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
   IF (ir == 0) THEN
      !write(*,*) 'This should be seen only when called from InitialState and ir is 0. ir:',ir

      ! =============================================================================
      ! === Translate inputs from ModelDailyState to variable names used in model ===
      ! =============================================================================

      ! Get id_prev from ModelDailyState
      id_prev = INT(ModelDailyState(Gridiv, cMDS_id_prev))

      SnowfallCum = ModelDailyState(Gridiv, cMDS_SnowfallCum)
      SnowAlb = ModelDailyState(Gridiv, cMDS_SnowAlb)

      porosity_id = ModelDailyState(Gridiv, cMDS_porosity)
      albDecTr_id = ModelDailyState(Gridiv, cMDS_albDecTr)
      albEveTr_id = ModelDailyState(Gridiv, cMDS_albEveTr)
      albGrass_id = ModelDailyState(Gridiv, cMDS_albGrass)
      DecidCap_id = ModelDailyState(Gridiv, cMDS_DecidCap)

      DecidCap_id_grids(Gridiv) = DecidCap_id
      albDecTr_id_grids(Gridiv) = albDecTr_id
      albEveTr_id_grids(Gridiv) = albEveTr_id
      albGrass_id_grids(Gridiv) = albGrass_id
      porosity_id_grids(Gridiv) = porosity_id

      ! ---- Phenology
      ! ---- LAI
      LAI_id = 0
      LAI_id(ivConif) = ModelDailyState(Gridiv, cMDS_LAIInitialEveTr)
      LAI_id(ivDecid) = ModelDailyState(Gridiv, cMDS_LAIInitialDecTr)
      LAI_id(ivGrass) = ModelDailyState(Gridiv, cMDS_LAIInitialGrass)

      ! GDD_id: GDD Values for one day
      GDD_id = ModelDailyState(Gridiv, cMDS_GDD1_0)
      ! SDD_id: SDD Values for one day
      SDD_id = ModelDailyState(Gridiv, cMDS_GDD2_0)
      ! Tmin, Tmax: daily minimum and maximum temperatures
      Tmin_id = ModelDailyState(Gridiv, cMDS_GDDMin)
      Tmax_id = ModelDailyState(Gridiv, cMDS_GDDMax)
      ! length of daylight
      lenDay_id = 0

      ! ---- Heating degree days, HDD
      ! HDD = 0
      ! HDD(id_prev,1)   = ModelDailyState(Gridiv,cMDS_HDD1)      ! 1 = Heating
      ! HDD(id_prev,2)   = ModelDailyState(Gridiv,cMDS_HDD2)      ! 2 = Cooling
      ! HDD(id_prev-3,3) = ModelDailyState(Gridiv,cMDS_TempCOld3) ! 3 will become average
      ! HDD(id_prev-2,3) = ModelDailyState(Gridiv,cMDS_TempCOld2)
      ! HDD(id_prev-1,3) = ModelDailyState(Gridiv,cMDS_TempCOld1)
      ! HDD(id_prev,3)   = ModelDailyState(Gridiv,cMDS_TempC)
      ! 4 = 5 day running mean
      ! 5 = daily precip total
      ! HDD(id_prev,6) = ModelDailyState(Gridiv,cMDS_DaysSinceRain) ! 6 = days since rain

      ! ---- Heating degree days, HDD_id: HDD Values for one day
      HDD_id(1:6) = 0

      ! HDD_id(1)=ModelDailyState(Gridiv,cMDS_HDD1)
      ! HDD_id(2)=ModelDailyState(Gridiv,cMDS_HDD2)
      ! HDD_id(3)=ModelDailyState(Gridiv,cMDS_TempC)
      ! ! 4 = 5 day running mean
      ! ! 5 = daily precip total
      ! HDD_id(6) = ModelDailyState(Gridiv,cMDS_DaysSinceRain)

      ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
      HDD_id_grids(:, Gridiv) = HDD_id(:)
      GDD_id_grids(:, Gridiv) = GDD_id(:)
      SDD_id_grids(:, Gridiv) = SDD_id(:)
      Tmin_id_grids(Gridiv) = Tmin_id
      Tmax_id_grids(Gridiv) = Tmax_id
      lenDay_id_grids(Gridiv) = lenDay_id
      LAI_id_grids(:, Gridiv) = LAI_id(:)

      ! daily water use
      WUDay_id = 0
      WUDay_id_grids(:, Gridiv) = WUDay_id(:)

      ! AlbDecTr_grids(:,Gridiv) = AlbDecTr(:)
      ! AlbEveTr_grids(:,Gridiv) = AlbEveTr(:)
      ! AlbGrass_grids(:,Gridiv) = AlbGrass(:)
      ! DecidCap_grids(:,Gridiv) = DecidCap(:)
      ! Porosity_grids(:,Gridiv) = Porosity(:)

      ! ---- Snow density of each surface
      SnowDens(1:nsurf) = ModelDailyState(Gridiv, cMDS_SnowDens(1:nsurf))

      ! =============================================================================
      ! === Translate inputs from ModelOutputData to variable names used in model ===
      ! =============================================================================
      ! ---- Above-ground state
      state_id(1:nsurf) = ModelOutputData(0, cMOD_State(1:nsurf), Gridiv)
      !     stateDay(0,Gridiv,1:nsurf) = ModelOutputData(0,cMOD_State(1:nsurf),Gridiv)
      ! ---- Below-ground
      soilstore_id(1:nsurf) = ModelOutputData(0, cMOD_SoilState(1:nsurf), Gridiv)
      !     soilmoistDay(0,Gridiv,1:nsurf) = ModelOutputData(0,cMOD_SoilState(1:nsurf),Gridiv)
      ! ---- Snow fraction
      SnowFrac(1:nsurf) = ModelOutputData(0, cMOD_SnowFrac(1:nsurf), Gridiv)
      ! ---- Snow water equivalent in SnowPack
      SnowPack(1:nsurf) = ModelOutputData(0, cMOD_SnowPack(1:nsurf), Gridiv)
      ! ---- Liquid (melted) water in SnowPack
      SnowWater(1:nsurf) = ModelOutputData(0, cMOD_SnowWaterState(1:nsurf), Gridiv)

   ENDIF  !ir = 0
   !=================================================================================

   !=========================== Write FileChoices.txt ===============================
   !=================================================================================
   ! Do once per grid per year (was in SUEWS_Initial.f95)
   IF (ir == 1 .AND. iMB == 1) THEN   !For first row of first block only
      !write(*,*) 'Writing to FileChoices for first chunk of met data per year per grid'
      FileChoices = TRIM(FileOutputPath)//TRIM(FileCode)//'_FileChoices.txt'
      OPEN (12, file=FileChoices, position='append')

      WRITE (grid_txt, '(I5)') INT(SurfaceChar(Gridiv, c_Grid))
      WRITE (year_txt, '(I4)') INT(SurfaceChar(Gridiv, c_Year))
      WRITE (SsG_YYYY, '(A12)') TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(ADJUSTL(year_txt))

      !write(12,*) '--------------------------------------------------------------------------------'
      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Surface characteristics'//' -----'
      ! Characteristics that apply to some or all surface types
      WRITE (12, '(8a10,a16)') 'Paved', 'Bldgs', 'EveTr', 'DecTr', 'Grass', 'BSoil', 'Water', 'Snow', ' SurfType'
      WRITE (12, 120) (sfr(iv), iv=1, nsurf), FCskip, ' SurfFr'
      WRITE (12, 120) FCskip, FCskip, IrrFracConif, IrrFracDecid, IrrFracGrass, FCskip, FCskip, FCskip, ' IrrFr'
      WRITE (12, 120) FCskip, FCskip, WUAreaEveTr_m2, WUAreaDecTr_m2, WUAreaGrass_m2, FCskip, FCskip, FCskip, ' WaterUseArea'
      WRITE (12, 120) FCskip, BldgH, EveTreeH, DecTreeH, FCskip, FCskip, FCskip, FCskip, ' H'
      WRITE (12, 120) FCskip, FAIBldg, FAIEveTree, FAIDecTree, FCskip, FCskip, FCskip, FCskip, ' FAI'
      WRITE (12, 120) FCskip, FCskip, AlbMin_EveTr, AlbMin_DecTr, AlbMin_Grass, FCskip, FCskip, SnowAlbMin, ' AlbedoMin'
      WRITE (12, 120) FCskip, FCskip, AlbMax_EveTr, AlbMax_DecTr, AlbMax_Grass, FCskip, FCskip, SnowAlbMax, ' AlbedoMax'
      !write(12,120) (alb(iv),iv=1,nsurf),SnowAlb, ' Albedo'   ! This is instantaneous value (not provided as input)
      WRITE (12, 120) (emis(iv), iv=1, nsurf), emis_snow, ' Emissivity'
      WRITE (12, 120) FCskip, FCskip, (baseT(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' BaseT'
      WRITE (12, 120) FCskip, FCskip, (baseTe(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' BaseTe'
      WRITE (12, 120) (StoreDrainPrm(1, iv), iv=1, nsurf), FCskip, ' StorageMin'
      WRITE (12, 120) (StoreDrainPrm(5, iv), iv=1, nsurf), FCskip, ' StorageMax'
      WRITE (12, 120) (WetThresh(iv), iv=1, nsurf), FCskip, ' WetThreshold'
      WRITE (12, 120) (StateLimit(iv), iv=1, nsurf), FCskip, ' StateLimit'
      WRITE (12, 120) (StoreDrainPrm(2, iv), iv=1, nsurf), FCskip, ' DrainageEq'    !real
      WRITE (12, 120) (StoreDrainPrm(3, iv), iv=1, nsurf), FCskip, ' DrainageCoef1'
      WRITE (12, 120) (StoreDrainPrm(4, iv), iv=1, nsurf), FCskip, ' DrainageCoef2'
      WRITE (12, 120) FCskip, FCskip, (GDDFull(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' GDDFull'
      WRITE (12, 120) FCskip, FCskip, (SDDFull(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' SDDFull'
      WRITE (12, 120) FCskip, FCskip, (LAImin(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' LAIMin'
      WRITE (12, 120) FCskip, FCskip, (LAImax(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' LAIMax'
      WRITE (12, 120) FCskip, FCskip, FCskip, PorMin_dec, FCSkip, FCskip, FCskip, FCskip, ' PorosityMin'
      WRITE (12, 120) FCskip, FCskip, FCskip, PorMax_dec, FCSkip, FCskip, FCskip, FCskip, ' PorosityMax'
      WRITE (12, '(2f10.3,3i10,  3f10.3,a16)') FCskip, FCskip, LAItype(1:nvegsurf), FCskip, FCskip, FCskip, ' LAIEq'   !integer
      WRITE (12, '(2f10.3,3f10.5,3f10.3,a16)') FCskip, FCskip, LAIPower(1, 1:nvegsurf), FCskip, FCskip, FCskip, ' LAI_LeafGP1'
      WRITE (12, '(2f10.3,3f10.5,3f10.3,a16)') FCskip, FCskip, LAIPower(2, 1:nvegsurf), FCskip, FCskip, FCskip, ' LAI_LeafGP2'
      WRITE (12, '(2f10.3,3f10.5,3f10.3,a16)') FCskip, FCskip, LAIPower(3, 1:nvegsurf), FCskip, FCskip, FCskip, ' LAI_LeafOP1'
      WRITE (12, '(2f10.3,3f10.5,3f10.3,a16)') FCskip, FCskip, LAIPower(4, 1:nvegsurf), FCskip, FCskip, FCskip, ' LAI_LeafOP2'
      WRITE (12, 120) FCskip, FCskip, (MaxConductance(iv), iv=1, nvegsurf), FCskip, FCskip, FCskip, ' MaxCond'
      WRITE (12, 120) (SoilDepth(iv), iv=1, (nsurf - 1)), FCskip, FCskip, ' SoilDepth'
      WRITE (12, 120) (SoilStoreCap(iv), iv=1, (nsurf - 1)), FCskip, FCskip, ' SoilStoreCap'
      WRITE (12, '(6f10.5,2f10.3,a16)') (SatHydraulicConduct(iv), iv=1, (nsurf - 1)), FCskip, FCskip, ' SatHydraulicConduct'
      ! Not currently coded, but add these later: SoilDensity, InfiltrationRate, OBS_SMDept, OBS_SMCap, OBS_SoilNotRocks
      WRITE (12, 120) (SnowPackLimit(iv), iv=1, (nsurf - 1)), FCskip, FCskip, ' SnowLimPatch'
      WRITE (12, 120) SnowLimPaved, SnowLimBldg, FCskip, FCskip, FCskip, FCskip, FCskip, FCskip, ' SnowLimRemove'
      WRITE (12, 120) (OHM_coef(1:nsurf, 1, 1)), OHM_coef(nsurf + 1, 1, 1), ' OHM_a1_Sum_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 2, 1)), OHM_coef(nsurf + 1, 2, 1), ' OHM_a1_Sum_Dry'
      WRITE (12, 120) (OHM_coef(1:nsurf, 3, 1)), OHM_coef(nsurf + 1, 3, 1), ' OHM_a1_Win_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 4, 1)), OHM_coef(nsurf + 1, 4, 1), ' OHM_a1_Win_Dry'
      WRITE (12, 120) (OHM_coef(1:nsurf, 1, 2)), OHM_coef(nsurf + 1, 1, 2), ' OHM_a2_Sum_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 2, 2)), OHM_coef(nsurf + 1, 2, 2), ' OHM_a2_Sum_Dry'
      WRITE (12, 120) (OHM_coef(1:nsurf, 3, 2)), OHM_coef(nsurf + 1, 3, 2), ' OHM_a2_Win_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 4, 2)), OHM_coef(nsurf + 1, 4, 2), ' OHM_a2_Win_Dry'
      WRITE (12, 120) (OHM_coef(1:nsurf, 1, 3)), OHM_coef(nsurf + 1, 1, 3), ' OHM_a3_Sum_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 2, 3)), OHM_coef(nsurf + 1, 2, 3), ' OHM_a3_Sum_Dry'
      WRITE (12, 120) (OHM_coef(1:nsurf, 3, 3)), OHM_coef(nsurf + 1, 3, 3), ' OHM_a3_Win_Wet'
      WRITE (12, 120) (OHM_coef(1:nsurf, 4, 3)), OHM_coef(nsurf + 1, 4, 3), ' OHM_a3_Win_Dry'
      WRITE (12, 120) (OHM_threshSW(1:nsurf)), OHM_threshSW(nsurf + 1), ' OHMthreshold_SW'
      WRITE (12, 120) (OHM_threshWD(1:nsurf)), OHM_threshWD(nsurf + 1), ' OHMthreshold_WD'

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Snow parameters'//' -----'
      WRITE (12, '(a12,11a10)') 'Grid', 'RadMeltF', 'TempMeltF', 'tau_a', 'tau_f', 'PLimAlb', 'SDensMin', 'SDensMax', &
         'tau_r', 'CRWMin', 'CRWMax', 'PLimSnow'
      WRITE (12, '(a12,11f10.4)') SsG_YYYY, RadMeltFact, TempMeltFact, tau_a, tau_f, PrecipLimitAlb, SnowDensMin, SnowDensMax, &
         tau_r, CRWmin, CRWmax, PrecipLimit

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Conductance parameters'//' -----'
      WRITE (12, '(a12,12a10)') 'Grid', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'TH', 'TL', 'S1', 'S2', 'Kmax', 'gsModel'
      WRITE (12, '(a12,11f10.3,i3)') SsG_YYYY, G1, G2, G3, G4, G5, G6, TH, TL, S1, S2, Kmax, gsModel

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Energy-use parameters'//' -----'
      WRITE (12, '(a12,11a10)') 'Grid', 'PopDensDaytime', 'BaseTHDD', 'QF_A_WD', 'QF_A_WE', 'QF_B_WD', 'QF_B_WE', 'QF_C_WD', &
         'QF_C_WE', 'AH_Min', 'AH_Slope', 'T_critic_Heating'
      WRITE (12, '(a12,11f10.3)') SsG_YYYY, PopDensDaytime, BaseTHDD, QF_A(1:2), QF_B(1:2), QF_C(1:2), &
         AH_Min, AH_Slope_Heating, T_critic_Heating

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Water-use parameters'//' -----'
      WRITE (12, '(a12,10a10)') 'Grid', 'IeStart', 'IeEnd', 'IntWatUse', 'Faut', &
         'Ie_a1', 'Ie_a2', 'Ie_a3', 'Ie_m1', 'Ie_m2', 'Ie_m3'
      WRITE (12, '(a12,2i10,8f10.3)') SsG_YYYY, Ie_start, Ie_end, InternalWaterUse_h, Faut, &
         Ie_a(1:3), Ie_m(1:3)

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Weekly profiles'//' -----'
      WRITE (12, '(a12,7a10,  a16)') 'Grid', '1_Sun', '2_Mon', '3_Tue', '4_Wed', '5_Thu', '6_Fri', '7_Sat', ' DayOfWeek'
      WRITE (12, '(a12,7f10.3,a16)') SsG_YYYY, DayWat(1:7), ' Irr allowed'
      WRITE (12, '(a12,7f10.3,a16)') SsG_YYYY, DayWatPer(1:7), ' Frac properties'

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Hourly profiles'//' -----'
      WRITE (12, '(a12,24i10,a20)') 'Grid', (iv, iv=0, 23), 'HourOfDay'
      WRITE (12, 121) SsG_YYYY, AHProf_24hr(0:23, 1), ' Anthrop heat WD'
      WRITE (12, 121) SsG_YYYY, AHProf_24hr(0:23, 2), ' Anthrop heat WE'
      WRITE (12, 121) SsG_YYYY, WUProfM_24hr(0:23, 1), ' Manual water use WD'
      WRITE (12, 121) SsG_YYYY, WUProfM_24hr(0:23, 2), ' Manual water use WE'
      WRITE (12, 121) SsG_YYYY, WUProfA_24hr(0:23, 1), ' Auto. water use WD'
      WRITE (12, 121) SsG_YYYY, WUProfA_24hr(0:23, 2), ' Auto. water use WE'
      WRITE (12, 121) SsG_YYYY, SnowProf_24hr(0:23, 1), ' Snow clearing WD'
      WRITE (12, 121) SsG_YYYY, SnowProf_24hr(0:23, 2), ' Snow clearing WE'

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Within-grid water distribution'//' -----'
      WRITE (12, '(9a10)') 'ToPaved', 'ToBldgs', 'ToEveTr', 'ToDecTr', 'ToGrass', 'ToBSoil', 'ToWater', 'ToROorSS'

      DO iv = 1, (nsurf - 1)
         WRITE (12, '(8f10.4)') (WaterDist(j, iv), j=1, nsurf + 1)
      ENDDO

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Other parameters'//' -----'
      WRITE (12, '(a12,7a10)') 'Grid', 'FlowChange', 'ROToWater', 'PipeCap', &   ! Water-related
         'DrRate', 'Cover', 'MaxRes', &   ! LUMPS-related
         'Trans'   ! NARP-related
      WRITE (12, '(a12,7f10.3)') SsG_YYYY, FlowChange, RunoffToWater, PipeCapacity, &
         DRAINRT, RAINCOVER, RAINMAXRES, &
         Trans_Site

      WRITE (12, *) '----- '//TRIM(ADJUSTL(SsG_YYYY))//' Site parameters'//' -----'
      WRITE (12, '(a12,9a10)') &
         'Grid', 'lat', 'lon', 'tz', 'alt', 'SurfA_ha', 'z', 'PopDensNighttime', 'z0_input', 'zd_input', 'StartDLS', 'EndDLS'
      WRITE (12, '(a12,4f10.4,f10.2,4f10.4,2i10)') &
         SsG_YYYY, lat, lng*(-1.0), timezone, alt, SurfaceArea_ha, z, PopDensNighttime, z0m_in, zdm_in, &
         startDLS, endDLS! DayLightSavingDay(1:2)

      WRITE (12, *) ''

      CLOSE (12)

      !==============================================================================
      ! Check input values are reasonable ===========================================

      ! Coefficients for anthropogenic heat models ----------------------------------
      IF (EmissionsMethod == 1) THEN   !Loridan et al. (2011) calculation
         IF (AH_min(1) == 0 .AND. Ah_slope_Heating(1) == 0 .AND. T_Critic_Heating(1) == 0) THEN
            CALL ErrorHint(53, 'Check QF calculation coefficients.', notUsed, notUsed, EmissionsMethod)
         ENDIF

      ELSEIF (EmissionsMethod == 2) THEN   !Jarvi et al. (2011) calculation
         IF (SUM(QF_A) == 0 .AND. SUM(QF_B) == 0 .AND. SUM(QF_C) == 0) THEN
            CALL ErrorHint(54, 'Check QF calculation coefficients.', notUsed, notUsed, EmissionsMethod)
         ENDIF
      ENDIF

      ! Morphometric parameters -----------------------------------------------------
      IF (RoughLenMomMethod == 1) THEN   !z0, zd values provided in input file
         ! Check z0m and zd are reasonable
         IF (z0m < 0.00001) CALL ErrorHint(1, 'z0 value provided is very small (RoughLenMomMethod=1).', z0m, notUsed, GridID)
         IF (zdm < 0.00001) CALL ErrorHint(1, 'zd value provided is very small (RoughLenMomMethod=1).', zdm, notUsed, GridID)
         zzd = z - zdm
      ELSEIF (RoughLenMomMethod == 3) THEN   !z0, zd calculated using FAI provided in input file
         ! Check FAIs reasonable
         IF (FAIBLdg < 0) CALL ErrorHint(1, &
                                         'FAI_Bldgs value provided is very small (RoughLenMomMethod=3)', &
                                         FAIBldg, notUsed, GridID)
         IF (FAITree < 0) CALL ErrorHint(1, &
                                         'FAI_EveTr/DecTr value provided is very small (RoughLenMomMethod=3)', &
                                         FAITree, notUsed, GridID)
      ENDIF

   ENDIF   !End for first row of first block only ===================================

   !=================================================================================
   !For each row of the met forcing file (ir), translate correct info for each grid
   ! into model variables
   IF (ir > 0) THEN
      ! =============================================================================
      ! === Translate met data from MetForcingData to variable names used in model ==
      ! =============================================================================
      iy = INT(MetForcingData(ir, 1, Gridiv)) !Integer variables
      id = INT(MetForcingData(ir, 2, Gridiv))
      it = INT(MetForcingData(ir, 3, Gridiv))
      imin = INT(MetForcingData(ir, 4, Gridiv))
      isec = 0 ! NOT used by SUEWS but by WRF-SUEWS via the cal_main interface
      qn1_obs = MetForcingData(ir, 5, Gridiv)      !Real values (kind(1d0))
      qh_obs = MetForcingData(ir, 6, Gridiv)
      qe_obs = MetForcingData(ir, 7, Gridiv)
      qs_obs = MetForcingData(ir, 8, Gridiv)
      qf_obs = MetForcingData(ir, 9, Gridiv)
      avu1 = MetForcingData(ir, 10, Gridiv)
      avrh = MetForcingData(ir, 11, Gridiv)
      Temp_C = MetForcingData(ir, 12, Gridiv)
      Press_hPa = MetForcingData(ir, 13, Gridiv)
      Precip = MetForcingData(ir, 14, Gridiv)
      avkdn = MetForcingData(ir, 15, Gridiv)
      snowFrac_obs = MetForcingData(ir, 16, Gridiv)
      ldown_obs = MetForcingData(ir, 17, Gridiv)
      fcld_obs = MetForcingData(ir, 18, Gridiv)
      wu_m3 = MetForcingData(ir, 19, Gridiv)
      xsmd = MetForcingData(ir, 20, Gridiv)
      LAI_obs = MetForcingData(ir, 21, Gridiv)
      kdiff = MetForcingData(ir, 22, Gridiv)
      kdir = MetForcingData(ir, 23, Gridiv)
      wdir = MetForcingData(ir, 24, Gridiv)

      ! get qn1 memory for previous time steps
      dqndt = dqndt_grids(Gridiv)
      qn1_av = qn1_av_grids(Gridiv)
      tair_av = tair_av_grids(Gridiv)

      IF (SnowUse == 1) THEN
         dqnsdt = dqnsdt_grids(Gridiv)
         qn1_s_av = qn1_s_av_grids(Gridiv)
      ENDIF

      ! added by TS 29 Jun 2018 to remove annual loops in main calculation
      GDD_id = GDD_id_grids(:, Gridiv)
      SDD_id = SDD_id_grids(:, Gridiv)
      Tmin_id = Tmin_id_grids(Gridiv)
      Tmax_id = Tmax_id_grids(Gridiv)
      lenDay_id = lenDay_id_grids(Gridiv)
      HDD_id = HDD_id_grids(:, Gridiv)
      LAI_id = LAI_id_grids(:, Gridiv)
      WUDay_id = WUDay_id_grids(:, Gridiv)

      DecidCap_id = DecidCap_id_grids(Gridiv)
      albDecTr_id = albDecTr_id_grids(Gridiv)
      albEveTr_id = albEveTr_id_grids(Gridiv)
      albGrass_id = albGrass_id_grids(Gridiv)
      porosity_id = porosity_id_grids(Gridiv)

      ! get met array for one grid used in AnOHM
      MetForcingData_grid = MetForcingData(:, :, Gridiv)

      ! Calculate dectime
      dectime = REAL(id - 1, KIND(1d0)) + REAL(it, KIND(1d0))/24 + REAL(imin, KIND(1d0))/(60*24)
      ! Create datetime stamp for error/warnings file
      WRITE (iy_text, '(i4)') iy
      WRITE (id_text, '(i3)') id
      WRITE (it_text, '(i2)') it
      WRITE (imin_text, '(i2)') imin
      !     datetime = TRIM(ADJUSTL(iy_text))//' '//TRIM(ADJUSTL(id_text))//' '//TRIM(ADJUSTL(it_text))//' '//TRIM(ADJUSTL(imin_text))
      WRITE (GridID_text, '(i10)') GridID

      ! =============================================================================
      ! === Translate values from ModelDailyState to variable names used in model ===
      ! =============================================================================
      ! porosity(id) = ModelDailyState(Gridiv,cMDS_porosity)
      ! albDecTr(id) = ModelDailyState(Gridiv,cMDS_albDecTr)
      ! albEveTr(id) = ModelDailyState(Gridiv,cMDS_albEveTr)
      ! albGrass(id) = ModelDailyState(Gridiv,cMDS_albGrass)
      ! DecidCap(id) = ModelDailyState(Gridiv,cMDS_DecidCap)

      ! SnowfallCum is instantaneous values and should be translated at each tstep, TS 17 Sep 2019
      SnowfallCum = ModelDailyState(Gridiv, cMDS_SnowfallCum)
      ! ---- Snow density of each surface
      SnowDens(1:nsurf) = ModelDailyState(Gridiv, cMDS_SnowDens(1:nsurf))
      ! ---- Snow albedo
      SnowAlb = ModelDailyState(Gridiv, cMDS_SnowAlb)

      ! =============================================================================
      ! === Translate values from ModelOutputData to variable names used in model ===
      ! =============================================================================
      ! ---- Above-ground state
      state_id(1:nsurf) = ModelOutputData(ir - 1, cMOD_State(1:nsurf), Gridiv)
      ! ---- Below-ground state
      soilstore_id(1:nsurf) = ModelOutputData(ir - 1, cMOD_SoilState(1:nsurf), Gridiv)
      ! ---- Snow fraction
      SnowFrac(1:nsurf) = ModelOutputData(ir - 1, cMOD_SnowFrac(1:nsurf), Gridiv)
      ! ---- Snow water equivalent in SnowPack
      SnowPack(1:nsurf) = ModelOutputData(ir - 1, cMOD_SnowPack(1:nsurf), Gridiv)
      ! ---- Liquid (melted) water in SnowPack
      SnowWater(1:nsurf) = ModelOutputData(ir - 1, cMOD_SnowWaterState(1:nsurf), Gridiv)

      ! ---- ice fraction
      IceFrac = IceFrac_grids(:, Gridiv)

      !Also translate ESTM forcing data
      IF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
         ! write(*,*) 'Translating ESTM forcing data'
         Ts5mindata(ir, 1:ncolsESTMdata) = ESTMForcingData(ir, 1:ncolsESTMdata, Gridiv)
         Ts5mindata_ir(1:ncolsESTMdata) = ESTMForcingData(ir, 1:ncolsESTMdata, Gridiv)
         CALL ESTM_translate(Gridiv)
      ENDIF

   ENDIF !ir>0   !===================================================================

   ! --------------------------------------------------------------------------------
   ! Check Initial Conditions are reasonable ----------------------------------------
   IF (ir == 1 .AND. iMB == 1) THEN   !For first row of first block only
      CALL CheckInitial
   ENDIF
   ! --------------------------------------------------------------------------------

   ! ======================================================================
   ! write out initial conditions for debugging supy
   if (ir == 1 .and. imb == 1) then
      FileStateInit = TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(ADJUSTL(year_txt))//'_state_init.txt'
      OPEN (12, file=FileStateInit, position='rewind')

      write (12, *) '&state_init'
      write (12, *) 'aerodynamicresistancemethod=', aerodynamicresistancemethod
      write (12, *) 'ah_min=', ah_min
      write (12, *) 'ahprof_24hr=', ahprof_24hr
      write (12, *) 'ah_slope_cooling=', ah_slope_cooling
      write (12, *) 'ah_slope_heating=', ah_slope_heating
      write (12, *) 'alb=', alb
      write (12, *) 'albmax_dectr=', albmax_dectr
      write (12, *) 'albmax_evetr=', albmax_evetr
      write (12, *) 'albmax_grass=', albmax_grass
      write (12, *) 'albmin_dectr=', albmin_dectr
      write (12, *) 'albmin_evetr=', albmin_evetr
      write (12, *) 'albmin_grass=', albmin_grass
      write (12, *) 'alpha_bioco2=', alpha_bioco2
      write (12, *) 'alpha_enh_bioco2=', alpha_enh_bioco2
      write (12, *) 'alt=', alt
      write (12, *) 'avkdn=', avkdn
      write (12, *) 'avrh=', avrh
      write (12, *) 'avu1=', avu1
      write (12, *) 'baset=', baset
      write (12, *) 'basete=', basete
      write (12, *) 'basethdd=', basethdd
      write (12, *) 'beta_bioco2=', beta_bioco2
      write (12, *) 'beta_enh_bioco2=', beta_enh_bioco2
      write (12, *) 'bldgh=', bldgh
      write (12, *) 'capmax_dec=', capmax_dec
      write (12, *) 'capmin_dec=', capmin_dec
      write (12, *) 'chanohm=', chanohm
      write (12, *) 'co2pointsource=', co2pointsource
      write (12, *) 'cpanohm=', cpanohm
      write (12, *) 'crwmax=', crwmax
      write (12, *) 'crwmin=', crwmin
      write (12, *) 'daywat=', daywat
      write (12, *) 'daywatper=', daywatper
      write (12, *) 'dectreeh=', dectreeh
      write (12, *) 'diagnose=', diagnose
      write (12, *) 'diagqn=', diagqn
      write (12, *) 'diagqs=', diagqs
      write (12, *) 'drainrt=', drainrt
      write (12, *) 'dt_since_start=', dt_since_start
      write (12, *) 'dqndt=', dqndt
      write (12, *) 'qn1_av=', qn1_av
      write (12, *) 'dqnsdt=', dqnsdt
      write (12, *) 'qn1_s_av=', qn1_s_av
      write (12, *) 'ef_umolco2perj=', ef_umolco2perj
      write (12, *) 'emis=', emis
      write (12, *) 'emissionsmethod=', emissionsmethod
      write (12, *) 'enef_v_jkm=', enef_v_jkm
      write (12, *) 'enddls=', enddls
      write (12, *) 'evetreeh=', evetreeh
      write (12, *) 'faibldg=', faibldg
      write (12, *) 'faidectree=', faidectree
      write (12, *) 'faievetree=', faievetree
      write (12, *) 'faut=', faut
      write (12, *) 'fcef_v_kgkm=', fcef_v_kgkm
      write (12, *) 'fcld_obs=', fcld_obs
      write (12, *) 'flowchange=', flowchange
      write (12, *) 'frfossilfuel_heat=', frfossilfuel_heat
      write (12, *) 'frfossilfuel_nonheat=', frfossilfuel_nonheat
      write (12, *) 'g1=', g1
      write (12, *) 'g2=', g2
      write (12, *) 'g3=', g3
      write (12, *) 'g4=', g4
      write (12, *) 'g5=', g5
      write (12, *) 'g6=', g6
      write (12, *) 'gdd_id=', gdd_id
      write (12, *) 'gddfull=', gddfull
      write (12, *) 'gridiv=', gridiv
      write (12, *) 'gsmodel=', gsmodel
      write (12, *) 'hdd_id=', hdd_id
      write (12, *) 'humactivity_24hr=', humactivity_24hr
      write (12, *) 'icefrac=', icefrac
      write (12, *) 'id=', id
      write (12, *) 'ie_a=', ie_a
      write (12, *) 'ie_end=', ie_end
      write (12, *) 'ie_m=', ie_m
      write (12, *) 'ie_start=', ie_start
      write (12, *) 'imin=', imin
      write (12, *) 'internalwateruse_h=', internalwateruse_h
      write (12, *) 'irrfracconif=', irrfracconif
      write (12, *) 'irrfracdecid=', irrfracdecid
      write (12, *) 'irrfracgrass=', irrfracgrass
      write (12, *) 'isec=', isec
      write (12, *) 'it=', it
      write (12, *) 'evapmethod=', evapmethod
      write (12, *) 'iy=', iy
      write (12, *) 'kkanohm=', kkanohm
      write (12, *) 'kmax=', kmax
      write (12, *) 'lai_id=', lai_id
      write (12, *) 'laicalcyes=', laicalcyes
      write (12, *) 'laimax=', laimax
      write (12, *) 'laimin=', laimin
      write (12, *) 'lai_obs=', lai_obs
      write (12, *) 'laipower=', laipower
      write (12, *) 'laitype=', laitype
      write (12, *) 'lat=', lat
      write (12, *) 'lenday_id=', lenday_id
      write (12, *) 'ldown_obs=', ldown_obs
      write (12, *) 'lng=', lng
      write (12, *) 'maxconductance=', maxconductance
      write (12, *) 'maxfcmetab=', maxfcmetab
      write (12, *) 'maxqfmetab=', maxqfmetab
      write (12, *) 'snowwater=', snowwater
      ! write (12, *) 'metforcingdata_grid=', metforcingdata_grid
      write (12, *) 'minfcmetab=', minfcmetab
      write (12, *) 'minqfmetab=', minqfmetab
      write (12, *) 'min_res_bioco2=', min_res_bioco2
      write (12, *) 'narp_emis_snow=', narp_emis_snow
      write (12, *) 'narp_trans_site=', narp_trans_site
      write (12, *) 'netradiationmethod=', netradiationmethod
      write (12, *) 'ohm_coef=', ohm_coef
      write (12, *) 'ohmincqf=', ohmincqf
      write (12, *) 'ohm_threshsw=', ohm_threshsw
      write (12, *) 'ohm_threshwd=', ohm_threshwd
      write (12, *) 'pipecapacity=', pipecapacity
      write (12, *) 'popdensdaytime=', popdensdaytime
      write (12, *) 'popdensnighttime=', popdensnighttime
      write (12, *) 'popprof_24hr=', popprof_24hr
      write (12, *) 'pormax_dec=', pormax_dec
      write (12, *) 'pormin_dec=', pormin_dec
      write (12, *) 'precip=', precip
      write (12, *) 'preciplimit=', preciplimit
      write (12, *) 'preciplimitalb=', preciplimitalb
      write (12, *) 'press_hpa=', press_hpa
      write (12, *) 'qf0_beu=', qf0_beu
      write (12, *) 'qf_a=', qf_a
      write (12, *) 'qf_b=', qf_b
      write (12, *) 'qf_c=', qf_c
      write (12, *) 'qn1_obs=', qn1_obs
      write (12, *) 'qh_obs=', qh_obs
      write (12, *) 'qs_obs=', qs_obs
      write (12, *) 'qf_obs=', qf_obs
      write (12, *) 'radmeltfact=', radmeltfact
      write (12, *) 'raincover=', raincover
      write (12, *) 'rainmaxres=', rainmaxres
      write (12, *) 'resp_a=', resp_a
      write (12, *) 'resp_b=', resp_b
      write (12, *) 'roughlenheatmethod=', roughlenheatmethod
      write (12, *) 'roughlenmommethod=', roughlenmommethod
      write (12, *) 'runofftowater=', runofftowater
      write (12, *) 's1=', s1
      write (12, *) 's2=', s2
      write (12, *) 'sathydraulicconduct=', sathydraulicconduct
      write (12, *) 'sddfull=', sddfull
      write (12, *) 'sdd_id=', sdd_id
      write (12, *) 'sfr=', sfr
      write (12, *) 'smdmethod=', smdmethod
      write (12, *) 'snowalb=', snowalb
      write (12, *) 'snowalbmax=', snowalbmax
      write (12, *) 'snowalbmin=', snowalbmin
      write (12, *) 'snowpacklimit=', snowpacklimit
      write (12, *) 'snowdens=', snowdens
      write (12, *) 'snowdensmax=', snowdensmax
      write (12, *) 'snowdensmin=', snowdensmin
      write (12, *) 'snowfallcum=', snowfallcum
      write (12, *) 'snowfrac=', snowfrac
      write (12, *) 'snowlimbldg=', snowlimbldg
      write (12, *) 'snowlimpaved=', snowlimpaved
      write (12, *) 'snowfrac_obs=', snowfrac_obs
      write (12, *) 'snowpack=', snowpack
      write (12, *) 'snowprof_24hr=', snowprof_24hr
      write (12, *) 'snowuse=', snowuse
      write (12, *) 'soildepth=', soildepth
      write (12, *) 'soilstore_id=', soilstore_id
      write (12, *) 'soilstorecap=', soilstorecap
      write (12, *) 'stabilitymethod=', stabilitymethod
      write (12, *) 'startdls=', startdls
      write (12, *) 'state_id=', state_id
      write (12, *) 'statelimit=', statelimit
      write (12, *) 'storageheatmethod=', storageheatmethod
      write (12, *) 'storedrainprm=', storedrainprm
      write (12, *) 'surfacearea=', surfacearea
      write (12, *) 'tair_av=', tair_av
      write (12, *) 'tau_a=', tau_a
      write (12, *) 'tau_f=', tau_f
      write (12, *) 'tau_r=', tau_r
      write (12, *) 'tmax_id=', tmax_id
      write (12, *) 'tmin_id=', tmin_id
      write (12, *) 't_critic_cooling=', t_critic_cooling
      write (12, *) 't_critic_heating=', t_critic_heating
      write (12, *) 'temp_c=', temp_c
      write (12, *) 'tempmeltfact=', tempmeltfact
      write (12, *) 'th=', th
      write (12, *) 'theta_bioco2=', theta_bioco2
      write (12, *) 'timezone=', timezone
      write (12, *) 'tl=', tl
      write (12, *) 'trafficrate=', trafficrate
      write (12, *) 'trafficunits=', trafficunits
      write (12, *) 'traffprof_24hr=', traffprof_24hr
      ! write (12, *) 'ts5mindata_ir=', ts5mindata_ir
      write (12, *) 'tstep=', tstep
      write (12, *) 'tstep_prev=', tstep_prev
      write (12, *) 'veg_type=', veg_type
      write (12, *) 'waterdist=', waterdist
      write (12, *) 'waterusemethod=', waterusemethod
      write (12, *) 'wetthresh=', wetthresh
      write (12, *) 'wu_m3=', wu_m3
      write (12, *) 'wuday_id=', wuday_id
      write (12, *) 'decidcap_id=', decidcap_id
      write (12, *) 'albdectr_id=', albdectr_id
      write (12, *) 'albevetr_id=', albevetr_id
      write (12, *) 'albgrass_id=', albgrass_id
      write (12, *) 'porosity_id=', porosity_id
      write (12, *) 'wuprofa_24hr=', wuprofa_24hr
      write (12, *) 'wuprofm_24hr=', wuprofm_24hr
      write (12, *) 'xsmd=', xsmd
      write (12, *) 'z=', z
      write (12, *) 'z0m_in=', z0m_in
      write (12, *) 'zdm_in=', zdm_in
      write (12, *) '/'

      WRITE (12, *) ''

      CLOSE (12)

   end if

   ! ======================================================================

   RETURN

120 FORMAT(8f10.3, a16)  !format (10g10.2)
121 FORMAT(a12, 24f10.4, a20)

END SUBROUTINE SUEWS_Translate
!===================================================================================

!SUEWS_TranslateBack
!Translates model variables to arrays for each grid
!Runs at the end of SUEWS_Calculations to store correct info for each grid
!Made by HW Nov 2014
!-----------------------------------------------------------------------------------
!Last modified:LJ 14 Sep 2015
!              HCW 28 Nov 2014
!===================================================================================
SUBROUTINE SUEWS_TranslateBack(Gridiv, ir, irMax)

   USE allocateArray
   USE ColNamesInputFiles
   USE ColNamesModelDailyState
   USE data_in
   USE defaultnotUsed
   USE gis_data
   USE Initial
   USE mod_z
   USE resist
   USE snowMod
   USE sues_data
   USE time

   IMPLICIT NONE

   INTEGER::Gridiv, & ! Index of the analysed grid (Gridcounter)
             ir, & ! Meteorological forcing file index (set to zero if SUEWS_Translate called from InitialState)
             irMax        ! Last row in current chunk of met data

   ! =============================================================================
   ! === Translate values from variable names used in model to ModelDailyState ===
   ! =============================================================================

   ! ModelDailyState(Gridiv,cMDS_porosity)    = porosity(id)
   ! ModelDailyState(Gridiv,cMDS_albDecTr)    = albDecTr(id)
   ! ModelDailyState(Gridiv,cMDS_albEveTr)    = albEveTr(id)
   ! ModelDailyState(Gridiv,cMDS_albGrass)    = albGrass(id)
   ! ModelDailyState(Gridiv,cMDS_DecidCap)    = DecidCap(id)
   ModelDailyState(Gridiv, cMDS_SnowfallCum) = SnowfallCum

   ModelDailyState(Gridiv, cMDS_porosity) = porosity_id
   ModelDailyState(Gridiv, cMDS_albDecTr) = albDecTr_id
   ModelDailyState(Gridiv, cMDS_albEveTr) = albEveTr_id
   ModelDailyState(Gridiv, cMDS_albGrass) = albGrass_id
   ModelDailyState(Gridiv, cMDS_DecidCap) = DecidCap_id

   ! Save required DailyState variables for the current grid (HCW 27 Nov 2014)
   ! HDD_grids(:,:,Gridiv)    = HDD(:,:)
   ! GDD_grids(:,:,Gridiv)    = GDD(:,:)
   ! LAI_grids(:,:,Gridiv)    = LAI(:,:)
   ! WUDay_grids(:,:,Gridiv)  = WUDay(:,:)
   ! AlbDecTr_grids(:,Gridiv) = AlbDecTr(:)
   ! AlbEveTr_grids(:,Gridiv) = AlbEveTr(:)
   ! AlbGrass_grids(:,Gridiv) = AlbGrass(:)
   ! DecidCap_grids(:,Gridiv) = DecidCap(:)
   ! Porosity_grids(:,Gridiv) = Porosity(:)

   ! ! update qn1 memory with values of current time step
   ! qn1_store(:,Gridiv)      = qn1_store_grid(:)
   ! qn1_av_store(:,Gridiv)   = qn1_av_store_grid(:)
   ! IF (SnowUse == 1) THEN
   !    qn1_S_store(:,Gridiv)    = qn1_S_store_grid(:)
   !    qn1_S_av_store(:,Gridiv) = qn1_S_av_store_grid(:)
   ! ENDIF

   ! update averaged qn1 memory
   dqndt_grids(Gridiv) = dqndt
   qn1_av_grids(Gridiv) = qn1_av
   tair_av_grids(Gridiv) = tair_av
   IF (SnowUse == 1) THEN
      dqnsdt_grids(Gridiv) = dqnsdt
      qn1_s_av_grids(Gridiv) = qn1_s_av
   ENDIF

   ! added by TS 29 Jun 2018 to remove annual loops in main calculation
   GDD_id_grids(:, Gridiv) = GDD_id
   SDD_id_grids(:, Gridiv) = SDD_id
   Tmin_id_grids(Gridiv) = Tmin_id
   Tmax_id_grids(Gridiv) = Tmax_id
   lenDay_id_grids(Gridiv) = lenDay_id
   HDD_id_grids(:, Gridiv) = HDD_id
   LAI_id_grids(:, Gridiv) = LAI_id
   WUDay_id_grids(:, Gridiv) = WUDay_id

   DecidCap_id_grids(Gridiv) = DecidCap_id
   albDecTr_id_grids(Gridiv) = albDecTr_id
   albEveTr_id_grids(Gridiv) = albEveTr_id
   albGrass_id_grids(Gridiv) = albGrass_id
   porosity_id_grids(Gridiv) = porosity_id

   ! ---- Snow density of each surface
   ModelDailyState(Gridiv, cMDS_SnowDens(1:nsurf)) = SnowDens(1:nsurf)
   ModelDailyState(Gridiv, cMDS_SnowAlb) = SnowAlb

   ! =============================================================================
   ! === Translate values from variable names used in model to ModelOutputData ===
   ! =============================================================================

   ModelOutputData(ir, cMOD_State(1:nsurf), Gridiv) = state_id(1:nsurf)
   ModelOutputData(ir, cMOD_SoilState(1:nsurf), Gridiv) = soilstore_id(1:nsurf)
   ModelOutputData(ir, cMOD_SnowFrac(1:nsurf), Gridiv) = SnowFrac(1:nsurf)
   ModelOutputData(ir, cMOD_SnowPack(1:nsurf), Gridiv) = SnowPack(1:nsurf)
   ModelOutputData(ir, cMOD_SnowWaterState(1:nsurf), Gridiv) = SnowWater(1:nsurf)

   ! ---- ice fraction
   IceFrac_grids(:, Gridiv) = IceFrac

   IF (ir == irMax) THEN   !Store variables ready for next chunk of met data
      ModelOutputData(0, cMOD_State(1:nsurf), Gridiv) = state_id(1:nsurf)
      ModelOutputData(0, cMOD_SoilState(1:nsurf), Gridiv) = soilstore_id(1:nsurf)
      ModelOutputData(0, cMOD_SnowFrac(1:nsurf), Gridiv) = SnowFrac(1:nsurf)
      ModelOutputData(0, cMOD_SnowPack(1:nsurf), Gridiv) = SnowPack(1:nsurf)
      ModelOutputData(0, cMOD_SnowWaterState(1:nsurf), Gridiv) = SnowWater(1:nsurf)
   ENDIF

   RETURN
endsubroutine SUEWS_TranslateBack
!===================================================================================
