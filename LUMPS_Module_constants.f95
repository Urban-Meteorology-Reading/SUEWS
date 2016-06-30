! sg feb 2012 - added some comments
! sg feb 2012 - changed number of surfaces to allocatable array
! lj jun 2012 - snow part added
! HW, LJ Oct 2014 - fixes to the structure
! HCW 03 Mar 2015 - tidied
! HCW 03 Jul 2015 - added albedo max & min to SUEWS_NonVeg.txt (min not used), SUEWS_Veg.txt, SUEWS_Water.txt (min not used)
! LJ 06 Jul 2015 - changed alb_snow, albsnowmin and albsnowmax to SnowAlb, SnowAlbMin and SnowAlbMax (to be systematic with
!                   other variables). Similarly denssnow changed to SnowDens. cMDS_SnowAlb=29 added.
! HCW 10 Mar 2016 - variable vsmd added for soil moisture of vegetated surfaces
! TS 14 Mar 2016 - multiple addtions for AnOHM
! HCW 14 Jun 2015 - updated columns for ESTM and column names for AnOHM

!==================================================================================================
MODULE allocateArray

  IMPLICIT NONE

  INTEGER:: ESTMOne
  
  ! ---- Set parameters for reading in data ------------------------------------------------------
  INTEGER, PARAMETER:: MaxNumberOfGrids=2000   !Max no. grids   !HCW changed to 2000 from 10000 so prog can run on windows (2GB lim)
  INTEGER, PARAMETER:: MaxLinesMet=8640        !Max no. lines to read in one go (for all grids, ie MaxLinesMet/NumberOfGrids each)

  ! ---- Set number of columns in input files ----------------------------------------------------
  INTEGER, PARAMETER:: ncolumnsSiteSelect=97        !SUEWS_SiteSelect.txt
  INTEGER, PARAMETER:: ncolumnsNonVeg=22            !SUEWS_NonVeg.txt
  INTEGER, PARAMETER:: ncolumnsVeg=33               !SUEWS_Veg.txt
  INTEGER, PARAMETER:: ncolumnsWater=19             !SUEWS_Water.txt
  INTEGER, PARAMETER:: ncolumnsSnow=23              !SUEWS_Snow.txt
  INTEGER, PARAMETER:: ncolumnsSoil=9               !SUEWS_Soil.txt
  INTEGER, PARAMETER:: ncolumnsConductance=12       !SUEWS_Conductance.txt
  INTEGER, PARAMETER:: ncolumnsOHMCoefficients=4    !SUEWS_OHMCoefficients.txt
  INTEGER, PARAMETER:: ncolumnsESTMCoefficients=52  !SUEWS_ESTMCoefficients.txt ! S.O. 04 Feb 2016
  INTEGER, PARAMETER:: ncolumnsAnthropogenicHeat=11 !SUEWS_AnthropogenicHeat.txt
  INTEGER, PARAMETER:: ncolumnsIrrigation=25        !SUEWS_Irrigation.txt
  INTEGER, PARAMETER:: ncolumnsProfiles=25          !SUEWS_Profiles.txt
  INTEGER, PARAMETER:: ncolumnsWGWaterDist=10       !SUEWS_WithinGridWaterDist.txt
  INTEGER, PARAMETER:: ncolumnsMetForcingData=24    !Meteorological forcing file (_data.txt)
  INTEGER, PARAMETER:: ncolsESTMdata=13           !ESTM input file (_ESTM_Ts_data.txt))
  
  ! ---- Set number of columns in output files ---------------------------------------------------
  INTEGER, PARAMETER:: ncolumnsDataOut=72,&    !Main output file (_5.txt). DataOut created in SUEWS_Calculations.f95
       ncolumnsDataOutSnow=102

  ! ---- Define input file headers ---------------------------------------------------------------
  CHARACTER(len=20),DIMENSION(ncolumnsSiteSelect)::        HeaderSiteSelect_File          !Header for SiteSelect.txt
  CHARACTER(len=20),DIMENSION(ncolumnsNonVeg)::            HeaderNonVeg_File              !Header for the nonveg surface
  CHARACTER(len=20),DIMENSION(ncolumnsNonVeg)::            HeaderNonVeg_Reqd              !Expected header for the nonveg surface
  CHARACTER(len=20),DIMENSION(ncolumnsVeg)::               HeaderVeg_File                 !Header for the veg surface
  CHARACTER(len=20),DIMENSION(ncolumnsVeg)::               HeaderVeg_Reqd                 !Expected header for the veg surface
  CHARACTER(len=20),DIMENSION(ncolumnsWater)::             HeaderWater_File               !Header for water surface
  CHARACTER(len=20),DIMENSION(ncolumnsWater)::             HeaderWater_Reqd               !Expected header for water surface
  CHARACTER(len=20),DIMENSION(ncolumnsSnow)::              HeaderSnow_File                !Header for Snow surface
  CHARACTER(len=20),DIMENSION(ncolumnsSnow)::              HeaderSnow_Reqd                !Expected header for Snow surface
  CHARACTER(len=20),DIMENSION(ncolumnsSoil)::              HeaderSoil_File                !Header for soils
  CHARACTER(len=20),DIMENSION(ncolumnsSoil)::              HeaderSoil_Reqd                !Expected header for soils
  CHARACTER(len=20),DIMENSION(ncolumnsConductance)::       HeaderCond_File                !Header for conductances
  CHARACTER(len=20),DIMENSION(ncolumnsConductance)::       HeaderCond_Reqd                !Expected header for conductances
  CHARACTER(len=20),DIMENSION(ncolumnsOHMCoefficients)::   HeaderOHMCoefficients_File     !Header for soils
  CHARACTER(len=20),DIMENSION(ncolumnsOHMCoefficients)::   HeaderOHMCoefficients_Reqd     !Expected header for soils
  CHARACTER(len=20),DIMENSION(ncolumnsESTMCoefficients)::  HeaderESTMCoefficients_File    !Header for soils            ! S.O. 04 Feb 2016
  CHARACTER(len=20),DIMENSION(ncolumnsESTMCoefficients)::  HeaderESTMCoefficients_Reqd    !Expected header for soils   ! S.O. 04 Feb 2016
  CHARACTER(len=20),DIMENSION(ncolumnsAnthropogenicHeat):: HeaderAnthropogenicHeat_File   !Header for QF
  CHARACTER(len=20),DIMENSION(ncolumnsAnthropogenicHeat):: HeaderAnthropogenicHeat_Reqd   !Expected header for QF
  CHARACTER(len=20),DIMENSION(ncolumnsIrrigation)::        HeaderIrrigation_File          !Header for Irrigation
  CHARACTER(len=20),DIMENSION(ncolumnsIrrigation)::        HeaderIrrigation_Reqd          !Expected header for Irrigation
  CHARACTER(len=20),DIMENSION(ncolumnsProfiles)::          HeaderProfiles_File            !Header for Profiles
  CHARACTER(len=20),DIMENSION(ncolumnsProfiles)::          HeaderProfiles_Reqd            !Expected header for Profiles
  CHARACTER(len=20),DIMENSION(ncolumnsWGWaterDist)::       HeaderWGWaterDist_File         !Header for Profiles
  CHARACTER(len=20),DIMENSION(ncolumnsWGWaterDist)::       HeaderWGWaterDist_Reqd         !Expected header for Profiles

  ! ---- Define arrays to store input information from SiteInfo spreadsheet ----------------------
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::SiteSelect                !Stores info from SiteSelect.txt
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::NonVeg_Coeff              !Coefficients for the nonveg surfaces
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Veg_Coeff                 !Coefficients for the veg surfaces
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Water_Coeff               !Coefficients for the water surface
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Snow_Coeff                !Coefficients for snow
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Soil_Coeff                !Coefficients for soil
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Conductance_Coeff         !Coefficients for conductances
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::OHMCoefficients_Coeff     !Coefficients for OHMCoefficients
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::ESTMCoefficients_Coeff    !Coefficients for ESTMCoefficients   ! S.O. 04 Feb 2016
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::AnthropogenicHeat_Coeff   !Coefficients for AnthropogenicHeat
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Irrigation_Coeff          !Coefficients for Irrigation
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::Profiles_Coeff            !Coefficients for Profiles
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::WGWaterDist_Coeff         !Coefficients for WithinGridWaterDist

  ! ---- Define arrays for model calculations ----------------------------------------------------
  INTEGER(KIND(1d0)),DIMENSION(:), ALLOCATABLE:: GridIDmatrix         !Array containing GridIDs in SiteSelect
  REAL(KIND(1d0)),DIMENSION(:,:),  ALLOCATABLE:: SurfaceChar          !Array for surface characteristics
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: MetForcingData       !Array for meteorological forcing data
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: ESTMForcingData      !Array for ESTM forcing data
  REAL(KIND(1d0)),DIMENSION(:,:),  ALLOCATABLE:: ModelDailyState      !DailyState array
  REAL(KIND(1d0)),DIMENSION(:),    ALLOCATABLE:: DailyStateFirstOpen
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: ModelOutputData      !Output data matrix
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: dataOut              !Main data output matrix
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: dataOutBL            !CBL output matrix
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: dataOutSOL           !SOLWEIG POI output matrix
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: dataOutSnow          !Main data output matrix
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: dataOutESTM          !ESTM output matrix

  ! ---- Define array for hourly profiles interpolated to tstep ----------------------------------
  REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE:: TstepProfiles
  REAL(KIND(1d0)),DIMENSION(:,:),  ALLOCATABLE:: AHProf_tstep
  REAL(KIND(1d0)),DIMENSION(:,:),  ALLOCATABLE:: WUProfM_tstep, WUProfA_tstep
  
  ! ---- For ESTM
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)::  Ts5mindata     !surface temperature input data
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:) ::   Tair24HR

  ! Column numbers for TstepProfiles
  INTEGER:: cTP_EnUseWD  = 1,&
            cTP_EnUseWE  = 2,&
            cTP_WUManuWD = 3,&
            cTP_WUManuWE = 4,&
            cTP_WUAutoWD = 5,&
            cTP_WUAutoWE = 6,&
            cTP_SnowCWD  = 7,&
            cTP_SnowCWE  = 8

  !-----------------------------------------------------------------------------------------------

  ! ---- Surface types ---------------------------------------------------------------------------
  INTEGER, PARAMETER:: nsurf=7                !Total number of surfaces
  INTEGER, PARAMETER:: NVegSurf=3             !Number of surfaces that are vegetated
  INTEGER, PARAMETER:: nsurfIncSnow=nsurf+1   !Number of surfaces + snow

  INTEGER:: PavSurf   = 1,&   !When all surfaces considered together (1-7)
       BldgSurf  = 2,&
       ConifSurf = 3,&
       DecidSurf = 4,&
       GrassSurf = 5,&   !New surface classes: Grass = 5th/7 surfaces
       BSoilSurf = 6,&   !New surface classes: Bare soil = 6th/7 surfaces
       WaterSurf = 7,&
       ExcessSurf= 8,&   !Runoff or subsurface soil in WGWaterDist
       NSurfDoNotReceiveDrainage=0,&   !Number of surfaces that do not receive drainage water (green roof)
       ivConif = 1,&     !When only vegetated surfaces considered (1-3)
       ivDecid = 2,&
       ivGrass = 3

  REAL(KIND(1d0)),DIMENSION(nsurf):: sfr   !Surface fractions [-]

  ! ---- Water balance for each surface  ---------------------------------------------------------
  !These variables are expressed as depths [mm] over each surface(is); the depth therefore varies with sfr(is)
  REAL(KIND(1d0)),DIMENSION(nsurf):: AddWater       !Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: AddWaterRunoff !Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
  ! N.B. this is not an amount; drain(is)*AddWaterRunoff(is) is the amount [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: chang          !Change in state [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: drain          !Drainage of each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: evap           !Evaporation from each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: runoff         !Runoff from each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: runoffSoil     !Soil runoff from each soil sub-surface [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: smd_nsurf      !Soil moisture deficit of each sub-surface [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: soilmoist      !Soil moisture of each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: soilmoistOld   !Soil moisture of each surface type from previous timestep [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: state          !Wetness status of each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: stateOld       !Wetness status of each surface type from previous timestep [mm]

  REAL(KIND(1d0)),DIMENSION(nsurf):: WetThresh      !When State > WetThresh, rs=0 limit in SUEWS_evap [mm] (specified in input files)
  REAL(KIND(1d0)),DIMENSION(nsurf):: StateLimit     !Limit for state of each surface type [mm] (specified in input files)

  ! ---- Soil characteristics specified in input files -------------------------------------------
  REAL(KIND(1d0)),DIMENSION(nsurf):: SatHydraulicConduct !Saturated hydraulic conductivity for each soil subsurface [mm s-1]
  REAL(KIND(1d0)),DIMENSION(nsurf):: SoilDepth           !Depth of sub-surface soil store for each surface [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf):: SoilStoreCap        !Capacity of soil store for each surface [mm]

  ! ---- Within-grid water distribution matrix ---------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(nsurf+1,nsurf-1)::WaterDist !Within-grid water distribution to other surfaces and runoff/soil store [-]

  ! ---- Drainage characteristics ----------------------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(6,nsurf):: surf   !Storage capacities and drainage equation info for each surface
  ! 1 - min storage capacity [mm]
  ! 2 - Drainage equation to use
  ! 3 - Drainage coeff 1 [units depend on choice of eqn]
  ! 4 - Drainage coeff 2 [units depend on choice of eqn]
  ! 5 - max storage capacity [mm]
  ! 6 - current storage capacity [mm]
  !-----------------------------------------------------------------------------------------------

  ! ---- Define arrays at daily timestep ---------------------------------------------------------
  INTEGER, PARAMETER:: ndays = 366   !Max no. days in a year used to specify size of daily arrays
  !! Could delete NDays and allocate these elsewhere once no. days is known
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5):: GDD          !Growing Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, 6):: HDD          !Heating Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 9):: WU_Day       !Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf):: LAI   !LAI for each veg surface [m2 m-2]

  ! Seasonality of deciduous trees accounted for by the following variables which change with time
  REAL(KIND(1d0)),DIMENSION( 0:ndays):: DecidCap   !Storage capacity of deciduous trees [mm]
  REAL(KIND(1d0)),DIMENSION( 0:ndays):: porosity   !Porosity of deciduous trees [-]

  REAL(KIND(1d0)),DIMENSION( 0:ndays):: albDecTr     !Albedo of deciduous trees [-]
  REAL(KIND(1d0)),DIMENSION( 0:ndays):: albEveTr     !Albedo of evergreen trees [-]
  REAL(KIND(1d0)),DIMENSION( 0:ndays):: albGrass   !Albedo of grass[-]

  REAL(KIND(1d0)):: AlbMin_DecTr,&   !Min albedo for deciduous trees [-]
       AlbMax_DecTr,&   !Max albedo for deciduous trees [-]
       CapMin_dec,&   !Min storage capacity for deciduous trees [mm] (from input information)
       CapMax_dec,&   !Max storage capacity for deciduous trees [mm] (from input information)
       PorMin_dec,&  !Min porosity for deciduous trees
       PorMax_dec,&    !Max porosity for deciduous trees
       AlbMin_EveTr,&   !Min albedo for evergreen trees [-]
       AlbMax_EveTr,&   !Max albedo for evergreen trees [-]
       AlbMin_Grass,&   !Min albedo for grass [-]
       AlbMax_Grass     !Max albedo for grass [-]

  ! Replicate arrays needed for DailyState, adding dimension to identify the grid, HCW 27 Nov 2014
  !! Could delete MaxNumberOfGrids and allocate these elsewhere once NumberOfGrids is known
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5,MaxNumberOfGrids):: GDD_grids
  REAL(KIND(1d0)),DIMENSION(-4:ndays, 6,MaxNumberOfGrids):: HDD_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 9,MaxNumberOfGrids):: WU_Day_grids
  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf,MaxNumberOfGrids):: LAI_grids

  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: albDecTr_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: DecidCap_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: porosity_grids

  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: albEveTr_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: albGrass_grids

  ! AnOHM related: added by TS 01 Mar 2016
  ! store AnOHM coef. of all sfc. by TS 09 Apr 2016
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: Bo_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: mAH_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: a1AnOHM_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: a2AnOHM_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids):: a3AnOHM_grids
  REAL(KIND(1d0)),DIMENSION( 0:ndays,MaxNumberOfGrids,nsurf,3):: a123AnOHM_gs
  !! store water states for AnOHM iteration, by TS 13 Apr 2016
  !REAL(KIND(1d0)),DIMENSION(0:ndays,MaxNumberOfGrids,nsurf):: soilmoistDay   !Soil moisture of each surface type at the end of a day [mm], 13 Apr 2016 TS
  !REAL(KIND(1d0)),DIMENSION(0:ndays,MaxNumberOfGrids,nsurf):: stateDay       !Wetness status of each existing surface type at the end of a day [mm], 13 Apr 2016 TS


  ! Day of week, month and season (used for water use and energy use calculations, and in OHM)
  INTEGER,DIMENSION(0:ndays,3)::DayofWeek   !1 - day of week; 2 - month; 3 - season
  !-----------------------------------------------------------------------------------------------

  ! --- Vegetation phenology ---------------------------------------------------------------------
  ! Parameters provided in input information for each vegetation surface (SUEWS_Veg.txt)
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: BaseT            !Base temperature for growing degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: BaseTe           !Base temperature for senescence degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: GDDFull          !Growing degree days needed for full capacity [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: SDDFull          !Senescence degree days needed to initiate leaf off [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: LaiMin           !Min LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: LaiMax           !Max LAI  [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: MaxConductance   !Max conductance [mm s-1]
  REAL(KIND(1d0)),DIMENSION(4)       :: LaiPower         !Coeffs for LAI equation: 1,2 - leaf growth; 3,4 - leaf off
  !! N.B. currently DecTr only, although input provided for all veg types
  INTEGER:: LAIType                                      !LAI equation to use: original (0) or new (1)
  !real(kind(1d0))::GDDmax,SDDMax                         ! Max GDD and SDD across all veg types [degC] (removed HCW 03 Mar 2015)

  !No longer used (removed HCW 27 Nov 2014)
  !real(kind(1d0)),dimension(0:23)::runT           ! running average T for the day
  !real(kind(1d0)),dimension(0:23)::runP           ! running total Precip for the day
  !real (kind(1d0))::avT_h, totP_h                 ! daily running average Temp, Total precip
  !-----------------------------------------------------------------------------------------------

  ! ---- Variables related to NARP ---------------------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(nsurf):: alb    !Albedo of each surface type [-]
  REAL(KIND(1d0)),DIMENSION(nsurf):: emis   !Emissivity of each surface type [-]

  ! Radiation balance components for different surfaces
  REAL(KIND(1d0)),DIMENSION(nsurf):: Tsurf_ind,&        !Surface temperature for each surface [degC]
       Tsurf_ind_snow,&   !Snow surface temperature for each surface [degC]
       Tsurf_ind_nosnow
  REAL(KIND(1d0)),DIMENSION(nsurf):: kup_ind,&          !Outgoing shortwave radiation for each surface [W m-2]
       kup_ind_snow,&     !Outgoing shortwave radiation for each snow surface [W m-2]
       kup_ind_nosnow
  REAL(KIND(1d0)),DIMENSION(nsurf):: lup_ind,&          !Outgoing longwave radiation for each surface [W m-2]
       lup_ind_snow,&     !Outgoing longwave radiation for each snow surface [W m-2]
       lup_ind_nosnow
  REAL(KIND(1d0)),DIMENSION(nsurf):: qn1_ind,&          !Net all-wave radiation for each surface [W m-2]
       qn1_ind_snow,&     !Net all-wave radiation for each snow surface [W m-2]
       qn1_ind_nosnow

  ! ---- NARP-specific parameters ----------------------------------------------------------------
  REAL(KIND(1d0))             :: NARP_LAT,NARP_LONG,NARP_YEAR,NARP_TZ,&
       NARP_ALB_SNOW,NARP_EMIS_SNOW,NARP_TRANS_SITE
  REAL(KIND(1D0))             :: NARP_G(365)   !!Should this be NDays?? - HCW
  INTEGER                     :: NARP_NPERHOUR
  REAL(KIND(1D0)),ALLOCATABLE :: NARP_KDOWN_HR(:)
  ! Constants required
  REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292,&
       RAD2DEG=57.29577951,&
       SIGMA_SB=5.67E-8
  !-----------------------------------------------------------------------------------------------

  ! ---- OHM coefficients ------------------------------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(9,4,3):: OHM_coef   !Array for OHM coefficients
  REAL(KIND(1d0)):: a1,a2,a3   !OHM coefficients, a1 [-]; a2 [h]; a3 [W m-2]
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids):: q1_grids,q2_grids,q3_grids,&   !For rate of change for OHM for each grid
       r1_grids,r2_grids,r3_grids     !Same for snow
  !-----------------------------------------------------------------------------------------------

  ! ---- Snow-related variables ------------------------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(nsurf):: changSnow,&       !Change in snowpack in mm
       maxSnowVol,&      !! Maximum snow volume
       MeltWaterStore,&  !!Liquid water in the snow pack of ith surface
       ev_snow,&          !!Evaporation from snowpack in mm
       mw_ind,&           !Melt water from individual surface in mm
       mw_indDay,&        !!Melt water per day from each surface type in m3
       runoffSnow,&       !!Runoff from snowpack in mm and in m3
       SnowDens,&        !Density of snow
       SnowFrac,&         !!Surface fraction of snow cover
       iceFrac,&
       snowInit,&
       snowDepth,&       !Depth of snow in cm
       SnowToSurf,&      !Meltwater flowing from snow to surface
       volSWE,&
       StateFraction,&   !Fraction of state that can freeze
       freezMelt,&       !Amount of freezing meltwater in mm for the ith surface area
       Qm_freezState,&   !Heat by freezing of surface state
       freezState,&      !Amount of freezing state in mm for the ith surface area
       FreezStateVol,&
       Qm_melt,&         !Heat consumption by snow melt
       Qm_rain,&         !Heat by rain falling on snow
       rainOnSnow,&      !Liquid precipitation falling on snow ()
       snowD,&
       deltaQi

  REAL(KIND(1d0)),DIMENSION(nsurf):: snowPack,&        !Amount of snow on each surface in mm
       snowPackOld
  INTEGER,DIMENSION(nsurf):: heiG,&                    !snow layer height
       snowCoverForms,&
       snowCalcSwitch=0          !Defines if snow related balance is made
  !-----------------------------------------------------------------------------------------------

  ! ---- Grid connections ------------------------------------------------------------------------
  !! Grid connections needs coding, currently no water transfer between grids
  ! Added HCW 14 Nov 2014
  INTEGER,PARAMETER:: nconns = 8   !Number of grids for between-grid connections
  REAL(KIND(1d0)),DIMENSION(nconns):: GridToFrac   !Fraction of water moving to the grid specified in GridTo [-]
  REAL(KIND(1d0)),DIMENSION(nconns):: GridTo       !Grid that water moves to
  !!character(len=15),dimension(2,MaxNumberOfGrids)::GridConnections  !List of different grid corrections
  !!real (kind(1d0)),dimension(MaxNumberOfGrids)::   GridConnectionsFrac   !Fraction of water moving between the different grids
  !-----------------------------------------------------------------------------------------------

  ! ---- AnOHM related variable, added by TS, 01 Mar 2016 ---------------------------------------------------------------
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids) :: a1AnOHM,a2AnOHM,a3AnOHM   ! OHM coefficients, a1 [-]; a2 [h]; a3 [W m-2]
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids) :: mAHAnOHM           ! daily mean AH [W m-2]
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids) :: BoAnOHMStart  ! initial Bo for interation [-]
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids) :: BoAnOHMEnd  ! final Bo for interation [-]
  REAL(KIND(1d0)),DIMENSION(nsurf):: cpAnOHM      ! heat capacity [??]
  REAL(KIND(1d0)),DIMENSION(nsurf):: kkAnOHM       ! heat conductivity [??]
  REAL(KIND(1d0)),DIMENSION(nsurf):: chAnOHM      ! bulk transfer coef. [-]
  !-----------------------------------------------------------------------------------------------

  ! ESTM variables for SUEWS surfaces
  REAL(KIND(1d0)),DIMENSION(5,nsurfIncSnow):: zSurf_SUEWSsurfs, &
                                              kSurf_SUEWSsurfs, &
                                              rSurf_SUEWSsurfs
  !-----------------------------------------------------------------------------------------------
  !---------------------------------- Column numbers ---------------------------------------------

  ! ---- Set column numbering for SurfaceChar ----------------------------------------------------
  ! Columns 1:80 are the same as in SiteSelect.txt and defined below
  INTEGER:: cc    !Column counter
  INTEGER,PARAMETER:: ccEndSI=ncolumnsSiteSelect

  ! Applicable to each surface
  INTEGER,DIMENSION(nsurf):: c_AlbMin     = (/(cc, cc=ccEndSI+ 0*nsurf+1,ccEndSI+ 0*nsurf+nsurf, 1)/)  !Min. albedo
  INTEGER,DIMENSION(nsurf):: c_AlbMax     = (/(cc, cc=ccEndSI+ 1*nsurf+1,ccEndSI+ 1*nsurf+nsurf, 1)/)  !Max. albedo
  INTEGER,DIMENSION(nsurf):: c_Emis       = (/(cc, cc=ccEndSI+ 2*nsurf+1,ccEndSI+ 2*nsurf+nsurf, 1)/)  !Emissivity
  INTEGER,DIMENSION(nsurf):: c_StorMin    = (/(cc, cc=ccEndSI+ 3*nsurf+1,ccEndSI+ 3*nsurf+nsurf, 1)/)  !Min. storage capacity (canopy)
  INTEGER,DIMENSION(nsurf):: c_StorMax    = (/(cc, cc=ccEndSI+ 4*nsurf+1,ccEndSI+ 4*nsurf+nsurf, 1)/)  !Max. storage capacity (canopy)
  INTEGER,DIMENSION(nsurf):: c_WetThresh  = (/(cc, cc=ccEndSI+ 5*nsurf+1,ccEndSI+ 5*nsurf+nsurf, 1)/)  !Threshold for wet evaporation [mm]
  INTEGER,DIMENSION(nsurf):: c_StateLimit = (/(cc, cc=ccEndSI+ 6*nsurf+1,ccEndSI+ 6*nsurf+nsurf, 1)/)  !Limit for surface state [mm]
  INTEGER,DIMENSION(nsurf):: c_DrEq       = (/(cc, cc=ccEndSI+ 7*nsurf+1,ccEndSI+ 7*nsurf+nsurf, 1)/)  !Drainage equation
  INTEGER,DIMENSION(nsurf):: c_DrCoef1    = (/(cc, cc=ccEndSI+ 8*nsurf+1,ccEndSI+ 8*nsurf+nsurf, 1)/)  !Drainage coef. 1
  INTEGER,DIMENSION(nsurf):: c_DrCoef2    = (/(cc, cc=ccEndSI+ 9*nsurf+1,ccEndSI+ 9*nsurf+nsurf, 1)/)  !Drainage coef. 2
  INTEGER,DIMENSION(nsurf):: c_SoilTCode  = (/(cc, cc=ccEndSI+10*nsurf+1,ccEndSI+10*nsurf+nsurf, 1)/)  !Soil type code

  ! N.B. not included in SUEWS_Water.txt
  INTEGER,DIMENSION(nsurf):: c_SnowLimPat =(/(cc, cc=ccEndSI+11*nsurf+1,ccEndSI+11*nsurf+nsurf, 1)/) !Snow limit for patchiness
  ! N.B. currently only in SUEWS_NonVeg.txt
  INTEGER,DIMENSION(nsurf):: c_SnowLimRem =(/(cc, cc=ccEndSI+12*nsurf+1,ccEndSI+12*nsurf+nsurf, 1)/) !Snow limit for removal
  ! AnOHM TS
  INTEGER,DIMENSION(nsurf):: c_CpAnOHM = (/(cc, cc=ccEndSI+13*nsurf+1,ccEndSI+13*nsurf+nsurf, 1)/) !heat capacity, AnOHM TS
  INTEGER,DIMENSION(nsurf):: c_KkAnOHM = (/(cc, cc=ccEndSI+14*nsurf+1,ccEndSI+14*nsurf+nsurf, 1)/) !heat conductivity, AnOHM TS
  INTEGER,DIMENSION(nsurf):: c_ChAnOHM = (/(cc, cc=ccEndSI+15*nsurf+1,ccEndSI+15*nsurf+nsurf, 1)/) !bulk transfer coef., AnOHM TS
 
  ! Find current column number
  INTEGER,PARAMETER:: ccEndI = (ccEndSI+15*nsurf+nsurf) !add columns for AnOHM, AnOHM TS

  ! Applicable to vegetated surfaces only
  INTEGER,DIMENSION(NVegSurf):: c_BaseT   =(/(cc, cc=ccEndI+ 0*nvegsurf+1,ccEndI+ 0*nvegsurf+nvegsurf, 1)/) !Base temp. for leaf-on
  INTEGER,DIMENSION(NVegSurf):: c_BaseTe  =(/(cc, cc=ccEndI+ 1*nvegsurf+1,ccEndI+ 1*nvegsurf+nvegsurf, 1)/) !Base temp. for leaf-off
  INTEGER,DIMENSION(NVegSurf):: c_GDDFull =(/(cc, cc=ccEndI+ 2*nvegsurf+1,ccEndI+ 2*nvegsurf+nvegsurf, 1)/) !GDD for full LAI
  INTEGER,DIMENSION(NVegSurf):: c_SDDFull =(/(cc, cc=ccEndI+ 3*nvegsurf+1,ccEndI+ 3*nvegsurf+nvegsurf, 1)/) !SDD for start of leaf-fall
  INTEGER,DIMENSION(NVegSurf):: c_LAIMin  =(/(cc, cc=ccEndI+ 4*nvegsurf+1,ccEndI+ 4*nvegsurf+nvegsurf, 1)/) !Min. LAI
  INTEGER,DIMENSION(NVegSurf):: c_LAIMax  =(/(cc, cc=ccEndI+ 5*nvegsurf+1,ccEndI+ 5*nvegsurf+nvegsurf, 1)/) !Max. LAI
  INTEGER,DIMENSION(NVegSurf):: c_GsMax   =(/(cc, cc=ccEndI+ 6*nvegsurf+1,ccEndI+ 6*nvegsurf+nvegsurf, 1)/) !Max. conductance
  INTEGER,DIMENSION(NVegSurf):: c_LAIEq   =(/(cc, cc=ccEndI+ 7*nvegsurf+1,ccEndI+ 7*nvegsurf+nvegsurf, 1)/) !LAI equation
  INTEGER,DIMENSION(NVegSurf):: c_LeafGP1 =(/(cc, cc=ccEndI+ 8*nvegsurf+1,ccEndI+ 8*nvegsurf+nvegsurf, 1)/) !Leaf growth power 1
  INTEGER,DIMENSION(NVegSurf):: c_LeafGP2 =(/(cc, cc=ccEndI+ 9*nvegsurf+1,ccEndI+ 9*nvegsurf+nvegsurf, 1)/) !Leaf growth power 2
  INTEGER,DIMENSION(NVegSurf):: c_LeafOP1 =(/(cc, cc=ccEndI+10*nvegsurf+1,ccEndI+10*nvegsurf+nvegsurf, 1)/) !Leaf-off power 1
  INTEGER,DIMENSION(NVegSurf):: c_LeafOP2 =(/(cc, cc=ccEndI+11*nvegsurf+1,ccEndI+11*nvegsurf+nvegsurf, 1)/) !Leaf-off power 2

  ! Find current column number
  INTEGER,PARAMETER:: ccEndP = (ccEndI+11*nvegsurf+nvegsurf)

  ! Applicable to snow only
  INTEGER:: c_SnowRMFactor = (ccEndP+ 1)
  INTEGER:: c_SnowTMFactor = (ccEndP+ 2)
  INTEGER:: c_SnowAlbMin   = (ccEndP+ 3)
  INTEGER:: c_SnowAlbMax   = (ccEndP+ 4)
  !integer:: c_SnowAlb      = (ccEndP+ 5) (ccEndP free if needed somewhere)
  INTEGER:: c_SnowEmis     = (ccEndP+ 6)
  INTEGER:: c_Snowtau_a    = (ccEndP+ 7)
  INTEGER:: c_Snowtau_f    = (ccEndP+ 8)
  INTEGER:: c_SnowPLimAlb  = (ccEndP+ 9)
  INTEGER:: c_SnowSDMin    = (ccEndP+10)
  INTEGER:: c_SnowSDMax    = (ccEndP+11)
  INTEGER:: c_Snowtau_r    = (ccEndP+12)
  INTEGER:: c_SnowCRWMin   = (ccEndP+13)
  INTEGER:: c_SnowCRWMax   = (ccEndP+14)
  INTEGER:: c_SnowPLimSnow = (ccEndP+15)

  ! Find current column number
  INTEGER,PARAMETER:: ccEndSn = (ccEndP+15)

  ! Soil information
  INTEGER,DIMENSION(nsurf):: c_SoilDepth    = (/(cc, cc=ccEndSn+ 0*nsurf+1,ccEndSn+ 0*nsurf+nsurf, 1)/)  ! Volumetric SM capacity
  INTEGER,DIMENSION(nsurf):: c_SoilStCap = (/(cc, cc=ccEndSn+ 1*nsurf+1,ccEndSn+ 1*nsurf+nsurf, 1)/)  ! Volumetric SM capacity
  INTEGER,DIMENSION(nsurf):: c_KSat         = (/(cc, cc=ccEndSn+ 2*nsurf+1,ccEndSn+ 2*nsurf+nsurf, 1)/)  ! Saturated hydraulic conductivity
  INTEGER,DIMENSION(nsurf):: c_SoilDens     = (/(cc, cc=ccEndSn+ 3*nsurf+1,ccEndSn+ 3*nsurf+nsurf, 1)/)  ! Soil Density
  INTEGER,DIMENSION(nsurf):: c_SoilInfRate  = (/(cc, cc=ccEndSn+ 4*nsurf+1,ccEndSn+ 4*nsurf+nsurf, 1)/)  ! Soil infiltration rate
  INTEGER,DIMENSION(nsurf):: c_ObsSMDepth   = (/(cc, cc=ccEndSn+ 5*nsurf+1,ccEndSn+ 5*nsurf+nsurf, 1)/)  ! Depth of SM obs
  INTEGER,DIMENSION(nsurf):: c_ObsSMMax     = (/(cc, cc=ccEndSn+ 6*nsurf+1,ccEndSn+ 6*nsurf+nsurf, 1)/)  ! Obs maximum SM [kg kg-1 OR m3 m-3]
  INTEGER,DIMENSION(nsurf):: c_ObsSNRFrac   = (/(cc, cc=ccEndSn+ 7*nsurf+1,ccEndSn+ 7*nsurf+nsurf, 1)/)  ! Obs fraction of soil without rocks

  ! Find current column number
  INTEGER,PARAMETER:: ccEndSo = (ccEndSn+ 7*nsurf+nsurf)

  ! Surface conductance
  INTEGER:: c_GsG1  = (ccEndSo+ 1)
  INTEGER:: c_GsG2  = (ccEndSo+ 2)
  INTEGER:: c_GsG3  = (ccEndSo+ 3)
  INTEGER:: c_GsG4  = (ccEndSo+ 4)
  INTEGER:: c_GsG5  = (ccEndSo+ 5)
  INTEGER:: c_GsG6  = (ccEndSo+ 6)
  INTEGER:: c_GsTH  = (ccEndSo+ 7)
  INTEGER:: c_GsTL  = (ccEndSo+ 8)
  INTEGER:: c_GsS1  = (ccEndSo+ 9)
  INTEGER:: c_GsS2  = (ccEndSo+10)
  INTEGER:: c_GsKmax  = (ccEndSo+11)

  ! Find current column number
  INTEGER,PARAMETER:: ccEndGs = (ccEndSo+11)

  ! OHM codes
  INTEGER,DIMENSION(nsurfIncSnow):: c_OHMCode_SWet  =(/(cc, cc=ccEndGs+ 0*nsurfIncSnow+1,&
       ccEndGs+ 0*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (summer wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_OHMCode_SDry  =(/(cc, cc=ccEndGs+ 1*nsurfIncSnow+1,&
       ccEndGs+ 1*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (summer dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_OHMCode_WWet  =(/(cc, cc=ccEndGs+ 2*nsurfIncSnow+1,&
       ccEndGs+ 2*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (winter wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_OHMCode_WDry  =(/(cc, cc=ccEndGs+ 3*nsurfIncSnow+1,&
       ccEndGs+ 3*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (winter dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a1_SWet       =(/(cc, cc=ccEndGs+ 4*nsurfIncSnow+1,&
       ccEndGs+ 4*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (summer wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a2_SWet       =(/(cc, cc=ccEndGs+ 5*nsurfIncSnow+1,&
       ccEndGs+ 5*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (summer wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a3_SWet       =(/(cc, cc=ccEndGs+ 6*nsurfIncSnow+1,&
       ccEndGs+ 6*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (summer wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a1_SDry       =(/(cc, cc=ccEndGs+ 7*nsurfIncSnow+1,&
       ccEndGs+ 7*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (summer dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a2_SDry       =(/(cc, cc=ccEndGs+ 8*nsurfIncSnow+1,&
       ccEndGs+ 8*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (summer dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a3_SDry       =(/(cc, cc=ccEndGs+ 9*nsurfIncSnow+1,&
       ccEndGs+ 9*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (summer dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a1_WWet       =(/(cc, cc=ccEndGs+10*nsurfIncSnow+1,&
       ccEndGs+10*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (winter wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a2_WWet       =(/(cc, cc=ccEndGs+11*nsurfIncSnow+1,&
       ccEndGs+11*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (winter wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a3_WWet       =(/(cc, cc=ccEndGs+12*nsurfIncSnow+1,&
       ccEndGs+12*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (winter wet)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a1_WDry       =(/(cc, cc=ccEndGs+13*nsurfIncSnow+1,&
       ccEndGs+13*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (winter dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a2_WDry       =(/(cc, cc=ccEndGs+14*nsurfIncSnow+1,&
       ccEndGs+14*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (winter dry)
  INTEGER,DIMENSION(nsurfIncSnow):: c_a3_WDry       =(/(cc, cc=ccEndGs+15*nsurfIncSnow+1,&
       ccEndGs+15*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (winter dry)

  ! ESTM code for each surface inclduing snow
  INTEGER,DIMENSION(nsurfIncSnow):: c_ESTMCode      = (/(cc, cc=ccEndGs+16*nsurfIncSnow+1,&
                                                             ccEndGs+16*nsurfIncSnow+nsurfIncSnow, 1)/)  !ESTM code
       
  ! Find current column number
  INTEGER,PARAMETER:: ccEndO = (ccEndGs+16*nsurfIncSnow+nsurfIncSnow)

  ! Anthropogenic heat
  INTEGER :: c_BaseTHDD  = (ccEndO+ 1)
  INTEGER :: c_QF_A1    = (ccEndO+ 2)
  INTEGER :: c_QF_B1     = (ccEndO+ 3)
  INTEGER :: c_QF_C1     = (ccEndO+ 4)
  INTEGER :: c_QF_A2     = (ccEndO+ 5)
  INTEGER :: c_QF_B2    = (ccEndO+ 6)
  INTEGER :: c_QF_C2     = (ccEndO+ 7)
  INTEGER :: c_AHMin     = (ccEndO+ 8)
  INTEGER :: c_AHSlope   = (ccEndO+ 9)
  INTEGER :: c_TCritic   = (ccEndO+10)

  ! Find current column number
  INTEGER,PARAMETER:: ccEndA = (ccEndO+10)

  ! Irrigation
  INTEGER :: c_IeStart    = (ccEndA+ 1)
  INTEGER :: c_IeEnd    = (ccEndA+ 2)
  INTEGER :: c_IntWU    = (ccEndA+ 3)
  INTEGER :: c_Faut    = (ccEndA+ 4)
  INTEGER,DIMENSION(3):: c_Ie_a      = (/(cc, cc=ccEndA+4+ 0*3+1, ccEndA+4 + 0*3+3, 1)/)  ! Automatic irrigation coeffs
  INTEGER,DIMENSION(3):: c_Ie_m      = (/(cc, cc=ccEndA+4+ 1*3+1, ccEndA+4 + 1*3+3, 1)/)  ! Manual irrigation coeffs
  INTEGER,DIMENSION(7):: c_DayWat    = (/(cc, cc=ccEndA+10+ 0*7+1,ccEndA+10+ 0*7+7, 1)/)  ! Irrigation allowed on each day
  INTEGER,DIMENSION(7):: c_DayWatPer = (/(cc, cc=ccEndA+10+ 1*7+1,ccEndA+10+ 1*7+7, 1)/)  ! Fraction properties using irrigation allowed on each day

  ! Find current column number
  INTEGER,PARAMETER:: ccEndIr = (ccEndA+10+ 1*7+7)

  ! Hourly profiles
  INTEGER,DIMENSION(24):: c_HrProfEnUseWD  = (/(cc, cc=ccEndIr+ 0*24+1, ccEndIr+ 0*24+24, 1)/)  ! Energy use, weekdays
  INTEGER,DIMENSION(24):: c_HrProfEnUseWE  = (/(cc, cc=ccEndIr+ 1*24+1, ccEndIr+ 1*24+24, 1)/)  ! Energy use, weekends
  INTEGER,DIMENSION(24):: c_HrProfWUManuWD = (/(cc, cc=ccEndIr+ 2*24+1, ccEndIr+ 2*24+24, 1)/)  ! Water use, manual, weekdays
  INTEGER,DIMENSION(24):: c_HrProfWUManuWE = (/(cc, cc=ccEndIr+ 3*24+1, ccEndIr+ 3*24+24, 1)/)  ! Water use, manual, weekends
  INTEGER,DIMENSION(24):: c_HrProfWUAutoWD = (/(cc, cc=ccEndIr+ 4*24+1, ccEndIr+ 4*24+24, 1)/)  ! Water use, automatic, weekdays
  INTEGER,DIMENSION(24):: c_HrProfWUAutoWE = (/(cc, cc=ccEndIr+ 5*24+1, ccEndIr+ 5*24+24, 1)/)  ! Water use, automatic, weekends
  INTEGER,DIMENSION(24):: c_HrProfSnowCWD  = (/(cc, cc=ccEndIr+ 6*24+1, ccEndIr+ 6*24+24, 1)/)  ! Snow clearing, weekdays
  INTEGER,DIMENSION(24):: c_HrProfSnowCWE  = (/(cc, cc=ccEndIr+ 7*24+1, ccEndIr+ 7*24+24, 1)/)  ! Snow clearing, weekends

  ! Find current column number
  INTEGER,PARAMETER:: ccEndPr = (ccEndIr+ 7*24+24)

  ! Within-grid water distribution (for each surface)
  INTEGER,DIMENSION(nsurf):: c_WGToPaved = (/(cc, cc=ccEndPr+ 0*nsurf+1,ccEndPr+ 0*nsurf+nsurf, 1)/) !Water dist to Paved
  INTEGER,DIMENSION(nsurf):: c_WGToBldgs = (/(cc, cc=ccEndPr+ 1*nsurf+1,ccEndPr+ 1*nsurf+nsurf, 1)/) !Water dist to Bldgs
  INTEGER,DIMENSION(nsurf):: c_WGToEveTr = (/(cc, cc=ccEndPr+ 2*nsurf+1,ccEndPr+ 2*nsurf+nsurf, 1)/) !Water dist to EveTr
  INTEGER,DIMENSION(nsurf):: c_WGToDecTr = (/(cc, cc=ccEndPr+ 3*nsurf+1,ccEndPr+ 3*nsurf+nsurf, 1)/) !Water dist to DecTr
  INTEGER,DIMENSION(nsurf):: c_WGToGrass = (/(cc, cc=ccEndPr+ 4*nsurf+1,ccEndPr+ 4*nsurf+nsurf, 1)/) !Water dist to Grass
  INTEGER,DIMENSION(nsurf):: c_WGToBSoil = (/(cc, cc=ccEndPr+ 5*nsurf+1,ccEndPr+ 5*nsurf+nsurf, 1)/) !Water dist to BSoil
  INTEGER,DIMENSION(nsurf):: c_WGToWater = (/(cc, cc=ccEndPr+ 6*nsurf+1,ccEndPr+ 6*nsurf+nsurf, 1)/) !Water dist to Water
  INTEGER,DIMENSION(nsurf):: c_WGToRunoff    = (/(cc, cc=ccEndPr+ 7*nsurf+1,ccEndPr+ 7*nsurf+nsurf, 1)/) !Water dist to runoff
  INTEGER,DIMENSION(nsurf):: c_WGToSoilStore = (/(cc, cc=ccEndPr+ 8*nsurf+1,ccEndPr+ 8*nsurf+nsurf, 1)/) !Water dist to sub-surface soil

  ! Find current column number
  INTEGER,PARAMETER:: ccEndWG = (ccEndPr+ 8*nsurf+nsurf)

  !ESTM
  ! Roof/surface characteristics for all surfaces including snow
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_thick1  = (/(cc, cc=ccEndWG+ 0*nsurfIncSnow+1,ccEndWG+ 0*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_k1      = (/(cc, cc=ccEndWG+ 1*nsurfIncSnow+1,ccEndWG+ 1*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_rhoCp1  = (/(cc, cc=ccEndWG+ 2*nsurfIncSnow+1,ccEndWG+ 2*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_thick2  = (/(cc, cc=ccEndWG+ 3*nsurfIncSnow+1,ccEndWG+ 3*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_k2      = (/(cc, cc=ccEndWG+ 4*nsurfIncSnow+1,ccEndWG+ 4*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_rhoCp2  = (/(cc, cc=ccEndWG+ 5*nsurfIncSnow+1,ccEndWG+ 5*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_thick3  = (/(cc, cc=ccEndWG+ 6*nsurfIncSnow+1,ccEndWG+ 6*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_k3      = (/(cc, cc=ccEndWG+ 7*nsurfIncSnow+1,ccEndWG+ 7*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_rhoCp3  = (/(cc, cc=ccEndWG+ 8*nsurfIncSnow+1,ccEndWG+ 8*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_thick4  = (/(cc, cc=ccEndWG+ 9*nsurfIncSnow+1,ccEndWG+ 9*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_k4      = (/(cc, cc=ccEndWG+10*nsurfIncSnow+1,ccEndWG+10*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_rhoCp4  = (/(cc, cc=ccEndWG+11*nsurfIncSnow+1,ccEndWG+11*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_thick5  = (/(cc, cc=ccEndWG+12*nsurfIncSnow+1,ccEndWG+12*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_k5      = (/(cc, cc=ccEndWG+13*nsurfIncSnow+1,ccEndWG+13*nsurfIncSnow+nsurfIncSnow, 1)/)
  INTEGER,DIMENSION(nsurfIncSnow):: c_Surf_rhoCp5  = (/(cc, cc=ccEndWG+14*nsurfIncSnow+1,ccEndWG+14*nsurfIncSnow+nsurfIncSnow, 1)/)
  ! Find current column number
  INTEGER,PARAMETER:: ccEndESTMB = (ccEndWG+14*nsurfIncSnow+nsurfIncSnow)
  ! Other ESTM characteristics are for built surfaces only
  INTEGER:: c_Wall_thick1  = (ccEndESTMB+ 1)
  INTEGER:: c_Wall_k1      = (ccEndESTMB+ 2)
  INTEGER:: c_Wall_rhoCp1  = (ccEndESTMB+ 3)
  INTEGER:: c_Wall_thick2  = (ccEndESTMB+ 4)
  INTEGER:: c_Wall_k2      = (ccEndESTMB+ 5)
  INTEGER:: c_Wall_rhoCp2  = (ccEndESTMB+ 6)
  INTEGER:: c_Wall_thick3  = (ccEndESTMB+ 7)
  INTEGER:: c_Wall_k3      = (ccEndESTMB+ 8)
  INTEGER:: c_Wall_rhoCp3  = (ccEndESTMB+ 9)
  INTEGER:: c_Wall_thick4  = (ccEndESTMB+10)
  INTEGER:: c_Wall_k4      = (ccEndESTMB+11)
  INTEGER:: c_Wall_rhoCp4  = (ccEndESTMB+12)
  INTEGER:: c_Wall_thick5  = (ccEndESTMB+13)
  INTEGER:: c_Wall_k5      = (ccEndESTMB+14)
  INTEGER:: c_Wall_rhoCp5  = (ccEndESTMB+15)
  INTEGER:: c_Internal_thick1  = (ccEndESTMB+16)
  INTEGER:: c_Internal_k1      = (ccEndESTMB+17)
  INTEGER:: c_Internal_rhoCp1  = (ccEndESTMB+18)
  INTEGER:: c_Internal_thick2  = (ccEndESTMB+19)
  INTEGER:: c_Internal_k2      = (ccEndESTMB+20)
  INTEGER:: c_Internal_rhoCp2  = (ccEndESTMB+21)
  INTEGER:: c_Internal_thick3  = (ccEndESTMB+22)
  INTEGER:: c_Internal_k3      = (ccEndESTMB+23)
  INTEGER:: c_Internal_rhoCp3  = (ccEndESTMB+24)
  INTEGER:: c_Internal_thick4  = (ccEndESTMB+25)
  INTEGER:: c_Internal_k4      = (ccEndESTMB+26)
  INTEGER:: c_Internal_rhoCp4  = (ccEndESTMB+27)
  INTEGER:: c_Internal_thick5  = (ccEndESTMB+28)
  INTEGER:: c_Internal_k5      = (ccEndESTMB+29)
  INTEGER:: c_Internal_rhoCp5  = (ccEndESTMB+30)  
  INTEGER:: c_nroom      =  (ccEndESTMB+31)  
  INTEGER:: c_alb_ibld   =  (ccEndESTMB+32)  
  INTEGER:: c_em_ibld    =  (ccEndESTMB+33)  
  INTEGER:: c_CH_iwall   =  (ccEndESTMB+34)  
  INTEGER:: c_CH_iroof   =  (ccEndESTMB+35)  
  INTEGER:: c_CH_ibld    =  (ccEndESTMB+36)  
  ! Find current column number
  INTEGER,PARAMETER:: ccEndESTMM = (ccEndESTMB+36)
  ! For Paved surfaces, there are 3 possible ESTM classes (with _Surf characteristics only)
  INTEGER,DIMENSION(3):: c_Surf_thick1_Paved  = (/(cc, cc=ccEndESTMM+ 0*3+1,ccEndESTMM+ 0*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_k1_Paved      = (/(cc, cc=ccEndESTMM+ 1*3+1,ccEndESTMM+ 1*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_rhoCp1_Paved  = (/(cc, cc=ccEndESTMM+ 2*3+1,ccEndESTMM+ 2*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_thick2_Paved  = (/(cc, cc=ccEndESTMM+ 3*3+1,ccEndESTMM+ 3*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_k2_Paved      = (/(cc, cc=ccEndESTMM+ 4*3+1,ccEndESTMM+ 4*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_rhoCp2_Paved  = (/(cc, cc=ccEndESTMM+ 5*3+1,ccEndESTMM+ 5*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_thick3_Paved  = (/(cc, cc=ccEndESTMM+ 6*3+1,ccEndESTMM+ 6*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_k3_Paved      = (/(cc, cc=ccEndESTMM+ 7*3+1,ccEndESTMM+ 7*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_rhoCp3_Paved  = (/(cc, cc=ccEndESTMM+ 8*3+1,ccEndESTMM+ 8*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_thick4_Paved  = (/(cc, cc=ccEndESTMM+ 9*3+1,ccEndESTMM+ 9*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_k4_Paved      = (/(cc, cc=ccEndESTMM+10*3+1,ccEndESTMM+10*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_rhoCp4_Paved  = (/(cc, cc=ccEndESTMM+11*3+1,ccEndESTMM+11*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_thick5_Paved  = (/(cc, cc=ccEndESTMM+12*3+1,ccEndESTMM+12*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_k5_Paved      = (/(cc, cc=ccEndESTMM+13*3+1,ccEndESTMM+13*3+3, 1)/)
  INTEGER,DIMENSION(3):: c_Surf_rhoCp5_Paved  = (/(cc, cc=ccEndESTMM+14*3+1,ccEndESTMM+14*3+3, 1)/)
  ! Find current column number
  INTEGER,PARAMETER:: ccEndESTMMP = (ccEndESTMM+14*3+3)
  ! For Bldgs surfaces, there are 5 possible ESTM classes (all characteristics)
  INTEGER,DIMENSION(5):: c_Surf_thick1_Bldgs  = (/(cc, cc=ccEndESTMMP+ 0*5+1,ccEndESTMMP+ 0*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_k1_Bldgs      = (/(cc, cc=ccEndESTMMP+ 1*5+1,ccEndESTMMP+ 1*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_rhoCp1_Bldgs  = (/(cc, cc=ccEndESTMMP+ 2*5+1,ccEndESTMMP+ 2*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_thick2_Bldgs  = (/(cc, cc=ccEndESTMMP+ 3*5+1,ccEndESTMMP+ 3*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_k2_Bldgs      = (/(cc, cc=ccEndESTMMP+ 4*5+1,ccEndESTMMP+ 4*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_rhoCp2_Bldgs  = (/(cc, cc=ccEndESTMMP+ 5*5+1,ccEndESTMMP+ 5*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_thick3_Bldgs  = (/(cc, cc=ccEndESTMMP+ 6*5+1,ccEndESTMMP+ 6*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_k3_Bldgs      = (/(cc, cc=ccEndESTMMP+ 7*5+1,ccEndESTMMP+ 7*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_rhoCp3_Bldgs  = (/(cc, cc=ccEndESTMMP+ 8*5+1,ccEndESTMMP+ 8*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_thick4_Bldgs  = (/(cc, cc=ccEndESTMMP+ 9*5+1,ccEndESTMMP+ 9*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_k4_Bldgs      = (/(cc, cc=ccEndESTMMP+10*5+1,ccEndESTMMP+10*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_rhoCp4_Bldgs  = (/(cc, cc=ccEndESTMMP+11*5+1,ccEndESTMMP+11*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_thick5_Bldgs  = (/(cc, cc=ccEndESTMMP+12*5+1,ccEndESTMMP+12*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_k5_Bldgs      = (/(cc, cc=ccEndESTMMP+13*5+1,ccEndESTMMP+13*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Surf_rhoCp5_Bldgs  = (/(cc, cc=ccEndESTMMP+14*5+1,ccEndESTMMP+14*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_thick1_Bldgs  = (/(cc, cc=ccEndESTMMP+15*5+1,ccEndESTMMP+15*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_k1_Bldgs      = (/(cc, cc=ccEndESTMMP+16*5+1,ccEndESTMMP+16*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_rhoCp1_Bldgs  = (/(cc, cc=ccEndESTMMP+17*5+1,ccEndESTMMP+17*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_thick2_Bldgs  = (/(cc, cc=ccEndESTMMP+18*5+1,ccEndESTMMP+18*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_k2_Bldgs      = (/(cc, cc=ccEndESTMMP+19*5+1,ccEndESTMMP+19*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_rhoCp2_Bldgs  = (/(cc, cc=ccEndESTMMP+20*5+1,ccEndESTMMP+20*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_thick3_Bldgs  = (/(cc, cc=ccEndESTMMP+21*5+1,ccEndESTMMP+21*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_k3_Bldgs      = (/(cc, cc=ccEndESTMMP+22*5+1,ccEndESTMMP+22*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_rhoCp3_Bldgs  = (/(cc, cc=ccEndESTMMP+23*5+1,ccEndESTMMP+23*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_thick4_Bldgs  = (/(cc, cc=ccEndESTMMP+24*5+1,ccEndESTMMP+24*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_k4_Bldgs      = (/(cc, cc=ccEndESTMMP+25*5+1,ccEndESTMMP+25*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_rhoCp4_Bldgs  = (/(cc, cc=ccEndESTMMP+26*5+1,ccEndESTMMP+26*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_thick5_Bldgs  = (/(cc, cc=ccEndESTMMP+27*5+1,ccEndESTMMP+27*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_k5_Bldgs      = (/(cc, cc=ccEndESTMMP+28*5+1,ccEndESTMMP+28*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Wall_rhoCp5_Bldgs  = (/(cc, cc=ccEndESTMMP+29*5+1,ccEndESTMMP+29*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_thick1_Bldgs  = (/(cc, cc=ccEndESTMMP+30*5+1,ccEndESTMMP+30*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_k1_Bldgs      = (/(cc, cc=ccEndESTMMP+31*5+1,ccEndESTMMP+31*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_rhoCp1_Bldgs  = (/(cc, cc=ccEndESTMMP+32*5+1,ccEndESTMMP+32*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_thick2_Bldgs  = (/(cc, cc=ccEndESTMMP+33*5+1,ccEndESTMMP+33*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_k2_Bldgs      = (/(cc, cc=ccEndESTMMP+34*5+1,ccEndESTMMP+34*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_rhoCp2_Bldgs  = (/(cc, cc=ccEndESTMMP+35*5+1,ccEndESTMMP+35*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_thick3_Bldgs  = (/(cc, cc=ccEndESTMMP+36*5+1,ccEndESTMMP+36*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_k3_Bldgs      = (/(cc, cc=ccEndESTMMP+37*5+1,ccEndESTMMP+37*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_rhoCp3_Bldgs  = (/(cc, cc=ccEndESTMMP+38*5+1,ccEndESTMMP+38*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_thick4_Bldgs  = (/(cc, cc=ccEndESTMMP+39*5+1,ccEndESTMMP+39*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_k4_Bldgs      = (/(cc, cc=ccEndESTMMP+40*5+1,ccEndESTMMP+40*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_rhoCp4_Bldgs  = (/(cc, cc=ccEndESTMMP+41*5+1,ccEndESTMMP+41*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_thick5_Bldgs  = (/(cc, cc=ccEndESTMMP+42*5+1,ccEndESTMMP+42*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_k5_Bldgs      = (/(cc, cc=ccEndESTMMP+43*5+1,ccEndESTMMP+43*5+5, 1)/)
  INTEGER,DIMENSION(5):: c_Internal_rhoCp5_Bldgs  = (/(cc, cc=ccEndESTMMP+44*5+1,ccEndESTMMP+44*5+5, 1)/) 
  INTEGER,DIMENSION(5):: c_nroom_Bldgs      =  (ccEndESTMMP+44*5+5+ 1)  
  INTEGER,DIMENSION(5):: c_alb_ibld_Bldgs   =  (ccEndESTMMP+44*5+5+ 2)  
  INTEGER,DIMENSION(5):: c_em_ibld_Bldgs    =  (ccEndESTMMP+44*5+5+ 3)  
  INTEGER,DIMENSION(5):: c_CH_iwall_Bldgs   =  (ccEndESTMMP+44*5+5+ 4)  
  INTEGER,DIMENSION(5):: c_CH_iroof_Bldgs   =  (ccEndESTMMP+44*5+5+ 5)  
  INTEGER,DIMENSION(5):: c_CH_ibld_Bldgs    =  (ccEndESTMMP+44*5+5+ 6)  
  
  !Last column number for SurfaceChar array
  INTEGER,PARAMETER:: MaxNCols_c = (ccEndESTMMP+44*5+5+ 6)
  !-----------------------------------------------------------------------------------------------

  ! ---- Set column numbering for ModelOutputData ------------------------------------------------
  ! Applicable to each surface
  INTEGER,PARAMETER:: ccMOD = 32
  INTEGER,DIMENSION(nsurf):: cMOD_State          =(/(cc, cc=ccMOD+ 0*nsurf+1,ccMOD+ 0*nsurf+nsurf, 1)/)  !Above ground state
  INTEGER,DIMENSION(nsurf):: cMOD_SoilState      =(/(cc, cc=ccMOD+ 1*nsurf+1,ccMOD+ 1*nsurf+nsurf, 1)/)  !Below ground state (soil store)
  INTEGER,DIMENSION(nsurf):: cMOD_SnowWaterState =(/(cc, cc=ccMOD+ 2*nsurf+1,ccMOD+ 2*nsurf+nsurf, 1)/)  !Liquid (melted) water
  INTEGER,DIMENSION(nsurf):: cMOD_SnowPack       =(/(cc, cc=ccMOD+ 3*nsurf+1,ccMOD+ 3*nsurf+nsurf, 1)/)  !SWE
  INTEGER,DIMENSION(nsurf):: cMOD_SnowFrac       =(/(cc, cc=ccMOD+ 4*nsurf+1,ccMOD+ 4*nsurf+nsurf, 1)/)  !Snow fraction
  INTEGER,DIMENSION(nsurf):: cMOD_SnowDens       =(/(cc, cc=ccMOD+ 5*nsurf+1,ccMOD+ 5*nsurf+nsurf, 1)/)  !Snow density

  !Last column number for ModelOutputData array
  INTEGER,PARAMETER:: MaxNCols_cMOD = ccMOD+ 5*nsurf+nsurf
  !-----------------------------------------------------------------------------------------------

  ! ---- Set column numbering for ModelDailyState ------------------------------------------------
  ! Applicable to each surface
  INTEGER,PARAMETER:: ccMDS = 30
  INTEGER,DIMENSION(nsurf):: cMDS_SnowDens       =(/(cc, cc=ccMDS+ 0*nsurf+1,ccMDS+ 0*nsurf+nsurf, 1)/)  !Snow density

  !Last column number for ModelDailyState array
  INTEGER,PARAMETER:: MaxNCols_cMDS = ccMDS+ 0*nsurf+nsurf
  !-----------------------------------------------------------------------------------------------

  ! ---- Set column numbering for ESTM_Ts_data input file ===-------------------------------------
  ! HCW 15 June 2016
  INTEGER, PARAMETER:: cTs_iy = 1
  INTEGER, PARAMETER:: cTs_id = 2
  INTEGER, PARAMETER:: cTs_it = 3
  INTEGER, PARAMETER:: cTs_imin = 4
  INTEGER, PARAMETER:: cTs_Tiair = 5
  INTEGER, PARAMETER:: cTs_Tsurf = 6
  INTEGER, PARAMETER:: cTs_Troof = 7
  INTEGER, PARAMETER:: cTs_Troad = 8
  INTEGER, PARAMETER:: cTs_Twall = 9
  INTEGER, PARAMETER:: cTs_Twall_n = 10
  INTEGER, PARAMETER:: cTs_Twall_e = 11
  INTEGER, PARAMETER:: cTs_Twall_s = 12
  INTEGER, PARAMETER:: cTs_Twall_w = 13
    
  
END MODULE allocateArray
!==================================================================================================

!==================================================================================================
MODULE Initial

  IMPLICIT NONE

  INTEGER::FirstYear,&          !First year to run (specified in SiteSelect.txt)
       LastYear,&           !Last year to run  (specified in SiteSelect.txt)
       FirstGrid,&    !First grid to run (as in SiteSelect)
       LastGrid,&            !Last grid to run  (as in SiteSelect)
       NumberOfGrids,&      !Number of grids
       GridCounter,&        !Counter for grids (i.e. from 1 to NumberOfGrids)

       ReadBlocksMetData,&  !Number of blocks of met data to read (for each grid, for each year)
       ReadLinesMetData,&   !Number of lines of met data in each block (for each grid)

       nlinesMetData,&            !Number of lines in Met Forcing file
       nlinesESTMdata,&           !Number of lines in ESTM Forcing file
       nlinesSiteSelect,&         !Number of lines in SUEWS_SiteSelect.txt
       nlinesNonVeg,&             !Number of lines in SUEWS_NonVeg.txt
       nlinesVeg,&                !Number of lines in SUEWS_Veg.txt
       nlinesWater,&             !Number of lines in SUEWS_Water.txt
       nlinesSnow,&               !Number of lines in SUEWS_Snow.txt
       nlinesSoil,&               !Number of lines in SUEWS_Soil.txt
       nlinesConductance,&        !Number of lines in SUEWS_Conductance.txt
       nlinesOHMCoefficients,&    !Number of lines in SUEWS_OHMCoefficients.txt
       nlinesESTMCoefficients,&   !Number of lines in SUEWS_ESTMCoefficients.txt
       nlinesAnthropogenicHeat,&  !Number of lines in SUEWS_AnthropogenicHeat.txt
       nlinesIrrigation,&         !Number of lines in SUEWS_Irrigation.txt
       nlinesProfiles,&           !Number of lines in SUEWS_Profiles.txt
       nlinesWGWaterDist,&        !Number of lines in SUEWS_WGWaterDist.txt
       nlines,&                   !Number of lines in different files
       SkippedLines,&           !Number of lines to skip over before reading each block of met data
       iv5            !Counter for code matching.

END MODULE Initial
!==================================================================================================

!==================================================================================================
MODULE data_in

  IMPLICIT NONE

  CHARACTER (len=90)::progname='SUEWS V2016b'  !<<<<<<<<<<<<<<<<<<

  ! ---- Run information ------------------------------------------------------------------------
  CHARACTER (len=20)::  FileCode   !Set in RunControl
  CHARACTER (len=150):: FileInputPath,&   !Filepath for input files (set in RunControl)
       FileOutputPath    !Filepath for output files (set in RunControl)
  ! ---- File names -----------------------------------------------------------------------------
  CHARACTER (len=150):: FileOut,&         !Output file name
       FileChoices,&     !Run characteristics file name
       FileMet,&         !Meteorological forcing file name
       FileDaily,&       !Daily State output file name
       FileESTMTs,&      !ESTM input file name   
       SOLWEIGpoiOut,&   !SOLWEIG poi file name
       BLout             !CLB output file name

  ! ---- Model options set in RunControl --------------------------------------------------------
  INTEGER:: AnthropHeatChoice,&    !QF in met file (0); Loridan et al. 2010 (1); Jarvi et al. 2011 (2)
       CBLuse,&               !CBL slab model used (1) or not used (0)
       MultipleMetFiles,&     !Indicates whether a single met file is used for all grids (0) or one for each grid (1)
       MultipleESTMFiles,&    !Indicates whether a single ESTM input data file is used for all grids (0) or one for each grid (1)
       KeepTstepFilesIn,&     !Delete (0) or keep (1) input met files at resolution of tstep (used by python, not fortran)
       KeepTstepFilesOut,&    !Delete (0) or keep (1) output files at resolution of tstep (used by python, not fortran)
       WriteSurfsFile,&       !Write output file containing variables for each surface (1) or not (0). Not currently used!!
       gsChoice,&             !Options for surface conductance calculation (1 - Ja11, 2 - adjusted method)
       NetRadiationChoice,&   !Options for net all-wave radiation calculation
       OHMIncQF,&             !OHM calculation uses Q* only (0) or Q*+QF (1)
       QSChoice,&             !OHM (1); QS in met file (2); AnOHM(3); ESTM(4)
       SkipHeaderSiteInfo,&   !Number of header lines to skip in SiteInfo files
       SkipHeaderMet,&        !Number of header lines to skip in met file input
       SNOWuse,&              !Snow part used (1) or not used (0)
       SOLWEIGuse,&           !SOLWEIG part used (calculates Tmrt and other fluxes on a grid, FL)
       smd_choice,&           !Use modelled (0) or observed(1,2) soil moisture
       WU_choice,&            !Use modelled (0) or observed (1) water use
       z0_method              !Defines method for calculating z0 & zd

  ! ---- Model options currently set in model, but may be moved to RunControl at a later date
  INTEGER:: AlbedoChoice,&         !No additional albedo varaition (0); zenith angle calculation (1)
                                !Currently set to 0 in SUEWS_Initial
       InputMetFormat,&       !Defines format for met input data: LUMPS format(1) or SUEWS format(10)
                                !Currently set to 10 in SUEWS_Initial
       ity,&                  !Evaporation calculated according to Rutter (1) or Shuttleworth (2)
                                !Currently set to 2 in OverallRunControl
       LAIcalcYes,&           !Use observed (0) or modelled (1) LAI
                                !Currently set to 1 in OverallRunControl
       WriteDailyState        !Daily state file written (1)
  !Currently set to 1 in SUEWS_Initial

  ! ---- Other options used within model --------------------------------------------------------
  INTEGER:: ldown_option           !Parameterisation used for downward longwave radiation (1/2/3)

  ! ---- Output file numbers --------------------------------------------------------------------
  INTEGER:: lfnout,&               !Error Output write units
       lfnoutC,&              !Clean output write units
       lfnOld                 !!Was used for GridConnections

  LOGICAL:: finish,once

  ! ---- Other options set in RunControl --------------------------------------------------------
  REAL (KIND(1d0)):: timezone      !Timezone (GMT=0)

  ! ---- Variables in alphabetical order --------------------------------------------------------
  !! Add units
  REAL (KIND(1d0)):: AH_MIN,&    !Minimum anthropogenic heat flux (AnthropHeatChoice = 1)
       AH_SLOPE,&  !Slope of the antrhropogenic heat flux calculation (AnthropHeatChoice = 1)
       alpha_qhqe,& !Alpha parameter used in LUMPS QH and QE calculations [-]
       avdens,&    !Average air density
       avkdn,&     !Average downwelling shortwave radiation
       avrh,&      !Average relative humidity
       avts,&      !Average surface temperature
       avu1,&      !Average wind speed
       azimuth,&   !Sun azimuth in degrees
       BaseTHDD,&  !Base temperature for QF
       E_mod,&     !Modelled latent heat flux with LUMPS  [W m-2]
       emis_snow,& !Emissivity of snow
       fcld,&      !Cloud fraction modelled
       fcld_obs,&  !Cloud fraction observed
       h_mod,&     !Modelled sensible heat flux with LUMPS [W m-2]
       kclear,&    !Theoretical downward shortwave radiation
       kdiff,&     !Diffuse shortwave radiation
       kdir,&      !Direct shortwave radiation
       kup,&       !Upward shortwave radiation
       lai_obs,&   !LAI for study area provided in met forcing file
       lat,&       !Latitude
       ldown, &    !Downward longwave radiation
       ldown_obs,& !Downwelling longwave radiation
       lng,&       !Longitude
       lup,&       !Upward longwave radiation
       NumCapita,& !Number of people in the study area per hectare [ha-1]
       PopDensDaytime,&   ! Daytime population density [ha-1] (i.e. workers)
       PopDensNighttime,& ! Nighttime population density [ha-1] (i.e. residents)
       Precip,&    !Precipitation per timestep [mm]
       Precip_hr,&    !Precipitation [mm hr-1]
       Press_hPa,&  !Station air pressure in hPa
       Pres_kPa,&   !Station air pressure in kPa
       qe,&        !Observed latent heat flux
       qe_obs,&
       qf,&        !Observed anthropogenic heat flux
       QF_SAHP,&    !Anthropogenic heat flux calculated by SAHP
       qh,&        !Observed sensible heat flux
       qh_obs,&
       qn1,&       !Net all-wave radiation for the study area
       qn1_bup,&
       qn1_obs,&   !Observed new all-wave radiation
       qn1_S,&     !Total net all-wave radiation for the snowpack
       qn1_SF,&    !Total net all-wave radiation for the snowfree surface
       qs,&        !Observed storage heat flux
       QSanOHM,&   !Simulated storage heat flux by AnOHM, TS 30 May 2016
       QSestm,&    !Simulated storage heat flux by ESTM, TS 30 May 2016
       snow,&      !snow cover
       snow_obs,&  !Observed snow cover
       T_CRITIC,& !Critical temperature
       Temp_C,&    !Air temperature
       trans_site,&  !Atmospheric transmittivity
       tsurf,&   !Surface temperature
       wdir,&      ! Wind direction
       wu_m3,&     !Water use provided in met forcing file [m3]
       xsmd,&      !Measured soil moisture deficit
       year,&      !Year of the measurements
       zenith_deg  !Sun zenith angle in degrees

  REAL(KIND(1d0)),DIMENSION(2)::Qf_A,Qf_B,Qf_C   !Qf coefficients
  REAL(KIND(1d0)),DIMENSION(0:23,2):: AHPROF     !Anthropogenic heat profiles for (1)weekdays / (2)weekends

  INTEGER,DIMENSION(2)::DayLightSavingDay   !DOY when daylight saving changes

  INTEGER::nCBLstep  !number of time steps of Runge-kutta methods in one hour

  !---------Water bucket (see B. Offerle's PhD)----------------------------------
  REAL (KIND(1D0)):: DRAINRT,&      !Drainage rate of the water bucket [mm hr-1]
       RAINBUCKET,&   !RAINFALL RESERVOIR [mm]
       RAINCOVER,&
       RAINMAXRES,&   !Maximum water bucket reservoir [mm]
       RAINRES,&      ! [mm]
       TEMPVEG        !TEMPORARY VEGETATIVE SURFACE FRACTION ADJUSTED BY RAINFALL

  !---------SOLWEIG variables---------------------------------------------------
  REAL(KIND(1D0)):: absL,&             ! Absorption coefficient of longwave radiation of a person
       absK,&             ! Absorption coefficient of shortwave radiation of a person
       heightgravity,&    ! Centre of gravity for a standing person
       TransMin,&         ! Tranmissivity of K through decidious vegetation (leaf on)
       TransMax           ! Tranmissivity of K through decidious vegetation (leaf off)

  INTEGER:: Posture,&                ! 1.Standing, 2.Sitting
       usevegdem,&          ! With vegetation (1)
       row,&                    ! Y coordinate for point of interest
       col,&                    ! X coordinate for point of interest
       onlyglobal,&             ! if no diffuse and direct SW, then =1
       SOLWEIGpoi_out,&         ! write output variables at point of interest
       Tmrt_out,&               ! write output Tmrt grid
       Lup2d_out,&              ! write output Lup grid
       Ldown2d_out,&            ! write output Ldown grid
       Kup2d_out,&              ! write output Kup grid
       Kdown2d_out,&            ! write output Kdown grid
       GVF_out,&                ! write output GroundViewFActor grid
       SOLWEIG_ldown,&          ! 1= use SOLWEIG code to estimate Ldown, 0=use SEUWS
       OutInterval,&            ! Output interval in minutes
       RunForGrid               ! If only one grid should be run. All grids -999

  CHARACTER (len=150):: DSMPath,&    ! Path to DSMs
       DSMname,&    ! Ground and building DSM
       CDSMname,&   ! Canopy DSM
       TDSMname,&   ! Trunk zone DSM
       SVFPath,&    ! Path to SVFs
       SVFsuffix,&  !
       buildingsname! Boolean matrix for locations of building pixels

  !--------- AnOHM related variables----------------------------------
  ! to be added here

END MODULE data_in
!==================================================================================================

!======================================================================================================
MODULE cbl_MODULE

  INTEGER::EntrainmentType,&  ! Entrainment type choice
       CO2_included,&     ! CO2 included
       InitialData_use,&  ! 1 read initial data, 0 do not
       qh_choice,&        ! selection of qh use to drive CBL growth 1=Suews 2=lumps 3=obs
       sondeflag      ! 1 read sonde or vertical profile data in 0 do not

  INTEGER,DIMENSION(366)::cblday=0

  CHARACTER (len=200), DIMENSION(366)::FileSonde=""
  CHARACTER (len=200)::InitialDataFileName
  REAL(KIND(1D0)):: wsb       ! subsidence velocity
  REAL(KIND(1d0)),DIMENSION(1:10):: cbldata
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::IniCBLdata

  !Parameters in CBL code
  INTEGER::zmax,&
       nEqn=4,&
       iCBLcount,&
       nlineInData
  REAL(KIND(1d0))::C2K=273.16

  REAL (KIND(1D0)):: usbl,ftbl,fqbl,fcbl,gamt,gamq,gamc,tpp,qpp,cp0!,tk

  REAL(KIND(1D0))::alpha3,&
       blh_m,&    ! Boundary layer height(m)
       blh1_m,&
       cm,&       ! CO2 concentration in CBL
                                !cp0,gamc,& !
       gamt_Km,&  ! Vertical gradient of theta (K/m)
       gamq_gkgm,&! Vertical gradient of specific humidity (g/kg/m)
       gamq_kgkgm,&! Vertical gradient of specific humidity (kg/kg/m)
                                !fcbl,&
       tm_C,&     ! Potential temperature in CBL (degree Celsius)
       tm_K,&     ! Potential temperature in CBL (K)
       tmp_K,&
       tp_C,&     ! Potential temperature just above Boundary layer height(degree Celsius)
       tp_K,&     ! Potential temperature just above Boundary layer height(K)
       tpp_K,&
       febl_kgkgms,&! Kinematic latent heat flux((kg/kg)*m/s)
       fhbl_Kms,&   ! Kinematic sensible heat flux(K*m/s)
       qm_gkg,&   ! Specific humidity in CBL(g/kg)
       qm_kgkg,&  ! Specific humidity in CBL(kg/kg)
       qp_gkg,&   ! Specific humidity above Boundary layer height(g/kg)
       qp_kgkg,&  ! Specific humidity above Boundary layer height(kg/kg)
       qpp_kgkg


  REAL (KIND(1D0)), DIMENSION (0:500,2):: gtheta,ghum ! Vertical gradient of theta and specific humidity from sonde data
  REAL (KIND(1D0)), DIMENSION(4)::y

END   MODULE cbl_MODULE
!===================================================================================

MODULE snowMod
  IMPLICIT NONE

  REAL (KIND(1D0))::AdjMeltFact,&    !Factor between melt and freezing factors
       CumSnowfall,&        !Cumulative snowfall
       fwh,&              !Weighted freezing water
       lvS_J_kg,&         !Latent heat of sublimation in J/kg
       mwh,&              !Weighted hourly water melt
       MwStore,&              !Meltwater storage
       PrecipLimit,&      !Temperature limit when precipitation occurs as snow
       PrecipLimitAlb,&   !Precipitation limit for albedo change (in mm)
       Qm,&               !Snow melt associated heat flux
       QmFreez,&          !Energy released in freezing of meltwater or surface state
       QmRain,&
       qn1_snow,&         !Net all-wave radiation of snowpack
       qn1_nosnow,&       !Same for the snow free surface
       RadMeltFact,&      !Radiation melt factor
       SnowAlb,&          !Snow albedo
       SnowAlbMin,&       !Minimum snow albedo
       SnowAlbMax,&       !Maximum snow albedo
       SnowDensMin,&      !Minimum density of snow
       SnowDensMax,&      !Maximum density of snow
       SnowLimBuild,&     !Snow removal limits for roofs in mm)
       SnowLimPaved,&     !Snow removal limits for paved surfaces in mm)
       swe,&        !Weighted snow water equivalent (in mm)
       tau_a,&            !Time constans related to albedo change
       tau_f,&
       tau_r,&            !Time constant for density increase.
       TempMeltFact,&     !Temperature melt factor
       volDay,&           !Volume of the melted water per day
       zf,&
       WaterHoldCapFrac,& !Water holding capacity factor
       CRWmin,& !Free water holding capacity of deep snowpack
       CRWmax  !Free water holding capacity of shallow snowpack

  REAL(KIND(1D0)), DIMENSION(2)::  SnowRemoval=0 ! Removal of snow in mm
  REAL(KIND(1d0)), DIMENSION(0:23,2):: snowProf  ! Timing of snow removal (0 or 1) Hourly, WD/WE

  INTEGER::SnowFractionChoice   !Choice how fraction of snow is calculated

END MODULE snowMod
!===================================================================================

!==================================================================================================
MODULE defaultNotUsed
  IMPLICIT NONE
  REAL (KIND(1d0)):: notUsed=-55.55,reall,NAN=-999,pNAN=999
  INTEGER:: notUsedI=-55, ios_out,errorChoice  !errorChoice defines if the problemfile is opened for the first time
END MODULE defaultNotUsed
!==================================================================================================

!==================================================================================================
MODULE time
  INTEGER:: iy,&            !Year
       id,&            !Day of year
       it,&            !Hour
       imin,&          !Minutes
       iostat_var,&      !File status from reading data (should not be here)
       DLS                            !day lightsavings =1 + 1h) =0

  REAL(KIND(1d0)):: dectime        !Decimal time
  REAL(KIND(1d0)):: halftimestep   !In decimal time based on interval
  REAL (KIND(1d0)):: tstepcount    !Count number of timesteps in this day
  INTEGER:: nofDaysThisYear        !Based on whether leap year or not

  INTEGER:: iy_prev_t, id_prev_t   !Value of iy and id at previous timestep

END MODULE time
!==================================================================================================

!===================================================================================
MODULE mod_grav
  REAL (KIND(1d0)):: grav=9.80665  !g - gravity - physics today august 1987
END MODULE mod_grav

!===================================================================================
MODULE mod_k
  REAL(KIND(1d0)) :: k=0.4,&             !Von Karman's contant
       k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000 !Limit for neutral stability
END MODULE mod_k

!===================================================================================
MODULE Thresh
  REAL(KIND(1d0)) :: IPThreshold_mmhr = 10   !Threshold for intense precipitation [mm hr-1]

END MODULE Thresh


!===================================================================================
MODULE gas
  !   press (mb) ea (mb)
  IMPLICIT NONE
  REAL (KIND(1d0))::  comp=0.9995
  REAL (KIND(1d0))::  epsil=0.62197   !ratio molecular weight of water vapor/dry air (kg/mol/kg/mol)
  REAL (KIND(1d0))::  epsil_gkg=621.97   !ratio molecular weight of water vapor/dry air in g/kg
  REAL (KIND(1d0))::  dry_gas=8.31451 !Dry gas constant (J/k/mol)
  REAL (KIND(1d0))::  gas_ct_wat=461.05 !Gas constant for water (J/kg/K)
  REAL (KIND(1d0))::  molar=0.028965 !Dry air molar fraction in kg/mol
  REAL (KIND(1d0))::  molar_wat_vap=0.0180153 !Molar fraction of water vapor in kg/mol
  REAL (KIND(1d0))::  gas_ct_dry=8.31451/0.028965 !j/kg/k=dry_gas/molar
  REAL (KIND(1d0))::  gas_ct_wv=8.31451/0.0180153 !j/kg/kdry_gas/molar_wat_vap
END MODULE gas

!**********************************************
MODULE mod_z
  REAL (KIND(1d0)) :: zzd,&  !Active measurement height (meas. height-displac. height)
       z0m,&  !Aerodynamic roughness length
       zdm,&  !Displacement height
       z      !Windspeed height
  REAL(KIND(1E10))::z0V      !Roughness length for vapour
END MODULE mod_z

!**********************************************
MODULE resist  !Variables related surface resistance calculations (P. 1744 in G&O1991)
  IMPLICIT NONE
  REAL (KIND(1d0)):: th,&             !Maximum temperature limit
       tl,&             !Minimum temperature limit
       Kmax,&           !Annual maximum hourly solar radiation
       g1,g2,g3,g4,&    !Fitted parameters related to
       g5,g6,s1,s2,&    !surface res. calculations
       tc,&             !Temperature parameter 1
       tc2              !Temperature parameter 2
END MODULE resist

!**********************************************
MODULE moist
  IMPLICIT NONE

  REAL (KIND(1d0))::avcp,&        !Specific heat capacity
       dens_dry,&    !Dry air density kg m-3
       dq,&          !Specific humidity deficit
       Ea_hPa,&      !Water vapour pressure in hPa
       Es_hPa,&      !Saturation vapour pressure in hPa
       lv_J_kg,&     !Latent heat of vaporization in [J kg-1]
       tlv,&         !Latent heat of vaporization per timestep
                                ![J kg-1 s-1] (tlv=lv_J_kg/tstep_real)
       psyc_hPa,&    !Psychometric constant in hPa
       psycIce_hPa,& !Psychometric constant in hPa for snow
       s_Pa,&        !Vapour pressure versus temperature slope in Pa
       s_hpa,&       !Vapour pressure versus temperature slope in hPa
       sIce_hpa,&    !Vapour pressure versus temperature slope in hPa above ice/snow
       vpd_hPa,&     !Vapour pressure deficit in hPa
       vpd_pa,&      !Vapour pressure deficit in Pa
       waterDens=999.8395 !Density of water in 0 cel deg

END MODULE moist
!**********************************************

MODULE gis_data
  IMPLICIT NONE

  REAL(KIND(1d0)):: Alt,&                        !Altitude in m
       areaunir,&                   !Unirrigated area
       areair,&                     !Irrigated area
       bldgH,&                      !Mean building height
       FAIbldg,&                    !Frontal area fraction of buildings
       FAItree,&                    !Frontal area fraction of trees
       FAIEveTree,&                    !Frontal area fraction of evergreen trees
       FAIDecTree,&                    !Frontal area fraction of deciduous trees
       grassfractionirrigated,&     !Irrigated grass fraction for LUMPS
       pavedfractionirrigated,&     !Irrigated paved area fraction for LUMPS
       TreeH,&                      !Mean tree height
       EveTreeH,&                     !Height of evergreen trees
       DecTreeH,&                     !Height of deciduous trees
       treefractionirrigated,&      !Irrigated tree fraction for LUMPS
       veg_fr,&                     !Vegetation fraction from land area
                                !- For LUMPS - dependent on user choice    & water
       VegFraction, &               ! sum of vegetation -not including water
       ImpervFraction,&             ! sum of surface cover fractions for impervious surfaces
       PervFraction,&               ! sum of surface cover fractions for pervious surfaces
       NonWaterFraction,&           ! sum of surface cover fractions for all except water surfaces
       areaZh                       !=(sfr(BldgSurf)+sfr(ConifSurf)+sfr(DecidSurf)) !Total area of buildings and trees

  INTEGER:: idgis,&      !Time integers used in the code
       itgis,&      !
       Veg_type=1    !Defines how vegetation is calculated for LUMPS

END MODULE gis_data

!************************************************************
MODULE sues_data
  IMPLICIT NONE

  INTEGER:: tstep,&    !Timestep [s] at which the model is run (set in RunControl)
       nsh,&      !Number of timesteps per hour
       t_interval   !Number of seconds in an hour [s] (now set in OverallRunControl)

  REAL(KIND(1d0)):: nsh_real,&   !nsh cast as a real for use in calculations
       tstep_real   !tstep cast as a real for use in calculations

  !Options for model setup (switches, etc) mainly set in RunControl
  INTEGER:: StabilityMethod,&   !Defines stability functions used (set in RunControl)
       RoughLen_heat     !Defines method for calculating roughness length for heat (set in RunControl)


  INTEGER:: in
  INTEGER:: is      !Integer to count over surface types

  !These are variables which currently have been removed from SuesInput.nml
  INTEGER::AerodynamicResistanceMethod=2 !The method used to calculate aerodynamic resistance

  INTEGER::Ie_start,&   !Starting time of water use (DOY)
       Ie_end       !Ending time of water use (DOY)

  REAL(KIND(1d0)),DIMENSION(2):: SurplusEvap !Surplus for evaporation in 5 min timestep
  ! sg -- need to determine size

  !Variables listed in SuesInput.nml
  REAL (KIND(1d0))::FlowChange,&         !Difference between the input and output flow in the water body
       PipeCapacity,&       !Capacity of pipes to transfer water
       RunoffToWater,&      !Fraction of surface runoff going to water body
       SmCap,&              !Volumetric/gravimetric soil moisture capacity
       SoilDensity,&        !Bulk density of soil
       SoilDepthMeas,&      !Soil depth of the measured soil moisture
       SoilRocks,&          !Fraction of rocks in soil
       SurfaceArea,&        !Surface area of the study area [m2]
       SurfaceArea_ha,&     !Surface area of the study area [ha]
       WaterBodyType,&      !If water body type is pond/lake (=1) or river (=2)
       WaterState,&         !State of the water body
       WaterStorCap,&       !Capacity of water body when surface is wet
       WUAreaEveTr_m2,&     !Water use area (evergreen trees) [m2]
       WUAreaDecTr_m2,&     !Water use area (deciduous trees) [m2]
       WUAreaGrass_m2,&     !Water use area (grass) [m2]
       WUAreaTotal_m2,&     !Water use area (total) [m2]
       wu_EveTr,&              !Water use for evergreen trees/shrubs [mm]
       wu_DecTr,&              !Water use for deciduous trees/shrubs [mm]
       wu_Grass                !Water use for grass [mm]

  !Other related to SUES
  REAL (KIND(1d0))::AdditionalWater,&     !Water flow from other grids
       ch_per_interval,&     !Change in state per interval
       chSnow_per_interval,& !Change in snow state per interval
       dI_dt,&               !Water flow between two stores
       dr_per_interval,&     !Drainage per interval
       ev_per_interval,&     !Evaporation per interval
       surf_chang_per_tstep,& !Change in surface state per timestep [mm] (for whole surface)
       tot_chang_per_tstep,&  !Change in surface and soilstate per timestep [mm] (for whole surface)
       NWstate_per_tstep,&     !State per timestep [mm] (for whole surface, excluding water body)
       state_per_tstep,&     !State per timestep [mm] (for whole surface)
       drain_per_tstep,&     !Drainage per timestep [mm] (for whole surface, excluding water body)
       runoff_per_tstep,&     !Runoff per timestep [mm] (for whole surface)
       runoffSoil_per_tstep,& !Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)
       ev_per_tstep,&     !Evaporation per timestep [mm] (for whole surface)
       qe_per_tstep,&     !QE [W m-2] (for whole surface)
       p_mm,&                !Inputs to surface water balance
       pin,&                 !Rain per time interval
       planF,&               !Areally weighted frontal area fraction
       rb,&                  !Boundary layer resistance
                                ! Water leaving each grid for grid-to-grid connectivity
       runoffAGimpervious,&     !Above ground runoff from impervious surface [mm] for whole surface area
       runoffAGveg,&            !Above ground runoff from vegetated surfaces [mm] for whole surface area
       runoffWaterBody,&        !Above ground runoff from water surface [mm] for whole surface area
       runoffPipes,&            !Runoff in pipes [mm] for whole surface area
       runoffAGimpervious_m3,&  !Volume of above ground runoff from impervious surface [m3]
       runoffAGveg_m3,&         !Volume of above ground runoff from vegetated surfaces [m3]
       runoffWaterBody_m3,&     !Volume of above ground runoff from water surface [m3]
       runoffPipes_m3,&         !Volume of runoff in pipes [m3]
       runoff_per_interval,&
                                ! Total water transported to each grid for grid-to-grid connectivity
       addImpervious,&      !Water from impervious surfaces of other grids [mm] for whole surface area
       addVeg,&             !Water from vegetated surfaces of other grids [mm] for whole surface area
       addWaterbody,&       !Water from water surface of other grids [mm] for whole surface area
       addPipes,&           !Water in pipes from other grids [mm] for whole surface area

       runoffSoil_per_interval,&
       qe_per_interval,&     !latent heat per interval
       soilmoistcap,&
       soilstate,&        !Area-averaged soil moisture [mm] for whole surface
       st_per_interval,&!Surface state per interval
       surplusWaterBody,&  !Extra runoff that goes to water body [mm] as specified by RunoffToWater
       tlv_sub,&
       overuse=0,&
       Zh               !Areally weighted roughness element height

  !Calculation of u*,stability and aerodynamic resistance
  REAL (KIND(1d0))::H,&          !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
       l_mod,&      !Monin-Obukhov length (either measured or modelled)
       psim,&       !Stability function of momentum
       psyh,&       !Stability function of heat
       RA,&         !Aerodynamic resistance
       RAsnow,&     !Aerodynamic resistance over snow
       tstar,&      !T*
       ustar,&      !Friction velocity
       z0_gis       !Roughness length for momentum from gis input file

  !Surface resistance related variables
  REAL (KIND(1d0))::resistsurf,& !Surface resistance
       gdq,&        !G(dq)
       qnm,&        !QMAX/(QMAX+G2)
       gq,&         !G(Q*)
       gtemp,&      !G(T)
       gl,&         !G(LAI)
       sdp,&        !S1/G6+S2
       smd,&        !Soil moisture deficit of the soil surface layer
       vsmd,&       !Soil moisture deficit for vegetated surfaces only (what about BSoil?)
       gs,&         !G(Soil moisture deficit)
       gsc          !Surface Layer Conductance

  !SUES latent heat flux related variables
  REAL (KIND(1d0))::  vdrc,&     !Second term up in calculation of E
       numPM,&    !Numerator of PM equation
       sp,&       !Term in calculation of E
       sae,&      !Same
       ev,&       !Evaporation
       rst,&      !Flag in SUEWS_Evap (gets set to 1 if surface dry; 0 if surface wet)
       qeph,&       !Latent heat flux (W m^-2)
       qeOut      !Latent heat flux [W m-2]


  !Water use related variables
  REAL (KIND(1d0)):: ext_wu,&         !External water use for the model timestep [mm] (over whole study area)
       Faut,&           !Fraction of irrigated area using automatic irrigation
       int_wu,&         !Internal water use for the model timestep [mm] (over whole study area)
       IrrFracConif,&  !Fraction of evergreen trees which are irrigated
       IrrFracDecid,&  !Fraction of deciduous trees which are irrigated
       IrrFracGrass,&  !Fraction of grass which is irrigated
       InternalWaterUse_h !Internal water use [mm h-1]

  ! 7 - number of days in week
  REAL(KIND(1d0)),DIMENSION(7)::DayWatPer,&  !% of houses following daily water
       DayWat       !Days of watering allowed
  REAL(KIND(1d0)),DIMENSION(0:23,2):: WUProfM,&   !Hourly profiles for water use (manual irrigation)
       WUProfA   !Hourly profiles for water use (automatic irrigation)


  REAL (KIND(1d0)),DIMENSION(3)::Ie_a,Ie_m   !Coefficients for automatic and manual irrigation models

END MODULE sues_data

!**********************************************
!===================================================================================
MODULE VegPhenogy
  IMPLICIT NONE
  REAL (KIND(1d0)):: VegPhenLumps
END MODULE VegPhenogy

MODULE filename
  CHARACTER (len=90)::  smithfile     !file for narp
END MODULE filename


MODULE InitialCond

  REAL (KIND(1d0))::LAIinitialEveTr,&
       LAIinitialDecTr,&
       LAIinitialGrass,&
       porosity0,&
       DecidCap0,&
       albDecTr0,&
       albEveTr0,&
       albGrass0,&
       Temp_C0,&
       GDD_1_0,&
       GDD_2_0,&
       SoilStorePavedState,&
       SoilStoreBldgsState,&
       SoilStoreEveTrState,&
       SoilStoreDecTrstate,&
       SoilStoreGrassState,&
       SoilStoreBSoilState,&
       SnowWaterPavedState,&
       SnowWaterBldgsState,&
       SnowWaterEveTrState,&
       SnowWaterDecTrState,&
       SnowWaterGrassState,&
       SnowWaterBSoilState,&
       SnowWaterWaterstate,&
       SnowPackPaved,&
       SnowPackBldgs,&
       SnowPackEveTr,&
       SnowPackDecTr,&
       SnowPackGrass,&
       SnowPackBSoil,&
       SnowPackWater,&
       SnowAlb0

  INTEGER::ID_Prev

END MODULE InitialCond

!-------------------------------------------------
!New modules for the column numbers


!-------------------------------------------------------------------------
MODULE ColNamesModelDailyState

  IMPLICIT NONE

  !========== Columns for ModelDailyState array =========================

  INTEGER:: cMDS_id_prev    = 3, &
       cMDS_HDD1            = 4, &
       cMDS_HDD2            = 5, &
       cMDS_TempC           = 6, &
       cMDS_TempCRM         = 7, &
       cMDS_Precip          = 8, &
       cMDS_DaysSinceRain   = 9, &
       cMDS_TempCOld1       = 10,&
       cMDS_TempCOld2       = 11,&
       cMDS_TempCOld3       = 12,&
       cMDS_GDDMin          = 13,&
       cMDS_GDDMax          = 14,&
       cMDS_GDD1_0          = 15,&
       cMDS_GDD2_0          = 16,&
       cMDS_LAIInitialEveTr = 17,&
       cMDS_LAIInitialDecTr = 18,&
       cMDS_LAIInitialGrass = 19,&
       cMDS_porosity        = 20,&
       cMDS_albEveTr        = 21,&
       cMDS_albDecTr        = 22,&
       cMDS_albGrass        = 23,&
       cMDS_DecidCap        = 24,&
       cMDS_CumSnowfall     = 25,&
       cMDS_LAIEveTr        = 26,&
       cMDS_LAIDecTr        = 27,&
       cMDS_LAIGrass        = 28,&
       cMDS_SnowAlb         = 29,&
       cMDS_BoRatio         = 30,& ! noontime Bowen ratio, added by TS
       cMDS_a1AnOHM         = 31,& ! a1 of AnOHM, added by TS
       cMDS_a2AnOHM         = 32,& ! a2 of AnOHM, added by TS
       cMDS_a3AnOHM         = 33 ! a3 of AnOHM, added by TS


END MODULE ColNamesModelDailyState


!-------------------------------------------------------------------------
MODULE ColNamesInputFiles

  IMPLICIT NONE
  
  INTEGER:: ccc    !Column counter

  ! Column names and numbers must match the input files

  !========== Columns for SUEWS_SiteSelect.txt ==========================
  ! Columns 1:97 are the same for SurfaceChar
  INTEGER::c_Grid     = 1,&
       c_Year     = 2,&
       c_StartDLS = 3,&
       c_EndDLS   = 4,&
                                ! Site info
       c_lat  = 5,&
       c_lng  = 6,&
       c_Area = 7,&
       c_Alt  = 8,&
                                ! Time info
       c_id   = 9,&
       c_it   = 10,&
       c_imin = 11,&
                                ! Surface fractions
       c_FrPaved = 12,&
       c_FrBldgs = 13,&
       c_FrEveTr = 14,&
       c_FrDecTr = 15,&
       c_FrGrass = 16,&
       c_FrBSoil = 17,&
       c_FrWater = 18,&
                                ! Irrigated fractions
       c_IrrEveTrFrac = 19,&
       c_IrrDecTrFrac = 20,&
       c_IrrGrassFrac = 21,&
                                ! Height information
       c_HBldgs   = 22,&
       c_HEveTr   = 23,&
       c_HDecTr   = 24,&
       c_z0m      = 25,&
       c_zdm      = 26,&
       c_FAIBldgs = 27,&
       c_FAIEveTr = 28,&
       c_FAIDecTr = 29,&
                                ! Population
       c_PopDensDay   = 30,&
       c_PopDensNight = 31,&
                                ! Codes for different surfaces
       c_PavedCode = 32,&  ! Links characteristics in SUEWS_NonVeg.txt
       c_BldgsCode = 33,&  ! Links characteristics in SUEWS_NonVeg.txt
       c_EveTrCode = 34,&  ! Links characteristics in SUEWS_Veg.txt
       c_DecTrCode = 35,&    ! Links characteristics in SUEWS_Veg.txt
       c_GrassCode = 36,&     ! Links characteristics in SUEWS_Veg.txt
       c_BSoilCode = 37,&  ! Links characteristics in SUEWS_Veg.txt
       c_WaterCode = 38,&       ! Links characteristics in SUEWS_Water.txt
                                ! LUMPS info
       c_LUMPSDr     = 39,&
       c_LUMPSCover  = 40,&
       c_LUMPSMaxRes = 41,&
                                ! NARP info
       c_NARPTrans      =42,&
                                ! Code for conductances
       c_CondCode      =43,&       ! Links characteristics in SUEWS_Conductance.txt
                                ! Code for snow
       c_SnowCode      =44,&    ! Links characteristics in SUEWS_Snow.txt
                                ! Codes for human impacts on energy, water and snow
       c_SnowProfWD= 45,&  ! Snow-clearing profile in SUEWS_Profile.txt (weekdays)
       c_SnowProfWE   = 46,&  ! Snow-clearing profile in SUEWS_Profile.txt (weekends)
       c_QFCode       = 47,&       ! Links anthropogenic heat info in SUEWS_AnthropogenicHeat.txt
       c_EnProfWD     = 48,&  ! Links to energy-use profile in SUEWS_Profile.txt (weekdays)
       c_EnProfWE     = 49,&  ! Links to energy-use profile in SUEWS_Profile.txt (weekends)
       c_IrrCode      = 50,&       ! Links irrigation info in SUEWS_Irrigation.txt
       c_WProfManuWD  = 51,&  ! Links to water-use profile in SUEWS_Profile.txt (manual irrigation, weekdays)
       c_WProfManuWE  = 52,&  ! Links to water-use profile in SUEWS_Profile.txt (manual irrigation, weekends)
       c_WProfAutoWD  = 53,&  ! Links to water-use profile in SUEWS_Profile.txt (automatic irrigation, weekdays)
       c_WProfAutoWE  = 54,&  ! Links to water-use profile in SUEWS_Profile.txt (automatic irrigation, weekends)
                                ! Flow information
       c_FlowChange    =55,&  ! Difference in input & output flows for water surface
       c_RunoffToWater =56,&    ! Fraction of above-ground runoff flowing to water surface
       c_PipeCapacity  =57,&  ! Pipe capacity [mm]
                                ! Runoff (to 8 adjacent grids)
       c_GridConnection1of8 = 58,&
       c_Fraction1of8       = 59,&
       c_GridConnection2of8 = 60,&
       c_Fraction2of8       = 61,&
       c_GridConnection3of8 = 62,&
       c_Fraction3of8       = 63,&
       c_GridConnection4of8 = 64,&
       c_Fraction4of8       = 65,&
       c_GridConnection5of8 = 66,&
       c_Fraction5of8       = 67,&
       c_GridConnection6of8 = 68,&
       c_Fraction6of8       = 69,&
       c_GridConnection7of8 = 70,&
       c_Fraction7of8       = 71,&
       c_GridConnection8of8 = 72,&
       c_Fraction8of8       = 73,&
                                ! Runoff within grid (for each surface type)
       c_WGPavedCode = 74,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGBldgsCode = 75,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGEveTrCode = 76,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGDecTrCode = 77,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGGrassCode = 78,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGBSoilCode = 79,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
       c_WGWaterCode = 80,&   ! Links to SUEWS_WaterDistibuteWithinGrid.txt
                                 ! Additional info for ESTM
       c_AreaWall = 81   ! Wall surface fraction (Awall/Agridcell)
       
       INTEGER,DIMENSION(3) :: c_Fr_ESTMClass_Paved = (/(ccc,ccc=82,84,1)/) ! Fraction of Paved surface with ESTM Class 1-5
       INTEGER,DIMENSION(3) :: c_Code_ESTMClass_Paved = (/(ccc,ccc=85,87,1)/) ! Code for Paved surface ESTM Class 1-5
       INTEGER,DIMENSION(5) :: c_Fr_ESTMClass_Bldgs = (/(ccc,ccc=88,92,1)/) ! Fraction of Bldgs surface with ESTM Class 1-5
       INTEGER,DIMENSION(5) :: c_Code_ESTMClass_Bldgs = (/(ccc,ccc=93,97,1)/) ! Code for Bldgs surface ESTM Class 1-5
       
  !========== Columns for SUEWS_NonVeg.txt ==========================
  INTEGER :: ci_Code   = 1, &
       ci_AlbMin       = 2, &
       ci_AlbMax       = 3, &
       ci_Emis         = 4, &
       ci_StorMin      = 5, &
       ci_StorMax      = 6, &
       ci_WetThresh    = 7, &
       ci_StateLimit   = 8, &
       ci_DrEq         = 9, &
       ci_DrCoef1      = 10,&
       ci_DrCoef2      = 11,&
       ci_SoilTCode    = 12,&
       ci_SnowLimPat   = 13,&
       ci_SnowLimRem   = 14,&
       ci_OHMCode_SWet = 15,&
       ci_OHMCode_SDry = 16,&
       ci_OHMCode_WWet = 17,&
       ci_OHMCode_WDry = 18,&
       ci_ESTMCode     = 19,& ! ESTM code for each surface (if 0 use codes in SiteSelect instead)
       ci_CpAnOHM      = 20,& ! heat capacity, added by TS AnOHM
       ci_KkAnOHM      = 21,& ! heat conductivity, added by TS AnOHM
       ci_ChAnOHM      = 22 ! bulk transfer coef., added by TS AnOHM

  !========== Columns for SUEWS_Veg.txt ============================
  INTEGER :: cp_Code   = 1, &
       cp_AlbMin       = 2, &
       cp_AlbMax       = 3, &
       cp_Emis         = 4, &
       cp_StorMin      = 5, &
       cp_StorMax      = 6, &
       cp_WetThresh    = 7, &
       cp_StateLimit   = 8, &
       cp_DrEq         = 9, &
       cp_DrCoef1      = 10,&
       cp_DrCoef2      = 11,&
       cp_SoilTCode    = 12,&
       cp_SnowLimPat   = 13,&
       cp_BaseT        = 14,&
       cp_BaseTe       = 15,&
       cp_GDDFull      = 16,&
       cp_SDDFull      = 17,&
       cp_LAIMin       = 18,&
       cp_LAIMax       = 19,&
       cp_GsMax        = 20,&
       cp_LAIEq        = 21,&
       cp_LeafGP1      = 22,&
       cp_LeafGP2      = 23,&
       cp_LeafOP1      = 24,&
       cp_LeafOP2      = 25,&
       cp_OHMCode_SWet = 26,&
       cp_OHMCode_SDry = 27,&
       cp_OHMCode_WWet = 28,&
       cp_OHMCode_WDry = 29,&
       cp_ESTMCode     = 30,& 
       cp_CpAnOHM      = 31,& ! heat capacity, added by TS AnOHM
       cp_KkAnOHM      = 32,& ! heat conductivity, added by TS AnOHM
       cp_ChAnOHM      = 33 ! bulk transfer coef., added by TS AnOHM


  !========== Columns for SUEWS_Water.txt ===============================
  INTEGER :: cw_Code   = 1, &
       cw_AlbMin       = 2, &
       cw_AlbMax       = 3, &
       cw_Emis         = 4, &
       cw_StorMin      = 5, &
       cw_StorMax      = 6, &
       cw_WetThresh    = 7, &
       cw_StateLimit   = 8, &
       cw_DrEq         = 9, &
       cw_DrCoef1      = 10,&
       cw_DrCoef2      = 11,&
       cw_OHMCode_SWet = 12,&
       cw_OHMCode_SDry = 13,&
       cw_OHMCode_WWet = 14,&
       cw_OHMCode_WDry = 15,&
       cw_ESTMCode     = 16,&
       cw_CpAnOHM      = 17,& ! heat capacity, added by TS AnOHM
       cw_KkAnOHM      = 18,& ! heat conductivity, added by TS AnOHM
       cw_ChAnOHM      = 19 ! bulk transfer coef., added by TS AnOHM

  !========== Columns for SUEWS_Snow.txt ================================
  INTEGER :: cs_Code   = 1, &
       cs_SnowRMFactor = 2, &
       cs_SnowTMFactor = 3, &
       cs_SnowAlbMin   = 4, &
       cs_SnowAlbMax   = 5, &
       cs_SnowEmis     = 6, &
       cs_Snowtau_a    = 7, &
       cs_Snowtau_f    = 8, &
       cs_SnowPLimAlb  = 9, &
       cs_SnowSDMin    = 10,&
       cs_SnowSDMax    = 11,&
       cs_Snowtau_r    = 12,&
       cs_SnowCRWMin   = 13,&
       cs_SnowCRWMax   = 14,&
       cs_SnowPLimSnow = 15,&
       cs_OHMCode_SWet = 16,&
       cs_OHMCode_SDry = 17,&
       cs_OHMCode_WWet = 18,&
       cs_OHMCode_WDry = 19,&
       cs_ESTMCode     = 20,&
       cs_CpAnOHM      = 21,& ! heat capacity, added by TS
       cs_KkAnOHM      = 22,& ! heat conductivity, added by TS
       cs_ChAnOHM      = 23 ! bulk transfer coef., added by TS

  !========== Columns for SUEWS_Soil.txt ================================
  INTEGER :: cSo_Code        =  1,&
       cSo_SoilDepth   =  2,&
       cSo_SoilStCap   =  3,&
       cSo_KSat         =  4,&
       cSo_SoilDens    =  5,&
       cSo_SoilInfRate =  6,&
       cSo_ObsSMDepth  =  7,&
       cSo_ObsSMMax    =  8,&
       cSo_ObsSNRFrac  =  9

  !========== Columns for SUEWS_Conductance.txt =========================
  INTEGER :: cc_Code   =  1,&
       cc_GsG1   =  2,&
       cc_GsG2   =  3,&
       cc_GsG3   =  4,&
       cc_GsG4   =  5,&
       cc_GsG5   =  6,&
       cc_GsG6   =  7,&
       cc_GsTH   =  8,&
       cc_GsTL   =  9,&
       cc_GsS1   =  10,&
       cc_GsS2   =  11,&
       cc_GsKmax =  12

  !========== Columns for SUEWS_OHMCoefficients.txt =====================
  INTEGER :: cO_Code = 1,&
       cO_a1  = 2,&
       cO_a2  = 3,&
       cO_a3  = 4

  !========== Columns for SUEWS_ESTMCoefficients.txt =====================!   ! S.O. 04 Feb 2016
  INTEGER::cE_Code    = 1, &
       cE_Surf_thick1 = 2, &  !Characteristics for 5x roof/surface layers
       cE_Surf_k1     = 3, &
       cE_Surf_rhoCp1 = 4, &
       cE_Surf_thick2 = 5, &
       cE_Surf_k2     = 6, &
       cE_Surf_rhoCp2 = 7, &
       cE_Surf_thick3 = 8, &
       cE_Surf_k3     = 9, &
       cE_Surf_rhoCp3 = 10, &
       cE_Surf_thick4 = 11, &
       cE_Surf_k4     = 12, &
       cE_Surf_rhoCp4 = 13, &
       cE_Surf_thick5 = 14, &
       cE_Surf_k5     = 15, &
       cE_Surf_rhoCp5 = 16, &
       cE_Wall_thick1 = 17, &   ! Characteristics for 5x external wall layers (used for Bldgs surfaces only)
       cE_Wall_k1     = 18, &
       cE_Wall_rhoCp1 = 19, &
       cE_Wall_thick2 = 20, &
       cE_Wall_k2     = 21, &
       cE_Wall_rhoCp2 = 22, &
       cE_Wall_thick3 = 23, &
       cE_Wall_k3     = 24, &
       cE_Wall_rhoCp3 = 25, &
       cE_Wall_thick4 = 26, &
       cE_Wall_k4     = 27, &
       cE_Wall_rhoCp4 = 28, &
       cE_Wall_thick5 = 29, &
       cE_Wall_k5     = 30, &
       cE_Wall_rhoCp5 = 31, &
       cE_Internal_thick1 = 32, &   ! Characteristics for 5x internal wall layers (used for Bldgs surfaces only)
       cE_Internal_k1     = 33, &
       cE_Internal_rhoCp1 = 34, &
       cE_Internal_thick2 = 35, &
       cE_Internal_k2     = 36, &
       cE_Internal_rhoCp2 = 37, &
       cE_Internal_thick3 = 38, &
       cE_Internal_k3     = 39, &
       cE_Internal_rhoCp3 = 40, &
       cE_Internal_thick4 = 41, &
       cE_Internal_k4     = 42, &
       cE_Internal_rhoCp4 = 43, &
       cE_Internal_thick5 = 44, &
       cE_Internal_k5     = 45, &
       cE_Internal_rhoCp5 = 46, &
       cE_nroom    = 47,&
       cE_alb_ibld = 48,&
       cE_em_ibld  = 49,&
       cE_CH_iwall = 50,&
       cE_CH_iroof = 51,&
       cE_CH_ibld  = 52

  !========== Columns for SUEWS_AnthropogenicHeat.txt ===================
  INTEGER ::   cA_Code    = 1,&
       cA_BaseTHDD = 2,&
       cA_QF_A1    = 3,&   !Weekday
       cA_QF_B1    = 4,&   !Weekday
       cA_QF_C1    = 5,&   !Weekday
       cA_QF_A2    = 6,&   !Weekend
       cA_QF_B2    = 7,&   !Weekend
       cA_QF_C2    = 8,&   !Weekend
       cA_AHMin    = 9,&
       cA_AHSlope  = 10,&
       cA_TCritic  = 11

  !========== Columns for SUEWS_Irrigation.txt ==========================

  INTEGER ::  cIr_Code     = 1,&
       cIr_IeStart     = 2,&
       cIr_IeEnd       = 3,&
       cIr_IntWU       = 4,&
       cIr_Faut        = 5,&
       cIr_Ie_a1       = 6,&
       cIr_Ie_a2       = 7,&
       cIr_Ie_a3       = 8,&
       cIr_Ie_m1       = 9,&
       cIr_Ie_m2       = 10,&
       cIr_Ie_m3       = 11,&
       cIr_DayWat1     = 12,&
       cIr_DayWat2     = 13,&
       cIr_DayWat3     = 14,&
       cIr_DayWat4     = 15,&
       cIr_DayWat5     = 16,&
       cIr_DayWat6     = 17,&
       cIr_DayWat7     = 18,&
       cIr_DayWatPer1  = 19,&
       cIr_DayWatPer2  = 20,&
       cIr_DayWatPer3  = 21,&
       cIr_DayWatPer4  = 22,&
       cIr_DayWatPer5  = 23,&
       cIr_DayWatPer6  = 24,&
       cIr_DayWatPer7  = 25

  !========== Columns for SUEWS_Profile.txt =============================

  INTEGER:: cc   !Column counter

  INTEGER:: cPr_Code = 1
  INTEGER,DIMENSION(24):: cPr_Hours = (/(cc, cc=2,25, 1)/)  ! Hourly profile data

  !========== Columns for SUEWS_WithinGridWaterDist.txt =================

  INTEGER::   cWG_Code     = 1,&
       cWG_ToPaved     = 2,&
       cWG_ToBldgs     = 3,&
       cWG_ToEveTr     = 4,&
       cWG_ToDecTr     = 5,&
       cWG_ToGrass     = 6,&
       cWG_ToBSoil     = 7,&
       cWG_ToWater     = 8,&
       cWG_ToRunoff    = 9,&
       cWG_ToSoilStore = 10

END MODULE ColNamesInputFiles

!----------------------------------------------------------------------------------------
MODULE ESTM_data !S.O. and FO

  ! =======ESTMinput.nml==============================
  INTEGER ::        evolveTibld,&
       TsurfChoice,&
       ibldCHmod

  REAL(KIND(1d0)):: LBC_soil,&        !Lowest boundary condition in soil
       THEAT_ON,&
       THEAT_OFF,&
       THEAT_fix,&
       ivf_iw,&    !Internal view factors : im
       ivf_ir,&
       ivf_ii,&
       ivf_if,&
       ivf_ww,&    !Internal view factors : wall
       ivf_wr,&
       ivf_wi,&
       ivf_wf,&
       ivf_rw,&    !Internal view factors : roof
       ivf_ri,&
       ivf_rf,&
       ivf_fw,&    !Internal view factors : floor
       ivf_fr,&
       ivf_fi


  !=======ESTMcoefficients.txt=================================
  INTEGER:: Nibld,&           !Number of layers in an internal element in buildings, calculated when the file is read.
       Nwall,&           !Number of layers in external wall
       Nroof,&           !Number of layers in roof
       Nground           !Number of layers in ground

  REAL(KIND(1d0)),DIMENSION(5):: zibld,&    !Thickness of layers in internal building
       zwall,&    !Thickness of layers in external wall
       zroof,&    !Thickness of layers in roof
       zground,&  !Thickness of layers in ground
       kibld,&    !Thermal conductivity of layers in internal building
       kwall,&    !Thermal conductivity of layers in external wall
       kroof,&    !Thermal conductivity of layers in roof
       kground,&  !Thermal conductivity of layers in ground
       ribld,&    !Volumetric heat capacity of layers in internal building
       rwall,&    !Volumetric heat capacity of layers in external wall
       rroof,&    !Volumetric heat capacity of layers in roof
       rground    !Volumetric heat capacity of layers in ground
  
  ! Paved and Bldgs surfaces can include 3 and 5 classes respectively
  !For the 3x Paved surfaces
  REAL(KIND(1d0)),DIMENSION(5,3):: zSurf_Paved  
  REAL(KIND(1d0)),DIMENSION(5,3):: kSurf_Paved
  REAL(KIND(1d0)),DIMENSION(5,3):: rSurf_Paved
  !For the 5x Bldgs surfaces
  REAL(KIND(1d0)),DIMENSION(5,5):: zSurf_Bldgs  
  REAL(KIND(1d0)),DIMENSION(5,5):: kSurf_Bldgs
  REAL(KIND(1d0)),DIMENSION(5,5):: rSurf_Bldgs
  REAL(KIND(1d0)),DIMENSION(5,5):: zwall_Bldgs  
  REAL(KIND(1d0)),DIMENSION(5,5):: kwall_Bldgs
  REAL(KIND(1d0)),DIMENSION(5,5):: rwall_Bldgs
  REAL(KIND(1d0)),DIMENSION(5,5):: zibld_Bldgs  
  REAL(KIND(1d0)),DIMENSION(5,5):: kibld_Bldgs
  REAL(KIND(1d0)),DIMENSION(5,5):: ribld_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: nroom_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: alb_ibld_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: em_ibld_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: CH_iwall_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: CH_iroof_Bldgs
  REAL(KIND(1d0)),DIMENSION(5):: CH_ibld_Bldgs
  
  REAL(KIND(1d0))::   nroom,&      !Number of rooms in internal building  (changed from integer to real HCW 16 Jun 2016)
       alb_ibld,& !albedo value of internal elements
       em_ibld,&  !emissivity of internal elements
       CH_iroof,& !bulk transfer coefficient of internal roof
       CH_iwall,& !bulk transfer coefficient of internal wall
       CH_ibld,&  !bulk transfer coefficient of internal building element
       fwall,&    !fraction of wall
       AreaWall   ! Area of wall within grid [m2]

  !=======ESTM Ts input=================================
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:)  ::  Tibld,Twall,Troof,Tground
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)::  Tw_4

  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  ::  Tibld_grids,Twall_grids,Troof_grids,Tground_grids
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:,:)::  Tw_4_grids
  
  !=======variables and parameters created in ESTM=============================
  REAL(KIND(1d0))                           ::  alb_avg,&
       alb_ground,&   !albedo value of ground
       alb_roof,&     !albedo value of roof
       alb_veg,&      !albedo value of veg
       CHAIR,&
       CHR,&
       em_ground,&    !emissivity of ground
       em_roof,&      !emissivity of roof
       em_veg,&       !emissivity of veg
       em_r,&         !emissivity of roof inside building
       em_w,&         !emissivity of internal wall
       em_i,&         !emissivity of ???
       em_f,&         !emissivity of floor
       fair,&         !fraction of air (or ratio of outdoor air volume to indoor air volume)
       fground,&      !fraction of ground
       fibld,&        !fraction of internal elements (?)
       finternal,&    !sum of froof, fibld and fwall
       froof,&        !fraction of roof
       fveg,&         !fraction of veg
       HW,&           !Height Width ratio
       LUP_ground,&
       LUP_ROOF,&
       LUP_VEG,&
       LUP_WALL,&
       minshc_airbld,&
       Pcoeff(5),&
       Qsground,&     !Storage heat flux into ground
       Qsroof,&       !Storage heat flux into roof
       Qswall,&       !Storage heat flux into wall
       Qs_4(4),&      !Storage heat flux into each external wall (N,E,S and W direction)
       Qsair,&        !Storage heat flux into air
       Qsibld,&       !Storage heat flux into internal building elements
       RVF_ground,&
       RVF_WALL,&
       RVF_ROOF,&
       RVF_CANYON,&
       RVF_VEG,&
       SHC_air,&
       SVF_ground,&   !Sky view factor from ground
       SVF_wall,&     !Sky view factor from wall
       SVF_roof,&     !Sky view factor from roof
       TANZENITH,&    !
       Tair1,&
       Tair2,&
       Tairday,&      !24hour average air temperature
       Tfloor,&
       Tievolve,&
       TN_roof,&
       TN_wall,&
       T0_wall,&
       T0_roof,&
       T0_ground,&
       T0_ibld,&
       WS,&           !Wind speed used in ESTM
       xvf_wall,&
       ZREF,&         !local scale reference height
       zvf_ground,&   !wall view factor from ground
       zvf_WALL     !wall view factor from ground

  ! Arrays to store variables for each grid
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: Tair2_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: lup_ground_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: lup_wall_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: lup_roof_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: Tievolve_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: T0_wall_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: T0_roof_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: T0_ground_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: T0_ibld_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: TN_roof_grids
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: TN_wall_grids
  
  ! Surface fractions for ESTM classes   
  REAL(KIND(1d0)),DIMENSION(3):: ESTMsfr_Paved     
  REAL(KIND(1d0)),DIMENSION(5):: ESTMsfr_Bldgs
  
  LOGICAL             ::bctype(2),&
       CFLfail=.FALSE.,&
       diagnoseTi=.FALSE.,&
       first,&
       HVAC=.FALSE.,&
       SPINDONE=.FALSE.

  REAL(KIND(1d0)),PARAMETER::alb_wall=0.23,em_wall=0.9  ! used only when radforce = T but radforce is always set to F.
  INTEGER, PARAMETER::        maxiter=100
  REAL(KIND(1d0)),PARAMETER:: conv=0.0001

  !=============variables maybe will be removed=======================================
  INTEGER             ::nalb,&
       nemis
  REAL(KIND(1d0))     ::sumalb,&
       sumemis

END MODULE ESTM_data

!----------------------------------------------------------------------------------
MODULE MathConstants

  REAL (KIND(1d0)),PARAMETER ::pi=3.14159265359
  REAL (KIND(1d0)),PARAMETER ::dtr=0.0174532925, rtd=57.2957795

END MODULE MathConstants

!----------------------------------------------------------------------------------
MODULE PhysConstants

  REAL (KIND(1d0)),PARAMETER :: C2K = 273.15   !Celsius to Kelvin
  REAL (KIND(1d0)),PARAMETER :: SBConst = 5.67051e-8   !Stefan Boltzmann constant [W m-2 K-4]

END MODULE PhysConstants

!C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95
