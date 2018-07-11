!========================================================================================
! a mini version of SUEWS to be coupled with WRF
! TS 22 Apr 2018: initial
! TS 11 Jun 2018: modified according to recent SUEWS development


MODULE SuMin_Module
  USE SUEWS_Driver,ONLY:SUEWS_cal_Main,&
       PavSurf,BldgSurf,ConifSurf,DecidSurf,GrassSurf,BSoilSurf,WaterSurf,&
       ivConif,ivDecid,ivGrass,&
       ncolumnsDataOutSUEWS,ncolumnsDataOutSnow,&
       ncolumnsDataOutESTM,ncolumnsDataOutDailyState

  IMPLICIT NONE

CONTAINS

  ! a mini version of SUEWS
  SUBROUTINE SuMin(&
       snowUse,EmissionsMethod,NetRadiationMethod,RoughLenHeatMethod,&! model options
       RoughLenMomMethod,StorageHeatMethod,AerodynamicResistanceMethod,OHMIncQF,&! model options
       iy,id,it,imin,isec,dt_since_start,tstep,tstep_prev,startDLS,endDLS,&! time-related input
       alt,lat,lng,Z,timezone,SurfaceArea,sfr,&! site-specific geographical settings
       z0m_in,zdm_in,&! roughness related settings
       alb,emis,SnowAlb,OHM_coef,WaterDist,&
       AHProf_24hr,HumActivity_24hr,PopProf_24hr,TraffProf_24hr,WUProfA_24hr,WUProfM_24hr,&
       qn1_av,dqndt,qn1_s_av,dqnsdt,&
       surf,DecidCap_id,albDecTr_id,albEveTr_id,albGrass_id,porosity_id,&
       GDD_id,HDD_id,LAI_id,WUDay_id,soilmoist_id,state_id,MeltWaterStore,&
       avkdn,avRh,avU1,Press_hPa,Temp_C,Precip,& ! forcing variables
       qh,qe,qsfc,tsk,CHKLOWQ)!output

    ! model configurations
    INTEGER,INTENT(in) ::snowUse
    INTEGER,INTENT(in) ::EmissionsMethod
    INTEGER,INTENT(in) ::NetRadiationMethod
    INTEGER,INTENT(IN) ::RoughLenHeatMethod
    INTEGER,INTENT(IN) ::RoughLenMomMethod
    INTEGER,INTENT(IN) ::StorageHeatMethod
    INTEGER,INTENT(IN) ::AerodynamicResistanceMethod
    INTEGER,INTENT(IN) ::OHMIncQF  !OHM calculation uses Q* only (0) or Q*+QF (1)

    ! time-related input
    INTEGER,INTENT(IN) ::iy
    INTEGER,INTENT(IN) ::id
    INTEGER,INTENT(IN) ::it
    INTEGER,INTENT(IN) ::imin
    INTEGER,INTENT(in) ::isec
    INTEGER,INTENT(in) ::dt_since_start ! time since simulation starts [s]
    INTEGER,INTENT(IN) ::tstep
    INTEGER,INTENT(IN) ::tstep_prev ! tstep size of the previous step
    INTEGER,INTENT(IN) ::startDLS ! start of daylight saving (inclusive)
    INTEGER,INTENT(IN) ::endDLS ! end of daylight saving (inclusive)


    ! site-specific geographical settings
    REAL(KIND(1D0)),INTENT(IN)              ::alt
    REAL(KIND(1D0)),INTENT(IN)              ::lat
    REAL(KIND(1D0)),INTENT(IN)              ::lng
    REAL(KIND(1D0)),INTENT(IN)              ::Z
    REAL(KIND(1D0)),INTENT(IN)              ::timezone
    REAL(KIND(1D0)),INTENT(IN)              ::SurfaceArea
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN) ::sfr


    ! roughness related settings
    REAL(KIND(1D0)),INTENT(in) ::z0m_in
    REAL(KIND(1D0)),INTENT(in) ::zdm_in

    ! radiation related settings:
    REAL(KIND(1D0)),DIMENSION(7),INTENT(INOUT) ::alb
    REAL(KIND(1D0)),DIMENSION(7),INTENT(IN)    ::emis
    REAL(KIND(1D0)),INTENT(INOUT) ::SnowAlb

    ! OHM coeffcients
    REAL(KIND(1D0)),DIMENSION(7+1,4,3),INTENT(IN) ::OHM_coef

    ! hydrology related settings
    REAL(KIND(1D0)),DIMENSION(7+1,7-1),INTENT(IN) ::WaterDist

    ! profiles at 24 hours
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: AHProf_24hr
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: HumActivity_24hr
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: PopProf_24hr
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: TraffProf_24hr
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: WUProfA_24hr
    REAL(KIND(1D0)),DIMENSION(0:23,2),INTENT(in) :: WUProfM_24hr


    ! daily states, also initial conditions

    REAL(KIND(1d0)),INTENT(INOUT) ::qn1_av
    REAL(KIND(1d0)),INTENT(INOUT) ::dqndt
    REAL(KIND(1d0)),INTENT(INOUT) ::qn1_s_av !qn1_av for snow
    REAL(KIND(1d0)),INTENT(INOUT) ::dqnsdt !dqndt for snow
    REAL(KIND(1d0)),INTENT(INOUT) ::DecidCap_id
    REAL(KIND(1d0)),INTENT(INOUT) ::albDecTr_id
    REAL(KIND(1d0)),INTENT(INOUT) ::albEveTr_id
    REAL(KIND(1d0)),INTENT(INOUT) ::albGrass_id
    REAL(KIND(1d0)),INTENT(INOUT) ::porosity_id
    REAL(KIND(1d0)),DIMENSION(5),INTENT(INOUT)   ::GDD_id     !Growing Degree Days (see SUEWS_DailyState.f95)
    REAL(KIND(1d0)),DIMENSION(6,2),INTENT(INOUT)   ::HDD_id     !Growing Degree Days (see SUEWS_DailyState.f95)
    REAL(KIND(1d0)),DIMENSION(3),INTENT(INOUT)   ::LAI_id     !LAI for each veg surface [m2 m-2]
    REAL(KIND(1d0)),DIMENSION(9),INTENT(INOUT)   ::WUDay_id
    REAL(KIND(1D0)),DIMENSION(7),INTENT(INOUT)   ::soilmoist_id
    REAL(KIND(1D0)),DIMENSION(7),INTENT(INOUT)   ::state_id
    REAL(KIND(1D0)),DIMENSION(7),INTENT(INOUT)   ::MeltWaterStore
    REAL(KIND(1D0)),DIMENSION(6,7),INTENT(INOUT) ::surf

    ! forcing variables
    REAL(KIND(1D0)),INTENT(IN)::avkdn
    REAL(KIND(1D0)),INTENT(IN)::avRh
    REAL(KIND(1D0)),INTENT(IN)::avU1
    REAL(KIND(1D0)),INTENT(IN)::Press_hPa
    REAL(KIND(1D0)),INTENT(IN)::Temp_C
    REAL(KIND(1D0)),INTENT(IN)::Precip

    ! output for WRF
    REAL(KIND(1D0)),INTENT(out)::qh
    REAL(KIND(1D0)),INTENT(out)::qe
    REAL(KIND(1D0)),INTENT(out)::qsfc
    REAL(KIND(1D0)),INTENT(out)::tsk
    REAL(KIND(1D0)),INTENT(out)::CHKLOWQ


    ! fixed settings in SuMin
    INTEGER,PARAMETER ::veg_type        = 1
    INTEGER,PARAMETER ::gsModel         = 2
    INTEGER,PARAMETER ::StabilityMethod = 3
    INTEGER,PARAMETER ::SMDMethod       = 0
    INTEGER,PARAMETER ::DiagQN          = 0
    INTEGER,PARAMETER ::DiagQS          = 0
    INTEGER,PARAMETER ::Diagnose        = 0
    INTEGER,PARAMETER ::ity             = 2
    INTEGER,PARAMETER ::LAICalcYes      = 1
    INTEGER,PARAMETER ::WaterUseMethod  = 0

    REAL(KIND(1D0)),PARAMETER:: LAI_obs   = 0
    REAL(KIND(1D0)),PARAMETER:: ldown_obs = 0
    REAL(KIND(1D0)),PARAMETER:: fcld_obs  = 0
    REAL(KIND(1D0)),PARAMETER:: snow_obs  = 0
    REAL(KIND(1D0)),PARAMETER:: qn1_obs   = 0
    REAL(KIND(1D0)),PARAMETER:: qh_obs    = 0
    REAL(KIND(1D0)),PARAMETER:: qf_obs    = 0
    REAL(KIND(1D0)),PARAMETER:: qs_obs    = 0

    REAL(KIND(1D0)),DIMENSION(7),PARAMETER ::SoilDepth = 0.2


    ! local variables not used for WRF coupling
    INTEGER::Gridiv
    INTEGER::Ie_end
    INTEGER::Ie_start

    ! parameters used in SUEWS for now:
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER:: SoilStoreCap=[150., 150., 150., 150., 150., 150., 0.]        !Capacity of soil store for each surface [mm]

    REAL(KIND(1d0)),PARAMETER:: AlbMin_DecTr=0.12   !Min albedo for deciduous trees [-]
    REAL(KIND(1d0)),PARAMETER:: AlbMax_DecTr=0.18   !Max albedo for deciduous trees [-]
    REAL(KIND(1d0)),PARAMETER:: AlbMin_EveTr=0.11   !Min albedo for evergreen trees [-]
    REAL(KIND(1d0)),PARAMETER:: AlbMax_EveTr=0.12   !Max albedo for evergreen trees [-]
    REAL(KIND(1d0)),PARAMETER:: AlbMin_Grass=0.18   !Min albedo for grass [-]
    REAL(KIND(1d0)),PARAMETER:: AlbMax_Grass=0.21    !Max albedo for grass [-]

    REAL(KIND(1d0)),PARAMETER:: CapMin_dec=0.3   !Min storage capacity for deciduous trees [mm] (from input information)
    REAL(KIND(1d0)),PARAMETER:: CapMax_dec=0.8   !Max storage capacity for deciduous trees [mm] (from input information)
    REAL(KIND(1d0)),PARAMETER:: PorMin_dec=0.2   !Min porosity for deciduous trees
    REAL(KIND(1d0)),PARAMETER:: PorMax_dec=0.6   !Max porosity for deciduous trees

    REAL(KIND(1d0)),PARAMETER:: FAIbldg    = 0. !Frontal area fraction of buildings
    REAL(KIND(1d0)),PARAMETER:: FAIEveTree = 0. !Frontal area fraction of evergreen trees
    REAL(KIND(1d0)),PARAMETER:: FAIDecTree = 0. !Frontal area fraction of deciduous trees

    REAL (KIND(1d0)),PARAMETER :: bldgH    = 10 !Mean building height
    REAL (KIND(1d0)),PARAMETER :: EveTreeH = 10 !Height of evergreen trees
    REAL (KIND(1d0)),PARAMETER :: DecTreeH = 10 !Height of deciduous trees

    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: BaseT          = [5,5,5]          !Base temperature for growing degree days [degC]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: BaseTe         = [11,11,11]       !Base temperature for senescence degree days [degC]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: GDDFull        = [300,300,300]    !Growing degree days needed for full capacity [degC]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: SDDFull        = [-450,-450,-450] !Senescence degree days needed to initiate leaf off [degC]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: LaiMin         = [4.,1.,1.6]      !Min LAI [m2 m-2]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: LaiMax         = [5.1,5.5,5.9]    !Max LAI [m2 m-2]
    REAL(KIND(1d0)),DIMENSION(3),PARAMETER:: MaxConductance = [7.4,11.7,30.1]  !Max conductance [mm s-1]

    REAL(KIND(1d0)),DIMENSION(4,3),PARAMETER:: LaiPower = RESHAPE(& !Coeffs for LAI equation: 1,2 - leaf growth; 3,4 - leaf off
         [[0.03,0.03,0.03],                                       &
         [0.0005,0.0005,0.0005],                                  &
         [0.03,0.03,0.03],                                        &
         [0.0005,0.0005,0.0005]],                                 &
         [4,3])

    INTEGER,DIMENSION(3),PARAMETER:: LAIType=0     !LAI equation to use: original (0) or new (1)


    REAL (KIND(1D0)),PARAMETER ::DRAINRT       = 0.25 !Drainage rate of the water bucket [mm hr-1]
    REAL (KIND(1D0)),PARAMETER ::RAINCOVER     = 1
    REAL (KIND(1D0)),PARAMETER ::RAINMAXRES    = 10   !Maximum water bucket reservoir [mm]
    REAL (KIND(1d0)),PARAMETER ::FlowChange    = 0    !Difference between the input and output flow in the water body
    REAL (KIND(1d0)),PARAMETER ::PipeCapacity  = 100  !Capacity of pipes to transfer water
    REAL (KIND(1d0)),PARAMETER ::RunoffToWater = 0.1  !Fraction of surface runoff going to water body


    REAL(KIND(1d0)),DIMENSION(7),PARAMETER:: StateLimit = [0.48, 0.25, 1.3, 0.8, 1.9, 1.0, 30000.] !Limit for state of each surface type [mm] (specified in input files)
    REAL(KIND(1d0)),DIMENSION(7),PARAMETER:: WetThresh  = [0.48, 0.25, 1.3, 0.8, 1.9, 1., 0.5]     !When State > WetThresh, rs=0 limit in SUEWS_evap [mm] (specified in input files)

    ! ---- Drainage characteristics ----------------------------------------------------------------
    ! 1 - min storage capacity [mm]
    ! 2 - Drainage equation to use
    ! 3 - Drainage coeff 1 [units depend on choice of eqn]
    ! 4 - Drainage coeff 2 [units depend on choice of eqn]
    ! 5 - max storage capacity [mm]
    REAL(KIND(1d0)),DIMENSION(5,7),PARAMETER:: surf_attr = RESHAPE(& ! variable to store the above five properties
         [[ 0.48 , 0.25 , 1.3 , 0.3 , 1.9 , 0.8 , 0.5 ],           &
         [ 3. , 3. , 2. , 2. , 2. , 3. , 0. ],                     &
         [10. , 10. , 0.013, 0.013, 0.013, 10. , 0. ],             &
         [ 3. , 3. , 1.71 , 1.71 , 1.71 , 3. , 0. ],               &
         [ 0.48 , 0.25 , 1.3 , 0.8 , 1.9 , 0.8 , 0.5 ]],           &
         [5,7])


    ! these will be assigned locally as data
    ! use gsModel=2 as in Ward et al. (2016)
    REAL (KIND(1d0)),PARAMETER::th   = 40   !Maximum temperature limit
    REAL (KIND(1d0)),PARAMETER::tl   = -10  !Minimum temperature limit
    REAL (KIND(1d0)),PARAMETER::Kmax = 1200 !Annual maximum hourly solar radiation
    REAL (KIND(1d0)),PARAMETER::g1   = 3.5  !Fitted parameters related to
    REAL (KIND(1d0)),PARAMETER::g2   = 200
    REAL (KIND(1d0)),PARAMETER::g3   = 0.1
    REAL (KIND(1d0)),PARAMETER::g4   = 0.7
    REAL (KIND(1d0)),PARAMETER::g5   = 30
    REAL (KIND(1d0)),PARAMETER::g6   = 0.05 !Fitted parameters related to
    REAL (KIND(1d0)),PARAMETER::s1   = 5.56
    REAL (KIND(1d0)),PARAMETER::s2   = 0    !surface res. calculations



    REAL(KIND(1d0)),DIMENSION(7+1),PARAMETER:: OHM_threshSW = [10,10,10,10,10,10,10,10]         !Arrays for OHM thresholds
    REAL(KIND(1d0)),DIMENSION(7+1),PARAMETER:: OHM_threshWD = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9] !Arrays for OHM thresholds

    REAL (KIND(1d0)),PARAMETER::  BaseTHDD=18.9  !Base temperature for QF


    ! local variables
    REAL(KIND(1D0))::CRWmax
    REAL(KIND(1D0))::CRWmin

    REAL(KIND(1D0))::EF_umolCO2perJ
    REAL(KIND(1D0))::EnEF_v_Jkm

    REAL(KIND(1D0))::Faut
    REAL(KIND(1D0))::FcEF_v_kgkm


    REAL(KIND(1D0))::FrFossilFuel_Heat
    REAL(KIND(1D0))::FrFossilFuel_NonHeat

    REAL(KIND(1D0))::InternalWaterUse_h
    REAL(KIND(1D0))::IrrFracConif
    REAL(KIND(1D0))::IrrFracDecid
    REAL(KIND(1D0))::IrrFracGrass


    REAL(KIND(1D0))::MaxQFMetab
    REAL(KIND(1D0))::MinQFMetab
    REAL(KIND(1D0))::NARP_EMIS_SNOW
    REAL(KIND(1D0))::NARP_TRANS_SITE
    REAL(KIND(1D0))::NumCapita

    REAL(KIND(1D0))::PopDensDaytime
    REAL(KIND(1D0))::PopDensNighttime

    REAL(KIND(1D0))::PrecipLimit
    REAL(KIND(1D0))::PrecipLimitAlb
    REAL(KIND(1D0))::RadMeltFact

    REAL(KIND(1D0))::SnowAlbMax
    REAL(KIND(1D0))::SnowAlbMin
    REAL(KIND(1D0))::SnowDensMax
    REAL(KIND(1D0))::SnowDensMin
    REAL(KIND(1D0))::SnowLimBuild
    REAL(KIND(1D0))::SnowLimPaved

    REAL(KIND(1D0))::tau_a
    REAL(KIND(1D0))::tau_f
    REAL(KIND(1D0))::tau_r

    REAL(KIND(1D0))::TempMeltFact

    REAL(KIND(1D0))::TrafficUnits
    REAL(KIND(1D0))::xsmd




    REAL(KIND(1D0)),DIMENSION(2)               ::AH_MIN
    REAL(KIND(1D0)),DIMENSION(2)               ::AH_SLOPE_Cooling
    REAL(KIND(1D0)),DIMENSION(2)               ::AH_SLOPE_Heating
    REAL(KIND(1D0)),DIMENSION(2)               ::QF0_BEU
    REAL(KIND(1D0)),DIMENSION(2)               ::Qf_A
    REAL(KIND(1D0)),DIMENSION(2)               ::Qf_B
    REAL(KIND(1D0)),DIMENSION(2)               ::Qf_C
    REAL(KIND(1D0)),DIMENSION(2)               ::T_CRITIC_Cooling
    REAL(KIND(1D0)),DIMENSION(2)               ::T_CRITIC_Heating
    REAL(KIND(1D0)),DIMENSION(2)               ::TrafficRate
    REAL(KIND(1D0)),DIMENSION(3)               ::Ie_a
    REAL(KIND(1D0)),DIMENSION(3)               ::Ie_m

    REAL(KIND(1D0)),DIMENSION(7)               ::DayWat
    REAL(KIND(1D0)),DIMENSION(7)               ::DayWatPer

    ! AnOHM related: not used
    REAL(KIND(1D0)),DIMENSION(7)               ::chAnOHM
    REAL(KIND(1D0)),DIMENSION(7)               ::cpAnOHM
    REAL(KIND(1D0)),DIMENSION(7)               ::kkAnOHM
    REAL(KIND(1D0)),DIMENSION(:,:),ALLOCATABLE ::MetForcingData_grid

    REAL(KIND(1D0)),DIMENSION(7)               ::SatHydraulicConduct

    REAL(KIND(1D0)),DIMENSION(7)               ::snowD


    REAL(KIND(1D0)),DIMENSION(3)               ::alpha_bioCO2
    REAL(KIND(1D0)),DIMENSION(3)               ::alpha_enh_bioCO2

    REAL(KIND(1D0)),DIMENSION(3)               ::beta_bioCO2
    REAL(KIND(1D0)),DIMENSION(3)               ::beta_enh_bioCO2

    REAL(KIND(1D0)),DIMENSION(3)               ::min_res_bioCO2
    REAL(KIND(1D0)),DIMENSION(3)               ::resp_a
    REAL(KIND(1D0)),DIMENSION(3)               ::resp_b

    REAL(KIND(1D0)),DIMENSION(0:23,2)          ::snowProf_24hr
    REAL(KIND(1D0)),DIMENSION(3)               ::theta_bioCO2
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE   ::Ts5mindata_ir !TODO:allocatable array can't serve as argument?


    REAL(KIND(1D0))                            ::SnowfallCum
    REAL(KIND(1d0)),DIMENSION(24*3600/tstep)   ::Tair24HR


    REAL(KIND(1D0)),DIMENSION(7)               ::IceFrac
    REAL(KIND(1D0)),DIMENSION(7)               ::SnowDens
    REAL(KIND(1D0)),DIMENSION(7)               ::snowFrac
    REAL(KIND(1D0)),DIMENSION(7)               ::SnowPack


    REAL(KIND(1D0)),DIMENSION(5)                           ::datetimeLine
    REAL(KIND(1D0)),DIMENSION(ncolumnsDataOutSUEWS-5)      ::dataOutLineSUEWS
    REAL(KIND(1D0)),DIMENSION(ncolumnsDataOutSnow-5)       ::dataOutLineSnow
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutESTM-5)       ::dataOutLineESTM
    REAL(KIND(1d0)),DIMENSION(ncolumnsDataOutDailyState-5) ::DailyStateLine


    PRINT*,''
    PRINT*, 'soilmoist_id',soilmoist_id
    ! soilmoist_id=MERGE(soilmoist_id,soilmoist_id*0,soilmoist_id>0)
    ! print*, 'soilmoist_id modified',soilmoist_id
    PRINT*, 'state_id',state_id


    CALL SUEWS_cal_Main(&
         AerodynamicResistanceMethod,AH_MIN,AHProf_24hr,AH_SLOPE_Cooling,& ! input&inout in alphabetical order
         AH_SLOPE_Heating,&
         alb,AlbMax_DecTr,AlbMax_EveTr,AlbMax_Grass,&
         AlbMin_DecTr,AlbMin_EveTr,AlbMin_Grass,&
         alpha_bioCO2,alpha_enh_bioCO2,alt,avkdn,avRh,avU1,BaseT,BaseTe,&
         BaseTHDD,beta_bioCO2,beta_enh_bioCO2,bldgH,CapMax_dec,CapMin_dec,&
         chAnOHM,cpAnOHM,CRWmax,CRWmin,DayWat,DayWatPer,&
         DecTreeH,Diagnose,DiagQN,DiagQS,DRAINRT,&
         dt_since_start,dqndt,qn1_av,dqnsdt,qn1_s_av,&
         EF_umolCO2perJ,emis,EmissionsMethod,EnEF_v_Jkm,endDLS,EveTreeH,FAIBldg,&
         FAIDecTree,FAIEveTree,Faut,FcEF_v_kgkm,fcld_obs,FlowChange,&
         FrFossilFuel_Heat,FrFossilFuel_NonHeat,G1,G2,G3,G4,G5,G6,GDD_id,&
         GDDFull,Gridiv,gsModel,HDD_id,HumActivity_24hr,&
         IceFrac,id,Ie_a,Ie_end,Ie_m,Ie_start,imin,&
         InternalWaterUse_h,IrrFracConif,IrrFracDecid,IrrFracGrass,isec,it,ity,&
         iy,kkAnOHM,Kmax,LAI_id,LAICalcYes,LAIMax,LAIMin,LAI_obs,&
         LAIPower,LAIType,lat,ldown_obs,lng,MaxConductance,MaxQFMetab,&
         MeltWaterStore,MetForcingData_grid,MinQFMetab,min_res_bioCO2,&
         NARP_EMIS_SNOW,NARP_TRANS_SITE,NetRadiationMethod,&
         NumCapita,OHM_coef,OHMIncQF,OHM_threshSW,&
         OHM_threshWD,PipeCapacity,PopDensDaytime,&
         PopDensNighttime,PopProf_24hr,PorMax_dec,PorMin_dec,&
         Precip,PrecipLimit,PrecipLimitAlb,Press_hPa,&
         QF0_BEU,Qf_A,Qf_B,Qf_C,&
         qn1_obs,qh_obs,qs_obs,qf_obs,&
         RadMeltFact,RAINCOVER,RainMaxRes,resp_a,resp_b,&
         RoughLenHeatMethod,RoughLenMomMethod,RunoffToWater,S1,S2,&
         SatHydraulicConduct,SDDFull,sfr,SMDMethod,SnowAlb,SnowAlbMax,&
         SnowAlbMin,snowD,SnowDens,SnowDensMax,SnowDensMin,SnowfallCum,snowFrac,&
         SnowLimBuild,SnowLimPaved,snow_obs,SnowPack,SnowProf_24hr,snowUse,SoilDepth,&
         soilmoist_id,soilstoreCap,StabilityMethod,startDLS,state_id,StateLimit,&
         StorageHeatMethod,surf,SurfaceArea,Tair24HR,tau_a,tau_f,tau_r,&
         T_CRITIC_Cooling,T_CRITIC_Heating,Temp_C,TempMeltFact,TH,&
         theta_bioCO2,timezone,TL,TrafficRate,TrafficUnits,&
         TraffProf_24hr,Ts5mindata_ir,tstep,tstep_prev,veg_type,&
         WaterDist,WaterUseMethod,WetThresh,&
         WUDay_id,&
         DecidCap_id,&
         albDecTr_id,&
         albEveTr_id,&
         albGrass_id,&
         porosity_id,&
         WUProfA_24hr,&
         WUProfM_24hr,xsmd,Z,z0m_in,zdm_in,&
         datetimeLine,dataOutLineSUEWS,dataOutLineSnow,dataOutLineESTM,&!output
         DailyStateLine)!output


    qh=dataOutLineSUEWS(9)
    qe=dataOutLineSUEWS(10)
    qsfc=dataOutLineSUEWS(16)
    tsk=dataOutLineSUEWS(5)+273.15
    CHKLOWQ=0.02

    PRINT*,''
    PRINT*, 'avkdn,kup,ldown,lup,tsurf'
    PRINT*, dataOutLineSUEWS(1:5)
    PRINT*, 'qn1,qf,qs,qh,qe'
    PRINT*, dataOutLineSUEWS(6:10)
    PRINT*,''
    ! IF ( ABS(qe)>1000 ) THEN
    !    zdm_in=0.
    !    PRINT*, 10./zdm_in
    ! END IF

  END SUBROUTINE SuMin


END MODULE SuMin_Module
