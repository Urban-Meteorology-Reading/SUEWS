! sg feb 2012 - added some comments
! sg feb 2012 - changed number of surfaces to allocatable array
! lj jun 2012 - snow part added

!===================================================================================
 module allocateArray
   implicit none

   integer::PavSurf=1,& ! when all surfaces    BareSoilSurf= treated separately below
   			BldgSurf=2,&
            ConifSurf=3,&
            DecidSurf=4,&
            GrassISurf=5,&
            GrassUSurf=6,&  
            WaterSurf=7, &
            ExcessSurf=8, &   ! runoff or soil
            NSurfDoNotReceiveDrainage=0,& ! Number of surfaces do not receive drainage water - green roof
            ivConif=1,&   ! order when just vegetated surfaces
            ivDecid=2,&
            ivGrassI=3,&
            ivGrassU=4,&
            CreateAnnual,& 
            alldays=1,&
            NumberDailyVars=25  !Number of columns in daily/monthly output file
            
   integer, parameter:: nSurf=7      ! total number of surfaces
   integer, parameter:: NVegSurf=4   ! number of surfaces that are vegetated
   integer, parameter:: maxGrid=100  ! maximum number of surface grids -- could be much larger e.g 1000
   integer, parameter:: NrowOhm =50  ! number of rows in the Ohm coefficient file
   integer, parameter:: NCoeffOhm=4  ! number of coefficients per surface to use (4)
   integer, parameter:: Ndays =366   ! sg -- could read this in -- and allocate arrays
          
   integer:: NPeriodsPerYear         ! calculated based on interval -- which is read in
   integer:: numMonth=12             ! sg -- could read this in -- and allocate arrays
   integer, parameter:: NumberSnowLayer = 1000  !Number of snow layers  
   integer,parameter:: sideG=90        !Number of horizontal grids to which snow melting is calculated

   !DailyState
   real (kind(1d0)),dimension(0:NDays, 5):: GDD ! growing degree days - for LAI
   
   ! check this may not be needed to go to -4 -- as 5 day mean only seems to be used in -OHMnew 
   real (kind(1d0)),dimension(-4:NDays, 6):: HDD ! heating dd - for QF
                                                !Matrix for heating and cooling degree days and Tair 
    										 !Last column contains the 5-day running mean								

   real (kind(1d0)),dimension(-4:NDays,NVegSurf)::LAI ! LAI
   real (kind(1d0)),dimension(NVegSurf):: LAI0         ! initial LAI


   real (kind(1d0)),dimension(nsurf):: AddWater       !Additional water from other surfaces
   real (kind(1d0)),dimension(nsurf):: AddWaterRunoff !Part of outflows going to runoff/soil
   real (kind(1d0)),dimension(nsurf):: chang          !Change in surface stores
   real (kind(1d0)),dimension(nsurf):: drain          !Drainage of each surface type
   real (kind(1d0)),dimension(nsurf):: runoff        !Runoff from each surface type
   real (kind(1d0)),dimension(nsurf):: runoffOut     !Runoff of existing surface types
   real (kind(1d0)),dimension(nsurf):: runoffSoil    !Soil runoff from each sub-surface
   real (kind(1d0)),dimension(nsurf):: runoffSoilOut !Soil runoff from existing soil stores
   real (kind(1d0)),dimension(nsurf):: SatHydraulicConduct ! Saturated Hydraulic conductivity 
   real (kind(1d0)),dimension(nsurf):: areasfr       !Area of each subsurface
  
   real (kind(1d0)),dimension(nsurf)::  smd_nsurf     !Soil moisture deficit of each sub-surface
   real (kind(1d0)),dimension(nsurf)::  smd_nsurfOut  !Soil moisture deficit in existing soil stores
   real (kind(1d0)),dimension(nsurf)::  SoilmoistOld  !Soil moisture from pervius timestep
   real (kind(1d0)),dimension(nsurf)::  soilmoist0    !Inital soil moisture
   real (kind(1d0)),dimension(nsurf)::  Soilmoist     !Soil moisture of each surface type
   real (kind(1d0)),dimension(nsurf)::  soilstorecap !Maximum capasity soil storage can hold for each surface
   real (kind(1d0)),dimension(nsurf)::  VolSoilMoistCap ! Maximum volumetric soil moisture capacity of each surface
   real (kind(1d0)),dimension(nsurf):: SoilStateOld,StateOld
 
   !========Variables related to NARP================================
   !Radiation balance components for different surfaces
   real (kind(1d0)),dimension(nsurf)::  Tsurf_ind,&
                                        Tsurf_ind_snow,&
                                        Tsurf_ind_nosnow 
   real (kind(1d0)),dimension(nsurf)::  kup_ind,&
                                        kup_ind_snow,&
                                        kup_ind_nosnow
   real (kind(1d0)),dimension(nsurf)::  lup_ind,&
                                        lup_ind_snow,&
                                        lup_ind_nosnow
   real (kind(1d0)),dimension(nsurf)::  qn1_ind,&
                                        qn1_ind_snow,&
                                        qn1_ind_nosnow
   real (kind(1d0)),dimension(nsurf)::  emis       ! emissivity
   real (kind(1d0)),dimension(nsurf)::  alb        !albedo
  
   !GIS information
   real (kind(1d0)),dimension(nsurf)::  sfr     ! surface fractions + bare soil
   real (kind(1d0))::BareSoilSurfFraction
 
   !surface states
   real (kind(1d0)),dimension(nsurf):: state      !Wetness status of each surface type in mm
   real (kind(1d0)),dimension(nsurf):: stateAll   !Wetness status of each surface type in mm 
   real (kind(1d0)),dimension(nsurf):: StateOut   !Wetness state of existing surfaces  
   real (kind(1d0)),dimension(nsurf):: state0     !inital Wetness status of each surface   state
 
   real (kind(1d0)),dimension(nsurf+1,nsurf-1)::WaterDist !Table which determines distribution of water in the canopy
   
 !===================Snow related===================================
   real (kind(1d0)),dimension(nsurf)::changSnow,&       !Change in snowpack in mm
   									  maxSnowVol,&      !! Maximum snow volume							  
									  MeltWaterStore,&  !!Liquid water in the snow pack of ith surface
                                      ev_snow,&        	!!Evaporation from snowpack in mm
                                      mw_ind,&         	!Melt water from individual surface in mm
                                      mw_indDay,&      	!!Melt water per day from each surface type in m3
									  runoffSnow,&     	!!Runoff from snowpack in mm and in m3
                                      densSnow,&        !Density of snow
                                      SnowDensInit,&
                                      snowFrac,&       	!!Surface fraction of snow cover
                                      iceFrac,&
                                      snowInit,&
                                      snowDepth,&      !Depth of snow in cm
                                      SnowToSurf,&     !Meltwater flowing from snow to surface
                                      volSWE,&
                                      StateFraction,&  !Fraction of state that can freeze                                      
                                      freezMelt,&      !Amount of freezing meltwater in mm for the ith surface area
                                      Qm_freezState,&  !Heat by freezing of surface state
                                      freezState,&     !Amount of freezing state in mm for the ith surface area
                                      FreezStateVol,&
                                      Qm_melt,&        !Heat consumption by snow melt
                                      Qm_rain,&        !Heat by rain falling on snow
                                      rainOnSnow,&     !Liquid precipitation falling on snow ()
                                      snowD,&
                                      deltaQi
 
   real (kind(1d0)),dimension(nsurf)::snowPack,&       !Amount of snow on each surface in mm
                                      snowPackOld
   integer,dimension(nsurf):: heiG,&                    !snow layer height
                              snowCoverForms,&
                              snowCalcSwitch=0          !Defines if snow related balance is made
                             
   real(kind(1d0)),dimension (NumberSnowLayer,sideG,nsurf):: snowMat=0 !Matrix of snow, size limited
   
   !==================Outer programs================================= 
   character(len=15),dimension(2,maxGrid)::GridConnections  !List of different grid corrections
   real (kind(1d0)),dimension(maxGrid)::   GridConnectionsFrac   !Fraction of water moving between the different grids
   integer,dimension(maxGrid):: laiID   !Day of year of LAI from previous year
   
   real (kind(1d0)),dimension(maxGrid,nVegSurf)::laiGrids  !Initializing LAI for vegetation surfaces    
   real(kind(1d0)),dimension(maxGrid,nsurf)::stateGrids   !Initializing states (max grid amount 1000)
   real(kind(1d0)),dimension(maxGrid,nsurf):: soilmoistGrids   
         
   real (kind(1d0)), dimension(0:ndays)::albDec    !Changing albedo for deciduous trees
   real (kind(1d0)), dimension(0:ndays)::DecidCap  !Max capacity of deciduous trees
   real (kind(1d0)), dimension(0:ndays):: porosity !Porosity of deciduous trees
   real (kind(1d0)), dimension(0:ndays,6):: WU_Day !Daily water use:  total,automatic,manual
   real (kind(1d0)):: AlbMin_dec=0.15,&
   					  				AlbMax_dec,&
                      CapMin_dec,&
                      CapMax_dec                                
   real (kind(1d0)), dimension(5,nsurf):: surf     !Surface information matrix defined in LUMPS_Initial.f90
       
   ! OHMNew
   real (kind(1d0)),dimension(4,nsurf+2)::co2use
   real (kind(1d0)),dimension(NrowOhm,NCoeffOhm):: co 
   
   ! LUMPS vegetation phenology
   integer:: jj1,jj2,jj3,jj4, changed,changeInit,LAItype
   real (kind(1d0)):: kk1,kk2,dmax,dmin,dx1,dx2    ! used to calculate Vegetation phenology for LUMPS
   real(kind(1d0)),dimension(0:23)::runT           ! running average T for the day
   real(kind(1d0)),dimension(0:23)::runP           ! running total Precip for the day
   real (kind(1d0))::avT_h, totP_h                 ! daily running average Temp, Total precip
   real(kind(1d0)),dimension(nvegsurf):: baseT     ! Base Temperature for growing degree days
   real(kind(1d0)),dimension(nvegsurf):: BaseTe    ! Base Temperature for senescence
   real(kind(1d0)),dimension(nvegsurf):: GDDFull   !Canopy growth parameters
   real(kind(1d0))::GDDmax,SDDMax                  ! limit across vegetation types
   real(kind(1d0)),dimension(nvegsurf):: LaiMin    ! minimum LAI
   real(kind(1d0)),dimension(nvegsurf):: LaiMax    ! maximum LAI
  ! real(kind(1d0)):: MaxLaiMax                     ! maximum LAI of all  sg 12 nov 13 -- no longer used
   real(kind(1d0)),dimension(nvegsurf):: SDDFull   ! Vegetation parameters 
   real(kind(1d0)),dimension(nvegsurf+2):: MaxConductance      ! maximum conductance for surface resistance
                   ! +2 so can use surface position in array
                   
   real(kind(1d0)),dimension(4):: laiPower
                    
   ! used in water Use, anthropogenic heat and storage calculation                   
   integer,dimension(0:Ndays,3)::DayofWeek  ! day of week, month, season
   real (kind(1d0)), dimension(3)::WU_prof      !% of water use profile
       
 end MODULE allocateArray
 
!======================================================================================================
MODULE cbl_MODULE

	integer::EntrainmentType,&  ! Entrainment type choice
			 CO2_included,&     ! CO2 included
             InitialData_use,&  ! 1 read initial data, 0 do not        
             qh_choice,&        ! selection of qh use to drive CBL growth 1=Suews 2=lumps 3=obs
             sondeflag      ! 1 read sonde or vertical profile data in 0 do not
       
    integer,dimension(366)::cblday=0   
      
    character (len=200), dimension(366)::FileSonde=""
    character (len=200)::InitialDataFileName    
    real(kind(1D0)):: wsb       ! subsidence velocity    
    real(kind(1d0)),dimension(0:1,0:9):: cbldata
    real(kind(1d0)),dimension(0:9)::cbld
    real(kind(1d0)),dimension(:,:),allocatable::IniCBLdata
    
  !Parameters in CBL code         
    integer::tstep_s,&
			 which_day,& 
             zmax,&            
             start1,&
             start2,&
             icount,&
             jday 
    integer::C2K=273.16,&
             nEqn=4      
             
               
 real (kind(1D0)):: usbl,ftbl,fqbl,fcbl,gamt,gamq,gamc,tpp,qpp,cp0!,tk
 
 real(kind(1D0))::alpha3,&
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

                                      
    real (kind(1D0)), dimension (0:500,2):: gtheta,ghum ! Vertical gradient of theta and specific humidity from sonde data
    real (kind(1D0)), dimension(4)::y
 
   END   MODULE cbl_MODULE
 !===================================================================================
 module snowMod
     implicit none
 
     real (kind(1D0))::AdjMeltFact,&	  !Factor between melt and freezing factors
     				   ALB_SNOW,&
                       albSnowMin,&       !Minimum snow albedo
                       albSnowMax,&       !Maximum snow albedo
                       CumSnowfall,&        !Cumulative snowfall
                       densSnowMin,&      !Minimum density of snow
                       densSnowMax,&      !Maximum density of snow
                       fwh,&              !Weighted freezing water
                       lvS_J_kg,&         !Latent heat of sublimation in J/kg
                       mwh,&              !Weighted hourly water melt
                       MwStore,&        		  !Meltwater storage	
                       PrecipLimit,&      !Temperature limit when precipitation occurs as snow 
                       PrecipLimitAlb,&   !Precipitation limit for albedo change (in mm)
                       Qm,&               !Snow melt associated heat flux
                       QmFreez,&          !Energy released in freezing of meltwater or surface state
                       QmRain,&
                       qn1_snow,&         !Net all-wave radiation of snowpack
                       qn1_nosnow,&       !Same for the snow free surface
                       RadMeltFact,&      !Radiation melt factor
                       SnowLimBuild,&     !Snow removal limits for roofs in mm)
					   SnowLimPaved,&     !Snow removal limits for paved surfaces in mm)
                       swe,&			  !Weighted snow water equivalent (in mm)
                       tau_a,&            !Time constans related to albedo change
                       tau_f,&
                       tau_r,&            !Time constant for density increase.
                       TempMeltFact,&     !Temperature melt factor
                       volDay,&           !Volume of the melted water per day 
                       zf,&
                       WaterHoldCapFrac,& !Water holding capacity factor
                       CRWmin,& !Free water holding capacity of deep snowpack
                       CRWmax  !Free water holding capacity of shallow snowpack
                        
     real(kind(1D0)), dimension(2)::  SnowRemoval=0 ! Removal of snow in mm 
     real(kind(1d0)), dimension(0:23):: snowProf  ! Timing of snow removal (hourly)
     
     integer::SnowFractionChoice   !Choice how fraction of snow is calculated
 
 end module snowMod

 !===================================================================================
 Module defaultNotUsed
 	implicit none
 	real (kind(1d0)):: notUsed=-55.55,reall,NAN=-999,pNAN=999
 	integer:: notUsedI=-55, ios_out,errorChoice  !errorChoice defines if the problemfile is opened for the first time
 end Module defaultNotUsed
 
 !===================================================================================
 MODULE data_in
 implicit  none
 
 !In alphabethical order    
     
 real (kind(1d0)):: AH_MIN,&   !Minimum anthropogenic heat flux
                    AH_SLOPE,& !Slope of the antrhropogenic heat flux calculation
                    alpha_qhqe,& !Alpha parameter used in LUMPS QH and QE calculations
                    avdens,&    !Average air density
                    avkdn,&     !Average downwelling shortwave radiation
                    avrh,&      !Average relative humidity
                    avts,&      !Average surface temperature
                    avu1,&      !Average wind speed
                    azimuth,&   !Sun azimuth in degrees
                    BaseTHDD,&             !Base temperature for QF               
                    defaultQf,&             !Default anthropogenic heat flux
                    defaultQs,&             !Default storage heat flux
                    E_mod,&      !Modelled latent heat flux with LUMPS
                    emis_snow,&            !Emissivity of snow
                    fcld,&       !Cloud fraction modelled
                    fcld_obs,&  !Cloud fraction observed
                    h_mod,&      !Modelled sensible heat flux with LUMPS
                    kclear,&  !Theoretical downward shortwave radiation
                    kdiff,&   !Diffuse shortwave radiation
                    kdir,&    !Direct shortwave radiation
                    kup,&     !Upward 
                    lai_hr,&    !LAI                
                    lat,&                  !Latitude
                    ldown, &  !Downward longwave radiation
                    ldown_obs,& !Downwelling longwave radiation
                    lng,&       !Longitude
                    lup,&     !Upward longwave radiation             
                    numCapita,&  !Number of people in the study area
                    Precip_hr,&  !Hourly precipitation (mm)
                    Press_hPa,&  !Station air pressure in hPa     
                    Pres_kPa,&   !Station air pressure in kPa
                    q1_noqf,&
                    qe,&        !Observed latent heat flux
                    qe_obs,&
                    qf,&        !Observed antrhropogeni heat flux
                    QF_SAHP,&    !Anthropogenic heat flux calculated by SAHP
                    qh,&        !Observed sensible heat flux
                    qh_obs,&
                    qn1,&       !Net all-wave radiation for the study area
                    qn1_bup,&
                    qn1_obs,&   !Observed new all-wave radiation
                    qn1_S,&     !Total net all-wave radiation for the snowfree surface
                    qn1_SF,&    !Total net all-wave radiation for the snowpack
                    qs,&        !Observed storage heat flux
                    snow,&      !snow cover
                    snow_obs,&  !Observed snow cover
                    T_CRITIC,& !Critical temperature
                    Temp_C,&    !Air temperature
                    timezone,&  !Timezone (GMT=0)
                    trans_site,&  !Atmospheric transmittivity
                    tsurf,&   !Surface temperature
                    wdir,&      ! Wind direction
                    wuh,&       !Hourly water use (mm)
                    xsmd,&      !Measured soil moisture deficit
                    year,&      !Year of the measurements
                    zenith_deg    !Sun zenith angle in degrees

  integer::FirstYear, nlines  !Number of lines in met forcing file.
              
  real (kind(1d0)),dimension(366,25)::day,season,month,yr_tot,all_tot !daily matrixes
  integer,dimension(:,:), allocatable:: dataMet1              !Meteorological input matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataMet2              !Meteorological input matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataOut1             !Main output matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataOut2             !NARP output matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataOut3             !Snow output matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataOut5min
  real(kind(1d0)),dimension(:,:), allocatable:: dataOutBL    !CBL output matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataOutSOL   !SOLWEIG POI output

  !Testing different annual reading
  integer,dimension(:,:), allocatable:: AnnualFileIN

  integer::AlbedoChoice,&         !If albedos dependency on zenith angle is taken into account
           AnthropHeatChoice,&    !Is anthropogenic heat calculated
           CBLuse,&               !s.o.
           commonchoiceallsites,& !Determines if multiple sites are considered  impacts OHM sub-model
           gisinputtype, &        !GisInputType -- 1/2/3/4
           Inputmetformat,&       !Defines format for met input data
           Interval,&             !Measurement interval in seconds
           ldown_option,&         !What parameterization is used for downward longwave radiation 1-2-3
           lfnout,&               !Error Output write units
           lfnoutC,&              !Clean output write units
           lfnOld,&
           lfnSAHP,&              !Number of SAHP file
           NARPOutput,&           !Defines if radiation components are separatley printed out
           netradiationchoice,&   !Is net all-wave radiation modeled (=2) or measured (=1)
           qschoice,&             !Defines if QS is calculated
           qual,&                 !Counter for writing 
           SkipHeaderGis,&        !Number of lines in gis file that are skipped
           SkipHeaderMet,&        !Number of lines to skip in met file input
           SNOWuse,&
           SOLWEIGout,&           !Calculates Tmrt and other fluxes on a grid, FL 
           write5min,&            !Defines if 5-min output is printed
           writedailyState=1
                   
      character (len=150)::FileInputPath,&   !Filepath for input file
                           FileOutputPath,&  !Filepath for output file
                           fileout,&         !Output file name   
                           FileErrorInf,&    !Output error file name
                           filechoices,&     !File with output of run characteristics
                           filemet,&         !Meteorological input filename
                           filegis,&         !GIS input datafile
                           fileOHM,&         !File name with ohm coefficients
                           FileOHMChoices,&  !File name where possible a1-a3 coefficients are located
                           fileWU,&          !WU filename
                           fileMonthly,&     !Monthly output file name
                           fileDaily,&       !Daily output filename
                           file5min,&        !5min output filename
                           NARPOut,&         !File name for NARP output
                           SnowOut,&         !File name for snow output file
                           SOLWEIGpoiOut,&   !File name for SOLWEIG poi file
                           BLOut,&           !File name for BL output file
                           FileSAHP          !SAHP file name
                           
      character (len=20)::FileCode,FileCodeO
      character (len=90),dimension(14)::keepheader
      character (len=90)::progname='SUEWS V2014b'  !!!!<<<<<<<<<<<<<<<<<<
        
      logical:: finish,once
 
      real(kind(1d0)),dimension(2)::Qf_A,Qf_B,Qf_C !Qf coefficients
      
      integer,dimension(2)::DayLightSavingDay    !The date when it is changed to daylight saving (in DOY)
                          
      real(kind(1d0)), dimension(0:23,2):: AHPROF !Heat profile -- ?? check why not 0:23
 
      real(kind(1d0))::nCBLstep  !number of time steps of Runge-kutta methods in one hour 
         
  !---------Water bucket (see B. Offerle's PhD)----------------------------------
      real (KIND(1D0)):: DRAINRT,&    !Drainage rate of the water bucket
                            RAINBUCKET,& !RAINFALL RESERVOIR
                            RAINCOVER,&  
                            RAINMAXRES,& !Maximum water bucket reservoir
                            RAINRES,& 
                            TEMPVEG ! TEMPORARY VEGETATIVE SURFACE FRACTION ADJUSTED BY RAINFALL
      
  !---------SOLWEIG variables---------------------------------------------------
     real(kind(1D0))::absL,&           ! Absorption coefficient of longwave radiation of a person         
                      absK,&           ! Absorption coefficient of shortwave radiation of a person
                      heightgravity,&  ! Centre of gravity for a standing person
                      TransMin,&         ! Tranmissivity of K through decidious vegetation (leaf on)
                      TransMax           ! Tranmissivity of K through decidious vegetation (leaf off)

     integer::Posture,&                ! 1.Standing, 2.Sitting
              usevegdem,& 	           ! With vegetation (1)
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
              SOLWEIG_ldown            ! 1= use SOLWEIG code to estimate Ldown, 0=use SEUWS 
              
     character (len=150)::DSMPath,&    ! Path to DSMs
                          DSMname,&    ! Ground and building DSM
                          CDSMname,&   ! Canopy DSM
                          TDSMname,&   ! Trunk zone DSM
                          SVFPath,&    ! Path to SVFs
                          SVFsuffix,&  !
                          buildingsname! Boolean matrix for locations of building pixels
     
             


 end module data_in
 !===================================================================================
 !**********************************************
 module time
   integer :: id,&            !Day of year
              it,&            !Hour
              iostat_var,&      !File status from reading data (should not be here)
              lastTimeofDAY=23, firstTimeofDay=0,&
              DLS                            !- day lightsavings =1 + 1h) =0  
 
             ! check what lasttime of day should be
   real (kind(1d0)):: dectime !Time is decimals
   real (kind(1d0)):: halftimestep !in decimal time based on interval in runcontrol
   real (kind(1d0)):: NPeriodsPerDay,&  ! number of time intervals
   hrcount ! count number of hours in this day
 integer::nofDaysThisYear   ! BASED ON WHETHEr leapyear or no
 end module time
 !===================================================================================
 module ohm_calc !OHM related variables
   implicit none
   real (kind(1d0)):: q1,q2,q3,& !Time derivatives of OHM coefficients
                      a1,a2,a3,&   !Coefficients for OHM
                      r1,r2,r3
                      ! ?? add what the dimensions are for this array
   real(kind(1d0)),dimension (9,4,3):: OHM_coef                   
   character (len=2)::ohm_txt   !Abbreviation concerning OHM calculation printed in output
 end module ohm_calc
 
 !===================================================================================
 module mod_grav    
   real (kind(1d0)):: grav=9.80665  !g - gravity - physics today august 1987
 end module mod_grav
 
 !===================================================================================
 module mod_k
   real(kind(1d0)) :: k=0.4,&             !Von Karman's contant
                      k2=0.16,&           !Power of Van Karman's contant
                      neut_limit=0.001000 !Limit for neutral stability
 end module mod_k
 
 !===================================================================================
 module gas
   !   press (mb) ea (mb)
   implicit none
   real (kind(1d0))::  comp=0.9995     
   real (kind(1d0))::  epsil=0.62197   !ratio molecular weight of water vapor/dry air (kg/mol/kg/mol)
   real (kind(1d0))::  epsil_gkg=621.97   !ratio molecular weight of water vapor/dry air in g/kg
   real (kind(1d0))::  dry_gas=8.31451 !Dry gas constant (J/k/mol) 
   real (kind(1d0))::  gas_ct_wat=461.05 !Gas constant for water (J/kg/K)
   real (kind(1d0))::  molar=0.028965 !Dry air molar fraction in kg/mol
   real (kind(1d0))::  molar_wat_vap=0.0180153 !Molar fraction of water vapor in kg/mol  
   real (kind(1d0))::  gas_ct_dry=8.31451/0.028965 !j/kg/k=dry_gas/molar
   real (kind(1d0))::  gas_ct_wv=8.31451/0.0180153 !j/kg/kdry_gas/molar_wat_vap 
 end module gas
 
 !**********************************************
 module mod_z
   real (kind(1d0)) :: zzd,&  !Active measurement height (meas. height-displac. height)
                       z0m,&  !Aerodynamic roughness length
                       zdm,&  !Displacement height
                       z      !Windspeed height
   real(kind(1E10))::z0V      !Roughness length for vapour
 end module mod_z
 
 !**********************************************
 module resist  !Variables related surface resistance calculations (P. 1744 in G&O1991)
 implicit none
 real (kind(1d0)):: th,&             !Maximum temperature limit
                        tl,&             !Minimum temperature limit
                        Kmax,&           !Annual maximum hourly solar radiation
                        g1,g2,g3,g4,&    !Fitted parameters related to 
                        g5,g6,s1,s2,&    !surface res. calculations
                        tc,&             !Temperature parameter 1 
                        tc2              !Temperature parameter 2
 end module resist
 
 !**********************************************
 module moist
   implicit none
   
   real (kind(1d0))::avcp,&        !Specific heat capacity                        
   					 dens_dry,&    !Dry air density kg m-3
                     dq,&          !Specific humidity deficit
                     Ea_hPa,&      !Water vapour pressure in hPa 
                     Es_hPa,&      !Saturation vapour pressure in hPa   
                     lv_J_kg,&     !Latent heat of vaporization in J/kg
                     psyc_hPa,&    !Psychometric constant in hPa
                     psycIce_hPa,& !Psychometric constant in hPa for snow
                     s_Pa,&        !Vapour pressure versus temperature slope in Pa
                     s_hpa,&       !Vapour pressure versus temperature slope in hPa
                     sIce_hpa,&    !Vapour pressure versus temperature slope in hPa above ice/snow
                     vpd_hPa,&     !Vapour pressure deficit in hPa
                     vpd_pa,&      !Vapour pressure deficit in Pa
                     waterDens=999.8395 !Density of water in 0 cel deg

 end module moist
 !**********************************************
 
 module gis_data
   implicit none                             
                                         
   real(kind(1d0)):: Alt,&                        !Altitude in m
   		             areaunir,&                   !Unirrigated area
                     areair,&                     !Irrigated area
                     bldgH,&                      !Mean building height
                     FAIbldg,&                    !Frontal area fraction of buildings
                     FAItree,&                    !Frontal area fraction of trees                   
                     grassfractionirrigated,&     !Irrigated grass fraction for LUMPS                                   
                     pavedfractionirrigated,&     !Irrigated paved area fraction for LUMPS
                     TreeH,&                      !Mean tree height
                     treefractionirrigated,&      !Irrigated tree fraction for LUMPS                    
                     veg_fr,&                     !Vegetation fraction from land area 
                                                  !- For LUMPS - dependent on user choice    & water
                     VegFraction, &               ! sum of vegetation -not including water
                     areaZh                       !=(sfr(BldgSurf)+sfr(ConifSurf)+sfr(DecidSurf)) !Total area of buildings and trees
 
   integer:: idgis,&      !Time integers used in the code    
             itgis,&      !
             Veg_type=1    !Defines how vegetation is calculated for LUMPS
   
 end module gis_data
 
 !************************************************************
 module sues_data
   implicit none                                  
       
   real(kind(1d0)),dimension(2):: SurPlus_evap !Surplus for evaporation in 5 min timestep
    ! sg -- need to determine size                 
                                  
   !Variables listed in SuesInput.nml                                      
   real (kind(1d0))::addImpervious,&      !Water from other grids of paved area
                     addPipes,&           !Water from other grids in pipes
                     addVeg,&             !Water from other grids of vegetation area
                     addWaterbody,&       !Water from from other grids in water body
					 DecidStorCapMin,&    !Storage capacity for deciduos trees in off-leaf period
                     FlowChange,&         !Difference between the input and output flow in the water body
                     PipeCapacity,&       !Capacity of pipes to transfer water
                     RunoffToWater,&      !Fraction of surface runoff goinf to water body
                     DecidStorCapMax,&    !Storage capacity for deciduos trees in on-leaf period
                     SmCap,&              !Volumetric/gravimetric soil moisture capacity
                     SoilDensity,&        !Bulk density of soil
                     SoilDepth,&          !Depth of the soil layer
                     SoilDepthMeas,&      !Soil depth of the measured soil moisture
                     SoilRocks,&          !Fraction of rocks in soil
                     SurfaceArea,&        !Surface area of the wateruse area
                     tstep,&              !Time step for SUES
                     WaterBodyType,&      !If water body type is pond/lake (=1) or river (=2)
                     WaterState,&         !State of the water body
                     WaterStorCap,&       !Capacity of water body when surface is wet
                     wateruseareaGrass,&  !Water use area
                     wateruseareaTrees,&
                     wuhTrees             !Water for trees/shrubs
                                  
   integer:: ity,&              !Evaporation calculated according to Shuttleworth (=2) or Rutter (=1)
             RoughLen_heat,&    !Method of calculating roughness length for heat
             smd_choice !Defines if soil moisture calculated or  observed (and type of input)
   
             
                    
   !Other related to SUES   
   real (kind(1d0))::AdditionalWater,&     !Water flow from other surfaces
                     ch_per_interval,&     !Change in state per interval
                     chSnow_per_interval,& !Change in snow state per interval
                     dI_dt,&               !Water flow between two stores
                     dr_per_interval,&     !Drainage per interval
                     ev_per_interval,&     !Evaporation per interval
                     evap_5min,&           !Evaporation per 5 minute
                     P,&
                     pin,&                 !Rain per time interval
                     planF,&               !Areally weighted frontal area fraction
                     rb,&                  !Boundary layer resistance
                     runoffAGimpervious,&  !Above ground runoff from impervious surface
                     runoffAGveg,&         !Above ground runoff from vegetated surfaces
                     runoffPipes,&         !Runoff in pipes
                     runoffWaterBody,&
                     runoff_per_interval,&
                     runoffSoil_per_interval,&
                     qe_per_interval,&     !latent heat per interval
                     soilmoist_state=0,&
                     soilmoistcap,&
                     soilstate,&
                     st_per_interval,&!Surface state per interval
                     stph1,&          !State per interval ??
                     surplusWaterBody,&
                     tlv,&            !Vaporization per timestep
                     tlv_sub,&
                     overuse=0,&
                     Zh               !Areally weighted roughness element height
 
   !Calculation of u*,stability and aerodynamic resistance  
   real (kind(1d0))::h,&          !Sensible heat flux used to calculate friction velocity
                     l_mod,&      !Monin-Obukhov length (either measured or modelled)
                     psim,&       !Stability function of momentum
                     psyh,&       !Stability function of heat
                     RA,&         !Aerodynamic resistance
                     RAsnow,&     !Aerodynamic resistance over snow
                     tstar,&      !T*
                     ustar,&      !Friction velocity
                     z0_gis       !Roughness length for momentum from gis input file
                     
   !Surface resistance related variables
   real (kind(1d0))::resistsurf,& !Surface resistance
                     gdq,&        !G(dq)
                     qnm,&        !QMAX/(QMAX+G2)
                     gq,&         !G(Q*)
                     gtemp,&      !G(T)
                     gl,&         !G(LAI)
                     sdp,&        !S1/G6+S2
                     smd,&        !Soil moisture deficit of the soil surface layer
                     gs,&         !G(Soil moisture deficit)
                     gsc          !Surface Layer Conductance      
                     
  !SUES latent heat flux related variables 
  real (kind(1d0))::  vdrc,&     !Second term up in calculation of E
                      e,&        !Upper part of the equation and e used as emissivity in NARP 
                      sp,&       !Term in calculation of E
                      sae,&      !Same
                      ev,&       !Evaporation
                      rst,&      !Some switch ??
                      qeph       !Latent heat flux (W m^-2)
                      
 
 !Water use related variables
  real (kind(1d0)):: ext_wuh,&        !External water use 
                     ext_wuhP,&
                     ext_wu,&         !
                     Faut,&           !Fraction of irrigated area using automatic irrigation
                     int_wuh,&        !Internal water use
                     int_wu,& 
                     IrrFractionTrees,& !Fraction of the surafce area of irrigated trees
                     IrrTrees,&         !Surface area fraction of irrigated trees
                     internalwateruse !Internal water use
  
 ! 7 - number of days in week                   
  real (kind(1d0)), dimension(7)::DayWatPer,&  !% of houses following daily water
                                  DayWat       !Days of watering allowed
  real (kind(1d0)), dimension(0:23)::HourWat !Hours of watering allowed+ % using
 
  real (kind(1d0)), dimension(3)::Ie_a,&
                                  Ie_m
  
  integer ::HourResChoice,&
            WU_choice,&  !Defines if water use is calculated
            Ie_start,&   !Starting time of water use (DOY)
            Ie_end       !Ending time of water use (DOY)     
 
   integer :: nsh,&
              nmin,&
              is,&             ! which surface
              in,&
              LAIcalcYes,&       ! 1 YES Calc LAI 0 do not
              z0_method,&  
              stabilitymethod !Method to calculate stability
                     
   !These are variables which currently have been removed from SuesInput.nml     
   integer::AerodynamicResistanceMethod=2 !The method used to calculate aerodynamic resistance
           
     
 end module sues_data
 
 !**********************************************
 !===================================================================================
 Module VegPhenogy
     implicit none
     real (kind(1d0)):: VegPhenLumps
 end Module VegPhenogy
 
 Module filename
     character (len=90)::  smithfile     !file for narp
 end Module filename

 
 Module InitialCond

   real (Kind(1d0))::LAIinitialET,&
                     LAIinitialDT,&
                     LAIinitialUG,&
                     LAIinitialIG,&
                     porosity0,&
                     DecidCap0,&
                     albDec0,&
                     Temp_C0,&
                     GDD_1_0,&
                     GDD_2_0,&
                     SnowWaterBldgState,&
                     SnowWaterETstate,&
                     SnowWaterDTState,&
                     SnowWaterIGState,&
                     SnowWaterUGState,&
                     SnowWaterPavstate,&
                     SnowWaterWaterstate,&
                     SnowPackBldg,&
                     SnowPackET,&
                     SnowPackDT,&
                     SnowPackIG,&
                     SnowPackUG,&
                     SnowPackWater,&
                     SnowPackPav
     integer::ID_Prev

 end Module InitialCond

 !C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95
 
 
 
 
 
 
 