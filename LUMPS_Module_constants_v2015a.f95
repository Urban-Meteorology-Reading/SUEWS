! sg feb 2012 - added some comments
! sg feb 2012 - changed number of surfaces to allocatable array
! lj jun 2012 - snow part added
! HW, LJ Oct 2014 - fixes to the structure

!===================================================================================
 module allocateArray

   IMPLICIT NONE
   
   ! Set number of columns in input files
   integer, parameter:: ncolumnsSiteSelect=80	!SUEWS_SiteSelect.txt
   integer, parameter:: ncolumnsImpervious=16   !SUEWS_Impervious.txt
   integer, parameter:: ncolumnsPervious=27     !SUEWS_Pervious.txt
   integer, parameter:: ncolumnsWater=12	!SUEWS_Water.txt
   integer, parameter:: ncolumnsSnow=20		!SUEWS_Snow.txt
   integer, parameter:: ncolumnsSoil=8		!SUEWS_Soil.txt
   integer, parameter:: ncolumnsConductance=12	!SUEWS_Conductance.txt
   integer, parameter:: ncolumnsOHMCoefficients=4 	!SUEWS_OHMCoefficients.txt
   integer, parameter:: ncolumnsAnthropogenicHeat=11 	!SUEWS_AnthropogenicHeat.txt
   integer, parameter:: ncolumnsIrrigation=25		!SUEWS_Irrigation.txt
   integer, parameter:: ncolumnsProfiles=25 		!SUEWS_Profiles.txt
   integer, parameter:: ncolumnsWGWaterDist=10 		!SUEWS_WithinGridWaterDist.txt
  
   ! For input file headers               
   character(len=20),dimension(ncolumnsSiteSelect) :: HeaderSiteSelect_File     !Header for SiteSelect.txt                
   character(len=20),dimension(ncolumnsImpervious) :: HeaderImp_File  !Header for the impervious surface
   character(len=20),dimension(ncolumnsImpervious) :: HeaderImp_Reqd !Expected Header for the impervious surface 	
   character(len=20),dimension(ncolumnsPervious)   :: HeaderPer_File  !Header for the pervious surface
   character(len=20),dimension(ncolumnsPervious)   :: HeaderPer_Reqd  !Expected Header for the pervious surface 	
   character(len=20),dimension(ncolumnsWater)      :: HeaderWater_File  !Header for water surface
   character(len=20),dimension(ncolumnsWater)      :: HeaderWater_Reqd  !Expected Header for water surface 	
   character(len=20),dimension(ncolumnsSnow)       :: HeaderSnow_File  !Header for Snow surface
   character(len=20),dimension(ncolumnsSnow)       :: HeaderSnow_Reqd    !Expected Header for Snow surface 	
   character(len=20),dimension(ncolumnsSoil)       :: HeaderSoil_File  !Header for soils
   character(len=20),dimension(ncolumnsSoil)	   :: HeaderSoil_Reqd    !Expected Header for soils
   character(len=20),dimension(ncolumnsConductance):: HeaderCond_File !Header for conductances 
   character(len=20),dimension(ncolumnsConductance):: HeaderCond_Reqd  !Expected Header for conductances 
   character(len=20),dimension(ncolumnsOHMCoefficients)    :: HeaderOHMCoefficients_File  !Header for soils
   character(len=20),dimension(ncolumnsOHMCoefficients)	   :: HeaderOHMCoefficients_Reqd    !Expected Header for soils   
   character(len=20),dimension(ncolumnsAnthropogenicHeat)    :: HeaderAnthropogenicHeat_File  !Header for QF
   character(len=20),dimension(ncolumnsAnthropogenicHeat)    :: HeaderAnthropogenicHeat_Reqd    !Expected Header for QF 
   character(len=20),dimension(ncolumnsIrrigation):: HeaderIrrigation_File !Header for Irrigation
   character(len=20),dimension(ncolumnsIrrigation):: HeaderIrrigation_Reqd  !Expected Header for Irrigation
   character(len=20),dimension(ncolumnsProfiles):: HeaderProfiles_File !Header for Profiles
   character(len=20),dimension(ncolumnsProfiles):: HeaderProfiles_Reqd  !Expected Header for Profiles
   character(len=20),dimension(ncolumnsWGWaterDist):: HeaderWGWaterDist_File !Header for Profiles
   character(len=20),dimension(ncolumnsWGWaterDist):: HeaderWGWaterDist_Reqd  !Expected Header for Profiles
   
   real(kind(1d0)),dimension(:,:),allocatable::SiteSelect           !Matrix of SiteSelect.txt
   real(kind(1d0)),dimension(:,:),allocatable::Impervious_Coeff   !Coefficients for the impervious surfaces
   real(kind(1d0)),dimension(:,:),allocatable::Pervious_Coeff     !Coefficients for the pervious surfaces
   real(kind(1d0)),dimension(:,:),allocatable::Water_Coeff        !Coefficients for the water surface
   real(kind(1d0)),dimension(:,:),allocatable::Snow_Coeff	  !Coefficients for snow
   real(kind(1d0)),dimension(:,:),allocatable::Soil_Coeff	  !Coefficients for soil
   real(kind(1d0)),dimension(:,:),allocatable::Conductance_Coeff  !Coefficients for conductances
   real(kind(1d0)),dimension(:,:),allocatable::OHMCoefficients_Coeff   !Coefficients for OHMCoefficients 
   real(kind(1d0)),dimension(:,:),allocatable::AnthropogenicHeat_Coeff !Coefficients for AnthropogenicHeat
   real(kind(1d0)),dimension(:,:),allocatable::Irrigation_Coeff  !Coefficients for Irrigation
   real(kind(1d0)),dimension(:,:),allocatable::Profiles_Coeff  !Coefficients for Profiles
   real(kind(1d0)),dimension(:,:),allocatable::WGWaterDist_Coeff  !Coefficients for WithinGridWaterDist
      
      
   real(kind(1d0)),dimension(:,:,:),allocatable:: MetForcingData  !Meteorological forcing data matrix kept
                                                                  !in the program at once
   real(kind(1d0)),dimension(:,:,:),allocatable::ModelOutputData  !Same for tsep output data matrix
   real(kind(1d0)),dimension(:,:),  allocatable::ModelDailyState    !DailyState matrix
   real(kind(1d0)),dimension(:),    allocatable::DailyStateFirstOpen
   real(kind(1d0)),dimension(:,:),  allocatable::SurfaceChar        !Matrix for the surface characteristics
   
   real(kind(1d0)),dimension(:,:,:),allocatable:: TstepProfiles   !Array for hourly profiles interpolated to Tstep
   ! Columns
   integer:: cTP_EnUseWD  = 1,&
   	     cTP_EnUseWE  = 2,&
       	     cTP_WUManuWD = 3,&
       	     cTP_WUManuWE = 4,&
       	     cTP_WUAutoWD = 5,&
       	     cTP_WUAutoWE = 6,&
   	     cTP_SnowCWD  = 7,&
   	     cTP_SnowCWE  = 8
    
   real(kind(1d0)),dimension(:,:),allocatable:: AHProf_tstep 
   real(kind(1d0)),dimension(:,:),allocatable:: WUProfM_tstep, WUProfA_tstep

   !--------------------------------------------------------------
   ! Surface types
   integer:: PavSurf   = 1,&   !When all surfaces considered together (1-7)
   	     BldgSurf  = 2,&
             ConifSurf = 3,&
             DecidSurf = 4,&          
             GrassSurf = 5,&   !New surface classes: Grass = 5th/7 surfaces
             BSoilSurf = 6,&   !New surface classes: Bare soil = 6th/7 surfaces
             WaterSurf = 7,&
             ExcessSurf= 8,&   ! runoff or soil
             NSurfDoNotReceiveDrainage=0,& ! Number of surfaces do not receive drainage water - green roof
             ivConif = 1,&     !When only vegetated (now pervious) surfaces considered (1-4)
             ivDecid = 2,&
             ivGrass = 3,&     !New surface classes: Grass = 3rd/4 veg (pervious) surfaces
             ivBSoil = 4,&     !New surface classes: BSoil = 4th/4 veg (pervious) surfaces
            CreateAnnual,& 
            alldays=1,&
            NumberDailyVars=25  !Number of columns in daily/monthly output file
            
   integer, parameter:: nSurf=7      ! total number of surfaces
   integer, parameter:: NVegSurf=3   ! number of surfaces that are vegetated
   integer, parameter:: nSurfIncSnow=nsurf+1 !Number of surfaces + snow
   integer, parameter:: maxGrid=100  ! maximum number of surface grids -- could be much larger e.g 1000
   integer, parameter:: NCoeffOhm=4  ! number of coefficients per surface to use (4)
   integer, parameter:: Ndays =366   ! sg -- could read this in -- and allocate arrays
          
   integer:: numMonth=12             ! sg -- could read this in -- and allocate arrays
   integer,parameter:: sideG=90        !Number of horizontal grids to which snow melting is calculated

   !DailyState
   real (kind(1d0)),dimension(0:NDays, 5):: GDD ! growing degree days - for LAI
   
   ! check this may not be needed to go to -4 -- as 5 day mean only seems to be used in -OHM subroutine 
   real (kind(1d0)),dimension(-4:NDays, 6):: HDD ! heating dd - for QF
                                                !Matrix for heating and cooling degree days and Tair 
    										 !Last column contains the 5-day running mean								

   real (kind(1d0)),dimension(-4:NDays,NVegSurf)::LAI 	! LAI
      
   ! Replicate variables needed for DailyState, adding dimension to identify the grid, HCW 27 Nov 2014
   !! Could delete maxGrid and allocate these elsewhere once NumberOfGrids is known
   real (kind(1d0)),dimension(-4:NDays,6,maxGrid):: HDD_grids
   real (kind(1d0)),dimension( 0:NDays,5,maxGrid):: GDD_grids
   real (kind(1d0)),dimension(-4:NDays,NVegSurf,maxGrid):: LAI_grids
   real (kind(1d0)),dimension( 0:NDays,9,maxGrid):: WU_Day_grids
   real (kind(1d0)),dimension( 0:NDays,maxGrid):: albDec_grids
   real (kind(1d0)),dimension( 0:NDays,maxGrid):: DecidCap_grids
   real (kind(1d0)),dimension( 0:NDays,maxGrid):: porosity_grids
   
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
   
   !========NARP SPECIFIC PARAMETERS================================
   
   REAL(KIND(1D0))             :: NARP_LAT,NARP_LONG,NARP_YEAR,NARP_TZ,&
                                  NARP_ALB_SNOW,NARP_EMIS_SNOW,NARP_TRANS_SITE
   !REAL(KIND(1d0)),allocatable:: NARP_ALB(:),NARP_EMIS(:)                             
   ! check everywhere else 366 days
   REAL(KIND(1D0))             :: NARP_G(365)
   INTEGER                     :: NARP_NPERHOUR
   REAL(KIND(1D0)),ALLOCATABLE :: NARP_KDOWN_HR(:)
   
   REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292,RAD2DEG=57.29577951,&
                                  SIGMA_SB=5.67E-8
  
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
    
   !==================Outer programs================================= 
   character(len=15),dimension(2,maxGrid)::GridConnections  !List of different grid corrections
   real (kind(1d0)),dimension(maxGrid)::   GridConnectionsFrac   !Fraction of water moving between the different grids
   integer,dimension(maxGrid):: laiID   !Day of year of LAI from previous year
         
   real (kind(1d0)), dimension(0:ndays)::albDec    !Changing albedo for deciduous trees
   real (kind(1d0)), dimension(0:ndays)::DecidCap  !Max capacity of deciduous trees
   real (kind(1d0)), dimension(0:ndays):: porosity !Porosity of deciduous trees
   real (kind(1d0)), dimension(0:ndays,9):: WU_Day !Daily water use (total, automatic, manual) for EveTr,DecTr,Grass
   real (kind(1d0)):: AlbMin_dec=0.15,&
		      AlbMax_dec,&
                      CapMin_dec,&
                      CapMax_dec                                
   real (kind(1d0)), dimension(5,nsurf):: surf     !Surface information matrix defined in LUMPS_Initial.f90
          
   ! LUMPS vegetation phenology
   integer:: jj1,jj2,jj3,jj4, changed,changeInit,LAItype
   real (kind(1d0)):: kk1,kk2,dmax,dmin,dx1,dx2    ! used to calculate Vegetation phenology for LUMPS
   !real(kind(1d0)),dimension(0:23)::runT           ! running average T for the day   !Not used (HCW 27 Nov 2014) 
   !real(kind(1d0)),dimension(0:23)::runP           ! running total Precip for the day   !Not used (HCW 27 Nov 2014) 
   !real (kind(1d0))::avT_h, totP_h                 ! daily running average Temp, Total precip   !Not used (HCW 27 Nov 2014) 
   real(kind(1d0)),dimension(nvegsurf):: baseT     ! Base Temperature for growing degree days
   real(kind(1d0)),dimension(nvegsurf):: BaseTe    ! Base Temperature for senescence
   real(kind(1d0)),dimension(nvegsurf):: GDDFull   !Canopy growth parameters
   real(kind(1d0))::GDDmax,SDDMax                  ! limit across vegetation types
   real(kind(1d0)),dimension(nvegsurf):: LaiMin    ! minimum LAI
   real(kind(1d0)),dimension(nvegsurf):: LaiMax    ! maximum LAI
   !real(kind(1d0)):: MaxLaiMax                     ! maximum LAI of all  sg 12 nov 13 -- no longer used
   real(kind(1d0)),dimension(nvegsurf):: SDDFull   ! Vegetation parameters 
   real(kind(1d0)),dimension(nvegsurf+2):: MaxConductance      ! maximum conductance for surface resistance
                   ! +2 so can use surface position in array
                   
   real(kind(1d0)),dimension(4):: laiPower
                    
   ! used in water Use, anthropogenic heat and storage calculation                   
   integer,dimension(0:Ndays,3)::DayofWeek  ! day of week, month, season
   real (kind(1d0)), dimension(3)::WU_prof      !% of water use profile
   
   
   ! ==== Grid connections ==== ! Added HCW 14 Nov 2014
   integer,parameter:: NConns = 8 	!Number of grids for between-grid connections
   real(kind(1d0)), dimension(NConns) :: GridToFrac
   real(kind(1d0)), dimension(NConns) :: GridTo
       
   !======== Column numbers for SurfaceChar =============================    
   ! Columns 1:80 are the same as in SiteSelect.txt and defined below
   integer:: cc	 !Column counter
   integer,parameter:: ccEndSI=ncolumnsSiteSelect
   
   ! Applicable to each surface
   integer,dimension(nsurf):: c_Alb	  =(/(cc, cc=ccEndSI+ 0*nsurf+1,ccEndSI+ 0*nsurf+nsurf, 1)/)  !Albedo
   integer,dimension(nsurf):: c_Emis	  =(/(cc, cc=ccEndSI+ 1*nsurf+1,ccEndSI+ 1*nsurf+nsurf, 1)/)  !Emissivity
   integer,dimension(nsurf):: c_StorMin	  =(/(cc, cc=ccEndSI+ 2*nsurf+1,ccEndSI+ 2*nsurf+nsurf, 1)/)  !Min. storage capacity (canopy)
   integer,dimension(nsurf):: c_StorMax	  =(/(cc, cc=ccEndSI+ 3*nsurf+1,ccEndSI+ 3*nsurf+nsurf, 1)/)  !Max. storage capacity (canopy)
   integer,dimension(nsurf):: c_DrEq	  =(/(cc, cc=ccEndSI+ 4*nsurf+1,ccEndSI+ 4*nsurf+nsurf, 1)/)  !Drainage equation
   integer,dimension(nsurf):: c_DrCoef1	  =(/(cc, cc=ccEndSI+ 5*nsurf+1,ccEndSI+ 5*nsurf+nsurf, 1)/)  !Drainage coef. 1    
   integer,dimension(nsurf):: c_DrCoef2   =(/(cc, cc=ccEndSI+ 6*nsurf+1,ccEndSI+ 6*nsurf+nsurf, 1)/)  !Drainage coef. 2
   integer,dimension(nsurf):: c_SoilStCap  =(/(cc, cc=ccEndSI+ 7*nsurf+1,ccEndSI+ 7*nsurf+nsurf, 1)/)  !Soil storage capacity (below surface)	
   integer,dimension(nsurf):: c_SoilTCode  =(/(cc, cc=ccEndSI+ 8*nsurf+1,ccEndSI+ 8*nsurf+nsurf, 1)/)  !Soil type code
   ! N.B. not included in SUEWS_Water.txt
   integer,dimension(nsurf):: c_SnowLimPat =(/(cc, cc=ccEndSI+ 9*nsurf+1,ccEndSI+ 9*nsurf+nsurf, 1)/)  !Snow limit for patchiness
   ! N.B. currently only in SUEWS_Impervious.txt
   integer,dimension(nsurf):: c_SnowLimRem =(/(cc, cc=ccEndSI+10*nsurf+1,ccEndSI+10*nsurf+nsurf, 1)/)  !Snow limit for removal
   
   ! Find current column number	
   integer,parameter:: ccEndI = (ccEndSI+10*nsurf+nsurf)
   
   ! Applicable to pervious surfaces only
   integer,dimension(NVegSurf):: c_BaseT   =(/(cc, cc=ccEndI+ 0*nvegsurf+1,ccEndI+ 0*nvegsurf+nvegsurf, 1)/) !Base temp. for leaf-on
   integer,dimension(NVegSurf):: c_BaseTe  =(/(cc, cc=ccEndI+ 1*nvegsurf+1,ccEndI+ 1*nvegsurf+nvegsurf, 1)/) !Base temp. for leaf-off
   integer,dimension(NVegSurf):: c_GDDFull =(/(cc, cc=ccEndI+ 2*nvegsurf+1,ccEndI+ 2*nvegsurf+nvegsurf, 1)/) !GDD for full LAI
   integer,dimension(NVegSurf):: c_SDDFull =(/(cc, cc=ccEndI+ 3*nvegsurf+1,ccEndI+ 3*nvegsurf+nvegsurf, 1)/) !SDD for start of leaf-fall
   integer,dimension(NVegSurf):: c_LAIMin  =(/(cc, cc=ccEndI+ 4*nvegsurf+1,ccEndI+ 4*nvegsurf+nvegsurf, 1)/) !Min. LAI
   integer,dimension(NVegSurf):: c_LAIMax  =(/(cc, cc=ccEndI+ 5*nvegsurf+1,ccEndI+ 5*nvegsurf+nvegsurf, 1)/) !Max. LAI
   integer,dimension(NVegSurf):: c_GsMax   =(/(cc, cc=ccEndI+ 6*nvegsurf+1,ccEndI+ 6*nvegsurf+nvegsurf, 1)/) !Max. conductance
   integer,dimension(NVegSurf):: c_LAIEq   =(/(cc, cc=ccEndI+ 7*nvegsurf+1,ccEndI+ 7*nvegsurf+nvegsurf, 1)/) !LAI equation
   integer,dimension(NVegSurf):: c_LeafGP1 =(/(cc, cc=ccEndI+ 8*nvegsurf+1,ccEndI+ 8*nvegsurf+nvegsurf, 1)/) !Leaf growth power 1
   integer,dimension(NVegSurf):: c_LeafGP2 =(/(cc, cc=ccEndI+ 9*nvegsurf+1,ccEndI+ 9*nvegsurf+nvegsurf, 1)/) !Leaf growth power 2
   integer,dimension(NVegSurf):: c_LeafOP1 =(/(cc, cc=ccEndI+10*nvegsurf+1,ccEndI+10*nvegsurf+nvegsurf, 1)/) !Leaf-off power 1
   integer,dimension(NVegSurf):: c_LeafOP2 =(/(cc, cc=ccEndI+11*nvegsurf+1,ccEndI+11*nvegsurf+nvegsurf, 1)/) !Leaf-off power 2
     
   ! Find current column number	
   integer,parameter:: ccEndP = (ccEndI+11*nvegsurf+nvegsurf)
   
   ! Applicable to snow only
   integer:: c_SnowRMFactor = (ccEndP+ 1)
   integer:: c_SnowTMFactor = (ccEndP+ 2)  
   integer:: c_SnowAlbMin   = (ccEndP+ 3)      
   integer:: c_SnowAlbMax   = (ccEndP+ 4)
   integer:: c_SnowAlb      = (ccEndP+ 5)
   integer:: c_SnowEmis     = (ccEndP+ 6)      
   integer:: c_Snowtau_a    = (ccEndP+ 7)
   integer:: c_Snowtau_f    = (ccEndP+ 8)
   integer:: c_SnowPLimAlb  = (ccEndP+ 9)
   integer:: c_SnowSDMin    = (ccEndP+10)
   integer:: c_SnowSDMax    = (ccEndP+11)
   integer:: c_Snowtau_r    = (ccEndP+12)
   integer:: c_SnowCRWMin   = (ccEndP+13)
   integer:: c_SnowCRWMax   = (ccEndP+14)
   integer:: c_SnowPLimSnow = (ccEndP+15)
   
   ! Find current column number	
   integer,parameter:: ccEndSn = (ccEndP+15)
      
   ! Soil information
   integer,dimension(nsurf):: c_VolSMCap    = (/(cc, cc=ccEndSn+ 0*nsurf+1,ccEndSn+ 0*nsurf+nsurf, 1)/)  ! Volumetric SM capacity
   integer,dimension(nsurf):: c_KSat        = (/(cc, cc=ccEndSn+ 1*nsurf+1,ccEndSn+ 1*nsurf+nsurf, 1)/)  ! Saturated hydraulic conductivity
   integer,dimension(nsurf):: c_SoilDens    = (/(cc, cc=ccEndSn+ 2*nsurf+1,ccEndSn+ 2*nsurf+nsurf, 1)/)  ! Soil Density
   integer,dimension(nsurf):: c_SoilInfRate = (/(cc, cc=ccEndSn+ 3*nsurf+1,ccEndSn+ 3*nsurf+nsurf, 1)/)  ! Soil infiltration rate
   integer,dimension(nsurf):: c_ObsSMDepth  = (/(cc, cc=ccEndSn+ 4*nsurf+1,ccEndSn+ 4*nsurf+nsurf, 1)/)  ! Depth of SM obs
   integer,dimension(nsurf):: c_ObsSMMax    = (/(cc, cc=ccEndSn+ 5*nsurf+1,ccEndSn+ 5*nsurf+nsurf, 1)/)  ! Obs maximum SM [kg kg-1 OR m3 m-3]
   integer,dimension(nsurf):: c_ObsSNRFrac  = (/(cc, cc=ccEndSn+ 6*nsurf+1,ccEndSn+ 6*nsurf+nsurf, 1)/)  ! Obs fraction of soil without rocks
   
   ! Find current column number	
   integer,parameter:: ccEndSo = (ccEndSn+ 6*nsurf+nsurf)
   
   ! Surface conductance
   integer:: c_GsG1	= (ccEndSo+ 1)
   integer:: c_GsG2	= (ccEndSo+ 2)
   integer:: c_GsG3	= (ccEndSo+ 3)
   integer:: c_GsG4	= (ccEndSo+ 4)
   integer:: c_GsG5	= (ccEndSo+ 5)
   integer:: c_GsG6	= (ccEndSo+ 6)
   integer:: c_GsTH	= (ccEndSo+ 7)
   integer:: c_GsTL	= (ccEndSo+ 8)
   integer:: c_GsS1	= (ccEndSo+ 9)
   integer:: c_GsS2	= (ccEndSo+10)
   integer:: c_GsKmax	= (ccEndSo+11)
       
   ! Find current column number	
   integer,parameter:: ccEndGs = (ccEndSo+11)
   
   ! OHM codes
   integer,dimension(nsurfIncSnow):: c_OHMCode_SWet  =(/(cc, cc=ccEndGs+ 0*nsurfIncSnow+1,&
   								ccEndGs+ 0*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (summer wet)
   integer,dimension(nsurfIncSnow):: c_OHMCode_SDry  =(/(cc, cc=ccEndGs+ 1*nsurfIncSnow+1,&
   								ccEndGs+ 1*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (summer dry)
   integer,dimension(nsurfIncSnow):: c_OHMCode_WWet  =(/(cc, cc=ccEndGs+ 2*nsurfIncSnow+1,&
   								ccEndGs+ 2*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (winter wet)
   integer,dimension(nsurfIncSnow):: c_OHMCode_WDry  =(/(cc, cc=ccEndGs+ 3*nsurfIncSnow+1,&
   								ccEndGs+ 3*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM code (winter dry)
   integer,dimension(nsurfIncSnow):: c_a1_SWet       =(/(cc, cc=ccEndGs+ 4*nsurfIncSnow+1,&
   								ccEndGs+ 4*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (summer wet)
   integer,dimension(nsurfIncSnow):: c_a2_SWet       =(/(cc, cc=ccEndGs+ 5*nsurfIncSnow+1,&
   								ccEndGs+ 5*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (summer wet)
   integer,dimension(nsurfIncSnow):: c_a3_SWet	     =(/(cc, cc=ccEndGs+ 6*nsurfIncSnow+1,&
   								ccEndGs+ 6*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (summer wet)
   integer,dimension(nsurfIncSnow):: c_a1_SDry       =(/(cc, cc=ccEndGs+ 7*nsurfIncSnow+1,&
   								ccEndGs+ 7*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (summer dry)
   integer,dimension(nsurfIncSnow):: c_a2_SDry       =(/(cc, cc=ccEndGs+ 8*nsurfIncSnow+1,&
   								ccEndGs+ 8*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (summer dry)
   integer,dimension(nsurfIncSnow):: c_a3_SDry       =(/(cc, cc=ccEndGs+ 9*nsurfIncSnow+1,&
   								ccEndGs+ 9*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (summer dry)
   integer,dimension(nsurfIncSnow):: c_a1_WWet       =(/(cc, cc=ccEndGs+10*nsurfIncSnow+1,&
   								ccEndGs+10*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (winter wet)
   integer,dimension(nsurfIncSnow):: c_a2_WWet       =(/(cc, cc=ccEndGs+11*nsurfIncSnow+1,&
   								ccEndGs+11*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (winter wet)
   integer,dimension(nsurfIncSnow):: c_a3_WWet       =(/(cc, cc=ccEndGs+12*nsurfIncSnow+1,&
   								ccEndGs+12*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (winter wet)
   integer,dimension(nsurfIncSnow):: c_a1_WDry       =(/(cc, cc=ccEndGs+13*nsurfIncSnow+1,&
   								ccEndGs+13*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a1 (winter dry)
   integer,dimension(nsurfIncSnow):: c_a2_WDry       =(/(cc, cc=ccEndGs+14*nsurfIncSnow+1,&
   								ccEndGs+14*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a2 (winter dry)
   integer,dimension(nsurfIncSnow):: c_a3_WDry       =(/(cc, cc=ccEndGs+15*nsurfIncSnow+1,&
   								ccEndGs+15*nsurfIncSnow+nsurfIncSnow, 1)/)  !OHM a3 (winter dry)   

   ! Find current column number	
   integer,parameter:: ccEndO = (ccEndGs+15*nsurfIncSnow+nsurfIncSnow)
   
   ! Anthropogenic heat
   integer :: c_BaseTHDD  = (ccEndO+ 1)
   integer :: c_QF_A1	  = (ccEndO+ 2)
   integer :: c_QF_B1     = (ccEndO+ 3)
   integer :: c_QF_C1     = (ccEndO+ 4)
   integer :: c_QF_A2 	  = (ccEndO+ 5)
   integer :: c_QF_B2	  = (ccEndO+ 6)
   integer :: c_QF_C2 	  = (ccEndO+ 7)
   integer :: c_AHMin     = (ccEndO+ 8)
   integer :: c_AHSlope   = (ccEndO+ 9)
   integer :: c_TCritic   = (ccEndO+10)
 
   ! Find current column number	
   integer,parameter:: ccEndA = (ccEndO+10)
    
   ! Irrigation 
   integer :: c_IeStart	  = (ccEndA+ 1)
   integer :: c_IeEnd	  = (ccEndA+ 2)
   integer :: c_IntWU	  = (ccEndA+ 3)
   integer :: c_Faut	  = (ccEndA+ 4)
   integer,dimension(3):: c_Ie_a      = (/(cc, cc=ccEndA+4+ 0*3+1, ccEndA+4 + 0*3+3, 1)/)  ! Automatic irrigation coeffs
   integer,dimension(3):: c_Ie_m      = (/(cc, cc=ccEndA+4+ 1*3+1, ccEndA+4 + 1*3+3, 1)/)  ! Manual irrigation coeffs
   integer,dimension(7):: c_DayWat    = (/(cc, cc=ccEndA+10+ 0*7+1,ccEndA+10+ 0*7+7, 1)/)  ! Irrigation allowed on each day
   integer,dimension(7):: c_DayWatPer = (/(cc, cc=ccEndA+10+ 1*7+1,ccEndA+10+ 1*7+7, 1)/)  ! Fraction properties using irrigation allowed on each day
   
   ! Find current column number	
   integer,parameter:: ccEndIr = (ccEndA+10+ 1*7+7)
       
   ! Hourly profiles
   integer,dimension(24):: c_HrProfEnUseWD  = (/(cc, cc=ccEndIr+ 0*24+1, ccEndIr+ 0*24+24, 1)/)  ! Energy use, weekdays
   integer,dimension(24):: c_HrProfEnUseWE  = (/(cc, cc=ccEndIr+ 1*24+1, ccEndIr+ 1*24+24, 1)/)  ! Energy use, weekends
   integer,dimension(24):: c_HrProfWUManuWD = (/(cc, cc=ccEndIr+ 2*24+1, ccEndIr+ 2*24+24, 1)/)  ! Water use, manual, weekdays
   integer,dimension(24):: c_HrProfWUManuWE = (/(cc, cc=ccEndIr+ 3*24+1, ccEndIr+ 3*24+24, 1)/)  ! Water use, manual, weekends
   integer,dimension(24):: c_HrProfWUAutoWD = (/(cc, cc=ccEndIr+ 4*24+1, ccEndIr+ 4*24+24, 1)/)  ! Water use, automatic, weekdays
   integer,dimension(24):: c_HrProfWUAutoWE = (/(cc, cc=ccEndIr+ 5*24+1, ccEndIr+ 5*24+24, 1)/)  ! Water use, automatic, weekends
   integer,dimension(24):: c_HrProfSnowCWD  = (/(cc, cc=ccEndIr+ 6*24+1, ccEndIr+ 6*24+24, 1)/)  ! Snow clearing, weekdays
   integer,dimension(24):: c_HrProfSnowCWE  = (/(cc, cc=ccEndIr+ 7*24+1, ccEndIr+ 7*24+24, 1)/)  ! Snow clearing, weekends
 
   ! Find current column number	
   integer,parameter:: ccEndPr = (ccEndIr+ 5*24+24)
       
   ! Within-grid water distribution (for each surface)
   integer,dimension(nsurf):: c_WGToPaved = (/(cc, cc=ccEndPr+ 0*nsurf+1,ccEndPr+ 0*nsurf+nsurf, 1)/) !Water dist to Paved
   integer,dimension(nsurf):: c_WGToBuilt = (/(cc, cc=ccEndPr+ 1*nsurf+1,ccEndPr+ 1*nsurf+nsurf, 1)/) !Water dist to Built
   integer,dimension(nsurf):: c_WGToEveTr = (/(cc, cc=ccEndPr+ 2*nsurf+1,ccEndPr+ 2*nsurf+nsurf, 1)/) !Water dist to EveTr
   integer,dimension(nsurf):: c_WGToDecTr = (/(cc, cc=ccEndPr+ 3*nsurf+1,ccEndPr+ 3*nsurf+nsurf, 1)/) !Water dist to DecTr
   integer,dimension(nsurf):: c_WGToGrass = (/(cc, cc=ccEndPr+ 4*nsurf+1,ccEndPr+ 4*nsurf+nsurf, 1)/) !Water dist to Grass
   integer,dimension(nsurf):: c_WGToBSoil = (/(cc, cc=ccEndPr+ 5*nsurf+1,ccEndPr+ 5*nsurf+nsurf, 1)/) !Water dist to BSoil
   integer,dimension(nsurf):: c_WGToWater = (/(cc, cc=ccEndPr+ 6*nsurf+1,ccEndPr+ 6*nsurf+nsurf, 1)/) !Water dist to Water
   integer,dimension(nsurf):: c_WGToRunoff    = (/(cc, cc=ccEndPr+ 7*nsurf+1,ccEndPr+ 7*nsurf+nsurf, 1)/) !Water dist to runoff 
   integer,dimension(nsurf):: c_WGToSoilStore = (/(cc, cc=ccEndPr+ 8*nsurf+1,ccEndPr+ 8*nsurf+nsurf, 1)/) !Water dist to sub-surface soil
 
   !Last column number for SurfaceChar array
   integer,parameter:: MaxNCols_c = ccEndPr+ 8*nsurf+nsurf
 
   !----------------------------------------------------------------------
   ! ---- For ModelOutputData ----
   ! Applicable to each surface
   integer,parameter:: ccMOD = 32
   integer,dimension(nsurf):: cMOD_State          =(/(cc, cc=ccMOD+ 0*nsurf+1,ccMOD+ 0*nsurf+nsurf, 1)/)  !Above ground state
   integer,dimension(nsurf):: cMOD_SoilState      =(/(cc, cc=ccMOD+ 1*nsurf+1,ccMOD+ 1*nsurf+nsurf, 1)/)  !Below ground state (soil store)
   integer,dimension(nsurf):: cMOD_SnowWaterState =(/(cc, cc=ccMOD+ 2*nsurf+1,ccMOD+ 2*nsurf+nsurf, 1)/)  !Liquid (melted) water
   integer,dimension(nsurf):: cMOD_SnowPack       =(/(cc, cc=ccMOD+ 3*nsurf+1,ccMOD+ 3*nsurf+nsurf, 1)/)  !SWE
   integer,dimension(nsurf):: cMOD_SnowFrac       =(/(cc, cc=ccMOD+ 4*nsurf+1,ccMOD+ 4*nsurf+nsurf, 1)/)  !Snow fraction
   !integer,dimension(nsurf):: cMOD_SnowDens      =(/(cc, cc=ccMOD+ 5*nsurf+1,ccMOD+ 5*nsurf+nsurf, 1)/)  !Snow density
   !SnowDens in ModelDailyState instead - why??
   
   !Last column number for ModelOutputData array
   integer,parameter:: MaxNCols_cMOD = ccMOD+ 5*nsurf+nsurf
   
   !----------------------------------------------------------------------
  
   !----------------------------------------------------------------------
   ! ---- For ModelDailyState ----
   ! Applicable to each surface
   integer,parameter:: ccMDS = 30
   integer,dimension(nsurf):: cMDS_SnowDens       =(/(cc, cc=ccMDS+ 0*nsurf+1,ccMDS+ 0*nsurf+nsurf, 1)/)  !Snow density
   
   !Last column number for ModelDailyState array
   integer,parameter:: MaxNCols_cMDS = ccMDS+ 0*nsurf+nsurf
   
   !----------------------------------------------------------------------
      
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
!New module for model initialization

 !-----------------------------------------------------------------------------
  MODULE Initial
 
     IMPLICIT NONE

     integer::FirstYear,&          !First year
     	      LastYear,&           !Last year
     	      FirstGrid,&	   !First grid (as in SiteSelect)
     	      LastGrid,&	   !Last grid (as in SiteSelect)
     	      NumberOfGrids,&      !Number of grids
              GridCounter,&        !Counter for grids (i.e. from 1 to NumberOfGrids)
              ReadlinesMetdata,&   !Number of lines in each block of met data read at once
              skippedLines,& 	   !Number of lines to skip over before reading each block of met data
              nlinesMetdata,&      !Total number of lines in met forcing file
              nlinesSiteSelect,&   !Number of lines in SUEWS_SiteSelect.txt
              nlinesImpervious,&   !Number of lines in SUEWS_Impervious.txt
              nlinesPervious,&     !Number of lines in SUEWS_Pervious.txt
              nlinesWater,&   	   !Number of lines in SUEWS_Water.txt
              nlinesSnow,&		
              nlinesSoil,&
              nlinesConductance,&
              nlinesOHMCoefficients,&
              nlinesAnthropogenicHeat,&
              nlinesIrrigation,&
              nlinesProfiles,&
              nlinesWGWaterDist,&
              nlines,&             !Number of lines in different files
              iv5		   !Counter for code matching

     character (len=150):: FileMet !Meteorological input filename. Change location

  END MODULE Initial
 !-----------------------------------------------------------------------------

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
     real(kind(1d0)), dimension(0:23,2):: snowProf  ! Timing of snow removal (0 or 1) Hourly, WD/WE
     
     integer::SnowFractionChoice   !Choice how fraction of snow is calculated
 
 end module snowMod

 !===================================================================================
 Module defaultNotUsed
 	implicit none
 	real (kind(1d0)):: notUsed=-55.55,reall,NAN=-999,pNAN=999
 	integer:: notUsedI=-55, ios_out,errorChoice  !errorChoice defines if the problemfile is opened for the first time
 end Module defaultNotUsed
 
 !===================================================================================
 
 module data_in
 
 IMPLICIT NONE
 
 ! In alphabetical order
 real (kind(1d0)):: AH_MIN,&    !Minimum anthropogenic heat flux (AnthropHeatChoice = 1)
                    AH_SLOPE,&  !Slope of the antrhropogenic heat flux calculation (AnthropHeatChoice = 1)
                    alpha_qhqe,& !Alpha parameter used in LUMPS QH and QE calculations
                    avdens,&    !Average air density
                    avkdn,&     !Average downwelling shortwave radiation
                    avrh,&      !Average relative humidity
                    avts,&      !Average surface temperature
                    avu1,&      !Average wind speed
                    azimuth,&   !Sun azimuth in degrees
                    BaseTHDD,&  !Base temperature for QF               
                    defaultQf,& !Default anthropogenic heat flux
                    defaultQs,& !Default storage heat flux
                    E_mod,&     !Modelled latent heat flux with LUMPS
                    emis_snow,& !Emissivity of snow
                    fcld,&      !Cloud fraction modelled
                    fcld_obs,&  !Cloud fraction observed
                    h_mod,&     !Modelled sensible heat flux with LUMPS
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
                    qn1_S,&     !Total net all-wave radiation for the snowpack
                    qn1_SF,&    !Total net all-wave radiation for the snowfree surface
                    qs,&        !Observed storage heat flux
                    snow,&      !snow cover
                    snow_obs,&  !Observed snow cover
                    T_CRITIC,& !Critical temperature
                    Temp_C,&    !Air temperature
                    timezone,&  !Timezone (GMT=0)
                    trans_site,&  !Atmospheric transmittivity
                    tsurf,&   !Surface temperature
                    wdir,&      ! Wind direction
                    wu_m3,&     !Water use provided in met forcing file [m3]
                    xsmd,&      !Measured soil moisture deficit
                    year,&      !Year of the measurements
                    zenith_deg  !Sun zenith angle in degrees


              
  real (kind(1d0)),dimension(366,25)::day,season,month,yr_tot,all_tot !daily matrixes
  integer,dimension(:,:), allocatable:: dataMet1              !Meteorological input matrix
  real(kind(1d0)),dimension(:,:), allocatable:: dataMet2              !Meteorological input matrix
  real(kind(1d0)),dimension(:,:,:), allocatable:: dataOut             !Main output matrix
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
           ldown_option,&         !What parameterization is used for downward longwave radiation 1-2-3
           lfnout,&               !Error Output write units
           lfnoutC,&              !Clean output write units
           lfnOld,&
           lfnSAHP,&              !Number of SAHP file
           NARPOutput,&           !Defines if radiation components are separatley printed out
           netradiationchoice,&   !Is net all-wave radiation modeled (=2) or measured (=1)
           qschoice,&             !Defines if QS is calculated
           SkipHeaderSiteInfo,&   !Number of header lines to skip in SiteInfo files 
           SkipHeaderMet,&        !Number of header lines to skip in met file input
           SNOWuse,&
           SOLWEIGout,&           !Calculates Tmrt and other fluxes on a grid, FL
           write5min,&            !Defines if 5-min output is printed
           writedailyState=1
                   
      character (len=150)::FileInputPath,&   !Filepath for input file
                           FileOutputPath,&  !Filepath for output file
                           fileout,&         !Output file name   
                           FileErrorInf,&    !Output error file name
                           filechoices,&     !File with output of run characteristics
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
      character (len=90)::progname='SUEWS V2014c'  !!!!<<<<<<<<<<<<<<<<<<
        
      logical:: finish,once
 
      real(kind(1d0)),dimension(2)::Qf_A,Qf_B,Qf_C !Qf coefficients
      
      integer,dimension(2)::DayLightSavingDay    !The date when it is changed to daylight saving (in DOY)
                          
      real(kind(1d0)), dimension(0:23,2):: AHPROF !Anthropogenic heat profiles for (1)weekdays / (2)weekends
 
      integer::nCBLstep  !number of time steps of Runge-kutta methods in one hour
         
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
   integer :: iy,&            !Year
              id,&            !Day of year
              it,&            !Hour
              imin,&          !Minutes
              iostat_var,&      !File status from reading data (should not be here)
              lastTimeofDAY=23, firstTimeofDay=0,&
              DLS                            !- day lightsavings =1 + 1h) =0  
 
             ! check what lasttime of day should be
   real (kind(1d0)):: dectime !Time is decimals
   real (kind(1d0)):: halftimestep !in decimal time based on interval 
   real (kind(1d0)):: tstepcount ! count number of timesteps in this day
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
 module Thresh
   real(kind(1d0)) :: IPThreshold_mmhr = 10   !Threshold for intense precipitation [mm hr-1]
   
 end module Thresh
 
 
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
   
   integer:: tstep,&    !Timestep [s] at which the model is run (set in RunControl)
             nsh,&      !Number of timesteps per hour
             interval   !Number of seconds in an hour [s] (now set in OverallRunControl)
   
   real(kind(1d0)):: nsh_real,&   !nsh cast as a real for use in calculations                   
                     tstep_real   !tstep cast as a real for use in calculations                   
             
   !Options for model setup (switches, etc) mainly set in RunControl 
   integer:: LAIcalcYes,&        !Use observed(0) or modelled(1) LAI  !Set to 1 in OverallRunControl
             ity,&               !Evaporation calculated according to Rutter(1) or Shuttleworth(2)  !Set to 2 in OverallRunControl
             z0_method,&         !Defines method for calculating z0 & zd (set in RunControl)
             WU_choice,&         !Use modelled(0) or observed(1) water use (set in RunControl)
             SMD_choice,&        !Defines if soil moisture is modelled(0) or observed(1,2) (set in RunControl)
             StabilityMethod,&   !Defines stability functions used (set in RunControl)
             RoughLen_heat       !Defines method for calculating roughness length for heat (set in RunControl)
                     
   integer:: in                     
   integer:: is      !Integer to count over surface types

   !These are variables which currently have been removed from SuesInput.nml     
   integer::AerodynamicResistanceMethod=2 !The method used to calculate aerodynamic resistance
   
   integer::Ie_start,&   !Starting time of water use (DOY)
            Ie_end       !Ending time of water use (DOY)     
       
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
                      rst,&      !Flag in SUEWS_Evap (gets set to 1 if surface dry; 0 if surface wet)
                      qeph       !Latent heat flux (W m^-2)
                      
 
 !Water use related variables
  real (kind(1d0)):: ext_wu,&         !External water use for the model timestep [mm] (over whole study area)
                     Faut,&           !Fraction of irrigated area using automatic irrigation
                     int_wu,&         !Internal water use for the model timestep [mm] (over whole study area)
                     IrrFractionTrees,& !Fraction of the surafce area of irrigated trees
                     IrrTrees,&         !Surface area fraction of irrigated trees
                     IrrFracConif,&	!Fraction of evergreen trees which are irrigated
		     IrrFracDecid,&	!Fraction of deciduous trees which are irrigated
		     IrrFracGrass,&	!Fraction of grass which is irrigated
                     InternalWaterUse_h !Internal water use [mm h-1]
  
 ! 7 - number of days in week                   
  real(kind(1d0)),dimension(7)::DayWatPer,&  !% of houses following daily water
                                  DayWat       !Days of watering allowed
  real(kind(1d0)),dimension(0:23,2):: WUProfM,&   !Hourly profiles for water use (manual irrigation) 
   				      WUProfA   !Hourly profiles for water use (automatic irrigation) 
 
 
  real (kind(1d0)),dimension(3)::Ie_a,&
                                  Ie_m           
     
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

   real (Kind(1d0))::LAIinitialEveTr,&
                     LAIinitialDecTr,&
                     LAIinitialGrass,&
                     porosity0,&
                     DecidCap0,&
                     albDec0,&
                     Temp_C0,&
                     GDD_1_0,&
                     GDD_2_0,&
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
                     SnowPackWater
                     
     integer::ID_Prev

 end Module InitialCond

!-------------------------------------------------
!New modules for the column numbers


!-------------------------------------------------------------------------
 Module ColNamesModelDailyState
 
   IMPLICIT NONE
    
   !========== Columns for ModelDailyState array =========================
 
   integer:: cMDS_id_prev	= 3,& 
   	     cMDS_HDD1		= 4,&
             cMDS_HDD2		= 5,&
             cMDS_TempC		= 6,&
             cMDS_TempCRM  	= 7,&
             cMDS_Precip	= 8,&
             cMDS_DaysSinceRain	= 9,&
             cMDS_TempCOld1	=10,&
	     cMDS_TempCOld2	=11,&
             cMDS_TempCOld3	=12,&
                        
             cMDS_GDDMin	=13,&
             cMDS_GDDMax	=14,&
             cMDS_GDD1_0	=15,&
             cMDS_GDD2_0	=16,&
             cMDS_LAIInitialEveTr =17,&
             cMDS_LAIInitialDecTr =18,&
             cMDS_LAIInitialGrass =19,&
                          
             cMDS_porosity	=20,&
             cMDS_albDec	=21,&
	     cMDS_DecidCap	=22,&
             cMDS_CumSnowfall	=23,&

	     cMDS_LAIEveTr	=24,&
	     cMDS_LAIDecTr	=25,&
	     cMDS_LAIGrass	=26
	     
 end Module ColNamesModelDailyState


!-------------------------------------------------------------------------
 Module ColNamesInputFiles
 
   IMPLICIT NONE
   
   ! Column names and numbers must match the input files
   
   !========== Columns for SUEWS_SiteSelect.txt ==========================
   ! Columns 1:80 are the same for SurfaceChar
   integer::c_Grid	    = 1,&
            c_Year	    = 2,&
            c_StartDLS	    = 3,&
            c_EndDLS	    = 4,&
            ! Site info
            c_lat	    = 5,&
            c_lng	    = 6,&
            c_Area	    = 7,&
            c_Alt	    = 8,&
            ! Time info
            c_id	    = 9,&
            c_it	    =10,&
            c_imin	    =11,&
            ! Surface fractions
            c_FrPaved	    =12,&
            c_FrBuilt	    =13,&
            c_FrEveTr	    =14,&
            c_FrDecTr	    =15,&
            c_FrGrass	    =16,&
            c_FrBSoil	    =17,&
            c_FrWater	    =18,&
            ! Irrigated fractions
            c_IrrEveTrFrac  =19,&
            c_IrrDecTrFrac  =20,&
            c_IrrGrassFrac  =21,&
            ! Height information
            c_HBuilt	    =22,&
            c_HEveTr	    =23,&
            c_HDecTr	    =24,&
            c_z0m	    =25,&
            c_zdm	    =26,&
            c_FAIBuilt	    =27,&
            c_FAIEveTr	    =28,&
            c_FAIDecTr	    =29,&
            ! Population
	    c_PopDensDay    =30,&
            c_PopDensNight  =31,&
	    ! Codes for different surfaces
            c_PavedCode	    =32,&	! Links characteristics in SUEWS_Impervious.txt
            c_BuiltCode	    =33,&	! Links characteristics in SUEWS_Impervious.txt
            c_EveTrCode	    =34,&	! Links characteristics in SUEWS_Pervious.txt
            c_DecTrCode	    =35,&  	! Links characteristics in SUEWS_Pervious.txt
            c_GrassCode	    =36,&   	! Links characteristics in SUEWS_Pervious.txt
            c_BSoilCode	    =37,&	! Links characteristics in SUEWS_Pervious.txt
            c_WaterCode	    =38,&       ! Links characteristics in SUEWS_Water.txt
	    ! LUMPS info
	    c_LUMPSDr	    =39,&
	    c_LUMPSCover    =40,&
	    c_LUMPSMaxRes   =41,&
	    ! NARP info
	    c_NARPTrans	    =42,&
	    ! Code for conductances
	    c_CondCode	    =43,&       ! Links characteristics in SUEWS_Conductance.txt
	    ! Code for snow		
            c_SnowCode      =44,&  	! Links characteristics in SUEWS_Snow.txt
            ! Codes for human impacts on energy, water and snow
            c_SnowProfWD    =45,&	! Snow-clearing profile in SUEWS_Profile.txt (weekdays)
	    c_SnowProfWE    =46,&	! Snow-clearing profile in SUEWS_Profile.txt (weekends)
	    c_QFCode        =47,&     	! Links anthropogenic heat info in SUEWS_AnthropogenicHeat.txt
	    c_EnProfWD	    =48,&	! Links to energy-use profile in SUEWS_Profile.txt (weekdays)
	    c_EnProfWE 	    =49,&	! Links to energy-use profile in SUEWS_Profile.txt (weekends)
	    c_IrrCode	    =50,&     	! Links irrigation info in SUEWS_Irrigation.txt
	    c_WProfManuWD   =51,&	! Links to water-use profile in SUEWS_Profile.txt (manual irrigation, weekdays)
	    c_WProfManuWE   =52,&	! Links to water-use profile in SUEWS_Profile.txt (manual irrigation, weekends)
	    c_WProfAutoWD   =53,&	! Links to water-use profile in SUEWS_Profile.txt (automatic irrigation, weekdays)
	    c_WProfAutoWE   =54,&	! Links to water-use profile in SUEWS_Profile.txt (automatic irrigation, weekends)
	    ! Flow information 	
	    c_FlowChange    =55,&	! Difference in input & output flows for water surface
	    c_RunoffToWater =56,&  	! Fraction of above-ground runoff flowing to water surface
	    c_PipeCapacity  =57,&	! Pipe capacity [mm]	
 	    ! Runoff (to 8 adjacent grids)		
            c_GridConnection1of8 =58,&
            c_Fraction1of8 	 =59,&
            c_GridConnection2of8 =60,&
            c_Fraction2of8 	 =61,&
            c_GridConnection3of8 =62,&
            c_Fraction3of8 	 =63,&
            c_GridConnection4of8 =64,&
            c_Fraction4of8 	 =65,&
            c_GridConnection5of8 =66,&
            c_Fraction5of8 	 =67,&
            c_GridConnection6of8 =68,&
            c_Fraction6of8 	 =69,&
            c_GridConnection7of8 =70,&
            c_Fraction7of8 	 =71,&
            c_GridConnection8of8 =72,&
            c_Fraction8of8 	 =73,&            
 	    ! Runoff within grid (for each surface type)
 	    c_WGPavedCode   =74,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
 	    c_WGBuiltCode   =75,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
 	    c_WGEveTrCode   =76,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
 	    c_WGDecTrCode   =77,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
 	    c_WGGrassCode   =78,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
     	    c_WGBSoilCode   =79,& 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		
 	    c_WGWaterCode   =80 	! Links to SUEWS_WaterDistibuteWithinGrid.txt		 	     	    
 	 
   !========== Columns for SUEWS_Impervious.txt ==========================
   integer :: ci_Code	      =  1,&
   	      ci_Alb 	      =  2,&
   	      ci_Emis         =  3,&
	      ci_StorMin      =  4,&
	      ci_StorMax      =  5,&
	      ci_DrEq         =  6,&
	      ci_DrCoef1      =  7,&
	      ci_DrCoef2      =  8,&
	      ci_SoilStCap    =  9,&
	      ci_SoilTCode    = 10,&
	      ci_SnowLimPat   = 11,&
              ci_SnowLimRem   = 12,&
	      ci_OHMCode_SWet = 13,&
	      ci_OHMCode_SDry = 14,&
	      ci_OHMCode_WWet = 15,&
	      ci_OHMCode_WDry = 16              
            
   !========== Columns for SUEWS_Pervious.txt ============================
   integer :: cp_Code	    =  1,&
   	      cp_Alb 	    =  2,&
   	      cp_Emis       =  3,&
	      cp_StorMin    =  4,&
	      cp_StorMax    =  5,&
	      cp_DrEq       =  6,&
	      cp_DrCoef1    =  7,&
	      cp_DrCoef2    =  8,&
	      cp_SoilStCap  =  9,&
	      cp_SoilTCode  = 10,&
	      cp_SnowLimPat = 11,&
	      cp_BaseT 	    = 12,&
	      cp_BaseTe	    = 13,&
	      cp_GDDFull    = 14,&
	      cp_SDDFull    = 15,&
	      cp_LAIMin	    = 16,&
	      cp_LAIMax	    = 17,&
	      cp_GsMax	    = 18,&
	      cp_LAIEq      = 19,&
	      cp_LeafGP1    = 20,&
	      cp_LeafGP2    = 21,&
	      cp_LeafOP1    = 22,&
	      cp_LeafOP2    = 23,&
	      cp_OHMCode_SWet = 24,&
	      cp_OHMCode_SDry = 25,&
	      cp_OHMCode_WWet = 26,&
	      cp_OHMCode_WDry = 27	      
   
   !========== Columns for SUEWS_Water.txt ===============================
   integer :: cw_Code	    =  1,&
   	      cw_Alb 	    =  2,&
   	      cw_Emis       =  3,&
	      cw_StorMin    =  4,&
	      cw_StorMax    =  5,&
	      cw_DrEq       =  6,&
	      cw_DrCoef1    =  7,&
	      cw_DrCoef2    =  8,&
	      cw_OHMCode_SWet = 9,&
	      cw_OHMCode_SDry = 10,&
	      cw_OHMCode_WWet = 11,&
	      cw_OHMCode_WDry = 12
	      
   !========== Columns for SUEWS_Snow.txt ================================	      
   integer :: cs_Code	      =  1,&
   	      cs_SnowRMFactor =  2,&
   	      cs_SnowTMFactor =  3,&
   	      cs_SnowAlbMin   =  4,&
   	      cs_SnowAlbMax   =  5,&
   	      cs_SnowAlb      =  6,&
   	      cs_SnowEmis     =  7,&
   	      cs_Snowtau_a    =  8,&
   	      cs_Snowtau_f    =  9,&
	      cs_SnowPLimAlb  = 10,&   	      
	      cs_SnowSDMin    = 11,&   	      
	      cs_SnowSDMax    = 12,&   	      
	      cs_Snowtau_r    = 13,&
      	      cs_SnowCRWMin   = 14,&
   	      cs_SnowCRWMax   = 15,&
   	      cs_SnowPLimSnow = 16,&     
	      cs_OHMCode_SWet = 17,&
	      cs_OHMCode_SDry = 18,&
	      cs_OHMCode_WWet = 19,&
	      cs_OHMCode_WDry = 20
	      
   !========== Columns for SUEWS_Soil.txt ================================
   integer :: cSo_Code        =  1,&
   	      cSo_VolSmCap    =  2,&
   	      cSo_KSat 	      =  3,&
   	      cSo_SoilDens    =  4,&
	      cSo_SoilInfRate =  5,&
	      cSo_ObsSMDepth  =  6,&
	      cSo_ObsSMMax    =  7,&
	      cSo_ObsSNRFrac  =  8
	      
   !========== Columns for SUEWS_Conductance.txt =========================
   integer :: cc_Code   =  1,&
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
   integer :: cO_Code   =  1,&
   	      cO_a1	=  2,&
   	      cO_a2   	=  3,&
   	      cO_a3   	=  4   
   	      
   !========== Columns for SUEWS_AnthropogenicHeat.txt ===================
   integer :: cA_Code 	  =  1,&
   	      cA_BaseTHDD =  2,&
   	      cA_QF_A1	  =  3,&   !Weekday
   	      cA_QF_B1    =  4,&   !Weekday
   	      cA_QF_C1    =  5,&   !Weekday
   	      cA_QF_A2 	  =  6,&   !Weekend 
   	      cA_QF_B2	  =  7,&   !Weekend
   	      cA_QF_C2 	  =  8,&   !Weekend
   	      cA_AHMin    =  9,&
   	      cA_AHSlope  = 10,&
   	      cA_TCritic  = 11

   !========== Columns for SUEWS_Irrigation.txt ==========================   	      
 	       	      
   integer ::	cIr_Code	=  1,&
   		cIr_IeStart 	=  2,&
   		cIr_IeEnd	=  3,&
   		cIr_IntWU	=  4,&
   		cIr_Faut	=  5,&
   		cIr_Ie_a1    	=  6,&
		cIr_Ie_a2    	=  7,&   				
		cIr_Ie_a3    	=  8,&
		cIr_Ie_m1    	=  9,&
		cIr_Ie_m2    	= 10,&
		cIr_Ie_m3    	= 11,&
   		cIr_DayWat1	= 12,&		
   		cIr_DayWat2	= 13,&		
   		cIr_DayWat3	= 14,&		
   		cIr_DayWat4	= 15,&		
   		cIr_DayWat5	= 16,&		
   		cIr_DayWat6	= 17,&		
   		cIr_DayWat7	= 18,&		
   		cIr_DayWatPer1	= 19,&		
   		cIr_DayWatPer2	= 20,&		
   		cIr_DayWatPer3	= 21,&		
   		cIr_DayWatPer4	= 22,&		
   		cIr_DayWatPer5	= 23,&		
   		cIr_DayWatPer6	= 24,&		
   		cIr_DayWatPer7	= 25		
   		
   !========== Columns for SUEWS_Profile.txt =============================   	         		
   
   integer:: cc	 !Column counter
   
   integer:: cPr_Code = 1
   integer,dimension(24):: cPr_Hours = (/(cc, cc=2,25, 1)/)  ! Hourly profile data

   !========== Columns for SUEWS_WithinGridWaterDist.txt =================
   
   integer:: 	cWG_Code 	= 1,&
   		cWG_ToPaved	= 2,&
   		cWG_ToBuilt	= 3,&
   		cWG_ToEveTr	= 4,&
   		cWG_ToDecTr	= 5,&
   		cWG_ToGrass	= 6,&
   		cWG_ToBSoil	= 7,&
   		cWG_ToWater	= 8,&
   		cWG_ToRunoff	= 9,&
   		cWG_ToSoilStore	= 10

 end Module ColNamesInputFiles



 !C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95C:\Users\sue\Dropbox\BLUEWS\2012av\LUMPS_Module_constants_v6_0.f95
 
 
 
 
 
 
 