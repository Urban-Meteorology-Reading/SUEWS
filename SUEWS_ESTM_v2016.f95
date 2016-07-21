!===============================================================================
SUBROUTINE ESTM_v2016(QSnet,Gridiv,ir)
  ! HCW Questions:
  !                - should TFloor be set in namelist instead of hard-coded here?
  !                - zref used for radiation calculation and fair is set to 2*BldgH here. For compatibility with the rest of the
  !                  SUEWS model, should this be the (wind speed) measurement height z specified in RunControl.nml?
  !                - In SUEWS_translate, fwall=AreaWall/SurfaceArea. Is this correct?
  !                - If froof=1 (i.e. whole grid is building), is HW=0 correct?  
  !                - Then is an IF(Fground ==0) option needed?  
  !                - alb_wall=0.23 and em_wall=0.9 are set in LUMPS_module_constants. Shouldn't these be provided as input?
  !                - Do the LUP calculations here need to be compatible with those in LUMPS_NARP?
  !                - File opening rewritten using existing error handling in SUEWS - can delete mod_error module from SUEWS_ESTM_functions
  !                - In SUEWS_ESTM_v2016, the first row is set to -999. This may be acceptable at the
  !                  start of the run but should be handled properly between blocks of met data? - need to check what's actually happening here.
  !                - Many duplicate functions in SUEWS_ESTM_functions need changing to the existing SUEWS versions.  
  !                - Are the following correctly initialised? T0_ibld, T0_ground, T0_wall, T0_roof, TN_wall, TN_roof, Tground, Twall, Troof, Tibld    
  !                - What are Nalb, sumalb, Nemis and Sumemis for? Are they used correctly?  
  !
    
  !SUEWS_ESTM_v2016
  ! Last modified HCW 14 Jun 2016
  !               HCW 27 Jun 2016 Corrected iESTMcount bug - now increases for all grids together
  !               HCW 30 Jun 2016 Major changes to handle grids and met blocks  
    
  !Contains calculation for each time step
  !Calculate local scale heat storage from single building and surroundings observations
  !OFferle, May 2003
  !
  !MODIFICATION HISTORY
  !             15 DECEMBER 2003
  !             (1) CHANGED AIR EXCHANGE RATE TO ALSO BE DEPENDENT ON OUTSIDE AIR TEMPERATURE
  !             (2) ADDED SPECIFIC HEAT CALCULATION FOR AIR
  !             (3) RH ADDED AS VAR(14) IN INPUT FILE
  !  12 JANUARY 2004
  !             (1) ADDED ESTIMATED AVG NADIR LOOKING EXTERNAL SURFACE TEMPERATURE TO OUTPUT
  !                 WEIGHTED BY SURFACE FRACTION AND SVF. SVF_ROOF=1, SVF_WALLS = 1-SVF_ground
  !             (2) ADDED NET RADIATION (RN) CALCULATIONS FOR ALL SURFACES AND AVG RN TO OUTPUT BASED ON
  !                 RADIOMETER VIEW FROM ZREF
  !       14 JANUARY 2004
  !             (1) MOVED GRID SPECIFIC PARAMETERS OUT OF NAMELIST, INTO PARAMETER FILE E.G. ALB,EM,F
  !                 T0 MAKE CHANGES IN LOCATIONS EASIER.
  !
  !       23 JANUARY 2004 - Version 2
  !             (1) PUT ALL TEMPERATURES INTO K
  !             (2) added interpolation routine to run on shorter timesteps.
  !                 tested so that it doesn't change solution for forced temperatures.
  !                 however in version 2 the energy balance at the surface isn't correctly solved
  !
  !       25 JANUARY 2004 - Version 3
  !             (1) added solution to energy balance at surface
  !             (2) removed some extraneous code
  !             (3) added vegetation fractions for future development
  !             (4) some changes to input namelist.
  !             (5) need to add wind speed dependence for exchange coefficients
  !             (6) changed the way Rn_net is calculated. also radiometer view factor relationships
  !             (7) added a calculation for heat loss/gain to outside air (that going into building mass is storage)
  !                 this is labelled QFBLD which it is in a sense.
  !       6 FEBRUARY 2004
  !             (1) ADDED INTERNAL VIEW FACTOR FILE FOR INTERNAL GEOMETRY, INCLUDING FIXED FLOOR TEMPERATURE
  !             (2) added MeanU, MinWS to config, and U to ILOC.
  !
  !      11 FEBRUARY 2004
  !             (1) added site lat, long, elevation to inputs
  !             (2) need zenith angle for wall direct radiation interception
  !
  !      15 JUNE 2004
  !             (1) CORRECTED OUTPUT OF T0 FOR FORCED SURFACE TEMPERATURES
  !             (2) MADE SOME CHANGES TO RADIATION, COMPUTATION OF AVERAGE ALBEDO, ADDED AVERAGE EMISSIVITY
  !             (3) ADDED OUTPUT FILE FOR RADIATION COMPONENTS
  !             (4) CHANGED INPUT CONFIG SO CONVECTIVE EXCHANGE COEFFS ARE NOT USER SELECTABLE
  !             (5) MAY STILL BE PROBLEMS WITH WALL RADIATIVE EXCHANGE AND AVERAGE ALBEDO
  !             (6) CHANGED THE WAY HEAT STORAGE IS COMPUTED IN HEATCONDUCTION MODULE BUT THIS SHOULD NOT CHANGE RESULTS
  !
  !      16 NOVEMBER 2004
  !             (1) CHANGED AIR EXCHANGE RATE TO BE BASES ON DAILY TEMPERATURE CHANGES.
  !             (2) WRITES OUT FINAL LAYER TEMPERATURES AND INCLUDES OPTION TO READ AT BEGINNING.
  !
  !DESCRIPTION: uses explicit time differencing to compute heat conduction through roof, walls,
  !             internal mass, and grounds (elements). Air heat storage is computed from average air temperature
  !             Boundary conditions are determined by measured surface temperature(s) or computed
  !             from energy balance at the surface.
  !             Internal air temperature can either be fixed or allowed to evolve in response
  !             to air mass exchanges and convective heating from internal surfaces.
  !
  !INPUT:
  !FORCING DATA: DTIME,KDOWN,LDOWN,TSURF,TAIR_OUT,TAIR_IN,TROOF,TWALL_AVG,TWALL_N,TWALL_E,TWALL_S,TWALL_W,Tground,RH,U
  !NAMELIST    : HEATSTORAGE_vFO.NML
  !    &config
  !    ifile=FORCING DATA
  !    ofile=OUTPUT FILE
  !    pfile=HEAT STORAGE PARAMETER FILE
  !    Nibld = INTERNAL MASS LAYERS
  !    Nwall = EXTERNAL WALL LAYERS
  !    Nroof = ROOF LAYERS
  !    Nground = ground/SOIL LAYERS
  !    LBC_soil = LOWER BOUNDARY CONDITION FOR ground/SOIL
  !    iloc= INPUT COLUMNS IN DATA FILE
  !    evolveTibld= USE DIAGNOSTIC VALUE FOR INTERNAL BUILDING TEMPERATURE
  !                        0: don't use, use measured
  !                        1: TURN ON USE when temp goess ABOVE TINT_ON, off when temp is below TINT_OFF
  !                        2: always use diagnostic
  !    THEAT_ON= TEMPERATURE AT WHICH HEAT CONTROL IS TURNED ON
  !    THEAT_OFF= TEMPERATURE AT WHICH HEAT CONTROL IS TURNED OFF
  !       THEAT_FIX = Fixed internal temperature for climate control
  !    oneTsurf= USE SINGLE SURFACE TEMPERATURE TO DRIVE ALL LAYERS
  !    radforce= USE RADIATIVE ENERGY BALANCE TO DRIVE EXTERNAL TEMPERATURES
  !    maxtimestep=302, maximum time step in s for filling data gaps.
  !       Alt = STATION HEIGHT (m) FOR PRESSURE CALCULATION
  !       SPINUP = NUMBER OF LINES TO USE FOR SPINUP (REPEATS THESE LINES BUT ONLY OUTPUTS THE 2ND TIME)
  !    INITTEMP = if TRUE INITIALIZES TEMPERATURES TO THOSE IN FINALTEMP.TXT FILE
  !    CH_ibld = INTERNAL BUILDING CONVECTIVE EXCHANGE COEFFICIENT
  !       **** THESE SHOULD DEPEND ON WIND SPEED BUT CURRENTLY DO NOT ****
  !       CHAIR = CONVECTIVE EXCHANGE COEFFICIENT FOR ROOF
  !       chair_ground = ... FOR ground
  !       chair_wall = ... FOR WALL
  !    /
  ! ***************** PARAMETER FILE VARIABLES
  !               fveg = FRACTION OF ground SURFACE COVERED WITH VEG
  !               zveg = VEGETATION HEIGHT
  !               alb_veg = VEGETATION ALBEDO
  !               em_veg = VEGETATION EMISSIVITY
  !               ZREF = REFERENCE HEIGHT FOR FLUX CALCULATION
  !               BldgH    = mean building height
  !               HW    = CANYON ASPECT RATION
  !               f_X   = FRACTION OF X WHERE X IS INTERNAL, WALL, ROOF, ground
  !               Alb_x = ALBEDO OF X
  !               em_ibld = EMISSIVITY OF X
  !               TX    = INITIAL LAYER TEMPERATURES
  !               zX    = LAYER THICKNESS
  !               kX    = LAYER THERMAL CONDUCTIVITY
  !               ribld = LAYER VOLUMETRIC HEAT CAPACITY
  !
  !****************** INTERNAL VIEW FACTOR FILE
  !OUTPUT:      fixed format text with single header, heatstorage for all elements, and temperatures
  !             for each element-layer.
  !===============================================================================
  
  USE initial
  USE meteo                                                               !!FO!! :METEOMOD.f95
  !USE mod_error
  USE mod_interp                                                          !!FO!! :mod_interp.f95
  USE mod_solver                                                          !!FO!! :mod_solver.f95
  USE modSolarCalc                                                        !!FO!! :modsolarcalc.f95
  USE MathConstants                                                       !!FO!! :MathConstants_module.f95
  USE PhysConstants
  USE heatflux
  USE ESTM_data
  USE mod_z
  USE data_in
  USE sues_data
  USE gis_data
  USE allocateArray
  USE time
  USE defaultNotUsed

  IMPLICIT NONE


  !Output to SUEWS
  REAL(KIND(1d0)),INTENT(out)::QSnet
  !Input from SUEWS, corrected as Gridiv by TS 09 Jun 2016
  INTEGER,INTENT(in)::Gridiv, ir
  
  !Use only in this subroutine
  INTEGER::i, ii
  INTEGER:: Tair2Set=0
  REAL(KIND(1d0))::AIREXHR, AIREXDT
  REAL(KIND(1d0)),DIMENSION(2)::bc
  REAL(KIND(1d0))::chair_ground,chair_wall
  REAL(KIND(1d0))::EM_EQUIV
  REAL(KIND(1d0))::kdz
  REAL(KIND(1d0))::kup_estm,LUP_net,kdn_estm
  REAL(KIND(1d0))::QHestm
  REAL(KIND(1d0))::QFBld !Anthropogenic heat from HVAC
  REAL(KIND(1d0))::shc_airbld
  REAL(KIND(1d0))::sw_hor,sw_vert
  REAL(KIND(1d0))::T0
  REAL(KIND(1d0))::Tinternal,Tsurf_all,Troof_in,Troad,Twall_all,Tw_n,Tw_e,Tw_s,Tw_w
  REAL(KIND(1d0))::Twallout(5),Troofout(5),Tibldout(5),Tgroundout(5)
  REAL(KIND(1d0))::Tadd,Tveg
  REAL(KIND(1d0))::Tairmix
  REAL(KIND(1d0))::RN
  REAL(KIND(1d0))::Rs_roof,Rl_roof,RN_ROOF
  REAL(KIND(1d0))::Rs_wall,Rl_wall,RN_WALL
  REAL(KIND(1d0))::Rs_ground,Rl_ground,RN_ground
  REAL(KIND(1d0))::Rs_ibld,Rl_ibld
  REAL(KIND(1d0))::Rs_iroof,Rl_iroof
  REAL(KIND(1d0))::Rs_iwall,Rl_iwall
  REAL(KIND(1d0))::zenith_rad
  REAL(KIND(1d0))::dum(50)
  REAL(KIND(1d0)),PARAMETER::WSmin=0.1  ! Check why there is this condition. S.O.
  LOGICAL::radforce, groundradforce
  
  radforce       = .FALSE.
  groundradforce = .FALSE. !Close the radiation scheme in original ESTM S.O.O.

  ! Set -999s for first row
  dum=(/-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
       -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
       -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
       -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
       -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999./)

  !External bulk exchange coefficients - set these somewhere more sensible***
  CHR=0.005
  CHAIR=CHR
  CHAIR_ground=CHAIR
  CHAIR_WALL=CHAIR
        
  !Get met data for use in ESTM subroutine     
  kdn_estm=avkdn
  WS=avu1
  IF (WS<WSMin) WS=WSmin
  Tair1=Temp_C+C2K
  ! Set initial value of Tair2 to air temp
  IF(Gridiv == 1) Tair2Set = Tair2Set+1
  IF(Tair2Set==1) THEN
    Tair2=Temp_C+C2K 
  ELSE
    Tair2 = Tair2_grids(Gridiv)
    ! Also get other variables for this grid
    Tievolve = Tievolve_grids(Gridiv)   
    lup_ground = lup_ground_grids(Gridiv)
    lup_wall = lup_wall_grids(Gridiv)
    lup_roof = lup_roof_grids(Gridiv)  
    T0_ibld = T0_ibld_grids(Gridiv)
    T0_ground = T0_ground_grids(Gridiv)
    T0_wall = T0_wall_grids(Gridiv)
    T0_roof = T0_roof_grids(Gridiv)
    TN_wall = TN_wall_grids(Gridiv)
    TN_roof = TN_roof_grids(Gridiv)
    Tground(:) = Tground_grids(:,Gridiv)
    Twall(:) = Twall_grids(:,Gridiv)
    Troof(:) = Troof_grids(:,Gridiv)
    Tibld(:) = Tibld_grids(:,Gridiv)
    Tw_4 = Tw_4_grids(:,:,Gridiv)  
    
  ENDIF
  
  ! Get Ts from Ts5min data array   
  Tinternal  = Ts5mindata(ir,cTs_Tiair)
  Tsurf_all  = Ts5mindata(ir,cTs_Tsurf)
  Troof_in   = Ts5mindata(ir,cTs_Troof)
  Troad      = Ts5mindata(ir,cTs_Troad)
  Twall_all  = Ts5mindata(ir,cTs_Twall)

  Tw_n       = Ts5mindata(ir,cTs_Twall_n)
  Tw_e       = Ts5mindata(ir,cTs_Twall_e)
  Tw_s       = Ts5mindata(ir,cTs_Twall_s)
  Tw_w       = Ts5mindata(ir,cTs_Twall_w)
 
  !    if (any(isnan(Ts5mindata))) then                                    !!FO!! can't use data when time gap is too big (or neg.) or data is NaN
  !        if (spindone) then                                                  !!FO!! writes a line of NaNs
  !            write(20,'(1F8.4,I6,100f10.1)') dectime,it,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
  !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
  !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
  !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.
  !        endif
  !        return ! changed from cycle            !!FO!! returns to beginning of do-loop, i.e. doesn't make any heat calculations
  !    endif
  !!FO!! this loop is used to run through the calculations (SPINUP number of times) to achieve better numerical stability
  !!FO!! when if condition is met the program starts over from the beginning of the input file and the calculations are performed
  !!FO!! ...once again, but this time saved to the output file
  !!FO!! it's because of this if statement arrangement that the output result starts at input time #2 (if not SPINUP=0 as originally in heatstorage_vFO.nml)
  !    IF (NLINESREAD>datalines.AND..NOT.SPINDONE) THEN                        !!FO!! at this stage spindone = .false.
  !        NLINESREAD=0
  !        SPINDONE=.TRUE.
  !PRINT*, "SPUNUP"
  !        return ! changed from cycle
  !    ENDIF

  !! Write first row of each met block as -999  
  !IF (first) THEN  !Set to true in ESTM_initials
  !!   !Tair2=Temp_C+C2K !This is now set in SUEWS_translate for ir=0 only
  !!   ! first=.FALSE.
  !   IF(Gridiv == NumberOfGrids) first=.FALSE.  !Set to false only after all grids have run
  !   dataOutESTM(ir,1:32,Gridiv)=(/REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),&
  !        REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,(dum(ii),ii=1,27)/)
  !   RETURN
  !ENDIF
  
  ! What are these constants? - Need defining somewhere
  zenith_rad=zenith_deg/180*PI
  IF (zenith_rad>0.AND.zenith_rad<PI/2.-HW) THEN  !ZENITH MUST BE HIGHER THAN BUILDINGS FOR DIRECT INTERCEPTION
     tanzenith = MIN(TAN(zenith_rad),5.67) !LIMITS TO ANGLES LESS THAN 80 EVEN FOR LOW HW
     tanzenith = tanzenith*kdn_estm/(1370*COS(zenith_rad)) !REDUCTION FACTOR FOR MAXIMUM
  ELSE
     tanzenith = 0.
  ENDIF

  SHC_air=HEATCAPACITY_AIR(Tair1,avrh,Press_hPa)   ! Use SUEWS version
  Tair24HR=EOSHIFT(Tair24hr, 1, Tair1, 1) !!!***Check this
  Tairday=SUM(Tair24HR)/(24*nsh)

    
  !Evolution of building temperature from heat added by convection
  SELECT CASE(evolvetibld)   !EvolveTiBld specifies which internal building temperature approach to use
  CASE(0); diagnoseTi=.FALSE.; HVAC=.FALSE. !use data in file                                !!FO!! use measured indoor temperature (Tref in Lodz2002HS.txt)
  CASE(1);                                                                                   !!FO!! use of HVAC to counteract T changes
     diagnoseTi=.TRUE.
     IF (Tievolve>THEAT_OFF) THEN   !THEAT_OFF now converted to Kelvin in ESTM_initials - HCW 15 Jun 2016
     !IF (Tievolve>THEAT_OFF+C2K) THEN
        HVAC=.FALSE.
     ELSEIF (Tievolve<THEAT_ON) THEN   !THEAT_OFF now converted to Kelvin in ESTM_initials - HCW 15 Jun 2016
     !ELSEIF (Tievolve<THEAT_ON+C2K) THEN
        HVAC=.TRUE.
     ENDIF
  CASE(2); diagnoseTi=.TRUE.                                                                 !!FO!! convection between ibld and inside of external walls(?)
  END SELECT

  !ASSUME AIR MIXES IN PROPORTION TO # OF EXCHANGES
  IF (Tairday>20.+C2K.AND.Tievolve>25.+C2K.AND.TAIR1<Tievolve.AND..NOT.HVAC) THEN
     AIREXHR = 2.0  !Windows or exterior doors on 3 sides (ASHRAE 1981 22.8)
  ELSEIF (Tairday<17.+C2K.OR.HVAC) THEN
     AIREXHR = 0.5 !No window or exterior doors, storm sash or weathertripped (ASHRAE 1981 22.8)
  ELSE
     AIREXHR = 1.0
  ENDIF

  AIREXDT=AIREXHR*(Tstep/3600.0)
  shc_airbld=HEATCAPACITY_AIR(TiEVOLVE,avrh,Press_hPa)
  IF (shc_airbld<minshc_airbld) minshc_airbld=shc_airbld

  !internal convective exchange coefficients                         !!FO!! ibldCHmod = 0 originally
  !iBldCHmod specifies method for convective exchange coeffs
  IF (ibldCHmod==1) THEN       !ASHRAE 2001
     CH_ibld  = 1.31*(ABS(T0_ibld-Tievolve))**0.25/shc_airbld
     CH_iwall = 1.31*(ABS(TN_wall-Tievolve))**0.25/shc_airbld
     CH_iroof = 1.52*(ABS(TN_roof-Tievolve))**0.25/shc_airbld
     IF (ABS(TN_roof-Tievolve)>0) CH_iroof=CH_iroof*0.39 !effect of convection is weaker downward
  ELSEIF (ibldCHmod==2) THEN   !Awbi, H.B. 1998, Energy and Buildings 28: 219-227
     CH_ibld  = 1.823*(ABS(T0_ibld-Tievolve))**0.293/shc_airbld
     CH_iwall = 1.823*(ABS(TN_wall-Tievolve))**0.293/shc_airbld
     CH_iroof = 2.175*(ABS(TN_roof-Tievolve))**0.308/shc_airbld
     IF (ABS(TN_roof-Tievolve)>0) CH_iroof=0.704*(ABS(TN_roof-Tievolve))**0.133/shc_airbld !effect of convection is weaker downward
  ENDIF

  !Evolving T = (Previous Temp + dT from Sensible heat flux) mixed with outside air
  !ASSUMES THE CH_BLD INCLUDES THE EFFECT OF VENTILATION RATE IN m/s (e.g. if a normal CH is .005 and
  !the value here is .003 the assumed ventilation is 0.6 m/s                                 !!FO!! CH_ibld=0.0015 from heatstorage_Barbican.nml => ventilation=0.3 m/s
  Tairmix =  (Tievolve + TAIR1*AIREXDT)/(1.0+AIREXDT)
  QFBld= froof*(Tievolve-Tairmix)*shc_airbld*BldgH/Tstep !heat added or lost, requires cooling or heating if HVAC on

  !!FO!! CH_xxxx has unit [m/s]  !!**HCW what is going on with tstep here??
  Tievolve = Tairmix+Tstep/BldgH/finternal* &                                                                         !!FO!! finternal(=froof+fibld+fwall) => normalisation of fractions
       (CH_ibld*fibld*(T0_ibld-Tievolve)+CH_iroof*froof*(TN_roof-Tievolve)+CH_iwall*fwall*(TN_wall-Tievolve))      !!FO!! [K] = [K] + [s/m]*([m/s]*([K]))

  IF (.NOT.diagnoseTi) Tievolve=Tinternal+C2K
  IF (HVAC) THEN !Run up/down to set point +/- 1 degree with adjustment of 90% per hour
     Tadd=(SIGN(-1.0d0,THEAT_fix-Tievolve)+THEAT_fix-Tievolve)*MIN(4.*Tstep/3600.0,0.9) !!**HCW check??
     Tievolve=Tievolve+Tadd
  ENDIF

  
  !========>RADIATION<================
  IF (kdn_estm<0) kdn_estm=0. !set non-zero shortwave to zero  !Should this be moved up to line 183/4?

  !external components, no diffuse
  !for reflections complete absorption is assumed
  !for shortwave these are net values
  !for longwave these are incoming only
  !MUST DIVIDE SHORTWAVE INTO DIRECT AND DIFFUSE
  sw_hor =kdn_estm           !incoming solar on horizontal surface
  sw_vert=kdn_estm*tanzenith !incoming solar on vertical surface = kdown(obs)*sin(zenith)/cos(zenith)

  Rs_roof=svf_roof*(1.0-alb_roof)*sw_hor
  Rl_roof=svf_roof*em_roof*ldown

  Rs_ground=svf_ground*(1.-alb_ground)*sw_hor+&
       zvf_ground*svf_wall*alb_wall*sw_vert*(1-alb_ground)+&
       zvf_ground*svf_ground*alb_ground*sw_hor*xvf_wall*alb_wall

  Rl_ground=svf_ground*ldown*em_ground+zvf_ground*(lup_wall+svf_wall*ldown*(1-em_wall))*em_ground

  Rs_wall=svf_wall*(1.-alb_wall)*sw_vert+&
       zvf_wall*svf_wall*alb_wall*sw_vert*(1.+zvf_wall*alb_wall)+&
       xvf_wall*svf_ground*alb_ground*sw_hor*(1-alb_wall)+&
       zvf_ground*xvf_wall*svf_ground*alb_ground*sw_hor*alb_wall

  !wall to wall exchange handled simultaneously with seb calc
  Rl_wall=svf_wall*ldown*em_wall+zvf_wall*svf_wall*ldown*(1-em_wall)*em_wall+&
       xvf_wall*(lup_ground+svf_ground*ldown*(1-em_ground))*em_wall

  !DIFFICULT TO DETERMINE WHAT THIS IS EXACTLY, DONT INCLUDE WALLS
  kup_estm=kdn_estm-RVF_ROOF*Rs_roof-(RVF_ground+RVF_WALL)*Rs_ground/svf_ground-RVF_VEG*ALB_VEG*kdn_estm
  IF (kdn_estm > 10 .AND. kup_estm > 0) THEN
     alb_avg = kup_estm/kdn_estm
     sumalb  = sumalb+alb_avg
     Nalb    = Nalb+1
  ENDIF

  
  !internal components
  Rs_ibld=0 ! This could change if there are windows (need solar angles or wall svf * fraction glazing * transmissivity)
  !internal incoming longwave terms do not include the view factors for its own surface e.g. for ibld and walls
  !added floor view factors
  Rl_ibld=SBConst*(ivf_iw*em_w*TN_wall**4 +&
       ivf_ir*em_r*TN_roof**4 +&
       ivf_if*em_f*Tfloor**4)
  Rs_iwall=0
  Rl_iwall=SBConst*(ivf_wi*em_i*T0_ibld**4 +&
       ivf_wr*em_r*TN_roof**4 +&
       ivf_wf*em_f*Tfloor**4)
  Rs_iroof=0
  Rl_iroof=SBConst*(ivf_ri*em_i*T0_ibld**4 +&
       ivf_rw*em_w*TN_wall**4 +&
       ivf_rf*em_f*Tfloor**4)

  !========>INTERNAL<================
  bctype=.FALSE.
  kdz=2*kibld(1)/zibld(1)
  Pcoeff=(/em_ibld*SBConst*(1-ivf_ii*em_ibld),0.0d0,0.0d0,kdz+shc_airbld*CH_ibld,&
       -kdz*Tibld(1)-shc_airbld*CH_ibld*Tievolve-Rs_ibld-Rl_ibld/)
  T0_ibld=NewtonPolynomial(T0_ibld,Pcoeff,conv,maxiter)
  bc(1)=T0_ibld                                                       !!FO!! this leads to Tibld(1) = Tibld(3) , i.e. ...
  bc(2)=bc(1)                                                         !!FO!! temperature equal on both sides of inside wall
  CALL heatcond1d(Tibld,Qsibld,zibld(1:Nibld),REAL(Tstep,KIND(1d0)),kibld(1:Nibld),ribld(1:Nibld),bc,bctype)

  !========>WALLS<================
  bctype=.FALSE.
  kdz=2*kwall(nwall)/zwall(nwall)
  Pcoeff=(/em_ibld*SBConst*(1-ivf_ww*em_ibld),0.0d0,0.0d0,kdz+shc_airbld*CH_iwall,&
       -kdz*Twall(nwall)-shc_airbld*CH_iwall*Tievolve-Rs_iwall-Rl_iwall/)
  TN_wall=NewtonPolynomial(TN_wall,Pcoeff,conv,maxiter)
  bc(2)=TN_wall                                                       !!FO!! boundary condition #2 = inner surface Twall, originally from lodz_parms_ltm.txt or finaltemp.txt

  IF (TsurfChoice<2 .OR.radforce) THEN
     IF (radforce) THEN                                              !!FO!! 1st prio: radforce
        kdz=2*kwall(1)/zwall(1)
        Pcoeff=(/em_wall*SBConst*(1-zvf_wall*em_wall),0.0d0,0.0d0,kdz+shc_air*chair_wall*WS,&
             -kdz*Twall(1)-shc_air*chair_wall*WS*Tair1-Rs_wall-Rl_wall/)
        T0_wall=NewtonPolynomial(T0_wall,Pcoeff,conv,maxiter)
        bc(1)=T0_wall                                               !!FO!! boundary condition #1 = outer surface Twall, originally from lodz_parms_ltm.txt or finaltemp.txt
     ELSEIF (TsurfChoice==0) THEN
        bc(1)=Tsurf_all+C2K; T0_wall=bc(1)
     ELSEIF (TsurfChoice==1) THEN
        bc(1)=Twall_all+C2K; T0_wall=bc(1)
     ENDIF                                                           !!FO!! Tsoil in Lodz2002HS.txt NB => Lodz2002HS.txt doesn't work with onewall = TRUE

     CALL heatcond1d(Twall,Qswall,zwall(1:nwall),REAL(Tstep,KIND(1d0)),kwall(1:nwall),rwall(1:nwall),bc,bctype)     !!FO!! new set of Twalls are calculated from heat conduction through wall

  ELSEIF(TsurfChoice==2) THEN!SPECIAL FOR 4 WALLS
     T0_wall=0.
     DO i=1,4 !do 4 walls
        bc(1)=Tw_n+Tw_e+Tw_s+Tw_w+C2K; T0_wall=T0_wall+bc(1)
        CALL heatcond1d(Tw_4(:,i),Qs_4(i),zwall(1:nwall),REAL(Tstep,KIND(1d0)),kwall(1:nwall),rwall(1:nwall),bc,bctype)
     ENDDO
     !Take average of 4 wall values
     T0_wall=T0_wall/4.
     Qswall = SUM(Qs_4)/4.
     Twall = SUM(Tw_4,2)/4.
  ENDIF

  !========>ROOF<================
  bctype=.FALSE.
  kdz=2*kroof(nroof)/zroof(nroof)
  Pcoeff=(/em_ibld*SBConst,0.0d0,0.0d0,kdz+shc_airbld*CH_iroof,&
       -kdz*Troof(nroof)-shc_airbld*CH_iroof*Tievolve-Rs_iroof-Rl_iroof/)
  TN_roof=NewtonPolynomial(TN_roof,Pcoeff,conv,maxiter)
  bc(2)=TN_roof

  IF (radforce) THEN
     kdz=2*kroof(1)/zroof(1)
     Pcoeff=(/em_roof*SBConst,0.0d0,0.0d0,kdz+shc_air*chair*WS,&
          -kdz*Troof(1)-shc_air*chair*WS*Tair1-Rs_roof-Rl_roof/)
     T0_roof=NewtonPolynomial(T0_roof,Pcoeff,conv,maxiter)
     bc(1)=T0_roof
  ELSEIF (TsurfChoice==0) THEN
     bc(1)=Tsurf_all+C2K; T0_roof=bc(1)
  ELSE
     bc(1)=Troof_in+C2K; T0_roof=bc(1)
  ENDIF

  CALL heatcond1d(Troof,Qsroof,zroof(1:nroof),REAL(Tstep,KIND(1d0)),kroof(1:nroof),rroof(1:nroof),bc,bctype)
  

  !========>ground<================
  bctype=.FALSE.
  kdz=2*kground(1)/zground(1)

  IF (radforce.OR.groundradforce) THEN
     Pcoeff=(/em_ground*SBConst,0.0d0,0.0d0,kdz+shc_air*chair_ground*WS,&
          -kdz*Tground(1)-shc_air*chair_ground*WS*Tair1-Rs_ground-Rl_ground/)
     T0_ground=NewtonPolynomial(T0_ground,Pcoeff,conv,maxiter)
     bc(1)=T0_ground
  ELSEIF (TsurfChoice==0) THEN
     bc(1)=Tsurf_all+C2K; T0_ground=bc(1)
  ELSE
     bc(1)=Troad+C2K; T0_ground=bc(1)
  ENDIF

  bc(2)=LBC_soil+C2K
  !     bc(2)=0.; bctype(2)=.t.

  IF ( fground/=0. )   THEN   ! check fground==0 scenario to avoid division-by-zero error, TS 21 Jul 2016
     CALL heatcond1d(Tground,Qsground,zground(1:Nground),REAL(Tstep,KIND(1d0)),kground(1:Nground),rground(1:Nground),bc,bctype)
  ELSE
     Qsground=NAN
  END IF

  Qsair = fair*SHC_air*(Tair1-Tair2)/Tstep
  Qsibld = Qsibld*fibld
  Qswall = Qswall*fwall
  Qsroof = Qsroof*froof
  Qsground = Qsground*fground
  Qsnet = Qsibld + Qswall + Qsroof + Qsground                              !!FO!! QSair not included; called QSNET in output file (column #10)

  !write(*,*) Qsair, QSibld, Qswall, Qsroof, Qsground, Qsnet
  
  !========>Radiation<================
  !note that the LUP for individual components does not include reflected
  LUP_ground = SBConst*EM_ground*T0_ground**4
  LUP_WALL   = SBConst*EM_WALL*T0_WALL**4
  LUP_ROOF   = SBConst*EM_ROOF*T0_ROOF**4
  TVEG       = TAIR1
  LUP_VEG    = SBConst*EM_VEG*TVEG**4
  T0         = RVF_ground*T0_ground+RVF_WALL*T0_WALL+RVF_ROOF*T0_ROOF+RVF_VEG*TVEG
  LUP_net    = RVF_ground*LUP_ground+RVF_WALL*LUP_WALL+RVF_ROOF*LUP_ROOF+RVF_VEG*LUP_VEG
  EM_EQUIV   = LUP_net/(SBConst*T0**4) !!FO!! apparent emissivity of the atmosphere [cloudless sky: >� Ldown from gases in the lowest 100 m] calculated from surface at T0
  RN_ground  = rs_ground+rl_ground-lup_ground
  RN_ROOF    = rs_roof+rl_roof-lup_roof
  RN_WALL    = rs_wall+rl_wall-lup_wall*(1-zvf_wall*em_wall)
  RN         = kdn_estm-kup_estm+ldown*EM_EQUIV-lup_net !!FO!! average net radiation (at z > zref ????) = shortwave down - shortwave up + [longwave down * apparent emissivity] - longwave up
  QHestm     = (T0-Tair1)*CHair*SHC_air*WS
  sumemis    = sumemis+EM_EQUIV
  nemis      = nemis+1

  ! IF (SPINDONE) THEN                                                      !!FO!! only the last set of values in the time interpolation loop is written to file

  IF (Nwall<5)THEN
     Twallout  =(/Twall,(dum(ii),  ii=1,(5-Nwall))/)
  ELSE
     Twallout=Twall
  ENDIF

  IF (Nroof<5) THEN
     Troofout  =(/Troof,(dum(ii),  ii=1,(5-Nroof))/);
  ELSE
     Troofout=Troof
  ENDIF

  IF (Nground<5)THEN
     Tgroundout=(/Tground,(dum(ii),ii=1,(5-Nground))/)
  ELSE
     Tgroundout=Tground
  ENDIF

  IF (Nibld<5)THEN
     Tibldout  =(/Tibld,(dum(ii),  ii=1,(5-Nibld))/)
  ELSE
     Tibldout=Tibld
  ENDIF
  
  dataOutESTM(ir,1:32,Gridiv)=(/REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),&
       REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,Qsnet,Qsair,Qswall,Qsroof,Qsground,Qsibld,&!11
       Twallout,Troofout,Tgroundout,Tibldout,Tievolve/)!21
  !kdn_estm,kup_estm,ldown,lup_net,RN,& !10
  !   Qsnet,Qsair,QHestm,QFBld,T0,Qswall,Qsroof,Qsground,Qsibld,RN_WALL,RN_ROOF,RN_ground,&   !12
  !  Twallout,TN_Wall,Troofout,TN_roof,Tgroundout,Tibldout,Tievolve,zenith_deg/)!8+XX
  ! WRITE(30,'(1F8.4,10F10.1)') dectime, kdn_estm, rs_wall, rs_roof, rs_ground, ldown,rl_wall, rl_roof, rl_ground
  !endif

  Tair2=Tair1

  ! Save variables for this grid
  Tair2_grids(Gridiv)=Tair1
  lup_ground_grids(Gridiv) = lup_ground
  lup_wall_grids(Gridiv) = lup_wall
  lup_roof_grids(Gridiv) = lup_roof  
  Tievolve_grids(Gridiv) = Tievolve
  T0_ibld_grids(Gridiv) = T0_ibld
  T0_ground_grids(Gridiv) = T0_ground
  T0_wall_grids(Gridiv) = T0_wall
  T0_roof_grids(Gridiv) = T0_roof
  TN_wall_grids(Gridiv) = TN_wall
  TN_roof_grids(Gridiv) = TN_roof
  Tground_grids(:,Gridiv) = Tground(:)
  Twall_grids(:,Gridiv) = Twall(:)
  Troof_grids(:,Gridiv) = Troof(:)
  Tibld_grids(:,Gridiv) = Tibld(:)
  Tw_4_grids(:,:,Gridiv) = Tw_4(:,:)
  
END SUBROUTINE ESTM_v2016
