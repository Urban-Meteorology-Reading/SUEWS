!This subroutine does the actual calculations of the SUEWS code (mostly old SUEWS_Temporal).
!Made by LJ and HW Oct 2014
!Gives in the grid ID (Gridiv) and number of line in the met forcing data to be analyzed (ir)
!Last modification
!
!Last modification:
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
!                 - Replaced soilmoist_state variable with soilstate (as seems to be duplicate)
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


  IMPLICIT NONE

  INTEGER        :: Gridiv,ir,i,ih,iMB
  LOGICAL        :: debug=.FALSE.
  REAL(KIND(1d0)):: idectime
  !real(kind(1d0)):: SnowDepletionCurve  !for SUEWS_Snow - not needed here (HCW 24 May 2016)
  REAL(KIND(1d0)):: lai_wt,qsatf
  INTEGER        :: irMax

  !==================================================================
  !==================================================================


  !Translate all data to the variables used in the model calculations
  IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_Translate...'
  CALL SUEWS_Translate(Gridiv,ir,iMB)
  ! ! load the final water states of the previous day to keep water balance, by TS 13 Apr 2016
  ! IF ( ir==1 ) THEN
  !      PRINT*, '********************************'
  !      PRINT*, 'starting state of', id,it,imin
  !    state(1:nsurf)     = stateDay(id-1,Gridiv,1:nsurf)
  !    soilmoist(1:nsurf) = soilmoistDay(id-1,Gridiv,1:nsurf)
  !      PRINT*, 'state:', state
  !      PRINT*, 'soilmoist', soilmoist
  ! END IF
  IF(Diagnose==1) WRITE(*,*) 'Calling RoughnessParameters...'
  CALL RoughnessParameters(Gridiv) ! Added by HCW 11 Nov 2014


  !=============Get data ready for the qs calculation====================
  IF(NetRadiationMethod==0) THEN !Radiative components are provided as forcing
     !avkdn=NAN                  !Needed for resistances for SUEWS.
     ldown     = NAN
     lup       = NAN
     kup       = NAN
     tsurf     = NAN
     lup_ind   = NAN
     kup_ind   = NAN
     tsurf_ind = NAN
     qn1_ind   = NAN
     Fcld      = NAN
  ENDIF

  IF(ldown_option==1) THEN
     Fcld = NAN
  ENDIF
  !=====================================================================
  ! Initialisation for OAF's water bucket scheme
  ! LUMPS only (Loridan et al. (2012))
  RAINRES = 0.
  RAINBUCKET = 0.
  E_MOD=0 !RAIN24HR = 0.;

  !=====================================================================

  !Initialize variables calculated at each 5-min timestep
  runoffAGveg             = 0
  runoffAGimpervious      = 0
  runoffWaterBody         = 0
  runoffSoil_per_tstep    = 0
  runoffSoil_per_interval = 0
  chSnow_per_interval     = 0
  mwh                     = 0 !Initialize snow melt and heat related to snowmelt
  fwh                     = 0
  Qm                      = 0 !Heat related to melting/freezing
  QmFreez                 = 0
  QmRain                  = 0
  Mw_ind                  = 0
  SnowDepth               = 0
  zf                      = 0
  deltaQi                 = 0
  swe                     = 0
  MwStore                 = 0
  WaterHoldCapFrac        = 0

  qh = -999 ! Added HCW 26 Feb 2015
  H  = -999 ! Added HCW 26 Feb 2015

  ! Calculate sun position
  idectime=dectime-halftimestep! sun position at middle of timestep before
  IF(Diagnose==1) WRITE(*,*) 'Calling sun_position...'
  CALL sun_position(year,idectime,timezone,lat,lng,alt,azimuth,zenith_deg)
  !write(*,*) DateTime, timezone,lat,lng,alt,azimuth,zenith_deg


  IF(CBLuse>=1)THEN ! If CBL is used, calculated Temp_C and RH are replaced with the obs.
     IF(Diagnose==1) WRITE(*,*) 'Calling CBL...'
     CALL CBL(ir,iMB,Gridiv)   !ir=1 indicates first row of each met data block
  ENDIF

  !Call the dailystate routine to get surface characteristics ready
  IF(Diagnose==1) WRITE(*,*) 'Calling DailyState...'
  CALL DailyState(Gridiv)

  IF(LAICalcYes==0)THEN
     lai(id-1,:)=lai_obs ! check -- this is going to be a problem as it is not for each vegetation class
  ENDIF

  !Calculation of density and other water related parameters
  IF(Diagnose==1) WRITE(*,*) 'Calling atmos_moist_lumps...'
  CALL atmos_moist_lumps(avdens)


  !======== Calculate soil moisture =========
  SoilMoistCap=0   !Maximum capacity of soil store [mm] for whole surface
  soilstate=0      !Area-averaged soil moisture [mm] for whole surface
  DO is=1,nsurf-1   !No water body included
     soilmoistCap=soilMoistCap+(soilstoreCap(is)*sfr(is)/NonWaterFraction)
     soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
  ENDDO

  !If loop removed HCW 26 Feb 2015
  !if (ir==1) then  !Calculate initial smd
  smd=soilmoistCap-soilstate
  !endif

  ! Calculate soil moisture for vegetated surfaces only (for use in surface conductance)
  vsmd=0
  DO is=ConifSurf,GrassSurf  !Vegetated surfaces only
     IF ( sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) ==0 ) THEN
        vsmd=0
     ELSE
        vsmd=vsmd+(soilstoreCap(is) - soilmoist(is))*sfr(is)/(sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))
     END IF
     !write(*,*) is, vsmd, smd
  ENDDO

  ! ===================NET ALLWAVE RADIATION================================
  IF(NetRadiationMethod>0)THEN

     IF (snowUse==0) snowFrac=snow_obs

     IF(ldown_option==1) THEN !Observed ldown provided as forcing
        ldown=ldown_obs
     ELSE
        ldown=-9              !to be filled in NARP
     ENDIF

     IF(ldown_option==2) THEN !observed cloud fraction provided as forcing
        fcld=fcld_obs
     ENDIF

     !write(*,*) DecidCap(id), id, it, imin, 'Calc - near start'

     ! Update variables that change daily and represent seasonal variability
     alb(DecidSurf)    = albDecTr(id) !Change deciduous albedo
     surf(6,DecidSurf) = DecidCap(id) !Change current storage capacity of deciduous trees
     ! Change EveTr and Grass albedo too
     alb(ConifSurf) = albEveTr(id)
     alb(GrassSurf) = albGrass(id)

     IF(Diagnose==1) WRITE(*,*) 'Calling NARP...'
     CALL NARP(SnowAlb,qn1_SF,qn1_S)
     !Temp_C,kclear,fcld,dectime,avkdn,avRH,qn1,kup,ldown,lup,tsurf,&
     !AlbedoChoice,ldown_option,Press_hPa,Ea_hPa,qn1_obs,&
     !zenith_deg,NetRadiationMethod,
  ELSE
     snowFrac = snow_obs
     qn1      = qn1_obs
     qn1_sf   = qn1_obs
     qn1_s    = qn1_obs
  ENDIF

  ! ===================SOLWEIG OUTPUT ========================================
  IF (SOLWEIGuse==1) THEN
     IF (OutInterval==imin) THEN
        IF (RunForGrid==-999) THEN
           IF(Diagnose==1) WRITE(*,*) 'Calling SOLWEIG_2014a_core...'
           CALL SOLWEIG_2014a_core(iMB)
           SolweigCount=SolweigCount+1
        ELSE
           IF (Gridiv == RunForGrid) THEN
              IF(Diagnose==1) WRITE(*,*) 'Calling SOLWEIG_2014a_core...'
              CALL SOLWEIG_2014a_core(iMB)
              SolweigCount=SolweigCount+1
           ENDIF
        ENDIF
     ENDIF
  ELSE
     SOLWEIGpoi_out=0
  ENDIF

  ! ===================ANTHROPOGENIC HEAT FLUX================================
  ih=it-DLS
  IF(ih<0) ih=23

  IF(AnthropHeatMethod==1) THEN
     IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_1...'
     CALL SAHP_1(qf_sahp,QF_SAHP_base,QF_SAHP_heat,id,ih,imin)
     qn1_bup=qn1
     qn1=qn1+QF_SAHP
  ELSEIF(AnthropHeatMethod==2) THEN
     IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_2...'
     CALL SAHP_2(qf_sahp,QF_SAHP_base,QF_SAHP_heat,id,ih,imin)
     qn1_bup=qn1
     qn1=qn1+QF_SAHP
  ELSEIF(AnthropHeatMethod==3) THEN
     IF(Diagnose==1) WRITE(*,*) 'Calling SAHP_3...'
     CALL SAHP_3(qf_sahp,id,ih,imin)
     qn1_bup=qn1
     qn1=qn1+QF_SAHP
  ELSE
     qn1_bup=qn1
     qn1=qn1+qf
  ENDIF

  ! -- qn1 is now QSTAR+QF (net all-wave radiation + anthropogenic heat flux)
  ! -- qn1_bup is QSTAR only
  IF(AnthropHeatMethod>=1) THEN
     qf=QF_SAHP
  ENDIF

  ! Calculate CO2 fluxes from anthropogenic components
  IF(Diagnose==1) WRITE(*,*) 'Calling CO2_anthro...'
  CALL CO2_anthro(id,ih,imin)

  ! =================STORAGE HEAT FLUX=======================================

  IF(StorageHeatMethod==1) THEN           !Use OHM to calculate QS
     IF(OHMIncQF == 1) THEN      !Calculate QS using QSTAR+QF
        IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
        CALL OHM(Gridiv)
     ELSEIF(OHMIncQF == 0) THEN  !Calculate QS using QSTAR
        qn1=qn1_bup
        IF(Diagnose==1) WRITE(*,*) 'Calling OHM...'
        CALL OHM(Gridiv)
     ENDIF
  ENDIF

  ! use AnOHM to calculate QS, TS 14 Mar 2016
  IF (StorageHeatMethod==3) THEN
     IF ( OHMIncQF == 1 ) THEN    !Calculate QS using QSTAR+QF
        IF(Diagnose==1) WRITE(*,*) 'Calling AnOHM...'
        CALL AnOHM(Gridiv)
     ELSEIF(OHMIncQF == 0) THEN   !Calculate QS using QSTAR
        qn1=qn1_bup
        IF(Diagnose==1) WRITE(*,*) 'Calling AnOHM...'
        CALL AnOHM(Gridiv)
     END IF
  END IF

  !Calculate QS using ESTM
  IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
     !CALL ESTM_v2016(QSestm,iMB)
     IF(Diagnose==1) WRITE(*,*) 'Calling ESTM...'
     CALL ESTM_v2016(QSestm,Gridiv,ir)  ! iMB corrected to Gridiv, TS 09 Jun 2016
     QS=QSestm   ! Use ESTM qs
  ENDIF

  ! don't use these composite QS options at the moment, TS 10 Jun 2016
  ! IF (StorageHeatMethod>=10)THEN ! Chose which QS will be used in SUEWS and output file
  !    !write(800,*)id,it,QS,QSanOHM,QSestm
  !    IF(StorageHeatMethod==14)THEN
  !       QS=QSestm
  !    ELSEIF(StorageHeatMethod==13)THEN
  !       QS=QSanOHM
  !    ENDIF
  ! ENDIF



  ! For the purpose of turbulent fluxes, remove QF from the net all-wave radiation
  qn1=qn1_bup  !Remove QF from QSTAR

  ! -- qn1 is now QSTAR only



  !==================Energy related to snow melting/freezing processes=======
  IF (snowUse==1)  THEN

     IF(Diagnose==1) WRITE(*,*) 'Calling MeltHeat'
     CALL MeltHeat

     ! If snow on ground, no irrigation, so veg_fr same in each case
     !New fraction of vegetation.
     !IF(veg_type==1)THEN         ! area vegetated
     veg_fr = sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+&
          sfr(GrassSurf)*(1-snowFrac(GrassSurf))+sfr(BSoilSurf)*(1-snowFrac(BSoilSurf))+&
          sfr(WaterSurf)*(1-snowFrac(WaterSurf))

     !ELSEIF(veg_type==2)THEN     ! area irrigated
     !   !!!veg_fr=sfr(GrassISurf)*(1-snowFrac(GrassUSurf))
     !END IF

  END IF

  !==========================Turbulent Fluxes================================

  tlv=lv_J_kg/tstep_real !Latent heat of vapourisation per timestep

  IF(Diagnose==1) WRITE(*,*) 'Calling LUMPS_QHQE...'
  CALL LUMPS_QHQE !Calculate QH and QE from LUMPS
  IF(debug)WRITE(*,*)press_Hpa,psyc_hPA,i

  IF(Diagnose==1) WRITE(*,*) 'Calling WaterUse...'
  CALL WaterUse !Gives the external and internal water uses per timestep

  IF(Precip>0) THEN   !Initiate rain data [mm]
     pin=Precip
  ELSE
     pin=0
  ENDIF


  ! Get first estimate of sensible heat flux. Modified by HCW 26 Feb 2015
  ! Calculate kinematic heat flux (w'T') from sensible heat flux [W m-2] from observed data (if available) or LUMPS
  IF(qh_obs/=NAN) THEN   !if(qh_obs/=NAN) qh=qh_obs   !Commented out by HCW 04 Mar 2015
     H=qh_obs/(avdens*avcp)  !Use observed value
  ELSE
     IF(h_mod/=NAN) THEN
        H = h_mod/(avdens*avcp)   !Use LUMPS value
     ELSE
        H=(qn1*0.2)/(avdens*avcp)   !If LUMPS has had a problem, we still need a value
        CALL ErrorHint(38,'LUMPS unable to calculate realistic value for H_mod.',h_mod, dectime, notUsedI)
     ENDIF
  ENDIF

  !------------------------------------------------------------------

  IF(Diagnose==1) WRITE(*,*) 'Calling STAB_lumps...'
  CALL STAB_lumps(H,StabilityMethod,ustar,L_mod) !u* and Obukhov length out

  IF(Diagnose==1) WRITE(*,*) 'Calling AerodynamicResistance...'
  CALL AerodynamicResistance(RA,AerodynamicResistanceMethod,StabilityMethod,RoughLenHeatMethod,&
       ZZD,z0m,k2,AVU1,L_mod,Ustar,VegFraction,psyh)      !RA out

  IF (snowUse==1) THEN
     IF(Diagnose==1) WRITE(*,*) 'Calling AerodynamicResistance...'
     CALL AerodynamicResistance(RAsnow,AerodynamicResistanceMethod,StabilityMethod,3,&
          ZZD,z0m,k2,AVU1,L_mod,Ustar,VegFraction,psyh)      !RA out
  ENDIF

  IF(Diagnose==1) WRITE(*,*) 'Calling SurfaceResistance...'
  CALL SurfaceResistance(id,it)   !qsc and surface resistance out
  IF(Diagnose==1) WRITE(*,*) 'Calling BoundaryLayerResistance...'
  CALL BoundaryLayerResistance


  ! Calculate CO2 fluxes from biogenic components
  IF(Diagnose==1) WRITE(*,*) 'Calling CO2_biogen...'
  CALL CO2_biogen
  ! Sum anthropogenic and biogenic CO2 flux components to find overall CO2 flux
  Fc = Fc_anthro + Fc_biogen

  sae   = s_hPa*(qn1_SF+qf-qs)    !s_haPa - slope of svp vs t curve. qn1 changed to qn1_SF, lj in May 2013
  vdrc  = vpd_hPa*avdens*avcp
  sp    = s_hPa/psyc_hPa
  numPM = sae+vdrc/ra

  !write(*,*) numPM, sae, vdrc/ra, s_hPA+psyc_hPa, NumPM/(s_hPA+psyc_hPa)

  !=====================================================================
  !========= Water balance calculations ================================
  ! Needs to run at small timesteps (i.e. minutes)
  ! Previously, v2014b switched to 60/NSH min intervals here
  ! Now whole model runs at a resolution of tstep

  ! Initialise water balance variables
  qe               = 0
  ev               = 0
  ev_snow          = 0
  SurplusEvap      = 0
  evap             = 0
  chang            = 0
  runoff           = 0
  runoffSoil       = 0
  surplusWaterBody = 0

  ! Added by HCW 13 Feb 2015
  qe_per_tstep         = 0     ![W m-2]
  ev_per_tstep         = 0
  drain_per_tstep      = 0
  surf_chang_per_tstep = 0
  tot_chang_per_tstep  = 0
  state_per_tstep      = 0
  NWstate_per_tstep    = 0
  runoff_per_tstep     = 0

  ! Retain previous surface state and soil moisture state
  stateOld     = state     !State of each surface [mm] for the previous timestep
  soilmoistOld = soilmoist !Soil moisture of each surface [mm] for the previous timestep

  !============= Grid-to-grid runoff =============
  ! Calculate additional water coming from other grids
  ! i.e. the variables addImpervious, addVeg, addWaterBody, addPipes
  !call RunoffFromGrid(GridFromFrac)  !!Need to code between-grid water transfer

  ! Sum water coming from other grids (these are expressed as depths over the whole surface)
  AdditionalWater = addPipes+addImpervious+addVeg+addWaterBody  ![mm]

  ! Initialise runoff in pipes
  runoffPipes         = addPipes !Water flowing in pipes from other grids. No need for scaling??
  !! CHECK p_i
  runoff_per_interval = addPipes !pipe plor added to total runoff.

  !================== Drainage ===================
  ! Calculate drainage for each soil subsurface (excluding water body)
  IF(Diagnose==1) WRITE(*,*) 'Calling Drainage...'
  DO is=1,nsurf-1
     CALL Drainage(surf(6,is),surf(2,is),surf(3,is),surf(4,is))
     !HCW added and changed to surf(6,is) here 20 Feb 2015
     drain_per_tstep=drain_per_tstep+(drain(is)*sfr(is)/NonWaterFraction)   !No water body included
  ENDDO

  drain(WaterSurf) = 0  ! Set drainage from water body to zero

  ! Distribute water within grid, according to WithinGridWaterDist matrix (Cols 1-7)
  IF(Diagnose==1) WRITE(*,*) 'Calling ReDistributeWater...'
  CALL ReDistributeWater   !Calculates AddWater(is)

  !======== Evaporation and surface state ========
  IF(Diagnose==1) WRITE(*,*) 'Calling evap_SUEWS and SoilStore...'
  DO is=1,nsurf   !For each surface in turn
     IF (snowCalcSwitch(is)==1) THEN
        IF (sfr(is)/=0) THEN
           IF(Diagnose==1) WRITE(*,*) 'Calling SnowCalc...'
           CALL snowCalc
        ELSE
           snowFrac(is) = 0
           SnowDens(is) = 0
           SnowPack(is) = 0
        ENDIF
     ELSE
        CALL Evap_SUEWS   !Calculates ev [mm]
        rss_nsurf(is) = rss !Store rss for each surface
        CALL soilstore    !Surface water balance and soil store updates (can modify ev, updates state)
        ! if ( id>160 ) then
        !   stop "stop after Evap_SUEWS"
        ! end if

        evap(is)     = ev !Store ev for each surface
        ! Sum evaporation from different surfaces to find total evaporation [mm]
        ev_per_tstep = ev_per_tstep+evap(is)*sfr(is)

        ! Sum change from different surfaces to find total change to surface state
        surf_chang_per_tstep = surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)
        ! Sum runoff from different surfaces to find total runoff
        runoff_per_tstep     = runoff_per_tstep+runoff(is)*sfr(is)
        ! Calculate total state (including water body)
        state_per_tstep      = state_per_tstep+(state(is)*sfr(is))
        ! Calculate total state (excluding water body)
        IF (is.NE.WaterSurf) NWstate_per_tstep=NWstate_per_tstep+(state(is)*sfr(is)/NonWaterFraction)

        ChangSnow(is)  = 0
        runoffSnow(is) = 0

     ENDIF
  ENDDO  !end loop over surfaces


  ! Convert evaporation to latent heat flux [W m-2]
  qe_per_tstep = ev_per_tstep*tlv
  qeOut        = qe_per_tstep

  !============ Sensible heat flux ===============
  ! Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
  qh=(qn1+qf+QmRain)-(qeOut+qs+Qm+QmFreez)     !qh=(qn1+qf+QmRain+QmFreez)-(qeOut+qs+Qm)

  ! Calculate QH using resistance method (for testing HCW 06 Jul 2016)
  IF(ra/=0) THEN
     qh_r = avdens*avcp*(tsurf-Temp_C)/ra
  ELSE
     qh_r=NAN
  ENDIF


  !============ surface-level diagonostics ===============
  ! wind speed:
  CALL diagSfc(0.,0.,ustar,avU10_ms,0)
  ! temperature:
  CALL diagSfc(tsurf,qh,ustar,t2_C,1)
  ! humidity:
  CALL diagSfc(qsatf(tsurf,Press_hPa)*1000,& ! Saturation specific humidity at surface in g/kg
       qeOut,ustar,q2_gkg,2)
  !============ surface-level diagonostics end ===============


  !write(*,*) Gridiv, qn1, qf, qh, qeOut, qs, qn1+qf-qs
  !if(ir > 155 .and. ir <165) pause
  !if((qn1+qf-qs) - (qeOut) < -1)  then
  !    write(*,*) '!!', (qn1+qf-qs),(qeOut)
  !    pause
  !endif
  !if(ir > 600) pause
  !if(Gridiv == 3) write(*,*) ''
  !if(Gridiv == 3 .and. ir == 100) pause
  !!if(Gridiv == 3 .and.ir == 200) pause
  !!if(Gridiv == 3 .and.ir == 300) pause
  !if(Gridiv == 3 .and.ir == 400) pause
  !if(Gridiv == 3 .and. ir == irMax ) pause

  ! Calculate volume of water that will move between grids
  ! Volume [m3] = Depth relative to whole area [mm] / 1000 [mm m-1] * SurfaceArea [m2]
  ! Need to use these volumes when converting back to addImpervious, AddVeg and AddWater
  runoffAGimpervious_m3 = runoffAGimpervious/1000 *SurfaceArea
  runoffAGveg_m3        = runoffAGveg/1000 *SurfaceArea
  runoffWaterBody_m3    = runoffWaterBody/1000 *SurfaceArea
  runoffPipes_m3        = runoffPipes/1000 *SurfaceArea


  !=== Horizontal movement between soil stores ===
  ! Now water is allowed to move horizontally between the soil stores
  IF(Diagnose==1) WRITE(*,*) 'Calling HorizontalSoilWater...'
  CALL HorizontalSoilWater

  !========== Calculate soil moisture ============
  soilstate=0       !Area-averaged soil moisture [mm] for whole surface
  DO is=1,nsurf-1   !No water body included
     soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
     IF (soilstate<0) THEN
        CALL ErrorHint(62,'SUEWS_Calculations: total soilstate < 0 (just added surface is) ',soilstate,NotUsed,is)
     ELSEIF (soilstate>SoilMoistCap) THEN
        CALL ErrorHint(62,'SUEWS_Calculations: total soilstate > capacity (just added surface is) ',soilstate,NotUsed,is)
        !SoilMoist_state=soilMoistCap !What is this LJ 10/2010 - SM exceeds capacity, but where does extra go?HCW 11/2014
     ENDIF
  ENDDO  !end loop over surfaces
  ! Calculate soil moisture deficit
  smd=soilMoistCap-soilstate   !One value for whole surface
  smd_nsurf=SoilstoreCap-soilmoist   !smd for each surface

  ! Soil stores can change after horizontal water movements
  ! Calculate total change in surface and soil state
  tot_chang_per_tstep = surf_chang_per_tstep   !Change in surface state
  DO is=1,(nsurf-1)   !No soil for water surface (so change in soil moisture is zero)
     tot_chang_per_tstep = tot_chang_per_tstep + ((SoilMoist(is)-SoilMoistOld(is))*sfr(is))   !Add change in soil state
  ENDDO

  !=====================================================================
  !====================== Prepare data for output ======================

  ! Check surface composition (HCW 21 Jul 2016)
  ! if totally water surface, set output to -999 for columns that do not exist
  IF(sfr(WaterSurf)==1) THEN
     NWstate_per_tstep = NAN  !no non-water surface
     smd = NAN   !no soil store beneath water surface
     soilstate = NAN
     runoffSoil_per_tstep = NAN
     drain_per_tstep = NAN   !no drainage from water surf
  ENDIF

  ! Remove non-existing surface type from surface and soil outputs   ! Added back in with NANs by HCW 24 Aug 2016
  DO is=1,nsurf
     IF (sfr(is)<0.00001) THEN
        stateOut(is)= NAN
        smd_nsurfOut(is)= NAN
        !runoffOut(is)= NAN
        !runoffSoilOut(is)= NAN
     ELSE
        stateOut(is)=state(is)
        smd_nsurfOut(is)=smd_nsurf(is)
        !runoffOut(is)=runoff(is)
        !runoffSoilOut(is)=runoffSoil(is)
     ENDIF
  ENDDO

  ! Remove negative state   !ErrorHint added, then commented out by HCW 16 Feb 2015 as should never occur
  !if(st_per_interval<0) then
  !   call ErrorHint(63,'SUEWS_Calculations: st_per_interval < 0',st_per_interval,NotUsed,NotUsedI)
  !   !st_per_interval=0
  !endif

  !Set limits on output data to avoid formatting issues ---------------------------------
  ! Set limits as +/-9999, depending on sign of value
  ! errorHints -> warnings commented out as writing these slow down the code considerably when there are many instances
  IF(ResistSurf > 9999) THEN
     !CALL errorHint(6,'rs set to 9999 s m-1 in output; calculated value > 9999 s m-1',ResistSurf,notUsed,notUsedI)
     ResistSurf=9999
  ENDIF

  IF(l_mod > 9999) THEN
     !CALL errorHint(6,'Lob set to 9999 m in output; calculated value > 9999 m',L_mod,notUsed,notUsedI)
     l_mod=9999
  ELSEIF(l_mod < -9999) THEN
     !CALL errorHint(6,'Lob set to -9999 m in output; calculated value < -9999 m',L_mod,notUsed,notUsedI)
     l_mod=-9999
  ENDIF


  ! Set NA values   !!Why only these variables??  !!ErrorHints here too - error hints can be very slow here
  IF(ABS(qh)>pNAN) qh=NAN
  IF(ABS(qh_r)>pNAN) qh_r=NAN
  IF(ABS(qeOut)>pNAN) qeOut=NAN
  IF(ABS(qs)>pNAN) qs=NAN
  IF(ABS(ch_per_interval)>pNAN) ch_per_interval=NAN
  IF(ABS(surf_chang_per_tstep)>pNAN) surf_chang_per_tstep=NAN
  IF(ABS(tot_chang_per_tstep)>pNAN) tot_chang_per_tstep=NAN
  IF(ABS(soilstate)>pNAN) soilstate=NAN
  IF(ABS(smd)>pNAN) smd=NAN

  ! If measured smd is used, set components to -999 and smd output to measured one
  IF (SMDMethod>0) THEN
     smd_nsurf=NAN
     smd_nsurfOut=NAN
     smd=xsmd
  ENDIF

  ! Calculate areally-weighted LAI
  IF(iy == (iy_prev_t+1) .AND. (id-1) == 0) THEN   !Check for start of next year and avoid using lai(id-1) as this is at the start of the year
     lai_wt=0
     DO is=1,nvegsurf
        lai_wt=lai_wt+lai(id_prev_t,is)*sfr(is+2)
     ENDDO
  ELSE
     lai_wt=0
     DO is=1,nvegsurf
        lai_wt=lai_wt+lai(id-1,is)*sfr(is+2)
     ENDDO
  ENDIF

  ! Calculate areally-weighted albedo
  bulkalbedo = 0
  DO is=1,nsurf
     bulkalbedo = bulkalbedo + alb(is)*sfr(is)
  ENDDO

  ! Save qh and qe for CBL in next iteration
  IF(Qh_choice==1) THEN   !use QH and QE from SUEWS
     qhforCBL(Gridiv) = qh
     qeforCBL(Gridiv) = qeOut
  ELSEIF(Qh_choice==2)THEN   !use QH and QE from LUMPS
     qhforCBL(Gridiv) = h_mod
     qeforCBL(Gridiv) = e_mod
  ELSEIF(qh_choice==3)THEN  !use QH and QE from OBS
     qhforCBL(Gridiv) = qh_obs
     qeforCBL(Gridiv) = qe_obs
     IF(qh_obs<-900.OR.qe_obs<-900)THEN  ! observed data has a problem
        CALL ErrorHint(22,'Unrealistic observed qh or qe_value.',qh_obs,qe_obs,qh_choice)
     ENDIF
  ENDIF



  !=====================================================================
  !====================== Write out files ==============================
  !Define the overall output matrix to be printed out step by step
  dataOut(ir,1:ncolumnsDataOut,Gridiv)=(/REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,&   !5
       avkdn,kup,ldown,lup,tsurf,&
       qn1,qf,qs,qh,qeOut,&
       h_mod,e_mod,qh_r,&
       precip,ext_wu,ev_per_tstep,runoff_per_tstep,tot_chang_per_tstep,&
       surf_chang_per_tstep,state_per_tstep,NWstate_per_tstep,drain_per_tstep,smd,&
       FlowChange/nsh_real,AdditionalWater,&
       runoffSoil_per_tstep,runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,&
       int_wu,wu_EveTr,wu_DecTr,wu_Grass,&
       (smd_nsurfOut(is),is=1,nsurf-1),&
       (stateOut(is),is=1,nsurf),&
       zenith_deg,azimuth,bulkalbedo,Fcld,&
       lai_wt,z0m,zdm,&
       ustar,l_mod,ra,ResistSurf,&
       Fc,&
       Fc_photo,Fc_respi,Fc_metab,Fc_traff,Fc_build,&
       qn1_SF,qn1_S,SnowAlb,&
       Qm,QmFreez,QmRain,swe,mwh,MwStore,chSnow_per_interval,&
       (SnowRemoval(is),is=1,2),&
       t2_C,q2_gkg,avU10_ms& ! surface-level diagonostics
       /)

  IF (snowUse==1) THEN
     dataOutSnow(ir,1:ncolumnsDataOutSnow,Gridiv)=(/REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),&               !2
          REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,&                                        !5
          SnowPack(1:nsurf),mw_ind(1:nsurf),Qm_melt(1:nsurf),&                                     !26
          Qm_rain(1:nsurf),Qm_freezState(1:nsurf),snowFrac(1:(nsurf-1)),&                          !46
          rainOnSnow(1:nsurf),&                                                                    !53
          qn1_ind_snow(1:nsurf),kup_ind_snow(1:nsurf),freezMelt(1:nsurf),&                         !74
          MeltWaterStore(1:nsurf),SnowDens(1:nsurf),&                                              !88
          snowDepth(1:nsurf),Tsurf_ind_snow(1:nsurf)/)                                             !102
  ENDIF

  !Calculate new snow fraction used in the next timestep if snowUse==1
  !Calculated only at end of each hour.
  !if (SnowFractionChoice==2.and.snowUse==1.and.it==23.and.imin==(nsh_real-1)/nsh_real*60) then
  !   do is=1,nsurf-1
  !      if ((snowPack(is)>0.and.mw_ind(is)>0)) then
  !         write(*,*) is,snowPack(is),snowD(is),mw_ind(is),snowFrac(is)!

  !         snowFrac(is)=SnowDepletionCurve(is,snowPack(is),snowD(is))
  !         write(*,*) snowFrac(is)
  !         pause
  !      elseif (snowPack(is)==0) then
  !         snowFrac(is)=0
  !      endif
  !   enddo
  !endif


  !write(*,*) DecidCap(id), id, it, imin, 'Calc - before translate back'
  !write(*,*) iy, id, it, imin, 'Calc - before translate back'
  !if(Gridiv==1)  write(*,*) iy, id, it, imin, HDD(id-1,5), HDD(id,5), HDD(id-1,6), HDD(id,6)
  !if(id==12) pause
  !write(*,*) ' '

  IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_TranslateBack...'
  CALL SUEWS_TranslateBack(Gridiv,ir,irMax)

  !  ! store water balance states of the day, by TS 13 Apr 2016
  !  IF ( ir==irMax ) THEN
  !       PRINT*, 'ending end of', id,it,imin
  !     stateDay(id,Gridiv,:)     = state(:)
  !     soilmoistDay(id,Gridiv,:) = soilmoist(:)
  !       PRINT*, 'state:', state
  !       PRINT*, 'soilmoist', soilmoist
  !       PRINT*, '********************************'
  !  END IF

  ! if ( id>10 ) then
  !   stop "stop to test state"
  ! end if

!!!if((id <=3 .or. id > 364).and. it == 0 .and. imin == 0) pause
!!!if((id <=3 .or. id > 364).and. it == 23 .and. imin == 55) pause
!!!if((id >=100 .or. id < 103).and. it == 0 .and. imin == 0) pause

  ! write(*,*) '------------'

END SUBROUTINE SUEWS_Calculations
