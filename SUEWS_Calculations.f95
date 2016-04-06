!This subroutine does the actual calculations of the SUEWS code (mostly old SUEWS_Temporal).
!Made by LJ and HW Oct 2014
!Gives in the grid ID (Gridiv) and number of line in the met forcing data to be analyzed (ir)
!Last modification
! TS 09 Mar 2016
!  Added AnOHM subroutine to calculate heat storage
!Last modification:
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

  use data_in
  use time
  use NARP_MODULE
  use defaultNotUsed
  use allocateArray
  use sues_data
  use snowMod
  use gis_data
  use initial
  use moist
  use mod_z
  use mod_k
  use solweig_module


  IMPLICIT NONE

  integer:: Gridiv,ir,i,ih,iMB
  logical:: debug=.false.
  real(kind(1d0)):: idectime
  real(kind(1d0)):: SnowDepletionCurve
  real(kind(1d0)):: lai_wt
  integer:: irMax
    
!==================================================================
!==================================================================
 
  !Translate all data to the variables used in the model calculations
  call SUEWS_Translate(Gridiv,ir,iMB)
  call RoughnessParameters ! Added by HCW 11 Nov 2014

!=============Get data ready for the qs calculation====================
 if(NetRadiationChoice==0) then !Radiative components are provided as forcing
   !avkdn=NAN                  !Needed for resistances for SUEWS.
   ldown=NAN
   lup=NAN
   kup=NAN
   tsurf=NAN
   lup_ind=NAN
   kup_ind=NAN
   tsurf_ind=NAN
   qn1_ind=NAN
   Fcld=NAN
 endif

 if(ldown_option==1) then
   Fcld=NAN
 endif
!=====================================================================
 ! Initialisation for OAF's water bucket scheme
 ! LUMPS only (Loridan et al. (2012))
 RAINRES = 0.
 RAINBUCKET = 0.
 E_MOD=0 !RAIN24HR = 0.;

!=====================================================================
 
 !Initialize variables calculated at each 5-min timestep
 runoffAGveg=0
 runoffAGimpervious=0
 runoffWaterBody=0
 runoffSoil_per_tstep=0
 runoffSoil_per_interval=0
 chSnow_per_interval=0
 mwh = 0       !Initialize snow melt and heat related to snowmelt
 fwh = 0
 Qm = 0        !Heat related to melting/freezing
 QmFreez = 0
 QmRain =0
 Mw_ind = 0
 SnowDepth = 0
 zf = 0
 deltaQi = 0
 swe = 0
 MwStore = 0
 WaterHoldCapFrac=0

 qh = -999 ! Added HCW 26 Feb 2015
 H = -999 ! Added HCW 26 Feb 2015
 
 ! Calculate sun position
 idectime=dectime-halftimestep! sun position at middle of timestep before
 call sun_position(year,idectime,timezone,lat,lng,alt,azimuth,zenith_deg)

 if(CBLuse>=1)then ! If CBL is used, calculated Temp_C and RH are replaced with the obs.  
   call CBL(ir,iMB)   !ir=1 indicates first row of each met data block
 endif
 
 !Call the dailystate routine to get surface characteristics ready
 call DailyState(Gridiv)
  
 if(LAICalcYes==0)then
   lai(id-1,:)=lai_obs ! check -- this is going to be a problem as it is not for each vegetation class
 endif

 !Calculation of density and other water related parameters
 call atmos_moist_lumps(avdens)


 !======== Calculate soil moisture =========
 SoilMoistCap=0   !Maximum capacity of soil store [mm] for whole surface
 soilstate=0      !Area-averaged soil moisture [mm] for whole surface
 do is=1,nsurf-1   !No water body included   
   soilmoistCap=soilMoistCap+(soilstoreCap(is)*sfr(is)/NonWaterFraction)
   soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
 enddo

 !If loop removed HCW 26 Feb 2015
 !if (ir==1) then  !Calculate initial smd
    smd=soilmoistCap-soilstate
 !endif

 ! Calculate soil moisture for vegetated surfaces only (for use in surface conductance)
 vsmd=0
 do is=ConifSurf,GrassSurf  !Vegetated surfaces only
    vsmd=vsmd+(soilstoreCap(is) - soilmoist(is))*sfr(is)/(sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))
    !write(*,*) is, vsmd, smd
 enddo

 ! ===================NET ALLWAVE RADIATION================================
 if(NetRadiationChoice>0)then

   if (snowUse==0) snowFrac=snow_obs

   if(ldown_option==1) then !Observed ldown provided as forcing
     ldown=ldown_obs
   else
     ldown=-9              !to be filled in NARP
   endif

   if(ldown_option==2) then !observed cloud fraction provided as forcing
     fcld=fcld_obs
   endif

   !write(*,*) DecidCap(id), id, it, imin, 'Calc - near start'
   
   ! Update variables that change daily and represent seasonal variability
   alb(DecidSurf)=albDec(id) !Change deciduous albedo
   surf(6,DecidSurf)=DecidCap(id)  !Change current storage capacity of deciduous trees
   ! Change EveTr and Grass albedo too
   alb(ConifSurf)=albEveTr(id)
   alb(GrassSurf)=albGrass(id)
   
   call NARP(SnowAlb,qn1_SF,qn1_S)
   !Temp_C,kclear,fcld,dectime,avkdn,avRH,qn1,kup,ldown,lup,tsurf,&
   !AlbedoChoice,ldown_option,Press_hPa,Ea_hPa,qn1_obs,&
   !zenith_deg,netRadiationChoice,
 else
   snowFrac = snow_obs
   qn1=qn1_obs
   qn1_sf=qn1_obs
   qn1_s=qn1_obs
 endif

! ===================SOLWEIG OUTPUT ========================================
 if (SOLWEIGuse==1) then
     if (OutInterval==imin) then
         if (RunForGrid==-999) then
             call SOLWEIG_2014a_core(iMB)
             SolweigCount=SolweigCount+1
         else
             if (Gridiv == RunForGrid) then
                 call SOLWEIG_2014a_core(iMB)
                 SolweigCount=SolweigCount+1
             endif             
         endif
     endif
 else
   SOLWEIGpoi_out=0
 endif

 ! ===================ANTHROPOGENIC HEAT FLUX================================
 ih=it-DLS
 if(ih<0) ih=23
 
 if(AnthropHeatChoice==1) then
    call SAHP_1_v2015(qf_sahp,id,ih,imin)
    qn1_bup=qn1
    qn1=qn1+QF_SAHP
 elseif(AnthropHeatChoice==2) then
    call SAHP_2_v2015(qf_sahp,id,ih,imin)
    qn1_bup=qn1
    qn1=qn1+QF_SAHP
 else
    qn1_bup=qn1
    qn1=qn1+qf
 endif
 
 ! -- qn1 is now QSTAR+QF (net all-wave radiation + anthropogenic heat flux)
 ! -- qn1_bup is QSTAR only
 if(AnthropHeatChoice>=1) then
    qf=QF_SAHP
 endif
 ! =================STORAGE HEAT FLUX=======================================
 if(QSChoice==1) then           !Use OHM to calculate QS
    if(OHMIncQF == 1) then      !Calculate QS using QSTAR+QF
      call OHM_v2015(Gridiv)
    elseif(OHMIncQF == 0) then  !Calculate QS using QSTAR
      qn1=qn1_bup  
      call OHM_v2015(Gridiv)
    endif   
 endif   
  
 if (QSChoice==3) then 	! use AnOHM to calculate QS
 	call AnOHM_v2016(Gridiv)
 end if   
 
 ! For the purpose of turbulent fluxes, remove QF from the net all-wave radiation
 qn1=qn1_bup  !Remove QF from QSTAR
 
 ! -- qn1 is now QSTAR only
 


 !==================Energy related to snow melting/freezing processes=======
 IF (snowUse==1)  then

   call MeltHeat

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
 
 call LUMPS_QHQE !Calculate QH and QE from LUMPS
 if(debug)write(*,*)press_Hpa,psyc_hPA,i
 
 call WaterUse !Gives the external and internal water uses per timestep

 if(Precip>0) then   !Initiate rain data [mm]
   pin=Precip
 else
   pin=0
 endif

 
 ! Get first estimate of sensible heat flux. Modified by HCW 26 Feb 2015  
 ! Calculate kinematic heat flux (w'T') from sensible heat flux [W m-2] from observed data (if available) or LUMPS 
 if(qh_obs/=NAN) then   !if(qh_obs/=NAN) qh=qh_obs   !Commented out by HCW 04 Mar 2015
    H=qh_obs/(avdens*avcp)  !Use observed value
 else
    if(h_mod/=NAN) then 
       H=h_mod/(avdens*avcp)   !Use LUMPS value
    else   
       H=(qn1*0.2)/(avdens*avcp)   !If LUMPS has had a problem, we still need a value
       call ErrorHint(38,'LUMPS unable to calculate realistic value for H_mod.',h_mod, dectime, notUsedI)
    endif   
 endif
  
 !------------------------------------------------------------------

 call STAB_lumps(H,StabilityMethod,ustar,L_mod) !u* and monin-obukhov length out

 call AerodynamicResistance(RA,AerodynamicResistanceMethod,StabilityMethod,RoughLen_heat,&
                            ZZD,z0m,k2,AVU1,L_mod,Ustar,VegFraction,psyh)      !RA out

 if (snowUse==1) then
    call AerodynamicResistance(RAsnow,AerodynamicResistanceMethod,StabilityMethod,3,&
                               ZZD,z0m,k2,AVU1,L_mod,Ustar,VegFraction,psyh)      !RA out
 endif

 call SurfaceResistance(id,it)   !qsc and surface resistance out
 call BoundaryLayerResistance


 sae=s_hPa*(qn1_SF+qf-qs)    !s_haPa - slope of svp vs t curve. qn1 changed to qn1_SF, lj in May 2013
 vdrc=vpd_hPa*avdens*avcp
 sp=s_hPa/psyc_hPa
 numPM=sae+vdrc/ra

 !write(*,*) numPM, sae, vdrc/ra, s_hPA+psyc_hPa, NumPM/(s_hPA+psyc_hPa)

 !=====================================================================
 !========= Water balance calculations ================================
 ! Needs to run at small timesteps (i.e. minutes)
 ! Previously, v2014b switched to 60/NSH min intervals here
 ! Now whole model runs at a resolution of tstep
 
 ! Initialise water balance variables
 qe=0
 ev=0           
 ev_snow=0
 SurplusEvap=0
 evap=0
 chang=0
 runoff=0
 runoffSoil=0
 surplusWaterBody=0    
 
 ! Added by HCW 13 Feb 2015
 qe_per_tstep=0     ![W m-2]
 ev_per_tstep=0 
 drain_per_tstep=0
 surf_chang_per_tstep=0
 tot_chang_per_tstep=0
 state_per_tstep=0
 NWstate_per_tstep=0
 runoff_per_tstep=0
 
 ! Retain previous surface state and soil moisture state
 stateOld = state           !State of each surface [mm] for the previous timestep
 soilmoistOld = soilmoist   !Soil moisture of each surface [mm] for the previous timestep
  
 !============= Grid-to-grid runoff =============
 ! Calculate additional water coming from other grids 
 ! i.e. the variables addImpervious, addVeg, addWaterBody, addPipes
 !call RunoffFromGrid(GridFromFrac)  !!Need to code between-grid water transfer

 ! Sum water coming from other grids (these are expressed as depths over the whole surface)
 AdditionalWater=addPipes+addImpervious+addVeg+addWaterBody  ![mm]
 
 ! Initialise runoff in pipes
 runoffPipes=addPipes   !Water flowing in pipes from other grids. No need for scaling??
 !! CHECK p_i
 runoff_per_interval=addPipes !pipe plor added to total runoff.
  
 !================== Drainage ===================
 ! Calculate drainage for each soil subsurface (excluding water body)
 do is=1,nsurf-1
    call Drainage(surf(6,is),surf(2,is),surf(3,is),surf(4,is))
    !HCW added and changed to surf(6,is) here 20 Feb 2015

    drain_per_tstep=drain_per_tstep+(drain(is)*sfr(is)/NonWaterFraction)   !No water body included
 enddo

 drain(WaterSurf) = 0  ! Set drainage from water body to zero
 
 ! Distribute water within grid, according to WithinGridWaterDist matrix (Cols 1-7)
 call ReDistributeWater   !Calculates AddWater(is)
  
 !======== Evaporation and surface state ========
 do is=1,nsurf   !For each surface in turn      
    if (snowCalcSwitch(is)==1) then
       if (sfr(is)/=0) then
          call snowCalc
       else
          snowFrac(is)=0
          SnowDens(is)=0
          SnowPack(is)=0
       endif
    else
       call Evap_SUEWS   !Calculates ev [mm]             
       call soilstore    !Surface water balance and soil store updates (can modify ev, updates state)

       evap(is) = ev  !Store ev for each surface
       
       ! Sum evaporation from different surfaces to find total evaporation [mm]
       ev_per_tstep=ev_per_tstep+evap(is)*sfr(is)

       ! Sum change from different surfaces to find total change to surface state
       surf_chang_per_tstep=surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)
       ! Sum runoff from different surfaces to find total runoff
       runoff_per_tstep=runoff_per_tstep+runoff(is)*sfr(is)      
       ! Calculate total state (including water body)
       state_per_tstep=state_per_tstep+(state(is)*sfr(is))
       ! Calculate total state (excluding water body)
       if (is.ne.WaterSurf) NWstate_per_tstep=NWstate_per_tstep+(state(is)*sfr(is)/NonWaterFraction)
              
       ChangSnow(is)=0
       runoffSnow(is)=0
       
    endif                  
 enddo  !end loop over surfaces

 
 ! Convert evaporation to latent heat flux [W m-2]
 qe_per_tstep=ev_per_tstep*tlv
 qeOut=qe_per_tstep
 
 !============ Sensible heat flux ===============
 ! Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
 qh=(qn1+qf+QmRain)-(qeOut+qs+Qm+QmFreez)     !qh=(qn1+qf+QmRain+QmFreez)-(qeOut+qs+Qm)
 
 
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
 runoffAGveg_m3 = runoffAGveg/1000 *SurfaceArea
 runoffWaterBody_m3 = runoffWaterBody/1000 *SurfaceArea
 runoffPipes_m3 = runoffPipes/1000 *SurfaceArea
 
  
 !=== Horizontal movement between soil stores ===
 ! Now water is allowed to move horizontally between the soil stores
 call HorizontalSoilWater
           
 !========== Calculate soil moisture ============
 soilstate=0       !Area-averaged soil moisture [mm] for whole surface
 do is=1,nsurf-1   !No water body included   
    soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
    if (soilstate<0) then
       call ErrorHint(62,'SUEWS_Calculations: total soilstate < 0 (just added surface is) ',soilstate,NotUsed,is)
    elseif (soilstate>SoilMoistCap) then
       call ErrorHint(62,'SUEWS_Calculations: total soilstate > capacity (just added surface is) ',soilstate,NotUsed,is)   
       !SoilMoist_state=soilMoistCap !What is this LJ 10/2010 - SM exceeds capacity, but where does extra go?HCW 11/2014
    endif
 enddo  !end loop over surfaces          
 ! Calculate soil moisture deficit
 smd=soilMoistCap-soilstate   !One value for whole surface
 smd_nsurf=SoilstoreCap-soilmoist   !smd for each surface
 
 ! Soil stores can change after horizontal water movements
 ! Calculate total change in surface and soil state
 tot_chang_per_tstep = surf_chang_per_tstep   !Change in surface state
 do is=1,(nsurf-1)   !No soil for water surface (so change in soil moisture is zero)
    tot_chang_per_tstep = tot_chang_per_tstep + ((SoilMoist(is)-SoilMoistOld(is))*sfr(is))   !Add change in soil state
 enddo
 
 !=====================================================================
 !====================== Prepare data for output ======================
 
 ! Remove non-existing surface type from surface and soil outputs   !Commented out by HCW 04 Mar 2015 as should not occur
 !do is=1,nsurf
 !   if (sfr(is)<0.00001) then
 !      stateOut(is)=0
 !      smd_nsurfOut(is)=0
 !      runoffOut(is)=0
 !      runoffSoilOut(is)=0
 !   else
 !      stateOut(is)=state(is)
 !      smd_nsurfOut(is)=smd_nsurf(is)
 !      runoffOut(is)=runoff(is)
 !      runoffSoilOut(is)=runoffSoil(is)
 !   endif
 !enddo
 
 ! Remove negative state   !ErrorHint added, then commented out by HCW 16 Feb 2015 as should never occur
 !if(st_per_interval<0) then
 !   call ErrorHint(63,'SUEWS_Calculations: st_per_interval < 0',st_per_interval,NotUsed,NotUsedI)   
 !   !st_per_interval=0 
 !endif   
 
 ! Set limit on ResistSurf
 if(ResistSurf>9999) ResistSurf=9999
     
 ! Set NA values   !!Why only these variables??  !!ErrorHints here too??
 if(abs(qh)>pNAN) qh=NAN
 if(abs(qeOut)>pNAN) qeOut=NAN
 if(abs(qs)>pNAN) qs=NAN
 if(abs(ch_per_interval)>pNAN) ch_per_interval=NAN
 if(abs(surf_chang_per_tstep)>pNAN) surf_chang_per_tstep=NAN
 if(abs(tot_chang_per_tstep)>pNAN) tot_chang_per_tstep=NAN
 if(abs(soilstate)>pNAN) soilstate=NAN
 if(abs(smd)>pNAN) smd=NAN

 ! If measured smd is used, set components to -999 and smd output to measured one
 if (smd_choice>0) then
     smd_nsurf=NAN
    !smd_nsurfOut=NAN
    smd=xsmd
 endif
 
 ! Calculate areally-weighted LAI
 if(iy == (iy_prev_t+1) .and. (id-1) == 0) then   !Check for start of next year and avoid using lai(id-1) as this is at the start of the year 
   lai_wt=0  
   do is=1,nvegsurf
       lai_wt=lai_wt+lai(id_prev_t,is)*sfr(is+2) 
   enddo        
 else
   lai_wt=0
   do is=1,nvegsurf
       lai_wt=lai_wt+lai(id-1,is)*sfr(is+2) 
   enddo
 endif

 
 !=====================================================================
 !====================== Write out files ==============================
 !Define the overall output matrix to be printed out step by step
 dataOut(ir,1:ncolumnsDataOut,Gridiv)=(/real(iy,kind(1D0)),real(id,kind(1D0)),real(it,kind(1D0)),real(imin,kind(1D0)),dectime,&   !5
        avkdn,kup,ldown,lup,tsurf,qn1,h_mod,e_mod,qs,qf,qh,qeOut,&                                       !17
        precip,ext_wu,ev_per_tstep,drain_per_tstep,&                                                     !21
        state_per_tstep,NWstate_per_tstep,surf_chang_per_tstep,tot_chang_per_tstep,&                     !25
        runoff_per_tstep,runoffSoil_per_tstep,runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,&  !31
        AdditionalWater,FlowChange/nsh_real,int_wu,wu_EveTr,wu_DecTr,wu_Grass,&                          !37
        ra,ResistSurf,ustar,l_mod,Fcld,&                                                                 !42
        soilstate,smd,(smd_nsurf(is),is=1,nsurf-1),(state(is),is=1,nsurf),&                              !57
        lai_wt,z0m,zdm,&                                                                                 !60
        qn1_SF,qn1_S,Qm,QmFreez,QmRain,swe,mwh,MwStore,(SnowRemoval(is),is=1,2),chSnow_per_interval,&    !71
        SnowAlb/)                                                                                        !72

 if (snowUse==1) then
    dataOutSnow(ir,1:ncolumnsDataOutSnow,Gridiv)=(/real(iy,kind(1D0)),real(id,kind(1D0)),&               !2
                real(it,kind(1D0)),real(imin,kind(1D0)),dectime,&                                        !5
                SnowPack(1:nsurf),mw_ind(1:nsurf),Qm_melt(1:nsurf),&                                     !26
                Qm_rain(1:nsurf),Qm_freezState(1:nsurf),snowFrac(1:(nsurf-1)),&                          !46
                rainOnSnow(1:nsurf),&                                                                    !53
                qn1_ind_snow(1:nsurf),kup_ind_snow(1:nsurf),freezMelt(1:nsurf),&                         !74
                MeltWaterStore(1:nsurf),SnowDens(1:nsurf),&                                              !88
                snowDepth(1:nsurf),Tsurf_ind_snow(1:nsurf)/)                                             !102
 endif
 
 !write(*,*) DecidCap(id), id, it, imin, 'Calc - before translate back'
 !write(*,*) iy, id, it, imin, 'Calc - before translate back'
 !if(Gridiv==1)  write(*,*) iy, id, it, imin, HDD(id-1,5), HDD(id,5), HDD(id-1,6), HDD(id,6)
 !if(id==12) pause
 !write(*,*) ' '

 call SUEWS_TranslateBack(Gridiv,ir,irMax)

 !!!if((id <=3 .or. id > 364).and. it == 0 .and. imin == 0) pause
 !!!if((id <=3 .or. id > 364).and. it == 23 .and. imin == 55) pause
 !!!if((id >=100 .or. id < 103).and. it == 0 .and. imin == 0) pause
 
 ! write(*,*) '------------'

 END SUBROUTINE SUEWS_Calculations
