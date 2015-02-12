!This subroutine does the actual calculations of the SUEWS code (mostly old SUEWS_Temporal).
!Made by LJ and HW Oct 2014
!Gives in the grid ID (Gridiv) and number of line in the met forcing data to be analyzed (ir)
!Last modification
! HCW 15 Jan 2015
! Added switch OHMIncQF to calculate QS with (1, default) or without (0) QF added to QSTAR
! To do 
!      - add iy and imin to output files (may impact LUMPS_RunoffFromGrid)
!      - move OHMIncQF to RunControl.nml
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
  use moist
  use mod_z
  use mod_k
  use solweig_module


  IMPLICIT NONE

  integer :: Gridiv,ir,i,ih,iMB
  logical:: debug=.false.
  !real(kind(1d0)),DIMENSION(4)::GridFromFrac !Not maybe needed in the future
  real(kind(1d0))::idectime
  integer:: imon,iday,iyr,iseas
  integer:: reset=1!,i,iv,ih,id_in,it_in,&
            !SunriseTime,SunsetTime,errFileYes,ind5min=1
  real(kind(1d0)):: SnowDepletionCurve
  real(kind(1d0))::lai_wt
  integer::ind5min=1
  integer::irMax
  
  integer::OHMIncQF=1 	!Move to RunControl.nml

  real(kind(1d0)):: ggg

!==================================================================
!==================================================================

  !Translate all data to the variables used in the model calculations
  !write(*,*) '**Calling SUEWS_Translate from SUEWS_Calculations'
  !write(*,*) Gridiv, ir, iMB
  call SUEWS_Translate(Gridiv,ir,iMB)
  call RoughnessParameters ! Added by HCW 11 Nov 2014
  
  !NARP Config now done in SUEWS_Translate

  !Water flow from other grids needs to be considered likely here ??
  
!=============Get data ready for the qs calculation====================
 STPH1=0
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
 !! Should these be initialised in SUEWS_Translate instead??
 !Initialization for OAF's water bucket scheme
 ! LUMPS only (Loridan et al. (2012)
 !RAINRES = 0.
 !RAINBUCKET = 0.
 !E_MOD=0 !RAIN24HR = 0.;

!=====================================================================

 !! Does flow from other grids, temperature and LAI initialization need to happen here??
 !!calls to SOLWEIG_initial and CBL_initial?? - these have been moved to SUEWS_Program, but are not necessarily completed.
 
 !INITIALIZE VARIABLE FOR THE LOOP
 runoffAGveg=0
 runoffAGimpervious=0
 runoffSoil_per_interval=0
 runoffWaterBody=0
 ev_per_interval=0
 qe_per_interval=0
 st_per_interval=0
 dr_per_interval=0
 ch_per_interval=0
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

 ! Calculate sun position
 idectime=dectime-halftimestep! sun position at middle of timestep before
 call sun_position(year,idectime,timezone,lat,lng,alt,azimuth,zenith_deg)

 if(CBLuse>=1)then ! If CBL is used, calculated Temp_C and RH are replaced with the obs.  
   call CBL(i,iMB)
 endif
 
 !Call the dailystate routine to get surface characteristics ready
 call DailyState(Gridiv)

 if(LAICalcYes==0)then
   ! check -- this is going to be a problem as it is not for each vegtation class
   lai(id-1,:)=lai_obs
 endif

 !Calculation of density and other water related parameters
 call atmos_moist_lumps(avdens)


 !========Calculate water storage capacity in soil=========
 
 SoilMoistCap=0
 soilstate=0
 do is=1,nsurf-1 !No water body included
   soilmoistCap=soilstoreCap(is)*sfr(is)+soilMoistCap
   soilstate=soilmoist(is)*sfr(is)+soilstate
 enddo

 !write(*,*) 'ir =',ir
 if (ir==1) then  !Calculate initial smd
    smd=soilmoistcap-soilstate
 endif

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

   ALB(DecidSurf)=albDec(ID-1) !Change deciduous albedo

   call narp(alb_snow,qn1_SF,qn1_S)
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
 if (SOLWEIGout==1) then
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

 ! =================STORAGE HEAT FLUX=======================================
 if(QSChoice==1) then		!Use OHM to calculate QS
    if(OHMIncQF == 1) then	!Calculate QS using QSTAR+QF
       call OHM_v2015		
    elseif(OHMIncQF == 0) then  !Calculate QS using QSTAR 
      qn1=qn1_bup  
      call OHM_v2015
    endif   
 endif   
 
 ! For the purpose of turbulent fluxes, remove QF from the net all-wave radiation
 qn1=qn1_bup  ! remove QF from QSTAR
 
 ! -- qn1 is now QSTAR only
 
 if(AnthropHeatChoice>=1) then
    qf=QF_SAHP
 endif

 !==================Energy related to snow melting/freezing processes=======
 IF (snowUse==1)  then

   call MeltHeat(i)

   ! If snow on ground, no irrigation, so veg_fr same in each case
   !New fraction of vegetation.
   !IF(veg_type==1)THEN         ! area vegetated
      veg_fr=sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+&
                 sfr(GrassSurf)*(1-snowFrac(GrassSurf))+sfr(BSoilSurf)*(1-snowFrac(BSoilSurf))+&
                 sfr(WaterSurf)*(1-snowFrac(WaterSurf))

   !ELSEIF(veg_type==2)THEN     ! area irrigated
   !   !!!veg_fr=sfr(GrassISurf)*(1-snowFrac(GrassUSurf))
   !END IF

 END IF

 !==========================Turbulent Fluxes================================

 call LUMPS_QHQE !Calculate QH and QE from LUMPS

 if(debug)write(*,*)press_Hpa,psyc_hPA,i

 call WaterUse !Gives the external and internal water uses per timestep

 if(Precip>0) then   !Initiate rain data [mm]
   pin=Precip
 else
   pin=0
 endif

 ! Kinematic Heat flux (w'T'). Changed by LJ
 if(i==1) qh=qh_obs
 if ((qh==NAN.or.qh==0).and.h_mod/=0) then
  !if (h_mod/=NAN.and.h_mod<1000) then
   H=h_mod/(avdens*avcp)
 else
 ! if LUMPS has had a problem then we need a value
   H=(qn1*0.2)/(avdens*avcp)
 !call ErrorHint(37,'LUMPS unable to calculate realistic value.',h_mod, dectime, notUsedI) !comment out by SO
 endif

 !------------------------------------------------------------------

 call STAB_lumps(H,StabilityMethod,ustar,L_mod) !u* and monin-obukhov length out

 call AerodynamicResistance(RA,AerodynamicResistanceMethod,StabilityMethod,RoughLen_heat,&
                            ZZD,z0m,k2,AVU1,L_mod,Ustar,veg_fr,psyh)      !RA out

 if (snowUse==1) then
      call AerodynamicResistance(RAsnow,AerodynamicResistanceMethod,StabilityMethod,3,&
                                 ZZD,z0m,k2,AVU1,L_mod,Ustar,veg_fr,psyh)      !RA out
 endif

 call SurfaceResistance(id,it)   !qsc and surface resistance out
 call BoundaryLayerResistance


 sae=s_hPa*(qn1_SF+qf-qs)    !s - slope of svp vs t curve
 !qn1 changed to qn1_SF, lj in May 2013

 vdrc=vpd_hPa*avdens*avcp
 sp=s_hPa/psyc_hPa
 tlv=lv_J_kg/tstep_real
 e=sae+vdrc/ra

 ! write(*,*)e,sae,vdrc,ra,vpd_hPa,avdens,avcp,s_Hpa,qn1,qs,qf
 ! pause
 !if(debug)write(*,*)dectime,press_HPa,i
 !#######End of water vapour calculations############################

 !!Need to do between-grid water transfer
 !!Additional water from other grids
 !call RunoffFromGrid(GridFromFrac)

 !Initiate pipe water and runoff including pipes
 runoffPipes=addPipes
 runoff_per_interval=addPipes
 
 
 ! v2014b switched to 60/NSH min intervals here; v2014c already runs at this timestep
 
 !surplusWaterBody=0  !!This was set to zero outside the 60/nsh timestep loop - adjust later!!
 
 Surplus_evap=0  !Initiate evaporation surplus
 evap_5min=0
 ev=0
 qe=0
 ev_snow=0

 SoilStateOld=soilMoist !Initialize storages
 StateOld=state

 ! Calculate drainage for each soil subsurface (excluding water body)
 do is=1,nsurf-1
    call Drainage(surf(1,is),surf(2,is),surf(3,is),surf(4,is)) !per interval
 enddo
 ! Set drainage from water body to zero
 Drain(WaterSurf)=0  

 ! Distribute water within grid
 call ReDistributeWater
 
 
 !! NEED TO CHECK THROUGH ALL OF THIS. SNOWCALC SUBROUTINE MAY ALSO NEED ADAPTING!!
 !===============EVAPORATION and STORAGE for each surface===============

 do is=1,nsurf  
         
    if (snowCalcSwitch(is)==1) then
       call snowCalc(i) 
    else
       call Evap_SUEWS  !qe and ev out
                 
       call soilstore(surf(1,is))   !Soil store updates
              
    !! Correct to comment this out now??
    !!if(in==nsh) then   	! if requirement added by HCW 30/07/2014
       if (is.ne.WaterSurf) st_per_interval=st_per_interval+state(is)*sfr(is)!State on non-water area (LJ 10/2010)
    !!endif
                      
    ! Add evaporation from different surfaces to total evaporation
    if (is==BldgSurf.or.is==PavSurf) then
       ev_per_interval=ev_per_interval+((ev-Surplus_evap(is))*sfr(is))
       qe_per_interval=qe_per_interval+((ev-Surplus_evap(is))*sfr(is))*lv_J_kg
       evap_5min=evap_5min+((ev-Surplus_evap(is))*sfr(is))
    else
       ev_per_interval=ev_per_interval+(ev*sfr(is))
       qe_per_interval=qe_per_interval+(ev*sfr(is))*lv_J_kg
       evap_5min=evap_5min+(ev*sfr(is))
    endif
       
    ch_per_interval=ch_per_interval+(state(is)-stateOld(is))*sfr(is) !Update surface storage change. LJ
    
    ChangSnow(is)=0
    runoffSnow(is)=0
               
    endif           
        
 enddo  !end loop over surfaces
 
  
 
 ! At this point water has moved between the canopy and soilstorages of each surface.
 ! ==================================================================================  
 ! Now water is allowed to move between the surface stores
 call HorizontalSoilWater
           
 !Calculate soilmoisture state
 Soilmoist_state=0  	!Initial value 
 do is=1,nsurf-1 	!Summing all together (excluding water body)
    soilmoist_state=soilmoist(is)*sfr(is)+soilmoist_state
    if(soilmoist_state<0)then
       call errorHint(29,'subroutine SUEWS_Calculations[soilmoist_state<0],dectime,soilmoist_state,sfr(is)',&
                      dectime,soilmoist_state,int(sfr(is)))
       call errorHint(29,'subroutine SUEWS_Calculations[soilmoist_state<0],dectime,soilmoist(is),sfr(is)',&
                      dectime,soilmoist(is),int(sfr(is)))
 
    endif
 enddo  !end loop over surfaces
           
 if (SoilMoist_state>soilMoistCap) then  !What is this LJ 10/2010 - SM exceeds capacity, but where does extra go?HCW 11/2014
    SoilMoist_state=soilMoistCap
 endif
 
 ! Calculate soil moisture deficit for each time-interval.
 smd=soilMoistCap-soilmoist_state
 smd_nsurf=SoilstoreCap-soilmoist
 
 ! Update changes in soil storages to overall ch_per_interval. 
 ! - Needed as soil stores can change after horizontal water movements. 
 do is=1,nsurf
    ch_per_interval = ch_per_interval+(SoilMoist(is)-SoilStateOld(is))*sfr(is)
 enddo
           
 !Remove non-existing surface type from surface and soil outputs
 do is=1,nsurf
    if (sfr(is)<0.00001)then
       StateOut(is)=0
       smd_nsurfOut(is)=0
       runoffOut(is)=0
       runoffSoilOut(is)=0
    else
       StateOut(is)=State(is)
       smd_nsurfOut(is)=smd_nsurf(is)
       runoffOut(is)=runoff(is)
       runoffSoilOut(is)=runoffSoil(is)
    endif
 enddo  !end loop over surfaces
 
 !write(*,*) dectime
 
 !open(13,file='TestingFiveMin.txt',position="append")
 !write(13,*) id,in,dectime
 !close(13)
 
 !!Write out 5-min file - fix this later (needs allocating??)
 !Modified HCW 30/07/2014 so that output columns and header match (was 69 cols with extra columns for water)
 if(write5min==1) then           !   Save 5 min results to a file
 !   dataOut5min(ind5min,1:64)=(/real(id,kind(1D0)),real(in,kind(1D0)),dectime,pin,ext_wu,ev_per_interval,&
 !                               stateOut(1:nsurf),smd_nsurfOut(1:(nsurf-1)),drain(1:(nsurf-1)),runoffOut(1:(nsurf-1)),&
 !                               runoffsoilOut(1:(nsurf-1)),runoffSnow(1:(nsurf-1)),snowPack(1:nsurf),&
 !                               ChangSnow(1:nsurf),mw_ind(1:nsurf)/)
    ind5min = ind5min+1
 endif

 ! end of 60/NSH timestep used to be here (v2014b) 


 !======FINAL STEPS BEFORE WRITING OUT====================================
 !! CHECKSurfs What is this - why no Paved surf?? (HCW, 26 Jan 2015) 
 AdditionalWater=addWaterBody*sfr(WaterSurf)+addPipes+addImpervious*sfr(BldgSurf)+addveg*  &
                 (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))

 qeph=qe_per_interval/t_Interval * t_Interval/tstep_real !Calculate evaporation per interval          
           
 !Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
 !qh=(qn1+qf+QmRain+QmFreez)-(qeph+qs+Qm) 
 qh=(qn1+qf+QmRain)-(qeph+qs+Qm+QmFreez) 
     
 ! Deleted HCW 26 Jan 2015 (as ext_wu already for the whole study area)
 !!Calcuate external water use           
 !ext_wuhP=0
 !if(IrrFracConif>0.or.IrrFracDecid>0.or.IrrFracGrass>0)then                 
 !   ext_wuhP=ext_wu*sfr(GrassSurf)+wuhTrees*IrrTrees       
 !endif
      
 !Remove negative state     
 if(st_per_interval<0) st_per_interval=0 

 !     write hour results
 if(ResistSurf>9999) ResistSurf=9999

 if(abs(qh)>pNAN) qh=NAN
 if(abs(qe)>pNAN) qe=NAN
 if(abs(qs)>pNAN) qs=NAN
 if(abs(qeph)>pNAN) qeph=NAN
 if(abs(ch_per_interval)>pNAN) ch_per_interval=NAN
 if(abs(soilmoist_state)>pNAN) soilmoist_state=NAN
 if(abs(smd)>pNAN) smd=NAN

 !If measured smd is used, set components to -999 and smd output to measured one
 if (smd_choice>0) then
    smd_nsurfOut=NAN
    smd=xsmd
 endif
 
 ! Calculate areally-weighted LAI
 lai_wt=0
 do is=1,nvegsurf
    lai_wt=lai_wt+lai(id-1,is)*sfr(is+2) 
 enddo
 
 ! Calculate snowdepth from SWE
 do is=1,nsurf
    if (densSnow(is)/=0) then  
       SnowDepth(is) = SnowPack(is)*waterDens/densSnow(is) 
    endif
    ! Calculate overall snow water equivalent
    swe = swe + SnowDepth(is)*sfr(is)*snowFrac(is)
    MwStore = MwStore + MeltWaterStore(is)*sfr(is)*snowFrac(is)
 enddo

 !���������FILE WRITE SECTION���������������������������
 !Define the overall output matrix to be printed out step by step
 !! N.B. ext_wu now appears multiple times here - need to tidy columns
 !! Need to match up with LUMPS_RunoffFromGrid !!
 dataOut(ir,1:192,Gridiv)=(/real(iy,kind(1D0)),real(id,kind(1D0)),real(it,kind(1D0)),real(imin,kind(1D0)),dectime,avkdn,kup,& !7
                           ldown,lup,tsurf,qn1,h_mod,e_mod,qs,qf,qh,qeph,Precip,ext_wu,ev_per_interval,dr_per_interval,&     !21
                           ch_per_interval,st_per_interval,runoffSoil_per_interval,runoff_per_interval,runoffPipes,&           !26
                           runoffAGimpervious,runoffAGveg,runoffWaterBody,ra,ResistSurf,ustar,l_mod,&                          !33
                           (smd_nsurfOut(is),is=1,nsurf-1),(stateOut(is),is=1,nsurf),Fcld,soilmoist_state,smd,lai_wt,&         !50
                           FlowChange*sfr(WaterSurf),AdditionalWater,ext_wu,ext_wu,int_wu,qn1_SF,qn1_S,Qm,QmFreez,&        !59
                           QmRain,swe,MwStore,(SnowRemoval(is),is=1,2),chSnow_per_interval,&                                   !65
                           kup_ind(1:nsurf),lup_ind(1:nsurf),tsurf_ind(1:nsurf),qn1_ind(1:nsurf),&                             !66-93
                           SnowPack(1:nsurf),mwh,mw_ind(1:nsurf),Qm_melt(1:nsurf),&                                            !94-115
                           Qm_rain(1:nsurf),Qm_freezState(1:nsurf),snowFrac(1:(nsurf-1)),alb_snow,rainOnSnow(1:nsurf),&        !116-143
                           qn1_ind_snow(1:nsurf),kup_ind_snow(1:nsurf),freezMelt(1:nsurf),MeltWaterStore(1:nsurf),&            !144-171
                           densSnow(1:nsurf),snowDepth(1:nsurf),Tsurf_ind_snow(1:nsurf)/)                                      !172-192

  !Calculate new snow fraction used in the next timestep if snowUse==1
  !This is done each hour?? Location likely needs to be changes
 
 if (SnowFractionChoice==2.and.snowUse==1) then
    do is=1,nsurf-1
       if (snowPack(is)>0.and.mw_ind(is)>0) then
          snowFrac(is)=SnowDepletionCurve(is,snowPack(is),snowD(is),SnowLimPaved,SnowLimBuild)
       elseif (snowPack(is)==0) then
          snowFrac(is)=0
       endif
    enddo
 endif      
 
 call SUEWS_TranslateBack(Gridiv,ir,iMB,irMax) 
 
 !write(*,*) 'In SUEWS_Calculations'
 !write(*,*) 'imin',imin             
 !write(*,*) 'it',it             
 !write(*,*) 'id',id
 !write(*,*) 'iy',iy
 
 !write(*,*) SnowPack
 !write(*,*) SnowFrac
 !write(*,*) densSnow
 !write(*,*) State
 !write(*,*) SoilMoist
 !write(*,*) ModelOutputData(ir,cMOD_SnowPack(1), Gridiv)
 !write(*,*) ModelOutputData(ir,cMOD_SnowPack(2), Gridiv) 
 !write(*,*) ModelOutputData(ir,cMOD_SnowPack(4), Gridiv)  
  
 write(*,*) '------------'


!GridFromFrac- close files??

!!accumulate??

!SnowPack_grid=SnowPack   !!is this needed??

 END SUBROUTINE SUEWS_Calculations