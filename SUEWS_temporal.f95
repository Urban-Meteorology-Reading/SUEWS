! vanu_holt_hanna.f   sg June 94
! sg aug 99 F90
!    
! PROGRAM to compute input for van Ulden & Holtslag (1983) JAM
! for qh & qe  alpha and  beta
! http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html
!
! V4.0 SUES added (Grimmond and Oke 1991, RR)
    
! VERSION 4.1 
! Thomas Loridan 05/2008
! NARP_MODULE modified to allow for cloud fraction parameterization at night

! 4.2: change output to 15 and 30 min also  --- sg june 2008
! 5.0: create more subroutines 
! add water balance model (Grimmond and Oke 1986 WRR)
!
! V5.1
! Different longwave down options (ldown_option) - T. Loridan - June 2009
! 1 - Ldown Observed (add data as last column in met file)  
! 2 - Ldown modelled from observed FCLD (add data as last column in met file)
! 3 - Ldown modelled from FCLD(RH,TA)
! 4 - Ldown modelled from FCLD(Kdown); i.e. day FCLD only 
!     cloud fraction is kept constant throught the night (Offerle et al. 2003, JAM) 
! 5 - Option 3 at night and 4 during the day (might cause discontinuities in Ldown)
!
! Added a new met format (6) to read and convert IUMIP data D. Young - Aug 2009

!V5.1 Update to GIS format to include only land cover fraction
!Output to different periods has been removed for the time being 

! December 2009 - Simple Anthropogenic Heat Parameterization (SAHP) routine added.
!                 input parameters are to be specified in the .sahp file
!
!               - Simple water bucket from Offerle (PhD thesis, 2003) also added
!V5.2 
!This version is a simplified version of LUMPS allowing only default meteorological
!input format (InputMetformat=0). Also SUES and anthropogenic heat emissions 
!are excluded. L. Jarvi (Jan 2010)

!V6.0=SUES-2: 
!SUES and UWB parts are now working. 
!NARP calculated from different surfaces, Simple WU model, soil storage flows
!L.Jarvi (11 Jun)

!Now subroutine for main program SUEWS_v1_0. L. Jarvi (Nov 2010)
!	Flows from other grids, 7th subsurface included, LAI calculation modified, 
!
!sg feb 2012 - number of surfaces - allocateArray  module added
! start at midnight
! daily calculations done each day - 
!lj June 2012 - development of snow to the model has started
!fl May 2014 - coupling with SOLWEIG has started
! hcw Jul 2014 - state on non-water area is now calculated only for the last timestep  
! hcw Jul 2014 - Corrected 5-min output file columns to match header (64 cols)

!----------------------------------------------------------------------------------
 
subroutine SUEWS_temporal(GridName,GridFrom,GridFromFrac,iyr,errFileYes,SnowPack_grid)
  use allocateArray  !module_LUMPS_constants,f90   
  use gas           ! module_LUMPS_constants,f90
  use mod_grav	    ! module_LUMPS_constants,f90
  use mod_z	        ! module_LUMPS_constants,f90
  use mod_k	        ! module_LUMPS_constants,f90
  use moist
  use ohm_calc
  use data_in	    ! module_LUMPS_constants,f90  
  use gis_data      ! LUMPS_gis_read.f90
  use SUES_data     ! LUMPS_Inital.f95
  use run_info      ! run_control_v2_1.f90
  use NARP_MODULE   ! LUMPS_NARP6_0.f90
  use time          ! LUMPS_Module_constants.f95
  use defaultnotUsed
  use VegPhenogy
  use snowMod
  use solweig_module
  
  implicit none 
  
  !Grid information 
  real(kind(1d0)),DIMENSION(4)::GridFromFrac
  character(len=15)::GridName
  character(len=15),dimension(2,4)::GridFrom
  character(len=15)::year_txt

  real (kind(1d0)),dimension(nsurf)::SnowPack_grid
  
  !Other variables 
  logical:: debug=.false.                            
  integer:: imon,iday,iyr,iseas,reset=1,i,iv,ih,id_in,it_in, SunriseTime,SunsetTime,errFileYes,ind5min=1
            
  real(kind(1d0))::lai_wt,dectime_nsh,SnowDepletionCurve,idectime
                   
  character(len=100)::FileNameOld,str2

  !Variables related to NARP
   !Radiation balance components for different surfaces
  real(kind(1D0))::NARP_ALB_is,NARP_EMIS_is,snowFracTot!,qn1_cum,kup_cum,lup_cum,tsurf_cum

  
 !Initialize the model (reading nml files, defining filepaths, printing filechoices)
 call OHMinitialize
 !===========================NARP CONFIG=============================================
 if(NetRadiationChoice>0)then ! I don't think this is needed anymore (FL)
    	call NARP_CONFIG(LAT,LNG,YEAR,TIMEZONE,ALB_SNOW,EMIS_SNOW,TRANS_SITE,Interval,ldown_option)
  		!This for the snow cover fractions
 		!Initiate NARP anyway in order to get surface temperatures 
 		!call NARP_CONFIG(LAT,LNG,YEAR,TIMEZONE,ALB,EMIS,ALB_SNOW,EMIS_SNOW,TRANS_SITE,INTERVAL,ldown_option)
 endif
 finish=.false.


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
!=====================================================================================
 !Initialization for OAF's water bucket scheme
 ! LUMPS only (Loridan et al. (2012)
 RAINRES = 0.
 RAINBUCKET = 0.
 E_MOD=0 !RAIN24HR = 0.;

!======Open files of other grids if flow exists=======================================
!Made by LJ 10/2010
!Read the additional water from surroundings grids here.4 possible options.
do is=1,4
   if (GridFromFrac(is)/=0) then!If runoff from other surfaces exists read data

   !File identifier
   if (is==1) then
       lfnOld=486
   elseif (is==2) then
       lfnOld=487
   elseif (is==3) then
       lfnOld=488
   elseif (is==4) then
       lfnOld=489
   endif
      
   Write(str2,'(i2)') interval/60
   Write(year_txt,'(i4)') int(year)
           
   FileNameOld=trim(FileOutputPath)//trim(GridFrom(1,is))//trim(year_txt)//'_'//trim(adjustl(str2))//'.txt'

   open(lfnOld,file=trim(FileNameOld),err=319,iostat=ios_out)
   call skipHeader(lfnOld,5)
   endif
 enddo
! CHECK Where this needs to actually occur at the moment just put here
  
! ! temperature  and LAI initialization
! need to think about what happens if LAI not calcualated


  if((CBLuse==1).or.(CBLuse==2)) call CBL_initial
  if(SOLWEIGout==1) call SOLWEIG_initial
    
 !=================================================================================
 !=========INTERVAL LOOP WHERE THE ACTUAL MODEL CALCULATIONS HAPPEN================
 !=================================================================================
 ! assume starting at midnight
 ! first couple of time steps will not be correct for Qs 
 
  do i=1,nlines  !367*48 !want more than it will be old: 4000 24 h days or 1000 96 *15 days!NPeriodsPerYear
   
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
    mw_ind = 0    
    SnowDepth = 0
    zf = 0
    deltaQi = 0
    swe = 0
    MwStore = 0
    WaterHoldCapFrac=0

    call ConvertMetData(i) !Get correct forcing data for each timestep

    if(finish) then
      write(*,*)id,it 
      if(finish)exit
    endif
    
    ! Calculate sun position
    idectime=dectime-halftimestep! sun position at middle of timestep before
    call sun_position(year,idectime,timezone,lat,lng,alt,azimuth,zenith_deg)
    

    if(CBLuse>=1)then
        call CBL(i)
    endif


    ! If GISInput Varies
    if(GISInputType==4)then
      id_in=id   ! sg  - gis data can be missing - so now check 
      it_in=it
       call read_gis(finish)
       if(id/=id_in.or.it/=it_in) call ErrorHint(21,FileGIS,real(id,kind(1d0)),real(it,kind(1d0)),it_in)
       if(finish)exit
    endif
    if(z0_method>1) call RoughnessParameters(id-1)   ! this will vary with porosity even with fixed sFr
    
    !=============DAILY LAI AND WATER USE====================================================
    !Calculate LAI and temperature related 
    call DailyState   

    if(LAICalcYes==0)then
      ! check -- this is going to be a problem as it is not for each vegtation class
       lai(id-1,:)=lai_hr        
    endif 
    
	!Calculation of density and other water related parameters
    call atmos_moist_lumps(avdens)

    
    !========Calculate water storage capasity in soil=========
    SoilMoistCap=0
    soilstate=0
  
    do is=1,nsurf-1 !No water body included
       soilmoistCap=soilstoreCap(is)*sfr(is)+soilMoistCap
       soilstate=soilmoist(is)*sfr(is)+soilstate
    enddo
    
    if (i==1) then  !Calculate initial smd
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
        call SOLWEIG_2014a_core(i)
    else
        SOLWEIGpoi_out=0
    endif  
    
    ! ===================ANTHROPOGENIC HEAT FLUX================================
    ih=it-DLS
 
    if(ih<0)then
         ih=23
    endif
    
    if(AnthropHeatChoice==1)then
       call SAHP(qF_sahp,ih,id)        
       qn1_bup=qn1
       qn1=qn1+QF_SAHP
       
    elseif(AnthropHeatChoice==2)then
       call SAHP_2(qf_sahp,ih,id)
        
       qn1_bup=qn1
       qn1=qn1+QF_SAHP
    else
       qn1_bup=qn1
       qn1=qn1+qf
    endif
    
    ! =================STORAGE HEAT FLUX=======================================
    if(QSChoice==1) call OHMnew
    
    
    !For the purpose of turbulent fluxes, remove QF from the net allwave radiation
    qn1=qn1_bup  ! remove QF from QSTAR
    if(AnthropHeatChoice>=1) then 
       qf=QF_SAHP
    endif

    !==================Energy related to snow melting/freezing processes=======
    if (snowUse==1)  then
    
      call MeltHeat(i)

      !New fraction of vegetation
      IF(veg_type==1)THEN         ! area vegetated
         veg_fr=sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+&
                sfr(GrassUSurf)*(1-snowFrac(GrassUSurf))+sfr(GrassISurf)*(1-snowFrac(GrassISurf))+&
                sfr(WaterSurf)*(1-snowFrac(WaterSurf))

      ELSEIF(veg_type==2)THEN     ! area irrigated
            veg_fr=sfr(GrassISurf)*(1-snowFrac(GrassUSurf))
      END IF
    endif
   
    !==========================Turbulent Fluxes================================

    call LUMPS_QHQE !Calculate QH and QE from LUMPS
    
    if(debug)write(*,*)press_Hpa,psyc_hPA,i

    call WaterUse !Gives the external and internal water uses per timestep

    if(Precip_hr>0)then  !Initiate rain data (per 300 s)
       pin=Precip_hr/nsh 
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

      !Calculate available energy times s
      
      sae=s_hPa*(qn1_SF+qf-qs)    !s - slope of svp vs t curve
                                  !qn1 changed to qn1_SF, lj in May 2013
    
      vdrc=vpd_hPa*avdens*avcp 
      sp=s_hPa/psyc_hPa
      tlv=lv_J_kg/tstep
      e=sae+vdrc/ra
     
      ! write(*,*)e,sae,vdrc,ra,vpd_hPa,avdens,avcp,s_Hpa,qn1,qs,qf
      ! pause
      !if(debug)write(*,*)dectime,press_HPa,i
      !#######End of water vapour calculations############################
   
      !Additional water from other grids
      call RunoffFromGrid(GridFromFrac)
      
      !Initiate pipe water and runoff including pipes
      runoffPipes=addPipes
      runoff_per_interval=addPipes
      
      
      !=======do calculation for the hour in 60/nsh min intervals==========
      ! loop witithin the hour for no. steps per interval: nsh
      surplusWaterBody=0
      
      do  in=1,nsh
          
          SurPlus_evap=0  !Initiate evaporation surplus
          evap_5min=0
          ev=0
          qe=0
          ev_snow=0
         
          SoilStateOld=soilMoist !Initialize storages
          StateOld=state
          
          !Calculate DRAINAGE for each subsurface excluding water body
          do is=1,nsurf-1
              call drainage(surf(1,is),surf(2,is),surf(3,is),surf(4,is)) !per interval
          enddo
          drain(WaterSurf)=0 !Set drainage from water body to zero

          !Define water gains on different surfaces
          call ReDistributeWater

          
          do is=1,nsurf  
          !===============EVAPORATION and STORAGE for each surface===============
            
            if (snowCalcSwitch(is)==1) then
              call snowCalc(i) 
                
            !Original one, no snow calculations are made!  
            else
               call Evap_SUEWS(surf(1,is))  !qe and ev out
               call soilstore(surf(1,is)) !Soil store updates
               
               if(in==nsh) then   	! if requirement added by HCW 30/07/2014
              		 if (is.ne.WaterSurf) st_per_interval=st_per_interval+state(is)*sfr(is)!State on non-water area (LJ 10/2010)
        	   endif
               
               !Add evaporation to total one
               if (is==BldgSurf.or.is==PavSurf) then
                   ev_per_interval=ev_per_interval+((ev-SurPlus_evap(is))*sfr(is))
                   qe_per_interval=qe_per_interval+((ev-SurPlus_evap(is))*sfr(is))*lv_J_kg
                   
                   evap_5min=evap_5min+((ev-SurPlus_evap(is))*sfr(is))
               else
                   ev_per_interval=ev_per_interval+(ev*sfr(is))
                   qe_per_interval=qe_per_interval+(ev*sfr(is))*lv_J_kg
                   
                   evap_5min=evap_5min+(ev*sfr(is))
               endif

               ch_per_interval=ch_per_interval+(state(is)-stateOld(is))*sfr(is) !Update surface storage change. LJ

               ChangSnow(is)=0
               runoffSnow(is)=0
            endif
          enddo !nsurf
         
          !At this point water has moved between the canopy and soilstorages of each surface.
          !==================================================================================  
          !Now water is allowed to move between the surface stores
          call HorizontalSoilWater
          
          !Calculate soilmoisture state
          Soilmoist_state=0   !Initial value 
          do is=1,nsurf-1 !Summing all together (excluding water body)
              soilmoist_state=soilmoist(is)*sfr(is)+soilmoist_state

              if(soilmoist_state<0)then
                call errorHint(29,'subroutine SUEWS_temporal[soilmoist_state<0],dectime,soilmoist_state,sfr(is)', &
                 dectime,soilmoist_state,int(sfr(is)))
                call errorHint(29,'subroutine SUEWS_temporal[soilmoist_state<0],dectime,soilmoist(is),sfr(is)', &
                dectime,soilmoist(is),int(sfr(is)))

              endif
          enddo
          
          if (SoilMoist_state>soilMoistCap) then!What is this LJ 10/2010
              SoilMoist_state=soilMoistCap
          endif

          !Calculate soil moisture deficit for each time-interval.
          smd=soilMoistCap-soilmoist_state
          smd_nsurf=SoilstoreCap-soilmoist

          !Update changes in soil storages to overall ch_per_interval. Needed as soil stores can change after
          !horizontal water movements. 
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
          enddo

          dectime_nsh = dectime + 1.0*(in-1)/nsh/24

         !Modified hcw 30/07/2014 so that output columns and header match (was 69 cols with extra columns for water)
          if(write5min==1) then           !   Save 5 min results to a file
            dataOut5min(ind5min,1:64)=(/real(id,kind(1D0)),real(in,kind(1D0)),dectime_nsh,pin,ext_wu,ev_per_interval,&
                              stateOut(1:nsurf),smd_nsurfOut(1:(nsurf-1)),drain(1:(nsurf-1)),runoffOut(1:(nsurf-1)),&
                              runoffsoilOut(1:(nsurf-1)),&
                              runoffSnow(1:(nsurf-1)),snowPack(1:nsurf),ChangSnow(1:nsurf),mw_ind(1:nsurf)/)
            ind5min = ind5min+1
          endif
         
      enddo !in=1,nsh (LJ)
      
      !======FINAL STEPS BEFORE WRITING OUT====================================
       
      AdditionalWater=addWaterBody*sfr(WaterSurf)+addPipes+addImpervious*sfr(BldgSurf)+addveg*  &
                     (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassISurf)+sfr(GrassUSurf))

      qeph=qe_per_interval/Interval !Calculate evaporation per interval
      
      !Calculate sensible heat flux as a residual (Modified by LJ in Nov 2012)
      !qh=(qn1+qf+QmRain+QmFreez)-(qeph+qs+Qm) 
      qh=(qn1+qf+QmRain)-(qeph+qs+Qm+QmFreez) 
     
      ext_wuhP=0
      if(GrassISurf>0.or.IrrTrees>0)then                 !Calcuate external water use       
          ext_wuhP=ext_wuh*sfr(GrassISurf)+wuhTrees*IrrTrees       
      end if
      
      if(st_per_interval<0) st_per_interval=0 !Remove negative state

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
    
      lai_wt=0
      do iv=1,NvegSurf
 		lai_wt=lai_wt+lai(id-1,iv)*sfr(iv+2)!Areally weighted LAI
      enddo

      !Calculate snowdepth from SWE
      
      do is = 1,Nsurf
          if (densSnow(is)/=0) then  
          	 SnowDepth(is) = SnowPack(is)*waterDens/densSnow(is) 
          endif
          !Calculate overall snow water equivalent
          swe = swe + SnowDepth(is)*sfr(is)*snowFrac(is)
          MwStore = MwStore + MeltWaterStore(is)*sfr(is)*snowFrac(is)

      enddo
      
     !¤¤¤¤¤¤¤¤¤FILE WRITE SECTION¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
     ! if this changes it will impact LUMPS_RunoffFromGrid
     dataOut1(i,1:62)=(/real(id,kind(1D0)),real(it,kind(1D0)),dectime,avkdn,kup,ldown,lup,tsurf,&  !8
                       qn1,h_mod,e_mod,qs,qf,qh,qeph,Precip_hr,ext_wuhP,ev_per_interval,& !18
                       dr_per_interval,ch_per_interval,st_per_interval,runoffSoil_per_interval,runoff_per_interval,& !23
                       runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,ra,ResistSurf,ustar,& !30
                       l_mod,(smd_nsurfOut(is),is=1,nsurf-1),(stateOut(is),is=1,nsurf),Fcld,soilmoist_state,smd,lai_wt,& !48
                       FlowChange*sfr(WaterSurf),AdditionalWater,ext_wuh,wuhTrees,qn1_SF,qn1_S,Qm,QmFreez,QmRain,& !57
                       swe,MwStore,(SnowRemoval(is),is=1,2),chSnow_per_interval/)

     dataOut2(i,1:30)=(/real(id,kind(1D0)),dectime,kup_ind(1:7),lup_ind(1:7),tsurf_ind(1:7),qn1_ind(1:7)/)
     if (snowUse==1)then!Shiho: This condition is needed when snowUse=0
     dataOut3(i,1:106)=(/real(id,kind(1D0)),real(it,kind(1D0)),dectime,SnowPack(1:7),SnowRemoval(1:2),mwh,mw_ind(1:7),&
                        Qm,Qm_melt(1:7),Qm_rain(1:7),Qm_freezState(1:7),snowFrac(1:6),alb_snow,rainOnSnow(1:7),&
                        qn1_ind_snow(1:7),kup_ind_snow(1:7),freezMelt(1:7),MeltWaterStore(1:7),densSnow(1:7),&
                        snowDepth(1:7),Tsurf_ind_snow(1:7),QmFreez/)
     endif


     !Writes to monthly and daily out
     !First day is not necessarily 1 Jan. Headers are only written with the first line
     !if(id==1.and.reset==1)then

     if(reset==1)then
        iyr=iyr+1
        reset=0
     
        !call OutputHeaders(ProgName,lfnOut,lfnOutC,text,veg_type,ldown_option,1,0) !Removed as unnecessary LJ Nov 2013
        call out_accumulate(366,day,14,0)
        call out_accumulate(12,month,15,0)
        write(15,*)'% season year-',iyr-1
        call out_accumulate(2,season,15,0)
        write(15,*)'% Year:'
        
        call out_accumulate(iyr,yr_tot,61,int(year))
        call accum_zero
        
     !elseif(id==2)then
     !   reset=1
     endif


     call day2month(id,imon,iday,iseas,year,lat)
     call accumulate(id,day)
     call accumulate(imon,month)
     call accumulate(iseas,season)
     
     call accumulate(iyr,yr_tot)
     call accumulate(1,all_tot)

  
     !Writing formats
     !301 format(2i4,f9.4,5f8.2,7f12.4,12f10.4,2f9.1,f10.4,g15.5,27f12.4)!s.o. 7F11.4 ->9F11.4
     !New format Nov 2013 LJ
     301 format(2(i3,1X),f8.4,4(f8.2,1X),(f7.2,1X),7(f8.2,1X),12(f8.3,1X),2(f7.1,1X),& !until 29 col
     (f7.2,1X),(g14.5,1X),12(f8.3,1X),(f10.3,1X),3(f8.3,1X),(f9.2,1X),4(f9.3,1X),5(f8.2,1X),5(f9.3,1X))!s.o. 7F11.4 ->9F11.4
  
    !Calculate new snow fraction used in the next timestep if snowUse==1
    if (SnowFractionChoice==2.and.snowUse==1) then
       do is=1,nsurf-1
          if (snowPack(is)>0.and.mw_ind(is)>0) then
             snowFrac(is)=SnowDepletionCurve(is,snowPack(is),snowD(is),SnowLimPaved,SnowLimBuild)
          elseif (snowPack(is)==0) then
             snowFrac(is)=0
          endif
       enddo
    endif
     
 end do !i=1,10000 (LJ)

 !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 !Save data to files
 call SUEWS_Output


 !Variables passed to future years. Currently: Surface and soil stores and LAI - occurs in nextInitial
 !if(errFileYes==1)close (lfnout)
 !close (lfnoutC)
 close (7)

 do is=1,4
  if (GridFromFrac(is)/=0) then!If runoff from other surfaces exists read data

  !File identifier
  if (is==1) then
      lfnOld=486
  elseif (is==2) then
      lfnOld=487
  elseif (is==3) then
      lfnOld=488
  elseif (is==4) then
      lfnOld=489
  endif
  close(lfnOld)
  endif
 enddo

 call OutputHeaders(ProgName,lfnOutC,text,veg_type,ldown_option,1)
 call out_accumulate(id,day,14,0)
 call out_accumulate(12,month,15,0)
 write(15,*)'% Season year-',iyr
 call out_accumulate(2,season,15,0)
 write(15,*)'% Year:'
 Call out_accumulate(1,all_tot,15,0)
 if(CreateAnnual==1)then
	call out_accumulate(iyr,yr_tot,61,int(year))
 endif
  
 call NextInitial(GridName) ! write grid state information for the next year by grid 
 SnowPack_grid=SnowPack
 
 ! print*,id
 close(15)
 close(14)
 close(140)
 close(118)
!-------------Switched off for SUEWS LJ Sep 2010----------------------------
! Creates 15 and 30 min data from 60 min data, or 30 and 60 min from 15 min data
! or 15 and 60 min from 30 min data
!Works currently only without Suesdata
!if (SuewsStatus==0) then!
!	call DifferentTimeIntervals
!endif
!---------------------------------------------------------------------------


 return

314 	call errorHint(11,trim(filemet),notUsed,notUsed,ios_out)
319 	call errorHint(11,trim(FileNameOld),notUsed,notUsed,ios_out)


 end subroutine SUEWS_temporal
