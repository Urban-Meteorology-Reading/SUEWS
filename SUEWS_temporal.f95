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
  
  implicit none 
  
  !Grid information 
  real(kind(1d0)),DIMENSION(4)::GridFromFrac
  character(len=15)::GridName
  character(len=15),dimension(2,4)::GridFrom
  character(len=15)::year_txt

  real (kind(1d0)),dimension(nsurf)::SnowPack_grid
  
  !Other variables 
  logical:: debug=.false.                            
  integer:: imon,iday,iyr,iseas,reset=1,i,iv,ih,id_in,it_in, SunriseTime,SunsetTime,errFileYes
            
  real(kind(1d0))::lai_wt,dectime_nsh,SnowDepletionCurve
                   
  character(len=100)::FileNameOld,str2

  !Variables related to NARP. Radiation balance components for different surfaces
  real(kind(1D0))::NARP_ALB_is,NARP_EMIS_is,qn1_cum,kup_cum,lup_cum,tsurf_cum,snowFracTot

  
 !Initialize the model (reading nml files, defining filepaths, printing filechoices)
 call OHMinitialize
 !===========================NARP CONFIG=============================================
 if (NARPOutput==1) then
        open(7,file=NARPOut) 
        write(7,110)           
110 	format('%id  dectime   kup_pav   kup_blg kup_everg   kup_dec kup_Irrgr    kup_Gr   kup_wtr  ', &
          ' lup_pav   lup_blg  lup_everg  lup_dec lup_Irrgr   lup_Gr    lup_wtr    Ts_pav    Ts_blg   Ts_everg',&
          '  Ts_dec   Ts_Irrgr     Ts_Gr    Ts_wtr    qn_pav   qn_blg   qn_everg   qn_dec   qn_Irrgr    qn_Gr   qn_wtr') 

 endif
 if(NetRadiationChoice>0)then 
    	call NARP_CONFIG(LAT,LNG,YEAR,TIMEZONE,ALB_SNOW,EMIS_SNOW,TRANS_SITE,Interval,ldown_option)
  		!This for the snow cover fractions
 		!Initiate NARP anyway in order to get surface temperatures 
 		!call NARP_CONFIG(LAT,LNG,YEAR,TIMEZONE,ALB,EMIS,ALB_SNOW,EMIS_SNOW,TRANS_SITE,INTERVAL,ldown_option)
 endif
 finish=.false.

 
!Snow outputfile
if (snowUse==1) then
        open(8,file=SnowOut)
        write(8,111)           
111 	format('%doy  it dectime  SWE_pav SWE_bldg SWE_evergr SWE_dec SWE_irrGr SWE_Gr SWE_water ',&
               ' SnowRem_pav SnowRem_bldg Mw Mw_pav Mw_bldg Mw_evergr  Mw_dec Mw_irrGr    Mw_Gr ',&
               'Mw_water     Qm   Qm_pav Qm_bldg Qm_evergr Qm_dec Qm_irrGr Qm_Gr Qm_water ',&
               'Qa_pav Qa_bldg Qa_evergr Qa_dec Qa_irrGr Qa_Gr Qa_water QmFr_pav ',&
               'QmFr_bldg QmFr_evergr QmFr_dec QmFr_irrGr QmFr_Gr QmFr_water ',&
               'fr_pav fr_bldg fr_evergr fr_dec fr_irrGr fr_Gr alb_snow  rainOnSnow_pav ',&
               'rainOnSnow_bldg rainOnSnow_evergr rainOnSnow_dec rainOnSnow_irrGr rainOnSnow_Gr',&
               'rainOnSnow_water Qn_pavSnow Qs_blgSnow Qs_evergrSnow Qs_decSnow Qs_irrGrSnow Qs_GrSnow ',&
               'Qs_wtrSnow kup_pavSnow kup_blgSnow kup_evergrSnow kup_decSnow kup_irrGrSnow kup_GrSnow ',&
               'kup_wtrSnow ') !Not yet final!
 endif


 !Anthropogenic Heat coefs
 ! this location allows change between years and locations
 

 !=============Get data ready for the qs calculation====================
 write(12,*) '# Met file: ', trim(fileMet), finish

 open(1,file=trim(fileMet),status='old',err=314,iostat=ios_out,position='rewind')
 call skipHeader(1,SkipHeaderMet)

 call SUEWS_Read(1)

 !Allocate output files. Ca be moved elsewhere
 allocate(dataOut1(nlines,62))
 allocate(dataOut2(nlines,30))
 allocate(dataOut3(nlines,106))

 STPH1=0
 if(NetRadiationChoice==0) then !Radiative components are provided as forcing 
    !avkdn=NAN     !Needed for resistances for SUEWS.
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


  if(CBLuse==1) call CBL_initial
    
 !=================================================================================
 !=========INTERVAL LOOP WHERE THE ACTUAL MODEL CALCULATIONS HAPPEN================
 !=================================================================================
 ! assume starting at midnight
 ! first couple of time steps will not be correct for Qs 
 
 
 
  do i=1,nlines !367*48 !want more than it will be old: 4000 24 h days or 1000 96 *15 days!NPeriodsPerYear
   
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


    call ConvertMetData(i)

    !call MetRead(i) !READ FORCING DATA
    
    if(finish) then
      write(*,*)id,it 
      if(finish)exit
    endif
    
    if(CBLuse==1)call CBL

    ! If GISInput Varies
    if(GISInputType==4)then
      id_in=id   ! sg  - gis data can be missing - so now check 
      it_in=it
       call read_gis(finish)
       if(id/=id_in.or.it/=it_in) call ErrorHint(21,FileGIS,float(id),float(it),it_in)
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

        if (snowUse==0) snowFrac=SNOW 
       
        if(ldown_option==1) then !Observed ldown provided as forcing
           ldown=ldown_obs
        else
           ldown=-9              !to be filled in NARP
        endif
        
        if(ldown_option==2) then !observed cloud fraction provided as forcing
              fcld=fcld_obs
        endif
        
        ALB(DecidSurf)=albDec(ID-1) !Change deciduous albedo

        call narp(qn1,kclear,kup,ldown,lup,fcld,tsurf,dectime,avkdn,Temp_C,avRH,&
                  Ea_hPa,Press_hPa,ldown_option,AlbedoChoice,qn1_obs,&
                  netRadiationChoice,alb_snow,qn1_SF,qn1_S)
    else
        snowFrac = snow
        qn1=qn1_obs
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
               
               if (is.ne.WaterSurf) st_per_interval=st_per_interval+state(is)*sfr(is)!State on non-water area (LJ 10/2010)
        
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
                 dectime,soilmoist_state,sfr(is))
                call errorHint(29,'subroutine SUEWS_temporal[soilmoist_state<0],dectime,soilmoist(is),sfr(is)', &
                dectime,soilmoist(is),sfr(is))

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

          if(write5min==1) then           !   write 5 min results       

            write(16,36)id,in,dectime_nsh,pin,ext_wu,ev_per_interval,(stateOut(is),is=1,nsurf),& !1-13
                  (smd_nsurfOut(is),is=1,nsurf-1),(drain(is),is=1,nsurf-1),(runoffOut(is),is=1,nsurf-1),& !14-19,20-25,26-31
                  (runoffsoilOut(is),is=1,nsurf-1),(runoffSnow(is),is=1,nsurf-1),(snowPack(is),is=1,nsurf),& !32-37,38-43,44-50
                  (ChangSnow(is),is=1,nsurf),(mw_ind(is),is=1,nsurf)
                   
36        format(2i3,f9.4,61f9.3)     !format(i3,i3,3f12.4,6f10.4, f14.3,25f10.4)   
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
      
     !中中中中了ILE WRITE SECTION中中中中中中中中中中中中中�
     ! if this changes it will impact LUMPS_RunoffFromGrid
     

    write(lfnoutC,301)id,it,&
          dectime,avkdn,kup,ldown,lup,tsurf,&  !8
          qn1,h_mod,e_mod,qs,qf,qh,qeph,Precip_hr,ext_wuhP,ev_per_interval,& !18
          dr_per_interval,ch_per_interval,st_per_interval,runoffSoil_per_interval,runoff_per_interval,& !23
          runoffPipes,runoffAGimpervious,runoffAGveg,runoffWaterBody,ra,ResistSurf,ustar,& !30
          l_mod,(smd_nsurfOut(is),is=1,nsurf-1),(stateOut(is),is=1,nsurf),Fcld,soilmoist_state,smd,lai_wt,& !48
          FlowChange*sfr(WaterSurf),AdditionalWater,ext_wuh,wuhTrees,qn1_SF,qn1_S,Qm,QmFreez,QmRain,& !57
          swe,MwStore,(SnowRemoval(is),is=1,2),chSnow_per_interval !62

          
     !Radiation output
     if (NARPOutput==1) then
     	write(7,117)id,dectime,(kup_ind(is),is=1,nsurf),(lup_ind(is),is=1,nsurf),&
          	(tsurf_ind(is),is=1,nsurf),(qn1_ind(is),is=1,nsurf)       
        117 format(i3,f9.4,28f10.3)
 	 endif 

     !Snow output
     if (snowUse==1) then 
        write(8,118)id,it,dectime,(SnowPack(is),is=1,nsurf),(SnowRemoval(is),is=1,2),mwh,& !1-13
             (mw_ind(is),is=1,nsurf),Qm,(Qm_melt(is),is=1,nsurf),& !14-28
             (Qm_rain(is),is=1,nsurf),(Qm_freezState(is),is=1,nsurf),& !29-42
             (snowFrac(is),is=1,nsurf-1),alb_snow,(rainOnSnow(is),is=1,nsurf),&!43-56
             (qn1_ind_snow(is),is=1,nsurf),(kup_ind_snow(is),is=1,nsurf),& !57-70
             (freezMelt(is),is=1,nsurf),(MeltWaterStore(is),is=1,nsurf),& !71-84
             (densSnow(is),is=1,nsurf),(snowDepth(is),is=1,nsurf),& !85 - 91, 92 - 98
             (Tsurf_ind_snow(is),is=1,nsurf),QmFreez !99 -105, 106
             
        118 format(2(i3,1X),(f8.4,1X),17(f8.3,1X),22(f7.2,1X),7(f7.2,1X),7(f7.3,1X),14(f7.2,1X),&
        14(f7.3,1X),21(f7.2,1X),(f14.3,1X))     
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

 !Variables passed to future years. Currently: Surface and soil stores and LAI - occurs in nextInitial 

 if(errFileYes==1)close (lfnout)
 close (lfnoutC)
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

 call OutputHeaders(ProgName,lfnOut,lfnOutC,text,veg_type,ldown_option,1,0)
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

 ! stop 'finished'
 return

314 	call errorHint(11,trim(filemet),notUsed,notUsed,ios_out)
319 	call errorHint(11,trim(FileNameOld),notUsed,notUsed,ios_out)


end subroutine SUEWS_temporal
