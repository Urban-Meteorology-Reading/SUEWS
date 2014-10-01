!This subroutine makes snow related calculations on hourly time step. Needed for the 
!available energy in LUMPS and SUEWS. Made by LJ in Dec 2012
!mod by lj in May 2013 - calculation of the energy balance for the snowpack was modified
!                        to use qn1_ind_snow(surf)


subroutine MeltHeat(i)

  use allocateArray
  use snowMod
  use data_in
  use snowMod
  use time
  use moist
  use sues_data
  use gis_data
  
  implicit none

  REAL(KIND(1d0)):: Watfreeze,& !State of snow on ith surface, total energy exchange (not exported!)
                    ci=2090,cw=4190       !Specific heat capacity of water
 
  INTEGER :: i 
 
  !Initialize snow variables

  snowCalcSwitch=0     !Initialize switches
  StateFraction=0
  snowCoverForms=0
      
  Qm_melt=0  
  Qm_freezState=0 
  Qm_rain=0   
  FreezMelt=0 
  FreezState=0
  FreezStateVol=0
  rainOnSnow = 0
  
  !Initialize mechanical snow removal
  SnowRemoval(1)=0
  SnowRemoval(2)=0
  
  !Go each surface through if surface type existing
  do is=1,nsurf
    if (sfr(is)/=0) then
      
        !If snowpack existing, these calculations are made
        if (SnowPack(is)>0) then

             SnowDepth(is) = (SnowPack(is)/1000)*waterDens/densSnow(is) !in m

             !Calculate meltwater related water flows from degree-day method
             
             !Melt equations
             if (Temp_C>=0) then
                 if (qn1_ind_snow(is)<0) then
                 	mw_ind(is) = TempMeltFact*Temp_C  ! in mm h-1
                 else
                    mw_ind(is) = RadMeltFact*(qn1_ind_snow(is))
                endif

             else  !Freezing equations
                 AdjMeltFact=1
                 mw_ind(is) = TempMeltFact*Temp_C*AdjMeltFact ! in mm h-1
                
             endif 

             if (mw_ind(is)>SnowPack(is)) then !Limited by the previous timestep snowpack
                mw_ind(is) = SnowPack(is)
             endif

             
             !Hourly heat consumed to snowmelt/refreezing. Converted from mm h-1 to mm s-1 and to m s-1
             Qm_melt(is) = waterDens*((mw_ind(is)/3600)/1000)*(lvS_J_kg-lv_J_kg)

             !If melt is negative this means freezing water in the snowpack
             if (mw_ind(is)<0) then
            
                FreezMelt(is) = -mw_ind(is)
                mw_ind(is)=0
                
                !Freezing water cannot exceed the meltwaterstore
                if (FreezMelt(is)>Meltwaterstore(is)) FreezMelt(is) = Meltwaterstore(is)
             
                !Recalculate melt related energy
                Qm_melt(is) = waterDens*((-FreezMelt(is)/3600)/1000)*(lvS_J_kg-lv_J_kg)

             endif
            
    
            !------If air temperature is above zero, precipitation causes advective heat to the
            !------snowpack. Eq (23) in Sun et al., 1999
            if (Temp_C>=PrecipLimit.and.Precip_hr>0) then
               Qm_rain(is) = waterDens*cw*(Temp_C-PrecipLimit)*(Precip_hr*0.001/3600)  !in W m-2
               
               if (Qm_rain(is)<0) then !Can only be positive
                  Qm_rain(is) = 0  
               else
                  rainOnSnow(is) = Precip_hr
               endif
            endif
   
        endif !End if snowpack
         
        !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
        !Freeze surface state if cold enough and there is water, this will freeze.
        if (Tsurf_ind(is)<0.and.state(is)>0) then
    
           snowCalcSwitch(is)=1 !If water on ground this forms ice and snow calculations are made
    
           !no water surface
           if (is/=WaterSurf) then
             
              if (snowpack(is)==0.or.snowfrac(is)==0) then !Snowpack forms. need to be fixed at some point LJ
                  FreezState(is) = state(is)
                  FreezStateVol(is) = FreezState(is)
              else                                !There is snow already on ground
                  FreezState(is) = state(is)
                  FreezStateVol(is) = state(is)*(1-snowFrac(is))/snowFrac(is)
              endif
              
              Qm_freezState(is) = -waterDens*(FreezState(is)/3600/1000)*(lvS_J_kg-lv_J_kg)
               
           else       !Water surface 
             !Calculate average value how much water can freeze above the water areas
             Watfreeze = 100*(0-Temp_C)/(waterDens*(lvS_J_kg-lv_J_kg)/3600)
            
             if (Watfreeze<=state(is)) then !Not all state freezes
                FreezState(is) = Watfreeze  
                Qm_freezState(is) = -waterDens*(Watfreeze/3600/1000)*(lvS_J_kg-lv_J_kg)
             else
                FreezState(is) = state(is)  !All state freezes
                Qm_freezState(is) = -waterDens*(state(is)/3600/1000)*(lvS_J_kg-lv_J_kg)
             endif
           endif
    
         endif
         
         !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
         
         !Define if any snowmelt calculations are made
         !snowpack existing,freezing occuring on ground or from precip
         if (is/=WaterSurf) then  
            if (snowPack(is)>0.or.(Precip_hr>0.and.Tsurf_ind(is)<0)) then
               snowCalcSwitch(is)=1
            endif
         else       !Water surface separately
            if (snowPack(WaterSurf)>0.or.FreezState(WaterSurf)>0) then
               snowCalcSwitch(WaterSurf)=1
            endif
         endif

         !Update snow density of each surface
         if (Precip_hr>0.and.Tsurf_ind(is)<0) then
         	densSnow(is) = densSnow(is)*snowPack(is)/(snowPack(is)+Precip_hr)+densSnowMin*Precip_hr/(snowPack(is)+Precip_hr)
  		 endif

         !Weighted variables for the whole area
         mwh = mwh + mw_ind(is)*sfr(is)*snowFrac(is)  !Snowmelt
         fwh = fwh + FreezMelt(is)*sfr(is)*snowFrac(is) !Freezing water 
         
         Qm = Qm + Qm_melt(is)*sfr(is)*snowFrac(is)   !Energy consumed to the melt/freezing.
    
         QmRain = QmRain + Qm_rain(is)*sfr(is)*snowFrac(is)
         !QmFreez = QmFreez + deltaQi(is)*sfr(is)*snowFrac(is)+Qm_freezState(is)*sfr(is)*(1-snowFrac(is))
         QmFreez = QmFreez + deltaQi(is)*sfr(is)*snowFrac(is)+Qm_freezState(is)*sfr(is)*(1-snowFrac(is))
     endif
     
    enddo !End surface type

    !Update snow albedo and density
    if (Precip_hr>0.and.sum(snowPack)>0.and.Temp_C<0) then
       CumSnowfall=CumSnowfall + Precip_hr
       if (CumSnowfall>PrecipLimitAlb)  alb_snow=albSnowMax
    else
        CumSnowfall=0
    endif

   
 end subroutine MeltHeat



!===============================================================================================
!===============================================================================================
subroutine SnowCalc(i)
  !Calculation of snow and water balance on 5 min timestep. Treats snowfree and snow covered 
  !areas separately. Weighting is taken into account in the overall values.
  !Last modified LJ in 24 May 2013
  !-------------------------------- 
  
  use allocateArray
  use defaultNotUsed
  use data_in
  use snowMod
  use time
  use moist
  use sues_data
  use gis_data
  use mod_k
  use mod_z
  
  implicit none
             
  REAL(KIND(1d0)):: Evap_SUEWS_Snow,&
                    MeltExcess,&  !Excess melt water that needs to leave snowpack
                    snowTotInit,&
                    EvPart,&
                    runoffTest,&
                    FWC           !Water holding capacity of snow in mm
                  
                    
  REAL(KIND(1d0))::SnowDepletionCurve !Function
  
  integer::i
  
  runoffSnow(is)=0 !Initialize for runoff caused by snowmelting
  runoff(is)=0
  runoffSoil(is)=0
  chang(is)=0
  changSnow(is)=0
  runoffTest=0
  SnowToSurf(is)=0
  
  !Initial snowpack + meltwater in it  
  snowTotInit=snowPack(is)+MeltWaterStore(is) 
  
 !Calculate water holding capacity (Jin et al. 1999)
  if (densSnow(is)>=200) then
     WaterHoldCapFrac=CRWmin
  else
     WaterHoldCapFrac=CRWmin+(CRWmax-CRWmin)*(200-densSnow(is))/200  
  endif

  
  !===========Initialization of evaporations====================
  EvPart=0
  ev=0

  !Calculate evaporation from snowpack and snow free surfaces (in mm)
  if (snowFrac(is)<1) call Evap_SUEWS(surf(1,is)) 
  
  if (snowFrac(is)>0) then
      ev_snow(is) = Evap_SUEWS_Snow(qn1_S,qf,qs,Qm_Melt(is),Qm_rain(is),lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RA,&
                                      psyc_hPa,ity,tstep,avcp,sIce_hPa,deltaQi(is))                        
  endif

  !If not enough water for evaporation in impervious surfaces,
  !evaporation is taken from pervious surfaces
  if (is>2) then
     if  ((VegFraction+sfr(WaterSurf))/=0) then
        EvPart=(SurPlus_evap(PavSurf)+SurPlus_evap(BldgSurf))*(sfr(is)/(VegFraction+sfr(WaterSurf)))
     endif
  endif
 
  
  !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
  !First if the surface is fully covered with snow or the snowpack forms 
  !(equally distributed) on the current hour
  
  if (is==WaterSurf.and.sfr(WaterSurf)>0) GO TO 606  !Water surface is treated separately

    
   
  !First is whole surface is covered with snow or if snowpack forms at this timestep
  if (snowFrac(is)==1.or.(snowPack(is)==0.and.pin>0.and.Temp_C<0.and.Tsurf_ind(is)<0)&
    .or.(freezState(is)>0)) then
        
     !------Snow pack exists-----------------------------------------
     if (SnowPack(is)>0) then
        ev_snow(is)=ev_snow(is)+EvPart !Evaporataion surplus

        
       !(Snowfall per interval+freezing of water) - (meltwater+evaporation from snowpack)
        changSnow(is)=(pin+freezMelt(is)/nsh)-(mw_ind(is)/nsh+ev_snow(is)) !Calculate change in the snowpack
        
        if (freezState(is)>0) then !snowfraction is not 1 but the snow free surface water freezes
           FreezStateVol(is)=FreezState(is)
           changSnow(is) = changSnow(is)+freezStateVol(is)/nsh
           ev=0

           if (in==1) iceFrac(is)=1-snowFrac(is)
           snowFrac(is)=1 
           
        endif
        
        !Let runoff happen if rain on snow event
        if (rainOnSnow(is)>0) then
           changSnow(is)=changSnow(is)-pin
           MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)/nsh
        endif
        
     !------Snowpack forms at this timestep from precipitation or freezing state---   
     else 
        ev=ev+EvPart
        changSnow(is)=(pin+freezStateVol(is)/nsh)-ev 
        snowFrac(is)=1
        iceFrac(is)=0
        densSnow(is)=densSnowMin
     endif
     
     !--------------------------------------------------------------------------------------
     SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack
     state(is)=state(is)-freezState(is)/nsh   !Update state if water has frozen
     

     !If snowpack exists after the state calculations
     if (snowPack(is)>0) then  
    
        !Add melted water to meltstore and freeze water according to freezMelt(is)
        MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is)/nsh - freezMelt(is)/nsh
      
        !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
        FWC = WaterHoldCapFrac*snowPack(is)
     
        !If FWC is exceeded, add meltwater to runoff 
        !if (MeltWaterStore(is)>=FWC.and.Temp_C>PrecipLimit) then
        if (MeltWaterStore(is)>=FWC) then
 			runoffSnow(is)=runoffSnow(is)+(MeltWaterStore(is)-FWC)
            MeltWaterStore(is) = FWC
        endif
        
        !At the end of it (hour) calculate possible snow removal
        !if (it>=6.and.it<=18.and.in==12.and.is<3) call snowRem 
        if (snowProf(it)==1.and.in==12.and.is<3) call snowRem
   
     !If snowPack is negative, it melts at this timestep
     elseif (SnowPack(is)<0) then

        !If freezing meltwater inside this hour, remove it from the MeltWaterStore
        MeltWaterStore(is)=MeltWaterStore(is)-freezMelt(is)/nsh+mw_ind(is)/nsh+SnowPack(is) 
        SnowPack(is)=0
        snowFrac(is)=0
  
        if (MeltWaterStore(is)<0) then !Not enough water in the meltwater store, 
                            
            if (ev_snow(is)/=0) then              
               ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased. ZZ Possible problem
               if (ev_snow(is)<0) ev_snow(is)=0
               changSnow(is)=changSnow(is)+MeltWaterStore(is)
               MeltWaterStore(is)=0
            else
              ev=ev+MeltWaterStore(is)
              if (ev<0) ev=0
              changSnow(is)=changSnow(is)+MeltWaterStore(is)
              MeltWaterStore(is)=0
            endif
            
        else
            chang(is)=MeltWaterStore(is)  !Meltwater goes to surface state as no snow exists anymore
            state(is)=state(is)+chang(is) 
            MeltWaterStore(is)=0
            
        endif
      
     endif !snowpack negative or positive

     !------Change state of snow and surface
     chSnow_per_interval=chSnow_per_interval+((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)
     ch_per_interval=ch_per_interval+(state(is)-stateOld(is))*sfr(is)

     !------Evaporation
     if (is==BldgSurf.or.is==PavSurf) then 
        ev_per_interval=ev_per_interval+(ev-SurPlus_evap(is))*sfr(is)+ev_snow(is)*sfr(is)
        qe_per_interval=qe_per_interval+ev_snow(is)*lvS_J_kg*sfr(is)+(ev-SurPlus_evap(is))*lv_J_kg*sfr(is)
       evap_5min=evap_5min+(ev-SurPlus_evap(is))*sfr(is)+ev_snow(is)*sfr(is)
     else  
        ev_per_interval=ev_per_interval+ev*sfr(is)+ev_snow(is)*sfr(is)
        qe_per_interval=qe_per_interval+ev_snow(is)*lvS_J_kg*sfr(is)+ev*lv_J_kg*sfr(is)
        evap_5min=evap_5min+ev*sfr(is)+ev_snow(is)*sfr(is)
     endif

     !==========RUNOFF===================
     runoffPipes=runoffPipes+runoffSnow(is)*sfr(is) !Runoff to pipes 
 
    !If pipe capacity is full, surface runoff occurs
    if (runoffPipes>PipeCapacity) then
         if (is==PavSurf.or.is==BldgSurf) then
  
             if (sfr(WaterSurf)>0.0000001) then
                 runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)*(1-RunoffToWater) !if pipes are full, water will stay above ground
                 surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
             else
                 runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)
             endif
             runoffPipes=PipeCapacity                                         !and pipes are in their max. capacity
                                          
         elseif (is>=ConifSurf.and.is<=GrassUSurf) then
         
             if (sfr(WaterSurf)>0.0000001) then
                 runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)*(1-RunoffToWater)       !if pipes are full, water will stay above ground
                 surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
             else
                  runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity) 
             endif
             runoffPipes=PipeCapacity                                 !and pipes are in their max. capacity
             
         endif
     endif
     
     runoff_per_interval=runoff_per_interval+runoffSnow(is)*sfr(is)
  

  !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
  !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
   
  !If both snow and snow free surface exists 
  else
   
    !Snowpack water balance for the whole surface area. In reality snow depth = snowPack/snowFrac(is)  
    changSnow(is)=(pin+freezMelt(is)/nsh)-(mw_ind(is)/nsh+ev_snow(is))

    if (rainOnSnow(is)>0) then  ! Precipitation is routed to meltwater store
        changSnow(is)=changSnow(is)-pin 
        MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)/nsh
    endif
  
    SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack
    
    !If snowpack exists at this point
    if (snowPack(is)>0) then  

      !Add melted water to meltstore and freeze water according to freezMelt(is)
      MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is)/nsh - freezMelt(is)/nsh
      
      !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
      FWC = WaterHoldCapFrac*snowPack(is)!(snowPack(is)+MeltWaterStore(is))
     
      !If FWC is exceeded, add meltwater to runoff
      !if (MeltWaterStore(is)>=FWC.and.Temp_C>PrecipLimit) then
      if (MeltWaterStore(is)>=FWC) then
         MeltExcess=0
         MeltExcess = MeltWaterStore(is)-FWC  
         
         if (snowFrac(is)>0.9) then
         	runoffSnow(is)=runoffSnow(is) + MeltExcess*snowFrac(is)
         	!stateExcess = stateExcess + MeltExcess*snowFrac(is)/(1-snowFrac(is))*(1-snowFrac(is))
         	SnowToSurf(is) =  SnowToSurf(is) + MeltExcess*snowFrac(is)
         else
            SnowToSurf(is) =  SnowToSurf(is) + MeltExcess*snowFrac(is)/(1-snowFrac(is))
         endif
         
         MeltWaterStore(is) = FWC
      endif
      
      !At the end of it (hour) calculate possible snowremoval
      !if (it>=6.and.it<=18.and.in==12.and.is<3) call snowRem 
      if (snowProf(it)==1.and.in==12.and.is<3) call snowRem 
        
    !If snowPack is negative
    elseif (SnowPack(is)<0) then

      !If freezing meltwater inside this hour, remove it from the MeltWaterStore
      MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is)/nsh - freezMelt(is)/nsh
      MeltWaterStore(is) = MeltWaterStore(is) + SnowPack(is) 
      SnowPack(is)=0
      snowFrac(is)=0
      
      if (MeltWaterStore(is)<0) then !Not enough water in the meltwater store, 
                                        
          ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased
          if (ev_snow(is)<0)  ev_snow(is)=0
          changSnow(is)=changSnow(is)+MeltWaterStore(is)
          MeltWaterStore(is)=0
          
      else
          chang(is)=MeltWaterStore(is)*snowFrac(is)/(1-snowFrac(is)) !Meltwater goes to surface state 
          SnowToSurf(is)=SnowToSurf(is)+chang(is) 
          MeltWaterStore(is)=0
      endif
      
    endif !snowpack negative or positive

   
    !----------Snow free surface-----------
 
    !pin=pin  !Add water from other surfaces (in volumes)
    
    if (is==PavSurf.or.is==BldgSurf) then  !Impervious surfaces (paved, buildings)
  
        !Surface store update       
        if (pin>10) then
           runoff(is)=runoff(is)+(pin+SnowToSurf(is)+AddWater(is)-10)
           chang(is)=10-(drain(is)+ev+freezState(is)/nsh) 
        else 
           chang(is)=pin+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is)/nsh)   !Add precip and water from other surfaces
        endif                                                !remove drainage, evap and freezing of state
     
        state(is)=state(is)+chang(is) !Change in state (for whole surface area areasfr(is))
     
        !Add water from neighbouring grids 
        if (is==PavSurf) state(is)=state(is)+(addImpervious/nsh)
        
        runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is) !Drainage (not flowing to other surfaces) goes to runoff
        
        if(state(is)<0.0) then  !Surface state cannot be negative
           SurPlus_Evap(is)=abs(state(is)) !take evaporation from other surfaces in mm 
           state(is)=0.0
        endif
       

    elseif(is>=3) then ! Pervious surfaces (conif, decid, grass unirr, grass irr)
     
       ev=ev+EvPart
     
       !Additional water input from other grids
       if (VegFraction/=0) then
        	pin=pin+addVeg/nsh*(sfr(is)/VegFraction)
       endif
       
        
       !Change in water stores 
       if (pin>10) then !if 5min precipitation is larger than 10 mm
           runoff(is)=runoff(is)+(pin+SnowToSurf(is)+AddWater(is)-10)
           chang(is)=10-(drain(is)+ev+freezState(is)/nsh)
       else
           chang(is)=pin+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is)/nsh)    
       endif
     
       state(is)=state(is)+chang(is)
        
       !Add water in soil store only if ground is not frozen
       !if (Temp_C>0) then
       soilmoist(is)=soilmoist(is)+Drain(is)*AddWaterRunoff(is)*(1-snowFrac(is))
       !else
       !	 runoff(is)=runoff(is)+Drain(is)*AddWaterRunoff(is)
       !endif		
     
       !If state of the surface is negative, remove water from soilstore
       if(state(is)<0.0) then  
   
          if ((soilmoist(is)+state(is))>=0.and.Temp_C>0) then !If water in soilstore, water is removed
            
              soilmoist(is)=soilmoist(is)+state(is)*(1-snowFrac(is))
              state(is)=0.0

          else !If not water in the soilstore evaporation does not occur
              chang(is)=chang(is)+state(is)
              ev=ev+state(is)
              state(is)=0.0
          endif
       endif !state is negative
       
       !If soilstorage is full at this point, excess will go to surface runoff
       if (soilmoist(is)>soilstoreCap(is)) then  
          runoffTest=runoffTest+(soilmoist(is)-soilstoreCap(is))
          soilmoist(is)=soilstoreCap(is)
       elseif (soilmoist(is)<0) then
          soilmoist(is)=0
       endif
    
     !State of non-water area in mm 
      st_per_interval=st_per_interval+state(is)*sfr(is)*(1-snowFrac(is))
   
   endif !Surface type

   !Calculate change in snowpack and state
   ch_per_interval=ch_per_interval+(state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is))
   chSnow_per_interval=chSnow_per_interval+((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)*snowFrac(is)
   
   !Add evaporation to total 
   if (is==BldgSurf.or.is==PavSurf) then
       ev_per_interval=ev_per_interval+(ev-SurPlus_evap(is))*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
       qe_per_interval=qe_per_interval+ev_snow(is)*lvS_J_kg*sfr(is)*snowFrac(is)&
                       +(ev-SurPlus_evap(is))*lv_J_kg*sfr(is)*(1-snowFrac(is))
       evap_5min=evap_5min+(ev-SurPlus_evap(is))*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
   else
       ev_per_interval=ev_per_interval+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
       qe_per_interval=qe_per_interval+ev_snow(is)*lvS_J_kg*sfr(is)*snowFrac(is)+ev*lv_J_kg*sfr(is)*(1-snowFrac(is))
       evap_5min=evap_5min+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
   endif
  
  !========RUNOFF=======================
   
  !!Add runoff to pipes 
  runoffPipes=runoffPipes+runoffSnow(is)*sfr(is)*snowFrac(is)+runoff(is)*sfr(is)*(1-snowFrac(is))+runoffTest*sfr(is)
  
    !If pipe capacity is full, surface runoff occurs
    if (runoffPipes>PipeCapacity) then
         if (is==PavSurf.or.is==BldgSurf) then
  
             if (sfr(WaterSurf)>0.0000001) then
                 runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)*(1-RunoffToWater) !if pipes are full, water will stay above ground
                 surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
             else
                 runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)
             endif
             runoffPipes=PipeCapacity                                         !and pipes are in their max. capacity
                                          
         elseif (is>=ConifSurf.and.is<=GrassUSurf) then
         
             if (sfr(WaterSurf)>0.0000001) then
                 runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)*(1-RunoffToWater)       !if pipes are full, water will stay above ground
                 surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
             else
                  runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity) 
             endif
             runoffPipes=PipeCapacity                                 !and pipes are in their max. capacity
             
         endif
     endif
     
     runoff_per_interval=runoff_per_interval+runoffSnow(is)*sfr(is)*snowFrac(is)+runoff(is)*sfr(is)*(1-snowFrac(is))&
                         +runoffTest*sfr(is)
     !runoff_per_interval=runoff_per_interval+runoffSnow(is)*sfr(is)+runoff(is)*sfr(is)+runoffTest*sfr(is)                   


  endif !Ground both snow covered and snow free    


  return

  !いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
  !WATERBODY is treated separately as state always below ice if ice existing
  !Calculate change in snowpack
606 changSnow(WaterSurf)=(pin+freezMelt(WaterSurf)/nsh+freezState(WaterSurf)/nsh)-&
                         (mw_ind(WaterSurf)/nsh+ev_snow(WaterSurf))
                         
    SnowPack(WaterSurf)=SnowPack(WaterSurf)+changSnow(WaterSurf) !Update snowpack
    state(WaterSurf)=state(WaterSurf)+FlowChange/nsh-freezState(WaterSurf)/nsh  !Update state below ice
     
    !If snowpack exists
    if (snowPack(WaterSurf)>0) then  
    
       !Add melted water to meltstore and freeze water according to freezMelt(is)
       MeltWaterStore(WaterSurf)=MeltWaterStore(WaterSurf)+mw_ind(WaterSurf)/nsh-freezMelt(WaterSurf)/nsh
      
       !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
       FWC = WaterHoldCapFrac*SnowPack(WaterSurf)
     
       !If FWC is exceeded, add meltwater to state
       if (MeltWaterStore(WaterSurf)>=FWC) then
          state(WaterSurf)=state(WaterSurf)+(MeltWaterStore(WaterSurf)-FWC)
          MeltWaterStore(WaterSurf) = FWC
       endif

     !If snowPack is negative, it melts at this timestep
     elseif (SnowPack(is)<0) then
     
        !Add water to the meltwater store
        !If freezing meltwater inside this hour, remove it from the MeltWaterStore
        MeltWaterStore(WaterSurf) = MeltWaterStore(WaterSurf)-freezMelt(WaterSurf)/nsh &
                                    + mw_ind(WaterSurf)/nsh
       
        state(WaterSurf)=state(WaterSurf)+MeltWaterStore(WaterSurf)+SnowPack(WaterSurf) !Add meltwater to state
        SnowPack(WaterSurf)=0
        if (state(WaterSurf)<0) ev_snow(WaterSurf)=ev_snow(WaterSurf)+state(WaterSurf)
 
     endif !snowpack negative or positive

     !Check water state separately
     if (state(WaterSurf)>Surf(5,WaterSurf)) then
         runoff(WaterSurf)=runoff(WaterSurf)+(state(WaterSurf)-Surf(5,WaterSurf))
         state(WaterSurf)=Surf(5,WaterSurf)
         runoffWaterBody=runoffWaterBody+runoff(WaterSurf)*sfr(WaterSurf)
     else  
         state(WaterSurf)=state(WaterSurf)+surplusWaterBody
        
         if (state(WaterSurf)>Surf(5,WaterSurf)) then
             runoffWaterBody=runoffWaterBody+(state(WaterSurf)-Surf(5,WaterSurf))*sfr(WaterSurf)
             state(WaterSurf)=Surf(5,WaterSurf)
         endif
     endif
     

     !Change state of snow and surface
     chSnow_per_interval=chSnow_per_interval+((snowPack(WaterSurf)+MeltWaterstore(WaterSurf))-snowTotInit)*sfr(WaterSurf)
     ch_per_interval=ch_per_interval+(state(WaterSurf)-stateOld(WaterSurf))*sfr(WaterSurf)
     
     !Evaporation
     ev_per_interval=ev_per_interval+ev*sfr(WaterSurf)+ev_snow(WaterSurf)*sfr(WaterSurf)
     qe_per_interval=qe_per_interval+ev_snow(WaterSurf)*lvS_J_kg*sfr(WaterSurf)+ev*lv_J_kg*sfr(WaterSurf)
     evap_5min=evap_5min+ev*sfr(is)+ev_snow(WaterSurf)*sfr(WaterSurf)
 
     runoff_per_interval=runoff_per_interval+(runoff(is)*sfr(is)) !The total runoff from the area
     
     if (snowPack(WaterSurf)>0) then     !Fraction only 1 or 0
         snowFrac(WaterSurf)=1
     else
         snowFrac(WaterSurf)=0
     endif

  end subroutine SnowCalc
  
  
  !==========================================================================
  !==========================================================================
  
 
  !==========================================================================
  FUNCTION Evap_SUEWS_Snow(qn1,qf,qs,Qm,QP,lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RAsnow,psyc_hPa,ity,&
                           tstep,avcp,sIce_hPa,Qfreez) RESULT(ev_snow)
  !Calculates evaporation from snow surface (ev_snow). 
  !INPUT:  
  

   implicit none
   
   real (Kind(1d0))::sae_snow,lvS_J_kg,e_snow,qe_snow,avdens,Temp_C,RA,Qfreez,&
                     qn1,qf,qs,Qm,ev_snow,vdrcIce,esIce_hPa,EaIce_hPa,avRh,Press_hPa,&
                     psyc_hPa,tstep,tlv_sub,avcp,sIce_hPa,QP,RAsnow
                     
   real (Kind(1d0)):: sat_vap_pressIce !Function
   integer:: ity,from=1
   
   !s * (Available energy)
   !sae_snow=sIce_hPa*(qn1+Qp-Qm-Qfreez)
   sae_snow=sIce_hPa*(Qp-Qm)


   !Saturation vapor pressure over ice
   esIce_hPa= sat_vap_pressIce(Temp_C,Press_hPa,from)
   EaIce_hPa=avRh/100*esIce_hPa  
   
   vdrcIce=(esIce_hPa-eaIce_hpa)*avdens*avcp
   tlv_sub=lvS_J_kg/tstep                  !Latent heat for sublimation!
  
   e_snow=sae_snow+vdrcIce/RAsnow

   
   qe_snow=e_snow/(sIce_hPa+psyc_hPa)!Latent heat (W/m^2)

   ev_snow=qe_snow/tlv_sub                !Evaporation (in mm)
  
   return
   
  END FUNCTION Evap_SUEWS_Snow
  
  !==========================================================================
  
  subroutine snowRem
  ! Calculates mechanical removal of snow from roofs ans roads 
   use allocateArray 
   use snowMod
   use sues_data
   use time
  
  implicit none

 
  if (is==PavSurf) then
     if (SnowPack(PavSurf)>SnowLimPaved) then
         SnowRemoval(PavSurf) = (SnowPack(PavSurf)-SnowLimPaved)*sfr(PavSurf)*snowfrac(PavSurf)
         SnowPack(PavSurf)=SnowLimPaved
         !snowPack(PavSurf)=snowPack(PavSurf)/snowFrac(PavSurf)
     endif
  endif
  if (is==BldgSurf)then
     if (SnowPack(BldgSurf)>SnowLimBuild) then
         SnowRemoval(2) = (SnowPack(BldgSurf)-SnowLimBuild)*sfr(BldgSurf)*snowfrac(BldgSurf)
         SnowPack(BldgSurf)=SnowLimBuild
         !snowPack(BldgSurf)=snowPack(BldgSurf)/snowFrac(BldgSurf)
     endif  
  endif
 
 
  end subroutine snowRem
 
  !-----------------------------------------------
  !-----------------------------------------------
  FUNCTION SnowDepletionCurve(is,swe,sweD,SnowLimPaved,SnowLimBuild) RESULT(asc)
  !This function calculates surface coverage of snow according to the
  !depletion curves in Valeo and Ho (2004). Done once per hour.
  
   use allocateArray
  
   IMPLICIT  NONE
  
   INTEGER::is
   REAL (KIND(1d0))::asc,sweD,swe,SnowLimPaved,SnowLimBuild
  
   !Impervious surface 
   if (is==PavSurf) then
       
     if (swe<=sweD) then      !Snow water equivalent below threshold
        asc=((swe/sweD))**2
     else
        asc=1
     endif
     
   !Built surface  
   elseif (is==BldgSurf) then
   
     if (swe<=sweD) then
        if ((swe/sweD)<0.9) then
          asc=(swe/sweD)*0.5
        else
          asc=(swe/sweD)**8
        endif  
     else
        asc=1
     endif
   elseif (is==WaterSurf) then
      if (swe>0) asc=1 
   
   !Vegetion surfaces    
   else
     if (swe<=sweD) then
        
        asc=1-((1/3.1416)*acos(2*(swe/sweD)-1))**1.7
     else
        asc=1
     endif
     
   endif

   !asc=real(int(10000.*asc))/10000  !4 decimal precision 

   RETURN
  END FUNCTION SnowDepletionCurve