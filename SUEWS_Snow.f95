!This subroutine makes snow related calculations at the model time step. Needed for the 
!available energy in LUMPS and SUEWS. Made by LJ in Dec 2012
!SUBROUTINES:
!  MeltHeat - Calculation of snow related energy processes
!  SnowCalc - Calculation of snow and soil storages
!  Evap_SUEWS_Snow - Calculation of evaporation from the snowpack
!  snowRem - Removal of snow my snow clearing
!  SnowDepletionCurve - Calculation of snow fractions
!Last modified
!  LJ 3 Feb 2016  - Changed so that not all surface water freezes in 5-min timestep.
!                    Re-organization of the snow routine due to this change
!                    Calculation of albedo moved from MeltHeat to SnowCalc
!  LJ 27 Jan 2016  - Tabs removed, cleaning of the code
!  HCW 08 Dec 2015 - Added check for no Paved surfaces
!  LJ 14 July 2015 - Code fixed to work with tstep.
!  HCW 06 Mar 2015 - Unused variable 'i' removed.
!  LJ Jan 2015     - Change the calculation from hourly timestep to timestep defined by nsh
!  LJ May 2013     - Calculation of the energy balance for the snowpack was modified
!                        to use qn1_ind_snow(surf)
!=======================================================================================

 subroutine MeltHeat

  use allocateArray
  use snowMod
  use data_in
  use snowMod
  use time
  use moist
  use sues_data
  use gis_data
  
  IMPLICIT NONE

  REAL(KIND(1d0)):: Watfreeze,&        !State of snow on ith surface, total energy exchange (not exported!)
                    cw=4190!,ci=2090   !Specific heat capacity of water
    
  !Initialize snow variables
  snowCalcSwitch=0
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
  SnowRemoval(PavSurf)=0
  SnowRemoval(BldgSurf)=0

  !=========================================================================================
  do is=1,nsurf  !Go each surface type through
    if (sfr(is)/=0) then  !If surface type existing,
      if (SnowPack(is)>0) then  !If snowpack existing, calculate meltwater related water flows

         SnowDepth(is) = (SnowPack(is)/1000)*waterDens/SnowDens(is) !Snow depth in m

         !Calculate meltwater related water flows with hourly degree-day method.

         !These are for snow melting
         if (Temp_C>=0) then
            if (qn1_ind_snow(is)<0) then
               mw_ind(is) = TempMeltFact*Temp_C             !(mm C−1 h−1)*(C) = in mm h-1
            else
               mw_ind(is) = RadMeltFact*(qn1_ind_snow(is))  !(mm m2 W−1 h−1)*(W m-2)= mm h-1 ??
            endif

         else  !Freezing equations
            AdjMeltFact=1  !Relationship between the temperature melt and freezing factors
            mw_ind(is) = TempMeltFact*Temp_C*AdjMeltFact ! in mm h-1
         endif

         !Previous equation give the hourly values, divide these with the timestep number 
         mw_ind(is) = mw_ind(is)/nsh_real

         if (mw_ind(is)>SnowPack(is)) mw_ind(is) = SnowPack(is)!Limited by the previous timestep snowpack

         !-----------------------------------------------------
         ! Heat consumed to snowmelt/refreezing within Tstep.
         ! Converted from mm nsh-1 to mm nsh-1 and to m s-1
         Qm_melt(is) = waterDens*((mw_ind(is)/tstep_real)/1000)*(lvS_J_kg-lv_J_kg)

         !If melt is negative this means freezing water in the snowpack
         if (mw_ind(is)<0) then
            
            FreezMelt(is) = -mw_ind(is) !Save this to variable FreezMelt
            mw_ind(is) = 0
                
            !Freezing water cannot exceed meltwater store
            if (FreezMelt(is)>Meltwaterstore(is)) FreezMelt(is) = Meltwaterstore(is)
             
            !Recalculate melt related energy
            Qm_melt(is) = waterDens*((-FreezMelt(is)/tstep_real)/1000)*(lvS_J_kg-lv_J_kg)
          endif

          !-----------------------------------------------------
          ! If air temperature is above zero, precipitation causes advective heat to the
          ! snowpack. Eq (23) in Sun et al., 1999
  		  ! Calculation done at resolution of the model timestep 
          if (Temp_C>=PrecipLimit.and.Precip>0) then
            Qm_rain(is) = waterDens*cw*(Temp_C-PrecipLimit)*(Precip*0.001/tstep_real)  !in W m-2
            if (Qm_rain(is)<0) then !Can only be positive
               Qm_rain(is) = 0
            else
               rainOnSnow(is) = Precip !Save information on the rain on snow event
            endif
          endif
   
        endif !End if snowpack
        
        !=================================================================
         
        !Freeze surface water state if cold enough.
        if (Tsurf_ind(is)<0.and.state(is)>0) then
    
           snowCalcSwitch(is)=1 !If water on ground this forms ice and snow calculations are made
    
           !Other surfaces than water treated first
           if (is/=WaterSurf) then

              !FreezState(is) = state(is) 
              !Previously all state could freeze in 5-min timestep. Now we calculate how much water
              !can freeze in timestep (see watersurf exaplanation)
              FreezState(is) = 100*(0-Temp_C)/(waterDens*(lvS_J_kg-lv_J_kg))

              !The amount of freezing water cannot be greater than the surface state
              if (FreezState(is)>state(is)) FreezState(is) = state(is)

              if (snowpack(is)==0.or.snowfrac(is)==0) then !Snowpack forms
                   FreezStateVol(is) = FreezState(is)
              else                                         !There is snow already on ground
                   FreezStateVol(is) = FreezState(is)*(1-snowFrac(is))/snowFrac(is)
              endif

              ! If the amount of freezing water is very small and there is state left to the ground
              ! no freezing of water will take place
              if (FreezStateVol(is)<0.00000000001.and.FreezState(is)<state(is)) then
                 FreezState(is) = 0
                 FreezStateVol(is) = 0
              endif

              !Calculate the heat exchange in W m-2
              Qm_freezState(is) = -waterDens*(FreezState(is)/tstep_real/1000)*(lvS_J_kg-lv_J_kg)

           !Water surface separately
           else
              !Calculate average value how much water can freeze above the water areas
              !Equation is -hA(T-T0) = rhoV(Cp+dT +Lf) in 5-min timestep
              !h=convective heat trasnfer,A, area of water,rwo water density,V volume, dT temperature difference
              !before and end of the 5-min period. dT equals zero, h=100 and when multiplied with Area, the equation
              !simplyfies to the this. LJ 14 July 2015
              Watfreeze = 100*(0-Temp_C)/(waterDens*(lvS_J_kg-lv_J_kg))
              FreezState(is) = Watfreeze
              Qm_freezState(is) = -waterDens*(Watfreeze/tstep_real/1000)*(lvS_J_kg-lv_J_kg)
           endif
    
        endif

       !======================================================================
       ! Define if any snowmelt calculations are made: snowpack existing,
       ! freezing occuring on ground or from precip
        if (is/=WaterSurf) then
          if (snowPack(is)>0.or.(Precip>0.and.Tsurf_ind(is)<0)) then
            snowCalcSwitch(is)=1
          endif
        else       !Water surface separately
          if (snowPack(WaterSurf)>0.or.FreezState(WaterSurf)>0) then
            snowCalcSwitch(WaterSurf)=1
          endif
        endif

        !Update snow density of each surface
        if (Precip>0.and.Tsurf_ind(is)<0.and.SnowPack(is)>0) then
            !write(*,*) SnowDens(is),snowPack(is),Precip,SnowDensMin
            SnowDens(is) = SnowDens(is)*snowPack(is)/(snowPack(is)+Precip)+SnowDensMin*Precip/(snowPack(is)+Precip)
            !write(*,*) SnowDens(is)
            !pause
        endif

        !Weighted variables for the whole area
        mwh = mwh + mw_ind(is)*sfr(is)*snowFrac(is)        !Snowmelt
        fwh = fwh + FreezMelt(is)*sfr(is)*snowFrac(is)     !Freezing water
        Qm = Qm + Qm_melt(is)*sfr(is)*snowFrac(is)         !Energy consumed to the melt/freezing.
        QmRain = QmRain + Qm_rain(is)*sfr(is)*snowFrac(is) !Rain on snow
        QmFreez=QmFreez+deltaQi(is)*sfr(is)*snowFrac(is)+Qm_freezState(is)*sfr(is)*(1-snowFrac(is)) !Freezing water
    endif
     
  enddo !End surface type

  !Update snow albedo to its maximum value if precipitation exists
  if (Precip>0.and.sum(snowPack)>0.and.Temp_C<0) then

    CumSnowfall=CumSnowfall + Precip

    if (CumSnowfall>PrecipLimitAlb) then
    !write(*,*) CumSnowfall,Precip,PrecipLimitAlb,PrecipLimitAlb/nsh_real
   !pause
      SnowAlb=SnowAlbMax
      CumSnowfall=0
    endif
  else
    CumSnowfall=0
  endif


 end subroutine MeltHeat


!===============================================================================================
!===============================================================================================
 subroutine SnowCalc
  !Calculation of snow and water balance on 5 min timestep. Treats snowfree and snow covered
  !areas separately. Weighting is taken into account in the overall values.
  !Last modified:
  !  LJ in 6 May 2015 - Modified to run with timestep
  !  HCW 06 Mar 2015 - Unused variable 'i' removed.
  !  HCW 26 Jan 2015 - Added weekday/weekend option for snow clearing profiles
  !  LJ in 24 May 2013
  !========================================================================
  
  use allocateArray
  use defaultNotUsed
  use data_in    
  use gis_data
  use mod_k
  use mod_z
  use moist
  use snowMod
  use sues_data
  use thresh
  use time
  
  IMPLICIT NONE
             
  REAL(KIND(1d0)):: Evap_SUEWS_Snow,&
                    MeltExcess,&      !Excess melt water that needs to leave snowpack
                    snowTotInit,&
                    EvPart,&
                    runoffTest,&
                    snowFracFresh,&   !Snow fraction for newly formed snowpack.
                    snowFracOld,&
                    FWC               !Water holding capacity of snow in mm
  REAL(KIND(1d0)):: SnowDepletionCurve

  integer:: iu                        !1=weekday OR 2=weekend
  !========================================================================
  !Initialize variables for the calculation of water storages and evaporation

  ! Use weekday or weekend snow clearing profile
  iu=1     !Set to 1=weekday
  if(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) iu=2  !Set to 2=weekend

  runoffSnow(is)=0 !Initialize for runoff caused by snowmelting
  runoff(is)=0
  runoffSoil(is)=0
  chang(is)=0
  changSnow(is)=0
  runoffTest=0
  SnowToSurf(is)=0
  EvPart=0
  ev=0
  snowFracFresh=0
  snowFracOld=0

  !Initial snowpack + meltwater in it  
  snowTotInit=snowPack(is)+MeltWaterStore(is) 
  
  !Calculate water holding capacity (Jin et al. 1999)
  if (SnowDens(is)>=200) then
     WaterHoldCapFrac=CRWmin
  else
     WaterHoldCapFrac=CRWmin+(CRWmax-CRWmin)*(200-SnowDens(is))/200  
  endif

  !======================================================================
  ! Calculate evaporation from snowpack and snow free surfaces (in mm)
  if (snowFrac(is)<1) call Evap_SUEWS !ev and qe for snow free surface out
  
  if (snowFrac(is)>0) then
      ev_snow(is) = Evap_SUEWS_Snow(Qm_Melt(is),Qm_rain(is),lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RAsnow,&
                                    psyc_hPa,tstep,avcp,sIce_hPa)
  endif

  !If not enough water for evaporation in impervious surfaces,
  !evaporation is taken from pervious surfaces
  if (is>2) then
     if  (PervFraction/=0) then
       EvPart=(SurplusEvap(PavSurf)*sfr(PavSurf)+SurplusEvap(BldgSurf)*sfr(BldgSurf))/PervFraction
     endif
  endif


  !============================================================================
  !1) Surface is fully covered with snow or the snowpack forms (equally distributed) 
  !   on the current Tstep
  !2) Both snow and snow free surface exists
  !3) Water surface is treated separately
  if (is==WaterSurf.and.sfr(WaterSurf)>0) GO TO 606

  !The calculations are divided into 3 parts
  ! 1) Surface is fully covered with snow
  ! 2) Surface is not fully covered with snow


 ! if (id==313.and.it==22.and.imin==35) then
 !    write(*,*) is,SnowPack(is), snowFrac(is),ev_snow(is),ev,EvPart,rainOnSnow(is),mw_ind(is)
 ! endif

  !1)------------------------------------------------------------------
  if (SnowPack(is)>0.and.snowFrac(is)==1) then

     ev_snow(is)=ev_snow(is)+EvPart !Evaporation surplus

     !(Snowfall per interval+freezing of melt water and surface state) - (meltwater+evaporation from snowpack)
     changSnow(is)=(Precip+freezMelt(is))-(mw_ind(is)+ev_snow(is)) !Calculate change in snowpack (in mm)

     !If rain on snow event, add this water to meltwaterstore
     if (rainOnSnow(is)>0) then
       changSnow(is)=changSnow(is)-Precip
       MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)
     endif

     SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack

     !---------If snowpack exists after the state calculations
     if (snowPack(is)>0) then

     !Add melted water to meltstore and freeze water according to freezMelt(is)
     MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is) - freezMelt(is)

     !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
     FWC = WaterHoldCapFrac*snowPack(is)

     !If FWC is exceeded, excess meltwater (MeltExcess) will leave from the snowpack
     if (MeltWaterStore(is)>=FWC) then
        MeltExcess = 0                      !Initialize the excess meltwater
        MeltExcess = MeltWaterStore(is)-FWC !Calculate the exceess water
        MeltWaterStore(is) = FWC            !Update the meltwaterstore to the maximum it can hold
        runoffSnow(is) = runoffSnow(is) + MeltExcess
     endif

     !At the end of the hour calculate possible snow removal
     if (SnowProf(it,iu)==1.and.is<3.and.(imin==(nsh_real-1)/nsh_real*60))  call snowRem

     !----------If snowPack is negative, it melts at this timestep
   elseif (SnowPack(is)<0) then

     !If freezing meltwater inside this timestep, remove it from the MeltWaterStore
     MeltWaterStore(is)=MeltWaterStore(is)-freezMelt(is)+mw_ind(is)+SnowPack(is)
     SnowPack(is)=0.0   !Set the snow pack and snow
     snowFracOld=1
     snowFrac(is)=0
     if (id==313.and.it==22.and.imin==35) then
       write(*,*) is,SnowPack(is), snowFrac(is),ev_snow(is),rainOnSnow(is),mw_ind(is),MeltWaterStore(is)
      ! pause
     endif
     if (MeltWaterStore(is)<0) then !Not enough water in the meltwater store,
        ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased.??
        if (ev_snow(is)<0) ev_snow(is)=0
        changSnow(is)=changSnow(is)+MeltWaterStore(is)
        MeltWaterStore(is)=0
     else
        chang(is)=MeltWaterStore(is)  !Meltwater goes to surface state as no snow exists anymore
        state(is)=state(is)+chang(is)
        MeltWaterStore(is)=0
     endif
   endif !snowpack negative or positive

   !---------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------
   !2) Surface not completely covered with snow
  elseif (snowFrac(is)<1) then

   !Snow calculations: snowpack can exist.
   if (SnowPack(is)>0) then
      ev_snow(is)=ev_snow(is)+EvPart !Evaporation surplus

      !----Snowpack water balance for the whole surface area. In reality snow depth = snowPack/snowFrac(is)
      !(Snowfall per interval+freezing of melt water and surface state) - (meltwater+evaporation from snowpack)
      changSnow(is)=(Precip+freezMelt(is)+freezStateVol(is))-(mw_ind(is)+ev_snow(is)) !Calculate change in snowpack (in mm)

      !If rain on snow event, add this water to meltwaterstore
      if (rainOnSnow(is)>0) then
         changSnow(is)=changSnow(is)-Precip
         MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)
      endif
      SnowPack(is)=SnowPack(is)+changSnow(is)

      !The fraction of snow will get a value of 1 (ie full snow cover):
      !Surface state is dry but precipitation occurs, no precipitation but all state freezes at a single timestep,
      !There is both precipitation and all surface state freezes
      !In this case no snow free state calculations are made as snowFrac(is)==1
      if ((Precip>0.and.state(is)==0).or.(Precip==0.and.FreezState(is)==state(is)).or.&
         (Precip>0.and.FreezState(is)==state(is))) then
         snowFrac(is)=1
         SnowDens(is)=SnowDensMin
         state(is)=0
         changSnow(is)=changSnow(is)-ev-drain(is)
      endif
   endif

   !Snowpack can also form at the current timestep (2). If this forms purely from snowfall or/and all water at surface freezes,
   !the whole surface will be covered with snow. If there is water on ground this snowfall can immediately melt
   !and in this case the snow fraction is not necessarily 1 but its information is saved to snowFracFresh that
   !is taken into account in snow fraction after calculation of state.
   if (snowpack(is)==0.and.Tsurf_ind(is)<0) then

    !The fraction of snow will get a value of 1 (ie full snow cover):
    !Surface state is dry but precipitation occurs, no precipitation but all state freezes at a single timestep,
    !There is both precipitation and all surface state freezes
    if ((Precip>0.and.state(is)==0).or.(Precip==0.and.FreezState(is)==state(is)).or.&
       (Precip>0.and.FreezState(is)==state(is))) then

       ev=ev+EvPart
       changSnow(is)=Precip+FreezStateVol(is)-ev-drain(is)
       SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack
       snowFrac(is)=1                           !If precipitation occurs, snow fraction = 1
       !iceFrac(is)=0
       SnowDens(is)=SnowDensMin
       state(is)=0
    endif

    if (FreezState(is)>0.and.FreezState(is)<state(is)) then
      changSnow(is)=(Precip+freezStateVol(is))
      SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack

      snowFracFresh=SnowDepletionCurve(is,snowPack(is),snowD(is))
      iceFrac(is)=1
      SnowDens(is)=SnowDensMin
      if (snowFracFresh<0.001) snowFracFresh=0.001

    endif
   endif


   !---------If snowpack exists after the state calculations
   if (snowPack(is)>0) then

      !Add melted water to meltstore and freeze water according to freezMelt(is)
      MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is) - freezMelt(is)

      !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
      FWC = WaterHoldCapFrac*snowPack(is)

      !If FWC is exceeded, excess meltwater (MeltExcess) will leave from the snowpack
      if (MeltWaterStore(is)>=FWC) then
         MeltExcess = 0                      !Initialize the excess meltwater
         MeltExcess = MeltWaterStore(is)-FWC !Calculate the exceess water
         MeltWaterStore(is) = FWC            !Update the meltwaterstore to the maximum it can hold

         !If the fraction of snow is greater than 0.8 or if the surface is is buildings,
         !the excess water will directly go to runoff. Otherwise it will flow to the
         !snow free area via SnowToSurf(is)
         if ((snowFrac(is)>0.8.and.is/=BldgSurf).or.(is==BldgSurf)) then
            runoffSnow(is) = runoffSnow(is) + MeltExcess
         else
            SnowToSurf(is) = SnowToSurf(is) + MeltExcess*snowFrac(is)/(1-snowFrac(is))
         endif
       endif

       !At the end of the hour calculate possible snow removal
       if (SnowProf(it,iu)==1.and.is<3.and.(imin==(nsh_real-1)/nsh_real*60))  call snowRem

       !----------If snowPack is negative, it melts at this timestep
    elseif (SnowPack(is)<0) then

     !If freezing meltwater inside this timestep, remove it from the MeltWaterStore
     MeltWaterStore(is)=MeltWaterStore(is)-freezMelt(is)+mw_ind(is)+SnowPack(is)

     SnowPack(is)=0.0   !Set the snow pack and snow
     snowFrac(is)=0.0

     if (MeltWaterStore(is)<0) then !Not enough water in the meltwater store,
        ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased.??
        if (ev_snow(is)<0) ev_snow(is)=0
        changSnow(is)=changSnow(is)+MeltWaterStore(is)
        MeltWaterStore(is)=0
     else
        chang(is)=MeltWaterStore(is)*snowFrac(is)/(1-snowFrac(is))  !Meltwater goes to surface state as no snow exists anymore
        SnowToSurf(is)=SnowToSurf(is)+chang(is)
        MeltWaterStore(is)=0
     endif
   endif !snowpack negative or positive


  !--------
  !Next the snow free surface (3). Calculations only done if snowfraction is smaller than 1
  if ((is==PavSurf.or.is==BldgSurf).and.snowFrac(is)<1) then  !Impervious surfaces (paved, buildings)

     !Surface store update. If precipitation is greater than the threshold, the exceeding water
     !goes directly to runoff
     if (precip>IPThreshold_mmhr/nsh_real) then
        !runoff = runoff + (precipitation+water from the snow surface+water from other surfaces-the thereshold limit)
        runoff(is)=runoff(is)+(Precip+SnowToSurf(is)+AddWater(is)-IPThreshold_mmhr/nsh_real)
        chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev+freezState(is))
     else
        !Add precip and water from other surfaces and remove drainage, evap and freezing of state
        chang(is)=Precip+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is))
     endif

     state(is)=state(is)+chang(is) !Change in state (for whole surface area areasfr(is))

     !Add water from impervious grids
     ! Check sfr/=0 added HCW 08 Dec 2015
     if (is==PavSurf.and.sfr(PavSurf)/=0) state(is)=state(is)+(addImpervious)/sfr(PavSurf)

     runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is) !Drainage (not flowing to other surfaces) goes to runoff

     if(state(is)<0.0) then  !Surface state cannot be negative
        SurplusEvap(is)=abs(state(is)) !take evaporation from other surfaces in mm
        ev = ev-SurplusEvap(is)
        state(is)=0.0
     endif

  elseif(is>=3.and.snowFrac(is)<1) then ! Pervious surfaces (conif, decid, grass unirr, grass irr)

     ev=ev+EvPart

     !Change in water stores
     if (Precip+addVeg*(sfr(is)/VegFraction)>(IPThreshold_mmhr/nsh_real)) then !if 5min precipitation is larger than 10 mm
        runoff(is)=runoff(is)+(Precip+addVeg*(sfr(is)/VegFraction)+SnowToSurf(is)+AddWater(is)-(IPThreshold_mmhr/nsh_real))
        chang(is)=(IPThreshold_mmhr/nsh_real)-(drain(is)+ev+freezState(is))
     else
        chang(is)=Precip+addVeg*(sfr(is)/VegFraction)+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is))
     endif

     state(is)=state(is)+chang(is)

     !Add water in soil store only if ground is not frozen
     if (Temp_C>0) then
        soilmoist(is)=soilmoist(is)+Drain(is)*AddWaterRunoff(is)*(1-snowFrac(is))
     else
        runoff(is)=runoff(is)+Drain(is)*AddWaterRunoff(is)
     endif

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

  endif !Surface fraction


  !-------------------------------------------------------------------------------------------------------------------
  !!Update the fraction of snow in the case of falling snow
  if (snowFracFresh>0) snowFrac(is)=snowFracFresh

  !Calculate change in snowpack and state for the respective surface areas
  !ch_per_interval=ch_per_interval+(state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is))
  surf_chang_per_tstep=surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is))
  chSnow_per_interval=chSnow_per_interval+((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)*max(snowFrac(is),snowfracOld)

  !Add evaporation to total
  if (is==BldgSurf.or.is==PavSurf) then
    ev_per_tstep=ev_per_tstep+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
    qe_per_tstep=qe_per_tstep+ev_snow(is)*lvS_J_kg*sfr(is)*snowFrac(is)&
                       +ev*lv_J_kg*sfr(is)*(1-snowFrac(is))
  else
    ev_per_tstep=ev_per_tstep+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*snowFrac(is)
    qe_per_tstep=qe_per_tstep+ev_snow(is)*lvS_J_kg*sfr(is)*snowFrac(is)+ev*lv_J_kg*sfr(is)*(1-snowFrac(is))
  endif
  
  !========RUNOFF=======================
   
  !Add runoff to pipes
  runoffPipes=runoffPipes+runoffSnow(is)*sfr(is)*snowFrac(is)+runoff(is)*sfr(is)*(1-snowFrac(is))+runoffTest*sfr(is)
  call updateFlood
  runoff_per_tstep=runoff_per_tstep+runoffSnow(is)*sfr(is)*snowFrac(is)+runoff(is)*sfr(is)*(1-snowFrac(is))&
                         +runoffTest*sfr(is)

  !===Update snow depth, weighted SWE, and Mwstore
  if (SnowDens(is)/=0) then
     SnowDepth(is) = SnowPack(is)*waterDens/SnowDens(is)
  endif

  ! Calculate overall snow water equivalent
  swe = swe + SnowPack(is)*sfr(is)*snowFrac(is)
  MwStore = MwStore + MeltWaterStore(is)*sfr(is)*snowFrac(is)


  !Calculate new snow fraction here
  if (SnowFractionChoice==2.and.it==23.and.imin==(nsh_real-1)/nsh_real*60) then
     !if (snowPack(is)>0.and.mw_ind(is)>0) then
     if (snowPack(is)>0) then
        snowFrac(is) = SnowDepletionCurve(is,snowPack(is),snowD(is))
        if (snowFrac(is)<0.001) snowFrac(is)=0.001  !The snow fraction minimum is 1% of the surface
     elseif (snowPack(is)>0.and.Precip>0.and.Tsurf_ind(is)<0.and.state(is)==0) then
        snowFrac(is) = 1
     endif

  elseif (snowPack(is)==0) then
    snowFrac(is)=0
  endif

 return

 !==========================================================================
 !WATERBODY is treated separately as state always below ice if ice existing
 !Calculate change in snowpack
606 changSnow(WaterSurf)=(Precip+freezMelt(WaterSurf)+freezState(WaterSurf))-&
                         (mw_ind(WaterSurf)+ev_snow(WaterSurf))
                         
    SnowPack(WaterSurf)=SnowPack(WaterSurf)+changSnow(WaterSurf) !Update snowpack
    state(WaterSurf)=state(WaterSurf)+FlowChange-freezState(WaterSurf)  !Update state below ice
     
    !If snowpack exists
    if (snowPack(WaterSurf)>0) then  
    
       !Add melted water to meltstore and freeze water according to freezMelt(is)
       MeltWaterStore(WaterSurf)=MeltWaterStore(WaterSurf)+mw_ind(WaterSurf)-freezMelt(WaterSurf)
      
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
        MeltWaterStore(WaterSurf) = MeltWaterStore(WaterSurf)-freezMelt(WaterSurf) &
                                    + mw_ind(WaterSurf)
       
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
     !ch_per_interval=ch_per_interval+(state(WaterSurf)-stateOld(WaterSurf))*sfr(WaterSurf)
     surf_chang_per_tstep=surf_chang_per_tstep+(state(WaterSurf)-stateOld(WaterSurf))*sfr(WaterSurf)
     
     !Evaporation
     ev_per_tstep=ev_per_tstep+ev*sfr(WaterSurf)+ev_snow(WaterSurf)*sfr(WaterSurf)
     qe_per_tstep=qe_per_tstep+ev_snow(WaterSurf)*lvS_J_kg*sfr(WaterSurf)+ev*lv_J_kg*sfr(WaterSurf) 
     runoff_per_tstep=runoff_per_tstep+(runoff(is)*sfr(is)) !The total runoff from the area
     
     if (snowPack(WaterSurf)>0) then     !Fraction only 1 or 0
         snowFrac(WaterSurf)=1
     else
         snowFrac(WaterSurf)=0
     endif

  end subroutine SnowCalc
  
  
  !==========================================================================
  !==========================================================================
  !Calculates evaporation from snow surface (ev_snow).

  FUNCTION Evap_SUEWS_Snow(Qm,QP,lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RAsnow,psyc_hPa,&
                           tstep,avcp,sIce_hPa) RESULT(ev_snow)
    
   IMPLICIT NONE

   !INPUT
   real (Kind(1d0))::Qm,QP,&        !melt heat, advect. heat
                     lvS_J_kg,avdens,avRh,&   !latent heat of sublimation, air density,relative humidity,
                     Press_hPa,Temp_C,&       !air pressure, air temperature
                     RAsnow,psyc_hPa,&        !aerodyn res snow, psychometric constant, type of evaporation calculation
                     avcp,sIce_hPa            !spec. heat, satured curve on snow

   !OTHER VARIABLES
   real (Kind(1d0))::e_snow,&     !PM equation obe line
                     sae_snow,&   !s * (Available energy)
                     qe_snow,&    !Latent heat flux
                     ev_snow,&    !Evaporation
                     vdrcIce,&    !Vapour pressure deficit
                     esIce_hPa,&  !Saturation vapor pressure over ice
                     EaIce_hPa,&  !Vapour pressure
                     tlv_sub,&    !Latent heat for sublimation
                     tstep_real   !timestep as real

   real (Kind(1d0)):: sat_vap_pressIce !Function

   integer:: tstep,from=1
   !-----------------------------------------------------
   
   tstep_real = real(tstep,kind(1d0))

   sae_snow=sIce_hPa*(Qp-Qm)   !Calculate the driving parameter in calculation of evaporation. Järvi et al. (2015)

   esIce_hPa= sat_vap_pressIce(Temp_C,Press_hPa,from) !Saturation vapor pressure over ice
   EaIce_hPa=avRh/100*esIce_hPa                       !Vapour pressure of water
   vdrcIce=(esIce_hPa-eaIce_hpa)*avdens*avcp          !Vapour pressure deficit
   tlv_sub=lvS_J_kg/tstep_real                        !Latent heat for sublimation
  
   e_snow=sae_snow+vdrcIce/RAsnow                     !PM equation
   qe_snow=e_snow/(sIce_hPa+psyc_hPa)                 !Latent heat (W/m^2)
   ev_snow=qe_snow/tlv_sub                            !Evaporation (in mm)
  
   return
   
  END FUNCTION Evap_SUEWS_Snow
  
  !==========================================================================
  !==========================================================================
  ! Calculates mechanical removal of snow from roofs ans roads
  subroutine snowRem

   use allocateArray 
   use snowMod
   use sues_data
   use time
  
   IMPLICIT NONE

  !write(*,*) is, SnowPack(is),SnowLimPaved,SnowLimBuild

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
  !write(*,*) is, SnowPack(is),SnowLimPaved,SnowLimBuild
 !pause
  end subroutine snowRem
 
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  FUNCTION SnowDepletionCurve(is,swe,sweD) RESULT(asc)
  !This function calculates surface coverage of snow according to the
  !depletion curves in Valeo and Ho (2004).
  !INPUT: is   Surface type number
  !       swe  Snow water content
  !       sweD Limit for
  
   use allocateArray
  
   IMPLICIT  NONE
  
   INTEGER::is
   REAL (KIND(1d0))::asc,sweD,swe


   !Impervious surface
   if (is==PavSurf) then
       
     if (swe<=sweD) then      !Snow water equivalent below threshold
        asc=((swe/sweD))**2
     else
        asc=1
     endif
     
   !Bldgs surface  
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