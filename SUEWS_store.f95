 subroutine soilstore
!------------------------------------------------------------------------------
!Calculation of storage change
!LJ 6 May 2015
! - Calculations of the piperunoff exceedings moved to separate subroutine updateFlood.
!   Now also called from snow subroutine
!Evaporation is modified using EvapPart
! - when no water on impervious surfaces, evap occurs above pervious surfaces instead
!Rewritten by HCW 12 Feb 2015
! Old variable 'p' for water input to the surface renamed to 'p_mm'
! All water now added to p_mm first, before threshold checks or other calculations
! Water from other grids now added to p_mm (instead of state for impervious surfaces)
! Removed division of runoff by nsh, as whole model now runs at the same timestep
! Adjusted transfer of ev between surfaces to conserve mass (not depth)
! Volumes used for water transport between grids to account for SurfaceArea changing between grids
! Added threshold check for state(WaterSurf) - was going negative
!Last modified HCW 09 Feb 2015
! Removed StorCap input because it is provided by module allocateArray
! Tidied and commented code
!Modified by LJ in November 2012: 
! P>10 was not taken into account for impervious surfaces - Was fixed.
! Above impervious surfaces possibility of the state to exceed max capacity was limited
!  although this should be possible - was fixed
!Modified by LJ 10/2010
!Rewritten mostly by LJ in 2010
! To do:
!	- Finish area normalisation for RG2G & finish coding GridConnections
!	- What is the 10 mm hr-1 threshold for?
! 	- Decide upon and correct storage capacities here & in evap subroutine
! 	- FlowChange units should be mm hr-1 - need to update everywhere
!	- Add SurfaceFlood(is)?
!	- What happens if sfr(is) = 0 or 1?
!	- Consider how irrigated trees actually works...
!------------------------------------------------------------------------------
  
  use allocateArray
  use data_in
  use defaultNotUsed
  use gis_data  
  use sues_data
  use thresh
  use time
  
  IMPLICIT NONE
  
  !real(kind(1d0)),dimension(nsurf):: SurfaceFlood   !Stores flood water when surface state exceeds storage capacity [mm]
  
  real(kind(1d0)):: EvPart   !Extra evaporation [mm] from impervious surfaces which cannot happen due to lack of water
  ! Initialise extra evaporation to zero
  EvPart=0   
   
  !SurfaceFlood(is) = 0 !!This probably needs to be carried over between timesteps, but reset for now
  
  !==================================================================
  ! Combine water inputs to the current surface
    
  ! ---- If surface is irrigated, add external water use to precip ----
  ! Add external water use for each surface type
  if(is==ConifSurf) then
     p_mm=pin+wu_EveTr
  elseif(is==DecidSurf) then
     p_mm=pin+wu_DecTr
  elseif(is==GrassSurf) then  
     p_mm=pin+wu_Grass
  else
     p_mm=pin
  endif
    
  ! ---- Add water from other surfaces within the same grid (RS2S) ----
  ! AddWater is the water supplied to the current surface from other surfaces 
  !  i.e. drain*WaterDist (see SUEWS_ReDistributeWater)
  p_mm=p_mm+AddWater(is)


 
  !==== Impervious surfaces (Paved, Buildings) ======================
  if(is==PavSurf.or.is==BldgSurf) then
     
     ! ---- Add water from neighbouring grids (RG2G) ----
     ! Add to PavSurf only, as water cannot flow onto buildings
     if (is==PavSurf) p_mm=p_mm+addImpervious/sfr(PavSurf)

     ! Calculate change in surface state (inputs - outputs)
     chang(is)=p_mm-(drain(is)+ev)	   
     
     ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
     !  and does not affect state
     if(p_mm>IPThreshold_mmhr/nsh_real) then
        runoff(is)=runoff(is)+(p_mm-IPThreshold_mmhr/nsh_real)
        chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev)	
     endif
     
     ! Calculate updated state using chang
     state(is)=state(is)+chang(is)

     ! Check state is within physical limits between zero (dry) and max. storage capacity
     if(state(is)<0.0) then   ! Cannot have a negative surface state
        ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
        ! Allow evaporation only until surface is dry (state(is)=0); additional evaporation -> evaporation surplus
        SurplusEvap(is)=abs(state(is))   !Surplus evaporation is that which tries to evaporate non-existent water
 	    ev = ev-SurplusEvap(is)	  !Limit evaporation according to water availability
        state(is)=0.0			  !Now surface is dry
    ! elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
    !    !! If state exceeds the storage capacity, then the excess goes to surface flooding
    !    !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
    !    runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
    !    state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
     endif     
       
     ! Recalculate change in surface state from difference with previous timestep
     chang(is) = state(is)-stateOld(is)
     
     ! Runoff -------------------------------------------------------
     ! For impervious surfaces, some of drain(is) becomes runoff
     runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is)   !Drainage (that is not flowing to other surfaces) goes to runoff 
          
     !So, up to this point, runoff(is) can have contributions if
     ! p_mm > ipthreshold (water input too fast)
     ! state > surf(6,is) (net water exceeds storage capacity)
     ! WaterDist specifies some fraction of drain(is) -> runoff
 
  !==== Pervious surfaces (Conif, Decid, Grass, BSoil, Water) =======
  elseif(is>=3) then

     ! Transfer evaporation surplus from impervious surfaces to pervious surfaces
     if(PervFraction/=0) then   !If pervious surfaces exist 
        EvPart=(SurplusEvap(PavSurf)*sfr(PavSurf)+SurplusEvap(BldgSurf)*sfr(BldgSurf))/PervFraction
     else   !If no pervious surface, SurplusEvap cannot be transferred and this evap cannot happen (will increase QH instead)
        EvPart=0
     endif
     ! Add surplus evaporation to ev for pervious surfaces
     ev=ev+EvPart

     !==== For Conif, Decid, Grass, BSoil surfaces ==================
     if (is/=WaterSurf) then
        
        ! ---- Add water from neighbouring grids (RG2G) ----
        ! Add to Grass and BSoil only, as water cannot flow onto trees
        if (is==GrassSurf.or.is==BSoilSurf) then 
           if ((sfr(GrassSurf)+sfr(BSoilSurf))/=0) then
              p_mm=p_mm+addVeg/(sfr(GrassSurf)+sfr(BSoilSurf))
           endif
        endif
        
        ! Calculate change in surface state (inputs - outputs)
	    chang(is)=p_mm-(drain(is)+ev)

	    ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
        !  and does not affect state
        if (p_mm>IPThreshold_mmhr/nsh_real) then   
           runoff(is)=runoff(is)+(p_mm-IPThreshold_mmhr/nsh_real)
           chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev)   
	    endif
	
	    ! Calculate updated state using chang
     	state(is)=state(is)+chang(is)

        ! Check state is within physical limits between zero (dry) and max. storage capacity
	    if(state(is)<0.0) then   ! Cannot have a negative surface state
	    ! If there is not sufficient water on the surface, then remove water from soilstore
 	    ! Allow evaporation until soilmoist is depleted and surface is dry
           if((soilmoist(is)+state(is))>=0) then
              soilmoist(is)=soilmoist(is)+state(is)
              state(is)=0.0
           else   
	      ! If there is not sufficient water on the surface or soilstore, then don't allow this evaporation to happen
              ev=ev-abs(state(is))   !Limit evaporation according to water availability
              state(is)=0.0	     !Now surface is dry	
	   endif
	   !! What about if there is some water in soilstore, but not enough to provide all the water for evaporation??	   
	   !! Is this saying water can be evaporated from the soilstore as easily as from the surface??
  	   !elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
	    !   !! If state exceeds the storage capacity, then the excess goes to surface flooding
	   !   !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
	   !   runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
	   !   state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
	endif     
	       
	! Recalculate change in surface state from difference with previous timestep
    chang(is) = state(is)-stateOld(is)

    !Where should this go? Used to be before previous part!!
	! Soilmoist -------------------------------------------------
	! For pervious surfaces (not water), some of drain(is) goes to soil storage
	soilmoist(is)=soilmoist(is)+drain(is)*AddWaterRunoff(is) !Drainage (that is not flowing to other surfaces) goes to soil storages
		
        ! If soilstore is full, the excess will go to runoff
     	if(soilmoist(is)>soilstoreCap(is)) then			!! Should this also go to flooding of some sort?
           runoff(is)=runoff(is)+(soilmoist(is)-soilstoreCap(is))
           soilmoist(is)=soilstoreCap(is)
     	elseif (soilmoist(is)<0) then   !!But where does this lack of water go? !!Can this really happen here??
           call ErrorHint(62,'SUEWS_store: soilmoist(is) < 0 ',soilmoist(is),NotUsed,is)
           ! Code this properly - soilmoist(is) < 0 shouldn't happen given the above loops
           !soilmoist(is)=0   !Groundwater / deeper soil should kick in 
     	endif
                
     
     !==== Water surface ========================================
     elseif (is==WaterSurf) then
       
        if(sfr(WaterSurf)/=0)then
           
           ! ---- Add water from neighbouring grids (RG2G) ----
           p_mm=p_mm+addWaterBody/sfr(WaterSurf)
     	 
     	   ! Calculate change in surface state (inputs - outputs) 
     	   ! No drainage for water surface
     	   ! FlowChange is the difference in input and output flows [mm hr-1]   !!Should this really be a constant??
           chang(is)=p_mm+FlowChange/nsh_real-(ev) 
        	 
	   ! Calculate updated state using chang
           state(is)=state(is)+chang(is)			        
			        
           ! Check state is within physical limits between zero (dry) and max. storage capacity        	 
           if(state(is)<0.0) then   ! Cannot have a negative surface state
              ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
              ev=ev-abs(state(is))   !Limit evaporation according to water availability
              state(is)=0.0	     !Now surface is dry	
 	   !elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
           !   !! If state exceeds the storage capacity, then the excess goes to surface flooding
           !   !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
           !   runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
           !   state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
           endif     
       
           ! Recalculate change in surface state from difference with previous timestep
           chang(is) = state(is)-stateOld(is)        	  
        	 
           ! If state exceeds limit, then excess goes to runoff (currently applies to water surf only)
           if (state(WaterSurf)>StateLimit(WaterSurf)) then
              runoff(WaterSurf)=runoff(WaterSurf)+(state(WaterSurf)-StateLimit(WaterSurf))
              state(WaterSurf)=StateLimit(WaterSurf)
              runoffWaterBody=runoffWaterBody+runoff(WaterSurf)*sfr(WaterSurf)
           else  
              state(WaterSurf)=state(WaterSurf)+surplusWaterBody
              if (state(WaterSurf)>StateLimit(WaterSurf)) then
                 runoffWaterBody=runoffWaterBody+(state(WaterSurf)-StateLimit(WaterSurf))*sfr(WaterSurf)
                 state(WaterSurf)=StateLimit(WaterSurf)
              endif
           endif
           
           ! Recalculate change in surface state from difference with previous timestep
           chang(is) = state(is)-stateOld(is)        
        endif
     
     endif   !end of WaterSurf
   
  endif   !end of different surfaces
  
  !==================================================================
  !==== RUNOFF ======================================================
 
  ! Need to consider areas here - SurfaceArea may vary between grids too
  ! - also implement where water for next surface is calculated (RunoffFromGrid subroutine)
  ! Calculations of the piperunoff exceedensances moved to separate subroutine so that from snow same 
  ! calculations can be made. LJ in May 2015
   
  if(is<WaterSurf) then   !Not for water body

     ! Add runoff to pipes 
     runoffPipes=runoffPipes+(runoff(is)*sfr(is))
     call updateFlood
  endif

  runoff_per_interval=runoff_per_interval+(runoff(is)*sfr(is)) !The total runoff from the area !!Check (HCW)
  
 end subroutine soilstore
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 subroutine updateFlood

  use allocateArray
  use sues_data

  IMPLICIT NONE

  ! If pipe capacity is full, surface runoff occurs
  ! N.B. this will happen each loop (replicates pipes filling up)
  if(runoffPipes>PipeCapacity) then

    !------Paved and building surface
    if(is==PavSurf.or.is==BldgSurf) then
        if(sfr(WaterSurf)>0.0000001) then
           ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
           ! RunoffToWater is specified in SUEWS_SiteSelect.txt
           runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)*(1-RunoffToWater)
           surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
        else
           ! Otherwise, all flood water must go to runoff
           runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)
        endif
                   !------other surfaces
     elseif(is>=ConifSurf.and.is<=BSoilSurf) then
        if(sfr(WaterSurf)>0.0000001) then
          ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
          runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)*(1-RunoffToWater)
          surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
        else
          ! Otherwise, all flood water must go to runoff
          runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)
        endif
      endif

      runoffPipes=PipeCapacity   !Pipes are at their max capacity

  endif   !If runoff exceed pipe capacity

 end subroutine updateFlood
