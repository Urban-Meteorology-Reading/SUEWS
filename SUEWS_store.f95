subroutine soilstore(StorCap)
!Rewritten mostly by LJ in 2010
!Modified by LJ in November 2012: 
!P>10 was not taken into account for impervious surfaces - Was fixed.
!Above impervious surfaces possibility of the state to exceed max capacity was limited
!although this should be possible - was fixed
!Calculation of storage change updated in this subroutine
!------------------------------------------------------------
  
  use data_in
  use SUES_data
  use gis_data
  use time
  use allocateArray
  
  implicit none
  real (Kind(1d0))::StorCap,EvPart

    
  !====If surface is irrigated grass, add external water use========
  if(is==GrassISurf) then  
     p=pin+ext_wu    
  elseif(is==ConifSurf.or.is==DecidSurf) then
     p=pin+wuhTrees/nsh 
  else
     p=pin
  endif
  !=================================================================
  !Add water from other surfaces at this point
  p=p+AddWater(is) 

  ! Reset runoff ffrom previous time interval
  runoff(is)=0.0
  runoffsoil(is)=0.0
    
  EvPart=0 !Extra evaporation if no water above impervious surfaces

  if(is==PavSurf.or.is==BldgSurf) then  !Impervious surfaces (paved, buildings)
   
     !Change water in the surface store 
     !If P too large!
     if (P>10) then
        runoff(is)=runoff(is)+(P-10)
        chang(is)=10-(drain(is)+ev)
     else
        chang(is)=P-(drain(is)+ev)
     endif
     
     state(is)=state(is)+chang(is)
     
     !Add water from neighbouring grids.
     !check - this is Bldg Surf originally - but I think should be paved (sg has changed)
     if (is==PavSurf) state(is)=state(is)+addImpervious/nsh
     
     runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is)!Drainage (not flowing to other surfaces) goes to runoff 
     
     if(state(is)<0.0) then  !Surface state cannot be negative
        SurPlus_Evap(is)=abs(state(is)) !This correspond addition to ev_per_interval
        state(is)=0.0
     !elseif (state(is)>StorCap) then !If state is higher than the max capacity??
     !   runoff(is)=runoff(is)+(state(is)-StorCap)
     !   state(is)=StorCap
  	 endif
       
 elseif(is>=3) then ! Pervious surfaces (conif, decid, grass unirr, grass irr, water body)
     
     !Add water to the surface store (is) 
     if  ((VegFraction+sfr(WaterSurf))/=0) then
     	EvPart=(SurPlus_evap(PavSurf)+SurPlus_evap(BldgSurf))*(sfr(is)/(VegFraction+sfr(WaterSurf)))
     else
        EvPart=0
     endif
     
	 ev=ev+EvPart
     
     if (is/=WaterSurf) then !Not for water body
        
        !Additional water input from other grids
        if (VegFraction/=0) then
        	P=P+addVeg/nsh*(sfr(is)/VegFraction)
        endif
        
        !Change in water stores !Modified by LJ in 4 Aug 2011.
        
        if (P>10) then !if 5min precipitation is larger than 10 mm
            runoff(is)=runoff(is)+(P-10)
            chang(is)=10-(drain(is)+ev)
		else
        	chang(is)=(P)-(drain(is)+ev)	
	    endif
		state(is)=state(is)+chang(is)
        
        !Add water in soil store. 
        soilmoist(is)=soilmoist(is)+Drain(is)*AddWaterRunoff(is) !Drainage goes to soil storages
   
		!If state of the surface is negative, remove water from soilstore
		if(state(is)<0.0) then  
        
       		if ((soilmoist(is)+state(is))>=0) then !If water in soilstore, water is removed
     			soilmoist(is)=soilmoist(is)+state(is)
        		state(is)=0.0
        	else !If not water in the soilstore evaporation does not occur
          		chang(is)=chang(is)-state(is)
          		ev=ev+state(is)
          		state(is)=0.0
			endif
        
        !elseif (state(is)>StorCap) then !!In original SUES, this was commented out
       !								!!Should be tested how change of S affects!!
       !    runoff(is)=runoff(is)+(state(is)-StorCap)
       !	   state(is)=StorCap
  	 	endif !state is negative

        !If soilstorage is full at this point, excess will go to surface runoff
     	if(soilmoist(is)>soilstoreCap(is))then  !Runoff occurs only through soilstorages
        	runoff(is)=runoff(is)+(soilmoist(is)-soilstoreCap(is))
        	soilmoist(is)=soilstoreCap(is)
     	elseif (soilmoist(is)<0) then
     		soilmoist(is)=0
     	endif
     !==============For water body===========================================
     elseif (is==WaterSurf) then
       
        !Add water from previous grid
   
        If(addWaterBody>0.and.sfr(watersurf)<0.000001)then
          print*,'subroutine soilstore'
          print*,'addwaterbody=',addwaterbody, 'but watersurf=',sfr(watersurf)
         stop
         endif   
           !If water body exists add water from surface
         if(sfr(watersurf)/=0)then
              ! check  what happens if no water body in this grid what would happen
              !  check whether this should be the maximum or minimm watersurf store (5, or (1
        !Do the changes in water body in certain order
             P=P+addWaterBody/nsh
     		chang(WaterSurf)=(P+FlowChange/nsh)-(ev)
			state(WaterSurf)=state(WaterSurf)+chang(WaterSurf)

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
        endif
     endif

     
 endif!Different surfaces

 
 !========RUNOFF=======================
 if (is<WaterSurf)then   !No water body

 	!Add runoff to pipes 
 	runoffPipes=runoffPipes+(runoff(is)*sfr(is))
    
    !If pipe capacity is full, surface runoff occurs
 	if (runoffPipes>PipeCapacity) then
   		if (is==PavSurf.or.is==BldgSurf) then

            if (sfr(WaterSurf)>0.0000001) then
   				runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)*(1-RunoffToWater) !if pipes are full, water will stay above ground
        		surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
            else
                runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)
            endif
            runoffPipes=PipeCapacity										 !and pipes are in their max. capacity
            							 
   		elseif (is>=ConifSurf.and.is<=GrassUSurf) then
        
        	if (sfr(WaterSurf)>0.0000001) then
        		runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)*(1-RunoffToWater)       !if pipes are full, water will stay above ground
                surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
            else
                 runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)	
        	endif
			runoffPipes=PipeCapacity								 !and pipes are in their max. capacity
   		    
		endif
 	endif
 endif
   
 runoff_per_interval=runoff_per_interval+(runoff(is)*sfr(is)) !The total runoff from the area
 
end subroutine soilstore
