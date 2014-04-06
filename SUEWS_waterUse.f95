subroutine waterUse
    !Funtion calculates water use per time interval reducing the internal 
    !water use if it has not yet been reduced.
    !Modified by LJ (10 Sep 2010)
    !--------------------------------------------------------------
	use SUES_data
	use data_in
    use time
    use allocateArray
    use defaultNotUsed
    implicit none
   
    real (kind(1d0))::WuFr=1
   integer::ih ! hour corrected for daylightsavings
    !--------------------------------------------------------------------
    !If measured water use is used divide wuh [m3] with water use area  
    if (WU_choice==1) then!
      if (wuh==NAN) then
         wuh=0
         wuhTrees=wuh/(WaterUseAreaTrees*10)
      else 
         wuh=wuh/(WaterUseAreaGrass*10)
         if(WaterUseAreaTrees>0)then
         	wuhTrees=wuh/(WaterUseAreaTrees*10)
         else
            wuhTrees=0
         endif
      endif
    
    !Else if water use is modelled, calculate hourly wu
    elseif (WU_choice==0) then  
        !If no hourly water use restrictions exist
         
        ! times all standard time - adjust for DLS -day light saving (calculated in daily state) if DLS (substract 1 h)
        ih=it-DLS
        if(ih<0)then
           ih=23
        endif
        
        !automatic then manual
        wuh=HourWat(ih)*WU_Day(id-1,2)
        wuhTrees=HourWat(ih)*WU_Day(id-1,5)

        !! if actually raining -- reduce manual fraction of water use
 
        WUFr=1 !Initialize this
        if(HDD(id,5)>2.and.WU_Day(id-1,3)>0)then
 			WUFr=0  ! check what this should be 1= all applied
  		endif
        
        wuh=wuh+(WuFr*HourWat(ih)*WU_Day(id-1,3))
        wuhTrees=wuhTrees+(WuFr*HourWat(ih)*WU_Day(id-1,6))
  
    endif !End WU choice

     
  
	! remove internal water use from grass only for now. Need to be fixed!!
    
	ext_wuh=wuh-((InternalWaterUse/24)+OverUse) !Hourly value
	
    if(ext_wuh<0)then
  		overUse=abs(ext_wuh)
  		ext_wuh=0
	else
   		OverUse=0
	endif
	int_wuh=wuh-ext_wuh


	!convert to part of an hour
	if(Int_wuh>0) then
  		int_wu=Int_wuh/nsh
	else
   		int_wu=0
	endif
      
	if(ext_wuh>0) then
   		ext_wu=ext_wuh/nsh
	else
   		ext_wu=0
	endif

end subroutine waterUse

