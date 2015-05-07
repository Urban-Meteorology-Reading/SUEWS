! Conversion of water use (irrigation) 
! Last modified by HCW 12 Feb 2015
!  - Water use [mm] now inidcates the amount of water supplied for each surface
! Last modified by HCW 26 Jan 2015
!  - Water use [mm] is the same for each surface at the moment and indicates the 
! amount of water supplied for each irrigated area
! To Do:
!	- Add functionality for water on paved surfaces (street cleaning, fountains)
!===================================================================================
 subroutine WaterUse
  
  use allocateArray
  use data_in
  use defaultNotUsed
  use sues_data
  use time
  
  IMPLICIT NONE
  
  real(kind(1d0)):: wu   		!Water use for the model timestep [mm]
  real(kind(1d0)):: InternalWaterUse    !Internal water use for the model timestep [mm]
  real(kind(1d0)):: WuFr=1
  integer:: ih   !Hour corrected for Daylight savings
  integer:: iu   !1=weekday OR 2=weekend    	    
    
  ! --------------------------------------------------------------------------------
  ! If water used is observed and provided in the met forcing file, units are m3
  ! Divide observed water use (in m3) by water use area to find water use (in mm)
  if (WU_choice==1) then   !If water use is observed
     ! Calculate water use area [m2] for each surface type
     WUAreaEveTr_m2 = IrrFracConif*sfr(ConifSurf)*SurfaceArea
     WUAreaDecTr_m2 = IrrFracDecid*sfr(DecidSurf)*SurfaceArea   
     WUAreaEveTr_m2 = IrrFracGrass*sfr(GrassSurf)*SurfaceArea   
     WUAreaTotal_m2 = WUAreaEveTr_m2 + WUAreaDecTr_m2 + WUAreaGrass_m2  

     !Set water use [mm] for each surface type to zero initially
     wu_EveTr=0
     wu_DecTr=0
     wu_Grass=0  
     if(wu_m3==NAN.or.wu_m3==0) then   !If no water use
        wu_m3=0
        wu=wu_m3
     else 			       !If water use
        if (WUAreaTotal_m2>0) then
           wu = (wu_m3/WUAreaTotal_m2/1000)     !Water use in mm for the whole irrigated area
           if (WUAreaEveTr_m2>0) then
              wu_EveTr=wu   			!Water use for Irr EveTr in mm - these are all the same at the moment
              wu_EveTr=wu_EveTr*IrrFracConif	!Water use for EveTr in mm
           endif
           if (WUAreaDecTr_m2>0) then
	      wu_DecTr=wu   			!Water use for Irr DecTr in mm - these are all the same at the moment
	      wu_DecTr=wu_DecTr*IrrFracDecid	!Water use for DecTr in mm
           endif
           if (WUAreaGrass_m2>0) then
              wu_Grass=wu   			!Water use for Irr Grass in mm - these are all the same at the moment
              wu_Grass=wu_Grass*IrrFracGrass	!Water use for Grass in mm
           endif            
           wu = (wu_m3/SurfaceArea/1000)        !Water use for the whole study area in mm 
        endif
     endif   
      
  ! --------------------------------------------------------------------------------  
  ! If water use is modelled, calculate at timestep of model resolution [mm]
  elseif (WU_choice==0) then   !If water use is modelled
 
     ! Account for Daylight saving
     ih=it-DLS
     if (ih<0) ih=23

     ! Weekday or weekend profile 	 
     iu=1     !Set to 1=weekday
     if(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) then  
        iu=2  !Set to 2=weekend
     endif
     
     !write(*,*) (NSH*(ih+1-1)+imin*NSH/60+1)
        
     ! ---- Automatic irrigation ----
     wu_EveTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,2)   !Automatic evergreen trees
     wu_DecTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,5)   !Automatic deciduous trees
     wu_Grass = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,8)   !Automatic grass
     
     ! ---- Manual irrigation ----
     WuFr=1 !Initialize WuFr to 1, but if raining, reduce manual fraction of water use
     ! If cumulative daily precipitation exceeds 2 mm
     if(HDD(id,5)>2) then    !.and.WU_Day(id-1,3)>0) then !Commented out HCW 23/01/2015
        WuFr=0   ! 0 -> No manual irrigation if raining 
     endif
     
     ! Add manual to automatic to find total irrigation
     wu_EveTr = wu_EveTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,3)) !Manual evergreen trees
     wu_DecTr = wu_DecTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,6)) !Manual deciduous trees
     wu_Grass = wu_Grass + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day(id-1,9)) !Manual grass

     ! Added HCW 12 Feb 2015.
     wu_EveTr=wu_EveTr*IrrFracConif	!Water use for EveTr [mm]
     wu_DecTr=wu_DecTr*IrrFracDecid	!Water use for DecTr [mm]
     wu_Grass=wu_Grass*IrrFracGrass	!Water use for Grass [mm]
     
     ! Total water use for the whole study area [mm]
     wu = wu_EveTr*sfr(ConifSurf) + wu_DecTr*sfr(DecidSurf) + wu_Grass*sfr(GrassSurf)
  
  endif   !End WU_choice
  ! --------------------------------------------------------------------------------

  ! Internal water use is supplied in SUEWS_Irrigation in mm h-1
  ! Convert to mm for the model timestep
  InternalWaterUse = InternalWaterUse_h/nsh_real
  
  ! Remove InternalWaterUse from the total water use
  ext_wu = wu-(InternalWaterUse+OverUse)
  ! Check ext_wu cannot be negative
  if (ext_wu<0) then
     overUse=abs(ext_wu)
     ext_wu=0
  else
     OverUse=0
  endif
	
  int_wu = wu-ext_wu
  
  ! Decrease the water use for each surface by the same proportion
  if(ext_wu/=0.and.wu/=0) then
     wu_EveTr = wu_EveTr*ext_wu/wu
     wu_DecTr = wu_DecTr*ext_wu/wu
     wu_Grass = wu_Grass*ext_wu/wu
  endif 
  
endsubroutine WaterUse
!===================================================================================
