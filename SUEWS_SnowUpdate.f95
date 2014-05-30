!In this subroutine the snow properties for the next timestep are updated

subroutine SnowUpdate
  use snowMod
  use data_in
  use sues_data
  use allocateArray
  
  implicit none

  real (Kind(1d0)):: tau_1=24*60*60!Function


  !Calculation of snow alebdo aging in an hourly timestep 
  !Lemonsu et al. 2010 (org: Verseghy (1991)&Baker et al. 1990)
  if (sum(SnowPack)>0) then 

      if (Temp_C<0) then
      	alb_snow = alb_snow-tau_a*(60*60)/tau_1
  	  else        
      	alb_snow = (alb_snow-albSnowMin)*exp(-tau_f*(60*60)/tau_1)+albSnowMin 
      endif	
    
      if (alb_snow<albSnowMin) alb_snow=albSnowMin 

  else
     alb_snow = 0
  endif

  !Update snow density
  do is=1,nsurf
     if (snowpack(is)>0) densSnow(is) = (densSnow(is)-densSnowMax)*exp(-tau_r*(60*60)/tau_1)+densSnowMax
     if (densSnow(is)>densSnowMax) densSnow(is)=densSnowMax
  enddo
  
end subroutine SnowUpdate

!====================================================================================
!====================================================================================
  