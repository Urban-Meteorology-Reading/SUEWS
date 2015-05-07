!In this subroutine the snow properties are updated. The aging functions work for hourly air
!temperature (dt=1). As the code timestep is typically smaller than one hour, the new albedo
!and density are calculated using the timestep air temperature and then these are scaled to timestep
!according to NSH.  Made by LJ in sprint 2014
!Last update:
!  LJ in May 2015 - Changed to work with shorter timestep. Cleaning of the code made
!
!========================================================================

 subroutine SnowUpdate(Temp_C_hr)

  use allocateArray
  use snowMod
  use sues_data

  IMPLICIT NONE

  REAL(KIND(1D0))::alb_change,&     !Change in snow albedo
                   dens_change,&    !Change in snow density
                   tau_1,& !Number of seconds in a day
                   Temp_C_hr      !Air temperature

  !Initialize
  alb_change=0
  dens_change=0
  tau_1=24*60*60

  !Calculation of snow albedo by Lemonsu et al. 2010 (org: Verseghy (1991)&Baker et al.(1990))

  if (sum(SnowPack)>0) then !Check if snow on any of the surfaces

      if (Temp_C_hr<0) then
        alb_change = tau_a*(60*60)/tau_1
      	alb_snow = alb_snow-alb_change/NSH_real
  	  else
        alb_change = exp(-tau_f*(60*60)/tau_1)
      	alb_snow = (alb_snow-albSnowMin)*alb_change/NSH_real+albSnowMin
      endif	
    
      if (alb_snow<albSnowMin) alb_snow=albSnowMin !Albedo cannot be smaller than the min albedo

  else
     alb_snow = 0
  endif

  !Update snow density
  do is=1,nsurf
     dens_change = exp(-tau_r*(60*60)/tau_1)
     if (snowpack(is)>0) densSnow(is) = (densSnow(is)-densSnowMax)*dens_change/NSH_real+densSnowMax
     if (densSnow(is)>densSnowMax) densSnow(is)=densSnowMax
  enddo


  
 end subroutine SnowUpdate

!====================================================================================
!====================================================================================
  