!In this subroutine the snow properties are updated. The aging functions work for hourly air
!temperature (dt=1). As the code timestep is typically smaller than one hour, the new albedo
!and density are calculated using the timestep air temperature and then these are scaled to timestep
!according to NSH.  Made by LJ in sprint 2014
!Last update:
!  LJ in 7 July 2015 - Changed to work with shorter timestep: defined by tstep. Cleaning of the code.
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

  !==========================================================
  !Calculation of snow albedo by Lemonsu et al. 2010 
  !(org: Verseghy (1991)&Baker et al.(1990))
  if (sum(SnowPack)>0) then !Check if snow on any of the surfaces

      if (Temp_C_hr<0) then
        !alb_change = tau_a*(60*60)/tau_1
        alb_change = tau_a*(tstep)/tau_1
      	snowAlb = SnowAlb-alb_change
  	  else
        !alb_change = exp(-tau_f*(60*60)/tau_1)
        alb_change = exp(-tau_f*(tstep)/tau_1)
      	SnowAlb = (SnowAlb-SnowAlbMin)*alb_change+SnowAlbMin
      endif	

      if (SnowAlb<SnowAlbMin) SnowAlb=SnowAlbMin !Albedo cannot be smaller than the min albedo

  else
     SnowAlb = 0
  endif

  !Update snow density: There is a mistake in Järvi et al. (2014): tau_h should be tau_1
  do is=1,nsurf

    !If snowPack existing
    if (snowPack(is)>0) then
       dens_change = exp(-tau_r*(tstep)/tau_1)
       if (snowpack(is)>0) SnowDens(is) = (SnowDens(is)-SnowDensMax)*dens_change+SnowDensMax
       if (SnowDens(is)>SnowDensMax) SnowDens(is)=SnowDensMax
    else
       SnowDens(is) = SnowDensMin
    endif
  enddo

 end subroutine SnowUpdate

!====================================================================================
!====================================================================================
  