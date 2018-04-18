!In this subroutine the snow properties are updated. The aging functions work for hourly air
!temperature (dt=1). As the code timestep is typically smaller than one hour, the new albedo
!and density are calculated using the timestep air temperature and then these are scaled to timestep
!according to NSH.  Made by LJ in sprint 2014
!Last update:
! TS 17 Sep 2017 - Improve the explicit interface
! LJ 7 July 2015 - Changed to work with shorter timestep: defined by tstep. Cleaning of the code.
!
!========================================================================

SUBROUTINE SnowUpdate(&
     nsurf,tstep,&!input
     Temp_C_hr,&
     tau_a,&
     tau_f,&
     tau_r,&
     SnowDensMax,&
     SnowDensMin,&
     SnowAlbMin,&
     SnowPack,&
     SnowAlb,&!inout
     SnowDens)


  IMPLICIT NONE

  INTEGER,INTENT(in)::nsurf
  INTEGER,INTENT(in)::tstep

  REAL(KIND(1D0)),INTENT(in)::Temp_C_hr        !Air temperature
  REAL(KIND(1D0)),INTENT(in)::tau_a
  REAL(KIND(1D0)),INTENT(in)::tau_f
  REAL(KIND(1D0)),INTENT(in)::tau_r
  REAL(KIND(1D0)),INTENT(in)::SnowDensMax
  REAL(KIND(1D0)),INTENT(in)::SnowDensMin
  REAL(KIND(1D0)),INTENT(in)::SnowAlbMin

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::SnowPack

  REAL(KIND(1d0)),INTENT(inout)::SnowAlb

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowDens


  INTEGER::is
  REAL(KIND(1D0))::alb_change,&     !Change in snow albedo
       dens_change,&    !Change in snow density
       tau_1         !Number of seconds in a day

  !Initialize
  alb_change=0
  dens_change=0
  tau_1=24*60*60

  !==========================================================
  !Calculation of snow albedo by Lemonsu et al. 2010
  !(org: Verseghy (1991)&Baker et al.(1990))
  IF (SUM(SnowPack)>0) THEN !Check if snow on any of the surfaces
     IF (Temp_C_hr<0) THEN
        !alb_change = tau_a*(60*60)/tau_1
        alb_change = tau_a*(tstep)/tau_1
        SnowAlb = SnowAlb-alb_change
     ELSE
        !alb_change = exp(-tau_f*(60*60)/tau_1)
        alb_change = EXP(-tau_f*(tstep)/tau_1)
        SnowAlb = (SnowAlb-SnowAlbMin)*alb_change+SnowAlbMin
     ENDIF
     IF (SnowAlb<SnowAlbMin) SnowAlb=SnowAlbMin !Albedo cannot be smaller than the min albedo
  ELSE
     SnowAlb = 0
  ENDIF

  !Update snow density: There is a mistake in Järvi et al. (2014): tau_h should be tau_1
  DO is=1,nsurf

     !If SnowPack existing
     IF (SnowPack(is)>0) THEN
        dens_change = EXP(-tau_r*(tstep)/tau_1)
        IF (SnowPack(is)>0) SnowDens(is) = (SnowDens(is)-SnowDensMax)*dens_change+SnowDensMax
        IF (SnowDens(is)>SnowDensMax) SnowDens(is)=SnowDensMax
     ELSE
        SnowDens(is) = SnowDensMin
     ENDIF
  ENDDO

  ! write(*,*) SnowAlb

END SUBROUTINE SnowUpdate

!====================================================================================
!====================================================================================
