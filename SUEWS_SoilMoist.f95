SUBROUTINE soilMoist_update(&
     nsurf,ConifSurf,DecidSurf,GrassSurf,&!input
     NonWaterFraction,&
     soilstoreCap,sfr,soilmoist,&
     soilmoistCap,soilstate,&!output
     vsmd,smd)
  IMPLICIT NONE

  INTEGER,INTENT(in)::nsurf,ConifSurf,DecidSurf,GrassSurf
  REAL(KIND(1d0)),INTENT(in)::NonWaterFraction
  REAL(KIND(1d0)),INTENT(in),DIMENSION(nsurf)::soilstoreCap,sfr,soilmoist

  REAL(KIND(1d0)),INTENT(out)::soilmoistCap,soilstate
  REAL(KIND(1d0)),INTENT(out)::vsmd,smd

  INTEGER :: is


  SoilMoistCap=0   !Maximum capacity of soil store [mm] for whole surface
  soilstate=0      !Area-averaged soil moisture [mm] for whole surface

  IF (NonWaterFraction/=0) THEN !Soil states only calculated if soil exists. LJ June 2017
     DO is=1,nsurf-1   !No water body included
        soilmoistCap=soilMoistCap+(soilstoreCap(is)*sfr(is)/NonWaterFraction)
        soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
     ENDDO
  ENDIF

  !If loop removed HCW 26 Feb 2015
  !if (ir==1) then  !Calculate initial smd
  smd=soilmoistCap-soilstate
  !endif

  ! Calculate soil moisture for vegetated surfaces only (for use in surface conductance)
  vsmd=0
  DO is=ConifSurf,GrassSurf  !Vegetated surfaces only
     IF ( sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) ==0 ) THEN
        vsmd=0
     ELSE
        vsmd=vsmd+(soilstoreCap(is) - soilmoist(is))*sfr(is)/(sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))
     END IF
     !write(*,*) is, vsmd, smd
  ENDDO


END SUBROUTINE soilMoist_update
