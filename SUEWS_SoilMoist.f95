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


!========== Calculate soil moisture ============
SUBROUTINE SUEWS_cal_SoilMoist(&
  nsurf,&
  NonWaterFraction,SoilMoistCap,SoilStoreCap,surf_chang_per_tstep,&
  soilmoist,soilmoistOld,sfr,smd,smd_nsurf,tot_chang_per_tstep,&
  soilstate)

  IMPLICIT NONE

  INTEGER,INTENT(in) ::nsurf
  REAL(KIND(1d0)),INTENT(in)::NonWaterFraction
  REAL(KIND(1d0)),INTENT(in)::SoilMoistCap
  REAL(KIND(1d0)),INTENT(in)::surf_chang_per_tstep
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::soilmoist
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::soilmoistOld
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::SoilStoreCap        !Capacity of soil store for each surface [mm]

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::smd_nsurf
  
  REAL(KIND(1d0)),INTENT(out)::soilstate
  REAL(KIND(1d0)),INTENT(out)::smd
  REAL(KIND(1d0)),INTENT(out)::tot_chang_per_tstep

  REAL(KIND(1d0)),PARAMETER::NotUsed=-999
  INTEGER :: is

  soilstate=0       !Area-averaged soil moisture [mm] for whole surface
  IF (NonWaterFraction/=0) THEN !Fixed for water surfaces only
     DO is=1,nsurf-1   !No water body included
        soilstate=soilstate+(soilmoist(is)*sfr(is)/NonWaterFraction)
        IF (soilstate<0) THEN
           CALL ErrorHint(62,'SUEWS_Calculations: total soilstate < 0 (just added surface is) ',soilstate,NotUsed,is)
        ELSEIF (soilstate>SoilMoistCap) THEN
           CALL ErrorHint(62,'SUEWS_Calculations: total soilstate > capacity (just added surface is) ',soilstate,NotUsed,is)
           !SoilMoist_state=soilMoistCap !What is this LJ 10/2010 - SM exceeds capacity, but where does extra go?HCW 11/2014
        ENDIF
     ENDDO  !end loop over surfaces
  ENDIF

  ! Calculate soil moisture deficit
  smd=soilMoistCap-soilstate   !One value for whole surface
  smd_nsurf=SoilstoreCap-soilmoist   !smd for each surface

  ! Soil stores can change after horizontal water movements
  ! Calculate total change in surface and soil state
  tot_chang_per_tstep = surf_chang_per_tstep   !Change in surface state
  DO is=1,(nsurf-1)   !No soil for water surface (so change in soil moisture is zero)
     tot_chang_per_tstep = tot_chang_per_tstep + ((SoilMoist(is)-soilmoistOld(is))*sfr(is))   !Add change in soil state
  ENDDO

END SUBROUTINE SUEWS_cal_SoilMoist
