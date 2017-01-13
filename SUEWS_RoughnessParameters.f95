SUBROUTINE RoughnessParameters
  ! Last modified by HCW 03 Mar 2015
  ! Get surface covers and frontal area fractions (LJ 11/2010)
  ! sg feb 2012 -- made separate subroutine
  !INPUT: id  Day index
  !--------------------------------------------------------------------------------
  USE data_in
  USE gis_data
  USE sues_data
  USE allocateArray
  USE mod_z
  USE defaultNotUsed
  USE time

  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !If total area of buildings and trees is larger than zero
  IF (areaZh/=0)THEN
     !Zh=bldgH*sfr(BldgSurf)/areaZh+treeH*sfr(ConifSurf)/areaZh+treeH*(1-porosity(id))*sfr(DecidSurf)/areaZh
     Zh=bldgH*sfr(BldgSurf)/areaZh + EveTreeH*sfr(ConifSurf)/areaZh + DecTreeH*(1-porosity(id))*sfr(DecidSurf)/areaZh
  ENDIF

  !Calculate Z0m and Zdm depending on the Z0 method
  IF(RoughLenMomMethod==2) THEN  !Rule of thumb (G&O 1999)
     Z0m=0.1*Zh
     Zdm=0.7*Zh
  ELSEIF(RoughLenMomMethod==3)THEN !MacDonald 1998
     IF (areaZh/=0)THEN  !Plan area fraction
        !planF=FAIBldg*sfr(BldgSurf)/areaZh+FAItree*sfr(ConifSurf)/areaZh+FAItree*(1-porosity(id))*sfr(DecidSurf)/areaZh
        planF=FAIBldg*sfr(BldgSurf)/areaZh + FAIEveTree*sfr(ConifSurf)/areaZh + FAIDecTree*(1-porosity(id))*sfr(DecidSurf)/areaZh
     ELSE
        planF=0.00001
        Zh=1
     ENDIF

     Zdm=(1+4.43**(-sfr(BldgSurf))*(sfr(BldgSurf)-1))*Zh
     Z0m=((1-Zdm/Zh)*EXP(-(0.5*1.0*1.2/0.4**2*(1-Zdm/Zh)*planF)**(-0.5)))*Zh
  ENDIF

  ZZD=Z-zdm

  ! Error messages if aerodynamic parameters negative
  IF(z0m<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, z0 < 0 m.',z0m,notUsed,notUsedI)
  IF(zzd<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, zd < 0 m.',zzd,notUsed,notUsedI)
END SUBROUTINE RoughnessParameters
