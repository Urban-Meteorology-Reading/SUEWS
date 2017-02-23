SUBROUTINE RoughnessParameters
  ! Get surface covers and frontal area fractions (LJ 11/2010)
  ! Last modified by HCW 08 Feb 2017 - fixed bug in Zh between grids, added default z0m, zdm
  !                  HCW 03 Mar 2015
  !                  sg feb 2012 - made separate subroutine
  !--------------------------------------------------------------------------------
  USE data_in
  USE gis_data
  USE sues_data
  USE allocateArray
  USE mod_z
  USE defaultNotUsed
  USE time

  IMPLICIT NONE

  REAL(KIND(1D0)):: z0m4Paved,z0m4Grass,z0m4BSoil,z0m4Water   !Default values for roughness lengths [m]
  
  ! Set default values (using Moene & van Dam 2013, Atmos-Veg-Soil Interactions, Table 3.3)
  Z0m4Paved = 0.003 !estimate 
  Z0m4Grass = 0.02
  Z0m4BSoil = 0.002
  Z0m4Water = 0.0005
  
  !------------------------------------------------------------------------------
  !If total area of buildings and trees is larger than zero, use tree heights and building heights to calculate zH
  IF (areaZh/=0) THEN
     Zh=bldgH*sfr(BldgSurf)/areaZh + EveTreeH*sfr(ConifSurf)/areaZh + DecTreeH*(1-porosity(id))*sfr(DecidSurf)/areaZh
  ELSE
     Zh=0   !Set Zh to zero if areaZh = 0
  ENDIF
  
  IF(Zh/=0)THEN
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
  ELSEIF(Zh==0)THEN   !If zh calculated to be zero, set default roughness length and displacement height
     IF(areaZh /= 0) CALL ErrorHint(15,'In SUEWS_RoughnessParameters.f95, zh = 0 m but areaZh > 0',zh,areaZh,notUsedI)
     !Estimate z0 and zd using default values and surfaces that do not contribute to areaZh
     IF(areaZh /= 1)THEN
        z0m = (z0m4Paved*sfr(PavSurf) + z0m4Grass*sfr(GrassSurf) + z0m4BSoil*sfr(BSoilSurf) + z0m4Water*sfr(WaterSurf))/(1-areaZh)
        zdm = 0
        CALL ErrorHint(15,'Setting z0m and zdm using default values',z0m,zdm,notUsedI)
     ELSEIF(areaZh==1)THEN  !If, for some reason, Zh = 0 and areaZh == 1, assume height of 10 m and use rule-of-thumb
        z0m = 1
        zdm = 7
        CALL ErrorHint(15,'Assuming mean height = 10 m, Setting z0m and zdm to default value',z0m,zdm,notUsedI)   
     ENDIF 
  ENDIF
     
  ZZD=Z-zdm

  ! Error messages if aerodynamic parameters negative
  IF(z0m<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, z0 < 0 m.',z0m,notUsed,notUsedI)
  IF(zdm<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, zd < 0 m.',zdm,notUsed,notUsedI)
  IF(zzd<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, (z-zd) < 0 m.',zzd,notUsed,notUsedI)
END SUBROUTINE RoughnessParameters
