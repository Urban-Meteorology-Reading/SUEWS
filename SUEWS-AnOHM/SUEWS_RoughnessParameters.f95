SUBROUTINE SUEWS_cal_RoughnessParameters(&
     RoughLenMomMethod,&! input:
     sfr,&! surface fractions
     bldgH,&
     EveTreeH,&
     DecTreeH,&
     porosity_id,&
     FAIBldg,FAIEveTree,FAIDecTree,Z,z0,zd,&
     planF,&! output:
     Zh,Z0m,Zdm,ZZD)
  ! Get surface covers and frontal area fractions (LJ 11/2010)
  ! Last modified:
  ! TS  18 Sep 2017 - added explicit interface
  ! HCW 08 Feb 2017 - fixed bug in Zh between grids, added default Z0m, Zdm
  ! HCW 03 Mar 2015
  ! sg feb 2012 - made separate subroutine
  !--------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,PARAMETER:: nsurf     = 7 ! number of surface types
  INTEGER,PARAMETER:: PavSurf   = 1 !When all surfaces considered together (1-7)
  INTEGER,PARAMETER:: BldgSurf  = 2
  INTEGER,PARAMETER:: ConifSurf = 3
  INTEGER,PARAMETER:: DecidSurf = 4
  INTEGER,PARAMETER:: GrassSurf = 5 !New surface classes: Grass = 5th/7 surfaces
  INTEGER,PARAMETER:: BSoilSurf = 6 !New surface classes: Bare soil = 6th/7 surfaces
  INTEGER,PARAMETER:: WaterSurf = 7

  INTEGER, INTENT(in) ::RoughLenMomMethod

  REAL(KIND(1d0)), DIMENSION(nsurf),INTENT(in) ::sfr! surface fractions


  REAL(KIND(1d0)), INTENT(in) ::bldgH
  REAL(KIND(1d0)), INTENT(in) ::EveTreeH
  REAL(KIND(1d0)), INTENT(in) ::DecTreeH
  REAL(KIND(1d0)), INTENT(in) ::porosity_id
  REAL(KIND(1d0)), INTENT(in) ::FAIBldg,FAIEveTree,FAIDecTree
  REAL(KIND(1d0)), INTENT(in) ::Z,Z0,Zd


  REAL(KIND(1d0)), INTENT(out) ::planF
  REAL(KIND(1d0)), INTENT(out) ::Zh
  REAL(KIND(1d0)), INTENT(out) ::Z0m
  REAL(KIND(1d0)), INTENT(out) ::Zdm
  REAL(KIND(1d0)), INTENT(out) ::ZZD



  REAL(KIND(1d0)) ::areaZh
  INTEGER, PARAMETER :: notUsedI=-55
  REAL(KIND(1d0)),PARAMETER:: notUsed=-55.5
  REAL(KIND(1D0)):: z0m4Paved,z0m4Grass,z0m4BSoil,z0m4Water   !Default values for roughness lengths [m]

  areaZh =(sfr(BldgSurf)+sfr(ConifSurf)+sfr(DecidSurf)) !Total area of buildings and trees

  ! Set default values (using Moene & van Dam 2013, Atmos-Veg-Soil Interactions, Table 3.3)
  Z0m4Paved = 0.003 !estimate
  Z0m4Grass = 0.02
  Z0m4BSoil = 0.002
  Z0m4Water = 0.0005
  ! test initialisation
  Z0m=z0*1.
  ZZD=0.
  Zdm=zd*1.

  !------------------------------------------------------------------------------
  !If total area of buildings and trees is larger than zero, use tree heights and building heights to calculate zH
  IF (areaZh/=0) THEN
     Zh=bldgH*sfr(BldgSurf)/areaZh + EveTreeH*sfr(ConifSurf)/areaZh + DecTreeH*(1-porosity_id)*sfr(DecidSurf)/areaZh
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
           !planF=FAIBldg*sfr(BldgSurf)/areaZh+FAItree*sfr(ConifSurf)/areaZh+FAItree*(1-porosity_id)*sfr(DecidSurf)/areaZh
           planF=FAIBldg*sfr(BldgSurf)/areaZh + FAIEveTree*sfr(ConifSurf)/areaZh + FAIDecTree*(1-porosity_id)*sfr(DecidSurf)/areaZh
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
        Z0m = (z0m4Paved*sfr(PavSurf) + z0m4Grass*sfr(GrassSurf) + z0m4BSoil*sfr(BSoilSurf) + z0m4Water*sfr(WaterSurf))/(1-areaZh)
        Zdm = 0
        CALL ErrorHint(15,'Setting Z0m and Zdm using default values',Z0m,Zdm,notUsedI)
     ELSEIF(areaZh==1)THEN  !If, for some reason, Zh = 0 and areaZh == 1, assume height of 10 m and use rule-of-thumb
        Z0m = 1
        Zdm = 7
        CALL ErrorHint(15,'Assuming mean height = 10 m, Setting Z0m and Zdm to default value',Z0m,Zdm,notUsedI)
     ENDIF
  ENDIF


  ZZD=Z-Zdm

  IF(RoughLenMomMethod==1) THEN  !use z0, zd values provided in input file
     Z0m=z0
     Zdm=zd
     ZZD=Z-zd
  ENDIF

  ! Error messages if aerodynamic parameters negative
  IF(Z0m<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, z0 < 0 m.',Z0m,notUsed,notUsedI)
  IF(Zdm<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, zd < 0 m.',Zdm,notUsed,notUsedI)
  IF(zzd<0) CALL ErrorHint(14,'In SUEWS_RoughnessParameters.f95, (z-zd) < 0 m.',zzd,notUsed,notUsedI)
END SUBROUTINE SUEWS_cal_RoughnessParameters
