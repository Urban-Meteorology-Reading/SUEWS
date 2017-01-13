!========================================================================================
 SUBROUTINE CO2_biogen
! Created by HCW Aug 2016 to calculate biogenic component of CO2 flux.
! -ve Fc for uptake; +ve Fc for emission      
!
! To Do:
!   - Add different options
!========================================================================================

  USE AllocateArray   
  USE data_in
  USE PhysConstants   
  USE time
  
  IMPLICIT NONE
  
  INTEGER:: iv ! counter
  INTEGER:: BiogenCO2Choice = 1  !Move to RunControl later 
  
  REAL(KIND(1d0)):: PAR_umolm2s1
  REAL(KIND(1d0)):: Schmid2000_Pho, Flanag2002_Pho, Schmid2000_Res
    
  REAL(KIND(1d0)):: KdnToPAR ! Conversion from Kdn to PAR 
  REAL(KIND(1d0)):: F02_AMax, F02_alpha ! Coefficients for Flanagan et al. (2002) model
  
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: active_veg_fr
  
  ! Define coefficients (could be moved elsewhere) ------------
  KdnToPAR = 0.473  !Papaioannou et al. (1993) (mean annual value)
  F02_AMax  = -16.3     !Flanagan et al. (2002) Mean of summer 1998, 1999, 2000
  F02_alpha = -0.0205   !Flanagan et al. (2002) Mean of summer 1998, 1999, 2000
  
  ! Calculate PAR from Kdown ----------------------------------
  PAR_umolm2s1 = JtoumolPAR * KdnToPAR * avKdn
  
  ! Get active vegetation surface -----------------------------  
  DO iv=ivConif,ivGrass   !For vegetated surfaces
     !active_veg_fr(iv) = (sfr(iv+2)*(1-snowFrac(iv+2)))*(lai(id-1,iv)-LaiMin(iv))/(LaiMax(iv)-LaiMin(iv)) !LJ - snow?
     active_veg_fr(iv) = (sfr(iv+2))*(lai(id-1,iv)-LaiMin(iv))/(LaiMax(iv)-LaiMin(iv)) !LJ - snow?
  ENDDO  
   
  ! Calculate carbon uptake due to photosynthesis -------------
  ! Eq 6 Schmid et al. (2000) empirical model for hardwood forest (GEP)
  Schmid2000_Pho = -(-1.4 + 35*(PAR_umolm2s1/(590 + PAR_umolm2s1)))  !umol m-2 s-1
  IF(avKdn <= 0 .or. Schmid2000_Pho > 0) THEN   !uptake in daytime only, and remove any 
    Schmid2000_Pho = 0  
  ENDIF
  ! Flanagan et al. (2002) empirical model for temperate grassland (GPP)
  Flanag2002_Pho = F02_AMax*F02_alpha*PAR_umolm2s1/(F02_alpha*PAR_umolm2s1 + F02_AMax)  !umol m-2 s-1
  IF(avKdn <= 0 .or. Flanag2002_Pho > 0) THEN   !uptake in daytime only, and remove any 
    Flanag2002_Pho = 0
  ENDIF
    
  ! Calculate ecosystem respiration (for natural surfaces) ----
  ! Eq 5 Schmid et al. (2000)
  Schmid2000_Res = 1.08*exp(0.064*Temp_C)   !umol m-2 s-1   !!!Switch to using soil temp?
  
  IF(BiogenCO2Choice == 1) THEN   ! Use Schmid2000 for trees; Flanag2002 for grass
     Fc_photo = Schmid2000_Pho*active_veg_fr(ConifSurf-2)+ & 
             Schmid2000_Pho*active_veg_fr(DecidSurf-2)+ & 
             Flanag2002_Pho*active_veg_fr(GrassSurf-2)             
     Fc_respi = Schmid2000_Res*(sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))        
  ENDIF
  
  ! Combine to find biogenic CO2 flux
  Fc_biogen = Fc_photo + Fc_respi
  
  RETURN
  
 ENDSUBROUTINE CO2_biogen
!========================================================================================

 
!========================================================================================
 SUBROUTINE CO2_anthro(id, ih, imin)
! Created by HCW Aug 2016 to calculate anthropogenic component of CO2 flux.
! -ve Fc for uptake; +ve Fc for emission      
!
! To Do:
!   - Add specific diurnal profiles for traffic and building energy use for AnthropCO2Method=3
!      (currently the same anthropogenic heat profile is applied to both traffic and buildings)
!   - Add anthropogenic latent heat too?
!========================================================================================

  USE AllocateArray   
  USE data_in
  USE DefaultNotUsed
  USE sues_data
    
  IMPLICIT NONE
  
  INTEGER:: iu ! 1=weekday or 2=weekend
  INTEGER:: ih ! hour accounting for DLS
  INTEGER:: id, imin
    
  REAL(KIND(1d0)):: QF_metab, QF_traff, QF_build    !W m-2
  REAL(KIND(1d0)):: DorNorT   !Daytime, night-time or transition time
    
  !!!Move to input file later!!! Rename SUEWS_AnthropogenicHeat.txt to SUEWS_Anthropogenic.txt and add these there?
  REAL(KIND(1d0)):: FcEF_v_kgkm, EnEF_v_Jkm, EF_umolCO2perJ, FracFossilFuel
  
  !AnthropCO2Method
  !1 - CO2 emissions based on QF calculated according to AnthropHeatMethod=2;
  !2 - CO2 emissions based on mean traffic rate and building energy use specified in SiteSelect 
  
  ! Define coefficients ---------------------------------------   ! Move to inputs?
  ! CO2 emission factors
  FcEF_v_kgkm = 0.2069  !Average car, UK (DECC) [kg km-1]   
  !FcEF_v_kgkm = 0.295   !Moriwaki & Kanda (2004)  [kg km-1]
  ! Energy emission factors
  EnEF_v_Jkm = 3.97e6  ! [J km−1] Sailor & Lu (2004)
  ! Emission factors for fuels
  EF_umolCO2perJ =  51.0/44   !Natural gas (Moriwaki & Kanda 2004)
  !EF_umolCO2perJ =  67.6/44   !Kerosene (Moriwaki & Kanda 2004)
  FracFossilFuel = 0.70   !Proportion of building energy use from fossil fuels rather than electricity (0.6-0.8 for Sw)
                          !Some variation with season
  
  ! Establish whether weekday or weekend ---------------------- 
  iu=1     !Set to 1=weekday
  IF(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) THEN  
     iu=2  !Set to 2=weekend
  ENDIF
  
  ! Calculate CO2 emissions from human metabolism -------------
  ! (Pop densities in ha-1 -> m-2)
  DorNorT = CO2m_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)   !1=night, 2=day, 1-2=transition
  IF(iu == 1) THEN !If weekday, switch day/night populations
     Fc_metab = 120*PopDensNighttime/10000*(2-DorNorT) + 280*PopDensDaytime/10000*(DorNorT-1) !umol m-2 s-1
     QF_metab = 75*PopDensNighttime/10000*(2-DorNorT) + 175*PopDensDaytime/10000*(DorNorT-1) !W m-2   from Sailor & Lu (2004)
  ELSEIF(iu == 2) THEN  !If weekend, use night population always
     Fc_metab = 120*PopDensNighttime/10000*(2-DorNorT) + 280*PopDensNighttime/10000*(DorNorT-1) !umol m-2 s-1
     QF_metab = 75*PopDensNighttime/10000*(2-DorNorT) + 175*PopDensNighttime/10000*(DorNorT-1) !W m-2
  ENDIF
    
  ! Calculate CO2 emissions from traffic ----------------------
  IF(AnthropCO2Method == 2) THEN
     ! Assume temperature independent part of QF is traffic + metabolism
     QF_traff = QF_SAHP_base - QF_metab
     IF(QF_traff < 0) THEN
         CALL ErrorHint(69,'QF metab exceeds base QF in Fc_anthro.',QF_metab,QF_SAHP_base,notUsedI)
         QF_traff = 0
     ENDIF
     Fc_traff = QF_traff / EnEF_v_Jkm * FcEF_v_kgkm*1e3*1e6/44   !Divide QF by energy emission factor and multiply by CO2 factor
  ELSEIF(AnthropCO2Method == 3) THEN
     ! Calculate using mean traffic rate [veh km m-2 s-1] * emission factor [kg km-1] * 1e3 g kg-1 /44 g mol-1 * 1e6 umol mol-1
     Fc_traff = TrafficRate * FcEF_v_kgkm*1e3*1e6/44 * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
     QF_traff = TrafficRate * EnEF_v_Jkm * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)  !Also calculate QF 
  ENDIF
     
  ! Calculate CO2 emissions from building energy use ----------
  IF(AnthropCO2Method == 2) THEN
     ! Assume temperature independent part of QF is traffic + metabolism, 
     ! CDD part is electric A/C (no local CO2 emissions)
     ! HDD part is building energy use, split between electric (no local emissions CO2) and combustion (CO2) heating
     Fc_build = QF_SAHP_heat * EF_umolCO2perJ * FracFossilFuel
  ELSEIF(AnthropCO2Method == 3) THEN
     ! Calculate using building energy use [W m-2]
     Fc_build = BuildEnergyUse * EF_umolCO2perJ * FracFossilFuel * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
     QF_build = BuildEnergyUse * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)  !Also calculate QF 
  ENDIF   
    
  ! Combine to find anthropogenic CO2 flux
  Fc_anthro = Fc_metab + Fc_traff + Fc_build
  
  IF(AnthropHeatMethod == 3) THEN !!! N.B. need to implement this QF in QS according to OHMIncQF !!!
     QF = QF_metab + QF_traff + QF_build 
  ENDIF
  
  
  RETURN
  
 ENDSUBROUTINE CO2_anthro
!======================================================================================== 