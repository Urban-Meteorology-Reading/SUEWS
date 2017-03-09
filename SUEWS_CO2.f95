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
    
  REAL(KIND(1d0)):: ActDorNorT, PopDorNorT   !Daytime, night-time or transition time
    
  REAL(KIND(1d0)):: QF_build, QF_traff, QF_metab
  
  !!!Move to input file later!!! Rename SUEWS_AnthropogenicHeat.txt to SUEWS_Anthropogenic.txt and add these there?
  REAL(KIND(1d0)):: FcEF_v_kgkm, EnEF_v_Jkm, EF_umolCO2perJ, FracFossilFuel, QF0_NB
  
  !AnthropCO2Method
  !2 - CO2 emissions based on QF calculated according to AnthropHeatMethod=2;
  !3 - CO2 emissions based on mean traffic rate and building energy use specified in SiteSelect 
  
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
  QF0_NB = 0.25
  
  ! Establish whether weekday or weekend ---------------------- 
  iu=1     !Set to 1=weekday
  IF(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) THEN  
     iu=2  !Set to 2=weekend
  ENDIF
  
  ! Calculate CO2 emissions from human metabolism -------------
  ! (Pop densities in ha-1 -> m-2)
  PopDorNorT = HumActivity_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu) !!! Set separately later !!!
  ActDorNorT = HumActivity_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)   !1=night, 2=day, 1-2=transition
  Fc_metab = (PopDensNighttime*(2-PopDorNorT) + PopDensDaytime*(PopDorNorT-1))/10000 * (120*(2-ActDorNorT) + 280*(ActDorNorT-1)) !umol m-2 s-1
  QF_metab = (PopDensNighttime*(2-PopDorNorT) + PopDensDaytime*(PopDorNorT-1))/10000 * (75*(2-ActDorNorT) + 175*(ActDorNorT-1)) !W m-2 
  
  ! Calculate CO2 emissions from traffic ----------------------
  IF(AnthropCO2Method == 2) THEN
     ! Assume fraction QF0_NB of temperature independent part of QF is traffic + metabolism
     QF_traff = QF0_NB*QF_SAHP_base - QF_metab   ! need to calculate QF_metab here again!!!
     IF(QF_traff < 0) THEN
         CALL ErrorHint(69,'QF metab exceeds base QF in Fc_anthro.',QF_metab,QF_SAHP_base,notUsedI)
         QF_traff = 0
     ENDIF
     Fc_traff = QF_traff / EnEF_v_Jkm * FcEF_v_kgkm*1e3*1e6/44   !Divide QF by energy emission factor and multiply by CO2 factor
  ELSEIF(AnthropCO2Method == 3) THEN
     ! Calculate using mean traffic rate [veh km cap-1 day-1] * emission factor [kg km-1] * 1e3 g kg-1 /44 g mol-1 * 1e6 umol mol-1
     ! Which popdens? !!!
      Fc_traff = PopDensNighttime/10000 * TrafficRate/(60*60*24) * FcEF_v_kgkm*1e3*1e6/44 * &
                   AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
  ENDIF
     
  ! Calculate CO2 emissions from building energy use ----------
  IF(AnthropCO2Method == 2) THEN
     ! Assume temperature independent part of QF is traffic + metabolism, 
     ! CDD part is electric A/C (no local CO2 emissions)
     ! HDD part is building energy use, split between electric (no local emissions CO2) and combustion (CO2) heating
     Fc_build = QF_SAHP_heat * EF_umolCO2perJ * FracFossilFuel
  ELSEIF(AnthropCO2Method == 3) THEN
     ! Calculate using building energy use [W cap-1]
     Fc_build = PopDensNighttime/10000 * BuildEnergyUse * EF_umolCO2perJ * FracFossilFuel * &
                 AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)   !Need to work out how to incorporate temp dependence !!!
     !Fc_build = Fc_build * 0.7 !Allen et al. (2010)
     !! Use daily mean temperature from previous day
     !IF(HDD(id-1,3) < 12.0) THEN      
     !   Fc_build = Fc_build * (1.0 + min(12.0-HDD(id-1,3),12.0--4)*0.05)   ! 5% increase in e (N.B. this is for heating only in Allen et al. 2010)!!!
     !   write(*,*) (1.0 + min(12.0-HDD(id-1,3),12.0--4)*0.05)
     !   ! Do not apply increase at high temp as this is AC (assume electric) N.B. slightly contradictory with FracFossilFuel
     !!ELSEIF(HDD(id-1,3) > 12.0) THEN
     !!   Fc_build = Fc_build * (1.0 + min(HDD(id-1,3)-12.0,35-12.0)*0.03)   ! 3% increase in e (N.B. this is for AC only in Allen et al. 2010)!!!
     !ENDIF                 
      
     ! Pigeon et al. (2007) (approx coeffs)
     IF(HDD(id-1,3) < 15.0) Fc_build = Fc_build * (1.0 + (15.0-HDD(id-1,3))*0.5)
     !write(*,*) (1.0 + (15.0-HDD(id-1,3))*0.5)
                 
  ENDIF   
    
  ! Sum components to give anthropogenic CO2 flux
  Fc_anthro = Fc_metab + Fc_traff + Fc_build
  
  !write(*,*) Fc_anthro, Fc_metab, Fc_traff, Fc_build
  !write(*,*) QF_SAHP_heat+QF_SAHP_base, QF_metab, QF_traff, QF_SAHP_heat
  
  RETURN
  
 ENDSUBROUTINE CO2_anthro
!======================================================================================== 