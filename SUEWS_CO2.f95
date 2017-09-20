!========================================================================================
 SUBROUTINE CO2_biogen
! Created by HCW Aug 2016 to calculate biogenic component of CO2 flux.
! -ve Fc for uptake; +ve Fc for emission
!
! Last modified:
! HCW 11 Apr 2017 - Tidied and merged with LJ code
! LJ   6 Apr 2017 - Minimum limit for soil respiration with BiogenCO2Choice = 2 was set to 0.6 umol m-2 s-1
!                 - Choice for non-rectancular hyperbola to calculate the biogenic CO2 flux added (BiogenCO2Choice = 2)
!                    (Bellucco et al. 2017, Agric. Forest. Met. 236, 113-122).
!                 - Both the "local Helsinki model" (BiogenCO2Choice = 2)) and "general model" BiogenCO2Choice = 3) are implemented
!                 - Snow fraction added to the calculation of active vegetation fraction and the soil respiration
!
! To Do:
!  - Active vegetation goes to zero with LAI minimum, but this needs to be changed so some minimum value
!    especially in the case of evergreentrees
!  - Move some of the parameters to input files
!  - Now on weekend nighttime population is used throughout the day. Do we need extra column in SiteSelect for daytime weekend population?
!========================================================================================

  USE AllocateArray
  USE data_in
  USE PhysConstants
  USE time

  IMPLICIT NONE

  INTEGER:: iv ! counter
  INTEGER:: BiogenCO2Choice = 1  !Move to RunControl later
  ! 1 - Rectangular hyperbola (Ruimy, Schmid, Flanagan)
  ! 2 - Non-rectangular hyperbola, Helsinki (Bellucco et al. 2016)
  ! 3 - Non-rectangular hyperbola, general  (Bellucco et al. 2016)

  REAL(KIND(1d0)):: PAR_umolm2s1

  REAL(KIND(1d0)):: KdnToPAR ! Conversion from Kdn to PAR

  REAL(KIND(1d0)):: Min_respi ! Minimum soil respiration rate (for cold-temperature limit)

  ! Coefficients for rectangular hyperbola light response curve (for each surface)
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: RecHyp_AMax    ! umol CO2 m-2 s-1
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: RecHyp_alpha   ! umol CO2 (umol photons)-1

  ! Coefficients for non-rectangular hyperbola (Bellucco et al. 2016)
  REAL(KIND(1d0)):: alpha_NRH, beta_NRH, theta_NRH

  REAL(KIND(1d0)):: Bellucco2016_Pho, Bellucco2016_Res

  REAL(KIND(1d0)),DIMENSION(nvegsurf):: active_veg_fr
  REAL(KIND(1d0)),DIMENSION(nvegsurf):: Fc_photo_surf   ! Photosynthesis for each vegetated surface

  ! Define coefficients (could be moved elsewhere) ------------
  !!! Tidy these up (nomenclature and signs)
  !KdnToPAR = 0.473  !Papaioannou et al. (1993) (mean annual value)
  KdnToPAR  = 0.46   !From Tsubo and Walker, 2005: PAR is on average 0.46 Kdown

  Min_respi = 0.6   ! LJ & UHEL (0.6 umol m-2 s-1 estimated from Hyytiala forest site)

  ! Calculate PAR from Kdown ----------------------------------
  PAR_umolm2s1 = JtoumolPAR * KdnToPAR * avKdn

  ! Calculate active vegetation surface -----------------------
  !Now this is zero always when LAI in its minimum values. This needs to vary between the summer time maximum and some minimum value
  !especially in the case of evergreen trees (i.e. in early March LAI can be in its minimum value, but air temperature and radiation
  !such that uptake can take place)
  DO iv=ivConif,ivGrass   !For vegetated surfaces. Snow included although quite often LAI will be in its minimum when snow on ground
     active_veg_fr(iv) = (sfr(iv+2)*(1-snowFrac(iv+2)))*(LAI(id-1,iv)-LAIMin(iv))/(LAIMax(iv)-LAIMin(iv))
  ENDDO


  IF(BiogenCO2Choice == 1) THEN   ! Rectangular hyperbola
     ! Set these in input files
     !! Schmid et al. (2000) empirical model for hardwood forest (Eq 6. with dark respiration rate removed)
     !RecHyp_AMax(:)  = 35
     !RecHyp_alpha(:) = 0.0593
     !! Flanagan et al. (2002) empirical model for temperate grassland (mean of summer 1998, 1999, 2000)
     !RecHyp_AMax(:)  = 16.3
     !RecHyp_alpha(:) = 0.0205
     ! Ruimy et al. (1995) review for plant canopies
     RecHyp_AMax(:)  = 43.35
     RecHyp_alpha(:) = 0.044

     ! Calculate carbon uptake due to photosynthesis -------------
     Fc_photo = 0
     DO iv=ivConif,ivGrass
        Fc_photo_surf(iv) = -RecHyp_AMax(iv)*RecHyp_alpha(iv)*PAR_umolm2s1/(RecHyp_alpha(iv)*PAR_umolm2s1 + RecHyp_AMax(iv))
        ! For active vegetation fraction only
        Fc_photo = Fc_photo + Fc_photo_surf(iv)*active_veg_fr(iv)  !umol m-2 s-1
     ENDDO

     ! Calculate ecosystem respiration ---------------------------
     ! Eq 5 Schmid et al. (2000), with minimum value from UHEL
     Fc_respi = MAX(Min_respi, 1.08*exp(0.064*Temp_C))   !umol m-2 s-1   !!!Switch to using soil temp?
     ! For natural surfaces only (with no snow)
     Fc_respi = Fc_respi* (sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+ &
                                sfr(GrassSurf)*(1-snowFrac(GrassSurf))+sfr(BSoilSurf)*(1-snowFrac(BSoilSurf)))

  ELSEIF(BiogenCO2Choice == 2 .OR. BiogenCO2Choice == 3) THEN   !Bellucco et al. (2016)
     IF(BiogenCO2Choice == 2) THEN   ! Local Helsinki model
        alpha_NRH = 0.031  ! umol CO2 umol photons-1
        beta_NRH  = 17.793 ! umol m^-2 s^-1
        theta_NRH = 0.723  ! -

        ! Respiration calculated from Helsinki (Järvi et al. 2012), 0.6 limit estimated from Hyytiala forest site
        Bellucco2016_Res = MAX(Min_respi, 3.229*exp(0.0329*Temp_C))

     ELSEIF (BiogenCO2Choice == 3) THEN ! General model
        ! Not currently recommended as includes also some anthropogenic impacts. Should maybe be separate from other biogen choices?

        ! Alpha and beta calculated as a function of vegetation cover fraction
        alpha_NRH = 0.005 + 0.016*(sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))   !umol CO2 umol photons-1
        beta_NRH  = -8.474 + 33.454*(sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf)) !umol m^-2 s^-1

        theta_NRH = 0.96   ! -

        ! Respiration set as fixed value
        Bellucco2016_Res = 2.43   ! umol m^-2 s^-1

     ENDIF

     ! Photosynthesis calculated using NRH
     Bellucco2016_Pho = -(1/(2*theta_NRH)*(alpha_NRH*PAR_umolm2s1+beta_NRH- &
                           sqrt((alpha_NRH*PAR_umolm2s1+beta_NRH)**2-4*alpha_NRH*beta_NRH*theta_NRH*PAR_umolm2s1)))

     Fc_photo = Bellucco2016_Pho*active_veg_fr(ConifSurf-2)+ &
                Bellucco2016_Pho*active_veg_fr(DecidSurf-2)+ &
                Bellucco2016_Pho*active_veg_fr(GrassSurf-2)

     Fc_respi = Bellucco2016_Res * (sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+ &
                                    sfr(GrassSurf)*(1-snowFrac(GrassSurf))+sfr(BSoilSurf)*(1-snowFrac(BSoilSurf)))

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
  IF (iu==2) PopDorNorT = 1   !!!
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
     ! Assume temperature independent part of QF is traffic + metabolism, !!! need QB_0_NB here!!!
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
