!========================================================================================
! Created by HCW Aug 2016 to calculate biogenic component of CO2 flux.
! This subroutine is still under development and in the equations there might be bugs and
! the code is not well commented.
!
! Last modified:
! MH 20 Jun 2017  - Tidied and renamed from SUEWS_CO2.f95 to SUEWS_CO2Biogen.f95
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

! EmissionsMethod:
!  11-16 - Rectangular hyperbola (Ruimy, Schmid, Flanagan)
!  21-26 - Non-rectangular hyperbola, Helsinki (Bellucco et al. 2017)
!  31-36 - Non-rectangular hyperbola, general  (Bellucco et al. 2017)
!========================================================================================
SUBROUTINE CO2_biogen(EmissionsMethod,&
     ivConif,ivDecid,ivGrass,ConifSurf,DecidSurf,GrassSurf,BSoilSurf,&
     snowFrac,nsurf,NVegSurf,avkdn,Temp_C,sfr,LAI_day,LAIMin,LAIMax,&
     alpha_bioCO2,beta_bioCO2,theta_bioCO2,alpha_enh_bioCO2,beta_enh_bioCO2,&
     resp_a,resp_b,min_res_bioCO2,Fc_biogen,Fc_respi,Fc_photo,&
     notUsed,notUsedI)

  IMPLICIT NONE
  INTEGER,INTENT(in):: EmissionsMethod
  INTEGER,INTENT(in)::&
       ivConif,ivDecid,ivGrass,ConifSurf,DecidSurf,GrassSurf,BSoilSurf,&
       nsurf,nvegSurf,&
       notUsedI
  REAL(KIND(1d0)),INTENT(in)::&
       avkdn,&
       Temp_C,&
       notUsed

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::&
       sfr,&   !Surface fractions [-]
       snowFrac
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in):: LAI_day
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::&
       LAIMin, LAIMax,&      ! [m2 m-2]
                                !  BiogenCO2Code,&       !Biogenic CO2 Code for SUEWS_BiogenCO2.txt
       alpha_bioCO2,&
       beta_bioCO2,&
       theta_bioCO2,&
       alpha_enh_bioCO2,&
       beta_enh_bioCO2,&
       resp_a,&
       resp_b,&
       min_res_bioCO2
  REAL(KIND(1D0)),INTENT(out):: &
       Fc_biogen,&
       Fc_respi,Fc_photo

  INTEGER:: iv ! counter

  REAL(KIND(1d0)):: &
       PAR_umolm2s1,&
       Bellucco2017_Pho,&     ! Photosynthesis (Bellucco et al. 2016)
       Bellucco2017_Res,&     ! Respiration (Bellucco et al. 2016)
       Bellucco2017_Res_surf,&! Respiration for each vegetated surface
       VegFracSum             ! Sum of vegetation fractions without water. Could be moved elsewhere later.

  REAL(KIND(1d0)),DIMENSION(nvegsurf)::&
       active_veg_fr,&         ! Active vegetation fraction
       Fc_photo_surf,&         ! Photosynthesis for each vegetated surface
       Bellucco2017_Pho_surf ! Photosynthesis for each vegetated surface

  REAL(KIND(1d0)),DIMENSION(nvegsurf)::&
       alpha_bioCO2_v2,&
       beta_bioCO2_v2,&
       theta_bioCO2_v2

  REAL(KIND(1d0)),PARAMETER :: &
       JtoumolPAR = 4.6,&
       KdntoPAR = 0.46

  !-----------------------------------------------------------------------

  ! Calculate PAR from Kdown ----------------------------------
  PAR_umolm2s1 = JtoumolPAR * KdntoPAR * avKdn

  VegFracSum = sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf)

  ! Calculate active vegetation surface -----------------------
  !Now this is zero always when LAI in its minimum values. This needs to vary between the summer time maximum and some minimum value
  !especially in the case of evergreen trees (i.e. in early March LAI can be in its minimum value, but air temperature and radiation
  !such that uptake can take place)
  DO iv=ivConif,ivGrass   !For vegetated surfaces. Snow included although quite often LAI will be in its minimum when snow on ground
     active_veg_fr(iv) = (sfr(iv+2)*(1-snowFrac(iv+2)))*(LAI_day(iv)-LAIMin(iv))/(LAIMax(iv)-LAIMin(iv))
  ENDDO
  !DO iv=ivConif,ivGrass   !For vegetated surfaces. Snow is not included?
  !   active_veg_fr(iv) = sfr(iv+2)*(LAI(id-1,iv)-LAIMin(iv))/(LAIMax(iv)-LAIMin(iv))
  !ENDDO


  IF(EmissionsMethod>=11 .AND. EmissionsMethod<=16) THEN   ! Rectangular hyperbola

     ! Calculate carbon uptake due to photosynthesis -------------
     Fc_photo = 0
     DO iv=ivConif,ivGrass
        Fc_photo_surf(iv) = -beta_bioCO2(iv)*alpha_bioCO2(iv)*PAR_umolm2s1/(alpha_bioCO2(iv)*PAR_umolm2s1 + beta_bioCO2(iv))
        ! For active vegetation fraction only
        Fc_photo = Fc_photo + Fc_photo_surf(iv)*active_veg_fr(iv)  !umol m-2 s-1
     ENDDO

  ELSEIF(EmissionsMethod>=21 .AND. EmissionsMethod<=26) THEN  !Local model, Bellucco et al. (2017)

     ! Calculate carbon uptake due to photosynthesis -------------
     Bellucco2017_Pho = 0
     DO iv=ivConif,ivGrass
        Bellucco2017_Pho_surf(iv) = -(1/(2*theta_bioCO2(iv))*(alpha_bioCO2(iv)*PAR_umolm2s1+beta_bioCO2(iv)- &
             SQRT((alpha_bioCO2(iv)*PAR_umolm2s1+beta_bioCO2(iv))**2-4* &
             alpha_bioCO2(iv)*beta_bioCO2(iv)*theta_bioCO2(iv)*PAR_umolm2s1)))
        ! For active vegetation fraction only
        Bellucco2017_Pho = Bellucco2017_Pho + Bellucco2017_Pho_surf(iv)*active_veg_fr(iv)
     ENDDO

     !Fc_photo = Bellucco2017_Pho*active_veg_fr(ConifSurf-2)+ &
     !           Bellucco2017_Pho*active_veg_fr(DecidSurf-2)+ &
     !           Bellucco2017_Pho*active_veg_fr(GrassSurf-2)

     Fc_photo = Bellucco2017_Pho

  ELSEIF(EmissionsMethod>=31 .AND. EmissionsMethod<=36) THEN  !General model, Bellucco et al. (2017)
     ! Not currently recommended as includes also some anthropogenic impacts. Should maybe be separate from other biogen choices?
     ! Alpha and beta calculated as a function of vegetation cover fraction

     !Different alpha, beta and theta vegetation cover values should be same in BiogenCO2Method = 3.
     IF(alpha_bioCO2(ivConif) == alpha_bioCO2(ivDecid) .AND. alpha_bioCO2(ivConif) == alpha_bioCO2(ivGrass) .AND. &
          beta_bioCO2(ivConif) == beta_bioCO2(ivDecid) .AND. beta_bioCO2(ivConif) == beta_bioCO2(ivGrass) .AND. &
          theta_bioCO2(ivConif) == theta_bioCO2(ivDecid) .AND. theta_bioCO2(ivConif) == theta_bioCO2(ivGrass)) THEN

        !Because different alpha, beta and theta values are same - only one vegetation type is needed.
        alpha_bioCO2_v2(ivConif) = alpha_bioCO2(ivConif) + alpha_enh_bioCO2(ivConif)* &
             (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))      !umol CO2 umol photons-1
        beta_bioCO2_v2(ivConif) = -beta_bioCO2(ivConif) + beta_enh_bioCO2(ivConif)* &
             (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))     !umol m^-2 s^-1

        !Photosynthesis
        !Bellucco2017_Pho = -(1/(2*theta_bioCO2(ivConif))*(alpha_bioCO2(ivConif)*PAR_umolm2s1+beta_bioCO2(ivConif)- &
        !                   sqrt((alpha_bioCO2(ivConif)*PAR_umolm2s1+beta_bioCO2(ivConif))**2-4* &
        !                   alpha_bioCO2(ivConif)*beta_bioCO2(ivConif)*theta_bioCO2(ivConif)*PAR_umolm2s1)))

        Bellucco2017_Pho = -(1/(2*theta_bioCO2(ivConif))*(alpha_bioCO2_v2(ivConif)*PAR_umolm2s1+beta_bioCO2_v2(ivConif)- &
             SQRT((alpha_bioCO2_v2(ivConif)*PAR_umolm2s1+beta_bioCO2_v2(ivConif))**2-4* &
             alpha_bioCO2_v2(ivConif)*beta_bioCO2_v2(ivConif)*theta_bioCO2(ivConif)*PAR_umolm2s1)))

     ELSE !If values are not same, then weighted average is calculated.
        CALL ErrorHint(74,'Check values in SUEWS_BiogenCO2.txt: ',notUsed,notUsed,notUsedI)

        ! Weighted averages
        alpha_bioCO2_v2(ivConif) = (alpha_bioCO2(ivConif)*sfr(ConifSurf)/VegFracSum + &
             alpha_bioCO2(ivDecid)*sfr(DecidSurf)/VegFracSum&
             + alpha_bioCO2(ivGrass)*sfr(GrassSurf)/VegFracSum) &
             / (alpha_bioCO2(ivConif) + alpha_bioCO2(ivDecid) + alpha_bioCO2(ivGrass))
        beta_bioCO2_v2(ivConif)  = (beta_bioCO2(ivConif)*sfr(ConifSurf)/VegFracSum + &
             beta_bioCO2(ivDecid)*sfr(DecidSurf)/VegFracSum&
             + beta_bioCO2(ivGrass)*sfr(GrassSurf)/VegFracSum )/ (beta_bioCO2(ivConif) + &
             beta_bioCO2(ivDecid) + beta_bioCO2(ivGrass))
        theta_bioCO2_v2(ivConif) = (theta_bioCO2(ivConif)*sfr(ConifSurf)/VegFracSum  + &
             theta_bioCO2(ivDecid)*sfr(DecidSurf)/VegFracSum&
             + theta_bioCO2(ivGrass)*sfr(GrassSurf)/VegFracSum )/( theta_bioCO2(ivConif) + &
             theta_bioCO2(ivDecid) + theta_bioCO2(ivGrass))

        ! Using weighted average values to calculate alpha and beta as a function of vegetation cover fraction
        alpha_bioCO2_v2(ivConif) = alpha_bioCO2_v2(ivConif) + alpha_enh_bioCO2(ivConif)* &
             (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))     !umol CO2 umol photons-1
        beta_bioCO2_v2(ivConif) = -beta_bioCO2_v2(ivConif) + beta_enh_bioCO2(ivConif)* &
             (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))     !umol m^-2 s^-1

        !Photosynthesis
        !Bellucco2017_Pho = -(1/(2*theta_bioCO2(ivConif))*(alpha_bioCO2(ivConif)*PAR_umolm2s1+beta_bioCO2(ivConif)- &
        !                   sqrt((alpha_bioCO2(ivConif)*PAR_umolm2s1+beta_bioCO2(ivConif))**2-4* &
        !                   alpha_bioCO2(ivConif)*beta_bioCO2(ivConif)*theta_bioCO2(ivConif)*PAR_umolm2s1)))
        Bellucco2017_Pho = -(1/(2*theta_bioCO2_v2(ivConif))*(alpha_bioCO2_v2(ivConif)*PAR_umolm2s1+beta_bioCO2_v2(ivConif)- &
             SQRT((alpha_bioCO2_v2(ivConif)*PAR_umolm2s1+beta_bioCO2_v2(ivConif))**2-4* &
             alpha_bioCO2_v2(ivConif)*beta_bioCO2_v2(ivConif)*theta_bioCO2_v2(ivConif)*PAR_umolm2s1)))

     ENDIF

     ! Calculate carbon uptake due to photosynthesis -------------
     Fc_photo = Bellucco2017_Pho*active_veg_fr(ConifSurf-2)+ &
          Bellucco2017_Pho*active_veg_fr(DecidSurf-2)+ &
          Bellucco2017_Pho*active_veg_fr(GrassSurf-2)

  ENDIF

  ! Calculate carbon uptake due to respiration -------------
  Bellucco2017_Res = 0.0
  Bellucco2017_Res_surf = 0.0
  IF (VegFracSum>0.01) THEN
     DO iv=ivConif,ivGrass
        IF (sfr(2+iv)>0.005) THEN
           Bellucco2017_Res_surf = MAX(min_res_bioCO2(iv), resp_a(iv)*EXP(resp_b(iv)*Temp_C))
           ! For active vegetation fraction only
           ! Bellucco2017_Res = Bellucco2017_Res + Bellucco2017_Res_surf(iv)*active_veg_fr(iv)
           Bellucco2017_Res = Bellucco2017_Res + Bellucco2017_Res_surf * sfr(2+iv)/VegFracSum
        ENDIF
     ENDDO
  ENDIF
  !Fc_respi = Bellucco2017_Res * (sfr(ConifSurf)*(1-snowFrac(ConifSurf))+sfr(DecidSurf)*(1-snowFrac(DecidSurf))+ &
  !                                  sfr(GrassSurf)*(1-snowFrac(GrassSurf))+sfr(BSoilSurf)*(1-snowFrac(BSoilSurf)))
  Fc_respi = Bellucco2017_Res * (sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(BSoilSurf))

  ! Combine to find biogenic CO2 flux
  Fc_biogen = Fc_photo + Fc_respi


  RETURN

ENDSUBROUTINE CO2_biogen
!========================================================================================
