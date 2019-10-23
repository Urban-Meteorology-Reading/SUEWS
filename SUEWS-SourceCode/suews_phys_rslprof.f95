module rsl_module
   implicit none

contains

   SUBROUTINE RSLProfile( &
      UStar, L_MOD, sfr, Zh, planF, StabilityMethod, &
      avcp,lv_J_kg, &
      Temp_C, avRH, Press_hPa, zMeas, qh, qe, &  ! input
      dataoutLineRSL) ! output
      !-----------------------------------------------------
      ! calculates windprofiles using MOST with a RSL-correction
      ! based on Harman & Finnigan 2007
      !
      ! last modified by:
      ! NT 16 Mar 2019
      ! TS 16 Oct 2019: improved consistency in parameters/varaibles within SUEWS
      !
      !-----------------------------------------------------
      USE AtmMoistStab_module, ONLY: cal_Stab, stab_psi_mom, stab_psi_heat, stab_phi_mom, stab_phi_heat
      USE meteo, ONLY: RH2qa

      IMPLICIT NONE
      INTEGER, PARAMETER:: nsurf = 7 ! number of surface types
      INTEGER, PARAMETER:: BldgSurf = 2
      INTEGER, PARAMETER:: ConifSurf = 3
      INTEGER, PARAMETER:: DecidSurf = 4

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in) ::sfr! surface fractions
      REAL(KIND(1d0)), INTENT(in):: zMeas  ! Air temperature/ moisture forcing height [m]
      REAL(KIND(1d0)), INTENT(in):: Temp_C ! Air temperature at forcing height [C]
      REAL(KIND(1d0)), INTENT(in):: avRH   ! relative humidity at forcing height [-]
      REAL(KIND(1d0)), INTENT(in):: Press_hPa ! pressure at forcing height [hPa]
      REAL(KIND(1d0)), INTENT(in):: UStar  ! Friction velocity [m s-1]
      REAL(KIND(1d0)), INTENT(in):: L_MOD  ! Obukhov length [m]
      REAL(KIND(1d0)), INTENT(in):: avcp  ! specific heat capacity [J kg-1 K-1]
      REAL(KIND(1d0)), INTENT(in):: qh  ! sensible heat flux [W m-2]
      ! REAL(KIND(1d0)), INTENT(in):: TStar  ! Temperature scale [K]
      REAL(KIND(1d0)), INTENT(in):: lv_J_kg  ! Latent heat of vaporization in [J kg-1]
      REAL(KIND(1d0)), INTENT(in):: qe     ! Latent heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(in):: Zh     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: planF  ! Frontal area index [-]
      INTEGER, INTENT(in)::StabilityMethod

      REAL(KIND(1d0)), PARAMETER:: cd_tree = 1.2, & ! drag coefficient tree canopy !!!!needs adjusting!!!
                                   a_tree = 0.05, & ! the foliage area per unit volume !!!!needs adjusting!!!
                                   kappa = 0.40, &! von karman constant
                                 !   lv_J_kg = 2.5E6, &! latent heat for water vapor!!! make consistant with rest of code
                                   beta_N = 0.40, &  ! H&F beta coefficient in neutral conditions from Theeuwes et al., 2019 BLM
                                   pi = 4.*ATAN(1.0), r = 0.1, &
                                   a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1. ! constraints to determine beta

      INTEGER, PARAMETER :: nz = 30   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy

      ! REAL(KIND(1d0)), INTENT(out), DIMENSION(nz*3):: zarrays ! Height array
      REAL(KIND(1d0)), INTENT(out), DIMENSION(nz*4):: dataoutLineRSL  ! Variables array (/U,T,q/)
      REAL(KIND(1d0)), DIMENSION(nz):: dif, dif2, psihat_z, psihath_z, zarray, &
                                       dataoutLineURSL, & ! wind speed array [m s-1]
                                       dataoutLineTRSL, & ! Temperature array [C]
                                       dataoutLineqRSL    ! Specific humidity array [g kg-1]

      REAL(KIND(1d0)):: zd, & ! displacement height
                        Lc_build, Lc_tree, Lc, & ! canopy drag length scale
                        Scc, & ! Schmidt number for temperature and humidity
                        dz, & ! height steps
                        phim, psimz, psimZh, psimz0, phi_hatmZh, phi_hathZh, phimzp, phimz, phihzp, phihz, psihz, psihza, &  ! stability function for momentum
                        betaHF, betaNL, beta, betaN2, &  ! beta coefficient from Harman 2012
                        elm, & ! mixing length
                        xx1, xx1_2, xxh1, xxh1_2, err, z01, dphi, dphih, &  ! dummy variables for stability functions
                        z0, &  ! roughness length from H&F
                        f, cm, c2, ch, c2h, & ! H&F'07 and H&F'08 'constants'
                        t_h, q_h, & ! H&F'08 canopy corrections
                        TStar,&
                        qa_gkg, qStar ! specific humidity scale
      INTEGER :: I, z, it, idx_can, idx_za
      !
      ! Step 1: Calculate grid-cel dependent constants
      ! Step 2: Calculate Beta (crucial for H&F method)
      ! Step 3: calculate the stability dependent H&F constants
      ! Step 4: determine psihat at levels above the canopy
      ! Step 5: Calculate z0 iteratively
      ! Step 6: Calculate mean variables above canopy
      ! Step 7: Calculate mean variables in canopy
      !
      ! Step 1
      ! Start setting up the parameters
      ! calculate Lc for tree grid fraction using eq 1 H&F'07 and rest of grid using C&B'04
      !
      Lc_build = (1.-sfr(BldgSurf))/planF*Zh  ! Coceal and Belcher 2004 assuming Cd = 2
      Lc_tree = 1./(cd_tree*a_tree)
      Lc = (1.-(sfr(BldgSurf) + sfr(ConifSurf) + sfr(ConifSurf)))/planF*Zh
      Scc = 0.5 + 0.3*TANH(2.*Lc/L_MOD)  ! Schmidt number Harman and Finnigan 2008: assuming the same for heat and momemntum
      f = 0.5*((1.+4.*r*Scc)**0.5) - 0.5
      !
      ! Define the height array
      !
      IF ((3.*Zh) < 10.) THEN
         dz = 1./3.      ! if canopy height is small use steps of 0.33333 m to get to 10 m
         zarray = (/(I, I=1, nz)/)*dz
      ELSE
         dz = Zh/10.
         zarray = (/(I, I=1, nz)/)*dz
      ENDIF

      DO z = 1, nz
         dif(z) = ABS(zarray(z) - Zh)
      ENDDO
      idx_can = MINLOC(dif, DIM=1)
      phim = stab_phi_mom(StabilityMethod, Lc/L_MOD)
      !
      ! Step 2:
      ! Parameterise beta according to Harman 2012 with upper limit of 0.5
      ! betaN for trees found to be 0.3 and for urban 0.4 linearly interpolate between the two using surface fractions
      betaN2 = 0.30 + (1.-sfr(ConifSurf) - sfr(ConifSurf))*0.1

      betaHF = betaN2/phim
      betaNL = (kappa/2.)/phim

      IF (Lc/L_MOD > a2) THEN
         beta = betaHF
      ELSE
         beta = betaNL + ((betaHF - betaNL)/(1.+a1*abs(Lc/L_MOD - a2)**a3))
      ENDIF

      IF (beta > 0.5) THEN
         beta = 0.5
      ENDIF
      zd = Zh - (beta**2.)*Lc
      elm = 2.*beta**3*Lc

      DO z = 1, nz
         dif2(z) = ABS(zarray(z) - (zMeas - zd))
      ENDDO
      idx_za = MINLOC(dif2, DIM=1)

      !
      ! Step 3:
      !
      psimZh = stab_psi_mom(StabilityMethod, (Zh - zd)/L_MOD)

      ! calculate phihatM according to H&F '07 and H&F '08 for heat and humidity
      xx1 = stab_phi_mom(StabilityMethod, (Zh - zd)/L_MOD)
      xx1_2 = stab_phi_mom(StabilityMethod, (Zh - zd + 1.)/L_MOD)

      xxh1 = stab_phi_heat(StabilityMethod, (Zh - zd)/L_MOD)
      xxh1_2 = stab_phi_heat(StabilityMethod, (Zh - zd + 1.)/L_MOD)

      phi_hatmZh = kappa/(2.*beta*xx1)
      phi_hathZh = kappa*Scc/(2.*beta*xxh1)
      dphi = xx1_2 - xx1
      dphih = xxh1_2 - xxh1
      IF (phi_hatmZh > 1.) THEN
         c2 = 0.5 ! more stable, but less correct
         c2h = 0.5
      ELSE
         c2 = (kappa*(3.-(2.*beta**2.*Lc/xx1*dphi)))/(2.*beta*xx1 - kappa)  ! if very unstable this might cause some high values of psihat_z
         c2h = (kappa*Scc*(2.+f - (dphih*2.*beta**2.*Lc/xxh1)))/(2.*beta*xxh1 - kappa*Scc)
      ENDIF
      cm = (1.-phi_hatmZh)*EXP(c2/2.)
      ch = (1.-phi_hathzh)*EXP(c2h/2.)

      !
      ! Step 4:
      !
      psihat_z = 0.*zarray
      DO z = nz - 1, idx_can - 1, -1
         phimz = stab_phi_mom(StabilityMethod, (zarray(z) - zd)/L_MOD)
         phimzp = stab_phi_mom(StabilityMethod, (zarray(z + 1) - zd)/L_MOD)
         phihz = stab_phi_heat(StabilityMethod, (zarray(z) - zd)/L_MOD)
         phihzp = stab_phi_heat(StabilityMethod, (zarray(z + 1) - zd)/L_MOD)

         psihat_z(z) = psihat_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
                       /(zarray(z + 1) - zd)
         psihat_z(z) = psihat_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - zd)/elm)) &
                       /(zarray(z) - zd)
         psihath_z(z) = psihath_z(z + 1) + dz/2.*phihzp*(ch*EXP(-1.*c2h*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
                        /(zarray(z + 1) - zd)
         psihath_z(z) = psihath_z(z) + dz/2.*phihz*(ch*EXP(-1.*c2h*beta*(zarray(z) - zd)/elm)) &
                        /(zarray(z) - zd)
      ENDDO
      !
      ! Step 5
      ! calculate z0 iteratively
      !
      z0 = 0.5  !first guess
      err = 10.
      psimz0 = 0.5
      it = 1
      DO WHILE ((err > 0.001) .AND. (it < 10))
         psimz0 = stab_psi_mom(StabilityMethod, z0/L_MOD)
         z01 = z0
         z0 = (Zh - zd)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0)*EXP(psihat_z(idx_can))
         err = ABS(z01 - z0)
         it = it + 1
      ENDDO

      psimz0 = stab_psi_mom(StabilityMethod, z0/L_MOD)
      psihza = stab_psi_heat(StabilityMethod, (zMeas - zd)/L_MOD)
      TStar = -1.*(qh/(avcp))/UStar
      qStar = -1.*(qe/lv_J_kg)/UStar
      qa_gkg = RH2qa(avRH/100, Press_hPa, Temp_c)
      !
      ! Step 6
      ! calculate above canopy wind speed
      !
      DO z = idx_can, nz
         psimz = stab_psi_mom(StabilityMethod, (zarray(z) - zd)/L_MOD)
         psihz = stab_psi_heat(StabilityMethod, (zarray(z) - zd)/L_MOD)
         dataoutLineURSL(z) = (LOG((zarray(z) - zd)/z0) - psimz + psimz0 - psihat_z(z) + psihat_z(idx_can))/kappa
         dataoutLineTRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za - 1))/kappa
         dataoutLineqRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za - 1))/kappa
      ENDDO
      !
      ! Step 7
      ! calculate in canopy wind speed
      !
      t_h = Scc*TStar/(beta*f)
      q_h = Scc*qStar/(beta*f)
      DO z = 1, idx_can
         dataoutLineURSL(z) = dataoutLineURSL(idx_can)*EXP(beta*(zarray(z) - Zh)/elm)
         dataoutLineTRSL(z) = ((dataoutLineTRSL(idx_can)*TStar) + t_h*EXP(beta*f*(zarray(z) - Zh)/elm) - t_h)/TStar
         dataoutLineqRSL(z) = ((dataoutLineqRSL(idx_can)*qStar) + q_h*EXP(beta*f*(zarray(z) - Zh)/elm) - q_h)/qStar
      ENDDO

      dataoutLineURSL = dataoutLineURSL*UStar
      dataoutLineTRSL = dataoutLineTRSL*TStar + Temp_C
      dataoutLineqRSL = (dataoutLineqRSL*qStar + qa_gkg/1000.)*1000.

      dataoutLineRSL = (/zarray, dataoutLineURSL, dataoutLineTRSL, dataoutLineqRSL/)
      ! zarrays = (/zarray, zarray, zarray/)

! print *, 'Wind speed', dataoutLineURSL
      ! DO z = 1, nz
      !    print *, dataoutLineTRSL(z)
      ! ENDDO
      ! DO z = 1, nz
      !    print *, zarray(z)
      ! ENDDO
! print *, 'Temperature', Temp_C, dataoutLineTRSL
! print *, 'qStar', qStar, qe
! print *, 'humidity' , qa_gkg, dataoutLineqRSL*1000.
   END SUBROUTINE RSLProfile

end module rsl_module
