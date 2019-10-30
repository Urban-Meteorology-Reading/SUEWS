module rsl_module
   USE AtmMoistStab_module, ONLY: cal_Stab, stab_psi_mom, stab_psi_heat, stab_phi_mom, stab_phi_heat
   USE meteo, ONLY: RH2qa, qa2RH
   USE allocateArray, ONLY: &
      nsurf, BldgSurf, ConifSurf, DecidSurf
   implicit none

contains

   SUBROUTINE RSLProfile( &
      Zh, z0m, zdm, &
      UStar, L_MOD, sfr, planF, StabilityMethod, &
      avcp, lv_J_kg, &
      Temp_C, avRH, Press_hPa, zMeas, qh, qe, &  ! input
      T2_C, q2_gkg, U10_ms, RH2, &!output
      z0, zd, psihatm_z, psihath_z, & !output
      dataoutLineRSL) ! output
      !-----------------------------------------------------
      ! calculates windprofiles using MOST with a RSL-correction
      ! based on Harman & Finnigan 2007
      !
      ! last modified by:
      ! NT 16 Mar 2019: initial version
      ! TS 16 Oct 2019: improved consistency in parameters/varaibles within SUEWS
      !
      !-----------------------------------------------------

      IMPLICIT NONE

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in) ::sfr! surface fractions
      REAL(KIND(1d0)), INTENT(in):: zMeas  ! Air temperature/ moisture forcing height [m]
      REAL(KIND(1d0)), INTENT(in):: Temp_C ! Air temperature at forcing height [C]
      REAL(KIND(1d0)), INTENT(in):: avRH   ! relative humidity at forcing height [-]
      REAL(KIND(1d0)), INTENT(in):: Press_hPa ! pressure at forcing height [hPa]
      REAL(KIND(1d0)), INTENT(in):: UStar  ! Friction velocity [m s-1]
      REAL(KIND(1d0)), INTENT(in):: L_MOD  ! Obukhov length [m]
      REAL(KIND(1d0)), INTENT(in):: avcp  ! specific heat capacity [J kg-1 K-1]
      REAL(KIND(1d0)), INTENT(in):: lv_J_kg  ! Latent heat of vaporization in [J kg-1]
      REAL(KIND(1d0)), INTENT(in):: qh  ! sensible heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(in):: qe     ! Latent heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(in):: Zh     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: z0m     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: zdm     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: planF  ! Frontal area index [-]

      INTEGER, INTENT(in)::StabilityMethod

      REAL(KIND(1d0)), INTENT(out):: T2_C ! Air temperature at 2 m [C]
      REAL(KIND(1d0)), INTENT(out):: q2_gkg ! Air specific humidity at 2 m [g kg-1]
      REAL(KIND(1d0)), INTENT(out):: U10_ms ! wind speed at 10 m [m s-1]
      REAL(KIND(1d0)), INTENT(out):: RH2 ! Air relative humidity [-]

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
      REAL(KIND(1d0)), DIMENSION(nz), INTENT(out):: psihatm_z
      REAL(KIND(1d0)), DIMENSION(nz), INTENT(out):: psihath_z
      REAL(KIND(1d0)), DIMENSION(nz):: dif
      ! REAL(KIND(1d0)), DIMENSION(nz):: psihatm_z, psihath_z
      REAL(KIND(1d0)), DIMENSION(nz):: zarray
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineURSL ! wind speed array [m s-1]
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineTRSL ! Temperature array [C]
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineqRSL ! Specific humidity array [g kg-1]

      REAL(KIND(1d0)), INTENT(out)::z0  ! roughness length from H&F
      REAL(KIND(1d0)), INTENT(out)::zd ! zero-plane displacement

      REAL(KIND(1d0))::Lc_build, Lc_tree, Lc ! canopy drag length scale
      REAL(KIND(1d0))::Scc ! Schmidt number for temperature and humidity
      REAL(KIND(1d0))::phim, psimz, psimZh, psimz0, phi_hatmZh, phi_hathZh, phimzp, phimz, phihzp, phihz, psihz, psihza  ! stability function for momentum
      REAL(KIND(1d0))::betaHF, betaNL, beta, betaN2  ! beta coefficient from Harman 2012
      REAL(KIND(1d0))::elm ! mixing length
      REAL(KIND(1d0))::xxm1, xxm1_2, xxh1, xxh1_2, err, z01, dphi, dphih ! dummy variables for stability functions
      REAL(KIND(1d0))::f, cm, c2, ch, c2h ! H&F'07 and H&F'08 'constants'
      REAL(KIND(1d0))::t_h, q_h ! H&F'08 canopy corrections
      REAL(KIND(1d0))::TStar ! temperature scale
      REAL(KIND(1d0))::sfr_zh ! land cover fraction of bluff bodies: buildings and trees
      REAL(KIND(1d0))::sfr_tr ! land cover fraction of trees
      REAL(KIND(1d0))::L_MOD_RSL ! Obukhov length used in RSL module with thresholds applied
      REAL(KIND(1d0))::zH_RSL ! mean canyon height used in RSL module with thresholds applied
      REAL(KIND(1d0))::psihatm_zp
      REAL(KIND(1d0))::psihath_zp
      REAL(KIND(1d0))::dz! initial height step
      REAL(KIND(1d0)), parameter::Zh_min = 0.15! limit for minimum canyon height used in RSL module
      REAL(KIND(1d0)), parameter::ratio_dz = 1.618! ratio between neighbouring height steps

      REAL(KIND(1d0))::qa_gkg, qStar ! specific humidity scale
      INTEGER :: I, z, it, idx_can, idx_za, idx_2m, idx_10m
      !
      ! Step 1: Calculate grid-cell dependent constants
      ! Step 2: Calculate Beta (crucial for H&F method)
      ! Step 3: calculate the stability dependent H&F constants
      ! Step 4: determine psihat at levels above the canopy
      ! Step 5: Calculate z0 iteratively
      ! Step 6: Calculate mean variables above canopy
      ! Step 7: Calculate mean variables in canopy
      !
      ! ! Step 1
      ! ! Start setting up the parameters
      ! ! calculate Lc for tree grid fraction using eq 1 H&F'07 and rest of grid using C&B'04
      ! !
      ! ! under stable conditions, set a threshold for L_MOD to avoid numerical issues. TS 28 Oct 2019
      ! L_MOD_RSL = merge(L_MOD, max(300., L_MOD), L_MOD < 0)

      ! ! zH_RSL
      ! zH_RSL = max(zh, Zh_min)

      ! ! land cover fraction of bluff bodies
      ! sfr_zh = sum(sfr([BldgSurf, ConifSurf, DecidSurf]))
      ! ! set a threshold of 0.99 for sfr_zh to avoid numerical difficulties
      ! sfr_zh = min(sfr_zh, 0.99)

      ! ! land cover fraction of trees
      ! sfr_tr = sum(sfr([ConifSurf, DecidSurf]))

      ! ! height scale for buildings !not used? why?
      ! Lc_build = (1.-sfr(BldgSurf))/planF*Zh_RSL  ! Coceal and Belcher 2004 assuming Cd = 2

      ! ! height scale for tress
      ! Lc_tree = 1./(cd_tree*a_tree) ! not used? why?

      ! ! height scale for bluff bodies
      ! Lc = (1.-sfr_zh)/planF*Zh_RSL
      ! Lc = cal_lc(zh_RSL, planF, sfr)

      ! ! phim at Lc
      ! phim = stab_phi_mom(StabilityMethod, Lc/L_MOD_RSL)

      ! Scc = 0.5 + 0.3*TANH(2.*Lc/L_MOD_RSL)  ! Schmidt number Harman and Finnigan 2008: assuming the same for heat and momemntum
      ! f = 0.5*((1.+4.*r*Scc)**0.5) - 0.5
      ! !

      ! !
      ! ! Step 2:
      ! ! Parameterise beta according to Harman 2012 with upper limit of 0.5
      ! ! betaN for trees found to be 0.3 and for urban 0.4 linearly interpolate between the two using surface fractions
      ! ! betaN2 = 0.30 + (1.-sfr(ConifSurf) - sfr(ConifSurf))*0.1
      ! if (sfr_zh > 0) then
      !    betaN2 = 0.30*sfr_tr/sfr_zh + (sfr_zh - sfr_tr)/sfr_zh*0.4
      ! else
      !    betaN2 = 0.35 ! what about if sfr_zh==0?
      ! end if

      ! betaHF = betaN2/phim
      ! betaNL = (kappa/2.)/phim

      ! IF (Lc/L_MOD_RSL > a2) THEN
      !    beta = betaHF
      ! ELSE
      !    beta = betaNL + ((betaHF - betaNL)/(1.+a1*abs(Lc/L_MOD_RSL - a2)**a3))
      ! ENDIF

      ! IF (beta > 0.5) THEN
      !    beta = 0.5
      ! ENDIF
      ! ! zd = merge(Zh_RSL - (beta**2.)*Lc, zdm, Zh_RSL>zh_min) ! if Zh>0, calculate z0 using RSL method; use SUEWS-based zdm otherwise
      ! zd = Zh_RSL - (beta**2.)*Lc
      ! elm = 2.*beta**3*Lc

      call RSL_cal_prms( &
         StabilityMethod, zh, L_MOD, sfr, planF, &!input
         L_MOD_RSL,zH_RSL,Lc, beta, zd, elm,Scc,f)

      ! Define the height array
      !
      ! IF ((3.*Zh_RSL) < 10.) THEN
      !    dz = min(1./3.,Zh_RSL/2.)      ! if canopy height is small use steps of 0.33333 m to get to 10 m
      !    zarray = (/(I, I=1, nz)/)*dz
      ! ELSE
      !    dz = Zh_RSL/10.
      !    zarray = (/(I, I=1, nz)/)*dz
      ! ENDIF
      IF (Zh_RSL <= 2) THEN
         ! add key heights for within canyon calculations
         zarray(1) = Zh_RSL/4.
         zarray(2) = Zh_RSL/2.
         zarray(3) = merge(0.99*Zh_RSL, Zh_RSL, Zh_RSL == 2.)
         ! define heights for above canyon calculations
         zarray(4) = 2.
      ELSE IF (Zh_RSL <= 10) THEN
         ! add key heights for within canyon calculations
         zarray(1) = 2.
         zarray(2) = (2 + Zh_RSL)/2.
         zarray(3) = merge(0.99*Zh_RSL, Zh_RSL, Zh_RSL == 10.)
         ! define heights for above canyon calculations
         zarray(4) = 10.
      else
         ! add key heights for within canyon calculations
         zarray(1) = 2.
         zarray(2) = 10.
         zarray(3) = (10 + Zh_RSL)/2.
         ! define heights for above canyon calculations
         zarray(4) = Zh_RSL
      ENDIF
      ! fill up other heights in zarray
      dz = (zMeas - zarray(4))/(nz - 4.)
      do i = 5, nz
         zarray(i) = zarray(4) + (i - 4)*dz
      end do

      ! add key heights (2m and 10m) to zarray
      ! 2m:
      DO z = 1, nz
         dif(z) = ABS(zarray(z) - 2)
      ENDDO
      idx_2m = MINLOC(dif, DIM=1)
      zarray(idx_2m) = 2
      ! 10m:
      DO z = 1, nz
         dif(z) = ABS(zarray(z) - 10)
      ENDDO
      idx_10m = MINLOC(dif, DIM=1)
      zarray(idx_10m) = 10

      ! determine index at the canyon top
      DO z = 1, nz
         dif(z) = ABS(zarray(z) - Zh_RSL)
      ENDDO
      idx_can = MINLOC(dif, DIM=1)
      zarray(idx_can) = Zh_RSL

      ! determine index at measurement height
      DO z = 1, nz
         ! dif2(z) = ABS(zarray(z) - (zMeas - zd))
         dif(z) = ABS(zarray(z) - zMeas)
      ENDDO
      idx_za = MINLOC(dif, DIM=1)
      ! zarray(idx_za) = zMeas - zd
      zarray(idx_za) = zMeas

      !
      ! Step 3:
      !
      psimZh = stab_psi_mom(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)

      ! calculate phihatM according to H&F '07 and H&F '08 for heat and humidity
      xxm1 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)
      xxm1_2 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd + 1.)/L_MOD_RSL)

      xxh1 = stab_phi_heat(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)
      xxh1_2 = stab_phi_heat(StabilityMethod, (Zh_RSL - zd + 1.)/L_MOD_RSL)

      phi_hatmZh = kappa/(2.*beta*xxm1)
      phi_hathZh = kappa*Scc/(2.*beta*xxh1)
      dphi = xxm1_2 - xxm1
      dphih = xxh1_2 - xxh1
      IF (phi_hatmZh > 1.) THEN
         c2 = 0.5 ! more stable, but less correct
         c2h = 0.5
      ELSE
         c2 = (kappa*(3.-(2.*beta**2.*Lc/xxm1*dphi)))/(2.*beta*xxm1 - kappa)  ! if very unstable this might cause some high values of psihat_z
         c2h = (kappa*Scc*(2.+f - (dphih*2.*beta**2.*Lc/xxh1)))/(2.*beta*xxh1 - kappa*Scc)
      ENDIF
      cm = (1.-phi_hatmZh)*EXP(c2/2.)
      ch = (1.-phi_hathzh)*EXP(c2h/2.)

      !
      ! Step 4:
      !
      psihatm_z = 0.*zarray
      psihath_z = 0.*zarray
      psihatm_zp = 0.
      psihath_zp = 0
      DO z = nz - 1, idx_can - 1, -1
         phimz = stab_phi_mom(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
         phimzp = stab_phi_mom(StabilityMethod, (zarray(z + 1) - zd)/L_MOD_RSL)
         phihz = stab_phi_heat(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
         phihzp = stab_phi_heat(StabilityMethod, (zarray(z + 1) - zd)/L_MOD_RSL)

         psihatm_z(z) = psihatm_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
                        /(zarray(z + 1) - zd)
         psihatm_z(z) = psihatm_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - zd)/elm)) &
                        /(zarray(z) - zd)
         ! psihatm_zp=psihatm_z(z)
         psihath_z(z) = psihath_z(z + 1) + dz/2.*phihzp*(ch*EXP(-1.*c2h*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
                        /(zarray(z + 1) - zd)
         psihath_z(z) = psihath_z(z) + dz/2.*phihz*(ch*EXP(-1.*c2h*beta*(zarray(z) - zd)/elm)) &
                        /(zarray(z) - zd)
         ! psihath_zp=psihath_z(z)
      ENDDO
      !
      ! Step 5
      ! calculate z0 iteratively
      !
      ! z0 = 0.5  !first guess
      z0 = z0m  !first guess
      ! if (Zh_RSL>Zh_min) then
      ! if Zh>Zh_min, calculate z0 using RSL method
      err = 10.
      psimz0 = 0.5
      it = 1
      DO WHILE ((err > 0.001) .AND. (it < 10))
         psimz0 = stab_psi_mom(StabilityMethod, z0/L_MOD_RSL)
         z01 = z0
         z0 = (Zh_RSL - zd)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0)*EXP(psihatm_z(idx_can))
         err = ABS(z01 - z0)
         it = it + 1
      ENDDO
      ! endif

      psimz0 = stab_psi_mom(StabilityMethod, z0/L_MOD_RSL)
      psihza = stab_psi_heat(StabilityMethod, (zMeas - zd)/L_MOD_RSL)
      TStar = -1.*(qh/(avcp))/UStar
      qStar = -1.*(qe/lv_J_kg)/UStar
      qa_gkg = RH2qa(avRH/100, Press_hPa, Temp_c)
      !
      ! Step 6
      ! calculate above canopy variables
      !
      DO z = idx_can, nz
         psimz = stab_psi_mom(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
         psihz = stab_psi_heat(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
         ! dataoutLineURSL(z) = (LOG((zarray(z) - zd)/z0) - psimz + psimz0 - psihatm_z(z) + psihatm_z(idx_can))/kappa ! this is different from Theeuwes et al. (2019 BLM)
         dataoutLineURSL(z) = (LOG((zarray(z) - zd)/z0) - psimz + psimz0 - psihatm_z(z))/kappa ! eqn. 3 in Theeuwes et al. (2019 BLM)
         ! dataoutLineTRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za - 1))/kappa
         ! dataoutLineqRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za - 1))/kappa
         dataoutLineTRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za))/kappa ! eqn. 4 in Theeuwes et al. (2019 BLM)
         dataoutLineqRSL(z) = (LOG((zarray(z) - zd)/(zMeas - zd)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za))/kappa
      ENDDO
      !
      ! Step 7
      ! calculate in canopy variables
      !
      t_h = Scc*TStar/(beta*f)
      q_h = Scc*qStar/(beta*f)
      DO z = 1, idx_can
         dataoutLineURSL(z) = dataoutLineURSL(idx_can)*EXP(beta*(zarray(z) - Zh_RSL)/elm)
         dataoutLineTRSL(z) = ((dataoutLineTRSL(idx_can)*TStar) + t_h*EXP(beta*f*(zarray(z) - Zh_RSL)/elm) - t_h)/TStar
         dataoutLineqRSL(z) = ((dataoutLineqRSL(idx_can)*qStar) + q_h*EXP(beta*f*(zarray(z) - Zh_RSL)/elm) - q_h)/qStar
      ENDDO

      dataoutLineURSL = dataoutLineURSL*UStar
      dataoutLineTRSL = dataoutLineTRSL*TStar + Temp_C
      dataoutLineqRSL = (dataoutLineqRSL*qStar + qa_gkg/1000.)*1000.

      dataoutLineRSL = (/zarray, dataoutLineURSL, dataoutLineTRSL, dataoutLineqRSL/)
      ! zarrays = (/zarray, zarray, zarray/)

      !
      ! Step 8
      ! retrieve the diagnostics at key heights
      !
      T2_C = dataoutLineTRSL(idx_2m)
      q2_gkg = dataoutLineqRSL(idx_2m)
      U10_ms = dataoutLineURSL(idx_10m)
      ! get relative humidity:
      RH2 = qa2RH(q2_gkg, press_hPa, T2_C)

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

   recursive function cal_psim_hat(StabilityMethod, z, zh_RSL, L_MOD_RSL, beta, Lc) result(psim_hat_z)
      ! calculate surface/skin temperature
      ! TS, 23 Oct 2019
      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) :: z ! height of interest [m]
      real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
      real(KIND(1D0)), intent(in) ::  L_MOD_RSL ! Obukhov length [m]

      ! output
      real(KIND(1D0)) ::psim_hat_z ! psim_hat at height of interest

      ! internal variables
      real(KIND(1D0)) ::zp ! a height above z used for iterative calculations
      real(KIND(1D0)) ::zd ! displacement height used in RSL
      real(KIND(1D0)) ::phim_lc ! displacement height used in RSL
      real(KIND(1D0)) ::phim_z ! displacement height used in RSL
      real(KIND(1D0)) ::phim_zp ! displacement height used in RSL
      real(KIND(1D0)) ::psim_hat_zp ! displacement height used in RSL
      real(KIND(1D0)) ::elm ! displacement height used in RSL
      real(KIND(1D0)) ::xxm1 ! displacement height used in RSL
      real(KIND(1D0)) ::xxm1_2 ! displacement height used in RSL
      real(KIND(1D0)) ::dphi ! displacement height used in RSL
      real(KIND(1D0)) ::phi_hatmZh ! displacement height used in RSL
      real(KIND(1D0)) ::cm ! displacement height used in RSL
      real(KIND(1D0)) ::c2 ! displacement height used in RSL

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

      zd = Zh_RSL - (beta**2.)*Lc
      elm = 2.*beta**3*Lc

      ! phim at Lc
      phim_lc = stab_phi_mom(StabilityMethod, Lc/L_MOD_RSL)

      phim_z = stab_phi_mom(StabilityMethod, (z - zd)/L_MOD_RSL)
      phim_zp = stab_phi_mom(StabilityMethod, (zp - zd)/L_MOD_RSL)

      psim_hat_zp = cal_psim_hat(StabilityMethod, zp, zh_RSL, L_MOD_RSL, beta, Lc)

      xxm1 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)
      xxm1_2 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd + 0.01)/L_MOD_RSL)

      dphi = (xxm1_2 - xxm1)/0.01

      phi_hatmZh = kappa/(2.*beta*xxm1)


      IF (phi_hatmZh > 1.) THEN
         ! more stable, but less correct
         c2 = 0.5
      ELSE
         ! if very unstable this might cause some high values of psihat_z
         c2 = (kappa*(3.-(2.*beta**2.*Lc/xxm1*dphi)))/(2.*beta*xxm1 - kappa)
      ENDIF

      cm = (1.-phi_hatmZh)*EXP(c2/2.)

      !Taylor's approximation for integral
      psim_hat_z = psim_hat_zp + dz/2.*phim_zp*(cm*EXP(-1.*c2*beta*(zp - zd)/elm))/(zp - zd)
      psim_hat_z = psim_hat_z + dz/2.*phim_z*(cm*EXP(-1.*c2*beta*(z - zd)/elm))/(z - zd)

   end function cal_psim_hat

   function cal_Lc(zH_RSL, planF, sfr) result(Lc)
      ! calculate surface/skin temperature
      ! TS, 23 Oct 2019
      implicit none
      real(KIND(1D0)), intent(in) ::  zH_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  planF ! height scale for bluff bodies [m]
      real(KIND(1D0)), DIMENSION(nsurf), intent(in) ::  sfr ! land cover fractions

      ! output
      real(KIND(1D0)) ::Lc

      ! internal variables
      real(KIND(1D0)) ::sfr_zh
      real(KIND(1D0)) ::sfr_tr

      ! REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      ! REAL(KIND(1d0)), PARAMETER::r = 0.1
      ! REAL(KIND(1d0)), PARAMETER::a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1.

      ! land cover fraction of bluff bodies
      sfr_zh = sum(sfr([BldgSurf, ConifSurf, DecidSurf]))
      ! set a threshold of 0.99 for sfr_zh to avoid numerical difficulties
      sfr_zh = min(sfr_zh, 0.99)

      ! land cover fraction of trees
      sfr_tr = sum(sfr([ConifSurf, DecidSurf]))

      ! height scale for buildings !not used? why?
      ! Lc_build = (1.-sfr(BldgSurf))/planF*Zh_RSL  ! Coceal and Belcher 2004 assuming Cd = 2

      ! height scale for tress
      ! Lc_tree = 1./(cd_tree*a_tree) ! not used? why?

      ! height scale for bluff bodies
      Lc = (1.-sfr_zh)/planF*Zh_RSL

   end function cal_Lc

   subroutine RSL_cal_prms( &
      StabilityMethod, zh, L_MOD, sfr, planF, &!input
      L_MOD_RSL,zH_RSL,Lc, beta, zd, elm,Scc,f)!output
      ! calculate surface/skin temperature
      ! TS, 23 Oct 2019
      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) ::  zh ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  planF ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  L_MOD ! Obukhov length [m]
      real(KIND(1D0)), DIMENSION(nsurf), intent(in) ::  sfr ! land cover fractions

      ! output
      real(KIND(1D0)), intent(out) ::L_MOD_RSL
      real(KIND(1D0)), intent(out) ::zH_RSL
      real(KIND(1D0)), intent(out) ::Lc
      real(KIND(1D0)), intent(out) ::beta ! psim_hat at height of interest
      real(KIND(1D0)), intent(out) ::zd
      real(KIND(1D0)), intent(out) ::elm

      ! internal variables
      real(KIND(1D0)) ::sfr_zh
      real(KIND(1D0)) ::sfr_tr
      real(KIND(1D0)) ::phim
      real(KIND(1D0)) ::Scc
      real(KIND(1D0)) ::f
      real(KIND(1D0)) ::betaN2
      real(KIND(1D0)) ::betaHF
      real(KIND(1D0)) ::betaNL

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      REAL(KIND(1d0)), PARAMETER::r = 0.1
      REAL(KIND(1d0)), PARAMETER::a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1.
      REAL(KIND(1d0)), parameter::Zh_min = 0.15! limit for minimum canyon height used in RSL module

      ! under stable conditions, set a threshold for L_MOD to avoid numerical issues. TS 28 Oct 2019
      L_MOD_RSL = merge(L_MOD, max(300., L_MOD), L_MOD < 0)

      ! zH_RSL
      zH_RSL = max(zh, Zh_min)

      ! land cover fraction of bluff bodies
      sfr_zh = sum(sfr([BldgSurf, ConifSurf, DecidSurf]))
      ! set a threshold of 0.99 for sfr_zh to avoid numerical difficulties
      sfr_zh = min(sfr_zh, 0.99)

      ! land cover fraction of trees
      sfr_tr = sum(sfr([ConifSurf, DecidSurf]))

      ! height scale for buildings !not used? why?
      ! Lc_build = (1.-sfr(BldgSurf))/planF*Zh_RSL  ! Coceal and Belcher 2004 assuming Cd = 2

      ! height scale for tress
      ! Lc_tree = 1./(cd_tree*a_tree) ! not used? why?

      ! height scale for bluff bodies
      Lc = (1.-sfr_zh)/planF*Zh_RSL

      ! phim at Lc
      phim = stab_phi_mom(StabilityMethod, Lc/L_MOD_RSL)

      Scc = 0.5 + 0.3*TANH(2.*Lc/L_MOD_RSL)  ! Schmidt number Harman and Finnigan 2008: assuming the same for heat and momemntum
      f = 0.5*((1.+4.*r*Scc)**0.5) - 0.5

      !
      ! Step 2:
      ! Parameterise beta according to Harman 2012 with upper limit of 0.5
      ! betaN for trees found to be 0.3 and for urban 0.4 linearly interpolate between the two using surface fractions
      ! betaN2 = 0.30 + (1.-sfr(ConifSurf) - sfr(ConifSurf))*0.1
      if (sfr_zh > 0) then
         betaN2 = 0.30*sfr_tr/sfr_zh + (sfr_zh - sfr_tr)/sfr_zh*0.4
      else
         betaN2 = 0.35 ! what about if sfr_zh==0?
      end if

      betaHF = betaN2/phim
      betaNL = (kappa/2.)/phim

      IF (Lc/L_MOD_RSL > a2) THEN
         beta = betaHF
      ELSE
         beta = betaNL + ((betaHF - betaNL)/(1.+a1*abs(Lc/L_MOD_RSL - a2)**a3))
      ENDIF

      IF (beta > 0.5) THEN
         beta = 0.5
      ENDIF

      zd = Zh_RSL - (beta**2.)*Lc
      elm = 2.*beta**3*Lc

   end subroutine RSL_cal_prms

end module rsl_module
