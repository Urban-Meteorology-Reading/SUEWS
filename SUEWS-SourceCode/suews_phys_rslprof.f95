module rsl_module
   USE AtmMoistStab_module, ONLY: cal_Stab, stab_psi_mom, stab_psi_heat, stab_phi_mom, stab_phi_heat
   USE meteo, ONLY: RH2qa, qa2RH
   use resist_module, only: SUEWS_cal_RoughnessParameters
   USE allocateArray, ONLY: &
      nsurf, BldgSurf, ConifSurf, DecidSurf, ncolumnsDataOutRSL
   implicit none

   INTEGER, PARAMETER :: nz = 30   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy

contains

   SUBROUTINE RSLProfile( &
      Zh, z0m, zdm, &
      L_MOD, sfr, FAI, StabilityMethod, &
      avcp, lv_J_kg, avdens, &
      avU1, Temp_C, avRH, Press_hPa, zMeas, qh, qe, &  ! input
      T2_C, q2_gkg, U10_ms, RH2, &!output
      dataoutLineRSL) ! output
      !-----------------------------------------------------
      ! calculates windprofiles using MOST with a RSL-correction
      ! based on Harman & Finnigan 2007
      !
      ! last modified by:
      ! NT 16 Mar 2019: initial version
      ! TS 16 Oct 2019: improved consistency in parameters/varaibles within SUEWS
      ! TODO how to improve the speed of this code
      !
      !-----------------------------------------------------

      IMPLICIT NONE

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in) ::sfr! surface fractions [-]
      REAL(KIND(1d0)), INTENT(in):: zMeas  ! height of atmospheric forcing [m]
      REAL(KIND(1d0)), INTENT(in):: avU1   ! Wind speed at forcing height [m s-1]
      REAL(KIND(1d0)), INTENT(in):: Temp_C ! Air temperature at forcing height [C]
      REAL(KIND(1d0)), INTENT(in):: avRH   ! relative humidity at forcing height [-]
      REAL(KIND(1d0)), INTENT(in):: Press_hPa ! pressure at forcing height [hPa]
      REAL(KIND(1d0)), INTENT(in):: L_MOD  ! Obukhov length [m]
      REAL(KIND(1d0)), INTENT(in):: avcp  ! specific heat capacity [J kg-1 K-1]
      REAL(KIND(1d0)), INTENT(in):: lv_J_kg  ! Latent heat of vaporization in [J kg-1]
      REAL(KIND(1d0)), INTENT(in):: avdens  ! air density [kg m-3]
      REAL(KIND(1d0)), INTENT(in):: qh  ! sensible heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(in):: qe     ! Latent heat flux [W m-2]
      REAL(KIND(1d0)), INTENT(in):: Zh     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: z0m     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: zdm     ! Mean building height [m]
      REAL(KIND(1d0)), INTENT(in):: FAI  ! Frontal area index [-]

      INTEGER, INTENT(in)::StabilityMethod

      REAL(KIND(1d0)), INTENT(out):: T2_C ! Air temperature at 2 m [C]
      REAL(KIND(1d0)), INTENT(out):: q2_gkg ! Air specific humidity at 2 m [g kg-1]
      REAL(KIND(1d0)), INTENT(out):: U10_ms ! wind speed at 10 m [m s-1]
      REAL(KIND(1d0)), INTENT(out):: RH2 ! Air relative humidity [-]

      INTEGER, PARAMETER :: nz = 30   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy

      REAL(KIND(1d0)), PARAMETER:: cd_tree = 1.2, & ! drag coefficient tree canopy !!!!needs adjusting!!!
                                   a_tree = 0.05, & ! the foliage area per unit volume !!!!needs adjusting!!!
                                   kappa = 0.40, &! von karman constant
                                   !   lv_J_kg = 2.5E6, &! latent heat for water vapor!!! make consistant with rest of code
                                   beta_N = 0.40, &  ! H&F beta coefficient in neutral conditions from Theeuwes et al., 2019 BLM
                                   pi = 4.*ATAN(1.0), r = 0.1, &
                                   a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1. ! constraints to determine beta

      ! Variables array [z,U,T,q, 12 debug vars]
      ! z: height array
      ! U,T,q: wind speed, air temp, specific humidity at z;
      ! debug vars: see dataoutLineRSL
      REAL(KIND(1d0)), INTENT(out), DIMENSION(ncolumnsDataOutRSL - 5):: dataoutLineRSL
      REAL(KIND(1d0)), DIMENSION(nz):: psihatm_z
      REAL(KIND(1d0)), DIMENSION(nz):: psihath_z
      REAL(KIND(1d0)), DIMENSION(nz):: dif
      ! REAL(KIND(1d0)), DIMENSION(nz):: psihatm_z, psihath_z
      REAL(KIND(1d0)), DIMENSION(nz):: zarray
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineURSL ! wind speed array [m s-1]
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineTRSL ! Temperature array [C]
      REAL(KIND(1d0)), DIMENSION(nz):: dataoutLineqRSL ! Specific humidity array [g kg-1]

      REAL(KIND(1d0))::z0_RSL  ! roughness length from H&F
      REAL(KIND(1d0))::zd_RSL ! zero-plane displacement

      ! REAL(KIND(1d0))::Lc_build, Lc_tree, Lc ! canopy drag length scale
      REAL(KIND(1d0))::Lc ! canopy drag length scale
      REAL(KIND(1d0))::Scc ! Schmidt number for temperature and humidity
      REAL(KIND(1d0))::psimz, psimz0, psimza, phimzp, phimz, phihzp, phihz, psihz, psihza  ! stability function for momentum
      ! REAL(KIND(1d0))::betaHF, betaNL, beta, betaN2  ! beta coefficient from Harman 2012
      REAL(KIND(1d0))::beta  ! beta coefficient from Harman 2012
      REAL(KIND(1d0))::elm ! mixing length
      ! REAL(KIND(1d0))::xxm1, xxm1_2, xxh1, xxh1_2, dphi, dphih ! dummy variables for stability functions
      REAL(KIND(1d0))::f, cm, c2, ch, c2h ! H&F'07 and H&F'08 'constants'
      REAL(KIND(1d0))::t_h, q_h ! H&F'08 canopy corrections
      REAL(KIND(1d0))::TStar_RSL ! temperature scale
      REAL(KIND(1d0))::UStar_RSL ! friction velocity used in RSL
      REAL(KIND(1d0))::PAI ! plan area index, including areas of roughness elements: buildings and trees
      ! REAL(KIND(1d0))::sfr_tr ! land cover fraction of trees
      ! REAL(KIND(1d0))::L_MOD ! Obukhov length used in RSL module with thresholds applied
      REAL(KIND(1d0))::zH_RSL ! mean canyon height used in RSL module with thresholds applied
      REAL(KIND(1d0))::dz! initial height step
      REAL(KIND(1d0))::phi_hatmZh, phim_zh
      ! REAL(KIND(1d0)), parameter::zH_min = 8! limit for minimum canyon height used in RSL module
      REAL(KIND(1d0)), parameter::ratio_dz = 1.618! ratio between neighbouring height steps

      REAL(KIND(1d0))::qa_gkg, qStar_RSL ! specific humidity scale
      INTEGER :: I, z, idx_can, idx_za, idx_2m, idx_10m
      INTEGER :: nz_can ! number of heights in canyon

      LOGICAL:: flag_RSL ! whether RSL correction is used

      ! CHARACTER(len=1024) :: Errmessage
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

      call RSL_cal_prms( &
         StabilityMethod, &!input
         zh, L_MOD, sfr, FAI, &!input
         zH_RSL, Lc, beta, zd_RSL, z0_RSL, elm, Scc, f, PAI)

      ! Define the height array with consideration of key heights
      ! set number of heights within canopy
      IF (Zh_RSL <= 2) THEN
         nz_can = 3
      ELSE IF (Zh_RSL <= 10) THEN
         nz_can = 10
      else
         nz_can = 15
      ENDIF
      ! fill up heights in canopy
      dz = Zh_RSL/nz_can
      do i = 1, nz_can
         zarray(i) = dz*i
      end do

      ! guaranttee 2 m is within the zarray
      if (dz > 2) zarray(1) = 2.

      zarray(nz_can) = Zh_RSL
      ! fill up heights above canopy
      dz = (zMeas - Zh_RSL)/(nz - nz_can)
      do i = nz_can + 1, nz
         zarray(i) = Zh_RSL + (i - nz_can)*dz
      end do

      ! add key heights (2m and 10m) to zarray
      ! 2m:
      ! DO z = 1, nz
      !    dif(z) = ABS(zarray(z) - 2)
      ! ENDDO
      ! idx_2m = MINLOC(dif, DIM=1)
      ! zarray(idx_2m) = 2
      idx_2m = 2
      ! 10m:
      ! DO z = 1, nz
      !    dif(z) = ABS(zarray(z) - 10)
      ! ENDDO
      ! idx_10m = MINLOC(dif, DIM=1)
      ! zarray(idx_10m) = 10
      idx_10m = 4

      ! determine index at the canyon top
      DO z = 1, nz
         dif(z) = ABS(zarray(z) - Zh_RSL)
      ENDDO
      idx_can = MINLOC(dif, DIM=1)
      zarray(idx_can) = Zh_RSL

      ! determine index at measurement height
      ! DO z = 1, nz
      !    dif(z) = ABS(zarray(z) - zMeas)
      ! ENDDO
      idx_za = nz
      zarray(idx_za) = zMeas

      ! see Fig 1 of Grimmond and Oke (1999) for the range for 'real cities'
      ! PAI ~ [0.1,.61], FAI ~ [0.05,0.45]
      flag_RSL = (1.-PAI)/FAI <= 18 .and. (1.-PAI)/FAI>.87

      if (flag_RSL) then
         ! use RSL approach to calculate correction factors
         ! Step 3: calculate the stability dependent H&F constants

         call cal_ch(StabilityMethod, zh_RSL, zd_RSL, Lc, beta, L_MOD, Scc, f, c2h, ch)
         call cal_cm(StabilityMethod, zH_RSL, zd_RSL, Lc, beta, L_MOD, c2, cm, phi_hatmZh, phim_zh)

         ! Step 4: determine psihat at levels above the canopy
         psihatm_z = 0.*zarray
         psihath_z = 0.*zarray
         DO z = nz - 1, idx_can, -1
            phimz = stab_phi_mom(StabilityMethod, (zarray(z) - zd_RSL)/L_MOD)
            phimzp = stab_phi_mom(StabilityMethod, (zarray(z + 1) - zd_RSL)/L_MOD)
            phihz = stab_phi_heat(StabilityMethod, (zarray(z) - zd_RSL)/L_MOD)
            phihzp = stab_phi_heat(StabilityMethod, (zarray(z + 1) - zd_RSL)/L_MOD)

            psihatm_z(z) = psihatm_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - zd_RSL)/elm)) &  !Taylor's approximation for integral
                           /(zarray(z + 1) - zd_RSL)
            psihatm_z(z) = psihatm_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - zd_RSL)/elm)) &
                           /(zarray(z) - zd_RSL)
            psihath_z(z) = psihath_z(z + 1) + dz/2.*phihzp*(ch*EXP(-1.*c2h*beta*(zarray(z + 1) - zd_RSL)/elm)) &  !Taylor's approximation for integral
                           /(zarray(z + 1) - zd_RSL)
            psihath_z(z) = psihath_z(z) + dz/2.*phihz*(ch*EXP(-1.*c2h*beta*(zarray(z) - zd_RSL)/elm)) &
                           /(zarray(z) - zd_RSL)
         ENDDO

      else

         ! correct parameters if RSL approach doesn't apply for scenario of isolated flows
         ! see Fig 1 of Grimmond and Oke (1999)
         ! when isolated flow or skimming flow, implying RSL doesn't apply, force RSL correction to zero
         psihatm_z = 0
         psihath_z = 0
         ! beta = 1.e6

         !correct RSL-based using SUEWS system-wide values
         z0_RSL = z0m
         zd_RSL = zdm
         if (zh_rsl <= zd_RSL) then
            ! this may happen as only building height is considered in calculation of zd
            zd_RSL = 0.99*zh_rsl
         end if

         ! correct elm uisng suggested valid thresholds by Harman and Finnigan (2007)

         if (L_MOD>0) then
            ! eqn 25 in HF07, stable condition:
            Lc=2.2*kappa/beta*L_MOD
         else
            ! eqn 26 in HF07, stable condition:
            Lc=-2/beta**2*L_MOD
         endif
         elm=cal_elm_RSL(beta, Lc)

         ! then MOST recovers from RSL correction
      end if

      ! Step 6: Calculate mean variables above canopy
      !
      psimz0 = stab_psi_mom(StabilityMethod, z0_RSL/L_MOD)
      psimza = stab_psi_mom(StabilityMethod, (zMeas - zd_RSL)/L_MOD)
      psihza = stab_psi_heat(StabilityMethod, (zMeas - zd_RSL)/L_MOD)
      UStar_RSL = avU1*kappa/(LOG((zMeas - zd_RSL)/z0_RSL) - psimza + psimz0 + psihatm_z(nz))

      ! set a lower limit for ustar to improve numeric stability
      UStar_RSL = merge(UStar_RSL, 0.15d0, UStar_RSL > 0.15d0)

      TStar_RSL = -1.*(qh/(avcp*avdens))/UStar_RSL
      qStar_RSL = -1.*(qe/lv_J_kg*avdens)/UStar_RSL
      qa_gkg = RH2qa(avRH/100, Press_hPa, Temp_c)

      DO z = idx_can, nz
         psimz = stab_psi_mom(StabilityMethod, (zarray(z) - zd_RSL)/L_MOD)
         psihz = stab_psi_heat(StabilityMethod, (zarray(z) - zd_RSL)/L_MOD)
         dataoutLineURSL(z) = (LOG((zarray(z) - zd_RSL)/z0_RSL) - psimz + psimz0 + psihatm_z(z))/kappa ! eqn. 3 in Theeuwes et al. (2019 BLM)
         dataoutLineTRSL(z) = (LOG((zarray(z) - zd_RSL)/(zMeas - zd_RSL)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za))/kappa ! eqn. 4 in Theeuwes et al. (2019 BLM)
         dataoutLineqRSL(z) = (LOG((zarray(z) - zd_RSL)/(zMeas - zd_RSL)) - psihz + psihza + psihath_z(z) - psihath_z(idx_za))/kappa
      ENDDO
      !
      ! Step 7: calculate in canopy variables
      !
      if (idx_can > 1) then
         t_h = Scc*TStar_RSL/(beta*f)
         q_h = Scc*qStar_RSL/(beta*f)
         DO z = 1, idx_can - 1
            dataoutLineURSL(z) = dataoutLineURSL(idx_can)*EXP(beta*(zarray(z) - Zh_RSL)/elm)
            dataoutLineTRSL(z) = dataoutLineTRSL(idx_can) + (t_h*EXP(beta*f*(zarray(z) - Zh_RSL)/elm) - t_h)/TStar_RSL
            dataoutLineqRSL(z) = dataoutLineqRSL(idx_can) + (q_h*EXP(beta*f*(zarray(z) - Zh_RSL)/elm) - q_h)/qStar_RSL
         ENDDO
      end if

      dataoutLineURSL = dataoutLineURSL*UStar_RSL
      dataoutLineTRSL = dataoutLineTRSL*TStar_RSL + Temp_C
      dataoutLineqRSL = (dataoutLineqRSL*qStar_RSL + qa_gkg/1000.)*1000.

      dataoutLineRSL = [zarray, dataoutLineURSL, dataoutLineTRSL, dataoutLineqRSL, &
                        !information for debugging
                        L_MOD, zH_RSL, Lc, beta, zd_RSL, z0_RSL, elm, Scc, f, UStar_RSL, FAI, PAI, merge(1.d0, 0.d0, flag_RSL) &
                        ]

      !
      ! Step 8
      ! retrieve the diagnostics at key heights
      !
      T2_C = interp_z(2d0, zarray, dataoutLineTRSL)
      q2_gkg = interp_z(2d0, zarray, dataoutLineqRSL)
      U10_ms = interp_z(10d0, zarray, dataoutLineURSL)
      ! get relative humidity:
      RH2 = qa2RH(q2_gkg, press_hPa, T2_C)

   END SUBROUTINE RSLProfile

   function interp_z(z_x, z, v) result(v_x)

      real(KIND(1D0)), intent(in) ::  z_x ! height to interpolate at
      real(KIND(1D0)), dimension(nz), intent(in) ::  z ! heights
      real(KIND(1D0)), dimension(nz), intent(in) ::  v ! values associated with heights

      ! output
      real(KIND(1D0)) ::v_x ! zd used in RSL

      ! local variables
      real(KIND(1D0)) ::slope! slope
      real(KIND(1D0)) ::dz! slope
      real(KIND(1D0)), dimension(nz) ::dif! slope
      integer :: idx_low! vertical index lower than z_x
      integer :: idx_x! vertical index lower than z_x
      integer :: idx_high! vertical index higher than z_x
      integer :: idx! vertical index higher than z_x
      integer, PARAMETER::nz = 30! vertical index higher than z_x

      ! initialise variables
      idx_x = -999

      dif = z - z_x
      idx_x = MAXLOC(dif, 1, dif == 0)
      idx_low = MAXLOC(dif, 1, dif < 0)
      idx_high = MINLOC(dif, 1, dif > 0)

      if (idx > 0) then
         ! z_x is one of zarray elements
         v_x = v(idx_x)
      else
         ! linear interpolation is performed
         dz = z(idx_high) - z(idx_low)
         slope = (v(idx_high) - v(idx_low))/dz
         v_x = v(idx_low) + (z_x - z(idx_low))*slope
      endif

   end function interp_z

   function cal_elm_RSL(beta, Lc) result(elm)

      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL

      ! output
      real(KIND(1D0)) ::elm ! zd used in RSL

      elm = 2.*beta**3*Lc

   end function cal_elm_RSL

   recursive function cal_psim_hat(StabilityMethod, z, zh_RSL, L_MOD, beta, Lc) result(psim_hat_z)
      ! calculate psi_hat for momentum
      ! TS, 23 Oct 2019
      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) :: z ! height of interest [m]
      real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
      real(KIND(1D0)), intent(in) ::  L_MOD ! Obukhov length [m]

      ! output
      real(KIND(1D0)) ::psim_hat_z ! psim_hat at height of interest

      ! internal variables
      real(KIND(1D0)) ::zp ! a height above z used for iterative calculations
      real(KIND(1D0)) ::zd_RSL ! displacement height used in RSL
      real(KIND(1D0)) ::phim_lc ! displacement height used in RSL
      real(KIND(1D0)) ::phim_z ! displacement height used in RSL
      real(KIND(1D0)) ::phim_zp ! displacement height used in RSL
      real(KIND(1D0)) ::phim_hat_zp ! displacement height used in RSL
      real(KIND(1D0)) ::phim_hat_z ! displacement height used in RSL
      real(KIND(1D0)) ::psim_hat_zp ! displacement height used in RSL
      real(KIND(1D0)) ::elm ! displacement height used in RSL
      ! real(KIND(1D0)) ::xxm1 ! displacement height used in RSL
      ! real(KIND(1D0)) ::xxm1_2 ! displacement height used in RSL
      ! real(KIND(1D0)) ::dphi ! displacement height used in RSL
      ! real(KIND(1D0)) ::phi_hatmZh ! displacement height used in RSL
      ! real(KIND(1D0)) ::cm
      ! real(KIND(1D0)) ::c2
      ! real(KIND(1D0)) ::phi_hatmZh, phim_zh

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

      if (z > 100) then
         psim_hat_z = 0.
         return
      end if

      zp = 1.01*z ! a height above z

      zd_RSL = cal_zd_RSL(Zh_RSL, beta, lc)
      elm = cal_elm_RSL(beta, lc)

      ! phim at Lc
      phim_lc = stab_phi_mom(StabilityMethod, Lc/L_MOD)

      phim_z = stab_phi_mom(StabilityMethod, (z - zd_RSL)/L_MOD)
      phim_zp = stab_phi_mom(StabilityMethod, (zp - zd_RSL)/L_MOD)

      psim_hat_zp = cal_psim_hat(StabilityMethod, zp, zh_RSL, L_MOD, beta, Lc)

      !Taylor's approximation for integral
      ! psim_hat_z = psim_hat_zp + dz/2.*phim_zp*(cm*EXP(-1.*c2*beta*(zp - zd_RSL)/elm))/(zp - zd_RSL)
      ! psim_hat_z = psim_hat_z + dz/2.*phim_z*(cm*EXP(-1.*c2*beta*(z - zd_RSL)/elm))/(z - zd_RSL)
      phim_hat_zp = cal_phim_hat(StabilityMethod, zp, zh_RSL, L_MOD, beta, lc)
      phim_hat_z = cal_phim_hat(StabilityMethod, z, zh_RSL, L_MOD, beta, lc)

      psim_hat_z = psim_hat_zp + dz/2.*phim_zp*(1 - phim_hat_zp)/(zp - zd_RSL)
      psim_hat_z = psim_hat_z + dz/2.*phim_z*(1 - phim_hat_z)/(z - zd_RSL)

   end function cal_psim_hat

   function cal_phim_hat(StabilityMethod, z, zh_RSL, L_MOD, beta, lc) result(phim_hat)
      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) :: z
      real(KIND(1D0)), intent(in) :: zh_RSL
      real(KIND(1D0)), intent(in) :: L_MOD
      real(KIND(1D0)), intent(in) :: beta
      real(KIND(1D0)), intent(in) :: lc
      real(KIND(1D0)) :: phim_hat
      real(KIND(1D0)) :: zd_RSL

      real(KIND(1D0)) :: elm
      real(KIND(1D0)) :: c2
      real(KIND(1D0)) :: cm, phi_hatmZh, phim_zh

      elm = cal_elm_RSL(beta, lc)

      zd_RSL = cal_zd_RSL(Zh_RSL, beta, lc)

      call cal_cm(StabilityMethod, zh_RSL, zd_RSL, Lc, beta, L_MOD, c2, cm, phi_hatmZh, phim_zh)

      phim_hat = 1 - cm*EXP(-1.*c2*beta*(z - zd_RSL)/elm)

   end function cal_phim_hat

   subroutine cal_cm(StabilityMethod, zh_RSL, zd_RSL, Lc, beta, L_MOD, c2, cm, phi_hatmZh, phim_zh)

      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      ! real(KIND(1D0)), intent(in) :: z ! height of interest [m]
      real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  zd_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
      real(KIND(1D0)), intent(in) ::  L_MOD ! Obukhov length [m]

      ! output
      real(KIND(1D0)), intent(out) ::c2
      real(KIND(1D0)), intent(out) ::cm
      real(KIND(1D0)), intent(out) ::phi_hatmZh
      real(KIND(1D0)), intent(out) ::phim_zh

      ! internal variables
      ! real(KIND(1D0)) ::phim_zh
      real(KIND(1D0)) ::phim_zhdz
      real(KIND(1D0)) ::dphi
      ! real(KIND(1D0)) ::phi_hatmZh

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

      phim_zh = stab_phi_mom(StabilityMethod, (Zh_RSL - zd_RSL)/L_MOD)
      phim_zhdz = stab_phi_mom(StabilityMethod, (Zh_RSL - zd_RSL + dz)/L_MOD)

      dphi = (phim_zhdz - phim_zh)/dz
      if (phim_zh /= 0.) then
         phi_hatmZh = kappa/(2.*beta*phim_zh)
      else
         ! neutral condition
         phi_hatmZh = 1.
      end if

      IF (phi_hatmZh >= 1.) THEN
         ! more stable, but less correct
         c2 = 0.5
         phi_hatmZh = 1.
      ELSE
         ! if very unstable this might cause some high values of psihat_z
         c2 = (kappa*(3.-(2.*beta**2.*Lc/phim_zh*dphi)))/(2.*beta*phim_zh - kappa)
      ENDIF
      ! force c2 to 0.5 for better stability. TS 14 Jul 2020
      ! TODO: a more proper threshold needs to be determined
      c2 = 0.5

      cm = (1.-phi_hatmZh)*EXP(c2/2.)

   end subroutine cal_cm

   subroutine cal_ch(StabilityMethod, zh_RSL, zd_RSL, Lc, beta, L_MOD, Scc, f, c2h, ch)

      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      ! real(KIND(1D0)), intent(in) :: z ! height of interest [m]
      real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  zd_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  Scc !
      real(KIND(1D0)), intent(in) ::  f !
      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
      real(KIND(1D0)), intent(in) ::  L_MOD ! Obukhov length [m]

      ! output
      real(KIND(1D0)), intent(out) ::ch
      real(KIND(1D0)), intent(out) ::c2h ! displacement height used in RSL

      ! internal variables
      real(KIND(1D0)) ::phih_zh ! displacement height used in RSL
      real(KIND(1D0)) ::phih_zhdz ! displacement height used in RSL
      real(KIND(1D0)) ::dphih ! displacement height used in RSL
      real(KIND(1D0)) ::phi_hathZh ! displacement height used in RSL

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

      phih_zh = stab_phi_heat(StabilityMethod, (Zh_RSL - zd_RSL)/L_MOD)
      phih_zhdz = stab_phi_heat(StabilityMethod, (Zh_RSL - zd_RSL + 1.)/L_MOD)

      dphih = phih_zhdz - phih_zh
      if (phih_zh /= 0.) then
         phi_hathZh = kappa*Scc/(2.*beta*phih_zh)
      else
         phi_hathZh = 1.
      end if

      IF (phi_hathZh >= 1.) THEN
         ! more stable, but less correct
         c2h = 0.5
         phi_hathZh = 1.
      ELSE
         ! if very unstable this might cause some high values of psihat_z
         c2h = (kappa*Scc*(2.+f - (dphih*2.*beta**2.*Lc/phih_zh)))/(2.*beta*phih_zh - kappa*Scc)
      ENDIF
      ! force c2h to 0.5 for better stability. TS 14 Jul 2020
      ! TODO: a more proper threshold needs to be determined
      c2h = 0.5

      ch = (1.-phi_hathZh)*EXP(c2h/2.)

   end subroutine cal_ch

   ! function cal_psihatm_z(StabilityMethod, nz, zarray, L_MOD_RSL, zH_RSL, Lc, beta, zd, elm) result(psihatm_z)

   !    ! calculate psi_hat for momentum
   !    ! TS, 23 Oct 2019
   !    implicit none
   !    integer, intent(in) :: StabilityMethod ! stability method
   !    integer, intent(in) :: nz ! number of vertical layers
   !    real(KIND(1D0)), DIMENSION(nz), intent(in) :: zarray ! height of interest [m]
   !    real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
   !    real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
   !    real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
   !    real(KIND(1D0)), intent(in) ::  L_MOD_RSL ! Obukhov length [m]
   !    real(KIND(1D0)), intent(in) ::zd ! displacement height used in RSL
   !    real(KIND(1D0)), intent(in) ::elm ! displacement height used in RSL

   !    ! output
   !    real(KIND(1D0)), DIMENSION(nz) ::psihatm_z ! psim_hat at height of interest

   !    ! internal variables
   !    ! real(KIND(1D0)) ::zp ! a height above z used for iterative calculations
   !    ! real(KIND(1D0)) ::phim_lc ! displacement height used in RSL
   !    real(KIND(1D0)) ::phimz ! displacement height used in RSL
   !    real(KIND(1D0)) ::phimzp ! displacement height used in RSL
   !    ! real(KIND(1D0)) ::psim_hat_zp ! displacement height used in RSL
   !    real(KIND(1D0)) ::xxm1 ! displacement height used in RSL
   !    real(KIND(1D0)) ::xxm1_2 ! displacement height used in RSL
   !    real(KIND(1D0)) ::dphi ! displacement height used in RSL
   !    real(KIND(1D0)) ::phi_hatmZh ! displacement height used in RSL
   !    real(KIND(1D0)) ::cm ! displacement height used in RSL
   !    real(KIND(1D0)) ::c2 ! displacement height used in RSL
   !    REAL(KIND(1d0)), DIMENSION(nz):: dif

   !    REAL(KIND(1d0)), PARAMETER::kappa = 0.40
   !    REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

   !    integer::z, idx_can

   !    psihatm_z = 0.*zarray

   !    ! determine index at the canyon top
   !    DO z = 1, nz
   !       dif(z) = ABS(zarray(z) - Zh_RSL)
   !    ENDDO
   !    idx_can = MINLOC(dif, DIM=1)
   !    ! zarray(idx_can) = Zh_RSL

   !    ! calculate phihatM according to H&F '07 and H&F '08 for heat and humidity
   !    xxm1 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)
   !    xxm1_2 = stab_phi_mom(StabilityMethod, (Zh_RSL - zd + 1.)/L_MOD_RSL)

   !    phi_hatmZh = kappa/(2.*beta*xxm1)
   !    dphi = xxm1_2 - xxm1

   !    IF (phi_hatmZh > 1.) THEN
   !       c2 = 0.5 ! more stable, but less correct
   !       ! c2h = 0.5
   !    ELSE
   !       c2 = (kappa*(3.-(2.*beta**2.*Lc/xxm1*dphi)))/(2.*beta*xxm1 - kappa)  ! if very unstable this might cause some high values of psihat_z
   !       ! c2h = (kappa*Scc*(2.+f - (dphih*2.*beta**2.*Lc/xxh1)))/(2.*beta*xxh1 - kappa*Scc)
   !    ENDIF
   !    cm = (1.-phi_hatmZh)*EXP(c2/2.)
   !    ! ch = (1.-phi_hathzh)*EXP(c2h/2.)

   !    DO z = nz - 1, idx_can - 1, -1
   !       phimz = stab_phi_mom(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
   !       phimzp = stab_phi_mom(StabilityMethod, (zarray(z + 1) - zd)/L_MOD_RSL)

   !       psihatm_z(z) = psihatm_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
   !                      /(zarray(z + 1) - zd)
   !       psihatm_z(z) = psihatm_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - zd)/elm)) &
   !                      /(zarray(z) - zd)

   !    ENDDO

   ! end function cal_psihatm_z

   ! function cal_psihath_z(StabilityMethod, nz, zarray, L_MOD_RSL, zH_RSL, Lc, beta, zd, elm, Scc, f) result(psihath_z)

   !    ! calculate psi_hat for momentum
   !    ! TS, 23 Oct 2019
   !    implicit none
   !    integer, intent(in) :: StabilityMethod ! stability method
   !    integer, intent(in) :: nz ! number of vertical layers

   !    real(KIND(1D0)), DIMENSION(nz), intent(in) :: zarray ! height of interest [m]
   !    real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
   !    real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
   !    real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL
   !    real(KIND(1D0)), intent(in) ::  Scc ! parameter in RSL
   !    real(KIND(1D0)), intent(in) ::  f ! parameter in RSL
   !    real(KIND(1D0)), intent(in) ::  L_MOD_RSL ! Obukhov length [m]
   !    real(KIND(1D0)), intent(in) ::  elm ! displacement height used in RSL
   !    real(KIND(1D0)), intent(in) ::zd ! displacement height used in RSL

   !    ! output
   !    real(KIND(1D0)), DIMENSION(nz) ::psihath_z ! psim_hat at height of interest

   !    ! internal variables
   !    ! real(KIND(1D0)) ::zp ! a height above z used for iterative calculations
   !    ! real(KIND(1D0)) ::phim_lc ! displacement height used in RSL
   !    real(KIND(1D0)) ::phihz ! displacement height used in RSL
   !    real(KIND(1D0)) ::phihzp ! displacement height used in RSL
   !    ! real(KIND(1D0)) ::psim_hat_zp ! displacement height used in RSL
   !    real(KIND(1D0)) ::xxh1 ! displacement height used in RSL
   !    real(KIND(1D0)) ::xxh1_2 ! displacement height used in RSL
   !    real(KIND(1D0)) ::dphih ! displacement height used in RSL
   !    real(KIND(1D0)) ::phi_hathZh ! displacement height used in RSL
   !    real(KIND(1D0)) ::ch ! displacement height used in RSL
   !    real(KIND(1D0)) ::c2h ! displacement height used in RSL
   !    REAL(KIND(1d0)), DIMENSION(nz):: dif

   !    REAL(KIND(1d0)), PARAMETER::kappa = 0.40
   !    REAL(KIND(1d0)), PARAMETER::dz = 0.1 !height step

   !    integer::z, idx_can

   !    psihath_z = 0.*zarray

   !    ! determine index at the canyon top
   !    DO z = 1, nz
   !       dif(z) = ABS(zarray(z) - Zh_RSL)
   !    ENDDO
   !    idx_can = MINLOC(dif, DIM=1)
   !    ! zarray(idx_can) = Zh_RSL

   !    ! calculate phihatM according to H&F '07 and H&F '08 for heat and humidity
   !    xxh1 = stab_phi_heat(StabilityMethod, (Zh_RSL - zd)/L_MOD_RSL)
   !    xxh1_2 = stab_phi_heat(StabilityMethod, (Zh_RSL - zd + 1.)/L_MOD_RSL)

   !    phi_hathZh = kappa*Scc/(2.*beta*xxh1)
   !    dphih = xxh1_2 - xxh1

   !    IF (phi_hathZh > 1.) THEN
   !       ! c2 = 0.5 ! more stable, but less correct
   !       c2h = 0.5
   !    ELSE
   !       ! c2 = (kappa*(3.-(2.*beta**2.*Lc/xxm1*dphi)))/(2.*beta*xxm1 - kappa)  ! if very unstable this might cause some high values of psihat_z
   !       c2h = (kappa*Scc*(2.+f - (dphih*2.*beta**2.*Lc/xxh1)))/(2.*beta*xxh1 - kappa*Scc)
   !    ENDIF
   !    ! cm = (1.-phi_hatmZh)*EXP(c2/2.)
   !    ch = (1.-phi_hathzh)*EXP(c2h/2.)

   !    DO z = nz - 1, idx_can - 1, -1
   !       phihz = stab_phi_heat(StabilityMethod, (zarray(z) - zd)/L_MOD_RSL)
   !       phihzp = stab_phi_heat(StabilityMethod, (zarray(z + 1) - zd)/L_MOD_RSL)

   !       psihath_z(z) = psihath_z(z + 1) + dz/2.*phihzp*(ch*EXP(-1.*c2h*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
   !                      /(zarray(z + 1) - zd)
   !       psihath_z(z) = psihath_z(z) + dz/2.*phihz*(ch*EXP(-1.*c2h*beta*(zarray(z) - zd)/elm)) &
   !                      /(zarray(z) - zd)

   !    ENDDO

   ! end function cal_psihath_z

   function cal_zd_RSL(zh_RSL, beta, Lc) result(zd_RSL)

      real(KIND(1D0)), intent(in) ::  zh_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(in) ::  beta ! parameter in RSL

      ! output
      real(KIND(1D0)) ::zd_RSL ! zd used in RSL

      zd_RSL = Zh_RSL - (beta**2.)*Lc
      !correct negative values using rule of thumb, TS 24 Jun 2020
      ! if (zd_RSL < 0) zd_RSL = 0.7*Zh_RSL

   end function cal_zd_RSL

   function cal_z0_RSL(StabilityMethod, zH_RSL, zd_RSL, beta, L_MOD_RSL, Lc) result(z0_RSL)
      ! calculate z0 iteratively
      ! TS, 23 Oct 2019
      implicit none
      integer, intent(in) ::StabilityMethod
      real(KIND(1D0)), intent(in) ::  zH_RSL ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  zd_RSL ! displacement height [m]
      real(KIND(1D0)), intent(in) ::  L_MOD_RSL ! Monin Obukhov length[m]
      real(KIND(1D0)), intent(in) ::  Lc ! canyon length scale [m]
      real(KIND(1D0)), intent(in) ::  beta ! height scale for bluff bodies [m]

      ! output
      real(KIND(1D0)) ::z0_RSL

      ! internal variables
      real(KIND(1D0)) ::psimZh, psimz0, z0_RSL_x, psihatm_Zh
      real(KIND(1D0)) ::err
      integer ::it

      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      ! REAL(KIND(1d0)), PARAMETER::r = 0.1
      ! REAL(KIND(1d0)), PARAMETER::a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1.

      psimZh = stab_psi_mom(StabilityMethod, (Zh_RSL - zd_RSL)/L_MOD_RSL)
      psihatm_Zh = cal_psim_hat(StabilityMethod, Zh_RSL, zh_RSL, L_MOD_RSL, beta, Lc)

      !first guess
      z0_RSL = 0.1*Zh_RSL
      err = 10.
      ! psimz0 = 0.5
      it = 1
      DO WHILE ((err > 0.001) .AND. (it < 10))
         z0_RSL_x = z0_RSL
         psimz0 = stab_psi_mom(StabilityMethod, z0_RSL_x/L_MOD_RSL)
         z0_RSL = (Zh_RSL - zd_RSL)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0)*EXP(psihatm_Zh)
         err = ABS(z0_RSL_x - z0_RSL)
         it = it + 1
      ENDDO

      ! set limit on z0_RSL for numeric stability
      z0_RSL=merge(z0_RSL,1d-2,z0_RSL<1d-2)

   end function cal_z0_RSL

   subroutine RSL_cal_prms( &
      StabilityMethod, zh, L_MOD, sfr, FAI, &!input
      zH_RSL, Lc, beta, zd_RSL, z0_RSL, elm, Scc, f, PAI)!output

      implicit none
      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) ::  zh ! canyon depth [m]
      real(KIND(1D0)), intent(in) ::  FAI ! frontal area index
      real(KIND(1D0)), intent(in) ::  L_MOD ! Obukhov length [m]
      real(KIND(1D0)), DIMENSION(nsurf), intent(in) ::  sfr ! land cover fractions

      ! output
      ! real(KIND(1D0)), intent(out) ::L_MOD ! Obukhov length used in RSL module with thresholds applied
      real(KIND(1D0)), intent(out) ::zH_RSL ! mean canyon height used in RSL module with thresholds applied
      real(KIND(1D0)), intent(out) ::Lc ! height scale for bluff bodies [m]
      real(KIND(1D0)), intent(out) ::beta ! psim_hat at height of interest
      real(KIND(1D0)), intent(out) ::zd_RSL ! displacement height to prescribe if necessary [m]
      real(KIND(1D0)), intent(out) ::z0_RSL ! roughness length [m]
      real(KIND(1D0)), intent(out) ::elm ! length scale used in RSL
      real(KIND(1D0)), intent(out) ::Scc ! parameter in RSL
      real(KIND(1D0)), intent(out) ::f ! parameter in RSL
      real(KIND(1D0)), intent(out) ::PAI ! plan area index inlcuding area of trees

      ! internal variables
      real(KIND(1D0)) ::sfr_tr
      real(KIND(1D0)) ::lc_over_L
      real(KIND(1D0)) ::betaHF
      real(KIND(1D0)) ::betaNL

      REAL(KIND(1d0)), PARAMETER::planF_low = 1e-6
      REAL(KIND(1d0)), PARAMETER::kappa = 0.40
      ! REAL(KIND(1d0)), PARAMETER::z0m= 0.40
      REAL(KIND(1d0)), PARAMETER::r = 0.1
      REAL(KIND(1d0)), PARAMETER::a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1.
      REAL(KIND(1d0)), parameter::Zh_min = 0.4! limit for minimum canyon height used in RSL module

      ! under stable conditions, set a threshold for L_MOD to avoid numerical issues. TS 28 Oct 2019
      ! L_MOD = merge(L_MOD, 300.d1, L_MOD < 300.)

      ! zH_RSL
      zH_RSL = max(zh, Zh_min)

      ! land cover fraction of bluff bodies
      PAI = sum(sfr([BldgSurf, ConifSurf, DecidSurf]))
      ! set a threshold for sfr_zh to avoid numerical difficulties
      PAI = min(PAI, 0.8)

      ! land cover fraction of trees
      sfr_tr = sum(sfr([ConifSurf, DecidSurf]))

      ! height scale for buildings !not used? why?
      ! Lc_build = (1.-sfr(BldgSurf))/FAI*Zh_RSL  ! Coceal and Belcher 2004 assuming Cd = 2

      ! height scale for tress
      ! Lc_tree = 1./(cd_tree*a_tree) ! not used? why?

      ! height scale for bluff bodies
      Lc = (1.-PAI)/FAI*Zh_RSL

      ! a normalised scale with a physcially valid range between [-2,2] (Harman 2012, BLM)
      lc_over_L = lc/L_MOD
      if (lc_over_L > 0) then
         lc_over_L = min(2., lc_over_L)
      else
         lc_over_L = max(-2., lc_over_L)
      end if

      ! Schmidt number Harman and Finnigan 2008: assuming the same for heat and momemntum
      Scc = 0.5 + 0.3*TANH(2.*lc_over_L)
      f = 0.5*((1.+4.*r*Scc)**0.5) - 0.5

      ! Step 2:
      ! Parameterise beta according to Harman 2012 with upper limit of 0.5

      call cal_beta_RSL(StabilityMethod, PAI, sfr_tr, lc_over_L, beta, betaHF, betaNL)

      zd_RSL = cal_zd_RSL(Zh_RSL, beta, lc)

      elm = cal_elm_RSL(beta, Lc)

      ! calculate z0 iteratively
      z0_RSL = cal_z0_RSL(StabilityMethod, zh_RSL, zd_RSL, beta, L_MOD, lc)

   end subroutine RSL_cal_prms

   subroutine cal_beta_RSL(StabilityMethod, PAI, sfr_tr, lc_over_L, beta, betaHF, betaNL)
      ! Step 2:
      ! Parameterise beta according to Harman 2012 with upper limit of 0.5
      implicit none

      integer, intent(in) :: StabilityMethod ! stability method
      real(KIND(1D0)), intent(in) :: PAI
      real(KIND(1D0)), intent(in) :: sfr_tr
      real(KIND(1D0)), intent(in) :: lc_over_L

      ! output
      real(KIND(1D0)), intent(out):: beta
      real(KIND(1D0)), intent(out):: betaHF
      real(KIND(1D0)), intent(out):: betaNL

      real(KIND(1D0)), PARAMETER :: kappa = 0.4
      REAL(KIND(1d0)), PARAMETER::a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1.
      ! real(KIND(1D0)) :: phim_hat
      ! real(KIND(1D0)) :: zd_RSL

      real(KIND(1D0)) :: betaN2, betaHF_0, betaNL_0
      INTEGER :: it
      real(KIND(1D0)) :: err
      real(KIND(1D0)) :: phim

      !

      ! betaN for trees found to be 0.3 and for urban 0.4 linearly interpolate between the two using surface fractions
      ! betaN2 = 0.30 + (1.-sfr(ConifSurf) - sfr(ConifSurf))*0.1
      if (PAI > 0) then
         betaN2 = 0.30*sfr_tr/PAI + (PAI - sfr_tr)/PAI*0.4
      ELSE
         betaN2 = 0.35
      endif

      ! betaHF
      ! it = 1
      ! phim = 1
      ! DO WHILE ((err > 0.001) .AND. (it < 10))
      !    betaHF_0 = betaN2/phim
      !    phim = stab_phi_mom(StabilityMethod, (betaHF_0**2)*Lc/L_MOD)
      !    betaHF = betaN2/phim
      !    err = ABS(betaHF - betaHF_0)
      !    it = it + 1
      ! ENDDO

      betaHF = cal_beta_lc(stabilityMethod, betaN2, lc_over_L)

      ! betaNL
      ! it = 1
      ! phim = 1
      ! DO WHILE ((err > 0.001) .AND. (it < 10))
      !    betaNL_0 = (kappa/2.)/phim
      !    phim = stab_phi_mom(StabilityMethod, (betaNL_0**2)*Lc/L_MOD)
      !    betaNL = (kappa/2.)/phim
      !    err = ABS(betaNL_0 - betaNL)
      !    it = it + 1
      ! ENDDO

      betaNL = cal_beta_lc(stabilityMethod, kappa/2., lc_over_L)

      IF (lc_over_L > a2) THEN
         beta = betaHF
      ELSE
         beta = betaNL + ((betaHF - betaNL)/(1.+a1*abs(lc_over_L - a2)**a3))
      ENDIF

      IF (beta > 0.5) THEN
         beta = 0.5
      ENDIF

   end subroutine cal_beta_RSL

   function cal_beta_lc(stabilityMethod, beta0, lc_over_l) result(beta_x)
      ! TS, 03 Aug 2020:
      ! iterative determination of beta depending on Lc/L
      ! ref: eqn 10 & 11 in Harman (2012, BLM)
      implicit none
      integer, intent(in) :: StabilityMethod
      real(KIND(1D0)), intent(in) :: beta0
      real(KIND(1D0)), intent(in) ::lc_over_l
      real(KIND(1D0)) :: beta_x

      real(KIND(1D0)) :: phim, err, beta_x0

      integer::it

      it = 1
      phim = 1
      err = 1
      ! print *, '***********************'
      ! print *, 'beta0:', beta0
      ! print *, 'Lc/L_MOD:', lc_over_l
      DO WHILE ((err > 0.001) .AND. (it < 20))
         beta_x0 = beta0/phim
         phim = stab_phi_mom(StabilityMethod, (beta_x0**2)*lc_over_l)
         beta_x = beta0/phim
         err = ABS(beta_x - beta_x0)
         ! print *, it, err, beta_x0, beta_x, phim, lc_over_l
         it = it + 1

      ENDDO
      ! print *, ''

   end function cal_beta_lc

end module rsl_module
