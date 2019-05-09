MODULE AtmMoistStab_module
   IMPLICIT NONE

CONTAINS
   !.c!! For Lumps Version 2 - no stability calculations
   ! Latent heat of sublimation when air temperature below zero added. LJ Nov 2012
   ! explict interface added to all subroutines, TS 08 Aug 2017
   SUBROUTINE LUMPS_cal_AtmMoist( &
      Temp_C, Press_hPa, avRh, dectime, &! input:
      lv_J_kg, lvS_J_kg, &! output:
      es_hPa, Ea_hPa, VPd_hpa, VPD_Pa, dq, dens_dry, avcp, air_dens)

      USE meteo, ONLY: &
         sat_vap_press_x, spec_heat_beer, &
         lat_vap, lat_vapSublim, spec_hum_def

      IMPLICIT NONE
      REAL(KIND(1d0))::vap_dens

      REAL(KIND(1d0)), INTENT(in):: &
         Temp_C, &
         Press_hPa, &
         avRh, dectime
      REAL(KIND(1d0)), INTENT(out):: &
         lv_J_kg, &!Latent heat of vaporization in [J kg-1]
         lvS_J_kg, &!Latent heat of sublimation in J/kg
         es_hPa, &!Saturation vapour pressure over water in hPa
         Ea_hPa, &!Vapour pressure of water in hPa
         VPd_hpa, & !vapour pressure deficit in hPa
         VPD_Pa, & !vapour pressure deficit in Pa
         dq, &!Specific humidity deficit in g/kg
         dens_dry, & !Vap density or absolute humidity         (kg/m3)
         avcp, &!specific heat capacity in J kg-1 K-1
         air_dens!Air density in kg/m3

      REAL(KIND(1d0)), PARAMETER:: &
         !  comp          = 0.9995, &
         !  epsil         = 0.62197,&           !ratio molecular weight of water vapor/dry air (kg/mol/kg/mol)
         !  epsil_gkg     = 621.97, &           !ratio molecular weight of water vapor/dry air in g/kg
         !  dry_gas       = 8.31451,&           !Dry gas constant (J/k/mol)
         !  gas_ct_wat    = 461.05,&            !Gas constant for water (J/kg/K)
         !  molar         = 0.028965,&          !Dry air molar fraction in kg/mol
         !  molar_wat_vap = 0.0180153,&         !Molar fraction of water vapor in kg/mol
         gas_ct_dry = 8.31451/0.028965, &  !j/kg/k=dry_gas/molar
         gas_ct_wv = 8.31451/0.0180153 !j/kg/kdry_gas/molar_wat_vap
      !  waterDens     = 999.8395            !Density of water in 0 cel deg
      INTEGER::from = 1

      !Saturation vapour pressure over water in hPa
      es_hPa = sat_vap_press_x(Temp_C, Press_hPa, from, dectime) ! dectime is more or less unnecessary here

      !Vapour pressure of water in hPa
      Ea_hPa = avRh/100*es_hPa

      ! if(debug.and.dectime>55.13.and.dectime<55.2)write(35,*)'%',Temp_C

      VPd_hpa = es_hPa - ea_hpa           !vapour pressure deficit in hPa
      VPD_Pa = (es_hPa*100) - (Ea_hPa*100)!vapour pressure deficit in Pa

      dq = spec_hum_def(vpd_hPa, Press_hPa) !Specific humidity deficit in g/kg

      !Vap density or absolute humidity         (kg/m3)
      vap_dens = (Ea_hPa*100/((Temp_C + 273.16)*gas_ct_wv))

      !density Dry Air Beer(1990)        kg/m3
      dens_dry = ((Press_hPa - Ea_hPa)*100)/(gas_ct_dry*(273.16 + Temp_C))

      !Air density in kg/m3
      air_dens = (Press_hPa*100)/(gas_ct_dry*(Temp_C + 273.16))

      !Calculate specific heat capacity in J kg-1 K-1
      avcp = spec_heat_beer(Temp_C, avRh, vap_dens, dens_dry)

      !Latent heat of vaporization in [J kg-1]
      lv_J_kg = lat_vap(Temp_C, Ea_hPa, Press_hPa, avcp, dectime)

      !Latent heat of sublimation in J/kg
      IF (Temp_C < 0.000) THEN
         lvS_J_kg = lat_vapSublim(Temp_C, Ea_hPa, Press_hPa, avcp)
      ELSE
         lvS_J_kg = 2834000
      ENDIF
      !if(debug)write(*,*)lv_J_kg,Temp_C,'lv2'
      IF (press_hPa < 900) THEN
         CALL ErrorHint(46, 'Function LUMPS_cal_AtmMoist', press_hPa, -55.55d0, -55)
      ENDIF
      RETURN
   END SUBROUTINE LUMPS_cal_AtmMoist

   !.c!! For Lumps Version 2 - no stability calculations
   !==========================================================
   !     Last change:
   !     TS   08 Aug 2017: added explicit interface
   !     TS   13 Jun 2017: corrections to the integral of stability functions
   !     MH   12 Apr 2017: Stable limit to exit do-loop
   !     LJ   25 Nov 2014: Limits for L
   !     LJ   19 Feb 2010
   !     SG   27 Mar 2000    4:44 pm
   !     ust - friction velocity
   !     L - monin obukhov stability length
   !       Van Ulden & Holtslag (1985) JCAM: 24: 1196-1207

   SUBROUTINE STAB_lumps( &
      ! input
      StabilityMethod, &
      dectime, & !Decimal time
      zzd, &     !Active measurement height (meas. height-displac. height)
      z0m, &     !Aerodynamic roughness length
      zdm, &     !Displacement height
      avU1, &    !Average wind speed
      Temp_C, &  !Air temperature
      H_init, & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
      ! output:
      L_MOD, & !Obukhov length
      TStar, & !T*
      UStar, & !Friction velocity
      zL)!Stability scale

      IMPLICIT NONE
      INTEGER, INTENT(in):: StabilityMethod

      REAL(KIND(1d0)), INTENT(in)::dectime !Decimal time
      REAL(KIND(1d0)), INTENT(in)::zzd     !Active measurement height (meas. height-displac. height)
      REAL(KIND(1d0)), INTENT(in)::z0m     !Aerodynamic roughness length
      REAL(KIND(1d0)), INTENT(in)::zdm     !Displacement height
      REAL(KIND(1d0)), INTENT(in)::avU1    !Average wind speed
      REAL(KIND(1d0)), INTENT(in)::Temp_C    !Air temperature
      REAL(KIND(1d0)), INTENT(in)::H_init    !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity

      REAL(KIND(1d0)), INTENT(out)::L_MOD!Obukhov length
      REAL(KIND(1d0)), INTENT(out)::TStar!T*
      REAL(KIND(1d0)), INTENT(out)::UStar!Friction velocity
      REAL(KIND(1d0)), INTENT(out)::zL ! Stability scale
      ! REAL(KIND(1d0)),INTENT(out)::psim   !Stability function of momentum

      REAL(KIND(1d0))::G_T_K, &
                        KUZ, &
                        LOLD, &
                        psim, &
                        z0l, &
                        psimz0, &
                        h
      REAL(KIND(1d0)), PARAMETER :: &
         k = 0.4, &             !Von Karman's contant
         grav = 9.80665  !g - gravity - physics today august 1987
      INTEGER, PARAMETER ::   notUsedI = -55

      INTEGER :: i

      LOGICAL :: debug = .FALSE.

      IF (debug) WRITE (*, *) StabilityMethod, z0m, avU1, H_init, UStar, L_MOD
      G_T_K = (Grav/(Temp_C + 273.16))*k !gravity constant/(Temperature*Von Karman Constant)
      KUZ = k*AvU1                     !Von Karman constant*mean wind speed
   IF (zzd < 0) CALL ErrorHint(32, 'Windspeed Ht too low relative to zdm [Stability calc]- values [z-zdm, zdm]', Zzd, zdm, notUsedI)

      UStar = KUZ/LOG(Zzd/z0m)      !Initial setting of u* and calc. of L_MOD (neutral situation)
      IF (ABS(H_init) < 0.001) THEN    ! prevent zero TStar
         h = 0.001
      ELSE
         h = H_init
      END IF
      TStar = (-H/UStar)
      L_MOD = (UStar**2)/(G_T_K*TStar)

      IF (LOG(zzd/z0m) < 0.001000) THEN
         ! PRINT*, 1/(z0m-z0m)
         CALL ErrorHint(17, 'In stability subroutine, (z-zd) < z0.', zzd, z0m, notUsedI)
      ENDIF
      print *,'QH = ', H
      DO i = 1, 330 !Iteration starts
         LOLD = L_MOD
         zL = zzd/L_MOD
         z0L = z0m/L_MOD  !z0m roughness length

         ! IF (z  L>2)THEN
         !    CALL ErrorHint(73,'LUMPS_atmos_functions_stab.f95: stability scale (z/L), UStar',zL,UStar,notUsedI)
         !    RETURN !MO-theory not necessarily valid above zL>2. Still causes problematic UStar values and correct limit might be 0.3.
         !    !Needs more investigations.
         ! END IF

         psim = stab_fn_mom(StabilityMethod, zL, zL)
         psimz0 = stab_fn_mom(StabilityMethod, zL, z0L)

         UStar = KUZ/(LOG(Zzd/z0m) - psim + psimz0) !Friction velocity in non-neutral situation

         TStar = (-H/UStar)
         L_MOD = (UStar**2)/(G_T_K*TStar)

         ! IF(UStar<0.001000)THEN       !If u* too small
         !    UStar=KUZ/(LOG(Zzd/z0m))
         !    PRINT*, 'UStar info',UStar,KUZ,(LOG(Zzd/z0m)),Zzd,z0m
         !    CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] zl,dectime',zl,dectime,notUsedI)
         !    CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] z0l,UStar',z0l,UStar,notUsedI)
         !    CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] psim,psimz0',psim,psimz0,notUsedI)
         !    CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] AVU1,log(zzd/z0m)',AVU1,LOG(zzd/z0m),notUsedI)
         !
         !    ! RETURN
         ! ENDIF

         IF (ABS(LOLD - L_MOD) < 0.01) THEN
            IF (ABS(L_MOD) > 1e6) L_MOD = L_MOD/ABS(L_MOD)*1e6
              EXIT
            CONTINUE
         ENDIF
      ENDDO

      ! limit zL to be within [-5,2]
      IF (zL < -5 .OR. zL > 2) THEN
         zL = MIN(2., MAX(-5., zL))
         ! limit other output variables as well as z/L
         L_MOD = zzd/zL
         z0L = z0m/L_MOD
         psim = stab_fn_mom(StabilityMethod, zL, zL)
         psimz0 = stab_fn_mom(StabilityMethod, zL, z0L)
         ! TS 01 Aug 2018: set a low limit at 0.15 m/s (Schumann 1987, BLM)
         ! to prevent potential issues in other stability-related calcualtions
         UStar = MAX(0.15, KUZ/(LOG(Zzd/z0m) - psim + psimz0))
         TStar = (-H/UStar)
      END IF

      ! TS: limit UStar and TStar to reasonable values
      ! 02 Aug 2018: set a low limit at 0.15 m/s (Schumann 1987, BLM)
      UStar = MAX(0.15, UStar)
      TStar = (-H/UStar)
      IF (UStar < 0.0001) THEN       !If u* still too small after iteration, then force quit simulation and write out error info
         ! UStar=KUZ/(LOG(Zzd/z0m))
         PRINT *, 'UStar', UStar, KUZ, (LOG(Zzd/z0m)), Zzd, z0m
         CALL ErrorHint(30, 'SUBROUTINE STAB_lumps:[ u*< 0.0001] zl,dectime', zl, dectime, notUsedI)
         CALL ErrorHint(30, 'SUBROUTINE STAB_lumps:[ u*< 0.0001] z0l,UStar', z0l, UStar, notUsedI)
         CALL ErrorHint(30, 'SUBROUTINE STAB_lumps:[ u*< 0.0001] psim,psimz0', psim, psimz0, notUsedI)
         CALL ErrorHint(30, 'SUBROUTINE STAB_lumps:[ u*< 0.0001] AVU1,log(zzd/z0m)', AVU1, LOG(zzd/z0m), notUsedI)

         ! RETURN
      ENDIF

      ! if ( L_MOD<-990 ) then
      !   print*, 'L_MOD too low',L_MOD
      !   print*, 1/(L_MOD-L_MOD)
      !
      ! end if

   END SUBROUTINE STAB_lumps

   !==================================================================

   FUNCTION stab_fn_mom(StabilityMethod, ZL, zl_f) RESULT(psym)
      !     StabilityMethod = 1-4 -
      !     PSYM - stability FUNCTION for momentum
      !Modified by LJ Mar 2010
      !Input:Used stability method, stability (z-d)/L, zeta (either (z-d)/L or z0/L)

      ! USE mod_z
      ! USE mod_k

      IMPLICIT NONE
      REAL(KIND(1d0)), PARAMETER :: &
         !  k=0.4,&             !Von Karman's contant
         !  k2=0.16,&           !Power of Van Karman's contant
         neut_limit = 0.001000 !Limit for neutral stability
      !  notUsedI=-55

      REAL(KIND(1d0)):: piover2, psym, zl, zl_f, x, x2
      INTEGER ::StabilityMethod

      PIOVER2 = ACOS(-1.)/2.
      !PRINT*,StabilityMethod,zl,"stab_fn_mom:"
      IF (ABS(zL) < neut_limit) THEN
         psym = 0
      ELSEIF (zL < -neut_limit) THEN    !Unstable

         IF (StabilityMethod == 1) THEN     !    Jensen et al 1984 - Van Ulden & Holtslag (1985) p 1206&
            psym = ((1.-16.*zl_f)**0.25) - 1
         ELSEIF (StabilityMethod == 2) THEN !Dyer (1974)(1-16z/L)**.25' k=0.41  mod. Hogstrom (1988)v15.2
            X = (1.-(15.2*zl_f))**0.25
            X2 = LOG((1 + (X**2.))/2.)
            PSYM = (2.*LOG((1 + X)/2.)) + X2 - (2.*ATAN(X)) + PIOVER2
         ELSEIF (StabilityMethod == 3) THEN ! Kondo (1975) adopted by Campbell & Norman eqn 7.26 p 97
            psym = 0.6*(2)*LOG((1 + (1 - 16*zl_f)**0.5)/2)
         ELSEIF (StabilityMethod == 4) THEN !Businger et al (1971) modifed  Hogstrom (1988)
            x = (1 - 19.3*zl_f)**(-0.25)
            X2 = LOG((1 + (X**2.))/2.)
            PSYM = (2.*LOG((1 + X)/2.)) + X2 - (2.*ATAN(X)) + PIOVER2
         ELSEIF (StabilityMethod == 7) THEN ! Dyer & Bradley (1982) (1-28z/L)**.25' k=0.4
            X = (1 - (28.*zl_f))**0.25  ! NT: changed + to - (bug -> checked reference)
            X2 = LOG((1 + X**2.)/2.)
            PSYM = (2.*LOG((1 + X)/2.)) + X2 - (2.*ATAN(X)) + PIOVER2
         ELSEIF (StabilityMethod == 5) THEN ! Zilitinkevich & Chalikov (1968) modified Hogstrom (1988)
            IF (zl_f >= -0.16) THEN
               x = 1 + 1.38*zl_f
            ELSE
               x = 0.42*(-1)*zl_f**0.333
            ENDIF
            X2 = LOG((1 + (X**2.))/2.)
            PSYM = (2.*LOG((1 + X)/2.)) + X2 - (2.*ATAN(X)) + PIOVER2

         ELSEIF (StabilityMethod == 6) THEN !     Foken and Skeib (1983)
            IF (zl_f >= 0.06) THEN
               x = 1
            ELSE
               x = ((-1)*zl_f/0.06)**0.25
            ENDIF
            X2 = LOG((1 + (X**2.))/2.)
            PSYM = (2.*LOG((1 + X)/2.)) + X2 - (2.*ATAN(X)) + PIOVER2
         ENDIF

      ELSEIF (zL > neut_limit) THEN            !Stable

         IF (StabilityMethod == 1) THEN         !Dyer (1974) k=0.35 x=1+5*zl Mod. Hogstrom (1988)
            psym = (-4.8)*zl_f
         ELSEIF (StabilityMethod == 2) THEN     !Van Ulden & Holtslag (1985) p 1206
            IF (zl_f > 1000.) THEN
               zl_f = 1000.
            END IF
            PSYM = (-17.*(1.-EXP(-0.29*zl_f)))
         ELSEIF (StabilityMethod == 4) THEN ! Businger et al (1971) modifed  Hogstrom (1988)
            ! psym=1+6*zl_f  ! this is NOT the integral form but the stability function, TS 13 Jun 2017
            psym = (-6)*zl_f   ! this is the integral form, TS 13 Jun 2017
         ELSEIF (StabilityMethod == 3) THEN ! Kondo (1975) adopted by Campbell & Norman eqn 7.26 p 97
            psym = (-6)*LOG(1 + zl_f)

         ENDIF
      ENDIF
      RETURN
   END FUNCTION stab_fn_mom

   FUNCTION stab_phi_mom(StabilityMethod, ZL, zl_f) RESULT(phim)
      !     StabilityMethod = 1-4 -
      !     phi - stability FUNCTION for momentum
      !Modified by NT May 2019 !!!!! check if all are correct!
      !Input:Used stability method, stability (z-d)/L, zeta (either (z-d)/L or z0/L)

      IMPLICIT NONE
      REAL(KIND(1d0)), PARAMETER :: &
         !  k=0.4,&             !Von Karman's contant
         !  k2=0.16,&           !Power of Van Karman's contant
         neut_limit = 0.001000 !Limit for neutral stability
      !  notUsedI=-55

      REAL(KIND(1d0)):: phim, zl, zl_f
      INTEGER ::StabilityMethod

      IF (ABS(zL) < neut_limit) THEN
         phim = 0
      ELSEIF (zL < -neut_limit) THEN    !Unstable

         IF (StabilityMethod == 1) THEN     !    Jensen et al 1984 - Van Ulden & Holtslag (1985) p 1206&
            phim = ((1.-16.*zl_f)**(-0.25))
         ELSEIF (StabilityMethod == 2) THEN !Dyer (1974)(1-16z/L)**.25' k=0.41  mod. Hogstrom (1988)v15.2
            phim = (1.-(15.2*zl_f))**(-0.25)
         ELSEIF (StabilityMethod == 3) THEN ! Kondo (1975) adopted by Campbell & Norman eqn 7.26 p 97
            phim = ((1.-16.*zl_f)**(-0.25))
         ELSEIF (StabilityMethod == 4) THEN !Businger et al (1971) modifed  Hogstrom (1988)
            phim = (1. - 19.*zl_f)**(-0.25)
         ELSEIF (StabilityMethod == 7) THEN ! Dyer & Bradley (1982) (1-28z/L)**.25' k=0.4
            phim = (1. - (28.*zl_f))**(-0.25)
         ELSEIF (StabilityMethod == 5) THEN ! Zilitinkevich & Chalikov (1968) modified Hogstrom (1988)
            IF (zl_f >= -0.16) THEN
               phim = 1 + 1.38*zl_f
            ELSE
               phim = 0.42*(-1)*zl_f*(-0.333)
            ENDIF
         ELSEIF (StabilityMethod == 6) THEN !     Foken and Skeib (1983)
            IF (zl_f >= 0.06) THEN
               phim = 1
            ELSE
               phim = ((-1)*zl_f/0.06)**(-0.25)
            ENDIF
         ENDIF

      ELSEIF (zL > neut_limit) THEN            !Stable

         IF (StabilityMethod == 1) THEN         !Dyer (1974) k=0.35 x=1+5*zl Mod. Hogstrom (1988)
            phim = 1.+(4.8)*zl_f
         ELSEIF (StabilityMethod == 2) THEN     !Van Ulden & Holtslag (1985) p 1206 ! NT: have no function for phim 
            phim = 1.+(4.8)*zl_f
         ELSEIF (StabilityMethod == 4) THEN ! Businger et al (1971) modifed  Hogstrom (1988)
            phim=1+6*zl_f 
         ELSEIF (StabilityMethod == 3) THEN ! Kondo (1975) adopted by Campbell & Norman eqn 7.26 p 97  !!NT: checked 
            phim = 1. + 6.* zl_f/(1. + zl_f)  !!NT: checked reference and updated 
         ENDIF
      ENDIF
      RETURN
   END FUNCTION stab_phi_mom  
   !_______________________________________________________________
   !
   ! PSYH - stability function for heat
   FUNCTION stab_fn_heat(StabilityMethod, ZL, zl_f) RESULT(psyh)
      ! USE mod_k
      IMPLICIT NONE
      REAL(KIND(1d0)), PARAMETER :: &
         !  k=0.4,&             !Von Karman's contant
         !  k2=0.16,&           !Power of Van Karman's contant
         neut_limit = 0.001000 !Limit for neutral stability
      !  notUsedI=-55

      REAL(KIND(1d0)):: zl, zl_f, psyh, x
      INTEGER :: StabilityMethod

      IF (ABS(zl) < neut_limit) THEN      !Neutral
         psyh = 0
      ELSEIF (zL < -neut_limit) THEN     ! Unstable
         IF (StabilityMethod == 3) THEN
            !campbell & norman eqn 7.26
            psyh = 0.6*(2)*LOG((1 + (1 - 16*zl_f)**0.5)/2)
         ELSE

            IF (StabilityMethod == 4) THEN ! Businger et al (1971) modifed  Hogstrom (1988)
               x = 0.95*(1.-11.6*zl_f)**(-0.5)
            ELSEIF (StabilityMethod == 7) THEN
               x = (1 - (28.*ZL))**0.25
            ELSEIF (StabilityMethod == 2) THEN ! Dyer 1974 X=(1.-(16.*ZL))**(0.5)modified Hosgstrom
               x = 0.95*(1.-15.2*zl_f)**0.5
            ENDIF
            PSYH = 2*LOG((1 + x**2)/2)
         ENDIF

      ELSE IF (zL > neut_limit) THEN    !Stable
         IF (zL <= 1) THEN ! weak/moderate stable
            IF (StabilityMethod == 4) THEN !Businger et al (1971) modifed  Hogstrom (1988)
               ! psyh=0.95+(7.8*zl_f) ! this is NOT the integral form but the stability function, TS 13 Jun 2017
               psyh = (-7.8)*zl_f   ! this is the integral form, TS 13 Jun 2017
            ELSE !Dyer (1974)  PSYH=(-5)*ZL        modifed  Hogstrom (1988)
               PSYH = (-4.5)*Zl_f
            ENDIF
         ELSE !zL>1, very stable. otherwise psyh would be too large. TS 13 Jun 2017
            ! adopt the form as Brutasert (1982) eqn 4.58. but following the coeffs. of the above eqns
            IF (StabilityMethod == 4) THEN !Businger et al (1971) modifed  Hogstrom (1988)
               psyh = (-7.8)*(1 + LOG(zl_f))
            ELSE !Dyer (1974)  PSYH=(-5)*ZL        modifed  Hogstrom (1988)
               PSYH = (-4.5)*(1 + LOG(zl_f))
            ENDIF
         END IF

      ENDIF

      RETURN
   END FUNCTION stab_fn_heat
   !--------------------------------------------------------------------------------
   ! psys - roughness sublayer correction psi_*
   !
   !     Garratt (1980) QJRMS Appendix 1 p 815/816

   FUNCTION stab_fn_rou(z, zstar) RESULT(psys)
      IMPLICIT NONE
      REAL(KIND(1d0))::alpha, zeta, z, psys, zstar, alpha1
      !     z wind speed height - z_d
      !     zstar height of the roughness sublayer
      !     eqn (a4) using alpha=0.5 alpha1=0.7
      alpha = 0.5
      alpha1 = 0.7
      zeta = z/zstar
      psys = (alpha - 1)*LOG(zeta) - (alpha*alpha1)*(1 - zeta) - (alpha*alpha1**2) &
             *(1 - zeta**2)/6.-(alpha*alpha1**3)*(1 - zeta**3)/24.
      RETURN
   END FUNCTION stab_fn_rou

END MODULE AtmMoistStab_module
