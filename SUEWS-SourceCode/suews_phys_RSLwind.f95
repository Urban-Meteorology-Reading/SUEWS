SUBROUTINE WindProfile( &
   UStar, L_MOD, sfr, Zh, planF, StabilityMethod, &  ! input
   zarray, dataoutLineURSL) ! output
   !-----------------------------------------------------
   ! calculates windprofiles using MOST with a RSL-correction
   ! based on Harman & Finnigan 2007
   !
   ! last modified by:
   ! NT 03/2019
   !
   !-----------------------------------------------------
   USE AtmMoistStab_module, ONLY: STAB_lumps, stab_fn_mom, stab_phi_mom

   IMPLICIT NONE
   INTEGER, PARAMETER:: nsurf = 7 ! number of surface types
   INTEGER, PARAMETER:: BldgSurf = 2
   INTEGER, PARAMETER:: ConifSurf = 3
   INTEGER, PARAMETER:: DecidSurf = 4

   REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in) ::sfr! surface fractions
   REAL(KIND(1d0)), INTENT(in):: UStar ! Friction velocity [m s-1]
   REAL(KIND(1d0)), INTENT(in):: L_MOD  ! Obukhov length [m]
   REAL(KIND(1d0)), INTENT(in):: Zh    ! Mean building height [m]
   REAL(KIND(1d0)), INTENT(in):: planF ! Frontal area index [-]
   INTEGER, INTENT(in)::StabilityMethod

   REAL(KIND(1d0)), PARAMETER:: cd_tree = 1.2, & ! drag coefficient tree canopy !!!!needs adjusting!!!
                                a_tree = 0.05, & ! the foliage area per unit volume !!!!needs adjusting!!!
                                kappa = 0.40, &! von karman constant
                                beta_N = 0.40, &  ! H&F beta coefficient in neutral conditions from Theeuwes et al., 2019 BLM
                                pi = 4.*ATAN(1.0), &
                                a1 = 4., a2 = -0.1, a3 = 1.5, a4 = -1. ! constraints to determine beta
   INTEGER, PARAMETER :: nz = 30   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy

   REAL(KIND(1d0)), INTENT(out), DIMENSION(nz):: zarray ! Height array
   REAL(KIND(1d0)), INTENT(out), DIMENSION(nz):: dataoutLineURSL ! Wind speed array

   REAL(KIND(1d0)), DIMENSION(nz)::dif, psihat_z

   REAL(KIND(1d0)):: zd, & ! displacement height
                     Lc_build, Lc_tree, Lc, & ! canopy drag length scale
                     dz, & ! height steps
                     phim, psimz, psimZh, psimz0, phi_hatmZh, phimzp, phimz, &  ! stability function for momentum
                     betaHF, betaNL, beta, betaN2, &  ! beta coefficient from Harman 2012
                     elm, & ! mixing length
                     xx1, xx1_2, err, z01, dphi, &  ! dummy variables for stability functions
                     z0, &  ! roughness length from H&F
                     cm, c2 ! H&F'07 'constants'
   INTEGER :: I, z, it, idx_can

   ! Start setting up the parameters
   ! calculate Lc for tree grid fraction using eq 1 H&F'07 and rest of grid using C&B'04
   Lc_build = (1.-sfr(BldgSurf))/planF*Zh  ! Coceal and Belcher 2004 assuming Cd = 2
   Lc_tree = 1./(cd_tree*a_tree)
   Lc = (1.-(sfr(BldgSurf) + sfr(ConifSurf) + sfr(ConifSurf)))/planF*Zh
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

   phim = stab_phi_mom(StabilityMethod, Lc/L_MOD, Lc/L_MOD)

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
   ! start calculations for above roof height
   ! start with stability at canopy top for z0 and phihat
   psimZh = stab_fn_mom(StabilityMethod, (Zh - zd)/L_MOD, (Zh - zd)/L_MOD)

   ! calculate phihatM according to H&F '07
   xx1 = stab_phi_mom(StabilityMethod, (Zh - zd)/L_MOD, (Zh - zd)/L_MOD)
   xx1_2 = stab_phi_mom(StabilityMethod, (Zh - zd + 1.)/L_MOD, (Zh - zd + 1.)/L_MOD)

   phi_hatmZh = kappa/(2.*beta*xx1)
   dphi = xx1_2 - xx1
   IF (phi_hatmZh > 1.) THEN
      c2 = 0.5 ! more stable but less correct
   ELSE
      c2 = (kappa*(3.-(2.*beta**2.*Lc/xx1*dphi)))/(2.*beta*xx1 - kappa)  ! if very unstable this might cause some high values of psihat_z
   ENDIF

   cm = (1.-phi_hatmZh)*EXP(c2/2.)

   psihat_z = 0.*zarray
   DO z = idx_can, nz - 1
      phimz = stab_phi_mom(StabilityMethod, (zarray(z) - Zh + zd)/L_MOD, (zarray(z) - Zh + zd)/L_MOD)
      phimzp = stab_phi_mom(StabilityMethod, (zarray(z + 1) - Zh + zd)/L_MOD, (zarray(z + 1) - Zh + zd)/L_MOD)

      psihat_z(z) = psihat_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - Zh + zd)/elm)) &  !Taylor's approximation for integral
                    /(zarray(z + 1) - Zh + zd)
      psihat_z(z) = psihat_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - Zh + zd)/elm)) &
                    /(zarray(z) - Zh + zd)
   ENDDO

   ! calculate z0 iteratively
   z0 = 0.5  !first guess
   err = 10.
   it = 1
   DO WHILE ((err > 0.001) .AND. (it < 10))
      psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)
      z01 = z0
      z0 = (Zh - zd)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0)*EXP(psihat_z(idx_can))
      err = ABS(z01 - z0)
      it = it + 1
   ENDDO

   psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)

   ! calculate above canopy wind speed
   DO z = idx_can, nz
      psimz = stab_fn_mom(StabilityMethod, (zarray(z) - zd)/L_MOD, (zarray(z) - zd)/L_MOD)
      dataoutLineURSL(z) = UStar/kappa*(LOG((zarray(z) - zd)/z0) - psimz + psimz0 - psihat_z(z) + psihat_z(idx_can))
   ENDDO
   ! calculate in canopy wind speed
   DO z = 1, idx_can
      dataoutLineURSL(z) = dataoutLineURSL(idx_can)*EXP(beta*(zarray(z) - Zh)/elm)
   ENDDO

END SUBROUTINE WindProfile
