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
   USE AtmMoistStab_module, ONLY: STAB_lumps, stab_fn_mom

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
                                a1 = 4., a2 = -0.1, a3 = 1.5, a4 = 1. ! constraints to determine beta
   INTEGER, PARAMETER :: nz = 30   ! number of levels 10 levels in canopy plus 20 (3 x Zh) above the canopy

   REAL(KIND(1d0)), INTENT(out), DIMENSION(nz):: zarray ! Height array
   REAL(KIND(1d0)), INTENT(out), DIMENSION(nz):: dataoutLineURSL ! Wind speed array

   REAL(KIND(1d0)), DIMENSION(nz)::psihat_z, psihat_z0

   REAL(KIND(1d0)):: zd, & ! displacement height
                     Lc_build, Lc_tree, Lc, & ! canopy drag length scale
                     dz, & ! height steps
                     phim, psimz, psimZh, psimz0, phi_hatmZh, phimzp, phimz, &  ! stability function for momentum
                     betaHF, betaNL, beta, &  ! beta coefficient from Harman 2012
                     elm, & ! mixing length
                     x1, xx1, xx1_2, err, z01, dphi, &  ! dummy variables for stability functions
                     z0, &  ! roughness length from H&F
                     cm, c2 ! H&F'07 'constants'
   INTEGER :: I, z, it

   ! Start setting up the parameters
   ! calculate Lc for tree grid fraction using eq 1 H&F'07 and rest of grid using C&B'04
   Lc_build = (1.-sfr(BldgSurf))/planF*Zh  ! Coceal and Belcher 2004 assuming Cd = 2
   Lc_tree = 1./(cd_tree*a_tree)
   Lc = (1.-(sfr(BldgSurf) + sfr(ConifSurf) + sfr(ConifSurf)))/planF*Zh
   print *, 'Lc', Lc
   print *, 'Z_H',Zh
   print *, 'FAI',planF
   print *, 'PAI',(sfr(BldgSurf) + sfr(ConifSurf) + sfr(ConifSurf))
   dz = Zh/10
   zarray = (/(I, I=1, nz)/)*dz

   !Method of determining beta from Harman 2012, BLM
   IF (Lc/L_MOD < 0.) THEN
      phim = (1.-16*Lc/L_MOD)**0.25
   ELSE
      phim = 1.+5.*Lc/L_MOD
   ENDIF

   betaHF = beta_N/phim
   betaNL = (kappa/2.)/phim

   IF (Lc/L_MOD > a2) THEN
      beta = betaHF
   ELSEIF (Lc/L_MOD < a4) THEN
      beta = beta_N
   ELSE
      beta = betaNL + ((betaHF - betaNL)/(1.+a1*abs(Lc/L_MOD - a2)**a3))
   ENDIF
   zd = Zh - (beta**2.)*Lc
   elm = 2.*beta**3*Lc

   ! start calculations for above roof height
   ! start with stability at canopy top for z0 and phihat
   IF ((Zh - zd)/L_MOD < 0.) THEN
      psimZh = stab_fn_mom(StabilityMethod, (Zh - zd)/L_MOD, (Zh - zd)/L_MOD)

      ! calculate phihatM according to H&F '07
      xx1 = (1.-16.*(Zh - zd)/L_MOD)**(-0.25)
      xx1_2 = (1.-16.*(Zh - zd + 1.)/L_MOD)**(-0.25)
      phi_hatmZh = kappa/(2.*beta*xx1)
      dphi = xx1_2 - xx1
      c2 = (kappa*(3.-(2.*beta**2.*Lc/xx1*dphi)))/(2.*beta*xx1 - kappa)
      cm = (1.-phi_hatmZh)*EXP(c2/2.)

      psihat_z = 0.*zarray
      DO z = 9, nz - 1
         phimz = (1.-16.*(zarray(z) - zd)/L_MOD)**(-0.25)
         print *, 'phimz here',phimz, L_MOD
         phimzp = (1.-16.*(zarray(z + 1) - zd)/L_MOD)**(-0.25)
         print *, 'phimzp here',phimzp
         psihat_z(z) = psihat_z(z + 1) + dz/2.*phimzp*(cm*EXP(-1.*c2*beta*(zarray(z + 1) - zd)/elm)) &  !Taylor's approximation for integral
                       /(zarray(z + 1) - zd)
         psihat_z(z) = psihat_z(z) + dz/2.*phimz*(cm*EXP(-1.*c2*beta*(zarray(z) - zd)/elm)) &
                       /(zarray(z) - zd)
      ENDDO

      print *, zarray(10), psihat_z(10)
      ! calculate z0 iteratively
      z0 = 0.5  !first guess
      err = 10.
      DO it = 1, 10
         psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)
         z01 = z0
         z0 = (Zh - zd)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0)*EXP(psihat_z(10))
         err = ABS(z01 - z0)
         IF (err < 0.001) EXIT
      ENDDO
      psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)

   ELSE   ! stable
      psimZh = stab_fn_mom(StabilityMethod, (Zh - zd)/L_MOD, (Zh - zd)/L_MOD)

      ! calculate phihatM according to H&F '07
      xx1 = (1.-5.*(Zh - zd)/L_MOD)
      xx1_2 = (1.-5.*(Zh - zd + 1.)/L_MOD)
      phi_hatmZh = kappa/(2.*beta*xx1)
      dphi = xx1_2 - xx1
      c2 = (kappa*(3.-(2.*beta**2.*Lc/xx1*dphi)))/(2.*beta*xx1 - kappa)
      cm = (1.-phi_hatmZh)*EXP(c2/2.)

      psihat_z = 0.*zarray
      DO z = 9, nz - 1
         phimz = (1.+5.*(zarray(z) - zd)/L_MOD)
         phimzp = (1.+5.*(zarray(z + 1) - zd)/L_MOD)
         psihat_z(z) = psihat_z(z + 1) + dz/2.*phimzp*(cm* &    !Taylor's approximation for integral
                                                       EXP(-1.*c2*beta*(zarray(z + 1) - zd)/elm))/(zarray(z + 1) - zd)
         psihat_z(z) = psihat_z(z) + dz/2.*phimz*(cm* &
                                                  EXP(-1.*c2*beta*(zarray(z) - zd)/elm))/(zarray(z) - zd)
      ENDDO

      ! calculate z0 iteratively
      z0 = 0.5  !first guess
      err = 10.
      DO it = 1, 10
         psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)
         z01 = z0
         z0 = (Zh - zd)*EXP(-1.*kappa/beta)*EXP(-1.*psimZh + psimz0) !* np.exp(PSIhat_Zh) add RSL correction term
         err = ABS(z01 - z0)
         IF (err < 0.001) EXIT
      ENDDO
      psimz0 = stab_fn_mom(StabilityMethod, z0/L_MOD, z0/L_MOD)

   ENDIF
!    print *, z0, it
   ! calculate above canopy wind speed
   DO z = 10, nz
      psimz = stab_fn_mom(StabilityMethod, (zarray(z) - zd)/L_MOD, (zarray(z) - zd)/L_MOD)
      dataoutLineURSL(z) = UStar/kappa*(LOG((zarray(z) - zd)/z0) - psimz + psimz0 - psihat_z(z) + psihat_z(10))
   ENDDO
   PRINT *, 'Uarray', dataoutLineURSL
   ! calculate in canopy wind speed
   DO z = 1, 10
      dataoutLineURSL(z) = dataoutLineURSL(10)*EXP(beta*(zarray(z) - Zh)/elm)
   ENDDO

!    DO z = 1, nz
!       print *, dataoutLineURSL(z)
!    ENDDO
!    DO z = 1, nz
!       print *, zarray(z)
!    ENDDO
END SUBROUTINE WindProfile
