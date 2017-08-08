!.c!! For Lumps Version 2 - no stability calculations
!==========================================================
!     Last change:
!     TS   13 Jun 2017: corrections to the integral of stability functions
!     MH   12 Apr 2017: Stable limit to exit do-loop
!     LJ   25 Nov 2014: Limits for L
!     LJ   19 Feb 2010
!     SG   27 Mar 2000    4:44 pm
!     ust - friction velocity
!     L - monin obukhov stability length
!       Van Ulden & Holtslag (1985) JCAM: 24: 1196-1207

SUBROUTINE STAB_lumps(&

                                ! input
     StabilityMethod,&
     dectime,& !Decimal time
     zzd,&     !Active measurement height (meas. height-displac. height)
     z0M,&     !Aerodynamic roughness length
     zdm,&     !Displacement height
     avU1,&    !Average wind speed
     Temp_C,&  !Air temperature
                                ! output:
         L,& !Obukhov length
     Tstar,& !T*
     USTAR,& !Friction velocity
     h,    & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
     psim  & !Stability function of momentum
     )
  ! USE mod_k
  ! USE mod_z
  ! USE mod_grav
  ! USE data_in
  ! USE time
  ! USE defaultnotUsed
  ! USE WhereWhen
  IMPLICIT NONE
  INTEGER,INTENT(in):: StabilityMethod

  REAL(KIND(1d0)),INTENT(in)::&
       dectime,& !Decimal time
       zzd,&     !Active measurement height (meas. height-displac. height)
       z0M,&     !Aerodynamic roughness length
       zdm,&     !Displacement height
       avU1,&    !Average wind speed
       Temp_C    !Air temperature


  REAL(KIND(1d0)),INTENT(out)::&
       L,    &!Obukhov length
       Tstar,&!T*
       USTAR,&!Friction velocity
       h,&    !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
       psim   !Stability function of momentum

  REAL(KIND(1d0))::stab_fn_mom,&
       G_T_k,&
       KUZ,&
       LOLD,&
       zL,&
       z0l,&
       stab_f,&
       psimz0
  REAL(KIND(1d0)),PARAMETER :: &
       k=0.4,&             !Von Karman's contant
       k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000,& !Limit for neutral stability
       grav=9.80665,&  !g - gravity - physics today august 1987
       notUsedI=-55

  INTEGER :: i

  LOGICAL :: debug=.FALSE.

  IF(debug) WRITE(*,*)StabilityMethod,z0M,avU1,h,USTAR,L
  G_T_k=(Grav/(Temp_C+273.16))*k !gravity constant/(Temperature*Von Karman Constant)
  KUZ=k*AvU1                     !Von Karman constant*mean wind speed
  IF(zzd<0) CALL ErrorHint(32,'Windspeed Ht too low relative to zdm [Stability calc]- values [z-zdm, zdm]',Zzd,zdm,notUsedI)

  USTAR=KUZ/LOG(Zzd/Z0M)      !Initial setting of u* and calc. of L (neutral situation)
  IF ( ABS(H)<0.001 ) THEN    ! prevent zero Tstar
     H=0.001
  END IF
  Tstar=(-H/ustar)
  L=(USTAR**2)/(G_T_K*Tstar)


  IF(LOG(zzd/z0M)<0.001000) CALL ErrorHint(17,'In stability subroutine, (z-zd) < z0.',zzd,z0m,notUsedI)
  DO i=1,330 !Iteration starts
     LOLD=L
     zL=zzd/L
     z0L=z0M/L  !z0M roughness length

     IF (zL>2)THEN
        CALL ErrorHint(73,'LUMPS_atmos_functions_stab.f95: stability parameter, ustar',zL,USTAR,notUsedI)
        RETURN !MO-theory not necessarily valid above zL>2. Still causes problematic ustar values and correct limit might be 0.3.
        !Needs more investigations.
     END IF

     psim=stab_fn_mom(StabilityMethod,zL,zL)
     psimz0=stab_fn_mom(StabilityMethod,zL,z0L)


     USTAR=KUZ/(LOG(Zzd/Z0M)-PSIM+psimz0) !Friction velocity in non-neutral situation

     IF(ustar<0.001000)THEN       !If u* too small
        USTAR=KUZ/(LOG(Zzd/Z0M))
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] zl,dectime',zl,dectime,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] z0l,ustar',z0l,ustar,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] psim,psimz0',psim,psimz0,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] AVU1,log(zzd/z0m)',AVU1,LOG(zzd/z0m),notUsedI)

        RETURN
     ENDIF

     tstar=(-H/ustar)
     L=(Ustar**2)/(G_T_K*Tstar)

     IF(ABS(LOLD-L)<0.01)THEN
        IF (ABS(L)>1e6) L = L/ABS(L)*1e6
        RETURN
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE stab_lumps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION stab_fn_mom(StabilityMethod,ZL,zl_f) RESULT(psym)
  !     StabilityMethod = 1-4 -
  !     PSYM - stability FUNCTION for momentum
  !Modified by LJ Mar 2010
  !Input:Used stability method, stability (z-d)/L, zeta (either (z-d)/L or z0/L)

  ! USE mod_z
  ! USE mod_k

  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER :: &
       k=0.4,&             !Von Karman's contant
       k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000,& !Limit for neutral stability
       notUsedI=-55

  REAL (KIND(1d0)):: piover2,psym,zl,zl_f,x,x2
  INTEGER ::StabilityMethod

  PIOVER2=ACOS(-1.)/2.
  !PRINT*,StabilityMethod,zl,"stab_fn_mom:"
  IF(ABS(zL)<neut_limit) THEN
     psym=0
  ELSEIF(zL<-neut_limit) THEN    !Unstable

     IF(StabilityMethod==1)THEN     !    Jensen et al 1984 - Van Ulden & Holtslag (1985) p 1206&
        psym=((1.-16.*zl_f)**0.25)-1
     ELSEIF(StabilityMethod==2) THEN !Dyer (1974)(1-16z/L)**.25' k=0.41  mod. Hogstrom (1988)v15.2
        X=(1.-(15.2*zl_f))**0.25
        X2=LOG((1+(X**2.))/2.)
        PSYM=(2.*LOG((1+X)/2.))+X2-(2.*ATAN(X))+PIOVER2
     ELSEIF(StabilityMethod==3)THEN !     campbell & norman eqn 7.26
        psym=0.6*(2)*LOG((1+(1-16*zl_f)**0.5)/2)
     ELSEIF(StabilityMethod==4) THEN !Businger et al (1971) modifed  Hogstrom (1988)
        x=(1-19.3*zl_f)**(-0.25)
        X2=LOG((1+(X**2.))/2.)
        PSYM=(2.*LOG((1+X)/2.))+X2-(2.*ATAN(X))+PIOVER2
     ELSEIF(StabilityMethod==7) THEN ! Dyer & Bradley (1982) (1-28z/L)**.25' k=0.4
        X=(1+(28.*zl_f))**0.25
        X2=LOG((1+X**2.)/2.)
        PSYM=(2.*LOG((1+X)/2.))+X2-(2.*ATAN(X))+PIOVER2
     ELSEIF(StabilityMethod==5)THEN ! Zilitinkevich & Chalikov (1968) modified Hogstrom (1988)
        IF(zl_f>=-0.16)THEN
           x=1+1.38*zl_f
        ELSE
           x=0.42*(-1)*zl_f**0.333
        ENDIF
        X2=LOG((1+(X**2.))/2.)
        PSYM=(2.*LOG((1+X)/2.))+X2-(2.*ATAN(X))+PIOVER2

     ELSEIF(StabilityMethod==6)THEN !     Foken and Skeib (1983)
        IF(zl_f>=0.06)THEN
           x=1
        ELSE
           x=((-1)*zl_f/0.06)**0.25
        ENDIF
        X2=LOG((1+(X**2.))/2.)
        PSYM=(2.*LOG((1+X)/2.))+X2-(2.*ATAN(X))+PIOVER2
     ENDIF

  ELSEIF(zL>neut_limit) THEN            !Stable

     IF(StabilityMethod==1)THEN         !Dyer (1974) k=0.35 x=1+5*zl Mod. Hogstrom (1988)
        psym=(-4.8)*zl_f
     ELSEIF(StabilityMethod==2)THEN     !Van Ulden & Holtslag (1985) p 1206
        IF ( zl_f >1000. ) THEN
           zl_f=1000.
        END IF
        PSYM=(-17.*(1.-EXP(-0.29*zl_f)))
     ELSEIF(StabilityMethod==4)THEN ! Businger et al (1971) modifed  Hogstrom (1988)
        ! psym=1+6*zl_f  ! this is NOT the integral form but the stability function, TS 13 Jun 2017
        psym=(-6)*zl_f   ! this is the integral form, TS 13 Jun 2017
     ELSEIF(StabilityMethod==3)THEN ! campbell & norman eqn 7.27 p 97
        psym=(-6)*LOG(1+zl_f)

     ENDIF
  ENDIF
  RETURN
END FUNCTION stab_fn_mom

!_______________________________________________________________
!
! PSYH - stability function for heat
FUNCTION stab_fn_heat(StabilityMethod,ZL,zl_f) RESULT (psyh)
  ! USE mod_k
  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER :: &
       k=0.4,&             !Von Karman's contant
       k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000,& !Limit for neutral stability
       notUsedI=-55

  REAL (KIND(1d0)):: zl,zl_f,psyh,x
  INTEGER :: StabilityMethod

  IF(ABS(zl)<neut_limit)THEN      !Neutral
     psyh=0
  ELSEIF(zL<-neut_limit) THEN     ! Unstable
     IF(StabilityMethod==3)THEN
        !campbell & norman eqn 7.26
        psyh=0.6*(2)*LOG((1+(1-16*zl_f)**0.5)/2)
     ELSE

        IF(StabilityMethod==4)THEN ! Businger et al (1971) modifed  Hogstrom (1988)
           x=0.95*(1.-11.6*zl_f)**(-0.5)
        ELSEIF(StabilityMethod==7) THEN
           x=(1-(28.*ZL))**0.25
        ELSEIF(StabilityMethod==2)THEN ! Dyer 1974 X=(1.-(16.*ZL))**(0.5)modified Hosgstrom
           x=0.95*(1.-15.2*zl_f)**0.5
        ENDIF
        PSYH=2*LOG((1+x**2)/2)
     ENDIF

  ELSE IF (zL>neut_limit) THEN    !Stable
     IF ( zL<=1 ) THEN ! weak/moderate stable
        IF(StabilityMethod==4)THEN !Businger et al (1971) modifed  Hogstrom (1988)
           ! psyh=0.95+(7.8*zl_f) ! this is NOT the integral form but the stability function, TS 13 Jun 2017
           psyh=(-7.8)*zl_f   ! this is the integral form, TS 13 Jun 2017
        ELSE !Dyer (1974)  PSYH=(-5)*ZL	modifed  Hogstrom (1988)
           PSYH=(-4.5)*Zl_f
        ENDIF
     ELSE !zL>1, very stable. otherwise psyh would be too large. TS 13 Jun 2017
        ! adopt the form as Brutasert (1982) eqn 4.58. but following the coeffs. of the above eqns
        IF(StabilityMethod==4)THEN !Businger et al (1971) modifed  Hogstrom (1988)
           psyh=(-7.8)*(1+LOG(zl_f))
        ELSE !Dyer (1974)  PSYH=(-5)*ZL	modifed  Hogstrom (1988)
           PSYH=(-4.5)*(1+LOG(zl_f))
        ENDIF
     END IF

  ENDIF

  RETURN
END FUNCTION stab_fn_heat
!--------------------------------------------------------------------------------
! psys - roughness sublayer correction psi_*
!
!     Garratt (1980) QJRMS Appendix 1 p 815/816

FUNCTION stab_fn_rou(z,zstar) RESULT (psys)
  IMPLICIT NONE
  REAL(KIND(1d0))::alpha,zeta,z,psys,zstar,alpha1
  !     z wind speed height - z_d
  !     zstar height of the roughness sublayer
  !     eqn (a4) using alpha=0.5 alpha1=0.7
  alpha=0.5
  alpha1=0.7
  zeta=z/zstar
  psys=(alpha-1)* LOG(zeta)-(alpha*alpha1)*(1-zeta)-(alpha*alpha1**2) &
       *(1-zeta**2)/6.-(alpha*alpha1**3)*(1-zeta**3)/24.
  RETURN
END FUNCTION stab_fn_rou
