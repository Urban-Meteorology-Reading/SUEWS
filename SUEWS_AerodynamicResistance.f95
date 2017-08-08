SUBROUTINE AerodynamicResistance(&

     ! input:
     ZZD,&
     z0m,&
     AVU1,&
     L_mod,&
     Ustar,&
     VegFraction,&
     AerodynamicResistanceMethod,&
     StabilityMethod,&
     RoughLenHeatMethod,&
     ! output:
     RA)

  ! Returns Aerodynamic resistance (RA) to the main program SUEWS_Calculations
  ! All ra equations reported in Thom & Oliver (1977)
  ! Modified by TS - interface modified
  ! Modified by LJ
  !   -Removal of tabs and cleaning the code
  ! Modified by HCW 03 Dec 2015 - changed lower limit on ra from 2 s m-1 to 10 s m-1 (to avoid unrealistically high evaporation rates)
  ! Modified by LJ in 12 April to not to be used with global variables
  ! To Do:
  !       - Check whether the thresholds 2-200 s m-1 are suitable over a range of z0!! HCW 04 Mar 2015
  ! OUTPUT: RA - Aerodynamic resistance [s m^-1]
  ! INPUT:  AerodynamicResistanceMethod = Method to calculate RA
  !         StabilityMethod = defines the method to calculate atmospheric stability
  !         RoughLenHeatMethod = Method to calculate heat roughness length
  !         *Measurement height minus* Displacement height (m) (was incorrectly labelled, corrected HCW 25 May 2016
  !         z0m = Aerodynamic roughness length (m)
  !         k2 = Power of Van Karman's constant (= 0.16 = 0.4^2)
  !         AVU1 = Mean wind speed
  !         L_mod = Obukhov length (m)
  !         Ustar = Friction velocity (m s-1)
  !         VegFraction = Fraction of vegetation
  !               (changed from veg_fr which also includes water surface by HCW 05 Nov 2015)


  ! USE DefaultNotUsed

  IMPLICIT NONE

  REAL(KIND(1d0)),INTENT(in)::&
       ZZD,&      !Active measurement height (meas. height-displac. height)
       z0m,&      !Aerodynamic roughness length
       AVU1,&     !Average wind speed
       L_mod,&    !Monin-Obukhov length (either measured or modelled)
       Ustar,&    !Friction velocity
       VegFraction!Fraction of vegetation
  INTEGER,INTENT(in)::&
       AerodynamicResistanceMethod,&
       StabilityMethod,&
       RoughLenHeatMethod

  REAL(KIND(1d0)),INTENT(out)::RA !Aerodynamic resistance [s m^-1]

  INTEGER, PARAMETER :: notUsedI=-55
  REAL(KIND(1d0)), PARAMETER :: &
       notUsed=-55.5,&
       k2=0.16,& !Power of Van Karman's constant (= 0.16 = 0.4^2)
       muu=1.46e-5 !molecular viscosity
  REAL(KIND(1d0))::stab_fn_heat,stab_fn_mom,&
       psym,&
       psyh,z0V


  !1)Monteith (1965)-neutral stability
  IF(AerodynamicResistanceMethod==1) THEN
     RA=(LOG(ZZD/z0m)**2)/(k2*AVU1)

     !2) Non-neutral stability
     !    PSIM - stability function for momentum
     !     PSYH - stability function for heat
     !    assuming stability functions the same for heat and water
  ELSEIF(AerodynamicResistanceMethod==2) THEN  !Dyer (1974)

     psym=stab_fn_mom(StabilityMethod,ZZD/L_mod,zzd/L_mod)
     psyh=stab_fn_heat(StabilityMethod,ZZD/L_mod,zzd/L_mod)

     !Z0V roughness length for vapour
     IF (RoughLenHeatMethod==1) THEN !Brutasert (1982) Z0v=z0/10(see Grimmond & Oke, 1986)
        z0V=Z0m/10
     ELSEIF (RoughLenHeatMethod==2) THEN ! Kawai et al. (2007)
       	!z0V=Z0m*exp(2-(1.2-0.9*veg_fr**0.29)*(Ustar*Z0m/muu)**0.25)
        ! Changed by HCW 05 Nov 2015 (veg_fr includes water; VegFraction = veg + bare soil)
        z0V=Z0m*EXP(2-(1.2-0.9*VegFraction**0.29)*(Ustar*Z0m/muu)**0.25)
     ELSEIF (RoughLenHeatMethod==3) THEN
        z0V=Z0m*EXP(-20.) ! Voogt and Grimmond, JAM, 2000
     ELSEIF (RoughLenHeatMethod==4) THEN
        z0V=Z0m*EXP(2-1.29*(Ustar*Z0m/muu)**0.25) !See !Kanda and Moriwaki (2007),Loridan et al. (2010)
     ENDIF

     IF(Zzd/L_mod==0.OR.Ustar==0) THEN
        RA=(LOG(ZZD/z0m)*LOG(ZZD/z0V))/(k2*AVU1) !Use neutral equation
     ELSE
        RA=((LOG(ZZD/z0m)-PSYM)*(LOG(ZZD/z0V)-PSYH))/(K2*AVU1)
     ENDIF

     !3) Thom and Oliver (1977)
  ELSEIF(AerodynamicResistanceMethod==3) THEN
     RA=(4.72*LOG(ZZD/z0m)**2)/(1 + 0.54*AVU1)
  ENDIF

  !If ra outside permitted range, adjust extreme values !!Check whether these thresholds are suitable over a range of z0
  IF(RA>200) THEN   !was 175
     CALL errorHint(7,'In AerodynamicResistance.f95, calculated ra > 200 s m-1; ra set to 200 s m-1',RA,notUsed,notUsedI)
     RA=200
  ELSEIF(RA<10) THEN   !found  By Shiho - fix Dec 2012  !Threshold changed from 2 to 10 s m-1 (HCW 03 Dec 2015)
     CALL errorHint(7,'In AerodynamicResistance.f95, calculated ra < 10 s m-1; ra set to 10 s m-1',RA,notUsed,notUsedI)
     RA=10
     ! RA=(log(ZZD/z0m))**2/(k2*AVU1)
     IF(avu1<0) WRITE(*,*) avu1,ra
  ENDIF

  RETURN
END SUBROUTINE AerodynamicResistance
