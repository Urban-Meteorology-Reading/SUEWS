SUBROUTINE diagSfc(xSurf,xFlux,us,VegFraction,xDiag,opt)
  ! Ting Sun 20 May 2017: calculate surface-level diagonostics

  USE mod_k
  USE mod_z
  USE sues_data
  USE data_in
  USE moist

  IMPLICIT NONE

  REAL(KIND(1d0)),INTENT(in) :: xSurf,xFlux,us,VegFraction
  REAL(KIND(1d0)),INTENT(out):: xDiag
  INTEGER,INTENT(in)         :: opt ! 0 for momentum, 1 for temperature, 2 for humidity

  REAL(KIND(1d0))            :: &
       psymz2,psymz10,psymz0,psyhz2,psyhz0,& ! stability correction functions
       z0h,& ! Roughness length for heat
       z2zd,z10zd,&
       muu=1.46e-5,& !molecular viscosity
       stab_fn_mom,stab_fn_heat !stability correction functions

  !***************************************************************
  ! log-law based stability corrections:
  ! Roughness length for heat
  IF (RoughLenHeatMethod==1) THEN !Brutasert (1982) z0h=z0/10(see Grimmond & Oke, 1986)
     z0h=z0m/10
  ELSEIF (RoughLenHeatMethod==2) THEN ! Kawai et al. (2007)
     !z0h=z0m*exp(2-(1.2-0.9*veg_fr**0.29)*(us*z0m/muu)**0.25)
     ! Changed by HCW 05 Nov 2015 (veg_fr includes water; VegFraction = veg + bare soil)
     z0h=z0m*EXP(2-(1.2-0.9*VegFraction**0.29)*(us*z0m/muu)**0.25)
  ELSEIF (RoughLenHeatMethod==3) THEN
     z0h=z0m*EXP(-20.) ! Voogt and Grimmond, JAM, 2000
  ELSEIF (RoughLenHeatMethod==4) THEN
     z0h=z0m*EXP(2-1.29*(us*z0m/muu)**0.25) !See !Kanda and Moriwaki (2007),Loridan et al. (2010)
  ENDIF

  ! z0h=z0m/5


  ! zX-z0
  z2zd=2+z0h   ! set lower limit as z0h to prevent arithmetic error
  z10zd=10+z0m ! set lower limit as z0m to prevent arithmetic error

  ! stability correction functions
  ! momentum:
  psymz10=stab_fn_mom(StabilityMethod,z10zd/L_mod,z10zd/L_mod)
  psymz2=stab_fn_mom(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
  psymz0=stab_fn_mom(StabilityMethod,z0m/L_mod,z0m/L_mod)

  ! heat and vapor: assuming both are the same

  psyhz2=stab_fn_heat(StabilityMethod,z2zd/L_mod,z2zd/L_mod)
  psyhz0=stab_fn_heat(StabilityMethod,z0h/L_mod,z0h/L_mod)
  !***************************************************************

  SELECT CASE (opt)
  CASE (0) ! wind (momentum) at 10 m
     xDiag=us/k*(LOG(z10zd/z0m)-psymz10+psymz0)

  CASE (1) ! temperature at 2 m
     xDiag=xSurf-xFlux/(k*us*avdens*avcp)*(LOG(z2zd/z0h)-psyhz2+psyhz0)
     IF ( ABS((LOG(z2zd/z0h)-psyhz2+psyhz0))>10 ) THEN
        PRINT*, '#####################################'
        PRINT*, 'xSurf',xSurf
        PRINT*, 'xFlux',xFlux
        PRINT*, 'k*us*avdens*avcp',k*us*avdens*avcp
        PRINT*, 'k',k
        PRINT*, 'us',us
        PRINT*, 'avdens',avdens
        PRINT*, 'avcp',avcp
        PRINT*, 'xFlux/X',xFlux/(k*us*avdens*avcp)
        PRINT*, 'stab',(LOG(z2zd/z0h)-psyhz2+psyhz0)
        PRINT*, 'LOG(z2zd/z0h)',LOG(z2zd/z0h)
        PRINT*, 'z2zd',z2zd,'L_mod',L_mod,'z0h',z0h
        PRINT*, 'z2zd/L_mod',z2zd/L_mod
        PRINT*, 'psyhz2',psyhz2
        PRINT*, 'psyhz0',psyhz0
        PRINT*, 'psyhz2-psyhz0',psyhz2-psyhz0
        PRINT*, 'xDiag',xDiag
        PRINT*, '*************************************'
     END IF


  CASE (2) ! humidity at 2 m
     xDiag=xSurf-xFlux/(k*us*avdens*tlv)*(LOG(z2zd/z0h)-psyhz2+psyhz0)

  END SELECT

END SUBROUTINE diagSfc
