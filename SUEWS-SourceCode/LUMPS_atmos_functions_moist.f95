MODULE AtmMoist_module
  IMPLICIT NONE

CONTAINS
  !.c!! For Lumps Version 2 - no stability calculations
  ! Latent heat of sublimation when air temperature below zero added. LJ Nov 2012
  ! explict interface added to all subroutines, TS 08 Aug 2017
  SUBROUTINE LUMPS_cal_AtmMoist(&
       Temp_C,Press_hPa,avRh,dectime,&! input:
       lv_J_kg,lvS_J_kg,&! output:
       es_hPa,Ea_hPa,VPd_hpa,VPD_Pa,dq,dens_dry,avcp,air_dens)

    IMPLICIT NONE
    REAL(KIND(1d0))::vap_dens

    REAL(KIND(1d0)),INTENT(in)::&
         Temp_C,&
         Press_hPa,&
         avRh,dectime
    REAL(KIND(1d0)),INTENT(out)::&
         lv_J_kg,&!Latent heat of vaporization in [J kg-1]
         lvS_J_kg,&!Latent heat of sublimation in J/kg
         es_hPa,&!Saturation vapour pressure over water in hPa
         Ea_hPa,&!Vapour pressure of water in hPa
         VPd_hpa,& !vapour pressure deficit in hPa
         VPD_Pa,& !vapour pressure deficit in Pa
         dq,&!Specific humidity deficit in g/kg
         dens_dry,& !Vap density or absolute humidity	 (kg/m3)
         avcp,&!specific heat capacity in J kg-1 K-1
         air_dens!Air density in kg/m3


    REAL (KIND(1d0)),PARAMETER:: &
                                !  comp          = 0.9995, &
                                !  epsil         = 0.62197,&           !ratio molecular weight of water vapor/dry air (kg/mol/kg/mol)
                                !  epsil_gkg     = 621.97, &           !ratio molecular weight of water vapor/dry air in g/kg
                                !  dry_gas       = 8.31451,&           !Dry gas constant (J/k/mol)
                                !  gas_ct_wat    = 461.05,&            !Gas constant for water (J/kg/K)
                                !  molar         = 0.028965,&          !Dry air molar fraction in kg/mol
                                !  molar_wat_vap = 0.0180153,&         !Molar fraction of water vapor in kg/mol
         gas_ct_dry    = 8.31451/0.028965,&  !j/kg/k=dry_gas/molar
         gas_ct_wv     = 8.31451/0.0180153 !j/kg/kdry_gas/molar_wat_vap
    !  waterDens     = 999.8395            !Density of water in 0 cel deg
    INTEGER::from=1

    !Saturation vapour pressure over water in hPa
    es_hPa = sat_vap_press(Temp_C,Press_hPa,from,dectime) ! dectime is more or less unnecessary here

    !Vapour pressure of water in hPa
    Ea_hPa=avRh/100*es_hPa

    ! if(debug.and.dectime>55.13.and.dectime<55.2)write(35,*)'%',Temp_C

    VPd_hpa=es_hPa-ea_hpa           !vapour pressure deficit in hPa
    VPD_Pa=(es_hPa*100)-(Ea_hPa*100)!vapour pressure deficit in Pa

    dq=(spec_hum_def(vpd_hPa,Press_hPa)) !Specific humidity deficit in g/kg

    !Vap density or absolute humidity	 (kg/m3)
    vap_dens=(Ea_hPa*100/((Temp_C+273.16)*gas_ct_wv))

    !density Dry Air Beer(1990)	kg/m3
    dens_dry=((Press_hPa-Ea_hPa)*100)/(gas_ct_dry*(273.16+Temp_C))

    !Air density in kg/m3
    air_dens=(Press_hPa*100)/(gas_ct_dry*(Temp_C+273.16))

    !Calculate specific heat capacity in J kg-1 K-1
    avcp=spec_heat_beer(Temp_C,avRh,vap_dens,dens_dry)

    !Latent heat of vaporization in [J kg-1]
    lv_J_kg=lat_vap(Temp_C,Ea_hPa,Press_hPa,avcp,dectime)

    !Latent heat of sublimation in J/kg
    IF(Temp_C<0.000) THEN
       lvS_J_kg=lat_vapSublim(Temp_C,Ea_hPa,Press_hPa,avcp)
    ELSE
       lvS_J_kg = 2834000
    ENDIF
    !if(debug)write(*,*)lv_J_kg,Temp_C,'lv2'
    IF(press_hPa<900) THEN
       CALL ErrorHint(46, 'Function LUMPS_cal_AtmMoist',press_hPa,-55.55, -55)
    ENDIF
    RETURN
  END SUBROUTINE LUMPS_cal_AtmMoist

  !=====================================================================
  ! sg sept 99 f90
  ! Uses eqns from Buck (1981) JAM 20, 1527-1532
  ! units T (deg C) e (mb) P (mb)
  ! f corrects for the fact that we are not dealing with pure water
  ! LJ Feb 2010
  !Changed to use the updated version (Buck research manual, 1996) from Buck (1981)
  !For water different equations in cold and warm temperatures

  FUNCTION sat_vap_press(Temp_c,PRESS_hPa,from,dectime) RESULT(es_hPa)
    ! USE time
    ! USE defaultnotUsed
    IMPLICIT NONE

    REAL(KIND(1d0))::temp_C,press_hpa,dectime!,pw
    REAL(KIND(1d0))::e_mb,f,press_kpa,es_hPA
    INTEGER:: from,iv
    INTEGER,PARAMETER::notUsedI=-55

    !If air temperature between -0.001 -
    IF(ABS(temp_C)<0.001000)THEN
       IF(from==1) THEN  ! not from determining Tw
          iv=INT(press_Hpa)
          CALL errorHint(29,'Function sat_vap_press: temp_C, dectime,press_Hpa = ',temp_C, dectime,iv)

       ENDIF
       temp_C=0.001000
    ENDIF

    Press_kPa=Press_hPa/10

    IF(Temp_C<50.AND.Temp_C>-40)THEN
       !e_mb=6.1121*EXP(((18.729-Temp_C/227.3)*Temp_C)/(Temp_C+257.87)) !Old one
       !f=1.00072+Press_hPa*(3.2E-6+5.9D-10*Temp_C**2)

       IF (Temp_C>=0.001000) THEN
          e_mb=6.1121*EXP(((18.678-Temp_C/234.5)*Temp_C)/(Temp_C+257.14))
          f=1.00072+Press_kPa*(3.2E-6+5.9E-10*Temp_C**2)
          es_hPa=e_mb*f

       ELSEIF (Temp_C<=-0.001000) THEN
          e_mb=6.1115*EXP(((23.036-Temp_C/333.7)*Temp_C)/(Temp_C+279.82))
          f=1.00022+Press_kPa*(3.83E-6+6.4E-10*Temp_C**2)
          es_hPa=e_mb*f
       ENDIF

    ELSE
       CALL ErrorHint(28,'FUNCTION sat_vap_press: [Temperature is out of range], Temp_C,dectime',Temp_C,dectime,notUsedI)

    ENDIF

    RETURN
  END FUNCTION sat_vap_press


  FUNCTION sat_vap_pressIce(Temp_c,PRESS_hPa,from,dectime) RESULT(es_hPa)
    ! USE time
    ! USE defaultnotUsed
    IMPLICIT NONE

    REAL(KIND(1d0))::e_mb,f,temp_C,press_hpa,press_kpa,es_hPA,dectime!,pw
    INTEGER:: from,iv
    INTEGER,PARAMETER::notUsedI=-55

    !If air temperature between -0.001 -
    IF(ABS(temp_C)<0.001000)THEN
       IF(from==1) THEN  ! not from determining Tw
          iv=INT(press_Hpa)
          CALL errorHint(29,'Function sat_vap_press: temp_C, dectime,press_Hpa = ',temp_C, dectime,iv)

       ENDIF
       temp_C=0.001000
    ENDIF

    Press_kPa=Press_hPa/10

    IF(Temp_C<50.AND.Temp_C>-40)THEN
       e_mb=6.1115*EXP(((23.036-Temp_C/333.7)*Temp_C)/(Temp_C+279.82))
       f=1.00022+Press_kPa*(3.83E-6+6.4E-10*Temp_C**2) !In hPa
       es_hPa=e_mb*f

    ELSE
       CALL ErrorHint(28,'FUNCTION sat_vap_press: [Temperature is out of range], Temp_C,dectime',Temp_C,dectime,notUsedI)

    ENDIF

    RETURN
  END FUNCTION sat_vap_pressIce

  !==========================================================
  !Output: specific humidity deficit in g/kg
  !Input: Dry air density and air pressure in hPa
  FUNCTION spec_hum_def(vpd_hPa,press_hPa) RESULT(dq)
    ! USE gas
    IMPLICIT NONE
    REAL(KIND(1d0))           :: press_hPa,vpd_hPa,dq
    REAL(KIND(1d0)),PARAMETER :: epsil_gkg = 621.97 !ratio molecular weight of water vapor/dry air in g/kg
    dq=epsil_gkg*vpd_hPa/press_hPa ! Phd Thesis II.13 p 196
  END FUNCTION spec_hum_def

  ! ==============================================================================
  FUNCTION spec_heat_beer(Temp_C,rh,rho_v,rho_d) RESULT (cp)
    ! Input: Air temperature, relative humidity, water vapour and dry air densities
    ! Output: heat capacity in units J kg-1 K-1
    ! Reference: Tom Beer, CSIRO, 1990. Applied Environmetrics Meteorological Tables.
    ! Can be found from SG:s office from Atmmos Moist map
    !-------------------------------------------------------------------------------

    ! USE defaultnotUsed
    IMPLICIT NONE

    REAL(KIND(1d0))::cp,cpd,cpm,rho_v,rho_d,rh,temp_C

    !Garratt equation a20 (1992)
    CPd = 1005.0+((Temp_C+23.16)**2)/3364.0 !Changed from 23.15 to 23.16

    !Beer (1990) for water vapor
    cpm = 1859 + 0.13*rH+ (19.3+0.569*rH)*(Temp_C/100.) + &
         (10.+0.5*rH)*(Temp_C/100.)**2

    IF(ABS(rho_d)<0.000100.OR.ABS(rho_v)<0.000100.OR.ABS(rho_d+rho_v)<0.000100)THEN
       CALL ErrorHint(42,'spec-heat_beer',rho_v,rho_d,INT(Temp_C))
    ENDIF

    cp=cpd*(rho_d/(rho_d+rho_v))+cpm*(rho_v/(rho_d+rho_v))

    !   print*,"cp: ",cp,cpd,rho_d,rho_v,cpm,rh
  END FUNCTION spec_heat_beer

  !==========================================================
  !Latent_heat.f sg nov 96
  !sg sep 99 converted f90 FUNCTION
  !Added calcualation of latent heat of sublimation, LJ June 2012

  FUNCTION Lat_vap(Temp_C,Ea_hPa,Press_hPa,cp,dectime) RESULT (lv_J_kg)
    !Input: Air temperature, Water vapour pressure, Air pressure, heat capacity
    !Output: latent heat of vaporization

    ! USE time
    ! USE SnowMod
    ! USE defaultnotUsed

    IMPLICIT NONE
    REAL(KIND(1d0))::cp,lv_J_kg,ea_fix,tw,&
         incr,es_tw,psyc,ea_est,press_hPa,ea_HPa, temp_C,dectime!,Temp_K
    ! REAL(KIND(1d0))::sat_vap_press,psyc_const ! functions

    LOGICAL:: switch1=.FALSE.,switch2=.FALSE.!,debug=.true.
    INTEGER:: ii,from=2
    REAL(KIND(1d0)),PARAMETER::notUsed=-55.55

    ea_fix=ea_hPa
    !if(debug) write(*,*)Temp_C, 'LV'
    !Temp_K=temp_C+273.16

    !lv=1.91846E6*(Temp_K/(Temp_K-33.91))**2

    lv_J_kg=(2500.25-2.365*temp_C)*1000  !First guess for lv in units J/kg


    tw=Temp_C/2.  !First estimate for wet bulb temperature
    incr=3.
    DO ii=1,100
       IF(Press_hPa<900) THEN
          CALL ErrorHint(45,'function Lat_vap',Press_hPA,notUsed,ii)
       ENDIF

       ! if(debug.and.dectime>55.13.and.dectime<55.2)write(35,*)'% 1',Tw

       es_tw=sat_vap_press(Tw,Press_hPa,from,dectime)  !Calculate saturation vapour pressure in hPa

       !if(debug.and.dectime>55.13.and.dectime<55.2)write(35,*)'% 2',Tw

       IF(Press_hPa<900) THEN
          CALL ErrorHint(45,'function Lat_vap - 2',Press_hPA,notUsed,ii)
       ENDIF

       psyc=psyc_const(cp,Press_hPa,lv_J_kg) !in units hPa/K

       IF(Press_hPa<900) THEN
          CALL ErrorHint(45,'function Lat _vap -31',Press_hPA,notUsed,ii)
       ENDIF

       ea_est=es_tw-psyc*(temp_C-tw)

       lv_J_kg=(2500.25-2.365*tw)*1e3

       IF(switch1.AND.switch2)THEN
          incr=incr/10.
          switch1=.FALSE.
          switch2=.FALSE.
       ENDIF
       IF(ABS(ea_est-ea_fix)<0.001000)THEN
          RETURN
       ELSEIF(ea_est > ea_fix)THEN
          tw=tw-incr
          switch1=.TRUE.
       ELSEIF(ea_est< ea_fix)THEN
          tw=tw+incr
          switch2=.TRUE.
       ENDIF
    ENDDO

    RETURN
  END FUNCTION Lat_vap


  FUNCTION Lat_vapSublim(Temp_C,Ea_hPa,Press_hPa,cp) RESULT (lvS_J_kg)
    !Input: Air temperature, Water vapour pressure, Air pressure, heat capacity
    !Output: latent heat of sublimation in units J/kg

    ! USE time

    IMPLICIT NONE

    REAL(KIND(1d0))::lvS_J_kg,temp_C,tw,incr,Ea_hPa,Press_hPa,cp
   ! REAL(KIND(1d0))::ea_fix,es_tw,psyc,ea_est,Temp_K
   ! REAL(KIND(1d0))::sat_vap_pressIce,psyc_const ! functions
   ! LOGICAL:: switch1=.FALSE.,switch2=.FALSE.!,debug=.true.
   ! INTEGER:: ii,from=2

    !Latent heat for sublimation
    !From Rogers&Yau (A short course in cloud physics), Wikipedia

   ! ea_fix=ea_hPa

    lvS_J_kg=(2834.1-0.29*temp_C)*1e3 !First guess for Ls in J/kg

    tw=Temp_C/2.  !First estimate for wet bulb temperature
    incr=3.
    Press_hPa=Press_hPa
    Ea_hPa=Ea_hPa
    cp=cp

    !DO ii=1,100

   !    es_tw=sat_vap_pressIce(Tw,Press_hPa,from)  !Calculate saturation vapour pressure in hPa

     ! psyc=psyc_const(cp,Press_hPa,lv_J_kg)

  !   ea_est=es_tw-psyc*(temp_C-tw)
  !  lvS_J_kg=(2834.1-0.29*tw)*1e3

  !   IF(switch1.AND.switch2)THEN
  !      incr=incr/10.
  !     switch1=.FALSE.
  !     switch2=.FALSE.
  !    ENDIF

  !   IF(ABS(ea_est-ea_fix)<0.001)THEN
  !     RETURN
  !   ELSEIF(ea_est > ea_fix)THEN
  !      tw=tw-incr
  !      switch1=.TRUE.
  !   ELSEIF(ea_est< ea_fix)THEN
  !      tw=tw+incr
  !      switch2=.TRUE.
  !    ENDIF
  !   ENDDO

   ! RETURN
  END FUNCTION Lat_vapSublim



  !=====================================================================
  !psyc_const.f   sg   nov 96
  !sg june 99 f90
  !calculate psyc - psychrometic constant Fritschen and Gay (1979)

  FUNCTION psyc_const(cp,Press_hPa,lv_J_kg) RESULT(psyc_hPa) !In units hPa/K
    USE gas

    IMPLICIT NONE
    REAL (KIND(1d0))::cp,lv_J_kg,press_hPa,psyc_hpa

    ! cp for moist air (shuttleworth p 4.13)
    IF(cp*press_hPa<900.OR.lv_J_kg<10000)THEN
       CALL errorHint(19,'in psychrometric constant calculation:  cp [J kg-1 K-1], p [hPa], Lv [J kg-1]',cp,Press_hPa,INT(lv_J_kg))
    ENDIF

    psyc_hPa=(cp*press_hPa)/(epsil*lv_J_kg)
    !    if(debug)write(*,*)psyc_hpa, 'g',cp,press_HPa,lv
    ! LV MJKg-1
    !www.cimis.water.ca.gov/infoEotPmEquation.jsp
    !psyc_hPa=(0.00163*(Press_hPa/10)/LV)*10
    ! write(*,*)psyc_hpa
    !psyc=psyc*100.! convert into Pa
  END FUNCTION psyc_const

  !==========================================================

  FUNCTION dewpoint(ea_hPa) RESULT(Temp_C_dew)
    ! ea = vapor pressure (hPa)
    ! td = dewpoint (oC)
    !calculates dewpoint in degC from
    ! http://www.atd.ucar.edu/weather_fl/dewpoint.html
    !     dewpoint = (237.3 * ln(e_vp/6.1078)) / (17.27 - (ln(e_vp/6.1078)))

    REAL(KIND(1d0))::ea_hPa,temp_C_dew
    Temp_C_dew = (237.3 * LOG(ea_hPa/6.1078)) / (17.27 - (LOG(ea_hPa/6.1078)))
  END FUNCTION dewpoint
  !===============================================================================
  FUNCTION slope_svp(temp_C) RESULT(s_hPa)
    !COEFFICENTS FOR CALCULATING desat/dT
    !Slope of the saturation vapor pressure vst air temperature curve

    IMPLICIT  NONE

    REAL (KIND(1d0)):: b1,b2,b3,b4,b5,b6,b7,S_hPa,temp_C
    B1=4.438099984D-1
    B2=2.857002636D-2
    B3=7.938054040D-4
    B4=1.215215065D-5
    B5=1.036561403D-7
    B6=3.532421810D-10
    B7=-7.090244804D-13

    !     s - slope of saturation vapour pressure curve - Lowe (1977) -T (K)
    !     mb /K
    S_hPa=B1+temp_C*(B2+temp_C*(B3+temp_C*(B4+temp_C*(B5+temp_C*(B6+B7*temp_C)))))
    ! write(*,*)'s',s_hpa,temp_C
    !s_Pa=s_Pa*100  ! Pa/K
    !www.cimis.water.ca.gov/infoEotPmEquation.jsp
    ! s_hPa=(((4099 *(es_hPa/10))/(Temp_C+273.3)**2))*10
    ! if(debug)write(*,*)s_hpa
    RETURN
  END FUNCTION slope_svp

  !===============================================================================
  FUNCTION slopeIce_svp(temp_C) RESULT(s_hPa)
    !COEFFICENTS FOR CALCULATING desat/dT
    !Slope of the saturation vapor pressure vst air temperature curve

    IMPLICIT  NONE

    REAL (KIND(1d0)):: b1,b2,b3,b4,b5,b6,b7,S_hPa,temp_C

    B1=5.030305237D-1
    B2=3.773255020D-2
    B3=1.267995369D-3
    B4=2.477563108D-5
    B5=3.005693132D-7
    B6=2.158542548D-9
    B7=7.131097725D-12

    ! s - slope of saturation vapour pressure curve - Lowe (1977) -T (K)
    ! mb /K
    S_hPa=B1+temp_C*(B2+temp_C*(B3+temp_C*(B4+temp_C*(B5+temp_C*(B6+B7*temp_C)))))

    RETURN
  END FUNCTION slopeIce_svp

  !------------------------------------------------------------------------
  FUNCTION qsatf(T,PMB) RESULT(qsat)
    !       MRR, 1987
    ! AT TEMPERATURE T (DEG C) AND PRESSURE PMB (MB), GET SATURATION SPECIFIC
    !       HUMIDITY (KG/KG) FROM TETEN FORMULA

    REAL (KIND(1D0))::T,es,qsat,PMB

    REAL (KIND(1D0)),PARAMETER::&

                                !Teten coefficients
         A=6.106,&
         B=17.27,&
         C=237.3,&
         molar=0.028965,& !Dry air molar fraction in kg/mol
         molar_wat_vap=0.0180153 !Molar fraction of water vapor in kg/mol


    IF(t.GT.55)THEN
       CALL ErrorHint(34,'Function qsatf',T,0.00D0,-55)
    ENDIF

    ES = A*dEXP(B*T/(C+T))
    qsat = (molar_wat_vap/molar)*ES/PMB!(rmh2o/rmair)*ES/PMB
  END FUNCTION qsatf

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

SUBROUTINE STAB_lumps(&

                                ! input
     StabilityMethod,&
     dectime,& !Decimal time
     zzd,&     !Active measurement height (meas. height-displac. height)
     z0M,&     !Aerodynamic roughness length
     zdm,&     !Displacement height
     avU1,&    !Average wind speed
     Temp_C,&  !Air temperature
     h_init,    & !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity
                                ! output:
     L_MOD,& !Obukhov length
     Tstar,& !T*
     UStar,& !Friction velocity
     psim)!Stability function of momentum

  IMPLICIT NONE
  INTEGER,INTENT(in):: StabilityMethod


  REAL(KIND(1d0)),INTENT(in)::dectime !Decimal time
  REAL(KIND(1d0)),INTENT(in)::zzd     !Active measurement height (meas. height-displac. height)
  REAL(KIND(1d0)),INTENT(in)::z0M     !Aerodynamic roughness length
  REAL(KIND(1d0)),INTENT(in)::zdm     !Displacement height
  REAL(KIND(1d0)),INTENT(in)::avU1    !Average wind speed
  REAL(KIND(1d0)),INTENT(in)::Temp_C    !Air temperature
  REAL(KIND(1d0)),INTENT(in)::h_init    !Kinematic sensible heat flux [K m s-1] used to calculate friction velocity


  REAL(KIND(1d0)),INTENT(out)::L_MOD!Obukhov length
  REAL(KIND(1d0)),INTENT(out)::Tstar!T*
  REAL(KIND(1d0)),INTENT(out)::UStar!Friction velocity
  REAL(KIND(1d0)),INTENT(out)::psim   !Stability function of momentum

  REAL(KIND(1d0))::G_T_k,&
       KUZ,&
       LOLD,&
       zL,&
       z0l,&
       psimz0,&
       h
  REAL(KIND(1d0)),PARAMETER :: &
       k=0.4,&             !Von Karman's contant
       grav=9.80665,&  !g - gravity - physics today august 1987
       notUsedI=-55

  INTEGER :: i

  LOGICAL :: debug=.FALSE.

  IF(debug) WRITE(*,*)StabilityMethod,z0M,avU1,h_init,UStar,L_MOD
  G_T_k=(Grav/(Temp_C+273.16))*k !gravity constant/(Temperature*Von Karman Constant)
  KUZ=k*AvU1                     !Von Karman constant*mean wind speed
  IF(zzd<0) CALL ErrorHint(32,'Windspeed Ht too low relative to zdm [Stability calc]- values [z-zdm, zdm]',Zzd,zdm,notUsedI)

  UStar=KUZ/LOG(Zzd/Z0M)      !Initial setting of u* and calc. of L_MOD (neutral situation)
  IF ( ABS(h_init)<0.001 ) THEN    ! prevent zero Tstar
     h=0.001
  ELSE
     h=h_init
  END IF
  Tstar=(-H/UStar)
  L_MOD=(UStar**2)/(G_T_K*Tstar)


  IF(LOG(zzd/z0M)<0.001000) CALL ErrorHint(17,'In stability subroutine, (z-zd) < z0.',zzd,z0m,notUsedI)
  DO i=1,330 !Iteration starts
     LOLD=L_MOD
     zL=zzd/L_MOD
     z0L=z0M/L_MOD  !z0M roughness length

     IF (zL>2)THEN
        CALL ErrorHint(73,'LUMPS_atmos_functions_stab.f95: stability parameter, UStar',zL,UStar,notUsedI)
        RETURN !MO-theory not necessarily valid above zL>2. Still causes problematic UStar values and correct limit might be 0.3.
        !Needs more investigations.
     END IF

     psim=stab_fn_mom(StabilityMethod,zL,zL)
     psimz0=stab_fn_mom(StabilityMethod,zL,z0L)


     UStar=KUZ/(LOG(Zzd/Z0M)-PSIM+psimz0) !Friction velocity in non-neutral situation

     IF(UStar<0.001000)THEN       !If u* too small
        UStar=KUZ/(LOG(Zzd/Z0M))
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] zl,dectime',zl,dectime,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] z0l,UStar',z0l,UStar,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] psim,psimz0',psim,psimz0,notUsedI)
        CALL ErrorHint(30,'SUBROUTINE STAB_lumps:[ u*< 0.001] AVU1,log(zzd/z0m)',AVU1,LOG(zzd/z0m),notUsedI)

        RETURN
     ENDIF

     tstar=(-H/UStar)
     L_MOD=(UStar**2)/(G_T_K*Tstar)

     IF(ABS(LOLD-L_MOD)<0.01)THEN
        IF (ABS(L_MOD)>1e6) L_MOD = L_MOD/ABS(L_MOD)*1e6
        RETURN
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE STAB_lumps

!==================================================================

FUNCTION stab_fn_mom(StabilityMethod,ZL,zl_f) RESULT(psym)
  !     StabilityMethod = 1-4 -
  !     PSYM - stability FUNCTION for momentum
  !Modified by LJ Mar 2010
  !Input:Used stability method, stability (z-d)/L, zeta (either (z-d)/L or z0/L)

  ! USE mod_z
  ! USE mod_k

  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER :: &
                                !  k=0.4,&             !Von Karman's contant
                                !  k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000 !Limit for neutral stability
  !  notUsedI=-55

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
                                !  k=0.4,&             !Von Karman's contant
                                !  k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000 !Limit for neutral stability
  !  notUsedI=-55

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


END MODULE AtmMoist_module
