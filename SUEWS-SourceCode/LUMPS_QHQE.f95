SUBROUTINE LUMPS_cal_QHQE(&
     veg_type,& !input
     snowUse,id,qn1,qf,qs,Qm,Temp_C,Veg_Fr,avcp,Press_hPa,lv_J_kg,&
     tstep_real,DRAINRT,nsh_real,&
     Precip,RainMaxRes,RAINCOVER,sfr,LAI,LAImax,LAImin,&
     H_mod,& !output
     E_mod,psyc_hPa,s_hPa,sIce_hpa,TempVeg,VegPhenLumps)
  !Calculates QH and QE for LUMPS. See Loridan et al. (2011)
  ! ref: Grimmond and Oke (2002) JAM and references within that
  !      Offerle (2003) -- add water bucket
  ! ref: Loridan et al. (2011) JAMC dynamic water & vegetation
  ! Last modified:
  ! LJ 27 Jan 2016  - Removal of tabs, cleaning the code
  ! HCW 04 Mar 2015 - Modified to account for model timestep (rather than hourly resolution)
  ! LJ Feb 2014     - The bug related to VegMax has been fixed (cannot divide by zero)
  ! LJ/SG May 2012  - Changed phenology to be consistent with SUEWS LAI. No longer Loridan et al. (2011)
  ! LJ June 2012    - Modified to work with snow (Qm added in the equations!)
  ! SG Feb 2012     - added some comments
  ! --------------------------------------------------------------
  USE meteo,ONLY:psyc_const,slope_svp,slopeice_svp

  IMPLICIT NONE
  INTEGER,PARAMETER::ndays=366
  INTEGER,PARAMETER::NSurf=7
  INTEGER,PARAMETER::NVegSurf=3
  INTEGER,PARAMETER::ivConif=1
  INTEGER,PARAMETER::ivGrass=3

  INTEGER,INTENT(in) :: veg_type  !Defines how vegetation is calculated for LUMPS
  INTEGER,INTENT(in) :: snowUse ! option of snow module
  INTEGER,INTENT(in) :: id ! day of year

  REAL(KIND(1d0)),INTENT(in) :: qn1! net all-wave radiation
  REAL(KIND(1d0)),INTENT(in) :: qf! anthropogenic heat flux
  REAL(KIND(1d0)),INTENT(in) :: qs! storage heat flux
  REAL(KIND(1d0)),INTENT(in) :: Qm!Snow melt associated heat flux
  REAL(KIND(1d0)),INTENT(in) :: Temp_C!air temperature in degC
  REAL(KIND(1d0)),INTENT(in) :: Veg_Fr!Vegetation fraction from land area
  REAL(KIND(1d0)),INTENT(in) :: avcp!Specific heat capacity
  REAL(KIND(1d0)),INTENT(in) :: Press_hPa!Station air pressure in hPa
  REAL(KIND(1d0)),INTENT(in) :: lv_J_kg!Latent heat of vaporization in [J kg-1]
  REAL(KIND(1d0)),INTENT(in) :: tstep_real ! time step in REAL
  REAL(KIND(1d0)),INTENT(in) :: DRAINRT!Drainage rate of the water bucket [mm hr-1]
  REAL(KIND(1d0)),INTENT(in) :: nsh_real! real cast of Number of timesteps per hour
  REAL(KIND(1d0)),INTENT(in) :: Precip!Precipitation per timestep [mm]
  REAL(KIND(1d0)),INTENT(in) :: RainMaxRes!Maximum water bucket reservoir [mm]
  REAL(KIND(1d0)),INTENT(in) :: RAINCOVER! LUMPS Limit when surface totally wet [mm]

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in) :: sfr! veg surface fractions [-]
  REAL(KIND(1D0)),DIMENSION(-4:NDAYS,NVEGSURF),INTENT(in) :: LAI! LAI(id-1,iv), LAI at the beginning of today
  REAL(KIND(1d0)),DIMENSION(3),INTENT(in) :: LAImax!Max LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(3),INTENT(in) :: LAImin    !Min LAI [m2 m-2]

  REAL(KIND(1d0)),INTENT(out) ::H_mod
  REAL(KIND(1d0)),INTENT(out) ::E_mod !turbulent fluxes: QH, QE
  REAL(KIND(1d0)),INTENT(out) ::psyc_hPa !Psychometric constant in hPa
  REAL(KIND(1d0)),INTENT(out) ::s_hPa!Vapour pressure versus temperature slope in hPa
  REAL(KIND(1d0)),INTENT(out) ::sIce_hpa!Vapour pressure versus temperature slope in hPa above ice/snow
  REAL(KIND(1d0)),INTENT(out) ::TempVeg !TEMPORARY VEGETATIVE SURFACE FRACTION ADJUSTED BY RAINFALL
  REAL(KIND(1d0)),INTENT(out) ::VegPhenLumps
  ! REAL(KIND(1d0)),INTENT(inout) ::RainBucket !RAINFALL RESERVOIR [mm]
  ! INTEGER::iv
  REAL(KIND(1d0)),DIMENSION(3) :: sfrVeg! veg surface fractions [-]                             !,start
  REAL(KIND(1d0)),DIMENSION(3) :: LAIDay! LAI(id-1,iv), LAI at the beginning of today
  REAL(KIND(1d0))::VegPhen,VegMax,VegMin,&   !Vegetation phenology for LUMPS
       psyc_s,&       !Psychometric constant
       alpha_sl,alpha_in,&    	  !Parameters used in LUMPS QH and QE calculations
       beta,&                      !Beta parameter used in LUMPS QH and QE calculations [W m-2]
       alpha_qhqe,RAINRES,RainBucket,tlv

  tlv=lv_J_kg/tstep_real !Latent heat of vapourisation per timestep
  ! initialize VegPhenLumps to output
  VegPhenLumps=0
  ! initialize rain-related variables
  RainBucket=0.

  sfrVeg=sfr(ivConif+2:ivGrass+2)

  LAIDay= LAI(id-1,veg_type)

  ! Calculate slope of the saturation vapour pressure vs air temp.
  s_hPa=slope_svp(Temp_C)
  psyc_hPa=psyc_const(avcp,Press_hPa,lv_J_kg)
  psyc_s=psyc_hPa/s_hPa

  !Calculate also sublimation ones if snow calculations are made.
  !Used also for LUMPS
  IF (snowUse==1) THEN
     IF (Temp_C<=0) THEN
        sIce_hpa=slopeIce_svp(Temp_C)
     ELSE
        sIce_hpa=slope_svp(Temp_C)
     ENDIF
     psyc_s=psyc_hPa/sIce_hPa   !Psychometric constant divided by the slope
  ENDIF

  ! replaced by sinusoidal vegetation formulation
  !alpha=gis(idgis,itgis,1)*alpha_sl+alpha_in

  !THE FOLLOWING ADJUSTS THE ALPHA and BETA PARAMETERs FOR RAINFALL.
  !ASSUMES THE SURFACE IS VEGETATION COVERED WITH RAIN > RAINCOVER mm/DAY
  !OTHERWISE INCREASES VEGETATION LINEAR WITH AMOUNT OF RAIN.

  ! !IF (E_MOD>0.) RainBucket=RainBucket-E_MOD*1.44E-3 !1.44E-3 MM/(W/M^2)/HR (i.e. 3600/(lv_J_kg))
  ! IF (E_MOD>0.) RainBucket=RainBucket-E_MOD/tlv   !Adjusted for per model timestep instead of per hour HCW 04 Mar 2015
  ! IF (Temp_C>0.) RainBucket=RainBucket - DRAINRT/nsh_real  !DRAINRT is specified in mm h-1
  ! IF (RainBucket<0.) RainBucket=0.
  ! IF (Precip>0) RainBucket=MIN(RainMaxRes,RainBucket+Precip)
  !
  ! RAINRES = RainBucket
  ! IF (RAINRES>RAINCOVER) RAINRES=RAINCOVER

  !--------Calculate vegetation phenology for LUMPS------------------------
  ! VegPhen=0
  ! VegMax=0
  ! VegMin=0
  VegPhen=DOT_PRODUCT(sfrVeg,LAIDay)
  VegMax=DOT_PRODUCT(sfrVeg,LAImax)
  VegMin=DOT_PRODUCT(sfrVeg,LAImin)

  ! DO iv=ivConif,ivGrass   !Normalized LAI for vegetation
  !    VegPhen = sfr(iv+2)*LAI(id-1,iv) + VegPhen
  !    VegMax  = sfr(iv+2)*LAImax(iv) + VegMax
  !    VegMin  = sfr(iv+2)*LAImax(iv) + VegMin
  ! ENDDO

  IF(VegMax<=0.01000) THEN   !If max vegetation is very small, TempVeg = 0;
     TempVeg=0
  ELSE
     VegPhenLumps=(VegPhen)/(VegMax)
     TempVeg=Veg_Fr*VegPhenLumps   !Now this is veg_fraction in general
  ENDIF

  IF (TempVeg>0.9000) THEN   !If vegetation fraction is larger than 0.7 (0.9?)
     beta = (20-3)*TempVeg+3
     alpha_qhqe=TempVeg*0.8+0.2
  ELSE
     beta=3
     IF(veg_type==1) THEN   !Area vegetated, including bare soil and water
        alpha_sl=0.686
        alpha_in=0.189
     ELSEIF(veg_type==2) THEN   !Area irrigated vegetation
        alpha_sl=0.610
        alpha_in=0.222
     ENDIF
     alpha_qhqe=TempVeg*alpha_sl+alpha_in
  ENDIF

  ! Calculate the actual heat fluxes
  H_mod= ((1-alpha_qhqe)+psyc_s)/(1+psyc_s)*(qn1+qf-qs-Qm)-beta   !Eq 3, Grimmond & Oke (2002)
  E_mod= (alpha_qhqe/(1+psyc_s)*(qn1+qf-qs-Qm))+beta              !Eq 4, Grimmond & Oke (2002)

  ! adjust RAINRES after E_mod calculation is done: ! moved here from above. TS, 13 Jan 2018
  !IF (E_MOD>0.) RainBucket=RainBucket-E_MOD*1.44E-3 !1.44E-3 MM/(W/M^2)/HR (i.e. 3600/(lv_J_kg))
  IF (E_MOD>0.) RainBucket=RainBucket-E_MOD/tlv   !Adjusted for per model timestep instead of per hour HCW 04 Mar 2015
  IF (Temp_C>0.) RainBucket=RainBucket - DRAINRT/nsh_real  !DRAINRT is specified in mm h-1
  IF (RainBucket<0.) RainBucket=0.
  IF (Precip>0) RainBucket=MIN(RainMaxRes,RainBucket+Precip)

  RAINRES = RainBucket
  IF (RAINRES>RAINCOVER) RAINRES=RAINCOVER


  RETURN

END SUBROUTINE LUMPS_cal_QHQE
