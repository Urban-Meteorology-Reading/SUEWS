SUBROUTINE LUMPS_QHQE(&
                                !input:
     qn1,qf,qs,Qm,Temp_C,avcp,Press_hPa,lv_J_kg,tlv,&
     DRAINRT,Precip,&
     RainBucket,RainMaxRes,RAINCOVER,&
     Veg_Fr,sfrVeg,laiDay,LAImax,LAImin,&
     nsh_real,veg_type,snowUse,&
                                !output:
     H_mod,E_mod,&
     psyc_hPa,s_hPa,sIce_hpa,TempVeg)
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

  ! USE allocateArray
  ! USE data_in
  ! USE defaultNotUsed
  ! USE gis_data
  ! USE moist
  ! USE snowMod
  ! USE sues_data
  ! use time
  ! USE VegPhenogy

  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(in)   ::&
       qn1,&        ! net all-wave radiation
       qf,&         ! anthropogenic heat flux
       qs,&         ! storage heat flux
       Qm,&         !Snow melt associated heat flux
       Temp_C,&     !air temperature in degC
       Veg_Fr,&     !Vegetation fraction from land area
       avcp,&       !Specific heat capacity
       Press_hPa,&  !Station air pressure in hPa
       lv_J_kg,&    !Latent heat of vaporization in [J kg-1]
       tlv,&        !Latent heat of vaporization per timestep
       DRAINRT,&    !Drainage rate of the water bucket [mm hr-1]
       nsh_real,&   ! real cast of Number of timesteps per hour
       Precip,&     !Precipitation per timestep [mm]
       RainMaxRes,& !Maximum water bucket reservoir [mm]
       RAINCOVER,&  ! LUMPS Limit when surface totally wet [mm]
       sfrVeg(3),&  ! veg surface fractions [-]
       laiDay(3),&  ! lai(id-1,iv), LAI at the beginning of today
       LAImax(3),&  !Max LAI [m2 m-2]
       LAImin(3)    !Min LAI [m2 m-2]
  REAL(KIND(1d0)),INTENT(inout) ::RainBucket !RAINFALL RESERVOIR [mm]
  INTEGER,INTENT(in) ::&
       veg_type,&  !Defines how vegetation is calculated for LUMPS
       snowUse ! option of snow module
  REAL(KIND(1d0)),INTENT(out)   ::&
       H_mod,E_mod,& !turbulent fluxes: QH, QE
       psyc_hPa,& !Psychometric constant in hPa
       s_hPa,& !Vapour pressure versus temperature slope in hPa
       sIce_hpa,&!Vapour pressure versus temperature slope in hPa above ice/snow
       TempVeg !TEMPORARY VEGETATIVE SURFACE FRACTION ADJUSTED BY RAINFALL


  INTEGER::iv                                 !,start
  REAL(KIND(1d0))::VegPhen,VegMax,VegMin,&   !Vegetation phenology for LUMPS
       slope_svp, slopeIce_svp,& !Slope of the saturation vapour pressure curve above watre and ice
       psyc_s,psyc_const,&       !Psychometric constant
       alpha_sl,alpha_in,&    	  !Parameters used in LUMPS QH and QE calculations
       beta,&                      !Beta parameter used in LUMPS QH and QE calculations [W m-2]
       alpha_qhqe,VegPhenLumps,RAINRES

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

  !IF (E_MOD>0.) RainBucket=RainBucket-E_MOD*1.44E-3 !1.44E-3 MM/(W/M^2)/HR (i.e. 3600/(lv_J_kg))
  IF (E_MOD>0.) RainBucket=RainBucket-E_MOD/tlv   !Adjusted for per model timestep instead of per hour HCW 04 Mar 2015
  IF (Temp_C>0.) RainBucket=RainBucket - DRAINRT/nsh_real  !DRAINRT is specified in mm h-1
  IF (RainBucket<0.) RainBucket=0.
  IF (Precip>0) RainBucket=MIN(RainMaxRes,RainBucket+Precip)

  RAINRES = RainBucket
  IF (RAINRES>RAINCOVER) RAINRES=RAINCOVER

  !--------Calculate vegetation phenology for LUMPS------------------------
  ! VegPhen=0
  ! VegMax=0
  ! VegMin=0
  VegPhen=DOT_PRODUCT(sfrVeg,laiDay)
  VegMax=DOT_PRODUCT(sfrVeg,LAImax)
  VegMin=DOT_PRODUCT(sfrVeg,LAImin)
  ! DO iv=ivConif,ivGrass   !Normalized LAI for vegetation
  !    VegPhen = sfr(iv+2)*lai(id-1,iv) + VegPhen
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

  RETURN

END SUBROUTINE LUMPS_QHQE
