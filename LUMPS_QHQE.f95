 subroutine LUMPS_QHQE
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

  use allocateArray
  use data_in
  use defaultNotUsed
  use gis_data      
  use moist
  use snowMod
  use sues_data
  use time
  use VegPhenogy
  
  IMPLICIT NONE
 
  integer::iv                                 !,start
  real (kind(1d0))::VegPhen,VegMax,VegMin,&   !Vegetation phenology for LUMPS
                    slope_svp, slopeIce_svp,& !Slope of the saturation vapour pressure curve above watre and ice
                    psyc_s,psyc_const,&       !Psychometric constant
                    alpha_sl,alpha_in,&    	  !Parameters used in LUMPS QH and QE calculations
                    beta                      !Beta parameter used in LUMPS QH and QE calculations [W m-2]

  ! Calculate slope of the saturation vapour pressure vs air temp.
  s_hPa=slope_svp(Temp_C) 
  psyc_hPa=psyc_const(avcp,Press_hPa,lv_J_kg)
  psyc_s=psyc_hPa/s_hPa

  !Calculate also sublimation ones if snow calculations are made. 
  !Used also for LUMPS
  if (snowUse==1) then
      if (Temp_C<=0) then
         sIce_hpa=slopeIce_svp(Temp_C) 
      else
         sIce_hpa=slope_svp(Temp_C)
      endif
      psyc_s=psyc_hPa/sIce_hPa   !Psychometric constant divided by the slope
  endif

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
  VegPhen=0
  VegMax=0
  VegMin=0
  do iv=ivConif,ivGrass   !Normalized LAI for vegetation
     VegPhen = sfr(iv+2)*lai(id-1,iv) + VegPhen
     VegMax  = sfr(iv+2)*LAImax(iv) + VegMax
     VegMin  = sfr(iv+2)*LAImin(iv) + VegMin
  enddo
  
  if(VegMax<=0.01000) then   !If max vegetation is very small, TempVeg = 0;
     TempVeg=0
  else
     VegPhenLumps=(VegPhen)/(VegMax)
     TempVeg=Veg_Fr*VegPhenLumps   !Now this is veg_fraction in general
  endif
 
  if (TempVeg>0.9000) then   !If vegetation fraction is larger than 0.7 (0.9?)
     beta = (20-3)*TempVeg+3
     alpha_qhqe=TempVeg*0.8+0.2
  else  
     beta=3
     if(veg_type==1) then   !Area vegetated, including bare soil and water
        alpha_sl=0.686 
        alpha_in=0.189
     elseif(veg_type==2) then   !Area irrigated vegetation
        alpha_sl=0.610 
        alpha_in=0.222
     endif
     alpha_qhqe=TempVeg*alpha_sl+alpha_in
  endif
     
  ! Calculate the actual heat fluxes
  H_mod= ((1-alpha_qhqe)+psyc_s)/(1+psyc_s)*(qn1+qf-qs-Qm)-beta   !Eq 3, Grimmond & Oke (2002)
  E_mod= (alpha_qhqe/(1+psyc_s)*(qn1+qf-qs-Qm))+beta              !Eq 4, Grimmond & Oke (2002)
 
  return

 end subroutine LUMPS_QHQE



