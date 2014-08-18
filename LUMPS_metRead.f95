subroutine MetRead(i)
!Reading met data
!Input: hour number i
!Code changed in Feb 2012 (LJ). Input fluxes qh and qe changed _obs as well as qn1_obs ending

  use data_in
  use gis_data
  use time
  use sues_data
  use defaultNotUsed
  implicit none
  integer::lfn=1, i

  real(kind(1D0)):: dir30
      
  finish=.false.

  
  if (InputMetFormat==0) then !This is the default format using LUMPS only
    
  	READ(lfn,*,iostat=iostat_var)id,it,dectime,qn1_obs,avu1,avrh,& 
             Temp_C,dir30,Pres_kPa,Precip_hr,avkdn,snow_obs,ldown_obs,fcld_obs
    
    qf=-999
    qs=-999
    qh_obs=-999
    qe_obs=-999
    xsmd=-99999
    
  elseif (InputMetFormat==10) then !SUEWS reading
      READ(lfn,*,iostat=iostat_var)id,it,dectime,qn1_obs,qh_obs,qe_obs,qs,qf,avu1,avrh,&
               Temp_C,Pres_kPa,Precip_hr,avkdn,snow_obs,ldown_obs,fcld_obs,&
               wuh,xsmd,lai_hr,kdiff,kdir,wdir
               
      ! check with LJ -- lai_hr only one veg type
  
      !Calculate soil moisture deficits from either volumetric or gravimetric soilstates
      if (smd_choice==1.and.xsmd/=-999) then !Soil moisture - volumetric 
         xsmd=(SmCap-xsmd)*SoilDepthMeas*SoilRocks
      elseif (smd_choice==2.and.xsmd/=-999) then !Soil moisture -gravimetric
         xsmd=(SmCap-xsmd)*SoilDensity*SoilDepthMeas*SoilRocks
      else
         xsmd=-999  
      endif
             
  else
     write(12,*)'HeaderInput.nml error InputMetFormat not usable ',InputmetFormat
  endif
  
  !################Meteorological variables reading done###########################
  Press_hPa=Pres_kPa*10. ! convert to hPa

  if (qs==-999.0) qs=defaultQs  
       
  if (AnthropHeatChoice==0.and.qf==-999)then
      call ErrorHint(30,'subroutine MetRead: [Qf default value going to be used],qf,dectime',qf,dectime,notUsedI)
      qf=defaultQf
  endif

  if(it==24) then
     id=id+1
     it=it+1
  endif

  if (int(dectime)/=id) then !Added my LJ in Feb 2013
    call ErrorHint(35,'Met Data: decimal time does not match with day of year',real(id,kind(1d0)),dectime,notUsedI)
  endif
    
  
  if(iostat_var<0)THEN
     iostat_var=0
     CLOSE(lfn)
     finish=.TRUE.
     RETURN
  ENDIF

  qual=0
  if(AvKdn<0) then
    call ErrorHint(27,'Met Data: avKdn - needed for Surf. resistance, If present, check file not tab delimited',&
                      avkdn,dectime, notUsedI)
     !sg removed this is causing the problems with resistances 
     !  AvKdn=0 !Solar radiation cannot be lower than 1
  endif
 
  
  if((ldown_option==1).and.(ldown_obs<0))then
     call ErrorHint(27,'Met Data: LWdn (ldown_obs) - impact Q* calc', ldown_obs,dectime, notUsedI)
    
  elseif(ldown_option==2) then
     if(fcld_obs==-999.0.or.fcld_obs<0.or.fcld_obs>1) then
        call ErrorHint(27,'Met Data: flcd_obs - impacts LW & Q* radiation', fcld_obs,dectime, notUsedI)
     endif	
  endif
  
  if(qn1_obs==-999.and.NetRadiationChoice==0) then  !If measured Q* is used and it is -999 
     call ErrorHint(27,'Met Data: Q* - will impact everything', qn1,dectime, notUsedI)
  endif
    
  if(avu1<=0) then !If wind speed is negative
    call ErrorHint(27,'Met Data: avU1 - impacts aeroydnamic resistances', avU1,dectime, notUsedI)
  endif
     

  if(Temp_C<-50.or.Temp_C>60)then !If temperature unrealistic
    call ErrorHint(27,'Met Data: Temp_C - beyond what is expected', Temp_C,dectime, notUsedI)
  endif
  
  if(avrh>100.or.avrh<1)then !If relative humidity larger than 100%
      call ErrorHint(27,'Met Data: avRH - beyond what is expected', avRH,dectime, notUsedI)
  endif
  
  if(Pres_kPa<90)then  !If pressure too low
    call ErrorHint(27,'Met Data: Pres_kPa - too low - this could be fixed in model',Pres_kPa ,dectime, notUsedI)
  endif

  if (Precip_hr<0) then  !If rain in negative, set it to zero
     call ErrorHint(27,'Met Data: Precip_hr - less than 0',Precip_hr ,dectime, notUsedI)
  endif

  if (snow_obs==NAN) snow_obs=0
    
  if (snowUse==0.and.(snow_obs<0.or.snow_obs>1)) then
     call ErrorHint(27,'Met Data: snow not between [0  1]',snow_obs ,dectime, notUsedI)
  endif

  if (xsmd<0.and.smd_choice==1) then  !If soil moisture deficit is zero
    call ErrorHint(27,'Met Data: xsmd - less than 0',xsmd ,dectime, notUsedI)
  endif
   
  !Check initial conditions
  if (i==1) call CheckInitial 

  RETURN	
 END SUBROUTINE MetRead

