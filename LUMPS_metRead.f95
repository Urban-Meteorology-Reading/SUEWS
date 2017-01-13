!Reading one line of meteorological forcing data in.
!Latest change:
!  Feb 2012, LJ:  Input fluxes qh and qe changed _obs as well as qn1_obs ending
!  Oct 2014, LJ:  Variables changed only be used in this part of code and these are passed to calling
!                 function in MetArray.
!  Jan 2015, HCW: Precip_hr, wuh and lai_hr changed for generic timesteps
!  Jan 2016, LJ:  Removal of tabs  
! To Do: 
!       - Check observed SM calculation
!---------------------------------------------------------------------------------------------------
  subroutine MetRead(MetArray,InputmetFormat,ldown_option,NetRadiationMethod,&
               snowUse,SMDMethod,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)

  use defaultNotUsed

  IMPLICIT NONE

  !INPUT
  real (kind(1d0)),dimension(24)::MetArray !Array leaving the subroutine within 
                                           !each INTERVAL (defined in RunControl.nml)
                                           ! - Met data now provided at a resolution of tstep, HCW Jan 2015
                                           ! so MetArray could be bypassed??
  
  real (kind(1d0))::SmCap,&
                    SoilDepthMeas,&        !Measured soil depth
                    SoilRocks,&            !Rocks on ground
                    SoilDensity            !Density of soil

  integer::InputmetFormat,&     !Format of the meteorological forcing file
           ldown_option,&       !Method of calculating Ldown
           NetRadiationMethod,& !Method of calculating Q*
           SMDMethod,&         !Method of measured soil moisture
           snowUse

  ! Variables read in
  real (kind(1d0))::avkdn,&     !Average downwelling shortwave radiation
                    avrh,&      !Average relative humidity
                    avu1,&      !Average wind speed
                    dectime,&   !Decimal time
                    fcld_obs,&  !Cloud fraction observed
                    iy,&        !Year
                    id,&        !Day
                    it,&        !Hour
                    imin,&      !Minute
                    kdiff,&     !Diffuse shortwave radiation
                    kdir,&      !Direct shortwave radiation
                    lai_obs,&   !Overall LAI of the study area
                    ldown_obs,& !Downwelling longwave radiation
                    Precip,& !Rainfall [mm]
                    Pres_hPa,&  !Station air pressure in hPa
                    Pres_kPa,&  !Station air pressure in kPa
                    snow_obs,&  !Observed surface fraction of snow (between 0 and 1)
                    qe_obs,&    !Observed latent heat flux
                    qf_obs,&    !Observed antrhropogeni heat flux
                    qh_obs,&    !Observed sensible heat flux
                    qn1_obs,&   !Observed net all-wave radiation
                    qs_obs,&    !Observed storage heat flux
                    Temp_C,&    !Air temperature
                    wdir,&      !Wind direction
                    wu_m3,&     !Water use provided in met forcing file [m3]
                    xsmd        !Measured soil moisture deficit

  integer::iostat_var,lfn=1
 
  !-----------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------
  
  if (InputMetFormat==0) then   !Default format using LUMPS only

    READ(lfn,*,iostat=iostat_var)iy,id,it,imin,qn1_obs,avu1,avrh,&
             Temp_C,wdir,Pres_kPa,Precip,avkdn,snow_obs,ldown_obs,fcld_obs

    !Set other variables needed while running SUEWS to zero
    qf_obs=NaN
    qs_obs=NaN
    qh_obs=NaN
    qe_obs=NaN
    xsmd=-99999
    kdiff=NaN
    kdir=NaN
    wdir=NaN

  elseif (InputMetFormat==10) then !SUEWS reading
      READ(lfn,*,iostat=iostat_var) iy,id,it,imin,qn1_obs,qh_obs,qe_obs,qs_obs,qf_obs,avu1,avrh,&
                                    Temp_C,Pres_kPa,Precip,avkdn,snow_obs,ldown_obs,fcld_obs,&
                                    wu_m3,xsmd,lai_obs,kdiff,kdir,wdir
                                   

  !write(*,*) 'In LUMPS_MetRead (1)'
  !write(*,*) 'imin',imin             
  !write(*,*) 'it',it             
  !write(*,*) 'id',id
  !write(*,*) 'iy',iy                        
                                   

      !Calculate observed soil moisture deficits from either volumetric or gravimetric soilstates
      if (SMDMethod==1.and.xsmd/=-999) then !Soil moisture - volumetric
         xsmd=(SmCap-xsmd)*SoilDepthMeas*SoilRocks
      elseif (SMDMethod==2.and.xsmd/=-999) then !Soil moisture -gravimetric
         xsmd=(SmCap-xsmd)*SoilDensity*SoilDepthMeas*SoilRocks
      else
         xsmd=-999
      endif
             
  else
     call ErrorHint(55,'RunControl.nml, InputMetFormat not usable.',notUsed,notUsed,InputmetFormat)
  endif
  
  !===============Meteorological variables reading done==========================
  Pres_hPa=Pres_kPa*10. ! convert to hPa

  !If hour is 23, change this to following day
  if(it==24) then
     id=id+1
     it=it+1
  endif

  if(iostat_var<0)THEN
     iostat_var=0
     CLOSE(lfn)
     RETURN
  ENDIF

  if(AvKdn<0) then
    call ErrorHint(27,'Met Data: avKdn - needed for Surf. resistance, If present, check file not tab delimited',&
                      avkdn,dectime,notUsedI)
     !sg removed this is causing the problems with resistances 
     !  AvKdn=0 !Solar radiation cannot be lower than 1
  endif

  if((ldown_option==1).and.(ldown_obs<0))then
     call ErrorHint(27,'Met Data: LWdn (ldown_obs) - impact Q* calc',ldown_obs,dectime,notUsedI)
    
  elseif(ldown_option==2) then
     if(fcld_obs==-999.0.or.fcld_obs<0.or.fcld_obs>1) then
        call ErrorHint(27,'Met Data: flcd_obs - impacts LW & Q* radiation',fcld_obs,dectime,notUsedI)
     endif  
  endif
  
  if(qn1_obs==-999.and.NetRadiationMethod==0) then  !If measured Q* is used and it is -999 
     call ErrorHint(27,'Met Data: Q* - will impact everything', qn1_obs,dectime, notUsedI)
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

  if (Precip<0) then  !If rain in negative, set it to zero
     call ErrorHint(27,'Met Data: Precip - less than 0',Precip ,dectime, notUsedI)
  endif

  if (snow_obs==NAN) snow_obs=0
    
  if (snowUse==0.and.(snow_obs<0.or.snow_obs>1)) then
     call ErrorHint(27,'Met Data: snow not between [0  1]',snow_obs ,dectime, notUsedI)
  endif

  if (xsmd<0.and.SMDMethod==1) then  !If soil moisture deficit is zero
    call ErrorHint(27,'Met Data: xsmd - less than 0',xsmd ,dectime, notUsedI)
  endif
   
  !Create an array to be printed out.
  MetArray(1:24)=(/iy,id,it,imin,qn1_obs,qh_obs,qe_obs,qs_obs,qf_obs,avu1,&
                   avrh,Temp_C,Pres_hPa,Precip,avkdn,snow_obs,ldown_obs,&
                   fcld_obs,wu_m3,xsmd,lai_obs,kdiff,kdir,wdir/)
                   
  !write(*,*) 'In LUMPS_MetRead (2)'
  !write(*,*) 'imin',imin             
  !write(*,*) 'it',it             
  !write(*,*) 'id',id
  !write(*,*) 'iy',iy        
  
  RETURN

 END SUBROUTINE MetRead

