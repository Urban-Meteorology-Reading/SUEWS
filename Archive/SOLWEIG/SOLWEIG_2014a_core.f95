! This is the core function of the SOLWEIG model
! 2013-10-27
! Fredrik Lindberg, fredrikl@gvc.gu.se
! G�teborg Urban Climate Group
! Gothenburg University
! Last modified:
!   LJ 27 Jan 2016  Renoval of tabs and fixing int-real conversions
!   HCW 02 Dec 2014 DEG2RAD and RAD2DEG commented out as now defined in AllocateArray


SUBROUTINE Solweig_2014a_core(iMBi)

  USE matsize
  USE solweig_module
  USE data_in
  USE gis_data
  USE time
  USE allocateArray

  IMPLICIT NONE
  INTEGER         :: DOY,hour,first,second,j,dfm,iMBi!,ith!onlyglobal,usevegdem,x,y,i
  REAL(KIND(1d0)) :: albedo_b,albedo_g,eground,ewall!absK,absL,,Fside,Fup
  REAL(KIND(1d0)) :: t,Tstart,height,psi!,timezone,lat,lng,alt,amaxvalue
  REAL(KIND(1d0)) :: altitude,zen!scale,azimuth,zenith
  REAL(KIND(1d0)) :: CI,CI_Tg,c,I0,Kt,Tg,Tgamp,Tw,Ktc,weight1
  REAL(KIND(1d0)) :: Ta,RH,P,radG,radD,radI,radI0!,idectime,tdectime!dectime,
  REAL(KIND(1d0)) :: corr,I0et,CIuncorr,s!,lati
  REAL(KIND(1d0)) :: SNDN,SNUP,DEC,DAYL!,timestepdec,YEAR
  REAL(KIND(1d0)) :: msteg,esky,ea
  ! Internal grids
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: tmp,Knight,svfbuveg,Tgmap0!,Tgmap
  !Search directions for Ground View Factors (GVF)
  REAL(KIND(1d0)) :: azimuthA(1:18)=(/ (j*(360.0/18.0),j=0,17) /)
  ! temporary parameters and variables for testing
  REAL(KIND(1d0)),PARAMETER   :: pi=3.141592653589793
  REAL(KIND(1d0)),PARAMETER   :: SBC=5.67051e-8
  !REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292,RAD2DEG=57.29577951 !Now defined in AllocateArray HCW 02 Dec 2014

  INTEGER:: firstTimeofDay=0  !!Needs updating for new model timestep

!!!!!! Begin program !!!!!!
  ! internal grids
  ALLOCATE(tmp(sizey,sizex))
  ALLOCATE(Knight(sizey,sizex))
  ALLOCATE(Tgmap0(sizey,sizex))
  ALLOCATE(svfbuveg(sizey,sizex))
  !allocate(Tgmap0(sizey,sizex))

  ! external grids
  IF (ALLOCATED(Kdown2d)) DEALLOCATE(Kdown2d); ALLOCATE(Kdown2d(sizey,sizex))
  IF (ALLOCATED(Kup2d)) DEALLOCATE(Kup2d); ALLOCATE(Kup2d(sizey,sizex))
  IF (ALLOCATED(Knorth)) DEALLOCATE(Knorth); ALLOCATE(Knorth(sizey,sizex))
  IF (ALLOCATED(Kwest)) DEALLOCATE(Kwest); ALLOCATE(Kwest(sizey,sizex))
  IF (ALLOCATED(Ksouth)) DEALLOCATE(Ksouth); ALLOCATE(Ksouth(sizey,sizex))
  IF (ALLOCATED(Keast)) DEALLOCATE(Keast); ALLOCATE(Keast(sizey,sizex))
  IF (ALLOCATED(Ldown2d)) DEALLOCATE(Ldown2d); ALLOCATE(Ldown2d(sizey,sizex))
  IF (ALLOCATED(Lup2d)) DEALLOCATE(Lup2d); ALLOCATE(Lup2d(sizey,sizex))
  IF (ALLOCATED(LNorth)) DEALLOCATE(LNorth); ALLOCATE(LNorth(sizey,sizex))
  IF (ALLOCATED(Lwest)) DEALLOCATE(Lwest); ALLOCATE(Lwest(sizey,sizex))
  IF (ALLOCATED(Lsouth)) DEALLOCATE(Lsouth); ALLOCATE(Lsouth(sizey,sizex))
  IF (ALLOCATED(Least)) DEALLOCATE(Least); ALLOCATE(Least(sizey,sizex))
  IF (ALLOCATED(gvf)) DEALLOCATE(gvf); ALLOCATE(gvf(sizey,sizex))
  IF (ALLOCATED(Sstr)) DEALLOCATE(Sstr); ALLOCATE(Sstr(sizey,sizex))
  IF (ALLOCATED(Tmrt)) DEALLOCATE(Tmrt); ALLOCATE(Tmrt(sizey,sizex))
  IF (ALLOCATED(shadow)) DEALLOCATE(shadow); ALLOCATE(shadow(sizey,sizex))
  IF (ALLOCATED(sos)) DEALLOCATE(sos); ALLOCATE(sos(sizey,sizex))
  IF (ALLOCATED(F_sh)) DEALLOCATE(F_sh); ALLOCATE(F_sh(sizey,sizex))
  IF (ALLOCATED(svfalfa)) DEALLOCATE(svfalfa); ALLOCATE(svfalfa(sizey,sizex))

  ! These variables should change name in the future...
  P=Press_hPa
  Ta=Temp_C
  RH=avrh
  radG=avkdn
  DOY=id
  hour=it
  height=heightgravity
  psi=trans
  albedo_b=alb(2) ! taken from Bldg (FunctionalTypes)
  albedo_g=alb(1) ! taken from Paved (FunctionalTypes)
  ewall=emis(2) ! taken from Bldg (FunctionalTypes)
  eground=emis(1) ! taken from Paved (FunctionalTypes)
  radD=kdiff
  radI=kdir

  ! Transmissivity of shortwave radiation through vegetation based on decid LAI
  IF (it==firstTimeofDay) THEN
     trans=TransMin+(LAImax(2)-LAI(id-1,2))*transperLAI
  ENDIF

  ! Radiation sensor setup offset in degrees
  t=0

  !Surface temperature difference at sunrise
  Tstart=3.41

  !Initialization of maps
  Knight=0.0

  tmp=1-(svf+svfveg-1)
  WHERE (tmp<=0) tmp=0.000000001 ! avoiding log(0)
  svfalfa=ASIN(EXP(LOG(tmp)/2))

  !Parameterization for Lup
  first=INT(ANINT(height)) !Radiative surface influence, Rule of thumb by Schmid et al. (1990).
  IF (first==0) THEN
     first=1
  END IF
  second=INT(ANINT(height*20))

  ! SVF combines for buildings and vegetation
  svfbuveg=(svf-(1-svfveg)*(1-psi))

  ! Sun position related things
  CALL DAYLEN(DOY,lat,DAYL,DEC,SNDN,SNUP)
  zen=zenith_deg*DEG2RAD
  altitude=90-zenith_deg

  !Determination of clear-sky emissivity from Prata (1996)
  ea=6.107*10**((7.5*Ta)/(237.3+Ta))*(RH/100)!Vapor pressure
  msteg=46.5*(ea/(Ta+273.15))
  esky=(1-(1+msteg)*EXP(-((1.2+3.0*msteg)**0.5)))-0.04

!!! DAYTIME !!!
  IF (altitude>0) THEN

     !Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
     !factor for low sun elevations after Lindberg et al. (2008)
     CALL clearnessindex_2013b(zen,DOY,Ta,RH/100,radG,lat,P,I0,CI,Kt,I0et,CIuncorr)
     IF (CI>1) CI=1
     CIlatenight=CI

     !Estimation of radD and radI if not measured after Reindl et al. (1990)
     IF (onlyglobal == 1) THEN
        CALL diffusefraction(radG,altitude,Kt,Ta,RH,radI,radD)
     END IF

     !Shadow images
     IF (usevegdem==1) THEN ! with vegetation
        CALL shadowingfunction_20(azimuth,altitude,scale,amaxvalue)
        shadow=(sh-(1-vegsh)*(1-psi))
     ELSE ! without vegetation
        CALL shadowingfunction_10(azimuth,altitude,scale)
        vegsh = 1.0D0
        shadow=sh
     END IF

     !Ground View Factors based on shadow patterns and sunlit walls
     gvf=0.0D0
     CALL wallinsun_veg(azimuth)
     DO j=1,SIZE(azimuthA)
        CALL sunonsurface_veg(azimuthA(j),scale,first,second,psi)
        gvf=gvf+sos
     END DO
     gvf=gvf/SIZE(azimuthA)+(buildings*(-1)+1)

     !Building height angle from svf
     CALL cylindric_wedge(zen) !Fraction shadow on building walls based on sun altitude and svf
     !F_sh(isnan(F_sh))=0.5 FIXTHIS

!!! Calculation of shortwave daytime radiative fluxes !!!
     Kdown2d=radI*shadow*SIN(altitude*DEG2RAD)+radD*svfbuveg+ &
          radG*albedo_b*(1-svfbuveg)*(1-F_sh) !*sin(altitude(i)*(pi/180));

     Kup2d=albedo_g*(radI*gvf*SIN(altitude*DEG2RAD)+radD*svfbuveg+ &
          radG*albedo_b*(1-svfbuveg)*(1-F_sh))

     CALL Kside_veg_v24(radI,radG,radD,azimuth,altitude,psi,t,albedo_b)

!!!! Surface temperature parameterisation during daytime !!!!
     dfm=ABS(172-DOY) !Day from midsommer
     Tgamp=0.000006*dfm**3-0.0017*dfm**2+0.0127*dfm+17.084+Tstart !sinus function for daily surface temperature wave
     !Tg=Tgamp*sin(((hour-rise)/(15-rise))*pi/2)-Tstart ! check if this should be 15 (3 pm)
     Tg=Tgamp*SIN(((dectime-DOY-SNUP/24)/(15./24-SNUP/24))*pi/2)-Tstart !new sunrise time 2014a
     IF (Tg<0) Tg=0! temporary for removing low Tg during morning 20140513, from SOLWEIG1D
     s=(dectime-DOY-SNUP/24)
     !New estimation of Tg reduction for non-clear situation based on Reindl et al. 1990
     Ktc=1.0
     CALL diffusefraction(I0,altitude,Ktc,Ta,RH,radI0,s)
     corr=0.1473*LOG(90.-(zen*RAD2DEG))+0.3454 ! 20070329 temporary correction of latitude from Lindberg et al. 2008
     CI_Tg=(radI/radI0)+(1.-corr)

     IF (CI_Tg>1) CI_Tg=1  !!FIX THIS?? .and. CI_Tg<inf then CI_Tg=1

     Tg=Tg*CI_Tg !new estimation
     Tw=Tg

!!!! Lup, daytime !!!!
     !Surface temperature wave delay - new as from 2014a
     Tgmap0=gvf*Tg+Ta ! current timestep
     IF (firstdaytime==1) THEN !"first in morning"
        Tgmap1=Tgmap0
     END IF
     IF (timeadd>=(59/1440.)) THEN !more or equal to 59 min
        weight1=EXP(-33.27*timeadd) !surface temperature delay function - 1 step
        Tgmap1=Tgmap0*(1-weight1)+Tgmap1*weight1
        Lup2d=SBC*eground*((Tgmap1+273.15)**4)
        IF (timestepdec>(59/1440.)) THEN
           timeadd=timestepdec
        ELSE
           timeadd=0
        END IF
     ELSE
        timeadd=timeadd+timestepdec
        weight1=EXP(-33.27*timeadd) !surface temperature delay function - 1 step
        Lup2d=SBC*eground*((Tgmap0*(1-weight1)+Tgmap1*weight1+273.15)**4)
     END IF
     firstdaytime=0

  ELSE !!!!!!! NIGHTTIME !!!!!!!!

     !Nocturnal cloud fraction from Offerle et al. 2003
     IF (dectime<(DOY+0.5) .AND. dectime>DOY .AND. altitude<1.0) THEN  !! THIS NEED SOME THOUGHT 20150211
        !j=0
        !do while (dectime<(DOY+SNUP/24))
        !!    call ConvertMetData(ith+j) ! read data at sunrise ??
        !    j=j+1
        !end do
        !call NARP_cal_SunPosition(year,dectime,timezone,lat,lng,alt,azimuth,zenith_deg)!this is not good
        !zen=zenith_deg*DEG2RAD
        !call clearnessindex_2013b(zen,DOY,Temp_C,RH/100,avkdn,lat,Press_hPa,I0,CI,Kt,I0et,CIuncorr)
        !!call ConvertMetData(ith) ! read data at current timestep again ??
        !call NARP_cal_SunPosition(year,dectime,timezone,lat,lng,alt,azimuth,zenith_deg)!this is not good

        CI=1.0
     ELSE
        IF (SolweigCount==1) THEN
           CI=1.0
        ELSE
           CI=CIlatenight
        END IF
     END IF

     Tw=0.0
     Tg=0.0

     !Nocturnal Kfluxes set to 0
     Kdown2d=0.0
     Kwest=0.0
     Kup2d=0.0
     Keast=0.0
     Ksouth=0.0
     Knorth=0.0
     shadow=0.0
     gvf=0.0

!!! Lup !!!
     Lup2d=SBC*eground*((Knight+Ta+Tg+273.15)**4)
     firstdaytime=1
  END IF

!!! Ldown !!!
  IF (SOLWEIG_ldown==1) THEN   ! fixed for non-clear situations 20140701
     Ldown2d=(svf+svfveg-1)*esky*SBC*((Ta+273.15)**4)+(2-svfveg-svfaveg)*ewall*SBC*((Ta+273.15)**4)+ &
          (svfaveg-svf)*ewall*SBC*((Ta+273.15+Tw)**4)+(2-svf-svfveg)*(1-ewall)*esky*SBC*((Ta+273.15)**4) !Jonsson et al. (2006)
     !Ldown2d=Ldown2d-25 ! Shown by Jonsson et al. (2006) and Duarte et al. (2006) ! Removed as from 20140701
     IF (CI>1) CI=1  !!FIX THIS?? .and. CI<inf) CI=1
     !if (CI < 0.95) then !non-clear conditions
     c=1-CI
     !Ldown2d=Ldown2d*(1-c)+c*SBC*((Ta+273.15)**4)
     Ldown2d=Ldown2d*(1-c)+c*((svf+svfveg-1)*SBC*((Ta+273.15)**4)+(2-svfveg-svfaveg)*ewall*SBC*((Ta+273.15)**4)+ &
          (svfaveg-svf)*ewall*SBC*((Ta+273.15+Tw)**4)+(2-svf-svfveg)*(1-ewall)*SBC*((Ta+273.15)**4))
     !end if
  ELSE
     Ldown2d=(svf+svfveg-1)*ldown+(2-svfveg-svfaveg)*ewall*SBC*((Ta+273.15)**4)+ &
          (svfaveg-svf)*ewall*SBC*((Ta+273.15+Tw)**4)+(2-svf-svfveg)*(1-ewall)*ldown !Jonsson et al. (2006)
  END IF

!!! Lside !!!
  CALL Lside_veg_v2(altitude,Ta,Tw,SBC,ewall,esky,t,CI)!azimuth,

!!! Calculation of radiant flux density and Tmrt !!!
  Sstr=absK*(Kdown2d*Fup+Kup2d*Fup+Knorth*Fside+Keast*Fside+Ksouth*Fside+Kwest*Fside) &
       +absL*(Ldown2d*Fup+Lup2d*Fup+Lnorth*Fside+Least*Fside+Lsouth*Fside+Lwest*Fside)
  Tmrt=SQRT(SQRT((Sstr/(absL*SBC))))-273.2

  IF (SOLWEIGpoi_out==1) THEN
     dataOutSOL(SolweigCount,1:4,iMBi)=[iy,id,it,imin]
     dataOutSOL(SolweigCount,5,iMBi)=dectime
     dataOutSOL(SolweigCount,6:ncolumnsdataOutSOL,iMBi)=[&
          azimuth,altitude,radG,radD,radI,&
          Kdown2d(row,col),Kup2d(row,col),Ksouth(row,col),Kwest(row,col),Knorth(row,col),Keast(row,col),&
          Ldown2d(row,col),Lup2d(row,col),Lsouth(row,col),Lwest(row,col),Lnorth(row,col),Least(row,col),&
          Tmrt(row,col),I0,CI,gvf(row,col),shadow(row,col),svf(row,col),svfbuveg(row,col),Ta,Tg]
  END IF

  CALL SaveGrids

  DEALLOCATE(tmp)
  DEALLOCATE(Knight)
  DEALLOCATE(Tgmap0)
  DEALLOCATE(svfbuveg)

  ! external grids
  DEALLOCATE(Kdown2d)
  DEALLOCATE(Kup2d)
  DEALLOCATE(Knorth)
  DEALLOCATE(Kwest)
  DEALLOCATE(Ksouth)
  DEALLOCATE(Keast)
  DEALLOCATE(Ldown2d)
  DEALLOCATE(Lup2d)
  DEALLOCATE(Lnorth)
  DEALLOCATE(Lwest)
  DEALLOCATE(Lsouth)
  DEALLOCATE(Least)

END SUBROUTINE Solweig_2014a_core
