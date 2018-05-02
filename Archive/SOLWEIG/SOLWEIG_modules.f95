!Last modified: LJ 27 Jan 2016 - Removal of tabs
 module matsize

    IMPLICIT NONE

    integer    :: sizex,sizey ! number of rows and cols of grid
    real(kind(1d0)), allocatable, dimension(:,:)  :: a,sh,vbshvegsh,vegsh
    real(kind(1d0)), allocatable, dimension(:,:)  :: bush,vegdem,vegdem2,tempgrid
    real(kind(1d0)), allocatable, dimension(:,:)  :: buildings,svf,svfE,svfS,svfW,svfN
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
    real(kind(1d0)), allocatable, dimension(:,:)  :: Kdown2d,Keast,Knorth,Ksouth,Kup2d,Kwest
    real(kind(1d0)), allocatable, dimension(:,:)  :: Ldown2d,Least,Lnorth,Lsouth,Lup2d,Lwest
    real(kind(1d0)), allocatable, dimension(:,:)  :: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfalfa,sos,Tgmap1
    real(kind(1d0)), allocatable, dimension(:,:)  :: viktveg,viktsky,viktrefl,viktwall,savegrid

 end module matsize

!svfvegbu,,viktaveg,viktonlywall,svfalfaE,svfalfaN,svfalfaS,svfalfaW,alfaB,betaB,betasun,Lground,Lrefl,Lsky,Lveg,Lwallsh,Lwallsun
! alfa,xa,ha,hkil,Ai,phi,qa,Za,tempsh,tempbu,tempbub,tempwallsun,weightsumwall,weightsumsh,gvf1,gvf2
!!,svfbuveg
!Lnight,,tempvegdem,tempvegdem2,fabovea,gabovea,tempbush,firstvegdem,vegsh2,tempgrid
!svfviktbuveg
!,stopbuild,stopveg,g,bushplant,Tgmap
!,temp,f,tmp,Knight,,tempb
!Tgmapgvf

module solweig_module
    IMPLICIT NONE

    real(kind(1d0)) :: timestepdec,& !time step in decimal time
                       CIlatenight,&
                       timeadd,&
                       firstdaytime,& ! if new day starts, =1 else =0
                       Fside,& ! fraction of a person seen from each cardinal point
                       Fup,& ! fraction of a person seen from down and up
                       scale,&
                       amaxvalue,&
                       trans,&
                       transperLAI,&
                       xllcorner,&
                       yllcorner,&
                       NoData,&
                       cellsize
    integer         :: SolweigCount

    integer    :: sizex,sizey ! number of rows and cols of grid
    real(kind(1d0)), allocatable, dimension(:,:)  :: a,sh,vbshvegsh,vegsh
    real(kind(1d0)), allocatable, dimension(:,:)  :: bush,vegdem,vegdem2,tempgrid
    real(kind(1d0)), allocatable, dimension(:,:)  :: buildings,svf,svfE,svfS,svfW,svfN
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
    real(kind(1d0)), allocatable, dimension(:,:)  :: Kdown2d,Keast,Knorth,Ksouth,Kup2d,Kwest
    real(kind(1d0)), allocatable, dimension(:,:)  :: Ldown2d,Least,Lnorth,Lsouth,Lup2d,Lwest
    real(kind(1d0)), allocatable, dimension(:,:)  :: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
    real(kind(1d0)), allocatable, dimension(:,:)  :: svfalfa,sos,Tgmap1
    real(kind(1d0)), allocatable, dimension(:,:)  :: viktveg,viktsky,viktrefl,viktwall,savegrid
CONTAINS
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

!>
 subroutine clearnessindex_2013b(zen,jday,Ta,RH,radG,lat,P,I0,CI,Kt,I0et,CIuncorr)
 !Last modified:
 !LJ 27 Jan 2016 - Removal of tabs

 implicit none
    ! Use somemodule

    real(kind(1d0))                 :: CI,CIuncorr,I0,I0et,Kt
    integer                         :: jday
    real(kind(1d0))                 :: P,radG,RH,Ta,zen,lat,iG,Itoa
    real(kind(1d0)), dimension(4)   :: G
    real(kind(1d0)), parameter          :: pi=3.141592653589793

    ! Variable declarations
    real*8 :: a2      !>
    !logical :: b      !>
    real*8 :: b2      !>
    real*8 :: corr      !>
    real*8 :: D      !>
    real*8 :: m      !>
    real*8 :: Tar      !>
    real*8 :: Td      !>
    real*8 :: Trpg      !>
    real*8 :: Tw      !>
    real*8 :: u      !>
    !

    ! Clearness Index at the Earth's surface calculated from Crawford and Duchon 1999

    if (P==-999) then
        p=1013 !Pressure in millibars
    else
        p=P*10 !Convert from hPa to millibars
    end if
    Itoa=1370 !Effective solar constant
    !call solar_ESdist(jday,D)
    call sun_distance(jday,D)
    !D=sun_distance(jday) !irradiance differences due to Sun-Earth distances
    m=35.*cos(zen)*((1224.*(cos(zen)**2)+1.)**(-1./2.)) !optical air mass at p=1013
    Trpg=1.021-0.084*(m*(0.000949*p+0.051))**0.5 !Transmission coefficient for Rayliegh scattering and permanent gases

    ! empirical constant depending on latitude
    if (lat<10) then
    G=(/3.37,2.85,2.80,2.64/)
    else if (lat>=10 .and. lat<20) then
    G=(/2.99,3.02,2.70,2.93/)
    else if (lat>=20 .and. lat<30) then
    G=(/3.60,3.00,2.98,2.93/)
    else if (lat>=30 .and. lat<40) then
    G=(/3.04,3.11,2.92,2.94/)
    else if (lat>=40 .and. lat<50) then
    G=(/2.70,2.95,2.77,2.71/)
    else if (lat>=50 .and. lat<60) then
    G=(/2.52,3.07,2.67,2.93/)
    else if (lat>=60 .and. lat<70) then
    G=(/1.76,2.69,2.61,2.61/)
    else if (lat>=70 .and. lat<80) then
    G=(/1.60,1.67,2.24,2.63/)
    else if (lat>=80 .and. lat<90) then
    G=(/1.11,1.44,1.94,2.02/)
    end if
    if (jday > 335 .or. jday <= 60) then
        iG=G(1)
    else if (jday > 60 .and. jday <= 152) then
        iG=G(2)
    else if (jday > 152 .and. jday <= 244) then
        iG=G(3)
    else if (jday > 244 .and. jday <= 335) then
        iG=G(4)
    end if
    !dewpoint calculation
    a2=17.27
    b2=237.7
    Td=(b2*(((a2*Ta)/(b2+Ta))+log(RH)))/(a2-(((a2*Ta)/(b2+Ta))+log(RH)))
    Td=(Td*1.8)+32 !Dewpoint (�F)
    u=exp(0.1133-log(iG+1)+0.0393*Td) !Precipitable water
    Tw=1-0.077*((u*m)**0.3) !Transmission coefficient for water vapor
    Tar=0.935**m !Transmission coefficient for aerosols

    I0=Itoa*cos(zen)*Trpg*Tw*D*Tar

    !!! This needs to be checked !!!
    !b=I0==abs(zen)>pi/2
    !I0(b==1)=0
    !clear b
    !if (not(isreal(I0))) then
    !    I0=0
    !end if

    corr=0.1473*log(90-(zen/pi*180))+0.3454 ! 20070329

    CIuncorr=radG/I0
    CI=CIuncorr+(1-corr)
    I0et=Itoa*cos(zen)*D !extra terrestial solar radiation
    Kt=radG/I0et

    !!!! This needs to be checked !!!
    !if (isnan(CI)) then
    !    CI=NaN
    !end if
    end subroutine clearnessindex_2013b

    !===============================================================================

 subroutine sun_distance(jday,D)

 ! Calculates solar irradiance variation based on mean earth sun distance
 ! with day of year as input.
 ! Partridge and Platt, 1975

    INTEGER          ::jday
    REAL(KIND(1d0))  ::b,D

    b = 2*3.141592654*jday/365
    D = sqrt(1.00011 + 0.034221 * cos(b) + 0.001280 * sin(b) + 0.000719 * cos(2*b) + 0.000077 * sin(2*b))

 end subroutine sun_distance

subroutine cylindric_wedge(zen)

use matsize

! Fraction of sunlit walls based on sun altitude and svf wieghted building angles

    implicit none

    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0))         :: zen      !>
    real(kind(1d0))         :: beta      !>
    real(kind(1d0)),allocatable,dimension(:,:) :: alfa,xa,ha,hkil,ba
    real(kind(1d0)),allocatable,dimension(:,:) :: Ai,phi,qa,Za
    real(kind(1d0)),allocatable,dimension(:,:) :: ukil,Ssurf
    !real(kind(1d0)), dimension(sizey,sizex) ::
    !real(kind(1d0)), dimension(sizey,sizex) ::

    allocate(alfa(sizey,sizex))
    allocate(ba(sizey,sizex))
    allocate(ha(sizey,sizex))
    allocate(xa(sizey,sizex))
    allocate(qa(sizey,sizex))
    allocate(Za(sizey,sizex))
    allocate(phi(sizey,sizex))
    allocate(ukil(sizey,sizex))
    allocate(Ai(sizey,sizex))
    allocate(Ssurf(sizey,sizex))
    allocate(hkil(sizey,sizex))

    beta=zen
    alfa=svfalfa

    xa=1.-2./(tan(alfa)*tan(beta))
    ha=2./(tan(alfa)*tan(beta))
    ba=(1./tan(alfa))
    hkil=2.*ba*ha


    qa = 0.0D0

    where (xa<0) !qa(xa<0)=tan(beta)/2
        qa=tan(beta)/2
    end where


    Za = 0.0D0

    phi = 0.0D0

    Ai = 0.0D0

    ukil = 0.0D0
    where (xa<0)
        !Za(xa<0)=((ba(xa<0).**2)-((qa(xa<0).**2)./4)).**0.5
        Za=(ba**2-qa**2/4.)**0.5
        !phi(xa<0)=atan(Za(xa<0)./qa(xa<0))
        phi=atan(Za/qa)
        !A1(xa<0)=(sin(phi(xa<0))-phi(xa<0).*cos(phi(xa<0)))./(1-cos(phi(xa<0)))
        Ai=(sin(phi)-phi*cos(phi))/(1-cos(phi))
        !ukil(xa<0)=2*ba(xa<0).*xa(xa<0).*A1(xa<0)
        ukil=2*ba*xa*Ai
    end where

    Ssurf=hkil+ukil

    F_sh=(2*pi*ba-Ssurf)/(2*pi*ba) !Xa

    deallocate(alfa)
    deallocate(ba)
    deallocate(ha)
    deallocate(xa)
    deallocate(qa)
    deallocate(Za)
    deallocate(phi)
    deallocate(ukil)
    deallocate(Ai)
    deallocate(Ssurf)
    deallocate(hkil)

end subroutine cylindric_wedge

! This subroutine estimates diffuse and directbeam radiation according to
! Reindl et al (1990), Solar Energy 45:1
subroutine diffusefraction(radG,altitude,Kt,Ta,RH,radI,radD)
    implicit none

    real(kind(1d0))                 :: radG,altitude,Kt,Ta,RH,radD,radI,alfa
    REAL(KIND(1D0)),PARAMETER       :: DEG2RAD=0.017453292,RAD2DEG=57.29577951  !!Already defined in AllocateArray module. Delete??

    alfa=altitude*DEG2RAD

    if (Ta<=-99 .or. RH<=-99) then !.or. isnan(Ta) .or. isnan(RH)) then
        if (Kt<=0.3) then
            radD=radG*(1.020-0.248*Kt)
        else if (Kt>0.3 .and. Kt<0.78) then
            radD=radG*(1.45-1.67*Kt)
        else if (Kt>=0.78) then
            radD=radG*0.147
        end if
    else
        !RH=RH/100
        if (Kt<=0.3) then
            radD=radG*(1-0.232*Kt+0.0239*sin(alfa)-0.000682*Ta+0.0195*(RH/100))
        else if (Kt>0.3 .and. Kt<0.78) then
            radD=radG*(1.329-1.716*Kt+0.267*sin(alfa)-0.00357*Ta+0.106*(RH/100))
        else if (Kt>=0.78) then
            radD=radG*(0.426*Kt-0.256*sin(alfa)+0.00349*Ta+0.0734*(RH/100))
        end if
    end if
    radI=(radG-radD)/(sin(alfa))

    !! Corrections for low sun altitudes (20130307)
    if (radI<0) then
    radI=0
    end if

    if (altitude<1 .and. radI>radG) then
    radI=radG
    end if

    if (radD>radG) then
    radD=radG
    end if
end subroutine diffusefraction
 ! This subroutine loads a ESRIACSII grid as a 2D array
 ! Last modified:
 ! LJ 27 Jan 2016-removal of tabs
 !-----------------------------------------------------
 subroutine LoadEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)
    use matsize

    implicit none
    real(kind(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    integer                           :: col,row
    character(len=100)                :: GridPath,GridName,GridFile,n
    !real(kind(1d0)),intent(out)       :: tempgrid
    !real(kind(1d0)),dimension(sizey,sizex),intent(out)       :: tempgrid!,allocatable

    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=trim(GridPath)//trim(GridName)
    open(99,File=GridFile,status='old')

	! Read Header
    read(99,*) n,sizex
    read(99,*) n,sizey
    read(99,*) n,xllcornerlo
    read(99,*) n,yllcornerlo
    read(99,*) n,cellsizelo
    read(99,*) n,NoDatalo

    allocate(tempgrid(sizex,sizey))

	! Read Matrix
    do row=1,sizey
       read(99,*) (tempgrid(row,col),col=1,sizex)
    end do
    close(99)

    return
    end subroutine LoadEsriAsciiGrid


    ! This subroutine saves a 2D array as an ESRIACSII grid
    subroutine SaveEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)
    use matsize

    implicit none
    real(kind(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    integer                           :: col,row
    character(len=100)                :: GridPath,GridName,GridFile
    !integer                           :: sizey,sizex!,intent(in)
    !real(kind(1d0)), allocatable, dimension(:,:):: grid
   ! real(kind(1d0)),dimension(sizey,sizex)       :: grid!,allocatable

    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=trim(GridPath)//trim(GridName)
    open(94,File=GridFile,status='unknown')

    ! Read Header
    write(94,"(A5,1x,I0)") 'ncols',sizex
    write(94,"(A5,1x,I0)") 'nrows',sizey
    write(94,"(A9,1x,F0.2)") 'xllcorner',xllcornerlo
    write(94,"(A9,1x,F0.2)") 'yllcorner',yllcornerlo
    write(94,"(A8,1x,F0.2)") 'cellsize',cellsizelo
    write(94,"(A12,1x,F0.2)") 'NODATA_value',NoDatalo

	! write Matrix
    do row=1,sizey
        write(94,100) (savegrid(row,col),col=1,sizex)
    end do
    close(94)
100 format(200(f6.2,1x))

    return
 end subroutine SaveEsriAsciiGrid
 ! setup for SOLWEIG
 ! FL may 2014
 subroutine SOLWEIG_Initial
 !Last modified LJ 27 Jan 2016 - Removal of Tabs and real-int fixes

  use matsize         ! All allocatable grids and related variables used in SOLWEIG
  use InitialCond
  use allocateArray
  use data_in
  use sues_data
  use defaultNotUsed
  use InitialCond
  use solweig_module
  use time

  implicit none

 character(len=100)  :: Path,GridFile,GridFolder
 real(kind(1d0))                 :: vegmax
 character(len=100),dimension(5) :: svfname
 character(len=100),dimension(10):: svfvegname
 logical                         :: exist
 integer                         :: firstday

namelist/SOLWEIGinput/Posture,&    ! 1.Standing, 2.Sitting
    absL,&            ! Absorption coefficient of longwave radiation of a person
    absK,&            ! Absorption coefficient of shortwave radiation of a person
    heightgravity,&   ! Center of gravity for a standing person
    usevegdem,&       ! With vegetation (1)
    DSMPath,&         ! Path to DSMs
    DSMname,&         ! Ground and building DSM
    CDSMname,&        ! Canopy DSM
    TDSMname,&        ! Trunk zone DSM
    TransMin,&        ! Transmissivity of K through decidious vegetation (leaf on)
    TransMax,&        ! Transmissivity of K through decidious vegetation (leaf off)
    SVFPath,&         ! Path to SVFs
    SVFsuffix,&       !
    buildingsname,&   ! Boolean matrix for locations of building pixels
    row,&             ! X coordinate for point of interest
    col,&             ! Y coordinate for point of interest
    onlyglobal,&      ! if no diffuse and direct, then =1
    SOLWEIGpoi_out,&  ! write output variables at point of interest
    Tmrt_out,&        ! write Tmrt grid to file
    Lup2d_out,&       ! write Lup grid to file
    Ldown2d_out,&     ! write Ldown grid to file
    Kup2d_out,&       ! write Kup grid to file
    Kdown2d_out,&     ! write Kdown grid to file
    GVF_out,&         ! write GroundViewFactor grid to file
    SOLWEIG_ldown,&   ! 1= use SOLWEIG code to estimate Ldown, 0=use SEUWS
    OutInterval,&     ! Output interval in minutes
    RunForGrid        ! If only one grid should be run. All grids -999

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Read in the SOLWEIGinput.nml file
    open(52,file=trim(FileInputPath)//'SOLWEIGinput.nml',err=274,status='old')
    read(52,nml=SOLWEIGinput)
    close(52)

    SolweigCount=1

    if (OutInterval == 60) then
        OutInterval=0
    endif

    if (Posture==1) then
        Fside=0.22
        Fup=0.06
    else
        Fside=0.1666666667
        Fup=0.166666667
    endif

    !timestepdec=real(t_interval)/(3600*24)
    timestepdec=real(OutInterval)/(real(t_interval)*24.)

	!!! Loading DSM !!!
    Path=trim(FileInputPath)//trim(DSMPath)
    call LoadEsriAsciiGrid(Path,DSMName,xllcorner,yllcorner,cellsize,NoData)
    allocate(a(sizey,sizex))
    a=tempgrid
    deallocate(tempgrid)

    scale=1/cellsize

    GridFolder=trim(FileOutputPath)//'Grids'

    ! Create grid folder !! This does not work in windows. needs to be done in python
    !inquire(file=GridFolder, exist=exist)
    !if (exist) then
    !else
    !    makedirectory = 'mkdir ' // trim(GridFolder)
    !    call system(makedirectory)
    !end if

    !!! Set up for vegetation scheme, or not !!!
    if (usevegdem==1) then
        ! Calculating transmissivity of short wave radiation  through vegetation based on decid LAI
        transperLAI=(TransMax-TransMin)/(LAImax(2)-LAImin(2))
        firstday = int(MetForcingData(1,2,1))
        trans=TransMin+(LAImax(2)-LAI(firstday-1,2))*transperLAI

	    ! Loading vegDSM (SDSM)
        Path=trim(FileInputPath)//trim(DSMPath)
        call LoadEsriAsciiGrid(Path,CDSMname,xllcorner,yllcorner,cellsize,NoData)
        allocate(vegdem(sizey,sizex))
        vegdem=tempgrid
        deallocate(tempgrid)

        ! Loading trunkDSM (TDSM)
        Path=trim(FileInputPath)//trim(DSMPath)
        call LoadEsriAsciiGrid(Path,TDSMname,xllcorner,yllcorner,cellsize,NoData)
        allocate(vegdem2(sizey,sizex))
        vegdem2=tempgrid
        deallocate(tempgrid)

    	! amaxvalue (used in calculation of vegetation shadows)
        vegmax=maxval(vegdem)
        amaxvalue=maxval(a)-minval(a)
        amaxvalue=max(amaxvalue,vegmax)

    	! Elevation vegdems if buildingDSM includes ground heights
        vegdem=vegdem+a
        where (vegdem==a)
            vegdem=0.0
        end where
        vegdem2=vegdem2+a;
        where (vegdem2==a)
            vegdem2=0.0
        end where
        ! Bush separation
        allocate(bush(sizex,sizey))
        where ((vegdem>0) .and. (vegdem2==0))
            bush=vegdem
        elsewhere
            bush=0.0
        end where
    else
        trans=1.00;
    endif

    !!! Loading/creating SVFs !!!
    Path=trim(FileInputPath)//trim(SVFPath)//trim(SVFsuffix)
    svfname=(/'svf.asc ','svfE.asc','svfN.asc','svfW.asc','svfS.asc'/)
    svfvegname=(/'svfveg.asc  ','svfEveg.asc ','svfNveg.asc ','svfWveg.asc ','svfSveg.asc ',&
        'svfaveg.asc ','svfEaveg.asc','svfNaveg.asc','svfWaveg.asc','svfSaveg.asc'/)
    ! SVFs, Should be done as a loop... ! How to change variable in a loop???
    call LoadEsriAsciiGrid(Path,svfname(1),xllcorner,yllcorner,cellsize,NoData)
    allocate(svf(sizey,sizex)); svf=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(2),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfE(sizey,sizex)); svfE=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(3),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfN(sizey,sizex)); svfN=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(4),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfW(sizey,sizex)); svfW=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(5),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfS(sizey,sizex)); svfS=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(1),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfveg(sizey,sizex)); svfveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(2),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfEveg(sizey,sizex)); svfEveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(3),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfNveg(sizey,sizex)); svfNveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(4),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfWveg(sizey,sizex)); svfWveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(5),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfSveg(sizey,sizex)); svfSveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(6),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfaveg(sizey,sizex)); svfaveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(7),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfEaveg(sizey,sizex)); svfEaveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(8),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfNaveg(sizey,sizex)) ; svfNaveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(9),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfWaveg(sizey,sizex)); svfWaveg=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(10),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfSaveg(sizey,sizex)); svfSaveg=tempgrid; deallocate(tempgrid)

    !!! Loading buildings grid !!!
    Path=trim(FileInputPath)//trim(DSMPath)
    GridFile=trim(Path)//trim(buildingsname)
    inquire(file=GridFile, exist=exist)
    if (exist) then
        call LoadEsriAsciiGrid(Path,buildingsname,xllcorner,yllcorner,cellsize,NoData)
        allocate(buildings(sizey,sizex))
        buildings=tempgrid
        deallocate(tempgrid)
    else
        !!! Not ready, should return error. Also for the other grids
    endif

    ! Time related info
    !timestepdec=t_INTERVAL/(1440.*60.)
    timeadd=0.00

    ! Initiate map for surface temperature delay
    allocate(Tgmap1(sizey,sizex))
    Tgmap1=0.0

    return

274 call ErrorHint(40,trim("SOLWEIGinput.nml FileCode is missing"),notUsed,notUsed,notUsedI)

end subroutine SOLWEIG_Initial
subroutine Kside_veg_v24(radI,radG,radD,azimuth,altitude,psi,t,albedo)
use matsize

    implicit none

    real(kind(1d0)), parameter     :: pi=3.141592653589793
    real(kind(1D0)) :: vikttot,aziE,aziN,aziS,aziW
    real(kind(1D0)) :: radI,radG,radD
    real(kind(1D0)) :: azimuth,altitude,psi,t,albedo

    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: svfviktbuveg

    allocate(svfviktbuveg(sizex,sizey))
    ! New reflection equation 2012-05-25
    vikttot=4.4897
    aziE=azimuth+t
    aziS=azimuth-90+t
    aziW=azimuth-180+t
    aziN=azimuth-270+t
    ! sunw=cos(altitude*(pi/180)); ! anngle to walls
    !! Kside with weights
    call Kvikt_veg( svfE,svfEveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (360-t) .or. azimuth <= (180-t)) then
        Keast=radI*shadow*cos(altitude*(pi/180))*sin(aziE*(pi/180))+ &
        radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Keast=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if

    call Kvikt_veg( svfS,svfSveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (90-t) .and. azimuth <= (270-t)) then
    Ksouth=radI*shadow*cos(altitude*(pi/180))*sin(aziS*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Ksouth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if

    call Kvikt_veg( svfW,svfWveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (180-t) .and. azimuth <= (360-t)) then
    Kwest=radI*shadow*cos(altitude*(pi/180))*sin(aziW*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Kwest=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if

    call Kvikt_veg( svfN,svfNveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth <= (90-t) .or. azimuth > (270-t)) then
    Knorth=radI*shadow*cos(altitude*(pi/180))*sin(aziN*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Knorth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if

    deallocate(svfviktbuveg)
    end subroutine Kside_veg_v24


subroutine Kvikt_veg(isvf, isvfveg, vikttot)
use matsize

    implicit none
    real(kind(1D0)) :: vikttot
    real(kind(1d0)), dimension(sizey,sizex) :: isvf
    real(kind(1d0)), dimension(sizey,sizex) :: isvfveg
    real(kind(1d0)), dimension(sizey,sizex) :: svfvegbu

    !! Least
    viktwall=(vikttot-(63.227*isvf**6-161.51*isvf**5+156.91*isvf**4-70.424*isvf**3+16.773*isvf**2-0.4863*isvf))/vikttot

    svfvegbu=(isvfveg+isvf-1) ! Vegetation plus buildings
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2&
         &-0.4863*svfvegbu))/vikttot
    viktveg=viktveg-viktwall
end subroutine Kvikt_veg
subroutine Lside_veg_v2(altitude,Ta,Tw,SBC,ewall,esky,t,CI)!azimuth,
! This m-file is the current one that estimates L from the four cardinal points 20100414
use matsize
use data_in

    implicit none
    real(kind(1D0))                             :: vikttot,aziE,aziN,aziS,aziW
    real(kind(1D0))                             :: altitude,Ta,Tw,SBC,ewall,esky,t,CI,c!azimuth,
    real(kind(1d0)),allocatable,dimension(:,:)  :: svfalfaE,svfalfaS,svfalfaW,svfalfaN
    real(kind(1d0)),allocatable,dimension(:,:)  :: alfaB,betaB,betasun
    real(kind(1d0)),allocatable,dimension(:,:)  :: Lground,Lrefl,Lsky,Lsky_allsky,Lveg,Lwallsh,Lwallsun
    real(kind(1d0)),allocatable,dimension(:,:)  :: viktonlywall,viktaveg,svfvegbu
    real(kind(1d0)),allocatable,dimension(:,:)  :: oneminussvfE,oneminussvfS,oneminussvfW,oneminussvfN
    real(kind(1d0)),parameter                   :: pi=3.141592653589793

    allocate(oneminussvfE(sizex,sizey))
    allocate(oneminussvfS(sizex,sizey))
    allocate(oneminussvfW(sizex,sizey))
    allocate(oneminussvfN(sizex,sizey))
    allocate(svfalfaE(sizex,sizey))
    allocate(svfalfaS(sizex,sizey))
    allocate(svfalfaW(sizex,sizey))
    allocate(svfalfaN(sizex,sizey))
    allocate(alfaB(sizex,sizey))
    allocate(betaB(sizex,sizey))
    allocate(betasun(sizex,sizey))
    allocate(Lground(sizex,sizey))
    allocate(Lrefl(sizex,sizey))
    allocate(Lsky(sizex,sizey))
    allocate(Lsky_allsky(sizex,sizey))
    allocate(Lveg(sizex,sizey))
    allocate(Lwallsh(sizex,sizey))
    allocate(Lwallsun(sizex,sizey))
    allocate(viktonlywall(sizex,sizey))
    allocate(viktaveg(sizex,sizey))
    allocate(svfvegbu(sizex,sizey))

    if (allocated(viktwall)) deallocate(viktwall); allocate(viktwall(sizey,sizex))
    if (allocated(viktsky)) deallocate(viktsky); allocate(viktsky(sizey,sizex))
    if (allocated(viktveg)) deallocate(viktveg); allocate(viktveg(sizey,sizex))
    if (allocated(viktrefl)) deallocate(viktrefl); allocate(viktrefl(sizey,sizex))

    !Building height angle from svf
    oneminussvfE=1.-svfE; where (oneminussvfE<=0) oneminussvfE=0.000000001 ! avoiding log(0)
    oneminussvfS=1.-svfS; where (oneminussvfS<=0) oneminussvfS=0.000000001 ! avoiding log(0)
    oneminussvfW=1.-svfW; where (oneminussvfW<=0) oneminussvfW=0.000000001 ! avoiding log(0)
    oneminussvfN=1.-svfN; where (oneminussvfN<=0) oneminussvfN=0.000000001 ! avoiding log(0)
    svfalfaE=asin(exp((log(oneminussvfE))/2))
    svfalfaS=asin(exp((log(oneminussvfS))/2))
    svfalfaW=asin(exp((log(oneminussvfW))/2))
    svfalfaN=asin(exp((log(oneminussvfN))/2))

    vikttot=4.4897
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t

    ! F_sh=cylindric_wedge(zen,svfalfa);!Fraction shadow on building walls based on sun altitude and svf
    ! F_sh(isnan(F_sh))=0.5;
    F_sh=2.*F_sh-1. !(cylindric_wedge scaled 0-1)

     if (SOLWEIG_ldown==1) then
        c=1-CI
        Lsky_allsky=esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4)
    else
        Lsky_allsky=ldown
    end if

    !! Least
    call Lvikt_veg(svfE,svfEveg,svfEaveg,vikttot)

    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaE)
        betaB=atan(tan((svfalfaE)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (180-t) .and. azimuth <= (360-t)) then
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziE*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if
    Lsky=((svfE+svfEveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Least=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lsouth
    call Lvikt_veg(svfS,svfSveg,svfSaveg,vikttot)

    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaS)
        betaB=atan(tan((svfalfaS)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth <= (90-t) .or. azimuth > (270-t)) then
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziS*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if
    Lsky=((svfS+svfSveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lsouth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lwest
    call Lvikt_veg(svfW,svfWveg,svfWaveg,vikttot)

    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaW)
        betaB=atan(tan((svfalfaW)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (360-t) .or. azimuth <= (180-t)) then
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziW*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if
    Lsky=((svfW+svfWveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lwest=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lnorth
    call Lvikt_veg(svfN,svfNveg,svfNaveg,vikttot)

    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaN)
        betaB=atan(tan((svfalfaN)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (90-t) .and. azimuth <= (270-t)) then
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziN*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if
    Lsky=((svfN+svfNveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lnorth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    deallocate(svfalfaE)
    deallocate(svfalfaS)
    deallocate(svfalfaW)
    deallocate(svfalfaN)
    deallocate(alfaB)
    deallocate(betaB)
    deallocate(betasun)
    deallocate(Lground)
    deallocate(Lrefl)
    deallocate(Lsky)
    deallocate(Lsky_allsky)
    deallocate(Lveg)
    deallocate(Lwallsh)
    deallocate(Lwallsun)
    deallocate(viktonlywall)
    deallocate(viktaveg)
    deallocate(svfvegbu)

end subroutine Lside_veg_v2
subroutine Lvikt_veg(isvf,isvfveg,isvfaveg,vikttot)
use matsize

    implicit none
    real(kind(1D0)) :: vikttot
    real(kind(1d0)), dimension(sizey,sizex) :: isvf
    real(kind(1d0)), dimension(sizey,sizex) :: isvfveg
    real(kind(1d0)), dimension(sizey,sizex) :: isvfaveg
    real(kind(1d0)), dimension(sizey,sizex) :: viktonlywall
    real(kind(1d0)), dimension(sizey,sizex) :: viktaveg
    !real(kind(1d0)), dimension(sizey,sizex) :: viktwall
    real(kind(1d0)), dimension(sizey,sizex) :: svfvegbu

    !allocate(svfalfaE(sizex,sizey))
    !allocate(svfalfaS(sizex,sizey))
    !allocate(svfalfaW(sizex,sizey))
    !allocate(svfalfaN(sizex,sizey))
    !allocate(alfaB(sizex,sizey))
    !allocate(betaB(sizex,sizey))

    !! Least
    viktonlywall=(vikttot-(63.227*isvf**6-161.51*isvf**5+156.91*isvf**4-70.424*isvf**3+16.773*isvf**2-0.4863*isvf))/vikttot

    viktaveg=(vikttot-(63.227*isvfaveg**6-161.51*isvfaveg**5+156.91*isvfaveg**4-70.424*isvfaveg**3+16.773*isvfaveg**2&
         &-0.4863*isvfaveg))/vikttot

    viktwall=viktonlywall-viktaveg

    svfvegbu=(isvfveg+isvf-1) ! Vegetation plus buildings
    viktsky=(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu)/vikttot
    viktrefl=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863&
         &*svfvegbu))/vikttot
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863&
         &*svfvegbu))/vikttot
    viktveg=viktveg-viktwall

end subroutine Lvikt_veg
!  FUNCTION TO RETURN 0 IF IX=0, 1 IF 0<IX<MAXPOS+1,-1 OTHERWISE.
!  MAXPOS is given as the maximum positive Integer.
subroutine issign(IX,MAXPOS,ISIGNM)
      real(kind(1d0)) IX,MAXPOS,ISIGNM
      ISIGNM=1.0
      IF(IX.LT.0.OR.IX.GT.MAXPOS)ISIGNM=-1
      IF(IX.EQ.0) ISIGNM=0
      RETURN
end subroutine issign


!=====================================================
! fuction to convert interger to string
character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    end function str


subroutine SaveGrids

use matsize
use solweig_module
use data_in
use time

implicit none
character(len=100)       :: GridPath,GridName,GridPath2
character(len=4)       ::doy,hour
!real(kind(1d0)), allocatable, dimension(:,:):: savegrid

    allocate(savegrid(sizey,sizex))

    if (Tmrt_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Tmrt_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Tmrt
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Lup2d_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Tmrt_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Tmrt
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Ldown2d_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Ldown2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Ldown2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Kup2d_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Kup2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Kup2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Kdown2d_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Kdown2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Kdown2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (GVF_out==1) then
        Gridpath2='Grids/'
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id
        write(hour,'(i2)') it
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='GVF_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=gvf
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if

    deallocate(savegrid)

end subroutine SaveGrids
	!------------------------------------------------!
	! Shadow casting algorithm, ground and buildings !
    !------------------------------------------------!
    subroutine shadowingfunction_10(azimuth,altitude,scale)
    use matsize
	! This m.file calculates shadows on a DSM
    ! Code originate from Carlo Rattis thesis
    ! This code is translated from Matlab by
    ! Fredrik Lindberg, Gothenburg University
    ! Last modified LJ 27 Jan 2016 - Removal of tabs and fixing real-int conversions

    implicit none
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0)), parameter          :: maxpos=10000000000.0
    !real, allocatable, dimension(:,:)   :: a,f,temp !! Defined in matsize
    real(kind(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    real(kind(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    integer                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,row,col !,sizex,sizey
	! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: temp,f

    !special cases
    if (altitude==90) then
        altitude=altitude-0.0001
    end if
    if (azimuth==0) then
        azimuth=azimuth-0.0001
    end if

    ! conversion
    degrees=pi/180
    azi=azimuth*degrees
    alt=altitude*degrees

    allocate(f(sizex,sizey))
    allocate(temp(sizex,sizey))

    if (allocated(sh)) deallocate(sh)
    allocate(sh(sizex,sizey))

    ! initialise parameters
    f=a
    dx=0
    dy=0
    dz=0
    temp=a*0.0
    index=1

    ! other loop parameters
    amaxvalue=maxval(a)
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=sin(azi)
    cosazimuth=cos(azi)
    tanazimuth=tan(azi)
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=abs(1./sinazimuth)
    dscos=abs(1./cosazimuth)
    tanaltitudebyscale=tan(alt)/scale


    ! main loop
    DO WHILE (amaxvalue>=dz .and. abs(dx)<=sizex .and. abs(dy)<=sizey)

        IF ((pibyfour <= azi .and. azi < threetimespibyfour) .or. (fivetimespibyfour <= azi .and. azi < seventimespibyfour)) THEN
            dy=signsinazimuth*index
            dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
            ds=dssin
        ELSE
            dy=signsinazimuth*abs(nint(index*tanazimuth))
            dx=-1.*signcosazimuth*index
            ds=dscos
        END IF

        dz=ds*index*tanaltitudebyscale
        temp=temp*0

        absdx=abs(dx)
        absdy=abs(dy)

        xc1=int((dx+absdx)/2)+1
        xc2=(sizex+int((dx-absdx)/2))
        yc1=int((dy+absdy)/2)+1
        yc2=(sizey+int((dy-absdy)/2))
        xp1=-int((dx-absdx)/2)+1
        xp2=(sizex-int((dx+absdx)/2))
        yp1=-int((dy-absdy)/2)+1
        yp2=(sizey-int((dy+absdy)/2))

        temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz

        f=max(f,temp)
        index=index+1

    END DO

    f=f-a
    where (f>0)
        f=-1
    end where
    sh=f+1
    !sh=f ! invert as in shadowingfunctionglobalradiation
    deallocate(f)
    deallocate(temp)

end subroutine shadowingfunction_10
    !------------------------------------------------!
    ! Shadow casting algorithm, vegetation           !
    !------------------------------------------------!

subroutine shadowingfunction_20(azimuth,altitude,scale,amaxvalue)
use matsize
    ! This m.file calculates shadows on a DSM and for vegetation units
    ! This code is translated from Matlab by Fredrik Lindberg, Gothenburg University

    implicit none
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0)), parameter          :: maxpos=10000000000.0
    real(kind(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    real(kind(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    integer                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,test!,row,col !,sizex,sizey
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: f,temp,tmp,stopbuild,stopveg,g,bushplant,tempvegdem,tempvegdem2
    real(kind(1d0)),allocatable,dimension(:,:) :: fabovea,gabovea,tempbush,firstvegdem,vegsh2

    !real 		  				  	    :: start_time,end_time

    !special case
    if (altitude==90) then
        altitude=altitude-0.0001
    end if
    if (azimuth==0) then
        azimuth=azimuth-0.0001
    end if

    ! conversion
    degrees=pi/180;
    azi=azimuth*degrees;
    alt=altitude*degrees;

    if (allocated(sh)) deallocate(sh)
    allocate(sh(sizex,sizey))
    if (allocated(vegsh)) deallocate(vegsh)
    allocate(vegsh(sizex,sizey))
    if (allocated(vbshvegsh)) deallocate(vbshvegsh)
    allocate(vbshvegsh(sizex,sizey))

	! allocation of grids
    allocate(f(sizex,sizey))
    allocate(temp(sizex,sizey))
    allocate(tmp(sizex,sizey))
    allocate(stopbuild(sizex,sizey))
    allocate(stopveg(sizex,sizey))
    allocate(g(sizex,sizey))
    allocate(bushplant(sizex,sizey))
    allocate(tempvegdem(sizex,sizey))
    allocate(tempvegdem2(sizex,sizey))
    allocate(fabovea(sizex,sizey))
    allocate(gabovea(sizex,sizey))
    allocate(firstvegdem(sizex,sizey))
    allocate(tempbush(sizex,sizey))
    allocate(vegsh2(sizex,sizey))

    ! initialise parameters
    f=a
    dx=0
    dy=0
    dz=0
    temp=a*0.0
    sh=temp
    vegsh=sh
    stopbuild=sh
    stopveg=sh
    vbshvegsh=sh
    g=sh
    bushplant=temp
    where (bush>1)
        bushplant=1
    end where

    index=1
	!test=0
    ! other loop parameters
    !amaxvalue=maxval(a)
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour;
    fivetimespibyfour=5.*pibyfour;
    seventimespibyfour=7.*pibyfour;
    sinazimuth=sin(azi);
    cosazimuth=cos(azi);
    tanazimuth=tan(azi);
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth);
    !signcosazimuth=cosazimuth/abs(cosazimuth);
    dssin=abs(1./sinazimuth);
    dscos=abs(1./cosazimuth);
    tanaltitudebyscale=tan(alt)/scale;

    DO WHILE (amaxvalue>=dz .and. abs(dx)<=sizex .and. abs(dy)<=sizey)

        IF ((pibyfour <= azi .and. azi < threetimespibyfour) .or. (fivetimespibyfour <= azi .and. azi < seventimespibyfour)) THEN
            dy=signsinazimuth*index
            dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
            ds=dssin
        ELSE
            dy=signsinazimuth*abs(nint(index*tanazimuth))
            dx=-1.*signcosazimuth*index
            ds=dscos
        END IF

        dz=ds*index*tanaltitudebyscale
        temp=temp*0
        tempvegdem=temp
        tempvegdem2=temp

        absdx=abs(dx)
        absdy=abs(dy)

        xc1=int(((dx+absdx)/2))+1
        xc2=(sizex+int((dx-absdx)/2))
        yc1=int((dy+absdy)/2)+1
        yc2=(sizey+int((dy-absdy)/2))
        xp1=-int((dx-absdx)/2)+1
        xp2=(sizex-int((dx+absdx)/2))
        yp1=-int((dy-absdy)/2)+1
        yp2=(sizey-int((dy+absdy)/2))

        temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz
        tempvegdem(xp1:xp2,yp1:yp2)=vegdem(xc1:xc2,yc1:yc2)-dz
        tempvegdem2(xp1:xp2,yp1:yp2)=vegdem2(xc1:xc2,yc1:yc2)-dz

        f=max(f,temp)
        where (f>a) !sh(f>a)=1;sh(f<=a)=0; !Moving building shadow
            sh=1
        elsewhere
            sh=0
        end where
        where (tempvegdem>a) !fabovea=tempvegdem>a; !vegdem above DEM
            fabovea=1
        elsewhere
            fabovea=0
        end where
        where (tempvegdem2>a) !gabovea=tempvegdem2>a; !vegdem2 above DEM
            gabovea=1
        elsewhere
            gabovea=0
        end where
        vegsh2=fabovea-gabovea
        vegsh=max(vegsh,vegsh2)
        where ((vegsh*sh)>0) !vegsh(vegsh.*sh>0)=0;! removing shadows 'behind' buildings
           vegsh=0
        end where
        vbshvegsh=vegsh+vbshvegsh

    	! vegsh at high sun altitudes
        if (index==1) then
            firstvegdem=tempvegdem-temp
            where (firstvegdem<=0)!firstvegdem(firstvegdem<=0)=1000;
                firstvegdem=1000
            end where
            where (firstvegdem<dz)!vegsh(firstvegdem<dz)=1;
                vegsh=1
            end where
            tmp=temp*0.0
            where (vegdem2>a)!vegsh=vegsh.*(vegdem2>a);
                tmp=1
            end where
            vegsh=vegsh*tmp
            vbshvegsh=temp*0.0 !vbshvegsh=zeros(sizex,sizey);
        end if

    	! Bush shadow on bush plant
        tmp=fabovea*bush
        if ((maxval(bush)>0) .and. (maxval(tmp)>0)) then
            tempbush=temp*0.0
            tempbush(xp1:xp2,yp1:yp2)=bush(xc1:xc2,yc1:yc2)-dz
            g=max(g,tempbush)
            g=bushplant*g
        end if

        index=index+1
    END DO

    sh=1-sh
    where (vbshvegsh>0)!vbshvegsh(vbshvegsh>0)=1;
        vbshvegsh=1
    end where
    vbshvegsh=vbshvegsh-vegsh;

    if (maxval(bush)>0) then
        g=g-bush
        where (g>0)!g(g>0)=1;g(g<0)=0;
            g=1
        elsewhere
            g=0
        end where
        vegsh=vegsh-bushplant+g
        where (vegsh<0)!vegsh(vegsh<0)=0;
           vegsh=0
        end where
    end if

    where (vegsh>0)!vegsh(vegsh>0)=1;
        vegsh=1
    end where
    vegsh=1-vegsh
    vbshvegsh=1-vbshvegsh

    !deallocation of grids
    deallocate(f)
    deallocate(temp)
    deallocate(tmp)
    deallocate(stopbuild)
    deallocate(stopveg)
    deallocate(g)
    deallocate(bushplant)
    deallocate(tempvegdem)
    deallocate(tempvegdem2)
    deallocate(fabovea)
    deallocate(gabovea)
    deallocate(firstvegdem)
    deallocate(tempbush)
    deallocate(vegsh2)

end subroutine shadowingfunction_20
 subroutine sunonsurface_veg(iazimuthA, scale, first, second, psi)
 ! This m-file creates a boolean image of sunlit walls.
 ! Shadows from both buildings and vegetation is accounted for
 ! moving building in the direction of the sun
 ! Last modified:
 ! LJ 27 Jan 2016 - Removal of tabs and fixing real-int conversions

 use matsize

    implicit none
    real(kind(1d0))             :: iazimuthA,iazimuth,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))             :: scale
    integer                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2,n,first,second
    real(kind(1d0))             :: dx,dy,ds,absdx,absdy,psi !,dz
    real(kind(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))             :: seventimespibyfour
    real(kind(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos
    real(kind(1d0)),allocatable,dimension(:,:)  ::weightsumwall,weightsumsh,gvf1,gvf2
    real(kind(1d0)),allocatable,dimension(:,:)  ::f,tempsh,tempbu,tempbub,tempwallsun,tempb,sh1

    real(kind(1d0)), parameter  :: pi=3.141592653589793
    real(kind(1d0)), parameter  :: maxpos=10000000000.0

    allocate(weightsumwall(sizex,sizey))
    allocate(weightsumsh(sizex,sizey))
    allocate(gvf1(sizex,sizey))
    allocate(gvf2(sizex,sizey))
    allocate(f(sizex,sizey))
    allocate(tempsh(sizex,sizey))
    allocate(tempbu(sizex,sizey))
    allocate(tempbub(sizex,sizey))
    allocate(tempwallsun(sizex,sizey))
    allocate(tempb(sizex,sizey))
    allocate(sh1(sizex,sizey))

    iazimuth=iazimuthA*(pi/180)
    !special cases
    if (iazimuth==0) then
        iazimuth=iazimuth+0.000001
    end if
    ! loop parameters
    index=0
    f=buildings
    sh1=sh-(1-vegsh)*(1-psi)
    dx=0
    dy=0
    ds=0

    tempsh = 0.0D0
    tempbu = 0.0D0
    tempbub = 0.0D0
    tempwallsun = 0.0D0
    !sh = 0.0D0
    weightsumsh = 0.0D0
    weightsumwall = 0.0D0

    first=int(real(first,kind(1d0))*scale)  !Int added around the equation as first and second are both integers
    second=int(real(second,kind(1d0))*scale)

    ! other loop parameters
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=sin(iazimuth)
    cosazimuth=cos(iazimuth)
    tanazimuth=tan(iazimuth)
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=abs(1./sinazimuth)
    dscos=abs(1./cosazimuth)

    !! The Shadow casting algorithm
    do n=1,second
        IF ((pibyfour <= iazimuth .and. iazimuth < threetimespibyfour) .or. (fivetimespibyfour&
             & <= iazimuth .and. iazimuth < seventimespibyfour)) THEN
            dy=signsinazimuth*index
            dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
            ds=dssin
        ELSE
            dy=signsinazimuth*abs(nint(index*tanazimuth))
            dx=-1.*signcosazimuth*index
            ds=dscos
        END IF

        absdx=abs(dx)
        absdy=abs(dy)

        xc1=int((dx+absdx)/2)+1
        xc2=(sizex+int((dx-absdx)/2))
        yc1=int((dy+absdy)/2)+1
        yc2=(sizey+int((dy-absdy)/2))
        xp1=-int((dx-absdx)/2)+1
        xp2=(sizex-int((dx+absdx)/2))
        yp1=-int((dy-absdy)/2)+1
        yp2=(sizey-int((dy+absdy)/2))

        tempbu(xp1:xp2,yp1:yp2)=buildings(xc1:xc2,yc1:yc2) !moving building

        tempsh(xp1:xp2,yp1:yp2)=sh1(xc1:xc2,yc1:yc2) !moving shadow image
        f=min(f,tempbu) !utsmetning of buildings

        weightsumsh=weightsumsh+tempsh*f

        tempwallsun(xp1:xp2,yp1:yp2)=sunwall(xc1:xc2,yc1:yc2) !moving building wall in sun image
        tempb=tempwallsun*f
        where ((tempb+tempbub)>0) !tempbub=(tempb+tempbub)>0==1
            tempbub=1.
        end where

        weightsumwall=weightsumwall+tempbub

        if (index*scale==first) then
            gvf1=(weightsumwall+weightsumsh)/first
            where (gvf1>1)
                gvf1=1.
            end where
        end if
        index=index+1

    end do
    gvf2=(weightsumsh+weightsumwall)/second
    where (gvf2>1)
        gvf2=1.
    end where

    ! Weighting
    sos=(gvf1*0.5+gvf2*0.4)/0.9

    deallocate(weightsumwall)
    deallocate(weightsumsh)
    deallocate(gvf1)
    deallocate(gvf2)
    deallocate(f)
    deallocate(tempsh)
    deallocate(tempbu)
    deallocate(tempbub)
    deallocate(tempwallsun)
    deallocate(tempb)
    deallocate(sh1)

end subroutine sunonsurface_veg

 !>
 subroutine wallinsun_veg(azimuth)
 ! This m-file creates a boolean image of sunlit walls.
 ! Shadows from both buildings and vegetation is accounted for
 ! moving building in the direction of the sun
 ! Last modified:
 !  LJ 27 Jan 2017 - Change of equations xc1...yp2 to account for the change from real to integer
 !---------------------------------------------------------------------------------

 use matsize
    implicit none
    real(kind(1d0))             :: azimuth,iazimuth
    integer                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2
    real(kind(1d0))             :: dx,dy,dz,ds,absdx,absdy
    real(kind(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))             :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos,azi
    real(kind(1d0)), parameter  :: pi=3.141592653589793
    real(kind(1d0)), parameter  :: maxpos=10000000000.0
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: temp,sh1

    allocate(temp(sizex,sizey))
    allocate(sh1(sizex,sizey))
    !allocate(vegsh(sizex,sizey))
    !allocate(a(sizex,sizey))
    !allocate(buildings(sizex,sizey))

    if (allocated(sunwall)) deallocate(sunwall); allocate(sunwall(sizey,sizex))

    iazimuth=azimuth+180
    if (iazimuth>=360) then
        iazimuth=iazimuth-360
    end if
    !special cases
    if (iazimuth==0) then
        iazimuth=iazimuth+0.00001
    end if
    ! conversion into radians
    azi=iazimuth*(pi/180)

    index=1
    dx=0.
    dy=0.
    dz=0.
    ds=0.

    temp = 0.0D0

    ! other loop parameters
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=sin(azi)
    cosazimuth=cos(azi)
    tanazimuth=tan(azi)
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=abs(1./sinazimuth)
    dscos=abs(1./cosazimuth)

    sh1=vegsh+sh-1.
    !! The Shadow casting algoritm
    IF ((pibyfour <= azi .and. azi < threetimespibyfour) .or. (fivetimespibyfour <= azi .and. azi < seventimespibyfour)) THEN
        dy=signsinazimuth*index
        dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
        ds=dssin
    ELSE
        dy=signsinazimuth*abs(nint(index*tanazimuth))
        dx=-1.*signcosazimuth*index
        ds=dscos
    END IF

    ! note: dx and dy represent absolute values while ds is an incremental value

    absdx=abs(dx)
    absdy=abs(dy)

    xc1=int((dx+absdx)/2)+1  !LJ added int to the equation to account for the conversion from real to int
    xc2=(sizex+int((dx-absdx)/2))
    yc1=int((dy+absdy)/2)+1
    yc2=(sizey+int((dy-absdy)/2))
    xp1=-int((dx-absdx)/2)+1
    xp2=(sizex-int((dx+absdx)/2))
    yp1=-int((dy-absdy)/2)+1
    yp2=(sizey-int((dy+absdy)/2))

    temp(xp1:xp2,yp1:yp2)= buildings(xc1:xc2,yc1:yc2)

    sunwall=temp-buildings
    where (sunwall==1) !f1(f1==1)=0
          sunwall=0
    end where
    where (sunwall==-1) !f1(f1==-1)=1
          sunwall=1
    end where
    sunwall=sh1*sunwall

    deallocate(temp)
    deallocate(sh1)

end subroutine wallinsun_veg


end module solweig_module
