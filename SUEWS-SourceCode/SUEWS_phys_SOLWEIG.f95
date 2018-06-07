!Last modified: LJ 27 Jan 2016 - Removal of tabs
! module matsize
!
!    IMPLICIT NONE
!
!    integer    :: sizex,sizey ! number of rows and cols of grid
!    real(kind(1d0)), allocatable, dimension(:,:)  :: a,sh,vbshvegsh,vegsh
!    real(kind(1d0)), allocatable, dimension(:,:)  :: bush,vegdem,vegdem2,tempgrid
!    real(kind(1d0)), allocatable, dimension(:,:)  :: buildings,svf,svfE,svfS,svfW,svfN
!    real(kind(1d0)), allocatable, dimension(:,:)  :: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
!    real(kind(1d0)), allocatable, dimension(:,:)  :: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
!    real(kind(1d0)), allocatable, dimension(:,:)  :: Kdown2d,Keast,Knorth,Ksouth,Kup2d,Kwest
!    real(kind(1d0)), allocatable, dimension(:,:)  :: Ldown2d,Least,Lnorth,Lsouth,Lup2d,Lwest
!    real(kind(1d0)), allocatable, dimension(:,:)  :: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
!    real(kind(1d0)), allocatable, dimension(:,:)  :: svfalfa,sos,Tgmap1
!    real(kind(1d0)), allocatable, dimension(:,:)  :: viktveg,viktsky,viktrefl,viktwall,savegrid
!
! end module matsize

!svfvegbu,,viktaveg,viktonlywall,svfalfaE,svfalfaN,svfalfaS,svfalfaW,alfaB,betaB,betasun,Lground,Lrefl,Lsky,Lveg,Lwallsh,Lwallsun
! alfa,xa,ha,hkil,Ai,phi,qa,Za,tempsh,tempbu,tempbub,tempwallsun,weightsumwall,weightsumsh,gvf1,gvf2
!!,svfbuveg
!Lnight,,tempvegdem,tempvegdem2,fabovea,gabovea,tempbush,firstvegdem,vegsh2,tempgrid
!svfviktbuveg
!,stopbuild,stopveg,g,bushplant,Tgmap
!,temp,f,tmp,Knight,,tempb
!Tgmapgvf

MODULE solweig_module
  USE data_in
  USE gis_data
  USE time
  USE allocateArray
  USE InitialCond
  USE sues_data
  USE defaultNotUsed

  IMPLICIT NONE

  REAL(KIND(1d0)) :: timestepdec,& !time step in decimal time
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
  INTEGER         :: SolweigCount

  INTEGER    :: sizex,sizey ! number of rows and cols of grid
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: a,sh,vbshvegsh,vegsh
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: bush,vegdem,vegdem2,tempgrid
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: buildings,svf,svfE,svfS,svfW,svfN
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: Kdown2d,Keast,Knorth,Ksouth,Kup2d,Kwest
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: Ldown2d,Least,Lnorth,Lsouth,Lup2d,Lwest
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: svfalfa,sos,Tgmap1
  REAL(KIND(1d0)), ALLOCATABLE, DIMENSION(:,:)  :: viktveg,viktsky,viktrefl,viktwall,savegrid
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

    ! USE matsize
    ! USE solweig_module


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
  SUBROUTINE clearnessindex_2013b(zen,jday,Ta,RH,radG,lat,P,I0,CI,Kt,I0et,CIuncorr)
    !Last modified:
    !LJ 27 Jan 2016 - Removal of tabs

    IMPLICIT NONE
    ! Use somemodule

    REAL(KIND(1d0))                 :: CI,CIuncorr,I0,I0et,Kt
    INTEGER                         :: jday
    REAL(KIND(1d0))                 :: P,radG,RH,Ta,zen,lat,iG,Itoa
    REAL(KIND(1d0)), DIMENSION(4)   :: G
    REAL(KIND(1d0)), PARAMETER          :: pi=3.141592653589793

    ! Variable declarations
    REAL*8 :: a2      !>
    !logical :: b      !>
    REAL*8 :: b2      !>
    REAL*8 :: corr      !>
    REAL*8 :: D      !>
    REAL*8 :: m      !>
    REAL*8 :: Tar      !>
    REAL*8 :: Td      !>
    REAL*8 :: Trpg      !>
    REAL*8 :: Tw      !>
    REAL*8 :: u      !>
    !

    ! Clearness Index at the Earth's surface calculated from Crawford and Duchon 1999

    IF (P==-999) THEN
       p=1013 !Pressure in millibars
    ELSE
       p=P*10 !Convert from hPa to millibars
    END IF
    Itoa=1370 !Effective solar constant
    !call solar_ESdist(jday,D)
    CALL sun_distance(jday,D)
    !D=sun_distance(jday) !irradiance differences due to Sun-Earth distances
    m=35.*COS(zen)*((1224.*(COS(zen)**2)+1.)**(-1./2.)) !optical air mass at p=1013
    Trpg=1.021-0.084*(m*(0.000949*p+0.051))**0.5 !Transmission coefficient for Rayliegh scattering and permanent gases

    ! empirical constant depending on latitude
    IF (lat<10) THEN
       G=(/3.37,2.85,2.80,2.64/)
    ELSE IF (lat>=10 .AND. lat<20) THEN
       G=(/2.99,3.02,2.70,2.93/)
    ELSE IF (lat>=20 .AND. lat<30) THEN
       G=(/3.60,3.00,2.98,2.93/)
    ELSE IF (lat>=30 .AND. lat<40) THEN
       G=(/3.04,3.11,2.92,2.94/)
    ELSE IF (lat>=40 .AND. lat<50) THEN
       G=(/2.70,2.95,2.77,2.71/)
    ELSE IF (lat>=50 .AND. lat<60) THEN
       G=(/2.52,3.07,2.67,2.93/)
    ELSE IF (lat>=60 .AND. lat<70) THEN
       G=(/1.76,2.69,2.61,2.61/)
    ELSE IF (lat>=70 .AND. lat<80) THEN
       G=(/1.60,1.67,2.24,2.63/)
    ELSE IF (lat>=80 .AND. lat<90) THEN
       G=(/1.11,1.44,1.94,2.02/)
    END IF
    IF (jday > 335 .OR. jday <= 60) THEN
       iG=G(1)
    ELSE IF (jday > 60 .AND. jday <= 152) THEN
       iG=G(2)
    ELSE IF (jday > 152 .AND. jday <= 244) THEN
       iG=G(3)
    ELSE IF (jday > 244 .AND. jday <= 335) THEN
       iG=G(4)
    END IF
    !dewpoint calculation
    a2=17.27
    b2=237.7
    Td=(b2*(((a2*Ta)/(b2+Ta))+LOG(RH)))/(a2-(((a2*Ta)/(b2+Ta))+LOG(RH)))
    Td=(Td*1.8)+32 !Dewpoint (�F)
    u=EXP(0.1133-LOG(iG+1)+0.0393*Td) !Precipitable water
    Tw=1-0.077*((u*m)**0.3) !Transmission coefficient for water vapor
    Tar=0.935**m !Transmission coefficient for aerosols

    I0=Itoa*COS(zen)*Trpg*Tw*D*Tar

!!! This needs to be checked !!!
    !b=I0==abs(zen)>pi/2
    !I0(b==1)=0
    !clear b
    !if (not(isreal(I0))) then
    !    I0=0
    !end if

    corr=0.1473*LOG(90-(zen/pi*180))+0.3454 ! 20070329

    CIuncorr=radG/I0
    CI=CIuncorr+(1-corr)
    I0et=Itoa*COS(zen)*D !extra terrestial solar radiation
    Kt=radG/I0et

!!!! This needs to be checked !!!
    !if (isnan(CI)) then
    !    CI=NaN
    !end if
  END SUBROUTINE clearnessindex_2013b

  !===============================================================================

  SUBROUTINE sun_distance(jday,D)

    ! Calculates solar irradiance variation based on mean earth sun distance
    ! with day of year as input.
    ! Partridge and Platt, 1975

    INTEGER          ::jday
    REAL(KIND(1d0))  ::b,D

    b = 2*3.141592654*jday/365
    D = SQRT(1.00011 + 0.034221 * COS(b) + 0.001280 * SIN(b) + 0.000719 * COS(2*b) + 0.000077 * SIN(2*b))

  END SUBROUTINE sun_distance

  SUBROUTINE cylindric_wedge(zen)
    ! Fraction of sunlit walls based on sun altitude and svf wieghted building angles

    IMPLICIT NONE

    REAL(KIND(1d0)), PARAMETER          :: pi=3.141592653589793
    REAL(KIND(1d0))         :: zen      !>
    REAL(KIND(1d0))         :: beta      !>
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: alfa,xa,ha,hkil,ba
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: Ai,phi,qa,Za
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: ukil,Ssurf
    !real(kind(1d0)), dimension(sizey,sizex) ::
    !real(kind(1d0)), dimension(sizey,sizex) ::

    ALLOCATE(alfa(sizey,sizex))
    ALLOCATE(ba(sizey,sizex))
    ALLOCATE(ha(sizey,sizex))
    ALLOCATE(xa(sizey,sizex))
    ALLOCATE(qa(sizey,sizex))
    ALLOCATE(Za(sizey,sizex))
    ALLOCATE(phi(sizey,sizex))
    ALLOCATE(ukil(sizey,sizex))
    ALLOCATE(Ai(sizey,sizex))
    ALLOCATE(Ssurf(sizey,sizex))
    ALLOCATE(hkil(sizey,sizex))

    beta=zen
    alfa=svfalfa

    xa=1.-2./(TAN(alfa)*TAN(beta))
    ha=2./(TAN(alfa)*TAN(beta))
    ba=(1./TAN(alfa))
    hkil=2.*ba*ha


    qa = 0.0D0

    WHERE (xa<0) !qa(xa<0)=tan(beta)/2
       qa=TAN(beta)/2
    END WHERE


    Za = 0.0D0

    phi = 0.0D0

    Ai = 0.0D0

    ukil = 0.0D0
    WHERE (xa<0)
       !Za(xa<0)=((ba(xa<0).**2)-((qa(xa<0).**2)./4)).**0.5
       Za=(ba**2-qa**2/4.)**0.5
       !phi(xa<0)=atan(Za(xa<0)./qa(xa<0))
       phi=ATAN(Za/qa)
       !A1(xa<0)=(sin(phi(xa<0))-phi(xa<0).*cos(phi(xa<0)))./(1-cos(phi(xa<0)))
       Ai=(SIN(phi)-phi*COS(phi))/(1-COS(phi))
       !ukil(xa<0)=2*ba(xa<0).*xa(xa<0).*A1(xa<0)
       ukil=2*ba*xa*Ai
    END WHERE

    Ssurf=hkil+ukil

    F_sh=(2*pi*ba-Ssurf)/(2*pi*ba) !Xa

    DEALLOCATE(alfa)
    DEALLOCATE(ba)
    DEALLOCATE(ha)
    DEALLOCATE(xa)
    DEALLOCATE(qa)
    DEALLOCATE(Za)
    DEALLOCATE(phi)
    DEALLOCATE(ukil)
    DEALLOCATE(Ai)
    DEALLOCATE(Ssurf)
    DEALLOCATE(hkil)

  END SUBROUTINE cylindric_wedge

  ! This subroutine estimates diffuse and directbeam radiation according to
  ! Reindl et al (1990), Solar Energy 45:1
  SUBROUTINE diffusefraction(radG,altitude,Kt,Ta,RH,radI,radD)
    IMPLICIT NONE

    REAL(KIND(1d0))                 :: radG,altitude,Kt,Ta,RH,radD,radI,alfa
    REAL(KIND(1D0)),PARAMETER       :: DEG2RAD=0.017453292,RAD2DEG=57.29577951  !!Already defined in AllocateArray module. Delete??

    alfa=altitude*DEG2RAD

    IF (Ta<=-99 .OR. RH<=-99) THEN !.or. isnan(Ta) .or. isnan(RH)) then
       IF (Kt<=0.3) THEN
          radD=radG*(1.020-0.248*Kt)
       ELSE IF (Kt>0.3 .AND. Kt<0.78) THEN
          radD=radG*(1.45-1.67*Kt)
       ELSE IF (Kt>=0.78) THEN
          radD=radG*0.147
       END IF
    ELSE
       !RH=RH/100
       IF (Kt<=0.3) THEN
          radD=radG*(1-0.232*Kt+0.0239*SIN(alfa)-0.000682*Ta+0.0195*(RH/100))
       ELSE IF (Kt>0.3 .AND. Kt<0.78) THEN
          radD=radG*(1.329-1.716*Kt+0.267*SIN(alfa)-0.00357*Ta+0.106*(RH/100))
       ELSE IF (Kt>=0.78) THEN
          radD=radG*(0.426*Kt-0.256*SIN(alfa)+0.00349*Ta+0.0734*(RH/100))
       END IF
    END IF
    radI=(radG-radD)/(SIN(alfa))

    !! Corrections for low sun altitudes (20130307)
    IF (radI<0) THEN
       radI=0
    END IF

    IF (altitude<1 .AND. radI>radG) THEN
       radI=radG
    END IF

    IF (radD>radG) THEN
       radD=radG
    END IF
  END SUBROUTINE diffusefraction
  ! This subroutine loads a ESRIACSII grid as a 2D array
  ! Last modified:
  ! LJ 27 Jan 2016-removal of tabs
  !-----------------------------------------------------
  SUBROUTINE LoadEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)

    IMPLICIT NONE
    REAL(KIND(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    INTEGER                           :: col,row
    CHARACTER(len=100)                :: GridPath,GridName,GridFile,n
    !real(kind(1d0)),intent(out)       :: tempgrid
    !real(kind(1d0)),dimension(sizey,sizex),intent(out)       :: tempgrid!,allocatable

    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=TRIM(GridPath)//TRIM(GridName)
    OPEN(99,File=GridFile,status='old')

    ! Read Header
    READ(99,*) n,sizex
    READ(99,*) n,sizey
    READ(99,*) n,xllcornerlo
    READ(99,*) n,yllcornerlo
    READ(99,*) n,cellsizelo
    READ(99,*) n,NoDatalo

    ALLOCATE(tempgrid(sizex,sizey))

    ! Read Matrix
    DO row=1,sizey
       READ(99,*) (tempgrid(row,col),col=1,sizex)
    END DO
    CLOSE(99)

    RETURN
  END SUBROUTINE LoadEsriAsciiGrid


  ! This subroutine saves a 2D array as an ESRIACSII grid
  SUBROUTINE SaveEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)


    IMPLICIT NONE
    REAL(KIND(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    INTEGER                           :: col,row
    CHARACTER(len=100)                :: GridPath,GridName,GridFile
    !integer                           :: sizey,sizex!,intent(in)
    !real(kind(1d0)), allocatable, dimension(:,:):: grid
    ! real(kind(1d0)),dimension(sizey,sizex)       :: grid!,allocatable

    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=TRIM(GridPath)//TRIM(GridName)
    OPEN(94,File=GridFile,status='unknown')

    ! Read Header
    WRITE(94,"(A5,1x,I0)") 'ncols',sizex
    WRITE(94,"(A5,1x,I0)") 'nrows',sizey
    WRITE(94,"(A9,1x,F0.2)") 'xllcorner',xllcornerlo
    WRITE(94,"(A9,1x,F0.2)") 'yllcorner',yllcornerlo
    WRITE(94,"(A8,1x,F0.2)") 'cellsize',cellsizelo
    WRITE(94,"(A12,1x,F0.2)") 'NODATA_value',NoDatalo

    ! write Matrix
    DO row=1,sizey
       WRITE(94,100) (savegrid(row,col),col=1,sizex)
    END DO
    CLOSE(94)
100 FORMAT(200(f6.2,1x))

    RETURN
  END SUBROUTINE SaveEsriAsciiGrid

  ! setup for SOLWEIG
  ! FL may 2014
  SUBROUTINE SOLWEIG_Initial
    !Last modified LJ 27 Jan 2016 - Removal of Tabs and real-int fixes

    ! USE matsize         ! All allocatable grids and related variables used in SOLWEIG


    IMPLICIT NONE

    CHARACTER(len=100)  :: Path,GridFile,GridFolder
    REAL(KIND(1d0))                 :: vegmax
    CHARACTER(len=100),DIMENSION(5) :: svfname
    CHARACTER(len=100),DIMENSION(10):: svfvegname
    LOGICAL                         :: exist
    INTEGER                         :: firstday

    NAMELIST/SOLWEIGinput/Posture,&    ! 1.Standing, 2.Sitting
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
    OPEN(52,file=TRIM(FileInputPath)//'SOLWEIGinput.nml',err=274,status='old')
    READ(52,nml=SOLWEIGinput)
    CLOSE(52)

    SolweigCount=1

    IF (OutInterval == 60) THEN
       OutInterval=0
    ENDIF

    IF (Posture==1) THEN
       Fside=0.22
       Fup=0.06
    ELSE
       Fside=0.1666666667
       Fup=0.166666667
    ENDIF

    !timestepdec=real(t_interval)/(3600*24)
    timestepdec=REAL(OutInterval)/(REAL(t_interval)*24.)

!!! Loading DSM !!!
    Path=TRIM(FileInputPath)//TRIM(DSMPath)
    CALL LoadEsriAsciiGrid(Path,DSMName,xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(a(sizey,sizex))
    a=tempgrid
    DEALLOCATE(tempgrid)

    scale=1/cellsize

    GridFolder=TRIM(FileOutputPath)//'Grids'

    ! Create grid folder !! This does not work in windows. needs to be done in python
    !inquire(file=GridFolder, exist=exist)
    !if (exist) then
    !else
    !    makedirectory = 'mkdir ' // trim(GridFolder)
    !    call system(makedirectory)
    !end if

!!! Set up for vegetation scheme, or not !!!
    IF (usevegdem==1) THEN
       ! Calculating transmissivity of short wave radiation  through vegetation based on decid LAI
       transperLAI=(TransMax-TransMin)/(LAImax(2)-LAImin(2))
       firstday = INT(MetForcingData(1,2,1))
       trans=TransMin+(LAImax(2)-LAI(firstday-1,2))*transperLAI

       ! Loading vegDSM (SDSM)
       Path=TRIM(FileInputPath)//TRIM(DSMPath)
       CALL LoadEsriAsciiGrid(Path,CDSMname,xllcorner,yllcorner,cellsize,NoData)
       ALLOCATE(vegdem(sizey,sizex))
       vegdem=tempgrid
       DEALLOCATE(tempgrid)

       ! Loading trunkDSM (TDSM)
       Path=TRIM(FileInputPath)//TRIM(DSMPath)
       CALL LoadEsriAsciiGrid(Path,TDSMname,xllcorner,yllcorner,cellsize,NoData)
       ALLOCATE(vegdem2(sizey,sizex))
       vegdem2=tempgrid
       DEALLOCATE(tempgrid)

       ! amaxvalue (used in calculation of vegetation shadows)
       vegmax=MAXVAL(vegdem)
       amaxvalue=MAXVAL(a)-MINVAL(a)
       amaxvalue=MAX(amaxvalue,vegmax)

       ! Elevation vegdems if buildingDSM includes ground heights
       vegdem=vegdem+a
       WHERE (vegdem==a)
          vegdem=0.0
       END WHERE
       vegdem2=vegdem2+a;
       WHERE (vegdem2==a)
          vegdem2=0.0
       END WHERE
       ! Bush separation
       ALLOCATE(bush(sizex,sizey))
       WHERE ((vegdem>0) .AND. (vegdem2==0))
          bush=vegdem
       ELSEWHERE
          bush=0.0
       END WHERE
    ELSE
       trans=1.00;
    ENDIF

!!! Loading/creating SVFs !!!
    Path=TRIM(FileInputPath)//TRIM(SVFPath)//TRIM(SVFsuffix)
    svfname=(/'svf.asc ','svfE.asc','svfN.asc','svfW.asc','svfS.asc'/)
    svfvegname=(/'svfveg.asc  ','svfEveg.asc ','svfNveg.asc ','svfWveg.asc ','svfSveg.asc ',&
         'svfaveg.asc ','svfEaveg.asc','svfNaveg.asc','svfWaveg.asc','svfSaveg.asc'/)
    ! SVFs, Should be done as a loop... ! How to change variable in a loop???
    CALL LoadEsriAsciiGrid(Path,svfname(1),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svf(sizey,sizex)); svf=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfname(2),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfE(sizey,sizex)); svfE=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfname(3),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfN(sizey,sizex)); svfN=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfname(4),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfW(sizey,sizex)); svfW=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfname(5),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfS(sizey,sizex)); svfS=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(1),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfveg(sizey,sizex)); svfveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(2),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfEveg(sizey,sizex)); svfEveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(3),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfNveg(sizey,sizex)); svfNveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(4),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfWveg(sizey,sizex)); svfWveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(5),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfSveg(sizey,sizex)); svfSveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(6),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfaveg(sizey,sizex)); svfaveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(7),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfEaveg(sizey,sizex)); svfEaveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(8),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfNaveg(sizey,sizex)) ; svfNaveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(9),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfWaveg(sizey,sizex)); svfWaveg=tempgrid; DEALLOCATE(tempgrid)
    CALL LoadEsriAsciiGrid(Path,svfvegname(10),xllcorner,yllcorner,cellsize,NoData)
    ALLOCATE(svfSaveg(sizey,sizex)); svfSaveg=tempgrid; DEALLOCATE(tempgrid)

!!! Loading buildings grid !!!
    Path=TRIM(FileInputPath)//TRIM(DSMPath)
    GridFile=TRIM(Path)//TRIM(buildingsname)
    INQUIRE(file=GridFile, exist=exist)
    IF (exist) THEN
       CALL LoadEsriAsciiGrid(Path,buildingsname,xllcorner,yllcorner,cellsize,NoData)
       ALLOCATE(buildings(sizey,sizex))
       buildings=tempgrid
       DEALLOCATE(tempgrid)
    ELSE
!!! Not ready, should return error. Also for the other grids
    ENDIF

    ! Time related info
    !timestepdec=t_INTERVAL/(1440.*60.)
    timeadd=0.00

    ! Initiate map for surface temperature delay
    ALLOCATE(Tgmap1(sizey,sizex))
    Tgmap1=0.0

    RETURN

274 CALL ErrorHint(40,TRIM("SOLWEIGinput.nml FileCode is missing"),notUsed,notUsed,notUsedI)

  END SUBROUTINE SOLWEIG_Initial

  SUBROUTINE Kside_veg_v24(radI,radG,radD,azimuth,altitude,psi,t,albedo)

    IMPLICIT NONE

    REAL(KIND(1d0)), PARAMETER     :: pi=3.141592653589793
    REAL(KIND(1D0)) :: vikttot,aziE,aziN,aziS,aziW
    REAL(KIND(1D0)) :: radI,radG,radD
    REAL(KIND(1D0)) :: azimuth,altitude,psi,t,albedo

    ! Internal grids
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: svfviktbuveg

    ALLOCATE(svfviktbuveg(sizex,sizey))
    ! New reflection equation 2012-05-25
    vikttot=4.4897
    aziE=azimuth+t
    aziS=azimuth-90+t
    aziW=azimuth-180+t
    aziN=azimuth-270+t
    ! sunw=cos(altitude*(pi/180)); ! anngle to walls
    !! Kside with weights
    CALL Kvikt_veg( svfE,svfEveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    IF (azimuth > (360-t) .OR. azimuth <= (180-t)) THEN
       Keast=radI*shadow*COS(altitude*(pi/180))*SIN(aziE*(pi/180))+ &
            radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    ELSE
       Keast=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    END IF

    CALL Kvikt_veg( svfS,svfSveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    IF (azimuth > (90-t) .AND. azimuth <= (270-t)) THEN
       Ksouth=radI*shadow*COS(altitude*(pi/180))*SIN(aziS*(pi/180))+ &
            radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    ELSE
       Ksouth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    END IF

    CALL Kvikt_veg( svfW,svfWveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    IF (azimuth > (180-t) .AND. azimuth <= (360-t)) THEN
       Kwest=radI*shadow*COS(altitude*(pi/180))*SIN(aziW*(pi/180))+ &
            radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    ELSE
       Kwest=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    END IF

    CALL Kvikt_veg( svfN,svfNveg,vikttot)
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    IF (azimuth <= (90-t) .OR. azimuth > (270-t)) THEN
       Knorth=radI*shadow*COS(altitude*(pi/180))*SIN(aziN*(pi/180))+ &
            radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    ELSE
       Knorth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    END IF

    DEALLOCATE(svfviktbuveg)
  END SUBROUTINE Kside_veg_v24


  SUBROUTINE Kvikt_veg(isvf, isvfveg, vikttot)


    IMPLICIT NONE
    REAL(KIND(1D0)) :: vikttot
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: isvf
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: isvfveg
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: svfvegbu

    !! Least
    viktwall=(vikttot-(63.227*isvf**6-161.51*isvf**5+156.91*isvf**4-70.424*isvf**3+16.773*isvf**2-0.4863*isvf))/vikttot

    svfvegbu=(isvfveg+isvf-1) ! Vegetation plus buildings
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2&
         &-0.4863*svfvegbu))/vikttot
    viktveg=viktveg-viktwall
  END SUBROUTINE Kvikt_veg

  SUBROUTINE Lside_veg_v2(altitude,Ta,Tw,SBC,ewall,esky,t,CI)!azimuth,
    ! This m-file is the current one that estimates L from the four cardinal points 20100414
    IMPLICIT NONE
    REAL(KIND(1D0))                             :: vikttot,aziE,aziN,aziS,aziW
    REAL(KIND(1D0))                             :: altitude,Ta,Tw,SBC,ewall,esky,t,CI,c!azimuth,
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  :: svfalfaE,svfalfaS,svfalfaW,svfalfaN
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  :: alfaB,betaB,betasun
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  :: Lground,Lrefl,Lsky,Lsky_allsky,Lveg,Lwallsh,Lwallsun
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  :: viktonlywall,viktaveg,svfvegbu
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  :: oneminussvfE,oneminussvfS,oneminussvfW,oneminussvfN
    REAL(KIND(1d0)),PARAMETER                   :: pi=3.141592653589793

    ALLOCATE(oneminussvfE(sizex,sizey))
    ALLOCATE(oneminussvfS(sizex,sizey))
    ALLOCATE(oneminussvfW(sizex,sizey))
    ALLOCATE(oneminussvfN(sizex,sizey))
    ALLOCATE(svfalfaE(sizex,sizey))
    ALLOCATE(svfalfaS(sizex,sizey))
    ALLOCATE(svfalfaW(sizex,sizey))
    ALLOCATE(svfalfaN(sizex,sizey))
    ALLOCATE(alfaB(sizex,sizey))
    ALLOCATE(betaB(sizex,sizey))
    ALLOCATE(betasun(sizex,sizey))
    ALLOCATE(Lground(sizex,sizey))
    ALLOCATE(Lrefl(sizex,sizey))
    ALLOCATE(Lsky(sizex,sizey))
    ALLOCATE(Lsky_allsky(sizex,sizey))
    ALLOCATE(Lveg(sizex,sizey))
    ALLOCATE(Lwallsh(sizex,sizey))
    ALLOCATE(Lwallsun(sizex,sizey))
    ALLOCATE(viktonlywall(sizex,sizey))
    ALLOCATE(viktaveg(sizex,sizey))
    ALLOCATE(svfvegbu(sizex,sizey))

    IF (ALLOCATED(viktwall)) DEALLOCATE(viktwall); ALLOCATE(viktwall(sizey,sizex))
    IF (ALLOCATED(viktsky)) DEALLOCATE(viktsky); ALLOCATE(viktsky(sizey,sizex))
    IF (ALLOCATED(viktveg)) DEALLOCATE(viktveg); ALLOCATE(viktveg(sizey,sizex))
    IF (ALLOCATED(viktrefl)) DEALLOCATE(viktrefl); ALLOCATE(viktrefl(sizey,sizex))

    !Building height angle from svf
    oneminussvfE=1.-svfE; WHERE (oneminussvfE<=0) oneminussvfE=0.000000001 ! avoiding log(0)
    oneminussvfS=1.-svfS; WHERE (oneminussvfS<=0) oneminussvfS=0.000000001 ! avoiding log(0)
    oneminussvfW=1.-svfW; WHERE (oneminussvfW<=0) oneminussvfW=0.000000001 ! avoiding log(0)
    oneminussvfN=1.-svfN; WHERE (oneminussvfN<=0) oneminussvfN=0.000000001 ! avoiding log(0)
    svfalfaE=ASIN(EXP((LOG(oneminussvfE))/2))
    svfalfaS=ASIN(EXP((LOG(oneminussvfS))/2))
    svfalfaW=ASIN(EXP((LOG(oneminussvfW))/2))
    svfalfaN=ASIN(EXP((LOG(oneminussvfN))/2))

    vikttot=4.4897
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t

    ! F_sh=cylindric_wedge(zen,svfalfa);!Fraction shadow on building walls based on sun altitude and svf
    ! F_sh(isnan(F_sh))=0.5;
    F_sh=2.*F_sh-1. !(cylindric_wedge scaled 0-1)

    IF (SOLWEIG_ldown==1) THEN
       c=1-CI
       Lsky_allsky=esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4)
    ELSE
       Lsky_allsky=ldown
    END IF

    !! Least
    CALL Lvikt_veg(svfE,svfEveg,svfEaveg,vikttot)

    IF (altitude>0) THEN ! daytime
       alfaB=ATAN(svfalfaE)
       betaB=ATAN(TAN((svfalfaE)*F_sh))
       betasun=((alfaB-betaB)/2)+betaB
       IF (azimuth > (180-t) .AND. azimuth <= (360-t)) THEN
          Lwallsun=SBC*ewall*((Ta+273.15+Tw*SIN(aziE*(pi/180)))**4)* &
               viktwall*(1-F_sh)*COS(betasun)*0.5
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
       ELSE
          Lwallsun=0
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
       END IF
    ELSE !nighttime
       Lwallsun=0
       Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    END IF
    Lsky=((svfE+svfEveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Least=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lsouth
    CALL Lvikt_veg(svfS,svfSveg,svfSaveg,vikttot)

    IF (altitude>0) THEN ! daytime
       alfaB=ATAN(svfalfaS)
       betaB=ATAN(TAN((svfalfaS)*F_sh))
       betasun=((alfaB-betaB)/2)+betaB
       IF (azimuth <= (90-t) .OR. azimuth > (270-t)) THEN
          Lwallsun=SBC*ewall*((Ta+273.15+Tw*SIN(aziS*(pi/180)))**4)* &
               viktwall*(1-F_sh)*COS(betasun)*0.5
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
       ELSE
          Lwallsun=0
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
       END IF
    ELSE !nighttime
       Lwallsun=0
       Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    END IF
    Lsky=((svfS+svfSveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lsouth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lwest
    CALL Lvikt_veg(svfW,svfWveg,svfWaveg,vikttot)

    IF (altitude>0) THEN ! daytime
       alfaB=ATAN(svfalfaW)
       betaB=ATAN(TAN((svfalfaW)*F_sh))
       betasun=((alfaB-betaB)/2)+betaB
       IF (azimuth > (360-t) .OR. azimuth <= (180-t)) THEN
          Lwallsun=SBC*ewall*((Ta+273.15+Tw*SIN(aziW*(pi/180)))**4)* &
               viktwall*(1-F_sh)*COS(betasun)*0.5
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
       ELSE
          Lwallsun=0
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
       END IF
    ELSE !nighttime
       Lwallsun=0
       Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    END IF
    Lsky=((svfW+svfWveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lwest=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    !! Lnorth
    CALL Lvikt_veg(svfN,svfNveg,svfNaveg,vikttot)

    IF (altitude>0) THEN ! daytime
       alfaB=ATAN(svfalfaN)
       betaB=ATAN(TAN((svfalfaN)*F_sh))
       betasun=((alfaB-betaB)/2)+betaB
       IF (azimuth > (90-t) .AND. azimuth <= (270-t)) THEN
          Lwallsun=SBC*ewall*((Ta+273.15+Tw*SIN(aziN*(pi/180)))**4)* &
               viktwall*(1-F_sh)*COS(betasun)*0.5
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
       ELSE
          Lwallsun=0
          Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
       END IF
    ELSE !nighttime
       Lwallsun=0
       Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    END IF
    Lsky=((svfN+svfNveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lnorth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    DEALLOCATE(svfalfaE)
    DEALLOCATE(svfalfaS)
    DEALLOCATE(svfalfaW)
    DEALLOCATE(svfalfaN)
    DEALLOCATE(alfaB)
    DEALLOCATE(betaB)
    DEALLOCATE(betasun)
    DEALLOCATE(Lground)
    DEALLOCATE(Lrefl)
    DEALLOCATE(Lsky)
    DEALLOCATE(Lsky_allsky)
    DEALLOCATE(Lveg)
    DEALLOCATE(Lwallsh)
    DEALLOCATE(Lwallsun)
    DEALLOCATE(viktonlywall)
    DEALLOCATE(viktaveg)
    DEALLOCATE(svfvegbu)

  END SUBROUTINE Lside_veg_v2

  SUBROUTINE Lvikt_veg(isvf,isvfveg,isvfaveg,vikttot)

    IMPLICIT NONE
    REAL(KIND(1D0)) :: vikttot
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: isvf
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: isvfveg
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: isvfaveg
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: viktonlywall
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: viktaveg
    !real(kind(1d0)), dimension(sizey,sizex) :: viktwall
    REAL(KIND(1d0)), DIMENSION(sizey,sizex) :: svfvegbu

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

  END SUBROUTINE Lvikt_veg
  !  FUNCTION TO RETURN 0 IF IX=0, 1 IF 0<IX<MAXPOS+1,-1 OTHERWISE.
  !  MAXPOS is given as the maximum positive Integer.
  SUBROUTINE issign(IX,MAXPOS,ISIGNM)
    REAL(KIND(1d0)) IX,MAXPOS,ISIGNM
    ISIGNM=1.0
    IF(IX<0.OR.IX>MAXPOS)ISIGNM=-1
    IF(IX==0) ISIGNM=0
    RETURN
  END SUBROUTINE issign


  !=====================================================
  ! fuction to convert interger to string
  CHARACTER(len=20) FUNCTION str(k)
    INTEGER, INTENT(in) :: k
    WRITE (str, *) k
    str = ADJUSTL(str)
  END FUNCTION str


  SUBROUTINE SaveGrids

    IMPLICIT NONE
    CHARACTER(len=100)       :: GridPath,GridName,GridPath2
    CHARACTER(len=4)       ::doy,hour
    !real(kind(1d0)), allocatable, dimension(:,:):: savegrid

    ALLOCATE(savegrid(sizey,sizex))

    IF (Tmrt_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='Tmrt_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=Tmrt
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF
    IF (Lup2d_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='Tmrt_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=Tmrt
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF
    IF (Ldown2d_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='Ldown2d_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=Ldown2d
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF
    IF (Kup2d_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='Kup2d_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=Kup2d
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF
    IF (Kdown2d_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='Kdown2d_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=Kdown2d
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF
    IF (GVF_out==1) THEN
       Gridpath2='Grids/'
       GridPath=TRIM(FileOutputPath)//TRIM(GridPath2)
       WRITE(doy,'(i3)') id
       WRITE(hour,'(i2)') it
       hour=ADJUSTL(hour)
       doy=ADJUSTL(doy)
       GridName='GVF_'//TRIM(doy)//'_'//TRIM(hour)//'.txt'
       savegrid=gvf
       CALL SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    END IF

    DEALLOCATE(savegrid)

  END SUBROUTINE SaveGrids
  !------------------------------------------------!
  ! Shadow casting algorithm, ground and buildings !
  !------------------------------------------------!
  SUBROUTINE shadowingfunction_10(azimuth,altitude,scale)
    ! This m.file calculates shadows on a DSM
    ! Code originate from Carlo Rattis thesis
    ! This code is translated from Matlab by
    ! Fredrik Lindberg, Gothenburg University
    ! Last modified LJ 27 Jan 2016 - Removal of tabs and fixing real-int conversions

    IMPLICIT NONE
    REAL(KIND(1d0)), PARAMETER          :: pi=3.141592653589793
    REAL(KIND(1d0)), PARAMETER          :: maxpos=10000000000.0
    !real, allocatable, dimension(:,:)   :: a,f,temp !! Defined in matsize
    REAL(KIND(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    REAL(KIND(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    REAL(KIND(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    REAL(KIND(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    INTEGER                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,row,col !,sizex,sizey
    ! Internal grids
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: temp,f

    !special cases
    IF (altitude==90) THEN
       altitude=altitude-0.0001
    END IF
    IF (azimuth==0) THEN
       azimuth=azimuth-0.0001
    END IF

    ! conversion
    degrees=pi/180
    azi=azimuth*degrees
    alt=altitude*degrees

    ALLOCATE(f(sizex,sizey))
    ALLOCATE(temp(sizex,sizey))

    IF (ALLOCATED(sh)) DEALLOCATE(sh)
    ALLOCATE(sh(sizex,sizey))

    ! initialise parameters
    f=a
    dx=0
    dy=0
    dz=0
    temp=a*0.0
    index=1

    ! other loop parameters
    amaxvalue=MAXVAL(a)
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=SIN(azi)
    cosazimuth=COS(azi)
    tanazimuth=TAN(azi)
    CALL issign(sinazimuth,maxpos,signsinazimuth)
    CALL issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=ABS(1./sinazimuth)
    dscos=ABS(1./cosazimuth)
    tanaltitudebyscale=TAN(alt)/scale


    ! main loop
    DO WHILE (amaxvalue>=dz .AND. ABS(dx)<=sizex .AND. ABS(dy)<=sizey)

       IF ((pibyfour <= azi .AND. azi < threetimespibyfour) .OR. (fivetimespibyfour <= azi .AND. azi < seventimespibyfour)) THEN
          dy=signsinazimuth*index
          dx=-1.*signcosazimuth*ABS(NINT(index/tanazimuth))
          ds=dssin
       ELSE
          dy=signsinazimuth*ABS(NINT(index*tanazimuth))
          dx=-1.*signcosazimuth*index
          ds=dscos
       END IF

       dz=ds*index*tanaltitudebyscale
       temp=temp*0

       absdx=ABS(dx)
       absdy=ABS(dy)

       xc1=INT((dx+absdx)/2)+1
       xc2=(sizex+INT((dx-absdx)/2))
       yc1=INT((dy+absdy)/2)+1
       yc2=(sizey+INT((dy-absdy)/2))
       xp1=-INT((dx-absdx)/2)+1
       xp2=(sizex-INT((dx+absdx)/2))
       yp1=-INT((dy-absdy)/2)+1
       yp2=(sizey-INT((dy+absdy)/2))

       temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz

       f=MAX(f,temp)
       index=index+1

    END DO

    f=f-a
    WHERE (f>0)
       f=-1
    END WHERE
    sh=f+1
    !sh=f ! invert as in shadowingfunctionglobalradiation
    DEALLOCATE(f)
    DEALLOCATE(temp)

  END SUBROUTINE shadowingfunction_10
  !------------------------------------------------!
  ! Shadow casting algorithm, vegetation           !
  !------------------------------------------------!

  SUBROUTINE shadowingfunction_20(azimuth,altitude,scale,amaxvalue)

    ! This m.file calculates shadows on a DSM and for vegetation units
    ! This code is translated from Matlab by Fredrik Lindberg, Gothenburg University

    IMPLICIT NONE
    REAL(KIND(1d0)), PARAMETER          :: pi=3.141592653589793
    REAL(KIND(1d0)), PARAMETER          :: maxpos=10000000000.0
    REAL(KIND(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    REAL(KIND(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    REAL(KIND(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    REAL(KIND(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    INTEGER                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,test!,row,col !,sizex,sizey
    ! Internal grids
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: f,temp,tmp,stopbuild,stopveg,g,bushplant,tempvegdem,tempvegdem2
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: fabovea,gabovea,tempbush,firstvegdem,vegsh2

    !real 		  				  	    :: start_time,end_time

    !special case
    IF (altitude==90) THEN
       altitude=altitude-0.0001
    END IF
    IF (azimuth==0) THEN
       azimuth=azimuth-0.0001
    END IF

    ! conversion
    degrees=pi/180;
    azi=azimuth*degrees;
    alt=altitude*degrees;

    IF (ALLOCATED(sh)) DEALLOCATE(sh)
    ALLOCATE(sh(sizex,sizey))
    IF (ALLOCATED(vegsh)) DEALLOCATE(vegsh)
    ALLOCATE(vegsh(sizex,sizey))
    IF (ALLOCATED(vbshvegsh)) DEALLOCATE(vbshvegsh)
    ALLOCATE(vbshvegsh(sizex,sizey))

    ! allocation of grids
    ALLOCATE(f(sizex,sizey))
    ALLOCATE(temp(sizex,sizey))
    ALLOCATE(tmp(sizex,sizey))
    ALLOCATE(stopbuild(sizex,sizey))
    ALLOCATE(stopveg(sizex,sizey))
    ALLOCATE(g(sizex,sizey))
    ALLOCATE(bushplant(sizex,sizey))
    ALLOCATE(tempvegdem(sizex,sizey))
    ALLOCATE(tempvegdem2(sizex,sizey))
    ALLOCATE(fabovea(sizex,sizey))
    ALLOCATE(gabovea(sizex,sizey))
    ALLOCATE(firstvegdem(sizex,sizey))
    ALLOCATE(tempbush(sizex,sizey))
    ALLOCATE(vegsh2(sizex,sizey))

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
    WHERE (bush>1)
       bushplant=1
    END WHERE

    index=1
    !test=0
    ! other loop parameters
    !amaxvalue=maxval(a)
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour;
    fivetimespibyfour=5.*pibyfour;
    seventimespibyfour=7.*pibyfour;
    sinazimuth=SIN(azi);
    cosazimuth=COS(azi);
    tanazimuth=TAN(azi);
    CALL issign(sinazimuth,maxpos,signsinazimuth)
    CALL issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth);
    !signcosazimuth=cosazimuth/abs(cosazimuth);
    dssin=ABS(1./sinazimuth);
    dscos=ABS(1./cosazimuth);
    tanaltitudebyscale=TAN(alt)/scale;

    DO WHILE (amaxvalue>=dz .AND. ABS(dx)<=sizex .AND. ABS(dy)<=sizey)

       IF ((pibyfour <= azi .AND. azi < threetimespibyfour) .OR. (fivetimespibyfour <= azi .AND. azi < seventimespibyfour)) THEN
          dy=signsinazimuth*index
          dx=-1.*signcosazimuth*ABS(NINT(index/tanazimuth))
          ds=dssin
       ELSE
          dy=signsinazimuth*ABS(NINT(index*tanazimuth))
          dx=-1.*signcosazimuth*index
          ds=dscos
       END IF

       dz=ds*index*tanaltitudebyscale
       temp=temp*0
       tempvegdem=temp
       tempvegdem2=temp

       absdx=ABS(dx)
       absdy=ABS(dy)

       xc1=INT(((dx+absdx)/2))+1
       xc2=(sizex+INT((dx-absdx)/2))
       yc1=INT((dy+absdy)/2)+1
       yc2=(sizey+INT((dy-absdy)/2))
       xp1=-INT((dx-absdx)/2)+1
       xp2=(sizex-INT((dx+absdx)/2))
       yp1=-INT((dy-absdy)/2)+1
       yp2=(sizey-INT((dy+absdy)/2))

       temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz
       tempvegdem(xp1:xp2,yp1:yp2)=vegdem(xc1:xc2,yc1:yc2)-dz
       tempvegdem2(xp1:xp2,yp1:yp2)=vegdem2(xc1:xc2,yc1:yc2)-dz

       f=MAX(f,temp)
       WHERE (f>a) !sh(f>a)=1;sh(f<=a)=0; !Moving building shadow
          sh=1
       ELSEWHERE
          sh=0
       END WHERE
       WHERE (tempvegdem>a) !fabovea=tempvegdem>a; !vegdem above DEM
          fabovea=1
       ELSEWHERE
          fabovea=0
       END WHERE
       WHERE (tempvegdem2>a) !gabovea=tempvegdem2>a; !vegdem2 above DEM
          gabovea=1
       ELSEWHERE
          gabovea=0
       END WHERE
       vegsh2=fabovea-gabovea
       vegsh=MAX(vegsh,vegsh2)
       WHERE ((vegsh*sh)>0) !vegsh(vegsh.*sh>0)=0;! removing shadows 'behind' buildings
          vegsh=0
       END WHERE
       vbshvegsh=vegsh+vbshvegsh

       ! vegsh at high sun altitudes
       IF (index==1) THEN
          firstvegdem=tempvegdem-temp
          WHERE (firstvegdem<=0)!firstvegdem(firstvegdem<=0)=1000;
             firstvegdem=1000
          END WHERE
          WHERE (firstvegdem<dz)!vegsh(firstvegdem<dz)=1;
             vegsh=1
          END WHERE
          tmp=temp*0.0
          WHERE (vegdem2>a)!vegsh=vegsh.*(vegdem2>a);
             tmp=1
          END WHERE
          vegsh=vegsh*tmp
          vbshvegsh=temp*0.0 !vbshvegsh=zeros(sizex,sizey);
       END IF

       ! Bush shadow on bush plant
       tmp=fabovea*bush
       IF ((MAXVAL(bush)>0) .AND. (MAXVAL(tmp)>0)) THEN
          tempbush=temp*0.0
          tempbush(xp1:xp2,yp1:yp2)=bush(xc1:xc2,yc1:yc2)-dz
          g=MAX(g,tempbush)
          g=bushplant*g
       END IF

       index=index+1
    END DO

    sh=1-sh
    WHERE (vbshvegsh>0)!vbshvegsh(vbshvegsh>0)=1;
       vbshvegsh=1
    END WHERE
    vbshvegsh=vbshvegsh-vegsh;

    IF (MAXVAL(bush)>0) THEN
       g=g-bush
       WHERE (g>0)!g(g>0)=1;g(g<0)=0;
          g=1
       ELSEWHERE
          g=0
       END WHERE
       vegsh=vegsh-bushplant+g
       WHERE (vegsh<0)!vegsh(vegsh<0)=0;
          vegsh=0
       END WHERE
    END IF

    WHERE (vegsh>0)!vegsh(vegsh>0)=1;
       vegsh=1
    END WHERE
    vegsh=1-vegsh
    vbshvegsh=1-vbshvegsh

    !deallocation of grids
    DEALLOCATE(f)
    DEALLOCATE(temp)
    DEALLOCATE(tmp)
    DEALLOCATE(stopbuild)
    DEALLOCATE(stopveg)
    DEALLOCATE(g)
    DEALLOCATE(bushplant)
    DEALLOCATE(tempvegdem)
    DEALLOCATE(tempvegdem2)
    DEALLOCATE(fabovea)
    DEALLOCATE(gabovea)
    DEALLOCATE(firstvegdem)
    DEALLOCATE(tempbush)
    DEALLOCATE(vegsh2)

  END SUBROUTINE shadowingfunction_20
  SUBROUTINE sunonsurface_veg(iazimuthA, scale, first, second, psi)
    ! This m-file creates a boolean image of sunlit walls.
    ! Shadows from both buildings and vegetation is accounted for
    ! moving building in the direction of the sun
    ! Last modified:
    ! LJ 27 Jan 2016 - Removal of tabs and fixing real-int conversions
    IMPLICIT NONE
    REAL(KIND(1d0))             :: iazimuthA,iazimuth,sinazimuth,cosazimuth,tanazimuth
    REAL(KIND(1d0))             :: scale
    INTEGER                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2,n,first,second
    REAL(KIND(1d0))             :: dx,dy,ds,absdx,absdy,psi !,dz
    REAL(KIND(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    REAL(KIND(1d0))             :: seventimespibyfour
    REAL(KIND(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  ::weightsumwall,weightsumsh,gvf1,gvf2
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)  ::f,tempsh,tempbu,tempbub,tempwallsun,tempb,sh1

    REAL(KIND(1d0)), PARAMETER  :: pi=3.141592653589793
    REAL(KIND(1d0)), PARAMETER  :: maxpos=10000000000.0

    ALLOCATE(weightsumwall(sizex,sizey))
    ALLOCATE(weightsumsh(sizex,sizey))
    ALLOCATE(gvf1(sizex,sizey))
    ALLOCATE(gvf2(sizex,sizey))
    ALLOCATE(f(sizex,sizey))
    ALLOCATE(tempsh(sizex,sizey))
    ALLOCATE(tempbu(sizex,sizey))
    ALLOCATE(tempbub(sizex,sizey))
    ALLOCATE(tempwallsun(sizex,sizey))
    ALLOCATE(tempb(sizex,sizey))
    ALLOCATE(sh1(sizex,sizey))

    iazimuth=iazimuthA*(pi/180)
    !special cases
    IF (iazimuth==0) THEN
       iazimuth=iazimuth+0.000001
    END IF
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

    first=INT(REAL(first,KIND(1d0))*scale)  !Int added around the equation as first and second are both integers
    second=INT(REAL(second,KIND(1d0))*scale)

    ! other loop parameters
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=SIN(iazimuth)
    cosazimuth=COS(iazimuth)
    tanazimuth=TAN(iazimuth)
    CALL issign(sinazimuth,maxpos,signsinazimuth)
    CALL issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=ABS(1./sinazimuth)
    dscos=ABS(1./cosazimuth)

    !! The Shadow casting algorithm
    DO n=1,second
       IF ((pibyfour <= iazimuth .AND. iazimuth < threetimespibyfour) .OR. (fivetimespibyfour&
            & <= iazimuth .AND. iazimuth < seventimespibyfour)) THEN
          dy=signsinazimuth*index
          dx=-1.*signcosazimuth*ABS(NINT(index/tanazimuth))
          ds=dssin
       ELSE
          dy=signsinazimuth*ABS(NINT(index*tanazimuth))
          dx=-1.*signcosazimuth*index
          ds=dscos
       END IF

       absdx=ABS(dx)
       absdy=ABS(dy)

       xc1=INT((dx+absdx)/2)+1
       xc2=(sizex+INT((dx-absdx)/2))
       yc1=INT((dy+absdy)/2)+1
       yc2=(sizey+INT((dy-absdy)/2))
       xp1=-INT((dx-absdx)/2)+1
       xp2=(sizex-INT((dx+absdx)/2))
       yp1=-INT((dy-absdy)/2)+1
       yp2=(sizey-INT((dy+absdy)/2))

       tempbu(xp1:xp2,yp1:yp2)=buildings(xc1:xc2,yc1:yc2) !moving building

       tempsh(xp1:xp2,yp1:yp2)=sh1(xc1:xc2,yc1:yc2) !moving shadow image
       f=MIN(f,tempbu) !utsmetning of buildings

       weightsumsh=weightsumsh+tempsh*f

       tempwallsun(xp1:xp2,yp1:yp2)=sunwall(xc1:xc2,yc1:yc2) !moving building wall in sun image
       tempb=tempwallsun*f
       WHERE ((tempb+tempbub)>0) !tempbub=(tempb+tempbub)>0==1
          tempbub=1.
       END WHERE

       weightsumwall=weightsumwall+tempbub

       IF (index*scale==first) THEN
          gvf1=(weightsumwall+weightsumsh)/first
          WHERE (gvf1>1)
             gvf1=1.
          END WHERE
       END IF
       index=index+1

    END DO
    gvf2=(weightsumsh+weightsumwall)/second
    WHERE (gvf2>1)
       gvf2=1.
    END WHERE

    ! Weighting
    sos=(gvf1*0.5+gvf2*0.4)/0.9

    DEALLOCATE(weightsumwall)
    DEALLOCATE(weightsumsh)
    DEALLOCATE(gvf1)
    DEALLOCATE(gvf2)
    DEALLOCATE(f)
    DEALLOCATE(tempsh)
    DEALLOCATE(tempbu)
    DEALLOCATE(tempbub)
    DEALLOCATE(tempwallsun)
    DEALLOCATE(tempb)
    DEALLOCATE(sh1)

  END SUBROUTINE sunonsurface_veg

  !>
  SUBROUTINE wallinsun_veg(azimuth)
    ! This m-file creates a boolean image of sunlit walls.
    ! Shadows from both buildings and vegetation is accounted for
    ! moving building in the direction of the sun
    ! Last modified:
    !  LJ 27 Jan 2017 - Change of equations xc1...yp2 to account for the change from real to integer
    !---------------------------------------------------------------------------------


    IMPLICIT NONE
    REAL(KIND(1d0))             :: azimuth,iazimuth
    INTEGER                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2
    REAL(KIND(1d0))             :: dx,dy,dz,ds,absdx,absdy
    REAL(KIND(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    REAL(KIND(1d0))             :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    REAL(KIND(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos,azi
    REAL(KIND(1d0)), PARAMETER  :: pi=3.141592653589793
    REAL(KIND(1d0)), PARAMETER  :: maxpos=10000000000.0
    ! Internal grids
    REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:) :: temp,sh1

    ALLOCATE(temp(sizex,sizey))
    ALLOCATE(sh1(sizex,sizey))
    !allocate(vegsh(sizex,sizey))
    !allocate(a(sizex,sizey))
    !allocate(buildings(sizex,sizey))

    IF (ALLOCATED(sunwall)) DEALLOCATE(sunwall); ALLOCATE(sunwall(sizey,sizex))

    iazimuth=azimuth+180
    IF (iazimuth>=360) THEN
       iazimuth=iazimuth-360
    END IF
    !special cases
    IF (iazimuth==0) THEN
       iazimuth=iazimuth+0.00001
    END IF
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
    sinazimuth=SIN(azi)
    cosazimuth=COS(azi)
    tanazimuth=TAN(azi)
    CALL issign(sinazimuth,maxpos,signsinazimuth)
    CALL issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=ABS(1./sinazimuth)
    dscos=ABS(1./cosazimuth)

    sh1=vegsh+sh-1.
    !! The Shadow casting algoritm
    IF ((pibyfour <= azi .AND. azi < threetimespibyfour) .OR. (fivetimespibyfour <= azi .AND. azi < seventimespibyfour)) THEN
       dy=signsinazimuth*index
       dx=-1.*signcosazimuth*ABS(NINT(index/tanazimuth))
       ds=dssin
    ELSE
       dy=signsinazimuth*ABS(NINT(index*tanazimuth))
       dx=-1.*signcosazimuth*index
       ds=dscos
    END IF

    ! note: dx and dy represent absolute values while ds is an incremental value

    absdx=ABS(dx)
    absdy=ABS(dy)

    xc1=INT((dx+absdx)/2)+1  !LJ added int to the equation to account for the conversion from real to int
    xc2=(sizex+INT((dx-absdx)/2))
    yc1=INT((dy+absdy)/2)+1
    yc2=(sizey+INT((dy-absdy)/2))
    xp1=-INT((dx-absdx)/2)+1
    xp2=(sizex-INT((dx+absdx)/2))
    yp1=-INT((dy-absdy)/2)+1
    yp2=(sizey-INT((dy+absdy)/2))

    temp(xp1:xp2,yp1:yp2)= buildings(xc1:xc2,yc1:yc2)

    sunwall=temp-buildings
    WHERE (sunwall==1) !f1(f1==1)=0
       sunwall=0
    END WHERE
    WHERE (sunwall==-1) !f1(f1==-1)=1
       sunwall=1
    END WHERE
    sunwall=sh1*sunwall

    DEALLOCATE(temp)
    DEALLOCATE(sh1)

  END SUBROUTINE wallinsun_veg


END MODULE solweig_module
