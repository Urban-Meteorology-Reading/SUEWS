MODULE NARP_MODULE
  !==============================================================================
  !NET ALL WAVE RADIATION PARAMETERIZATION ROUTINES
  !B. OFFERLE
  !DEPT OF GEOGRAPHY
  !INDIANA UNIVERSITY
  !bofferle@indiana.edu
  !
  !MODIFIED: 19 DEC 2002
  !CURRENTLY THE SMITH GRID IS ONLY VALID FOR THE N. HEMISPHERE
  !
  !Thomas Loridan, May 2008
  !4.1: MODIFICATION FOR CLOUD FRACTION PARAMETERIZATION AT NIGHT USING THE RATE OF COOLING.
  !     EOSHIFT INTRINSIC FUNCTION WAS ALSO REMOVED BECAUSE IT IS COMPILER DEPENDENT.
  !
  !6.0  T. Loridan - June 2009
  !     Different longwave down options (ldown_option)
  ! 1 - Ldown Observed (add data as last column in met file)
  ! 2 - Ldown modelled from observed FCLD (add data as last column in met file)
  ! 3 - Ldown modelled from FCLD(RH,TA)
  ! 4 - Ldown modelled from FCLD(Kdown); i.e. day FCLD only
  !     cloud fraction is kept constant throught the night (Offerle et al. 2003, JAM)
  ! 5 - Option 3 at night and 4 during the day (might cause discontinuities in Ldown)

  !SUEWS   L. Jarvi - Oct 2010
  !Currently Ldown options 4 and 5 commented out in order to reduce input files.
  !Last modified:
  ! TS 06 Aug 2017 - interface modified to receive explict input and output arguments
  ! LJ 27 Jan 2016 - Removal of tabs, cleaning of the code
  ! FL July 2014 - Variables are moved to modules in NARP subroutine. Snow related should also in future.
  ! FL Nov 2013 - A new sun postion algorithm added
  ! LJ May 2013 - Main program NARP changed to take subsurfaces and snow into account here and not
  ! in the main program
  ! LJ Oct 2012 - Zenith angle change in the calculation of albedo added
  ! sg feb 2012 - Allocatable array module added

  !==============================================================================================
  ! USE allocateArray

  IMPLICIT NONE

CONTAINS

  !==============================================================================
  SUBROUTINE NARP(&
       nsurf,sfr,snowFrac,alb,emis,IceFrac,&! input:
       NARP_G,NARP_TRANS_SITE,NARP_EMIS_SNOW,&
       DTIME,ZENITH_deg,kdown,Temp_C,RH,Press_hPa,qn1_obs,&
       SnowAlb,&
       AlbedoChoice,ldown_option,NetRadiationMethod,DiagQN,&
       QSTARall,QSTAR_SF,QSTAR_S,kclear,KUPall,LDown,LUPall,fcld,TSURFall)! output:
    !KCLEAR,FCLD,DTIME,KDOWN,QSTARall,KUPall,LDOWN,LUPall,TSURFall,&
    !AlbedoChoice,ldown_option,Temp_C,Press_hPa,Ea_hPa,qn1_obs,RH,&
    !,zenith_degnetRadiationChoice,

    !SUBROUTINE NARP(QSTAR,DTIME,KDOWN,LDOWN,T,RH,PRESS,FCLD,SNOW)
    !returns estimate of Q* given the meteorological fields and prior
    !configuration call.

    !OUTPUT FIELDS
    !QSTARall (W/m2) = net all wave radiation
    !KCLEAR          = clear sky incoming solar radiation
    !KUPall          = reflect solar radiation
    !LDOWN           = incoming longwave radiation (observed or modelled depending on ldown_option)
    !LUPall          = outgoing longwaver radiation
    !FCLD            = estimated cloud fraction (used only for emissivity estimate)
    !                  FCLD USED AS INPUT ALSO
    !TSURFall (DEG C)= estimated surface temperature
    !QSTAR_SF        = net all wave radiation for snow free surface
    !QSTAR_S         = net all wave radiation for SnowPack

    !INPUT FIELDS
    !DTIME (days) = local time, not daylight savings
    !KDOWN (W/m2) = incoming solar radiation
    !T (DEG C)    = air temperature near model height
    !RH (%)       = relative humidity near model height
    !PRESS (mb)   = station pressure, use estimate if unavailable
    !qn1_obs      = Observed Q*


    !INTERNAL FIELDS
    ! TemP_K = air temperature in K
    ! ES_hPs = saturation vapor pressure (hPa)
    ! EA_Pa = vapor pressure (hPa)
    ! TD = dewpoint (C)
    ! ZENITH = solar zenith angle
    ! ALB0 = surface albedo
    ! EMIS0 = surface emissivity
    ! EMIS_A = atmospheric emissivity
    ! TRANS = atmospheric transmissivity
    ! LUPCORR = correction for surface heating by kdown (W/m2)
    ! SIGMATK4 = energy flux density
    ! KDOWN_HR = hourly average insolation
    ! DOY = day of year

    !Modified by LJ to calcuate snow free and SnowPack components (May 2013)
    !Modified to define variables in data_in module
    !-------------------------------------------------------------------------------
    ! USE allocateArray
    ! use gis_data
    ! use data_in ! Included 20140701, FL
    ! use moist   ! Included 20140701, FL
    ! use time    ! Included 20140701, FL
    REAL(KIND(1D0)),DIMENSION(nsurf),INTENT(in) ::sfr
    REAL(KIND(1D0)),DIMENSION(nsurf),INTENT(in) ::snowFrac
    REAL(KIND(1D0)),DIMENSION(nsurf),INTENT(in) ::alb
    REAL(KIND(1D0)),DIMENSION(nsurf),INTENT(in) ::emis
    REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in) ::IceFrac
    REAL(KIND(1D0)),DIMENSION(365),INTENT(in) ::NARP_G

    REAL(KIND(1D0)),INTENT(in) ::DTIME
    REAL(KIND(1D0)),INTENT(in) ::ZENITH_deg
    REAL(KIND(1D0)),INTENT(in) ::kdown
    REAL(KIND(1D0)),INTENT(in) ::Temp_C
    REAL(KIND(1D0)),INTENT(in) ::RH
    REAL(KIND(1D0)),INTENT(in) ::Press_hPa
    REAL(KIND(1D0)),INTENT(in) ::qn1_obs
    REAL(KIND(1D0)),INTENT(in) ::SnowAlb
    REAL(KIND(1D0)),INTENT(in) ::NARP_TRANS_SITE
    REAL(KIND(1D0)),INTENT(in) ::NARP_EMIS_SNOW

    INTEGER,INTENT(in) ::nsurf
    INTEGER,INTENT(in) ::AlbedoChoice
    INTEGER,INTENT(in) ::ldown_option
    INTEGER,INTENT(in) ::NetRadiationMethod
    INTEGER,INTENT(in) ::DiagQN

    REAL(KIND(1D0)),INTENT(out) ::QSTARall
    REAL(KIND(1D0)),INTENT(out) ::QSTAR_SF
    REAL(KIND(1D0)),INTENT(out) ::QSTAR_S
    REAL(KIND(1D0)),INTENT(out) ::kclear
    REAL(KIND(1D0)),INTENT(out) ::KUPall
    REAL(KIND(1D0)),INTENT(out) ::ldown
    REAL(KIND(1D0)),INTENT(out) ::LUPall
    REAL(KIND(1D0)),INTENT(out) ::fcld
    REAL(KIND(1D0)),INTENT(out) ::TSURFall

    REAL(KIND(1d0)),DIMENSION(nsurf) ::qn1_ind
    REAL(KIND(1d0)),DIMENSION(nsurf) ::kup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf) ::lup_ind
    REAL(KIND(1d0)),DIMENSION(nsurf) ::Tsurf_ind

    REAL(KIND(1d0)),DIMENSION(nsurf) ::qn1_ind_nosnow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::kup_ind_nosnow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::lup_ind_nosnow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::Tsurf_ind_nosnow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::qn1_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::kup_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::lup_ind_snow
    REAL(KIND(1d0)),DIMENSION(nsurf) ::Tsurf_ind_snow

    REAL(KIND(1D0)) ::Temp_K,TD,ZENITH,QSTAR,QSTAR_SNOW,KUP_SNOW,LUP_SNOW,TSURF_SNOW,KUP,LUP,TSURF
    REAL(KIND(1D0)) ::ALB0,EMIS0,EMIS_A,TRANS!,RH,DTIME,KDOWN
    REAL(KIND(1D0)) ::LUPCORR,SIGMATK4,KDOWN_HR=0.
    INTEGER         ::DOY, is

    REAL(KIND(1D0))::qn1_cum,kup_cum,lup_cum,tsurf_cum,&   !Cumulative radiation components
         qn1_is,kup_is,lup_is,tsurf_is,&       !Sub-surface radiation components
         SF_all,ALB1

    REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292,&
        !  RAD2DEG=57.29577951,&
         SIGMA_SB=5.67E-8

    !Initialize variables
    ! RH=avrh
    ! DTIME=dectime
    ! KDOWN=avkdn
    Temp_K=Temp_C+273.16
    SIGMATK4=SIGMA_SB*Temp_K**4
    TD=DEWPOINT(Temp_C,RH)
    ! Sun postition is now calculated in the main loop, FL
    !ZENITH=SOLAR_ZENITH(NARP_LAT,NARP_LONG,NARP_TZ,DTIME)
    !call sun_position(NARP_YEAR,DTIME,NARP_TZ,NARP_LAT,NARP_LONG,Alt,AZIMUTH,ZENITH)
    ZENITH=ZENITH_deg*DEG2RAD
    DOY=INT(DTIME)
    IF(DOY==366)doy=365

    !===================================================================================
    !Calculate radiation for each sub-surface
    qn1_cum=0
    kup_cum=0
    lup_cum=0
    tsurf_cum=0

    QSTAR_SF=0
    QSTAR_S=0

    !Total snowfree surface fraction
    SF_all=0
    DO is = 1,nsurf
       IF (sfr(is)/=0) SF_all = SF_all + sfr(is)*(1-snowFrac(is))
    ENDDO

    DO is=1,nsurf

       EMIS_A=PRATA_EMIS(Temp_K,Press_hPa)

       !--------------------------------------------------
       !-------SNOW FREE SURFACE--------------------------

       IF (AlbedoChoice==1.AND.180*ZENITH/ACOS(0.0)<90) THEN
          ALB0=ALB(is)+0.5e-16*(180*ZENITH/ACOS(0.0))**8 !AIDA 1982
       ELSE
          ALB0=ALB(is)
       ENDIF
       EMIS0=EMIS(is)

       !Downward longwave radiation
       IF((ldown_option==4) .OR. (ldown_option==5)) THEN !Estimate FCLD from Kdown (Offerle et al. 2003)
          IF (ZENITH<1.5) THEN !DAYTIME CALCULATIONS
             TRANS=TRANSMISSIVITY(Press_hPa,TD,NARP_G(DOY),ZENITH)
             KCLEAR=ISURFACE(DOY,ZENITH)*TRANS*NARP_TRANS_SITE
             IF (KCLEAR>50.) THEN
                FCLD=CLOUD_FRACTION(KDOWN,KCLEAR)
                EMIS_A=EMIS_CLOUD_SQ(EMIS_A,FCLD)
             ELSE
                IF(ldown_option==5) THEN ! Use RH when Kdown can not be used
                   FCLD=WC_fraction(RH,Temp_C)
                   EMIS_A=EMIS_CLOUD(EMIS_A,FCLD)
                ELSE
                   !FCLD is left to the latest calculable value
                   EMIS_A=EMIS_CLOUD_SQ(EMIS_A,FCLD)
                ENDIF
             ENDIF
          ELSE !NIGHT TIME CALCULATIONS
             IF(ldown_option==4) THEN
                !FCLD is left to the latest calculable value
                EMIS_A=EMIS_CLOUD_SQ(EMIS_A,FCLD)
             ELSEIF((ldown_option==5)) THEN ! Use RH
                FCLD=WC_fraction(RH,Temp_C)
                EMIS_A=EMIS_CLOUD(EMIS_A,FCLD)
             ENDIF
          ENDIF
       ELSEIF(ldown_option==3) THEN !Always use RH
          FCLD=WC_fraction(RH,Temp_C)
          EMIS_A=EMIS_CLOUD(EMIS_A,FCLD)
       ELSEIF(ldown_option==2) THEN ! FCLD obs, came from input
          EMIS_A=EMIS_CLOUD(EMIS_A,FCLD)
       ENDIF

       IF(ldown_option>1) THEN ! Forcing available if ldown_option=1, model otherwise
          LDOWN=EMIS_A*SIGMATK4
       ENDIF

       !----------------------------------------------------------------------------
       !Note that this is not averaged over the hour for cases where time step < 1hr
       KDOWN_HR=KDOWN
       IF (KDOWN_HR>0) THEN
          LUPCORR=(1-ALB0)*(0.08*KDOWN_HR)
       ELSE
          LUPCORR=0.
       ENDIF

       KUP=ALB0*KDOWN
       TSURF=((EMIS0*SIGMATK4+LUPCORR)/(EMIS0*SIGMA_SB))**0.25 !Eqs. (14) and (15),

       LUP=EMIS0*SIGMATK4+LUPCORR+(1-EMIS0)*LDOWN              !Eq (16) in Offerle et al. (2002)
       QSTAR=KDOWN-KUP+LDOWN-LUP
       TSURF=TSURF-273.16

       !======================================================================
       !Snow related parameters if snow pack existing
       IF (snowFrac(is)>0) THEN
          IF (AlbedoChoice==1.AND.180*ZENITH/ACOS(0.0)<90) THEN
             ALB1=SnowAlb+0.5e-16*(180*ZENITH/ACOS(0.0))**8 !AIDA 1982
          ELSE
             ALB1=SnowAlb
          ENDIF

          KUP_SNOW = (ALB1*(snowFrac(is)-snowFrac(is)*IceFrac(is))+ALB0*snowFrac(is)*IceFrac(is))*KDOWN
          TSURF_SNOW=((NARP_EMIS_SNOW*SIGMATK4)/(NARP_EMIS_SNOW*SIGMA_SB))**0.25 !Snow surface temperature

          !IF (TSURF_SNOW>273.16) TSURF_SNOW=min(273.16,Temp_K)!Set this to 2 degrees (melted water on top)
          !open(34,file='TestingSnowFrac.txt',position='append')
          !write(34,*) dectime,is,alb1,alb0,snowFrac(is),IceFrac(is),KDOWN,KUP_snow
          !close(34)

          LUP_SNOW = NARP_EMIS_SNOW*SIGMA_SB*TSURF_SNOW**4+(1-NARP_EMIS_SNOW)*LDOWN
          QSTAR_SNOW = KDOWN-KUP_SNOW+LDOWN-LUP_SNOW
          TSURF_SNOW = TSURF_SNOW-273.16

       ELSE
          KUP_SNOW = 0
          LUP_SNOW = 0
          TSURF_SNOW = 0
          QSTAR_SNOW = 0
          !QSTAR_ICE = 0
          !KUP_ICE = 0
       ENDIF

       qn1_ind_nosnow(is)=QSTAR          !Define sub-surface radiation components
       kup_ind_nosnow(is)=KUP
       lup_ind_nosnow(is)=LUP
       Tsurf_ind_nosnow(is)=TSURF

       qn1_ind_snow(is)=QSTAR_SNOW        !Define snow sub-surface radiation components
       kup_ind_snow(is)=KUP_SNOW
       lup_ind_snow(is)=LUP_SNOW
       Tsurf_ind_snow(is)=TSURF_SNOW


       IF (SF_all/=0)THEN
          QSTAR_SF = QSTAR_SF + QSTAR*sfr(is)*(1-snowFrac(is))/SF_all
       ELSE
          QSTAR_SF = QSTAR_SF + QSTAR*sfr(is)*(1-snowFrac(is))
       ENDIF

       IF ((1-SF_all)/=0)THEN
          QSTAR_S = QSTAR_S + QSTAR_SNOW*sfr(is)*snowFrac(is)/(1-SF_all)
       ELSE
          QSTAR_S = QSTAR_S + QSTAR_SNOW*sfr(is)*snowFrac(is)
       ENDIF

       !---------------------------------------------------------------------
       !Calculate weighted variables for each subsurface
       qn1_is = QSTAR*(1-snowFrac(is))+QSTAR_SNOW*snowFrac(is)
       kup_is = KUP*(1-snowFrac(is))+KUP_SNOW*snowFrac(is)
       lup_is = LUP*(1-snowFrac(is))+LUP_SNOW*snowFrac(is)
       tsurf_is = TSURF*(1-snowFrac(is))+TSURF_SNOW*snowFrac(is)

       qn1_cum=qn1_cum+(qn1_is*sfr(is))  !Calculate cumulative radiation components
       kup_cum=kup_cum+(kup_is*sfr(is))
       lup_cum=lup_cum+(lup_is*sfr(is))
       tsurf_cum=tsurf_cum+(tsurf_is*sfr(is))

       qn1_ind(is)=qn1_is                !Define sub-surface radiation components
       kup_ind(is)=kup_is
       lup_ind(is)=lup_is
       Tsurf_ind(is)=tsurf_is

    ENDDO !End of the surface types



    !Set overall radiation components
    IF (NetRadiationMethod/=3000) THEN !Observed Q* used and snow is modeled
       QSTARall=qn1_cum
    ELSE
       QSTARall=qn1_obs
    ENDIF

    KUPall=kup_cum
    LUPall=lup_cum
    TSURFall=tsurf_cum

    ! qn1=QSTARall
    ! kup=KUPall
    !LDOWN has same name
    ! lup=LUPall
    tsurf=TSURFall
    !if (kup>500) then
    ! write(*,*) Kdown, kup, kup_ind(1),kup_ind(2),kup_ind(3),kup_ind(4),kup_ind(5),kup_ind(6),SnowAlb
    ! pause
    !endif

    IF(DiagQN==1) WRITE(*,*) 'kdown: ',kdown,'kup:',kup,'ldown: ',ldown,'lup: ',lup

  END SUBROUTINE NARP

  !==============================================================================
  !OTHER REQUIRED

  !==============================================================================
  FUNCTION dewpoint(Temp_C,rh) RESULT(td)
    ! ea = vapor pressure (hPa)
    ! td = dewpoint (oC)
    ! calculates dewpoint in degC from
    ! http://www.atd.ucar.edu/weather_fl/dewpoint.html
    ! dewpoint = (237.3 * ln(e_vp/6.1078)) / (17.27 - (ln(e_vp/6.1078)))

    REAL(KIND(1d0))::rh,td,Temp_C,g
    !http://en.wikipedia.org/wiki/Dew_point
    g=((17.27*Temp_C)/(237.7+Temp_C))+LOG(rh/100)
    Td=(237.7*g)/(17.27-g)
    !td = (237.3 * LOG(ea_hPa/6.1078)) / (17.27 - (LOG(ea_hPa/6.1078)))
  END FUNCTION dewpoint
  !===============================================================================
  FUNCTION PRATA_EMIS(Temp_K,EA_hPa) RESULT(EMIS_A)
    ! clear sky emissivity function Prata 1996
    REAL(KIND(1d0))::Temp_K,ea_hPa,EMIS_A
    REAL(KIND(1d0))::W

    W=46.5*(ea_hPa/Temp_K)
    EMIS_A=1.-(1.+W)*EXP(-SQRT(1.2+3.*W))
  END FUNCTION PRATA_EMIS
  !===============================================================================
  FUNCTION EMIS_CLOUD(EMIS_A,FCLD) RESULT(em_adj)
    !calculates adjusted emissivity due to clouds
    REAL(KIND(1d0))::EMIS_A,FCLD,em_adj
    !T. Loridan, removing the square for FCLD in the emissivity correction
    !em_adj=EMIS_A+(1.-EMIS_A)*FCLD*FCLD
    em_adj=EMIS_A+(1.-EMIS_A)*FCLD
  END FUNCTION EMIS_CLOUD
  !===============================================================================
  FUNCTION EMIS_CLOUD_SQ(EMIS_A,FCLD) RESULT(em_adj)
    !calculates adjusted emissivity due to clouds
    REAL(KIND(1d0))::EMIS_A,FCLD,em_adj
    em_adj=EMIS_A+(1.-EMIS_A)*FCLD*FCLD
  END FUNCTION EMIS_CLOUD_SQ
  !===============================================================================
  FUNCTION cloud_fraction(KDOWN,KCLEAR) RESULT(FCLD)
    REAL(KIND(1d0))::KDOWN,KCLEAR,FCLD

    FCLD=1.-KDOWN/KCLEAR
    IF(FCLD>1.) FCLD=1.
    IF(FCLD<0.) FCLD=0.
  END FUNCTION cloud_fraction
  !===============================================================================
  FUNCTION WC_fraction(RH,Temp) RESULT(FWC)
    ! Thomas Loridan, King's College London: June 2009
    ! Parameterisation of fraction water content using the relative humidity
    REAL(KIND(1d0)),INTENT(in)   :: RH      !Relative Humidity in %
    REAL(KIND(1d0)),INTENT(in)   :: Temp    !Temperature in degre C

    REAL(KIND(1d0))              :: FWC     !Fraction water content between 0 and 1
    REAL(KIND(1d0))              :: A, B    !Parameters in the expo

    !Parameters
    !A=0.078
    !B=0.026

    A=0.185
    B=0.00019*Temp+0.015

    !FWC parameterization
    FWC=A * (EXP(B * RH)-1)
    IF(FWC>1.) FWC=1.
    IF(FWC<0.) FWC=0.
  END FUNCTION WC_fraction
  !===============================================================================
  !FUNCTION solar_zenith(lat,lng,timezone,dectime) RESULT(zenith)
  !  !Stull, 1989
  !  !returns zenith in radians
  !  !lat, lng in RADS
  !  REAL(KIND(1d0)) ::lat,lng,timezone,dectime,zenith, eqtime
  !  REAL(KIND(1d0)) ::ha, decl,tst, time_offset,gamma
  !
  !  gamma=2.*3.141592654/365.25463*dectime
  !  eqtime=229.18*(7.5e-5+1.868e-3*COS(gamma)-0.032077*SIN(gamma)&
  !       -0.014615*COS(2.*gamma)-0.040849*SIN(2.*gamma))
  !  decl=6.918e-3-0.399912*COS(gamma)+0.070257*SIN(gamma)&
  !       -0.006758*COS(2.*gamma)+9.07e-4*SIN(2.*gamma)-2.697e-3*COS(3.*gamma)&
  !       +1.48e-3*SIN(3.*gamma)
  !  time_offset=eqtime-4.*lng*RAD2DEG-60.*timezone
  !  tst=(dectime-FLOOR(dectime))*1440.+time_offset
  !  ha=(tst/4.)-180.
  !  ha=ha*DEG2RAD
  !  zenith=ACOS(SIN(lat)*SIN(decl)+COS(lat)*COS(decl)*COS(ha))
  !END FUNCTION solar_zenith

  !===============================================================================
  FUNCTION ISURFACE(doy,zenith) RESULT(Isurf)
    ! Calculates ground level solar irradiance clear sky
    ! assuming transmissivity = 1
    ! let it report zero if zenith >= 90
    REAL(KIND(1d0))::zenith,Isurf
    INTEGER::doy
    REAL(KIND(1d0))::Rmean, Rse, cosZ,Itoa
    REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292

    Rmean = 149.6                 !Stull 1998
    Rse=solar_ESdist(doy)
    IF(zenith<90.*DEG2RAD) THEN
       cosZ = COS(zenith)
       Itoa = 1370.*(Rmean/Rse)**2  !top of the atmosphere
       Isurf = Itoa*cosZ            !ground level solar irradiance in W/m2
    ELSE
       Isurf = 0.
    ENDIF

  END FUNCTION ISURFACE

  !===============================================================================
  FUNCTION solar_ESdist(doy) RESULT(Rse)
    !from Stull, 1998   Keep! called from SOLWEIG_clearnessindex_2013b
    INTEGER          ::doy
    REAL(KIND(1d0))             ::Rse
    REAL(KIND(1d0)) ::MA,nu,e,a

    e = 0.0167
    a = 146.457

    MA = 2.*3.141592654*(doy-3)/365.25463 !Mean anomaly
    nu=MA+0.0333988*SIN(MA)+.0003486*SIN(2.*MA)+5e-6*SIN(3.*MA) !true anomaly
    Rse = a*(1-e*e)/(1+e*COS(nu))

  END FUNCTION solar_ESdist

  !===============================================================================
  FUNCTION SmithLambda(lat) RESULT(G)
    USE FileName
    USE defaultnotUsed
    !read kriged data based on Smith 1966 (JAM)
    INTEGER :: lat,ios,ilat
    REAL(KIND(1d0)),DIMENSION(365):: G

    !open(99,file="Smith1966.grd",access="direct",action="read",recl=365*4,iostat=ios)
    !read(99,rec=lat+1,iostat=ios) G
    OPEN(99,file=smithFile,iostat=ios)
    DO ilat=1,lat
       READ(99,*)
    ENDDO
    READ(99,*,iostat=ios)ilat, G
    IF (ios/=0) THEN
       CALL  ErrorHint(11,'reading Smith1966.grd (ios).',notUsed,notUsed,ios)
    ENDIF
    CLOSE(99)
  END FUNCTION SmithLambda

  !===============================================================================
  FUNCTION transmissivity(Press_hPa,Temp_C_dew,G,zenith) RESULT(trans)
    ! bulk atmospheric transmissivity (Crawford and Duchon, 1999)
    ! P = pressure (hPa)
    ! Td = dewpoint (C)
    ! G parameter is empirical value from Smith 1966 (JAM)
    ! zenith in radians
    ! if zenith > 80 use the value for 80.

    REAL(KIND(1d0)) ::Press_hPa,TemP_C_dew,zenith,G,trans
    REAL(KIND(1d0))::m,TrTpg,u,Tw,Ta,cosZ
    REAL(KIND(1d0))::Tdf
    REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292


    IF (zenith>80.*DEG2RAD) THEN
       cosZ=COS(80.*DEG2RAD)
    ELSE
       cosZ=COS(zenith)
    ENDIF

    Tdf = TemP_C_dew*1.8+32. !celsius to fahrenheit
    !	Transmission coefficients
    m = 35*cosZ/SQRT(1224.*cosZ*cosZ+1)            !optical air mass at p=1013 mb
    !Rayleigh & permanent gases
    TrTpg = 1.021-0.084*SQRT(m*(0.000949*Press_hPa+0.051)) !first two trans coeff
    u = EXP(0.113-LOG(G+1)+0.0393*Tdf)             !precipitable water
    Tw = 1-0.077*(u*m)**0.3             !vapor transmission coe3ff.
    Ta = 0.935**m                       !4th trans coeff
    trans = TrTpg*Tw*Ta                 !bulk atmospheric transmissivity
  END FUNCTION transmissivity
  !===============================================================================
END MODULE NARP_MODULE
