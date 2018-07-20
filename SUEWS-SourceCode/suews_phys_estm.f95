!==============================================================================
MODULE mod_interp
  !     Created on Thu Jan 22 00:06:32 2004
  IMPLICIT NONE

CONTAINS
  ELEMENTAL FUNCTION interp1d(x1,x2,y1,y2,xi) RESULT(yi)
    REAL(KIND(1d0)),INTENT(in) ::x1,x2,xi
    REAL(KIND(1d0)),INTENT(in) ::y1,y2
    REAL (KIND(1D0))::b0,b1
    REAL(KIND(1d0))         ::yi
    !integer         ::ny                 !!!!!FO!!!!!
    b1=(y2-y1)/(x2-x1)
    b0=y1-b1*x1
    yi=b0+b1*xi
  END FUNCTION interp1d
END MODULE mod_interp

!=============================================================================
MODULE mod_solver
  !     Created on Thu Jan 22 07:50:01 2004
  !     Copyright (c) 2001 MyCompany. All rights reserved.

  IMPLICIT NONE

CONTAINS

  FUNCTION NewtonPolynomial(x0,Pcoeff,conv,maxiter) RESULT(x)
    !Solves Newton's Method for a polynomial of the form
    !f(x)=Pcoeff(1)*x^n+Pcoeff(2)*x^n-1+...+Pcoeff(n+1)
    !                f(x(i))
    ! x(i+1) = x(i)- -------
    !                f'(x(i))
    !
    !conv is the level required for convergence
    !maxiter is the maximum allowed iterations
    !----------------------------------------------------
    REAL(KIND(1d0)) ::x0,x,conv
    REAL(KIND(1d0)) ::Pcoeff(:)
    REAL(KIND(1d0)) ::e, xprev
    REAL (KIND(1D0))::f,fp
    INTEGER ::maxiter
    INTEGER ::niter
    LOGICAL ::converged=.FALSE.
    INTEGER ::n,i,j

    e=HUGE(1.)
    n=SIZE(Pcoeff)
    x=x0
    DO i=1,maxiter
       IF (ABS(e)<conv) THEN
          converged=.TRUE.
          EXIT
       ENDIF
       f=0; fp = 0
       DO j=1,n-1
          f = f+Pcoeff(j)*x**(n-j)
          fp = fp + Pcoeff(j)*(n-j)*x**(n-j-1)                              !!FO!! derivative
       ENDDO

       f = f + Pcoeff(n)
       xprev = x
       IF (fp==0.) fp=TINY(1.)
       x = xprev - f/fp
       e = x-xprev
    ENDDO
    niter=i-1
    IF (.NOT.converged) THEN
       PRINT*, "Solution did not converge. Niter=", niter, " Error=", e
       x=x0
    ENDIF
  END FUNCTION NewtonPolynomial
END MODULE mod_solver
!==============================================================================

!==============================================================================
MODULE modSolarCalc
  USE MathConstants
  IMPLICIT NONE

CONTAINS
  !=======================================================
  FUNCTION min_zenith(lat,doy) RESULT(zmin)
    !returns max zenith
    !returns zenith in radians for lat, lng in degrees
    REAL(KIND(1d0)) ::lat,dectime,zmin
    REAL(KIND(1d0)) ::latr,decl
    INTEGER :: doy
    dectime=float(doy)
    latr=lat*dtr
    decl=0.409*COS(2*pi*(dectime-173)/365.25)
    zmin=pi/2.-ASIN(SIN(latr)*SIN(decl)-COS(latr)*COS(decl)*(-1))
  END FUNCTION min_zenith
  !=======================================================

  !=======================================================
  FUNCTION Local_apparent_time(lng,dectime) RESULT(la_time)
    !Oke, 1989, equation of time elsewhere
    REAL(KIND(1d0)) ::lng,dectime,la_time
    REAL(KIND(1d0)) ::gamma,eqtime,lmst

    lmst=dectime-4.*lng/60./1440.
    gamma=2.*pi/365.*(lmst-1.)
    eqtime=229.18*(7.5e-5+1.868e-3*COS(gamma)-0.032077*SIN(gamma)&
         &    -0.014615*COS(2.*gamma)-0.040849*SIN(2.*gamma))
    la_time=lmst+eqtime/1440.
  END FUNCTION Local_apparent_time
  !=======================================================

  SUBROUTINE Solar_angles(lat,lng,timezone,dectime,decl,zenith,azimuth)

    REAL, INTENT(in)  ::lat,lng,timezone,dectime
    INTEGER                 ::doy,hour,mn
    REAL(KIND(1d0)), INTENT(out)  ::decl,zenith,azimuth
    REAL (KIND(1d0))  ::ha,latr,eqtime,tst,&
         time_offset,gamma           !!!!!FO!!!!! lngr, phi, theta

    latr=lat*pi/180.
    doy=FLOOR(dectime)
    hour=FLOOR((dectime-doy)*24.)
    mn=FLOOR((dectime-doy-hour/24.)*60.)   !!!Check this

    gamma=2.*pi/365.25463*(doy-1.+(hour-12.)/24.)
    eqtime=229.18*(7.5e-5+1.868e-3*COS(gamma)-0.032077*SIN(gamma)&
         &    -0.014615*COS(2.*gamma)-0.040849*SIN(2.*gamma))
    decl=6.918e-3-0.399912*COS(gamma)+0.070257*SIN(gamma)&
         &    -0.006758*COS(2.*gamma)+9.07e-4*SIN(2.*gamma)-2.697e-3*COS(3.*gamma)&
         &    +1.48e-3*SIN(3.*gamma)
    time_offset=eqtime-4.*lng+60.*timezone
    tst=hour*60.+mn+time_offset
    ha=(tst/4.)-180.
    ha=ha*pi/180.

    zenith=ACOS(SIN(latr)*SIN(decl)+COS(latr)*COS(decl)*COS(ha))
    azimuth=pi+ACOS((SIN(latr)*COS(zenith)-SIN(decl))/(COS(latr)*SIN(zenith)))

    RETURN
  END SUBROUTINE Solar_angles

  !=========================== SolarTimes ==================================
  SUBROUTINE Solar_Times(lat,lng,timezone,dectime,sunrise,sunset,snoon)
    !  for sunrise and sunset ha = ha(zenith=90)
    !  timezone is offset to GMT e.g. -5 for EST

    REAL(KIND(1d0)), INTENT(in)  ::lat,lng,timezone,dectime
    INTEGER                 ::doy
    REAL(KIND(1d0)), INTENT(out)  ::sunrise, sunset, snoon
    REAL(KIND(1d0))  :: ha, latr, eqtime, gamma, zenith, decl
    latr = lat*dtr
    zenith=90.833*dtr
    doy=FLOOR(dectime)
    gamma=2.*pi/365.*(float(doy)-0.5) !fractional year
    eqtime=229.18*(7.5e-5+1.868e-3*COS(gamma)-0.032077*SIN(gamma)&
         &    -0.014615*COS(2.*gamma)-0.040849*SIN(2.*gamma))
    decl=6.918e-3-0.399912*COS(gamma)+0.070257*SIN(gamma)&
         &    -0.006758*COS(2.*gamma)+9.07e-4*SIN(2.*gamma)-2.697e-3*COS(3.*gamma)&
         &    +1.48e-3*SIN(3.*gamma)
    ha=ACOS(COS(zenith)/(COS(latr)*COS(decl))-TAN(latr)*TAN(decl))
    ha=ha*rtd
    sunrise = (720.-4.*(lng-ha)-eqtime)/60.-timezone
    sunset = (720.-4.*(lng+ha)-eqtime)/60.-timezone
    snoon = (720.-4.*lng-eqtime)/60.-timezone
    RETURN
  END SUBROUTINE Solar_Times
  !=======================================================
  FUNCTION kdown_surface(doy,zenith) RESULT(Isurf)
    ! Calculates ground level solar irradiance clear sky
    ! assuming transmissivity = 1
    ! let it report zero if zenith >= 90
    REAL(KIND(1d0))    ::zenith,Isurf
    INTEGER    ::doy
    REAL (KIND(1d0))::Rmean, Rse, cosZ,Itoa

    Rmean = 149.6   !Stull 1998
    Rse=solar_ESdist(doy)
    IF(zenith<pi/2.) THEN
       cosZ = COS(zenith)
       Itoa = 1370*(Rmean/Rse)**2    !top of the atmosphere
       Isurf = Itoa*cosZ      !ground level solar irradiance in W/m2
    ELSE
       Isurf = 0.
    ENDIF

  END FUNCTION kdown_surface

  !=======================================================
  FUNCTION SmithLambda(lat) RESULT(G)
    !read kriged data based on Smith 1966 (JAM)
    INTEGER :: lat,ios
    REAL,DIMENSION(365):: G

    OPEN(99,file="Smith1966.grd",access="direct",action="read",recl=365*4,iostat=ios)
    IF (ios/=0) THEN
       PRINT*, "Iostat=",ios," reading Smith1966.grd"
       STOP
    ENDIF
    READ(99,rec=lat+1,iostat=ios) G
    IF (ios/=0) PRINT*, "Iostat=", ios, " reading Smith1966.grd"
    CLOSE(99)
  END FUNCTION SmithLambda
  !=======================================================
  FUNCTION transmissivity_CD(P,Td,G,zenith) RESULT(trans)           !!!!!FO!!!!! ,doy
    ! bulk atmospheric transmissivity (Crawford and Duchon, 1999)
    ! P = pressure (hPa)
    ! Td = dewpoint (C)
    ! G parameter is empirical value from Smith 1966 (JAM)
    ! zenith in radians


    !        integer         ::doy           !!!!!FO!!!!!
    REAL(KIND(1d0))    ::P,Td,zenith,G,trans
    REAL (KIND(1d0))::m,TrTpg,u,Tw,Ta,cosZ
    REAL (KIND(1d0))::Tdf

    IF (zenith>80.*dtr) THEN
       cosZ=COS(80.*dtr)
    ELSE
       cosZ=COS(zenith)
    ENDIF
    Tdf = Td*1.8+32. !celsius to fahrenheit
    !  Transmission coefficients
    m = 35*cosZ/SQRT(1224.*cosZ*cosZ+1) !optical air mass at p=1013 mb
    TrTpg = 1.021-0.084*SQRT(m*(0.000949*P+0.051)) !first two trans coeff
    u = EXP(0.113-LOG(G+1)+0.0393*Tdf) !precipitable water
    Tw = 1-0.077*(u*m)**0.3    !vapor transmission coe3ff.
    Ta = 0.935**m        !4th trans coeff
    trans = TrTpg*Tw*Ta              !bulk atmospherics transmissivity
  END FUNCTION transmissivity_CD

  !   !=======================================================
  !   ! NB:this FUNCTION is problematic:
  !   ! these variables are NEVER initialized/calculated: tr,tg,tw,ta
  !   ! but used for calcuLAItng Sdir!
  !   FUNCTION kdown_niemala(S0,vap_press,Tk) RESULT(kdown)           !!!!!FO!!!!! ,albedo
  !     !from Niemala et al. 2001. Atmospheric Research 58:141-154
  !     ! using generalized formula only, no empirical data from site
  !     !11 October 2001
  !     ! S0=base solar insolation at the surface
  !     ! albedo = surface albedo (required for multireflection diffuse)
  !     ! vap_press=screen level vapor pressure (Pa)
  !     ! Tk=screen level air temperature (K)
  !     ! ta=aerosol transmissivity
  !     ! tR=Rayleigh scattering transmissivity
  !     ! tg=uniformly mixed gas transmissivity
  !     ! tw=water vapor transmissivity
  !     ! toz=ozone transmissivity
  ! !!!!!!!INCOMPLETE
  !     REAL  ::S0,vap_press,Tk,kdown           !!!!!FO!!!!! ,albedo
  !     REAL  ::Sdir,Diffuse,Diffuse_R,Diffuse_a,Diffuse_m
  !     REAL  ::tr,tg,tw,ta,toz, theta
  !     CALL transmissivity(vap_press,Tk,theta,tr,tg,tw,ta,toz)
  !     Sdir=0.9751*S0*tr*tg*tw*ta*toz
  !     Diffuse_R=0.
  !     Diffuse_a=0.
  !     Diffuse_m=0.
  !     Diffuse=Diffuse_R+Diffuse_a+Diffuse_m
  !     kdown=Sdir+Diffuse
  !   END FUNCTION kdown_niemala
  !   !=======================================================
  !   SUBROUTINE transmissivity(vap_press,Tk,theta,tr,tg,tw,ta,toz)
  !     !calculates atmospheric transmissivities for
  !     ! ta=aerosol transmissivity
  !     ! tR=Rayleigh scattering transmissivity
  !     ! tg=uniformly mixed gas transmissivity
  !     ! tw=water vapor transmissivity
  !     ! toz=ozone transmissivity
  !     ! from Iqbal (1983) and Niemala(2001)
  !     ! vap_press (Pa)
  !     ! Tk air temp (K)
  !     ! w=precipitable water content (cm)
  !     !==========INCOMPLETE
  !     REAL  ::Tk, vap_press,theta
  !     REAL  ::tr,tg,tw,ta,toz,w
  !     w=0.493*vap_press/Tk
  !     ta=0.59+0.012*theta-1.336e-4*theta*theta
  !
  !   END SUBROUTINE transmissivity

  !=======================================================
  FUNCTION solar_ESdist(doy) RESULT(Rse)
    !from Stull, 1998
    INTEGER  ::doy
    REAL(KIND(1d0))  ::Rse
    REAL (KIND(1d0)) ::MA,nu,e,a

    e = 0.0167
    a = 146.457

    MA = 2.*pi*(doy-3)/365.25463 !Mean anomaly
    nu=MA+0.0333988*SIN(MA)+.0003486*SIN(2.*MA)+5e-6*SIN(3.*MA) !true anomaly
    Rse = a*(1-e*e)/(1+e*COS(nu))

  END FUNCTION solar_ESdist

END MODULE modSolarCalc

!=====================================================================================
!=====================================================================================
MODULE heatflux
  IMPLICIT NONE
CONTAINS

  SUBROUTINE heatcond1d(T,Qs,dx,dt,k,rhocp,bc,bctype)
    REAL(KIND(1d0)),INTENT(inout)::T(:)
    REAL(KIND(1d0)),INTENT(in)::dx(:),dt,k(:),rhocp(:),bc(2)
    REAL(KIND(1d0)),INTENT(out)::Qs
    LOGICAL,INTENT(in)::bctype(2)
    INTEGER         ::i,n!,j       !!!!!FO!!!!!
    REAL(KIND(1d0)),ALLOCATABLE::w(:),a(:),T1(:)
    n=SIZE(T)
    ALLOCATE(w(0:n),a(n),T1(n))
    !w = interface tempea
    w(1:n)=T
    w(0)=bc(1); w(n)=bc(2)
    !convert from flux to equivalent temperature, not exact
    ! F = k dT/dX => dx*F/k + Ti = Ti
    IF (bctype(1)) w(0)=bc(1)*0.5*dx(1)/k(1)+w(1)
    IF (bctype(2)) w(n)=bc(2)*0.5*dx(n)/k(n)+w(n)

    a=k/dx
    DO i=1,n-1
       w(i)=(T(i+1)*a(i+1)+T(i)*a(i))/(a(i)+a(i+1))
    ENDDO
    !!FO!! print*, 'w: ', w
    DO i=1,n
       T1(i) = (dt/rhocp(i))*(w(i-1)-2*T(i) + w(i))*2*a(i)/dx(i) + T(i)
    ENDDO
    !!FO!! print*, 'T1: ', T1
    !for storage the internal distribution of heat should not be important
    Qs = (w(0)-T(1))*2*a(1) + (w(n)-T(n))*2*a(n)                           !!FO!! k*d(dT/dx)/dx = rhoCp*(dT/dt) -- rhoCp*(dT/dt)*dx = dQs -- dQs = k*d(dT/dx)
    ! Qs=sum((T1-T)*rhocp*dx)/dt!
    T=T1
  END SUBROUTINE heatcond1d
END MODULE heatflux

!==================================================================================


MODULE ESTM_module
  !===============================================================================
  ! revision history:
  ! TS 09 Oct 2017: re-organised ESTM subroutines into a module
  !===============================================================================
  IMPLICIT NONE


CONTAINS

  !======================================================================================
  ! Subroutine to read in ESTM data in the same way as met data (SUEWS_InitializeMetData)
  ! HCW 30 Jun 2016
  SUBROUTINE SUEWS_GetESTMData(lunit)
    USE allocateArray, ONLY:ncolsestmdata, estmforcingdata
    USE data_in, ONLY: fileestmts, skipheadermet
    USE sues_data, ONLY: tstep_real, tstep
    USE defaultnotUsed, ONLY: notused, ios_out
    USE Initial, ONLY: skippedlines, readlinesmetdata, gridcounter

    IMPLICIT NONE

    INTEGER,INTENT(in)::lunit
    INTEGER::i,iyy !,RunNumber,NSHcounter
    INTEGER :: iostat_var
    REAL(KIND(1d0)),DIMENSION(ncolsESTMdata):: ESTMArray
    REAL(KIND(1d0)):: imin_prev, ih_prev, iday_prev, tstep_estm   !For checks on temporal resolution of estm data

    !---------------------------------------------------------------

    !Open the file for reading and read the actual data
    !write(*,*) FileESTMTs
    OPEN(lunit,file=TRIM(FileESTMTs),status='old',err=315)
    CALL skipHeader(lunit,SkipHeaderMet)

    ! Skip to the right place in the ESTM file, depending on how many chunks have been read already
    IF (skippedLines>0) THEN
       DO iyy=1,skippedLines
          READ(lunit,*)
       ENDDO
    ENDIF

    ! Read in next chunk of ESTM data and fill ESTMForcingData array with data for every timestep
    DO i=1,ReadlinesMetdata
       READ(lunit,*,iostat=iostat_var) ESTMArray
       ESTMForcingData(i,1:ncolsESTMdata,GridCounter) = ESTMArray
       ! Check timestamp of met data file matches TSTEP specified in RunControl
       IF(i==1) THEN
          imin_prev = ESTMArray(4)
          ih_prev   = ESTMArray(3)
          iday_prev = ESTMArray(2)
       ELSEIF(i==2) THEN
          tstep_estm = ((ESTMArray(4)+60*ESTMArray(3)) - (imin_prev+60*ih_prev))*60   !tstep in seconds
          IF(tstep_estm/=tstep_real.AND.ESTMArray(2)==iday_prev) THEN
             CALL ErrorHint(39,'TSTEP in RunControl does not match TSTEP of ESTM data (DOY).',REAL(tstep,KIND(1d0)),tstep_estm,&
                  INT(ESTMArray(2)))
          ENDIF
       ENDIF
    ENDDO

    CLOSE(lunit)

    RETURN

315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,ios_out)

  END SUBROUTINE SUEWS_GetESTMData
  !======================================================================================


  !======================================================================================
  SUBROUTINE ESTM_initials

    ! Last modified HCW 30 Jun 2016 - reading in now done by SUEWS_GetESTMData subroutine.
    !                                 ESTM_initials now only runs once per run at the very start.
    ! Last modified HCW 15 Jun 2016 - code now reads in 5-min file (interpolation done beforehand, outside of SUEWS itself)

    USE defaultNotUsed
    USE PhysConstants, ONLY: c2k
    USE ESTM_data
    USE allocateArray
    USE data_in, ONLY: fileinputpath
    USE Initial, ONLY: numberofgrids

    IMPLICIT NONE

    !=====Read ESTMinput.nml================================
    NAMELIST/ESTMinput/TsurfChoice,&
         evolveTibld,              &
         ibldCHmod,                &
         LBC_soil,                 &
         THEAT_ON,                 &
         THEAT_OFF,                &
         THEAT_fix

    OPEN(511,file=TRIM(FileInputPath)//'ESTMinput.nml',status='old')
    READ(511,nml=ESTMinput)
    CLOSE(511)

    !Convert specified temperatures to Kelvin
    THEAT_ON=THEAT_ON+C2K
    THEAT_OFF=THEAT_OFF+C2K
    THEAT_fix=THEAT_fix+C2K

    ALLOCATE(Tair2_grids(NumberOfGrids))
    ALLOCATE(lup_ground_grids(NumberOfGrids))
    ALLOCATE(lup_wall_grids(NumberOfGrids))
    ALLOCATE(lup_roof_grids(NumberOfGrids))
    ALLOCATE(Tievolve_grids(NumberOfGrids))
    ALLOCATE(T0_ibld_grids(NumberOfGrids))
    ALLOCATE(T0_ground_grids(NumberOfGrids))
    ALLOCATE(T0_wall_grids(NumberOfGrids))
    ALLOCATE(T0_roof_grids(NumberOfGrids))
    ALLOCATE(TN_wall_grids(NumberOfGrids))
    ALLOCATE(TN_roof_grids(NumberOfGrids))

  END SUBROUTINE ESTM_initials
  !======================================================================================


  SUBROUTINE ESTM_translate(Gridiv)
    ! HCW 30 Jun 2016

    USE defaultNotUsed,ONLY: nan
    USE PhysConstants, ONLY: c2k, sbconst
    USE ESTM_data
    USE allocateArray
    USE gis_data, ONLY: bldgh
    USE Initial, ONLY: numberofgrids

    IMPLICIT NONE
    INTEGER :: i
    !REAL(KIND(1d0)) :: CFLval
    !REAL(KIND(1d0)) :: t5min
    REAL(KIND(1d0))::W,WB
    !CHARACTER (len=20)::FileCodeX
    !CHARACTER (len=150):: FileFinalTemp
    !LOGICAL:: inittemps=.FALSE.
    INTEGER:: ESTMStart=0
    INTEGER:: Gridiv

    !Set initial values at the start of each run for each grid
    IF(Gridiv == 1) ESTMStart = ESTMStart+1
    IF(ESTMStart==1) THEN

       !write(*,*) ' ESTMStart: ',ESTMStart, 'initialising ESTM for grid no. ', Gridiv

       TFLOOR=20.0 ! This is used only when radforce =T  !TODO:  should be put in the namelist
       TFLOOR=TFLOOR+C2K

       ! Initial values
       Tievolve=20.0 + C2K
       SHC_air=1230.0
       minshc_airbld=1300

       ! ---- Internal view factors ----
       !constant now but should be calculated in the future
       IVF_IW =   0.100000
       IVF_IR =   0.000000
       IVF_II =   0.900000
       IVF_IF =   0.000000
       IVF_WW =   0.050000
       IVF_WR =   0.000000
       IVF_WI =   0.950000
       IVF_WF =   0.000000
       IVF_RW =   0.050000
       IVF_RI =   0.950000
       IVF_RF =   0.000000
       IVF_FW =   0.050000
       IVF_FR =   0.000000
       IVF_FI =   0.950000

       Tair24HR=C2K

       !Ts5mindata(1,ncolsESTMdata) = -999
       ! !Fill Ts5mindata for current grid and met block - this is done in SUEWS_translate
       Ts5mindata(1,1:ncolsESTMdata) = ESTMForcingData(1,1:ncolsESTMdata,Gridiv)


       ! ---- Initialization of variables and parameters for first row of run for each grid ----
       ! N layers are calculated in SUEWS_translate
       IF ( .NOT. ALLOCATED(Tibld) ) THEN
          ! print*, "Nibld",Nibld
          ! print*, "Nwall",Nwall
          ! print*, "Nroof",Nroof
          ! print*, "Nground",Nground
          ALLOCATE(Tibld(Nibld),Twall(Nwall),Troof(Nroof),Tground(Nground),Tw_4(Nwall,4))
          ALLOCATE(Tibld_grids(Nibld,NumberOfGrids), &
               Twall_grids(Nwall,NumberOfGrids), &
               Troof_grids(Nroof,NumberOfGrids), &
               Tground_grids(Nground,NumberOfGrids), &
               Tw_4_grids(Nwall,4,NumberOfGrids))
       ENDIF

       ! Transfer variables from Ts5mindata to variable names
       ! N.B. column numbers here for the following file format - need to change if input columns change!
       ! dectime iy id it imin Tiair Tsurf Troof Troad Twall Twall_n Twall_e Twall_s Twall_w
       !        1  2  3  4    5     6     7     8     9     10      11      12      13       !new
       !
       ! Calculate temperature of each layer in Kelvin
       !  QUESTION: what if (Nground/Nwall/Nroof-1)==0? TS 21 Oct 2017
       DO i=1,Nground
          Tground(i)=(LBC_soil-Ts5mindata(1,cTs_Troad))*(i-1)/(Nground-1)+Ts5mindata(1,cTs_Troad)+C2K
       ENDDO
       DO i=1,Nwall
          Twall(i)=(Ts5mindata(1,cTs_Tiair)-Ts5mindata(1,cTs_Twall))*(i-1)/(Nwall-1)+Ts5mindata(1,cTs_Twall)+C2K
       ENDDO
       DO i=1,Nroof
          Troof(i)=(Ts5mindata(1,cTs_Tiair)-Ts5mindata(1,cTs_Troof))*(i-1)/(Nroof-1)+Ts5mindata(1,cTs_Troof)+C2K
       ENDDO
       Tibld(1:Nibld)=Ts5mindata(1,cTs_Tiair)+C2K

    ENDIF  !End of loop run only at start (for each grid)

    ! ---- Parameters related to land surface characteristics ----
    ! QUESTION: Would Zref=z be more appropriate?
    ZREF=2.0*BldgH      !!FO!! BldgH: mean bulding hight, zref: local scale reference height (local: ~ 10^2 x 10^2 -- 10^3 x 10^3 m^2)

    svf_ground=1.0
    svf_roof=1.0

    ! ==== roof (i.e. Bldgs)
    !froof=sfr(BldgSurf)   ! Moved to SUEWS_translate HCW 16 Jun 2016
    alb_roof=alb(BldgSurf)
    em_roof=emis(BldgSurf)

    ! ==== vegetation (i.e. EveTr, DecTr, Grass)
    !fveg=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)  ! Moved to SUEWS_translate HCW 16 Jun 2016
    IF(fveg/=0) THEN
       alb_veg=(alb(ConifSurf)*sfr(ConifSurf) + alb(DecidSurf)*sfr(DecidSurf) + alb(GrassSurf)*sfr(GrassSurf))/fveg
       em_veg=(emis(ConifSurf)*sfr(ConifSurf) + emis(DecidSurf)*sfr(DecidSurf) + emis(GrassSurf)*sfr(GrassSurf))/fveg
    ELSE ! check fveg==0 scenario to avoid division-by-zero error, TS 21 Oct 2017
       alb_veg=NAN
       em_veg=NAN
    ENDIF

    ! ==== ground (i.e. Paved, EveTr, DecTr, Grass, BSoil, Water - all except Bldgs)
    !fground=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(PavSurf)+sfr(BsoilSurf)+sfr(WaterSurf) ! Moved to SUEWS_translate HCW 16 Jun 2016
    IF(fground/=0) THEN
       alb_ground=(alb(ConifSurf)*sfr(ConifSurf)+alb(DecidSurf)*sfr(DecidSurf)&
            +alb(GrassSurf)*sfr(GrassSurf)+alb(PavSurf)*sfr(PavSurf)&
            +alb(BsoilSurf)*sfr(BsoilSurf)+alb(WaterSurf)*sfr(WaterSurf))/fground
       em_ground=(emis(ConifSurf)*sfr(ConifSurf)+emis(DecidSurf)*sfr(DecidSurf)&
            +emis(GrassSurf)*sfr(GrassSurf)+emis(PavSurf)*sfr(PavSurf)&
            +emis(BsoilSurf)*sfr(BsoilSurf)+emis(WaterSurf)*sfr(WaterSurf))/fground
    ELSE ! check fground==0 scenario to avoid division-by-zero error, TS 21 Jul 2016
       alb_ground=NAN
       em_ground=NAN
    ENDIF

    IF(froof<1.0) THEN
       HW=fwall/(2.0*(1.0-froof))
    ELSE
       !  HW=0  !HCW if only roof, no ground
       HW=0.00001  ! to avoid zero-HW scenario TS 21 Oct 2017

    END IF
    HW=MAX(0.00001,HW)! to avoid zero-HW scenario TS 27 Oct 2017

    IF (Fground==1.0) THEN   !!FO!! if only ground, i.e. no houses
       W=1
       WB=0
       SVF_ground=1.
       zvf_WALL=0.
       SVF_WALL=0.
       SVF_ROOF=1.
       zvf_ground=0.
       xvf_wall=0.
       RVF_CANYON=1.
       RVF_ground=1.-FVEG
       RVF_ROOF=0
       RVF_WALL=0
       RVF_VEG=FVEG
    ELSE IF ( Fground==0.0 ) THEN !check fground==0 (or HW==0) scenario to avoid division-by-zero error, TS 21 Jul 2016
       ! the following values are calculated given HW=0
       W=0
       WB=1
       zvf_WALL= 0 !COS(ATAN(2/HW))  when HW=0                                 !!FO!! wall view factor for wall
       HW=0
       SVF_ground=MAX(COS(ATAN(2*HW)),0.00001)!!FO!! sky view factor for ground ! to avoid zero-division scenario TS 21 Oct 2017
       SVF_WALL=(1-zvf_WALL)/2                                                 !!FO!! sky view factor for wall
       zvf_ground=1-svf_ground                                                 !!FO!! wall view factor for ground
       xvf_wall=svf_wall                                                       !!FO!! ground view factor
       !   RVF_CANYON=COS(ATAN(2*ZREF/W))
       !   RVF_ROOF=1-RVF_CANYON
       !   RVF_WALL=(COS(ATAN(2*(ZREF-BldgH)/W))-RVF_CANYON)*RVF_CANYON
       !   RVF_ground=RVF_CANYON-RVF_WALL
       RVF_ground=(fground-fveg)*SVF_ground
       RVF_veg=fveg*SVF_ground
       RVF_ROOF=froof
       RVF_Wall=1-RVF_ROOF-RVF_ground-RVF_VEG
    ELSE
       W=BldgH/HW   !What about if HW = 0, need to add IF(Fground ==0) option ! fixed by setting a small number, TS 21 Oct 2017
       WB=W*SQRT(FROOF/Fground)
       SVF_ground=COS(ATAN(2*HW))                                              !!FO!! sky view factor for ground
       zvf_WALL=COS(ATAN(2/HW))                                                !!FO!! wall view factor for wall
       SVF_WALL=(1-zvf_WALL)/2                                                 !!FO!! sky view factor for wall
       zvf_ground=1-svf_ground                                                 !!FO!! wall view factor for ground
       xvf_wall=svf_wall                                                       !!FO!! ground view factor
       !   RVF_CANYON=COS(ATAN(2*ZREF/W))
       !   RVF_ROOF=1-RVF_CANYON
       !   RVF_WALL=(COS(ATAN(2*(ZREF-BldgH)/W))-RVF_CANYON)*RVF_CANYON
       !   RVF_ground=RVF_CANYON-RVF_WALL
       RVF_ground=(fground-fveg)*SVF_ground
       RVF_veg=fveg*SVF_ground
       RVF_ROOF=froof
       RVF_Wall=1-RVF_ROOF-RVF_ground-RVF_VEG
    ENDIF

    alb_avg=alb_ground*RVF_ground + alb_wall*RVF_WALL + alb_roof*RVF_ROOF + alb_veg*RVF_VEG

    sumalb=0.; nalb=0
    sumemis=0.; nemis=0

    !set emissivity for ceiling, wall and floor inside of buildings
    em_r = em_ibld; em_w=em_ibld; em_i=em_ibld; em_f=em_ibld

    !internal elements
    IF (nroom==0) THEN
       fibld = (FLOOR(BldgH/3.1-0.5)-1)*froof
    ELSE
       fibld = (2.-2./nroom)*fwall + (FLOOR(BldgH/3.1-0.5)-1)*froof
    ENDIF

    IF (fibld==0) fibld=0.00001 !this just ensures a solution to radiation
    finternal = froof+fibld+fwall
    IF (finternal==0) finternal=0.00001 ! to avoid zero-devision error TS 21 Oct 2017
    fair=zref-BldgH*froof
    IF (fair==0) fair=0.00001 ! to avoid zero-devision error TS 21 Oct 2017
    !ivf_ii=1.-ivf_iw-ivf_ir-ivf_if    !S.O. I do not know these are should be calculated or read from input files
    !ivf_ww=1.-ivf_wi-ivf_wr-ivf_wf
    !ivf_rw=1.-ivf_ri-ivf_rf;
    !ivf_fr=ivf_rf;

    IF ((ivf_ii+ivf_iw+ivf_ir+ivf_if > 1.0001) .OR. &
         (ivf_wi+ivf_ww+ivf_wr+ivf_wf > 1.0001) .OR. &
         (ivf_ri+ivf_rw+ivf_rf > 1.0001) .OR. &
         (ivf_fi+ivf_fw+ivf_fr > 1.0001) .OR. &
         (ivf_ii+ivf_iw+ivf_ir+ivf_if < 0.9999) .OR. &
         (ivf_wi+ivf_ww+ivf_wr+ivf_wf < 0.9999) .OR. &
         (ivf_ri+ivf_rw+ivf_rf < 0.9999) .OR. &
         (ivf_fi+ivf_fw+ivf_fr < 0.9999)) THEN
       PRINT*, "At least one internal view factor <> 1. Check ivf in ESTMinput.nml"
    ENDIF

    !=======Initial setting==============================================
    !! Rewritten by HCW 15 Jun 2016 to use existing SUEWS error handling
    !IF(inittemps) THEN
    !   write(*,*) 'inittemps:',inittemps
    !   FileFinalTemp=TRIM(FileOutputPath)//TRIM(FileCodeX)//'_ESTM_finaltemp.txt'
    !   OPEN(99,file=TRIM(FileFinalTemp),status='old',err=316)  ! Program stopped if error opening file
    !   READ(99,*) Twall,Troof,Tground,Tibld                    ! Twall, Troof, Tground & Tibld get new values
    !   CLOSE(99)
    !ENDIF
    !
    !!IF (inittemps) THEN                                                        !!FO!! inittemps=.true. set in nml file
    !!   OPEN(99,file='outputfiles/finaltemp.txt',status='old',iostat=ios)       !!FO!! has to exist
    !!
    !!   IF (ios/=0) CALL error('outputfiles/finaltemp.txt',ios,1)               !!FO!! calls mod_error.f95, writes that the opening failed and stops prg
    !!   IF (ios/=0) THEN
    !!      Twall   = (/273., 285., 291./)
    !!      Troof   = (/273., 285., 291./)
    !!      Tground = (/273., 275., 280., 290./)
    !!      Tibld   = (/293., 293., 293./)
    !!   ELSE
    !!      READ(99,*) Twall,Troof,Tground,Tibld                             !!FO!! if finaltemp.txt exists Twall[3], Troof[3], Tground[4] & Tibld[3] get new values
    !!      CLOSE(99)
    !!   ENDIF
    !!ENDIF

    !where (isnan(Twall))
    !    Twall = 273
    !endwhere
    !where (isnan(Troof))
    !    Troof = 273
    !endwhere
    !where (isnan(Tground))
    !    Tground = 281
    !endwhere
    !where (isnan(Tibld))
    !    Tibld = 293
    !endwhere

    IF(ESTMStart==1) THEN
       DO i=1,4
          Tw_4(:,i) = Twall  !!FO!! Tw_4 holds three differnet temp:s for each wall layer but the same set for all points of the compass
       ENDDO

       !initialize surface temperatures
       T0_ground=Tground(1)
       T0_wall=Twall(1)
       T0_roof=Troof(1)
       T0_ibld=Tibld(1)
       TN_roof=Troof(nroof)
       TN_wall=Twall(nwall)

       !initialize outgoing longwave   !QUESTION: HCW - Are these calculations compatible with those in LUMPS_NARP?
       LUP_ground=SBConst*EM_ground*T0_ground**4
       LUP_WALL=SBConst*EM_WALL*T0_WALL**4
       LUP_ROOF=SBConst*EM_ROOF*T0_ROOF**4

       !  PRINT*,"W,WB= ",W,WB
       !  PRINT*,'SVF_ground ','SVF_WALL ','zvf_WALL ','HW '
       !  PRINT*,SVF_ground,SVF_WALL,zvf_WALL,HW
       !  PRINT*,'RVF_ground ','RVF_WALL ','RVF_ROOF ','RVF_VEG'
       !  PRINT*,RVF_ground,RVF_WALL,RVF_ROOF,RVF_VEG
       !  print*,'Alb_avg (VF)=',alb_avg
       !  print*,'z0m, Zd', z0m, ZD


    ENDIF

    first=.TRUE.

    !======Courant-Friedrichs-Lewy condition=================================
    !This is comment out by S.O. for now
    ! NB: should be recovered for calculation stability, FO and TS, 11 Oct 2017
    !   CFLval = minval(0.5*zibld*zibld*ribld/kibld)   !!FO!! z*z*r/k => unit [s]
    !   if (Tstep>CFLval) then !CFL condition   !!FO!! CFL condition:  Courant�Friedrichs�Lewy condition is a necessary condition for convergence while solving
    !      write(*,*) "IBLD: CFL condition: Tstep=",Tstep,">",CFLval !!FO!! certain partial differential equations numerically by the method of finite differences (like eq 5 in Offerle et al.,2005)
    !      CFLfail=.TRUE.
    !   endif
    !   CFLval = minval(0.5*zroof*zroof*rroof/kroof)
    !   if (Tstep>CFLval) then !CFL condition
    !      write(*,*) "ROOF: CFL condition: Tstep=",Tstep,">",CFLval
    !      CFLfail=.TRUE.
    !   endif
    !   CFLval = minval(0.5*zwall*zwall*rwall/kwall)
    !   if (Tstep>CFLval) then !CFL condition
    !      write(*,*) "WALL: CFL condition: Tstep=",Tstep,">",CFLval
    !      CFLfail=.TRUE.
    !   endif
    !   CFLval = minval(0.5*zground*zground*rground/kground)
    !   if (Tstep>CFLval) then !CFL condition
    !      write(*,*) "ground: CFL condition: Tstep=",Tstep,">",CFLval
    !      CFLfail=.TRUE.
    !   endif
    !   if (CFLfail) then
    !      write(*,*) "Increase dX or decrease maxtimestep. Hit any key to continue"
    !      read (*,*)
    !   endif


    ! Tiaircyc = (1+(LondonQSJune_Barbican.Tair-Tiair)./(5*Tiair)).*(Tiair + 0.4*sin(LondonQSJune_Barbican.HOUR*2*pi/24-10/24*2*pi))    !!FO!! outdoor temp affected



    RETURN

    !     315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,NotUsedI)
    ! 316 CALL errorHint(11,TRIM(fileFinalTemp),notUsed,notUsed,NotUsedI)

  END SUBROUTINE ESTM_translate

  !===============================================================================
  SUBROUTINE ESTM(&
       Gridiv,&!input
       nsh,tstep,&
       avkdn, avu1, temp_c, zenith_deg, avrh, press_hpa, ldown,&
       bldgh,Ts5mindata_ir,&
       Tair24HR,&!inout
       dataOutLineESTM,QS)!output
    ! NB: HCW Questions:
    !                - should TFloor be set in namelist instead of hard-coded here?
    !                - zref used for radiation calculation and fair is set to 2*BldgH here. For compatibility with the rest of the
    !                  SUEWS model, should this be the (wind speed) measurement height z specified in RunControl.nml?
    !                - In SUEWS_translate, fwall=AreaWall/SurfaceArea. Is this correct?
    !                - If froof=1 (i.e. whole grid is building), is HW=0 correct?
    !                - Then is an IF(Fground ==0) option needed?
    !                - alb_wall=0.23 and em_wall=0.9 are set in LUMPS_module_constants. Shouldn't these be provided as input?
    !                - Do the LUP calculations here need to be compatible with those in LUMPS_NARP?
    !                - File opening rewritten using existing error handling in SUEWS - can delete mod_error module from SUEWS_ESTM_functions
    !                - In SUEWS_ESTM_v2016, the first row is set to -999. This may be acceptable at the
    !                  start of the run but should be handled properly between blocks of met data? - need to check what's actually happening here.
    !                - Many duplicate functions in SUEWS_ESTM_functions need changing to the existing SUEWS versions.
    !                - Are the following correctly initialised? T0_ibld, T0_ground, T0_wall, T0_roof, TN_wall, TN_roof, Tground, Twall, Troof, Tibld
    !                - What are Nalb, sumalb, Nemis and Sumemis for? Are they used correctly?
    !

    !SUEWS_ESTM_v2016
    ! revision:
    ! HCW 14 Jun 2016
    ! HCW 27 Jun 2016 Corrected iESTMcount bug - now increases for all grids together
    ! HCW 30 Jun 2016 Major changes to handle grids and met blocks
    ! TS  09 Oct 2017 Added explicit interface

    !Contains calculation for each time step
    !Calculate local scale heat storage from single building and surroundings observations
    !OFferle, May 2003
    !
    !MODIFICATION HISTORY
    !             15 DECEMBER 2003
    !             (1) CHANGED AIR EXCHANGE RATE TO ALSO BE DEPENDENT ON OUTSIDE AIR TEMPERATURE
    !             (2) ADDED SPECIFIC HEAT CALCULATION FOR AIR
    !             (3) RH ADDED AS VAR(14) IN INPUT FILE
    !  12 JANUARY 2004
    !             (1) ADDED ESTIMATED AVG NADIR LOOKING EXTERNAL SURFACE TEMPERATURE TO OUTPUT
    !                 WEIGHTED BY SURFACE FRACTION AND SVF. SVF_ROOF=1, SVF_WALLS = 1-SVF_ground
    !             (2) ADDED NET RADIATION (RN) CALCULATIONS FOR ALL SURFACES AND AVG RN TO OUTPUT BASED ON
    !                 RADIOMETER VIEW FROM ZREF
    !       14 JANUARY 2004
    !             (1) MOVED GRID SPECIFIC PARAMETERS OUT OF NAMELIST, INTO PARAMETER FILE E.G. ALB,EM,F
    !                 T0 MAKE CHANGES IN LOCATIONS EASIER.
    !
    !       23 JANUARY 2004 - Version 2
    !             (1) PUT ALL TEMPERATURES INTO K
    !             (2) added interpolation routine to run on shorter timesteps.
    !                 tested so that it doesn't change solution for forced temperatures.
    !                 however in version 2 the energy balance at the surface isn't correctly solved
    !
    !       25 JANUARY 2004 - Version 3
    !             (1) added solution to energy balance at surface
    !             (2) removed some extraneous code
    !             (3) added vegetation fractions for future development
    !             (4) some changes to input namelist.
    !             (5) need to add wind speed dependence for exchange coefficients
    !             (6) changed the way Rn_net is calculated. also radiometer view factor relationships
    !             (7) added a calculation for heat loss/gain to outside air (that going into building mass is storage)
    !                 this is labelled QFBLD which it is in a sense.
    !       6 FEBRUARY 2004
    !             (1) ADDED INTERNAL VIEW FACTOR FILE FOR INTERNAL GEOMETRY, INCLUDING FIXED FLOOR TEMPERATURE
    !             (2) added MeanU, MinWS to config, and U to ILOC.
    !
    !      11 FEBRUARY 2004
    !             (1) added site lat, long, elevation to inputs
    !             (2) need zenith angle for wall direct radiation interception
    !
    !      15 JUNE 2004
    !             (1) CORRECTED OUTPUT OF T0 FOR FORCED SURFACE TEMPERATURES
    !             (2) MADE SOME CHANGES TO RADIATION, COMPUTATION OF AVERAGE ALBEDO, ADDED AVERAGE EMISSIVITY
    !             (3) ADDED OUTPUT FILE FOR RADIATION COMPONENTS
    !             (4) CHANGED INPUT CONFIG SO CONVECTIVE EXCHANGE COEFFS ARE NOT USER SELECTABLE
    !             (5) MAY STILL BE PROBLEMS WITH WALL RADIATIVE EXCHANGE AND AVERAGE ALBEDO
    !             (6) CHANGED THE WAY HEAT STORAGE IS COMPUTED IN HEATCONDUCTION MODULE BUT THIS SHOULD NOT CHANGE RESULTS
    !
    !      16 NOVEMBER 2004
    !             (1) CHANGED AIR EXCHANGE RATE TO BE BASES ON DAILY TEMPERATURE CHANGES.
    !             (2) WRITES OUT FINAL LAYER TEMPERATURES AND INCLUDES OPTION TO READ AT BEGINNING.
    !
    !DESCRIPTION: uses explicit time differencing to compute heat conduction through roof, walls,
    !             internal mass, and grounds (elements). Air heat storage is computed from average air temperature
    !             Boundary conditions are determined by measured surface temperature(s) or computed
    !             from energy balance at the surface.
    !             Internal air temperature can either be fixed or allowed to evolve in response
    !             to air mass exchanges and convective heating from internal surfaces.
    !
    !INPUT:
    !FORCING DATA: DTIME,KDOWN,LDOWN,TSURF,TAIR_OUT,TAIR_IN,TROOF,TWALL_AVG,TWALL_N,TWALL_E,TWALL_S,TWALL_W,Tground,RH,U
    !NAMELIST    : HEATSTORAGE_vFO.NML
    !    &config
    !    ifile=FORCING DATA
    !    ofile=OUTPUT FILE
    !    pfile=HEAT STORAGE PARAMETER FILE
    !    Nibld = INTERNAL MASS LAYERS
    !    Nwall = EXTERNAL WALL LAYERS
    !    Nroof = ROOF LAYERS
    !    Nground = ground/SOIL LAYERS
    !    LBC_soil = LOWER BOUNDARY CONDITION FOR ground/SOIL
    !    iloc= INPUT COLUMNS IN DATA FILE
    !    evolveTibld= USE DIAGNOSTIC VALUE FOR INTERNAL BUILDING TEMPERATURE
    !                        0: don't use, use measured
    !                        1: TURN ON USE when temp goess ABOVE TINT_ON, off when temp is below TINT_OFF
    !                        2: always use diagnostic
    !    THEAT_ON= TEMPERATURE AT WHICH HEAT CONTROL IS TURNED ON
    !    THEAT_OFF= TEMPERATURE AT WHICH HEAT CONTROL IS TURNED OFF
    !       THEAT_FIX = Fixed internal temperature for climate control
    !    oneTsurf= USE SINGLE SURFACE TEMPERATURE TO DRIVE ALL LAYERS
    !    radforce= USE RADIATIVE ENERGY BALANCE TO DRIVE EXTERNAL TEMPERATURES
    !    maxtimestep=302, maximum time step in s for filling data gaps.
    !       Alt = STATION HEIGHT (m) FOR PRESSURE CALCULATION
    !       SPINUP = NUMBER OF LINES TO USE FOR SPINUP (REPEATS THESE LINES BUT ONLY OUTPUTS THE 2ND TIME)
    !    INITTEMP = if TRUE INITIALIZES TEMPERATURES TO THOSE IN FINALTEMP.TXT FILE
    !    CH_ibld = INTERNAL BUILDING CONVECTIVE EXCHANGE COEFFICIENT
    !       **** THESE SHOULD DEPEND ON WIND SPEED BUT CURRENTLY DO NOT ****
    !       CHAIR = CONVECTIVE EXCHANGE COEFFICIENT FOR ROOF
    !       chair_ground = ... FOR ground
    !       chair_wall = ... FOR WALL
    !    /
    ! ***************** PARAMETER FILE VARIABLES
    !               fveg = FRACTION OF ground SURFACE COVERED WITH VEG
    !               zveg = VEGETATION HEIGHT
    !               alb_veg = VEGETATION ALBEDO
    !               em_veg = VEGETATION EMISSIVITY
    !               ZREF = REFERENCE HEIGHT FOR FLUX CALCULATION
    !               BldgH    = mean building height
    !               HW    = CANYON ASPECT RATION
    !               f_X   = FRACTION OF X WHERE X IS INTERNAL, WALL, ROOF, ground
    !               Alb_x = ALBEDO OF X
    !               em_ibld = EMISSIVITY OF X
    !               TX    = INITIAL LAYER TEMPERATURES
    !               zX    = LAYER THICKNESS
    !               kX    = LAYER THERMAL CONDUCTIVITY
    !               ribld = LAYER VOLUMETRIC HEAT CAPACITY
    !
    !****************** INTERNAL VIEW FACTOR FILE
    !OUTPUT:      fixed format text with single header, heatstorage for all elements, and temperatures
    !             for each element-layer.
    !===============================================================================

    USE meteo, ONLY: pi, heatcapacity_air
    USE mod_solver
    USE modSolarCalc                                                        !!FO!! :modsolarcalc.f95
    USE MathConstants                                                       !!FO!! :MathConstants_module.f95
    USE PhysConstants
    USE heatflux
    USE ESTM_data

    IMPLICIT NONE
    INTEGER, PARAMETER:: ncolsESTMdata=13
    ! INTEGER, PARAMETER:: ncolumnsDataOutESTM=32
    INTEGER, PARAMETER:: cTs_Tiair = 5
    INTEGER, PARAMETER:: cTs_Tsurf = 6
    INTEGER, PARAMETER:: cTs_Troof = 7
    INTEGER, PARAMETER:: cTs_Troad = 8
    INTEGER, PARAMETER:: cTs_Twall = 9
    INTEGER, PARAMETER:: cTs_Twall_n = 10
    INTEGER, PARAMETER:: cTs_Twall_e = 11
    INTEGER, PARAMETER:: cTs_Twall_s = 12
    INTEGER, PARAMETER:: cTs_Twall_w = 13
    REAL(KIND(1d0)),PARAMETER::NAN=-999

    INTEGER,INTENT(in)::Gridiv
    INTEGER,INTENT(in)::nsh,tstep
    ! INTEGER,INTENT(in)::iy !Year
    ! INTEGER,INTENT(in)::id !Day of year
    ! INTEGER,INTENT(in)::it !Hour
    ! INTEGER,INTENT(in)::imin          !Minutes

    REAL(KIND(1d0)),INTENT(in)::avkdn
    REAL(KIND(1d0)),INTENT(in)::avu1
    REAL(KIND(1d0)),INTENT(in)::temp_c
    REAL(KIND(1d0)),INTENT(in)::zenith_deg
    REAL(KIND(1d0)),INTENT(in)::avrh
    REAL(KIND(1d0)),INTENT(in)::press_hpa
    REAL(KIND(1d0)),INTENT(in)::ldown
    REAL(KIND(1d0)),INTENT(in)::bldgh
    ! REAL(KIND(1d0)),INTENT(in):: dectime        !Decimal time
    REAL(KIND(1d0)),DIMENSION(ncolsESTMdata),INTENT(in)::  Ts5mindata_ir     !surface temperature input data
    REAL(KIND(1d0)),DIMENSION(24*nsh),INTENT(inout) ::   Tair24HR ! may be replaced with MetForcingData by extracting the Tiar part

    REAL(KIND(1d0)),DIMENSION(27),INTENT(out):: dataOutLineESTM
    !Output to SUEWS
    REAL(KIND(1d0)),INTENT(out)::QS
    !Input from SUEWS, corrected as Gridiv by TS 09 Jun 2016


    !Use only in this subroutine
    INTEGER::i, ii
    INTEGER:: Tair2Set=0
    REAL(KIND(1d0))::AIREXHR, AIREXDT
    REAL(KIND(1d0)),DIMENSION(2)::bc
    REAL(KIND(1d0))::chair_ground,chair_wall
    REAL(KIND(1d0))::EM_EQUIV
    REAL(KIND(1d0))::kdz
    REAL(KIND(1d0))::kup_estm,LUP_net,kdn_estm
    REAL(KIND(1d0))::QHestm
    REAL(KIND(1d0))::QFBld !Anthropogenic heat from HVAC
    REAL(KIND(1d0))::shc_airbld
    REAL(KIND(1d0))::sw_hor,sw_vert
    REAL(KIND(1d0))::T0
    REAL(KIND(1d0))::Tinternal,Tsurf_all,Troof_in,Troad,Twall_all,Tw_n,Tw_e,Tw_s,Tw_w
    REAL(KIND(1d0))::Twallout(5),Troofout(5),Tibldout(5),Tgroundout(5)
    REAL(KIND(1d0))::Tadd,Tveg
    REAL(KIND(1d0))::Tairmix
    REAL(KIND(1d0))::RN
    REAL(KIND(1d0))::Rs_roof,Rl_roof,RN_ROOF
    REAL(KIND(1d0))::Rs_wall,Rl_wall,RN_WALL
    REAL(KIND(1d0))::Rs_ground,Rl_ground,RN_ground
    REAL(KIND(1d0))::Rs_ibld,Rl_ibld
    REAL(KIND(1d0))::Rs_iroof,Rl_iroof
    REAL(KIND(1d0))::Rs_iwall,Rl_iwall
    REAL(KIND(1d0))::zenith_rad
    REAL(KIND(1d0))::dum(50)
    REAL(KIND(1d0))::bldgHX ! local bldgHX=max(bldgH,0.001)
    REAL(KIND(1d0)),PARAMETER::WSmin=0.1  ! Check why there is this condition. S.O.
    LOGICAL::radforce, groundradforce

    radforce       = .FALSE.
    groundradforce = .FALSE. !Close the radiation scheme in original ESTM S.O.O.

    bldgHX=MAX(bldgH,0.001) ! this is to prevent zero building height

    ! Set -999s for first row
    ! dum=(/-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
    !      -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
    !      -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
    !      -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,&
    !      -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999./)
    dum=[(-999,i=1,50)]

    !External bulk exchange coefficients - set these somewhere more sensible***
    CHR=0.005
    CHAIR=CHR
    CHAIR_ground=CHAIR
    CHAIR_WALL=CHAIR

    !Get met data for use in ESTM subroutine
    kdn_estm=avkdn
    WS=avu1
    IF (WS<WSMin) WS=WSmin
    Tair1=Temp_C+C2K
    ! Set initial value of Tair2 to air temp
    IF(Gridiv == 1) Tair2Set = Tair2Set+1
    IF(Tair2Set==1) THEN
       Tair2=Temp_C+C2K
    ELSE
       Tair2 = Tair2_grids(Gridiv)
       ! Also get other variables for this grid
       Tievolve = Tievolve_grids(Gridiv)
       lup_ground = lup_ground_grids(Gridiv)
       lup_wall = lup_wall_grids(Gridiv)
       lup_roof = lup_roof_grids(Gridiv)
       T0_ibld = T0_ibld_grids(Gridiv)
       T0_ground = T0_ground_grids(Gridiv)
       T0_wall = T0_wall_grids(Gridiv)
       T0_roof = T0_roof_grids(Gridiv)
       TN_wall = TN_wall_grids(Gridiv)
       TN_roof = TN_roof_grids(Gridiv)
       Tground(:) = Tground_grids(:,Gridiv)
       Twall(:) = Twall_grids(:,Gridiv)
       Troof(:) = Troof_grids(:,Gridiv)
       Tibld(:) = Tibld_grids(:,Gridiv)
       Tw_4 = Tw_4_grids(:,:,Gridiv)

    ENDIF

    ! Get Ts from Ts5min data array
    ! Tinternal  = Ts5mindata(ir,cTs_Tiair)
    ! Tsurf_all  = Ts5mindata(ir,cTs_Tsurf)
    ! Troof_in   = Ts5mindata(ir,cTs_Troof)
    ! Troad      = Ts5mindata(ir,cTs_Troad)
    ! Twall_all  = Ts5mindata(ir,cTs_Twall)
    !
    ! Tw_n       = Ts5mindata(ir,cTs_Twall_n)
    ! Tw_e       = Ts5mindata(ir,cTs_Twall_e)
    ! Tw_s       = Ts5mindata(ir,cTs_Twall_s)
    ! Tw_w       = Ts5mindata(ir,cTs_Twall_w)
    Tinternal  = Ts5mindata_ir(cTs_Tiair)
    Tsurf_all  = Ts5mindata_ir(cTs_Tsurf)
    Troof_in   = Ts5mindata_ir(cTs_Troof)
    Troad      = Ts5mindata_ir(cTs_Troad)
    Twall_all  = Ts5mindata_ir(cTs_Twall)

    Tw_n       = Ts5mindata_ir(cTs_Twall_n)
    Tw_e       = Ts5mindata_ir(cTs_Twall_e)
    Tw_s       = Ts5mindata_ir(cTs_Twall_s)
    Tw_w       = Ts5mindata_ir(cTs_Twall_w)

    !    if (any(isnan(Ts5mindata))) then                                    !!FO!! can't use data when time gap is too big (or neg.) or data is NaN
    !        if (spindone) then                                                  !!FO!! writes a line of NaNs
    !            write(20,'(1F8.4,I6,100f10.1)') dectime,it,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
    !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
    !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,&
    !                -0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.,-0./0.
    !        endif
    !        return ! changed from cycle            !!FO!! returns to beginning of do-loop, i.e. doesn't make any heat calculations
    !    endif
    !!FO!! this loop is used to run through the calculations (SPINUP number of times) to achieve better numerical stability
    !!FO!! when if condition is met the program starts over from the beginning of the input file and the calculations are performed
    !!FO!! ...once again, but this time saved to the output file
    !!FO!! it's because of this if statement arrangement that the output result starts at input time #2 (if not SPINUP=0 as originally in heatstorage_vFO.nml)
    !    IF (NLINESREAD>datalines.AND..NOT.SPINDONE) THEN                        !!FO!! at this stage spindone = .false.
    !        NLINESREAD=0
    !        SPINDONE=.TRUE.
    !PRINT*, "SPUNUP"
    !        return ! changed from cycle
    !    ENDIF

    !! Write first row of each met block as -999
    !IF (first) THEN  !Set to true in ESTM_initials
    !!   !Tair2=Temp_C+C2K !This is now set in SUEWS_translate for ir=0 only
    !!   ! first=.FALSE.
    !   IF(Gridiv == NumberOfGrids) first=.FALSE.  !Set to false only after all grids have run
    !   dataOutESTM(ir,1:32,Gridiv)=(/REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),&
    !        REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)),dectime,(dum(ii),ii=1,27)/)
    !   RETURN
    !ENDIF

    ! What are these constants? - Need defining somewhere
    zenith_rad=zenith_deg/180*PI
    IF (zenith_rad>0.AND.zenith_rad<PI/2.-HW) THEN  !ZENITH MUST BE HIGHER THAN BUILDINGS FOR DIRECT INTERCEPTION
       tanzenith = MIN(TAN(zenith_rad),5.67) !LIMITS TO ANGLES LESS THAN 80 EVEN FOR LOW HW
       tanzenith = tanzenith*kdn_estm/(1370*COS(zenith_rad)) !REDUCTION FACTOR FOR MAXIMUM
    ELSE
       tanzenith = 0.
    ENDIF

    SHC_air=HEATCAPACITY_AIR(Tair1,avrh,Press_hPa)   ! Use SUEWS version
    Tair24HR=EOSHIFT(Tair24HR, 1, Tair1, 1) !!!*** NB: Check this. and is this the tair of past 24 hrs? TS 10 Oct 2017
    Tairday=SUM(Tair24HR)/(24*nsh)


    !Evolution of building temperature from heat added by convection
    SELECT CASE(evolvetibld)   !EvolveTiBld specifies which internal building temperature approach to use
    CASE(0)
       diagnoseTi=.FALSE.
       HVAC=.FALSE. !use data in file                                !!FO!! use measured indoor temperature (Tref in Lodz2002HS.txt)
    CASE(1)                                                                                   !!FO!! use of HVAC to counteract T changes
       diagnoseTi=.TRUE.
       IF (Tievolve>THEAT_OFF) THEN   !THEAT_OFF now converted to Kelvin in ESTM_initials - HCW 15 Jun 2016
          !IF (Tievolve>THEAT_OFF+C2K) THEN
          HVAC=.FALSE.
       ELSEIF (Tievolve<THEAT_ON) THEN   !THEAT_OFF now converted to Kelvin in ESTM_initials - HCW 15 Jun 2016
          !ELSEIF (Tievolve<THEAT_ON+C2K) THEN
          HVAC=.TRUE.
       ENDIF
    CASE(2)
       diagnoseTi=.TRUE.                                                                 !!FO!! convection between ibld and inside of external walls(?)
    END SELECT

    !ASSUME AIR MIXES IN PROPORTION TO # OF EXCHANGES
    IF (Tairday>20.+C2K.AND.Tievolve>25.+C2K.AND.TAIR1<Tievolve.AND..NOT.HVAC) THEN
       AIREXHR = 2.0  !Windows or exterior doors on 3 sides (ASHRAE 1981 22.8)
    ELSEIF (Tairday<17.+C2K.OR.HVAC) THEN
       AIREXHR = 0.5 !No window or exterior doors, storm sash or weathertripped (ASHRAE 1981 22.8)
    ELSE
       AIREXHR = 1.0
    ENDIF

    AIREXDT=AIREXHR*(Tstep/3600.0)
    shc_airbld=MAX(HEATCAPACITY_AIR(TiEVOLVE,avrh,Press_hPa),0.00001) ! to avoid zero-division scenario TS 21 Oct 2017
    IF (shc_airbld<minshc_airbld) minshc_airbld=shc_airbld

    !internal convective exchange coefficients                         !!FO!! ibldCHmod = 0 originally
    !iBldCHmod specifies method for convective exchange coeffs
    IF (ibldCHmod==1) THEN       !ASHRAE 2001
       CH_ibld  = 1.31*(ABS(T0_ibld-Tievolve))**0.25/shc_airbld
       CH_iwall = 1.31*(ABS(TN_wall-Tievolve))**0.25/shc_airbld
       CH_iroof = 1.52*(ABS(TN_roof-Tievolve))**0.25/shc_airbld
       IF (ABS(TN_roof-Tievolve)>0) CH_iroof=CH_iroof*0.39 !effect of convection is weaker downward
    ELSEIF (ibldCHmod==2) THEN   !Awbi, H.B. 1998, Energy and Buildings 28: 219-227
       CH_ibld  = 1.823*(ABS(T0_ibld-Tievolve))**0.293/shc_airbld
       CH_iwall = 1.823*(ABS(TN_wall-Tievolve))**0.293/shc_airbld
       CH_iroof = 2.175*(ABS(TN_roof-Tievolve))**0.308/shc_airbld
       IF (ABS(TN_roof-Tievolve)>0) CH_iroof=0.704*(ABS(TN_roof-Tievolve))**0.133/shc_airbld !effect of convection is weaker downward
    ENDIF

    !Evolving T = (Previous Temp + dT from Sensible heat flux) mixed with outside air
    !ASSUMES THE CH_BLD INCLUDES THE EFFECT OF VENTILATION RATE IN m/s (e.g. if a normal CH is .005 and
    !the value here is .003 the assumed ventilation is 0.6 m/s                                 !!FO!! CH_ibld=0.0015 from heatstorage_Barbican.nml => ventilation=0.3 m/s
    Tairmix =  (Tievolve + TAIR1*AIREXDT)/(1.0+AIREXDT)
    QFBld= froof*(Tievolve-Tairmix)*shc_airbld*bldgHX/Tstep !heat added or lost, requires cooling or heating if HVAC on

    !!FO!! CH_xxxx has unit [m/s]  !!**HCW what is going on with tstep here??
    Tievolve = Tairmix+Tstep/bldgHX/finternal* &                                                                         !!FO!! finternal(=froof+fibld+fwall) => normalisation of fractions
         (CH_ibld*fibld*(T0_ibld-Tievolve)+CH_iroof*froof*(TN_roof-Tievolve)+CH_iwall*fwall*(TN_wall-Tievolve))      !!FO!! [K] = [K] + [s/m]*([m/s]*([K]))

    IF (.NOT.diagnoseTi) Tievolve=Tinternal+C2K
    IF (HVAC) THEN !Run up/down to set point +/- 1 degree with adjustment of 90% per hour
       Tadd=(SIGN(-1.0d0,THEAT_fix-Tievolve)+THEAT_fix-Tievolve)*MIN(4.*Tstep/3600.0,0.9) !!**HCW check??
       Tievolve=Tievolve+Tadd
    ENDIF


    !========>RADIATION<================
    IF (kdn_estm<0) kdn_estm=0. !set non-zero shortwave to zero  !Should this be moved up to line 183/4?

    !external components, no diffuse
    !for reflections complete absorption is assumed
    !for shortwave these are net values
    !for longwave these are incoming only
    !MUST DIVIDE SHORTWAVE INTO DIRECT AND DIFFUSE
    sw_hor =kdn_estm           !incoming solar on horizontal surface
    sw_vert=kdn_estm*tanzenith !incoming solar on vertical surface = kdown(obs)*sin(zenith)/cos(zenith)

    Rs_roof=svf_roof*(1.0-alb_roof)*sw_hor
    Rl_roof=svf_roof*em_roof*ldown

    Rs_ground=svf_ground*(1.-alb_ground)*sw_hor+&
         zvf_ground*svf_wall*alb_wall*sw_vert*(1-alb_ground)+&
         zvf_ground*svf_ground*alb_ground*sw_hor*xvf_wall*alb_wall

    Rl_ground=svf_ground*ldown*em_ground+zvf_ground*(lup_wall+svf_wall*ldown*(1-em_wall))*em_ground

    Rs_wall=svf_wall*(1.-alb_wall)*sw_vert+&
         zvf_wall*svf_wall*alb_wall*sw_vert*(1.+zvf_wall*alb_wall)+&
         xvf_wall*svf_ground*alb_ground*sw_hor*(1-alb_wall)+&
         zvf_ground*xvf_wall*svf_ground*alb_ground*sw_hor*alb_wall

    !wall to wall exchange handled simultaneously with seb calc
    Rl_wall=svf_wall*ldown*em_wall+zvf_wall*svf_wall*ldown*(1-em_wall)*em_wall+&
         xvf_wall*(lup_ground+svf_ground*ldown*(1-em_ground))*em_wall

    !DIFFICULT TO DETERMINE WHAT THIS IS EXACTLY, DONT INCLUDE WALLS
    kup_estm=kdn_estm-RVF_ROOF*Rs_roof-(RVF_ground+RVF_WALL)*Rs_ground/svf_ground-RVF_VEG*ALB_VEG*kdn_estm
    IF (kdn_estm > 10 .AND. kup_estm > 0) THEN
       alb_avg = kup_estm/kdn_estm
       sumalb  = sumalb+alb_avg
       Nalb    = Nalb+1
    ENDIF


    !internal components
    Rs_ibld=0 ! This could change if there are windows (need solar angles or wall svf * fraction glazing * transmissivity)
    !internal incoming longwave terms do not include the view factors for its own surface e.g. for ibld and walls
    !added floor view factors
    Rl_ibld=SBConst*(ivf_iw*em_w*TN_wall**4 +&
         ivf_ir*em_r*TN_roof**4 +&
         ivf_if*em_f*Tfloor**4)
    Rs_iwall=0
    Rl_iwall=SBConst*(ivf_wi*em_i*T0_ibld**4 +&
         ivf_wr*em_r*TN_roof**4 +&
         ivf_wf*em_f*Tfloor**4)
    Rs_iroof=0
    Rl_iroof=SBConst*(ivf_ri*em_i*T0_ibld**4 +&
         ivf_rw*em_w*TN_wall**4 +&
         ivf_rf*em_f*Tfloor**4)

    !========>INTERNAL<================
    bctype=.FALSE.
    kdz=2*kibld(1)/zibld(1)
    Pcoeff=(/em_ibld*SBConst*(1-ivf_ii*em_ibld),0.0d0,0.0d0,kdz+shc_airbld*CH_ibld,&
         -kdz*Tibld(1)-shc_airbld*CH_ibld*Tievolve-Rs_ibld-Rl_ibld/)
    T0_ibld=NewtonPolynomial(T0_ibld,Pcoeff,conv,maxiter)
    bc(1)=T0_ibld                                                       !!FO!! this leads to Tibld(1) = Tibld(3) , i.e. ...
    bc(2)=bc(1)                                                         !!FO!! temperature equal on both sides of inside wall
    CALL heatcond1d(Tibld,Qsibld,zibld(1:Nibld),REAL(Tstep,KIND(1d0)),kibld(1:Nibld),ribld(1:Nibld),bc,bctype)

    !========>WALLS<================
    bctype=.FALSE.
    kdz=2*kwall(nwall)/zwall(nwall)
    Pcoeff=(/em_ibld*SBConst*(1-ivf_ww*em_ibld),0.0d0,0.0d0,kdz+shc_airbld*CH_iwall,&
         -kdz*Twall(nwall)-shc_airbld*CH_iwall*Tievolve-Rs_iwall-Rl_iwall/)
    TN_wall=NewtonPolynomial(TN_wall,Pcoeff,conv,maxiter)
    bc(2)=TN_wall                                                       !!FO!! boundary condition #2 = inner surface Twall, originally from lodz_parms_ltm.txt or finaltemp.txt

    IF (TsurfChoice<2 .OR.radforce) THEN
       IF (radforce) THEN                                              !!FO!! 1st prio: radforce
          kdz=2*kwall(1)/zwall(1)
          Pcoeff=(/em_wall*SBConst*(1-zvf_wall*em_wall),0.0d0,0.0d0,kdz+shc_air*chair_wall*WS,&
               -kdz*Twall(1)-shc_air*chair_wall*WS*Tair1-Rs_wall-Rl_wall/)
          T0_wall=NewtonPolynomial(T0_wall,Pcoeff,conv,maxiter)
          bc(1)=T0_wall                                               !!FO!! boundary condition #1 = outer surface Twall, originally from lodz_parms_ltm.txt or finaltemp.txt
       ELSEIF (TsurfChoice==0) THEN
          bc(1)=Tsurf_all+C2K; T0_wall=bc(1)
       ELSEIF (TsurfChoice==1) THEN
          bc(1)=Twall_all+C2K; T0_wall=bc(1)
       ENDIF                                                           !!FO!! Tsoil in Lodz2002HS.txt NB => Lodz2002HS.txt doesn't work with onewall = TRUE

       CALL heatcond1d(Twall,Qswall,zwall(1:nwall),REAL(Tstep,KIND(1d0)),kwall(1:nwall),rwall(1:nwall),bc,bctype)     !!FO!! new set of Twalls are calculated from heat conduction through wall

    ELSEIF(TsurfChoice==2) THEN!SPECIAL FOR 4 WALLS
       T0_wall=0.
       DO i=1,4 !do 4 walls
          bc(1)=Tw_n+Tw_e+Tw_s+Tw_w+C2K; T0_wall=T0_wall+bc(1)
          CALL heatcond1d(Tw_4(:,i),Qs_4(i),zwall(1:nwall),REAL(Tstep,KIND(1d0)),kwall(1:nwall),rwall(1:nwall),bc,bctype)
       ENDDO
       !Take average of 4 wall values
       T0_wall=T0_wall/4.
       Qswall = SUM(Qs_4)/4.
       Twall = SUM(Tw_4,2)/4.
    ENDIF

    !========>ROOF<================
    bctype=.FALSE.
    kdz=2*kroof(nroof)/zroof(nroof)
    Pcoeff=(/em_ibld*SBConst,0.0d0,0.0d0,kdz+shc_airbld*CH_iroof,&
         -kdz*Troof(nroof)-shc_airbld*CH_iroof*Tievolve-Rs_iroof-Rl_iroof/)
    TN_roof=NewtonPolynomial(TN_roof,Pcoeff,conv,maxiter)
    bc(2)=TN_roof

    IF (radforce) THEN
       kdz=2*kroof(1)/zroof(1)
       Pcoeff=(/em_roof*SBConst,0.0d0,0.0d0,kdz+shc_air*chair*WS,&
            -kdz*Troof(1)-shc_air*chair*WS*Tair1-Rs_roof-Rl_roof/)
       T0_roof=NewtonPolynomial(T0_roof,Pcoeff,conv,maxiter)
       bc(1)=T0_roof
    ELSEIF (TsurfChoice==0) THEN
       bc(1)=Tsurf_all+C2K; T0_roof=bc(1)
    ELSE
       bc(1)=Troof_in+C2K; T0_roof=bc(1)
    ENDIF

    CALL heatcond1d(Troof,Qsroof,zroof(1:nroof),REAL(Tstep,KIND(1d0)),kroof(1:nroof),rroof(1:nroof),bc,bctype)


    !========>ground<================
    bctype=.FALSE.
    kdz=2*kground(1)/zground(1)

    IF (radforce.OR.groundradforce) THEN
       Pcoeff=(/em_ground*SBConst,0.0d0,0.0d0,kdz+shc_air*chair_ground*WS,&
            -kdz*Tground(1)-shc_air*chair_ground*WS*Tair1-Rs_ground-Rl_ground/)
       T0_ground=NewtonPolynomial(T0_ground,Pcoeff,conv,maxiter)
       bc(1)=T0_ground
    ELSEIF (TsurfChoice==0) THEN
       bc(1)=Tsurf_all+C2K; T0_ground=bc(1)
    ELSE
       bc(1)=Troad+C2K; T0_ground=bc(1)
    ENDIF

    bc(2)=LBC_soil+C2K
    !     bc(2)=0.; bctype(2)=.t.

    IF ( fground/=0. )   THEN   ! check fground==0 scenario to avoid division-by-zero error, TS 21 Jul 2016
       CALL heatcond1d(Tground,Qsground,zground(1:Nground),REAL(Tstep,KIND(1d0)),kground(1:Nground),rground(1:Nground),bc,bctype)
    ELSE
       Qsground=NAN
    END IF

    Qsair = fair*SHC_air*(Tair1-Tair2)/Tstep
    Qsibld = Qsibld*fibld
    Qswall = Qswall*fwall
    Qsroof = Qsroof*froof
    Qsground = Qsground*fground
    QS = Qsibld + Qswall + Qsroof + Qsground                              !!FO!! QSair not included; called QS in output file (column #10)


    !write(*,*) Qsair, QSibld, Qswall, Qsroof, Qsground, QS

    !========>Radiation<================
    !note that the LUP for individual components does not include reflected
    LUP_ground = SBConst*EM_ground*T0_ground**4
    LUP_WALL   = SBConst*EM_WALL*T0_WALL**4
    LUP_ROOF   = SBConst*EM_ROOF*T0_ROOF**4
    TVEG       = TAIR1
    LUP_VEG    = SBConst*EM_VEG*TVEG**4
    T0         = RVF_ground*T0_ground+RVF_WALL*T0_WALL+RVF_ROOF*T0_ROOF+RVF_VEG*TVEG
    LUP_net    = RVF_ground*LUP_ground+RVF_WALL*LUP_WALL+RVF_ROOF*LUP_ROOF+RVF_VEG*LUP_VEG
    EM_EQUIV   = LUP_net/(SBConst*T0**4) !!FO!! apparent emissivity of the atmosphere [cloudless sky: >� Ldown from gases in the lowest 100 m] calculated from surface at T0
    RN_ground  = rs_ground+rl_ground-lup_ground
    RN_ROOF    = rs_roof+rl_roof-lup_roof
    RN_WALL    = rs_wall+rl_wall-lup_wall*(1-zvf_wall*em_wall)
    RN         = kdn_estm-kup_estm+ldown*EM_EQUIV-lup_net !!FO!! average net radiation (at z > zref ????) = shortwave down - shortwave up + [longwave down * apparent emissivity] - longwave up
    QHestm     = (T0-Tair1)*CHair*SHC_air*WS
    sumemis    = sumemis+EM_EQUIV
    nemis      = nemis+1

    ! IF (SPINDONE) THEN                                                      !!FO!! only the last set of values in the time interpolation loop is written to file

    IF (Nwall<5)THEN
       Twallout  =(/Twall,(dum(ii),  ii=1,(5-Nwall))/)
    ELSE
       Twallout=Twall
    ENDIF

    IF (Nroof<5) THEN
       Troofout  =(/Troof,(dum(ii),  ii=1,(5-Nroof))/);
    ELSE
       Troofout=Troof
    ENDIF

    IF (Nground<5)THEN
       Tgroundout=(/Tground,(dum(ii),ii=1,(5-Nground))/)
    ELSE
       Tgroundout=Tground
    ENDIF

    IF (Nibld<5)THEN
       Tibldout  =(/Tibld,(dum(ii),  ii=1,(5-Nibld))/)
    ELSE
       Tibldout=Tibld
    ENDIF

    ! dataOutESTM(ir,1:ncolumnsDataOutESTM,Gridiv)=[&
    !      REAL(iy,KIND(1D0)),REAL(id,KIND(1D0)),REAL(it,KIND(1D0)),REAL(imin,KIND(1D0)), dectime,&!5
    !      QS,Qsair,Qswall,Qsroof,Qsground,Qsibld,&!11
    !      Twallout,Troofout,Tgroundout,Tibldout,Tievolve]!32 !NB. These all have 5 elements except Tievolve (1).
    dataOutLineESTM=[&
         QS,Qsair,Qswall,Qsroof,Qsground,Qsibld,&!6
         Twallout,Troofout,Tgroundout,Tibldout,Tievolve]!27 !NB. These all have 5 elements except Tievolve (1).
    ! set invalid values to nan
    dataOutLineESTM=set_nan(dataOutLineESTM)
    ! call r8vec_print(ncolumnsDataOutESTM-5,dataOutESTM(ir,6:ncolumnsDataOutESTM,Gridiv),'dataOutESTM')


    Tair2=Tair1

    ! Save variables for this grid
    Tair2_grids(Gridiv)=Tair1
    lup_ground_grids(Gridiv) = lup_ground
    lup_wall_grids(Gridiv) = lup_wall
    lup_roof_grids(Gridiv) = lup_roof
    Tievolve_grids(Gridiv) = Tievolve
    T0_ibld_grids(Gridiv) = T0_ibld
    T0_ground_grids(Gridiv) = T0_ground
    T0_wall_grids(Gridiv) = T0_wall
    T0_roof_grids(Gridiv) = T0_roof
    TN_wall_grids(Gridiv) = TN_wall
    TN_roof_grids(Gridiv) = TN_roof
    Tground_grids(:,Gridiv) = Tground(:)
    Twall_grids(:,Gridiv) = Twall(:)
    Troof_grids(:,Gridiv) = Troof(:)
    Tibld_grids(:,Gridiv) = Tibld(:)
    Tw_4_grids(:,:,Gridiv) = Tw_4(:,:)

  END SUBROUTINE ESTM

  !===============set variable of invalid value to NAN====================================
  ELEMENTAL FUNCTION set_nan(x) RESULT(xx)
    IMPLICIT NONE
    REAL(KIND(1d0)),PARAMETER::pNAN=9999
    REAL(KIND(1d0)),PARAMETER::NAN=-999
    REAL(KIND(1d0)),INTENT(in)::x
    REAL(KIND(1d0))::xx

    IF(ABS(x)>pNAN) THEN
       xx=NAN
    ELSE
       xx=x
    ENDIF

  END FUNCTION set_nan
  !========================================================================

END MODULE ESTM_module
