!========================================================================================
!> AnOHM: Analytical Objective Hysteresis Model
!> @author
!> Ting Sun, ting.sun@reading.ac.uk
!> @brief
!> calculate heat storage.
!> @ref model details refer to
!> https://doi.org/10.5194/gmd-2016-300
! history:
! 20160301: initial version
! 20170109: updated dqndt calculation in accordance with SUEWS_OHM.f95 (HCW)
! 20170810: revamped structure
! 20170825: improved Bowen calculation
!========================================================================================
MODULE AnOHM_module

  IMPLICIT NONE
CONTAINS

  !========================================================================================
  !> High level wrapper for AnOHM calculation
  !! @brief
  !! calculate heat storage based within AnOHM framework.
  !! @returns
  !! -# grid ensemble heat storage:
  !! QS = a1*(Q*)+a2*(dQ*/dt)+a3
  !! -# grid ensemble OHM coefficients: a1, a2 and a3
  SUBROUTINE AnOHM(&
       qn1,qn1_store,qn1_av_store,&
       MetForcingData_grid,moist_surf,&
       alb, emis, cp, kk, ch,&! input
       sfr,nsurf,nsh,AnthropHeatMethod,id,Gridiv,&
       a1,a2,a3,qs)! output

    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:,:)::MetForcingData_grid !< met forcing array of grid

    REAL(KIND(1d0)),INTENT(in):: qn1               !< net all-wave radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: sfr(nsurf)        !< surface fraction (0-1) [-]
    REAL(KIND(1d0)),INTENT(in):: moist_surf(nsurf) !< non-dimensional surface wetness status (0-1) [-]

    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::alb  !< albedo [-]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::emis !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::cp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::kk   !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::ch   !< bulk transfer coef [J m-3 K-1]

    INTEGER,INTENT(in):: id                !< day of year [-]
    INTEGER,INTENT(in):: Gridiv            !< grid id [-]
    INTEGER,INTENT(in):: AnthropHeatMethod !< AnthropHeat option [-]
    INTEGER,INTENT(in):: nsurf             !< number of surfaces [-]
    INTEGER,INTENT(in):: nsh               !< number of timesteps in one hour [-]

    REAL(KIND(1d0)),INTENT(inout)::qn1_store(nsh) !< stored qn1 [W m-2]
    REAL(KIND(1d0)),INTENT(inout)::qn1_av_store(2*nsh+1) !< average net radiation over previous hour [W m-2]

    REAL(KIND(1d0)),INTENT(out):: a1 !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out):: a2 !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out):: a3 !< AnOHM coefficients of grid [W m-2]
    REAL(KIND(1d0)),INTENT(out):: qs !< storage heat flux [W m-2]

    INTEGER :: is,xid !< @var qn1 net all-wave radiation
    REAL(KIND(1d0)),PARAMETER::NotUsed=-55.5!< @var qn1 net all-wave radiation
    INTEGER,PARAMETER::notUsedI=-55!< @var qn1 net all-wave radiation
    LOGICAL :: idQ ! whether id contains enough data

    REAL(KIND(1d0))                  :: dqndt       !< rate of change of net radiation [W m-2 h-1] at t-2
    ! REAL(KIND(1d0))                  :: surfrac     !< surface fraction accounting for SnowFrac if appropriate
    REAL(KIND(1d0)),DIMENSION(nsurf) :: xa1,xa2,xa3 !< temporary AnOHM coefs.
    ! REAL(KIND(1d0))                  :: qn1_av      ! average net radiation over previous hour [W m-2]
    ! REAL(KIND(1d0))                  :: nsh_nna     ! number of timesteps per hour with non -999 values (used for spinup)

    ! initialize the coefficients
    xa1 = 0.1
    xa2 = 0.2
    xa3 = 10

    ! to test if the current met block contains enough data for AnOHM
    ! TODO: more robust selection should be implemented
    ! daylight hours >= 6
    idQ=COUNT(MetForcingData_grid(:,2)==id .AND. & ! day of year
         MetForcingData_grid(:,4)==0 .AND. & ! minutes
         MetForcingData_grid(:,15)>0) & ! Sd
         .GE. 6

    ! PRINT*, idQ
    IF ( idQ ) THEN
       ! given enough data, calculate coefficients of day `id`
       xid=id
    ELSE
       ! otherwise calculate coefficients of yesterday: `id-1`
       xid=id-1
    END IF

    DO is=1,nsurf
       !  IF ( sfr(is) > .001 ) THEN
       !   call AnOHM to calculate the coefs.
       CALL AnOHM_coef(is,xid,Gridiv,MetForcingData_grid,moist_surf,AnthropHeatMethod,& !input
            alb, emis, cp, kk, ch,&! input
            xa1(is),xa2(is),xa3(is))                         ! output
       ! print*, 'AnOHM_coef are: ',xa1,xa2,xa3
       !  ELSE
       !     xa1(is  )=0.1
       !     xa2(is)=0.1
       !     xa3(is)=1
       !  END IF

    ENDDO

    !   calculate the areally-weighted OHM coefficients
    a1=DOT_PRODUCT(xa1,sfr)
    a2=DOT_PRODUCT(xa2,sfr)
    a3=DOT_PRODUCT(xa3,sfr)


    !   Calculate radiation part ------------------------------------------------------------
    qs=-999          !qs  = Net storage heat flux  [W m-2]
    IF(qn1>-999) THEN   !qn1 = Net all-wave radiation [W m-2]

       ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
       CALL OHM_dqndt_cal(nsh,qn1,qn1_store,qn1_av_store,dqndt)

       ! Calculate net storage heat flux
       CALL OHM_QS_cal(qn1,dqndt,a1,a2,a3,qs)

    ELSE
       CALL ErrorHint(21,'SUEWS_AnOHM.f95: bad value for qn found during qs calculation.',qn1,NotUsed,notUsedI)
    ENDIF

  END SUBROUTINE AnOHM
  !========================================================================================

  !========================================================================================
  !> High level wrapper for AnOHM coefficients calculation
  !! @brief
  !! calculate OHM coefficients within AnOHM framework.
  !! @returns
  !! -# OHM coefficients of a given surface type: a1, a2 and a3
  SUBROUTINE AnOHM_coef(&
       sfc_typ,xid,xgrid,&!input
       MetForcingData_grid,moist,AnthropHeatMethod,& !input
       alb, emis, cp, kk, ch,&! input
       xa1,xa2,xa3)                         ! output

    IMPLICIT NONE

    ! input
    INTEGER,INTENT(in):: sfc_typ           !< surface type [-]
    INTEGER,INTENT(in):: xid               !< day of year [-]
    INTEGER,INTENT(in):: xgrid             !< grid id [-]
    INTEGER,INTENT(in):: AnthropHeatMethod !< AnthropHeat option [-]

    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: alb                 !< albedo [-]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: emis                !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: cp                  !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: kk                  !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: ch                  !< bulk transfer coef [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)   :: moist               !< surface wetness status [-]
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:,:) :: MetForcingData_grid !< met forcing array of grid

    ! output
    REAL(KIND(1d0)),INTENT(out) :: xa1 !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out) :: xa2 !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out) :: xa3 !< AnOHM coefficients of grid [W m-2]

    ! surface temperature related scales:
    REAL(KIND(1d0)):: ATs   !< daily amplitude of surface temperature [K]
    REAL(KIND(1d0)):: mTs   !< daily mean of surface temperature [K]
    REAL(KIND(1d0)):: gamma !< phase difference between Ts and Sd [rad]

    !   forcing scales
    REAL(KIND(1d0))::ASd,mSd,tSd !< solar radiation
    REAL(KIND(1d0))::ATa,mTa,tTa !< air temperature
    REAL(KIND(1d0))::tau         !< phase lag between Sd and Ta (Ta-Sd)
    REAL(KIND(1d0))::mWS,mWF,mAH !< mean values of WS, WF and AH


    !   forcings:
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::Sd   ! incoming solar radiation [W m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::Ta   ! air temperature [degC]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::RH   ! relative humidity [%]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::pres ! air pressure [hPa]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::WS   ! wind speed [m s-1]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::WF   ! water flux density [m3 s-1 m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::AH   ! anthropogenic heat [W m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE::tHr  ! time [hr]

    !   sfc. properties:
    REAL(KIND(1d0)) ::xalb   ! albedo,
    REAL(KIND(1d0)) ::xemis  ! emissivity,
    REAL(KIND(1d0)) ::xcp    ! heat capacity,
    REAL(KIND(1d0)) ::xk     ! thermal conductivity,
    REAL(KIND(1d0)) ::xch    ! bulk transfer coef.
    REAL(KIND(1d0)) ::xBo    ! Bowen ratio
    REAL(KIND(1d0)) ::xeta   ! effective absorption coefficient
    REAL(KIND(1d0)) ::xmu    ! effective absorption fraction
    REAL(KIND(1d0)) ::xmoist ! surface wetness



    ! locally saved variables:
    ! if coefficients have been calculated, just reload them
    ! otherwise, do the calculation
    INTEGER,SAVE :: id_save,grid_save
    REAL(KIND(1d0)), SAVE:: coeff_grid_day(7,3)=-999.

    ! INTEGER :: WaterSurf=7

    ! PRINT*, 'xid,id_save',xid,id_save
    ! PRINT*, 'xgrid,grid_save',xgrid,grid_save
    ! PRINT*, 'sfc_typ',sfc_typ
    ! PRINT*, 'coeff_grid_day',coeff_grid_day(sfc_typ,:)
    IF ( xid==id_save .AND. xgrid ==grid_save) THEN
       ! if coefficients have been calculated, just reload them
       !  print*, 'here no repetition'
       xa1=coeff_grid_day(sfc_typ,1)
       xa2=coeff_grid_day(sfc_typ,2)
       xa3=coeff_grid_day(sfc_typ,3)
    ELSE
      !  PRINT*, ''
      !  PRINT*, 'surface:',sfc_typ
      !  PRINT*, 'xid',xid,id_save


       ! load forcing characteristics:
       CALL AnOHM_Fc(&
            xid,MetForcingData_grid,AnthropHeatMethod,& ! input
            ASd,mSd,tSd,ATa,mTa,tTa,tau,mWS,mWF,mAH)    ! output

       ! load forcing variables:
       CALL AnOHM_FcLoad(&
            xid,MetForcingData_grid,AnthropHeatMethod,& ! input
            Sd,Ta,RH,pres,WS,WF,AH,tHr)                 ! output

       ! load sfc. properties:
       xalb   = alb(sfc_typ)
       xemis  = emis(sfc_typ)
       xcp    = cp(sfc_typ)
       xk     = kk(sfc_typ)
       xch    = ch(sfc_typ)
       xmoist = moist(sfc_typ)

       !  PRINT*, 'xBo before:',xBo
       ! calculate Bowen ratio:
       CALL AnOHM_Bo_cal(&
            sfc_typ,&
            Sd,Ta,RH,pres,tHr,               & ! input: forcing
            ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH, & ! input: forcing
            xalb,xemis,xcp,xk,xch,xmoist,    & ! input: sfc properties
            tSd,                             & ! input: peaking time of Sd in hour
            xBo)                               ! output: Bowen ratio
       !  PRINT*, 'xBo after:',xBo

       !  calculate AnOHM coefficients
       SELECT CASE (sfc_typ)
       CASE (1:6) ! land surfaces
          CALL  AnOHM_coef_land_cal(&
               ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,& ! input: forcing
               xalb,xemis,xcp,xk,xch,xBo,      & ! input: sfc properties
               xa1,xa2,xa3,ATs,mTs,gamma)                    ! output: surface temperature related scales by AnOHM

       CASE (7) ! water surface
          ! NB:give fixed values for the moment
          xeta = 0.3
          xmu  = 0.2
          CALL  AnOHM_coef_water_cal(&
               ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
               xalb,xemis,xcp,xk,xch,xBo,xeta,xmu,&   ! input: sfc properties
               xa1,xa2,xa3,ATs,mTs,gamma)            ! output

          ! save variables for marking status as done
          id_save =xid
          grid_save=xgrid

       END SELECT
       ! save variables for reusing values of the same day
       coeff_grid_day(sfc_typ,:)=(/xa1,xa2,xa3/)
    END IF

  END SUBROUTINE AnOHM_coef
  !========================================================================================

  !========================================================================================
  !> calculate the surface temperature related parameters (ATs, mTs, gamma)
  !> based on forcings and sfc. conditions
  !> @return
  !> @b xTs surface temperature at local time
  SUBROUTINE AnOHM_xTs(&
       sfc_typ,& !input: surface type
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
       xalb,xemis,xcp,xk,xch,xBo,&   ! input: sfc properties
       tSd,& !input: peaking time of Sd in hour
       xTHr,&! input: time in hour
       xTs)! output: surface temperature

    IMPLICIT NONE
    ! input:
    INTEGER,INTENT(in):: sfc_typ !< surface type (land: 1-6, water: 7)

    ! input: forcing scales
    REAL(KIND(1d0)),INTENT(in):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: ATa !< daily amplitude of air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: mTa !< daily mean air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(in):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(in):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(in):: mAH !< daily mean anthropogenic heat flux [W m-2]

    ! input: sfc properties
    REAL(KIND(1d0)),INTENT(in):: xalb  !< albedo [-]
    REAL(KIND(1d0)),INTENT(in):: xemis !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in):: xcp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xk    !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in):: xch   !< bulk transfer coef [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xBo   !< Bowen ratio [-]
    REAL(KIND(1d0)):: xeta  !< effective absorption fraction [-]
    REAL(KIND(1d0)):: xmu   !< effective absorption coefficient [m-1]

    ! input: temporal-related
    REAL(KIND(1d0)),INTENT(in):: tSd  !< local peaking time of Sd, hour
    REAL(KIND(1d0)),INTENT(in):: xTHr !< local time to calculate Ts, hour

    ! output:
    REAL(KIND(1d0)),INTENT(out) :: xTs !< surface temperature at xTHr(hr)

    !   local
    REAL(KIND(1d0)) :: &
         xa1,xa2,xa3,&!coefficients
         ATs,mTs,gamma !surface temperature related scales by AnOHM

    ! constant:
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth

    SELECT CASE (sfc_typ)
    CASE (1:6) ! land surfaces
       CALL  AnOHM_coef_land_cal(&
            ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,& ! input: forcing
            xalb,xemis,xcp,xk,xch,xBo,      & ! input: sfc properties
            xa1,xa2,xa3,ATs,mTs,gamma)                    ! output: surface temperature related scales by AnOHM


    CASE (7) ! water surface
       ! !   NB:give fixed values for the moment
       xeta  = 0.3
       xmu   = 0.2
       CALL  AnOHM_coef_water_cal(&
            ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
            xalb,xemis,xcp,xk,xch,xBo,xeta,xmu,&   ! input: sfc properties
            xa1,xa2,xa3,ATs,mTs,gamma)            ! output

    END SELECT

    ! for a local time xTHr (in hour):
    xTs=ATs*SIN(OMEGA*(xTHr-tSd+6)*3600-gamma)+mTs

  END SUBROUTINE AnOHM_xTs
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_coef_land_cal(&
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
       xalb,xemis,xcp,xk,xch,xBo,&   ! input: sfc properties
       xa1,xa2,xa3,ATs,mTs,gamma)            ! output

    IMPLICIT NONE

    ! input: forcing scales
    REAL(KIND(1d0)),INTENT(in):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: ATa !< daily amplitude of air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: mTa !< daily mean air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(in):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(in):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(in):: mAH !< daily mean anthropogenic heat flux [W m-2]
    ! input: sfc properties
    REAL(KIND(1d0)),INTENT(in):: xalb  !< albedo [-]
    REAL(KIND(1d0)),INTENT(in):: xemis !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in):: xcp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xk    !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in):: xch   !< bulk transfer coef [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xBo   !< Bowen ratio [-]

    ! output
    REAL(KIND(1d0)),INTENT(out) :: xa1   !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out) :: xa2   !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out) :: xa3   !< AnOHM coefficients of grid [W m-2]
    REAL(KIND(1d0)),INTENT(out) :: ATs   !< daily amplitude of surface temperature [K]
    REAL(KIND(1d0)),INTENT(out) :: mTs   !< daily mean of surface temperature [K]
    REAL(KIND(1d0)),INTENT(out) :: gamma !< phase difference between Ts and Sd [K]


    !   constant
    REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth

    !   local variables:
    REAL(KIND(1d0)) :: beta              ! inverse Bowen ratio
    REAL(KIND(1d0)) :: f,fL,fT           ! energy redistribution factors
    REAL(KIND(1d0)) :: lambda            ! thermal diffusivity
    REAL(KIND(1d0)) :: delta,m,n         ! water flux related variables
    REAL(KIND(1d0)) :: ms,ns             ! m, n related
    REAL(KIND(1d0)) :: ceta,cphi         ! phase related temporary variables
    REAL(KIND(1d0)) :: eta,phi,xlag      ! phase related temporary variables
    REAL(KIND(1d0)) :: xx1,xx2,xx3,xchWS ! temporary use
    ! LOGICAL         :: flagGood = .TRUE. ! quality flag, T for good, F for bad



    !   give fixed values for test
    !   properties
    !   xalb  = .2
    !   xemis = .9
    !   xcp   = 1e6
    !   xk    = 1.2
    !   xch   = 4
    !   xBo    = 2
    !   forcings
    !   ASd = 400
    !   mSd = 100
    !   ATa = 23
    !   mTa = 23+C2K
    !   tau = PI/6
    !   WS  = 1
    !   AH  = 0
    !   Wf  = 0

    ! PRINT*, '********sfc_typ: ',sfc_typ,' start********'
    !   initialize flagGood as .TRUE.
    ! flagGood = .TRUE.
    ! PRINT*, flagGood

    !   calculate sfc properties related parameters:
    xchWS = xch*mWS
    ! xchWS = 0.5*xchWS

    !     xch    = CRA/RA !update bulk trsf. coeff. with RA (aerodyn. res.)
    beta   = 1/xBo
    ! PRINT*, 'beta:', beta
    f      = ((1+beta)*xchWS+4*SIGMA*xemis*mTa**3)
    !     print*, 'xch',xch,'mTa',mTa,'dLu',4*SIGMA*xemis*mTa**3
    fL     = 4*SIGMA*xemis*mTa**3
    !     print*, 'fL',fL
    fT     = (1+beta)*xchWS
    !     print*, 'fT',fT
    lambda = xk/xcp
    !     print*, 'lambda',lambda
    delta  = SQRT(.5*(mWF**2+SQRT(mWF**4+16*lambda**2*OMEGA**2)))
    !     print*, 'delta',delta
    m      = (2*lambda)/(delta+mWF)
    n      = delta/OMEGA
    !     print*, 'm',m,'n',n


    !   calculate surface temperature related parameters:
    mTs   = (mSd*(1-xalb)/f)+mTa
    ms    = 1+xk/(f*m)
    ns    = xk/(f*n)
    !     print*, 'ms,ns:', ms,ns
    xx1   = f*SIN(tau)*ATa
    !     print*, 'xx1',xx1
    xx2   = (1-xalb)*ASd+f*COS(tau)*ATa
    !     print*, 'xx2',xx2
    gamma = ATAN(ns/ms)+ATAN(xx1/xx2)
    !     print*, 'gamma:', gamma
    ATs   = -(SIN(tau)*ATa)/(ns*COS(gamma)-ms*SIN(gamma))
    !     print*,  'ATs:', ATs

    !   calculate net radiation related parameters:
    xx1  = (ns*COS(gamma)+SIN(gamma)-ms*SIN(gamma))*SIN(tau)*ATa*fL
    xx2  = (xalb-1)*(ns*COS(gamma)-ms*SIN(gamma))*ASd
    xx3  = (-ms*COS(tau)*SIN(tau)+COS(gamma)*(ns*COS(tau)+SIN(tau)))*ATa*fL
    xx2  = xx2-xx3
    phi  = ATAN(xx1/xx2)
    xx3  = (ns*COS(gamma)-ms*SIN(gamma))
    xx1  = (1+SIN(gamma)/xx3)*SIN(tau)*ATa*fL
    xx2  = (xalb-1)*ASd-(COS(tau)+COS(gamma)*SIN(tau)/xx3)*ATa*fL
    cphi = SQRT(xx1**2+xx2**2)

    !   calculate heat storage related parameters:
    xx1  = m*COS(gamma)-n*SIN(gamma)
    xx2  = m*SIN(gamma)+n*COS(gamma)
    eta  = ATAN(xx1/xx2)
    !     if ( eta<0 ) print*, 'lambda,delta,m,n,gamma:', lambda,delta,m,n,gamma
    xx1  = xk**2*(m**2+n**2)*ATs**2
    xx2  = m**2*n**2
    ceta = SQRT(xx1/xx2)

    !   key phase lag:
    xlag = eta-phi

    !   calculate the OHM coeffs.:
    !   a1:
    xa1 = (ceta*COS(xlag))/cphi

    !   a2:
    xa2 = (ceta*SIN(xlag))/(OMEGA*cphi)
    xa2 = xa2/3600 ! convert the unit from s-1 to h-1

    !   a3:
    xa3  = -xa1*(fT/f)*(mSd*(1-xalb))-mAH

    ! !   quality checking:
    ! !   quality checking of forcing conditions
    ! ! IF ( ASd < 0 .OR. ATa < 0 .OR. ATs < 0 .OR. tau<-4.0/12*Pi) flagGood = .FALSE.
    ! !   quality checking of a1
    ! IF ( .NOT. (xa1>0 .AND. xa1<0.7)) THEN
    !   !  flagGood = .FALSE.
    !    IF (xa1 >0.7) xa1=MAX(0.7,xa1)
    ! ENDIF
    ! !   quality checking of a2
    ! IF ( .NOT. (xa2>-0.5 .AND. xa2<0.5)) THEN
    !   !  flagGood = .FALSE.
    !    !  IF ( xa2>0.5) xa2 = 0.5
    !    !  IF (xa2<-0.5) xa2 = -0.5
    ! ENDIF
    ! !   quality checking of a3
    ! ! IF ( .NOT. (xa3<0)) flagGood = .FALSE.

    !     print*, 'ceta,cphi', ceta,cphi
    !     print*, 'tau,eta,phi,xlag in deg:',tau/pi*180,eta/pi*180,phi/pi*180,xlag/pi*180

    ! PRINT*, '********sfc_typ: ',sfc_typ,' end********'

  END SUBROUTINE AnOHM_coef_land_cal
  !========================================================================================

  !========================================================================================
  !> a wrapper for retrieving AnOHM coefficients of water body
  !> @returns
  !! -# OHM coefficients of a given surface type: a1, a2 and a3
  !> @Caution: NB: this SUBROUTINE hasn't been well tested.
  SUBROUTINE AnOHM_coef_water_cal(&
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
       xalb,xemis,xcp,xk,xch,xBo,xeta,xmu,&   ! input: sfc properties
       xa1,xa2,xa3,ATs,mTs,gamma)            ! output


    IMPLICIT NONE

    ! input: forcing scales
    REAL(KIND(1d0)),INTENT(in):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: ATa !< daily amplitude of air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: mTa !< daily mean air temperature [K]
    REAL(KIND(1d0)),INTENT(in):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(in):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(in):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(in):: mAH !< daily mean anthropogenic heat flux [W m-2]

    ! input: sfc properties
    REAL(KIND(1d0)),INTENT(in):: xalb  !< albedo [-]
    REAL(KIND(1d0)),INTENT(in):: xemis !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in):: xcp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xk    !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in):: xch   !< bulk transfer coef [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in):: xBo   !< Bowen ratio [-]
    REAL(KIND(1d0)),INTENT(in):: xeta  !< effective absorption fraction [-]
    REAL(KIND(1d0)),INTENT(in):: xmu   !< effective absorption coefficient [m-1]

    ! output
    REAL(KIND(1d0)),INTENT(out) :: xa1   !< AnOHM coefficients of grid [-]
    REAL(KIND(1d0)),INTENT(out) :: xa2   !< AnOHM coefficients of grid [h]
    REAL(KIND(1d0)),INTENT(out) :: xa3   !< AnOHM coefficients of grid [W m-2]
    REAL(KIND(1d0)),INTENT(out) :: ATs   !< daily amplitude of surface temperature [K]
    REAL(KIND(1d0)),INTENT(out) :: mTs   !< daily mean of surface temperature [K]
    REAL(KIND(1d0)),INTENT(out) :: gamma !< phase difference between Ts and Sd [K]

    !   constant
    REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth

    !   local variables:
    REAL(KIND(1d0)) :: beta                   ! inverse Bowen ratio
    REAL(KIND(1d0)) :: f,fL,fT                ! energy redistribution factors
    REAL(KIND(1d0)) :: lambda,calb            ! temporary use
    REAL(KIND(1d0)) :: delta             ! sfc. temperature related variables
    ! REAL(KIND(1d0)) :: delta,m,n              ! sfc. temperature related variables
    REAL(KIND(1d0)) :: xm,xn                  ! m, n related
    ! REAL(KIND(1d0)) :: gamma              ! phase lag scale
    REAL(KIND(1d0)) :: phi              ! phase lag scale
    ! REAL(KIND(1d0)) :: ATs,mTs                ! surface temperature amplitude
    REAL(KIND(1d0)) :: czeta,ctheta           ! phase related temporary variables
    REAL(KIND(1d0)) :: zeta,theta,xlag           ! phase related temporary variables
    REAL(KIND(1d0)) :: xx1,xx2,xx3            ! temporary use
    REAL(KIND(1d0)) :: kappa                  ! temporary use
    REAL(KIND(1d0)) :: dtau,dpsi,dphi         ! temporary use
    REAL(KIND(1d0)) :: cdtau,cdpsi,cdphi      ! temporary use
    REAL(KIND(1d0)) :: xxT,xxkappa,xxdltphi,xchWS   ! temporary use
    ! LOGICAL :: flagGood = .TRUE.  ! quality flag, T for good, F for bad

    ! ====not used====
    REAL(KIND(1d0)) :: dummy
    dummy=mah+mwf
    ! ====not used====

    !   calculate sfc properties related parameters:
    xm     = xk*xmu**2
    xn     = xcp*OMEGA
    phi    = ATAN(xn/xm)
    kappa  = SQRT(xcp*OMEGA/(2*xk))

    ! mWS=SUM(WS, dim=1, mask=(WS>0))/SIZE(WS, dim=1)
    xchWS    = xch*mWS
    beta   = 1/xBo
    f      = ((1+beta)*xchWS+4*SIGMA*xemis*mTa**3)
    fL     = 4*SIGMA*xemis*mTa**3
    fT     = (1+beta)*xchWS

    calb=1-xalb


    lambda= SQRT(xm**2+xn**2)

    dtau=ATAN(xk*kappa/(f+xk*kappa))
    dpsi=ATAN((xk-xmu)/kappa)
    dphi=ATAN(kappa*(f+xk*xmu)/(f*(kappa-xmu)+xk*kappa*(2*kappa-xmu)))
    cdtau=SQRT((xk*kappa)**2+(f+xk*kappa)**2)
    cdpsi=SQRT((xk-xmu)**2+kappa**2)
    cdphi=cdtau*cdpsi

    !   calculate surface temperature related parameters:
    ! daily mean:
    mTs   = (mSd*(1-xalb+xeta)/f)+mTa
    ! amplitude:
    xx1 = (xk*xeta*xmu*calb*ASd*cdpsi)**2
    xx2 = 2*lambda*SQRT(xx1)*(calb*ASd*SIN(phi-dpsi)+f*ATa*SIN(tau+phi-dpsi))
    xx3 = lambda**2*((calb*ASd+COS(tau)*f*ATa)**2+(SIN(tau)*f*ATa)**2)
    ATs = 1/(cdtau*lambda)*SQRT(xx1+xx2+xx3)
    ! phase lag:
    xx1 = (xk*kappa*calb*ASd+cdtau*f*ATa*SIN(tau+dtau))*lambda    &
         +(xk*xeta*xmu*calb*ASd*cdphi*SIN(phi+dphi))
    xx2 = ((f+xk*kappa)*calb*ASd-cdtau*f*ATa*COS(tau+dtau))*lambda&
         -xk*xeta*xmu*calb*ASd*cdphi*COS(phi+dphi)
    delta = ATAN(xx1/xx2)
    gamma = delta

    ! calculate net radiation related parameters:
    ! phase lag:
    xx1    = fL*(ATs*SIN(delta)-ATa*SIN(tau))
    xx2    = calb*ASd-fL*(ATs*COS(delta)-ATa*COS(tau))
    theta  = ATAN(xx1/xx2)
    ! amplitude:
    ctheta = SQRT(xx1**2+xx2**2)

    !   calculate heat storage related parameters:
    ! scales:
    xxT      = SQRT(2.)*kappa*lambda*ATs
    xxkappa  = cdpsi*xeta*xmu*ASd
    xxdltphi = COS(delta)*SIN(dpsi)*COS(phi)-SIN(delta)*COS(dpsi)*SIN(phi)
    ! phase lag:
    xx1  = xxT*SIN(PI/4-delta)+xxkappa*SIN(phi+dpsi)
    xx2  = xxT*SIN(PI/4+delta)-xxkappa*SIN(PHI-dpsi)
    zeta = ATAN(xx1/xx2)
    ! amplitude:
    xx1   = 2*SQRT(2.)*xxkappa*xxT*xxdltphi
    xx2   = (1-COS(2*dpsi)*COS(2*phi))*xxkappa**2
    xx3   = xxT**2
    czeta = xk/lambda*SQRT(xx1+xx2+xx3)


    !   calculate the OHM coeffs.:
    xlag = zeta-theta
    ! a1:
    xa1  = (czeta*COS(xlag))/ctheta

    !   write(*,*) 'ceta,xlag,cphi:', ceta,xlag,cphi
    !   a2:
    xa2  = (czeta*SIN(xlag))/(OMEGA*ctheta)
    xa2  = xa2/3600 ! convert the unit from s-1 to h-1

    !   a3:
    xa3  = mSd*(xalb-1)*(xeta+(fT-fL*xeta)/f*xa1)


    !   quality checking:
    !   quality checking of forcing conditions
    ! IF ( ASd < 0 .OR. ATa < 0 .OR. ATs < 0 .OR. tau<-4.0/12*Pi) flagGood = .FALSE.
    ! !   quality checking of a1
    ! IF ( .NOT. (xa1>0 .AND. xa1<0.7)) THEN
    !    flagGood = .FALSE.
    !    IF (xa1 >0.7) xa1=MAX(0.7,xa1)
    ! ENDIF
    ! !   quality checking of a2
    ! IF ( .NOT. (xa2>-0.5 .AND. xa2<0.5)) THEN
    !    flagGood = .FALSE.
    !    !  IF ( xa2>0.5) xa2 = 0.5
    !    !  IF (xa2<-0.5) xa2 = -0.5
    ! ENDIF
    ! !   quality checking of a3
    ! IF ( .NOT. (xa3<0)) flagGood = .FALSE.

    !   skip the first day for quality checking
    ! IF ( xid == 1 ) flagGood = .TRUE.

    ! print*,  'sfc_typ_water:', sfc_typ
    ! print*, 'a1,a2,a3:', xa1,xa2,xa3

  END SUBROUTINE AnOHM_coef_water_cal
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_Fc(&
       xid,MetForcingData_grid,AnthropHeatMethod,& ! input
       ASd,mSd,tSd,ATa,mTa,tTa,tau,mWS,mWF,mAH)    ! output

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in):: xid
    INTEGER,INTENT(in):: AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:,:) ::MetForcingData_grid

    !   output
    REAL(KIND(1d0)),INTENT(out):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(out):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(out):: tSd !< local peaking time of solar radiation [hr]
    REAL(KIND(1d0)),INTENT(out):: ATa !< daily amplitude of air temperature [degC]
    REAL(KIND(1d0)),INTENT(out):: mTa !< daily mean air temperature [degC]
    REAL(KIND(1d0)),INTENT(out):: tTa !< local peaking time of air temperature [hour]
    REAL(KIND(1d0)),INTENT(out):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(out):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(out):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(out):: mAH !< daily mean anthropogenic heat flux [W m-2]

    !   forcings:
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: Sd   !< incoming solar radiation [W m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: Ta   !< air temperature [degC]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: RH   !< relative humidity [%]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: pres !< Atmospheric pressure [hPa]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: WS   ! wind speed [m s-1]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: WF   ! water flux density [m3 s-1 m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: AH   ! anthropogenic heat [W m-2]
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: tHr  ! local time [hr]


    ! load forcing variables:
    CALL AnOHM_FcLoad(xid,MetForcingData_grid,AnthropHeatMethod,Sd,Ta,RH,pres,WS,WF,AH,tHr)
    ! calculate forcing scales for AnOHM:
    CALL AnOHM_FcCal(Sd,Ta,WS,WF,AH,tHr,ASd,mSd,tSd,ATa,mTa,tTa,tau,mWS,mWF,mAH)

    ! CALL r8vec_print(SIZE(sd, dim=1),sd,'Sd')
    ! PRINT*, ASd,mSd,tSd

    ! CALL r8vec_print(SIZE(ta, dim=1),ta,'Ta')
    ! PRINT*, ATa,mTa,tTa

  END SUBROUTINE AnOHM_Fc
  !========================================================================================

  !========================================================================================
  !> load forcing series for AnOHM_FcCal
  SUBROUTINE AnOHM_FcLoad(&
       xid,MetForcingData_grid,AnthropHeatMethod,& ! input
       Sd,Ta,RH,pres,WS,WF,AH,tHr) ! output

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in):: xid !< day of year
    INTEGER,INTENT(in):: AnthropHeatMethod !< AnthropHeat option
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:,:) ::MetForcingData_grid !< met forcing array of grid

    !   output
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: Sd  !< incoming solar radiation [W m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: Ta  !< air temperature [degC]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: RH  !< relative humidity [%]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: pres!< atmospheric pressure [mbar]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: WS  !< wind speed [m s-1]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: WF  !< water flux density [m3 s-1 m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: AH  !< anthropogenic heat [W m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: tHr !< local time  [hr]


    !   local variables:
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE :: subMet ! subset array of daytime series

    INTEGER :: err
    INTEGER :: lenMetData,nVar


    LOGICAL, ALLOCATABLE :: metMask(:)

    ! IF ( xid==73 ) THEN
    !    PRINT*, 'MetForcingData_grid in AnOHM_FcLoad, shape:',SHAPE(MetForcingData_grid)
    !    PRINT*, MetForcingData_grid(UBOUND(MetForcingData_grid, dim=1),:)
    ! END IF


    ! determine the length of subset
    lenMetData = COUNT(&
         MetForcingData_grid(:,2)==xid & ! day=xid
         .AND. MetForcingData_grid(:,4)==0)! tmin=0

    ! construct mask
    IF (ALLOCATED(metMask)) DEALLOCATE(metMask, stat=err)
    ALLOCATE(metMask(lenMetData))
    metMask=(MetForcingData_grid(:,2)==xid & ! day=xid
         .AND. MetForcingData_grid(:,4)==0)! tmin=0

    ! construct array for time and met variables
    nVar=8! number of variables to retrieve
    ! print*, 'good 1'
    ! allocate subMet:
    IF (ALLOCATED(subMet)) DEALLOCATE(subMet, stat=err)
    ALLOCATE(subMet(lenMetData,nVar))
    subMet=RESHAPE(PACK(MetForcingData_grid(:,(/3,& !time: hour
         15,12,11,13,10,12,9/)),&! met: Sd, Ta, RH, pres, WS, WF, AH
                                ! NB: WF not used: a placeholder
         SPREAD(metMask, dim=2, ncopies=nVar)),& ! replicate mask vector to 2D array
         (/lenMetData,nVar/)) ! convert to target shape


    ! re-allocate arrays as their sizes may change during passing
    IF (ALLOCATED(tHr)) DEALLOCATE(tHr, stat=err)
    ALLOCATE(tHr(lenMetData))
    IF (ALLOCATED(Sd)) DEALLOCATE(Sd, stat=err)
    ALLOCATE(Sd(lenMetData))
    IF (ALLOCATED(Ta)) DEALLOCATE(Ta, stat=err)
    ALLOCATE(Ta(lenMetData))
    IF (ALLOCATED(RH)) DEALLOCATE(RH, stat=err)
    ALLOCATE(RH(lenMetData))
    IF (ALLOCATED(pres)) DEALLOCATE(pres, stat=err)
    ALLOCATE(pres(lenMetData))
    IF (ALLOCATED(WS)) DEALLOCATE(WS, stat=err)
    ALLOCATE(WS(lenMetData))
    IF (ALLOCATED(WF)) DEALLOCATE(WF, stat=err)
    ALLOCATE(WF(lenMetData))
    IF (ALLOCATED(AH)) DEALLOCATE(AH, stat=err)
    ALLOCATE(AH(lenMetData))

    ! load the sublist into forcing variables
    tHr  = subMet(:,1)! time in hour
    Sd   = subMet(:,2)
    Ta   = subMet(:,3)
    RH   = subMet(:,4)
    pres = subMet(:,5)
    WS   = subMet(:,6)
    WF   = 0          ! set as 0 for the moment
    IF ( AnthropHeatMethod == 0 ) THEN
       AH = subMet(:,8)    ! read in from MetForcingData_grid,
    ELSE
       AH = 0 ! temporarily change to zero; TODO: chenge back to modelled value
       !  AH = mAH_grids(xid-1,xgrid)
    END IF

    ! IF ( xid==73 ) THEN
    !    PRINT*, 'in  AnOHM_FcLoad'
    !    PRINT*, 'id:',xid
    !    PRINT*, 'Sd:', Sd
    !    PRINT*, 'Ta:', Ta
    !    PRINT*, 'RH:', RH
    !    PRINT*, 'pres:', pres
    !    PRINT*, 'WS:', WS
    !    PRINT*, 'WF:', WF
    !    PRINT*, 'AH:', AH
    !    PRINT*, 'tHr:', tHr
    ! END IF

  END SUBROUTINE AnOHM_FcLoad
  !========================================================================================

  !========================================================================================
  !> calculate the key parameters of a sinusoidal curve for AnOHM forcings
  !> i.e., a, b, c in a*Sin(Pi/12*t+b)+c
  SUBROUTINE AnOHM_FcCal(&
       Sd,Ta,WS,WF,AH,tHr,&                     ! input
       ASd,mSd,tSd,ATa,mTa,tTa,tau,mWS,mWF,mAH) ! output
    IMPLICIT NONE

    ! input
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: Sd  !< incoming shortwave radiation [W m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: Ta  !< air temperature [degC]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: WS  !< wind speed [m s-1]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: WF  !< water flux density [m3 s-1 m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: AH  !< anthropogenic heat [W m-2]
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: tHr !< time [hr]

    ! output
    REAL(KIND(1d0)),INTENT(out):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(out):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(out):: tSd !< local peaking time of solar radiation [hr]
    REAL(KIND(1d0)),INTENT(out):: ATa !< daily amplitude of air temperature [degC]
    REAL(KIND(1d0)),INTENT(out):: mTa !< daily mean air temperature [degC]
    REAL(KIND(1d0)),INTENT(out):: tTa !< local peaking time of air temperature [hr]
    REAL(KIND(1d0)),INTENT(out):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(out):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(out):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(out):: mAH !< daily mean anthropogenic heat flux [W m-2]

    !   constant
    REAL(KIND(1d0)), PARAMETER :: PI  = ATAN(1.0)*4 ! Pi
    REAL(KIND(1d0)), PARAMETER :: C2K = 273.15      ! degC to K

    !   local variables:
    REAL(KIND(1d0)),ALLOCATABLE :: tHrDay(:) ! daytime tHr when Sd>0
    REAL(KIND(1d0)),ALLOCATABLE :: selX(:)   ! daytime sublist of met variable when Sd>0


    ! REAL(KIND(1d0)) :: xx              ! temporary use
    INTEGER :: err,lenDay
    LOGICAL,DIMENSION(:),ALLOCATABLE::SdMask


    lenDay=COUNT(Sd>5, dim=1)
    ALLOCATE(SdMask(lenDay), stat=err)
    IF ( err/= 0) PRINT *, "SdMask: Allocation request denied"
    SdMask=Sd>5

    ! CALL r8vec_print(24,Sd,'Sd')
    ! CALL r8vec_print(24,tHr,'tHr')
    ALLOCATE(tHrDay(lenDay), stat=err)
    IF ( err/= 0) PRINT *, "tHrDay: Allocation request denied"
    tHrDay=PACK(tHr, mask=SdMask)


    ALLOCATE(selX(lenDay), stat=err)
    IF ( err/= 0) PRINT *, "selX: Allocation request denied"

    !   calculate sinusoidal scales of Sd:
    ! PRINT*, 'Calc. Sd...'
    selX=PACK(Sd, mask=SdMask)
    ASd=(MAXVAL(selX)-MINVAL(selX))/2
    mSd=SUM(selX)/lenDay
    tSd=12
    CALL AnOHM_ShapeFit(tHrDay,selX,ASd,mSd,tSd)
    ! CALL r8vec_print(lenDay,selX,'Sd Day:')
    !   modify ill-shaped days to go through
    ! IF ( ASd < 0 .OR. tSd > 15) THEN
    !    !         ASd = abs(ASd)
    !    !         tSd = 12 ! assume Sd peaks at 12:00LST
    !    CALL r8vec_print(lenDay,tHrDay,'tHrDay:')
    !    CALL r8vec_print(lenDay,selX,'Sd Day:')
    !    PRINT*, 'ASd:', ASd
    !    PRINT*, 'mSd:', mSd
    !    PRINT*, 'tSd:', tSd
    ! END IF
    ! PRINT*, 'Sd Day:', selX

    !   calculate sinusoidal scales of Ta:
    ! PRINT*, 'Calc. Ta...'
    selX=PACK(Ta, mask=SdMask)
    ATa=(MAXVAL(selX)-MINVAL(selX))/2
    mTa=SUM(selX)/lenDay
    tTa=12
    CALL AnOHM_ShapeFit(tHrDay,selX,ATa,mTa,tTa)
    ! CALL r8vec_print(lenDay,selX,'Ta Day:')
    IF ( mTa < 60 ) mTa = mTa+C2K ! correct the Celsius to Kelvin
    !   modify ill-shaped days to go through
    IF ( ATa < 0 ) THEN
       !         ATa = abs(ATa)
       !         tTa = 14 ! assume Ta peaks at 14:00LST
       CALL r8vec_print(lenDay,selX,'Ta Day:')
       PRINT*, 'ATa:', ATa
       PRINT*, 'mTa:', mTa
       PRINT*, 'tTa:', tTa
    END IF
    ! PRINT*, 'Ta:', Ta(10:16)


    !   calculate the phase lag between Sd and Ta:
    tau = (tTa-tSd)/24*2*PI
    ! PRINT*, 'tau:', tau

    !   calculate the mean values:
    selX=PACK(WS, mask=SdMask)
    mWS = SUM(selX)/lenDay  ! mean value of WS

    selX=PACK(WF, mask=SdMask)
    mWF = SUM(selX)/lenDay  ! mean value of WF

    selX=PACK(AH, mask=SdMask)
    mAH = SUM(selX)/lenDay  ! mean value of AH
    ! PRINT*, 'mWS:', mWS
    !     print*, 'mWF:', mWF
    !     print*, 'mAH:', mAH
    IF (ALLOCATED(SdMask)) DEALLOCATE(SdMask, stat=err)
    IF ( err/= 0) PRINT *, "SdMask: Deallocation request denied"

    IF (ALLOCATED(tHrDay)) DEALLOCATE(tHrDay, stat=err)
    IF ( err/= 0) PRINT *, "tHrDay: Deallocation request denied"

    IF (ALLOCATED(selX)) DEALLOCATE(selX, stat=err)
    IF ( err/= 0) PRINT *, "selX: Deallocation request denied"

  END SUBROUTINE AnOHM_FcCal
  !========================================================================================

  !========================================================================================
  !> calculate the key parameters of a sinusoidal curve for AnOHM forcings
  !> i.e., a, b, c in a*Sin(Pi/12*t+b)+c, where t is in hour
  SUBROUTINE AnOHM_ShapeFit(&
       tHr,obs,&      !input
       amp,mean,tpeak)!output
    IMPLICIT NONE

    !   input
    REAL(KIND(1d0)),INTENT(in) :: tHr(:)  !< time in hour
    REAL(KIND(1d0)),INTENT(in) :: obs(:)  !< observation

    !   output
    REAL(KIND(1d0)),INTENT(out) :: amp     !< amplitude
    REAL(KIND(1d0)),INTENT(out) :: mean    !< average
    REAL(KIND(1d0)),INTENT(out) :: tpeak   !< peaking time (h)

    INTEGER :: m,n,info,err         ! temporary use

    ! EXTERNAL fSin
    REAL ( KIND(1d0) ),ALLOCATABLE:: fvec(:),x(:)

    REAL ( KIND(1d0) ):: tol= 0.00001D+00 ! tolerance

    n=3 ! number of parameters to determine
    m=SIZE(tHr, dim=1) ! number of observation pairs
    ! PRINT*, 'm',m


    ALLOCATE(fvec(m), stat=err)
    IF ( err/= 0) PRINT *, "fvec: Allocation request denied"

    ALLOCATE(x(n), stat=err)
    IF ( err/= 0) PRINT *, "x: Allocation request denied"


    ! initial guess for fitting
    x=(/mean,amp,tpeak/)

    ! PRINT*, 'x',x
    ! iflag=1
    ! CALL fSin(m,n,x,tHr,obs,fvec,iflag)
    ! CALL r8vec_print(3,x,'x')

    ! use minpack subroutine 'lmstr1' to fit the sinusoidal form
    CALL lmdif1( fSin, m, n, x, tHr, obs, fvec, tol, info )

    ! x   = (/mean,amp,delta/), see subroutine 'fSin' for the meaning of this vector
    mean  = x(1)
    amp   = x(2)
    tpeak = x(3)+6 ! when t = delta + 6 (hr of LST), t is tpeak.

    ! adjust fitted parameters to physical range:
    ! amp>0
    IF ( amp<0 ) THEN
       !  PRINT*, ''
       !  PRINT*, 'before:'
       !  PRINT*, amp,tpeak,mean
       !  CALL r8vec_print(SIZE(tHR, dim=1),tHR,'tHR')
       !  CALL r8vec_print(SIZE(obs, dim=1),obs,'obs')
       amp=ABS(amp)
       tpeak=x(3)-12+6+24
       tpeak=MOD(tpeak,24.)
       !  PRINT*, 'after:'
       !  PRINT*, amp,tpeak,mean
    END IF
    tpeak=MOD(tpeak,24.)

  END SUBROUTINE AnOHM_ShapeFit
  !========================================================================================

  !========================================================================================
  !> sinusoidal function f(t) for fitting:
  !> f(t) = mean+amp*Sin(Pi/12(t-delta))
  !> x    = (/mean,amp,delta/) contains the fitting parameters
  SUBROUTINE fSin(m, n, x, xdat, ydat, fvec, iflag)

    IMPLICIT NONE
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    ! REAL ( kind = 8 ) fjrow(n)
    REAL ( kind = 8 ) fvec(m),& ! residual vector
         xdat(m),ydat(m) ! (x,y) observations for fitting

    INTEGER ( kind = 4 ) iflag,i
    REAL ( kind = 8 ) x(n)
    REAL(KIND(1d0)), PARAMETER :: PI  = ATAN(1.0)*4      ! Pi

    IF ( iflag == 0 ) THEN

       WRITE ( *, '(a)' ) ''
       DO i = 1, n
          WRITE ( *, '(g14.6)' ) x(i)
       END DO

    ELSE  IF ( iflag == 1 ) THEN
       fvec(1:m) = x(1) + x(2) * SIN(2*PI/24*(xdat(1:m)-x(3))) - ydat(1:m)

    END IF
    RETURN

  END SUBROUTINE fSin
  !========================================================================================

  !========================================================================================
  !> estimate daytime Bowen ratio for calculation of AnOHM coefficients
  SUBROUTINE AnOHM_Bo_cal(&
       sfc_typ,& ! surface type
       Sd,Ta,RH,pres,tHr,              & ! input: forcing
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,& ! input: forcing
       xalb,xemis,xcp,xk,xch,xSM,      & ! input: sfc properties
       tSd,                            & ! input: peaking time of Sd in hour
       xBo)                              ! output: Bowen ratio


    IMPLICIT NONE

    ! input:
    INTEGER, INTENT(in) :: sfc_typ ! unknown Bowen ratio

    ! input: daytime series
    REAL(kind = 8), INTENT(in), DIMENSION(:) ::Sd  !< incoming solar radiation [W m-2]
    REAL(kind = 8), INTENT(in), DIMENSION(:) ::Ta  !< air temperature [degC]
    REAL(kind = 8), INTENT(in), DIMENSION(:) ::RH  !< relative humidity [%]
    REAL(kind = 8), INTENT(in), DIMENSION(:) ::pres!< Atmospheric pressure [mbar]
    REAL(kind = 8), INTENT(in), DIMENSION(:) ::tHr !< local time  [hr]

    ! input: forcing scales
    REAL(KIND(1d0)),INTENT(in):: ASd !< daily amplitude of solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: mSd !< daily mean solar radiation [W m-2]
    REAL(KIND(1d0)),INTENT(in):: tSd !< local peaking time of solar radiation [hr]
    REAL(KIND(1d0)),INTENT(in):: ATa !< daily amplitude of air temperature [degC]
    REAL(KIND(1d0)),INTENT(in):: mTa !< daily mean air temperature [degC]
    REAL(KIND(1d0)),INTENT(in):: tau !< phase lag between Sd and Ta (Ta-Sd) [rad]
    REAL(KIND(1d0)),INTENT(in):: mWS !< daily mean wind speed [m s-1]
    REAL(KIND(1d0)),INTENT(in):: mWF !< daily mean underground moisture flux [m3 s-1 m-2]
    REAL(KIND(1d0)),INTENT(in):: mAH !< daily mean anthropogenic heat flux [W m-2]

    ! input: surface properties
    REAL(KIND(1d0)),INTENT(in) :: xalb  !< albedo [-]
    REAL(KIND(1d0)),INTENT(in) :: xemis !< emissivity [-]
    REAL(KIND(1d0)),INTENT(in) :: xcp   !< heat capacity [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in) :: xk    !< thermal conductivity [W m-1 K-1]
    REAL(KIND(1d0)),INTENT(in) :: xch   !< bulk transfer coef [J m-3 K-1]
    REAL(KIND(1d0)),INTENT(in) :: xSM   !< surface moisture status [-]

    ! output:
    REAL(kind=8), INTENT(out) :: xBo ! unknown Bowen ratio [-]

    ! EXTERNAL fcnBo

    REAL(kind=8), ALLOCATABLE :: x(:),fvec(:),prms(:)

    INTEGER :: lenDay,n,m,info,err,nVar,nPrm
    LOGICAL, DIMENSION(:), ALLOCATABLE :: maskDay
    REAL(kind=8) :: tol=1E-20

    ! LOGICAL, ALLOCATABLE,DIMENSION(:) :: metMask

    ! daytime mask
    ALLOCATE(maskDay(SIZE(sd)),stat=err)
    maskDay=sd>5
    ! length of daytime series
    lenDay=SIZE(PACK(sd, mask=maskDay), dim=1)
    ! ! daytime mask
    ! IF (ALLOCATED(metMask)) DEALLOCATE(metMask, stat=err)
    ! ALLOCATE(metMask(lenDay))
    ! metMask=(Sd>0)


    ! assign x vector:
    n=1
    ALLOCATE(x(n), stat=err)
    IF ( err/= 0) PRINT *, "x: Allocation request denied"
    ALLOCATE(fvec(n), stat=err)
    IF ( err/= 0) PRINT *, "fvec: Allocation request denied"

    ! pass initial Bowen ratio:
    xBo=10
    x(1)= xBo

    ! NB: set these numbers properly if any changes made to this subroutine:
    ! number of parameters to pass:
    nPrm=16
    ! number of met variables to pass:
    nVar=5
    ! length of x vector that holds unknown and parameters
    m=nPrm+nVar*lenDay

    ALLOCATE(prms(m), stat=err)
    IF ( err/= 0) PRINT *, "prms: Allocation request denied"

    ! pass forcing scales:
    prms(1)=ASd
    prms(2)=mSd
    prms(3)=ATa
    prms(4)=mTa
    prms(5)=tau
    prms(6)=mWS
    prms(7)=mWF
    prms(8)=mAH

    ! pass sfc. property scales:
    prms(9)=xalb
    prms(10)=xemis
    prms(11)=xcp
    prms(12)=xk
    prms(13)=xch
    prms(14)=xSM

    ! pass tSd:
    prms(15)=tSd

    ! pass tSd:
    prms(16)=sfc_typ*1.0


    ! extract daytime series
    prms(nPrm+1:m)=PACK((/Sd,Ta,RH,pres,tHr/), &
         mask=PACK(SPREAD(maskDay, dim=2, ncopies=nVar),.TRUE.))

    ! PRINT*, 'xBo before solve:',x(1)
    ! PRINT*, 'fvec before solve:',fvec(1)
    ! PRINT*, 'xSM:',xSM
    ! solve nonlinear equation fcnBo(x)=0
    CALL hybrd1(fcnBo,n,x,fvec,tol,info,m,prms)
    xBo=x(1)
    ! PRINT*, 'xBo after solve: ',x(1)
    ! PRINT*, 'fvec after solve:',fvec(1)

    IF (ALLOCATED(x)) DEALLOCATE(x, stat=err)
    IF ( err/= 0) PRINT *, "x: Deallocation request denied"
    IF (ALLOCATED(fvec)) DEALLOCATE(fvec, stat=err)
    IF ( err/= 0) PRINT *, "fvec: Deallocation request denied"
    IF (ALLOCATED(prms)) DEALLOCATE(prms, stat=err)
    IF ( err/= 0) PRINT *, "prms: Deallocation request denied"

  END SUBROUTINE AnOHM_Bo_cal
  !========================================================================================

  !========================================================================================
  !> this fucntion will construct an equaiton for Bo calculation
  SUBROUTINE fcnBo( n, x, fvec, iflag, m, prms )

    IMPLICIT NONE
    !    Input, external FCN, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !      subroutine fcn ( n, x, fvec, iflag )
    INTEGER ( kind = 4 ):: n
    INTEGER ( kind = 4 ):: m
    INTEGER ( kind = 4 ):: iflag
    REAL ( kind = 8 ):: fvec(n)
    REAL ( kind = 8 ):: x(n) ! x(1) as unknown Bo
    REAL ( kind = 8 ):: prms(m) ! prms(i) used for passing parameters

    ! the unknow: Bowen ratio
    REAL(kind = 8):: xBo

    ! forcing scales
    REAL(kind = 8):: ASd
    REAL(kind = 8):: mSd
    REAL(kind = 8):: ATa
    REAL(kind = 8):: mTa
    REAL(kind = 8):: tau
    REAL(kind = 8):: mWS
    REAL(kind = 8):: mWF
    REAL(kind = 8):: mAH

    ! sfc. property scales
    REAL(kind = 8):: xalb
    REAL(kind = 8):: xemis
    REAL(kind = 8):: xcp
    REAL(kind = 8):: xk
    REAL(kind = 8):: xch

    ! peaking time of Sd in hour
    REAL(kind = 8)::tSd

    ! surface moisture status [-]
    REAL(kind = 8)::xSM

    ! surface type
    INTEGER :: sfc_typ


    ! array of daytime series
    REAL(kind=8), DIMENSION(:,:),ALLOCATABLE :: dayArray

    ! daytime series
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: Sd  !< incoming solar radiation, W m-2
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: Ta  !< air temperature, degC
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: RH  !< relative humidity, %
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: pres!< Atmospheric pressure, mbar
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: tHr !< time, hour
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: Ts  ! surface temperature, degC
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: qa  ! specific humidity, kg kg-1
    REAL(kind = 8), DIMENSION(:),ALLOCATABLE :: qs  ! Saturated specific humidity at surface,kg kg-1

    ! conversion constant:
    REAL(kind = 8),PARAMETER:: C2K=273.15 ! degC to K

    ! thermodynamic values at standard temperature (293.15 K) and pressure (101325 pascals):
    REAL(kind = 8),PARAMETER:: cp_air = 1006.38    ! specific heat capacity of air, J kg-1 K-1
    REAL(kind = 8),PARAMETER:: Lv_air = 2264.705E3 ! latent heat of vaporization of water, J kg-1

    INTEGER :: lenDay,i,err,nVar,nPrm


    IF ( iflag == 0 ) THEN
       ! only print out the x vector

       WRITE ( *, '(a)' ) ''
       CALL r8vec_print(n,x,'x in fcnBo')

    ELSE  IF ( iflag == 1 ) THEN
       ! calculate fvec at x:
       ! pass unknown:
       xBo=x(1)

       ! pass forcing scales:
       ASd = prms(1)
       mSd = prms(2)
       ATa = prms(3)
       mTa = prms(4)
       tau = prms(5)
       mWS = prms(6)
       mWF = prms(7)
       mAH = prms(8)

       ! pass sfc. property scales:
       xalb  = prms(9)
       xemis = prms(10)
       xcp   = prms(11)
       xk    = prms(12)
       xch   = prms(13)
       xSM   = MIN(prms(14),1.0)

       ! pass tSd:
       tSd=prms(15)

       ! pass tSd:
       sfc_typ=INT(prms(16))

       ! set number of parameters according to the above code
       nPrm=16

       ! number of met variables to pass:
       nVar=5
       ! length of daytime series
       lenDay=(m-nPrm)/nVar
       ! allocate daytime series
       ALLOCATE(dayArray(nVar,lenDay), stat=err)
       IF ( err/= 0) PRINT *, "dayArray: Allocation request denied"
       dayArray=RESHAPE(prms(nPrm+1:SIZE(prms)), shape=(/nVar,lenDay/),order=(/2,1/))

       ! pass daytime series
       ! Sd:
       ALLOCATE(Sd(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "Sd: Allocation request denied"
       Sd(:)=dayArray(1,:)
       ! Ta:
       ALLOCATE(Ta(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "Ta: Allocation request denied"
       Ta(:)=dayArray(2,:)
       ! RH:
       ALLOCATE(RH(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "RH: Allocation request denied"
       RH(:)=dayArray(3,:)
       ! pres:
       ALLOCATE(pres(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "pres: Allocation request denied"
       pres(:)=dayArray(4,:)
       ! tHr:
       ALLOCATE(tHr(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "tHr: Allocation request denied"
       tHr(:)=dayArray(5,:)


       !  PRINT*, 'n',n
       !  PRINT*, 'lenDay',lenDay
       !  PRINT*, 'nVar',nVar
       !  PRINT*, 'shape of dayArray',SHAPE(dayArray)
       !  PRINT*, 'shape of met variables',SHAPE(Ta),SHAPE(RH),SHAPE(pres),SHAPE(tHr)

       ALLOCATE(Ts(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "Ts: Allocation request denied"
       ALLOCATE(qs(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "qs: Allocation request denied"
       ALLOCATE(qa(lenDay), stat=err)
       IF ( err/= 0) PRINT *, "qa: Allocation request denied"


       ! calculate sums of QH and QE series in the daytime
       IF ( xSM==0 ) THEN
          xBo=1000 ! extremely arid
       ELSE
          ! PRINT*, 'lenDay',lenDay
          DO i = 1, lenDay, 1
             ! calculate surface temperature
             CALL AnOHM_xTs(&
                  sfc_typ,&
                  ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
                  xalb,xemis,xcp,xk,xch,xBo,&   ! input: sfc properties
                  tSd,& !input: peaking time of Sd in hour
                  tHr(i),&! input: time in hour
                  Ts(i))! output: surface temperature, K

             ! convert K to degC
             Ts(i)=Ts(i)-C2K

             ! calculate saturation specific humidity
             qs(i)=qsat_fn(Ts(i),pres(i))

             ! calculate specific humidity
             qa(i)=qa_fn(Ta(i),RH(i),pres(i))

             ! PRINT*,''
             ! PRINT*, 'tHr',tHr(i)
             ! PRINT*, 'Sd',Sd(i)
             ! PRINT*, 'Ts',Ts(i)
             ! PRINT*, 'pres',pres(i)
             ! PRINT*, 'qs',qs(i)
             ! PRINT*, 'Ta',Ta(i)
             ! PRINT*, 'RH',RH(i)
             ! PRINT*, 'pres',pres(i)
             ! PRINT*, 'qa',qa(i)

          END DO

          ! below for testing:
          ! rho_air=1.293, air density, kg m-3
          ! ra=60, nominal aerodynamic resistence, s m-1
          ! PRINT*,''
          ! PRINT*, 'QH:',SUM(cp_air*(Ts-Ta)*1.293/60)/lenDay
          ! PRINT*, 'xSM:',xSM
          ! PRINT*, 'QE:',SUM(xSM*Lv_air*(qs-qa)*1.293/60)/lenDay
          ! PRINT*,''

          xBo = SUM(cp_air*(Ts-Ta))/& ! sum(QH)
               SUM(xSM*Lv_air*(qs-qa))! sum(QE)

       END IF

       ! f(Bo)-Bo==0
       fvec(1)=x(1)-xBo

       ! deallocate arrays
       IF (ALLOCATED(dayArray)) DEALLOCATE(dayArray, stat=err)
       IF ( err/= 0) PRINT *, "dayArray: Deallocation request denied"
       IF (ALLOCATED(Sd)) DEALLOCATE(Sd, stat=err)
       IF ( err/= 0) PRINT *, "Sd: Deallocation request denied"
       IF (ALLOCATED(Ta)) DEALLOCATE(Ta, stat=err)
       IF ( err/= 0) PRINT *, "Ta: Deallocation request denied"
       IF (ALLOCATED(RH)) DEALLOCATE(RH, stat=err)
       IF ( err/= 0) PRINT *, "RH: Deallocation request denied"
       IF (ALLOCATED(pres)) DEALLOCATE(pres, stat=err)
       IF ( err/= 0) PRINT *, "pres: Deallocation request denied"
       IF (ALLOCATED(tHr)) DEALLOCATE(tHr, stat=err)
       IF ( err/= 0) PRINT *, "tHr: Deallocation request denied"
       IF (ALLOCATED(Ts)) DEALLOCATE(Ts, stat=err)
       IF ( err/= 0) PRINT *, "Ts: Deallocation request denied"
       IF (ALLOCATED(qa)) DEALLOCATE(qa, stat=err)
       IF ( err/= 0) PRINT *, "qa: Deallocation request denied"
       IF (ALLOCATED(qs)) DEALLOCATE(qs, stat=err)
       IF ( err/= 0) PRINT *, "qs: Deallocation request denied"

    END IF
    RETURN

  END SUBROUTINE fcnBo
  !========================================================================================

  !========================================================================================
  !> calculate saturation vapor pressure (es)
  !> at air temperature (Ta)
  !> (MRR, 1987)
  FUNCTION esat_fn(Ta) RESULT(esat)
    REAL (KIND(1D0))::Ta  !< air temperature [degC]
    REAL (KIND(1D0))::esat  !< saturation vapor pressure [hPa]

    REAL (KIND(1D0)),PARAMETER::A=6.106 !< Teten coefficients
    REAL (KIND(1D0)),PARAMETER::B=17.27 !< Teten coefficients
    REAL (KIND(1D0)),PARAMETER::C=237.3 !< Teten coefficients

    esat = A*EXP(B*Ta/(C+Ta))
  END FUNCTION esat_fn
  !========================================================================================

  !========================================================================================
  !> calculate saturation specific humidity (qsat)
  !> at air temperature (Ta) and atmospheric pressure (pres)
  !> (MRR, 1987)
  FUNCTION qsat_fn(Ta,pres) RESULT(qsat)
    REAL (KIND(1D0))::Ta  !< air temperature [degC]
    REAL (KIND(1D0))::es  !< saturation vapor pressure [hPa]
    REAL (KIND(1D0))::qsat!< saturation specific humidity [kg kg-1]
    REAL (KIND(1D0))::pres!< atmospheric pressure [hPa]

    REAL (KIND(1D0)),PARAMETER::molar         = 0.028965  !< Dry air molar fraction [kg mol-1]
    REAL (KIND(1D0)),PARAMETER::molar_wat_vap = 0.0180153 !< Molar fraction of water vapor [kg mol-1]

    es = esat_fn(Ta)
    qsat = (molar_wat_vap/molar)*es/pres!(rmh2o/rmair)*ES/PMB
  END FUNCTION qsat_fn
  !========================================================================================

  !========================================================================================
  !> convert relative humidity (RH) to specific humidity (qa)
  !> at air temperature (Ta) and atmospheric pressure (pres)
  FUNCTION qa_fn(Ta,RH,pres) RESULT(qa)

    REAL (KIND(1D0))::Ta  !< air temperature [degC]
    REAL (KIND(1D0))::RH  !< relative humidity [%]
    REAL (KIND(1D0))::ea  !< vapor pressure [hPa]
    REAL (KIND(1D0))::es  !< saturation vapor pressure [hPa]
    REAL (KIND(1D0))::pres!< atmospheric pressure [hPa]
    REAL (KIND(1D0))::qa  !< saturation specific humidity [kg kg-1]

    REAL (KIND(1D0)),PARAMETER::molar         = 0.028965  !< Dry air molar fraction [kg mol-1]
    REAL (KIND(1D0)),PARAMETER::molar_wat_vap = 0.0180153 !< Molar fraction of water vapor [kg mol-1]

    es = esat_fn(Ta)
    ea = es*RH/100
    qa = (molar_wat_vap/molar)*ea/pres!(rmh2o/rmair)*ES/PMB
  END FUNCTION qa_fn
  !========================================================================================

END MODULE AnOHM_module
