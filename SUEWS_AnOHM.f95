MODULE AnOHM_module

  IMPLICIT NONE
CONTAINS

  !========================================================================================
  SUBROUTINE AnOHM(qn1,qn1_store,qn1_av_store,&
       MetForcingData_grid,&
       alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
       sfr,nsurf,nsh,AnthropHeatMethod,id,&
       a1,a2,a3,qs)! output
    ! author: Ting Sun
    !
    ! purpose:
    ! calculate heat storage based within AnOHM framework.
    !
    ! TODO modify the notes
    ! input:
    ! 1) Gridiv: grid number, a global variable.
    ! with Gridiv, required met forcings and sfc. properties can be loaded for calculation.
    !
    ! output:
    ! 1) grid ensemble heat storage.
    ! QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! https://doi.org/10.5194/gmd-2016-300
    !
    ! history:
    ! 20160301: initial version
    ! 20170109  updated dqndt calculation in accordance with SUEWS_OHM.f95 (HCW)
    !========================================================================================

    ! USE allocateArray
    ! USE data_in
    ! USE defaultNotUsed
    ! USE gis_data
    ! USE sues_data
    ! USE time
    ! USE InitialCond

    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:,:)::&
         MetForcingData_grid!,&! met forcing array of grid
    REAL(KIND(1d0)),INTENT(in)::&
         qn1,& ! net all-wave radiation
         sfr(nsurf)! surface
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::&
         alb, & ! albedo,
         emis,    & ! emissivity,
         cpAnOHM, & ! heat capacity,
         kkAnOHM, & ! thermal conductivity
         chAnOHM    ! bulk transfer coef.

    INTEGER,INTENT(in)::&
         id,&
         AnthropHeatMethod,&!AnthropHeat option
         nsurf,& ! number of surfaces
         nsh ! number of timesteps in one hour

    REAL(KIND(1d0)),INTENT(inout)::&
         qn1_store(nsh),qn1_av_store(2*nsh+1)

    REAL(KIND(1d0)),INTENT(out):: &
         a1,a2,a3,&!AnOHM coefficients of grid
         qs! storage heat flux

    INTEGER :: Gridiv,is,xid
    REAL(KIND(1d0)),PARAMETER::NotUsed=-55.5
    INTEGER,PARAMETER::notUsedI=-55
    LOGICAL :: idQ ! whether id contains enough data

    REAL(KIND(1d0))                  :: dqndt       ! rate of change of net radiation [W m-2 h-1] at t-2
    ! REAL(KIND(1d0))                  :: surfrac     ! surface fraction accounting for SnowFrac if appropriate
    REAL(KIND(1d0)),DIMENSION(nsurf) :: xa1,xa2,xa3 ! temporary AnOHM coefs.
    ! REAL(KIND(1d0))                  :: qn1_av      ! average net radiation over previous hour [W m-2]
    ! REAL(KIND(1d0))                  :: nsh_nna     ! number of timesteps per hour with non -999 values (used for spinup)

    ! INTERFACE
    !    SUBROUTINE AnOHM_coef(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,& !input
    !         xa1,xa2,xa3)                         ! output
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: sfc_typ,xid, xgrid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),INTENT(out):: xa1,xa2,xa3
    !    END SUBROUTINE AnOHM_coef
    ! END INTERFACE

    ! ------AnOHM coefficients --------
    !   print*, 'n surf:', nsurf
    !   Loop through surface types at the beginning of a day------------------------
    ! IF ( it==0 .AND. imin==5 ) THEN
    !       print*, '----- AnOHM called -----'
    !       print*, 'Grid@id:', Gridiv, id
    ! ------Set to zero initially------
    !  a1AnOHM(Gridiv) = 0   ![-]
    !  a2AnOHM(Gridiv) = 0   ![h]
    !  a3AnOHM(Gridiv) = 0   ![W m-2]
    !----------------------------------
    xa1 = 0.1
    xa2 = 0.2
    xa3 = 10
    ! PRINT*, 'MetForcingData_grid in AnOHM, shape:',SHAPE(MetForcingData_grid)
    ! PRINT*, MetForcingData_grid(1,1)
    ! PRINT*, 'qn1_store info in AnOHM, shape',SHAPE(qn1_store)
    ! PRINT*, 'qn1_av_store info in AnOHM, shape',SHAPE(qn1_av_store)

    ! to test if the current met block contains enough data for AnOHM
    ! daylight hours >= 8
    idQ=COUNT(MetForcingData_grid(:,2)==id .AND. & ! day of year
         MetForcingData_grid(:,4)==0 .AND. & ! minutes
         MetForcingData_grid(:,15)>0) & ! Sd
         .GE. 8

    ! PRINT*, idQ
    IF ( idQ ) THEN
       ! given enough data, calculate coefficients of day `id`
       xid=id
    ELSE
       ! otherwise calculate coefficients of yesterday: `id-1`
       xid=id-1
    END IF

    DO is=1,nsurf

       !   call AnOHM to calculate the coefs.
       CALL AnOHM_coef(is,xid,Gridiv,MetForcingData_grid,AnthropHeatMethod,& !input
            alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
            xa1(is),xa2(is),xa3(is))                         ! output
       ! print*, 'AnOHM_coef are: ',xa1,xa2,xa3

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
  SUBROUTINE AnOHM_coef(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,& !input
       alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
       xa1,xa2,xa3)                         ! output

    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
    ! generic subroutine to call AnOHM_coef_land and AnOHM_coef_water based on surface types.
    !
    ! input:
    ! 1) sfc_typ: surface type.
    ! these properties will be loaded:
    ! xemis: emissivity, 1
    ! xcp: heat capacity, J m-3
    ! xk: thermal conductivity, W m K-1
    ! xch: bulk turbulent transfer coefficient,
    ! Bo: Bowen ratio (i.e. QH/QE), 1
    ! 2) xid: day of year
    ! will be used to retrieve forcing diurnal cycles of ixd.
    !
    ! output:
    ! a1, a2, and a3
    ! in the relationship:
    ! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! https://doi.org/10.5194/gmd-2016-300
    !
    ! history:
    ! 20160222: initial version
    !========================================================================================


    ! USE allocateArray
    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in):: sfc_typ, xid, xgrid,AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:) ::&
         alb,  &
         emis, &
         cpAnOHM,  &
         kkAnOHM,  &
         chAnOHM
    REAL(KIND(1d0)),INTENT(in) :: MetForcingData_grid(:,:)

    !   output
    REAL(KIND(1d0)),INTENT(out) :: xa1, xa2, xa3


    ! locally saved variables:
    ! if coefficients have been calculated, just reload them
    ! otherwise, do the calculation
    INTEGER,SAVE :: id_save=-999,grid_save= -999
    REAL(KIND(1d0)), SAVE:: coeff_grid_day(7,3)=-999.

    INTEGER :: WaterSurf=7
    ! INTERFACE
    !    SUBROUTINE AnOHM_coef_land(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
    !         xa1,xa2,xa3)
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: sfc_typ,xid, xgrid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),INTENT(out):: xa1,xa2,xa3
    !    END SUBROUTINE AnOHM_coef_land
    !    SUBROUTINE AnOHM_coef_water(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
    !         xa1,xa2,xa3)
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: sfc_typ,xid, xgrid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),INTENT(out):: xa1,xa2,xa3
    !    END SUBROUTINE AnOHM_coef_water
    ! END INTERFACE
    ! PRINT*, 'xid,id_save',xid,id_save
    ! PRINT*, 'xgrid,grid_save',xgrid,grid_save
    ! PRINT*, 'sfc_typ',sfc_typ
    ! PRINT*, 'coeff_grid_day',coeff_grid_day(sfc_typ,:)
    IF ( xid==id_save .AND. xgrid ==grid_save) THEN
       ! if coefficients have been calculated, just reload them
       xa1=coeff_grid_day(sfc_typ,1)
       xa2=coeff_grid_day(sfc_typ,2)
       xa3=coeff_grid_day(sfc_typ,3)
    ELSE
       ! otherwise, do the calculation
       IF ( sfc_typ<WaterSurf ) THEN
          CALL AnOHM_coef_land(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
               alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
               xa1,xa2,xa3)            ! output
       ELSE
          CALL AnOHM_coef_water(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
               alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
               xa1,xa2,xa3)            ! output

          ! mark the day and grid as finished when the last surface (i.e., water) is done
          id_save =xid
          grid_save=xgrid
       END IF
       coeff_grid_day(sfc_typ,:)=(/xa1,xa2,xa3/)
    END IF
    ! PRINT*, 'coeff_grid_day',coeff_grid_day(sfc_typ,:)

  END SUBROUTINE AnOHM_coef
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_coef_land(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
       alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
       xa1,xa2,xa3)            ! output
    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
    !
    ! input:
    ! 1) sfc_typ: surface type.
    ! these properties will be loaded:
    ! xemis: emissivity, 1
    ! xcp: heat capacity, J m-3
    ! xk: thermal conductivity, W m K-1
    ! xch: bulk turbulent transfer coefficient,
    ! Bo: Bowen ratio (i.e. QH/QE), 1
    ! 2) xid: day of year
    ! will be used to retrieve forcing diurnal cycles of ixd.
    !
    ! output:
    ! a1, a2, and a3
    ! in the relationship:
    ! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! https://doi.org/10.5194/gmd-2016-300
    !
    ! history:
    ! 20160222: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE data_in
    ! USE defaultNotUsed
    ! USE gis_data
    ! USE sues_data
    ! USE time

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in):: sfc_typ, xid,xgrid, AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:) ::&
         alb,  &
         emis, &
         cpAnOHM,  &
         kkAnOHM,  &
         chAnOHM
    REAL(KIND(1d0)),INTENT(in):: MetForcingData_grid(:,:)

    !   output
    REAL(KIND(1d0)),INTENT(out) :: xa1, xa2, xa3

    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    ! REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K
    ! REAL(KIND(1d0)), PARAMETER :: CRA   = 915.483          ! converting RA (aerodyn. res.) to bulk trsf. coeff., [kg s-3]


    !   forcing scales
    REAL(KIND(1d0)):: &
         ASd,mSd,    & ! solar radiation
         ATa,mTa,    & ! air temperature
         tau,        & ! phase lag between Sd and Ta (Ta-Sd)
         mWS,mWF,mAH ! mean values of WS, WF and AH


    !   sfc. properties:
    REAL(KIND(1d0)) :: &
         xalb,   &    !  albedo,
         xemis,  &    !  emissivity,
         xcp,    &    !  heat capacity,
         xk,     &    !  thermal conductivity,
         xch,    &    !  bulk transfer coef.
         xBo          !  Bowen ratio
    ! INTERFACE
    !    SUBROUTINE AnOHM_Fc(xid,MetForcingData_grid,AnthropHeatMethod,ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: xid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),INTENT(out):: ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH
    !    END SUBROUTINE AnOHM_Fc
    ! END INTERFACE
    ! PRINT*, 'MetForcingData_grid before AnOHM_Fc, shape:',SHAPE(MetForcingData_grid)
    ! PRINT*, MetForcingData_grid(1,1)
    ! load forcing characteristics:
    CALL AnOHM_Fc(xid,MetForcingData_grid,AnthropHeatMethod,& !input
         ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)! output

    !   load sfc. properties:
    CALL AnOHM_SfcLoad(sfc_typ,xid,xgrid,&            ! input
         alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
         xalb,xemis,xcp,xk,xch,xBo) ! output

    ! calculate coefficients:
    CALL AnOHM_coef_land_cal(ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&! input: forcing
         xalb,xemis,xcp,xk,xch,xBo,&   ! input: sfc properties
         xa1,xa2,xa3) ! output

  END SUBROUTINE AnOHM_coef_land
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_coef_land_cal(&
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
       xalb,xemis,xcp,xk,xch,xBo,&   ! input: sfc properties
       xa1,xa2,xa3)            ! output

    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
    !
    ! input:
    ! 1) sfc_typ: surface type.
    ! these properties will be loaded:
    ! xemis: emissivity, 1
    ! xcp: heat capacity, J m-3
    ! xk: thermal conductivity, W m K-1
    ! xch: bulk turbulent transfer coefficient,
    ! Bo: Bowen ratio (i.e. QH/QE), 1
    ! 2) xid: day of year
    ! will be used to retrieve forcing diurnal cycles of ixd.
    !
    ! output:
    ! a1, a2, and a3
    ! in the relationship:
    ! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! https://doi.org/10.5194/gmd-2016-300
    !
    ! history:
    ! 20160222: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE data_in
    ! USE defaultNotUsed
    ! USE gis_data
    ! USE sues_data
    ! USE time

    IMPLICIT NONE

    !   input
    REAL(KIND(1d0)),INTENT(in):: &
                                ! input: forcing scales
         ASd,mSd,    & ! solar radiation
         ATa,mTa,    & ! air temperature
         tau,        & ! phase lag between Sd and Ta (Ta-Sd)
         mWS,mWF,mAH,& ! mean values of WS, WF and AH
                                ! input: sfc properties
         xalb,& ! albedo,
         xemis,&! emissivity,
         xcp,&  ! heat capacity,
         xk,&   ! thermal conductivity,
         xch,&  ! bulk transfer coef.
         xBo  ! Bowen ratio


    !   output
    REAL(KIND(1d0)),INTENT(out) :: xa1, xa2, xa3

    !   constant
    REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K
    ! REAL(KIND(1d0)), PARAMETER :: CRA   = 915.483          ! converting RA (aerodyn. res.) to bulk trsf. coeff., [kg s-3]

    !   local variables:
    REAL(KIND(1d0))    :: beta               ! inverse Bowen ratio
    REAL(KIND(1d0))    :: f,fL,fT            ! energy redistribution factors
    REAL(KIND(1d0))    :: lambda             ! thermal diffusivity
    REAL(KIND(1d0))    :: delta,m,n          ! water flux related variables
    REAL(KIND(1d0))    :: ms,ns              ! m, n related
    REAL(KIND(1d0))    :: gamma              ! phase lag scale
    REAL(KIND(1d0))    :: ATs,mTs            ! amplitude and mean value of surface temperature
    REAL(KIND(1d0))    :: ceta,cphi          ! phase related temporary variables
    REAL(KIND(1d0))    :: eta,phi,xlag       ! phase related temporary variables
    REAL(KIND(1d0))    :: xx1,xx2,xx3,xchWS  ! temporary use
    ! REAL(KIND(1d0))    :: solTs              ! surface temperature
    LOGICAL :: flagGood = .TRUE.  ! quality flag, T for good, F for bad
    ! INTEGER :: lenDay=24,err


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
    flagGood = .TRUE.
    ! PRINT*, flagGood




    !   calculate sfc properties related parameters:
    xchWS = xch*mWS
    xchWS = 0.5*xchWS

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
    !     xa1 = min(xa1,0.7)

    !   a2:
    xa2 = (ceta*SIN(xlag))/(OMEGA*cphi)
    xa2 = xa2/3600 ! convert the unit from s-1 to h-1
    !     xa2 = 0.5*xa2  ! modify it to be more reasonable
    !     xa2 = max(min(xa2,0.7),-0.7)

    !   a3:
    xa3  = -xa1*(fT/f)*(mSd*(1-xalb))-mAH

    !   quality checking:
    !   quality checking of forcing conditions
    IF ( ASd < 0 .OR. ATa < 0 .OR. ATs < 0 .OR. tau<-4.0/12*Pi) flagGood = .FALSE.
    !   quality checking of a1
    IF ( .NOT. (xa1>0 .AND. xa1<0.7)) THEN
       flagGood = .FALSE.
       IF (xa1 >0.7) xa1=MAX(0.7,xa1)
    ENDIF
    !   quality checking of a2
    IF ( .NOT. (xa2>-0.5 .AND. xa2<0.5)) THEN
       flagGood = .FALSE.
       !  IF ( xa2>0.5) xa2 = 0.5
       !  IF (xa2<-0.5) xa2 = -0.5
    ENDIF
    !   quality checking of a3
    IF ( .NOT. (xa3<0)) flagGood = .FALSE.

    !     print*, 'ceta,cphi', ceta,cphi
    !     print*, 'tau,eta,phi,xlag in deg:',tau/pi*180,eta/pi*180,phi/pi*180,xlag/pi*180

    ! PRINT*, '********sfc_typ: ',sfc_typ,' end********'

  END SUBROUTINE AnOHM_coef_land_cal
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_coef_water(sfc_typ,xid,xgrid,MetForcingData_grid,AnthropHeatMethod,&   ! input
       alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
       xa1,xa2,xa3)            ! output
    ! author: Ting Sun
    ! date: 20161124
    !
    ! purpose:
    ! designed for water surface
    ! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
    ! NB: this SUBROUTINE hasn't been well tested.
    !
    ! input:
    ! 1) sfc_typ: surface type.
    ! these properties will be loaded:
    ! xemis: emissivity, 1
    ! xcp: heat capacity, J m-3
    ! xk: thermal conductivity, W m K-1
    ! xch: bulk turbulent transfer coefficient,
    ! Bo: Bowen ratio (i.e. QH/QE), 1
    ! xeta: effective absorption coefficient, 1
    ! xmu: effective absorption fraction, m-1
    ! 2) xid: day of year
    ! will be used to retrieve forcing diurnal cycles of ixd.
    !
    ! output:
    ! a1, a2, and a3
    ! in the relationship:
    ! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! the AnOHM-water paper to be added.
    !
    ! history:
    ! 20161124: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE data_in
    ! USE defaultNotUsed
    ! USE gis_data
    ! USE sues_data
    ! USE time

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in):: sfc_typ, xid, xgrid,AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in)::&
         MetForcingData_grid(:,:)
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:)::&
         alb, & ! albedo,
         emis,    & ! emissivity,
         cpAnOHM, & ! heat capacity,
         kkAnOHM, & ! thermal conductivity
         chAnOHM    ! bulk transfer coef.

    !   output
    REAL(KIND(1d0)),INTENT(out) :: xa1, xa2, xa3

    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    ! REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K


    !   forcing scales
    REAL(KIND(1d0)):: &
         ASd,mSd,    & ! solar radiation
         ATa,mTa,    & ! air temperature
         tau,        & ! phase lag between Sd and Ta (Ta-Sd)
         mWS,mWF,mAH ! mean values of WS, WF and AH

    !   sfc. properties:
    REAL(KIND(1d0)) :: &
         xalb, & ! albedo,
         xemis,& ! emissivity,
         xcp,   & ! heat capacity,
         xk,     & ! thermal conductivity,
         xch,   & ! bulk transfer coef.
         xBo,          & ! Bowen ratio
         xeta,         & ! effective absorption coefficient
         xmu ! effective absorption fraction
    ! INTERFACE
    !    SUBROUTINE AnOHM_Fc(xid,MetForcingData_grid,AnthropHeatMethod,ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: xid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),INTENT(out):: ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH
    !    END SUBROUTINE AnOHM_Fc
    ! END INTERFACE

    ! load forcing characteristics:
    CALL AnOHM_Fc(xid,MetForcingData_grid,AnthropHeatMethod,& !input
         ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)! output
    !   write(*,*) ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH

    !   load sfc. properties:
    CALL AnOHM_SfcLoad(sfc_typ,xid,xgrid,&            ! input
         alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
         xalb,xemis,xcp,xk,xch,xBo) ! output
    !   write(*,*) 'here the properties:'
    !   write(*,*) xalb,xemis,xcp,xk,xch,xBo
    ! !   give fixed values for the moment
    xeta  = 0.3
    xmu   = 0.2

    ! calculate coefficients:
    CALL AnOHM_coef_water_cal(ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&
         xalb,xemis,xcp,xk,xch,xBo,xeta,xmu,&
         xa1,xa2,xa3)

  END SUBROUTINE AnOHM_coef_water
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_coef_water_cal(&
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH,&   ! input: forcing
       xalb,xemis,xcp,xk,xch,xBo,xeta,xmu,&   ! input: sfc properties
       xa1,xa2,xa3)            ! output
    ! author: Ting Sun
    ! date: 20161124
    !
    ! purpose:
    ! designed for water surface
    ! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
    ! NB: this SUBROUTINE hasn't been well tested.
    !
    ! input: see the code below
    !
    ! output:
    ! a1, a2, and a3
    ! in the relationship:
    ! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
    !
    ! ref:
    ! the AnOHM-water paper to be added.
    !
    ! history:
    ! 20161124: initial version
    !========================================================================================
    USE allocateArray
    USE data_in
    USE defaultNotUsed
    USE gis_data
    USE sues_data
    USE time

    IMPLICIT NONE

    !   input
    REAL(KIND(1d0)),INTENT(in):: &
                                ! input: forcing scales
         ASd,mSd,    & ! solar radiation
         ATa,mTa,    & ! air temperature
         tau,        & ! phase lag between Sd and Ta (Ta-Sd)
         mWS,mWF,mAH,& ! mean values of WS, WF and AH
                                ! input: sfc properties
         xalb,& ! albedo,
         xemis,&! emissivity,
         xcp,&  ! heat capacity,
         xk,&   ! thermal conductivity,
         xch,&  ! bulk transfer coef.
         xBo,&  ! Bowen ratio
         xeta,& ! effective absorption coefficient
         xmu    ! effective absorption fraction

    !   output
    REAL(KIND(1d0)),INTENT(out):: xa1, xa2, xa3

    !   constant
    REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K


    !   local variables:
    REAL(KIND(1d0)) :: beta                   ! inverse Bowen ratio
    REAL(KIND(1d0)) :: f,fL,fT                ! energy redistribution factors
    REAL(KIND(1d0)) :: lambda,calb            ! temporary use
    REAL(KIND(1d0)) :: delta             ! sfc. temperature related variables
    ! REAL(KIND(1d0)) :: delta,m,n              ! sfc. temperature related variables
    REAL(KIND(1d0)) :: xm,xn                  ! m, n related
    ! REAL(KIND(1d0)) :: gamma              ! phase lag scale
    REAL(KIND(1d0)) :: phi              ! phase lag scale
    REAL(KIND(1d0)) :: ATs,mTs                ! surface temperature amplitude
    REAL(KIND(1d0)) :: czeta,ctheta           ! phase related temporary variables
    REAL(KIND(1d0)) :: zeta,theta,xlag           ! phase related temporary variables
    REAL(KIND(1d0)) :: xx1,xx2,xx3            ! temporary use
    REAL(KIND(1d0)) :: kappa                  ! temporary use
    REAL(KIND(1d0)) :: dtau,dpsi,dphi         ! temporary use
    REAL(KIND(1d0)) :: cdtau,cdpsi,cdphi      ! temporary use
    REAL(KIND(1d0)) :: xxT,xxkappa,xxdltphi,xchWS   ! temporary use
    LOGICAL :: flagGood = .TRUE.  ! quality flag, T for good, F for bad

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
    xx1   = (xk*xeta*xmu*calb*ASd*cdpsi)**2
    xx2=2*lambda*SQRT(xx1)*(calb*ASd*SIN(phi-dpsi)+f*ATa*SIN(tau+phi-dpsi))
    xx3=lambda**2*((calb*ASd+COS(tau)*f*ATa)**2+(SIN(tau)*f*ATa)**2)
    ATs=1/(cdtau*lambda)*SQRT(xx1+xx2+xx3)
    ! phase lag:
    xx1=(xk*kappa*calb*ASd+cdtau*f*ATa*SIN(tau+dtau))*lambda&
         +(xk*xeta*xmu*calb*ASd*cdphi*SIN(phi+dphi))
    xx2=((f+xk*kappa)*calb*ASd-cdtau*f*ATa*COS(tau+dtau))*lambda&
         -xk*xeta*xmu*calb*ASd*cdphi*COS(phi+dphi)
    delta=ATAN(xx1/xx2)

    !   calculate net radiation related parameters:
    ! phase lag:
    xx1=fL*(ATs*SIN(delta)-ATa*SIN(tau))
    xx2=calb*ASd-fL*(ATs*COS(delta)-ATa*COS(tau))
    theta=ATAN(xx1/xx2)
    ! amplitude:
    ctheta=SQRT(xx1**2+xx2**2)

    !   calculate heat storage related parameters:
    ! scales:
    xxT=SQRT(2.)*kappa*lambda*ATs
    xxkappa=cdpsi*xeta*xmu*ASd
    xxdltphi=COS(delta)*SIN(dpsi)*COS(phi)-SIN(delta)*COS(dpsi)*SIN(phi)
    ! phase lag:
    xx1=xxT*SIN(PI/4-delta)+xxkappa*SIN(phi+dpsi)
    xx2=xxT*SIN(PI/4+delta)-xxkappa*SIN(PHI-dpsi)
    zeta=ATAN(xx1/xx2)
    ! amplitude:
    xx1=2*SQRT(2.)*xxkappa*xxT*xxdltphi
    xx2=(1-COS(2*dpsi)*COS(2*phi))*xxkappa**2
    xx3=xxT**2
    czeta=xk/lambda*SQRT(xx1+xx2+xx3)


    !   calculate the OHM coeffs.:
    xlag=zeta-theta
    !   a1:
    xa1  = (czeta*COS(xlag))/ctheta

    !   write(*,*) 'ceta,xlag,cphi:', ceta,xlag,cphi
    !   a2:
    xa2  = (czeta*SIN(xlag))/(OMEGA*ctheta)
    xa2  = xa2/3600 ! convert the unit from s-1 to h-1

    !   a3:
    xa3  = mSd*(xalb-1)*(xeta+(fT-fL*xeta)/f*xa1)


    !   quality checking:
    !   quality checking of forcing conditions
    IF ( ASd < 0 .OR. ATa < 0 .OR. ATs < 0 .OR. tau<-4.0/12*Pi) flagGood = .FALSE.
    !   quality checking of a1
    IF ( .NOT. (xa1>0 .AND. xa1<0.7)) THEN
       flagGood = .FALSE.
       IF (xa1 >0.7) xa1=MAX(0.7,xa1)
    ENDIF
    !   quality checking of a2
    IF ( .NOT. (xa2>-0.5 .AND. xa2<0.5)) THEN
       flagGood = .FALSE.
       !  IF ( xa2>0.5) xa2 = 0.5
       !  IF (xa2<-0.5) xa2 = -0.5
    ENDIF
    !   quality checking of a3
    IF ( .NOT. (xa3<0)) flagGood = .FALSE.

    !   skip the first day for quality checking
    ! IF ( xid == 1 ) flagGood = .TRUE.

    ! print*,  'sfc_typ_water:', sfc_typ
    ! print*, 'a1,a2,a3:', xa1,xa2,xa3

  END SUBROUTINE AnOHM_coef_water_cal
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_Fc(xid,MetForcingData_grid,AnthropHeatMethod,& ! input
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)    ! output
    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the key parameters of a sinusoidal curve for AnOHM forcings
    ! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
    !
    ! input:
    ! 1) xid : day of year
    ! 2) xgrid: grid
    !
    ! output:
    ! 1) ASd, mSd: amplitude and mean value of Sd
    ! 2) ATa, mTa: amplitude and mean value of Ta
    ! 3) tau: phase lag between Sd and Ta
    ! 4) mWS: mean value of WS
    ! 5) mWF: mean value of WF
    ! 6) mAH: mean value of AH
    !
    ! history:
    ! 20160226: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE ColNamesInputFiles
    ! USE data_in
    ! USE defaultNotUsed
    ! USE initial
    ! USE sues_data
    ! USE time

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in) :: xid,AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)

    !   output
    REAL(KIND(1d0)),INTENT(out) :: ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH

    !   forcings:
    REAL(KIND(1d0)), DIMENSION(:),ALLOCATABLE:: &
         Sd,& !   incoming solar radiation
         Ta,& !   air temperature
         WS,& !   wind speed
         WF,& !   water flux density
         AH,& !   anthropogenic heat
         tHr  !   time in hour

    ! PRINT*, 'MetForcingData_grid in AnOHM_Fc, shape:',SHAPE(MetForcingData_grid)
    ! PRINT*, MetForcingData_grid(1,1),MetForcingData_grid(1,5)
    ! INTERFACE
    !    SUBROUTINE AnOHM_FcLoad(xid,MetForcingData_grid,AnthropHeatMethod,Sd,Ta,WS,WF,AH,tHr)
    !      IMPLICIT NONE
    !      INTEGER,INTENT(in) :: xid,AnthropHeatMethod
    !      REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)
    !      REAL(KIND(1d0)),DIMENSION(:),INTENT(out),ALLOCATABLE :: Sd,Ta,WS,WF,AH,tHr
    !    END SUBROUTINE AnOHM_FcLoad
    !    SUBROUTINE AnOHM_FcCal(Sd,Ta,WS,WF,AH,tHr,& ! input
    !         ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)       ! output
    !      REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: Sd,Ta,WS,WF,AH,tHr
    !      REAL(KIND(1d0)),INTENT(out) ::&
    !           ASd,mSd,&   ! Sd scales
    !           ATa,mTa,&   ! Ta scales
    !           tau,&       ! phase lag between Sd and Ta
    !           mWS,&       ! mean WS
    !           mWF,&       ! mean WF
    !           mAH         ! mean AH
    !    END SUBROUTINE AnOHM_FcCal
    ! END INTERFACE

    ! load forcing variables:
    CALL AnOHM_FcLoad(xid,MetForcingData_grid,AnthropHeatMethod,Sd,Ta,WS,WF,AH,tHr)
    ! calculate forcing scales for AnOHM:
    CALL AnOHM_FcCal(Sd,Ta,WS,WF,AH,tHr,ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)

  END SUBROUTINE AnOHM_Fc
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_FcLoad(xid,MetForcingData_grid,AnthropHeatMethod,& ! input
       Sd,Ta,WS,WF,AH,tHr) ! output
    ! author: Ting Sun
    !
    ! purpose:
    ! load forcing series for AnOHM_FcCal
    !
    ! input:
    ! 1) xid : day of year
    ! 2) xgrid: grid
    !
    ! output:
    ! hourly values between 00:00 and 23:00 (local time, inclusive):
    ! 1) Sd: incoming shortwave radiation,  W m-2
    ! 2) Ta: air temperature, K
    ! 3) WS: wind speed, m s-1
    ! 4) WF: water flux density, ???
    ! 5) AH: anthropogenic heat, W m-2
    ! 6) tHr: time in hour, hr
    !
    ! history:
    ! 20170720: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE ColNamesInputFiles
    ! USE data_in
    ! USE defaultNotUsed
    ! USE initial
    ! USE sues_data
    ! USE time

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in) :: xid,AnthropHeatMethod
    REAL(KIND(1d0)),INTENT(in) ::MetForcingData_grid(:,:)

    !   output
    REAL(KIND(1d0)),DIMENSION(:),INTENT(out), ALLOCATABLE:: Sd,Ta,WS,WF,AH,tHr

    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    ! REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K

    !   local variables:
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE :: subMet !

    ! INTEGER :: i,xx(8),err            ! temporary use
    INTEGER :: err            ! temporary use
    INTEGER :: lenMetData

    ! INTEGER :: irRange(2)   ! row range in MetData containing id values


    ! INTEGER, ALLOCATABLE :: lSub(:) ! array to retrieve data of the specific day (id)
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
    ALLOCATE(metMask(lenMetData))
    metMask=(MetForcingData_grid(:,2)==xid & ! day=xid
         .AND. MetForcingData_grid(:,4)==0)! tmin=0

    ! construct array for time and met variables
    ALLOCATE(subMet(lenMetData,6))
    subMet=RESHAPE(PACK(MetForcingData_grid(:,(/3,& !time: hour
         15,12,10,12,9& ! met: Sd, Ta, WS, WF, AH
         /)),&
         SPREAD(metMask, dim=2, ncopies=6)),& ! replicate mask vector to 2D array
         (/lenMetData,6/)) ! convert to target shape

    ! re-allocate arrays as their sizes may change during passing
    IF (ALLOCATED(tHr)) DEALLOCATE(tHr, stat=err)
    ALLOCATE(tHr(lenMetData))
    IF (ALLOCATED(Sd)) DEALLOCATE(Sd, stat=err)
    ALLOCATE(Sd(lenMetData))
    IF (ALLOCATED(Ta)) DEALLOCATE(Ta, stat=err)
    ALLOCATE(Ta(lenMetData))
    IF (ALLOCATED(WS)) DEALLOCATE(WS, stat=err)
    ALLOCATE(WS(lenMetData))
    IF (ALLOCATED(WF)) DEALLOCATE(WF, stat=err)
    ALLOCATE(WF(lenMetData))
    IF (ALLOCATED(AH)) DEALLOCATE(AH, stat=err)
    ALLOCATE(AH(lenMetData))

    ! load the sublist into forcing variables
    tHr = subMet(:,1)! time in hour:
    Sd  = subMet(:,2)
    Ta  = subMet(:,3)
    WS  = subMet(:,4)
    WF  = 0 ! set as 0 for the moment
    IF ( AnthropHeatMethod == 0 ) THEN
       AH = subMet(:,6)    ! read in from MetForcingData_grid,
    ELSE
       AH = 0 ! temporarily change to zero; TODO: chenge back to modelled value
       !  AH = mAH_grids(xid-1,xgrid)
    END IF

    ! IF ( xid==73 ) THEN
    !    PRINT*, 'in  AnOHM_FcLoad'
    !    PRINT*, 'id:',xid
    !    PRINT*, 'Sd:', Sd
    !    PRINT*, 'Ta:', Ta
    !    PRINT*, 'WS:', WS
    !    PRINT*, 'WF:', WF
    !    PRINT*, 'AH:', AH
    !    PRINT*, 'tHr:', tHr
    ! END IF


  END SUBROUTINE AnOHM_FcLoad
  !========================================================================================

  !========================================================================================
  SUBROUTINE AnOHM_FcCal(Sd,Ta,WS,WF,AH,tHr,& ! input
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)    ! output
    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the key parameters of a sinusoidal curve for AnOHM forcings
    ! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
    !
    ! input:
    ! hourly values between 00:00 and 23:00 (local time, inclusive):
    ! 1) Sd: incoming shortwave radiation,  W m-2
    ! 2) Ta: air temperature, K
    ! 3) WS: wind speed, m s-1
    ! 4) WF: water flux density, ???
    ! 5) AH: anthropogenic heat, W m-2
    ! 6) tHr: time in hour, hr
    !
    ! output:
    ! 1) ASd, mSd: amplitude and mean value of Sd
    ! 2) ATa, mTa: amplitude and mean value of Ta
    ! 3) tau: phase lag between Sd and Ta
    ! 4) mWS: mean value of WS
    ! 5) mWF: mean value of WF
    ! 6) mAH: mean value of AH
    !
    ! ref:
    !
    ! history:
    ! 20160224: initial version
    !========================================================================================
    IMPLICIT NONE

    !   input
    REAL(KIND(1d0)),DIMENSION(:),INTENT(in):: Sd,Ta,WS,WF,AH,tHr

    !   output
    REAL(KIND(1d0)),INTENT(out) ::&
         ASd,mSd,&   ! Sd scales
         ATa,mTa,&   ! Ta scales
         tau,&       ! phase lag between Sd and Ta
         mWS,&       ! mean WS
         mWF,&       ! mean WF
         mAH         ! mean AH

    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K

    !   local variables:
    REAL(KIND(1d0)) :: tSd,tTa                ! peaking timestamps
    REAL(KIND(1d0)),ALLOCATABLE :: &
         tHrDay(:),& ! daytime tHr when Sd>0
         selX(:)     ! daytime sublist of met variable when Sd>0

    ! REAL(KIND(1d0)) :: xx              ! temporary use
    INTEGER :: err,lenDay
    LOGICAL,DIMENSION(:),ALLOCATABLE::SdMask

    ! INTERFACE
    !    SUBROUTINE AnOHM_ShapeFit(tHr,obs,amp,mean,tpeak)
    !      IMPLICIT NONE
    !      !   input
    !      REAL(KIND(1d0)),INTENT(in) :: tHr(:)  ! time in hour
    !      REAL(KIND(1d0)),INTENT(in) :: obs(:)  ! observation
    !      !   output
    !      REAL(KIND(1d0)),INTENT(out) :: amp     ! amplitude
    !      REAL(KIND(1d0)),INTENT(out) :: mean    ! average
    !      REAL(KIND(1d0)),INTENT(out) :: tpeak   ! peaking time (h)
    !
    !    END SUBROUTINE AnOHM_ShapeFit
    !
    ! END INTERFACE


    lenDay=COUNT(Sd>5, dim=1)
    ALLOCATE(SdMask(lenDay), stat=err)
    IF ( err/= 0) PRINT *, "SdMask: Allocation request denied"
    SdMask=Sd>0.1

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
    !   modify ill-shaped days to go through
    IF ( ASd < 0 .OR. tSd > 15) THEN
       !         ASd = abs(ASd)
       !         tSd = 12 ! assume Sd peaks at 12:00LST
       CALL r8vec_print(lenDay,tHrDay,'tHrDay:')
       CALL r8vec_print(lenDay,selX,'Sd Day:')
       PRINT*, 'ASd:', ASd
       PRINT*, 'mSd:', mSd
       PRINT*, 'tSd:', tSd
    END IF
    ! PRINT*, 'Sd Day:', selX

    !   calculate sinusoidal scales of Ta:
    ! PRINT*, 'Calc. Ta...'
    selX=PACK(Ta, mask=SdMask)
    ATa=(MAXVAL(selX)-MINVAL(selX))/2
    mTa=SUM(selX)/lenDay
    tTa=12
    CALL AnOHM_ShapeFit(tHrDay,selX,ATa,mTa,tTa)
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
  SUBROUTINE AnOHM_ShapeFit(tHr,obs,amp,mean,tpeak)
    ! author: Ting Sun
    !
    ! purpose:
    ! calculate the key parameters of a sinusoidal curve for AnOHM forcings
    ! i.e., a, b, c in a*Sin(Pi/12*t+b)+c, where t is in hour
    !
    ! input:
    ! obs: hourly values between 10:00 and 16:00 (local time, inclusive)
    !
    ! output:
    ! 1) amp  : amplitude
    ! 2) mean : mean value
    ! 3) tpeak: daily peaking hour (h)
    !
    ! ref:
    !
    ! history:
    ! 20160224: initial version
    ! 20170719: use minpack to approximate the sinusoidal shape
    !========================================================================================
    IMPLICIT NONE

    !   input
    REAL(KIND(1d0)),INTENT(in) :: tHr(:)  ! time in hour
    REAL(KIND(1d0)),INTENT(in) :: obs(:)  ! observation

    !   output
    REAL(KIND(1d0)),INTENT(out) :: amp     ! amplitude
    REAL(KIND(1d0)),INTENT(out) :: mean    ! average
    REAL(KIND(1d0)),INTENT(out) :: tpeak   ! peaking time (h)

    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    ! REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K

    !   local variables:
    ! REAL(KIND(1d0))    :: coefs(3,7)           ! coefficients for least squre solution
    ! REAL(KIND(1d0))    :: aCosb,aSinb,a,b,c    ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
    ! REAL(KIND(1d0))    :: xx,mbias,mbiasObs    ! temporary use
    INTEGER :: m,n,info,err         ! temporary use

    ! INTEGER ( kind = 4 ), PARAMETER :: m = 11
    ! INTEGER ( kind = 4 ), PARAMETER :: n = 3
    ! INTEGER ( kind = 4 ) :: ldfjac = m

    ! EXTERNAL fSin
    REAL ( KIND(1d0) ),ALLOCATABLE:: fvec(:),x(:)

    ! INTEGER ( kind = 4 ) iflag
    ! INTEGER ( kind = 4 ) info
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
  SUBROUTINE fSin(m, n, x, xdat, ydat, fvec, iflag)
    ! purpose:
    ! sinusoidal function f(t) for fitting:
    ! f(t) = mean+amp*Sin(Pi/12(t-delta))
    ! x    = (/mean,amp,delta/) contains the fitting parameters


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
  SUBROUTINE AnOHM_SfcLoad(sfc_typ,xid,xgrid,&            ! input
       alb, emis, cpAnOHM, kkAnOHM, chAnOHM,&! input
       xalb,xemis,xcp,xk,xch,xBo) ! output

    ! author: Ting Sun
    !
    ! purpose:
    ! load surface properties.
    !
    ! input:
    ! 1) sfc : surface type
    ! 2) xid : day of year
    ! 3) grid: grid number
    !
    ! output:
    ! surface properties:
    ! 1) alb,    ! albedo,
    ! 2) emiss,  ! emissivity,
    ! 3) cp,     ! heat capacity,
    ! 4) k,      ! thermal conductivity,
    ! 5) ch,     ! bulk transfer coef.
    ! 6) Bo      ! Bowen ratio
    !
    !
    ! history:
    ! 20160228: initial version
    !========================================================================================
    ! USE allocateArray
    ! USE ColNamesInputFiles
    ! USE data_in
    ! USE defaultNotUsed
    ! USE initial
    ! USE sues_data
    ! USE time
    ! USE InitialCond

    IMPLICIT NONE

    !   input
    INTEGER,INTENT(in) :: sfc_typ, xid, xgrid
    REAL(KIND(1d0)),INTENT(in),DIMENSION(:) ::&
         alb,  &
         emis, &
         cpAnOHM,  &
         kkAnOHM,  &
         chAnOHM

    !   output:
    !   sfc_typ. properties:
    REAL(KIND(1d0)),INTENT(out) ::&
         xalb,       &    !  albedo,
         xemis,      &    !  emissivity,
         xcp,        &    !  heat capacity,
         xk,         &    !  thermal conductivity,
         xch,        &    !  bulk transfer coef.
         xBo              !  Bowen ratio
    INTEGER :: xxx
    !   constant
    ! REAL(KIND(1d0)), PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
    ! REAL(KIND(1d0)), PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
    ! REAL(KIND(1d0)), PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
    ! REAL(KIND(1d0)), PARAMETER :: C2K   = 273.15           ! degC to K


    !   load properties from global variables
    xalb  = alb(sfc_typ)
    xemis = emis(sfc_typ)
    xcp   = cpAnOHM(sfc_typ)
    xk    = kkAnOHM(sfc_typ)
    xch   = chAnOHM(sfc_typ)

    !   print*, 'xalb:', xalb
    !   print*, 'xemis:', xemis
    !   print*, 'sfc_typ:', sfc_typ
    !   print*, 'xcp:', xcp
    !   print*, 'xk:', xk
    !   print*, 'xch:', xch

    ! suppress warnings
    xxx=xid
    xxx=xgrid


    !   load Bowen ratio of yesterday from DailyState calculation
    ! IF ( BoInit==NAN ) THEN
    !    CALL ErrorHint(68,'No initial Bowen ratio found in InitialConditions.',BoInit,notUsed,notUsedI)
    ! ELSE
    !    xBo = BoInit ! get Initial value from initial condition namelist
    ! END IF
    ! set as previous day's Bo, 20160708 TS
    ! xBo = Bo_grids(xid-1,xgrid) ! default Bo will be read in when xid = 0
    xBo=0
    ! IF ( xid==1 ) THEN
    !
    !    !  xBo = Bo_grids(xid-1,xgrid)
    ! ELSE
    !    !  IF ( BoAnOHMEnd(xgrid) > BoAnOHMStart(xgrid) ) THEN
    !    !     xBo = BoAnOHMEnd(xgrid)/2+BoAnOHMStart(xgrid)/2
    !    !  ELSE
    !    !     xBo = (BoAnOHMStart(xgrid)+BoAnOHMEnd(xgrid))/2
    !    !  END IF
    !
    ! END IF
    !     if ( xBo>30 .or. xBo<-1 ) write(unit=*, fmt=*) "Bo",xBo

    !     xBo = 2.
    xBo = MAX(MIN(xBo,30.),0.1)

    ! if ( xid==1 ) then
    !   print*, 'xBo:', xBo
    ! end if
    !   print*, 'xBo:', xBo

  END SUBROUTINE AnOHM_SfcLoad
  !========================================================================================



END MODULE AnOHM_module
