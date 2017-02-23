!========================================================================================
SUBROUTINE AnOHM_v2016(Gridiv)
  ! author: Ting Sun
  !
  ! purpose:
  ! calculate heat storage based within AnOHM framework.
  !
  ! input:
  ! 1) Gridiv: grid number, a global variable.
  ! with Gridiv, required met forcings and sfc. properties can be loaded for calculation.
  !
  ! output:
  ! 1) grid ensemble heat storage.
  ! QS = a1*(Q*)+a2*(dQ*/dt)+a3
  !
  ! ref:
  ! the AnOHM paper to be added.
  !
  ! history:
  ! 20160301: initial version
  ! 20170109  updated dqndt calculation in accordance with SUEWS_OHM.f95 (HCW) 
  !========================================================================================

  USE allocateArray
  USE data_in
  USE defaultNotUsed
  USE gis_data
  USE sues_data
  USE time
  USE InitialCond

  IMPLICIT NONE

  !INTEGER :: i,ii
  INTEGER :: Gridiv

  REAL    :: dqndt        ! rate of change of net radiation [W m-2 h-1] at t-2
  REAL    :: surfrac      ! surface fraction accounting for SnowFrac if appropriate
  REAL    :: xa1,xa2,xa3  ! temporary AnOHM coefs.
  REAL(KIND(1d0))    :: qn1_av       ! average net radiation over previous hour [W m-2]
  REAL(KIND(1d0))    :: nsh_nna      ! number of timesteps per hour with non -999 values (used for spinup)

  ! ------AnOHM coefficients --------
  !   print*, 'n surf:', nsurf
  !   Loop through surface types at the beginning of a day------------------------
  IF ( it==0 .AND. imin==5 ) THEN
     !       print*, '----- AnOHM called -----'
     !       print*, 'Grid@id:', Gridiv, id
     ! ------Set to zero initially------
     a1AnOHM(Gridiv) = 0   ![-]
     a2AnOHM(Gridiv) = 0   ![h]
     a3AnOHM(Gridiv) = 0   ![W m-2]
     !----------------------------------

     DO is=1,nsurf
        surfrac=sfr(is)
        !   print*, 'surfrac of ', is, 'is: ',surfrac
        !   initialize the coefs.
        xa1 = 0.1
        xa2 = 0.2
        xa3 = 10
        !   call AnOHM to calculate the coefs.
        IF ( is<nsurf ) THEN
           !  land surfaces
           CALL AnOHM_coef(is,id,Gridiv,xa1,xa2,xa3)
        ELSE
           !  water surface
           CALL AnOHM_coef_water(nsurf,id,Gridiv,xa1,xa2,xa3)
        END IF
        !   calculate the areally-weighted OHM coefficients
        a1AnOHM(Gridiv) = a1AnOHM(Gridiv)+surfrac*xa1
        a2AnOHM(Gridiv) = a2AnOHM(Gridiv)+surfrac*xa2
        a3AnOHM(Gridiv) = a3AnOHM(Gridiv)+surfrac*xa3

     ENDDO

     PRINT*, 'grid coeff:'
     PRINT*, a1AnOHM(Gridiv),a2AnOHM(Gridiv),a3AnOHM(Gridiv)
     !   end of loop over surface types -----------------------------------------
     !  IF ( id>365 ) THEN
     !     PRINT*, '----- OHM coeffs -----'
     !     WRITE(*,'(3f10.4)') a1AnOHM(Gridiv),a2AnOHM(Gridiv),a3AnOHM(Gridiv)
     !  END IF

  END IF

  !   Calculate radiation part ------------------------------------------------------------
  qs=NAN          !qs  = Net storage heat flux  [W m-2]
  IF(qn1>-999) THEN   !qn1 = Net all-wave radiation [W m-2]
     ! Old calculations
     !dqndt = 0.5*(qn1-q2_grids(Gridiv))*nsh_real/2       !gradient at t-1
     !!       Calculate net storage heat flux
     !qs = qn1*a1AnOHM(Gridiv)+dqndt*a2AnOHM(Gridiv)+a3AnOHM(Gridiv)          !Eq 4, Grimmond et al. 1991
     !
     !q1_grids(Gridiv) = q2_grids(Gridiv) !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
     !q2_grids(Gridiv) = q3_grids(Gridiv) !q2 = net radiation at t-1
     !q3_grids(Gridiv) = qn1              !q3 = net radiation at t (at t-1 when q3 used in next timestep)
     
     ! New calculations (HCW Dec 2016)
     ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
     if(nsh > 1) then
        qn1_store(1:(nsh-1),Gridiv) = qn1_store(2:nsh,Gridiv)    
        qn1_store(nsh,Gridiv) = qn1
        nsh_nna = sum(spread(1,1,nsh), mask=qn1_store(:,Gridiv) /= -999) !Find how many are not -999s  !bug fixed HCW 08 Feb 2017    
        qn1_av = sum(qn1_store(:,Gridiv), mask=qn1_store(:,Gridiv) /= -999)/nsh_nna
     elseif(nsh==1) then
         qn1_store(:,Gridiv) = qn1
         qn1_av = qn1
     endif         
     ! Store hourly average values (calculated every timestep) for previous 2 hours
     if(nsh > 1) then
        qn1_av_store(1:(2*nsh),Gridiv) = qn1_av_store(2:(2*nsh+1),Gridiv)    
        qn1_av_store(2*nsh+1,Gridiv) = qn1_av
     elseif(nsh==1) then
        qn1_av_store(:,Gridiv) = qn1_av
     endif    
     ! Calculate dQ* per dt for 60 min (using running mean Q* at t hours and (t-2) hours)
     if(any(qn1_av_store == -999)) then
        dqndt=0  ! Set dqndt term to zero for spinup
     else
        dqndt=0.5*(qn1_av_store((2*nsh+1),Gridiv)-qn1_av_store(1,Gridiv))
     endif
     
     ! Calculate net storage heat flux
     qs = qn1*a1AnOHM(Gridiv)+dqndt*a2AnOHM(Gridiv)+a3AnOHM(Gridiv)
  
  ELSE
     CALL ErrorHint(21,'SUEWS_AnOHM.f95: bad value for qn found during qs calculation.',qn1,NotUsed,notUsedI)
  ENDIF

END SUBROUTINE AnOHM_v2016
!========================================================================================


!========================================================================================
SUBROUTINE AnOHM_coef(sfc_typ,xid,xgrid,&   ! input
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
  ! the AnOHM paper to be added.
  !
  ! history:
  ! 20160222: initial version
  !========================================================================================
  USE allocateArray
  USE data_in
  USE defaultNotUsed
  USE gis_data
  USE sues_data
  USE time

  IMPLICIT NONE

  !   input
  INTEGER:: sfc_typ, xid, xgrid

  !   output
  REAL :: xa1, xa2, xa3

  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K
  REAL, PARAMETER :: CRA   = 915.483          ! converting RA (aerodyn. res.) to bulk trsf. coeff., [kg s-3]

  !   sfc. properties:
  REAL :: xalb,   &    !  albedo,
       xemis,  &    !  emissivity,
       xcp,    &    !  heat capacity,
       xk,     &    !  thermal conductivity,
       xch,    &    !  bulk transfer coef.
       xBo          !  Bowen ratio

  !   forcings:
  REAL, DIMENSION(24) :: Sd,& !   incoming solar radiation
       Ta,& !   air temperature
       WS,& !   wind speed
       WF,& !   anthropogenic heat
       AH   !   water flux density


  !   local variables:
  REAL(KIND(1d0))    :: beta               ! inverse Bowen ratio
  REAL(KIND(1d0))    :: f,fL,fT            ! energy redistribution factors
  REAL(KIND(1d0))    :: lambda             ! thermal diffusivity
  REAL(KIND(1d0))    :: delta,m,n          ! water flux related variables
  REAL(KIND(1d0))    :: ms,ns              ! m, n related
  REAL(KIND(1d0))    :: gamma              ! phase lag scale
  REAL(KIND(1d0))    :: ASd,mSd            ! solar radiation
  REAL(KIND(1d0))    :: ATa,mTa            ! air temperature
  REAL(KIND(1d0))    :: tau                ! phase lag between Sd and Ta (Ta-Sd)
  REAL(KIND(1d0))    :: ATs!,mTs            ! surface temperature amplitude
  REAL(KIND(1d0))    :: ceta,cphi          ! phase related temporary variables
  REAL(KIND(1d0))    :: eta,phi,xlag       ! phase related temporary variables
  REAL(KIND(1d0))    :: mWS,mWF,mAH        ! mean values of WS, WF and AH
  REAL(KIND(1d0))    :: xx1,xx2,xx3        ! temporary use
  !REAL(KIND(1d0))    :: solTs              ! surface temperature
  LOGICAL :: flagGood = .TRUE.  ! quality flag, T for good, F for bad

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

  !   load met. forcing data:
  CALL AnOHM_FcLoad(sfc_typ,xid,xgrid,&   ! input
       Sd,Ta,WS,WF,AH)         ! output
  !   print*, 'here the forcings:'
  !   print*, Sd,Ta,WS,WF,AH

  !   load forcing characteristics:
  CALL AnOHM_FcCal(Sd,Ta,WS,WF,AH,&               ! input
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)  ! output
  !   print*, ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH
  ! if ( sfc_typ==1 ) then
  ! print*, 'AH of sfc',sfc_typ,':', mAH
  ! end if
  !   load sfc. properties:
  CALL AnOHM_SfcLoad(sfc_typ,xid,xgrid,&          ! input
       xalb,xemis,xcp,xk,xch,xBo)  ! output
  !   print*, 'here the properties:'
  !   print*, xalb,xemis,xcp,xk,xch,xBo


  !   initial Bowen ratio
  BoAnOHMStart(xgrid) = xBo

  !   calculate sfc properties related parameters:
  xch = xch*mWS
  xch = 0.5*xch

  !     xch    = CRA/RA !update bulk trsf. coeff. with RA (aerodyn. res.)
  beta   = 1/xBo
  ! PRINT*, 'beta:', beta
  f      = ((1+beta)*xch+4*SIGMA*xemis*mTa**3)
  !     print*, 'xch',xch,'mTa',mTa,'dLu',4*SIGMA*xemis*mTa**3
  fL     = 4*SIGMA*xemis*mTa**3
  !     print*, 'fL',fL
  fT     = (1+beta)*xch
  !     print*, 'fT',fT
  lambda = xk/xcp
  !     print*, 'lambda',lambda
  delta  = SQRT(.5*(mWF**2+SQRT(mWF**4+16*lambda**2*OMEGA**2)))
  !     print*, 'delta',delta
  m      = (2*lambda)/(delta+mWF)
  n      = delta/OMEGA
  !     print*, 'm',m,'n',n


  !   calculate surface temperature related parameters:
  !   mTs   = (mSd*(1-xalb)/f)+mTa
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
     IF ( xa2>0.5) xa2 = 0.5
     IF (xa2<-0.5) xa2 = -0.5
  ENDIF
  !   quality checking of a3
  IF ( .NOT. (xa3<0)) flagGood = .FALSE.

  !   skip the first day for quality checking
  IF ( xid == 1 ) flagGood = .TRUE.


  !     if ( flagGood ) then
  !         print*, 'a1,a2,a3:', xa1,xa2,xa3
  !     else
  !         print*, 'a1,a2,a3:'
  !         print*, 'today:', xa1,xa2,xa3
  !         print*, 'yesterday:', a123AnOHM_gs(xid-1,xgrid,sfc_typ,:)
  ! !         if ( xid>20 ) stop 'bad values met!'
  !     end if



  !   update good values of AnOHM coefs.:
  ! IF ( .NOT. flagGood ) THEN
  !       replace values of today of yesterday
  ! xa1 = a123AnOHM_gs(xid-1,xgrid,sfc_typ,1)
  ! xa2 = a123AnOHM_gs(xid-1,xgrid,sfc_typ,2)
  ! xa3 = a123AnOHM_gs(xid-1,xgrid,sfc_typ,3)
  ! ENDIF

  !   update today's values with good values
  a123AnOHM_gs(xid,xgrid,sfc_typ,1) = xa1
  a123AnOHM_gs(xid,xgrid,sfc_typ,2) = xa2
  a123AnOHM_gs(xid,xgrid,sfc_typ,3) = xa3



  !     print*, 'ceta,cphi', ceta,cphi
  !     print*, 'tau,eta,phi,xlag in deg:',tau/pi*180,eta/pi*180,phi/pi*180,xlag/pi*180

  ! PRINT*, '********sfc_typ: ',sfc_typ,' end********'

END SUBROUTINE AnOHM_coef
!========================================================================================

!========================================================================================
SUBROUTINE AnOHM_coef_water(sfc_typ,xid,xgrid,&   ! input
     xa1,xa2,xa3)            ! output
  ! author: Ting Sun
  ! date: 20161124
  !
  ! purpose:
  ! designed for water surface
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
  ! the AnOHM paper to be added.
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
  INTEGER:: sfc_typ, xid, xgrid

  !   output
  REAL :: xa1, xa2, xa3

  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K

  !   sfc. properties:
  REAL :: xalb,   &    !  albedo,
       xemis,  &    !  emissivity,
       xcp,    &    !  heat capacity,
       xk,     &    !  thermal conductivity,
       xch,    &    !  bulk transfer coef.
       xBo,    &    !  Bowen ratio
       xeta,  &     ! effective absorption coefficient
       xmu         ! effective absorption fraction

  !   forcings:
  REAL, DIMENSION(24) :: Sd,& !   incoming solar radiation
       Ta,& !   air temperature
       WS,& !   wind speed
       WF,& !   anthropogenic heat
       AH   !   water flux density


  !   local variables:
  REAL :: beta                   ! inverse Bowen ratio
  REAL :: f,fL,fT                ! energy redistribution factors
  REAL :: lambda,calb            ! temporary use
  REAL :: delta,m,n              ! sfc. temperature related variables
  REAL :: xm,xn                  ! m, n related
  REAL :: gamma,phi              ! phase lag scale
  REAL :: ASd,mSd                ! solar radiation
  REAL :: ATa,mTa                ! air temperature
  REAL :: tau                    ! phase lag between Sd and Ta (Ta-Sd)
  REAL :: ATs,mTs                ! surface temperature amplitude
  REAL :: czeta,ctheta           ! phase related temporary variables
  REAL :: zeta,theta,xlag           ! phase related temporary variables
  REAL :: mWS,mWF,mAH            ! mean values of WS, WF and AH
  REAL :: xx1,xx2,xx3            ! temporary use
  REAL :: kappa                  ! temporary use
  REAL :: dtau,dpsi,dphi         ! temporary use
  REAL :: cdtau,cdpsi,cdphi      ! temporary use
  REAL :: xxT,xxkappa,xxdltphi   ! temporary use
  LOGICAL :: flagGood = .TRUE.  ! quality flag, T for good, F for bad


  ! !   give fixed values for test
  !
  ! !   properties
  ! xalb  = .05
  ! xemis = .95
  ! xcp   = 2e6
  ! xk    = 1.5
  ! xch   = 2
  ! xBo   = 1/3.
  xeta  = 0.3
  xmu   = 0.2
  !
  !
  ! !   forcings
  ! ASd = 400
  ! mSd = 200
  ! ATa = 5
  ! mTa = 25+C2K
  ! tau = PI/6
  ! WS  = 2.1
  ! AH  = 0
  ! Wf  = 0

  !   load met. forcing data:
  CALL AnOHM_FcLoad(sfc_typ,xid,xgrid,&   ! input
       Sd,Ta,WS,WF,AH)         ! output
  !   write(*,*) 'here the forcings:'
  !   write(*,*) Sd,Ta,WS,WF,AH

  !   write(unit=*, fmt=*) 'DOY:', xid

  ! load forcing characteristics:
  CALL AnOHM_FcCal(Sd,Ta,WS,WF,AH,&               ! input
       ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)  ! output
  !   write(*,*) ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH
  !   load sfc. properties:
  CALL AnOHM_SfcLoad(sfc_typ,xid,xgrid,&          ! input
       xalb,xemis,xcp,xk,xch,xBo)  ! output
  !   write(*,*) 'here the properties:'
  !   write(*,*) xalb,xemis,xcp,xk,xch,xBo

  !   initial Bowen ratio
  BoAnOHMStart(xgrid) = xBo

  !   calculate sfc properties related parameters:
  xm     = xk*xmu**2
  xn     = xcp*OMEGA
  phi    = ATAN(xn/xm)
  kappa  = SQRT(xcp*OMEGA/(2*xk))

  mWS=SUM(WS, dim=1, mask=(WS>0))/SIZE(WS, dim=1)
  xch    = xch*mWS
  beta   = 1/xBo
  f      = ((1+beta)*xch+4*SIGMA*xemis*mTa**3)
  fL     = 4*SIGMA*xemis*mTa**3
  fT     = (1+beta)*xch

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
     IF ( xa2>0.5) xa2 = 0.5
     IF (xa2<-0.5) xa2 = -0.5
  ENDIF
  !   quality checking of a3
  IF ( .NOT. (xa3<0)) flagGood = .FALSE.

  !   skip the first day for quality checking
  IF ( xid == 1 ) flagGood = .TRUE.

  WRITE(*,*) 'sfc_typ_water:', sfc_typ
  WRITE(*,*) 'a1,a2,a3:', xa1,xa2,xa3

END SUBROUTINE AnOHM_coef_water
!========================================================================================

!========================================================================================
SUBROUTINE AnOHM_FcCal(Sd,Ta,WS,WF,AH,&                 ! input
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
  !
  ! output:
  ! 1) ASd, mSd: amplitude and mean value of Sd
  ! 2) ATa, mTa: amplitude and mean value of Ta
  ! 3) tau: phase lag between Sd and Ta
  ! 4) mWS: mean value of WS
  ! 4) mWF: mean value of WF
  ! 4) mAH: mean value of AH
  !
  ! ref:
  !
  ! history:
  ! 20160224: initial version
  !========================================================================================
  IMPLICIT NONE

  !   input
  REAL :: Sd(24),&    !
       Ta(24),&    !
       WS(24),&    !
       Wf(24),&    !
       AH(24)      !

  !   output
  REAL :: ASd,mSd,&   ! Sd scales
       ATa,mTa,&   ! Ta scales
       tau,&       ! phase lag between Sd and Ta
       mWS,&       ! mean WS
       mWF,&       ! mean WF
       mAH         ! mean AH

  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K

  !   local variables:
  REAL(KIND(1d0)) :: tSd,tTa         ! peaking timestamps
  ! REAL(KIND(1d0)) :: aCosb,aSinb,b,c ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
  ! REAL(KIND(1d0)) :: xx              ! temporary use


  !   calculate sinusoidal scales of Sd:
  !     print*, 'Calc. Sd...'
  CALL AnOHM_ShapeFit(Sd(10:16),ASd,mSd,tSd)
  !   modify ill-shaped days to go through
  !     if ( ASd < 0 ) then
  !         ASd = abs(ASd)
  !         tSd = 12 ! assume Sd peaks at 12:00LST
  !     end if
  !   print*, 'ASd:', ASd
  !   print*, 'mSd:', mSd
  !   print*, 'tSd:', tSd

  !   calculate sinusoidal scales of Ta:
  !     print*, 'Calc. Ta...'
  CALL AnOHM_ShapeFit(Ta(10:16),ATa,mTa,tTa)
  IF ( mTa < 60 ) mTa = mTa+C2K ! correct the Celsius to Kelvin
  !   modify ill-shaped days to go through
  !     if ( ATa < 0 ) then
  !         ATa = abs(ATa)
  !         tTa = 14 ! assume Ta peaks at 14:00LST
  !     end if
  !   print*, 'ATa:', ATa
  !   print*, 'mTa:', mTa
  !   print*, 'tTa:', tTa

  !   calculate the phase lag between Sd and Ta:
  tau = (tTa-tSd)/24*2*PI
  !   print*, 'tau:', tau

  !   calculate the mean values:
  mWS = SUM(WS(10:16))/7  ! mean value of WS
  mWF = SUM(WF(10:16))/7  ! mean value of WF
  mAH = SUM(AH(10:16))/7  ! mean value of AH
  !     print*, 'mWS:', mWS
  !     print*, 'mWF:', mWF
  !     print*, 'mAH:', mAH


END SUBROUTINE AnOHM_FcCal
!========================================================================================

!========================================================================================
SUBROUTINE AnOHM_ShapeFit(obs,amp,mean,tpeak)
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
  !========================================================================================
  IMPLICIT NONE

  !   input
  REAL :: obs(7)  ! observation (daytime 10:00–16:00)

  !   output
  REAL :: amp     ! amplitude
  REAL :: mean    ! average
  REAL :: tpeak   ! peaking time (h)

  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K

  !   local variables:
  REAL    :: coefs(3,7)           ! coefficients for least squre solution
  REAL    :: aCosb,aSinb,a,b,c    ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
  REAL    :: xx,mbias,mbiasObs    ! temporary use
  INTEGER :: i                    ! temporary use

  !   c coeffs.:
  !   exact values:
  !   coefs(1,:) = (/2*(-615 - 15*Sqrt(2.) + 211*Sqrt(3.) + 80*Sqrt(6.)), &
  !                 3*(-414 - 12*Sqrt(2.) + 140*Sqrt(3.) + 75*Sqrt(6.)),      &
  !                 -1437 + 156*Sqrt(2.) + 613*Sqrt(3.) + 92*Sqrt(6.),        &
  !                 2*(-816 + 177*Sqrt(2.) + 337*Sqrt(3.) + 13*Sqrt(6.)), &
  !                 -1437 + 156*Sqrt(2.) + 613*Sqrt(3.) + 92*Sqrt(6.),        &
  !                 3*(-414 - 12*Sqrt(2.) + 140*Sqrt(3.) + 75*Sqrt(6.)),      &
  !                 2*(-615 - 15*Sqrt(2.) + 211*Sqrt(3.) + 80*Sqrt(6.))/)
  !   coefs(1,:) = coefs(1,:)/(-9450 + 534*Sqrt(2.) + 3584*Sqrt(3.) + 980*Sqrt(6.))
  !   apprx. values:
  coefs(1,:) = (/1.72649, 0.165226, -0.816223, -1.15098, -0.816223, 0.165226, 1.72649/)
  c          = dot_PRODUCT(coefs(1,:),obs)

  !   aCos coeffs.:
  !   exact values:
  !   coefs(2,:) = (/249*(-6 + Sqrt(2.)) + Sqrt(3.)*(124 + 347*Sqrt(2.)), &
  !                  -3*(-263 - 372*Sqrt(2.) + Sqrt(3.) + 325*Sqrt(6.)),      &
  !                  81 - 48*Sqrt(2.) - 247*Sqrt(3.) + 174*Sqrt(6.),          &
  !                  3*(27 - 545*Sqrt(2.) - 45*Sqrt(3.) + 340*Sqrt(6.)),      &
  !                  1740 - 303*Sqrt(2.) - 632*Sqrt(3.) - 73*Sqrt(6.),        &
  !                  -207 + 1368*Sqrt(2.) - 505*Sqrt(3.) - 338*Sqrt(6.),      &
  !                  -990 - 747*Sqrt(2.) + 1398*Sqrt(3.) - 155*Sqrt(6.)/)
  !   coefs(2,:) = coefs(2,:)/(-9450 + 534*Sqrt(2.) + 3584*Sqrt(3.) + 980*Sqrt(6.))
  !   apprx. values:
  coefs(2,:) = (/0.890047, 0.302243, -0.132877, -0.385659, -0.438879, -0.288908, 0.054033/)
  aCosb      = dot_PRODUCT(coefs(2,:),obs)

  !   cSin coeffs.:
  !   exact values:
  !   coefs(3,:) = (/-28 - 13*Sqrt(2.) - 5*Sqrt(3.)*(-6 + Sqrt(2.)),  &
  !                  -15 + Sqrt(2.) - 9*Sqrt(3.) + 12*Sqrt(6.),       &
  !                  63 - 7*Sqrt(2.) - 21*Sqrt(3.) - 5*Sqrt(6.),      &
  !                  -9 - 7*Sqrt(3.) + 11*Sqrt(6.),                   &
  !                  -32 - 7*Sqrt(2.) + 28*Sqrt(3.) - Sqrt(6.),       &
  !                  -3 + 27*Sqrt(2.) - 5*Sqrt(3.) - 11*Sqrt(6.), &
  !                  24 - Sqrt(2.) - 16*Sqrt(3.) - Sqrt(6.)/)
  !   coefs(3,:) = coefs(3,:)/(-224 + 36*Sqrt(2.) + 58*Sqrt(3.) + 28*Sqrt(6.))
  !   apprx. values:
  coefs(3,:) = (/1.64967, -0.0543156, -1.10791, -1.4393, -1.02591, 0.104083, 1.87368/)
  aSinb      = dot_PRODUCT(coefs(3,:),obs)

  !   calculate b and a:
  b = ATAN(aSinb/aCosb)
  IF ( b>0 ) b = b-Pi     ! handle over Pi/2 phase lag
  a = aSinb/SIN(b)

  mbias    = 0.
  mbiasObs = 0.
  !     write(*, *) '      obs        ','     sim       ','     diff     '
  DO i = 1, 7, 1
     !     print out the sim and obs pairs
     xx       = a*SIN(Pi/12*(9+i)+b)+c
     mbias    = mbias+ABS((xx-obs(i)))
     mbiasObs = mbiasObs+ABS(obs(i))
     !         write(*, *) obs(i),xx,xx-obs(i)
  END DO
  mbias = mbias/mbiasObs*100
  !     write(*, *) 'mean bias (%): ', mbias


  !   get results:
  amp=a                   ! amplitude
  mean=c                  ! mean value
  tpeak=(PI/2-b)/PI*12    ! convert to timestamp (h)

  IF ( a<0 ) THEN
     !         print*, 'mean bias (%): ', mbias
     !         print*, 'obs: ', obs
     !         stop
  ENDIF


END SUBROUTINE AnOHM_ShapeFit
!========================================================================================

!========================================================================================
SUBROUTINE AnOHM_FcLoad(sfc,xid,xgrid,& ! input
     Sd,Ta,WS,WF,AH) ! output
  ! author: Ting Sun
  !
  ! purpose:
  ! calculate the key parameters of a sinusoidal curve for AnOHM forcings
  ! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
  !
  ! input:
  ! 1) sfc: surface type
  ! 2) xid : day of year
  !
  ! output:
  ! hourly values between 00:00 and 23:00 (local time, inclusive):
  ! 1) Sd: incoming shortwave radiation,  W m-2
  ! 2) Ta: air temperature, K
  ! 3) WS: wind speed, m s-1
  ! 4) WF: water flux density, ???
  ! 5) AH: anthropogenic heat, W m-2
  !
  ! ref:
  !
  !
  ! todo:
  ! check the completeness of forcings of given day (i.e., xid)
  !
  ! history:
  ! 20160226: initial version
  !========================================================================================
  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE initial
  USE sues_data
  USE time

  IMPLICIT NONE

  !   input
  INTEGER :: sfc, xid, xgrid

  !   output
  REAL :: Sd(24),&    !
       Ta(24),&    !
       WS(24),&    !
       WF(24),&    !
       AH(24)      !

  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K

  !   local variables:
  ! REAL :: tSd,tTa         ! peaking timestamps
  ! REAL :: aCosb,aSinb,b,c ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
  ! REAL :: xx              ! temporary use

  INTEGER :: i            ! temporary use


  INTEGER :: irRange(2)   ! row range in MetData containing id values
  INTEGER :: lenMetData

  INTEGER, ALLOCATABLE :: lSub(:) ! array to retrieve data of the specific day (id)

  lenMetData = SIZE(MetForcingData(:,2,xgrid))

  ! determine the positions of met forcings of xid in the whole met chunk
  ALLOCATE(lSub(lenMetData))
  lSub = (/(i,i=1,lenMetData)/)
  ! select pieces of xid
  WHERE (MetForcingData(:,2,xgrid) == xid)
     lSub = 1
  ELSEWHERE
     lSub = 0
  END WHERE
  irRange = (/MAXLOC(lSub),MAXLOC(lSub)+nsd-1/)
  ! PRINT*, 'irRange:', irRange

  ! IF ( sfc==1 .AND. xid ==366 ) THEN
  !    PRINT*, 'irRange:', irRange
  !    PRINT*, 'id:',xid
  !   !  print*, 'id:', MetForcingData(:,2,xgrid)
  !
  ! END IF



  !   load the sublist into forcings:
  Sd = MetForcingData(irRange(1):irRange(2):nsh,15,xgrid)     ! nsh is timesteps per hour, a global var.
  Ta = MetForcingData(irRange(1):irRange(2):nsh,12,xgrid)
  WS = MetForcingData(irRange(1):irRange(2):nsh,10,xgrid)
  WF = MetForcingData(irRange(1):irRange(2):nsh,12,xgrid)*0   ! set as 0 for debug
  IF ( AnthropHeatMethod == 0 ) THEN
     AH = MetForcingData(irRange(1):irRange(2):nsh,9,xgrid)*0    ! read in from MetForcingData,
  ELSE
     AH = mAH_grids(xid-1,xgrid)
  END IF

  !   print*, 'Sd:', Sd
  !   print*, 'Ta:', Ta
  !   print*, 'WS:', WS
  !   print*, 'WF:', WF
  !   print*, 'AH:', AH

END SUBROUTINE AnOHM_FcLoad
!========================================================================================

!========================================================================================
SUBROUTINE AnOHM_SfcLoad(sfc,xid,xgrid,&            ! input
     xalb,xemis,xcp,xk,xch,xBo)  ! output

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
  ! ref:
  !
  ! history:
  ! 20160228: initial version
  !========================================================================================
  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE initial
  USE sues_data
  USE time
  USE InitialCond

  IMPLICIT NONE

  !   input
  INTEGER :: sfc, xid, xgrid

  !   output:
  !   sfc. properties:
  REAL :: xalb,       &    !  albedo,
       xemis,      &    !  emissivity,
       xcp,        &    !  heat capacity,
       xk,         &    !  thermal conductivity,
       xch,        &    !  bulk transfer coef.
       xBo              !  Bowen ratio
  !   constant
  REAL, PARAMETER :: SIGMA = 5.67e-8          ! Stefan-Boltzman
  REAL, PARAMETER :: PI    = ATAN(1.0)*4      ! Pi
  REAL, PARAMETER :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
  REAL, PARAMETER :: C2K   = 273.15           ! degC to K


  !   load properties from global variables
  xalb  = alb(sfc)
  xemis = emis(sfc)
  xcp   = cpAnOHM(sfc)
  xk    = kkAnOHM(sfc)
  xch   = chAnOHM(sfc)

  !   print*, 'xalb:', xalb
  !   print*, 'xemis:', xemis
  !   print*, 'sfc:', sfc
  !   print*, 'xcp:', xcp
  !   print*, 'xk:', xk
  !   print*, 'xch:', xch


  !   load Bowen ratio of yesterday from DailyState calculation
  IF ( BoInit==NAN ) THEN
     CALL ErrorHint(68,'No initial Bowen ratio found in InitialConditions.',BoInit,notUsed,notUsedI)
  ELSE
     xBo = BoInit ! get Initial value from initial condition namelist
  END IF
  ! set as previous day's Bo, 20160708 TS
  xBo = Bo_grids(xid-1,xgrid) ! default Bo will be read in when xid = 0
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
