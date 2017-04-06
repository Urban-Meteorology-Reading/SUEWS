! Note: INTERVAL is now set to 3600 s in Initial (it is no longer set in RunControl) HCW 29 Jan 2015
! Last modified:
!  HCW 29 Mar 2017 - Changed third dimension of dataOutBL to Gridiv (was previously iMB which seems incorrect)
!  NT 6 Apr 2017 - include top of the CBL variables in RKUTTA scheme + add flag to include or exclude subsidence
!  LJ 27 Jan 2016 - Removal of tabs

SUBROUTINE CBL(ifirst,iMB,Gridiv)

  USE mod_z
  USE mod_k
  USE gas
  USE time
  USE data_in
  USE sues_data
  USE moist
  USE allocateArray
  USE defaultNotUsed
  USE cbl_module
  USE gis_data
  USE WhereWhen
  
  IMPLICIT NONE

  REAL(KIND(1d0))::sat_vap_press
  REAL(KIND(1d0))::qh_use,qe_use,tm_K_zm,qm_gkg_zm
  REAL(KIND(1d0))::Temp_C1,avrh1,es_hPa1
  REAL(KIND(1d0))::secs0,secs1,Lv
  INTEGER::idoy,ifirst,iMB,Gridiv,startflag
  REAL(KIND(1d0)), PARAMETER::pi=3.141592653589793d+0,d2r=pi/180.

  ! Reset iCBLcount at start of each metblock (HCW added 29/03/2017)
  IF(ifirst == 1) THEN
     iCBLcount = 0
  ENDIF
  
  !Skip first loop and unspecified days
  !IF((ifirst==1 .AND. iMB==1) .OR. CBLday(id)==0) THEN   !HCW modified condition to check for first timestep of the model run
  IF(ifirst==1 .OR. CBLday(id)==0) THEN   !HCW modified 29/03/2017
     iCBLcount=iCBLcount+1
  !write(*,*) 'ifirst or nonCBLday', DateTime, iCBLcount   
     dataOutBL(iCBLcount,1:ncolumnsdataOutBL,Gridiv)=(/REAL(iy,8),REAL(id,8),REAL(it,8),REAL(imin,8),dectime, &
                                                    (NAN,is=6,ncolumnsdataOutBL)/)
     RETURN
  ELSEIF(avkdn<5)THEN
     CALL CBL_initial(qh_use,qe_use,tm_K_zm,qm_gkg_zm,startflag,iMb, Gridiv)
     RETURN
  ENDIF

  IF(startflag==0)THEN !write down initial values in previous time step
  !write(*,*) 'startflag', DateTime, iCBLcount   
     dataOutBL(iCBLcount,1:ncolumnsdataOutBL,Gridiv)=(/REAL(iy,8),REAL(id,8),REAL(it,8),REAL(imin,8),dectime,blh_m,tm_K, &
               qm_kgkg*1000,tp_K,qp_kgkg*1000,(NAN,is=11,20),gamt_Km,gamq_kgkgm/)
     startflag=1
  ENDIF

  qh_use=qhforCBL(Gridiv)   !HCW 21 Mar 2017
  qe_use=qeforCBL(Gridiv)
  IF(qh_use<-900.OR.qe_use<-900)THEN  ! observed data has a problem
     CALL ErrorHint(22,'Unrealistic qh or qe_value for CBL.',qh_use,qe_use,qh_choice)
  ENDIF
  !!Heat flux choices - these are now made in SUEWS_Calculations for qhforCBL and qeCBL, rather than here
  !IF(Qh_choice==1) THEN   !from SUEWS
  !  !qh_use=qh
  !   !qe_use=qeph
  !   qh_use=qhforCBL(Gridiv)   !HCW 21 Mar 2017
  !   qe_use=qeforCBL(Gridiv)
  !ELSEIF(qh_choice==2)THEN !from LUMPS
  !   qh_use=H_mod
  !   qe_use=E_mod
  !ELSEIF(qh_choice==3)THEN  !from OBS
  !   IF(qh_obs<-900.OR.qe_obs<-900)THEN  ! observed data has a problem
  !      CALL ErrorHint(22,'Unrealistic observed qh or qe_value.',qh_obs,qe_obs,qh_choice)
  !   ENDIF
  !   qh_use=qh_obs
  !   qe_use=qe_obs
  !ENDIF

  !-------Main loop of CBL calculation--------------------------------------
  !-------------------------------------------------------------------------

  cbldata(1)=float(it)+float(imin)/60.
  cbldata(2)=qh_use
  cbldata(3)=qe_use
  cbldata(4)=avdens
  cbldata(5)=lv_J_kg
  cbldata(6)=avcp
  cbldata(7)=avu1
  cbldata(8)=ustar
  cbldata(9)=Press_hPa
  cbldata(10)=psyh

  secs0=cbldata(1)*3600.
  secs1=secs0+float(tstep) ! time in seconds
  ! Kinematic fluxes
  fhbl_Kms    = cbldata(2)/ (cbldata(4)*cbldata(6))  !qh_use/(avdens*avcp)      ! units: degK * m/s
  febl_kgkgms = cbldata(3)/ (cbldata(4)*cbldata(5))  !qe_use/(avdens*lv_J_kg)   ! units: kg/kg * m/s
  IF(CO2_included==1)THEN
     fcbl = 0!fc(i)/(rmco2/volm)      ! units: mol/mol * m/s
  ELSE
     cm=NAN
  ENDIF

!   tpp_K=tp_K
!   qpp_kgkg=qp_kgkg

  IF(sondeflag.EQ.1) THEN
     CALL gamma_sonde
  ENDIF
  !     	set up array for Runge-Kutta call
  blh1_m=blh_m
  y(1)=blh_m ! integrate h, t, q, c from time s(i-1)
  y(2)=tm_K  ! to time s(i) using runge-kutta solution
  y(3)=qm_kgkg   ! of slab CBL equations
  y(4)=cm
  y(5)=tp_K
  y(6)=qp_kgkg

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL rkutta(neqn,secs0,secs1,y,1)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++
  blh_m   =y(1)
  tm_K    =y(2)  ! potential temperature, units: deg C  <-NT: shouldn't this be K?
  qm_kgkg =y(3)  ! specific humidity, units: kg/kg
  cm      =y(4)  ! co2 concentration,units: mol/mol
  tp_K    =y(5)  ! potential temperature top of CBL: K
  qp_kgkg =y(6)  ! specific humidity top of CBL: kg/kg
  !     	compute derived quantities for this time step

!NT: now included in rkutta
!   tp_K   = tpp_K     + (gamt_Km*(blh_m-blh1_m))
!   qp_kgkg = qpp_kgkg + (gamq_kgkgm*(blh_m-blh1_m))
! 
!   IF (tp_K.LT.tm_K) THEN
!      tp_K = tm_K
!   ENDIF

  !		th = tm_K - (grav/cbldata(5))*blh_m			 ! actual temp just below z=h
  !		dh = qsatf(th,cbldata(8)) - qm_kgkg           ! deficit just below z=h

  tp_C=tp_K-C2K
  tm_C=tm_K-C2K

  !	 delt = tp_K - tm_K ! temp
  !	 delq = qp_kgkg - qm_kgkg ! humidity
  !deltv = (tp_K - tm_K) + 0.61*tm_k*(qp_kgkg - qm_kgkg)  ! pot virtual temp

  qm_gkg=qm_kgkg*1000 !humidities: kg/kg -> g/kg

  !Output time correction
  idoy=id
  !If(it==0 .and. imin==55) idoy=id-1
  IF(it==0 .AND. imin==(nsh_real-1)/nsh_real*60) idoy=id-1  !Modified by HCW 04 Mar 2015 in case model timestep is not 5-min


  IF((qh_choice==1).OR.(qh_choice==2))THEN !BLUEWS or BLUMPS
     !Stability correction
     !tm_K_zm=tm_K+cbldata(10)*cbldata(2)/(k*cbldata(8)*cbldata(6)*cbldata(4))
     Temp_C=tm_K/((1000/cbldata(9))**(gas_ct_dry/cbldata(6)))-C2K
     es_hPa=sat_vap_press(Temp_C,cbldata(9),1)
     lv=(2500.25-2.365*Temp_C)*1000
     !qm_gkg_zm=qm_gkg+cbldata(10)*cbldata(3)/(k*cbldata(8)*cbldata(4)*lv)
     avrh=100*((qm_gkg*cbldata(9)/(622+qm_gkg))/es_hPa) !check pressure
     IF(avrh>100)THEN
        CALL errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbldata(1)/24.0,avrh,100)
        avrh=100
     ENDIF
     iCBLcount=iCBLcount+1
     !write(*,*) 'qh1or2', DateTIme, iCBLcount   
     dataOutBL(iCBLcount,1:ncolumnsdataOutBL,Gridiv)=(/REAL(iy,8),REAL(id,8),REAL(it,8),REAL(imin,8),dectime,blh_m,tm_K, & 
                qm_kgkg*1000, tp_K,qp_kgkg*1000,&
          Temp_C,avrh,cbldata(2),cbldata(3),cbldata(9),cbldata(7),cbldata(8),cbldata(4),cbldata(5),cbldata(6),&
          gamt_Km,gamq_kgkgm/)
  ELSEIF(qh_choice==3)THEN ! CBL
     !tm_K_zm=tm_K+cbldata(10)*cbldata(2)/(k*cbldata(8)*cbldata(6)*cbldata(4))
     Temp_C1=tm_K/((1000/cbldata(9))**(gas_ct_dry/cbldata(6)))-C2K
     es_hPa1=sat_vap_press(Temp_C1,cbldata(9),1)
     lv=(2500.25-2.365*Temp_C1)*1000
     !qm_gkg_zm=qm_gkg+cbldata(10)*cbldata(3)/(k*cbldata(8)*cbldata(4)*lv)
     avrh1=100*((qm_gkg*cbldata(8)/(622+qm_gkg))/es_hPa1) !check pressure
     IF(avrh1>100)THEN
        CALL errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbldata(1)/24.0,avrh,100)
        avrh1=100
     ENDIF
     iCBLcount=iCBLcount+1
     !write(*,*) 'qh3', DateTIme, iCBLcount   
     dataOutBL(iCBLcount,1:ncolumnsdataOutBL,Gridiv)=(/REAL(iy,8),REAL(id,8),REAL(it,8),REAL(imin,8),dectime,blh_m,tm_K, &
          qm_kgkg*1000,tp_K,qp_kgkg*1000,&
          Temp_C1,avrh1,cbldata(2),cbldata(3),cbldata(9),cbldata(7),cbldata(8),cbldata(4),cbldata(5),cbldata(6),&
          gamt_Km,gamq_kgkgm/)
  ENDIF

  RETURN

END SUBROUTINE CBL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE CBL_ReadInputData
  USE allocateArray
  USE data_in
  USE sues_data
  USE cbl_module
  USE initial
  USE WhereWhen

  IMPLICIT NONE

  INTEGER::i, ios
  REAL(KIND(1d0))::l

  NAMELIST/CBLInput/EntrainmentType,&
       QH_choice,&
       CO2_included,&
       cblday,&
       wsb,&
       InitialData_use,&
       InitialDataFileName,&
       sondeflag,&
       FileSonde

  OPEN(51,file=TRIM(FileInputPath)//'CBLInput.nml',status='old', err=24)
  READ(51,nml=CBLInput,err=24)
  CLOSE(51)

  !Read initial values if it's needed
  IF(InitialData_use==1 .OR. InitialData_use==2)THEN
     OPEN(52,file=TRIM(FileInputPath)//TRIM(InitialDataFileName),status='old', err=25)
     READ(52,*)
     nlineInData = 0   !Initialise nlines
     DO
        READ(52,*, iostat=ios) l
        IF(ios<0 .or. l == -9) EXIT   !IF (l == -9) EXIT
        nlineInData = nlineInData + 1
     ENDDO
     CLOSE(52)
    
     ALLOCATE(IniCBLdata(1:nlineInData,1:8))
     OPEN(52,file=TRIM(FileInputPath)//TRIM(InitialDataFileName),status='old', err=25)
     READ(52,*)
     DO i=1,nlineInData
        READ(52,*)IniCBLdata(i,1:8)
     ENDDO
     CLOSE(52)
  ENDIF

  IF(CO2_included==0)THEN
     fcbl=0       ! hard-wire no CO2
  ENDIF

  iCBLcount=0

  RETURN

24 CALL ErrorHint(24,'CBLInput.nml',0.00D0,0.000D0,0)
25 CALL ErrorHint(24,TRIM(FileInputPath)//TRIM(InitialDataFileName),0.00D0,0.00D0,0)

END SUBROUTINE CBL_ReadInputData

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE CBL_initial(qh_use,qe_use,tm_K_zm,qm_gkg_zm,startflag,iMB, Gridiv)

  USE mod_z
  USE mod_k
  USE gas
  USE time
  USE data_in
  USE sues_data
  USE moist
  USE allocateArray
  USE defaultNotUsed
  USE cbl_module
  USE gis_data
  USE WhereWhen
  
  IMPLICIT NONE

  REAL(KIND(1d0))::qh_use,qe_use,tm_K_zm,qm_gkg_zm
  REAL(KIND(1d0))::qsatf,sat_vap_press,lv
  INTEGER::i,nLineDay,iMB,Gridiv,startflag

  
  qh_use=qhforCBL(Gridiv)   !HCW 21 Mar 2017
  qe_use=qeforCBL(Gridiv)
  IF(qh_use<-900.OR.qe_use<-900)THEN  ! observed data has a problem
     CALL ErrorHint(22,'Unrealistic qh or qe_value for CBL.',qh_use,qe_use,qh_choice)
  ENDIF
  !!Heat flux choices - these are now made in SUEWS_Calculations for qhforCBL and qeCBL, rather than here
  !IF(Qh_choice==1) THEN   !from SUEWS
  !   !qh_use=qh
  !   !qe_use=qeph
  !   qh_use=qhforCBL(Gridiv)   !HCW 21 Mar 2017
  !   qe_use=qeforCBL(Gridiv)
  !ELSEIF(qh_choice==2)THEN !from LUMPS
  !   qh_use=H_mod
  !   qe_use=E_mod
  !ELSEIF(qh_choice==3)THEN  !from OBS
  !   IF(qh_obs<-900.OR.qe_obs<-900)THEN  ! observed data has a problem
  !      CALL ErrorHint(22,'Unrealistic observed qh or qe_value.',qh_obs,qe_obs,qh_choice)
  !   ENDIF
  !   qh_use=qh_obs
  !   qe_use=qe_obs
  !ENDIF
  

  blh_m=NAN
  iCBLcount=iCBLcount+1
  !write(*,*) 'cblinitial', DateTIme, iCBLcount   
  dataOutBL(iCBLcount,1:ncolumnsdataOutBL,Gridiv)=(/REAL(iy,8),REAL(id,8),REAL(it,8),REAL(imin,8),dectime, &
                                                 (NAN,is=6,ncolumnsdataOutBL)/)

  nLineDay=0
  DO i=1,nlineInData
     IF (INT(IniCBLdata(i,1))<=id)THEN
        nLineDay=nLineDay+1
     ENDIF
  ENDDO
 
  IF(InitialData_use==2) THEN
     blh_m=IniCBLdata(nLineDay,2)
     gamt_Km=IniCBLdata(nLineDay,3)
     gamq_gkgm=IniCBLdata(nLineDay,4)
     tp_K=IniCBLdata(nLineDay,5)
     qp_gkg=IniCBLdata(nLineDay,6)
     tm_K=IniCBLdata(nLineDay,7)
     qm_gkg=IniCBLdata(nLineDay,8)
  ELSEIF(InitialData_use==1 .AND. IniCBLdata(nlineDay,1)==id)THEN   ! Changed from i to nlineDay, HCW 29 March 2017
     blh_m=IniCBLdata(nLineDay,2)
     gamt_Km=IniCBLdata(nLineDay,3)
     gamq_gkgm=IniCBLdata(nLineDay,4)
     tm_K_zm=(Temp_C+C2K)*((1000/Press_hPa)**(gas_ct_dry/avcp))
     tm_K=tm_K_zm-psyh*qh_use/(k*ustar*avcp*avdens)
     es_hPa=sat_vap_press(Temp_C,Press_hPa,1)
     qm_gkg_zm=622*avrh/(100*Press_hPa/es_hPa-avrh)
     lv=(2500.25-2.365*temp_C)*1000
     qm_gkg=qm_gkg_zm-psyh*qe_use/(k*ustar*avdens*lv)
     tp_K=tm_K
     qp_gkg=qm_gkg
  ELSEIF(InitialData_use==0)THEN
     blh_m=241.5
     gamt_Km=0.043
     gamq_gkgm=0.0092
     tm_K_zm=(Temp_C+C2K)*((1000/Press_hPa)**(gas_ct_dry/avcp))
     tm_K=tm_K_zm-psyh*qh_use/(k*ustar*avcp*avdens)
     es_hPa=sat_vap_press(Temp_C,Press_hPa,1)
     qm_gkg_zm=622*avrh/(100*Press_hPa/es_hPa-avrh)
     lv=(2500.25-2.365*temp_C)*1000
     qm_gkg=es_hPa-psyh*qe_use/(k*ustar*avdens*lv)
     tp_K=tm_K
     qp_gkg=qm_gkg
  ENDIF

  gamq_kgkgm=gamq_gkgm/1000.
  qp_kgkg=qp_gkg/1000    !humidities: g/kg -> kg/kg   q+
  qm_kgkg=qm_gkg/1000    !conc at mixing layer height h
  tp_C=tp_K-C2K
  tm_C=tm_K-C2K

  IF(sondeflag==1 .AND. cblday(id)==1) THEN
     !if gamma theta varies with z (selected by setting gthetaflag=1)
     !if gamma q varies with z (selected by setting ghumflag=1)
     CALL sonde(id)
     gamt_Km=0
     gamq_kgkgm=0
  ENDIF

  !adjusting qp and pm in case of saturation
  IF(qp_kgkg.GT.qsatf(tp_C,Press_hPa).OR.qp_kgkg.LT.0)THEN
     qp_kgkg = qsatf(tp_C,Press_hPa)
  ENDIF
  IF(qm_kgkg.GT.qsatf(tm_C,Press_hPa).OR.qm_kgkg.LT.0) THEN
     qm_kgkg = qsatf(tm_C,Press_hPa)
  ENDIF

  !    if((CBLuse==2).and.(zenith_deg>=90))then
  !    blh_m=188
  !    endif
  startflag=0


END SUBROUTINE CBL_initial

!------------------------------------------------------------------------
!------------------------------------------------------------------------
FUNCTION qsatf(T,PMB) RESULT(qsat)
  !       MRR, 1987
  ! AT TEMP T (DEG C) AND PRESSURE PMB (MB), GET SATURATION SPECIFIC
  !       HUMIDITY (KG/KG) FROM TETEN FORMULA
  !$$$$$$     use mod_parameter
  !$$$$$$     use mod_teten
  USE gas
  !$$$$$$     real (kind(1D0)):: T,es,qsat,pmb
  REAL (KIND(1D0))::T,es,qsat,PMB
  REAL (KIND(1D0))::A=6.106, B=17.27, C=237.3  !Teten coefficients

  IF(t.GT.55)THEN
     CALL ErrorHint(34,'Function qsatf',T,0.00D0,notUsedI)
  ENDIF

  ES = A*dEXP(B*T/(C+T))
  qsat = (molar_wat_vap/molar)*ES/PMB!(rmh2o/rmair)*ES/PMB
END FUNCTION qsatf
