! Note: INTERVAL is now set to 3600 s in Initial (it is no longer set in RunControl) HCW 29 Jan 2015

subroutine CBL(ifirst,iMB)

    use mod_z     
    use mod_k     
    use gas       
    use time      
    use data_in     
    use sues_data 
    use moist     
    use allocateArray   
    use defaultNotUsed
    use cbl_module
    use gis_data
    implicit none
   
    real(Kind(1d0))::sat_vap_press
    real(Kind(1d0))::qh_use,qe_use,tm_K_zm,qm_gkg_zm
    real(Kind(1d0))::Temp_C1,avrh1,es_hPa1   
    real(Kind(1d0))::secs0,secs1,Lv
    integer::idoy,ifirst,iMB,startflag
    real(Kind(1d0)), parameter::pi=3.141592653589793d+0,d2r=pi/180.

    
    !Skip first loop and unspecified days 
    if((ifirst==1 .and. iMB==1) .or. CBLday(id)==0) then   !HCW modified condition to check for first timestep of the model run
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,(NAN,is=6,22)/)                 
        return
    elseif(avkdn<5)then
        call CBL_initial(qh_use,qe_use,tm_K_zm,qm_gkg_zm,startflag,iMb)
        return
    endif

    if(startflag==0)then !write down initial values in previous time step
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,blh_m,tm_K,qm_kgkg*1000,&
        tp_K,qp_kgkg*1000,(NAN,is=11,20),gamt_Km,gamq_kgkgm/)
        startflag=1
    endif
    
    !Heat flux choices           
    if(Qh_choice==1) then   !from SUEWS
        qh_use=qh
        qe_use=qeph
    elseif(qh_choice==2)then !from LUMPS
        qh_use=H_mod
        qe_use=E_mod
    elseif(qh_choice==3)then  !from OBS
        if(qh_obs<-900.or.qe_obs<-900)then  ! observed data has a problem
            call ErrorHint(22,'Problem in observed Qh/Qe_value',qh_obs,qe_obs,qh_choice)
        endif
        qh_use=qh_obs
        qe_use=qe_obs
    endif

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
    if(CO2_included==1)then
        fcbl = 0!fc(i)/(rmco2/volm)      ! units: mol/mol * m/s
    else
        cm=NAN  
    endif

    tpp_K=tp_K
    qpp_kgkg=qp_kgkg                                    

    if(sondeflag.eq.1) then
            call gamma_sonde       
    endif
!     	set up array for Runge-Kutta call
    blh1_m=blh_m
    y(1)=blh_m ! integrate h, t, q, c from time s(i-1)
    y(2)=tm_K  ! to time s(i) using runge-kutta solution
    y(3)=qm_kgkg   ! of slab CBL equations
    y(4)=cm

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call rkutta(neqn,secs0,secs1,y,1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    blh_m   =y(1)
    tm_K    =y(2)  ! potential temperature, units: deg C 
    qm_kgkg =y(3)  ! specific humidity, units: kg/kg
    cm      =y(4)  ! co2 concentration,units: mol/mol   

!     	compute derived quantities for this time step

    tp_K   = tpp_K     + (gamt_Km*(blh_m-blh1_m))
    qp_kgkg = qpp_kgkg + (gamq_kgkgm*(blh_m-blh1_m)) 

    if (tp_K.lt.tm_K) then
            tp_K = tm_K
    endif

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
    if(it==0 .and. imin==(nsh_real-1)/nsh_real*60) idoy=id-1  !Modified by HCW 04 Mar 2015 in case model timestep is not 5-min

    
    if((qh_choice==1).or.(qh_choice==2))then !BLUEWS or BLUMPS
        !Stability correction
        !tm_K_zm=tm_K+cbldata(10)*cbldata(2)/(k*cbldata(8)*cbldata(6)*cbldata(4))
        Temp_C=tm_K/((1000/cbldata(9))**(gas_ct_dry/cbldata(6)))-C2K 
        es_hPa=sat_vap_press(Temp_C,cbldata(9),1)           
        lv=(2500.25-2.365*Temp_C)*1000
        !qm_gkg_zm=qm_gkg+cbldata(10)*cbldata(3)/(k*cbldata(8)*cbldata(4)*lv)
        avrh=100*((qm_gkg*cbldata(9)/(622+qm_gkg))/es_hPa) !check pressure
        if(avrh>100)then
            call errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbldata(1)/24.0,avrh,100)
            avrh=100     
        endif
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,blh_m,tm_K,qm_kgkg*1000,&
        tp_K,qp_kgkg*1000,&
        Temp_C,avrh,cbldata(2),cbldata(3),cbldata(9),cbldata(7),cbldata(8),cbldata(4),cbldata(5),cbldata(6),&
        gamt_Km,gamq_kgkgm/) 
    elseif(qh_choice==3)then ! CBL
        !tm_K_zm=tm_K+cbldata(10)*cbldata(2)/(k*cbldata(8)*cbldata(6)*cbldata(4))
        Temp_C1=tm_K/((1000/cbldata(9))**(gas_ct_dry/cbldata(6)))-C2K   
        es_hPa1=sat_vap_press(Temp_C1,cbldata(9),1)
        lv=(2500.25-2.365*Temp_C1)*1000
        !qm_gkg_zm=qm_gkg+cbldata(10)*cbldata(3)/(k*cbldata(8)*cbldata(4)*lv)
        avrh1=100*((qm_gkg*cbldata(8)/(622+qm_gkg))/es_hPa1) !check pressure
        if(avrh1>100)then
            call errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbldata(1)/24.0,avrh,100)
            avrh1=100     
        endif
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,blh_m,tm_K,qm_kgkg*1000,&
        tp_K,qp_kgkg*1000,&
        Temp_C1,avrh1,cbldata(2),cbldata(3),cbldata(9),cbldata(7),cbldata(8),cbldata(4),cbldata(5),cbldata(6),&
        gamt_Km,gamq_kgkgm/)
    endif

    return 

end subroutine CBL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine CBL_ReadInputData
	use allocateArray
	use data_in
	use sues_data
	use cbl_module
        use initial
	implicit none
        integer::i
        real(kind(1d0))::l

	namelist/CBLInput/EntrainmentType,&
                            QH_choice,&
                            CO2_included,&
                            cblday,&
                            wsb,&
                            InitialData_use,&
                            InitialDataFileName,&
                            sondeflag,&
                            FileSonde

	open(51,file=trim(FileInputPath)//'CBLInput.nml',status='old', err=24) 
        read(51,nml=CBLInput,err=24)
	close(51)

        !Read initial values if it's needed
	if(InitialData_use==1 .or. InitialData_use==2)then
            open(52,file=trim(FileInputPath)//trim(InitialDataFileName),status='old', err=25)
            read(52,*)
            nlineInData = 0   !Initialise nlines
            do
                read(52,*) l
                if (l == -9) exit
                nlineInData = nlineInData + 1
            enddo            
            close(52)

            allocate(IniCBLdata(1:nlineInData,1:8))
            open(52,file=trim(FileInputPath)//trim(InitialDataFileName),status='old', err=25)
            read(52,*)
            do i=1,nlineInData
              read(52,*)IniCBLdata(i,1:8)
            enddo
            close(52)
	endif

	if(CO2_included==0)then
		fcbl=0      ! hard-wire no CO2
	endif

        iCBLcount=0

	return

24  call ErrorHint(24,'CBLInput.nml',0.00D0,0.000D0,0) 
25  call ErrorHint(24,trim(FileInputPath)//trim(InitialDataFileName),0.00D0,0.00D0,0)
   
end subroutine CBL_ReadInputData

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine CBL_initial(qh_use,qe_use,tm_K_zm,qm_gkg_zm,startflag,iMB)
  
    use mod_z     
    use mod_k     
    use gas       
    use time      
    use data_in     
    use sues_data 
    use moist     
    use allocateArray   
    use defaultNotUsed
    use cbl_module
    use gis_data
    implicit none
    
    real(Kind(1d0))::qh_use,qe_use,tm_K_zm,qm_gkg_zm
    real(Kind(1d0))::qsatf,sat_vap_press,lv
    integer::i,nLineDay,iMB,startflag
    
    !Heat flux choices           
    if(Qh_choice==1) then   !from SUEWS
        qh_use=qh
        qe_use=qeph
    elseif(qh_choice==2)then !from LUMPS
        qh_use=H_mod
        qe_use=E_mod
    elseif(qh_choice==3)then  !from OBS
        if(qh_obs<-900.or.qe_obs<-900)then  ! observed data has a problem
            call ErrorHint(22,'Problem in observed Qh/Qe_value',qh_obs,qe_obs,qh_choice)
        endif
        qh_use=qh_obs
        qe_use=qe_obs
    endif
    
    
    blh_m=NAN
    iCBLcount=iCBLcount+1
    dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,(NAN,is=6,22)/)   
    
    nLineDay=0
    do i=1,nlineInData
        if (int(IniCBLdata(i,1))<=id)then
        nLineDay=nLineDay+1
        endif
    enddo
    
    if(InitialData_use==2) then
        blh_m=IniCBLdata(nLineDay,2)
        gamt_Km=IniCBLdata(nLineDay,3)
        gamq_gkgm=IniCBLdata(nLineDay,4)
        tp_K=IniCBLdata(nLineDay,5)
        qp_gkg=IniCBLdata(nLineDay,6)
        tm_K=IniCBLdata(nLineDay,7)
        qm_gkg=IniCBLdata(nLineDay,8)
    elseif(InitialData_use==1 .and. IniCBLdata(i,1)==id)then
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
    elseif(InitialData_use==0)then
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
    endif

    gamq_kgkgm=gamq_gkgm/1000.
    qp_kgkg=qp_gkg/1000    !humidities: g/kg -> kg/kg   q+
    qm_kgkg=qm_gkg/1000    !conc at mixing layer height h
    tp_C=tp_K-C2K
    tm_C=tm_K-C2K
    
    if(sondeflag==1 .and. cblday(id)==1) then 
        !if gamma theta varies with z (selected by setting gthetaflag=1)
        !if gamma q varies with z (selected by setting ghumflag=1)      
        call sonde(id)
        gamt_Km=0
        gamq_kgkgm=0
    endif

    !adjusting qp and pm in case of saturation         
    if(qp_kgkg.gt.qsatf(tp_C,Press_hPa).or.qp_kgkg.lt.0)then  
        qp_kgkg = qsatf(tp_C,Press_hPa)
    endif
    if(qm_kgkg.gt.qsatf(tm_C,Press_hPa).or.qm_kgkg.lt.0) then 
        qm_kgkg = qsatf(tm_C,Press_hPa)      
    endif
    
!    if((CBLuse==2).and.(zenith_deg>=90))then
!    blh_m=188
!    endif
    startflag=0

    
end subroutine CBL_initial

!------------------------------------------------------------------------
!------------------------------------------------------------------------
FUNCTION qsatf(T,PMB) result(qsat)
!       MRR, 1987
! AT TEMP T (DEG C) AND PRESSURE PMB (MB), GET SATURATION SPECIFIC
!       HUMIDITY (KG/KG) FROM TETEN FORMULA
!$$$$$$     use mod_parameter
!$$$$$$     use mod_teten
    use gas
!$$$$$$     real (kind(1D0)):: T,es,qsat,pmb
    real (kind(1D0))::T,es,qsat,PMB
    real (kind(1D0))::A=6.106, B=17.27, C=237.3  !Teten coefficients

	if(t.gt.55)then
       call ErrorHint(34,'Function qsatf',T,0.00D0,notUsedI)	   
	endif

	ES = A*dEXP(B*T/(C+T))
	qsat = (molar_wat_vap/molar)*ES/PMB!(rmh2o/rmair)*ES/PMB
END function qsatf
