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
    real(Kind(1d0))::qsatf,qh_use,qe_use,tm_K_zm,qm_gkg_zm
    real(Kind(1d0))::Temp_C1,avrh1,es_hPa1   
    real(Kind(1d0))::secs0,secs1,Lv
    character(len=6)::iday
    integer::i,j,idoy,ifirst,ncol,iMB
    real(Kind(1d0)), parameter::pi=3.141592653589793d+0,d2r=pi/180.

    !Skip CBL for first loop
    if(ifirst==1 .or. CBLday(id)==0)then
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,(NAN,is=6,22)/)                 
        return
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

    if(CBLuse==1)then	! CBL model is used only daytime 
        if(zenith_deg>=85) then 
            start1=0
            blh_m=NAN
            cbldata(1,1)=float(it)+float(imin)/60.
            cbldata(1,2)=qh_use
            cbldata(1,3)=qe_use
            cbldata(1,4)=avdens
            cbldata(1,5)=lv_J_kg
            cbldata(1,6)=avcp
            cbldata(1,7)=avu1
            cbldata(1,8)=ustar
            cbldata(1,9)=Press_hPa
            cbldata(1,10)=psyh 
            iCBLcount=iCBLcount+1
            dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,(NAN,is=6,22)/)             
            return
        endif      
    elseif((CBLuse==2).and.(start2==0)) then! CBL model is used fulltime but not working yet
        start2=1
        blh_m=NAN
        cbldata(1,1)=float(it)+float(imin)/60.
        cbldata(1,2)=qh_use
        cbldata(1,3)=qe_use
        cbldata(1,4)=avdens
        cbldata(1,5)=lv_J_kg
        cbldata(1,6)=avcp
        cbldata(1,7)=avu1
        cbldata(1,8)=ustar
        cbldata(1,9)=Press_hPa
        cbldata(1,10)=psyh  
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,(NAN,is=6,22)/) 
        return    
    endif  

    if(start1.eq.0)then         
        !m-  potential temperature and specific humidity in the well-mixed portion of the CBL
        !p (+)  potential temperature and specific humidity immediately above zi
        jday=jday+1
        if(InitialData_use==2) then
            which_day=int(IniCBLdata(jday,1))
            blh_m=IniCBLdata(jday,2)
            gamt_Km=IniCBLdata(jday,3)
            gamq_gkgm=IniCBLdata(jday,4)
            tp_K=IniCBLdata(jday,5)
            qp_gkg=IniCBLdata(jday,6)
            tm_K=IniCBLdata(jday,7)
            qm_gkg=IniCBLdata(jday,8)
        elseif(InitialData_use==1)then
            which_day=int(IniCBLdata(jday,1))
            blh_m=IniCBLdata(jday,2)
            gamt_Km=IniCBLdata(jday,3)
            gamq_gkgm=IniCBLdata(jday,4)
            tm_K_zm=(Temp_C+C2K)*((1000/Press_hPa)**(gas_ct_dry/avcp))
            tm_K=tm_K_zm-psyh*qh_use/(k*ustar*avcp*avdens)             
            es_hPa=sat_vap_press(Temp_C,Press_hPa,1)
            qm_gkg_zm=622*avrh/(100*Press_hPa/es_hPa-avrh)
            lv=(2500.25-2.365*temp_C)*1000
            qm_gkg=qm_gkg_zm-psyh*qe_use/(k*ustar*avdens*lv) 
            tp_K=tm_K
            qp_gkg=qm_gkg
        elseif(InitialData_use==0)then
            which_day=id
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

        if(sondeflag.eq.1) then 
            !if gamma theta varies with z (selected by setting gthetaflag=1)
            !if gamma q varies with z (selected by setting ghumflag=1)      
            call sonde
            gamt_Km=0
            gamq_kgkgm=0
        endif

        qp_kgkg=qp_gkg/1000    !humidities: g/kg -> kg/kg   q+
        qm_kgkg=qm_gkg/1000    !conc at mixing layer height h
        tp_C=tp_K-C2K
        tm_C=tm_K-C2K

        !adjusting qp and pm in case of saturation         
        if(qp_kgkg.gt.qsatf(tp_C,Press_hPa).or.qp_kgkg.lt.0)then  
            qp_kgkg = qsatf(tp_C,Press_hPa)
        endif
        if(qm_kgkg.gt.qsatf(tm_C,Press_hPa).or.qm_kgkg.lt.0) then 
            qm_kgkg = qsatf(tm_C,Press_hPa)      
        endif

        write(iday,'(i3)')which_day

        start1=1
    endif
      
    if((CBLuse==2).and.(zenith_deg>=90))then
      blh_m=188
    endif

    idoy=id
    If(it==0 .and. imin==55) idoy=id-1     

    cbldata(2,1)=float(it)+float(imin)/60.
    cbldata(2,2)=qh_use
    cbldata(2,3)=qe_use
    cbldata(2,4)=avdens
    cbldata(2,5)=lv_J_kg
    cbldata(2,6)=avcp
    cbldata(2,7)=avu1
    cbldata(2,8)=ustar
    cbldata(2,9)=Press_hPa
    cbldata(2,10)=psyh

    do i=1,nCBLstep

        !Unit of calculation time is decimal hour
        cbld(1)=float((i-1)*tstep_s)/3600.+cbldata(1,1)

        secs0=cbld(1)*3600.
        secs1=secs0+float(tstep_s) ! time in seconds

        do j=2,10
            cbld(j)=(i*(cbldata(2,j)-cbldata(1,j))/float(nCBLstep))+cbldata(1,j)
        enddo 

        ! Kinematic fluxes   
        fhbl_Kms    = cbld(2)/ (cbld(4)*cbld(6))  !qh_use/(avdens*avcp)      ! units: degK * m/s 
        febl_kgkgms = cbld(3)/ (cbld(4)*cbld(5))  !qe_use/(avdens*lv_J_kg)   ! units: kg/kg * m/s 

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
        call rkutta(neqn,secs0,secs1,y,10)
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

!		th = tm_K - (grav/cbld(5))*blh_m			 ! actual temp just below z=h
!		dh = qsatf(th,cbld(8)) - qm_kgkg           ! deficit just below z=h

        tp_C=tp_K-C2K      
        tm_C=tm_K-C2K

        !	 delt = tp_K - tm_K ! temp
        !	 delq = qp_kgkg - qm_kgkg ! humidity
        !deltv = (tp_K - tm_K) + 0.61*tm_k*(qp_kgkg - qm_kgkg)  ! pot virtual temp

        qm_gkg=qm_kgkg*1000 !humidities: kg/kg -> g/kg

    enddo

    if((qh_choice==1).or.(qh_choice==2))then !BLUEWS or BLUMPS
        !Stability correction
        !tm_K_zm=tm_K+cbld(10)*cbld(2)/(k*cbld(8)*cbld(6)*cbld(4))
        Temp_C=tm_K/((1000/cbld(9))**(gas_ct_dry/cbld(6)))-C2K 
        es_hPa=sat_vap_press(Temp_C,cbld(9),1)           
        lv=(2500.25-2.365*Temp_C)*1000
        !qm_gkg_zm=qm_gkg+cbld(10)*cbld(3)/(k*cbld(8)*cbld(4)*lv)
        avrh=100*((qm_gkg*cbld(9)/(622+qm_gkg))/es_hPa) !check pressure
        if(avrh>100)then
            call errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbld(1)/24.0,avrh,100)
            avrh=100     
        endif
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,blh_m,tm_K,qm_kgkg*1000,&
        tp_K,qp_kgkg*1000,&
        Temp_C,avrh,cbld(2),cbld(3),cbld(9),cbld(7),cbld(8),cbld(4),cbld(5),cbld(6),&
        gamt_Km,gamq_kgkgm/) 
    else ! CBL
        !tm_K_zm=tm_K+cbld(10)*cbld(2)/(k*cbld(8)*cbld(6)*cbld(4))
        Temp_C1=tm_K/((1000/cbld(9))**(gas_ct_dry/cbld(6)))-C2K   
        es_hPa1=sat_vap_press(Temp_C1,cbld(9),1)
        lv=(2500.25-2.365*Temp_C1)*1000
        !qm_gkg_zm=qm_gkg+cbld(10)*cbld(3)/(k*cbld(8)*cbld(4)*lv)
        avrh1=100*((qm_gkg*cbld(8)/(622+qm_gkg))/es_hPa1) !check pressure
        if(avrh1>100)then
            call errorHint(34,'subroutine CBL dectime, relative humidity',idoy+cbld(1)/24.0,avrh,100)
            avrh1=100     
        endif
        iCBLcount=iCBLcount+1
        dataOutBL(iCBLcount,1:22,iMB)=(/real(iy,8),real(id,8),real(it,8),real(imin,8),dectime,blh_m,tm_K,qm_kgkg*1000,&
        tp_K,qp_kgkg*1000,&
        Temp_C1,avrh1,cbld(2),cbld(3),cbld(9),cbld(7),cbld(8),cbld(4),cbld(5),cbld(6),&
        gamt_Km,gamq_kgkgm/)
    endif

    cbldata(1,1)=float(it)+float(imin)/60.
    cbldata(1,2)=qh_use
    cbldata(1,3)=qe_use
    cbldata(1,4)=avdens
    cbldata(1,5)=lv_J_kg
    cbldata(1,6)=avcp
    cbldata(1,7)=avu1
    cbldata(1,8)=ustar
    cbldata(1,9)=Press_hPa
    cbldata(1,10)=psyh

    return 

end subroutine CBL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine CBL_initial
	use allocateArray
	use data_in
	use sues_data
	use cbl_module
        use initial
	implicit none
        integer::i,j,nlineInData
        real(kind(1d0))::l
        real(kind(1d0)),dimension(nlines,8)::indata

	namelist/CBLInput/EntrainmentType,&
                            QH_choice,&
                            CO2_included,&
                            cblday,&
                            wsb,&
                            tstep_s,&
                            InitialData_use,&
                            InitialDataFileName,&
                            sondeflag,&
                            FileSonde
    
	open(51,file=trim(FileInputPath)//'CBLInput.nml',status='old', err=24) 
        read(51,nml=CBLInput,err=24)
	close(51)
	nCBLstep=(t_Interval)/tstep_s

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

	start1=0
        start2=0
        jday=0
        iCBLcount=0

	return

24  	call ErrorHint(24,'CBLInput.nml',0.00D0,0.000D0,0) 
25 	call ErrorHint(24,trim(FileInputPath)//trim(InitialDataFileName),0.00D0,0.00D0,0)
   
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
