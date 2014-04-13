subroutine CBL
	! LUMPS_Module_constants_v6_0.f95
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
	IMPLICIT None
   
	real(Kind(1d0))::sat_vap_press,qsatf,deltv,qh_use,qe_use,tm_K_zm
	real(Kind(1d0))::sol_decl,Temp_C1,avrh1,es_hPa1   
	real(Kind(1d0))::secs0,secs1
	character(len=6)::iday
	integer::SunriseTime,SunsetTime,i,j
	real, parameter::pi=3.141592653589793d+0,d2r=pi/180.

	if(CBLday(id)==0) then   
		blh_m=NAN
		return
	endif
   
	!Calculate Sunrisetime & Sunsettime
	if(CBLday(id)==1)then
		if((mod(year,100.).ne.0.).and.(mod(year,4.).eq.0.))then  !leap year
			sol_decl=0.409*cos(2.*pi*float(id-174)/float(366))     
		elseif(mod(year,400.).eq.0.)then  !leap year
			sol_decl=0.409*cos(2.*pi*float(id-174)/float(366))
		else !not leap year 
			sol_decl=0.409*cos(2.*pi*float(id-173)/float(365)) 
		endif
		if(lng.ge.0)then
     	SunriseTime=24./(2.*pi)*(lng*d2r+acos((sin(lat*d2r)*sin(sol_decl)-sin(0.))/(cos(lat*d2r)*cos(sol_decl))))+Timezone+1 
     	SunsetTime =24./(2.*pi)*(lng*d2r-acos((sin(lat*d2r)*sin(sol_decl)-sin(0.))/(cos(lat*d2r)*cos(sol_decl))))+Timezone+24  
		else
        SunriseTime=24./(2.*pi)*(-lng*d2r+acos((sin(lat*d2r)*sin(sol_decl)-sin(0.))/(cos(lat*d2r)*cos(sol_decl))))+Timezone+1
     	SunsetTime =24./(2.*pi)*(-lng*d2r-acos((sin(lat*d2r)*sin(sol_decl)-sin(0.))/(cos(lat*d2r)*cos(sol_decl))))+Timezone+24  
     	endif
	endif

	!Heat flux choices                  
	if(Qh_choice==1) then   !SUEWS
      	qh_use=qh
      	qe_use=qeph
	elseif(qh_choice==2)then ! LUMPS
      	qh_use=H_mod
      	qe_use=E_mod
	elseif(qh_choice==3)then  !obs
		if(qh_obs<-900.or.qe_obs<-900)then  ! observed data has a problem
        	call ErrorHint(22,'Qh_choice_from_CBLinital',qh_obs,qe_obs,qh_choice)
		endif
		qh_use=qh_obs
		qe_use=qe_obs
	endif  
   
	!Take the initial values at a time step before starting CBL model and return
	if((it.le.SunriseTime).or.(it.gt.SunsetTime)) then ! To do -- need to put another alternative in  
    	start1=0
        blh_m=NAN
		if(interval==60)then
			cbldata(0,0)=it
			cbldata(0,1)=qh_use
			cbldata(0,2)=qe_use
			cbldata(0,3)=avdens
			cbldata(0,4)=lv_J_kg
			cbldata(0,5)=avcp
			cbldata(0,6)=avu1
			cbldata(0,7)=ustar
			cbldata(0,8)=Press_hPa    
		endif
		return
	elseif((qh_use<0).and.(start1==0)) then!If qh_use is still below zero, CBL model is not used and skipped to next loop.
		blh_m=NAN
		if(interval==60)then
        	cbldata(0,0)=it
        	cbldata(0,1)=qh_use
        	cbldata(0,2)=qe_use
        	cbldata(0,3)=avdens
        	cbldata(0,4)=lv_J_kg
        	cbldata(0,5)=avcp
        	cbldata(0,6)=avu1
        	cbldata(0,7)=ustar
        	cbldata(0,8)=Press_hPa
		endif
		return  
	endif

	!! Start CBL model
	if(start1.eq.0)then         
		!m-  potential temperature and specific humidity in the well-mixed portion of the CBL
		!p (+)  potential temperature and specific humidity immediately above zi
        
        !!choose how to get initial data
        if(InitialData_use==2) then!read intial data from a file
			read(52,*)which_day,blh_m,gamt_Km,gamq_gkgm,tp_K,qp_gkg,tm_K,qm_gkg
        elseif(InitialData_use==1)then
        	read(52,*)which_day,blh_m,gamt_Km,gamq_gkgm!,tp_K,qp_gkg,tm_K,qm_gkg
            tm_K_zm=(Temp_C+C2K)*((1000/Press_hPa)**(gas_ct_dry/avcp))
            tm_K=tm_K_zm-psyh*qh_use/(k*ustar*avcp*avdens)             
			es_hPa=sat_vap_press(Temp_C,Press_hPa,1)
			!avrh=100*((qm_gkg*Press_hPa/(622+qm_gkg))/es_hPa)
            qm_gkg=622*avrh/(100*Press_hPa/es_hPa-avrh)
            !!this needs to be check
            tp_K=tm_K
            qp_gkg=qm_gkg
        endif
		write(12,'(i3,f5.0,6g14.6)')which_day,blh_m,gamt_Km,gamq_gkgm,tp_K,qp_gkg,tm_K,qm_gkg

		gamq_kgkgm=gamq_gkgm/1000.
       
		if(sondeflag.eq.1) then 
        	!if gamma theta varies with z (selected by setting gthetaflag=1)
        	!if gamma q varies with z (selected by setting ghumflag=1)      
        	call sonde
			gamt_Km=0
			gamq_kgkgm=0
     	endif
               
		qp_kgkg=qp_gkg/1000    !humidities: g/kg -> kg/kg   q+
		qm_kgkg=qm_gkg/1000    ! conc at mixing layer height h
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
		open(53,file=trim(FileOutputPath)//'CBL_'//trim(adjustl(iday))//'.txt',status='unknown')
        WRITE(53,*)'which_day blh_m gamt_Km gamq_gkgm tp_K qp_gkg tm_K qm_gkg'
		write(53,'(i3,f5.0,6g14.6)')which_day,blh_m,gamt_Km,gamq_kgkgm,tp_K,qp_gkg,tm_K,qm_gkg 
		write(53, 102)
102  	format('id    it   dectime    z         theta       q       theta+      q+       Temp_C',&
               '    rh        QH_use    QE_use  Press_hPa    avu1      ustar     avdens',&
               '    lv_J_kg avcp         gamt       gamq') 
     	start1=start1+1
!$$$$$$      avblh=blh_m
!$$$$$$      orig_gamt_Km=gamt_Km
!$$$$$$      orig_gamq_kgkgm=gamq_kgkgm
	endif
      
!$$$$$$   if(CBLday(id)==1)then
!$$$$$$     if((it.lt.(SunriseTime+0.01)).or.(it.ge.SunsetTime+2)) then
!$$$$$$       blh_m=avblh
!$$$$$$     endif
!$$$$$$     if(it==(SunsetTime+1))then
!$$$$$$       blh_m=(blh_m-avblh)/2+avblh
!$$$$$$     endif
!$$$$$$   endif    

  
!$$$$$$   if((it==SunsetTime))then
!$$$$$$     gamt_Km=-gamt_Km
!$$$$$$   endif
!$$$$$$   if((it==SunriseTime))then
!$$$$$$     gamt_Km=-gamt_Km
!$$$$$$   endif
!$$$$$$         if(it>=(SunsetTime+1).or.it<(Sunrisetime))then     
!$$$$$$         gamt_Km=-orig_gamt_Km*100
!$$$$$$         gamq_kgkgm=-orig_gamq_kgkgm*100
!$$$$$$         write(*,*)it,gamt_Km
!$$$$$$         !pause
!$$$$$$         else
!$$$$$$         gamt_Km=orig_gamt_Km
!$$$$$$         gamq_kgkgm=orig_gamq_kgkgm
!$$$$$$         endif

		cbldata(1,0)=it
		cbldata(1,1)=qh_use
        cbldata(1,2)=qe_use
        cbldata(1,3)=avdens
        cbldata(1,4)=lv_J_kg
        cbldata(1,5)=avcp
        cbldata(1,6)=avu1
        cbldata(1,7)=ustar
        cbldata(1,8)=Press_hPa

	do i=0,nCblStep-1

		cbld(0)=(float(i*tstep_s)/3600.)+cbldata(0,0)

		secs0=cbld(0)*3600.
		secs1=secs0+float(tstep_s) ! time in seconds
      
		do j=1,8
			cbld(j)=(i*(cbldata(1,j)-cbldata(0,j))/nCBLstep)+cbldata(0,j)
		enddo 
           
		! Kinematic fluxes   
		fhbl_Kms    = cbld(1)/ (cbld(3)*cbld(5))  !qh_use/(avdens*avcp)      ! units: degK * m/s 
		febl_kgkgms = cbld(2)/ (cbld(3)*cbld(4))  !qe_use/(avdens*lv_J_kg)   ! units: kg/kg * m/s 

		if(CO2_included==1)then
			fcbl = 0!fc(i)/(rmco2/volm)      ! units: mol/mol * m/s
		else
			cm=NAN  ! s.o.
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
		 
		deltv = (tp_K - tm_K) + 0.61*tm_k*(qp_kgkg - qm_kgkg)  ! pot virtual temp
      
		qm_gkg=qm_kgkg*1000 !humidities: kg/kg -> g/kg

		if((qh_choice==1).or.(qh_choice==2))then !BLUEWS or BLUMPS
         	!Stability correction
			tm_K_zm=tm_K+psyh*cbld(1)/(k*cbld(7)*cbld(5)*cbld(3))
			Temp_C=tm_K_zm/((1000/cbld(8))**(gas_ct_dry/cbld(5)))-C2K
			es_hPa=sat_vap_press(Temp_C,cbld(8),1)
			avrh=100*((qm_gkg*cbld(8)/(622+qm_gkg))/es_hPa) !check pressure
			if(avrh>100)then
				call errorHint(34,'subroutine CBL dectime, relative humidity',id+cbld(0)/24.0,avrh,100)
         		avrh=100     
			endif   
	        write(53,'(i4,f6.2,f9.4,13(f10.3,1x),2(f15.7,1X),2(g15.7,1x))')id,cbld(0),id+cbld(0)/24.,blh_m,tm_K,&
            qm_kgkg*1000,tp_K,qp_kgkg*1000,Temp_C,avrh,cbld(1),cbld(2),cbld(8),cbld(6),cbld(7),cbld(3),&
            cbld(4),cbld(5),gamt_Km,gamq_kgkgm

		else ! CBL
         	tm_K_zm=tm_K+psyh*cbld(1)/(k*cbld(7)*cbld(5)*cbld(3))
         	Temp_C1=tm_K_zm/((1000/cbld(8))**(gas_ct_dry/cbld(5)))-C2K
      	 	es_hPa1=sat_vap_press(Temp_C1,cbld(8),1)
      	 	avrh1=100*((qm_gkg*cbld(8)/(622+qm_gkg))/es_hPa1) !check pressure
			if(avrh1>100)then
         		call errorHint(34,'subroutine CBL dectime, relative humidity',id+cbld(0)/24.0,avrh,100)
         		avrh1=100     
			endif     
	        write(53,'(i4,f6.2,f9.4,13(f10.3,1x),2(f15.7,1X),2(g15.7,1x))')id,cbld(0),id+cbld(0)/24.,blh_m,tm_K,&
            qm_kgkg*1000,tp_K,qp_kgkg*1000,Temp_C,avrh,cbld(1),cbld(2),cbld(8),cbld(6),cbld(7),cbld(3),&
            cbld(4),cbld(5),gamt_Km,gamq_kgkgm
		endif            
	enddo

	cbldata(0,0)=it
	cbldata(0,1)=qh_use
	cbldata(0,2)=qe_use
	cbldata(0,3)=avdens
	cbldata(0,4)=lv_J_kg
 	cbldata(0,5)=avcp
	cbldata(0,6)=avu1
	cbldata(0,7)=ustar
	cbldata(0,8)=Press_hPa
        
	return 

end subroutine CBL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Subroutine CBL_initial
	use allocateArray
	use data_in
	use cbl_module
	implicit none

	namelist/CBLInput/EntrainmentType,wsb,QH_choice,sondeflag,CO2_included,&
   		FileSonde,InitialDataFileName,tstep_s,cblday,InitialData_use
    
	open(51,file=trim(FileInputPath)//'CBLInput.nml',status='old', err=24) 
	read(51,nml=CBLInput,err=24)
	close(51)
	nCBLstep=(Interval*60)/tstep_s
	write(12,*)'--------------CBLml---------------------------'
	write(12,nml=CBLInput)
	write(12,*)'day blh_m   gamt_Km     gamq_gkgm     tp_K          qp_gkg        tm_K         qm_gkg'

	if(InitialData_use==1 .or. InitialData_use==2)then
   		open(52,file=trim(FileInputPath)//trim(InitialDataFileName),status='old', err=25)
   		read(52,*)
	endif
   
	if(CO2_included==0)then
		fcbl=0      ! hard-wire no CO2
	endif

	start1=0
 
	return

24  	call ErrorHint(24,'CBLInput.nml',0.00D0,0.000D0,0) 
25 		call ErrorHint(24,trim(FileInputPath)//trim(InitialDataFileName),0.00D0,0.00D0,0)
   
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
	   print*,t, " Function qsatf: Temperature doesn't appear to be in deg C"
	   pause
	endif

	ES = A*dEXP(B*T/(C+T))
	qsat = (molar_wat_vap/molar)*ES/PMB!(rmh2o/rmair)*ES/PMB
END function qsatf
