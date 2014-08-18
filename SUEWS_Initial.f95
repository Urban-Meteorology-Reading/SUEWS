! sg feb 2012
! only run once at start - fixed for all grids and all years
 subroutine OverallRunControl
  use resist        ! LUMPS_metRead.f95
  use data_in      ! LUMPS_metRead.f95
  use ohm_calc
  use run_info     ! run_control_v2_1.f90
  use SUES_data    !
  USE GIS_data     ! LUMPS_gis_read.f95
  use mod_z	       ! module_LUMPS_constants,f90
  use FileName
  use allocateArray  ! module_LUMPS_constants,f90
  use defaultNotUsed
  use time
  use snowMod
  
  IMPLICIT NONE
  
  real (Kind(1d0)):: skip, NSH_real
  integer::iv,i !sh (1=NH)
  character (len=50)::ProbLine
  
  namelist/RunControl/AnthropHeatChoice,&
        CBLuse,& !s.o.
        GISInputType,&
        NetRadiationChoice,&
        RoughLen_heat,&
        QSChoice,&
        smd_choice,&							
        StabilityMethod,&
        write5min,&
        writedailystate,&	
        WU_choice,&
        z0_method,&  
        defaultQf,&   
        defaultQs,&
        FileCode,&
        FileInputPath,&
        FileOutputPath,&
        SkipHeaderGIS,&
        SkipHeaderMet,&
        SnowFractionChoice,&
        SNOWuse,&
        SOLWEIGout,&
        Interval,&          
        TIMEZONE,& 
        Tstep,& 
        Z
        
  FileCode='none'
  !smithFile='Smith1966.grd'
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Read in the runcontrol.nml file
  open(55,File='RunControl.nml',err=200,status='old') !Change with needs
  read(55,nml=RunControl,err=201)
  close(55)

  !Problems with filecode?
  if (FileCode=='none') call ErrorHint(26,trim("RunControl.nml FileCode is missing"),notUsed,notUsed,notUsedI)
  
  !--------------------------------------------------------------------------
 
  !Read in FileChoices.txt
  FileChoices=trim(FileOutputPath)//trim(FileCode)//'FileChoices.txt'
  open(12,file=FileChoices,err=203)
  write(12,nml=RunControl)

  !Determines what should be done with respect to radiation
  AlbedoChoice=0
  NARPOutput=0
  ldown_option=0
  if(netRadiationChoice==0)then     !Observation from the input file is used
    
    if(snowUse==1) then             !If snow is modelled, NARP is needed for surface temperature
       netRadiationChoice=3000
       ldown_option=3   !Ldown will be modeled
       !NetRadiationChoice=NetRadiationChoice/1000
    endif
    
  elseif(netRadiationChoice>0)then  !Modelled is used
       AlbedoChoice=-9
       NARPOutput=-9
       if(NetRadiationChoice<10) then
            AlbedoChoice=0
            NARPOutput=0
            if(NetRadiationChoice==1)ldown_option=1
            if(NetRadiationChoice==2)ldown_option=2
            if(NetRadiationChoice==3)ldown_option=3
              
       elseif(NetRadiationChoice>9.and.NetRadiationChoice<100) then  
            AlbedoChoice=0  
            NARPOutput=1
            if(NetRadiationChoice==10)ldown_option=1
            if(NetRadiationChoice==20)ldown_option=2
            if(NetRadiationChoice==30)ldown_option=3
            NetRadiationChoice=NetRadiationChoice/10
            
       elseif(NetRadiationChoice>=100.and.NetRadiationChoice<1000) then
            AlbedoChoice=1  
            NARPOutput=1
            if(NetRadiationChoice==100)ldown_option=1
            if(NetRadiationChoice==200)ldown_option=2
            if(NetRadiationChoice==300)ldown_option=3
           	NetRadiationChoice=NetRadiationChoice/100
       endif
       
       !If bad netRadiationChoice value
	   if(netRadiationChoice>3.or. AlbedoChoice==-9.or.NARPOutput==-9)then
           write(*,*) 'NetRadiationChoice=',  NetRadiationChoice
           write(*,*) 'Value not usable'
           stop
       endif
  endif   
  ! %%%%%%%%%%%%%%%%%%%%SUEWS_FunctionalTypes.txt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Open SUEWS_FunctionalTypes
  
  open(55,file=trim(FileInputPath)//trim('SUEWS_FunctionalTypes.txt'),err=204,status='old')
  ProbLine='header'
  read(55,*,err=205)  ! skip header
  ProbLine='alb'
  read(55,*,err=205) (alb(iv),iv=1,nsurf), alb_snow                  ! 1
  
  alBMax_dec=alb(DecidSurf)
  ProbLine='emis'
  read(55,*,err=205) (emis(iv),iv=1,nsurf), emis_snow                ! 2
  ProbLine='BaseT'
  read(55,*,err=205) skip, skip, (baseT(iv),iv=1,nVegsurf),skip,skip   ! 3
  ProbLine='BaseTe'
  read(55,*) skip, skip, (baseTe(iv),iv=1,nVegsurf),skip,skip   ! 4
  ProbLine='! storage capacity minimum'
  read(55,*,err=205) (Surf(1,iv),iv=1,nsurf), skip         ! 5   ! storage capacity minimum/default
  ProbLine='! 6    ! storage capacity maximum'
  read(55,*,err=205) (Surf(5,iv),iv=1,nsurf), skip         ! 6    ! storage capacity maximum
  ProbLine='! 7   ! drainage equation'
  read(55,*,err=205) (Surf(2,iv),iv=1,nsurf), skip         ! 7   ! drainage equation
  ProbLine='! 8    ! drainage coef1'
  read(55,*,err=205) (Surf(3,iv),iv=1,nsurf), skip         ! 8    ! drainage coef1
  ProbLine='! 9    ! drainage coef2'
  read(55,*,err=205) (Surf(4,iv),iv=1,nsurf), skip         ! 9    ! drainage coef2
  ProbLine='!10 GDD full'
  read(55,*,err=205) skip, skip, (GDDFull(iv),iv=1,nVegsurf),skip,skip  ! 10
  ProbLine='!11 SDD full'
  read(55,*,err=205) skip, skip,(SDDFull(iv),iv=1,nVegsurf),skip,skip   ! 11
  ProbLine='!12  LAI min'
  read(55,*,err=205) skip, skip, (LAImin(iv),iv=1,nVegsurf),skip,skip ! 12  ! LAI min
  ProbLine='!13  LAI max'
  read(55,*,err=205) skip, skip, (LAImax(iv),iv=1,nVegsurf),skip,skip  ! 13 ! LAI max
  ProbLine='! 14 Maximum conductance'
  read(55,*,err=205) skip, skip, (MaxConductance(iv),iv=1,nVegsurf),skip,skip  ! 14 Maximum conductance
  ProbLine='!15 Soil store capacity'
  read(55,*,err=205) (soilstoreCap(iv),iv=1,nsurf), skip    !15 Soil store capacity
  ProbLine='!16 Volumetric soil moisture capacity'
  read(55,*,err=205) (VolSoilMoistCap(iv),iv=1,nsurf), skip !16 Volumetric soil moisture capacity
  ProbLine='!17 %Hydraulic conductivity of saturated soil'
  read(55,*,err=205) (SatHydraulicConduct(iv),iv=1, nsurf), skip !17 %Hydraulic conductivity of saturated soil
  ProbLine= 'header'
  read(55,*,err=205)  ! skip header
  ProbLine= '! 19 conductance parameters'
  read(55,*,err=205) G1,G2,G3,G4,G5,G6   ! 19 conductance parameters
  ProbLine='! 20 conductance parameters'
  read(55,*,err=205) TH,TL,S1,S2, Kmax   ! 20 conductance parameters
  ProbLine= 'header'
  read(55,*,err=205)  ! skip header Soil related (measurements)
  ProbLine='Soil related (measurements)'
  read(55,*,err=205) SoilDensity, SoilDepthMeas, SoilRocks,SmCap  !22
  ProbLine= 'header'
  read(55,*,err=205)  ! skip header LUMPS related
  ProbLine='! 23LUMPS related '
  read(55,*,err=205)  DRAINRT, RAINCOVER, RAINMAXRES ! 24 LUMPS related (1)drainage rate (2),for adjusting alpha/beta for wet surface= , (3)Maximum water bucket reservoir
  ProbLine= 'header'
  read(55,*,err=205) ! skip header snow related
  ProbLine='!25 Snow Related'
  read(55,*,err=205) RadMeltFact,TempMeltFact,albSnowMin,albSnowMax,tau_a,tau_f,PrecipLimitAlb !26
  read(55,*,err=205) densSnowMin,densSnowMax,tau_r,CRWmin,CRWmax,PrecipLimit                   !27
  read(55,*,err=205) (snowD(iv),iv=1,nsurf-1)!28
  ProbLine='!28 Snow water hold capacity'
  read(55,*,err=205)  ! skip header !29
  ProbLine=' Narp Related'
  read(55,*,err=205) TRANS_SITE !30
  read(55,*,err=205) ! skip header snow related !31
  ProbLine=' 30 LAI related' ! skip header
  read(55,*,err=205) LAItype, (laiPower(iv),iv=1,4)
  
  !!sg_09Nov13 - remove this
  !Maxima LAI's
 ! MaxLaiMax=max(Laimax(1),LaiMax(2))
!  do iv=3,nVegsurf
!      MaxLaiMax=max(MaxLaiMax,LaiMax(iv))
!  enddo
      
  CapMin_dec=surf(1,DecidSurf)
  CapMax_dec=surf(5,DecidSurf)
  
  GDDmax=0
  SDDmax=0
  do iv=1,NVegSurf
      GDDmax=max(GDDFull(iv),GDDmax)
      SDDmax=min(SDDFull(iv),SDDmax)
  enddo
  close (55) 
  file_qs=.false.
    
 
  !----------------------------------------------------------------------
  !SUEWS run information
  Inputmetformat=10
  hrcount=0
  LAICalcYes=1
  ITY=2
  NPeriodsPerDay= 24 ! number of time periods per day

  !----------------------------------------------------------------------
  write(*,*)'-------------------------------------------------'
  write(*,*)"LUMPS/Suews 2012a - relevant references"
  write(*,*)"LUMPS - Grimmond and Oke (2002) JAM, 41, 79-810"
  write(*,*)"OHM - Grimmond and Oke (1999) JAM, 38, 922-940"
  write(*,*)"NARP - Offerle et al. (2003) JAM"
  write(*,*)"SUES - Evaporation Grimmond & Oke (1991) WRR"
  write(*,*)"Water Balance Model Grimmond et al. (1986) WRR"
  write(*,*)"NARP - long wave improvements (Loridan et al. 2011 JAMC)"
  write(*,*)"SUEWS - anthropogenic heat, etc (Jarvi et al. 2011 JH)"
  write(*,*)'-------------------------------------------------'

  !----------------------------------------------------------------------
  !Write nml files to FileChoices.txt
  
  write(12,*)'--------SUEWS_FunctionalTypes.txt---------------------------- '
  write(12,*)'!Paved  Bldg   Conif  Decid  GrassI   GrassU  Water  Snow        -9 not applicable  '
  write(12,120) (alb(iv),iv=1,nsurf), alb_snow,' albedo -1'             ! 1
  write(12,120) (emis(iv),iv=1,nsurf), emis_snow, 'emis -2'               ! 2
  write(12,120) skip, skip, (baseT(iv),iv=1,nVegsurf),skip,skip ,'BaseT'  ! 3
  write(12,120) skip, skip, (baseTe(iv),iv=1,nVegsurf),skip,skip, 'BaseTe'   ! 4
  write(12,120) (Surf(1,iv),iv=1,nsurf), skip ,'storage capacity minimum/default' !5
  write(12,120) (Surf(5,iv),iv=1,nsurf), skip ,'storage capacity maximum'       ! 
  write(12,120) (Surf(2,iv),iv=1,nsurf), skip ,'drain equation' 
  write(12,120) (Surf(3,iv),iv=1,nsurf), skip ,'dr coef1'     ! 7    ! 
  write(12,120) (Surf(4,iv),iv=1,nsurf), skip ,'dr coef2'      ! 8    !  
  write(12,'(8f8.0, 2g10.4)') skip, skip, (GDDFull(iv),iv=1,nVegsurf),skip,skip, 'GDDFull '  ! 10
  write(12,'(8f8.0, 2g10.4)') skip, skip,(SDDFull(iv),iv=1,nVegsurf),skip,skip,'SDDFull'  ! 11
  write(12,120) skip, skip, (LAImin(iv),iv=1,nVegsurf),skip,skip,'LAI min'! 12  !
  write(12,120) skip, skip, (LAImax(iv),iv=1,nVegsurf),skip,skip,'LAI max'  ! 13 ! 
  write(12,120) skip, skip, (MaxConductance(iv),iv=1,nVegsurf),skip,skip,'MaxCond'  ! 14
  write(12,'(8f8.0, 2g10.4)') (soilstoreCap(iv),iv=1,nsurf), skip, 'soilstoreCap'     ! 15
  write(12,120) skip, skip, (VolSoilMoistCap(iv),iv=1,nVegsurf),skip,skip,'VolSoilMoistCap'  ! 16
  write(12,120) (SatHydraulicConduct(iv),iv=1,nsurf),skip,'SatHydraulicConduct'  ! 17
  write(12,'(5g10.4, g10.5, g20.0)') G1,G2,G3,G4,G5,G6,' conductance parameters'   ! 17 
  write(12,'(4g10.2,g10.5,g20.0)') TH,TL,S1,S2, Kmax,' conductance parameters'   ! 18 
  write(12,*)  ! skip header Soil rlated
  write(12,'(g12.6,g10.2,f6.0,3g8.3,g8.1)')SoilDensity,SoilDepthMeas, SoilRocks,SmCap,'Soil'  !24
  write(12,*)  ! skip header LUMPS related
  write(12,'(10g8.2)')DRAINRT,RAINCOVER,RAINMAXRES,'LUMPS (1)drainage rate,adjust alpha/beta wet surface(3)Max water bucket'  ! 26
  write(12,*)! skip header snow related  
  write(12,'(10g8.2)')RadMeltFact,TempMeltFact,albSnowMin,albSnowMax,tau_a,tau_f,PrecipLimitAlb
  write(12,'(10g8.2)')densSnowMin,densSnowMax,tau_r,CRWmin,CRWmax,PrecipLimit    
  write(12,*)! skip header NARP related  
  write(12,'(10g8.2)') TRANS_SITE, 'tran_site'
  write(12,*)! skip header LAI related  
  write(12,'(10g8.2)') LAItype, (laiPower(iv),iv=1,2)
  
  120	 format (10g10.2)
 
  !Interval (sec),TStep (s),NSH number of steps per interval
  
  NSH_real = interval/TSTEP 

  if (NSH_real<2) then    !If NSH_real smaller that two, only one period is run
    NSH=1
    TSTEP=real(interval,kind(1d0))
  else
    if (NSH_real==int(NSH_real)) then
    	NSH=interval/TSTEP    !NSH needs to be an integer
        
    else
      	call ErrorHint(39,'File: RunControl',TSTEP,real(interval,kind(1d0)),notUsedI)
      
    endif
  endif
  
  NMIN=TSTEP/60
  
  NPeriodsPerYear= NPeriodsPerDay*NDays ! assume 366 max but could be less  
  
  halftimestep=real(interval)/2/(24*3600) ! used in sun_position to get sunpos in the middle of timestep
 
  ! ---------------OHMChoices-----------------------------------------------
  ! coefficients for storage heat flux (will select by grid later)
  fileOHMChoices=trim(FileInputPath)//'OHM_Coefficients.txt'
  co=0
  ! ---------------OHMChoices-----------------------------------------------
  open(8,file=trim(fileOHMChoices),status='old',err=202)
  do i=1,NrowOhm
	 READ(8,*,iostat=iostat_var) (co(i,iv),iv=1,NCoeffOhm)
  enddo
  close (8)
  write(12,*) '---Number of rows in OHM_Coefficients.txt =',i-1

  return
  !-------problems---------------------------------------------------------------------------------------
200		call ProblemsText('RunControl.nml -file not found')
		call PauseStop
201  	call ProblemsText('RunControl.nml - problem in file')  
		write(500,nml=RunControl)
        call PauseStop
202		call ProblemsText(trim(fileOHMChoices))
 		call PauseStop
203     call ProblemsText(trim(FileChoices))
 		call PauseStop
204     call ProblemsText(trim(FileInputPath)//trim('SUEWS_FunctionalTypes.txt'))
 		call PauseStop       
205		call ProblemsText(trim(FileInputPath)//trim('SUEWS_FunctionalTypes.txt')// ProbLine)
		call PauseStop
219     call ErrorHint(12,trim('SnowControl'),notUsed,notUsed,ios_out)
        call PauseStop         
     
end subroutine OverallRunControl
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine OpenAnnualFiles(Grid)
	!This subroutine opens the annual files and saves header to those
    use allocateArray
    use data_in
    
	IMPLICIT NONE

    character(len=15)::grid

    if(CreateAnnual==1)then
           if(int(year)==FirstYear)then
               open(61,file=trim(FileOutputPath)//trim(Grid)//'_annual.txt',err=200)
               write(61,140)
           else
               open(61,file=trim(FileOutputPath)//trim(Grid)//'_annual.txt',access="append",err=200)
           endif

        endif
        if(alldays==1) then
           if(int(year)==FirstYear)then
               open(90,file=trim(FileOutputPath)//trim(Grid)//'_alldays.txt',err=201)
               write(90,141)
           else
               open(90,file=trim(FileOutputPath)//trim(Grid)//'_alldays.txt',access="append",err=201)
           endif
        endif
        
        if(writedailyState==1)then
          if(int(year)==FirstYear)then
               open(60, file=trim(FileOutputPath)//trim(Grid)//'DailyState.txt',err=202)
               write(60,142)
          else     
               open(60, file=trim(FileOutputPath)//trim(Grid)//'DailyState.txt',access="append",err=202)   
          endif
        endif

140 	format('%yr counter      qn      qs       qf     qe_S       pp      ext_Ie',&
        '     int_Ie   tot_ie      E_S     Change     R_Soil      R         Fw ',&
        '    addWater   QH_S      Qm   delta_QSI   Qrain    SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i') 

141 	format('%yr     day     qn      qs       qf     qe_S       pp      ext_Ie',&
        '     int_Ie   tot_ie      E_S     Change     R_Soil      R         Fw ',&
        '    addWater   QH_S      Qm   delta_QSI   Qrain    SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i')

142      format('%year id    HD1h  HDD2c  HDD3m  HDT5d   Prec   DaSR  GDD1g GDD2s  GDmn   GDmx ',&
       ' dayLG   LAIc   LAId  LAIgI  LAIgU  DEcap    Por Albdec WUgr(1) WUgr(2) WUgr(3)',&
       ' WUtr(1) WUtr(2) WUtr(3) LAIch LAIlumps alb_snow dens_snow_pav dens_snow_bldg',&
       ' dens_snow_evergr dens_snow_dec dens_snow_Irrgr dens_snow_Gr dens_snow_water')
       
		return

200 	call ProblemsText(trim(FileOutputPath)//trim(Grid)//'_annual.txt')
		call PauseStop 

201 	call ProblemsText(trim(FileOutputPath)//trim(Grid)//'_alldays.txt')
        call PauseStop 

202 	call ProblemsText(trim(FileOutputPath)//trim(Grid)//'DailyState.txt')
		call PauseStop 
end subroutine OpenAnnualFiles

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine RunControlByGridByYear    
 use data_in       ! LUMPS_metRead.f95
 use ohm_calc
 use run_info      ! run_control_v2_1.f90
 use SUES_data     !
 USE GIS_data      ! LUMPS_gis_read.f95
 use mod_z     ! module_LUMPS_constants,f90
 use FileName
 use allocateArray  ! module_LUMPS_constants,f90
 use time
 use defaultNotUsed
 use snowMod
 
 IMPLICIT NONE
 integer:: which,iv,j,ic
 real(kind(1d0))::rate(9)
 character (len=150):: CanopyName,FileHr !Name of the input files
 
 namelist/SiteSpecificParam/lat,&
 							lng,& 
                            SurfaceArea,&
                            RunoffToWater,&  
                            WaterUseAreaGrass,&
                            WaterUseAreaTrees,&
                            FlowChange,& 
                            PipeCapacity,&
 							Faut,&          !Irrigation
                            Ie_end,&  
                            Ie_start,&
                            Ie_a,&
                            Ie_m,&
                            DayWat,& 
                            DayWatPer,&   
                            InternalWaterUse,&  ! in mm (/24)    
                            SnowLimBuild,&
				            SnowLimPaved				                
                         
                        
 FileWU=trim(FileInputPath)//'SiteSpecificParam'//trim(FileCode)//'.nml'
 FileHr=trim(FileInputPath)//'HourlyProfile'//trim(FileCode)//'.txt'

 open(7,File=FileHr,err=251,status='old')
 read(7,*,err=251,iostat=ios_out) !skip header
 ic=1
 read(7,*,err=250,iostat=ios_out)(HourWat(j),j=0,23)
 ic=2
 read(7,*,err=250,iostat=ios_out)(AhProF(j,1),j=0,23)
 ic=3
 read(7,*,err=250,iostat=ios_out)(AhProF(j,2),j=0,23)
 ic=4
 read(7,*,err=250,iostat=ios_out)(snowProf(j),j=0,23)
 close (7)

 open(7,File=FileWU,err=200,status='old') !Change with needs
 read(7,nml=SiteSpecificParam,err=211, iostat=ios_out)
 close(7)

 write(12,*)'----------',trim(FileHr),'-hr:0 to 23-- after HourResChoice--------'
 write(12,'(24f6.2, a10)')HourWat, ' Hour Wat'
 write(12,'(24f6.2,a20)')(AhProF(j,1),j=0,23),'AhProf weekday -1 '
 write(12,'(24f6.2,a20)')(AhProF(j,2),j=0,23),'AhProf weekend -2 '
   
 if(abs(sum(HourWat)-1.000)>0.001)then
  if(abs(sum(HourWat)-1.000)<0.02)then
    	HourWat=HourWat*(1/sum(HourWat))
  else
  		call ErrorHint(7,trim(FileHr),sum(HourWat),notUsed,notUsedI)
  endif
 endif
 SurfaceArea=SurfaceArea*10000!Change surface area from ha to m^2
 write(12,*)'----------',trim(FileWU),'----------'
 write(12,nml=SiteSpecificParam)
 
 !------------------------------------------------------------------------------------------                    
 CanopyName=trim(FileInputPath)//'WaterDist'//trim(FileCode)//'.txt'  
 open(6,File=trim(CanopyName),err=201,status='old') 
 read(6,*,iostat=ios_out) ! skip header
 if(ios_out<0)call ErrorHint(11,trim(canopyname),notUsed,notUsed,ios_out)
 
! initialize to all zero
 WaterDist=0 !Table which determines distribution of water in the canopy
 do iv=1, Nsurf-1 
 	read(6,*,iostat=ios_out,end=501)(rate(j),j=1,9),which
       
    if(ios_out/=0)then
		call ProblemsText(trim(CanopyName))
       	write(500,*)iv,' number of lines read but expecting ',nsurf-1
      	write(500,*)ios_Out, " [if this is -1 EOF; -2 EOR]"
        call PauseStop	  
 	endif
 501 	if(rate(which)/=0) call ErrorHint(8,trim(CanopyName),rate(which),notUsed,notUsedI)
    	
  ! check LJ why can this not go to soil and runof
  	if(rate(8)/=0.and.rate(9)/=0)  call ErrorHint(9,trim(CanopyName),rate(8),rate(9),notUsedI)
    if(sum(rate)>1.0000001.or.sum(rate)<0.9999999)call ErrorHint(10,trim(CanopyName),sum(rate),notUsed,which)
      	
	 do j=1,Nsurf  
		 WaterDist(j,which)=rate(j) 
	 enddo
	 if(rate(8)/=0) then
 		WaterDist(8,which)=rate(8)
	 else
 	  	WaterDist(8,which)=rate(9)
	 endif
  enddo
  
   write(12,*)'waterDistribution'
   write(12,*)'To !Paved Build  Conif  Decid GrassIrr GrassUnirr  Water Runoff  Soil SurfaceFrom   '
  
   do iv=1,nSurf-1
        write(12,'(i4, 8f6.2)')iv,(WaterDist(j,iv),j=1,nsurf+1)
   enddo
        
   close(6)                                
return
200		call ProblemsText(trim(FileWU))
		call PauseStop 
201		call ProblemsText(trim(CanopyName))
		call PauseStop         
251		call ErrorHint(11,trim(FileHr),notUsed,notUsed,ios_out)
250     reall=real(ios_out)
		call ErrorHint(13,trim(FileHr),reall,notUsed,ic)
211 	write(*,nml=SiteSpecificParam) 
 		call ErrorHint(12,trim(FileWU),daywatper(1),daywatper(2),ios_out)
		write(500,nml=SiteSpecificParam) 
        call PauseStop 

end subroutine RunControlByGridByYear          
!----------------------------------------------------------------------------------------------
! sg feb 2012 - 
!----------------------------------------------------------------------------------------------

 subroutine InitialState(GridName,errFileYes)
  use resist        ! LUMPS_metRead.f95
  use data_in       ! LUMPS_metRead.f95
  use ohm_calc
  use run_info      ! run_control_v2_1.f90
  use SUES_data     !
  USE GIS_data      ! LUMPS_gis_read.f95
  use mod_z         ! module_LUMPS_constants,f90
  use FileName
  use allocateArray  ! module_LUMPS_constants,f90
  use time
  use defaultNotUsed
  use snowMod
  use InitialCond
  
  IMPLICIT NONE
  
  character (len=15)::GridName
  character (len=10)::str2 !Variables related to filepaths 
  real (KIND(1d0))::BldgState,ConifState,DecidState,GrassIState,GrassUState,PavState,&
              	soilstoreDecState,soilstoreGrassUnirState,soilstoreConifstate,&
          		soilstoreGrassIrrState,soilstorePavstate,soilstoreBldgState,SnowFracBldg,&
                SnowFracConif, SnowFracDec, SnowFracGrassIrr, SnowFracGrassUnir,&
                SnowFracPav, SnowFracWater,SnowDensBldg,SnowDensConif,SnowDensDec,&
                SnowDensGrassIrr,SnowDensGrassUnir,SnowDensPav,SnowDensWater  
                    
  character (len=4)::year_txt   
  character(len=150)::fileInit 
  integer:: j,DaysSinceRain,wd,seas,date,mb,year_int,switch=0,id_next,calc,errFileYes

  !These relate to InitialConditions.nml and are only used in this part of the code

  namelist/InitialConditions/DaysSinceRain,&
                  Temp_C0,&
                  ID_Prev,&
                  GDD_1_0,&
                  GDD_2_0,&
                  BldgState,&
                  ConifState,&
                  DecidState,&
                  PavState,&
                  GrassIState,&
                  GrassUState,&
                  LAIinitialConif,&            ! 
                  LAIinitialDecid,&
                  LAIinitialGrassU,&
                  LAIinitialGrassI,&
                  porosity0,&
                  DecidCap0,&
                  albDec0,&
                  soilstoreBldgState,&
                  soilstoreConifstate,&
                  soilstoreDecState,&
                  soilstoreGrassUnirState,&
                  soilstoreGrassIrrState,&
                  soilstorePavstate,&
                  WaterState,&
				  SnowWaterBldgState,&
				  SnowWaterConifstate,&
                  SnowWaterDecState,&
                  SnowWaterGrassIrrState,&
                  SnowWaterGrassUnirState,&
                  SnowWaterPavstate,&
                  SnowWaterWaterstate,&
				  SnowPackBldg,& 
                  SnowPackConif,&
                  SnowPackDec,&
                  SnowPackGrassIrr,&
                  SnowPackGrassUnir,&
                  SnowPackWater,&
                  SnowPackPav,&
                  SnowFracBldg,&
                  SnowFracConif,&
                  SnowFracDec,&
                  SnowFracGrassIrr,&
                  SnowFracGrassUnir,&
                  SnowFracPav,&
                  SnowFracWater,&
                  SnowDensBldg,&
                  SnowDensConif,&
                  SnowDensDec,&
                  SnowDensGrassIrr,&
                  SnowDensGrassUnir,&
                  SnowDensPav,&
                  SnowDensWater

  year_int=int(year)
  write(year_txt,'(I4)')year_int

  FileInit=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//trim(year_txt)//'.nml'
  call ErrorHint(44,FileInit,notUsed,notUsed,notUsedI)
   
  open(55,File=trim(FileInit),err=200,status='old') !Change with needs
  read(55,iostat=ios_out,nml=InitialConditions,err=203)
  close(55)
  write(12,nml=InitialConditions)
  
  if(id_prev>=364)id_prev=0
    
  changed=0 ! use for LUMPS phenology
  changeInit=0
  ! check -- this needs to be set for passing between years
  !'-2' because array nvegsurf
  lai=0
  lai(id_prev,ivConif)=LAIinitialConif
  lai(id_prev,ivdecid)=LAIinitialDecid
  lai(id_prev,ivGrassI)=LAIinitialGrassI
  lai(id_prev,ivGrassU)=LAIinitialGrassU
  
  ! growing degree days
  GDD(:,1)=0
  GDD(:,2)=0
  GDD(:,5)=0
  GDD(:,3)=90 ! going to check for minimum
  GDD(:,4)=-90 ! going to check for maximum
  GDD(id_prev,1)=GDD_1_0
  GDD(id_prev,2)=GDD_2_0
  
  jj1=0
  jj2=0
  jj3=0
  jj4=0
  HDD=0
  if(AnthropHeatChoice>=0)then !!BaseT is used always.??Check
    ! sg 18/04/12 - file name changed
    FileSAHP=trim(FileInputPath)//trim(FileCode)//'SAHP.nml'  
    call SAHP_Coefs(temp_C0,id_prev)
  endif
  ! assume that the temperature has been the same for the previous days
  HDD(id_prev-3,3)=Temp_C0
  HDD(id_prev-2,3)=Temp_C0
  HDD(id_prev-1,3)=Temp_C0
  HDD(id_prev,3)=Temp_C0
  HDD(id_prev,6)=daysSinceRain

  !Initialize hourly temperature and precipitation + 5 day mean used 
  !to define thermal growing season
  runT(0:23)=Temp_C0
  runP(0:23)=0
  
  ! initialize all initally to the previous day value (i.e. day before run starts)
  porosity=porosity0
  AlbDec=albDec0
  decidCap=DecidCap0
  CumSnowfall=0
  
  !------DEFINE INPUT FILE PATHS AND FILE NAMES------------------------------
  FileMet=trim(FileInputPath)//trim(FileCode)//'_data.txt'
  FileGIS=trim(FileInputPath)//trim(FileCode)//'.gis'
  FileOHM=trim(FileInputPath)//'SelectOHM'//trim(FileCode)//'.txt'

            !===============READ GIS DATA  if not varying ==================================================================
  open(3,file=trim(fileGIS),status='old',err=313)
  call skipHeader(3,SkipHeaderGIS)
  finish=.false.
   ! GISInputType=4   Varies each time step
   ! GISInputType=3   Stays the same each hour
   
   IF(gisInputType==3) then
      call read_gis(finish)
          if(z0_method==1) then
               if(z0m<0.00001) call ErrorHint(5,trim(fileGIS),z0m,notUsed,notUsedI)
               if(zdm<0.00001) call ErrorHint(6,trim(fileGIS),zdm,notUsed,notUsedI)
                zzd=z-zdm
          elseif(z0_method==3)then
                if(FAIBLdg<0) call ErrorHint(1,trim(fileGIS),Faibldg,notUsed,notUsedI)
                if(FAITree<0)call ErrorHint(2,trim(fileGIS),faitree,notUsed,notUsedI)
          else   
               call RoughnessParameters(id_prev)
          endif
      close(3)
   elseIF(gisInputType==4) then
      call read_gis(finish)
          if(z0_method==1) then
               if(z0m<0.00001) call ErrorHint(5,trim(fileGIS),z0m,notUsed,notUsedI)
               if(zdm<0.00001) call ErrorHint(6,trim(fileGIS),zdm,notUsed,notUsedI)
                zzd=z-zdm
          elseif(z0_method==3)then
                if(FAIBLdg<0) call ErrorHint(1,trim(fileGIS),Faibldg,notUsed,notUsedI)
                if(FAITree<0)call ErrorHint(2,trim(fileGIS),faitree,notUsed,notUsedI)
          else   
               call RoughnessParameters(id_prev)
          endif
      backspace(3)
 
   endif

if(id_prev==0)then
  year=year-1
  call LeapYearCalc (year,id_prev)
  switch=1
!  print*,'switch'
endif 
  
call day2month(id_prev,mb,date,seas,year,lat)!Calculate real date from doy
call Day_of_Week(date,mb,year,wd)!Calculate weekday (1=Sun,...)

if(switch==1)then
  year=year+1
  id_prev=0
  switch=0
endif
  
dayofWeek(id_prev,1)=wd  ! day of week
dayofWeek(id_prev,2)=mb  ! month
dayofweek(id_prev,3)=seas ! season

! in case first hour not 0
id_next=id_prev+1
if(id_next>nofDaysThisYear) then
  id_next=1
  year=year+1
  switch=1
  call ErrorHint(43,'switch- years',notUsed,notUsed,notUsedI)
endif
call day2month(id_next,mb,date,seas,year,lat)!Calculate real date from doy
call Day_of_Week(date,mb,year,wd)!Calculate weekday (1=Sun,...)
if(switch==1)then
  year=year-1
  switch=0
endif
dayofWeek(id_next,1)=wd  ! day of week
dayofWeek(id_next,2)=mb  ! month
dayofweek(id_next,3)=seas ! season


if(switch==1)then
  year=year+1
  id_prev=0
  switch=0
endif

WU_day=0
id=id_prev
it=LastTimeOfDay
if (id_prev>=DayLightSavingDay(1).and.id_prev<=DayLightSavingDay(2)) then!Summer time
    DLS=1
else
    DLS=0
endif
        

!Calculate daily water use if this is modelled
IrrTrees =0
if (WU_choice==0) then
  ! check if
  calc=0
  if (DayWat(wd)==1.0) then      !if=1 - then this is a day that has watering
      if (lat>=0)then             !Northern Hemisphere
          if (id>=Ie_start.and.id<=Ie_end) calc=1 !if day between irrigation period               
      else                        !Southern Hemisphere
          calc=1
          if (id>=Ie_end.and.id<=Ie_start) calc=0 !if day between irrigation period                       
      endif
      if(calc==1) then              
          ! HDD(id,6) -- daysSincerain, HDD(id3)- mean airT
          ! automatic

          WU_day(id,2)=Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(GrassISurf)*DayWatPer(wd) 
          ! manual
          WU_day(id,3)=(1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(GrassISurf)*DayWatPer(wd)
                      !aut=(-84.535+9.959*TempWU_Daily(i)+3.674*DaysSinceRain(i))*Faut
                      !man=(-23.36+2.988*TempWU_Daily(i)+1.102*DaysSinceRain(i))*(1-Faut)
          if (WU_Day(id,2)<0) WU_Day(id,2)=0 !If modelled WU is negative -> 0
          if (WU_Day(id,3)<0) WU_Day(id,3)=0 !If modelled WU is negative -> 0
                              
          WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))
          
          !Calculate the fraction for irrigated trees/shrubs. Added by LJ in 9 September 2013
          IrrTrees = sfr(ConifSurf)*IrrFractionTrees+sfr(DecidSurf)*IrrFractionTrees

          WU_day(id,5)=Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*IrrTrees*DayWatPer(wd) 
          WU_day(id,6)=(1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*IrrTrees*DayWatPer(wd)
                           
            
          if (WU_Day(id,5)<0) WU_Day(id,5)=0 !If modelled WU is negative -> 0
          if (WU_Day(id,6)<0) WU_Day(id,6)=0 !If modelled WU is negative -> 0
                              
          WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))
 
        else
           WU_Day(id,1)=0
           WU_Day(id,2)=0
           WU_Day(id,3)=0
           WU_Day(id,4)=0
           WU_Day(id,5)=0
           WU_Day(id,6)=0
       endif
    endif
 endif

 
     !Wetness status of each surface  

!Above ground state
State(PavSurf)=PavState
State(BldgSurf)=BldgState
State(ConifSurf)=ConifState
State(DecidSurf)=DecidState
State(GrassISurf)=GrassIState
State(GrassUSurf)=GrassUState
State(WaterSurf)=WaterState   !State of water body

!Maximum capasity soil storage can hold for each surface type

!::::Below ground state::::::::::::::::::::::::::::::::::::::::::
!Max capacity of the different soil levels

 !Soil moisture of each surface type
   
!Initial conditions for soilstores
soilmoist(PavSurf)=soilstorePavstate
soilmoist(BldgSurf)=soilstoreBldgState
soilmoist(ConifSurf)=soilstoreConifstate
soilmoist(DecidSurf)=soilstoreDecState
soilmoist(GrassISurf)=soilstoreGrassIrrState
soilmoist(GrassUSurf)=soilstoreGrassUnirState
soilmoist(WaterSurf)=0                 !No soil layer for water body


!!SNOWPACK PARAMETERS
!Initial parameters for snowpack 
SnowPack(PavSurf)=SnowPackPav
SnowPack(BldgSurf)=SnowPackBldg
SnowPack(ConifSurf)=SnowPackConif
SnowPack(DecidSurf)=SnowPackDec
SnowPack(GrassISurf)=SnowPackGrassIrr
SnowPack(GrassUSurf)=SnowPackGrassUnir
SnowPack(WaterSurf)=SnowPackWater

!Amount of liquid (melted water) in the snowpack
MeltWaterStore(PavSurf) = SnowWaterPavstate
MeltWaterStore(BldgSurf) = SnowWaterBldgState
MeltWaterStore(ConifSurf) = SnowWaterConifstate
MeltWaterStore(DecidSurf) = SnowWaterDecState
MeltWaterStore(GrassISurf) = SnowWaterGrassIrrState
MeltWaterStore(GrassUSurf) = SnowWaterGrassUnirState
MeltWaterStore(WaterSurf) = SnowWaterWaterstate

!Initial fraction of snow
snowFrac(PavSurf)=SnowFracPav
snowFrac(BldgSurf)=SnowFracBldg
snowFrac(ConifSurf)=SnowFracConif
snowFrac(DecidSurf)=SnowFracDec
snowFrac(GrassISurf)=SnowFracGrassIrr
snowFrac(GrassUSurf)=SnowFracGrassUnir
snowFrac(WaterSurf)=SnowFracWater

iceFrac=0.2


densSnow(PavSurf) = SnowDensPav    !Initialize snow density
densSnow(BldgSurf) = SnowDensBldg 
densSnow(ConifSurf) = SnowDensConif
densSnow(DecidSurf) = SnowDensDec
densSnow(GrassISurf) = SnowDensGrassIrr
densSnow(GrassUSurf) = SnowDensGrassUnir
densSnow(WaterSurf) = SnowDensWater

! for first 3 h Q* just assume at night
! check ?? this should not happen at the start of the year if continuing

q1=-101
q2=-100
q3=-99

r1=-101
r2=-100
r3=-99


 !%%%%%%%%%%OUTPUT FILENAMES AND THEIR PATHS%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Write(str2,'(i2)') Interval/60
 FileOut=trim(FileOutputPath)//trim(FileCode)//'_'//trim(adjustl(str2))//'.txt'
 FileErrorInf=trim(FileOutputPath)//trim(FileCode)//'_ErrorFile.txt'
 NARPOut=trim(FileOutputPath)//trim(FileCode)//'_NARPOut.txt'
 fileMonthly=trim(FileOutputPath)//trim(FileCode)//'_MonthlyFile.txt'
 fileDaily=trim(FileOutputPath)//trim(FileCode)//'_DailyFile.txt'
 file5min=trim(FileOutputPath)//trim(FileCode)//'_5min.txt'
 SnowOut=trim(FileOutputPath)//trim(FileCode)//'_SnowOut.txt'
 SOLWEIGpoiOut=trim(FileOutputPath)//trim(FileCode)//'_SOLWEIGpoiOut.txt'
 BLOut=trim(FileOutputPath)//trim(FileCode)//'_BL.txt'

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !DEFINE DIFFERENT INITIALIZATION PARAMETERS
 once=.true.

 lfn_us=10
 lfnOut=35  ! Output error file

 
!==================OUTPUT FILE OPTIONS=============================================
!Open daily and monthly
   open(14,file=fileDaily,err=201)
   open(15,file=fileMonthly,err=202)


 ! zero arrays------------------------------------------------------------
 !sg -- 16 variables written out

 do j=1,NumberDailyVars
   all_tot(1,j)=0  
 enddo
call accum_zero    
! zero arrays------------------------------------------------------------
    
return
200		call ProblemsText(trim(FileInit))
		call PauseStop
201		call ProblemsText(trim(fileDaily))
		call PauseStop
202		call ProblemsText(trim(fileMonthly))
		call PauseStop
203 	call ErrorHint(16,trim(FileInit),notUsed,notUsed,ios_out)
      	write(500,nml=InitialConditions)       
        call PauseStop
111 	call ProblemsText(trim(FileErrorInf))
		call PauseStop
313 	call ProblemsText(trim(filegis))
		call PauseStop

 end subroutine InitialState

!--------------------------------------------------------------------------

  
subroutine accum_zero
use data_in
implicit none
   day=0
   month=0
   yr_tot=0
   season=0  

return
end subroutine accum_zero

!--------------------------------------------------------------------------

subroutine NextInitial(GridName)
use allocateArray  ! module_LUMPS_constants,f90
use time
use data_in
use snowMod

IMPLICIT NONE
character (len=15)::GridName
character (len=4)::year_txt
!real(Kind(1d0))::year
integer:: year_int

if (id>360) then 
  year_int=int(year+1)
  write(year_txt,'(I4)')year_int
  open(55,File=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//trim(adjustl(year_txt))//'.nml',err=200) 

else
  year_int=int(year)
  write(year_txt,'(I4)')year_int
  open(55,File=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//trim(adjustl(year_txt))//'end'//'.nml',err=201) 
endif

write(55,*)'&InitialConditions'
write(55,*)'DaysSinceRain=',int(HDD(id,6))
if(it/=LastTimeofday)then
  ! need to do this otherwise it will not be completed
	id=id-1  
endif

write(55,*)'Temp_C0=',HDD(id,3)  
write(55,*)'ID_Prev=',id
write(55,*)'GDD_1_0=',GDD(id,1)
write(55,*)'GDD_2_0=',GDD(id,2)
write(55,*)'BldgState=',State(BldgSurf)
write(55,*)'ConifState=',State(ConifSurf)
write(55,*)'DecidState=',State(DecidSurf)
write(55,*)'PavState=',State(PavSurf)
write(55,*)'GrassIState=',State(GrassISurf)
write(55,*)'GrassUState=',State(GrassUSurf)
write(55,*)'LAIinitialConif=',lai(id,ivConif)         ! 
write(55,*)'LAIinitialDecid=',lai(id,ivdecid)
write(55,*)'LAIinitialGrassU=',lai(id,ivGrassI)
write(55,*)'LAIinitialGrassI=',lai(id,ivGrassU)
write(55,*)'porosity0=',porosity(id)
write(55,*)'DecidCap0=',decidCap(id)
write(55,*)'albDec0=',AlbDec(id)
write(55,*)'soilstoreBldgState=',soilmoist(BldgSurf)
write(55,*)'soilstoreConifstate=',soilmoist(ConifSurf)
write(55,*)'soilstoreDecState=',soilmoist(DecidSurf)
write(55,*)'soilstoreGrassUnirState=',soilmoist(GrassISurf)
write(55,*)'soilstoreGrassIrrState=',soilmoist(GrassUSurf)
write(55,*)'soilstorePavstate=',soilmoist(PavSurf)
write(55,*)'WaterState=',State(WaterSurf)
write(55,*)'SnowWaterBldgState=',MeltWaterStore(BldgSurf)
write(55,*)'SnowWaterConifstate=',MeltWaterStore(ConifSurf)
write(55,*)'SnowWaterDecState=',MeltWaterStore(DecidSurf)
write(55,*)'SnowWaterGrassIrrState=',MeltWaterStore(GrassISurf)
write(55,*)'SnowWaterGrassUnirState=',MeltWaterStore(GrassUSurf)
write(55,*)'SnowWaterPavstate=',MeltWaterStore(PavSurf)
write(55,*)'SnowWaterWaterstate=',MeltWaterStore(WaterSurf)
write(55,*)'SnowPackBldg=',SnowPack(BldgSurf)
write(55,*)'SnowPackConif=',SnowPack(ConifSurf)
write(55,*)'SnowPackDec=',SnowPack(DecidSurf)
write(55,*)'SnowPackGrassIrr=',SnowPack(GrassISurf)
write(55,*)'SnowPackGrassUnir=',SnowPack(GrassUSurf)
write(55,*)'SnowPackPav=',SnowPack(PavSurf)
write(55,*)'SnowPackWater=',SnowPack(WaterSurf)
write(55,*)'SnowFracBldg=',SnowFrac(BldgSurf)
write(55,*)'SnowFracConif=',SnowFrac(ConifSurf)
write(55,*)'SnowFracDec=',SnowFrac(DecidSurf)
write(55,*)'SnowFracGrassIrr=',SnowFrac(GrassISurf)
write(55,*)'SnowFracGrassUnir=',SnowFrac(GrassUSurf)
write(55,*)'SnowFracPav=',SnowFrac(PavSurf)
write(55,*)'SnowFracWater=',SnowFrac(WaterSurf)
write(55,*)'SnowDensBldg=',densSnow(BldgSurf)
write(55,*)'SnowDensConif=',densSnow(ConifSurf)
write(55,*)'SnowDensDec=',densSnow(DecidSurf)
write(55,*)'SnowDensGrassIrr=',densSnow(GrassISurf)
write(55,*)'SnowDensGrassUnir=',densSnow(GrassUSurf)
write(55,*)'SnowDensPav=',densSnow(PavSurf)
write(55,*)'SnowDensWater=',densSnow(WaterSurf)
write(55,*)'/'
close(55)
return

200 	call ProblemsText(trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//trim(adjustl(year_txt))//'.nml')
		call PauseStop
201     call ProblemsText(trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//trim(adjustl(year_txt))//'end'//'.nml')
		call PauseStop
        
end subroutine NextInitial


 !¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
 subroutine CheckInitial  
 !Check that the parameters in InitialConditions file.
 !Called from metRead. Added by LJ in 8/2/2013
 
 use data_in
 use InitialCond
 use defaultNotUsed
 use time
 use allocateArray
 use snowMod

 
 implicit none

 if (Temp_C0<(Temp_C-10).or.Temp_C0>(Temp_C+10)) then
     call ErrorHint(36,'InitialCond: Check temperature', Temp_C0, Temp_C, notUsedI)
 endif
 if (ID_Prev/=id-1) then
   call ErrorHint(36,'InitialCond: Check previous day', real(ID_Prev,kind(1d0)), real(id,kind(1d0)), notUsedI)
 endif

 !Check if LAI values are OK. Need to treat different hemispheres as well as 
 !tropics separately.
 if (lat>40) then
   if ((LAIinitialConif>LAImin(ConifSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialConif<LAImax(ConifSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialConif in InitialConditions file', LAIinitialConif, LAImin(ConifSurf), notUsedI)
   endif
   if ((LAIinitialDecid>LAImin(DecidSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialDecid<LAImax(DecidSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialDecid in InitialConditions file', LAIinitialDecid, LAImin(DecidSurf), notUsedI)
   endif
   if ((LAIinitialGrassU>LAImin(GrassUSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialGrassU<LAImax(GrassUSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrassU in InitialConditions file', LAIinitialGrassU, LAImin(GrassUSurf), notUsedI)
   endif
   if ((LAIinitialGrassI>LAImin(GrassISurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialGrassI<LAImax(GrassISurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrassI in InitialConditions file', LAIinitialGrassI, LAImin(GrassISurf), notUsedI)
   endif    
   
 elseif (lat<-40) then
   if ((LAIinitialConif<LAImax(ConifSurf-2)-1.and.(id<60.or.id>330)).or.&
     (LAIinitialConif>LAImin(ConifSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialConif in InitialConditions file', LAIinitialConif, LAImax(ConifSurf), notUsedI)
   endif
   if ((LAIinitialDecid>LAImax(DecidSurf-2)-1.and.(id<60.or.id>330)).or.&
     (LAIinitialDecid>LAImin(DecidSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialDecid in InitialConditions file', LAIinitialDecid, LAImax(DecidSurf), notUsedI)
   endif
   if ((LAIinitialGrassU<LAImax(GrassUSurf-2)-1.and.(id<60.or.id>330)).or.&
     (LAIinitialGrassU>LAImin(GrassUSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrassU in InitialConditions file', LAIinitialGrassU, LAImax(GrassUSurf), notUsedI)
   endif
   if ((LAIinitialGrassI<LAImax(GrassISurf-2)-1.and.(id<60.or.id>330)) .or.&
     (LAIinitialGrassI>LAImin(GrassISurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrassI in InitialConditions file', LAIinitialGrassI, LAImax(GrassISurf), notUsedI)
   endif    
 
 elseif (lat<10.and.lat>-10) then
 
   if (LAIinitialConif<LAImax(ConifSurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialConif in InitialConditions file', LAIinitialConif, LAImax(ConifSurf), notUsedI)
   endif
   if (LAIinitialDecid<LAImax(DecidSurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialDecid in InitialConditions file', LAIinitialDecid, LAImax(DecidSurf), notUsedI)
   endif
   if (LAIinitialGrassU<LAImax(GrassUSurf-2)-0.5) then
      call ErrorHint(37,'ICheck LAIinitialGrassU in InitialConditions file', LAIinitialGrassU, LAImax(GrassUSurf), notUsedI)
   endif
   if (LAIinitialGrassI<LAImax(GrassISurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialGrassI in InitialConditions file', LAIinitialGrassI, LAImax(GrassISurf), notUsedI)
   endif
  
 endif

 !Soilstore check
 if (soilmoist(BldgSurf)>soilstoreCap(BldgSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of building soil store.',&
                   soilmoist(BldgSurf), soilstoreCap(BldgSurf), notUsedI)
 endif
 if (soilmoist(PavSurf)>soilstoreCap(PavSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of paved soil store.',&
                   soilmoist(PavSurf), soilstoreCap(PavSurf), notUsedI)
 endif
 if (soilmoist(ConifSurf)>soilstoreCap(ConifSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of conif soil store.',&
                   soilmoist(ConifSurf), soilstoreCap(ConifSurf), notUsedI)
 endif
 if (soilmoist(DecidSurf)>soilstoreCap(DecidSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of deciduous soil store.',&
                   soilmoist(DecidSurf), soilstoreCap(DecidSurf), notUsedI)
 endif
 if (soilmoist(GrassUSurf)>soilstoreCap(GrassUSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of unirrigated soil store.',&
                   soilmoist(GrassUSurf), soilstoreCap(GrassUSurf), notUsedI)
 endif
 if (soilmoist(GrassISurf)>soilstoreCap(GrassISurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of irrigated soil store.',&
                   soilmoist(GrassISurf), soilstoreCap(GrassISurf), notUsedI)
 endif

 !Snow stuff
 if (snowUse==1) then
    if (SnowWaterBldgState>CRWmax*SnowPackBldg) then
       call ErrorHint(36,'InitialCond: SnowWaterBldgState', SnowWaterBldgState, SnowPackBldg, notUsedI)
    endif 
    if (SnowWaterPavstate>CRWmax*SnowPackPav) then
       call ErrorHint(36,'InitialCond: SnowWaterPavState', SnowWaterPavstate, SnowPackPav, notUsedI)
    endif 
    if (SnowWaterConifstate>CRWmax*SnowPackConif) then
       call ErrorHint(36,'InitialCond: SnowWaterConifstate', SnowWaterConifstate, SnowPackConif, notUsedI)
    endif 
    if (SnowWaterDecState>CRWmax*SnowPackDec) then
       call ErrorHint(36,'InitialCond: SnowWaterDecState', SnowWaterDecState, SnowPackDec, notUsedI)
    endif 
    if (SnowWaterGrassIrrState>CRWmax*SnowPackGrassIrr) then
       call ErrorHint(36,'InitialCond: SnowWaterGrassIrrState', SnowWaterGrassIrrState, SnowPackGrassIrr, notUsedI)
    endif 
    if (SnowWaterGrassUnirState>CRWmax*SnowPackGrassUnir) then
       call ErrorHint(36,'InitialCond: SnowWaterGrassUnirState', SnowWaterGrassUnirState, SnowPackGrassUnir, notUsedI)
    endif 
 endif

!SnowWaterWaterstate,& ??
!SnowPackWater,& ??

 





 end subroutine CheckInitial  

  