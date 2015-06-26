!=========================================================================
! sg feb 2012
! only run once at start - fixed for all grids and all years
 
 SUBROUTINE OverallRunControl
! Last modified by HCW 06 Mar 2015 - Removed options 10,20,30 (NARPOutput) for NetRadiationChoice
! Last modified by HCW 06 Feb 2015
!  File ID numbers changed so they are unique
! Last modified by HCW 19 Dec 2014
! To Do: 
! 	- Holidays.txt input file needs to be read in and coded into model
!	- Add column header checks for SiteSelect
!-------------------------------------------------------------------------
  
  use allocateArray 
  use ColNamesInputFiles
  use data_in      
  use defaultNotUsed
  use FileName
  use initial
  use gis_data     
  use mod_z	   
  use resist       
  use snowMod
  use sues_data    
  use time
   
  IMPLICIT NONE
  
  integer:: iv,i,ii,SkipCounter            !iv and i, ii are integers used in do loops
  character(len=50):: FileN
  
  ! ---- Namelist for RunControl.nml ----
  namelist/RunControl/AnthropHeatChoice,&
        CBLuse,& !s.o.
        gsChoice,&
        NetRadiationChoice,&
        RoughLen_heat,&
        QSChoice,&
        OHMIncQF,&
        smd_choice,&							
        StabilityMethod,&
        WU_choice,&
        z0_method,&  
        FileCode,&
        FileInputPath,&
        FileOutputPath,&
        SkipHeaderSiteInfo,&
        SkipHeaderMet,&
        MultipleMetFiles,&
        KeepTstepFilesIn,&
        KeepTstepFilesOut,&
        WriteSurfsFile,&
        SnowFractionChoice,&
        SNOWuse,&
        SOLWEIGuse,&
        TIMEZONE,& 
        Tstep,& 
        Z
  ! -------------------------------------
        
  FileCode='none'
  !smithFile='Smith1966.grd'
   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Read in the RunControl.nml file
  open(55,File='RunControl.nml',err=200,status='old') !Change with needs
  read(55,nml=RunControl,err=201)
  close(55)
  
  !Check for problems with FileCode
  if (FileCode=='none') call ErrorHint(26,trim("RunControl.nml FileCode is missing"),notUsed,notUsed,notUsedI)
   
  !-----------------------------------------------------------------------
 
  !Write RunControl information to FileChoices.txt
  FileChoices=trim(FileOutputPath)//trim(FileCode)//'_FileChoices.txt'
  open (12,file=FileChoices,err=203)
  write(12,nml=RunControl)
  close(12)
  
  !Determine what should be done with respect to radiation
  AlbedoChoice=0
  ldown_option=0
  if(netRadiationChoice==0)then     !Observed Q* from the met input file will be used
    if(snowUse==1) then             !If snow is modelled, NARP is needed for surface temperature
       netRadiationChoice=3000
       ldown_option=3   	    !Ldown will be modelled
       !NetRadiationChoice=NetRadiationChoice/1000
    endif
    
  elseif(netRadiationChoice>0)then  !Modelled Q* is used (NARP)
       AlbedoChoice=-9
       if(NetRadiationChoice<10) then
            AlbedoChoice=0
            if(NetRadiationChoice==1)ldown_option=1
            if(NetRadiationChoice==2)ldown_option=2
            if(NetRadiationChoice==3)ldown_option=3
                   
       elseif(NetRadiationChoice>=100.and.NetRadiationChoice<1000) then
            AlbedoChoice=1  
            if(NetRadiationChoice==100)ldown_option=1
            if(NetRadiationChoice==200)ldown_option=2
            if(NetRadiationChoice==300)ldown_option=3
               NetRadiationChoice=NetRadiationChoice/100
       endif
       
       !If bad NetRadiationChoice value
	   if(netRadiationChoice>3.or. AlbedoChoice==-9)then
           write(*,*) 'NetRadiationChoice=',NetRadiationChoice
           write(*,*) 'Value not usable'
           stop
       endif
  endif
  
  !------------------------------------------------------------------
  !Print run information on the screen
  write(*,*)'--------------------------------------------------------'
  write(*,*)"LUMPS/Suews 2012a - relevant references"
  write(*,*)"LUMPS - Grimmond and Oke (2002) JAM, 41, 79-810"
  write(*,*)"OHM - Grimmond and Oke (1999) JAM, 38, 922-940"
  write(*,*)"NARP - Offerle et al. (2003) JAM"
  write(*,*)"SUES - Evaporation Grimmond & Oke (1991) WRR"
  write(*,*)"Water Balance Model Grimmond et al. (1986) WRR"
  write(*,*)"NARP - Long wave improvements (Loridan et al. 2011 JAMC)"
  write(*,*)"SUEWS - Anthropogenic heat, etc (Jarvi et al. 2011 JH)"
  write(*,*)"SUEWS - Snow module included (Jarvi et al. 2014 GMD)"
  write(*,*)'--------------------------------------------------------'    
  

  !=======================================================================
  !======================== Read input files ============================= 	
  ! This part reads the input files derived from the SiteInfo spreadsheet
  
  write(*,*) 'Reading the following input files:'

  !=======================SUEWS_SiteSelect.txt============================
  FileN='SUEWS_SiteSelect.txt'
  call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesSiteSelect=nlines
  allocate(SiteSelect(nlinesSiteSelect,ncolumnsSiteSelect))
  !Read input file 
  open(21,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
  do SkipCounter=1,(SkipHeaderSiteInfo-1)
     read(21,*) 	!Skip lines before header
  enddo
  read(21,*) (HeaderSiteSelect_File(iv),iv=1,ncolumnsSiteSelect) !Get header
    
  do i=1,nlinesSiteSelect
     read(21,*) (SiteSelect(i,iv),iv=1,ncolumnsSiteSelect)
     !write(*,*) (SiteSelect(i,iv),iv=1,ncolumnsSiteSelect)
  enddo
  close(21)
  
  !call InputHeaderCheck(FileN) !! Need to add column checks for SiteSelect.txt
  
  !=======================SUEWS_NonVeg.txt============================
  FileN='SUEWS_NonVeg.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesNonVeg=nlines
    allocate(NonVeg_Coeff(nlinesNonVeg,ncolumnsNonVeg))
    !Read input file 
    open(22,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(22,*) 	!Skip lines before header
    enddo
    read(22,*) (HeaderNonVeg_File(iv),iv=1,ncolumnsNonVeg) !Get header
      
    do i=1,nlinesNonVeg
       read(22,*) (NonVeg_Coeff(i,iv),iv=1,ncolumnsNonVeg)
       !write(*,*) (NonVeg_Coeff(i,iv),iv=1,ncolumnsNonVeg)
    enddo
  close(22)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesNonVeg
     do ii=1,nlinesNonVeg
        if(NonVeg_Coeff(i,ci_Code)==NonVeg_Coeff(ii,ci_Code) .and. i/=ii) then
           write(*,*) 'Code',NonVeg_Coeff(i,ci_Code),'in Impervious.txt not unique!'
           call ErrorHint(60,FileN,NonVeg_Coeff(i,ci_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   
  
  !=======================SUEWS_Veg.txt==============================
  FileN='SUEWS_Veg.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesVeg=nlines
    allocate(Veg_Coeff(nlinesVeg,ncolumnsVeg))
    !Read input file 
    open(23,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(23,*) 	!Skip lines before header
    enddo
    read(23,*) (HeaderVeg_File(iv),iv=1,ncolumnsVeg) !Get header
      
    do i=1,nlinesVeg
       read(23,*) (Veg_Coeff(i,iv),iv=1,ncolumnsVeg)
       !write(*,*) (Veg_Coeff(i,iv),iv=1,ncolumnsVeg)
    enddo
  close(23)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesVeg
     do ii=1,nlinesVeg
        if(Veg_Coeff(i,cp_Code)==Veg_Coeff(ii,cp_Code) .and. i/=ii) then
           write(*,*) 'Code',Veg_Coeff(i,cp_Code),'in Pervious.txt not unique!'
           call ErrorHint(60,FileN,Veg_Coeff(i,cp_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   

  !=======================SUEWS_Water.txt=================================
  FileN='SUEWS_Water.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesWater=nlines
    allocate(Water_Coeff(nlinesWater,ncolumnsWater))
    !Read input file 
    open(24,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(24,*) 	!Skip lines before header
    enddo
    read(24,*) (HeaderWater_File(iv),iv=1,ncolumnsWater) !Get header
      
    do i=1,nlinesWater
       read(24,*) (Water_Coeff(i,iv),iv=1,ncolumnsWater)
       !write(*,*) (Water_Coeff(i,iv),iv=1,ncolumnsWater)
    enddo
  close(24)
  
  call InputHeaderCheck(FileN)
  
  ! Check codes are unique
  do i=1,nlinesWater
     do ii=1,nlinesWater
        if(Water_Coeff(i,cw_Code)==Water_Coeff(ii,cw_Code) .and. i/=ii) then
           write(*,*) 'Code',Water_Coeff(i,cw_Code),'in Water.txt not unique!'
           call ErrorHint(60,FileN,Water_Coeff(i,cw_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   


  !=======================SUEWS_Snow.txt==================================
  FileN='SUEWS_Snow.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesSnow=nlines
    allocate(Snow_Coeff(nlinesSnow,ncolumnsSnow))
    !Read input file 
    open(25,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(25,*) 	!Skip lines before header
    enddo
    read(25,*) (HeaderSnow_File(iv),iv=1,ncolumnsSnow) !Get header
      
    do i=1,nlinesSnow
       read(25,*) (Snow_Coeff(i,iv),iv=1,ncolumnsSnow)
       !write(*,*) (Snow_Coeff(i,iv),iv=1,ncolumnsSnow)
    enddo
  close(25)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesSnow
     do ii=1,nlinesSnow
        if(Snow_Coeff(i,cs_Code)==Snow_Coeff(ii,cs_Code) .and. i/=ii) then
           write(*,*) 'Code',Snow_Coeff(i,cs_Code),'in Snow.txt not unique!'
           call ErrorHint(60,FileN,Snow_Coeff(i,cs_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   

  
  !=======================SUEWS_Soil.txt==================================
  FileN='SUEWS_Soil.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesSoil=nlines
    allocate(Soil_Coeff(nlinesSoil,ncolumnsSoil))
    !Read input file 
    open(26,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(26,*) 	!Skip lines before header
    enddo
    read(26,*) (HeaderSoil_File(iv),iv=1,ncolumnsSoil) !Get header
      
    do i=1,nlinesSoil
       read(26,*) (Soil_Coeff(i,iv),iv=1,ncolumnsSoil)
       !write(*,*) (Soil_Coeff(i,iv),iv=1,ncolumnsSoil)
    enddo
  close(26)
  
  call InputHeaderCheck(FileN)
  
  ! Check codes are unique
  do i=1,nlinesSoil
     do ii=1,nlinesSoil
        if(Soil_Coeff(i,cSo_Code)==Soil_Coeff(ii,cSo_Code) .and. i/=ii) then
           write(*,*) 'Code',Soil_Coeff(i,cSo_Code),'in Soil.txt not unique!'
           call ErrorHint(60,FileN,Soil_Coeff(i,cSo_Code),notUsed,notUsedI)
        endif
     enddo
  enddo     

  !===================SUEWS_Conductance.txt===============================
  FileN='SUEWS_Conductance.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesConductance=nlines
    allocate(Conductance_Coeff(nlinesConductance,ncolumnsConductance))
    !Read input file 
    open(27,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(27,*) 	!Skip lines before header
    enddo
    read(27,*) (HeaderCond_File(iv),iv=1,ncolumnsConductance) !Get header
      
    do i=1,nlinesConductance
       read(27,*) (Conductance_Coeff(i,iv),iv=1,ncolumnsConductance)
       !write(*,*) (Conductance_Coeff(i,iv),iv=1,ncolumnsConductance)
    enddo
  close(27)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesConductance
     do ii=1,nlinesConductance
        if(Conductance_Coeff(i,cc_Code)==Conductance_Coeff(ii,cc_Code) .and. i/=ii) then
           write(*,*) 'Code',Conductance_Coeff(i,cc_Code),'in Conductance.txt not unique!'
           call ErrorHint(60,FileN,Conductance_Coeff(i,cc_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   
  
  !===================SUEWS_OHMCoefficients.txt===========================
  FileN='SUEWS_OHMCoefficients.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesOHMCoefficients=nlines
    allocate(OHMCoefficients_Coeff(nlinesOHMCoefficients,ncolumnsOHMCoefficients))
    !Read input file 
    open(28,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(28,*) 	!Skip lines before header
    enddo
    read(28,*) (HeaderOHMCoefficients_File(iv),iv=1,ncolumnsOHMCoefficients) !Get header
  
    do i=1,nlinesOHMCoefficients
       read(28,*) (OHMCoefficients_Coeff(i,iv),iv=1,ncolumnsOHMCoefficients)
       !write(*,*) (OHMCoefficients_Coeff(i,iv),iv=1,ncolumnsOHMCoefficients)
    enddo
  close(28)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesOHMCoefficients
     do ii=1,nlinesOHMCoefficients
        if(OHMCoefficients_Coeff(i,cO_Code)==OHMCoefficients_Coeff(ii,cO_Code) .and. i/=ii) then
           write(*,*) 'Code',OHMCoefficients_Coeff(i,cO_Code),'in OHMCoefficients.txt not unique!'
           call ErrorHint(60,FileN,OHMCoefficients_Coeff(i,cO_Code),notUsed,notUsedI)
        endif
     enddo
  enddo     

  !================SUEWS_AnthropogenicHeat.txt============================
  FileN='SUEWS_AnthropogenicHeat.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesAnthropogenicHeat=nlines
    allocate(AnthropogenicHeat_Coeff(nlinesAnthropogenicHeat,ncolumnsAnthropogenicHeat))
    !Read input file 
    open(29,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(29,*) 	!Skip lines before header
    enddo
    read(29,*) (HeaderAnthropogenicHeat_File(iv),iv=1,ncolumnsAnthropogenicHeat) !Get header
      
    do i=1,nlinesAnthropogenicHeat
       read(29,*) (AnthropogenicHeat_Coeff(i,iv),iv=1,ncolumnsAnthropogenicHeat)
       !write(*,*) (AnthropogenicHeat_Coeff(i,iv),iv=1,ncolumnsAnthropogenicHeat)
    enddo
  close(29)

  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesAnthropogenicHeat
     do ii=1,nlinesAnthropogenicHeat
        if(AnthropogenicHeat_Coeff(i,cA_Code)==AnthropogenicHeat_Coeff(ii,cA_Code) .and. i/=ii) then
           write(*,*) 'Code',AnthropogenicHeat_Coeff(i,cA_Code),'in AnthropogenicHeat.txt not unique!'
           call ErrorHint(60,FileN,AnthropogenicHeat_Coeff(i,cA_Code),notUsed,notUsedI)
        endif
     enddo
  enddo     

  !================SUEWS_Irrigation.txt===================================
  FileN='SUEWS_Irrigation.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesIrrigation=nlines
    allocate(Irrigation_Coeff(nlinesIrrigation,ncolumnsIrrigation))
    !Read input file 
    open(30,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(30,*) 	!Skip lines before header
    enddo
    read(30,*) (HeaderIrrigation_File(iv),iv=1,ncolumnsIrrigation) !Get header
      
    do i=1,nlinesIrrigation
       read(30,*) (Irrigation_Coeff(i,iv),iv=1,ncolumnsIrrigation)
       !write(*,*) (Irrigation_Coeff(i,iv),iv=1,ncolumnsIrrigation)
    enddo
  close(30)        
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesIrrigation
     do ii=1,nlinesIrrigation
        if(Irrigation_Coeff(i,cIr_Code)==Irrigation_Coeff(ii,cIr_Code) .and. i/=ii) then
           write(*,*) 'Code',Irrigation_Coeff(i,cIr_Code),'in Irrigation.txt not unique!'
           call ErrorHint(60,FileN,Irrigation_Coeff(i,cIr_Code),notUsed,notUsedI)
        endif
     enddo
  enddo     
  
  !================SUEWS_Profiles.txt=====================================
  FileN='SUEWS_Profiles.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesProfiles=nlines
    allocate(Profiles_Coeff(nlinesProfiles,ncolumnsProfiles))
    !Read input file 
    open(31,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(31,*) 	!Skip lines before header
    enddo
    read(31,*) (HeaderProfiles_File(iv),iv=1,ncolumnsProfiles) !Get header
      
    do i=1,nlinesProfiles
       read(31,*) (Profiles_Coeff(i,iv),iv=1,ncolumnsProfiles)
       !write(*,*) (Profiles_Coeff(i,iv),iv=1,ncolumnsProfiles)
    enddo
  close(31)
  
  call InputHeaderCheck(FileN)

  ! Check codes are unique
  do i=1,nlinesProfiles
     do ii=1,nlinesProfiles
        if(Profiles_Coeff(i,cPr_Code)==Profiles_Coeff(ii,cPr_Code) .and. i/=ii) then
           write(*,*) 'Code',Profiles_Coeff(i,cPr_Code),'in Profiles.txt not unique!'
           call ErrorHint(60,FileN,Profiles_Coeff(i,cPr_Code),notUsed,notUsedI)
        endif
     enddo
  enddo   

 !================SUEWS_WithinGridWaterDist.txt===========================
  FileN='SUEWS_WithinGridWaterDist.txt'
    call NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
    nlinesWGWaterDist=nlines
    allocate(WGWaterDist_Coeff(nlinesWGWaterDist,ncolumnsWGWaterDist))
    !Read input file 
    open(32,file=trim(FileInputPath)//trim(FileN),err=300,status='old')
    do SkipCounter=1,(SkipHeaderSiteInfo-1)
       read(32,*) 	!Skip lines before header
    enddo
    read(32,*) (HeaderWGWaterDist_File(iv),iv=1,ncolumnsWGWaterDist) !Get header
      
    do i=1,nlinesWGWaterDist
       read(32,*) (WGWaterDist_Coeff(i,iv),iv=1,ncolumnsWGWaterDist)
       !write(*,*) (WGWaterDist_Coeff(i,iv),iv=1,ncolumnsWGWaterDist)
    enddo
  close(32)  
  
  call InputHeaderCheck(FileN)
  
  ! Check codes are unique
  do i=1,nlinesWGWaterDist
     do ii=1,nlinesWGWaterDist
        if(WGWaterDist_Coeff(i,cWG_Code)==WGWaterDist_Coeff(ii,cWG_Code) .and. i/=ii) then
           write(*,*) 'Code',WGWaterDist_Coeff(i,cWG_Code),'in WithinGridWaterDist.txt not unique!'
           call ErrorHint(60,FileN,WGWaterDist_Coeff(i,cWG_Code),notUsed,notUsedI)
        endif
     enddo
  enddo    
    
  !=======================================================================
  !=======================================================================
    
  !-----------------------------------------------------------------------
  !SUEWS run information
  InputMetFormat=10	!Input met data file in LUMPS format(1) or SUEWS format(10)
  LAICalcYes=1		!Use observed(0) or modelled(1) LAI
  ity=2			!Evaporation calculated according to Rutter(1) or Shuttleworth(2)
  WriteDailyState = 1   !Daily state file written
  tstepcount=0
   
  t_INTERVAL = 3600   !Number of seconds in an hour  
   
  !Calculate nsh (number of steps per hour) from model timestep (tstep) set in in RunControl
  nsh_real = t_INTERVAL/real(tstep,kind(1d0))
  
  ! Check nsh is an integer	
  if(nsh_real==int(nsh_real)) then
     nsh = int(nsh_real)   
  else
     call ErrorHint(39,'TSTEP must divide into t_INTERVAL exactly.',real(tstep,kind(1d0)),real(t_INTERVAL,kind(1d0)),notUsedI)
  endif 
  
  ! Check nsh is reasonable
  if(nsh_real<6.or.nsh_real>60) then   
     call ErrorHint(39,'TSTEP is too small or too large.',real(tstep,kind(1d0)),real(t_INTERVAL,kind(1d0)),notUsedI)        
  endif
  
  ! Cast integer nsh as nsh_real for use in calculations
  nsh_real = real(nsh,kind(1d0))
  ! Cast integer tstep as tstep_real for use in calculations
  tstep_real = real(tstep,kind(1d0))
    
  !! Check this is still valid for v2015a  
  HalfTimeStep=real(tstep_real)/2/(24*3600)   !Used in sun_position to get sunpos in the middle of timestep
  
  return

  !-------Possible problems-----------------------------------------------
  200 call ErrorHint(47,'RunControl.nml',notUsed,notUsed,notUsedI)
  201 call ErrorHint(48,'RunControl.nml',notUsed,notUsed,notUsedI)
 
  203 call ErrorHint(47,trim(FileChoices),notUsed,notUsed,notUsedI)
 
  300 call ErrorHint(48,trim(FileN),notUsed,notUsed,notUsedI)
  !-----------------------------------------------------------------------
   
  !pause
  
  END SUBROUTINE OverallRunControl
!=========================================================================



!=========================================================================
!This subroutine finds the number of rows in each input file
!INPUT: FileN      		Name of the input file 
!       SkipHeaderLines 	Number of header rows to skip
!Made by LJ/HCW in Oct 2014
 
 subroutine NumberRows(FileN,SkipHeaderLines)

   use data_in
   use DefaultNotUsed
   use Initial

   IMPLICIT NONE

   character(len=50):: FileN
   integer:: SkipHeaderLines, RunNumber
   integer:: SkipCounter
      
   write(*,*) FileN
   open(39,file=trim(FileInputPath)//trim(FileN),err=204,status='old')

   if(SkipHeaderLines > 0) then
     do SkipCounter=1,SkipHeaderLines
        read(39,*,err=205) 
        !write(*,*) SkipCounter, SkipHeaderLines
     enddo
   endif

   nlines = 0 !Initialize nlines
   do
      read(39,*) RunNumber
       if(RunNumber==-9) exit
       nlines = nlines + 1
   end do
   !write(*,*) 'nlines read: ',nlines
   close(39)
   
  return

204 call ErrorHint(47,trim(FileInputPath)//trim(FileN),notUsed,notUsed,notUsedI)
205 call ErrorHint(48,trim(FileInputPath)//trim(FileN),notUsed,notUsed,notUsedI)

 end subroutine NumberRows
!=========================================================================



!----------------------------------------------------------------------------------------------
!Moves input information corresponding to row rr of SiteSelect into row Gridiv of SurfaceChar
! - currently done once per year for each grid
!Made by HW&LJ Oct 2014.
!Last modified: HCW 26 Jan 2015
! Interpolated hourly energy use profiles to resolution of model timestep.
! Interpolated hourly water use profiles to resolution of model timestep.
! Normalised energy use & water use profiles as required.
! Hourly snow clearing profile is a 0-1 switch, and thus is not interpolated nor normalised.
! To Do 
!	- Rename profiles 1-24 rather than 0-23?
!----------------------------------------------------------------------------------------------
 SUBROUTINE InitializeSurfaceCharacteristics(Gridiv,rr)

   use allocateArray
   use ColNamesInputFiles
   use data_in
   use defaultNotUsed
   use Initial
   use sues_data
   
   IMPLICIT NONE

   integer:: Gridiv,&    !Row of SurfaceChar where input information will be stored
             rr          !Row of SiteSelect that matches current grid and year
                   
   !-------------------------------------------------------------------------------------------
   
   ! Initialise row of SurfaceChar
   SurfaceChar(Gridiv,:) = -999
         	
   ! Transfer data in SiteSelect to SurfaceChar
   SurfaceChar(Gridiv,1:ncolumnsSiteSelect) = SiteSelect(rr,1:ncolumnsSiteSelect) !Cols in same order as in SiteSelect.txt
     
   ! ======== Retrieve information from other input files via codes ========
          
   ! ---- Find code for Paved surface (Impervious) ----
   call CodeMatchNonVeg(rr,c_PavedCode)
   ! Transfer characteristics to SurfaceChar for Paved surface
   SurfaceChar(gridiv,c_Alb(PavSurf)) 	       = NonVeg_Coeff(iv5,ci_Alb)
   SurfaceChar(gridiv,c_Emis(PavSurf))         = NonVeg_Coeff(iv5,ci_Emis)
   SurfaceChar(gridiv,c_StorMin(PavSurf))      = NonVeg_Coeff(iv5,ci_StorMin)
   SurfaceChar(gridiv,c_StorMax(PavSurf))      = NonVeg_Coeff(iv5,ci_StorMax)
   SurfaceChar(gridiv,c_WetThresh(PavSurf))   = NonVeg_Coeff(iv5,ci_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(PavSurf))   = NonVeg_Coeff(iv5,ci_StateLimit)
   SurfaceChar(gridiv,c_DrEq(PavSurf)) 	       = NonVeg_Coeff(iv5,ci_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(PavSurf))      = NonVeg_Coeff(iv5,ci_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(PavSurf))      = NonVeg_Coeff(iv5,ci_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(PavSurf))    = NonVeg_Coeff(iv5,ci_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(PavSurf))   = NonVeg_Coeff(iv5,ci_SnowLimPat)
   SurfaceChar(gridiv,c_SnowLimRem(PavSurf))   = NonVeg_Coeff(iv5,ci_SnowLimRem)
   SurfaceChar(gridiv,c_OHMCode_SWet(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WDry)
   ! Use SoilCode for Paved to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(PavSurf))
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(PavSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(PavSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(PavSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(PavSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(PavSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(PavSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(PavSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(PavSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac) 
   
   ! Get OHM characteristics for Paved
   call CodeMatchOHM(Gridiv,PavSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,PavSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,PavSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,PavSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   
   ! Get water distribution (within grid) for Paved
   call CodeMatchDist(rr,c_WGPavedCode,cWG_ToPaved)
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(PavSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(PavSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)   
       
   ! ---- Find code for Bldgs surface (Impervious) ----
   call CodeMatchNonVeg(rr,c_BldgsCode) 
   ! Transfer characteristics to SurfaceChar for Bldgs surface
   SurfaceChar(gridiv,c_Alb(BldgSurf))          = NonVeg_Coeff(iv5,ci_Alb)
   SurfaceChar(gridiv,c_Emis(BldgSurf))         = NonVeg_Coeff(iv5,ci_Emis)
   SurfaceChar(gridiv,c_StorMin(BldgSurf))      = NonVeg_Coeff(iv5,ci_StorMin)
   SurfaceChar(gridiv,c_StorMax(BldgSurf))      = NonVeg_Coeff(iv5,ci_StorMax)
   SurfaceChar(gridiv,c_WetThresh(BldgSurf))   = NonVeg_Coeff(iv5,ci_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(BldgSurf))   = NonVeg_Coeff(iv5,ci_StateLimit)
   SurfaceChar(gridiv,c_DrEq(BldgSurf))         = NonVeg_Coeff(iv5,ci_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(BldgSurf))      = NonVeg_Coeff(iv5,ci_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(BldgSurf))      = NonVeg_Coeff(iv5,ci_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(BldgSurf))    = NonVeg_Coeff(iv5,ci_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(BldgSurf))   = NonVeg_Coeff(iv5,ci_SnowLimPat)
   SurfaceChar(gridiv,c_SnowLimRem(BldgSurf))   = NonVeg_Coeff(iv5,ci_SnowLimRem)
   SurfaceChar(gridiv,c_OHMCode_SWet(BldgSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(BldgSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(BldgSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(BldgSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WDry)   
   ! Use SoilCode for Bldgs to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(BldgSurf))
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(BldgSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(BldgSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(BldgSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(BldgSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(BldgSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(BldgSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(BldgSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(BldgSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac)
   !Get OHM characteristics for Bldgs
   call CodeMatchOHM(Gridiv,BldgSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,BldgSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,BldgSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,BldgSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)      
   
   ! Get water distribution (within grid) for Bldgs
   call CodeMatchDist(rr,c_WGBldgsCode,cWG_ToBldgs)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(BldgSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(BldgSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)   
 
   ! ---- Find code for EveTr surface (Pervious) ----
   call CodeMatchVeg(rr,c_EveTrCode)       
   ! Transfer characteristics to SurfaceChar for EveTr surface
   ! All surfaces (1-nsurf)
   SurfaceChar(gridiv,c_Alb(ConifSurf))        = Veg_Coeff(iv5,cp_Alb)
   SurfaceChar(gridiv,c_Emis(ConifSurf))       = Veg_Coeff(iv5,cp_Emis)
   SurfaceChar(gridiv,c_StorMin(ConifSurf))    = Veg_Coeff(iv5,cp_StorMin)
   SurfaceChar(gridiv,c_StorMax(ConifSurf))    = Veg_Coeff(iv5,cp_StorMax)
   SurfaceChar(gridiv,c_WetThresh(ConifSurf)) = Veg_Coeff(iv5,cp_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(ConifSurf)) = Veg_Coeff(iv5,cp_StateLimit)
   SurfaceChar(gridiv,c_DrEq(ConifSurf))       = Veg_Coeff(iv5,cp_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(ConifSurf))    = Veg_Coeff(iv5,cp_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(ConifSurf))    = Veg_Coeff(iv5,cp_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(ConifSurf))  = Veg_Coeff(iv5,cp_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(ConifSurf)) = Veg_Coeff(iv5,cp_SnowLimPat)
   ! Veg surfaces only (1-nvegsurf)
   SurfaceChar(gridiv,c_BaseT(ivConif))        = Veg_Coeff(iv5,cp_BaseT)
   SurfaceChar(gridiv,c_BaseTe(ivConif))       = Veg_Coeff(iv5,cp_BaseTe)
   SurfaceChar(gridiv,c_GDDFull(ivConif))      = Veg_Coeff(iv5,cp_GDDFull)
   SurfaceChar(gridiv,c_SDDFull(ivConif))      = Veg_Coeff(iv5,cp_SDDFull)
   SurfaceChar(gridiv,c_LAIMin(ivConif))       = Veg_Coeff(iv5,cp_LAIMin)
   SurfaceChar(gridiv,c_LAIMax(ivConif))       = Veg_Coeff(iv5,cp_LAIMax)
   SurfaceChar(gridiv,c_GsMax(ivConif))        = Veg_Coeff(iv5,cp_GsMax)
   SurfaceChar(gridiv,c_LAIEq(ivConif))        = Veg_Coeff(iv5,cp_LAIEq)
   SurfaceChar(gridiv,c_LeafGP1(ivConif))      = Veg_Coeff(iv5,cp_LeafGP1)
   SurfaceChar(gridiv,c_LeafGP2(ivConif))      = Veg_Coeff(iv5,cp_LeafGP2)
   SurfaceChar(gridiv,c_LeafOP1(ivConif))      = Veg_Coeff(iv5,cp_LeafOP1)
   SurfaceChar(gridiv,c_LeafOP2(ivConif))      = Veg_Coeff(iv5,cp_LeafOP2)
   ! OHM codes
   SurfaceChar(gridiv,c_OHMCode_SWet(ConifSurf)) = Veg_Coeff(iv5,cp_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(ConifSurf)) = Veg_Coeff(iv5,cp_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(ConifSurf)) = Veg_Coeff(iv5,cp_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(ConifSurf)) = Veg_Coeff(iv5,cp_OHMCode_WDry)     
   ! Use SoilCode for EveTr to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(ConifSurf))     
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(ConifSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(ConifSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(ConifSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(ConifSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(ConifSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(ConifSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(ConifSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(ConifSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac)
   !Get OHM characteristics for Conif
   call CodeMatchOHM(Gridiv,ConifSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,ConifSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,ConifSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,ConifSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)      
   
   ! Get water distribution (within grid) for EveTr
   call CodeMatchDist(rr,c_WGEveTrCode,cWG_ToEveTr)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(ConifSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(ConifSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)     
     
   ! ---- Find code for DecTr surface (Pervious) ----
   call CodeMatchVeg(rr,c_DecTrCode)        
   ! Transfer characteristics to SurfaceChar for DecTr surface
   ! All surfaces (1-nsurf)
   SurfaceChar(gridiv,c_Alb(DecidSurf))        = Veg_Coeff(iv5,cp_Alb)
   SurfaceChar(gridiv,c_Emis(DecidSurf))       = Veg_Coeff(iv5,cp_Emis)
   SurfaceChar(gridiv,c_StorMin(DecidSurf))    = Veg_Coeff(iv5,cp_StorMin)
   SurfaceChar(gridiv,c_StorMax(DecidSurf))    = Veg_Coeff(iv5,cp_StorMax)
   SurfaceChar(gridiv,c_WetThresh(DecidSurf)) = Veg_Coeff(iv5,cp_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(DecidSurf)) = Veg_Coeff(iv5,cp_StateLimit)
   SurfaceChar(gridiv,c_DrEq(DecidSurf))       = Veg_Coeff(iv5,cp_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(DecidSurf))    = Veg_Coeff(iv5,cp_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(DecidSurf))    = Veg_Coeff(iv5,cp_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(DecidSurf))  = Veg_Coeff(iv5,cp_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(DecidSurf)) = Veg_Coeff(iv5,cp_SnowLimPat)
   ! Veg surfaces only (1-nvegsurf)
   SurfaceChar(gridiv,c_BaseT(ivDecid))        = Veg_Coeff(iv5,cp_BaseT)
   SurfaceChar(gridiv,c_BaseTe(ivDecid))       = Veg_Coeff(iv5,cp_BaseTe)
   SurfaceChar(gridiv,c_GDDFull(ivDecid))      = Veg_Coeff(iv5,cp_GDDFull)
   SurfaceChar(gridiv,c_SDDFull(ivDecid))      = Veg_Coeff(iv5,cp_SDDFull)
   SurfaceChar(gridiv,c_LAIMin(ivDecid))       = Veg_Coeff(iv5,cp_LAIMin)
   SurfaceChar(gridiv,c_LAIMax(ivDecid))       = Veg_Coeff(iv5,cp_LAIMax)
   SurfaceChar(gridiv,c_GsMax(ivDecid))        = Veg_Coeff(iv5,cp_GsMax)
   SurfaceChar(gridiv,c_LAIEq(ivDecid))        = Veg_Coeff(iv5,cp_LAIEq)
   SurfaceChar(gridiv,c_LeafGP1(ivDecid))      = Veg_Coeff(iv5,cp_LeafGP1)
   SurfaceChar(gridiv,c_LeafGP2(ivDecid))      = Veg_Coeff(iv5,cp_LeafGP2)
   SurfaceChar(gridiv,c_LeafOP1(ivDecid))      = Veg_Coeff(iv5,cp_LeafOP1)
   SurfaceChar(gridiv,c_LeafOP2(ivDecid))      = Veg_Coeff(iv5,cp_LeafOP2)
   ! OHM codes
   SurfaceChar(gridiv,c_OHMCode_SWet(DecidSurf)) = Veg_Coeff(iv5,cp_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(DecidSurf)) = Veg_Coeff(iv5,cp_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(DecidSurf)) = Veg_Coeff(iv5,cp_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(DecidSurf)) = Veg_Coeff(iv5,cp_OHMCode_WDry)      
   ! Use SoilCode for DecTr to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(DecidSurf))           
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(DecidSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(DecidSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(DecidSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(DecidSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(DecidSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(DecidSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(DecidSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(DecidSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac)   
   !Get OHM characteristics for Decid
   call CodeMatchOHM(Gridiv,DecidSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,DecidSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,DecidSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,DecidSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)         
  
  ! Get water distribution (within grid) for DecTr
   call CodeMatchDist(rr,c_WGDecTrCode,cWG_ToDecTr)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(DecidSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(DecidSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)       
   
   ! ---- Find code for Grass surface (Pervious) ----
   call CodeMatchVeg(rr,c_GrassCode)       
   ! Transfer characteristics to SurfaceChar for Grass surface
   ! All surfaces (1-nsurf)
   SurfaceChar(gridiv,c_Alb(GrassSurf))        = Veg_Coeff(iv5,cp_Alb)
   SurfaceChar(gridiv,c_Emis(GrassSurf))       = Veg_Coeff(iv5,cp_Emis)
   SurfaceChar(gridiv,c_StorMin(GrassSurf))    = Veg_Coeff(iv5,cp_StorMin)
   SurfaceChar(gridiv,c_StorMax(GrassSurf))    = Veg_Coeff(iv5,cp_StorMax)
   SurfaceChar(gridiv,c_WetThresh(GrassSurf)) = Veg_Coeff(iv5,cp_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(GrassSurf)) = Veg_Coeff(iv5,cp_StateLimit)
   SurfaceChar(gridiv,c_DrEq(GrassSurf))       = Veg_Coeff(iv5,cp_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(GrassSurf))    = Veg_Coeff(iv5,cp_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(GrassSurf))    = Veg_Coeff(iv5,cp_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(GrassSurf))  = Veg_Coeff(iv5,cp_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(GrassSurf)) = Veg_Coeff(iv5,cp_SnowLimPat)
   ! Veg surfaces only (1-nvegsurf)
   SurfaceChar(gridiv,c_BaseT(ivGrass))        = Veg_Coeff(iv5,cp_BaseT)
   SurfaceChar(gridiv,c_BaseTe(ivGrass))       = Veg_Coeff(iv5,cp_BaseTe)
   SurfaceChar(gridiv,c_GDDFull(ivGrass))      = Veg_Coeff(iv5,cp_GDDFull)
   SurfaceChar(gridiv,c_SDDFull(ivGrass))      = Veg_Coeff(iv5,cp_SDDFull)
   SurfaceChar(gridiv,c_LAIMin(ivGrass))       = Veg_Coeff(iv5,cp_LAIMin)
   SurfaceChar(gridiv,c_LAIMax(ivGrass))       = Veg_Coeff(iv5,cp_LAIMax)
   SurfaceChar(gridiv,c_GsMax(ivGrass))        = Veg_Coeff(iv5,cp_GsMax)
   SurfaceChar(gridiv,c_LAIEq(ivGrass))        = Veg_Coeff(iv5,cp_LAIEq)
   SurfaceChar(gridiv,c_LeafGP1(ivGrass))      = Veg_Coeff(iv5,cp_LeafGP1)
   SurfaceChar(gridiv,c_LeafGP2(ivGrass))      = Veg_Coeff(iv5,cp_LeafGP2)
   SurfaceChar(gridiv,c_LeafOP1(ivGrass))      = Veg_Coeff(iv5,cp_LeafOP1)
   SurfaceChar(gridiv,c_LeafOP2(ivGrass))      = Veg_Coeff(iv5,cp_LeafOP2)
   ! OHM codes
   SurfaceChar(gridiv,c_OHMCode_SWet(GrassSurf)) = Veg_Coeff(iv5,cp_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(GrassSurf)) = Veg_Coeff(iv5,cp_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(GrassSurf)) = Veg_Coeff(iv5,cp_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(GrassSurf)) = Veg_Coeff(iv5,cp_OHMCode_WDry)         
   ! Use SoilCode for Grass to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(GrassSurf))     
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(GrassSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(GrassSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(GrassSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(GrassSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(GrassSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(GrassSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(GrassSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(GrassSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac)   
   !Get OHM characteristics for Grass
   call CodeMatchOHM(Gridiv,GrassSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,GrassSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,GrassSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,GrassSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)      
   
  ! Get water distribution (within grid) for Grass
   call CodeMatchDist(rr,c_WGGrassCode,cWG_ToGrass)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(GrassSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(GrassSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)        
        
   ! ---- Find code for BSoil surface (Impervious) ----
   call CodeMatchNonVeg(rr,c_BSoilCode)       
   ! Transfer characteristics to SurfaceChar for BSoil surface
   ! All surfaces (1-nsurf)
   SurfaceChar(gridiv,c_Alb(BSoilSurf))        = NonVeg_Coeff(iv5,ci_Alb)
   SurfaceChar(gridiv,c_Emis(BSoilSurf))       = NonVeg_Coeff(iv5,ci_Emis)
   SurfaceChar(gridiv,c_StorMin(BSoilSurf))    = NonVeg_Coeff(iv5,ci_StorMin)
   SurfaceChar(gridiv,c_StorMax(BSoilSurf))    = NonVeg_Coeff(iv5,ci_StorMax)
   SurfaceChar(gridiv,c_WetThresh(BSoilSurf)) = NonVeg_Coeff(iv5,ci_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(BSoilSurf)) = NonVeg_Coeff(iv5,ci_StateLimit)
   SurfaceChar(gridiv,c_DrEq(BSoilSurf))       = NonVeg_Coeff(iv5,ci_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(BSoilSurf))    = NonVeg_Coeff(iv5,ci_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(BSoilSurf))    = NonVeg_Coeff(iv5,ci_DrCoef2)
   SurfaceChar(gridiv,c_SoilTCode(BSoilSurf))  = NonVeg_Coeff(iv5,ci_SoilTCode)
   SurfaceChar(gridiv,c_SnowLimPat(BSoilSurf)) = NonVeg_Coeff(iv5,ci_SnowLimPat)
   SurfaceChar(gridiv,c_SnowLimRem(BSoilSurf))   = NonVeg_Coeff(iv5,ci_SnowLimRem)
   SurfaceChar(gridiv,c_OHMCode_SWet(BSoilSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(BSoilSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(BSoilSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(BSoilSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WDry)         
   ! Use SoilCode for BSoil to find code for soil characteristics
   call CodeMatchSoil(Gridiv,c_SoilTCode(BSoilSurf))         
   ! Transfer soil characteristics to SurfaceChar
   SurfaceChar(gridiv,c_SoilDepth(BSoilSurf))    = Soil_Coeff(iv5,cSo_SoilDepth)
   SurfaceChar(gridiv,c_SoilStCap(BSoilSurf))    = Soil_Coeff(iv5,cSo_SoilStCap)
   SurfaceChar(gridiv,c_KSat(BSoilSurf))        = Soil_Coeff(iv5,cSo_KSat)
   SurfaceChar(gridiv,c_SoilDens(BSoilSurf))    = Soil_Coeff(iv5,cSo_SoilDens)
   SurfaceChar(gridiv,c_SoilInfRate(BSoilSurf)) = Soil_Coeff(iv5,cSo_SoilInfRate)
   SurfaceChar(gridiv,c_ObsSMDepth(BSoilSurf))  = Soil_Coeff(iv5,cSo_ObsSMDepth)
   SurfaceChar(gridiv,c_ObsSMMax(BSoilSurf))    = Soil_Coeff(iv5,cSo_ObsSMMax)
   SurfaceChar(gridiv,c_ObsSNRFrac(BSoilSurf))  = Soil_Coeff(iv5,cSo_ObsSNRFrac)   
   ! Get OHM characteristics for BSoil
   call CodeMatchOHM(Gridiv,BSoilSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,BSoilSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,BSoilSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,BSoilSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)      
   
   ! Get water distribution (within grid) for Bare soil
   call CodeMatchDist(rr,c_WGBSoilCode,cWG_ToBSoil)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(BSoilSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(BSoilSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)    
   

   ! ---- Find code for Water surface (Water) ----
   call CodeMatchWater(rr,c_WaterCode)              
   ! Transfer characteristics to SurfaceChar for Water surface
   ! All surfaces (1-nsurf)
   SurfaceChar(gridiv,c_Alb(WaterSurf))        = Water_Coeff(iv5,cw_Alb)
   SurfaceChar(gridiv,c_Emis(WaterSurf))       = Water_Coeff(iv5,cw_Emis)
   SurfaceChar(gridiv,c_StorMin(WaterSurf))    = Water_Coeff(iv5,cw_StorMin)
   SurfaceChar(gridiv,c_StorMax(WaterSurf))    = Water_Coeff(iv5,cw_StorMax)
   SurfaceChar(gridiv,c_WetThresh(WaterSurf)) = Water_Coeff(iv5,cw_WetThresh)
   SurfaceChar(gridiv,c_StateLimit(WaterSurf)) = Water_Coeff(iv5,cw_StateLimit)
   SurfaceChar(gridiv,c_DrEq(WaterSurf))       = Water_Coeff(iv5,cw_DrEq)
   SurfaceChar(gridiv,c_DrCoef1(WaterSurf))    = Water_Coeff(iv5,cw_DrCoef1)
   SurfaceChar(gridiv,c_DrCoef2(WaterSurf))    = Water_Coeff(iv5,cw_DrCoef2)
   ! OHM codes 
   SurfaceChar(gridiv,c_OHMCode_SWet(WaterSurf)) = Water_Coeff(iv5,cw_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(WaterSurf)) = Water_Coeff(iv5,cw_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(WaterSurf)) = Water_Coeff(iv5,cw_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(WaterSurf)) = Water_Coeff(iv5,cw_OHMCode_WDry)    
   ! Get OHM characteristics for Water   
   call CodeMatchOHM(Gridiv,WaterSurf,'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,WaterSurf,'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,WaterSurf,'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,WaterSurf,'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)         
   
   ! Get water distribution (within grid) for Water
   call CodeMatchDist(rr,c_WGWaterCode,cWG_ToWater)  
   ! Transfer distribution to SurfaceChar
   SurfaceChar(Gridiv,c_WGToPaved(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToPaved)   
   SurfaceChar(Gridiv,c_WGToBldgs(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBldgs)   
   SurfaceChar(Gridiv,c_WGToEveTr(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToEveTr)   
   SurfaceChar(Gridiv,c_WGToDecTr(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToDecTr)   
   SurfaceChar(Gridiv,c_WGToGrass(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToGrass)   
   SurfaceChar(Gridiv,c_WGToBSoil(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToBSoil)   
   SurfaceChar(Gridiv,c_WGToWater(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToWater)   
   SurfaceChar(Gridiv,c_WGToRunoff(WaterSurf))  = WGWaterDist_Coeff(iv5,cWG_ToRunoff)   
   SurfaceChar(Gridiv,c_WGToSoilStore(WaterSurf))    = WGWaterDist_Coeff(iv5,cWG_ToSoilStore)     
    

   ! ---- Find code for Snow surface (Snow) ----
   call CodeMatchSnow(rr,c_SnowCode)     
   ! Transfer characteristics to SurfaceChar for Snow surface
   SurfaceChar(gridiv,c_SnowRMFactor) = Snow_Coeff(iv5,cs_SnowRMFactor)
   SurfaceChar(gridiv,c_SnowTMFactor) = Snow_Coeff(iv5,cs_SnowTMFactor)
   SurfaceChar(gridiv,c_SnowAlbMin)   = Snow_Coeff(iv5,cs_SnowAlbMin)
   SurfaceChar(gridiv,c_SnowAlbMax)   = Snow_Coeff(iv5,cs_SnowAlbMax)
   SurfaceChar(gridiv,c_SnowAlb)      = Snow_Coeff(iv5,cs_SnowAlb)
   SurfaceChar(gridiv,c_SnowEmis)     = Snow_Coeff(iv5,cs_SnowEmis)
   SurfaceChar(gridiv,c_Snowtau_a)    = Snow_Coeff(iv5,cs_Snowtau_a)
   SurfaceChar(gridiv,c_Snowtau_f)    = Snow_Coeff(iv5,cs_Snowtau_f)
   SurfaceChar(gridiv,c_SnowPLimAlb)  = Snow_Coeff(iv5,cs_SnowPLimAlb)
   SurfaceChar(gridiv,c_SnowSDMin)    = Snow_Coeff(iv5,cs_SnowSDMin)
   SurfaceChar(gridiv,c_SnowSDMax)    = Snow_Coeff(iv5,cs_SnowSDMax)
   SurfaceChar(gridiv,c_Snowtau_r)    = Snow_Coeff(iv5,cs_Snowtau_r)
   SurfaceChar(gridiv,c_SnowCRWMin)   = Snow_Coeff(iv5,cs_SnowCRWMin)
   SurfaceChar(gridiv,c_SnowCRWMax)   = Snow_Coeff(iv5,cs_SnowCRWMax)
   SurfaceChar(gridiv,c_SnowPLimSnow) = Snow_Coeff(iv5,cs_SnowPLimSnow)
   SurfaceChar(gridiv,c_OHMCode_SWet(nsurf+1)) = Snow_Coeff(iv5,cs_OHMCode_SWet)
   SurfaceChar(gridiv,c_OHMCode_SDry(nsurf+1)) = Snow_Coeff(iv5,cs_OHMCode_SDry)
   SurfaceChar(gridiv,c_OHMCode_WWet(nsurf+1)) = Snow_Coeff(iv5,cs_OHMCode_WWet)
   SurfaceChar(gridiv,c_OHMCode_WDry(nsurf+1)) = Snow_Coeff(iv5,cs_OHMCode_WDry)   
   ! Get OHM characteristics for Snow
   call CodeMatchOHM(Gridiv,(nsurf+1),'SWet')  !Summer wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,(nsurf+1),'SDry')  !Summer dry	
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)   
   call CodeMatchOHM(Gridiv,(nsurf+1),'WWet')  !Winter wet
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)
   call CodeMatchOHM(Gridiv,(nsurf+1),'WDry')  !Winter dry
   ! Transfer OHM characteristics to SurfaceChar
   SurfaceChar(Gridiv,c_a1_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
   SurfaceChar(Gridiv,c_a2_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
   SurfaceChar(Gridiv,c_a3_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)         
         
   ! ---- Find code for Surface conductances ----
   call CodeMatchConductance(rr,c_CondCode)
   ! Transfer conductance characteristics to SurfaceChar
   SurfaceChar(gridiv,c_GsG1) 	    = Conductance_Coeff(iv5,cc_GsG1)
   SurfaceChar(gridiv,c_GsG2) 	    = Conductance_Coeff(iv5,cc_GsG2)
   SurfaceChar(gridiv,c_GsG3) 	    = Conductance_Coeff(iv5,cc_GsG3)
   SurfaceChar(gridiv,c_GsG4) 	    = Conductance_Coeff(iv5,cc_GsG4)
   SurfaceChar(gridiv,c_GsG5) 	    = Conductance_Coeff(iv5,cc_GsG5)
   SurfaceChar(gridiv,c_GsG6) 	    = Conductance_Coeff(iv5,cc_GsG6)
   SurfaceChar(gridiv,c_GsTH) 	    = Conductance_Coeff(iv5,cc_GsTH)
   SurfaceChar(gridiv,c_GsTL) 	    = Conductance_Coeff(iv5,cc_GsTL)
   SurfaceChar(gridiv,c_GsS1) 	    = Conductance_Coeff(iv5,cc_GsS1)
   SurfaceChar(gridiv,c_GsS2) 	    = Conductance_Coeff(iv5,cc_GsS2)
   SurfaceChar(gridiv,c_GsKmax)     = Conductance_Coeff(iv5,cc_GsKmax)
       
   ! ---- Find code for Anthropogenic heat ----
   call CodeMatchAnthropogenicHeat(rr,c_QFCode)    
   ! Transfer Anthropogenic heat characteristics to SurfaceChar
   SurfaceChar(gridiv,c_BaseTHDD)   = AnthropogenicHeat_Coeff(iv5,cA_BaseTHDD)
   SurfaceChar(gridiv,c_QF_A1) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_A1)
   SurfaceChar(gridiv,c_QF_B1) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_B1)
   SurfaceChar(gridiv,c_QF_C1) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_C1)
   SurfaceChar(gridiv,c_QF_A2) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_A2)
   SurfaceChar(gridiv,c_QF_B2) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_B2)
   SurfaceChar(gridiv,c_QF_C2) 	    = AnthropogenicHeat_Coeff(iv5,cA_QF_C2)
   SurfaceChar(gridiv,c_AHMin) 	    = AnthropogenicHeat_Coeff(iv5,cA_AHMin)
   SurfaceChar(gridiv,c_AHSlope)    = AnthropogenicHeat_Coeff(iv5,cA_AHSlope)
   SurfaceChar(gridiv,c_TCritic)    = AnthropogenicHeat_Coeff(iv5,cA_TCritic)
       
   ! ---- Find code for Irrigation ----
   call CodeMatchIrrigation(rr,c_IrrCode)    
   ! Transfer Irrigation characteristics to SurfaceChar
   SurfaceChar(gridiv,c_IeStart)     = Irrigation_Coeff(iv5,cIr_IeStart)
   SurfaceChar(gridiv,c_IeEnd)       = Irrigation_Coeff(iv5,cIr_IeEnd)
   SurfaceChar(gridiv,c_IntWU)       = Irrigation_Coeff(iv5,cIr_IntWU)
   SurfaceChar(gridiv,c_Faut)        = Irrigation_Coeff(iv5,cIr_Faut)
   SurfaceChar(gridiv,c_Ie_a)        = Irrigation_Coeff(iv5,cIr_Ie_a1:cIr_Ie_a3)
   SurfaceChar(gridiv,c_Ie_m)        = Irrigation_Coeff(iv5,cIr_Ie_m1:cIr_Ie_m3)
   SurfaceChar(gridiv,c_DayWat)      = Irrigation_Coeff(iv5,cIr_DayWat1:cIr_DayWat7)
   SurfaceChar(gridiv,c_DayWatPer)   = Irrigation_Coeff(iv5,cIr_DayWatPer1:cIr_DayWatPer7)
       
   ! ---- Find code for Hourly Profiles ----
   ! Energy use (weekdays)
   call CodeMatchProf(rr,c_EnProfWD)
   SurfaceChar(gridiv,c_HrProfEnUseWD) = Profiles_Coeff(iv5,cPr_Hours)
   ! Energy use (weekends)
   call CodeMatchProf(rr,c_EnProfWE)
   SurfaceChar(gridiv,c_HrProfEnUseWE) = Profiles_Coeff(iv5,cPr_Hours)
   ! Water use profile (manual, weekdays)
   call CodeMatchProf(rr,c_WProfManuWD)
   SurfaceChar(gridiv,c_HrProfWUManuWD)  = Profiles_Coeff(iv5,cPr_Hours)
   ! Water use profile (manual, weekends)
   call CodeMatchProf(rr,c_WProfManuWE)
   SurfaceChar(gridiv,c_HrProfWUManuWE)  = Profiles_Coeff(iv5,cPr_Hours)
   ! Water use profile (automatic, weekdays)
   call CodeMatchProf(rr,c_WProfAutoWD)
   SurfaceChar(gridiv,c_HrProfWUAutoWD) = Profiles_Coeff(iv5,cPr_Hours)
   ! Water use profile (automatic, weekends)
   call CodeMatchProf(rr,c_WProfAutoWE)
   SurfaceChar(gridiv,c_HrProfWUAutoWE) = Profiles_Coeff(iv5,cPr_Hours)
   ! Snow clearing profile (weekdays)
   call CodeMatchProf(rr,c_SnowProfWD)
   SurfaceChar(gridiv,c_HrProfSnowCWD) = Profiles_Coeff(iv5,cPr_Hours)
   ! Snow clearing profile (weekends)
   call CodeMatchProf(rr,c_SnowProfWE)
   SurfaceChar(gridiv,c_HrProfSnowCWE) = Profiles_Coeff(iv5,cPr_Hours)
   
   ! ---- Interpolate Hourly Profiles to model timestep and normalise
   TstepProfiles(Gridiv,:,:) = -999   !Initialise TstepProfiles
   ! Energy use
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_EnUseWD,c_HrProfEnUseWD)
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_EnUseWE,c_HrProfEnUseWE)
   ! For energy use, normalise so the AVERAGE of the multipliers is equal to 1
   TstepProfiles(Gridiv,cTP_EnUseWD,:) = TstepProfiles(Gridiv,cTP_EnUseWD,:) / sum(TstepProfiles(Gridiv,cTP_EnUseWD,:))*24*nsh_real
   TstepProfiles(Gridiv,cTP_EnUseWE,:) = TstepProfiles(Gridiv,cTP_EnUseWE,:) / sum(TstepProfiles(Gridiv,cTP_EnUseWE,:))*24*nsh_real
     
   ! Water use
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUManuWD,c_HrProfWUManuWD)
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUManuWE,c_HrProfWUManuWE)
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUAutoWD,c_HrProfWUAutoWD)
   call SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUAutoWE,c_HrProfWUAutoWE)
   ! For water use, normalise so the SUM of the multipliers is equal to 1 (profile is multiplied by daily water use)
   TstepProfiles(Gridiv,cTP_WUManuWD,:) = TstepProfiles(Gridiv,cTP_WUManuWD,:) / sum(TstepProfiles(Gridiv,cTP_WUManuWD,:))
   TstepProfiles(Gridiv,cTP_WUManuWE,:) = TstepProfiles(Gridiv,cTP_WUManuWE,:) / sum(TstepProfiles(Gridiv,cTP_WUManuWE,:))
   TstepProfiles(Gridiv,cTP_WUAutoWD,:) = TstepProfiles(Gridiv,cTP_WUAutoWD,:) / sum(TstepProfiles(Gridiv,cTP_WUAutoWD,:))
   TstepProfiles(Gridiv,cTP_WUAutoWE,:) = TstepProfiles(Gridiv,cTP_WUAutoWE,:) / sum(TstepProfiles(Gridiv,cTP_WUAutoWE,:))
      
 end subroutine InitializeSurfaceCharacteristics


!----------------------------------------------------------------------------------------------
!Calculates the initial conditions for each grid on the first year of run
!Made by sg feb 2012 -
!Latest modified:
!20 Oct 2014, LJ: Saves to slot 0 of the output matrix
!06 Nov 2014 HCW
!----------------------------------------------------------------------------------------------

 
!-------------------------------------------------------------------------
 SUBROUTINE InitialState(GridName,year_int,Gridiv,year_txt)
 ! Last modified by HCW 03 Dec 2014
 ! To do:
 !	 - Check running means (5-day temperature)
 !------------------------------------------------------------------------
 
  use allocateArray 
  use data_in 
  use ColNamesModelDailyState
  use defaultNotUsed
  use FileName
  use gis_data
  use InitialCond
  use mod_z    
  use resist  
  use snowMod
  use sues_data 
  use time
  
  IMPLICIT NONE
  
  character(len=20):: GridName    !Name of the evaluated grid
  character(len=10):: str2        !Variables related to filepaths
  character(len=150):: fileInit   !Initial conditions filename
  character(len=4):: year_txt     !year in txt format
  integer::DaysSinceRain,Gridiv,& !number of days since rain, grid number,
           gamma1,gamma2          !switches related to cooling and heating degree days
  integer::wd,seas,date,mb,&      !weekday information, season, date, month
           year_int,switch=0,&    !year as an integer, switch related to previous day
           id_next,calc           !next day,counter in irrigation calculations
   
  real (KIND(1d0))::PavedState,BldgsState,EveTrState,DecTrState,GrassState,BSoilState,&
              	    SnowFracPaved,SnowFracBldgs,SnowFracEveTr,SnowFracDecTr,&
                    SnowFracGrass,SnowFracBSoil,SnowFracWater,&
                    SnowDensPaved,SnowDensBldgs,SnowDensEveTr,SnowDensDecTr,&
                    SnowDensGrass,SnowDensBSoil,SnowDensWater

  !-----------------------------------------------------------------------

  namelist/InitialConditions/DaysSinceRain,&
                  Temp_C0,&
                  ID_Prev,&
                  GDD_1_0,&
                  GDD_2_0,&
                  LAIinitialEveTr,&            
                  LAIinitialDecTr,&
                  LAIinitialGrass,&
                  porosity0,&
                  albDec0,&
                  decidCap0,&
                  PavedState,&
                  BldgsState,&
	          EveTrState,&
	          DecTrState,&
	          GrassState,&
                  BSoilState,&
                  WaterState,&
                  SoilStorePavedState,&
                  SoilStoreBldgsState,&
                  SoilStoreEveTrState,&
                  SoilStoreDecTrState,&
                  SoilStoreGrassState,&
                  SoilStoreBSoilState,&
		  SnowWaterPavedState,&
		  SnowWaterBldgsState,&
		  SnowWaterEveTrState,&
                  SnowWaterDecTrState,&
                  SnowWaterGrassState,&
                  SnowWaterBSoilState,&
                  SnowWaterWaterState,&
                  SnowPackPaved,&
		  SnowPackBldgs,&
                  SnowPackEveTr,&
                  SnowPackDecTr,&
                  SnowPackGrass,&
                  SnowPackBSoil,&
                  SnowPackWater,&
                  SnowFracPaved,&
                  SnowFracBldgs,&
                  SnowFracEveTr,&
                  SnowFracDecTr,&
                  SnowFracGrass,&
                  SnowFracBSoil,&
                  SnowFracWater,&
                  SnowDensPaved,&
                  SnowDensBldgs,&
                  SnowDensEveTr,&
                  SnowDensDecTr,&
                  SnowDensGrass,&
                  SnowDensBSoil,&
                  SnowDensWater

  ! Define InitialConditions file ----------------------------------------
  FileInit=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//'.nml'

  ! Open, read and close InitialConditions file --------------------------
  open(56,File=trim(FileInit),err=600,status='old') !Change with needs
  read(56,iostat=ios_out,nml=InitialConditions,err=601)
  close(56)
  
  ! Write InitialConditions to FileChoices -------------------------------
  FileChoices=trim(FileOutputPath)//trim(FileCode)//'_FileChoices.txt'
  open(12,file=FileChoices,position='append')
  write(12,*)'----------',trim(FileInit),'----------'
  write(12,nml=InitialConditions)
  close(12)

  !-----------------------------------------------------------------------

  ! Previous day DOY number (needed in file allocations)
  if(id_prev>=364) id_prev=0  !If previous day is larger than 364, set this to zero
  
  ! -- Save id_prev to ModelDailyState array --
  ModelDailyState(Gridiv,cMDS_id_prev) = id_prev

  ! -- Save phenology info in InitialConditions to ModelDailyState array --
  ModelDailyState(Gridiv,cMDS_LAIInitialEveTr) = LAIInitialEveTr
  ModelDailyState(Gridiv,cMDS_LAIInitialDecTr) = LAIInitialDecTr
  ModelDailyState(Gridiv,cMDS_LAIInitialGrass) = LAIInitialGrass
  ModelDailyState(Gridiv,cMDS_GDD1_0) = GDD_1_0
  ModelDailyState(Gridiv,cMDS_GDD2_0) = GDD_2_0
  ModelDailyState(Gridiv,cMDS_GDDMin) =  90   !Going to check for minimum GDD - ??
  ModelDailyState(Gridiv,cMDS_GDDMax) = -90   !Going to check for maximum GDD - ??

  ! Initialize to the previous day's value (i.e. day before run starts)
  ModelDailyState(Gridiv,cMDS_porosity) = porosity0
  ModelDailyState(Gridiv,cMDS_albDec)   = albDec0
  ModelDailyState(Gridiv,cMDS_DecidCap) = DecidCap0
  ModelDailyState(Gridiv,cMDS_CumSnowfall) = 0 !!Check this

  ! -- Anthropogenic heat flux initializations -- 
  ! Need to get BaseTHDD from SurfaceChar, as info not transferred until SUEWS_Translate called 
  BaseTHDD = SurfaceChar(Gridiv,c_BaseTHDD)
    
  if(AnthropHeatChoice>=0) then    
     !Calculations related to heating and cooling degree days (BaseT is used always)   
     if ((Temp_C0-BaseTHDD)>=0) then   !Cooling
        gamma2=1
     else
        gamma2=0
     endif
     if ((BaseTHDD-Temp_C0)>=0) then   !Heating
        gamma1=1
     else
        gamma1=0
     endif
     ModelDailyState(Gridiv,cMDS_HDD1) = gamma1*(BaseTHDD-Temp_C0) ! Heating
     ModelDailyState(Gridiv,cMDS_HDD2) = gamma2*(Temp_C0-BaseTHDD) ! Cooling
  endif

  ModelDailyState(Gridiv,cMDS_TempC) = Temp_C0
  ModelDailyState(Gridiv,cMDS_DaysSinceRain) = real(DaysSinceRain,kind(1d0))

  ! Assume that the temperature has been the same for the previous days
  ModelDailyState(Gridiv,cMDS_TempCOld1) = Temp_C0
  ModelDailyState(Gridiv,cMDS_TempCOld2) = Temp_C0
  ModelDailyState(Gridiv,cMDS_TempCOld3) = Temp_C0
  
  !! These variables don't seem to be needed (commented out HCW 27 Nov 2014)
  !! If required, they will need updating for a non-hourly timestep
  !! Initialize hourly temperature and precipitation + 5 day mean used to define thermal growing season
  !runT(0:23)=Temp_C0 
  !runP(0:23)=0
    
  finish=.false.  
  
  ! -- Save snow density info in InitialConditions to ModelDailyState array --  
  ModelDailyState(Gridiv,cMDS_SnowDens(PavSurf))    = SnowDensPaved   
  ModelDailyState(Gridiv,cMDS_SnowDens(BldgSurf))   = SnowDensBldgs
  ModelDailyState(Gridiv,cMDS_SnowDens(ConifSurf))  = SnowDensEveTr
  ModelDailyState(Gridiv,cMDS_SnowDens(DecidSurf))  = SnowDensDecTr
  ModelDailyState(Gridiv,cMDS_SnowDens(GrassSurf))  = SnowDensGrass
  ModelDailyState(Gridiv,cMDS_SnowDens(BSoilSurf))  = SnowDensBSoil
  ModelDailyState(Gridiv,cMDS_SnowDens(WaterSurf))  = SnowDensWater
 
  !! Where is this from??   
  IceFrac=0.2   !Estimated fraction of ice. Should be improved in the future

  ! ==============================================================
  ! ============ Save states to ModelOutputData ==================    

     
  ! -- Initial wetness status of each surface (above ground) --
  ModelOutputData(0,cMOD_State(PavSurf),   Gridiv) = PavedState
  ModelOutputData(0,cMOD_State(BldgSurf),  Gridiv) = BldgsState
  ModelOutputData(0,cMOD_State(ConifSurf), Gridiv) = EveTrState
  ModelOutputData(0,cMOD_State(DecidSurf), Gridiv) = DecTrState
  ModelOutputData(0,cMOD_State(GrassSurf), Gridiv) = GrassState
  ModelOutputData(0,cMOD_State(BSoilSurf), Gridiv) = BSoilState
  ModelOutputData(0,cMOD_State(WaterSurf), Gridiv) = WaterState
  
  ! -- Initial soil stores for each surface (below ground) --
  ModelOutputData(0,cMOD_SoilState(PavSurf),   Gridiv) = SoilStorePavedState
  ModelOutputData(0,cMOD_SoilState(BldgSurf),  Gridiv) = SoilStoreBldgsState
  ModelOutputData(0,cMOD_SoilState(ConifSurf), Gridiv) = SoilStoreEveTrstate
  ModelOutputData(0,cMOD_SoilState(DecidSurf), Gridiv) = SoilStoreDecTrState
  ModelOutputData(0,cMOD_SoilState(GrassSurf), Gridiv) = SoilStoreGrassState
  ModelOutputData(0,cMOD_SoilState(BSoilSurf), Gridiv) = SoilStoreBSoilState
  ModelOutputData(0,cMOD_SoilState(WaterSurf), Gridiv) = 0 ! No soil layer for water surface
  
  ! -- Initial snow water equivalent for each surface --
  ModelOutputData(0,cMOD_SnowPack(PavSurf),   Gridiv) = SnowPackPaved
  ModelOutputData(0,cMOD_SnowPack(BldgSurf),  Gridiv) = SnowPackBldgs
  ModelOutputData(0,cMOD_SnowPack(ConifSurf), Gridiv) = SnowPackEveTr
  ModelOutputData(0,cMOD_SnowPack(DecidSurf), Gridiv) = SnowPackDecTr
  ModelOutputData(0,cMOD_SnowPack(GrassSurf), Gridiv) = SnowPackGrass
  ModelOutputData(0,cMOD_SnowPack(BSoilSurf), Gridiv) = SnowPackBSoil
  ModelOutputData(0,cMOD_SnowPack(WaterSurf), Gridiv) = SnowPackWater
  
  ! -- Initial liquid (melted) water for each surface --
  ModelOutputData(0,cMOD_SnowWaterState(PavSurf),   Gridiv) = SnowWaterPavedState
  ModelOutputData(0,cMOD_SnowWaterState(BldgSurf),  Gridiv) = SnowWaterBldgsState
  ModelOutputData(0,cMOD_SnowWaterState(ConifSurf), Gridiv) = SnowWaterEveTrstate
  ModelOutputData(0,cMOD_SnowWaterState(DecidSurf), Gridiv) = SnowWaterDecTrState
  ModelOutputData(0,cMOD_SnowWaterState(GrassSurf), Gridiv) = SnowWaterGrassState
  ModelOutputData(0,cMOD_SnowWaterState(BSoilSurf), Gridiv) = SnowWaterBSoilState
  ModelOutputData(0,cMOD_SnowWaterState(WaterSurf), Gridiv) = SnowWaterWaterState
 
  ! -- Initial fraction of snow on each surface --
  ModelOutputData(0,cMOD_SnowFrac(PavSurf),   Gridiv) = SnowFracPaved
  ModelOutputData(0,cMOD_SnowFrac(BldgSurf),  Gridiv) = SnowFracBldgs
  ModelOutputData(0,cMOD_SnowFrac(ConifSurf), Gridiv) = SnowFracEveTr
  ModelOutputData(0,cMOD_SnowFrac(DecidSurf), Gridiv) = SnowFracDecTr
  ModelOutputData(0,cMOD_SnowFrac(GrassSurf), Gridiv) = SnowFracGrass
  ModelOutputData(0,cMOD_SnowFrac(BSoilSurf), Gridiv) = SnowFracBSoil
  ModelOutputData(0,cMOD_SnowFrac(WaterSurf), Gridiv) = SnowFracWater 
     
  !At this point translate arrays to variables (needed for RoughnessParameters)
  call SUEWS_Translate(Gridiv,0,0)  
    
  !Calculation of roughness parameters (N.B. uses porosity)
  call RoughnessParameters

 !=============================================================================
 ! If the run start day is at previous year, then calculate the number of days
 ! in that year.

 !First we need to know if the previous day given in initial conditions (id_prev) is
 !on previous year as this is needed in the initialization of dayofWeek matrix.
 !In this case switch is set to one for date calculations.
 if(id_prev==0)then                     !If id_prev = 0, means that the first modelled day is 1 Jan
   year_int=year_int-1                  !1) find the number of days on that year
   call LeapYearCalc (year_int,id_prev) !2) set switch to 1 so that the code knows to change back to current year (switch=0)
   switch=1
 endif

 call day2month(id_prev,mb,date,seas,year_int,lat) !Calculate date information (mb = month, date = day,...)
 call Day_of_Week(date,mb,year_int,wd)             !Calculate weekday of the previous day (wd) (1=Sun, ..., 7=Sat)

 !After the day in previous year switch is changed back to zero: 
 !    ie not previous day anymore
 !Also the size of dayofweek is from 0:NdaysinYear meaning 
 !that in zero slot is the previous day information
 if(switch==1)then
   year_int=year_int+1
   id_prev=0
   switch=0
 endif
  
 dayofWeek(id_prev,1)=wd   ! day of week
 dayofWeek(id_prev,2)=mb   ! month
 dayofweek(id_prev,3)=seas ! season (summer=1, winter=2) needed for accumulation

 ! in case next day goes to next year calculate again the date information for dayofweek matrix.
 id_next=id_prev+1
 if(id_next>nofDaysThisYear) then
   id_next=1
   year_int=year_int+1
   switch=1
   call ErrorHint(43,'switch- years',notUsed,notUsed,notUsedI)
 endif

 call day2month(id_next,mb,date,seas,year_int,lat) !Calculate real date from doy
 call Day_of_Week(date,mb,year_int,wd)             !Calculate weekday (1=Sun, ..., 7=Sat)
 
 if(switch==1)then
   iy=iy-1
   switch=0
 endif

 dayofWeek(id_next,1)=wd  ! day of week
 dayofWeek(id_next,2)=mb  ! month
 dayofweek(id_next,3)=seas ! season

 !=============================================================================
 
 !id=id_prev
 !it= 23 !!LastTimeOfDay
 if (id_prev>=DayLightSavingDay(1).and.id_prev<=DayLightSavingDay(2)) then  !Summertime
    DLS=1
 else
    DLS=0
 endif
          
 ! -----------------------------------------------------------------------
 ! Calculate daily water use if modelled (i.e. if WU_choice = 0). 
 ! Calculated from previous day information given in InitialConditions file
 
 WU_day=0                !Initialize WU_day
 if (WU_choice==0) then  !Model water use
    calc=0

    if (DayWat(wd)==1.0) then !if DayWat(wd)=1.0 (irrigation occurs on this day)
       if (lat>=0) then            !Northern Hemisphere
          if (id>=Ie_start.and.id<=Ie_end) calc=1 !if day between irrigation period               
       else                        !Southern Hemisphere
          calc=1
          if (id>=Ie_end.and.id<=Ie_start) calc=0 !if day between irrigation period                       
       endif
       if(calc==1) then                     
          ! Model daily water use based on HDD(id,6)(days since rain) and HDD(id,3)(average temp)
          ! ---- Automatic irrigation (evergreen trees) ----
	      WU_day(id,2) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(ConifSurf)*IrrFracConif*DayWatPer(wd)
	      if (WU_Day(id,2)<0) WU_Day(id,2)=0   !If modelled WU is negative -> 0
	      ! ---- Manual irrigation (evergreen trees) ----
	      WU_day(id,3) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(ConifSurf)*IrrFracConif*DayWatPer(wd)
	      if (WU_Day(id,3)<0) WU_Day(id,3)=0   !If modelled WU is negative -> 0
	      ! ---- Total evergreen trees water use (automatic + manual) ----
	      WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))
	                              
	      ! ---- Automatic irrigation (deciduous trees) ----
	      WU_day(id,5) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(DecidSurf)*IrrFracDecid*DayWatPer(wd)
	      if (WU_Day(id,5)<0) WU_Day(id,5)=0   !If modelled WU is negative -> 0
	      ! ---- Manual irrigation (deciduous trees) ----
	      WU_day(id,6) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(DecidSurf)*IrrFracDecid*DayWatPer(wd)
	      if (WU_Day(id,6)<0) WU_Day(id,6)=0   !If modelled WU is negative -> 0
	      ! ---- Total deciduous trees water use (automatic + manual) ----
	      WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))
	                              
	      ! ---- Automatic irrigation (grass) ----
	      WU_day(id,8) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(GrassSurf)*IrrFracGrass*DayWatPer(wd)
	      if (WU_Day(id,8)<0) WU_Day(id,8)=0   !If modelled WU is negative -> 0
	      ! ---- Manual irrigation (grass) ----
	      WU_day(id,9) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(GrassSurf)*IrrFracGrass*DayWatPer(wd)
	      if (WU_Day(id,9)<0) WU_Day(id,9)=0   !If modelled WU is negative -> 0
	      ! ---- Total grass water use (automatic + manual) ----
          WU_Day(id,7)=(WU_day(id,8)+WU_day(id,9))
       else
          WU_Day(id,1)=0
          WU_Day(id,2)=0
          WU_Day(id,3)=0
          WU_Day(id,4)=0
          WU_Day(id,5)=0
          WU_Day(id,6)=0
          WU_Day(id,7)=0
          WU_Day(id,8)=0
          WU_Day(id,9)=0
       endif
    endif
 endif

 ! -----------------------------------------------------------------------
   
 !Initialise rates of change variables for OHM calculation 
 ! check ?? this should not happen at the start of the year if continuing

 q1_grids(Gridiv)=-101
 q2_grids(Gridiv)=-100
 q3_grids(Gridiv)=-99

 r1_grids(Gridiv)=-101
 r2_grids(Gridiv)=-100
 r3_grids(Gridiv)=-99


 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !DEFINE DIFFERENT INITIALIZATION PARAMETERS
 !Are still needed ??
 !once=.true.
 return

600 call ErrorHint(47,trim(FileInit),notUsed,notUsed,notUsedI)
601 call ErrorHint(48,trim(FileInit),notUsed,notUsed,ios_out)

 end subroutine InitialState

!-------------------------------------------------------------------------
 subroutine NextInitial(GridName,year_int)
 ! Modified by HCW 21 Nov 2014
 ! Last day of year is not anymore the number of days on that year, but rather
 ! id == 1. Thus nofDaysThisYear was changed to 1. LJ 9/4/2015
 !------------------------------------------------------------------------
  use allocateArray 
  use ColNamesInputFiles
  use ColNamesModelDailyState
  use data_in
  use defaultNotUsed
  use Initial
  use sues_data
  use snowMod
  use time
  
  IMPLICIT NONE

  character (len=15)::GridName
  character (len=4)::year_txt2
  integer:: year_int2
  integer:: year_int

  year=year_int   !HCW added 21 Nov 2014

  if (id==1) then  !nofDaysThisYear changed to 1
     year_int2=int(year+1)
     write(year_txt2,'(I4)')year_int2
     open(57,File=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//'.nml',err=200)
  else
        year_int2=int(year)
        write(year_txt2,'(I4)')year_int2
        open(57,File=trim(FileInputPath)//trim("InitialConditions")//trim(GridName)//'end.nml',err=201)
     endif
  !endif

  write(57,*)'&InitialConditions'
  write(57,*)'DaysSinceRain=',int(HDD(id,6))
   
  !! If last time of day, then DailyState variables will have been updated so can write out arrays for id rather than id-1
  !if(it==23 .and. imin == (nsh_real-1)/nsh_real*60) then  !!LastTimeofday
  !   id=id+1
  !endif
  
  write(57,*)'Temp_C0=',HDD(nofDaysThisYear,3)
  write(57,*)'ID_Prev=',nofDaysThisYear
  write(57,*)'GDD_1_0=',GDD(nofDaysThisYear,1)
  write(57,*)'GDD_2_0=',GDD(nofDaysThisYear,2)
  write(57,*)'PavedState=',State(PavSurf)
  write(57,*)'BldgsState=',State(BldgSurf)
  write(57,*)'EveTrState=',State(ConifSurf)
  write(57,*)'DecTrState=',State(DecidSurf)
  write(57,*)'GrassState=',State(GrassSurf)
  write(57,*)'BSoilState=',State(BSoilSurf)
  write(57,*)'WaterState=',State(WaterSurf)
  write(57,*)'LAIinitialEveTr=',lai(nofDaysThisYear,ivConif)
  write(57,*)'LAIinitialDecTr=',lai(nofDaysThisYear,ivDecid)
  write(57,*)'LAIinitialGrass=',lai(nofDaysThisYear,ivGrass)
  write(57,*)'porosity0=',porosity(id)
  write(57,*)'DecidCap0=',decidCap(id)
  write(57,*)'albDec0=',AlbDec(id)
  write(57,*)'soilstorePavedState=',soilmoist(PavSurf)
  write(57,*)'soilstoreBldgsState=',soilmoist(BldgSurf)
  write(57,*)'soilstoreEveTrState=',soilmoist(ConifSurf)
  write(57,*)'soilstoreDecTrState=',soilmoist(DecidSurf)
  write(57,*)'soilstoreGrassState=',soilmoist(GrassSurf)
  write(57,*)'soilstoreBSoilState=',soilmoist(BSoilSurf)
  write(57,*)'SnowWaterPavedstate=',MeltWaterStore(PavSurf)
  write(57,*)'SnowWaterBldgsState=',MeltWaterStore(BldgSurf)
  write(57,*)'SnowWaterEveTrState=',MeltWaterStore(ConifSurf)
  write(57,*)'SnowWaterDecTrState=',MeltWaterStore(DecidSurf)
  write(57,*)'SnowWaterGrassState=',MeltWaterStore(GrassSurf)
  write(57,*)'SnowWaterBSoilState=',MeltWaterStore(BSoilSurf)
  write(57,*)'SnowWaterWaterState=',MeltWaterStore(WaterSurf)
  write(57,*)'SnowPackPaved=',SnowPack(PavSurf)
  write(57,*)'SnowPackBldgs=',SnowPack(BldgSurf)
  write(57,*)'SnowPackEveTr=',SnowPack(ConifSurf)
  write(57,*)'SnowPackDecTr=',SnowPack(DecidSurf)
  write(57,*)'SnowPackGrass=',SnowPack(GrassSurf)
  write(57,*)'SnowPackBSoil=',SnowPack(BSoilSurf)
  write(57,*)'SnowPackWater=',SnowPack(WaterSurf)
  write(57,*)'SnowFracPaved=',SnowFrac(PavSurf)
  write(57,*)'SnowFracBldgs=',SnowFrac(BldgSurf)
  write(57,*)'SnowFracEveTr=',SnowFrac(ConifSurf)
  write(57,*)'SnowFracDecTr=',SnowFrac(DecidSurf)
  write(57,*)'SnowFracGrass=',SnowFrac(GrassSurf)
  write(57,*)'SnowFracBSoil=',SnowFrac(BSoilSurf)
  write(57,*)'SnowFracWater=',SnowFrac(WaterSurf)
  write(57,*)'SnowDensPaved=',densSnow(PavSurf)
  write(57,*)'SnowDensBldgs=',densSnow(BldgSurf)
  write(57,*)'SnowDensEveTr=',densSnow(ConifSurf)
  write(57,*)'SnowDensDecTr=',densSnow(DecidSurf)
  write(57,*)'SnowDensGrass=',densSnow(GrassSurf)
  write(57,*)'SnowDensBSoil=',densSnow(BSoilSurf)
  write(57,*)'SnowDensWater=',densSnow(WaterSurf)
  write(57,*)'/'
  close(57)
  
  if(it==23 .and. imin == (nsh_real-1)/nsh_real*60) then
     id=id-1  
  endif
  
  return

200 call ErrorHint(49,trim("InitialConditions")//trim(GridName)// &
		   '_'//trim(adjustl(year_txt2))//'.nml',notUsed,notUsed,notUsedI)
201 call ErrorHint(49,trim("InitialConditions")//trim(GridName)// &
		   '_'//trim(adjustl(year_txt2))//'end.nml',notUsed,notUsed,notUsedI)

 end subroutine NextInitial
!------------------------------------------------------------------------- 



 !=======================================================================
 !=======================================================================
 !This subroutine prepares a meteorological forcing file.

 subroutine SUEWS_InitializeMetData(lunit)

   use allocateArray
   use data_in
   use sues_data
   use time
   use defaultnotUsed
   use Initial

   IMPLICIT NONE

   integer::lunit,i,iyy,RunNumber!,NSHcounter
   real (kind(1d0)),dimension(24)::MetArray
   real(kind(1d0)):: imin_prev, ih_prev, iday_prev, tstep_met   !For checks on temporal resolution of met data

   !---------------------------------------------------------------

   !Open the file for reading and read the actual data
   !write(*,*) fileMet
   open(lunit,file=trim(fileMet),status='old',err=314)
   call skipHeader(lunit,SkipHeaderMet)

   ! Skip to the right place in the met file, depending on how many chunks have been read already
   if (skippedLines>0) then
       do iyy=1,skippedLines
       read(lunit,*)
       enddo
   endif

   ! Read in next chunk of met data and fill MetForcingData array with data for every timestep
   !NSHcounter = 1
   DO i=1,ReadlinesMetdata
      call MetRead(MetArray,InputmetFormat,ldown_option,NetRadiationChoice,&
                   snowUse,smd_choice,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)
      !DO iv=1,NSH
      !    MetForcingData(NSHcounter,1:24,GridCounter) = MetArray
      !   NSHcounter = NSHcounter + 1
      !ENDDO
      MetForcingData(i,1:24,GridCounter) = MetArray
      ! Check timestamp of met data file matches TSTEP specified in RunControl
      if(i==1) then
         imin_prev = MetArray(4)
         ih_prev   = MetArray(3)
         iday_prev   = MetArray(2)
      elseif(i==2) then
         tstep_met = ((MetArray(4)+60*MetArray(3)) - (imin_prev+60*ih_prev))*60   !tstep in seconds
         if(tstep_met.ne.tstep_real.and.MetArray(2)==iday_prev) then
            call ErrorHint(39,'TSTEP in RunControl does not match TSTEP of met data (DOY).',real(tstep,kind(1d0)),tstep_met,&
                           int(MetArray(2)))        
         endif    
      endif   
        
   ENDDO

   CLOSE(lunit)

   return

314 call errorHint(11,trim(fileMet),notUsed,notUsed,ios_out)


 end subroutine SUEWS_InitializeMetData
!----------------------------------------------------------------------------------------------


 !====================================================================================
 subroutine CheckInitial  
 !Check the parameters in InitialConditions file.
 !Modified by HCW 04 Mar 2014, changed soilmoist(is) checks to use names given in InitialConditions 
 !Added by LJ in 8/2/2013
 
 use allocateArray
 use data_in
 use defaultNotUsed
 use InitialCond
 use snowMod
 use time

 implicit none

 real(kind(1d0)):: pTol   !Precision tolerance for range checks
 
 if (Temp_C0<(Temp_C-10).or.Temp_C0>(Temp_C+10)) then
     call ErrorHint(36,'InitialCond: Check temperature', Temp_C0, Temp_C, notUsedI)
 endif
 
 if (ID_Prev/=id-1) then
   call ErrorHint(36,'InitialCond: Check previous day', real(ID_Prev,kind(1d0)), real(id,kind(1d0)), notUsedI)
 endif

  !!This part currently does not work for multiple grids as Initial conditions values get overwritten.
  ! Simple checks that Initial Conditions are within the specified ranges (within a precision tolerance)
  pTol = 0.00001  
  !write(*,*) LAIInitialEveTr,LAImin(ConifSurf-2),LAImax(ConifSurf-2)
  if(LAIInitialEveTr < (LAImin(ConifSurf-2)-pTol)) then
     call ErrorHint(36,'Intial LAI for EveTr < min value in SUEWS_Veg.txt!', LAIMin(ConifSurf-2), LAIInitialEveTr, notUsedI)
  endif
  if(LAIInitialEveTr > (LAImax(ConifSurf-2)+pTol)) then
     call ErrorHint(36,'Intial LAI for EveTr > max value in SUEWS_Veg.txt!', LAIMax(ConifSurf-2), LAIInitialEveTr, notUsedI)
  endif
  !write(*,*) LAIInitialDecTr,LAImin(DecidSurf-2)
  if(LAIInitialDecTr < (LAImin(DecidSurf-2)-pTol)) then
     call ErrorHint(36,'Intial LAI for DecTr < min value in SUEWS_Veg.txt!', LAIMin(DecidSurf-2), LAIInitialDecTr, notUsedI)
  endif
  if(LAIInitialDecTr > (LAImax(DecidSurf-2)+pTol)) then
     call ErrorHint(36,'Intial LAI for DecTr > max value in SUEWS_Veg.txt!', LAIMax(DecidSurf-2), LAIInitialDecTr, notUsedI)
  endif
  if(LAIInitialGrass < (LAImin(GrassSurf-2)-pTol)) then
     call ErrorHint(36,'Intial LAI for Grass < min value in SUEWS_Veg.txt!', LAIMin(GrassSurf-2), LAIInitialGrass, notUsedI)
  endif
  if(LAIInitialGrass > (LAImax(GrassSurf-2)+pTol)) then
    call ErrorHint(36,'Intial LAI for Grass > max value in SUEWS_Veg.txt!', LAIMax(GrassSurf-2), LAIInitialGrass, notUsedI)
  endif

  !write(*,*) AlbDec0, albmin_dec, albmax_dec
  if(AlbDec0 < (AlbMin_dec-pTol)) then
     !call ErrorHint(36,'Intial albedo for DecTr < min value in SUEWS_Veg.txt!', AlbMin_dec, AlbDec0, notUsedI)
     call ErrorHint(36,'Intial albedo for DecTr < min value!', AlbMin_dec, AlbDec0, notUsedI)
  endif
  if(AlbDec0 > (AlbMax_dec+pTol)) then
     !call ErrorHint(36,'Intial albedo for DecTr > max value in SUEWS_Veg.txt!', AlbMax_dec, AlbDec0, notUsedI)
     call ErrorHint(36,'Intial albedo for DecTr < max value!', AlbMax_dec, AlbDec0, notUsedI)
  endif 
  !DecidCap0, Porosity0...
 
 !Check more thoroughly if LAI values are OK. Need to treat different hemispheres as well as 
 !tropics separately.
 if (lat>40) then
   if ((LAIinitialEveTr>LAImin(ConifSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialEveTr<LAImax(ConifSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImin(ConifSurf), notUsedI)
   endif
   if ((LAIinitialDecTr>LAImin(DecidSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialDecTr<LAImax(DecidSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImin(DecidSurf), notUsedI)
   endif
   if ((LAIinitialGrass>LAImin(GrassSurf-2)+1.and.(id<60.or.id>330)).or.&
     (LAIinitialGrass<LAImax(GrassSurf-2)-1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImin(GrassSurf), notUsedI)
   endif    
   
 elseif (lat<-40) then
   if ((LAIinitialEveTr<LAImax(ConifSurf-2)-1.and.(id<60.or.id>330)).or.&
     (LAIinitialEveTr>LAImin(ConifSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImax(ConifSurf), notUsedI)
   endif
   if ((LAIinitialDecTr>LAImax(DecidSurf-2)-1.and.(id<60.or.id>330)).or.&
     (LAIinitialDecTr>LAImin(DecidSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImax(DecidSurf), notUsedI)
   endif
   if ((LAIinitialGrass<LAImax(GrassSurf-2)-1.and.(id<60.or.id>330)) .or.&
     (LAIinitialGrass>LAImin(GrassSurf-2)+1.and.(id>130.and.id<244))) then
      call ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImax(GrassSurf), notUsedI)
   endif    
 
 elseif (lat<10.and.lat>-10) then
 
   if (LAIinitialEveTr<LAImax(ConifSurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImax(ConifSurf), notUsedI)
   endif
   if (LAIinitialDecTr<LAImax(DecidSurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImax(DecidSurf), notUsedI)
   endif
   if (LAIinitialGrass<LAImax(GrassSurf-2)-0.5) then
      call ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImax(GrassSurf), notUsedI)
   endif
  
 endif

 !Soilstore check
 if (SoilstoreBldgsState>soilstoreCap(BldgSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of building soil store.',&
                   SoilstoreBldgsState, soilstoreCap(BldgSurf), notUsedI)
 endif
 if (SoilstorePavedState>soilstoreCap(PavSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of paved soil store.',&
                   SoilstorePavedState, soilstoreCap(PavSurf), notUsedI)
 endif
 if (SoilstoreEveTrState>soilstoreCap(ConifSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of conif soil store.',&
                   SoilstoreEveTrState, soilstoreCap(ConifSurf), notUsedI)
 endif
 if (SoilstoreDecTrState>soilstoreCap(DecidSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of deciduous soil store.',&
                   SoilstoreDecTrState, soilstoreCap(DecidSurf), notUsedI)
 endif
 if (SoilstoreBSoilState>soilstoreCap(BSoilSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of bare soil soil store.',&
                   SoilstoreBSoilState, soilstoreCap(BSoilSurf), notUsedI)
 endif
 if (SoilstoreGrassState>soilstoreCap(GrassSurf)) then
    call ErrorHint(36,'InitialCond: Check initial condition of grass soil store.',&
                   SoilstoreGrassState, soilstoreCap(GrassSurf), notUsedI)
 endif

 !Snow stuff
 if (snowUse==1) then
    if (SnowWaterBldgsState>CRWmax*SnowPackBldgs) then
       call ErrorHint(36,'InitialCond: SnowWaterBldgsState', SnowWaterBldgsState, SnowPackBldgs, notUsedI)
    endif 
    if (SnowWaterPavedState>CRWmax*SnowPackPaved) then
       call ErrorHint(36,'InitialCond: SnowWaterPavedState', SnowWaterPavedState, SnowPackPaved, notUsedI)
    endif 
    if (SnowWaterEveTrState>CRWmax*SnowPackEveTr) then
       call ErrorHint(36,'InitialCond: SnowWaterEveTrstate', SnowWaterEveTrstate, SnowPackEveTr, notUsedI)
    endif 
    if (SnowWaterDecTrState>CRWmax*SnowPackDecTr) then
       call ErrorHint(36,'InitialCond: SnowWaterDecTrState', SnowWaterDecTrState, SnowPackDecTr, notUsedI)
    endif 
    if (SnowWaterGrassState>CRWmax*SnowPackGrass) then
       call ErrorHint(36,'InitialCond: SnowWaterGrassState', SnowWaterGrassState, SnowPackGrass, notUsedI)
    endif 
    if (SnowWaterBSoilState>CRWmax*SnowPackBSoil) then
       call ErrorHint(36,'InitialCond: SnowWaterGrassUnirState', SnowWaterBSoilState, SnowPackBSoil, notUsedI)
    endif 
 endif

!SnowWaterWaterstate,& ??
!SnowPackWater,& ??

 
 end subroutine CheckInitial  


  
