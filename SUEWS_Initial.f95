!=========================================================================
! sg feb 2012
! only run once at start - fixed for all grids and all years

SUBROUTINE OverallRunControl
  ! Last modified:
  ! LJ 27 Jan 2016  - Removal of tabs, cleaning of the code
  ! HCW 06 Mar 2015 - Removed options 10,20,30 (NARPOutput) for NetRadiationChoice
  ! HCW 06 Feb 2015 - File ID numbers changed so they are unique
  ! HCW 19 Dec 2014
  ! To Do:
  !   - Holidays.txt input file needs to be read in and coded into model
  !  - Add column header checks for SiteSelect
  !-------------------------------------------------------------------------

  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE FileName
  USE initial
  USE gis_data
  USE mod_z
  USE resist
  USE snowMod
  USE sues_data
  USE time

  IMPLICIT NONE

  INTEGER:: iv,i,ii,SkipCounter            !iv and i, ii are integers used in do loops
  CHARACTER(len=50):: FileN

  ! ---- Namelist for RunControl.nml ----
  NAMELIST/RunControl/AnthropHeatChoice,&
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
  OPEN(55,File='RunControl.nml',err=200,status='old') !Change with needs
  READ(55,nml=RunControl,err=201)
  CLOSE(55)

  !Check for problems with FileCode
  IF (FileCode=='none') CALL ErrorHint(26,TRIM("RunControl.nml FileCode is missing"),notUsed,notUsed,notUsedI)

  !-----------------------------------------------------------------------

  !Write RunControl information to FileChoices.txt
  FileChoices=TRIM(FileOutputPath)//TRIM(FileCode)//'_FileChoices.txt'
  OPEN (12,file=FileChoices,err=203)
  WRITE(12,nml=RunControl)
  CLOSE(12)

  !Determine what should be done with respect to radiation
  AlbedoChoice=0
  ldown_option=0
  IF(netRadiationChoice==0)THEN    !Observed Q* from the met input file will be used
     IF(snowUse==1) THEN            !If snow is modelled, NARP is needed for surface temperature
        netRadiationChoice=3000
        ldown_option=3              !Ldown will be modelled
        !NetRadiationChoice=NetRadiationChoice/1000
     ENDIF

  ELSEIF(netRadiationChoice>0)THEN  !Modelled Q* is used (NARP)
     AlbedoChoice=-9
     IF(NetRadiationChoice<10) THEN
        AlbedoChoice=0
        IF(NetRadiationChoice==1)ldown_option=1
        IF(NetRadiationChoice==2)ldown_option=2
        IF(NetRadiationChoice==3)ldown_option=3

     ELSEIF(NetRadiationChoice>=100.AND.NetRadiationChoice<1000) THEN
        AlbedoChoice=1
        IF(NetRadiationChoice==100)ldown_option=1
        IF(NetRadiationChoice==200)ldown_option=2
        IF(NetRadiationChoice==300)ldown_option=3
        NetRadiationChoice=NetRadiationChoice/100
     ENDIF

     !If bad NetRadiationChoice value
     IF(netRadiationChoice>3.OR. AlbedoChoice==-9)THEN
        WRITE(*,*) 'NetRadiationChoice=',NetRadiationChoice
        WRITE(*,*) 'Value not usable'
        STOP
     ENDIF
  ENDIF
  IF (gsChoice>2)THEN
     WRITE(*,*) 'gsChoice=', gsChoice
     WRITE(*,*) 'value too large'
     STOP
  END IF

  !------------------------------------------------------------------
  !Print run information on the screen
  WRITE(*,*)'--------------------------------------------------------'
  WRITE(*,*)"LUMPS/Suews 2016 - relevant references"
  WRITE(*,*)"LUMPS - Grimmond and Oke (2002) JAM, 41, 79-810"
  WRITE(*,*)"OHM - Grimmond and Oke (1999) JAM, 38, 922-940"
  WRITE(*,*)"NARP - Offerle et al. (2003) JAM"
  WRITE(*,*)"SUES - Evaporation Grimmond & Oke (1991) WRR"
  WRITE(*,*)"Water Balance Model Grimmond et al. (1986) WRR"
  WRITE(*,*)"NARP - Long wave improvements (Loridan et al. 2011 JAMC)"
  WRITE(*,*)"SUEWS - Anthropogenic heat, etc (Jarvi et al. 2011 JH)"
  WRITE(*,*)"SUEWS - Snow module included (Jarvi et al. 2014 GMD)"
  WRITE(*,*)'--------------------------------------------------------'


  !=======================================================================
  !======================== Read input files =============================
  ! This part reads the input files derived from the SiteInfo spreadsheet

  WRITE(*,*) 'Reading the following input files:'

  !=======================SUEWS_SiteSelect.txt============================
  FileN='SUEWS_SiteSelect.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesSiteSelect=nlines
  ALLOCATE(SiteSelect(nlinesSiteSelect,ncolumnsSiteSelect))
  !Read input file
  OPEN(21,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(21,*)   !Skip lines before header
  ENDDO
  READ(21,*) (HeaderSiteSelect_File(iv),iv=1,ncolumnsSiteSelect) !Get header

  DO i=1,nlinesSiteSelect
     READ(21,*) (SiteSelect(i,iv),iv=1,ncolumnsSiteSelect)
     !write(*,*) (SiteSelect(i,iv),iv=1,ncolumnsSiteSelect)
  ENDDO
  CLOSE(21)

  !call InputHeaderCheck(FileN) !! Need to add column checks for SiteSelect.txt

  !=======================SUEWS_NonVeg.txt============================
  FileN='SUEWS_NonVeg.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesNonVeg=nlines
  ALLOCATE(NonVeg_Coeff(nlinesNonVeg,ncolumnsNonVeg))
  !Read input file
  OPEN(22,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(22,*)   !Skip lines before header
  ENDDO
  READ(22,*) (HeaderNonVeg_File(iv),iv=1,ncolumnsNonVeg) !Get header

  DO i=1,nlinesNonVeg
     READ(22,*) (NonVeg_Coeff(i,iv),iv=1,ncolumnsNonVeg)
     !write(*,*) (NonVeg_Coeff(i,iv),iv=1,ncolumnsNonVeg)
  ENDDO
  CLOSE(22)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesNonVeg
     DO ii=1,nlinesNonVeg
        IF(NonVeg_Coeff(i,ci_Code)==NonVeg_Coeff(ii,ci_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',NonVeg_Coeff(i,ci_Code),'in Impervious.txt not unique!'
           CALL ErrorHint(60,FileN,NonVeg_Coeff(i,ci_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !=======================SUEWS_Veg.txt==============================
  FileN='SUEWS_Veg.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesVeg=nlines
  ALLOCATE(Veg_Coeff(nlinesVeg,ncolumnsVeg))
  !Read input file
  OPEN(23,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(23,*)   !Skip lines before header
  ENDDO
  READ(23,*) (HeaderVeg_File(iv),iv=1,ncolumnsVeg) !Get header

  DO i=1,nlinesVeg
     READ(23,*) (Veg_Coeff(i,iv),iv=1,ncolumnsVeg)
     !write(*,*) (Veg_Coeff(i,iv),iv=1,ncolumnsVeg)
  ENDDO
  CLOSE(23)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesVeg
     DO ii=1,nlinesVeg
        IF(Veg_Coeff(i,cp_Code)==Veg_Coeff(ii,cp_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Veg_Coeff(i,cp_Code),'in Pervious.txt not unique!'
           CALL ErrorHint(60,FileN,Veg_Coeff(i,cp_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !=======================SUEWS_Water.txt=================================
  FileN='SUEWS_Water.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesWater=nlines
  ALLOCATE(Water_Coeff(nlinesWater,ncolumnsWater))
  !Read input file
  OPEN(24,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(24,*)   !Skip lines before header
  ENDDO
  READ(24,*) (HeaderWater_File(iv),iv=1,ncolumnsWater) !Get header

  DO i=1,nlinesWater
     READ(24,*) (Water_Coeff(i,iv),iv=1,ncolumnsWater)
     !write(*,*) (Water_Coeff(i,iv),iv=1,ncolumnsWater)
  ENDDO
  CLOSE(24)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesWater
     DO ii=1,nlinesWater
        IF(Water_Coeff(i,cw_Code)==Water_Coeff(ii,cw_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Water_Coeff(i,cw_Code),'in Water.txt not unique!'
           CALL ErrorHint(60,FileN,Water_Coeff(i,cw_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO


  !=======================SUEWS_Snow.txt==================================
  FileN='SUEWS_Snow.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesSnow=nlines
  ALLOCATE(Snow_Coeff(nlinesSnow,ncolumnsSnow))
  !Read input file
  OPEN(25,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(25,*)   !Skip lines before header
  ENDDO
  READ(25,*) (HeaderSnow_File(iv),iv=1,ncolumnsSnow) !Get header

  DO i=1,nlinesSnow
     READ(25,*) (Snow_Coeff(i,iv),iv=1,ncolumnsSnow)
     !write(*,*) (Snow_Coeff(i,iv),iv=1,ncolumnsSnow)
  ENDDO
  CLOSE(25)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesSnow
     DO ii=1,nlinesSnow
        IF(Snow_Coeff(i,cs_Code)==Snow_Coeff(ii,cs_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Snow_Coeff(i,cs_Code),'in Snow.txt not unique!'
           CALL ErrorHint(60,FileN,Snow_Coeff(i,cs_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO


  !=======================SUEWS_Soil.txt==================================
  FileN='SUEWS_Soil.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesSoil=nlines
  ALLOCATE(Soil_Coeff(nlinesSoil,ncolumnsSoil))
  !Read input file
  OPEN(26,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(26,*)   !Skip lines before header
  ENDDO
  READ(26,*) (HeaderSoil_File(iv),iv=1,ncolumnsSoil) !Get header

  DO i=1,nlinesSoil
     READ(26,*) (Soil_Coeff(i,iv),iv=1,ncolumnsSoil)
     !write(*,*) (Soil_Coeff(i,iv),iv=1,ncolumnsSoil)
  ENDDO
  CLOSE(26)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesSoil
     DO ii=1,nlinesSoil
        IF(Soil_Coeff(i,cSo_Code)==Soil_Coeff(ii,cSo_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Soil_Coeff(i,cSo_Code),'in Soil.txt not unique!'
           CALL ErrorHint(60,FileN,Soil_Coeff(i,cSo_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !===================SUEWS_Conductance.txt===============================
  FileN='SUEWS_Conductance.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesConductance=nlines
  ALLOCATE(Conductance_Coeff(nlinesConductance,ncolumnsConductance))
  !Read input file
  OPEN(27,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(27,*)   !Skip lines before header
  ENDDO
  READ(27,*) (HeaderCond_File(iv),iv=1,ncolumnsConductance) !Get header

  DO i=1,nlinesConductance
     READ(27,*) (Conductance_Coeff(i,iv),iv=1,ncolumnsConductance)
     !write(*,*) (Conductance_Coeff(i,iv),iv=1,ncolumnsConductance)
  ENDDO
  CLOSE(27)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesConductance
     DO ii=1,nlinesConductance
        IF(Conductance_Coeff(i,cc_Code)==Conductance_Coeff(ii,cc_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Conductance_Coeff(i,cc_Code),'in Conductance.txt not unique!'
           CALL ErrorHint(60,FileN,Conductance_Coeff(i,cc_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !===================SUEWS_OHMCoefficients.txt===========================
  FileN='SUEWS_OHMCoefficients.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesOHMCoefficients=nlines
  ALLOCATE(OHMCoefficients_Coeff(nlinesOHMCoefficients,ncolumnsOHMCoefficients))
  !Read input file
  OPEN(28,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(28,*)   !Skip lines before header
  ENDDO
  READ(28,*) (HeaderOHMCoefficients_File(iv),iv=1,ncolumnsOHMCoefficients) !Get header

  DO i=1,nlinesOHMCoefficients
     READ(28,*) (OHMCoefficients_Coeff(i,iv),iv=1,ncolumnsOHMCoefficients)
     !write(*,*) (OHMCoefficients_Coeff(i,iv),iv=1,ncolumnsOHMCoefficients)
  ENDDO
  CLOSE(28)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesOHMCoefficients
     DO ii=1,nlinesOHMCoefficients
        IF(OHMCoefficients_Coeff(i,cO_Code)==OHMCoefficients_Coeff(ii,cO_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',OHMCoefficients_Coeff(i,cO_Code),'in OHMCoefficients.txt not unique!'
           CALL ErrorHint(60,FileN,OHMCoefficients_Coeff(i,cO_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !===================SUEWS_ESTMCoefficients.txt===========================
  FileN='SUEWS_ESTMCoefficients.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesESTMCoefficients=nlines
  ALLOCATE(ESTMCoefficients_Coeff(nlinesESTMCoefficients,ncolumnsESTMCoefficients))
  !Read input file
  OPEN(28,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(28,*)   !Skip lines before header
  ENDDO
  READ(28,*) (HeaderESTMCoefficients_File(iv),iv=1,ncolumnsESTMCoefficients) !Get header

  DO i=1,nlinesESTMCoefficients
     READ(28,*) (ESTMCoefficients_Coeff(i,iv),iv=1,ncolumnsESTMCoefficients)
  ENDDO
  CLOSE(28)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesESTMCoefficients
     DO ii=1,nlinesESTMCoefficients
        IF(ESTMCoefficients_Coeff(i,cE_Code)==ESTMCoefficients_Coeff(ii,cE_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',ESTMCoefficients_Coeff(i,cE_Code),'in ESTMCoefficients.txt not unique!'
           CALL ErrorHint(60,FileN,ESTMCoefficients_Coeff(i,cE_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !================SUEWS_AnthropogenicHeat.txt============================
  FileN='SUEWS_AnthropogenicHeat.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesAnthropogenicHeat=nlines
  ALLOCATE(AnthropogenicHeat_Coeff(nlinesAnthropogenicHeat,ncolumnsAnthropogenicHeat))
  !Read input file
  OPEN(29,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(29,*)   !Skip lines before header
  ENDDO
  READ(29,*) (HeaderAnthropogenicHeat_File(iv),iv=1,ncolumnsAnthropogenicHeat) !Get header

  DO i=1,nlinesAnthropogenicHeat
     READ(29,*) (AnthropogenicHeat_Coeff(i,iv),iv=1,ncolumnsAnthropogenicHeat)
     !write(*,*) (AnthropogenicHeat_Coeff(i,iv),iv=1,ncolumnsAnthropogenicHeat)
  ENDDO
  CLOSE(29)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesAnthropogenicHeat
     DO ii=1,nlinesAnthropogenicHeat
        IF(AnthropogenicHeat_Coeff(i,cA_Code)==AnthropogenicHeat_Coeff(ii,cA_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',AnthropogenicHeat_Coeff(i,cA_Code),'in AnthropogenicHeat.txt not unique!'
           CALL ErrorHint(60,FileN,AnthropogenicHeat_Coeff(i,cA_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !================SUEWS_Irrigation.txt===================================
  FileN='SUEWS_Irrigation.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesIrrigation=nlines
  ALLOCATE(Irrigation_Coeff(nlinesIrrigation,ncolumnsIrrigation))
  !Read input file
  OPEN(30,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(30,*)   !Skip lines before header
  ENDDO
  READ(30,*) (HeaderIrrigation_File(iv),iv=1,ncolumnsIrrigation) !Get header

  DO i=1,nlinesIrrigation
     READ(30,*) (Irrigation_Coeff(i,iv),iv=1,ncolumnsIrrigation)
     !write(*,*) (Irrigation_Coeff(i,iv),iv=1,ncolumnsIrrigation)
  ENDDO
  CLOSE(30)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesIrrigation
     DO ii=1,nlinesIrrigation
        IF(Irrigation_Coeff(i,cIr_Code)==Irrigation_Coeff(ii,cIr_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Irrigation_Coeff(i,cIr_Code),'in Irrigation.txt not unique!'
           CALL ErrorHint(60,FileN,Irrigation_Coeff(i,cIr_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !================SUEWS_Profiles.txt=====================================
  FileN='SUEWS_Profiles.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesProfiles=nlines
  ALLOCATE(Profiles_Coeff(nlinesProfiles,ncolumnsProfiles))
  !Read input file
  OPEN(31,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(31,*)   !Skip lines before header
  ENDDO
  READ(31,*) (HeaderProfiles_File(iv),iv=1,ncolumnsProfiles) !Get header

  DO i=1,nlinesProfiles
     READ(31,*) (Profiles_Coeff(i,iv),iv=1,ncolumnsProfiles)
     !write(*,*) (Profiles_Coeff(i,iv),iv=1,ncolumnsProfiles)
  ENDDO
  CLOSE(31)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesProfiles
     DO ii=1,nlinesProfiles
        IF(Profiles_Coeff(i,cPr_Code)==Profiles_Coeff(ii,cPr_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',Profiles_Coeff(i,cPr_Code),'in Profiles.txt not unique!'
           CALL ErrorHint(60,FileN,Profiles_Coeff(i,cPr_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !================SUEWS_WithinGridWaterDist.txt===========================
  FileN='SUEWS_WithinGridWaterDist.txt'
  CALL NumberRows(FileN,SkipHeaderSiteInfo)     !Find number of rows in input file
  nlinesWGWaterDist=nlines
  ALLOCATE(WGWaterDist_Coeff(nlinesWGWaterDist,ncolumnsWGWaterDist))
  !Read input file
  OPEN(32,file=TRIM(FileInputPath)//TRIM(FileN),err=300,status='old')
  DO SkipCounter=1,(SkipHeaderSiteInfo-1)
     READ(32,*)    !Skip lines before header
  ENDDO
  READ(32,*) (HeaderWGWaterDist_File(iv),iv=1,ncolumnsWGWaterDist) !Get header

  DO i=1,nlinesWGWaterDist
     READ(32,*) (WGWaterDist_Coeff(i,iv),iv=1,ncolumnsWGWaterDist)
     !write(*,*) (WGWaterDist_Coeff(i,iv),iv=1,ncolumnsWGWaterDist)
  ENDDO
  CLOSE(32)

  CALL InputHeaderCheck(FileN)

  ! Check codes are unique
  DO i=1,nlinesWGWaterDist
     DO ii=1,nlinesWGWaterDist
        IF(WGWaterDist_Coeff(i,cWG_Code)==WGWaterDist_Coeff(ii,cWG_Code) .AND. i/=ii) THEN
           WRITE(*,*) 'Code',WGWaterDist_Coeff(i,cWG_Code),'in WithinGridWaterDist.txt not unique!'
           CALL ErrorHint(60,FileN,WGWaterDist_Coeff(i,cWG_Code),notUsed,notUsedI)
        ENDIF
     ENDDO
  ENDDO

  !=======================================================================
  !=======================================================================

  !-----------------------------------------------------------------------
  !SUEWS run information
  InputMetFormat=10    !Input met data file in LUMPS format(1) or SUEWS format(10)
  LAICalcYes=1         !Use observed(0) or modelled(1) LAI
  ity=2                !Evaporation calculated according to Rutter(1) or Shuttleworth(2)
  WriteDailyState = 1  !Daily state file written
  tstepcount=0

  t_INTERVAL = 3600   !Number of seconds in an hour

  !Calculate nsh (number of steps per hour) from model timestep (tstep) set in in RunControl
  nsh_real = t_INTERVAL/REAL(tstep,KIND(1d0))

  ! Check nsh is an integer
  IF(nsh_real==INT(nsh_real)) THEN
     nsh = INT(nsh_real)
  ELSE
     CALL ErrorHint(39,'TSTEP must divide into t_INTERVAL exactly.',REAL(tstep,KIND(1d0)),REAL(t_INTERVAL,KIND(1d0)),notUsedI)
  ENDIF

  ! Check nsh is reasonable
  IF(nsh_real<6.OR.nsh_real>60) THEN
     CALL ErrorHint(39,'TSTEP is too small or too large.',REAL(tstep,KIND(1d0)),REAL(t_INTERVAL,KIND(1d0)),notUsedI)
  ENDIF

  ! Cast integer nsh as nsh_real for use in calculations
  nsh_real = REAL(nsh,KIND(1d0))
  ! Cast integer tstep as tstep_real for use in calculations
  tstep_real = REAL(tstep,KIND(1d0))

  !! Check this is still valid for v2016a
  HalfTimeStep=REAL(tstep_real)/2/(24*3600)   !Used in sun_position to get sunpos in the middle of timestep

  RETURN

  !-------Possible problems-----------------------------------------------
200 CALL ErrorHint(47,'RunControl.nml',notUsed,notUsed,notUsedI)
201 CALL ErrorHint(48,'RunControl.nml',notUsed,notUsed,notUsedI)

203 CALL ErrorHint(47,TRIM(FileChoices),notUsed,notUsed,notUsedI)

300 CALL ErrorHint(48,TRIM(FileN),notUsed,notUsed,notUsedI)
  !-----------------------------------------------------------------------

  !pause

END SUBROUTINE OverallRunControl
!=========================================================================



!=========================================================================
!This subroutine finds the number of rows in each input file
!INPUT: FileN          Name of the input file
!       SkipHeaderLines   Number of header rows to skip
!Made by LJ/HCW in Oct 2014

SUBROUTINE NumberRows(FileN,SkipHeaderLines)

  USE data_in
  USE DefaultNotUsed
  USE Initial

  IMPLICIT NONE

  CHARACTER(len=50):: FileN
  INTEGER:: SkipHeaderLines, RunNumber
  INTEGER:: SkipCounter

  WRITE(*,*) FileN
  OPEN(39,file=TRIM(FileInputPath)//TRIM(FileN),err=204,status='old')

  IF(SkipHeaderLines > 0) THEN
     DO SkipCounter=1,SkipHeaderLines
        READ(39,*,err=205)
        !write(*,*) SkipCounter, SkipHeaderLines
     ENDDO
  ENDIF

  nlines = 0 !Initialize nlines
  DO
     READ(39,*) RunNumber
     IF(RunNumber==-9) EXIT
     nlines = nlines + 1
  END DO
  !write(*,*) 'nlines read: ',nlines
  CLOSE(39)

  RETURN

204 CALL ErrorHint(47,TRIM(FileInputPath)//TRIM(FileN),notUsed,notUsed,notUsedI)
205 CALL ErrorHint(48,TRIM(FileInputPath)//TRIM(FileN),notUsed,notUsed,notUsedI)

END SUBROUTINE NumberRows
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
!  - Rename profiles 1-24 rather than 0-23?
!----------------------------------------------------------------------------------------------
SUBROUTINE InitializeSurfaceCharacteristics(Gridiv,rr)

  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE Initial
  USE sues_data

  IMPLICIT NONE

  INTEGER:: Gridiv,&    !Row of SurfaceChar where input information will be stored
       rr          !Row of SiteSelect that matches current grid and year

  !-------------------------------------------------------------------------------------------

  ! Initialise row of SurfaceChar
  SurfaceChar(Gridiv,:) = -999

  ! Transfer data in SiteSelect to SurfaceChar
  SurfaceChar(Gridiv,1:ncolumnsSiteSelect) = SiteSelect(rr,1:ncolumnsSiteSelect) !Cols in same order as in SiteSelect.txt

  ! ======== Retrieve information from other input files via codes ========

  ! ---- Find code for Paved surface (Impervious) ----
  CALL CodeMatchNonVeg(rr,c_PavedCode)
  ! Transfer characteristics to SurfaceChar for Paved surface
  SurfaceChar(gridiv,c_AlbMin(PavSurf))       = NonVeg_Coeff(iv5,ci_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(PavSurf))       = NonVeg_Coeff(iv5,ci_AlbMax)
  SurfaceChar(gridiv,c_Emis(PavSurf))         = NonVeg_Coeff(iv5,ci_Emis)
  SurfaceChar(gridiv,c_StorMin(PavSurf))      = NonVeg_Coeff(iv5,ci_StorMin)
  SurfaceChar(gridiv,c_StorMax(PavSurf))      = NonVeg_Coeff(iv5,ci_StorMax)
  SurfaceChar(gridiv,c_WetThresh(PavSurf))    = NonVeg_Coeff(iv5,ci_WetThresh)
  SurfaceChar(gridiv,c_StateLimit(PavSurf))   = NonVeg_Coeff(iv5,ci_StateLimit)
  SurfaceChar(gridiv,c_DrEq(PavSurf))         = NonVeg_Coeff(iv5,ci_DrEq)
  SurfaceChar(gridiv,c_DrCoef1(PavSurf))      = NonVeg_Coeff(iv5,ci_DrCoef1)
  SurfaceChar(gridiv,c_DrCoef2(PavSurf))      = NonVeg_Coeff(iv5,ci_DrCoef2)
  SurfaceChar(gridiv,c_SoilTCode(PavSurf))    = NonVeg_Coeff(iv5,ci_SoilTCode)
  SurfaceChar(gridiv,c_SnowLimPat(PavSurf))   = NonVeg_Coeff(iv5,ci_SnowLimPat)
  SurfaceChar(gridiv,c_SnowLimRem(PavSurf))   = NonVeg_Coeff(iv5,ci_SnowLimRem)
  SurfaceChar(gridiv,c_OHMCode_SWet(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SWet)
  SurfaceChar(gridiv,c_OHMCode_SDry(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_SDry)
  SurfaceChar(gridiv,c_OHMCode_WWet(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WWet)
  SurfaceChar(gridiv,c_OHMCode_WDry(PavSurf)) = NonVeg_Coeff(iv5,ci_OHMCode_WDry)
  SurfaceChar(Gridiv,c_CpAnOHM(PavSurf))           = NonVeg_Coeff(iv5,ci_CpAnOHM)   ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(PavSurf))           = NonVeg_Coeff(iv5,ci_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(PavSurf))           = NonVeg_Coeff(iv5,ci_ChAnOHM)  ! bulk transfer coef., AnOHM TS

  ! Use SoilCode for Paved to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(PavSurf))
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
  CALL CodeMatchOHM(Gridiv,PavSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,PavSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,PavSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,PavSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(PavSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  !Get ESTM parameters for ground
  ! Transfer ESTM characteristics to SurfaceChar
  CALL CodeMatchESTM(Gridiv,PavSurf,'Grnd')
  !roof
  SurfaceChar(Gridiv,c_thick1_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_r)
  SurfaceChar(Gridiv,c_k1_r(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_r)
  SurfaceChar(Gridiv,c_rhoCp1_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_r)
  SurfaceChar(Gridiv,c_thick2_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_r)
  SurfaceChar(Gridiv,c_k2_r(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_r)
  SurfaceChar(Gridiv,c_rhoCp2_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_r)
  SurfaceChar(Gridiv,c_thick3_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_r)
  SurfaceChar(Gridiv,c_k3_r(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_r)
  SurfaceChar(Gridiv,c_rhoCp3_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_r)
  SurfaceChar(Gridiv,c_thick4_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_r)
  SurfaceChar(Gridiv,c_k4_r(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_r)
  SurfaceChar(Gridiv,c_rhoCp4_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_r)
  SurfaceChar(Gridiv,c_thick5_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_r)
  SurfaceChar(Gridiv,c_k5_r(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_r)
  SurfaceChar(Gridiv,c_rhoCp5_r(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_r)   !parameters below are empty
  SurfaceChar(Gridiv,c_thick1_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_e)
  SurfaceChar(Gridiv,c_k1_e(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_e)
  SurfaceChar(Gridiv,c_rhoCp1_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_e)
  SurfaceChar(Gridiv,c_thick2_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_e)
  SurfaceChar(Gridiv,c_k2_e(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_e)
  SurfaceChar(Gridiv,c_rhoCp2_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_e)
  SurfaceChar(Gridiv,c_thick3_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_e)
  SurfaceChar(Gridiv,c_k3_e(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_e)
  SurfaceChar(Gridiv,c_rhoCp3_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_e)
  SurfaceChar(Gridiv,c_thick4_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_e)
  SurfaceChar(Gridiv,c_k4_e(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_e)
  SurfaceChar(Gridiv,c_rhoCp4_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_e)
  SurfaceChar(Gridiv,c_thick5_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_e)
  SurfaceChar(Gridiv,c_k5_e(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_e)
  SurfaceChar(Gridiv,c_rhoCp5_e(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_e)
  SurfaceChar(Gridiv,c_thick1_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_i)
  SurfaceChar(Gridiv,c_k1_i(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_i)
  SurfaceChar(Gridiv,c_rhoCp1_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_i)
  SurfaceChar(Gridiv,c_thick2_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_i)
  SurfaceChar(Gridiv,c_k2_i(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_i)
  SurfaceChar(Gridiv,c_rhoCp2_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_i)
  SurfaceChar(Gridiv,c_thick3_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_i)
  SurfaceChar(Gridiv,c_k3_i(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_i)
  SurfaceChar(Gridiv,c_rhoCp3_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_i)
  SurfaceChar(Gridiv,c_thick4_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_i)
  SurfaceChar(Gridiv,c_k4_i(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_i)
  SurfaceChar(Gridiv,c_rhoCp4_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_i)
  SurfaceChar(Gridiv,c_thick5_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_i)
  SurfaceChar(Gridiv,c_k5_i(PavSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_i)
  SurfaceChar(Gridiv,c_rhoCp5_i(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_i)
  SurfaceChar(Gridiv,c_nroom(PavSurf))       = ESTMCoefficients_Coeff(iv5,cE_nroom)
  SurfaceChar(Gridiv,c_alb_ibld(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_alb_ibld)
  SurfaceChar(Gridiv,c_em_ibld(PavSurf))     = ESTMCoefficients_Coeff(iv5,cE_em_ibld)
  SurfaceChar(Gridiv,c_CH_iwall(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_CH_iwall)
  SurfaceChar(Gridiv,c_CH_iroof(PavSurf))    = ESTMCoefficients_Coeff(iv5,cE_CH_iroof)
  SurfaceChar(Gridiv,c_CH_ibld(PavSurf))     = ESTMCoefficients_Coeff(iv5,cE_CH_ibld)
  SurfaceChar(Gridiv,c_fwall(PavSurf))       = ESTMCoefficients_Coeff(iv5,cE_fwall)

  ! Get water distribution (within grid) for Paved
  CALL CodeMatchDist(rr,c_WGPavedCode,cWG_ToPaved)
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
  CALL CodeMatchNonVeg(rr,c_BldgsCode)
  ! Transfer characteristics to SurfaceChar for Bldgs surface
  SurfaceChar(gridiv,c_AlbMin(BldgSurf))       = NonVeg_Coeff(iv5,ci_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(BldgSurf))       = NonVeg_Coeff(iv5,ci_AlbMax)
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
  SurfaceChar(Gridiv,c_CpAnOHM(BldgSurf))            = NonVeg_Coeff(iv5,ci_CpAnOHM)   ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(BldgSurf))            = NonVeg_Coeff(iv5,ci_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(BldgSurf))            = NonVeg_Coeff(iv5,ci_ChAnOHM)  ! bulk transfer coef., AnOHM TS

  ! Use SoilCode for Bldgs to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(BldgSurf))
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
  CALL CodeMatchOHM(Gridiv,BldgSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BldgSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BldgSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BldgSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(BldgSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  !Get ESTM parameters for Bldgs (roof, external wall and internal element)
  ! Transfer ESTM characteristics to SurfaceChar
  CALL CodeMatchESTM(Gridiv,BldgSurf,'Bldg')
  !roof
  SurfaceChar(Gridiv,c_thick1_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_r)
  SurfaceChar(Gridiv,c_k1_r(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_r)
  SurfaceChar(Gridiv,c_rhoCp1_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_r)
  SurfaceChar(Gridiv,c_thick2_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_r)
  SurfaceChar(Gridiv,c_k2_r(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_r)
  SurfaceChar(Gridiv,c_rhoCp2_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_r)
  SurfaceChar(Gridiv,c_thick3_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_r)
  SurfaceChar(Gridiv,c_k3_r(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_r)
  SurfaceChar(Gridiv,c_rhoCp3_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_r)
  SurfaceChar(Gridiv,c_thick4_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_r)
  SurfaceChar(Gridiv,c_k4_r(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_r)
  SurfaceChar(Gridiv,c_rhoCp4_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_r)
  SurfaceChar(Gridiv,c_thick5_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_r)
  SurfaceChar(Gridiv,c_k5_r(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_r)
  SurfaceChar(Gridiv,c_rhoCp5_r(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_r)
  !external wall
  SurfaceChar(Gridiv,c_thick1_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_e)
  SurfaceChar(Gridiv,c_k1_e(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_e)
  SurfaceChar(Gridiv,c_rhoCp1_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_e)
  SurfaceChar(Gridiv,c_thick2_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_e)
  SurfaceChar(Gridiv,c_k2_e(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_e)
  SurfaceChar(Gridiv,c_rhoCp2_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_e)
  SurfaceChar(Gridiv,c_thick3_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_e)
  SurfaceChar(Gridiv,c_k3_e(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_e)
  SurfaceChar(Gridiv,c_rhoCp3_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_e)
  SurfaceChar(Gridiv,c_thick4_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_e)
  SurfaceChar(Gridiv,c_k4_e(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_e)
  SurfaceChar(Gridiv,c_rhoCp4_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_e)
  SurfaceChar(Gridiv,c_thick5_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_e)
  SurfaceChar(Gridiv,c_k5_e(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_e)
  SurfaceChar(Gridiv,c_rhoCp5_e(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_e)
  !internal element in buildings
  SurfaceChar(Gridiv,c_thick1_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick1_i)
  SurfaceChar(Gridiv,c_k1_i(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k1_i)
  SurfaceChar(Gridiv,c_rhoCp1_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP1_i)
  SurfaceChar(Gridiv,c_thick2_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick2_i)
  SurfaceChar(Gridiv,c_k2_i(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k2_i)
  SurfaceChar(Gridiv,c_rhoCp2_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP2_i)
  SurfaceChar(Gridiv,c_thick3_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick3_i)
  SurfaceChar(Gridiv,c_k3_i(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k3_i)
  SurfaceChar(Gridiv,c_rhoCp3_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP3_i)
  SurfaceChar(Gridiv,c_thick4_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick4_i)
  SurfaceChar(Gridiv,c_k4_i(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k4_i)
  SurfaceChar(Gridiv,c_rhoCp4_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP4_i)
  SurfaceChar(Gridiv,c_thick5_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_thick5_i)
  SurfaceChar(Gridiv,c_k5_i(BldgSurf))        = ESTMCoefficients_Coeff(iv5,cE_k5_i)
  SurfaceChar(Gridiv,c_rhoCp5_i(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_rhoCP5_i)
  SurfaceChar(Gridiv,c_nroom(BldgSurf))       = ESTMCoefficients_Coeff(iv5,cE_nroom)
  SurfaceChar(Gridiv,c_alb_ibld(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_alb_ibld)
  SurfaceChar(Gridiv,c_em_ibld(BldgSurf))     = ESTMCoefficients_Coeff(iv5,cE_em_ibld)
  SurfaceChar(Gridiv,c_CH_iwall(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_CH_iwall)
  SurfaceChar(Gridiv,c_CH_iroof(BldgSurf))    = ESTMCoefficients_Coeff(iv5,cE_CH_iroof)
  SurfaceChar(Gridiv,c_CH_ibld(BldgSurf))     = ESTMCoefficients_Coeff(iv5,cE_CH_ibld)
  SurfaceChar(Gridiv,c_fwall(BldgSurf))       = ESTMCoefficients_Coeff(iv5,cE_fwall)

  ! Get water distribution (within grid) for Bldgs
  CALL CodeMatchDist(rr,c_WGBldgsCode,cWG_ToBldgs)
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
  CALL CodeMatchVeg(rr,c_EveTrCode)
  ! Transfer characteristics to SurfaceChar for EveTr surface
  ! All surfaces (1-nsurf)
  SurfaceChar(gridiv,c_AlbMin(ConifSurf))     = Veg_Coeff(iv5,cp_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(ConifSurf))     = Veg_Coeff(iv5,cp_AlbMax)
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
  ! AnOHM TS
  SurfaceChar(Gridiv,c_CpAnOHM(ConifSurf))           = Veg_Coeff(iv5,cp_CpAnOHM) ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(ConifSurf))           = Veg_Coeff(iv5,cp_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(ConifSurf))           = Veg_Coeff(iv5,cp_ChAnOHM)  ! bulk transfer coef., AnOHM TS

  ! Use SoilCode for EveTr to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(ConifSurf))
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
  CALL CodeMatchOHM(Gridiv,ConifSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,ConifSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,ConifSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,ConifSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(ConifSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! Get water distribution (within grid) for EveTr
  CALL CodeMatchDist(rr,c_WGEveTrCode,cWG_ToEveTr)
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
  CALL CodeMatchVeg(rr,c_DecTrCode)
  ! Transfer characteristics to SurfaceChar for DecTr surface
  ! All surfaces (1-nsurf)
  SurfaceChar(gridiv,c_AlbMin(DecidSurf))     = Veg_Coeff(iv5,cp_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(DecidSurf))     = Veg_Coeff(iv5,cp_AlbMax)
  SurfaceChar(gridiv,c_Emis(DecidSurf))       = Veg_Coeff(iv5,cp_Emis)
  SurfaceChar(gridiv,c_StorMin(DecidSurf))    = Veg_Coeff(iv5,cp_StorMin)
  SurfaceChar(gridiv,c_StorMax(DecidSurf))    = Veg_Coeff(iv5,cp_StorMax)
  SurfaceChar(gridiv,c_WetThresh(DecidSurf))  = Veg_Coeff(iv5,cp_WetThresh)
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
  ! AnOHM TS
  SurfaceChar(Gridiv,c_CpAnOHM(DecidSurf))           = Veg_Coeff(iv5,cp_CpAnOHM) ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(DecidSurf))           = Veg_Coeff(iv5,cp_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(DecidSurf))           = Veg_Coeff(iv5,cp_ChAnOHM)  ! bulk transfer coef., AnOHM TS
  ! Use SoilCode for DecTr to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(DecidSurf))
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
  CALL CodeMatchOHM(Gridiv,DecidSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,DecidSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,DecidSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,DecidSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(DecidSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! Get water distribution (within grid) for DecTr
  CALL CodeMatchDist(rr,c_WGDecTrCode,cWG_ToDecTr)
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
  CALL CodeMatchVeg(rr,c_GrassCode)
  ! Transfer characteristics to SurfaceChar for Grass surface
  ! All surfaces (1-nsurf)
  SurfaceChar(gridiv,c_AlbMin(GrassSurf))     = Veg_Coeff(iv5,cp_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(GrassSurf))     = Veg_Coeff(iv5,cp_AlbMax)
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
  ! AnOHM TS
  SurfaceChar(Gridiv,c_CpAnOHM(GrassSurf))           = Veg_Coeff(iv5,cp_CpAnOHM) ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(GrassSurf))           = Veg_Coeff(iv5,cp_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(GrassSurf))           = Veg_Coeff(iv5,cp_ChAnOHM)  ! bulk transfer coef., AnOHM TS
  ! Use SoilCode for Grass to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(GrassSurf))
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
  CALL CodeMatchOHM(Gridiv,GrassSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,GrassSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,GrassSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,GrassSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(GrassSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! Get water distribution (within grid) for Grass
  CALL CodeMatchDist(rr,c_WGGrassCode,cWG_ToGrass)
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
  CALL CodeMatchNonVeg(rr,c_BSoilCode)
  ! Transfer characteristics to SurfaceChar for BSoil surface
  ! All surfaces (1-nsurf)
  SurfaceChar(gridiv,c_AlbMin(BSoilSurf))     = NonVeg_Coeff(iv5,ci_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(BSoilSurf))     = NonVeg_Coeff(iv5,ci_AlbMax)
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
  SurfaceChar(Gridiv,c_CpAnOHM(BSoilSurf))           = NonVeg_Coeff(iv5,ci_CpAnOHM) ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(BSoilSurf))           = NonVeg_Coeff(iv5,ci_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(BSoilSurf))           = NonVeg_Coeff(iv5,ci_ChAnOHM)  ! bulk transfer coef., AnOHM TS
  ! Use SoilCode for BSoil to find code for soil characteristics
  CALL CodeMatchSoil(Gridiv,c_SoilTCode(BSoilSurf))
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
  CALL CodeMatchOHM(Gridiv,BSoilSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BSoilSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BSoilSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,BSoilSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(BSoilSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! Get water distribution (within grid) for Bare soil
  CALL CodeMatchDist(rr,c_WGBSoilCode,cWG_ToBSoil)
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
  CALL CodeMatchWater(rr,c_WaterCode)
  ! Transfer characteristics to SurfaceChar for Water surface
  ! All surfaces (1-nsurf)
  SurfaceChar(gridiv,c_AlbMin(WaterSurf))     = Water_Coeff(iv5,cw_AlbMin)
  SurfaceChar(gridiv,c_AlbMax(WaterSurf))     = Water_Coeff(iv5,cw_AlbMax)
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
  ! AnOHM TS
  SurfaceChar(Gridiv,c_CpAnOHM(WaterSurf))           = Water_Coeff(iv5,cw_CpAnOHM) ! heat capacity, AnOHM TS
  SurfaceChar(Gridiv,c_KkAnOHM(WaterSurf))           = Water_Coeff(iv5,cw_KkAnOHM)  ! heat conductivity, AnOHM TS
  SurfaceChar(Gridiv,c_ChAnOHM(WaterSurf))           = Water_Coeff(iv5,cw_ChAnOHM)  ! bulk transfer coef., AnOHM TS
  ! Get OHM characteristics for Water
  CALL CodeMatchOHM(Gridiv,WaterSurf,'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,WaterSurf,'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,WaterSurf,'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,WaterSurf,'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(WaterSurf))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! Get water distribution (within grid) for Water
  CALL CodeMatchDist(rr,c_WGWaterCode,cWG_ToWater)
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
  CALL CodeMatchSnow(rr,c_SnowCode)
  ! Transfer characteristics to SurfaceChar for Snow surface
  SurfaceChar(gridiv,c_SnowRMFactor) = Snow_Coeff(iv5,cs_SnowRMFactor)
  SurfaceChar(gridiv,c_SnowTMFactor) = Snow_Coeff(iv5,cs_SnowTMFactor)
  SurfaceChar(gridiv,c_SnowAlbMin)   = Snow_Coeff(iv5,cs_SnowAlbMin)
  SurfaceChar(gridiv,c_SnowAlbMax)   = Snow_Coeff(iv5,cs_SnowAlbMax)
  !SurfaceChar(gridiv,c_SnowAlb)      = Snow_Coeff(iv5,cs_SnowAlb)
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
  !    SurfaceChar(Gridiv,c_CpAnOHM(nsurf+1))           = Snow_Coeff(iv5,cs_CpAnOHM)   ! heat capacity, AnOHM TS
  !    SurfaceChar(Gridiv,c_KkAnOHM(nsurf+1))           = Snow_Coeff(iv5,cs_KkAnOHM)  ! heat conductivity, AnOHM TS
  !    SurfaceChar(Gridiv,c_ChAnOHM(nsurf+1))           = Snow_Coeff(iv5,cs_ChAnOHM)  ! bulk transfer coef., AnOHM TS
  ! Get OHM characteristics for Snow
  CALL CodeMatchOHM(Gridiv,(nsurf+1),'SWet')  !Summer wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,(nsurf+1),'SDry')  !Summer dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_SDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,(nsurf+1),'WWet')  !Winter wet
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WWet(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)
  CALL CodeMatchOHM(Gridiv,(nsurf+1),'WDry')  !Winter dry
  ! Transfer OHM characteristics to SurfaceChar
  SurfaceChar(Gridiv,c_a1_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a1)
  SurfaceChar(Gridiv,c_a2_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a2)
  SurfaceChar(Gridiv,c_a3_WDry(nsurf+1))    = OHMCoefficients_Coeff(iv5,cO_a3)

  ! ---- Find code for Surface conductances ----
  CALL CodeMatchConductance(rr,c_CondCode)
  ! Transfer conductance characteristics to SurfaceChar
  SurfaceChar(gridiv,c_GsG1)       = Conductance_Coeff(iv5,cc_GsG1)
  SurfaceChar(gridiv,c_GsG2)       = Conductance_Coeff(iv5,cc_GsG2)
  SurfaceChar(gridiv,c_GsG3)       = Conductance_Coeff(iv5,cc_GsG3)
  SurfaceChar(gridiv,c_GsG4)       = Conductance_Coeff(iv5,cc_GsG4)
  SurfaceChar(gridiv,c_GsG5)       = Conductance_Coeff(iv5,cc_GsG5)
  SurfaceChar(gridiv,c_GsG6)       = Conductance_Coeff(iv5,cc_GsG6)
  SurfaceChar(gridiv,c_GsTH)       = Conductance_Coeff(iv5,cc_GsTH)
  SurfaceChar(gridiv,c_GsTL)       = Conductance_Coeff(iv5,cc_GsTL)
  SurfaceChar(gridiv,c_GsS1)       = Conductance_Coeff(iv5,cc_GsS1)
  SurfaceChar(gridiv,c_GsS2)       = Conductance_Coeff(iv5,cc_GsS2)
  SurfaceChar(gridiv,c_GsKmax)     = Conductance_Coeff(iv5,cc_GsKmax)

  ! ---- Find code for Anthropogenic heat ----
  CALL CodeMatchAnthropogenicHeat(rr,c_QFCode)
  ! Transfer Anthropogenic heat characteristics to SurfaceChar
  SurfaceChar(gridiv,c_BaseTHDD)   = AnthropogenicHeat_Coeff(iv5,cA_BaseTHDD)
  SurfaceChar(gridiv,c_QF_A1)      = AnthropogenicHeat_Coeff(iv5,cA_QF_A1)
  SurfaceChar(gridiv,c_QF_B1)      = AnthropogenicHeat_Coeff(iv5,cA_QF_B1)
  SurfaceChar(gridiv,c_QF_C1)      = AnthropogenicHeat_Coeff(iv5,cA_QF_C1)
  SurfaceChar(gridiv,c_QF_A2)      = AnthropogenicHeat_Coeff(iv5,cA_QF_A2)
  SurfaceChar(gridiv,c_QF_B2)      = AnthropogenicHeat_Coeff(iv5,cA_QF_B2)
  SurfaceChar(gridiv,c_QF_C2)      = AnthropogenicHeat_Coeff(iv5,cA_QF_C2)
  SurfaceChar(gridiv,c_AHMin)      = AnthropogenicHeat_Coeff(iv5,cA_AHMin)
  SurfaceChar(gridiv,c_AHSlope)    = AnthropogenicHeat_Coeff(iv5,cA_AHSlope)
  SurfaceChar(gridiv,c_TCritic)    = AnthropogenicHeat_Coeff(iv5,cA_TCritic)

  ! ---- Find code for Irrigation ----
  CALL CodeMatchIrrigation(rr,c_IrrCode)
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
  CALL CodeMatchProf(rr,c_EnProfWD)
  SurfaceChar(gridiv,c_HrProfEnUseWD) = Profiles_Coeff(iv5,cPr_Hours)
  ! Energy use (weekends)
  CALL CodeMatchProf(rr,c_EnProfWE)
  SurfaceChar(gridiv,c_HrProfEnUseWE) = Profiles_Coeff(iv5,cPr_Hours)
  ! Water use profile (manual, weekdays)
  CALL CodeMatchProf(rr,c_WProfManuWD)
  SurfaceChar(gridiv,c_HrProfWUManuWD)  = Profiles_Coeff(iv5,cPr_Hours)
  ! Water use profile (manual, weekends)
  CALL CodeMatchProf(rr,c_WProfManuWE)
  SurfaceChar(gridiv,c_HrProfWUManuWE)  = Profiles_Coeff(iv5,cPr_Hours)
  ! Water use profile (automatic, weekdays)
  CALL CodeMatchProf(rr,c_WProfAutoWD)
  SurfaceChar(gridiv,c_HrProfWUAutoWD) = Profiles_Coeff(iv5,cPr_Hours)
  ! Water use profile (automatic, weekends)
  CALL CodeMatchProf(rr,c_WProfAutoWE)
  SurfaceChar(gridiv,c_HrProfWUAutoWE) = Profiles_Coeff(iv5,cPr_Hours)
  ! Snow clearing profile (weekdays)
  CALL CodeMatchProf(rr,c_SnowProfWD)
  SurfaceChar(gridiv,c_HrProfSnowCWD) = Profiles_Coeff(iv5,cPr_Hours)
  ! Snow clearing profile (weekends)
  CALL CodeMatchProf(rr,c_SnowProfWE)
  SurfaceChar(gridiv,c_HrProfSnowCWE) = Profiles_Coeff(iv5,cPr_Hours)

  ! ---- Interpolate Hourly Profiles to model timestep and normalise
  TstepProfiles(Gridiv,:,:) = -999   !Initialise TstepProfiles
  ! Energy use
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_EnUseWD,c_HrProfEnUseWD)
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_EnUseWE,c_HrProfEnUseWE)
  ! For energy use, normalise so the AVERAGE of the multipliers is equal to 1
  TstepProfiles(Gridiv,cTP_EnUseWD,:) = TstepProfiles(Gridiv,cTP_EnUseWD,:) / SUM(TstepProfiles(Gridiv,cTP_EnUseWD,:))*24*nsh_real
  TstepProfiles(Gridiv,cTP_EnUseWE,:) = TstepProfiles(Gridiv,cTP_EnUseWE,:) / SUM(TstepProfiles(Gridiv,cTP_EnUseWE,:))*24*nsh_real

  ! Water use
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUManuWD,c_HrProfWUManuWD)
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUManuWE,c_HrProfWUManuWE)
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUAutoWD,c_HrProfWUAutoWD)
  CALL SUEWS_InterpHourlyProfiles(Gridiv,cTP_WUAutoWE,c_HrProfWUAutoWE)
  ! For water use, normalise so the SUM of the multipliers is equal to 1 (profile is multiplied by daily water use)
  TstepProfiles(Gridiv,cTP_WUManuWD,:) = TstepProfiles(Gridiv,cTP_WUManuWD,:) / SUM(TstepProfiles(Gridiv,cTP_WUManuWD,:))
  TstepProfiles(Gridiv,cTP_WUManuWE,:) = TstepProfiles(Gridiv,cTP_WUManuWE,:) / SUM(TstepProfiles(Gridiv,cTP_WUManuWE,:))
  TstepProfiles(Gridiv,cTP_WUAutoWD,:) = TstepProfiles(Gridiv,cTP_WUAutoWD,:) / SUM(TstepProfiles(Gridiv,cTP_WUAutoWD,:))
  TstepProfiles(Gridiv,cTP_WUAutoWE,:) = TstepProfiles(Gridiv,cTP_WUAutoWE,:) / SUM(TstepProfiles(Gridiv,cTP_WUAutoWE,:))

END SUBROUTINE InitializeSurfaceCharacteristics


!----------------------------------------------------------------------------------------------
!Calculates the initial conditions for each grid on the first year of run
!Made by sg feb 2012 -
!Latest modified:
!20 Oct 2014, LJ: Saves to slot 0 of the output matrix
!06 Nov 2014 HCW
!----------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------
SUBROUTINE InitialState(GridName,year_int,Gridiv,year_txt)
  ! Last modified HCW 03 Jul 2015
  ! Added initial conditions albEveTr0 and albGrass0
  ! Last modified by HCW 03 Dec 2014
  ! To do:
  !   - Check running means (5-day temperature)
  !------------------------------------------------------------------------

  USE allocateArray
  USE data_in
  USE ColNamesModelDailyState
  USE defaultNotUsed
  USE FileName
  USE gis_data
  USE InitialCond
  USE mod_z
  USE resist
  USE snowMod
  USE sues_data
  USE time

  IMPLICIT NONE

  CHARACTER(len=20):: GridName    !Name of the evaluated grid
  !character(len=10):: str2        !Variables related to filepaths
  CHARACTER(len=150):: fileInit   !Initial conditions filename
  CHARACTER(len=4):: year_txt     !year in txt format
  INTEGER::DaysSinceRain,Gridiv,& !number of days since rain, grid number,
       gamma1,gamma2          !switches related to cooling and heating degree days
  INTEGER::wd,seas,date,mb,&      !weekday information, season, date, month
       year_int,switch=0,&    !year as an integer, switch related to previous day
       id_next,calc           !next day,counter in irrigation calculations

  REAL (KIND(1d0))::PavedState,BldgsState,EveTrState,DecTrState,GrassState,BSoilState,&
       SnowFracPaved,SnowFracBldgs,SnowFracEveTr,SnowFracDecTr,          &
       SnowFracGrass,SnowFracBSoil,SnowFracWater,                        &
       SnowDensPaved,SnowDensBldgs,SnowDensEveTr,SnowDensDecTr,          &
       SnowDensGrass,SnowDensBSoil,SnowDensWater

  !-----------------------------------------------------------------------

  NAMELIST/InitialConditions/DaysSinceRain,&
       Temp_C0,&
       ID_Prev,&
       GDD_1_0,&
       GDD_2_0,&
       LAIinitialEveTr,&
       LAIinitialDecTr,&
       LAIinitialGrass,&
       AlbEveTr0,&
       albDec0,&
       AlbGrass0,&
       decidCap0,&
       porosity0,&
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
       SnowDensWater,&
       SnowAlb0

  ! Define InitialConditions file ----------------------------------------
  FileInit=TRIM(FileInputPath)//TRIM("InitialConditions")//TRIM(GridName)//'.nml'

  ! Open, read and close InitialConditions file --------------------------
  OPEN(56,File=TRIM(FileInit),err=600,status='old') !Change with needs
  READ(56,iostat=ios_out,nml=InitialConditions,err=601)
  CLOSE(56)

  ! Write InitialConditions to FileChoices -------------------------------
  FileChoices=TRIM(FileOutputPath)//TRIM(FileCode)//'_FileChoices.txt'
  OPEN(12,file=FileChoices,position='append')
  WRITE(12,*)'----------',TRIM(FileInit),'----------'
  WRITE(12,nml=InitialConditions)
  CLOSE(12)

  !-----------------------------------------------------------------------

  ! Previous day DOY number (needed in file allocations)
  IF(id_prev>=364) id_prev=0  !If previous day is larger than 364, set this to zero

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
  ModelDailyState(Gridiv,cMDS_albEveTr)   = albEveTr0
  ModelDailyState(Gridiv,cMDS_albGrass)   = albGrass0

  ! -- Anthropogenic heat flux initializations --
  ! Need to get BaseTHDD from SurfaceChar, as info not transferred until SUEWS_Translate called
  BaseTHDD = SurfaceChar(Gridiv,c_BaseTHDD)

  IF(AnthropHeatChoice>=0) THEN
     !Calculations related to heating and cooling degree days (BaseT is used always)
     IF ((Temp_C0-BaseTHDD)>=0) THEN   !Cooling
        gamma2=1
     ELSE
        gamma2=0
     ENDIF
     IF ((BaseTHDD-Temp_C0)>=0) THEN   !Heating
        gamma1=1
     ELSE
        gamma1=0
     ENDIF
     ModelDailyState(Gridiv,cMDS_HDD1) = gamma1*(BaseTHDD-Temp_C0) ! Heating
     ModelDailyState(Gridiv,cMDS_HDD2) = gamma2*(Temp_C0-BaseTHDD) ! Cooling
  ENDIF

  ModelDailyState(Gridiv,cMDS_TempC) = Temp_C0
  ModelDailyState(Gridiv,cMDS_DaysSinceRain) = REAL(DaysSinceRain,KIND(1d0))

  ! Assume that the temperature has been the same for the previous days
  ModelDailyState(Gridiv,cMDS_TempCOld1) = Temp_C0
  ModelDailyState(Gridiv,cMDS_TempCOld2) = Temp_C0
  ModelDailyState(Gridiv,cMDS_TempCOld3) = Temp_C0

  !! These variables don't seem to be needed (commented out HCW 27 Nov 2014)
  !! If required, they will need updating for a non-hourly timestep
  !! Initialize hourly temperature and precipitation + 5 day mean used to define thermal growing season
  !runT(0:23)=Temp_C0
  !runP(0:23)=0

  finish=.FALSE.

  ! -- Save snow density and snow albedo info in InitialConditions to ModelDailyState array --
  ModelDailyState(Gridiv,cMDS_SnowDens(PavSurf))    = SnowDensPaved
  ModelDailyState(Gridiv,cMDS_SnowDens(BldgSurf))   = SnowDensBldgs
  ModelDailyState(Gridiv,cMDS_SnowDens(ConifSurf))  = SnowDensEveTr
  ModelDailyState(Gridiv,cMDS_SnowDens(DecidSurf))  = SnowDensDecTr
  ModelDailyState(Gridiv,cMDS_SnowDens(GrassSurf))  = SnowDensGrass
  ModelDailyState(Gridiv,cMDS_SnowDens(BSoilSurf))  = SnowDensBSoil
  ModelDailyState(Gridiv,cMDS_SnowDens(WaterSurf))  = SnowDensWater
  ModelDailyState(Gridiv,cMDS_SnowAlb)  = SnowAlb0

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
  CALL SUEWS_Translate(Gridiv,0,0)
  PRINT*, 'id_prev in initial after translate',id_prev

  !Calculation of roughness parameters (N.B. uses porosity)
  CALL RoughnessParameters


  !=============================================================================
  ! If the run start day is at previous year, then calculate the number of days
  ! in that year.

  !First we need to know if the previous day given in initial conditions (id_prev) is
  !on previous year as this is needed in the initialization of dayofWeek matrix.
  !In this case switch is set to one for date calculations.
  IF(id_prev==0)THEN                     !If id_prev = 0, means that the first modelled day is 1 Jan
     year_int=year_int-1                  !1) find the number of days on that year
     CALL LeapYearCalc (year_int,id_prev) !2) set switch to 1 so that the code knows to change back to current year (switch=0)
     switch=1
  ENDIF

  CALL day2month(id_prev,mb,date,seas,year_int,lat) !Calculate date information (mb = month, date = day,...)
  CALL Day_of_Week(date,mb,year_int,wd)             !Calculate weekday of the previous day (wd) (1=Sun, ..., 7=Sat)

  !After the day in previous year switch is changed back to zero:
  !    ie not previous day anymore
  !Also the size of dayofweek is from 0:NdaysinYear meaning
  !that in zero slot is the previous day information
  IF(switch==1)THEN
     year_int=year_int+1
     id_prev=0
     switch=0
  ENDIF

  dayofWeek(id_prev,1)=wd   ! day of week
  dayofWeek(id_prev,2)=mb   ! month
  dayofweek(id_prev,3)=seas ! season (summer=1, winter=2) needed for accumulation

  ! in case next day goes to next year calculate again the date information for dayofweek matrix.
  id_next=id_prev+1
  IF(id_next>nofDaysThisYear) THEN
     id_next=1
     year_int=year_int+1
     switch=1
     CALL ErrorHint(43,'switch- years',notUsed,notUsed,notUsedI)
  ENDIF

  CALL day2month(id_next,mb,date,seas,year_int,lat) !Calculate real date from doy
  CALL Day_of_Week(date,mb,year_int,wd)             !Calculate weekday (1=Sun, ..., 7=Sat)

  IF(switch==1)THEN
     iy=iy-1
     switch=0
  ENDIF

  dayofWeek(id_next,1)=wd  ! day of week
  dayofWeek(id_next,2)=mb  ! month
  dayofweek(id_next,3)=seas ! season

  !=============================================================================

  !id=id_prev
  !it= 23 !!LastTimeOfDay
  IF (id_prev>=DayLightSavingDay(1).AND.id_prev<=DayLightSavingDay(2)) THEN  !Summertime
     DLS=1
  ELSE
     DLS=0
  ENDIF

  ! -----------------------------------------------------------------------
  ! Calculate daily water use if modelled (i.e. if WU_choice = 0).
  ! Calculated from previous day information given in InitialConditions file

  WU_day=0                !Initialize WU_day
  IF (WU_choice==0) THEN  !Model water use
     calc=0

     IF (DayWat(wd)==1.0) THEN !if DayWat(wd)=1.0 (irrigation occurs on this day)
        IF (lat>=0) THEN            !Northern Hemisphere
           IF (id>=Ie_start.AND.id<=Ie_end) calc=1 !if day between irrigation period
        ELSE                        !Southern Hemisphere
           calc=1
           IF (id>=Ie_end.AND.id<=Ie_start) calc=0 !if day between irrigation period
        ENDIF
        IF(calc==1) THEN
           ! Model daily water use based on HDD(id,6)(days since rain) and HDD(id,3)(average temp)

           ! ---- Automatic irrigation (evergreen trees) ----
           WU_day(id,2) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(ConifSurf)*IrrFracConif*DayWatPer(wd)
           IF (WU_Day(id,2)<0) WU_Day(id,2)=0   !If modelled WU is negative -> 0

           ! ---- Manual irrigation (evergreen trees) ----
           WU_day(id,3) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(ConifSurf)*IrrFracConif*DayWatPer(wd)
           IF (WU_Day(id,3)<0) WU_Day(id,3)=0   !If modelled WU is negative -> 0

           ! ---- Total evergreen trees water use (automatic + manual) ----
           WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))

           ! ---- Automatic irrigation (deciduous trees) ----
           WU_day(id,5) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(DecidSurf)*IrrFracDecid*DayWatPer(wd)
           IF (WU_Day(id,5)<0) WU_Day(id,5)=0   !If modelled WU is negative -> 0

           ! ---- Manual irrigation (deciduous trees) ----
           WU_day(id,6) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(DecidSurf)*&
                IrrFracDecid*DayWatPer(wd)
           IF (WU_Day(id,6)<0) WU_Day(id,6)=0   !If modelled WU is negative -> 0

           ! ---- Total deciduous trees water use (automatic + manual) ----
           WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))

           ! ---- Automatic irrigation (grass) ----
           WU_day(id,8) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(GrassSurf)*&
                IrrFracGrass*DayWatPer(wd)
           IF (WU_Day(id,8)<0) WU_Day(id,8)=0   !If modelled WU is negative -> 0
           ! ---- Manual irrigation (grass) ----
           WU_day(id,9) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(GrassSurf)*&
                IrrFracGrass*DayWatPer(wd)
           IF (WU_Day(id,9)<0) WU_Day(id,9)=0   !If modelled WU is negative -> 0
           ! ---- Total grass water use (automatic + manual) ----
           WU_Day(id,7)=(WU_day(id,8)+WU_day(id,9))
        ELSE
           WU_Day(id,1)=0
           WU_Day(id,2)=0
           WU_Day(id,3)=0
           WU_Day(id,4)=0
           WU_Day(id,5)=0
           WU_Day(id,6)=0
           WU_Day(id,7)=0
           WU_Day(id,8)=0
           WU_Day(id,9)=0
        ENDIF
     ENDIF
  ENDIF

  ! -----------------------------------------------------------------------

  !Initialise rates of change variables for OHM calculation
  ! check ?? this should not happen at the start of the year if continuing

  !q1_grids(Gridiv)=-20
  !q2_grids(Gridiv)=-22
  !q3_grids(Gridiv)=-24

  !r1_grids(Gridiv)=-20
  !r2_grids(Gridiv)=-22
  !r3_grids(Gridiv)=-24

  q1_grids(Gridiv)=-101
  q2_grids(Gridiv)=-100
  q3_grids(Gridiv)=-99

  r1_grids(Gridiv)=-101
  r2_grids(Gridiv)=-100
  r3_grids(Gridiv)=-99

  ! ---- AnOHM TS ---------------------
  ! initialize Bowen ratio
  Bo_grids(0,:)=1.
  mAH_grids(0,:)=25.

  ! -----------------------------------

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !DEFINE DIFFERENT INITIALIZATION PARAMETERS
  !Are still needed ??
  !once=.true.
  RETURN

600 CALL ErrorHint(47,TRIM(FileInit),notUsed,notUsed,notUsedI)
601 CALL ErrorHint(48,TRIM(FileInit),notUsed,notUsed,ios_out)

END SUBROUTINE InitialState

!-------------------------------------------------------------------------
SUBROUTINE NextInitial(GridName,year_int)
  ! Last modified LJ 06 Jul 2015
  ! Initial conditions of SnowAlb added, densSnow changed to SnowDens
  ! Last modified HCW 03 Jul 2015
  ! Added initial conditions albEveTr0 and albGrass0
  ! Modified by HCW 21 Nov 2014
  ! Last day of year is not anymore the number of days on that year, but rather
  ! id == 1. Thus nofDaysThisYear was changed to 1. LJ 9/4/2015
  !------------------------------------------------------------------------
  USE allocateArray
  USE ColNamesInputFiles
  USE ColNamesModelDailyState
  USE data_in
  USE defaultNotUsed
  USE Initial
  USE sues_data
  USE snowMod
  USE time

  IMPLICIT NONE

  CHARACTER (len=15)::GridName
  CHARACTER (len=4)::year_txt2
  INTEGER:: year_int2
  INTEGER:: year_int

  year=year_int   !HCW added 21 Nov 2014

  IF (id==1) THEN  !nofDaysThisYear changed to 1
     year_int2=INT(year+1)
     WRITE(year_txt2,'(I4)')year_int2
     OPEN(57,File=TRIM(FileInputPath)//TRIM("InitialConditions")//TRIM(GridName)//'.nml',err=200)
  ELSE
     year_int2=INT(year)
     WRITE(year_txt2,'(I4)')year_int2
     OPEN(57,File=TRIM(FileInputPath)//TRIM("InitialConditions")//TRIM(GridName)//'end.nml',err=201)
  ENDIF
  !endif

  WRITE(57,*)'&InitialConditions'
  WRITE(57,*)'DaysSinceRain=',INT(HDD(id,6))

  !! If last time of day, then DailyState variables will have been updated so can write out arrays for id rather than id-1
  !if(it==23 .and. imin == (nsh_real-1)/nsh_real*60) then  !!LastTimeofday
  !   id=id+1
  !endif

  WRITE(57,*)'Temp_C0=',HDD(nofDaysThisYear,3)
  WRITE(57,*)'ID_Prev=',nofDaysThisYear
  WRITE(57,*)'GDD_1_0=',GDD(nofDaysThisYear,1)
  WRITE(57,*)'GDD_2_0=',GDD(nofDaysThisYear,2)
  WRITE(57,*)'PavedState=',State(PavSurf)
  WRITE(57,*)'BldgsState=',State(BldgSurf)
  WRITE(57,*)'EveTrState=',State(ConifSurf)
  WRITE(57,*)'DecTrState=',State(DecidSurf)
  WRITE(57,*)'GrassState=',State(GrassSurf)
  WRITE(57,*)'BSoilState=',State(BSoilSurf)
  WRITE(57,*)'WaterState=',State(WaterSurf)
  WRITE(57,*)'LAIinitialEveTr=',lai(nofDaysThisYear,ivConif)
  WRITE(57,*)'LAIinitialDecTr=',lai(nofDaysThisYear,ivDecid)
  WRITE(57,*)'LAIinitialGrass=',lai(nofDaysThisYear,ivGrass)
  WRITE(57,*)'albEveTr0=',AlbEveTr(id)
  WRITE(57,*)'albDec0=',AlbDec(id)
  WRITE(57,*)'albGrass0=',AlbGrass(id)
  WRITE(57,*)'DecidCap0=',decidCap(id)
  WRITE(57,*)'porosity0=',porosity(id)
  WRITE(57,*)'soilstorePavedState=',soilmoist(PavSurf)
  WRITE(57,*)'soilstoreBldgsState=',soilmoist(BldgSurf)
  WRITE(57,*)'soilstoreEveTrState=',soilmoist(ConifSurf)
  WRITE(57,*)'soilstoreDecTrState=',soilmoist(DecidSurf)
  WRITE(57,*)'soilstoreGrassState=',soilmoist(GrassSurf)
  WRITE(57,*)'soilstoreBSoilState=',soilmoist(BSoilSurf)
  WRITE(57,*)'SnowWaterPavedstate=',MeltWaterStore(PavSurf)
  WRITE(57,*)'SnowWaterBldgsState=',MeltWaterStore(BldgSurf)
  WRITE(57,*)'SnowWaterEveTrState=',MeltWaterStore(ConifSurf)
  WRITE(57,*)'SnowWaterDecTrState=',MeltWaterStore(DecidSurf)
  WRITE(57,*)'SnowWaterGrassState=',MeltWaterStore(GrassSurf)
  WRITE(57,*)'SnowWaterBSoilState=',MeltWaterStore(BSoilSurf)
  WRITE(57,*)'SnowWaterWaterState=',MeltWaterStore(WaterSurf)
  WRITE(57,*)'SnowPackPaved=',SnowPack(PavSurf)
  WRITE(57,*)'SnowPackBldgs=',SnowPack(BldgSurf)
  WRITE(57,*)'SnowPackEveTr=',SnowPack(ConifSurf)
  WRITE(57,*)'SnowPackDecTr=',SnowPack(DecidSurf)
  WRITE(57,*)'SnowPackGrass=',SnowPack(GrassSurf)
  WRITE(57,*)'SnowPackBSoil=',SnowPack(BSoilSurf)
  WRITE(57,*)'SnowPackWater=',SnowPack(WaterSurf)
  WRITE(57,*)'SnowFracPaved=',SnowFrac(PavSurf)
  WRITE(57,*)'SnowFracBldgs=',SnowFrac(BldgSurf)
  WRITE(57,*)'SnowFracEveTr=',SnowFrac(ConifSurf)
  WRITE(57,*)'SnowFracDecTr=',SnowFrac(DecidSurf)
  WRITE(57,*)'SnowFracGrass=',SnowFrac(GrassSurf)
  WRITE(57,*)'SnowFracBSoil=',SnowFrac(BSoilSurf)
  WRITE(57,*)'SnowFracWater=',SnowFrac(WaterSurf)
  WRITE(57,*)'SnowDensPaved=',SnowDens(PavSurf)
  WRITE(57,*)'SnowDensBldgs=',SnowDens(BldgSurf)
  WRITE(57,*)'SnowDensEveTr=',SnowDens(ConifSurf)
  WRITE(57,*)'SnowDensDecTr=',SnowDens(DecidSurf)
  WRITE(57,*)'SnowDensGrass=',SnowDens(GrassSurf)
  WRITE(57,*)'SnowDensBSoil=',SnowDens(BSoilSurf)
  WRITE(57,*)'SnowDensWater=',SnowDens(WaterSurf)
  WRITE(57,*)'SnowAlb0=',SnowAlb
  WRITE(57,*)'/'
  CLOSE(57)

  IF(it==23 .AND. imin == (nsh_real-1)/nsh_real*60) THEN
     id=id-1
  ENDIF

  RETURN

200 CALL ErrorHint(49,TRIM("InitialConditions")//TRIM(GridName)// &
       '_'//TRIM(ADJUSTL(year_txt2))//'.nml',notUsed,notUsed,notUsedI)
201 CALL ErrorHint(49,TRIM("InitialConditions")//TRIM(GridName)// &
       '_'//TRIM(ADJUSTL(year_txt2))//'end.nml',notUsed,notUsed,notUsedI)

END SUBROUTINE NextInitial
!-------------------------------------------------------------------------



!=======================================================================
!=======================================================================
!This subroutine prepares a meteorological forcing file.

SUBROUTINE SUEWS_InitializeMetData(lunit)

  USE allocateArray
  USE data_in
  USE sues_data
  USE time
  USE defaultnotUsed
  USE Initial

  IMPLICIT NONE

  INTEGER::lunit,i,iyy !,RunNumber,NSHcounter
  REAL (KIND(1d0)),DIMENSION(24)::MetArray
  REAL(KIND(1d0)):: imin_prev, ih_prev, iday_prev, tstep_met   !For checks on temporal resolution of met data

  !---------------------------------------------------------------

  !Open the file for reading and read the actual data
  !write(*,*) fileMet
  OPEN(lunit,file=TRIM(fileMet),status='old',err=314)
  CALL skipHeader(lunit,SkipHeaderMet)

  ! Skip to the right place in the met file, depending on how many chunks have been read already
  IF (skippedLines>0) THEN
     DO iyy=1,skippedLines
        READ(lunit,*)
     ENDDO
  ENDIF

  ! Read in next chunk of met data and fill MetForcingData array with data for every timestep
  !NSHcounter = 1
  DO i=1,ReadlinesMetdata
     CALL MetRead(MetArray,InputmetFormat,ldown_option,NetRadiationChoice,&
          snowUse,smd_choice,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)
     !DO iv=1,NSH
     !    MetForcingData(NSHcounter,1:24,GridCounter) = MetArray
     !   NSHcounter = NSHcounter + 1
     !ENDDO
     MetForcingData(i,1:24,GridCounter) = MetArray
     ! Check timestamp of met data file matches TSTEP specified in RunControl
     IF(i==1) THEN
        imin_prev = MetArray(4)
        ih_prev   = MetArray(3)
        iday_prev   = MetArray(2)
     ELSEIF(i==2) THEN
        tstep_met = ((MetArray(4)+60*MetArray(3)) - (imin_prev+60*ih_prev))*60   !tstep in seconds
        IF(tstep_met.NE.tstep_real.AND.MetArray(2)==iday_prev) THEN
           CALL ErrorHint(39,'TSTEP in RunControl does not match TSTEP of met data (DOY).',REAL(tstep,KIND(1d0)),tstep_met,&
                INT(MetArray(2)))
        ENDIF
     ENDIF

  ENDDO

  CLOSE(lunit)

  RETURN

314 CALL errorHint(11,TRIM(fileMet),notUsed,notUsed,ios_out)


END SUBROUTINE SUEWS_InitializeMetData
!----------------------------------------------------------------------------------------------


!====================================================================================
SUBROUTINE CheckInitial
  !Check the parameters in InitialConditions file.
  !Modified by HCW 04 Mar 2014, changed soilmoist(is) checks to use names given in InitialConditions
  !Added by LJ in 8/2/2013

  USE allocateArray
  USE data_in
  USE defaultNotUsed
  USE InitialCond
  USE snowMod
  USE time

  IMPLICIT NONE

  !real(kind(1d0)):: pTol   !Precision tolerance for range checks

  IF (Temp_C0<(Temp_C-10).OR.Temp_C0>(Temp_C+10)) THEN
     CALL ErrorHint(36,'InitialCond: Check temperature', Temp_C0, Temp_C, notUsedI)
  ENDIF

  IF (ID_Prev/=id-1) THEN
     CALL ErrorHint(36,'InitialCond: Check previous day', REAL(ID_Prev,KIND(1d0)), REAL(id,KIND(1d0)), notUsedI)
  ENDIF

!!!This part currently does not work for multiple grids as Initial conditions values get overwritten.
  !! Simple checks that Initial Conditions are within the specified ranges (within a precision tolerance)
  !pTol = 0.00001
  !!write(*,*) LAIInitialEveTr,LAImin(ConifSurf-2),LAImax(ConifSurf-2)
  !if(LAIInitialEveTr < (LAImin(ConifSurf-2)-pTol)) then
  !  call ErrorHint(36,'Intial LAI for EveTr < min value in SUEWS_Veg.txt!', LAIMin(ConifSurf-2), LAIInitialEveTr, notUsedI)
  !endif
  !if(LAIInitialEveTr > (LAImax(ConifSurf-2)+pTol)) then
  !   call ErrorHint(36,'Intial LAI for EveTr > max value in SUEWS_Veg.txt!', LAIMax(ConifSurf-2), LAIInitialEveTr, notUsedI)
  !endif
  !!write(*,*) LAIInitialDecTr,LAImin(DecidSurf-2)
  !if(LAIInitialDecTr < (LAImin(DecidSurf-2)-pTol)) then
  !   call ErrorHint(36,'Intial LAI for DecTr < min value in SUEWS_Veg.txt!', LAIMin(DecidSurf-2), LAIInitialDecTr, notUsedI)
  !endif
  !if(LAIInitialDecTr > (LAImax(DecidSurf-2)+pTol)) then
  !   call ErrorHint(36,'Intial LAI for DecTr > max value in SUEWS_Veg.txt!', LAIMax(DecidSurf-2), LAIInitialDecTr, notUsedI)
  !endif
  !if(LAIInitialGrass < (LAImin(GrassSurf-2)-pTol)) then
  !   call ErrorHint(36,'Intial LAI for Grass < min value in SUEWS_Veg.txt!', LAIMin(GrassSurf-2), LAIInitialGrass, notUsedI)
  !endif
  !if(LAIInitialGrass > (LAImax(GrassSurf-2)+pTol)) then
  !  call ErrorHint(36,'Intial LAI for Grass > max value in SUEWS_Veg.txt!', LAIMax(GrassSurf-2), LAIInitialGrass, notUsedI)
  !endif

  !!write(*,*) AlbDec0, albmin_DecTr, albmax_DecTr
  !if(AlbDec0 < (AlbMin_DecTr-pTol)) then
  !   !call ErrorHint(36,'Intial albedo for DecTr < min value in SUEWS_Veg.txt!', AlbMin_DecTr, AlbDec0, notUsedI)
  !   call ErrorHint(36,'Intial albedo for DecTr < min value!', AlbMin_DecTr, AlbDec0, notUsedI)
  !endif
  !if(AlbDec0 > (AlbMax_dec+pTol)) then
  !   !call ErrorHint(36,'Intial albedo for DecTr > max value in SUEWS_Veg.txt!', AlbMax_dec, AlbDec0, notUsedI)
  !   call ErrorHint(36,'Intial albedo for DecTr < max value!', AlbMax_dec, AlbDec0, notUsedI)
  !endif
  !!DecidCap0, Porosity0...

  !Check more thoroughly if LAI values are OK. Need to treat different hemispheres as well as
  !tropics separately.
  IF (lat>40) THEN
     IF ((LAIinitialEveTr>LAImin(ConifSurf-2)+1.AND.(id<60.OR.id>330)).OR.&
          (LAIinitialEveTr<LAImax(ConifSurf-2)-1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImin(ConifSurf), notUsedI)
     ENDIF
     IF ((LAIinitialDecTr>LAImin(DecidSurf-2)+1.AND.(id<60.OR.id>330)).OR.&
          (LAIinitialDecTr<LAImax(DecidSurf-2)-1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImin(DecidSurf), notUsedI)
     ENDIF
     IF ((LAIinitialGrass>LAImin(GrassSurf-2)+1.AND.(id<60.OR.id>330)).OR.&
          (LAIinitialGrass<LAImax(GrassSurf-2)-1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImin(GrassSurf), notUsedI)
     ENDIF

  ELSEIF (lat<-40) THEN
     IF ((LAIinitialEveTr<LAImax(ConifSurf-2)-1.AND.(id<60.OR.id>330)).OR.&
          (LAIinitialEveTr>LAImin(ConifSurf-2)+1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImax(ConifSurf), notUsedI)
     ENDIF
     IF ((LAIinitialDecTr>LAImax(DecidSurf-2)-1.AND.(id<60.OR.id>330)).OR.&
          (LAIinitialDecTr>LAImin(DecidSurf-2)+1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImax(DecidSurf), notUsedI)
     ENDIF
     IF ((LAIinitialGrass<LAImax(GrassSurf-2)-1.AND.(id<60.OR.id>330)) .OR.&
          (LAIinitialGrass>LAImin(GrassSurf-2)+1.AND.(id>130.AND.id<244))) THEN
        CALL ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImax(GrassSurf), notUsedI)
     ENDIF

  ELSEIF (lat<10.AND.lat>-10) THEN

     IF (LAIinitialEveTr<LAImax(ConifSurf-2)-0.5) THEN
        CALL ErrorHint(37,'Check LAIinitialEveTr in InitialConditions file', LAIinitialEveTr, LAImax(ConifSurf), notUsedI)
     ENDIF
     IF (LAIinitialDecTr<LAImax(DecidSurf-2)-0.5) THEN
        CALL ErrorHint(37,'Check LAIinitialDecTr in InitialConditions file', LAIinitialDecTr, LAImax(DecidSurf), notUsedI)
     ENDIF
     IF (LAIinitialGrass<LAImax(GrassSurf-2)-0.5) THEN
        CALL ErrorHint(37,'Check LAIinitialGrass in InitialConditions file', LAIinitialGrass, LAImax(GrassSurf), notUsedI)
     ENDIF

  ENDIF

  !Soilstore check
  IF (SoilstoreBldgsState>soilstoreCap(BldgSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of building soil store.',&
          SoilstoreBldgsState, soilstoreCap(BldgSurf), notUsedI)
  ENDIF
  IF (SoilstorePavedState>soilstoreCap(PavSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of paved soil store.',&
          SoilstorePavedState, soilstoreCap(PavSurf), notUsedI)
  ENDIF
  IF (SoilstoreEveTrState>soilstoreCap(ConifSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of conif soil store.',&
          SoilstoreEveTrState, soilstoreCap(ConifSurf), notUsedI)
  ENDIF
  IF (SoilstoreDecTrState>soilstoreCap(DecidSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of deciduous soil store.',&
          SoilstoreDecTrState, soilstoreCap(DecidSurf), notUsedI)
  ENDIF
  IF (SoilstoreBSoilState>soilstoreCap(BSoilSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of bare soil soil store.',&
          SoilstoreBSoilState, soilstoreCap(BSoilSurf), notUsedI)
  ENDIF
  IF (SoilstoreGrassState>soilstoreCap(GrassSurf)) THEN
     CALL ErrorHint(36,'InitialCond: Check initial condition of grass soil store.',&
          SoilstoreGrassState, soilstoreCap(GrassSurf), notUsedI)
  ENDIF

  !Snow stuff
  IF (snowUse==1) THEN
     IF (SnowWaterBldgsState>CRWmax*SnowPackBldgs) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterBldgsState', SnowWaterBldgsState, SnowPackBldgs, notUsedI)
     ENDIF
     IF (SnowWaterPavedState>CRWmax*SnowPackPaved) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterPavedState', SnowWaterPavedState, SnowPackPaved, notUsedI)
     ENDIF
     IF (SnowWaterEveTrState>CRWmax*SnowPackEveTr) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterEveTrstate', SnowWaterEveTrstate, SnowPackEveTr, notUsedI)
     ENDIF
     IF (SnowWaterDecTrState>CRWmax*SnowPackDecTr) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterDecTrState', SnowWaterDecTrState, SnowPackDecTr, notUsedI)
     ENDIF
     IF (SnowWaterGrassState>CRWmax*SnowPackGrass) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterGrassState', SnowWaterGrassState, SnowPackGrass, notUsedI)
     ENDIF
     IF (SnowWaterBSoilState>CRWmax*SnowPackBSoil) THEN
        CALL ErrorHint(36,'InitialCond: SnowWaterGrassUnirState', SnowWaterBSoilState, SnowPackBSoil, notUsedI)
     ENDIF
  ENDIF

  !SnowWaterWaterstate,& ??
  !SnowPackWater,& ??


END SUBROUTINE CheckInitial
