!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!Main program of SUEWS version 1.0  -------------------------------------------
! Model calculations are made in two stages:
! (1) initialise the run for each block of met data (iblock from 1 to ReadBlocksMetData)
! (2) perform the actual model calculations (SUEWS_Calculations)
! After reading in all input information from SiteInfo, the code loops
!  - over years
!  - then over blocks (met data read in and stored for each block)
!  - then over rows
!  - then over grids
!Last modified by TS 10 Apr 2017  - conditional compilation blocks added for netCDF adaptation
!Last modified by LJ 6 Apr 2017   - Snow initialisation, allocation and deallocation added
!Last modified by HCW 10 Feb 2017 - Disaggregation of met forcing data
!Last modified by HCW 12 Jan 2017 - Changes to InitialConditions
!Last modified by HCW 26 Aug 2016 - CO2 flux added
!Last modified by HCW 04 Jul 2016 - GridID can now be up to 10 digits long
!Last modified by HCW 29 Jun 2016 - Reversed over-ruling of ReadLinesMetdata so this is not restricted here to one day
!Last modified by HCW 27 Jun 2016 - Re-corrected grid number for output files. N.B. Gridiv seems to have been renamed iGrid
!                                 - Met file no longer has grid number attached if same met data used for all grids
!Last modified by HCW 24 May 2016 - InitialConditions file naming altered
!                                   Unused year_txt argument removed from InitialState
!                 LJ  30 Mar 2016 - Grid run order changed from linear to non-linear
!Last modified by TS 14 Mar 2016  - Include AnOHM daily iteration
!Last modified by HCW 25 Jun 2015 - Fixed bug in LAI calculation at year change
!Last modified by HCW 12 Mar 2015
!Last modified by HCW 26 Feb 2015
!Last modified by HCW 03 Dec 2014
!
! To do:
!   - Water movement between grids (GridConnections) not yet coded
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
PROGRAM SUEWS_Program

   USE AllocateArray
   USE ColNamesInputFiles
   USE Data_in
   USE DefaultNotUsed
   USE Initial
   USE MetDisagg
   USE Sues_Data
   USE Time
   USE WhereWhen
   USE ctrl_output
   USE ESTM_module, ONLY: SUEWS_GetESTMData, ESTM_initials
   USE BLUEWS_module, ONLY: CBL_ReadInputData

   IMPLICIT NONE

   CHARACTER(len=4) :: year_txt  !Year as a text string

   CHARACTER(len=20):: FileCodeX, & !Current file code
                       FileCodeXWY, & !File code without year
                       FileCodeXWG !File code without grid
   CHARACTER(len=20):: grid_txt, &    !Grid number as a text string (from FirstGrid to LastGrid)
                       tstep_txt, &   !Model timestep (in minutes) as a text string (minutes)
                       ResIn_txt, ResInESTM_txt     !Resolution of Original met/ESTM forcing file as a text string (minutes)

   INTEGER:: nlinesLimit, &   !Max number of lines that can be read in one go for each grid
             NumberOfYears   !Number of years to be run

   INTEGER:: year_int, & ! Year as an integer (from SiteSelect rather than met forcing file)
             igrid, &     !Grid number (from 1 to NumberOfGrids)
             iblock, &    !Block number (from 1 to ReadBlocksMetData)
             ir, irMax, &  !Row number within each block (from 1 to irMax)
             rr !Row of SiteSelect corresponding to current year and grid

   REAL::  timeStart, timeFinish ! profiling use, AnOHM TS

   !==========================================================================

   ! Start counting cpu time
   CALL CPU_TIME(timeStart)

   ! initialise simulation time
   dt_since_start = 0

   WRITE (*, *) '========================================================'
   WRITE (*, *) 'Running ', progname

   ! Initialise error file (0 -> problems.txt file will be newly created)
   errorChoice = 0
   ! Initialise error file (0 -> warnings.txt file will be newly created)
   warningChoice = 0
   ! Initialise OutputFormats to 1 so that output format is written out only once per run
   OutputFormats = 1

   ! Initialise WhereWhen variables for error handling
   GridID_text = '00000'
   datetime = '00000'

   ! Read RunControl.nml and all .txt input files from SiteSelect spreadsheet
   CALL overallRunControl

   WRITE (tstep_txt, '(I5)') tstep/60  !Get tstep (in minutes) as a text string
   WRITE (ResIn_txt, '(I5)') ResolutionFilesIn/60  !Get ResolutionFilesIn (in minutes) as a text string
   WRITE (ResInESTM_txt, '(I5)') ResolutionFilesInESTM/60

   ! Find first and last year of the current run
   FirstYear = MINVAL(INT(SiteSelect(:, c_Year)))
   LastYear = MAXVAL(INT(SiteSelect(:, c_Year)))

   NumberOfYears = LastYear - FirstYear + 1 !Find the number of years to run

   !Find the the number of grids within each year in SUEWS_SiteSelect.txt
   ! N.B. need to have the same grids for each year
   NumberOfGrids = INT(nlinesSiteSelect/NumberOfYears)

   !! Find the first and last grid numbers (N.B. need to have the same grids for each year)
   !FirstGrid = minval(int(SiteSelect(:,c_Grid)))
   !LastGrid  = maxval(int(SiteSelect(:,c_Grid)))
   IF (NumberOfGrids > MaxNumberOfGrids) THEN
      CALL ErrorHint(64, &
                     'No. of grids exceeds max. possible no. of grids.', &
                     REAL(MaxNumberOfGrids, KIND(1d0)), NotUsed, NumberOfGrids)
   ENDIF

   ALLOCATE (GridIDmatrix(NumberOfGrids)) !Get the nGrid numbers correctly
   ALLOCATE (GridIDmatrix0(NumberOfGrids)) !Get the nGrid numbers correctly

   DO igrid = 1, NumberOfGrids
      GridIDmatrix(igrid) = INT(SiteSelect(igrid, c_Grid))
   ENDDO

#ifdef nc
   ! sort grid matrix to conform the geospatial layout as in QGIS, TS 14 Dec 2016
   IF (ncMode == 1) THEN
      GridIDmatrix0 = GridIDmatrix
      CALL sortGrid(GridIDmatrix0, GridIDmatrix, nRow, nCol)
   ENDIF
   ! GridIDmatrix0 stores the grid ID in the original order
#endif

   ! GridIDmatrix=GridIDmatrix0
   WRITE (*, *) '--------------------------------------------'
   WRITE (*, *) 'Years identified:', FirstYear, 'to', LastYear
   WRITE (*, *) 'No. grids identified:', NumberOfGrids, 'grids'
   WRITE (*, *) 'Maximum No. grids allowed:', MaxNumberOfGrids, 'grids'

   ! Set limit on number of lines to read
   nlinesLimit = INT(FLOOR(MaxLinesMet/REAL(NumberOfGrids, KIND(1d0))))  !Uncommented HCW 29 Jun 2016
   !nlinesLimit = 24*nsh  !Commented out HCW 29 Jun 2016

   ! ---- Allocate arrays ----------------------------------------------------
   ! Daily state needs to be outside year loop to transfer states between years
   ALLOCATE (ModelDailyState(NumberOfGrids, MaxNCols_cMDS))   !DailyState
   ALLOCATE (DailyStateFirstOpen(NumberOfGrids))             !Initialisation for header

   ! ---- Initialise arrays --------------------------------------------------
   ModelDailyState(:, :) = -999
   DailyStateFirstOpen(:) = 1

   ! -------------------------------------------------------------------------
   ! Initialise ESTM (reads ESTM nml, should only run once)
   IF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
      IF (Diagnose == 1) WRITE (*, *) 'Calling ESTM_initials...'
      CALL ESTM_initials
   ENDIF

   !==========================================================================
   DO year_int = FirstYear, LastYear   !Loop through years

      WRITE (*, *) ' '
      WRITE (year_txt, '(I4)') year_int  !Get year as a text string

      ! Find number of days in the current year
      CALL LeapYearCalc(year_int, nofDaysThisYear)

      ! Prepare to disaggregate met data to model time-step (if required) ------
      ! Find number of model time-steps per resolution of original met forcing file
      Nper_real = ResolutionFilesIn/REAL(Tstep, KIND(1d0))
      Nper = INT(Nper_real)

      IF (Nper /= Nper_real) THEN
         CALL ErrorHint(2, 'Problem in SUEWS_Program: check resolution of met forcing data (ResolutionFilesIn)'// &
                        'and model time-step (Tstep).', &
                        REAL(Tstep, KIND(1d0)), NotUsed, ResolutionFilesIn)
      ELSEIF (Nper > 1) THEN
         WRITE (*, *) 'Resolution of met forcing data: ', TRIM(ADJUSTL(ResIn_txt)), ' min;', &
            ' model time-step: ', TRIM(ADJUSTL(tstep_txt)), ' min', ' -> SUEWS will perform disaggregation.'
         IF (Diagnose == 1) WRITE (*, *) 'Getting information for met disaggregation'
         ! Get names of original met forcing file(s) to disaggregate (using first grid)
         WRITE (grid_txt, '(I10)') GridIDmatrix(1)  !Get grid as a text string

         ! Get met file name for this grid: SSss_YYYY_data_RR.txt
         FileOrigMet = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)//'_data_' &
                       //TRIM(ADJUSTL(ResIn_txt))//'.txt'
         ! But if each grid has the same met file, met file name does not include grid number
         IF (MultipleMetFiles /= 1) THEN
            FileOrigMet = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_data_' &
                          //TRIM(ADJUSTL(ResIn_txt))//'.txt'
         ENDIF

         nlinesOrigMetdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
         nlinesOrigMetdata = count_lines(TRIM(FileOrigMet))
         ! WRITE(*,*) 'nlinesOrigMetdata', nlinesOrigMetdata

         ReadLinesOrigMetData = nlinesOrigMetdata   !Initially set limit as the size of  file
         IF (nlinesOrigMetData*Nper > nlinesLimit) THEN   !But restrict if this limit exceeds memory capacity
            ReadLinesOrigMetData = INT(nlinesLimit/Nper)
         ENDIF
         ! make sure the metblocks read in consists of complete diurnal cycles
         nsdorig = nsd/Nper
         ReadLinesOrigMetData = INT(MAX(nsdorig*(ReadLinesOrigMetData/nsdorig), nsdorig))
         !WRITE(*,*) 'ReadlinesOrigMetdata', ReadlinesOrigMetdata
         WRITE (*, *) 'Original met data will be read in chunks of ', ReadlinesOrigMetdata, 'lines.'

         ReadBlocksOrigMetData = INT(CEILING(REAL(nlinesOrigMetData, KIND(1d0))/REAL(ReadLinesOrigMetData, KIND(1d0))))

         ! Set ReadLinesMetdata and ReadBlocksMetData
         ReadLinesMetdata = ReadLinesOrigMetdata*Nper
         ReadBlocksMetData = INT(CEILING(REAL(nlinesOrigMetData*Nper, KIND(1d0))/REAL(ReadLinesMetdata, KIND(1d0))))
         WRITE (*, *) 'Processing current year in ', ReadBlocksMetData, 'blocks.'

         nlinesMetdata = nlinesOrigMetdata*Nper

      ELSEIF (Nper == 1) THEN
         WRITE (*, *) 'ResolutionFilesIn = Tstep: no disaggregation needed for met data.'

         !-----------------------------------------------------------------------
         ! Find number of lines in met forcing file for current year (nlinesMetdata)
         ! Need to know how many lines will be read each iteration
         ! Use first grid as an example as the number of lines is the same for all grids
         ! within one year
         WRITE (grid_txt, '(I10)') GridIDmatrix(1)  !Get grid as a text string

         ! Get met file name for this year for this grid
         FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
         FileMet = TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
         !If each grid has the same met file, met file name does not include grid number (HCW 27 Jun 2016)
         IF (MultipleMetFiles /= 1) THEN
            FileCodeXWG = TRIM(FileCode)//'_'//TRIM(year_txt) !File code without grid
            FileMet = TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
         ENDIF

         nlinesMetdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
         nlinesMetdata = count_lines(TRIM(FileMet))
         !-----------------------------------------------------------------------

         ! To conserve memory, read met data in blocks
         ! Find number of lines that can be read in each block (i.e. read in at once)
         ReadLinesMetdata = nlinesMetData   !Initially set limit as the size of the met file (N.B.solves problem with Intel fortran)
         IF (nlinesMetData > nlinesLimit) THEN   !But restrict if this limit exceeds memory capacity
            ReadLinesMetdata = nlinesLimit
         ENDIF
         ! make sure the metblocks read in consists of complete diurnal cycles, TS 08 Jul 2016
         ReadLinesMetdata = INT(MAX(nsd*(ReadLinesMetdata/nsd), nsd))

         WRITE (*, *) 'Met data will be read in blocks of ', ReadLinesMetdata, 'lines.'

         ! Find number of blocks of met data
         ReadBlocksMetData = INT(CEILING(REAL(nlinesMetData, KIND(1d0))/REAL(ReadLinesMetdata, KIND(1d0))))
         WRITE (*, *) 'Processing current year in ', ReadBlocksMetData, 'blocks.'

      ENDIF

      !write(*,*) ReadBlocksMetData, ReadBlocksOrigMetData

      ! ---- Allocate arrays--------------------------------------------------
      IF (Diagnose == 1) WRITE (*, *) 'Allocating arrays in SUEWS_Program.f95...'
      ALLOCATE (SurfaceChar(NumberOfGrids, MaxNCols_c))                                   !Surface characteristics
      ALLOCATE (MetForcingData(ReadLinesMetdata, ncolumnsMetForcingData, NumberOfGrids))   !Met forcing data
      ALLOCATE (MetForcingData_grid(ReadLinesMetdata, ncolumnsMetForcingData))   !Met forcing data
      ALLOCATE (ModelOutputData(0:ReadLinesMetdata, MaxNCols_cMOD, NumberOfGrids))         !Data at model timestep
      ALLOCATE (dataOutSUEWS(ReadLinesMetdata, ncolumnsDataOutSUEWS, NumberOfGrids))                 !Main output array
      dataOutSUEWS = NaN ! initialise Main output array
      ALLOCATE (dataOutDailyState(ndays, ncolumnsDataOutDailyState, NumberOfGrids))                 !DailyState array
      dataOutDailyState = NaN ! initialise DailyState
      IF (SOLWEIGuse == 1) ALLOCATE (dataOutSOL(ReadLinesMetdata, ncolumnsdataOutSOL, NumberOfGrids))     !SOLWEIG POI output
      IF (CBLuse >= 1) ALLOCATE (dataOutBL(ReadLinesMetdata, ncolumnsdataOutBL, NumberOfGrids))       !CBL output
      ! IF (SnowUse == 1) THEN
      IF (.NOT. ALLOCATED(dataOutSnow)) ALLOCATE (dataOutSnow(ReadLinesMetdata, ncolumnsDataOutSnow, NumberOfGrids))   !Snow output
      ! ALLOCATE(qn1_S_store(NSH,NumberOfGrids))
      ! ALLOCATE(qn1_S_store_grid(NSH))
      ! ALLOCATE(qn1_S_av_store(2*NSH+1,NumberOfGrids))
      ! ALLOCATE(qn1_S_av_store_grid(2*NSH+1))
      IF (.NOT. ALLOCATED(qn1_s_av_grids)) ALLOCATE (qn1_s_av_grids(NumberOfGrids))
      IF (.NOT. ALLOCATED(dqnsdt_grids)) ALLOCATE (dqnsdt_grids(NumberOfGrids))
      ! qn1_S_store(:,:) = NAN
      ! qn1_S_av_store(:,:) = NaN
      ! qn1_S_store_grid(:) = NAN
      ! qn1_S_av_store_grid(:) = NaN
      qn1_s_av_grids = 0 ! Initialise to 0
      dqnsdt_grids = 0 ! Initialise to 0
      ! ENDIF
      ! IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
      IF (.NOT. ALLOCATED(dataOutESTM)) ALLOCATE (dataOutESTM(ReadlinesMetdata, 32, NumberOfGrids)) !ESTM output
      ! ENDIF
      ! ALLOCATE(TstepProfiles(NumberOfGrids,14,24*NSH))   !Hourly profiles interpolated to model timestep
      ! ALLOCATE(AHProf_tstep(24*NSH,2))                   !Anthropogenic heat profiles at model timestep
      ! ALLOCATE(WUProfM_tstep(24*NSH,2))                  !Manual water use profiles at model timestep
      ! ALLOCATE(WUProfA_tstep(24*NSH,2))                  !Automatic water use profiles at model timestep
      ! ALLOCATE(HumActivity_tstep(24*NSH,2))              !Human activity profiles at model timestep
      ! ALLOCATE(TraffProf_tstep(24*NSH,2))                !Traffic profiles at model timestep
      ! ALLOCATE(PopProf_tstep(24*NSH,2))                  !Population profiles at model timestep
      ! ALLOCATE(qn1_store(NSH,NumberOfGrids))
      ! ALLOCATE(qn1_av_store(2*NSH+1,NumberOfGrids))
      ! ALLOCATE(qn1_store_grid(NSH))
      ! ALLOCATE(qn1_av_store_grid(2*NSH+1))
      ! ALLOCATE(qhforCBL(NumberOfGrids))
      ! ALLOCATE(qeforCBL(NumberOfGrids))
      ALLOCATE (qn1_av_grids(NumberOfGrids))
      ALLOCATE (dqndt_grids(NumberOfGrids))
      !! QUESTION: Add snow clearing (?)

      ! qn1_store(:,:) = NAN ! Initialise to -999
      ! qn1_av_store(:,:) = NAN ! Initialise to -999
      ! qn1_store_grid(:)=NAN
      ! qn1_av_store_grid(:)=NAN
      qn1_av_grids = 0 ! Initialise to 0
      dqndt_grids = 0 ! Initialise to 0

      ! qhforCBL(:) = NAN
      ! qeforCBL(:) = NAN

      ! QUESTION: Initialise other arrays here?

      IF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
         ! Prepare to disaggregate ESTM data to model time-step (if required) ------
         ! Find number of model time-steps per resolution of original met forcing file
         NperESTM_real = ResolutionFilesInESTM/REAL(Tstep, KIND(1d0))
         NperESTM = INT(NperESTM_real)
         IF (NperESTM /= NperESTM_real) THEN
            CALL ErrorHint(2, 'Problem in SUEWS_Program: check resolution of ESTM forcing data (ResolutionFilesInESTM)'// &
                           'and model time-step (Tstep).', &
                           REAL(Tstep, KIND(1d0)), NotUsed, ResolutionFilesInESTM)
         ELSEIF (NperESTM > 1) THEN
            WRITE (*, *) 'Resolution of ESTM forcing data: ', TRIM(ADJUSTL(ResInESTM_txt)), ' min;', &
               ' model time-step: ', TRIM(ADJUSTL(tstep_txt)), ' min', ' -> SUEWS will perform disaggregation.'
            IF (Diagnose == 1) WRITE (*, *) 'Getting information for ESTM disaggregation'
            ! Get names of original met forcing file(s) to disaggregate (using first grid)
            WRITE (grid_txt, '(I10)') GridIDmatrix(1)  !Get grid as a text string

            ! Get met file name for this grid
            FileESTMTs = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)//'_ESTM_Ts_data_' &
                         //TRIM(ADJUSTL(ResInESTM_txt))//'.txt'
            ! But if each grid has the same ESTM file, file name does not include grid number
            IF (MultipleESTMFiles /= 1) THEN
               FileESTMTs = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_ESTM_Ts_data_' &
                            //TRIM(ADJUSTL(ResInESTM_txt))//'.txt'
            ENDIF

            ! Find number of lines in orig ESTM file
            nlinesOrigESTMdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
            nlinesOrigESTMdata = count_lines(TRIM(FileESTMTs))

            ! Check ESTM data and met data will have the same length (so that ESTM file can be read in same blocks as met data)
            IF (nlinesOrigESTMdata*NperESTM /= nlinesMetData) THEN
               CALL ErrorHint(66, &
                              'Downscaled ESTM and met input files will have different lengths', REAL(nlinesMetdata, KIND(1d0)), &
                              NotUsed, nlinesESTMdata*NperESTM)
            ENDIF

            !write(*,*) 'nlinesOrigESTMdata', nlinesOrigESTMdata
            ! Set number of lines to read from original ESTM file using met data blocks
            ReadLinesOrigESTMData = ReadLinesMetdata/NperESTM
            !WRITE(*,*) 'ReadlinesOrigESTMdata', ReadlinesOrigESTMdata
            WRITE (*, *) 'Original ESTM data will be read in chunks of ', ReadlinesOrigESTMdata, 'lines.'

            nlinesESTMdata = nlinesOrigESTMdata*NperESTM

         ELSEIF (NperESTM == 1) THEN
            WRITE (*, *) 'ResolutionFilesInESTM = Tstep: no disaggregation needed for met data.'

            !-----------------------------------------------------------------------
            ! Find number of lines in ESTM forcing file for current year (nlinesESTMdata)
            WRITE (grid_txt, '(I10)') GridIDmatrix(1)  !Get grid as a text string (use first grid as example)
            ! Get file name for this year for this grid
            FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
            FileESTMTs = TRIM(FileInputPath)//TRIM(FileCodeX)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
            !If each grid has the same met file, met file name does not include grid number
            IF (MultipleESTMFiles /= 1) THEN
               FileCodeXWG = TRIM(FileCode)//'_'//TRIM(year_txt) !File code without grid
               FileESTMTs = TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
            ENDIF

            ! Find number of lines in ESTM file
            nlinesESTMdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
            nlinesESTMdata = count_lines(TRIM(FileESTMTs))
            !-----------------------------------------------------------------------

            ! Check ESTM data and met data are same length (so that ESTM file can be read in same blocks as met data)
            IF (nlinesESTMdata /= nlinesMetdata) THEN
               CALL ErrorHint(66, 'ESTM input file different length to met forcing file', REAL(nlinesMetdata, KIND(1d0)), &
                              NotUsed, nlinesESTMdata)
            ENDIF

         ENDIF

         ! Allocate arrays to receive ESTM forcing data
         ! IF (.NOT. ALLOCATED(ESTMForcingData)) ALLOCATE(ESTMForcingData(1:ReadLinesMetdata,ncolsESTMdata,NumberOfGrids))
         ! IF (.NOT. ALLOCATED(Ts5mindata)) ALLOCATE(Ts5mindata(1:ReadLinesMetdata,ncolsESTMdata))
         ! IF (.NOT. ALLOCATED(Ts5mindata_ir)) ALLOCATE(Ts5mindata_ir(ncolsESTMdata))
         !
         ! IF (.NOT. ALLOCATED(Tair24HR)) ALLOCATE(Tair24HR(24*nsh))
         ! if ( /= 0) print *, ": Deallocation request denied"
         ! ALLOCATE(Tair24HR(24*nsh))

      ENDIF
      IF (.NOT. ALLOCATED(ESTMForcingData)) ALLOCATE (ESTMForcingData(1:ReadLinesMetdata, ncolsESTMdata, NumberOfGrids))
      IF (.NOT. ALLOCATED(Ts5mindata)) ALLOCATE (Ts5mindata(1:ReadLinesMetdata, ncolsESTMdata))
      IF (.NOT. ALLOCATED(Ts5mindata_ir)) ALLOCATE (Ts5mindata_ir(ncolsESTMdata))
      IF (.NOT. ALLOCATED(Tair24HR)) ALLOCATE (Tair24HR(24*nsh))
      ! ------------------------------------------------------------------------

      !-----------------------------------------------------------------------
      SkippedLines = 0  !Initialise lines to be skipped in met forcing file
      SkippedLinesOrig = 0  !Initialise lines to be skipped in original met forcing file
      SkippedLinesOrigESTM = 0  !Initialise lines to be skipped in original met forcing file

      DO iblock = 1, ReadBlocksMetData   !Loop through blocks of met data

         ! WRITE(*,*) iblock,'/',ReadBlocksMetData

         ! Model calculations are made in two stages:
         ! (1) initialise the run for each block of met data (iblock from 1 to ReadBlocksMetData)
         ! (2) perform the actual model calculations (SUEWS_Calculations)

         GridCounter = 1   !Initialise counter for grids in each year
         DO igrid = 1, NumberOfGrids   !Loop through grids

            GridID = GridIDmatrix(igrid)   !store grid here for referencing error codes
            WRITE (grid_txt, '(I10)') GridIDmatrix(igrid)   !Get grid ID as a text string

            ! (1) First stage: initialise run if this is the first iteration this year
            ! (1a) Initialise surface characteristics
            IF (iblock == 1) THEN
               IF (Diagnose == 1) WRITE (*, *) 'First block of data - doing initialisation'
               ! (a) Transfer characteristics from SiteSelect to correct row of SurfaceChar
               DO rr = 1, nlinesSiteSelect
                  !Find correct grid and year
                  IF (Diagnose == 1) WRITE (*, *) 'grid found:', SiteSelect(rr, c_Grid), 'grid needed:', GridIDmatrix(igrid)
                  IF (Diagnose == 1) WRITE (*, *) 'year found:', SiteSelect(rr, c_Year), 'year needed:', year_int
                  IF (SiteSelect(rr, c_Grid) == GridIDmatrix(igrid) .AND. SiteSelect(rr, c_Year) == year_int) THEN
                     IF (Diagnose == 1) WRITE (*, *) 'Match found (grid and year) for rr = ', rr, 'of', nlinesSiteSelect
                     CALL InitializeSurfaceCharacteristics(GridCounter, rr)
                     EXIT
                  ELSEIF (rr == nlinesSiteSelect) THEN
                     WRITE (*, *) 'Program stopped! Year', year_int, 'and/or grid', igrid, 'not found in SiteSelect.txt.'
                     CALL ErrorHint(59, 'Cannot find year and/or grid in SiteSelect.txt', REAL(igrid, KIND(1d0)), NotUsed, year_int)
                  ENDIF
               ENDDO
            ENDIF   !end first block of met data

            ! (1b) Initialise met data
            IF (Nper > 1) THEN
               ! Disaggregate met data ---------------------------------------------------

               ! Set maximum value for ReadLinesOrigMetData to handle end of file (i.e. small final block)
               IF (iblock == ReadBlocksMetData) THEN   !For last block of data in file
                  ReadLinesOrigMetDataMAX = nlinesOrigMetdata - (iblock - 1)*ReadLinesOrigMetdata
               ELSE
                  ReadLinesOrigMetDataMAX = ReadLinesOrigMetData
               ENDIF
               !write(*,*) ReadLinesOrigMetDataMAX, ReadLinesOrigMetData
               ! Get names of original met forcing file(s) to disaggregate
               ! Get met file name for this grid: SSss_YYYY_data_RR.txt
               IF (MultipleMetFiles == 1) THEN   !If each grid has its own met file
                  FileOrigMet = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)//'_data_' &
                                //TRIM(ADJUSTL(ResIn_txt))//'.txt'
                  ! Also set file name for downscaled file
                  FileDscdMet = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)//'_data_' &
                                //TRIM(ADJUSTL(tstep_txt))//'.txt'
                  ! Disaggregate met data
                  CALL DisaggregateMet(iblock, igrid)
               ELSE
                  ! If each grid has the same met file, met file name does not include grid number, and only need to disaggregate once
                  FileOrigMet = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_data_' &
                                //TRIM(ADJUSTL(ResIn_txt))//'.txt'
                  FileDscdMet = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_data_' &
                                //TRIM(ADJUSTL(tstep_txt))//'.txt'
                  IF (igrid == 1) THEN       !Disaggregate for the first grid only
                     CALL DisaggregateMet(iblock, igrid)
                  ELSE                    !Then for subsequent grids simply copy data
                     MetForcingData(1:ReadLinesMetdata, 1:24, GridCounter) = MetForcingData(1:ReadLinesMetdata, 1:24, 1)
                  ENDIF
               ENDIF

            ELSEIF (Nper == 1) THEN
               ! Get met forcing file name for this year for the first grid
               ! Can be something else than 1
               FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
               FileCodeXWG = TRIM(FileCode)//'_'//TRIM(year_txt) !File code without grid
               !  IF(iblock==1) WRITE(*,*) 'Current FileCode: ', FileCodeX

               ! For every block of met data ------------------------------------
               ! Initialise met forcing data into 3-dimensional matrix
               !write(*,*) 'Initialising met data for block',iblock
               IF (MultipleMetFiles == 1) THEN   !If each grid has its own met file
                  FileMet = TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                  CALL SUEWS_InitializeMetData(1)
               ELSE                             !If one met file used for all grids
                  !FileMet=TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                  ! If one met file used for all grids, look for met file with no grid code (FileCodeXWG)
                  FileMet = TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                  IF (igrid == 1) THEN       !Read for the first grid only
                     CALL SUEWS_InitializeMetData(1)
                  ELSE                          !Then for subsequent grids simply copy data
                     MetForcingData(1:ReadLinesMetdata, 1:24, GridCounter) = MetForcingData(1:ReadLinesMetdata, 1:24, 1)
                  ENDIF
               ENDIF
            ENDIF   !end of nper statement

            ! Only for the first block of met data, read initial conditions (moved from above, HCW 12 Jan 2017)
            IF (iblock == 1) THEN
               FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
               !write(*,*) ' Now calling InitialState'
               CALL InitialState(FileCodeX, year_int, GridCounter, NumberOfGrids)
            ENDIF

            ! Initialise ESTM if required, TS 05 Jun 2016; moved inside grid loop HCW 27 Jun 2016
            IF (StorageHeatMethod == 4 .OR. StorageHeatMethod == 14) THEN
               IF (NperESTM > 1) THEN
                  ! Disaggregate ESTM data --------------------------------------------------
                  ! Set maximum value for ReadLinesOrigESTMData to handle end of file (i.e. small final block)
                  IF (iblock == ReadBlocksMetData) THEN   !For last block of data in file
                     ReadLinesOrigESTMDataMAX = nlinesOrigESTMdata - (iblock - 1)*ReadLinesOrigESTMdata
                  ELSE
                     ReadLinesOrigESTMDataMAX = ReadLinesOrigESTMData
                  ENDIF
                  !write(*,*) ReadLinesOrigESTMDataMAX, ReadLinesOrigESTMData
                  ! Get names of original ESTM forcing file(s) to disaggregate
                  ! Get ESTM file name for this grid: SSss_YYYY_ESTM_Ts_data_RR.txt
                  IF (MultipleESTMFiles == 1) THEN   !If each grid has its own ESTM file
                     FileOrigESTM = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt) &
                                    //'_ESTM_Ts_data_'//TRIM(ADJUSTL(ResInESTM_txt))//'.txt'
                     ! Also set file name for downscaled file
                     FileDscdESTM = TRIM(FileInputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt) &
                                    //'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                     ! Disaggregate ESTM data
                     CALL DisaggregateESTM(iblock)
                  ELSE
                     ! If each grid has the same ESTM file, ESTM file name does not include grid number, and only need to disaggregate once
                     FileOrigESTM = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_ESTM_Ts_data_' &
                                    //TRIM(ADJUSTL(ResInESTM_txt))//'.txt'
                     FileDscdESTM = TRIM(FileInputPath)//TRIM(FileCode)//'_'//TRIM(year_txt)//'_ESTM_Ts_data_' &
                                    //TRIM(ADJUSTL(tstep_txt))//'.txt'
                     IF (igrid == 1) THEN       !Disaggregate for the first grid only
                        CALL DisaggregateESTM(iblock)
                     ELSE                    !Then for subsequent grids simply copy data
                        ESTMForcingData(1:ReadLinesMetdata, 1:ncolsESTMdata, GridCounter) = ESTMForcingData(1:ReadLinesMetdata, &
                                                                                                            1:ncolsESTMdata, 1)
                     ENDIF
                  ENDIF

               ELSEIF (NperESTM == 1) THEN
                  ! Get ESTM forcing file name for this year for the first grid
                  FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
                  FileCodeXWG = TRIM(FileCode)//'_'//TRIM(year_txt) !File code without grid
                  ! For every block of ESTM data ------------------------------------
                  ! Initialise ESTM forcing data into 3-dimensional matrix
                  !write(*,*) 'Initialising ESTM data for block',iblock
                  IF (MultipleESTMFiles == 1) THEN   !If each grid has its own met file
                     FileESTMTs = TRIM(FileInputPath)//TRIM(FileCodeX)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                     !write(*,*) 'Calling GetESTMData...', FileCodeX, iblock, igrid
                     CALL SUEWS_GetESTMData(101)
                  ELSE                             !If one ESTM file used for all grids
                     FileESTMTs = TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                     !write(*,*) 'Calling GetESTMData...', FileCodeX, iblock, igrid
                     IF (igrid == 1) THEN       !Read for the first grid only
                        CALL SUEWS_GetESTMData(101)
                     ELSE                          !Then for subsequent grids simply copy data
                        ESTMForcingData(1:ReadLinesMetdata, 1:ncolsESTMdata, GridCounter) = ESTMForcingData(1:ReadLinesMetdata, &
                                                                                                            1:ncolsESTMdata, 1)
                     ENDIF
                  ENDIF
               ENDIF   !end of nperESTM statement
            ENDIF

            GridCounter = GridCounter + 1   !Increase GridCounter by 1 for next grid

         ENDDO !end loop over grids
         skippedLines = skippedLines + ReadLinesMetdata   !Increase skippedLines ready for next block
         skippedLinesOrig = skippedLinesOrig + ReadlinesOrigMetdata   !Increase skippedLinesOrig ready for next block
         skippedLinesOrigESTM = skippedLinesOrigESTM + ReadlinesOrigESTMdata   !Increase skippedLinesOrig ready for next block
         !write(*,*) iblock
         !write(*,*) ReadLinesMetdata, readlinesorigmetdata
         !write(*,*) skippedLines, skippedLinesOrig, skippedLinesOrig*Nper

         ! Initialise the modules on the first day
         ! Initialise CBL and SOLWEIG parts if required
         IF (iblock == 1) THEN
            IF ((CBLuse == 1) .OR. (CBLuse == 2)) CALL CBL_ReadInputData
         ENDIF

         ! NB: SOLWEIG is disabled in v2018a
         ! QUESTION: not sure if solweig should be within iblock if statement too?
         ! IF(SOLWEIGuse==1) CALL SOLWEIG_initial

         !write(*,*) 'Initialisation done'
         ! First stage: initialisation done ----------------------------------

         ! (2) Second stage: do calculations at 5-min time-steps -------------
         ! First set maximum value of ir
         IF (iblock == ReadBlocksMetData) THEN   !For last block of data in file
            irMax = nlinesMetdata - (iblock - 1)*ReadLinesMetdata
         ELSE
            irMax = ReadLinesMetdata
         ENDIF

         DO ir = 1, irMax   !Loop through rows of current block of met data
            GridCounter = 1    !Initialise counter for grids in each year

            DO igrid = 1, NumberOfGrids   !Loop through grids
               IF (Diagnose == 1) WRITE (*, *) 'Row (ir):', ir, '/', irMax, 'of block (iblock):', iblock, '/', ReadBlocksMetData, &
                  'Grid:', GridIDmatrix(igrid)

               ! Call model calculation code
               !  IF(ir==1) WRITE(*,*) 'Now running block ',iblock,'/',ReadBlocksMetData,' of year ',year_int,'...'
               WRITE (grid_txt, '(I10)') GridIDmatrix(igrid)   !Get grid ID as a text string
               FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
               IF (ir == 1 .AND. igrid == 1) THEN
                  WRITE (*, *) TRIM(ADJUSTL(FileCodeX)), &
                     ': Now running block ', iblock, '/', ReadBlocksMetData, ' of ', TRIM(year_txt), '...'
               ENDIF
               IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_Calculations...'
               ! print*, 'before cal:',sum(Tair24HR)
               CALL SUEWS_Calculations(GridCounter, ir, iblock, irMax)
               ! print*, 'after cal:',sum(Tair24HR)
               IF (Diagnose == 1) WRITE (*, *) 'SUEWS_Calculations finished...'

               ! Record iy and id for current time step to handle last row in yearly files (YYYY 1 0 0)
               !  IF(GridCounter == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
               IF (igrid == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
                  iy_prev_t = iy
                  id_prev_t = id
               ENDIF

               ! Write state information to new InitialConditions files
               IF (ir == irMax) THEN              !If last row...
                  IF (iblock == ReadBlocksMetData) THEN    !...of last block of met data
                     WRITE (grid_txt, '(I10)') GridIDmatrix(igrid)
                     !  WRITE(year_txtNext,'(I4)') year_int+1  !Get next year as a string format
                     !  FileCodeX     = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
                     !  FileCodeXNext = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txtNext)
                     !  CALL NextInitial(FileCodeXNext,year_int)
                     FileCodeXwy = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt)) !File code without year (HCW 24 May 2016)
                     IF (Diagnose == 1) WRITE (*, *) 'Calling NextInitial...'
                     CALL NextInitial(FileCodeXwy, year_int)

                  ENDIF
               ENDIF

               GridCounter = GridCounter + 1   !Increase GridCounter by 1 for next grid
            ENDDO !end loop over grids

            ! update simulation time since start
            dt_since_start = dt_since_start + tstep

            !! TODO: water movements between the grids needs to be taken into account here

         ENDDO !end loop over rows of met data

         ! Write output files in blocks --------------------------------
#ifdef nc
         IF (ncMode == 0) THEN
#endif
            DO igrid = 1, NumberOfGrids
               IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_Output...'
               CALL SUEWS_Output(irMax, iblock, igrid, year_int)
            ENDDO
#ifdef nc
         ENDIF

         IF (ncMode == 1) THEN
            ! write resulst in netCDF
            IF (Diagnose == 1) WRITE (*, *) 'Calling SUEWS_Output_nc...'
            CALL SUEWS_Output(irMax)
            ! write input information in netCDF as well for future development
            !  IF ( iblock==1 ) THEN
            !     CALL SiteSelect_txt2nc
            !  ENDIF
         ENDIF
#endif

      ENDDO !end loop over blocks of met data
      !-----------------------------------------------------------------------

      ! ---- Decallocate arrays ----------------------------------------------
      IF (Diagnose == 1) WRITE (*, *) 'Deallocating arrays in SUEWS_Program.f95...'
      DEALLOCATE (SurfaceChar)
      DEALLOCATE (MetForcingData)
      DEALLOCATE (MetForcingData_grid)
      DEALLOCATE (ModelOutputData)
      DEALLOCATE (dataOutSUEWS)
      DEALLOCATE (dataOutDailyState)
      ! IF (SnowUse == 1) THEN
      DEALLOCATE (dataOutSnow)
      ! DEALLOCATE(qn1_S_store)
      ! DEALLOCATE(qn1_S_av_store)
      ! DEALLOCATE(qn1_S_store_grid)
      ! DEALLOCATE(qn1_S_av_store_grid)
      DEALLOCATE (qn1_s_av_grids)
      DEALLOCATE (dqnsdt_grids)
      ! ENDIF
      ! IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
      DEALLOCATE (dataOutESTM) !ESTM output
      DEALLOCATE (ESTMForcingData)
      DEALLOCATE (Ts5mindata)
      DEALLOCATE (Ts5mindata_ir)
      DEALLOCATE (Tair24HR)
      ! ENDIF
      ! DEALLOCATE(TstepProfiles)
      ! DEALLOCATE(AHProf_tstep)
      ! DEALLOCATE(WUProfM_tstep)
      ! DEALLOCATE(WUProfA_tstep)
      ! DEALLOCATE(HumActivity_tstep)
      ! DEALLOCATE(TraffProf_tstep)
      ! DEALLOCATE(PopProf_tstep)
      ! DEALLOCATE(qn1_store)
      ! DEALLOCATE(qn1_av_store)
      ! DEALLOCATE(qn1_store_grid)
      ! DEALLOCATE(qn1_av_store_grid)
      ! DEALLOCATE(qhforCBL)
      ! DEALLOCATE(qeforCBL)
      DEALLOCATE (qn1_av_grids)
      DEALLOCATE (dqndt_grids)

      ! ----------------------------------------------------------------------

   ENDDO  !end loop over years

   ! ---- Decallocate array --------------------------------------------------
   ! Daily state needs to be outside year loop to transfer states between years
   IF (ALLOCATED(ModelDailyState)) DEALLOCATE (ModelDailyState)
   ! Also needs to happen at the end of the run
   IF (ALLOCATED(UseColumnsDataOut)) DEALLOCATE (UseColumnsDataOut)
   ! -------------------------------------------------------------------------

   ! get cpu time consumed
   CALL CPU_TIME(timeFinish)
   WRITE (*, *) "Time = ", timeFinish - timeStart, " seconds."

   !Write to problems.txt that run has completed
   IF (errorChoice == 0) THEN  !if file has not been opened previously
      OPEN (500, file='problems.txt')
      errorChoice = 1
   ELSE
      OPEN (500, file='problems.txt', position="append")
   ENDIF
   !Writing of the problem file
   WRITE (500, *) '--------------'
   WRITE (500, *) 'Run completed.'
   WRITE (500, *) '0'  ! Write out error code 0 if run completed
   CLOSE (500)

   ! Also print to screen
   WRITE (*, *) "----- SUEWS run completed -----"

   STOP 'finished'

   ! 313 CALL errorHint(11,TRIM(FileOrigMet),notUsed,notUsed,ios_out)
   ! 314 CALL errorHint(11,TRIM(FileMet),notUsed,notUsed,ios_out)
   ! 315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,NotUsedI)

END PROGRAM SUEWS_Program
