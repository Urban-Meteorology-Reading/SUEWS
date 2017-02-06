!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!Main program of SUEWS version 1.0  -------------------------------------------
! Model calculations are made in two stages:
! (1) initialise the run for each block of met data (iv from 1 to ReadBlocksMetData)
! (2) perform the actual model calculations (SUEWS_Calculations)
! After reading in all input information from SiteInfo, the code loops
!  - over years
!  - then over blocks (met data read in and stored for each block)
!  - then over rows
!  - then over grids
!
!Last modified by HCW 12 Jan 2017 - Changes to InitialConditions
!Last modified by HCW 26 Aug 2016 - CO2 flux added
!Last modified by HCW 04 Jul 2016 - GridID can now be up to 10 digits long
!Last modified by HCW 29 Jun 2016 - Reversed over-ruling of ReadLinesMetData so this is not restricted here to one day
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

  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE initial
  USE sues_data
  USE time

  IMPLICIT NONE

  CHARACTER(len = 4) :: year_txt  !Year as a text string

  CHARACTER(len = 20):: FileCodeX,& !Current file code
       FileCodeXWY,& !File code without year
       FileCodeXWG !File code without grid
  CHARACTER(len = 20):: grid_txt, & !Grid number as a text string (from FirstGrid to LastGrid)
       tstep_txt   !Model timestep (in minutes) as a text string

  INTEGER:: nlinesLimit,&   !Max number of lines that can be read in one go for each grid
       NumberOfYears   !Number of years to be run

  INTEGER:: igrid,& !Grid number (from FirstGrid to LastGrid)
       iv,       & !Block number (from 1 to ReadBlocksMetData)
       ir,irMax, & !Row number within each block (from 1 to irMax)
       year_int, & ! Year as an integer (from SiteSelect rather than met forcing file)
       rr !Row of SiteSelect corresponding to current year and grid

  LOGICAL:: PrintPlace = .FALSE.   !Prints row, block, and grid number to screen if TRUE;

  REAL::  timeStart, timeFinish ! profiling use, AnOHM TS
  INTEGER(KIND(1d0)),ALLOCATABLE :: GridIDmatrix0(:)
  ! integer :: ncMode = 1 ! if the output should be written in netCDF, TS, 08 Dec 201616
  ! REAL :: xErr      ! error in Bo iteration, AnOHM TS 20160331
  ! LOGICAL, ALLOCATABLE :: flagRerunAnOHM(:)   ! iteration run to make Bo converge,AnOHM TS

  !==========================================================================

  ! Start counting cpu time
  CALL cpu_TIME(timeStart)

  WRITE(*,*) '========================================================'
  WRITE(*,*) 'Running ',progname

  ! Initialise error file (0 -> problems.txt file is created)
  errorChoice=0
  ! Initialise OutputFormats to 1 so that output format is written out only once per run
  OutputFormats = 1
  ! Set Diagnose switch to off (0). If Diagnose = 1 is set in RunControl, model progress will be printed
  Diagnose = 0  

  ! Read RunControl.nml and all .txt input files from SiteSelect spreadsheet

  CALL overallRunControl

  ! Find first and last year of the current run
  FirstYear = MINVAL(INT(SiteSelect(:,c_Year)))
  LastYear  = MAXVAL(INT(SiteSelect(:,c_Year)))

  NumberOfYears = LastYear-FirstYear+1 !Find the number of years to run

  !Find the the number of grids within each year in SUEWS_SiteSelect.txt
  ! N.B. need to have the same grids for each year
  NumberOfGrids = INT(nlinesSiteSelect/NumberOfYears)

  !! Find the first and last grid numbers (N.B. need to have the same grids for each year)
  !FirstGrid = minval(int(SiteSelect(:,c_Grid)))
  !LastGrid  = maxval(int(SiteSelect(:,c_Grid)))
  IF(NumberOfGrids > MaxNumberOfGrids) THEN
     CALL ErrorHint(64,'No. of grids exceeds max. possible no. of grids.',REAL(MaxNumberOfGrids,KIND(1d0)),NotUsed,NumberOfGrids)
  ENDIF

  ALLOCATE (GridIDmatrix(NumberOfGrids)) !Get the nGrid numbers correctly
  ALLOCATE (GridIDmatrix0(NumberOfGrids)) !Get the nGrid numbers correctly

  DO igrid=1,NumberOfGrids
     GridIDmatrix(igrid) = INT(SiteSelect(igrid,c_Grid))
     GridIDmatrix0(igrid) =INT(SiteSelect(igrid,c_Grid))
  ENDDO

  ! sort grid matrix to conform the geospatial layout as in QGIS
  IF (ncMode==1) CALL sortGrid(GridIDmatrix0,GridIDmatrix,nRow,nCol)
  ! GridIDmatrix0 stores the grid ID in the original order


  ! GridIDmatrix=GridIDmatrix0
  WRITE(*,*) '--------------------------------------------'
  WRITE(*,*) 'Years identified:',FirstYear,'to',LastYear
  WRITE(*,*) 'No. grids identified:',NumberOfGrids,'grids'

  ! ---- Allocate arrays ----------------------------------------------------
  ! Daily state needs to be outside year loop to transfer states between years
  ALLOCATE(ModelDailyState(NumberOfGrids,MaxNCols_cMDS))   !DailyState
  ALLOCATE(DailyStateFirstOpen(NumberOfGrids))             !Initialisation for header
  !ALLOCATE(flagRerunAnOHM(NumberOfGrids))                  !flag for rerun AnOHM
  !allocate(BoAnOHMStart(NumberOfGrids))                    !initial Bo
  !allocate(BoAnOHMEnd(NumberOfGrids))                      !final Bo

  ! print*, 'good 1'
  ! ---- Initialise arrays --------------------------------------------------
  ModelDailyState(:,:) = -999
  DailyStateFirstOpen(:) = 1
  !flagRerunAnOHM(:) = .TRUE.

  ! -------------------------------------------------------------------------

  ! Initialise ESTM (reads ESTM nml, should only run once)
  IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
     IF(Diagnose==1) write(*,*) 'Calling ESTM_initials...' 
     CALL ESTM_initials
  ENDIF

  !==========================================================================
  DO year_int=FirstYear,LastYear   !Loop through years

     WRITE(*,*) ' '
     WRITE(year_txt,'(I4)') year_int  !Get year as a text string

     ! Find number of days in the current year
     CALL LeapYearCalc (year_int,nofDaysThisYear)

     !-----------------------------------------------------------------------
     ! Find number of lines in met forcing file for current year (nlinesMetdata)
     ! Need to know how many lines will be read each iteration
     ! Use first grid as an example as the number of lines is the same for all grids
     ! within one year
     WRITE(grid_txt,'(I10)') GridIDmatrix(1)  !Get grid as a text string
     WRITE(tstep_txt,'(I5)') tstep/60  !Get tstep (in minutes) as a text string

     ! Get met file name for this year for this grid
     FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
     FileMet   = TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
     !If each grid has the same met file, met file name does not include grid number (HCW 27 Jun 2016)
     IF(MultipleMetFiles /= 1) THEN
        FileCodeXWG=TRIM(FileCode)//'_'//TRIM(year_txt) !File code without grid
        FileMet=TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
     ENDIF

     ! Open this example met file
     OPEN(10,file=TRIM(FileMet),status='old',err=314)
     CALL skipHeader(10,SkipHeaderMet)  !Skip header

     ! Find number of lines in met file
     nlinesMetdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
     DO
        READ(10,*) iv
        IF (iv == -9) EXIT
        nlinesMetdata = nlinesMetdata + 1
     ENDDO
     CLOSE(10)
     !-----------------------------------------------------------------------

     ! To conserve memory, read met data in blocks
     ! Find number of lines that can be read in each block (i.e. read in at once)
     ReadLinesMetData = nlinesMetData   !Initially set limit as the size of the met file (N.B.solves problem with Intel fortran)
     nlinesLimit = INT(FLOOR(MaxLinesMet/REAL(NumberOfGrids,KIND(1d0))))  !Uncommented HCW 29 Jun 2016
     !nlinesLimit = 24*nsh  !Commented out HCW 29 Jun 2016
     IF(nlinesMetData > nlinesLimit) THEN   !But restrict if this limit exceeds memory capacity
        ReadLinesMetData = nlinesLimit
     ENDIF
     ! make sure the metblocks read in consists of complete diurnal cycles, TS 08 Jul 2016
     ReadLinesMetData = INT(MAX(nsd*(ReadLinesMetData/nsd), nsd))

     WRITE(*,*) 'Met data will be read in chunks of ',ReadlinesMetdata,'lines.'

     ! Find number of blocks of met data
     ReadBlocksMetData = INT(CEILING(REAL(nlinesMetData,KIND(1d0))/REAL(ReadLinesMetData,KIND(1d0))))
     WRITE(*,*) 'Processing current year in ',ReadBlocksMetData,'blocks.'

     ! ---- Allocate arrays--------------------------------------------------
     IF(Diagnose==1) write(*,*) 'Allocating arrays in SUEWS_Program.f95...' 
     ALLOCATE(SurfaceChar(NumberOfGrids,MaxNCols_c))                                   !Surface characteristics
     ALLOCATE(MetForcingData(ReadlinesMetdata,ncolumnsMetForcingData,NumberOfGrids))   !Met forcing data
     ALLOCATE(ModelOutputData(0:ReadlinesMetdata,MaxNCols_cMOD,NumberOfGrids))         !Data at model timestep
     ALLOCATE(dataOut(ReadlinesMetdata,ncolumnsDataOut,NumberOfGrids))                 !Main output array
     IF (SOLWEIGuse == 1) ALLOCATE(dataOutSOL(ReadlinesMetdata,ncolumnsdataOutSOL,NumberOfGrids))     !SOLWEIG POI output
     IF (CBLuse >= 1)     ALLOCATE(dataOutBL(ReadlinesMetdata,ncolumnsdataOutBL,NumberOfGrids))       !CBL output
     IF (SnowUse == 1)    ALLOCATE(dataOutSnow(ReadlinesMetdata,ncolumnsDataOutSnow,NumberOfGrids))   !Snow output
     IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14) ALLOCATE(dataOutESTM(ReadlinesMetdata,32,NumberOfGrids)) !ESTM output
     ALLOCATE(TstepProfiles(NumberOfGrids,10,24*NSH))   !Hourly profiles interpolated to model timestep
     ALLOCATE(AHProf_tstep(24*NSH,2))                   !Anthropogenic heat profiles at model timestep
     ALLOCATE(WUProfM_tstep(24*NSH,2))                  !Manual water use profiles at model timestep
     ALLOCATE(WUProfA_tstep(24*NSH,2))                  !Automatic water use profiles at model timestep
     ALLOCATE(CO2m_tstep(24*NSH,2))
     ALLOCATE(qn1_store(NSH,NumberOfGrids))
     ALLOCATE(qn1_av_store(2*NSH+1,NumberOfGrids))
     !! Add snow clearing (?)

     qn1_store(:,:) = NAN ! Initialise to -999
     qn1_av_store(:,:) = NAN ! Initialise to -999
     ! Initialise other arrays here???


     ! Check number of lines in ESTM forcing file (nlinesESTMdata)
     IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
        IF(MultipleESTMFiles  == 1) THEN  !if separate ESTM files for each grid
           FileESTMTs=TRIM(FileInputPath)//TRIM(FileCodeX)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
           !write(*,*) 'Calling ESTM initials...', FileCodeX, iv, igrid
           !CALL ESTM_initials(FileCodeX)
        ELSE   !if one ESMT file for all grids
           FileESTMTs=TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
           !write(*,*) 'Calling ESTM initials...', FileCodeX, iv, igrid
           !CALL ESTM_initials(FileCodeX)
        ENDIF
        ! Open this example met file
        OPEN(101,file=TRIM(FileESTMTs),status='old',action='read',err=315)
        CALL skipHeader(101,SkipHeaderMet)  !Skip header
        ! Find number of lines in ESTM data file
        nlinesESTMdata = 0
        DO
           READ(101,*) iv
           IF (iv == -9) EXIT
           nlinesESTMdata = nlinesESTMdata + 1
        ENDDO
        CLOSE(101)

        ! Check ESTM data and met data are same length (so that ESTM file can be read same blocks as met data
        IF(nlinesESTMdata /= nlinesMetdata) THEN
           CALL ErrorHint(66,'ESTM input file different length to met forcing file',REAL(nlinesESTMdata,KIND(1d0)), &
                NotUsed,nlinesMetdata)
        ENDIF
        ! Allocate array to receive ESTM forcing data
        ALLOCATE(ESTMForcingData(1:ReadlinesMetdata,ncolsESTMdata,NumberOfGrids))
        ALLOCATE(Ts5mindata(1:ReadlinesMetdata,ncolsESTMdata))
        ALLOCATE(Tair24HR(24*nsh))
     ENDIF
     ! ----------------------------------------------------------------------

     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     SkippedLines=0  !Initialise lines to be skipped in met forcing file

     DO iv=1,ReadBlocksMetData   !Loop through blocks of met data
        WRITE(*,*) iv,'/',ReadBlocksMetData,'ReadBlocksMetData (iv loop)'

        ! Model calculations are made in two stages:
        ! (1) initialise the run for each block of met data (iv from 1 to ReadBlocksMetData)
        ! (2) perform the actual model calculations (SUEWS_Calculations)

        ! (1) First stage: initialise run -----------------------------------
        GridCounter=1   !Initialise counter for grids in each year
        DO igrid=1,NumberOfGrids   !Loop through grids

           WRITE(grid_txt,'(I10)') GridIDmatrix(igrid)   !Get grid ID as a text string
           ! Get met forcing file name for this year for the first grid
           ! Can be something else than 1
           FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
           !  IF(iv==1) WRITE(*,*) 'Current FileCode: ', FileCodeX

           ! For the first block of met data --------------------------------
           IF(iv == 1) THEN
              !write(*,*) 'First block of data - doing initialisation'
              ! (a) Transfer characteristics from SiteSelect to correct row of SurfaceChar
              DO rr=1,nlinesSiteSelect
                 !Find correct grid and year
                 IF(SiteSelect(rr,c_Grid)==GridIDmatrix(igrid).AND.SiteSelect(rr,c_Year)==year_int) THEN
                    !write(*,*) 'Match found (grid and year) for rr = ', rr
                    CALL InitializeSurfaceCharacteristics(GridCounter,rr)
                    EXIT
                 ELSEIF(rr == nlinesSiteSelect) THEN
                    WRITE(*,*) 'Program stopped! Year',year_int,'and/or grid',igrid,'not found in SiteSelect.txt.'
                    CALL ErrorHint(59,'Cannot find year and/or grid in SiteSelect.txt',REAL(igrid,KIND(1d0)),NotUsed,year_int)
                 ENDIF
              ENDDO
              !! (b) get initial conditions (moved to below, HCW 12 Jan 2017)
              !write(*,*) 'Calling InitialState'
              !CALL InitialState(FileCodeX,year_int,GridCounter,NumberOfGrids)
           ENDIF   !end first block of met data

           ! For every block of met data ------------------------------------
           ! Initialise met forcing data into 3-dimensional matrix
           !write(*,*) 'Initialising met data for block',iv
           IF(MultipleMetFiles == 1) THEN   !If each grid has its own met file
              FileMet=TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
              CALL SUEWS_InitializeMetData(1)
           ELSE                             !If one met file used for all grids
              !FileMet=TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
              ! If one met file used for all grids, look for met file with no grid code (FileCodeXWG)
              FileMet=TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
              IF(igrid == 1) THEN       !Read for the first grid only
                 CALL SUEWS_InitializeMetData(1)
              ELSE                          !Then for subsequent grids simply copy data
                 MetForcingData(1:ReadlinesMetdata,1:24,GridCounter) = MetForcingData(1:ReadlinesMetdata,1:24,1)
              ENDIF
           ENDIF

           ! Only for the first block of met data, read initial conditions (moved from above, HCW 12 Jan 2017)
           IF(iv == 1) THEN
              !write(*,*) ' Now calling InitialState'
              CALL InitialState(FileCodeX,year_int,GridCounter,NumberOfGrids)
           ENDIF

           ! Initialise ESTM if required, TS 05 Jun 2016; moved inside grid loop HCW 27 Jun 2016
           IF(StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
              IF(MultipleESTMFiles  == 1) THEN  !if separate ESTM files for each grid
                 FileESTMTs=TRIM(FileInputPath)//TRIM(FileCodeX)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                 !write(*,*) 'Calling GetESTMData...', FileCodeX, iv, igrid
                 CALL SUEWS_GetESTMData(101)
              ELSE   !if one ESMT file for all grids
                 FileESTMTs=TRIM(FileInputPath)//TRIM(FileCodeXWG)//'_ESTM_Ts_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
                 !write(*,*) 'Calling GetESTMData...', FileCodeX, iv, igrid
                 IF(igrid == 1) THEN       !Read for the first grid only
                    CALL SUEWS_GetESTMData(101)
                 ELSE                          !Then for subsequent grids simply copy data
                    ESTMForcingData(1:ReadlinesMetdata,1:ncolsESTMdata,GridCounter) = ESTMForcingData(1:ReadlinesMetdata, &
                         1:ncolsESTMdata,1)
                 ENDIF
              ENDIF
           ENDIF

           GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid

        ENDDO !end loop over grids
        skippedLines = skippedLines + ReadlinesMetdata   !Increase skippedLines ready for next block

        ! Initialise the modules on the first day
        ! if ( iv==1 ) then
        ! Initialise CBL and SOLWEIG parts if required
        IF((CBLuse==1).OR.(CBLuse==2)) CALL CBL_ReadInputData
        IF(SOLWEIGuse==1) CALL SOLWEIG_initial

        !write(*,*) 'Initialisation done'
        ! First stage: initialisation done ----------------------------------

        ! (2) Second stage: do calculations at 5-min timesteps --------------
        ! First set maximum value of ir
        IF(iv == ReadBlocksMetData) THEN   !For last block of data in file
           irMax = nlinesMetdata - (iv-1)*ReadLinesMetdata
        ELSE
           irMax = ReadLinesMetdata
        ENDIF

        ! ! iteration for AnOHM running by do-while, 12 Mar 2016 TS ------------
        ! iter       = 0
        ! BoAnOHMEnd = NAN

        ! flagRerunAnOHM = .TRUE.

        ! DO WHILE ( ANY(flagRerunAnOHM) .AND. iter < 20 )
        !    iter = iter+1
        !  PRINT*, 'iteration:',iter

        DO ir=1,irMax   !Loop through rows of current block of met data
           GridCounter=1    !Initialise counter for grids in each year

           DO igrid=1,NumberOfGrids   !Loop through grids
              IF(PrintPlace) WRITE(*,*) 'Row (ir):', ir,'/',irMax,'of block (iv):', iv,'/',ReadBlocksMetData,&
                   'Grid:',GridIDmatrix(igrid)
              IF(Diagnose==1) WRITE(*,*) 'Row (ir):', ir,'/',irMax,'of block (iv):', iv,'/',ReadBlocksMetData,&
                   'Grid:',GridIDmatrix(igrid)     

              !  ! Translate daily state back so as to keep water balance at beginning of a day
              !  IF ( StorageHeatMethod==3 .AND. ir==1) THEN
              !     CALL SUEWS_Translate(igrid,0,iv)
              !  END IF

              ! Call model calculation code
              !  IF(ir==1) WRITE(*,*) 'Now running block ',iv,'/',ReadBlocksMetData,' of year ',year_int,'...'
              WRITE(grid_txt,'(I10)') GridIDmatrix(igrid)   !Get grid ID as a text string
              FileCodeX=TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
              IF(ir==1) THEN
                 WRITE(*,*) TRIM(ADJUSTL(FileCodeX)),': Now running block ',iv,'/',ReadBlocksMetData,' of ',TRIM(year_txt),'...'
              ENDIF
              IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_Calculations...'
              CALL SUEWS_Calculations(GridCounter,ir,iv,irMax)
              IF(Diagnose==1) WRITE(*,*) 'SUEWS_Calculations finished...'
              
              ! Record iy and id for current time step to handle last row in yearly files (YYYY 1 0 0)
              !  IF(GridCounter == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
              IF(igrid == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
                 iy_prev_t = iy
                 id_prev_t = id
              ENDIF

              ! Write state information to new InitialConditions files
              IF(ir == irMax) THEN              !If last row...
                 IF(iv == ReadBlocksMetData) THEN    !...of last block of met data
                    WRITE(grid_txt,'(I10)') GridIDmatrix(igrid)
                    !  WRITE(year_txtNext,'(I4)') year_int+1  !Get next year as a string format
                    !  FileCodeX     = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
                    !  FileCodeXNext = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txtNext)
                    !  CALL NextInitial(FileCodeXNext,year_int)
                    FileCodeXwy=TRIM(FileCode)//TRIM(ADJUSTL(grid_txt)) !File code without year (HCW 24 May 2016)
                    IF(Diagnose==1) WRITE(*,*) 'Calling NextInitial...'
                    CALL NextInitial(FileCodeXwy,year_int)

                 ENDIF
              ENDIF


              ! print AnOHM coeffs. info.:
              !  IF ( it == 0 .AND. imin == 5 )  THEN
              !     WRITE(*, '(a13,2f10.4)') 'Start: a1, a2',a1AnOHM_grids(id,igrid),a2AnOHM_grids(id,igrid)
              !  ENDIF

              ! IF ( ir == irMax-10) THEN
              !    xErr = ABS(a1AnOHM(igrid)-a1AnOHM_grids(id,igrid))/ABS(a1AnOHM(igrid))+&
              !         ABS(a2AnOHM(igrid)-a2AnOHM_grids(id,igrid))/ABS(a2AnOHM(igrid))
              !    xErr = ABS(xErr/2)
              !    !     WRITE(*, '(a13,2f10.4)') 'End: a1, a2',a1AnOHM(igrid),a2AnOHM(igrid)
              !    !     WRITE(*, '(a10,f10.4,2x,i2,x,i2)') 'Error(%)',xErr*100,it,imin
              ! ENDIF



              ! IF ( StorageHeatMethod == 3 .AND. ir == irMax .AND. xErr < 0.1) THEN
              !    flagRerunAnOHM(igrid)=.FALSE.
              !    ! WRITE(unit=*, fmt=*) '*********'
              !    ! WRITE(unit=*, fmt=*) 'converged:'
              !    ! WRITE(*, '(a8,f10.6)') 'a1 Start',a1AnOHM_grids(id,igrid)
              !    ! WRITE(*, '(a8,f10.6)') 'a1 End',a1AnOHM(igrid)
              !    ! WRITE(*, '(a8,f10.6)') 'diff.',ABS(a1AnOHM(igrid)-a1AnOHM_grids(id,igrid))
              !    ! WRITE(unit=*, fmt=*) '*********'
              ! ENDIF

              ! bypass the do-while loop for converge checking
              ! IF ( StorageHeatMethod /= 3 ) flagRerunAnOHM(igrid)=.FALSE.


              GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid
           ENDDO !end loop over grids

           !!water movements between the grids needs to be taken into account here ??

        ENDDO !end loop over rows of met data
        ! print*, 'finish running:', iv

        ! ENDDO ! do-while loop for AnOHM end --------------------------


        ! Write output files in blocks --------------------------------
        IF ( ncMode .EQ. 0 ) THEN
           DO igrid=1,NumberOfGrids
IF(Diagnose==1) WRITE(*,*) 'Calling SUEWS_Output...'
              CALL SUEWS_Output(igrid,year_int,iv,irMax,GridIDmatrix(igrid))  !GridIDmatrix required for correct naming of output files
              !CALL SUEWS_Output(igrid,year_int,iv,irMax)
           ENDDO
        ELSE
           ! write resulst in netCDF
           CALL SUEWS_Output_nc(year_int,iv,irMax)
           ! write input information in netCDF as well for future development
           IF ( iv==1 ) THEN
              CALL SiteSelect_txt2nc
           ENDIF
        ENDIF
        ! print*, 'finish output:',iv

     ENDDO !end loop over blocks of met data
     !-----------------------------------------------------------------------

     ! ---- Decallocate arrays ----------------------------------------------
     IF(Diagnose==1) WRITE(*,*) 'Deallocating arrays in SUEWS_Program.f95...'
     DEALLOCATE(SurfaceChar)
     DEALLOCATE(MetForcingData)
     DEALLOCATE(ModelOutputData)
     DEALLOCATE(dataOut)
     IF (SnowUse == 1) DEALLOCATE(dataOutSnow)
     DEALLOCATE(TstepProfiles)
     DEALLOCATE(AHProf_tstep)
     DEALLOCATE(WUProfM_tstep)
     DEALLOCATE(WUProfA_tstep)
     DEALLOCATE(CO2m_tstep)
     DEALLOCATE(qn1_store)
     DEALLOCATE(qn1_av_store)
     ! ----------------------------------------------------------------------

  ENDDO  !end loop over years

  ! ---- Decallocate array --------------------------------------------------
  ! Daily state needs to be outside year loop to transfer states between years
  DEALLOCATE(ModelDailyState)
  ! Also needs to happen at the end of the run
  DEALLOCATE(UseColumnsDataOut)
  ! -------------------------------------------------------------------------

  ! get cpu time consumed
  CALL cpu_TIME(timeFinish)
  WRITE(*,*) "Time = ",timeFinish-timeStart," seconds."

  !Write to problems.txt that run has completed
  IF (errorChoice==0) THEN  !if file has not been opened previously
     OPEN(500,file='problems.txt')
     errorChoice=1
  ELSE
     OPEN(500,file='problems.txt',position="append")
  ENDIF
  !Writing of the problem file
  WRITE(500,*) '--------------'
  WRITE(500,*) 'Run completed.'
  WRITE(500,*) '0'  ! Write out error code 0 if run completed

  ! Also print to screen
  WRITE(*,*) "----- SUEWS run completed -----"

  STOP 'finished'


314 CALL errorHint(11,TRIM(FileMet),notUsed,notUsed,ios_out)
315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,NotUsedI)

END PROGRAM SUEWS_Program
