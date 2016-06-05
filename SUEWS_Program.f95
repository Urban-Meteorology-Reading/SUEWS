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
!
!Last modified by TS 14 Mar 2016  - Include AnOHM daily interation
!Last modified by HCW 25 Jun 2015 - Fixed bug in LAI calculation at year change
!Last modified by HCW 12 Mar 2015
!Last modified by HCW 26 Feb 2015
!Last modified by HCW 03 Dec 2014
!Last modified by LJ  30 Mar 2016 - Grid run order changed from linear to non-linear
!                 HCW 25 Jun 2015 - Fixed bug in LAI calculation at year change
!                 HCW 12 Mar 2015
!                 HCW 26 Feb 2015
!                 HCW 03 Dec 2014
!
! To do:
!   - Snow modules need updating for water balance part
!   - Water movement between grids (GridConnections) not yet coded
!   - Check all arrays are allocated/deallocated and initialised in the correct place
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

  CHARACTER(len = 4) ::year_txt, & !Year as a text string
       year_txtNext !Following year as a text string (used for NextInitial)
  CHARACTER(len = 20)::FileCodeX,& !Current file code
       FileCodeXNext !File code for the following year
  CHARACTER(len = 20)::grid_txt, & !Grid number as a text string (from FirstGrid to LastGrid)
       tstep_txt !Model timestep (in minutes) as a text string

  INTEGER:: nlinesLimit,&   !Max number of lines that can be read in one go for each grid
       NumberOfYears   !Number of years to be run

  INTEGER::igrid,& !Grid number (from FirstGrid to LastGrid)
       iv,       & !Block number (from 1 to ReadBlocksMetData)
       ir,irMax, & !Row number within each block (from 1 to irMax)
       year_int, & ! Year as an integer (from SiteSelect rather than met forcing file)
       iter,     & ! iteraion counter, AnOHM TS
       rr !Row of SiteSelect corresponding to current year and grid

  LOGICAL:: PrintPlace = .FALSE.   !Prints row, block, and grid number to screen if TRUE;

  REAL :: timeStart, timeFinish ! profiling use, AnOHM TS
  REAL :: xErr      ! error in Bo iteration, AnOHM TS 20160331
  LOGICAL, ALLOCATABLE :: flagRerunAnOHM(:)   ! iteration run to make Bo converge,AnOHM TS

  !  ! ---- Allocate arrays--------------------------------------------------
  !  allocate(SurfaceChar(NumberOfGrids,MaxNCols_c))   !Surface characteristics
  !  allocate(MetForcingData(1:ReadlinesMetdata,ncolumnsMetForcingData,NumberOfGrids))   !Met forcing data
  !  allocate(ModelOutputData(0:ReadlinesMetdata,MaxNCols_cMOD,NumberOfGrids))           !Data at model timestep
  !  allocate(dataOut(1:ReadlinesMetdata,ncolumnsDataOut,NumberOfGrids))                 !Main output array
  !  if (SOLWEIGuse == 1) then
  !     allocate(dataOutSOL(1:ReadlinesMetdata,28,NumberOfGrids))                        !SOLWEIG POI output
  !  endif
  !  if (CBLuse >= 1) then
  !     allocate(dataOutBL(1:ReadlinesMetdata,22,NumberOfGrids))                         !CBL output
  !  endif
  !  if (SnowUse == 1) then
  !     allocate(dataOutSnow(1:ReadlinesMetdata,ncolumnsDataOutSnow,NumberOfGrids))      !Snow output array
  !  endif
  !  if(QSChoice==4 .or. QSChoice==14) then
  !      allocate(dataOutESTM(1:ReadlinesMetdata,32,NumberOfGrids))
  !  endif
  !
  !  allocate(TstepProfiles(NumberOfGrids,6,24*NSH))  !Hourly profiles interpolated to model timestep
  !  allocate(AHProf_tstep(24*NSH,2))                 !Anthropogenic heat profiles at model timestep
  !  allocate(WUProfM_tstep(24*NSH,2))                !Manual water use profiles at model timestep
  !  allocate(WUProfA_tstep(24*NSH,2))                !Automatic water use profiles at model timestep
  !  !! Add snow clearing (?)
  !  ! ----------------------------------------------------------------------

  ! ---- Initialise arrays  !! Does this need to happen here??

  !==========================================================================

  ! start counting cpu time
  CALL cpu_TIME(timeStart)


  ! Initialise error file (0 -> problems.txt file is created)
  errorChoice=0

  ! Read RunControl.nml and all input files from SiteSelect spreadsheet.
  ! This is saved to SiteSelect datamatrix
  CALL overallRunControl

  ! First find first and last year of the current run
  FirstYear = MINVAL(INT(SiteSelect(:,c_Year)))
  LastYear  = MAXVAL(INT(SiteSelect(:,c_Year)))

  NumberOfYears = LastYear-FirstYear+1 !Find the number of years to run

  !Find the the number of grids within each year in SiteSelect and GridIDs
  !(N.B. need to have the same grids for each year)
  NumberOfGrids = INT(nlinesSiteSelect/NumberOfYears)

  !! Find the first and last grid numbers (N.B. need to have the same grids for each year)
  !FirstGrid = minval(int(SiteSelect(:,c_Grid)))
  !LastGrid  = maxval(int(SiteSelect(:,c_Grid)))
  IF(NumberOfGrids > MaxNumberOfGrids) THEN
     CALL ErrorHint(64,'No. of grids exceeds max. possible no. of grids.',REAL(MaxNumberOfGrids,KIND(1d0)),NotUsed,NumberOfGrids)
  ENDIF

  ALLOCATE (GridIDmatrix(NumberOfGrids)) !Get the nGrid numbers correctly
  DO igrid=1,NumberOfGrids
     GridIDmatrix(igrid) = INT(SiteSelect(igrid,c_Grid))
  ENDDO

  WRITE(*,*) '--------------------------------------------'
  WRITE(*,*) 'Years identified:',FirstYear,'to',LastYear
  WRITE(*,*) 'Grids identified:',NumberOfGrids,'grids'

  ! ---- Allocate arrays ----------------------------------------------------
  ! Daily state needs to be outside year loop to transfer states between years
  ALLOCATE(ModelDailyState(NumberOfGrids,MaxNCols_cMDS))   !DailyState
  ALLOCATE(DailyStateFirstOpen(NumberOfGrids))             !Initialization for header
  ALLOCATE(flagRerunAnOHM(NumberOfGrids))                  !flag for rerun AnOHM
  !   allocate(BoAnOHMStart(NumberOfGrids))       ! initial Bo
  !   allocate(BoAnOHMEnd(NumberOfGrids))         ! final Bo



  ! ---- Initialise arrays --------------------------------------------------
  ModelDailyState(:,:) = -999
  DailyStateFirstOpen(:) = 1

  flagRerunAnOHM(:) = .TRUE.
  ! -------------------------------------------------------------------------

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
     WRITE(grid_txt,'(I5)') GridIDmatrix(1)  !Get grid as a text string
     WRITE(tstep_txt,'(I5)') tstep/60  !Get tstep (in minutes) as a text string

     ! Get met file name for this year for this grid
     FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
     FileMet   = TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'

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
     !        nlinesLimit = int(floor(MaxLinesMet/real(NumberOfGrids,kind(1d0))))
     nlinesLimit = 24*nsh
     IF(nlinesMetData > nlinesLimit) THEN   !But restrict if this limit exceeds memory capacity
        ReadLinesMetData = nlinesLimit
     ENDIF
     !write(*,*) 'Met data will be read in chunks of',ReadlinesMetdata,'lines.'

     ! Find number of blocks of met data
     ReadBlocksMetData = INT(CEILING(REAL(nlinesMetData,KIND(1d0))/REAL(ReadLinesMetData,KIND(1d0))))
     !write(*,*) 'Met data will be read in',ReadBlocksMetData,'blocks.'

     ! ---- Allocate arrays--------------------------------------------------
     ALLOCATE(SurfaceChar(NumberOfGrids,MaxNCols_c))                                               !Surface characteristics
     ALLOCATE(MetForcingData(1:ReadlinesMetdata,ncolumnsMetForcingData,NumberOfGrids))             !Met forcing data
     ALLOCATE(ModelOutputData(0:ReadlinesMetdata,MaxNCols_cMOD,NumberOfGrids))                     !Data at model timestep
     ALLOCATE(dataOut(1:ReadlinesMetdata,ncolumnsDataOut,NumberOfGrids))                           !Main output array
     IF (SOLWEIGuse == 1) ALLOCATE(dataOutSOL(1:ReadlinesMetdata,28,NumberOfGrids))                !SOLWEIG POI output
     IF (CBLuse >= 1)  ALLOCATE(dataOutBL(1:ReadlinesMetdata,22,NumberOfGrids))                    !CBL output
     IF (SnowUse == 1) ALLOCATE(dataOutSnow(1:ReadlinesMetdata,ncolumnsDataOutSnow,NumberOfGrids)) !Snow output array
     IF (QSChoice==4 .OR. QSChoice==14) ALLOCATE(dataOutESTM(1:ReadlinesMetdata,32,NumberOfGrids)) !ESTM output array, TS 05 Jun 2016
     ALLOCATE(TstepProfiles(NumberOfGrids,6,24*NSH))                        !Hourly profiles interpolated to model timestep
     ALLOCATE(AHProf_tstep(24*NSH,2))                                       !Anthropogenic heat profiles at model timestep
     ALLOCATE(WUProfM_tstep(24*NSH,2))                                      !Manual water use profiles at model timestep
     ALLOCATE(WUProfA_tstep(24*NSH,2))                                      !Automatic water use profiles at model timestep
     !! Add snow clearing (?)

    !  ! test ESTM initialisation
    !  IF(QSChoice==4 .OR. QSChoice==14) THEN
    !     PRINT*, 'day:', iv
     !
    !     !  if ( iv>1 ) then
    !     CALL ESTM_initials(FileCodeX)
     !
    !     !  end if
    !  ENDIF


     ! ----------------------------------------------------------------------

     ! ---- Initialise arrays  !! Does this need to happen here??

     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     SkippedLines=0  !Initialise lines to be skipped in met forcing file

     DO iv=1,ReadBlocksMetData   !Loop through blocks of met data
        !write(*,*) iv,'/',ReadBlocksMetData,'ReadBlocksMetData (iv loop)'

        ! Model calculations are made in two stages:
        ! (1) initialise the run for each block of met data (iv from 1 to ReadBlocksMetData)
        ! (2) perform the actual model calculations (SUEWS_Calculations)


        ! (1) First stage: initialise run -----------------------------------
        GridCounter=1   !Initialise counter for grids in each year
        DO igrid=1,NumberOfGrids   !Loop through grids

           WRITE(grid_txt,'(I5)') GridIDmatrix(igrid)   !Get grid ID as a text string
           ! Get met forcing file name for this year for the first grid
           ! Can be something else than 1
           FileCodeX = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
           IF(iv==1) WRITE(*,*) 'Current FileCode: ', FileCodeX

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
              ! (b) get initial conditions
              CALL InitialState(FileCodeX,year_int,GridCounter,year_txt)
           ENDIF   !end first block of met data

           ! For every block of met data ------------------------------------
           ! Initialise met forcing data into 3-dimensional matrix
           !write(*,*) 'Initialising met data for block',iv

           IF(MultipleMetFiles == 1) THEN   !If each grid has its own met file
              FileMet=TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
              CALL SUEWS_InitializeMetData(1)
           ELSE                             !If one met file used for all grids
              FileMet=TRIM(FileInputPath)//TRIM(FileCodeX)//'_data_'//TRIM(ADJUSTL(tstep_txt))//'.txt'
              IF(igrid == 1) THEN       !Read for the first grid only
                 CALL SUEWS_InitializeMetData(1)
              ELSE                          !Then for subsequent grids simply copy data
                 MetForcingData(1:ReadlinesMetdata,1:24,GridCounter) = MetForcingData(1:ReadlinesMetdata,1:24,1)
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

        ! Initialise ESTM if required, TS 05 Jun 2016
        ! print*, "before call ESTM_initials:", FileCodeX
        IF(QSChoice==4 .OR. QSChoice==14) THEN
           PRINT*, 'day:', iv

          !  if ( iv>1 ) then
             CALL ESTM_initials(FileCodeX)

          !  end if
        ENDIF

        ! end if


        !write(*,*) 'Initialisation done'
        ! First stage: initialisation done ----------------------------------


        ! (2) Second stage: do calculations at 5-min timesteps --------------
        ! First set maximum value of ir
        IF(iv == ReadBlocksMetData) THEN   !For last block of data in file
           irMax = nlinesMetdata - (iv-1)*ReadLinesMetdata
        ELSE
           irMax = ReadLinesMetdata
        ENDIF

        ! iteration for AnOHM running by do-while, 12 Mar 2016 TS ------------
        iter       = 0
        BoAnOHMEnd = NAN

        flagRerunAnOHM = .TRUE.

        DO WHILE ( ANY(flagRerunAnOHM) .AND. iter < 20 )
           iter = iter+1
           !  PRINT*, 'iteration:',iter


           DO ir=1,irMax   !Loop through rows of current block of met data
              GridCounter=1    !Initialise counter for grids in each year


              DO igrid=1,NumberOfGrids   !Loop through grids
                 IF(PrintPlace) WRITE(*,*) 'Row (ir):', ir,'/',irMax,'of block (iv):', iv,'/',ReadBlocksMetData,&
                      'Grid:',GridIDmatrix(igrid)

                 !  ! Translate daily state back so as to keep water balance at beginning of a day
                 !  IF ( QSChoice==3 .AND. ir==1) THEN
                 !     CALL SUEWS_Translate(igrid,0,iv)
                 !  END IF

                 ! Call model calculation code
                 IF(ir==1) WRITE(*,*) 'Now running block ',iv,'/',ReadBlocksMetData,' of year ',year_int,'...'
                 CALL SUEWS_Calculations(GridCounter,ir,iv,irMax)


                 ! Record iy and id for current time step to handle last row in yearly files (YYYY 1 0 0)
                 !  IF(GridCounter == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
                 IF(igrid == NumberOfGrids) THEN   !Adjust only when the final grid has been run for this time step
                    iy_prev_t = iy
                    id_prev_t = id
                 ENDIF

                 ! Write state information to new InitialConditions files
                 IF(ir == irMax) THEN              !If last row...
                    IF(iv == ReadBlocksMetData) THEN    !...of last block of met data
                       WRITE(grid_txt,'(I5)') GridIDmatrix(igrid)
                       WRITE(year_txtNext,'(I4)') year_int+1  !Get next year as a string format
                       FileCodeX     = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txt)
                       FileCodeXNext = TRIM(FileCode)//TRIM(ADJUSTL(grid_txt))//'_'//TRIM(year_txtNext)
                       CALL NextInitial(FileCodeXNext,year_int)
                    ENDIF
                 ENDIF


                 ! print AnOHM coeffs. info.:
                 !  IF ( it == 0 .AND. imin == 5 )  THEN
                 !     WRITE(*, '(a13,2f10.4)') 'Start: a1, a2',a1AnOHM_grids(id,igrid),a2AnOHM_grids(id,igrid)
                 !  ENDIF

                 IF ( ir == irMax-10) THEN
                    xErr = ABS(a1AnOHM(igrid)-a1AnOHM_grids(id,igrid))/ABS(a1AnOHM(igrid))+&
                         ABS(a2AnOHM(igrid)-a2AnOHM_grids(id,igrid))/ABS(a2AnOHM(igrid))
                    xErr = ABS(xErr/2)
                    !     WRITE(*, '(a13,2f10.4)') 'End: a1, a2',a1AnOHM(igrid),a2AnOHM(igrid)
                    !     WRITE(*, '(a10,f10.4,2x,i2,x,i2)') 'Error(%)',xErr*100,it,imin
                 ENDIF



                 IF ( QSChoice == 3 .AND. ir == irMax .AND. xErr < 0.1) THEN
                    flagRerunAnOHM(igrid)=.FALSE.
                    ! WRITE(unit=*, fmt=*) '*********'
                    ! WRITE(unit=*, fmt=*) 'converged:'
                    ! WRITE(*, '(a8,f10.6)') 'a1 Start',a1AnOHM_grids(id,igrid)
                    ! WRITE(*, '(a8,f10.6)') 'a1 End',a1AnOHM(igrid)
                    ! WRITE(*, '(a8,f10.6)') 'diff.',ABS(a1AnOHM(igrid)-a1AnOHM_grids(id,igrid))
                    ! WRITE(unit=*, fmt=*) '*********'
                 ENDIF

                 ! bypass the do-while loop for converge checking
                 IF ( QSChoice /= 3 ) flagRerunAnOHM(igrid)=.FALSE.


                 GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid
              ENDDO !end loop over grids

              !!water movements between the grids needs to be taken into account here ??

           ENDDO !end loop over rows of met data

        ENDDO ! do-while loop for AnOHM end --------------------------


        ! Write output files in blocks --------------------------------
        DO igrid=1,NumberOfGrids
           !              call SUEWS_Output(igrid,year_int,iv,irMax,GridIDmatrix(igrid))
           CALL SUEWS_Output(igrid,year_int,iv,irMax)
        ENDDO

     ENDDO !end loop over blocks of met data
     !-----------------------------------------------------------------------

     ! ---- Decallocate arrays ----------------------------------------------
     DEALLOCATE(SurfaceChar)
     DEALLOCATE(MetForcingData)
     DEALLOCATE(ModelOutputData)
     DEALLOCATE(dataOut)
     IF (SnowUse == 1) DEALLOCATE(dataOutSnow)
     DEALLOCATE(TstepProfiles)
     DEALLOCATE(AHProf_tstep)
     DEALLOCATE(WUProfM_tstep)
     DEALLOCATE(WUProfA_tstep)
     ! ----------------------------------------------------------------------

  ENDDO  !end loop over years

  ! ---- Decallocate array --------------------------------------------------
  ! Daily state needs to be outside year loop to transfer states between years
  DEALLOCATE(ModelDailyState)
  ! -------------------------------------------------------------------------

  ! get cpu time consumed
  CALL cpu_TIME(timeFinish)
  WRITE(*,*) "Time = ",timeFinish-timeStart," seconds."

  STOP 'finished'


314 CALL errorHint(11,TRIM(FileMet),notUsed,notUsed,ios_out)

END PROGRAM SUEWS_Program
