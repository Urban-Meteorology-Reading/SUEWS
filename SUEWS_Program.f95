!Main program of SUEWS version 1.0
!Last modified by HCW 03 Dec 2014
! To do:
!	- CBL and SOLWEIG (calls to _initial subroutines) not tested
! 	- Snow modules not tested since v2014b -> v2014c
! 	- water movement between grids (GridConnections) not yet coded
!	- check all arrays are intialised and initialised in the correct place
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 program SUEWS_Program

    use allocateArray
    use ColNamesInputFiles
    use data_in
    use defaultNotUsed   
    use initial
    use sues_data
    use time
    
    IMPLICIT NONE

    integer:: year_int			!Year as an integer
    character(len=4)::  year_txt,&	!Year as a text string
    			year_txtNext	!Following year as a text string (used for NextInitial)
    character(len=20):: grid_txt	!Grid as a text string (FirstGrid to LastGrid)
    character(len=20):: FileCodeX,&	!Has format SsGGGGG_YYYY
    			FileCodeXNext	  
    integer:: ReadMetTimes		!Number of blocks of met data
    integer:: errFileYes
    integer:: i,&	   ! FirstGrid to LastGrid (i.e. grid number)
    		  iv,&	   ! 1 to ReadMetTimes (i.e. met data block number)
    		  ir,irMax,& ! Row of met data within each chunk of met data
    		  rr 	   ! Row of SiteSelect corresponding to current year and grid
    		   
    !==========================================================================
    
    ! Initialise error file (0 -> problems.txt file is created)
    errorChoice=0

    ! Read RunControl.nml and input files from SiteInfo spreadsheet
    call overallRunControl

    ! First find first and last year of the current run
    FirstYear = minval(int(SiteSelect(:,c_Year)))
    LastYear  = maxval(int(SiteSelect(:,c_Year)))
    
    ! Find the first and last grid numbers (N.B. need to have the same grids for each year)
    FirstGrid = minval(int(SiteSelect(:,c_Grid))) 
    LastGrid  = maxval(int(SiteSelect(:,c_Grid))) 
    NumberOfGrids = LastGrid-FirstGrid+1   !Number of grids of the current run
    
    write(*,*) '--------------------------------------------'
    write(*,*) 'Years identified:',FirstYear,'to',LastYear
    write(*,*) 'Grids identified:',FirstGrid,'to',LastGrid

    ! ---- Allocate array -----------------------------------------------------
    ! Daily state needs to be outside year loop to transfer states between years
    allocate(ModelDailyState(NumberOfGrids,MaxNCols_cMDS))   !DailyState      
    allocate(DailyStateFirstOpen(NumberOfGrids))             !Initialization for header
    ! -------------------------------------------------------------------------
 
    ! ---- Initialise arrays
    ModelDailyState(:,:) = -999
    DailyStateFirstOpen(:) = 1

    !==========================================================================
    DO year_int=FirstYear,LastYear   !Loop through years
    
       write(*,*) ' '
       write(*,*) 'Now running year',year_int
      
       write(year_txt,'(I4)') year_int  !Get year as a text string

       ! Find number of days in the current year
       call LeapYearCalc (year_int,nofDaysThisYear)
       
       !-----------------------------------------------------------------------
       ! Find number of lines in met forcing file for current year (nlinesMetdata)
       !  Need to know how many lines will be read each iteration
       !  Use FirstGrid as an example
       write(grid_txt,'(I5)') FirstGrid  !Get grid as a text string
       ! Get met file name for this year for this grid
       FileCodeX=trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(year_txt)
       FileMet=trim(FileInputPath)//trim(FileCodeX)//'_data.txt'

       ! Open this example met file
       open(10,file=trim(FileMet),status='old',err=314)
       call skipHeader(10,SkipHeaderMet)  !Skip header
       ! Find number of lines in met file
       nlinesMetdata = 0   !Initialise nlines
       do
          read(10,*) iv
          if (iv == -9) exit
          nlinesMetdata = nlinesMetdata + 1
       enddo
       close(10)
       !-----------------------------------------------------------------------

       ! To conserve memory, read met data in blocks of ReadlinesMetdata
       ReadlinesMetdata = int(floor(10000/real(NumberOfGrids,kind(1d0))))
       !write(*,*) 'Met data will be read in chunks of',ReadlinesMetdata,'lines.'

       ! Number of blocks of met data
       ReadMetTimes = int(ceiling(real(nlinesMetdata,kind(1d0))/real(ReadlinesMetdata,kind(1d0))))
       !write(*,*) 'Met data will be read in',ReadMetTimes,'blocks.'

       ! ---- Allocate arrays--------------------------------------------------
       ! N.B. literal values here must exceed the number of required columns!
       allocate(SurfaceChar(NumberOfGrids,MaxNCols_c))	!Surface characteristics	
       allocate(TstepProfiles(NumberOfGrids,6,24*NSH))	!Hourly profiles interpolated to model timestep
       allocate(AHProf_tstep(24*NSH,2))			!AHProf at model timestep
       allocate(WUProfM_tstep(24*NSH,2))		!WUProfM at model timestep
       allocate(WUProfA_tstep(24*NSH,2))		!WUProfA at model timestep
       !! Add water use & snow clearing
       allocate(MetForcingData(1:ReadlinesMetdata,24,NumberOfGrids))   	      !Met forcing data 
       allocate(ModelOutputData(0:ReadlinesMetdata,MaxNCols_cMOD,NumberOfGrids))  !5-min output
       allocate(dataOut(1:ReadlinesMetdata,192,NumberOfGrids))  	!Main output array
       ! ----------------------------------------------------------------------
          
       ! ---- Initialise arrays
       !! Does this need to happen here??
       ! N.B. Each row of SurfaceChar is set to -999 in InitializeSurfaceCharacteristics
       
       !-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       SkippedLines=0  !Initialise lines to be skipped in met forcing file
	
       DO iv=1,ReadMetTimes   !Loop through blocks of met data
 	  
 	  !write(*,*) iv,'/',ReadMetTimes,'ReadMetTimes (iv loop)'

          ! Model calculations are made in two stages: 
          ! (1) initialise the run for each block of met data (iv from 1 to ReadMetTimes)
          ! (2) perform the actual model calculations (SUEWS_Calculations)

          ! (1) First stage: initialise run ------------------------------
          GridCounter=1   !Initialise counter for grids in each year
          DO i=FirstGrid,LastGrid   !Loop through grids
	    	
  	     !write(*,*) i,'/',NumberOfGrids, 'grids (i loop).'
             write(grid_txt,'(I5)') i   !Get grid as a text string
             ! Get met file name for this year for this grid
             FileCodeX=trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(year_txt)
             FileMet=trim(FileInputPath)//trim(FileCodeX)//'_data.txt'
             write(*,*) 'Current FileCode', FileCodeX      
                          	     
  	     ! For the first block of met data ---------------------------
  	     if(iv == 1) then  
	        write(*,*) 'First block of data - doing initialisation'
                ! (a) Transfer characteristics from SiteSelect to correct row of SurfaceChar
                do rr=1,nlinesSiteSelect
                   if(SiteSelect(rr,c_Grid)==i.and.SiteSelect(rr,c_Year)==year_int)then
                      !write(*,*) 'Match found (grid and year) for rr = ', rr
                      call InitializeSurfaceCharacteristics(GridCounter,rr)
		      exit
		   elseif(rr == nlinesSiteSelect) then 
		      write(*,*) 'Program stopped! Year',year_int,'and/or grid',i,'not found in SiteSelect.txt.'
		      !! HCW 20 Nov 2014 - I can't get valueI to write out to the error file correctly
		      call ErrorHint(59,'Cannot find year and/or grid in SiteSelect.txt',i,NotUsed,year_int)
                   endif     
	        enddo
	        ! (b) get initial conditions
                call InitialState(FileCodeX,errFileYes,year_int,GridCounter,year_txt)
             endif   !end first block of met data
                          
             ! For every block of met data -------------------------------
             ! Initialise met forcing data into 3-dimensional matrix
             write(*,*) 'Initialising met data for block',iv
             call SUEWS_InitializeMetData(1)   !Here 1 refers to FileMet unit 
             
             GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid

          ENDDO !end loop over grids
          skippedLines = skippedLines + ReadlinesMetdata   !Increase skippedLines ready for next block

          !Test reading of met data
          !open(78,file='TestingMetData.txt',position="append")
          !do ir=1,(ReadlinesMetdata*NSH)
          !   write(78,*) MetForcingData(ir,1:10,1)
          !enddo
          !close(78)
              
          if((CBLuse==1).or.(CBLuse==2)) call CBL_initial  !These need to be fixed??
          if(SOLWEIGout==1) call SOLWEIG_initial	   !These need to be fixed??

	  write(*,*) 'Initialisation done'
	  ! First stage: initialisation done -----------------------------
      
          ! (2) Second stage: do calculations at 5-min timesteps ---------
          ! First set maximum value of ir
          if(iv == ReadMetTimes) then   !For last block of data in file
             !irMax = nlinesMetdata*NSH - (iv-1)*ReadlinesMetdata*NSH
             irMax = nlinesMetdata - (iv-1)*ReadlinesMetdata
          else 				
           !irMax = ReadlinesMetdata*NSH
           irMax = ReadlinesMetdata
          endif   
          DO ir=1,irMax   !Loop through rows of current block of met data
             GridCounter=1    !Initialise counter for grids in each year
             DO i=FirstGrid,LastGrid   !Loop through grids
           
 	      write(*,*) 'Row (ir):', ir, '/',irMax,'of block (iv):',iv,'Grid:',i	           
 
              ! Call model calculation code
              call SUEWS_Calculations(GridCounter,ir,iv,irMax)
              
              ! Write state information to new InitialConditions files
              if(ir == irMax) then              !If last row...
                 if(iv == ReadMetTimes) then    !...of last block of met data 
                    write(grid_txt,'(I5)') i
                    write(year_txtNext,'(I4)') year_int+1  !Get next year as a string format
                    FileCodeX    =trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(year_txt)
                    FileCodeXNext=trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(year_txtNext)
                    call NextInitial(FileCodeX,year_int,FileCodeXNext)
	         endif
              endif
           
              GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid
              ENDDO !end loop over grids
           
              !!water movements between the grids needs to be taken into account here ??

          ENDDO !end loop over rows of met data

          ! Write output files in blocks --------------------------------
          DO i=FirstGrid,LastGrid
             call SUEWS_Output(i,year_int,iv,irMax)
	     write(*,*) 'Output files written', id,year_int,i
          ENDDO
       !!pause
           
       ENDDO !end loop over blocks of met data
       !-----------------------------------------------------------------------
        
       ! ---- Decallocate arrays ----------------------------------------------
       deallocate(SurfaceChar)
       deallocate(MetForcingData)
       deallocate(ModelOutputData)
       deallocate(dataOut)
       ! ----------------------------------------------------------------------

    ENDDO  !end loop over years

    ! ---- Decallocate array --------------------------------------------------
    ! Daily state needs to be outside year loop to transfer states between years
    deallocate(ModelDailyState)
    ! -------------------------------------------------------------------------

    stop 'finished'

 314 call errorHint(11,trim(FileMet),notUsed,notUsed,ios_out)

 end program SUEWS_Program
