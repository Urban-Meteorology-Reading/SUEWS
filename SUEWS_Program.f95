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
!Last modified by HCW 03 Mar 2015
!Last modified by HCW 26 Feb 2015
!Last modified by HCW 03 Dec 2014
!
! To do:
! 	- Snow modules need updating for water balance part
! 	- Water movement between grids (GridConnections) not yet coded
!	- Check all arrays are allocated/deallocated and initialised in the correct place
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

    character(len=4)::  year_txt,&	!Year as a text string
    			year_txtNext	!Following year as a text string (used for NextInitial)
    character(len=20):: FileCodeX,&	!Current file code
    			FileCodeXNext	!File code for the following year  
    character(len=20):: grid_txt	!Grid number as a text string (from FirstGrid to LastGrid)
    
    integer:: nlinesLimit   !Max number of lines that can be read in one go for each grid
                               
    integer:: 	i,&	   ! Grid number (from FirstGrid to LastGrid)
    		iv,&	   ! Block number (from 1 to ReadBlocksMetData)
    		ir,irMax,& ! Row number within each block (from 1 to irMax)
    		rr,& 	   ! Row of SiteSelect corresponding to current year and grid
                year_int   ! Year as an integer (from SiteSelect rather than met forcing file) 
    
    integer:: errFileYes            
    
    logical:: PrintPlace=.true.   !Prints row, block, and grid number to screen if TRUE
    		     
    !==========================================================================
    
    ! Initialise error file (0 -> problems.txt file is created)
    errorChoice=0

    ! Read RunControl.nml and all input files from SiteInfo spreadsheet
    call overallRunControl

    ! First find first and last year of the current run
    FirstYear = minval(int(SiteSelect(:,c_Year)))
    LastYear  = maxval(int(SiteSelect(:,c_Year)))
    
    ! Find the first and last grid numbers (N.B. need to have the same grids for each year)
    FirstGrid = minval(int(SiteSelect(:,c_Grid))) 
    LastGrid  = maxval(int(SiteSelect(:,c_Grid))) 
    NumberOfGrids = LastGrid-FirstGrid+1   !Number of grids of the current run
    if(NumberOfGrids > MaxNumberOfGrids) then
       call ErrorHint(64,'No. of grids exceeds max. possible no. of grids.',real(MaxNumberOfGrids,kind(1d0)),NotUsed,NumberOfGrids)
    endif    
    
    write(*,*) '--------------------------------------------'
    write(*,*) 'Years identified:',FirstYear,'to',LastYear
    write(*,*) 'Grids identified:',FirstGrid,'to',LastGrid

    ! ---- Allocate arrays ----------------------------------------------------
    ! Daily state needs to be outside year loop to transfer states between years
    allocate(ModelDailyState(NumberOfGrids,MaxNCols_cMDS))   !DailyState        
    allocate(DailyStateFirstOpen(NumberOfGrids))             !Initialization for header
    ! ---- Initialise arrays --------------------------------------------------
    ModelDailyState(:,:) = -999
    DailyStateFirstOpen(:) = 1
    ! -------------------------------------------------------------------------
    
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
       nlinesMetdata = 0   !Initialise nlinesMetdata (total number of lines in met forcing file)
       do
          read(10,*) iv
          if (iv == -9) exit
          nlinesMetdata = nlinesMetdata + 1
       enddo
       close(10)
       !-----------------------------------------------------------------------
       
       ! To conserve memory, read met data in blocks
       ! Find number of lines that can be read in each block (i.e. read in at once)
       ReadLinesMetData = nlinesMetData   !Initially set limit as the size of the met file (N.B.solves problem with Intel fortran)
       nlinesLimit = int(floor(MaxLinesMet/real(NumberOfGrids,kind(1d0))))
       if(nlinesMetData > nlinesLimit) then   !But restrict if this limit exceeds memory capacity
          ReadLinesMetData = nlinesLimit
       endif    
       !write(*,*) 'Met data will be read in chunks of',ReadlinesMetdata,'lines.'

       ! Find number of blocks of met data
       ReadBlocksMetData = int(ceiling(real(nlinesMetData,kind(1d0))/real(ReadLinesMetData,kind(1d0))))
       !write(*,*) 'Met data will be read in',ReadBlocksMetData,'blocks.'

       ! ---- Allocate arrays--------------------------------------------------
       allocate(SurfaceChar(NumberOfGrids,MaxNCols_c))   !Surface characteristics
       allocate(MetForcingData(1:ReadlinesMetdata,ncolumnsMetForcingData,NumberOfGrids))   !Met forcing data 
       allocate(ModelOutputData(0:ReadlinesMetdata,MaxNCols_cMOD,NumberOfGrids))           !Data at model timestep
       allocate(dataOut(1:ReadlinesMetdata,ncolumnsDataOut,NumberOfGrids))  	           !Main output array
       if (SOLWEIGuse == 1) then
          allocate(dataOutSOL(1:ReadlinesMetdata,28,NumberOfGrids))                        !SOLWEIG POI output
       endif
       if (CBLuse >= 1) then
          allocate(dataOutBL(1:ReadlinesMetdata,22,NumberOfGrids))  	                   !CBL output
       endif
       
       allocate(TstepProfiles(NumberOfGrids,6,24*NSH))	!Hourly profiles interpolated to model timestep
       allocate(AHProf_tstep(24*NSH,2))			!Anthropogenic heat profiles at model timestep
       allocate(WUProfM_tstep(24*NSH,2))		!Manual water use profiles at model timestep
       allocate(WUProfA_tstep(24*NSH,2))		!Automatic water use profiles at model timestep
       !! Add snow clearing (?)      
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
          DO i=FirstGrid,LastGrid   !Loop through grids
	    	
  	     !write(*,*) i,'/',NumberOfGrids, 'grids (i loop).'
             write(grid_txt,'(I5)') i   !Get grid as a text string
             ! Get met forcing file name for this year for this grid
             FileCodeX=trim(FileCode)//trim(adjustl(grid_txt))//'_'//trim(year_txt)
             write(*,*) 'Current FileCode', FileCodeX      
                          	     
  	     ! For the first block of met data --------------------------------
  	     if(iv == 1) then  
	        !write(*,*) 'First block of data - doing initialisation'
                ! (a) Transfer characteristics from SiteSelect to correct row of SurfaceChar
                do rr=1,nlinesSiteSelect
                   if(SiteSelect(rr,c_Grid)==i.and.SiteSelect(rr,c_Year)==year_int)then
                      !write(*,*) 'Match found (grid and year) for rr = ', rr
                      call InitializeSurfaceCharacteristics(GridCounter,rr)
		      exit
		   elseif(rr == nlinesSiteSelect) then 
		      write(*,*) 'Program stopped! Year',year_int,'and/or grid',i,'not found in SiteSelect.txt.'
		      call ErrorHint(59,'Cannot find year and/or grid in SiteSelect.txt',real(i,kind(1d0)),NotUsed,year_int)
                   endif     
	        enddo
	        ! (b) get initial conditions
                call InitialState(FileCodeX,errFileYes,year_int,GridCounter,year_txt)
             endif   !end first block of met data
                          
             ! For every block of met data ------------------------------------
             ! Initialise met forcing data into 3-dimensional matrix
             write(*,*) 'Initialising met data for block',iv
             
             if(MultipleMetFiles == 1) then   !If each grid has its own met file
                FileMet=trim(FileInputPath)//trim(FileCodeX)//'_data.txt'
                call SUEWS_InitializeMetData(1)   
             else                             !If one met file used for all grids  
                FileMet=trim(FileInputPath)//trim(FileCode)//'_'//trim(year_txt)//'_data.txt'
                if(i == FirstGrid) then       !Read for the first grid only  
                   call SUEWS_InitializeMetData(1)
                else                          !Then for subsequent grids simply copy data  
                   MetForcingData(1:ReadlinesMetdata,1:24,GridCounter) = MetForcingData(1:ReadlinesMetdata,1:24,1)
                endif
             endif

             GridCounter = GridCounter+1   !Increase GridCounter by 1 for next grid

          ENDDO !end loop over grids
          skippedLines = skippedLines + ReadlinesMetdata   !Increase skippedLines ready for next block

          ! Initialise CBL and SOLWEIG parts if required
          if((CBLuse==1).or.(CBLuse==2)) call CBL_ReadInputData
          if(SOLWEIGuse==1) call SOLWEIG_initial

	  !write(*,*) 'Initialisation done'
	  ! First stage: initialisation done ----------------------------------
      
          ! (2) Second stage: do calculations at 5-min timesteps --------------
          ! First set maximum value of ir
          if(iv == ReadBlocksMetData) then   !For last block of data in file
             irMax = nlinesMetdata - (iv-1)*ReadLinesMetdata
          else 				
             irMax = ReadLinesMetdata
          endif   
          DO ir=1,irMax   !Loop through rows of current block of met data
             GridCounter=1    !Initialise counter for grids in each year
             DO i=FirstGrid,LastGrid   !Loop through grids
           
 	      if(PrintPlace) write(*,*) 'Row (ir):', ir,'/',irMax,'of block (iv):', iv,'/',ReadBlocksMetData,'Grid:',i	           
 
              ! Call model calculation code
              call SUEWS_Calculations(GridCounter,ir,iv,irMax)
              
              ! Write state information to new InitialConditions files
              if(ir == irMax) then              !If last row...
                 if(iv == ReadBlocksMetData) then    !...of last block of met data 
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
	     write(*,*) 'Output files written (reached DOY, year, Grid)', id,year_int,i
          ENDDO
                  
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
