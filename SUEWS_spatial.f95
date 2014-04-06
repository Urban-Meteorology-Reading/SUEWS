!The spatial part of SUEWS: LJ in Nov 2010
!Last modified lj/sg May 2012
!----------------------------------------------------------------------------------------------------------
subroutine SUEWS_spatial(year_txt,iyr)

	  use allocateArray
      use data_in
      
      IMPLICIT NONE

      integer::iostat_var,NroGrids,i,ii,jj,iyr ,errFileYes=0
      character(len=15),dimension(2,4)::GridFrom                    !Grid connections for each grid    
      real (kind(1d0)),DIMENSION(4)::GridFromFrac                   !Fraction of water moving between these
      character(len=15)::grid
      character(len=4)::year_txt
      character(len=100)::GridConnectionsName
       
      real (kind(1d0)),DIMENSION(100,7)::SnowPackG !Max numer of grids set to 100!!

      
      GridConnectionsFrac=0     
      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      !Open Spatial SUEWS information file GridConnectionsYYYY.txt
	  GridConnectionsName=trim(FileInputPath)//'GridConnections'//trim(adjustl(year_txt))//'.txt'
      open(89,file=trim(GridConnectionsName),status='old',err=317)
	  READ(89,*,iostat=iostat_var) NroGrids  !Read number of grids
      
      
      !Grid names and connections are read in
      do i=1,NroGrids
      	 READ(89,*,iostat=iostat_var) GridConnections(1,i),GridConnectionsFrac(i),GridConnections(2,i)
      enddo
      close (89)
      
	  !=======Call the program itself grid by grid========================
      Grid='none' !Initiate Grid number
      
      do i=1,NroGrids
         if (trim(GridConnections(1,i))/=trim(Grid)) then
             Grid=GridConnections(1,i) !Define grid
             GridFromFrac=0
                
			 !Check what grids get water from other grids 
             jj=1
             do ii=1, NroGrids 
                if (trim(GridConnections(2,ii))==trim(Grid)) then
                   GridFrom(1,jj)=GridConnections(1,ii)
                   GridFromFrac(jj)=GridConnectionsFrac(ii)
                   jj=jj+1
                endif
             enddo

    		 call OpenAnnualFiles(Grid) !Open annual and dailystate output files
             
             !Define new filecode for each grid separately and save to FileChoices
             FileCode=trim(Grid)//trim(year_txt)
             write(12,*) "*************",FileCode,"*************"
             
			 call RunControlByGridByYear  !Call grid specific runcontrol        
	         Call InitialState(Grid,errFileYes)      !Initial state of the run            
             call SUEWS_temporal(Grid,GridFrom,GridFromFrac,iyr,errFileYes,SnowPackG(i,:)) !Call the code grid by grid
             
          endif
	  enddo
   
      return
317	  call ProblemsText(trim( GridConnectionsName))
      call PauseStop	


end subroutine SUEWS_spatial