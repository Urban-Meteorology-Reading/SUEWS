!Main program of SUEWS version 1.0
!Last modified by LJ Nov/2010
!::::::::::::::::::::::::::::::::::::::::::::::::::::
program SUEWS_Program
	  use allocateArray
    use data_In
    use time
    use defaultNotUsed
	 
	  implicit none
    integer::NroYears,i,year_int,iyr                              !Initializing variables used
	  character(len=4)::year_txt

  	stateGrids=NAN 
    soilmoistGrids=NAN
    laiGrids=NAN

    
    call overallRunControl   ! Reads RunControl.nml and FunctionalTypes.nml located in SUEWS_initial
	    
    !Open file including yearly information
    open(98,file=trim(FileInputPath)//trim('ModelledYears.txt'),status='old',err=317)
    READ(98,*,iostat=iostat_var) NroYears  !Read number of years
    
    !Use daylighsavings to determine if NH or SH
    READ(98,*,iostat=iostat_var)year,DayLightSavingDay(1),DayLightSavingDay(2)
    
    FirstYear=int(year) !Define the first year
    
    rewind(98)
    READ(98,*,iostat=iostat_var) NroYears  !Read number of years
    if(NroYears>1)then
	    createAnnual=1
    else
        createAnnual=0
    endif

	FileCodeO=FileCode
    iyr=0
	do i=1,NroYears
         READ(98,*,iostat=iostat_var)year,DayLightSavingDay(1),DayLightSavingDay(2)
         year_int=int(year)
         write(year_txt,'(I4)')year_int
         FileCode=trim(FileCodeO)//trim(adjustl(year_txt))
   
       	 write(12,*)'================== ',trim(fileCode),' =========================='	
         write(*,*)'================== ',trim(fileCode),' =========================='
         call LeapYearCalc (year,nofDaysThisYear)
         call SUEWS_spatial(year_txt,iyr) ! via allocatearray->>stateGrids,soilmoistGrids,laiGrids,laiID

    enddo
    close(12) 
    stop 'finished'

317		call ProblemsText(trim('ModelledYears.txt'))
      	call PauseStop

end program SUEWS_Program
