!Main program of SUEWS version 1.0
!Last modified by LJ Nov/2010
!::::::::::::::::::::::::::::::::::::::::::::::::::::
 program SUEWS_Program
    use allocateArray
    use data_In
    use time
    use defaultNotUsed
	 
    IMPLICIT NONE

    integer::NroYears,i,year_int,iyr,iyy                              !Initializing variables used
    character(len=4)::year_txt
    !----------------------------------------------------------
  	stateGrids=NAN 
    soilmoistGrids=NAN
    laiGrids=NAN
    errorChoice=0 !Initialize error file that is printed out

    call overallRunControl   ! Reads RunControl.nml and FunctionalTypes.nml located in SUEWS_initial

    !Open file including yearly and daylighsaving information (needed for NH or SH).
    !Read only the number of years at this point.
    open(98,file=trim(FileInputPath)//trim('ModelledYears.txt'),status='old',err=317)
    READ(98,*,iostat=iostat_var) NroYears  !Read in number of years
    close(98)

    !If only one year is run, no annual file will be created.
    if(NroYears>1)then
	    createAnnual=1
    else
        createAnnual=0
    endif


	FileCodeO=FileCode !Initialize file code from the RunControl.nml
    iyr=0
	do i=1,NroYears
         !Open the file again
         open(98,file=trim(FileInputPath)//trim('ModelledYears.txt'),status='old',err=317)

         !Read empty lines on those years that have already been read in
         do iyy=1,i
          read(98,*)
         enddo
         READ(98,*,iostat=iostat_var)year,DayLightSavingDay(1),DayLightSavingDay(2)
         close(98)

         if (i==1) FirstYear=year
         year_int=int(year)

         write(year_txt,'(I4)') year_int
         FileCode=trim(FileCodeO)//trim(adjustl(year_txt))

         write(*,*)'================== ',trim(fileCode),' =========================='
         call LeapYearCalc (year_int,nofDaysThisYear)
         call SUEWS_spatial(year_int,year_txt,iyr) ! via allocatearray->>stateGrids,soilmoistGrids,laiGrids,laiID

    enddo
    !close(98)
    stop 'finished'

317		call ProblemsText(trim('ModelledYears.txt'))
      	call PauseStop

 end program SUEWS_Program
