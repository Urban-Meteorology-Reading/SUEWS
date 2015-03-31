 subroutine ErrorHint(errh,ProblemFile,value,value2,valueI) ! real
!errh        -- Create a numbered code for the situation so get a unique message to help solve the problem
!ProblemFile -- Filename where the problem occurs
!value       -- Real number with correct type
!value2      -- Real number with correct type
!valueI      -- Integer 2
!Last modified LJ 8 Feb 2013
! sg 29/7/14 - close (500)
! LJ 2/10/2014 - addition of comments
!-------------------------------------------------------------------------------------------------

 use defaultNotUsed

 IMPLICIT NONE

 real(kind(1d0)):: value,value2

 character (len=*)::ProblemFile                !Name of the problem file
 character (len=150)::text1='unknown problem'   !Initialization of text
 integer:: errh,ValueI,ValueI2					! v7,v8 initialised as false, HCW 28/10/2014
 logical:: v1=.false.,v2=.false.,v3=.false.,v4=.false.,v5=.false.,v6=.false.,v7=.false.,v8=.false.
 logical:: returnTrue=.false.

 ! Initialise returnTrue as false (HCW 29/10/2014)
 ! - need to do this in Fortran as values assigned in declarations are not applied 
 ! on subsequent calling of the subroutine
 returnTrue=.false.
 ! Initialise v1-v8 as false
 v1=.false.
 v2=.false.
 v3=.false.
 v4=.false.
 v5=.false.
 v6=.false.
 v7=.false.
 v8=.false.

 call ProblemsText(ProblemFile)                 !Call the subroutine that opens the problem.txt file

 !The list of knows possible problems of the code:
 !  text1 is the error message written to the ProblemFile.
 !  v1 -v7 are different possibilities for what numbers will be written out
 !  ReturnTrue is true if the model run can continue. (Comment modified by HCW 17/10/2014) 
 if(errh==1)then
    text1='for FAIBLDG - for z0_method selected - this value looks inappropriate'
    v1=.true.
 elseif(errh==2)then
    text1='for FAITree - for z0_method selected - this value looks inappropriate'
    v1=.true.
 elseif(errh==3) then
    text1='sdec1=sdec2 - check/adjust GDD and SDDFull & DailyState'
    v3=.true.
 elseif(errh==4) then
    text1='sdec3=sdec=4 - check/adjust GDD and SDDFull & DailyState'
    v3=.true.
 elseif(errh==5)then
    text1='Value for z0 in SiteSelect looks inappropriate'
    v1=.true.
 elseif(errh==6)then
    text1='Value for zd in SiteSelect looks inappropriate'
    v1=.true.
 elseif(errh==7) then
    text1='RA value exceeds thresholds set in model'
    v1=.true.
    returnTrue=.true.
 elseif(errh==8) then
    text1='Should be zero (water cannot move from one surface to the same surface).'
    v1=.true.
 elseif(errh==9) then
    text1= 'One of these (water to Runoff or to SoilStore) should be zero.'
    v2=.true.
 elseif(errh==10) then
    text1='Should sum to 1.'
    v1=.true.
 elseif(errh==11) then
    text1=' this file not found- ios_out'
    v3=.true.
 elseif(errh==12) then
    text1= 'daywatper values - need a space in front after = if less than 1 '
    v3=.true.
    returnTrue=.true.
 elseif(errh==13) then
    text1= ' problem with line, ios_out error: '
    v6=.true.
 elseif(errh==14) then
    text1= ' Above subroutine has calculated of z0m - this value looks inappropriate'
    v1=.true.
 elseif(errh==15) then
    text1= ' Above subroutine has calculated of zd - this value looks inappropriate'
    v1=.true.
 elseif(errh==16) then
    text1=' Check gridConnectionsYYYY.txt file to see if year embedded in name, ios_out'
    v3=.true.
    returnTrue=.true.
 elseif(errh==17) then
    text1= 'LOG(zzd/z0M)<0.001) problem: zzd, z0m'
    v2=.true.
    returnTrue=.true.
 elseif(errh==18) then
    text1='check the soil depth relative to: SoilStoreCap(is), soilmoist(is), whichsurface no'
    v4=.true.
    !probably will stop but need to check this is really the problem
    returnTrue=.true.
 elseif(errh==19)then
    text1=' cp, press , lv - not stopping after pause'
    v4=.true.
    returnTrue=.true.
 elseif(errh==20)then
    text1=' skip lines, ios_out If running multiple grids, check the order of the grids.'
    v5=.true.
 elseif(errh==21)then
    text1='Check qn, qn1_S, qf.'
    v1=.true.
          !text1=' missing times in GIS file: id- GIS it-gis, it-met'  !gis file no longer used - errh 21 recycled
	  !v4=.true.	  
 elseif(errh==22)then
    text1=' QH_observed, QE_observed, QH_choice: '
    v4=.true.
 elseif(errh==23)then
    text1='CBL-sonde -need to increase size of izm:zmax,izm'
    v5=.true.
 elseif(errh==24) then
    text1='CBL file problem - opening'
    v8=.true.
 elseif(errh==25) then
    text1='CBL file problem -- reading sonde data, line:'
    v3=.true.
 elseif(errh==26) then
    text1='Check that FileOutputPath and FileCode are specified and have double quotes around them'
    v8=.true.
 elseif(errh==27)then
    text1='Problems with Met data -forcing data: variable value, dectime'
    v2=.true.  ! 2 real
 elseif(errh==28) then
    text1='Processing in subroutine indicated has a problem, variables'
    returntrue=.true.
    v3=.true.  ! 1 integer
 elseif(errh==29) then
    text1='Processing in subroutine indicated has a problem, time, variables'
    returntrue=.true.
    v7=.true.  ! 1 real, 2 integers
 elseif(errh==30) then
    text1='Processing in subroutine indicated has a problem, time, variables'
    returntrue=.true.
    v2=.true.  ! 2 real
 elseif(errh==31) then
    text1='Processing in subroutine indicated has a problem, time, variables'
    returntrue=.true.
    v1=.true.  ! 1 real
 elseif(errh==32) then
    text1='Model applicable to local scale, z<z0d'
    v2=.true.  ! 2 real
 elseif(errh==33) then
    text1 = 'Number of snow layers too large.'
    v1=.true.  ! 1 real
 elseif(errh==34) then
    text1= 'Air temperature > 55 C -- will keep running'
    v1=.true.  ! 1 real
    returntrue=.true.
 elseif(errh==35) then
    text1 = 'Problems with Met data -forcing data: doy, dectime'
    v2 = .true.  ! 2 real
 elseif(errh==36) then
    text1 = 'Problems in the InitialConditions file: Initial value, compared value'
    v2 = .true.  !2 real
 elseif(errh==37) then
    text1 = 'Problems in the InitialConditions file: Initial value, compared value'
    returntrue=.true.
    v2 = .true.  !2 real
 elseif(errh==38) then
    text1 = 'H=(qn*0.2)/(avdens*avcp)'
    returntrue=.true.
    v1 = .true.  !2 real
 elseif(errh==39) then
    text1 = 'Different value of TSTEP needed (300 s recommended). Resolution of met data must match TSTEP set in RunControl.'
    v2 = .true.  !2 real
 elseif(errh==40) then
    text1='SOLWEIG file problem - opening'
    v8=.true.
 elseif(errh==41) then
    text1= ' addwaterbody= Error1-- but watersurf=  Error 2'
    v2=.true. !2 real
 elseif(errh==42)then
    text1= 'abs(rho_d)<0.001000.OR.abs(rho_v)<0.001000.OR.abs(rho_d+rho_v)<0.001000) rho_v,rho_d, T'
    returntrue=.true.
    v4=.true. !2 real, temperature as an integer
 elseif(errh==43) then
    text1='Switching Years - will keep running'
    returntrue=.true.
    v8=.true.
 elseif(errh==44)then
    text1='Initial File Name - will keep going'
    returntrue=.true.
    v8=.true.
 elseif(errh==45)then
    text1='Pressure < 900 hPa, Loop Number'
    returntrue=.true.
    v5 = .true.
 elseif(errh==46)then
    text1 = 'Pressure < 900 hPa'
    returntrue = .true.
    v1 = .true.
 elseif(errh==47)then
    text1 = 'File missing'
    returntrue = .true.
 elseif(errh==48)then
    text1 = 'Something wrong in the rows of the file'
    returntrue = .true.
 elseif(errh==49)then
    text1 = 'Problems in saving to InitialConditionsYYYY.nml'
 elseif(errh==50)then
    text1 = 'Wrong number of lines read: nsurf, [-1 EOF; -2 EOR]'
    v1 = .true.
 elseif(errh==51)then
    text1 = 'Problems in opening the file'
    write(*,*) ProblemFile
 elseif(errh==52)then
    text1 = 'Problems in opening the output file'
 elseif(errh==53)then
    text1 = 'AH_min=0.and.Ah_slope=0.and.T_Critic=0, AnthropHeatChoice='
    returntrue = .true.
    v3 = .true.
 elseif(errh==54)then
    text1 = 'QF_A=0.and.QF_B=0.and.QF_C=0, AnthropHeatChoice='
    returntrue = .true.
    v3 = .true.
 elseif(errh==55)then
    text1 = 'InputmetFormat='
    returntrue = .true.
    v3 = .true.
 elseif(errh==56)then
    text1 = 'Check input files against manual (N.B. Case sensitive).'
    v8 = .true.
 elseif(errh==57)then
    text1 = 'not found. Check input files.'
    v1 = .true.
 elseif(errh==58)then
    text1 = 'File header not specified in model code.'
    v8 = .true.
 elseif(errh==59)then
    text1 = 'not found. Check SUEWS_SiteSelect.txt.'
    v6 = .true.
 elseif(errh==60)then
    text1 = 'non-unique code.'
    v1 = .true.
 elseif(errh==61)then
    text1 = 'Check coefficients and drainage equation specified in input files.'
    v4 = .true.
    returntrue = .true.
 elseif(errh==62)then
    text1 = 'Problem with soil moisture calculation.'
    v5 = .true.
 elseif(errh==63)then
    text1 = 'Problem with calculation.'
    v1 = .true.
 elseif(errh==64)then
    text1 = 'SUEWS cannot currently handle this many grids.'
    v6 = .true.
 endif
 !---------------------------------------------------------------------
 !This part of the code determines how the error message is written out
    
 if(v1) then ! 1 real
    write(500,*)'ERROR value: =', value
 elseif(v2) then ! 2 real
    write(500,*)'ERROR values: =', value, value2
 elseif(v3) then ! 1 integer
    write(500,*)'ERROR value: =', valueI
 elseif(v4) then ! 2 real, 1 integer
    write(500,*)'ERROR values: =', value, value2, valueI
 elseif(v5) then ! 1 real 1 integer
    write(500,*)'ERROR values: =', value, valueI
 elseif(v6) then ! 2 integer
    valueI2=int(value)
    write(500,*)'ERROR values: =', valueI, valueI2
 elseif(v7) then
    valueI2=int(value2)
    write(500,*)'ERROR values: =', value, valueI2, valueI
 elseif(v8) then
    ! no error values
 endif
     
     
 !Write the actual comment the problems file
 write(500,*) trim(text1)

 close(500)

 !When returnTrue=true, then the program can continue despite the errors seen
 if(returnTrue) then
     !write(*,*)'Problems.txt has been closed and overwritten if other errors occur'
    return  !Continue program
 endif

 call PauseStop(ProblemFile)        !Stop the program
      
 return
 end subroutine ErrorHint

!=============================================================

 subroutine ProblemsText(ProblemFile)

    use defaultNotUsed
    IMPLICIT NONE

    character (len=*):: ProblemFile

    !Opening problems.txt file: First option is selected if the file is opened for the first time
    !Second option for later points
    if (errorChoice==0) then
        open(500,file='problems.txt')
        write(*,*) 'See problems.txt for possible issues in the run.'
        errorChoice=1
    else
        open(500,file='problems.txt',position="append")
    endif

    !Writing of the problem file
    write(500,*)'problem file: ',trim(ProblemFile)

    return
 end subroutine ProblemsText


 subroutine PauseStop(ProblemFile)

   IMPLICIT NONE
   character (len=*):: ProblemFile

   write(*,*)'problem file: ',trim(ProblemFile)
   write(*,*)'see problems.txt'

   !pause
   stop
 end subroutine PauseStop
