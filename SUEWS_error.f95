SUBROUTINE ErrorHint(errh,ProblemFile,VALUE,value2,valueI) ! real
  !errh        -- Create a numbered code for the situation so get a unique message to help solve the problem
  !ProblemFile -- Filename where the problem occurs
  !value       -- Real number with correct type
  !value2      -- Real number with correct type
  !valueI      -- Integer 2
  !Last modified LJ 8 Feb 2013
  ! sg 29/7/14 - close (500)
  ! LJ 2/10/2014 - addition of comments
  ! HCW 25 May 2016 Added warning/error labels to distinguish serious errors (that stop program)

  !-------------------------------------------------------------------------------------------------

  USE defaultNotUsed

  IMPLICIT NONE

  REAL(KIND(1d0)):: VALUE,value2

  CHARACTER (len=*)::ProblemFile                 ! Name of the problem file
  CHARACTER (len=150)::text1='unknown problem'   ! Initialization of text
  INTEGER:: errh,ValueI,ValueI2                  ! v7,v8 initialised as false, HCW 28/10/2014
  LOGICAL:: v1=.FALSE.,v2=.FALSE.,v3=.FALSE.,v4=.FALSE.,v5=.FALSE.,v6=.FALSE.,v7=.FALSE.,v8=.FALSE.
  LOGICAL:: returnTrue=.FALSE.

  ! Initialise returnTrue as false (HCW 29/10/2014)
  ! - need to do this in Fortran as values assigned in declarations are not applied
  ! on subsequent calling of the subroutine
  returnTrue=.FALSE.
  ! Initialise v1-v8 as false
  v1=.FALSE.
  v2=.FALSE.
  v3=.FALSE.
  v4=.FALSE.
  v5=.FALSE.
  v6=.FALSE.
  v7=.FALSE.
  v8=.FALSE.

  CALL ProblemsText(ProblemFile)                 !Call the subroutine that opens the problem.txt file

  !The list of knows possible problems of the code:
  !  text1 is the error message written to the ProblemFile.
  !  v1 -v7 are different possibilities for what numbers will be written out
  !  ReturnTrue is true if the model run can continue. (Comment modified by HCW 17/10/2014)
  IF(errh==1)THEN
     text1='for FAIBLDG - for z0_method selected - this value looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==2)THEN
     text1='for FAITree - for z0_method selected - this value looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==3) THEN
     text1='sdec1=sdec2 - check/adjust GDD and SDDFull & DailyState'
     v3=.TRUE.
  ELSEIF(errh==4) THEN
     text1='sdec3=sdec=4 - check/adjust GDD and SDDFull & DailyState'
     v3=.TRUE.
  ELSEIF(errh==5)THEN
     text1='Value for z0 in SiteSelect looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==6)THEN
     text1='Value for zd in SiteSelect looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==7) THEN
     text1='RA value exceeds thresholds set in model'
     v1=.TRUE.
     returnTrue=.TRUE.
  ELSEIF(errh==8) THEN
     text1='Should be zero (water cannot move from one surface to the same surface).'
     v1=.TRUE.
  ELSEIF(errh==9) THEN
     text1= 'One of these (water to Runoff or to SoilStore) should be zero.'
     v2=.TRUE.
  ELSEIF(errh==10) THEN
     text1='Should sum to 1.'
     v1=.TRUE.
  ELSEIF(errh==11) THEN
     text1=' this file not found- ios_out'
     v3=.TRUE.
  ELSEIF(errh==12) THEN
     text1= 'daywatper values - need a space in front after = if less than 1 '
     v3=.TRUE.
     returnTrue=.TRUE.
  ELSEIF(errh==13) THEN
     text1= ' problem with line, ios_out error: '
     v6=.TRUE.
  ELSEIF(errh==14) THEN
     text1= ' Above subroutine has calculated of z0m - this value looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==15) THEN
     text1= ' Above subroutine has calculated of zd - this value looks inappropriate'
     v1=.TRUE.
  ELSEIF(errh==16) THEN
     text1=' Check gridConnectionsYYYY.txt file to see if year embedded in name, ios_out'
     v3=.TRUE.
     returnTrue=.TRUE.
  ELSEIF(errh==17) THEN
     text1= 'LOG(zzd/z0M)<0.001) problem: zzd, z0m'
     v2=.TRUE.
     returnTrue=.TRUE.
  ELSEIF(errh==18) THEN
     text1='check the soil depth relative to: SoilStoreCap(is), soilmoist(is), whichsurface no'
     v4=.TRUE.
     !probably will stop but need to check this is really the problem
     returnTrue=.TRUE.
  ELSEIF(errh==19)THEN
     text1=' cp, press , lv - not stopping after pause'
     v4=.TRUE.
     returnTrue=.TRUE.
  ELSEIF(errh==20)THEN
     text1=' skip lines, ios_out If running multiple grids, check the order of the grids.'
     v5=.TRUE.
  ELSEIF(errh==21)THEN
     text1='Check qn, qn1_S, qf.'
     v1=.TRUE.
     !text1=' missing times in GIS file: id- GIS it-gis, it-met'  !gis file no longer used - errh 21 recycled
     !v4=.true.
  ELSEIF(errh==22)THEN
     text1=' QH_observed, QE_observed, QH_choice: '
     v4=.TRUE.
  ELSEIF(errh==23)THEN
     text1='CBL-sonde -need to increase size of izm:zmax,izm'
     v5=.TRUE.
  ELSEIF(errh==24) THEN
     text1='CBL file problem - opening'
     v8=.TRUE.
  ELSEIF(errh==25) THEN
     text1='CBL file problem -- reading sonde data, line:'
     v3=.TRUE.
  ELSEIF(errh==26) THEN
     text1='Check that FileOutputPath and FileCode are specified and have double quotes around them'
     v8=.TRUE.
  ELSEIF(errh==27)THEN
     text1='Problems with Met data -forcing data: variable value, dectime'
     v2=.TRUE.  ! 2 real
  ELSEIF(errh==28) THEN
     text1='Processing in subroutine indicated has a problem, variables'
     returntrue=.TRUE.
     v3=.TRUE.  ! 1 integer
  ELSEIF(errh==29) THEN
     text1='Processing in subroutine indicated has a problem, time, variables'
     returntrue=.TRUE.
     v7=.TRUE.  ! 1 real, 2 integers
  ELSEIF(errh==30) THEN
     text1='Processing in subroutine indicated has a problem, time, variables'
     returntrue=.TRUE.
     v2=.TRUE.  ! 2 real
  ELSEIF(errh==31) THEN
     text1='Processing in subroutine indicated has a problem, time, variables'
     returntrue=.TRUE.
     v1=.TRUE.  ! 1 real
  ELSEIF(errh==32) THEN
     text1='Model applicable to local scale, z<z0d'
     v2=.TRUE.  ! 2 real
  ELSEIF(errh==33) THEN
     text1 = 'Number of snow layers too large.'
     v1=.TRUE.  ! 1 real
  ELSEIF(errh==34) THEN
     text1= 'Air temperature > 55 C -- will keep running'
     v1=.TRUE.  ! 1 real
     returntrue=.TRUE.
  ELSEIF(errh==35) THEN
     text1 = 'Problems with Met data -forcing data: doy, dectime'
     v2 = .TRUE.  ! 2 real
  ELSEIF(errh==36) THEN
     text1 = 'Problems in the InitialConditions file: Initial value, compared value'
     v2 = .TRUE.  !2 real
  ELSEIF(errh==37) THEN
     text1 = 'Problems in the InitialConditions file: Initial value, compared value'
     returntrue=.TRUE.
     v2 = .TRUE.  !2 real
  ELSEIF(errh==38) THEN
     text1 = 'H=(qn*0.2)/(avdens*avcp)'
     returntrue=.TRUE.
     v1 = .TRUE.  !2 real
  ELSEIF(errh==39) THEN
     text1 = 'Different value of TSTEP needed (300 s recommended). Resolution of forcing data must match TSTEP set in RunControl.'
     v4 = .TRUE.  !2 real, 1 int
  ELSEIF(errh==40) THEN
     text1='SOLWEIG file problem - opening'
     v8=.TRUE.
  ELSEIF(errh==41) THEN
     text1= ' addwaterbody= Error1-- but watersurf=  Error 2'
     v2=.TRUE. !2 real
  ELSEIF(errh==42)THEN
     text1= 'abs(rho_d)<0.001000.OR.abs(rho_v)<0.001000.OR.abs(rho_d+rho_v)<0.001000) rho_v,rho_d, T'
     returntrue=.TRUE.
     v4=.TRUE. !2 real, temperature as an integer
  ELSEIF(errh==43) THEN
     text1='Switching Years - will keep running'
     returntrue=.TRUE.
     v8=.TRUE.
  ELSEIF(errh==44)THEN
     text1='Initial File Name - will keep going'
     returntrue=.TRUE.
     v8=.TRUE.
  ELSEIF(errh==45)THEN
     text1='Pressure < 900 hPa, Loop Number'
     returntrue=.TRUE.
     v5 = .TRUE.
  ELSEIF(errh==46)THEN
     text1 = 'Pressure < 900 hPa'
     returntrue = .TRUE.
     v1 = .TRUE.
  ELSEIF(errh==47)THEN
     text1 = 'File missing'
     !returntrue = .TRUE.
  ELSEIF(errh==48)THEN
     text1 = 'Something wrong in the rows of the file'
     !returntrue = .TRUE.
  ELSEIF(errh==49)THEN
     text1 = 'Problems in saving to InitialConditionsYYYY.nml'
  ELSEIF(errh==50)THEN
     text1 = 'Wrong number of lines read: nsurf, [-1 EOF; -2 EOR]'
     v1 = .TRUE.
  ELSEIF(errh==51)THEN
     text1 = 'Problems in opening the file'
     WRITE(*,*) ProblemFile
  ELSEIF(errh==52)THEN
     text1 = 'Problems in opening the output file'
  ELSEIF(errh==53)THEN
     text1 = 'AH_min=0.and.Ah_slope=0.and.T_Critic=0, AnthropHeatChoice='
     returntrue = .TRUE.
     v3 = .TRUE.
  ELSEIF(errh==54)THEN
     text1 = 'QF_A=0.and.QF_B=0.and.QF_C=0, AnthropHeatChoice='
     returntrue = .TRUE.
     v3 = .TRUE.
  ELSEIF(errh==55)THEN
     text1 = 'InputmetFormat='
     returntrue = .TRUE.
     v3 = .TRUE.
  ELSEIF(errh==56)THEN
     text1 = 'Check input files against manual (N.B. Case sensitive).'
     v8 = .TRUE.
  ELSEIF(errh==57)THEN
     text1 = 'not found. Check input files.'
     v1 = .TRUE.
  ELSEIF(errh==58)THEN
     text1 = 'File header not specified in model code.'
     v8 = .TRUE.
  ELSEIF(errh==59)THEN
     text1 = 'not found. Check SUEWS_SiteSelect.txt.'
     v6 = .TRUE.
  ELSEIF(errh==60)THEN
     text1 = 'non-unique code.'
     v1 = .TRUE.
  ELSEIF(errh==61)THEN
     text1 = 'Check coefficients and drainage equation specified in input files.'
     v4 = .TRUE.
     returntrue = .TRUE.
  ELSEIF(errh==62)THEN
     text1 = 'Problem with soil moisture calculation.'
     v5 = .TRUE.
  ELSEIF(errh==63)THEN
     text1 = 'Problem with calculation.'
     v1 = .TRUE.
  ELSEIF(errh==64)THEN
     text1 = 'SUEWS cannot currently handle this many grids.'
     v6 = .TRUE.
  ELSEIF(errh==65) THEN
     text1='Negative gs calculated! Check suitability of parameters in Conductance.txt.'
     returntrue=.TRUE.
     v7=.TRUE.  ! 1 real, 2 integers
  ELSEIF(errh==66)THEN
     text1 = 'Different number of lines in ESTM forcing and Met forcing files.'
     v6 = .TRUE.
  ELSEIF(errh==67) THEN
     text1='ESTMClass1 automatically set to 100%.'
     returntrue=.TRUE.
     v1=.TRUE.
  ELSEIF(errh==68) THEN
     text1='Initial Bowen ratio automatically set to 1.'
     returntrue=.TRUE.
     v1=.TRUE.
  ELSEIF(errh==69) THEN
     text1='Setting QF_traff to zero. Check input data.'
     returntrue=.TRUE.
     v2=.TRUE.
  ELSEIF(errh==70) THEN
     text1='Specify profile values between 1 (night) and 2 (day).'
     v8=.TRUE.
  ENDIF
  !---------------------------------------------------------------------
  !This part of the code determines how the error message is written out

  IF(v1) THEN ! 1 real
     WRITE(500,*)'Error value: =', VALUE
  ELSEIF(v2) THEN ! 2 real
     WRITE(500,*)'Error values: =', VALUE, value2
  ELSEIF(v3) THEN ! 1 integer
     WRITE(500,*)'Error value: =', valueI
  ELSEIF(v4) THEN ! 2 real, 1 integer
     WRITE(500,*)'Error values: =', VALUE, value2, valueI
  ELSEIF(v5) THEN ! 1 real 1 integer
     WRITE(500,*)'Error values: =', VALUE, valueI
  ELSEIF(v6) THEN ! 2 integer
     valueI2=INT(VALUE)
     WRITE(500,*)'Error values: =', valueI, valueI2
  ELSEIF(v7) THEN
     valueI2=INT(value2)
     WRITE(500,*)'Error values: =', VALUE, valueI2, valueI
  ELSEIF(v8) THEN
     ! no error values
  ENDIF


  ! Write comment to problems.txt
  IF(returnTrue) THEN
     WRITE(500,*) TRIM(text1),' - WARNING'
  ELSE
     WRITE(500,*) 'ERROR! Program stopped: ',TRIM(text1)
  ENDIF
  !!Write the actual comment the problems file
  !write(500,*) trim(text1)

  CLOSE(500)

  !When returnTrue=true, then the program can continue despite the errors seen
  IF(returnTrue) THEN
     !write(*,*)'Problems.txt has been closed and overwritten if other errors occur'
     RETURN  !Continue program
  ENDIF

  CALL PauseStop(ProblemFile)        !Stop the program

  RETURN
END SUBROUTINE ErrorHint

!=============================================================

 SUBROUTINE ProblemsText(ProblemFile)

    USE defaultNotUsed
    IMPLICIT NONE

    CHARACTER (len=*):: ProblemFile

    !Opening problems.txt file: First option is selected if the file is opened for the first time
    !Second option for later points
    IF (errorChoice==0) THEN
        OPEN(500,file='problems.txt')
        WRITE(*,*) 'See problems.txt for possible issues in the run.'
        errorChoice=1
    ELSE
        OPEN(500,file='problems.txt',position="append")
    ENDIF

    !Writing of the problem file
    WRITE(500,*)'problem file: ',TRIM(ProblemFile)

    RETURN
 END SUBROUTINE ProblemsText


 SUBROUTINE PauseStop(ProblemFile)

   IMPLICIT NONE
   CHARACTER (len=*):: ProblemFile

   WRITE(*,*)'problem file: ',TRIM(ProblemFile)
   WRITE(*,*)'see problems.txt'

   !pause
   STOP
 END SUBROUTINE PauseStop
