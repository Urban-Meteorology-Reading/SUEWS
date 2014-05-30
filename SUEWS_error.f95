subroutine ErrorHint(errh,ProblemFile,value,value2,valueI) ! real
!errh        -- Create a numbered code for the situation so get a unique message to help solve the problem
!ProblemFile -- Filename where the problem occurs
!value       -- Real number with correct type
!value2      -- Real number with correct type
!valueI      -- Integer 2
!Last modified LJ 8 Feb 2013
!-------------------------------------------------------------------------------------------------

Use defaultNotUsed
real(kind(1d0)):: value,value2

character (len=*)::  ProblemFile
character (len=150)::  text1='unknown problem'
integer:: errh,ValueI,ValueI2
logical:: v1=.false.,v2=.false.,v3=.false.,v4=.false.,v5=.false.,v6=.false.,returnTrue=.false.,v7=.true.
  
	 call ProblemsText(ProblemFile) 
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
         text1='for z0m - GIS Input selected for z0_method - this value looks inappropriate'
          v1=.true.
     elseif(errh==6)then     	 
         text1='for zd - GIS Input selected for z0_method - this value looks inappropriate'
          v1=.true.
     elseif(errh==7) then
		 text1='HourWat should add to 1 in this file'
         v1=.true.
     elseif(errh==8) then
     	 text1='this should be zero - water Distribution'
         v1=.true.
	 elseif(errh==9) then
        text1= ' printed - one of these should be zero - - water Distribution'
        v2=.true.
 	 elseif(errh==10) then
        text1=' should add to 1  in above file'
        v2=.true.
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
       ! probably will stop but need to check this is really the problem
       returnTrue=.true.
     elseif(errh==19)then
      	text1=' cp, press , lv - not stopping after pause'
        v4=.true.
        returnTrue=.true.
     elseif(errh==20)then
     	text1=' skip lines, ios_out If running multiple grids, check the order of the grids.'
	    v5=.true.            
     elseif(errh==21)then
     	text1=' missing times in GIS file: id- GIS it-gis, it-met'
	    v4=.true.
      elseif(errh==22)then
     	text1=' QH_observed, QE_observed, QH_choice: '
	    v4=.true.         
	  elseif(errh==23)then
     	text1='CBL-sonde -need to increase size of izm:zmax,izm'
	    v5=.true. 
      elseif(errh==24) then
      	text1='CBL file problem - opening'
        v7=.true.
      elseif(errh==25) then
      	text1='CBL file problem -- reading sonde data, line:'
        v3=.true.  
      elseif(errh==26) then
        text1='Check that FileOutputPath and FileCode are specified and have double quotes around them'
        v7=.true.
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
        v4=.true.  ! 2 real, 1 integer 
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
        returntrue=.true.
        !Reserved for CBL
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
        text1 = 'Tstep and INTERVAL do not match!'
        v2 = .true.  !2 real    
      endif
    
     if(v1) then ! 1 real
        	write(500,*)'ERROR value: =', value
  		   write(*,*)'ERROR value: =', value
     elseif(v2) then ! 2 real
        	write(500,*)'ERROR values: =', value, value2
  		    write(*,*)'ERROR values: =', value, value2
     elseif(v3) then ! 1 integer
       	    write(500,*)'ERROR value: =', valueI
  	   	    write(*,*)'ERROR value: =', valueI
     elseif(v4) then ! 2 real, 1 integer 
            write(500,*)'ERROR values: =', value, value2, valueI
  		    write(*,*)'ERROR values: =', value, value2,valueI
     elseif(v5) then ! 1 real 1 integer
     	   write(500,*)'ERROR values: =', value, valueI
  		   write(*,*)'ERROR values: =', value, valueI
     elseif(v6) then ! 2 integer
      	   valueI2=int(value)     
     	   write(500,*)'ERROR values: =', valueI, valueI2
  		   write(*,*)'ERROR values: =', valueI, valueI2
    elseif(v7) then
    ! no error values  
    endif
    
    write(500,*)trim(text1)
    write(*,*)trim(text1)
    if(returnTrue)return  !Continue program 
    call PauseStop        !Stop the program
    
return
end subroutine ErrorHint



subroutine ProblemsText(ProblemFile)
character (len=*):: ProblemFile
	open(500,file='problems.txt')
    write(500,*)'problem file: ',trim(ProblemFile)
    write(*,*)'problem file: ',trim(ProblemFile)
    write(*,*)'see problems.txt'
    return
end subroutine ProblemsText

subroutine PauseStop
pause
stop
end subroutine PauseStop
