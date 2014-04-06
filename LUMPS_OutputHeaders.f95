subroutine OutputHeaders(ProgName,lfnOut,lfnOutC,text,veg_type,ldown_option,selectHeader,errFileYes)
  
CHARACTER (len=90),DIMENSION(14)::text
integer :: Veg_type,lfnOut,lfnOutC, selectHeader,errFileYes,ldown_option
character(len=90)::ProgName
 
101 format("% common",a30,a30,/,'%',a30,a30)


if(selectHeader==2)then

    if (errFileYes==1) then
		!Modify these for Version 6.1!!(LJ)
    	write(lfnOut,*)"% Version=", trim( ProgName)
		write(lfnOut,101)trim(text(1)),trim(text(2)),trim(text(3)),trim(text(4))
    	write(lfnOut,*) '% veg_type:', veg_type,' NetRadiation/ldown_option:', ldown_option
    	write(lfnOut,30)

    endif
        
    30 format('%id it  dectime    kdown     kup    ldown     lup     Tsurf     qn      h_mod    e_mod',&
            '     qs        QF       QH       QE      P/i      Ie/i     E/i      DR/i    Ch/i', &
            '     ST/i    ROsoil/i  RO/i    ROpipe   ROpav     ROveg   ROwater    RA      RS     ustar',&
            '    L_mod SoilSt_pav SoilSt_blg SoilSt_everg SoilSt_dec SoilSt_Irrgr SoilSt_Gr St_pav ',&
            'St_blg St_everg St_dec St_Irrgr St_Gr St_water  Fcld    SoilState    smd       LAI',&
            '       Fw     addWater Iegrass/i Ietrees/i qn1_SF    qn1_S     Qm ',&
            '   delta_QS    Qrain     SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i ')
            
    !Clean data files
    write(lfnOutC,'(a)', advance = 'no') "% Version=" !formatting beginning of the row removes extra white 
    write(lfnOutC,*) trim( ProgName)  !here one can still use automatic formatting
	write(lfnOutC,101)trim(text(1)),trim(text(2)),trim(text(3)),trim(text(4))
    write(lfnOutC,'(a)', advance = 'no') '% veg_type:'
    write(lfnOutC,*) veg_type,' ldown_option:', ldown_option
    write(lfnOutC,30)
              
endif

!For daily and monthly files
if (selectHeader==1)then
	write(14,140)
	140 format('%day counter   qn      qs       qf     qe_S       pp      ext_Ie',&
        '     int_Ie   tot_ie      E_S     Change     R_Soil      R         Fw ',&
        '    addWater   QH_S      Qm   delta_QSI   Qrain    SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i')
	write(15,140)
endif

!LUMPS only header
if (selectHeader==3)then

  if (errFileYes==1) then
    write(lfnOut,*)"% Version=", trim( ProgName)
	write(lfnOut,101)trim(text(1)),trim(text(2)),trim(text(3))
    write(lfnOut,*) '% veg_type:', veg_type,' ldown_option:', ldown_option

    write(lfnout,40)
  endif
    
    40 format('%id   it  dectime    kdown    kup     ldown     lup     Tsurf        qn       h_mod       e_mod',&
              '     qs           qf        FCLD        V')
	
    write(lfnOutC,*)"% Version=", trim( ProgName)
	write(lfnOutC,101)trim(text(1)),trim(text(2)),trim(text(3))
    write(lfnOutC,*) '% veg_type:', veg_type,' ldown_option:', ldown_option

    write(lfnoutC,40)
 endif

 return
 end subroutine OutputHeaders
