!In this subroutine the output files will be opened and the output matrixes will be printed out.
!Currently works only with hourly, NARP and snow output files (for singe grid)
!Last change by LJ in 13 April 2014
!-----------------------------------------------------------------------------------------------
 subroutine SUEWS_Output!(DataOut1,DataOut2,DataOut3,NARPOutput,snowUse,)
 !INPUT: DataOut1 = Main data output matrix
 !       DataOut2 = NARP output matrix
 !       DataOut3 = Snow output matrix

	use sues_data
    use data_in
    use allocateArray
    use gis_data
    use time

    IMPLICIT NONE

    CHARACTER (len=90),DIMENSION(14)::text
    integer::i,j !,lfnOutC

    !First open output files and print headers to them
    if (NARPOutput==1) then
        open(7,file=NARPOut)
        write(7,110)
        110 format('%id  dectime   kup_pav   kup_blg kup_everg   kup_dec kup_Irrgr    kup_Gr   kup_wtr  ', &
                   ' lup_pav   lup_blg  lup_everg  lup_dec lup_Irrgr   lup_Gr    lup_wtr    Ts_pav    Ts_blg   Ts_everg',&
                   '  Ts_dec   Ts_Irrgr     Ts_Gr    Ts_wtr    qn_pav   qn_blg   qn_everg   qn_dec   qn_Irrgr    qn_Gr   qn_wtr')
    endif

    !Actual output file
    lfnOutC=39 !Clean output file
    open(lfnOutC,file=trim(FileOut),err=112)

    !Header output is printedll
    text(1)="_FileChoices.txt: options selected"
    text(2)="DailyState.txt:LAI,HDD etc"
    text(3)="% "
    text(4)="% "
    call OutputHeaders(ProgName,lfnOutC,text,veg_type,ldown_option,2) ! LUMPS_OutputHeaders.f95
    do j=1,4
      keepHeader(j)=text(j)
    end do

    !Snow outputfile
    if (snowUse==1) then
        open(8,file=SnowOut)
        write(8,111)
        111  format('%doy  it dectime  SWE_pav SWE_bldg SWE_evergr SWE_dec SWE_irrGr SWE_Gr SWE_water ',&
                ' SnowRem_pav SnowRem_bldg Mw Mw_pav Mw_bldg Mw_evergr  Mw_dec Mw_irrGr    Mw_Gr ',&
                'Mw_water     Qm   Qm_pav Qm_bldg Qm_evergr Qm_dec Qm_irrGr Qm_Gr Qm_water ',&
                'Qa_pav Qa_bldg Qa_evergr Qa_dec Qa_irrGr Qa_Gr Qa_water QmFr_pav ',&
                'QmFr_bldg QmFr_evergr QmFr_dec QmFr_irrGr QmFr_Gr QmFr_water ',&
                'fr_pav fr_bldg fr_evergr fr_dec fr_irrGr fr_Gr alb_snow  rainOnSnow_pav ',&
                'rainOnSnow_bldg rainOnSnow_evergr rainOnSnow_dec rainOnSnow_irrGr rainOnSnow_Gr',&
                'rainOnSnow_water Qn_pavSnow Qs_blgSnow Qs_evergrSnow Qs_decSnow Qs_irrGrSnow Qs_GrSnow ',&
                'Qs_wtrSnow kup_pavSnow kup_blgSnow kup_evergrSnow kup_decSnow kup_irrGrSnow kup_GrSnow ',&
                'kup_wtrSnow ')

    endif

    if(write5min==1) then                ! if going to write 5 min data out
       open(16,file=file5min,err=204)
       write(16,161)

 161   	format('%id 5min dectime    pp       Ie        E     St_pav   St_blg   St_everg St_dec ' ,&
               ' St_IrrGr  St_Gr   St_water SoilSt_pav SoilSt_blg SoilSt_everg SoilSt_dec',&
               ' SoilSt_IrrGr SoilSt_Gr D_pav D_blg  D_everg    D_dec    D_IrrGr     D_Gr     r_pav',&
               '     r_blg    r_everg    r_dec    r_IrrGr    r_Gr  soilr_pav soilr_bldg soilr_everg soilr_dec',&
               ' soilr_IrrGr soilr_Gr snowr_pav snowr_bldg snowr_everg snowr_dec snowr_IrrGr snowr_Gr',&
               ' SWE_pav SWE_bldg SWE_everg SWE_dec SWE_IrrGr SWE_Gr SWE_water snowCh_pav snowCh_bldg',&
               ' snowCh_everg snowCh_dec snowCh_IrrGr snowCh_Gr snowCh_water mwh_pav mwh_bldg mwh_everg',&
               ' mwh_dec mwh_IrrGr mwh_Gr mwh_water')
     endif

    !Actual data writings
    do i=1,nlines
        write(lfnoutC,301) int(dataOut1(i,1)),int(dataOut1(i,2)),(dataOut1(i,is),is = 3,62)

        if (NARPOutput==1)  write(7,117) int(dataOut2(i,1)),(dataOut2(i,is),is = 2,30)
        if (snowUse==1)  write(8,118) int(dataOut3(i,1)),int(dataOut1(i,2)),(dataOut3(i,is),is = 3,106)
    enddo

    if(write5min==1) then
      do i=1,nlines*int(INTERVAL/Tstep)
         write(16,36) int(dataOut5min(i,1)),int(dataOut5min(i,2)),(dataOut5min(i,is),is = 3,69)

      enddo
    endif

!Writing formats
!301 format(2i4,f9.4,5f8.2,7f12.4,12f10.4,2f9.1,f10.4,g15.5,27f12.4)!s.o. 7F11.4 ->9F11.4
!New format Nov 2013 LJ
301 format(2(i3,1X),f8.4,4(f8.2,1X),(f7.2,1X),7(f8.2,1X),12(f8.3,1X),2(f7.1,1X),&
           (f7.2,1X),(g14.5,1X),16(f8.3,1X),(f9.2,1X),4(f9.3,1X),5(f8.2,1X),5(f9.3,1X))!
118 format(2(i3,1X),(f8.4,1X),17(f8.3,1X),22(f7.2,1X),7(f7.2,1X),7(f7.3,1X),14(f7.2,1X),14(f7.3,1X),21(f7.2,1X),(f14.3,1X))
117 format(i3,f9.4,28f10.3)
36  format(2i3,f9.4,66f9.3)     !format(i3,i3,3f12.4,6f10.4, f14.3,25f10.4)

 close (lfnoutC)
 close (8)
 close (7)

 return

!Error commands
112 call ProblemsText(trim(FileOut))
    call PauseStop

204 call ProblemsText(trim(file5min))
    call PauseStop

 end subroutine