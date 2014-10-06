!In this subroutine the output files will be opened and the output matrixes will be printed out.
!Currently works only with hourly, NARP and snow output files (for singe grid)
!Last change by LJ in 13 April 2014
!LAst change by FL in 10 June 2014
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
    use defaultNotUsed


    IMPLICIT NONE

    CHARACTER (len=90),DIMENSION(14)::text
    integer::i,j !,lfnOutC

    !First open output files and print headers to them
    if (NARPOutput==1) then
        open(7,file=NARPOut)
        write(7,110)
        110 format('%id  dectime   kup_pav   kup_blgs   kup_ET    kup_dT    kup_IG    kup_UG   kup_wtr  ', &
                   ' lup_pav  lup_bldg    lup_ET    lup_DT   lup_IG    lup_UG    lup_wtr    Ts_pav    Ts_bldg    Ts_ET',&
                   '     Ts_DT     Ts_IG     Ts_UG    Ts_wtr    qn_pav   qn_bldg     qn_ET     qn_DT     qn_IG    qn_UG     qn_wtr')
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
        111  format('%id   it dectime  SWE_pav SWE_bldg   SWE_ET   SWE_DT   SWE_IG   SWE_UG  SWE_wtr ',&
                'SnowRem_pav SnowRem_bldg Mw  Mw_pav  Mw_bldg    Mw_ET    Mw_DT    Mw_IG    Mw_UG ',&
                '  Mw_wtr     Qm   Qm_pav Qm_bldg   Qm_ET   Qm_DT   Qm_IG   Qm_UG  Qm_wtr ',&
                ' Qa_pav Qa_bldg   Qa_ET   Qa_DT   Qa_IG   Qa_UG  Qa_wtr QmFr_pav ',&
                'QmFr_bldg QmFr_ET QmFr_DT QmFr_IG QmFr_UG QmFr_wtr ',&
                'fr_pav fr_bldg fr_ET  fr_DT   fr_IG   fr_UG alb_snow RainSn_pav ',&
                'RainSn_bldg RainSn_ET RainSn_DT RainSn_IG RainSn_UG ',&
                'RainSn_wtr Qn_pavSnow Qn_blgsSnow Qn_ETSnow Qn_DTSnow Qn_IGSnow Qn_UG ',&
                'Qs_wtrSnow kup_pavSnow kup_blgsSnow kup_ETSnow kup_DTSnow kup_IGSnow kup_UGSnow ',&
                'kup_wtrSnow frMelt_pav frMelt_bldg frMelt_ET frMelt_DT frMelt_IG frMelt_UG frMelt_wtr',&
                'MwStore_pav MwStore_bldg MwStore_ET MwStore_DT MwStore_IG MwStore_UG MwStore_wtr',&
                'densSnow_pav densSnow_bldg densSnow_ET densSnow_DT densSnow_IG densSnow_UG densSnow_wtr',&
                'Sd_pav Sd_bldg Sd_ET Sd_DT Sd_IG Sd_UG Sd_water Tsnow_pav Tsnow_bldg Tsnow_ET Tsnow_DT',&
                'Tsnow_IG Tsnow_UG Tsnow_wtr delta_Qi')
    endif

    !SOLWEIG outputfile
    if (SOLWEIGpoi_out==1) then
        open(9,file=SOLWEIGpoiOut)
        write(9,113)
113     format('%doy dectime  azimuth altitude GlobalRad DiffuseRad DirectRad ',&
                ' Kdown2d    Kup2d    Ksouth     Kwest    Knorth     Keast ',&
                ' Ldown2d    Lup2d    Lsouth     Lwest    Lnorth     Least ',&
                '   Tmrt       I0       CI        gvf      shadow    svf    svfbuveg    Ta    Tg')
    endif    
    
    if(write5min==1) then                ! if going to write 5 min data out
       open(16,file=file5min,err=204)
       write(16,161)

 161   	format('%id 5min dectime    pp       Ie        E     St_pav   St_bldg   St_ET    St_DT ' ,&
               '   St_IG    St_UG St_wtr SoilSt_pav SoilSt_bldg SoilSt_ET SoilSt_DT',&
               ' SoilSt_IG SoilSt_UG D_pav D_bldg  D_ET    D_DT    D_IG     D_UG     r_pav',&
               '     r_bldg    r_ET    r_DT    r_IG    r_UG  soilr_pav soilr_bldg soilr_ET soilr_DT',&
               ' soilr_IG soilr_UG snowr_pav snowr_bldg snowr_ET snowr_DT snowr_IG snowr_UG',&
               ' SWE_pav SWE_bldg SWE_ET SWE_DT   SWE_IG   SWE_UG  SWE_wtr snowCh_pav snowCh_bldg',&
               ' snowCh_ET snowCh_DT snowCh_IG snowCh_UG snowCh_wtr mwh_pav mwh_bldg mwh_ET',&
               ' mwh_DT mwh_IG mwh_UG mwh_wtr')
    endif
    
     !BL ouputfile
     if (CBLuse>=1)then
        open(53,file=BLOut,status='unknown')
		write(53, 102)
102  	format('id    it   dectime    z         theta       q       theta+      q+       Temp_C',&
               '    rh        QH_use    QE_use  Press_hPa    avu1      ustar     avdens',&
               '    lv_J_kg avcp         gamt       gamq') 
     endif

    !Actual data writings
    do i=1,nlines
        write(lfnoutC,301) int(dataOut1(i,1)),int(dataOut1(i,2)),(dataOut1(i,is),is = 3,62)
        if (NARPOutput==1)  write(7,117) int(dataOut2(i,1)),(dataOut2(i,is),is = 2,30)
        if (snowUse==1)  write(8,118) int(dataOut3(i,1)),int(dataOut1(i,2)),(dataOut3(i,is),is = 3,106)
        if (SOLWEIGpoi_out==1) write(9,119) int(dataOutSOL(i,1)),(dataOutSOL(i,is),is = 2,28) 
    enddo

    if(write5min==1) then
      do i=1,nlines*int(INTERVAL/Tstep)
         write(16,36) int(dataOut5min(i,1)),int(dataOut5min(i,2)),(dataOut5min(i,is),is = 3,64)

      enddo
    endif

	if(CBLuse>=1) then
    	do i=1,(nlines-1)*nCBLstep
        write(*,*)i,int(dataOutBL(i,1)),dataOutBL(i,2),dataOutBL(i,3)
    		write(53,114)int(dataOutBL(i,1)),(dataOutBL(i,is),is=2,20) 
    	enddo  
    endif  

!Writing formats
!301 format(2i4,f9.4,5f8.2,7f12.4,12f10.4,2f9.1,f10.4,g15.5,27f12.4)!s.o. 7F11.4 ->9F11.4
!New format Nov 2013 LJ
301 format(2(i3,1X),f8.4,4(f8.2,1X),(f7.2,1X),7(f8.2,1X),12(f8.3,1X),2(f7.1,1X),&
           (f7.2,1X),(g14.5,1X),16(f8.3,1X),(f9.2,1X),4(f9.3,1X),5(f8.2,1X),5(f9.3,1X))!
118 format(2(i3,1X),(f8.4,1X),17(f8.3,1X),22(f7.2,1X),7(f7.2,1X),7(f7.3,1X),14(f7.2,1X),14(f7.3,1X),21(f7.2,1X),(f14.3,1X))
119 format(1(i3,1X),4(f8.4,1X),23(f9.3,1X))
117 format(i3,f9.4,28f10.3)
36  format(2i3,f9.4,61f9.3)     !format(i3,i3,3f12.4,6f10.4, f14.3,25f10.4)
114 format((i3,1X),(f5.2,1X),(f8.4,1X),14(f15.7,1x),3(f15.7,1X))

!Close all output files
 close (lfnoutC)
 close (8)
 close (7)
 close (9)

 return

 !Error commands
 112 call ErrorHint(52,trim(fileOut),notUsed,notUsed,notUsedI)
 204 call ErrorHint(52,trim(file5min),notUsed,notUsed,notUsedI)


 end subroutine