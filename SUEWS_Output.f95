!In this subroutine the output files will be opened and the output matrices will be printed out.
!
!Last change:
! LJ in 13 April 2014
! FL in 10 June 2014
! HCW 18 Nov 2014
! LJ 5 Jan 2015: code cleaned, daily and monthly filesaving added
!-----------------------------------------------------------------------------------------------
 subroutine SUEWS_Output(Gridiv, year_int, iv, irMax) !(DataOut1,DataOut2,DataOut3,NARPOutput,snowUse,)
 !INPUT: Gridiv = Grid number
 !       year_int = Year as a integer
 !       iv = Block number of met data
 !       irMax = Maximum number of rows in met data

  use sues_data
  use data_in
  use allocateArray
  use gis_data
  use time
  use defaultNotUsed
  use initial

  IMPLICIT NONE

  character (len=90),dimension(14)::text
  integer::i,j !,lfnOutC
  integer:: Gridiv, year_int, iv, irMax
  character(len=10):: str2, grstr2, yrstr2
  character(len=100):: rawpath    
    
  !================DEFINE OUTPUT FILENAME AND ITS PATH================
  write(str2,'(i2)') TSTEP/60
  write(grstr2,'(i2)') Gridiv
  write(yrstr2,'(i4)') year_int

  rawpath=trim(FileOutputPath)//trim(FileCode)//trim(adjustl(grstr2))//'_'//trim(adjustl(yrstr2))
  FileOut=trim(rawpath)//'_'//trim(adjustl(str2))//'.txt'

  !================OPEN OUTPUT FILE AND PRINT HEADER================

  ! Hourly output file
  lfnOutC=39  !Output file code

  if(iv == 1) then
    open(lfnOutC,file=trim(FileOut),err=112)
    write(lfnOutC,110)

    110 format('%iy id it imin  dectime    kdown     kup    ldown     lup     Tsurf     qn      h_mod    e_mod',&
        '     qs        QF       QH       QE      P/i      Ie/i     E/i      DR/i    Ch/i', &
        '     ST/i    ROsoil/i  RO/i    ROpipe   ROpav     ROveg   ROwater    RA      RS     ustar',&
        '    L_mod SoilSt_pav SoilSt_bldg SoilSt_ET SoilSt_DT SoilSt_IG SoilSt_UG St_pav ',&
        'St_bldg   St_ET    St_DT    St_IG    St_UG   St_wtr     Fcld SoilState    smd       LAI',&
        '       Fw     addWater Ie/i Ie/i qn1_SF    qn1_S     Qm ',&
        '   delta_QS    Qrain     SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i kup_pav   kup_blgs',&
        'kup_ET    kup_dT    kup_IG    kup_UG   kup_wtr  ', &
        ' lup_pav  lup_bldg    lup_ET    lup_DT   lup_IG    lup_UG    lup_wtr    Ts_pav    Ts_bldg    Ts_ET',&
        '     Ts_DT     Ts_IG     Ts_UG    Ts_wtr    qn_pav   qn_bldg     qn_ET     qn_DT     qn_IG    qn_UG     qn_wtr',&
        'SWE_pav SWE_bldg   SWE_ET   SWE_DT   SWE_IG   SWE_UG  SWE_wtr ',&
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
        'Tsnow_IG Tsnow_UG Tsnow_wtr')
  else
    open(lfnOutC,file=trim(FileOut),position='append')!,err=112)
  endif

 !================ACTUAL DATA WRITING================

  do i=1,irMax
        write(lfnoutC,301) int(dataOut(i,1,Gridiv)),int(dataOut(i,2,Gridiv)),int(dataOut(i,3,Gridiv)),&
                           int(dataOut(i,4,Gridiv)),(dataOut(i,is,Gridiv),is = 5,192)
  enddo

  !================WRITING FORMAT================
  301 format((i4,1X),3(i3,1X),(f8.4,1X),(f8.2,1X),3(f8.2,1X),(f7.2,1X),7(f8.2,1X),12(f8.3,1X),2(f7.1,1X),&
          (f7.2,1X),(g14.5,1X),16(f8.3,1X),(f9.2,1X),4(f9.3,1X),5(f8.2,1X),5(f9.3,1X),&
          28f10.3,17(f8.3,1X),22(f8.2,1X),7(f8.2,1X),7(f8.3,1X),14(f8.2,1X),&
          14(f8.3,1X),21(f8.2,1X),(f14.3,1X))

  !304 format(1(i3,1X),4(f8.4,1X),23(f9.3,1X))                                                   !Solweig output
  !305 format((i3,1X),(f5.2,1X),(f8.4,1X),14(f15.7,1x),3(f15.7,1X))                              !CBL output


  !================CLOSE OUTPUTFILE================
  close (lfnoutC)

  return

  !Error commands
  112 call ErrorHint(52,trim(fileOut),notUsed,notUsed,notUsedI)
  204 call ErrorHint(52,trim(file5min),notUsed,notUsed,notUsedI)


 end subroutine