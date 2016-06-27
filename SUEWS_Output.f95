!In this subroutine the output files will be opened and the output matrices will be printed out.
!
!Last change:
! HCW 25 May 2016 - changed QF, QH QE in output file to qf, qh, qe to match input file. Also e.g. RA, RS, L_mod -> ra, rs, L_Ob
! HCW 12 Nov 2015
! Added z0m and zdm to output file
! HCW 27 Apr 2015
! Increased output resolution of water balance components (N.B. model time steps < 5 min may require higher precision)
! LJ in 13 April 2014
! FL in 10 June 2014
! HCW 18 Nov 2014
! LJ 5 Jan 2015: code cleaned, daily and monthly filesaving added
!-----------------------------------------------------------------------------------------------
SUBROUTINE SUEWS_Output(Gridiv, year_int, iv, irMax,  GridID)
  !INPUT: Gridiv = Grid number
  !       year_int = Year as a integer
  !       iv = Block number of met data
  !       irMax = Maximum number of rows in met data

  USE sues_data
  USE data_in
  USE allocateArray
  USE gis_data
  USE time
  USE defaultNotUsed
  USE initial
  USE solweig_module
  USE cbl_module
  USE ESTM_data


  IMPLICIT NONE

  INTEGER::i ,azero!,lfnOutC
  INTEGER:: Gridiv, year_int, iv, irMax, GridID
  CHARACTER(len=10):: str2, grstr2, yrstr2
  CHARACTER(len=100):: rawpath, SnowOut,ESTMOut

  !================DEFINE OUTPUT FILENAME AND ITS PATH================
  WRITE(str2,'(i2)') TSTEP/60
  WRITE(grstr2,'(i5)') GridID
  WRITE(yrstr2,'(i4)') year_int

  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grstr2))//'_'//TRIM(ADJUSTL(yrstr2))
  FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(str2))//'.txt'
  SOLWEIGpoiOut=TRIM(rawpath)//'_SOLWEIGpoiOut.txt'
  ESTMOut=TRIM(rawpath)//'_ESTM_5.txt'
  BLOut=TRIM(rawpath)//'_BL.txt'
  SnowOut=TRIM(rawpath)//'_snow_5.txt'

  !================OPEN OUTPUT FILE AND PRINT HEADER================

  ! Hourly output file
  lfnOutC=39  !Output file code

  IF(iv == 1) THEN
     OPEN(lfnOutC,file=TRIM(FileOut),err=112)
     WRITE(lfnOutC,110)

110  FORMAT('%iy id it imin dectime ',&
     'kdown kup ldown lup Tsurf qn h_mod e_mod qs qf qh qe ',&
     'P/i Ie/i E/i Dr/i ',&
     'St/i NWSt/i surfCh/i totCh/i ',&
     'St/i NWSt/i surfCh/i totCh/i ',&
     'St/i NWSt/i surfCh/i totCh/i ',&
     'RO/i ROsoil/i ROpipe ROpav ROveg ROwater ',&
     'AdditionalWater FlowChange WU_int WU_EveTr WU_DecTr WU_Grass ',&
     'ra rs ustar L_Ob Fcld ',&
     'SoilSt smd smd_Paved smd_Bldgs smd_EveTr smd_DecTr smd_Grass smd_BSoil ',&
     'St_Paved St_Bldgs St_EveTr St_DecTr St_Grass St_BSoil St_Water ',&
     'LAI z0m zdm ',&
     'qn1_SF qn1_S Qm QmFreez QmRain SWE Mw MwStore snowRem_Paved snowRem_Bldgs ChSnow/i ',&
     'SnowAlb ')

     !These belon to NARP ouput file
     ! 'kup_Paved kup_Bldgs kup_EveTr kup_DecTr kup_Grass kup_BSoil kup_Water ',&
     ! 'lup_Paved lup_Bldgs lup_EveTr lup_DecTr lup_Grass lup_BSoil lup_Water ',&
     ! 'Ts_Paved Ts_Bldgs Ts_EveTr Ts_DecTr Ts_Grass Ts_BSoil Ts_Water ',&
     ! 'qn_Paved qn_Bldgs qn_EveTr qn_DecTr qn_Grass qn_BSoil qn_Water ',&

  ELSE
     OPEN(lfnOutC,file=TRIM(FileOut),position='append')!,err=112)
  ENDIF

  !SOLWEIG outputfile
  IF (SOLWEIGpoi_out==1) THEN
     OPEN(9,file=SOLWEIGpoiOut)
     WRITE(9,113)
113  FORMAT('%doy dectime  azimuth altitude GlobalRad DiffuseRad DirectRad ',&
          ' Kdown2d    Kup2d    Ksouth     Kwest    Knorth     Keast ',&
          ' Ldown2d    Lup2d    Lsouth     Lwest    Lnorth     Least ',&
          '   Tmrt       I0       CI        gvf      shadow    svf    svfbuveg    Ta    Tg')
  ENDIF

  !BL ouputfile
  IF (CBLuse>=1)THEN
     OPEN(53,file=BLOut,status='unknown')
     WRITE(53, 102)
102  FORMAT('iy  id   it imin dectime         z            theta          q',&
          '               theta+          q+              Temp_C          rh',&
          '              QH_use          QE_use          Press_hPa       avu1',&
          '            ustar           avdens          lv_J_kg         avcp',&
          '            gamt            gamq')
  ENDIF

  !Snow outputfile
  IF (SnowUse>=1) THEN
     IF(iv == 1) THEN
        OPEN(54,file=SnowOut)
        WRITE(54, 114)
114     FORMAT('%iy  id   it imin dectime ',&
             'SWE_Paved SWE_Bldgs SWE_EveTr SWE_DecTr SWE_Grass SWE_BSoil SWE_Water ',&
             'Mw_Paved Mw_Bldgs Mw_EveTr Mw_DecTr Mw_Grass Mw_BSoil Mw_Water ',&
             'Qm_Paved Qm_Bldgs Qm_EveTr Qm_DecTr Qm_Grass Qm_BSoil Qm_Water ',&
             'Qa_Paved Qa_Bldgs Qa_EveTr Qa_DecTr Qa_Grass Qa_BSoil Qa_Water ',&
             'QmFr_Paved QmFr_Bldgs QmFr_EveTr QmFr_DecTr QmFr_Grass QmFr_BSoil QmFr_Water ',&
             'fr_Paved fr_Bldgs fr_EveTr fr_DecTr fr_Grass fr_BSoil ',&
             'RainSn_Paved RainSn_Bldgs RainSn_EveTr RainSn_DecTr RainSn_Grass RainSn_BSoil RainSn_Water ',&
             'Qn_PavedSnow Qn_BldgsSnow Qn_EveTrSnpw Qn_DecTrSnow Qn_GrassSnpw Qn_BSoilSnow Qn_WaterSnow ',&
             'kup_PavedSnow kup_BldgsSnow kup_EveTrSnpw kup_DecTrSnow kup_GrassSnpw kup_BSoilSnow kup_WaterSnow ',&
             'frMelt_Paved frMelt_Bldgs frMelt_EveTr frMelt_DecTr frMelt_Grass frMelt_BSoil frMelt_Water ',&
             'MwStore_Paved MwStore_Bldgs MwStore_EveTr MwStore_DecTr MwStore_Grass MwStore_BSoil MwStore_Water ',&
             'SnowDens_Paved SnowDens_Bldgs SnowDens_EveTr SnowDens_DecTr SnowDens_Grass SnowDens_BSoil SnowDens_Water ',&
             'Sd_Paved Sd_Bldgs Sd_EveTr Sd_DecTr Sd_Grass Sd_BSoil Sd_Water ',&
             'Tsnow_Paved Tsnow_Bldgs Tsnow_EveTr Tsnow_DecTr Tsnow_Grass Tsnow_BSoil Tsnow_Water')
     ELSE
        OPEN(54,file=TRIM(SnowOut),position='append')
     ENDIF
  ENDIF

  !ESTM ouputfile
  IF (QSChoice==4 .OR. QSChoice==14)THEN
     OPEN(58,file=ESTMOut,status='unknown')
     WRITE(58, 115)
115  FORMAT('%iy id it imin dectime ',&
          'QSNET QSAIR QSWALL QSROOF QSGROUND QSIBLD ',&
          'TWALL1 TWALL2 TWALL3 TWALL4 TWALL5 ',&       !T0_WALL TWALL1 TWALL2 TWALL3 TN_WALL
          'TROOF1 TROOF2 TROOF3 TROOF4 TROOF5 ',&
          'TGROUND1 TGROUND2 TGROUND3 TGROUND4 TGROUND5 ',&
          'TiBLD1 TiBLD2 TiBLD3 TiBLD4 TiBLD5 TaBLD ')
  ENDIF

  !================ACTUAL DATA WRITING================
  IF (SOLWEIGpoi_out==1) THEN
     DO i=1,SolweigCount-1
        WRITE(9,304) INT(dataOutSOL(i,1,Gridiv)),(dataOutSOL(i,is,Gridiv),is = 2,28)
     ENDDO
  ENDIF

  DO i=1,irMax
     WRITE(lfnoutC,301) INT(dataOut(i,1,Gridiv)),INT(dataOut(i,2,Gridiv)),INT(dataOut(i,3,Gridiv)),INT(dataOut(i,4,Gridiv)),&
          dataOut(i,5:ncolumnsDataOut,Gridiv)
  ENDDO

  IF(CBLuse>=1) THEN
     DO i=1,iCBLcount
        WRITE(53,305)(INT(dataOutBL(i,is,Gridiv)),is=1,4),(dataOutBL(i,is,Gridiv),is=5,22)
     ENDDO
  ENDIF

  IF(SnowUse>=1) THEN
     DO i=1,irmax
        WRITE(54,306)(INT(dataOutSnow(i,is,Gridiv)),is=1,4),(dataOutSnow(i,is,Gridiv),is=5,ncolumnsDataOutSnow)
     ENDDO
  ENDIF
  IF (QSChoice==4 .OR. QSChoice==14)THEN
     DO i=1,iESTMcount
        WRITE(58, 307)(INT(dataOutESTM(i,is,Gridiv)),is=1,4),(dataOutESTM(i,is,Gridiv),is=5,32)
     ENDDO
  ENDIF

  !================WRITING FORMAT================
  ! Main output file at model timestep
  ! Do NOT change from 301 here - read by python wrapper
  ! 301_Format
301 FORMAT((i4,1X),3(i3,1X),(f8.4,1X),&       !time parameters 5
       5(f9.4,1X),7(f9.4,1X),&            !17
       4(f10.6,1X),&                      !21
       1(f10.5,1X),3(f10.6,1X),&          !25
       6(f10.6,1X),&                      !31
       2(f9.3,1X),4(f9.4,1X),&            !37
       3(f10.5,1X),(g14.7,1X),(f10.5,1X),& !42
       2(f10.4,1X),6(f10.5,1X),7(f10.4,1X),& !57
       3(f10.4,1X),&                       !60 LAI z0m zdm
       5(f10.4,1X),6(f10.6,1X),&           !71
       1(f8.4,1X))                        !72 albedo snow

  !==================== This part read by python wrapper ======================
  ! Update to match output columns, header and format
  ! Average, sum, or use last value to go from model timestep to 60-min output
  ! 301_Instructions
  ! TimeCol = [1,2,3,4,5]
  ! AvCol  = [6,7,8,9,10,11,12,13,14,15,16,17,  38,39,40,41,42,  61,62,63,64,65]
  ! SumCol = [18,19,20,21,  24,25, 26,27,28,29,30,31,  32,33,34,35,36,37,  69,70,71]
  ! LastCol  = [22,23,  43,44,45,46,47,48,49,50,  51,52,53,54,55,56,57,  58,59,60,  66,67,68,  72]

304 FORMAT(1(i3,1X),4(f8.4,1X),23(f9.3,1X))          !Solweig output
305 FORMAT((i4,1X),3(i3,1X),(f8.4,1X),17(f15.7,1x))  !CBL output
306 FORMAT((i4,1X),3(i3,1X),(f8.4,1X)&               !Snow out
       7(f10.4,1X),7(f10.4,1X),7(f10.4,1X),&
       7(f10.4,1X),7(f12.4,1X),6(f10.4,1X),&
       7(f10.4,1X),&
       7(f10.4,1X),7(f10.4,1X),7(f10.4,1X),&
       7(f10.4,1X),7(f10.4,1X),7(f10.4,1X),7(f10.4,1X))
307 FORMAT((i4,1X),3(i3,1X),(f8.4,1X),27(f12.4,1X))  !ESTM out

  !================CLOSE OUTPUTFILE================
  CLOSE (lfnoutC)
  CLOSE (9)
  CLOSE (53)
  CLOSE (54)
  CLOSE (58)
  RETURN

  !Error commands
112 CALL ErrorHint(52,TRIM(fileOut),notUsed,notUsed,notUsedI)


END SUBROUTINE SUEWS_Output
