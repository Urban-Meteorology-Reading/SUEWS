!In this subroutine the output files will be opened and the output matrices will be printed out.
!
!Last change:
! HCW 12 Dec 2016 - Restructured writing of output files and introduced families of output variables
! HCW 04 Jul 2016 - GridID can now be up to 10 digits long. If file not found or not read correctly, program stops
! HCW 29 Jun 2016 - Fixed bug in output file (4 columns were repeated twice)
!                 - Changed ESTM output file to be written for more than one met block
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
SUBROUTINE SUEWS_Output(Gridiv, year_int, iv, irMax, CurrentGrid)
  !INPUT: Gridiv   = Grid number
  !       year_int = Year as a integer
  !       iv       = Block number of met data
  !       irMax    = Maximum number of rows in met data
  !       CurrentGrid   = Grid ID (according to SiteSelect.txt)

  USE allocateArray
  USE cbl_module
  USE data_in
  USE defaultNotUsed
  USE ESTM_data
  USE gis_data
  USE initial
  USE SetupOutput
  USE solweig_module
  USE sues_data 
  USE time

  IMPLICIT NONE

  INTEGER:: Gridiv, year_int, iv, irMax, CurrentGrid !inputs
  INTEGER:: i  
  
  CHARACTER(len=10):: str2, grstr2, yrstr2
  CHARACTER(len=100):: rawpath, SnowOut,ESTMOut, FileOutFormat

  CHARACTER(len=10),DIMENSION(nColumnsDataOut):: HeaderAll, FormatAll  !Header and formats for all output variables
  CHARACTER(len=12),DIMENSION(nColumnsDataOut):: UnitsAll              !Units and formats for all output variables
  CHARACTER(len=12*nColumnsDataOut):: HeaderOut, FormatOut, UnitsOut   !Header and format for selected output variables (untrimmed)
  CHARACTER(len=1),DIMENSION(nColumnsDataOut):: AggAll                !Aggregation method required
  CHARACTER(len=5*nColumnsDataOut),ALLOCATABLE:: AggOut 
  CHARACTER(len=3*nColumnsDataOut),ALLOCATABLE:: ColNos
  CHARACTER(len=10):: fy, ft, fd, f94, f104, f106   !Useful formats
  CHARACTER(len= 1):: aT, aA, aS, aL   !Useful formats
  CHARACTER(len= 5):: itext

  ! Define useful formats here  
  fy   = '(i0004,1X)'   !4 digit integer for year
  ft   = '(i0003,1X)'   !3 digit integer for id, it, imin
  fd   = '(f08.4,1X)'   !3 digits + 4 dp for dectime
  f94  = '(f09.4,1X)'   !standard output format: 4 dp + 4 digits
  f104 = '(f10.4,1X)'   !standard output format: 4 dp + 5 digits
  f106 = '(f10.6,1X)'   !standard output format: 6 dp + 3 digits
  ! Define aggregation methods here (for wrapper)
  aT = '0'   !time columns
  aA = '1'   !average
  aS = '2'   !sum
  aL = '3'   !last value
  
  !========== Set file path and file names ==========
  WRITE(str2,'(i2)') TSTEP/60
  WRITE(grstr2,'(i10)') CurrentGrid
  WRITE(yrstr2,'(i4)') year_int

  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grstr2))//'_'//TRIM(ADJUSTL(yrstr2))
  FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(str2))//'.txt'
  SOLWEIGpoiOut=TRIM(rawpath)//'_SOLWEIGpoiOut.txt'
  ESTMOut=TRIM(rawpath)//'_ESTM_5.txt'
  BLOut=TRIM(rawpath)//'_BL.txt'
  SnowOut=TRIM(rawpath)//'_snow_5.txt'
  FileOutFormat=TRIM(FileOutputPath)//TRIM(FileCode)//'_YYYY_'//TRIM(ADJUSTL(str2))//'_OutputFormat.txt'
  
  
  !========== Get headers and write out output info to file ==========
  ! To add extra columns, change all these (Header, Units, Format, Agg) together
  ! Could change to read from external file later
  IF(OutputFormats==1) THEN   !Once per run
  
     ! Set all output variables here. This must agree with dataOut (see SUEWS_Calculations.f95) and FormatAll
     HeaderAll(:) = (/ '      Year','       DOY','      Hour','       Min','   Dectime', &   !datetime info (1-5)
                       '     Kdown','       Kup','     Ldown','       Lup','     Tsurf', &   !radiation components (6-10)
                       '        QN','        QF','        QS','        QH','        QE', &   !energy fluxes (11-15)
                       '   QHlumps','   QElumps','   QHresis', &                             !energy fluxes (other approaches) (16-18)
                       '      Rain','       Irr','      Evap','        RO','     TotCh', &   !water balance components (19-23)
                       '    SurfCh','     State',' NWtrState','  Drainage','       SMD', &   !water balance components cont. (24-28)
                       '    FlowCh','  AddWater', &                                          !water balance components cont. (29-30)
                       '    ROSoil','    ROPipe','     ROImp','     ROVeg','   ROWater', &   !runoff components (31-35)
                       '     WUInt','   WUEveTr','   WUDecTr','   WUGrass', &                !water use (36-39)
                       '  SMDPaved','  SMDBldgs','  SMDEveTr','  SMDDecTr','  SMDGrass','  SMDBSoil', &   !smd for each surface (40-45)
                       '   StPaved','   StBldgs','   StEveTr','   StDecTr','   StGrass','   StBSoil','   StWater',&   !states (46-52)
                       '    Zenith','   Azimuth','   AlbBulk','      Fcld', &                ! extra radiation info (53-56)
                       '       LAI','       z0m','       zdm', &                             ! extra surface info (57-59)
                       '     ustar','       Lob','        ra','        rs', &                ! turbulence (60-63)
                       '        Fc', &                                                       ! CO2 flux (64)
                       '   FcPhoto','   FcRespi','   FcMetab','   FcTraff','   FcBuild', &   ! CO2 flux components (65-69)
                       '  QNSnowFr','    QNSnow','   AlbSnow', &                             ! snow-related (radiation) (70-72)
                       '        QM','  QMFreeze','    QMRain','       SWE',' MeltWater','MeltWStore','    SnowCh', &   !snow (73-79)
                       'SnowRPaved','SnowRBldgs' /)                                          !snow-related (removal) (80-81)

     UnitsAll(:) = (/  '        YYYY','         DOY','          HH','          MM','         day', &   !datetime info (1-5)
                       '       W_m-2','       W_m-2','       W_m-2','       W_m-2','        degC', &   !radiation components (6-10)
                       '       W_m-2','       W_m-2','       W_m-2','       W_m-2','       W_m-2', &   !energy fluxes (11-15)
                       '       W_m-2','       W_m-2','       W_m-2', &                                 !energy fluxes (other approaches) (16-18)
                       '          mm','          mm','          mm','          mm','          mm', &   !water balance components (19-23)
                       '          mm','          mm','          mm','          mm','          mm', &   !water balance components cont. (24-28)
                       '          mm','          mm', &                                                !water balance components cont. (29-30)
                       '          mm','          mm','          mm','          mm','          mm', &   !runoff components (31-35)
                       '          mm','          mm','          mm','          mm', &                  !water use (36-39)
                       '          mm','          mm','          mm','          mm','          mm','          mm', &   !smd for each surface (40-45)
                       '          mm','          mm','          mm','          mm','          mm','          mm','          mm',&   !states (46-52)
                       '         deg','         deg','           -','           -', &                  ! extra radiation info (53-56)
                       '      m2_m-2','           m','           m', &                                 ! extra surface info (57-59)
                       '       m_s-1','           m','       s_m-1','       s_m-1', &                  ! turbulence (60-63)
                       'umol_m-2_s-1', &                                                               ! CO2 flux (64)
                       'umol_m-2_s-1','umol_m-2_s-1','umol_m-2_s-1','umol_m-2_s-1','umol_m-2_s-1', &   ! CO2 flux components (65-69)
                       '       W_m-2','       W_m-2','           -', &                                 ! snow-related (radiation) (70-72)
                       '       W_m-2','       W_m-2','       W_m-2','          mm','          mm','          mm','          mm', &   !snow (73-79)
                       '          mm','          mm' /)                                                !snow-related (removal) (80-81)                    

     ! Use SPREAD function to duplicate repeat formats for groups of similar variables   
     FormatAll(:) = (/ fy,SPREAD(ft,1,3),fd, &      !date time info
                       SPREAD(f94,1,5), &           !radiation components
                       SPREAD(f94,1,5), &           !energy fluxes
                       SPREAD(f94,1,3), &           !energy fluxes (other approaches)
                       SPREAD(f106,1,5), &          !water balance components
                       f106,f104,f106,f106,f94, &   !water balance components cont.
                       f104,f104, &                 !water balance components cont.
                       SPREAD(f106,1,5), &          !runoff components
                       SPREAD(f94,1,4), &           !water use
                       SPREAD(f94,1,6), &           !smd for each surface
                       SPREAD(f94,1,6),f104, &      !state for each surface
                       SPREAD(f94,1,4), &           !extra radiation info
                       SPREAD(f94,1,3), &           !extra surface info
                       f94,f104,f94,f94, &          !turbulence
                       f94, &                       !CO2 flux
                       SPREAD(f94,1,5), &           !CO2 flux components
                       SPREAD(f94,1,3), &           !snow-related (radiation)
                       SPREAD(f106,1,7), &          !snow-related (radiation)
                       SPREAD(f94,1,2) /)           !snow-related (removal)
                       
     ! Set type of aggregation required by wrapper   
     AggAll(:) = (/ SPREAD(aT,1,5), &         !date time info
                    SPREAD(aA,1,5), &         !radiation components
                    SPREAD(aA,1,5), &         !energy fluxes
                    SPREAD(aA,1,3), &         !energy fluxes (other approaches)
                    SPREAD(aS,1,5), &         !water balance components
                    aS,aL,aL,aS,aL, &         !water balance components cont.
                    aS,aS, &                  !water balance components cont.
                    SPREAD(aS,1,5), &         !runoff components
                    SPREAD(aS,1,4), &         !water use
                    SPREAD(aL,1,6), &         !smd for each surface
                    SPREAD(aL,1,7), &         !state for each surface
                    aL,aL,aA,aA, &            !extra radiation info
                    SPREAD(aA,1,3), &         !extra surface info
                    SPREAD(aA,1,4), &         !turbulence
                    aA, &                     !CO2 flux
                    SPREAD(aA,1,5), &         !CO2 flux components
                    SPREAD(aA,1,3), &         !snow-related (radiation)
                    aA,aA,aA,aS,aS,aS,aS, &   !snow-related (radiation)
                    SPREAD(aS,1,2) /)         !snow-related (removal)                       

     ! Select variables to be written out                  
     IF(WriteOutOption == 0) THEN   !all (not snow-related)
        UsecolumnsDataOut = (/ (i, i=1,69, 1) /)
     ELSEIF(WriteOutOption == 1) THEN   !all plus snow-related 
        UsecolumnsDataOut = (/ (i, i=1,nColumnsDataOut, 1) /)
     ELSEIF(WriteOutOption == 2) THEN   !minimal output
        UsecolumnsDataOut = (/ (i, i=1,15, 1),(i, i=19,28, 1), 53,54,55,56, 57, 60,61, 64 /) 
     ELSE
        write(*,*) 'RunControl: WriteOutOption code not recognised, so writing out all variables.'    
        UsecolumnsDataOut = (/ (i, i=1,69, 1) /)
     ENDIF

     ! Create subset of HeaderAll and FormatAll for selected variables only   
     HeaderOut=''
     FormatOut=''
     UnitsOut=''
     AggOut=''
     ColNos=''
     DO i=1,SIZE(UseColumnsDataOut)
        HeaderOut=trim(HeaderOut)//' '//adjustl(HeaderAll(UsecolumnsDataOut(i)))
        !write(*,*) HeaderOut
        UnitsOut=trim(UnitsOut)//' '//adjustl(UnitsAll(UsecolumnsDataOut(i)))
        !write(*,*) UnitsOut
        FormatOut=trim(FormatOut)//' '//adjustl(FormatAll(UsecolumnsDataOut(i)))
        !write(*,*) FormatOut
        AggOut=trim(AggOut)//' '//adjustl(AggAll(UsecolumnsDataOut(i)))
        !write(*,*) AggOut
        write(itext,'(i3)') i
        ColNos=trim(ColNos)//' '//adjustl(itext)
     ENDDO  
     !HeaderUse=trim(adjustl(HeaderOut))//' ' !with extra space at end of header row
     ALLOCATE(CHARACTER(LEN(trim(adjustl(HeaderOut)))):: HeaderUse)
     ALLOCATE(CHARACTER(LEN(trim(adjustl(UnitsOut)))):: UnitsUse)
     ALLOCATE(CHARACTER(LEN(trim(adjustl(FormatOut)))):: FormatUse)
     ALLOCATE(CHARACTER(LEN(trim(adjustl(AggOut)))):: AggUse)
     ALLOCATE(CHARACTER(LEN(trim(adjustl(ColNos)))):: ColNosUse)
     HeaderUse=trim(adjustl(HeaderOut))
     UnitsUse=trim(adjustl(UnitsOut))
     FormatUse='('//trim(adjustl(FormatOut))//')'     
     AggUse=trim(adjustl(AggOut))
     ColNosUse=trim(adjustl(ColNos))
     !write(*,*) '||',HeaderUse,'||'
     !write(*,*) '||',FormatUse,'||'  
     !write(*,*) '||',ColNosUse,'||'
         

    !=========== Write output format info to file ===========  
     OPEN(50,file=TRIM(FileOutFormat),err=111)
     WRITE(50,'(a)') ColNosUse
     WRITE(50,'(a)') HeaderUse
     WRITE(50,'(a)') UnitsUse
     WRITE(50,'(a)') FormatUse(2:(LEN(FormatUse)-1))   !also write formats to output file (without outer brackets)
     WRITE(50,'(a)') AggUse
     CLOSE (50)
     OutputFormats = 0
  ENDIF
  
  !========== Open output file (and first time print header) ==========
  
  ! Main output file --------------------------------------------------
  lfnOutC=39  !Output file code
  IF (iv==1) THEN
     OPEN(lfnOutC,file=TRIM(FileOut),err=112)
     WRITE(lfnOutC,'(a)') HeaderUse
  ELSE
     OPEN(lfnOutC,file=TRIM(FileOut),position='append')!,err=112)
  ENDIF
  
  ! SOLWEIG output file -----------------------------------------------
  IF (SOLWEIGpoi_out==1) THEN
     OPEN(9,file=SOLWEIGpoiOut)
     WRITE(9,113)
113  FORMAT('%doy dectime  azimuth altitude GlobalRad DiffuseRad DirectRad ', &
            ' Kdown2d    Kup2d    Ksouth     Kwest    Knorth     Keast ', &
            ' Ldown2d    Lup2d    Lsouth     Lwest    Lnorth     Least ', &
            '   Tmrt       I0       CI        gvf      shadow    svf    svfbuveg    Ta    Tg')
  ENDIF
  
  ! BL ouput file -----------------------------------------------------
  IF (CBLuse>=1) THEN
     OPEN(53,file=BLOut,status='unknown')
     WRITE(53, 102)
102  FORMAT('iy  id   it imin dectime         z            theta          q',&
            '               theta+          q+              Temp_C          rh',&
            '              QH_use          QE_use          Press_hPa       avu1',&
            '            ustar           avdens          lv_J_kg         avcp',&
            '            gamt            gamq')
  ENDIF

  ! Snow output file --------------------------------------------------
  IF (SnowUse>=1) THEN
     IF (iv==1) THEN
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

  ! ESTM ouput file ---------------------------------------------------
  IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14) THEN
     IF (iv==1) THEN    
        OPEN(58,file=TRIM(ESTMOut),status='unknown')
        WRITE(58, 115)
115     FORMAT('%iy id it imin dectime ',&
               'QSNET QSAIR QSWALL QSROOF QSGROUND QSIBLD ',&
               'TWALL1 TWALL2 TWALL3 TWALL4 TWALL5 ',&       !T0_WALL TWALL1 TWALL2 TWALL3 TN_WALL
               'TROOF1 TROOF2 TROOF3 TROOF4 TROOF5 ',&
               'TGROUND1 TGROUND2 TGROUND3 TGROUND4 TGROUND5 ',&
               'TiBLD1 TiBLD2 TiBLD3 TiBLD4 TiBLD5 TaBLD ')
     ELSE
        OPEN(58,file=TRIM(ESTMOut),position='append')
     ENDIF
  ENDIF

  !These belong to NARP ouput file
     ! 'kup_Paved kup_Bldgs kup_EveTr kup_DecTr kup_Grass kup_BSoil kup_Water ',&
     ! 'lup_Paved lup_Bldgs lup_EveTr lup_DecTr lup_Grass lup_BSoil lup_Water ',&
     ! 'Ts_Paved Ts_Bldgs Ts_EveTr Ts_DecTr Ts_Grass Ts_BSoil Ts_Water ',&
     ! 'qn_Paved qn_Bldgs qn_EveTr qn_DecTr qn_Grass qn_BSoil qn_Water ',&
  
  
  !========== Write out data ==========
  DO i=1,irMax
      WRITE(lfnoutC,FormatUse) INT(dataOut(i,PACK(UseColumnsDataOut, UsecolumnsDataOut < 5),Gridiv)),&
           dataOut(i,PACK(UseColumnsDataOut, UsecolumnsDataOut >= 5),Gridiv)
     !WRITE(lfnoutC,301) (INT(dataOut(i,is,Gridiv)),is=1,4),&      
     !      dataOut(i,5:ncolumnsDataOut,Gridiv)
  
  ENDDO
  
  
  IF (SOLWEIGpoi_out==1) THEN
     DO i=1,SolweigCount-1
        WRITE(9,304) INT(dataOutSOL(i,1,Gridiv)),(dataOutSOL(i,is,Gridiv),is=2,ncolumnsdataOutSOL)
     ENDDO
  ENDIF

  IF(CBLuse>=1) THEN
     DO i=1,iCBLcount
        WRITE(53,305)(INT(dataOutBL(i,is,Gridiv)),is=1,4),(dataOutBL(i,is,Gridiv),is=5,ncolumnsdataOutBL)
     ENDDO
  ENDIF

  IF(SnowUse>=1) THEN
     DO i=1,irmax
        WRITE(54,306)(INT(dataOutSnow(i,is,Gridiv)),is=1,4),(dataOutSnow(i,is,Gridiv),is=5,ncolumnsDataOutSnow)
     ENDDO
  ENDIF
  
  IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14)THEN
     DO i=1,irMax
        WRITE(58, 307)(INT(dataOutESTM(i,is,Gridiv)),is=1,4),(dataOutESTM(i,is,Gridiv),is=5,32)
     ENDDO
  ENDIF

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
111 CALL ErrorHint(52,TRIM(fileOutFormat),notUsed,notUsed,notUsedI)
112 CALL ErrorHint(52,TRIM(fileOut),notUsed,notUsed,notUsedI)


END SUBROUTINE SUEWS_Output
