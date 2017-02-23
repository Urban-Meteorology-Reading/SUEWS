 MODULE MetDisagg
 !========================================================================================    
 ! Disaggregation of meteorological forcing data    
 !  Code to disaggregate met forcing data from resolution provided to the model time-step
 ! Created: HCW 10 Feb 2017    
 !   
 ! ---- Key for MetDisaggMethod settings ----
 !10  ->  linear disaggregation for timestamps representing the end of the averaging period 
 !20  ->  linear disaggregation for instantaneous variables (e.g. air temperature, humidity, pressure, wind speed in WFDEI dataset)
 !100 -> evenly distribute rainfall among all subintervals in a rainy interval
 !101 -> evenly distribute rainfall among RainAmongN subintervals in a rainy interval 
 !        - requires RainAmongN to be set in RunControl.nml
 ! If KdownZen = 1 -> include additional zenith check in kdown disaggregation
 ! 
 ! N.B. wdir downscaling is currently not implemented   
 !
 !========================================================================================

 USE AllocateArray
 USE ColNamesInputFiles
 USE Data_In
 USE Initial

 IMPLICIT NONE

 CONTAINS
 
  !======================================================================================
  SUBROUTINE DisaggregateMet(iBlock, igrid)
  ! Subroutine to disaggregate met forcing data to model time-step     
  ! HCW 10 Feb 2017 
  !======================================================================================
 
    USE Sues_Data
    USE DefaultNotUsed

    IMPLICIT NONE

    INTEGER:: lunit = 100
    INTEGER:: tdiff   !Time difference (in minutes) between first and second rows of original met forcing file
    INTEGER:: i,ii  !counter
    INTEGER:: iBlock, igrid
    INTEGER,DIMENSION(Nper):: seq1Nper
    INTEGER,DIMENSION(nsd):: seq1nsd
    INTEGER,DIMENSION(nColumnsMetForcingData):: MetDisaggMethod   ! Stores method to use for disaggregating met data
    REAL(KIND(1d0)),DIMENSION(nColumnsMetForcingData):: MetArrayOrig
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigMetData*Nper,ncolumnsMetForcingData):: Met_tt    
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigMetData*Nper):: Met_tt_kdownAdj    
    CHARACTER(LEN=9),DIMENSION(ncolumnsMetForcingData):: HeaderMet
    CHARACTER(LEN=10*ncolumnsMetForcingData):: HeaderMetOut
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigMetData):: dectimeOrig 
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigMetData*Nper):: dectimeDscd, dectimeFast 
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigMetData*Nper):: idectime ! sun position at middle of time-step before

    INTEGER, DIMENSION(Nper):: temp_iy, temp_id, temp_ih, temp_im, temp_ihm

    ! Allocate and initialise arrays to receive original forcing data --------------------
    ALLOCATE(MetForDisagg(ReadLinesOrigMetData,nColumnsMetForcingData))
    ALLOCATE(MetForDisaggPrev(nColumnsMetForcingData))
    ALLOCATE(MetForDisaggNext(nColumnsMetForcingData))
    MetForDisagg(:,:) = -999
    MetForDisaggPrev(:) = -999
    MetForDisaggNext(:) = -999
    ! Initialise array to receive disaggregated data
    Met_tt = -999
    ! Intialise disaggregation method
    MetDisaggMethod = -999

    ! Generate useful sequences
    seq1Nper = (/(i, i=1,Nper, 1)/)
    seq1nsd = (/(i, i=1,nsd, 1)/)    

    ! Get methods to use for disaggregation from RunControl
    IF(DiagnoseDisagg==1) write(*,*) 'DisaggMethod: ',DisaggMethod, 'RainDisaggMethod:',RainDisaggMethod, 'RainAmongN:',RainAmongN
    IF(DisaggMethod == 1) THEN
       MetDisaggMethod(:) = 10   !linear disaggregation of averages    
    ELSEIF(DisaggMethod == 2) THEN
       MetDisaggMethod(:) = 20   !linear disaggregation of instantaneous values
    ELSEIF(DisaggMethod == 3) THEN   !WFDEI set up, where T, Q, pres, U are instantaneous
       MetDisaggMethod(:) = 10   !linear disaggregation of averages    
       MetDisaggMethod(10:13) = 20   !linear disagg instantaneous values for U, RH, Tair, pres
    ELSE
       CALL errorHint(2,'Problem in SUEWS_MetDisagg: DisaggMethod value should be 1, 2, or 3', &
                        NotUsed,NotUsed,DisaggMethod)
    ENDIF    
    ! Set rainfall    
    MetDisaggMethod(14) = RainDisaggMethod


    ! Read data ---------------------------------------------------------------------
    IF(DiagnoseDisagg==1) write(*,*) 'Reading file: ', TRIM(FileOrigMet)
    OPEN(lunit,file=TRIM(FileOrigMet),status='old')
    ! CALL skipHeader(lunit,SkipHeaderMet)  !Skip header -> read header instead
    READ(lunit,*) HeaderMet
    !write(*,*) HeaderMet
    ! Skip over lines that have already been read and downscaled
    IF (SkippedLinesOrig>0) THEN
       DO i=1,skippedLinesOrig-1   ! minus 1 here because last line of last block needs to be read again
          READ(lunit,*)
       ENDDO
       ! Read in last line of previous block
       CALL MetRead(lunit,MetArrayOrig,InputmetFormat,ldown_option,NetRadiationMethod,&
                     snowUse,SMDMethod,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)
       MetForDisaggPrev(1:ncolumnsMetForcingData) = MetArrayOrig
    ENDIF    
    ! Read in current block
    DO i=1, ReadLinesOrigMetDataMax
       CALL MetRead(lunit,MetArrayOrig,InputmetFormat,ldown_option,NetRadiationMethod,&
                   snowUse,SMDMethod,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)
       MetForDisagg(i,1:ncolumnsMetForcingData) = MetArrayOrig
    ENDDO
    ! Read in first line of next block (except for last block)
    IF(iBlock/=ReadBlocksOrigMetData) THEN
       CALL MetRead(lunit,MetArrayOrig,InputmetFormat,ldown_option,NetRadiationMethod,&
                      snowUse,SMDMethod,SoilDepthMeas,SoilRocks,SoilDensity,SmCap)
       MetForDisaggNext(1:ncolumnsMetForcingData) = MetArrayOrig
    ENDIF
    CLOSE(lunit)

    ! Check resolution of original met forcing data -------------------------------------
    ! Find time difference (in minutes) between first and second row
    tdiff = INT(MetForDisagg(2,4)-MetForDisagg(1,4))   !Try using minutes
    IF(tdiff == 0) tdiff = INT(MetForDisagg(2,3)-MetForDisagg(1,3))*60   !If no difference in minutes, try using hours
    IF(tdiff < 0) THEN   !If time difference is negative (e.g. change of day), instead use second and third row
       tdiff = INT(MetForDisagg(3,4)-MetForDisagg(2,4))
       IF(tdiff == 0) tdiff = INT(MetForDisagg(3,3)-MetForDisagg(2,3))*60   !If no difference in minutes, try using hours
    ENDIF    
    ! Check actual resolution matches specified input resolution
    IF(tdiff /= ResolutionFilesIn/60) THEN
       CALL errorHint(2,'Problem in SUEWS_MetDisagg: timestamps in met forcing file inconsistent with ResolutionFilesIn', &
                        REAL(ResolutionFilesIn,KIND(1d0)),NotUsed,tdiff*60)
    ENDIF

    ! Disaggregate time columns ---------------------------------------------------------
    write(*,*) 'Disaggregating met forcing data (',TRIM(FileOrigMet),') to model time-step...'
    ! Convert to dectime
    dectimeOrig = MetForDisagg(:,2) + MetForDisagg(:,3)/24.0 + MetForDisagg(:,4)/(60.0*24.0)

    DO i=1,ReadLinesOrigMetDataMax
       ! Downscale dectime using dectimeOrig(i) [becomes timestamp of last subinterval]
       dectimeDscd(Nper*(i-1)+Seq1Nper) = dectimeOrig(i) - (tstep/60.0)/(60.0*24.0)*(/(ii, ii=(Nper-1),0, -1)/)  
       ! Convert to required formats
       temp_iy   = INT(MetForDisagg(i,1))   !Copy year
       temp_id   = FLOOR(dectimeDscd(Nper*(i-1)+Seq1Nper))   !DOY
       ! To avoid precision errors, round here
       !  - this should not be a problem as a difference of 1 = 1 min, so a difference of 0.001 << 1 min
       temp_ihm  = NINT(((dectimeDscd(Nper*(i-1)+Seq1Nper) - temp_id/1.0)*60.0*24.0)*1000.0)/1000   !Minutes of the day (1440 max)
       temp_ih = (temp_ihm-MOD(temp_ihm,60))/60   !Hours             
       temp_im = MOD(temp_ihm,60)   !Minutes      

       IF(dectimeOrig(i) == 1.0000 .and. i > 1) THEN   !If year changes and it is not the beginning of the dataset
         write(*,*) 'Year change encountered: ', dectimeOrig(i), dectimeOrig(i-1)
         ! Re-downscale dectime using dectimeOrig(i-1)
         dectimeDscd(Nper*(i-1)+Seq1Nper) = dectimeOrig(i-1) + (tstep/60.0)/(60.0*24.0)*Seq1Nper  
         ! Convert to required formats
         temp_iy   = INT(MetForDisagg(i,1))   !Copy year
         temp_id   = FLOOR(dectimeDscd(Nper*(i-1)+Seq1Nper))   !DOY
         temp_ihm  = NINT(((dectimeDscd(Nper*(i-1)+Seq1Nper) - temp_id/1.0)*60.0*24.0)*1000.0)/1000   !Mins of the day (1440 max)
         temp_ih = (temp_ihm-MOD(temp_ihm,60))/60   !Hours             
         temp_im = MOD(temp_ihm,60)   !Minutes      
         ! Adjust year and DOY to account for year change
         temp_iy(1:(Nper-1)) = temp_iy(1:(Nper-1)) - 1  !Subtract 1 from year for all except final timestamp
         temp_id(Nper) = 1  !Set day for final timestamp to 1
       ENDIF       

       !IF(i==1 .or. i== ReadlinesOrigMetDataMax) THEN
       !   write(*,*) temp_iy
       !   write(*,*) temp_id
       !   !write(*,*) temp_ihm
       !   write(*,*) temp_ih
       !   write(*,*) temp_im
       !ENDIF   

       ! Copy to Met_tt array
       Met_tt(Nper*(i-1)+Seq1Nper,1) = temp_iy
       Met_tt(Nper*(i-1)+Seq1Nper,2) = temp_id
       Met_tt(Nper*(i-1)+Seq1Nper,3) = temp_ih
       Met_tt(Nper*(i-1)+Seq1Nper,4) = temp_im

    ENDDO

    ! Disaggregate other columns --------------------------------------------------------
    DO ii=5,ncolumnsMetForcingData
    IF(ii == 14) THEN  !Do something different for rainfall and snowfall (if present)
      IF(MetDisaggMethod(14) == 100) THEN
         Met_tt(:,14) = DisaggP_amongN(MetForDisagg(:,14),Nper,Nper,ReadLinesOrigMetData,ReadLinesOrigMetDataMax)
         IF(ALL(MetForDisagg(:,16)==-999)) THEN
            Met_tt(:,16) = -999
         ELSE
            Met_tt(:,16) = DisaggP_amongN(MetForDisagg(:,16),Nper,Nper,ReadLinesOrigMetData,ReadLinesOrigMetDataMax) 
         ENDIF   
      ELSEIF(MetDisaggMethod(14) == 101) THEN
         IF(RainAmongN == -999) THEN
            CALL ErrorHint(2,'Problem in SUEWS_MetDisagg: RainDisaggMethod requires RainAmongN', &
                               REAL(RainAmongN,KIND(1d0)),NotUsed,RainDisaggMethod)
         ELSEIF(RainAmongN > Nper) THEN
            CALL ErrorHint(2,'Problem in SUEWS_MetDisagg: RainAmongN > Nper',REAL(Nper,KIND(1d0)),NotUsed,RainAmongN)
         ELSE
            Met_tt(:,14) = DisaggP_amongN(MetForDisagg(:,14),RainAmongN, Nper,ReadLinesOrigMetData,ReadLinesOrigMetDataMax)
            IF(ALL(MetForDisagg(:,16)==-999)) THEN
               Met_tt(:,16) = -999
            ELSE
               Met_tt(:,16) = DisaggP_amongN(MetForDisagg(:,16),RainAmongN,Nper,ReadLinesOrigMetData,ReadLinesOrigMetDataMax)
            ENDIF
         ENDIF      
      ELSE
         write(*,*) 'Disaggregation code for rain not recognised'    
      ENDIF   
    ELSEIF(ii == 24) THEN  !wind direction disaggregation not coded yet...
      IF(ANY(MetForDisagg(:,ii)/=-999)) THEN
         write(*,*) 'Disaggregation of wind direction not currently implemented!'
      ENDIF
    ELSE   
      IF(ALL(MetForDisagg(:,ii)==-999)) THEN
         !IF(DiagnoseDisagg==1) write(*,*) 'No data for col.', ii  
         Met_tt(:,ii) = -999
      ELSE
         Met_tt(:,ii) = Disagg_Lin(MetForDisagg(:,ii),MetForDisaggPrev(ii),MetForDisaggNext(ii),MetDisaggMethod(ii), & 
                                      Nper,ReadLinesOrigMetData,ReadLinesOrigMetDataMax,iBlock)
      ENDIF
    ENDIF
    ENDDO

    ! Adjust kdown disaggregation using zenith angle
    IF(KdownZen == 1) THEN
    IF(DiagnoseDisagg==1) write(*,*) 'Adjusting disaggregated kdown using zenith angle' 
    Met_tt_kdownAdj(:) = Met_tt(:,15) 
    ! Translate location data from SurfaceChar to find solar angles
    lat = SurfaceChar(igrid,c_lat)
    lng = SurfaceChar(igrid,c_lng)*(-1.0)  !HCW switched sign of lng 12 Dec 2016. Input should now be -ve for W, +ve for E
    timezone = SurfaceChar(igrid,c_tz)
    alt = SurfaceChar(igrid,c_Alt)
    ! Calculate dectime at downscaled time-step
    dectimeFast(:) = Met_tt(:,2) + Met_tt(:,3)/24.0 + Met_tt(:,4)/(60.0*24.0)
    idectime=dectimeFast-halftimestep! sun position at middle of timestep before
    DO i=1,(ReadLinesOrigMetDataMax*Nper)
      CALL sun_position(Met_tt(i,2),idectime(i),timezone,lat,lng,alt,azimuth,zenith_deg)
      ! If sun below horizon, set disaggregated kdown to zero
      IF(zenith_deg > 90) THEN
         !write(*,*) Met_tt(i,1:4)
         Met_tt_kdownAdj(i) = 0.0
      ENDIF    
    ENDDO
    ! Redistribute kdown over each day
    DO i=1,(ReadLinesOrigMetDataMax*Nper/nsd) ! Loop over each day
      Met_tt_kdownAdj((i-1)*nsd+seq1nsd) = Met_tt_kdownAdj( (i-1)*nsd+seq1nsd) * &
                 SUM(Met_tt((i-1)*nsd+seq1nsd,15 ))/SUM(Met_tt_kdownAdj((i-1)*nsd+seq1nsd))   
    ENDDO
    ! Copy adjusted kdown back to Met_tt array
    Met_tt(:,15) = Met_tt_kdownAdj(:)
    ENDIF

    ! Copy disaggregated data to MetForcingDataArray
    MetForcingData(:,1:24,GridCounter) = Met_tt(:,1:24) 

    ! If snow is -999, set to zero (also in LUMPS_metRead.f95)
    IF(ALL(MetForcingData(:,16,GridCounter) == -999)) MetForcingData(:,16,GridCounter)=0

    ! Undo pressure conversion again for writing out
    Met_tt(:,13) = Met_tt(:,13)/10.0

    ! Write out disaggregated file ------------------------------------------------------
    IF(KeepTstepFilesIn == 1) THEN
    IF (iBlock==1) THEN
      ! Prepare header     
      DO i=1,ncolumnsMetForcingData
         IF(i==1) THEN
            HeaderMetOut=ADJUSTL(HeaderMet(i))
         ELSE
            HeaderMetOut=TRIM(HeaderMetOut)//' '//ADJUSTL(HeaderMet(i))
         ENDIF
      ENDDO
      ! Write out header
      OPEN(78,file=TRIM(FileDscdMet),err=112)
      WRITE(78,'(a)') HeaderMetOut
    ELSE
      OPEN(78,file=TRIM(FileDscdMet),position='append')!,err=112)
    ENDIF
    ! Write out data
    DO i=1,(ReadLinesOrigMetDataMax*Nper)
      WRITE(78,303) (INT(Met_tt(i,ii)), ii=1,4), Met_tt(i,5:ncolumnsMetForcingData)   
    ENDDO
    IF(iBlock == ReadBlocksOrigMetData) THEN
     WRITE(78,'(i2)') -9
     WRITE(78,'(i2)') -9
    ENDIF
    CLOSE (78)   !Close output file
    ENDIF 


    303 FORMAT((i4,1X), 3(i3,1X), 9(f9.2,1X), (f9.4,1X), 10(f9.2,1X))  !Allows 4 dp for rainfall

    ! Deallocate arrays -----------------------------------------------------------------
    DEALLOCATE(MetForDisagg)
    DEALLOCATE(MetForDisaggPrev)
    DEALLOCATE(MetForDisaggNext)

    RETURN      

    112 CALL ErrorHint(52,TRIM(FileDscdMet),notUsed,notUsed,notUsedI)

  ENDSUBROUTINE DisaggregateMet
!======================================================================================        

  !======================================================================================
  SUBROUTINE DisaggregateESTM(iBlock, igrid)
  ! Subroutine to disaggregate met forcing data to model time-step     
  ! HCW 10 Feb 2017 
  !======================================================================================
 
    USE Sues_Data
    USE DefaultNotUsed

    IMPLICIT NONE

    INTEGER:: lunit = 101
    INTEGER:: tdiff   !Time difference (in minutes) between first and second rows of original met forcing file
    INTEGER:: i,ii  !counter
    INTEGER:: iBlock, igrid
    INTEGER,DIMENSION(NperESTM):: seq1NperESTM
    INTEGER,DIMENSION(nsd):: seq1nsd
    INTEGER,DIMENSION(ncolsESTMdata):: ESTMDisaggMethod   ! Stores method to use for disaggregating met data
    REAL(KIND(1d0)),DIMENSION(ncolsESTMdata):: ESTMArrayOrig
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigESTMData*NperESTM,ncolsESTMdata):: ESTM_tt    
    CHARACTER(LEN=9),DIMENSION(ncolsESTMdata):: HeaderESTM
    CHARACTER(LEN=10*ncolsESTMdata):: HeaderESTMOut
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigESTMData):: dectimeOrig
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigESTMData*NperESTM):: dectimeDscd, dectimeFast 
    REAL(KIND(1d0)),DIMENSION(ReadLinesOrigESTMData*NperESTM):: idectime ! sun position at middle of time-step before
    INTEGER::iostat_var
    
    INTEGER, DIMENSION(NperESTM):: temp_iy, temp_id, temp_ih, temp_im, temp_ihm

    ! Allocate and initialise arrays to receive original forcing data --------------------
    ALLOCATE(ESTMForDisagg(ReadLinesOrigESTMData,ncolsESTMdata))
    ALLOCATE(ESTMForDisaggPrev(ncolsESTMdata))
    ALLOCATE(ESTMForDisaggNext(ncolsESTMdata))
    ESTMForDisagg(:,:) = -999
    ESTMForDisaggPrev(:) = -999
    ESTMForDisaggNext(:) = -999
    ! Initialise array to receive disaggregated data
    ESTM_tt = -999
    ! Intialise disaggregation method
    ESTMDisaggMethod = -999

    ! Generate useful sequences
    seq1NperESTM = (/(i, i=1,NperESTM, 1)/)
    seq1nsd = (/(i, i=1,nsd, 1)/)    

    ! Get methods to use for disaggregation from RunControl 
    ! (N.B.DisaggMethodESTM is set as 1 or 2 in RunControl; ESTMDisaggMethod is array of ncolsESTMdata used here)
    IF(DisaggMethodESTM == 1) THEN
       ESTMDisaggMethod(:) = 10   !linear disaggregation of averages    
    ELSEIF(DisaggMethodESTM == 2) THEN
       ESTMDisaggMethod(:) = 20   !linear disaggregation of instantaneous values
    ELSE
       CALL errorHint(2,'Problem in SUEWS_ESTMDisagg: DisaggMethodESTM value should be 1 or 2', &
                        NotUsed,NotUsed,DisaggMethodESTM)
    ENDIF    
    
    ! Read data ---------------------------------------------------------------------
    IF(DiagnoseDisaggESTM==1) write(*,*) 'Reading file: ', TRIM(FileOrigESTM)
    OPEN(lunit,file=TRIM(FileOrigESTM),status='old')
    ! CALL skipHeader(lunit,SkipHeaderMet)  !Skip header -> read header instead
    READ(lunit,*) HeaderESTM
    !write(*,*) HeaderMet
    ! Skip over lines that have already been read and downscaled
    IF (SkippedLinesOrigESTM>0) THEN
       DO i=1,skippedLinesOrigESTM-1   ! minus 1 here because last line of last block needs to be read again
          READ(lunit,*)
       ENDDO
       ! Read in last line of previous block
       READ(lunit,*,iostat=iostat_var) ESTMArrayOrig
       ESTMForDisaggPrev(1:ncolsESTMdata) = ESTMArrayOrig
    ENDIF    
    ! Read in current block
    DO i=1, ReadLinesOrigESTMDataMax
       READ(lunit,*,iostat=iostat_var) ESTMArrayOrig
       ESTMForDisagg(i,1:ncolsESTMdata) = ESTMArrayOrig
    ENDDO
    ! Read in first line of next block (except for last block)
    IF(iBlock/=ReadBlocksOrigMetData) THEN
       READ(lunit,*,iostat=iostat_var) ESTMArrayOrig
       ESTMForDisaggNext(1:ncolsESTMdata) = ESTMArrayOrig
    ENDIF
    CLOSE(lunit)

    ! Check resolution of original met forcing data -------------------------------------
    ! Find time difference (in minutes) between first and second row
    tdiff = INT(ESTMForDisagg(2,4)-ESTMForDisagg(1,4))   !Try using minutes
    IF(tdiff == 0) tdiff = INT(ESTMForDisagg(2,3)-ESTMForDisagg(1,3))*60   !If no difference in minutes, try using hours
    IF(tdiff < 0) THEN   !If time difference is negative (e.g. change of day), instead use second and third row
       tdiff = INT(ESTMForDisagg(3,4)-ESTMForDisagg(2,4))
       IF(tdiff == 0) tdiff = INT(ESTMForDisagg(3,3)-ESTMForDisagg(2,3))*60   !If no difference in minutes, try using hours
    ENDIF    
    ! Check actual resolution matches specified input resolution
    IF(tdiff /= ResolutionFilesInESTM/60) THEN
       CALL errorHint(2,'Problem in SUEWS_ESTMDisagg: timestamps in ESTM forcing file inconsistent with ResolutionFilesInESTM', &
                        REAL(ResolutionFilesInESTM,KIND(1d0)),NotUsed,tdiff*60)
    ENDIF

    ! Disaggregate time columns ---------------------------------------------------------
    write(*,*) 'Disaggregating ESTM forcing data (',TRIM(FileOrigESTM),') to model time-step...'
    ! Convert to dectime
    dectimeOrig = ESTMForDisagg(:,2) + ESTMForDisagg(:,3)/24.0 + ESTMForDisagg(:,4)/(60.0*24.0)

    DO i=1,ReadLinesOrigESTMDataMax
       ! Downscale dectime using dectimeOrig(i) [becomes timestamp of last subinterval]
       dectimeDscd(NperESTM*(i-1)+Seq1NperESTM) = dectimeOrig(i) - (tstep/60.0)/(60.0*24.0)*(/(ii, ii=(NperESTM-1),0, -1)/)  
       ! Convert to required formats
       temp_iy   = INT(ESTMForDisagg(i,1))   !Copy year
       temp_id   = FLOOR(dectimeDscd(NperESTM*(i-1)+Seq1NperESTM))   !DOY
       ! To avoid precision errors, round here
       !  - this should not be a problem as a difference of 1 = 1 min, so a difference of 0.001 << 1 min
       temp_ihm  = NINT(((dectimeDscd(NperESTM*(i-1)+Seq1NperESTM) - temp_id/1.0)*60.0*24.0)*1000.0)/1000   !Minutes of the day (1440 max)
       temp_ih = (temp_ihm-MOD(temp_ihm,60))/60   !Hours             
       temp_im = MOD(temp_ihm,60)   !Minutes      

       IF(dectimeOrig(i) == 1.0000 .and. i > 1) THEN   !If year changes and it is not the beginning of the dataset
         write(*,*) 'Year change encountered: ', dectimeOrig(i), dectimeOrig(i-1)
         ! Re-downscale dectime using dectimeOrig(i-1)
         dectimeDscd(NperESTM*(i-1)+Seq1NperESTM) = dectimeOrig(i-1) + (tstep/60.0)/(60.0*24.0)*Seq1NperESTM  
         ! Convert to required formats
         temp_iy   = INT(ESTMForDisagg(i,1))   !Copy year
         temp_id   = FLOOR(dectimeDscd(NperESTM*(i-1)+Seq1NperESTM))   !DOY
         temp_ihm  = NINT(((dectimeDscd(NperESTM*(i-1)+Seq1NperESTM) - temp_id/1.0)*60.0*24.0)*1000.0)/1000   !Mins of the day (1440 max)
         temp_ih = (temp_ihm-MOD(temp_ihm,60))/60   !Hours             
         temp_im = MOD(temp_ihm,60)   !Minutes      
         ! Adjust year and DOY to account for year change
         temp_iy(1:(NperESTM-1)) = temp_iy(1:(NperESTM-1)) - 1  !Subtract 1 from year for all except final timestamp
         temp_id(NperESTM) = 1  !Set day for final timestamp to 1
       ENDIF       

       !IF(i==1 .or. i== ReadlinesOrigESTMDataMax) THEN
       !   write(*,*) temp_iy
       !   write(*,*) temp_id
       !   !write(*,*) temp_ihm
       !   write(*,*) temp_ih
       !   write(*,*) temp_im
       !ENDIF   

       ! Copy to ESTM_tt array
       ESTM_tt(NperESTM*(i-1)+Seq1NperESTM,1) = temp_iy
       ESTM_tt(NperESTM*(i-1)+Seq1NperESTM,2) = temp_id
       ESTM_tt(NperESTM*(i-1)+Seq1NperESTM,3) = temp_ih
       ESTM_tt(NperESTM*(i-1)+Seq1NperESTM,4) = temp_im
       
    ENDDO

    
    ! Disaggregate other columns --------------------------------------------------------
    ! All other columns are temperatures
    DO ii=5,ncolsESTMdata
       IF(ALL(ESTMForDisagg(:,ii)==-999)) THEN
          !IF(DiagnoseDisaggESTM==1) write(*,*) 'No data for col.', ii  
          ESTM_tt(:,ii) = -999
       ELSE
          ESTM_tt(:,ii) = Disagg_Lin(ESTMForDisagg(:,ii),ESTMForDisaggPrev(ii),ESTMForDisaggNext(ii),ESTMDisaggMethod(ii), &
                                        NperESTM,ReadLinesOrigESTMData,ReadLinesOrigESTMDataMax,iBlock)
      ENDIF
    ENDDO
    
    ! Copy disaggregated data to MetForcingDataArray
    ESTMForcingData(:,1:ncolsESTMdata,GridCounter) = ESTM_tt(:,1:ncolsESTMdata) 

    ! Write out disaggregated file ------------------------------------------------------
    IF(KeepTstepFilesIn == 1) THEN
       IF (iBlock==1) THEN
         ! Prepare header     
         DO i=1,ncolsESTMdata
            IF(i==1) THEN
               HeaderESTMOut=ADJUSTL(HeaderESTM(i))
            ELSE
               HeaderESTMOut=TRIM(HeaderESTMOut)//' '//ADJUSTL(HeaderESTM(i))
            ENDIF
         ENDDO
         ! Write out header
         OPEN(78,file=TRIM(FileDscdESTM),err=113)
         WRITE(78,'(a)') HeaderESTMOut
       ELSE
          OPEN(78,file=TRIM(FileDscdESTM),position='append')!,err=113)
       ENDIF
       ! Write out data
       DO i=1,(ReadLinesOrigESTMDataMax*NperESTM)
          WRITE(78,304) (INT(ESTM_tt(i,ii)), ii=1,4), ESTM_tt(i,5:ncolsESTMdata)   
       ENDDO
       IF(iBlock == ReadBlocksOrigMetData) THEN
          WRITE(78,'(i2)') -9
          WRITE(78,'(i2)') -9
       ENDIF
       CLOSE (78)   !Close output file
    ENDIF 


    304 FORMAT((i4,1X), 3(i3,1X), 9(f9.4,1X)) 

    ! Deallocate arrays -----------------------------------------------------------------
    DEALLOCATE(ESTMForDisagg)
    DEALLOCATE(ESTMForDisaggPrev)
    DEALLOCATE(ESTMForDisaggNext)

    RETURN      

    113 CALL ErrorHint(52,TRIM(FileDscdESTM),notUsed,notUsed,notUsedI)

  ENDSUBROUTINE DisaggregateESTM
!======================================================================================          
  

! Define functions here:
FUNCTION Disagg_Lin(Slow,SlowPrev,SlowNext,DisaggType,Nper_loc,ReadLinesOrig_loc,ReadLinesOrigMax_loc,iBlock) RESULT(Fast)

  USE DefaultNotUsed
  USE sues_data

  IMPLICIT NONE

  INTEGER:: DisaggType   !Type of disaggregation: 10 for averaged variables; 20 for instantaneous variables
  INTEGER:: Nper_loc     !Number of subintervals per interval (local Nper)
  INTEGER:: ReadLinesOrig_loc,ReadLinesOrigMax_loc   !Number of lines to read in original file (local)
  INTEGER:: iBlock
  REAL(KIND(1d0)),DIMENSION(ReadLinesOrig_loc*Nper_loc):: Fast  !Array to receive disaggregated data
  REAL(KIND(1d0)),DIMENSION(ReadLinesOrig_loc):: Slow   !Array to disaggregate
  REAL(KIND(1d0)):: SlowPrev, SlowNext
  INTEGER,DIMENSION(Nper_loc):: FastRows   !Group of rows that are filled with each iteration
  INTEGER,DIMENSION(FLOOR(Nper_loc/2.0)):: FirstRows10   !Rows at the beginning that are not filled during iteration (for averages)
  INTEGER,DIMENSION(Nper_loc-FLOOR(Nper_loc/2.0)):: LastRows10    !Rows at the end that are not filled during iteration
  INTEGER,DIMENSION(Nper_loc):: FirstRows20   !Rows at the beginning that are not filled during iteration (for instantaneous)
  INTEGER,DIMENSION(Nper_loc):: seq1Nper_loc   !1 to Nper_loc
  INTEGER:: XNper_loc   !XNper_loc = 2 for even Nper_loc; XNper_loc=1 for odd Nper_loc    
  INTEGER:: i,ii   !counters

  ! Calculate XNper_loc (differentiates between disaggregations with odd and even Nper_loc)
  IF(MOD(Nper_loc,2)==0) XNper_loc=2
  IF(MOD(Nper_loc,2)==1) XNper_loc=1

  seq1Nper_loc = (/(i, i=1,Nper_loc, 1)/)

  ! Setup counters for iteration 
  IF(DisaggType==10) THEN
     FastRows = FLOOR(Nper_loc/2.0) + seq1Nper_loc  ! Rows to create at model time-step
     FirstRows10 = (/(i, i=1,(FastRows(1)-1), 1)/)   !For start of dataset
     LastRows10 =  (/(i, i=Nper_loc*(ReadLinesOrigMax_loc-1-1)+FastRows(Nper_loc)+1,(ReadLinesOrigMax_loc*Nper_loc),1)/)  ! For end of dataset
  ELSEIF(DisaggType==20) THEN
     FastRows = Nper_loc + seq1Nper_loc   !Rows to create at model time-step    
     FirstRows20 = (/(i, i=1,(FastRows(1)-1), 1)/)   !For start of dataset
  ENDIF

  ! Initialise fast array to -999
  Fast = -999
  ! Linearly disaggregate 
  IF(DisaggType==10) THEN   !Averaged variables
     IF(DiagnoseDisagg==1) write(*,*) 'Linearly disaggregating averaged variable'
     DO i=1,(ReadLinesOrigMax_loc-1) 
        Fast(Nper_loc*(i-1)+FastRows) = Slow(i) - &
                                           (Slow(i+1)-Slow(i))/(XNper_loc*Nper_loc) + &
                                              (Slow(i+1)-Slow(i))/Nper_loc*(/(ii, ii=1,Nper_loc, 1)/)                                              
     ENDDO

     ! For first few rows, use previous met block
     IF(iBlock==1) THEN
        Fast(FirstRows10) = Fast(FastRows(1))   !Use repeat values at the start of the year
     ELSE 
        Fast(FirstRows10) = SlowPrev - &
                             (Slow(1)-SlowPrev)/(XNper_loc*Nper_loc) + &
                               (Slow(1)-SlowPrev)/Nper_loc * & 
                                  (/(ii, ii=(Nper_loc-SIZE(FirstRows10)+1),Nper_loc, 1)/)
     ENDIF
     ! For last few rows, use next met block
     IF(iBlock==ReadBlocksOrigMetData) THEN
        Fast(LastRows10) = Fast(Nper_loc*(ReadLinesOrigMax_loc-1-1)+FastRows(Nper_loc))   !Use repeat values at the end of the year
     ELSE
        Fast(LastRows10) = Slow(ReadLinesOrigMax_loc) - &
                            (SlowNext-Slow(ReadLinesOrigMax_loc))/(XNper_loc*Nper_loc) + &
                              (SlowNext-Slow(ReadLinesOrigMax_loc))/Nper_loc * &
                                (/(ii, ii=1,SIZE(LastRows10), 1)/)
     ENDIF
  ELSEIF(DisaggType==20) THEN   !Instantaneous variables
     IF(DiagnoseDisagg==1) write(*,*) 'Linearly disaggregating instantaneous variable'
     DO i=1,(ReadLinesOrigMax_loc-1) 
        Fast(Nper_loc*(i-1)+FastRows) = (Slow(i) + &
                                          (Slow(i+1)-Slow(i))/Nper_loc*2*(seq1Nper_loc-1) + &
                                              Slow(i))/2
     ENDDO
     ! For first few rows, use previous met block
     IF(iBlock==1) THEN
        Fast(FirstRows20) = Fast(FastRows(1))   !Use repeat values at the start of the year
     ELSE
        Fast(FirstRows20) = (SlowPrev + &
                              (Slow(1)-SlowPrev)/Nper_loc*2 * &
                                ((/(ii, ii=(Nper_loc-SIZE(FirstRows20)+1),Nper_loc, 1)/)-1) + &
                                   SlowPrev)/2
     ENDIF
     !! Last few rows are already filled for the instantaneous value disaggregation
     !IF(iBlock==ReadBlocksOrigMetData) THEN
     !   Fast(LastRows20) = Fast(Nper_loc*(ReadLinesOrigMax_loc-1-1)+FastRows(Nper_loc))   !Use repeat values at the end of the year
     !ELSE
     !   Fast(LastRows20) = (Slow(ReadLinesOrigMax_loc) + &
     !                     (SlowNext-Slow(ReadLinesOrigMax_loc))/Nper_loc*2 * &
     !                        ((/(ii, ii=1,SIZE(LastRows20), 1)/)-1) + &
     !                       Slow(ReadLinesOrigMax_loc))/2
     !ENDIF
  ENDIF

  IF(ANY(Fast(1:ReadLinesOrigMax_loc*Nper_loc) == -999)) THEN
     WRITE(*,*) 'Problem: -999s (',COUNT(Fast(1:ReadLinesOrigMax_loc*Nper_loc)==-999),') in disaggregated data.'
     CALL errorHint(13,'Problem in SUEWS_MetDisagg: -999 values in disaggregated data.',NotUsed,NotUsed,NotUsedI)
  ENDIF    

ENDFUNCTION Disagg_Lin  
!======================================================================================

!======================================================================================
FUNCTION DisaggP_amongN(Slow,amongN, Nper_loc, ReadLinesOrig_loc, ReadLinesOrigMax_loc) RESULT(Fast)
! Subroutine to disaggregate precipitation by evenly distributing among N subintervals
!  (i.e. equal intensity in N subintervals) 
! See Ward et al. (in review), meanN, 0.5N or 0.25N approach    
! HCW 10 Feb 2017 
!======================================================================================

  USE DefaultNotUsed
  USE sues_data

  IMPLICIT NONE

  INTEGER:: MetCol       !Column for met variable to disaggregate, type of disaggregation 
  INTEGER:: amongN       !Number of subintervals over which rain will be distributed
  INTEGER:: Nper_loc     !Number of subintervals per interval (local Nper)
  INTEGER:: ReadLinesOrig_loc,ReadLinesOrigMax_loc   !Number of lines to read in original file (local)
  REAL(KIND(1d0)),DIMENSION(ReadLinesOrig_loc*Nper_loc):: Fast  !Array to receive disaggregated data  
  REAL(KIND(1d0)),DIMENSION(ReadLinesOrig_loc):: Slow   !Array to disaggregate
  INTEGER,DIMENSION(:),ALLOCATABLE:: Subintervals  !Array of subintervals that contain rain
  INTEGER,DIMENSION(Nper_loc):: seq1Nper_loc   !1 to Nper_loc
  INTEGER:: i

  ! For each averaging period, get subintervals which will receive rain
  ALLOCATE(Subintervals(amongN))
  Subintervals(:) = -999

  seq1Nper_loc = (/(i, i=1,Nper_loc, 1)/)

  IF(DiagnoseDisagg==1) write(*,*) 'Distributing over ',amongN,' subintervals for variable'
  ! If all subintervals are to contain rain, don't need to generate random numbers
  IF(amongN == Nper_loc) THEN
     Subintervals(:) = seq1Nper_loc
  ENDIF
  IF(amongN > Nper_loc) CALL errorHint(2,'Problem in SUEWS_MetDisagg: no. of rainy periods cannot exceed number of subintervals', &
                                            REAL(Nper_loc,KIND(1d0)),NotUsed,amongN)


  ! Initialise fast array to -999
  Fast = -999
  DO i=1,ReadLinesOrigMax_loc 
     Fast(Nper_loc*(i-1)+seq1Nper_loc) = 0   !Fill all subintervals with zeros initially
     IF(Slow(i) > 0) THEN   !If there is some rainfall during this interval...
        IF(amongN < Nper_loc) THEN
           Subintervals(:) = -999   
           Subintervals = RandomSamples(amongN,Nper_loc)
        ENDIF
        Fast(Nper_loc*(i-1)+SubIntervals) = Slow(i)/amongN
     ENDIF
  ENDDO

  IF(ANY(Fast(1:ReadLinesOrigMax_loc*Nper_loc) == -999)) THEN
     write(*,*) 'Problem: -999s (',COUNT(Fast(1:ReadLinesOrigMax_loc*Nper_loc) == -999),') in disaggregated data'
     CALL errorHint(13,'Problem in SUEWS_MetDisagg: -999 values in disaggregated data.',NotUsed,NotUsed,NotUsedI)
  ENDIF

ENDFUNCTION DisaggP_amongN
!======================================================================================

!======================================================================================
FUNCTION RandomSamples(N,OutOf) RESULT(Samples)
! Generates N/OutOf random samples without repeats
!   e.g. for N = 3 and OutOf = 12, a possibility for Samples = 7,3,11   
! HCW 10 Feb 2017    
!======================================================================================

IMPLICIT NONE   

INTEGER:: i   !counter
INTEGER:: N   !number of samples to return
INTEGER:: OutOf   !number to sample from
INTEGER:: X   !next sample to be added
REAL(KIND(1D0)):: r   !random number
INTEGER,DIMENSION(:),ALLOCATABLE:: Samples   !Array to receive random samples

! Allocate and initialise Samples
ALLOCATE(Samples(N))
Samples(:) = -999

! Generate random sample (no repeats)
i=0 !Set counter to zero initially
DO WHILE (any(Samples == -999))
CALL random_number(r)
X = int(r*OutOf)+1
!write(*,*) X
!write(*,*) COUNT(Samples == X)
IF(COUNT(Samples == X) == 0) THEN
  ! Only keep if this subinterval has not already been selected
  i=i+1
  Samples(i)=X
ENDIF   
!write(*,*) Samples
ENDDO       

ENDFUNCTION RandomSamples
!======================================================================================



ENDMODULE MetDisagg
!========================================================================================