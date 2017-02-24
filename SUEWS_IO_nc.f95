!========================================================================================
! netCDF conversion subroutines for SUEWS
! author: Ting Sun
!
! disclamier:
!     This code employs the netCDF Fortran 90 API.
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
!     Part of the work is under the help of examples provided by the documentation.
!
! purpose:
! these subroutines write out the results of SUEWS in netCDF format.
!
!
! history:
! 20161209: initial version
! 20161213: standalise the txt2nc procedure
!========================================================================================


SUBROUTINE SiteSelect_txt2nc
  !
  USE allocateArray
  USE ColNamesInputFiles
  USE data_in
  USE defaultNotUsed
  USE FileName
  USE initial
  USE gis_data
  USE mod_z
  USE resist
  USE snowMod
  USE sues_data
  USE time
  USE netcdf

  IMPLICIT NONE

  CHARACTER(len=100):: rawpath

  ! We are writing 2D data, {y, x}
  INTEGER, PARAMETER :: NDIMS = 2, iVarStart=1
  INTEGER :: NX,NY
  !
  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  INTEGER :: ncID, varID, dimids(NDIMS),varIDx,varIDy
  INTEGER :: x_dimid, y_dimid,iVar
  REAL(KIND(1d0)), ALLOCATABLE :: varOut(:,:), varSeq(:),varSeq0(:),&
       varX(:,:),varY(:,:),xLat(:,:),xLon(:,:)
  INTEGER :: idVar(iVarStart:ncolumnsSiteSelect)
  CHARACTER(len=25):: nameVarList(iVarStart:ncolumnsSiteSelect),ivarStr2

  ! variable names:
  nameVarList=(/&
       "Grid                     ",  "Year                     ",  "StartDLS                 ",&
       "EndDLS                   ",  "lat                      ",  "lng                      ",&
       "Timezone                 ",  "SurfaceArea              ",  "Alt                      ",&
       "z                        ",  "id                       ",  "ih                       ",&
       "imin                     ",  "Fr_Paved                 ",  "Fr_Bldgs                 ",&
       "Fr_EveTr                 ",  "Fr_DecTr                 ",  "Fr_Grass                 ",&
       "Fr_Bsoil                 ",  "Fr_Water                 ",  "IrrFr_EveTr              ",&
       "IrrFr_DecTr              ",  "IrrFr_Grass              ",  "H_Bldgs                  ",&
       "H_EveTr                  ",  "H_DecTr                  ",  "z0                       ",&
       "zd                       ",  "FAI_Bldgs                ",  "FAI_EveTr                ",&
       "FAI_DecTr                ",  "PopDensDay               ",  "PopDensNight             ",&
       "TrafficRate              ",  "BuildEnergyUse           ",  "Code_Paved               ",&
       "Code_Bldgs               ",  "Code_EveTr               ",  "Code_DecTr               ",&
       "Code_Grass               ",  "Code_Bsoil               ",  "Code_Water               ",&
       "LUMPS_DrRate             ",  "LUMPS_Cover              ",  "LUMPS_MaxRes             ",&
       "NARP_Trans               ",  "CondCode                 ",  "SnowCode                 ",&
       "SnowClearingProfWD       ",  "SnowClearingProfWE       ",  "AnthropogenicCode        ",&
       "EnergyUseProfWD          ",  "EnergyUseProfWE          ",  "ActivityProfWD           ",&
       "ActivityProfWE           ",  "IrrigationCode           ",  "WaterUseProfManuWD       ",&
       "WaterUseProfManuWE       ",  "WaterUseProfAutoWD       ",  "WaterUseProfAutoWE       ",&
       "FlowChange               ",  "RunoffToWater            ",  "PipeCapacity             ",&
       "GridConnection1of8       ",  "Fraction1of8             ",  "GridConnection2of8       ",&
       "Fraction2of8             ",  "GridConnection3of8       ",  "Fraction3of8             ",&
       "GridConnection4of8       ",  "Fraction4of8             ",  "GridConnection5of8       ",&
       "Fraction5of8             ",  "GridConnection6of8       ",  "Fraction6of8             ",&
       "GridConnection7of8       ",  "Fraction7of8             ",  "GridConnection8of8       ",&
       "Fraction8of8             ",  "WithinGridPavedCode      ",  "WithinGridBldgsCode      ",&
       "WithinGridEveTrCode      ",  "WithinGridDecTrCode      ",  "WithinGridGrassCode      ",&
       "WithinGridUnmanBSoilCode ",  "WithinGridWaterCode      ",  "AreaWall                 ",&
       "Fr_ESTMClass_Paved1      ",  "Fr_ESTMClass_Paved2      ",  "Fr_ESTMClass_Paved3      ",&
       "Code_ESTMClass_Paved1    ",  "Code_ESTMClass_Paved2    ",  "Code_ESTMClass_Paved3    ",&
       "Fr_ESTMClass_Bldgs1      ",  "Fr_ESTMClass_Bldgs2      ",  "Fr_ESTMClass_Bldgs3      ",&
       "Fr_ESTMClass_Bldgs4      ",  "Fr_ESTMClass_Bldgs5      ",  "Code_ESTMClass_Bldgs1    ",&
       "Code_ESTMClass_Bldgs2    ",  "Code_ESTMClass_Bldgs3    ",  "Code_ESTMClass_Bldgs4    ",&
       "Code_ESTMClass_Bldgs5    " /)

  !=======================convert txt to netCDF============================
  ! file names
  rawpath=TRIM(FileInputPath)
  FileOut=TRIM(rawpath)//'SUEWS_SiteSelect.nc'
  ! PRINT*, 'GOOD 1'

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  CALL check( nf90_create(FileOut, NF90_CLOBBER, ncID) )
  ! print*, FileOut
  ! PRINT*, 'GOOD 2'

  ! define the dimension of spatial array/frame in the output
  nX=nCol
  nY=nRow

  ! allocate arrays
  ALLOCATE(varOut(nX,nY))
  ALLOCATE(varSeq0(nX*nY))
  ALLOCATE(varSeq(nX*nY))
  ALLOCATE(xLon(nX,nY))
  ALLOCATE(xLat(nX,nY))
  ALLOCATE(varY(nX,nY))
  ALLOCATE(varX(nX,nY))

  ! latitude:
  varSeq0=SiteSelect(1:nX*nY,5)
  CALL sortSeqReal(varSeq0,varSeq,nY,nX)
  xLat = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
  ! PRINT*, 'before flipping:',xLat(1:5,1)
  xLat =xLat(:,nY:1:-1)
  ! PRINT*, 'after flipping:',xLat(1:5,1)

  ! longitude:
  varSeq0=SiteSelect(1:nX*nY,6)
  CALL sortSeqReal(varSeq0,varSeq,nY,nX)
  xLon = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )


  ! pass values to coordinate variables
  varY = xLat
  varX = xLon


  ! Define the dimensions. NetCDF will hand back an ID for each.
  ! nY = ncolumnsDataOut-4
  ! nx = NumberOfGrids
  ! CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )
  CALL check( nf90_def_dim(ncID, "west_east", NX, x_dimid) )
  CALL check( nf90_def_dim(ncID, "south_north", NY, y_dimid) )
  ! PRINT*, 'GOOD 3'

  dimids =  (/x_dimid, y_dimid/)

  ! define 2D variables:

  ! define coordinate variables:
  CALL check( nf90_def_var(ncID,'xLon', NF90_REAL, (/x_dimid, y_dimid/), varIDx))
  CALL check( nf90_put_att(ncID,varIDx,'units','degree_east') )

  CALL check( nf90_def_var(ncID,'xLat', NF90_REAL, (/x_dimid, y_dimid/), varIDy))
  CALL check( nf90_put_att(ncID,varIDy,'units','degree_north') )

  ! define variables in SiteSelect
  DO iVar = iVarStart, ncolumnsSiteSelect, 1
     ! define variable name
     ivarStr2=nameVarList(iVar)

     ! Define the variable. The type of the variable in this case is
     ! NF90_REAL.
     CALL check( nf90_def_var(ncID,TRIM(ivarStr2), NF90_REAL, dimids, varID) )
     CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
     !  print*, 'put att good'
     !  CALL check( nf90_put_att(ncID,varID,'units',TRIM(ADJUSTL(iunitStr2))) )
     !  print*, 'put unit good'
     idVar(iVar)=varID
  END DO
  CALL check( nf90_enddef(ncID) )
  ! End define mode. This tells netCDF we are done defining metadata.

  ! put coordinate variables:
  CALL check( nf90_put_var(ncID, varIDx, varX) )
  CALL check( nf90_put_var(ncID, varIDy, varY) )
  CALL check( NF90_SYNC(ncID) )

  ! put 2D variables in
  DO iVar = iVarStart, ncolumnsSiteSelect, 1
     !  PRINT*, 'dim1:', SIZE(dataOut(1:irMax,iVar,:), dim=1)
     !  PRINT*, 'dim2:',SIZE(dataOut(1:irMax,iVar,:), dim=2)
     ! reshape dataOut to be aligned in checker board form
     varSeq0=SiteSelect(1:nX*nY,iVar)
     CALL sortSeqReal(varSeq0,varSeq,nY,nX)
     varOut = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
     varOut = varOut(:,nY:1:-1)
     !  get the variable id
     varID= idVar(iVar)
     !  PRINT*, 'GOOD 5',iVar

     CALL check( nf90_put_var(ncID, varID, varOut) )
     CALL check(NF90_SYNC(ncID))
  END DO
  IF (ALLOCATED(varOut)) DEALLOCATE(varOut)
  IF (ALLOCATED(varSeq0)) DEALLOCATE(varSeq0)
  IF (ALLOCATED(varSeq)) DEALLOCATE(varSeq)
  IF (ALLOCATED(xLon)) DEALLOCATE(xLon)
  IF (ALLOCATED(xLat)) DEALLOCATE(xLat)
  IF (ALLOCATED(varY)) DEALLOCATE(varY)
  IF (ALLOCATED(varX)) DEALLOCATE(varX)

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  CALL check( nf90_close(ncID) )

END SUBROUTINE SiteSelect_txt2nc

! SUBROUTINE InitialConditions_txt2nc
!   !
!   USE allocateArray
!   USE ColNamesInputFiles
!   USE data_in
!   USE defaultNotUsed
!   USE FileName
!   USE initial
!   USE gis_data
!   USE mod_z
!   USE resist
!   USE snowMod
!   USE sues_data
!   USE time
!   USE netcdf
!
!   IMPLICIT NONE
!
!
!   !
!   CHARACTER(len=100):: rawpath
!   !
!   ! We are writing 2D data, {y, x}
!   INTEGER, PARAMETER :: NDIMS = 2, iVarStart=1
!   INTEGER :: NX,NY
!   !
!   ! When we create netCDF files, variables and dimensions, we get back
!   ! an ID for each one.
!   INTEGER :: ncID, varID, dimids(NDIMS)
!   INTEGER :: x_dimid, y_dimid,iVar
!   REAL(KIND(1d0)), ALLOCATABLE :: varOut(:,:), varSeq(:),varSeq0(:)
!   INTEGER :: idVar(iVarStart:ncolumnsSiteSelect)
!   CHARACTER(len=25):: nameVarList(iVarStart:ncolumnsSiteSelect),ivarStr2
!
!   ! variable names:
!   nameVarList=(/&
!
!        /)
!
!   !=======================convert txt to netCDF============================
!   ! file names
!   rawpath=TRIM(FileInputPath)
!   FileOut=TRIM(rawpath)//'SUEWS_InitialCondition.nc'
!   ! PRINT*, 'GOOD 1'
!
!   ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
!   ! overwrite this file, if it already exists.
!   CALL check( nf90_create(FileOut, NF90_CLOBBER, ncID) )
!   ! print*, FileOut
!   ! PRINT*, 'GOOD 2'
!
!
!   nX=nCol
!   nY=nRow
!   ! write out each variable
!   ALLOCATE(varOut(nX,nY))
!   ALLOCATE(varSeq(nlinesSiteSelect))
!   ALLOCATE(varSeq0(nlinesSiteSelect))
!
!
!   ! Define the dimensions. NetCDF will hand back an ID for each.
!   ! nY = ncolumnsDataOut-4
!   ! nx = NumberOfGrids
!   ! CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )
!   CALL check( nf90_def_dim(ncID, "west_east", NX, x_dimid) )
!   CALL check( nf90_def_dim(ncID, "south_north", NY, y_dimid) )
!   ! PRINT*, 'GOOD 3'
!
!   dimids =  (/x_dimid, y_dimid/)
!
!   ! define 2D variables:
!   DO iVar = iVarStart, ncolumnsSiteSelect, 1
!      ! define variable name
!      ivarStr2=nameVarList(iVar)
!
!      ! Define the variable. The type of the variable in this case is
!      ! NF90_REAL.
!      CALL check( nf90_def_var(ncID,TRIM(ivarStr2), NF90_REAL, dimids, varID) )
!      idVar(iVar)=varID
!   END DO
!   CALL check( nf90_enddef(ncID) )
!   ! End define mode. This tells netCDF we are done defining metadata.
!
!   ! then other 3D variables
!   DO iVar = iVarStart, ncolumnsSiteSelect, 1
!      !  PRINT*, 'dim1:', SIZE(dataOut(1:irMax,iVar,:), dim=1)
!      !  PRINT*, 'dim2:',SIZE(dataOut(1:irMax,iVar,:), dim=2)
!      ! reshape dataOut to be aligned in checker board form
!      varSeq0=SiteSelect(1:nX*nY,iVar)
!      CALL sortSeqReal(varSeq0,varSeq,nY,nX)
!      varOut = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
!      varOut = varOut(:,nY:1:-1)
!      !  get the variable id
!      varID= idVar(iVar)
!      !  PRINT*, 'GOOD 5',iVar
!
!      CALL check( nf90_put_var(ncID, varID, varOut) )
!      CALL check(NF90_SYNC(ncID))
!   END DO
!   IF (ALLOCATED(varOut)) DEALLOCATE(varOut)
!
!   ! Close the file. This frees up any internal netCDF resources
!   ! associated with the file, and flushes any buffers.
!   CALL check( nf90_close(ncID) )
!
! END SUBROUTINE InitialConditions_txt2nc



!===========================================================================!
! write the output of final SUEWS results in netCDF
!   with spatial layout of QGIS convention
! the spatial matrix arranges successive rows down the page (i.e., north to south)
!   and succesive columns across (i.e., west to east)
! the output file frequency is the same as metblocks in the main SUEWS loop
!===========================================================================!
SUBROUTINE SUEWS_Output_nc(year_int,iblock,irMax)
  !INPUT: year_int = Year as a integer
  !       iblock = Block number of met data
  !       irMax = Maximum number of rows in met data
  USE allocateArray
  USE cbl_module
  USE data_in
  USE defaultNotUsed
  USE ESTM_data
  USE gis_data
  USE initial
  USE solweig_module
  USE sues_data
  USE time
  USE strings
  ! USE SetupOutput
  USE netCDF


  IMPLICIT NONE


  INTEGER:: year_int, iblock, irMax, CurrentGrid !inputs
  INTEGER:: i, j, nlinesOut,nsize
  REAL(KIND(1d0)),ALLOCATABLE:: dataOutProc0(:,:,:),dataOutProc(:,:),dataOutProcX(:,:,:)

  CHARACTER(len=10):: str2, str2_tt, grstr2, yrstr2,iblockStr2
  CHARACTER(len=100):: rawpath, SnowOut,ESTMOut, FileOutFormat


  ! N.B. if change lengths here, also adjust in MODULE AllocateArray accordingly
  CHARACTER(len=10),DIMENSION(nColumnsDataOut):: HeaderAll, FormatAll  !Header and formats for all output variables
  CHARACTER(len=14*nColumnsDataOut):: HeaderOut, FormatOut             !Header and format for selected output variables (untrimmed)
  CHARACTER(len=12*nColumnsDataOut):: FormatOutNoSep, HeaderOutNoSep   !Format for selected output variables (not comma sep)
  CHARACTER(len=12),DIMENSION(nColumnsDataOut):: UnitsAll              !Units for all output variables
  CHARACTER(len=14*nColumnsDataOut):: UnitsOut                         !Units for selected output variables
  CHARACTER(len=50),DIMENSION(nColumnsDataOut):: LongNmAll             !LongName for all output variables
  CHARACTER(len=52*nColumnsDataOut):: LongNmOut                        !LongName for selected output variables (untrimmed)
  CHARACTER(len= 1),DIMENSION(nColumnsDataOut):: AggregAll             !Aggregation method required for all output variables
  CHARACTER(len= 3*nColumnsDataOut):: AggregOut                        !Aggregation method required for selected output variables
  CHARACTER(len= 4*nColumnsDataOut):: ColNos
  CHARACTER(len= 1),ALLOCATABLE:: AggregUseX(:)                        !Aggregation method array required for selected output variables

  CHARACTER(len=10):: fy, ft, fd, f94, f104, f106   !Useful formats
  CHARACTER(len= 1):: aT, aA, aS, aL   !Useful formats
  CHARACTER(len= 3):: itext

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
  WRITE(str2_tt,'(i4)') tstep/60
  WRITE(str2,'(i4)') ResolutionFilesOut/60
  WRITE(grstr2,'(i10)') CurrentGrid
  WRITE(yrstr2,'(i4)') year_int
  WRITE(iblockStr2,'(i4)') iblock

  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//'_'//TRIM(ADJUSTL(yrstr2)) ! output resolution added, TS 9 Feb 2017
  ! For files at specified output resolution
  FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(iblockStr2))//'_'//TRIM(ADJUSTL(str2))//'.nc'
  SOLWEIGpoiOut=TRIM(rawpath)//'_SOLWEIGpoiOut.nc'
  ESTMOut=TRIM(rawpath)//'_ESTM_'//TRIM(ADJUSTL(str2))//'.nc' ! output resolution added, TS 10 Feb 2017
  BLOut=TRIM(rawpath)//'_BL.nc'
  SnowOut=TRIM(rawpath)//'_snow_5.nc'
  ! keep the OutputFormat file in plain text
  FileOutFormat=TRIM(FileOutputPath)//TRIM(FileCode)//'_YYYY_'//TRIM(ADJUSTL(str2))//'_OutputFormat.txt'

  ! For main data file at model time-step (may not be used if KeepTstepFilesOut is 0)
  FileOut_tt=TRIM(rawpath)//'_'//TRIM(ADJUSTL(iblockStr2))//'_'//TRIM(ADJUSTL(str2_tt))//'.nc'

  !========== Get headers and write out output info to file ==========
  ! To add extra columns, change all these (Header, Units, LongNm, Format, Agg) together
  ! Could change to read from external file later
  IF(OutputFormats==1) THEN   !Once per run
     ! Set all output variables here. This must agree with dataOut (see SUEWS_Calculations.f95)
     HeaderAll(:) = '-'   !Initialise
     UnitsAll (:) = '-'
     FormatAll(:) = '-'
     AggregAll(:) = '-'
     LongNmAll(:) = '-'

     HeaderAll(1:5) = (/'   Year','    DOY','   Hour','    Min','Dectime'/)   !datetime info
     UnitsAll (1:5) = (/'   YYYY','    DOY','     HH','     MM','      -'/)
     FormatAll(1:5) = (/fy,ft,ft,ft,fd/)
     AggregAll(1:5) = aT
     LongNmAll(1:5) = (/'        Year',' Day of Year','        Hour','      Minute','Decimal time'/)

     HeaderAll(6:9) = (/'Kdown','  Kup','Ldown','  Lup'/)  !radiation components
     UnitsAll (6:9) = 'W_m-2'
     FormatAll(6:9) = f94
     AggregAll(6:9) = aA
     LongNmAll(6:9) = (/'Incoming shortwave radiation','Outgoing shortwave radiation', &
          ' Incoming longwave radiation',' Outgoing longwave radiation'/)

     HeaderAll(10) = 'Tsurf'
     UnitsAll (10) = 'degC'
     FormatAll(10) = f94
     AggregAll(10) = aA
     LongNmAll(10) = 'Bulk surface temperature'

     HeaderAll(11:15) = (/'QN','QF','QS','QH','QE'/)  !energy fluxes
     UnitsAll (11:15) = 'W_m-2'
     FormatAll(11:15) = f94
     AggregAll(11:15) = aA
     LongNmAll(11:15) = (/' Net all-wave radiation','Anthropogenic heat flux','  Net storage heat flux', &
          '     Sensible heat flux','       Latent heat flux'/)

     HeaderAll(16:18) = (/'QHlumps','QElumps','QHresis'/)   !energy fluxes (other approaches)
     UnitsAll (16:18)  = 'W_m-2'
     FormatAll(16:18) = f94
     AggregAll(16:18) = aA
     LongNmAll(16:18) = (/ '      Sensible heat flux (using LUMPS)','        Latent heat flux (using LUMPS)', &
          'Sensible heat flux (resistance method)'/)

     HeaderAll(19:23) = (/' Rain','  Irr',' Evap','   RO','TotCh'/)   !water balance components
     UnitsAll (19:23) = 'mm'
     FormatAll(19:23) = f106
     AggregAll(19:23) = aS
     LongNmAll(19:23) = (/'                            Rain','                      Irrigation', &
          '                     Evaporation','                          Runoff', &
          'Surface and soil moisture change'/)

     HeaderAll(24:28) = (/'    SurfCh','     State',' NWtrState','  Drainage','       SMD'/)   !water balance components cont.
     UnitsAll (24:28) = 'mm'
     FormatAll(24:28) = (/ f106,f104,f106,f106,f94 /)
     AggregAll(24:28) = (/aS,aL,aL,aS,aL/)
     LongNmAll(24:28) = (/'                   Surface moisture change','                     Surface wetness state', &
          'Surface wetness state (non-water surfaces)','                                  Drainage', &
          '                     Soil moisture deficit'/)


     HeaderAll(29:30) = (/'  FlowCh','AddWater'/)                                          !water balance components cont.
     UnitsAll (29:30) = 'mm'
     FormatAll(29:30) = f104
     AggregAll(29:30) = aS
     LongNmAll(29:30) = (/' Additional flow into water body','Addtional water from other grids'/)

     HeaderAll(31:35) = (/' ROSoil',' ROPipe','  ROImp','  ROVeg','ROWater'/)   !runoff components
     UnitsAll (31:35) = 'mm'
     FormatAll(31:35) = f106
     AggregAll(31:35) = aS
     LongNmAll(31:35) = (/'                 Runoff to soil','                Runoff to pipes', &
          'Runoff over impervious surfaces',' Runoff over vegetated surfaces', &
          '       Runoff for water surface'/)

     HeaderAll(36:39) = (/'  WUInt','WUEveTr','WUDecTr','WUGrass'/)                !water use
     UnitsAll (36:39) = 'mm'
     FormatAll(36:39) = f94
     AggregAll(36:39) = aS
     LongNmAll(36:39) = (/'             InternalWaterUse', &
          'Water use for evergreen trees','Water use for deciduous trees','          Water use for grass'/)

     HeaderAll(40:45) = (/'SMDPaved','SMDBldgs','SMDEveTr','SMDDecTr','SMDGrass','SMDBSoil'/)  !smd for each surface
     UnitsAll (40:45) = 'mm'
     FormatAll(40:45) = f94
     AggregAll(40:45) = aL
     LongNmAll(40:45) = (/'         Soil moisture deficit for paved surface', &
          '      Soil moisture deficit for building surface', &
          'Soil moisture deficit for evergreen tree surface', &
          'Soil moisture deficit for deciduous tree surface', &
          '         Soil moisture deficit for grass surface', &
          '     Soil moisture deficit for bare soil surface'/)

     HeaderAll(46:52) = (/'StPaved','StBldgs','StEveTr','StDecTr','StGrass','StBSoil','StWater'/)   !states
     UnitsAll (46:52) = 'mm'
     FormatAll(46:52) = (/SPREAD(f94,1,6),f104/)
     AggregAll(46:52) = aL
     LongNmAll(46:52) = (/'         Surface wetness state for paved surface', &
          '      Surface wetness state for building surface', &
          'Surface wetness state for evergreen tree surface', &
          'Surface wetness state for deciduous tree surface', &
          '         Surface wetness state for grass surface', &
          '     Surface wetness state for bare soil surface', &
          '         Surface wetness state for water surface'/)

     HeaderAll(53:54) = (/' Zenith','Azimuth'/)! solar angles
     UnitsAll (53:54) = 'deg'
     FormatAll(53:54) = f94
     AggregAll(53:54) = aL
     LongNmAll(53:54) = (/' Solar zenith angle','Solar azimuth angle'/)

     HeaderAll(55:56) = (/'AlbBulk','   Fcld'/)                ! extra radiation info
     UnitsAll (55:56) = '-'
     FormatAll(55:56) = f94
     AggregAll(55:56) = aA
     LongNmAll(55:56) = (/'   Bulk albedo','Cloud fraction'/)

     HeaderAll(57)    = 'LAI'   ! extra surface info
     UnitsAll (57)    = 'm2_m-2'
     FormatAll(57)    = f94
     AggregAll(57)    = aA
     LongNmAll(57)    = 'Leaf area index'

     HeaderAll(58:59) = (/'z0m','zdm'/)
     UnitsAll (58:59) = 'm'
     FormatAll(58:59) = f94
     AggregAll(58:59) = aA
     LongNmAll(58:59) = (/' Roughness length for momentum','Zero-plane displacement height'/)

     HeaderAll(60:63) = (/'ustar','  Lob','   ra','   rs'/)                ! turbulence
     UnitsAll (60:63) = (/'m_s-1','    m','s_m-1','s_m-1'/)
     FormatAll(60:63) = (/f94,f104,f94,f94/)
     AggregAll(60:63) = aA
     LongNmAll(60:63) = (/'     Friction velocity','        Obukhov length', &
          'Aerodynamic resistance','    Surface resistance'/)

     HeaderAll(64:69) = (/'     Fc','FcPhoto','FcRespi','FcMetab','FcTraff','FcBuild'/)   ! CO2 flux & components
     UnitsAll (64:69) = 'umol_m-2_s-1'
     FormatAll(64:69) = f94
     AggregAll(64:69) = aA
     LongNmAll(64:69) = (/'                    CO2 flux', &
          'CO2 flux from photosynthesis','   CO2 flux from respiration', &
          '    CO2 flux from metabolism','       CO2 flux from traffic','     CO2 flux from buildings'/)

     HeaderAll(70:72) = (/'QNSnowFr','  QNSnow',' AlbSnow'/)                             ! snow-related (radiation)
     UnitsAll (70:72) = (/'W_m-2','W_m-2','    -'/)
     FormatAll(70:72) = f94
     AggregAll(70:72) = aA
     LongNmAll(70:72) = (/'Net all-wave radiation for non-snow area','    Net all-wave radiation for snow area', &
          '                             Snow albedo'/)

     HeaderAll(73:75) = (/'        QM','  QMFreeze','    QMRain'/)
     UnitsAll (73:75) = 'W_m-2'
     FormatAll(73:75) = f106
     AggregAll(73:75) = aA
     LongNmAll(73:75) = (/'   Snow-related heat exchange','       Internal energy change','Heat released by rain on snow'/)

     HeaderAll(76:79) = (/'       SWE',' MeltWater','MeltWStore','    SnowCh'/)   !snow
     UnitsAll (76:79) = 'mm'
     FormatAll(76:79) = f106
     AggregAll(76:79) = aS
     LongNmAll(76:79) = (/'Snow water equivalent','            Meltwater','      Meltwater store','  Change in snow pack' /)

     HeaderAll(80:81) = (/'SnowRPaved','SnowRBldgs'/)   !snow-related (removal)
     UnitsAll (80:81) = 'mm'
     FormatAll(80:81) = f94
     AggregAll(80:81) = aS
     LongNmAll(80:81) = (/'   Snow removed from paved surface','Snow removed from building surface' /)


     ! Select variables to be written out
     !write(*,*) 'WriteOutOption:', WriteOutOption
     IF(WriteOutOption == 0) THEN   !all (not snow-related)
        ALLOCATE(UseColumnsDataOut(69))
        UseColumnsDataOut = (/ (i, i=1,69, 1) /)
     ELSEIF(WriteOutOption == 1) THEN   !all plus snow-related
        ALLOCATE(UseColumnsDataOut(nColumnsDataOut))
        UseColumnsDataOut = (/ (i, i=1,nColumnsDataOut, 1) /)
     ELSEIF(WriteOutOption == 2) THEN   !minimal output
        ALLOCATE(UseColumnsDataOut(33))
        UseColumnsDataOut = (/ (i, i=1,15, 1),(i, i=19,28, 1), 53,54,55,56, 57, 60,61, 64 /)
     ELSE
        WRITE(*,*) 'RunControl: WriteOutOption code not recognised, so writing out all variables.'
        ALLOCATE(UseColumnsDataOut(69))
        UseColumnsDataOut = (/ (i, i=1,69, 1) /)
     ENDIF

     ! Create subset of HeaderAll and FormatAll for selected variables only
     DO i=1,SIZE(UseColumnsDataOut)
        WRITE(itext,'(i3)') i
        IF(i==1) THEN
           HeaderOut=ADJUSTL(HeaderAll(UseColumnsDataOut(i)))
           HeaderOutNoSep=ADJUSTL(HeaderAll(UseColumnsDataOut(i)))
           UnitsOut=ADJUSTL(UnitsAll(UseColumnsDataOut(i)))
           FormatOut=ADJUSTL(FormatAll(UseColumnsDataOut(i)))
           FormatOutNoSep=ADJUSTL(FormatAll(UseColumnsDataOut(i)))
           LongNmOut=ADJUSTL(LongNmAll(UseColumnsDataOut(i)))
           AggregOut=ADJUSTL(AggregAll(UseColumnsDataOut(i)))
           ColNos=ADJUSTL(itext)
        ELSE
           HeaderOut=TRIM(HeaderOut)//';'//ADJUSTL(HeaderAll(UseColumnsDataOut(i)))
           HeaderOutNoSep=TRIM(HeaderOutNoSep)//' '//ADJUSTL(HeaderAll(UseColumnsDataOut(i)))
           !write(*,*) HeaderOut
           UnitsOut=TRIM(UnitsOut)//';'//ADJUSTL(UnitsAll(UseColumnsDataOut(i)))
           !write(*,*) UnitsOut
           FormatOut=TRIM(FormatOut)//';'//ADJUSTL(FormatAll(UseColumnsDataOut(i)))
           FormatOutNoSep=TRIM(FormatOutNoSep)//' '//ADJUSTL(FormatAll(UseColumnsDataOut(i)))
           !write(*,*) FormatOut
           LongNmOut=TRIM(LongNmOut)//';'//ADJUSTL(LongNmAll(UseColumnsDataOut(i)))
           !write(*,*) LongNmOut
           AggregOut=TRIM(AggregOut)//';'//ADJUSTL(AggregAll(UseColumnsDataOut(i)))
           !write(*,*) AggregOut
           ColNos=TRIM(ColNos)//';'//ADJUSTL(itext)
        ENDIF
     ENDDO
     !HeaderUse=trim(adjustl(HeaderOut))//' ' !with extra space at end of header row
     !  PRINT*, 'mem start'

     !ALLOCATE(CHARACTER(LEN(trim(adjustl(HeaderOut)))):: HeaderUse)
     !  PRINT*, 'mem',1

     !ALLOCATE(CHARACTER(LEN(trim(adjustl(UnitsOut)))):: UnitsUse)
     !ALLOCATE(CHARACTER(LEN(trim(adjustl(LongNmOut)))):: LongNmUse)

     !  PRINT*, LEN(TRIM(ADJUSTL(FormatOut)))
     !  PRINT*, FormatOut
     !ALLOCATE(CHARACTER(LEN(trim(adjustl(FormatOut)))):: FormatUse)
     !  PRINT*, 'mem',3

     !  PRINT*, LEN(TRIM(ADJUSTL(AggregOut)))
     !  PRINT*, AggregOut
     !ALLOCATE(CHARACTER(LEN(trim(adjustl(AggregOut)))):: AggregUse)
     !  PRINT*, 'mem',4

     !ALLOCATE(CHARACTER(LEN(trim(adjustl(ColNos)))):: ColNosUse)
     HeaderUse=TRIM(ADJUSTL(HeaderOut))

     HeaderUseNoSep=TRIM(ADJUSTL(HeaderOutNoSep))
     UnitsUse=TRIM(ADJUSTL(UnitsOut))
     LongNmUse=TRIM(ADJUSTL(LongNmOut))
     FormatUse=TRIM(ADJUSTL(FormatOut))
     FormatUseNoSep='('//TRIM(ADJUSTL(FormatOutNoSep))//')'
     AggregUse=TRIM(ADJUSTL(AggregOut))
     ColNosUse=TRIM(ADJUSTL(ColNos))
     !write(*,*) '||',TRIM(HeaderUse),'||'
     !write(*,*) '||',TRIM(FormatUse),'||'
     !write(*,*) '||',TRIM(ColNosUse),'||'


     !=========== Write output format info to file ===========
     OPEN(50,file=TRIM(FileOutFormat),err=111)
     WRITE(50,'(a)') TRIM(ColNosUse)
     WRITE(50,'(a)') TRIM(HeaderUse)
     WRITE(50,'(a)') TRIM(LongNmUse)
     WRITE(50,'(a)') TRIM(UnitsUse)
     WRITE(50,'(a)') TRIM(FormatUse)   !also write formats to output file (without outer brackets)
     WRITE(50,'(a)') TRIM(AggregUse)
     CLOSE (50)
     OutputFormats = 0
  ENDIF

  ALLOCATE(AggregUseX(SIZE(UseColumnsDataOut)))
  CALL parse(AggregUse,';',AggregUseX,SIZE(UseColumnsDataOut))
  ! PRINT*, 'good 1'
  ! PRINT*, AggregUseX


  !========== Write out data ==========
  IF ( ResolutionFilesOut == Tstep .OR. KeepTstepFilesOut == 1) THEN ! output frequency same as input, or specify to keep raw output files (HCW 20 Feb 2017)
     ! original output

     ! write out processed/aggregated data
     CALL SUEWS_Write_nc(TRIM(FileOut_tt),HeaderUse,LongNmUse,UnitsUse,&
          dataOut(1:irMax,1:SIZE(UseColumnsDataOut),1:NumberOfGrids),irMax)

     !  IF (SOLWEIGpoi_out==1) THEN
     !     DO i=1,SolweigCount-1
     !        WRITE(9,304) INT(dataOutSOL(i,1,Gridiv)),(dataOutSOL(i,is,Gridiv),is=2,ncolumnsdataOutSOL)
     !     ENDDO
     !  ENDIF
     !
     !  IF(CBLuse>=1) THEN
     !     DO i=1,iCBLcount
     !        WRITE(53,305)(INT(dataOutBL(i,is,Gridiv)),is=1,4),(dataOutBL(i,is,Gridiv),is=5,ncolumnsdataOutBL)
     !     ENDDO
     !  ENDIF
     !
     !  IF(SnowUse>=1) THEN
     !     DO i=1,irmax
     !        WRITE(54,306)(INT(dataOutSnow(i,is,Gridiv)),is=1,4),(dataOutSnow(i,is,Gridiv),is=5,ncolumnsDataOutSnow)
     !     ENDDO
     !  ENDIF
     !
     !  IF (StorageHeatMethod==4 .OR. StorageHeatMethod==14)THEN
     !     ! write out processed/aggregated data
     !     CALL SUEWS_Write_nc(TRIM(ESTMOut),HeaderUse,LongNmUse,UnitsUse,&
     !          dataOutESTM(1:irMax,1:UseColumnsDataOut,1:NumberOfGrids),irMax)
     !
     !     DO i=1,irMax
     !        WRITE(58, 307)(INT(dataOutESTM(i,is,Gridiv)),is=1,4),(dataOutESTM(i,is,Gridiv),is=5,32)
     !     ENDDO
     !  ENDIF

  ENDIF

  IF ( ResolutionFilesOut /= Tstep ) THEN ! if output frequency different from input, TS 09 Feb 2017
     ! write out every nlinesOut, 60.*60/ResolutionFilesOut = output frequency per hour
     nlinesOut=INT(nsh/(60.*60/ResolutionFilesOut))

     ! Main output file --------------------------------------------------
     !  lfnOutC=39  !Output file code
     !  IF (iv==1) THEN
     !     OPEN(lfnOutC,file=TRIM(FileOut),err=112)
     !     WRITE(lfnOutC,'(a)') HeaderUseNoSep
     !  ELSE
     !     OPEN(lfnOutC,file=TRIM(FileOut),position='append')!,err=112)
     !  ENDIF
     nsize=0
     DO i=nlinesOut,irMax,nlinesOut
        nsize=nsize+1
     ENDDO
     ALLOCATE(dataOutProcX(nsize,SIZE(UseColumnsDataOut),NumberOfGrids))

     nsize=0
     DO i=nlinesOut,irMax,nlinesOut
        nsize=nsize+1
        ALLOCATE(dataOutProc0(nlinesOut,SIZE(UseColumnsDataOut),NumberOfGrids))
        ALLOCATE(dataOutProc(SIZE(UseColumnsDataOut),NumberOfGrids))

        dataOutProc0=dataOut(i-nlinesOut+1:i,1:SIZE(UseColumnsDataOut),1:NumberOfGrids)

        DO j = 1, SIZE(AggregUseX), 1
           ! print*, i,j
           ! aggregating different variables
           SELECT CASE (AggregUseX(j))
           CASE ('0') !time columns, aT
              dataOutProc(j,:)=dataOutProc0(nlinesOut,j,:)
           CASE ('1') !average, aA
              dataOutProc(j,:)=SUM(dataOutProc0(:,j,:))/nlinesOut
           CASE ('2') !sum, aS
              dataOutProc(j,:)=SUM(dataOutProc0(:,j,:))
           CASE ('3') !last value,aL
              dataOutProc(j,:)=dataOutProc0(nlinesOut,j,:)
           END SELECT
           IF ( Diagnose==1 .AND. i==irMax ) THEN
              PRINT*, 'raw data of ',j,':'
              PRINT*, dataOutProc0(:,j,1)
              PRINT*, 'aggregated with method: ',AggregUseX(j)
              PRINT*, dataOutProc(j,1)
              PRINT*, ''
           END IF
        END DO
        dataOutProcX(nsize,:,:)=dataOutProc(:,:)




        IF (ALLOCATED(dataOutProc0)) DEALLOCATE(dataOutProc0)
        IF (ALLOCATED(dataOutProc)) DEALLOCATE(dataOutProc)

     ENDDO
     ! write out processed/aggregated data
     CALL SUEWS_Write_nc(fileOut,HeaderUse,LongNmUse,UnitsUse,dataOutProcX,nsize)
     IF (ALLOCATED(dataOutProcX)) DEALLOCATE(dataOutProcX)


  ENDIF


  IF (ALLOCATED(AggregUseX)) DEALLOCATE(AggregUseX)

  !================CLOSE OUTPUTFILE================
  RETURN

  !Error commands
  ! 110 CALL ErrorHint(52,TRIM(fileOut_tt),notUsed,notUsed,notUsedI)
111 CALL ErrorHint(52,TRIM(fileOutFormat),notUsed,notUsed,notUsedI)
  ! 112 CALL ErrorHint(52,TRIM(fileOut),notUsed,notUsed,notUsedI)


END SUBROUTINE SUEWS_Output_nc

SUBROUTINE SUEWS_Write_nc(fileOutput,header,longNm,units,dataOutput,irMax)
  ! generic subroutine to write out data in netCDF format

  USE allocateArray
  USE cbl_module
  USE data_in
  USE defaultNotUsed
  USE ESTM_data
  USE gis_data
  USE initial
  USE solweig_module
  USE sues_data
  USE time
  USE strings
  ! USE SetupOutput
  USE netCDF


  IMPLICIT NONE

  CHARACTER(len=*):: fileOutput
  CHARACTER(len=14*SIZE(UseColumnsDataOut)):: header
  CHARACTER(len=50*SIZE(UseColumnsDataOut)):: longNm
  CHARACTER(len=14*SIZE(UseColumnsDataOut)):: units

  ! We are writing 3D data, {time, y, x}
  INTEGER, PARAMETER :: NDIMS = 3, iVarStart=6
  INTEGER :: NX,NY,irMax

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  INTEGER :: ncID, varID, dimids(NDIMS),varIDGrid
  INTEGER :: x_dimid, y_dimid,time_dimid,iVar,varIDx,varIDy
  REAL(KIND(1d0)) :: dataOutput(1:irMax,1:SIZE(UseColumnsDataOut),1:NumberOfGrids)
  REAL(KIND(1d0)), ALLOCATABLE :: varOut(:,:,:),&
       varX(:,:),varY(:,:),&
       xLat(:,:),xLon(:,:),&
       varSeq0(:),varSeq(:)
  INTEGER :: idVar(iVarStart:SIZE(UseColumnsDataOut))
  CHARACTER(len=50):: nameVarList(SIZE(UseColumnsDataOut)),ivarStr2,&
       longNmList(SIZE(UseColumnsDataOut)),ilongNmStr2,&
       unitList(SIZE(UseColumnsDataOut)),iunitStr2


  !================GET THE VARIABLE-RELAED PROPERTIES================
  CALL parse(header,';',nameVarList,SIZE(UseColumnsDataOut))
  CALL parse(longNm,';',longNmList,SIZE(UseColumnsDataOut))
  CALL parse(units,';',unitList,SIZE(UseColumnsDataOut))

  ! define the dimension of spatial array/frame in the output
  nX   = nCol
  nY   = nRow

  ALLOCATE(varSeq0(nX*nY))
  ALLOCATE(varSeq(nX*nY))
  ALLOCATE(xLon(nX,nY))
  ALLOCATE(xLat(nX,nY))
  ALLOCATE(varY(nX,nY))
  ALLOCATE(varX(nX,nY))

  ! latitude:
  varSeq0=SiteSelect(1:nX*nY,5)
  CALL sortSeqReal(varSeq0,varSeq,nY,nX)
  xLat = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
  ! PRINT*, 'before flipping:',xLat(1:2,1)
  xLat =xLat(:,nY:1:-1)
  ! PRINT*, 'after flipping:',xLat(1:2,1)

  ! longitude:
  varSeq0=SiteSelect(1:nX*nY,6)
  CALL sortSeqReal(varSeq0,varSeq,nY,nX)
  xLon = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )


  ! pass values to coordinate variables
  varY = xLat
  varX = xLon
  ! PRINT*, 'size x dim 1:',SIZE(varX, dim=1)
  ! PRINT*, 'size x dim 2:',SIZE(varX, dim=2)

  ! file names
  ! rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(yrStr2))//'_'//TRIM(ADJUSTL(iblockStr2))
  ! FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(tstepStr2))//'.nc'

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  PRINT*, 'writing file:',TRIM(fileOutput)
  CALL check( nf90_create(TRIM(fileOutput), NF90_CLOBBER, ncID) )

  ! Define the dimensions. NetCDF will hand back an ID for each.
  ! nY = ncolumnsDataOut-4
  ! nx = NumberOfGrids
  CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )
  CALL check( nf90_def_dim(ncID, "west_east", NX, x_dimid) )
  CALL check( nf90_def_dim(ncID, "south_north", NY, y_dimid) )
  ! PRINT*, 'good define dim'

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/x_dimid, y_dimid, time_dimid/)

  ! write out each variable
  ALLOCATE(varOut(nX,nY,irMax))

  ! define all variables
  ! define coordinate variables:
  CALL check( nf90_def_var(ncID,'xLon', NF90_REAL, (/x_dimid, y_dimid/), varIDx))
  CALL check( nf90_put_att(ncID,varIDx,'units','degree_east') )

  CALL check( nf90_def_var(ncID,'xLat', NF90_REAL, (/x_dimid, y_dimid/), varIDy))
  CALL check( nf90_put_att(ncID,varIDy,'units','degree_north') )

  ! PRINT*, 'good define var'

  ! define grid_ID:
  CALL check( nf90_def_var(ncID,'grid_ID', NF90_INT, (/x_dimid, y_dimid/), varID))
  CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
  varIDGrid=varID

  ! define other 3D variables:
  DO iVar = iVarStart, SIZE(UseColumnsDataOut), 1
     ! define variable name
     ivarStr2=nameVarList(iVar)
     iunitStr2=unitList(iVar)
     ilongNmStr2=longNmList(iVar)
     !  PRINT*, iunitStr2

     ! Define the variable. The type of the variable in this case is
     ! NF90_REAL.
     !  PRINT*, TRIM(ADJUSTL(ivarStr2))
     CALL check( nf90_def_var(ncID,TRIM(ADJUSTL(ivarStr2)), NF90_REAL, dimids, varID) )
     !  PRINT*, 'define good'
     CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
     !  PRINT*, 'put coordinates good'
     CALL check( nf90_put_att(ncID,varID,'units',TRIM(ADJUSTL(iunitStr2))) )
     !  PRINT*, 'put unit good'
     CALL check( nf90_put_att(ncID,varID,'longname',TRIM(ADJUSTL(ilongNmStr2))) )
     !  PRINT*, 'put longname good'
     idVar(iVar)=varID
  END DO
  CALL check( nf90_enddef(ncID) )
  ! End define mode. This tells netCDF we are done defining metadata.

  ! put all variable values into netCDF datasets
  ! put coordinate variables:
  CALL check( nf90_put_var(ncID, varIDx, varX) )
  CALL check( nf90_put_var(ncID, varIDy, varY) )
  CALL check( NF90_SYNC(ncID) )
  ! PRINT*, 'good put var'


  ! put grid_ID:
  CALL check( nf90_put_var(ncID, varIDGrid, RESHAPE(GridIDmatrix,(/nX,nY/),order = (/1,2/))) )
  ! PRINT*, 'good put varIDGrid',varIDGrid

  CALL check( NF90_SYNC(ncID) )

  ! then other 3D variables
  DO iVar = iVarStart, SIZE(UseColumnsDataOut), 1
     !  PRINT*, 'dim1:', SIZE(dataOut(1:irMax,iVar,:), dim=1)
     !  PRINT*, 'dim2:',SIZE(dataOut(1:irMax,iVar,:), dim=2)
     ! reshape dataOut to be aligned in checker board form
     varOut = RESHAPE(dataOutput(1:irMax,iVar,:),(/nX,nY,irMax/),order = (/3,1,2/) )
     varOut = varOut(:,nY:1:-1,:)
     !  get the variable id
     varID= idVar(iVar)
     !  PRINT*, 'good put iVar',iVar
     CALL check( nf90_put_var(ncID, varID, varOut) )
     !  PRINT*, 'good put var',varID
     CALL check(NF90_SYNC(ncID))
  END DO

  IF (ALLOCATED(varOut)) DEALLOCATE(varOut)
  IF (ALLOCATED(varSeq0)) DEALLOCATE(varSeq0)
  IF (ALLOCATED(varSeq)) DEALLOCATE(varSeq)
  IF (ALLOCATED(xLon)) DEALLOCATE(xLon)
  IF (ALLOCATED(xLat)) DEALLOCATE(xLat)
  IF (ALLOCATED(varY)) DEALLOCATE(varY)
  IF (ALLOCATED(varX)) DEALLOCATE(varX)

  ! IF (ALLOCATED(UseColumnsDataOut)) DEALLOCATE(UseColumnsDataOut)


  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  CALL check( nf90_close(ncID) )

  ! PRINT*, "*** SUCCESS writing netCDF file:"
  ! PRINT*, FileOut


END SUBROUTINE SUEWS_Write_nc

! SUBROUTINE SUEWS_DailyStateOut_nc(year_int, igrid)
!   USE sues_data
!   USE data_in
!   USE allocateArray
!   USE gis_data
!   USE time
!   USE defaultNotUsed
!   USE initial
!   USE solweig_module
!   USE cbl_module
!   USE ESTM_data
!   USE VegPhenogy
!   USE snowMod
!   USE netCDF
!
!   IMPLICIT NONE
!
!   INTEGER :: year_int, igrid
!   INTEGER,PARAMETER :: lenVarList=42
!   INTEGER :: idVar(lenVarList),dimIDs(3)
!   CHARACTER(len=20):: varNameList(lenVarList),varName,yrStr2,ivarStr2
!   CHARACTER(len=200)::rawpath
!   REAL:: ModelDailyStateOut(lenVarList),varOut
!   INTEGER :: x_dimid, y_dimid,time_dimid,iVar,varID,ncID,NX,NY
!
!   ! variable names
!   varNameList= (/&
!        "HDD1_h              ","HDD2_c              ","HDD3_Tmean          ",&
!        "HDD4_T5d            ","P/day               ","DaysSR              ",&
!        "GDD1_g              ","GDD2_s              ","GDD3_Tmin           ",&
!        "GDD4_Tmax           ","GDD5_DayLHrs        ","LAI_EveTr           ",&
!        "LAI_DecTr           ","LAI_Grass           ","DecidCap            ",&
!        "Porosity            ","AlbEveTr            ","AlbDecTr            ",&
!        "AlbGrass            ","WU_EveTr(1)         ","WU_EveTr(2)         ",&
!        "WU_EveTr(3)         ","WU_DecTr(1)         ","WU_DecTr(2)         ",&
!        "WU_DecTr(3)         ","WU_Grass(1)         ","WU_Grass(2)         ",&
!        "WU_Grass(3)         ","deltaLAI            ","LAIlumps            ",&
!        "AlbSnow             ","Dens_Snow_Paved     ","Dens_Snow_Bldgs     ",&
!        "Dens_Snow_EveTr     ","Dens_Snow_DecTr     ","Dens_Snow_Grass     ",&
!        "Dens_Snow_BSoil     ","Dens_Snow_Water     ","BoAnOHMEnd          ",&
!        "a1AnOHM             ","a2AnOHM             ","a3AnOHM             "/)
!   ! variable values
!   ModelDailyStateOut=&
!        (/HDD(id,1:6),GDD(id,1:5),&
!        LAI(id,1:nvegsurf),&
!        DecidCap(id),Porosity(id),AlbEveTr(id),AlbDecTr(id),AlbGrass(id),&
!        WU_day(id-1,1:9),&
!        deltaLAI,VegPhenLumps,SnowAlb,SnowDens(1:7),&
!        BoAnOHMEnd(igrid),a1AnOHM(igrid),a2AnOHM(igrid),a3AnOHM(igrid)/)
!
!
!   ! file names
!   WRITE(yrStr2,'(i4)') year_int
!   rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(yrStr2))
!   FileOut=TRIM(rawpath)//'_DailyState.nc'
!
!   ! if first day and first grid
!   ! create dailystate nc file
!   IF ( id==1 ) THEN
!      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
!      ! overwrite this file, if it already exists.
!      CALL check( nf90_create(FileOut, NF90_CLOBBER, ncID) )
!
!      ! Define the dimensions. NetCDF will hand back an ID for each.
!      ! nY = ncolumnsDataOut-4
!      ! nx = NumberOfGrids
!      CALL check( nf90_def_dim(ncID, "x", NX, x_dimid) )
!      CALL check( nf90_def_dim(ncID, "y", NY, y_dimid) )
!      CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )
!
!      ! The dimids array is used to pass the IDs of the dimensions of
!      ! the variables. Note that in fortran arrays are stored in
!      ! column-major format.
!      dimids =  (/x_dimid, y_dimid, time_dimid/)
!      ! define all variables
!      DO iVar = 1, lenVarList
!         ! define variable name
!         ! WRITE(ivarStr2,'(i10)') iVar
!         ivarStr2=varNameList(iVar)
!
!         ! Define the variable. The type of the variable in this case is
!         ! NF90_REAL.
!         CALL check( nf90_def_var(ncID,TRIM(ivarStr2), NF90_REAL, dimids, varID) )
!         idVar(iVar)=varID
!      END DO
!      CALL check( nf90_enddef(ncID) )
!      ! End define mode. This tells netCDF we are done defining metadata.
!
!      ! write initial DailyState into nc file
!      DO iVar = 1, lenVarList
!         WRITE(ivarStr2,'(i10)') iVar
!         ! write out each variable
!         varOut = ModelDailyStateOut(iVar)
!         !  get the variable id
!         varID= idVar(iVar)
!         CALL check( nf90_put_var(ncID, varID, varOut) )
!         CALL check( NF90_SYNC(ncID) )
!      END DO
!
!   ELSE
!      ! otherwise append DailyState to existing files
!      ! day 2 and ongoing
!      ! open DailyState nc file and check out variable information
!      CALL check( nf90_open(FileOut, nf90_nowrite, ncID) )
!
!      ! write initial DailyState into nc file
!      DO iVar = 1, lenVarList
!         ! get the variable id
!         varName=varNameList(iVar)
!         CALL check( nf90_inq_varid(ncID, varName, varID) )
!         ! write out each variable
!         varOut = ModelDailyStateOut(iVar)
!         CALL check( nf90_put_var(ncID, varID, varOut, start=(/id/)) )
!         CALL check( NF90_SYNC(ncID) )
!      END DO
!   END IF
!
! END SUBROUTINE SUEWS_DailyStateOut_nc

!===========================================================================!
! convert a vector of grids to a matrix
! the grid IDs in seqGrid2Sort follow the QGIS convention
! the spatial matrix arranges successive rows down the page (i.e., north to south)
!   and succesive columns across (i.e., west to east)
! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector
!===========================================================================!
SUBROUTINE grid2mat(seqGrid2Sort, seqGridSorted, matGrid, nRow, nCol)

  IMPLICIT NONE

  INTEGER,DIMENSION(nRow*nCol) :: seqGrid2Sort,seqGridSorted
  INTEGER,DIMENSION(nRow,nCol) :: matGrid
  INTEGER :: nRow, nCol,i,j,loc

  CALL sortGrid(seqGrid2Sort, seqGridSorted, nRow, nCol)
  PRINT*, 'old:'
  PRINT*, seqGrid2Sort(1:5)
  PRINT*, 'sorted:'
  PRINT*, seqGridSorted(1:5)
  PRINT*, ''
  DO i = 1, nRow
     DO j = 1, nCol
        loc=(i-1)*nCol+j
        ! PRINT*, i,j,loc
        ! PRINT*, seqGridSorted(loc)
        matGrid(i,j)=seqGridSorted(loc)
     END DO
  END DO

END SUBROUTINE grid2mat

!===========================================================================!
! convert sequence of REAL values to a matrix
! the grid IDs in seqGrid2Sort follow the QGIS convention
! the spatial matrix arranges successive rows down the page (i.e., north to south)
!   and succesive columns across (i.e., west to east)
! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector
!===========================================================================!
SUBROUTINE seq2mat(seq2Sort, seqSorted, matGrid, nRow, nCol)

  IMPLICIT NONE

  REAL(KIND(1d0)),DIMENSION(nRow*nCol) :: seq2Sort,seqSorted
  REAL(KIND(1d0)),DIMENSION(nRow,nCol) :: matGrid
  INTEGER :: nRow, nCol,i,j,loc

  CALL sortSeqReal(seq2Sort, seqSorted, nRow, nCol)
  PRINT*, 'old:'
  PRINT*, seq2Sort(1:5)
  PRINT*, 'sorted:'
  PRINT*, seqSorted(1:5)
  PRINT*, ''
  DO i = 1, nRow
     DO j = 1, nCol
        loc=(i-1)*nCol+j
        ! PRINT*, i,j,loc
        ! PRINT*, seqGridSorted(loc)
        matGrid(i,j)=seqSorted(loc)
     END DO
  END DO

END SUBROUTINE seq2mat

!===========================================================================!
! sort a sequence of LONG values into the specially aligned sequence per QGIS
!===========================================================================!
SUBROUTINE sortGrid(seqGrid2Sort, seqGridSorted, nRow, nCol)
  USE qsort_c_module
  ! convert a vector of grids to a matrix
  ! the grid IDs in seqGrid2Sort follow the QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector

  IMPLICIT NONE
  INTEGER :: nRow, nCol,i=1,j=1,xInd,len

  INTEGER,DIMENSION(nRow*nCol),INTENT(in) :: seqGrid2Sort
  INTEGER,DIMENSION(nRow*nCol),INTENT(out) :: seqGridSorted
  INTEGER,DIMENSION(nRow*nCol) :: locSorted
  INTEGER :: loc
  REAL:: ind(nRow*nCol,2)
  REAL :: seqGridSortedReal(nRow*nCol),val

  ! number of grids
  len=nRow*nCol
  ! PRINT*, 'after input:'
  ! PRINT*, 'seqGrid2Sort:'
  ! PRINT*, seqGrid2Sort(1:5)
  ! PRINT*, '****'
  !
  !
  ! PRINT*, 'after input:'
  ! PRINT*, 'seqGridSorted:'
  ! PRINT*, seqGridSorted(1:5)
  ! PRINT*, '****'


  ! fill in an nRow*nCol array with values to determine sequence
  xInd=1
  DO i = 1, nRow
     DO j = 1, nCol
        !  {row, col, value for sorting, index in new sequence}
        ind(xInd,:)=(/i+j+i/(nRow+1.),xInd*1./)
        xInd=xInd+1
     END DO
  END DO
  ! PRINT*, 'old after sorting:'
  ! PRINT*, seqGridSorted(1:5)
  ! PRINT*, 'ind:'
  ! PRINT*, ind(:,1)
  ! PRINT*, 'ind. seq:'
  ! PRINT*, ind(:,2)

  ! then sorted ind(:,3) will have the same order as seqGrid2Sort
  ! sort ind(:,3)
  seqGridSortedReal=ind(:,1)*1.
  CALL QsortC(seqGridSortedReal)
  ! print*, 'sorted real:'
  ! print*, seqGridSortedReal

  ! get index of each element of old sequence in the sorted sequence
  DO i=1,len
     ! value in old sequence
     !  val=ind(i,3)*1.
     val=seqGridSortedReal(i)
     DO j=1,len
        IF ( val .EQ. ind(j,1)*1.) THEN
           ! location in sorted sequence
           locSorted(i)=j
        END IF
     END DO
  END DO

  ! put elements of old sequence in the sorted order
  DO i = 1, len
     loc=locSorted(i)
     seqGridSorted(loc)=seqGrid2Sort(i)
  END DO
  seqGridSorted=seqGridSorted(len:1:-1)
  ! PRINT*, 'loc sorted:'
  ! PRINT*, locSorted

  ! PRINT*, 'sorted:'
  ! PRINT*, seqGridSorted(1:5)
  ! PRINT*, 'sort subroutine end!'

END SUBROUTINE sortGrid

!===========================================================================!
! sort a sequence of REAL values into the specially aligned sequence per QGIS
!===========================================================================!
SUBROUTINE sortSeqReal(seqReal2Sort, seqRealSorted, nRow, nCol)
  USE qsort_c_module
  ! convert a vector of grids to a matrix
  ! the grid IDs in seqReal2Sort follow the QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! seqRealSorted stores the grid IDs as aligned in matGrid but squeezed into a vector

  IMPLICIT NONE
  INTEGER :: nRow, nCol,i=1,j=1,xInd,len

  REAL(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(in) :: seqReal2Sort
  REAL(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(out) :: seqRealSorted
  INTEGER(KIND(1d0)),DIMENSION(nRow*nCol) :: locSorted
  INTEGER(KIND(1d0)) :: loc
  REAL:: ind(nRow*nCol,2)
  REAL :: seqRealSortedReal(nRow*nCol),val

  ! number of grids
  len=nRow*nCol
  ! PRINT*, 'after input:'
  ! PRINT*, 'seqReal2Sort:'
  ! PRINT*, seqReal2Sort(1:5)
  ! PRINT*, '****'
  !
  !
  ! PRINT*, 'after input:'
  ! PRINT*, 'seqRealSorted:'
  ! PRINT*, seqRealSorted(1:5)
  ! PRINT*, '****'


  ! fill in an nRow*nCol array with values to determine sequence
  xInd=1
  DO i = 1, nRow
     DO j = 1, nCol
        !  {row, col, value for sorting, index in new sequence}
        ind(xInd,:)=(/i+j+i/(nRow+1.),xInd*1./)
        xInd=xInd+1
     END DO
  END DO
  ! PRINT*, 'old after sorting:'
  ! PRINT*, seqRealSorted(1:5)
  ! PRINT*, 'ind:'
  ! PRINT*, ind(:,1)
  ! PRINT*, 'ind. seq:'
  ! PRINT*, ind(:,2)

  ! then sorted ind(:,3) will have the same order as seqReal2Sort
  ! sort ind(:,3)
  seqRealSortedReal=ind(:,1)*1.
  CALL QsortC(seqRealSortedReal)
  ! print*, 'sorted real:'
  ! print*, seqRealSortedReal

  ! get index of each element of old sequence in the sorted sequence
  DO i=1,len
     ! value in old sequence
     !  val=ind(i,3)*1.
     val=seqRealSortedReal(i)
     DO j=1,len
        IF ( val .EQ. ind(j,1)*1.) THEN
           ! location in sorted sequence
           locSorted(i)=j
        END IF
     END DO
  END DO

  ! put elements of old sequence in the sorted order
  DO i = 1, len
     loc=locSorted(i)
     seqRealSorted(loc)=seqReal2Sort(i)
  END DO
  seqRealSorted=seqRealSorted(len:1:-1)
  ! PRINT*, 'loc sorted:'
  ! PRINT*, locSorted

  ! PRINT*, 'sorted:'
  ! PRINT*, seqRealSorted(1:5)
  ! PRINT*, 'sort subroutine end!'

END SUBROUTINE sortSeqReal

!===========================================================================!
! a wrapper for checking netCDF status
!===========================================================================!
SUBROUTINE check(status)
  USE netcdf

  INTEGER, INTENT ( in) :: status

  IF(status /= nf90_noerr) THEN
     PRINT *, TRIM(nf90_strerror(status))
     STOP "Stopped"
  END IF
END SUBROUTINE check
