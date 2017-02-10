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
  varSeq0=SiteSelect(:,5)
  CALL sortSeqReal(varSeq0,varSeq,nY,nX)
  xLat = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
  ! PRINT*, 'before flipping:',xLat(1:5,1)
  xLat =xLat(:,nY:1:-1)
  ! PRINT*, 'after flipping:',xLat(1:5,1)

  ! longitude:
  varSeq0=SiteSelect(:,6)
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
     varSeq0=SiteSelect(:,iVar)
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
!      varSeq0=SiteSelect(:,iVar)
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
SUBROUTINE SUEWS_Output_nc(year_int,iv,irMax)
  !INPUT: year_int = Year as a integer
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
  ! USE SetupOutput
  USE netCDF


  IMPLICIT NONE

  INTEGER:: year_int, iv, irMax, i
  CHARACTER(len=10):: tstepStr2, ivStr2, yrStr2
  CHARACTER(len=100):: rawpath

  ! We are writing 3D data, {time, y, x}
  INTEGER, PARAMETER :: NDIMS = 3, iVarStart=6
  INTEGER :: NX,NY

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  INTEGER :: ncID, varID, dimids(NDIMS),varIDGrid
  INTEGER :: x_dimid, y_dimid,time_dimid,iVar,varIDx,varIDy
  REAL(KIND(1d0)), ALLOCATABLE :: varOut(:,:,:),&
       varX(:,:),varY(:,:),&
       xLat(:,:),xLon(:,:),&
       varSeq0(:),varSeq(:)
  INTEGER :: idVar(iVarStart:ncolumnsDataOut)
  CHARACTER(len=20):: nameVarList(iVarStart:ncolumnsDataOut),&
       unitVarList(iVarStart:ncolumnsDataOut),&
       ivarStr2,iunitStr2
  CHARACTER(len=12*nColumnsDataOut):: HeaderOut, UnitsOut   !Header and units for selected output variables (untrimmed)
  CHARACTER(len=12*nColumnsDataOut):: HeaderUsed, UnitsUsed   !Header and units for selected output variables (untrimmed)


  !================DEFINE OUTPUT VARIABLES AND UNITS================
  ! variable names:
  nameVarList=(/&
       '     kdown','       kup','     ldown','       lup','     Tsurf', &   !radiation components (6-10)
       '        qn','        qf','        qs','        qh','        qe', &   !energy fluxes (11-15)
       '  qh_LUMPS','  qe_LUMPS','      qh_r', &                             !energy fluxes (other approaches) (16-18)
       '       P_i','      Ie_i','       E_i','      RO_i','   totCh_i', &   !water balance components (19-23)
       '  surfCh_i','      St_i','    NWSt_i','      Dr_i','       smd', &   !water balance components cont. (24-28)
       '    FlowCh','  AddWater', &                                          !water balance components cont. (29-30)
       '    ROsoil','    ROpipe','    RO_imp','    RO_veg','    RO_wat', &   !runoff components (31-35)
       '    wu_int','  wu_EveTr','  wu_DecTr','  wu_Grass', &                !water use (36-39)
       ' smd_Paved',' smd_Bldgs',' smd_EveTr',' smd_DecTr',' smd_Grass',' smd_BSoil', &   !smd for each surface (40-45)
       '  St_Paved','  St_Bldgs','  St_EveTr','  St_DecTr','  St_Grass','  St_BSoil','  St_Water',&   !states (46-52)
       '    zenith','   azimuth','  alb_bulk','      Fcld', &                ! extra radiation info (53-56)
       '       LAI','       z0m','       zdm', &                             ! extra surface info (57-59)
       '     ustar','       Lob','        ra','        rs', &                ! turbulence (60-63)
       '        Fc', &                                                       ! CO2 flux (64)
       '  Fc_photo','  Fc_respi','  Fc_metab','  Fc_traff','  Fc_build', &   ! CO2 flux components (65-69)
       '   qn_SnFr','  qn_Sn   ','  alb_snow', &                             ! snow-related (radiation) (70-72)
       '        qm','  qmFreeze','    qmRain','       SWE','        Mw','   MwStore','    SnCh_i', &   !snow (73-79)
       ' SnR_Paved',' SnR_Bldgs'&                                            !snow-related (removal) (80-81)
       /)

  unitVarList=(/&
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

  ! Select variables to be written out
  IF(WriteOutOption == 0) THEN   !all (not snow-related)
     IF ( .NOT. ALLOCATED(UseColumnsDataOut) ) THEN
        ALLOCATE(UseColumnsDataOut(69-iVarStart+1))
     END IF
     UsecolumnsDataOut = (/ (i, i=iVarStart,69, 1) /)
  ELSEIF(WriteOutOption == 1) THEN   !all plus snow-related
     IF ( .NOT. ALLOCATED(UseColumnsDataOut) ) THEN
        ALLOCATE(UseColumnsDataOut(nColumnsDataOut-iVarStart+1))
     END IF

     UsecolumnsDataOut = (/ (i, i=iVarStart,nColumnsDataOut, 1) /)
  ELSEIF(WriteOutOption == 2) THEN   !minimal output
     IF ( .NOT. ALLOCATED(UseColumnsDataOut) ) THEN
        ALLOCATE(UseColumnsDataOut(33-iVarStart+1))
     END IF
     UsecolumnsDataOut = (/ (i, i=iVarStart,15, 1),(i, i=19,28, 1), 53,54,55,56, 57, 60,61, 64 /)
  ELSE
     WRITE(*,*) 'RunControl: WriteOutOption code not recognised, so writing out all variables.'
     IF ( .NOT. ALLOCATED(UseColumnsDataOut) ) THEN
        ALLOCATE(UseColumnsDataOut(69-iVarStart+1))
     END IF
     UsecolumnsDataOut = (/ (i, i=iVarStart,69, 1) /)
  ENDIF

  ! Create subset of HeaderAll and FormatAll for selected variables only
  HeaderOut=''
  UnitsOut=''
  ! FormatOut=''
  ! AggOut=''
  ! ColNos=''
  DO i=1,SIZE(UseColumnsDataOut)
     HeaderOut=TRIM(HeaderOut)//' '//ADJUSTL(nameVarList(UsecolumnsDataOut(i)))
     !write(*,*) HeaderOut
     UnitsOut=TRIM(UnitsOut)//' '//ADJUSTL(unitVarList(UsecolumnsDataOut(i)))
     !write(*,*) UnitsOut
     !  FormatOut=TRIM(FormatOut)//' '//ADJUSTL(FormatAll(UsecolumnsDataOut(i)))
     !write(*,*) FormatOut
     !  AggOut=TRIM(AggOut)//' '//ADJUSTL(AggAll(UsecolumnsDataOut(i)))
     !write(*,*) AggOut
     !  WRITE(itext,'(i3)') i
     !  ColNos=TRIM(ColNos)//' '//ADJUSTL(itext)
  ENDDO

  ! PRINT*, 'mem start'
  ! ALLOCATE(CHARACTER(LEN(TRIM(ADJUSTL(HeaderOut)))):: HeaderUse)
  ! PRINT*, 'mem', 1
  ! ALLOCATE(CHARACTER(LEN(TRIM(ADJUSTL(UnitsOut)))):: UnitsUse)
  ! PRINT*, 'mem', 2
  ! ALLOCATE(CHARACTER(LEN(trim(adjustl(FormatOut)))):: FormatUse)
  ! ALLOCATE(CHARACTER(LEN(trim(adjustl(AggOut)))):: AggUse)
  ! ALLOCATE(CHARACTER(LEN(trim(adjustl(ColNos)))):: ColNosUse)
  HeaderUsed=TRIM(ADJUSTL(HeaderOut))
  UnitsUsed=TRIM(ADJUSTL(UnitsOut))
  ! FormatUse='('//trim(adjustl(FormatOut))//')'
  ! AggUse=trim(adjustl(AggOut))
  ! ColNosUse=trim(adjustl(ColNos))


  !================DEFINE OUTPUT FILENAME AND ITS PATH================
  WRITE(tstepStr2,'(i2)') TSTEP/60
  WRITE(ivStr2,'(i10)') iv
  WRITE(yrStr2,'(i4)') year_int

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
  ! PRINT*, 'size x dim 1:',SIZE(varX, dim=1)
  ! PRINT*, 'size x dim 2:',SIZE(varX, dim=2)

  ! file names
  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(yrStr2))//'_'//TRIM(ADJUSTL(ivStr2))
  FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(tstepStr2))//'.nc'

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  CALL check( nf90_create(FileOut, NF90_CLOBBER, ncID) )

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
  DO iVar = iVarStart, ncolumnsDataOut, 1
     ! define variable name
     ivarStr2=nameVarList(iVar)
     iunitStr2=unitVarList(iVar)
     !  PRINT*, iunitStr2

     ! Define the variable. The type of the variable in this case is
     ! NF90_REAL.
     !  print*, TRIM(ADJUSTL(ivarStr2))
     CALL check( nf90_def_var(ncID,TRIM(ADJUSTL(ivarStr2)), NF90_REAL, dimids, varID) )
     !  print*, 'define good'
     CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
     !  print*, 'put att good'
     CALL check( nf90_put_att(ncID,varID,'units',TRIM(ADJUSTL(iunitStr2))) )
     !  print*, 'put unit good'
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
  CALL check( NF90_SYNC(ncID) )

  ! then other 3D variables
  DO iVar = iVarStart, ncolumnsDataOut, 1
     !  PRINT*, 'dim1:', SIZE(dataOut(1:irMax,iVar,:), dim=1)
     !  PRINT*, 'dim2:',SIZE(dataOut(1:irMax,iVar,:), dim=2)
     ! reshape dataOut to be aligned in checker board form
     varOut = RESHAPE(dataOut(1:irMax,iVar,:),(/nX,nY,irMax/),order = (/3,1,2/) )
     varOut = varOut(:,nY:1:-1,:)
     !  get the variable id
     varID= idVar(iVar)
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

  ! IF (ALLOCATED(UseColumnsDataOut)) DEALLOCATE(UseColumnsDataOut)


  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  CALL check( nf90_close(ncID) )

  ! PRINT*, "*** SUCCESS writing netCDF file:"
  ! PRINT*, FileOut


END SUBROUTINE SUEWS_Output_nc


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

  INTEGER(KIND(1d0)),DIMENSION(nRow*nCol) :: seqGrid2Sort,seqGridSorted
  INTEGER(KIND(1d0)),DIMENSION(nRow,nCol) :: matGrid
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

  INTEGER(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(in) :: seqGrid2Sort
  INTEGER(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(out) :: seqGridSorted
  INTEGER(KIND(1d0)),DIMENSION(nRow*nCol) :: locSorted
  INTEGER(KIND(1d0)) :: loc
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
