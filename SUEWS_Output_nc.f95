!In this subroutine the output files will be opened and the output matrices will be printed out.
!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_wr.f90,v 1.7 2006/12/09 18:44:58 russ Exp $
!-----------------------------------------------------------------------------------------------
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
  USE netCDF


  IMPLICIT NONE

  ! INTEGER::i ,azero!,lfnOutC
  INTEGER:: year_int, iv, irMax!, GridID
  CHARACTER(len=10):: tstepStr2, ivStr2, yrStr2,ivarStr2
  CHARACTER(len=100):: rawpath!, SnowOut,ESTMOut

  ! We are writing 3D data, {time, y, x}
  INTEGER, PARAMETER :: NDIMS = 3
  INTEGER :: NX,NY

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  INTEGER :: ncID, varID, dimids(NDIMS)
  INTEGER :: x_dimid, y_dimid,time_dimid,iVar
  REAL, ALLOCATABLE :: varOut(:,:,:)
  INTEGER :: idVar(5:ncolumnsDataOut)

  ! This is the data array we will write. It will just be filled with
  ! a progression of integers for this example.
  ! integer :: data_out(NY, NX)

  ! Loop indexes, and error handling.
  ! integer :: x, y

  !================DEFINE OUTPUT FILENAME AND ITS PATH================
  WRITE(tstepStr2,'(i2)') TSTEP/60
  WRITE(ivStr2,'(i10)') iv
  WRITE(yrStr2,'(i4)') year_int

  ! define the dimension of spatial array/frame in the output
  nX=nCol
  nY=nRow

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
  CALL check( nf90_def_dim(ncID, "x", NX, x_dimid) )
  CALL check( nf90_def_dim(ncID, "y", NY, y_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/x_dimid, y_dimid, time_dimid/)

  ! write out each variable
  ALLOCATE(varOut(nX,nY,irMax))

  ! define all variables
  DO iVar = 5, ncolumnsDataOut, 1
     ! define variable name
     WRITE(ivarStr2,'(i10)') iVar

     ! Define the variable. The type of the variable in this case is
     ! NF90_REAL.
     CALL check( nf90_def_var(ncID,'var'//TRIM(ivarStr2), NF90_REAL, dimids, varID) )
     idVar(iVar)=varID
  END DO
  CALL check( nf90_enddef(ncID) )
  ! End define mode. This tells netCDF we are done defining metadata.


  ! put all variable values into netCDF datasets
  DO iVar = 5, ncolumnsDataOut, 1
     WRITE(ivarStr2,'(i10)') iVar
     varOut = RESHAPE(dataOut(1:irMax,iVar,:),(/nX,nY,irMax/),order = (/2,1,3/) )
     !  get the variable id
     varID= idVar(iVar)
     CALL check( nf90_put_var(ncID, varID, varOut) )
     CALL check(NF90_SYNC(ncID))
  END DO
  IF (ALLOCATED(varOut)) DEALLOCATE(varOut)

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  CALL check( nf90_close(ncID) )

  ! PRINT*, "*** SUCCESS writing netCDF file:"
  ! PRINT*, FileOut

END SUBROUTINE SUEWS_Output_nc


SUBROUTINE SUEWS_DailyStateOut_nc(year_int, igrid)
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
  USE VegPhenogy
  USE snowMod
  USE netCDF

  INTEGER :: year_int, igrid
  INTEGER,PARAMETER :: lenVarList=42
  INTEGER :: idVar(lenVarList),dimIDs(3)
  CHARACTER(len=20):: varNameList(lenVarList),varName,yrStr2,ivarStr2
  CHARACTER(len=200)::rawpath
  REAL:: ModelDailyStateOut(lenVarList),varOut
  INTEGER :: x_dimid, y_dimid,time_dimid,iVar,varID




  ! variable names
  varNameList= &
       (/"HDD1_h              ","HDD2_c              ","HDD3_Tmean          ",&
       "HDD4_T5d            ","P/day               ","DaysSR              ",&
       "GDD1_g              ","GDD2_s              ","GDD3_Tmin           ",&
       "GDD4_Tmax           ","GDD5_DayLHrs        ","LAI_EveTr           ",&
       "LAI_DecTr           ","LAI_Grass           ","DecidCap            ",&
       "Porosity            ","AlbEveTr            ","AlbDecTr            ",&
       "AlbGrass            ","WU_EveTr(1)         ","WU_EveTr(2)         ",&
       "WU_EveTr(3)         ","WU_DecTr(1)         ","WU_DecTr(2)         ",&
       "WU_DecTr(3)         ","WU_Grass(1)         ","WU_Grass(2)         ",&
       "WU_Grass(3)         ","deltaLAI            ","LAIlumps            ",&
       "AlbSnow             ","Dens_Snow_Paved     ","Dens_Snow_Bldgs     ",&
       "Dens_Snow_EveTr     ","Dens_Snow_DecTr     ","Dens_Snow_Grass     ",&
       "Dens_Snow_BSoil     ","Dens_Snow_Water     ","BoAnOHMEnd          ",&
       "a1AnOHM             ","a2AnOHM             ","a3AnOHM             "/)
  ! variable values
  ModelDailyStateOut=&
       (/HDD(id,1:6),GDD(id,1:5),&
       LAI(id,1:nvegsurf),&
       DecidCap(id),Porosity(id),AlbEveTr(id),AlbDecTr(id),AlbGrass(id),&
       WU_day(id-1,1:9),&
       deltaLAI,VegPhenLumps,SnowAlb,SnowDens(1:7),&
       BoAnOHMEnd(igrid),a1AnOHM(igrid),a2AnOHM(igrid),a3AnOHM(igrid)/)


  ! file names
  WRITE(yrStr2,'(i4)') year_int
  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(yrStr2))
  FileOut=TRIM(rawpath)//'_DailyState.nc'

  ! if first day and first grid
  ! create dailystate nc file
  IF ( id==1 ) THEN
     ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
     ! overwrite this file, if it already exists.
     CALL check( nf90_create(FileOut, NF90_CLOBBER, ncID) )

     ! Define the dimensions. NetCDF will hand back an ID for each.
     ! nY = ncolumnsDataOut-4
     ! nx = NumberOfGrids
     CALL check( nf90_def_dim(ncID, "x", NX, x_dimid) )
     CALL check( nf90_def_dim(ncID, "y", NY, y_dimid) )
     CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )

     ! The dimids array is used to pass the IDs of the dimensions of
     ! the variables. Note that in fortran arrays are stored in
     ! column-major format.
     dimids =  (/x_dimid, y_dimid, time_dimid/)
     ! define all variables
     DO iVar = 1, lenVarList
        ! define variable name
        ! WRITE(ivarStr2,'(i10)') iVar
        ivarStr2=varNameList(iVar)

        ! Define the variable. The type of the variable in this case is
        ! NF90_REAL.
        CALL check( nf90_def_var(ncID,TRIM(ivarStr2), NF90_REAL, dimids, varID) )
        idVar(iVar)=varID
     END DO
     CALL check( nf90_enddef(ncID) )
     ! End define mode. This tells netCDF we are done defining metadata.

     ! write initial DailyState into nc file
     DO iVar = 1, lenVarList
        WRITE(ivarStr2,'(i10)') iVar
        ! write out each variable
        varOut = ModelDailyStateOut(iVar)
        !  get the variable id
        varID= idVar(iVar)
        CALL check( nf90_put_var(ncID, varID, varOut) )
        CALL check( NF90_SYNC(ncID) )
     END DO

  ELSE
     ! otherwise append DailyState to existing files
     ! day 2 and ongoing
     ! open DailyState nc file and check out variable information
     CALL check( nf90_open(FileOut, nf90_nowrite, ncID) )

     ! write initial DailyState into nc file
     DO iVar = 1, lenVarList
        ! get the variable id
        varName=varNameList(iVar)
        CALL check( nf90_inq_varid(ncID, varName, varID) )
        ! write out each variable
        varOut = ModelDailyStateOut(iVar)
        CALL check( nf90_put_var(ncID, varID, varOut, start=(/id/)) )
        CALL check( NF90_SYNC(ncID) )
     END DO
  END IF

END SUBROUTINE SUEWS_DailyStateOut_nc

SUBROUTINE grid2mat(seqGridOld, seqGridSorted, matGrid, nRow, nCol)
  ! convert a vector of grids to a matrix
  ! the grid IDs in seqGridOld follow the QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector


  INTEGER,DIMENSION(nRow*nCol) :: seqGridOld,seqGridSorted
  INTEGER,DIMENSION(nRow,nCol) :: matGrid




END SUBROUTINE grid2mat


SUBROUTINE check(status)
  USE netcdf

  INTEGER, INTENT ( in) :: status

  IF(status /= nf90_noerr) THEN
     PRINT *, TRIM(nf90_strerror(status))
     STOP "Stopped"
  END IF
END SUBROUTINE check
