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
  USE netCDF


  IMPLICIT NONE

  INTEGER::i ,azero!,lfnOutC
  INTEGER:: Gridiv, year_int, iv, irMax, GridID
  CHARACTER(len=10):: tstepStr2, ivStr2, yrStr2
  CHARACTER(len=100):: rawpath, SnowOut,ESTMOut



  ! This is the name of the data file we will create.
  ! character (len = *), parameter :: FILE_NAME = "simple_xy.nc"

  ! We are writing 3D data, {time, y, x}
  integer, parameter :: NDIMS = 3
  INTEGER :: NX, NY

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  INTEGER :: ncid, varid, dimids(NDIMS)
  INTEGER :: x_dimid, y_dimid,time_dimid

  ! This is the data array we will write. It will just be filled with
  ! a progression of integers for this example.
  ! integer :: data_out(NY, NX)

  ! Loop indexes, and error handling.
  ! integer :: x, y

  !================DEFINE OUTPUT FILENAME AND ITS PATH================
  WRITE(tstepStr2,'(i2)') TSTEP/60
  WRITE(ivStr2,'(i10)') iv
  WRITE(yrStr2,'(i4)') year_int


  rawpath=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(yrStr2))//'_'//TRIM(ADJUSTL(ivStr2))
  FileOut=TRIM(rawpath)//'_'//TRIM(ADJUSTL(tstepStr2))//'.nc'
  ! SOLWEIGpoiOut=TRIM(rawpath)//'_SOLWEIGpoiOut.txt'
  ! ESTMOut=TRIM(rawpath)//'_ESTM_5.txt'
  ! BLOut=TRIM(rawpath)//'_BL.txt'
  ! SnowOut=TRIM(rawpath)//'_snow_5.txt'

  ! Create some pretend data. If this wasn't an example program, we
  ! would have some real data to write, for example, model output.
  ! do x = 1, NX
  !    do y = 1, NY
  !       data_out(y, x) = (x - 1) * NY + (y - 1)
  !    end do
  ! end do

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  CALL check( nf90_create(FileOut, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each.
  nY=ncolumnsDataOut-4
  nx=NumberOfGrids
  CALL check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, time_dimid) )
  CALL check( nf90_def_dim(ncid, "x", NX, x_dimid) )
  CALL check( nf90_def_dim(ncid, "y", NY, y_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ y_dimid, x_dimid, time_dimid/)

  ! Define the variable. The type of the variable in this case is
  ! NF90_REAL.
  CALL check( nf90_def_var(ncid, "data", NF90_REAL, dimids, varid) )
  print*, 'GOOD 1'

  ! End define mode. This tells netCDF we are done defining metadata.
  CALL check( nf90_enddef(ncid) )
  print*, 'GOOD 2'


  ! CALL check( nf90_inquire_variable(ncid, varid, ndims = numDims) )
  ! dimids

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.

  ! CALL check( nf90_put_var(ncid, varid, dataOut(1:irMax,5:ncolumnsDataOut,:)) )

  CALL check( nf90_put_var(ncid, varid, reshape(dataOut(1:irMax,5:ncolumnsDataOut,:),(/ny, nx, irmax/) ) ))
  print*, 'GOOD 3'

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  CALL check( nf90_close(ncid) )

  PRINT *, "*** SUCCESS writing example file simple_xy.nc! "




END SUBROUTINE SUEWS_Output_nc



SUBROUTINE check(status)
  use netcdf

  INTEGER, INTENT ( in) :: status

  IF(status /= nf90_noerr) THEN
     PRINT *, TRIM(nf90_strerror(status))
     STOP "Stopped"
  END IF
END SUBROUTINE check
