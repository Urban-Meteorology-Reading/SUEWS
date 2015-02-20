!Interpolates hourly profiles provided in SUEWS_Profiles.txt
! to resolution of the model timestep
! HCW 06 Feb 2015
!===================================================================================
subroutine SUEWS_InterpHourlyProfiles(Gridiv,TstepP_ID,SurfChar_HrProf)

  use allocateArray
  use ColNamesInputFiles
  use sues_data

  IMPLICIT NONE

  integer:: i,j, ii   !Used to count over hours and sub-hourly timesteps
  integer:: Gridiv, TstepP_ID
  integer,dimension(24):: SurfChar_HrProf
  real(kind(1d0)):: deltaProf   !Change in hourly profiles per model timestep               

  ! Copy value for first hour
  TstepProfiles(Gridiv,TstepP_ID,1) = SurfaceChar(Gridiv,SurfChar_HrProf(1))
  do i=1,24
      j = (i+1)  
      if(i == 24) j = 1   !If last hour of day, loop round to first hour of day for interpolation
      deltaProf = ((SurfaceChar(Gridiv,SurfChar_HrProf(j)) - SurfaceChar(Gridiv,SurfChar_HrProf(i))))/nsh_real
      do ii=1,nsh
         if((nsh*(i-1)+ii+1) < (23*nsh+nsh+1))  then
            TstepProfiles(Gridiv,TstepP_ID,(nsh*(i-1)+ii+1)) = SurfaceChar(Gridiv,SurfChar_HrProf(i)) + deltaProf*ii              
         endif
      enddo         
   enddo
   
endsubroutine SUEWS_InterpHourlyProfiles
!===================================================================================