subroutine BoundaryLayerResistance
  use sues_data
  use data_in
  use mod_z	        ! module_LUMPS_constants,f90
  use mod_k	        ! module_LUMPS_constants,f90
  implicit none

  if(ustar<0.01) then
     ustar=avu1/log(zzd/z0m)*k
  end if
  
  rb=(1.1/ustar)+(5.6*(ustar**0.333333))!rb - boundary layer resistance shuttleworth
  
  return
end subroutine BoundaryLayerResistance
