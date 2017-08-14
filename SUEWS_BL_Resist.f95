SUBROUTINE BoundaryLayerResistance(&
     ! input:
     zzd,&     !Active measurement height (meas. height-displac. height)
     z0M,&     !Aerodynamic roughness length
     avU1,&    !Average wind speed

     ! input/output:
     USTAR,&
     ! output:
     rb)
  ! use sues_data
  ! use data_in
  ! use mod_z       ! module_LUMPS_constants,f90
  ! use mod_k       ! module_LUMPS_constants,f90
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(in)::&
       zzd,&     !Active measurement height (meas. height-displac. height)
       z0M,&     !Aerodynamic roughness length
       avU1    !Average wind speed

  REAL(KIND(1d0)),INTENT(inout)::&
       USTAR!Friction velocity

  REAL(KIND(1d0)),INTENT(out)::&
       rb   !boundary layer resistance shuttleworth

  REAL(KIND(1d0)),PARAMETER :: &
       k=0.4,&             !Von Karman's contant
       k2=0.16,&           !Power of Van Karman's contant
       neut_limit=0.001000,& !Limit for neutral stability
       grav=9.80665,&  !g - gravity - physics today august 1987
       notUsedI=-55

  IF(ustar<0.01) THEN
     ustar=avu1/LOG(zzd/z0m)*k
  END IF

  rb=(1.1/ustar)+(5.6*(ustar**0.333333))!rb - boundary layer resistance shuttleworth

  RETURN
END SUBROUTINE BoundaryLayerResistance
