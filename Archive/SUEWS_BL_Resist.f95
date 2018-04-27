SUBROUTINE BoundaryLayerResistance(&
     zzd,& ! input:    !Active measurement height (meas. height-displac. height)
     z0M,&     !Aerodynamic roughness length
     avU1,&    !Average wind speed
     UStar,&! input/output:
     rb)! output:

  IMPLICIT NONE

  REAL(KIND(1d0)),INTENT(in)::zzd     !Active measurement height (meas. height-displac. height)
  REAL(KIND(1d0)),INTENT(in)::z0M     !Aerodynamic roughness length
  REAL(KIND(1d0)),INTENT(in)::avU1    !Average wind speed

  REAL(KIND(1d0)),INTENT(inout)::UStar!Friction velocity

  REAL(KIND(1d0)),INTENT(out)::rb   !boundary layer resistance shuttleworth

  REAL(KIND(1d0)),PARAMETER :: k=0.4

  IF(UStar<0.01) THEN
     UStar=avu1/LOG(zzd/z0m)*k
  END IF

  rb=(1.1/UStar)+(5.6*(UStar**0.333333))!rb - boundary layer resistance shuttleworth

  RETURN
END SUBROUTINE BoundaryLayerResistance
