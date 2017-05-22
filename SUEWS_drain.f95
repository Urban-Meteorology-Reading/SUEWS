 subroutine drainage(StorCap,DrainEq,DrainCoef1,DrainCoef2)
!Calculation of drainage for each land surface. 
!INPUT: Storage capacity, type of drainage equation used, drainage coefficients
!       used in the equation
!Modified by HCW 16 Feb 2015
!  Removed option of Eq 4 (calculation needs to be checked before re-implementing).
!  Code writes an error if calculated drainage exceeds surface state (but code continues).
!  This may indicate inappropriate drainage equation, storage capacities or model tstep.
!Modified by LJ in Aug 2011. Drainage cannot exceed the surface storage. 
!Modified LJ in 10/2010
!------------------------------------------------------------------------------    
  
 use allocateArray
 use gis_data
 use sues_data
 use time
  
 IMPLICIT NONE
  
 real (kind(1d0))::StorCap,DrainCoef1,DrainCoef2,DrainEq
  
 !If surface is dry, no drainage occurs
 if(state(is)<0.000000001) then
    drain(is)=0.0
 else
    if(int(DrainEq)==1) then   !Falk and Niemczynowicz (1978): Drainage equation for paved, buildings and irrigated grass

       if (state(is)<StorCap) then
         drain(is)=0   !No drainage if state is less than storage capacity
       else
         drain(is)=(DrainCoef1*(state(is)-StorCap)**DrainCoef2)/nsh_real
       endif
      
    elseif(int(DrainEq)==2) then   !Rutter eqn corrected for c=0, see Eq 9 of Calder & Wright 1986
       drain(is)=(DrainCoef1*(exp(DrainCoef2*state(is))-1))/nsh_real
        ! N.B. -1 is correct here but brackets are wrong in G&O 1991 Eq 5 & Ja11 Eq 18.   
      
    elseif(int(DrainEq)==3) then   !Falk and Niemczynowicz (1978)
        drain(is)=(DrainCoef1*(state(is)**DrainCoef2))/nsh_real
            
     ! Option 4 removed by HCW 16 Feb 2015, as it is not used and appears to be problematic
     !elseif(int(DrainEq)==4) then    !Rutter eqn not corrected for c=0
     !   drain(is)=DrainCoef1*exp(DrainCoef2*(state(is)-StorCap))
     !   drain(is)=drain(is)*tstep_real/60 !i.e. multiply by no. mins per timestep  !Is this correct?? Why not divide by nsh_real?
    endif

    ! Check value obtained is physically reasonable
    ! More water cannot drain than is in the surface state
    ! although high initial rate of drainage tries to drain more water than is in state within tstep
    ! May indicate shorter tstep needed, or a more suitable equation
    if (drain(is)>state(is)) then
      !write(*,*) 'Drainage:', is, drain(is), state(is), drain(is)-state(is), DrainEq, DrainCoef1, DrainCoef2, nsh_real
      call ErrorHint(61,'SUEWS_drain: drain(is) > state(is) for surface is ',drain(is),state(is),is)
      drain(is)=state(is)   !All water in state is drained (but no more)
    elseif(drain(is)<0.0001) then
      drain(is)=0
    endif
 endif

 return

 end subroutine drainage
!------------------------------------------------------------------------------   