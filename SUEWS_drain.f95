SUBROUTINE drainage(&

     ! input:
     is,&
     state_is,&
     StorCap,&
     DrainEq,&
     DrainCoef1,&
     DrainCoef2,&

     nsh_real,&
     ! output:
     drain_is)
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

  ! use allocateArray
  ! use gis_data
  ! use sues_data
  ! use time

  IMPLICIT NONE
  INTEGER,INTENT(in)::&
       is ! surface type number
  REAL (KIND(1d0)),INTENT(in)::&
       state_is,  &!Wetness status of surface type "is" [mm]
       StorCap,   &!current storage capacity [mm]
       DrainCoef1,&!Drainage coeff 1 [units depend on choice of eqn]
       DrainCoef2,&!Drainage coeff 2 [units depend on choice of eqn]
       DrainEq,   &!Drainage equation to use
       nsh_real    !nsh cast as a real for use in calculations
  REAL (KIND(1d0)),INTENT(out)::&
       drain_is!Drainage of surface type "is" [mm]


  !If surface is dry, no drainage occurs
  IF(state_is<0.000000001) THEN
     drain_is=0.0
  ELSE
     IF(INT(DrainEq)==1) THEN   !Falk and Niemczynowicz (1978): Drainage equation for paved, buildings and irrigated grass

        IF (state_is<StorCap) THEN
           drain_is=0   !No drainage if state is less than storage capacity
        ELSE
           drain_is=(DrainCoef1*(state_is-StorCap)**DrainCoef2)/nsh_real
        ENDIF

     ELSEIF(INT(DrainEq)==2) THEN   !Rutter eqn corrected for c=0, see Eq 9 of Calder & Wright 1986
        drain_is=(DrainCoef1*(EXP(DrainCoef2*state_is)-1))/nsh_real
        ! N.B. -1 is correct here but brackets are wrong in G&O 1991 Eq 5 & Ja11 Eq 18.

     ELSEIF(INT(DrainEq)==3) THEN   !Falk and Niemczynowicz (1978)
        drain_is=(DrainCoef1*(state_is**DrainCoef2))/nsh_real

        ! Option 4 removed by HCW 16 Feb 2015, as it is not used and appears to be problematic
        !elseif(int(DrainEq)==4) then    !Rutter eqn not corrected for c=0
        !   drain(is)=DrainCoef1*exp(DrainCoef2*(state(is)-StorCap))
        !   drain(is)=drain(is)*tstep_real/60 !i.e. multiply by no. mins per timestep  !Is this correct?? Why not divide by nsh_real?
     ENDIF

     ! Check value obtained is physically reasonable
     ! More water cannot drain than is in the surface state
     ! although high initial rate of drainage tries to drain more water than is in state within tstep
     ! May indicate shorter tstep needed, or a more suitable equation
     IF (drain_is>state_is) THEN
        !write(*,*) 'Drainage:', is, drain(is), state(is), drain(is)-state(is), DrainEq, DrainCoef1, DrainCoef2, nsh_real
        CALL ErrorHint(61,'SUEWS_drain: drain_is > state_is for surface is ',drain_is,state_is,is)
        drain_is=state_is   !All water in state is drained (but no more)
     ELSEIF(drain_is<0.0001) THEN
        drain_is=0
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE drainage
!------------------------------------------------------------------------------
