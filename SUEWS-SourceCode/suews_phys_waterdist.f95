MODULE WaterDist_module
   USE allocateArray, ONLY: nsurf, &
                            PavSurf, BldgSurf, &
                            ConifSurf, DecidSurf, GrassSurf, &
                            BSoilSurf, WaterSurf, ExcessSurf

   IMPLICIT NONE
   ! INTEGER, PARAMETER :: nsurf = 7
   ! INTEGER, PARAMETER :: PavSurf = 1
   ! INTEGER, PARAMETER :: BldgSurf = 2
   ! INTEGER, PARAMETER :: ConifSurf = 3
   ! INTEGER, PARAMETER :: DecidSurf = 4
   ! INTEGER, PARAMETER :: GrassSurf = 5
   ! INTEGER, PARAMETER :: BSoilSurf = 6
   ! INTEGER, PARAMETER :: WaterSurf = 7
   ! INTEGER, PARAMETER :: ExcessSurf = 8
CONTAINS

   !------------------------------------------------------------------------------
   SUBROUTINE drainage( &
      is, & !input
      state_is, &
      StorCap, &
      DrainEq, &
      DrainCoef1, &
      DrainCoef2, &
      nsh_real, &
      drain_is)!output

      !Calculation of drainage for each land surface.
      !INPUT: Storage capacity, type of drainage equation used, drainage coefficients
      !       used in the equation
      !Modified by HCW 16 Feb 2015
      !  Removed option of Eq 4 (calculation needs to be checked before re-implementing).
      !  Code writes an error if calculated drainage exceeds surface state_id (but code continues).
      !  This may indicate inappropriate drainage equation, storage capacities or model tstep.
      !Modified by LJ in Aug 2011. Drainage cannot exceed the surface storage.
      !Modified LJ in 10/2010
      !------------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(in):: is ! surface type number

      REAL(KIND(1d0)), INTENT(in)::state_is !Wetness status of surface type "is" [mm]
      REAL(KIND(1d0)), INTENT(in)::StorCap !current storage capacity [mm]
      REAL(KIND(1d0)), INTENT(in)::DrainCoef1 !Drainage coeff 1 [units depend on choice of eqn]
      REAL(KIND(1d0)), INTENT(in)::DrainCoef2 !Drainage coeff 2 [units depend on choice of eqn]
      REAL(KIND(1d0)), INTENT(in)::DrainEq !Drainage equation to use
      REAL(KIND(1d0)), INTENT(in)::nsh_real    !nsh cast as a real for use in calculations
      REAL(KIND(1d0)), INTENT(out):: drain_is!Drainage of surface type "is" [mm]

      !If surface is dry, no drainage occurs
      IF (state_is < 0.000000001) THEN
         drain_is = 0.0
      ELSE
         IF (INT(DrainEq) == 1) THEN   !Falk and Niemczynowicz (1978): Drainage equation for paved, buildings and irrigated grass

            IF (state_is < StorCap) THEN
               drain_is = 0   !No drainage if state_id is less than storage capacity
            ELSE
               drain_is = (DrainCoef1*(state_is - StorCap)**DrainCoef2)/nsh_real
            ENDIF

         ELSEIF (INT(DrainEq) == 2) THEN   !Rutter eqn corrected for c=0, see Eq 9 of Calder & Wright 1986
            drain_is = (DrainCoef1*(EXP(DrainCoef2*state_is) - 1))/nsh_real
            ! N.B. -1 is correct here but brackets are wrong in G&O 1991 Eq 5 & Ja11 Eq 18.

         ELSEIF (INT(DrainEq) == 3) THEN   !Falk and Niemczynowicz (1978)
            drain_is = (DrainCoef1*(state_is**DrainCoef2))/nsh_real

         ENDIF

         ! Check value obtained is physically reasonable
         ! More water cannot drain than is in the surface state_id
         ! although high initial rate of drainage tries to drain more water than is in state_id within tstep
         ! May indicate shorter tstep needed, or a more suitable equation
         IF (drain_is > state_is) THEN
            !write(*,*) 'Drainage:', is, drain(is), state_id(is), drain(is)-state_id(is), DrainEq, DrainCoef1, DrainCoef2, nsh_real
            CALL ErrorHint(61, 'SUEWS_drain: drain_is > state_is for surface is ', drain_is, state_is, is)
            drain_is = state_is   !All water in state_id is drained (but no more)
         ELSEIF (drain_is < 0.0001) THEN
            drain_is = 0
         ENDIF
      ENDIF

      RETURN

   END SUBROUTINE drainage
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   SUBROUTINE soilstore( &
      is, sfr, PipeCapacity, RunoffToWater, pin, & ! input:
      wu_EveTr, wu_DecTr, wu_Grass, drain, AddWater, addImpervious, nsh_real, stateOld, AddWaterRunoff, &
      PervFraction, addVeg, SoilStoreCap, addWaterBody, FlowChange, StateLimit, &
      runoffAGimpervious, surplusWaterBody, & ! inout:
      runoffAGveg, runoffPipes, ev, soilstore_id, SurplusEvap, runoffWaterBody, &
      p_mm, chang, runoff, state_id)!output:
      !------------------------------------------------------------------------------
      !Calculation of storage change
      ! LJ 27 Jan 2016
      !   -Removed tabs and cleaned the code
      ! HCW 08 Dec 2015
      !   -Added if-loop check for no Paved surfaces
      ! LJ 6 May 2015
      !   - Calculations of the piperunoff exceedings moved to separate subroutine updateFlood.
      !   - Now also called from snow subroutine
      !   - Evaporation is modified using EvapPart
      !   - when no water on impervious surfaces, evap occurs above pervious surfaces instead
      ! Rewritten by HCW 12 Feb 2015
      !   - Old variable 'p' for water input to the surface renamed to 'p_mm'
      !   - All water now added to p_mm first, before threshold checks or other calculations
      !   - Water from other grids now added to p_mm (instead of state_id for impervious surfaces)
      !   - Removed division of runoff by nsh, as whole model now runs at the same timestep
      !   - Adjusted transfer of ev between surfaces to conserve mass (not depth)
      !   - Volumes used for water transport between grids to account for SurfaceArea changing between grids
      !   - Added threshold check for state_id(WaterSurf) - was going negative
      ! Last modified HCW 09 Feb 2015
      !   - Removed StorCap input because it is provided by module allocateArray
      !   - Tidied and commented code
      ! Modified by LJ in November 2012:
      !   - P>10 was not taken into account for impervious surfaces - Was fixed.
      !   - Above impervious surfaces possibility of the state_id to exceed max capacity was limited
      !     although this should be possible - was fixed
      ! Modified by LJ 10/2010
      ! Rewritten mostly by LJ in 2010
      ! To do:
      !   - Finish area normalisation for RG2G & finish coding GridConnections
      !   - What is the 10 mm hr-1 threshold for?
      !  - Decide upon and correct storage capacities here & in evap subroutine
      !  - FlowChange units should be mm hr-1 - need to update everywhere
      !   - Add SurfaceFlood(is)?
      !   - What happens if sfr(is) = 0 or 1?
      !   - Consider how irrigated trees actually works...
      !------------------------------------------------------------------------------

      IMPLICIT NONE

      !Stores flood water when surface state_id exceeds storage capacity [mm]
      INTEGER, INTENT(in)::is ! surface type

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr! surface fractions
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::AddWater!Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::stateOld!Wetness status of each surface type from previous timestep [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::AddWaterRunoff!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::SoilStoreCap!Capacity of soil store for each surface [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::StateLimit!Limit for state_id of each surface type [mm] (specified in input files)

      REAL(KIND(1d0)), INTENT(in)::PipeCapacity!Capacity of pipes to transfer water
      REAL(KIND(1d0)), INTENT(in)::RunoffToWater!Fraction of surface runoff going to water body
      REAL(KIND(1d0)), INTENT(in)::pin!Rain per time interval
      REAL(KIND(1d0)), INTENT(in)::wu_EveTr!Water use for evergreen trees/shrubs [mm]
      REAL(KIND(1d0)), INTENT(in)::wu_DecTr!Water use for deciduous trees/shrubs [mm]
      REAL(KIND(1d0)), INTENT(in)::wu_Grass!Water use for grass [mm]
      REAL(KIND(1d0)), INTENT(in)::addImpervious!Water from impervious surfaces of other grids [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(in)::nsh_real!nsh cast as a real for use in calculations
      REAL(KIND(1d0)), INTENT(in)::PervFraction! sum of surface cover fractions for impervious surfaces
      REAL(KIND(1d0)), INTENT(in)::addVeg!Water from vegetated surfaces of other grids [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(in)::addWaterBody!Water from water surface of other grids [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(in)::FlowChange!Difference between the input and output flow in the water body

      REAL(KIND(1d0)), INTENT(inout)::runoffAGimpervious!Above ground runoff from impervious surface [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(inout)::surplusWaterBody!Extra runoff that goes to water body [mm] as specified by RunoffToWater
      REAL(KIND(1d0)), INTENT(inout)::runoffAGveg!Above ground runoff from vegetated surfaces [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(inout)::runoffPipes!Runoff in pipes [mm] for whole surface area
      REAL(KIND(1d0)), INTENT(inout)::ev!Evaporation
      REAL(KIND(1d0)), INTENT(inout)::runoffWaterBody!Above ground runoff from water surface [mm] for whole surface area

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout)::soilstore_id  !Soil moisture of each surface type [mm]
      REAL(KIND(1d0)), DIMENSION(2), INTENT(inout)    ::SurplusEvap!Surplus for evaporation in 5 min timestep

      REAL(KIND(1d0)), INTENT(out)::p_mm!Inputs to surface water balance

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::chang !Change in state_id [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::runoff!Runoff from each surface type [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in) ::drain !Drainage of each surface type [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::state_id !Wetness status of each surface type [mm]

      !Extra evaporation [mm] from impervious surfaces which cannot happen due to lack of water
      REAL(KIND(1d0)):: EvPart
      REAL(KIND(1d0)), PARAMETER:: NotUsed = -55.5
      REAL(KIND(1d0)), PARAMETER:: IPThreshold_mmhr = 10 ! NB:this should be an input and can be specified. SG 25 Apr 2018

      !Initialise extra evaporation to zero
      EvPart = 0

      !SurfaceFlood(is) = 0 !!This probably needs to be carried over between timesteps, but reset for now

      !==================================================================
      ! Combine water inputs to the current surface
      ! Add external water use for each surface type
      SELECT CASE (is)
      CASE (ConifSurf)
         p_mm = pin + wu_EveTr
      CASE (DecidSurf)
         p_mm = pin + wu_DecTr
      CASE (GrassSurf)
         p_mm = pin + wu_Grass
      CASE default
         p_mm = pin
      END SELECT

      ! Add water from other surfaces within the same grid (RS2S) ----
      ! AddWater is the water supplied to the current surface from other surfaces
      !  i.e. drain*WaterDist (see SUEWS_ReDistributeWater)
      p_mm = p_mm + AddWater(is)
      !==================================================================

      !========surface-specific calculation=========================
      SELECT CASE (is)
      CASE (PavSurf, BldgSurf)
         !==== Impervious surfaces (Paved, Buildings) ======================
         ! Add water from neighbouring grids (RG2G)
         ! Add to PavSurf only, as water cannot flow onto buildings
         IF (is == PavSurf) THEN
            IF (sfr(PavSurf) /= 0) THEN   ! If loop added HCW 08 Dec 2015
               p_mm = p_mm + addImpervious/sfr(PavSurf)
            ENDIF
         ENDIF

         ! Calculate change in surface state_id (inputs - outputs)
         chang(is) = p_mm - (drain(is) + ev)

         ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
         ! and does not affect state_id
         IF (p_mm > IPThreshold_mmhr/nsh_real) THEN
            runoff(is) = runoff(is) + (p_mm - IPThreshold_mmhr/nsh_real)
            chang(is) = IPThreshold_mmhr/nsh_real - (drain(is) + ev)
         ENDIF

         ! Calculate updated state_id using chang
         state_id(is) = stateOld(is) + chang(is)

         ! Check state_id is within physical limits between zero (dry) and max. storage capacity
         IF (state_id(is) < 0.0) THEN   ! Cannot have a negative surface state_id
            ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
            ! Allow evaporation only until surface is dry (state_id(is)=0); additional evaporation -> evaporation surplus
            SurplusEvap(is) = ABS(state_id(is))   !Surplus evaporation is that which tries to evaporate non-existent water
            ev = ev - SurplusEvap(is)          !Limit evaporation according to water availability
            state_id(is) = 0.0                    !Now surface is dry
            ! elseif (state_id(is)>StoreDrainPrm(6,is)) then   !!This should perhaps be StateLimit(is)
            !    !! If state_id exceeds the storage capacity, then the excess goes to surface flooding
            !    !SurfaceFlood(is)=SurfaceFlood(is)+(state_id(is)-StoreDrainPrm(6,is))   !!Need to deal with this properly
            !    runoff(is)=runoff(is)+(state_id(is)-StoreDrainPrm(6,is))   !!needs to go to flooding
            !    state_id(is)=StoreDrainPrm(6,is)              !Now surface state_id is at max (storage) capacity
         ENDIF

         ! Recalculate change in surface state_id from difference with previous timestep
         chang(is) = state_id(is) - stateOld(is)

         ! Runoff -------------------------------------------------------
         ! For impervious surfaces, some of drain(is) becomes runoff
         runoff(is) = runoff(is) + drain(is)*AddWaterRunoff(is)   !Drainage (that is not flowing to other surfaces) goes to runoff

         !So, up to this point, runoff(is) can have contributions if
         ! p_mm > ipthreshold (water input too fast)
         ! state_id > StoreDrainPrm(6,is) (net water exceeds storage capacity)
         ! WaterDist specifies some fraction of drain(is) -> runoff

      CASE (ConifSurf:BSoilSurf)
         !==== For Conif, Decid, Grass, BSoil surfaces ==================
         ! Transfer evaporation surplus from impervious surfaces to pervious surfaces
         EvPart = MERGE( &
                  DOT_PRODUCT(SurplusEvap(PavSurf:BldgSurf), sfr(PavSurf:BldgSurf)/PervFraction), &
                  0d0, &
                  PervFraction /= 0)

         ! Add surplus evaporation to ev for pervious surfaces
         ev = ev + EvPart

         ! ---- Add water from neighbouring grids (RG2G) ----
         ! Add to Grass and BSoil only, as water cannot flow onto trees
         IF (is == GrassSurf .OR. is == BSoilSurf) THEN
            IF ((sfr(GrassSurf) + sfr(BSoilSurf)) /= 0) THEN
               p_mm = p_mm + addVeg/(sfr(GrassSurf) + sfr(BSoilSurf))
            ENDIF
         ENDIF

         ! Calculate change in surface state_id (inputs - outputs)
         chang(is) = p_mm - (drain(is) + ev)

         ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
         !  and does not affect state_id
         IF (p_mm > IPThreshold_mmhr/nsh_real) THEN
            runoff(is) = runoff(is) + (p_mm - IPThreshold_mmhr/nsh_real)
            chang(is) = IPThreshold_mmhr/nsh_real - (drain(is) + ev)
         ENDIF

         ! Calculate updated state_id using chang
         state_id(is) = stateOld(is) + chang(is)

         ! Check state_id is within physical limits between zero (dry) and max. storage capacity
         IF (state_id(is) < 0.0) THEN   ! Cannot have a negative surface state_id
            ! If there is not sufficient water on the surface, then remove water from soilstore
            ! Allow evaporation until soilstore_id is depleted and surface is dry
            IF ((soilstore_id(is) + state_id(is)) >= 0) THEN
               soilstore_id(is) = soilstore_id(is) + state_id(is)
               state_id(is) = 0.0
               ! If there is not sufficient water on the surface or soilstore, then don't allow this evaporation to happen
            ELSE
               ev = ev - ABS(state_id(is))   !Limit evaporation according to water availability
               state_id(is) = 0.0          !Now surface is dry
            ENDIF

            !elseif (state_id(is)>StoreDrainPrm(6,is)) then   !!This should perhaps be StateLimit(is)
            !   !! If state_id exceeds the storage capacity, then the excess goes to surface flooding
            !   !SurfaceFlood(is)=SurfaceFlood(is)+(state_id(is)-StoreDrainPrm(6,is))   !!Need to deal with this properly
            !   runoff(is)=runoff(is)+(state_id(is)-StoreDrainPrm(6,is))   !!needs to go to flooding
            !   state_id(is)=StoreDrainPrm(6,is)              !Now surface state_id is at max (storage) capacity
         ENDIF

         ! Recalculate change in surface state_id from difference with previous timestep
         chang(is) = state_id(is) - stateOld(is)

         !Where should this go? Used to be before previous part!!
         ! soilstore_id -------------------------------------------------
         ! For pervious surfaces (not water), some of drain(is) goes to soil storage
         ! Drainage (that is not flowing to other surfaces) goes to soil storages
         soilstore_id(is) = soilstore_id(is) + drain(is)*AddWaterRunoff(is)

         ! If soilstore is full, the excess will go to runoff
         IF (soilstore_id(is) > SoilStoreCap(is)) THEN  ! TODO: this should also go to flooding of some sort
            runoff(is) = runoff(is) + (soilstore_id(is) - SoilStoreCap(is))
            soilstore_id(is) = SoilStoreCap(is)
         ELSEIF (soilstore_id(is) < 0) THEN   !! QUESTION: But where does this lack of water go? !!Can this really happen here?
            CALL ErrorHint(62, 'SUEWS_store: soilstore_id(is) < 0 ', soilstore_id(is), NotUsed, is)
            ! Code this properly - soilstore_id(is) < 0 shouldn't happen given the above loops
            !soilstore_id(is)=0   !Groundwater / deeper soil should kick in
         ENDIF

      CASE (WaterSurf)
         IF (sfr(WaterSurf) /= 0) THEN

            ! ---- Add water from neighbouring grids (RG2G) ----
            p_mm = p_mm + addWaterBody/sfr(WaterSurf)

            ! Calculate change in surface state_id (inputs - outputs)
            ! No drainage for water surface
            ! FlowChange is the difference in input and output flows [mm hr-1]
            chang(is) = p_mm + FlowChange/nsh_real - ev

            ! Calculate updated state_id using chang
            state_id(is) = stateOld(is) + chang(is)

            ! Check state_id is within physical limits between zero (dry) and max. storage capacity
            IF (state_id(is) < 0.0) THEN   ! Cannot have a negative surface state_id
               ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
               ev = ev - ABS(state_id(is))   !Limit evaporation according to water availability
               state_id(is) = 0.0          !Now surface is dry
               !elseif (state_id(is)>StoreDrainPrm(6,is)) then   !!This should perhaps be StateLimit(is)
               !   !! If state_id exceeds the storage capacity, then the excess goes to surface flooding
               !   !SurfaceFlood(is)=SurfaceFlood(is)+(state_id(is)-StoreDrainPrm(6,is))   !!Need to deal with this properly
               !   runoff(is)=runoff(is)+(state_id(is)-StoreDrainPrm(6,is))   !!needs to go to flooding
               !   state_id(is)=StoreDrainPrm(6,is)              !Now surface state_id is at max (storage) capacity
            ENDIF

            ! Recalculate change in surface state_id from difference with previous timestep
            chang(is) = state_id(is) - stateOld(is)

            ! If state_id exceeds limit, then excess goes to runoff (currently applies to water StoreDrainPrm only)
            IF (state_id(WaterSurf) > StateLimit(WaterSurf)) THEN
               runoff(WaterSurf) = runoff(WaterSurf) + (state_id(WaterSurf) - StateLimit(WaterSurf))
               state_id(WaterSurf) = StateLimit(WaterSurf)
               runoffWaterBody = runoffWaterBody + runoff(WaterSurf)*sfr(WaterSurf)
            ELSE
               state_id(WaterSurf) = state_id(WaterSurf) + surplusWaterBody
               IF (state_id(WaterSurf) > StateLimit(WaterSurf)) THEN
                  runoffWaterBody = runoffWaterBody + (state_id(WaterSurf) - StateLimit(WaterSurf))*sfr(WaterSurf)
                  state_id(WaterSurf) = StateLimit(WaterSurf)
               ENDIF
            ENDIF

            ! Recalculate change in surface state_id from difference with previous timestep
            chang(is) = state_id(is) - stateOld(is)
         ENDIF
      END SELECT
      !==================================================================

      !==== RUNOFF ======================================================
      ! TODO: to consider areas here - SurfaceArea may vary between grids too
      ! - also implement where water for next surface is calculated (RunoffFromGrid subroutine)
      ! Calculations of the piperunoff exceedensances moved to separate subroutine so that from snow same
      ! calculations can be made. LJ in May 2015

      IF (is < WaterSurf) THEN   !Not for water body
         !  CALL updateFlood
         CALL updateFlood( &
            is, runoff, &! input:
            sfr, PipeCapacity, RunoffToWater, &
            runoffAGimpervious, surplusWaterBody, runoffAGveg, runoffPipes)! inout:
      ENDIF

   END SUBROUTINE soilstore
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   SUBROUTINE updateFlood( &
      is, runoff, &! input:
      sfr, PipeCapacity, RunoffToWater, &
      runoffAGimpervious, surplusWaterBody, runoffAGveg, runoffPipes)! inout:

      IMPLICIT NONE

      INTEGER, INTENT(in) :: is
      REAL(KIND(1d0)), INTENT(in) :: sfr(nsurf), runoff(nsurf), PipeCapacity, RunoffToWater
      REAL(KIND(1d0)), INTENT(inout) :: runoffAGimpervious, surplusWaterBody, runoffAGveg, runoffPipes

      ! Add runoff to pipes
      runoffPipes = runoffPipes + (runoff(is)*sfr(is))

      ! If pipe capacity is full, surface runoff occurs
      ! N.B. this will happen each loop (replicates pipes filling up)
      IF (runoffPipes > PipeCapacity) THEN

         !------Paved and building surface
         IF (is == PavSurf .OR. is == BldgSurf) THEN
            IF (sfr(WaterSurf) > 0.0000001) THEN
               ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
               ! RunoffToWater is specified in SUEWS_SiteSelect.txt
               runoffAGimpervious = runoffAGimpervious + (runoffPipes - PipeCapacity)*(1 - RunoffToWater)
               surplusWaterBody = surplusWaterBody + (runoffPipes - PipeCapacity)*RunoffToWater
            ELSE
               ! Otherwise, all flood water must go to runoff
               runoffAGimpervious = runoffAGimpervious + (runoffPipes - PipeCapacity)
            ENDIF
            !------other surfaces
         ELSEIF (is >= ConifSurf .AND. is <= BSoilSurf) THEN
            IF (sfr(WaterSurf) > 0.0000001) THEN
               ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
               runoffAGveg = runoffAGveg + (runoffPipes - PipeCapacity)*(1 - RunoffToWater)
               surplusWaterBody = surplusWaterBody + (runoffPipes - PipeCapacity)*RunoffToWater
            ELSE
               ! Otherwise, all flood water must go to runoff
               runoffAGveg = runoffAGveg + (runoffPipes - PipeCapacity)
            ENDIF
         ENDIF

         runoffPipes = PipeCapacity   !Pipes are at their max capacity

      ENDIF   !If runoff exceed pipe capacity

   END SUBROUTINE updateFlood
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   SUBROUTINE ReDistributeWater( &
      snowUse, WaterDist, sfr, Drain, &! input:
      AddWaterRunoff, AddWater)! output:
      !Drainage moves into different parts defined by WaterDistSS_YYYY.txt. LJ 2010
      !AddWater(is) is that amount of water that is gained for each surface
      !Latest update takes snow into account. 22/03/2013 LJ
      !-------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(in)::snowUse!Snow part used (1) or not used (0)

      REAL(KIND(1d0)), INTENT(in)::WaterDist(nsurf + 1, nsurf - 1) !Within-grid water distribution to other surfaces and runoff/soil store [-]
      REAL(KIND(1d0)), INTENT(in)::sfr(nsurf)                !Surface fractions [-]
      REAL(KIND(1d0)), INTENT(in)::Drain(nsurf)               !Drainage of each surface type [mm]

      REAL(KIND(1d0)), INTENT(out)::AddWaterRunoff(nsurf)!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
      REAL(KIND(1d0)), INTENT(out)::AddWater(nsurf)        !Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]

      INTEGER::ii, jj
      INTEGER::NSurfDoNotReceiveDrainage = 0!Number of surfaces that do not receive drainage water (green roof)

      !Fractions that go to runoff from each surface
      DO ii = 1, nsurf - 1   !not water in the calculation
         AddWaterRunoff(ii) = WaterDist(8, ii)
      ENDDO
      AddWaterRunoff(WaterSurf) = 0
      AddWater = 0

      DO ii = 1, nsurf - NSurfDoNotReceiveDrainage !go through surfaces from 1 to 7. These gain water through drainage
         DO jj = 1, nsurf - (NSurfDoNotReceiveDrainage + 1) !From where surface ii can gain water - can't gain water from itself

            IF (sfr(ii) /= 0) THEN !Water movement takes place only if surface fraction exists

               !No snow calculations!
               IF (snowUse == 0) THEN
                  AddWater(ii) = AddWater(ii) + (Drain(jj)*sfr(jj)/sfr(ii))*WaterDist(ii, jj) !Original

                  !Snow included, This needs to be fixed at some point. LJ Mar 2013
               ELSE
                  AddWaterRunoff(jj) = AddWaterRunoff(jj) + WaterDist(ii, jj) !No receiving surface -> runoff
               ENDIF

            ELSE
               AddWaterRunoff(jj) = AddWaterRunoff(jj) + WaterDist(ii, jj) !If no receiving surface exists,
               !water fraction goes to AddWaterRunoff
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE ReDistributeWater
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   SUBROUTINE SUEWS_update_SoilMoist( &
      NonWaterFraction, &!input
      SoilStoreCap, sfr, soilstore_id, &
      SoilMoistCap, SoilState, &!output
      vsmd, smd)
      IMPLICIT NONE

      ! INTEGER,INTENT(in)::nsurf,ConifSurf,DecidSurf,GrassSurf
      REAL(KIND(1d0)), INTENT(in)::NonWaterFraction
      REAL(KIND(1d0)), INTENT(in), DIMENSION(nsurf)::SoilStoreCap, sfr, soilstore_id

      REAL(KIND(1d0)), INTENT(out)::SoilMoistCap, SoilState
      REAL(KIND(1d0)), INTENT(out)::vsmd, smd

      INTEGER :: is

      SoilMoistCap = 0   !Maximum capacity of soil store [mm] for whole surface
      SoilState = 0      !Area-averaged soil moisture [mm] for whole surface

      IF (NonWaterFraction /= 0) THEN !Soil states only calculated if soil exists. LJ June 2017
         DO is = 1, nsurf - 1   !No water body included
            SoilMoistCap = SoilMoistCap + (SoilStoreCap(is)*sfr(is)/NonWaterFraction)
            SoilState = SoilState + (soilstore_id(is)*sfr(is)/NonWaterFraction)
         ENDDO
      ENDIF

      !If loop removed HCW 26 Feb 2015
      !if (ir==1) then  !Calculate initial smd
      smd = SoilMoistCap - SoilState
      !endif

      ! Calculate soil moisture for vegetated surfaces only (for use in surface conductance)
      vsmd = 0
      DO is = ConifSurf, GrassSurf  !Vegetated surfaces only
         IF (sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf) == 0) THEN
            vsmd = 0
         ELSE
            vsmd = vsmd + (SoilStoreCap(is) - soilstore_id(is))*sfr(is)/(sfr(ConifSurf) + sfr(DecidSurf) + sfr(GrassSurf))
         END IF
         !write(*,*) is, vsmd, smd
      ENDDO

   END SUBROUTINE SUEWS_update_SoilMoist
   !------------------------------------------------------------------------------

   !========== Calculate soil moisture ============
   SUBROUTINE SUEWS_cal_SoilState( &
      SMDMethod, xsmd, NonWaterFraction, SoilMoistCap, &!input
      SoilStoreCap, surf_chang_per_tstep, &
      soilstore_id, soilstoreOld, sfr, &
      smd, smd_nsurf, tot_chang_per_tstep, SoilState)!output

      IMPLICIT NONE
      INTEGER, PARAMETER :: nsurf = 7

      INTEGER, INTENT(in) ::SMDMethod
      REAL(KIND(1d0)), INTENT(in)::xsmd
      REAL(KIND(1d0)), INTENT(in)::NonWaterFraction
      REAL(KIND(1d0)), INTENT(in)::SoilMoistCap

      REAL(KIND(1d0)), INTENT(in)::surf_chang_per_tstep
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::soilstore_id !Soil moisture of each surface type [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::soilstoreOld !Soil moisture of each surface type from previous timestep [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::sfr
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(in)::SoilStoreCap        !Capacity of soil store for each surface [mm]

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(out)::smd_nsurf !smd for each surface
      REAL(KIND(1d0)), INTENT(out)::SoilState !Area-averaged soil moisture [mm] for whole surface
      REAL(KIND(1d0)), INTENT(out)::smd !One value for whole surface
      REAL(KIND(1d0)), INTENT(out)::tot_chang_per_tstep !Change in surface state_id

      REAL(KIND(1d0)), PARAMETER::NotUsed = -999
      REAL(KIND(1d0)), PARAMETER::NAN = -999
      INTEGER :: is

      SoilState = 0       !Area-averaged soil moisture [mm] for whole surface
      IF (NonWaterFraction /= 0) THEN !Fixed for water surfaces only
         DO is = 1, nsurf - 1   !No water body included
            SoilState = SoilState + (soilstore_id(is)*sfr(is)/NonWaterFraction)
            IF (SoilState < 0) THEN
               CALL ErrorHint(62, 'SUEWS_Calculations: total SoilState < 0 (just added surface is) ', SoilState, NotUsed, is)
            ELSEIF (SoilState > SoilMoistCap) THEN
               CALL ErrorHint(62, 'SUEWS_Calculations: total SoilState > capacity (just added surface is) ', SoilState, NotUsed, is)
               !SoilMoist_state=SoilMoistCap !What is this LJ 10/2010 - QUESTION: SM exceeds capacity, but where does extra go?HCW 11/2014
            ENDIF
         ENDDO  !end loop over surfaces
         ! SoilState = DOT_PRODUCT(soilstore_id(1:nsurf - 1), sfr(1:nsurf - 1))/NonWaterFraction
         ! IF (SoilState < 0) THEN
         !    CALL ErrorHint(62, 'SUEWS_Calculations: total SoilState < 0 (just added surface is) ', SoilState, NotUsed, is)
         ! ELSEIF (SoilState > SoilMoistCap) THEN
         !    CALL ErrorHint(62, 'SUEWS_Calculations: total SoilState > capacity (just added surface is) ', SoilState, NotUsed, is)
         !    !SoilMoist_state=SoilMoistCap !What is this LJ 10/2010 - QUESTION: SM exceeds capacity, but where does extra go?HCW 11/2014
         ! ENDIF
      ENDIF

      ! Calculate soil moisture deficit
      smd = SoilMoistCap - SoilState   !One value for whole surface
      smd_nsurf = SoilStoreCap - soilstore_id   !smd for each surface

      ! Soil stores can change after horizontal water movements
      ! Calculate total change in surface and soil state_id
      tot_chang_per_tstep = surf_chang_per_tstep   !Change in surface state_id
      DO is = 1, (nsurf - 1)   !No soil for water surface (so change in soil moisture is zero)
         tot_chang_per_tstep = tot_chang_per_tstep + ((soilstore_id(is) - soilstoreOld(is))*sfr(is))   !Add change in soil state_id
      ENDDO

      IF (SMDMethod > 0) THEN ! use observed value
         !  smd_nsurf=NAN
         smd_nsurf = NAN
         smd = xsmd
      ENDIF

   END SUBROUTINE SUEWS_cal_SoilState
   !===================================================================================

   SUBROUTINE SUEWS_cal_HorizontalSoilWater( &
      sfr, &! input: ! surface fractions
      SoilStoreCap, &!Capacity of soil store for each surface [mm]
      SoilDepth, &!Depth of sub-surface soil store for each surface [mm]
      SatHydraulicConduct, &!Saturated hydraulic conductivity for each soil subsurface [mm s-1]
      SurfaceArea, &!Surface area of the study area [m2]
      NonWaterFraction, &! sum of surface cover fractions for all except water surfaces
      tstep_real, & !tstep cast as a real for use in calculations
      soilstore_id, &! inout: !Soil moisture of each surface type [mm]
      runoffSoil, &!Soil runoff from each soil sub-surface [mm]
      runoffSoil_per_tstep &!  output:!Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)
      )
      !Transfers water in soil stores of land surfaces LJ (2010)
      !Change the model to use varying hydraulic conductivity instead of constant value LJ (7/2011)
      !If one of the surface's soildepth is zero, no water movement is considered
      ! LJ  15/06/2017 Modification:   - Moved location of runoffSoil_per_tstep within previous if-loop to avoid dividing with zero with 100% water surface
      ! HCW 22/02/2017 Modifications:  - Minor bug fixed in VWC1/B_r1 comparison - if statements reversed
      ! HCW 13/08/2014 Modifications:  - Order of surfaces reversed (for both is and jj loops)
      !                                - Number of units (e.g. properties) added to distance calculation
      ! HCW 12/08/2014 Modifications:  - Distance changed from m to mm in dI_dt calculation
      !                                - dI_dt [mm s-1] multiplied by no. seconds in timestep -> dI [mm]
      !                                - if MatPot is set to max. value (100000 mm), Km set to 0 mm s-1
      !                                - Provide parameters for residual volumetric soil moisture [m3 m-3]
      !                                   (currently hard coded as 0.1 m3 m-3 for testing)
      !
      !------------------------------------------------------
      ! use SUES_data
      ! use gis_data
      ! use time
      ! use allocateArray

      IMPLICIT NONE

      REAL(KIND(1d0)), INTENT(in) ::sfr(nsurf)! surface fractions
      REAL(KIND(1d0)), INTENT(in) ::SoilStoreCap(nsurf)!Capacity of soil store for each surface [mm]
      REAL(KIND(1d0)), INTENT(in) ::SoilDepth(nsurf)!Depth of sub-surface soil store for each surface [mm]
      REAL(KIND(1d0)), INTENT(in) ::SatHydraulicConduct(nsurf)!Saturated hydraulic conductivity for each soil subsurface [mm s-1]
      REAL(KIND(1d0)), INTENT(in) ::SurfaceArea!Surface area of the study area [m2]
      REAL(KIND(1d0)), INTENT(in) ::NonWaterFraction! sum of surface cover fractions for all except water surfaces
      REAL(KIND(1d0)), INTENT(in) ::tstep_real !tstep cast as a real for use in calculations

      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout) ::soilstore_id!Soil moisture of each surface type [mm]
      REAL(KIND(1d0)), DIMENSION(nsurf), INTENT(inout) ::runoffSoil!Soil runoff from each soil sub-surface [mm]

      REAL(KIND(1d0)), INTENT(out) :: runoffSoil_per_tstep!Runoff to deep soil per timestep [mm] (for whole surface, excluding water body)

      INTEGER::jj, is
      REAL(KIND(1d0)):: &
         DimenWaterCon1, DimenWaterCon2, &
         SoilMoistCap_Vol1, &
         SoilMoist_vol1, &
         SoilMoistCap_Vol2, &
         SoilMoist_vol2, &
         B_r1, MatPot1, Km1, &
         B_r2, MatPot2, Km2, &
         Distance, KmWeight, dI, &
         dI_dt!Water flow between two stores

      REAL(KIND(1d0)), PARAMETER:: &
         alphavG = 0.0005, &  !Set alphavG to match value in van Genuchten (1980) [mm-1]
         NUnits = 1   !Can change to represent plot/base unit size

      ! SoilMoist_vol1,2     = Volumetric soil moisture [m3 m-3]
      ! SoilMoistCap_vol1,2  = Volumetric soil moisture capacity [m3 m-3] (from FunctionalTypes)
      ! MatPot1,2            = Water potential (i.e. pressure head) of store [mm]
      ! DimenWaterCon1,2     = Dimensionless water content, or relative saturation [-]
      ! Distance             = Distance between two stores [m]
      ! B_r1,2               = Residual volumetric soil moisture [m3 m-3]
      ! Km1,2                = Hydraulic conductivity of store [mm s-1]
      ! KmWeight             = Weighted hydraulic conductivity [mm s-1]
      ! alphavG              = Parameter (could depend on soil texture) [mm-1]
      ! dI                   = Water flow between stores [mm] dI = dI_dt * no. secs in each timestep
      !                         if dI > 0, first surface gains water, second surface loses water
      ! NUnits               = Number of repeating units (e.g. properties, blocks) for distance calculation [-]

      runoffSoil_per_tstep = 0

      DO is = 1, nsurf - 1 !nsurf-1,1,-1  !Loop through each surface, excluding water surface (runs backwards as of 13/08/2014, HCW)

         IF (sfr(is) /= 0 .AND. SoilStoreCap(is) > 0) THEN  !If particular surface area exists
            ! and is capable of storing water (SoilStoreCap [mm])
            DO jj = is + 1, nsurf - 1 !is-1,1,-1  !Sub-loop through remaining surfaces (runs backwards as of 13/08/2014, HCW)

               IF (sfr(jj) /= 0 .AND. SoilStoreCap(jj) > 0) THEN  !If other surface area exists
                  ! and is capable of storing water

                  ! ---- For surface 1 -----------------------------------------------------
                  ! Calculate non-saturated VWC
                  SoilMoistCap_Vol1 = SoilStoreCap(is)/SoilDepth(is) !Volumetric soil moisture capacity [m3 m-3] (i.e. saturated VWC)
                  SoilMoist_vol1 = soilstore_id(is)/SoilDepth(is) !Volumetric soil moisture [m3 m-3]

                  !B_r1=SoilMoistCap_Vol1-SoilMoist_vol1  !Residual soil moisture content [m3 m-3]
                  B_r1 = 0.1 !HCW 12/08/2014 Temporary fix
                  ! Need to add residual soil moisture values to FunctionalTypes
                  !B_r1=VolSoilMoistRes(is) !Residual soil moisture content [m3 m-3]

                  !Order of if statements reversed HCW 22 Feb 2017
                  !If soil moisture less than or equal to residual value, set MatPot to max and Km to 0 to suppress water movement
                  IF (B_r1 >= SoilMoist_vol1) THEN
                     MatPot1 = 100000
                     Km1 = 0 !Added by LJ in Nov 2013
                     ! Otherwise, there should be enough water in the soil to allow horizontal transfer
                  ELSE
                     DimenWaterCon1 = (SoilMoist_vol1 - B_r1)/(SoilMoistCap_Vol1 - B_r1) !Dimensionless water content [-]

                     ! If very large or very small, adjust for calculation of MatPot and Km
                     IF (DimenWaterCon1 > 0.99999) THEN
                        DimenWaterCon1 = DimenWaterCon1 - 0.0001 !This cannot equal 1
                     ENDIF

                     IF (DimenWaterCon1 < 0.00000005) THEN
                        DimenWaterCon1 = DimenWaterCon1 + 0.0000001   !Added HCW 22 Feb 2017
                     ENDIF

                     !van Genuchten (1980), with n=2 and m = 1-1/n = 1/2
                     !Water potential of first store [mm] (van Genuchten 1980, Eq 3 rearranged)
                     MatPot1 = SQRT(1/DimenWaterCon1**2 - 1)/alphavG

                     !Hydraulic conductivity of first store [mm s-1] (van Genuchten 1980, Eq 8)
                     Km1 = SatHydraulicConduct(is)*SQRT(DimenWaterCon1)*(1 - (1 - DimenWaterCon1**2)**0.5)**2

                     !Check this value (HCW 12/08/2014)
                     IF (MatPot1 > 100000) THEN
                        MatPot1 = 100000  !Max. potential is 100000 mm (van Genuchten 1980)
                        Km1 = 0   !Added by HCW 12/08/2014
                     ENDIF

                  ENDIF

                  ! ---- For surface 2 -----------------------------------------------------
                  ! Calculate non-saturated VWC
                  SoilMoistCap_Vol2 = SoilStoreCap(jj)/SoilDepth(jj) !Volumetric soil moisture capacity [m3 m-3] (i.e. saturated VWC)
                  SoilMoist_vol2 = soilstore_id(jj)/SoilDepth(jj) !Volumetric soil moisture [m3 m-3]

                  !B_r2=SoilMoistCap_Vol2-SoilMoist_vol2  !Residual soil moisture content [m3 m-3]
                  B_r2 = 0.1 !HCW 12/08/2014 Temporary fix
                  ! Need to add residual soil moisture values to FunctionalTypes
                  !B_r2=VolSoilMoistRes(jj) !Residual soil moisture content [m3 m-3]

                  !If soil moisture below residual value, set MatPot to maximum
                  IF (B_r2 >= SoilMoist_vol2) THEN
                     MatPot2 = 100000
                     Km2 = 0 !Added by LJ in Nov 2013
                  ELSE
                     DimenWaterCon2 = (SoilMoist_vol2 - B_r2)/(SoilMoistCap_Vol2 - B_r2) !Dimensionless water content [-]

                     IF (DimenWaterCon2 > 0.99999) THEN
                        DimenWaterCon2 = DimenWaterCon2 - 0.0001 !This cannot equal 1
                     ENDIF

                     IF (DimenWaterCon2 < 0.00000005) THEN
                        DimenWaterCon2 = DimenWaterCon2 + 0.0000001   !Added HCW 22 Feb 2017
                     ENDIF

                     !van Genuchten (1980), with n=2 and m = 1-1/n = 1/2
                     !Water potential of second store [mm] (van Genuchten 1980, Eq 3 rearranged)
                     MatPot2 = SQRT(1/DimenWaterCon2**2 - 1)/alphavG

                     !Hydraulic conductivity of second store [mm s-1] (van Genuchten 1980, Eq 8)
                     Km2 = SatHydraulicConduct(jj)*SQRT(DimenWaterCon2)*(1 - (1 - DimenWaterCon2**2)**0.5)**2

                     IF ((MatPot2) > 100000) THEN
                        MatPot2 = 100000 !Max. potential is 100000 mm (van Genuchten 1980)
                        Km2 = 0   !Added by HCW 12/08/2014
                     ENDIF

                  ENDIF

                  ! ------------------------------------------------------------------------

                  !Find distance between the two stores (see Jarvi et al. 2011)
                  !SurfaceArea in m2 (changed from ha to m2 n SUEWS_Initial), so Distance in m
                  Distance = (SQRT(sfr(is)*SurfaceArea/NUnits) + SQRT(sfr(jj)*SurfaceArea/NUnits))/2

                  !Calculate areally-weighted hydraulic conductivity [mm s-1]
                  KmWeight = (sfr(is)*Km1 + sfr(jj)*Km2)/(sfr(is) + sfr(jj))

                  !Find water flow between the two stores [mm s-1] (Green-Ampt equation, Hillel 1971)
                  !Multiply Distance by 1000 to convert m to mm (HCW 12/08/2014)
                  dI_dt = -(KmWeight)*(-MatPot1 + MatPot2)/(Distance*1000)

                  !Multiply dI_dt by number of seconds in timestep to convert mm s-1 to mm
                  !Added by HCW 12/08/2014
                  dI = dI_dt*tstep_real  !Use dI instead of dI_dt in the following calculations

                  !Move water (in mm) ------------------------------------------------------
                  !Water moves only if (i) there is sufficient water to move and (ii) there is space to move it

                  ! If there is sufficient water in both surfaces, allow movement of dI to occur
                  IF ((soilstore_id(jj) >= dI*sfr(is)/sfr(jj)) .AND. ((soilstore_id(is) + dI) >= 0)) THEN
                     soilstore_id(is) = soilstore_id(is) + dI
                     soilstore_id(jj) = soilstore_id(jj) - dI*sfr(is)/sfr(jj)  !Check (HCW 13/08/2014) - QUESTION: why adjust for jj and not is?

                     ! If insufficient water in first surface to move dI, instead move as much as possible
                  ELSEIF ((soilstore_id(is) + dI) < 0) THEN
                     soilstore_id(jj) = soilstore_id(jj) + soilstore_id(is)*sfr(is)/sfr(jj) !HCW 12/08/2014 switched order of these two lines
                     soilstore_id(is) = 0    !Check (HCW 13/08/2014) - QUESTION: can SM actually go to zero, or is this inconsistent with SMres?

                     ! If insufficient water in second surface to move dI, instead move as much as possible
                  ELSE
                     soilstore_id(is) = soilstore_id(is) + soilstore_id(jj)*sfr(jj)/sfr(is)
                     soilstore_id(jj) = 0
                  ENDIF

                  !If soil moisture exceeds capacity, excess goes to soil runoff (first surface)
                  IF (soilstore_id(is) > SoilStoreCap(is)) THEN
                     runoffSoil(is) = runoffSoil(is) + (soilstore_id(is) - SoilStoreCap(is))
                     soilstore_id(is) = SoilStoreCap(is)
                     !elseif (soilstore_id(is)<0) then  !HCW 13/08/2014 commented out as should never be true here anyway...
                     !   soilstore_id(is)=0             ! ... and if so, need to do more here (i.e. account for other water too)
                  ENDIF

                  !If soil moisture exceeds capacity, excess goes to soil runoff (second surface)
                  IF (soilstore_id(jj) > SoilStoreCap(jj)) THEN
                     runoffSoil(jj) = runoffSoil(jj) + (soilstore_id(jj) - SoilStoreCap(jj))
                     soilstore_id(jj) = SoilStoreCap(jj)
                     !elseif (soilstore_id(jj)<0) then  !HCW 13/08/2014 commented out (as above)
                     !         soilstore_id(jj)=0
                  ENDIF

               ENDIF  !end if second surface exists and is capable of storing water

            ENDDO  !end jj loop over second surface

            runoffSoil_per_tstep = runoffSoil_per_tstep + (runoffSoil(is)*sfr(is)/NonWaterFraction)  !Excludes water body. Moved here as otherwise code crashed when NonWaterFraction=0

         ENDIF  !end if first surface exists and is capable of storing water

         !runoffSoil_per_tstep=runoffSoil_per_tstep+(runoffSoil(is)*sfr(is)/NonWaterFraction)  !Excludes water body

      ENDDO !is loop over first surface

   END SUBROUTINE SUEWS_cal_HorizontalSoilWater
   !===================================================================================

   !===================================================================================
   SUBROUTINE SUEWS_cal_WaterUse( &
      nsh_real, & ! input:
      wu_m3, SurfaceArea, sfr, &
      IrrFracConif, IrrFracDecid, IrrFracGrass, &
      DayofWeek_id, WUProfA_24hr, WUProfM_24hr, &
      InternalWaterUse_h, HDD_id, WUDay_id, &
      WaterUseMethod, NSH, it, imin, DLS, &
      WUAreaEveTr_m2, WUAreaDecTr_m2, & ! output:
      WUAreaGrass_m2, WUAreaTotal_m2, &
      wu_EveTr, wu_DecTr, wu_Grass, int_wu, ext_wu)
      ! Conversion of water use (irrigation)
      ! Last modified:
      !  TS 30 Oct 2018  - fixed a bug in external water use
      !  TS 08 Aug 2017  - addded explicit interface
      !  LJ  6 Apr 2017  - WUchoice changed to WaterUseMethod
      !  TK 14 Mar 2017  - Corrected the variable name WUAreaEveTr_m2 -> WUAreaGrass_m2 (row 35)
      !                    Corrected conversion from m to mm /1000 -> *1000 (row 47 and 60)
      !  LJ 27 Jan 2016  - Removing Tab:s and cleaning the code
      !  HCW 12 Feb 2015 - Water use [mm] now inidcates the amount of water supplied for each surface
      !  HCW 26 Jan 2015 - Water use [mm] is the same for each surface at the moment and indicates the
      !                    amount of water supplied for each irrigated area
      !
      ! To Do:
      !        - Add functionality for water on paved surfaces (street cleaning, fountains)

      IMPLICIT NONE
      INTEGER, PARAMETER :: nsurf = 7

      REAL(KIND(1d0)), INTENT(in)::nsh_real
      REAL(KIND(1d0)), INTENT(in)::wu_m3 ! external water input (e.g., irrigation) in m^3
      REAL(KIND(1d0)), INTENT(in)::SurfaceArea !Surface area of the study area [m2]
      REAL(KIND(1d0)), INTENT(in)::sfr(nsurf)!Surface fractions [-]
      REAL(KIND(1d0)), INTENT(in)::IrrFracConif!Fraction of evergreen trees which are irrigated
      REAL(KIND(1d0)), INTENT(in)::IrrFracDecid!Fraction of deciduous trees which are irrigated
      REAL(KIND(1d0)), INTENT(in)::IrrFracGrass!Fraction of grass which is irrigated
      REAL(KIND(1d0)), INTENT(in)::InternalWaterUse_h !Internal water use [mm h-1]
      ! WUProfA_tstep(24*NSH,2),& !Automatic water use profiles at model timestep
      ! WUProfM_tstep(24*NSH,2),& !Manual water use profiles at model timestep
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::WUProfA_24hr !Automatic water use profiles at hourly scales
      REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in)::WUProfM_24hr !Manual water use profiles at hourly scales

      REAL(KIND(1d0)), DIMENSION(12), INTENT(in)::HDD_id !HDD(id-1), Heating Degree Days (see SUEWS_DailyState.f95)
      REAL(KIND(1d0)), DIMENSION(9), INTENT(in)::WUDay_id!WUDay(id-1), Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)

      INTEGER, INTENT(in):: &
         DayofWeek_id(3), & !DayofWeek(id) 1 - day of week; 2 - month; 3 - season
         WaterUseMethod, & !Use modelled (0) or observed (1) water use
         NSH, &!Number of timesteps per hour
         it, & !Hour
         imin, & !Minutes
         DLS !day lightsavings =1 + 1h) =0
      !  nsurf

      REAL(KIND(1d0)), INTENT(out):: &
         WUAreaEveTr_m2, &
         WUAreaDecTr_m2, &
         WUAreaGrass_m2, &
         WUAreaTotal_m2, &
         wu_EveTr, &
         wu_DecTr, &
         wu_Grass, &
         int_wu, &
         ext_wu

      REAL(KIND(1d0)):: &
         InternalWaterUse, &    !Internal water use for the model timestep [mm]
         WuFr = 1, &
         wu = 0!Water use for the model timestep [mm]
      INTEGER:: ih   !Hour corrected for Daylight savings
      INTEGER:: iu   !1=weekday OR 2=weekend
      INTEGER :: tstep ! timestep in second
      REAL(KIND(1d0)), PARAMETER::NAN = -999.
      REAL(KIND(1d0)):: OverUse
      REAL(KIND(1d0)):: rain_cum_daily ! accumulated daily rainfall

      REAL(KIND(1d0)):: get_Prof_SpecTime_sum

      ! NB: set OverUse as 0 as done module_constants, TS 22 Oct 2017
      ! and the logic for calculating OverUse to be determined
      OverUse = 0

      ! timestep in second
      tstep = INT(3600/NSH)

      ! accumulated daily rainfall
      rain_cum_daily = HDD_id(11)

      ! --------------------------------------------------------------------------------
      ! If water used is observed and provided in the met forcing file, units are m3
      ! Divide observed water use (in m3) by water use area to find water use (in mm)
      IF (WaterUseMethod == 1) THEN   !If water use is observed
         ! Calculate water use area [m2] for each surface type
         WUAreaEveTr_m2 = IrrFracConif*sfr(ConifSurf)*SurfaceArea
         WUAreaDecTr_m2 = IrrFracDecid*sfr(DecidSurf)*SurfaceArea
         WUAreaGrass_m2 = IrrFracGrass*sfr(GrassSurf)*SurfaceArea
         WUAreaTotal_m2 = WUAreaEveTr_m2 + WUAreaDecTr_m2 + WUAreaGrass_m2

         !Set water use [mm] for each surface type to zero initially
         wu_EveTr = 0
         wu_DecTr = 0
         wu_Grass = 0
         IF (wu_m3 == NAN .OR. wu_m3 == 0) THEN !If no water use
            ! wu_m3=0
            wu = 0
         ELSE                            !If water use
            IF (WUAreaTotal_m2 > 0) THEN
               wu = (wu_m3/WUAreaTotal_m2*1000)  !Water use in mm for the whole irrigated area
               IF (WUAreaEveTr_m2 > 0) THEN
                  wu_EveTr = wu                    !Water use for Irr EveTr in mm - these are all the same at the moment
                  wu_EveTr = wu_EveTr*IrrFracConif !Water use for EveTr in mm
               ENDIF
               IF (WUAreaDecTr_m2 > 0) THEN
                  wu_DecTr = wu                        !Water use for Irr DecTr in mm - these are all the same at the moment
                  wu_DecTr = wu_DecTr*IrrFracDecid     !Water use for DecTr in mm
               ENDIF
               IF (WUAreaGrass_m2 > 0) THEN
                  wu_Grass = wu                    !Water use for Irr Grass in mm - these are all the same at the moment
                  wu_Grass = wu_Grass*IrrFracGrass !Water use for Grass in mm
               ENDIF
               wu = (wu_m3/SurfaceArea*1000)     !Water use for the whole study area in mm
            ENDIF
         ENDIF

         ! --------------------------------------------------------------------------------
         ! If water use is modelled, calculate at timestep of model resolution [mm]
      ELSEIF (WaterUseMethod == 0) THEN   !If water use is modelled

         ! Account for Daylight saving
         ih = it - DLS
         IF (ih < 0) ih = 23

         ! Weekday or weekend profile
         iu = 1     !Set to 1=weekday
         !  IF(DayofWeek(id,1)==1.OR.DayofWeek(id,1)==7) THEN
         IF (DayofWeek_id(1) == 1 .OR. DayofWeek_id(1) == 7) THEN
            iu = 2  !Set to 2=weekend
         ENDIF

         !write(*,*) (NSH*(ih+1-1)+imin*NSH/60+1)

         ! ---- Automatic irrigation ----
         ! wu_EveTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(2)   !Automatic evergreen trees
         ! wu_DecTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(5)   !Automatic deciduous trees
         ! wu_Grass = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(8)   !Automatic grass
         wu_EveTr = get_Prof_SpecTime_sum(ih, imin, 0, WUProfA_24hr(:, iu), tstep)*WUDay_id(2)   !Automatic evergreen trees
         wu_DecTr = get_Prof_SpecTime_sum(ih, imin, 0, WUProfA_24hr(:, iu), tstep)*WUDay_id(5)   !Automatic deciduous trees
         wu_Grass = get_Prof_SpecTime_sum(ih, imin, 0, WUProfA_24hr(:, iu), tstep)*WUDay_id(8)   !Automatic grass

         ! PRINT*, ''
         ! PRINT*, 'WUDay_id(2) ',WUDay_id(2)
         ! PRINT*, 'profile ',get_Prof_SpecTime_sum(ih,imin,0,WUProfA_24hr(:,iu),tstep)
         ! PRINT*, 'manual:'
         ! PRINT*, 'wu_EveTr',wu_EveTr
         ! PRINT*, 'wu_DecTr',wu_DecTr
         ! PRINT*, 'wu_Grass',wu_Grass

         ! ---- Manual irrigation ----
         WuFr = 1 !Initialize WuFr to 1, but if raining, reduce manual fraction of water use
         ! If cumulative daily precipitation exceeds 2 mm
         IF (rain_cum_daily > 2) THEN    !.and.WUDay(id-1,3)>0) then !Commented out HCW 23/01/2015
            WuFr = 0   ! 0 -> No manual irrigation if raining
         ENDIF

         ! Add manual to automatic to find total irrigation
         ! wu_EveTr = wu_EveTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(3)) !Manual evergreen trees
         ! wu_DecTr = wu_DecTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(6)) !Manual deciduous trees
         ! wu_Grass = wu_Grass + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WUDay_id(9)) !Manual grass
         wu_EveTr = wu_EveTr + (get_Prof_SpecTime_sum(ih, imin, 0, WUProfM_24hr(:, iu), tstep)*WuFr*WUDay_id(3)) !Manual evergreen trees
         wu_DecTr = wu_DecTr + (get_Prof_SpecTime_sum(ih, imin, 0, WUProfM_24hr(:, iu), tstep)*WuFr*WUDay_id(6)) !Manual deciduous trees
         wu_Grass = wu_Grass + (get_Prof_SpecTime_sum(ih, imin, 0, WUProfM_24hr(:, iu), tstep)*WuFr*WUDay_id(9)) !Manual grass

         ! PRINT*, 'auto:'
         ! PRINT*, 'wu_EveTr',wu_EveTr
         ! PRINT*, 'wu_DecTr',wu_DecTr
         ! PRINT*, 'wu_Grass',wu_Grass
         ! Added HCW 12 Feb 2015.
         !wu_EveTr=wu_EveTr*sfr(ConifSurf)*IrrFracConif        !Water use for EveTr [mm]
         !wu_DecTr=wu_DecTr*sfr(DecidSurf)*IrrFracDecid        !Water use for DecTr [mm]
         !wu_Grass=wu_Grass*sfr(GrassSurf)*IrrFracGrass        !Water use for Grass [mm]
         wu_EveTr = wu_EveTr*IrrFracConif  !Water use for EveTr [mm]
         wu_DecTr = wu_DecTr*IrrFracDecid  !Water use for DecTr [mm]
         wu_Grass = wu_Grass*IrrFracGrass  !Water use for Grass [mm]

         ! PRINT*, 'auto:'
         ! PRINT*, 'IrrFracConif',IrrFracConif
         ! PRINT*, 'IrrFracDecid',IrrFracDecid
         ! PRINT*, 'IrrFracGrass',IrrFracGrass

         ! Total water use for the whole study area [mm]
         wu = wu_EveTr*sfr(ConifSurf) + wu_DecTr*sfr(DecidSurf) + wu_Grass*sfr(GrassSurf)

      ENDIF   !End WU_choice
      ! --------------------------------------------------------------------------------

      ! Internal water use is supplied in SUEWS_Irrigation in mm h-1
      ! Convert to mm for the model timestep
      InternalWaterUse = InternalWaterUse_h/nsh_real

      ! Remove InternalWaterUse from the total water use
      ext_wu = wu - (InternalWaterUse + OverUse)
      ! Check ext_wu cannot be negative
      IF (ext_wu < 0) THEN
         overUse = ABS(ext_wu)
         ext_wu = 0
      ELSE
         OverUse = 0
      ENDIF

      int_wu = wu - ext_wu

      ! Decrease the water use for each surface by the same proportion
      IF (ext_wu /= 0 .AND. wu /= 0) THEN
         wu_EveTr = wu_EveTr*ext_wu/wu
         wu_DecTr = wu_DecTr*ext_wu/wu
         wu_Grass = wu_Grass*ext_wu/wu
      ENDIF

   END SUBROUTINE SUEWS_cal_WaterUse
   !===================================================================================

END MODULE WaterDist_module
