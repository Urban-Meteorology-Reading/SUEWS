SUBROUTINE soilstore(&
     is,& ! input: ! surface type
     sfr,&! surface fractions
     PipeCapacity,&!Capacity of pipes to transfer water
     RunoffToWater,&!Fraction of surface runoff going to water body
     pin,&!Rain per time interval
     wu_EveTr,&!Water use for evergreen trees/shrubs [mm]
     wu_DecTr,&!Water use for deciduous trees/shrubs [mm]
     wu_Grass,&!Water use for grass [mm]
     AddWater,&!Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
     addImpervious,&!Water from impervious surfaces of other grids [mm] for whole surface area
     nsh_real,&!nsh cast as a real for use in calculations
     stateOld,&!Wetness status of each surface type from previous timestep [mm]
     AddWaterRunoff,&!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
     PervFraction,&! sum of surface cover fractions for impervious surfaces
     addVeg,&!Water from vegetated surfaces of other grids [mm] for whole surface area
     soilstoreCap,&!Capacity of soil store for each surface [mm]
     addWaterBody,&!Water from water surface of other grids [mm] for whole surface area
     FlowChange,&!Difference between the input and output flow in the water body
     StateLimit,&!Limit for state of each surface type [mm] (specified in input files)
     runoffAGimpervious,&!  inout:!Above ground runoff from impervious surface [mm] for whole surface area
     surplusWaterBody,&!Extra runoff that goes to water body [mm] as specified by RunoffToWater
     runoffAGveg,&!Above ground runoff from vegetated surfaces [mm] for whole surface area
     runoffPipes,&!Runoff in pipes [mm] for whole surface area
     ev,&!Evaporation
     soilmoist,&!Soil moisture of each surface type [mm]
     SurplusEvap,&!Surplus for evaporation in 5 min timestep
     runoffWaterBody,&!Above ground runoff from water surface [mm] for whole surface area
     runoff_per_interval,&! Total water transported to each grid for grid-to-grid connectivity
     p_mm,&!output: !Inputs to surface water balance
     chang,&!Change in state [mm]
     runoff,&!Runoff from each surface type [mm]
     drain,&!Drainage of each surface type [mm]
     state&!Wetness status of each surface type [mm]
     )
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
  !   - Water from other grids now added to p_mm (instead of state for impervious surfaces)
  !   - Removed division of runoff by nsh, as whole model now runs at the same timestep
  !   - Adjusted transfer of ev between surfaces to conserve mass (not depth)
  !   - Volumes used for water transport between grids to account for SurfaceArea changing between grids
  !   - Added threshold check for state(WaterSurf) - was going negative
  ! Last modified HCW 09 Feb 2015
  !   - Removed StorCap input because it is provided by module allocateArray
  !   - Tidied and commented code
  ! Modified by LJ in November 2012:
  !   - P>10 was not taken into account for impervious surfaces - Was fixed.
  !   - Above impervious surfaces possibility of the state to exceed max capacity was limited
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

  ! USE allocateArray
  ! USE data_in
  ! USE defaultNotUsed
  ! USE gis_data
  ! USE sues_data
  ! USE thresh
  ! USE time

  IMPLICIT NONE

  !Stores flood water when surface state exceeds storage capacity [mm]
  !real(kind(1d0)),dimension(nsurf):: SurfaceFlood
  INTEGER, PARAMETER:: nsurf=7                !Total number of surfaces
  ! INTEGER, PARAMETER:: NVegSurf=3             !Number of surfaces that are vegetated
  ! INTEGER, PARAMETER:: nsurfIncSnow=nsurf+1   !Number of surfaces + snow

  INTEGER:: PavSurf   = 1,&   !When all surfaces considered together (1-7)
       BldgSurf  = 2,&
       ConifSurf = 3,&
       DecidSurf = 4,&
       GrassSurf = 5,&   !New surface classes: Grass = 5th/7 surfaces
       BSoilSurf = 6,&   !New surface classes: Bare soil = 6th/7 surfaces
       WaterSurf = 7
      !  ExcessSurf= 8,&   !Runoff or subsurface soil in WGWaterDist
      !  NSurfDoNotReceiveDrainage=0,&   !Number of surfaces that do not receive drainage water (green roof)
      !  ivConif = 1,&     !When only vegetated surfaces considered (1-3)
      !  ivDecid = 2,&
      !  ivGrass = 3

  INTEGER,INTENT(in)::is ! surface type


  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr(nsurf)! surface fractions
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::AddWater(nsurf)!Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::stateOld(nsurf)!Wetness status of each surface type from previous timestep [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::AddWaterRunoff(nsurf)!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::soilstoreCap(nsurf)!Capacity of soil store for each surface [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::StateLimit(nsurf)!Limit for state of each surface type [mm] (specified in input files)

  REAL(KIND(1d0)),INTENT(in)::PipeCapacity!Capacity of pipes to transfer water
  REAL(KIND(1d0)),INTENT(in)::RunoffToWater!Fraction of surface runoff going to water body
  REAL(KIND(1d0)),INTENT(in)::pin!Rain per time interval
  REAL(KIND(1d0)),INTENT(in)::wu_EveTr!Water use for evergreen trees/shrubs [mm]
  REAL(KIND(1d0)),INTENT(in)::wu_DecTr!Water use for deciduous trees/shrubs [mm]
  REAL(KIND(1d0)),INTENT(in)::wu_Grass!Water use for grass [mm]
  REAL(KIND(1d0)),INTENT(in)::addImpervious!Water from impervious surfaces of other grids [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(in)::nsh_real!nsh cast as a real for use in calculations
  REAL(KIND(1d0)),INTENT(in)::PervFraction! sum of surface cover fractions for impervious surfaces
  REAL(KIND(1d0)),INTENT(in)::addVeg!Water from vegetated surfaces of other grids [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(in)::addWaterBody!Water from water surface of other grids [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(in)::FlowChange!Difference between the input and output flow in the water body


  REAL(KIND(1d0)),INTENT(inout)::runoffAGimpervious!Above ground runoff from impervious surface [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(inout)::surplusWaterBody!Extra runoff that goes to water body [mm] as specified by RunoffToWater
  REAL(KIND(1d0)),INTENT(inout)::runoffAGveg!Above ground runoff from vegetated surfaces [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(inout)::runoffPipes!Runoff in pipes [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(inout)::ev!Evaporation
  REAL(KIND(1d0)),INTENT(inout)::runoffWaterBody!Above ground runoff from water surface [mm] for whole surface area
  REAL(KIND(1d0)),INTENT(inout)::runoff_per_interval! Total water transported to each grid for grid-to-grid connectivity

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::soilmoist  !Soil moisture of each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(2),INTENT(inout)    ::SurplusEvap!Surplus for evaporation in 5 min timestep

  REAL(KIND(1d0)),INTENT(out)::p_mm!Inputs to surface water balance

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::chang !Change in state [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoff!Runoff from each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::drain !Drainage of each surface type [mm]
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::state !Wetness status of each surface type [mm]

  !Extra evaporation [mm] from impervious surfaces which cannot happen due to lack of water
  REAL(KIND(1d0)):: EvPart
  REAL(KIND(1d0)),PARAMETER:: NotUsed=-55.5,&
       IPThreshold_mmhr=10

  !Initialise extra evaporation to zero
  EvPart=0

  !SurfaceFlood(is) = 0 !!This probably needs to be carried over between timesteps, but reset for now

  !==================================================================
  ! Combine water inputs to the current surface
  ! Add external water use for each surface type
  IF(is==ConifSurf) THEN
     p_mm=pin+wu_EveTr
  ELSEIF(is==DecidSurf) THEN
     p_mm=pin+wu_DecTr
  ELSEIF(is==GrassSurf) THEN
     p_mm=pin+wu_Grass
  ELSE
     p_mm=pin
  ENDIF

  ! Add water from other surfaces within the same grid (RS2S) ----
  ! AddWater is the water supplied to the current surface from other surfaces
  !  i.e. drain*WaterDist (see SUEWS_ReDistributeWater)
  p_mm=p_mm+AddWater(is)

  !==== Impervious surfaces (Paved, Buildings) ======================
  IF(is==PavSurf.OR.is==BldgSurf) THEN

     ! Add water from neighbouring grids (RG2G)
     ! Add to PavSurf only, as water cannot flow onto buildings
     IF (is==PavSurf) THEN
        IF(sfr(PavSurf)/=0) THEN   ! If loop added HCW 08 Dec 2015
           p_mm=p_mm+addImpervious/sfr(PavSurf)
        ENDIF
     ENDIF

     ! Calculate change in surface state (inputs - outputs)
     chang(is)=p_mm-(drain(is)+ev)

     ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
     ! and does not affect state
     IF(p_mm>IPThreshold_mmhr/nsh_real) THEN
        runoff(is)=runoff(is)+(p_mm-IPThreshold_mmhr/nsh_real)
        chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev)
     ENDIF

     ! Calculate updated state using chang
     state(is)=state(is)+chang(is)

     ! Check state is within physical limits between zero (dry) and max. storage capacity
     IF(state(is)<0.0) THEN   ! Cannot have a negative surface state
        ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
        ! Allow evaporation only until surface is dry (state(is)=0); additional evaporation -> evaporation surplus
        SurplusEvap(is)=ABS(state(is))   !Surplus evaporation is that which tries to evaporate non-existent water
        ev = ev-SurplusEvap(is)          !Limit evaporation according to water availability
        state(is)=0.0                    !Now surface is dry
        ! elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
        !    !! If state exceeds the storage capacity, then the excess goes to surface flooding
        !    !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
        !    runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
        !    state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
     ENDIF

     ! Recalculate change in surface state from difference with previous timestep
     chang(is) = state(is)-stateOld(is)

     ! Runoff -------------------------------------------------------
     ! For impervious surfaces, some of drain(is) becomes runoff
     runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is)   !Drainage (that is not flowing to other surfaces) goes to runoff

     !So, up to this point, runoff(is) can have contributions if
     ! p_mm > ipthreshold (water input too fast)
     ! state > surf(6,is) (net water exceeds storage capacity)
     ! WaterDist specifies some fraction of drain(is) -> runoff

     !==== Pervious surfaces (Conif, Decid, Grass, BSoil, Water) =======
  ELSEIF(is>=3) THEN

     ! Transfer evaporation surplus from impervious surfaces to pervious surfaces
     IF(PervFraction/=0) THEN   !If pervious surfaces exist
        EvPart=(SurplusEvap(PavSurf)*sfr(PavSurf)+SurplusEvap(BldgSurf)*sfr(BldgSurf))/PervFraction
     ELSE         !If no pervious surface, SurplusEvap cannot be transferred and this evap cannot
        EvPart=0  !happen (will increase QHinstead)
     ENDIF

     ! Add surplus evaporation to ev for pervious surfaces
     ev=ev+EvPart

     !==== For Conif, Decid, Grass, BSoil surfaces ==================
     IF (is/=WaterSurf) THEN

        ! ---- Add water from neighbouring grids (RG2G) ----
        ! Add to Grass and BSoil only, as water cannot flow onto trees
        IF (is==GrassSurf.OR.is==BSoilSurf) THEN
           IF ((sfr(GrassSurf)+sfr(BSoilSurf))/=0) THEN
              p_mm=p_mm+addVeg/(sfr(GrassSurf)+sfr(BSoilSurf))
           ENDIF
        ENDIF

        ! Calculate change in surface state (inputs - outputs)
        chang(is)=p_mm-(drain(is)+ev)

        ! If p_mm is too large, excess goes to runoff (i.e. the rate of water supply is too fast)
        !  and does not affect state
        IF (p_mm>IPThreshold_mmhr/nsh_real) THEN
           runoff(is)=runoff(is)+(p_mm-IPThreshold_mmhr/nsh_real)
           chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev)
        ENDIF

        ! Calculate updated state using chang
        state(is)=state(is)+chang(is)

        ! Check state is within physical limits between zero (dry) and max. storage capacity
        IF(state(is)<0.0) THEN   ! Cannot have a negative surface state
           ! If there is not sufficient water on the surface, then remove water from soilstore
           ! Allow evaporation until soilmoist is depleted and surface is dry
           IF((soilmoist(is)+state(is))>=0) THEN
              soilmoist(is)=soilmoist(is)+state(is)
              state(is)=0.0
              ! If there is not sufficient water on the surface or soilstore, then don't allow this evaporation to happen
           ELSE
              ev=ev-ABS(state(is))   !Limit evaporation according to water availability
              state(is)=0.0          !Now surface is dry
           ENDIF
           !! What about if there is some water in soilstore, but not enough to provide all the water for evaporation??
           !! Is this saying water can be evaporated from the soilstore as easily as from the surface??
           !elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
           !   !! If state exceeds the storage capacity, then the excess goes to surface flooding
           !   !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
           !   runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
           !   state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
        ENDIF

        ! Recalculate change in surface state from difference with previous timestep
        chang(is) = state(is)-stateOld(is)

        !Where should this go? Used to be before previous part!!
        ! Soilmoist -------------------------------------------------
        ! For pervious surfaces (not water), some of drain(is) goes to soil storage
        ! Drainage (that is not flowing to other surfaces) goes to soil storages
        soilmoist(is)=soilmoist(is)+drain(is)*AddWaterRunoff(is)

        ! If soilstore is full, the excess will go to runoff
        IF(soilmoist(is)>soilstoreCap(is)) THEN              !! Should this also go to flooding of some sort?
           runoff(is)=runoff(is)+(soilmoist(is)-soilstoreCap(is))
           soilmoist(is)=soilstoreCap(is)
        ELSEIF (soilmoist(is)<0) THEN   !!But where does this lack of water go? !!Can this really happen here??
           CALL ErrorHint(62,'SUEWS_store: soilmoist(is) < 0 ',soilmoist(is),NotUsed,is)
           ! Code this properly - soilmoist(is) < 0 shouldn't happen given the above loops
           !soilmoist(is)=0   !Groundwater / deeper soil should kick in
        ENDIF

        !==== Water surface ========================================
     ELSEIF (is==WaterSurf) THEN

        IF(sfr(WaterSurf)/=0)THEN

           ! ---- Add water from neighbouring grids (RG2G) ----
           p_mm=p_mm+addWaterBody/sfr(WaterSurf)

           ! Calculate change in surface state (inputs - outputs)
           ! No drainage for water surface
           ! FlowChange is the difference in input and output flows [mm hr-1]   !!Should this really be a constant??
           chang(is)=p_mm+FlowChange/nsh_real-(ev)

           ! Calculate updated state using chang
           state(is)=state(is)+chang(is)

           ! Check state is within physical limits between zero (dry) and max. storage capacity
           IF(state(is)<0.0) THEN   ! Cannot have a negative surface state
              ! If there is not sufficient water on the surface, then don't allow this evaporation to happen
              ev=ev-ABS(state(is))   !Limit evaporation according to water availability
              state(is)=0.0          !Now surface is dry
              !elseif (state(is)>surf(6,is)) then   !!This should perhaps be StateLimit(is)
              !   !! If state exceeds the storage capacity, then the excess goes to surface flooding
              !   !SurfaceFlood(is)=SurfaceFlood(is)+(state(is)-surf(6,is))   !!Need to deal with this properly
              !   runoff(is)=runoff(is)+(state(is)-surf(6,is))   !!needs to go to flooding
              !   state(is)=surf(6,is)              !Now surface state is at max (storage) capacity
           ENDIF

           ! Recalculate change in surface state from difference with previous timestep
           chang(is) = state(is)-stateOld(is)

           ! If state exceeds limit, then excess goes to runoff (currently applies to water surf only)
           IF (state(WaterSurf)>StateLimit(WaterSurf)) THEN
              runoff(WaterSurf)=runoff(WaterSurf)+(state(WaterSurf)-StateLimit(WaterSurf))
              state(WaterSurf)=StateLimit(WaterSurf)
              runoffWaterBody=runoffWaterBody+runoff(WaterSurf)*sfr(WaterSurf)
           ELSE
              state(WaterSurf)=state(WaterSurf)+surplusWaterBody
              IF (state(WaterSurf)>StateLimit(WaterSurf)) THEN
                 runoffWaterBody=runoffWaterBody+(state(WaterSurf)-StateLimit(WaterSurf))*sfr(WaterSurf)
                 state(WaterSurf)=StateLimit(WaterSurf)
              ENDIF
           ENDIF

           ! Recalculate change in surface state from difference with previous timestep
           chang(is) = state(is)-stateOld(is)
        ENDIF

     ENDIF   !end of WaterSurf

  ENDIF   !end of different surfaces

  !==================================================================
  !==== RUNOFF ======================================================

  ! Need to consider areas here - SurfaceArea may vary between grids too
  ! - also implement where water for next surface is calculated (RunoffFromGrid subroutine)
  ! Calculations of the piperunoff exceedensances moved to separate subroutine so that from snow same
  ! calculations can be made. LJ in May 2015

  IF(is<WaterSurf) THEN   !Not for water body

     ! Add runoff to pipes
     runoffPipes=runoffPipes+(runoff(is)*sfr(is))
     !  CALL updateFlood
     CALL updateFlood(&
                                ! input:
          nsurf,is,PavSurf,BldgSurf,WaterSurf,ConifSurf,BSoilSurf,&
          sfr,PipeCapacity,RunoffToWater,&
                                ! inout:
          runoffAGimpervious,surplusWaterBody,runoffAGveg,runoffPipes&
          )
  ENDIF

  runoff_per_interval=runoff_per_interval+(runoff(is)*sfr(is)) !The total runoff from the area !!Check (HCW)

END SUBROUTINE soilstore
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 SUBROUTINE updateFlood(&

! input:
   nsurf,is,PavSurf,BldgSurf,WaterSurf,ConifSurf,BSoilSurf,&
   sfr,PipeCapacity,RunoffToWater,&
! inout:
runoffAGimpervious,surplusWaterBody,runoffAGveg,runoffPipes&
   )

  ! USE allocateArray
  ! USE sues_data

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nsurf,is,PavSurf,BldgSurf,WaterSurf,ConifSurf,BSoilSurf
  REAL(KIND(1d0)), INTENT(in) :: sfr(nsurf),PipeCapacity,RunoffToWater
  REAL(KIND(1d0)), INTENT(inout) :: runoffAGimpervious,surplusWaterBody,runoffAGveg,runoffPipes

  ! If pipe capacity is full, surface runoff occurs
  ! N.B. this will happen each loop (replicates pipes filling up)
  IF(runoffPipes>PipeCapacity) THEN

    !------Paved and building surface
    IF(is==PavSurf.OR.is==BldgSurf) THEN
        IF(sfr(WaterSurf)>0.0000001) THEN
           ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
           ! RunoffToWater is specified in SUEWS_SiteSelect.txt
           runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)*(1-RunoffToWater)
           surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
        ELSE
           ! Otherwise, all flood water must go to runoff
           runoffAGimpervious=runoffAGimpervious+(runoffPipes-PipeCapacity)
        ENDIF
                   !------other surfaces
     ELSEIF(is>=ConifSurf.AND.is<=BSoilSurf) THEN
        IF(sfr(WaterSurf)>0.0000001) THEN
          ! If there is some water present, the water surface will take some of the flood water (fraction RunoffToWater)
          runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)*(1-RunoffToWater)
          surplusWaterBody=surplusWaterBody+(runoffPipes-PipeCapacity)*RunoffToWater
        ELSE
          ! Otherwise, all flood water must go to runoff
          runoffAGveg=runoffAGveg+(runoffPipes-PipeCapacity)
        ENDIF
      ENDIF

      runoffPipes=PipeCapacity   !Pipes are at their max capacity

  ENDIF   !If runoff exceed pipe capacity

 END SUBROUTINE updateFlood
