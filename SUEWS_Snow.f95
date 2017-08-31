!This subroutine makes snow related calculations at the model time step. Needed for the
!available energy in LUMPS and SUEWS. Made by LJ in Dec 2012
!SUBROUTINES:
!  MeltHeat - Calculation of snow related energy processes
!  SnowCalc - Calculation of snow and soil storages
!  Evap_SUEWS_Snow - Calculation of evaporation from the snowpack
!  snowRem - Removal of snow my snow clearing
!  SnowDepletionCurve - Calculation of snow fractions
!Last modified
!  TS 31 Aug 2017 - fixed the incomplete explicit interfaces
!  LJ 24 Aug 2017 - added explicit interfaces
!  LJ 3 May 2016  - Changed so that not all surface water freezes in 5-min timestep.
!                    Re-organization of the snow routine due to this change
!                    Calculation of albedo moved from MeltHeat to SnowCalc
!  LJ 27 Jan 2016  - Tabs removed, cleaning of the code
!  HCW 08 Dec 2015 - Added check for no Paved surfaces
!  LJ 14 July 2015 - Code fixed to work with tstep.
!  HCW 06 Mar 2015 - Unused variable 'i' removed.
!  LJ Jan 2015     - Change the calculation from hourly timestep to timestep defined by nsh
!  LJ May 2013     - Calculation of the energy balance for the snowpack was modified
!                        to use qn1_ind_snow(surf)
!=======================================================================================

SUBROUTINE MeltHeat(&!input
     nsurf,&
     bldgsurf,&
     PavSurf,&
     WaterSurf,&
     nsh_real,&
     lvS_J_kg,&
     lv_J_kg,&
     tstep_real,&
     RadMeltFact,&
     TempMeltFact,&
     SnowAlb,&
     SnowDens,&
     SnowAlbMax,&
     SnowDensMin,&
     Qm_melt,&
     Qm_freezState,&
     Qm_rain,&
     deltaQi,&
     FreezMelt,&
     FreezState,&
     FreezStateVol,&
     rainOnSnow,&
     Tsurf_ind,&
     state,&
     sfr,&
     Temp_C,&
     Precip,&
     PrecipLimit,&
     PrecipLimitAlb,&
     qn1_ind_snow,&
     SnowDepth,&
     Meltwaterstore,&
     waterdens,&
     SnowPack,&      !updated state
     CumSnowfall,&
     snowFrac,&
     mwh,&
     fwh,&
     Qm,&
     QmFreez,&
     QmRain,&
     snowCalcSwitch,&
     mw_ind)       !output



  IMPLICIT NONE

  INTEGER,INTENT(in)::bldgsurf,&
       nsurf,&
       PavSurf,&
       WaterSurf

  INTEGER,DIMENSION(nsurf),INTENT(out)::snowCalcSwitch

  !These are input to the module
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::sfr,&
       Tsurf_ind,&
       state,&
       qn1_ind_snow,&
       Meltwaterstore,&
       deltaQi


  !Output, dimension nsurf
  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::Qm_melt,&
       Qm_freezState,&
       Qm_rain,&
       FreezMelt,&
       FreezState,&
       FreezStateVol,&
       rainOnSnow,&
       SnowDepth

  !Input and output as this is updated in this subroutine
  REAL(KIND(1d0)),DIMENSION(nsurf)::SnowPack,&
       snowFrac,&
       mw_ind,&
       SnowAlb,&
       SnowDens


  REAL(KIND(1d0)),INTENT(in)::&
       lvS_J_kg,&
       lv_J_kg,&
       tstep_real,&
       RadMeltFact,&
       TempMeltFact,&
       SnowAlbMax,&
       SnowDensMin,&
       Temp_C,&
       Precip,&
       PrecipLimit,&
       PrecipLimitAlb,&
       nsh_real,&
       waterdens

  REAL(KIND(1d0)),INTENT(out)::&
       mwh,&
       fwh,&
       Qm,&
       QmFreez,&
       QmRain,&
       CumSnowfall

  REAL(KIND(1d0))::AdjMeltFact,&
       cw=4190,&  !,ci=2090   !Specific heat capacity of water
       Watfreeze

  INTEGER::is

  !Initialize snow variables
  snowCalcSwitch=0
  Qm_melt=0
  Qm_freezState=0
  Qm_rain=0
  FreezMelt=0
  FreezState=0
  FreezStateVol=0
  rainOnSnow = 0


  !=========================================================================================
  DO is=1,nsurf  !Go each surface type through
     IF (sfr(is)/=0) THEN  !If surface type existing,

        IF (SnowPack(is)>0) THEN  !If snowpack existing, calculate meltwater related water flows

           SnowDepth(is) = (SnowPack(is)/1000)*waterDens/SnowDens(is) !Snow depth in m

           !Calculate meltwater related water flows with hourly degree-day method.

           !These are for snow melting
           IF (Temp_C>=0) THEN
              IF (qn1_ind_snow(is)<0) THEN
                 mw_ind(is) = TempMeltFact*Temp_C             !(mm C−1 h−1)*(C) = in mm h-1
              ELSE
                 mw_ind(is) = RadMeltFact*(qn1_ind_snow(is))  !(mm m2 W−1 h−1)*(W m-2)= mm h-1 ??
              ENDIF

           ELSE  !Freezing equations
              AdjMeltFact=1  !Relationship between the temperature melt and freezing factors
              mw_ind(is) = TempMeltFact*Temp_C*AdjMeltFact ! in mm h-1
           ENDIF

           !Previous equation give the hourly values, divide these with the timestep number
           mw_ind(is) = mw_ind(is)/nsh_real

           IF (mw_ind(is)>SnowPack(is)) mw_ind(is) = SnowPack(is)!Limited by the previous timestep snowpack

           !-----------------------------------------------------
           ! Heat consumed to snowmelt/refreezing within Tstep.
           ! Converted from mm nsh-1 to mm nsh-1 and to m s-1
           Qm_melt(is) = waterDens*((mw_ind(is)/tstep_real)/1000)*(lvS_J_kg-lv_J_kg)

           !If melt is negative this means freezing water in the snowpack
           IF (mw_ind(is)<0) THEN

              FreezMelt(is) = -mw_ind(is) !Save this to variable FreezMelt
              mw_ind(is) = 0

              !Freezing water cannot exceed meltwater store
              IF (FreezMelt(is)>Meltwaterstore(is)) FreezMelt(is) = Meltwaterstore(is)

              !Recalculate melt related energy
              Qm_melt(is) = waterDens*((-FreezMelt(is)/tstep_real)/1000)*(lvS_J_kg-lv_J_kg)
           ENDIF

           !-----------------------------------------------------
           ! If air temperature is above zero, precipitation causes advective heat to the
           ! snowpack. Eq (23) in Sun et al., 1999
           ! Calculation done at resolution of the model timestep
           IF (Temp_C>=PrecipLimit.AND.Precip>0) THEN
              Qm_rain(is) = waterDens*cw*(Temp_C-PrecipLimit)*(Precip*0.001/tstep_real)  !in W m-2
              IF (Qm_rain(is)<0) THEN !Can only be positive
                 Qm_rain(is) = 0
              ELSE
                 rainOnSnow(is) = Precip !Save information on the rain on snow event
              ENDIF
           ENDIF

        ENDIF !End if snowpack

        !=================================================================

        !Freeze surface water state if cold enough.
        IF (Tsurf_ind(is)<0.AND.state(is)>0) THEN

           snowCalcSwitch(is)=1 !If water on ground this forms ice and snow calculations are made

           !Other surfaces than water treated first
           IF (is/=WaterSurf) THEN

              !FreezState(is) = state(is)
              !Previously all state could freeze in 5-min timestep. Now we calculate how much water
              !can freeze in a timestep based on the same temperature freezing fraction.
              FreezState(is) = -TempMeltFact*Tsurf_ind(is)/nsh_real

              !The amount of freezing water cannot be greater than the surface state
              IF (FreezState(is)>state(is)) FreezState(is) = state(is)

              IF (snowpack(is)==0.OR.snowfrac(is)==0) THEN !Snowpack forms
                 FreezStateVol(is) = FreezState(is)
              ELSE                                         !There is snow already on ground
                 FreezStateVol(is) = FreezState(is)*(1-snowFrac(is))/snowFrac(is)
              ENDIF

              ! If the amount of freezing water is very small and there is state left to the ground
              ! no freezing of water will take place
              IF (FreezStateVol(is)<0.00000000001.AND.FreezState(is)<state(is)) THEN
                 FreezState(is) = 0
                 FreezStateVol(is) = 0
              ENDIF

              !Calculate the heat exchange in W m-2
              Qm_freezState(is) = -waterDens*(FreezState(is)/tstep_real/1000)*(lvS_J_kg-lv_J_kg)

              !Water surface separately
           ELSE
              !Calculate average value how much water can freeze above the water areas
              !Equation is -hA(T-T0) = rhoV(Cp+dT +Lf) in 5-min timestep
              !h=convective heat trasnfer,A, area of water,rwo water density,V volume, dT temperature difference
              !before and end of the 5-min period. dT equals zero, h=100 and when multiplied with Area, the equation
              !simplyfies to the this. LJ 14 July 2015
              Watfreeze = 100*(0-Temp_C)/(waterDens*(lvS_J_kg-lv_J_kg))
              FreezState(is) = Watfreeze
              Qm_freezState(is) = -waterDens*(Watfreeze/tstep_real/1000)*(lvS_J_kg-lv_J_kg)
           ENDIF

        ENDIF

        !======================================================================
        ! Define if any snowmelt calculations are made: snowpack existing,
        ! freezing occuring on ground or from precip
        IF (is/=WaterSurf) THEN
           IF (snowPack(is)>0.OR.(Precip>0.AND.Tsurf_ind(is)<0)) THEN
              snowCalcSwitch(is)=1
           ENDIF
        ELSE       !Water surface separately
           IF (snowPack(WaterSurf)>0.OR.FreezState(WaterSurf)>0) THEN
              snowCalcSwitch(WaterSurf)=1
           ENDIF
        ENDIF

        !Update snow density of each surface
        IF (Precip>0.AND.Tsurf_ind(is)<0.AND.SnowPack(is)>0) THEN
           SnowDens(is) = SnowDens(is)*snowPack(is)/(snowPack(is)+Precip)+SnowDensMin*Precip/(snowPack(is)+Precip)
        ENDIF

        !Weighted variables for the whole area
        mwh = mwh + mw_ind(is)*sfr(is)*snowFrac(is)        !Snowmelt
        fwh = fwh + FreezMelt(is)*sfr(is)*snowFrac(is)     !Freezing water
        Qm = Qm + Qm_melt(is)*sfr(is)*snowFrac(is)         !Energy consumed to the melt/freezing.
        QmRain = QmRain + Qm_rain(is)*sfr(is)*snowFrac(is) !Rain on snow
        QmFreez=QmFreez+deltaQi(is)*sfr(is)*snowFrac(is)+Qm_freezState(is)*sfr(is)*(1-snowFrac(is)) !Freezing water
     ENDIF

  ENDDO !End surface type

  !Update snow albedo to its maximum value if precipitation exists
  IF (Precip>0.AND.SUM(snowPack)>0.AND.Temp_C<0) THEN

     CumSnowfall=CumSnowfall + Precip

     IF (CumSnowfall>PrecipLimitAlb) THEN
        SnowAlb=SnowAlbMax
        CumSnowfall=0
     ENDIF
  ELSE
     CumSnowfall=0
  ENDIF


END SUBROUTINE MeltHeat


!===============================================================================================
!===============================================================================================
SUBROUTINE SnowCalc(&
     ity,&
     id,&  !input
     nsurf,&
     tstep,&
     imin,&
     it,&
     is,&
     snowfractionchoice,&
     nsh_real,&
     DayofWeek,&
     CRWmin,&
     CRWmax,&
     lvS_J_kg,&
     lv_j_kg,&
     avdens,&
     waterdens,&
     avRh,&
     Press_hPa,&
     precip,&
     Temp_C,&
     RAsnow,&
     psyc_hPa,&
     avcp,&
     sIce_hPa,&
     surf,&
     ConifSurf,&
     BSoilSurf,&
     BldgSurf,&
     PavSurf,&
     WaterSurf,&
     PervFraction,&
     VegFraction,&
     addimpervious,&
     soilstorecap,&
     sfr,&
     SnowDens,&
     snowdensmin,&
     Qm_Melt,&
     Qm_rain,&
     Tsurf_ind,&
     drain,&
     SurplusEvap,&
     snowPack,&
     snowFrac,&
     mw_ind,&
     rainonsnow,&
     freezmelt,&
     freezstate,&
     freezstatevol,&
     MeltWaterStore,&
     runoffPipes,&
     mwstore,&
     runoffwaterbody,&
     iceFrac,&
     addwater,&
     addwaterrunoff,&
     SnowDepth,&
     WetThresh,&
     numPM,&
     s_hPa,&
     ResistSurf,&
     sp,&
     ra,&
     rb,&
     tlv,&
     runoffSnow,&   !Output
     runoff,&
     runoffSoil,&
     chang,&
     changSnow,&
     SnowToSurf,&
     addVeg,&
     ev,&
     ev_snow,&
     state,&
     stateold,&
     swe,&
     chSnow_per_interval,&
     ev_per_tstep,&
     qe_per_tstep,&
     runoff_per_tstep,&
     surf_chang_per_tstep,&
     PipeCapacity,&
     RunoffToWater,&
     runoffAGimpervious,&
     runoffAGveg,&
     FlowChange,&
     surpluswaterbody,&
     snowprof,&
     soilmoist,&
     snowD)

  !Calculation of snow and water balance on 5 min timestep. Treats snowfree and snow covered
  !areas separately. Weighting is taken into account in the overall values.
  !Last modified:
  !  LJ in 6 May 2015 - Modified to run with timestep
  !  HCW 06 Mar 2015 - Unused variable 'i' removed.
  !  HCW 26 Jan 2015 - Added weekday/weekend option for snow clearing profiles
  !  LJ in 24 May 2013
  !========================================================================

  !USE allocateArray
  !USE defaultNotUsed
  !USE data_in
  !USE gis_data
  !USE mod_k
  !USE mod_z
  !USE moist
  !USE snowMod
  !USE sues_data
  !USE thresh
  !USE time

  IMPLICIT NONE



  INTEGER,INTENT(in)::id,&
       nsurf,&
       tstep,&
       imin,&
       it,&
       is,&
       snowfractionchoice,&
       ConifSurf,&
       BSoilSurf,&
       BldgSurf,&
       PavSurf,&
       WaterSurf
  INTEGER,INTENT(in) :: ity!Evaporation calculated according to Rutter (1) or Shuttleworth (2)

  REAL(KIND(1d0)),INTENT(in)::CRWmin,&
       CRWmax,&
       nsh_real,&
       lvS_J_kg,&
       lv_j_kg,&
       avdens,&
       waterdens,&
       avRh,&
       Press_hPa,&
       Temp_C,&
       RAsnow,&
       psyc_hPa,&
       avcp,&
       sIce_hPa,&
       PervFraction,&
       vegfraction,&
       addimpervious,&
       numPM,&
       s_hPa,&
       ResistSurf,&
       sp,&
       ra,&
       rb,&
       tlv,&
       snowdensmin,&
       precip,&
       PipeCapacity,&
       RunoffToWater,&
       runoffAGimpervious,&
       runoffAGveg,&
       surpluswaterbody

  REAL(KIND(1d0)),DIMENSION(366,2),INTENT(in)::DayofWeek
  REAL(KIND(1d0)),DIMENSION(6,nsurf),INTENT(in):: surf


  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(out)::runoffSnow,& !Initialize for runoff caused by snowmelting
       runoff,&
       runoffSoil,&
       chang,&
       changSnow,&
       SnowToSurf,&
       state,&
       snowD,&
       ev_snow,&
       soilmoist

  REAL(KIND(1d0)), DIMENSION(0:23,2),INTENT(out):: snowProf

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)::drain,&
       WetThresh,&
       stateold,&
       mw_ind,&
       soilstorecap,&
       rainonsnow,&
       freezmelt,&
       freezstate,&
       freezstatevol,&
       Qm_Melt,&
       Qm_rain,&
       Tsurf_ind,&
       sfr


  REAL(KIND(1d0)),INTENT(out)::swe,ev

  !Updated status: input and output
  REAL(KIND(1d0)),DIMENSION(nsurf)::snowPack,&
       snowFrac,&
       MeltWaterStore,&
       SnowDepth,&
       iceFrac,&
       addwater,&
       addwaterrunoff,&
       SnowDens

  REAL(KIND(1d0))::chSnow_per_interval,&
       ev_per_tstep,&
       qe_per_tstep,&
       runoff_per_tstep,&
       surf_chang_per_tstep,&
       runoffPipes,&
       mwstore,&
       runoffwaterbody,&
       FlowChange

  REAL(KIND(1d0)),DIMENSION(2):: SurplusEvap


  REAL(KIND(1d0))::&
       qe,&
       rss,&
       addVeg,&
       Evap_SUEWS_Snow,&
       MeltExcess,&      !Excess melt water that needs to leave snowpack
       snowTotInit,&
       EvPart,&
       runoffTest,&
       snowFracFresh1,&   !Snow fraction for newly formed snowpack
       snowFracFresh2,&   !Snow fraction for newly formed snowpack from state only
       snowFracOld,&
       WaterHoldCapFrac,&
       FWC                !Water holding capacity of snow in mm

  REAL(KIND(1d0)):: SnowDepletionCurve

  INTEGER:: iu                        !1=weekday OR 2=weekend
  REAL(KIND(1d0)),PARAMETER :: IPThreshold_mmhr = 10   !Threshold for intense precipitation [mm hr-1]
  !========================================================================
  !Initialize variables for the calculation of water storages and evaporation

  ! Use weekday or weekend snow clearing profile
  iu=1     !Set to 1=weekday
  IF(DayofWeek(id,1)==1.OR.DayofWeek(id,1)==7) iu=2  !Set to 2=weekend

  !write(*,*) is
  runoffSnow(is)=0 !Initialize for runoff caused by snowmelting
  runoff(is)=0
  runoffSoil(is)=0
  chang(is)=0
  changSnow(is)=0
  runoffTest=0
  SnowToSurf(is)=0
  EvPart=0
  ev=0
  snowFracFresh1=0
  snowFracFresh2=0
  snowFracOld=0

  !Initial snowpack + meltwater in it
  snowTotInit=snowPack(is)+MeltWaterStore(is)

  !Calculate water holding capacity (Jin et al. 1999)
  IF (SnowDens(is)>=200) THEN
     WaterHoldCapFrac=CRWmin
  ELSE
     WaterHoldCapFrac=CRWmin+(CRWmax-CRWmin)*(200-SnowDens(is))/200
  ENDIF

  !======================================================================
  ! Calculate evaporation from snowpack and snow free surfaces (in mm)
  ! IF (snowFrac(is)<1) CALL Evap_SUEWS !ev and qe for snow free surface out

  IF (snowFrac(is)<1) CALL Evap_SUEWS(&

                                ! input:
       ity,&!Evaporation calculated according to Rutter (1) or Shuttleworth (2)
       state(is),& ! wetness status
       WetThresh(is),&!When State > WetThresh, rs=0 limit in SUEWS_evap [mm] (specified in input files)
       surf(6,is),& ! = surf(is,6), current storage capacity [mm]
       numPM,&!numerator of P-M eqn
       s_hPa,&!Vapour pressure versus temperature slope in hPa
       psyc_hPa,&!Psychometric constant in hPa
       ResistSurf,&!Surface resistance
       sp,&!Term in calculation of E
       ra,&!Aerodynamic resistance
       rb,&!Boundary layer resistance
       tlv,&!Latent heat of vaporization per timestep [J kg-1 s-1], (tlv=lv_J_kg/tstep_real)

                                ! output:
       rss,&
       ev,&
       qe) ! latent heat flux [W m-2]

  IF (snowFrac(is)>0) THEN
     ev_snow(is) = Evap_SUEWS_Snow(Qm_Melt(is),Qm_rain(is),lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RAsnow,&
          psyc_hPa,tstep,avcp,sIce_hPa)
  ENDIF

  !If not enough water for evaporation in impervious surfaces,
  !evaporation is taken from pervious surfaces
  IF (is>2) THEN
     IF  (PervFraction/=0) THEN
        EvPart=(SurplusEvap(PavSurf)*sfr(PavSurf)+SurplusEvap(BldgSurf)*sfr(BldgSurf))/PervFraction
     ENDIF
  ENDIF


  !============================================================================
  !Water surface is treated separately
  IF (is==WaterSurf.AND.sfr(WaterSurf)>0) GO TO 606

  !The calculations are divided into 2 main parts
  ! 1) Surface is fully covered with snow at the beginning of the time step
  ! 2) Surface is not fully covered with snow but rather part is snow free OR
  !    surface not orginally covered with snow, but the snow forms at the current timestep

  !1)------------------------------------------------------------------
  !  ------------------------------------------------------------------
  IF (SnowPack(is)>0.AND.snowFrac(is)==1) THEN

     ev_snow(is)=ev_snow(is)+EvPart !Evaporation surplus

     !(Snowfall per interval+freezing of melt water and surface state) - (meltwater+evaporation from snowpack)
     changSnow(is)=(Precip+freezMelt(is))-(mw_ind(is)+ev_snow(is)) !Calculate change in snowpack (in mm)

     !If rain on snow event, add this water to meltwaterstore
     IF (rainOnSnow(is)>0) THEN
        changSnow(is)=changSnow(is)-Precip
        MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)
     ENDIF

     SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack

     !---------If snowpack exists after the state calculations
     IF (snowPack(is)>0) THEN

        !Add melted water to meltstore and freeze water according to freezMelt(is)
        MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is) - freezMelt(is)

        !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
        FWC = WaterHoldCapFrac*snowPack(is)

        !If FWC is exceeded, excess meltwater (MeltExcess) will leave from the snowpack
        IF (MeltWaterStore(is)>=FWC) THEN
           MeltExcess = 0                      !Initialize the excess meltwater
           MeltExcess = MeltWaterStore(is)-FWC !Calculate the exceess water
           MeltWaterStore(is) = FWC            !Update the meltwaterstore to the maximum it can hold
           runoffSnow(is) = runoffSnow(is) + MeltExcess
        ENDIF

        !At the end of the hour calculate possible snow removal
        IF (SnowProf(it,iu)==1.AND.is<3.AND.(imin==(nsh_real-1)/nsh_real*60))  CALL snowRem

        !----------If snowPack is negative, it melts at this timestep
     ELSEIF (SnowPack(is)<0) THEN

        !If freezing meltwater inside this timestep, remove it from the MeltWaterStore
        MeltWaterStore(is)=MeltWaterStore(is)-freezMelt(is)+mw_ind(is)+SnowPack(is)
        SnowPack(is)=0.0   !Set the snow pack and snow
        snowFracOld=1
        snowFrac(is)=0
        snowDens(is)=0

        IF (MeltWaterStore(is)<0) THEN !Not enough water in the meltwater store,
           ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased.??
           IF (ev_snow(is)<0) ev_snow(is)=0
           changSnow(is)=changSnow(is)+MeltWaterStore(is)
           MeltWaterStore(is)=0
        ELSE
           chang(is)=MeltWaterStore(is)  !Meltwater goes to surface state as no snow exists anymore
           state(is)=state(is)+chang(is)
           MeltWaterStore(is)=0
        ENDIF
     ENDIF !snowpack negative or positive


     !2)------Surface not fully covered with snow-------------------------------------------
     !  ------------------------------------------------------------------------------------
  ELSEIF (snowFrac(is)<1) THEN

     !Snow calculations: snowpack can either exist or form at the current timestep
     IF (SnowPack(is)>0) THEN
        ev_snow(is)=ev_snow(is)+EvPart !Evaporation surplus


        !----Snowpack water balance for the whole surface area. In reality snow depth = snowPack/snowFrac(is)
        !(Snowfall per interval+freezing of melt water and surface state) - (meltwater+evaporation from snowpack)
        changSnow(is)=(Precip+freezMelt(is)+freezStateVol(is))-(mw_ind(is)+ev_snow(is)) !Calculate change in snowpack (in mm)

        !If rain on snow event, add this water to meltwaterstore
        IF (rainOnSnow(is)>0) THEN
           changSnow(is)=changSnow(is)-Precip
           MeltWaterStore(is) = MeltWaterStore(is)+rainOnSnow(is)
        ENDIF
        SnowPack(is)=SnowPack(is)+changSnow(is)


        !The fraction of snow will update when:
        !a) Surface state is dry but precipitation occurs =1
        !b) There is both precipitation and all surface state freezes =1
        !c) No precipitation but all state freezes at a single timestep =2
        !d) Part of the surface freezes
        IF (Precip>0.AND.FreezState(is)==state(is)) THEN !both a) and b)
           snowFracFresh1=1
        ELSEIF (Precip==0.AND.FreezState(is)>0.AND.FreezState(is)==state(is)) THEN
           snowFracFresh1=1

           !snowFracFresh1=SnowDepletionCurve(is,snowPack(is),snowD(is))
           !if (snowFracFresh1<0.001) snowFracFresh1=0.001
        ELSEIF (FreezState(is)>0.AND.FreezState(is)<state(is)) THEN !This if not all water freezes
           snowFracFresh1=0.95 !Now this fraction set to something close to one. Should be improved in the future at some point
           !if (is==1)then
           ! write(*,*) id,it,imin,snowfrac(is),FreezState(is),state(is)
           ! pause
           !endif
        ENDIF

        !Snowpack can also form at the current timestep (2). If this forms purely from snowfall or/and all water at surface freezes,
        !the whole surface will be covered with snow. If there is water on ground this snowfall can immediately melt
        !and in this case the snow fraction is not necessarily 1 but its information is saved to snowFracFresh that
        !is taken into account in snow fraction after calculation of state.
     ELSEIF (snowpack(is)==0.AND.Tsurf_ind(is)<0) THEN

        !The fraction of snow will get a value of 1 (ie full snow cover):
        !Surface state is dry but precipitation occurs, no precipitation but all state freezes at a single timestep,
        !There is both precipitation and all surface state freezes
        IF ((Precip>0.AND.state(is)==0).OR.(Precip==0.AND.FreezState(is)==state(is)).OR.&
             (Precip>0.AND.FreezState(is)==state(is))) THEN

           !ev=ev+EvPart
           changSnow(is)=Precip+FreezStateVol(is)
           SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack

           snowFracFresh1=1
           iceFrac(is)=FreezState(is)/(FreezState(is)+Precip)
           SnowDens(is)=SnowDensMin
        ENDIF

        IF (FreezState(is)>0.AND.FreezState(is)<state(is)) THEN

           changSnow(is)=Precip+freezStateVol(is)
           SnowPack(is)=SnowPack(is)+changSnow(is)  !Update snowpack
           snowFracFresh2=0.95 !Now this fraction set to something close to one. Should be improved in the future at some point

           !snowFracFresh2=SnowDepletionCurve(is,snowPack(is),snowD(is))
           !if (snowFracFresh2<0.001) snowFracFresh2=0.001
           iceFrac(is)=1
           SnowDens(is)=SnowDensMin
           !write(*,*) 2,is,id,it,imin,snowfrac(is),FreezState(is),state(is),state(is)+Precip
           !pause

        ENDIF
     ENDIF

     !---------If snowpack exists after the state calculations
     IF (snowPack(is)>0) THEN

        !Add melted water to meltstore and freeze water according to freezMelt(is)
        MeltWaterStore(is) = MeltWaterStore(is) + mw_ind(is) - freezMelt(is)

        !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
        FWC = WaterHoldCapFrac*snowPack(is)

        !If FWC is exceeded, excess meltwater (MeltExcess) will leave from the snowpack
        IF (MeltWaterStore(is)>=FWC) THEN
           MeltExcess = 0                      !Initialize the excess meltwater
           MeltExcess = MeltWaterStore(is)-FWC !Calculate the exceess water
           MeltWaterStore(is) = FWC            !Update the meltwaterstore to the maximum it can hold

           !If the fraction of snow is greater than 0.8 or if the surface is is buildings,
           !the excess water will directly go to runoff. Otherwise it will flow to the
           !snow free area via SnowToSurf(is)
           IF ((snowFrac(is)>0.9.AND.is/=BldgSurf).OR.(is==BldgSurf)) THEN
              runoffSnow(is) = runoffSnow(is) + MeltExcess
           ELSE
              SnowToSurf(is) = SnowToSurf(is) + MeltExcess*snowFrac(is)/(1-snowFrac(is))
           ENDIF
        ENDIF

        !At the end of the hour calculate possible snow removal
        IF (SnowProf(it,iu)==1.AND.is<3.AND.(imin==(nsh_real-1)/nsh_real*60))  CALL snowRem

        !----------If snowPack is negative, it melts at this timestep
     ELSEIF (SnowPack(is)<0) THEN

        !If freezing meltwater inside this timestep, remove it from the MeltWaterStore
        MeltWaterStore(is)=MeltWaterStore(is)-freezMelt(is)+mw_ind(is)+SnowPack(is)

        SnowPack(is)=0.0   !Set the snow pack and snow
        snowFracFresh1=0
        snowFracFresh2=0
        snowDens(is)=0

        IF (MeltWaterStore(is)<0) THEN !Not enough water in the meltwater store,
           ev_snow(is)=ev_snow(is)+MeltWaterStore(is) !evaporation from snow is decreased.??
           IF (ev_snow(is)<0) ev_snow(is)=0
           changSnow(is)=changSnow(is)+MeltWaterStore(is)
           MeltWaterStore(is)=0
        ELSE
           SnowToSurf(is)=SnowToSurf(is)+MeltWaterStore(is)*snowFrac(is)/(1-snowFrac(is))
           MeltWaterStore(is)=0
        ENDIF
     ENDIF !snowpack negative or positive


     !--------
     !Next the snow free surface (3). Calculations only done if snowfraction is smaller than 1
     IF ((is==PavSurf.OR.is==BldgSurf).AND.snowFrac(is)<1) THEN  !Impervious surfaces (paved, buildings)

        !Surface store update. If precipitation is greater than the threshold, the exceeding water
        !goes directly to runoff
        IF (precip>IPThreshold_mmhr/nsh_real) THEN
           !runoff = runoff + (precipitation+water from the snow surface+water from other surfaces-the thereshold limit)
           runoff(is)=runoff(is)+(Precip+SnowToSurf(is)+AddWater(is)-IPThreshold_mmhr/nsh_real)
           chang(is)=IPThreshold_mmhr/nsh_real-(drain(is)+ev+freezState(is))
        ELSE
           !Add precip and water from other surfaces and remove drainage, evap and freezing of state
           chang(is)=Precip+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is))
        ENDIF

        state(is)=state(is)+chang(is) !Change in state (for whole surface area areasfr(is))

        !Add water from impervious grids
        ! Check sfr/=0 added HCW 08 Dec 2015
        IF (is==PavSurf.AND.sfr(PavSurf)/=0) state(is)=state(is)+(addImpervious)/sfr(PavSurf)

        runoff(is)=runoff(is)+drain(is)*AddWaterRunoff(is) !Drainage (not flowing to other surfaces) goes to runoff

        IF(state(is)<0.0) THEN  !Surface state cannot be negative
           SurplusEvap(is)=ABS(state(is)) !take evaporation from other surfaces in mm
           ev = ev-SurplusEvap(is)
           state(is)=0.0
        ENDIF

     ELSEIF(is>=3.AND.snowFrac(is)<1) THEN ! Pervious surfaces (conif, decid, grass unirr, grass irr)

        ev=ev+EvPart

        !Change in water stores
        IF (Precip+addVeg*(sfr(is)/VegFraction)>(IPThreshold_mmhr/nsh_real)) THEN !if 5min precipitation is larger than 10 mm
           runoff(is)=runoff(is)+(Precip+addVeg*(sfr(is)/VegFraction)+SnowToSurf(is)+AddWater(is)-(IPThreshold_mmhr/nsh_real))
           chang(is)=(IPThreshold_mmhr/nsh_real)-(drain(is)+ev+freezState(is))
        ELSE
           chang(is)=Precip+addVeg*(sfr(is)/VegFraction)+SnowToSurf(is)+AddWater(is)-(drain(is)+ev+freezState(is))
        ENDIF

        state(is)=state(is)+chang(is)

        !Add water in soil store only if ground is not frozen
        IF (Temp_C>0) THEN
           soilmoist(is)=soilmoist(is)+Drain(is)*AddWaterRunoff(is)*(1-snowFrac(is))
        ELSE
           runoff(is)=runoff(is)+Drain(is)*AddWaterRunoff(is)
        ENDIF

        !If state of the surface is negative, remove water from soilstore
        IF(state(is)<0.0) THEN

           IF ((soilmoist(is)+state(is))>=0.AND.Temp_C>0) THEN !If water in soilstore, water is removed

              soilmoist(is)=soilmoist(is)+state(is)*(1-snowFrac(is))
              state(is)=0.0

           ELSE !If not water in the soilstore evaporation does not occur
              chang(is)=chang(is)+state(is)
              ev=ev+state(is)
              state(is)=0.0
           ENDIF
        ENDIF !state is negative

        !If soilstorage is full at this point, excess will go to surface runoff
        IF (soilmoist(is)>soilstoreCap(is)) THEN
           runoffTest=runoffTest+(soilmoist(is)-soilstoreCap(is))
           soilmoist(is)=soilstoreCap(is)
        ELSEIF (soilmoist(is)<0) THEN
           soilmoist(is)=0
        ENDIF

     ENDIF !Surface type

  ENDIF !Surface fraction

  !-------------------------------------------------------------------------------------------------------------------

  !Calculate change in snowpack and state for the respective surface areas
  !Here the case where not all surface state freezes is handled
  IF (snowFracFresh2>0) THEN
     surf_chang_per_tstep=surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is))&
          -Precip*sfr(is)*(1-snowFracFresh2)
     chSnow_per_interval=chSnow_per_interval+((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)*(1-snowFrac(is))&
          -Precip*sfr(is)*snowFracFresh2
  ELSE
     surf_chang_per_tstep=surf_chang_per_tstep+(state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is))
     chSnow_per_interval=chSnow_per_interval+((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)*MAX(snowFrac(is),snowfracOld)
  ENDIF

  !Add evaporation to total
  IF (is==BldgSurf.OR.is==PavSurf) THEN
     ev_per_tstep=ev_per_tstep+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)
     qe_per_tstep=qe_per_tstep+ev_snow(is)*lvS_J_kg*sfr(is)*snowFrac(is)&
          +ev*lv_J_kg*sfr(is)*(1-snowFrac(is))
  ELSE
     ev_per_tstep=ev_per_tstep+ev*sfr(is)*(1-snowFrac(is))+ev_snow(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)
     qe_per_tstep=qe_per_tstep+ev_snow(is)*lvS_J_kg*sfr(is)*MAX(snowFrac(is),snowfracOld)+ev*lv_J_kg*sfr(is)*(1-snowFrac(is))
  ENDIF

  !========RUNOFF=======================

  !Add runoff to pipes
  runoffPipes=runoffPipes+runoffSnow(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)+runoff(is)*sfr(is)*(1-snowFrac(is))&
       +runoffTest*sfr(is)
  CALL updateFlood(&
                                ! input:
       nsurf,is,PavSurf,BldgSurf,WaterSurf,ConifSurf,BSoilSurf,&
       sfr,PipeCapacity,RunoffToWater,&
                                ! inout:
       runoffAGimpervious,surplusWaterBody,runoffAGveg,runoffPipes&
       )

  runoff_per_tstep=runoff_per_tstep+runoffSnow(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)+runoff(is)*sfr(is)*(1-snowFrac(is))&
       +runoffTest*sfr(is)

  !===Update snow depth, weighted SWE, and Mwstore
  IF (SnowDens(is)/=0) THEN
     SnowDepth(is) = SnowPack(is)*waterDens/SnowDens(is)
  ENDIF

  ! Calculate overall snow water equivalent
  swe = swe + SnowPack(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)
  MwStore = MwStore + MeltWaterStore(is)*sfr(is)*MAX(snowFrac(is),snowfracOld)

  !if (id==6.and.it==13.and.imin==20) then!
  !if (id==85.and.it==3.and.imin==10) then!
  ! if (id==92.and.it==21.and.imin==35) then!
  !  write(*,*)  ((snowPack(is)+MeltWaterstore(is))-snowTotInit)*sfr(is)*(1-snowFrac(is)),&
  !              runoff(is)*sfr(is)*(1-snowFrac(is)),&
  !              ev*sfr(is)*(1-snowFrac(is)),&
  !              (state(is)-stateOld(is))*sfr(is)*(1-snowFrac(is)),Precip*sfr(is)
  !  write(*,*)  changSnow(is),runoff(is),ev,chang(is),runoffTest,FreezState(is) !changSnow(is)-freezMelt(is)
  !  write(*,*)  is,Precip,runoff_per_tstep,ev_per_tstep,surf_chang_per_tstep,chSnow_per_interval
  !  write(*,*)  is,Precip-runoff_per_tstep-ev_per_tstep,surf_chang_per_tstep+chSnow_per_interval
  !  write(*,*)  is,snowFrac(is),sfr(is),sfr(is)*ev_snow(is)
  !  pause
  ! endif

  !Only now update the new snow fractions both in the case that snow existing already on ground
  !and snow forms at the current timestep
  IF (snowFracFresh1>0) snowFrac(is)=snowFracFresh1
  IF (snowFracFresh2>0) snowFrac(is)=snowFracFresh2

  !Calculate new snow fraction here.
  !Tässä ongelmana että snow fraction muuttuu vain kun on sulamisvettä ja on vika tunti.
  !Tämä ei juuri koskaan toteudu johtuen lämpötilan vuorokausisyklistä
  !Kokeile tässä ajaa kahdella tavalla 1) ei tarvita Mw:tä
  !                                    2) päivitys voi tapahtua millon vain
  !if (SnowFractionChoice==2.and.imin==(nsh_real-1)/nsh_real*60) then
  IF (SnowFractionChoice==2) THEN
     IF (snowPack(is)>0.AND.mw_ind(is)>0) THEN
        snowFrac(is) = SnowDepletionCurve(is,snowPack(is),snowD(is))
        IF (snowFrac(is)<0.001) snowFrac(is)=0.001  !The snow fraction minimum is 1% of the surface
     ELSEIF (snowPack(is)==0) THEN
        snowFrac(is)=0
     ENDIF
  ENDIF

  RETURN

  !==========================================================================
  !WATERBODY is treated separately as state always below ice if ice existing
  !Calculate change in snowpack
606 changSnow(WaterSurf)=(Precip+freezMelt(WaterSurf)+freezState(WaterSurf))-&
       (mw_ind(WaterSurf)+ev_snow(WaterSurf))

  SnowPack(WaterSurf)=SnowPack(WaterSurf)+changSnow(WaterSurf) !Update snowpack
  state(WaterSurf)=state(WaterSurf)+FlowChange-freezState(WaterSurf)  !Update state below ice

  !If snowpack exists
  IF (snowPack(WaterSurf)>0) THEN

     !Add melted water to meltstore and freeze water according to freezMelt(is)
     MeltWaterStore(WaterSurf)=MeltWaterStore(WaterSurf)+mw_ind(WaterSurf)-freezMelt(WaterSurf)

     !Calculate water holding capacity (FWC: Valeo and Ho, 2004) of the snowpack
     FWC = WaterHoldCapFrac*SnowPack(WaterSurf)

     !If FWC is exceeded, add meltwater to state
     IF (MeltWaterStore(WaterSurf)>=FWC.AND.Temp_C>=0) THEN
        state(WaterSurf)=state(WaterSurf)+(MeltWaterStore(WaterSurf)-FWC)
        MeltWaterStore(WaterSurf) = FWC
     ENDIF

     !If snowPack is negative, it melts at this timestep
  ELSEIF (SnowPack(is)<0) THEN

     !Add water to the meltwater store
     !If freezing meltwater inside this hour, remove it from the MeltWaterStore
     MeltWaterStore(WaterSurf) = MeltWaterStore(WaterSurf)-freezMelt(WaterSurf) &
          + mw_ind(WaterSurf)

     state(WaterSurf)=state(WaterSurf)+MeltWaterStore(WaterSurf)+SnowPack(WaterSurf) !Add meltwater to state
     SnowPack(WaterSurf)=0
     IF (state(WaterSurf)<0) ev_snow(WaterSurf)=ev_snow(WaterSurf)+state(WaterSurf)

  ENDIF !snowpack negative or positive

  !Check water state separately
  IF (state(WaterSurf)>Surf(5,WaterSurf)) THEN
     runoff(WaterSurf)=runoff(WaterSurf)+(state(WaterSurf)-Surf(5,WaterSurf))
     state(WaterSurf)=Surf(5,WaterSurf)
     runoffWaterBody=runoffWaterBody+runoff(WaterSurf)*sfr(WaterSurf)
  ELSE
     state(WaterSurf)=state(WaterSurf)+surplusWaterBody

     IF (state(WaterSurf)>Surf(5,WaterSurf)) THEN
        runoffWaterBody=runoffWaterBody+(state(WaterSurf)-Surf(5,WaterSurf))*sfr(WaterSurf)
        state(WaterSurf)=Surf(5,WaterSurf)
     ENDIF
  ENDIF


  !Change state of snow and surface
  chSnow_per_interval=chSnow_per_interval+((snowPack(WaterSurf)+MeltWaterstore(WaterSurf))-snowTotInit)*sfr(WaterSurf)
  !ch_per_interval=ch_per_interval+(state(WaterSurf)-stateOld(WaterSurf))*sfr(WaterSurf)
  surf_chang_per_tstep=surf_chang_per_tstep+(state(WaterSurf)-stateOld(WaterSurf))*sfr(WaterSurf)

  !Evaporation
  ev_per_tstep=ev_per_tstep+ev*sfr(WaterSurf)+ev_snow(WaterSurf)*sfr(WaterSurf)
  qe_per_tstep=qe_per_tstep+ev_snow(WaterSurf)*lvS_J_kg*sfr(WaterSurf)+ev*lv_J_kg*sfr(WaterSurf)
  runoff_per_tstep=runoff_per_tstep+(runoff(is)*sfr(is)) !The total runoff from the area

  IF (snowPack(WaterSurf)>0) THEN     !Fraction only 1 or 0
     snowFrac(WaterSurf)=1
  ELSE
     snowFrac(WaterSurf)=0
  ENDIF

END SUBROUTINE SnowCalc


!==========================================================================
!==========================================================================
!Calculates evaporation from snow surface (ev_snow).

FUNCTION Evap_SUEWS_Snow(Qm,QP,lvS_J_kg,avdens,avRh,Press_hPa,Temp_C,RAsnow,psyc_hPa,&
     tstep,avcp,sIce_hPa) RESULT(ev_snow)

  IMPLICIT NONE

  !INPUT
  REAL (KIND(1d0))::Qm,QP,&        !melt heat, advect. heat
       lvS_J_kg,avdens,avRh,&   !latent heat of sublimation, air density,relative humidity,
       Press_hPa,Temp_C,&       !air pressure, air temperature
       RAsnow,psyc_hPa,&        !aerodyn res snow, psychometric constant, type of evaporation calculation
       avcp,sIce_hPa            !spec. heat, satured curve on snow

  !OTHER VARIABLES
  REAL (KIND(1d0))::e_snow,&     !PM equation obe line
       sae_snow,&   !s * (Available energy)
       qe_snow,&    !Latent heat flux
       ev_snow,&    !Evaporation
       vdrcIce,&    !Vapour pressure deficit
       esIce_hPa,&  !Saturation vapor pressure over ice
       EaIce_hPa,&  !Vapour pressure
       tlv_sub,&    !Latent heat for sublimation
       tstep_real   !timestep as real

  REAL (KIND(1d0)):: sat_vap_pressIce !Function

  INTEGER:: tstep,from=1
  !-----------------------------------------------------

  tstep_real = REAL(tstep,KIND(1d0))

  sae_snow=sIce_hPa*(Qp-Qm)   !Calculate the driving parameter in calculation of evaporation. Järvi et al. (2015)

  esIce_hPa= sat_vap_pressIce(Temp_C,Press_hPa,from) !Saturation vapor pressure over ice
  EaIce_hPa=avRh/100*esIce_hPa                       !Vapour pressure of water
  vdrcIce=(esIce_hPa-eaIce_hpa)*avdens*avcp          !Vapour pressure deficit
  tlv_sub=lvS_J_kg/tstep_real                        !Latent heat for sublimation

  e_snow=sae_snow+vdrcIce/RAsnow                     !PM equation
  qe_snow=e_snow/(sIce_hPa+psyc_hPa)                 !Latent heat (W/m^2)
  ev_snow=qe_snow/tlv_sub                            !Evaporation (in mm)

  RETURN

END FUNCTION Evap_SUEWS_Snow

!==========================================================================
!==========================================================================
! Calculates mechanical removal of snow from roofs ans roads
SUBROUTINE snowRem

  USE allocateArray
  USE snowMod
  USE sues_data
  USE time

  IMPLICIT NONE

  !write(*,*) is, SnowPack(is),SnowLimPaved,SnowLimBuild

  IF (is==PavSurf) THEN
     IF (SnowPack(PavSurf)>SnowLimPaved) THEN
        SnowRemoval(PavSurf) = (SnowPack(PavSurf)-SnowLimPaved)*sfr(PavSurf)*snowfrac(PavSurf)
        SnowPack(PavSurf)=SnowLimPaved
        !snowPack(PavSurf)=snowPack(PavSurf)/snowFrac(PavSurf)
     ENDIF
  ENDIF
  IF (is==BldgSurf)THEN
     IF (SnowPack(BldgSurf)>SnowLimBuild) THEN
        SnowRemoval(2) = (SnowPack(BldgSurf)-SnowLimBuild)*sfr(BldgSurf)*snowfrac(BldgSurf)
        SnowPack(BldgSurf)=SnowLimBuild
        !snowPack(BldgSurf)=snowPack(BldgSurf)/snowFrac(BldgSurf)
     ENDIF
  ENDIF
  !write(*,*) is, SnowPack(is),SnowLimPaved,SnowLimBuild
  !pause
END SUBROUTINE snowRem

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION SnowDepletionCurve(is,swe,sweD) RESULT(asc)
  !This function calculates surface coverage of snow according to the
  !depletion curves in Valeo and Ho (2004).
  !INPUT: is   Surface type number
  !       swe  Snow water content
  !       sweD Limit for

  USE allocateArray

  IMPLICIT  NONE

  INTEGER::is
  REAL (KIND(1d0))::asc,sweD,swe


  !Impervious surface
  IF (is==PavSurf) THEN

     IF (swe<=sweD) THEN      !Snow water equivalent below threshold
        asc=((swe/sweD))**2
     ELSE
        asc=1
     ENDIF

     !Bldgs surface
  ELSEIF (is==BldgSurf) THEN

     IF (swe<=sweD) THEN
        IF ((swe/sweD)<0.9) THEN
           asc=(swe/sweD)*0.5
        ELSE
           asc=(swe/sweD)**8
        ENDIF
     ELSE
        asc=1
     ENDIF
  ELSEIF (is==WaterSurf) THEN
     IF (swe>0) asc=1

     !Vegetion surfaces
  ELSE
     IF (swe<=sweD) THEN

        asc=1-((1/3.1416)*ACOS(2*(swe/sweD)-1))**1.7
     ELSE
        asc=1
     ENDIF

  ENDIF

  !asc=real(int(10000.*asc))/10000  !4 decimal precision

  RETURN
END FUNCTION SnowDepletionCurve
