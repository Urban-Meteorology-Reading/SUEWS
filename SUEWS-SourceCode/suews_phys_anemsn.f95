!===================================================================================
!Simple Anthropogenic Heat and Carbon Dioxide Parameterization routines
!This subroutine is still under development and in the equations concerning CO2 fluxes
!there are bugs and missing comments.
!Last modified
! TS 01 Jul 2018 - replaced annual HDD array with a simplied daily array HDD_id_use
! MH 29 Jun 2017 -  Finalised the code to calculate the anthropogenic emissions of heat and CO2
! HCW 21 Apr 2017 - renamed from SUEWS_SAHP.f95. Now includes CO2 fluxes as well as QF.
!                   Tidied code. Combined three subroutines into 1 to avoid repetition.
! HCW 24 Feb 2017 - Added new anthropogenic heat flux calculation (AnthropEmissionsMethod=3)
! HCW 25 Aug 2016 - Outputs base QF (part without temp. dependence)
! LJ 27 Jan 2016  - Removal of Tabs
! HCW 20 Jan 2015 - v2015 applies a profile at each model timestep
!                   these have been interpolated from the hourly profile input data (SUEWS_Profiles)
!                   vthey are now normalised (sum to 1) in InitializeSurfaceCharacteristics
!                   N.B. previous versions were not applying weekday/weekend profiles correctly
!
! EmissionsMethod = 0 - Use values in met forcing file, or default QF
! EmissionsMethod = 1 - Method according to Loridan et al. (2011) : SAHP
! EmissionsMethod = 2 - Method according to Jarvi et al. (2011)   : SAHP_2
!
!===================================================================================

SUBROUTINE AnthropogenicEmissions( &
   EmissionsMethod, &
   id, it, imin, DLS, nsh, DayofWeek_id, &
   EF_umolCO2perJ, FcEF_v_kgkm, EnEF_v_Jkm, TrafficUnits, &
   FrFossilFuel_Heat, FrFossilFuel_NonHeat, &
   MinQFMetab, MaxQFMetab, &
   PopDensDaytime, PopDensNighttime, &
   Temp_C, HDD_id, Qf_A, Qf_B, Qf_C, &
   AH_MIN, AH_SLOPE_Heating, AH_SLOPE_Cooling, &
   T_CRITIC_Heating, T_CRITIC_Cooling, &
   TrafficRate, &
   QF0_BEU, QF_SAHP, &
   Fc_anthro, Fc_metab, Fc_traff, Fc_build, &
   AHProf_24hr, HumActivity_24hr, TraffProf_24hr, PopProf_24hr)
   ! Simple anthropogenic heat parameterisation and co2 calculation
   ! Calculates QF_SAHP and Fc_anthro

   IMPLICIT NONE

   INTEGER, INTENT(in):: EmissionsMethod
   INTEGER, INTENT(in):: &
      id, & !Hour
      it, & !Hour
      imin, & !Minutes
      DLS, & !day lightsavings =1 +1h =0
      nsh !Number of timesteps per hour

   INTEGER, DIMENSION(3), INTENT(in)::DayofWeek_id   !1 - day of week; 2 - month; 3 - season

   REAL(KIND(1d0)), DIMENSION(12), INTENT(in):: HDD_id !Heating Degree Days (see SUEWS_DailyState.f95)

   REAL(KIND(1d0)), DIMENSION(2), INTENT(in):: &
      Qf_A, Qf_B, Qf_C, &    !Qf coefficients
      AH_MIN, &            !Minimum anthropogenic heat flux (AnthropHeatMethod = 1)
      AH_SLOPE_Heating, &  !Slope of the antrhropogenic heat flux calculation (AnthropHeatMethod = 1)
      AH_SLOPE_Cooling, &
      T_CRITIC_Heating, &  !Critical temperature
      T_CRITIC_Cooling, &  !Critical cooling temperature
      TrafficRate, &       !Traffic rate
      QF0_BEU

   REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in):: AHProf_24hr
   REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in):: HumActivity_24hr
   REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in):: TraffProf_24hr
   REAL(KIND(1d0)), DIMENSION(0:23, 2), INTENT(in):: PopProf_24hr

   REAL(KIND(1D0)), INTENT(in):: &
      EF_umolCO2perJ, &
      FcEF_v_kgkm, &       !CO2 Emission factor
      EnEF_v_Jkm, &        !Energy Emission factor
      TrafficUnits, &      !Traffic Units choice
      FrFossilFuel_Heat, & !Fraction of Fossil Fuel heat
      FrFossilFuel_NonHeat, &!Fraction of Fossil Fuel non heat
      MinQFMetab, &        !Minimum QF Metabolism
      MaxQFMetab, &        !Maximum QF Metabolism
      
      PopDensDaytime, &    !Daytime population density [ha-1] (i.e. workers)
      PopDensNighttime, &  !Nighttime population density [ha-1] (i.e. residents)
      Temp_C           !Air temperature

   REAL(KIND(1D0)), INTENT(out):: &
      QF_SAHP, &
      Fc_anthro, &
      Fc_metab, Fc_traff, Fc_build

   INTEGER:: &
      iu, &               !1=weekday OR 2=weekend
      ih

   REAL(KIND(1D0)):: &
      get_Prof_SpecTime_inst, &
      get_Prof_SpecTime_mean, & !external function to get profile value at sepcified timestamp
      HDD_daily, & !daily HDD
      CDD_daily, & !daily HDD
      Tair_avg_daily, & !daily mean air temperature
      DP_x_RhoPop, DP_x_RhoPop_traff, &
      MinFcMetab, MaxFcMetab, &
      QF_build, QF_metab, QF_traff, &
      QF_SAHP_base, &     !Anthropogenic heat flux calculated by SAHP (temp independent part)
      QF_SAHP_heat, &     !Anthropogenic heat flux calculated by SAHP (heating part only)
      QF_SAHP_ac, &       !AC contribution
      NumCapita, &         !Number of people in the study area per hectare [ha-1]
      PopDorNorT, &       ! Population
      ActDorNorT, &       ! Human activity
      TraffDorNorT, &     ! Traffic
      AHDorNorT          ! Anthropogenic heat
      

   ! transfer HDD values to local explict variables
   HDD_daily = HDD_id(7)
   CDD_daily = HDD_id(8)
   ! NB: temporarily use 5-day running mean to test the performance
   Tair_avg_daily = HDD_id(10)

   ! Tair_avg_daily= HDD_id_use(3) ! this is daily

   ! use average of both daytime and nighttime population values 
   ! TS 27 Dec 2018: moved from `translate` here to simplify the interface 
   NumCapita = (PopDensDaytime + PopDensNighttime)/2

   !-----------------------------------------------------------------------
   ! Account for Daylight saving
   ih = it - DLS
   IF (ih < 0) ih = 23

   ! Set weekday/weekend counter
   iu = 1   !Set to 1=weekday
   IF (DayofWeek_id(1) == 1 .OR. DayofWeek_id(1) == 7) iu = 2  !Set to 2=weekend

   ! Calculate energy emissions and CO2 from human metabolism -------------
   ! Pop dens (cap ha-1 -> cap m-2) x activity level (W cap-1)
   ! PopDorNorT   = PopProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)      ! 1=night, 2=day, 1-2=transition
   ! ActDorNorT   = HumActivity_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)  ! 1=night, 2=day, 1-2=transition
   ! TraffDorNorT = TraffProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)    ! normalise so the AVERAGE of the multipliers is equal to 1
   ! AHDorNorT    = AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)       ! normalise so the AVERAGE of the multipliers is equal to 1
   ! PRINT*, ''
   ! PRINT*, 'AHDorNorT old:',AHDorNorT
   PopDorNorT = get_Prof_SpecTime_inst(ih, imin, 0, PopProf_24hr(:, iu))
   ActDorNorT = get_Prof_SpecTime_mean(ih, imin, 0, HumActivity_24hr(:, iu))
   TraffDorNorT = get_Prof_SpecTime_mean(ih, imin, 0, TraffProf_24hr(:, iu))
   AHDorNorT = get_Prof_SpecTime_mean(ih, imin, 0, AHProf_24hr(:, iu))
   ! PRINT*, 'AHDorNorT new:',AHDorNorT
   ! Diurnal profile times population density [cap ha-1]
   DP_x_RhoPop = AHDorNorT*NumCapita

   MinFcMetab = 8/5.*MinQFMetab
   MaxFcMetab = 8/5.*MaxQFMetab ! Sailor & Lu (2004)

   !QF_metab = (PopDensNighttime*(2-PopDorNorT) + PopDensDaytime*(PopDorNorT-1))/10000 * &
   !              (MinQFMetab*(2-ActDorNorT) + MaxQFMetab*(ActDorNorT-1)) !W m-2
   !Fc_metab = (PopDensNighttime*(2-PopDorNorT) + PopDensDaytime*(PopDorNorT-1))/10000 * &
   !           (MinFcMetab*(2-ActDorNorT) + MinFcMetab*(ActDorNorT-1)) !umol m-2 s-1

   !QF_metab = (PopDensNighttime*MinQFMetab*((2-ActDorNorT)+(2-PopDorNorT))/2 + &
   !           PopDensDaytime*MaxQFMetab*((ActDorNorT-1)+(PopDorNorT-1))/2)/10000 !W m-2

   !Fc_metab = (PopDensNighttime*MinFcMetab*((2-ActDorNorT)+(2-PopDorNorT))/2 + &
   !           PopDensDaytime*MaxFcMetab*((ActDorNorT-1)+(PopDorNorT-1))/2)/10000 !umol m-2 s-1

   !----------------------------------------------------------------------------------------------------
   !Diffferent methods to calculate the anthopogenic heat and CO2 emissions.
   !1-3: CO2 emission is calculated from QF which can be calculated with three methods
   !41-43: CO2 emission is calculated with local information. QF methods are used for housing and human metabolism

   IF (EmissionsMethod == 1 .OR. EmissionsMethod == 4 .OR. EmissionsMethod == 11 .OR. EmissionsMethod == 14 .OR. &
       EmissionsMethod == 21 .OR. EmissionsMethod == 24 .OR. EmissionsMethod == 31 .OR. EmissionsMethod == 34) THEN   ! (formerly SAHP_1 subroutine)
      ! Loridan et al. (2011) JAMC Eq 13: linear relation with air temperature
      ! Weekday/weekend differences due to profile only
      ! Now scales with population density

      IF (Temp_C < T_CRITIC_Heating(iu)) THEN
         QF_SAHP = (AH_MIN(iu) + AH_SLOPE_Heating(iu)*(T_CRITIC_Heating(iu) - Temp_C))*AHDorNorT
      ELSE
         QF_SAHP = AH_MIN(iu)*AHDorNorT
      ENDIF

      ! Need to be checked later, not recommended to use
      ! QF_SAHP_base = AH_MIN(iu) * DP_x_RhoPop     ! Temperature-independent contribution
      QF_SAHP_base = AH_MIN(iu)*AHDorNorT         ! Temperature-independent contribution
      QF_SAHP_heat = QF_SAHP - QF_SAHP_base         ! Heating contribution
      QF_SAHP_ac = 0                              ! No AC contribution with this method

   ELSEIF (EmissionsMethod == 2 .OR. EmissionsMethod == 5 .OR. EmissionsMethod == 12 .OR. EmissionsMethod == 15 .OR. &
           EmissionsMethod == 22 .OR. EmissionsMethod == 25 .OR. EmissionsMethod == 32 .OR. EmissionsMethod == 35) THEN   ! (formerly SAHP_2 subroutine)
      ! Jarvi et al. (2011) JH Eq 3 using HDD and CDD
      ! Weekday/weekend differences due to profile and coefficients QF_a,b,c
      ! Scales with population density
      QF_SAHP = (Qf_a(iu) + Qf_b(iu)*CDD_daily + Qf_c(iu)*HDD_daily)*DP_x_RhoPop  !This contains QF from all three sources: buildings, metabolism and traffic!
      QF_SAHP_base = (Qf_a(iu))*DP_x_RhoPop                ! Temperature-independent contribution from buildings, traffic and human metabolism
      QF_SAHP_heat = (Qf_c(iu)*HDD_daily)*DP_x_RhoPop    ! Heating contribution
      QF_SAHP_ac = (Qf_b(iu)*CDD_daily)*DP_x_RhoPop    ! Cooling (AC) contribution

   ELSEIF (EmissionsMethod == 3 .OR. EmissionsMethod == 6 .OR. EmissionsMethod == 13 .OR. EmissionsMethod == 16 .OR. &
           EmissionsMethod == 23 .OR. EmissionsMethod == 26 .OR. EmissionsMethod == 33 .OR. EmissionsMethod == 36) THEN
      ! Updated Loridan et al. (2011) method using daily (not instantaneous) air temperature (HDD_id_use(3))
      ! Linear relation with air temperature
      ! Weekday/weekend differences due to profile only
      ! Scales with population density

      ! Need to be checked later, not recommended to use
      ! QF_SAHP_base = AH_MIN(iu) * DP_x_RhoPop       ! Temperature-independent contribution
      QF_SAHP_base = AH_MIN(iu)*AHDorNorT           ! Temperature-independent contribution

      IF (Tair_avg_daily < T_CRITIC_Heating(iu)) THEN     ! Heating
         QF_SAHP = (AH_MIN(iu) + AH_SLOPE_Heating(iu)*(T_CRITIC_Heating(iu) - Tair_avg_daily))*AHDorNorT
         QF_SAHP_heat = QF_SAHP - QF_SAHP_base        ! Heating contribution
         QF_SAHP_ac = 0

      ELSEIF (Tair_avg_daily > T_CRITIC_Cooling(iu)) THEN ! Air-conditioning
         QF_SAHP = (AH_MIN(iu) + AH_SLOPE_Cooling(iu)*(Tair_avg_daily - T_CRITIC_Cooling(iu)))*AHDorNorT
         QF_SAHP_heat = 0
         QF_SAHP_ac = QF_SAHP - QF_SAHP_base          ! AC contribution

      ELSE
         QF_SAHP = AH_MIN(iu)*AHDorNorT
      ENDIF

      ! ELSEIF(EmissionsMethod==0) THEN
      !    QF_SAHP   = 0
      !    Fc_anthro = 0
      !    Fc_metab  = 0
      !    Fc_traff  = 0
      !    Fc_build  = 0
   ENDIF

   RETURN
ENDSUBROUTINE AnthropogenicEmissions
!========================================================================================
