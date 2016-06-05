!-------------------------------------------------------------------------
SUBROUTINE InputHeaderCheck(FileName)
  ! Checks columns in input files match the columns expected by model code
  ! Model code columns are defined here
  ! Latest update:
  !   TS 02 Mar 2016  - AnOHM related variables added
  !   LJ 27 Jan 2016  - Removal of tabs
  !   LJ 07 July 2015 - snow albedo removed
  !   HCW 12 Nov 2014
  !-------------------------------------------------------------------------

  USE allocateArray
  USE ColNamesInputFiles
  USE defaultNotUsed

  IMPLICIT NONE

  CHARACTER (len=50):: FileName

  ! ========== Define expected column names here ==========
  ! =======================================================

  ! ========== SUEWS_NonVeg.txt =============
  HeaderNonVeg_Reqd(ci_Code)         = "Code"
  HeaderNonVeg_Reqd(ci_AlbMin)       = "AlbedoMin"
  HeaderNonVeg_Reqd(ci_AlbMax)       = "AlbedoMax"
  HeaderNonVeg_Reqd(ci_Emis)         = "Emissivity"
  HeaderNonVeg_Reqd(ci_StorMin)      = "StorageMin"
  HeaderNonVeg_Reqd(ci_StorMax)      = "StorageMax"
  HeaderNonVeg_Reqd(ci_WetThresh)    = "WetThreshold"
  HeaderNonVeg_Reqd(ci_StateLimit)   = "StateLimit"
  HeaderNonVeg_Reqd(ci_DrEq)         = "DrainageEq"
  HeaderNonVeg_Reqd(ci_DrCoef1)      = "DrainageCoef1"
  HeaderNonVeg_Reqd(ci_DrCoef2)      = "DrainageCoef2"
  HeaderNonVeg_Reqd(ci_SoilTCode)    = "SoilTypeCode"
  HeaderNonVeg_Reqd(ci_SnowLimPat)   = "SnowLimPatch"
  HeaderNonVeg_Reqd(ci_SnowLimRem)   = "SnowLimRemove"
  HeaderNonVeg_Reqd(ci_OHMCode_SWet) = "OHMCode_SummerWet"
  HeaderNonVeg_Reqd(ci_OHMCode_SDry) = "OHMCode_SummerDry"
  HeaderNonVeg_Reqd(ci_OHMCode_WWet) = "OHMCode_WinterWet"
  HeaderNonVeg_Reqd(ci_OHMCode_WDry) = "OHMCode_WinterDry"
  HeaderNonVeg_Reqd(ci_cpAnOHM)           = "cp"		! AnOHM TS
  HeaderNonVeg_Reqd(ci_kkAnOHM)           = "kk"		! AnOHM TS
  HeaderNonVeg_Reqd(ci_chAnOHM)           = "ch"		! AnOHM TS

  ! ========== SUEWS_Veg.txt ===============
  HeaderVeg_Reqd(cp_Code)         = "Code"
  HeaderVeg_Reqd(cp_AlbMin)       = "AlbedoMin"
  HeaderVeg_Reqd(cp_AlbMax)       = "AlbedoMax"
  HeaderVeg_Reqd(cp_Emis)         = "Emissivity"
  HeaderVeg_Reqd(cp_StorMin)      = "StorageMin"
  HeaderVeg_Reqd(cp_StorMax)      = "StorageMax"
  HeaderVeg_Reqd(cp_WetThresh)    = "WetThreshold"
  HeaderVeg_Reqd(cp_StateLimit)   = "StateLimit"
  HeaderVeg_Reqd(cp_DrEq)         = "DrainageEq"
  HeaderVeg_Reqd(cp_DrCoef1)      = "DrainageCoef1"
  HeaderVeg_Reqd(cp_DrCoef2)      = "DrainageCoef2"
  HeaderVeg_Reqd(cp_SoilTCode)    = "SoilTypeCode"
  HeaderVeg_Reqd(cp_SnowLimPat)   = "SnowLimPatch"
  HeaderVeg_Reqd(cp_BaseT)        = "BaseT"
  HeaderVeg_Reqd(cp_BaseTe)       = "BaseTe"
  HeaderVeg_Reqd(cp_GDDFull)      = "GDDFull"
  HeaderVeg_Reqd(cp_SDDFull)      = "SDDFull"
  HeaderVeg_Reqd(cp_LAIMin)       = "LAIMin"
  HeaderVeg_Reqd(cp_LAIMax)       = "LAIMax"
  HeaderVeg_Reqd(cp_GsMax)        = "MaxConductance"
  HeaderVeg_Reqd(cp_LAIEq)        = "LAIEq"
  HeaderVeg_Reqd(cp_LeafGP1)      = "LeafGrowthPower1"
  HeaderVeg_Reqd(cp_LeafGP2)      = "LeafGrowthPower2"
  HeaderVeg_Reqd(cp_LeafOP1)      = "LeafOffPower1"
  HeaderVeg_Reqd(cp_LeafOP2)      = "LeafOffPower2"
  HeaderVeg_Reqd(cp_OHMCode_SWet) = "OHMCode_SummerWet"
  HeaderVeg_Reqd(cp_OHMCode_SDry) = "OHMCode_SummerDry"
  HeaderVeg_Reqd(cp_OHMCode_WWet) = "OHMCode_WinterWet"
  HeaderVeg_Reqd(cp_OHMCode_WDry) = "OHMCode_WinterDry"
  HeaderVeg_Reqd(cp_cpAnOHM)           = "cp"		! AnOHM TS
  HeaderVeg_Reqd(cp_kkAnOHM)           = "kk"		! AnOHM TS
  HeaderVeg_Reqd(cp_chAnOHM)           = "ch"		! AnOHM TS

  ! ========== SUEWS_Water.txt ==================
  HeaderWater_Reqd(cw_Code)         = "Code"
  HeaderWater_Reqd(cw_AlbMin)       = "AlbedoMin"
  HeaderWater_Reqd(cw_AlbMax)       = "AlbedoMax"
  HeaderWater_Reqd(cw_Emis)         = "Emissivity"
  HeaderWater_Reqd(cw_StorMin)      = "StorageMin"
  HeaderWater_Reqd(cw_StorMax)      = "StorageMax"
  HeaderWater_Reqd(cw_WetThresh)    = "WetThreshold"
  HeaderWater_Reqd(cw_StateLimit)   = "StateLimit"
  HeaderWater_Reqd(cw_DrEq)         = "DrainageEq"
  HeaderWater_Reqd(cw_DrCoef1)      = "DrainageCoef1"
  HeaderWater_Reqd(cw_DrCoef2)      = "DrainageCoef2"
  HeaderWater_Reqd(cw_OHMCode_SWet) = "OHMCode_SummerWet"
  HeaderWater_Reqd(cw_OHMCode_SDry) = "OHMCode_SummerDry"
  HeaderWater_Reqd(cw_OHMCode_WWet) = "OHMCode_WinterWet"
  HeaderWater_Reqd(cw_OHMCode_WDry) = "OHMCode_WinterDry"
  HeaderWater_Reqd(cw_cpAnOHM)           = "cp"		! AnOHM TS
  HeaderWater_Reqd(cw_kkAnOHM)           = "kk"		! AnOHM TS
  HeaderWater_Reqd(cw_chAnOHM)           = "ch"		! AnOHM TS

  ! ========== SUEWS_Snow.txt ===================
  HeaderSnow_Reqd(cs_Code)         = "Code"
  HeaderSnow_Reqd(cs_SnowRMFactor) = "RadMeltFactor"
  HeaderSnow_Reqd(cs_SnowTMFactor) = "TempMeltFactor"
  HeaderSnow_Reqd(cs_SnowAlbMin)   = "AlbedoMin"
  HeaderSnow_Reqd(cs_SnowAlbMax)   = "AlbedoMax"
  HeaderSnow_Reqd(cs_SnowEmis)     = "Emissivity"
  HeaderSnow_Reqd(cs_Snowtau_a)    = "tau_a"
  HeaderSnow_Reqd(cs_Snowtau_f)    = "tau_f"
  HeaderSnow_Reqd(cs_SnowPLimAlb)  = "PrecipLimAlb"
  HeaderSnow_Reqd(cs_SnowSDMin)    = "SnowDensMin"
  HeaderSnow_Reqd(cs_SnowSDMax)    = "SnowDensMax"
  HeaderSnow_Reqd(cs_Snowtau_r)    = "tau_r"
  HeaderSnow_Reqd(cs_SnowCRWMin)   = "CRWMin"
  HeaderSnow_Reqd(cs_SnowCRWMax)   = "CRWMax"
  HeaderSnow_Reqd(cs_SnowPLimSnow) = "PrecipLimSnow"
  HeaderSnow_Reqd(cs_OHMCode_SWet) = "OHMCode_SummerWet"
  HeaderSnow_Reqd(cs_OHMCode_SDry) = "OHMCode_SummerDry"
  HeaderSnow_Reqd(cs_OHMCode_WWet) = "OHMCode_WinterWet"
  HeaderSnow_Reqd(cs_OHMCode_WDry) = "OHMCode_WinterDry"
  HeaderSnow_Reqd(cs_cpAnOHM)           = "cp"		! AnOHM TS
  HeaderSnow_Reqd(cs_kkAnOHM)           = "kk"		! AnOHM TS
  HeaderSnow_Reqd(cs_chAnOHM)           = "ch"		! AnOHM TS


  ! ========== SUEWS_Soil.txt ===================
  HeaderSoil_Reqd(cSo_Code)        = "Code"
  HeaderSoil_Reqd(cSo_SoilDepth)   = "SoilDepth"
  HeaderSoil_Reqd(cSo_SoilStCap)   = "SoilStoreCap"
  HeaderSoil_Reqd(cSo_KSat)        = "SatHydraulicCond"
  HeaderSoil_Reqd(cSo_SoilDens)    = "SoilDensity"
  HeaderSoil_Reqd(cSo_SoilInfRate) = "InfiltrationRate"
  HeaderSoil_Reqd(cSo_ObsSMDepth)  = "OBS_SMDepth"
  HeaderSoil_Reqd(cSo_ObsSMMax)    = "OBS_SMCap"
  HeaderSoil_Reqd(cSo_ObsSNRFrac)  = "OBS_SoilNotRocks"

  ! ========== SUEWS_Conductance.txt ============
  HeaderCond_Reqd(cc_Code)         = "Code"
  HeaderCond_Reqd(cc_GsG1)         = "G1"
  HeaderCond_Reqd(cc_GsG2)         = "G2"
  HeaderCond_Reqd(cc_GsG3)         = "G3"
  HeaderCond_Reqd(cc_GsG4)         = "G4"
  HeaderCond_Reqd(cc_GsG5)         = "G5"
  HeaderCond_Reqd(cc_GsG6)         = "G6"
  HeaderCond_Reqd(cc_GsTH)         = "TH"
  HeaderCond_Reqd(cc_GsTL)         = "TL"
  HeaderCond_Reqd(cc_GsS1)         = "S1"
  HeaderCond_Reqd(cc_GsS2)         = "S2"
  HeaderCond_Reqd(cc_GsKmax)       = "Kmax"

  ! ========== SUEWS_OHMCoefficients.txt ========
  HeaderOHMCoefficients_Reqd(cO_Code) = "Code"
  HeaderOHMCoefficients_Reqd(cO_a1)   = "a1"
  HeaderOHMCoefficients_Reqd(cO_a2)   = "a2"
  HeaderOHMCoefficients_Reqd(cO_a3)   = "a3"

  ! ========== SUEWS_ESTMCoefficients.txt ========
  HeaderESTMCoefficients_Reqd(cE_Code)       = "Code"
  HeaderESTMCoefficients_Reqd(cE_thick1_r)   = "thick1"
  HeaderESTMCoefficients_Reqd(cE_k1_r)       = "k1"
  HeaderESTMCoefficients_Reqd(cE_rhoCp1_r)   = "rhoCp1"
  HeaderESTMCoefficients_Reqd(cE_thick2_r)   = "thick2"
  HeaderESTMCoefficients_Reqd(cE_k2_r)       = "k2"
  HeaderESTMCoefficients_Reqd(cE_rhoCp2_r)   = "rhoCp2"
  HeaderESTMCoefficients_Reqd(cE_thick3_r)   = "thick3"
  HeaderESTMCoefficients_Reqd(cE_k3_r)       = "k3"
  HeaderESTMCoefficients_Reqd(cE_rhoCp3_r)   = "rhoCp3"
  HeaderESTMCoefficients_Reqd(cE_thick4_r)   = "thick4"
  HeaderESTMCoefficients_Reqd(cE_k4_r)       = "k4"
  HeaderESTMCoefficients_Reqd(cE_rhoCp4_r)   = "rhoCp4"
  HeaderESTMCoefficients_Reqd(cE_thick5_r)   = "thick5"
  HeaderESTMCoefficients_Reqd(cE_k5_r)       = "k5"
  HeaderESTMCoefficients_Reqd(cE_rhoCp5_r)   = "rhoCp5"
  HeaderESTMCoefficients_Reqd(cE_thick1_e)   = "thick1"
  HeaderESTMCoefficients_Reqd(cE_k1_e)       = "k1"
  HeaderESTMCoefficients_Reqd(cE_rhoCp1_e)   = "rhoCp1"
  HeaderESTMCoefficients_Reqd(cE_thick2_e)   = "thick2"
  HeaderESTMCoefficients_Reqd(cE_k2_e)       = "k2"
  HeaderESTMCoefficients_Reqd(cE_rhoCp2_e)   = "rhoCp2"
  HeaderESTMCoefficients_Reqd(cE_thick3_e)   = "thick3"
  HeaderESTMCoefficients_Reqd(cE_k3_e)       = "k3"
  HeaderESTMCoefficients_Reqd(cE_rhoCp3_e)   = "rhoCp3"
  HeaderESTMCoefficients_Reqd(cE_thick4_e)   = "thick4"
  HeaderESTMCoefficients_Reqd(cE_k4_e)       = "k4"
  HeaderESTMCoefficients_Reqd(cE_rhoCp4_e)   = "rhoCp4"
  HeaderESTMCoefficients_Reqd(cE_thick5_e)   = "thick5"
  HeaderESTMCoefficients_Reqd(cE_k5_e)       = "k5"
  HeaderESTMCoefficients_Reqd(cE_rhoCp5_e)   = "rhoCp5"
  HeaderESTMCoefficients_Reqd(cE_thick1_i)   = "thick1"
  HeaderESTMCoefficients_Reqd(cE_k1_i)       = "k1"
  HeaderESTMCoefficients_Reqd(cE_rhoCp1_i)   = "rhoCp1"
  HeaderESTMCoefficients_Reqd(cE_thick2_i)   = "thick2"
  HeaderESTMCoefficients_Reqd(cE_k2_i)       = "k2"
  HeaderESTMCoefficients_Reqd(cE_rhoCp2_i)   = "rhoCp2"
  HeaderESTMCoefficients_Reqd(cE_thick3_i)   = "thick3"
  HeaderESTMCoefficients_Reqd(cE_k3_i)       = "k3"
  HeaderESTMCoefficients_Reqd(cE_rhoCp3_i)   = "rhoCp3"
  HeaderESTMCoefficients_Reqd(cE_thick4_i)   = "thick4"
  HeaderESTMCoefficients_Reqd(cE_k4_i)       = "k4"
  HeaderESTMCoefficients_Reqd(cE_rhoCp4_i)   = "rhoCp4"
  HeaderESTMCoefficients_Reqd(cE_thick5_i)   = "thick5"
  HeaderESTMCoefficients_Reqd(cE_k5_i)       = "k5"
  HeaderESTMCoefficients_Reqd(cE_rhoCp5_i)   = "rhoCp5"
  HeaderESTMCoefficients_Reqd(cE_nroom)      = "nroom"
  HeaderESTMCoefficients_Reqd(cE_alb_ibld)   = "alb_ibld"
  HeaderESTMCoefficients_Reqd(cE_em_ibld)    = "em_ibld"
  HeaderESTMCoefficients_Reqd(cE_CH_iwall)   = "CH_iwall"
  HeaderESTMCoefficients_Reqd(cE_CH_iroof)   = "CH_iroof"
  HeaderESTMCoefficients_Reqd(cE_CH_ibld)    = "CH_ibld"
  HeaderESTMCoefficients_Reqd(cE_fwall)      = "fwall"

  ! ========== SUEWS_AnthropogenicHeat.txt ======
  HeaderAnthropogenicHeat_Reqd(cA_Code)     = "Code"
  HeaderAnthropogenicHeat_Reqd(cA_BaseTHDD) = "BaseTHDD"
  HeaderAnthropogenicHeat_Reqd(cA_QF_A1)    = "QF_A_Weekday"
  HeaderAnthropogenicHeat_Reqd(cA_QF_B1)    = "QF_B_Weekday"
  HeaderAnthropogenicHeat_Reqd(cA_QF_C1)    = "QF_C_Weekday"
  HeaderAnthropogenicHeat_Reqd(cA_QF_A2)    = "QF_A_Weekend"
  HeaderAnthropogenicHeat_Reqd(cA_QF_B2)    = "QF_B_Weekend"
  HeaderAnthropogenicHeat_Reqd(cA_QF_C2)    = "QF_C_Weekend"
  HeaderAnthropogenicHeat_Reqd(cA_AHMin)    = "AHMin"
  HeaderAnthropogenicHeat_Reqd(cA_AHSlope)  = "AHSlope"
  HeaderAnthropogenicHeat_Reqd(cA_TCritic)  = "TCritic"

  ! ========== SUEWS_Irrigation.txt =============
  HeaderIrrigation_Reqd(cIr_Code)         = "Code"
  HeaderIrrigation_Reqd(cIr_IeStart)      = "Ie_start"
  HeaderIrrigation_Reqd(cIr_IeEnd)        = "Ie_end"
  HeaderIrrigation_Reqd(cIr_IntWU)        = "InternalWaterUse"
  HeaderIrrigation_Reqd(cIr_Faut)         = "Faut"
  HeaderIrrigation_Reqd(cIr_Ie_a1)        = "Ie_a1"
  HeaderIrrigation_Reqd(cIr_Ie_a2)        = "Ie_a2"
  HeaderIrrigation_Reqd(cIr_Ie_a3)        = "Ie_a3"
  HeaderIrrigation_Reqd(cIr_Ie_m1)        = "Ie_m1"
  HeaderIrrigation_Reqd(cIr_Ie_m2)        = "Ie_m2"
  HeaderIrrigation_Reqd(cIr_Ie_m3)        = "Ie_m3"
  HeaderIrrigation_Reqd(cIr_DayWat1)      = "DayWat(1)"
  HeaderIrrigation_Reqd(cIr_DayWat2)      = "DayWat(2)"
  HeaderIrrigation_Reqd(cIr_DayWat3)      = "DayWat(3)"
  HeaderIrrigation_Reqd(cIr_DayWat4)      = "DayWat(4)"
  HeaderIrrigation_Reqd(cIr_DayWat5)      = "DayWat(5)"
  HeaderIrrigation_Reqd(cIr_DayWat6)      = "DayWat(6)"
  HeaderIrrigation_Reqd(cIr_DayWat7)      = "DayWat(7)"
  HeaderIrrigation_Reqd(cIr_DayWatPer1)   = "DayWatPer(1)"
  HeaderIrrigation_Reqd(cIr_DayWatPer2)   = "DayWatPer(2)"
  HeaderIrrigation_Reqd(cIr_DayWatPer3)   = "DayWatPer(3)"
  HeaderIrrigation_Reqd(cIr_DayWatPer4)   = "DayWatPer(4)"
  HeaderIrrigation_Reqd(cIr_DayWatPer5)   = "DayWatPer(5)"
  HeaderIrrigation_Reqd(cIr_DayWatPer6)   = "DayWatPer(6)"
  HeaderIrrigation_Reqd(cIr_DayWatPer7)   = "DayWatPer(7)"

  ! ========== SUEWS_Profiles.txt ===============
  HeaderProfiles_Reqd(cPr_Code)      = "Code"
  HeaderProfiles_Reqd(cPr_Hours( 1)) = "0"
  HeaderProfiles_Reqd(cPr_Hours( 2)) = "1"
  HeaderProfiles_Reqd(cPr_Hours( 3)) = "2"
  HeaderProfiles_Reqd(cPr_Hours( 4)) = "3"
  HeaderProfiles_Reqd(cPr_Hours( 5)) = "4"
  HeaderProfiles_Reqd(cPr_Hours( 6)) = "5"
  HeaderProfiles_Reqd(cPr_Hours( 7)) = "6"
  HeaderProfiles_Reqd(cPr_Hours( 8)) = "7"
  HeaderProfiles_Reqd(cPr_Hours( 9)) = "8"
  HeaderProfiles_Reqd(cPr_Hours(10)) = "9"
  HeaderProfiles_Reqd(cPr_Hours(11)) = "10"
  HeaderProfiles_Reqd(cPr_Hours(12)) = "11"
  HeaderProfiles_Reqd(cPr_Hours(13)) = "12"
  HeaderProfiles_Reqd(cPr_Hours(14)) = "13"
  HeaderProfiles_Reqd(cPr_Hours(15)) = "14"
  HeaderProfiles_Reqd(cPr_Hours(16)) = "15"
  HeaderProfiles_Reqd(cPr_Hours(17)) = "16"
  HeaderProfiles_Reqd(cPr_Hours(18)) = "17"
  HeaderProfiles_Reqd(cPr_Hours(19)) = "18"
  HeaderProfiles_Reqd(cPr_Hours(20)) = "19"
  HeaderProfiles_Reqd(cPr_Hours(21)) = "20"
  HeaderProfiles_Reqd(cPr_Hours(22)) = "21"
  HeaderProfiles_Reqd(cPr_Hours(23)) = "22"
  HeaderProfiles_Reqd(cPr_Hours(24)) = "23"

  ! ========== SUEWS_WithinGridWaterDist.txt ====
  HeaderWGWaterDist_Reqd(cWG_Code)        = "Code"
  HeaderWGWaterDist_Reqd(cWG_ToPaved)     = "ToPaved"
  HeaderWGWaterDist_Reqd(cWG_ToBldgs)     = "ToBldgs"
  HeaderWGWaterDist_Reqd(cWG_ToEveTr)     = "ToEveTr"
  HeaderWGWaterDist_Reqd(cWG_ToDecTr)     = "ToDecTr"
  HeaderWGWaterDist_Reqd(cWG_ToGrass)     = "ToGrass"
  HeaderWGWaterDist_Reqd(cWG_ToBSoil)     = "ToBSoil"
  HeaderWGWaterDist_Reqd(cWG_ToWater)     = "ToWater"
  HeaderWGWaterDist_Reqd(cWG_ToRunoff)    = "ToRunoff"
  HeaderWGWaterDist_Reqd(cWG_ToSoilStore) = "ToSoilStore"

  ! =======================================================


  !write(*,*) 'Checking header for ', FileName
  ! Check columns in input files match model code

  IF(FileName == 'SUEWS_NonVeg.txt') THEN
     IF(ANY(HeaderNonVeg_File /= HeaderNonVeg_Reqd)) THEN
        WRITE(*,*) HeaderNonVeg_File == HeaderNonVeg_Reqd
        WRITE(*,*) HeaderNonVeg_File
        WRITE(*,*) HeaderNonVeg_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_NonVeg.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Veg.txt') THEN
     IF(ANY(HeaderVeg_File /= HeaderVeg_Reqd)) THEN
        WRITE(*,*) HeaderVeg_File == HeaderVeg_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Veg.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Water.txt') THEN
     IF(ANY(HeaderWater_File /= HeaderWater_Reqd)) THEN
        WRITE(*,*) HeaderWater_File == HeaderWater_Reqd
        WRITE(*,*) HeaderWater_File
        WRITE(*,*) HeaderWater_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Water.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Snow.txt') THEN
     IF(ANY(HeaderSnow_File /= HeaderSnow_Reqd)) THEN
        WRITE(*,*) HeaderSnow_File == HeaderSnow_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Snow.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Soil.txt') THEN
     IF(ANY(HeaderSoil_File /= HeaderSoil_Reqd)) THEN
        WRITE(*,*) HeaderSoil_File == HeaderSoil_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Soil.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Conductance.txt') THEN
     IF(ANY(HeaderCond_File /= HeaderCond_Reqd)) THEN
        WRITE(*,*) HeaderCond_File == HeaderCond_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Cond.txt does not match model code.',notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_OHMCoefficients.txt') THEN
     IF(ANY(HeaderOHMCoefficients_File /= HeaderOHMCoefficients_Reqd)) THEN
        WRITE(*,*) HeaderOHMCoefficients_File == HeaderOHMCoefficients_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_OHMCoefficients.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_ESTMCoefficients.txt') THEN
     IF(ANY(HeaderESTMCoefficients_File /= HeaderESTMCoefficients_Reqd)) THEN
        WRITE(*,*) HeaderESTMCoefficients_File == HeaderESTMCoefficients_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_ESTMCoefficients.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_AnthropogenicHeat.txt') THEN
     IF(ANY(HeaderAnthropogenicHeat_File /= HeaderAnthropogenicHeat_Reqd)) THEN
        WRITE(*,*) HeaderAnthropogenicHeat_File == HeaderAnthropogenicHeat_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_AnthropogenicHeat.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Irrigation.txt') THEN
     IF(ANY(HeaderIrrigation_File /= HeaderIrrigation_Reqd)) THEN
        WRITE(*,*) HeaderIrrigation_File == HeaderIrrigation_Reqd
        WRITE(*,*) HeaderIrrigation_File
        WRITE(*,*) HeaderIrrigation_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Irrigation.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_Profiles.txt') THEN
     IF(ANY(HeaderProfiles_File /= HeaderProfiles_Reqd)) THEN
        WRITE(*,*) HeaderProfiles_File == HeaderProfiles_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_Profiles.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSEIF(FileName == 'SUEWS_WithinGridWaterDist.txt') THEN
     IF(ANY(HeaderWGWaterDist_File /= HeaderWGWaterDist_Reqd)) THEN
        WRITE(*,*) HeaderWGWaterDist_File == HeaderWGWaterDist_Reqd
        CALL ErrorHint(56,'Names or order of columns in SUEWS_WithinGridWaterDist.txt does not match model code.',&
             notUsed,notUsed,notUsedI)
     ENDIF

  ELSE
     WRITE(*,*) 'Problem in subroutine InputHeaderCheck. File header not specified in model code for ',FileName
     CALL ErrorHint(58,FileName,notUsed,notUsed,notUsedI)

  ENDIF

ENDSUBROUTINE InputHeaderCheck

!-------------------------------------------------------------------------
