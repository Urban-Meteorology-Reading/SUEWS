!-------------------------------------------------------------------------
SUBROUTINE InputHeaderCheck(FileName)
! Checks columns in input files match the columns expected by model code 
! Model code columns are defined here
! HCW 12 Nov 2014
!-------------------------------------------------------------------------

  use allocateArray 
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  character (len=50):: FileName
      
  ! ========== Define expected column names here ==========
  ! =======================================================
 
  ! ========== SUEWS_NonVeg.txt =============
  HeaderNonVeg_Reqd(ci_Code)	  = "Code"
  HeaderNonVeg_Reqd(ci_AlbMin)	  = "AlbedoMin"
  HeaderNonVeg_Reqd(ci_AlbMax)	  = "AlbedoMax"
  HeaderNonVeg_Reqd(ci_Emis)	  = "Emissivity"
  HeaderNonVeg_Reqd(ci_StorMin)	  = "StorageMin"
  HeaderNonVeg_Reqd(ci_StorMax)	  = "StorageMax"
  HeaderNonVeg_Reqd(ci_WetThresh)= "WetThreshold"
  HeaderNonVeg_Reqd(ci_StateLimit)= "StateLimit"
  HeaderNonVeg_Reqd(ci_DrEq)      = "DrainageEq"
  HeaderNonVeg_Reqd(ci_DrCoef1)   = "DrainageCoef1"
  HeaderNonVeg_Reqd(ci_DrCoef2)   = "DrainageCoef2"
  HeaderNonVeg_Reqd(ci_SoilTCode)    = "SoilTypeCode"
  HeaderNonVeg_Reqd(ci_SnowLimPat)   = "SnowLimPatch"
  HeaderNonVeg_Reqd(ci_SnowLimRem)   = "SnowLimRemove"     
  HeaderNonVeg_Reqd(ci_OHMCode_SWet) = "OHMCode_SummerWet"     
  HeaderNonVeg_Reqd(ci_OHMCode_SDry) = "OHMCode_SummerDry"     
  HeaderNonVeg_Reqd(ci_OHMCode_WWet) = "OHMCode_WinterWet"     
  HeaderNonVeg_Reqd(ci_OHMCode_WDry) = "OHMCode_WinterDry"     
   
  ! ========== SUEWS_Veg.txt ===============
  HeaderVeg_Reqd(cp_Code)	  = "Code"
  HeaderVeg_Reqd(cp_AlbMin)	  = "AlbedoMin"
  HeaderVeg_Reqd(cp_AlbMax)	  = "AlbedoMax"
  HeaderVeg_Reqd(cp_Emis)	  = "Emissivity"
  HeaderVeg_Reqd(cp_StorMin)	  = "StorageMin"
  HeaderVeg_Reqd(cp_StorMax)	  = "StorageMax"
  HeaderVeg_Reqd(cp_WetThresh)    = "WetThreshold"
  HeaderVeg_Reqd(cp_StateLimit)   = "StateLimit"
  HeaderVeg_Reqd(cp_DrEq)	  = "DrainageEq"
  HeaderVeg_Reqd(cp_DrCoef1)	  = "DrainageCoef1"
  HeaderVeg_Reqd(cp_DrCoef2)	  = "DrainageCoef2"
  HeaderVeg_Reqd(cp_SoilTCode)	  = "SoilTypeCode"
  HeaderVeg_Reqd(cp_SnowLimPat)   = "SnowLimPatch"
  HeaderVeg_Reqd(cp_BaseT)	  = "BaseT"      
  HeaderVeg_Reqd(cp_BaseTe)	  = "BaseTe"      
  HeaderVeg_Reqd(cp_GDDFull)	  = "GDDFull"      
  HeaderVeg_Reqd(cp_SDDFull)	  = "SDDFull"      
  HeaderVeg_Reqd(cp_LAIMin)	  = "LAIMin"      
  HeaderVeg_Reqd(cp_LAIMax)	  = "LAIMax"      
  HeaderVeg_Reqd(cp_GsMax)	  = "MaxConductance"      
  HeaderVeg_Reqd(cp_LAIEq)	  = "LAIEq"      
  HeaderVeg_Reqd(cp_LeafGP1)	  = "LeafGrowthPower1"      
  HeaderVeg_Reqd(cp_LeafGP2)	  = "LeafGrowthPower2"      
  HeaderVeg_Reqd(cp_LeafOP1)	  = "LeafOffPower1"      
  HeaderVeg_Reqd(cp_LeafOP2)	  = "LeafOffPower2"    
  HeaderVeg_Reqd(cp_OHMCode_SWet) = "OHMCode_SummerWet"     
  HeaderVeg_Reqd(cp_OHMCode_SDry) = "OHMCode_SummerDry"     
  HeaderVeg_Reqd(cp_OHMCode_WWet) = "OHMCode_WinterWet"     
  HeaderVeg_Reqd(cp_OHMCode_WDry) = "OHMCode_WinterDry"        
   
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
  HeaderWater_Reqd(cw_DrCoef2) 	    = "DrainageCoef2"   
  HeaderWater_Reqd(cw_OHMCode_SWet) = "OHMCode_SummerWet"     
  HeaderWater_Reqd(cw_OHMCode_SDry) = "OHMCode_SummerDry"     
  HeaderWater_Reqd(cw_OHMCode_WWet) = "OHMCode_WinterWet"     
  HeaderWater_Reqd(cw_OHMCode_WDry) = "OHMCode_WinterDry"       
   
  ! ========== SUEWS_Snow.txt ===================
  HeaderSnow_Reqd(cs_Code) 	   = "Code"
  HeaderSnow_Reqd(cs_SnowRMFactor) = "RadMeltFactor"
  HeaderSnow_Reqd(cs_SnowTMFactor) = "TempMeltFactor"
  HeaderSnow_Reqd(cs_SnowAlbMin)   = "AlbedoMin"
  HeaderSnow_Reqd(cs_SnowAlbMax)   = "AlbedoMax"
  HeaderSnow_Reqd(cs_SnowAlb)	   = "Albedo"
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
  
  ! ========== SUEWS_Soil.txt ===================
  HeaderSoil_Reqd(cSo_Code)        = "Code" 
  HeaderSoil_Reqd(cSo_SoilDepth)   = "SoilDepth" 
  HeaderSoil_Reqd(cSo_SoilStCap)   = "SoilStoreCap" 
  HeaderSoil_Reqd(cSo_KSat)  	   = "SatHydraulicCond" 
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
   
  if(FileName == 'SUEWS_NonVeg.txt') then
     if(ANY(HeaderNonVeg_File /= HeaderNonVeg_Reqd)) then
        write(*,*) HeaderNonVeg_File == HeaderNonVeg_Reqd
        write(*,*) HeaderNonVeg_File
        write(*,*) HeaderNonVeg_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_NonVeg.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif  
    
  elseif(FileName == 'SUEWS_Veg.txt') then
     if(ANY(HeaderVeg_File /= HeaderVeg_Reqd)) then
        write(*,*) HeaderVeg_File == HeaderVeg_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Veg.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif 
    
  elseif(FileName == 'SUEWS_Water.txt') then
     if(ANY(HeaderWater_File /= HeaderWater_Reqd)) then
        write(*,*) HeaderWater_File == HeaderWater_Reqd
        write(*,*) HeaderWater_File
        write(*,*) HeaderWater_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Water.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif 
  
  elseif(FileName == 'SUEWS_Snow.txt') then
     if(ANY(HeaderSnow_File /= HeaderSnow_Reqd)) then
        write(*,*) HeaderSnow_File == HeaderSnow_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Snow.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif 
   
  elseif(FileName == 'SUEWS_Soil.txt') then 
     if(ANY(HeaderSoil_File /= HeaderSoil_Reqd)) then
        write(*,*) HeaderSoil_File == HeaderSoil_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Soil.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif 
    
  elseif(FileName == 'SUEWS_Conductance.txt') then
     if(ANY(HeaderCond_File /= HeaderCond_Reqd)) then
         write(*,*) HeaderCond_File == HeaderCond_Reqd
         call ErrorHint(56,'Names or order of columns in SUEWS_Cond.txt does not match model code.',notUsed,notUsed,notUsedI)
     endif 
   
  elseif(FileName == 'SUEWS_OHMCoefficients.txt') then 
     if(ANY(HeaderOHMCoefficients_File /= HeaderOHMCoefficients_Reqd)) then
        write(*,*) HeaderOHMCoefficients_File == HeaderOHMCoefficients_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_OHMCoefficients.txt does not match model code.',&
                       notUsed,notUsed,notUsedI)
     endif 
   
  elseif(FileName == 'SUEWS_AnthropogenicHeat.txt') then
     if(ANY(HeaderAnthropogenicHeat_File /= HeaderAnthropogenicHeat_Reqd)) then
        write(*,*) HeaderAnthropogenicHeat_File == HeaderAnthropogenicHeat_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_AnthropogenicHeat.txt does not match model code.',&
                       notUsed,notUsed,notUsedI)
     endif    
   
  elseif(FileName == 'SUEWS_Irrigation.txt') then 
     if(ANY(HeaderIrrigation_File /= HeaderIrrigation_Reqd)) then
        write(*,*) HeaderIrrigation_File == HeaderIrrigation_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Irrigation.txt does not match model code.',&
                       notUsed,notUsed,notUsedI)
     endif    
   
  elseif(FileName == 'SUEWS_Profiles.txt') then
     if(ANY(HeaderProfiles_File /= HeaderProfiles_Reqd)) then
        write(*,*) HeaderProfiles_File == HeaderProfiles_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_Profiles.txt does not match model code.',&
                       notUsed,notUsed,notUsedI)
      endif    
   
  elseif(FileName == 'SUEWS_WithinGridWaterDist.txt') then 
     if(ANY(HeaderWGWaterDist_File /= HeaderWGWaterDist_Reqd)) then
        write(*,*) HeaderWGWaterDist_File == HeaderWGWaterDist_Reqd
        call ErrorHint(56,'Names or order of columns in SUEWS_WithinGridWaterDist.txt does not match model code.',&
                       notUsed,notUsed,notUsedI)
     endif
    
  else
     write(*,*) 'Problem in subroutine InputHeaderCheck. File header not specified in model code for ',FileName
     call ErrorHint(58,FileName,notUsed,notUsed,notUsedI)
  
  endif
   
ENDSUBROUTINE InputHeaderCheck

!-------------------------------------------------------------------------
