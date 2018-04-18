! Subroutines for matching codes in the input files
!  could re-write as a generic function later...

 SUBROUTINE CodeMatchOHM(Gridiv,is,SWWD) 
! Matches OHM coefficients a1, a2, a3 via OHM codes 
! for summer/winter wet/dry conditions 
! HCW 03 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: gridiv
  integer:: is
  character(len=4):: SWWD
 
  iv5=0 ! Reset iv5 to zero 
 
  if(SWWD == 'SWet') then
    
       do iv5=1,nlinesOHMCoefficients
          if (OHMCoefficients_Coeff(iv5,cO_Code)==SurfaceChar(gridiv,c_OHMCode_SWet(is))) then
          exit
          elseif(iv5 == nlinesOHMCoefficients) then 
          write(*,*) 'Program stopped! OHM code (summer wet)',SurfaceChar(gridiv,c_OHMCode_SWet(is)),&
         		   'not found in OHM_Coefficients.txt for surface',is,'.'
          call ErrorHint(57,'Cannot find OHM code (summer wet)',SurfaceChar(gridiv,c_OHMCode_SWet(is)),notUsed,notUsedI)
          endif
     enddo    
     
   elseif(SWWD == 'SDry') then
    
       do iv5=1,nlinesOHMCoefficients
          if (OHMCoefficients_Coeff(iv5,cO_Code)==SurfaceChar(gridiv,c_OHMCode_SDry(is))) then
          exit
          elseif(iv5 == nlinesOHMCoefficients) then 
          write(*,*) 'Program stopped! OHM code (summer dry)',SurfaceChar(gridiv,c_OHMCode_SDry(is)),&
         		   'not found in OHM_Coefficients.txt for surface',is,'.'
          call ErrorHint(57,'Cannot find OHM code (summer dry)',SurfaceChar(gridiv,c_OHMCode_SDry(is)),notUsed,notUsedI)
          endif
     enddo         
     
  elseif(SWWD == 'WWet') then
  
     do iv5=1,nlinesOHMCoefficients
        if (OHMCoefficients_Coeff(iv5,cO_Code)==SurfaceChar(gridiv,c_OHMCode_WWet(is))) then
        exit
        elseif(iv5 == nlinesOHMCoefficients) then 
        write(*,*) 'Program stopped! OHM code (winter wet)',SurfaceChar(gridiv,c_OHMCode_WWet(is)),&
       		   'not found in OHM_Coefficients.txt for surface',is,'.'
        call ErrorHint(57,'Cannot find OHM code (winter wet)',SurfaceChar(gridiv,c_OHMCode_WWet(is)),notUsed,notUsedI)
        endif
     enddo        

  elseif(SWWD == 'WDry') then
  
     do iv5=1,nlinesOHMCoefficients
        if (OHMCoefficients_Coeff(iv5,cO_Code)==SurfaceChar(gridiv,c_OHMCode_WDry(is))) then
        exit
        elseif(iv5 == nlinesOHMCoefficients) then 
        write(*,*) 'Program stopped! OHM code (winter dry)',SurfaceChar(gridiv,c_OHMCode_WDry(is)),&
       		   'not found in OHM_Coefficients.txt for surface',is,'.'
        call ErrorHint(57,'Cannot find OHM code (winter dry)',SurfaceChar(gridiv,c_OHMCode_WDry(is)),notUsed,notUsedI)
        endif
     enddo             
  
  else
     write(*,*) 'Problem with CodeMatchOHM (in SUEWS_CodeMatch.f95). ',SWWD,' not recognised. Needs to be one of: ',&
     	        'SWet = Summer Wet, SDry = Summer Dry, WWet = WinterWet, WDry = Winter Dry. N.B. Case sensitive. '
     stop
  endif
   
  return
ENDSUBROUTINE CodeMatchOHM
! ---------------------------------------------------------
    
 SUBROUTINE CodeMatchESTM(Gridiv,is) 
! Matches ESTM coefficients via ESTM code 
! Modified HCW 16 Jun 2016 - for SUEWS surface types
!                          - removed summer/winter wet/dry option   
! S.O. 04 Feb 2016
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: gridiv
  integer:: is
   
  iv5=0 ! Reset iv5 to zero
  
  do iv5=1,nlinesESTMCoefficients
     if (ESTMCoefficients_Coeff(iv5,cE_Code)==SurfaceChar(gridiv,c_ESTMCode(is))) then
        exit
     elseif(iv5 == nlinesESTMCoefficients) then 
        write(*,*) 'Program stopped! ESTM code',SurfaceChar(gridiv,c_ESTMCode(is)),&
                   'not found in ESTM_Coefficients.txt for surface',is,'.'
        call ErrorHint(57,'Cannot find ESTM code',SurfaceChar(gridiv,c_ESTMCode(is)),notUsed,notUsedI)
     endif
  enddo      
           
  return
 ENDSUBROUTINE CodeMatchESTM
! ---------------------------------------------------------    

 SUBROUTINE CodeMatchESTM_Class(Gridiv,is,ii) 
! Matches ESTM coefficients via ESTM codes in SiteSelect for Paved and Bldgs ESTM classes
! HCW 16 Jun 2016
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: gridiv
  integer:: is, ii
   
  iv5=0 ! Reset iv5 to zero
  
  if(is == BldgSurf) then      
     do iv5=1,nlinesESTMCoefficients
        if (ESTMCoefficients_Coeff(iv5,cE_Code)==SurfaceChar(gridiv,c_Code_ESTMClass_Bldgs(ii))) then
           exit
        elseif(iv5 == nlinesESTMCoefficients) then 
           write(*,*) 'Program stopped! ESTM code',SurfaceChar(gridiv,c_Code_ESTMClass_Bldgs(ii)),&
                      'not found in ESTM_Coefficients.txt for surface',is,'.'
           call ErrorHint(57,'Cannot find ESTM code',SurfaceChar(gridiv,c_Code_ESTMClass_Bldgs(ii)),notUsed,notUsedI)
        endif
     enddo      
  elseif(is == PavSurf) then         
     do iv5=1,nlinesESTMCoefficients
        if (ESTMCoefficients_Coeff(iv5,cE_Code)==SurfaceChar(gridiv,c_Code_ESTMClass_Paved(ii))) then
           exit
        elseif(iv5 == nlinesESTMCoefficients) then 
           write(*,*) 'Program stopped! ESTM code',SurfaceChar(gridiv,c_Code_ESTMClass_Paved(ii)),&
                      'not found in ESTM_Coefficients.txt for surface',is,'.'
           call ErrorHint(57,'Cannot find ESTM code',SurfaceChar(gridiv,c_Code_ESTMClass_Paved(ii)),notUsed,notUsedI)
        endif
     enddo      
  else
     write(*,*) 'Problem with CodeMatchESTM_Class (in SUEWS_CodeMatch.f95). ',is,' not correct. Needs to be either ',&
     	        '1 = Paved surfaced, 2 = Bldgs surfaces.'
     stop
  endif   
  return
 ENDSUBROUTINE CodeMatchESTM_Class
! ---------------------------------------------------------    

SUBROUTINE CodeMatchProf(Gridiv,SurfaceCharCodeCol)
! Matches Soil characteristics via codes in *SurfaceChar*
 ! for energy use/water use/snow clearing
! HCW 20 Nov 2014
 ! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: Gridiv
  integer:: SurfaceCharCodeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesProfiles
     if (Profiles_Coeff(iv5,cPr_Code)==SurfaceChar(Gridiv,SurfaceCharCodeCol)) then
     exit
     elseif(iv5 == nlinesProfiles) then 
     write(*,*) 'Program stopped! Profile code ',SurfaceChar(Gridiv,SurfaceCharCodeCol),'not found in SUEWS_Profiles.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Profiles.txt',SurfaceChar(Gridiv,SurfaceCharCodeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchProf 
! ---------------------------------------------------------

SUBROUTINE CodeMatchDist(rr,CodeCol,codeColSameSurf)
! Matches within-grid water distribution via codes 
! Checks water cannot flow from one surface to the same surface
! Checks water distribution fractions sum to 1.
! HCW 10 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol, codeColSameSurf
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesWGWaterDist
     if (WGWaterDist_Coeff(iv5,cWG_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesWGWaterDist) then 
     write(*,*) 'Program stopped! Within-grid water distribution code ',SiteSelect(rr,codeCol),&
     		'not found in SUEWS_WaterDistWithinGrid.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_WaterDistWithinGrid.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  ! Check water flow to same surface is zero (previously in RunControlByGridByYear in SUEWS_Initial.f95)
  if(WGWaterDist_Coeff(iv5,codeColSameSurf) /= 0) then
     call ErrorHint(8,'Diagonal elements should be zero as water cannot move from one surface to the same surface.', &
                    WGWaterDist_Coeff(iv5,codeColSameSurf),notUsed,notUsedI)
  endif 
  
  !! MODIFY THIS??
  ! Check water either moves to runoff or soilstore, but not to both
  ! Model returns an error if both ToRunoff and ToSoilStore are non-zero. 
  !! - Probably should remove this...
  !! - Also look at SUEWS_translate, as the non-zero value goes into WaterDist
  if(WGWaterDist_Coeff(iv5,cWG_ToRunoff)/=0.and.WGWaterDist_Coeff(iv5,cWG_ToSoilStore)/=0) then
     call ErrorHint(9,'One of these (ToRunoff,ToSoilStore) should be zero.', &
                    WGWaterDist_Coeff(iv5,cWG_ToRunoff),WGWaterDist_Coeff(iv5,cWG_ToSoilStore),notUsedI)
  endif
  
  !! Also do for water surface once implemented
  if(codeCol /= c_WGWaterCode) then   ! Except for Water surface
     ! Check total water distribution from each surface adds up to 1
     if(sum(WGWaterDist_Coeff(iv5,cWG_ToPaved:cWG_ToSoilStore)) > 1.0000001.or.sum(WGWaterDist_Coeff(iv5,&
             cWG_ToPaved:cWG_ToSoilStore)) < 0.9999999 ) then
        call ErrorHint(8,'Total water distribution from each surface should add up to 1.',&
                       sum(WGWaterDist_Coeff(iv5,cWG_ToPaved:cWG_ToSoilStore)),notUsed,notUsedI)
     endif     
  endif
  
  return
ENDSUBROUTINE CodeMatchDist 
! ---------------------------------------------------------


 SUBROUTINE CodeMatchNonVeg(rr,CodeCol)
! Matches Impervious characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero

  do iv5=1,nlinesNonVeg
     if (NonVeg_Coeff(iv5,ci_Code)==SiteSelect(rr,codeCol)) then
        exit
     elseif(iv5 == nlinesNonVeg) then 
        write(*,*) 'Program stopped! NonVeg code ',SiteSelect(rr,codeCol),'not found in SUEWS_NonVeg.txt.'
        call ErrorHint(57,'Cannot find code in SUEWS_NonVeg.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
 ENDSUBROUTINE CodeMatchNonVeg
! ---------------------------------------------------------   


 SUBROUTINE CodeMatchVeg(rr,CodeCol)
  ! Matches Pervious characteristics via codes in SiteSelect
  ! HCW 20 Nov 2014
  ! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesVeg
     if (Veg_Coeff(iv5,cp_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesVeg) then 
     write(*,*) 'Program stopped! Veg code ',SiteSelect(rr,codeCol),'not found in SUEWS_Vegs.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Veg.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchVeg 
! --------------------------------------------------------- 


SUBROUTINE CodeMatchWater(rr,CodeCol)
! Matches Water characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesWater
     if (Water_Coeff(iv5,cw_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesWater) then 
     write(*,*) 'Program stopped! Water code ',SiteSelect(rr,codeCol),'not found in SUEWS_Water.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Water.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchWater 
! --------------------------------------------------------- 


SUBROUTINE CodeMatchSnow(rr,CodeCol)
! Matches Snow characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesSnow
     if (Snow_Coeff(iv5,cs_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesSnow) then 
     write(*,*) 'Program stopped! Snow code ',SiteSelect(rr,codeCol),'not found in SUEWS_Snow.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Snow.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchSnow
! --------------------------------------------------------- 

SUBROUTINE CodeMatchConductance(rr,CodeCol)
! Matches Conductance characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesConductance
     if (Conductance_Coeff(iv5,cc_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesConductance) then 
     write(*,*) 'Program stopped! Conductance code ',SiteSelect(rr,codeCol),'not found in SUEWS_Conductance.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Conductance.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchConductance
! --------------------------------------------------------- 
 
 
SUBROUTINE CodeMatchAnthropogenic(rr,CodeCol)
! Matches AnthropogenicHeat characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! MH 21 Jun 2017
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesAnthropogenic
     if (Anthropogenic_Coeff(iv5,cA_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesAnthropogenic) then 
     write(*,*) 'Program stopped! Anthropogenic code ',SiteSelect(rr,codeCol),'not found in SUEWS_AnthropogenicHeat.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_AnthropogenicHeat.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchAnthropogenic
! --------------------------------------------------------- 


SUBROUTINE CodeMatchIrrigation(rr,CodeCol)
! Matches Irrigation characteristics via codes in SiteSelect
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: rr
  integer:: codeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesIrrigation
     if (Irrigation_Coeff(iv5,cIr_Code)==SiteSelect(rr,codeCol)) then
     exit
     elseif(iv5 == nlinesIrrigation) then 
     write(*,*) 'Program stopped! Irrigation code ',SiteSelect(rr,codeCol),'not found in SUEWS_Irrigation.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Irrigation.txt',SiteSelect(rr,codeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchIrrigation
! --------------------------------------------------------- 

SUBROUTINE CodeMatchSoil(Gridiv,SurfaceCharCodeCol)
! Matches Soil characteristics via codes in *SurfaceChar*
! HCW 20 Nov 2014
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: Gridiv
  integer:: SurfaceCharCodeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesSoil
     if (Soil_Coeff(iv5,cSo_Code)==SurfaceChar(Gridiv,SurfaceCharCodeCol)) then
     exit
     elseif(iv5 == nlinesSoil) then 
     write(*,*) 'Program stopped! Soil code ',SurfaceChar(Gridiv,SurfaceCharCodeCol),'not found in SUEWS_Soil.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_Soil.txt',SurfaceChar(Gridiv,SurfaceCharCodeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchSoil
! --------------------------------------------------------- 
	
SUBROUTINE CodeMatchBiogen(Gridiv,SurfaceCharCodeCol)
! Matches Biogen characteristics via codes in *SuraceChar*
! MH 16 Jun 2017
! ---------------------------------------------------------

  use allocateArray
  use Initial
  use ColNamesInputFiles
  use defaultNotUsed
  
  IMPLICIT NONE
  
  integer:: Gridiv
  integer:: SurfaceCharCodeCol
 
  iv5=0 ! Reset iv5 to zero
 
  do iv5=1,nlinesBiogen
     if (Biogen_Coeff(iv5,cB_Code)==SurfaceChar(Gridiv,SurfaceCharCodeCol)) then
     exit
     elseif(iv5 == nlinesBiogen) then 
     write(*,*) 'Program stopped! Biogen code ',SurfaceChar(Gridiv,SurfaceCharCodeCol),'not found in SUEWS_BiogenCO2.txt.'
     call ErrorHint(57,'Cannot find code in SUEWS_BiogenCO2.txt',SurfaceChar(Gridiv,SurfaceCharCodeCol),notUsed,notUsedI)
     endif
  enddo   
  
  return
ENDSUBROUTINE CodeMatchBiogen





! --------------------------------------------------------- 
