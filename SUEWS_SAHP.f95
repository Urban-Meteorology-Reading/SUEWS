
!===================================================================================
!Simple Anthropogenic Heat Parameterization routines
!Last modified
! HCW 24 Feb 2017 - Added new anthropogenic heat flux calculation (AnthropHeatMethod=3)
! HCW 25 Aug 2016 - Outputs base QF (part without temp. dependence)
! LJ 27 Jan 2016  - Removal of Tabs
! HCW 20 Jan 2015 - v2015 applies a profile at each model timestep
!                   these have been interpolated from the hourly profile input data (SUEWS_Profiles)
!                   vthey are now normalised (sum to 1) in InitializeSurfaceCharacteristics
!                   N.B. previous versions were not applying weekday/weekend profiles correctly
! 
! AnthropHeatMethod = 1 - Method according to Loridan et al. (2011) : SAHP
! AnthropHeatMethod = 2 - Method according to Jarvi et al. (2011)   : SAHP_2
! AnthropHeatMethod = 0 - Use values in met forcing file, or default QF
!
!===================================================================================

!-----------------------------------------------------------------------------------
 subroutine SAHP_1(QF_o,QF_o_base,QF_o_heat,id,ih,imin)
 ! Called if AnthropHeatMethod = 1
 ! Method according to Loridan et al. (2011)
 ! Weekday/weekend differences due to profile only

  use allocateArray
  use data_in
  use sues_data
  
  IMPLICIT NONE
  
  integer:: id,&   !Day
            ih,&   !Hour, with daylight saving accounted for (ih, not it)
            imin,& !Minute
            iu     !1=weekday OR 2=weekend
  	    
  real (kind(1d0)):: QF_o   !Output: modelled QF [W m-2]
  real (kind(1d0)):: QF_o_base, QF_o_heat  !Output: temperature-independent part and heating only part of modelled QF  [W m-2]

  iu=1     !Set to 1=weekday
  if(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) then  
     iu=2  !Set to 2=weekend
  endif
  
  ! Linear relation with air temperature from Loridan et al. (2011) JAMC Eq 13
  if(Temp_C.lt.T_CRITIC) then    
     QF_o = AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(AH_MIN + AH_SLOPE*(T_CRITIC-Temp_C))         
  else
     QF_o = AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*AH_MIN
  endif 
  
  QF_o_base = AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*AH_MIN
  QF_o_heat = QF_o - QF_o_base
     
 endsubroutine SAHP_1
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
subroutine SAHP_2(QF_o,QF_o_base,QF_o_heat,id,ih,imin)
! Called if AnthropHeatMethod = 2
! Method according to Jarvi et al. (2011)  
! Weekday/weekend differences due to profile and coefficients QF_a,b,c
  
  use allocateArray
  use data_in
  use sues_data
  
  IMPLICIT NONE
   
  integer:: id,&   !Day
            ih,&   !Hour, with daylight saving accounted for (ih, not it)
            imin,& !Minute
            iu     !1=weekday OR 2=weekend
   	    
  real (kind(1d0)):: QF_o   !Output: modelled QF  [W m-2]
  real (kind(1d0)):: QF_o_base, QF_o_heat  !Output: temperature-independent part and heating only part of modelled QF  [W m-2]

  iu=1     !Set to 1=weekday
  if(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) then  
     iu=2  !Set to 2=weekend
  endif

  ! Uses HDD and CDD, see Jarvi et al. (2011) JH Eq 3 x pop density

  !write(*,*) HDD(id-1,1), HDD(id-1,2)

  !write(*,*) '--------- QF ---------'
  !write(*,*) id, ih, imin
  !write(*,*) (NSH*(ih+1-1)+imin*NSH/60+1)
  !write(*,*) AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
  !write(*,*) (Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1))*numCapita

  
  QF_o = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1)))*numCapita

  QF_o_base = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_a(iu)))*numCapita
  QF_o_heat = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_c(iu)*HDD(id-1,1)))*numCapita

  !write(*,*) QF_o
  !write(*,*) '----------------------'
  !write(*,*) " " 

endsubroutine SAHP_2
!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------
SUBROUTINE SAHP_3(QF_o,id,ih,imin)
! Called if AnthropHeatMethod = 3
! Method according to !!!
! Created 24 Feb 2017, by HCW


  USE AllocateArray
  USE data_in
  USE sues_data
  
  IMPLICIT NONE
   
  INTEGER:: id,&   !Day
            ih,&   !Hour, with daylight saving accounted for (ih, not it)
            imin,& !Minute
            iu     !1=weekday OR 2=weekend
   	    
  REAL(KIND(1d0)):: QF_o   !Output: modelled QF  [W m-2]
 
  REAL(KIND(1d0)):: QF_metab, QF_traff, QF_build    !W m-2
  REAL(KIND(1d0)):: PopDorNorT, ActDorNorT  !Daytime (2), night-time (1) or transition time (1-2)
    
  !!!Move to input file later!!! Rename SUEWS_AnthropogenicHeat.txt to SUEWS_Anthropogenic.txt and add these there?
  REAL(KIND(1d0)):: EnEF_v_Jkm
  
  ! Define coefficients ---------------------------------------   ! Move to inputs? !!!
  ! Energy emission factors
  EnEF_v_Jkm = 3.97e6  ! [J km−1] Sailor & Lu (2004)
  
  ! Establish whether weekday or weekend ---------------------- 
  iu=1     !Set to 1=weekday
  IF(DayofWeek(id,1)==1 .or. DayofWeek(id,1)==7) THEN  
     iu=2  !Set to 2=weekend
  ENDIF
   
  ! Calculate energy release from human metabolism -------------
  ! from Sailor & Lu (2004) 
  PopDorNorT = HumActivity_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu) !!! Set separately later !!!
  ActDorNorT = HumActivity_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
  ! Pop dens (cap ha-1 -> cap m-2) x activity level (W cap-1) 
  QF_metab = (PopDensNighttime*(2-PopDorNorT) + PopDensDaytime*(PopDorNorT-1))/10000 * (75*(2-ActDorNorT) + 175*(ActDorNorT-1))
  
  ! Calculate energy released from traffic ---------------------
  ! Use mean traffic rate [veh km cap-1 day-1] * emission factor [J km-1]
  ! Which popdens? !!!
  QF_traff = PopDensNighttime/10000 * TrafficRate/(60*60*24) * EnEF_v_Jkm * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
  
  ! Calculate CO2 emissions from building energy use ----------
  ! Use building energy use [W cap-1]
  QF_build = PopDensNighttime/10000 * BuildEnergyUse * AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)   !Need to work out how to incorporate temp dependence !!!
  !!! See Allen et al (2010) --> Ruth & Lin
  !!! QF_o = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1)))*numCapita
  
  ! Sum components to give QF_o
  QF_o = QF_metab + QF_traff + QF_build 
   
  !write(*,*) QF_o, QF_metab, QF_traff, QF_build
  
  RETURN
  
ENDSUBROUTINE SAHP_3
!-----------------------------------------------------------------------------------