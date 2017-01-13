
!===================================================================================
!Simple Anthropogenic Heat Parameterization routines
!Last modified
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
 subroutine SAHP_1_v2015(QF_o,QF_o_base,QF_o_heat,id,ih,imin)
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
     
 endsubroutine SAHP_1_v2015
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
subroutine SAHP_2_v2015(QF_o,QF_o_base,QF_o_heat,id,ih,imin)
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
!write(*,*) (NSH*(ih+1-1)+imin*NSH/60+1)
!write(*,*) AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)
!write(*,*) (Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1))*numCapita
   
  QF_o = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1)))*numCapita

  QF_o_base = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_a(iu)))*numCapita
  QF_o_heat = (AHProf_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*(Qf_c(iu)*HDD(id-1,1)))*numCapita

  !write(*,*) QF_o
!write(*,*) '----------------------'
!write(*,*) " " 

endsubroutine SAHP_2_v2015
!-----------------------------------------------------------------------------------


