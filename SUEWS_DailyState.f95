! Calculation of daily state variables
! Updates each timestep, but correct values calculated only at the end of each day
! Rest of the code uses value from previous day
! Responds to what has happened in the past (temperature, rainfall, etc)
! N.B. If changes are made here, may need to update code in SUEWS_Initial accordingly
!
! Last modified HCW 1 Jun 2015
! Bug fix 05 Jun now fixed in a different way - DecidCap is now treated the same as DecidAlb 
!  so should cope with multiple grids.
! Last modified HCW 05 Jun 2015
! Bug fix - set all current storage capacities (surf(6,)) to min. value, then set for DecTr
! Last modified LJ 11 Mar 2015
! Removed switch as no longer necessary
! Last modified HCW 06 Mar 2015
!  iy used instead of year which does not have a value here
! Last modified HCW 20 Feb 2015
!  Added surf(6,is) for the current storage capacity
!  Updated and corrected DailyState output file
! Last modified LJ 05 Feb 2015
!  DailyState saving fixed. Now header is printed and the file closed and opened as suggested.
! N.B. Bug in daily Precip - needs fixing!!! - HCW thinks this is fixed 20 Feb 2015
! Last modified HCW 26 Jan 2015
!  sfr and IrrFracs deleted from WU_Day calculations, so that WU_Day is not spread over
!   the total area
! Modified HCW 22 Jan 2015
! WU_Day now has 9 columns (EveTr, DecTr, Grass; automatic, manual, total) HCW 23/01/2015 
! Handles values for different grids (Gridiv & ir arguments) HCW 27/11/2014
! Added the calculation of surface temperature
! Snow albedo aging and calculation of snow density added, LJ 22/02/2013
! Calculation of LAI senescence from previous day length added, LJ 22/07/2013
! sg feb 2012 - rewritten from LUMPS_LAI so done in real time
!
! To Do - Check SnowUpdate following adjustment to 5-min timestep (LJ to do)
!	- Account for change of year in 5-day running mean?
!	- Check LAI calcs (N/S hemisphere similarities; use of day length)
!       - Take out doy limits (140,170, etc) and code as parameters
!	- Could add different coefficients (Ie_m, Ie_a) for each vegetation type
!==============================================================================
 subroutine DailyState(Gridiv)

  use allocateArray
  use data_in      
  use defaultNotUsed
  use Initial
  use snowMod
  use sues_data
  use time
  use VegPhenogy
      
  IMPLICIT NONE
        
  integer:: Gridiv
  integer:: gamma1,gamma2  !Switch for heating and cooling degree days
  integer:: iv,&           !Loop over vegetation types
  	        jj,&           !Loop over previous 5 days
            calc,&         !Water use calculation is done when calc = 1
  	        wd,&           !Number of weekday (Sun=1,...Sat=7)
            mb,&           !Months
            seas,&         !Season (summer=1, winter=2)
            date,&         !Day
            critDays       !Limit for GDD when GDD or SDD is set to zero

  real(kind(1d0)):: capChange,porChange,albChange,deltaLAI
  real(kind(1d0)):: no,yes,indHelp   !Switches and checks for GDD
  character(len=10):: grstr2
     
  ! --------------------------------------------------------------------------------
  ! ------------- Key to daily arrays ----------------------------------------------
  ! HDD(,1) ---- Heating  		 [degC] ! GDD(,1) ---- Growing  	[degC]
  ! HDD(,2) ---- Cooling   		 [degC] ! GDD(,2) ---- Senescence	[degC]   
  ! HDD(,3) ---- Daily mean temp	 [degC]	! GDD(,3) ---- Daily min temp	[degC]
  ! HDD(,4) ---- 5-day running mean temp [degC] ! GDD(,4) ---- Daily max temp	[degC]
  ! HDD(,5) ---- Daily precip total	 [mm]	! GDD(,5) ---- Daytime hours    [h]
  ! HDD(,6) ---- Days since rain   	 [d]
  !
  ! LAI(,1:3) -- LAI for each veg surface [m2 m-2]
  !
  ! WU_Day(,1) - Daily water use total for Irr EveTr (automatic+manual) [mm]
  ! WU_Day(,2) - Automatic irrigation for Irr EveTr 		  	[mm]
  ! WU_Day(,3) - Manual irrigation for Irr EveTr 			[mm]
  ! WU_Day(,4) - Daily water use total for Irr DecTr (automatic+manual) [mm]
  ! WU_Day(,5) - Automatic irrigation for Irr DecTr 		  	[mm]
  ! WU_Day(,6) - Manual irrigation for Irr DecTr 			[mm]
  ! WU_Day(,7) - Daily water use total for Irr Grass (automatic+manual) [mm]
  ! WU_Day(,8) - Automatic irrigation for Irr Grass 		    	[mm]
  ! WU_Day(,9) - Manual irrigation for Irr Grass 			[mm]
  ! --------------------------------------------------------------------------------
      
     
  !! Initialization ----------------------------------------------------------------- 
  !! These variables don't seem to be needed (commented out HCW 27 Nov 2014)
  !! If required, they will need updating for a non-hourly timestep
  !runT(it)=Temp_C      !runT has been initialized in SUEWS_initial to the previous day average
  !avT_h=sum(runT)/24   !Average daily temperature
  !runP(it)=Precip   !Same for precipitation
  !totP_h=sum(runP)     !Daily sum for precipitation
    
  ! Daily min and max temp (these get updated through the day) ---------------------
  GDD(id,3) = min(Temp_C,GDD(id,3))     !Daily min T in column 3
  GDD(id,4) = max(Temp_C,GDD(id,4))     !Daily max T in column 4
  if (avkdn>10) then
     GDD(id,5) = GDD(id,5)+1/nsh_real   !Cumulate daytime hours !Divide by nsh (HCW 01 Dec 2014)
  endif
      
  ! Calculations related to heating and cooling degree days (HDD) ------------------
  ! See Sailor & Vasireddy (2006) EMS Eq 1,2 (theirs is hourly timestep)
  if ((BaseTHDD-Temp_C)>=0) then   !Heating
     gamma1=1
  else
     gamma1=0 
  endif
  
  if ((Temp_C-BaseTHDD)>=0) then   !Cooling
     gamma2=1
  else
     gamma2=0    
  endif
    
  if(Gridiv == 1) tstepcount=tstepcount+1   !Add 1 to tstepcount only once for all grids  
       
  HDD(id,1)=HDD(id,1) + gamma1*(BaseTHDD-Temp_C)   !Heating
  HDD(id,2)=HDD(id,2) + gamma2*(Temp_C-BaseTHDD)   !Cooling
  HDD(id,3)=HDD(id,3) + Temp_C                     !Will become daily average temperature
  ! 	 4 ------------------------------------!   !5-day running mean  
  HDD(id,5)=HDD(id,5) + Precip                     !Daily precip total  
  ! 	 6 ------------------------------------!   !Days since rain

  ! Update snow density, albedo surface fraction
  if (snowUse==1) call SnowUpdate(Temp_C)

  ! ================================================================================
  ! This next part occurs only on the first or last timestep of each day
    
  ! On first timestep of each day, define whether the day each a workday or weekend
  if (it==0.and.imin==0) then

     call day2month(id,mb,date,seas,iy,lat)	!Calculate real date from doy
     call Day_of_Week(date,mb,iy,wd)   	    !Calculate weekday (1=Sun, ..., 7=Sat)

     dayofWeek(id,1)=wd      !Day of week
     dayofWeek(id,2)=mb      !Month
     dayofweek(id,3)=seas    !Season

  ! On last timestep, perform the daily calculations -------------------------------
  ! Daily values not correct until end of each day, so main program uses day before
  elseif (it==23.and.imin==(nsh_real-1)/nsh_real*60) then
     !write(*,*) 'Last timestep of day'    
     
     ! Heating degree days (HDD) -------------
     HDD(id,1)=HDD(id,1)/tstepcount   !Heating
     HDD(id,2)=HDD(id,2)/tstepcount   !Cooling
     HDD(id,3)=HDD(id,3)/tstepcount   !Average temp
        
     if(Gridiv == NumberOfGrids) tstepcount=0  !Set to zero only after last grid has run
     
     ! Calculate 5-day running mean temp
     !! Need to deal with the previous year
     do jj=1,5         
        HDD(id,4)=HDD(id,4) + HDD(id-(jj-1),3)
     enddo
     HDD(id,4) = HDD(id,4)/5
          
     ! Calculate days since rain     
     if(HDD(id,5)>0) then        !Rain occurred
        HDD(id,6)=0
     else
        HDD(id,6)=HDD(id-1,6)+1  !Days since rain
     endif
      
     
     ! Calculate modelled daily water use ------------------------------------------
     if (WU_choice==0) then   !If water use is to be modelled (rather than observed)
        wd=dayofWeek(id,1)
        
        if (DayWat(wd)==1.0) then      !1 indicates watering permitted on this day
           calc=0
           if (lat>=0) then            !Northern Hemisphere
              if (id>=Ie_start-1.and.id<=Ie_end+1) calc=1   !Day between irrigation period               
           else                        !Southern Hemisphere
              calc=1
              if (id>=Ie_end.and.id<=Ie_start) calc=0       !Day between irrigation period                       
           endif
         
           if(calc==1) then   
              ! Model daily water use based on HDD(id,6)(days since rain) and HDD(id,3)(average temp)
              ! WU_Day is the amount of water [mm] per day, applied to each of the irrigated areas
              ! N.B. These are the same for each vegetation type at the moment

              ! ---- Automatic irrigation (evergreen trees) ----
	          WU_day(id,2) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
	          if (WU_Day(id,2)<0) WU_Day(id,2)=0   !If modelled WU is negative -> 0

              ! ---- Manual irrigation (evergreen trees) ----
	          WU_day(id,3) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
	          if (WU_Day(id,3)<0) WU_Day(id,3)=0   !If modelled WU is negative -> 0

	          ! ---- Total evergreen trees water use (automatic + manual) ----
              WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))
                            
              ! ---- Automatic irrigation (deciduous trees) ----
	          WU_day(id,5) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
	          if (WU_Day(id,5)<0) WU_Day(id,5)=0   !If modelled WU is negative -> 0

              ! ---- Manual irrigation (deciduous trees) ----
	          WU_day(id,6) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
	          if (WU_Day(id,6)<0) WU_Day(id,6)=0   !If modelled WU is negative -> 0

              ! ---- Total deciduous trees water use (automatic + manual) ----
              WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))
                            
              ! ---- Automatic irrigation (grass) ----
              WU_day(id,8) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd) 
              if (WU_Day(id,8)<0) WU_Day(id,8)=0   !If modelled WU is negative -> 0

              ! ---- Manual irrigation (grass) ----
              WU_day(id,9) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
              if (WU_Day(id,9)<0) WU_Day(id,9)=0   !If modelled WU is negative -> 0

              ! ---- Total grass water use (automatic + manual) ----      
              WU_Day(id,7)=(WU_day(id,8)+WU_day(id,9))

           else   !If no irrigation on this day
              WU_Day(id,1)=0
              WU_Day(id,2)=0
              WU_Day(id,3)=0
              WU_Day(id,4)=0
              WU_Day(id,5)=0
              WU_Day(id,6)=0
              WU_Day(id,7)=0
	          WU_Day(id,8)=0
              WU_Day(id,9)=0
           endif
        endif
     endif 
   
     !------------------------------------------------------------------------------
     ! Calculation of LAI from growing degree days
     ! This was revised and checked on 16 Feb 2014 by LJ 
     !------------------------------------------------------------------------------

     critDays=50   !Critical limit for GDD when GDD or SDD is set to zero 
    
     ! Loop through vegetation types (iv)
     do iv=1,NVegSurf
        ! Calculate GDD for each day from the minimum and maximum air temperature
        yes =((GDD(id,3)+GDD(id,4))/2-BaseT(iv))    !Leaf on
        no  =((GDD(id,3)+GDD(id,4))/2-BaseTe(iv))   !Leaf off
      
        indHelp = 0   !Help switch to allow GDD to go to zero in sprint-time !!What does this mean?? HCW
          
        if(yes<0) then   !GDD cannot be negative 
           indHelp=yes   !Amount of negative GDD
           yes=0
        endif
            
        if(no>0) no=0    !SDD cannot be positive 
               
        ! Calculate cumulative growing and senescence degree days 
        GDD(id,1) = GDD(id-1,1)+yes
        GDD(id,2) = GDD(id-1,2)+no
       
        ! Possibility for cold spring
        if(GDD(id,2)<=SDDFull(iv).and.indHelp<0) then
           GDD(id,1)=0
        endif
            
        if(GDD(id,1)>=GDDFull(iv)) then   !Start senescence
           GDD(id,1)=GDDFull(iv)          !Leaves should not grow so delete yes from earlier
           if(GDD(id,2)<-critDays) GDD(id,1)=0
        endif
          
        if (GDD(id,2)<=SDDFull(iv)) then   !After senescence now start growing leaves      
           GDD(id,2)=SDDFull(iv)           !Leaves off so add back earlier 
           if(GDD(id,1)>critDays) GDD(id,2)=0
        endif
          
        ! With these limits SDD, GDD is set to zero
        if(GDD(id,2)<-critDays.and.GDD(id,2)>SDDFull(iv))  GDD(id,1)=0    
        if(GDD(id,1)> critDays.and.GDD(id,1)<GDDFull(iv))  GDD(id,2)=0 
          
        ! Now calculate LAI itself
        if(lat>=0) then   !Northern hemispere
           if (id==140.and.GDD(id,2)/=0)  GDD(id,2)=0  !If SDD is not zero by mid May, this is forced
           ! Set SDD to zero in summer time
           if (GDD(id,1)> critDays.and.id<170) GDD(id,2)=0
           ! Set GDD zero in winter time
           if (GDD(id,2)<-critDays.and.id>170) GDD(id,1)=0
      
           if (LAItype < 0.5) then   !Original LAI type
              if(GDD(id,1)>0.and.GDD(id,1)<GDDFull(iv)) then       !Leaves can still grow
                 lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
              elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then   !Start senescence
                 lai(id,iv)=(lai(id-1,iv)**laiPower(3)*GDD(id,2)*laiPower(4))+lai(id-1,iv) 
              else
                 lai(id,iv)=lai(id-1,iv)
              endif
           elseif (LAItype>=0.5) then            
              if(GDD(id,1)>0.and.GDD(id,1)<GDDFull(iv)) then        !Leaves can still grow
                 lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
              !! Use day length to start senescence at high latitudes (N hemisphere)  
              elseif (GDD(id,5)<=12.and.GDD(id,2)>SDDFull(iv)) then !Start senescence	
                 lai(id,iv)=(lai(id-1,iv)*laiPower(3)*(1-GDD(id,2))*laiPower(4))+lai(id-1,iv)
              else
                 lai(id,iv)=lai(id-1,iv)
              endif
           endif
            
        elseif (lat<0) then   !Southern hemisphere !! N.B. not identical to N hemisphere - return to later
           if (id==300.and.GDD(id,2)/=0)  GDD(id,2)=0   !If SDD is not zero by late Oct, this is forced
           ! Set SDD to zero in summer time
           if (GDD(id,1)> critDays.and.id>250) GDD(id,2)=0
           ! Set GDD zero in winter time
           if (GDD(id,2)<-critDays.and.id<250) GDD(id,1)=0
        
           if (LAItype < 0.5) then   !Original LAI type
              if(GDD(id,1)>0.and.GDD(id,1)<GDDFull(iv)) then
                 lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
              elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then
                 lai(id,iv)=(lai(id-1,iv)**laiPower(3)*GDD(id,2)*laiPower(4))+lai(id-1,iv) 
              else
                 lai(id,iv)=lai(id-1,iv)
              endif
           else  
              if(GDD(id,1)>0.and.GDD(id,1)<GDDFull(iv)) then
                 lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
              !! Day length not used to start senescence in S hemisphere (not much land) 
              elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then
                 lai(id,iv)=(lai(id-1,iv)*laiPower(3)*(1-GDD(id,2))*laiPower(4))+lai(id-1,iv) 
              else
                 lai(id,iv)=lai(id-1,iv)
              endif
           endif
        endif   !N or S hemisphere
           
        ! Check LAI within limits; if not set to limiting value   
        if(lai(id,iv).gt.LAImax(iv))then
           lai(id,iv)=LAImax(iv)
        elseif(lai(id,iv).lt.LAImin(iv))then
           lai(id,iv)=LAImin(iv)
        endif    
     
     enddo   !End of loop over veg surfaces
     !------------------------------------------------------------------------------
       
     !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     ! Calculate the development of deciduous cover, albedo and porosity
     ! If only LUMPS is used, set deciduous capacities to 0
                   
     !assume porosity Change based on GO99- Heisler?
     ! max-min for 
     ! sg changed to increase or decrease appropriately for change in LAI
     iv=ivDecid
     deltaLAI=0
     CapChange=0
     porChange=0
     albChange=0
            
     if((LAI(ID,iv)-LAI(ID-1,iv))/=0) then
        deltaLAI=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LaiMin(iv))
        CapChange=(CapMin_dec-CapMax_dec)* deltaLAI
        albChange=(alBMax_dec-AlbMin_dec)* deltaLAI          
        porChange=(0.2-0.6)* deltaLAI     
     endif        
     DecidCap(id)=DecidCap(id-1)-CapChange
     porosity(id)=porosity(id-1)-porChange
     albDec(id)=albDec(id-1)+albChange
      
     ! -----------------------------------------------------------------------------
     
     !--------------------------------------------------------------------------
     !Write out DailyState file (1 row per day)

     if (writedailyState==1) then
        !Define filename
        write(grstr2,'(i2)') Gridiv      !Convert grid number for output file name
        FileDaily=trim(FileOutputPath)//trim(FileCode)//trim(adjustl(grstr2))//'_DailyState.txt'

        ! If first modelled day, open the file and save header
        if (DailyStateFirstOpen(Gridiv)==1) then
           open(60,file=FileDaily)
           write(60,142)
           142  format('%year id ',&                                          !2
                       'HDD1_h HDD2_c HDD3_Tmean HDD4_T5d P/day DaysSR ',&    !8
                       'GDD1_g GDD2_s GDD3_Tmin GDD4_Tmax GDD5_DayLHrs ',&    !13
                       'LAI_EveTr LAI_DecTr LAI_Grass ',&                     !16
                       'DecidCap Porosity AlbDec ',&                          !19
                       'WU_EveTr(1) WU_EveTr(2) WU_EveTr(3) ',&               !22
                       'WU_DecTr(1) WU_DecTr(2) WU_DecTr(3) ',&               !25
                       'WU_Grass(1) WU_Grass(2) WU_Grass(3) ',&               !28
                       'deltaLAI LAIlumps AlbSnow dens_snow_pav ',&           !32
                       'dens_snow_bldg dens_snow_EveTr dens_snow_DecTr',&     !35
                       'dens_snow_Grass dens_snow_Bares dens_snow_wtr')       !38
                       
           DailyStateFirstOpen(Gridiv)=0
        ! Otherwise open file to append
        else
           open(60,file=FileDaily,position='append')
        endif      
        
        ! Write actual data     
        write(60,601) iy,id,&
                      HDD(id,1:6),GDD(id,1:5),&
                      LAI(id,1:nvegsurf),&
                      DecidCap(id),Porosity(id),AlbDec(id),&
                      WU_day(id-1,1:9),&
                      deltaLAI,VegPhenLumps,alb_snow,densSnow(1:7)
            
        601 format(2(i4,1X),&
                   6(f6.1,1X), 5(f6.1,1X),&
                   3(f6.2,1X),&
                   3(f6.2,1X),&
                   9(f7.3,1X),&
                   2(f7.2,1X),8(f7.2,1X))
         
        ! Close the daily state file
        close(60)
     endif  !End writedailystate
     
     ! Set Daylight Saving ---------------------------------------------------------   
     if (id>=DayLightSavingDay(1).and.id<=DayLightSavingDay(2)) then   !Summer time
        DLS=1
     else
        DLS=0
     endif
     ! -----------------------------------------------------------------------------
         
  endif   !End of section done only at the end of each day (i.e. only once per day)
  ! ================================================================================    
                                
  return
endsubroutine DailyState
!=================================================================================== 