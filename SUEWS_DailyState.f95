!sg feb 2012
!only need to calculate some things each day - but they need to respond to what has happened in the past
!Added the calculation of surface temperature
!Snow albedo aging and calculation of snow density added, LJ 22/02/2013
!Calculation of LAI senescence from previous day length added, LJ 22/7/2013

Subroutine DailyState
    !==============================================================================
    !  Updates each timestep
    !sg feb 2012 - rewritten from LUMPS_LAI so done in real time
    !==============================================================================
    use sues_data
    use allocateArray
    use time
    use data_in        ! sdec1,sdec2,sdec3,sdec4
    use VegPhenogy
    use snowMod
    use defaultNotUsed
    
    implicit none
    
    integer::gamma1,gamma2,iv,jj,wd,seas,date,mb,switch=0,j,critDays,calc, seasonLAI
    real (kind(1d0)):: no,yes,CapChange,porChange,albChange,deltaLAI, indHelp
    
    !-------------------------------------------------------------------------------------
    
    !Initialization 
    runT(it)=Temp_C    !runT has been initialized in SUEWS_initial to the previous day average
    avT_h=sum(runT)/24 !Average daily temperature
    runP(it)=Precip_hr !Same for recipitation
    totP_h=sum(runP)   !Daily sum for precipitation
    
    !GDD-1 growing  2-- senesecence
    ! changes through the day 
    GDD(id,3)=min(Temp_C,GDD(id,3))     ! daily min T in column 3
    GDD(id,4)=max(Temp_C,GDD(id,4))     ! daily max T in column 4
    if (avkdn>10)  GDD(id,5)=GDD(id,5)+1 ! Cumulate daytime hours
    

    !Calculations related to heating and cooling degree days (HDD)
    
    if ((Temp_C-BaseTHDD)>=0) then !Cooling
       gamma2=1
    else
        gamma2=0    
    endif
    
    if ((BaseTHDD-Temp_C)>=0) then !Heating
       gamma1=1
    else
       gamma1=0 
    endif
    
    hrcount=hrcount+1         
        
    HDD(id,1)=HDD(id,1) + gamma1*(BaseTHDD-Temp_C) ! Heating
    HDD(id,2)=HDD(id,2) + gamma2*(Temp_C-BaseTHDD) ! Cooling
    HDD(id,3)=HDD(id,3) + Temp_C                   ! will become average 
    ! 4--------------------------------------------! 5 day running mean  
    HDD(id,5)=HDD(id,5) + Precip_hr                ! daily precip total
    ! 6 -------------------------------------------! days since rain
    

    !Calculation of snow albedo and density from previous timestep (delta t = 1 h).
    call SnowUpdate

    !=============================================================================
    !Occurs only on first or last hours of the day
    if(it==FirstTimeofDay) then
      if(id<=1)then
          year=year-1
          call LeapYearCalc (year,id)
          switch=1
          call ErrorHint(43,'switch- to last day of last year',notUsed,notUsed,notUsedI)
          
      endif 
  
      call day2month(id,mb,date,seas,year,lat)!Calculate real date from doy
      call Day_of_Week(date,mb,year,wd)!Calculate weekday (1=Sun,...)

      if(switch==1)then
        year=year+1
        id=1
        switch=0
      endif
     
      dayofWeek(id,1)=wd      ! day of week
      dayofWeek(id,2)=mb      ! month
      dayofweek(id,3)=seas    ! season



  
  !---------------Daily calculation -------------------
  ! use the day before in main programme
  !----------------------------------------------------
  elseif(it==LastTimeOfDay) then
 
  !---------------------------------------------------
      !Heating degree days (HDD)
  !---------------------------------------------------
      HDD(id,1)=HDD(id,1)/hrcount  !Heating
      HDD(id,2)=HDD(id,2)/hrcount  !Cooling
      HDD(id,3)=HDD(id,3)/hrcount  !Average
      
      !write(12,*)hdd(id,3),hrcount,id,it,year
      hrcount=0
      do jj=1,5  ! 5 day running mean    
        ! need to deal with the previous year - to do
          HDD(id,4)=HDD(id,4) + HDD(id-(jj-1),3)
      enddo
      HDD(id,4) = HDD(id,4)/5 ! 5 day mean
      if(HDD(id,5)>0) then  ! rain occurred
         HDD(id,6)=0
      else
         HDD(id,6)=HDD(id-1,6)+1     ! days since rain
      endif
  !------------------------------------------------------
  !       Calculate daily water use if this is modelled
  !------------------------------------------------------
    
      if (WU_choice==0) then
        
        wd=dayofWeek(id,1)
        if (DayWat(wd)==1.0) then      !if=1 - then this is a day that has watering
           calc=0
            if (lat>=0)then             !Northern Hemisphere
                if (id>=Ie_start-1.and.id<=Ie_end+1) calc=1 !if day between irrigation period               
            else                        !Southern Hemisphere
                calc=1
                if (id>=Ie_end.and.id<=Ie_start) calc=0
                       !if day between irrigation period                       
            endif
            
            if(calc==1) then             
                !
                ! HDD(id,6) -- daysSincerain, HDD(id3)- mean airT
                ! automatic
                WU_day(id,2)=Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*sfr(GrassISurf)*DayWatPer(wd) 
                ! manual
                WU_day(id,3)=(1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*sfr(GrassISurf)*DayWatPer(wd)
                           
            
                if (WU_Day(id,2)<0) WU_Day(id,2)=0 !If modelled WU is negative -> 0
                if (WU_Day(id,3)<0) WU_Day(id,3)=0 !If modelled WU is negative -> 0
                                    
                WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))

				!Calculate the fraction for irrigated trees/shrubs
                !Added by LJ in 2 September 2013
                IrrTrees = sfr(ConifSurf)*IrrFractionTrees+sfr(DecidSurf)*IrrFractionTrees

                WU_day(id,5)=Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*IrrTrees*DayWatPer(wd) 
                ! manual
                WU_day(id,6)=(1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*IrrTrees*DayWatPer(wd)
                           
            
                if (WU_Day(id,5)<0) WU_Day(id,5)=0 !If modelled WU is negative -> 0
                if (WU_Day(id,6)<0) WU_Day(id,6)=0 !If modelled WU is negative -> 0
                                    
                WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))
              
                
              else
                 WU_Day(id,1)=0
                 WU_Day(id,2)=0
                 WU_Day(id,3)=0
                 WU_Day(id,4)=0
                 WU_Day(id,5)=0
                 WU_Day(id,6)=0
             endif
          endif
       endif
   
  !------------------------------------------------------
  !    Calculation of LAI from growing degree days
  !    This was revised and checked on 16 Feb 2014 by LJ 
  !------------------------------------------------------
  
       critDays=50 !Critical limit for GDD when GDD or SDD is set to zero 
       
       !Go vegetation type (iv) through
       do iv=1,NVegSurf

          !Calculate GDD for each day from the minimum and maximum air temperature.
          yes=((GDD(id,3)+GDD(id,4))/2-BaseT(iv))   ! leaf on
          no=((GDD(id,3)+GDD(id,4))/2-BaseTe(iv))   ! leaf off
         
          indHelp = 0   !Help switch to allow GDD fo to zero in sprint-time
          
          if(yes<0) then !GDD cannot be negative 
            indHelp=yes  !Amount of negative GDD
            yes=0
          endif
            
          if(no>0) no=0 !SDD cannot be positive 
            
          !Calculate cumulative growing and senescence degree days 
          GDD(id,1)=GDD(id-1,1)+yes
          GDD(id,2)=GDD(id-1,2)+no
       
          !Possibility for cold spring
          if (GDD(id,2)<=SDDFull(iv).and.indHelp<0) then
             GDD(id,1)=0
          endif
            
          
          if(GDD(id,1)>=GDDFull(iv))then ! start senescence
              GDD(id,1)=GDDFull(iv)       ! leaves should not grow so delete yes from earlier
              if(GDD(id,2)<-critDays)GDD(id,1)=0
          endif
          
          if (GDD(id,2)<=SDDFull(iv))then  ! after senescence now start growing leaves      
               GDD(id,2)=SDDFull(iv)    ! leaves off so add back earlier 
               if(GDD(id,1)>critDays)GDD(id,2)=0
          endif
          
          !With these limits SDD of GDD is set to zero
          if(GDD(id,2)<-critDays.and.GDD(id,2)>SDDFull(iv))  GDD(id,1)=0    
          if(GDD(id,1)>critDays.and.GDD(id,1)<GDDFull(iv))  GDD(id,2)=0 
          
          !Now the LAI itself

          if(lat>=0) then  !Northern hemispere
             
             if (id==140.and.GDD(id,2)/=0)  GDD(id,2)=0 !If SDD (GDD(id,1)) is not zero by early may,this is forced
       
             !Set SDD to zero in summer time
             if (GDD(id,1)>critDays.and.id<170) GDD(id,2)=0
             !Set GDD zero in winter time
             if (GDD(id,2)<-critDays.and.id>170) GDD(id,1)=0
      
             
             
             if (LAItype < 0.5) then !Original LAI type
                if(gdd(id,1)>0.and.gdd(id,1)<gddfull(iv)) then      ! then leaves can still grow
                  lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
                elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then  ! start senescence
                   lai(id,iv)=(lai(id-1,iv)**laiPower(3)*GDD(id,2)*laiPower(4))+lai(id-1,iv) 
                else
                  lai(id,iv)=lai(id-1,iv)
                endif
             elseif (LAItype>=0.5) then
                if(gdd(id,1)>0.and.gdd(id,1)<gddfull(iv)) then      ! then leaves can still grow
                    lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
                  elseif (GDD(id,5)<=12.and.GDD(id,2)>SDDFull(iv)) then  ! start senescence
                    lai(id,iv)=(lai(id-1,iv)*laiPower(3)*(1-GDD(id,2))*laiPower(4))+lai(id-1,iv)
                  else
                    lai(id,iv)=lai(id-1,iv)
                  endif
              
              endif
            
          elseif (lat<0) then     !Southern hemisphere 

             if (id==300.and.GDD(id,2)/=0)  GDD(id,2)=0 !If SDD (GDD(id,1)) is not zero by early may,this is forced
       
             !Set SDD to zero in summer time
             if (GDD(id,1)>critDays.and.id>250) GDD(id,2)=0
             !Set GDD zero in winter time
             if (GDD(id,2)<-critDays.and.id<250) GDD(id,1)=0
      
             

             if (LAItype < 0.5) then !Original LAI type
               
                if(gdd(id,1)>0.and.gdd(id,1)<gddfull(iv)) then
                      lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
                elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then
                      lai(id,iv)=(lai(id-1,iv)**laiPower(3)*GDD(id,2)*laiPower(4))+lai(id-1,iv) 
                else
                      lai(id,iv)=lai(id-1,iv)
                endif

             else
                if(gdd(id,1)>0.and.gdd(id,1)<gddfull(iv)) then
                      lai(id,iv)=(lai(id-1,iv)**laiPower(1)*GDD(id,1)*laiPower(2))+lai(id-1,iv)
                elseif(GDD(id,2)<0.and.GDD(id,2)>SDDFull(iv)) then
                      lai(id,iv)=(lai(id-1,iv)*laiPower(3)*(1-GDD(id,2))*laiPower(4))+lai(id-1,iv) 
                else
                      lai(id,iv)=lai(id-1,iv)
                endif

             endif
        endif                      
                                 
        if(lai(id,iv).gt.LAImax(iv))then
           lai(id,iv)=LAImax(iv)
        elseif(lai(id,iv).lt.LAImin(iv))then
           lai(id,iv)=LAImin(iv)
        endif    
   
       enddo

       
       !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       !Calculate the development of deciduous cover, albedo and porosity
       !If only LUMPS is used, set deciduous capacities to 0
              
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
           !CapChange=(surf(1,DecidSurf)-surf(5,DecidSurf))* deltaLAI
           CapChange=(CapMin_dec-CapMax_dec)* deltaLAI
           albChange=(alBMax_dec-AlbMin_dec)* deltaLAI          
           porChange=(0.2-0.6)* deltaLAI     
       endif        
       DecidCap(id)=DecidCap(id-1)-CapChange
       porosity(id)=porosity(id-1)-porChange
       albDec(id)=albDec(id-1)+albChange
       SURF(1,DecidSurf)=DecidCap(id)  !Surface moisture capacity of  deciduous trees        
        
!---------------------------------------------------------------------------------------
       !Write to the dailystate file
       write(60,601)int(year),id,(HDD(id,j),j=1,6),(GDD(id,j),j=1,5),(lai(id,iv),iv=1,nvegsurf),&
                    decidcap(id),porosity(id),albdec(id), WU_day(id-1,1),WU_day(id-1,2),&
                    Wu_day(id-1,3),WU_day(id-1,4),WU_day(id-1,5),Wu_day(id-1,6),deltaLAI,&
                    VegPhenLumps,alb_snow,(densSnow(j),j=1,6)
       601 format(2(i4,1X),6(f6.1,1X),2(f5.0,1X),3(f6.1,1X),4(f6.2,1X),3(f6.2,1X),6(f7.3,1X),&
       17(f7.2,1X))  
       
       if (id>=DayLightSavingDay(1).and.id<=DayLightSavingDay(2)) then!Summer time
           DLS=1
       else
           DLS=0
       endif
   !  pause 'end of day'
    
    endif   ! just done once a day
 !write(12,'(2i4,4f6.1)')id,it,avT_H,hdd(id-1,3),totP_h,hdd(id,5)

  return
  end subroutine DailyState
  