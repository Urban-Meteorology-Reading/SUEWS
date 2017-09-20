! Calculation of daily state variables
! Responds to what has happened in the past (temperature, rainfall, etc)
! Updates each time step, but for many variables, correct values are calculated only at the end of each day!
! --> for these variables, the rest of the code MUST use values from the previous day
! N.B. Some of this code is repeated in SUEWS_Initial
! --> so if changes are made here, SUEWS_Initial may also need to be updated accordingly
! N.B. Currently, daily variables are calculated using 00:00-23:55 timestamps (for 5-min resolution); should use 00:05-00:00
!
! Last modified:
!  TS 18 Sep 2017  - Added explicit interface
!  TS 07 Jun 2017  - Improve the format of output with more friendly alignment
!  HCW 04 Jul 2016 - GridID can now be up to 10 digits long
!  HCW 25 May 2016 - Added extra columns to daily state file (albedo for EveTr and Grass)
!  HCW 24 May 2016 - Bug fixed in naming of DailyState file (now uses GridIDmatrix(Gridiv) rather than Gridiv)
!  LJ 27 Jan 2016  - Removal of tabs
!  HCW 20 Aug 2015 - Sign of the porosity change corrected so that porosity is greatest when LAI is smallest
!  HCW 03 Jul 2015 - Increased output resolution of P/day in DailyState file to avoid rounding errors.
!                    Albedo of EveTr and Grass now adjusted based on change in LAI for EveTr and Grass
!                    (rather than DecTr)
!  HCW 29 Jun 2015 - Added albChange for EveTr and Grass surfaces
!  HCW 11 Jun 2015 - Bug fix from 05 Jun now fixed in a different way -
!                    DecidCap is now treated the same as DecidAlb so should be able to cope with multiple grids.
!  HCW 05 Jun 2015 - Bug fix - set all current storage capacities (surf(6,)) to min. value, then set for DecTr
!  LJ 11 Mar 2015  - Removed switch as no longer necessary
!  HCW 06 Mar 2015 - iy used instead of year which does not have a value here
!  HCW 20 Feb 2015 - Added surf(6,is) for the current storage capacity
!  Updated and corrected DailyState output file
!  LJ 05 Feb 2015  - DailyState saving fixed. Now header is printed and the file closed and opened as suggested.
! N.B. Bug in daily Precip - needs fixing!!! - HCW thinks this is fixed 20 Feb 2015
!  HCW 26 Jan 2015 - sfr and IrrFracs deleted from WU_Day calculations, so that WU_Day is not spread over
!  the total area
!  HCW 23 Jan 2015 - WU_Day now has 9 columns (EveTr, DecTr, Grass; automatic, manual, total)
!  HCW 27 Nov 2014 - Handles values for different grids (Gridiv & ir arguments)
! Added the calculation of surface temperature
!  LJ 22 Feb 2013  - Snow albedo aging and calculation of snow density added,
!  LJ 22 Jul 2013  - Calculation of LAI senescence from previous day length added
! sg feb 2012 - rewritten from LUMPS_LAI so done in real time
!
! To Do
!   - Account for change of year in 5-day running mean?
!   - Check LAI calcs (N/S hemisphere similarities; use of day length)
!       - Take out doy limits (140,170, etc) and code as parameters
!   - Could add different coefficients (Ie_m, Ie_a) for each vegetation type
!==============================================================================
SUBROUTINE DailyState(&
     iy,id,it,imin,Gridiv,tstep,&!input
     WaterUseMethod,snowUse,Ie_start,Ie_end,&
     ReadLinesMetdata,&
     ncolumnsDataOut,&
     NumberOfGrids,&
     LAICalcYes,&
     LAIType,&
     GridIDmatrix,&
     nsh_real,&
     avkdn,&
     Temp_C,&
     Precip,&
     BaseTHDD,&
     lat,&
     Faut,&
     LAI_obs,&
     tau_a,&
     tau_f,&
     tau_r,&
     SnowDensMax,&
     SnowDensMin,&
     SnowAlbMin,&
     alBMax_DecTr,&
     alBMax_EveTr,&
     alBMax_Grass,&
     AlbMin_DecTr,&
     AlbMin_EveTr,&
     AlbMin_Grass,&
     CapMax_dec,&
     CapMin_dec,&
     PorMax_dec,&
     PorMin_dec,&
     VegPhenLumps,&
     Ie_a,&
     Ie_m,&
     DayWatPer,&
     DayWat,&
     SnowPack,&
     BaseT,&
     BaseTe,&
     GDDFull,&
     SDDFull,&
     LAIMin,&
     LAIMax,&
     LAIPower,&
     dataOut ,&
     FileCode,&
     FileOutputPath,&

     a1,& !inout
     a2,&
     a3,&
     tstepcount,&
     SnowAlb,&
     DecidCap,&
     albDecTr,&
     albEveTr,&
     albGrass,&
     porosity,&
     GDD,&
     HDD,&
     SnowDens,&
     LAI,&
     DailyStateFirstOpen,&

     DayofWeek,&!output
     WU_Day)

  ! USE allocateArray
  ! USE data_in
  ! USE defaultNotUsed
  ! USE Initial
  ! USE snowMod
  ! USE sues_data
  ! USE time
  ! USE VegPhenogy
  ! USE InitialCond

  IMPLICIT NONE
  INTEGER, PARAMETER:: ndays = 366   !Max no. days in a year used to specify size of daily arrays
  INTEGER, PARAMETER:: nsurf=7                !Total number of surfaces
  INTEGER, PARAMETER:: NVegSurf=3             !Number of surfaces that are vegetated
  ! INTEGER, PARAMETER:: nsurfIncSnow=nsurf+1   !Number of surfaces + snow
  INTEGER, PARAMETER:: MaxNumberOfGrids=2000   !Max no. grids   !HCW changed to 2000 from 10000 so prog can run on windows (2GB lim)
  ! INTEGER, PARAMETER:: MaxLinesMet=8640        !Max no. lines to read in one go (for all grids, ie MaxLinesMet/NumberOfGrids each)


  INTEGER,INTENT(in)::iy
  INTEGER,INTENT(in)::id
  INTEGER,INTENT(in)::it
  INTEGER,INTENT(in)::imin

  INTEGER,INTENT(in)::Gridiv
  INTEGER,INTENT(in)::tstep
  INTEGER,INTENT(in)::WaterUseMethod
  INTEGER,INTENT(in)::snowUse
  INTEGER,INTENT(in)::Ie_start   !Starting time of water use (DOY)
  INTEGER,INTENT(in)::Ie_end       !Ending time of water use (DOY)
  INTEGER,INTENT(in)::ReadLinesMetdata
  INTEGER,INTENT(in)::ncolumnsDataOut
  INTEGER,INTENT(in)::NumberOfGrids
  INTEGER,INTENT(in)::LAICalcYes


  INTEGER,DIMENSION(nvegsurf),INTENT(in):: LAIType                  !LAI equation to use: original (0) or new (1)

  INTEGER,DIMENSION(MaxNumberOfGrids),INTENT(in):: GridIDmatrix         !Array containing GridIDs in SiteSelect after sorting

  REAL(KIND(1d0)),INTENT(in)::nsh_real
  REAL(KIND(1d0)),INTENT(in)::avkdn
  REAL(KIND(1d0)),INTENT(in)::Temp_C
  REAL(KIND(1d0)),INTENT(in)::Precip
  REAL(KIND(1d0)),INTENT(in)::BaseTHDD
  REAL(KIND(1d0)),INTENT(in)::lat
  REAL(KIND(1d0)),INTENT(in)::Faut
  REAL(KIND(1d0)),INTENT(in)::LAI_obs
  REAL(KIND(1D0)),INTENT(in)::tau_a
  REAL(KIND(1D0)),INTENT(in)::tau_f
  REAL(KIND(1D0)),INTENT(in)::tau_r
  REAL(KIND(1D0)),INTENT(in)::SnowDensMax
  REAL(KIND(1D0)),INTENT(in)::SnowDensMin
  REAL(KIND(1D0)),INTENT(in)::SnowAlbMin
  REAL(KIND(1d0)),INTENT(in)::alBMax_DecTr
  REAL(KIND(1d0)),INTENT(in)::alBMax_EveTr
  REAL(KIND(1d0)),INTENT(in)::alBMax_Grass
  REAL(KIND(1d0)),INTENT(in)::AlbMin_DecTr
  REAL(KIND(1d0)),INTENT(in)::AlbMin_EveTr
  REAL(KIND(1d0)),INTENT(in)::AlbMin_Grass
  REAL(KIND(1d0)),INTENT(in)::CapMax_dec
  REAL(KIND(1d0)),INTENT(in)::CapMin_dec
  REAL(KIND(1d0)),INTENT(in)::PorMax_dec
  REAL(KIND(1d0)),INTENT(in)::PorMin_dec
  REAL(KIND(1d0)),INTENT(in) ::VegPhenLumps

  REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::Ie_a
  REAL(KIND(1d0)),DIMENSION(3),INTENT(in) ::Ie_m !Coefficients for automatic and manual irrigation models
  REAL(KIND(1d0)),DIMENSION(7),INTENT(in) ::DayWatPer !% of houses following daily water
  REAL(KIND(1d0)),DIMENSION(7),INTENT(in) ::DayWat !Days of watering allowed


  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(in)      ::SnowPack
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::BaseT !Base temperature for growing degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::BaseTe !Base temperature for senescence degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::GDDFull !Growing degree days needed for full capacity [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::SDDFull !Senescence degree days needed to initiate leaf off [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::LAIMin !Min LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)   ::LAIMax !Max LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(4,nvegsurf),INTENT(in) ::LAIPower !Coeffs for LAI equation: 1,2 - leaf growth; 3,4 - leaf off


  REAL(KIND(1d0)),DIMENSION(ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids),INTENT(in)::dataOut

  CHARACTER (len = 20),INTENT(in) :: FileCode       !Set in RunControl
  CHARACTER (len = 150),INTENT(in):: FileOutputPath !Filepath for output files (set in RunControl)

  REAL(KIND(1d0)),INTENT(inout)::a1
  REAL(KIND(1d0)),INTENT(inout)::a2
  REAL(KIND(1d0)),INTENT(inout)::a3
  REAL(KIND(1d0)),INTENT(inout)::tstepcount
  REAL(KIND(1d0)),INTENT(inout)::SnowAlb

  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::DecidCap
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albDecTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albEveTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albGrass
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::porosity
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5),INTENT(inout):: GDD !Growing Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, 6),INTENT(inout):: HDD          !Heating Degree Days (see SUEWS_DailyState.f95)

  REAL(KIND(1d0)),DIMENSION(nsurf),INTENT(inout)::SnowDens
  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf),INTENT(inout):: LAI !LAI for each veg surface [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids),INTENT(inout) ::DailyStateFirstOpen


  INTEGER,DIMENSION(0:ndays,3),INTENT(out)::DayofWeek

  REAL(KIND(1d0)),DIMENSION(0:ndays,9),INTENT(out):: WU_Day       !Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)

  INTEGER::date
  REAL(KIND(1d0))::deltaLAI

  REAL(KIND(1d0))::xBo
  REAL(KIND(1d0))::xmAH
  REAL(KIND(1d0)),DIMENSION(44) ::DailyStateLine
  ! INTEGER::gamma1
  ! INTEGER::gamma2

  ! --------------------------------------------------------------------------------
  ! ------------- Key to daily arrays ----------------------------------------------
  ! HDD(,1) ---- Heating         [degC] ! GDD(,1) ---- Growing      [degC]
  ! HDD(,2) ---- Cooling         [degC] ! GDD(,2) ---- Senescence   [degC]
  ! HDD(,3) ---- Daily mean temp     [degC] ! GDD(,3) ---- Daily min temp   [degC]
  ! HDD(,4) ---- 5-day running mean temp [degC] ! GDD(,4) ---- Daily max temp   [degC]
  ! HDD(,5) ---- Daily precip total  [mm]   ! GDD(,5) ---- Daytime hours    [h]
  ! HDD(,6) ---- Days since rain     [d]
  !
  ! LAI(,1:3) -- LAI for each veg surface [m2 m-2]
  !
  ! WU_Day(,1) - Daily water use total for Irr EveTr (automatic+manual) [mm]
  ! WU_Day(,2) - Automatic irrigation for Irr EveTr             [mm]
  ! WU_Day(,3) - Manual irrigation for Irr EveTr            [mm]
  ! WU_Day(,4) - Daily water use total for Irr DecTr (automatic+manual) [mm]
  ! WU_Day(,5) - Automatic irrigation for Irr DecTr             [mm]
  ! WU_Day(,6) - Manual irrigation for Irr DecTr            [mm]
  ! WU_Day(,7) - Daily water use total for Irr Grass (automatic+manual) [mm]
  ! WU_Day(,8) - Automatic irrigation for Irr Grass                 [mm]
  ! WU_Day(,9) - Manual irrigation for Irr Grass            [mm]
  ! --------------------------------------------------------------------------------

  CALL init_DailyState(&
       id,& !input
       Gridiv,&
       avkdn,&
       Temp_C,&
       Precip,&
       BaseTHDD,&
       nsh_real,&
       tstepcount,& !inout
       GDD,&
       HDD)
  ! !! Initialization -----------------------------------------------------------------
  ! !! These variables don't seem to be needed (commented out HCW 27 Nov 2014)
  ! !! If required, they will need updating for a non-hourly timestep
  ! !runT(it)=Temp_C      !runT has been initialized in SUEWS_initial to the previous day average
  ! !avT_h=sum(runT)/24   !Average daily temperature
  ! !runP(it)=Precip   !Same for precipitation
  ! !totP_h=sum(runP)     !Daily sum for precipitation
  !
  ! ! Daily min and max temp (these get updated through the day) ---------------------
  ! GDD(id,3) = MIN(Temp_C,GDD(id,3))     !Daily min T in column 3
  ! GDD(id,4) = MAX(Temp_C,GDD(id,4))     !Daily max T in column 4
  ! IF (avkdn>10) THEN
  !    GDD(id,5) = GDD(id,5)+1/nsh_real   !Cumulate daytime hours !Divide by nsh (HCW 01 Dec 2014)
  ! ENDIF
  !
  ! ! Calculations related to heating and cooling degree days (HDD) ------------------
  ! ! See Sailor & Vasireddy (2006) EMS Eq 1,2 (theirs is hourly timestep)
  ! IF ((BaseTHDD-Temp_C)>=0) THEN   !Heating
  !    gamma1=1
  ! ELSE
  !    gamma1=0
  ! ENDIF
  !
  ! IF ((Temp_C-BaseTHDD)>=0) THEN   !Cooling
  !    gamma2=1
  ! ELSE
  !    gamma2=0
  ! ENDIF
  !
  ! IF(Gridiv == 1) tstepcount=tstepcount+1   !Add 1 to tstepcount only once for all grids
  !
  ! HDD(id,1)=HDD(id,1) + gamma1*(BaseTHDD-Temp_C)   !Heating
  ! HDD(id,2)=HDD(id,2) + gamma2*(Temp_C-BaseTHDD)   !Cooling
  ! HDD(id,3)=HDD(id,3) + Temp_C                     !Will become daily average temperature
  ! !      4 ------------------------------------!   !5-day running mean
  ! HDD(id,5)=HDD(id,5) + Precip                     !Daily precip total
  ! !      6 ------------------------------------!   !Days since rain




  ! Update snow density, albedo surface fraction
  IF (snowUse==1) CALL SnowUpdate(&
       nsurf,tstep,Temp_C,tau_a,tau_f,tau_r,&!input
       SnowDensMax,SnowDensMin,SnowAlbMin,SnowPack,&
       SnowAlb,SnowDens)!inout

  ! ================================================================================
  ! This next part occurs only on the first or last timestep of each day ===========

  ! --------------------------------------------------------------------------------
  ! On first timestep of each day, define whether the day each a workday or weekend
  IF (it==0.AND.imin==0) THEN

     !  CALL day2m  onth(id,mb,date,seas,iy,lat) !Calculate real date from doy
     !  CALL Day_of_Week(date,mb,iy,wd)        !Calculate weekday (1=Sun, ..., 7=Sat)
     !
     !  dayofWeek(id,1)=wd      !Day of week
     !  dayofWeek(id,2)=mb      !Month
     !  dayofweek(id,3)=seas    !Season

     CALL SUEWS_Cal_DailyStateStart(id,date,iy,lat,dayofWeek(id,:))

     ! --------------------------------------------------------------------------------
     ! On last timestep, perform the daily calculations -------------------------------
     ! Daily values not correct until end of each day,
     !  so main program should use values from the previous day
  ELSEIF (it==23 .AND. imin==(nsh_real-1)/nsh_real*60) THEN
     !write(*,*) 'Last timestep of day'

     !  ! Heating degree days (HDD) -------------
     !  HDD(id,1)=HDD(id,1)/tstepcount   !Heating
     !  HDD(id,2)=HDD(id,2)/tstepcount   !Cooling
     !  HDD(id,3)=HDD(id,3)/tstepcount   !Average temp
     !
     !  ! Calculate 5-day running mean temp     !!Need to deal with the previous year - CHECK!!
     !  DO jj=1,5
     !     HDD(id,4)=HDD(id,4) + HDD(id-(jj-1),3)
     !  ENDDO
     !  HDD(id,4) = HDD(id,4)/5
     !
     !  ! Calculate number of days since rain
     !  IF(HDD(id,5)>0) THEN        !Rain occurred
     !     HDD(id,6)=0
     !  ELSE
     !     HDD(id,6)=HDD(id-1,6)+1  !Days since rain
     !  ENDIF
     CALL update_HDD(&
          id,tstepcount,& !input
          HDD) !inout

     IF(Gridiv == NumberOfGrids) tstepcount=0  !Set to zero only after last grid has run

     ! Calculate modelled daily water use ------------------------------------------
     CALL update_WaterUse(&
          id,WaterUseMethod,DayofWeek,lat,Faut,HDD,&!input
          Ie_a,Ie_m,Ie_start,Ie_end,DayWatPer,DayWat,&
          WU_Day) !output
     !  IF (WaterUseMethod==0) THEN   !If water use is to be modelled (rather than observed)
     !
     !     wd=dayofWeek(id,1)
     !
     !     IF (DayWat(wd)==1.0) THEN      !1 indicates watering permitted on this day
     !        calc=0
     !        IF (lat>=0) THEN            !Northern Hemisphere
     !           IF (id>=Ie_start-1.AND.id<=Ie_end+1) calc=1   !Day between irrigation period
     !        ELSE                        !Southern Hemisphere
     !           calc=1
     !           IF (id>=Ie_end.AND.id<=Ie_start) calc=0       !Day between irrigation period
     !        ENDIF
     !
     !        IF(calc==1) THEN
     !           ! Model daily water use based on HDD(id,6)(days since rain) and HDD(id,3)(average temp)
     !           ! WU_Day is the amount of water [mm] per day, applied to each of the irrigated areas
     !           ! N.B. These are the same for each vegetation type at the moment
     !
     !           ! ---- Automatic irrigation (evergreen trees) ----
     !           WU_day(id,2) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,2)<0) WU_Day(id,2)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Manual irrigation (evergreen trees) ----
     !           WU_day(id,3) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,3)<0) WU_Day(id,3)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Total evergreen trees water use (automatic + manual) ----
     !           WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))
     !
     !           ! ---- Automatic irrigation (deciduous trees) ----
     !           WU_day(id,5) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,5)<0) WU_Day(id,5)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Manual irrigation (deciduous trees) ----
     !           WU_day(id,6) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,6)<0) WU_Day(id,6)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Total deciduous trees water use (automatic + manual) ----
     !           WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))
     !
     !           ! ---- Automatic irrigation (grass) ----
     !           WU_day(id,8) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,8)<0) WU_Day(id,8)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Manual irrigation (grass) ----
     !           WU_day(id,9) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
     !           IF (WU_Day(id,9)<0) WU_Day(id,9)=0   !If modelled WU is negative -> 0
     !
     !           ! ---- Total grass water use (automatic + manual) ----
     !           WU_Day(id,7)=(WU_day(id,8)+WU_day(id,9))
     !
     !        ELSE   !If no irrigation on this day
     !           WU_Day(id,1)=0
     !           WU_Day(id,2)=0
     !           WU_Day(id,3)=0
     !           WU_Day(id,4)=0
     !           WU_Day(id,5)=0
     !           WU_Day(id,6)=0
     !           WU_Day(id,7)=0
     !           WU_Day(id,8)=0
     !           WU_Day(id,9)=0
     !        ENDIF
     !     ENDIF
     !  ENDIF

     !------------------------------------------------------------------------------
     ! Calculation of LAI from growing degree days
     ! This was revised and checked on 16 Feb 2014 by LJ
     !------------------------------------------------------------------------------
     CALL update_GDDLAI(&
          id,LAICalcYes,& !input
          lat,LAI_obs,&
          BaseT,&
          BaseTe,&
          GDDFull,&
          SDDFull,&
          LAIMin,&
          LAIMax,&
          LAIPower,LAIType,&
          GDD,LAI) !inout

     !  critDays=50   !Critical limit for GDD when GDD or SDD is set to zero
     !
     !  ! Loop through vegetation types (iv)
     !  DO iv=1,NVegSurf
     !     ! Calculate GDD for each day from the minimum and maximum air temperature
     !     yes =((GDD(id,3)+GDD(id,4))/2-BaseT(iv))    !Leaf on
     !     no  =((GDD(id,3)+GDD(id,4))/2-BaseTe(iv))   !Leaf off
     !
     !     indHelp = 0   !Help switch to allow GDD to go to zero in sprint-time !!What does this mean?? HCW
     !
     !     IF(yes<0) THEN   !GDD cannot be negative
     !        indHelp=yes   !Amount of negative GDD
     !        yes=0
     !     ENDIF
     !
     !     IF(no>0) no=0    !SDD cannot be positive
     !
     !     ! Calculate cumulative growing and senescence degree days
     !     GDD(id,1) = GDD(id-1,1)+yes
     !     GDD(id,2) = GDD(id-1,2)+no
     !
     !     ! Possibility for cold spring
     !     IF(GDD(id,2)<=SDDFull(iv).AND.indHelp<0) THEN
     !        GDD(id,1)=0
     !     ENDIF
     !
     !     IF(GDD(id,1)>=GDDFull(iv)) THEN   !Start senescence
     !        GDD(id,1)=GDDFull(iv)          !Leaves should not grow so delete yes from earlier
     !        IF(GDD(id,2)<-critDays) GDD(id,1)=0
     !     ENDIF
     !
     !     IF (GDD(id,2)<=SDDFull(iv)) THEN   !After senescence now start growing leaves
     !        GDD(id,2)=SDDFull(iv)           !Leaves off so add back earlier
     !        IF(GDD(id,1)>critDays) GDD(id,2)=0
     !     ENDIF
     !
     !     ! With these limits SDD, GDD is set to zero
     !     IF(GDD(id,2)<-critDays.AND.GDD(id,2)>SDDFull(iv))  GDD(id,1)=0
     !     IF(GDD(id,1)> critDays.AND.GDD(id,1)<GDDFull(iv))  GDD(id,2)=0
     !
     !     ! Now calculate LAI itself
     !     IF(lat>=0) THEN   !Northern hemispere
     !        IF (id==140.AND.GDD(id,2)/=0)  GDD(id,2)=0  !If SDD is not zero by mid May, this is forced
     !        ! Set SDD to zero in summer time
     !        IF (GDD(id,1)> critDays.AND.id<170) GDD(id,2)=0
     !        ! Set GDD zero in winter time
     !        IF (GDD(id,2)<-critDays.AND.id>170) GDD(id,1)=0
     !
     !        IF (LAItype(iv) < 0.5) THEN   !Original LAI type
     !           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN       !Leaves can still grow
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
     !           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN   !Start senescence
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(3,iv)*GDD(id,2)*LAIPower(4,iv))+LAI(id-1,iv)
     !           ELSE
     !              LAI(id,iv)=LAI(id-1,iv)
     !           ENDIF
     !        ELSEIF (LAItype(iv)>=0.5) THEN
     !           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN        !Leaves can still grow
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
     !              !! Use day length to start senescence at high latitudes (N hemisphere)
     !           ELSEIF (GDD(id,5)<=12.AND.GDD(id,2)>SDDFull(iv)) THEN !Start senescence
     !              LAI(id,iv)=(LAI(id-1,iv)*LAIPower(3,iv)*(1-GDD(id,2))*LAIPower(4,iv))+LAI(id-1,iv)
     !           ELSE
     !              LAI(id,iv)=LAI(id-1,iv)
     !           ENDIF
     !        ENDIF
     !
     !     ELSEIF (lat<0) THEN   !Southern hemisphere !! N.B. not identical to N hemisphere - return to later
     !        IF (id==300.AND.GDD(id,2)/=0)  GDD(id,2)=0   !If SDD is not zero by late Oct, this is forced
     !        ! Set SDD to zero in summer time
     !        IF (GDD(id,1)> critDays.AND.id>250) GDD(id,2)=0
     !        ! Set GDD zero in winter time
     !        IF (GDD(id,2)<-critDays.AND.id<250) GDD(id,1)=0
     !
     !        IF (LAItype(iv) < 0.5) THEN   !Original LAI type
     !           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
     !           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(3,iv)*GDD(id,2)*LAIPower(4,iv))+LAI(id-1,iv)
     !           ELSE
     !              LAI(id,iv)=LAI(id-1,iv)
     !           ENDIF
     !        ELSE
     !           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN
     !              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
     !              !! Day length not used to start senescence in S hemisphere (not much land)
     !           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN
     !              LAI(id,iv)=(LAI(id-1,iv)*LAIPower(3,iv)*(1-GDD(id,2))*LAIPower(4,iv))+LAI(id-1,iv)
     !           ELSE
     !              LAI(id,iv)=LAI(id-1,iv)
     !           ENDIF
     !        ENDIF
     !     ENDIF   !N or S hemisphere
     !
     !     ! Check LAI within limits; if not set to limiting value
     !     IF(LAI(id,iv).GT.LAImax(iv))THEN
     !        LAI(id,iv)=LAImax(iv)
     !     ELSEIF(LAI(id,iv).LT.LAImin(iv))THEN
     !        LAI(id,iv)=LAImin(iv)
     !     ENDIF
     !
     !  ENDDO   !End of loop over veg surfaces
     !
     !  IF(LAICalcYes==0)THEN ! moved to DailyState, TS 18 Sep 2017
     !     LAI(id-1,:)=LAI_obs ! check -- this is going to be a problem as it is not for each vegetation class
     !  ENDIF
     !  !------------------------------------------------------------------------------

     CALL update_Veg(&
          id,&!input
          LAImax,&
          LAIMin,&
          alBMax_DecTr,&
          alBMax_EveTr,&
          alBMax_Grass,&
          AlbMin_DecTr,&
          AlbMin_EveTr,&
          AlbMin_Grass,&
          CapMax_dec,&
          CapMin_dec,&
          PorMax_dec,&
          PorMin_dec,&
          DecidCap,&!inout
          albDecTr,&
          albEveTr,&
          albGrass,&
          porosity,&
          LAI,&
          deltaLAI)!output

     !  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     !  ! Calculate the development of vegetation cover
     !  ! Albedo changes with LAI for each vegetation type
     !  ! Storage capacity and porosity are updated based on DecTr LAI only (seasonal variation in Grass and EveTr assumed small)
     !  ! If only LUMPS is used, set deciduous capacities to 0
     !  ! Assume porosity Change based on GO99- Heisler??
     !  deltaLAI=0
     !  deltaLAIEveTr=0
     !  deltaLAIGrass=0
     !  CapChange=0
     !  porChange=0
     !  albChangeDecTr=0
     !  albChangeEveTr=0
     !  albChangeGrass=0
     !
     !  iv=ivDecid
     !  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     !     deltaLAI=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     !     albChangeDecTr=(alBMax_DecTr-AlbMin_DecTr)* deltaLAI
     !     CapChange=(CapMin_dec-CapMax_dec)* deltaLAI
     !     porChange=(PorMin_dec-PorMax_dec)* deltaLAI
     !  ENDIF
     !
     !  iv=ivConif
     !  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     !     deltaLAIEveTr=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     !     albChangeEveTr=(alBMax_EveTr-AlbMin_EveTr)* deltaLAIEveTr    !!N.B. Currently uses deltaLAI for deciduous trees only!!
     !  ENDIF
     !
     !  iv=ivGrass
     !  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     !     deltaLAIGrass=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     !     albChangeGrass=(alBMax_Grass-AlbMin_Grass)* deltaLAIGrass    !!N.B. Currently uses deltaLAI for deciduous trees only!!
     !  ENDIF
     !
     !  iv=ivDecid
     !
     !  !write(*,*) deltaLAI, deltaLAIEveTr, deltaLAIGrass
     !
     !  DecidCap(id) = DecidCap(id-1) - CapChange
     !  albDecTr(id) = albDecTr(id-1) + albChangeDecTr
     !  porosity(id) = porosity(id-1) + porChange !- changed to + by HCW 20 Aug 2015 (porosity greatest when LAI smallest)
     !  !Also update albedo of EveTr and Grass surfaces
     !  albEveTr(id) = albEveTr(id-1) + albChangeEveTr
     !  albGrass(id) = albGrass(id-1) + albChangeGrass

     CALL update_AnOHM(&
          Gridiv,id,& !input
          ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids,dataOut,&
          a1,a2,a3,&!inout
          xBo,xmAH) !output

     !  ! ------- AnOHM related, added by TS 01 Mar 2016 ------------------------------
     !  !   read in QH and QE of the current day
     !  ALLOCATE(lSub(SIZE(dataOut(:,2,Gridiv))))
     !  WHERE (dataOut(:,2,Gridiv) == id)
     !     lSub = 1
     !  ELSEWHERE
     !     lSub = 0
     !  END WHERE
     !  irRange = (/MAXLOC(lSub),MAXLOC(lSub)+SUM(lSub)/)
     !  ! write(*,*) 'good here: 424',irRange(2), irRange(1),sum(lSub),maxloc(lSub)
     !  !   calculate Bowen ratio
     !  !   check if the diurnal cycle of present day (id) is complete
     !  IF (irRange(2)- irRange(1) >= nsh*24-1) THEN
     !     !   load the sublist into forcings:
     !     !   QH and QE for Bowen ratio
     !
     !     ! needs improvement:
     !     ! TODO: include all timesteps as in the calculations
     !     xQH  = dataOut(irRange(1):irRange(2):nsh,14,Gridiv)
     !     xQE  = dataOut(irRange(1):irRange(2):nsh,15,Gridiv)
     !     mxQH = SUM(xQH(10:16))/7
     !     mxQE = SUM(xQE(10:16))/7
     !
     !
     !     ! handle extreme dry condition to prevent NAN
     !     IF ( ABS(mxQE) < 0.1 ) mxQE = 0.1
     !
     !     xBo  = mxQH/mxQE
     !     !   calculate daily mean AH
     !     xAH  = dataOut(irRange(1):irRange(2):nsh,15,Gridiv)
     !     xmAH  = SUM(xAH(:))/24
     !  ELSE
     !     !   give a default Bo as 2 if no enough flux data to update with:
     !     xBo  = 10.
     !     xmAH = 25.
     !  END IF
     !  !   handle NAN
     !  IF ( xBo/=xBo ) THEN
     !     xBo = 2.
     !  END IF
     !  IF ( xmAH/=xmAH ) THEN
     !     xmAH = 25.
     !  END IF
     !
     !
     !  BoAnOHMEnd(Gridiv)  = xBo
     !  Bo_grids(id,Gridiv) = BoAnOHMEnd(Gridiv)
     !  BoInit              = Bo_grids(id,Gridiv) ! BoInit will be written out to InitialConditions of next year
     !
     !  IF(debug) THEN
     !     IF ( ABS(xBo) > 30. ) WRITE(unit=*, fmt=*) "mean QH, QE:", mxQH,mxQE
     !     IF ( ABS(xBo) < 1/30. ) WRITE(unit=*, fmt=*) "mean QH, QE:", mxQH,mxQE
     !     IF ( HDD(id,5)>0 ) PRINT*, 'Rainfall',HDD(id,5)
     !     PRINT*, 'xBo is',xBo
     !  ENDIF
     !
     !  mAHAnOHM(Gridiv)     = xmAH
     !  mAH_grids(id,Gridiv) = mAHAnOHM(Gridiv)
     !  ! load current AnOHM coef.:
     !  !  IF ( StorageHeatMethod==3 ) THEN
     !  !     !  if ( id>364 ) then
     !  !     !    print*, 'test in DailyState'
     !  !     !    print*, a1AnOHM(Gridiv),a2AnOHM(Gridiv),a3AnOHM(Gridiv)
     !  !     !
     !  !     !  end if
     !  !     a1AnOHM(Gridiv) = a1AnOHM_grids(id,Gridiv)
     !  !     a2AnOHM(Gridiv) = a2AnOHM_grids(id,Gridiv)
     !  !     a3AnOHM(Gridiv) = a3AnOHM_grids(id,Gridiv)
     !  !  ELSE
     !  ! output the normal OHM coeff.:
     !  a1AnOHM(Gridiv) = a1
     !  a2AnOHM(Gridiv) = a2
     !  a3AnOHM(Gridiv) = a3
     !  !  END IF
     !
     !  ! --------- AnOHM related End ----------------------------------------------------
    !  PRINT*,'OK here'
     ! -----------------------------------------------------------------------------
print*, 'VegPhenLumps here 1',VegPhenLumps
     CALL update_DailyState(&
          iy,id,&!input
          GDD,HDD,LAI,&
          DecidCap,albDecTr,albEveTr,albGrass,porosity,&
          WU_Day,&
          deltaLAI,VegPhenLumps,&
          SnowAlb,SnowDens,&
          xBo,a1,a2,a3,&
          DailyStateLine)!output

     CALL output_DailyState(&
          Gridiv,GridIDmatrix,&!input
          FileCode,FileOutputPath,&
          DailyStateLine,&
          DailyStateFirstOpen)


  ENDIF   !End of section done only at the end of each day (i.e. only once per day)
  ! ================================================================================

  RETURN

END SUBROUTINE DailyState


SUBROUTINE init_DailyState(&
     id,& !input
     Gridiv,&
     avkdn,&
     Temp_C,&
     Precip,&
     BaseTHDD,&
     nsh_real,&
     tstepcount,& !inout
     GDD,&
     HDD)
  IMPLICIT NONE
  INTEGER,PARAMETER::ndays    = 366

  INTEGER,INTENT(in)::id
  INTEGER,INTENT(in)::Gridiv
  REAL(KIND(1d0)),INTENT(in)::avkdn
  REAL(KIND(1d0)),INTENT(in)::Temp_C
  REAL(KIND(1d0)),INTENT(in)::Precip
  REAL(KIND(1d0)),INTENT(in)::BaseTHDD
  REAL(KIND(1d0)),INTENT(in)::nsh_real

  REAL(KIND(1d0)),INTENT(inout)::tstepcount
  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5),INTENT(inout):: GDD !Growing Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, 6),INTENT(inout):: HDD          !Heating Degree Days (see SUEWS_DailyState.f95)

  INTEGER::gamma1
  INTEGER::gamma2
  !! Initialization -----------------------------------------------------------------
  !! These variables don't seem to be needed (commented out HCW 27 Nov 2014)
  !! If required, they will need updating for a non-hourly timestep
  !runT(it)=Temp_C      !runT has been initialized in SUEWS_initial to the previous day average
  !avT_h=sum(runT)/24   !Average daily temperature
  !runP(it)=Precip   !Same for precipitation
  !totP_h=sum(runP)     !Daily sum for precipitation

  ! Daily min and max temp (these get updated through the day) ---------------------
  GDD(id,3) = MIN(Temp_C,GDD(id,3))     !Daily min T in column 3
  GDD(id,4) = MAX(Temp_C,GDD(id,4))     !Daily max T in column 4
  IF (avkdn>10) THEN
     GDD(id,5) = GDD(id,5)+1/nsh_real   !Cumulate daytime hours !Divide by nsh (HCW 01 Dec 2014)
  ENDIF

  ! Calculations related to heating and cooling degree days (HDD) ------------------
  ! See Sailor & Vasireddy (2006) EMS Eq 1,2 (theirs is hourly timestep)
  IF ((BaseTHDD-Temp_C)>=0) THEN   !Heating
     gamma1=1
  ELSE
     gamma1=0
  ENDIF

  IF ((Temp_C-BaseTHDD)>=0) THEN   !Cooling
     gamma2=1
  ELSE
     gamma2=0
  ENDIF

  IF(Gridiv == 1) tstepcount=tstepcount+1   !Add 1 to tstepcount only once for all grids

  HDD(id,1)=HDD(id,1) + gamma1*(BaseTHDD-Temp_C)   !Heating
  HDD(id,2)=HDD(id,2) + gamma2*(Temp_C-BaseTHDD)   !Cooling
  HDD(id,3)=HDD(id,3) + Temp_C                     !Will become daily average temperature
  !      4 ------------------------------------!   !5-day running mean
  HDD(id,5)=HDD(id,5) + Precip                     !Daily precip total
  !      6 ------------------------------------!   !Days since rain

END SUBROUTINE init_DailyState


SUBROUTINE update_AnOHM(&
     Gridiv,id,& !input
     ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids,dataOut,&
     a1,a2,a3,&!inout
     xBo,xmAH) !output
  IMPLICIT NONE
  INTEGER,INTENT(in) :: Gridiv,id,ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids

  REAL(KIND(1d0)),DIMENSION(ReadLinesMetdata,ncolumnsDataOut,NumberOfGrids),INTENT(in)::dataOut
  REAL(KIND(1d0)),INTENT(inout)::a1,a2,a3
  REAL(KIND(1d0)),INTENT(out)::xBo,xmAH
  !   local variables:
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE :: subMet ! subset array of daytime series

  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE :: xQH,xQE,xAH
  REAL(KIND(1d0))::mxQH,mxQE,xa1,xa2,xa3

  INTEGER :: err
  INTEGER :: lenMetData,nVar


  LOGICAL, ALLOCATABLE :: metMask(:)


  ! determine the length of subset
  lenMetData = COUNT(&
       dataOut(:,2,Gridiv)==id & ! day=xid
       .AND. dataOut(:,4,Gridiv)==0)! tmin=0

  ! construct mask
  IF (ALLOCATED(metMask)) DEALLOCATE(metMask, stat=err)
  ALLOCATE(metMask(lenMetData))
  metMask=(dataOut(:,2,Gridiv)==id & ! day=xid
       .AND. dataOut(:,4,Gridiv)==0)! tmin=0

  ! construct array for time and met variables
  nVar=3! number of variables to retrieve
  ! print*, 'good 1'
  ! allocate subMet:
  IF (ALLOCATED(subMet)) DEALLOCATE(subMet, stat=err)
  ALLOCATE(subMet(lenMetData,nVar))
  subMet=RESHAPE(PACK(dataOut(:,(/14,15,16/),Gridiv),&! QH,QE,AH
       SPREAD(metMask, dim=2, ncopies=nVar)),& ! replicate mask vector to 2D array
       (/lenMetData,nVar/)) ! convert to target shape

  ! re-allocate arrays as their sizes may change during passing
  IF (ALLOCATED(xQH)) DEALLOCATE(xQH, stat=err)
  ALLOCATE(xQH(lenMetData))
  IF (ALLOCATED(xQE)) DEALLOCATE(xQE, stat=err)
  ALLOCATE(xQE(lenMetData))
  IF (ALLOCATED(xAH)) DEALLOCATE(xAH, stat=err)
  ALLOCATE(xAH(lenMetData))

  xQH  = subMet(:,1)
  xQE  = subMet(:,2)
  mxQH = SUM(xQH(10:16))/7
  mxQE = SUM(xQE(10:16))/7
  ! handle extreme dry condition to prevent NAN
  IF ( ABS(mxQE) < 0.1 ) mxQE = 0.1
  xBo  = mxQH/mxQE

  !   calculate daily mean AH
  xAH  = subMet(:,3)
  xmAH  = SUM(xAH(:))/24

  xa1=a1
  xa2=a2
  xa3=a3

  ! BoAnOHMEnd(Gridiv)  = xBo
  ! Bo_grids(id,Gridiv) = BoAnOHMEnd(Gridiv)
  ! BoInit              = Bo_grids(id,Gridiv) ! BoInit will be written out to InitialConditions of next year

  ! IF(debug) THEN
  !    IF ( ABS(xBo) > 30. ) WRITE(unit=*, fmt=*) "mean QH, QE:", mxQH,mxQE
  !    IF ( ABS(xBo) < 1/30. ) WRITE(unit=*, fmt=*) "mean QH, QE:", mxQH,mxQE
  !    IF ( HDD(id,5)>0 ) PRINT*, 'Rainfall',HDD(id,5)
  !    PRINT*, 'xBo is',xBo
  ! ENDIF

  ! mAHAnOHM(Gridiv)     = xmAH
  ! mAH_grids(id,Gridiv) = mAHAnOHM(Gridiv)
  ! load current AnOHM coef.:
  !  IF ( StorageHeatMethod==3 ) THEN
  !     !  if ( id>364 ) then
  !     !    print*, 'test in DailyState'
  !     !    print*, a1AnOHM(Gridiv),a2AnOHM(Gridiv),a3AnOHM(Gridiv)
  !     !
  !     !  end if
  !     a1AnOHM(Gridiv) = a1AnOHM_grids(id,Gridiv)
  !     a2AnOHM(Gridiv) = a2AnOHM_grids(id,Gridiv)
  !     a3AnOHM(Gridiv) = a3AnOHM_grids(id,Gridiv)
  !  ELSE
  ! output the normal OHM coeff.:
  ! a1AnOHM(Gridiv) = a1
  ! a2AnOHM(Gridiv) = a2
  ! a3AnOHM(Gridiv) = a3

END SUBROUTINE update_AnOHM



SUBROUTINE update_Veg(&
     id,&!input
     LAImax,&
     LAIMin,&
     alBMax_DecTr,&
     alBMax_EveTr,&
     alBMax_Grass,&
     AlbMin_DecTr,&
     AlbMin_EveTr,&
     AlbMin_Grass,&
     CapMax_dec,&
     CapMin_dec,&
     PorMax_dec,&
     PorMin_dec,&
     DecidCap,&!inout
     albDecTr,&
     albEveTr,&
     albGrass,&
     porosity,&
     LAI,&
     deltaLAI)!output

  IMPLICIT NONE
  INTEGER,PARAMETER::ndays    = 366
  INTEGER,PARAMETER::nvegsurf = 3
  INTEGER,PARAMETER::ivConif  = 1 !When only vegetated surfaces considered (1-3)
  INTEGER,PARAMETER::ivDecid  = 2
  INTEGER,PARAMETER::ivGrass  = 3

  INTEGER,INTENT(in)::id
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::LAImax
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)::LAIMin

  REAL(KIND(1d0)),INTENT(in)::alBMax_DecTr
  REAL(KIND(1d0)),INTENT(in)::alBMax_EveTr
  REAL(KIND(1d0)),INTENT(in)::alBMax_Grass
  REAL(KIND(1d0)),INTENT(in)::AlbMin_DecTr
  REAL(KIND(1d0)),INTENT(in)::AlbMin_EveTr
  REAL(KIND(1d0)),INTENT(in)::AlbMin_Grass
  REAL(KIND(1d0)),INTENT(in)::CapMax_dec
  REAL(KIND(1d0)),INTENT(in)::CapMin_dec
  REAL(KIND(1d0)),INTENT(in)::PorMax_dec
  REAL(KIND(1d0)),INTENT(in)::PorMin_dec

  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::DecidCap
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albDecTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albEveTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::albGrass
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(inout)::porosity

  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf),INTENT(inout)::LAI
  REAL(KIND(1d0)),INTENT(out)::deltaLAI

  INTEGER::iv

  REAL(KIND(1d0))::albChangeDecTr
  REAL(KIND(1d0))::albChangeEveTr
  REAL(KIND(1d0))::albChangeGrass
  REAL(KIND(1d0))::CapChange

  REAL(KIND(1d0))::deltaLAIEveTr
  REAL(KIND(1d0))::deltaLAIGrass
  REAL(KIND(1d0))::porChange
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ! Calculate the development of vegetation cover
  ! Albedo changes with LAI for each vegetation type
  ! Storage capacity and porosity are updated based on DecTr LAI only (seasonal variation in Grass and EveTr assumed small)
  ! If only LUMPS is used, set deciduous capacities to 0
  ! Assume porosity Change based on GO99- Heisler??
  deltaLAI=0
  deltaLAIEveTr=0
  deltaLAIGrass=0
  CapChange=0
  porChange=0
  albChangeDecTr=0
  albChangeEveTr=0
  albChangeGrass=0

  iv=ivDecid
  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     deltaLAI=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     albChangeDecTr=(alBMax_DecTr-AlbMin_DecTr)* deltaLAI
     CapChange=(CapMin_dec-CapMax_dec)* deltaLAI
     porChange=(PorMin_dec-PorMax_dec)* deltaLAI
  ENDIF

  iv=ivConif
  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     deltaLAIEveTr=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     albChangeEveTr=(alBMax_EveTr-AlbMin_EveTr)* deltaLAIEveTr    !!N.B. Currently uses deltaLAI for deciduous trees only!!
  ENDIF

  iv=ivGrass
  IF((LAI(ID,iv)-LAI(ID-1,iv))/=0) THEN
     deltaLAIGrass=(LAI(id,iv)-LAI(id-1,iv))/(LAImax(iv)-LAIMin(iv))
     albChangeGrass=(alBMax_Grass-AlbMin_Grass)* deltaLAIGrass    !!N.B. Currently uses deltaLAI for deciduous trees only!!
  ENDIF

  iv=ivDecid

  !write(*,*) deltaLAI, deltaLAIEveTr, deltaLAIGrass

  DecidCap(id) = DecidCap(id-1) - CapChange
  albDecTr(id) = albDecTr(id-1) + albChangeDecTr
  porosity(id) = porosity(id-1) + porChange !- changed to + by HCW 20 Aug 2015 (porosity greatest when LAI smallest)
  !Also update albedo of EveTr and Grass surfaces
  albEveTr(id) = albEveTr(id-1) + albChangeEveTr
  albGrass(id) = albGrass(id-1) + albChangeGrass

END SUBROUTINE update_Veg



SUBROUTINE update_GDDLAI(&
     id,LAICalcYes,& !input
     lat,LAI_obs,&
     BaseT,&
     BaseTe,&
     GDDFull,&
     SDDFull,&
     LAIMin,&
     LAIMax,&
     LAIPower,LAIType,&
     GDD,LAI) !inout
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  ! Calculation of LAI from growing degree days
  ! This was revised and checked on 16 Feb 2014 by LJ
  !------------------------------------------------------------------------------
  INTEGER, PARAMETER:: ndays = 366   !Max no. days in a year used to specify size of daily arrays
  ! INTEGER, PARAMETER:: nsurf=7                !Total number of surfaces
  INTEGER, PARAMETER:: NVegSurf=3             !Number of surfaces that are vegetated
  ! INTEGER, PARAMETER:: nsurfIncSnow=nsurf+1   !Number of surfaces + snow
  ! INTEGER, PARAMETER:: MaxNumberOfGrids=2000   !Max no. grids   !HCW changed to 2000 from 10000 so prog can run on windows (2GB lim)
  ! INTEGER, PARAMETER:: MaxLinesMet=8640        !Max no. lines to read in one go (for all grids, ie MaxLinesMet/NumberOfGrids each)

  INTEGER,INTENT(in)::id
  INTEGER,INTENT(in)::LAICalcYes

  REAL(KIND(1d0)),INTENT(in)::lat
  REAL(KIND(1d0)),INTENT(in)::LAI_obs

  ! --- Vegetation phenology ---------------------------------------------------------------------
  ! Parameters provided in input information for each vegetation surface (SUEWS_Veg.txt)
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: BaseT          !Base temperature for growing degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: BaseTe         !Base temperature for senescence degree days [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: GDDFull        !Growing degree days needed for full capacity [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: SDDFull        !Senescence degree days needed to initiate leaf off [degC]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: LAIMin         !Min LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(nvegsurf),INTENT(in)  :: LAIMax         !Max LAI [m2 m-2]
  REAL(KIND(1d0)),DIMENSION(4,nvegsurf),INTENT(in):: LAIPower       !Coeffs for LAI equation: 1,2 - leaf growth; 3,4 - leaf off
  !! N.B. currently DecTr only, although input provided for all veg types
  INTEGER,DIMENSION(nvegsurf),INTENT(in):: LAIType                  !LAI equation to use: original (0) or new (1)

  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5),INTENT(inout)       :: GDD !Growing Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf),INTENT(inout):: LAI !LAI for each veg surface [m2 m-2]

  REAL(KIND(1d0)):: no,yes,indHelp   !Switches and checks for GDD

  INTEGER :: critDays,iv


  critDays=50   !Critical limit for GDD when GDD or SDD is set to zero

  ! Loop through vegetation types (iv)
  DO iv=1,NVegSurf
     ! Calculate GDD for each day from the minimum and maximum air temperature
     yes =((GDD(id,3)+GDD(id,4))/2-BaseT(iv))    !Leaf on
     no  =((GDD(id,3)+GDD(id,4))/2-BaseTe(iv))   !Leaf off

     indHelp = 0   !Help switch to allow GDD to go to zero in sprint-time !!What does this mean?? HCW

     IF(yes<0) THEN   !GDD cannot be negative
        indHelp=yes   !Amount of negative GDD
        yes=0
     ENDIF

     IF(no>0) no=0    !SDD cannot be positive

     ! Calculate cumulative growing and senescence degree days
     GDD(id,1) = GDD(id-1,1)+yes
     GDD(id,2) = GDD(id-1,2)+no

     ! Possibility for cold spring
     IF(GDD(id,2)<=SDDFull(iv).AND.indHelp<0) THEN
        GDD(id,1)=0
     ENDIF

     IF(GDD(id,1)>=GDDFull(iv)) THEN   !Start senescence
        GDD(id,1)=GDDFull(iv)          !Leaves should not grow so delete yes from earlier
        IF(GDD(id,2)<-critDays) GDD(id,1)=0
     ENDIF

     IF (GDD(id,2)<=SDDFull(iv)) THEN   !After senescence now start growing leaves
        GDD(id,2)=SDDFull(iv)           !Leaves off so add back earlier
        IF(GDD(id,1)>critDays) GDD(id,2)=0
     ENDIF

     ! With these limits SDD, GDD is set to zero
     IF(GDD(id,2)<-critDays.AND.GDD(id,2)>SDDFull(iv))  GDD(id,1)=0
     IF(GDD(id,1)> critDays.AND.GDD(id,1)<GDDFull(iv))  GDD(id,2)=0

     ! Now calculate LAI itself
     IF(lat>=0) THEN   !Northern hemispere
        IF (id==140.AND.GDD(id,2)/=0)  GDD(id,2)=0  !If SDD is not zero by mid May, this is forced
        ! Set SDD to zero in summer time
        IF (GDD(id,1)> critDays.AND.id<170) GDD(id,2)=0
        ! Set GDD zero in winter time
        IF (GDD(id,2)<-critDays.AND.id>170) GDD(id,1)=0

        IF (LAItype(iv) < 0.5) THEN   !Original LAI type
           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN       !Leaves can still grow
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN   !Start senescence
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(3,iv)*GDD(id,2)*LAIPower(4,iv))+LAI(id-1,iv)
           ELSE
              LAI(id,iv)=LAI(id-1,iv)
           ENDIF
        ELSEIF (LAItype(iv)>=0.5) THEN
           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN        !Leaves can still grow
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
              !! Use day length to start senescence at high latitudes (N hemisphere)
           ELSEIF (GDD(id,5)<=12.AND.GDD(id,2)>SDDFull(iv)) THEN !Start senescence
              LAI(id,iv)=(LAI(id-1,iv)*LAIPower(3,iv)*(1-GDD(id,2))*LAIPower(4,iv))+LAI(id-1,iv)
           ELSE
              LAI(id,iv)=LAI(id-1,iv)
           ENDIF
        ENDIF

     ELSEIF (lat<0) THEN   !Southern hemisphere !! N.B. not identical to N hemisphere - return to later
        IF (id==300.AND.GDD(id,2)/=0)  GDD(id,2)=0   !If SDD is not zero by late Oct, this is forced
        ! Set SDD to zero in summer time
        IF (GDD(id,1)> critDays.AND.id>250) GDD(id,2)=0
        ! Set GDD zero in winter time
        IF (GDD(id,2)<-critDays.AND.id<250) GDD(id,1)=0

        IF (LAItype(iv) < 0.5) THEN   !Original LAI type
           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(3,iv)*GDD(id,2)*LAIPower(4,iv))+LAI(id-1,iv)
           ELSE
              LAI(id,iv)=LAI(id-1,iv)
           ENDIF
        ELSE
           IF(GDD(id,1)>0.AND.GDD(id,1)<GDDFull(iv)) THEN
              LAI(id,iv)=(LAI(id-1,iv)**LAIPower(1,iv)*GDD(id,1)*LAIPower(2,iv))+LAI(id-1,iv)
              !! Day length not used to start senescence in S hemisphere (not much land)
           ELSEIF(GDD(id,2)<0.AND.GDD(id,2)>SDDFull(iv)) THEN
              LAI(id,iv)=(LAI(id-1,iv)*LAIPower(3,iv)*(1-GDD(id,2))*LAIPower(4,iv))+LAI(id-1,iv)
           ELSE
              LAI(id,iv)=LAI(id-1,iv)
           ENDIF
        ENDIF
     ENDIF   !N or S hemisphere

     ! Check LAI within limits; if not set to limiting value
     IF(LAI(id,iv).GT.LAImax(iv))THEN
        LAI(id,iv)=LAImax(iv)
     ELSEIF(LAI(id,iv).LT.LAImin(iv))THEN
        LAI(id,iv)=LAImin(iv)
     ENDIF

  ENDDO   !End of loop over veg surfaces

  IF(LAICalcYes==0)THEN ! moved to DailyState, TS 18 Sep 2017
     LAI(id-1,:)=LAI_obs ! check -- this is going to be a problem as it is not for each vegetation class
  ENDIF
  !------------------------------------------------------------------------------

END SUBROUTINE update_GDDLAI



SUBROUTINE update_WaterUse(&
     id,WaterUseMethod,DayofWeek,lat,Faut,HDD,&!input
     Ie_a,Ie_m,Ie_start,Ie_end,DayWatPer,DayWat,&
     WU_Day) !output

  IMPLICIT NONE
  INTEGER,PARAMETER :: ndays=366

  INTEGER,INTENT(in) :: id
  INTEGER,INTENT(in) :: WaterUseMethod
  INTEGER,INTENT(in)::Ie_start   !Starting time of water use (DOY)
  INTEGER,INTENT(in)::Ie_end       !Ending time of water use (DOY)
  INTEGER,DIMENSION(0:ndays,3),INTENT(in)::DayofWeek

  REAL(KIND(1d0)),INTENT(in)::lat
  REAL(KIND(1d0)),INTENT(in)::Faut          !Fraction of irrigated area using automatic irrigation

  REAL(KIND(1d0)),DIMENSION(-4:366,6),INTENT(in)::HDD
  REAL(KIND(1d0)),DIMENSION(3),INTENT(in)::Ie_a
  REAL(KIND(1d0)),DIMENSION(3),INTENT(in)::Ie_m   !Coefficients for automatic and manual irrigation models
  REAL(KIND(1d0)),DIMENSION(7),INTENT(in)::DayWatPer  !% of houses following daily water
  REAL(KIND(1d0)),DIMENSION(7),INTENT(in)::DayWat       !Days of watering allowed


  REAL(KIND(1d0)),DIMENSION(0:ndays,9),INTENT(out):: WU_Day       !Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)

  INTEGER::wd,&
       calc        !Water use calculation is done when calc = 1






  IF (WaterUseMethod==0) THEN   !If water use is to be modelled (rather than observed)

     wd=dayofWeek(id,1)

     IF (DayWat(wd)==1.0) THEN      !1 indicates watering permitted on this day
        calc=0
        IF (lat>=0) THEN            !Northern Hemisphere
           IF (id>=Ie_start-1.AND.id<=Ie_end+1) calc=1   !Day between irrigation period
        ELSE                        !Southern Hemisphere
           calc=1
           IF (id>=Ie_end.AND.id<=Ie_start) calc=0       !Day between irrigation period
        ENDIF

        IF(calc==1) THEN
           ! Model daily water use based on HDD(id,6)(days since rain) and HDD(id,3)(average temp)
           ! WU_Day is the amount of water [mm] per day, applied to each of the irrigated areas
           ! N.B. These are the same for each vegetation type at the moment

           ! ---- Automatic irrigation (evergreen trees) ----
           WU_day(id,2) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,2)<0) WU_Day(id,2)=0   !If modelled WU is negative -> 0

           ! ---- Manual irrigation (evergreen trees) ----
           WU_day(id,3) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,3)<0) WU_Day(id,3)=0   !If modelled WU is negative -> 0

           ! ---- Total evergreen trees water use (automatic + manual) ----
           WU_Day(id,1)=(WU_day(id,2)+WU_day(id,3))

           ! ---- Automatic irrigation (deciduous trees) ----
           WU_day(id,5) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,5)<0) WU_Day(id,5)=0   !If modelled WU is negative -> 0

           ! ---- Manual irrigation (deciduous trees) ----
           WU_day(id,6) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,6)<0) WU_Day(id,6)=0   !If modelled WU is negative -> 0

           ! ---- Total deciduous trees water use (automatic + manual) ----
           WU_Day(id,4)=(WU_day(id,5)+WU_day(id,6))

           ! ---- Automatic irrigation (grass) ----
           WU_day(id,8) = Faut*(Ie_a(1)+Ie_a(2)*HDD(id,3)+Ie_a(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,8)<0) WU_Day(id,8)=0   !If modelled WU is negative -> 0

           ! ---- Manual irrigation (grass) ----
           WU_day(id,9) = (1-Faut)*(Ie_m(1)+Ie_m(2)*HDD(id,3)+Ie_m(3)*HDD(id,6))*DayWatPer(wd)
           IF (WU_Day(id,9)<0) WU_Day(id,9)=0   !If modelled WU is negative -> 0

           ! ---- Total grass water use (automatic + manual) ----
           WU_Day(id,7)=(WU_day(id,8)+WU_day(id,9))

        ELSE   !If no irrigation on this day
           WU_Day(id,1)=0
           WU_Day(id,2)=0
           WU_Day(id,3)=0
           WU_Day(id,4)=0
           WU_Day(id,5)=0
           WU_Day(id,6)=0
           WU_Day(id,7)=0
           WU_Day(id,8)=0
           WU_Day(id,9)=0
        ENDIF
     ENDIF
  ENDIF

END SUBROUTINE update_WaterUse



SUBROUTINE update_HDD(&
     id,tstepcount,& !input
     HDD) !inout
  IMPLICIT NONE
  INTEGER,INTENT(in)::id
  REAL(KIND(1d0)),INTENT(in)::tstepcount
  REAL(KIND(1d0)),DIMENSION(-4:366,6),INTENT(inout):: HDD

  INTEGER :: jj


  ! Heating degree days (HDD) -------------
  HDD(id,1)=HDD(id,1)/tstepcount   !Heating
  HDD(id,2)=HDD(id,2)/tstepcount   !Cooling
  HDD(id,3)=HDD(id,3)/tstepcount   !Average temp

  ! Calculate 5-day running mean temp     !!Need to deal with the previous year - CHECK!!
  DO jj=1,5
     HDD(id,4)=HDD(id,4) + HDD(id-(jj-1),3)
  ENDDO
  HDD(id,4) = HDD(id,4)/5

  ! Calculate number of days since rain
  IF(HDD(id,5)>0) THEN        !Rain occurred
     HDD(id,6)=0
  ELSE
     HDD(id,6)=HDD(id-1,6)+1  !Days since rain
  ENDIF

END SUBROUTINE update_HDD



SUBROUTINE SUEWS_Cal_DailyStateEnd

  IMPLICIT NONE

END SUBROUTINE SUEWS_Cal_DailyStateEnd



SUBROUTINE SUEWS_Cal_DailyStateStart(&
     id,date,iy,lat,&!input
     dayofWeek_id)
  IMPLICIT NONE
  INTEGER,INTENT(in) :: id,date,iy
  REAL(KIND(1d0)),INTENT(in) ::lat

  INTEGER,DIMENSION(3),INTENT(out) ::dayofWeek_id

  INTEGER ::wd,mb,seas

  CALL day2month(id,mb,date,seas,iy,lat) !Calculate real date from doy
  CALL Day_of_Week(date,mb,iy,wd)        !Calculate weekday (1=Sun, ..., 7=Sat)

  dayofWeek_id(1)=wd      !Day of week
  dayofWeek_id(2)=mb      !Month
  dayofweek_id(3)=seas    !Season

END SUBROUTINE SUEWS_Cal_DailyStateStart

! transfer results to a one-line output for DailyState
SUBROUTINE update_DailyState(&
     iy,id,&!input
     GDD,HDD,LAI,&
     DecidCap,albDecTr,albEveTr,albGrass,porosity,&
     WU_Day,&
     deltaLAI,VegPhenLumps,&
     SnowAlb,SnowDens,&
     xBo,a1,a2,a3,&
     DailyStateLine)!output

  IMPLICIT NONE
  INTEGER,PARAMETER::ndays=366
  INTEGER,PARAMETER::nvegsurf=3

  INTEGER,INTENT(in) ::iy,id

  REAL(KIND(1d0)),DIMENSION( 0:ndays, 5),INTENT(in):: GDD          !Growing Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, 6),INTENT(in):: HDD          !Heating Degree Days (see SUEWS_DailyState.f95)
  REAL(KIND(1d0)),DIMENSION(-4:ndays, nvegsurf),INTENT(in):: LAI   !LAI for each veg surface [m2 m-2]

  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(in) ::DecidCap
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(in) ::albDecTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(in) ::albEveTr
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(in) ::albGrass
  REAL(KIND(1d0)),DIMENSION( 0:ndays),INTENT(in) ::porosity
  REAL(KIND(1d0)),DIMENSION(0:ndays,9),INTENT(in):: WU_Day !Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)

  REAL(KIND(1d0)),INTENT(in) ::deltaLAI
  REAL(KIND(1d0)),INTENT(in) ::VegPhenLumps
  REAL(KIND(1d0)),INTENT(in) ::SnowAlb
  REAL(KIND(1d0)),DIMENSION(7),INTENT(in)::SnowDens
  REAL(KIND(1d0)),INTENT(in) ::xBo
  REAL(KIND(1d0)),INTENT(in) ::a1
  REAL(KIND(1d0)),INTENT(in) ::a2
  REAL(KIND(1d0)),INTENT(in) ::a3

  REAL(KIND(1d0)),DIMENSION(44),INTENT(out) :: DailyStateLine


  ! Write actual data
  DailyStateLine(1:2)   = [iy,id]
  DailyStateLine(3:8)   = HDD(id,1:6)
  DailyStateLine(9:13)  = GDD(id,1:5)
  DailyStateLine(14:16) = LAI(id,1:nvegsurf)
  DailyStateLine(17:21) = [DecidCap(id),Porosity(id),AlbEveTr(id),AlbDecTr(id),AlbGrass(id)]
  DailyStateLine(22:30) = WU_day(id-1,1:9)
  DailyStateLine(31)    = deltaLAI
  DailyStateLine(32)    = VegPhenLumps
  DailyStateLine(33:40) = [SnowAlb,SnowDens(1:7)]
  DailyStateLine(41:44) = [xBo,a1,a2,a3]


END SUBROUTINE update_DailyState


!===================================================================================
!TODO: to implement a generic SUBROUTINE for DailyState output
SUBROUTINE output_DailyState(&
     Gridiv,GridIDmatrix,&!input
     FileCode,FileOutputPath,&
     DailyStateLine,&
     DailyStateFirstOpen)!inout

  IMPLICIT NONE

  INTEGER,PARAMETER:: MaxNumberOfGrids=2000   !Max no. grids   !HCW changed to 2000 from 10000 so prog can run on windows (2GB lim)


  INTEGER,INTENT(in)::Gridiv
  INTEGER,DIMENSION(MaxNumberOfGrids),INTENT(in):: GridIDmatrix         !Array containing GridIDs in SiteSelect after sorting

  CHARACTER (len = 20),INTENT(in) :: FileCode       !Set in RunControl
  CHARACTER (len = 150),INTENT(in):: FileOutputPath !Filepath for output files (set in RunControl)

  REAL(KIND(1d0)),DIMENSION(44),INTENT(in) :: DailyStateLine
  REAL(KIND(1d0)),DIMENSION(MaxNumberOfGrids),INTENT(inout):: DailyStateFirstOpen

  CHARACTER(len=10):: grstr2

  CHARACTER (len=150)::FileDaily       !Daily State output file name



  CHARACTER(len=10):: headerDaily(60)='' ! headers of daily variables



  !write(*,*) 'writing out daily state for day id:',id
  !Define filename
  ! WRITE(grstr2,'(i5)') Gridiv      !Convert grid number for output file name
  WRITE(grstr2,'(i10)') GridIDmatrix(Gridiv)      !Bug fix HCW 24/05/2016 - name file with Grid as in SiteSelect

  FileDaily=TRIM(FileOutputPath)//TRIM(FileCode)//TRIM(ADJUSTL(grstr2))//'_DailyState.txt'

  ! If first modelled day, open the file and save header
  IF (DailyStateFirstOpen(Gridiv)==1) THEN
     headerDaily(1:44)=[CHARACTER(len=10) ::'Year','DOY',&                !2
          'HDD1_h','HDD2_c','HDD3_Tmean','HDD4_T5d','P/day','DaysSR',&    !8
          'GDD1_g','GDD2_s','GDD3_Tmin','GDD4_Tmax','GDD5_DLHrs',&    !13
          'LAI_EveTr','LAI_DecTr','LAI_Grass',&                     !16
          'DecidCap','Porosity','AlbEveTr','AlbDecTr','AlbGrass',&      !21
          'WU_EveTr1','WU_EveTr2','WU_EveTr3',&               !24
          'WU_DecTr1','WU_DecTr2','WU_DecTr3',&               !27
          'WU_Grass1','WU_Grass2','WU_Grass3',&               !30
          'deltaLAI','LAIlumps','AlbSnow','DSnowPvd',&        !34
          'DSnowBldgs','DSnowEveTr','DSnowDecTr',&    !37
          'DSnowGrass','DSnowBSoil','DSnowWater',&    !40
          'BoAnOHMEnd','a1AnOHM','a2AnOHM','a3AnOHM'] !44 TS AnOHM 05 Mar 2016

     OPEN(60,file=FileDaily)
     WRITE(60,602) ADJUSTR(headerDaily(1:44))

     DailyStateFirstOpen(Gridiv)=0
     ! Otherwise open file to append
  ELSE
     OPEN(60,file=FileDaily,position='append')
     ! Write actual data
     WRITE(60,601) INT(DailyStateLine(1)),INT(DailyStateLine(2)),&
          DailyStateLine(3:44)
  ENDIF





601 FORMAT(2(i10,1X),&
       4(f10.1,1X),1(f10.4,1X),1(f10.1,1X), 5(f10.1,1X),&
       3(f10.2,1X),&
       5(f10.2,1X),&
       9(f10.3,1X),&
       2(f10.2,1X),8(f10.2,1X),&
       4(f10.2,1X))
602 FORMAT(2(a10,1X),&
       4(a10,1X),1(a10,1X),1(a10,1X), 5(a10,1X),&
       3(a10,1X),&
       5(a10,1X),&
       9(a10,1X),&
       2(a10,1X),8(a10,1X),&
       4(a10,1X))
  ! Close the daily state file
  CLOSE(60)

END SUBROUTINE output_DailyState
