! Conversion of water use (irrigation)
! Last modified:
!  TS 08 Aug 2017  - addded explicit interface
!  LJ  6 Apr 2017  - WUchoice changed to WaterUseMethod
!  TK 14 Mar 2017  - Corrected the variable name WUAreaEveTr_m2 -> WUAreaGrass_m2 (row 35)
!                    Corrected conversion from m to mm /1000 -> *1000 (row 47 and 60)
!  LJ 27 Jan 2016  - Removing Tab:s and cleaning the code
!  HCW 12 Feb 2015 - Water use [mm] now inidcates the amount of water supplied for each surface
!  HCW 26 Jan 2015 - Water use [mm] is the same for each surface at the moment and indicates the
!                    amount of water supplied for each irrigated area
!
! To Do:
!	- Add functionality for water on paved surfaces (street cleaning, fountains)
!===================================================================================
SUBROUTINE WaterUse(&
     ! input:
     nsh_real,&
     SurfaceArea,&
     sfr,&
     IrrFracConif,&
     IrrFracDecid,&
     IrrFracGrass,&
     DayofWeek_id,&
     WUProfA_tstep,&
     WUProfM_tstep,&
     InternalWaterUse_h,&
     HDD_id,&
     WU_Day_id,&
     WaterUseMethod,&
     ConifSurf,&
     DecidSurf,&
     GrassSurf,&
     NSH,&
     it,imin,DLS,nsurf,&
     OverUse,&
     !  output:
     WUAreaEveTr_m2,&
     WUAreaDecTr_m2,&
     WUAreaGrass_m2,&
     WUAreaTotal_m2,&
     wu_EveTr,&
     wu_DecTr,&
     wu_Grass,&
     wu_m3,&
     int_wu,&
     ext_wu)

  ! USE allocateArray
  ! USE data_in
  ! USE defaultNotUsed
  ! USE sues_data
  ! USE time

  IMPLICIT NONE

  REAL(KIND(1d0)),INTENT(in):: &
       nsh_real,&
       SurfaceArea,& !Surface area of the study area [m2]
       sfr(nsurf),& !Surface fractions [-]
       IrrFracConif,&!Fraction of evergreen trees which are irrigated
       IrrFracDecid,&!Fraction of deciduous trees which are irrigated
       IrrFracGrass,&!Fraction of grass which is irrigated
       DayofWeek_id(3),& !DayofWeek(id) 1 - day of week; 2 - month; 3 - season
       WUProfA_tstep(24*NSH,2),& !Automatic water use profiles at model timestep
       WUProfM_tstep(24*NSH,2),& !Manual water use profiles at model timestep
       InternalWaterUse_h,& !Internal water use [mm h-1]
       HDD_id(6),& !HDD(id-1), Heating Degree Days (see SUEWS_DailyState.f95)
       WU_Day_id(9) !WU_Day(id-1), Daily water use for EveTr, DecTr, Grass [mm] (see SUEWS_DailyState.f95)

  INTEGER,INTENT(in):: &
       WaterUseMethod,& !Use modelled (0) or observed (1) water use
       ConifSurf,& !surface code
       DecidSurf,& !surface code
       GrassSurf,& !surface code
       NSH,&!Number of timesteps per hour
       it,& !Hour
       imin,& !Minutes
       DLS,& !day lightsavings =1 + 1h) =0
       nsurf

  REAL(KIND(1d0)),INTENT(inout):: &
       OverUse

  REAL(KIND(1d0)),INTENT(out):: &
       WUAreaEveTr_m2,&
       WUAreaDecTr_m2,&
       WUAreaGrass_m2,&
       WUAreaTotal_m2,&
       wu_EveTr,&
       wu_DecTr,&
       wu_Grass,&
       wu_m3,&
       int_wu,&
       ext_wu


  REAL(KIND(1d0)):: &
       InternalWaterUse,&    !Internal water use for the model timestep [mm]
       WuFr=1,&
       wu!Water use for the model timestep [mm]
  INTEGER:: ih   !Hour corrected for Daylight savings
  INTEGER:: iu   !1=weekday OR 2=weekend
  REAL(KIND(1d0)),PARAMETER::NAN=-999.

  ! --------------------------------------------------------------------------------
  ! If water used is observed and provided in the met forcing file, units are m3
  ! Divide observed water use (in m3) by water use area to find water use (in mm)
  IF (WaterUseMethod==1) THEN   !If water use is observed
     ! Calculate water use area [m2] for each surface type
     WUAreaEveTr_m2 = IrrFracConif*sfr(ConifSurf)*SurfaceArea
     WUAreaDecTr_m2 = IrrFracDecid*sfr(DecidSurf)*SurfaceArea
     WUAreaGrass_m2 = IrrFracGrass*sfr(GrassSurf)*SurfaceArea
     WUAreaTotal_m2 = WUAreaEveTr_m2 + WUAreaDecTr_m2 + WUAreaGrass_m2

     !Set water use [mm] for each surface type to zero initially
     wu_EveTr=0
     wu_DecTr=0
     wu_Grass=0
     IF(wu_m3==NAN.OR.wu_m3==0) THEN !If no water use
        wu_m3=0
        wu=wu_m3
     ELSE                            !If water use
        IF (WUAreaTotal_m2>0) THEN
           wu = (wu_m3/WUAreaTotal_m2*1000)  !Water use in mm for the whole irrigated area
           IF (WUAreaEveTr_m2>0) THEN
              wu_EveTr=wu                    !Water use for Irr EveTr in mm - these are all the same at the moment
              wu_EveTr=wu_EveTr*IrrFracConif !Water use for EveTr in mm
           ENDIF
           IF (WUAreaDecTr_m2>0) THEN
              wu_DecTr=wu                        !Water use for Irr DecTr in mm - these are all the same at the moment
              wu_DecTr=wu_DecTr*IrrFracDecid     !Water use for DecTr in mm
           ENDIF
           IF (WUAreaGrass_m2>0) THEN
              wu_Grass=wu                    !Water use for Irr Grass in mm - these are all the same at the moment
              wu_Grass=wu_Grass*IrrFracGrass !Water use for Grass in mm
           ENDIF
           wu = (wu_m3/SurfaceArea*1000)     !Water use for the whole study area in mm
        ENDIF
     ENDIF

     ! --------------------------------------------------------------------------------
     ! If water use is modelled, calculate at timestep of model resolution [mm]
  ELSEIF (WaterUseMethod==0) THEN   !If water use is modelled

     ! Account for Daylight saving
     ih=it-DLS
     IF (ih<0) ih=23

     ! Weekday or weekend profile
     iu=1     !Set to 1=weekday
     !  IF(DayofWeek(id,1)==1.OR.DayofWeek(id,1)==7) THEN
     IF(DayofWeek_id(1)==1.OR.DayofWeek_id(1)==7) THEN
        iu=2  !Set to 2=weekend
     ENDIF

     !write(*,*) (NSH*(ih+1-1)+imin*NSH/60+1)

     ! ---- Automatic irrigation ----
     wu_EveTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(2)   !Automatic evergreen trees
     wu_DecTr = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(5)   !Automatic deciduous trees
     wu_Grass = WUProfA_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(8)   !Automatic grass

     ! ---- Manual irrigation ----
     WuFr=1 !Initialize WuFr to 1, but if raining, reduce manual fraction of water use
     ! If cumulative daily precipitation exceeds 2 mm
     IF(HDD_id(5)>2) THEN    !.and.WU_Day(id-1,3)>0) then !Commented out HCW 23/01/2015
        WuFr=0   ! 0 -> No manual irrigation if raining
     ENDIF

     ! Add manual to automatic to find total irrigation
     wu_EveTr = wu_EveTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(3)) !Manual evergreen trees
     wu_DecTr = wu_DecTr + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(6)) !Manual deciduous trees
     wu_Grass = wu_Grass + (WuFr*WUProfM_tstep((NSH*(ih+1-1)+imin*NSH/60+1),iu)*WU_Day_id(9)) !Manual grass

     ! Added HCW 12 Feb 2015.
     !wu_EveTr=wu_EveTr*sfr(ConifSurf)*IrrFracConif	!Water use for EveTr [mm]
     !wu_DecTr=wu_DecTr*sfr(DecidSurf)*IrrFracDecid	!Water use for DecTr [mm]
     !wu_Grass=wu_Grass*sfr(GrassSurf)*IrrFracGrass	!Water use for Grass [mm]
     wu_EveTr=wu_EveTr*IrrFracConif  !Water use for EveTr [mm]
     wu_DecTr=wu_DecTr*IrrFracDecid  !Water use for DecTr [mm]
     wu_Grass=wu_Grass*IrrFracGrass  !Water use for Grass [mm]

     ! Total water use for the whole study area [mm]
     wu = wu_EveTr*sfr(ConifSurf) + wu_DecTr*sfr(DecidSurf) + wu_Grass*sfr(GrassSurf)

  ENDIF   !End WU_choice
  ! --------------------------------------------------------------------------------

  ! Internal water use is supplied in SUEWS_Irrigation in mm h-1
  ! Convert to mm for the model timestep
  InternalWaterUse = InternalWaterUse_h/nsh_real

  ! Remove InternalWaterUse from the total water use
  ext_wu = wu-(InternalWaterUse+OverUse)
  ! Check ext_wu cannot be negative
  IF (ext_wu<0) THEN
     overUse=ABS(ext_wu)
     ext_wu=0
  ELSE
     OverUse=0
  ENDIF

  int_wu = wu-ext_wu

  ! Decrease the water use for each surface by the same proportion
  IF(ext_wu/=0.AND.wu/=0) THEN
     wu_EveTr = wu_EveTr*ext_wu/wu
     wu_DecTr = wu_DecTr*ext_wu/wu
     wu_Grass = wu_Grass*ext_wu/wu
  ENDIF

endsubroutine WaterUse
!===================================================================================
