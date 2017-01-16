!======================================================================================
! Subroutine to read in ESTM data in the same way as met data (SUEWS_InitializeMetData)
! HCW 30 Jun 2016
SUBROUTINE SUEWS_GetESTMData(lunit)
    
  USE allocateArray
  USE data_in
  USE sues_data
  USE time
  USE defaultnotUsed
  USE Initial

  IMPLICIT NONE

  INTEGER::lunit,i,iyy !,RunNumber,NSHcounter
  integer::iostat_var
  REAL (KIND(1d0)),DIMENSION(ncolsESTMdata):: ESTMArray
  REAL(KIND(1d0)):: imin_prev, ih_prev, iday_prev, tstep_estm   !For checks on temporal resolution of estm data

  !---------------------------------------------------------------
 
  !Open the file for reading and read the actual data
  !write(*,*) FileESTMTs
  OPEN(lunit,file=TRIM(FileESTMTs),status='old',err=315)
  CALL skipHeader(lunit,SkipHeaderMet)

  ! Skip to the right place in the ESTM file, depending on how many chunks have been read already
  IF (skippedLines>0) THEN
     DO iyy=1,skippedLines
        READ(lunit,*)
     ENDDO
  ENDIF

  ! Read in next chunk of ESTM data and fill ESTMForcingData array with data for every timestep
  DO i=1,ReadlinesMetdata
     READ(lunit,*,iostat=iostat_var) ESTMArray
     ESTMForcingData(i,1:ncolsESTMdata,GridCounter) = ESTMArray
     ! Check timestamp of met data file matches TSTEP specified in RunControl
     IF(i==1) THEN
        imin_prev = ESTMArray(4)
        ih_prev   = ESTMArray(3)
        iday_prev = ESTMArray(2)
     ELSEIF(i==2) THEN
        tstep_estm = ((ESTMArray(4)+60*ESTMArray(3)) - (imin_prev+60*ih_prev))*60   !tstep in seconds
        IF(tstep_estm.NE.tstep_real.AND.ESTMArray(2)==iday_prev) THEN
           CALL ErrorHint(39,'TSTEP in RunControl does not match TSTEP of ESTM data (DOY).',REAL(tstep,KIND(1d0)),tstep_estm,&
                INT(ESTMArray(2)))
        ENDIF
     ENDIF
  ENDDO

  CLOSE(lunit)

  RETURN

315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,ios_out)

END SUBROUTINE SUEWS_GetESTMData
!======================================================================================

!======================================================================================
SUBROUTINE ESTM_initials
  
  ! Last modified HCW 30 Jun 2016 - reading in now done by SUEWS_GetESTMData subroutine.
  !                                 ESTM_initials now only runs once per run at the very start.
  ! Last modified HCW 15 Jun 2016 - code now reads in 5-min file (interpolation done beforehand, outside of SUEWS itself)    
    
  USE defaultNotUsed  
  USE heatflux
  USE meteo                                                               !!FO!! :METEOMOD.f95
  USE mod_interp                                                          !!FO!! :mod_interp.f95
  USE mod_solver                                                          !!FO!! :mod_solver.f95
  USE modSolarCalc                                                        !!FO!! :modsolarcalc.f95
  USE MathConstants                                                       !!FO!! :MathConstants_module.f95
  USE PhysConstants
  USE ESTM_data
  USE allocateArray
  USE gis_data
  USE sues_data
  USE data_in
  USE Initial

  IMPLICIT NONE
    
  !=====Read ESTMinput.nml================================
  NAMELIST/ESTMinput/TsurfChoice,&
       evolveTibld,              &
       ibldCHmod,                &
       LBC_soil,                 &
       THEAT_ON,                 &
       THEAT_OFF,                &
       THEAT_fix

  OPEN(511,file=TRIM(FileInputPath)//'ESTMinput.nml',status='old')
  READ(511,nml=ESTMinput)
  CLOSE(511)
 
  !Convert specified temperatures to Kelvin
  THEAT_ON=THEAT_ON+C2K
  THEAT_OFF=THEAT_OFF+C2K
  THEAT_fix=THEAT_fix+C2K
  
  ALLOCATE(Tair2_grids(NumberOfGrids))
  ALLOCATE(lup_ground_grids(NumberOfGrids))
  ALLOCATE(lup_wall_grids(NumberOfGrids))
  ALLOCATE(lup_roof_grids(NumberOfGrids))
  ALLOCATE(Tievolve_grids(NumberOfGrids))  
  ALLOCATE(T0_ibld_grids(NumberOfGrids))  
  ALLOCATE(T0_ground_grids(NumberOfGrids))  
  ALLOCATE(T0_wall_grids(NumberOfGrids))  
  ALLOCATE(T0_roof_grids(NumberOfGrids))  
  ALLOCATE(TN_wall_grids(NumberOfGrids))  
  ALLOCATE(TN_roof_grids(NumberOfGrids))  
  
 END SUBROUTINE ESTM_initials
!======================================================================================
  

 SUBROUTINE ESTM_translate(Gridiv)  
  ! HCW 30 Jun 2016 
    
  USE defaultNotUsed  
  USE heatflux
  USE meteo                                                               !!FO!! :METEOMOD.f95
  USE mod_interp                                                          !!FO!! :mod_interp.f95
  USE mod_solver                                                          !!FO!! :mod_solver.f95
  !USE mod_error
  USE modSolarCalc                                                        !!FO!! :modsolarcalc.f95
  USE MathConstants                                                       !!FO!! :MathConstants_module.f95
  USE PhysConstants
  USE ESTM_data
  USE allocateArray
  USE gis_data
  USE sues_data
  USE data_in
  USE Initial

  IMPLICIT NONE
  INTEGER :: i
  !REAL(KIND(1d0)) :: CFLval
  !REAL(KIND(1d0)) :: t5min
  REAL(KIND(1d0))::W,WB
  !CHARACTER (len=20)::FileCodeX
  !CHARACTER (len=150):: FileFinalTemp
  !LOGICAL:: inittemps=.FALSE.
  INTEGER:: ESTMStart=0
  INTEGER:: Gridiv

  !Set initial values at the start of each run for each grid
  IF(Gridiv == 1) ESTMStart = ESTMStart+1
  IF(ESTMStart==1) THEN     
    
    !write(*,*) ' ESTMStart: ',ESTMStart, 'initialising ESTM for grid no. ', Gridiv    
      
    TFLOOR=20.0 ! This is used only when radforce =T  !HCW should this be put in the namelist instead?
    TFLOOR=TFLOOR+C2K
    
    ! Initial values
    Tievolve=20.0 + C2K    
    SHC_air=1230.0
    minshc_airbld=1300
     
    ! ---- Internal view factors ----
    !constant now but should be calculated in the future
    IVF_IW =   0.100000
    IVF_IR =   0.000000
    IVF_II =   0.900000
    IVF_IF =   0.000000
    IVF_WW =   0.050000
    IVF_WR =   0.000000
    IVF_WI =   0.950000
    IVF_WF =   0.000000
    IVF_RW =   0.050000
    IVF_RI =   0.950000
    IVF_RF =   0.000000
    IVF_FW =   0.050000
    IVF_FR =   0.000000
    IVF_FI =   0.950000
   
    Tair24HR=C2K
  
    !Ts5mindata(1,ncolsESTMdata) = -999
  ! !Fill Ts5mindata for current grid and met block - this is done in SUEWS_translate
    Ts5mindata(1,1:ncolsESTMdata) = ESTMForcingData(1,1:ncolsESTMdata,Gridiv) 
 
    
    ! ---- Initialization of variables and parameters for first row of run for each grid ----
    ! N layers are calculated in SUEWS_translate
    IF ( .NOT. ALLOCATED(Tibld) ) THEN
      ! print*, "Nibld",Nibld
      ! print*, "Nwall",Nwall
      ! print*, "Nroof",Nroof
      ! print*, "Nground",Nground  
       ALLOCATE(Tibld(Nibld),Twall(Nwall),Troof(Nroof),Tground(Nground),Tw_4(Nwall,4))
       ALLOCATE(Tibld_grids(Nibld,NumberOfGrids), &
                Twall_grids(Nwall,NumberOfGrids), &
                Troof_grids(Nroof,NumberOfGrids), &
                Tground_grids(Nground,NumberOfGrids), &
                Tw_4_grids(Nwall,4,NumberOfGrids)) 
    ENDIF
    
    ! Transfer variables from Ts5mindata to variable names
    ! N.B. column numbers here for the following file format - need to change if input columns change!
    !dectime iy id it imin Tiair Tsurf Troof Troad Twall Twall_n Twall_e Twall_s Twall_w
    !        1  2  3  4    5     6     7     8     9     10      11      12      13       !new
     
    ! Calculate temperature of each layer in Kelvin
    DO i=1,Nground
       Tground(i)=(Ts5mindata(1,cTs_Tiair)-Ts5mindata(1,cTs_Troad))*(i-1)/(Nground-1)+Ts5mindata(1,cTs_Troad)+C2K   
    ENDDO
    DO i=1,Nwall
       Twall(i)=(Ts5mindata(1,cTs_Tiair)-Ts5mindata(1,cTs_Twall))*(i-1)/(Nwall-1)+Ts5mindata(1,cTs_Twall)+C2K  
    ENDDO
    DO i=1,Nroof
       Troof(i)=(Ts5mindata(1,cTs_Tiair)-Ts5mindata(1,cTs_Troof))*(i-1)/(Nroof-1)+Ts5mindata(1,cTs_Troof)+C2K
    ENDDO
    Tibld(1:Nibld)=Ts5mindata(1,cTs_Tiair)+C2K
    
  ENDIF  !End of loop run only at start (for each grid)
        
  ! ---- Parameters related to land surface characteristics ----
  ZREF=2.0*BldgH   !Would Zref=z be more appropriate?                              !!FO!! BldgH: mean bulding hight, zref: local scale reference height (local: ~ 10^2 x 10^2 -- 10^3 x 10^3 m^2)
  
  svf_ground=1.0
  svf_roof=1.0
  
  ! ==== roof (i.e. Bldgs)
  !froof=sfr(BldgSurf)   ! Moved to SUEWS_translate HCW 16 Jun 2016
  alb_roof=alb(BldgSurf)
  em_roof=emis(BldgSurf)

  ! ==== vegetation (i.e. EveTr, DecTr, Grass) 
  !fveg=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)  ! Moved to SUEWS_translate HCW 16 Jun 2016
  IF(fveg/=0) THEN
     alb_veg=(alb(ConifSurf)*sfr(ConifSurf) + alb(DecidSurf)*sfr(DecidSurf) + alb(GrassSurf)*sfr(GrassSurf))/fveg
     em_veg=(emis(ConifSurf)*sfr(ConifSurf) + emis(DecidSurf)*sfr(DecidSurf) + emis(GrassSurf)*sfr(GrassSurf))/fveg
  ENDIF 
  
  ! ==== ground (i.e. Paved, EveTr, DecTr, Grass, BSoil, Water - all except Bldgs)
  !fground=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(PavSurf)+sfr(BsoilSurf)+sfr(WaterSurf) ! Moved to SUEWS_translate HCW 16 Jun 2016
  IF(fground/=0) THEN
     alb_ground=(alb(ConifSurf)*sfr(ConifSurf)+alb(DecidSurf)*sfr(DecidSurf)&
          +alb(GrassSurf)*sfr(GrassSurf)+alb(PavSurf)*sfr(PavSurf)&
          +alb(BsoilSurf)*sfr(BsoilSurf)+alb(WaterSurf)*sfr(WaterSurf))/fground
     em_ground=(emis(ConifSurf)*sfr(ConifSurf)+emis(DecidSurf)*sfr(DecidSurf)&
          +emis(GrassSurf)*sfr(GrassSurf)+emis(PavSurf)*sfr(PavSurf)&
          +emis(BsoilSurf)*sfr(BsoilSurf)+emis(WaterSurf)*sfr(WaterSurf))/fground
  ELSE ! check fground==0 scenario to avoid division-by-zero error, TS 21 Jul 2016
     alb_ground=NAN
     em_ground=NAN
  ENDIF

  IF(froof<1.0) THEN
     HW=fwall/(2.0*(1.0-froof))
  ELSE
     HW=0  !HCW if only roof, no ground
  ENDIF

  IF (Fground==1.0) THEN   !!FO!! if only ground, i.e. no houses
     W=1
     WB=0
     SVF_ground=1.
     zvf_WALL=0.
     SVF_WALL=0.
     SVF_ROOF=1.
     zvf_ground=0.
     xvf_wall=0.
     RVF_CANYON=1.
     RVF_ground=1.-FVEG
     RVF_ROOF=0
     RVF_WALL=0
     RVF_VEG=FVEG
  ELSE IF ( Fground==0.0 ) THEN !check fground==0 (or HW==0) scenario to avoid division-by-zero error, TS 21 Jul 2016
    ! the following values are calculated given HW=0
     W=0
     WB=1
     zvf_WALL= 0 !COS(ATAN(2/HW))  when HW=0                                 !!FO!! wall view factor for wall
     HW=0
     SVF_ground=COS(ATAN(2*HW))                                              !!FO!! sky view factor for ground
     SVF_WALL=(1-zvf_WALL)/2                                                 !!FO!! sky view factor for wall
     zvf_ground=1-svf_ground                                                 !!FO!! wall view factor for ground
     xvf_wall=svf_wall                                                       !!FO!! ground view factor
     !   RVF_CANYON=COS(ATAN(2*ZREF/W))
     !   RVF_ROOF=1-RVF_CANYON
     !   RVF_WALL=(COS(ATAN(2*(ZREF-BldgH)/W))-RVF_CANYON)*RVF_CANYON
     !   RVF_ground=RVF_CANYON-RVF_WALL
     RVF_ground=(fground-fveg)*SVF_ground
     RVF_veg=fveg*SVF_ground
     RVF_ROOF=froof
     RVF_Wall=1-RVF_ROOF-RVF_ground-RVF_VEG
  ELSE
     W=BldgH/HW   !What about if HW = 0 ! need to add IF(Fground ==0) option?
     WB=W*SQRT(FROOF/Fground)
     SVF_ground=COS(ATAN(2*HW))                                              !!FO!! sky view factor for ground
     zvf_WALL=COS(ATAN(2/HW))                                                !!FO!! wall view factor for wall
     SVF_WALL=(1-zvf_WALL)/2                                                 !!FO!! sky view factor for wall
     zvf_ground=1-svf_ground                                                 !!FO!! wall view factor for ground
     xvf_wall=svf_wall                                                       !!FO!! ground view factor
     !   RVF_CANYON=COS(ATAN(2*ZREF/W))
     !   RVF_ROOF=1-RVF_CANYON
     !   RVF_WALL=(COS(ATAN(2*(ZREF-BldgH)/W))-RVF_CANYON)*RVF_CANYON
     !   RVF_ground=RVF_CANYON-RVF_WALL
     RVF_ground=(fground-fveg)*SVF_ground
     RVF_veg=fveg*SVF_ground
     RVF_ROOF=froof
     RVF_Wall=1-RVF_ROOF-RVF_ground-RVF_VEG
  ENDIF

  alb_avg=alb_ground*RVF_ground + alb_wall*RVF_WALL + alb_roof*RVF_ROOF + alb_veg*RVF_VEG

  sumalb=0.; nalb=0
  sumemis=0.; nemis=0

  !set emissivity for ceiling, wall and floor inside of buildings
  em_r = em_ibld; em_w=em_ibld; em_i=em_ibld; em_f=em_ibld

  !internal elements
  IF (nroom==0) THEN
     fibld = (FLOOR(BldgH/3.1-0.5)-1)*froof
  ELSE
     fibld = (2.-2./nroom)*fwall + (FLOOR(BldgH/3.1-0.5)-1)*froof
  ENDIF

  IF (fibld==0) fibld=0.00001 !this just ensures a solution to radiation
  finternal = froof+fibld+fwall
  fair=zref-BldgH*froof
  !ivf_ii=1.-ivf_iw-ivf_ir-ivf_if    !S.O. I do not know these are should be calculated or read from input files
  !ivf_ww=1.-ivf_wi-ivf_wr-ivf_wf
  !ivf_rw=1.-ivf_ri-ivf_rf;
  !ivf_fr=ivf_rf;

  IF ((ivf_ii+ivf_iw+ivf_ir+ivf_if > 1.0001) .OR. &
       (ivf_wi+ivf_ww+ivf_wr+ivf_wf > 1.0001) .OR. &
       (ivf_ri+ivf_rw+ivf_rf > 1.0001) .OR. &
       (ivf_fi+ivf_fw+ivf_fr > 1.0001) .OR. &
       (ivf_ii+ivf_iw+ivf_ir+ivf_if < 0.9999) .OR. &
       (ivf_wi+ivf_ww+ivf_wr+ivf_wf < 0.9999) .OR. &
       (ivf_ri+ivf_rw+ivf_rf < 0.9999) .OR. &
       (ivf_fi+ivf_fw+ivf_fr < 0.9999)) THEN
     PRINT*, "At least one internal view factor <> 1. Check ivf in ESTMinput.nml"
  ENDIF

  !!!=======Initial setting==============================================
  !! Rewritten by HCW 15 Jun 2016 to use existing SUEWS error handling
  !IF(inittemps) THEN
  !   write(*,*) 'inittemps:',inittemps
  !   FileFinalTemp=TRIM(FileOutputPath)//TRIM(FileCodeX)//'_ESTM_finaltemp.txt' 
  !   OPEN(99,file=TRIM(FileFinalTemp),status='old',err=316)  ! Program stopped if error opening file
  !   READ(99,*) Twall,Troof,Tground,Tibld                    ! Twall, Troof, Tground & Tibld get new values
  !   CLOSE(99)
  !ENDIF
  !
  !!IF (inittemps) THEN                                                        !!FO!! inittemps=.true. set in nml file
  !!   OPEN(99,file='outputfiles/finaltemp.txt',status='old',iostat=ios)       !!FO!! has to exist
  !!
  !!   IF (ios/=0) CALL error('outputfiles/finaltemp.txt',ios,1)               !!FO!! calls mod_error.f95, writes that the opening failed and stops prg
  !!   IF (ios/=0) THEN
  !!      Twall   = (/273., 285., 291./)
  !!      Troof   = (/273., 285., 291./)
  !!      Tground = (/273., 275., 280., 290./)
  !!      Tibld   = (/293., 293., 293./)
  !!   ELSE
  !!      READ(99,*) Twall,Troof,Tground,Tibld                             !!FO!! if finaltemp.txt exists Twall[3], Troof[3], Tground[4] & Tibld[3] get new values
  !!      CLOSE(99)
  !!   ENDIF
  !!ENDIF
  
  !where (isnan(Twall))
  !    Twall = 273
  !endwhere
  !where (isnan(Troof))
  !    Troof = 273
  !endwhere
  !where (isnan(Tground))
  !    Tground = 281
  !endwhere
  !where (isnan(Tibld))
  !    Tibld = 293
  !endwhere
     
  IF(ESTMStart==1) THEN 
     DO i=1,4
        Tw_4(:,i) = Twall  !!FO!! Tw_4 holds three differnet temp:s for each wall layer but the same set for all points of the compass
     ENDDO
  
     !initialize surface temperatures
     T0_ground=Tground(1)
     T0_wall=Twall(1)
     T0_roof=Troof(1)
     T0_ibld=Tibld(1)
     TN_roof=Troof(nroof)
     TN_wall=Twall(nwall)

     !initialize outgoing longwave   !HCW - Are these calculations compatible with those in LUMPS_NARP?
     LUP_ground=SBConst*EM_ground*T0_ground**4
     LUP_WALL=SBConst*EM_WALL*T0_WALL**4
     LUP_ROOF=SBConst*EM_ROOF*T0_ROOF**4
 
     !  PRINT*,"W,WB= ",W,WB
     !  PRINT*,'SVF_ground ','SVF_WALL ','zvf_WALL ','HW '
     !  PRINT*,SVF_ground,SVF_WALL,zvf_WALL,HW
     !  PRINT*,'RVF_ground ','RVF_WALL ','RVF_ROOF ','RVF_VEG'
     !  PRINT*,RVF_ground,RVF_WALL,RVF_ROOF,RVF_VEG
     !  print*,'Alb_avg (VF)=',alb_avg
     !  print*,'Z0m, Zd', Z0M, ZD

  
  ENDIF 
    
  first=.TRUE.
       
     !======Courant�Friedrichs�Lewy condition=================================
     !This is comment out by S.O. for now
     !   CFLval = minval(0.5*zibld*zibld*ribld/kibld)   !!FO!! z*z*r/k => unit [s]
     !   if (Tstep>CFLval) then !CFL condition   !!FO!! CFL condition:  Courant�Friedrichs�Lewy condition is a necessary condition for convergence while solving
     !      write(*,*) "IBLD: CFL condition: Tstep=",Tstep,">",CFLval !!FO!! certain partial differential equations numerically by the method of finite differences (like eq 5 in Offerle et al.,2005)
     !      CFLfail=.TRUE.
     !   endif
     !   CFLval = minval(0.5*zroof*zroof*rroof/kroof)
     !   if (Tstep>CFLval) then !CFL condition
     !      write(*,*) "ROOF: CFL condition: Tstep=",Tstep,">",CFLval
     !      CFLfail=.TRUE.
     !   endif
     !   CFLval = minval(0.5*zwall*zwall*rwall/kwall)
     !   if (Tstep>CFLval) then !CFL condition
     !      write(*,*) "WALL: CFL condition: Tstep=",Tstep,">",CFLval
     !      CFLfail=.TRUE.
     !   endif
     !   CFLval = minval(0.5*zground*zground*rground/kground)
     !   if (Tstep>CFLval) then !CFL condition
     !      write(*,*) "ground: CFL condition: Tstep=",Tstep,">",CFLval
     !      CFLfail=.TRUE.
     !   endif
     !   if (CFLfail) then
     !      write(*,*) "Increase dX or decrease maxtimestep. Hit any key to continue"
     !      read (*,*)
     !   endif


     ! Tiaircyc = (1+(LondonQSJune_Barbican.Tair-Tiair)./(5*Tiair)).*(Tiair + 0.4*sin(LondonQSJune_Barbican.HOUR*2*pi/24-10/24*2*pi))    !!FO!! outdoor temp affected
      
    
      
     RETURN 
     
!     315 CALL errorHint(11,TRIM(fileESTMTs),notUsed,notUsed,NotUsedI)
!     316 CALL errorHint(11,TRIM(fileFinalTemp),notUsed,notUsed,NotUsedI)
     
 END SUBROUTINE ESTM_translate
