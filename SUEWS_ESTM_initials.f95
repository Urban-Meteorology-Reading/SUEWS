SUBROUTINE ESTM_initials(FileCodeX)
  USE heatflux
  USE meteo                                                               !!FO!! :METEOMOD.f95
  USE mod_interp                                                          !!FO!! :mod_interp.f95
  USE mod_solver                                                          !!FO!! :mod_solver.f95
  USE mod_error
  USE modSolarCalc                                                        !!FO!! :modsolarcalc.f95
  USE MathConstants                                                       !!FO!! :MathConstants_module.f95
  USE ESTM_data
  USE allocateArray
  USE gis_data
  USE sues_data
  USE data_in

  IMPLICIT NONE
  INTEGER :: i,ii,ios,ios1,ios2
  INTEGER :: datalines,datalines5min
  REAL(KIND(1d0)) :: CFLval
  REAL(KIND(1d0)) :: t5min
  REAL(KIND(1d0)),ALLOCATABLE,DIMENSION(:,:)::Ts15mindata
  REAL(KIND(1d0))::W,WB
  CHARACTER (len=20),INTENT(in)::FileCodeX
  CHARACTER (len=150)::FileESTMTs,FileESTMoutput
  LOGICAL:: inittemps=.FALSE.


  !=====Reading ESTMinput.nml================================
  NAMELIST/ESTMinput/TsurfChoice,&
       evolveTibld,              &
       ibldCHmod,                &
       LBC_soil,                 &
       THEAT_ON,                 &
       THEAT_OFF,                &
       THEAT_fix

  OPEN(51,file=TRIM(FileInputPath)//'ESTMinput.nml',status='old')
  READ(51,nml=ESTMinput)
  CLOSE(51)

  dtperday=86400.0/Tstep
  ALLOCATE(Tair24HR(dtperday))
  Tair24HR=C2K

  !=======read input file===============================================================

  FileESTMTs=TRIM(FileInputPath)//TRIM(FileCodeX)//'_ESTM_Ts_data.txt'
  OPEN(10,file=TRIM(FileESTMTs),status='old',action='read',iostat=ios1)! Read the file of Ts forcing
  IF (ios1/=0) CALL error(FileESTMTs,ios1)

  READ(10,*) !skip header line
  datalines=0
  DO WHILE(ios1==0)  ! interpolate 15 min to 5min
     READ(10,*,iostat=ios1)
     datalines=datalines+1
  ENDDO
  REWIND(10)

  ALLOCATE(Ts15mindata(datalines,10))
  ALLOCATE(Ts5mindata(((datalines-1)*3+1),10))

  datalines5min=0
  READ(10,*) !skip header line
  DO i=1,datalines-1  ! interpolate 15 min to 5min
     READ(10,*) Ts15mindata(i,:)
     IF (i==1) THEN
        datalines5min=datalines5min+1
        Ts5mindata(datalines5min,:)=Ts15mindata(i,:)
     ELSE
        DO ii=1,3
           datalines5min=datalines5min+1
           t5min=Ts15mindata(i-1,1)+(Ts15mindata(i,1)-Ts15mindata(i-1,1))*ii/3
           Ts5mindata(datalines5min,:)=interp1d(Ts15mindata(i-1,1),Ts15mindata(i,1),Ts15mindata(i-1,:),Ts15mindata(i,:),t5min)
        ENDDO
     ENDIF
  ENDDO

  !=====Initialization of variables and paramters===================================

  ALLOCATE(Tibld(Nibld),Twall(Nwall),Troof(Nroof),Tground(Nground),Tw_4(Nwall,4))

  !CONVERT ALL TEMPS TO KELVIN
  DO i=1,Nground
     Tground(i)=(Ts5mindata(1,2)-Ts5mindata(1,5))*(i-1)/(Nground-1)+Ts5mindata(1,5)+C2K
  ENDDO
  DO i=1,Nwall
     Twall(i)=(Ts5mindata(1,2)-Ts5mindata(1,6))*(i-1)/(Nwall-1)+Ts5mindata(1,6)+C2K
  ENDDO
  DO i=1,Nroof
     Troof(i)=(Ts5mindata(1,2)-Ts5mindata(1,4))*(i-1)/(Nroof-1)+Ts5mindata(1,4)+C2K
  ENDDO

  Tibld(1:Nibld)=Ts5mindata(1,2)+C2K

  THEAT_fix=THEAT_fix+C2K
  Tfloor=20. ! This is used only when radforce =T
  TFLOOR=TFLOOR+C2K

  !=====Internal view factors=====================================================
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

  !=====Parameters related to land surface characteristics==========================
  ZREF=2.*BldgH                                                          !!FO!! BldgH: mean bulding hight, zref: local scale reference height (local: ~ 10^2 x 10^2 -- 10^3 x 10^3 m^2)
  svf_ground = 1.
  svf_ROOF=1.
  Tievolve = 20.+C2K;

  !roof
  froof=sfr(BldgSurf)
  alb_roof=alb(BldgSurf)
  em_roof=emis(BldgSurf)

  !vegetation
  fveg=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)
  alb_veg=(alb(ConifSurf)*sfr(ConifSurf)+alb(DecidSurf)*sfr(DecidSurf)+alb(GrassSurf)*sfr(GrassSurf))/fveg
  em_veg=(emis(ConifSurf)*sfr(ConifSurf)+emis(DecidSurf)*sfr(DecidSurf)+emis(GrassSurf)*sfr(GrassSurf))/fveg

  !ground
  fground=sfr(ConifSurf)+sfr(DecidSurf)+sfr(GrassSurf)+sfr(PavSurf)+sfr(BsoilSurf)+sfr(WaterSurf) !! S.O. This is calculated based on current version of ESTM but maybe will be changed.
  alb_ground=(alb(ConifSurf)*sfr(ConifSurf)+alb(DecidSurf)*sfr(DecidSurf)&
       +alb(GrassSurf)*sfr(GrassSurf)+alb(PavSurf)*sfr(PavSurf)&
       +alb(BsoilSurf)*sfr(BsoilSurf)+alb(WaterSurf)*sfr(WaterSurf))/fground
  em_ground=(emis(ConifSurf)*sfr(ConifSurf)+emis(DecidSurf)*sfr(DecidSurf)&
       +emis(GrassSurf)*sfr(GrassSurf)+emis(PavSurf)*sfr(PavSurf)&
       +emis(BsoilSurf)*sfr(BsoilSurf)+emis(WaterSurf)*sfr(WaterSurf))/fground

  HW=fwall/(2.*(1.-froof))

  IF (Fground==1.) THEN                                                         !!FO!! if only ground, i.e. no houses
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
  ELSE
     W=BldgH/HW
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

  alb_avg=alb_ground*RVF_ground+alb_wall*RVF_WALL+alb_roof*RVF_ROOF+alb_veg*RVF_VEG

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

  !!=======Initial setting==============================================

  IF (inittemps) THEN                                                        !!FO!! inittemps=.true. set in nml file
     OPEN(99,file='outputfiles/finaltemp.txt',status='old',iostat=ios)       !!FO!! has to exist

     IF (ios/=0) CALL error('outputfiles/finaltemp.txt',ios,1)               !!FO!! calls mod_error.f95, writes that the opening failed and stops prg
     IF (ios/=0) THEN
        Twall   = (/273., 285., 291./)
        Troof   = (/273., 285., 291./)
        Tground = (/273., 275., 280., 290./)
        Tibld   = (/293., 293., 293./)
     ELSE
        READ(99,*) Twall,Troof,Tground,Tibld                             !!FO!! if finaltemp.txt exists Twall[3], Troof[3], Tground[4] & Tibld[3] get new values
        CLOSE(99)
     ENDIF
  ENDIF
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

  DO i=1,4; Tw_4(:,i) = Twall; ENDDO                                          !!FO!! Tw_4 holds three differnet temp:s for each wall layer but the same set for all points of the compass

     !initialize surface temperatures
     T0_ground=Tground(1); T0_wall=Twall(1); T0_roof=Troof(1); T0_ibld=Tibld(1);
     TN_roof=Troof(nroof); TN_wall=Twall(nwall);

     !initialize outgoing longwave
     LUP_ground=SIGMA*EM_ground*T0_ground**4
     LUP_WALL=SIGMA*EM_WALL*T0_WALL**4
     LUP_ROOF=SIGMA*EM_ROOF*T0_ROOF**4

     !  PRINT*,"W,WB= ",W,WB
     !  PRINT*,'SVF_ground ','SVF_WALL ','zvf_WALL ','HW '
     !  PRINT*,SVF_ground,SVF_WALL,zvf_WALL,HW
     !  PRINT*,'RVF_ground ','RVF_WALL ','RVF_ROOF ','RVF_VEG'
     !  PRINT*,RVF_ground,RVF_WALL,RVF_ROOF,RVF_VEG
     !  print*,'Alb_avg (VF)=',alb_avg
     !print*,'Z0m, Zd', Z0M, ZD

     SHC_air=1230.
     minshc_airbld=1300
     first=.TRUE.
     iESTMcount=0

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
   END SUBROUTINE ESTM_initials
