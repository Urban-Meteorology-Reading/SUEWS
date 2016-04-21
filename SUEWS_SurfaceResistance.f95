SUBROUTINE SurfaceResistance(id,it)
  ! Last modified by HCW 18 Jun 2015
  ! Alternative gs calculation added using different functional forms and new coefficients
  ! Error hints added
  ! Last modified by LJ in 24 April 2013
  ! Added impact of snow fraction in LAI and in soil moisture deficit
  ! HCW 31/07/2014 Modified condition on g6 part to select meas/mod smd
  ! HCW 01/03/2016 SM dependence is now on modelled smd for vegetated surfaces only (vsmd) (Note:  obs smd still not operational!)

  USE allocateArray
  USE data_in
  USE gis_data
  USE moist
  USE resist
  USE sues_data

  IMPLICIT NONE

  INTEGER::id,it,iv
  REAL(KIND(1d0))::id_real

  !gsChoice = 1 - original Jarvis 76 and Jarvi 2011 method
  !gsChoice = 2 - new method

  id_real = REAL(id) !Day of year in real number

  IF(gsChoice == 1) THEN
     ! ResistSurf - Surface resistance --------------------------------------------
     ! At nighttime gsc set at arbitrary low value: gsc=0.1 mm/s (Shuttleworth, 1988b)
     IF(avkdn<=0)THEN  !Is this limit good here? ZZ
        gsc=0.1
     ELSE
        QNM=Kmax/(Kmax+G2)
        !gq=(qn1/(g2+qn1))/qnm !With net all-wave radiation

        gq=(avkdn/(G2+avkdn))/QNM !With Kdown

        !Specific humidity deficit
        IF(dq<G4) THEN
           gdq=1-G3*dq
        ELSE
           gdq=1-G3*G4
        ENDIF

        TC=(TH-G5)/(G5-TL)
        TC2=(G5-TL)*(TH-G5)**TC

        !If air temperature below TL or above TH,
        !fit it to TL+0.1/TH-0.1
        IF (Temp_C<=tl) THEN
           gtemp=(tl+0.1-tl)*(th-(tl+0.1))**tc/tc2

           !Call error only if no snow on ground
           IF (MIN(snowFrac(1),snowFrac(2),snowFrac(3),snowFrac(4),snowFrac(5),snowFrac(6))/=1) THEN
              CALL errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TL=0.1,Temp_c,id,it',&
                   REAL(Temp_c,KIND(1d0)),id_real,it)
           ENDIF

        ELSEIF (Temp_C>=th) THEN
           gtemp=((th-0.1)-tl)*(th-(th-0.1))**tc/tc2
           CALL errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TH=39.9,Temp_c,id,it',&
                REAL(Temp_c,KIND(1d0)),id_real,it)
        ELSE
           gtemp=(Temp_C-tl)*(th-Temp_C)**tc/tc2  ! temp
        ENDIF

        !   soil mosture deficit
        !   write(*,*)'smd:',lai(id,2),smd,xsmd

        sdp=S1/G6+S2

        IF(smd_choice>0)THEN         !Modified from ==1 to > 0 by HCW 31/07/2014
           gs=1-EXP(g6*(xsmd-sdp))  !Measured soil moisture deficit is used
        ELSE
           gs=1-EXP(g6*(vsmd-sdp))   !Modelled is used
        ENDIF

        gs = gs*(1-SUM(snowFrac(1:6))/6)

        IF(gs<0)THEN
           CALL errorHint(65,'SUEWS_SurfaceResistance (gsChoice=1): gs < 0 calculated, setting to 0.0001',gs,id_real,it)
           gs=0.0001
        ENDIF

        !LAI
        !Original way
        !gl=((lai(id,2)*areaunir/lm)+areair)/(areair+areaunir)

        !New way
        gl=0    !First initialize
        ! vegetated surfaces
        ! check basis for values koe - maximum surface conductance
        !  print*,id,it,sfr
        DO iv=ivConif,ivGrass
           !    print*,iv,laimax(iv),MaxConductance(iv),sfr(iv+2),lai(id-1,iv)
           !     pause
           ! check
           ! sg feb 2012 -- changed to use previous day LAI value
           !  and uses the maximum LAI read in
           !gl=gl+sfr(iv+2)*lai(id-1,iv)/MaxLaiMax*MaxConductance(iv)
           !gl=gl+(sfr(iv+2)*(1-snowFrac(iv+2)))*lai(id-1,iv)/MaxLaiMax*MaxConductance(iv)
           !sg 12nov13 -- removed maxLAImax - so just use the veg by class
           gl=gl+(sfr(iv+2)*(1-snowFrac(iv+2)))*lai(id-1,iv)/LaiMax(iv)*MaxConductance(iv)
        ENDDO

        gsc=(g1*gq*gdq*gtemp*gs*gl) !Original

        IF(gsc<=0) THEN
           CALL errorHint(65,'SUEWS_SurfaceResistance (gsChoice=1): gsc <= 0, setting to 0.1 mm s-1',gsc,id_real,it)
           gsc=0.1
        ENDIF

     ENDIF

  ELSEIF(gsChoice == 2) THEN
     !At nighttime, currently set gsc to arbitrary low value 0.1 mm s-1 (Shuttleworth, 1988b)   !!Address this later - HCW!!
     IF(avkdn<=0) THEN
        gsc=0.1
     ELSE   !Calculate components of surface conductance
        ! ---- g(kdown)----
        QNM = Kmax/(Kmax+G2)
        gq = (avkdn/(avkdn+G2))/QNM
        IF(avkdn >= Kmax) THEN   !! Add proper error handling later - HCW!!
           WRITE(*,*) 'Kmax exceeds Kdn setting to g(Kdn) to 1'
           gq = 1
        ENDIF
        ! ---- g(delq) ----
        gdq = G3 + (1-G3)*(G4**dq)   !Ogink-Hendriks (1995) Eq 12 (using G3 as Kshd and G4 as r)
        ! ---- g(Tair) ----
        Tc=(TH-G5)/(G5-TL)
        Tc2=(G5-TL)*(TH-G5)**Tc
        ! If air temperature below TL or above TH, then use value for TL+0.1 or TH-0.1
        IF (Temp_C <= TL) THEN
           gtemp=(TL+0.1-TL)*(TH-(TL+0.1))**Tc/Tc2
           ! Call error only if no snow on ground
           IF (MIN(snowFrac(1),snowFrac(2),snowFrac(3),snowFrac(4),snowFrac(5),snowFrac(6))/=1) THEN
              CALL errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TL+0.1,Temp_C,id,it',&
                   REAL(Temp_c,KIND(1d0)),id_real,it)
           ENDIF
        ELSEIF (Temp_C >= TH) THEN
           gtemp=((TH-0.1)-TL)*(TH-(TH-0.1))**Tc/Tc2
           CALL errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TH-0.1,Temp_C,id,it',&
                REAL(Temp_c,KIND(1d0)),id_real,it)
        ELSE
           gtemp=(Temp_C-TL)*(TH-Temp_C)**Tc/Tc2
        ENDIF
        ! ---- g(smd) ----
        sdp=S1/G6+S2
        IF(smd_choice>0) THEN                           !Modified from ==1 to > 0 by HCW 31/07/2014
           gs=(1-EXP(g6*(xsmd-sdp)))/(1-EXP(g6*(-sdp))) !Use measured smd
        ELSE
           gs=1-EXP(g6*(vsmd-sdp))   !Use modelled smd
           gs=(1-EXP(g6*(vsmd-sdp)))/(1-EXP(g6*(-sdp)))
        ENDIF

        gs = gs*(1-SUM(snowFrac(1:6))/6)

        IF(gs<0)THEN
           CALL errorHint(65,'SUEWS_SurfaceResistance (gsChoice=2): gs < 0 calculated, setting to 0.0001',gs,id_real,it)
           gs=0.0001
        ENDIF

        ! ---- g(LAI) ----
        gl=0    !Initialise
        DO iv=ivConif,ivGrass   !For vegetated surfaces
           gl=gl+(sfr(iv+2)*(1-snowFrac(iv+2)))*lai(id-1,iv)/LaiMax(iv)*MaxConductance(iv)
        ENDDO

        ! Multiply parts together
        gsc=(G1*gq*gdq*gtemp*gs*gl)

        IF(gsc<=0) THEN
           CALL errorHint(65,'SUEWS_SurfaceResistance (gsChoice=2): gsc <= 0, setting to 0.1 mm s-1',gsc,id_real,it)
           gsc=0.1
        ENDIF

     ENDIF

  ENDIF

  ResistSurf=1/(gsc/1000)  ![s m-1]

  RETURN
END SUBROUTINE SurfaceResistance
