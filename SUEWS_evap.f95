SUBROUTINE Evap_SUEWS
  !------------------------------------------------------------------------------
  !-Calculates evaporation for each surface from modified Penman-Monteith eqn
  !-State determines whether each surface type is dry or wet (wet/transition)
  !-Wet surfaces below storage capacity are in transition
  ! and QE depends on the state and storage capacity (i.e. varies with surface);
  ! for wet or dry surfaces QE does not vary between surface types
  !-See Sect 2.4 of Jarvi et al. (2011) Ja11
  !
  !Last modified:
  !  HCW 06 Jul 2016
  !   Moved rss declaration to LUMPS_Module_Constants so it can be written out
  !  HCW 11 Jun 2015
  !   Added WetThresh to distinguish wet/partially wet surfaces from the storage capacities used in SUEWS_drain
  !  HCW 30 Jan 2015
  !   Removed StorCap input because it is provided by module allocateArray
  !   Tidied and commented code
  !  LJ 10/2010
  !------------------------------------------------------------------------------

  USE allocateArray
  USE data_in
  USE defaultNotUsed

  USE moist
  USE sues_data

  IMPLICIT NONE

  REAL(KIND(1d0)):: rbsg,&  !Boundary-layer resistance x (slope/psychrometric const + 1) [s m-1]
       rsrbsg,&  !rs + rbsg [s m-1]
       W,&  !Depends on the amount of water on the canopy [-]
       r,&
       x

  ! Use Penman-Monteith eqn modified for urban areas (Eq6, Jarvi et al. 2011)
  ! Calculation independent of surface characteristics
  ! Uses value of rs for whole area (calculated based on LAI of veg surfaces in SUEWS_SurfaceResistance)

  ! PRINT*, 'is',is,'SMOIS',state(is)
  ! Dry surface ---------------------------------------------------------------
  IF(state(is)<=0.001) THEN
     qe  = numPM/(s_hPa+psyc_hPa*(1+ResistSurf/ra))  !QE [W m-2] (numPM = numerator of P-M eqn)
     ev  = qe/tlv              !Ev [mm]    (qe[W m-2]/tlv[J kg-1 s-1]*1/density_water[1000 kg m-3])
     W   = NAN                  !W not needed for dry surfaces (set to -999)
     rst = 1                  !Set flag indicating dry surface(1)

     ! Wet surface ---------------------------------------------------------------
  ELSE
     rst=0   !Set flag=0 indicating wet surface(0)

     ! Evaporation calculated according to Rutter(ity=1) or Shuttleworth(ity=2).
     !Set in SUEWS_initial (so not an input to the model)
     IF(ity==2) THEN   !-- Shuttleworth (1978) --
        rbsg   = rb*(sp+1)           !Boundary-layer resistance x (slope/psychro + 1)
        rsrbsg = ResistSurf+rbsg   !rs + rsbg

        ! If surface is completely wet, set rs to zero -------------------
        !if(state(is)>=surf(6,is).or.ResistSurf<25) then   !If at storage capacity or rs is small
        IF(state(is)>=WetThresh(is).OR.ResistSurf<25) THEN   !If at storage capacity or rs is small
           W=1                                            !So that rs=0 (Eq7, Jarvi et al. 2011)
           ! If surface is in transition, use rss ---------------------------
        ELSE   !if((state(is)<StorCap).and.(state(is)>0.001).or.(ResistSurf<50)) then
           r = (ResistSurf/ra)*(ra-rb)/rsrbsg
           W = (r-1)/(r-(WetThresh(is)/state(is)))
        ENDIF

        ! PRINT*, 'r',r
        ! PRINT*, 'W',W

        rss=(1/((W/rbsg)+((1-W)/rsrbsg)))-rbsg !Redefined surface resistance for wet
        ! PRINT*, 'resistances:',rbsg,rsrbsg,rss
        !surfaces (zero if W=1). Eq7, Jarvi et al. (2011)
        qe = numPM/(s_hPa+psyc_hPa*(1+rss/ra))   !QE [W m-2]
        ev = qe/tlv                              !Ev [mm]
        ! PRINT*, 'numPM',numPM
        ! PRINT*, 'qe',qe

     ELSEIF(ity==1) THEN   !-- Rutter --
        qe = numPM/(s_hPa+psyc_hPa)
        ev = qe/tlv

        IF(state(is) >= surf(6,is)) THEN
           x = 1.0
        ELSE
           x = state(is)/surf(6,is)
        ENDIF
        ev = ev*x     !QE [W m-2]
        qe = ev*tlv   !Ev [mm]
     ENDIF   !Rutter/Shuttleworth calculation
  ENDIF   !Wet/dry surface



  ! IF ( id>190 ) THEN
  !    STOP "stop in Evap_SUEWS"
  !
  ! END IF
END SUBROUTINE Evap_SUEWS
