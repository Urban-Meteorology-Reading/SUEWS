module evap_module
   implicit none

contains
   SUBROUTINE Evap_SUEWS( &
      EvapMethod, state_is, WetThresh_is, capStore_is, &!input
      vpd_hPa, avdens, avcp, qn_e, s_hPa, psyc_hPa, ResistSurf, RA, rb, tlv, &
      rss, ev, qe) !output
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

      IMPLICIT NONE
      INTEGER, INTENT(in) :: EvapMethod!Evaporation calculated according to Rutter (1) or Shuttleworth (2)

      REAL(KIND(1d0)), INTENT(in)::state_is ! wetness status
      REAL(KIND(1d0)), INTENT(in)::WetThresh_is!When State > WetThresh, RS=0 limit in SUEWS_evap [mm] (specified in input files)
      REAL(KIND(1d0)), INTENT(in)::capStore_is ! = StoreDrainPrm(6,is), current storage capacity [mm]

      REAL(KIND(1d0)), INTENT(in)::vpd_hPa ! vapour pressure deficit [hPa]
      REAL(KIND(1d0)), INTENT(in)::avdens ! air density
      REAL(KIND(1d0)), INTENT(in)::avcp ! air heat capacity
      REAL(KIND(1d0)), INTENT(in)::qn_e !net available energy for evaporation
      REAL(KIND(1d0)), INTENT(in)::s_hPa!Vapour pressure versus temperature slope in hPa
      REAL(KIND(1d0)), INTENT(in)::psyc_hPa!Psychometric constant in hPa
      REAL(KIND(1d0)), INTENT(in)::ResistSurf!Surface resistance
      ! REAL(KIND(1d0)),INTENT(in)::sp!Term in calculation of E
      REAL(KIND(1d0)), INTENT(in)::RA!Aerodynamic resistance
      REAL(KIND(1d0)), INTENT(in)::rb!Boundary layer resistance
      REAL(KIND(1d0)), INTENT(in)::tlv!Latent heat of vaporization per timestep [J kg-1 s-1], (tlv=lv_J_kg/tstep_real)

      REAL(KIND(1d0)), INTENT(out)::rss !Redefined surface resistance for wet
      REAL(KIND(1d0)), INTENT(out)::ev ! evapotranspiration [mm]
      REAL(KIND(1d0)), INTENT(out)::qe ! latent heat flux [W m-2]

      REAL(KIND(1d0))::numPM!numerator of P-M eqn
      REAL(KIND(1d0))::rbsg  !Boundary-layer resistance x (slope/psychrometric const + 1) [s m-1]
      REAL(KIND(1d0))::rsrbsg  !RS + rbsg [s m-1]
      REAL(KIND(1d0))::rst
      REAL(KIND(1d0))::W  !Depends on the amount of water on the canopy [-]
      REAL(KIND(1d0))::x
      REAL(KIND(1d0))::r

      REAL(KIND(1d0)), PARAMETER::  NAN = -999

      ! Use Penman-Monteith eqn modified for urban areas (Eq6, Jarvi et al. 2011)
      ! Calculation independent of surface characteristics
      ! Uses value of RS for whole area (calculated based on LAI of veg surfaces in SUEWS_SurfaceResistance.f95)

      ! PRINT*, 'is',is,'SMOIS',state(is)
      ! PRINT*, 'SMOIS',state_is,state_is<=0.001
      ! PRINT*, 'EvapMethod',EvapMethod

      !numerator of P-M eqn, refer to Eq6, Jarvi et al. 2011
      numPM = s_hPa*qn_e + vpd_hPa*avdens*avcp/RA !s_haPa - slope of svp vs t curve.

      ! Dry surface ---------------------------------------------------------------
      IF (state_is <= 0.001) THEN
         qe = numPM/(s_hPa + psyc_hPa*(1 + ResistSurf/RA))  !QE [W m-2] (numPM = numerator of P-M eqn)
         ev = qe/tlv !Ev [mm] (qe[W m-2]/tlv[J kg-1 s-1]*1/density_water[1000 kg m-3])
         W = NAN    !W not needed for dry surfaces (set to -999)
         rst = 1      !Set flag indicating dry surface(1)

         ! Wet surface ---------------------------------------------------------------
      ELSE
         rst = 0   !Set flag=0 indicating wet surface(0)

         ! Evaporation calculated according to Rutter(EvapMethod=1) or Shuttleworth(EvapMethod=2).
         !Set in SUEWS_initial (so not an input to the model)
         IF (EvapMethod == 2) THEN   !-- Shuttleworth (1978) --
            rbsg = rb*(s_hPa/psyc_hPa + 1)           !Boundary-layer resistance x (slope/psychro + 1)
            rsrbsg = ResistSurf + rbsg   !RS + rsbg

            ! If surface is completely wet, set RS to zero -------------------
            !if(state(is)>=StoreDrainPrm(6,is).or.ResistSurf<25) then   !If at storage capacity or RS is small
            IF (state_is >= WetThresh_is .OR. ResistSurf < 25) THEN   !If at storage capacity or RS is small
               W = 1                                            !So that RS=0 (Eq7, Jarvi et al. 2011)
               ! If surface is in transition, use rss ---------------------------
            ELSE   !if((state(is)<StorCap).and.(state(is)>0.001).or.(ResistSurf<50)) then
               r = (ResistSurf/RA)*(RA - rb)/rsrbsg
               W = (r - 1)/(r - (WetThresh_is/state_is))
            ENDIF

            ! PRINT*, 'r',r
            ! PRINT*, 'W',W

            rss = (1/((W/rbsg) + ((1 - W)/rsrbsg))) - rbsg !Redefined surface resistance for wet
            ! PRINT*, 'resistances:',rbsg,rsrbsg,rss
            !surfaces (zero if W=1). Eq7, Jarvi et al. (2011)
            qe = numPM/(s_hPa + psyc_hPa*(1 + rss/RA))   !QE [W m-2]
            ev = qe/tlv                              !Ev [mm]
            ! PRINT*, 'numPM',numPM
            ! PRINT*, 'qe',qe

         ELSEIF (EvapMethod == 1) THEN   !-- Rutter --
            qe = numPM/(s_hPa + psyc_hPa)
            ev = qe/tlv

            x = MERGE(1d0, state_is/capStore_is, state_is > capStore_is)
            ev = ev*x !QE [W m-2]
            qe = ev*tlv !Ev [mm]
         ENDIF   !Rutter/Shuttleworth calculation
      ENDIF   !Wet/dry surface

      ! IF ( id>190 ) THEN
      !    STOP "stop in Evap_SUEWS_new"
      !
      ! END IF
   ENDSUBROUTINE Evap_SUEWS

end module evap_module
