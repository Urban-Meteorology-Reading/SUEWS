subroutine Evap_SUEWS
!------------------------------------------------------------------------------
!Calcualates evaporation for each surface from modified Penman-Monteith eqn
!State determines whether each surface type is dry or wet (wet/transition)
!Wet surfaces below storage capacity are in transition 
! and QE depends on the state and storage capacity (i.e. varies with surface);
! for wet or dry surfaces QE does not vary between surface types
!See Sect 2.4 of Jarvi et al. (2011) Ja11
!
!Last modified HCW 30 Jan 2015
! Removed StorCap input because it is provided by module allocateArray
! Tidied and commented code
!Last modified LJ 10/2010
!------------------------------------------------------------------------------      
  
  use allocateArray
  use data_in
  use defaultNotUsed
  use moist
  use sues_data
    
  IMPLICIT NONE
  
  real(kind(1d0)):: rss,&	!Redefined surface resistance for transition [s m-1]
  		    rbsg,&	!Boundary-layer resistance x (slope/psychrometric const + 1) [s m-1]
  		    rsrbsg,&	!rs + rbsg [s m-1]
  		    W,&		!Depends on the amount of water on the canopy [-]
  		    r,x  
        
  ! Dry surface ---------------------------------------------------------------
  ! Use Penman-Monteith eqn modified for urban areas (Eq6, Jarvi et al. 2011)      
  ! Calculation independent of surface characteristics
  ! Uses value of rs for whole area (calculated based on LAI of veg surfaces in SUEWS_SurfaceResistance)
  if(state(is)<=0.001) then    
     qe=numPM/(s_hPa+psyc_hPa*(1+ResistSurf/ra))	!QE [W m-2] (numPM = numerator of P-M eqn)
     ev=qe/tlv					!Ev [mm]    (qe[W m-2]/tlv[J kg-1 s-1]*1/density_water[1000 kg m-3])	
     W=NAN   !W not needed for dry surfaces (set to -999)
     rst=1   !Set flag indicating dry surface(1)	
  
  else 
  ! Wet surface ---------------------------------------------------------------
     rst=0   !Set flag indicating wet surface(0)	
     
     ! Evaporation calculated according to Rutter(ity=1) or Shuttleworth(ity=2)
  
     if(ity==2) then   !-- Shuttleworth (1978) --
        rbsg=rb*(sp+1)           !Boundary-layer resistance x (slope/psychro + 1)
        rsrbsg=ResistSurf+rbsg   !rs + rsbg
        ! If surface is completely wet, set rs to zero -------------------
        if(state(is)>=surf(6,is).or.ResistSurf<25) then   !If at storage capacity or rs is small
           W=1   !So that rs=0 (Eq7, Jarvi et al. 2011)
        ! If surface is in transition, use rss ---------------------------
        else   !if((state(is)<StorCap).and.(state(is)>0.001).or.(ResistSurf<50)) then
           r=(ResistSurf/ra)*(ra-rb)/rsrbsg
	   W=(r-1)/(r-(surf(6,is)/state(is)))
        endif
        ! Calculate redefined surface resistance for wet surfaces (zero if W=1)
        ! Eq7, Jarvi et al. 2011
        rss=(1/((W/rbsg)+((1-W)/rsrbsg)))-rbsg   
        qe=numPM/(s_hPa+psyc_hPa*(1+rss/ra))   !QE [W m-2]
        ev=qe/tlv 				 !Ev [mm]		

     elseif(ity==1) then   !-- Rutter --
        qe=numPM/(s_hPa+psyc_hPa)
        ev=qe/tlv
        if(state(is)>=surf(6,is)) then
           x=1.0
        else
           x=state(is)/surf(6,is)
        endif
        ev=ev*x	    !QE [W m-2]
        qe=ev*tlv   !Ev [mm]
     endif   !Rutter/Shuttleworth calculation
  endif   !Wet/dry surface

end subroutine Evap_SUEWS
!------------------------------------------------------------------------------