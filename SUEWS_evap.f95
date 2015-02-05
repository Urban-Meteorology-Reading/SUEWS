subroutine Evap_SUEWS(StorCap)
!Input: storage capacity of each surface cover  
!Last time modified LJ 10/2010
!------------------------------------------------------------------      
  use data_in
  use moist
  use SUES_data
  use allocateArray
  use defaultNotUsed
  implicit none
  real (Kind(1d0))::StorCap,RSS,R,RSRBSG,RBSG,x,w
     
if(state(is)<=0.001) then    !Dry surface
   qe=e/(s_hPa+psyc_hPa*(1+ResistSurf/ra))
   ev=qe/tlv
   W=NAN
   rst=1

!write(*,*) 'dry', is, ev, qe, state(is), surf(1,is)

else 
  ! surface wet--------------------------------------------------------- 
   rst=0
   if(ity==2) then   !Shuttleworth (1978)      
      rbsg=rb*(sp+1) !Boundary_layer_resistance*(slope/psychom.+1)
      rsrbsg=ResistSurf+rbsg !Multiplyed by surface resistance      
      if((state(is)>=StorCap).or.(ResistSurf<25)) then !Wet surface
         w=1
      else   !if((state(is)<StorCap).and.(state(is)>0.001).or.(ResistSurf<50)) then !Surface in transition
         r=(ResistSurf/ra)*(ra-rb)/rsrbsg
	 w=(r-1)/(r-(surf(1,is)/state(is)))
      endif
      
      rss=(1/((w/rbsg)+((1-w)/rsrbsg)))-rbsg
                        
      qe=e/(s_hPa+psyc_hPa*(rss/ra+1))  !Latent heat (W/m^2)      
      ev=qe/tlv 						!Evaporation (in)

!write(*,*) 'wet', r, w
!write(*,*) 'wet', is, ev, qe, state(is), surf(1,is)          

    elseif(ity==1) then                 !  rutter              
       qe=e/(s_hPa+psyc_hPa)
       ev=qe/tlv
       if(state(is)>=StorCap) then
          x=1.0
       else
        x=state(is)/StorCap
      endif
     ev=ev*x
    	qe=ev*tlv
    endif

!write(*,*) 'wet', is, ev, qe, state(is), surf(1,is)    
endif

!write(*,*) is, ev, qe, state(is), surf(1,is)

end subroutine Evap_SUEWS