subroutine SurfaceResistance(id,it)
!Last modified by LJ in 24 April 2013
!Added impact of snow fraction in LAI and in soil moisture deficit


  use sues_data
  use moist
  use data_in
  use resist        ! lumps_MetRead.f95
  use gis_data      ! LUMPS_gis_read.f90  
  use allocateArray

  implicit none

  integer::id,it,iv
  real(kind(1d0))::id_real

  id_real = real(id) !Day of year in real number

  ! ResistSurf - Surface resistance --------------------------------------------
  ! At nighttime gsc set at arbitrary low value: gsc=0.1 mm/s (Shuttleworth, 1988b)
  if(avkdn<=0)then  !Is this limit good here? ZZ 
      gsc=0.1
  else 
     
     QNM=Kmax/(Kmax+G2)
    !gq=(qn1/(g2+qn1))/qnm !With net all-wave radiation
      
     gq=(avkdn/(G2+avkdn))/QNM !With Kdown
      
     !Specific humidity deficit
     if(dq<G4) then
        gdq=1-G3*dq
     else
        gdq=1-G3*G4
     endif
  
     TC=(TH-G5)/(G5-TL)     
     TC2=(G5-TL)*(TH-G5)**TC    
      
     !If air temperature below TL or above TH,
     !fit it to TL+0.1/TH-0.1
     if (Temp_C<=tl) then
        gtemp=(tl+0.1-tl)*(th-(tl+0.1))**tc/tc2
  
        !Call error only if no snow on ground
        if (min(snowFrac(1),snowFrac(2),snowFrac(3),snowFrac(4),snowFrac(5),snowFrac(6))/=1) then
           call errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TL=0.1,Temp_c,id,it',&
                real(Temp_c,kind(1d0)),id_real,it)
        endif
     
     elseif (Temp_C>=th) then
        gtemp=((th-0.1)-tl)*(th-(th-0.1))**tc/tc2
        call errorHint(29,'subroutine SurfaceResistance:T changed to fit limits TH=39.9,Temp_c,id,it',&
             real(Temp_c,kind(1d0)),id_real,it)
     else
        gtemp=(Temp_C-tl)*(th-Temp_C)**tc/tc2  ! temp
     endif
      
  !   soil mosture deficit
  !   write(*,*)'smd:',lai(id,2),smd,xsmd
   
      sdp=S1/G6+S2
      
      if(smd_choice==1)then
          gs=1-exp(g6*(xsmd-sdp))!Measured soil moisture deficit is used
      else
          gs=1-exp(g6*(smd-sdp))!Modelled is used
      endif

      
      gs = gs*(1-sum(snowFrac(1:6))/6)
 
   
      if(gs<0)then
          gs=0.0001
      endif
   
      !LAI
      !Original way
      !gl=((lai(id,2)*areaunir/lm)+areair)/(areair+areaunir) 
      
      !New way
      gl=0    !First initialize 
      ! vegetated surfaces
      ! check basis for values koe - maximum surface conductance
    !  print*,id,it,sfr
      do iv=ivConif,ivGrassU
    !    print*,iv,laimax(iv),MaxConductance(iv),sfr(iv+2),lai(id-1,iv)
    !     pause
          ! check 
          ! sg feb 2012 -- changed to use previous day LAI value
          !  and uses the maximum LAI read in
          !gl=gl+sfr(iv+2)*lai(id-1,iv)/MaxLaiMax*MaxConductance(iv)
        !!  gl=gl+(sfr(iv+2)*(1-snowFrac(iv+2)))*lai(id-1,iv)/MaxLaiMax*MaxConductance(iv)
        !! sg 12nov13 -- removed maxLAImax - so just use the veg by class
             gl=gl+(sfr(iv+2)*(1-snowFrac(iv+2)))*lai(id-1,iv)/LaiMax(iv)*MaxConductance(iv)
          
      enddo
      
      gsc=(g1*gq*gdq*gtemp*gs*gl) !Original
      
      IF(gsc<=0) gsc=0.1 
     
 endif
  
 ResistSurf=1/(gsc/1000)

  return
end subroutine SurfaceResistance
