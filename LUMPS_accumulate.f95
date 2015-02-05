!Subroutines included:
!	subroutine accumulate(id,total)
!	subroutine out_accumulate(lim,total,unit)
!
subroutine accumulate(id,total)
!Last modified by LJ 9/2010
!sg feb 2012 -- allocate arrays each time
  use data_in
  use gis_data
  use SUES_data
  use allocateArray
  use snowMod
  
  implicit none
  
  integer:: id
  real (kind(1d0)),dimension(ndays,NumberDailyVars)::total
  
  ! check why is this 390 rather than 366?
  
 
  total(id,1)=id
  total(id,2)=total(id,2)+1 ! counter
  total(id,3)=total(id,3)+qn1
  total(id,4)=total(id,4)+qs
  total(id,5)=total(id,5)+qf
  total(id,6)=total(id,6)+qeph
  total(id,7)=total(id,7)+Precip
  total(id,8)=total(id,8)+ext_wu
  total(id,9)=total(id,9)+int_wu
  total(id,10)=total(id,10)+int_wu+ext_wu
  total(id,11)=total(id,11)+ev_per_interval
  total(id,12)=total(id,12)+ch_per_interval
 !total(id,13)=soilmoist_state
  total(id,13)=total(id,13)+runoffSoil_per_interval
  total(id,14)=total(id,14)+runoff_per_interval
  total(id,15)=total(id,15)+FlowChange*sfr(WaterSurf)
  total(id,16)=total(id,16)+AdditionalWater
  total(id,17)=total(id,17)+QH  !SUEWS sensible heat flux
  
  total(id,18)=total(id,18)+Qm  !Heat related to snow melt
  total(id,19)=total(id,19)+QmFreez !Freezing energy change
  total(id,20)=total(id,20)+QmRain  !Internal energy change
  total(id,21)=total(id,21)+swe !snow water equivalent 
  total(id,22)=total(id,22)+MwStore !Meltwater store
  total(id,23)=total(id,23)+SnowRemoval(1) !Snow removal
  total(id,24)=total(id,24)+SnowRemoval(2) !Snow removal
  total(id,25)=total(id,25)+chSnow_per_interval

  return
end subroutine accumulate
! ==================================================================
subroutine out_accumulate(lim,total,unit,sel)
 use data_in
 use SUES_data
 use  allocateArray
 
 implicit none

 integer:: lim,unit,i,j,sel
 ! check why is this 390 rather than 366?
  real (kind(1d0)),dimension(366,NumberDailyVars)::total

 do i=1,lim 
  if(total(i,2)>0)then
    if(sel>0)then
	  total(i,1)=sel
	endif
      
   total(i,3)=total(i,3)/total(i,2)  
   total(i,4)=total(i,4)/total(i,2)
   total(i,5)=total(i,5)/total(i,2)
   total(i,6)=total(i,6)/total(i,2)
   total(i,17)=total(i,17)/total(i,2) !QH
   total(i,18)=total(i,18)/total(i,2) !QM
   total(i,19)=total(i,19)/total(i,2) !delta_Qsi
   total(i,20)=total(i,20)/total(i,2) !QP
   total(i,21)=total(i,21)/total(i,2) !snowpack
    
    !Write data to daily files.
    write(unit,140)(total(i,j),j=1,NumberDailyVars) !When snow allowed
    !write(unit,140)(total(i,j),j=1,17)
    
    if(unit==14.and.allDays==1.and.total(i,1)>0)then	
    	write(90,140)year, total(i,1),(total(i,j),j=3,NumberDailyVars)
        !write(90,140)year, total(i,1),(total(i,j),j=3,17)
    endif
    
 
   endif 
enddo
!140 format(F5.0,F6.0,17f13.5)
140 format((F5.0,1X),(F5.0,1X),4(f8.2,1X),10(f9.3,1X),4(f8.2,1X),7(f8.3,1X))

return
end subroutine out_accumulate