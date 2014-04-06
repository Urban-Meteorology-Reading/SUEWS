subroutine drainage(StorCap,DrainEq,DrainCoef1,DrainCoef2)
!Calculation of drainage for each land surface .Last modified LJ in 10/2010
use SUES_data
use gis_data
use time
use allocateArray
implicit none
real (Kind(1d0))::StorCap,DrainCoef1,DrainCoef2,DrainEq
!Modified by LJ in Aug 2011. Drainage cannot exceed the surface storage 

if(state(is)<0.000000001) then 
   drain(is)=0.0
else
  
   if(int(DrainEq)==1) then   !Falk and niemczynowicz (1978): Drainage equation for paved, buildings and irrigated grass
   	  if(state(is)<StorCap) then
         drain(is)=0
      else
         drain(is)=(DrainCoef1*(state(is)-StorCap)**DrainCoef2)/nsh
      endif
      
   elseif(int(DrainEq)==2) then   !     rutter eqn corrected for c=0
      drain(is)=(DrainCoef1*((exp(DrainCoef2*state(is)))-1))/nsh
      
   elseif(int(DrainEq)==3) then   !    falk and niemczynowicz (1978) equation
      drain(is)=(DrainCoef1*(state(is)**DrainCoef2))/nsh
      
   elseif(int(DrainEq)==4) then    !    rutter eqn not corrected for c=0
      drain(is)=DrainCoef1*exp(DrainCoef2*(state(is)-StorCap))
      
      drain(is)=drain(is)*nmin
   endif
   !drain(is)=drain(is)*sfr(is)
    !if(drain(is)>25)then
    !  write(*,*)'drainage very large'
    !  write(*,*)dectime,drain(is),state(is),int(draineq)
    !  pause
   !endif
   
   !Drainage cannot be larger than the surface state
   if (drain(is)>state(is)) then
       drain(is)=state(is)
   elseif(drain(is)<0.0001) then
       drain(is)=0
   endif
   
   dr_per_interval=dr_per_interval+drain(is)*sfr(is)
   

   !If surface area of type is 
   
!$$$$$$     if(drain(is)>0)then
!$$$$$$         write(*,*)'dr',drain(is),dr_per_interval,state(is)
!$$$$$$       pause
!$$$$$$     endif
endif

return
end subroutine drainage
