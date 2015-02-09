
!-------------------------------------------------------------------------
!Snowfraction added as a comment text. LJ 15 Jan 2013

SUBROUTINE OHMnew
  use SUES_data     
  USE GIS_data      ! LUMPS_gis_read.f95
  use data_in      
  use ohm_calc	
  use time
 use  allocateArray    ! sg feb 2012 - add allocated arrays
  use defaultNotUsed
 implicit none
 integer::i,ii,isel
 
 real(kind(1d0))::dqn,surfrac
  !Define what coefficients from the OHM_coef file are used	
    a1=0
    a2=0
    a3=0    
    
    !Find if the day belongs to summer (limit is 5 degrees)
    ! check could change by latitude i.e. 5 deg C not always appropriate or ref for basis
    
    if (HDD(id-1,4)>=5) then    !summer half of year
       ii=0
    else
       ii=2
    endif
    
    if(BareSoilSurfFraction>0)then    
       is=GrassISurf
       if (state(is)>0) then   !wet surface  ! should add soil moisture
            i=ii+1
       else                    !dry surface
            i=ii+2
       endif   
       !Snow cover added
       a1=a1+BareSoilSurfFraction*OHM_coef(is,i,1)*(1-snowFrac(is))  !The actual, areally weighted coefficients
       a2=a2+BareSoilSurfFraction*OHM_coef(is,i,2)*(1-snowFrac(is))
       a3=a3+BareSoilSurfFraction*OHM_coef(is,i,3)*(1-snowFrac(is))      
	endif

    do isel=1,nsurf  !+1  each surface type through,3d not possible at the moment!                        
        if(isel==GrassUSurf)then                  	
            is=GrassISurf
            surfrac=sfr(GrassUSurf)-BareSoilSurfFraction
        else
            is=isel
            surfrac=sfr(is)
        endif  
        
        if (state(is)>0) then   !wet surface  ! should add soil moisture
            i=ii+1
        else                   !dry surface
            i=ii+2
        	if(is>BldgSurf.and.is/=WaterSurf)then    ! wet soil  
       			 if(soilmoist(is)/soilstoreCap(is)>0.9)then   ! ?? check what the best value for this should be 
          			i=ii+1
                 endif
            endif             		 
        endif    

       !If snow on ground,this affects the surface fraction 
       !Buildings needs to be taken into account, what about paved ??
       !if (is>2.and.is/=WaterSurf.and.snowUse==1) then
       if (is/=BldgSurf.and.is/=WaterSurf.and.snowUse==1) then
           surfrac=surfrac*(1-snowFrac(is))
       endif
      
       a1=a1+surfrac*OHM_coef(is,i,1)  !The actual, areally weighted coefficients
       a2=a2+surfrac*OHM_coef(is,i,2)
       a3=a3+surfrac*OHM_coef(is,i,3)      
      
    enddo  

  if(qn1>-333)then
  	  qs=NAN  
  	  if(q1>-333.and.q3>-333)then
     	 dqn=(q3-q1)/2         !Difference between the future and past timesteps
     	 qs=QN1*a1+dqn*a2+a3
  	  endif
  	  q1=q2
  	  q2=q3
  	  q3=qn1
  endif    

  !Do snow separately is calculations are made! Added by LJ in August 2013
  if (qn1_S>-333.and.snowUse==1) then
 	  deltaQi=NAN  
  	  if(r1>-333.and.r3>-333)then
     	 dqn=(r3-r1)/2         !Difference between the future and past timesteps
     	 deltaQi=qn1_S*OHM_coef(nsurf+2,1,1)+dqn*OHM_coef(nsurf+2,1,2)+OHM_coef(nsurf+2,1,3)
  	  endif
  	  r1=r2
  	  r2=r3
  	  r3=qn1_S
  endif
  
      
  !Wind speed is not now embedded in the model
  !Seasons may need to be fixed
 return
ENDSUBROUTINE OHMnew


!! ---- No longer needed after v2014b (HCW) ----
! SUBROUTINE OHMinitialize
!!This subroutine deals with the new OHM calculation where
!!surface state and wind speed are taken into account
!!Made by LJ in 30 Aug 2011 according to the matlab codes by SG
!!---------------------------------------------------------------
!  USE GIS_data      ! LUMPS_gis_read.f95
!  use data_in      !
!  use ohm_calc
!  use allocateArray   ! sg feb 2012 - add allocated arrays
!  use time
!  use defaultNotUsed
!
!  IMPLICIT NONE
!
!  integer:: i, ii, iii, rr
!  ! check this is not being used -subvalues
!
!  !-------OPEN OHM FILES-------------------------
!  !The length of the user defined file is fixed
!  open(7,file=trim(fileOHM),status='old',err=200)
!  !  +1 is  -- soil
!  write(12,*)'-------',trim(fileOHM),'----------------------'
!  ! a1,a2,a3, code for the line of coefficients
!  do i=1,4
!	READ(7,*,iostat=iostat_var) (co2use(i,rr),rr=1,nsurf+2)
!    write(12,'(8g12.4)')(co2use(i,rr),rr=1,nsurf+2)
!  enddo
!  close(7)
!
! !Calculate the actual coefficients
! OHM_coef=0
!
! do ii=1,nsurf+2   !cols
!    do i=1,4 ! rows      (seasons/state)
!        if(co2use(i,ii)>0) then !If value is realistic start to compare with co
!			do iii=1,NrowOhm !Go each row in the co file through
!                if (co(iii,4)==co2use(i,ii)) then               
!					OHM_coef(ii,i,1) = co(iii,1) 
!                    OHM_coef(ii,i,2) = co(iii,2) 
!            		OHM_coef(ii,i,3) = co(iii,3)                    
!                endif
!            enddo
!        endif
!    enddo
! enddo
! ! this would have bare soil coefficient under UnIrrigated Grass -- but this should be fixed by water state
! ! canyon - in nsurf+1
! write(12,*)' OHM coefficients----------------------'
! do i=1,4
!    do ii=1, nsurf+1
!    	write(12,'(2i4,3g10.3)')ii,i, (OHM_coef(ii,i,iii),iii=1,3)
!    enddo
! enddo
!
! close(12)
! return
!
!200 call ErrorHint(47,trim(fileOHM),notUsed,notUsed,notUsedI)
!
! ENDSUBROUTINE OHMinitialize