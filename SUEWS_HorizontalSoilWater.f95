
subroutine HorizontalSoilWater
!Program tranfer water in soil stores of land surfaces LJ (2010)
!Change the model to use varying hydraulic conductivity instead of constant value LJ (7/2011)
!If one of the surface's soildepth is zero, not water movement is considered 
!------------------------------------------------------
use SUES_data
use gis_data
use time
use allocateArray
 
implicit none

integer::jj

real(kind(1d0)) ::SoilMoist_vol1,SoilMoistCap_vol1,SoilMoist_vol2,SoilMoistCap_vol2,&
 		MatPot1,DimenWaterCon1,MatPot2,DimenWaterCon2,Distance,B_r1,B_r2,Km1,Km2,KmWeight 
        
 runoffSoil=0.0
 do is=1,nsurf-1  ! exclude water surface
      
    if (sfr(is)/=0.and.SoilStoreCap(is)>0) then !If particular surface area exists
      
    	do jj=is+1,nsurf-1

           if (sfr(jj)/=0.and.SoilStoreCap(jj)>0) then !If other surface area exists
             
				!Calculate Soil depth in order to calculate non-saturated volumetric water content
                SoilDepth=SoilStoreCap(is)/VolSoilMoistCap(is)
                
   				!SoilMoist_vol1=SoilMoist(is)/SoilDepth       !Volumetric soil moisture
    			!SoilMoistCap_Vol1=SoilStoreCap(is)/SoilDepth !Volumetric soil moisture capacity
                
                SoilMoistCap_Vol1=VolSoilMoistCap(is)  !Volumetric soil moisture capacity
                SoilMoist_vol1=SoilMoist(is)/SoilDepth !Volumetric soil moisture
                
                B_r1=SoilMoistCap_Vol1-SoilMoist_vol1
  
                if (B_r1<SoilMoist_vol1) then 
    				DimenWaterCon1=(SoilMoist_vol1-B_r1)/(SoilMoistCap_Vol1-B_r1)!Dimensionless water content
                    
                    if (DimenWaterCon1>0.99999) DimenWaterCon1=DimenWaterCon1-0.0001 !This cannot equal 1!!!
                    ! write(*,*)DimenWaterCon1,SoilMoist(is),SoilDepth ,SoilStoreCap(is), B_r1,is
                    ! pause
                    if(DimenWaterCon1<0.00000005) then
                      call ErrorHint(18,'HorizontalSoilWater - check soil moisture capacity & volumetric water content',&
                      				 SoilStoreCap(is), soilmoist(is),is)
                    endif
                    MatPot1=sqrt(1/DimenWaterCon1**2-1)/0.0005 !Water potential of the first store
                    
                    !Hydraulic conductivity of the first store (eq 8 in van Genuchten 1980)
                    !m = 1-1/n where n=2 meaning that m=0.5
                    Km1=SatHydraulicConduct(is)*sqrt(DimenWaterCon1)*(1-(1-DimenWaterCon1**2)**0.5)**2
                    
                    if (MatPot1>100000) MatPot1=100000! Max. potential is 100 000 (van Genuchten 1980)
           
                else
                    MatPot1 = 100000
                    Km1 = 0 !Added by LJ in Nov 2013
                endif
                
                !-------------------------------------------------------------------------
                SoilDepth=SoilStoreCap(jj)/VolSoilMoistCap(jj) ! Soil depth of the othe surface
                
                SoilMoistCap_Vol2=VolSoilMoistCap(jj)
        		SoilMoist_vol2=SoilMoist(jj)/SoilDepth !Volumetric soil moisture
    			
                 
                
                B_r2=SoilMoistCap_Vol2-SoilMoist_vol2
                
                if (B_r2<SoilMoist_vol2) then 
    				DimenWaterCon2=(SoilMoist_vol2-B_r2)/(SoilMoistCap_Vol2-B_r2) !Dimensionless water content
                
                	if (DimenWaterCon2>0.99999) DimenWaterCon2=DimenWaterCon2-0.0001 !This cannot equal 1!!!
                
    				MatPot2=sqrt(1/DimenWaterCon2**2-1)/0.0005 !Water potential of surface two
 
					!Hydraulic conductivity of the second store (eq 8 in van Genuchten 1980)
                    !m = 1-1/n where n=2 meaning that m=0.5
                    Km2=SatHydraulicConduct(jj)*sqrt(DimenWaterCon2)*(1-(1-DimenWaterCon2**2)**0.5)**2
                    
				    if ((MatPot2)>100000) MatPot2=100000 ! Max. potential is 100 000 (van Genuchten 1980)
					
                else
                    MatPot2=100000
                    Km2 = 0 !Added by LJ in Nov 2013
                endif
                
                Distance=(sqrt(sfr(is)*SurfaceArea)+sqrt(sfr(jj)*SurfaceArea))/2!Distance of the two stores
                
				!Areally weighted new Km
                KmWeight=(sfr(is)*Km1+sfr(jj)*Km2)/(sfr(is)+sfr(jj))
                
				dI_dt=-(KmWeight)*(-MatPot1+MatPot2)/Distance !Water flow between the stores

                !Water moves only if there is space for water
                if ((SoilMoist(jj)>=dI_dt*sfr(is)/sfr(jj)).and.((SoilMoist(is)+dI_dt)>=0)) then
                	SoilMoist(is)=SoilMoist(is)+dI_dt
        			SoilMoist(jj)=SoilMoist(jj)-dI_dt*sfr(is)/sfr(jj)
                    
                elseif ((SoilMoist(is)+dI_dt)<0) then
                    SoilMoist(is)=0
                    SoilMoist(jj)=SoilMoist(jj)+SoilMoist(is)*sfr(is)/sfr(jj)

                else
                    SoilMoist(is)=SoilMoist(is)+SoilMoist(jj)*sfr(jj)/sfr(is)
                    SoilMoist(jj)=0
                endif
                

        		if (SoilMoist(is)>SoilStoreCap(is)) then
            		runoffSoil(is)=runoffSoil(is)+(SoilMoist(is)-SoilStoreCap(is))
					SoilMoist(is)=SoilStoreCap(is)
        		elseif (SoilMoist(is)<0) then
        			SoilMoist(is)=0
        		endif

        		if (SoilMoist(jj)>SoilStoreCap(jj)) then
            		runoffSoil(jj)=runoffSoil(jj)+(SoilMoist(jj)-SoilStoreCap(jj))		
					SoilMoist(jj)=SoilStoreCap(jj)
        		elseif (SoilMoist(jj)<0) then
                 call errorHint(30,'Subroutine HorizontalWater: [soilmoisture(is)<0- set to 0],dectime,SoilMoist(is),is',  &
                 dectime,SoilMoist(is),is)                		
        			SoilMoist(jj)=0
        		endif
           endif
          
    enddo
    endif 
 runoffSoil_per_interval=runoffSoil_per_interval+(runoffSoil(is)*sfr(is))
     
 enddo

end subroutine