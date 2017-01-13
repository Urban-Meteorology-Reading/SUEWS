
 subroutine HorizontalSoilWater
 !Transfers water in soil stores of land surfaces LJ (2010)
 !Change the model to use varying hydraulic conductivity instead of constant value LJ (7/2011)
 !If one of the surface's soildepth is zero, no water movement is considered
 ! HCW 12/08/2014 Modifications:  - Distance changed from m to mm in dI_dt calculation
 !                                - dI_dt [mm s-1] multiplied by no. seconds in timestep -> dI [mm]
 !                                - if MatPot is set to max. value (100000 mm), Km set to 0 mm s-1
 !                                - Provide parameters for residual volumetric soil moisture [m3 m-3]
 !                                   (currently hard coded as 0.1 m3 m-3 for testing)
 ! HCW 13/08/2014 Modifications:  - Order of surfaces reversed (for both is and jj loops)
 !                                - Number of units (e.g. properties) added to distance calculation
 !
 !------------------------------------------------------
 use SUES_data
 use gis_data
 use time
 use allocateArray
 
 IMPLICIT NONE

 integer::jj

 real(kind(1d0))::SoilMoist_vol1,SoilMoistCap_vol1,SoilMoist_vol2,SoilMoistCap_vol2,&
                  MatPot1,DimenWaterCon1,MatPot2,DimenWaterCon2,Distance,B_r1,B_r2,Km1,Km2,KmWeight,&
                  alphavG,dI, NUnits
                         
 ! SoilMoist_vol1,2     = Volumetric soil moisture [m3 m-3]        
 ! SoilMoistCap_vol1,2  = Volumetric soil moisture capacity [m3 m-3] (from FunctionalTypes)
 ! MatPot1,2            = Water potential (i.e. pressure head) of store [mm]
 ! DimenWaterCon1,2     = Dimensionless water content, or relative saturation [-]
 ! Distance             = Distance between two stores [m]
 ! B_r1,2               = Residual volumetric soil moisture [m3 m-3]
 ! Km1,2                = Hydraulic conductivity of store [mm s-1]
 ! KmWeight             = Weighted hydraulic conductivity [mm s-1]
 ! alphavG              = Parameter (could depend on soil texture) [mm-1]
 ! dI                   = Water flow between stores [mm] dI = dI_dt * no. secs in each timestep
 !                         if dI > 0, first surface gains water, second surface loses water  
 ! NUnits               = Number of repeating units (e.g. properties, blocks) for distance calculation [-]
               
 alphavG=0.0005  !Set alphavG to match value in van Genuchten (1980) [mm-1]                
 runoffSoil=0.0  !Set soil runoff to zero [mm]
 
 NUnits = 1   !Can change to represent plot/base unit size
 
 do is=1,nsurf-1 !nsurf-1,1,-1  !Loop through each surface, excluding water surface (runs backwards as of 13/08/2014, HCW)
    
    if (sfr(is)/=0.and.SoilStoreCap(is)>0) then  !If particular surface area exists 
                                                 ! and is capable of storing water (SoilStoreCap [mm])
        do jj=is+1,nsurf-1 !is-1,1,-1  !Sub-loop through remaining surfaces (runs backwards as of 13/08/2014, HCW)
           
           if (sfr(jj)/=0.and.SoilStoreCap(jj)>0) then  !If other surface area exists 
                                                        ! and is capable of storing water
              
              ! ---- For surface 1 -----------------------------------------------------
              ! Calculate non-saturated VWC
              SoilMoistCap_Vol1=SoilStoreCap(is)/SoilDepth(is) !Volumetric soil moisture capacity [m3 m-3] (i.e. saturated VWC)
              SoilMoist_vol1=SoilMoist(is)/SoilDepth(is) !Volumetric soil moisture [m3 m-3]
                
              !B_r1=SoilMoistCap_Vol1-SoilMoist_vol1  !Residual soil moisture content [m3 m-3]   
              B_r1=0.1 !HCW 12/08/2014 Temporary fix 
              ! Need to add residual soil moisture values to FunctionalTypes
              !B_r1=VolSoilMoistRes(is) !Residual soil moisture content [m3 m-3]
              
              if(B_r1<SoilMoist_vol1) then 
                 DimenWaterCon1=(SoilMoist_vol1-B_r1)/(SoilMoistCap_Vol1-B_r1) !Dimensionless water content [-]
                    
                 if(DimenWaterCon1>0.99999) then
                    DimenWaterCon1=DimenWaterCon1-0.0001 !This cannot equal 1
                 endif
                 ! write(*,*)DimenWaterCon1,SoilMoist(is),SoilDepth(is),SoilStoreCap(is), B_r1,is
                 ! pause
                 
                 !Check this (HCW 12/08/2014) - is this needed? Why is it not done for other surface?
                 if(DimenWaterCon1<0.00000005) then
                    call ErrorHint(18,'In HorizontalSoilWater.f95, check soil moisture capacity & soil moisture for this surface.',&
                        SoilStoreCap(is), soilmoist(is),is)
                 endif
                 
                 !van Genuchten (1980), with n=2 and m = 1-1/n = 1/2
                 !Water potential of first store [mm] (van Genuchten 1980, Eq 3 rearranged)
                 MatPot1=sqrt(1/DimenWaterCon1**2-1)/alphavG 
                    
                 !Hydraulic conductivity of first store [mm s-1] (van Genuchten 1980, Eq 8)
                 Km1=SatHydraulicConduct(is)*sqrt(DimenWaterCon1)*(1-(1-DimenWaterCon1**2)**0.5)**2
                    
                 !Check this value (HCW 12/08/2014)
                 if(MatPot1>100000) then
                    MatPot1 = 100000  !Max. potential is 100000 mm (van Genuchten 1980)
                    Km1 = 0   !Added by HCW 12/08/2014
                 endif
                
              else  !If soil moisture below residual value, set MatPot to maximum
                 MatPot1 = 100000
                 Km1 = 0 !Added by LJ in Nov 2013
              endif
              
              ! ---- For surface 2 -----------------------------------------------------
              ! Calculate non-saturated VWC
              SoilMoistCap_Vol2=SoilStoreCap(jj)/SoilDepth(jj) !Volumetric soil moisture capacity [m3 m-3] (i.e. saturated VWC)             
              SoilMoist_vol2=SoilMoist(jj)/SoilDepth(jj) !Volumetric soil moisture [m3 m-3]
    			
              !B_r2=SoilMoistCap_Vol2-SoilMoist_vol2  !Residual soil moisture content [m3 m-3]   
              B_r2=0.1 !HCW 12/08/2014 Temporary fix 
              ! Need to add residual soil moisture values to FunctionalTypes
              !B_r2=VolSoilMoistRes(jj) !Residual soil moisture content [m3 m-3]                 
                
              if(B_r2<SoilMoist_vol2) then 
                 DimenWaterCon2=(SoilMoist_vol2-B_r2)/(SoilMoistCap_Vol2-B_r2) !Dimensionless water content [-]
                
                 if(DimenWaterCon2>0.99999) then
                    DimenWaterCon2=DimenWaterCon2-0.0001 !This cannot equal 1
                 endif 
      
                 !van Genuchten (1980), with n=2 and m = 1-1/n = 1/2
                 !Water potential of second store [mm] (van Genuchten 1980, Eq 3 rearranged)           
                 MatPot2=sqrt(1/DimenWaterCon2**2-1)/alphavG
 
                 !Hydraulic conductivity of second store [mm s-1] (van Genuchten 1980, Eq 8)					
                 Km2=SatHydraulicConduct(jj)*sqrt(DimenWaterCon2)*(1-(1-DimenWaterCon2**2)**0.5)**2
                    
                 if((MatPot2)>100000) then
                    MatPot2=100000 !Max. potential is 100000 mm (van Genuchten 1980)
                    Km2 = 0   !Added by HCW 12/08/2014
                 endif
					
              else  !If soil moisture below residual value, set MatPot to maximum
                 MatPot2=100000
                 Km2 = 0 !Added by LJ in Nov 2013
              endif
              
              ! ------------------------------------------------------------------------
                
              !Find distance between the two stores (see Jarvi et al. 2011)
              !SurfaceArea in m2 (changed from ha to m2 n SUEWS_Initial), so Distance in m
              Distance=(sqrt(sfr(is)*SurfaceArea/NUnits)+sqrt(sfr(jj)*SurfaceArea/NUnits))/2
                                
			  !Calculate areally-weighted hydraulic conductivity [mm s-1]
              KmWeight=(sfr(is)*Km1+sfr(jj)*Km2)/(sfr(is)+sfr(jj))
                
			  !Find water flow between the two stores [mm s-1] (Green-Ampt equation, Hillel 1971)
			  !Multiply Distance by 1000 to convert m to mm (HCW 12/08/2014)
              dI_dt=-(KmWeight)*(-MatPot1+MatPot2)/(Distance*1000)
                
              !Multiply dI_dt by number of seconds in timestep to convert mm s-1 to mm
              !Added by HCW 12/08/2014
              dI=dI_dt*tstep_real  !Use dI instead of dI_dt in the following calculations
                          
              !Move water (in mm) ------------------------------------------------------
              !Water moves only if (i) there is sufficient water to move and (ii) there is space to move it
              
              ! If there is sufficient water in both surfaces, allow movement of dI to occur
              if ((SoilMoist(jj)>=dI*sfr(is)/sfr(jj)).and.((SoilMoist(is)+dI)>=0)) then
                 SoilMoist(is)=SoilMoist(is)+dI                   
                 SoilMoist(jj)=SoilMoist(jj)-dI*sfr(is)/sfr(jj)  !Check (HCW 13/08/2014) - why adjust for jj and not is?
              
              ! If insufficient water in first surface to move dI, instead move as much as possible
              elseif ((SoilMoist(is)+dI)<0) then
                 SoilMoist(jj)=SoilMoist(jj)+SoilMoist(is)*sfr(is)/sfr(jj) !HCW 12/08/2014 switched order of these two lines
                 SoilMoist(is)=0    !Check (HCW 13/08/2014) - can SM actually go to zero, or is this inconsistent with SMres?
                   
              ! If insufficient water in second surface to move dI, instead move as much as possible
              else
                 SoilMoist(is)=SoilMoist(is)+SoilMoist(jj)*sfr(jj)/sfr(is)
                 SoilMoist(jj)=0
              endif
                
              !If soil moisture exceeds capacity, excess goes to soil runoff (first surface) 
              if (SoilMoist(is)>SoilStoreCap(is)) then
                 runoffSoil(is)=runoffSoil(is)+(SoilMoist(is)-SoilStoreCap(is))
                 SoilMoist(is)=SoilStoreCap(is)
        	  !elseif (SoilMoist(is)<0) then  !HCW 13/08/2014 commented out as should never be true here anyway... 
        	  !   SoilMoist(is)=0             ! ... and if so, need to do more here (i.e. account for other water too)
              endif

              !If soil moisture exceeds capacity, excess goes to soil runoff (second surface) 
              if (SoilMoist(jj)>SoilStoreCap(jj)) then
                 runoffSoil(jj)=runoffSoil(jj)+(SoilMoist(jj)-SoilStoreCap(jj))
                 SoilMoist(jj)=SoilStoreCap(jj)
        	  !elseif (SoilMoist(jj)<0) then  !HCW 13/08/2014 commented out (as above)
        	  !	 SoilMoist(jj)=0
              endif
           
           endif  !end if second surface exists and is capable of storing water
          
        enddo  !end jj loop over second surface
    
    endif  !end if first surface exists and is capable of storing water
 
 runoffSoil_per_tstep=runoffSoil_per_tstep+(runoffSoil(is)*sfr(is)/NonWaterFraction)  !Excludes water body
     
 enddo !is loop over first surface

end subroutine