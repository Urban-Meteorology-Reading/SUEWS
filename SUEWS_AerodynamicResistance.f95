 subroutine AerodynamicResistance(RA,AerodynamicResistanceMethod,StabilityMethod,RoughLen_heat,&
            ZZD,z0m,k2,AVU1,L_mod,Ustar,VegFraction,psyh) ! psyh is added. shiho
        
 ! Returns Aerodynamic resistance (RA) to the main program SUEWS_Calculations
 ! All ra equations reported in Thom & Oliver (1977)
 ! Modified by LJ
 !   -Removal of tabs and cleaning the code
 ! Modified by HCW 03 Dec 2015 - changed lower limit on ra from 2 s m-1 to 10 s m-1 (to avoid unrealistically high evaporation rates)
 ! Modified by LJ in 12 April to not to be used with global variables
 ! To Do:
 !       - Check whether the thresholds 2-200 s m-1 are suitable over a range of z0!! HCW 04 Mar 2015
 ! OUTPUT: RA - Aerodynamic resistance [s m^-1]
 ! INPUT:  AerodynamicResistanceMethod = Method to calculate RA
 !         StabilityMethod = defines the method to calculate atmospheric stability
 !         RoughLen_heat = Method to calculate heat roughness length
 !         ZZD = Displacement height (m)
 !         z0m = Aerodynamic roughness length (m)
 !         k2 = Power of Van Karman's constant (= 0.16 = 0.4^2)
 !         AVU1 = Mean wind speed
 !         L_mod = Obukhov length (m)
 !         Ustar = Friction velocity (m s-1)
 !         VegFraction = Fraction of vegetation
 !               (changed from veg_fr which also includes water surface by HCW 05 Nov 2015)
         

 use DefaultNotUsed
            
 IMPLICIT NONE

 real (kind(1d0))::psym,psyh,stab_fn_heat,stab_fn_mom,ZZD,z0m,k2,AVU1,L_mod,Ustar,RA,z0V,VegFraction, &
 muu=1.46e-5 !molecular viscosity
 integer::AerodynamicResistanceMethod,StabilityMethod,RoughLen_heat

            
 !1)Monteith (1965)-neutral stability
 if(AerodynamicResistanceMethod==1) then
     RA=(log(ZZD/z0m)**2)/(k2*AVU1)
     
  !2) Non-neutral stability
  !    PSIM - stability function for momentum
  !     PSYH - stability function for heat
  !    assuming stability functions the same for heat and water
 elseif(AerodynamicResistanceMethod==2) then  !Dyer (1974)
       
     psym=stab_fn_mom(StabilityMethod,ZZD/L_mod,zzd/L_mod)
     psyh=stab_fn_heat(StabilityMethod,ZZD/L_mod,zzd/L_mod)

     !Z0V roughness length for vapour
     if (RoughLen_heat==1) then !Brutasert (1982) Z0v=z0/10(see Grimmond & Oke, 1986)
        z0V=Z0m/10
         elseif (RoughLen_heat==2) then ! Kawai et al. (2007)
       	!z0V=Z0m*exp(2-(1.2-0.9*veg_fr**0.29)*(Ustar*Z0m/muu)**0.25) 
        ! Changed by HCW 05 Nov 2015 (veg_fr includes water; VegFraction = veg + bare soil)     
        z0V=Z0m*exp(2-(1.2-0.9*VegFraction**0.29)*(Ustar*Z0m/muu)**0.25)
         elseif (RoughLen_heat==3) then
        z0V=Z0m*exp(-20.) ! Voogt and Grimmond, JAM, 2000
     elseif (RoughLen_heat==4) then     
        z0V=Z0m*exp(2-1.29*(Ustar*Z0m/muu)**0.25) !See !Kanda and Moriwaki (2007),Loridan et al. (2010)
     endif
       
     if(Zzd/L_mod==0.or.Ustar==0) then
        RA=(log(ZZD/z0m)*log(ZZD/z0V))/(k2*AVU1) !Use neutral equation
     else
        RA=((log(ZZD/z0m)-PSYM)*(log(ZZD/z0V)-PSYH))/(K2*AVU1)
     endif
     
     !3) Thom and Oliver (1977)
  elseif(AerodynamicResistanceMethod==3) then
     RA=(4.72*log(ZZD/z0m)**2)/(1 + 0.54*AVU1)
  endif

  !If Ra too large  ! this was 175 (??check)  !!Check whether these thresholds are suitable over a range of z0
  
  if(RA>200) then           
     call errorHint(7,'RA > 200 s m-1; value set to 200 s m-1',RA,notUsed,notUsedI)
     RA=200
  elseif(RA<10)then   ! found  By Shiho - fix Dec 2012  !Threshold changed from 2 to 10 s m-1 (HCW 03 Dec 2015) 
     call errorHint(7,'RA < 10 s m-1; value set to 10 s m-1',RA,notUsed,notUsedI)
     RA=10
     ! RA=(log(ZZD/z0m))**2/(k2*AVU1)
     if(avu1<0) write(*,*) avu1,ra
  endif
  
  return
end subroutine AerodynamicResistance
