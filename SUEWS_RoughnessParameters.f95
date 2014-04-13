Subroutine RoughnessParameters(id)
use gis_data
use data_in
use sues_data
use allocateArray
use mod_z
use defaultNotUsed
implicit none
integer::id


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Get surface covers and frontal area fractions (LJ 11/2010)
!  areaZh=(sfr(BldgSurf)+sfr(ConifSurf)+sfr(DecidSurf)) !Total area of buildings and trees    
! sg feb 2012 -- made separate subroutine
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::       
     if (areaZh/=0)then
     Zh=bldgH*sfr(BldgSurf)/areaZh+treeH*sfr(ConifSurf)/areaZh+treeH*(1-porosity(id))*sfr(DecidSurf)/areaZh
     endif
     if(z0_Method==2) then  !Rule of thumb (G&O 1999)
        Z0m=0.1*Zh 
        Zdm=0.7*Zh
     elseif(z0_Method==3)then !MacDonald 1998
         if (areaZh/=0)then  !Plan area fraction
             planF=FAIBldg*sfr(BldgSurf)/areaZh+FAItree*sfr(ConifSurf)/areaZh+FAItree*(1-porosity(id))*sfr(DecidSurf)/areaZh           
         else
             planF=0.00001
             Zh=1
         endif
           
     	Zdm=(1+4.43**(-sfr(BldgSurf))*(sfr(BldgSurf)-1))*Zh
        Z0m=((1-Zdm/Zh)*exp(-(0.5*1.0*1.2/0.4**2*(1-Zdm/Zh)*planF)**(-0.5)))*Zh
     endif
     
     ZZD=Z-zdm

     ! error messages if aerodynamic parameters negative
     if(z0m<0) call ErrorHint(14,'SUEWS_RoughnessParameters',z0m,notUsed,notUsedI)
     if(zzd<0)call ErrorHint(15,'SUEWS_RoughnessParameters',zzd,notUsed,notUsedI)
end subroutine RoughnessParameters
