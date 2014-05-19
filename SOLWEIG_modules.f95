! 
module matsize
    
    IMPLICIT NONE
    
    integer								        	:: sizex,sizey
    real(kind(1d0)), allocatable, dimension(:,:) 	:: a,sh,vbshvegsh,vegsh
    real(kind(1d0)), allocatable, dimension(:,:) 	:: bush,vegdem,vegdem2,tempgrid
    real(kind(1d0)), allocatable, dimension(:,:) 	:: buildings,svf,svfE,svfS,svfW,svfN 
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
    real(kind(1d0)), allocatable, dimension(:,:) 	:: Kdown,Keast,Knorth,Ksouth,Kup,Kwest
    real(kind(1d0)), allocatable, dimension(:,:) 	:: Ldown,Least,Lnorth,Lsouth,Lup,Lwest
    real(kind(1d0)), allocatable, dimension(:,:) 	:: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfalfa,sos,Tgmap1!Tgmapgvf
    real(kind(1d0)), allocatable, dimension(:,:) 	:: viktveg,viktsky,viktrefl,viktwall
    
    real(kind(1d0)) :: CIlatenight,timeadd,firstdaytime

!svfvegbu,,viktaveg,viktonlywall,svfalfaE,svfalfaN,svfalfaS,svfalfaW,alfaB,betaB,betasun,Lground,Lrefl,Lsky,Lveg,Lwallsh,Lwallsun
! alfa,xa,ha,hkil,Ai,phi,qa,Za,tempsh,tempbu,tempbub,tempwallsun,weightsumwall,weightsumsh,gvf1,gvf2
!!,svfbuveg
!Lnight,,tempvegdem,tempvegdem2,fabovea,gabovea,tempbush,firstvegdem,vegsh2,tempgrid
!svfviktbuveg
!,stopbuild,stopveg,g,bushplant,Tgmap
!,temp,f,tmp,Knight,,tempb
    
end module matsize
