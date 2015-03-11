! 
module matsize
    
    IMPLICIT NONE
    
    integer				        	:: sizex,sizey ! number of rows and cols of grid
    real(kind(1d0)), allocatable, dimension(:,:) 	:: a,sh,vbshvegsh,vegsh
    real(kind(1d0)), allocatable, dimension(:,:) 	:: bush,vegdem,vegdem2,tempgrid
    real(kind(1d0)), allocatable, dimension(:,:) 	:: buildings,svf,svfE,svfS,svfW,svfN 
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfveg,svfEveg,svfSveg,svfWveg,svfNveg
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg,last
    real(kind(1d0)), allocatable, dimension(:,:) 	:: Kdown2d,Keast,Knorth,Ksouth,Kup2d,Kwest
    real(kind(1d0)), allocatable, dimension(:,:) 	:: Ldown2d,Least,Lnorth,Lsouth,Lup2d,Lwest
    real(kind(1d0)), allocatable, dimension(:,:) 	:: gvf,Tmrt,shadow,Sstr,F_sh,sunwall
    real(kind(1d0)), allocatable, dimension(:,:) 	:: svfalfa,sos,Tgmap1
    real(kind(1d0)), allocatable, dimension(:,:) 	:: viktveg,viktsky,viktrefl,viktwall,savegrid
    
end module matsize

!svfvegbu,,viktaveg,viktonlywall,svfalfaE,svfalfaN,svfalfaS,svfalfaW,alfaB,betaB,betasun,Lground,Lrefl,Lsky,Lveg,Lwallsh,Lwallsun
! alfa,xa,ha,hkil,Ai,phi,qa,Za,tempsh,tempbu,tempbub,tempwallsun,weightsumwall,weightsumsh,gvf1,gvf2
!!,svfbuveg
!Lnight,,tempvegdem,tempvegdem2,fabovea,gabovea,tempbush,firstvegdem,vegsh2,tempgrid
!svfviktbuveg
!,stopbuild,stopveg,g,bushplant,Tgmap
!,temp,f,tmp,Knight,,tempb
!Tgmapgvf
    
module solweig_module
    IMPLICIT NONE
    
    real(kind(1d0)) :: timestepdec,& !time step in decimal time
                       CIlatenight,&
                       timeadd,&
                       firstdaytime,& ! if new day starts, =1 else =0
                       Fside,& ! fraction of a person seen from each cardinal point
                       Fup,& ! fraction of a person seen from down and up
                       scale,&
                       amaxvalue,&
                       trans,&
                       transperlai,&
                       xllcorner,&
                       yllcorner,&
                       NoData,&
                       cellsize
    integer         :: SolweigCount
    
end module solweig_module
