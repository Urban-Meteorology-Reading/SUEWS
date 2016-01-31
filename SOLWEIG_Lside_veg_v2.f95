subroutine Lside_veg_v2(altitude,Ta,Tw,SBC,ewall,esky,t,CI)!azimuth,
! This m-file is the current one that estimates L from the four cardinal points 20100414
use matsize
use data_in

    implicit none 
    real(kind(1D0))                             :: vikttot,aziE,aziN,aziS,aziW
    real(kind(1D0))                             :: altitude,Ta,Tw,SBC,ewall,esky,t,CI,c!azimuth,
    real(kind(1d0)),allocatable,dimension(:,:)  :: svfalfaE,svfalfaS,svfalfaW,svfalfaN
    real(kind(1d0)),allocatable,dimension(:,:)  :: alfaB,betaB,betasun
    real(kind(1d0)),allocatable,dimension(:,:)  :: Lground,Lrefl,Lsky,Lsky_allsky,Lveg,Lwallsh,Lwallsun
    real(kind(1d0)),allocatable,dimension(:,:)  :: viktonlywall,viktaveg,svfvegbu
    real(kind(1d0)),allocatable,dimension(:,:)  :: oneminussvfE,oneminussvfS,oneminussvfW,oneminussvfN
    real(kind(1d0)),parameter                   :: pi=3.141592653589793

    allocate(oneminussvfE(sizex,sizey))
    allocate(oneminussvfS(sizex,sizey))
    allocate(oneminussvfW(sizex,sizey))    
    allocate(oneminussvfN(sizex,sizey))    
    allocate(svfalfaE(sizex,sizey))
    allocate(svfalfaS(sizex,sizey))
    allocate(svfalfaW(sizex,sizey))    
    allocate(svfalfaN(sizex,sizey))    
    allocate(alfaB(sizex,sizey))
    allocate(betaB(sizex,sizey))
    allocate(betasun(sizex,sizey))    
    allocate(Lground(sizex,sizey))    
    allocate(Lrefl(sizex,sizey))
    allocate(Lsky(sizex,sizey))
    allocate(Lsky_allsky(sizex,sizey))
    allocate(Lveg(sizex,sizey))    
    allocate(Lwallsh(sizex,sizey))    
    allocate(Lwallsun(sizex,sizey))
    allocate(viktonlywall(sizex,sizey))
    allocate(viktaveg(sizex,sizey))    
    allocate(svfvegbu(sizex,sizey))

    if (allocated(viktwall)) deallocate(viktwall); allocate(viktwall(sizey,sizex))
    if (allocated(viktsky)) deallocate(viktsky); allocate(viktsky(sizey,sizex))
    if (allocated(viktveg)) deallocate(viktveg); allocate(viktveg(sizey,sizex))
    if (allocated(viktrefl)) deallocate(viktrefl); allocate(viktrefl(sizey,sizex))
    
    !Building height angle from svf
    oneminussvfE=1.-svfE; where (oneminussvfE<=0) oneminussvfE=0.000000001 ! avoiding log(0)
    oneminussvfS=1.-svfS; where (oneminussvfS<=0) oneminussvfS=0.000000001 ! avoiding log(0)
    oneminussvfW=1.-svfW; where (oneminussvfW<=0) oneminussvfW=0.000000001 ! avoiding log(0)
    oneminussvfN=1.-svfN; where (oneminussvfN<=0) oneminussvfN=0.000000001 ! avoiding log(0)
    svfalfaE=asin(exp((log(oneminussvfE))/2))
    svfalfaS=asin(exp((log(oneminussvfS))/2))
    svfalfaW=asin(exp((log(oneminussvfW))/2))
    svfalfaN=asin(exp((log(oneminussvfN))/2))
   
    vikttot=4.4897
    aziW=azimuth+t
    aziN=azimuth-90+t
    aziE=azimuth-180+t
    aziS=azimuth-270+t
   
    ! F_sh=cylindric_wedge(zen,svfalfa);!Fraction shadow on building walls based on sun altitude and svf 
    ! F_sh(isnan(F_sh))=0.5; 
    F_sh=2.*F_sh-1. !(cylindric_wedge scaled 0-1)
    
     if (SOLWEIG_ldown==1) then 
        c=1-CI
        Lsky_allsky=esky*SBC*((Ta+273.15)**4)*(1-c)+c*SBC*((Ta+273.15)**4) 
    else
        Lsky_allsky=ldown
    end if
    
    !! Least 
    call Lvikt_veg(svfE,svfEveg,svfEaveg,vikttot) 
   
    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaE)
        betaB=atan(tan((svfalfaE)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (180-t) .and. azimuth <= (360-t)) then 
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziE*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if 
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if 
    Lsky=((svfE+svfEveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Least=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
   
    !! Lsouth 
    call Lvikt_veg(svfS,svfSveg,svfSaveg,vikttot) 
   
    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaS)
        betaB=atan(tan((svfalfaS)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth <= (90-t) .or. azimuth > (270-t)) then 
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziS*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if 
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if 
    Lsky=((svfS+svfSveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lsouth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
   
    !! Lwest 
    call Lvikt_veg(svfW,svfWveg,svfWaveg,vikttot) 
   
    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaW)
        betaB=atan(tan((svfalfaW)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (360-t) .or. azimuth <= (180-t)) then 
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziW*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if 
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if 
    Lsky=((svfW+svfWveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lwest=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl
   
    !! Lnorth 
    call Lvikt_veg(svfN,svfNveg,svfNaveg,vikttot) 
   
    if (altitude>0) then ! daytime
        alfaB=atan(svfalfaN)
        betaB=atan(tan((svfalfaN)*F_sh))
        betasun=((alfaB-betaB)/2)+betaB
        if (azimuth > (90-t) .and. azimuth <= (270-t)) then 
            Lwallsun=SBC*ewall*((Ta+273.15+Tw*sin(aziN*(pi/180)))**4)* &
            viktwall*(1-F_sh)*cos(betasun)*0.5
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*F_sh*0.5
        else
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
        end if 
        else !nighttime
            Lwallsun=0
            Lwallsh=SBC*ewall*((Ta+273.15)**4)*viktwall*0.5
    end if 
    Lsky=((svfN+svfNveg-1)*Lsky_allsky)*viktsky*0.5
    Lveg=SBC*ewall*((Ta+273.15)**4)*viktveg*0.5
    Lground=Lup2d*0.5
    Lrefl=(Ldown2d+Lup2d)*(viktrefl)*(1-ewall)*0.5
    Lnorth=Lsky+Lwallsun+Lwallsh+Lveg+Lground+Lrefl

    deallocate(svfalfaE)
    deallocate(svfalfaS)
    deallocate(svfalfaW)    
    deallocate(svfalfaN)    
    deallocate(alfaB)
    deallocate(betaB)
    deallocate(betasun)    
    deallocate(Lground)    
    deallocate(Lrefl)
    deallocate(Lsky)
    deallocate(Lsky_allsky)
    deallocate(Lveg)    
    deallocate(Lwallsh)    
    deallocate(Lwallsun)
    deallocate(viktonlywall)
    deallocate(viktaveg)    
    deallocate(svfvegbu)
    
end subroutine Lside_veg_v2 
