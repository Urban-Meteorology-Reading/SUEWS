subroutine Kside_veg_v24(radI,radG,radD,azimuth,altitude,psi,t,albedo)
use matsize 
    
    implicit none
    
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1D0))	:: vikttot,aziE,aziN,aziS,aziW 
    real(kind(1D0))	:: radI,radG,radD
    real(kind(1D0))	:: azimuth,altitude,psi,t,albedo
    
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: svfviktbuveg
  
    allocate(svfviktbuveg(sizex,sizey))
    ! New reflection equation 2012-05-25 
    vikttot=4.4897
    aziE=azimuth+t
    aziS=azimuth-90+t
    aziW=azimuth-180+t
    aziN=azimuth-270+t
    ! sunw=cos(altitude*(pi/180)); ! anngle to walls 
    !! Kside with weights 
    call Kvikt_veg( svfE,svfEveg,vikttot) 
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (360-t) .or. azimuth <= (180-t)) then 
        Keast=radI*shadow*cos(altitude*(pi/180))*sin(aziE*(pi/180))+ &
        radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Keast=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if 

    call Kvikt_veg( svfS,svfSveg,vikttot) 
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (90-t) .and. azimuth <= (270-t)) then 
    Ksouth=radI*shadow*cos(altitude*(pi/180))*sin(aziS*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Ksouth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if 

    call Kvikt_veg( svfW,svfWveg,vikttot) 
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth > (180-t) .and. azimuth <= (360-t)) then 
    Kwest=radI*shadow*cos(altitude*(pi/180))*sin(aziW*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Kwest=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if 

    call Kvikt_veg( svfN,svfNveg,vikttot) 
    svfviktbuveg=(viktwall+(viktveg)*(1-psi))
    if (azimuth <= (90-t) .or. azimuth > (270-t)) then 
    Knorth=radI*shadow*cos(altitude*(pi/180))*sin(aziN*(pi/180))+ &
    radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    else
    Knorth=radD*(1-svfviktbuveg)+radG*albedo*svfviktbuveg*(1-F_sh) !*sin(altitude*(pi/180));
    end if 
    
    deallocate(svfviktbuveg)
    end subroutine Kside_veg_v24 

    
subroutine Kvikt_veg(isvf, isvfveg, vikttot)
use matsize

    implicit none
    real(kind(1D0))	:: vikttot     
    real(kind(1d0)), dimension(sizey,sizex) :: isvf
    real(kind(1d0)), dimension(sizey,sizex) :: isvfveg
    real(kind(1d0)), dimension(sizey,sizex) :: svfvegbu 
    
    !! Least 
    viktwall=(vikttot-(63.227*isvf**6-161.51*isvf**5+156.91*isvf**4-70.424*isvf**3+16.773*isvf**2-0.4863*isvf))/vikttot

    svfvegbu=(isvfveg+isvf-1) ! Vegetation plus buildings
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2&
         &-0.4863*svfvegbu))/vikttot
    viktveg=viktveg-viktwall
end subroutine Kvikt_veg