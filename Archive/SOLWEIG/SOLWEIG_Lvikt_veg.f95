subroutine Lvikt_veg(isvf,isvfveg,isvfaveg,vikttot)
use matsize 
   
    implicit none
    real(kind(1D0)) :: vikttot
    real(kind(1d0)), dimension(sizey,sizex) :: isvf
    real(kind(1d0)), dimension(sizey,sizex) :: isvfveg
    real(kind(1d0)), dimension(sizey,sizex) :: isvfaveg
    real(kind(1d0)), dimension(sizey,sizex) :: viktonlywall
    real(kind(1d0)), dimension(sizey,sizex) :: viktaveg
    !real(kind(1d0)), dimension(sizey,sizex) :: viktwall
    real(kind(1d0)), dimension(sizey,sizex) :: svfvegbu    

    !allocate(svfalfaE(sizex,sizey))
    !allocate(svfalfaS(sizex,sizey))
    !allocate(svfalfaW(sizex,sizey))    
    !allocate(svfalfaN(sizex,sizey))    
    !allocate(alfaB(sizex,sizey))
    !allocate(betaB(sizex,sizey))    
    
    !! Least 
    viktonlywall=(vikttot-(63.227*isvf**6-161.51*isvf**5+156.91*isvf**4-70.424*isvf**3+16.773*isvf**2-0.4863*isvf))/vikttot
   
    viktaveg=(vikttot-(63.227*isvfaveg**6-161.51*isvfaveg**5+156.91*isvfaveg**4-70.424*isvfaveg**3+16.773*isvfaveg**2&
         &-0.4863*isvfaveg))/vikttot
   
    viktwall=viktonlywall-viktaveg
   
    svfvegbu=(isvfveg+isvf-1) ! Vegetation plus buildings
    viktsky=(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863*svfvegbu)/vikttot
    viktrefl=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863&
         &*svfvegbu))/vikttot
    viktveg=(vikttot-(63.227*svfvegbu**6-161.51*svfvegbu**5+156.91*svfvegbu**4-70.424*svfvegbu**3+16.773*svfvegbu**2-0.4863&
         &*svfvegbu))/vikttot
    viktveg=viktveg-viktwall
   
end subroutine Lvikt_veg 
