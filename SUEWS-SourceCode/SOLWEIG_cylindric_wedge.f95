subroutine cylindric_wedge(zen)

use matsize

! Fraction of sunlit walls based on sun altitude and svf wieghted building angles 

    implicit none 
    
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0))         :: zen      !>  
    real(kind(1d0))         :: beta      !> 
    real(kind(1d0)),allocatable,dimension(:,:) :: alfa,xa,ha,hkil,ba
    real(kind(1d0)),allocatable,dimension(:,:) :: Ai,phi,qa,Za
    real(kind(1d0)),allocatable,dimension(:,:) :: ukil,Ssurf 
    !real(kind(1d0)), dimension(sizey,sizex) :: 
    !real(kind(1d0)), dimension(sizey,sizex) :: 

    allocate(alfa(sizey,sizex))
    allocate(ba(sizey,sizex))
    allocate(ha(sizey,sizex))
    allocate(xa(sizey,sizex))
    allocate(qa(sizey,sizex))
    allocate(Za(sizey,sizex))
    allocate(phi(sizey,sizex))
    allocate(ukil(sizey,sizex))
    allocate(Ai(sizey,sizex))
    allocate(Ssurf(sizey,sizex))
    allocate(hkil(sizey,sizex))
    
    beta=zen
    alfa=svfalfa

    xa=1.-2./(tan(alfa)*tan(beta))
    ha=2./(tan(alfa)*tan(beta))
    ba=(1./tan(alfa))
    hkil=2.*ba*ha


    qa = 0.0D0 

    where (xa<0) !qa(xa<0)=tan(beta)/2
        qa=tan(beta)/2
    end where
    
   
    Za = 0.0D0
    
    phi = 0.0D0 
    
    Ai = 0.0D0
    
    ukil = 0.0D0 
    where (xa<0) 
        !Za(xa<0)=((ba(xa<0).**2)-((qa(xa<0).**2)./4)).**0.5
        Za=(ba**2-qa**2/4.)**0.5
        !phi(xa<0)=atan(Za(xa<0)./qa(xa<0))
        phi=atan(Za/qa)
        !A1(xa<0)=(sin(phi(xa<0))-phi(xa<0).*cos(phi(xa<0)))./(1-cos(phi(xa<0)))
        Ai=(sin(phi)-phi*cos(phi))/(1-cos(phi))
        !ukil(xa<0)=2*ba(xa<0).*xa(xa<0).*A1(xa<0)
        ukil=2*ba*xa*Ai
    end where

    Ssurf=hkil+ukil

    F_sh=(2*pi*ba-Ssurf)/(2*pi*ba) !Xa
    
    deallocate(alfa) 
    deallocate(ba) 
    deallocate(ha) 
    deallocate(xa) 
    deallocate(qa) 
    deallocate(Za) 
    deallocate(phi)
    deallocate(ukil) 
    deallocate(Ai) 
    deallocate(Ssurf)
    deallocate(hkil)

end subroutine cylindric_wedge 
