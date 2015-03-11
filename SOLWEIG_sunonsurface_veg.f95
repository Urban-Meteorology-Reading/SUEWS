subroutine sunonsurface_veg(iazimuthA, scale, first, second, psi)
! This m-file creates a boolean image of sunlit walls. 
! Shadows from both buildings and vegetation is accounted for 
! moving building in the direction of the sun

use matsize
    implicit none
    real(kind(1d0))             :: iazimuthA,iazimuth,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))             :: scale 
    integer                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2,n,first,second
    real(kind(1d0))             :: dx,dy,ds,absdx,absdy,psi !,dz
    real(kind(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))             :: seventimespibyfour
    real(kind(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos
    real(kind(1d0)),allocatable,dimension(:,:)  ::weightsumwall,weightsumsh,gvf1,gvf2
    real(kind(1d0)),allocatable,dimension(:,:)  ::f,tempsh,tempbu,tempbub,tempwallsun,tempb,sh1

    real(kind(1d0)), parameter  :: pi=3.141592653589793
    real(kind(1d0)), parameter  :: maxpos=10000000000.0
    
    allocate(weightsumwall(sizex,sizey))    
    allocate(weightsumsh(sizex,sizey))
    allocate(gvf1(sizex,sizey))
    allocate(gvf2(sizex,sizey))
    allocate(f(sizex,sizey))
    allocate(tempsh(sizex,sizey))    
    allocate(tempbu(sizex,sizey))
    allocate(tempbub(sizex,sizey))
    allocate(tempwallsun(sizex,sizey))
    allocate(tempb(sizex,sizey))    
    allocate(sh1(sizex,sizey))
 
    iazimuth=iazimuthA*(pi/180)
    !special cases
    if (iazimuth==0) then
        iazimuth=iazimuth+0.000001
    end if
    ! loop parameters     
    index=0
    f=buildings
    sh1=sh-(1-vegsh)*(1-psi)
    dx=0
    dy=0
    ds=0 

    tempsh = 0.0D0 
    tempbu = 0.0D0
    tempbub = 0.0D0
    tempwallsun = 0.0D0
    !sh = 0.0D0
    weightsumsh = 0.0D0
    weightsumwall = 0.0D0

    first=first*scale
    second=second*scale

    ! other loop parameters
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=sin(iazimuth)
    cosazimuth=cos(iazimuth)
    tanazimuth=tan(iazimuth)
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=abs(1./sinazimuth)
    dscos=abs(1./cosazimuth)
    
    !! The Shadow casting algorithm 
    do n=1,second 
        IF ((pibyfour <= iazimuth .and. iazimuth < threetimespibyfour) .or. (fivetimespibyfour&
             & <= iazimuth .and. iazimuth < seventimespibyfour)) THEN
            dy=signsinazimuth*index
            dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
            ds=dssin
        ELSE
            dy=signsinazimuth*abs(nint(index*tanazimuth))
            dx=-1.*signcosazimuth*index
            ds=dscos
        END IF 

        absdx=abs(dx)
        absdy=abs(dy)
   
        xc1=((dx+absdx)/2)+1
        xc2=(sizex+(dx-absdx)/2)
        yc1=((dy+absdy)/2)+1
        yc2=(sizey+(dy-absdy)/2)
        xp1=-((dx-absdx)/2)+1
        xp2=(sizex-(dx+absdx)/2)
        yp1=-((dy-absdy)/2)+1
        yp2=(sizey-(dy+absdy)/2)

        tempbu(xp1:xp2,yp1:yp2)=buildings(xc1:xc2,yc1:yc2) !moving building

        tempsh(xp1:xp2,yp1:yp2)=sh1(xc1:xc2,yc1:yc2) !moving shadow image
        f=min(f,tempbu) !utsmetning of buildings

        weightsumsh=weightsumsh+tempsh*f

        tempwallsun(xp1:xp2,yp1:yp2)=sunwall(xc1:xc2,yc1:yc2) !moving building wall in sun image
        tempb=tempwallsun*f
        where ((tempb+tempbub)>0) !tempbub=(tempb+tempbub)>0==1
            tempbub=1.
        end where
        
        weightsumwall=weightsumwall+tempbub

        if (index*scale==first) then 
            gvf1=(weightsumwall+weightsumsh)/first
            where (gvf1>1)
                gvf1=1.
            end where
        end if 
        index=index+1

    end do 
    gvf2=(weightsumsh+weightsumwall)/second
    where (gvf2>1)
        gvf2=1.
    end where

    ! Weighting 
    sos=(gvf1*0.5+gvf2*0.4)/0.9
    
    deallocate(weightsumwall)    
    deallocate(weightsumsh)
    deallocate(gvf1)
    deallocate(gvf2)
    deallocate(f)
    deallocate(tempsh)    
    deallocate(tempbu)
    deallocate(tempbub)
    deallocate(tempwallsun)
    deallocate(tempb)
    deallocate(sh1)  

end subroutine  sunonsurface_veg 


