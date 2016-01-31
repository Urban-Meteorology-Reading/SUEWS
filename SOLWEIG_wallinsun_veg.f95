 !>
 subroutine wallinsun_veg(azimuth)
 ! This m-file creates a boolean image of sunlit walls.
 ! Shadows from both buildings and vegetation is accounted for
 ! moving building in the direction of the sun
 ! Last modified: 
 !  LJ 27 Jan 2017 - Change of equations xc1...yp2 to account for the change from real to integer
 !---------------------------------------------------------------------------------

 use matsize
    implicit none
    real(kind(1d0))             :: azimuth,iazimuth 
    integer                     :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2
    real(kind(1d0))             :: dx,dy,dz,ds,absdx,absdy
    real(kind(1d0))             :: pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))             :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))             :: signsinazimuth,signcosazimuth,dssin,dscos,azi
    real(kind(1d0)), parameter  :: pi=3.141592653589793
    real(kind(1d0)), parameter  :: maxpos=10000000000.0
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: temp,sh1
    
    allocate(temp(sizex,sizey))
    allocate(sh1(sizex,sizey))    
    !allocate(vegsh(sizex,sizey))
    !allocate(a(sizex,sizey))
    !allocate(buildings(sizex,sizey))
    
    if (allocated(sunwall)) deallocate(sunwall); allocate(sunwall(sizey,sizex))
    
    iazimuth=azimuth+180
    if (iazimuth>=360) then 
        iazimuth=iazimuth-360
    end if 
    !special cases
    if (iazimuth==0) then
        iazimuth=iazimuth+0.00001
    end if
    ! conversion into radians 
    azi=iazimuth*(pi/180)

    index=1
    dx=0.
    dy=0.
    dz=0.
    ds=0.

    temp = 0.0D0 

    ! other loop parameters
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour
    fivetimespibyfour=5.*pibyfour
    seventimespibyfour=7.*pibyfour
    sinazimuth=sin(azi)
    cosazimuth=cos(azi)
    tanazimuth=tan(azi)
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth)
    !signcosazimuth=cosazimuth/abs(cosazimuth)
    dssin=abs(1./sinazimuth)
    dscos=abs(1./cosazimuth)
    
    sh1=vegsh+sh-1.
    !! The Shadow casting algoritm 
    IF ((pibyfour <= azi .and. azi < threetimespibyfour) .or. (fivetimespibyfour <= azi .and. azi < seventimespibyfour)) THEN
        dy=signsinazimuth*index
        dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
        ds=dssin
    ELSE
        dy=signsinazimuth*abs(nint(index*tanazimuth))
        dx=-1.*signcosazimuth*index
        ds=dscos
    END IF 

    ! note: dx and dy represent absolute values while ds is an incremental value 

    absdx=abs(dx)
    absdy=abs(dy)
   
    xc1=int((dx+absdx)/2)+1  !LJ added int to the equation to account for the conversion from real to int
    xc2=(sizex+int((dx-absdx)/2))
    yc1=int((dy+absdy)/2)+1
    yc2=(sizey+int((dy-absdy)/2))
    xp1=-int((dx-absdx)/2)+1
    xp2=(sizex-int((dx+absdx)/2))
    yp1=-int((dy-absdy)/2)+1
    yp2=(sizey-int((dy+absdy)/2))

    temp(xp1:xp2,yp1:yp2)= buildings(xc1:xc2,yc1:yc2)

    sunwall=temp-buildings
    where (sunwall==1) !f1(f1==1)=0
          sunwall=0
    end where
    where (sunwall==-1) !f1(f1==-1)=1
          sunwall=1
    end where    
    sunwall=sh1*sunwall
    
    deallocate(temp)
    deallocate(sh1)
    
end subroutine  wallinsun_veg 
