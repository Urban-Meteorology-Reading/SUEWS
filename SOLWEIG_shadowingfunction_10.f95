	!------------------------------------------------!
	! Shadow casting algorithm, ground and buildings !
    !------------------------------------------------!
    subroutine shadowingfunction_10(azimuth,altitude,scale)
    use matsize
	! This m.file calculates shadows on a DSM
    ! Code originate from Carlo Rattis thesis
    ! This code is translated from Matlab by
    ! Fredrik Lindberg, Gothenburg University
    ! Last modified LJ 27 Jan 2016 - Removal of tabs and fixing real-int conversions
    
    implicit none
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0)), parameter          :: maxpos=10000000000.0
    !real, allocatable, dimension(:,:)   :: a,f,temp !! Defined in matsize
    real(kind(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    real(kind(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    integer                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,row,col !,sizex,sizey
	! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: temp,f
    
    !special cases
    if (altitude==90) then
        altitude=altitude-0.0001
    end if
    if (azimuth==0) then
        azimuth=azimuth-0.0001
    end if
        
    ! conversion
    degrees=pi/180
    azi=azimuth*degrees
    alt=altitude*degrees

    allocate(f(sizex,sizey))    
    allocate(temp(sizex,sizey))
    
    if (allocated(sh)) deallocate(sh)
    allocate(sh(sizex,sizey))
       
    ! initialise parameters
    f=a
    dx=0
    dy=0
    dz=0
    temp=a*0.0
    index=1

    ! other loop parameters
    amaxvalue=maxval(a)
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
    tanaltitudebyscale=tan(alt)/scale


    ! main loop
    DO WHILE (amaxvalue>=dz .and. abs(dx)<=sizex .and. abs(dy)<=sizey)

        IF ((pibyfour <= azi .and. azi < threetimespibyfour) .or. (fivetimespibyfour <= azi .and. azi < seventimespibyfour)) THEN
            dy=signsinazimuth*index
            dx=-1.*signcosazimuth*abs(nint(index/tanazimuth))
            ds=dssin
        ELSE
            dy=signsinazimuth*abs(nint(index*tanazimuth))
            dx=-1.*signcosazimuth*index
            ds=dscos
        END IF

        dz=ds*index*tanaltitudebyscale
        temp=temp*0
    
        absdx=abs(dx)
        absdy=abs(dy)
   
        xc1=int((dx+absdx)/2)+1
        xc2=(sizex+int((dx-absdx)/2))
        yc1=int((dy+absdy)/2)+1
        yc2=(sizey+int((dy-absdy)/2))
        xp1=-int((dx-absdx)/2)+1
        xp2=(sizex-int((dx+absdx)/2))
        yp1=-int((dy-absdy)/2)+1
        yp2=(sizey-int((dy+absdy)/2))
    
        temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz
    
        f=max(f,temp)
        index=index+1
        
    END DO
 
    f=f-a
    where (f>0)
        f=-1
    end where
    sh=f+1
    !sh=f ! invert as in shadowingfunctionglobalradiation
    deallocate(f)    
    deallocate(temp)

end subroutine shadowingfunction_10