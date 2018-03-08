    !------------------------------------------------!
    ! Shadow casting algorithm, vegetation           !
    !------------------------------------------------!

subroutine shadowingfunction_20(azimuth,altitude,scale,amaxvalue) 
use matsize
    ! This m.file calculates shadows on a DSM and for vegetation units
    ! This code is translated from Matlab by Fredrik Lindberg, Gothenburg University
    
    implicit none
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    real(kind(1d0)), parameter          :: maxpos=10000000000.0
    real(kind(1d0))                     :: degrees,azi,alt,dx,dy,dz,ds,absdx,absdy,azimuth,altitude
    real(kind(1d0))                     :: amaxvalue,pibyfour,threetimespibyfour,fivetimespibyfour
    real(kind(1d0))                     :: seventimespibyfour,sinazimuth,cosazimuth,tanazimuth
    real(kind(1d0))                     :: signsinazimuth,signcosazimuth,dssin,dscos,tanaltitudebyscale,scale
    integer                             :: index,xc1,xc2,yc1,yc2,xp1,xp2,yp1,yp2!,test!,row,col !,sizex,sizey
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: f,temp,tmp,stopbuild,stopveg,g,bushplant,tempvegdem,tempvegdem2
    real(kind(1d0)),allocatable,dimension(:,:) :: fabovea,gabovea,tempbush,firstvegdem,vegsh2   
    
    !real 		  				  	    :: start_time,end_time
    
    !special case
    if (altitude==90) then
        altitude=altitude-0.0001
    end if
    if (azimuth==0) then
        azimuth=azimuth-0.0001
    end if
        
    ! conversion
    degrees=pi/180;
    azi=azimuth*degrees;
    alt=altitude*degrees;

    if (allocated(sh)) deallocate(sh)
    allocate(sh(sizex,sizey))
    if (allocated(vegsh)) deallocate(vegsh)
    allocate(vegsh(sizex,sizey))
    if (allocated(vbshvegsh)) deallocate(vbshvegsh)
    allocate(vbshvegsh(sizex,sizey))
    
	! allocation of grids
    allocate(f(sizex,sizey)) 
    allocate(temp(sizex,sizey))
    allocate(tmp(sizex,sizey))
    allocate(stopbuild(sizex,sizey))
    allocate(stopveg(sizex,sizey))
    allocate(g(sizex,sizey))
    allocate(bushplant(sizex,sizey))
    allocate(tempvegdem(sizex,sizey))
    allocate(tempvegdem2(sizex,sizey))
    allocate(fabovea(sizex,sizey))
    allocate(gabovea(sizex,sizey))
    allocate(firstvegdem(sizex,sizey))
    allocate(tempbush(sizex,sizey))
    allocate(vegsh2(sizex,sizey))
    
    ! initialise parameters
    f=a
    dx=0
    dy=0
    dz=0
    temp=a*0.0
    sh=temp
    vegsh=sh
    stopbuild=sh
    stopveg=sh
    vbshvegsh=sh
    g=sh
    bushplant=temp
    where (bush>1)
        bushplant=1
    end where
    
    index=1
	!test=0
    ! other loop parameters
    !amaxvalue=maxval(a)
    pibyfour=pi/4.
    threetimespibyfour=3.*pibyfour;
    fivetimespibyfour=5.*pibyfour;
    seventimespibyfour=7.*pibyfour;
    sinazimuth=sin(azi);
    cosazimuth=cos(azi);
    tanazimuth=tan(azi);
    call issign(sinazimuth,maxpos,signsinazimuth)
    call issign(cosazimuth,maxpos,signcosazimuth)
    !signsinazimuth=sinazimuth/abs(sinazimuth);
    !signcosazimuth=cosazimuth/abs(cosazimuth);
    dssin=abs(1./sinazimuth);
    dscos=abs(1./cosazimuth);
    tanaltitudebyscale=tan(alt)/scale;

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
        tempvegdem=temp
        tempvegdem2=temp
        
        absdx=abs(dx)
        absdy=abs(dy)
    
        xc1=int(((dx+absdx)/2))+1
        xc2=(sizex+int((dx-absdx)/2))
        yc1=int((dy+absdy)/2)+1
        yc2=(sizey+int((dy-absdy)/2))
        xp1=-int((dx-absdx)/2)+1
        xp2=(sizex-int((dx+absdx)/2))
        yp1=-int((dy-absdy)/2)+1
        yp2=(sizey-int((dy+absdy)/2))

        temp(xp1:xp2,yp1:yp2)= a(xc1:xc2,yc1:yc2)-dz
        tempvegdem(xp1:xp2,yp1:yp2)=vegdem(xc1:xc2,yc1:yc2)-dz
        tempvegdem2(xp1:xp2,yp1:yp2)=vegdem2(xc1:xc2,yc1:yc2)-dz

        f=max(f,temp)
        where (f>a) !sh(f>a)=1;sh(f<=a)=0; !Moving building shadow
            sh=1
        elsewhere
            sh=0
        end where
        where (tempvegdem>a) !fabovea=tempvegdem>a; !vegdem above DEM
            fabovea=1
        elsewhere
            fabovea=0
        end where
        where (tempvegdem2>a) !gabovea=tempvegdem2>a; !vegdem2 above DEM
            gabovea=1
        elsewhere
            gabovea=0
        end where        
        vegsh2=fabovea-gabovea
        vegsh=max(vegsh,vegsh2)
        where ((vegsh*sh)>0) !vegsh(vegsh.*sh>0)=0;! removing shadows 'behind' buildings
           vegsh=0
        end where
        vbshvegsh=vegsh+vbshvegsh
       
    	! vegsh at high sun altitudes
        if (index==1) then
            firstvegdem=tempvegdem-temp
            where (firstvegdem<=0)!firstvegdem(firstvegdem<=0)=1000;
                firstvegdem=1000
            end where
            where (firstvegdem<dz)!vegsh(firstvegdem<dz)=1;
                vegsh=1
            end where
            tmp=temp*0.0
            where (vegdem2>a)!vegsh=vegsh.*(vegdem2>a);
                tmp=1
            end where
            vegsh=vegsh*tmp
            vbshvegsh=temp*0.0 !vbshvegsh=zeros(sizex,sizey);
        end if
    
    	! Bush shadow on bush plant
        tmp=fabovea*bush
        if ((maxval(bush)>0) .and. (maxval(tmp)>0)) then
            tempbush=temp*0.0
            tempbush(xp1:xp2,yp1:yp2)=bush(xc1:xc2,yc1:yc2)-dz
            g=max(g,tempbush)
            g=bushplant*g
        end if
      
        index=index+1
    END DO

    sh=1-sh
    where (vbshvegsh>0)!vbshvegsh(vbshvegsh>0)=1;
        vbshvegsh=1
    end where
    vbshvegsh=vbshvegsh-vegsh;

    if (maxval(bush)>0) then
        g=g-bush
        where (g>0)!g(g>0)=1;g(g<0)=0;
            g=1
        elsewhere
            g=0
        end where
        vegsh=vegsh-bushplant+g
        where (vegsh<0)!vegsh(vegsh<0)=0;
           vegsh=0
        end where
    end if

    where (vegsh>0)!vegsh(vegsh>0)=1;
        vegsh=1
    end where
    vegsh=1-vegsh
    vbshvegsh=1-vbshvegsh

    !deallocation of grids
    deallocate(f)
    deallocate(temp)
    deallocate(tmp)
    deallocate(stopbuild)
    deallocate(stopveg)
    deallocate(g)
    deallocate(bushplant)
    deallocate(tempvegdem)
    deallocate(tempvegdem2)
    deallocate(fabovea)
    deallocate(gabovea)
    deallocate(firstvegdem)
    deallocate(tempbush)
    deallocate(vegsh2)
    
end subroutine shadowingfunction_20
