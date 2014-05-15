! This is the core function of the SOLWEIG model 
! 2013-10-27 
! Fredrik Lindberg, fredrikl@gvc.gu.se 
! Göteborg Urban Climate Group 
! Gothenburg University

subroutine Solweig_2014a_core(scale,x,y,albedo_b,albedo_g,absK,absL,ewall,eground,Fside,Fup,lat,lng,alt,timezone, &
    year,DOY,hour,dectime,Ta,RH,P,radG,radD,radI,usevegdem,onlyglobal,height,amaxvalue,psi,ith,timestepdec)

use matsize
 
    implicit none
    integer         :: usevegdem,onlyglobal,DOY,hour,first,second,i,j,x,y,dfm,ith
    real(kind(1d0)) :: absK,absL,albedo_b,albedo_g,eground,ewall,Fside,Fup 
    real(kind(1d0)) :: t,Tstart,lat,lng,alt,timezone,height,amaxvalue,psi
    real(kind(1d0)) :: scale,azimuth,altitude,zen,zenith
    real(kind(1d0)) :: CI,CI_Tg,c,I0,Kt,Tg,Tgamp,Tw,Ktc,weight1
    real(kind(1d0)) :: dectime,Ta,RH,P,radG,radD,radI,radI0,idectime,tdectime
    real(kind(1d0)) :: corr,I0et,CIuncorr,s
    real(kind(1d0)) :: tTa,tP,tRH,tradG,thour    
	real(kind(1d0)) :: SNDN,SNUP,DEC,DAYL,YEAR,timestepdec
    ! Internal grids
    real(kind(1d0)),allocatable,dimension(:,:) :: tmp,Knight,svfbuveg,Tgmap0!,Tgmap
    !Search directions for Ground View Factors (GVF)
    real(kind(1d0)) :: azimuthA(1:18)=(/ (i*(360.0/18.0),i=0,17) /)   
    ! temporary parameters and variables for testing
    real(kind(1d0)),parameter   :: pi=3.141592653589793
    real(kind(1d0)),parameter   :: SBC=5.67051e-8
    REAL(KIND(1D0)),PARAMETER   :: DEG2RAD=0.017453292,RAD2DEG=57.29577951
    REAL(KIND(1D0))             :: msteg,esky,ea

    !!!!!! Begin program !!!!!!
    ! internal grids
    allocate(tmp(sizey,sizex))
    allocate(Knight(sizey,sizex))
    allocate(Tgmap0(sizey,sizex))
    allocate(svfbuveg(sizey,sizex))
    !allocate(Tgmap0(sizey,sizex))
    
    ! external grids
    if (allocated(Kdown)) deallocate(Kdown); allocate(Kdown(sizey,sizex))
    if (allocated(Kup)) deallocate(Kup); allocate(Kup(sizey,sizex))
    if (allocated(Knorth)) deallocate(Knorth); allocate(Knorth(sizey,sizex))
    if (allocated(Kwest)) deallocate(Kwest); allocate(Kwest(sizey,sizex))
    if (allocated(Ksouth)) deallocate(Ksouth); allocate(Ksouth(sizey,sizex))
    if (allocated(Keast)) deallocate(Keast); allocate(Keast(sizey,sizex))
    if (allocated(Ldown)) deallocate(Ldown); allocate(Ldown(sizey,sizex))
    if (allocated(Lup)) deallocate(Lup); allocate(Lup(sizey,sizex))
    if (allocated(LNorth)) deallocate(LNorth); allocate(LNorth(sizey,sizex))
    if (allocated(Lwest)) deallocate(Lwest); allocate(Lwest(sizey,sizex))
    if (allocated(Lsouth)) deallocate(Lsouth); allocate(Lsouth(sizey,sizex))
    if (allocated(Least)) deallocate(Least); allocate(Least(sizey,sizex))
    if (allocated(gvf)) deallocate(gvf); allocate(gvf(sizey,sizex))
    if (allocated(Sstr)) deallocate(Sstr); allocate(Sstr(sizey,sizex))
    if (allocated(Tmrt)) deallocate(Tmrt); allocate(Tmrt(sizey,sizex))
    if (allocated(shadow)) deallocate(shadow); allocate(shadow(sizey,sizex))
    if (allocated(sos)) deallocate(sos); allocate(sos(sizey,sizex))
    if (allocated(F_sh)) deallocate(F_sh); allocate(F_sh(sizey,sizex))
    if (allocated(svfalfa)) deallocate(svfalfa); allocate(svfalfa(sizey,sizex))
        
    ! Radiation sensor setup offset in degrees 
    t=0

    !Surface temperature difference at sunrise 
    Tstart=3.41

    !Initialization of maps 
    Knight=0.0
    
    tmp=1-(svf+svfveg-1) 
    where (tmp<=0) tmp=0.000000001 ! avoiding log(0)
    svfalfa=asin(exp(log(tmp)/2))

    !Parameterisarion for Lup 
    first=anint(height) !Radiative surface influence, Rule of thumb by Schmid et al. (1990).
    if (first==0) then 
        first=1
    end if 
    second=anint(height*20)
    
    ! SVF combines for buildings and vegetation    
    svfbuveg=(svf-(1-svfveg)*(1-psi))
    
	! Sun position related things
    idectime=dectime-timestepdec/2! sun position at middle of timestep before
    call sun_position(year,idectime,timezone,lat,lng,alt,azimuth,zenith)
    call DAYLEN(DOY,LAT*RAD2DEG,DAYL,DEC,SNDN,SNUP)
    zen=zenith*DEG2RAD
    altitude=90-zenith
    
    !Determination of clear-sky emissivity from Prata (1996), (This part might be calculated outside (NARP) later) 
    ea=6.107*10**((7.5*Ta)/(237.3+Ta))*(RH/100)!Vapor pressure
    msteg=46.5*(ea/(Ta+273.15))
    esky=(1-(1+msteg)*exp(-((1.2+3.0*msteg)**0.5)))-0.04
    
    
    !!! DAYTIME !!!
    if (altitude>0) then 
    
        !Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction 
        !factor for low sun elevations after Lindberg et al. (2008) 
        call clearnessindex_2013b(zen,DOY,Ta,RH/100,radG,lat,P,I0,CI,Kt,I0et,CIuncorr) 
        if (CI>1) CI=1  !!FIX THIS?? .and. CI<inf) CI=1
        CIlatenight=CI
    
        !Estimation of radD and radI if not measured after Reindl et al. (1990) 
        if (onlyglobal == 1) then 
            call diffusefraction(radG,altitude,Kt,Ta,RH,radI,radD) 
        end if 
    
        !Shadow images 
        if (usevegdem==1) then ! with vegetation
            call shadowingfunction_20(azimuth,altitude,scale,amaxvalue) 
            shadow=(sh-(1-vegsh)*(1-psi))
        else ! without vegetation
            call shadowingfunction_10(azimuth,altitude,scale)
            vegsh = 1.0D0 
            shadow=sh
        end if 
    
        !Ground View Factors based on shadowpatterns and sunlit walls 
        gvf=0.0D0
        call wallinsun_veg(azimuth) 
        do j=1,size(azimuthA) 
            call sunonsurface_veg(azimuthA(j),scale,first,second,psi) 
            gvf=gvf+sos
        end do 
        gvf=gvf/size(azimuthA)+(buildings*(-1)+1)
        
        !Building height angle from svf 
        call cylindric_wedge(zen) !Fraction shadow on building walls based on sun altitude and svf
        !F_sh(isnan(F_sh))=0.5 FIXTHIS
        
        !!! Calculation of shortwave daytime radiative fluxes !!! 
        Kdown=radI*shadow*sin(altitude*DEG2RAD)+radD*svfbuveg+ &
            radG*albedo_b*(1-svfbuveg)*(1-F_sh) !*sin(altitude(i)*(pi/180));
        
        Kup=albedo_g*(radI*gvf*sin(altitude*DEG2RAD)+radD*svfbuveg+ &
            radG*albedo_b*(1-svfbuveg)*(1-F_sh))
        
        call Kside_veg_v24(radI,radG,radD,azimuth,altitude,psi,t,albedo_b)
        
        !!!! Surface temperature parameterisation during daytime !!!! 
        dfm=abs(172-DOY) !Day from midsommer
        Tgamp=0.000006*dfm**3-0.0017*dfm**2+0.0127*dfm+17.084+Tstart !sinus function for daily surface temperature wave
        !Tg=Tgamp*sin(((hour-rise)/(15-rise))*pi/2)-Tstart ! check if this should be 15 (3 pm)
        Tg=Tgamp*sin(((dectime-DOY-SNUP/24)/(15/24-SNUP/24))*pi/2)-Tstart !new sunrise time 2014a
        if (Tg<0) Tg=0! temporary for removing low Tg during morning 20140513, from SOLWEIG1D
        
        !New estimation of Tg reduction for non-clear situation based on Reindl et al. 1990 
        Ktc=1.0
        call diffusefraction(I0,altitude,Ktc,Ta,RH,radI0,s) 
        corr=0.1473*log(90.-(zen*RAD2DEG))+0.3454 ! 20070329 temporary correction of latitude from Lindberg et al. 2008
        CI_Tg=(radI/radI0)+(1.-corr)
        
        if (CI_Tg>1) CI_Tg=1  !!FIX THIS?? .and. CI_Tg<inf then CI_Tg=1
        
        Tg=Tg*CI_Tg !new estimation
        Tw=Tg
        
        !!!! Lup, daytime !!!! 
        !Surface temperature wave delay 
        !if (hour>rise) then 
        !    Tgmapgvf1=((gvf-Tgmapgvf)*0.75*(-1.0))*Tg
         !   Tgmap=gvf*Tg+Ta
        !    Tgmap=Tgmap+Tgmapgvf1
        !    Lup=SBC*eground*((Tgmap+273.15)**4)
        !    Tgmapgvf=gvf
        !else
        !    Lup=SBC*eground*((gvf*Tg+Ta+273.15)**4)
        !    Tgmapgvf=gvf
        !end if 
        
        !Surface temperature wave delay - new as from 2014a
        Tgmap0=gvf*Tg+Ta ! current timestep
        if (firstdaytime==1) then !"first in morning"
            Tgmap1=Tgmap0
        end if
        if (timeadd>=(59/1440)) then !more or equal to 59 min
            weight1=exp(-33.27*timeadd) !surface temperature delay function - 1 step
            Tgmap1=Tgmap0*(1-weight1)+Tgmap1*weight1
            Lup=SBC*eground*((Tgmap1+273.15)**4)
            if (timestepdec>(59/1440)) then
                timeadd=timestepdec
            else
                timeadd=0
            end if
         else
            timeadd=timeadd+timestepdec
            weight1=exp(-33.27*timeadd) !surface temperature delay function - 1 step
            Lup=SBC*eground*((Tgmap0*(1-weight1)+Tgmap1*weight1+273.15)**4)
         end if
         firstdaytime=0
        
    else !!!!!!! NIGHTTIME !!!!!!!!
    
        !Nocturnal cloudfraction from Offerle et al. 2003 
        if (dectime<(DOY+0.5) .and. dectime>DOY .and. altitude<0) then 
            i=0
            do while (tdectime<(DOY+SNUP/24))
                READ(67,*) s,thour,tdectime,s,s,s,s,s,s,tRH,&
                       tTa,tP,s,tradG,s,s,s,&
                       s,s,s,s,s
                i=i+1
            end do
            do j=1,i
                backspace 67
            end do
        call sun_position(year,tdectime,timezone,lat,lng,alt,s,zenith)
        zen=zenith*DEG2RAD
        call clearnessindex_2013b(zen,DOY,tTa,tRH/100,tradG,lat,tP,I0,CI,Kt,I0et,CIuncorr) 
        else
            if (ith==0) then
                CI=1.0
            else
                CI=CIlatenight
            end if
        end if   

        Tw=0.0
        Tg=0.0
    
        !Nocturnal Kfluxes set to 0 
        Kdown=0.0
        Kwest=0.0
        Kup=0.0
        Keast=0.0
        Ksouth=0.0
        Knorth=0.0
        shadow=0.0
    
        !!! Lup !!! 
        Lup=SBC*eground*((Knight+Ta+Tg+273.15)**4)
        firstdaytime=1
    end if 
    
    !!! Ldown !!! 
    Ldown=(svf+svfveg-1)*esky*SBC*((Ta+273.15)**4)+(2-svfveg-svfaveg)*ewall*SBC*((Ta+273.15)**4)+ &
        (svfaveg-svf)*ewall*SBC*((Ta+273.15+Tw)**4)+(2-svf-svfveg)*(1-ewall)*esky*SBC*((Ta+273.15)**4) !Jonsson et al. (2006)
    Ldown=Ldown-25 ! Shown by Jonsson et al. (2006) and Duarte et al. (2006)
    if (CI>1) CI=1  !!FIX THIS?? .and. CI<inf) CI=1
	if (CI < 0.95) then !non-clear conditions
        c=1-CI
        Ldown=Ldown*(1-c)+c*SBC*((Ta+273.15)**4)
    end if 
    
    !!! Lside !!!
    call Lside_veg_v2(azimuth,altitude,Ta,Tw,SBC,ewall,esky,t)
    
    !!! Calculation of radiant flux density and Tmrt !!! 
    Sstr=absK*(Kdown*Fup+Kup*Fup+Knorth*Fside+Keast*Fside+Ksouth*Fside+Kwest*Fside) &
        +absL*(Ldown*Fup+Lup*Fup+Lnorth*Fside+Least*Fside+Lsouth*Fside+Lwest*Fside)
    Tmrt=sqrt(sqrt((Sstr/(absL*SBC))))-273.2
    
	!write(*,*) esky
	!write(*,*) Tmrt(y,x)
	!write(*,*) Kdown(y,x)
    !write(*,*) Kup(y,x)
    !write(*,*) Ksouth(y,x)
    !write(*,*) Kwest(y,x)
    !write(*,*) Knorth(y,x)
    !write(*,*) Keast(y,x)
	!write(*,*) Ldown(y,x)
    !write(*,*) Lup(y,x)
    !write(*,*) Lsouth(y,x)
    !write(*,*) Lwest(y,x)
    !write(*,*) Lnorth(y,x)
    !write(*,*) Least(y,x)
    !write(*,*) 
    
    deallocate(tmp)
    deallocate(Knight)
    deallocate(Tgmap0)
    deallocate(svfbuveg)
        
    ! external grids
    deallocate(Kdown)
    deallocate(Kup)
    deallocate(Knorth)
    deallocate(Kwest)
    deallocate(Ksouth)
    deallocate(Keast)
    deallocate(Ldown)
    deallocate(Lup)
    deallocate(Lnorth)
    deallocate(Lwest)
    deallocate(Lsouth)
    deallocate(Least)
    
end subroutine Solweig_2014a_core 

