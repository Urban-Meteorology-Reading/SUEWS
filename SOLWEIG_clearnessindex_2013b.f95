!> 
subroutine clearnessindex_2013b(zen,jday,Ta,RH,radG,lat,P,I0,CI,Kt,I0et,CIuncorr)
implicit none 
    ! Use somemodule 

    real(kind(1d0))                 :: CI,CIuncorr,I0,I0et,Kt
    integer                         :: jday 
    real(kind(1d0))                 :: P,radG,RH,Ta,zen,lat,iG,Itoa
    real(kind(1d0)), dimension(4)   :: G 
    real(kind(1d0)), parameter          :: pi=3.141592653589793
    
    ! Variable declarations 
    real*8 :: a2      !>  
    !logical :: b      !>  
    real*8 :: b2      !>  
    real*8 :: corr      !>  
    real*8 :: D      !> 
    real*8 :: m      !>  
    real*8 :: Tar      !>  
    real*8 :: Td      !>  
    real*8 :: Trpg      !>  
    real*8 :: Tw      !>  
    real*8 :: u      !>  
    ! 

    ! Clearness Index at the Earth's surface calculated from Crawford and Duchon 1999 

    if (P==-999) then 
        p=1013 !Pressure in millibars
    else
        p=P*10 !Convert from hPa to millibars
    end if 
    Itoa=1370 !Effective solar constant
    !call solar_ESdist(jday,D)
    call sun_distance(jday,D)
    !D=sun_distance(jday) !irradiance differences due to Sun-Earth distances
    m=35.*cos(zen)*((1224.*(cos(zen)**2)+1.)**(-1./2.)) !optical air mass at p=1013
    Trpg=1.021-0.084*(m*(0.000949*p+0.051))**0.5 !Transmission coefficient for Rayliegh scattering and permanent gases

    ! empirical constant depending on latitude 
    if (lat<10) then 
    G=(/3.37,2.85,2.80,2.64/) 
    else if (lat>=10 .and. lat<20) then 
    G=(/2.99,3.02,2.70,2.93/) 
    else if (lat>=20 .and. lat<30) then 
    G=(/3.60,3.00,2.98,2.93/) 
    else if (lat>=30 .and. lat<40) then 
    G=(/3.04,3.11,2.92,2.94/) 
    else if (lat>=40 .and. lat<50) then 
    G=(/2.70,2.95,2.77,2.71/) 
    else if (lat>=50 .and. lat<60) then 
    G=(/2.52,3.07,2.67,2.93/) 
    else if (lat>=60 .and. lat<70) then 
    G=(/1.76,2.69,2.61,2.61/) 
    else if (lat>=70 .and. lat<80) then 
    G=(/1.60,1.67,2.24,2.63/) 
    else if (lat>=80 .and. lat<90) then 
    G=(/1.11,1.44,1.94,2.02/) 
    end if 
    if (jday > 335 .or. jday <= 60) then
        iG=G(1)
    else if (jday > 60 .and. jday <= 152) then
        iG=G(2)
    else if (jday > 152 .and. jday <= 244) then 
        iG=G(3)
    else if (jday > 244 .and. jday <= 335) then 
        iG=G(4)
    end if 
    !dewpoint calculation 
    a2=17.27
    b2=237.7
    Td=(b2*(((a2*Ta)/(b2+Ta))+log(RH)))/(a2-(((a2*Ta)/(b2+Ta))+log(RH)))
    Td=(Td*1.8)+32 !Dewpoint (°F)
    u=exp(0.1133-log(iG+1)+0.0393*Td) !Precipitable water
    Tw=1-0.077*((u*m)**0.3) !Transmission coefficient for water vapor
    Tar=0.935**m !Transmission coefficient for aerosols

    I0=Itoa*cos(zen)*Trpg*Tw*D*Tar
    
    !!! This needs to be checked !!!
    !b=I0==abs(zen)>pi/2
    !I0(b==1)=0
    !clear b
    !if (not(isreal(I0))) then 
    !    I0=0
    !end if 
    
    corr=0.1473*log(90-(zen/pi*180))+0.3454 ! 20070329

    CIuncorr=radG/I0
    CI=CIuncorr+(1-corr)
    I0et=Itoa*cos(zen)*D !extra terrestial solar radiation
    Kt=radG/I0et
    
    !!!! This needs to be checked !!!
    !if (isnan(CI)) then
    !    CI=NaN
    !end if 
    end subroutine clearnessindex_2013b 
    
    !===============================================================================

subroutine sun_distance(jday,D)
 
! Calculates solar irradiance variation based on mean earth sun distance
! with day of year as input.
! Partridge and Platt, 1975

    INTEGER          	::jday
    REAL(KIND(1d0)) 	::b,D

	b = 2*3.141592654*jday/365
	D = sqrt(1.00011 + 0.034221 * cos(b) + 0.001280 * sin(b) + 0.000719 * cos(2*b) + 0.000077 * sin(2*b))

end subroutine sun_distance
