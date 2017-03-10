 subroutine sun_position(year,idectime,UTC,locationlatitude,locationlongitude,locationaltitude,sunazimuth,sunzenith)
    implicit none 
 
    integer :: month,day,hour,min,seas,dayofyear
    REAL(KIND(1D0)) :: sec,year,idectime,UTC
    REAL(KIND(1D0)) :: juliancentury,julianday,julianephemeris_century,julianephemeris_day,&
                       julianephemeris_millenium
    REAL(KIND(1D0)) :: earth_heliocentric_positionlatitude,earth_heliocentric_positionlongitude,&
                       earth_heliocentric_positionradius
    REAL(KIND(1D0)) :: sun_geocentric_positionlatitude, sun_geocentric_positionlongitude    
    REAL(KIND(1D0)) :: nutationlongitude,nutationobliquity    
    REAL(KIND(1D0)) :: corr_obliquity    
    REAL(KIND(1D0)) :: aberration_correction    
    REAL(KIND(1D0)) :: apparent_sun_longitude    
    REAL(KIND(1D0)) :: apparent_stime_at_greenwich    
    REAL(KIND(1D0)) :: sun_rigth_ascension    
    REAL(KIND(1D0)) :: sun_geocentric_declination    
    REAL(KIND(1D0)) :: locationlongitude, observer_local_hour    
    REAL(KIND(1D0)) :: topocentric_sun_positionrigth_ascension ,topocentric_sun_positionrigth_ascension_parallax
    REAL(KIND(1D0)) :: topocentric_sun_positiondeclination,locationlatitude,locationaltitude    
    REAL(KIND(1D0)) :: topocentric_local_hour    
    REAL(KIND(1D0)) :: sunazimuth,sunzenith

    ! This function compute the sun position (zenith and azimuth angle (in degrees) at the observer 
    ! location) as a function of the observer local time and position. 
    !
    ! Input lat and lng should be in degrees, alt in meters. 
    ! 
    ! It is an implementation of the algorithm presented by Reda et Andreas in: 
    ! Reda, I., Andreas, A. (2003) Solar position algorithm for solar 
    ! radiation application. National Renewable Energy Laboratory (NREL) 
    ! Technical report NREL/TP-560-34302. 
    ! This document is avalaible at www.osti.gov/bridge  
    ! Code is translated from matlab code by Fredrik Lindberg (fredrikl@gvc.gu.se)
    ! Last modified: LJ 27 Jan 2016 - Tabs removed


    ! Convert to timevectors from dectime and year
    call dectime_to_timevec(idectime,hour,min,sec)
    dayofyear=floor(idectime)
    call day2month(dayofyear,month,day,seas,year,locationlatitude)

    ! 1. Calculate the Julian Day, and Century. Julian Ephemeris day, century 
    ! and millenium are calculated using a mean delta_t of 33.184 seconds. 
    call julian_calculation(year,month,day,hour,min,sec,UTC,juliancentury,julianday,julianephemeris_century,&
         &julianephemeris_day,julianephemeris_millenium)
    
    ! 2. Calculate the Earth heliocentric longitude, latitude, and radius 
    ! vector (L, B, and R) 
    call earth_heliocentric_position_calculation(julianephemeris_millenium,earth_heliocentric_positionlatitude,&
         &earth_heliocentric_positionlongitude,earth_heliocentric_positionradius)

    ! 3. Calculate the geocentric longitude and latitude 
    call sun_geocentric_position_calculation(earth_heliocentric_positionlongitude,earth_heliocentric_positionlatitude,&
         & sun_geocentric_positionlatitude, sun_geocentric_positionlongitude)

    ! 4. Calculate the nutation in longitude and obliquity (in degrees). 
    call nutation_calculation(julianephemeris_century,nutationlongitude,nutationobliquity)

    ! 5. Calculate the true obliquity of the ecliptic (in degrees). 
    call corr_obliquity_calculation(julianephemeris_millenium, nutationobliquity, corr_obliquity)

    ! 6. Calculate the aberration correction (in degrees) 
    call abberation_correction_calculation(earth_heliocentric_positionradius, aberration_correction)
    
    ! 7. Calculate the apparent sun longitude in degrees) 
    call apparent_sun_longitude_calculation(sun_geocentric_positionlongitude, nutationlongitude,&
         & aberration_correction, apparent_sun_longitude)
    
    ! 8. Calculate the apparent sideral time at Greenwich (in degrees) 
    call apparent_stime_at_greenwich_calculation(julianday,juliancentury, nutationlongitude, &
         &corr_obliquity, apparent_stime_at_greenwich)

    ! 9. Calculate the sun rigth ascension (in degrees) 
    call sun_rigth_ascension_calculation(apparent_sun_longitude, corr_obliquity, sun_geocentric_positionlatitude, &
         &sun_rigth_ascension)

    ! 10. Calculate the geocentric sun declination (in degrees). Positive or 
    ! negative if the sun is north or south of the celestial equator. 
    call sun_geocentric_declination_calculation(apparent_sun_longitude, corr_obliquity, sun_geocentric_positionlatitude, &
         &sun_geocentric_declination)
    
    ! 11. Calculate the observer local hour angle (in degrees, westward from south). 
    call observer_local_hour_calculation(apparent_stime_at_greenwich, locationlongitude, sun_rigth_ascension, observer_local_hour)

    ! 12. Calculate the topocentric sun position (rigth ascension, declination and 
    ! rigth ascension parallax in degrees) 
    call topocentric_sun_position_calculate(topocentric_sun_positionrigth_ascension,&
         &topocentric_sun_positionrigth_ascension_parallax,topocentric_sun_positiondeclination,locationaltitude,&
         &locationlatitude,observer_local_hour,sun_rigth_ascension,sun_geocentric_declination,&
         &earth_heliocentric_positionradius)
    
    ! 13. Calculate the topocentric local hour angle (in degrees) 
    call topocentric_local_hour_calculate(observer_local_hour, topocentric_sun_positionrigth_ascension_parallax,&
         & topocentric_local_hour)       

    ! 14. Calculate the topocentric zenith and azimuth angle (in degrees) 
    call sun_topocentric_zenith_angle_calculate(locationlatitude , topocentric_sun_positiondeclination,&
         & topocentric_local_hour, sunazimuth,sunzenith)


!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Subfunction definitions ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!! 
CONTAINS

subroutine julian_calculation(year,month,day,hour,min,sec,UTC,juliancentury,julianday,julianephemeris_century&
     &,julianephemeris_day,julianephemeris_millenium)
    implicit none 

    REAL(KIND(1D0)) :: A,B,D,delta_t      
    REAL(KIND(1D0)) :: juliancentury      
    REAL(KIND(1D0)) :: julianday      
    REAL(KIND(1D0)) :: julianephemeris_century   
    REAL(KIND(1D0)) :: julianephemeris_day   
    REAL(KIND(1D0)) :: julianephemeris_millenium     
    REAL(KIND(1D0)) :: M,sec,year,UTC 
    integer :: day,hour,min,month
    !REAL(KIND(1D0)) :: time      !>   
    REAL(KIND(1D0)) :: ut_time ,Y !tt,
    ! 
    ! This function compute the julian day and julian century from the local 
    ! time and timezone information. Ephemeris are calculated with a delta_t=0 
    ! seconds. 

    if (month == 1 .or. month == 2) then 
    Y = year - 1.
    M = month + 12
    else
    Y = year
    M = month
    end if 
    ut_time = ((float(hour) - UTC)/24.) + (float(min)/(60.*24.)) + (sec/(60.*60.*24.)) ! time of day in UT time.
    D = day + ut_time ! Day of month in decimal time, ex. 2sd day of month at 12:30:30UT, D=2.521180556

    ! In 1582, the gregorian calendar was adopted 
    if (year == 1582.) then 
    if (month == 10) then 
    if (day <= 4) then ! The Julian calendar ended on October 4, 1582
        B = 0
    else if (day >= 15) then ! The Gregorian calendar started on October 15, 1582
        A = floor(Y/100)
        B = 2 - A + floor(A/4)
    else
        !disp('This date never existed!. Date automatically set to October 4, 1582')
        month = 10
        day = 4
        B = 0
    end if 
    else if (month<10) then ! Julian calendar
        B = 0
    else ! Gregorian calendar
        A = floor(Y/100)
        B = 2 - A + floor(A/4)
    end if 

    else if (year<1582.) then ! Julian calendar
    B = 0
    else
    A = floor(Y/100) ! Gregorian calendar
    B = 2 - A + floor(A/4)
    end if 

    julianday = floor(365.25*(Y+4716.)) + floor(30.6001*(M+1)) + D + B - 1524.5

    delta_t = 0. ! 33.184;
    julianephemeris_day = julianday + (delta_t/86400)

    juliancentury = (julianday - 2451545.) / 36525.

    julianephemeris_century = (julianephemeris_day - 2451545.) / 36525.

    julianephemeris_millenium = julianephemeris_century / 10.
    
end subroutine  julian_calculation 
 
subroutine earth_heliocentric_position_calculation(julianephemeris_millenium,earth_heliocentric_positionlatitude&
     &,earth_heliocentric_positionlongitude,earth_heliocentric_positionradius)
implicit none 

    REAL(KIND(1D0)) :: julianephemeris_millenium      !>  

    REAL(KIND(1D0)),dimension(64) :: A0      !>  
    REAL(KIND(1D0)),dimension(34) :: A1      !>  
    REAL(KIND(1D0)),dimension(20) :: A2      !>  
    REAL(KIND(1D0)),dimension(7) :: A3      !>  
    REAL(KIND(1D0)),dimension(3) :: A4      !>  
    REAL(KIND(1D0)) :: A5      !>  
    REAL(KIND(1D0)),dimension(64) :: B0      !>  
    REAL(KIND(1D0)),dimension(34) :: B1      !>  
    REAL(KIND(1D0)),dimension(20) :: B2      !>  
    REAL(KIND(1D0)),dimension(7) :: B3      !>  
    REAL(KIND(1D0)),dimension(3) :: B4      !>  
    REAL(KIND(1D0)) :: B5      !>  
    REAL(KIND(1D0)),dimension(64) :: C0      !>  
    REAL(KIND(1D0)),dimension(34) :: C1      !>  
    REAL(KIND(1D0)),dimension(20) :: C2      !>  
    REAL(KIND(1D0)),dimension(7) :: C3      !>  
    REAL(KIND(1D0)),dimension(3) :: C4      !>  
    REAL(KIND(1D0)) :: C5     
    REAL(KIND(1D0)),dimension(40) :: A0j      !>  
    REAL(KIND(1D0)),dimension(10) :: A1j     !>  
    REAL(KIND(1D0)),dimension(6) :: A2j      !>  
    REAL(KIND(1D0)),dimension(2) :: A3j      !>  
    REAL(KIND(1D0)) :: A4j      !>  
    REAL(KIND(1D0)),dimension(40) :: B0j      !>   
    REAL(KIND(1D0)),dimension(10) :: B1j      !>   
    REAL(KIND(1D0)),dimension(6) :: B2j      !>  
    REAL(KIND(1D0)),dimension(2) :: B3j      !>  
    REAL(KIND(1D0)) :: B4j      !>  
    REAL(KIND(1D0)),dimension(40) :: C0j      !>  
    REAL(KIND(1D0)),dimension(10) :: C1j      !>  
    REAL(KIND(1D0)),dimension(6) :: C2j      !>  
    REAL(KIND(1D0)),dimension(2) :: C3j      !>  
    REAL(KIND(1D0)) :: C4j      !>  
    REAL(KIND(1D0)),dimension(5) ::  A0i
    REAL(KIND(1D0)),dimension(5) ::  B0i 
    REAL(KIND(1D0)),dimension(5) ::  C0i
    REAL(KIND(1D0)),dimension(2) ::  A1i
    REAL(KIND(1D0)),dimension(2) ::  B1i
    REAL(KIND(1D0)),dimension(2) ::   C1i
    REAL(KIND(1D0)) :: earth_heliocentric_positionlatitude      !>  
    REAL(KIND(1D0)) :: earth_heliocentric_positionlongitude      !>  
    REAL(KIND(1D0)) :: earth_heliocentric_positionradius      !>  
    REAL(KIND(1D0)) :: JME      !>  
    REAL(KIND(1D0)) :: L0      !>  
    !REAL(KIND(1D0)) :: L0_terms      !>  
    REAL(KIND(1D0)) :: L1      !>  
    !REAL(KIND(1D0)) :: L1_terms      !>  
    REAL(KIND(1D0)) :: L2      !>  
    !REAL(KIND(1D0)) :: L2_terms      !>  
    REAL(KIND(1D0)) :: L3      !>  
    !REAL(KIND(1D0)) :: L3_terms      !>  
    REAL(KIND(1D0)) :: L4      !>  
    !REAL(KIND(1D0)) :: L4_terms      !>  
    REAL(KIND(1D0)) :: L5      !>  
    !REAL(KIND(1D0)), dimension(3) :: L5_terms      !>  
    !REAL(KIND(1D0)) :: R0_terms      !>  
    !REAL(KIND(1D0)) :: R1_terms      !>  
    !REAL(KIND(1D0)) :: R2_terms      !>  
    !REAL(KIND(1D0)) :: R3_terms      !>  
    !REAL(KIND(1D0)), dimension(3) :: R4_terms      !>  
    REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0
    
    ! This function compute the earth position relative to the sun, using 
    ! tabulated values. 

    A0=(/175347046,3341656,34894,3497,3418,3136,2676,2343,1324,1273,1199,990,902,857,780,753,505,&
         &492,357,317,284,271,243,206,205,202,156,132,126,115,103,102,102,99,98,86,85,85,80,79,71,&
         &74,74,70,62,61,57,56,56,52,52,51,49,41,41,39,37,37,36,36,33,30,30,25/)
    B0 = (/0.,4.669256800,4.626100,2.744100,2.828900,3.627700,4.418100,6.135200,0.7425000,2.037100,&
         &1.109600,5.233000,2.045000,3.508000,1.179000,2.533000,4.583000,4.205000,2.92,5.849000,&
         &1.899000,0.315,0.345,4.806000,1.869000,2.445800,0.833,3.411000,1.083000,0.645,0.636,0.976,&
         &4.267000,6.21,0.68,5.98,1.3,3.67,1.81,3.04,1.76,3.5,4.68,0.83,3.98,1.82,2.78,4.39,3.47,0.19,&
         &1.33,0.28,0.49,5.37,2.4,6.17,6.04,2.57,1.71,1.78,0.59,0.44,2.74,3.16/)
    C0 = (/0.,6283.075850,12566.15170,5753.384900,3.523100,77713.77150,7860.419400,3930.209700,&
         &11506.76980,529.6910,1577.343500,5884.927,26.29800,398.1490,5223.694,5507.553,&
         &18849.22800,775.5230,0.067,11790.62900,796.2980,10977.07900,5486.778,2544.314,&
         &5573.143,6069.777,213.2990,2942.463,20.77500,0.98,4694.003,15720.83900,7.114000,&
         &2146.170,155.4200,161000.6900,6275.960,71430.70,17260.15,12036.46,5088.630,3154.690,&
         &801.8200,9437.760,8827.390,7084.900,6286.600,14143.50,6279.550,12139.55,1748.020,&
         &5856.480,1194.450,8429.240,19651.05,10447.39,10213.29,1059.380,2352.870,6812.770,&
         &17789.85,83996.85,1349.870,4690.480/)
    A1 = (/628331966747.000,206059.,4303.,425.,119.,109.,93.,72.,68.,67.,59.,56.,45.,36.,29.,21.,19.,19.,17.,16.,&
         &16.,15.,12.,12.,12.,12.,11.,10.,10.,9.,9.,8.,6.,6./)
    B1 =(/ 0.,2.678235,2.635100,1.59,5.796000,2.966000,2.59,1.14,1.87,4.41,2.89,2.17,0.40,0.47,&
         &2.65,5.34,1.85,4.97,2.99,0.030,1.43,1.21,2.83,3.26,5.27,2.08,0.77,1.3,4.24,2.7,5.64,&
         &5.3,2.65,4.67/)
    C1 =(/ 0.,6283.075850,12566.15170,3.523000,26.29800,1577.344,18849.23,529.6900,398.1500,&
         &5507.550,5223.690,155.4200,796.3000,775.5200,7.11,0.98,5486.780,213.3000,6275.960,&
         &2544.310,2146.170,10977.08,1748.020,5088.630,1194.450,4694.,553.5700,3286.600,&
         &1349.870,242.7300,951.7200,2352.870,9437.760,4690.480/)
    A2 =(/ 52919,8720,309,27,16,16,10,9,7,5,4,4,3,3,3,3,3,3,2,2/)
    B2 = (/0.,1.072100,0.867,0.050,5.19,3.68,0.76,2.06,0.83,4.66,1.03,3.44,5.14,6.05,1.19,&
         &6.12,0.31,2.28,4.38,3.75/)
    C2 =(/ 0.,6283.075800,12566.15200,3.52,26.3,155.4200,18849.23,77713.77,775.5200,1577.340,&
         &7.11,5573.140,796.3000,5507.550,242.7300,529.6900,398.1500,553.5700,5223.690,0.98/)
    A3 = (/289,35,17,3,1,1,1/)
    B3 = (/5.8440,0.,5.4900,5.2000,4.7200,5.3000,5.9700/)
    C3 = (/6283.076,0.,12566.15,155.4200,3.52,18849.23,242.7300/)
    A4 = (/114,8,1/)
    B4 = (/3.1420,4.1300,3.8400/) 
    C4 =  (/0.,6283.08,12566.15/)
    A5 =1.
    B5 =3.1400
    C5 =0.
    
    JME = julianephemeris_millenium

    ! Compute the Earth Heliochentric longitude from the tabulated values. 
    L0 = sum(A0 * cos(B0 + (C0 * JME)))
    L1 = sum(A1 * cos(B1 + (C1 * JME)))
    L2 = sum(A2 * cos(B2 + (C2 * JME)))
    L3 = sum(A3 * cos(B3 + (C3 * JME)))
    L4 = sum(A4 * cos(B4 + (C4 * JME)))
    L5 = A5 * cos(B5 + (C5 * JME))

    earth_heliocentric_positionlongitude = &
         &(L0 + (L1 * JME) + (L2 * JME**2) + (L3 * JME**3) + (L4 * JME**4) + (L5 * JME**5)) / 1e8
    ! Convert the longitude to degrees. 
    earth_heliocentric_positionlongitude = earth_heliocentric_positionlongitude * 180./pi
    ! Limit the range to [0,360[; 
    earth_heliocentric_positionlongitude=set_to_range(earth_heliocentric_positionlongitude)

    A0i = (/280,102,80,44,32/)
    B0i = (/3.19900000000000,5.42200000000000,3.88000000000000,3.70000000000000,4./)
    C0i = (/84334.6620000000,5507.55300000000,5223.69000000000,2352.87000000000,1577.34000000000/)
    A1i = (/9,6/)
    B1i =(/3.90000000000000,1.73000000000000/)
    C1i = (/5507.55000000000,5223.69000000000/)

    L0 = sum(A0i * cos(B0i + (C0i * JME)))
    L1 = sum(A1i * cos(B1i + (C1i * JME)))

    earth_heliocentric_positionlatitude = (L0 + (L1 * JME)) / 1e8
    ! Convert the latitude to degrees. 
    earth_heliocentric_positionlatitude = earth_heliocentric_positionlatitude * 180/pi
    ! Limit the range to [0,360]; 
    earth_heliocentric_positionlatitude=set_to_range(earth_heliocentric_positionlatitude)

    A0j = (/100013989,1670700,13956,3084,1628,1576,925,542,472,346,329,307,243,212,186,175,110,&
         &98,86,86,85,63,57,56,49,47,45,43,39,38,37,37,36,35,33,32,32,28,28,26/)
    B0j = (/0.,3.09846350,3.05525000,5.1985,1.1739,2.8469,5.453,4.564,3.661,0.964,5.90,0.299,&
         &4.273,5.847,5.022,3.012,5.055,0.890,5.69,1.27,0.270,0.920,2.01,5.24,3.25,2.58,5.54,&
         &6.01,5.36,2.39,0.830,4.90,1.67,1.84,0.240,0.180,1.78,1.21,1.90,4.59/)
    C0j = (/0.,6283.07585,12566.1517,77713.7715,5753.38490,7860.41940,11506.7700,3930.21000,&
         &5884.92700,5507.55300,5223.69400,5573.14300,11790.6290,1577.34400,10977.0790,18849.2280,&
         &5486.77800,6069.78000,15720.8400,161000.690,17260.1500,529.69,83996.8500,71430.7000,&
         &2544.31000,775.52,9437.76000,6275.96000,4694.,8827.39000,19651.0500,12139.5500,&
         &12036.4600,2942.46000,7084.9,5088.63000,398.15,6286.6,6279.55000,10447.3900/)
    A1j = (/103019,1721,702,32,31,25,18,10,9,9/)
    B1j = (/1.10749000,1.0644,3.142,1.02,2.84,1.32,1.42,5.91,1.42,0.270/)
    C1j = (/6283.07585,12566.1517,0.,18849.2300,5507.55000,5223.69000,1577.34000,10977.0800,&
         &6275.96000,5486.78000/)
    A2j = (/4359,124,12,9,6,3/)
    B2j = (/5.7846,5.579,3.14,3.63,1.87,5.47/)
    C2j = (/6283.07580,12566.1520,0.,77713.7700,5573.14000,18849./)
    A3j = (/145,7/)
    B3j = (/4.273,3.92/)
    C3j = (/6283.07600,12566.1500/)
    A4j = 4
    B4j = 2.56
    C4j = 6283.08000

    ! Compute the Earth heliocentric radius vector 
    L0 = sum(A0j * cos(B0j + (C0j * JME)))
    L1 = sum(A1j * cos(B1j + (C1j * JME)))
    L2 = sum(A2j * cos(B2j + (C2j * JME)))
    L3 = sum(A3j * cos(B3j + (C3j * JME)))
    L4 = A4j * cos(B4j + (C4j * JME))

    ! Units are in AU 
    earth_heliocentric_positionradius = &
         &(L0 + (L1 * JME) + (L2 * JME**2) + (L3 * JME**3) + (L4 * JME**4)) / 1e8

end subroutine  earth_heliocentric_position_calculation 

subroutine sun_geocentric_position_calculation(earth_heliocentric_positionlongitude,&
     &earth_heliocentric_positionlatitude, sun_geocentric_positionlatitude, &
     &sun_geocentric_positionlongitude)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: earth_heliocentric_positionlongitude      !>  
    REAL(KIND(1D0)), intent(in) :: earth_heliocentric_positionlatitude      !> 
    REAL(KIND(1D0)) :: sun_geocentric_positionlatitude      !>  
    REAL(KIND(1D0)) :: sun_geocentric_positionlongitude      !>  
    
    ! This function compute the sun position relative to the earth. 

    sun_geocentric_positionlongitude = earth_heliocentric_positionlongitude + 180.0
    ! Limit the range to [0,360]; 
    sun_geocentric_positionlongitude=set_to_range(sun_geocentric_positionlongitude)

    sun_geocentric_positionlatitude = -earth_heliocentric_positionlatitude
    ! Limit the range to [0,360] 
    sun_geocentric_positionlatitude=set_to_range(sun_geocentric_positionlatitude)
end subroutine  sun_geocentric_position_calculation 

subroutine nutation_calculation(julianephemeris_century,nutationlongitude,nutationobliquity)
    implicit none 
 
    REAL(KIND(1D0)), intent(in) :: julianephemeris_century      !>  
    REAL(KIND(1D0)), dimension(63) :: delta_longitude      !>  
    REAL(KIND(1D0)), dimension(63) :: delta_obliquity      !>  
    REAL(KIND(1D0)) :: JCE      !>  
    REAL(KIND(1D0)) :: nutationlongitude      !>  
    REAL(KIND(1D0)) :: nutationobliquity      !>  
    REAL(KIND(1D0)) , dimension(4) :: p0,p1,p2,p3,p4 
    REAL(KIND(1D0)), dimension(63) ::tabulated_argument      !>  
    REAL(KIND(1D0)) :: X0      !>  
    REAL(KIND(1D0)) :: X1      !>  
    REAL(KIND(1D0)) :: X2      !>  
    REAL(KIND(1D0)) :: X3      !>  
    REAL(KIND(1D0)) :: X4      !>  
    REAL(KIND(1D0)), dimension(5) :: Xi      !>  
    integer, dimension(315) :: Y_terms1     !>  
    integer, dimension(5,63) ::Y_terms
    REAL(KIND(1D0)), dimension(252) :: nutation_terms1     !>  
    REAL(KIND(1D0)), dimension(4,63) ::nutation_terms
    integer :: i
    REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0 
    
    ! This function compute the nutation in longtitude and in obliquity, in 
    ! degrees. 

    ! All Xi are in degrees. 
    JCE = julianephemeris_century

    ! 1. Mean elongation of the moon from the sun 
    p0 = (/ (1/189474.),-0.0019142,445267.11148,297.85036 /) 
    ! X0 = polyval(p, JCE); 
    X0 = p0(1) * JCE**3 + p0(2) * JCE**2 + p0(3) * JCE + p0(4) ! This is faster than polyval...

    ! 2. Mean anomaly of the sun (earth) 
    p1 = (/ -(1/300000.),-0.0001603,35999.05034,357.52772 /)
    ! X1 = polyval(p, JCE); 
    X1 = p1(1) * JCE**3 + p1(2) * JCE**2 + p1(3) * JCE + p1(4)

    ! 3. Mean anomaly of the moon 
    p2 = (/(1/56250.),0.0086972,477198.867398,134.96298 /)
    ! X2 = polyval(p, JCE); 
    X2 = p2(1) * JCE**3 + p2(2) * JCE**2 + p2(3) * JCE + p2(4)

    ! 4. Moon argument of latitude 
    p3 = (/ (1/327270.),-0.0036825,483202.017538,93.27191 /) 
    ! X3 = polyval(p, JCE); 
    X3 = p3(1) * JCE**3 + p3(2) * JCE**2 + p3(3) * JCE + p3(4)

    ! 5. Longitude of the ascending node of the moon's mean orbit on the 
    ! ecliptic, measured from the mean equinox of the date 
    p4 = (/ (1/450000.),0.0020708,-1934.136261,125.04452 /)
    ! X4 = polyval(p, JCE); 
    X4 = p4(1) * JCE**3 + p4(2) * JCE**2 + p4(3) * JCE + p4(4)

    ! Y tabulated terms from the original code 
    Y_terms1 =  (/0,0,0,0,1,-2,0,0,2,2,0,0,0,2,2,0,0,0,0,2,0,1,0,0,0,0,0,1,0,0,-2,1,0,2,2,0,0,0,2,1, &
                0,0,1,2,2,-2,-1,0,2,2,-2,0,1,0,0,-2,0,0,2,1,0,0,-1,2,2,2,0,0,0,0,0,0,1,0,1,2,0,-1,2,2,& 
                0,0,-1,0,1,0,0,1,2,1,-2,0,2,0,0,0,0,-2,2,1,2,0,0,2,2,0,0,2,2,2,0,0,2,0,0,-2,0,1,2,2,0,&
                0,0,2,0,-2,0,0,2,0,0,0,-1,2,1,0,2,0,0,0,2,0,-1,0,1,-2,2,0,2,2,0,1,0,0,1,-2,0,1,0,1,0,-1,0,0,1,0,0,&
                2,-2,0,2,0,-1,2,1,2,0,1,2,2,0,1,0,2,2,-2,1,1,0,0,0,-1,0,2,2,2,0,0,2,1,2,0,1,0,0,-2,0,2,2,2,-2,0,1,&
                2,1,2,0,-2,0,1,2,0,0,0,1,0,-1,1,0,0,-2,-1,0,2,1,-2,0,0,0,1,0,0,2,2,1,-2,0,2,0,1,-2,1,0,2,1,0,0,1,-2,&
                0,-1,0,1,0,0,-2,1,0,0,0,1,0,0,0,0,0,0,1,2,0,0,0,-2,2,2,-1,-1,1,0,0,0,1,1,0,0,0,-1,1,2,2,2,-1,-1,2,2,&
                0,0,3,2,2,2,-1,0,2,2/)
    Y_terms=reshape(Y_terms1,(/5,63/))
    nutation_terms1 = (/-171996.,-174.2,92025.,8.9,-13187.,-1.6,5736.,-3.1,-2274.,-0.2,977.,-0.5,2062.,0.2,-895.,0.5,&
                    1426.,-3.4,54.,-0.1,712.,0.1,-7.,0.,-517.,1.2,224.,-0.6,-386.,-0.4,200.,0.,-301.,0.,129.,-0.1,&
                    217.,-0.5,-95.,0.3,-158.,0.,0.,0.,129.,0.1,-70.,0.,123.,0.,-53.,0.,63.,0.,0.,0.,63.,0.1,-33.,0.,-59.,0.,26.,0.,&
                    -58.,-0.1,32.,0.,-51.,0.,27.,0.,48.,0.,0.,0.,46.,0.,-24.,0.,-38.,0.,16.,0.,-31.,0.,13.,0.,29.,0.,&
                    0.,0.,29.,0.,-12.,0.,26.,0.,0.,0.,-22.,0.,0.,0.,21.,0.,-10.,0.,17.,-0.1,0.,0.,16.,0.,-8.,0.,-16.,0.1,7.,0.,&
                    -15.,0.,9.,0.,-13.,0.,7.,0.,-12.,0.,6.,0.,11.,0.,0.,0.,-10.,0.,5.,0.,-8.,0.,3.,0.,7.,0.,-3.,0.,-7.,0.,0.,0.,&
                    -7.,0.,3.,0.,-7.,0.,3.,0.,6.,0.,0.,0.,6.,0.,-3.,0.,6.,0.,-3.,0.,-6.,0.,3.,0.,-6.,0.,3.,0.,5.,0.,0.,0.,-5.,0.,&
                    3.,0.,-5.,0.,3.,0.,-5.,0.,3.,0.,4.,0.,0.,0.,4.,0.,0.,0.,4.,0.,0.,0.,-4.,0.,0.,0.,-4.,0.,0.,0.,-4.,0.,0.,0.,&
                    3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0.,-3.,0.,0.,0./)
    nutation_terms=reshape(nutation_terms1,(/4,63/))
    ! Using the tabulated values, compute the delta_longitude and 
    ! delta_obliquity. 
    Xi = (/X0, X1, X2, X3, X4/)
    
    do i=1,63
        tabulated_argument(i)=&
             &((Y_terms(1,i)*Xi(1))+(Y_terms(2,i)*Xi(2))+(Y_terms(3,i)*Xi(3))+(Y_terms(4,i)*Xi(4))+(Y_terms(5,i)*Xi(5)))*pi/180
    end do
    
    delta_longitude = ((nutation_terms(1,:) + (nutation_terms(2,:) * JCE))) * sin(tabulated_argument)
    delta_obliquity = ((nutation_terms(3,:) + (nutation_terms(4,:) * JCE))) * cos(tabulated_argument)

    ! Nutation in longitude 
    nutationlongitude = sum(delta_longitude) / 36000000.0

    ! Nutation in obliquity 
    nutationobliquity = sum(delta_obliquity) / 36000000.0

end subroutine  nutation_calculation 

subroutine corr_obliquity_calculation(julianephemeris_millenium, nutationobliquity, corr_obliquity)
    implicit none 

    REAL(KIND(1D0)), intent(out) :: corr_obliquity      !>  
    REAL(KIND(1D0)), intent(in) :: julianephemeris_millenium     !>  
    REAL(KIND(1D0)), intent(in) :: nutationobliquity     !>  
    REAL(KIND(1D0)) :: mean_obliquity      !>  
    REAL(KIND(1D0)), dimension(11) :: p      !>  
    REAL(KIND(1D0)) :: U      !>  

    ! This function compute the true obliquity of the ecliptic. 


    p = (/ 2.45,5.79,27.87,7.12,-39.05,-249.67,-51.38,1999.25,-1.55,-4680.93,84381.448 /) 
    ! mean_obliquity = polyval(p, julian.ephemeris_millenium/10); 

    U = julianephemeris_millenium/10
    mean_obliquity =&
         &p(1)*U**10 + p(2)*U**9 + p(3)*U**8 + p(4)*U**7 + p(5)*U**6 + p(6)*U**5 + p(7)*U**4 + &
         &p(8)*U**3 + p(9)*U**2 + p(10)*U + p(11)

    corr_obliquity = (mean_obliquity/3600) + nutationobliquity
end subroutine  corr_obliquity_calculation 
 
subroutine abberation_correction_calculation(earth_heliocentric_positionradius, aberration_correction)
    implicit none 

    REAL(KIND(1D0)), intent(out) :: aberration_correction      !>  
    REAL(KIND(1D0)), intent(in) :: earth_heliocentric_positionradius     !>  
 
    ! This function compute the aberration_correction, as a function of the 
    ! earth-sun distance. 

    aberration_correction = -20.4898/(3600*earth_heliocentric_positionradius)

end subroutine  abberation_correction_calculation 

subroutine apparent_sun_longitude_calculation(sun_geocentric_positionlongitude, nutationlongitude,&
     & aberration_correction, apparent_sun_longitude)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: aberration_correction      !>  
    REAL(KIND(1D0)), intent(out) :: apparent_sun_longitude      !>  
    REAL(KIND(1D0)), intent(in) :: nutationlongitude      !>  
    REAL(KIND(1D0)), intent(in) :: sun_geocentric_positionlongitude      !>  

    ! This function compute the sun apparent longitude 

    apparent_sun_longitude = sun_geocentric_positionlongitude + nutationlongitude + aberration_correction
    
end subroutine  apparent_sun_longitude_calculation 

subroutine apparent_stime_at_greenwich_calculation(julianday,juliancentury, nutationlongitude,&
     & corr_obliquity, apparent_stime_at_greenwich)
    implicit none 

    REAL(KIND(1D0)), intent(out) :: apparent_stime_at_greenwich      !>  
    REAL(KIND(1D0)), intent(in) :: corr_obliquity      !>  
    REAL(KIND(1D0)), intent(in) :: julianday     !> 
    REAL(KIND(1D0)), intent(in) :: juliancentury     !>  
    REAL(KIND(1D0)), intent(in) :: nutationlongitude      !>  
    REAL(KIND(1D0)) :: JC      !>  
    REAL(KIND(1D0)) :: JD      !>  
    REAL(KIND(1D0)) :: mean_stime      !>  
    REAL(KIND(1D0)),PARAMETER       :: pi=3.14159265358979d+0 
    
    ! This function compute the apparent sideral time at Greenwich. 

    JD = julianday
    JC = juliancentury

    ! Mean sideral time, in degrees 
    mean_stime = 280.46061837d+0 + (360.98564736629d+0*(JD-2451545.0d+0)) + (0.000387933d+0*JC**2) - (JC**3/38710000.0d+0)

    ! Limit the range to [0-360]; 
    mean_stime=set_to_range(mean_stime)

    apparent_stime_at_greenwich = mean_stime + (nutationlongitude * cos(corr_obliquity * pi/180))
end subroutine  apparent_stime_at_greenwich_calculation 

subroutine sun_rigth_ascension_calculation(apparent_sun_longitude, corr_obliquity, &
     &sun_geocentric_positionlatitude, sun_rigth_ascension)
    implicit none 
 
    REAL(KIND(1D0)), intent(in) :: apparent_sun_longitude      !>  
    REAL(KIND(1D0)), intent(in) :: corr_obliquity      !>  
    REAL(KIND(1D0)), intent(in) :: sun_geocentric_positionlatitude      !>  
    REAL(KIND(1D0)), intent(out) :: sun_rigth_ascension      !>  
    REAL(KIND(1D0)) :: argument_denominator      !>  
    REAL(KIND(1D0)) :: argument_numerator      !>  
    REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0 
    
    ! This function compute the sun rigth ascension. 

    argument_numerator = (sin(apparent_sun_longitude * pi/180.0) * cos(corr_obliquity * pi/180.0)) - &
    (tan(sun_geocentric_positionlatitude * pi/180.0) * sin(corr_obliquity * pi/180.0))
    argument_denominator = cos(apparent_sun_longitude * pi/180.0)

    sun_rigth_ascension = atan2(argument_numerator, argument_denominator) * 180.0/pi
    ! Limit the range to [0,360]; 
    sun_rigth_ascension=set_to_range(sun_rigth_ascension)
end subroutine  sun_rigth_ascension_calculation 

subroutine sun_geocentric_declination_calculation(apparent_sun_longitude, corr_obliquity, &
     &sun_geocentric_positionlatitude, sun_geocentric_declination)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: apparent_sun_longitude      !>  
    REAL(KIND(1D0)), intent(in) :: corr_obliquity      !>  
    REAL(KIND(1D0)), intent(out) :: sun_geocentric_declination      !>  
    REAL(KIND(1D0)), intent(in) :: sun_geocentric_positionlatitude     !>  
    REAL(KIND(1D0)) :: argument      !>  
    REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0  

    argument = (sin(sun_geocentric_positionlatitude * pi/180.0) * cos(corr_obliquity * pi/180.0)) + &
    (cos(sun_geocentric_positionlatitude * pi/180.0) * sin(corr_obliquity * pi/180) * sin(apparent_sun_longitude * pi/180.0))

    sun_geocentric_declination = asin(argument) * 180.0/pi
end subroutine  sun_geocentric_declination_calculation 

subroutine observer_local_hour_calculation(apparent_stime_at_greenwich, locationlongitude, &
     &sun_rigth_ascension, observer_local_hour)
    implicit none 
 
    REAL(KIND(1D0)), intent(in) :: apparent_stime_at_greenwich      !>  
    REAL(KIND(1D0)), intent(in) :: locationlongitude     !>  
    REAL(KIND(1D0)), intent(out) :: observer_local_hour      !>  
    REAL(KIND(1D0)), intent(in) :: sun_rigth_ascension      !>  


    observer_local_hour = apparent_stime_at_greenwich + locationlongitude - sun_rigth_ascension
    ! Set the range to [0-360] 
    observer_local_hour=set_to_range(observer_local_hour)
end subroutine  observer_local_hour_calculation 

subroutine topocentric_sun_position_calculate(topocentric_sun_positionrigth_ascension &
     &,topocentric_sun_positionrigth_ascension_parallax,topocentric_sun_positiondeclination,&
     &locationaltitude,locationlatitude,observer_local_hour,sun_rigth_ascension,&
     &sun_geocentric_declination,earth_heliocentric_positionradius)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: earth_heliocentric_positionradius   
    REAL(KIND(1D0)), intent(in) :: locationlatitude      !>  
    REAL(KIND(1D0)), intent(in) :: locationaltitude
    REAL(KIND(1D0)), intent(in) :: observer_local_hour      !>  
    REAL(KIND(1D0)),  intent(in) :: sun_geocentric_declination      !>  
    REAL(KIND(1D0)),  intent(in) :: sun_rigth_ascension      !>  
    REAL(KIND(1D0)) :: denominator      !>  
    REAL(KIND(1D0)) :: eq_horizontal_parallax      !>  
    REAL(KIND(1D0)) :: nominator      !>  
    REAL(KIND(1D0)) :: sun_rigth_ascension_parallax      !>  
    REAL(KIND(1D0)) :: topocentric_sun_positiondeclination      !>  
    REAL(KIND(1D0)) :: topocentric_sun_positionrigth_ascension      !>  
    REAL(KIND(1D0)) :: topocentric_sun_positionrigth_ascension_parallax      !>  
    REAL(KIND(1D0)) :: u      !>  
    REAL(KIND(1D0)) :: x      !>  
    REAL(KIND(1D0)) :: y      !>
    REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0
    
    ! topocentric_sun_positionrigth_ascension_parallax
    ! This function compute the sun position (rigth ascension and declination) 
    ! with respect to the observer local position at the Earth surface. 

    ! Equatorial horizontal parallax of the sun in degrees 
    eq_horizontal_parallax = 8.794 / (3600 * earth_heliocentric_positionradius)

    ! Term u, used in the following calculations (in radians) 
    u = atan(0.99664719 * tan(locationlatitude * pi/180))

    ! Term x, used in the following calculations 
    x = cos(u) + ((locationaltitude/6378140) * cos(locationaltitude * pi/180))

    ! Term y, used in the following calculations 
    y = (0.99664719d+0 * sin(u)) + ((locationaltitude/6378140) * sin(locationlatitude * pi/180))

    ! Parallax in the sun rigth ascension (in radians) 
    nominator = -x * sin(eq_horizontal_parallax * pi/180.0) * sin(observer_local_hour * pi/180.0)
    denominator = cos(sun_geocentric_declination * pi/180.0) - &
    (x * sin(eq_horizontal_parallax * pi/180.0) * cos(observer_local_hour * pi/180.0))
    sun_rigth_ascension_parallax = atan2(nominator, denominator)
    ! Conversion to degrees. 
    topocentric_sun_positionrigth_ascension_parallax = sun_rigth_ascension_parallax * 180.0/pi

    ! Topocentric sun rigth ascension (in degrees) 
    topocentric_sun_positionrigth_ascension = sun_rigth_ascension + (sun_rigth_ascension_parallax * 180.0/pi)

    ! Topocentric sun declination (in degrees) 
    nominator = (sin(sun_geocentric_declination * pi/180.0) - (y*sin(eq_horizontal_parallax * pi/180.0)))&
         & * cos(sun_rigth_ascension_parallax)
    denominator = cos(sun_geocentric_declination * pi/180.0) - (y*sin(eq_horizontal_parallax * pi/180.0))&
         & * cos(observer_local_hour * pi/180.0)
    topocentric_sun_positiondeclination = atan2(nominator, denominator) * 180.0/pi
end subroutine  topocentric_sun_position_calculate 

subroutine topocentric_local_hour_calculate(observer_local_hour, topocentric_sun_positionrigth_ascension_parallax,&
     & topocentric_local_hour)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: observer_local_hour      !>  
    REAL(KIND(1D0)), intent(out) :: topocentric_local_hour      !>  
    REAL(KIND(1D0)), intent(in) :: topocentric_sun_positionrigth_ascension_parallax     !>  

    ! This function compute the topocentric local jour angle in degrees 

    topocentric_local_hour = observer_local_hour - topocentric_sun_positionrigth_ascension_parallax
end subroutine  topocentric_local_hour_calculate 
 
subroutine sun_topocentric_zenith_angle_calculate(locationlatitude , topocentric_sun_positiondeclination, &
     &topocentric_local_hour, sunazimuth,sunzenith)
    implicit none 

    REAL(KIND(1D0)), intent(in) :: locationlatitude     !>  
    REAL(KIND(1D0)), intent(in) :: topocentric_local_hour      !>  
    REAL(KIND(1D0)), intent(in) :: topocentric_sun_positiondeclination   
    REAL(KIND(1D0)) :: corr_elevation      !>  
    REAL(KIND(1D0)) :: apparent_elevation      !>  
    REAL(KIND(1D0)) :: argument      !>  
    REAL(KIND(1D0)) :: denominator      !>  
    REAL(KIND(1D0)) :: nominator      !>  
    REAL(KIND(1D0)) :: refraction_corr      !>  
    REAL(KIND(1D0)) :: sunazimuth      !>  
    REAL(KIND(1D0)) :: sunzenith      !>  
     REAL(KIND(1D0)),PARAMETER       :: pi=3.141592653589793d+0
    ! This function compute the sun zenith angle, taking into account the 
    ! atmospheric refraction. A default temperature of 283K and a 
    ! default pressure of 1010 mbar are used. 

    ! Topocentric elevation, without atmospheric refraction 
    argument = (sin(locationlatitude * pi/180.0) * sin(topocentric_sun_positiondeclination * pi/180.0)) + &
    (cos(locationlatitude * pi/180.0) * cos(topocentric_sun_positiondeclination * pi/180.0) * &
         &cos(topocentric_local_hour * pi/180.0))
    corr_elevation = asin(argument) * 180.0/pi

    ! Atmospheric refraction correction (in degrees) 
    argument = corr_elevation + (10.3/(corr_elevation + 5.11))
    refraction_corr = 1.02 / (60 * tan(argument * pi/180.0))

    ! For exact pressure and temperature correction, use this, 
    ! with P the pressure in mbar amd T the temperature in Kelvins: 
    ! refraction_corr = (P/1010) * (283/T) * 1.02 / (60 * tan(argument * pi/180)); 

    ! Apparent elevation 
    apparent_elevation = corr_elevation + refraction_corr

    sunzenith = 90.0 - apparent_elevation

    ! Topocentric azimuth angle. The +180 conversion is to pass from astronomer 
    ! notation (westward from south) to navigation notation (eastward from 
    ! north); 
    nominator = sin(topocentric_local_hour * pi/180.0)
    denominator = (cos(topocentric_local_hour * pi/180.0) * sin(locationlatitude * pi/180.0)) - &
    (tan(topocentric_sun_positiondeclination * pi/180.0) * cos(locationlatitude * pi/180.0))
    sunazimuth = (atan2(nominator, denominator) * 180.0/pi) + 180.0
    ! Set the range to [0-360] 
    sunazimuth=set_to_range(sunazimuth)

end subroutine  sun_topocentric_zenith_angle_calculate 

FUNCTION set_to_range(var) RESULT(vari)
! This function make sure the variable is in the specified range.     

    REAL(KIND(1D0)) :: max_interval      !>  
    REAL(KIND(1D0)) :: min_interval      !>  
    REAL(KIND(1D0)) :: var 
    REAL(KIND(1D0)) :: vari 
    ! 
    max_interval=360.0
    min_interval=0.0

    vari = var - max_interval * floor(var/max_interval)

    if (vari<min_interval) then 
    vari = vari + max_interval
    end if 

END FUNCTION  set_to_range 

end subroutine  sun_position 


