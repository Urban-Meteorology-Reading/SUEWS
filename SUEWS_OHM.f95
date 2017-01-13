!========================================================================================
 SUBROUTINE OHM(Gridiv)
! Made by HCW Jan 2015 to replace OHMnew (no longer needed).
! Calculates net storage heat flux (QS) from Eq 4, Grimmond et al. 1991, Atm Env.
! Accounts for variable timesteps in dQ*/dt term.
! BareSoilSurfFraction removed so bare soil is now handled like other surfaces.
! Snow part changed from summer wet to winter wet coefficients.
! Changed -333 checks to -999 checks and added error handling
! Gradient now calculated for t-1 (was previously calculated for t-2).
! Modified by HCW 25 Feb 2015
!  Adapted q1,q2,q3 & r1,r2,r3 for multiple grids
! Modified by HCW 14 Dec 2016
!  Thresholds for Summer/Winter and Wet/Dry now provided in input files
!  Calculation of dqndt now uses hourly running mean rather than instantaneous values   
! To Do:
!   - No canyons implemented at the moment [OHM_coef(nsurf+1,,)]
!========================================================================================

  use allocateArray
  use data_in
  use defaultNotUsed
  use gis_data
  use sues_data
  use time

  IMPLICIT NONE

  integer:: i,ii
  integer:: Gridiv
  
  real(kind(1d0)):: nsh_nna ! number of timesteps per hour with non -999 values (used for spinup)

  real(kind(1d0)):: dqndt    !Rate of change of net radiation [W m-2 h-1] at t-1
  real(kind(1d0)):: surfrac  !Surface fraction accounting for SnowFrac if appropriate

  real(kind(1d0)):: qn1_av, qn1_S_av    !Average net radiation over previous hour [W m-2]
    
  !These are now provided in SiteInfo (OHMthresh for Summer/Winter and Wet/Dry)
  !!real(kind(1d0)):: OHM_TForSummer = 5  !Use summer coefficients if 5-day Tair >= 5 degC
  !real(kind(1d0)):: OHM_TForSummer = 10  !Use summer coefficients if 5-day Tair >= 10 degC - modified for UK HCW 14 Dec 2015
  !real(kind(1d0)):: OHM_SMForWet = 0.9  !Use wet coefficients if SM close to soil capacity

  
  ! OHM coefficients --------
  ! Set to zero initially
  a1=0   ![-]
  a2=0   ![h]
  a3=0   ![W m-2]
  ! -------------------------

  ! Loop through surface types ----------------------------------------------------------
  do is=1,nsurf
     surfrac=sfr(is)

     ! Use 5-day running mean Tair to decide whether it is summer or winter ----------------
     if(HDD(id-1,4) >= OHM_threshSW(is)) then !Summer
        ii=0
     else          !Winter
        ii=2
     endif

     if(state(is) > 0) then     !Wet surface
           i=ii+1
        else                    !Dry surface
           i=ii+2
           ! If the surface is dry but SM is close to capacity, use coefficients for wet surfaces
           if(is>BldgSurf.and.is/=WaterSurf)then    !Wet soil (i.e. EveTr, DecTr, Grass, BSoil surfaces)
              if(soilmoist(is)/soilstoreCap(is) > OHM_threshWD(is) ) then
                 i=ii+1
              endif
           endif
        endif

     ! If snow, adjust surface fractions accordingly
     if(SnowUse==1.and.is/=BldgSurf.and.is/=WaterSurf) then   ! ?? Why is BldgSurf excluded here?
        surfrac=surfrac*(1-SnowFrac(is))
     endif

     ! Calculate the areally-weighted OHM coefficients
     a1 = a1+surfrac*OHM_coef(is,i,1)
     a2 = a2+surfrac*OHM_coef(is,i,2)
     a3 = a3+surfrac*OHM_coef(is,i,3)

  enddo  !end of loop over surface types ------------------------------------------------

  !write(*,*) '----- OHM coeffs -----'
  !write(*,*) a1,a2,a3

  ! Old OHM calculations (up to v2016a)
  !! Calculate radiation part ------------------------------------------------------------
  !qs=NAN              !qs  = Net storage heat flux  [W m-2]
  !if(qn1>-999) then   !qn1 = Net all-wave radiation [W m-2]
  !   !if(q1>-999.and.q3>-999) then
  !      !dqndt = 0.5*(q3-q1)*nsh_real                !gradient at t-2
  !      dqndt = 0.5*(qn1-q2_grids(Gridiv))*nsh_real   !gradient at t-1
  !
  !      !Calculate net storage heat flux
  !      qs = qn1*a1 + dqndt*a2 + a3   !Eq 4, Grimmond et al. 1991
  !   !endif
  !   !q1=q2  !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
  !   !q2=q3  !q2 = net radiation at t-1
  !   !q3=qn1  !q3 = net radiation at t   (at t-1 when q3 used in next timestep)
  !   q1_grids(Gridiv) = q2_grids(Gridiv) !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
  !   q2_grids(Gridiv) = q3_grids(Gridiv) !q2 = net radiation at t-1
  !   q3_grids(Gridiv) = qn1              !q3 = net radiation at t (at t-1 when q3 used in next timestep)
  !else
  !   call ErrorHint(21,'Bad value for qn1 found during OHM calculation',qn1,NotUsed,notUsedI)
  !endif

  ! New OHM calculations (v2017a onwards) using running mean (HCW Dec 2016)
  ! Calculate radiation part ------------------------------------------------------------
  qs=NAN              !qs  = Net storage heat flux  [W m-2]
  if(qn1>-999) then   !qn1 = Net all-wave radiation [W m-2]   
     ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
     if(nsh > 1) then
        qn1_store(1:(nsh-1),Gridiv) = qn1_store(2:nsh,Gridiv)    
        qn1_store(nsh,Gridiv) = qn1
        nsh_nna = sum(qn1_store(:,Gridiv)/qn1_store(:,Gridiv), mask=qn1_store(:,Gridiv) /= -999) !Find how many are not -999s      
        qn1_av = sum(qn1_store(:,Gridiv), mask=qn1_store(:,Gridiv) /= -999)/nsh_nna
     elseif(nsh==1) then
         qn1_store(:,Gridiv) = qn1
         qn1_av = qn1
     endif         
     ! Store hourly average values (calculated every timestep) for previous 2 hours
     if(nsh > 1) then
        qn1_av_store(1:(2*nsh),Gridiv) = qn1_av_store(2:(2*nsh+1),Gridiv)    
        qn1_av_store(2*nsh+1,Gridiv) = qn1_av
     elseif(nsh==1) then
        qn1_av_store(:,Gridiv) = qn1_av
     endif    
     ! Calculate dQ* per dt for 60 min (using running mean Q* at t hours and (t-2) hours)
     if(any(qn1_av_store == -999)) then
        dqndt=0  ! Set dqndt term to zero for spinup
     else
        dqndt=0.5*(qn1_av_store((2*nsh+1),Gridiv)-qn1_av_store(1,Gridiv))
     endif
     
     ! Calculate net storage heat flux
     qs = qn1*a1 + dqndt*a2 + a3   !Eq 4, Grimmond et al. 1991
     
  else
     call ErrorHint(21,'In SUEWS_OHM.f95: bad value for qn found during qs calculation.',qn1,NotUsed,notUsedI)
  endif

  !write(*,*) qs
  !write(*,*) '--------------------'

  ! Do snow calculations separately -----
  ! Added by LJ in August 2013
  if(snowUse==1) then
     deltaQi=NAN
     if(qn1_S>-999) then
        ! Old OHM calculations (commented out HCW Dec 2016)
        !!if(r1>-999.and.r3>-999) then
        !   !dqndt = 0.5*(r3-r1)*nsh_real    !gradient at t-2
        !   dqndt = 0.5*(qn1_S-r2_grids(Gridiv))*nsh_real     !gradient at t-1
        !   ! Calculate net storage heat flux for snow surface (winter wet conditions HCW 15/01/2015)
        !   deltaQi = qn1_S*OHM_coef(nsurf+2,3,1) + dqndt*OHM_coef(nsurf+2,3,2) + OHM_coef(nsurf+2,3,3)
        !!endif
        !r1_grids(Gridiv)=r2_grids(Gridiv)
        !r2_grids(Gridiv)=r3_grids(Gridiv)
        !r3_grids(Gridiv)=qn1_S
        ! New OHM calculations 
        ! Store instantaneous qn1 values for previous hour (qn1_store) and average (qn1_av)
        if(nsh > 1) then
           qn1_S_store(1:(nsh-1),Gridiv) = qn1_S_store(2:nsh,Gridiv)    
           qn1_S_store(nsh,Gridiv) = qn1_S
           nsh_nna = sum(qn1_S_store(:,Gridiv)/qn1_S_store(:,Gridiv), mask=qn1_S_store(:,Gridiv) /= -999) !Find how many are not -999s      
           qn1_S_av = sum(qn1_S_store(:,Gridiv), mask=qn1_S_store(:,Gridiv) /= -999)/nsh_nna
        elseif(nsh==1) then
           qn1_S_store(:,Gridiv) = qn1_S
           qn1_S_av = qn1_S
        endif         
        ! Store hourly average values (calculated every timestep) for previous 2 hours
        if(nsh > 1) then
           qn1_S_av_store(1:(2*nsh),Gridiv) = qn1_S_av_store(2:(2*nsh+1),Gridiv)    
           qn1_S_av_store(2*nsh+1,Gridiv) = qn1_S_av
        elseif(nsh==1) then
           qn1_S_av_store(:,Gridiv) = qn1_S_av
        endif    
        ! Calculate dQ* per dt for 60 min (using running mean Q* at t hours and (t-2) hours)
        if(any(qn1_S_av_store == -999)) then
           dqndt=0  ! Set dqndt term to zero for spinup
        else
           dqndt=0.5*(qn1_S_av_store((2*nsh+1),Gridiv)-qn1_S_av_store(1,Gridiv))
        endif
         
        ! Calculate net storage heat flux for snow surface (winter wet conditions)
        deltaQi = qn1_S*OHM_coef(nsurf+2,3,1) + dqndt*OHM_coef(nsurf+2,3,2) + OHM_coef(nsurf+2,3,3)   
        
     else
        call ErrorHint(21,'In SUEWS_OHM.f95: bad value for qn(snow) found during qs calculation.',qn1_S,NotUsed,notUsedI)
     endif

  endif

 return
ENDSUBROUTINE OHM
!========================================================================================
