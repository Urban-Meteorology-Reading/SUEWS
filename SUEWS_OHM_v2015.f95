!========================================================================================
 SUBROUTINE OHM_v2015(Gridiv)
! Made by HCW Jan 2015 to replace OHMnew (no longer needed).
! Calculates net storage heat flux (QS) from Eq 4, Grimmond et al. 1991, Atm Env.
! Accounts for variable timesteps in dQ*/dt term.
! BareSoilSurfFraction removed so bare soil is now handled like other surfaces.
! Snow part changed from summer wet to winter wet coefficients.
! Changed -333 checks to -999 checks and added error handling
! Gradient now calculated for t-1 (was previously calculated for t-2).
! Modified by HCW 25 Feb 2015
!  Adapted q1,q2,q3 & r1,r2,r3 for multiple grids
! To Do:
!   - Change OHM_TForSummer with latitude (5 degC not always appropriate) - add to inputs
!        - Is OHM_SMForWet = 0.9 appropriate?
!   - No canyons implemented at the moment [OHM_coef(nsurf+1,,)]
! ?? Why is BldgSurf treated differently in terms of SnowFrac? (HCW)
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

  real(kind(1d0)):: dqndt    !Rate of change of net radiation [W m-2 h-1] at t-2
  real(kind(1d0)):: surfrac  !Surface fraction accounting for SnowFrac if appropriate

  !real(kind(1d0)):: OHM_TForSummer = 5  !Use summer coefficients if 5-day Tair >= 5 degC
  real(kind(1d0)):: OHM_TForSummer = 10  !Use summer coefficients if 5-day Tair >= 10 degC - modified for UK HCW 14 Dec 2015
  real(kind(1d0)):: OHM_SMForWet = 0.9  !Use wet coefficients if SM close to soil capacity


  ! OHM coefficients --------
  ! Set to zero initially
  a1=0   ![-]
  a2=0   ![h]
  a3=0   ![W m-2]
  ! -------------------------

  ! Use 5-day running mean Tair to decide whether it is summer or winter ----------------
  if(HDD(id-1,4) >= OHM_TForSummer) then !Summer
     ii=0
  else          !Winter
     ii=2
  endif

  ! Loop through surface types ----------------------------------------------------------
  do is=1,nsurf
     surfrac=sfr(is)
     if(state(is) > 0) then     !Wet surface
           i=ii+1
        else                    !Dry surface
           i=ii+2
           ! If the surface is dry but SM is close to capacity, use coefficients for wet surfaces
           if(is>BldgSurf.and.is/=WaterSurf)then    !Wet soil (i.e. EveTr, DecTr, Grass, BSoil surfaces)
              if(soilmoist(is)/soilstoreCap(is) > OHM_SMForWet) then
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

  ! Calculate radiation part ------------------------------------------------------------
  qs=NAN              !qs  = Net storage heat flux  [W m-2]
  if(qn1>-999) then   !qn1 = Net all-wave radiation [W m-2]
     !if(q1>-999.and.q3>-999) then
        !dqndt = 0.5*(q3-q1)*nsh_real                !gradient at t-2
        dqndt = 0.5*(qn1-q2_grids(Gridiv))*nsh_real   !gradient at t-1

        !Calculate net storage heat flux
        qs = qn1*a1 + dqndt*a2 + a3   !Eq 4, Grimmond et al. 1991
     !endif
     !q1=q2  !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
     !q2=q3  !q2 = net radiation at t-1
     !q3=qn1  !q3 = net radiation at t   (at t-1 when q3 used in next timestep)
     q1_grids(Gridiv) = q2_grids(Gridiv) !q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
     q2_grids(Gridiv) = q3_grids(Gridiv) !q2 = net radiation at t-1
     q3_grids(Gridiv) = qn1              !q3 = net radiation at t (at t-1 when q3 used in next timestep)
  else
     call ErrorHint(21,'Bad value for qn1 found during OHM calculation',qn1,NotUsed,notUsedI)
  endif

  !write(*,*) qs
  !write(*,*) '--------------------'

  ! Do snow calculations separately -----
  ! Added by LJ in August 2013
  if(snowUse==1) then
     deltaQi=NAN
     if(qn1_S>-999) then
        !if(r1>-999.and.r3>-999) then
           !dqndt = 0.5*(r3-r1)*nsh_real    !gradient at t-2
           dqndt = 0.5*(qn1_S-r2_grids(Gridiv))*nsh_real     !gradient at t-1
           ! Calculate net storage heat flux for snow surface (winter wet conditions HCW 15/01/2015)
           deltaQi = qn1_S*OHM_coef(nsurf+2,3,1) + dqndt*OHM_coef(nsurf+2,3,2) + OHM_coef(nsurf+2,3,3)
        !endif
        r1_grids(Gridiv)=r2_grids(Gridiv)
        r2_grids(Gridiv)=r3_grids(Gridiv)
        r3_grids(Gridiv)=qn1_S
     else
     call ErrorHint(21,'Bad value for qn1_S found during OHM calculation',qn1,NotUsed,notUsedI)
     endif
  endif

 return
ENDSUBROUTINE OHM_v2015
!========================================================================================
