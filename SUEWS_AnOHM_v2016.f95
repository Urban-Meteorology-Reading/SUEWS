!========================================================================================
subroutine AnOHM_v2016(Gridiv)
! author: Ting Sun
! 
! purpose:
! calculate heat storage based within AnOHM framework.
! 
! input:
! 1) Gridiv: grid number, a global variable. 
! with Gridiv, required met forcings and sfc. properties can be loaded for calculation.
! 
! output:
! 1) grid ensemble heat storage.
! QS = a1*(Q*)+a2*(dQ*/dt)+a3
! 
! ref:
! the AnOHM paper to be added.
! 
! history:
! 20160301: initial version
!========================================================================================

	use allocateArray
	use data_in
	use defaultNotUsed
	use gis_data
	use sues_data     
	use time
	
	IMPLICIT NONE
	
	integer	:: i,ii
	integer	:: Gridiv
	
	real	:: dqndt 		! rate of change of net radiation [W m-2 h-1] at t-2
	real	:: surfrac  	! surface fraction accounting for SnowFrac if appropriate
	real	:: xa1,xa2,xa3 	! temporary AnOHM coefs.
  
   
   
! ------AnOHM coefficients --------
!	write(*,*) 'n surf:', nsurf
!	Loop through surface types at the beginning of a day------------------------
	if ( it==0 .and. imin==5 ) then
! 		write(*,*) '----- AnOHM called -----'
! 		write(*,*) 'Grid@id:', Gridiv, id
		! ------Set to zero initially------
			a1AnOHM(Gridiv) = 0.1   ![-]
			a2AnOHM(Gridiv) = 0.2   ![h]
			a3AnOHM(Gridiv) = 10   	![W m-2]  
		!----------------------------------
		do 	is=1,nsurf
			surfrac=sfr(is)
	! 		write(*,*) 'surfrac of ', is, 'is: ',surfrac
	!	initialize the coefs.
		xa1 = 0.1
		xa2 = 0.2
		xa3 = 10
	!	call AnOHM to calculate the coefs.	
		call AnOHM_coef(is,id,Gridiv,xa1,xa2,xa3)
     
	!	calculate the areally-weighted OHM coefficients
			a1AnOHM(Gridiv) = a1AnOHM(Gridiv)+surfrac*xa1  
			a2AnOHM(Gridiv) = a2AnOHM(Gridiv)+surfrac*xa2
			a3AnOHM(Gridiv) = a3AnOHM(Gridiv)+surfrac*xa3      
      
		enddo  
	!	end of loop over surface types -----------------------------------------
! 		write(*,*) '----- OHM coeffs -----'
! 		write(*,*) a1AnOHM(Gridiv),a2AnOHM(Gridiv),a3AnOHM(Gridiv)
	end if
	

	

!   Calculate radiation part ------------------------------------------------------------
  	qs=NAN  		!qs  = Net storage heat flux  [W m-2]
  	if(qn1>-999) then	!qn1 = Net all-wave radiation [W m-2]     
		dqndt = 0.5*(qn1-q2_grids(Gridiv))*nsh_real  	  !gradient at t-1
! 		Calculate net storage heat flux
  	   	qs = qn1*a1AnOHM(Gridiv)+dqndt*a2AnOHM(Gridiv)+a3AnOHM(Gridiv)   		!Eq 4, Grimmond et al. 1991
		
		q1_grids(Gridiv) = q2_grids(Gridiv)	!q1 = net radiation at t-2 (at t-3 when q1 used in next timestep)
		q2_grids(Gridiv) = q3_grids(Gridiv)	!q2 = net radiation at t-1
		q3_grids(Gridiv) = qn1				!q3 = net radiation at t   (at t-1 when q3 used in next timestep)
  	else
     	call ErrorHint(21,'Bad value for qn1 found during OHM calculation',qn1,NotUsed,notUsedI)
  	endif    

end subroutine AnOHM_v2016
!========================================================================================


!========================================================================================
subroutine AnOHM_coef(sfc_typ,xid,xgrid,& 	! input
					xa1,xa2,xa3)          	! output
! author: Ting Sun
! 
! purpose:
! calculate the OHM coefs. (a1, a2, and a3) based on forcings and sfc. conditions
! 
! input:
! 1) sfc_typ: surface type. 
! these properties will be loaded:
! xemis: emissivity, 1
! xcp: heat capacity, J m-3
! xk: thermal conductivity, W m K-1
! xch: bulk turbulent transfer coefficient,
! Bo: Bowen ratio (i.e. QH/QE), 1
! 2) xid: day of year
! will be used to retrieve forcing diurnal cycles of ixd.
! 
! output:
! a1, a2, and a3 
! in the relationship:
! delta_QS = a1*(Q*)+a2*(dQ*/dt)+a3
! 
! ref:
! the AnOHM paper to be added.
! 
! history:
! 20160222: initial version
!========================================================================================
	use allocateArray
	use data_in
	use defaultNotUsed
	use gis_data
	use sues_data     
	use time
	
	implicit none

! 	input
	integer:: sfc_typ, xid, xgrid

! 	output
	real :: xa1, xa2, xa3

! 	constant
	real, parameter :: SIGMA = 5.67e-8   		! Stefan-Boltzman
	real, parameter :: PI    = atan(1.0)*4 		! Pi
	real, parameter :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
	real, parameter :: C2K   = 273.15    		! degC to K

! 	sfc. properties:
	real :: xalb,	&    ! 	albedo, 
		 	xemis, 	&    ! 	emissivity, 
			xcp,	&    ! 	heat capacity, 
			xk, 	&    ! 	thermal conductivity,
			xch, 	&    ! 	bulk transfer coef.
			xBo	         ! 	Bowen ratio

! 	forcings: 	
	real, dimension(24) :: Sd,& ! 	incoming solar radiation
						   Ta,& ! 	air temperature
						   WS,& ! 	wind speed
						   WF,& ! 	anthropogenic heat
						   AH	! 	water flux density


!	local variables:
	real :: beta 			! inverse Bowen ratio
	real :: f,fL,fT 		! energy redistribution factors
	real :: lambda 			! thermal diffusivity
	real :: delta,m,n		! water flux related variables
	real :: ms,ns			! m, n related
	real :: gamma			! phase lag scale
	real :: ASd,mSd			! solar radiation
	real :: ATa,mTa			! air temperature
	real :: tau				! phase lag between Sd and Ta (Ta-Sd)
	real :: ATs,mTs			! surface temperature amplitude
	real :: ceta,cphi		! phase related temporary variables
	real :: eta,phi,xlag	! phase related temporary variables
	real :: mWS,mWF,mAH 	! mean values of WS, WF and AH
	real :: xx1,xx2,xx3		! temporary use
	real :: solTs			! surface temperature

! 	give fixed values for test
! 	properties
! 	xalb  = .2
! 	xemis = .9
! 	xcp   = 1e6
! 	xk    = 1.2
! 	xch   = 4
! 	xBo    = 2
! 	forcings
! 	ASd = 400
! 	mSd = 100
! 	ATa = 23
! 	mTa = 23+C2K
! 	tau = PI/6
! 	WS  = 1
! 	AH  = 0
! 	Wf  = 0

! 	load met. forcing data: 
	call AnOHM_FcLoad(sfc_typ,xid,xgrid,&	! input
					Sd,Ta,WS,WF,AH) 		! output
! 	write(*,*) 'here the forcings:'
! 	write(*,*) Sd,Ta,WS,WF,AH
	
! 	write(unit=*, fmt=*) 'DOY:', xid

!	load forcing characteristics:
	call AnOHM_FcCal(Sd,Ta,WS,WF,AH,&				! input
				  ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)	! output
! 	write(*,*) ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH
! 	load sfc. properties:
	call AnOHM_SfcLoad(sfc_typ,xid,xgrid,&			! input
						xalb,xemis,xcp,xk,xch,xBo)  ! output
! 	write(*,*) 'here the properties:'
! 	write(*,*) xalb,xemis,xcp,xk,xch,xBo

! 	initial Bowen ratio
	BoAnOHMStart(xgrid) = xBo
	
! 	calculate sfc properties related parameters:
	xch     = xch*mWS
	beta   = 1/xBo
	f      = ((1+beta)*xch+4*SIGMA*xemis*mTa**3)
	fL     = 4*SIGMA*xemis*mTa**3
	fT     = (1+beta)*xch
	lambda = xk/xcp
	delta  = sqrt(.5*(mWF**2+sqrt(mWF**4+16*lambda**2*OMEGA**2)))
	m      = (2*lambda)/(delta+mWF)
	n      = delta/OMEGA

	
! 	calculate surface temperature related parameters:
! 	mTs   = (mSd*(1-xalb)/f)+mTa
	ms    = 1+xk/(f*m)
	ns    = xk/(f*n)
	xx1   = f*sin(tau)*ATa
! 	write(*,*) 'ns,ms,f,tau,ATa:', ns,ms,f,tau,ATa
	xx2   = (1-xalb)*ASd+f*cos(tau)*ATa
! 	write(*,*) 'xalb,ASd,f:', xalb,ASd,f
	gamma = atan(ns/ms)+atan(xx1/xx2)
	
	ATs   = -(sin(tau)*ATa)/(ns*cos(gamma)-ms*sin(gamma))
	
! 	calculate net radiation related parameters:
	xx1  = (ns*cos(gamma)+sin(gamma)-ms*sin(gamma))*sin(tau)*ATa*fL
	xx2  = (xalb-1)*(ns*cos(gamma)-ms*sin(gamma))*ASd
	xx3  = (-ms*cos(tau)*sin(tau)+cos(gamma)*(ns*cos(tau)+sin(tau)))
	xx2  = xx2-xx3
	phi  = atan(xx1/xx2)
	xx3  = (ns*cos(gamma)-ms*sin(gamma))
	xx1  = (1+sin(gamma)/xx3)*sin(tau)*ATa*fL
	xx2  = (xalb-1)*ASd-(cos(tau)+cos(gamma)*sin(tau)/xx3)*ATa*fL
	cphi = sqrt(xx1**2+xx2**2)

! 	calculate heat storage related parameters:
	xx1  = m*cos(gamma)-n*sin(gamma)
	xx2  = m*sin(gamma)+n*cos(gamma)
	eta  = atan(xx1/xx2)
! 	write(*,*) 'm,n,gamma:', m,n,gamma	
	xx1  = xk**2*(m**2+n**2)*ATs**2
	xx2  = m**2*n**2
	ceta = sqrt(xx1/xx2)	
	
! 	calculate the OHM coeffs.:
	xlag = eta-phi
! 	write(*,*) 'eta,phi:', eta,phi	
	xa1  = (ceta*cos(xlag))/cphi
! 	write(*,*) 'ceta,xlag,cphi:', ceta,xlag,cphi
	
	xa2  = (ceta*sin(xlag))/(OMEGA*cphi)
	xa2  = xa2/3600 ! convert the unit from s-1 to h-1
	xa3  = (mSd*(xalb-1)*ceta*cos(xlag)*fT)/(f*cphi)-mAH
	
! 	write(*,*) 'sfc_typ:', sfc_typ
! 	write(*,*) 'a1,a2,a3:', xa1,xa2,xa3
	
end subroutine AnOHM_coef
!========================================================================================

!========================================================================================
subroutine AnOHM_FcCal(Sd,Ta,WS,WF,AH,&  				! input
					ASd,mSd,ATa,mTa,tau,mWS,mWF,mAH)	! output
! author: Ting Sun
! 
! purpose:
! calculate the key parameters of a sinusoidal curve for AnOHM forcings
! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
! 
! input:
! hourly values between 00:00 and 23:00 (local time, inclusive):
! 1) Sd: incoming shortwave radiation,  W m-2
! 2) Ta: air temperature, K
! 3) WS: wind speed, m s-1
! 4) WF: water flux density, ???
! 5) AH: anthropogenic heat, W m-2
! 
! output:
! 1) ASd, mSd: amplitude and mean value of Sd
! 2) ATa, mTa: amplitude and mean value of Ta
! 3) tau: phase lag between Sd and Ta
! 4) mWS: mean value of WS
! 4) mWF: mean value of WF
! 4) mAH: mean value of AH
! 
! ref:
! 
! history:
! 20160224: initial version
!========================================================================================
	implicit none

! 	input
	real :: Sd(24),&	! 
			Ta(24),&	! 
			WS(24),&	! 
			Wf(24),&	! 
			AH(24)		! 
	
! 	output	
	real :: ASd,mSd,&	! Sd scales
			ATa,mTa,&	! Ta scales
			tau,&		! phase lag between Sd and Ta
			mWS,&		! mean WS
			mWF,&		! mean WF
			mAH		 	! mean AH

! 	constant
	real, parameter :: SIGMA = 5.67e-8   		! Stefan-Boltzman
	real, parameter :: PI    = atan(1.0)*4 		! Pi
	real, parameter :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
	real, parameter :: C2K   = 273.15    		! degC to K

!	local variables:
	real :: tSd,tTa			! peaking timestamps
	real :: aCosb,aSinb,b,c ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
	real :: xx				! temporary use

! 	calculate sinusoidal scales of Sd:
	call AnOHM_ShapeFit(Sd(10:16),ASd,mSd,tSd)
! 	write(*,*) 'ASd:', ASd
! 	write(*,*) 'mSd:', mSd
! 	write(*,*) 'tSd:', tSd

! 	calculate sinusoidal scales of Ta:
	call AnOHM_ShapeFit(Ta(10:16),ATa,mTa,tTa)
! 	write(*,*) 'ATa:', ATa
! 	write(*,*) 'mTa:', mTa
! 	write(*,*) 'tTa:', tTa

! 	calculate the phase lag between Sd and Ta:
	tau = (tTa-tSd)/24*2*PI
! 	write(*,*) 'tau:', tau
	
! 	calculate the mean values:
	mWS = sum(WS(10:16))/7	! mean value of WS
	mWF = sum(WF(10:16))/7	! mean value of WF
	mAH = sum(AH(10:16))/7	! mean value of AH
! 	write(*,*) 'mWS:', mWS
! 	write(*,*) 'mWF:', mWF
! 	write(*,*) 'mAH:', mAH


end subroutine AnOHM_FcCal
!========================================================================================

!========================================================================================
subroutine AnOHM_ShapeFit(obs,amp,mean,tpeak)
! author: Ting Sun
! 
! purpose:
! calculate the key parameters of a sinusoidal curve for AnOHM forcings
! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
! 
! input:
! obs: hourly values between 10:00 and 16:00 (local time, inclusive)
! 
! output:
! 1) amp  : amplitude
! 2) mean : mean value
! 3) tpeak: daily peaking hour (h)
! 
! ref:
! 
! history:
! 20160224: initial version
!========================================================================================
	implicit none

! 	input
	real :: obs(7)	! observation (daytime 10:00–16:00)
	
! 	output	
	real :: amp		! amplitude
	real :: mean	! average
	real :: tpeak	! peaking time (h)

! 	constant
	real, parameter :: SIGMA = 5.67e-8   		! Stefan-Boltzman
	real, parameter :: PI    = atan(1.0)*4 		! Pi
	real, parameter :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
	real, parameter :: C2K   = 273.15    		! degC to K

!	local variables:
	real 	:: coefs(3,7) 			! coefficients for least squre solution
	real 	:: aCosb,aSinb,a,b,c 	! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
	real 	:: xx,mbias					! temporary use
	integer	:: i					! temporary use

! 	c coeffs.:
! 	exact values:
! 	coefs(1,:) = (/2*(-615 - 15*Sqrt(2.) + 211*Sqrt(3.) + 80*Sqrt(6.)),	&
! 				  3*(-414 - 12*Sqrt(2.) + 140*Sqrt(3.) + 75*Sqrt(6.)),		&
! 				  -1437 + 156*Sqrt(2.) + 613*Sqrt(3.) + 92*Sqrt(6.), 		&
! 				  2*(-816 + 177*Sqrt(2.) + 337*Sqrt(3.) + 13*Sqrt(6.)),	&
! 				  -1437 + 156*Sqrt(2.) + 613*Sqrt(3.) + 92*Sqrt(6.),		&
! 				  3*(-414 - 12*Sqrt(2.) + 140*Sqrt(3.) + 75*Sqrt(6.)),		&
! 				  2*(-615 - 15*Sqrt(2.) + 211*Sqrt(3.) + 80*Sqrt(6.))/)
! 	coefs(1,:) = coefs(1,:)/(-9450 + 534*Sqrt(2.) + 3584*Sqrt(3.) + 980*Sqrt(6.))
! 	apprx. values:
	coefs(1,:) = (/1.72649, 0.165226, -0.816223, -1.15098, -0.816223, 0.165226, 1.72649/)
	c          = dot_product(coefs(1,:),obs)

! 	aCos coeffs.:
! 	exact values:
! 	coefs(2,:) = (/249*(-6 + Sqrt(2.)) + Sqrt(3.)*(124 + 347*Sqrt(2.)),	&
! 				   -3*(-263 - 372*Sqrt(2.) + Sqrt(3.) + 325*Sqrt(6.)),		&
! 				   81 - 48*Sqrt(2.) - 247*Sqrt(3.) + 174*Sqrt(6.),			&
! 				   3*(27 - 545*Sqrt(2.) - 45*Sqrt(3.) + 340*Sqrt(6.)),		&
! 				   1740 - 303*Sqrt(2.) - 632*Sqrt(3.) - 73*Sqrt(6.),		&
! 				   -207 + 1368*Sqrt(2.) - 505*Sqrt(3.) - 338*Sqrt(6.),		&
! 				   -990 - 747*Sqrt(2.) + 1398*Sqrt(3.) - 155*Sqrt(6.)/)
! 	coefs(2,:) = coefs(2,:)/(-9450 + 534*Sqrt(2.) + 3584*Sqrt(3.) + 980*Sqrt(6.))
! 	apprx. values:
	coefs(2,:) = (/0.890047, 0.302243, -0.132877, -0.385659, -0.438879, -0.288908, 0.054033/)
	aCosb      = dot_product(coefs(2,:),obs)
	
!	cSin coeffs.:
! 	exact values: 
! 	coefs(3,:) = (/-28 - 13*Sqrt(2.) - 5*Sqrt(3.)*(-6 + Sqrt(2.)),	&
! 				   -15 + Sqrt(2.) - 9*Sqrt(3.) + 12*Sqrt(6.),		&
! 				   63 - 7*Sqrt(2.) - 21*Sqrt(3.) - 5*Sqrt(6.),		&
! 				   -9 - 7*Sqrt(3.) + 11*Sqrt(6.),					&
! 				   -32 - 7*Sqrt(2.) + 28*Sqrt(3.) - Sqrt(6.),		&
! 				   -3 + 27*Sqrt(2.) - 5*Sqrt(3.) - 11*Sqrt(6.),	&
! 				   24 - Sqrt(2.) - 16*Sqrt(3.) - Sqrt(6.)/)
! 	coefs(3,:) = coefs(3,:)/(-224 + 36*Sqrt(2.) + 58*Sqrt(3.) + 28*Sqrt(6.))
! 	apprx. values:
	coefs(3,:) = (/1.64967, -0.0543156, -1.10791, -1.4393, -1.02591, 0.104083, 1.87368/)
	aSinb      = dot_product(coefs(3,:),obs)
	
! 	calculate b and a:
	b = atan(aSinb/aCosb)
	if ( b>0 ) b = b-Pi 	! handle over Pi/2 phase lag
	a = aSinb/sin(b)
	
! 	mbias=0.
! 	write(*, *) '      obs        ','     sim       ','     diff     '
! 	do i = 1, 7, 1
! ! 	print out the sim and obs pairs
! 		xx    = a*sin(Pi/12*(9+i)+b)+c
! 		mbias = mbias+abs(xx-obs(i))
! 		write(*, *) obs(i),xx,xx-obs(i)
! 	end do
! 	mbias =mbias/7
! 	write(*, *) 'mean bias: ', mbias
	
	
! 	get results:
	amp=a					! amplitude
	mean=c 					! mean value
	tpeak=(PI/2-b)/PI*12 	! convert to timestamp (h)

end subroutine AnOHM_ShapeFit
!========================================================================================

!========================================================================================
subroutine AnOHM_FcLoad(sfc,xid,xgrid,&	! input
						Sd,Ta,WS,WF,AH) ! output
! author: Ting Sun
! 
! purpose:
! calculate the key parameters of a sinusoidal curve for AnOHM forcings
! i.e., a, b, c in a*Sin(Pi/12*t+b)+c
! 
! input:
! 1) sfc: surface type
! 2) xid : day of year
! 
! output:
! hourly values between 00:00 and 23:00 (local time, inclusive):
! 1) Sd: incoming shortwave radiation,  W m-2
! 2) Ta: air temperature, K
! 3) WS: wind speed, m s-1
! 4) WF: water flux density, ???
! 5) AH: anthropogenic heat, W m-2
! 
! ref:
! 
! 
! todo:
! check the completeness of forcings of given day (i.e., xid)
! 
! history:
! 20160226: initial version
!========================================================================================
	use allocateArray
	use ColNamesInputFiles
	use data_in
	use defaultNotUsed   
	use initial
	use sues_data
	use time

	implicit none

! 	input
	integer :: sfc, xid, xgrid 
		
! 	output	
	real :: Sd(24),&	!
			Ta(24),&	!
			WS(24),&	!
			WF(24),&	!
			AH(24)		!
			
! 	constant
	real, parameter :: SIGMA = 5.67e-8   		! Stefan-Boltzman
	real, parameter :: PI    = atan(1.0)*4 		! Pi
	real, parameter :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
	real, parameter :: C2K   = 273.15    		! degC to K

!	local variables:
	real :: tSd,tTa			! peaking timestamps
	real :: aCosb,aSinb,b,c ! parameters for fitted shape: a*Sin(Pi/12*t+b)+c
	real :: xx				! temporary use

	integer :: i			! temporary use
	
	
	integer :: irRange(2)	! row range in MetData containing id values
	integer :: nStepHour	! time steps in an hour
	integer :: lenMetData

	integer, allocatable :: lSub(:) ! array to retrieve data of the specific day (id)

	lenMetData = size(MetForcingData(:,2,xgrid))
	
	allocate(lSub(lenMetData))
	lSub = (/(i,i=1,lenMetData)/)

	where (MetForcingData(:,2,xgrid) == xid) 
		lSub = 1
		elsewhere
		lSub = 0	
	end where
	irRange = (/maxloc(lSub),maxloc(lSub)+sum(lSub)-1/)
! 	write(*,*) 'irRange:', irRange
	
	
! 	load the sublist into forcings:
	Sd = MetForcingData(irRange(1):irRange(2):nsh,15,xgrid) 	! nsh is timesteps per hour, a global var.
	Ta = MetForcingData(irRange(1):irRange(2):nsh,12,xgrid)
	WS = MetForcingData(irRange(1):irRange(2):nsh,10,xgrid)
	WF = MetForcingData(irRange(1):irRange(2):nsh,12,xgrid)*0	! set as 0 for debug
	if ( anthropheatchoice == 0 ) then
		AH = MetForcingData(irRange(1):irRange(2):nsh,9,xgrid)*0	! read in from MetForcingData,	
	else
		AH = mAH_grids(xid-1,xgrid)
	end if
	
! 	write(*,*) 'Sd:', Sd
! 	write(*,*) 'Ta:', Ta
! 	write(*,*) 'WS:', WS
! 	write(*,*) 'WF:', WF
! 	write(*,*) 'AH:', AH

end subroutine AnOHM_FcLoad
!========================================================================================

!========================================================================================
subroutine AnOHM_SfcLoad(sfc,xid,xgrid,&			! input
						xalb,xemis,xcp,xk,xch,xBo)	! output

! author: Ting Sun
! 
! purpose:
! load surface properties.
! 
! input:
! 1) sfc : surface type
! 2) xid : day of year
! 3) grid: grid number
! 
! output:
! surface properties:
! 1) alb,	 ! albedo, 
! 2) emiss,  ! emissivity, 
! 3) cp,	 ! heat capacity, 
! 4) k, 	 ! thermal conductivity,
! 5) ch, 	 ! bulk transfer coef.
! 6) Bo	     ! Bowen ratio	
! 
! ref:
! 
! history:
! 20160228: initial version
!========================================================================================
	use allocateArray
	use ColNamesInputFiles
	use data_in
	use defaultNotUsed   
	use initial
	use sues_data
	use time

	implicit none

! 	input
	integer :: sfc, xid, xgrid 
		
! 	output:
! 	sfc. properties:
	real :: xalb,		&    ! 	albedo, 
		 	xemis, 		&    ! 	emissivity, 
			xcp,		&    ! 	heat capacity, 
			xk, 		&    ! 	thermal conductivity,
			xch, 		&    ! 	bulk transfer coef.
			xBo 	         ! 	Bowen ratio		
! 	constant
	real, parameter :: SIGMA = 5.67e-8   		! Stefan-Boltzman
	real, parameter :: PI    = atan(1.0)*4 		! Pi
	real, parameter :: OMEGA = 2*Pi/(24*60*60)  ! augular velocity of Earth
	real, parameter :: C2K   = 273.15    		! degC to K

	
! 	load properties from global variables
	xalb  = alb(sfc)
	xemis = emis(sfc)
	xcp   = cpAnOHM(sfc)
	xk    = kkAnOHM(sfc)
	xch   = chAnOHM(sfc)
	
! 	write(*,*) 'xalb:', xalb
! 	write(*,*) 'xemis:', xemis
! 	write(*,*) 'sfc:', sfc
! 	write(*,*) 'xcp:', xcp
! 	write(*,*) 'xk:', xk
! 	write(*,*) 'xch:', xch
	
	
! 	load Bowen ratio of yesterday from DailyState calculation
	if ( BoAnOHMEnd(xgrid) == NAN ) then
		xBo = Bo_grids(xid-1,xgrid)
	else
		xBo = BoAnOHMEnd(xgrid)		
	end if
	
! 	write(*,*) 'xBo:', xBo

end subroutine AnOHM_SfcLoad
!========================================================================================
