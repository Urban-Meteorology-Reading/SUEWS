!-----------------------------------------------------------------------
! from CBL modelling Cleugh and Grimmond (2000) BLM
!-----------------------------------------------------------------------
	SUBROUTINE RKUTTA(neqn,XA,XB,Y,NSTEPS)
!       XA=s0
!       XB=s1
!       Y(1)=blh_m
!       Y(2)=tm_K
!       Y(3)=qm_kgkg
!       Y(4)=cm    
!       JOHN KNIGHT, 1985 (AMENDED BY MRR, 23-SEP-85)
!       EXPLICIT FOURTH-ORDER RUNGE-KUTTA METHOD FOR FIRST-ORDER ODE SYSTEM
!       OF NE EQUATIONS, WITH INITIAL VALUES SUPPLIED.
!       MEANING OF PARAMETERS:
!       NE     = NUMBER OF EQUATIONS (MAX 21)
!       XA     = VALUE OF INDEPENDENT VARIABLE ON ENTRY
!       XB     = VALUE OF INDEPENDENT VARIABLE AT END OF INTERVAL
!       Y(NE)  = ON ENTRY: INITIAL VALUES OF DEPENDENT VARIABLES AT XA
!       ON EXIT:  CALCULATED VALUES OF DEPENDENT VARIABLES AT XB
!       NSTEPS = NUMBER OF INTEGRATION STEPS OVER INTERVAL (XA,XB)
!       DIFF  = NAME OF USER-SUPPLIED SUBROUTINE TO CALCULATE DERIVATIVES
!       DYDX (DIFF MUST BE DECLARED EXTERNAL IN CALLING PROGRAM).
!       PARAMETERS IN SUBROUTINE DIFF(NE,X,Y,DYDX):
!       NEqn = NUMBER OF EQUATIONS
!       X = INDEPENDENT VARIABLE
!       Y = ARRAY (LENGTH NE) OF VALUES OF DEPENDENT VARIABLES
!       DYDX = ON EXIT, ARRAY (LENGTH NE) OF VALUES OF DERIVATIVES
!	IMPLICIT real*8 (A-H,O-Z)
	implicit none
    integer::ns,nsteps, nj,n,neqn
	real(kind(1D0)), dimension (neqn):: y
	real(kind(1D0)), dimension (21):: dydx,arg
	real(kind(1D0)), dimension (21,5):: rk
	real(kind(1D0)), dimension (4):: coef
    real (kind(1D0)):: XA,XB,step,X,xx
   
	coef(1)=1.0
	coef(2)=0.5
	coef(3)=0.5
	coef(4)=1.0
!	print*,"rk1: ",xa,xb,y
	STEP = (XB-XA)/NSTEPS 

	DO NS = 1,NSTEPS
	   DO  NJ = 1,nEqn
	      RK(NJ,1) = 0
	   enddo
	   X = XA+(NS-1)*STEP
	   DO N = 1,4
	      IF (N.EQ.1)then
		  XX = X 
	      elseIF (N.GT.1)then
		  XX = X + COEF(N)*STEP
	      endif
	      DO NJ = 1,nEqn
		     ARG(NJ) = Y(NJ) + COEF(N)*RK(NJ,N)
	      enddo

	      CALL DIFF(xx,ARG,DYDX)

	      DO NJ = 1,nEqn
		     RK(NJ,N+1) = STEP*DYDX(NJ)
	      enddo
	   enddo
         
	   DO  NJ = 1,nEqn
	      DO  N = 1,4
		  		Y(NJ) = Y(NJ) + RK(NJ,N+1)/(6*COEF(N))
	      enddo
	   enddo
	enddo
	
	RETURN
	END subroutine rkutta
!---------------------------------------------------------------------
!---------------------------------------------------------------------
	
	subroutine diff(s,y1,dyds)
    ! in y1,neqn
    ! out dyds
    
!       calculates derivatives for cbl slab model
!       y(1) = h = cbl depth(m)
!       y(2) = t = potential temp(K)
!       y(3) = q = specific humidity(kg/kg)
!       y(4) = c = CO2 concentration
    use data_in
    use sues_data
!    use allocateArray
     use time
     USE CBL_MODULE
    use defaultnotUsed
    use mod_grav	 
	implicit none
    real (kind(1D0)), dimension(neqn)::dyds,y1
	real(kind(1d0)) :: zero=0.0,s
    real(kind(1d0)) :: h1,t_K,q_kgkg,c,cp,ws
    real(kind(1D0)):: delt_K,delq_kgkg,delc
    real(kind(1D0)):: gamtv_Km,deltv_K,ftv_Kms
    real(kind(1D0)):: ftva_Kms,delb,qs2,qs3
    real(kind(1D0)):: dhds,dtds,dqds,dcds
    real(kind(1D0)):: conk,conn,cona,conc,cont

!	print*,"diff: timestamp:",s
!    pause
	h1     = y1(1)!m
	t_K    = y1(2)!K
	q_kgkg = y1(3)!kg/kg
	c      = y1(4)
    
!       find t, q, c above inversion, and jumps across inversion
!       tp = tp + gamt*h
!       qp = qp0 + gamq*h 

	cp        = 0 ! cp0 + gamc* h1   ! todo 
    
	delt_K    = tpp_K    - t_K 
	delq_kgkg = qpp_kgkg - q_kgkg
	delc      = cp - c 

!       find potential virtual temperature flux, gradient and jump
	ftv_Kms  = fhbl_Kms + 0.61 * tm_K * febl_kgkgms
	gamtv_Km = gamt_Km  + 0.61 * tm_K * gamq_kgkgm!/1000
	deltv_K  = delt_K   + 0.61 * tm_K * delq_kgkg 

!       find velocity scale ws
	ftva_Kms = max(ftv_Kms,zero) ! virtual heat flux
	ws = (h1*ftva_Kms*grav/tm_K)**0.3333333333

!       find dhds using one of 4 alternative schemes chosen by ient:
	if (EntrainmentType.eq.2) then
!       EntrainmentType=1: encroachment (as in McN and S 1986 eq 16))
	   dhds = ftva_Kms/(h1*gamtv_Km)
       
	else if (EntrainmentType.eq.1) then
!       EntrainmentType=2: Driedonks 1981 (as in McN and S 1986 eq 13)
	   if (deltv_K.le.0.01) then 
              dhds = ftva_Kms/(h1*gamtv_Km)        
              
	      call errorHint(30,"subroutine diff [CBL: Deltv_K<0.01 EntrainmentType=1], deltv_K,delt_K,",deltv_K,delt_K,notUsedI)
          call errorHint(30,"subroutine diff [CBL: Deltv_K<0.01 EntrainmentType=1], tm_K,TPP_K,y1",tm_K,TPP_K, notUsedI)
          call errorHint(31,"subroutine diff [CBL: Deltv_K<0.01 EntrainmentType=1], y1",y1,notUsed,notUsedI)
	   else
              delb = grav*deltv_K/tm_K
              conc = 0.2
              cona = 5.0
              dhds = (conc*ws**3 + cona*cbld(7)**3)/(h1*delb)
	   end if                          
   
	else if (EntrainmentType.eq.4) then
!       EntrainmentType=3: Tennekes 1973 (as in R 1991 eqs 3,4)
       alpha3=0.7
	   if (deltv_K.le.0.01) then 
	      dhds = ftva_Kms/(h1*gamtv_Km)
	      call ErrorHint(31, 'subroutine difflfnout: [CBL: deltv_K<0.01 EntrainmentType=4],deltv_K',&
          deltv_K,notUsed,notUsedI)
	   else    
	      dhds = alpha3*ftva_Kms/deltv_K
	   end if

!       write (4,*) tpp, gamq, dhds, deltv
	   
	else if (EntrainmentType.eq.3) then
!       EntrainmentType=4: Rayner and Watson 1991 eq 21
	   conn = 1.33
	   conk = 0.18
	   cont = 0.80
	   qs3 = ws**3 + (conn*cbld(7))**3
	   qs2 = qs3**(0.6666666667) 
	   
	   if (deltv_K.le.0.01) then 
              dhds = ftva_Kms/(h1*gamtv_Km)
         call ErrorHint(31, 'subroutine difflfnout: [CBL: deltv_K<0.01 EntrainmentType=3],deltv_K',&
         deltv_K,notUsed,notUsedI)

	    
	   else
              delb = grav*deltv_K/tm_K
              dhds = (conk*qs3) / (cont*qs2 + h1*delb)
	   end if

	else
	   call ErrorHint(24, 'BLUEWS_DIff- CBL- illegal alpha',notUsed,notUsed,notUsedI)
	end if
! find dtds, dqds, dc/ds:
!	wsb is the subsidence velocity. Try using: -0.01, -0.05, -0.1.   

	dtds = fhbl_Kms/h1    + delt_K    *(dhds-wsb)/h1
	dqds = febl_kgkgms/h1 + delq_kgkg *(dhds-wsb)/h1
	dcds = fcbl/h1        + delc      *(dhds-wsb)/h1
    
	dyds(1) = dhds
	dyds(2) = dtds
	dyds(3) = dqds
	dyds(4) = dcds
    

 	return
	end subroutine diff

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine sonde
! read sonde or vertical profile data - when available          
!use allocateArray
use data_in
use cbl_module
implicit none
integer::i,fn=101,izm=500,notUsedI=-9999
character (len=200)::FileN
real (kind(1d0)):: dxx
real (kind(1d0)),parameter::notUsed=-9999.99

	FileN=trim(FileInputPath)//trim(FileSonde(which_day))
    open(fn,file=FileN,status="old",err=24)
    ! todo gneric skip header
	read(fn,*)
    read(fn,*)
    read(fn,*)
       
	do i=1,1000
	   read(fn,*,end=900,err=25)gtheta(i,1),dxx,gtheta(i,2),ghum(i,1),dxx,ghum(i,2)      
       ghum(i,2) = ghum(i,2) 
	enddo
900	    zmax=i-1
		if(zmax.gt.izm)then
	         call ErrorHint(23,FileN,float(zmax),notUsed,izm)	
	    endif
        close(fn)
        return
24		call ErrorHint(24,FileN,notUsed,notUsed, notUsedI)
25		call ErrorHint(25,FileN,notUsed,notUsed,i)
return
end subroutine sonde
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine gamma_sonde
use cbl_module
!use allocateArray

implicit none
real(kind(1D0))::gamtt,gamqq
integer::j
! gtheta(i,1),dxx,gtheta(i,2),ghum(i,1),dxx,ghum(i,2) 
!search for correct gamma theta, depends on h(i-1), 
!               ie current value for mixed layer depth
 if (sondeflag.eq.1) then
      do j=2,zmax
          if (blh_m.ge.gtheta(j-1,1)) then
              gamtt = gtheta(j-1,2)
          endif 
          gamt_Km = gamtt    
      enddo
            
    do j=2,zmax
         if (blh_m.ge.ghum(j-1,1)) then
             gamqq = ghum(j-1,2)
         endif
         gamq_kgkgm = gamqq/1000.
    enddo
endif 
return

end subroutine gamma_sonde
 