
!==============================================================================
!Simple Anthropogenic Heat Parameterization routine
!Last modified LJ (9/2010). 
!		Addition of SAHP_2:Takes weekend and weekday profiles into account
!       Take daylightsaving into account

!========================================= 
SUBROUTINE SAHP_Coefs(temp_C0,id_prev)
use allocateArray
use data_in
implicit none


INTEGER        :: gamma1,gamma2,id_prev
real(kind(1d0)):: temp_C0
	
namelist/AnthropogenicHeat/NumCapita,&
	BaseTHDD,&
    QF_A,&
    QF_B,&
	QF_C,&
    AH_min,&
    AH_slope,&
    T_Critic

	QF_A=0
    QF_B=0
	QF_C=0
    AH_min=0
    AH_slope=0
    T_Critic=0
    
    lfnSAHP=41 !SAHP file

	open(lfnSAHP,file=trim(FileSAHP),status='old')
    
	read(lfnSAHP,nml=AnthropogenicHeat,err=200)
	write(12,nml=AnthropogenicHeat)
 	close(lfnSAHP)
 
if(AnthropHeatChoice==1) then
  if(AH_min==0.and.Ah_slope==0.and.T_Critic==0)then
     	call ProblemsText(trim(FileSAHP))    	
        write(500,*)AnthropHeatChoice, "=AnthropHeatChoice"
        write(500,*)"AH_min=0.and.Ah_slope=0.and.T_Critic=0"        
      	write(500,nml=AnthropogenicHeat)       	
        Call PauseStop
   endif
elseif(AnthropHeatChoice==2) then
  if(sum(QF_A)==0.and.sum(QF_B)==0.and.sum(QF_C)==0)then
    	call ProblemsText(trim(FileSAHP))  
        write(500,*)AnthropHeatChoice, "=AnthropHeatChoice"
        write(500,*)"QF_A=0.and.QF_B=0.and.QF_C=0"
      	write(500,nml=AnthropogenicHeat)
        Call PauseStop
   endif
endif

!Calculations related to heating and cooling degree days
if ((Temp_C0-BaseTHDD)>=0) then !Cooling
   gamma2=1
else
    gamma2=0	
endif

if ((BaseTHDD-Temp_C0)>=0) then !Heating
   gamma1=1
else
   gamma1=0	
endif

HDD(id_prev,1)=gamma1*(BaseTHDD-Temp_C0) ! Heating
HDD(id_prev,2)= gamma2*(Temp_C0-BaseTHDD) ! Cooling

 return
 200  	Call ProblemsText(trim(FileSAHP))
 		write(*,nml=AnthropogenicHeat)
      	write(500,nml=AnthropogenicHeat)
       	call PauseStop
      

END SUBROUTINE SAHP_Coefs

!============================================================================== 
SUBROUTINE SAHP(qf_o,ih,id)
use data_in
use allocateArray
implicit none
INTEGER::id,iu=2,ih ! iu 1=weekend 2=weekday
real (kind(1d0)):: qf_o

if(DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) then  ! weekend
  iu=1  
endif

!Linear by Thomas  
!Loridan et al. (2011) JAMC
  
if(Temp_C.lt.T_CRITIC) then
   QF_o=AHPROF(ih,iu)*(AH_MIN+AH_SLOPE*(T_CRITIC-Temp_C))
else
   QF_o=AHPROF(ih,iu)*AH_MIN
endif 

END SUBROUTINE SAHP

!============================================================================== 
SUBROUTINE SAHP_2(qf_o,ih,id)
use data_in
use allocateArray
implicit none
INTEGER::id,iu=2,ih ! iu 1=weekend 2=weekday
real (kind(1d0)):: QF_o

!New one by Leena
! Jarvi et al. (2011) JH

if (DayofWeek(id,1)==1.or.DayofWeek(id,1)==7) then  ! weekend
  iu=1  
endif

QF_o=(AHPROF(ih,iu)*(Qf_a(iu)+Qf_b(iu)*HDD(id-1,2)+Qf_c(iu)*HDD(id-1,1)))*numCapita !Change numCapita to densities  
!print*,dayofweek(id,1),id,iu

  
!print*,AHPROF_COEF(iu),Qf_a(iu),Qf_b(iu),HDD(id-1,2),Qf_c(iu),HDD(id-1,1),numCapita
!pause

END SUBROUTINE SAHP_2                           

!============================================================================== 

