! sg Feb 2012 -- moving all related suborutine together
! subroutine run_control - checks that value is reasonable
! subrotuine skipHeader  - jumps over header rows in input files
! Last modified:
! LJ 27 Jan 2016 - Removal of tabs, cleaning of the code
!------------------------------------------------------------

 !Information for the run
 MODULE run_info
   IMPLICIT NONE
   CHARACTER (len=90),DIMENSION(14)::text
   INTEGER::lim0=0,lim1=1,lim2=2,lim4=4,lim3=3,lim6=6,lim8=8,lim12=12,lfn_us
   LOGICAL ::file_qs
 END MODULE run_info

! run_control
! called from: LUMPS_initial

 SUBROUTINE run_control(eval,LowerLimit,Upperlimit)
   ! ver - determines if value to be read is an integer or real and returns the value
   ! if ver=-9 - then use integer
   USE run_info
   IMPLICIT NONE
   INTEGER::eval,i,lowerlimit,upperlimit
   CHARACTER (len=4)::check

   IF(file_qs)THEN
      101  READ(lfn_us,*)check
      write(*,*)check
      DO i=1,3
        IF(check(i:i)=="#")THEN
           !   write(*,*)check(i:i),i,"check"
           GOTO 101
        ELSE
           BACKSPACE(lfn_us)
           READ(lfn_us ,*)eval
           !write(*,*)eval, TEXT(1)
           exit
        ENDIF
      ENDDO
   ENDIF

   WRITE(12,120)eval,text(1)

   IF(eval<Lowerlimit.OR.eval>upperlimit)THEN
     WRITE(*,*)"Value out of range"
     WRITE(*,*)eval,text(1)
     STOP
   ENDIF

   WRITE(*,120)eval,text(1)
   120 FORMAT(i4,2x,a90)

   RETURN
 END SUBROUTINE run_control

! Skip Headers------------------------------
 SUBROUTINE SkipHeader(lfn,skip)
   use defaultnotUsed
   IMPLICIT NONE

   integer::skip,lfn,i
   do I=1,skip
     read(lfn,*,err=201,iostat=ios_out)
   end do

   RETURN

   201  reall=real(skip)
   call ErrorHint(20,'In SkipHeader subroutine.',reall,notUsed,ios_out)
 END SUBROUTINE SkipHeader
