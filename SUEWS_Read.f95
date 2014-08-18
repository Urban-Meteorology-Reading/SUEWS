
!Read first how many lines in the input file there is. Then read met forcing data into matrix datamet
 subroutine SUEWS_Read(lunit)

 use data_in
 use sues_data
 use time
 use defaultnotUsed


 IMPLICIT NONE

 integer::lunit, i

 !First read the number of lines in the forcing data
 open(lunit,file=trim(fileMet),status='old',err=314,position='rewind')
 call skipHeader(lunit,SkipHeaderMet)

 nlines = 0 !Initialize nlines
 DO
    READ (lunit,*, END=10)
    nlines = nlines + 1
 END DO
 10 CLOSE (lunit)

 !Allocate hourly input and ouutput matrixes
 allocate(dataMet1(nlines,2))
 allocate(dataMet2(nlines,21))
 allocate(dataOut1(nlines,62))
 allocate(dataOut2(nlines,30))
 allocate(dataOut3(nlines,106))
 allocate(dataOut5min(nlines*nsh,69))
 allocate(dataOutBL(nlines*4,20))
 allocate(dataOutSOL(nlines,28))

!Open the file again. Code later better
 open(1,file=trim(fileMet),status='old',err=314,position='rewind')
 call skipHeader(1,SkipHeaderMet)

 DO i=1,nlines

    call MetRead(i)
    dataMet1(i,1:2) = (/id,it/)
    dataMet2(i,1:21) = (/dectime,qn1_obs,qh_obs,qe_obs,qs,qf,avu1,avrh,Temp_C,Press_hPa,&
                       Precip_hr,avkdn,snow_obs,ldown_obs,fcld_obs,wuh,xsmd,lai_hr,kdiff,kdir,wdir/)
 ENDDO

 CLOSE(1)
 return

314 call errorHint(11,trim(filemet),notUsed,notUsed,ios_out)


 end subroutine


!#####################################
 subroutine ConvertMetData(i)

    use time
    use data_in

    IMPLICIT NONE

    integer::i


    id = dataMet1(i,1)   !Integer variables
    it = dataMet1(i,2)

    dectime = dataMet2(i,1)
    qn1_obs = dataMet2(i,2)
    qh_obs = dataMet2(i,3)
    qe_obs = dataMet2(i,4)
    qs = dataMet2(i,5)
    qf = dataMet2(i,6)
    avu1 = dataMet2(i,7)
    avrh = dataMet2(i,8)
    Temp_C = dataMet2(i,9)
    Press_hPa = dataMet2(i,10)
    Precip_hr = dataMet2(i,11)
    avkdn = dataMet2(i,12)
    snow_obs = dataMet2(i,13)
    ldown_obs = dataMet2(i,14)
    fcld_obs = dataMet2(i,15)
    wuh = dataMet2(i,16)
    xsmd = dataMet2(i,17)
    lai_hr = dataMet2(i,18)
    kdiff = dataMet2(i,19)
    kdir = dataMet2(i,20)
    wdir = dataMet2(i,21)

 end subroutine
