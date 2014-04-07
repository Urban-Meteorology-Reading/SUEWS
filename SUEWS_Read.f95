!In this subroutine the output files will be printed out.

 subroutine SUEWS_Read(lunit)

 use data_in
 use time

 IMPLICIT NONE

 integer::lunit, i

 nlines = 0
 DO
    READ (lunit,*, END=10)
    nlines = nlines + 1
 END DO

 10 CLOSE (lunit)

 !nlines = 744

 allocate(dataMet(nlines,20))


!Open again
 open(1,file=trim(fileMet),status='old',position='rewind')
 call skipHeader(1,SkipHeaderMet)

 DO i=1,nlines

    call MetRead(i)
    !call GisRead(i)

    dataMet(i,1) = real(id,kind(1d0))
    dataMet(i,2) = real(it,kind(1d0))
    dataMet(i,3) = dectime
    dataMet(i,4) = qn1_obs
    dataMet(i,5) = qh_obs
    dataMet(i,6) = qe_obs
    dataMet(i,7) = qs
    dataMet(i,8) = qf
    dataMet(i,9) = avu1
    dataMet(i,10) = avrh
    dataMet(i,11) = Temp_C
    dataMet(i,12) = Pres_kPa
    dataMet(i,13) = Precip_hr
    dataMet(i,14) = avkdn
    dataMet(i,15) = snow
    dataMet(i,16) = ldown_obs
    dataMet(i,17) = fcld_obs
    dataMet(i,18) = wuh
    dataMet(i,19) = xsmd
    dataMet(i,20) = lai_hr


 ENDDO

 CLOSE(1)

 end subroutine


 subroutine ConvertMetData(i)

    use time
    use data_in

    IMPLICIT NONE

    integer::i

    id = int(dataMet(i,1))
    it = int(dataMet(i,2))
    dectime = dataMet(i,3)
    qn1_obs = dataMet(i,4)
    qh_obs = dataMet(i,5)
    qe_obs = dataMet(i,6)
    qs = dataMet(i,7)
    qf = dataMet(i,8)
    avu1 = dataMet(i,9)
    avrh = dataMet(i,10)
    Temp_C = dataMet(i,11)
    Pres_kPa = dataMet(i,12)
    Precip_hr = dataMet(i,13)
    avkdn = dataMet(i,14)
    snow = dataMet(i,15)
    ldown_obs = dataMet(i,16)
    fcld_obs = dataMet(i,17)
    wuh = dataMet(i,18)
    xsmd = dataMet(i,19)
    lai_hr = dataMet(i,20)

 end subroutine
