!  FUNCTION TO RETURN 0 IF IX=0, 1 IF 0<IX<MAXPOS+1,-1 OTHERWISE. 
!  MAXPOS is given as the maximum positive Integer. 
subroutine issign(IX,MAXPOS,ISIGNM)
      real(kind(1d0)) IX,MAXPOS,ISIGNM 
      ISIGNM=1.0 
      IF(IX.LT.0.OR.IX.GT.MAXPOS)ISIGNM=-1 
      IF(IX.EQ.0) ISIGNM=0 
      RETURN 
end subroutine issign


!=====================================================
! fuction to convert interger to string
character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    end function str
    

subroutine SaveGrids

use matsize
use solweig_module
use data_in
use time

implicit none
character(len=100)       :: GridPath,GridName,GridPath2
character(len=4)       ::doy,hour
!real(kind(1d0)), allocatable, dimension(:,:):: savegrid

    allocate(savegrid(sizey,sizex))

    if (Tmrt_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Tmrt_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Tmrt
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Lup2d_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Tmrt_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Tmrt
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Ldown2d_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Ldown2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Ldown2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Kup2d_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Kup2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Kup2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (Kdown2d_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='Kdown2d_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=Kdown2d
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if
    if (GVF_out==1) then
        Gridpath2='Grids/'    
        GridPath=trim(FileOutputPath)//trim(GridPath2)
        write(doy,'(i3)') id 
        write(hour,'(i2)') it 
        hour=adjustl(hour)
        doy=adjustl(doy)
        GridName='GVF_'//trim(doy)//'_'//trim(hour)//'.txt'
        savegrid=gvf
        call SaveEsriAsciiGrid(GridPath,GridName,xllcorner,yllcorner,cellsize,NoData)
    end if    
    
    deallocate(savegrid)

end subroutine SaveGrids