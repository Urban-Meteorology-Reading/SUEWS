! This subroutine loads a ESRIACSII grid as a 2D array
    
    subroutine LoadEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)	
    use matsize
    
	implicit none
    real(kind(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    integer                           :: col,row
    character(len=100)                :: GridPath,GridName,GridFile,n
    !real(kind(1d0)),intent(out)       :: tempgrid
    !real(kind(1d0)),dimension(sizey,sizex),intent(out)       :: tempgrid!,allocatable
    
    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=trim(GridPath)//trim(GridName)
  	open(99,File=GridFile,status='old') 
    
	! Read Header
    read(99,*) n,sizex
    read(99,*) n,sizey
    read(99,*) n,xllcornerlo
    read(99,*) n,yllcornerlo
    read(99,*) n,cellsizelo
    read(99,*) n,NoDatalo

    allocate(tempgrid(sizex,sizey))     
    
	! Read Matrix
    do row=1,sizey
      	read(99,*) (tempgrid(row,col),col=1,sizex)
    end do
    close(99)

    return
    end subroutine LoadEsriAsciiGrid
    
    
! This subroutine saves a 2D array as an ESRIACSII grid
subroutine SaveEsriAsciiGrid(GridPath,GridName,xllcornerlo,yllcornerlo,cellsizelo,NoDatalo)	
    use matsize
    
    implicit none
    real(kind(1d0))                   :: xllcornerlo,yllcornerlo,cellsizelo,NoDatalo
    integer                           :: col,row
    character(len=100)                :: GridPath,GridName,GridFile
    !integer                           :: sizey,sizex!,intent(in)
    !real(kind(1d0)), allocatable, dimension(:,:):: grid
   ! real(kind(1d0)),dimension(sizey,sizex)       :: grid!,allocatable
    
    ! Loading DSM
    !GridPath='D:\SOLWEIG2013b_Fortran\Inputdata\'
    !GridName='kr_dem.asc'
    GridFile=trim(GridPath)//trim(GridName)
    open(94,File=GridFile,status='unknown') 
    
    ! Read Header
    write(94,"(A5,1x,I0)") 'ncols',sizex
    write(94,"(A5,1x,I0)") 'nrows',sizey
    write(94,"(A9,1x,F0.2)") 'xllcorner',xllcornerlo
    write(94,"(A9,1x,F0.2)") 'yllcorner',yllcornerlo
    write(94,"(A8,1x,F0.2)") 'cellsize',cellsizelo
    write(94,"(A12,1x,F0.2)") 'NODATA_value',NoDatalo

	! write Matrix
    do row=1,sizey
        write(94,100) (savegrid(row,col),col=1,sizex)
    end do
    close(94)
100 format((f6.2,1x))
 
    return
    end subroutine SaveEsriAsciiGrid