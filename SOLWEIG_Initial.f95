! setup for SOLWEIG
! FL may 2014
subroutine SOLWEIG_Initial
use matsize         ! All allocatable grids and related variables used in SOLWEIG
use InitialCond 
use allocateArray
use data_in  
use sues_data
use defaultNotUsed
use InitialCond
use solweig_module
use time

implicit none
    
character(len=100)  :: Path,GridFile,GridFolder
real(kind(1d0))                 :: vegmax
character(len=100),dimension(5) :: svfname
character(len=100),dimension(10):: svfvegname
logical                         :: exist
integer                         :: firstday

namelist/SOLWEIGinput/Posture,&    ! 1.Standing, 2.Sitting
    absL,&            ! Absorption coefficient of longwave radiation of a person  
    absK,&            ! Absorption coefficient of shortwave radiation of a person
    heightgravity,&   ! Center of gravity for a standing person
    usevegdem,&       ! With vegetation (1)
    DSMPath,&         ! Path to DSMs
    DSMname,&         ! Ground and building DSM
    CDSMname,&        ! Canopy DSM
    TDSMname,&        ! Trunk zone DSM
    TransMin,&        ! Transmissivity of K through decidious vegetation (leaf on)
    TransMax,&        ! Transmissivity of K through decidious vegetation (leaf off)
    SVFPath,&         ! Path to SVFs
    SVFsuffix,&       !
    buildingsname,&   ! Boolean matrix for locations of building pixels 
    row,&             ! X coordinate for point of interest
    col,&             ! Y coordinate for point of interest
    onlyglobal,&      ! if no diffuse and direct, then =1
    SOLWEIGpoi_out,&  ! write output variables at point of interest
    Tmrt_out,&        ! write Tmrt grid to file
    Lup2d_out,&       ! write Lup grid to file
    Ldown2d_out,&     ! write Ldown grid to file
    Kup2d_out,&       ! write Kup grid to file
    Kdown2d_out,&     ! write Kdown grid to file
    GVF_out,&         ! write GroundViewFactor grid to file
    SOLWEIG_ldown,&   ! 1= use SOLWEIG code to estimate Ldown, 0=use SEUWS
    OutInterval,&     ! Output interval in minutes
    RunForGrid        ! If only one grid should be run. All grids -999

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !Read in the SOLWEIGinput.nml file
    open(52,file=trim(FileInputPath)//'SOLWEIGinput.nml',err=274,status='old')
    read(52,nml=SOLWEIGinput)
    close(52)
    
    SolweigCount=1
    
    if (OutInterval == 60) then
        OutInterval=0
    endif
    
    if (Posture==1) then
        Fside=0.22
        Fup=0.06
    else
        Fside=0.1666666667
        Fup=0.166666667
    endif
   
    !timestepdec=real(t_interval)/(3600*24)
    timestepdec=real(OutInterval)/(real(t_interval)*24.)
 
	!!! Loading DSM !!!
    Path=trim(FileInputPath)//trim(DSMPath)
    call LoadEsriAsciiGrid(Path,DSMName,xllcorner,yllcorner,cellsize,NoData)
    allocate(a(sizey,sizex))
    a=tempgrid
    deallocate(tempgrid)
    
    scale=1/cellsize
    
    GridFolder=trim(FileOutputPath)//'Grids'
    
    ! Create grid folder !! This does not work in windows. needs to be done in python
    !inquire(file=GridFolder, exist=exist)
    !if (exist) then
    !else
    !    makedirectory = 'mkdir ' // trim(GridFolder)
    !    call system(makedirectory)
    !end if
    
    !!! Set up for vegetation scheme, or not !!!
    if (usevegdem==1) then 
        ! Calculating transmissivity of short wave radiation  through vegetation based on decid lai
        transperlai=(TransMax-TransMin)/(laimax(2)-laimin(2))
        firstday = MetForcingData(1,2,1) 
        trans=TransMin+(laimax(2)-lai(firstday-1,2))*transperlai
               	
	! Loading vegDSM (SDSM)
        Path=trim(FileInputPath)//trim(DSMPath)
    	call LoadEsriAsciiGrid(Path,CDSMname,xllcorner,yllcorner,cellsize,NoData)
        allocate(vegdem(sizey,sizex))
        vegdem=tempgrid
        deallocate(tempgrid)

        ! Loading trunkDSM (TDSM)
        Path=trim(FileInputPath)//trim(DSMPath)
    	call LoadEsriAsciiGrid(Path,TDSMname,xllcorner,yllcorner,cellsize,NoData)
        allocate(vegdem2(sizey,sizex))
        vegdem2=tempgrid
        deallocate(tempgrid)        
        
    	! amaxvalue (used in calculation of vegetation shadows)
    	vegmax=maxval(vegdem)
    	amaxvalue=maxval(a)-minval(a)
    	amaxvalue=max(amaxvalue,vegmax)
    
    	! Elevation vegdems if buildingDSM includes ground heights
    	vegdem=vegdem+a
    	where (vegdem==a)
        	vegdem=0.0
    	end where    
	vegdem2=vegdem2+a;
    	where (vegdem2==a)
        	vegdem2=0.0
    	end where
        ! Bush separation
        allocate(bush(sizex,sizey))
        where ((vegdem>0) .and. (vegdem2==0))
        	bush=vegdem
        elsewhere
        	bush=0.0
        end where
 	else
    	trans=1.00;
    endif    
    
    !!! Loading/creating SVFs !!!
    Path=trim(FileInputPath)//trim(SVFPath)//trim(SVFsuffix)
    svfname=(/'svf.asc ','svfE.asc','svfN.asc','svfW.asc','svfS.asc'/)
    svfvegname=(/'svfveg.asc  ','svfEveg.asc ','svfNveg.asc ','svfWveg.asc ','svfSveg.asc ',&
        'svfaveg.asc ','svfEaveg.asc','svfNaveg.asc','svfWaveg.asc','svfSaveg.asc'/)
    ! SVFs, Should be done as a loop... ! How to change variable in a loop???      
    call LoadEsriAsciiGrid(Path,svfname(1),xllcorner,yllcorner,cellsize,NoData)
    allocate(svf(sizey,sizex)); svf=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(2),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfE(sizey,sizex)); svfE=tempgrid; deallocate(tempgrid)        
    call LoadEsriAsciiGrid(Path,svfname(3),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfN(sizey,sizex)); svfN=tempgrid; deallocate(tempgrid)        
    call LoadEsriAsciiGrid(Path,svfname(4),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfW(sizey,sizex)); svfW=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfname(5),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfS(sizey,sizex)); svfS=tempgrid; deallocate(tempgrid)
    call LoadEsriAsciiGrid(Path,svfvegname(1),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfveg(sizey,sizex)); svfveg=tempgrid; deallocate(tempgrid)
   	call LoadEsriAsciiGrid(Path,svfvegname(2),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfEveg(sizey,sizex)); svfEveg=tempgrid; deallocate(tempgrid)        
   	call LoadEsriAsciiGrid(Path,svfvegname(3),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfNveg(sizey,sizex)); svfNveg=tempgrid; deallocate(tempgrid)        
  	call LoadEsriAsciiGrid(Path,svfvegname(4),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfWveg(sizey,sizex)); svfWveg=tempgrid; deallocate(tempgrid)
   	call LoadEsriAsciiGrid(Path,svfvegname(5),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfSveg(sizey,sizex)); svfSveg=tempgrid; deallocate(tempgrid)
   	call LoadEsriAsciiGrid(Path,svfvegname(6),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfaveg(sizey,sizex)); svfaveg=tempgrid; deallocate(tempgrid)
   	call LoadEsriAsciiGrid(Path,svfvegname(7),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfEaveg(sizey,sizex)); svfEaveg=tempgrid; deallocate(tempgrid)        
   	call LoadEsriAsciiGrid(Path,svfvegname(8),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfNaveg(sizey,sizex)) ; svfNaveg=tempgrid; deallocate(tempgrid)        
   	call LoadEsriAsciiGrid(Path,svfvegname(9),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfWaveg(sizey,sizex)); svfWaveg=tempgrid; deallocate(tempgrid)
   	call LoadEsriAsciiGrid(Path,svfvegname(10),xllcorner,yllcorner,cellsize,NoData)
    allocate(svfSaveg(sizey,sizex)); svfSaveg=tempgrid; deallocate(tempgrid)
    
    !!! Loading buildings grid !!!
    Path=trim(FileInputPath)//trim(DSMPath)
    GridFile=trim(Path)//trim(buildingsname)
    inquire(file=GridFile, exist=exist)
    if (exist) then
    	call LoadEsriAsciiGrid(Path,buildingsname,xllcorner,yllcorner,cellsize,NoData)
        allocate(buildings(sizey,sizex))
        buildings=tempgrid
        deallocate(tempgrid)
    else
        !!! Not ready, should return error. Also for the other grids
    endif
    
    ! Time related info
    !timestepdec=t_INTERVAL/(1440.*60.) 
    timeadd=0.00 
    
    ! Initiate map for surface temperature delay
    allocate(Tgmap1(sizey,sizex))
    Tgmap1=0.0
    
    return
    
274 call ErrorHint(40,trim("SOLWEIGinput.nml FileCode is missing"),notUsed,notUsed,notUsedI) 
    
end subroutine SOLWEIG_Initial