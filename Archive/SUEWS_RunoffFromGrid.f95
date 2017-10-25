subroutine RunoffFromGrid(GridFromFrac)
	!Take water flows from other grids into account
	!LJ in 10/2010
    
    use data_in
    use SUES_data
    
    IMPLICIT NONE

    !Define columns of th eoutput file (if changed, need to be taken into account)
    real(kind(1d0))::a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,&
                     a24,a25,a26,a27
    real(kind(1d0)),dimension(4)::GridFromFrac
    integer::iostat_var
  
! check - give source -- i.e. subrutine that is writing this et 

    !Initiate water body data
    addWaterBody=0
    addImpervious=0
    addVeg=0
    addPipes=0
    
    !Read the additional water from surroundings grids here. 4 possible options.
    do in=1,4

       if (GridFromFrac(in)/=0) then!If runoff from other surfaces exists read data

       		!File identifier
            if (in==1) then
                lfnOld=486
            elseif (in==2) then
                lfnOld=487
            elseif (in==3) then
                lfnOld=488
            elseif (in==4) then
                lfnOld=489
            endif
       
          !a24=runoffPipes		 !!Need to alter these to be the volumes [_m3] and divide by SurfaceArea
          !a25=runoffAGimpervious
          !a26=runoffAGveg
          !a27=runoffWaterBody

            READ(lfnOld,*,iostat=iostat_var)a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,&
                                            a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27
                          
            addWaterBody=addWaterBody+a27       
            addPipes=addPipes+a24*GridFromFrac(in)
            addImpervious=addImpervious+a25*GridFromFrac(in)
            addVeg=addVeg+a26*GridFromFrac(in)
            
            if(iostat_var<0)THEN
               iostat_var=0
               CLOSE(lfnOld)
               RETURN
            ENDIF
       endif

    enddo  

 return

 
end subroutine RunoffFromGrid
