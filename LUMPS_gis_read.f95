
 SUBROUTINE read_gis(finish)
 !Reading of GIS data. In SUEWS only options 3 and 4 in use. If using other options, check the reading.
 !Changes: Calculation of z0m and zv LJ (11/2010)
 !		  Irrigation fraction modified to match with SUEWS LJ (4 Aug 2011)
 ! sg feb 2012   - surface types allocated, and array allocated
 ! --------------- extra material removed
 ! check this could be cleaned up
 ! New surface irrigated trees/shrubs added to the gis file (LJ Sep 2013)
 !-------------------------------------------------------------------------
  USE gis_data
  USE sues_data
  use time
  USE mod_z
  use gas
  use allocateArray  !module_LUMPS_constants,f90
  use defaultnotUsed

  IMPLICIT NONE

  logical::finish
      								! only used in this subroutine
  real(kind(1d0)):: build,&         !Build fraction from land area
                    cany3d,&        !Fraction of canyons
                    con_sh,&        !Coniferous
                    dec_sh,&        !deciduous
                    grassUnIrr,&    !Unirrigated grass fraction
					grassIrr,&      !Fraction of irrigated grass
                    tree_sh,&       !Tree fraction from land area
                    unman,&         !Unman fraction from land area egbare soil
                    water,&         !Water fraction from land area
                    ximper          !Paved fraction from land area
                         

  !character (len=80)::SubroutineName='Read_GIS in LUMPS_gisread.f95'
  !Only 3 and 4 are used in SUEWS (LJ Oct 2010)
 
  READ(3,*,iostat=ios_out)id,it,build,ximper,unman,con_sh,dec_sh,grassUnirr,GrassIrr,water,&
     	BldgH,TreeH,FAIBldg,FAITree,z0m,zdm,Alt,IrrFractionTrees

  IF (ios_out<0) THEN
      CLOSE (3)
      finish=.TRUE.
      RETURN
  END IF
  !tree_sh=con_sh+dec_sh !In case if needed in future
    
  sfr(PavSurf)=ximper     			! pav
  sfr(BldgSurf)=build     			! building
  sfr(ConifSurf)=con_sh       		! coniferous 
  sfr(DecidSurf)=dec_sh       		! deciduous
  sfr(GrassISurf)=grassIrr        	! irigated grass   
  sfr(GrassUSurf)=grassUnirr+Unman  ! unirrigatedgrass
  sfr(WaterSurf)=water
  
  BareSoilSurfFraction=Unman
  VegFraction= sfr(ConifSurf)+ sfr(DecidSurf)+sfr(GrassISurf)+sfr(GrassUSurf)

  areaZh=(sfr(BldgSurf)+sfr(ConifSurf)+sfr(DecidSurf)) !Total area of buildings and trees    

  !Check do the fractions add up close to one
  if(sum(sfr)>1.001.or.sum(sfr)<0.999) call ErrorHint(10,'GIS File- surface fractions',sum(sfr),notUsed,notUsedI)
      	
  
  !How vegetation is defined in LUMPS
  IF(veg_type==1)THEN         ! area vegetated
     veg_fr=con_sh+dec_sh+grassUnIrr+grassIrr+water+unman

  ELSEIF(veg_type==2)THEN     ! area irrigated
     !veg_fr=tree_sh*TreeFractionIrrigated+gras*GrassFractionIrrigated+watr+ximper*PavedFractionIrrigated
  	 veg_fr=grassIrr+IrrFractionTrees*sfr(ConifSurf)+IrrFractionTrees*sfr(DecidSurf)  !Onlu irrigated grass taken into account
  END IF

 
  RETURN

 END SUBROUTINE read_gis




