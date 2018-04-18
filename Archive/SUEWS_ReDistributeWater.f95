SUBROUTINE ReDistributeWater(&
     ! input:
     nsurf,& ! surface type number
     WaterSurf,&
     snowUse,&
     WaterDist,  &
     sfr,   &!
     Drain,&
     ! output:
     AddWaterRunoff,&
     addWater&
     )
  !Drainage moves into different parts defined by WaterDistSS_YYYY.txt. LJ 2010
  !addWater(is) is that amount of water that is gained for each surface
  !Latest update takes snow into account. 22/03/2013 LJ
  !-------------------------------------------------------------------

  ! use allocateArray
  ! use data_in
  ! use gis_data
  ! use sues_data

  IMPLICIT NONE
  INTEGER,INTENT(in)::&
       nsurf,         & ! number of surface types
       WaterSurf,     &!=7, water surface code
       snowUse!Snow part used (1) or not used (0)
  REAL (KIND(1d0)),INTENT(in)::&
       WaterDist(nsurf+1,nsurf-1), &!Within-grid water distribution to other surfaces and runoff/soil store [-]
       sfr(nsurf),                 &!Surface fractions [-]
       Drain(nsurf)                 !Drainage of each surface type [mm]
  REAL (KIND(1d0)),INTENT(out)::&
       AddWaterRunoff(nsurf),&!Fraction of water going to runoff/sub-surface soil (WGWaterDist) [-]
       addWater(nsurf)        !Water from other surfaces (WGWaterDist in SUEWS_ReDistributeWater.f95) [mm]

  INTEGER::ii,jj,&
       NSurfDoNotReceiveDrainage=0!Number of surfaces that do not receive drainage water (green roof)

  !Fractions that go to runoff from each surface
  DO ii=1,nsurf-1   !not water in the calculation
     AddWaterRunoff(ii)=WaterDist(8,ii)
  ENDDO
  AddWaterRunoff(WaterSurf)=0
  addWater=0

  DO ii=1,nsurf-NSurfDoNotReceiveDrainage !go through surfaces from 1 to 7. These gain water through drainage
     DO jj=1,nsurf-(NSurfDoNotReceiveDrainage+1) !From where surface ii can gain water - can't gain water from itself

        IF (sfr(ii)/=0) THEN !Water movement takes place only if surface fraction exists

           !No snow calculations!
           IF (snowUse==0) THEN
              AddWater(ii)=AddWater(ii)+(Drain(jj)*sfr(jj)/sfr(ii))*WaterDist(ii,jj) !Original

              !Snow included, This needs to be fixed at some point. LJ Mar 2013
           ELSE
              AddWaterRunoff(jj)=AddWaterRunoff(jj)+WaterDist(ii,jj) !No receiving surface -> runoff
           ENDIF

        ELSE
           AddWaterRunoff(jj)=AddWaterRunoff(jj)+WaterDist(ii,jj) !If no receiving surface exists,
           !water fraction goes to AddWaterRunoff
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE ReDistributeWater
