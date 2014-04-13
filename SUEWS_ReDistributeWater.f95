subroutine ReDistributeWater
!Drainage moves into different parts defined by WaterDistSS_YYYY.txt. LJ 2010
!addWater(is) has that amount of water that is gained for each surface
!
!Latest update takes snow into account. 22/03/2013 LJ
!-------------------------------------------------------------------

  use allocateArray
  use Sues_data
  use gis_data
  use data_in
     
  implicit none
  
  integer::ii,jj
  real (kind(1d0)):: WaterVol
  
  !Fractions that go to runoff from each surface
  do ii=1,nsurf-1 ! not water in the calculation
      AddWaterRunoff(ii)=WaterDist(8,ii) 
  enddo
  AddWaterRunoff(WaterSurf)=0
  addWater=0

  do ii=1,nsurf-NSurfDoNotReceiveDrainage !go through surfaces from 1 to 7. These gain water through drainage
    do jj=1,nsurf-(NSurfDoNotReceiveDrainage+1) !From where surface ii can gain water - can't gain water from itself
       
        if (sfr(ii)/=0) then !Water movement takes place only if surface fraction exists

           !No snow calculations!
           if (snowUse==0) then
              AddWater(ii)=AddWater(ii)+(Drain(jj)*sfr(jj)/sfr(ii))*WaterDist(ii,jj) !Original

           !Snow included, This needs to be fixed at some point. LJ Mar 2013
           else
              AddWaterRunoff(jj)=AddWaterRunoff(jj)+WaterDist(ii,jj) !No receiving surface -> runoff 
           endif
                                                             
        else
           AddWaterRunoff(jj)=AddWaterRunoff(jj)+WaterDist(ii,jj) !If no receiving surface exists,
                                                                  !water fraction goes to AddWaterRunoff
        endif
    enddo
  enddo
end subroutine ReDistributeWater