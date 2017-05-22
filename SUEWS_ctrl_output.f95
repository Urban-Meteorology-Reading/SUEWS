MODULE ctrl_output
  !========================================================================================
  ! generic output functions for SUEWS
  ! authors: Ting Sun (ting.sun@reading.ac.uk)
  !
  ! disclamier:
  !     This code employs the netCDF Fortran 90 API.
  !     Full documentation of the netCDF Fortran 90 API can be found at:
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
  !     Part of the work is under the help of examples provided by the documentation.
  !
  ! purpose:
  ! these subroutines write out the results of SUEWS in netCDF format.
  !
  !
  ! history:
  ! TS 20161209: initial version of netcdf function
  ! TS 20161213: standalise the txt2nc procedure
  ! TS 20170414: generic output procedures
  !========================================================================================


  USE allocateArray
  USE cbl_module
  USE data_in
  USE defaultNotUsed
  USE ESTM_data
  USE gis_data
  USE initial
  USE solweig_module
  USE sues_data
  USE time
  USE strings

  IMPLICIT NONE


  INTEGER :: i

  CHARACTER(len=10),PARAMETER:: & !Define useful formats here
       fy   = '(i0004,1X)',& !4 digit integer for year
       ft   = '(i0004,1X)',& !3 digit integer for id, it, imin
       fd   = '(f08.4,1X)',& !3 digits + 4 dp for dectime
       f94  = '(f09.4,1X)',& !standard output format: 4 dp + 4 digits
       f104 = '(f10.4,1X)',& !standard output format: 4 dp + 5 digits
       f106 = '(f10.6,1X)',& !standard output format: 6 dp + 3 digits
       f146 = '(f14.6,1X)'   !standard output format: 6 dp + 7 digits

  CHARACTER(len= 1),PARAMETER:: & ! Define aggregation methods here
       aT = '0',&   !time columns
       aA = '1',&   !average
       aS = '2',&   !sum
       aL = '3'     !last value

  CHARACTER(len= 3):: itext

  ! define type: variable attributes
  TYPE varAttr
     CHARACTER(len = 15) :: header ! short name in headers
     CHARACTER(len = 12) :: unit   ! unit
     CHARACTER(len = 14) :: fmt    ! output format
     CHARACTER(len = 50) :: longNm ! long name for detailed description
     CHARACTER(len = 1)  :: aggreg ! aggregation method
     CHARACTER(len = 10) :: group  ! group: datetime, default, ESTM, Snow, etc.
     INTEGER             :: level  ! output priority level: 0 for highest (defualt output)
  END TYPE varAttr

  ! initialise valist
  TYPE(varAttr) :: varList(300)

  ! datetime:
  DATA(varList(i), i=1,5)/&
       varAttr('Year'    , 'YYYY' , fy , 'Year'         , aT , 'datetime' , 0),&
       varAttr('DOY'     , 'DOY'  , ft , 'Day of Year'  , aT , 'datetime' , 0),&
       varAttr('Hour'    , 'HH'   , ft , 'Hour'         , aT , 'datetime' , 0),&
       varAttr('Min'     , 'MM'   , ft , 'Minute'       , aT , 'datetime' , 0),&
       varAttr('Dectime' , '-'    , fd , 'Decimal time' , aT , 'datetime' , 0)&
       /

  ! defualt:
  DATA(varList(i), i=6,84)/&
       varAttr('Kdown'      , 'W_m-2'        , f94  , 'Incoming shortwave radiation'                     , aA , '' , 0)     , &
       varAttr('Kup'        , 'W_m-2'        , f94  , 'Outgoing shortwave radiation'                     , aA , '' , 0)     , &
       varAttr('Ldown'      , 'W_m-2'        , f94  , 'Incoming longwave radiation'                      , aA , '' , 0)     , &
       varAttr('Lup'        , 'W_m-2'        , f94  , 'Outgoing longwave radiation'                      , aA , '' , 0)     , &
       varAttr('Tsurf'      , 'degC'         , f94  , 'Bulk surface temperature'                         , aA , '' , 0)     , &
       varAttr('QN'         , 'W_m-2'        , f94  , 'Net all-wave radiation'                           , aA , '' , 0)     , &
       varAttr('QF'         , 'W_m-2'        , f94  , 'Anthropogenic heat flux'                          , aA , '' , 0)     , &
       varAttr('QS'         , 'W_m-2'        , f94  , 'Net storage heat flux'                            , aA , '' , 0)     , &
       varAttr('QH'         , 'W_m-2'        , f94  , 'Sensible heat flux'                               , aA , '' , 0)     , &
       varAttr('QE'         , 'W_m-2'        , f94  , 'Latent heat flux'                                 , aA , '' , 0)     , &
       varAttr('QHlumps'    , 'W_m-2'        , f94  , 'Sensible heat flux (using LUMPS)'                 , aA , '' , 1)     , &
       varAttr('QElumps'    , 'W_m-2'        , f94  , 'Latent heat flux (using LUMPS)'                   , aA , '' , 1)     , &
       varAttr('QHresis'    , 'W_m-2'        , f94  , 'Sensible heat flux (resistance method)'           , aA , '' , 1)     , &
       varAttr('Rain'       , 'mm'           , f106 , 'Rain'                                             , aS , '' , 0)     , &
       varAttr('Irr'        , 'mm'           , f106 , 'Irrigation'                                       , aS , '' , 0)     , &
       varAttr('Evap'       , 'mm'           , f106 , 'Evaporation'                                      , aS , '' , 0)     , &
       varAttr('RO'         , 'mm'           , f106 , 'Runoff'                                           , aS , '' , 0)     , &
       varAttr('TotCh'      , 'mm'           , f106 , 'Surface and soil moisture change'                 , aS , '' , 0)     , &
       varAttr('SurfCh'     , 'mm'           , f106 , 'Surface moisture change'                          , aS , '' , 0)     , &
       varAttr('State'      , 'mm'           , f104 , 'SurfaceWetnessState'                              , aL , '' , 0)     , &
       varAttr('NWtrState'  , 'mm'           , f106 , 'Surface wetness state (non-water surfaces)'       , aL , '' , 0)     , &
       varAttr('Drainage'   , 'mm'           , f106 , 'Drainage'                                         , aS , '' , 0)     , &
       varAttr('SMD'        , 'mm'           , f94  , 'SoilMoistureDeficit'                              , aL , '' , 0)     , &
       varAttr('FlowCh'     , 'mm'           , f104 , 'Additional flow into water body'                  , aS , '' , 1)     , &
       varAttr('AddWater'   , 'mm'           , f104 , 'Addtional water from other grids'                 , aS , '' , 1)     , &
       varAttr('ROSoil'     , 'mm'           , f106 , 'Runoff to soil'                                   , aS , '' , 1)     , &
       varAttr('ROPipe'     , 'mm'           , f106 , 'Runoff to pipes'                                  , aS , '' , 1)     , &
       varAttr('ROImp'      , 'mm'           , f106 , 'Runoff over impervious surfaces'                  , aS , '' , 1)     , &
       varAttr('ROVeg'      , 'mm'           , f106 , 'Runoff over vegetated surfaces'                   , aS , '' , 1)     , &
       varAttr('ROWater'    , 'mm'           , f106 , 'Runoff for water surface'                         , aS , '' , 1)     , &
       varAttr('WUInt'      , 'mm'           , f94  , 'InternalWaterUse'                                 , aS , '' , 1)     , &
       varAttr('WUEveTr'    , 'mm'           , f94  , 'Water use for evergreen trees'                    , aS , '' , 1)     , &
       varAttr('WUDecTr'    , 'mm'           , f94  , 'Water use for deciduous trees'                    , aS , '' , 1)     , &
       varAttr('WUGrass'    , 'mm'           , f94  , 'Water use for grass'                              , aS , '' , 1)     , &
       varAttr('SMDPaved'   , 'mm'           , f94  , 'Soil moisture deficit for paved surface'          , aL , '' , 1)     , &
       varAttr('SMDBldgs'   , 'mm'           , f94  , 'Soil moisture deficit for building surface'       , aL , '' , 1)     , &
       varAttr('SMDEveTr'   , 'mm'           , f94  , 'Soil moisture deficit for evergreen tree surface' , aL , '' , 1)     , &
       varAttr('SMDDecTr'   , 'mm'           , f94  , 'Soil moisture deficit for deciduous tree surface' , aL , '' , 1)     , &
       varAttr('SMDGrass'   , 'mm'           , f94  , 'Soil moisture deficit for grass surface'          , aL , '' , 1)     , &
       varAttr('SMDBSoil'   , 'mm'           , f94  , 'Soil moisture deficit for bare soil surface'      , aL , '' , 1)     , &
       varAttr('StPaved'    , 'mm'           , f94  , 'Surface wetness state for paved surface'          , aL , '' , 1)     , &
       varAttr('StBldgs'    , 'mm'           , f94  , 'Surface wetness state for building surface'       , aL , '' , 1)     , &
       varAttr('StEveTr'    , 'mm'           , f94  , 'Surface wetness state for evergreen tree surface' , aL , '' , 1)     , &
       varAttr('StDecTr'    , 'mm'           , f94  , 'Surface wetness state for deciduous tree surface' , aL , '' , 1)     , &
       varAttr('StGrass'    , 'mm'           , f94  , 'Surface wetness state for grass surface'          , aL , '' , 1)     , &
       varAttr('StBSoil'    , 'mm'           , f94  , 'Surface wetness state for bare soil surface'      , aL , '' , 1)     , &
       varAttr('StWater'    , 'mm'           , f104 , 'Surface wetness state for water surface'          , aL , '' , 1)     , &
       varAttr('Zenith'     , 'deg'          , f94  , 'Solar zenith angle'                               , aL , '' , 0)     , &
       varAttr('Azimuth'    , 'deg'          , f94  , 'Solar azimuth angle'                              , aL , '' , 0)     , &
       varAttr('AlbBulk'    , '-'            , f94  , 'Bulk albedo'                                      , aA , '' , 0)     , &
       varAttr('Fcld'       , '-'            , f94  , 'Cloud fraction'                                   , aA , '' , 0)     , &
       varAttr('LAI'        , 'm2_m-2'       , f94  , 'Leaf area index'                                  , aA , '' , 0)     , &
       varAttr('z0m'        , 'm'            , f94  , 'Roughness length for momentum'                    , aA , '' , 1)     , &
       varAttr('zdm'        , 'm'            , f94  , 'Zero-plane displacement height'                   , aA , '' , 1)     , &
       varAttr('ustar'      , 'm_s-1'        , f94  , 'Friction velocity'                                , aA , '' , 0)     , &
       varAttr('Lob'        , 'm'            , f104 , 'Obukhov length'                                   , aA , '' , 0)     , &
       varAttr('ra'         , 's_m-1'        , f94  , 'Aerodynamic resistance'                           , aA , '' , 1)     , &
       varAttr('rs'         , 's_m-1'        , f94  , 'Surface resistance'                               , aA , '' , 1)     , &
       varAttr('Fc'         , 'umol_m-2_s-1' , f94  , 'CO2 flux'                                         , aA , '' , 0)     , &
       varAttr('FcPhoto'    , 'umol_m-2_s-1' , f94  , 'CO2 flux from photosynthesis'                     , aA , '' , 1)     , &
       varAttr('FcRespi'    , 'umol_m-2_s-1' , f94  , 'CO2 flux from respiration'                        , aA , '' , 1)     , &
       varAttr('FcMetab'    , 'umol_m-2_s-1' , f94  , 'CO2 flux from metabolism'                         , aA , '' , 1)     , &
       varAttr('FcTraff'    , 'umol_m-2_s-1' , f94  , 'CO2 flux from traffic'                            , aA , '' , 1)     , &
       varAttr('FcBuild'    , 'umol_m-2_s-1' , f94  , 'CO2 flux from buildings'                          , aA , '' , 1)     , &
       varAttr('QNSnowFr'   , 'W_m-2'        , f94  , 'Net all-wave radiation for non-snow area'         , aA , '' , 2)     , &
       varAttr('QNSnow'     , 'W_m-2'        , f94  , 'Net all-wave radiation for snow area'             , aA , '' , 2)     , &
       varAttr('AlbSnow'    , '-'            , f94  , 'Snow albedo'                                      , aA , '' , 2)     , &
       varAttr('QM'         , 'W_m-2'        , f106 , 'Snow-related heat exchange'                       , aA , '' , 2)     , &
       varAttr('QMFreeze'   , 'W_m-2'        , f106 , 'Internal energy change'                           , aA , '' , 2)     , &
       varAttr('QMRain'     , 'W_m-2'        , f106 , 'Heat released by rain on snow'                    , aA , '' , 2)     , &
       varAttr('SWE'        , 'mm'           , f104 , 'Snow water equivalent'                            , aA , '' , 2)     , &
       varAttr('MeltWater'  , 'mm'           , f104 , 'Meltwater'                                        , aA , '' , 2)     , &
       varAttr('MeltWStore' , 'mm'           , f104 , 'Meltwater store'                                  , aA , '' , 2)     , &
       varAttr('SnowCh'     , 'mm'           , f104 , 'Change in snow pack'                              , aA , '' , 2)     , &
       varAttr('SnowRPaved' , 'mm'           , f94  , 'Snow removed from paved surface'                  , aS , '' , 2)     , &
       varAttr('SnowRBldg'  , 'mm'           , f94  , 'Snow removed from building surface'               , aS , '' , 2)     , &
       varAttr('T2'         , 'degC'         , f94  , 'Air temperature at 2 m'                           , aA , '' , 0)     , &
       varAttr('Q2'         , 'g_kg-1'       , f94  , 'Specific humidity at 2 m'                         , aA , '' , 0)     , &
       varAttr('U10'        , 'm_s-1'        , f94  , 'Wind speed at 10 m'                               , aA , '' , 0)   &

       /

  ! SOLWEIG:
  DATA(varList(i), i=85,110)/&
       varAttr('azimuth'    , 'to_add' , f106 , 'azimuth'    , aA , 'SOLWEIG' , 0)  , &
       varAttr('altitude'   , 'to_add' , f106 , 'altitude'   , aA , 'SOLWEIG' , 0)  , &
       varAttr('GlobalRad'  , 'to_add' , f106 , 'GlobalRad'  , aA , 'SOLWEIG' , 0)  , &
       varAttr('DiffuseRad' , 'to_add' , f106 , 'DiffuseRad' , aA , 'SOLWEIG' , 0)  , &
       varAttr('DirectRad'  , 'to_add' , f106 , 'DirectRad'  , aA , 'SOLWEIG' , 0)  , &
       varAttr('Kdown2d'    , 'to_add' , f106 , 'Kdown2d'    , aA , 'SOLWEIG' , 0)  , &
       varAttr('Kup2d'      , 'to_add' , f106 , 'Kup2d'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Ksouth'     , 'to_add' , f106 , 'Ksouth'     , aA , 'SOLWEIG' , 0)  , &
       varAttr('Kwest'      , 'to_add' , f106 , 'Kwest'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Knorth'     , 'to_add' , f106 , 'Knorth'     , aA , 'SOLWEIG' , 0)  , &
       varAttr('Keast'      , 'to_add' , f106 , 'Keast'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Ldown2d'    , 'to_add' , f106 , 'Ldown2d'    , aA , 'SOLWEIG' , 0)  , &
       varAttr('Lup2d'      , 'to_add' , f106 , 'Lup2d'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Lsouth'     , 'to_add' , f106 , 'Lsouth'     , aA , 'SOLWEIG' , 0)  , &
       varAttr('Lwest'      , 'to_add' , f106 , 'Lwest'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Lnorth'     , 'to_add' , f106 , 'Lnorth'     , aA , 'SOLWEIG' , 0)  , &
       varAttr('Least'      , 'to_add' , f106 , 'Least'      , aA , 'SOLWEIG' , 0)  , &
       varAttr('Tmrt'       , 'to_add' , f106 , 'Tmrt'       , aA , 'SOLWEIG' , 0)  , &
       varAttr('I0'         , 'to_add' , f106 , 'I0'         , aA , 'SOLWEIG' , 0)  , &
       varAttr('CI'         , 'to_add' , f106 , 'CI'         , aA , 'SOLWEIG' , 0)  , &
       varAttr('gvf'        , 'to_add' , f106 , 'gvf'        , aA , 'SOLWEIG' , 0)  , &
       varAttr('shadow'     , 'to_add' , f106 , 'shadow'     , aA , 'SOLWEIG' , 0)  , &
       varAttr('svf'        , 'to_add' , f106 , 'svf'        , aA , 'SOLWEIG' , 0)  , &
       varAttr('svfbuveg'   , 'to_add' , f106 , 'svfbuveg'   , aA , 'SOLWEIG' , 0)  , &
       varAttr('Ta'         , 'to_add' , f106 , 'Ta'         , aA , 'SOLWEIG' , 0)  , &
       varAttr('Tg'         , 'to_add' , f106 , 'Tg'         , aA , 'SOLWEIG' , 0)&
       /

  ! BL:
  DATA(varList(i), i=111,127)/&
       varAttr('z'         , 'to_add' , f104 , 'z'         , aA , 'BL' , 0)  , &
       varAttr('theta'     , 'to_add' , f104 , 'theta'     , aA , 'BL' , 0)  , &
       varAttr('q'         , 'to_add' , f104 , 'q'         , aA , 'BL' , 0)  , &
       varAttr('theta+'    , 'to_add' , f104 , 'theta+'    , aA , 'BL' , 0)  , &
       varAttr('q+'        , 'to_add' , f104 , 'q+'        , aA , 'BL' , 0)  , &
       varAttr('Temp_C'    , 'to_add' , f104 , 'Temp_C'    , aA , 'BL' , 0)  , &
       varAttr('rh'        , 'to_add' , f104 , 'rh'        , aA , 'BL' , 0)  , &
       varAttr('QH_use'    , 'to_add' , f104 , 'QH_use'    , aA , 'BL' , 0)  , &
       varAttr('QE_use'    , 'to_add' , f104 , 'QE_use'    , aA , 'BL' , 0)  , &
       varAttr('Press_hPa' , 'to_add' , f104 , 'Press_hPa' , aA , 'BL' , 0)  , &
       varAttr('avu1'      , 'to_add' , f104 , 'avu1'      , aA , 'BL' , 0)  , &
       varAttr('ustar'     , 'to_add' , f104 , 'ustar'     , aA , 'BL' , 0)  , &
       varAttr('avdens'    , 'to_add' , f104 , 'avdens'    , aA , 'BL' , 0)  , &
       varAttr('lv_J_kg'   , 'to_add' , f104 , 'lv_J_kg'   , aA , 'BL' , 0)  , &
       varAttr('avcp'      , 'to_add' , f104 , 'avcp'      , aA , 'BL' , 0)  , &
       varAttr('gamt'      , 'to_add' , f104 , 'gamt'      , aA , 'BL' , 0)  , &
       varAttr('gamq'      , 'to_add' , f104 , 'gamq'      , aA , 'BL' , 0)&
       /

  ! Snow:
  DATA(varList(i), i=128,224)/&
       varAttr('SWE_Paved'      , 'to_add' , f106 , 'SWE_Paved'      , aA , 'snow' , 0)  , &
       varAttr('SWE_Bldgs'      , 'to_add' , f106 , 'SWE_Bldgs'      , aA , 'snow' , 0)  , &
       varAttr('SWE_EveTr'      , 'to_add' , f106 , 'SWE_EveTr'      , aA , 'snow' , 0)  , &
       varAttr('SWE_DecTr'      , 'to_add' , f106 , 'SWE_DecTr'      , aA , 'snow' , 0)  , &
       varAttr('SWE_Grass'      , 'to_add' , f106 , 'SWE_Grass'      , aA , 'snow' , 0)  , &
       varAttr('SWE_BSoil'      , 'to_add' , f106 , 'SWE_BSoil'      , aA , 'snow' , 0)  , &
       varAttr('SWE_Water'      , 'to_add' , f106 , 'SWE_Water'      , aA , 'snow' , 0)  , &
       varAttr('Mw_Paved'       , 'to_add' , f106 , 'Mw_Paved'       , aA , 'snow' , 0)  , &
       varAttr('Mw_Bldgs'       , 'to_add' , f106 , 'Mw_Bldgs'       , aA , 'snow' , 0)  , &
       varAttr('Mw_EveTr'       , 'to_add' , f106 , 'Mw_EveTr'       , aA , 'snow' , 0)  , &
       varAttr('Mw_DecTr'       , 'to_add' , f106 , 'Mw_DecTr'       , aA , 'snow' , 0)  , &
       varAttr('Mw_Grass'       , 'to_add' , f106 , 'Mw_Grass'       , aA , 'snow' , 0)  , &
       varAttr('Mw_BSoil'       , 'to_add' , f106 , 'Mw_BSoil'       , aA , 'snow' , 0)  , &
       varAttr('Mw_Water'       , 'to_add' , f106 , 'Mw_Water'       , aA , 'snow' , 0)  , &
       varAttr('Qm_Paved'       , 'to_add' , f106 , 'Qm_Paved'       , aA , 'snow' , 0)  , &
       varAttr('Qm_Bldgs'       , 'to_add' , f106 , 'Qm_Bldgs'       , aA , 'snow' , 0)  , &
       varAttr('Qm_EveTr'       , 'to_add' , f106 , 'Qm_EveTr'       , aA , 'snow' , 0)  , &
       varAttr('Qm_DecTr'       , 'to_add' , f106 , 'Qm_DecTr'       , aA , 'snow' , 0)  , &
       varAttr('Qm_Grass'       , 'to_add' , f106 , 'Qm_Grass'       , aA , 'snow' , 0)  , &
       varAttr('Qm_BSoil'       , 'to_add' , f106 , 'Qm_BSoil'       , aA , 'snow' , 0)  , &
       varAttr('Qm_Water'       , 'to_add' , f106 , 'Qm_Water'       , aA , 'snow' , 0)  , &
       varAttr('Qa_Paved'       , 'to_add' , f106 , 'Qa_Paved'       , aA , 'snow' , 0)  , &
       varAttr('Qa_Bldgs'       , 'to_add' , f106 , 'Qa_Bldgs'       , aA , 'snow' , 0)  , &
       varAttr('Qa_EveTr'       , 'to_add' , f106 , 'Qa_EveTr'       , aA , 'snow' , 0)  , &
       varAttr('Qa_DecTr'       , 'to_add' , f106 , 'Qa_DecTr'       , aA , 'snow' , 0)  , &
       varAttr('Qa_Grass'       , 'to_add' , f106 , 'Qa_Grass'       , aA , 'snow' , 0)  , &
       varAttr('Qa_BSoil'       , 'to_add' , f106 , 'Qa_BSoil'       , aA , 'snow' , 0)  , &
       varAttr('Qa_Water'       , 'to_add' , f106 , 'Qa_Water'       , aA , 'snow' , 0)  , &
       varAttr('QmFr_Paved'     , 'to_add' , f106 , 'QmFr_Paved'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_Bldgs'     , 'to_add' , f106 , 'QmFr_Bldgs'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_EveTr'     , 'to_add' , f106 , 'QmFr_EveTr'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_DecTr'     , 'to_add' , f106 , 'QmFr_DecTr'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_Grass'     , 'to_add' , f106 , 'QmFr_Grass'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_BSoil'     , 'to_add' , f106 , 'QmFr_BSoil'     , aA , 'snow' , 0)  , &
       varAttr('QmFr_Water'     , 'to_add' , f106 , 'QmFr_Water'     , aA , 'snow' , 0)  , &
       varAttr('fr_Paved'       , 'to_add' , f106 , 'fr_Paved'       , aA , 'snow' , 0)  , &
       varAttr('fr_Bldgs'       , 'to_add' , f106 , 'fr_Bldgs'       , aA , 'snow' , 0)  , &
       varAttr('fr_EveTr'       , 'to_add' , f106 , 'fr_EveTr'       , aA , 'snow' , 0)  , &
       varAttr('fr_DecTr'       , 'to_add' , f106 , 'fr_DecTr'       , aA , 'snow' , 0)  , &
       varAttr('fr_Grass'       , 'to_add' , f106 , 'fr_Grass'       , aA , 'snow' , 0)  , &
       varAttr('fr_BSoil'       , 'to_add' , f106 , 'fr_BSoil'       , aA , 'snow' , 0)  , &
       varAttr('RainSn_Paved'   , 'to_add' , f146 , 'RainSn_Paved'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_Bldgs'   , 'to_add' , f146 , 'RainSn_Bldgs'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_EveTr'   , 'to_add' , f146 , 'RainSn_EveTr'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_DecTr'   , 'to_add' , f146 , 'RainSn_DecTr'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_Grass'   , 'to_add' , f146 , 'RainSn_Grass'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_BSoil'   , 'to_add' , f146 , 'RainSn_BSoil'   , aA , 'snow' , 0)  , &
       varAttr('RainSn_Water'   , 'to_add' , f146 , 'RainSn_Water'   , aA , 'snow' , 0)  , &
       varAttr('Qn_PavedSnow'   , 'to_add' , f146 , 'Qn_PavedSnow'   , aA , 'snow' , 0)  , &
       varAttr('Qn_BldgsSnow'   , 'to_add' , f146 , 'Qn_BldgsSnow'   , aA , 'snow' , 0)  , &
       varAttr('Qn_EveTrSnpw'   , 'to_add' , f146 , 'Qn_EveTrSnpw'   , aA , 'snow' , 0)  , &
       varAttr('Qn_DecTrSnow'   , 'to_add' , f146 , 'Qn_DecTrSnow'   , aA , 'snow' , 0)  , &
       varAttr('Qn_GrassSnpw'   , 'to_add' , f146 , 'Qn_GrassSnpw'   , aA , 'snow' , 0)  , &
       varAttr('Qn_BSoilSnow'   , 'to_add' , f146 , 'Qn_BSoilSnow'   , aA , 'snow' , 0)  , &
       varAttr('Qn_WaterSnow'   , 'to_add' , f146 , 'Qn_WaterSnow'   , aA , 'snow' , 0)  , &
       varAttr('kup_PavedSnow'  , 'to_add' , f146 , 'kup_PavedSnow'  , aA , 'snow' , 0)  , &
       varAttr('kup_BldgsSnow'  , 'to_add' , f146 , 'kup_BldgsSnow'  , aA , 'snow' , 0)  , &
       varAttr('kup_EveTrSnpw'  , 'to_add' , f146 , 'kup_EveTrSnpw'  , aA , 'snow' , 0)  , &
       varAttr('kup_DecTrSnow'  , 'to_add' , f146 , 'kup_DecTrSnow'  , aA , 'snow' , 0)  , &
       varAttr('kup_GrassSnpw'  , 'to_add' , f146 , 'kup_GrassSnpw'  , aA , 'snow' , 0)  , &
       varAttr('kup_BSoilSnow'  , 'to_add' , f146 , 'kup_BSoilSnow'  , aA , 'snow' , 0)  , &
       varAttr('kup_WaterSnow'  , 'to_add' , f146 , 'kup_WaterSnow'  , aA , 'snow' , 0)  , &
       varAttr('frMelt_Paved'   , 'to_add' , f146 , 'frMelt_Paved'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_Bldgs'   , 'to_add' , f146 , 'frMelt_Bldgs'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_EveTr'   , 'to_add' , f146 , 'frMelt_EveTr'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_DecTr'   , 'to_add' , f146 , 'frMelt_DecTr'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_Grass'   , 'to_add' , f146 , 'frMelt_Grass'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_BSoil'   , 'to_add' , f146 , 'frMelt_BSoil'   , aA , 'snow' , 0)  , &
       varAttr('frMelt_Water'   , 'to_add' , f146 , 'frMelt_Water'   , aA , 'snow' , 0)  , &
       varAttr('MwStore_Paved'  , 'to_add' , f146 , 'MwStore_Paved'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_Bldgs'  , 'to_add' , f146 , 'MwStore_Bldgs'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_EveTr'  , 'to_add' , f146 , 'MwStore_EveTr'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_DecTr'  , 'to_add' , f146 , 'MwStore_DecTr'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_Grass'  , 'to_add' , f146 , 'MwStore_Grass'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_BSoil'  , 'to_add' , f146 , 'MwStore_BSoil'  , aA , 'snow' , 0)  , &
       varAttr('MwStore_Water'  , 'to_add' , f146 , 'MwStore_Water'  , aA , 'snow' , 0)  , &
       varAttr('SnowDens_Paved' , 'to_add' , f146 , 'SnowDens_Paved' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_Bldgs' , 'to_add' , f146 , 'SnowDens_Bldgs' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_EveTr' , 'to_add' , f146 , 'SnowDens_EveTr' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_DecTr' , 'to_add' , f146 , 'SnowDens_DecTr' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_Grass' , 'to_add' , f146 , 'SnowDens_Grass' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_BSoil' , 'to_add' , f146 , 'SnowDens_BSoil' , aA , 'snow' , 0)  , &
       varAttr('SnowDens_Water' , 'to_add' , f146 , 'SnowDens_Water' , aA , 'snow' , 0)  , &
       varAttr('Sd_Paved'       , 'to_add' , f106 , 'Sd_Paved'       , aA , 'snow' , 0)  , &
       varAttr('Sd_Bldgs'       , 'to_add' , f106 , 'Sd_Bldgs'       , aA , 'snow' , 0)  , &
       varAttr('Sd_EveTr'       , 'to_add' , f106 , 'Sd_EveTr'       , aA , 'snow' , 0)  , &
       varAttr('Sd_DecTr'       , 'to_add' , f106 , 'Sd_DecTr'       , aA , 'snow' , 0)  , &
       varAttr('Sd_Grass'       , 'to_add' , f106 , 'Sd_Grass'       , aA , 'snow' , 0)  , &
       varAttr('Sd_BSoil'       , 'to_add' , f106 , 'Sd_BSoil'       , aA , 'snow' , 0)  , &
       varAttr('Sd_Water'       , 'to_add' , f106 , 'Sd_Water'       , aA , 'snow' , 0)  , &
       varAttr('Tsnow_Paved'    , 'to_add' , f146 , 'Tsnow_Paved'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_Bldgs'    , 'to_add' , f146 , 'Tsnow_Bldgs'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_EveTr'    , 'to_add' , f146 , 'Tsnow_EveTr'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_DecTr'    , 'to_add' , f146 , 'Tsnow_DecTr'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_Grass'    , 'to_add' , f146 , 'Tsnow_Grass'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_BSoil'    , 'to_add' , f146 , 'Tsnow_BSoil'    , aA , 'snow' , 0)  , &
       varAttr('Tsnow_Water'    , 'to_add' , f146 , 'Tsnow_Water'    , aA , 'snow' , 0)&
       /

! ESTM:
  DATA(varList(i), i=225,246)/&
       varAttr('QSIBLD'   , 'W_m-2' , f106 , 'QSIBLD'   , aA , 'Storage Internal building'                , 0)   , &
       varAttr('TWALL1'   , 'degK'  , f106 , 'TWALL1'   , aA , 'Temperature in wall layer 1'              , 0)   , &
       varAttr('TWALL2'   , 'degK'  , f106 , 'TWALL2'   , aA , 'Temperature in wall layer 2'              , 0)   , &
       varAttr('TWALL3'   , 'degK'  , f106 , 'TWALL3'   , aA , 'Temperature in wall layer 3'              , 0)   , &
       varAttr('TWALL4'   , 'degK'  , f106 , 'TWALL4'   , aA , 'Temperature in wall layer 4'              , 0)   , &
       varAttr('TWALL5'   , 'degK'  , f106 , 'TWALL5'   , aA , 'Temperature in wall layer 5'              , 0)   , &
       varAttr('TROOF1'   , 'degK'  , f106 , 'TROOF1'   , aA , 'Temperature in roof layer 1'              , 0)   , &
       varAttr('TROOF2'   , 'degK'  , f106 , 'TROOF2'   , aA , 'Temperature in roof layer 2'              , 0)   , &
       varAttr('TROOF3'   , 'degK'  , f106 , 'TROOF3'   , aA , 'Temperature in roof layer 3'              , 0)   , &
       varAttr('TROOF4'   , 'degK'  , f106 , 'TROOF4'   , aA , 'Temperature in roof layer 4'              , 0)   , &
       varAttr('TROOF5'   , 'degK'  , f106 , 'TROOF5'   , aA , 'Temperature in roof layer 5'              , 0)   , &
       varAttr('TGROUND1' , 'degK'  , f106 , 'TGROUND1' , aA , 'Temperature in ground layer 1'            , 0)   , &
       varAttr('TGROUND2' , 'degK'  , f106 , 'TGROUND2' , aA , 'Temperature in ground layer 2'            , 0)   , &
       varAttr('TGROUND3' , 'degK'  , f106 , 'TGROUND3' , aA , 'Temperature in ground layer 3'            , 0)   , &
       varAttr('TGROUND4' , 'degK'  , f106 , 'TGROUND4' , aA , 'Temperature in ground layer 4'            , 0)   , &
       varAttr('TGROUND5' , 'degK'  , f106 , 'TGROUND5' , aA , 'Temperature in ground layer 5'            , 0)   , &
       varAttr('TiBLD1'   , 'degK'  , f106 , 'TiBLD1'   , aA , 'Temperature in internal building layer 1' , 0)   , &
       varAttr('TiBLD2'   , 'degK'  , f106 , 'TiBLD2'   , aA , 'Temperature in internal building layer 2' , 0)   , &
       varAttr('TiBLD3'   , 'degK'  , f106 , 'TiBLD3'   , aA , 'Temperature in internal building layer 3' , 0)   , &
       varAttr('TiBLD4'   , 'degK'  , f106 , 'TiBLD4'   , aA , 'Temperature in internal building layer 4' , 0)   , &
       varAttr('TiBLD5'   , 'degK'  , f106 , 'TiBLD5'   , aA , 'Temperature in internal building layer 5' , 0)   , &
       varAttr('TaBLD'    , 'degK'  , f106 , 'TaBLD'    , aA , 'Indoor air temperature'                   , 0) &
       /

CONTAINS

  ! main output wrapper function
  SUBROUTINE SUEWS_Output_txt(iv,irMax,Gridiv)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: iv,irMax,Gridiv

    INTEGER :: xx,err,outLevel,i
    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varlistX
    CHARACTER(len=10) :: grpList0(5)
    CHARACTER(len=10),DIMENSION(:),ALLOCATABLE :: grpList

    ! determine outLevel
    SELECT CASE (WriteOutOption)
    CASE (0) !all (not snow-related)
       outLevel=1
    CASE (1) !all plus snow-related
       outLevel=2
    CASE (2) !minimal output
       outLevel=0
    END SELECT


    ! determine groups to output
    grpList0(1)=''
    grpList0(2)='SOLWEIG'
    grpList0(3)='BL'
    grpList0(4)='snow'
    grpList0(5)='ESTM'
    xx=COUNT((/.TRUE.,&
         SOLWEIGpoi_out==1,&
         CBLuse>=1,&
         SnowUse>=1,&
         StorageHeatMethod==4 .OR. StorageHeatMethod==14/))

    ! PRINT*, grpList0,xx

    ALLOCATE(grpList(xx), stat=err)
    IF ( err/= 0) PRINT *, "grpList: Allocation request denied"

    grpList=PACK(grpList0, &
         mask=((/.TRUE.,&
         SOLWEIGpoi_out==1,&
         CBLuse>=1,&
         SnowUse>=1,&
         StorageHeatMethod==4 .OR. StorageHeatMethod==14/)))

    ! PRINT*, grpList,SIZE(grpList, dim=1)

    ! loop over all groups
    DO i = 1, SIZE(grpList),1
       !PRINT*, 'i',i
       xx=COUNT(varlist%group == TRIM(grpList(i)), dim=1)
       !  PRINT*, 'number of variables:',xx
       ALLOCATE(varlistX(5+xx), stat=err)
       IF ( err/= 0) PRINT *, "varlistX: Allocation request denied"
       ! datetime
       varlistX(1:5)=varlist(1:5)
       ! variable
       varlistX(6:5+xx)=PACK(varlist, mask=(varlist%group == TRIM(grpList(i))))



       ! all output frequency option:
       ! as forcing:
       IF ( ResolutionFilesOut == Tstep .OR. KeepTstepFilesOut == 1 ) THEN
          CALL SUEWS_Output_txt_grp(iv,irMax,varlistX,Gridiv,outLevel,Tstep)
       ENDIF
       !  as specified ResolutionFilesOut:
       IF ( ResolutionFilesOut /= Tstep ) THEN
          CALL SUEWS_Output_txt_grp(iv,irMax,varlistX,Gridiv,outLevel,ResolutionFilesOut)
       ENDIF

       IF (ALLOCATED(varlistX)) DEALLOCATE(varlistX, stat=err)
       IF ( err/= 0) PRINT *, "varlistX: Deallocation request denied"
       !  PRINT*, 'i',i,'end'

    END DO
  END SUBROUTINE SUEWS_Output_txt


  ! output wrapper function for one group
  SUBROUTINE SUEWS_Output_txt_grp(iv,irMax,varlist,Gridiv,outLevel,outFreq_s)
    IMPLICIT NONE

    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: iv,irMax,Gridiv,outLevel,outFreq_s

    REAL(KIND(1d0))::dataOutX(irMax,SIZE(varlist))
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::dataOutX_agg


    ! determine dataout array according to variable group
    SELECT CASE (TRIM(varlist(SIZE(varlist))%group))
    CASE ('') !default
       dataOutX=dataout(1:irMax,1:SIZE(varlist),Gridiv)

    CASE ('SOLWEIG') !SOLWEIG
       ! todo: inconsistent data structure
       dataOutX=dataOutSOL(1:irMax,1:SIZE(varlist),Gridiv)

    CASE ('BL') !BL
       dataOutX=dataOutBL(1:irMax,1:SIZE(varlist),Gridiv)

    CASE ('snow')    !snow
       dataOutX=dataOutSnow(1:irMax,1:SIZE(varlist),Gridiv)

    CASE ('ESTM')    !ESTM
       dataOutX=dataOutESTM(1:irMax,1:SIZE(varlist),Gridiv)

    END SELECT

    ! PRINT*, 'n of varlistX: ',SIZE(varlist)
    ! PRINT*, 'varlistX: ',varlist%header
    ! PRINT*, 'varlistX group: ',varlist%group


    ! aggregation:
    CALL SUEWS_Output_Agg(dataOutX_agg,dataOutX,varlist,irMax,outFreq_s)

    ! output:
    ! initialise file when processing first metblock
    IF ( iv == 1 ) CALL SUEWS_Output_Init(dataOutX_agg,varlist,Gridiv,outLevel)

    ! append the aggregated data to the specific txt file
    CALL SUEWS_Write_txt(dataOutX_agg,varlist,Gridiv,outLevel)

  END SUBROUTINE SUEWS_Output_txt_grp

  ! initialise an output file with file name and headers
  SUBROUTINE SUEWS_Output_Init(dataOut,varlist,Gridiv,outLevel)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::dataOut
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: Gridiv,outLevel

    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varlistSel
    INTEGER :: xx,err,fn,i,nargs
    CHARACTER(len=100) :: FileOut
    CHARACTER(len=3) :: itext
    CHARACTER(len=6) :: args(5)
    CHARACTER(len=16*SIZE(varlist)) :: FormatOut
    CHARACTER(len=16) :: formatX,headerX
    CHARACTER(len=16), DIMENSION(:), ALLOCATABLE:: headerOut

    ! select variables to output
    xx=COUNT((varList%level<= outLevel), dim=1)
    WRITE(itext,'(i3)') xx
    ALLOCATE(varlistSel(xx), stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Allocation request denied"
    varlistSel=PACK(varList, mask=(varList%level<= outLevel))


    ! generate file name
    CALL filename_gen(dataOut,varList,Gridiv,FileOut)

    ! store right-aligned headers
    ALLOCATE(headerOut(xx), stat=err)
    IF ( err/= 0) PRINT *, "headerOut: Allocation request denied"

    ! create format string:
    DO i = 1, SIZE(varlistSel)
       CALL parse(varlistSel(i)%fmt,'if.,',args,nargs)
       formatX=ADJUSTL('(a'//TRIM(args(2))//',1x)')
       ! adjust headers to right-aligned
       WRITE(headerOut(i),formatX) ADJUSTR(TRIM(ADJUSTL(varlistSel(i)%header)))
       IF ( i==1 ) THEN
          FormatOut=ADJUSTL(TRIM(formatX))
       ELSE
          FormatOut=TRIM(FormatOut)//' '//ADJUSTL(TRIM(formatX))
       END IF
    END DO
    FormatOut='('//TRIM(ADJUSTL(FormatOut))//')'

    ! create file
    fn=9
    OPEN(fn,file=TRIM(ADJUSTL(FileOut)),status='unknown')

    ! write out headers
    WRITE(fn, FormatOut) headerOut
    CLOSE(fn)

    ! write out format file
    CALL formatFile_gen(dataOut,varlist,Gridiv,outLevel)

    ! clean up
    IF (ALLOCATED(varlistSel)) DEALLOCATE(varlistSel, stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Deallocation request denied"
    IF (ALLOCATED(headerOut)) DEALLOCATE(headerOut, stat=err)
    IF ( err/= 0) PRINT *, "headerOut: Deallocation request denied"

  END SUBROUTINE SUEWS_Output_Init

  ! generate output format file
  SUBROUTINE formatFile_gen(dataOut,varlist,Gridiv,outLevel)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::dataOut
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: Gridiv,outLevel

    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varlistSel
    INTEGER :: xx,err,fn
    CHARACTER(len=100) :: FileOut
    CHARACTER(len=50*300) :: str_cat
    CHARACTER(len=50) :: str_x=''
    CHARACTER(len=3) :: itext

    ! get filename
    CALL filename_gen(dataOut,varList,Gridiv,FileOut,1)

    !select variables to output
    xx=COUNT((varList%level<= outLevel), dim=1)
    ALLOCATE(varlistSel(xx), stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Allocation request denied"
    varlistSel=PACK(varList, mask=(varList%level<= outLevel))

    ! create file
    fn=9
    OPEN(fn,file=TRIM(ADJUSTL(FileOut)),status='unknown')

    ! write out format strings
    ! column number:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       WRITE(itext,'(i3)') i
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(itext))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(itext)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! header:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       str_x=varlistSel(i)%header
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(str_x))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(str_x)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! long name:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       str_x=varlistSel(i)%longNm
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(str_x))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(str_x)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! unit:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       str_x=varlistSel(i)%unit
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(str_x))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(str_x)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! format:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       str_x=varlistSel(i)%fmt
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(str_x))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(str_x)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! aggregation method:
    str_cat=''
    DO i = 1, SIZE(varlistSel)
       str_x=varlistSel(i)%aggreg
       IF ( i==1 ) THEN
          str_cat=TRIM(ADJUSTL(str_x))
       ELSE
          str_cat=TRIM(str_cat)//';'//ADJUSTL(str_x)
       ENDIF
    END DO
    WRITE(fn,'(a)') TRIM(str_cat)

    ! close file
    CLOSE(fn)

    ! clean up
    IF (ALLOCATED(varlistSel)) DEALLOCATE(varlistSel, stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Deallocation request denied"

  END SUBROUTINE formatFile_gen

  ! aggregate data to specified resolution
  SUBROUTINE SUEWS_Output_Agg(dataOut_agg,dataOut,varlist,irMax,outFreq_s)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::dataOut
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: irMax,outFreq_s
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,INTENT(out)::dataOut_agg

    INTEGER ::  nlinesOut,i,j,x
    REAL(KIND(1d0))::dataOut_aggX(1:SIZE(varlist))
    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::dataOut_agg0
    nlinesOut=INT(nsh/(60.*60/outFreq_s))
    ! nGrid=SIZE(dataOut, dim=3)

    ALLOCATE(dataOut_agg(INT(irMax/nlinesOut),SIZE(varlist)))
    ALLOCATE(dataOut_agg0(nlinesOut,SIZE(varlist)))


    DO i=nlinesOut,irMax,nlinesOut
       x=i/nlinesOut
       dataOut_agg0=dataOut(i-nlinesOut+1:i,:)
       DO j = 1, SIZE(varlist), 1
          ! aggregating different variables
          SELECT CASE (varlist(j)%aggreg)
          CASE (aT) !time columns, aT
             dataOut_aggX(j)=dataOut_agg0(nlinesOut,j)
          CASE (aA) !average, aA
             dataOut_aggX(j)=SUM(dataOut_agg0(:,j))/nlinesOut
          CASE (aS) !sum, aS
             dataOut_aggX(j)=SUM(dataOut_agg0(:,j))
          CASE (aL) !last value, aL
             dataOut_aggX(j)=dataOut_agg0(nlinesOut,j)
          END SELECT

          IF ( Diagnose==1 .AND. i==irMax ) THEN
             ! IF ( i==irMax ) THEN
             PRINT*, 'raw data of ',j,':'
             PRINT*, dataOut_agg0(:,j)
             PRINT*, 'aggregated with method: ',varlist(j)%aggreg
             PRINT*, dataOut_aggX(j)
             PRINT*, ''
          END IF
       END DO
       dataOut_agg(x,:)=dataOut_aggX
    END DO

  END SUBROUTINE SUEWS_Output_Agg


  ! append output data to the specific file at the specified outLevel
  SUBROUTINE SUEWS_Write_txt(dataOut,varlist,Gridiv,outLevel)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::dataOut
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: Gridiv,outLevel

    REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE::dataOutSel
    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varlistSel
    CHARACTER(len=100) :: FileOut
    INTEGER :: fn,i,xx,err
    CHARACTER(len=12*SIZE(varlist)) :: FormatOut

    !select variables to output
    xx=COUNT((varList%level<= outLevel), dim=1)
    ALLOCATE(varlistSel(xx), stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Allocation request denied"
    varlistSel=PACK(varList, mask=(varList%level<= outLevel))

    ! copy data accordingly
    ALLOCATE(dataOutSel(SIZE(dataOut, dim=1),xx), stat=err)
    IF ( err/= 0) PRINT *, "dataOutSel: Allocation request denied"
    ! print*, SIZE(varList%level),PACK((/(i,i=1,SIZE(varList%level))/), varList%level <= outLevel)
    ! print*, irMax,shape(dataOut)
    dataOutSel=dataOut(:,PACK((/(i,i=1,SIZE(varList%level))/), varList%level <= outLevel))


    ! create format string:
    DO i = 1, SIZE(varlistSel)
       IF ( i==1 ) THEN
          FormatOut=ADJUSTL(varlistSel(i)%fmt)
       ELSE
          FormatOut=TRIM(FormatOut)//' '//ADJUSTL(varlistSel(i)%fmt)
       END IF
    END DO
    FormatOut='('//TRIM(ADJUSTL(FormatOut))//')'

    ! get filename
    CALL filename_gen(dataOutSel,varlistSel,Gridiv,FileOut)

    ! write out data
    fn=50
    OPEN(fn,file=TRIM(fileout),position='append')!,err=112)
    DO i=1,SIZE(dataOutSel,dim=1)
       WRITE(fn,FormatOut) &
            INT(dataOutSel(i,1:4)),&
            dataOutSel(i,5:SIZE(varlistSel))
    ENDDO
    CLOSE (fn)

    IF (ALLOCATED(varlistSel)) DEALLOCATE(varlistSel, stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Deallocation request denied"

    IF (ALLOCATED(dataOutSel)) DEALLOCATE(dataOutSel, stat=err)
    IF ( err/= 0) PRINT *, "dataOutSel: Deallocation request denied"

  END SUBROUTINE SUEWS_Write_txt


  SUBROUTINE filename_gen(dataOut,varList,Gridiv,FileOut,opt_fmt)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:),INTENT(in)::dataOut ! to determine year & output frequency
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist ! to determine output group
    INTEGER,INTENT(in) :: Gridiv ! to determine grid name as in SiteSelect
    INTEGER,INTENT(in),OPTIONAL :: opt_fmt ! to determine if a format file
    CHARACTER(len=100),INTENT(out) :: FileOut ! the output file name

    CHARACTER(len=10):: str_out_min,str_grid,&
         str_date,str_year,str_DOY,str_grp,str_sfx
    INTEGER :: year_int,DOY_int,val_fmt

    IF( PRESENT(opt_fmt) ) val_fmt = opt_fmt

    ! date:
    year_int=INT(dataOut(1,1))
    DOY_int=INT(dataOut(1,2))
    WRITE(str_year,'(i4)') year_int
    WRITE(str_DOY,'(i3.3)') DOY_int
    str_date='_'//TRIM(ADJUSTL(str_year))
#ifdef nc
    ! add DOY as a specifier
    IF (ncMode==1) str_date=TRIM(ADJUSTL(str_date))//TRIM(ADJUSTL(str_DOY))
#endif


    ! output frequency in minute:
    WRITE(str_out_min,'(i4)') &
         INT(dataOut(2,3)-dataOut(1,3))*60& ! hour
         +INT(dataOut(2,4)-dataOut(1,4))     !minute
    str_out_min='_'//TRIM(ADJUSTL(str_out_min))

    ! group: output type
    str_grp=varList(6)%group
    IF ( LEN(TRIM(str_grp)) > 0 ) str_grp='_'//TRIM(ADJUSTL(str_grp))

    ! grid name:
    WRITE(str_grid,'(i10)') GridIDmatrix(Gridiv)
#ifdef nc
    IF (ncMode==1) str_grid='' ! grid name not needed by nc files
#endif

    ! suffix:
    str_sfx='.txt'
#ifdef nc
    IF ( ncMode==1 ) str_sfx='.nc'
#endif

    ! filename: fileout
    FileOut=TRIM(FileOutputPath)//&
         TRIM(FileCode)//&
         TRIM(ADJUSTL(str_grid))//&
         TRIM(ADJUSTL(str_date))//&
         TRIM(ADJUSTL(str_grp))//&
         TRIM(ADJUSTL(str_out_min))//&
         TRIM(ADJUSTL(str_sfx))

    ! filename: format
    IF ( val_fmt==1 ) THEN
       FileOut=TRIM(FileOutputPath)//&
            TRIM(FileCode)//&
            TRIM(ADJUSTL(str_grp))//&
            '_OutputFormat.txt'
    END IF


  END SUBROUTINE filename_gen

  !========================================================================================
  ! netCDF conversion subroutines for SUEWS
  ! author: Ting Sun
  !
  ! disclamier:
  !     This code employs the netCDF Fortran 90 API.
  !     Full documentation of the netCDF Fortran 90 API can be found at:
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
  !     Part of the work is under the help of examples provided by the documentation.
  !
  ! purpose:
  ! these subroutines write out the results of SUEWS in netCDF format.
  !
  !
  ! history:
  ! 20161209: initial version
  ! 20161213: standalise the txt2nc procedure
  !========================================================================================


  !===========================================================================!
  ! write the output of final SUEWS results in netCDF
  !   with spatial layout of QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! the output file frequency is the same as metblocks in the main SUEWS loop
  !===========================================================================!

#ifdef nc
  SUBROUTINE SUEWS_Output_nc(irMax)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: irMax

    INTEGER :: xx,err,outLevel
    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varlistX
    CHARACTER(len=10) :: grpList0(5)
    CHARACTER(len=10),DIMENSION(:),ALLOCATABLE :: grpList

    ! determine outLevel
    SELECT CASE (WriteOutOption)
    CASE (0) !all (not snow-related)
       outLevel=1
    CASE (1) !all plus snow-related
       outLevel=2
    CASE (2) !minimal output
       outLevel=0
    END SELECT


    ! determine groups to output
    ! todo: needs to be smarter, automate this filtering
    grpList0(1)=''
    grpList0(2)='SOLWEIG'
    grpList0(3)='BL'
    grpList0(4)='snow'
    grpList0(5)='ESTM'
    xx=COUNT((/.TRUE.,&
         SOLWEIGpoi_out==1,&
         CBLuse>=1,&
         SnowUse>=1,&
         StorageHeatMethod==4 .OR. StorageHeatMethod==14/))

    ! PRINT*, grpList0,xx

    ALLOCATE(grpList(xx), stat=err)
    IF ( err/= 0) PRINT *, "grpList: Allocation request denied"

    grpList=PACK(grpList0, &
         mask=((/.TRUE.,&
         SOLWEIGpoi_out==1,&
         CBLuse>=1,&
         SnowUse>=1,&
         StorageHeatMethod==4 .OR. StorageHeatMethod==14/)))

    ! PRINT*, grpList

    ! loop over all groups
    DO i = 1, SIZE(grpList)
       xx=COUNT(varlist%group == TRIM(grpList(i)), dim=1)
       !  PRINT*, 'number of variables:',xx
       ALLOCATE(varlistX(5+xx), stat=err)
       IF ( err/= 0) PRINT *, "varlistX: Allocation request denied"
       ! datetime
       varlistX(1:5)=varlist(1:5)
       ! variable
       varlistX(6:5+xx)=PACK(varlist, mask=(varlist%group == TRIM(grpList(i))))

       ! all output frequency option:
       ! as forcing:
       IF ( ResolutionFilesOut == Tstep .OR. KeepTstepFilesOut == 1 ) THEN
          CALL SUEWS_Output_nc_grp(irMax,varlistX,outLevel,Tstep)
       ENDIF
       !  as specified ResolutionFilesOut:
       IF ( ResolutionFilesOut /= Tstep ) THEN
          CALL SUEWS_Output_nc_grp(irMax,varlistX,outLevel,ResolutionFilesOut)
       ENDIF
       IF (ALLOCATED(varlistX)) DEALLOCATE(varlistX, stat=err)
       IF ( err/= 0) PRINT *, "varlistX: Deallocation request denied"
    END DO

  END SUBROUTINE SUEWS_Output_nc


  SUBROUTINE SUEWS_Output_nc_grp(irMax,varlist,outLevel,outFreq_s)
    IMPLICIT NONE

    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: irMax,outLevel,outFreq_s

    REAL(KIND(1d0))::dataOutX(irMax,SIZE(varlist),NumberOfGrids)
    REAL(KIND(1d0)),ALLOCATABLE::dataOutX_agg(:,:,:),dataOutX_agg0(:,:)
    INTEGER :: iGrid,err


    ! determine dataout array according to variable group
    SELECT CASE (TRIM(varlist(SIZE(varlist))%group))
    CASE ('') !default
       dataOutX=dataout(1:irMax,1:SIZE(varlist),:)

    CASE ('SOLWEIG') !SOLWEIG
       ! todo: inconsistent data structure
       dataOutX=dataOutSOL(1:irMax,1:SIZE(varlist),:)

    CASE ('BL') !BL
       dataOutX=dataOutBL(1:irMax,1:SIZE(varlist),:)

    CASE ('snow')    !snow
       dataOutX=dataOutSnow(1:irMax,1:SIZE(varlist),:)

    CASE ('ESTM')    !ESTM
       dataOutX=dataOutESTM(1:irMax,1:SIZE(varlist),:)

    END SELECT

    ! aggregation:
    DO iGrid = 1, NumberOfGrids
       CALL SUEWS_Output_Agg(dataOutX_agg0,dataOutX(:,:,iGrid),varlist,irMax,outFreq_s)
       IF (.NOT. ALLOCATED(dataOutX_agg)) THEN
          ALLOCATE(dataOutX_agg(SIZE(dataOutX_agg0, dim=1),SIZE(varlist),NumberOfGrids), stat=err)
          IF ( err/= 0) PRINT *, ": Allocation request denied"
       ENDIF
       dataOutX_agg(:,:,iGrid)=dataOutX_agg0
    END DO


    ! write out data
    CALL SUEWS_Write_nc(dataOutX_agg,varlist,outLevel)
  END SUBROUTINE SUEWS_Output_nc_grp


  SUBROUTINE SUEWS_Write_nc(dataOut,varlist,outLevel)
    ! generic subroutine to write out data in netCDF format
    USE netCDF

    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:,:,:),INTENT(in)::dataOut
    TYPE(varAttr),DIMENSION(:),INTENT(in)::varlist
    INTEGER,INTENT(in) :: outLevel

    CHARACTER(len=100):: fileOut
    REAL(KIND(1d0)),DIMENSION(:,:,:),ALLOCATABLE::dataOutSel
    TYPE(varAttr),DIMENSION(:),ALLOCATABLE::varListSel

    ! We are writing 3D data, {time, y, x}
    INTEGER, PARAMETER :: NDIMS = 3, iVarStart=6
    INTEGER :: NX,NY,nTime,nVar,err

    ! When we create netCDF files, variables and dimensions, we get back
    ! an ID for each one.
    INTEGER :: ncID, varID, dimids(NDIMS),varIDGrid
    INTEGER :: x_dimid,y_dimid,time_dimid,iVar,varIDx,varIDy,varIDt
    REAL(KIND(1d0)), ALLOCATABLE :: varOut(:,:,:),&
         varX(:,:),varY(:,:),&
         xLat(:,:),xLon(:,:),&
         varSeq0(:),varSeq(:),xTime(:)

    INTEGER :: idVar(iVarStart:SIZE(varlist))
    CHARACTER(len=50):: header_str,longNm_str,unit_str
    CHARACTER(len = 4)  :: yrStr2
    CHARACTER(len = 40) :: startStr2

    ! determine number of times
    nTime=SIZE(dataOut, dim=1)

    !select variables to output
    nVar=COUNT((varList%level<= outLevel), dim=1)
    ALLOCATE(varlistSel(nVar), stat=err)
    IF ( err/= 0) PRINT *, "varlistSel: Allocation request denied"
    varlistSel=PACK(varList, mask=(varList%level<= outLevel))

    ! copy data accordingly
    ALLOCATE(dataOutSel(nTime,nVar,NumberOfGrids), stat=err)
    IF ( err/= 0) PRINT *, "dataOutSel: Allocation request denied"
    ! print*, SIZE(varList%level),PACK((/(i,i=1,SIZE(varList%level))/), varList%level <= outLevel)
    ! print*, nTime,shape(dataOut)
    dataOutSel=dataOut(:,PACK((/(i,i=1,SIZE(varList))/), varList%level <= outLevel),:)

    ! determine filename
    CALL filename_gen(dataOutSel(:,:,1),varlistSel,1,FileOut)

    ! set year string
    WRITE(yrStr2,'(i4)') INT(dataOut(1,1,1))
    ! get start for later time unit creation
    startStr2=TRIM(yrStr2)//'-01-01 00:00:00'

    ! define the dimension of spatial array/frame in the output
    nX = nCol
    nY = nRow

    ALLOCATE(varSeq0(nX*nY))
    ALLOCATE(varSeq(nX*nY))
    ALLOCATE(xLon(nX,nY))
    ALLOCATE(xLat(nX,nY))
    ALLOCATE(varY(nX,nY))
    ALLOCATE(varX(nX,nY))

    ! latitude:
    varSeq0=SiteSelect(1:nX*nY,5)
    CALL sortSeqReal(varSeq0,varSeq,nY,nX)
    xLat = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )
    ! PRINT*, 'before flipping:',xLat(1:2,1)
    xLat =xLat(:,nY:1:-1)
    ! PRINT*, 'after flipping:',xLat(1:2,1)

    ! longitude:
    varSeq0=SiteSelect(1:nX*nY,6)
    CALL sortSeqReal(varSeq0,varSeq,nY,nX)
    xLon = RESHAPE(varSeq,(/nX,nY/),order = (/1,2/) )


    ! pass values to coordinate variables
    varY = xLat
    varX = xLon
    ! PRINT*, 'size x dim 1:',SIZE(varX, dim=1)
    ! PRINT*, 'size x dim 2:',SIZE(varX, dim=2)


    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    PRINT*, 'writing file:',TRIM(fileOut)
    CALL check( nf90_create(TRIM(fileOut), NF90_CLOBBER, ncID) )

    ! Define the dimensions. NetCDF will hand back an ID for each.
    ! nY = ncolumnsDataOut-4
    ! nx = NumberOfGrids
    CALL check( nf90_def_dim(ncID, "time", NF90_UNLIMITED, time_dimid) )
    CALL check( nf90_def_dim(ncID, "west_east", NX, x_dimid) )
    CALL check( nf90_def_dim(ncID, "south_north", NY, y_dimid) )
    ! PRINT*, 'good define dim'

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/x_dimid, y_dimid, time_dimid/)

    ! write out each variable
    ALLOCATE(varOut(nX,nY,nTime))

    ! define all variables
    ! define time variable:
    CALL check( nf90_def_var(ncID,'time', NF90_REAL, time_dimid, varIDt))
    CALL check( nf90_put_att(ncID,varIDt,'units','minutes since '//startStr2 ) )

    ! define coordinate variables:
    CALL check( nf90_def_var(ncID,'xLon', NF90_REAL, (/x_dimid, y_dimid/), varIDx))
    CALL check( nf90_put_att(ncID,varIDx,'units','degree_east') )

    CALL check( nf90_def_var(ncID,'xLat', NF90_REAL, (/x_dimid, y_dimid/), varIDy))
    CALL check( nf90_put_att(ncID,varIDy,'units','degree_north') )

    ! PRINT*, 'good define var'

    ! define grid_ID:
    CALL check( nf90_def_var(ncID,'grid_ID', NF90_INT, (/x_dimid, y_dimid/), varID))
    CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
    varIDGrid=varID

    ! define other 3D variables:
    DO iVar = iVarStart, nVar
       ! define variable name
       header_str = varListSel(iVar)%header
       unit_str   = varListSel(iVar)%unit
       longNm_str = varListSel(iVar)%longNm
       !  PRINT*, unit_str

       ! Define the variable. The type of the variable in this case is
       ! NF90_REAL.
       !  PRINT*, TRIM(ADJUSTL(header_str))
       CALL check( nf90_def_var(ncID,TRIM(ADJUSTL(header_str)), NF90_REAL, dimids, varID) )
       !  PRINT*, 'define good'
       CALL check( nf90_put_att(ncID,varID,'coordinates','xLon xLat') )
       !  PRINT*, 'put coordinates good'
       CALL check( nf90_put_att(ncID,varID,'units',TRIM(ADJUSTL(unit_str))) )
       !  PRINT*, 'put unit good'
       CALL check( nf90_put_att(ncID,varID,'longname',TRIM(ADJUSTL(longNm_str))) )
       !  PRINT*, 'put longname good'
       idVar(iVar)=varID
    END DO
    CALL check( nf90_enddef(ncID) )
    ! End define mode. This tells netCDF we are done defining metadata.

    ! put all variable values into netCDF datasets
    ! put time variable in minute:
    xTime=(dataOutSel(1:nTime,2,1)-1)*24*60+dataOutSel(1:nTime,3,1)*60+dataOutSel(1:nTime,4,1)
    CALL check( nf90_put_var(ncID, varIDt, xTime) )

    ! put coordinate variables:
    CALL check( nf90_put_var(ncID, varIDx, varX) )
    CALL check( nf90_put_var(ncID, varIDy, varY) )
    CALL check( NF90_SYNC(ncID) )
    ! PRINT*, 'good put var'


    ! put grid_ID:
    CALL check( nf90_put_var(ncID, varIDGrid, RESHAPE(GridIDmatrix,(/nX,nY/),order = (/1,2/))) )
    ! PRINT*, 'good put varIDGrid',varIDGrid

    CALL check( NF90_SYNC(ncID) )

    ! then other 3D variables
    DO iVar = iVarStart, nVar
       !  PRINT*, 'dim1:', SIZE(dataOut(1:nTime,iVar,:), dim=1)
       !  PRINT*, 'dim2:',SIZE(dataOut(1:nTime,iVar,:), dim=2)
       ! reshape dataOut to be aligned in checker board form
       varOut = RESHAPE(dataOutSel(1:nTime,iVar,:),(/nX,nY,nTime/),order = (/3,1,2/) )
       varOut = varOut(:,nY:1:-1,:)
       !  get the variable id
       varID= idVar(iVar)
       !  PRINT*, 'good put iVar',iVar
       CALL check( nf90_put_var(ncID, varID, varOut) )
       !  PRINT*, 'good put var',varID
       CALL check(NF90_SYNC(ncID))
    END DO

    IF (ALLOCATED(varOut)) DEALLOCATE(varOut)
    IF (ALLOCATED(varSeq0)) DEALLOCATE(varSeq0)
    IF (ALLOCATED(varSeq)) DEALLOCATE(varSeq)
    IF (ALLOCATED(xLon)) DEALLOCATE(xLon)
    IF (ALLOCATED(xLat)) DEALLOCATE(xLat)
    IF (ALLOCATED(varY)) DEALLOCATE(varY)
    IF (ALLOCATED(varX)) DEALLOCATE(varX)

    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    CALL check( nf90_close(ncID) )

    ! PRINT*, "*** SUCCESS writing netCDF file:"
    ! PRINT*, FileOut
  END SUBROUTINE SUEWS_Write_nc

  !===========================================================================!
  ! convert a vector of grids to a matrix
  ! the grid IDs in seqGrid2Sort follow the QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector
  !===========================================================================!
  SUBROUTINE grid2mat(seqGrid2Sort, seqGridSorted, matGrid, nRow, nCol)

    IMPLICIT NONE

    INTEGER,DIMENSION(nRow*nCol) :: seqGrid2Sort,seqGridSorted
    INTEGER,DIMENSION(nRow,nCol) :: matGrid
    INTEGER :: nRow, nCol,i,j,loc

    CALL sortGrid(seqGrid2Sort, seqGridSorted, nRow, nCol)
    PRINT*, 'old:'
    PRINT*, seqGrid2Sort(1:5)
    PRINT*, 'sorted:'
    PRINT*, seqGridSorted(1:5)
    PRINT*, ''
    DO i = 1, nRow
       DO j = 1, nCol
          loc=(i-1)*nCol+j
          ! PRINT*, i,j,loc
          ! PRINT*, seqGridSorted(loc)
          matGrid(i,j)=seqGridSorted(loc)
       END DO
    END DO
  END SUBROUTINE grid2mat

  !===========================================================================!
  ! convert sequence of REAL values to a matrix
  ! the grid IDs in seqGrid2Sort follow the QGIS convention
  ! the spatial matrix arranges successive rows down the page (i.e., north to south)
  !   and succesive columns across (i.e., west to east)
  ! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector
  !===========================================================================!
  SUBROUTINE seq2mat(seq2Sort, seqSorted, matGrid, nRow, nCol)

    IMPLICIT NONE

    REAL(KIND(1d0)),DIMENSION(nRow*nCol) :: seq2Sort,seqSorted
    REAL(KIND(1d0)),DIMENSION(nRow,nCol) :: matGrid
    INTEGER :: nRow, nCol,i,j,loc

    CALL sortSeqReal(seq2Sort, seqSorted, nRow, nCol)
    PRINT*, 'old:'
    PRINT*, seq2Sort(1:5)
    PRINT*, 'sorted:'
    PRINT*, seqSorted(1:5)
    PRINT*, ''
    DO i = 1, nRow
       DO j = 1, nCol
          loc=(i-1)*nCol+j
          ! PRINT*, i,j,loc
          ! PRINT*, seqGridSorted(loc)
          matGrid(i,j)=seqSorted(loc)
       END DO
    END DO
  END SUBROUTINE seq2mat

  !===========================================================================!
  ! sort a sequence of LONG values into the specially aligned sequence per QGIS
  !===========================================================================!
  SUBROUTINE sortGrid(seqGrid2Sort, seqGridSorted, nRow, nCol)
    USE qsort_c_module
    ! convert a vector of grids to a matrix
    ! the grid IDs in seqGrid2Sort follow the QGIS convention
    ! the spatial matrix arranges successive rows down the page (i.e., north to south)
    !   and succesive columns across (i.e., west to east)
    ! seqGridSorted stores the grid IDs as aligned in matGrid but squeezed into a vector

    IMPLICIT NONE
    INTEGER :: nRow, nCol,i=1,j=1,xInd,len

    INTEGER,DIMENSION(nRow*nCol),INTENT(in) :: seqGrid2Sort
    INTEGER,DIMENSION(nRow*nCol),INTENT(out) :: seqGridSorted
    INTEGER,DIMENSION(nRow*nCol) :: locSorted
    INTEGER :: loc
    REAL:: ind(nRow*nCol,2)
    REAL :: seqGridSortedReal(nRow*nCol),val

    ! number of grids
    len=nRow*nCol
    ! PRINT*, 'after input:'
    ! PRINT*, 'seqGrid2Sort:'
    ! PRINT*, seqGrid2Sort(1:5)
    ! PRINT*, '****'
    !
    !
    ! PRINT*, 'after input:'
    ! PRINT*, 'seqGridSorted:'
    ! PRINT*, seqGridSorted(1:5)
    ! PRINT*, '****'


    ! fill in an nRow*nCol array with values to determine sequence
    xInd=1
    DO i = 1, nRow
       DO j = 1, nCol
          !  {row, col, value for sorting, index in new sequence}
          ind(xInd,:)=(/i+j+i/(nRow+1.),xInd*1./)
          xInd=xInd+1
       END DO
    END DO
    ! PRINT*, 'old after sorting:'
    ! PRINT*, seqGridSorted(1:5)
    ! PRINT*, 'ind:'
    ! PRINT*, ind(:,1)
    ! PRINT*, 'ind. seq:'
    ! PRINT*, ind(:,2)

    ! then sorted ind(:,3) will have the same order as seqGrid2Sort
    ! sort ind(:,3)
    seqGridSortedReal=ind(:,1)*1.
    CALL QsortC(seqGridSortedReal)
    ! print*, 'sorted real:'
    ! print*, seqGridSortedReal

    ! get index of each element of old sequence in the sorted sequence
    DO i=1,len
       ! value in old sequence
       !  val=ind(i,3)*1.
       val=seqGridSortedReal(i)
       DO j=1,len
          IF ( val .EQ. ind(j,1)*1.) THEN
             ! location in sorted sequence
             locSorted(i)=j
          END IF
       END DO
    END DO

    ! put elements of old sequence in the sorted order
    DO i = 1, len
       loc=locSorted(i)
       seqGridSorted(loc)=seqGrid2Sort(i)
    END DO
    seqGridSorted=seqGridSorted(len:1:-1)
    ! PRINT*, 'loc sorted:'
    ! PRINT*, locSorted

    ! PRINT*, 'sorted:'
    ! PRINT*, seqGridSorted(1:5)
    ! PRINT*, 'sort subroutine end!'
  END SUBROUTINE sortGrid

  !===========================================================================!
  ! sort a sequence of REAL values into the specially aligned sequence per QGIS
  !===========================================================================!
  SUBROUTINE sortSeqReal(seqReal2Sort, seqRealSorted, nRow, nCol)
    USE qsort_c_module
    ! convert a vector of grids to a matrix
    ! the grid IDs in seqReal2Sort follow the QGIS convention
    ! the spatial matrix arranges successive rows down the page (i.e., north to south)
    !   and succesive columns across (i.e., west to east)
    ! seqRealSorted stores the grid IDs as aligned in matGrid but squeezed into a vector

    IMPLICIT NONE
    INTEGER :: nRow, nCol,i=1,j=1,xInd,len

    REAL(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(in) :: seqReal2Sort
    REAL(KIND(1d0)),DIMENSION(nRow*nCol),INTENT(out) :: seqRealSorted
    INTEGER(KIND(1d0)),DIMENSION(nRow*nCol) :: locSorted
    INTEGER(KIND(1d0)) :: loc
    REAL:: ind(nRow*nCol,2)
    REAL :: seqRealSortedReal(nRow*nCol),val

    ! number of grids
    len=nRow*nCol
    ! PRINT*, 'after input:'
    ! PRINT*, 'seqReal2Sort:'
    ! PRINT*, seqReal2Sort(1:5)
    ! PRINT*, '****'
    !
    !
    ! PRINT*, 'after input:'
    ! PRINT*, 'seqRealSorted:'
    ! PRINT*, seqRealSorted(1:5)
    ! PRINT*, '****'


    ! fill in an nRow*nCol array with values to determine sequence
    xInd=1
    DO i = 1, nRow
       DO j = 1, nCol
          !  {row, col, value for sorting, index in new sequence}
          ind(xInd,:)=(/i+j+i/(nRow+1.),xInd*1./)
          xInd=xInd+1
       END DO
    END DO
    ! PRINT*, 'old after sorting:'
    ! PRINT*, seqRealSorted(1:5)
    ! PRINT*, 'ind:'
    ! PRINT*, ind(:,1)
    ! PRINT*, 'ind. seq:'
    ! PRINT*, ind(:,2)

    ! then sorted ind(:,3) will have the same order as seqReal2Sort
    ! sort ind(:,3)
    seqRealSortedReal=ind(:,1)*1.
    CALL QsortC(seqRealSortedReal)
    ! print*, 'sorted real:'
    ! print*, seqRealSortedReal

    ! get index of each element of old sequence in the sorted sequence
    DO i=1,len
       ! value in old sequence
       !  val=ind(i,3)*1.
       val=seqRealSortedReal(i)
       DO j=1,len
          IF ( val .EQ. ind(j,1)*1.) THEN
             ! location in sorted sequence
             locSorted(i)=j
          END IF
       END DO
    END DO

    ! put elements of old sequence in the sorted order
    DO i = 1, len
       loc=locSorted(i)
       seqRealSorted(loc)=seqReal2Sort(i)
    END DO
    seqRealSorted=seqRealSorted(len:1:-1)
    ! PRINT*, 'loc sorted:'
    ! PRINT*, locSorted

    ! PRINT*, 'sorted:'
    ! PRINT*, seqRealSorted(1:5)
    ! PRINT*, 'sort subroutine end!'
  END SUBROUTINE sortSeqReal

  !===========================================================================!
  ! a wrapper for checking netCDF status
  !===========================================================================!

  SUBROUTINE check(status)
    USE netcdf
    IMPLICIT NONE

    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT *, TRIM(nf90_strerror(status))
       STOP "Stopped"
    END IF
  END SUBROUTINE check
#endif

END MODULE ctrl_output
