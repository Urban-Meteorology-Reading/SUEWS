from SUEWS_driver import suews_driver as sd
import numpy as np
# print sd.__doc__
# print sd.suews_cal_main.__doc__
dir(sd)
len(dir(sd))

xp = 1233455
sd.set_nan([99999, 12])
# print sd.set_nan.__doc__
# print sd.square.__doc__
sd.square(1.9)


# print sd.suews_cal_qn.__doc__

ncolumnsdataout = 10  # : input int, optional #: shape(dataout,1)
nsh = 10  # : input int, optional#: (shape(ahprof_tstep,0))/(24)
numberofgrids = 10  # : input int, optional#: shape(dataout,2)
readlinesmetdata = 10  # : input int, optional#: shape(dataout,0)

a1 = np.float(1.3)  # : in/output rank-0 array(float,'d')
a2 = np.float(1.3)  # : in/output rank-0 array(float,'d')
a3 = np.float(1.3)  # : in/output rank-0 array(float,'d')
addimpervious = 1.3  # : input float
addpipes = 1.3  # : input float
addveg = 1.3  # : input float
addwater = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
addwaterbody = 1.3  # : input float
# : in/output rank-1 array('d') with bounds (7)
addwaterrunoff = 1.3 * np.ones(7)
aerodynamicresistancemethod = 1  # : input int
ah_min = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
# : input rank-2 array('d') with bounds (24 * nsh,2)
ahprof_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
# : input rank-1 array('d') with bounds (2)
ah_slope_cooling = 1.3 * np.ones(2)
# : input rank-1 array('d') with bounds (2)
ah_slope_heating = 1.3 * np.ones(2)
alb = .3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
# : in/output rank-1 array('d') with bounds (367)
albdectr = 1.3 * np.ones(367)
albedochoice = 1  # : input int
# : in/output rank-1 array('d') with bounds (367)
albevetr = 1.3 * np.ones(367)
# : in/output rank-1 array('d') with bounds (367)
albgrass = 1.3 * np.ones(367)
albmax_dectr = 1.3  # : input float
albmax_evetr = 1.3  # : input float
albmax_grass = 1.3  # : input float
albmin_dectr = 1.3  # : input float
albmin_evetr = 1.3  # : input float
albmin_grass = 1.3  # : input float
alpha_bioco2 = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
# : input rank-1 array('d') with bounds (3)
alpha_enh_bioco2 = 1.3 * np.ones(3)
alt = 1.3  # : input float
areazh = 1.3  # : input float
avdens = 1.3  # : input float
avkdn = 1000  # : input float
avrh = 1.3  # : input float
avu1 = 1.3  # : input float
baset = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
basete = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
basethdd = 1.3  # : input float
beta_bioco2 = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
beta_enh_bioco2 = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
biogenco2code = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
bldgh = 1.3  # : input float
capmax_dec = 1.3  # : input float
capmin_dec = 1.3  # : input float
chanohm = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
cpanohm = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
crwmax = 1.3  # : input float
crwmin = 1.3  # : input float
# : input rank-3 array('d') with bounds (readlinesmetdata,ncolumnsdataout,numberofgrids)
dataout = 1.3 * \
    np.ones((readlinesmetdata, ncolumnsdataout, numberofgrids), order='F')
# : in/output rank-3 array('d') with bounds (readlinesmetdata,32,numberofgrids)
dataoutestm = 1.3 * np.ones((readlinesmetdata, 32, numberofgrids), order='F')
# : in/output rank-2 array('i') with bounds (367,3)
dayofweek = 1 * np.ones((367, 3), order='F',dtype=np.int32)
daywat = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
daywatper = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
# : in/output rank-1 array('d') with bounds (367)
decidcap = 1.3 * np.ones(367)
dectime = 1.3  # : input float
dectreeh = 1.3  # : input float
diagnose = 1  # : input int
diagqn = 1  # : input int
diagqs = 1  # : input int
dls = 1  # : input int
drainrt = 1.3  # : input float
ef_umolco2perj = 1.3  # : input float
emis = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
emissionsmethod = 1  # : input int
enef_v_jkm = 1.3  # : input float
evetreeh = 1.3  # : input float
faibldg = 1.3  # : input float
faidectree = 1.3  # : input float
faievetree = 1.3  # : input float
faut = 1.3  # : input float
fcef_v_kgkm = 1.3  # : input float
fcld_obs = 1.3  # : input float
frfossilfuel_heat = 1.3  # : input float
frfossilfuel_nonheat = 1.3  # : input float
g1 = 1.3  # : input float
g2 = 1.3  # : input float
g3 = 1.3  # : input float
g4 = 1.3  # : input float
g5 = 1.3  # : input float
g6 = 1.3  # : input float
# : in/output rank-2 array('d') with bounds (367,5)
gdd = 1.3 * np.ones((367, 5), order='F')
gddfull = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
gridiv = 1  # : input int
gsmodel = 1  # : input int
halftimestep = 1.3  # : input float
# : in/output rank-2 array('d') with bounds (371,6)
hdd = 1.3 * np.ones((371, 6), order='F')
# : input rank-2 array('d') with bounds (24 * nsh,2)
humactivity_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
icefrac = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
id = 1  # : input int
ie_a = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
ie_end = 3  # : input int
ie_m = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
ie_start = 1  # : input int
imin = 1  # : input int
internalwateruse_h = 1.3  # : input float
ir = 1  # : input int
irrfracconif = 1.3  # : input float
irrfracdecid = 1.3  # : input float
irrfracgrass = 1.3  # : input float
it = 1  # : input int
ity = 1  # : input int
iy = 1  # : input int
k = 1.3  # : input float
kkanohm = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
kmax = 1.3  # : input float
# : in/output rank-2 array('d') with bounds (371,3)
lai = 1.3 * np.ones((371, 3), order='F')
laicalcyes = 3  # : input int
laimax = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
laimin = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
lai_obs = 1.3  # : input float
# : input rank-2 array('d') with bounds (4,3)
laipower = 1.3 * np.ones((4, 3), order='F')
laitype = 1 * np.ones(3)  # : input rank-1 array('i') with bounds (3)
lat = 1.3  # : input float
ldown_obs = 1.3  # : input float
ldown_option = 1  # : input int
lng = 1.3  # : input float
maxconductance = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
maxqfmetab = 1.3  # : input float
# : in/output rank-1 array('d') with bounds (7)
meltwaterstore = 1.3 * np.ones(7)
# : input rank-3 array('d') with bounds (f2py_metforcingdata_d0,f2py_metforcingdata_d1,f2py_metforcingdata_d2)
metforcingdata = 1.3 * np.ones((30, 30, 30), order='F')
minqfmetab = 1.3  # : input float
min_res_bioco2 = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
narp_emis_snow = 1.3  # : input float
narp_g = 1.3 * np.ones(365)  # : input rank-1 array('d') with bounds (365)
narp_trans_site = 1.3  # : input float
netradiationmethod = 1  # : input int
nonwaterfraction = 1.3  # : input float
nsh_real = 1.3  # : input float
numcapita = 1.3  # : input float
# : input rank-3 array('d') with bounds (9,4,3)
ohm_coef = 1.3 * np.ones((9, 4, 3), order='F')
ohmincqf = 0  # : input int
ohm_threshsw = 1.3 * np.ones(9)  # : input rank-1 array('d') with bounds (9)
ohm_threshwd = 1.3 * np.ones(9)  # : input rank-1 array('d') with bounds (9)
overuse = 1.3  # : in/output rank-0 array(float,'d')
pervfraction = 1.3  # : input float
pipecapacity = 1.3  # : input float
popdensdaytime = 1.3  # : input float
popdensnighttime = 1.3  # : input float
# : input rank-2 array('d') with bounds (24 * nsh,2)
popprof_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
pormax_dec = 1.3  # : input float
pormin_dec = 1.3  # : input float
# : in/output rank-1 array('d') with bounds (367)
porosity = 1.3 * np.ones(367)
precip = 1.3  # : input float
preciplimit = 1.3  # : input float
preciplimitalb = 1.3  # : input float
press_hpa = np.float64(921.3)  # : input float
qf0_beu = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
qf_a = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
qf_b = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
qf_c = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
qh_obs = 1.3  # : input float
# : in/output rank-1 array('d') with bounds (2 * nsh + 1)
qn1_av_store = 1.3 * np.ones((2 * nsh + 1), order='F')
qn1_obs = 1.3  # : input float
# : in/output rank-1 array('d') with bounds (2 * nsh + 1)
qn1_s_av_store = 1.3 * np.ones((2 * nsh + 1), order='F')
# : in/output rank-1 array('d') with bounds (nsh)
qn1_s_store = 1.3 * np.ones(nsh)
# : in/output rank-1 array('d') with bounds (nsh)
qn1_store = 1.3 * np.ones(nsh)
radmeltfact = 1.3  # : input float
raincover = 1.3  # : input float
rainmaxres = 1.3  # : input float
resp_a = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
resp_b = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
roughlenheatmethod = 1  # : input int
roughlenmommethod = 1  # : input int
runoff_per_interval = 1.3  # : in/output rank-0 array(float,'d')
runoffsoil = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
runofftowater = 1.3  # : input float
s1 = 1.3  # : input float
s2 = 1.3  # : input float
# : input rank-1 array('d') with bounds (7)
sathydraulicconduct = 1.3 * np.ones(7)
sddfull = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
sfr = 1.3 * np.ones(7)   # : input rank-1 array('d') with bounds (7)
smdmethod = 1  # : input int
snowalb = np.float(1.3)  # : in/output rank-0 array(float,'d')
snowalbmax = 1.3  # : input float
snowalbmin = 1.3  # : input float
snowdens = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
snowdensmax = 1.3  # : input float
snowdensmin = 1.3  # : input float
snowdepth = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
snowfrac = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
snowfractionchoice = 1  # : input int
snowlimbuild = 1.3  # : input float
snowlimpaved = 1.3  # : input float
snow_obs = 1.3  # : input float
snowpack = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
snowuse = 1  # : input int
soildepth = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
soilmoist = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
soilstorecap = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
stabilitymethod = 1  # : input int
state = 1.3 * np.ones(7)  # : in/output rank-1 array('d') with bounds (7)
statelimit = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
storageheatmethod = 1  # : input int
# : in/output rank-2 array('d') with bounds (6,7)
surf = 1.3 * np.ones((6, 7), order='F')
surfacearea = 1.3  # : input float
surplusevap = 1.3 * np.ones(2)  # : in/output rank-1 array('d') with bounds (2)
# : in/output rank-1 array('d') with bounds (24 * nsh)
tair24hr = 1.3 * np.ones(24 * nsh)
tau_a = 1.3  # : input float
tau_f = 1.3  # : input float
tau_r = 1.3  # : input float
# : input rank-1 array('d') with bounds (2)
t_critic_cooling = 1.3 * np.ones(2)
# : input rank-1 array('d') with bounds (2)
t_critic_heating = 1.3 * np.ones(2)
temp_c = 1.3  # : input float
tempmeltfact = 1.3  # : input float
th = 1.3  # : input float
theta_bioco2 = 1.3 * np.ones(3)  # : input rank-1 array('d') with bounds (3)
timezone = 1.3  # : input float
tl = 1.3  # : input float
trafficrate = 1.3 * np.ones(2)  # : input rank-1 array('d') with bounds (2)
trafficunits = 1.3  # : input float
# : input rank-2 array('d') with bounds (24 * nsh,2)
traffprof_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
# : input rank-1 array('d') with bounds (f2py_ts5mindata_ir_d0)
ts5mindata_ir = 1.3 * np.ones(20)
tstep = 5  # : input int
tstepcount = np.float(1.3)  # : in/output rank-0 array(float,'d')
tstep_real = 1.3  # : input float
tsurf_ind = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
vegfraction = 1.3  # : input float
veg_type = 1  # : input int
waterdens = 1.3  # : input float
# : input rank-2 array('d') with bounds (8,6)
waterdist = 1.3 * np.ones((8, 6), order='F')
waterusemethod = 1  # : input int
wetthresh = 1.3 * np.ones(7)  # : input rank-1 array('d') with bounds (7)
# : in/output rank-2 array('d') with bounds (367,9)
wu_day = 1.3 * np.ones((367, 9), order='F')
# : input rank-2 array('d') with bounds (24 * nsh,2)
wuprofa_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
# : input rank-2 array('d') with bounds (24 * nsh,2)
wuprofm_tstep = 1.3 * np.ones((24 * nsh, 2), order='F')
xsmd = 1.3  # : input float
year = 1.3  # : input float
z = 1.3  # : input float


# additionalwater,\
# avu10_ms,\
# azimuth,\
# chang,\
# changsnow,\
# chsnow_per_interval,\
# cumsnowfall,\
# dens_dry,\
# drain_per_tstep,\
# ea_hpa,\
# e_mod,\
# es_hpa,\
# ev,\
# evap,\
# ev_per_tstep,\
# ev_snow,\
# ext_wu,\
# fc,\
# fc_anthro,\
# fc_biogen,\
# fc_build,\
# fcld,\
# fc_metab,\
# fc_photo,\
# fc_respi,\
# fc_traff,\
# flowchange,\
# freezmelt,\
# fwh,\
# gsc,\
# h_mod,\
# int_wu,\
# kclear,\
# kup,\
# kup_ind_snow,\
# ldown,\
# l_mod,\
# lup,\
# mwh,\
# mw_ind,\
# mwstore,\
# nwstate_per_tstep,\
# planf,\
# p_mm,\
# psim,\
# q2_gkg,\
# qeout,\
# qe_per_tstep,\
# qf,\
# qf_sahp,\
# qh,\
# qh_r,\
# qm,\
# qmfreez,\
# qm_freezstate,\
# qm_melt,\
# qm_rain,\
# qmrain,\
# qn1,\
# qn1_ind_snow,\
# qn1_s,\
# qn1_sf,\
# qs,\
# ra,\
# rainonsnow,\
# resistsurf,\
# rss,\
# rss_nsurf,\
# runoff,\
# runoffagimpervious,\
# runoffagveg,\
# runoff_per_tstep,\
# runoffpipes,\
# runoffpipes_m3,\
# runoffsnow,\
# runoffsoil_per_tstep,\
# runoffwaterbody,\
# runoffwaterbody_m3,\
# smd,\
# smd_nsurf,\
# snowd,\
# snowprof,\
# snowremoval,\
# snowtosurf,\
# soilstate,\
# state_per_tstep,\
# surf_chang_per_tstep,\
# swe,\
# t2_c,\
# tempveg,\
# tot_chang_per_tstep,\
# tstar,\
# tsurf,\
# tsurf_ind_snow,\
# ustar,\
# vpd_pa,\
# wuareadectr_m2,\
# wuareaevetr_m2,\
# wuareagrass_m2,\
# wuareatotal_m2,\
# wu_dectr,\
# wu_evetr,\
# wu_grass,\
# wu_m3,\
# xbo,\
# z0m,\
# zdm,\
# zenith_deg,\
# zh
result= sd.suews_cal_main(
    a1,
    a2,
    a3,
    addimpervious,
    addpipes,
    addveg,
    addwater,
    addwaterbody,
    addwaterrunoff,
    aerodynamicresistancemethod,
    ah_min,
    ahprof_tstep,
    ah_slope_cooling,
    ah_slope_heating,
    alb,
    albdectr,
    albedochoice,
    albevetr,
    albgrass,
    albmax_dectr,
    albmax_evetr,
    albmax_grass,
    albmin_dectr,
    albmin_evetr,
    albmin_grass,
    alpha_bioco2,
    alpha_enh_bioco2,
    alt,
    areazh,
    avdens,
    avkdn,
    avrh,
    avu1,
    baset,
    basete,
    basethdd,
    beta_bioco2,
    beta_enh_bioco2,
    biogenco2code,
    bldgh,
    capmax_dec,
    capmin_dec,
    chanohm,
    cpanohm,
    crwmax,
    crwmin,
    dataout,
    dataoutestm,
    dayofweek,
    daywat,
    daywatper,
    decidcap,
    dectime,
    dectreeh,
    diagnose,
    diagqn,
    diagqs,
    dls,
    drainrt,
    ef_umolco2perj,
    emis,
    emissionsmethod,
    enef_v_jkm,
    evetreeh,
    faibldg,
    faidectree,
    faievetree,
    faut,
    fcef_v_kgkm,
    fcld_obs,
    frfossilfuel_heat,
    frfossilfuel_nonheat,
    g1,
    g2,
    g3,
    g4,
    g5,
    g6,
    gdd,
    gddfull,
    gridiv,
    gsmodel,
    halftimestep,
    hdd,
    humactivity_tstep,
    icefrac,
    id,
    ie_a,
    ie_end,
    ie_m,
    ie_start,
    imin,
    internalwateruse_h,
    ir,
    irrfracconif,
    irrfracdecid,
    irrfracgrass,
    it,
    ity,
    iy,
    k,
    kkanohm,
    kmax,
    lai,
    laicalcyes,
    laimax,
    laimin,
    lai_obs,
    laipower,
    laitype,
    lat,
    ldown_obs,
    ldown_option,
    lng,
    maxconductance,
    maxqfmetab,
    meltwaterstore,
    metforcingdata,
    minqfmetab,
    min_res_bioco2,
    narp_emis_snow,
    narp_g,
    narp_trans_site,
    netradiationmethod,
    nonwaterfraction,
    nsh_real,
    numcapita,
    ohm_coef,
    ohmincqf,
    ohm_threshsw,
    ohm_threshwd,
    overuse,
    pervfraction,
    pipecapacity,
    popdensdaytime,
    popdensnighttime,
    popprof_tstep,
    pormax_dec,
    pormin_dec,
    porosity,
    precip,
    preciplimit,
    preciplimitalb,
    press_hpa,
    qf0_beu,
    qf_a,
    qf_b,
    qf_c,
    qh_obs,
    qn1_av_store,
    qn1_obs,
    qn1_s_av_store,
    qn1_s_store,
    qn1_store,
    radmeltfact,
    raincover,
    rainmaxres,
    resp_a,
    resp_b,
    roughlenheatmethod,
    roughlenmommethod,
    runoff_per_interval,
    runoffsoil,
    runofftowater,
    s1,
    s2,
    sathydraulicconduct,
    sddfull,
    sfr,
    smdmethod,
    snowalb,
    snowalbmax,
    snowalbmin,
    snowdens,
    snowdensmax,
    snowdensmin,
    snowdepth,
    snowfrac,
    snowfractionchoice,
    snowlimbuild,
    snowlimpaved,
    snow_obs,
    snowpack,
    snowuse,
    soildepth,
    soilmoist,
    soilstorecap,
    stabilitymethod,
    state,
    statelimit,
    storageheatmethod,
    surf,
    surfacearea,
    surplusevap,
    tair24hr,
    tau_a,
    tau_f,
    tau_r,
    t_critic_cooling,
    t_critic_heating,
    temp_c,
    tempmeltfact,
    th,
    theta_bioco2,
    timezone,
    tl,
    trafficrate,
    trafficunits,
    traffprof_tstep,
    ts5mindata_ir,
    tstep,
    tstepcount,
    tstep_real,
    tsurf_ind,
    vegfraction,
    veg_type,
    waterdens,
    waterdist,
    waterusemethod,
    wetthresh,
    wu_day,
    wuprofa_tstep,
    wuprofm_tstep,
    xsmd,
    year,
    z,
    [ncolumnsdataout,
     nsh,
     numberofgrids,
     readlinesmetdata])



zenith_deg=1
diagqn=1
snowuse=0
snowfrac,ldown,fcld,qn1,qn1_sf,qn1_s,kclear,kup,lup,tsurf,qn1_ind_snow,kup_ind_snow,tsurf_ind_snow = sd.suews_cal_qn(netradiationmethod,snowuse,ldown_option,id,diagnose,snow_obs,ldown_obs,fcld_obs,dectime,zenith_deg,avkdn,temp_c,avrh,press_hpa,qn1_obs,snowalb,albedochoice,diagqn,narp_g,narp_trans_site,narp_emis_snow,icefrac,sfr,emis,alb,albdectr,decidcap,albevetr,albgrass,surf)
qn1
# print sd.suews_cal_qn.__doc__







#
