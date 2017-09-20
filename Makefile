CC = gfortran $(CFLAGS)          # compiler
CC_nc = gfortran $(CFLAGS_nc)    # compiler with netcdf support
TARGET = SUEWS_V2017b      # program name

# static flag for different OS to correctly link static libs
# so gfortran dependency can be relaxed
# but netCDF is still linked in the dynamic way as suggested by UCAR
ifeq ($(OS),Windows_NT)
  	STATIC = -static # mingw
		STATICLIB =
else
		STATIC = -static-libgfortran -static-libgcc # single -static won't work on macOS
		STATICLIB = libquadmath.a # force-ship this static lib
endif

CFLAGS = $(STATIC) -g -pg -Wall -Wtabs -fbounds-check -cpp -Wno-unused-dummy-argument -Wno-unused-variable \
				-fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal
CFLAGS_nc = $(CFLAGS) \
				-I`nc-config --includedir` -Dnc=1 # options for netcdf build

# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
MODULES = LUMPS_Module_constants.o  \
          LUMPS_metRead.o  \
		  		SUEWS_MetDisagg.o \
          SOLWEIG_modules.o  \
          SUEWS_Files_run_Control.o \
					precmod.o \
					stringmod.o \
					qsort_c_module.o\
					LUMPS_NARP_v3.o \
					SUEWS_AnOHM.o \
					SUEWS_OHM.o \
					SUEWS_driver.o \
					minpack.o

# Rest of the files including modules and functions which are independent
OTHERS =  BLUEWS_CBL.o   \
          BLUEWS_Diff.o  \
          SUEWS_translate.o \
          SUEWS_HorizontalSoilWater.o \
          LUMPS_atmos_functions_moist.o \
          LUMPS_atmos_functions_stab.o \
          SUEWS_ReDistributeWater.o \
          SUEWS_RoughnessParameters.o \
          SUEWS_RunoffFromGrid.o \
          NARP_sun_position_v2.o \
          SUEWS_SAHP.o \
          LUMPS_OutputHeaders.o \
          LUMPS_QHQE.o \
          SUEWS_AerodynamicResistance.o \
          SUEWS_store.o \
          SUEWS_BL_Resist.o \
          SUEWS_SurfaceResistance.o \
          SUEWS_drain.o \
          SUEWS_TimeRelatedSubroutines.o \
          SUEWS_error.o \
          SUEWS_waterUse.o \
          SUEWS_evap.o \
          SUEWS_CodeMatch.o \
          SUEWS_InputHeaders.o \
          SUEWS_InterpHourlyProfiles.o \
          SOLWEIG_2014a_core.o  \
					SOLWEIG_Initial.o  \
          SOLWEIG_clearnessindex_2013b.o  \
          SOLWEIG_cylindric_wedge.o  \
          SOLWEIG_diffusefraction.o  \
          SOLWEIG_EsriAsciiGrid.o  \
          SOLWEIG_Kside_veg_v24.o  \
          SOLWEIG_Lside_veg_v2.o  \
          SOLWEIG_Lvikt_veg.o  \
          SOLWEIG_misc.o  \
          SOLWEIG_shadowingfunction_10.o  \
          SOLWEIG_shadowingfunction_20.o  \
          SOLWEIG_sunonsurface_veg.o  \
          SOLWEIG_wallinsun_veg.o  \
          SUEWS_ESTM_functions.o \
          SUEWS_ESTM_initials.o \
          SUEWS_ESTM_v2016.o \
          SUEWS_CO2.o \
					SUEWS_Initial.o \
					SUEWS_SoilMoist.o \
					SUEWS_Snow.o \
					SUEWS_SnowUpdate.o \
					SUEWS_Diagnostics.o

# modules under rapid development
TEST =		SUEWS_Calculations.o \
					SUEWS_DailyState.o


# Build main program - main uses MODULES and OTHERS
main: SUEWS_Program.f95 $(MODULES) $(OTHERS) $(TEST)
	$(CC) SUEWS_ctrl_output.f95  -c ; \
	$(CC) SUEWS_Program.f95  -c ; \
	$(CC) SUEWS_Program.o $(MODULES) $(OTHERS) $(TEST) \
	$(STATICLIB) \
	SUEWS_ctrl_output.o	-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
netcdf: SUEWS_Program.f95 $(MODULES) $(OTHERS) $(TEST)
	$(CC_nc) SUEWS_ctrl_output.f95  -c ; \
	$(CC_nc) SUEWS_Program.f95  -c ; \
	$(CC_nc) SUEWS_Program.o $(MODULES) $(OTHERS) $(TEST) \
	$(STATICLIB) \
	SUEWS_ctrl_output.o \
	-L`nc-config --libdir` -lnetcdf -lnetcdff \
	-o $(TARGET)

# If OTHERS have changed, compile them again
$(OTHERS): $(MODULES) $(subst .o,.f95, $(OTHERS))
	$(CC) -c $(subst .o,.f95, $@)

# If MODULES have changed, compile them again
$(MODULES): $(subst .o,.f95, $(MODULES))
	$(CC) -c $(subst .o,.f95, $@)

# If TEST have changed, compile them again
$(TEST): $(subst .o,.f95, $(TEST))
	$(CC) -c $(subst .o,.f95, $@)

# If wanted, clean all *.o files after build
clean:
	-rm -rf *.o *.mod *.dSYM $(TARGET)
