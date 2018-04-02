CC = gfortran $(CFLAGS)          # compiler
TARGET = SUEWS_V2017c      # program name

# static flag for different OS to correctly link static libs
# so gfortran dependency can be relaxed
# but netCDF is still linked in the dynamic way as suggested by UCAR
ifeq ($(OS),Windows_NT)
			STATIC = -static -mcmodel=medium# mingw
			STATICLIB =
else
			UNAME_S := $(shell uname -s)
			ifeq ($(UNAME_S),Linux) # Linux
											STATIC = -static-libgfortran -static-libgcc # single -static won't work on macOS
											STATICLIB =
			endif
			ifeq ($(UNAME_S),Darwin) # macOS
											STATIC = -static-libgfortran -static-libgcc # single -static won't work on macOS
											STATICLIB = #libquadmath.a # force-ship this static lib
			endif
endif

CFLAGS = $(STATIC) -g -pg -Wall -Wtabs -fbounds-check -cpp -Wno-unused-dummy-argument -Wno-unused-variable \
									-fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal


# netCDF-related settings:
netcdf: NETCDFINC = `nc-config --includedir` # ordinary path for netCDF directories
netcdf:	NETCDFLIB = `nc-config --libdir` # ordinary path for netCDF directories

nc4fr: NETCDFINC = /home/xlinfr/apps/lib4cdo/include # path only valid with Fredrik's HPC
nc4fr: NETCDFLIB = /home/xlinfr/apps/lib4cdo/lib # path only valid with Fredrik's HPC

netcdf nc4fr:	CFLAGS += -Dnc=1 -I$(NETCDFINC) # options for netcdf build


# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
MODULES = LUMPS_Module_constants.o  \
					LUMPS_metRead.o  \
					LUMPS_NARP.o \
					SUEWS_MetDisagg.o \
					LUMPS_atmos_functions_moist.o \
					SOLWEIG_modules.o  \
					BLUEWS_module.o \
					SUEWS_Files_run_Control.o \
					stringmod.o \
					minpack.o\
					SUEWS_DailyState.o \
					qsort_c_module.o\
					SUEWS_ESTM.o \
					SUEWS_AnOHM.o \
					SUEWS_OHM.o \
					SUEWS_WaterDist.o \
					SUEWS_Snow.o \
					SUEWS_ctrl_output.o \
					SUEWS_driver.o


# Rest of the files including modules and functions which are independent
OTHERS =  SUEWS_translate.o \
					SUEWS_HorizontalSoilWater.o \
					LUMPS_atmos_functions_stab.o \
					SUEWS_RoughnessParameters.o \
					LUMPS_QHQE.o \
					SUEWS_AerodynamicResistance.o \
					SUEWS_BL_Resist.o \
					SUEWS_SurfaceResistance.o \
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
					SUEWS_Initial.o \
					SUEWS_AnthropogenicEmissions.o \
					SUEWS_BiogenCO2.o \
					SUEWS_SnowUpdate.o



# modules under rapid development
TEST =      SUEWS_Calculations.o


# Build main program - main uses MODULES and OTHERS
main: SUEWS_Program.f95 $(MODULES) $(OTHERS) $(TEST)
			$(CC) SUEWS_Program.f95  -c ; \
			$(CC) SUEWS_Program.o $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
netcdf: SUEWS_Program.f95 $(MODULES) $(OTHERS) $(TEST)
				$(CC) SUEWS_Program.f95  -c ; \
				$(CC) SUEWS_Program.o $(MODULES) $(OTHERS) $(TEST) \
				$(STATICLIB) \
				-L$(NETCDFLIB) -lnetcdf -lnetcdff \
				-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
nc4fr: SUEWS_Program.f95 $(MODULES) $(OTHERS) $(TEST)
			 $(CC) SUEWS_Program.f95  -c ; \
			 $(CC) SUEWS_Program.o $(MODULES) $(OTHERS) $(TEST) \
			 $(STATICLIB) \
			 -L$(NETCDFLIB) -Wl,--rpath -Wl,$(NETCDFLIB) -lnetcdff -lnetcdf \
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
