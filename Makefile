	

CC = gfortran $(CFLAGS)		        	# compiler
TARGET = SUEWS_V2015a			# program name
CFLAGS = -g -Wall -Wtabs 
# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
MODULES = LUMPS_Module_constants_v2015a.o  \
          LUMPS_metRead.o  \
          SOLWEIG_modules.o  \
          SUEWS_Files_run_Control.o
# Rest of the files including modules and functions which are independend
OTHERS =  BLUEWS_CBL.o   \
          LUMPS_NARP_v3.o \
          SUEWS_Calculations.o  \
          BLUEWS_Diff.o  \
          SUEWS_translate.o \
          SUEWS_HorizontalSoilWater.o \
          LUMPS_accumulate.o \
          SUEWS_Initial.o \
          LUMPS_atmos_functions_moist.o \
          SUEWS_OHM_v2015.o \
          LUMPS_atmos_functions_stab.o \
          SUEWS_ReDistributeWater.o \
          SUEWS_RoughnessParameters.o \
          SUEWS_RunoffFromGrid.o \
          NARP_sun_position_v2.o \
          SUEWS_SAHP_v2015.o \
          LUMPS_OutputHeaders.o \
          SUEWS_Snow.o \
          LUMPS_QHQE.o \
          SUEWS_AerodynamicResistance.o \
          SUEWS_store.o \
          SUEWS_BL_Resist.o \
          SUEWS_SurfaceResistance.o \
          SUEWS_DailyState.o \
          SUEWS_SnowUpdate.o \
          SUEWS_drain.o \
          SUEWS_TimeRelatedSubroutines.o \
          SUEWS_error.o \
          SUEWS_waterUse.o \
          SUEWS_evap.o \
          SUEWS_Output.o  \
          SUEWS_CodeMatch.o \
          SUEWS_InputHeaders.o \
          SUEWS_InterpHourlyProfiles.o \
          SOLWEIG_2014a_core.o  \
          SOLWEIG_clearnessindex_2013b.o  \
          SOLWEIG_cylindric_wedge.o  \
          SOLWEIG_diffusefraction.o  \
          SOLWEIG_EsriAsciiGrid.o  \
          SOLWEIG_Initial.o  \
          SOLWEIG_Kside_veg_v24.o  \
          SOLWEIG_Lside_veg_v2.o  \
          SOLWEIG_Lvikt_veg.o  \
          SOLWEIG_misc.o  \
          SOLWEIG_shadowingfunction_10.o  \
          SOLWEIG_shadowingfunction_20.o  \
          SOLWEIG_sunonsurface_veg.o  \
          SOLWEIG_wallinsun_veg.o  \

# Build main program - main uses MODULES and OTHERS
main: SUEWS_Program.f95 $(MODULES) $(OTHERS)
	$(CC) SUEWS_Program.f95 $(CFLAGS) -c; \
	$(CC) SUEWS_Program.o $(MODULES) $(OTHERS) -o $(TARGET)

# If OTHERS have changed, compile them again
$(OTHERS): $(MODULES) $(subst .o,.f95, $(OTHERS))
	$(CC) -c $(subst .o,.f95, $@)

# If MODULES have changed, compile them again
$(MODULES): $(subst .o,.f95, $(MODULES))
	$(CC) -c $(subst .o,.f95, $@)

# If wanted, clean all *.o files after build
clean: 
	    -rm *.o