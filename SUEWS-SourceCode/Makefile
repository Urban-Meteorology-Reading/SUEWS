
CC = gfortran $(CFLAGS)          # compiler

file       := SUEWS_const.f95
variable   := $(shell grep 'progname' ${file})
TARGET_str := $(lastword $(subst =, ,${variable})) # program name with quotes
TARGET     = $(subst $\', ,${TARGET_str}) # program name

# static flag for different OS to correctly link static libs
# so gfortran dependency can be relaxed
# but netCDF is still linked in the dynamic way as suggested by UCAR
SUEWS_build_base = ../ReleaseRepo/build/

ifeq ($(OS),Windows_NT)
			STATIC = -static -mcmodel=medium# mingw
			STATICLIB =
			SUEWS_build_dir = $(addprefix $(SUEWS_build_base), Win10x64)
else
			UNAME_S := $(shell uname -s)
			ifeq ($(UNAME_S),Linux) # Linux
				STATIC = -static-libgfortran -static-libgcc # single -static won't work on macOS
				STATICLIB =
				SUEWS_build_dir = $(addprefix $(SUEWS_build_base), Linux)
			endif
			ifeq ($(UNAME_S),Darwin) # macOS
				STATIC = -static-libgfortran -static-libgcc# single -static won't work on macOS
				STATICLIB = #libquadmath.a # force-ship this static lib
				SUEWS_build_dir = $(addprefix $(SUEWS_build_base), macOS)
			endif
endif

CFLAGS = $(STATIC) -g -pg -Wall -Wtabs -fbounds-check -cpp \
					-Wno-unused-dummy-argument -Wno-unused-variable \
					-fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal


# netCDF-related settings:
netcdf: NETCDFINC = `nc-config --includedir` # ordinary path for netCDF directories
netcdf:	NETCDFLIB = `nc-config --libdir` # ordinary path for netCDF directories

nc4fr: NETCDFINC = /home/xlinfr/apps/lib4cdo/include # path only valid with Fredrik's HPC
nc4fr: NETCDFLIB = /home/xlinfr/apps/lib4cdo/lib # path only valid with Fredrik's HPC

netcdf nc4fr:	CFLAGS += -Dnc=1 -I$(NETCDFINC) # options for netcdf build


# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
UTILS = stringmod.o \
				qsort_c_module.o\
				SUEWS_TimeRelatedSubroutines.o \
				SUEWS_const.o  \
				meteo.o \
				minpack.o

MODULES =	NARP_module.o \
					SUEWS_ctrl_input.o  \
					LUMPS_atmos.o \
					SOLWEIG_module.o  \
					BLUEWS_module.o \
					SUEWS_WaterDist.o \
					SUEWS_Snow.o \
					SUEWS_DailyState.o \
					SUEWS_ESTM.o \
					SUEWS_AnOHM.o \
					SUEWS_OHM.o \
					SUEWS_ctrl_output.o \
					SUEWS_driver.o


# Rest of the files including modules and functions which are independent
OTHERS =  SUEWS_translate.o \
					LUMPS_QHQE.o \
					SUEWS_Resist.o \
					SUEWS_error.o \
					SUEWS_evap.o \
					SUEWS_Initial.o \
					SUEWS_AnthropogenicEmissions.o \
					SUEWS_BiogenCO2.o


# modules under rapid development
TEST =      SUEWS_Calculations.o


# some test code
# $(info $$variable is [${variable}])
# $(info $$TARGET is [${TARGET}])

# Build main program - main uses MODULES and OTHERS
main: SUEWS_Program.f95 $(UTILS) $(MODULES) $(OTHERS) $(TEST)
			$(CC) SUEWS_Program.f95  -c ; \
			$(CC) SUEWS_Program.o $(UTILS) $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)
			mkdir -p $(SUEWS_build_dir)
			cp $(TARGET) $(SUEWS_build_dir)/.

# Build main program with NETCDF support - main uses MODULES and OTHERS
netcdf: SUEWS_Program.f95 $(UTILS) $(MODULES) $(OTHERS) $(TEST)
				$(CC) SUEWS_Program.f95  -c ; \
				$(CC) SUEWS_Program.o $(UTILS) $(MODULES) $(OTHERS) $(TEST) \
				$(STATICLIB) \
				-L$(NETCDFLIB) -lnetcdf -lnetcdff \
				-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
nc4fr: SUEWS_Program.f95 $(UTILS) $(MODULES) $(OTHERS) $(TEST)
			 $(CC) SUEWS_Program.f95  -c ; \
			 $(CC) SUEWS_Program.o $(UTILS) $(MODULES) $(OTHERS) $(TEST) \
			 $(STATICLIB) \
			 -L$(NETCDFLIB) -Wl,--rpath -Wl,$(NETCDFLIB) -lnetcdff -lnetcdf \
			 -o $(TARGET)


# If UTILS have changed, compile them again
$(UTILS): $(subst .o,.f95, $(UTILS))
						$(CC) -c $(subst .o,.f95, $@)

# If MODULES have changed, compile them again
$(MODULES): $(UTILS) $(subst .o,.f95, $(MODULES))
						$(CC) -c $(subst .o,.f95, $@)

# If OTHERS have changed, compile them again
$(OTHERS): $(UTILS) $(MODULES) $(subst .o,.f95, $(OTHERS))
					 $(CC) -c $(subst .o,.f95, $@)


# If TEST have changed, compile them again
$(TEST): $(subst .o,.f95, $(TEST))
				 $(CC) -c $(subst .o,.f95, $@)

# If wanted, clean all *.o files after build
clean:
			 -rm -rf *.o *.mod *.dSYM $(TARGET)
