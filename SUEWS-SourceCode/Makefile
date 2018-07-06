# compilers
FC = gfortran $(FFLAGS)
CC = gcc $(CFLAGS)

file       := SUEWS_ctrl_const.f95 # get program name form source code
variable   := $(shell grep 'progname' ${file})
TARGET_str := $(lastword $(subst =, ,${variable})) # program name with quotes
TARGET     := $(strip $(subst $\', ,${TARGET_str})) # program name

# static flag for different OS to correctly link static libs
# so gfortran dependency can be relaxed
# but netCDF is still linked in the dynamic way as suggested by UCAR
SUEWS_build_base = ../ReleaseRepo/build/
SUEWS_test_dir = ../BenchmarkTest/

ifeq ($(OS),Windows_NT)
			C_MODULE = strptime.o # this is needed by datetime_module: mingw lacks `strptime`
			STATIC = -static -mcmodel=medium# mingw
			STATICLIB =
			SUEWS_build_dir = $(addprefix $(SUEWS_build_base), Win10x64)
else
			C_MODULE =
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

FFLAGS = $(STATIC) -g -pg -Wall -Wtabs -fbounds-check -cpp \
					-Wno-unused-dummy-argument -Wno-unused-variable \
					-fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal
CFLAGS =


# netCDF-related settings:
netcdf: NETCDFINC = `nc-config --includedir` # ordinary path for netCDF directories
netcdf:	NETCDFLIB = `nc-config --libdir` # ordinary path for netCDF directories

nc4fr: NETCDFINC = /home/xlinfr/apps/lib4cdo/include # path only valid with Fredrik's HPC
nc4fr: NETCDFLIB = /home/xlinfr/apps/lib4cdo/lib # path only valid with Fredrik's HPC

netcdf nc4fr:	FFLAGS += -Dnc=1 -I$(NETCDFINC) # options for netcdf build

# coverage testing related settings:
cov: FFLAGS += -fprofile-arcs -ftest-coverage # add options for coverage testing
cov: suffix := _cov
cov: TARGET := $(strip $(TARGET))$(strip $(suffix))


clean: suffix :=
clean: TARGET := $(strip $(TARGET))$(strip $(suffix))


# All the files which include modules used by other modules (these therefore
# needs to be compiled first)

UTILS = SUEWS_ctrl_const.o  \
				SUEWS_util_stringmod.o \
				SUEWS_util_qsort.o\
				SUEWS_util_time.o \
				SUEWS_util_meteo.o \
				SUEWS_util_datetime.o \
				SUEWS_util_minpack.o

MODULES =	SUEWS_phys_NARP.o \
					SUEWS_ctrl_Input.o  \
					SUEWS_phys_AtmMoistStab.o \
					SUEWS_phys_BLUEWS.o \
					SUEWS_phys_WaterDist.o \
					SUEWS_phys_Snow.o \
					SUEWS_phys_DailyState.o \
					SUEWS_phys_ESTM.o \
					SUEWS_phys_AnOHM.o \
					SUEWS_phys_OHM.o \
					SUEWS_ctrl_Output.o \
					SUEWS_ctrl_Driver.o

# Rest of the files including modules and functions which are independent
OTHERS =  SUEWS_ctrl_Translate.o \
					SUEWS_phys_LUMPS.o \
					SUEWS_phys_Resist.o \
					SUEWS_ctrl_Error.o \
					SUEWS_phys_Evap.o \
					SUEWS_ctrl_Init.o \
					SUEWS_phys_AnEmsn.o \
					SUEWS_phys_BiogenCO2.o

# modules under rapid development
TEST =  SUEWS_ctrl_Calculations.o


# some test code
# $(info $$variable is [${variable}])
# $(info $$TARGET is [${TARGET}])

# Build main program - main uses MODULES and OTHERS
main: SUEWS_Program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			$(FC) SUEWS_Program.f95  -c ; \
			$(FC) SUEWS_Program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)
			mkdir -p $(SUEWS_build_dir)
			cp -f $(TARGET) $(SUEWS_build_dir)/.

# Build main program with NETCDF support - main uses MODULES and OTHERS
netcdf: SUEWS_Program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
				$(FC) SUEWS_Program.f95  -c ; \
				$(FC) SUEWS_Program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
				$(STATICLIB) \
				-L$(NETCDFLIB) -lnetcdf -lnetcdff \
				-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
nc4fr: SUEWS_Program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			 $(FC) SUEWS_Program.f95  -c ; \
			 $(FC) SUEWS_Program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			 $(STATICLIB) \
			 -L$(NETCDFLIB) -Wl,--rpath -Wl,$(NETCDFLIB) -lnetcdff -lnetcdf \
			 -o $(TARGET)

# Build main program with coverage checking options enabled
cov: SUEWS_Program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			$(FC) SUEWS_Program.f95  -c ; \
			$(FC) SUEWS_Program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)
			mkdir -p $(SUEWS_build_dir)
			cp -f $(TARGET) $(SUEWS_build_dir)/.


# Build main program and run test cases
check:
		$(MAKE) clean;
		$(MAKE) main;
		cd $(SUEWS_test_dir); python test_dev.py


# If UTILS have changed, compile them again
$(UTILS): $(subst .o,.f95, $(UTILS))
						$(FC) -c $(subst .o,.f95, $@)

$(C_MODULE): $(subst .o,.c, $(C_MODULE))
						$(CC) -c $(subst .o,.c, $@)

# If MODULES have changed, compile them again
$(MODULES): $(UTILS) $(subst .o,.f95, $(MODULES))
						$(FC) -c $(subst .o,.f95, $@)

# If OTHERS have changed, compile them again
$(OTHERS): $(UTILS) $(MODULES) $(subst .o,.f95, $(OTHERS))
					 $(FC) -c $(subst .o,.f95, $@)

# If TEST have changed, compile them again
$(TEST): $(subst .o,.f95, $(TEST))
				 $(FC) -c $(subst .o,.f95, $@)

# If wanted, clean all *.o files after build
clean:
	$(info $$TARGET is [${TARGET}])
	-rm -rf *.o *.mod *.dSYM *.gcda *.gcno ${TARGET}
