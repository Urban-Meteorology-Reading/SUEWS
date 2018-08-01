# compilers
FC = gfortran $(FFLAGS)
CC = gcc $(CFLAGS)

file       := suews_ctrl_const.f95 # get program name form source code
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

UTILS = suews_ctrl_const.o  \
				suews_util_stringmod.o \
				suews_util_qsort.o\
				suews_util_time.o \
				suews_util_meteo.o \
				suews_util_datetime.o \
				suews_util_minpack.o

MODULES =	suews_phys_narp.o \
					suews_ctrl_input.o  \
					suews_phys_atmmoiststab.o \
					suews_phys_bluews.o \
					suews_phys_waterdist.o \
					suews_phys_snow.o \
					suews_phys_dailystate.o \
					suews_phys_estm.o \
					suews_phys_anohm.o \
					suews_phys_ohm.o \
					suews_ctrl_output.o \
					suews_ctrl_driver.o

# Rest of the files including modules and functions which are independent
OTHERS =  suews_ctrl_translate.o \
					suews_phys_lumps.o \
					suews_phys_resist.o \
					suews_ctrl_error.o \
					suews_phys_evap.o \
					suews_ctrl_init.o \
					suews_phys_anemsn.o \
					suews_phys_biogenco2.o

# modules under rapid development
TEST =  suews_ctrl_calculations.o


# some test code
# $(info $$variable is [${variable}])
# $(info $$TARGET is [${TARGET}])

# Build main program - main uses MODULES and OTHERS
main: suews_program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			$(FC) suews_program.f95  -c ; \
			$(FC) suews_program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)
			mkdir -p $(SUEWS_build_dir)
			cp -f $(TARGET) $(SUEWS_build_dir)/.

# Build main program with NETCDF support - main uses MODULES and OTHERS
netcdf: suews_program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
				$(FC) suews_program.f95  -c ; \
				$(FC) suews_program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
				$(STATICLIB) \
				-L$(NETCDFLIB) -lnetcdf -lnetcdff \
				-o $(TARGET)

# Build main program with NETCDF support - main uses MODULES and OTHERS
nc4fr: suews_program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			 $(FC) suews_program.f95  -c ; \
			 $(FC) suews_program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			 $(STATICLIB) \
			 -L$(NETCDFLIB) -Wl,--rpath -Wl,$(NETCDFLIB) -lnetcdff -lnetcdf \
			 -o $(TARGET)

# Build main program with coverage checking options enabled
cov: suews_program.f95 $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST)
			$(FC) suews_program.f95  -c ; \
			$(FC) suews_program.o $(C_MODULE) $(UTILS) $(DATETIME) $(MODULES) $(OTHERS) $(TEST) \
			$(STATICLIB) \
			-o $(TARGET)
			mkdir -p $(SUEWS_build_dir)
			cp -f $(TARGET) $(SUEWS_build_dir)/.


# Build main program and run test cases
check:
		$(MAKE) clean;
		$(MAKE) main;
		cd $(SUEWS_test_dir)/code; python 1.test_dev.py


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
