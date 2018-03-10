# -*- makefile -*-

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON = /c/Users/sunt05/Anaconda2/python.exe
	F2PY_PY= /c/Users/sunt05/Anaconda2/Scripts/f2py.py
	F2PY_EXE = $(PYTHON) $(F2PY_PY)
	TARGET=$(MODULE).pyd
else
	UNAME_S := $(shell uname -s)
	TARGET=$(MODULE).so
	ifeq ($(UNAME_S),Linux) # Linux
		F2PY_EXE = f2py
	endif
	ifeq ($(UNAME_S),Darwin) # macOS
		F2PY_EXE = f2py
	endif
	PYTHON=python
endif

MODULE=SUEWS_driver
.PHONY:main,clean


LDSHARED = -static
F2PY = $(F2PY_EXE)
SUEWS_dir = ./SUEWS-AnOHM
# F2PY_FLAGS  = --quiet
F2PY_FLAGS  = --verbose
F2PY_F90FLAGS = --f90flags='-Wall -Wtabs -cpp' --opt='-O3' --build-dir build-f2py
# F2PY_FLAGS += --debug-capi
# F2PY_FLAGS += -DF2PY_REPORT_ON_ARRAY_COPY=1
# F2PY_FLAGS += -DF2PY_REPORT_ATEXIT

# All the files which include modules used by other modules (these therefore
# needs to be compiled first)
FILES = LUMPS_Module_constants.f95  \
				SUEWS_driver.f95

# ${MODULE}.so: ${MODULE}.pyf ${MODULE}.f90
# 	${F2PY} ${F2PY_FLAGS} -m ${MODULE} -c $< ${MODULE}.f90

main:
	# $(MAKE) -C $(SUEWS_dir) clean;
	# $(MAKE) -C $(SUEWS_dir) main; # make SUEWS with the `main` recipe
	-ln -sf $(SUEWS_dir)/*.o .; # link objects
	-ln -sf $(SUEWS_dir)/*.a .; # link objects
	-ln -sf $(SUEWS_dir)/*.mod .; # link objects
	-ln -sf $(SUEWS_dir)/*.f95 .; # link objects
	ar crs libSUEWS.a *.o *.a
	-rm $(subst .f95,.o, $(FILES));
	$(F2PY) $(F2PY_FLAGS) -m $(MODULE) -c $(F2PY_F90FLAGS) $(FILES) *.a;
	-mv $(TARGET) SuPy/.
	-mv libSUEWS.a SuPy/.
	-rm -rf *.o *.mod *.f95 *.a *.dSYM


# If wanted, clean all *.o files after build
clean:
			# $(info $$TARGET is [$(TARGET)])
			 -rm -rf *.o *.mod *.dSYM $(TARGET) SuPy/$(MODULE).*;
			 -$(PYTHON) setup.py clean --all

# clean all existing builds, rebuild f2py libs, build wheels and submit
pip:
			$(MAKE) clean;
			$(MAKE) main;
			$(PYTHON) setup.py bdist_wheel upload
