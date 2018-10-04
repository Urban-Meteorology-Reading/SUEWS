# -*- makefile -*-
.PHONY: main clean test pip supy

# OS-specific configurations
OSFLAG 				:=
ifeq ($(OS),Windows_NT)
	OSFLAG += -D WIN32
	PYTHON:=python
	ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
		OSFLAG += -D AMD64
	endif
	ifeq ($(PROCESSOR_ARCHITECTURE),x86)
		OSFLAG += -D IA32
	endif
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAG += -D LINUX
		pyver := 37
		PYTHON:=/opt/python/cp${pyver}-cp${pyver}m/bin/python
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAG += -D OSX
		PYTHON:=python
	endif
		UNAME_P := $(shell uname -p)
	ifeq ($(UNAME_P),x86_64)
		OSFLAG += -D AMD64
	endif
		ifneq ($(filter %86,$(UNAME_P)),)
	OSFLAG += -D IA32
		endif
	ifneq ($(filter arm%,$(UNAME_P)),)
		OSFLAG += -D ARM
	endif
endif

PYTHON_ver:= $(shell $(PYTHON) --version)


MODULE=suews_driver

SUEWS_dir = ../SUEWS-SourceCode

# SuPy_dir = supy-driver

# PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)

# $(info $$UNAME_S is ${UNAME_S})
# make suews_driver
main:
	$(MAKE) show
	$(PYTHON) -m pip install -r requirements.txt
	$(PYTHON) -m pip install auditwheel
	$(PYTHON) setup.py bdist_wheel


# make fortran exe
suews:
	$(MAKE) -C $(SUEWS_dir) -f Makefile.gfortran clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) -f Makefile.gfortran main; # make SUEWS with the `main` recipe
	-rm -rf *.o *.mod *.f95 *.a *.dSYM

show:
	@echo $(info $$OSFLAG is [${OSFLAG}])
	@echo $(info $$PYTHON is ${PYTHON_ver})


# If wanted, clean all *.o files after build
clean:
	$(MAKE) show
	rm -rf ${MODULE}/*.so ${MODULE}/*.pyc ${MODULE}/__pycache__ ${MODULE}/*.dSYM
	 -$(PYTHON) setup.py clean --all

# clean all existing builds, rebuild f2py libs, build wheels and test
test:
	$(MAKE) show;
	$(PYTHON) setup.py test

# upload wheels to pypi using twine
upload:
	$(MAKE) show;
	twine upload --skip-existing dist/*whl
	# $(shell twinup wheelhouse/*whl)
