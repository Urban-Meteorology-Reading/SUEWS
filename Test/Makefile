# -*- makefile -*-
.PHONY: main clean test pip supy

SUEWS_build_base = ./build/

# system specific settings
ifeq ($(OS),Windows_NT)
			SUEWS_build_dir = $(addprefix $(SUEWS_build_base), win64)
else
			UNAME_S := $(shell uname -s)
			ifeq ($(UNAME_S),Linux) # Linux
				SUEWS_build_dir = $(addprefix $(SUEWS_build_base), Linux)
			endif
			ifeq ($(UNAME_S),Darwin) # macOS
				SUEWS_build_dir = $(addprefix $(SUEWS_build_base), macOS)
			endif
endif



SUEWS_dir = ../SUEWS-SourceCode

makefile = Makefile.gfortran

SuPy_dir = supy-driver

PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)



# test the current version
test:
	# $(MAKE) -C $(SUEWS_dir) -f $(makefile) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) -f $(makefile) main; # make SUEWS with the `main` recipe
	cd code; python test_dev.py

clean:
	rm -rf *.zip

# include common settings
include $(SUEWS_dir)/include.common