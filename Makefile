# -*- makefile -*-
.PHONY: main clean test pip supy

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON_exe = python.exe
	# F2PY_PY= /c/Users/sunt05/Anaconda2/Scripts/f2py.py
	# F2PY_EXE = $(PYTHON) $(F2PY_PY)
	TARGET=$(MODULE).pyd
else
	UNAME_S := $(shell uname -s)
	TARGET=$(MODULE).so

	ifeq ($(UNAME_S),Linux) # Linux
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

	ifeq ($(UNAME_S),Darwin) # macOS
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

endif

MODULE=SUEWS_driver

SUEWS_dir = SUEWS-SourceCode

release_dir = ReleaseRepo

makefile = Makefile.gfortran

SuPy_dir = supy-driver

PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)



# make fortran exe
main:
	$(MAKE) -C $(SUEWS_dir) -f $(makefile) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) -f $(makefile) main; # make SUEWS with the `main` recipe
	-rm -rf *.o *.mod *.f95 *.a *.dSYM

# make fortran exe and run test cases
test:
	$(MAKE) -C $(SUEWS_dir) -f $(makefile) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) -f $(makefile) test; # make SUEWS with the `main` recipe

# make fortran exe, run test cases and pack release archive
release:
	$(MAKE) pip
	$(MAKE) main
	$(MAKE) -C $(release_dir) pack; # clean Fortran SUEWS build

# make supy dist
driver:
	$(info $$PYTHON is [${PYTHON}])
	$(MAKE) -C $(SuPy_dir) suews; # clean Fortran SUEWS build
	$(MAKE) -C $(SuPy_dir) main; # make SUEWS with the `main` recipe

pip:
	pip install pipreqs
	pipreqs $(release_dir) --savepath requirements.txt
	pip install -r requirements.txt
	rm -rf requirements.txt

# If wanted, clean all *.o files after build
clean:
	$(MAKE) -C $(SUEWS_dir) clean
	$(MAKE) -C $(SuPy_dir) clean
