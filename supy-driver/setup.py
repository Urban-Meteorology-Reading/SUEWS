# from setuptools import setup, Distribution
from setuptools import Distribution
from numpy.distutils.core import Extension, setup
import platform
import glob
import os
from pathlib import Path
import subprocess
import shutil


# wrap OS-specific `SUEWS_driver` libs
sysname = platform.system()
lib_basename = 'supy_driver'
if sysname == 'Windows':
    lib_name = lib_basename + '.pyd'
elif sysname == 'Darwin':
    lib_name = lib_basename + '.so'
elif sysname == 'Linux':
    lib_name = lib_basename + '.so'

# change compiler settings
if sysname == 'Windows':
    shutil.copyfile('win-setup.cfg', 'setup.cfg')

# load SUEWS Fortran source files
dir_f95 = '../SUEWS-SourceCode'
target_f95 = [
    os.path.join(dir_f95, f)
    for f in
    ['suews_ctrl_const.f95',
     'suews_ctrl_driver.f95']]
all_f95 = glob.glob(os.path.join(dir_f95, '*.f95'))
exclude_f95 = [
    os.path.join(dir_f95, f)
    for f in
    ['suews_c_wrapper.f95',
     'suews_ctrl_sumin.f95',
     'suews_program.f95']
]
other_f95 = list(
    set(all_f95)
    - set(target_f95)
    - set(exclude_f95)
)
other_obj = [f.replace('.f95', '.o') for f in other_f95]
if sysname == 'Windows':
    other_obj.append(os.path.join(dir_f95, 'strptime.o'))

src_f95 = target_f95 + other_f95
# for f in target_f95 + other_obj:
#     print(f)


def readme():
    with open('README.md') as f:
        return f.read()
# dir_source='SUEWS-SourceCode'
# path_source = Path(dir_source)
# str(path_source)


def get_suews_version(dir_source=dir_f95, ver_minor=2):
    path_source = Path(dir_source)
    path_makefile = (path_source / 'include.common')
    # identify `file` to retrieve version
    with open(str(path_makefile)) as fm:
        for line in fm:
            if 'file ' in line:
                file = line.split(':=')[-1].split('#')[0].strip()

    # get version from `file`
    path_constfile = (path_source / file)
    with open(str(path_constfile)) as fm:
        for line in fm:
            if 'progname' in line:
                ver = line.split('SUEWS_V')[-1].replace("'", '').strip()
                ver += str(ver_minor)

    # cast `ver` to the driver package
    path_pkg_init = Path('.')/lib_basename/'version.py'
    with open(str(path_pkg_init),'w') as fm:
        fm.write("__version__='{ver}'".format(ver=ver))

    return ver


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False


# print('will build', lib_name)


ext_modules = [
    Extension('supy_driver.suews_driver',
              target_f95,
              extra_f90_compile_args=['-cpp'],
              f2py_options=[
                  # '--quiet',
                  # '--debug-capi',
                  # ('-DF2PY_REPORT_ATEXIT' if sysname == 'Linux' else ''),
              ],
              extra_objects=other_obj,
              extra_link_args=[('' if sysname == 'Linux' else '-static')])]

setup(name='supy_driver',
      version=get_suews_version(ver_minor=24),
      description='the SUEWS driver driven by f2py',
      long_description=readme(),
      url='https://github.com/sunt05/SuPy',
      author='Ting Sun',
      author_email='ting.sun@reading.ac.uk',
      # license='GPL-V3.0',
      packages=['supy_driver'],
      package_data={
          'supy_driver': [
              # lib_name,
           # '*.json'
          ]
      },
      distclass=BinaryDistribution,
      ext_modules=ext_modules,
      python_requires='>=3.5',
      install_requires=[
          'numpy>=1.15.2'
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)


# use auditwheel to repair file name
path_file=Path(__file__).resolve().parent
print(path_file)
if sysname == 'Linux':
    fn_wheel = sorted(path_file.glob('dist/*whl'), key=os.path.getmtime)[-1]
    subprocess.call(["auditwheel", "repair", fn_wheel])


# change compiler settings
if sysname == 'Windows':
    os.remove('setup.cfg')
