from setuptools import Distribution
from numpy.distutils.core import Extension, setup
import platform
import glob
import os
from pathlib import Path
import subprocess
import shutil
from nonstopf2py import f2py
# from gen_suewsdrv import merge_source

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
    [
        'suews_ctrl_const.f95',
        'suews_ctrl_error.f95',
        'suews_phys_narp.f95',
        'suews_phys_atmmoiststab.f95',
        'suews_phys_resist.f95',
        'suews_phys_evap.f95',
        'suews_phys_snow.f95',
        'suews_phys_dailystate.f95',
        'suews_phys_lumps.f95',
        'suews_phys_anemsn.f95',
        'suews_phys_rslprof.f95',
        'suews_phys_biogenco2.f95',
        'suews_phys_ohm.f95',
        'suews_phys_waterdist.f95',
        'suews_ctrl_driver.f95',
    ]
]
all_f95 = glob.glob(os.path.join(dir_f95, '*.f95'))
exclude_f95 = [
    os.path.join(dir_f95, f)
    for f in
    [
        'suews_c_wrapper.f95',
        'suews_ctrl_sumin.f95',
        'suews_program.f95',
    ]
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

# # combine for files to use:
# file_all_f95 = 'suews_all.f95'
# # with open(file_all_f95, 'wb') as wfd:
# #     for f in src_f95:
# #         with open(f, 'rb') as fd:
# #             shutil.copyfileobj(fd, wfd)
# # directory of SUEWS source code
# # this dir is included as a git submodule so DON'T make ANY change there
# file_all_f95 = 'suews_all.f95'
# path_src_SUEWS = Path(dir_f95).resolve()
# print(path_src_SUEWS)
# # 4. generate SUEWS related source files from $dir_src_SUEWS and add them to $dir_WRF_SUEWS
# path_sf_suewsdrv = Path(file_all_f95)
# print(f'calling merge_source to generate {file_all_f95}')
# merge_source(path_src_SUEWS, path_sf_suewsdrv)
# # for f in target_f95:
# #     print(f)


def readme():
    f = """
    `supy_driver` is `F2PY`-based python binary package for `supy` with `SUEWS` as the computation core.
    """
    return f
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
    with open(str(path_pkg_init), 'w') as fm:
        fm.write("__version__='{ver}'".format(ver=ver))

    return ver


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False


# # print('will build', lib_name)
# for x in other_obj:
#     print(x)

ext_modules = [
    Extension('supy_driver.suews_driver',
              target_f95,
              extra_compile_args=['-D_POSIX_C_SOURCE=200809L'],
              extra_f90_compile_args=['-cpp'],
              f2py_options=[
                  # '--quiet',
                #   '--debug-capi',
                  # ('-DF2PY_REPORT_ATEXIT' if sysname == 'Linux' else ''),
              ],
              extra_objects=other_obj,
              extra_link_args=[('' if sysname == 'Linux' else '-static')])]

setup(name='supy_driver',
      # update version info here!
      version=get_suews_version(ver_minor=3),
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


# check latest build
path_dir_driver = Path(__file__).resolve().parent
list_wheels = [str(x) for x in path_dir_driver.glob('dist/*whl')]
fn_wheel = sorted(list_wheels, key=os.path.getmtime)[-1]
# print(list_wheels, fn_wheel)

# use auditwheel to repair file name for Linux
if sysname == 'Linux':
    # path_dir_driver = Path(__file__).resolve().parent
    # list_wheels = [str(x) for x in path_dir_driver.glob('dist/*whl')]
    # fn_wheel = sorted(list_wheels, key=os.path.getmtime)[-1]
    # print(list_wheels, fn_wheel)
    subprocess.call(["auditwheel", "repair", fn_wheel])
    subprocess.call(["ls", "-lrt"])


# change compiler settings
if sysname == 'Windows':
    os.remove('setup.cfg')
