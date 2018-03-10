# from setuptools import setup, Distribution
from setuptools import Distribution
from numpy.distutils.core import Extension, setup
import platform
import glob


def readme():
    with open('README.rst') as f:
        return f.read()


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False


# wrap OS-specific `SUEWS_driver` libs
sysname = platform.system()
if sysname == 'Windows':
    lib_name = 'SUEWS_driver.pyd'
elif sysname == 'Darwin':
    lib_name = 'SUEWS_driver.so'
else:
    lib_name = 'SUEWS_driver.so'

target_f95 = ['SUEWS-AnOHM/LUMPS_Module_constants.f95',
              'SUEWS-AnOHM/SUEWS_driver.f95']
other_f95 = set(glob.glob('SUEWS-AnOHM/*f95')) - \
    set(target_f95) - set(['SUEWS-AnOHM/SUEWS_C_wrapper.f95'])
other_f95 = list(other_f95)
other_obj = [f.replace('.f95', '.o') for f in other_f95]
src_f95 = target_f95 + other_f95
print src_f95,other_obj
ext_modules = [
    Extension('SUEWS_driver',
              target_f95,
              # src_f95,
              extra_f90_compile_args=['-cpp'],
              # library_dirs=['./SUEWS-AnOHM'],
              extra_objects=other_obj,
              extra_link_args=['-static'])]

setup(name='supy',
      version='0.1.3a7',
      description='the SUEWS model that speaks python',
      long_description=readme(),
      url='https://github.com/sunt05/SuPy',
      author='Ting Sun',
      author_email='ting.sun@reading.ac.uk',
      license='GPL-V3.0',
      packages=['supy'],
      package_data={
          'supy': [
              # lib_name,
           '*.json'
          ]
      },
      distclass=BinaryDistribution,
      ext_modules=ext_modules,
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'f90nml'
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
