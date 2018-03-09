from setuptools import setup, Distribution
import platform


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
          'supy': [lib_name, '*.json']
      },
      distclass=BinaryDistribution,
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'f90nml'
      ],
      include_package_data=True,
      zip_safe=False)
