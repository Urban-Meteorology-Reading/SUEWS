from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='supy',
      version='0.1.1-alpha',
      description='the SUEWS model that speaks python',
      long_description=readme(),
      url='https://github.com/sunt05/SuPy',
      author='Ting Sun',
      author_email='ting.sun@reading.ac.uk',
      license='GPL-V3.0',
      packages=['supy'],
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'f90nml'
      ],
      zip_safe=False)
