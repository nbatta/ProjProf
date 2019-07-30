from distutils.core import setup, Extension
import os

setup(name='projprof',
      version='0.1',
      description='Radial profile projection tool for Sunyeav-Zeldovich observations',
      url='https://github.com/nbatta/ProjProf',
      author='Nicholas Battaglia',
      author_email='nicholas.battaglia@gmail.com',
      license='BSD-2-Clause',
      packages=['ProjProf'],
      package_dir={'ProjProf':'src'},
      zip_safe=False)
