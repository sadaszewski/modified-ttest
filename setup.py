#
# Copyright (C) Stanislaw Adaszewski, 2023
# Email: s.adaszewski@adared.ch, s.adaszewski@gmail.com
# License: GPL v3
#

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension('modified_ttest.ext', ['modified_ttest/mod_ttest_c.c',
            'modified_ttest/mod_ttest.pyx'],
        include_dirs=[ numpy.get_include() ],
        libraries=[],
        library_dirs=[],
        #extra_compile_args=['-std=c++14'],
        #language='c++'
    )
]

setup(
    name='modified-ttest',
    version='1.0',
    author='Stanislaw Adaszewski',
    author_email='s.adaszewski@adared.ch',
    url='https://github.com/sadaszewski/modified-ttest',
    packages=['modified_ttest'],
    license='The GNU General Public License v3.0',
    description='Modified t-test reimplementation from R\'s SpatialPack (paper by Dutilleul)',
    long_description='Modified t-test reimplementation from R\'s SpatialPack (paper by Dutilleul)',
    scripts=[],
    install_requires=[
        'numpy',
    ],
    ext_modules=cythonize(extensions)
)