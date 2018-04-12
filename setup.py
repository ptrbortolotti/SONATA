# -*- coding: utf-8 -*-

from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='SONATA',
    version='0.0.1',
    description='Crosssectional Rotor Blade Preprocessor',
    author='Tobias Pflumm',
    author_email='tobias.pflumm@tum.de',
    url='https://gitlab.lrz.de/gu32kij/SONATA',
    ext_modules=cythonize("SONATA\cbm\mesh\mesh_byprojection.pyx"),
)
