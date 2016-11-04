# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='SONATA',
    version='0.0.1',
    description='Crosssectional Rotor Blade Preprocessor',
    long_description=readme,
    author='Tobias Pflumm',
    author_email='tobias.pflumm@tum.de',
    url='https://gitlab.lrz.de/gu32kij/SONATA',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
