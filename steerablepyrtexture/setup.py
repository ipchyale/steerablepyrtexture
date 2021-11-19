#! /usr/bin/env python

from setuptools import setup, Extension
import importlib
import os


setup(
    name='steerablepyrtexture',
    version='0.5',
    description=('Modifation of pyrtools to include texture analysis'),
    license=' ',
    url=' ',
    author='Nicholas Rogers',
    author_email='nickr7@gmail.com',
    keywords='multi-scale image-processing',
    packages=['steerablepyrtexture', 'steerablepyrtexture.pyramids', 'steerablepyrtexture.tools'],
    package_data={'': ['*.h', 'LICENSE']},
    install_requires=['numpy>=1.1',
                      'scipy>=0.18',
                      'matplotlib>=1.5',
                      'Pillow>=3.4',
                      'tqdm>=4.29',
                      'requests>=2.21'],
    )
