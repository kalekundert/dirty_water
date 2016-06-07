#!/usr/bin/env python3
# encoding: utf-8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import re
with open('dirty_water/__init__.py') as file:
    version_pattern = re.compile("__version__ = '(.*)'")
    version = version_pattern.search(file.read()).group(1)

with open('README.rst') as file:
    readme = file.read()

setup(
    name='dirty_water',
    version=version,
    author='Kale Kundert',
    author_email='kale.kundert@ucsf.edu',
    description='',
    long_description=readme,
    url='https://github.com/kalekundert/dirty_water',
    packages=[
        'dirty_water',
    ],
    include_package_data=True,
    install_requires=[
        'nonstdlib',
    ],
    license='MIT',
    zip_safe=False,
    keywords=[
        'dirty_water',
    ],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries',
    ],
)
