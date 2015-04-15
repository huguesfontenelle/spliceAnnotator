#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages
from os import path
import os
from pip.req import parse_requirements
import uuid

VERSION="0.1"

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements("requirements.txt", session=uuid.uuid1())

# reqs is a list of requirement
# e.g. ['django==1.5.1', 'mezzanine==1.4.6']
reqs = [str(ir.req) for ir in install_reqs]

setup(name="splice",
    version=VERSION,
    description="Splice annotation",
    long_description=open("README.md").read(),
    author="Hugues Fontenelle",
    author_email="hugues.fontenelle@medisin.uio.no",
    url="",
    license="Unknown",
    packages=find_packages(),
    include_package_data=True,
    install_requires=reqs,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'License :: Other/Proprietary License',
        'Natural Language :: Norwegian',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )

