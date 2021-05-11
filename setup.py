#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy',
                'os_command_py',
                'pdb_manip_py',
                'pdb2pqr-htmd-propka30']

test_requirements = ['pytest>=3', ]

setup(
    version='0.3.0',
    author="Samuel Murail",
    author_email='samuel.murail@u-paris.fr',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Docking_py is a python library allowing a simplified use of the Smina, vina, qvina2 and qvinaw docking software. Docking_py can be easily automatize and scripted.",
    entry_points={
        'console_scripts': [
            'docking_py=docking_py.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v2 (GPLv2)",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='docking_py',
    name='docking_py',
    packages=find_packages(include=['docking_py', 'docking_py.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/samuelmurail/docking_py',
    zip_safe=False,
)
