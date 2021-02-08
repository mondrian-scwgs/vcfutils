#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'click>=7.0',
    'pandas',
    'pyyaml',
    'numpy',
    'pytest',
    'Sphinx==1.8.5',
]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Shah Lab",
    author_email='todo@todo.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="todo: add description",
    entry_points={
        'console_scripts': [
            'vcfutils=vcfutils.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='vcfutils',
    name='vcfutils',
    packages=find_packages(include=['vcfutils', 'vcfutils.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mondrian-scwgs/vcfutils',
    version="v0.1.4",
    zip_safe=False,
)
