import os
from setuptools import setup, find_packages
import sys

try:
    import numpy
except ImportError:
    raise ImportError(
        "Please install NumPy first, or use the Anaconda Python Distribution "
        "(https://store.continuum.io/cshop/anaconda/) which comes with NumPy "
        "installed."
    )



version_py = os.path.join(os.path.dirname(__file__), 'metaseq', 'version.py')
version = open(version_py).read().split('=')[-1].strip().replace('"','')

long_description = open('README.rst').read()
setup(
        name='metaseq',
        version=version,
        description="Integrative analysis of high-thoughput sequencing data",
        #long_description=long_description,
        license="MIT",
        install_requires=['bx-python>=0.7.1', 'matplotlib>=1.3.1',
                          'pysam>=0.7', 'pandas>=0.13.1', 'pycurl',
                          'pybedtools>=0.6.4', 'gffutils>=0.8', 'urlgrabber',
                          'argparse', 'PyYAML',],
        packages=find_packages(),
        package_data={
            'metaseq':[
                'test/data/gdc*',
                'test/data/make_examples_from_pybedtools.py',
                'test/data/x.*',
            ]
        },
        scripts=[
            'metaseq/scripts/download_metaseq_example_data.py',
            'metaseq/scripts/metaseq-cli',
            'metaseq/scripts/speedtest.py',
        ],
        author='Ryan Dale',
        author_email='dalerr@niddk.nih.gov',
        url='http://github.com/daler/metaseq',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'Intended Audience :: System Administrators',
            'Operating System :: POSIX',
            'Operating System :: MacOS :: MacOS X',
            'Environment :: Console',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.', 
        ]
    )
