import ez_setup
ez_setup.use_setuptools()
from setuptools import setup
import sys
import os
import numpy

version_py = os.path.join(os.path.dirname(__file__), 'metaseq', 'version.py')
version = open(version_py).read().split('=')[-1].strip().replace('"','')

long_description = open('README.rst').read()
setup(
        name='metaseq',
        version=version,
        long_description=long_description,
        install_requires=['bx-python', 'numpy', 'HTSeq', 'matplotlib', 'scipy',
                          'scikits.learn', 'pysam', 'statsmodels',
                          'pybedtools', 'gffutils', 'argparse'],
        packages=['metaseq', 'metaseq.test', 
                  'metaseq.test.data',
                  'metaseq.integration'],
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
