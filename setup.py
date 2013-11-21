import ez_setup
ez_setup.use_setuptools()
from setuptools import setup
import sys
import os
import numpy


version_py = os.path.join(os.path.dirname(__file__), 'metaseq', 'version.py')
version = open(version_py).read().split('=')[-1].strip().replace('"','')


long_description = """
metaseq is a Python framework for genomic data analysis (primarily
high-throughput sequencing, but can be used for much more).  It ties
together other frameworks like BEDTools/pybedtools, samtools/pysam, bx-python,
HTSeq, gffutils, and matplotlib.
"""
setup(
        name='metaseq',
        version=version,
        long_description=long_description,
        install_requires=['bx-python', 'numpy', 'HTSeq', 'matplotlib', 'scipy',
                          'scikit-learn', 'pybedtools', 'gffutils',
                          'argparse'],
        packages=['metaseq', 'metaseq.test', 'metaseq.test.data', 'metaseq.integration'],
        package_data={'metaseq':['test/data/*']},
        scripts=[
            'metaseq/scripts/download_metaseq_example_data.py',
            'metaseq/scripts/metaseq-cli',
        ],
        author='Ryan Dale',
        author_email='dalerr@niddk.nih.gov',
        url='http://github.com/daler/metaseq',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
