from distutils.core import setup
import setuptools
import os
from pathlib import Path

scripts = []
script_dir = Path(__file__).resolve().parent.joinpath('ThackTech', 'scripts')
for s in script_dir.glob('*.py'):
    if s.name != '__init__.py':
        ep = s.name.replace('.py', '')
        scripts.append('{} = ThackTech.scripts.{}:main'.format(ep,ep))

setup(
    name='ThackTech',
    version='0.1',
    description='Libraries for process management, pipelineing, and other goodies.',
    url='https://github.com/quantumdot/ThackTech',
    author='Joshua K. Thackray',
    author_email='thackray@rutgers.edu',
    license='MIT',
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'metaseq',
        'intervaltree',
        'pysam',
        'pybedtools',
        'dill',
        'multiprocess', #this is different than the built in multiprocessing!
        'tabulate',
        'pdfrw',
        #'HTseq'
        #'glob', 
        #'sqlite3',
    ],
    entry_points={
        'console_scripts': scripts,
    },
    zip_safe=False
)
