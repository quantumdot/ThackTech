from setuptools import setup

setup(name='ThackTech',
      version='0.1',
      description='Libraries for process management, pipelineing, and other goodies.',
      url='https://github.com/quantumdot/ThackTech',
      author='Joshua K. Thackray',
      author_email='thackray@rutgers.edu',
      license='MIT',
      packages=['ThackTech'],
      install_requires=[
		  'numpy',
		  'scipy',
          'intervaltree',
		  'pysam',
          'pybedtools',
		  'glob',
		  'dill',
		  'multiprocess',
		  'tabulate',
		  'sqlite3',
      ],
      zip_safe=False)