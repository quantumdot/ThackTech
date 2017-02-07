from setuptools import setup

setup(name='ThackTech',
      version='0.1',
      description='Libraries for process management, pipelineing, and other goodies.',
      #url='http://github.com/storborg/funniest',
      author='Josh Thackray',
      author_email='thackray@rutgers.edu',
      license='MIT',
      packages=['funniest'],
      install_requires=[
		  'numpy',
		  'scipy',
          'intervaltree',
		  'pysam',
		  'glob',
		  'dill',
		  'multiprocess',
		  'tabulate',
		  'sqlite3',
      ],
      zip_safe=False)