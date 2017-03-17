from setuptools import setup


setup(  name='ThackTech',
        version='0.1',
        description='Libraries for process management, pipelineing, and other goodies.',
        url='https://github.com/quantumdot/ThackTech',
        author='Joshua K. Thackray',
        author_email='thackray@rutgers.edu',
        license='MIT',
        install_requires=[
            'numpy',
            'scipy',
            'intervaltree',
            'pysam',
            'pybedtools',
            'dill',
            'multiprocess', #this is different than the built in multiprocessing!
            'tabulate',
            #'glob', 
            #'sqlite3',
        ],
        scripts=[
            "scripts/bedToBedGraph.py",
            "scripts/wigToBedGraph.py",
            "scripts/cleanBedGraph.py",
            "scripts/profiler.py"
        ],
        zip_safe=False)