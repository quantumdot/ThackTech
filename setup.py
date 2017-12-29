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
            'pandas',
            'matplotlib',
            'metaseq',
            'intervaltree',
            'pysam',
            'pybedtools',
            'dill',
            'multiprocess', #this is different than the built in multiprocessing!
            'tabulate',
            #'HTseq'
            #'glob', 
            #'sqlite3',
        ],
        scripts=[
            "scripts/bedToBedGraph.py",
            "scripts/wigToBedGraph.py",
            "scripts/cleanBedGraph.py",
            "scripts/profiler.py",
            "scripts/correlator.py",
            "scripts/bedstats.py",
            "scripts/bowtie_alignment_pipeline.py",
            "scripts/macs_peakcall_pipeline.py",
            "scripts/epic_peakcall_pipeline.py",
            "scripts/subsample_fastq.py",
            "scripts/keggTermToGeneBed.py",
            "scripts/analyze_radial_intensity.py",
            "scripts/multi_venn.py",
            "scripts/designCNVprobes.py",
            "scripts/subsample_fastq.py",
            "scripts/tuxedo2_pipeline.py",
            "scripts/ginfo_inspect.py",
            "scripts/manifest_manipulate.py"
        ],
        zip_safe=False)