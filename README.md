# ThackTech


# Introduction
The `ThackTech` package contains several modules and subpackages for genomic-related analysis. A centerpiece of the packages is the `ThackTech.Pipelines` package, which offers portable pipline systems and standardized analysis modules. By portability, the system attempts to be system agnostic (with the limitation of being mostly linux-bound), but can be run on anything from a laptop, workstation, or even a SLURM-managed cluster.

# Install

## Framework Install
Installation should is simple:

`pip install -e git+ssh://git@github.com/quantumdot/ThackTech.git#egg=ThackTech`

You can also use the `--user` and `--no-deps` options. The repository is currently private, so see the owner for gaining access credentials.

## Configuration
Configuration options are set through ini-style files, and can be located in several locations with order precidence. In order, the configuration files are read from the following locations (with files read later having higher precidence):
1. `<Package_Install_Location>/conf`
2. `/etc/thacktech`
3. `~/.config/thacktech`
4. `$THACKTECH_CONF` environment variable
5. `$PWD` during invocation

There are several configuration types available that affect different parts of the framework (and are detailed below). For each of these, a config file of the name `<config_type>[.default].ini` is searched for. Configuration files can also be registered through the `ThackTech.conf` module.

### genomes
Genome configuration is setup to take advantage of the iGenomes reference setup, but index and other entries may also be setup manually. To specify a genome based on an iGenomes directory specify the following:

```
[mm9] #genome name
size=1.91e9	#effective_genome_size
goldenpath=/mnt/ref/reference/Mus_musculus/UCSC/mm9 #root of reference genome
```

### pipeline
Pipeline configuration elements go here. This mostly deals with pipeline runner configurations

```
[general]
shm_dir: /run/shm
runner: slurm

[slurm_runner]
partition: "main"
nodes: 1
threads: 8
time_limit: "1:00:00"

[parallel_runner]
threads: 8

[serial_runner]
threads: 8
```

### pipeline_modules
Configurations for pipeline modules go in these files, with the section header as the module name, and settings relating to module parameters. An example for the trimmomatic module:

```
[Trimmomatic]
trim_adapt_fa_SE=/mnt/ref/adapters/TruSeq3-SE.fa
trim_adapt_fa_PE=/mnt/ref/adapters/TruSeq3-PE-2.fa
leading=3
trailing=3
sliding_window_width=4
sliding_window_qthresh=15
min_length=25
```

# Pipelines
## Premade Pipelines
Pipeline-related infrastructure is located in the `ThackTech.Pipelines` package.  

There are a few pre-configured pipelines ready for use (these should be installed into your bin dir automagically):
* `bowtie_alignment_pipeline.py`
* `macs_peakcall_pipeline.py`
* `epic_peakcall_pipeline.py`

## Creating your own pipeline

Making your own pipeline is relatively simple:

```python
#import some framework classes
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner

#create a pipeline object
pipeline = AnalysisPipeline('Hello World')

#import an analysis module and add it to the pipeline
from ThackTech.Pipelines.PipelineModules import HelloWorld
pipeline.append_module(HelloWorld.HelloWorld())

#add an analysis module with parameters
from ThackTech.Pipelines.PipelineModules import Sleep
sleeper = Sleep.Sleep()
sleeper.set_parameter('sleep_time', 10) #number of seconds to sleep for
pipeline.append_module(sleeper)

#initialize a runner
from ThackTech.Pipelines import ParallelPipelineRunner
runner = ParallelPipelineRunner(pipeline, threads=4)

#run the pipeline over some samples
samples = [PipelineSample('sample_'+i, 'mm9', '/path/to/destination') for i in range(5)]
runner.run(samples)
```

## Pipeline Documentation
Documenting pipeline and analysis parameters is of the upmost importance. Therefore, `AnalysisPipeline` and modules derived from `PipelineModule` offer a `documentation()` method that returns a string documentation of the pipeline or module in question. This is also useful for interogating module parameters in and interactive session. Consider the output of the method call `sleeper.documentation()` in the context of the script above:

```
Sleep
----------------------------------------
CRITICAL:    False
PROCESSORS:  1
DESCRIPTION: Take a nap.
----------------------------------------
PARAMETERS:
Name        Type      Value    Default  Nullable    Choices    Description
----------  ------  -------  ---------  ----------  ---------  -----------------------------------------
sleep_time  int          10         10  False       None       Amount of time, in seconds, to sleep for.

RESOLVERS:
        No Resolvers Declared

TOOL VERSIONS:
        No Tools Declared
----------------------------------------
```

Consider the output of the method call `pipelien.documentation()` in the context of the script above:

```
========================================
Begin Pipeline: Hello World
========================================

 |  |  |  |  |  |  |  |  |  |  |  |  |
 V  V  V  V  V  V  V  V  V  V  V  V  V

HelloWorld
----------------------------------------
CRITICAL:    False
PROCESSORS:  1
DESCRIPTION: Say "Hello World".
----------------------------------------
PARAMETERS:
        No Parameters Declared

RESOLVERS:
        No Resolvers Declared

TOOL VERSIONS:
        No Tools Declared
----------------------------------------

 |  |  |  |  |  |  |  |  |  |  |  |  |
 V  V  V  V  V  V  V  V  V  V  V  V  V

Sleep
----------------------------------------
CRITICAL:    False
PROCESSORS:  1
DESCRIPTION: Take a nap.
----------------------------------------
PARAMETERS:
Name        Type      Value    Default  Nullable    Choices    Description
----------  ------  -------  ---------  ----------  ---------  -----------------------------------------
sleep_time  int          10         10  False       None       Amount of time, in seconds, to sleep for.

RESOLVERS:
        No Resolvers Declared

TOOL VERSIONS:
        No Tools Declared
----------------------------------------

 |  |  |  |  |  |  |  |  |  |  |  |  |
 V  V  V  V  V  V  V  V  V  V  V  V  V

========================================
End Pipeline: Hello World
========================================
```








