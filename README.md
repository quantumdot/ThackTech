# ThackTech


# Introduction
The `ThackTech` package contains several modules and subpackages for genomic-related analysis. A centerpiece of the packages is the `ThackTech.Pipelines` package, which offers portable pipline systems and standardized analysis modules. By portability, the system attempts to be system agnostic (with the limitation of being mostly linux-bound), but can be run on anything from a laptop, workstation, or even a SLURM-managed cluster.

Importantly, the pipelines are declerative in nature and decouple data from the logic of execution. 

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
goldenpath=/mnt/reference/Mus_musculus/UCSC/mm9 #root of reference genome
```

Genome configuration may also be specified more manually:
```
[mm9]
size=1.91e9
chrsize=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chrom.sizes
index.BowtieIndex=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome
index.Bowtie2Index=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
index.BWAIndex=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome
fasta.Genome=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
fasta.chr1=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr1.fa
fasta.chr2=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr2.fa
#...
fasta.chrX=/mnt/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chrX.fa
```
One interesting note: Specifying an "goldenpath" will first parse the directory specified looking for valid indicies/fasta files/other stuff, but then later directives may override the auto-discovered locations which may be useful for specifying, for example, a bowtie index with higer density of the suffix-array sample (`--offrate`). 

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

# Authoring Analysis Modules
Creating an analysis module is pretty straightforward.
### Create a new class that extends `ThackTech.Pipelines.PipelineModule`
```python
from ThackTech.Pipelines import PipelineModule
class HelloWorld(PipelineModule):
	def __init__(self, **kwargs):
		super_args = dict(name='HelloWorld', short_description='Say "Hello World".')
		super_args.update(**kwargs)
		super(HelloWorld, self).__init__(**super_args)
```
Notes:
* It is important to call the super constructor in the way shown so that consumers of the module may override parameters such as the module name

### Implement some methods
At the minimum, an analysis module must implement the `run()` method:
```python
def run(self, cxt):
	cxt.log.write("Hello World from:\n%s\n" % (platform.uname(),))
	cxt.log.flush()
#end run()
```
The run method is passed a context object derived from `ThackTech.Pipelines.ModuleRunContext` and contains references to the sample, logging utilities, and information regarding the pipeline state. All computational work should occur within the `run()` method or methods directly called by `run()`

### Declare Module Parameters
In all but the very simplest analysis modules, it is useful to declare parameters that affect how the module runs. `PipelineModule`'s contain a collection of `ModuleParameter` objects that provide type-safe decleration of options for a given module. This parameters may also have values assigned at runtime from installed configuration files.

Module parameters should be declared by overriding the `PipelineModule._declare_parameters()` method. Within the method decleration, call the `PipelineModule.add_parameter()` method, passing a `ModuleParameter` object:
```python
def _declare_parameters(self):
	self.add_parameter(ModuleParameter('samtools_path', str, 'samtools', desc="Path to samtools"))
#end _declare_parameters()
```
Here we define a string type parameter named `samtools_path` with a default value of 'samtools'. We can get the value of this parameter at runtime by using the ``PipelineModule.get_parameter_value()` method.

Other important information about `ModuleParameters`
* Value is type-safe, so on setting `value`, the passed value is coerced to the type set for the parameter
* Parameters can be nullable (i.e. accept `None` values)
* Parameters can define valid values (besides type) using the `choices` attribute

### Declare Module Resolvers
Most modules will require some input that can only be determined at runtime, for example, a file that is created by a module earlier in the pipeline. To facilitate a declaritive approach that decouples data from the computation, Modules offer resolvers are essentially named callables.

Resolvers should be declared by overriding the `PipelineModule._declare_resolvers()` method, and place in the body calls to `PipelineModule._name_resolver()`. 
```python
def _declare_resolvers(self):
	self._name_resolver('alignments')
#end _declare_resolvers()
```

Within the `run()` method, one can retrieve the runtime result of the resolver by calling `PipelineModule.resolve_input()` with the resolver name and the current module run context.

### Other Notes about Modules
* Set the number of processors dedicated to this module using the `processors` attribute.
* Set the criticalness of the module with the `critical` boolean attribute. Modules marked as critical will cause the entire pipeline to fail if the module fails, otherwise the module fault will be marked as a warning and pipeline will continue to execute.
* The `critical` and `processors` attributes may be passed to the module constructor.
* It may be sometimes useful to override the module name by passing a `name` kwarg to the constructor to differentiate multiple of the same module within a given pipeline.

