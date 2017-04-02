
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
# import subprocess


GLOBAL_MANAGER = mp.Manager()
try:
	CPU_COUNT = mp.cpu_count()
except:
	CPU_COUNT = 1
if CPU_COUNT is None:
	CPU_COUNT = 1

from ThackTech import Common

from ThackTech.Pipelines.GenomeInfo import GenomeInfo
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext
from ThackTech.Pipelines.Context import ModuleRunContext, BaseModuleContext
from ThackTech.Pipelines.PipelineSample import PipelineSample
from ThackTech.Pipelines.ModuleParameter import ModuleParameter
from ThackTech.Pipelines.PipelineModule import PipelineModule
from ThackTech.Pipelines.AnalysisPipeline import AnalysisPipeline
from ThackTech.Pipelines.PipelineRunner import PipelineRunner
from ThackTech.Pipelines.SerialPipelineRunner import SerialPipelineRunner
from ThackTech.Pipelines.ParallelPipelineRunner import ParallelPipelineRunner
from ThackTech.Pipelines.SlurmPipelineRunner import SlurmPipelineRunner


from ThackTech.Pipelines import PipelineModules

#from ThackTech.Processes import MultiStatusProgressBar, MultiStatusProgressItem

#from ThackTech.Pipelines import 
__all__ = [	'AnalysisPipeline',
			'GenomeInfo',
			'FileInfo',
			'FileContext',
			'PipelineSample',
			'ModuleParameter',
			'PipelineModule',
			'PipelineRunner',
			'SerialPipelineRunner',
			'ParallelPipelineRunner',
			'SlurmPipelineRunner',
			'ModuleRunContext',
			'BaseModuleContext',
			'GLOBAL_MANAGER'
		]
								
