# import os
# import sys
# import time
# import platform
# import traceback
# import re
# #from multiprocessing import Pool, Manager 
# #from multiprocessing import Manager
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
# import subprocess


GLOBAL_MANAGER = mp.Manager()

from ThackTech import Common

from ThackTech.Pipelines.GenomeInfo import GenomeInfo
from ThackTech.Pipelines.FileInfo import FileInfo
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
			'PipelineSample',
			'ModuleParameter',
			'PipelineModule',
			'PipelineRunner',
			'SerialPipelineRunner',
			'ParallelPipelineRunner',
			'SlurmPipelineRunner'
		]
								
