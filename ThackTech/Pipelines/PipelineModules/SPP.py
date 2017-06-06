import os
import sys
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.ModuleParameter import ModuleParameter


class SPP(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='SPP', short_description='Cross-correlation analysis using SPP')
		super_args.update(**kwargs)
		super(SPP, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('spp_path', str, 'run_spp.R', desc="Path to the run_spp.R script."))
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bam')
		pass
	#end __declare_resolvers()
	
	def run(self, cxt):
		spp_dir = os.path.join(cxt.sample.dest, 'spp')
		filetools.ensure_dir(spp_dir)
		spp_results = os.path.join(spp_dir, 'spp_results.tsv')
		sys.stdout.flush()
		
		bam = self.resolve_input('bam', cxt)
		
		#compute cross-correlation scores
		if cxt.sample.get_attribute('PE'):
			first_mate_bam = self.extract_first_mate(bam.fullpath)
		else:
			first_mate_bam = bam.fullpath
			
		spp_args = [
			'Rscript', self.get_parameter_value('spp_path'), 
			'-savp',
			'-savd',
			'-rf',
			'-s=0:1:400',
			'-c={}'.format(first_mate_bam), 
			'-out={}'.format(spp_results),
			'-odir={}'.format(spp_dir),
			'-p={}'.format(self.processors)
		]
		sys.stdout.write("\t-> Performing Quality and cross-correlation analysis using SPP......")
		sys.stdout.write("\n..............................................\n")
		sys.stdout.write(" ".join(spp_args))
		sys.stdout.write("\n..............................................\n")
		sys.stdout.flush()
		self._run_subprocess(spp_args)
		sys.stdout.write('\t-> Completed Cross-correlation analysis...\n')
		sys.stdout.flush()
		
		if cxt.sample.get_attribute('PE'):
			os.remove(first_mate_bam)
			os.remove(first_mate_bam+'.bai')
		
		return {
		
		}
	#end run()
	
	def extract_first_mate(self, bam):
		output = os.path.join(dest, os.path.splitext(os.path.basename(bam))[0]+'.firstmates.sam')
		cmd = 'samtools view -h -f 0x0040 %s > %s' % (bam, output)
		sys.stdout.write("\t-> Extracting first mates......")
		sys.stdout.write("\n..............................................\n")
		sys.stdout.write(cmd)
		sys.stdout.write("\n..............................................\n")
		sys.stdout.flush()
		self._run_subprocess(cmd, shell=True)
		return sam_to_bam(output, self.processors)
	#end extract_first_mate()
#end class SPP