import os
import sys
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class SPP(PipelineModule):
	
	def __init__(self, **kwargs):
		super(SPP, self).__init__('SPP', 'Cross-correlation analysis using SPP', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		spp_dir = os.path.join(cxt.sample.dest, 'spp')
		filetools.ensure_dir(spp_dir)
		spp_results = os.path.join(spp_dir, 'spp_results.tsv')
		sys.stdout.flush()
		
		#compute cross-correlation scores
		if cxt.sample.get_attribute('PE'):
			first_mate_bam = self.extract_first_mate()
			os.remove(first_mate_bam)
			os.remove(first_mate_bam+'.bai')
		else:
			first_mate_bam = cxt.sample.get_file()
			
		sys.stdout.write('\t-> Completed Cross-correlation analysis...\n')
		sys.stdout.flush()
		out_dir = os.path.dirname(outfile)
		spp_args = [
			'Rscript', '/home/josh/scripts/phantompeakqualtools/run_spp.R', 
			'-savp',
			'-savd',
			'-rf',
			'-s=0:1:400',
			('-c=%s' % (bam,)), 
			('-out=%s' % (outfile,)),
			('-odir=%s' % (out_dir,)),
			('-p=%d' % (slef.processors,))
		]
		sys.stdout.write("\t-> Performing Quality and cross-correlation analysis using SPP......")
		sys.stdout.write("\n..............................................\n")
		sys.stdout.write(" ".join(spp_args))
		sys.stdout.write("\n..............................................\n")
		sys.stdout.flush()
		self._run_subprocess(spp_args)
		
		return {
		
		}
	#end run()
	
	def extract_first_mate(self):
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