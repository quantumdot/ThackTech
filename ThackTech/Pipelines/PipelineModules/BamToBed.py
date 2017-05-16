import os
import subprocess
import gzip
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class BamToBed(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='BamToBed', short_description='Convert BAM to BED')
		super_args.update(**kwargs)
		super(BamToBed, self).__init__(**super_args)
		
		self.add_parameter(ModuleParameter('gzip', bool, True, desc="Gzip output"))
		
		self._name_resolver('bam')
	#end __init__()
	
	def tool_versions(self):
		return {
			'bamToBed': self._call_output("bamToBed -h 2>&1 | perl -ne 'if(m/Version: v([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def load_modules(self, cxt):
		pass
	#end load_modules()
	
	def run(self, cxt):
		cxt.log.write("\t-> Converting BAM to BED...\n")
		
		inputfile = self.resolve_input('bam', cxt)
		
		cmd = [
			'bamToBed',
			'-i', inputfile.fullpath
		]
		if cxt.sample.get_attribute('PE'):
			cxt.log.write("\t\t-> Conversion in PE mode...\n")
			cmd.append('-bedpe')
			outfile = os.path.join(cxt.sample.dest, inputfile.basename_with_ext('bedpe'))
		else:
			outfile = os.path.join(cxt.sample.dest, inputfile.basename_with_ext('bed'))
		
		if self.get_parameter_value('gzip'):
			cxt.log.write("\t\t-> Results will be gzipped...\n")
			outhandle = gzip.open(outfile+'.gz', 'wb')
		else:
			outhandle = open(outfile, 'w')
			
		p = subprocess.Popen(cmd, stdout=outhandle, stderr=cxt.log)
		p.communicate()
		outhandle.close()
		
		outbam = FileInfo(outfile, FileContext.from_module_context(cxt, "bed_from_bam"))
		return [outbam]
	#end run()
#end class SamToBam