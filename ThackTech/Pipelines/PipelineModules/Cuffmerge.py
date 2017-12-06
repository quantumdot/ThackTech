import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class Cuffmerge(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cuffmerge', short_description='Merge Transcript assemblies with Cuffmerge')
		super_args.update(**kwargs)
		super(Cuffmerge, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		cm_path = self.get_parameter_value_as_string('cuffmerge_path')
		return {
			'cuffmerge': self._call_output(cm_path+" --version 2>&1 | perl -ne 'if(m/.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('cuffmerge_path', str, 'cuffmerge', desc="Path to the cuffmerge executable."))
		
		self.add_parameter(ModuleParameter('min_isoform_fraction', float, 0.05, desc="Discard isoforms with abundance below this."))
		self.add_parameter(ModuleParameter('keep_temp', bool, False, desc="Keep all intermediate files during merge."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('assemblies')
		self._name_resolver('guide_gff')
	#end __declare_resolvers()
	
	def run(self, cxt):
		cm_args = [
			self.get_parameter_value_as_string('cuffmerge_path'),
			'-o', cxt.sample.dest,
			'--num-threads', str(self.processors),
			'--min-isoform-fraction', self.get_parameter_value_as_string('min_isoform_fraction')
			
			# @todo: implement these options:
			#-g/--ref-gtf
			#-s/--ref-sequence
		]
		if self.get_parameter_value('keep_temp'):
			cm_args.append('--keep-tmp')
		
		#add the input transcript assemblies that we are merging
		assemblies = self.resolve_input('assemblies', cxt)
		assembly_list_path = os.path.join(cxt.sample.dest, 'assembly_GTF_list.txt')
		cm_args.append(assembly_list_path)
		with open(assembly_list_path, 'w') as assembly_list:
			for a in assemblies:
				assembly_list.write('{}\n'.format(a.fullpath))
		
		
		
		#OK, we now have all the arguments setup, lets actually run cuffmerge
		cxt.log.write("\t->  Merging transcript assemblies with cuffmerge......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(cm_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(cm_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return [
			FileInfo(os.path.join(cxt.sample.dest, 'merged.gtf'), FileContext.from_module_context(cxt, 'merged_transcript_assembly'))
		]
		
	#end run()
#end class HelloWorld