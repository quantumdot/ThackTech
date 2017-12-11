import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class GffCompare(PipelineModule):
	'''
	STATUS: EXPERIMENTAL - Not Fully Implemented!
	'''
	
	def __init__(self, **kwargs):
		super_args = dict(name='GffCompare', short_description='Compare GFF files.')
		super_args.update(**kwargs)
		super(GffCompare, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		p = self.get_parameter_value_as_string('gffcompare_path')
		return {
			'gffcompare': self._call_output(p+" --version 2>&1 | perl -ne 'if(m/.*gffcompare.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('gffcompare_path', str, 'gffcompare', desc="Path to the gffcompare executable."))
		
		self.add_parameter(ModuleParameter('only_ref_trans', bool, False, desc="-R: consider only the reference transcripts that overlap any of the input transfrags (Sn correction)"))
		#self.add_parameter(ModuleParameter('', , , desc=""))
		#self.add_parameter(ModuleParameter('', , , desc=""))
		#self.add_parameter(ModuleParameter('', , , desc=""))
		#self.add_parameter(ModuleParameter('', , , desc=""))
		
		self.add_parameter(ModuleParameter('no_sample_output', bool, False, desc="-T: do not generate .tmap and .refmap files for each input file."))
		
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('reference_gtf')
		self._name_resolver('gtfs')
	#end __declare_resolvers()
	
	def run(self, cxt):		
		
		out_prefix = os.path.join(cxt.sample.dest, cxt.sample.name+'_gffcmp')
		out_files = []
		
		gtfs = self.resolve_input('gtfs', cxt)
		
		gffcmp_args = [
			self.get_parameter_value_as_string('gffcompare_path'),
			'-o', out_prefix,
			
			
		]
		
		if self.get_parameter_value('no_sample_output'):
			gffcmp_args.append('-T')
		else:
			for gtf in gtfs:
				gtf_base = gtf.filename_strip_all_ext
				out_files.append(FileInfo(gtf_base+'.refmap', FileContext.from_module_context(cxt, 'ref_trans_map')))
				out_files.append(FileInfo(gtf_base+'.tmap', FileContext.from_module_context(cxt, 'best_ref_trans')))
		
		#add the input GTFs as positional args
		gtfs = self.resolve_input('gtfs', cxt)
		gffcmp_args.extend([g.fullpath for g in gtfs])
		
		
		return [
			#always
			FileInfo(out_prefix+'.stats', FileContext.from_module_context(cxt, 'summary_stats')),
			FileInfo(out_prefix+'.tracking', FileContext.from_module_context(cxt, 'transfrags_tracking')),
			
			#depends on number of gtf files arguments
			FileInfo(out_prefix+'.combined.gtf', FileContext.from_module_context(cxt, 'transfrags_union')), #multiple
			#OR
			FileInfo(out_prefix+'.annotated.gtf', FileContext.from_module_context(cxt, 'annotated_gtf')), #single
			
			#for each gtf input:
			FileInfo('<gff_in>.refmap', FileContext.from_module_context(cxt, 'role')),
			FileInfo('<gff_in>.tmap', FileContext.from_module_context(cxt, 'role')),
			
			
		]
		
	#end run()

#end class GffCompare
