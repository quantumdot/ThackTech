import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class MergeBams(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MergeBams', short_description='Merge BAM files')
		super_args.update(**kwargs)
		super(MergeBams, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('samtools_path', str, 'samtools', desc="Path to samtools"))
		self.add_parameter(ModuleParameter('postfix', str, None, nullable=True, desc="Postfix to apply to merged bam filename"))
		self.add_parameter(ModuleParameter('sort', bool, True, desc="Sort the resulting merged BAM"))
		self.add_parameter(ModuleParameter('index', bool, True, desc="Index the resulting merged BAM. If sort is enabled, indexing occurs after sorting."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
	#end __declare_resolvers()

	
	def tool_versions(self):
		return {
			'samtools': self._call_output(self.get_parameter_value('samtools_path')+" 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		
		inbams = self.resolve_input('alignments', cxt)
		if inbams is None or (isinstance(inbams, list) and len(inbams) <= 0):
			cxt.log.write("No bam files provided, exiting module\n")
			return
			
		outdir = os.path.join(cxt.sample.dest, 'merged_bam')
		filetools.ensure_dir(outdir)
		outbase = os.path.join(outdir, cxt.sample.name)
		if self.get_parameter_value('postfix') is not None:
			outbase += "."+self.get_parameter_value_as_string('postfix')
		outfile = outbase+'.merged.bam'
		
		samtools_cmd = [ self.get_parameter_value('samtools_path'), 'merge', outfile ] + inbams
		
		cxt.log.write("-> Merging bam files with samtools......\n")
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(samtools_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(samtools_cmd, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		cxt.sample.add_file(FileInfo(outfile, FileContext.from_module_context(cxt, 'merged_bam')))
		
		if self.get_parameter_value('sort'):
			from ThackTech.Pipelines.PipelineModules import SortBam
			sorter = SortBam.SortBam(processors=self.processors)
			sorter.set_parameter('samtools_path', self.get_parameter_value('samtools_path'))
			sorter.set_parameter('overwrite', True)
			sorter.set_resolver('alignments', lambda c: c.sample.find_files(lambda f: f.fullpath == outfile))
			sorter.run(cxt)
			
		if self.get_parameter_value('index'):
			from ThackTech.Pipelines.PipelineModules import IndexBam
			indexer = IndexBam.IndexBam(processors=self.processors)
			indexer.set_parameter('samtools_path', self.get_parameter_value('samtools_path'))
			indexer.set_resolver('alignments', lambda c: c.sample.find_files(lambda f: f.fullpath == outfile))
			indexer.run(cxt)
		
		#return {
		#	'merged_bam': outfile
		#}
	#end run()
#end class GeneratePseudoreplicates