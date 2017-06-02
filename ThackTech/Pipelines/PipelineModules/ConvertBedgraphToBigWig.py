import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ConvertBedgraphToBigWig(PipelineModule):

	def __init__(self, **kwargs):
		super_args = dict(name='BedgraphToBigWig', short_description='Converting Bedgraph To BigWig')
		super_args.update(**kwargs)
		super(ConvertBedgraphToBigWig, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('bgtobw_path', str, 'bedGraphToBigWig'))
		self.add_parameter(ModuleParameter('bedtools_path', str, 'bedtools'))
		self.add_parameter(ModuleParameter('bedclip_path', str, 'bedClip'))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bedgraphs')
	#end __declare_resolvers()
	
	
	def run(self, cxt):
		bdgs = self.resolve_input('bedgraphs', cxt)
		
		cxt.log.write('Found {} bedgraph files for conversion.\n'.format(len(bdgs)))
		cxt.log.flush()
		procs = []
		output_files = {}
		for bdg in bdgs:
			if bdg.fullpath.endswith('.gz'):
				bdg_nogz = os.path.splitext(bdg.fullpath)[0]
			else:
				bdg_nogz = bdg.fullpath
			if os.path.splitext(bdg_nogz)[1].lower() not in ['.bdg', '.bedgraph']:
				cxt.log.write('Skipping file "{}" as it does not appear to be a bedgraph file!\n'.format(bdg.fullpath))
				cxt.log.flush()
				continue
			bw = os.path.splitext(bdg_nogz)[0]+'.bw'
			output_files[bdg.cxt.role] = bw
			#roundabout way to convert to BW as MACS can produce bdg entries that fall outside of valid chromosome ranges.
			#commands base on the script by Tao, but done (mostly) through pipes to eliminate intermediate files.
			#The bedgraphToBigwig program cannot take stdin however, as the program does two passes over the input file.
			#This is done to reduce memory usage, however if memory is not a concern, tham it is possible to use the 
			#wigToBigwig utility on the clipped bedgraph through stdin, but this uses a lot of memory!
			#see: https://gist.github.com/taoliu/2469050 
			cmd  = 'zcat -f "{in_bdg}" '										#read compressed bedgraph
			cmd += '| {bedtools_path} slop -i /dev/stdin -g "{chrsize}" -b 0 '	#convert bedgraph to bed
			cmd += '| {bedclip_path} /dev/stdin "{chrsize}" "{tmp_bdg}"'		#remove bed entries extending outside valid coordinates
			cmd += '; ' 														#critical that we pause until the clipped version is complete
			cmd += '{bgtobw_path} "{tmp_bdg}" "{chrsize}" "{out_bw}"'			#finally, convert to bigwig
			cmd += '; ' 														#critical that we pause until the conversion to bigwig is complete
			cmd += 'rm -f "{tmp_bdg}"' 											#make sure we cleanup our mess!
			
			cmd = cmd.format(bedtools_path=self.get_parameter_value('bedtools_path'),
							 bedclip_path=self.get_parameter_value('bedclip_path'),
							 bgtobw_path=self.get_parameter_value('bgtobw_path'),
							 chrsize=cxt.sample.genome.chrsize,
							 tmp_bdg=bdg_nogz+'.clip',
							 out_bw=bw,
							 in_bdg=bdg.fullpath)

			cxt.log.write('Converting Bedgraph to BigWig for {}....'.format(bdg.fullpath))
			cxt.log.write("\n..............................................\n")
			cxt.log.write(cmd)
			cxt.log.write("\n..............................................\n\n")
			cxt.log.flush()
			procs.append(subprocess.Popen(cmd, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log, shell=True))
		cxt.log.flush()
		for p in procs:
			p.communicate()
			
		return output_files
	#end run()
#end class ConvertBedgraphToBigWig