import glob
import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ConvertBedgraphToBigWig(PipelineModule):

	def __init__(self):
		PipelineModule.__init__(self, 'BedgraphToBigWig', 'Converting Bedgraph To BigWig')
		self.add_parameter(ModuleParameter('bgtobw_path', str,	'bedGraphToBigWig'))
		self.add_parameter(ModuleParameter('bedtools_path', str,	'bedtools'))
		self.add_parameter(ModuleParameter('bedclip_path', str,	'bedClip'))
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	
	def run(self, cxt):
		
		bdgs = {}#glob.glob(os.path.join(cxt.sample.dest, '*', '*', '*.bdg.gz'))
		if cxt.sample.has_file_group('MACS2'):
			bdgs['treatment_signal'] = cxt.sample.get_file('MACS2', 'treatment_signal')
			if cxt.sample.has_file('MACS2', 'control_signal'):
				bdgs['control_signal'] = cxt.sample.get_file('MACS2', 'control_signal')
				
		elif cxt.sample.has_file_group('MACS1'):
			bdgs['treatment_signal'] = cxt.sample.get_file('MACS1', 'treatment_signal')
			if cxt.sample.has_file('MACS1', 'control_signal'):
				bdgs['control_signal'] = cxt.sample.get_file('MACS1', 'control_signal')
		
		
		cxt.log.write('Found %d bedgraph files for conversion.\n' % (len(bdgs),))
		cxt.log.flush()
		procs = []
		output_files = {}
		for label, bdg in bdgs.iteritems():
			if bdg.endswith('.gz'):
				bdg_nogz = os.path.splitext(bdg)[0]
			else:
				bdg_nogz = bdg
			if os.path.splitext(bdg_nogz)[1].lower() not in ['.bdg', '.bedgraph']:
				cxt.log.write('Skipping file "%s" as it does not appear to be a bedgraph file!\n' % (bdg,))
				cxt.log.flush()
				continue
			bw = os.path.splitext(bdg_nogz)[0]+'.bw'
			output_files[label] = bw
			#roundabout way to convert to BW as MACS can produce bdg entries that fall outside of valid chromosome ranges.
			#commands base on the script by Tao, but done (mostly) through pipes to eliminate intermediate files.
			#The bedgraphToBigwig program cannot take stdin however, as the program does two passes over the input file.
			#This is done to reduce memory usage, however if memory is not a concern, tham it is possible to use the 
			#wigToBigwig utility on the clipped bedgraph through stdin, but this uses a lot of memory!
			#see: https://gist.github.com/taoliu/2469050 
			cmd  = 'zcat -f "'+bdg+'" '																					#read compressed bedgraph
			cmd += '| '+self.get_parameter_value('bedtools_path')+' slop -i /dev/stdin -g "'+cxt.sample.genome.chrsize+'" -b 0 '	#convert bedgraph to bed
			cmd += '| '+self.get_parameter_value('bedclip_path')+' /dev/stdin "'+cxt.sample.genome.chrsize+'" "'+bdg_nogz+'.clip"'#remove bed entries extending outside valid coordinates
			cmd += '; ' #critical that we pause until the clipped version is complete
			cmd += self.get_parameter_value('bgtobw_path')+' "'+bdg_nogz+'.clip" "'+cxt.sample.genome.chrsize+'" "'+bw+'"'		#finally, convert to bigwig
			cmd += '; rm -f "'+bdg_nogz+'.clip"' #make sure we cleanup our mess!

			cxt.log.write('Converting Bedgraph to BigWig for %s....' % (bdg,))
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