import os
import subprocess
import re
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class CEAS(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='CEAS', short_description='Cis-Regulatory Annotation System')
		super_args.update(**kwargs)
		super(CEAS, self).__init__(**super_args)
		
		self.add_parameter(ModuleParameter('ceas_wig_path', str, 'ceas'))
		self.add_parameter(ModuleParameter('ceas_bw_path', str, '/home/josh/scripts/cistrome/ceasBW'))
		self.add_parameter(ModuleParameter('additional_args', list, []))
		
		self._name_resolver('signal')
		self._name_resolver('peaks')
	#end __init__()
	
	def tool_versions(self):
		versions = {}
		tools = {
			'ceas_wig':		self.get_parameter_value('ceas_wig_path'),
			'ceas_bigwig':	self.get_parameter_value('ceas_bw_path')
		}
		for name, path in tools:
			vout = subprocess.check_output(path+' --version 2>&1 | tail -n 1', shell=True, stderr=subprocess.STDOUT)
			versions[name] = re.match(r"ceas[\s-]+([\d. \(\)\w]+)", vout).group(1)
		return versions
	#end tool_versions
	
	
	def show_version(self, handle=None, fancy=True, path=''):
		buff = ""
		if fancy:
			buff = 'CEAS version is....\n  -> '
		buff += subprocess.check_output(path+' --version 2>&1 | tail -n 1', shell=True, stderr=subprocess.STDOUT)
		if fancy:
			buff += '\n'
		if handle is not None:
			handle.write(buff)
			handle.flush()
		else:
			return buff
	#end show_version()
	
	def run(self, cxt):
		signal_file = self.resolve_input('signal', cxt)
		signal_compressed = False
		
		#check for compressed files
		if os.path.splitext(signal_file)[1] == '.gz':
			signal_compressed = True
			cxt.log.write("Decompressing signal file......\n")
			proc = subprocess.Popen(['gunzip', '--keep', signal_file], cwd=cxt.sample.dest)
			proc.communicate()
			signal_file = os.path.splitext(signal_file)[0]
		
		if os.path.splitext(signal_file)[1].lower() in ['.bw', '.bigwig']:
			ceas_exe = self.get_parameter_value('ceas_bw_path')
			self.get_parameter('additional_args').value + ['-l', cxt.sample.genome.chrsize]
		else:
			ceas_exe = self.get_parameter_value('ceas_wig_path')

		ceas_args = [
			ceas_exe,
			'-b', self.resolve_input('peaks', cxt),
			'-w',  signal_file,
			'-g', ('/mnt/ref/ceas/%s.refGene' % (cxt.sample.genome.name,)),
			'--name', cxt.sample.name
		] + self.get_parameter_value('additional_args')
			
		#Run CEAS
		cxt.log.write("Performing Cis-regulatory Element Annotation......\n")
		self.show_version(cxt.log, True, ceas_exe)
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(ceas_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()	

		self._run_subprocess(ceas_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		#generate the model figure
		fig_path = os.path.join(cxt.sample.dest, cxt.sample.name+'.R')
		if os.path.exists(fig_path):
			cxt.log.write('\t-> Generating PDF of results.....\n')
			cxt.log.flush()
			with open(fig_path, 'r') as rmodel:
				with open(os.devnull, 'w') as devnull:
					self._run_subprocess(['R', '--quiet', '--vanilla'], stdin=rmodel, stderr=subprocess.STDOUT, stdout=devnull, cwd=cxt.sample.dest)

		if signal_compressed:
			os.remove(signal_file)
		
		#return output_files
	#end run()
#end class ChanceAnalysis