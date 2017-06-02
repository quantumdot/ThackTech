import os
from ThackTech.Pipelines import PipelineModule



class OutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='OutputManifest', short_description='Generate Output Manifest')
		super_args.update(**kwargs)
		super(OutputManifest, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()

	
	def run(self, cxt):
		cxt.log.write("-> Writing sample output manifest....\n")
		cxt.log.flush()
		dest = os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv')
		if os.path.exists(dest):
			write_headers = False
		else:
			write_headers = True
		
		with open(dest, 'w') as output_manifest:
			if write_headers:
				output_manifest.write('id\tparent\tpipeline\tstep\tmodule\trole\tpath\tattributes\n')
			
			i = 0
			for f in cxt.sample.files:
				self.write_output(i, "None", f, output_manifest)
				
				if len(f.companions) > 0:
					j=1
					for fc in f.companions:
						self.write_output(i+j, i, fc, output_manifest)
						j += 1
					i += j
				
				i += 1
	#end run()
	
	def write_output(self, index, parent, f, out):
		out.write("{id}\t{parent}\t{pipeline}\t{step}\t{module}\t{role}\t{path}\t{attributes}\n".format(id=index,
																							  parent=parent,
																							  pipeline=f.cxt.pipeline,
																							  step=f.cxt.step,
																							  module=f.cxt.module,
																							  role=f.cxt.role,
																							  path=f.fullpath,
																							  attributes=';'.join("%s=%r" % (key,val) for (key,val) in f.attributes.iteritems())))
#end class OutputManifest

