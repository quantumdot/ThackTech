import os
import sys
from dateutil.relativedelta import relativedelta
import subprocess
import shlex
import pysam

def human_time_diff(start, end):
	'''returns a human readable string describing the elapsed time between start and end using units familiar to humans (i.e. months, days, hours, minutes... etc.) '''
	attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']
	delta = relativedelta(seconds=(end-start))
	return ', '.join(['%d %s' % (getattr(delta, attr), getattr(delta, attr) > 1 and attr or attr[:-1]) for attr in attrs if getattr(delta, attr)])
#end human_time_diff()

def ensure_dir(path):
	''' Creats the directories specified by path if they do not already exist '''
	if not os.path.exists(path):
		os.makedirs(path)
	return path
#end ensure_dir()

def which(name, all=False):
	results = []
	for path in os.getenv("PATH").split(os.path.pathsep):
		full_path = os.path.join(path, name)
		if os.path.exists(full_path):
			results.append(full_path)
			if not all:
				break
	return results
#end which()

def prepend_file(filename, data):
	'''Prepends the file specified by filename with data supplied by data '''
	with open(filename, 'r') as original:
		origional_data = original.read()
	with open(filename, 'w') as modified:
		modified.write(data + origional_data)
#end prepend_file()

def basename_noext(path, complete=False):
	'''
		Returns the path base name. If complete is True, then try to remove
		all extentions in the basename (until no '.' are remaining in the string
	'''
	base = os.path.splitext(os.path.basename(path))[0]
	if not complete:
		return base
	else:
		while True:
			new_base = os.path.splitext(os.path.basename(base))[0]
			if new_base == base:
				return new_base
			base = new_base
#end basename_noext()

def run_pipe(steps, outfile=None):
	#break this out into a recursive function
	#TODO:  capture stderr
	p = None
	p_next = None
	first_step_n = 1
	last_step_n = len(steps)
	for n,step in enumerate(steps, start=first_step_n):
		print "step %d: %s" %(n,step)
		if n == first_step_n:
			if n == last_step_n and outfile: #one-step pipeline with outfile
				with open(outfile, 'w') as fh:
					print "one step shlex: %s to file: %s" %(shlex.split(step), outfile)
					p = subprocess.Popen(shlex.split(step), stdout=fh)
				break
			print "first step shlex to stdout: %s" %(shlex.split(step))
			p = subprocess.Popen(shlex.split(step), stdout=subprocess.PIPE)
			#need to close p.stdout here?
		elif n == last_step_n and outfile: #only treat the last step specially if you're sending stdout to a file
			with open(outfile, 'w') as fh:
				print "last step shlex: %s to file: %s" %(shlex.split(step), outfile)
				p_last = subprocess.Popen(shlex.split(step), stdin=p.stdout, stdout=fh)
				p.stdout.close()
				p = p_last
		else: #handles intermediate steps and, in the case of a pipe to stdout, the last step
			print "intermediate step %d shlex to stdout: %s" %(n,shlex.split(step))
			p_next = subprocess.Popen(shlex.split(step), stdin=p.stdout, stdout=subprocess.PIPE)
			p.stdout.close()
			p = p_next
	out,err = p.communicate()
	return out,err
#end run_pipe()

def get_known_compression_ext():
	''' Return a list of known and supported compressed file extensions. Used by is_compressed() and extract() '''
	return ['.tar.bz2', '.tar.gz', '.tar.xz', '.bz2', '.gz', '.tar', '.tbz2', '.tgz', '.zip']
#end get_known_compression_ext()

def is_compressed(file):
	'''	Uses a pretty naive method of checking file extension to see if a file
		*looks* like it is compressed.'''
	for ext in get_known_compression_ext():
		if file.endswith(ext):
			return True
#end is_compressed

def extract(file, destdir, keeporigional=True, overwrite=False, wait=True):
	'''Extracts the file given by file to destination directory. Attempts to automagically detect the
	   compression format used. Returns the new file name and the subprocess process spawned.'''
	basename = os.path.basename(file)
	destfilename = basename
	extdepth = 0
	cmd = []
	
	if file.endswith('.tar.bz2'):
		extdepth = 2
		cmd = ['tar', 'xOjf']
	elif file.endswith('.tar.gz'):
		extdepth = 2
		cmd = ['tar', 'xOzf']
	elif file.endswith('.tar.xz'):
		extdepth = 2
		cmd = ['tar', 'xOJf']
	elif file.endswith('.bz2'):
		extdepth = 1
		cmd = ['bunzip2', '--decompress', '--stdout']
		if keeporigional:
			cmd.append('--keep')
	elif file.endswith('.gz'):
		extdepth = 1
		cmd = ['gunzip', '--stdout']
		#if keeporigional: #stdout should keep origional file
		#	cmd.append('--keep')
	elif file.endswith('.tar'):
		extdepth = 1
		cmd = ['tar', 'xOf']
	elif file.endswith('.tbz2'):
		extdepth = 1
		cmd = ['tar', 'xOjf']
	elif file.endswith('.tgz'):
		extdepth = 1
		cmd = ['tar', 'xOzf']
	elif file.endswith('.zip'):
		extdepth = 1
		cmd = ['unzip', '-p']
	#elif file.endswith('lzma'):
	#	extdepth = 1
	#	cmd = unlzma ./"$1"
	#elif file.endswith('rar'):
	#	extdepth = 1
	#	cmd = unrar x -ad ./"$1"
	#elif file.endswith('Z'):
	#	extdepth = 1
	#	cmd = uncompress ./"$1"
	#elif file.endswith('7z'):
	#	extdepth = 1
	#	cmd = 7z x ./"$1"
	#elif file.endswith('xz'):
	#	extdepth = 1
	#	#cmd = unxz ./"$1"
	#elif file.endswith('exe'):
	#	extdepth = 1
	#	#cmd = cabextract ./"$1"
	else:
		raise IOError('extract: "%s" - unknown archive method' % (file,))
	
	#remove the appropriate number of file extensions....
	for i in range(extdepth):
		destfilename = os.path.splitext(destfilename)[0]
	
	destination = os.path.join(destdir, destfilename)
	if os.path.exists(destination) and not overwrite:
		#we were told to not overwrite, and the file exists already!
		p = subprocess.Popen(["echo"])
		if wait:
			p.communicate()
	else:
		with open(destination, 'wb') as out:
			#sys.stderr.write(" ".join(cmd + [file])+"\n") 
			p = subprocess.Popen(cmd + [file], stdout=out)
			if wait:
				p.communicate()
	return (destination, p)
#end extract()



def is_bam_PE(bam_path):
	'''
		Find is a sam or bam file contains even one paired-end read
	'''
	with pysam.AlignmentFile(bam_path, "rb") as samfile:
		for read in samfile.fetch():
			if read.is_paired:
				return True
	return False
#end is_bam_PE()
