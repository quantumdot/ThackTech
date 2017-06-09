import os
import sys
from dateutil.relativedelta import relativedelta
import subprocess
import shlex
import collections, itertools

def human_time_diff(start, end):
	'''returns a human readable string describing the elapsed time between start and end using units familiar to humans (i.e. months, days, hours, minutes... etc.) '''
	attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']
	delta = relativedelta(seconds=(end-start))
	return ', '.join(['%d %s' % (getattr(delta, attr), getattr(delta, attr) > 1 and attr or attr[:-1]) for attr in attrs if getattr(delta, attr)])
#end human_time_diff()



def window(it, winsize, step=1):
	"""Sliding window iterator."""
	it=iter(it)  # Ensure we have an iterator
	l=collections.deque(itertools.islice(it, winsize))
	while 1:  # Continue till StopIteration gets raised.
		yield tuple(l)
		for i in range(step):
			l.append(next(it))
			l.popleft()


def which(name, all_results=False):
	results = []
	for path in os.getenv("PATH").split(os.path.pathsep):
		full_path = os.path.join(path, name)
		if os.path.exists(full_path):
			results.append(full_path)
			if not all_results:
				break
	return results
#end which()

def run_pipe(steps, outfile=None, stderr=None):
	#break this out into a recursive function
	#TODO:  capture stderr
	p = None
	p_next = None
	first_step_n = 1
	last_step_n = len(steps)
	for n,step in enumerate(steps, start=first_step_n):
		sys.stderr.write("step {}: {}".format(n, step))
		if n == first_step_n:
			if n == last_step_n and outfile: #one-step pipeline with outfile
				with open(outfile, 'w') as fh:
					sys.stderr.write("one step shlex: {} to file: {}".format(shlex.split(step), outfile))
					p = subprocess.Popen(shlex.split(step), stdout=fh, stderr=stderr)
				break
			sys.stderr.write("first step shlex to stdout: {}".format(shlex.split(step)))
			p = subprocess.Popen(shlex.split(step), stdout=subprocess.PIPE, stderr=stderr)
			#need to close p.stdout here?
		elif n == last_step_n and outfile: #only treat the last step specially if you're sending stdout to a file
			with open(outfile, 'w') as fh:
				sys.stderr.write("last step shlex: {} to file: {}".format(shlex.split(step), outfile))
				p_last = subprocess.Popen(shlex.split(step), stdin=p.stdout, stdout=fh, stderr=stderr)
				p.stdout.close()
				p = p_last
		else: #handles intermediate steps and, in the case of a pipe to stdout, the last step
			sys.stderr.write("intermediate step {} shlex to stdout: {}".format(n, shlex.split(step)))
			p_next = subprocess.Popen(shlex.split(step), stdin=p.stdout, stdout=subprocess.PIPE, stderr=stderr)
			p.stdout.close()
			p = p_next
	out,err = p.communicate()
	return out,err
#end run_pipe()

