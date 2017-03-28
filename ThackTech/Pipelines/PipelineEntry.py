#!/usr/bin/env python

import os
import argparse
import threading
import dill
import time
from ThackTech.Pipelines.PipelineRunner import _execute_pipeline_on_sample



if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('pipeline', help="pipeline pickle file")
	parser.add_argument('sample', help="sample pickle file")
	parser.add_argument('status', help="sample status pickle file")
	args = parser.parse_args()
	
	if not os.path.exists(args.pipeline):
		raise ValueError("Path %s not found!" % (args.pipeline,))
	if not os.path.exists(args.cxt.sample):
		raise ValueError("Path %s not found!" % (args.cxt.sample,))
	if not os.path.exists(args.status):
		raise ValueError("Path %s not found!" % (args.status,))
	
	try:
		with open(args.pipeline, 'rb') as pf:
			pipeline_pickle = dill.load(pf)
		with open(args.cxt.sample, 'rb') as sf:
			sample_pickle = dill.load(sf)
		with open(args.status, 'rb') as ssf:
			status_pickle = dill.load(ssf)
		
		tasks_statuses = dict()
		tasks_statuses[sample_pickle.name] = status_pickle
		
		thread = threading.Thread(target=_execute_pipeline_on_sample, args=(pipeline_pickle, sample_pickle, tasks_statuses))
		thread.start()
		
		while thread.is_alive():
			time.sleep(0.5)
			with open(args.status, 'wb', 0) as f:
				dill.dump(tasks_statuses[sample_pickle.name], f)
		thread.join()
	except Exception as e:
		print e
		raise
