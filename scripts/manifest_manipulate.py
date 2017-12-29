#!/usr/bin/env python

import os
import sys
import argparse
from ThackTech import filetools
from ThackTech.Pipelines import PipelineSample



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', nargs='+', help="Path to manifest(s) to manipulate.")
    action_choices = ['show', 'del']
    parser.add_argument('--action', action='store', choices=action_choices, default=action_choices[0], help="Action to perform")
    parser.add_argument('--nocommit', action='store_true', help="Do not commit any changes, just show what would be done.")
    
    filter_group = parser.add_argument_group('Filters')
    filter_group.add_argument('--pipeline')
    filter_group.add_argument('--step')
    filter_group.add_argument('--module')
    filter_group.add_argument('--role')
    filter_group.add_argument('--path')
    filter_group.add_argument('--attribute')
    args = parser.parse_args()
    
    changed = False
    for manifest in args.manifest:
        m_path = os.path.abspath(manifest)
        s = PipelineSample(filetools.basename_noext(m_path, True), 'mm9', os.path.dirname(m_path))
        s.read_file_manifest(m_path)
        
        filter_func = generate_filter(args)
        for f in s.files:
            if filter_func(f):
                if args.action == 'show':
                    print f
                elif args.action == 'del':
                    changed = True
                    sys.stderr.write('Removing file {} from manifest {}'.format(str(f), s.name))
                    s.remove_file(f)
                    
        if changed and not args.nocommit:
            sys.stderr.write('Writing out manifest for {}'.formt(m_path))
            s.write_file_manifest(m_path)
    
    
def generate_filter(args):
    
    def passes_filter(f):
        if args.pipeline is not None and f.cxt.pipeline not in args.pipeline:
            return False
        if args.step is not None and f.cxt.step not in args.step:
            return False
        if args.module is not None and f.cxt.module not in args.module:
            return False
        if args.role is not None and f.cxt.role not in args.role:
            return False
        if args.attribute is not None:
            for attr in args.attribute:
                if not f.has_attribute_value(*attr.split('=')):
                    return False
        return True
    #end passes_filter()
    return passes_filter
#end generate_filter()
    
    