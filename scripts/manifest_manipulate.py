#!/usr/bin/env python

import os
import sys
import argparse
from ThackTech import filetools
from ThackTech.Pipelines import PipelineSample



def main():
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest='action')
    show_cmd = subparsers.add_parser('show', parents=[parser], help="Show files matching filters")
    del_cmd = subparsers.add_parser('del', parents=[parser], help="Remove files matching filters")
    
    move_cmd = subparsers.add_parser('move', parents=[parser], help="Move files matching filters")
    move_cmd.add_argument('dest', help="destination for files that match. Use string formatting tokens for variables. i.e. {sname}. Tokens: [sdest, sname]")
    move_cmd.add_argument('--fs', action='store_true', help="move the file on the filesystem in addition to changing the manifest entry.")

    parser.add_argument('manifest', nargs='+', help="Path to manifest(s) to manipulate.")
    #action_choices = ['show', 'del']
    #parser.add_argument('--action', action='store', choices=action_choices, default=action_choices[0], help="Action to perform")
    parser.add_argument('--nocommit', action='store_true', help="Do not commit any changes, just show what would be done.")
    
    filter_group = parser.add_argument_group('Filters')
    filter_group.add_argument('--pipeline', action='append')
    filter_group.add_argument('--step', action='append')
    filter_group.add_argument('--module', action='append')
    filter_group.add_argument('--role', action='append')
    filter_group.add_argument('--path', action='append')
    filter_group.add_argument('--attribute', action='append')
    args = parser.parse_args()
    
    filter_func = generate_filter(args)
    action_func = globals()["action_"+args.action]
        
    
    for manifest in args.manifest:
        m_path = os.path.abspath(manifest)
        m_dir = os.path.dirname(m_path)
        m_base = filetools.basename_noext(m_path, True)
        m_name = m_base.replace("_output_manifest", "")
        sys.stderr.write("Processing manifest {}\n".format(m_path))
        s = PipelineSample(m_name, 'mm9', m_dir)
        s.read_file_manifest(m_path)
        sys.stderr.write("Inferred Information:\n")
        sys.stderr.write("Sample Name: {}\n".format(s.name))
        sys.stderr.write("Sample Dest: {}\n".format(s.dest))
        sys.stderr.write("-----------------------------------------\n")
        
        match_count = 0
        change_count = 0
        for f in s.files[:]:
            if filter_func(f):
                match_count += 1
                action_func(s, f, args)
        
        if match_count <= 0:
            sys.stderr.write('No items matched.\n')
        sys.stderr.write("-----------------------------------------\n")
        if match_count > 0:
            sys.stderr.write('Matched {} items\n'.format(match_count))
            
        if change_count >= 0 and not args.nocommit:
            sys.stderr.write('Writing {} changes to output manifest {}\n'.format(change_count, m_path))
            s.write_file_manifest(m_path)
        elif change_count >= 0 and args.nocommit:
            sys.stderr.write('Running in --nocommit mode, {} changes will not be saved.\n'.format(change_count))
        sys.stderr.write('\n\n')
#end main()
 
#########################################
# BEGIN Action Plugins
#########################################

def action_show(s, f, args):
    print f
    return False
#end action_show()

def action_del(s, f, args):
    sys.stderr.write('Removing file {}\n'.format(str(f)))
    s.remove_file(f)
    return True
#end action_del()

def action_move(s, f, args):
    d = args.dest.format({'sdest': s.dest, 'sname': s.name})
    d_fullpath = os.path.join(d, f.basename)
    s.remove_file(f)
    if args.fs:
        sys.stderr.write('Moving file {} -> {}'.format(f.fullpath, d_fullpath))
        f.move(d)
    else:
        sys.stderr.write('Changing path {} -> {}'.format(f.fullpath, d_fullpath))
        f.__fullpath = d_fullpath
    s.add_file(f)
    return True
#end action_move
 
#########################################
# END Action Plugins
#########################################
    
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
    







if __name__ == "__main__":
    main()




