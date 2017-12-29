#!/usr/bin/env python

import os
import sys
import argparse
from ThackTech import filetools
from ThackTech.Pipelines import PipelineSample



def main():
    oh = OptionHelper()
    args = oh.parse_args()
    
    filter_func = generate_filter(args)
    action_func = globals()["action_"+args.command]
        
    
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
                if action_func(s, f, args):
                    change_count +=1
                
        
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

def action_remove(s, f, args):
    sys.stderr.write('Removing file {}\n'.format(str(f)))
    s.remove_file(f)
    return True
#end action_del()

def action_move(s, f, args):
    d = args.dest.format(**{'sdest': s.dest, 'sname': s.name})
    d_fullpath = os.path.join(d, f.basename)
    s.remove_file(f)
    if args.fs:
        sys.stderr.write('Moving file {} -> {}\n'.format(f.fullpath, d_fullpath))
        f.move(d)
    else:
        sys.stderr.write('Changing path {} -> {}\n'.format(f.fullpath, d_fullpath))
        f._set_path(d_fullpath)
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
    

class OptionHelper(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Manipulate output manifest(s)',
            usage='''manifest_manipulate.py <command> [<args>]

The most commonly used commands are:
   show     Show files matching filters
   remove   Remove files matching filters
   move     Move files matching filters
''')
        self.parser.add_argument('command', help='Subcommand to run')
        
    def parse_args(self):
        args = self.parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print 'Unrecognized command'
            self.parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        cmd_parser = getattr(self, args.command)()
        self.add_filter_opts(cmd_parser)
        self.add_manifest_opts(cmd_parser)
        cmd_args = cmd_parser.parse_args(sys.argv[2:])
        cmd_args.command = args.command
        return cmd_args
        
    
    def add_filter_opts(self, parser):
        filter_group = parser.add_argument_group('Filters')
        filter_group.add_argument('--pipeline', action='append')
        filter_group.add_argument('--step', action='append')
        filter_group.add_argument('--module', action='append')
        filter_group.add_argument('--role', action='append')
        filter_group.add_argument('--path', action='append')
        filter_group.add_argument('--attribute', action='append')
        
    def add_manifest_opts(self, parser):
        parser.add_argument('manifest', nargs='+', help="Path to manifest(s) to manipulate.")
        parser.add_argument('--nocommit', action='store_true', help="Do not commit any changes, just show what would be done.")
    
    def show(self):
        parser = argparse.ArgumentParser(description='Show files matching filters')
        return parser

    def remove(self):
        parser = argparse.ArgumentParser(description='Show files matching filters')
        return parser
    
    def move(self):
        parser = argparse.ArgumentParser(description='Show files matching filters')
        parser.add_argument('dest', help="destination for files that match. Use string formatting tokens for variables. i.e. {sname}. Tokens: [sdest, sname]")
        parser.add_argument('--fs', action='store_true', help="move the file on the filesystem in addition to changing the manifest entry.")
        return parser





if __name__ == "__main__":
    main()




