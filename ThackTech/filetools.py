import os
import sys
import gzip
import errno
import tempfile
import subprocess




def make_str_filename_safe(unsafe_string):
    """Makes a string safe for use in filenames
    
    Removes any chars not alpha-numeric or in the set {-, _}, and replaces them with an empty string
    """
    keepcharacters = ('_','-')
    return "".join([c for c in unsafe_string if c.isalnum() or c in keepcharacters]).rstrip()
#end make_label_filename_safe()


def ensure_dir(path):
    """Creats the directories specified by path if they do not already exist.
    """
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as exception:
            #if the exception is raised because the directory already exits,
            #than our work is done and everything is OK, otherwise re-raise the error
            #THIS CAN OCCUR FROM A POSSIBLE RACE CONDITION!!!
            if exception.errno != errno.EEXIST:
                raise
    return path
#end ensure_dir()

def basename_noext(path, complete=False):
    '''Returns the path base name. If complete is True, then try to remove
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

def open_any(path, mode, **kwargs):
    if path.endswith('.gz'):
        return gzip.open(path, mode, **kwargs)
    else:
        return open(path, mode, **kwargs)
#end open_any()

def get_known_compression_ext():
    ''' Return a list of known and supported compressed file extensions. Used by is_compressed() and extract() '''
    return ['.tar.bz2', '.tar.gz', '.tar.xz', '.bz2', '.gz', '.tar', '.tbz2', '.tgz', '.zip']
#end get_known_compression_ext()

def is_compressed(filename):
    '''Uses a naive method of checking file extension to see if a file *looks* like it is compressed.'''
    for ext in get_known_compression_ext():
        if filename.endswith(ext):
            return True
#end is_compressed

def extract(filename, destdir, keeporigional=True, overwrite=False, wait=True, stderr=None):
    '''Extracts a file given by file to destination directory. Attempts to automagically detect the
       compression format used. Returns the new file name and the subprocess process spawned.'''
    basename = os.path.basename(filename)
    destfilename = basename
    extdepth = 0
    cmd = []
    
    if filename.endswith('.tar.bz2'):
        extdepth = 2
        cmd = ['tar', 'xOjf']
    elif filename.endswith('.tar.gz'):
        extdepth = 2
        cmd = ['tar', 'xOzf']
    elif filename.endswith('.tar.xz'):
        extdepth = 2
        cmd = ['tar', 'xOJf']
    elif filename.endswith('.bz2'):
        extdepth = 1
        cmd = ['bunzip2', '--decompress', '--stdout']
        if keeporigional:
            cmd.append('--keep')
    elif filename.endswith('.gz'):
        extdepth = 1
        cmd = ['gunzip', '--stdout']
        #if keeporigional: #stdout should keep origional file
        #    cmd.append('--keep')
    elif filename.endswith('.tar'):
        extdepth = 1
        cmd = ['tar', 'xOf']
    elif filename.endswith('.tbz2'):
        extdepth = 1
        cmd = ['tar', 'xOjf']
    elif filename.endswith('.tgz'):
        extdepth = 1
        cmd = ['tar', 'xOzf']
    elif filename.endswith('.zip'):
        extdepth = 1
        cmd = ['unzip', '-p']
    #elif filename.endswith('lzma'):
    #    extdepth = 1
    #    cmd = unlzma ./"$1"
    #elif filename.endswith('rar'):
    #    extdepth = 1
    #    cmd = unrar x -ad ./"$1"
    #elif filename.endswith('Z'):
    #    extdepth = 1
    #    cmd = uncompress ./"$1"
    #elif filename.endswith('7z'):
    #    extdepth = 1
    #    cmd = 7z x ./"$1"
    #elif filename.endswith('xz'):
    #    extdepth = 1
    #    #cmd = unxz ./"$1"
    #elif filename.endswith('exe'):
    #    extdepth = 1
    #    #cmd = cabextract ./"$1"
    else:
        raise IOError('extract: "%s" - unknown archive method' % (filename,))
    
    #remove the appropriate number of file extensions....
    for i in range(extdepth):
        destfilename = os.path.splitext(destfilename)[0]
    
    destination = os.path.join(destdir, destfilename)
    if os.path.exists(destination) and not overwrite:
        #we were told to not overwrite, and the file exists already!
        #spawn a dummy process so we can return it...
        p = subprocess.Popen(["echo"])
        if wait:
            p.communicate()
    else:
        with open(destination, 'wb') as out:
            #sys.stderr.write(" ".join(cmd + [file])+"\n") 
            p = subprocess.Popen(cmd + [filename], stdout=out, stderr=stderr)
            if wait:
                p.communicate()
    return (destination, p)
#end extract()


def prepend_file(filename, data):
    '''Prepends the file specified by filename with data supplied by data '''
    with open(filename, 'r') as original:
        origional_data = original.read()
    with open(filename, 'w') as modified:
        modified.write(data + origional_data)
#end prepend_file()


# class Tee(object):
#     """This object behaves like a file, but "tee-s" the data across multiple files-like objects
#     
#     """
#     def __init__(self, *args):
#         """Initialize this Tee object
#         
#         Parameters:
#             args: arbitrary number of file-like objects or strings. If string, it is assumed to be a file path and is opened in mode 'w'
#             
#         """
#         self.__innerhandle = tempfile.NamedTemporaryFile()
#         self.__handles = []
#         for arg in args:
#             if isinstance(arg, basestring):
#                 self.__handles.append(open(arg, 'w'))
#             else:
#                 self.__handles.append(arg)
#                 
#     def __poll(self):
#         pass
#                 
#     def release(self, filehandle):
#         self.__handles.remove(filehandle)
#         
#     def write(self, text):
#         """Write a string to the files held by this Tee
#         
#         Parameters:
#             text: String to write
#         """
#         for h in self.__handles:
#             h.write(text)
#             
#     def writelines(self, sequence):
#         """Write a sequence of strings to the files held by this Tee
#         
#         Parameters:
#             text: Sequence of strings to write
#         """
#         for text in sequence:
#             for h in self.__handles: 
#                 h.write(text)
# 
#     def flush(self):
#         """Flushes the output of all files held by this Tee
#         """
#         for h in self.__handles:
#             h.flush()
#     
#     def __del__(self):
#         for h in self.__handles:
#             if h != sys.stdout and h != sys.stderr:
#                 h.close()
# #end class Tee
