import os
import sys
try:
    import configparser
except:
    import ConfigParser as configparser



__config_dirs_to_search = [
    os.path.dirname(__file__), #built in default config
    os.path.expanduser("~/.config/thacktech"), #local per-user config
    "/etc/thacktech", #system-wide config
    os.environ.get("THACKTECH_CONF"), #dynamic settable config dir
    os.getcwd() #check the current working directory
]

__known_configs = {}

def get_config(name, force_reparse=False):
    if name not in __known_configs or force_reparse:
        config = configparser.SafeConfigParser()
        paths = []
        for loc in __config_dirs_to_search:
            sys.stderr.write("Looking for config files in {}\n".format(loc))
            if loc is not None:
                paths.append(os.path.join(loc, name+".default.ini"))
                paths.append(os.path.join(loc, name+".ini"))
        
        config.read(paths)
        __known_configs[name] = config
        
    return __known_configs[name]
#end get_config
        
        

    