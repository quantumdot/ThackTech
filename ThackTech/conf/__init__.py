import os
import sys
try:
    import configparser
except:
    import ConfigParser as configparser



__config_dirs_to_search = [
    os.path.dirname(__file__), #built in default config (conf directory within package)
    "/etc/thacktech", #system-wide config
    os.path.expanduser("~/.config/thacktech"), #local per-user config
    os.environ.get("THACKTECH_CONF") #dynamic settable config dir
]


try:
    #calling os.getcwd() sometimes fails with an exception
    #this was observed in an interactive python session
    #check the current working directory
    __config_dirs_to_search.append(os.getcwd())
except:
    pass


__known_configs = {}

def get_config(name, force_reparse=False, **kwargs):
    if name not in __known_configs or force_reparse:
        config = configparser.SafeConfigParser(**kwargs)
        paths = []
        for loc in __config_dirs_to_search:
            #sys.stderr.write("Looking for config files in {}\n".format(loc))
            if loc is not None:
                paths.append(os.path.join(loc, name+".default.ini"))
                paths.append(os.path.join(loc, name+".ini"))
        
        config.read(paths)
        __known_configs[name] = config
        
    return __known_configs[name]
#end get_config


    