import os
import sys
try:
    import configparser
except:
    import ConfigParser as configparser



__config_dirs_to_search = [
    os.path.dirname(__file__),
    os.path.expanduser("~/.config/thacktech"),
    "/etc/thacktech",
    os.environ.get("THACKTECH_CONF")
]


def get_config(name):
    config = configparser.SafeConfigParser()
    paths = []
    for loc in __config_dirs_to_search:
        if loc is not None:
            paths.append(os.path.join(loc, name+".default.ini"))
            paths.append(os.path.join(loc, name+".ini"))
    
    config.read(paths)    
    return config
#end get_config
        
        

    