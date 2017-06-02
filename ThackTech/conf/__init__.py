import os
import sys
try:
    import configparser
except:
    import ConfigParser as configparser



__config_search_paths = [
    os.path.dirname(__file__), #built in default config (conf directory within package)
    "/etc/thacktech", #system-wide config
    os.path.expanduser("~/.config/thacktech"), #local per-user config
    os.environ.get("THACKTECH_CONF") #dynamic settable config dir
]

try:
    #calling os.getcwd() sometimes fails with an exception
    #this was observed in an interactive python session
    #check the current working directory
    __config_search_paths.append(os.getcwd())
except:
    pass



__config_locations = {}
def register_config_location(conf_type, path):
    global __config_locations
    if conf_type not in __config_locations:
        __config_locations[conf_type] = []
        __register_default_locations(conf_type)
    
    paths_to_add = []
    if os.path.exists(path) and os.path.isdir(path):
        paths_to_add.append(os.path.join(path, conf_type+".default.ini"))
        paths_to_add.append(os.path.join(path, conf_type+".ini"))
    else:
        paths_to_add.append(path)

    for p in paths_to_add:
        if p not in __config_locations[conf_type]:
            __config_locations[conf_type].append(p)
#end register_config_location()

def __register_default_locations(conf_type):
    global __config_search_paths
    for loc in __config_search_paths:
        if loc is not None:
            register_config_location(conf_type, loc)
#end __register_default_locations()    

__known_configs = {}
def get_config(conf_type, conf_files=[], force_reparse=False, **kwargs):
    """Gets a ConfigParser object for the given configuration type
    
    The results of this operation are cached by default. One may force a reparse 
    by passing True to the force_reparse argument.
    
    If a list of file paths is passed to conf_files then these are appended to the
    list of config files collected from the default/built-in locations. In other words,
    configuration file order determines precedence of the declared key/value pairs such 
    that config files later in order have higher precedence (can override) than config
    files earlier in order.
    
    """
    global __known_configs
    global __config_locations
    if conf_type not in __known_configs or force_reparse:
        if conf_type not in __config_locations:
            __register_default_locations(conf_type)
        
        paths = __config_locations[conf_type] + conf_files
        
        config = configparser.SafeConfigParser(**kwargs)
        config.read(paths)
        __known_configs[conf_type] = config
        
    return __known_configs[conf_type]
#end get_config


    