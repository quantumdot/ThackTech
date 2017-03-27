import os
from ConfigParser import ConfigParser


__dirs_to_search = [
    
    
]

def get_config(name):
    
    config = ConfigParser()
    for d in __dirs_to_search:
        default_file = os.path.join(d, name+".default.ini")
        if os.path.exists(default_file):
            config.
    dbconf.readfp(open('default.ini'))
    if os.path.exists('environment.ini'):
        dbconf.readfp(open('environment.ini'))
    dbconf.get('database', 'server') # Returns 192.168.0.12