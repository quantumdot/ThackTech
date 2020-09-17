import urllib
import urllib2
import urlparse
import xml.etree.ElementTree as ET





class BioDas(object):
    
    def __init__(self, site="http://genome.ucsc.edu/cgi-bin/das"):
        url_parts = urlparse.urlparse(site, 'http')
        self.scheme = url_parts.scheme
        self.host = url_parts.netloc
        self.path = url_parts.path
        self.dabase = None
        
    def select_database(self, database):
        self.dabase = database
    
    def request(self, database, command, **data):
        if database is None:
            url = "{scheme}://{host}{path}/{database}/{command}".format(scheme=self.scheme, host=self.host, path=self.path, database=database, command=command)
        else:
            url = "{scheme}://{host}{path}/{command}".format(scheme=self.scheme, host=self.host, path=self.path, command=command)
        url += "?"+urllib.urlencode(data, True)
        
        response = urllib2.urlopen(url)
        with open(cache_path, 'w+') as cf:
            cf.write(f.read())
        f.close()

    
    def parse_response(self, response, output='json'):
        supported_output_types = ['xml', 'json']
        if output.lower() not in supported_output_types:
            raise ValueError("{} is not supported as a return type!")
        
        root = ET.fromstring(response.read())
        
        return root
    
    def sources(self):
        response = self.request(None, 'dsn')
        results = []
        for dsn in response.findall("DSN"):
            results.append({
                'id': dsn.find('SOURCE').get('id'),
                'version': dsn.find('SOURCE').get('version'),
                'name': dsn.find('SOURCE').text,
                'description': dsn.find('DESCRIPTION').text,
                'url': dsn.find('MAPMASTER').text,
            })
        return results
    #end sources()
        
    def entry_points(self, database):
        response = self.request(database, 'entry_points')
        #.....do something with response.....
    
    def sequence(self, segments, database):   
        response = self.request(database, 'command')
        root.findtext("./SEQUENCE/DNA").replace('\n','').replace('\r','')
        
        
        
        
        
        
        
        
        
        
        
        
        
        