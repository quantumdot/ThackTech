import os
import sys
import time
import subprocess
import pybedtools



class gfServer(object):
    
    def __init__(self, name, executible, host, port, genome2bit, **kwargs):
        """Initializes a new gfServer instance, but deos not start it.
        
        API:
            start() - starts the server, if it is not already running
            stop() - stops the server
            isPCR() - query the server for isPCR hits
        
        Parameters:
            name: (string) the name of this server, typically the genome assembly/build this server uses
            executable: (string) path to the gfServer executable. It is also assumed that the gfPcr executable is in the same directory as gfServer
            host: (string) Hostname to use in setting up the server
            port: (int) port number to assign to this server
            genome2bit: (string) path to a fasta file in 2bit format that you would like the server to index
            
        """
        self.name = name
        self.exe = executible
        self.host = host
        self.port = int(port)
        self.genome = genome2bit
        
        self.gfServe_params = {
			'tile_size': 11,
	        'step_size': 5,
	        'min_match': 2,
	        'max_gap': 2,
	        'translate': False,
	        'mask': False,
	        'rep_match': 1024,
	        'max_dna_hits': 100,
	        'max_trans_hits': 200,
	        'max_nt_size': 40000,
	        'max_as_size': 8000
		}
        self.gfServe_params.update(kwargs)
    #end __init__()
        
        
    def start(self):
        """Starts the gfServer
        
        """
        sys.stderr.write('Starting gfServer\n')
        if subprocess.call('{gfserver} status {gfhost} {gfport}'.format(gfserver=self.exe, gfhost=self.host, gfport=self.port), shell=True) == 0:
            sys.stderr.write('gfServer appears to be running. Skipping gfServer startup...\n\n\n')
            return
        
        cwd = os.getcwd()
        log_dir = os.path.join(cwd, 'gfserv')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'gfserver.{}.temp.log'.format(self.name))
        if os.path.exists(log_file):
            os.remove(log_file)   
        os.chdir(os.path.dirname(self.genome))
        try:
            time.sleep(5)
            cmd = '{gfserver} -canStop -log="{log}" '.format(gfserver=self.exe, log=log_file)
            cmd += '-tileSize={tile_size} -stepSize={step_size} -minMatch={min_match} -maxGap={max_gap} ' \
            	 + '-repMatch={rep_match} -maxDnaHits={max_dna_hits} -maxTransHits={max_trans_hits} -maxNtSize={max_nt_size} -maxAsSize={max_as_size} '.format(**self.gfServe_params)
            
            if self.gfServe_params['mask']:
                cmd += '-mask '
            if self.gfServe_params['translate']:
                cmd += '-trans '

            cmd += 'start {gfhost} {gfport} {genome} &> /dev/null &'.format(gfhost=self.host, gfport=self.port, genome=os.path.basename(self.genome))
             
            #print cmd
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Execution failed for gfServer:\n"+e)
            sys.exit(1)
        os.chdir(cwd)
        time.sleep(5) # sleep 5 sec
        gfready = 0
        while not gfready:
            gfready = self.__read_gfServer_log(log_file)
            if gfready == -1:
                sys.stderr.write('gfServer start error! Check '+log_file+' file, exit.')
                sys.exit(1)
            time.sleep(5)
        sys.stderr.write('gfServer ready\n\n\n')
    #end method start()   
        
    def stop(self):
        """Stops the gfServer
        
        """
        cmd = [
            self.exe,
            'stop', 
            self.host,
            str(self.port)
        ]
        p = subprocess.Popen(cmd)
        p.communicate()
    #end method stop()
        
    
    def __read_gfServer_log(self, infile):
        with open(infile, 'r') as f:
            for i in f:
                if 'Server ready for queries!' in i:
                    f.close()
                    return 1
                elif 'gfServer aborted' in i or ('error' in i and not 'getsockopt error' in i):
                    f.close()
                    return -1
        return 0
    #end method read_gfServer_log()
    
    def isPCR(self, forward, reverse, out='bed', maxsize=4000, minPerfect=15, minGood=15):
        """Performs a in-silico PCR search for amplicons resulting from primers `forward` and `reverse`
        
        Parameters:
            forward: (string) Forward primer to use
            reverse: (string) Reverse primer to use
            out: (string) Output format for results
            maxsize: (int) Maximum amplicon size allowed
            minPerfect: (int) see isPCR documentation
            minGood: (int) see isPCR documentation
            
        Returns:
            if `out` == 'bed' then data is parsed into a pybedtools.BedTool instance
            otherwise, the raw output data is returned
        """
        cmd = [
            os.path.join(os.path.dirname(self.exe), 'gfPcr'), 
            self.host,
            str(self.port),
            os.path.dirname(self.genome), 
            forward,
            reverse,
            'stdout',
            '-out=bed',
            '-maxSize={}'.format(maxsize),
            '-minPerfect={}'.format(minPerfect),
            '-minGood={}'.format(minGood),
            #'-name=%s' % (primerpair.name,) #using this option seems to make gfPcr throw a (harmless?) error freeing memory, we dont really need the name anyway so do not provide.
        ]
        
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdoutdata, stderrdata) = p.communicate()
        
        if out == 'bed':
            hits =  pybedtools.BedTool(stdoutdata, from_string=True)
            return hits
        else:
            return stdoutdata
    #end method isPCR()   
#end class gfServer   
    