import os
import sys
import subprocess
import getopt

LOAD = 0.8

flags1 = ""
flags2 = [ "node=" ]

try:
    opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
except getopt.GetoptError, err:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

node = "Tile"

for option,value in opts:
    if option in ('--node'):
        node = value

txbr_dir = os.path.abspath(".")

hostQuery = "rocks list host"

# First, get the list of the different hosts

prc = subprocess.Popen(hostQuery, shell=True, stdout=subprocess.PIPE)
hosts_details = prc.communicate()[0].split('\n')

all_hosts = [  h.split(":")[0] for h in hosts_details if h.find(":")!=-1 ]
hosts = [ h for h in all_hosts if h.find(node)!=-1 ]

# Second, get the number of CPUs on each host

cpu_dict = {}
for host in hosts:
   prc = subprocess.Popen("ssh -x %s \"%s\"" %(hosts[0], "python -c 'from multiprocessing import cpu_count; print cpu_count()'"), shell=True, stdout=subprocess.PIPE)
   n = int(prc.communicate()[0])
   print "host: %s    number of CPUs: %i" %(host,n)
   n = int(LOAD*n)
   cpu_dict[host] = n

if len(hosts)==0:
   prc = subprocess.Popen("python -c 'from multiprocessing import cpu_count; print cpu_count()'", shell=True, stdout=subprocess.PIPE)
   n = int(prc.communicate()[0])
   print "localhost:  number of CPUs: %i" %(n)
   n = int(LOAD*n)
   cpu_dict["localhost"] = n
   
# Third, create the configuration file

infofile = os.path.join(txbr_dir,'data','machinefile.cpu.txt')

try:
    f = open(infofile,'w')
    f.write(str(cpu_dict))
    f.close()
except:
    print "Error with file %s" %infofile


