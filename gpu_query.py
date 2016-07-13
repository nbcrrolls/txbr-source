import os
import sys
import subprocess
import getopt

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

os.system('make --directory=./txbr/txbr/bckprj/cuda')

hostQuery = "rocks list host"
deviceQuery = os.path.join(txbr_dir, "txbr", "txbr", "bckprj", "cuda", "deviceQuery")

# First, get the list of the different hosts

prc = subprocess.Popen(hostQuery, shell=True, stdout=subprocess.PIPE)
hosts_details = prc.communicate()[0].split('\n')

all_hosts = [  h.split(":")[0] for h in hosts_details if h.find(":")!=-1 ]
hosts = [ h for h in all_hosts if h.find(node)!=-1 ]

if len(hosts)==0:
    print "No node in %s" %(all_hosts)

# Second, get the number of GPUs on each host

gpu_dict = {}
for host in hosts:
   prc = subprocess.Popen("ssh -x %s %s -n" %(hosts[0],deviceQuery), shell=True, stdout=subprocess.PIPE)
   n = int(prc.communicate()[0])
   gpu_dict[host] = n
   print "host: %s    number of GPUs: %i" %(host,n)
   
# Third, create the configuration file

infofile = os.path.join(txbr_dir,'data','machinefile.gpu.txt')

try:
    f = open(infofile,'w')
    f.write(str(gpu_dict))
    f.close()
except:
    print "Error with file %s" %infofile


