import os.path

data_dir =  os.path.join(os.path.dirname(__file__),'./../../data')

numberOfCPUForNode = { 'localhost': 2 }
infoCPUBoard = []

infofile = os.path.join(data_dir,'config/machinefile.cpu.txt')

def loadInfo( file ):

    global numberOfCPUForNode

    if not os.path.exists(file): return

    try:
        f = open(file,'r')
        numberOfCPUForNode = eval(f.read().strip())
        f.close()
    except:
        return


def infoForCPUHosts():

    nodes = sorted(numberOfCPUForNode.keys())
        
    return ( nodes, numberOfCPUForNode )


def numberOfCPUs():
    '''Returns the total number of available CPU.'''
        
    values = numberOfCPUForNode.values()
    
    n = 0

    for value in values:
        n += value
        
    return n
    
    
loadInfo(infofile)
