import os.path

data_dir =  os.path.join(os.path.dirname(__file__),'./../../data')

numberOfGPUBoardsForNode = {}
infoGPUBoard = []

infofile = os.path.join(data_dir,'config/machinefile.gpu.txt')

def loadInfo(file):
    
    global numberOfGPUBoardsForNode
    
    if not os.path.exists(file): return
    
    try:
        f = open(file,'r')
        numberOfGPUBoardsForNode = eval(f.read().strip())
        f.close()
    except:
        return
    
    global infoGPUBoard
    
    del infoGPUBoard[:]
    
    nodes = sorted(numberOfGPUBoardsForNode.keys())
    
    for node in nodes:
        for deviceId in range(numberOfGPUBoardsForNode[node]):
            infoGPUBoard.append((node,deviceId))
    
    
def infoForGPUHosts():

    nodes = sorted(numberOfGPUBoardsForNode.keys())
        
    return ( nodes, numberOfGPUBoardsForNode, infoGPUBoard )


def numberOfGPUBoards():
    '''Returns the total number of available GPU boards.'''
        
    values = numberOfGPUBoardsForNode.values()
    
    n_gpu_boards = 0

    for value in values:
        n_gpu_boards += value
        
    return n_gpu_boards


def isGPUEnabled():
    
    return len(numberOfGPUBoardsForNode.keys())>0

    
loadInfo(infofile)
