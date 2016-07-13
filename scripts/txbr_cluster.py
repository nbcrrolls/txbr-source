#!/usr/bin/python

import sys
import os
import getopt
import numpy
import modl
import stats

from enthought.mayavi import mlab

BCC = 'bcc'
FCC = 'fcc'
RHB = 'rhombohedral'

def usage():

    print
    print 'Usage: %s.py -i file' %os.path.basename(sys.argv[0])
    print
    
    
def plotCorrelation( file ):
    
    if file==None or not os.path.exists(file): 
        return
    
    XYZ = modl.loadAllPoints(file)
    
    stats.plotCorrelation(*stats.getPairCorrelation(XYZ))
    
    
def plotLattice( lattice='bcc', parameters=None ):
    
    if lattice==BCC: XYZ = stats.bcc()
    if lattice==FCC: XYZ = stats.fcc()
    if lattice==BCC: XYZ = stats.rhombohedral(parameters[0])
        
    plot3D(XYZ)
    mlab.show()
    
    
def main():
    
    flags1 = "hi:"
    flags2 = [ "help" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-i'):
            file = value
            
    plotCorrelation(file)
    #filename = "/ncmirdata3/sphan/virus_array/wk_txbr/FHV102709-26_z_-146.0.out_mode1fl2_trim.3dpkmod"
    
    
main()