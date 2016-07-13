import sys
import numpy
import struct
import util

from ctypes import *
from txbr import log

def load2DRemap( filename ):
    
    log.info("Load Remap configuration from file %s." %filename)

    file = open(filename,'rb')

    local_magnification = struct.unpack("d",file.read(sizeof(c_double)))[0]
    
    (order,nunber_of_terms,pw_n1,pw_n2) = struct.unpack("iiii",file.read(4*sizeof(c_int)))

    power_orders = struct.unpack(''.join(['d' for i in range(pw_n1*pw_n2)]),file.read(pw_n1*pw_n2*sizeof(c_double)))
    power_orders = numpy.resize(power_orders,(pw_n1,pw_n2))

    (nproj,n2,n3,n1) = struct.unpack("iiii",file.read(4*sizeof(c_int)))
    
    map2D = struct.unpack(''.join(['d' for i in range(n1*n2*n3)]),file.read(n1*n2*n3*sizeof(c_double)))
    map2D = numpy.resize(map2D,(n1,n2,n3))

    file.close()

    log.info("Local Magnification: %f" %local_magnification)
    log.info("Order: %i" %order)
    log.info("Number of Terms: %i" %nunber_of_terms)
    log.info("Shape of the Power Order Matrix: (%i,%i)" %(pw_n1,pw_n2))
    log.info("Power Order Matrix:")
    log.info("%s" %power_orders)
    log.info("Number of Projections: %i" %nproj)
    log.info("Shape of the Remap Array: (%i,%i,%i)" %(n1,n2,n3))
    log.info("map2D:")
    log.info("%s" %map2D)
    
    return (local_magnification, order, power_orders, map2D)
    
    
def store2DRemap( map2D, order, filename, localMagnification = 1.0):

    (ntilts,n,number_of_terms) = map2D.shape
    power_order = util.powerOrder(order,dim=2)
    ntilts = map2D.shape[0]

    file = open(filename,'wb')
    
    file.write(struct.pack("d", localMagnification))
    file.write(struct.pack("iiii", order, number_of_terms, 2, number_of_terms))
    file.write(struct.pack(''.join(['d' for i in range(number_of_terms)]),*power_order[0]))
    file.write(struct.pack(''.join(['d' for i in range(number_of_terms)]),*power_order[1]))
    file.write(struct.pack("iiii", ntilts, map2D.shape[1], map2D.shape[2], map2D.shape[0]))
    file.write(struct.pack(''.join(['d' for i in range(map2D.size)]),*map2D.ravel()))
    
    file.close()


if __name__ == '__main__':
    
    if len(sys.argv)>1: 
        
        filename = sys.argv[1]

        (order, map2D, power_orders) = load2DRemap(filename)
    
        if ('--store' in sys.argv): store2DRemap( map2D, order, '%s.store' %filename )
        
    else:
        
        print "Provide the name of the binary remap file"