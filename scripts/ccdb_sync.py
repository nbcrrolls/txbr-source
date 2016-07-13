#!/usr/bin/python

import sys
import os
import os.path
import logging
import getopt
import ccdb

def usage():
    """Usage for the command *cdb_sync.py*."""

    print 'Usage: %s.py -d directory  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    --d directory"
    print "        The directory to copy"
    print "    --pid"
    print "        The project ID for the data"
    print "    --mpid"
    print "        The MPID of the data"
    print "    --check"
    print "        Check if it is possible to sync the data!"
    print "    -h or --help"
    print "        Help Information"

def main():
    '''Main routine for fixing names in the datajail'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hd:",["help","directory=","pid=","mpid=","check"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    keywords = {}

    keywords['directory'] = "."
    keywords['pid'] = None
    keywords['mpid'] = None
    check = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d","--directory"):
            keywords['directory'] = value_
        elif option_ in ("--pid"):
            try:
                keywords['pid'] = int(value_)
            except ValueError:
                print "Variable pid should be an integer!"
                sys.exit(0)
        elif option_ in ("--mpid"):
            try:
                keywords['mpid'] = int(value_)
            except ValueError:
                print "Variable mpid should be an integer!"
                sys.exit(0)
        elif option_ in ("--check"):
            check = True
        else:
            assert False, "unhandled option"

    if not os.path.exists(keywords['directory']):
        logging.warning("The source directory %s does not exists" %keywords['directory'])
        sys.exit(0)

    if keywords['mpid']==None:
        logging.warning("You must provide a MPID!")
        sys.exit(0)

    if keywords['pid']==None:
        logging.warning("You must provide a PID!")
        sys.exit(0)

    pathsrc = keywords['directory']
    pathdest = ccdb.mpidProcessDataDirectory(**keywords)
    
    if pathdest==None:
        logging.warning("The destination directory is not accessible")
        sys.exit(0)

    print "Sync data [%s] to [%s]" %(pathsrc,pathdest)

    sizeOfSource = ccdb.sizeOfDirectory(pathsrc,unit="Gb")
    diskUsage = ccdb.diskUsage(mountPoint=ccdb.getMountPoint(pathdest),unit="Gb")

    if sizeOfSource>2*diskUsage[1]:
        logging.warning("There is probably not enough space to sync the data")
        sys.exit(0)

    if check:
        rawdata = ccdb.mpidRawDataDirectory(**keywords)
        if rawdata:
            files = os.listdir(ccdb.mpidRawDataDirectory(**keywords))
            print files

    if not check:
        ccdb.sync( pathsrc, pathdest )


if __name__ == '__main__': main()
