#!/usr/bin/python

import sys
import os
import os.path
import glob
import re
import getopt

TITAN_SCOPE = "titan"
SPIRIT_SCOPE = "spirit"

JAIL_DIR = "/ncmirdata5/rawdata/data/%s/%s/CCDBID_%s" # %(scope,user,mpid)

def usage():
    """Usage for the command *txbr_sirt.py*."""

    print 'Usage: %s.py -b basename[,...]  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    --scope value"
    print "        The scope with which the data has been aquired"
    print "    --user value"
    print "        The CCDB that aquired the data."
    print "    --mpid"
    print "        The MPID of the data to fix"
    print "    -h or --help"
    print "        Help Information"


def fix_mpid( mpid, scope="*", user="*", doRename=False, **keywords ):
    """Consolidate the file naming convention (relatively to the mpid value) in the data-jail"""

    path = JAIL_DIR %(scope,user,mpid)

    files = glob.glob(path)

    if len(files)==0:
        print "There is nothing to fix..."
        sys.exit(0)

    if len(files)>1:
        print "There are more than one instance"
        sys.exit(0)

    path= files[0]

    print "root-path: %s" %(path)

    __walk__( path, mpid, doRename=doRename )


def __walk__( path, mpid, doRename=False ):

    if not os.access( path, os.W_OK) and doRename:
        print "No permission to write on %s" %path
        sys.exit(0)

    head,trail = os.path.split(path)
    pattern = "CCDBid_(\d*)"
    mpid_match = re.match(pattern,trail)
    if mpid_match:
        newtrail = re.sub(pattern,"CCDBid_%s"%mpid,trail)
        newpath = os.path.join(head,newtrail)
        print "%s: %s ->%s" %(mpid,trail,newtrail)
        if os.path.islink(path):
            link = os.readlink(path)
            newlink = re.sub(pattern,"CCDBid_%s"%mpid,link)
            if doRename and not os.path.lexists(newpath):
                os.remove(path)
                os.symlink(newlink,newpath)
                print "Relinked %s to %s" %(path, newlink)
        elif doRename and not os.path.lexists(newpath):
            os.rename( path, newpath )
            print "Renamed %s to %s" %(path, newpath)
        if doRename: path = newpath

    if os.path.isdir(path):
        files = os.listdir(path)
        for file in files:
            __walk__( os.path.join(path,file), mpid, doRename=doRename)
    else:
        pass


def main():
    '''Main routine for fixing names in the datajail'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hb:u:",["help","basename=","scope=","user=","mpid="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    keywords = {}

    keywords['basename'] = None
    keywords['scope'] = "*"
    keywords['user'] = "*"
    keywords['mpid'] = None

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-b","--basename"):
            keywords['basename'] = value_
        elif option_ in ("--scope"):
            keywords['scope'] = value_
        elif option_ in ("--user"):
            keywords['user'] = value_
        elif option_ in ("--mpid"):
            try:
                keywords['mpid'] = int(value_)
            except ValueError:
                print "Variable mpid should be an integer!"
                sys.exit(0)
        else:
            assert False, "unhandled option"

    if keywords['mpid']==None:
        print "You must provide a MPID!"
        sys.exit(0)

    fix_mpid( **keywords)


if __name__ == '__main__': main()
