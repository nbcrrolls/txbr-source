#!/home/sphan/usr/local/bin/python

import sys
import os.path
import getopt
import ccdb

def usage():

    print
    print 'Usage: %s.py [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -h or --help"
    print "        Help Information"

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hvi",['help','verbose','info'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    info = False

    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('--info'):
            info = True

    if info: ccdb.checkNodes()
    else:
        for host,info in ccdb.surveyNodes().iteritems():
            print "%s\n" %info


if __name__ == '__main__':
    main()
