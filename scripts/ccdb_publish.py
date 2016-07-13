#!/usr/bin/python

import sys
import os.path
import getopt
import txbr.utilities
import mrc
import ccdb

from txbr.txbrconf import PUB_DIR
from gluon.sql import DAL

db = DAL('sqlite://storage.sqlite',folder='/home/sphan/Desktop/www/web2py/applications/txbr/databases',auto_import=True)

def usage():

    print "Usage: %s -b basename -d directory" % sys.argv[0]
    print
    print "    --pid value"
    print "        The PID for which to create thumbnails"
    print "    -h or --help"
    print "        Help Information"


def __roamProject__( projectid ):
    """Roam a project for mpids"""

    projects = ccdb.loadProjects()

    trials = projects[projectid]['trials']
    mpids = {}

    for mpid,description in trials.items():
        microscopy = ccdb.loadProjectMicroscopy( projectid, mpid, link2View=False )
        dir = microscopy[0]
        if dir==None:
            print "There is no rawdata for %s in %s" %(mpid,dir)
            continue
        workdirectory = os.path.realpath(os.path.join(dir,"..","..","processed_data"))
        files = os.listdir(dir)
        files = [ file for file in files if not os.path.isdir(os.path.join(dir,file)) ]
        files = [ file for file in files if file.endswith(".st") or file.endswith(".mrc") ]
        files.sort()
        if len(files)==0:
            print "There is no rawdata for %s" %mpid
            continue
        basename,extension = os.path.splitext(files[0])
        mpids[mpid] = { 'basename':basename, 'd':dir, 'wd':workdirectory }

    return mpids


def createThumbnailForProject( projectid, overWrite=False, test=False ):
    """Create thumbnails (if possible) for all the datasets (MPID) within a given
    project specified by its PID

    :param: The project id.
    :overWrite: if True, regenerate a thumbnail even if it already exists.
    """

    mpids = __roamProject__( projectid ) # All the mpids for a project

    print

    for mpid,description in mpids.items():

        record = db.microscopy_thumbnail(microscopy_nb=mpid)

        print "Need tumbnail for mpid: %i   %s" %(mpid,record==None)

        dir = description['d']
        wd = description['wd']
        basename = description['basename']

        if not os.access(wd,os.W_OK):
            print "%s is not writable!" %wd
            continue
            
        if record!=None and not overWrite: continue

        src_file = os.path.join( dir,"%s.st" %basename )

        if not os.path.exists(src_file) or (os.system("header %s > /dev/null" %src_file)==0):
            files = os.listdir(dir)
            for file in files:
                path = os.path.join(dir,file)
                if os.path.isdir(path): continue
                ismrc = (os.system("header %s > /dev/null" %path)==0)
                if ismrc:
                    src_file = path
                    break

        if not os.path.lexists(src_file): continue

        print "Will generate thumbnail from %s" %src_file

        mount = ccdb.getMountPoint(wd)
        s = os.statvfs(mount)
        space = (s.f_bavail*s.f_frsize)
        if space==0:
            print "There is no free space on %s!" %mount

        if test: continue   # Nothing to do (thumbnails are not generated)

        try: # Generate the thumbnail

            pub_dir = os.path.join( wd, PUB_DIR )
            if not os.path.exists(pub_dir): os.makedirs(pub_dir,mode=0777)

            ccdb.resetStat( pub_dir )

            gallery_dir = os.path.join( pub_dir, "gallery" )
            if not os.path.exists(gallery_dir): os.makedirs(gallery_dir,mode=0777)

            mrcVolume = mrc.MRCFile(src_file)
            filename = os.path.join( pub_dir, "%s.thumbnail" %basename )
            txbr.utilities.saveAsPNG( mrcVolume.getMiddleSliceInZ(), filename , size=(150,150))

        except:
            
            pass


        try: # Upload the thumbnail into the database

            filename = os.path.join( pub_dir, "%s.thumbnail" %basename )
            
            filename = filename + ".gif"

            stream = open(filename,'rb')

            if record!=None:
                db(db.microscopy_thumbnail.microscopy_nb==mpid).update(data=db.microscopy_thumbnail.data.store(stream,filename))
            else:
                db.microscopy_thumbnail.insert(microscopy_nb=mpid,data=db.microscopy_thumbnail.data.store(stream,filename))

            stream.close()

            db.commit()

        except:

            pass



def main():
    '''Main routine for generating thumbnails in ccdb'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hd:b:",["help","wd=","pid="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("--wd"):
            work_directory = value_
        elif option_ in ("--pid"):
            pid = int(value_)
        else:
            assert False, "unhandled option"

    createThumbnailForProject( pid )


if __name__ == '__main__': main()

