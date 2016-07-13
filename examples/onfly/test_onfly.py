import time
import os.path
import glob
import shutil
import re
import txbr.onthefly

BACKUP = 'backup'

def testotf( directory, basename, extension ):

    print "Test the on-the-fly-routine"

    print "directory: %s" %directory
    print "basename: %s" %basename
    print "extension: %s" %extension

    # Preparation work - Grab the initial micrographs

    src_dir = os.path.join( directory, "WIB" )

    file_pattern = '*%s*.%s' %( basename, extension )

    tilts = glob.glob( os.path.join( src_dir, BACKUP, file_pattern ) )

    angles = []
    anglefile = {}

    for tilt in tilts:
        angle_match = re.match('\S*_([+-.\d]*).%s' %extension, tilt)
        if angle_match:
            angle = float(angle_match.group(1))
            angles.append( angle )
            anglefile[angle] = tilt

    angles.sort()
    
    inputs = [ anglefile[angle] for angle in angles ]
    outputs = [ os.path.join( src_dir, os.path.split(anglefile[angle])[1]) for angle in angles ]

    for fout in outputs:
        if os.path.exists(fout): os.remove(fout)

    # Load the image stack

    wd = os.path.join( directory, "txbr-processing" )
    stack = txbr.utilities.loadStack( directory, basename, workDirectory=wd, scope="fei_titan", bin=4 )

    # Grab the initial micrographs

    delta_t = 20    # update every delta_t seconds

    rec = txbr.onthefly.QuickReconstruction( stack )
    rec.clear()

    for index,(input,output) in enumerate(zip(inputs,outputs)):
        shutil.copyfile(input,output)
        if index%5==4 or index==len(inputs):
            rec.update()
            time.sleep(delta_t)

 #   rec.finalizeBeadAlignement()

if __name__ == '__main__':

    options1 = {'directory':"/ncmirdata3/sphan/otf/ju", 'basename':"Ju_cont2.mrc", 'extension':txbr.utilities.TIFF}
    options2 = {'directory':"/ncmirdata3/sphan/otf/ju2", 'basename':"CCDBid_74225_Ju_cont2.mrc", 'extension':txbr.utilities.JPG}

    testotf(**options1)

    