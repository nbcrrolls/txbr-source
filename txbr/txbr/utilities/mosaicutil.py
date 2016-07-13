import logging
import os.path
import re
import modl
import txbr

from txbr.utilities import FID_EXTENSION
from txbr.utilities import MARKER_EXTENSION
from txbr.utilities import TRACK_EXTENSION

MODEL_OFFSET = 1

log = logging.getLogger()

def mergeModelsFromConstraints( basename, extensions, constraints, keep_constraints_only=False ):
    '''Create marker files for series specified by its basename and extensions and a set of constraints'''
    
    displayConstraints(constraints)
         
    models = [ modl.Model("%s%s%s" %(basename, ext, FID_EXTENSION)) for ext in extensions ]
    
    for model in models: 
        model.loadFromFile()
        
    # Initialize the output models
    
    output_models = []
    
    for model,extension in zip(models,extensions):
        
        model_out_name = basename + extension + MARKER_EXTENSION
        model_out = modl.Model(model_out_name, template=model)
        
        output_models.append(model_out)
        
        trkPath = basename + extension + TRACK_EXTENSION
        if os.path.lexists(trkPath): os.unlink(trkPath)
        os.symlink(model_out_name,trkPath)
    
    # Populate the output models
    
    index_constraints_dict = {}
    markercount = 0

    for index,(model,extension) in enumerate(zip(models,extensions)):
        for object in model.objects:
            for indexOfContour,contour in enumerate(object.contours):
                key = "%s-%i" %(extension,indexOfContour+MODEL_OFFSET)
                if keep_constraints_only and not key in constraints:
                    continue
                if key in constraints and constraints[key] in index_constraints_dict:
                    indexOfConstraint = int(index_constraints_dict[constraints[key]])
                    object_out = output_models[index].objects[indexOfConstraint]
                else:
                    for index_out,model_out in enumerate(output_models):
                        if index==index_out:
                            object_out = model_out.addNewObject()
                            markercount += 1
                        else:
                            model_out.addNewObject()
                    if key in constraints:
                        index_constraints_dict[constraints[key]] = object_out.indexOfObject
                        object_out.name = constraints[key]
                for p in contour.points:
                    c = object_out.addNewContour()
                    c.addPoint(p[0],p[1],p[2])
    
    for model_out in output_models:
        model_out.save()
        
    if not os.path.lexists(txbr.ALIGN_DIR):
        os.mkdir(txbr.ALIGN_DIR)
        
    filename = os.path.join( txbr.ALIGN_DIR, '%s.mosaic.log' %(basename))
             
    f = open(filename,'w')
        
    log.info("Total of %i common markers!" %(len(index_constraints_dict)))
    for key,value in sorted(index_constraints_dict.iteritems(), key = lambda (k,v): (v,k)):
        c = "=".join([ pt for pt in constraints.keys() if constraints[pt]==key ])
        line = "%.4i -> %s (%-8s)" %(index_constraints_dict[key],c,key)
        log.info(line)
        f.write(line + '\n')
        
    f.close()
        
    log.info("Total number of markers: %i  (%i shared)" %(markercount, len(index_constraints_dict)))    
        
        
def displayConstraints( constraints ):
    '''Display the constraints'''
            
    log.info("Display constraints!")
    for key,value in sorted(constraints.iteritems(), key = lambda (k,v): (v,k)):
        log.info("%-8s: %s" %(key,value))

    
def extractMosaicConstraints( basename, extensions ):
    '''Extract constraints between fiducial files from a set of tilt series'''
    
    log.info("Extract the mosaic constraints from the fiducial files.")
        
    constraints = {}  # Just read the label of the contour...
        
    models = [ modl.Model("%s%s%s" %(basename, ext, FID_EXTENSION)) for ext in extensions ]
    
    for index,(model,extension) in enumerate(zip(models,extensions)):
        model.loadFromFile()
        for object in model.objects:
            for indexOfContour,contour in enumerate(object.contours):
                if contour.label!=None and len(contour.label)!=0:
                    key = "%s-%i" %(extension, indexOfContour+MODEL_OFFSET)
                    constraints[key]= contour.label
                    
    return constraints


def indicesIn( values ):
    
    output = []
    
    values = values.split(",")
    
    for value in values:
        item = value.split("-")
        if len(item)==1:
            output.append(int(item[0]))
        elif len(item)==2:
            obj_start = int(item[0])
            obj_end = int(item[1])
            n = obj_end - obj_start + 1
            output.extend([obj_start+i for i in range(n)])
            
    return output
            
            
def extractMosaicConstraintsFile( basename, extensions, mosaic_file ):
    '''Merge a list of models into one final model.'''
    
    log.info("Build a common fiducial file with a mosaic file constraint for:")
    for index,ext in enumerate(extensions):
        log.info("%i %s%s%s" %( index, basename, ext, FID_EXTENSION))
        
    constraints = {}

    p1 = re.compile("(\S+):(\S+)=(\S+):(\S+)")
    
    count = 0
        
    f = open(mosaic_file)
    
    for line in f.readlines():
        if line.startswith("#"): # Skip the constraint
            continue
        constraint1_match = p1.match(line)
        if constraint1_match:
            ext_a = constraint1_match.group(1)
            ext_b = constraint1_match.group(3)
            if not ext_a in extensions or not ext_b in extensions:
                print "NOT IN: %s %s" %(ext_a,ext_b)
                continue
            obj_a = indicesIn(constraint1_match.group(2))
            obj_b = indicesIn(constraint1_match.group(4))
            ext_a = [ ext_a for i in range(len(obj_a)) ]
            ext_b = [ ext_b for i in range(len(obj_b)) ]
            if len(obj_a)!=len(obj_b): 
                print "problem with %s" %line
                continue
        else:
            continue
        n = len(ext_a)
        for i in range(n):
            key_a = "%s-%i" %(ext_a[i],obj_a[i])
            key_b = "%s-%i" %(ext_b[i],obj_b[i])
            if not key_a in constraints and not key_b in constraints:
                constraints[key_a] = "constraint-%i" %count
                constraints[key_b] = "constraint-%i" %count
                count += 1
            elif key_a in constraints:
                constraints[key_b] = constraints[key_a]
            elif key_b in constraints:
                constraints[key_a] = constraints[key_b]
                 
    f.close()
    
    return constraints
        
    
def mergeMosaicModels( basename, extensions, mosaic_file=None, keep_constraints_only=False ):
    '''Merge a list of models into one final model. For the fiducial files'''
    
    if mosaic_file!=None and os.path.exists(mosaic_file):
        constraints = extractMosaicConstraintsFile( basename, extensions, mosaic_file )
    else:
        constraints = extractMosaicConstraints( basename, extensions )
        
    mergeModelsFromConstraints( basename, extensions, constraints, keep_constraints_only=keep_constraints_only )
    