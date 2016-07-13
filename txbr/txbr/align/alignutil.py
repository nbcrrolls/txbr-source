import contalign

def runOrthogonalAlignment( project, decimate=None, init=False ):
    """
    Run an orthogonal alignment for a given tomographic reconstruction project.
    """

    if project==None: return

    orthogonal = True
    multiple = True
    fix_angles = [ False, False, False ]
    fix_t = [ False, False, False ]

    keywords = {}

    keywords['model'] = contalign.__encode_model__( orthogonal, multiple, fix_angles, fix_t )
    keywords['decimate'] = decimate
    keywords['skip_ncg'] = False
    keywords['init'] = init

    alignment = contalign.TxBRContourAlign( project, **keywords )
    alignment.process()

    XYZ, mask = alignment.extract3DPointMarkers()
    projection = alignment.project.series[0].projection

    projmap = {}
    for index in range(projection.numberOfExposures()):
        projmap[index] = projection.getAffineTransformation(index)

    return projmap,XYZ

    

def runRegularAlignment( project, decimate=None, init=False ):
    """
    Run a regular alignment for a given tomographic reconstruction project.
    """

    orthogonal = False
    multiple = True
    fix_angles = [ False, False, False ]
    fix_t = [ False, False, False ]

    keywords = {}

    keywords['model'] = contalign.__encode_model__( orthogonal, multiple, fix_angles, fix_t )
    keywords['decimate'] = decimate
    keywords['skip_ncg'] = False
    keywords['init'] = init

    alignment = contalign.TxBRContourAlign( project, **keywords )
    alignment.process()

    XYZ, mask = alignment.extract3DPointMarkers()
    projection = alignment.project.series[0].projection

    projmap = {}
    for index in range(projection.numberOfExposures()):
        projmap[index] = projection.getAffineTransformation(index)

    return projmap,XYZ
