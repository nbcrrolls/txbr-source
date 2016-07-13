import os.path

TXBR_EXTENSION = ".txbr"
FID_EXTENSION = ".fid"
MARKER_EXTENSION = ".mrk"
TRACK_EXTENSION = ".trk"
TRAJ_EXTENSION = ".trj"
ST_EXTENSION = ".st"
MRC_EXTENSION = ".mrc"
PREALI_EXTENSION = ".preali"
ALI_EXTENSION = ".ali"
XG_EXTENSION = ".xg"
PREXG_EXTENSION = ".prexg"
RAWTLT_EXTENSION = ".rawtlt"
MOSAIC_EXTENSION = ".mosaic"
MDOC_EXTENSION = ".mdoc"
MOD_EXTENSION = ".mod"


def extract_series_name_from_directory( directory ):
    '''Extract the common basename given a directory.'''

    if os.path.lexists(directory): names = os.listdir(directory)
    else: names = []

    names = [ os.path.splitext(name)[0] for name in names if name.endswith(ST_EXTENSION) ]

    return extract_series_name( names )


def extract_series_name( names ):
    '''Extract the common basename among a list of names.'''

    if names==None or len(names)==0: 
        return ( None, [] )

    lmax = min([len(name) for name in names])

    basename = ''
    for index in range(lmax):
        basename += names[0][index]
        for i in range(1,len(names)):
            if basename!=names[i][:index+1]:
                 basename = basename[:index]
                 extensions = [name[index:] for name in names]
                 return [basename,extensions]

    extensions = [name[lmax:] for name in names]
    extensions.sort()
    
    return ( basename, extensions )
        
        
def extract_series_name_from_dir( basename, dir, extension=TXBR_EXTENSION ):
    '''Extract the series names given a prefix <basename> and a directory <dir>.
    Files whose starts with <basename>, belongs to <dir> and ends up with
    <extension> are taking into account.
    Variable <basename> can be a comma separated string list
    '''
        
    if os.path.exists(dir):
        n = len(extension)
        basenames = basename.split(",")
        # First pass for a common names
        files = []
        for b in basenames:
            files.extend([ file[:-n] for file in os.listdir(dir) \
                     if file.startswith(b) and file.endswith(extension) ])
        basename, extensions = extract_series_name(files)
        if len(extensions)>1:
            try:
                extensions.remove('')
            except ValueError:
                pass
        elif len(extensions)==0:
            print "No file with extension %s in %s" %(extension,dir)
        return [ basename + ext for ext in extensions ]
    else:
        return [ basename ]
        

def extract_series_name_on_basename_extension( names, dir, extension=TXBR_EXTENSION ):
    '''Find a partition ( basename, extensions ) for the list <names> such that there is
    a '.fid' file named afetr the variable <basename> concatenated with variable <extension> 
    in the directory <dir>'''
    
    basenames = extract_series_name_from_dir( ",".join(names), dir, extension=extension )
   
    basename, extensions = extract_series_name( basenames )
    
    # Evntually check if there is a mosaic files...
    files = [ file[:-len(extension)] for file in os.listdir(dir) \
                 if file.endswith(extension) ]

    for mf in files:
        if basename.startswith(mf):
            prefix_size = len(mf)
            extensions = [ basename[prefix_size:] + ext for ext in extensions]
            basename = basename[:prefix_size]
            break
        
    return ( basename, extensions )


