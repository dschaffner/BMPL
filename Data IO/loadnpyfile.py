import numpy as np
import os

def loadnpyfile(filename,supress=False):

    savefile = filename
    savefile = os.path.normpath(savefile)

    file = np.load(savefile)
    #if not supress:
    #    print 'Arrays loaded: '
    #    for arr in file.files:
    #        print arr    
    return file