#!/usr/bin/env python2.5
# Author: Brad McRae

import os
import sys
import shutil


try:
    import arcpy
    gp = arcpy.gp
    arcgisscripting = arcpy
    arc10 = True
except:
    arc10 = False
    import arcgisscripting
    gp = arcgisscripting.create()
    
gprint = gp.addmessage

_SCRIPT_NAME = "lm_delete_cwds"


def delete_cwd_dir():
    """Deletes cost-weighted distance directory and CWD rasters

    """               
    projectDir = sys.argv[1]
    cwdBaseDir = os.path.join(projectDir, "datapass\\cwd")
    try:
        if os.path.exists(cwdBaseDir):
            gp.delete_management(cwdBaseDir)
    except:
        try:
            if os.path.exists(cwdBaseDir):
                shutil.rmtree(cwdBaseDir)
        except:
            gp.AddError("Unable to delete cwd directory.  One of the rasters "
                        "might have been open in ArcMap.\n You may "
                        'need to re-start ArcMap to release the file lock.')
            for msg in range(0, gp.MessageCount):
                if gp.GetSeverity(msg) == 2:
                    gp.AddReturnMessage(msg)
                print gp.AddReturnMessage(msg)
            exit(0)        
        
    return

if __name__ == "__main__":
     delete_cwd_dir()