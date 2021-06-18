# Author: Brad McRae

import os
import sys
import shutil

import arcpy


_SCRIPT_NAME = "lm_delete_cwds"


def delete_cwd_dir(argv=None):
    """Deletes cost-weighted distance directory and CWD rasters

    """
    if argv is None:
        argv = sys.argv  # Get parameters from ArcGIS tool dialog

    projectDir = argv[1]
    cwdBaseDir = os.path.join(projectDir, "datapass\\cwd")
    try:
        if os.path.exists(cwdBaseDir):
            arcpy.Delete_management(cwdBaseDir)
    except Exception:
        try:
            if os.path.exists(cwdBaseDir):
                shutil.rmtree(cwdBaseDir)
        except Exception:
            arcpy.AddError("Unable to delete cwd directory.  One of the rasters "
                        "might have been open in ArcMap.\n You may "
                        'need to re-start ArcMap to release the file lock.')
            for msg in range(0, arcpy.GetMessageCount() - 1):
                if arcpy.GetSeverity(msg) == 2:
                    arcpy.AddReturnMessage(msg)
                print(arcpy.AddReturnMessage(msg))
            exit(0)

    return

if __name__ == "__main__":
     delete_cwd_dir()
