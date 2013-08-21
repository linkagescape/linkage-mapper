#!/usr/bin/env python2.5

"""Aggregates resistance rasters to coarser cell size using average resistance.
    Written by Brad McRae
Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension (toolbox is Arc 10).
Python 2.5


"""

__filename__ = "raster_aggregator.py"
__version__ = "2013_0610"

import os.path as path
import arcgisscripting
import os.path as path
import sys
import traceback
gp = arcgisscripting.create(9.3)
gp.CheckOutExtension("Spatial")
gp.OverwriteOutput = True

gprint = gp.addmessage

def str2bool(pstr):
    """Convert ESRI boolean string to Python boolean type"""
    return pstr == 'true'
    
def raster_aggregator():
    """Main function 

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.  

    """
    try:
        
        OUTPUTDIR = sys.argv[1]  # Output directory  
        AG_FACTOR =  int(sys.argv[2])
        METHOD = sys.argv[3]
        SMOOTH = str2bool(sys.argv[4])
        RESRAS = {}#list of resistance rasters       
        RESRAS[1] = sys.argv[5]
        RESRAS[2] = sys.argv[6]
        RESRAS[3] = sys.argv[7]
        RESRAS[4] = sys.argv[8]
        RESRAS[5] = sys.argv[9]
        
        if AG_FACTOR < 2 or AG_FACTOR > 99:
            msg = ('ERROR: Cell factor must be between 2 and 99.')
            gp.AddError(msg)
            exit(1)        
        
        gp.RefreshCatalog(OUTPUTDIR)

        OUTPUTGDBNAME = 'agregated_rasters.gdb'
        OUTPUTGDB = path.join(OUTPUTDIR, path.basename(OUTPUTGDBNAME))
        if not gp.exists(OUTPUTGDB):
            gp.createfilegdb(OUTPUTDIR, path.basename(OUTPUTGDBNAME))
        
        numRasters = 1
        if RESRAS[2] != '#':
            numRasters = 2
        if RESRAS[3] != '#':
            numRasters = 3
        if RESRAS[4] != '#':
            numRasters = 4
        if RESRAS[5] != '#':
            numRasters = 5
            
        gprint('\nThere are ' + str(numRasters) + ' rasters to aggregate.')
        gprint('\nCell factor is ' + str(AG_FACTOR))
        gprint('\nCell sizes will be multiplied by this amount')
        for rasterNum in range(1,numRasters+1):
            inputRaster = RESRAS[rasterNum]
            gp.SnapRaster = inputRaster
            oldCellSize = gp.Describe(inputRaster).MeanCellHeight
            dir,fileName = path.split(inputRaster)  
        
            if SMOOTH == True and METHOD == "MEAN":
                gprint('\nSmoothing cell values by taking mean of ' + str(AG_FACTOR) +'x' + str(AG_FACTOR) + ' neighborhood')
                if len(fileName)>10:
                    smoothRasterFN = fileName[0:10] + '_sm'
                else:
                    smoothRasterFN = fileName + '_sm'
                smoothRasterFull = path.join(OUTPUTDIR,smoothRasterFN)
                
                InNeighborhood = ("Rectangle " + str(AG_FACTOR) + " " + 
                                     str(AG_FACTOR) + " CELL")                              

                InNoDataOption = "NODATA"

                # FocalStatistics
                gp.FocalStatistics_sa(inputRaster, smoothRasterFull, InNeighborhood, "", InNoDataOption)

                inputRaster = smoothRasterFull

            if len(fileName)>10:
                outRasterFN = fileName[0:10] + '_ag'
            else:
                outRasterFN = fileName + '_ag'
            agRasterFull = path.join(OUTPUTDIR,outRasterFN)
            if gp.exists(agRasterFull):
                gp.delete_management(agRasterFull)
            
            gp.workspace = OUTPUTDIR
            
            gprint('\nAggregating raster "' + str(inputRaster) + '"')
            gprint('Old size was ' + str(oldCellSize))
            gprint('New cell size will be ' + str(oldCellSize * AG_FACTOR))
            gprint('Aggregation method: ' + METHOD)
            gp.Aggregate_sa(inputRaster, outRasterFN, AG_FACTOR, METHOD, 
                            "TRUNCATE", "NODATA")        

            if SMOOTH == True and METHOD == 'MEAN':
                finalRasterFN = fileName + '_CellFactor' + str(AG_FACTOR) + '_' + METHOD + '_sm'
            else:
                finalRasterFN = fileName + '_CellFactor' + str(AG_FACTOR) + '_' + METHOD     
                    
            gp.SnapRaster = agRasterFull           

            gp.workspace = OUTPUTGDB
            gprint("\nAggregated raster will be saved to: " + finalRasterFN)
            gprint("in geodatabase " + OUTPUTGDB)
            gp.CopyRaster_management(agRasterFull, finalRasterFN)  

            #clean up
            try:
                gp.delete_management(agRasterFull)
            except:
                pass
            if SMOOTH == True:           
                try:
                    gp.delete_management(smoothRasterFull)
                except:
                    pass
                
        gprint('Done.')
        
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        raise_python_error(__filename__)

        
def raise_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    gp.AddError("Geoprocessing error on **" + line + "** of " + filename +
                " :")
    dashline(1)
    for msg in range(0, gp.MessageCount):
        if gp.GetSeverity(msg) == 2:
            gp.AddReturnMessage(msg)
        # dashline(2)
        print gp.AddReturnMessage(msg)
        # dashline(2)
    exit(0)


def raise_python_error(filename):

    """Handle python errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]

    gp.AddError("Python error on **" + line + "** of " + filename)
    gp.AddError(err)
    # dashline(2)
    exit(0)
        

def dashline(lspace=0):
    """Output dashed line in tool output dialog.

       0 = No empty line
       1 = Empty line before
       2 = Empty line after

    """
    if lspace == 1:
        gp.addmessage('\n')
    gp.addmessage('---------------------------------')
    if lspace == 2:
        gp.addmessage('\n')
        
        
if __name__ == "__main__":
    raster_aggregator()
