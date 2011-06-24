##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_util.py 726 2010-12-02 15:33:49Z mcrae $
##

import ConfigParser, os, string, gzip
import numpy
import time

#NEEDED FOR GDAL ON BRAD'S LAPTOP
# import sys
# sys.path.append('/Library/Python/2.5/site-packages/GDAL-1.5.2-py2.5-macosx-10.5-i386.egg') 

##gdal_available = False #GDAL DISABLED UNTIL READ BUG RESOLVED
##try:
##    from osgeo import gdal_array, gdal
##    from osgeo.gdalconst import *
##except ImportError:
##    gdal_available = False

# Disable GDAL as it is error-prone for some cases for now. VS - 4/5/09
gdal_available = False
    
from string import split
from numpy import loadtxt, where
    
def readConfigFile(configFile):

    if os.path.isfile(configFile)==False:
        raise RuntimeError('File "'  + configFile + '" does not exist')

    config = ConfigParser.ConfigParser()
    config.read(configFile)
    options={}
    
    options['low_memory_mode']=False        
    options['use_mask']=False
    options['mask_file']='None' 
    options['use_included_pairs']=False
    options['included_pairs_file']='None' 
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None' 
    options['data_type']='raster' 
    options['version']='unknown'
    for section in config.sections():
        for option in config.options(section):
            try:
                options[option]=config.getboolean(section, option)
            except:
                options[option]=config.get(section, option)
    return options

def writeConfigFile(configFile, options):
   
    config = ConfigParser.ConfigParser()
 
    sections={}
    section='Version'
    sections['version']=section
    
    section='Connection scheme for raster habitat data'
    sections['connect_four_neighbors_only']=section
    sections['connect_using_avg_resistances']=section
    
    section='Short circuit regions (aka polygons)'
    sections['use_polygons']=section
    sections['polygon_file']=section
    
    section='Options for advanced mode'
    sections['source_file']=section
    sections['ground_file']=section
    sections['ground_file_is_resistances']=section
    sections['use_unit_currents']=section
    sections['use_direct_grounds']=section
    sections['remove_src_or_gnd']=section
    
    section='Calculation options'
    sections['solver']=section
    sections['print_timings']=section
    sections['low_memory_mode']=section    
    
    section='Output options'
    sections['output_file']=section
    sections['write_cur_maps']=section
    sections['write_cum_cur_map_only']=section
    sections['log_transform_maps']=section
    sections['write_volt_maps']=section
    sections['compress_grids']=section
    
    section='Mask file'
    sections['use_mask']=section
    sections['mask_file']=section  
    
    section='Options for pairwise and one-to-all and all-to-one modes'
    sections['use_included_pairs']=section
    sections['included_pairs_file']=section
    sections['point_file']=section
    sections['point_file_contains_polygons']=section
    
    section='Options for one-to-all and all-to-one modes'
    sections['use_variable_source_strengths']=section
    sections['variable_source_file']=section
   
    section='Habitat raster or graph'
    sections['habitat_file']=section
    sections['habitat_map_is_resistances']=section
    
    section="Circuitscape mode"
    sections['scenario']=section
    sections['data_type']=section

    if options['ground_file_is_resistances']=='not entered':
        options['ground_file_is_resistances'] = False
    if options['point_file_contains_polygons']=='not entered':
        options['point_file_contains_polygons'] = False
 
    for option in sections:
        try:
            config.add_section(sections[option])
        except:
            pass
    for option in sections:
        config.set(sections[option], option, options[option])

    f = open(configFile, 'w')
    config.write(f)
    f.close()
 


def setDefaultOptions():
    options = {}
    options['data_type']='raster' 
    options['version']='unknown'
    options['low_memory_mode']=False
    options['scenario']='not entered'
    options['habitat_file']='(Browse for a habitat map file)'
    options['habitat_map_is_resistances']=True
    options['point_file']='(Browse for file with locations of focal points or areas)'
    options['point_file_contains_polygons']=False
    options['connect_four_neighbors_only']=True
    options['connect_using_avg_resistances']=True
    options['use_polygons']=False
    options['polygon_file']='(Browse for a short-circuit region file)'
    options['source_file']='(Browse for a current source file)'
    options['ground_file']='(Browse for a ground point file)'
    options['ground_file_is_resistances']=True
    options['use_unit_currents']=False
    options['use_direct_grounds']=False
    options['remove_src_or_gnd']='not entered'
    options['output_file']='(Choose a base name for output files)'
    options['write_cur_maps']=False
    options['write_cum_cur_map_only']=False
    options['log_transform_maps']=False
    options['write_volt_maps']=False
    options['solver']='cg+amg'
    options['compress_grids']=False
    options['print_timings']=False
    options['use_mask']=False
    options['mask_file']='None' 
    options['use_included_pairs']=False
    options['included_pairs_file']='None' 
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None' 
    
    
    return options

def checkOptions(options):
    if options['scenario']=='not entered':
        all_options_entered=False
        message = 'Please choose a scenario' 
        
    elif options['habitat_file']=='(Browse for a habitat map file)':
        all_options_entered=False
        message = 'Please choose a raster habitat map file'
	    
    elif options['habitat_map_is_resistances']=='not entered':
        all_options_entered=False
        message = 'Please choose a habitat data type'
            
    elif options['scenario']=='pairwise'and options['point_file']=='(Browse for a file with focal points or regions)':
        all_options_entered=False
        message = 'Please choose a focal node file'

    elif options['scenario']=='one-to-all'and options['point_file']=='(Browse for a file with focal points or regions)':
        all_options_entered=False
        message = 'Please choose a focal node file'
        
    elif options['connect_four_neighbors_only']=='not entered':
        all_options_entered=False
        message = 'Please choose a cell connection scheme'
            
    elif options['connect_using_avg_resistances']=='not entered':
        all_options_entered=False 
        message = 'Please choose a cell connection calculation'
            
    elif options['scenario']=='advanced'and options['source_file']=='(Browse for a current source file)':
        all_options_entered=False 
        message = 'Please enter a current source file'
            
    elif options['scenario']=='advanced'and options['ground_file']=='(Browse for a ground point file)':
        all_options_entered=False
        message = 'Ground point file does not exist!'
            
    elif options['scenario']=='advanced'and options['ground_file_is_resistances']=='not entered':
        all_options_entered=False 
        message = 'Please choose a ground data type'
            
    elif options['use_polygons']==True and options['polygon_file']=='(Browse for a short-circuit region file)':
        all_options_entered=False
        message = 'Please enter a short-circuit region file or uncheck this option'
            
    elif options['output_file']=='(Choose an output file name)':
        all_options_entered=False
        message = 'Please choose an output file name'
            
    else:
      all_options_entered=True
      message='None'

    return all_options_entered, message      


def read_header(filename):
    if os.path.isfile(filename)==False:
        raise RuntimeError('File "'  + filename + '" does not exist')
    f = open(filename, 'r')
    try:
        [ign, ncols] = string.split(f.readline())
    except ValueError:
            raise  RuntimeError('Unable to read ASCII grid: "'  + filename + '". If file is a text list, please use .txt extension.')
    ncols = int(ncols)
    [ign, nrows] = string.split(f.readline())
    nrows = int(nrows)
    [ign, xllcorner] = string.split(f.readline())
    xllcorner = float(xllcorner)
    [ign, yllcorner] = string.split(f.readline())
    yllcorner = float(yllcorner)
    [ign, cellsize] = string.split(f.readline())
    cellsize = float(cellsize)
   
    try:
        [ign, nodata] = string.split(f.readline())
        try:
            nodata= int(nodata)
        except ValueError:
            nodata= float(nodata)
    except ValueError:
        nodata=False
  
    f.close()
 
    
    return ncols, nrows, xllcorner, yllcorner, cellsize, nodata 

def reader(filename, type):
    if os.path.isfile(filename)==False:      
        raise RuntimeError('File "'  + filename + '" does not exist')
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)

    if gdal_available == True:
        map = numpy.float64(gdal_array.LoadFile(filename))
        if nodata==True:    
            map = where(map==nodata, -9999, map)
    else:
        if nodata==False:
            map = loadtxt(filename, skiprows=5, dtype=type)
        else:
            map = loadtxt(filename, skiprows=6, dtype=type)
            map = where(map==nodata, -9999, map)
    if nrows==1:
        temp=numpy.zeros((1,map.size))
        temp[0,:]=map
        map=temp
    if ncols==1:
        temp=numpy.zeros((map.size,1))
        temp[:,0]=map
        map=temp       

    return map

def writer(file, data, state, compress):
#     if gdal_available == True:
#         driver = gdal.GetDriverByName ('MEM')
#         out = driver.Create (file, 
#                              xsize=data.shape[0],
#                              ysize=data.shape[1],
#                              bands=1, 
#                              eType=GDT_Float64)
#         out.SetGeoTransform([state['xllcorner'],  # top left x
#                              state['cellsize'],   # w-e pixel resolution
#                              0,                   # rotation
#                              state['yllcorner'],  # top left y
#                              0,                   # rotation
#                              state['cellsize'],   # n-s pixel resolution
#                              ])
#         band = out.GetRasterBand(1)
#         band.SetNoDataValue(state['nodata'])
#         gdal_array.BandWriteArray(band, data)

#         driver = gdal.GetDriverByName ('AAIGrid')
#         out_aai = driver.CreateCopy (file, out)
#     else:

    if gdal_available == True:
        format = "MEM"
        driver = gdal.GetDriverByName( format )
        dst_ds = driver.Create( file, len(data[0]), len(data),1,gdal.GDT_Float32)
        ull=state['yllcorner']+(state['cellsize'])*(len(data))
        dst_ds.SetGeoTransform([state['xllcorner'],  # left x
                             state['cellsize'],   # w-e pixel resolution
                             0,                   # rotation
                             ull,                 # upper left corner y
                             0,                   # rotation
                             state['cellsize'],   # n-s pixel resolution
                             ])
        dst_ds.GetRasterBand(1).WriteArray( data )
        format = 'AAIGrid'
        driver = gdal.GetDriverByName(format)
        dst_ds_new = driver.CreateCopy(file, dst_ds) #STILL GETTING LEADING SPACES.
        dst_ds = None
        
        
    else:
        f = False
        if compress == True:
            file = file + '.gz'
            f = gzip.open(file, 'w')
        else:
            f = open(file, 'w')

        f.write('ncols         ' + str(state['ncols']) + '\n')
        f.write('nrows         ' + str(state['nrows']) + '\n')
        f.write('xllcorner     ' + str(state['xllcorner']) + '\n')
        f.write('yllcorner     ' + str(state['yllcorner']) + '\n')
        f.write('cellsize      ' + str(state['cellsize']) + '\n')
        f.write('NODATA_value  ' + str(state['nodata']) + '\n')
        
        delimiter = ''
        fmt = ['%.6f ']*state['ncols']
        format = delimiter.join(fmt)
        for row in data:
            f.write(format % tuple(row) + '\n')

        f.close()


def elapsed_time(startTime): 
    now=time.time()
    elapsed=now-startTime
    secs=int(elapsed)
    mins=int(elapsed/60)
    hours=int(mins/60)
    mins=mins-hours*60
    secs=secs-mins*60-hours*3600
    return hours,mins,secs






    # try:
#         options['low_memory_mode']=config.getboolean("circuitscape options", "low_memory_mode")
#     except:
#         options['low_memory_mode']=False        
#     options['habitat_file']=config.get("circuitscape options", "habitat_file")
#     options['scenario']=config.get("circuitscape options", "scenario")
#     options['habitat_file']=config.get("circuitscape options", "habitat_file")
#     options['habitat_map_is_resistances']=config.getboolean("circuitscape options", "habitat_map_is_resistances")
#     options['point_file']=config.get("circuitscape options", "point_file")
#     try:
#         options['point_file_contains_polygons']=config.getboolean("circuitscape options", "point_file_contains_polygons")
#     except:
#         options['point_file_contains_polygons']=config.get("circuitscape options", "point_file_contains_polygons")        
#     options['connect_four_neighbors_only']=config.getboolean("circuitscape options", "connect_four_neighbors_only")
#     options['connect_using_avg_resistances']=config.getboolean("circuitscape options", "connect_using_avg_resistances")
#     options['use_polygons']=config.getboolean("circuitscape options", "use_polygons")   
#     options['polygon_file']=config.get("circuitscape options", "polygon_file")
# 
#     try:
#         options['use_mask']=config.getboolean("BETA options", "use_mask")
#     except:
#         options['use_mask']=False
#     try:
#         options['mask_file']=config.get("BETA options", "mask_file")    
#     except:    
#         options['mask_file']='None' 
#     try:
#         options['use_included_pairs']=config.getboolean("BETA options", "use_included_pairs")
#     except:
#         options['use_included_pairs']=False
# 
#     try:
#         options['included_pairs_file']=config.get("BETA options", "included_pairs_file")    
#     except:
#         options['included_pairs_file']='None' 
# 
# 
#     try:
#         options['use_variable_source_strengths']=config.getboolean("BETA options", "use_variable_source_strengths")
#     except:
#         options['use_variable_source_strengths']=False
# 
#     try:
#         options['variable_source_file']=config.get("BETA options", "variable_source_file")        
#     except:
#         options['variable_source_file']='None' 
#     
#     options['source_file']=config.get("circuitscape options", "source_file")
#     options['ground_file']=config.get("circuitscape options", "ground_file")
#     try:
#         options['ground_file_is_resistances']=config.getboolean("circuitscape options", "ground_file_is_resistances")
#     except:
#         options['ground_file_is_resistances']=config.get("circuitscape options", "ground_file_is_resistances")        
#     options['use_unit_currents']=config.getboolean("circuitscape options", "use_unit_currents")
#     options['use_direct_grounds']=config.getboolean("circuitscape options", "use_direct_grounds")
#     options['remove_src_or_gnd']=config.get("circuitscape options", "remove_src_or_gnd")
#     options['output_file']=config.get("circuitscape options", "output_file")
#     options['write_cur_maps']=config.getboolean("circuitscape options", "write_cur_maps")
#     options['write_cum_cur_map_only']=config.getboolean("circuitscape options", "write_cum_cur_map_only")
#     options['log_transform_maps']=config.getboolean("circuitscape options", "log_transform_maps")
#     options['write_volt_maps']=config.getboolean("circuitscape options", "write_volt_maps")
#     options['solver']=config.get("circuitscape options", "solver")
#     options['compress_grids']=config.getboolean("circuitscape options", "compress_grids")
#     options['print_timings'] = config.getboolean("circuitscape options", "print_timings")
#   




# 
#     config.add_section("circuitscape options")
#     config.add_section("BETA options")
    
        #Need following options to be boolean 
#     config.set("circuitscape options", "iterate", options['iterate'])        
#     if options['ground_file_is_resistances']=='not entered':
#         options['ground_file_is_resistances'] = False
#     if options['point_file_contains_polygons']=='not entered':
#         options['point_file_contains_polygons'] = False
#     config.set("circuitscape options", "low_memory_mode", options['low_memory_mode'])       
#     config.set("circuitscape options", "scenario", options['scenario'])
#     config.set("circuitscape options", "habitat_file", options['habitat_file'])
#     config.set("circuitscape options", "habitat_map_is_resistances", options['habitat_map_is_resistances'])
#     config.set("circuitscape options", "point_file", options['point_file'])
#     config.set("circuitscape options", "point_file_contains_polygons", options['point_file_contains_polygons'])
#     config.set("circuitscape options", "connect_four_neighbors_only", options['connect_four_neighbors_only'])
#     config.set("circuitscape options", "connect_using_avg_resistances", options['connect_using_avg_resistances'])
#     config.set("circuitscape options", "use_polygons", options['use_polygons'])
#     config.set("circuitscape options", "polygon_file", options['polygon_file'])
# 
#     config.set("circuitscape options", "source_file", options['source_file'])
#     config.set("circuitscape options", "ground_file", options['ground_file'])
#     config.set("circuitscape options", "ground_file_is_resistances", options['ground_file_is_resistances'])
#     config.set("circuitscape options", "use_unit_currents", options['use_unit_currents'])
#     config.set("circuitscape options", "use_direct_grounds", options['use_direct_grounds'])
#     config.set("circuitscape options", "remove_src_or_gnd", options['remove_src_or_gnd'])
#     config.set("circuitscape options", "output_file", options['output_file'])
#     config.set("circuitscape options", "write_cur_maps", options['write_cur_maps'])
#     config.set("circuitscape options", "write_cum_cur_map_only", options['write_cum_cur_map_only'])
#     config.set("circuitscape options", "log_transform_maps", options['log_transform_maps'])#
#     config.set("circuitscape options", "write_volt_maps", options['write_volt_maps'])
#     config.set("circuitscape options", "solver", options['solver'])
#     config.set("circuitscape options", "compress_grids", options['compress_grids'])
#     config.set("circuitscape options", "print_timings", options['print_timings'])
# 
#     config.set("BETA options", "use_mask", options['use_mask'])
#     config.set("BETA options", "mask_file", options['mask_file']) 
#     config.set("BETA options", "use_included_pairs", options['use_included_pairs'])
#     config.set("BETA options", "included_pairs_file", options['included_pairs_file']) 
#     config.set("BETA options", "use_variable_source_strengths", options['use_variable_source_strengths'])
#     config.set("BETA options", "variable_source_file", options['variable_source_file']) 
    
