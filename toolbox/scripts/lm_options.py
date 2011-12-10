#Candidate OPTIONS
#max width
#keep cwd

#!/usr/bin/env python2.5

"""Master script for linkage mapper.

Reguired Software:
Python 2.5
Numpy

"""

__filename__ = "lm_options.py"
__version__ = "0.6.5"


import arcgisscripting
import os
import sys
import traceback
import ConfigParser
import shutil

gp = arcgisscripting.create(9.3)
gprint = gp.addmessage

# def write_config_file(configFile, options):
    # """Write Linkage Mapper options to config file.

    # """
    # try:
        
    # except:
        # raise_python_error(__filename__)
        

def lm_options():
    """Set options for linkage mapper.

    """
    try:        
        options = {}
        PROJECTDIR = sys.argv[1]
        options['lcc_cutoff'] = sys.argv[2]
        configFile = os.path.join(PROJECTDIR, "lm_options.ini")
        config = ConfigParser.ConfigParser() 
        sections={}
        section='Linkage Mapper Options'
        sections['lcc_cutoff']=section
        
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
        gprint('Options written to: ' + configFile)
    except:
        gprint('****Failed to set Linkage Mapper options. Details follow.****')
        raise_python_error(__filename__)
         
    return


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

if __name__ == "__main__":
    lm_options()
    