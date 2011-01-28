import arcgisscripting

gp = arcgisscripting.create(9.3)

try:
    import numpy
except:    
    gp.AddError("Numpy is not installed.  Please get the Python 2.5-compatible version from http://sourceforge.net/projects/numpy/files/NumPy/1.4.1/numpy-1.4.1-win32-superpack-python2.5.exe/download ") 
    sys.exit(1) # end script process
gp.addmessage('\n------------------\nHOORAY, you already have Numpy installed!\n\n\n------------------')

    
