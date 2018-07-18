#!/usr/bin/env python2.5
# Code written by Peter Hoffmann.  Modified by Brad McRae for Linkage Mapper.

import sys
import traceback
import time
import arcgisscripting
import lm_util as lu
gp = arcgisscripting.create(9.3)
# gwarn = gp.addwarning
gprint = gp.addmessage

class retry(object):
    default_exceptions = (Exception,)
    def __init__(self, tries, exceptions=None, delay=1):
        """
        Decorator for retrying a function if exception occurs
        
        tries -- num tries 
        exceptions -- exceptions to catch
        delay -- wait between retries
        """
        self.tries = tries
        if exceptions is None:
            exceptions = retry.default_exceptions
        self.exceptions =  exceptions
        self.delay = delay

    def __call__(self, f):
        def fn(*args, **kwargs):
            exception = None
            for _ in range(self.tries):
                try:
                    return f(*args, **kwargs)
                except self.exceptions, e:
                    import traceback
                    tb = sys.exc_info()[2]  # get the traceback object
                    # tbinfo contains the error's line number and the code
                    tbinfo = traceback.format_tb(tb)[1]
                    line = tbinfo.split(", ")[1]
                    filename = tbinfo.split(", ")[0]
                    filename = filename.rsplit("File ")[1]
                    lu.warn('--------------------------------------------------')
                    msg = ("The following error is being reported "
                           "on " + line + " of " + filename + ":")                        
                    lu.warn(msg)
                    # lu.write_log(msg)
                    lu.warn(str(e))
                    # lu.write_log(str(e))
                    lu.print_drive_warning()

                    delaytime = self.delay*10*(_ + 1)
                    lu.warn("Will try again. ")
                    lu.warn('---------RETRY #' + str(_+1) + ' OUT OF ' + 
                                  str(self.tries) + ' IN ' +
                                  str(delaytime) + ' SECONDS---------\n')
                    lu.snooze(delaytime)
                    exception = e

            #if no success after tries, raise last exception
            raise exception
        return fn
