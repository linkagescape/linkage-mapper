"""Decorator for retrying a function if exception occurs.

Code written by Peter Hoffmann.  Modified by Brad McRae for Linkage Mapper.

http://peter-hoffmann.com/2010/retry-decorator-python.html

"""

import sys
import traceback

import lm_util as lu


class Retry(object):
    """Decorator for retrying a function if exception occurs."""

    default_exceptions = (Exception,)

    def __init__(self, tries, exceptions=None, delay=1):
        """Init retry decorator.

        tries -- num tries
        exceptions -- exceptions to catch
        delay -- wait between retries
        """
        self.tries = tries
        if exceptions is None:
            exceptions = Retry.default_exceptions
        self.exceptions = exceptions
        self.delay = delay

    def __call__(self, func_in):
        """Run decorator function."""
        def func_out(*args, **kwargs):
            """Run and if necessary retry running function."""
            exception = None
            for i in range(self.tries):
                try:
                    return func_in(*args, **kwargs)
                except self.exceptions as func_error:
                    tbobj = sys.exc_info()[2]  # Get the traceback object
                    # tbinfo contains the error's line number and the code
                    tbinfo = traceback.format_tb(tbobj)[1]
                    line = tbinfo.split(", ")[1]
                    filename = tbinfo.split(", ")[0]
                    filename = filename.rsplit("File ")[1]
                    lu.warn('------------------------------------------------'
                            '--')
                    msg = ("The following error is being reported "
                           "on " + line + " of " + filename + ":")
                    lu.warn(msg)
                    lu.warn(str(func_error))
                    lu.print_drive_warning()

                    delaytime = self.delay * 10 * (i + 1)
                    lu.warn("Will try again. ")
                    lu.warn('---------RETRY #' + str(i + 1) + ' OUT OF ' +
                            str(self.tries) + ' IN ' +
                            str(delaytime) + ' SECONDS---------\n')
                    lu.snooze(delaytime)
                    exception = func_error

            # If no success after tries, raise last exception
            raise exception

        return func_out
