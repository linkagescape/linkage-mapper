#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor utility module.

"""

import os
import time

from cc_config import cc_env


def mk_proj_dir(in_dir):
    """Create directory off project folder if it does not already exist"""
    new_dir = os.path.join(cc_env.proj_dir, in_dir)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    return new_dir
    
class Ctimer():
    def __enter__(self):
        self.start = time.time()
    def __exit__(self, *args): 
        print time.time() - self.start    
