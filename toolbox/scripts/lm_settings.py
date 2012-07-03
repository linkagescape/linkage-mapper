### USER SETTABLE VARIABLES
CALCNONNORMLCCS = False  # Mosiac non-normalized LCCs in step 5 (Boolean- set to True or False)
WRITETRUNCRASTER = True  # Truncate mosaic corridor raster at an upper limit, i.e., use width cutoff. (Boolean- set to True or False)
CWDTHRESH = 200000  # CWD corridor width cutoff to use in truncated raster (Integer)
MINCOSTDIST = None  # Minimum cost distance- any corridor shorter than this will not be mapped (Integer)
MINEUCDIST = None  # Minimum euclidean distance- any core areas closer than this will not be connected (Integer)
SAVENORMLCCS = False  # Save inidvidual normalized LCC grids, not just mosaic (Boolean- set to True or False)
SIMPLIFY_CORES = True  # Simplify cores before calculating distances (Boolean- set to True or False)
                       # This speeds up distance calculations in step 2,
                       # but Euclidean distances will be less precise.
