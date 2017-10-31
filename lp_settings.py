"""Linkage Priority user-configurable settings."""
# Authors: John Gallo and Randal Greene 2017


## Moved to GUI:
# CLOSEWEIGHT = 0.1 # relative closeness weight in CSP calculation (this plus following weights should sum to 1)
# PERMWEIGHT = 0.1 # relative permeability weight in CSP calculation
# CAVWEIGHT = 0.1 # relative core area value weight in CSP calculation
# ECIVWEIGHT = 0.0 # relative ECIV weight in CSP calculation
# CEDWEIGHT = 0.7 # relative climate envelope difference weight in CSP calculation (only used if climate envelope input provided)
# PROPCSPKEEP = 0.008 # proportion of top CSP values to keep (use 1 to keep all)
# TRUNCWEIGHT = 0.5  # truncated corridors weight in PERCI calculation (this plus following weight should sum to 1)
# ERCIWEIGHT = 0.5  # ERCI weight in PERCI calculation
# RESWEIGHT = 0.33  # resistance weight in CAV calculation (this plus following weights should sum to 1)
# SIZEWEIGHT = 0.33  # size weight in CAV calculation
# APWEIGHT = 0.34  # area/perimeter ratio weight in CAV calculation
# ECAVWEIGHT = 0.0  # ecav weight in CAV calculation (only used if ECAV input provided)
# CFCWEIGHT = 0.0  # relative Current Flow Centrality (CFC) weight in CAV calculation (only used if CF_Central present in Cores input)
# OCAVWEIGHT = 0.0  # ocav weight in CAV calculation (only used if OCAV input provided)


# Advanced Settings:
RELPERMNORMETH = 0  # relative permeability normalization method (use 0 for score range normalization; any other value for maximum value normalization)
RELCLOSENORMETH = 0  # relative closeness value normalization method (use 0 for score range normalization; any other value for maximum value normalization)
CALCLP = True  # calculate linkage priority
NORMCORRNORMETH = 0  # normalized corridor normalization method (use 0 for score range normalization; any other value for maximum value normalization)
RESNORMETH = 1  # resistance normalization method (use 0 for score range normalization; any other value for maximum value normalization)
SIZENORMETH = 1  # size normalization method (use 0 for score range normalization; any other value for maximum value normalization)
APNORMETH = 1  # area/perimeter ratio normalization method (use 0 for score range normalization; any other value for maximum value normalization)
ECAVNORMETH = 1  # ecav normalization method (use 0 for score range normalization; any other value for maximum value normalization)
CFCNORMETH = 1  # cfc normalization method (use 0 for score range normalization; any other value for maximum value normalization)
MINCPV = 0  # minimum corridor priority value (use 0 to keep all)
NORMALIZERCI = True  # normalize RCI
TRUNCNORMETH = 0  # truncated raster normalization method (use 0 for score range normalization; any other value for maximum value normalization)
CALCBP = True  # calculate blended priority (requires CALCLP above to also be True)
NORMALIZELP = True  # normalize Linkage Priority
NORMALIZEBP = True  # normalize Blended Priority
KEEPINTERMEDIATE = True # keep intermediate outputs for troubleshooting purposes
MAXCSPWEIGHT = 0.5 # relative max CSP value weight in CPV calculation
MEANCSPWEIGHT = 0.5 # relative mean CSP value weight in CPV calculation
