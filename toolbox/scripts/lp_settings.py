"""Linkage Priority user-configurable settings."""
# Authors: John Gallo and Randal Greene 2017


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
