"""Linkage Priority user-configurable settings."""
# Authors: John Gallo and Randal Greene 2017


CALCLP = True  # calculate linkage priority
CALCBP = True  # calculate blended priority (requires CALCLP to also be True)

NORMALIZERCI = True  # normalize RCI
NORMALIZELP = True  # normalize Linkage Priority
NORMALIZEBP = True  # normalize Blended Priority

RELPERMNORMETH = "SCORE_RANGE"  # relative permeability normalization method
RELCLOSENORMETH = "SCORE_RANGE"  # relative closeness value normalization method
NORMCORRNORMETH = "SCORE_RANGE"  # normalized corridor normalization method
RESNORMETH = "MAX_VALUE"  # resistance normalization method
SIZENORMETH = "MAX_VALUE"  # size normalization method
APNORMETH = "MAX_VALUE"  # area/perimeter ratio normalization method
ECAVNORMETH = "MAX_VALUE"  # ecav normalization method
CFCNORMETH = "MAX_VALUE"  # cfc normalization method
TRUNCNORMETH = "SCORE_RANGE"  # truncated raster normalization method

MINCPV = 0  # minimum corridor priority value (use 0 to keep all)
MAXCSPWEIGHT = 0.5  # relative max CSP value weight in CPV calculation
MEANCSPWEIGHT = 0.5  # relative mean CSP value weight in CPV calculation

KEEPINTERMEDIATE = True  # keep intermediate outputs for troubleshooting purposes
