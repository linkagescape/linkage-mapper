"""Linkage Priority user-configurable settings."""
# Authors: John Gallo and Randal Greene 2017

# Calculate Corridor Specific Value (CSP) or CSP & Blended Priority (BP)
CALCCSPBP = 2  # No_Calc=0, CSP=1, CSP_BP=2

RELPERMNORMETH = "SCORE_RANGE"  # relative permeability normalization method
RELCLOSENORMETH = "SCORE_RANGE"  # relative closeness value normalization method
RESNORMETH = "MAX_VALUE"  # resistance normalization method
SIZENORMETH = "MAX_VALUE"  # size normalization method
APNORMETH = "MAX_VALUE"  # area/perimeter ratio normalization method
ECAVNORMETH = "MAX_VALUE"  # ecav normalization method
CFCNORMETH = "MAX_VALUE"  # cfc normalization method
CANALOGNORMETH = "SCORE_RANGE"  # climate analog normalization method
CPREFERNORMETH = "SCORE_RANGE"  # climate preference normalization method
NORMCORRNORMETH = "SCORE_RANGE"  # normalized corridor normalization method

MINCPV = 0  # minimum corridor priority value (use 0 to keep all)
MAXCSPWEIGHT = 0.5  # relative max CSP value weight in CPV calculation
MEANCSPWEIGHT = 0.5  # relative mean CSP value weight in CPV calculation

HIGHERCE_COOLER = False  # higher climate envelop values are cooler

KEEPINTERMEDIATE = True  # keep intermediate outputs for troubleshooting purposes
