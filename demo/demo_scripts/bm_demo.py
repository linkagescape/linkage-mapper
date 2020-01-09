"""Script to run Barrier Mapper tool.

Assumes Linkage Mapper scripts and demo data are in their default folders.

"""

import sys
import os


def in_params(demo_dir):
    """Define input paramaters."""
    prj_dir = os.path.join(demo_dir,  # Folder containing demos
                           'demo_run',  # Container folder to hold model runs
                           'lm')  # Linkage Mapper model run folder

    return (
        [os.path.basename(__file__),  # Script Name ** Do not modify **
        prj_dir,  # Project Directory (1)
        os.path.join(demo_dir,
                     'demoData\\resistances'),  # Resistance Raster (2)
        400,  # Minimum Detection Radius (3)
        1200,  # Maximum Detection Radius (4)
        400,  # Radius Step Value (5)
        "Calculate both maximum and sum",  # Method for combining (6)
        "true",  # Save barrier rasters for each search radius (7)
        "false"]  # Calculate percent improvment scores (8)
        )


def main():
    """Set path and run model."""
    demo_path = (os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

    sys.path.append(os.path.join(demo_path, '..\\toolbox\\scripts'))
    import barrier_master
    barrier_master.bar_master(in_params(demo_path))


if __name__ == "__main__":
    main()
