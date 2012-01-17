#!/usr/bin/env python2.5

"""Script to run code checkers and profilers"""
import os
import sys
import glob
from datetime import datetime
from subprocess import call


def main():
    """Script to run code checkers and profilers"""
    # Set variables
    os.chdir(os.path.dirname(__file__))
    srcdir = os.path.join("..", "toolbox", "scripts")
    tday = datetime.today()
    outdir = os.path.join("..", "..", "review", tday.strftime("%Y%m%d"))
    fdatetime = tday.strftime(".%Y%m%d%H%M")
    mainm = "lm_run.py"
    pfile = os.path.join(outdir, os.path.splitext(mainm)[0] +
                         tday.strftime(".%Y%m%d%H%M") + ".pstats")
    filelist = glob.glob(os.path.join(srcdir, "*.py"))

    # Get selection option
    print
    if len(sys.argv) > 1:
        choice = abs(int(sys.argv[1]))
    else:
        print "Code checking and profiling tools"
        print "================================="
        print "1. pep8"
        print "2. pyflakes"
        print "3. pylint"
        print "4. pychecker"
        print "5. gprof2dot"
        print "6. pycallgraph"
        print "7. runsnake"
        print "8. epydoc"
        print "9. run all steps"
        print
        choice = raw_input("Please select one of the above: ")
        print "\n\n"
        choice = abs(int(choice))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if choice < 10:
        print "Files to be reviewed:"
        for i in range(0, len(filelist)):
            print filelist[i]
        print
        print "Output directory:"
        print os.path.abspath(outdir)
        print

        if choice == 1:
            rpep8(outdir, srcdir, fdatetime)
        elif choice == 2:
            rpyflakes(outdir, srcdir, fdatetime)
        elif choice == 3:
            rpylint(outdir, fdatetime, filelist)
        elif choice == 4:
            rpychecker(outdir, fdatetime, mainm, filelist)
        elif choice == 5:
            rcprofile(pfile, mainm)
            rgprof2dot(outdir, fdatetime, pfile)
        elif choice == 6:
            rpycallgraph(outdir, fdatetime, mainm)
        elif choice == 7:
            rcprofile(pfile, mainm)
            rrunsnake(pfile)
        elif choice == 8:
            repydoc(filelist)
        elif choice == 9:
            rpep8(outdir, srcdir, fdatetime)
            rpyflakes(outdir, srcdir, fdatetime)
            rpylint(outdir, fdatetime, filelist)
            rpychecker(outdir, fdatetime, mainm, filelist)
            rcprofile(pfile, mainm)
            rgprof2dot(outdir, fdatetime, pfile)
            rpycallgraph(outdir, fdatetime, mainm)
            rrunsnake(pfile)
            repydoc(filelist)

    print
    print "Exiting"


# Functions to run tools
def rpep8(p8odir, p8sdir, p8datetime):
    """"Run pep8"""
    print "Reviewing code with pep8"
    fpep8 = open(os.path.join(p8odir, "pep8" + p8datetime + ".txt"), "w")
    call(["pep8", "-r", "--statistics", os.path.join(p8sdir, "*.py")],
         stdout=fpep8)
    fpep8.close()


def rpyflakes(flkodir, flksrcdir, flkdatetime):
    """Run pyflakes"""
    print "Reviewing code with pyflakes"
    fpyflakes = open(os.path.join(flkodir, "pyflk" + flkdatetime + ".txt"),
                     "w")
    call(["pyflakes", flksrcdir], stdout=fpyflakes, shell=True)
    fpyflakes.close()


def rpylint(lntodir, lntdatetime, lntflist):
    """Run pylint"""
    print "Reviewing code with pylint"
    for i in range(0, len(lntflist)):
        fpylint = open(os.path.join(lntodir, "pylint." +
            os.path.basename(lntflist[i]) + lntdatetime + '.txt'), "w")
        call(["pylint", lntflist[i]], stdout=fpylint, shell=True)
        fpylint.close()


def rpychecker(chkodir, chkdatetime, chkmainm, chkflist):
    """Run pychecker"""
    print "Reviweing code with pychecker"
    print os.path.splitext(chkmainm)[0]
    fpychecker = open(os.path.join(chkodir, "pychk" +
                      os.path.splitext(chkmainm)[0] + chkdatetime + ".txt"),
                      "w")
    call("pychecker --only -# 1000 -q -t -v -g -n -a -I -8 -1 -A -G -m -f"
         " lm_argv.py " + ' '.join(chkflist), stdout=fpychecker,
         shell=True)
    fpychecker.close()


def rcprofile(cppfile, cpmainm):
    """Run cProfiler"""
    print "Profiling code"
    call("Python -m cProfile -o " + cppfile + " " + cpmainm, shell=True)


def rgprof2dot(g2dodir, g2ddatetime, g2dpfile):
    """"Run gprof2dot"""
    print "Creating call graph using gprof2dot"
    ogpro = "gprof" + g2ddatetime + ".txt"
    ogprof = open(os.path.join(g2dodir, ogpro), "w")
    call('gprof2dot.py -f pstats ' + g2dpfile, shell=True, stdout=ogprof,
         cwd=g2dodir)
    ogprof.close()
    odot = "gprof2dot" + g2ddatetime + ".png"
    call('dot -Tpng -o ' + odot + " " + ogpro, shell=True, cwd=g2dodir)


def rpycallgraph(cgodir, cgdatetime, cgmainm):
    """"Run pycallgraph"""
    print "Creating call graph using pycallgraph"
    opycall = os.path.join(cgodir, "pycallgraph" + cgdatetime + ".png")
    call("pycallgraph -o " + opycall + " " + cgmainm, shell=True)


def rrunsnake(rspfile):
    """"Run runsnake"""
    print "Startup Run Snake"
    call(["runsnake", rspfile])


def repydoc(epyflist):
    """Run epydoc"""
    print "epydoc API documentation"
    epydir = os.path.join("..", "..", "review", "epydoc")
    call('epydoc.py --html -o ' + epydir + ' --name'
         ' "Linkage Mapper" --graph  all ' + ' '.join(epyflist),
         shell=True)


if __name__ == "__main__":
    main()
