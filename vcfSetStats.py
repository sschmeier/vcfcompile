#!/usr/bin/env python
"""
NAME: vcfSetStats.py
====================

DESCRIPTION
===========

Read a vcf-file that is a result of "gatk3 CombineVariants" 
and compile statistics of the caller intersections.

INSTALLATION
============

Nothing special. Uses only standard libs.

USAGE
=====

python vcfSetStats.py test.vcf.gz


TODO
====

 - Make use of cyvcf for speed.


VERSION HISTORY
===============

0.0.1    2018911/06    Initial version.

LICENCE
=======
2019, copyright Sebastian Schmeier
s.schmeier@gmail.com // https://www.sschmeier.com

template version: 2.0 (2018/12/19)
"""
import sys
import os
import os.path
import argparse
import csv
import gzip
import bz2
import zipfile
import time
import re
import operator

csv.field_size_limit(sys.maxsize)

__version__ = "0.0.1"
__date__ = "2019/11/06"
__email__ = "s.schmeier@protonmail.com"
__author__ = "Sebastian Schmeier"

# For color handling on the shell
try:
    from colorama import init, Fore

    # INIT color
    # Initialise colours for multi-platform support.
    init()
    reset = Fore.RESET
    colors = {
        "success": Fore.GREEN,
        "error": Fore.RED,
        "warning": Fore.YELLOW,
        "info": "",
    }
except ImportError:
    sys.stderr.write(
        "colorama lib desirable. " + 'Install with "conda install colorama".\n\n'
    )
    reset = ""
    colors = {"success": "", "error": "", "warning": "", "info": ""}


def alert(atype, text, log, repeat=False):
    if repeat:
        textout = "{} [{}] {}\r".format(
            time.strftime("%Y%m%d-%H:%M:%S"), atype.rjust(7), text
        )
    else:
        textout = "{} [{}] {}\n".format(
            time.strftime("%Y%m%d-%H:%M:%S"), atype.rjust(7), text
        )

    log.write("{}{}{}".format(colors[atype], textout, reset))
    if atype == "error":
        sys.exit(1)


def success(text, log=sys.stderr):
    alert("success", text, log)


def error(text, log=sys.stderr):
    alert("error", text, log)


def warning(text, log=sys.stderr):
    alert("warning", text, log)


def info(text, log=sys.stderr, repeat=False):
    alert("info", text, log)


def parse_cmdline():
    """ Parse command-line args. """
    # parse cmd-line ----------------------------------------------------------
    description = 'Read a vcf-file that is a result of "gatk3 CombineVariants" and compile statistics of the caller intersections. Prints to standard out. Some stats go to standard error.'

    version = "version {}, date {}".format(__version__, __date__)
    epilog = "Copyright {} ({})".format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("--version", action="version", version="{}".format(version))
    parser.add_argument("file", metavar="FILE", help="vcf-file.")
    parser.add_argument(
        "--qual",
        dest="qual",
        metavar="NUMBER",
        type=float,
        default=0.0,
        help="Only consider variants with a QUAL value equal or greater than this value. [default = 0]",
    )
    parser.add_argument(
        "--snpeffType",
        metavar="TYPE",
        default=None,
        help='Only consider variants with this SnpEff effect annotation (HIGH, MODERATE, LOW, MODIFIER). [default: not considered"]',
    )

    # if no arguments supplied print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ["-", "stdin"]:
        filehandle = sys.stdin
    elif filename.split(".")[-1] == "gz":
        filehandle = gzip.open(filename, "rt")
    elif filename.split(".")[-1] == "bz2":
        filehandle = bz2.open(filename, "rt")
    elif filename.split(".")[-1] == "zip":
        filehandle = zipfile.ZipFile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    if args.snpeffType:
        reg_genes = re.compile("\|(HIGH|MODERATE|LOW|MODIFIER)\|(.+?)\|")

    try:
        fileobj = load_file(args.file)
    except IOError:
        error('Could not load file "{}". EXIT.'.format(f))

    reg_set = re.compile("set=(.+?)(?:;|$)")

    csv_reader_obj = csv.reader(fileobj, delimiter="\t", quoting=csv.QUOTE_NONE)
    i = 0  # number of variants in file
    iDroppedQual = 0
    iDroppedEff = 0
    iConsidered = 0
    callerSets = {}
    for a in csv_reader_obj:
        if a[0][0] == "#":  # comment
            continue

        i += 1

        if a[5] == ".":
            if args.qual > 0:
                iDroppedQual += 1
                continue
        else:
            qual = float(a[5])  # quality value
            if qual < args.qual:
                iDroppedQual += 1
                continue

        if args.snpeffType:
            res_genes = reg_genes.findall(a[7])
            if len(res_genes) == 0:
                iDroppedEff += 1
                continue

        iConsidered += 1
        res_set = reg_set.search(a[7])
        if not res_set:
            error("Could not extract set from line:\n{}\n".format("\t".join(a)))
        callers = res_set.groups()[0].split("-")
        callers.sort()
        callSet = "|".join(callers)
        callerSets[callSet] = callerSets.get(callSet, 0) + 1

    success("Variants in file: {}".format(i))
    success("Number of variants dropped due to QUAL: {}".format(iDroppedQual))
    success("Number of variants dropped due to EFF: {}".format(iDroppedEff))

    callerSets_sorted = sorted(callerSets.items(), key=operator.itemgetter(1))
    callerSets_sorted.reverse()

    outfileobj = sys.stdout
    # For printing to stdout
    # SIGPIPE is throwing exception when piping output to other tools
    # like head. => http://docs.python.org/library/signal.html
    # use a try - except clause to handle
    try:
        outfileobj.write("Set\tNumCallers\tNumVars\tPctVars\n")
        for t in callerSets_sorted:
            cset = t[0]
            numC = len(cset.split("|"))
            num = t[1]
            pct = num * 100.0 / iConsidered
            outfileobj.write("{}\t{}\t{}\t{}\n".format(cset, numC, num, pct))
        # flush output here to force SIGPIPE to be triggered
        # while inside this try block.
        sys.stdout.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shut-down
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE

    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == "__main__":
    sys.exit(main())
