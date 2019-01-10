#!/usr/bin/env python
"""
NAME: vcffilter.py
==================

DESCRIPTION
===========

Read vcf-file and filter variants.
Very crude method. Just for re-checking results.

INSTALLATION
============

Nothing special. Uses only standard libs.

USAGE
=====

python vcffilter.py *.vcf.gz


TODO
====

 - Make use of cyvcf for speed.


VERSION HISTORY
===============

0.0.1    20190110      Initial version.

LICENCE
=======
2018-2019, copyright Sebastian Schmeier
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
import logging


csv.field_size_limit(sys.maxsize)

__version__ = '0.0.1'
__date__ = '2019/01/10'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'

# For color handling on the shell
try:
    from colorama import init, Fore
    # INIT color
    # Initialise colours for multi-platform support.
    init()
    reset = Fore.RESET
    colors = {'success': Fore.GREEN,
              'error': Fore.RED,
              'warning': Fore.YELLOW,
              'info': ''}
except ImportError:
    sys.stderr.write('colorama lib desirable. ' +
                     'Install with "conda install colorama".\n\n')
    reset = ''
    colors = {'success': '', 'error': '', 'warning': '', 'info': ''}


def alert(atype, text, log, repeat=False):
    if repeat:
        textout = '{} [{}] {}\r'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)
    else:
        textout = '{} [{}] {}\n'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)

    log.write('{}{}{}'.format(colors[atype], textout, reset))
    if atype == 'error':
        sys.exit(1)


def success(text, log=sys.stderr):
    alert('success', text, log)


def error(text, log=sys.stderr):
    alert('error', text, log)


def warning(text, log=sys.stderr):
    alert('warning', text, log)


def info(text, log=sys.stderr, repeat=False):
    alert('info', text, log)


def parse_cmdline():
    """ Parse command-line args. """
    # parse cmd-line ----------------------------------------------------------
    description = 'Read vcf-file and filter based on annotation values. The defaults are from GATK.'
                  
    version = 'version {}, date {}'.format(__version__, __date__)
    epilog = 'Copyright {} ({})'.format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='{}'.format(version))
    parser.add_argument(
        'file',
        metavar='FILE',
        help='vcf-file.')
    parser.add_argument('--QD',
                        metavar='FLOAT',
                        type=float,
                        default=2.0,
                        help='Filter on QD > FLOAT. [default=2.].')
    parser.add_argument('--FS',
                        metavar='FLOAT',
                        type=float,
                        default=30.,
                        help='Filter on FS < FLOAT. [default=30.].')
    parser.add_argument('--DP',
                        metavar='FLOAT',
                        type=float,
                        default=10.,
                        help='Filter on DP > FLOAT. [default=10.].')
    parser.add_argument('--MQ',
                        metavar='FLOAT',
                        type=float,
                        default=40.,
                        help='Filter on MQ > FLOAT. [default=40.].')
    parser.add_argument('--MQRankSum',
                        metavar='FLOAT',
                        type=float,
                        default=-12.5,
                        help='Filter on MQRankSum > FLOAT. [default=-12.5].')
    parser.add_argument('--ReadPosRankSum',
                        metavar='FLOAT',
                        type=float,
                        default=-8.0,
                        help='Filter on ReadPosRankSum > FLOAT. [default=-8.0].')
    
    parser.add_argument('--warn',
        action="store_true",
        default=False,
        help='Do not throw an exception if the value could not be extracted '+ \
        ' from a vcf line. Instead only print warning to stderr.')
    parser.add_argument('--failed',
        metavar='FILE',
        type=str,
        default=None,
        help='vcf-File to store failed variants in. [default = None]')
    
    # if no arguments supplied print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename, 'rt')
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.open(filename, 'rt')
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.ZipFile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def main():
    """ The main funtion. """
    #logger = logging.getLogger(__name__)
    args, parser = parse_cmdline()

    outfileobj = sys.stdout

    if args.failed:
        if args.failed.split('.')[-1] == 'gz':
            outfileobj_failed = gzip.open(args.failed, 'wt')
        elif args.failed.split('.')[-1] == 'bz2':
            outfileobj_failed = bz2.open(args.failed, 'wt')
        else:
            outfileobj_failed = open(args.failed, "w") 
    

    reg_qd = re.compile(";QD=(.+?);")
    reg_dp = re.compile(";DP=(.+?);")
    reg_fs = re.compile(";FS=(.+?);")
    reg_mq = re.compile(";MQ=(.+?);")
    reg_rprs = re.compile(";ReadPosRankSum=(.+?);")
    reg_mqrs = re.compile(";MQRankSum=(.+?);")

    dict_regs = {"QD":reg_qd,
                 "DP":reg_dp,
                 "FS":reg_fs,
                 "MQ":reg_mq,
                 "ReadPosRankSum":reg_rprs,
                 "MQRankSum":reg_mqrs}

    dict_tests = {"QD":args.QD,
                  "DP": args.DP,
                  "FS": args.FS,
                  "MQ": args.MQ,
                  "ReadPosRankSum": args.ReadPosRankSum,
                  "MQRankSum": args.MQRankSum}
    
    try:
        fileobj = load_file(args.file)
    except IOError:
        error('Could not load file "{}". EXIT.'.format(args.file))        
       
    # delimited file handler
    csv_reader_obj = csv.reader(fileobj, delimiter="\t", quoting=csv.QUOTE_NONE)
    i = 0
    iYay = 0
    iNay = 0
    iV = 0
    iNotFound = 0
    for a in csv_reader_obj:
        i += 1

        if a[0][0] == "#":  # comment
            outfileobj.write("{}\n".format("\t".join(a)))
            if args.failed:
                outfileobj_failed.write("{}\n".format("\t".join(a)))
            continue

        fail = 0
        iV += 1
        for name,reg in dict_regs.items():
            res = reg.search(a[7])

            if not res:
                outstr = 'Could not find "{}" value. Removed variant.\n'.format(name) + \
                         'Line ({}): {}'.format(i,'\t'.join(a))
                if args.warn:
                    warning(outstr)
                    iNay += 1
                    iNotFound += 1
                    fail = 1
                else:
                    error(outstr)                       
            else:
                try:
                    value = float(res.group(1))
                except ValueError:
                    error("Could not convert {} to float.".format(res.group(1)))

                if name == "FS":
                    if value >= dict_tests[name]:
                        fail = 1
                        break
                else:
                    if value <= dict_tests[name]:
                        fail = 1
                        break
                        
        if not fail:
            iYay += 1
            # For printing to stdout
            # SIGPIPE is throwing exception when piping output to other tools
            # like head. => http://docs.python.org/library/signal.html
            # use a try - except clause to handle
            try:
                outfileobj.write("{}\n".format("\t".join(a)))
                # flush output here to force SIGPIPE to be triggered
                # while inside this try block.
                sys.stdout.flush()
            except BrokenPipeError:
                # Python flushes standard streams on exit; redirect remaining output
                # to devnull to avoid another BrokenPipeError at shut-down
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
                sys.exit(1)  # Python exits with error code 1 on EPIPE
        else:
            iNay += 1
            if args.failed:
                outfileobj_failed.write("{}\n".format("\t".join(a)))


    success("Variants in file: {}".format(iV))
    success("Variants passed all filters: {}".format(iYay))
    success("Variants failed at least one filter: {}".format(iNay))
    success("  Of those at least one filter could not been found for: {}".format(iNotFound))
    
    
    
    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())
