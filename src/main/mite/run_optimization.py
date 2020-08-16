#!/usr/bin/env python3

##################################################################################
##                                                                              ##
##   This module implements optimization routines based on the                  ##
##   multiresolution method.                                                    ##
##                                                                              ##
##   This file is part of the featsel program                                   ##
##   Copyright (C) 2019 Gustavo Mendes Maciel                                   ##
##                                                                              ##
##   This program is free software: you can redistribute it and/or modify       ##
##   it under the terms of the GNU General Public License as published by       ##
##   the Free Software Foundation, either version 3 of the License, or          ##
##   (at your option) any later version.                                        ##
##                                                                              ##
##   This program is distributed in the hope that it will be useful,            ##
##   but WITHOUT ANY WARRANTY; without even the implied warranty of             ##
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              ##
##   GNU General Public License for more details.                               ##
##                                                                              ##
##   You should have received a copy of the GNU General Public License          ##
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.      ##
##                                                                              ##
##################################################################################

import argparse
import logging
import math
import os
import subprocess
from datetime import datetime

from phyloproteomicanalysis import MiteToNexusWriter
import create_nexus_files as nexus
import multiresolution as mr
import run_cadm as cadm

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
log_path = dir_path + '../../log/'
reports_path = dir_path + '../../reports/report_mite.txt'
nexus_filename = 'mite.nex'
tree_filename = 'mite.nex.con.tre'
nexus_best_filename = 'best.nex'
tree_best_filename = 'best.nex.con.tre'

global_best = ([0], [])
best_history = []

iterations_count = 0

# Returns the elapsed time
def elapsed(start, end):
    return str(end - start)

# Constructs the NEXUS file from the partitioned ion maps.
def construct_nexus(mtnw, args, partition):
    start = datetime.now()
    complete_path = nexus.run(
        mtnw, args.w, args.h, args.binary, args.f, partition,
        '/' + dt_string + '_' + args.ptype + '_' + str(args.pargs[0]) + 'x' +
        str(args.pargs[1]) + '-' + str(args.pargs[2]) + 'x' +
        str(args.pargs[3]) + '-' + str(args.pargs[4])
    )
    end = datetime.now()
    logging.info('NEXUS creation elapsed time: ' + elapsed(start, end))

    return complete_path

# Runs MrBayes with the NEXUS file as input.
def run_mrbayes(nexus_path, iter_num, args):
    mblog_path = log_path + 'mb/' + dt_string + '_' + str(args.pargs[0]) + 'x'
    mblog_path += str(args.pargs[1]) + '-' + str(args.pargs[2]) + 'x'
    mblog_path += str(args.pargs[3]) + '-' + str(args.pargs[4]) + '/'

    if not os.path.exists(os.path.dirname(mblog_path)):
        os.makedirs(os.path.dirname(mblog_path))

    outfile = open(mblog_path + 'mb_log' + str(iter_num) + '.txt', "w")
    logging.info(
        'Starting to run MrBayes with ' + str(args.nproc) + ' processors'
    )
    start = datetime.now()

    try:
        subprocess.check_call(
            ['mb', nexus_path + '/' + nexus_filename],
            stdout=outfile,
            stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError as e:
        logging.exception(e.returncode, e.output)

    end = datetime.now()
    logging.info('MrBayes elapsed time: ' + elapsed(start, end))

# Runs the CADM test, comparing the constructed tree with the mtDNA tree.
def run_cadm(nexus_path):
    logging.info(
        'Starting to run CADM test (NEXUS file: ' + nexus_path + '/' +
        nexus_filename + ')'
    )
    stats = cadm.run(nexus_path + '/' + tree_filename, reports_path)
    logging.info("CADM results:\n" + str(stats))

    return stats 

# Assesses the results
def assess_results(stats, partition, nexus_path):
    global global_best
    global best_history

    if stats[0] > global_best[0][0]:
        logging.info('W = ' + str(stats[0]) + ' is the new global best')
        global_best = (stats, partition)
        best_history.append((stats, partition))
        os.rename(nexus_path + '/' + tree_filename, nexus_path + '/' + tree_best_filename)
        os.rename(nexus_path + '/' + nexus_filename, nexus_path + '/' + nexus_best_filename)
        logging.info('Created/updated ' + nexus_best_filename + ' and ' + tree_best_filename + ' files')

# Run optimization procedure
def run_basic_optimization(args):
    global iterations_count
    mtnw = MiteToNexusWriter(args.xml_path, args.nexus_path, args.binary)
    shape = mtnw.mites[0].matrix.shape
    partition = []
    lstartr = args.pargs[0]
    lstartc = args.pargs[1]
    lendr = args.pargs[2]
    lendc = args.pargs[3]
    lstep = args.pargs[4]

    logging.info(
        '\nXML_PATH=' + args.xml_path +
        '\nNEXUS_PATH=' + args.nexus_path +
        '\nWINDOW_WIDTH=' + str(args.w) + ', WINDOW_HEIGHT=' + str(args.h) +
        '\nBINARY=' + str(args.binary) + ', f=' + str(args.f) +
        '\nNPROC=' + str(args.nproc) +
        '\nPARTITION TYPE=' + args.ptype + ' (' + str(args.pargs[0]) + 'x' +
        str(args.pargs[1]) + ' to ' + str(args.pargs[2]) + 'x' +
        str(args.pargs[3]) + ' with step ' + str(args.pargs[4]) + ')'
    )

    n = 1
    for i, j in zip(range(lstartr, lendr + 1, lstep), range(lstartc, lendc + 1, lstep)):
        logging.info('Starting iteration ' + str(n))
        start = datetime.now()

        if args.ptype == 'S':
            partition = mr.generate_uniform_partition(shape, i, j)
        elif args.ptype == 's':
            partition = mr.generate_square_partition(shape, i, j)

        nexus_complete_path = construct_nexus(mtnw, args, partition)
        run_mrbayes(nexus_complete_path, n, args)
        stats = run_cadm(nexus_complete_path)
        assess_results(stats, partition, nexus_complete_path)
        iterations_count += 1
        end = datetime.now()
        logging.info('Iteration ' + str(n) + ' elapsed time: ' + elapsed(start, end))
        n += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "nexus_path", type=str, metavar="NEXUS_PATH", help="the path to store "
        "generated nexus files"
    )
    parser.add_argument(
        "--xml_path", type=str, help="the path to the XML files that will "
        "be used to construct the ion maps", metavar="PATH",
        default=dir_path+"/../../input/ion_map/xml/"
    )
    parser.add_argument(
        "--window_width", type=int, help="the width of the window used to "
        "reduce the MITEs", required=True, dest="w"
    )
    parser.add_argument(
        "--window_height", type=int, help="the height of the window used to "
        "reduce the MITEs", required=True, dest="h"
    )
    parser.add_argument(
        "-b", "--binary", action="store_true", help="use binary MITEs to create the "
        "NEXUS files"
    )
    parser.add_argument(
        "--f", type=float, help="value for the frequency of 1's necessary to "
        "get 1 in the new window"
    )
    parser.add_argument(
        "-p", "--parallel", type=int, metavar="NPROC", help="run MrBayes in "
        "parallel with %(metavar)s number of processors", default=1,
        dest="nproc"
    )
    parser.add_argument(
        "--partition_type", type=str, choices=['S', 's'], required=True,
        help="defines how the ion maps will be partitioned (S: per matrix "
        "size, s: per part size)", dest="ptype"
    )
    parser.add_argument(
        "--partition_args", nargs=5, type=int,
        metavar=("min_row_arg", "min_col_arg", "max_row_arg", "max_col_arg", "step"),
        help="if partition type is S, defines the min and max number of times"
        "the matrix will be divided; if s, defines the parts min and max size",
        dest="pargs", required=True
    )

    args = parser.parse_args()

    if (args.binary and args.f is None):
        parser.error("The --binary argument requires --f to be set")

    if (not args.binary and args.f is not None):
        parser.error("The argument --f requires the --binary flag to be set")

    global dt_string
    now = datetime.now()
    dt_string = now.strftime("%d%m%Y_%H:%M:%S")

    if not os.path.exists(os.path.dirname(log_path)):
        os.makedirs(os.path.dirname(log_path))

    logging.basicConfig(
        filename=(
            log_path + dt_string + '_' + args.ptype + '_' + str(args.pargs[0])
            + 'x' + str(args.pargs[1]) + '-' + str(args.pargs[2]) + 'x' +
            str(args.pargs[3]) + '-' + str(args.pargs[4]) + '.txt'
        ),
        level=logging.INFO,
        filemode="w",
        format='%(asctime)s %(levelname)s:%(message)s',
        datefmt='%d/%m/%Y %H:%M:%S'
    )
    logging.info(
        'Starting optimization routine. Partitioning type: ' + args.ptype +
        ' (' + str(args.pargs[0]) + 'x' + str(args.pargs[1]) + ' to ' +
        str(args.pargs[2]) + 'x' + str(args.pargs[3]) + ' with step ' +
        str(args.pargs[4]) + ')'
    )
    start = datetime.now()
    run_basic_optimization(args)
    end = datetime.now()
    logging.info('Total elapsed time: ' + elapsed(start, end))
    logging.info('Number of iterations: ' + str(iterations_count))
    logging.info('Global best: ' + str(global_best))
    logging.info('Best history: ' + str(best_history))
