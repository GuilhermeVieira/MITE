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

from MiteToNexusWriter import MiteToNexusWriter
import create_nexus_files as nexus
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

# Returns the time elapsed
def elapsed(start, end):
    return str(end - start)

# Returns a uniform partition that will be used in all ion maps.
def generate_uniform_partition(shape, pcount_row, pcount_col):
    logging.info(
        'Generating an uniform partition (' + str(pcount_row) + 'x' +
        str(pcount_col) + ' parts)'
    )
    partition = []
    step_row = int(math.ceil(shape[0] / pcount_row))
    step_col = int(math.ceil(shape[1] / pcount_col))

    for i in range(0, shape[0], step_row):
        for j in range(0, shape[1], step_col):
            slicex = slice(i, i + step_row)
            slicey = slice(j, j + step_col)
            partition.append((slicex, slicey))

    return partition

# Constructs the NEXUS file from the partitioned ion maps.
def construct_nexus(mtnw, args, partition):
    start = datetime.now()
    complete_path = nexus.run(
        mtnw, args.window_width, args.window_height, args.binary,
        args.f, partition
    )
    end = datetime.now()
    logging.info('NEXUS creation elapsed time: ' + elapsed(start, end))

    return complete_path

# Runs MrBayes with the NEXUS file as input.
def run_mrbayes(nexus_path, run_num, nproc):
    mblog_path = log_path + 'mb/' + dt_string + '/'

    if not os.path.exists(os.path.dirname(mblog_path)):
        os.makedirs(os.path.dirname(mblog_path))

    outfile = open(mblog_path + 'mb_log' + str(run_num) + '.txt', "w")
    logging.info('Starting to run MrBayes with ' + str(nproc) + ' processors')
    start = datetime.now()

    try:
        subprocess.check_call(
            ['mpirun', '-np', str(nproc), 'mb', nexus_path + '/' + nexus_filename],
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
    stats = cadm.run(nexus_path, reports_path)
    logging.info("CADM results: " + str(stats))

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
    mtnw = MiteToNexusWriter(args.xml_path, args.nexus_path, args.binary)
    shape = mtnw.mites[0].matrix.shape
    global iterations_count

    logging.info(
        'MIN_PCOUNT=' + str(args.min_pcount) +
        '\nMAX_PCOUNT=' + str(args.max_pcount) +
        '\nXML_PATH=' + args.xml_path +
        '\nNEXUS_PATH=' + args.nexus_path +
        '\nWINDOW_WIDTH=' + str(args.window_width) + ', WINDOW_HEIGHT=' + str(args.window_height) +
        '\nBINARY=' + str(args.binary) + ', f=' + str(args.f) +
        '\nNPROC=' + str(args.parallel)
    )

    for i, j in zip(range(min_pcount_row, max_pcount_row + 1), range(min_pcount_col, max_pcount_col + 1)):
        logging.info('Starting iteration ' + str(i))
        start = datetime.now()
        partition = generate_uniform_partition(shape, i, j)
        nexus_complete_path = construct_nexus(mtnw, args, partition)
        run_mrbayes(nexus_complete_path, i, args.parallel)
        stats = run_cadm(nexus_complete_path)
        assess_results(stats, partition, nexus_complete_path)
        iterations_count += 1
        end = datetime.now()
        logging.info('Iteration ' + str(i) + ' elapsed time: ' + elapsed(start, end))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "min_pcount", type=int, metavar="MIN_PCOUNT", help="initial value for partition count"
    )
    parser.add_argument(
        "max_pcount", type=int, metavar="MAX_PCOUNT", help="last value for partition count"
    )
    parser.add_argument(
        "nexus_path", type=str, metavar="NEXUS_PATH", help="the path to store "
        "generated nexus files"
    )
    parser.add_argument(
        "--xml_path", type=str, help="the path to the XML files that will "
        "be used to construct the ion maps",
        default=dir_path+"/../../input/ion_map/xml/"
    )
    parser.add_argument(
        "--window_width", type=int, help="the width of the window used to "
        "reduce the MITEs", required=True
    )
    parser.add_argument(
        "--window_height", type=int, help="the height of the window used to "
        "reduce the MITEs", required=True
    )
    parser.add_argument(
        "-b", "--binary", action="store_true", help="use binary MITEs to create the "
        "NEXUS files"
    )
    parser.add_argument(
        "-f", "--f", type=float, help="value for the frequency of 1's necessary to "
        "get 1 in the new window"
    )
    parser.add_argument(
        "-p", "--parallel", type=int, metavar="NPROC", help="Run MrBayes in "
        "parallel with %(metavar)s number of processors", default=1
    )

    args = parser.parse_args()

    if (args.binary and args.f is None):
        parser.error("The --binary argument requires --f to be set")

    if (not args.binary and args.f is not None):
        parser.error("The argument --f requires the --binary flag to be set")

    global dt_string
    global min_pcount_row, min_pcount_column, max_pcount_row, max_pcount_column
    now = datetime.now()
    dt_string = now.strftime("%d%m%Y_%H:%M:%S")
    min_pcount_row = min_pcount_col = args.min_pcount
    max_pcount_row = max_pcount_col = args.max_pcount

    if not os.path.exists(os.path.dirname(log_path)):
        os.makedirs(os.path.dirname(log_path))

    logging.basicConfig(
        filename=(log_path + dt_string + '.txt'),
        level=logging.INFO,
        filemode="w",
        format='%(asctime)s %(levelname)s:%(message)s',
        datefmt='%d/%m/%Y %H:%M:%S'
    )
    logging.info(
        'Starting optimization routine. Partitioning ion maps up to '
        + str(max_pcount_row) + 'x' + str(max_pcount_col) + ' parts'
    )
    start = datetime.now()
    run_basic_optimization(args)
    end = datetime.now()
    logging.info('Total elapsed time: ' + elapsed(start, end))
    logging.info('Number of iterations: ' + str(iterations_count))
    logging.info('Global best: ' + str(global_best))
    logging.info('Best history: ' + str(best_history))
