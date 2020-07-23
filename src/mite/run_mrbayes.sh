#!/bin/sh

##################################################################################
##                                                                              ##
##   This script executes the MrBayes.                                          ##
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

if [ -z "$1" ]
then
    echo "You have to specify the number of processors that will be used to run MrBayes!"
else
    for dir in ../../output/nexus/tres-especies; do
        printf '\n'
        echo $dir
        mpirun --use-hwthread-cpus -np $1 mb "$dir/mite.nex"
    done
fi
