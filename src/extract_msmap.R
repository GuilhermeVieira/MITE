##################################################################################
##                                                                              ##
##   R script that generates ion intensity maps (MITEs) in CSV files from mzXML ##
##   raw data files, given the mzXML data path, the output path, the MS level   ##
##   (1 or 2) and the m/z resolution.                                           ##
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

library("MSnbase")
library("optparse")

option_list = list(
    make_option("--dataDir",
                type="character",
                default=NULL,
                help="mzXML data path",
                metavar="path"),
    make_option("--outDir",
                type="character",
                default=NULL,
                help="output path",
                metavar="path"),
    make_option("--msLevel",
                type="double",
                default=NULL,
                help="ms level",
                metavar="double"),
    make_option("--mzRes",
                type="double",
                default=NULL,
                help="mz resolution",
                metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dataDir) || is.null(opt$outDir) || is.null(opt$msLevel)
    || is.null(opt$mzRes)) {
    print_help(opt_parser)
    stop("all arguments must be supplied", call.=FALSE)
}

files <- dir(path=opt$dataDir, pattern="\\.mzXML$")

out_dir <- paste(opt$outDir, opt$mzRes, sep="/")
dir.create(out_dir, showWarnings=FALSE)

labels <- c("DATA_DIR", "OUTPUT_DIR", "MS_LEVEL", "MZ_RESOLUTION")
values <- c(opt$dataDir, out_dir, opt$msLevel, opt$mzRes)
cat(sprintf("%-16s%s\n", labels, values), sep="")
cat("\n")
cat("FILES:", "\n")

for (f in files) {
    cat(sprintf("%-55s", paste(opt$dataDir, f, sep="/")))
    fname <- paste(tools::file_path_sans_ext(f), "csv", sep=".")
    rawdata <- readMSData(paste(opt$dataDir, f, sep="/"), mode="onDisk")
    c <- chromatogram(rawdata)
    lowMz <- mz(c)[1]
    highMz <- mz(c)[2]
    hd <- fData(rawdata)
    ms <- which(hd$msLevel == opt$msLevel)
    mite <- MSmap(rawdata, ms, lowMz, highMz, opt$mzRes, zeroIsNA=FALSE)
    write.csv(msMap(mite), file=paste(out_dir, fname, sep="/"))
    cat("OK", "\n")
}
