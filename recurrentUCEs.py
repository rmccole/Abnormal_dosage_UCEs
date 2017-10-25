#!/usr/bin/env python

"""
This script allows the user to check for recurrence of UCEs among one or more variant files. The UCEs should be drawn
from a master file which contains a unique identifier for each UCE, its coordinates and its type
(exonic/intronic/intergenic). The variant files should be normal 3-column interval files.

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

from __future__ import print_function
import argparse
import os
import sys


def getArgs(strInput=None):
    parser = argparse.ArgumentParser("Check one or more interval files for UCE recurrence, using a UCE master file to"
                                     "uniquely identify UCEs.")
    parser.add_argument('-u', '--uces', type=argparse.FileType('rU'), required=True,
                        help="An input master UCE file (ID, chr, stop, str, type, gene)")
    parser.add_argument('-a', '--against', type=argparse.FileType('rU'), required=True, nargs='+',
                        help="One or more uncollapsed interval files")
    parser.add_argument('-o', '--output',
                        help="The name of the output file (default: 'recurrent_UCEs.txt')")
    if strInput:
        print ("Given debug argument string: {0}".format(strInput))
        return parser.parse_args(strInput.split())
    return parser.parse_args()


def getBool(sPrompt, iDefault=None):
    while True:
        try:
            sLine = raw_input(sPrompt)
            # If no user entry and default value exists, return default
            if not sLine and iDefault is not None:
                return iDefault
            sLine = sLine.strip().lower()  # case-insensitive
            if sLine == 'y':
                return True
            elif sLine == 'n':
                return False
            else:
                raise ValueError
        except (ValueError, EOFError):
            print("Please enter [Y/n]")


def uceDict(fileobj, nFiles):
    hUCE = {}
    for line in fileobj:
        uce = line.strip().split('\t')
        ID = uce[0]
        interval = formatInt(uce[1:4])
        _type = uce[4]
        try:
            gene = uce[5]
        except IndexError:
            gene = ""
        # Put uces into dict, with ID as key and coords, type, and gene as value
        aVal = [interval, "{0}\t{1}".format(_type, gene)]
        for j in range(nFiles):
            count = 0
            aVal.append(count)
        hUCE[ID] = aVal
    return hUCE


def formatInt(aInterval):
    """ Format an 3-column interval correctly """
    return [aInterval[0], int(aInterval[1]), int(aInterval[2])]


def check(interval, hUCEs, n):
    col = n + 2  # Columns are from [2] onwards as [0] and [1] are UCE coords and type, respectively
    for UCE_ID in hUCEs:
        uce_interval = hUCEs[UCE_ID][0]  # Get interval corresponding to UCE ID
        if overlap(uce_interval, interval):
            hUCEs[UCE_ID][col] += 1  # Increment count
    return


def overlap(intervalA, intervalB):
    # Unpack intervals into chr, start, stop
    chrA, startA, stopA = intervalA
    chrB, startB, stopB = intervalB
    if chrA == chrB:
        if startB <= startA <= stopB:  # Check if interval A starts within interval B
            return True
        elif startB <= stopA <= stopB:  # Check if interval A ends within interval B
            return True
        elif startA <= startB and stopA >= stopB:  # Check if interval A surrounds interval B
            return True
        else:  # If not, then intervals don't overlap
            return False
    return False


def overwriteCheck(filename):
    if os.path.isfile(filename):
        if getBool("{0} already exists in current directory, overwrite? [Y/n]: ".format(filename)):
            return True
        else:
            return False
    return True


def write(hDict, outputArg, againstFiles):
    if outputArg:
        filename = outputArg
    else:
        filename = "recurrent_UCEs.txt"
    againstHeader = " Count\t".join(againstFiles)
    if overwriteCheck(filename):
        with open(filename, 'w') as fh:
            fh.write("UCE_ID\tChr\tStart\tStop\tType\tGene\t{0}\n".format(againstHeader))
            for key, val in hDict.items():
                fh.write("{0}\t{3}\t{4}\t{5}\t{1}\t{2}\n".format(key, val[1], "\t".join(map(str, val[2:])), *val[0]))
            print ("\nWrote to {0}".format(filename))
    else:
        print ("\nDid not write to {0}, exiting...".format(filename))
        sys.exit(1)


def main(args):
    nFiles = len(args.against)
    hUCEs = uceDict(args.uces, nFiles)
    for n, varfile in enumerate(args.against):
        print ("Checking for UCE reoccurence in {0}".format(os.path.basename(varfile.name)))
        for lino, line in enumerate(varfile, 1):
            print ("Checking interval {0}".format(lino), end='\r')
            interval = formatInt(line.strip().split('\t'))
            check(interval, hUCEs, n)
        print ("", end='\n')
    write(hUCEs, args.output, [os.path.basename(infile.name) for infile in args.against])


if __name__ == "__main__":
    args = getArgs()
    main(args)
