#!/usr/bin/env python
"""
Sum coordinates

Given one or more 1-based interval files, this script sums up the coordinates covered by
the intervals and returns the number of lines and number of bases covered.
Collapses file to prevent double-counting of intervals.

Chamith Fonseka
10 April 2013
"""
import argparse
import sys


def getArgs(strInput=None):
    parser = argparse.ArgumentParser(description="Returns total coverage and number of lines of one or more"
                                                 " 1-based interval files. Collapses overlapping intervals to avoid "
                                                 "double-counting")
    parser.add_argument('-u','--uncollapse', action='store_true',
                        help="Sum coordinates of given files without collapsing overlapping intervals")
    parser.add_argument("file", type=argparse.FileType('rU'), nargs='+',
                        help="One or more 3-column interval files")
    if strInput:
        print "Given debug argument string: {0}".format(strInput)
        return parser.parse_args(strInput.split())
    return parser.parse_args()


def counter(aIntervals):
    """ Count bases in every interval and sum across them """
    iBaseCoverage = iIntervalCount = 0
    for interval in aIntervals:
        # Add one to correct 0-based count
        iBaseCoverage += (interval[2] - interval[1] + 1)
        iIntervalCount += 1
    return iBaseCoverage, iIntervalCount


def collapse(aIntervals):
    # Initialize variables
    strChr = iStart = iStop = 0
    aOut = []
    for aInterval in aIntervals:
        # Test if an interval has been stored (always past first loop)
        if strChr:
            # Test if next interval is on a different chr OR if start of
            # next interval is larger than stop of previous interval 
            if strChr != aInterval[0] or int(aInterval[1]) > (iStop + 1):
                # Write interval and clear start/stop coordinates
                aOut.append([strChr, iStart, iStop])
                iStart = iStop = 0
        strChr = aInterval[0]
        # Advance to next interval if iStart is empty
        if not iStart:
            iStart = int(aInterval[1])
            # If next interval overlaps, adjust stop to larger coordinate
        if int(aInterval[2]) > iStop:
            iStop = int(aInterval[2])
        # Write last line
    aOut.append([strChr, iStart, iStop])
    return aOut


if __name__ == "__main__":
    args = getArgs()
    if not args.uncollapse:
        sys.stderr.write('Collapsing overlapping intervals...\n')
    for inFile in args.file:
        try:
            aIntervals = map(lambda x: [x[0], int(x[1]), int(x[2])],
                         [line.strip().split('\t') for line in inFile if line.strip()])
            aIntervals.sort(key=lambda x: (x[0], x[1], x[2]))
        except IndexError, ValueError:
            print "Unable to parse lines in {0}, exiting...".format(inFile.name)
            sys.exit(1)
        if not args.uncollapse:
            aCollapsed = collapse(aIntervals)
        else:
            aCollapsed = list(aIntervals)
        iCoverage, iCount = counter(aCollapsed)
        print "{0}\t{1}\t{2}".format(inFile.name, iCount, iCoverage)