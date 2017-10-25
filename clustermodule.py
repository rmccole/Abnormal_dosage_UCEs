#!/usr/bin/env python
'''
Module to implement the clustering feature in randomoverlaps3.py

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

'''

import sys

def cpartial_overlap(aIntervalA, aIntervalB):
    """
    Returns True if interval A overlaps for at least one bp with interval B, them returns the number of overlapping bp

    """
    # Check that both inputs are 3-column intervals
    if not len (aIntervalA) == len(aIntervalB) == 3:
        raise Exception("Regions could not be overlapped")
    # Check that both intervals are on the same chromosome, if not return False
    if aIntervalA[0] == aIntervalB[0]:
        # Unpack coordinates into integer variables
        iIntervalAStart, iIntervalAEnd = aIntervalA[1:]
        iIntervalBStart, iIntervalBEnd =  aIntervalB[1:]
        # Check if start coordinate of interval A lies within interval B
        if cpoint_checker(iIntervalAStart, aIntervalB[1:]):
            return True
        # If not, check if end coordinate of interval A lies within interval B
        elif cpoint_checker(iIntervalAEnd, aIntervalB[1:]):
            return True
        # If not, check if interval A surrounds interval B
        elif iIntervalAStart <= iIntervalBStart and iIntervalAEnd >= iIntervalBEnd:
            return True
        # If not, then intervals do not overlap
        else:
            return False
    else:
        return False

def cpoint_checker(iPoint, aInterval):
    ''' Returns True if the given point is larger than interval start and
    smaller than interval end'''
    # Ensure all comparisons are between integers
    iIntervalStart, iIntervalEnd = map(int, aInterval)
    iPoint = int(iPoint)
    if iPoint >= iIntervalStart and iPoint <= iIntervalEnd:
        return True

def ccollapse(aaIntervals):
    # Initialize variables
    strChr = iStart = iStop = 0
    aOut = []
    for aInterval in aaIntervals:
        # Test if an interval has been stored (always past first loop)
        if strChr:
            # Test if next interval is on a different chr OR if start of
            # next interval is larger than stop of previous interval 
            if strChr != aInterval[0] or aInterval[1] > (iStop + 1):
                # Write interval and clear start/stop coordinates
                aOut.append([strChr, iStart, iStop])
                iStart = iStop = 0
        strChr = aInterval[0]
        # Advance to next interval if iStart is empty
        if not iStart:
            iStart = int(aInterval[1])
        # If next interval overlaps, adjust stop to larger coordinate
        if aInterval[2] > iStop:
             iStop = aInterval[2]
    # Write last line
    aOut.append([strChr, iStart, iStop])
    return aOut

def cluster(aaUCEs, iClusterWidth, hChrEnds):
    '''Converts a list of intervals into a list of clusters. The cluster width
    should be given in kb.'''
    # Convert cluter width (given in kb) to bp, divided by 2 to create flanks
    iClusterWidth *= 500
    # Convert all UCEs to clusters and put into a new list
    aaClusters = []
    for aInterval in aaUCEs:
        strChr = str(aInterval[0])
        # Convert coordinates to integers
        iStart, iStop = map(int, aInterval[1:])
        # Extend start coordinate leftwards, adjusting start if necessary
        iStart -= iClusterWidth
        if iStart < 1:
            iStart = 1
        # Extend stop coordinate rightwards, adjusting end if it overruns
        # the end of the chromosome
        iStop += + iClusterWidth
        if hChrEnds.get(strChr) is None:
            raise Exception("Could not find chromosomes correctly. Ensure that "
                            "all files are based on the same genome.")
        else:
            iMaxStop = hChrEnds.get(strChr)
        if iStop > iMaxStop:
            iStop = iMaxStop
        aaClusters.append([strChr, iStart, iStop])
    # Collapse overlapping clusters
    aaClusters.sort(key=lambda x: (x[0], x[1], x[2]))
    aaCollapsedClusters = ccollapse(aaClusters)
    return aaCollapsedClusters
    
def c_trackuces(aaClusters, aaUCEs):
    aaClustersSizes = []
    aaClusters.sort(key=lambda x: (x[0], x[1], x[2]))
    aaUCEs.sort(key=lambda x: (x[0], x[1], x[2]))
    for aCluster in aaClusters:
        aUCESizes = []
        for aUCE in aaUCEs:
            # If UCE is not on same chr as cluster, skip
            if aUCE[0] != aCluster[0]:
                continue
            # if UCE stop coordinate is less than cluster start coordinate, skip
            if int(aUCE[2]) < int(aCluster[1]):
                continue
            # if UCE start coordinate is greater than cluster stop coordinate,
            # no more UCEs can be found, so skip to next cluster
            if int(aUCE[1]) > int(aCluster[2]):
                break
            # Add size of UCEs that overlap to list associated with cluster
            if cpartial_overlap(aUCE, aCluster) is True:
                # Append a tuple containing the length of the UCE and its offset
                # from the cluster start coordinate
                iLenUCE = int(aUCE[2]) - int(aUCE[1])
                iOffsetUCE = int(aUCE[1]) - int(aCluster[1])
                aUCESizes.append((iOffsetUCE, iLenUCE))
        if len(aUCESizes) < 1:
            print "No UCEs in cluster: " + "\t".join(map(str, aCluster))
            return(aCluster)
        aaClustersSizes.append((aCluster, aUCESizes))
    return(aaClustersSizes)

if __name__ == "__main__":
    print("This is a module designed to implement the clustering feature in "
          "the randomoverlaps3.py script. It is not meant to be run "
          "independently.")
    sys.exit("Exiting...")
