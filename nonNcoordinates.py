#!/usr/bin/env python

"""
Given a FASTA file, returns intervals of stretches that do not have any Ns. When
given the argument -n, it will also skip over lower-case bases, which usually
represent repeat-masked regions. This script can operate on two different types
of FASTA files. Given the -g argument, it expects a FASTA file with one entry
per chromosome, and assumes the start coordinate is 0. Otherwise, it will expect
multiple entries per chromosome and parse them accordingly, reading in the start
and end coordinates from the header.

Requirement: BioPython to parse FASTA files

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
try:
    import argparse
    import re
    import sys
    from Bio import SeqIO
    from Bio.Alphabet import IUPAC
except ImportError as err:
    print "Missing a required module: {0}".format(err)
    sys.exit("Exiting...")

# Catch command-line arguments 
parser = argparse.ArgumentParser(description="Returns coordinates of non-N stretches in the input FASTA file in "
                                             "interval format (1-based starts).")
parser.add_argument("input", type=argparse.FileType("rU"),
                    help="FASTA file")
parser.add_argument("-n", "--nonrep", action="store_true",
                    help="Convert repeat-masked bases to N and remove those regions from output intervals.")
parser.add_argument("-g", "--genomic", action="store_true",
                    help="Expect one FASTA entry per chromosme and header that is only the chromosome number "
                         "(ex. >chr1)")


# Compile patterns for sequence parsing
repmasked_bases = re.compile(r'[actg]')
unsequenced_bases = re.compile(r'[Nn]')


def single_FASTA_parser(FileIn, nonrep):
    # Write message before looping through FASTA file
    if nonrep:
        sys.stderr.write("Removing repeat-masked elements..." + "\n")
        # Loop over each FASTA entry
    for seq_record in SeqIO.parse(FileIn, "fasta",
                                  alphabet=IUPAC.unambiguous_dna):
        # Parse seq_record to give the chromosome name and sequence as a string
        # Expects that FASTA entry header is only the chromosome line (ex. chr1)
        strChr = seq_record.id
        strSeq = str(seq_record.seq)

        if nonrep:
            # Mask lower-case bases to Ns
            strSeq = repmasked_bases.sub("N", strSeq)

        # Reset counters and lists, assumes FASTA entry starts at the beginning
        # of chromosome
        abs_count = 0
        left = 0
        aCoordinates = []
        # Count bases, excluding Ns
        for base in strSeq:
            abs_count += 1
            if unsequenced_bases.match(str(base)):
                if left:
                    # Subtract 1 because abs_count is at an N
                    right = (abs_count - 1)
                    aCoordinates.append("{0}\t{1}\t{2}".format(strChr, left, right))
                    # Clear line and left/right counts
                    left = 0
                else:
                    pass
            else:
                if not left:
                    # Set left boundary to first non-N base
                    left = abs_count
                else:
                    pass
        if left:
            # If interval is open but not closed and parser has run through all
            # bases for that entry, close interval and append to output
            right = abs_count
            aCoordinates.append("{0}\t{1}\t{2}".format(strChr, left, right))
            # Write intervals to stdout
        stdout_writer(aCoordinates)


def multi_FASTA_parser(FileIn, nonrep):
    # Write message before looping through FASTA file
    if nonrep:
        sys.stdout.write("Removing repeat-masked elements..." + "\n")
        # Loop over each FASTA entry
    for seq_record in SeqIO.parse(FileIn, "fasta",
                                  alphabet=IUPAC.unambiguous_dna):
        # Parse seq_record to give the chromosome name, start, stop and
        # sequence as strings
        astrSeqRecord = str(seq_record.id).strip().split("_")
        strChr = astrSeqRecord[1]
        strRecStart = astrSeqRecord[2]
        strSeq = str(seq_record.seq)

        if nonrep:
            # Mask lower-case bases to Ns
            strSeq = repmasked_bases.sub("N", strSeq)

        # Reset counters and lists, assign interval start and end
        # Adjust abs_count so that it correctly counts bases
        abs_count = int(strRecStart) - 1
        left = 0
        aCoordinates = []
        # Count bases, excluding Ns
        for base in strSeq:
            abs_count += 1
            if unsequenced_bases.match(str(base)):
                if left:
                    # Subtract 1 because abs_count is at an N
                    right = (abs_count - 1)
                    aCoordinates.append("{0}\t{1}\t{2}".format(strChr, left, right))
                    # Clear line and left/right counts
                    left = 0
                else:
                    pass
            else:
                if not left:
                    # Set left boundary to first non-N base
                    left = abs_count
                else:
                    pass
        if left:
            # If interval is open but not closed and parser has run through all
            # bases for that entry, close interval and append to output
            right = abs_count
            aCoordinates.append("{0}\t{1}\t{2}".format(strChr, left, right))
            # Write intervals to stdout
        stdout_writer(aCoordinates)


def stdout_writer(aList):
    print '\n'.join(aList)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.genomic:
        single_FASTA_parser(args.input, args.nonrep)
    else:
        multi_FASTA_parser(args.input, args.nonrep)
