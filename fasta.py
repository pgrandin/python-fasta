#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A FASTA parser

This code implements a FASTA parser, combining the sequences from an input
file into the unique combination of all the sequences.

"""

__author__ = "Pierre Grandin"
__email__ = "linkedin@kazer.org"

import getopt
import sys
import logging
from fragment import Sequence

try:
    OPTS, ARGS = getopt.getopt(sys.argv[1:], 'i:d', ['inputfile=','debug'])
except getopt.GetoptError:
    print 'Usage: driver.py -i <inputfile> [--debug]'
    exit(2)


verbosity = logging.INFO

for opt, arg in OPTS:
    if opt in ("-i", "--inputfile"):
        inputfile = arg
    if opt in ("-d", "--debug"):
        verbosity = logging.DEBUG

logging.basicConfig(level=verbosity,
                    format='%(asctime)s %(levelname)s %(message)s') # %(name)s
LOGGER = logging.getLogger(__name__)
logging.getLogger("requests").setLevel(logging.WARNING)

def insert_fragment_in_place(sequences, fragment):
    """ Inserts a given fragment in the correct sequence in the array
    of sequences if possible, or create a new sequence and add it
    to the list of sequences if required
    """
    # First, check if the fragment fits in any of the existing sequences
    inserted = False
    for i in range(0, len(sequences)):
        sequence = sequences[i]
        LOGGER.debug("Looking at sequence #{} of size {}".format(i, sequence.length()))
        inserted = sequence.insert_if_overlaps(fragment)
        if inserted:
            # If the fragment was inserted in the sequence, we can move on
            break

    # If we were unable to insert in an existing sequence, create a new one
    if inserted == False:
        LOGGER.debug("Unable to fit this fragment in an existing sequence. Creating a new sequence")
        sequence = Sequence()
        sequence.append(fragment)
        sequences.append(sequence)
        LOGGER.debug("{} sequences now in the array".format(len(sequences)))
    else:
        LOGGER.debug("{} fragments now in the sequence".format(sequence.length()))

# We will store the sequences of fragments in an array
sequences = []
# We need to rebuild each fragment from the multi-line
# data from the FASTA file

current_fragment = ""
with open(inputfile, "r") as ins:
    for line in ins:
        if not line.startswith('>'):
            sanitized = line.rstrip()
            current_fragment = current_fragment + sanitized
            if len(sanitized) < 60:
                insert_fragment_in_place(sequences, current_fragment)
                current_fragment = ""


# At this point we have a bunch of sequences in the array, each containing
# 1 or more fragments. We need to parse the fragments again, now that we have
# loaded all the datas, to re-assemble the full sequence if possible.

# To do this, we only need to compare the first and last fragment of each sequence
# with the opposite fragment of another sequence to check if they match.

# This code assumes that the data have been deduplicated in the input file
# Deduplication would be easy to add at loading time if needed
previous_size = len(sequences)
while len(sequences) > 1:
    LOGGER.debug("I have {} sequences to assemble".format(len(sequences)))
    for i in range(0, len(sequences)-1):
        # Because we alter the length of the list of sequences,
        # we need to check if we haven't reached the end of the list after
        # it has been modified, as the range() will not be updated during the
        # processing
        if i > len(sequences)-1:
            break
        sequence = sequences[i]
        for j in range(i+1, len(sequences)):
            other = sequences[j]
            merged_list = sequence.merge_if_overlaps(other)
            if merged_list != None:
                sequences.remove(other)
                sequences.remove(sequence)
                sequences.append(merged_list)
                break
            else:
                LOGGER.debug("No match found between lists {} and {}".format(i,j))

    if previous_size == len(sequences):
        # We did not find any sequences to merge.
        # Effectively, we are stuck. It is possible that we are missing
        # a fragment to reconstruct one unique sequence
        LOGGER.error("No merge done, I'm stuck! We are probably missing a fragment")
        exit(2)
    else:
        previous_size = len(sequences)

# We have only one sequence left in the array, let's print it
print sequences[0].flatten(max_width=60)


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
