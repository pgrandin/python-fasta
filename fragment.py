"""A linked link implementation, with extra features to add fragments based upon string comparison

This module implements a singly-linked list, targeted at ordering DNS sequence from fragments,
and providing some methods to add fragments in the correct place in the sequence if they overlap
by at least 50% of their length.

"""

__author__ = "Pierre Grandin"
__email__ = "linkedin@kazer.org"


from suffixtree import SuffixTree

class Fragment(object):
    """ A simple fragment, to hold data from one specific fragment

    This class implements the methods to store information about one speficic fragment:
     - the DNA sequence using a string representation (such as T/C/G/A)
     - a pointer to the next fragment of the sequence
     - the offset at which this fragment overlaps from the previous fragment in the sequence
        (only used when reconstructing the full sequence)
    """

    def __init__(self, data=None, next_fragment=None):
        """ Initializes a new fragment, with optional 'data' and 'next_fragment' parameters """
        self.data = data
        self.next_fragment = next_fragment
        self.offset = -1

    def get_data(self):
        """ Returns the data of the current fragment """
        return self.data

    def get_next(self):
        """ Returns the next fragment """
        return self.next_fragment

    def set_next(self, new_next):
        """ Set the next fragment """
        self.next_fragment = new_next

    def set_offset(self, offset):
        """ Sets the offset of the current fragment """
        self.offset = offset

    def get_offset(self):
        """ Gets the offset of the current fragment """
        return self.offset

    def find_overlap(self, data):
        """ Determines if the data passed as argument overlap
        the data from the current fragment.

        If so, returns a dict with the offsets for both the current
        fragment and the data used for the comparison. This is useful to
        correctly set the offset when adding the fragment from the Sequence
        class

        """
        suffix_tree = SuffixTree()
        suffix_tree.append_string(self.data)
        suffix_tree.append_string(data)

        lcs = suffix_tree.find_longest_common_substrings()
        if(len(lcs)) == 1:
            substring_length = len(lcs[0])
            if substring_length / float(len(self.data)) > 0.5:
                pos_1 = self.data.find(lcs[0])
                pos_2 = data.find(lcs[0])
                return {'my_index' : pos_1, 'their_index' : pos_2}
        return None

class Sequence(object):
    """ A simple singly-linked list class to store and manipulate
    representations of the fragments.
    """
    def __init__(self, first_fragment=None):
        """ Initializes a sequence, with an optional first fragment """
        self.head = first_fragment
        self.tail = first_fragment
        if first_fragment != None:
            self._length = 1
        else:
            self._length = 0

    def get_first(self):
        """ Returns the first fragment of the sequence """
        return self.head

    def get_last(self):
        """ Returns the last fragment of the sequence """
        return self.tail

    def length(self):
        """ Returns the length (number of fragments) of the sequence """
        return self._length

    def append(self, data, offset=-1):
        """ Appends a new fragment to the sequence. """
        new_fragment = Fragment(data)
        new_fragment.set_offset(offset)

        # If we are appending to an empty list, also
        # update the head pointer to point to the new
        # fragment
        if self._length == 0:
            self.head = new_fragment
        else:
            # Add the new fragment after the current tail
            previous_tail = self.tail
            previous_tail.set_next(new_fragment)

        # Move the tail pointer to the new fragment
        self.tail = new_fragment

        # Lastly, increment the length counter
        self._length = self._length + 1

    def insert_if_overlaps(self, data):
        """ Inserts a fragment in the correct position in the sequence.
        The correct position can only be the first or the last fragment
        of the sequence, the only places where we should find an overlap
        """
        current = self.head
        previous = None
        has_at_least_one_overlap = False
        while current:
            if current.get_data() != data:
                overlap = current.find_overlap(data)
                if overlap != None:
                    has_at_least_one_overlap = True
                    new_fragment = Fragment(data)
                    if overlap['their_index'] > 0:
                        # This data goes before the existing set
                        if previous:
                            previous.set_next(new_fragment)
                        new_fragment.set_next(current)
                        new_fragment.set_offset(overlap['their_index'])
                        self.head = new_fragment
                        self._length = self._length + 1
                    elif overlap['my_index'] > 0:
                        new_fragment.set_next(current.get_next())
                        current.set_next(new_fragment)
                        current.set_offset(overlap['my_index'])
                        self._length = self._length + 1
                        self.tail = new_fragment
            # An item might overlap with two existing sequences in the chain
            previous = current
            current = current.get_next()
        return has_at_least_one_overlap

    @staticmethod
    def merge(first, second, offset):
        """ Merge two sequences. Technically this appends the second sequence
        at the end of the first one while preserving all the offsets
        """
        current = second.get_first()
        while current:
            first.append(current.get_data(), current.get_offset())
            current = current.get_next()
        return first

    def merge_if_overlaps(self, other_sequence):
        """ Look for a match between the current sequence and the other
        sequence passed in argument.
        Returns the sequences merged in the correct order if a match is found,
        and None otherwise
        """
        my_head = self.head
        ol_head = other_sequence.get_first()
        my_tail = self.tail
        ol_tail = other_sequence.get_last()

        overlap = my_head.find_overlap(ol_tail.get_data())
        if overlap != None:
            return self.merge(other_sequence, self, overlap['their_index'])

        overlap = my_tail.find_overlap(ol_head.get_data())
        if overlap != None:
            return self.merge(self, other_sequence, overlap['my_index'])

        return None

    def dump(self):
        """ Dumps the first 30 chars of each fragment of the sequence.
        Useful only for debugging.
        """
        current = self.head
        while current:
            print "{}".format(current.get_data()[0:30])
            current = current.get_next()


    def flatten(self, max_width=0, line_termination="\n"):
        """ Flatten the sequence, effectively returning the text representation
        of all the fragments correctly merged together.
        Optional arguments allow to specify the max width of the returned
        string, along with the line termination character(s) (defaults to
        '\n'
        """
        current = self.head
        out = ""
        while current:
            if current.get_offset() > 0:
                out = out + current.get_data()[:current.get_offset()]
            else:
                out = out + current.get_data()
            current = current.get_next()
        if max_width > 0:
            output = ""
            length = len(out)
            for index in range(0, length, max_width):
                output = output + out[index:index+max_width] + line_termination
        else:
            output = out

        return output

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
