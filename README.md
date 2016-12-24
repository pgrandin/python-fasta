fasta.py, a simple FASTA parser written in python.

This code implements a FASTA parser. The example code 'fasta.py' shows how to
leverage the Fragment and List classes to efficiently manage the fragments and
reconstruct the unique sequence based upon all of the input fragments.

Example usage :

$ python fasta.py -i example.txt
ATTAGACCTGCCGGAATAC

# Saving the output :
$ python fasta.py -i example.txt > output.txt

This package contains :
- fasta.py, the example implementation
- the fragment.py module containing the Fragment and List class
- the suffixtree.py module used for efficient string comparison

