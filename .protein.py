import os # get list of files in directory
from Bio.Alphabet import generic_protein #rebuilding Seq objects
from Bio.PDB.Polypeptide import *
import Bio
import Bio.PDB as PDB
from Bio.PDB.PDBParser import PDBParser # parsing PDB file
from collections import defaultdict
from sys import argv
import Bio.SeqIO as SeqIO # import read/write fasta functions
from Bio.SeqRecord import SeqRecord
import itertools # combine pairs of fragments when more than 2
                 # instances of the same fragment are found
import numpy as np

filename, folder = argv


def clean_seq(Seq_object):
    '''Solving parse(, pdb-seqres) issue that puts X
    every 13 residues when using legacy pdb files
    SEQRES fields must have 13 fields (old pdb files
    put the pdb id in field 14'''

    seq = Seq_object.seq
    alphabet = Seq_object.alphabet
    seq = seq.replace("X", "")
    
    new = Seq(seq, alphabet)
    return(new)
    

flen = 5

# Instantiate a CA builder with radius 4
ppb=CaPPBuilder(radius = 4)

k = 0
my_sequences = []
my_polypeptides = {}

for key, s in my_structures.items(): # Iterate through parsed structures
    # Build polypeptide list
    pp_list = ppb.build_peptides(s, aa_only = 1)
    my_polypeptides[key] = pp_list
    
    # Extract sequences for every polypeptide
    # (Usually just 1 pp, so 1 seq, but sometimes there are > 1 pps)
    # Generate a list where every element is the sequence of the structure
    # (Usually length 1)
    seq = map(lambda pp: str(pp.get_sequence()), pp_list)

    if len(seq) > 1:
        print key
        print seq
    # If there's more than 1 seq because there are gaps, join them with X
    # and return a single string
    # Rebuild the string object
    seq = SeqRecord(Seq("X".join(seq), generic_protein),
              id = my_files[k],
              description = "")
    #for pp in ppb.build_peptides(s): # Build polypeptide
    #    seq = pp.get_sequence()      # Extract seq
    #    my_sequences.append(seq)     # Append to out list
    my_sequences.append(seq)
    k += 1

print my_sequences[0]
print my_sequences[9]
#print len(my_files)
#print len(my_sequences)



## Export built sequences to multifasta_file
with open(folder + ".fasta", "w") as output_handle:
    SeqIO.write(my_sequences, output_handle, "fasta")


# Initalize fragments dict
my_fragments = {}
# Generate fragments of length n using list comprehension for all seqs
n = 5
k = 0
for seq in my_sequences:
    # save seq_id
    seq_id = my_files[k]
    # save fragment sequence and starting position in seq as a tuple
    # exclude fragments containing X
    fragments = [(seq[i:i+n], i) for i in range(0, len(seq) - n) if not "X" in seq[i:i+n]]
    
    for f in fragments: # for every tuple (f_seq, f_start)
        ## add it to the entry corresponding to that fragment
        ## if the entry does not exist, create it and store a list of length 1
        ## if the entry exists, add the tuple to the already existing list
        
        my_fragments.setdefault(str(f[0].seq), []).append((seq_id, f[1]))

    k += 1

# my_fragments is a dictionary counting how many times and where
# the key appears in the set of sequences

## export fragments to long csv data
handle = open("fragments.csv", "w+")
for f, value in my_fragments.items():
    for t in value:
        line = "%s,%s,%d\n" % (f, t[0], t[1])
        handle.write(line)


handle.close()




def extract_CA(seq_id, data, flen):

   seq = my_sequences[my_files.index(seq_id)]
   

def extract_X(key, data):

    '''Receives a tuple saying in which sequence the
     fragment occurs, and in which residue'''
    seq_id = data[0]
    print seq_id
    print key
    print data[1]
    
    X = []
    for i in range(5):
        # extract CA vector for residue from query_index to query_index + 4
        # (5 in total)
        print i
        r = residues[data[1] + i]
        print r
        vector = r["CA"].get_coord()
         
        X.append(vector)

    X = np.asmatrix(np.vstack(X))
    return X


def center_pair(pair):
    '''Receives a pair of coordinates matrices
    and returns a pair where one the matrices
    hass been moved so that both share the same center'''
    # extract X and Y
    X, Y = pair
 
    # center X and Y
    X_center = np.mean(X, axis = 0)
    Y_center = np.mean(Y, axis = 0)
    difference = Y_center - X_center
    Y = Y + difference
    return [X, Y]

def perform_SVD(pair):
     '''Receives a centered pair of matrices
     and returns U, Y and RMSD by singular value decomposition'''
     X, Y = pair

     n = X.shape[0]

     E0 = sum(np.linalg.norm(X, axis = 1) ^ 2 + np.linalg.norm(Y, axis = 1) ^ 2)

     R = Y * X.t 

     V, S, Wt = SVD(R)

     Z = np.diag([1, 1, -1])

     U = Wt.t * V.t

     if np.linalg.det(U) == 1:
         U = W * Z * V.t

     # Rotate Y using U
     Yr = U * Y
     
     Eu = E0 - 2 * sum(S) 
     RMSD = np.sqrt( 1/n *  Eu )
 
     result = {"U" : U, "Y" : Yr, "RMSD" : RMSD}
     return result
     
 


    
# access every pair, extract the CA vectors and perform SVD
for fragment, value in my_fragments.items():
    if len(value) > 1:
        match_n = len(value)
        matrices = []
        for match in value:
            m = extract_X(fragment, match)
            matrices.append(m)

        pairwise_combinations = itertools.combinations(matrices, 2)
        for pair in pairwise_combinations:
            pair = center_pair(pair)
            pair_svd = perform_SVD(pair)
            print pair_svd 

    

#structure = parser.get_structure("foo", "top100H/1aacH")
#for pp in ppb.build_peptides(structure):
#    pp.get_sequence()



