from Bio.PDB.PDBParser import PDBParser      # parsing PDB file
from Bio.PDB.Polypeptide import *            # build polypeptides
import os                                    # get list of files in dir
from Bio.Alphabet import generic_protein     # rebuilding Seq objects
import Bio
import Bio.PDB as PDB
from Bio.PDB.PDBParser import PDBParser      # parsing PDB file
#from collections import defaultdict
from sys import argv
import Bio.SeqIO as SeqIO                    # import read/write
                                             # fasta functions
from Bio.SeqRecord import SeqRecord          
import itertools                             # combine pairs of fragments
                                             # when more than 2 instances
                                             # of same fragm seq are found
import numpy as np                           # mathematical computations
from Bio.SeqUtils import IUPACData  as iupac # Required in real aminoacid
from Bio.PDB.Polypeptide import three_to_one # checks
import random # for selecting random pairs


aa_dict = iupac.protein_letters_3to1
aa_list = map(lambda aa: aa.upper(), aa_dict.keys())

def is_aa(residue):
    resname = residue.get_resname()
    return True if resname in aa_list and residue.has_id('CA') else False
    

def clean_seq(Seq_object):
    '''Solving parse(, pdb-seqres) issue that puts X
    every 13 residues when using legacy pdb files
    SEQRES fields must have 13 fields (old pdb files
    put the pdb id in field 14'''

    seq = str(Seq_object)
    alphabet = Seq_object.alphabet
    seq = seq.replace("X", "")

    new = Seq(seq, alphabet)
    return(seq)


def parse_structure(folder):
    '''Receives a folder with files whose filenames
    refer to pdb files in top100H folder.
    Returns a dictionary of structures
    where every key gives the structure id.'''
    print "Parsing structures"

    # Instantiate database
    my_files = sorted(os.listdir(folder))
    
    parser = PDBParser()

    # Initialize structures dict 
    my_structures = {}

    # Iterate through my_files and get_structure "foo" from file f
    nstr = 0
    for f in my_files:
        structure = parser.get_structure(f, "top/" + f)
        my_structures[f] = structure # Store structure with key given by fname
        nstr += 1
        #c_id = 0
        #selected_chains = []
        #for c in structure[0]:
        #    if check_chain(c, f, c_id):         
        #        selected_chains.append(c)
        #    c_id += 1
        #refined_structure = iter(selected_chains)

        #my_structures[f] = refined_structure # Store structure with key given by fname

    print "Parsed %d structures" % nstr
    return my_structures
  

def get_visible_chains(SeqRecord_list, structure):

    chains = list(structure[0])
    chains_id = []
    for c in chains:
        chains_id.append(c.get_id())

    visible = []
    for SeqRecord in SeqRecord_list:
        if SeqRecord.id in chains_id:
           visible.append(SeqRecord.id)

    return visible

def remove_invisible_seqs(rr_list, structure, key):
    '''Clean entries in SeqRecord_list that have
    no coordinates info'''
    seqs = []

    visible = get_visible_chains(rr_list, structure)
    newlist = []
    my_max = len(rr_list)
    j = 0
    i = 0
    while i + j < my_max:
        record = rr_list[j]
        # Extract id
        seq_id = record.id
        # If the id is not in the visible
        if not seq_id in visible:
            rr_list.pop(j)
            i += 1
        else:
            j += 1
 

    return rr_list 

       
 
def parse_sequence(folder):
    '''Receives a folder with files whose filenames
    refer to pdb files in top100H folder.
    A my_structures dictionary is supplied
    to select sequences that have their structural counterpart.
    Returns a dictionary of lists
    where every entry has a list of Seq objects as
    parsed by SeqIO.parse(file, "pdb-seqres")'''

    print "Parsing sequences from SEQRES in PDB file"
    acum = 0

    # Instantiate database
    my_files = sorted(os.listdir(folder))

    # Initalize sequences dict
    my_sequences = {}

    # For every file parse the sequences saved in SEQRES field
    # Debug legacy error that puts X after every 13 res
    for f in my_files:
        handle = open("top/" + f, "rU")
        clean_records = []
        for record in SeqIO.parse(handle, "pdb-seqres"):
            #record.seq = clean_seq(record.seq)
            clean_records.append(record)

        handle.close()
        my_sequences[f] = clean_records # Store Seq objects list
                                        # with key given by fname
        acum += 1

    print "Parsed %d lists of SeqRecord objects" % acum
    return my_sequences

## Refine sequences so that only sequences with coordinates are available
def refine_sequences(my_sequences, my_structures):
    '''Removes SeqRecords that have no spatial counterpart'''
    for key, value in my_sequences.items():
        new = remove_invisible_seqs(value, my_structures[key], key)
        my_sequences[key] = new

    return my_sequences

        
             

def build_polypeptide(s, r):
    '''Receives structure object and radius.
    Returns a polypeptide list for the structure object
    built using CaPPBuilder with radius r'''
    ppb = CaPPBuilder(radius = r)
    pp_list = ppb.build_peptides(s, aa_only = 1)

    # Process pp_list so that you get rid of heteroatoms
    
    return pp_list


def rewrite_seq(seq, my_chain, key):
    #print "Processing sequence in %s" % key
    i = 0 # Iterating throuh letters in seq coming from my_sequences
    j = 0 # Iterating through residues in my_chain
    my_chain = list(my_chain)
    final_seq = ""
    changed = 0
    while j < len(my_chain) and i < len(seq):
        if is_aa(my_chain[j]):
            r1 = three_to_one(my_chain[j].get_resname())
            r2 = seq[i]
            if r1 == r2:
               final_seq = final_seq + r1
               j += 1   # move to the next res in my_chain
            else:
               final_seq = final_seq + "X"
               changed = 1
               # don't move to the next residue in my_chain
               # seq has to catch up because there is a gap
               
            i += 1 # move to the next residue in seq
        else:
            j += 1 # move on anyway
        
    
    return [final_seq, changed]



def integrate_seq_str(my_sequences, my_structures):
    '''Refines my_sequences so that residues for which there is no
    coordinates info are marked with X'''
    #fragments = []  # UNCOMMENT ME IF THERE IS STH WRONG
    np_seq = 0
    for key, value in my_sequences.items():
        s = my_structures[key]
        c_ids = []
        for c in s[0]:
            c_ids.append(c.get_id())
         
        if len(c_ids) == 1:
            record = value[0]
            seq = str(record.seq)
            seq_id = record.id
            my_chain = list(s[0])[0]
            final_seq, changed = rewrite_seq(seq, my_chain, key)
            my_sequences[key] = [SeqRecord(Seq(final_seq, generic_protein),
                                             id = seq_id)]

            
        else:
            k = 0
            for record in value:
                seq = str(record.seq)
                seq_id = record.id
                c_idx = c_ids.index(seq_id)
                # Select chain with id equal to seq_id
                my_chain = list(s[0])[c_idx]
                i = 0
                final_seq, changed = rewrite_seq(seq, my_chain, key)
                my_sequences[key][k] = SeqRecord(Seq(final_seq, generic_protein),
                                                 id = seq_id)
                k += 1

        np_seq += changed
                    

    print "Processed %d sequences" % np_seq
    return my_sequences 


class Fragment(object):

    def __init__(self, seq, start, seq_id, chain_id, X, atoms):
        self.seq = seq
        self.start = start 
        self.seq_id = seq_id 
        self.chain_id = chain_id
        self.X = X 
        self.atoms = atoms

    def center(self):
#        if len(self.X) > 1:
#            print self.X
#            print type(self.X)
#            raise Exception("More than 1 matrix for fragment")
        return np.mean(self.X, axis = 0)

    def set_X(self, X):
        #del self.X
        self.X = X

    
def generate_fragments(my_sequences, my_structures, flen):

    my_fragments = {} # Initialize fragments dictionary
                      # fragment sequences will be used as keys
                      # value will be a list of fragment objects
                      # sharing the same 5 res seq with key
    acum2 = 0
    acum3 = 0
    nfragments = 0
    for key, value in my_sequences.items():    # for every list of SeqRecords
                                               # in each pdb file
        
        k = 0                                  # iterate in SeqRecords list
                                               # k = the index of the current
                                               # record in the SR list
        local = []                             # list that will store fragments
                                               # extracted from a given SR list
        for record in value:                   # for every SeqRecord
                                               # in the current list
            c = list(my_structures[key][0])[k] # extract the corresponding chain
                                               # from my_structures

            chain_id = c.get_id()              # extract chain id
            seq = str(record.seq)              # extract SRecord sequence


            ## List comprehension
            ## Put every length 5 substring in seq that does not contain X
            ## as an element of local_f
            ## Local_f stores whole fragments extracted from seq coming from
            ## value[k]
            
            local_f = [(seq[i:i + flen], i) for i in range(0, len(seq) - flen) if not "X" in seq[i:i+flen]]
            
            # Generate a fragment object for every subseq in local_f
            # Store 5 res subseq, starting pos in original seq, whole seqid,
            # chain id and CA_matrix.
            m = 0
            for f in local_f:
                reject_f = False
                my_matrix = []
                atoms = []
                residues = list(c)[m:m+flen]

                  
                acum = 0
                res_count = 0
                for r in residues:
                    if not is_aa(r):
                        reject_f = True 
                        break
                    else:
                        #print f[0][res_count], r.get_resname() checks that t
                        a = r['CA']
                        atoms.append(a)
                        coords = a.get_coord()
                        if len(coords) == 3:
                            acum += 1
                        my_matrix.append(coords)
    
                        CA_matrix = np.asmatrix(my_matrix)
                        res_count += 1
                if reject_f == True:
                    break

                if acum < flen:
                    #print list(residues)
                    acum2 += 1
               
                acum3 += 1        
                if CA_matrix.shape == (flen, 3): # in principle it should always 
                                              # at least in the exam
                    fragment = Fragment(f[0], f[1], key,            # Instantiate
                                        chain_id, CA_matrix, atoms) # Fragment class
                    nfragments += 1
                    my_fragments.setdefault(f[0], []).append(fragment)
                    m += 1 # move on to the next residue



            k += 1 # move on to next record in value (SR list)

    #print acum2
    #print acum3
    print "Generated %d fragments" % nfragments
    return my_fragments


def focus(f1, f2):
    '''Receives 2 Fragment objects and centers them around 0'''
    f1_center = f1.center()
    f2_center = f2.center()

    f1.set_X(f1.X - f1_center)
    f2.set_X(f2.X - f2_center)


    result = [f1, f2, f1_center, f2_center]
    return result


def rotate_pair(pair):
    '''Receives a centered pair of fragments and returns Y, U and S'''

    # Put them in 0, 0, 0
    f1, f2, f1_center, f2_center = focus(pair[0], pair[1])

    # Received as 5 x 3
    X = f1.X
    Y = f2.X

    R = Y.T * X # 3 x 5 * 5 x 3 = 3 x 3

    V, S, Wt = np.linalg.svd(R)

    #print "R"
    #print R
    #print "V"
    #print V
    #print "Wt"
    #print Wt
 
    #Z = np.diag([1, 1, -1])

    U = Wt.T * V.T

    if np.linalg.det(U) < 0:
        #print "Reflection catch"
        Wt[2] = -Wt[2]
        U = Wt.T * V.T

    # Rotate Y using U
    Yr = (Y * U.T) # + f1_center  # 5 x 3 * 3 x 3 = 5 x 3

    # Set f1 matrix to orignal place
    #f1.set_X(f1.X + f1_center)

    result = [X, Yr, U.T, S]
    return result

def compute_RMSD(pair, S):

    # Received as n x 3
    f1 = pair[0]
    f2 = pair[1]


    X = f1.X
    Y = f2.X

    # Transpose to 3 x 5
    X = X.T
    Y = Y.T
   
    #print "Tras pasar por rotate"
    #print f1.center() 
    #print f2.center() 

    n = float(X.shape[1])

#    sq_diff = np.square(X - Y) # 3 x 5
#
#    my_sum = np.sum(sq_diff, axis = 0) # suma por columnas
#                                       # becomes 1 x 5
#    print my_sum                         
#    rmsd = np.sqrt(1 / n * my_sum)
#    RMSD = np.sum(rmsd, axis = 1)[0, 0] # sum all rmsd values

    E0 = np.sum(np.square(np.linalg.norm(X, axis = 0)) + np.square(np.linalg.norm(Y, axis = 0)))
#    print E0

    RMSD = np.sqrt((1 / n) * (E0 - 2 * np.sum(S)))

    return RMSD

def signal(my_fragments):
    
    print "Working on fragment matches"
    my_RMSD = []
    #i = 0
    for key, value in my_fragments.items():
        pairwise_combinations = itertools.combinations(value, 2)
        for c in pairwise_combinations:
            f1 = c[0]
            f2 = c[1]
    #        before = f2.X
            X, Y, U, S = rotate_pair([f1, f2])
            f1.set_X(X)
            f2.set_X(Y)
    #        after = f2.X
            RMSD = compute_RMSD([f1, f2], S)

    #        if i < 5:
    #            print f1.X
    #            print after
    #            print RMSD
    #            i += 1
            
            result = (c[0].seq, c[0].seq_id, c[0].chain_id, c[0].start,
                      c[1].seq, c[1].seq_id, c[1].chain_id, c[1].start,
                      RMSD)
            my_RMSD.append(result)


    return my_RMSD


def random_pairs(my_fragments, l):
    
    print "Generating %d random pairs" % l
    my_RMSD = [] * l
    fragments = []
    for value in my_fragments.values():
        for f in value:
            fragments.append(f)

    b = len(fragments) - 1 #  last index in fragments
    #print b
    i = 0    
    while i < l: # generate as many random pairs as
                 # real matches are available

        first = random.randint(0, b) 
        second = random.randint(0, b) 

        f1 = fragments[first]
        f2 = fragments[second]

        X, Y, U, S = rotate_pair([f1, f2])
        f1.set_X(X)
        f2.set_X(Y)
        RMSD = compute_RMSD([f1, f2], S)

        result = (f1.seq, f1.seq_id, f1.chain_id, f1.start,
                  f2.seq, f2.seq_id, f2.chain_id, f2.start,
                  RMSD)
 
        my_RMSD.append(result)
        i += 1
    return my_RMSD
