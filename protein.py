from user_defined import *
#import Bio.PDB.Superimposer as Superimposer
#import pickle
from sys import argv

# Set fragment length
filename, flen = argv

flen = int(flen)

print "%d-mers chosen" % flen

out_folder = 'shiny/sb_exam/out/'

my_sequences = parse_sequence('top',)  # Parse sequences (fasta sequence)
my_structures = parse_structure('top') # Parse chains (CaPPBuild)

my_sequences = refine_sequences(my_sequences, my_structures)
#Keep sequences that have a chain counterpart


#### Export built sequences to multifasta_file
#fasta_proteins = {}
#for key, SR_list in my_sequences.items():
#    value = []
#    i = 0
#    for SR in SR_list:
#        new_SR = SeqRecord(Seq(str(SR.seq),
#                           generic_protein),
#                           id = key,
#                           description = SR.id)
#        value.append(new_SR)
#        i += 1
#
#
#    fasta_proteins[key] = value
#    
#
#with open('protein_sequences' + ".fasta", "w") as output_handle:
#    for value in fasta_proteins.values():
#        SeqIO.write(value, output_handle, "fasta")
###

                                       
my_sequences = integrate_seq_str(my_sequences, my_structures)
#Put X in chain breaks

my_fragments = generate_fragments(my_sequences, my_structures, flen)
#Generate fragments using the sequences taking into account chain breaks


n_matches = []

for key, value in my_fragments.items():
    i = len(value) 
    n_matches.append(i)

print "No matches %s" % n_matches.count(1)
print "One match %s" % n_matches.count(2)
print "Two matches %s" % n_matches.count(3)
print "Three matches %s" % n_matches.count(4)
print "Four matches %s" % n_matches.count(5)

# Export fragments
handle = open(out_folder + "%d-mer_fragments.csv" % flen, "w+")
handle.write("seq,seq_id,chain_id,start\n")
k = ["seq", "seq_id", "chain_id", "start"]
for key, value in my_fragments.items():
    for fragment in value:
        my_data = [fragment.__dict__.get(entry) for entry in k]
        #seq = fragment.seq
        #seq_id = fragment.seq_id
        #chain_id = fragment.chain_id
        #start = fragment.start
        #my_data = [seq, seq_id, chain_id, start]
        handle.write(','.join(map(str, my_data)) + '\n')

handle.close()


## Testing the SVD-RMSD algorithm
#f1 = my_fragments.values()[62][0]
#X = f1.X
##print "inicio"
##print f1.center()
#f2 = my_fragments.values()[72][0]
#
#
#X = f1.X
#Y = f2.X
#print "X and Y at start" # checked that it's the same for superimposers
#print X
#print Y
#
#X, Yr, U, S = rotate_pair([f1, f2])
#
#
#print "X and Y after centering and rotation"
#print X
#print Yr
#f1.set_X(X)
#f2.set_X(Yr)
#
#print "Is the center now equal?" # checked it's true
#if np.sum(abs(f1.center() - f2.center()), axis = 1)[0, 0] < 0.1:
#    print True
#else:
#    print f1.X - f2.X
#    print False
#
##print "X and Y after setting to the fragments."
##print "Should be just the same as above" # checked
##X = f1.X
##Y = f2.X 
##print X
##print Y
#
#
#
##MM, UU, SS = rotate_pair([f1, f2])
##print "Should be identity matrix" checked
##print UU
## UU = I proofs my algorithm thinks it can't be reduced
#
#rmsd = compute_RMSD([f1, f2], S)
#
#X = f1.X
#Y = f2.X
#
#
##print X
##print Y
##
#print "U"
#print U
#print rmsd
#
###
#sup = Superimposer()
## Specify the atom lists
## 'fixed' and 'moving' are lists of Atom objects
## The moving atoms will be put on the fixed atoms
#fixed = f1.atoms
#moving = f2.atoms
##set_atoms(fixed, moving)
#sup.set_atoms(fixed, moving)
## Print rotation/translation/rmsd
#SU, ST = map(np.asmatrix, sup.rotran)
#Srmsd = sup.rms
## Apply rotation/translation to the moving atoms
#sup.apply(moving)
#
#f2 = my_fragments.values()[72][0]
#SYr = np.dot((f2.X + ST), SU)
#
#print "SU"
#print SU
#print Srmsd


# Saving the objects:
#with open('objs.pickle', 'w') as f:  # Python 3: open(..., 'wb')
#    pickle.dump([X, Y, SYr, Yr, SU, U], f)



# Compute rmsd from matching pairs
my_RMSD = signal(my_fragments)
l = len(my_RMSD)

# Select l random pairs, center them and compute rmsd
my_rRMSD = random_pairs(my_fragments, l)


# Export data
header = 'f1.seq,f1.seq_id,f1.chain_id,f1.start,f2.seq,f2.seq_id,f2.chain_id,f2.start,rmsd\n'
handle = open(out_folder + '%d-mers_rmsd.txt' % flen, 'w+')

handle.write(header)
for entry in my_RMSD:
    handle.write(','.join(map(str, entry)) + '\n')
handle.close()

handle = open(out_folder + '%d-mers_random.txt' % flen, 'w+')
handle.write(header)
for entry in my_rRMSD:
    handle.write(','.join(map(str, entry)) + '\n')

handle.close()
