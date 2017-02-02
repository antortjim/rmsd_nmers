from protein_user_defined import *
#import Bio.PDB.Superimposer as Superimposer
#import pickle
from sys import argv

# Set fragment length
filename, flen = argv

flen = int(flen)

print "%d-mers chosen" % flen

out_folder = '../shiny/out/'
db = '../../top/'

my_sequences = parse_sequence(db,)  # Parse sequences (fasta sequence)
my_structures = parse_structure(db) # Parse chains (CaPPBuild)

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



# Compute rmsd from matching pairs
my_RMSD = signal(my_fragments)
l = len(my_RMSD)

# Select l random pairs, center them and compute rmsd
my_rRMSD = random_pairs(my_fragments, l)


# Export data to files for R posterior processing
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
