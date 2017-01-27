from rna_user_defined import *

# exam
seq1_seq = "AUCAGUUCUAGCAGGAGCUGUACUCAGAGACUCGGGAAAUUUUCCCGGAAUUUUACCCGGGUUUUUACGU"
seq1_str = "..(((((((....))))))).....(((((((((((.((..(((...)))..)).)))))))))))...."

seq2_seq = "AUCGGUUCCAGCAGGAACUGUACUCGGGGGCUCGGGAAACCCUCCCGGGGUUUUACCCGGGUUUUUACGU"
seq2_str = "..(((((((....))))))).(((((((..(((((((.....)))))))......)))))))........"

# RNAalifold sequence
ali_seq =  "AUCGGUUCCAGCAGGAACUGUACUCGGGGGCUCGGGAAACCCUCCCGGGGUUUUACCCGGGUUUUUACGU"
ali_str =  "..(((((((....))))))).(((((((..(((((((.....)))))))......)))))))........"


print "Hamming distance between structures"
print hamming_distance(seq1_str, seq2_str)

print "Hamming distance between sequences"
print hamming_distance(seq1_seq, seq2_seq)


print "Base pair distance between structures"
print bp_distance(seq1_str, seq2_str)

print "Alifold results"
print "Hamming distance"
print hamming_distance(seq1_str, ali_str)
print hamming_distance(seq2_str, ali_str)


print "BP distance"
print bp_distance(seq1_str, ali_str)
print bp_distance(seq2_str, ali_str)
