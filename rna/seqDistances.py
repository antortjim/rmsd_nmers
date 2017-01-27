import numpy as np


def hamming_distance(structure_1, structure_2):
    
    acum = 0 # number of mismatches
    for i in range(len(structure_1)): # for every position in structure_1
        if structure_1[i] != structure_2[i]: # if position i in str_1
                                             # is not equal to str_2[i]
            acum += 1                        # sum 1 to acum



    return acum

def check_bifurcation(structure):
    
    acum = 0
    bifurcation_id = 0
    bifurcations = []
    
    # iterate through the structure string and get index and value
    for idx, element in enumerate(structure):
        # if element is ( add 1 to acum.
        if element == "(":
            acum += 1
            # if acum is now 1, it means a new bifurcation has started
            if acum == 1:
                current_start = idx
        # if element is ) rest 1 to acum 
        elif element == ")":
            acum -= 1
            # if acum is now 0, it means the current bifurcation is over
            if acum == 0:
               bifurcations.append((current_start, idx))

    if acum == 0:
        return bifurcations
    elif acum > 0:
        return "Non valid structure. Missing closing base pairs"
    else:
        return "Non valid structure. Too many closing base pairs"

def extract_base_pairs(structure):
    
    bifurcations = check_bifurcation(structure)
    base_pairs = []
    for b in bifurcations:
        b_start = b[0] 
        b_end = b[1] + 1
        substring = structure[b_start:b_end]
        reverse = substring[::-1]
        acum = 0
       
        #opening = []
        #opening_rank = []
        #closing = []
        #closing_rank = []
      
        opening = {}
        closing = {}

        for idx, c in enumerate(substring):
            idx = idx + b_start # start to count from start of complete structure, not the substring
            if c == "(":
                acum += 1

                #opening.append(idx)
                #opening_rank.append(acum)
                
                try:
                    opening[str(acum)].append(idx)

                except KeyError:
                    opening[str(acum)] = [idx]

            if c == ")":
                acum -= 1
                #closing.append(idx)
                #closing_rank.append(acum + 1)

                try:
                    closing[str(acum + 1)].append(idx)

                except KeyError:
                    closing[str(acum + 1)] = [idx]


        bp = [] 
        for key in opening.keys(): 
            value = opening[key]
            for v in value:
                end = closing[key].pop(0)
                base_pairs.append((v, end))
                 
        base_pairs = base_pairs + bp        

    return base_pairs

def count_not_shared(list_1, list_2):

    acum = 0
 
    for l1 in list_1:
        if not l1 in list_2:
            acum += 1

    return acum
           

   
                
def bp_distance(structure_1, structure_2):


    bp1, bp2 = map(extract_base_pairs, [structure_1, structure_2])

    ns1 = count_not_shared(bp1, bp2)
    ns2 = count_not_shared(bp2, bp1)

    return ns1 + ns2


# exam
seq1_seq = "AUCAGUUCUAGCAGGAGCUGUACUCAGAGACUCGGGAAAUUUUCCCGGAAUUUUACCCGGGUUUUUACGU"
seq1_str = "..(((((((....))))))).....(((((((((((.((..(((...)))..)).)))))))))))...."

seq2_seq = "AUCGGUUCCAGCAGGAACUGUACUCGGGGGCUCGGGAAACCCUCCCGGGGUUUUACCCGGGUUUUUACGU"
seq2_str = "..(((((((....))))))).(((((((..(((((((.....)))))))......)))))))........"

ali_seq =  "AUCGGUUCCAGCAGGAACUGUACUCGGGGGCUCGGGAAACCCUCCCGGGGUUUUACCCGGGUUUUUACGU"
ali_str =  "..(((((((....))))))).(((((((..(((((((.....)))))))......)))))))........"


print "Hamming distance between structures"
print hamming_distance(seq1_str, seq2_str)

print "Base pair distance between structures"
print bp_distance(seq1_str, seq2_str)

print "Alifold results"
print "Hamming distance"
print hamming_distance(seq1_str, ali_str)
print hamming_distance(seq2_str, ali_str)


print "BP distance"
print bp_distance(seq1_str, ali_str)
print bp_distance(seq2_str, ali_str)
