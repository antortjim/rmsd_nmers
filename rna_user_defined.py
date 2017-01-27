import numpy as np


def hamming_distance(s1, s2):
    '''Returns hamming distance between s1 and s2'''
    acum = 0 # number of mismatches
    for i in range(len(s1)): # for every position in s1
        if structure_1[i] != s2[i]: # if position i in s1
                                             # is not equal to s2[i]
            acum += 1                        # sum 1 to acum



    return acum

def check_bifurcation(structure):
    '''Checks for bifurcations in the dot bracket notation.
    One bifurcation is represented as a tuple of 2 indices
    First one states where the bifurcation starts (the first "(")
    Las ones states where the bifurcation ends (last ")")
    If there are no bifurcations, the list of tuples has length 1
    and the indices are those of the first ( and the last ) in the
    whole structure'''

    # Initialize an acum that will store how many opening parentheresis
    # remain to be closed
    # Everytime acum reaches 0 after havin been positive counts
    # as one bifurcation (the end of one)
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
    '''Extracts the start and end indices of all giving base pairs
    in a dot bracket structure'''

    bifurcations = check_bifurcation(structure) # Checks for bifurcation
                                                # The code below is repeated for every
                                                # bifurc because they are independent folding units
                                                # and thus must be processed independtly
    base_pairs = []                             # Initialize list of tuples (start_id, end_id)
    for b in bifurcations:
        b_start = b[0] 
        b_end = b[1] + 1
        substring = structure[b_start:b_end]
        reverse = substring[::-1]
        acum = 0
       
     
        opening = {}
        closing = {}

        for idx, c in enumerate(substring):
            idx = idx + b_start # start to count from start of complete structure, not the substring
            if c == "(":
                acum += 1

                try:
                    opening[str(acum)].append(idx)

                except KeyError:
                    opening[str(acum)] = [idx]

            if c == ")":
                acum -= 1

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
    '''Counts which tuples are available in list_1 and not in list_2
    Input lists are extract_base_pairs output'''
    acum = 0
 
    for l1 in list_1:
        if not l1 in list_2:
            acum += 1

    return acum
           

                   
def bp_distance(structure_1, structure_2):
    '''Returns the bp_distance, computed as the return of count_not_shared
    in both directions'''

    bp1, bp2 = map(extract_base_pairs, [structure_1, structure_2])

    ns1 = count_not_shared(bp1, bp2)
    ns2 = count_not_shared(bp2, bp1)

    return ns1 + ns2
