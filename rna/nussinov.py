import pandas as pd
import numpy as np


def matrix_initalizer(seq_length):
    """
    Generate nussinov matrix required to run algorithm on seq
    """
    nussinov_table = np.zeros([seq_length, seq_length], dtype = int)
    return nussinov_table


def check_binding(letter_1, letter_2):
    """
    Check whether letter 2 is allowed to bind letter 1
    """
    allowed_binding = {'A' : 'U', 'C' : 'G', 'U' : ['A', 'G'] , 'G' : ['C', 'U']}

    result = letter_2 in allowed_binding[letter_1]
    
    return result


def nussinov_algorithm(seq, m, bifurcations):
     
    M = matrix_initalizer(len(seq))
    l = len(seq)

    D = matrix_initalizer(len(seq))

    #m = size of minimal loop in window, trivially 1
    #d = id of the diagonal we are filling.
    #its index is 0 based, so we start to fill the diagonal with id 1+m
    #that is, for m = 1, the diagonal filled first is #2.
    for d in range(m + 1, l):
        #i should go from 0 (the first row, up to l-d
        #the bigger d, the less rows we have to explore
        for i in range(0, l - d):
            j = i + d

            letter_1 = seq[i]
            letter_2 = seq[j]
 
            option_0 = M[i, j - 1]
            option_1 = M[i + 1, j]
            
            if check_binding(letter_1, letter_2):
                option_2 = M[i + 1, j - 1] + 1
            else:
                option_2 = 0

            options = [option_0, option_1, option_2]

            if bifurcations:
                best = 0
                best_k = i + 1
                for k in range(i + 1, j):
                    score = M[i, k] + M[k + 1, j]
                    best = score if score > best else best
                    best_k = k if score > best else best_k

                option_4 = best
                options.append(option_4)
          
            decision = options.index(max(options))
            print decision
            maximum = options[decision]
 
            M[i, j] = maximum
            D[i, j] = decision
            #if decision == 3:
            #    D[i, j] = -best_k

    return (M, D)

def traceback(D):
    seq_len = D.shape[0]
    result = "." * seq_len

    D[seq_len, seq_len]
