#Imports
import numpy as np

#Initialising an array of amino acid keys for reference
amino_array = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X']

#Reading BLOSUM62 substitution values from a text file
# and inputting into a 2d array
f = open('BLOSUM.txt', 'r')
g = f.read().splitlines() #splitting text file into lines

for i in range(22):
    g[i] = [int(j) for j in g[i].split()] #converting each line to an array of ints

BLOSUM62 = g
f.close()

#A function for calculating the value diagonal to a given element in the matrix
def diagonal(n1, n2, substitution):
    A1 = amino_array.index(n1) #finding the index of key in amino acid array
    A2 = amino_array.index(n2)
    return substitution[A1][A2] #return substitution value from BLOSUM62

#A function for calculating the direction of the path needed for backtracking
def pointers(di, ho, ve):
    pointer = max(di, ho, ve)

    if (di == pointer):
        return 'D'
    elif(ho == pointer):
        return 'H'
    elif (ve == pointer):
        return 'V'

#Constructing the matrices and calculating the final score
def nw_blosum62(s1, s2, gap, substitution) :
    gap_penalty = gap # Penalty value

    n = len(s1) + 1
    m = len(s2) + 1

    f_matrix = np.zeros((m,n), dtype = int) # Initialises an empty alignment matrix
    p_matrix = np.zeros((m,n), dtype = str) # Initialises an empty pointer matrix for backtracking

    # Fill first row and column of matrix with gap penalty
    for i in range(m):
        f_matrix[i][0] = gap_penalty * i
        p_matrix[i][0] = 'V'

    for j in range(n):
        f_matrix[0][j] = gap_penalty * j
        p_matrix[0][j] = 'H'

    # Fill the matrix
    p_matrix[0][0] = 0
    for i in range(1,m):
        for j in range(1,n):
            di = f_matrix[i-1][j-1] + diagonal(s1[j-1], s2[i-1], substitution)
            ho = f_matrix[i-1][j] + gap_penalty
            ve = f_matrix[i][j-1] + gap_penalty
            f_matrix[i][j] = max(di, ho, ve)
            p_matrix[i][j] = pointers(di, ho, ve)

    score = f_matrix[-1][-1]

    print("Alignment Score: {0}".format(score))
    print("\n" + "Alignment Matrix:")
    print(f_matrix)
    print("\n" + "Pointer Matrix:")
    print(p_matrix)

sequence1 = 'ARGCBHTSWYGHDRPFKL'
sequence2 = 'BRGQZHTWYYGHDRMIHB'

nw_blosum62(sequence1, sequence2, -4, BLOSUM62)
