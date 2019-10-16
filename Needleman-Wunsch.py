import numpy as np

def diagonal(n1, n2, pt):
    if (n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']

def pointers(di, ho, ve):
    pointer = max(di, ho, ve)

    if (di == pointer):
        return 'D'
    elif(ho == pointer):
        return 'H'
    elif (ve == pointer):
        return 'V'

def nw_simple(s1, s2, match, mismatch, gap) :
    penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap} # Dictionary for penalty values

    n = len(s1) + 1
    m = len(s2) + 1

    f_matrix = np.zeros((m,n), dtype = int) # Initialises an empty alignment matrix
    p_matrix = np.zeros((m,n), dtype = str) # Initialises an empty pointer matrix for backtracking

    # Fill first row and column of matrix with gap penalty
    for i in range(m):
        f_matrix[i][0] = penalty['GAP'] * i
        p_matrix[i][0] = 'V'

    for j in range(n):
        f_matrix[0][j] = penalty['GAP'] * j
        p_matrix[0][j] = 'H'

    # Fill the matrix
    p_matrix[0][0] = 0
    for i in range(1,m):
        for j in range(1,n):
            di = f_matrix[i-1][j-1] + diagonal(s1[j-1], s2[i-1], penalty)
            ho = f_matrix[i-1][j] + penalty['GAP']
            ve = f_matrix[i][j-1] + penalty['GAP']
            f_matrix[i][j] = max(di, ho, ve)
            p_matrix[i][j] = pointers(di, ho, ve)

    score = f_matrix[-1][-1]

    print("Alignment Score: {0}".format(score))
    print("\n" + "Alignment Matrix:")
    print(f_matrix)
    print("\n" + "Pointer Matrix:")
    print(p_matrix)

    sequence1 = 'TGCCA'
sequence2 = 'TCCA'

nw_simple(sequence1, sequence2, 1, -1, -2)
