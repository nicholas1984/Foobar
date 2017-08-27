"""
This is my solution to Google Foobar problem 3.2.  
"""

from fractions import Fraction

def split_states(transition_mat):
    """
    reads in a state transition matrix and returns two lists; one which 
    contains the terminal states the other with the transient states.
    
    example:
        transition = [[0,1,0,0,0,1],  [4,0,0,3,2,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0] ]
        split_states(transition) = ([2, 3, 4, 5], [0, 1])
    """
    terminal = []
    transient = []
    for index, array in enumerate(transition_mat):
        if sum(array) == 0:
            terminal.append(index)
        else:
            transient.append(index)
    return terminal, transient

def get_row_probabilities(transition_mat, row_index, states, matrix='R'):
    """
    transition_mat: matrix with entry (i,j) representing the number of oberved 
        transitions from state i to state j. 
    row_index: the specific row in transition_mat we want to convert to probabilities.
    states: list of either 1) transient state indices, or 2) terminal state indices.
    matirx: argument determining if we want rows of Ninv or R.

    example:
        transition_mat = [[0,1,0,0,0,1],  [4,0,0,3,2,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0] ]
        transient_row = 0
        term_states = [2, 3, 4, 5]
        get_row_probabilities(transition_mat, transient_row, term_states) 
            = [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 2)]
    """
    denominator = sum(transition_mat[row_index])
    if matrix == 'R':
        return [Fraction(transition_mat[row_index][col], denominator) for col in states]
    else:
        return [1-Fraction(transition_mat[row_index][col], denominator) if row_index == col else Fraction(-transition_mat[row_index][col], denominator) for col in states]

def get_R_and_Ninv(transition_mat, trans_states, term_states):
    """
    transition_mat: matrix with entry i,j representing the number of oberved 
        transitions from state i to state j.
    trans_states: list of indices for transient states.
    term_states: list of indices for absorbing states.
    
    Ninv: inverse of the fundamental matrix.
    R is a matrix with entriy (i,j) representing the probability of transitioning 
        from transitional state i to absorbing state j.
        
    example:
        transition_mat = [[0,1,0,0,0,1],  [4,0,0,3,2,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0] ]
        trans_states = [0, 1]
        term_states = [2, 3, 4, 5]
        get_R_and_Ninv(transition_mat, trans_states, term_states)
            R = [[Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 2)],[Fraction(0, 1), Fraction(1, 3), Fraction(2, 9), Fraction(0, 1)]]
            Ninv = [[Fraction(1, 1), Fraction(-1, 2)], [Fraction(-4, 9), Fraction(1, 1)]]
    """
    R = []
    Ninv = []
    for transient_row in trans_states:
        R.append(get_row_probabilities(transition_mat, transient_row, term_states))
        Ninv.append(get_row_probabilities(transition_mat, transient_row, trans_states, 'Ninv'))
    return R, Ninv

def elem_row_operation(A, row_to_scale, row_to_replace, Ainv=False):
    """
    A: matrix we are performing a row operation on.
    row_to_scale: index of row in A we will scale.
    row_to_replace: index of row in A we will replace.
    Ainv(optional): perform same elemantary row operation to compute A inverse.
    Performs an elementary row operation on matrix A and Ainv.
    
    example:
        A = [[Fraction(1,1), Fraction(1,2), Fraction(0,1)], 
             [Fraction(1,2), Fraction(1,5), Fraction(1,1)], 
             [Fraction(1,3), Fraction(1,1), Fraction(1,1)]
            ]
        row_to_scale = 0
        row_to_replace = 1
        Ainv = [ [Fraction(1,1), Fraction(0,1), Fraction(0,1) ], 
                 [Fraction(0,1), Fraction(1,1), Fraction(0,1)],
                 [Fraction(0,1), Fraction(0,1), Fraction(1,1)]
               ]
        
        elem_row_operation(A, row_to_scale, row_to_replace, Ainv)
        A = [[Fraction(1,1), Fraction(1,2), Fraction(0,1)], 
             [Fraction(0,1), Fraction(-1,20), Fraction(1,1)], 
             [Fraction(1,3), Fraction(1,1), Fraction(1,1)]
            ]
        Ainv = [ [Fraction(1,1), Fraction(0,1), Fraction(0,1) ], 
                 [Fraction(-1,2), Fraction(1,1), Fraction(0,1)],
                 [Fraction(0,1), Fraction(0,1), Fraction(1,1)]
        
    """
    if row_to_scale == row_to_replace:
        scalar = 1/A[row_to_scale][row_to_scale]
        A[row_to_replace] = [scalar*A[row_to_scale][col] for col in range(len(A[row_to_replace])) ]
        if Ainv:
            Ainv[row_to_replace] = [scalar*Ainv[row_to_scale][col] for col in range(len(Ainv[row_to_replace])) ]
    else: 
        scalar = -A[row_to_replace][row_to_scale]/A[row_to_scale][row_to_scale]
        A[row_to_replace] = [scalar*A[row_to_scale][col] + A[row_to_replace][col] for col in range(len(A[row_to_replace])) ]
        if Ainv:
            Ainv[row_to_replace] = [scalar*Ainv[row_to_scale][col] + Ainv[row_to_replace][col] for col in range(len(Ainv[row_to_replace])) ]

def swap_rows(A, row0, row1):
    """
    A: a matrix with more than 1 row
    row0: index of a row to interchange
    row1: index of another row to interchange
    
    example:
        A = [ [1,0], [0,1]]
        row0 = 0
        row1 = 1
        swap_rows(A, row0, row1) = [[0, 1], [1, 0]]
    """
    temp = A[row0]
    A[row0] = A[row1]
    A[row1] = temp
    
def identify_pivot(M, index):
    """
    M: a matrix such that M[index][index] = 0
    index: an index for M

    example:
        M = [ [0,1], [1,0] ]
        index = 0
        identify_pivot(M, index) = 1
    """
    for i in range(index, len(M)):
        if M[i][index] != 0:
            return i

def get_triangular(M, R, Minv=False):
    """
    M matrix to put in triangular form
    Minv matrix which we transform to the inverse of M
    R ordered list either range(len(M)) or the reverse.
    
    If R == range(len(M)) transform M to upper triangular by performing elementary row operations.
    If R == reverse(list(range(len(M)))) transform M to lower triangular by performing elementary row operations.
    
    example:
        M = [ [1, 1],
              [1, 2]
            ]
        Minv = [ [1, 0],
                 [0, 1]
               ]
        R = [0, 1]
        get_triangular(M, Minv, R)
        M = [ [1, 1],
              [0, 1]
            ]
        Minv = [ [1,  0],
                 [-1, 1]
               ]
    """
    for keep_index in R[0:len(M)-1]:
            for delete_index in R[1:len(M)]:
                if M[delete_index][keep_index] == 0:
                    continue
                else:
                    if M[keep_index][keep_index] == 0:
                        swap_index = identify_pivot(M, keep_index)
                        swap_rows(M, swap_index, keep_index)
                        if Minv:
                            swap_rows(Minv, swap_index, keep_index)
                    elem_row_operation(M, keep_index, delete_index, Minv)
            if M[keep_index][keep_index] == 1:
                continue
            else:
                elem_row_operation(M, keep_index, keep_index, Minv)

def invert_matrix(M):
    """
    Invert the matrix M using Gaussian elimination
    Note: we do not use pivoting because diagonal entries will be nonzero.
    
    example:
        M = [
             [Fraction(1,1), Fraction(1,2), Fraction(0,1)], 
             [Fraction(1,2), Fraction(1,5), Fraction(1,1)], 
             [Fraction(1,3), Fraction(1,1), Fraction(1,1)]
            ]
        
        invert_matrix(M) = [
                            [Fraction(48, 53), Fraction(30, 53), Fraction(-30, 53)],
                            [Fraction(10, 53), Fraction(-60, 53), Fraction(60, 53)],
                            [Fraction(-26, 53), Fraction(50, 53), Fraction(3, 53)]
                           ]
    """
    if len(M) == 2:
        det_inv = Fraction(1, M[0][0]*M[1][1] - M[0][1]*M[1][0])
        Minv = [ [det_inv*M[1][1], -det_inv*M[0][1] ], [ -det_inv*M[1][0], det_inv*M[0][0] ] ]
    else:
        R = range(len(M))
        Minv = [ [1 if i == j else 0 for j in R] for i in R]
        get_triangular(M, R, Minv)
        get_triangular(M, list(reversed(R)), Minv)
    return Minv
    
def mat_mul_row(A, B, row_index):
    """
    A and B are matrices to be multiplied.
    retrun (A*B)[row]
    """

    return [ sum(A[row_index][i]*B[i][j] for i in range(len(A[0]))) for j in range(len(B[0])) ]

            
        


def mat_mul(A, B):
    """
    A and B are square matrices
    Computes AxB
    
    example:
        A = [
             [Fraction(1,1), Fraction(1,2), Fraction(0,1)], 
             [Fraction(1,2), Fraction(1,5), Fraction(1,1)], 
             [Fraction(1,3), Fraction(1,1), Fraction(1,1)]
            ]
        B = [
             [Fraction(48, 53), Fraction(30, 53), Fraction(-30, 53)],
             [Fraction(10, 53), Fraction(-60, 53), Fraction(60, 53)],
             [Fraction(-26, 53), Fraction(50, 53), Fraction(3, 53)]
            ]
        mat_mul(A,B) = [
                         [Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)],
                         [Fraction(0, 1), Fraction(1, 1), Fraction(0, 1)],
                         [Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)]
                       ]
    """
    return [ mat_mul_row(A, B, index) for  index, row in enumerate(A) ]

def get_lcd_list(denominators):
    """
    calculate the lowest common denominator of a list of numbers
    """
    lcd = 1
    for val in denominators:
        lcd = get_lcd(lcd, val)
    return lcd

def get_lcd(a, b):
    """
    calculate lowest common denominator of a and b
    """
    return int(a*b/gcd(a,b))
    
def gcd(a, b):
    """
    calculate greates common divisor of a and b using Euclid's algorithm
    """
    while b != 0:
       t = b 
       b = a % b
       a = t
    return a    
    
def answer(transition_mat):
    """
    transition_mat is a square matrix with entry (i,j) representing the number 
    of times state i has transitioned to state j.  
    
    We return an array, 'ans', of length n where the probability of starting in state 0
    and ending up in the absorbing state of smallest index is ans[0]/ans[n-1], etc.
    
    To accomplish this we will form the components for the canonical form of a 
    Markov state transition matrix.  
    
    For further details see https://en.wikipedia.org/wiki/Absorbing_Markov_chain.
    
    example:
        transition = [[0,1,0,0,0,1],  [4,0,0,3,2,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0], ]
        solution(transition) = [0, 3, 2, 9, 14]
    """
    if len(transition_mat) == 1:
        return [1,1]
    term_states, trans_states = split_states(transition_mat)
    #print(term_states, trans_states)
    R, Ninv = get_R_and_Ninv(transition_mat, trans_states, term_states)   
    N = invert_matrix(Ninv)
    #print(N)
    #print(R)
    probabilities = mat_mul_row(N, R, 0)
    denom = max([val.denominator for val in probabilities])
    ans = [ int(val.numerator*(denom/val.denominator)) for val in probabilities ]
    ans.append(denom)
    return ans
    

#m = [[1,1], [0,0]]
#x = mat_mul(m,m)
#print(x)
#m = [[0,1,0,0,0,1],  [4,0,0,3,2,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0],  [0,0,0,0,0,0], ]
m = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
print(answer(m))



#eof