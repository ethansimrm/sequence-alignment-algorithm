"""
Algorithmic Thinking Project 4
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    This function returns a dictionary of dictionaries whose entries are 
    indexed by pairs of characters in alphabet plus ’-’. 
    The score for any entry indexed by one or more dashes is dash_score. 
    The score for the remaining diagonal entries is diag_score. 
    Finally, the score for the remaining off-diagonal entries is off_diag_score.
    """
    ans_dict = {}
    set_chars = set(alphabet)
    set_chars.add("-")
    for character_1 in set_chars:
        ans_dict[character_1] = {}
        for character_2 in set_chars:
            if character_1 == "-" or character_2 == "-":
                score = dash_score
            elif character_1 == character_2:
                score = diag_score
            else: 
                score = off_diag_score
            ans_dict[character_1][character_2] = score
    return ans_dict

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Takes as input two sequences seq_x and seq_y whose elements 
    share a common alphabet with scoring_matrix. The function computes 
    and returns the alignment matrix for seq_x and seq_y 
    as described in the Homework. If global_flag is True, 
    each entry of the alignment matrix is computed using the method 
    described in Question 8 of the Homework. If global_flag is False, 
    each entry is computed using the method described in Question 12 
    of the Homework.
    """
    row_num = len(seq_x) + 1
    col_num = len(seq_y) + 1
    #seq_x will index the rows, and seq_y will index the columns.
    #We need to use list comprehension (and not simple multiplication) to avoid reference sharing.
    ans_matrix = [[0 for dummy_col in range(col_num)] for dummy_row in range(row_num)]
    for row_idx in range(1, row_num):
        ans_matrix[row_idx][0] = ans_matrix[row_idx - 1][0] + scoring_matrix[seq_x[row_idx - 1]]["-"]
        #If global_flag is False, then we do local alignment - we want to update this change at every step, because
        #subsequent steps will depend on previous steps' values.
        if not global_flag:
            if ans_matrix[row_idx][0] < 0:
                ans_matrix[row_idx][0] = 0
    for col_idx in range(1, col_num):
        ans_matrix[0][col_idx] = ans_matrix[0][col_idx - 1] + scoring_matrix["-"][seq_y[col_idx - 1]]
        if not global_flag:
            if ans_matrix[0][col_idx] < 0:
                ans_matrix[0][col_idx] = 0
    for row_idx in range(1, row_num):
        for col_idx in range(1, col_num):
            #We take the maximum of three cases, namely the diagonal, left, and top boxes
            #of the matrix.
            case_diag = ans_matrix[row_idx - 1][col_idx - 1] + scoring_matrix[seq_x[row_idx - 1]][seq_y[col_idx - 1]]
            case_left = ans_matrix[row_idx - 1][col_idx] + scoring_matrix[seq_x[row_idx - 1]]["-"]
            case_top = ans_matrix[row_idx][col_idx - 1] + scoring_matrix["-"][seq_y[col_idx - 1]]
            ans_matrix[row_idx][col_idx] = max(case_diag, case_left, case_top)
            if not global_flag:
                if ans_matrix[row_idx][col_idx] < 0:
                    ans_matrix[row_idx][col_idx] = 0    
    return ans_matrix

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet 
    with scoring_matrix. This function computes a global alignment of seq_x and seq_y 
    using the global alignment matrix alignment_matrix.The function returns a tuple 
    of the form (score, align_x, align_y) where score is the score of the global alignment 
    align_x and align_y. Note that align_x and align_y should have the same length and 
    may include the padding character ’-’.
    """
    pointer_x = len(seq_x)
    pointer_y = len(seq_y)
    align_x = ""
    align_y = ""
    while pointer_x != 0 and pointer_y != 0:
        #Here, we traceback from the maximum value (lower right of matrix) to determine the previous maximum values which 
        #gave rise to it - this allows us to build back the alignment.
        #Why is the lower right value the maximum? Because global alignments must be of the same length - we trace paths
        #"towards" the lower right as we add characters (adding char, char means diag down/right) or dashes to either sequence,
        #until we hit that length. This is the alignment matrix. Of course, we penalise additions per the scoring matrix.
        if alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x - 1][pointer_y - 1] + scoring_matrix[seq_x[pointer_x - 1]][seq_y[pointer_y - 1]]:
            align_x = seq_x[pointer_x - 1] + align_x
            align_y = seq_y[pointer_y - 1] + align_y
            pointer_x -= 1
            pointer_y -= 1
        elif alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x - 1][pointer_y] + scoring_matrix[seq_x[pointer_x - 1]]["-"]:
            align_x = seq_x[pointer_x - 1] + align_x
            align_y = "-" + align_y
            pointer_x -= 1
        elif alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x][pointer_y - 1] + scoring_matrix["-"][seq_y[pointer_y - 1]]:
            align_x = "-" + align_x
            align_y = seq_y[pointer_y - 1] + align_y
            pointer_y -= 1
    while pointer_x != 0:
        align_x = seq_x[pointer_x - 1] + align_x
        align_y = "-" + align_y
        pointer_x -= 1
    while pointer_y != 0:
        align_x = "-" + align_x
        align_y = seq_y[pointer_y - 1] + align_y
        pointer_y -= 1
    return (alignment_matrix[len(seq_x)][len(seq_y)], align_x, align_y)

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Like compute_global_alignment, but alignment_matrix is a local alignment matrix.
    The function returns a tuple of the form (score, align_x, align_y) where score is 
    the score of the optimal local alignment align_x and align_y. Note that align_x 
    and align_y should have the same length and may include the padding character ’-’.
    """
    #First, scan the matrix for the maximum value
    max_value = 0
    max_row = 0
    max_col = 0
    for row_num in range(len(seq_x) + 1):
        for col_num in range(len(seq_y) + 1):
            if alignment_matrix[row_num][col_num] > max_value:
                max_value = alignment_matrix[row_num][col_num]
                max_row = row_num
                max_col = col_num
    pointer_x = max_row
    pointer_y = max_col
    align_x = ""
    align_y = ""
    #Now, we apply the technique from compute_global_alignment.
    #Note that we do not pad with dashes before and after, because this
    #would unnecessarily decrease the score. We just want a local alignment,
    #which only considers the shortest relevant section of the longer sequence
    #For example, if we align AA and TAAT we just want AA and AA to be returned.
    while pointer_x != 0 and pointer_y != 0:
        if alignment_matrix[pointer_x][pointer_y] == 0:
            break
        elif alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x - 1][pointer_y - 1] + scoring_matrix[seq_x[pointer_x - 1]][seq_y[pointer_y - 1]]:   
            align_x = seq_x[pointer_x - 1] + align_x
            align_y = seq_y[pointer_y - 1] + align_y
            pointer_x -= 1
            pointer_y -= 1
        elif alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x - 1][pointer_y] + scoring_matrix[seq_x[pointer_x - 1]]["-"]:
            align_x = seq_x[pointer_x - 1] + align_x
            align_y = "-" + align_y
            pointer_x -= 1
        elif alignment_matrix[pointer_x][pointer_y] == alignment_matrix[pointer_x][pointer_y - 1] + scoring_matrix["-"][seq_y[pointer_y - 1]]:
            align_x = "-" + align_x
            align_y = seq_y[pointer_y - 1] + align_y
            pointer_y -= 1
    return(max_value, align_x, align_y)  
                   
