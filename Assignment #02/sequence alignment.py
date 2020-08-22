import numpy as np


"""function for matching and mismatching scores"""


def matching_score(seq1, seq2, match, mismatch):
    if seq1 == seq2:
        return match
    else:
        return mismatch


""" main function """


def global_alignment(seq1, seq2, gap, match, mismatch):
    sum = 0
    alignment = []
    scores = []
    row = len(seq1)
    col = len(seq2)
    dp_matrix = np.zeros((row + 1, col + 1))

    """ filling the dynamic programming matrix"""

    for i in range(row + 1):
        dp_matrix[i, 0] = i * gap
    for j in range(col + 1):
        dp_matrix[0, j] = j * gap
    for i in range(row + 1):
        for j in range(col + 1):
            dp_matrix[i, j] = max(
                dp_matrix[i - 1, j - 1] + matching_score(seq1[i - 1], seq2[j - 1], match, mismatch),
                dp_matrix[i - 1, j] + gap, dp_matrix[i, j - 1] + gap)

    """" trace back for getting sequence layout"""

    while row > 0 and col > 0:

        if (dp_matrix[row - 1, col - 1] + matching_score(seq1[row - 1], seq2[col - 1], match, mismatch)) \
                == dp_matrix[row, col]:  # diagonal trace back

            # add matching or mismatching scores to scores list
            scores.append(matching_score(seq1[row - 1], seq2[col - 1], match, mismatch))
            alignment.append(seq1[row - 1])  # add base to alignment list
            row -= 1
            col -= 1

        elif dp_matrix[row - 1, col] + gap == dp_matrix[row, col]:  # up trace back
            scores.append(gap)  # add gap to scores list
            alignment.append('-')   # add '-' to alignment list
            row -= 1

        elif dp_matrix[row, col - 1] + gap == dp_matrix[row, col]:  # left trace back
            scores.append(gap)  # add gap to scores list
            alignment.append('-')
            col -= 1

    alignment.reverse()

    """ summation of all elements in the array to get the score of aligned sequence """

    for i in range(len(scores)):
        sum = sum + scores[i]
        i += 1

    return alignment, sum


seq1 = ['C', 'T', 'A', 'T', 'T', 'G', 'A', 'A', 'C', 'A', 'T']
seq2 = ['C', 'T', 'A', 'T', 'T', 'G', 'A', 'C', 'G', 'T', 'A', 'A', 'C', 'A', 'T']
match = 4
mismatch = -1
gap = -3
layout, sum = global_alignment(seq1, seq2, gap, match, mismatch)
print('Aligned sequence layout:', layout)
print('The score =', sum)
