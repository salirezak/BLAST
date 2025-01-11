import numpy as np
import pandas as pd

def local_alignment(seq1: str, seq2: str, open_gap_penalty: int, extend_gap_penalty: int, match_scoring_matrix: pd.DataFrame) -> tuple:
    score_matrix = np.zeros((len(seq2)+1, len(seq1)+1))
    dir_matrix = np.full((len(seq2)+1, len(seq1)+1), "", dtype=object)
    dir_matrix[0, 1:] = "R"
    dir_matrix[1:, 0] = "D"

    for i in range(1, len(seq2)+1):
        for j in range(1, len(seq1)+1):
            match_score = match_scoring_matrix[seq2[i-1]][seq1[j-1]]
            gap_penalty_d = extend_gap_penalty if "D" in dir_matrix[i-1, j] else open_gap_penalty
            gap_penalty_r = extend_gap_penalty if "R" in dir_matrix[i, j-1] else open_gap_penalty

            possible_scores = np.array([
                score_matrix[i-1, j-1] + match_score,
                score_matrix[i-1, j  ] - gap_penalty_d,
                score_matrix[i  , j-1] - gap_penalty_r,
                0
            ])

            score_matrix[i, j] = possible_scores.max()
            dir_matrix[i, j] = ["M", "D", "R", ""][possible_scores.argmax()]

    max_score = np.max(score_matrix)
    i_right, j_right = np.argwhere(score_matrix == max_score)[-1]


    align_1, align_2 = "", ""
    i, j = i_right, j_right
    while score_matrix[i, j] != 0:
        directions = dir_matrix[i, j]
        if "M" in directions:
            i -= 1
            j -= 1
            align_1 += seq1[j]
            align_2 += seq2[i]
        elif "D" in directions:
            i -= 1
            align_1 += "-"
            align_2 += seq2[i]
        elif "R" in directions:
            j -= 1
            align_1 += seq1[j]
            align_2 += "-"
    i_left, j_left = i, j

    align_1 = align_1[::-1]
    align_2 = align_2[::-1]

    return align_1, align_2, ((j_left, j_right), (i_left, i_right)), max_score