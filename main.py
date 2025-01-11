#blastp
import pandas as pd
from utils import generate_neighborhood, calculate_alignment_score
from fasta import read_fasta
from database import get_database
from local_alignment import local_alignment
from config import K, T, SCORING_MATRIX, SCORING_MATRIX_NAME, HSP_T, LOCAL_ALIGNMENT_MARGIN, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.style.use("ggplot")
# from rich import print



def merge_overlap_seqs(seed_hits: list) -> list:
    merged_seed_hits = []
    while True:
        tmp_len = len(seed_hits)
        while seed_hits:
            seed_hit = seed_hits.pop(0)
            i_range, j_range = range(seed_hit[0][0], seed_hit[0][1]+1), range(seed_hit[1][0], seed_hit[1][1]+1)
            for idx, next_seed_hit in enumerate(seed_hits):
                if next_seed_hit[0][0] in i_range and next_seed_hit[1][0] in j_range and \
                        next_seed_hit[0][1] - seed_hit[0][1] == next_seed_hit[1][1] - seed_hit[1][1]:
                    seed_hit = (seed_hit[0][0], next_seed_hit[0][1]), (seed_hit[1][0], next_seed_hit[1][1])
                    i_range, j_range = range(seed_hit[0][0], seed_hit[0][1]+1), range(seed_hit[1][0], seed_hit[1][1]+1)
                    seed_hits.pop(idx)
            merged_seed_hits.append(seed_hit)

        if len(merged_seed_hits) == tmp_len:
            break
        else:
            seed_hits = merged_seed_hits
            merged_seed_hits = []

    return merged_seed_hits

def find_seed_hits(seq: str, database: dict, k: int, threshold: int, scoring_matrix: pd.DataFrame) -> dict:
    neighborhood1 = generate_neighborhood(seq, k, threshold, scoring_matrix)
    seed_hits = {}
    for name, (seq2, neighborhood2) in database.items():
        current_seed_hits = []
        for i, neighbors1 in enumerate(neighborhood1):
            for j, neighbors2 in enumerate(neighborhood2):
                if neighbors1 & neighbors2:
                    seed_hit = ((i, i+k), (j, j+k))
                    current_seed_hits.append(seed_hit)

        seed_hits[name] = (seq2, merge_overlap_seqs(current_seed_hits))

    return seed_hits

def extend_hsps(seq: str, seed_hits: dict, threshold: int, scoring_matrix: pd.DataFrame) -> tuple:
    hsps = {}
    for name, (seq2, seed_hit_list) in seed_hits.items():
        current_hsps = []
        for seed_hit in seed_hit_list:
            (left_i, right_i), (left_j, right_j) = seed_hit[0], seed_hit[1]
            max_score = calculate_alignment_score(seq[left_i: right_i], seq2[left_j: right_j], scoring_matrix)

            current_score = max_score
            x, y = left_i, left_j
            while x > 0 and y > 0:
                current_score += scoring_matrix.loc[seq[x-1], seq2[y-1]]
                if current_score > max_score:
                    max_score = current_score
                    left_i, left_j = x-1, y-1
                if current_score < max_score - threshold:
                    break
                x, y = x-1, y-1

            current_score = max_score
            x, y = right_i, right_j
            while x < len(seq)-1 and y < len(seq2)-1:
                current_score += scoring_matrix.loc[seq[x], seq2[y]]
                if current_score > max_score:
                    max_score = current_score
                    right_i, right_j = x+1, y+1
                if current_score < max_score - threshold:
                    break
                x, y = x+1, y+1

            if max_score >= threshold:
                current_hsps.append(((left_i, right_i), (left_j, right_j)))

        hsps[name] = (seq2, merge_overlap_seqs(current_hsps))

    return hsps

def local_alignment_hsps(seq: str, hsps: dict, margin: int, scoring_matrix: pd.DataFrame, open_gap_penalty: int, extend_gap_penalty: int,):
    alignments = {}
    for name, (seq2, hsp_list) in hsps.items():
        current_alignments = ("", "", (), 0)
        max_score = 0
        for (i_left, i_right), (j_left, j_right) in hsp_list:
            sub_seq = seq[max(i_left-margin, 0): min(i_right+margin, len(seq))]
            sub_seq2 = seq2[max(j_left-margin, 0): min(j_right+margin, len(seq2))]

            align_1, align_2, pos, score = local_alignment(sub_seq, sub_seq2, open_gap_penalty, extend_gap_penalty, scoring_matrix)
            pos_ = (
                (
                    pos[0][0] + max(i_left-margin, 0),
                    pos[0][0] + max(i_left-margin, 0) + len(list(filter(lambda x: x != "-", list(align_1))))
                ),
                (
                    pos[1][0] + max(j_left-margin, 0),
                    pos[1][0] + max(j_left-margin, 0) + len(list(filter(lambda x: x != "-", list(align_2))))
                )
            )

            if score > max_score:
                max_score = score
                current_alignments = (align_1, align_2, pos_, score)

        alignments[name] = (seq2, current_alignments)

    return alignments

def plot(seq: str, seq_name: str, seed_hits: dict, hsps: dict, alignments: dict):
    gap_pp = {
        "label": "Gap",
        "c": "red",
        "lw": 0,
        "alpha": 1,
        "ls": "None",
        "marker": "o",
        "zorder": 10
    }
    match_pp = {
        "label": "Match",
        "c": "green",
        "lw": 0,
        "alpha": 1,
        "ls": "None",
        "marker": "o",
        "zorder": 10
    }
    mismatch_pp = {
        "label": "Mismatch",
        "c": "orange",
        "lw": 0,
        "alpha": 1,
        "ls": "None",
        "marker": "o",
        "zorder": 10
    }
    seed_hit_pp = {
        "label": "Seed Hit",
        "c": "pink",
        "lw": 3,
        "alpha": 1,
        "ls": "-",
        "marker": None,
        "zorder": 8
    }
    hsp_pp = {
        "label": "HSP",
        "c": "purple",
        "lw": 1,
        "alpha": 1,
        "ls": "-",
        "marker": None,
        "zorder": 9
    }

    for seq2_name, (seq2, seed_hit_list), (_, hsp_list), (_, (align1, align2, align_pos, ـ)) in zip(seed_hits.keys(), seed_hits.values(), hsps.values(), alignments.values()):
        fig, ax = plt.subplots(figsize=(15, 5))

        for i, j in seed_hit_list:
            ax.plot([j[0], j[1]], [i[0], i[1]], **seed_hit_pp)

        for i, j in hsp_list:
            ax.plot([j[0], j[1]], [i[0], i[1]], **hsp_pp)

        if align_pos:
            i_left, j_left = align_pos[0][0]-1, align_pos[1][0]-1
            for char1, char2 in zip(align1, align2):
                if char1 == "-":
                    j_left += 1
                    ax.scatter([j_left], [i_left], **gap_pp)
                elif char2 == "-":
                    i_left += 1
                    ax.scatter([j_left], [i_left], **gap_pp)
                elif char1 == char2:
                    i_left += 1
                    j_left += 1
                    ax.scatter([j_left], [i_left], **match_pp)
                else:
                    i_left += 1
                    j_left += 1
                    ax.scatter([j_left], [i_left], **mismatch_pp)

        ax.set_yticks(range(len(seq)))
        _ = len(str(len(seq)+1))
        ax.set_yticklabels(map(lambda x: f"{x[1]} {str(x[0]+1).zfill(_)}", enumerate(list(seq))), fontsize=9)
        ax.set_ylim(-0.5, len(seq) - 0.5)
        ax.set_ylabel(f"Q: {seq_name}", fontsize=9)

        ax.set_xticks(range(len(seq2)))
        _ = len(str(len(seq2)+1))
        ax.set_xticklabels(map(lambda x: f"{x[1]} {str(x[0]+1).zfill(_)}", enumerate(list(seq2))), rotation=45, fontsize=9)
        ax.set_xlim(-0.5, len(seq2) - 0.5)
        ax.set_xlabel(seq2_name, fontsize=9)

        handles = [Line2D([], [],  **pp) for pp in [seed_hit_pp, hsp_pp, gap_pp, match_pp, mismatch_pp]]
        labels = [pp["label"] for pp in [seed_hit_pp, hsp_pp, gap_pp, match_pp, mismatch_pp]]
        ax.legend(handles, labels)

        # ax.grid(True, "major", color="gray", alpha=0.5, zorder=1)

        fig.tight_layout()
        fig.savefig(f"Q_{seq_name.split()[0]}___D_{seq2_name.split()[0]}.png", dpi=600)
        plt.close(fig)

def result(seq: str, seq_name: str, alignments: dict) -> str:

    alignments_ = sorted(alignments.items(), key=lambda x: x[1][1][3], reverse=True)

    r = f"Query:\t{seq_name}\n\t{seq}\n\tLength: {len(seq)}"
    r += "\n" + "="*20 + "\n\n"

    for i, (seq2_name, (seq2, (a1, a2, pos, score))) in enumerate(alignments_):
        seq2 = "\n".join(["\t"+seq2[i:i+70] for i in range(0, len(seq2), 70)])
        r += f"\n{i+1} DB:\t{seq2_name}\n{seq2}\n\tLength: {len(seq2)}"
        r += "\n"
        r += f"\nMax Score: {score}"
        r += f"\nQ: {a1}"
        r += f"\nD: {a2}"
        r += "\n\n" + "-"*20 + ""

    return r

SEQ_NAME, SEQUENCE = list(read_fasta("query.fasta").items())[0]
SEQ_NAME = SEQ_NAME[1:]
DATABASE = get_database("database.fasta", K, T, SCORING_MATRIX, SCORING_MATRIX_NAME)

seed_hits = find_seed_hits(SEQUENCE, DATABASE, K, T, SCORING_MATRIX)
hsps = extend_hsps(SEQUENCE, seed_hits, HSP_T, SCORING_MATRIX)
alignments = local_alignment_hsps(SEQUENCE, hsps, LOCAL_ALIGNMENT_MARGIN, SCORING_MATRIX, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY)
r = result(SEQUENCE, SEQ_NAME, alignments)
plot(SEQUENCE, SEQ_NAME, seed_hits, hsps, alignments)

with open("result.txt", "w") as f:
    f.writelines(r)
