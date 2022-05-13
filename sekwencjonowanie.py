import numpy as np
from typing import *
from random import choice


def get_offsets_matrix(names: list) -> np.ndarray:
    """
    Calculate offset for each pair of oligonucleotides.
    Offset is a number of letters (nucleotides) we have to move in the first oligonucleotide to achieve the beginning of the second one.
    Eg. CGTT, GTTG -> 1, because we have to move by only one letter – C.
    Later it is computed only for the end of the first oligonucleotide (precisely: the same amount of letters we had in each oligonucleotide in the given data).
    See: append_new_oligonucleotide – offsets for the new oligo
    :param names: list of oligonucleotides
    :return: matrix of offsets for each pair of them
    """

    l = len(names)
    offsets = np.empty([l, l], dtype=int)
    for i in range(l):
        for j in range(l):
            flag = False
            for k in range(len(names[j])):
                if names[i][k:] == names[j][:len(names[i]) - k]:
                    offsets[i][j] = k
                    flag = True
                    break
            if not flag:
                offsets[i][j] = len(names[0])
    return offsets


def append_new_oligonucleotide(names: List[str], offsets: np.ndarray, old_row: int, old_col: int, o_le: int) -> Tuple[np.ndarray, List[str]]:
    """
    Collage a pair of oligonucleotides according to their offset.
    Remove them; Save the result on place of old_col nucleotide.
    :param names: all oligonucleotides
    :param offsets: offset for each par of oligonucleotides
    :param old_row: first oligonucleotide to collage
    :param old_col: first oligonucleotide to collage
    :param o_le: langth of basic oligonucleotide
    :return: changed offsets and names matrices
    """
    print("\nCollage")

    offset = offsets[old_row][old_col]

    # offsets for the new oligo
    offsets[:, old_col] = offsets[:, old_row]
    offsets[old_col][old_col] = 0

    offsets = np.delete(offsets, old_row, 0)
    offsets = np.delete(offsets, old_row, 1)

    # collage names
    print(names[old_row] + '\n' + " "*(len(names[old_row]) - o_le + offset) + names[old_col])
    names[old_col] = names[old_row][:len(names[old_row]) - o_le + offset] + names[old_col][:]
    print(names[old_col] + '\n')

    names.pop(old_row)

    return offsets, names


def solve_conflict(conflict_row: int, offsets: np.ndarray, offset: int, names: list) -> int:
    """
    Choose oligonucleotide to collage in the row with > 1 the same offsets
    :param conflict_row: index of row with conflicts
    :param offsets: offset for each par of oligonucleotides
    :param offset: current offset (lowest)
    :param names: all oligonucleotides
    :return: index of oligonucleotide to collage or -1 when cannot solve the conflict
    """
    print("\nsolve conflict")
    print(f"{offsets[conflict_row]}\nrow {conflict_row}, offset {offset}")
    offsets_cols = np.nonzero(offsets[conflict_row] == offset)[0]

    # dłuższy z konfliktowych oligonukleotydów
    names_lengths = np.array([len(names[a]) for a in offsets_cols])
    longest_oligo = np.nonzero(names_lengths == max(names_lengths))[0]

    if len(longest_oligo) == 1:
        return offsets_cols[longest_oligo[0]]

    # z lepszym dopasowaniem w kolumnie
    mins_in_cols = []

    for i in range(len(offsets_cols)):
        val = min(a for a in list(offsets[:, offsets_cols[i]]) if a > offset)
        mins_in_cols.append(val)

    collageable_oligos = np.nonzero(mins_in_cols == (min(mins_in_cols)))[0]
    if len(collageable_oligos) > 1:
        return offsets_cols
    else:
        return collageable_oligos[0]


def solve_conflict_randomly(repeated_offsets_rows, offsets, offset):
    chosen_row = choice(list(repeated_offsets_rows))
    offsets_cols = np.nonzero(offsets[chosen_row] == offset)[0]
    chosen_col = choice(offsets_cols)
    return chosen_row, chosen_col
