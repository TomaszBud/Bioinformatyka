import numpy as np


def get_offsets_matrix(names: list) -> np.ndarray:
    """
    :param names:
    :return:
    """
    l = len(names)
    offsets = np.empty([l, l])
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


def append_new_oligonucleotide(names: list, offsets: np.ndarray, old_row: int, old_col: int, offset: int, o_le: int) -> np.ndarray:
    """
    :param names:
    :param offsets:
    :param old_row:
    :param old_col:
    :param offset:
    :param o_le:
    :return:
    """

    offsets[:, old_col] = offsets[:, old_row]
    offsets[old_col][old_col] = 0
    offsets = np.delete(offsets, old_row, 0)
    offsets = np.delete(offsets, old_row, 1)

    print(names[old_row] + '\n' + " "*(len(names[old_row]) - o_le + offset) + names[old_col])
    names[old_col] = names[old_row][:len(names[old_row]) - o_le + offset] + names[old_col][:]
    print(names[old_col] + '\n')

    names.pop(old_row)

    return offsets


def solve_conflict(conflict_row: int, offsets: np.ndarray, lowest: int, names: list) -> int:
    print(conflict_row, offsets[conflict_row], lowest)
    offsets_cols = np.where(offsets[conflict_row] == lowest)[0]

    # dłuższy z konfliktowych oligonukl
    names_lengths = [len(names[a]) for a in offsets_cols]
    longest_oligo = np.where(names_lengths == max(names_lengths))   #tu jest błąd
    print(offsets_cols, names_lengths, longest_oligo,  max(names_lengths))
    if len(longest_oligo) == 1:
        return offsets_cols[longest_oligo[0]]

    # z lepszym dopasowaniem w kolumnie
    mins_in_cols = []

    for i in range(len(offsets_cols)):
        print("val ", i, " ", offsets_cols[i])
        val = min(a for a in list(offsets[:, offsets_cols[i]]) if a > lowest)
        mins_in_cols.append(val)

    collageable_oligos = np.where(mins_in_cols == (min(mins_in_cols)))[0]
    if len(collageable_oligos) > 1:
        return -1
    else:
        return collageable_oligos[0]
