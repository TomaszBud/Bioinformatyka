from sekwencjonowanie import *
from random import shuffle


def update_variables(oligos_left: int, offsets: np.ndarray, row: int, col: int, sth_changed: bool) -> Tuple[int, int, int, bool, set, set]:
    oligos_left = len(offsets)
    if row >= oligos_left:
        row = oligos_left - 1
    if col >= oligos_left:
        col = oligos_left - 1
    sth_changed = True

    return oligos_left, row, col, sth_changed, set(), set()


if __name__ == "__main__":
    names = []  # oligonucleotides

    with open("200-40-2", "r", encoding="UTF-8") as file:
        for line in file:
            names.append(line.strip("\n"))

    SEQ_LENGTH = 200
    shuffle(names)
    offsets = get_offsets_matrix(names)
    OLIGO_LENGTH = len(names[1])
    lowest = 1  # current offset
    row = 0
    col = 0
    oligos_left = len(offsets)
    repeated_offsets_rows = set()  # rows with > 1 element with current offset
    repeated_offsets_cols = set()
    print(names, offsets, sep='\n')
    sth_changed = False  # was this iteration successful

    # stop when the sequence has desired length or no more oligos left
    while (not any([len(x) >= SEQ_LENGTH for x in names])) and len(names) != 1:
        while row < oligos_left:
            if sum(offsets[row] == lowest) > 1:
                repeated_offsets_rows.add(row)
            else:
                # only one offset in row - process it
                while col < oligos_left:
                    if offsets[row][col] == lowest and row != col:
                        if sum(offsets[:,col] == lowest) > 1:
                            repeated_offsets_cols.add(col)
                        else:
                            offsets, names = append_new_oligonucleotide(names, offsets, row, col, OLIGO_LENGTH)
                            print(f"nowe offsety:\n{offsets}")
                            oligos_left, row, col, sth_changed, repeated_offsets_rows, repeated_offsets_cols = update_variables(oligos_left,                                                                         offsets, row, col, sth_changed)
                            break
                    col += 1
            row += 1
            col = 0

        if not any([y == lowest for y in offsets.flatten()]):
            lowest += 1

            print(f"New lowest: {lowest}")

        if not sth_changed:
            # solve a conflict
            for conflicted_row in repeated_offsets_rows:
                chosen_oligo = solve_conflict(conflicted_row, offsets, lowest, names)
                print(f"chosen oligo: {chosen_oligo}")
                if chosen_oligo > -1:
                    offsets, names = append_new_oligonucleotide(names, offsets, conflicted_row, chosen_oligo,
                                                                OLIGO_LENGTH)
                    # repeated_offsets_rows.remove(conflicted_row)
                    oligos_left, row, col, sth_changed, repeated_offsets_rows, repeated_offsets_cols = update_variables(oligos_left, offsets,
                                                                                                 row, col, sth_changed)
                    break

        if not sth_changed:
            for conflicted_col in repeated_offsets_cols:
                chosen_oligo = solve_conflict(conflicted_col, offsets.T, lowest, names)
                print(f"chosen oligo: {chosen_oligo}")
                if chosen_oligo > -1:
                    offsets, names = append_new_oligonucleotide(names, offsets, chosen_oligo, conflicted_col,
                                                                OLIGO_LENGTH)
                    # repeated_offsets_rows.remove(conflicted_row)
                    oligos_left, row, col, sth_changed, repeated_offsets_cols, repeated_offsets_rows = update_variables(oligos_left, offsets,
                                                                                                 row, col, sth_changed)
                    break

            if not sth_changed and len(repeated_offsets_rows):
                chosen_row, chosen_col = solve_conflict_randomly(repeated_offsets_rows, offsets, lowest)
                offsets, names = append_new_oligonucleotide(names, offsets, chosen_row, chosen_col,
                                                            OLIGO_LENGTH)
                # repeated_offsets_rows.remove(conflicted_row)
                oligos_left, row, col, sth_changed, repeated_offsets_rows, repeated_offsets_cols = update_variables(oligos_left, offsets, row,
                                                                                             col, sth_changed)

            if not sth_changed and len(repeated_offsets_cols):
                chosen_col, chosen_row = solve_conflict_randomly(repeated_offsets_cols, offsets.T, lowest)
                offsets, names = append_new_oligonucleotide(names, offsets, chosen_row, chosen_col,
                                                            OLIGO_LENGTH)
                # repeated_offsets_rows.remove(conflicted_row)
                oligos_left, row, col, sth_changed, repeated_offsets_rows, repeated_offsets_cols = update_variables(
                    oligos_left, offsets, row,
                    col, sth_changed)
        if lowest > OLIGO_LENGTH/4:
            offsets, names = try_to_collage(names, offsets, lowest, SEQ_LENGTH)
            oligos_left = len(names)

        row = 0
        sth_changed = False
        print(f"Oligos left: {len(names)}")

    print(f"\n---\nStan końcowy: {names}\n{[len(a) for a in names]}\n{max([len(a) for a in names])}\n\nPoszukiwana sekwencja to: {names[np.argmax([len(a) for a in names])]}")