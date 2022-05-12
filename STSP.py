from sekwencjonowanie import *


def update_variables(oligos_left: int, offsets: np.ndarray, row: int, col: int, sth_changed: bool) -> Tuple[int, int, int, bool]:
    oligos_left = len(offsets)
    if row >= oligos_left:
        row = oligos_left - 1
    if col >= oligos_left:
        col = oligos_left - 1
    sth_changed = True

    return oligos_left, row, col, sth_changed


if __name__ == "__main__":
    names = []  # oligonucleotides

    with open("krotki.txt", "r", encoding="UTF-8") as file:
        for line in file:
            names.append(line.strip("\n"))

    SEQ_LENGTH = 7
    offsets = get_offsets_matrix(names)
    OLIGO_LENGTH = len(names[1])
    lowest = 1  # current offset
    row = 0
    col = 0
    oligos_left = len(offsets)
    repeated_offsets_rows = set()  # rows with > 1 element with current offset
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
                        offsets, names = append_new_oligonucleotide(names, offsets, row, col, OLIGO_LENGTH)
                        print(f"nowe offsety:\n{offsets}")
                        oligos_left, row, col, sth_changed = update_variables(oligos_left, offsets, row, col, sth_changed)
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
                    offsets, names = append_new_oligonucleotide(names, offsets, conflicted_row, chosen_oligo, OLIGO_LENGTH)
                    repeated_offsets_rows.remove(conflicted_row)
                    oligos_left, row, col, sth_changed = update_variables(oligos_left, offsets, row, col, sth_changed)
                    break

        row = 0
        sth_changed = False
        print(f"Oligos left: {len(names)}")

    print(f"\n---\nStan końcowy: {names}")
