from sekwencjonowanie import *

if __name__ == "__main__":
    names = []
    with open("krotki.txt", "r", encoding="UTF-8") as file:
        for line in file:
            names.append(line.strip("\n"))

    SEQ_LENGTH = 7
    offsets = get_offsets_matrix(names)
    OLIGO_LENGTH = len(names[1])
    lowest = 1
    row = 0
    col = 0
    oligos_left = len(offsets)
    repeated_offsets_rows = []
    print(offsets, names)
    sth_changed = False

    while (not any([len(x) >= SEQ_LENGTH for x in names])) and len(names) != 1:
        while row < oligos_left:
            if sum(offsets[row] == lowest) > 1:
                repeated_offsets_rows.append(row)
            else:
                while col < oligos_left:
                    if offsets[row][col] == lowest and row != col:
                        offsets = append_new_oligonucleotide(names, offsets, row, col, lowest, OLIGO_LENGTH)
                        oligos_left = len(offsets)
                        if row >= oligos_left:
                            row = oligos_left - 1
                        if col >= oligos_left:
                            col = oligos_left - 1
                        sth_changed = True
                    # print("debug1")


                    col += 1
            row += 1
            col = 0
            # print("out")
        if not any([y == lowest for y in offsets.flatten()]):
            lowest += 1
            print(f"Lowest {lowest}")
        if not sth_changed:
            print(offsets)
            for i, conflicted_row in enumerate(repeated_offsets_rows):
                chosen_oligo = solve_conflict(conflicted_row, offsets, lowest, names)
                if chosen_oligo > -1:
                    append_new_oligonucleotide(names, offsets, conflicted_row, chosen_oligo, lowest, OLIGO_LENGTH)
                    del repeated_offsets_rows[i]

        row = 0
        sth_changed = False
        print(f"len(names) {len(names)}")
