from sekwencjonowanie import *

if __name__ == "__main__":
    names = []
    with open("krotki.txt", "r", encoding="UTF-8") as file:
        for line in file:
            names.append(line.strip("\n"))

    seq_length = 7
    offsets = get_offsets_matrix(names)
    oligo_length = len(names[1])
    lowest = 1
    row = 0
    col = 0
    oligos_left = len(offsets)
    repeated_offsets_rows = []
    print(offsets, names)
    sth_changed = False

    while (not any([len(x) >= seq_length for x in names])) and len(names) != 1:
        while row < oligos_left:
            if sum(offsets[row] == lowest) > 1:
                repeated_offsets_rows.append(row)
            else:
                while col < oligos_left:
                    if offsets[row][col] == lowest and row != col:
                        offsets = append_new_oligonucleotide(names, offsets, row, col, lowest, oligo_length)
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
            solve_conflict(repeated_offsets_rows.pop(), offsets, lowest)
        row = 0
        sth_changed = False
        print(len(names))
