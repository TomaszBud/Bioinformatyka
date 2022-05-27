from random import shuffle
import numpy as np
from typing import *
from random import choice


class Collager:
    def __init__(self):
        self.names = []
        self.offsets = None
        self.SEQ_LENGTH = None
        self.OLIGO_LENGTH = None
        self.col = 0
        self.row = 0
        self.oligos_left = None
        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()

    def read_instance_from_file(self, file_name):
        with open(file_name, "r", encoding="UTF-8") as file:
            for line in file:
                self.names.append(line.strip("\n"))

        self.SEQ_LENGTH = 200
        shuffle(self.names)
        self.get_offsets_matrix()
        self.OLIGO_LENGTH = len(self.names[1])
        self.oligos_left = len(self.names)

    def update_variables(self, sth_changed) -> bool:
        self.oligos_left = len(self.offsets)
        if self.row >= self.oligos_left:
            self.row = self.oligos_left - 1
        if self.col >= self.oligos_left:
            self.col = self.oligos_left - 1
        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()
        sth_changed = True

        return sth_changed

    def run_collager(self, long_collaging=False, forbidden_names=set()):
        lowest = 1
        sth_changed = False
        # stop when the sequence has desired length or no more oligos left
        while (not any([len(x) >= self.SEQ_LENGTH for x in self.names])) and len(self.names) != 1:
            while self.row < self.oligos_left and self.names[self.row] not in forbidden_names:
                if sum(self.offsets[self.row] == lowest) > 1:
                    self.repeated_offsets_rows.add(self.row)
                else:
                    # only one offset in row - process it
                    while self.col < self.oligos_left and self.names[self.col] not in forbidden_names:
                        if self.offsets[self.row][self.col] == lowest and self.row != self.col:
                            if sum(self.offsets[:, self.col] == lowest) > 1:
                                self.repeated_offsets_cols.add(self.col)
                            else:
                                self.append_new_oligonucleotide(self.row, self.col)
                                print(f"nowe offsety:\n{self.offsets}")
                                sth_changed = self.update_variables(sth_changed)
                                break
                        self.col += 1
                self.row += 1
                self.col = 0

            if not any([y == lowest for y in self.offsets.flatten()]):
                lowest += 1

                print(f"New lowest: {lowest}")

            if not sth_changed:
                # solve a conflict
                for conflicted_row in self.repeated_offsets_rows:
                    chosen_oligo = self.solve_conflict(conflicted_row, self.offsets, lowest)
                    print(f"chosen oligo: {chosen_oligo}")
                    if chosen_oligo > -1:
                        self.append_new_oligonucleotide(conflicted_row, chosen_oligo)
                        # repeated_offsets_rows.remove(conflicted_row)
                        sth_changed = self.update_variables( sth_changed)
                        break

            if not sth_changed:
                for conflicted_col in self.repeated_offsets_cols:
                    chosen_oligo = self.solve_conflict(conflicted_col, self.offsets.T, lowest)
                    print(f"chosen oligo: {chosen_oligo}")
                    if chosen_oligo > -1:
                        self.append_new_oligonucleotide(chosen_oligo, conflicted_col)
                        # repeated_offsets_rows.remove(conflicted_row)
                        sth_changed = self.update_variables( sth_changed)
                        break

                if not sth_changed and len(self.repeated_offsets_rows):
                    chosen_row, chosen_col = self.solve_conflict_randomly(self.offsets, lowest)
                    self.append_new_oligonucleotide(chosen_row, chosen_col)
                    # repeated_offsets_rows.remove(conflicted_row)
                    sth_changed = self.update_variables( sth_changed)

                if not sth_changed and len(self.repeated_offsets_cols):
                    chosen_col, chosen_row = self.solve_conflict_randomly(self.offsets.T, lowest)
                    self.append_new_oligonucleotide(chosen_row, chosen_col)
                    # repeated_offsets_rows.remove(conflicted_row)
                    sth_changed = self.update_variables(sth_changed)
            if lowest > self.OLIGO_LENGTH / 4 and not long_collaging:
                self.offsets, self.names = self.try_to_collage(self.offsets, lowest, self.SEQ_LENGTH)
                self.oligos_left = len(self.names)
                sth_changed = self.update_variables(sth_changed)

            self.row = 0
            sth_changed = False
            print(f"Oligos left: {len(self.names)}")
        return self.names

    def get_offsets_matrix(self):
        """
        Calculate offset for each pair of oligonucleotides.
        Offset is a number of letters (nucleotides) we have to move in the first oligonucleotide to achieve the beginning of the second one.
        Eg. CGTT, GTTG -> 1, because we have to move by only one letter – C.
        Later it is computed only for the end of the first oligonucleotide (precisely: the same amount of letters we had in each oligonucleotide in the given data).
        See: append_new_oligonucleotide – offsets for the new oligo
        :param names: list of oligonucleotides
        :return: matrix of offsets for each pair of them
        """

        l = len(self.names)
        self.offsets = np.empty([l, l], dtype=int)
        for i in range(l):
            for j in range(l):
                flag = False
                for k in range(len(self.names[j])):
                    if self.names[i][k:] == self.names[j][:len(self.names[i]) - k]:
                        self.offsets[i][j] = k
                        flag = True
                        break
                if not flag:
                    self.offsets[i][j] = len(self.names[0])

    def append_new_oligonucleotide(self, old_row: int, old_col: int):
        """
        Collage a pair of oligonucleotides according to their offset.
        Remove them; Save the result on place of old_col nucleotide.
        :param self.names: all oligonucleotides
        :param offsets: offset for each par of oligonucleotides
        :param old_row: first oligonucleotide to collage
        :param old_col: first oligonucleotide to collage
        :param o_le: langth of basic oligonucleotide
        :return: changed offsets and self.names matrices
        """
        print("\nCollage")

        offset = self.offsets[old_row][old_col]

        # self.offsets for the new oligo
        self.offsets[:, old_col] = self.offsets[:, old_row]
        self.offsets[old_col][old_col] = 0

        self.offsets = np.delete(self.offsets, old_row, 0)
        self.offsets = np.delete(self.offsets, old_row, 1)

        # collage self.names
        print(self.names[old_row] + '\n' + " " * (len(self.names[old_row]) - self.OLIGO_LENGTH + offset) + self.names[old_col])
        self.names[old_col] = self.names[old_row][:len(self.names[old_row]) - self.OLIGO_LENGTH + offset] + self.names[old_col][:]
        print(self.names[old_col] + '\n')

        self.names.pop(old_row)

    def solve_conflict(self, conflict_row: int, offsets: np.ndarray, offset: int) -> int:
        """
        Choose oligonucleotide to collage in the row with > 1 the same offsets
        :param conflict_row: index of row with conflicts
        :param offsets: offset for each par of oligonucleotides
        :param offset: current offset (lowest)
        :param self.names: all oligonucleotides
        :return: index of oligonucleotide to collage or -1 when cannot solve the conflict
        """
        print("\nsolve conflict")
        print(f"{offsets[conflict_row]}\nrow {conflict_row}, offset {offset}")
        offsets_cols = np.nonzero(offsets[conflict_row] == offset)[0]

        # dłuższy z konfliktowych oligonukleotydów
        names_lengths = np.array([len(self.names[a]) for a in offsets_cols])
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
            return -1
        else:
            return collageable_oligos[0]

    def solve_conflict_randomly(self, offsets, offset):
        chosen_row = choice(list(self.repeated_offsets_rows))
        offsets_cols = np.nonzero(offsets[chosen_row] == offset)[0]
        chosen_col = choice(offsets_cols)
        return chosen_row, chosen_col

    def try_to_collage(self, offsets, offset, length) -> Tuple[np.ndarray, list]:
        """
        Get the firsts longest oligos, which sum exceeds sequence lenght and (TODO: try to collage them).

        :param self.names: all oligonucleotides
        :param offsets: offset for each par of oligonucleotides
        :param offset: current offset (lowest)
        :param length: desired sequence lenght
        :return: slices of the offsets and self.names matrices
        """

        print("\n-----------------try to collage-------------")
        # sort self.names by their length
        names_order = sorted(range(len(self.names)), key=lambda i: len(self.names[i]), reverse=True)
        print(names_order)

        curr_len = 0
        how_many_oligos = 0
        flag = False

        # get the firsts longest oligos, which sum exceeds sequence lenght
        for i, name_index in enumerate(names_order):
            curr_len += len(self.names[name_index])
            if curr_len > length + offset * i:
                how_many_oligos += i + 1
                flag = True
                break


        if not flag:
            print(
                "ALARM!!! TU NIE POWINIEN BYć".upper())  # TODO: nie wejdzie, jeżeli będzie miał za mało oligo. Ma za mało, bo je usuwamy (patrz: slice matrices)
            how_many_oligos = len(self.names)
            return offsets, self.names

        # slice matrices
        print("wybrał", how_many_oligos)
        print("długości", [len(self.names[x]) for x in names_order[:how_many_oligos]])

        forbidden = set(np.delete(self.names.copy(), names_order[:how_many_oligos]))
        print("forbidden", forbidden)

        self.run_collager(True, forbidden)


if __name__ == "__main__":
    file_name = '200-40-2'
    collager = Collager()
    collager.read_instance_from_file(file_name)
    collager.run_collager()

    print(collager.names, collager.offsets, sep='\n')
    print(f"\n---\nStan końcowy: {collager.names}\n{[len(a) for a in collager.names]}\n{max([len(a) for a in collager.names])}\n\nPoszukiwana sekwencja to: {collager.names[np.argmax([len(a) for a in collager.names])]}")