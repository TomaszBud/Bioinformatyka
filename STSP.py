from random import shuffle
import numpy as np
from typing import *
from random import choice


class Collager:
    def __init__(self):
        self.names = []  # all oligonucleotides
        self.offsets = None  # offset for each par of oligonucleotides
        self.SEQ_LENGTH = None  # desired sequence length
        self.OLIGO_LENGTH = None  # length of basic oligonucleotide
        self.col = 0
        self.row = 0
        self.oligos_left = None
        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()

    def read_instance_from_file(self, file_name: str):
        with open(file_name, "r", encoding="UTF-8") as file:
            for line in file:
                self.names.append(line.strip("\n"))

        self.SEQ_LENGTH = int(file_name.split(".")[1].replace("-", ",").replace("+", ",").split(",")[0])

        shuffle(self.names)
        self.get_offsets_matrix()
        self.OLIGO_LENGTH = len(self.names[1])
        self.oligos_left = len(self.names)

    def update_variables(self) -> bool:
        print("czyszczenie")
        self.oligos_left = len(self.offsets)
        if self.row >= self.oligos_left:
            self.row = self.oligos_left - 1
        if self.col >= self.oligos_left:
            self.col = self.oligos_left - 1

        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()
        sth_changed = True
        return sth_changed

    def run_collager(self, extra_collaging=True, forbidden_names=set()):
        lowest = 1

        # stop when the sequence has desired length or no more oligos left
        while (not any([len(x) >= self.SEQ_LENGTH for x in self.names])) and len(self.names) != 1:
            sth_changed = False
            self.row = 0

            while self.row < self.oligos_left and self.names[self.row] not in forbidden_names:
                # look for conflicts
                if sum(self.offsets[self.row] == lowest) > 1:
                    self.repeated_offsets_rows.add(self.row)
                    print(f"Nowy konflikt w wierszu {self.row}")
                else:
                    # only one such offset in row
                    while self.col < self.oligos_left and self.names[self.col] not in forbidden_names:
                        if self.offsets[self.row][self.col] == lowest and self.row != self.col:
                            if sum(self.offsets[:, self.col] == lowest) > 1:
                                self.repeated_offsets_cols.add(self.col)
                                print(f"Nowy konflikt w kolumnie {self.col}")
                            else:
                                # only one such offset in col - process it
                                self.collage_oligonucleotides(self.row, self.col)
                                print(f"nowe offsety:\n{self.offsets}")
                                sth_changed = True
                                break
                        self.col += 1
                self.row += 1
                self.col = 0

            if not any([y == lowest for y in self.offsets.flatten()]):
                lowest += 1
                print(f"New lowest: {lowest}")

            # CONFLICTS SOLVING
            if not sth_changed:
                print(f"Solving conflicts: {self.repeated_offsets_rows}, {self.repeated_offsets_cols}")
                sth_changed = self.solve_conflicts(lowest, 'r')

                if not sth_changed:
                    sth_changed = self.solve_conflicts(lowest, 'c')

                    if not sth_changed and len(self.repeated_offsets_rows):
                        self.solve_conflict_randomly(lowest, 'r')
                        sth_changed = True

                    if not sth_changed and len(self.repeated_offsets_cols):
                        self.solve_conflict_randomly(lowest, 'c')
                        sth_changed = True

                    if not sth_changed:
                        print("Nic się nie zmieniło w tej iteracji")
                        lowest += 1
                        print(f"New lowest (oszukany): {lowest}")

            # if extra_collaging and lowest > self.OLIGO_LENGTH / 4:   #TODO: dobrać parametr
            #     print("extra_collaging")
            #     self.try_to_collage(lowest)
            #     sth_changed = self.update_variables()

            print(f"Oligos left: {len(self.names)}")
        return self.names

    def solve_conflicts(self, offset: int, axis='r'):
        """
        try to solve a conflict in from repeated offsets (in rows or columns)

        :param offset: current offset
        :param axis: 'r'/ 'c' – row or column
        :return:
        """
        if axis == 'r':
            conflicts = self.repeated_offsets_rows
            off = self.offsets
        else:
            conflicts = self.repeated_offsets_cols
            off = self.offsets.T

        for conflicted_item in conflicts:
            chosen_oligo = self.solve_conflict(conflicted_item, off, offset)
            print(f"chosen oligo: {chosen_oligo}")

            if chosen_oligo > -1:
                self.collage_oligonucleotides(conflicted_item, chosen_oligo)
                return True
        return False

    def get_offsets_matrix(self):
        """
        Calculate offset for each pair of oligonucleotides.
        Offset is a number of letters (nucleotides) we have to move in the first oligonucleotide to achieve the beginning of the second one.
        Eg. CGTT, GTTG -> 1, because we have to move by only one letter – C.
        Later it is computed only for the end of the first oligonucleotide (precisely: the same amount of letters we had in each oligonucleotide in the given data).
        See: append_new_oligonucleotide – offsets for the new oligo
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

    def solve_conflict(self, conflict_row: int, offsets: np.ndarray, offset: int) -> int:
        """
        Choose oligonucleotide to collage in the row with > 1 the same offsets
        :param conflict_row: index of row with conflicts
        :param offsets: offset for each par of oligonucleotides
        :param offset: current offset (lowest)
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

    def solve_conflict_randomly(self, offset: int, axis='r'):
        offsets = self.offsets if axis == 'r' else self.offsets.T
        chosen_row = choice(list(self.repeated_offsets_rows))
        offsets_cols = np.nonzero(offsets[chosen_row] == offset)[0]
        chosen_col = choice(offsets_cols)

        self.collage_oligonucleotides(chosen_row, chosen_col)

        return True

    def try_to_collage(self, offset) -> Tuple[np.ndarray, list]:
        """
        Get the firsts longest oligos, which sum exceeds sequence lenght and try to collage them.

        :param offset: current offset (lowest)
        :return:
        """

        print("\n-----------------try to collage-------------")
        # sort names by their length
        names_order = sorted(range(len(self.names)), key=lambda i: len(self.names[i]), reverse=True)
        print(names_order)

        curr_len = 0
        how_many_oligos = 0
        flag = False

        # get the firsts longest oligos, which sum exceeds sequence lenght
        for i, name_index in enumerate(names_order):
            curr_len += len(self.names[name_index])
            if curr_len > self.SEQ_LENGTH + offset * i:
                how_many_oligos += i + 1
                flag = True
                break

        if not flag:
            print("ALARM!!! TU NIE POWINIEN BYĆ")  # nie wejdzie, jeżeli będzie miał za mało oligo, żeby stworzyć pełną sekwencję
            # how_many_oligos = len(self.names)
            return None

        print(f"wybrał {how_many_oligos} oligo")
        print("Ich długości: ", [len(self.names[x]) for x in names_order[:how_many_oligos]])

        # tworzy listę nazw oligonukleotydów, których nie może użyć do składania sekwencji
        forbidden = set(np.delete(self.names.copy(), names_order[:how_many_oligos]))
        print("Forbidden", forbidden)

        self.run_collager(False, forbidden)

        #BUG: czasem są lepsze offsety wśród zakazanych oligo i nie zwiększa lowest, ale też nic nie może zrobić, bo wszystkie zbiory są puste


    def collage_oligonucleotides(self, old_row: int, old_col: int):
        """
        Collage a pair of oligonucleotides according to their offset.
        Remove them; Save the result on place of old_col nucleotide.
        :param old_row: first oligonucleotide to collage
        :param old_col: first oligonucleotide to collage
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
        print(self.names[old_row] + '\n' + " " * (len(self.names[old_row]) - self.OLIGO_LENGTH + offset) + self.names[
            old_col])
        self.names[old_col] = self.names[old_row][:len(self.names[old_row]) - self.OLIGO_LENGTH + offset] + self.names[
                                                                                                                old_col][
                                                                                                            :]
        print(self.names[old_col] + '\n')

        self.names.pop(old_row)

        self.update_variables()

        return True


if __name__ == "__main__":
    file_name = '9.200-40'
    collager = Collager()
    collager.read_instance_from_file(file_name)
    collager.run_collager()

    print(collager.names, collager.offsets, sep='\n')
    print(
        f"\n---\nStan końcowy: {collager.names}\n{[len(a) for a in collager.names]}\n{max([len(a) for a in collager.names])}\n\nPoszukiwana sekwencja to: {collager.names[np.argmax([len(a) for a in collager.names])]}")
