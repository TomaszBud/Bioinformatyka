import logging
from random import shuffle
import numpy as np
from typing import *
from random import choice
from Taboo import Taboo

class Collager:
    def __init__(self):
        self.names = []  # all oligonucleotides
        self.names_components = dict()
        self.BASE_NAMES = []
        self.offsets = None  # offset for each par of oligonucleotides
        self.BASE_OFFSETS = []
        self.SEQ_LENGTH = None  # desired sequence length
        self.OLIGO_LENGTH = None  # length of basic oligonucleotide
        self.col = 0
        self.row = 0
        self.oligos_left = None
        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()

    def read_instance_from_file(self, path, file_name: str):
        with open(f"{path}/{file_name}", "r", encoding="UTF-8") as file:
            for i, line in enumerate(file):
                oligo = line.strip("\n")
                self.BASE_NAMES.append(oligo)
        logging.info(self.names_components)


        shuffle(self.BASE_NAMES)
        self.names = self.BASE_NAMES.copy()
        for i, name in enumerate(self.BASE_NAMES):
            self.names_components[name] = [i]

        self.get_offsets_matrix()
        self.OLIGO_LENGTH = len(self.names[1])
        self.oligos_left = len(self.names)

        self.SEQ_LENGTH = int(file_name.split(".")[1].replace("-", ",").replace("+", ",").split(",")[0]) + self.OLIGO_LENGTH - 1


    def update_variables(self) -> bool:
        logging.debug("czyszczenie")
        self.oligos_left = len(self.offsets)
        if self.row >= self.oligos_left:
            self.row = self.oligos_left - 1
        if self.col >= self.oligos_left:
            self.col = self.oligos_left - 1

        self.repeated_offsets_cols = set()
        self.repeated_offsets_rows = set()
        sth_changed = True
        return sth_changed

    def run_collager(self, f_param=3, extra_collaging=True, forbidden_names=set()):
        lowest = 1

        # stop when the sequence has desired length or no more oligos left
        while (not any([len(x) >= self.SEQ_LENGTH for x in self.names])) and len(self.names) != 1 and len(self.names) > len(forbidden_names) + 1:
            sth_changed = False
            self.row = 0

            while self.row < self.oligos_left:
                # logging.debug( f"Row is forbidden: {self.names[self.row] in forbidden_names}. Row: {self.names[self.row]}.\nForbidden {forbidden_names}")
                if self.names[self.row] not in forbidden_names:
                    # is offset in this row?
                    lowest_in_row = (self.offsets[self.row] == lowest)
                    nr_of_lowest_in_row = sum(lowest_in_row)

                    if nr_of_lowest_in_row > 1:
                        self.repeated_offsets_rows.add(self.row)
                        logging.debug(f"Nowy konflikt w wierszu {self.row}")

                    elif nr_of_lowest_in_row == 1:
                        # only one such offset in row - check column
                        self.col = np.argwhere(lowest_in_row)[0][0]

                        if self.names[self.col] not in forbidden_names and self.row != self.col:
                            lowest_in_col = self.offsets[:, self.col] == lowest
                            nr_of_lowest_in_col = sum(lowest_in_col)

                            if nr_of_lowest_in_col > 1:
                                self.repeated_offsets_cols.add(self.col)
                                logging.debug(f"Nowy konflikt w kolumnie {self.col}")
                            else:
                                # only one such offset in col - process it
                                self.collage_oligonucleotides(self.row, self.col)
                                logging.debug(f"nowe offsety:\n{self.offsets}")
                                sth_changed = True
                                break

                self.row += 1
                self.col = 0

                """
                #NOTE: problem z tym rozwiązaniem:
                przy [[1,2,1],
                      [1,0,1]]
                wykryje tylko konflikty w wierszach. Czy to nam przeszkadza?
                """

            # CONFLICTS SOLVING
            if not sth_changed:
                sth_changed = self.solve_conflicts(lowest)
                if not sth_changed: logging.debug("Nic się nie zmieniło w tej iteracji")

            if extra_collaging and lowest > self.OLIGO_LENGTH / f_param:   #TODO: dobrać parametr
                logging.debug("extra_collaging")
                self.try_to_collage(lowest)
                self.update_variables()

            if not any([y == lowest for y in self.offsets.flatten()]) or not sth_changed:
                lowest += 1
                logging.debug(f"New offset: {lowest}")

            logging.debug(f"Oligos left: {len(self.names)}\n")

        return self.names

    def solve_conflicts(self, offset: int):
        logging.debug(f"Solving conflicts: {self.repeated_offsets_rows}, {self.repeated_offsets_cols}")
        if len(self.repeated_offsets_rows) > 0:
            if self.solve_conflicts_by_axis(offset, 'r'):
                return True
        if len(self.repeated_offsets_cols) > 0:
            if self.solve_conflicts_by_axis(offset, 'c'):
                return True
        if len(self.repeated_offsets_rows) > 0:
            self.solve_conflict_randomly(offset, 'r')
            return True
        if len(self.repeated_offsets_cols) > 0:
            self.solve_conflict_randomly(offset, 'c')
            return True

        return False

    def solve_conflicts_by_axis(self, offset: int, axis='r'):
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
            logging.debug(f"chosen oligo: {chosen_oligo}")

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
        self.BASE_OFFSETS = np.empty([l, l], dtype=int)
        for i in range(l):
            for j in range(l):
                flag = False
                for k in range(len(self.names[j])):
                    if self.names[i][k:] == self.names[j][:len(self.names[i]) - k]:
                        self.BASE_OFFSETS[i][j] = k
                        flag = True
                        break
                if not flag:
                    self.BASE_OFFSETS[i][j] = len(self.names[0])

        self.offsets = self.BASE_OFFSETS.copy()

    def solve_conflict(self, conflict_row: int, offsets: np.ndarray, offset: int) -> int:
        """
        Choose oligonucleotide to collage in the row with > 1 the same offsets
        :param conflict_row: index of row with conflicts
        :param offsets: offset for each par of oligonucleotides
        :param offset: current offset (offset)
        :return: index of oligonucleotide to collage or -1 when cannot solve the conflict
        """
        logging.debug("\nsolve conflict")
        logging.debug(f"{offsets[conflict_row]}\nrow {conflict_row}, offset {offset}")
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
        """
        Choose random conflict and solve it.
        :param offset: current offset (offset)
        :param axis: 'r'/ 'c' – row or column
        :return:
        """
        if axis == 'r':
            offsets = self.offsets
            repeated = self.repeated_offsets_rows
        else:
            offsets = self.offsets.T
            repeated = self.repeated_offsets_cols

        chosen_row = choice(list(repeated))
        offsets_cols = np.nonzero(offsets[chosen_row] == offset)[0]
        chosen_col = choice(offsets_cols)

        self.collage_oligonucleotides(chosen_row, chosen_col)

        return True

    def try_to_collage(self, offset) -> Tuple[np.ndarray, list]:
        """
        Get the firsts longest oligos, which sum exceeds sequence lenght and try to collage them.

        :param offset: current offset (offset)
        :return:
        """

        logging.debug("\n-----------------try to collage-------------")
        # sort names by their length
        names_order = sorted(range(len(self.names)), key=lambda i: len(self.names[i]), reverse=True)
        logging.debug(names_order)

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
            logging.error("ALARM!!! TU NIE POWINIEN BYĆ")  # nie wejdzie, jeżeli będzie miał za mało oligo, żeby stworzyć pełną sekwencję
            # how_many_oligos = len(self.names)
            return None

        logging.debug(f"wybrał {how_many_oligos} oligo")
        logging.debug(f"Ich długości: {[len(self.names[x]) for x in names_order[:how_many_oligos]]}")

        # tworzy listę nazw oligonukleotydów, których nie może użyć do składania sekwencji
        forbidden = set(np.delete(self.names.copy(), names_order[:how_many_oligos]))
        logging.debug(f"Forbidden {forbidden}")

        self.run_collager(extra_collaging=False, forbidden_names=forbidden)

        #BUG: czasem są lepsze offsety wśród zakazanych oligo i nie zwiększa offset, ale też nic nie może zrobić, bo wszystkie zbiory są puste

    def collage_oligonucleotides(self, old_row: int, old_col: int):
        """
        Collage a pair of oligonucleotides according to their offset.
        Remove them; Save the result on place of old_col nucleotide.
        :param old_row: first oligonucleotide to collage
        :param old_col: second oligonucleotide to collage
        :return: changed offsets and self.names matrices
        """

        logging.debug("\nCollage")

        # left_base_name = self.BASE_NAMES.index(self.names[old_row][-self.OLIGO_LENGTH:])

        left_name = self.names[old_row]
        right_name = self.names[old_col]

        offset = self.offsets[old_row][old_col]

        # self.offsets for the new oligo
        self.offsets[:, old_col] = self.offsets[:, old_row]
        self.offsets[old_col][old_col] = 0

        self.offsets = np.delete(self.offsets, old_row, 0)
        self.offsets = np.delete(self.offsets, old_row, 1)

        # collage self.names
        logging.debug(left_name + '\n' + " " * (len(left_name) - self.OLIGO_LENGTH + offset) + right_name)

        new_oligo_name = left_name[:len(left_name) - self.OLIGO_LENGTH + offset] + right_name[:]
        self.names[old_col] = new_oligo_name
        logging.debug(right_name + '\n')
        self.names.pop(old_row)


        # collage names components
        self.names_components[new_oligo_name] = self.names_components[left_name] + self.names_components[right_name]

        if len(left_name) > self.OLIGO_LENGTH:
            self.names_components.pop(left_name)
        if len(right_name) > self.OLIGO_LENGTH:
            self.names_components.pop(right_name)

        self.update_variables()
        if old_row < old_col:
            self.row = old_col - 1
        else:
            self.row = old_col

        return True




if __name__ == "__main__":
    # logging.basicConfig(format='%(message)s', level=logging.DEBUG) #DEBUG, INFO, WARNING, ERROR, CRITICAL
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    path = "./data"
    file_name = '9.200+20'
    best_solution = ([], 0)

    collager = Collager()
    collager.read_instance_from_file(path, file_name)
    collager.run_collager()

    res = collager.names[np.argmax([len(a) for a in collager.names])]
    indexes_in_solution = collager.names_components[res]
    seq_len = max([len(a) for a in collager.names])
    density = len(indexes_in_solution)/ len(res)

    logging.debug(f"{collager.names}\n{collager.offsets}")
    print(f"\n---\nStan po collagerze:",
          f"długość sekwencji:{seq_len}",
          f"zagęszczenie: {density}\n",
          f"Poszukiwana sekwencja: {res}\n---\n", sep="\n")

    if seq_len == collager.SEQ_LENGTH:
        best_solution = ([collager.BASE_NAMES[i] for i in indexes_in_solution], density)
        if len(collager.names) != 1:  # dla negatywnych
            exit(0)

    # print(collager.calc_density(collager.names_components[res]))

    taboo = Taboo(collager.BASE_OFFSETS, collager.BASE_NAMES, indexes_in_solution, collager.SEQ_LENGTH, best_solution)
    solution = taboo.run()
    # collager.taboo(collager.names_components[res])

    taboo_seq = taboo.collage_sequence_from_solution(solution)
    print(f"\n---\nStan po taboo:",
          f"długość sekwencji:{len(taboo_seq)}",
          f"zagęszczenie: {taboo.calc_density(solution, len(taboo_seq))}\n",
          f"Poszukiwana sekwencja: {taboo_seq}",
          f"Użyto {len(solution)} oligonukleotydów\n---\n", sep="\n")

    logging.info(f"odległość od pożądanego wyniku: {(collager.SEQ_LENGTH-collager.OLIGO_LENGTH+1) - len(solution)}")
