import pandas as pd
from random import shuffle


class ShrinkoExtenderoinator:
    def __init__(self, base_offsets, base_names, solution, seq_len, best_solution):
        self.BASE = pd.DataFrame(base_offsets, index=base_names, columns=base_names,
                                 dtype=int)  # offset matrix with cols and rows names = nucleotides names
        self.solution = [base_names[i] for i in solution]  # not joined oligo names included in solution
        self.OLIGO_LENGTH = len(base_names[0])  # length of basic oligo
        self.SEQ_LENGTH = seq_len  # desired seq lenght
        self.best_solution = best_solution
        self.best_density = (best_solution[1], len(best_solution[0]))

    def run(self) -> list:
        print("run")
        """
        :return: solution – names of oligonucleotides in sequence
        """
        max_iterations = 3

        density = 0
        for i in range(max_iterations):
            prev_density = density
            seq_len = self.shrink()

            # prepare offsets matrix
            offsets = self.BASE.copy()

            # delete banned oligos
            for t in self.solution[1:-1]:
                offsets.drop(t, axis=0, inplace=True)
                offsets.drop(t, axis=1, inplace=True)

            # offsets for the longest oligo
            longest = self.solution[0]

            offsets.loc[longest] = offsets.loc[self.solution[-1]]
            offsets[longest][longest] = 0
            offsets.drop(self.solution[-1], axis=1, inplace=True)
            offsets.drop(self.solution[-1], axis=0, inplace=True)

            seq_len = self.extend(longest, offsets)
            density = self.calc_density(self.solution, seq_len)

            if seq_len == self.SEQ_LENGTH and density > self.best_solution[1]:
                self.best_solution = (self.solution.copy(), density)  # remember it

            if density > self.best_density[0]:
                self.best_density = (density, seq_len)

        if self.best_solution[1] > density:
            return self.best_solution[0]  # back to better result
        else:
            return self.solution

    def shrink(self):
        """ delete oligos, which make density worse"""

        min_len = self.SEQ_LENGTH - self.OLIGO_LENGTH*2  # teraz pogarsza najlepsze rozwiązanie

        sth_changed = True
        seq_len = len(self.collage_sequence_from_solution(self.solution))
        current_density = 0

        while sth_changed:

            if len(self.solution) <= 1:
                print("Za krótkie rozwiązanie")
                break

            sth_changed = False
            current_density = self.calc_density(self.solution, seq_len)

            best_density = 0
            best_density_index = -1

            # find the one to eliminate
            x = [a for a in range(len(self.solution))]
            shuffle(x)
            for i in x:
                temp_seq_len = self.calc_seq_len(i, self.solution, seq_len)
                density = self.calc_density(self.solution[:i] + self.solution[i + 1:], temp_seq_len)
                if density >= best_density:
                    best_density = density
                    best_density_index = i

            # eliminate it
            if best_density > current_density or seq_len > min_len:  # or current_density < min_density
                seq_len = self.calc_seq_len(best_density_index, self.solution, seq_len)

                self.solution.pop(best_density_index)

                sth_changed = True
        # print(f"shrink {seq_len}")
        return seq_len

    def calc_seq_len(self, oligo_to_del_index, solution, seq_len):
        """
        Calc how solution length changes when oligo is deleted
        :param oligo_to_del_index:
        :param solution:
        :param seq_len:
        :return:
        """
        taboo_oligo = solution[oligo_to_del_index]

        if oligo_to_del_index == 0:
            next_oligo = self.solution[oligo_to_del_index + 1]
            return seq_len - self.BASE[next_oligo][taboo_oligo]
        elif oligo_to_del_index == len(solution) - 1:
            prev_oligo = self.solution[oligo_to_del_index - 1]
            return seq_len - self.BASE[taboo_oligo][prev_oligo]
        else:
            prev_oligo = self.solution[oligo_to_del_index - 1]
            next_oligo = self.solution[oligo_to_del_index + 1]
            seq_len -= self.BASE[taboo_oligo][prev_oligo] - (self.OLIGO_LENGTH - self.BASE[next_oligo][taboo_oligo])
            seq_len -= self.OLIGO_LENGTH - self.BASE[next_oligo][prev_oligo]
            return seq_len

    def extend(self, the_one: str, offsets: pd.DataFrame):
        """Add oligos to the_one to extend its lenght"""
        seq_len = len(self.collage_sequence_from_solution(self.solution))
        last_change = None

        # stop when the sequence has desired length or no more oligos left
        while (not seq_len >= self.SEQ_LENGTH) and offsets.shape[0] != 1:
            sth_changed = False
            row_without_main_diagonal = offsets.loc[the_one].drop(the_one, inplace=False)
            lowest_in_row = row_without_main_diagonal.idxmin()

            col_without_main_diagonal = offsets[the_one].drop(the_one, inplace=False)
            lowest_in_col = col_without_main_diagonal.idxmin()

            if offsets[the_one][lowest_in_col] < offsets[lowest_in_row][the_one]:
                # doklej offsets[the_one][lowest_in_col] na początek

                last_change = f'b{offsets[the_one][lowest_in_col]}'

                self.solution.insert(0, lowest_in_col)
                seq_len += offsets[the_one][lowest_in_col]
                self.change_offsets(lowest_in_col, the_one, offsets)


            else:
                # doklej offsets[lowest_in_row][the_one] na koniec

                last_change = f'e{offsets[lowest_in_row][the_one]}'
                self.solution.append(lowest_in_row)
                seq_len += offsets[lowest_in_row][the_one]

                self.change_offsets(the_one, lowest_in_row, offsets)
                the_one = lowest_in_row


        if seq_len > self.SEQ_LENGTH:
            if last_change[0] == 'b':
                self.solution.pop(0)
            if last_change[0] == 'e':
                self.solution.pop()
            seq_len -= int(last_change[1:])

        # TODO: usunąć
        if seq_len != len(self.collage_sequence_from_solution(self.solution)):
            print("znowu źle liczy długość...")

        # print(f"extend {seq_len}")
        return seq_len

    def collage_base_oligos(self, left: str, right: str):
        """collage two oligos based on their offset"""

        offset = self.BASE[right][left]
        new_oligo_name = left[:len(left) - self.OLIGO_LENGTH + offset] + right[:]

        return new_oligo_name

    def collage_sequence_from_solution(self, solution: list):
        """ having the list of names of oligonucleotides, collage them to one sequence """
        collaged_sequence = solution[0]

        for i in range(len(solution) - 1):
            new_fragment = self.collage_base_oligos(solution[i], solution[i + 1])
            collaged_sequence = collaged_sequence[:-self.OLIGO_LENGTH] + new_fragment

        return collaged_sequence

    def calc_density(self, solution: list, seq_len: int) -> float:
        """ how many oligos used for the sequence with this length """
        used_oligos = len(solution)
        collaged_sequence = seq_len

        return used_oligos / collaged_sequence

    @staticmethod
    def change_offsets(left: str, right: str, offsets: pd.DataFrame):
        """change offsets for the new oligo"""
        offsets[right] = offsets[left]
        offsets[right][right] = 0
        offsets.drop(left, axis=1, inplace=True)
        offsets.drop(left, axis=0, inplace=True)

        return offsets
