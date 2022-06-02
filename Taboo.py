import pandas as pd
import numpy as np
import logging
# x = pd.DataFrame(np.array(range(9)).reshape((3,3)), index=['a', 'b', 'c'], columns=['a', 'b', 'c'])

class Taboo:
    def __init__(self, base_offsets, base_names, solution, seq_len):
        self.BASE = pd.DataFrame(base_offsets, index=base_names, columns=base_names, dtype=int)
        self.solution = [base_names[i] for i in solution]  # nazwy oligo, niepołączone
        self.OLIGO_LENGTH = len(base_names[0])
        self.SEQ_LENGTH = seq_len
        self.taboo = []

        logging.basicConfig(format='%(message)s', level=logging.INFO)

        print(self.BASE)

    def run(self):
        self.shrink()
        # prepare offsets matrix
        offsets = self.BASE.copy()

        # delete banned oligos
        for t in self.solution[1:-1] + self.taboo:
            offsets.drop(t, axis=0, inplace=True)
            offsets.drop(t, axis=1, inplace=True)


        # offsets for the longest oligo
        longest = self.solution[0]

        offsets.loc[longest] = offsets.loc[self.solution[-1]]
        offsets.drop(self.solution[-1], axis=1, inplace=True)
        offsets.drop(self.solution[-1], axis=0, inplace=True)

        self.extend(longest, offsets)
        print("result:", self.calc_density(self.solution), offsets)

        return self.collage_sequence_from_solution(self.solution)

    def shrink(self):
        """ delete oligos, which make density worse"""

        print("shrink")
        sth_changed = True
        while sth_changed:
            sth_changed = False
            current_density = self.calc_density(self.solution)
            print(f"dens: {current_density}")
            best_density = current_density
            best_density_index = -1

            # find the one to eliminate
            for i in range(len(self.solution)):
                density = self.calc_density(self.solution[:i] + self.solution[i + 1:])
                if density >= best_density:
                    best_density = density
                    best_density_index = i

            # eliminate it
            if best_density_index != -1:
                taboo_oligo = self.solution[best_density_index]
                self.taboo.append(taboo_oligo)

                self.solution.pop(best_density_index)
                sth_changed = True

    def extend(self, the_one: str, offsets: pd.DataFrame):
        """Add oligos to the_one to extend its lenght"""
        lowest = 1
        seq_len = len(self.collage_sequence_from_solution(self.solution))

        # stop when the sequence has desired length or no more oligos left
        while (not seq_len >= self.SEQ_LENGTH) and offsets.shape[0] != 1:
            sth_changed = False

            row_without_main_diagonal = offsets.loc[the_one].drop(the_one)
            lowest_in_row = row_without_main_diagonal.idxmin()

            col_without_main_diagonal = offsets[the_one].drop(the_one)
            lowest_in_col = col_without_main_diagonal.idxmin()

            if offsets[the_one][lowest_in_col] > offsets[lowest_in_row][the_one]:
                # doklej na koniec
                self.solution.append(lowest_in_row)
                seq_len += offsets[lowest_in_row][the_one]
                self.change_offsets(the_one, lowest_in_row, offsets)


            else:
                # doklej na początek
                self.solution.insert(0, lowest_in_col)
                seq_len += offsets[the_one][lowest_in_col]
                self.change_offsets(lowest_in_col, the_one, offsets)

    def collage_base_oligos(self, left: str, right: str):
        """collage two oligos based on their offset"""
        offset = self.BASE[right][left]

        # logging.info('\n'+ left + '\n' + " " * (len(left) - self.OLIGO_LENGTH + offset) + right)

        new_oligo_name = left[:len(left) - self.OLIGO_LENGTH + offset] + right[:]

        # logging.info(new_oligo_name)
        return new_oligo_name

    def collage_sequence_from_solution(self, solution: list):
        """ having the list of names of oligonucleotides, collage them to one sequence """
        collaged_sequence = solution[0]

        for i in range(len(solution) - 1):
            new_fragment = self.collage_base_oligos(solution[i], solution[i + 1])
            collaged_sequence = collaged_sequence[:-self.OLIGO_LENGTH] + new_fragment

        return collaged_sequence

    def calc_density(self, solution: list) -> float:
        """ how many oligos used for the sequence with this length """
        # collage sequence by indexes

        used_oligos = len(solution)
        collaged_sequence = self.collage_sequence_from_solution(solution)

        # print(collaged_sequence + '\n' + " " * (len(collaged_sequence) - self.OLIGO_LENGTH) + new_fragment)

        return used_oligos / len(collaged_sequence)

    @staticmethod
    def change_offsets(left: str, right: str, offsets: pd.DataFrame):
        # self.offsets for the new oligo
        """change offsets and delete used oligos"""
        offsets.loc[right] = offsets.loc[left]
        offsets[right][right] = 0

        offsets.drop(columns=left, inplace=True)
        offsets.drop(rows=left, inplace=True)

        # offsets[:, right] = offsets[:, left]
        # offsets[right][right] = 0
        #
        # offsets = np.delete(offsets, left, 0)
        # offsets = np.delete(offsets, left, 1)
        # base_indexes_header.pop(left)

        return offsets
