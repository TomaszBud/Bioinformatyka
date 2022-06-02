import pandas as pd
import numpy as np
import logging

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


        # offsets for the longest one
        longest = self.solution[0]

        offsets.loc[longest] = offsets.loc[self.solution[-1]]
        offsets.drop(self.solution[-1], axis=1, inplace=True)
        offsets.drop(self.solution[-1], axis=0, inplace=True)

        # self.solution = self.extend(offsets, base_indexes_header, longest_index, self.solution)
        print(offsets)

        return self.collage_sequence_from_solution(self.solution)

    def shrink(self):
        print("shrink")
        sth_changed = True
        while sth_changed:
            sth_changed = False
            current_density = self.calc_density(self.solution)
            print(f"dens: {current_density}")
            best_density = current_density
            best_density_index = -1

            for i in range(len(self.solution)):
                density = self.calc_density(self.solution[:i] + self.solution[i + 1:])
                if density >= best_density:
                    best_density = density
                    best_density_index = i

            if best_density_index != -1:
                taboo_oligo = self.solution[best_density_index]
                self.taboo.append(taboo_oligo)

                self.solution.pop(best_density_index)
                sth_changed = True

    def collage_base_oligos(self, left: str, right: str):
        offset = self.BASE[right][left]

        # logging.info('\n'+ left + '\n' + " " * (len(left) - self.OLIGO_LENGTH + offset) + right)

        new_oligo_name = left[:len(left) - self.OLIGO_LENGTH + offset] + right[:]

        # logging.info(new_oligo_name)
        return new_oligo_name

    def collage_sequence_from_solution(self, solution):
        collaged_sequence = solution[0]

        for i in range(len(solution) - 1):
            new_fragment = self.collage_base_oligos(solution[i], solution[i + 1])
            collaged_sequence = collaged_sequence[:-self.OLIGO_LENGTH] + new_fragment

        return collaged_sequence

    def calc_density(self, solution: list) -> float:
        # collage sequence by indexes

        used_oligos = len(solution)
        collaged_sequence = self.collage_sequence_from_solution(solution)

        # print(collaged_sequence + '\n' + " " * (len(collaged_sequence) - self.OLIGO_LENGTH) + new_fragment)

        return used_oligos / len(collaged_sequence)




    def taboo_collage(self, left: int, right: int, offsets: np.ndarray, base_indexes_header: list):
        # self.offsets for the new oligo
        offsets[:, right] = offsets[:, left]
        offsets[right][right] = 0

        offsets = np.delete(offsets, left, 0)
        offsets = np.delete(offsets, left, 1)
        base_indexes_header.pop(left)

        return offsets


    def extend(self, offsets: np.ndarray, base_indexes_header: list, seq_index: int, solution):
        lowest = 1

        # stop when the sequence has desired length or no more oligos left
        while (not any([len(x) >= self.SEQ_LENGTH for x in self.names])) and len(self.names) != 1:
            sth_changed = False

            row_without_main_diagonal = np.hstack((offsets[seq_index][:seq_index], offsets[seq_index][seq_index + 1:]))
            lowest_in_row_index = np.argmin(row_without_main_diagonal)[0]

            col_without_main_diagonal = np.hstack((offsets[:seq_index][seq_index], offsets[seq_index + 1:][seq_index]))
            lowest_in_col_index = np.argmin(col_without_main_diagonal)[0]

            if offsets[lowest_in_col_index][seq_index] > offsets[seq_index][lowest_in_row_index]:
                # doklej na koniec
                solution.append(base_indexes_header[lowest_in_row_index])
                self.taboo_collage(seq_index, lowest_in_row_index, offsets, base_indexes_header)


            else:
                # doklej na początek
                solution.insert(0, base_indexes_header[lowest_in_col_index])
                self.taboo_collage(lowest_in_col_index, seq_index, offsets, base_indexes_header)

        return solution