import pandas as pd
import logging
import numpy as np
# x = pd.DataFrame(np.array(range(9)).reshape((3,3)), index=['a', 'b', 'c'], columns=['a', 'b', 'c'])


class Taboo:
    def __init__(self, base_offsets, base_names, solution, seq_len):
        self.BASE = pd.DataFrame(base_offsets, index=base_names, columns=base_names, dtype=int)  # offset matrix with cols and rows names = nucleotides names
        self.solution = [base_names[i] for i in solution]  # not joined oligo names included in solution
        self.OLIGO_LENGTH = len(base_names[0])  # length of basic oligo
        self.SEQ_LENGTH = seq_len   # desired seq lenght
        self.TABOO_TIME = 3     # how long can oligo be in taboo list
        self.taboo = [[], []]   # oligos' names and time left for them in taboo list

        logging.basicConfig(format='%(message)s', level=logging.DEBUG)

    def run(self):
        """
        taboo heuristics
        :return:
        """
        seq_len = 0

        while seq_len != self.SEQ_LENGTH or self.calc_density(self.solution, seq_len) < 0.70:  #TODO: co jak będzie dużo negatywnych?
            print("run")
            self.update_taboo()
            # print(self.taboo[1])

            seq_len = self.shrink()
            # print(self.taboo[1])

            # prepare offsets matrix
            offsets = self.BASE.copy()

            # delete banned oligos
            for t in self.solution[1:-1] + list(self.taboo[0]):
                offsets.drop(t, axis=0, inplace=True)
                offsets.drop(t, axis=1, inplace=True)


            # offsets for the longest oligo
            longest = self.solution[0]

            offsets.loc[longest] = offsets.loc[self.solution[-1]]
            offsets[longest][longest] = 0
            offsets.drop(self.solution[-1], axis=1, inplace=True)
            offsets.drop(self.solution[-1], axis=0, inplace=True)

            seq_len = self.extend(longest, offsets)
            print("result:", self.calc_density(self.solution, seq_len), "\nUsed oligos len: ", len(self.solution))

        return self.solution


    def update_taboo(self):
        """
        Decrease time left in taboo list for each oligo
        """
        i = 0
        while i < len(self.taboo[1]):
            self.taboo[1][i] -= 1
            if self.taboo[1][i] == 0:
                self.taboo[0].pop(i)
                self.taboo[1].pop(i)
            i += 1


    def shrink(self):
        """ delete oligos, which make density worse"""
        #TODO: shrink powinien bardziej skracać. Skrócił do 207, później extend do 214 i tak w kółko.
        min_len = self.SEQ_LENGTH-self.OLIGO_LENGTH # teraz pogarsza najlepsze rozwiązanie

        print("shrink")
        sth_changed = True
        seq_len = len(self.collage_sequence_from_solution(self.solution))
        while sth_changed or seq_len > min_len:
            sth_changed = False
            current_density = self.calc_density(self.solution, seq_len)

            best_density = 0
            best_density_index = -1

            # find the one to eliminate
            for i in range(len(self.solution)):
                temp_seq_len = self.calc_seq_len(i, self.solution, seq_len)
                density = self.calc_density(self.solution[:i] + self.solution[i + 1:], temp_seq_len)
                if density >= best_density:
                    best_density = density
                    best_density_index = i

            # eliminate it
            if best_density > current_density or seq_len > min_len:
                seq_len = self.calc_seq_len(best_density_index, self.solution, seq_len)
                self.taboo[0].append(self.solution[best_density_index])
                self.taboo[1].append(self.TABOO_TIME)

                self.solution.pop(best_density_index)

                sth_changed = True
        print(f"seq_len: {seq_len}, {len(self.collage_sequence_from_solution(self.solution))}")
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
        print("extend")
        seq_len = len(self.collage_sequence_from_solution(self.solution))

        # stop when the sequence has desired length or no more oligos left
        while (not seq_len >= self.SEQ_LENGTH) and offsets.shape[0] != 1:
            sth_changed = False
            # TODO: remember about last added oligo and delete it when seq_len exceeds SEQ_LENGTH
            row_without_main_diagonal = offsets.loc[the_one].drop(the_one, inplace=False)
            lowest_in_row = row_without_main_diagonal.idxmin()

            col_without_main_diagonal = offsets[the_one].drop(the_one, inplace=False)
            lowest_in_col = col_without_main_diagonal.idxmin()

            if offsets[the_one][lowest_in_col] < offsets[lowest_in_row][the_one]:
                # doklej offsets[the_one][lowest_in_col] na początek
                self.solution.insert(0, lowest_in_col)
                print("na poczatek", offsets[the_one][lowest_in_col], offsets[lowest_in_col][the_one])
                seq_len += offsets[the_one][lowest_in_col]
                self.change_offsets(lowest_in_col, the_one, offsets)
            else:
                # doklej offsets[lowest_in_row][the_one] na koniec
                self.solution.append(lowest_in_row)
                print("na koniec", offsets[lowest_in_row][the_one], offsets[the_one][lowest_in_row])
                seq_len += offsets[lowest_in_row][the_one]

                self.change_offsets(the_one, lowest_in_row, offsets)
                the_one = lowest_in_row


            print(seq_len, len(self.collage_sequence_from_solution(self.solution)))  #TODO: źle liczy
        return seq_len

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

        print(self.solution, collaged_sequence)
        return collaged_sequence

    def calc_density(self, solution: list, seq_len: int) -> float:
        """ how many oligos used for the sequence with this length """
        used_oligos = len(solution)
        collaged_sequence = seq_len

        return used_oligos / collaged_sequence

    @staticmethod
    def change_offsets(left: str, right: str, offsets: pd.DataFrame):
        # self.offsets for the new oligo
        """change offsets and delete used oligos"""
        offsets[right] = offsets[left]
        offsets[right][right] = 0
        offsets.drop(left, axis=1, inplace=True)
        offsets.drop(left, axis=0, inplace=True)

        return offsets
