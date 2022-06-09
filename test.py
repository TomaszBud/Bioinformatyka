from os import walk
from collager import Collager
import numpy as np
import logging
from Taboo import *


def run_algorithm(file_name, print=False):
    best_solution = ([], 0)

    collager = Collager()
    collager.read_instance_from_file(path, file_name)
    collager.run_collager()

    res = collager.names[np.argmax([len(a) for a in collager.names])]
    indexes_in_solution = collager.names_components[res]
    seq_len = max([len(a) for a in collager.names])
    density = len(indexes_in_solution) / len(res)

    if seq_len == collager.SEQ_LENGTH:
        best_solution = ([collager.BASE_NAMES[i] for i in indexes_in_solution], density)
        if len(collager.names) != 1:  # dla negatywnych
            return seq_len, density, 0

    taboo = Taboo(collager.BASE_OFFSETS, collager.BASE_NAMES, indexes_in_solution, collager.SEQ_LENGTH, best_solution)
    solution = taboo.run()

    taboo_seq = taboo.collage_sequence_from_solution(solution)

    return len(taboo_seq), taboo.calc_density(solution, len(taboo_seq)), (collager.SEQ_LENGTH - collager.OLIGO_LENGTH + 1) - len(solution)


if __name__ == "__main__":
    path = "./data"
    f_param = 1

    print("Działa")
    with open(f"results_{f_param}.csv", 'w', encoding='UTF-8') as out_file:
        out_file.write(f"z taboo")
        out_file.write("file_name;seq_len;density;difference;+/-\n")
        for _, _, files in walk(path, topdown=True):
            for i, file_name in enumerate(files):
                print(f"nr {i}. plik {file_name}")
                seq_len, density, difference = run_algorithm(file_name)
                print("Zapisano!")
                sign = '+' if file_name.find("+") > 0 else '-'

                out_file.write(f"{file_name};{seq_len};{density};{difference};{sign}\n")
