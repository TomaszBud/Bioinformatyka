from collager import Collager
from Taboo import Taboo
import time
import numpy as np
import logging


def run_algorithm(file_name, show=False):
    best_solution = ([], 0)
    path = "data/"

    collager = Collager()
    collager.read_instance_from_file(path, file_name)
    start = time.time()
    collager.run_collager()
    end = time.time()
    res = collager.names[np.argmax([len(a) for a in collager.names])]
    indexes_in_solution = collager.names_components[res]
    seq_len = max([len(a) for a in collager.names])
    density = len(indexes_in_solution) / len(res)
    collager_time = end-start

    if (pos := file_name.find("-")) > 0:
        desired_density = (collager.SEQ_LENGTH - collager.OLIGO_LENGTH + 1) - int(file_name[pos+1:]) / collager.SEQ_LENGTH
    else:
        desired_density = (collager.SEQ_LENGTH - collager.OLIGO_LENGTH + 1) / collager.SEQ_LENGTH

    if show:
        logging.debug(f"{collager.names}\n{collager.offsets}")
        print(f"\n---\nStan po collagerze:",
              f"długość sekwencji:{seq_len}",
              f"zagęszczenie: {density}\n",
              f"Poszukiwana sekwencja: {res}\n---\n", sep="\n")

    if seq_len == collager.SEQ_LENGTH:
        best_solution = ([collager.BASE_NAMES[i] for i in indexes_in_solution], density)
        if len(collager.names) != 1:  # dla negatywnych
            return seq_len, density, desired_density - density, len(indexes_in_solution), collager_time, 0

    taboo = Taboo(collager.BASE_OFFSETS, collager.BASE_NAMES, indexes_in_solution, collager.SEQ_LENGTH, best_solution)

    start = time.time()
    solution = taboo.run()
    end = time.time()

    taboo_seq = taboo.collage_sequence_from_solution(solution)
    density = taboo.calc_density(solution, len(taboo_seq))

    taboo_time = end - start

    if show:
        print(f"\n---\nStan po taboo:",
              f"długość sekwencji:{len(taboo_seq)}",
              f"zagęszczenie: {taboo.calc_density(solution, len(taboo_seq))}\n",
              f"Poszukiwana sekwencja: {taboo_seq}",
              f"Użyto {len(solution)} oligonukleotydów\n---\n", sep="\n")

        logging.info(
            f"odległość od pożądanego wyniku: {(collager.SEQ_LENGTH - collager.OLIGO_LENGTH + 1) - len(solution)}")

    return len(taboo_seq), density, desired_density - density, len(solution), collager_time, taboo_time

if __name__ == '__main__':
    run_algorithm('9.200+20', show=True)