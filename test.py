from os import walk
from main import run_algorithm
import numpy as np

if __name__ == "__main__":
    path = "./our_data"

    print("Działa")
    with open(f"results_five_our.csv", 'w', encoding='UTF-8') as out_file:
        out_file.write(f"z taboo\n")
        out_file.write("file_name;seq_len;density;best_density;difference_in_density;used_oligos;+/-;collager_time;se_time\n")
        for _, _, files in walk(path, topdown=True):
            for i, file_name in enumerate(files):
                for j in range(5):
                    print(f"nr {i}. plik {file_name}")
                    seq_len, density, difference_in_density, used_oligos, collager_time, se_time, tab_dens = run_algorithm(file_name, path="our_data/")
                    print("Zapisano!")
                    sign = '+' if file_name.find("+") > 0 else '-'
                    out_file.write(f"{file_name};{seq_len};{density};{tab_dens};{difference_in_density};{used_oligos};{sign};{collager_time};{se_time}\n")

    # with open(f"results_ten.csv", 'w', encoding='UTF-8') as out_file:
    #     print("Działa")
    #     out_file.write(f"z taboo\n")
    #     out_file.write(
    #         "file_name;seq_len;best_density;density_max;density_avg;density_std;difference_in_density;used_oligos;+/-;collager_time;se_time\n")
    #     for _, _, files in walk(path, topdown=True):
    #         for i, file_name in enumerate(files):
    #             seq_len = [0] * 10
    #             density = [0] * 10
    #             difference_in_density = [0] * 10
    #             used_oligos = [0] * 10
    #             collager_time = [0] * 10
    #             se_time = [0] * 10
    #             for j in range(10):
    #                 print(f"nr {i}. plik {file_name}")
    #                 seq_len[j], density[j], difference_in_density[j], used_oligos[j], collager_time[j], se_time[
    #                     j], tab_dens = run_algorithm(file_name)
    #                 print("Zapisano!")
    #             sign = '+' if file_name.find("+") > 0 else '-'
    #             density_max = max(density)
    #             density_avg = np.mean(density)
    #             density_std = np.std(density)
    #             out_file.write(f"{file_name};{max(seq_len)};{tab_dens};{density_max};{density_avg};{density_std};{min(difference_in_density)};{max(used_oligos)};{sign};{np.mean(collager_time)};{np.mean(se_time)}\n")