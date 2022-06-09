from os import walk
from main import run_algorithm
import numpy as np

if __name__ == "__main__":
    path = "./data"
    f_param = 1
    seq_len = [0]*10
    density = [0]*10
    difference_in_density = [0]*10
    used_oligos = [0]*10
    collager_time = [0]*10
    se_time = [0] * 10
    print("Działa")
    with open(f"results_final.csv", 'w', encoding='UTF-8') as out_file:
        out_file.write(f"z taboo\n")
        out_file.write("file_name;seq_len;best_density;density_max;density_avg;density_std;difference_in_density;used_oligos;+/-;collager_time;se_time\n")
        for _, _, files in walk(path, topdown=True):
            for i, file_name in enumerate(files):
                for j in range(10):
                    print(f"nr {i}. plik {file_name}")
                    seq_len[i], density[i], difference_in_density[i], used_oligos[i], collager_time[i], se_time[i], tab_dens = run_algorithm(file_name)
                    print("Zapisano!")
                sign = '+' if file_name.find("+") > 0 else '-'
                density_max = max(density)
                density_avg = np.avg(density)
                density_std = np.std(density)
                out_file.write(f"{file_name};{tab_dens};{max(seq_len)};{density_max};{density_avg};{density_std};{max(difference_in_density)};{max(used_oligos)};{sign};{min(collager_time)};{minse_time}\n")


