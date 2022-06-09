from os import walk
from main import run_algorithm




if __name__ == "__main__":
    path = "./data"
    f_param = 1

    print("Działa")
    with open(f"results_{f_param}.csv", 'w', encoding='UTF-8') as out_file:
        out_file.write(f"z taboo")
        out_file.write("file_name;seq_len;density;difference_in_density;used_oligos;+/-;collager_time;taboo_time\n")
        for _, _, files in walk(path, topdown=True):
            for i, file_name in enumerate(files):
                print(f"nr {i}. plik {file_name}")
                seq_len, density, difference_in_density, used_oligos, collager_time, taboo_time = run_algorithm(file_name)
                print("Zapisano!")
                sign = '+' if file_name.find("+") > 0 else '-'

                out_file.write(f"{file_name};{seq_len};{density};{difference_in_density};{used_oligos};{sign};{collager_time};{taboo_time}\n")


