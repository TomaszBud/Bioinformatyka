from os import walk
from STSP import Collager
import numpy as np
import logging

if __name__ == "__main__":
    path = "./data"
    f_param = 1

    print("Działa")
    with open(f"results_{f_param}.csv", 'w', encoding='UTF-8') as out_file:
        out_file.write(f"z poprawiającą, f_param: {f_param}, ustawianie wiersza na nowo sklejony oligo")
        out_file.write("file_name;oligos_left;length;difference;+/-\n")
        for _, _, files in walk(path, topdown=True):
            for i, file_name in enumerate(files):
                print(f"nr {i}. plik {file_name}")
                collager = Collager()
                collager.read_instance_from_file(path, file_name)
                collager.run_collager(f_param=f_param)

                logging.debug(f"{collager.names}\n{collager.offsets}")


                oligos_left = len([len(a) for a in collager.names])
                result = collager.names[np.argmax([len(a) for a in collager.names])]
                length = max([len(a) for a in collager.names])
                difference = length - collager.SEQ_LENGTH
                sign = '+' if file_name.find("+") > 0 else '-'

                out_file.write(f"{file_name};{oligos_left};{length};{difference};{sign}\n")
                print(
                    f"\n---\nStan końcowy: {[len(a) for a in collager.names]}\n{oligos_left}\n{length} {difference}\n\nPoszukiwana sekwencja to: {collager.names[np.argmax([len(a) for a in collager.names])]}")