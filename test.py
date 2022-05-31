from os import walk
from STSP import Collager
import numpy as np
import logging

if __name__ == "__main__":
    path = "./data"
    print("Działa")
    with open("results", 'w', encoding='UTF-8') as out_file:
        for _, _, files in walk(path, topdown=True):
            for file_name in files:
                print(file_name)
                collager = Collager()
                collager.read_instance_from_file(path, file_name)
                collager.run_collager()
    
                logging.debug(f"{collager.names}\n{collager.offsets}")

                oligos_left = len([len(a) for a in collager.names])
                result = collager.names[np.argmax([len(a) for a in collager.names])]
                length = max([len(a) for a in collager.names])
                difference = length - (collager.SEQ_LENGTH+collager.OLIGO_LENGTH-1)

                out_file.write(f"{file_name};{oligos_left};{result};{length};{difference}")
                print(
                    f"\n---\nStan końcowy: {collager.names}\n{oligos_left}\n{length} {difference}\n\nPoszukiwana sekwencja to: {collager.names[np.argmax([len(a) for a in collager.names])]}")