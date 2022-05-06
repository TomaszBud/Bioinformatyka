import numpy as np


def get_neighbour_matrix(matrix):
    l = len(matrix)
    neighbours = np.empty([l, l])
    for i in range(l):
        for j in range(l):
            flag = False
            for k in range(len(matrix[j])):
                if matrix[i][k:] == matrix[j][:len(matrix[i]) - k]:
                    neighbours[i][j] = k
                    flag = True
                    break
            if not flag:
                neighbours[i][j] = len(matrix[0])
    return neighbours


def append_new_oligonucleotide(names, neighbours, old_row: int, old_col: int, offset: int, o_le: int):
    neighbours[:, old_col] = neighbours[:, old_row]
    neighbours[old_col][old_col] = 0
    neighbours = np.delete(neighbours, old_row, 0)
    neighbours = np.delete(neighbours, old_row, 1)

    # names[old_col] = names[old_row][:offset + len(names[old_row]) - o_le] + names[old_col]

    print(names[old_row] + '\n' + " "*(len(names[old_row]) - o_le + offset) + names[old_col])

    # print(names[old_row][:offset], names[old_col][:])
    names[old_col] = names[old_row][:len(names[old_row]) - o_le + offset] + names[old_col][:]
    print(names[old_col] + '\n')

    names.pop(old_row)
    return neighbours


# matrix = []
# with open("krotki.txt", "r", encoding="UTF-8") as file:
#     for line in file:
#         matrix.append(line.strip("\n"))
#
# print(matrix)
#
# neighbours = get_neighbour_matrix(matrix)
# for i in range(len(neighbours)):
#     print(matrix[i], neighbours[i])
# print(len(neighbours[0]) == len(neighbours[1]))
#
#
# ones = []
# ones_c = 0
# sumki = 0
# # for i in range(len(neighbours)):
# #     ones.append([])
# #     for j in range(len(neighbours[i])):
# #         if neighbours[i][j] == 1:
# #             ones_c += 1
# #             ones[i].append(1)
# #         else:
# #             ones[i].append(0)
# #     # print(sum(ones[i]))
# #     print(np.mean(neighbours[i]))
# print(ones_c, sumki)