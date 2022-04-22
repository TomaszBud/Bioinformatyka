import numpy as np
def get_neighbour_matrix(matrix):
    neighbours = []
    for i in range(len(matrix)):
        neighbours.append([])
        for j in range(len(matrix)):
            flag = False
            for k in range(len(matrix[i])):
                if matrix[i][k:len(matrix[i])] == matrix[j][0:len(matrix[i]) - k]:
                    neighbours[i].append(k)
                    flag = True
                    break
            if not flag:
                neighbours[i].append(10)

    return neighbours


def append_new_oligonucletide(matrix, oligonucleotide, neighbours, old1, old2):
    for j in range(len(matrix)):
        flag = False
        for k in range(len(matrix[old1])):
            if matrix[old1][k:len(matrix[old1])] == matrix[j][0:len(matrix[old1]) - k]:
                neighbours[old1][old2] = k
                neighbours[old2][old1] = k
                flag = True
                break
        if not flag:
            neighbours[old1].append(10)





matrix = []
with open("500+50.txt", "r", encoding="UTF-8") as file:
    for line in file:
        matrix.append(line.strip("\n"))

print(matrix)

neighbours = get_neighbour_matrix(matrix)
for i in range(len(neighbours)):
    print(matrix[i], neighbours[i])
print(len(neighbours[0]) == len(neighbours[1]))


ones = []
ones_c = 0
sumki = 0
for i in range(len(neighbours)):
    ones.append([])
    for j in range(len(neighbours[i])):
        if neighbours[i][j] == 2:
            ones_c += 1
            ones[i].append(1)
        else:
            ones[i].append(0)
    # print(sum(ones[i]))
    print(np.mean(neighbours[i]))
print(ones_c, sumki)