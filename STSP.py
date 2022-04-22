from sekwencjonowanie import get_neighbour_matrix, append_new_oligonucletide

matrix = []
with open("500+50.txt", "r", encoding="UTF-8") as file:
    for line in file:
        matrix.append(line.strip("\n"))
# print(matrix)


length = 200
cost = 0
gain = 0
neighbours = get_neighbour_matrix(matrix)
o_le = len(matrix[1])
# print(o_le)
lowest = 1
while any(matrix) < length:
    for i in range(len(neighbours)):
        for j in range(len(neighbours[i])):
            if neighbours[i][j] == lowest and i != j:
                print(len(matrix[i])-lowest, len(matrix[i]), lowest, matrix[j])
                matrix[i] = matrix[i] + matrix[j][len(matrix[i])-lowest:]
                matrix.pop(j)
                for k in range(len(neighbours)):
                    neighbours[k].pop(j)
                neighbours.pop(j)
                append_new_oligonucletide(matrix, matrix[i], neighbours, i, j)

    if not any(neighbours) == lowest:
        lowest += 1


print(matrix)
for i in range(len(neighbours)):
    print(matrix[i], neighbours[i])
    
