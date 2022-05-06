from sekwencjonowanie import get_neighbour_matrix, append_new_oligonucleotide

if __name__ == "__main__":
    matrix = []
    with open("200-40-2", "r", encoding="UTF-8") as file:
        for line in file:
            matrix.append(line.strip("\n"))
# print(matrix)


    length = 209
    cost = 0
    gain = 0
    neighbours = get_neighbour_matrix(matrix)
    o_le = len(matrix[1])
    # print(o_le)
    lowest = 1
    i = 0
    j = 0
    leng = len(neighbours)
    print(neighbours, matrix)
    while (not any([len(x) >= length for x in matrix])) and len(matrix) != 1:
        while i < leng:
            while j < leng:
                if neighbours[i][j] == lowest and i != j:
                    # print(len(matrix[i])-lowest, len(matrix[i]), lowest, matrix[j])
                    neighbours = append_new_oligonucleotide(matrix, neighbours, i, j, lowest, o_le)
                    leng = len(neighbours)
                    if i >= leng:
                        i = leng - 1
                    if j >= leng:
                        j = leng - 1
                    # neighbours = get_neighbour_matrix(matrix)
                    # print(neighbours, "\n", matrix)
                j += 1
            i += 1
            j = 0
        if not any([y == lowest for y in neighbours.flatten()]):
            lowest += 1
            print(f"Lowest {lowest}")

        i = 0
        print(len(matrix))




    # [x for x in range(10)]

    # print(neighbours)


    # print(matrix)
    # for i in range(len(neighbours)):
    #     print(matrix[i], neighbours[i])

