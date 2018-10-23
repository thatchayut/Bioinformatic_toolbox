def avgFromList(list):
    #transpose matrix to calculate average of each column
    transposed_list = list.transpose()
    # print(transposed_list)
    number_of_sample = transposed_list.shape[1]
    list_sum = transposed_list.sum(axis=1)

    # print("Before calculating ...")
    # print(list_sum)
    for element in list_sum:
        element /= number_of_sample
    # print("After calculating ...")
    # print(list_sum)

    # transpose the matrix back to be the actual output
    list_sum = list_sum.transpose()
    # print(list_sum)
    return list_sum


