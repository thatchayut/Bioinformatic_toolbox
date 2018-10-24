import copy
import numpy as np

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

def meanCorrected(list, list_global_avg):
    matrix = copy.deepcopy(list)
    for row_index in range(0, len(matrix)):
        for col_index in range(0, len(matrix[row_index])):
            matrix[row_index][col_index] -= list_global_avg[col_index]
    # print("result = ")
    # print(matrix)
    return matrix
    
def covariance(list):
    matrix = copy.deepcopy(list)
    transpose_matrix = matrix.transpose()
    number_of_sample = matrix.shape[0]
    # print("number_of_sample = " + str(number_of_sample))
    # print("transpose_matrix")
    # print(transpose_matrix)
    # print(transpose_matrix.shape)
    # print("matrix")
    # print(matrix)
    # print(matrix.shape)
    result = np.matmul(transpose_matrix, matrix)
    # print("result")
    # print(result)
    # print(result.shape)

    for row_index in range(0, len(result)):
        for col_index in range(0, len(result[row_index])):
            result[row_index][col_index] /= number_of_sample
    return result

def poolCovariance(list_1, list_2, number_of_sample_all, number_of_sample_1, number_of_sample_2):
    matrix_1 = copy.deepcopy(list_1)
    matrix_2 = copy.deepcopy(list_2)
    multiplier_1 = number_of_sample_1 / number_of_sample_all
    multiplier_2 = number_of_sample_2 / number_of_sample_all
    result = np.zeros((matrix_1.shape[0], matrix_1.shape[1]))

    for row_index in range(0, matrix_1.shape[0]):
        print(row_index)
        for col_index in range(0, matrix_1.shape[1]):
            result[row_index, col_index] = (((multiplier_1) * matrix_1[row_index, col_index]) + ((multiplier_2) * matrix_2[row_index, col_index]))
    return result



