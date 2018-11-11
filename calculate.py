import copy
import numpy as np

def avgFromList(list):
    #transpose matrix to calculate average of each column
    transposed_list = list.transpose()
    number_of_sample = transposed_list.shape[1]
    list_sum = transposed_list.sum(axis=1)

    for element in list_sum:
        element /= number_of_sample

    # transpose the matrix back to be the actual output
    list_sum = list_sum.transpose()
    return list_sum

def meanCorrected(list, list_global_avg):
    matrix = copy.deepcopy(list)
    for row_index in range(0, len(matrix)):
        for col_index in range(0, len(matrix[row_index])):
            matrix[row_index][col_index] -= list_global_avg[col_index]
    return matrix
    
def covariance(list):
    matrix = copy.deepcopy(list)
    transpose_matrix = matrix.transpose()
    number_of_sample = matrix.shape[0]
    result = np.matmul(transpose_matrix, matrix)

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
        for col_index in range(0, matrix_1.shape[1]):
            result[row_index, col_index] = (((multiplier_1) * matrix_1[row_index, col_index]) + ((multiplier_2) * matrix_2[row_index, col_index]))
    return result

def inversed(matrix):
    inversed_matrix = np.linalg.inv(matrix)
    return inversed_matrix

#  If we do not know the prior probability, we just assume it is equal to total sample of each group divided by the total samples
def findPriorProb(number_of_sample_all, number_of_sample_1, number_of_sample_2):
    arg_1 = number_of_sample_1 / number_of_sample_all
    arg_2 = number_of_sample_2 / number_of_sample_all
    prior_prob = np.matrix([[arg_1], [arg_2]])
    return prior_prob

# calculate discriminative function
def findDiscriminative(matrix, avg_matrix, inversed_pool_covariance, prior_prob):
    list_result = []
    # print(prior_prob)
    for row_index in range(0, matrix.shape[0]):
        arg_1 = avg_matrix.dot(inversed_pool_covariance)
        arg_1 = arg_1.dot(matrix[row_index].transpose())
        arg_2 = np.matmul(avg_matrix, inversed_pool_covariance)
        arg_2 = np.matmul(arg_2, avg_matrix.transpose())
        result = arg_1 - ((1 / 2) * arg_2) + np.log(prior_prob)
        # print((np.log(prior_prob)))
        # print(result.shape)
        list_result.append(result)
    return list_result

def findOutput(f1, f2):
    list_f1 = []
    for element_index  in range(0, len(f1)):
        for row_index in range(0, (f1[element_index].shape[0])):
            for col_index in range(0, (f1[element_index].shape[1])):
                # print(f1[element_index][row_index, col_index])
                list_f1.append(f1[element_index][row_index, col_index])
    # print(list_f1)

    list_f2 = []
    for element_index  in range(0, len(f2)):
        for row_index in range(0, (f2[element_index].shape[0])):
            for col_index in range(0, (f1[element_index].shape[1])):
                # print(f2[element_index][row_index, col_index])
                list_f2.append(f2[element_index][row_index, col_index])
    # print(list_f2)

    # check which class is selected
    result =[]
    for element_index in range(0, len(list_f1)):
        if(list_f1[element_index] > list_f2[element_index]):
            result.append(1)
        elif(list_f1[element_index] < list_f2[element_index]):
            result.append(0)
        else:
            result.append("*")
    # print(result)
    return result

# chunk data(l) into n folds
def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

