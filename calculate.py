import copy
import numpy as np
import pandas as pd
import math

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
            # result.append("*")
            result.append(0.5)
    # print(result)
    return result

# split data l into n folds
def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

# check if 2 lists have equal size
def checkEqualListSize(list_1, list_2):
    check_valid = False
    num_of_chunks = None
    if (len(list_1) == len(list_2)):
        check_valid = True
        num_of_chunks = len(list_1)
    else:
        print("WARNING : # chunks in 1 st set is not equal to # chunks in 2nd")
    return check_valid, num_of_chunks

# calculating discriminative function
def lda(list_all_input, list_part_1, list_part_2):

    # create list with gene expression
    # list_all_input = []
    # list_part_1 = []
    # list_part_2 = []

    # for column in all_input:
    #     list_each_sample = []
    #     for element in all_input[column]:
    #         list_each_sample.append(element)
    #     list_all_input.append(list_each_sample)
    
    # for column in part_1:
    #     list_each_sample = []
    #     for element in part_1[column]:
    #         list_each_sample.append(element)
    #     list_part_1.append(list_each_sample)

    # for column in part_2:
    #     list_each_sample = []
    #     for element in part_2[column]:
    #         list_each_sample.append(element)
    #     list_part.append(list_each_sample)

    # create matrix using for calculating 
    matrix_all_input = np.matrix(list_all_input)
    matrix_part_1 = np.matrix(list_part_1)
    matrix_part_2 = np.matrix(list_part_2)
    # matrix_training_output = np.matrix(list_training_output).transpose()

    # print(matrix_training_output.transpose(1,0))
    # print("---------------------- Matrix for all input ---------------------")
    # print(matrix_all_input)
    # print("-------------------- Matrix for each feature --------------------")
    # print("Relapse ... ")
    # print(matrix_part_1)
    # print("NO relapses")
    # print(matrix_part_2)
    # print("-------------------- Matrix for output class --------------------")
    # print(matrix_training_output)
    # print("-----------------------------------------------------------------")
    # print()

    # calculate average of each feature
    avg_all_input = avgFromList(matrix_all_input)
    avg_part_1 = avgFromList(matrix_part_1)
    avg_part_2 = avgFromList(matrix_part_2)

    # calculate mean corrected data
    mean_corrected_part_1= meanCorrected(matrix_part_1, avg_all_input)
    mean_corrected_part_2 = meanCorrected(matrix_part_2, avg_all_input)

    # calculate covariance matrix
    covariance_part_1 = covariance(mean_corrected_part_1)
    covariance_part_2 = covariance(mean_corrected_part_2)

    # calculate pooled covariance matrix
    number_of_sample_part_1= matrix_part_1.shape[0]
    number_of_sample_part_2 = matrix_part_2.shape[0]
    number_of_sample = matrix_all_input.shape[0]
    pool_covariance = poolCovariance(covariance_part_1, covariance_part_2, number_of_sample, number_of_sample_part_1, \
                      number_of_sample_part_2)

    # calculate inversed matrix
    inversed_pool_covariance = inversed(pool_covariance)

    #  If we do not know the prior probability, we just assume it is equal to total sample of each group divided by the total samples
    prior_prob = findPriorProb(number_of_sample, number_of_sample_part_1, number_of_sample_part_2)

    # find output 
    f1 = findDiscriminative(matrix_all_input, avg_part_1, inversed_pool_covariance, prior_prob[0])
    f2 = findDiscriminative(matrix_all_input, avg_part_2, inversed_pool_covariance, prior_prob[1])

    actual_output = findOutput(f1, f2)

    return actual_output

def getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_id, rows_to_read_file_pathway):
    # prepare files to be used
    cols_to_read_file_to_convert = ["ID_REF", sample_id]
    # rows_to_read_file_pathway = 1329
    file_ref = pd.read_csv(file_ref_name)
    # For the last version, 'nrows' in file_to_convert has to be removed
    file_to_convert = pd.read_csv(file_to_convert_name, usecols = cols_to_read_file_to_convert, nrows = 50)
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # list all probe id
    list_probe_id = []
    for element in file_to_convert.loc[:, 'ID_REF']:
        list_probe_id.append(element)

    num_of_probe_id = len(list_probe_id)
    count_not_found = 0
    # scan probe id by its entrez id
    list_entrez_id = []
    print("Scanning in progress ...")
    print()
    for index in range(0, num_of_probe_id):
        entrez_id = file_ref.loc[file_ref['From'].isin([list_probe_id[index]])]
        # add to list if entrez_id is found
        if (entrez_id.empty):
            list_entrez_id.append("none")
            count_not_found += 1
        else:
            list_entrez_id.append(entrez_id.iloc[0][1])  
    num_of_available_data = num_of_probe_id - count_not_found

    # print(file_pathway)

    # create dictionary to collect each pathway
    pathways = {}
    for i in range(0, rows_to_read_file_pathway):
        key = i
        pathway = []
        pathway_name = file_pathway.iloc[i, 0]
        list_gene_expression = []
        for element in file_pathway.iloc[i, 2:-1]:
            if (str(element).isnumeric()):
                list_gene_name_with_expression = []
                # check if genes in each pathway have their expression value
                list_gene_same_entrez = []
                for list_entrez_id_index in range(0, len(list_entrez_id)):
                    if (element == list_entrez_id[list_entrez_id_index]):
                        # WARNING : NEED to adjust column in file_to_convert.iloc[i, 1] for use further
                        list_gene_same_entrez.append(file_to_convert.iloc[list_entrez_id_index, 1])
                # if gene expression is not found, assume it as 'zero'
                if not list_gene_same_entrez:
                    list_gene_same_entrez.append(0.0)
                
                # calculate genes expression of the same entrez id by using their average
                num_of_same_entrez = len(list_gene_same_entrez)
                avg_gene_expression = (sum(list_gene_same_entrez) / num_of_same_entrez)
                avg_gene_expression = round(avg_gene_expression, 7)
                
                # create a gene list in format [gene entrez id, gene expression]
                list_gene_name_with_expression.append(element)
                list_gene_name_with_expression.append(avg_gene_expression)
                list_gene_expression.append(list_gene_name_with_expression)
                # list_gene_expression.append(element)
                # list_gene_expression.append(avg_gene_expression)

        # combine pathway name with its gene expressions
        pathway.append(pathway_name)
        pathway.append(list_gene_expression)
        
        pathways[key] = pathway
    return pathways
    # print(pathways[1])
    # print(len(pathways))

def mean(list_input):
    size_of_list = len(list_input)
    sum_of_element = sum(list_input)
    mean_of_list = (sum_of_element / size_of_list)
    mean_of_list = round(mean_of_list, 6)
    return mean_of_list

def sd(list_input):
    # prepare required data to calculate sd
    size_of_list = len(list_input)
    mean_of_list = mean(list_input)

    # calculate sd
    sum_of_element = 0

    for element in list_input:
        sum_of_element += math.pow((element - mean_of_list), 2)
    
    sd_of_list = math.sqrt(sum_of_element / (size_of_list - 1))
    sd_of_list = round(sd_of_list, 6)
    return sd_of_list