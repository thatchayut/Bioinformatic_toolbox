import copy
import numpy as np
import pandas as pd
import math
from copy import deepcopy
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

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
    print("FIND OUTPUT ...")
    list_f1 = []
    for element_index  in range(0, len(f1)):
        for row_index in range(0, (f1[element_index].shape[0])):
            for col_index in range(0, (f1[element_index].shape[1])):
                # print(f1[element_index][row_index, col_index])
                list_f1.append(f1[element_index][row_index, col_index])
    # print("list_f1 : " + str(list_f1))

    list_f2 = []
    for element_index  in range(0, len(f2)):
        for row_index in range(0, (f2[element_index].shape[0])):
            for col_index in range(0, (f1[element_index].shape[1])):
                # print(f2[element_index][row_index, col_index])
                list_f2.append(f2[element_index][row_index, col_index])
    # print("list_f2 : " + str(list_f2))

    # find confidence score
    result = []
    for element_index in range(0, len(list_f1)):
        confidence_score = (list_f1[element_index] / list_f2[element_index])
        result.append(confidence_score)
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
        print(" WARNING : Number of chunks in the first dataset is not equal to Number of chunks in the second dataset.")
    return check_valid, num_of_chunks

# calculating discriminative function
def lda(list_all_input, list_part_1, list_part_2):

    # create matrix using for calculating 
    matrix_all_input = np.matrix(list_all_input)
    matrix_part_1 = np.matrix(list_part_1)
    matrix_part_2 = np.matrix(list_part_2)

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

    # if det is 0, the inverse matrix cannot be created. So, add a little noise to data if its det is 0.
    det_pool_covariance = np.linalg.det(pool_covariance)
    if (det_pool_covariance == 0):
        pool_covariance_row_size = pool_covariance.shape[0]
        pool_covariance_col_size = pool_covariance.shape[1]
        # noise = (0.00001 * (np.random.rand(pool_covariance_row_size, pool_covariance_col_size)))
        noise = np.random.rand(pool_covariance_row_size, pool_covariance_col_size)
        matrix_noise = np.matrix(noise)

        # add noise to pool_covariance
        pool_covariance += matrix_noise

    # calculate inversed matrix
    inversed_pool_covariance = inversed(pool_covariance)

    #  If we do not know the prior probability, we just assume it is equal to total sample of each group divided by the total samples
    prior_prob = findPriorProb(number_of_sample, number_of_sample_part_1, number_of_sample_part_2)

    # find output 
    f1 = findDiscriminative(matrix_all_input, avg_part_1, inversed_pool_covariance, prior_prob[0])
    f2 = findDiscriminative(matrix_all_input, avg_part_2, inversed_pool_covariance, prior_prob[1])

    actual_output = findOutput(f1, f2)

    return actual_output

# Old version
# def getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_id, rows_to_read_file_pathway, mean_of_data = 0, sd_of_data= 0, \
#                 max_of_data = 0, min_of_data = 0, method = "z_score"):
# function to create pathway of a given sample
def getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_id, rows_to_read_file_pathway, list_mean_sd_gene_expression_by_probe_id = None, normalize = True):

    # prepare files to be used
    cols_to_read_file_to_convert = ["ID_REF", sample_id]
    # rows_to_read_file_pathway = 1329
    file_ref = pd.read_csv(file_ref_name)
    # For the last version, 'nrows' in file_to_convert has to be removed
    file_to_convert = pd.read_csv(file_to_convert_name, usecols = cols_to_read_file_to_convert)
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # list all probe id
    list_probe_id = []
    for element in file_to_convert.loc[:, 'ID_REF']:
        list_probe_id.append(element)

    num_of_probe_id = len(list_probe_id)
    count_not_found = 0
    # scan probe id by its entrez id
    list_entrez_id = []
    print(" # Process : mapping sample " + str(sample_id) + " to pathways")
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

                        gene_probe_id = file_to_convert.iloc[list_entrez_id_index, 0]
                        gene_expression = file_to_convert.iloc[list_entrez_id_index, 1]
                        
                        # normalize gene using z-score 
                        if (normalize is True):
                            for gene_index in range(0, len(list_mean_sd_gene_expression_by_probe_id)):
                                gene_name = list_mean_sd_gene_expression_by_probe_id[gene_index][0]
                                if (gene_name == gene_probe_id):
                                    gene_mean = list_mean_sd_gene_expression_by_probe_id[gene_index][1]
                                    gene_sd = list_mean_sd_gene_expression_by_probe_id[gene_index][2]
                                    gene_expression_zscore = zscore(gene_expression, gene_mean, gene_sd)

                                    list_gene_same_entrez.append(gene_expression_zscore)
                        else:
                            list_gene_same_entrez.append(gene_expression)

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

        # combine pathway name with its gene expressions
        pathway.append(pathway_name)
        pathway.append(list_gene_expression)
        
        pathways[key] = pathway
    return pathways

def getPathwayLLR(file_ref_name, file_to_convert_name, file_pathway_name, sample_id, rows_to_read_file_pathway, list_gene_lambda_value):
    # prepare files to be used
    cols_to_read_file_to_convert = ["ID_REF", sample_id]
    # rows_to_read_file_pathway = 1329
    file_ref = pd.read_csv(file_ref_name)
    # For the last version, 'nrows' in file_to_convert has to be removed
    file_to_convert = pd.read_csv(file_to_convert_name, usecols = cols_to_read_file_to_convert)
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # list all probe id
    list_probe_id = []
    for element in file_to_convert.loc[:, 'ID_REF']:
        list_probe_id.append(element)

    num_of_probe_id = len(list_probe_id)
    count_not_found = 0

    # scan probe id by its entrez id
    # find matching between probe id and entrez id
    list_entrez_id = []
    print(" # Process : mapping sample " + str(sample_id) + " to pathways")
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

                        gene_probe_id = file_to_convert.iloc[list_entrez_id_index, 0]
                        gene_expression = file_to_convert.iloc[list_entrez_id_index, 1]
                        
                        for gene_index in range(0, len(list_gene_lambda_value)):
                            gene_name = list_gene_lambda_value[gene_index][0]

                            if (gene_name == gene_probe_id):
                                lambda_value = list_gene_lambda_value[gene_index][1]
                                list_gene_same_entrez.append(lambda_value)

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

        # combine pathway name with its gene expressions
        pathway.append(pathway_name)
        pathway.append(list_gene_expression)
        
        pathways[key] = pathway
    return pathways


def zscore(value, mean_value, sd_value):
    result = ((value - mean_value) / sd_value)
    return result

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

def sfs(list_pathway_name, list_desired_output, samples_relapse, samples_no_relapse, samples_test):
    # print("list_pathway_name size : " + str(len(list_pathway_name)))
    # print("list_desired_output : " + str(len(list_desired_output)))
    # print("samples_relapse size : " + str(len(samples_relapse)))
    # print("samples_no_relapse size : " + str(len(samples_no_relapse)))
    # print("samples_test size : " + str(len(samples_test)))
    # print("\n")

    check_finish = False
    check_improve_auc = True

    max_auc_score_over_all_features = 0

    feature_set_final = []
    
    # original 
    # list_pathway_selected = []

    # version 1
    # fix top pathway as the first index
    list_pathway_selected = [0]

    # version 2 
    # select 5% of top pathways
    # percent = 5
    # num_of_top_n_percent_pathway = math.ceil((len(list_pathway_name) * percent) / 100)
    # list_pathway_selected = []
    # for index in range(0, num_of_top_n_percent_pathway):
    #     list_pathway_selected.append(index)

    num_of_pathways = len(list_pathway_name)

    while (check_finish == False):

        if (check_improve_auc == True):
            max_auc_in_consider = 0
            list_pathway = []
            for pathway_index in range(0, num_of_pathways):
                list_pathway_to_consider = deepcopy(list_pathway_selected)
                
                # original 
                if (pathway_index not in list_pathway_selected):
                    list_pathway_to_consider.extend([pathway_index])
                    
                    input_relapse_to_test = []
                    for sample_index in range(0, len(samples_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_relapse in range(0, len(samples_relapse[sample_index])):
                            if (pathway_index_relapse in list_pathway_to_consider):
                                pathway_activity_to_test = samples_relapse[sample_index][pathway_index_relapse]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_relapse_to_test.append(list_pathway_each_sample_to_test)

                    input_no_relapse_to_test = []
                    for sample_index in range(0, len(samples_no_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_no_relapse in range(0, len(samples_no_relapse[sample_index])):
                            if (pathway_index_no_relapse in list_pathway_to_consider):
                                pathway_activity_to_test = samples_no_relapse[sample_index][pathway_index_no_relapse]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_no_relapse_to_test.append(list_pathway_each_sample_to_test)
                    
                    input_test_to_test = []
                    for sample_index in range(0, len(samples_test)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_test in range(0, len(samples_test[sample_index])):
                            if (pathway_index_test in list_pathway_to_consider):
                                pathway_activity_to_test = samples_test[sample_index][pathway_index_test]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_test_to_test.append(list_pathway_each_sample_to_test)

                    list_actual_output = lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)
                    auc_score = roc_auc_score(list_desired_output, list_actual_output)

                    # print("input_relapse_to_test size : " + str(len(input_relapse_to_test)))
                    # print("input_no_relapse_to_test size : " + str(len(input_no_relapse_to_test)))
                    # print("input_test_to_test : " + str(len(input_test_to_test)))
                    
                    list_actual_output = lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)
                    print(" Number of actual outputs : " +str(len(list_actual_output)))
                    print(" Number of desired outputs : " + str(len(list_desired_output)))
                    auc_score = roc_auc_score(list_desired_output, list_actual_output)
                    
                    print(" List of actual outputs : " + str(list_actual_output))
                    print(" List of desired outputs : " + str(list_desired_output))
                    print(" AUC score : " + str(auc_score))

                    if (auc_score >= max_auc_in_consider):
                        max_auc_in_consider = auc_score
                        list_pathway = deepcopy(list_pathway_to_consider)
                    print(" List of pathways to be considered : " + str(list_pathway_to_consider))
                    print(" Max AUC score of the list of pathways to be considered : " + str(max_auc_in_consider))

            if (max_auc_in_consider >= max_auc_score_over_all_features):
                max_auc_score_over_all_features = max_auc_in_consider
                list_pathway_selected = deepcopy(list_pathway)
            else:
                check_improve_auc = False
        else:
            check_finish = True

        list_pathway_selected_name = []
        for i in range(0, len(list_pathway_selected)):
            pathway = list_pathway_selected[i]
            list_pathway_selected_name.append(list_pathway_name[pathway])

        print(" Feature set : " + str(list_pathway_selected_name))
        print(" AUC score from this feature set : " + str(max_auc_score_over_all_features))

        feature_set_final = deepcopy(list_pathway_selected_name)

    return feature_set_final, max_auc_score_over_all_features

def sfsAdvance(list_pathway_name, list_desired_output, samples_relapse_input, samples_no_relapse_input, samples_test_input, \
    top_rank = None, top_rank_percent = None):
    # print("list_pathway_name size : " + str(len(list_pathway_name)))
    # print("list_desired_output : " + str(len(list_desired_output)))
    # print("samples_relapse_input size : " + str(len(samples_relapse_input)))
    # print("samples_no_relapse_input size : " + str(len(samples_no_relapse_input)))
    # print("samples_test_input size : " + str(len(samples_test_input)))
    # print("\n")

    # sort pathways in each sample based on list_pathway_name
    # for class 'relapse'
    samples_relapse = []
    for sample_index in range(0, len(samples_relapse_input)):
        sample = []
        for top_pathway_name in list_pathway_name:
            for pathway_index in range(0, len(samples_relapse_input[sample_index])):
                pathway_name = samples_relapse_input[sample_index][pathway_index][0]
                if (pathway_name == top_pathway_name):
                    pathway_activity = samples_relapse_input[sample_index][pathway_index][1]

                    sample.append(pathway_activity)
        samples_relapse.append(sample)

    # for class 'non-relapse'
    samples_no_relapse = []
    for sample_index in range(0, len(samples_no_relapse_input)):
        sample = []
        for top_pathway_name in list_pathway_name:
            for pathway_index in range(0, len(samples_no_relapse_input[sample_index])):
                pathway_name = samples_no_relapse_input[sample_index][pathway_index][0]
                if (pathway_name == top_pathway_name):
                    pathway_activity = samples_no_relapse_input[sample_index][pathway_index][1]

                    sample.append(pathway_activity)
        samples_no_relapse.append(sample)
    
    # for testing data
    samples_test= []
    for sample_index in range(0, len(samples_test_input)):
        sample = []
        for top_pathway_name in list_pathway_name:
            for pathway_index in range(0, len(samples_test_input[sample_index])):
                pathway_name = samples_test_input[sample_index][pathway_index][0]
                if (pathway_name == top_pathway_name):
                    pathway_activity = samples_test_input[sample_index][pathway_index][1]

                    sample.append(pathway_activity)
        samples_test.append(sample)


    # preparing for feature selection
    check_finish = False
    check_improve_auc = True

    max_auc_score_over_all_features = 0

    feature_set_final = []
    
    # original 
    # list_pathway_selected = []

    # # fix top pathway as the first index
    # list_pathway_selected = [0]

    # version 2 
    # select 5% of top pathways
    # percent = 5
    # num_of_top_n_percent_pathway = math.ceil((len(list_pathway_name) * percent) / 100)
    # list_pathway_selected = []
    # for index in range(0, num_of_top_n_percent_pathway):
    #     list_pathway_selected.append(index)
    

    # get the list to be used as an initial set
    list_pathway_selected = []

    # if use top n pathway
    if (top_rank is not None):
        num_of_top_pathway = top_rank
        for index in (0, num_of_top_pathway):
            list_pathway_selected.append(index)
    # if use top n percent pathway
    elif (top_rank_percent is not None):
        num_of_top_percent_pathway = math.ceil((len(list_pathway_name) * percent) / 100)
        for index in range(0, num_of_top_percent_pathway):
            list_pathway_selected.append(index)
    # if use only the highest rank
    else:
        list_pathway_selected = [0]



    num_of_pathways = len(list_pathway_name)

    while (check_finish == False):

        if (check_improve_auc == True):
            max_auc_in_consider = 0
            list_pathway = []
            for pathway_index in range(0, num_of_pathways):
                list_pathway_to_consider = deepcopy(list_pathway_selected)
                
                # original 
                if (pathway_index not in list_pathway_selected):
                # last_index_in_list_pathway_to_consider = len(list_pathway_selected) - 1
                # if (pathway_index not in list_pathway_selected) and (pathway_index > list_pathway_selected[last_index_in_list_pathway_to_consider]):
                    list_pathway_to_consider.extend([pathway_index])

                    input_relapse_to_test = []
                    for sample_index in range(0, len(samples_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_relapse in range(0, len(samples_relapse[sample_index])):
                            if (pathway_index_relapse in list_pathway_to_consider):
                                pathway_activity_to_test = samples_relapse[sample_index][pathway_index_relapse]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_relapse_to_test.append(list_pathway_each_sample_to_test)

                    input_no_relapse_to_test = []
                    for sample_index in range(0, len(samples_no_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_no_relapse in range(0, len(samples_no_relapse[sample_index])):
                            if (pathway_index_no_relapse in list_pathway_to_consider):
                                pathway_activity_to_test = samples_no_relapse[sample_index][pathway_index_no_relapse]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_no_relapse_to_test.append(list_pathway_each_sample_to_test)

                    input_test_to_test = []
                    for sample_index in range(0, len(samples_test)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index_test in range(0, len(samples_test[sample_index])):
                            if (pathway_index_test in list_pathway_to_consider):
                                pathway_activity_to_test = samples_test[sample_index][pathway_index_test]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_test_to_test.append(list_pathway_each_sample_to_test)

                    list_actual_output = lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)
                    auc_score = roc_auc_score(list_desired_output, list_actual_output)
                    
                    list_actual_output = lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)

                    auc_score = roc_auc_score(list_desired_output, list_actual_output)
                    
                    print(" List of actual outputs : " + str(list_actual_output))
                    print(" List of desired outputs : " + str(list_desired_output))
                    print(" AUC score : " + str(auc_score))

                    if (auc_score >= max_auc_in_consider):
                        max_auc_in_consider = auc_score
                        list_pathway = deepcopy(list_pathway_to_consider)

            if (max_auc_in_consider >= max_auc_score_over_all_features):
                max_auc_score_over_all_features = max_auc_in_consider
                list_pathway_selected = deepcopy(list_pathway)
            else:
                check_improve_auc = False
        else:
            check_finish = True

        list_pathway_selected_name = []
        for i in range(0, len(list_pathway_selected)):
            pathway = list_pathway_selected[i]
            list_pathway_selected_name.append(list_pathway_name[pathway])

        print(" Feature set : " + str(list_pathway_selected_name))
        print(" AUC score from this feature set : " + str(max_auc_score_over_all_features))

        feature_set_final = deepcopy(list_pathway_selected_name)

    return feature_set_final, max_auc_score_over_all_features





