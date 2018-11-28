#!/usr/bin/python
import pandas as pd
import math
import calculate
import random
from scipy import stats
import numpy as np
from sklearn.metrics import roc_auc_score
from copy import deepcopy
import time

def main():
    # record start time
    start_time = time.time()

    # prepare data
    # row_to_read = 22283
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    
    # version 1: consider only relapse and non-relapse within 5 years
    # file_training_output = pd.read_csv("mapping_sample_to_class_relapse.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])
    
    # version 2: consider non-relapse and relapse (not in specific period of time)
    file_training_output_relapse = pd.read_csv("mapping_sample_to_class_relapse.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])
    file_training_output_no_relapse = pd.read_csv("mapping_sample_to_class_no_relapse.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])
    file_training_output = pd.read_csv("mapping_sample_to_class_no_relapse.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # this will be used in calculating lda
    training_input = file_training_input
    training_output = file_training_output

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    # print(list_gene_name)

    # separate data into 2 classes

    # version 1: consider only relapse and non-relapse within 5 years
    # sample_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['1'])]
    # sample_no_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['0'])]

    # version 2: consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_training_output_relapse.loc[file_training_output_relapse['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_training_output_no_relapse.loc[file_training_output_no_relapse['relapse (1=True)'].isin(['0'])]
    # print(sample_no_relapse)
    
    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)
    # print(list_sample_relapse)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)

    # shuffle data to make each chunk does not depend on sample order
    random.shuffle(list_sample_relapse)
    print("list_sample_relapse SIZE = " + str(len(list_sample_relapse)))
    random.shuffle(list_sample_no_relapse)
    print("list_sample_no_relapse SIZE = " + str(len(list_sample_no_relapse)))

    # get number of folds
    while True:
        num_of_folds = input("Number of folds: ")
        if (num_of_folds.isnumeric() == False):
            print("WARNING : Invalid input must be numeric")
        elif(int(num_of_folds) > len(list_sample_relapse)):
            print("WARNING : Number of folds exceeds the size of the 1st dataset")
        elif(int(num_of_folds) > len(list_sample_no_relapse)):
            print("WARNING : Number of folds exceeds the size of the 2nd dataset")
        elif(int(num_of_folds) <= 1):
            print("WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)

    # get number of ranked gene to be shown
    while True:
        number_of_ranked_gene = input("Number of ranked feature: ")
        if ((number_of_ranked_gene.isnumeric() == False) or (int(number_of_ranked_gene) > row_to_read) or (int(number_of_ranked_gene) <= 0)):
            print("Invalid input...")
        else:
            break
    
    # get epoch
    while True:
        epoch = input("Epoch : ")
        if ((epoch.isnumeric() == False) or (int(epoch) <= 0)):
            print("Invalid input...")
        else:
            break

    # get output file's name
    file_name = input("Name of Output File : ")

    # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    for epoch_count in range(0, int(epoch)):
        start_epoch_time = time.time()
        result_file.write("#################################### Epoch : " + str(epoch_count + 1) + " ####################################\n")
        print("#################################### Epoch : " + str(epoch_count + 1) + " ####################################")
        # split data into k parts
        chunk_relapse_size = math.ceil(len(list_sample_relapse) / num_of_folds)
        chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse) / num_of_folds)
        # print(chunk_size)    

        chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_relapse_size))
        print("# chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

        chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_no_relapse_size))
        print("# chunks in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # list to collect maximun AUC in each fold
        list_max_auc = []

        # do only if number of chunks of both datasets are equal
        if (check_valid == True):
            for first_layer_test_index in range(0, num_of_chunks):
                feature_set = []
                feature_set_name = []
                top_n_genes_name_for_eval = []
                # keep testing data from each class
                first_layer_test_relapse = chunk_list_relapse[first_layer_test_index]
                first_layer_test_no_relapse = chunk_list_no_relapse[first_layer_test_index]
                print("\n------------------------------------------ K : " + str(first_layer_test_index + 1) + " of Epoch " + str(epoch_count + 1) + " --------------------------------")
                print("test relapse =" + str(first_layer_test_relapse))
                print("test no relapse = " + str(first_layer_test_no_relapse))
                print()
                # find training data
                # first layer 
                first_layer_train_relapse = []
                for first_layer_train_index in range(0, num_of_chunks):
                    if (chunk_list_relapse[first_layer_train_index] is not first_layer_test_relapse):
                        first_layer_train_relapse.append(chunk_list_relapse[first_layer_train_index])
                print("1st layer train relapse size = " + str(len(first_layer_train_relapse)))
                print("1st layer train relapse = " + str(first_layer_train_relapse))

                first_layer_train_no_relapse = []
                for first_layer_train_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse[first_layer_train_index] is not first_layer_test_no_relapse):
                        first_layer_train_no_relapse.append(chunk_list_no_relapse[first_layer_train_index])
                print("1st layer train no relapse size = " + str(len(first_layer_train_no_relapse)))
                print("1st layer train no relapse = " + str(first_layer_train_no_relapse))

                # merge all element in each list to be used in second layer
                print("\n##### merge remaining element in each training list to be used in the next step #####")
                second_list_sample_relapse = []
                for i in range(0, len(first_layer_train_relapse)):
                    second_list_sample_relapse.extend(first_layer_train_relapse[i])
                print("size of list sample relapse = " + str(len(second_list_sample_relapse)))
                print("list sample relapse for next step = " + str(second_list_sample_relapse))

                second_list_sample_no_relapse = []
                for i in range(0, len(first_layer_train_no_relapse)):
                    second_list_sample_no_relapse.extend(first_layer_train_no_relapse[i])
                print("size of list sample no relapse  = " + str(len(second_list_sample_no_relapse)))
                print("list sample no relapse for next step = " + str(second_list_sample_no_relapse))

                # splitting lists to use them as testing and training set
                # given that we use 4-fold cross validation in this level
                print("\n#### given that we use 4-fold cross validation in this level ####")
                second_num_of_fold = 4
                second_chunk_relapse_size = math.ceil(len(second_list_sample_relapse) / second_num_of_fold)
                second_chunk_no_relapse_size = math.ceil(len(second_list_sample_no_relapse) / second_num_of_fold)

                second_chunk_list_relapse = list(calculate.chunks(second_list_sample_relapse, second_chunk_relapse_size))
                print("# chunks in second_chunk_list_relapse = " + str(len(second_chunk_list_relapse)))
                second_chunk_list_no_relapse = list(calculate.chunks(second_list_sample_no_relapse, second_chunk_no_relapse_size))
                print("# chunks in second_chunk_list_no_relapse = " + str(len(second_chunk_list_no_relapse)))

                second_check_valid, second_num_of_chunks = calculate.checkEqualListSize(second_chunk_list_relapse, second_chunk_list_no_relapse)

                # do only if number of chunks of both datasets are equal
                if (second_check_valid == True):
                    # for second_layer_test_index in range(0, second_num_of_chunks):
                    second_layer_test_index = random.randint(0, second_num_of_chunks - 1)
                    # keep testing data from eacch class
                    second_layer_test_relapse =  second_chunk_list_relapse[second_layer_test_index]
                    second_layer_test_no_relapse = second_chunk_list_no_relapse[second_layer_test_index]
                    # print("\n------------------------------------------ L : " + str(second_layer_test_index + 1) + " --------------------------------")
                    print("second test relapse =" + str(second_layer_test_relapse))
                    print("second test no relapse = " + str(second_layer_test_no_relapse))
                    print()
                
                    # separate training dataset from testing dataset to use in t-test ranking
                    second_layer_train_relapse = []
                    for second_layer_train_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_relapse[second_layer_train_index] is not second_layer_test_relapse):
                            second_layer_train_relapse.append(second_chunk_list_relapse[second_layer_train_index])
                    print("2nd layer train relapse size = " + str(len(second_layer_train_relapse)))
                    print("2nd layer train relapse = " + str(second_layer_train_relapse))
                    
                    second_layer_train_no_relapse = []
                    for second_layer_train_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_no_relapse[second_layer_train_index] is not second_layer_test_no_relapse):
                            second_layer_train_no_relapse.append(second_chunk_list_no_relapse[second_layer_train_index])
                    print("2nd layer train no relapse size = " + str(len(second_layer_train_no_relapse)))
                    print("2nd layer train no relapse = " + str(second_layer_train_no_relapse))

                    # prepare dataset for t-test
                    # merge all samples in the same class
                    print("\n#### merge all samples in the same class to be used in t-test ####")
                    ttest_list_sample_relapse = []
                    for i in range(0, len(second_layer_train_relapse)):
                        ttest_list_sample_relapse.extend(second_layer_train_relapse[i])
                    print("size of ttest list sample relapse = " + str(len(ttest_list_sample_relapse)))
                    print("ttest list sample relapse = " + str(ttest_list_sample_relapse))

                    ttest_list_sample_no_relapse = []
                    for i in range(0, len(second_layer_train_no_relapse)):
                        ttest_list_sample_no_relapse.extend(second_layer_train_no_relapse[i])
                    print("size of ttest list sample no relapse = " + str(len(ttest_list_sample_no_relapse)))
                    print("ttest list sample no relapse = " + str(ttest_list_sample_no_relapse))

                    # get gene expression for each gene from samples with relapse
                    list_gene_exp_relapse = []
                    for i in range(0, row_to_read):
                        gene_exp_relapse = []
                        for column in  file_training_input.loc[i, ttest_list_sample_relapse]:
                            gene_exp_relapse.append(column)
                        list_gene_exp_relapse.append(gene_exp_relapse)
                    # print(list_gene_exp_relapse)

                    # get gene expression for each gene from samples with no relapse
                    list_gene_exp_no_relapse = []
                    for i in range(0, row_to_read):
                        gene_exp_no_relapse = []
                        for column in  file_training_input.loc[i, ttest_list_sample_no_relapse]:
                            gene_exp_no_relapse.append(column)
                        list_gene_exp_no_relapse.append(gene_exp_no_relapse)
                    # print(list_gene_exp_no_relapse)

                    # conducting t-test
                    ttest_result = []
                    for i in range(0, row_to_read):      
                        score = []
                        # print(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False))
                        # get absolute magnitude of t-test value
                        abs_ttest_value = math.fabs(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[0])
                        p_value = stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[1]
                        # add element with this format (gene_order_id, ttest_value)
                        score.append(i)
                        score.append(abs_ttest_value)
                        ttest_result.append(score)
                    # ranking elements using their t-test value in descending order
                    ttest_result.sort(key=lambda x: x[1], reverse=True)
                    # print(ttest_result)

                    # create list of ranked gene
                    ranked_gene = []
                    for i in range(0, len(ttest_result)):
                        gene_order_id = ttest_result[i][0]
                        # print(gene_order_id)
                        # print(list_gene_name[gene_order_id][1])
                        ranked_gene.append(list_gene_name[gene_order_id][1])
                    # print(ranked_gene)

                    # show top ranked feature
                    top_n_genes_name = []
                    print("#### t-test ranking ####")
                    for i in range(0, int(number_of_ranked_gene)):
                        top_n_genes_name.append(ranked_gene[i])
                        print(ranked_gene[i] + " => " + "t-test value : " + str(ttest_result[i][1]))

                    # rank gene id of each sample in training data
                    # print("\n#### sorting gene order by t-test ranking for each class ####")
                    # for class 'relapse'
                    # print("#### class 'Relapse' ####")
                    col_to_read_relapse = ["ID_REF"]
                    col_to_read_relapse.extend(second_layer_train_relapse[0])
                    # print(col_to_read_relapse)
                    file_training_input_relapse = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_relapse)
                    # print(file_training_input_relapse)
                    top_n_genes_relapse = file_training_input_relapse.loc[file_training_input_relapse['ID_REF'].isin(top_n_genes_name)]
                    # print(top_n_genes_relapse)
                    top_n_genes_relapse['gene_id'] = top_n_genes_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))
                    # print(top_n_genes_relapse)
                    top_n_genes_relapse_sorted  = top_n_genes_relapse.sort_values(by = ['gene_id'])
                    top_n_genes_relapse_sorted.drop(columns = 'gene_id', inplace = True)

                    top_n_genes_relapse_sorted_train = top_n_genes_relapse_sorted
                    top_n_genes_relapse_sorted_train.drop(columns = 'ID_REF', inplace = True)
                    # print(top_n_genes_relapse_sorted_train)

                    # for class 'no relapse'
                    # print("#### class 'no Relapse' ####")
                    col_to_read_no_relapse = ["ID_REF"]
                    col_to_read_no_relapse.extend(second_layer_train_no_relapse[0])
                    # print(col_to_read_no_relapse)
                    file_training_input_no_relapse = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_no_relapse)
                    # print(file_training_input_no_relapse)
                    top_n_genes_no_relapse = file_training_input_no_relapse.loc[file_training_input_no_relapse['ID_REF'].isin(top_n_genes_name)]
                    # print(top_n_genes_no_relapse)
                    top_n_genes_no_relapse['gene_id'] = top_n_genes_no_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))
                    # print(top_n_genes_no_relapse)
                    top_n_genes_no_relapse_sorted  = top_n_genes_no_relapse.sort_values(by = ['gene_id'])
                    top_n_genes_no_relapse_sorted.drop(columns = 'gene_id', inplace = True)

                    top_n_genes_no_relapse_sorted_train = top_n_genes_no_relapse_sorted
                    top_n_genes_no_relapse_sorted_train.drop(columns = 'ID_REF', inplace = True)
                    # print(top_n_genes_no_relapse_sorted_train)
                    
                    # Preparing testing data for feature selection
                    # print("#### Testing data relapse & no-relapse for feature selection ####")
                    second_layer_test_all = []
                    second_layer_test_all.extend(second_layer_test_relapse)
                    second_layer_test_all.extend(second_layer_test_no_relapse)   
                    # print(second_layer_test_all)  
                    # output for testing data
                    # second_layer_test_output = trainipseng_output.loc[training_output['GEO asscession number'].isin(second_layer_test_all)]
                    # print(second_layer_test_output)
                    # sort gene order of testing data
                    col_to_read_second_layer_test_gene = ["ID_REF"]
                    col_to_read_second_layer_test_gene.extend(second_layer_test_all)
                    second_layer_test_gene = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_second_layer_test_gene)
                    second_layer_top_n_test = second_layer_test_gene.loc[second_layer_test_gene['ID_REF'].isin(top_n_genes_name)]
                    second_layer_top_n_test['gene_id'] = second_layer_top_n_test['ID_REF'].apply(lambda name: top_n_genes_name.index(name))   
                    second_layer_top_n_test_sorted = second_layer_top_n_test.sort_values(by = ['gene_id'])
                    second_layer_top_n_test_sorted.drop(columns = 'gene_id', inplace = True)

                    top_n_test_sorted = second_layer_top_n_test_sorted
                    top_n_test_sorted.drop(columns = 'ID_REF', inplace = True)
                    # print(top_n_test_sorted)

                    # use top-rank feature as the first feature in lda classifier
                    # prepare list for input 
                    # list of all input data (testing data)
                    list_second_layer_top_n_test_sorted = []
                    for column in range(0, len(top_n_test_sorted)):
                        list_each_sample = []
                        for element in top_n_test_sorted.iloc[column]:
                            list_each_sample.append(element)
                            # list_each_sample = list(np.transpose(list_each_sample))
                            # print(list_each_sample)
                        list_second_layer_top_n_test_sorted.append(list_each_sample)
                    list_second_layer_top_n_test_sorted = list(np.transpose(list_second_layer_top_n_test_sorted))
                    # print(list_second_layer_top_n_test_sorted)

                    # output for testing data
                    second_layer_test_output = training_output.loc[training_output['GEO asscession number'].isin(second_layer_test_all)]
                    # sorting data according to its order in testing data
                    list_sample_to_read = list(second_layer_top_n_test_sorted.columns.values)
                    # print(list_sample_to_read)
                    second_layer_test_output['sample_id'] = second_layer_test_output['GEO asscession number'].apply(lambda name: list_sample_to_read.index(name))
                    second_layer_test_output = second_layer_test_output.sort_values(by = ['sample_id'])
                    second_layer_test_output.drop(columns = 'sample_id', inplace = True)
                    # create list of output
                    list_desired_output = []
                    for element in second_layer_test_output.loc[:, 'relapse (1=True)']:
                        list_desired_output.append(element)

                    # list of gene expression and sample of class 'relapse'
                    list_top_n_gene_relapse_sorted = []
                    for column in range(0, len(top_n_genes_relapse_sorted_train)):
                        list_each_sample = []
                        for element in top_n_genes_relapse_sorted_train.iloc[column]:
                            list_each_sample.append(element)
                        list_top_n_gene_relapse_sorted.append(list_each_sample)
                    list_top_n_gene_relapse_sorted = list(np.transpose(list_top_n_gene_relapse_sorted))
                    # print(list_top_n_gene_relapse_sorted)
                    
                    # list of gene expression and sample of class 'no relapse'
                    list_top_n_gene_no_relapse_sorted = []
                    for column in range(0, len(top_n_genes_no_relapse_sorted_train)):
                        list_each_sample = []
                        for element in top_n_genes_no_relapse_sorted_train.iloc[column]:
                            list_each_sample.append(element)
                        list_top_n_gene_no_relapse_sorted.append(list_each_sample)
                    list_top_n_gene_no_relapse_sorted = list(np.transpose(list_top_n_gene_no_relapse_sorted))
                    # print(list_top_n_gene_no_relapse_sorted)

                    # find set of genes to be used as a feature
                    check_finish = False
                    count_iteration = 1
                    gene_order = [0]
                    list_auc = []
                    while (check_finish == False):
                        if (count_iteration >= int(number_of_ranked_gene)):
                            check_finish = True
                        else:
                            max_auc_score = 0
                            gene_index_in_list = None
                            for i in range(0, int(number_of_ranked_gene)):
                                gene_order_test = deepcopy(gene_order)
                                gene_order_test.extend([i])
                                # select gene to be used in lda                       
                                input_relapse = []
                                for sample_index in range(0, len(list_top_n_gene_relapse_sorted)):
                                    list_each_sample = []
                                    for element_id in range(0, len(list_top_n_gene_relapse_sorted[sample_index])):
                                        if (element_id in gene_order_test):
                                            list_each_sample.append(list_top_n_gene_relapse_sorted[sample_index][element_id])
                                    input_relapse.append(list_each_sample)
                                # print(input_relapse)

                                input_no_relapse = []
                                for sample_index in range(0, len(list_top_n_gene_no_relapse_sorted)):
                                    list_each_sample = []
                                    for element_id in range(0, len(list_top_n_gene_no_relapse_sorted[sample_index])):
                                        if (element_id in gene_order_test):
                                            list_each_sample.append(list_top_n_gene_no_relapse_sorted[sample_index][element_id])
                                    input_no_relapse.append(list_each_sample)
                                # print(input_no_relapse)

                                input_testing_data = []
                                for sample_index in range(0, len(list_second_layer_top_n_test_sorted)):
                                    list_each_sample = []
                                    for element_id in range(0, len(list_second_layer_top_n_test_sorted[sample_index])):
                                        if (element_id in gene_order_test):
                                            list_each_sample.append(list_second_layer_top_n_test_sorted[sample_index][element_id])
                                    input_testing_data.append(list_each_sample)
                                # print(input_testing_data)
                                list_actual_output = calculate.lda(input_testing_data, input_relapse, input_no_relapse)
                                # print("actual output : " + str(list_actual_output))
                                # print(len(list_actual_output))
                                # print("desired output : " + str(second_layer_test_output))
                                # print("desired output : " + str(list_desired_output))
                                # print(len(list_desired_output))

                                # calculate AUC score
                                auc_score = roc_auc_score(list_desired_output, list_actual_output)
                                # print("auc_score = " + str(auc_score))

                                if (auc_score > max_auc_score):
                                    max_auc_score = auc_score
                                    gene_index_in_list = i
                                    print(max_auc_score)
                                    if max_auc_score not in list_auc:
                                        list_auc.append(max_auc_score)
                            # do not add gene that already exists in a feature
                            if (gene_index_in_list not in gene_order):
                                gene_order.extend([gene_index_in_list])                       
                            count_iteration += 1  

                    list_max_auc.append(max(list_auc))        
                    gene_order.sort()       
                    # print("gene_order used as feature : " + str(gene_order))

                    # get gene_name
                    gene_order_name = []
                    for element in gene_order:
                        gene_order_name.append(top_n_genes_name[element])
                    # print("feature : " + str(gene_order_name))

                    # copy required data to be used in evaluation
                    top_n_genes_name_for_eval = deepcopy(top_n_genes_name)
                    feature_set  = deepcopy(gene_order)
                    feature_set_name = deepcopy(gene_order_name)

                # preparing data for evaluation and creating classifier
                # for class 'relapse'
                print("#### class 'Relapse' for creating classifier ####")
                col_to_read_relapse_for_eval = ["ID_REF"]
                # col_to_read_relapse_for_eval.extend(first_layer_test_relapse)
                col_to_read_relapse_for_eval.extend(first_layer_train_relapse[0])
                # print(col_to_read_relapse)
                file_training_input_relapse_for_eval = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_relapse_for_eval)
                # print(file_training_input_relapse)
                top_n_genes_relapse_for_eval = file_training_input_relapse.loc[file_training_input_relapse['ID_REF'].isin(feature_set_name)]
                # print(top_n_genes_relapse)
                top_n_genes_relapse_for_eval['gene_id'] = top_n_genes_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))
                # print(top_n_genes_relapse)
                top_n_genes_relapse_sorted_for_eval  = top_n_genes_relapse_for_eval.sort_values(by = ['gene_id'])
                top_n_genes_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)
                # top_n_genes_relapse_sorted_train = top_n_genes_relapse_sorted
                top_n_genes_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)
                print(top_n_genes_relapse_sorted_for_eval)

                # for class 'no relapse'
                print("#### class 'no Relapse' for creating classifier ####")
                col_to_read_no_relapse_for_eval = ["ID_REF"]
                # col_to_read_no_relapse_for_eval.extend(first_layer_test_no_relapse)
                col_to_read_no_relapse_for_eval.extend(first_layer_train_no_relapse[0])
                # print(col_to_read_no_relapse)
                file_training_input_no_relapse_for_eval = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_no_relapse_for_eval)
                # print(file_training_input_no_relapse)
                top_n_genes_no_relapse_for_eval = file_training_input_no_relapse_for_eval.loc[file_training_input_no_relapse_for_eval['ID_REF'].isin(feature_set_name)]
                # print(top_n_genes_no_relapse)
                top_n_genes_no_relapse_for_eval['gene_id'] = top_n_genes_no_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))
                # print(top_n_genes_no_relapse)
                top_n_genes_no_relapse_sorted_for_eval  = top_n_genes_no_relapse_for_eval.sort_values(by = ['gene_id'])
                top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)
                # top_n_genes_no_relapse_sorted_train = top_n_genes_no_relapse_sorted
                top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)
                print(top_n_genes_no_relapse_sorted_for_eval)            

                print("#### Testing data relapse & no-relapse for evaluation ####")
                first_layer_test_all = []
                first_layer_test_all.extend(first_layer_test_relapse)
                first_layer_test_all.extend(first_layer_test_no_relapse)  
                print(first_layer_test_all)

                col_to_read_first_layer_test_gene = ["ID_REF"]
                col_to_read_first_layer_test_gene.extend(first_layer_test_all)
                first_layer_test_gene = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_first_layer_test_gene)
                first_layer_top_n_test = first_layer_test_gene.loc[first_layer_test_gene['ID_REF'].isin(feature_set_name)]
                first_layer_top_n_test['gene_id'] = first_layer_top_n_test['ID_REF'].apply(lambda name: feature_set_name.index(name))   
                first_layer_top_n_test_sorted = first_layer_top_n_test.sort_values(by = ['gene_id'])
                first_layer_top_n_test_sorted.drop(columns = 'gene_id', inplace = True)

                top_n_test_sorted_for_eval = first_layer_top_n_test_sorted
                top_n_test_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)
                # print(top_n_test_sorted_for_eval)

                # prepare list for input 
                # list of all input data (testing data)
                list_first_layer_top_n_test_sorted = []
                for column in range(0, len(top_n_test_sorted_for_eval)):
                    list_each_sample = []
                    for element in top_n_test_sorted_for_eval.iloc[column]:
                        list_each_sample.append(element)
                        # list_each_sample = list(np.transpose(list_each_sample))
                        # print(list_each_sample)
                    list_first_layer_top_n_test_sorted.append(list_each_sample)
                list_first_layer_top_n_test_sorted = list(np.transpose(list_first_layer_top_n_test_sorted))
                # print(top_n_test_sorted_for_eval)

                # output for testing data
                first_layer_test_output = training_output.loc[training_output['GEO asscession number'].isin(first_layer_test_all)]
                # sorting data according to its order in testing data
                list_sample_to_read_for_eval = list(first_layer_top_n_test_sorted.columns.values)
                # print(list_sample_to_read)
                first_layer_test_output['sample_id'] = first_layer_test_output['GEO asscession number'].apply(lambda name: list_sample_to_read_for_eval.index(name))
                first_layer_test_output = first_layer_test_output.sort_values(by = ['sample_id'])
                first_layer_test_output.drop(columns = 'sample_id', inplace = True)
                # create list of output
                list_desired_output_for_eval = []
                for element in first_layer_test_output.loc[:, 'relapse (1=True)']:
                    list_desired_output_for_eval.append(element)

                # list of gene expression and sample of class 'relapse' for evaluation
                list_top_n_gene_relapse_sorted_for_eval = []
                for column in range(0, len(top_n_genes_relapse_sorted_for_eval)):
                    list_each_sample = []
                    for element in top_n_genes_relapse_sorted_for_eval.iloc[column]:
                        list_each_sample.append(element)
                    list_top_n_gene_relapse_sorted_for_eval.append(list_each_sample)
                list_top_n_gene_relapse_sorted_for_eval = list(np.transpose(list_top_n_gene_relapse_sorted_for_eval))
                # print(list_top_n_gene_relapse_sorted_for_eval)

                # list of gene expression and sample of class 'no relapse' for evaluation
                list_top_n_gene_no_relapse_sorted_for_eval = []
                for column in range(0, len(top_n_genes_no_relapse_sorted_for_eval)):
                    list_each_sample = []
                    for element in top_n_genes_no_relapse_sorted_for_eval.iloc[column]:
                        list_each_sample.append(element)
                    list_top_n_gene_no_relapse_sorted_for_eval.append(list_each_sample)
                list_top_n_gene_no_relapse_sorted_for_eval = list(np.transpose(list_top_n_gene_no_relapse_sorted_for_eval))
                # print(list_top_n_gene_no_relapse_sorted)            

                # calculate lda to get actual output
                input_relapse_for_eval = []
                for sample_index in range(0, len(list_top_n_gene_relapse_sorted_for_eval)):
                    list_each_sample = []
                    for element_id in range(0, len(list_top_n_gene_relapse_sorted_for_eval[sample_index])):
                        if (element_id in feature_set):
                            list_each_sample.append(list_top_n_gene_relapse_sorted_for_eval[sample_index][element_id])
                    input_relapse_for_eval.append(list_each_sample)

                input_no_relapse_for_eval = []
                for sample_index in range(0, len(list_top_n_gene_no_relapse_sorted_for_eval)):
                    list_each_sample = []
                    for element_id in range(0, len(list_top_n_gene_no_relapse_sorted_for_eval[sample_index])):
                        if (element_id in feature_set):
                            list_each_sample.append(list_top_n_gene_no_relapse_sorted_for_eval[sample_index][element_id])
                    input_no_relapse_for_eval.append(list_each_sample)

                input_testing_data_for_eval = []
                for sample_index in range(0, len(list_first_layer_top_n_test_sorted)):
                    list_each_sample = []
                    for element_id in range(0, len(list_first_layer_top_n_test_sorted[sample_index])):
                        if (element_id in feature_set):
                            list_each_sample.append(list_first_layer_top_n_test_sorted[sample_index][element_id])
                    input_testing_data_for_eval.append(list_each_sample)

                list_actual_output_for_eval = calculate.lda(input_testing_data_for_eval, input_relapse_for_eval, input_no_relapse_for_eval)

                # calculate AUC score
                auc_score_for_eval = roc_auc_score(list_desired_output_for_eval, list_actual_output_for_eval)

                print("#### Evaluation of " + str(first_layer_test_index + 1) + " - fold ####")
                print("Feature Set : " + str(feature_set_name))
                print("Actual Output : " + str(list_actual_output_for_eval))
                print("Desired Output : " + str(list_desired_output_for_eval))
                print("AUC ROC score = " + str(auc_score_for_eval))

                # record ending time of this iteration
                end_epoch_time = time.time()
                time_elapse_epoch_second = end_epoch_time - start_epoch_time
                time_elapse_epoch_minute = time_elapse_epoch_second / 60
                time_elapse_epoch_hour = time_elapse_epoch_minute / 60

                # write output to an output file
                result_file.write("Fold : " + str(first_layer_test_index + 1) + "\n")
                result_file.write("Feature Set : " + str(feature_set_name) + "\n")
                result_file.write("Actual Output : " + str(list_actual_output_for_eval) + "\n")
                result_file.write("Desired Output : " + str(list_desired_output_for_eval) + "\n")
                result_file.write("AUC ROC Score : " + str(auc_score_for_eval) +  "\n")
                result_file.write("\n")
        result_file.write("Maximum AUC ROC score of feature in each fold = " + str(list_max_auc) + "\n")
        result_file.write("Time Elapse : " + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours)\n")
        print("Time Elapse : " + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours)\n")
        print("Maximum AUC ROC score of feature in each fold  = " + str(list_max_auc))
    # record end time
    end_time = time.time()
    time_elapse_second = end_time - start_time
    time_elapse_minute = time_elapse_second / 60
    time_elapse_hour = time_elapse_minute / 60
    print("Time Elapse : " + str(time_elapse_minute) + " minutes (" + str(time_elapse_hour) + " hours)")
    result_file.write("Total Time Elapse : " + str(time_elapse_minute) + " minutes (" + str(time_elapse_hour) + " hours)\n")

    result_file.close()
                         
if __name__ == '__main__':
    main()