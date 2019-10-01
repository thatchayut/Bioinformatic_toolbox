from scipy import stats
from sklearn.metrics import roc_auc_score
from copy import deepcopy

import math
import calculate
import random
import time
import data_handler
import pandas as pd
import numpy as np
import sys

import config as cf

def test():
    print("testestest")

def gene_based():
    # record start time
    start_time = time.time()

    # acquire required files
    # a file containing mapping between probes IDs and sample
    file_training_input_name = cf.gene_based["file_training_input_name"]
    # check whether this file is valid
    if data_handler.validateFile(file_training_input_name) is False :
        sys.exit(0)

    # number of rows of this file to be read
    row_to_read = cf.gene_based["row_to_read"]
    # check whether number of rows is valid
    if data_handler.validateRowToRead(row_to_read) is False :
        sys.exit(0)
    row_to_read = int(row_to_read)

    # name of a file containing mapping between samples and their class
    file_training_output_name = cf.gene_based["file_training_output_name"]
    # check whether this file is valid
    if data_handler.validateFile(file_training_output_name) is False :
        sys.exit(0)

    # acquire required information to conduct an experiment
    # number of epochs
    epochs = cf.gene_based["epochs"]
    # check whether number of epochs is valid
    if data_handler.validateEpochs(epochs) is False :
        sys.exit(0)
    epochs = int(epochs)

    # prepare data

    # row_to_read = 22283
    # file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_input = pd.read_csv(file_training_input_name, nrows = row_to_read)

    # consider non-relapse and relapse (not in specific period of time)
    file_training_output = pd.read_csv(file_training_output_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

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
    
    # separate data into 2 classes
    # consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['0'])]
    
    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)
    
    # shuffle data to make each chunk does not depend on sample order
    random.shuffle(list_sample_relapse) 
    random.shuffle(list_sample_no_relapse)

    # number of folds
    num_of_folds = cf.gene_based["num_of_folds"]
    # check whether number of folds is valid
    if data_handler.validateNumofFolds(num_of_folds, list_sample_relapse, list_sample_no_relapse) is False :
        sys.exit(0)
    num_of_folds = int(num_of_folds)

    # number of top-ranked features
    num_of_ranked_genes = cf.gene_based["num_of_ranked_genes"]
    # check whether number of top-ranked features is valid
    if data_handler.validateNumofRankedGenes(num_of_ranked_genes, row_to_read) is False :
        sys.exit(0)
    num_of_ranked_genes = int(num_of_ranked_genes)

    # prepare an output file
    output_file_name  = "gene_based_" + str(epochs) + "ep_" + str(num_of_folds) + "f_" + "top" + str(num_of_ranked_genes) + "_" + cf.gene_based['output_file_name']

    result_file = open("./result/" +str(output_file_name) + ".txt", "w+")

    # write file information 
    result_file.write("Dataset : " + str(file_training_input_name) +"\n")
    result_file.write("Number of epochs : " + str(epochs) + "\n")
    result_file.write("Number of folds : " + str(num_of_folds) + "\n")
    result_file.write("Number of top-ranked features : " + str(num_of_ranked_genes) + "\n")
    result_file.write("Number of samples in class relapse : " + str(len(list_sample_relapse)) + "\n")
    result_file.write("Number of samples in class non-relapse : " + str(len(list_sample_no_relapse)) + "\n")
    result_file.write("\n")

    print(" Number of samples in class relapse : " + str(len(list_sample_relapse)))
    print(" Number of samples in class non-relapse : " + str(len(list_sample_no_relapse)))

    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    # list to collect feature counter
    list_feature_counter = []

    # start an experiment
    for epoch_count in range(0, int(epochs)):
        # preparing for cross validation
        start_epoch_time = time.time()
        result_file.write("#################################### Epoch : " + str(epoch_count + 1) + " ####################################\n")
        print("#################################### Epoch : " + str(epoch_count + 1) + " ####################################\n")

        # split data into k parts
        chunk_relapse_size = math.ceil(len(list_sample_relapse) / num_of_folds)
        chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse) / num_of_folds)

        chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_relapse_size))
        chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_no_relapse_size))

        # check if number of chunks is valid
        check_valid_num_of_chunks, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)
        if check_valid_num_of_chunks is False:
            print(" WARNING : Number of chunks in both classes must be the same")
            sys.exit(0)

        # list to collect maximun AUC in each fold
        list_max_auc = []

        # list and variable to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        print(" # Process : Cross-validation")
        # do only if number of chunks of both datasets are equal
        # if (check_valid_num_of_chunks == True):
        for first_layer_test_index in range(0, num_of_chunks):
            feature_set = []
            feature_set_name = []
            top_n_genes_name_for_eval = []

            # get testing data from each class
            first_layer_test_relapse = chunk_list_relapse[first_layer_test_index]
            first_layer_test_no_relapse = chunk_list_no_relapse[first_layer_test_index]
            print("\n------------------------------------------ K : " + str(first_layer_test_index + 1) + " of Epoch " + str(epoch_count + 1) + " --------------------------------")

            # get training data
            # first layer 
            first_layer_train_relapse = []
            for first_layer_train_index in range(0, num_of_chunks):
                if (chunk_list_relapse[first_layer_train_index] is not first_layer_test_relapse):
                    first_layer_train_relapse.append(chunk_list_relapse[first_layer_train_index]

            first_layer_train_no_relapse = []
            for first_layer_train_index in range(0, num_of_chunks):
                if (chunk_list_no_relapse[first_layer_train_index] is not first_layer_test_no_relapse):
                    first_layer_train_no_relapse.append(chunk_list_no_relapse[first_layer_train_index])
            
            # merge all element in the same class
            second_list_sample_relapse = []
            for i in range(0, len(first_layer_train_relapse)):
                second_list_sample_relapse.extend(first_layer_train_relapse[i])
            print(" Samples in class relapse used as trainning set = " + str(second_list_sample_relapse) + "\n")

            second_list_sample_no_relapse = []
            for i in range(0, len(first_layer_train_no_relapse)):
                second_list_sample_no_relapse.extend(first_layer_train_no_relapse[i])
            print(" Samples in class non-relapse used as training set : " + str(second_list_sample_no_relapse) + "\n")

            # splitting lists to use them as marker evaluation set and feature selection set
            # given that we separate it into 3 parts
            # divide training set into 3 parts (2/3 for marker evaluation and 1/3 for feature selection)
            print(" Process : Feature selection")
            second_num_of_fold = 3
            second_chunk_relapse_size = math.ceil(len(second_list_sample_relapse) / second_num_of_fold)
            second_chunk_no_relapse_size = math.ceil(len(second_list_sample_no_relapse) / second_num_of_fold)

            second_chunk_list_relapse = list(calculate.chunks(second_list_sample_relapse, second_chunk_relapse_size))
            second_chunk_list_no_relapse = list(calculate.chunks(second_list_sample_no_relapse, second_chunk_no_relapse_size))

            # check if number of chunks is valid
            second_check_valid_num_of_chunks, second_num_of_chunks = calculate.checkEqualListSize(second_chunk_list_relapse, second_chunk_list_no_relapse)
            if second_check_valid_num_of_chunks is False:
                print(" WARNING : Number of chunks in both classes must be the same")
                sys.exit(0)

            # do only if number of chunks of both datasets are equal
            # if (second_check_valid == True):
            second_layer_test_index = random.randint(0, second_num_of_chunks - 1)

            # get testing data from each class
            second_layer_test_relapse =  second_chunk_list_relapse[second_layer_test_index]
            second_layer_test_no_relapse = second_chunk_list_no_relapse[second_layer_test_index]
            
            # separate training dataset from testing dataset to use in t-test ranking
            second_layer_train_relapse = []
            for second_layer_train_index in range(0, second_num_of_chunks):
                if (second_chunk_list_relapse[second_layer_train_index] is not second_layer_test_relapse):
                    second_layer_train_relapse.append(second_chunk_list_relapse[second_layer_train_index])
            
            second_layer_train_no_relapse = []
            for second_layer_train_index in range(0, second_num_of_chunks):
                if (second_chunk_list_no_relapse[second_layer_train_index] is not second_layer_test_no_relapse):
                    second_layer_train_no_relapse.append(second_chunk_list_no_relapse[second_layer_train_index])

            # prepare dataset for conducting t-test
            # merge all samples in the same class
            # samples in class relapse used as marker evaluation set
            ttest_list_sample_relapse = []
            for i in range(0, len(second_layer_train_relapse)):
                ttest_list_sample_relapse.extend(second_layer_train_relapse[i])

            # samples in class non-relapse used as marker evaluation set
            ttest_list_sample_no_relapse = []
            for i in range(0, len(second_layer_train_no_relapse)):
                ttest_list_sample_no_relapse.extend(second_layer_train_no_relapse[i])

            # get gene expression for each gene from samples with relapse
            list_gene_exp_relapse = []
            for i in range(0, row_to_read):
                gene_exp_relapse = []
                for column in  file_training_input.loc[i, ttest_list_sample_relapse]:
                    gene_exp_relapse.append(column)
                list_gene_exp_relapse.append(gene_exp_relapse)
            
            # get gene expression for each gene from samples with no relapse
            list_gene_exp_no_relapse = []
            for i in range(0, row_to_read):
                gene_exp_no_relapse = []
                for column in  file_training_input.loc[i, ttest_list_sample_no_relapse]:
                    gene_exp_no_relapse.append(column)
                list_gene_exp_no_relapse.append(gene_exp_no_relapse)
            
            # conducting t-test
            print(" # Process : Calculating t-score")
            ttest_result = []
            for i in range(0, row_to_read):      
                score = []
                # get absolute magnitude of t-test value
                abs_ttest_value = math.fabs(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[0])
                p_value = stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[1]
                # add element with this format (gene_order_id, ttest_value)
                score.append(i)
                score.append(abs_ttest_value)
                ttest_result.append(score)
            # ranking elements using their t-test value in descending order
            ttest_result.sort(key=lambda x: x[1], reverse=True)

            # create list of ranked gene
            ranked_gene = []
            for i in range(0, len(ttest_result)):
                gene_order_id = ttest_result[i][0]

                ranked_gene.append(list_gene_name[gene_order_id][1])
            
            # rank gene id of each sample in training data
            # for class 'relapse'
            col_to_read_relapse = ["ID_REF"]
            col_to_read_relapse.extend(ttest_list_sample_relapse)
            file_training_input_relapse = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_relapse)
            top_n_genes_relapse = file_training_input_relapse.loc[file_training_input_relapse['ID_REF'].isin(top_n_genes_name)]
            top_n_genes_relapse['gene_id'] = top_n_genes_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))
            top_n_genes_relapse_sorted  = top_n_genes_relapse.sort_values(by = ['gene_id'])
            top_n_genes_relapse_sorted.drop(columns = 'gene_id', inplace = True)

            top_n_genes_relapse_sorted_train = top_n_genes_relapse_sorted
            top_n_genes_relapse_sorted_train.drop(columns = 'ID_REF', inplace = True)

            # for class 'no relapse'
            col_to_read_no_relapse = ["ID_REF"]
            col_to_read_no_relapse.extend(ttest_list_sample_no_relapse)
            file_training_input_no_relapse = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_no_relapse)
            top_n_genes_no_relapse = file_training_input_no_relapse.loc[file_training_input_no_relapse['ID_REF'].isin(top_n_genes_name)]
            top_n_genes_no_relapse['gene_id'] = top_n_genes_no_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))
            top_n_genes_no_relapse_sorted  = top_n_genes_no_relapse.sort_values(by = ['gene_id'])
            top_n_genes_no_relapse_sorted.drop(columns = 'gene_id', inplace = True)

            top_n_genes_no_relapse_sorted_train = top_n_genes_no_relapse_sorted
            top_n_genes_no_relapse_sorted_train.drop(columns = 'ID_REF', inplace = True)

            # Preparing testing data for feature selection
            second_layer_test_all = []
            second_layer_test_all.extend(second_layer_test_relapse)
            second_layer_test_all.extend(second_layer_test_no_relapse)  

            # prepare expect output for testing data
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

            # use top-rank feature as the first feature in lda classifier
            # prepare list for input 
            # list of all input data (testing data)
            list_second_layer_top_n_test_sorted = []
            for column in range(0, len(top_n_test_sorted)):
                list_each_sample = []
                for element in top_n_test_sorted.iloc[column]:
                    list_each_sample.append(element)
                list_second_layer_top_n_test_sorted.append(list_each_sample)
            list_second_layer_top_n_test_sorted = list(np.transpose(list_second_layer_top_n_test_sorted))

            # desired output for testing data
            second_layer_test_output = training_output.loc[training_output['GEO asscession number'].isin(second_layer_test_all)]
            # sorting data according to its order in testing data
            list_sample_to_read = list(second_layer_top_n_test_sorted.columns.values)

            second_layer_test_output['sample_id'] = second_layer_test_output['GEO asscession number'].apply(lambda name: list_sample_to_read.index(name))
            second_layer_test_output = second_layer_test_output.sort_values(by = ['sample_id'])
            second_layer_test_output.drop(columns = 'sample_id', inplace = True)

            # create list of desired output
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

            # list of gene expression and sample of class 'no relapse'
            list_top_n_gene_no_relapse_sorted = []
            for column in range(0, len(top_n_genes_no_relapse_sorted_train)):
                list_each_sample = []
                for element in top_n_genes_no_relapse_sorted_train.iloc[column]:
                    list_each_sample.append(element)
                list_top_n_gene_no_relapse_sorted.append(list_each_sample)
            list_top_n_gene_no_relapse_sorted = list(np.transpose(list_top_n_gene_no_relapse_sorted))

            # find set of genes to be used as a feature using sequential feature selection
            print(" # Process : Sequential Forward Selection (SFS)")
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

                        input_no_relapse = []
                        for sample_index in range(0, len(list_top_n_gene_no_relapse_sorted)):
                            list_each_sample = []
                            for element_id in range(0, len(list_top_n_gene_no_relapse_sorted[sample_index])):
                                if (element_id in gene_order_test):
                                    list_each_sample.append(list_top_n_gene_no_relapse_sorted[sample_index][element_id])
                            input_no_relapse.append(list_each_sample)

                        input_testing_data = []
                        for sample_index in range(0, len(list_second_layer_top_n_test_sorted)):
                            list_each_sample = []
                            for element_id in range(0, len(list_second_layer_top_n_test_sorted[sample_index])):
                                if (element_id in gene_order_test):
                                    list_each_sample.append(list_second_layer_top_n_test_sorted[sample_index][element_id])
                            input_testing_data.append(list_each_sample)

                        list_actual_output = calculate.lda(input_testing_data, input_relapse, input_no_relapse)

                        # calculate AUC score
                        auc_score = roc_auc_score(list_desired_output, list_actual_output)

                        if (auc_score > max_auc_score):
                            max_auc_score = auc_score
                            gene_index_in_list = i
                            if max_auc_score not in list_auc:
                                list_auc.append(max_auc_score)

                    # do not add gene that already exists in a feature
                    if (gene_index_in_list not in gene_order):
                        gene_order.extend([gene_index_in_list])                       
                    count_iteration += 1  

            list_max_auc.append(max(list_auc))        
            gene_order.sort()       

            # get gene_name
            gene_order_name = []
            for element in gene_order:
                gene_order_name.append(top_n_genes_name[element])

            # copy required data to be used in evaluation
            top_n_genes_name_for_eval = deepcopy(top_n_genes_name)
            feature_set  = deepcopy(gene_order)
            feature_set_name = deepcopy(gene_order_name)

        # count feature frequency
        if (int(epochs) > 1):
            for feature_index in range(0, len(feature_set_name)):
                # if list feature counter is empty
                if not list_feature_counter:
                    feature_counter = []
                    feature_name = feature_set_name[feature_index]
                    feature_frequency = 1

                    feature_counter.append(feature_name)
                    feature_counter.append(feature_frequency)

                    list_feature_counter.append(feature_counter)
                else:
                    feature_name = feature_set_name[feature_index]

                    # check if this feature exist in the feature counter list
                    check_found = False
                    for feature_counter_index in range(0, len(list_feature_counter)):
                        feature_counter_name = list_feature_counter[feature_counter_index][0]

                        if (feature_name == feature_counter_name):
                            feature_frequency = list_feature_counter[feature_counter_index][1]
                            feature_frequency += 1

                            list_feature_counter[feature_counter_index][1] = feature_frequency
                            check_found = True
                    
                    # if this feature is not exist in a list feature counter
                    if (check_found == False):
                        feature_counter = []
                        feature_name = feature_set_name[feature_index]
                        feature_frequency = 1

                        feature_counter.append(feature_name)
                        feature_counter.append(feature_frequency)

                        list_feature_counter.append(feature_counter)

        # preparing data for evaluation and creating classifier
        print(" # Process : Prepare classifiers and testing data")

        # for class 'relapse'
        col_to_read_relapse_for_eval = ["ID_REF"]
        col_to_read_relapse_for_eval.extend(second_list_sample_relapse)
        file_training_input_relapse_for_eval = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_relapse_for_eval)
        top_n_genes_relapse_for_eval = file_training_input_relapse.loc[file_training_input_relapse['ID_REF'].isin(feature_set_name)]
        top_n_genes_relapse_for_eval['gene_id'] = top_n_genes_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))
        top_n_genes_relapse_sorted_for_eval  = top_n_genes_relapse_for_eval.sort_values(by = ['gene_id'])
        top_n_genes_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)
        top_n_genes_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

        # for class 'no relapse'
        col_to_read_no_relapse_for_eval = ["ID_REF"]
        col_to_read_no_relapse_for_eval.extend(second_list_sample_no_relapse)
        file_training_input_no_relapse_for_eval = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_no_relapse_for_eval)
        top_n_genes_no_relapse_for_eval = file_training_input_no_relapse_for_eval.loc[file_training_input_no_relapse_for_eval['ID_REF'].isin(feature_set_name)]
        top_n_genes_no_relapse_for_eval['gene_id'] = top_n_genes_no_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))
        top_n_genes_no_relapse_sorted_for_eval  = top_n_genes_no_relapse_for_eval.sort_values(by = ['gene_id'])
        top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)
        top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

        # get all testing data
        first_layer_test_all = []
        first_layer_test_all.extend(first_layer_test_relapse)
        first_layer_test_all.extend(first_layer_test_no_relapse)  

        # get only genes which are in a feature set
        col_to_read_first_layer_test_gene = ["ID_REF"]
        col_to_read_first_layer_test_gene.extend(first_layer_test_all)
        first_layer_test_gene = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read, usecols = col_to_read_first_layer_test_gene)
        first_layer_top_n_test = first_layer_test_gene.loc[first_layer_test_gene['ID_REF'].isin(feature_set_name)]
        first_layer_top_n_test['gene_id'] = first_layer_top_n_test['ID_REF'].apply(lambda name: feature_set_name.index(name))   
        first_layer_top_n_test_sorted = first_layer_top_n_test.sort_values(by = ['gene_id'])
        first_layer_top_n_test_sorted.drop(columns = 'gene_id', inplace = True)

        top_n_test_sorted_for_eval = first_layer_top_n_test_sorted
        top_n_test_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

        # prepare list for input 
        # list of all input data (testing data)
        list_first_layer_top_n_test_sorted = []
            for column in range(0, len(top_n_test_sorted_for_eval)):
                list_each_sample = []
                for element in top_n_test_sorted_for_eval.iloc[column]:
                    list_each_sample.append(element)

                list_first_layer_top_n_test_sorted.append(list_each_sample)
            list_first_layer_top_n_test_sorted = list(np.transpose(list_first_layer_top_n_test_sorted))
        
        # desired output for testing data
        first_layer_test_output = training_output.loc[training_output['GEO asscession number'].isin(first_layer_test_all)]

        # sorting data according to its order in testing data
        list_sample_to_read_for_eval = list(first_layer_top_n_test_sorted.columns.values)
        first_layer_test_output['sample_id'] = first_layer_test_output['GEO asscession number'].apply(lambda name: list_sample_to_read_for_eval.index(name))
        first_layer_test_output = first_layer_test_output.sort_values(by = ['sample_id'])
        first_layer_test_output.drop(columns = 'sample_id', inplace = True)

        # create list of desired output
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

        # list of gene expression and sample of class 'no relapse' for evaluation
        list_top_n_gene_no_relapse_sorted_for_eval = []
        for column in range(0, len(top_n_genes_no_relapse_sorted_for_eval)):
            list_each_sample = []
            for element in top_n_genes_no_relapse_sorted_for_eval.iloc[column]:
                list_each_sample.append(element)
            list_top_n_gene_no_relapse_sorted_for_eval.append(list_each_sample)
        list_top_n_gene_no_relapse_sorted_for_eval = list(np.transpose(list_top_n_gene_no_relapse_sorted_for_eval))    

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

        






