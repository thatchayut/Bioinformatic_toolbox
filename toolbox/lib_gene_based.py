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
            