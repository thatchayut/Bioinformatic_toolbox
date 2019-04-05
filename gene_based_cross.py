import pandas as pd
import math
import calculate
import random
import numpy as np
import time
import os
import add_ons
from scipy import stats
from sklearn.metrics import roc_auc_score
from copy import deepcopy



def main():
    # record start time
    start_time = time.time()

    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Method : Gene-Based Classification")
    print(" # Experiment : Cross-Dataset")
    print(" # This method requires 2 datasets.")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains mapping between gene probe IDs and samples of the first dataset")
    print(" #   [2] Number of rows of the file containing mapping between gene probe IDs and samples of the first dataset to be read")
    print(" #   [3] A file contains mapping between samples and their class of the first dataset")  
    print(" #   [4] A file contains mapping between gene probe IDs and samples of the second dataset")
    print(" #   [5] Number of rows of the file contains mapping between gene probe IDs and samples of the second dataset to be read")
    print(" #   [6] A file contains mapping samples and their class of the second dataset")
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked to provide required information to conduct an experiment including")
    print(" #   [1] Number of epochs")
    print(" #   [2] Number of folds")
    print(" #   [3] Number of top-ranked feature")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()

    # prepare variables
    file_gene_first_dataset_name = None
    file_gene_second_dataset_name = None

    row_to_read_file_gene_first_dataset = None
    row_to_read_file_gene_second_dataset = None
    
    file_output_first_dataset_name = None
    file_output_second_dataset_name = None

    epoch = None
    num_of_folds = None
    number_of_ranked_gene = None

    file_name = None


    print(" # Enter required information about the first dataset ")
    print(" 1. Enter name of a file containing mapping between probes IDs and samples of the first dataset ")
    file_gene_first_dataset_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_gene_first_dataset = input(" Number of rows : ")
        if (row_to_read_file_gene_first_dataset.isnumeric() == False):
            print(" WARNING : Number of rows must be numeric.")
        elif (int(row_to_read_file_gene_first_dataset) < 1):
            print("WARNING : Number of rows cannot be lower than 1.")
        else:
            break
    row_to_read_file_gene_first_dataset = int(row_to_read_file_gene_first_dataset)
    print()

    print(" 3. Enter name of file containing mapping between samples and their class of the first dataset")
    file_output_first_dataset_name = add_ons.checkFileValid()
    print()

    print(" # Enter required information about the second dataset ")
    print(" 1. Enter name of a file containing mapping between probes IDs and samples of the second dataset ")
    file_gene_second_dataset_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_gene_second_dataset = input(" Number of rows : ")
        if(row_to_read_file_gene_second_dataset.isnumeric() == False):
            print(" WARNING : Number of rows must be numeric.")
        elif (int(row_to_read_file_gene_second_dataset) < 1):
            print(" WARNING : Number of rows cannot be lower than 1.")
        else:
            break
    row_to_read_file_gene_second_dataset = int(row_to_read_file_gene_second_dataset)
    print()

    print(" 3. Enter name of a containing mapping between samples and their class of the second dataset")
    file_output_second_dataset_name = add_ons.checkFileValid()
    print()

    # prepare data
    # for 1st dataset
    # default row_to_read_file_gene_first_dataset for "GSE2034-22071 (edited).csv" = 22283
    file_gene_first_dataset = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset)

    # default file_output_first_dataset_name for "GSE2034-22071 (edited).csv" = "mapping_sample_to_class_gse2034.csv"
    file_output_first_dataset = pd.read_csv(file_output_first_dataset_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # for 2nd dataset
    # default row_to_read_file_second_dataset for "GSE3494_GPL96.csv" = 22283
    file_gene_second_dataset = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset) 

    # default file_output_second_dataset_name for "GSE3494_GPL96.csv" = "mapping_sample_to_class_gse3494.csv"
    file_output_second_dataset = pd.read_csv(file_output_second_dataset_name, usecols = ["GEO asscession number", "relapse (1=True)"])

    # get gene order id with its name
    list_gene_name_first_dataset = []
    for i in range(0, row_to_read_file_gene_first_dataset):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_gene_first_dataset.loc[i, "ID_REF"])
        list_gene_name_first_dataset.append(gene_name)
    
    # consider non-relapse and relapse (not in specific period of time)
    # for 1st dataset
    sample_relapse_first_dataset = file_output_first_dataset.loc[file_output_first_dataset['relapse (1=True)'].isin(['1'])]
    sample_no_relapse_first_dataset = file_output_first_dataset.loc[file_output_first_dataset['relapse (1=True)'].isin(['0'])]

    # for 2nd dataset
    sample_relapse_second_dataset = file_output_second_dataset.loc[file_output_second_dataset['relapse (1=True)'].isin(['1'])]
    sample_no_relapse_second_dataset = file_output_second_dataset.loc[file_output_second_dataset['relapse (1=True)'].isin(['0'])]

    # add GEO asscession number to each list
    # for 1st dataset
    list_sample_relapse_first_dataset = []
    for element in sample_relapse_first_dataset.loc[:, 'GEO asscession number']:
        list_sample_relapse_first_dataset.append(element)

    list_sample_no_relapse_first_dataset = []
    for element in sample_no_relapse_first_dataset.loc[:, 'GEO asscession number']:
        list_sample_no_relapse_first_dataset.append(element)
    
    # for 2nd dataset
    list_sample_relapse_second_dataset = []
    for element in sample_relapse_second_dataset.loc[:, 'GEO asscession number']:
        list_sample_relapse_second_dataset.append(element)
    
    list_sample_no_relapse_second_dataset = []
    for element in sample_no_relapse_second_dataset.loc[:, 'GEO asscession number']:
        list_sample_no_relapse_second_dataset.append(element)

    print(" # Enter required information to conduct an experiment")
    print(" 1. Enter number of epochs ")
    while True:
        epoch = input(" Epochs : ")

        if (epoch.isnumeric() == False):
            print(" WARNING : Number of epochs must be numeric.")
        elif (int(epoch) <= 0):
            print(" WARINING : Number of epochs must be greater than 0.")
        else:
            break
    print()

    print(" 2. Enter number of folds ")
    while True:
        num_of_folds = input(" Number of folds: ")
        if (num_of_folds.isnumeric() == False):
            print(" WARNING : Number of folds must be numeric")
        
        # these conditions are not available in mock-up
        elif(int(num_of_folds) > len(list_sample_relapse_second_dataset)):
            print("WARNING : Number of folds exceeds the size of samples in class relapse in the second dataset.")
        elif(int(num_of_folds) > len(list_sample_no_relapse_second_dataset)):
            print("WARNING : Number of folds exceeds the size of samples in clss non-relapse in the second dataset.")

        elif(int(num_of_folds) <= 1):
            print(" WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)    
    print()

    print(" 3. Enter number of top-ranked features")
    while True:
        number_of_ranked_gene = input(" Number of top-ranked features: ")
        if(number_of_ranked_gene.isnumeric() == False):
            print(" WARNING : Number of top-ranked features must be numeric.")

        # these conditions are not available in mock-up
        elif(int(number_of_ranked_gene) > row_to_read_file_gene_first_dataset):
            print(" WARINING : Number of top-ranked features must not exceed available genes from the first file.")

        elif(int(number_of_ranked_gene) <= 0):
            print(" WARNING : Number of top-ranked features must not be lower than or equal to 0.")    
        else:
            break
    print()

    file_name = input(" # Enter name of an output file : ")

    # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    # record dataset
    result_file.write(" The first dataset : " + str(file_gene_first_dataset_name) + "\n")
    result_file.write(" The second dataset : " + str(file_gene_second_dataset_name) + "\n")
    result_file.write("\n")
    
    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    for epoch_count in range(0, int(epoch)):
        start_epoch_time = time.time()
        result_file.write("#################################### Epoch : " + str(epoch_count + 1) + " ####################################\n")
        print("#################################### Epoch : " + str(epoch_count + 1) + " ####################################")

        # shuffle it to make it flexible for epoch changed
        # for 1st dataset
        random.shuffle(list_sample_relapse_first_dataset)
        random.shuffle(list_sample_no_relapse_first_dataset)

        # for 2nd dataset
        random.shuffle(list_sample_relapse_second_dataset)
        random.shuffle(list_sample_no_relapse_second_dataset)

        # conduct feature selection on 1st dataset
        # split data into 3 parts
        num_of_chunk_feature_selection = 3
        chunk_relapse_size = math.ceil(len(list_sample_relapse_first_dataset) / num_of_chunk_feature_selection)
        chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse_first_dataset) / num_of_chunk_feature_selection)

        chunk_list_relapse = list(calculate.chunks(list_sample_relapse_first_dataset, chunk_relapse_size))
        print()
        print(" Number of chunks in class relapse of the first dataset : " + str(len(chunk_list_relapse)) + "\n")

        chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse_first_dataset, chunk_no_relapse_size))
        print(" Number of chunks in class non-relapse of the first dataset : " + str(len(chunk_list_no_relapse)) + "\n")

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # variable to collect data from sfs
        feature_set = []
        feature_set_name = []
        top_n_genes_name_for_eval = []

        print()
        print(" # The first dataset is used in feature selection.\n")
        print(" # Process : Feature selection")
        # do only if number of chunks of both datasets are equal
        if (check_valid == True):
            # random a chunk of data to be use as a feature selection set
            feature_selection_index = random.randint(0, num_of_chunks - 1)

            # get feature selection set from each class
            feature_selection_relapse =  chunk_list_relapse[feature_selection_index]
            feature_selection_no_relapse = chunk_list_no_relapse[feature_selection_index]

            # separate marker evaluation dataset from feature selection dataset
            # for class "relapse"
            marker_evaluation_relapse = []
            for marker_evaluation_index in range(0, num_of_chunks):
                if (chunk_list_relapse[marker_evaluation_index] is not feature_selection_relapse):
                    marker_evaluation_relapse.append(chunk_list_relapse[marker_evaluation_index])

            # for class "non-relapse"
            marker_evaluation_no_relapse = []
            for marker_evaluation_index in range(0, num_of_chunks):
                if (chunk_list_no_relapse[marker_evaluation_index] is not feature_selection_no_relapse):
                    marker_evaluation_no_relapse.append(chunk_list_no_relapse[marker_evaluation_index])

            # merge all samples in the same class
            # for class "relapse"
            list_sample_relapse_marker_evaluation = []
            for i in range(0, len(marker_evaluation_relapse)):
                list_sample_relapse_marker_evaluation.extend(marker_evaluation_relapse[i])
            print(" Samples in class relapse used as marker evaluation set : " + str(list_sample_relapse_marker_evaluation) + "\n")

            # for class "non-relapse"
            list_sample_no_relapse_marker_evaluation = []
            for i in range(0, len(marker_evaluation_no_relapse)):
                list_sample_no_relapse_marker_evaluation.extend(marker_evaluation_no_relapse[i])
            print(" Samples in class non-relapse used as marker evaluation set : " + str(list_sample_no_relapse_marker_evaluation) + "\n")

            # get gene expression for each gene from samples with relapse
            # for class "relapse"
            list_gene_expression_relapse = []
            for i in range(0, row_to_read_file_gene_first_dataset):
                gene_expression_relapse = []
                for column in  file_gene_first_dataset.loc[i, list_sample_relapse_marker_evaluation]:
                    gene_expression_relapse.append(column)
                list_gene_expression_relapse.append(gene_expression_relapse)
            
            # for class "non-relapse"
            list_gene_expression_no_relapse = []
            for i in range(0, row_to_read_file_gene_first_dataset):
                gene_expression_no_relapse = []
                for column in  file_gene_first_dataset.loc[i, list_sample_no_relapse_marker_evaluation]:
                    gene_expression_no_relapse.append(column)
                list_gene_expression_no_relapse.append(gene_expression_no_relapse)
            
            print(" # Process : Calculating t-score")
            # conducting t-test
            ttest_result = []
            for i in range(0, row_to_read_file_gene_first_dataset):      
                score = []

                # get absolute magnitude of t-test value
                abs_ttest_value = math.fabs(stats.ttest_ind(list_gene_expression_relapse[i], list_gene_expression_no_relapse[i], equal_var = False)[0])
                p_value = stats.ttest_ind(list_gene_expression_relapse[i], list_gene_expression_no_relapse[i], equal_var = False)[1]
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
                ranked_gene.append(list_gene_name_first_dataset[gene_order_id][1])
            
            # show top ranked feature
            top_n_genes_name = []
            print(" #### t-score ranking ####")
            for i in range(0, int(number_of_ranked_gene)):
                top_n_genes_name.append(ranked_gene[i])
                print(" " + str(ranked_gene[i]) + " => " + " t-score : " + str(ttest_result[i][1]))
            
            # rank gene id of each sample in training data
            # for class 'relapse'
            col_to_read_relapse = ["ID_REF"]
            col_to_read_relapse.extend(list_sample_relapse_marker_evaluation)

            file_training_input_relapse = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset, usecols = col_to_read_relapse)
            
            top_n_genes_relapse = file_training_input_relapse.loc[file_training_input_relapse['ID_REF'].isin(top_n_genes_name)]
            top_n_genes_relapse['gene_id'] = top_n_genes_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))

            top_n_genes_relapse_sorted  = top_n_genes_relapse.sort_values(by = ['gene_id'])
            top_n_genes_relapse_sorted.drop(columns = 'gene_id', inplace = True)

            top_n_genes_relapse_sorted_train = top_n_genes_relapse_sorted
            top_n_genes_relapse_sorted.drop(columns = 'ID_REF', inplace = True)

            # for class 'no relapse'
            col_to_read_no_relapse = ["ID_REF"]
            col_to_read_no_relapse.extend(list_sample_no_relapse_marker_evaluation)

            file_training_input_no_relapse = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset, usecols = col_to_read_no_relapse)

            top_n_genes_no_relapse = file_training_input_no_relapse.loc[file_training_input_no_relapse['ID_REF'].isin(top_n_genes_name)]
            top_n_genes_no_relapse['gene_id'] = top_n_genes_no_relapse['ID_REF'].apply(lambda name: top_n_genes_name.index(name))

            top_n_genes_no_relapse_sorted  = top_n_genes_no_relapse.sort_values(by = ['gene_id'])
            top_n_genes_no_relapse_sorted.drop(columns = 'gene_id', inplace = True)

            top_n_genes_no_relapse_sorted_train = top_n_genes_no_relapse_sorted
            top_n_genes_no_relapse_sorted_train.drop(columns = 'ID_REF', inplace = True)

            # Preparing testing data for feature selection
            feature_selection_all = []
            feature_selection_all.extend(feature_selection_relapse)
            feature_selection_all.extend(feature_selection_no_relapse)   

            col_to_read_feature_selection = ["ID_REF"]
            col_to_read_feature_selection.extend(feature_selection_all)

            file_feature_selection = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset, usecols = col_to_read_feature_selection)
            file_feature_selection_top_n_test = file_feature_selection.loc[file_feature_selection['ID_REF'].isin(top_n_genes_name)]

            file_feature_selection_top_n_test['gene_id'] = file_feature_selection_top_n_test['ID_REF'].apply(lambda name: top_n_genes_name.index(name))
            file_feature_selection_top_n_test_sorted = file_feature_selection_top_n_test.sort_values(by = ['gene_id'])
            file_feature_selection_top_n_test_sorted.drop(columns = 'gene_id', inplace = True)

            top_n_test_sorted = file_feature_selection_top_n_test_sorted
            top_n_test_sorted.drop(columns = 'ID_REF', inplace = True)

            # use top-rank feature as the first feature in lda classifier
            # prepare list for input 
            # list of all input data (testing data)
            list_feature_selection_top_n_test_sorted = []
            for column in range(0, len(top_n_test_sorted)):
                list_each_sample = []
                for element in top_n_test_sorted.iloc[column]:
                    list_each_sample.append(element)
                list_feature_selection_top_n_test_sorted.append(list_each_sample)
            list_feature_selection_top_n_test_sorted = list(np.transpose(list_feature_selection_top_n_test_sorted))

            # output for feature selection
            file_feature_selection_output = file_output_first_dataset.loc[file_output_first_dataset['GEO asscession number'].isin(feature_selection_all)]

            # sorting data according to its order in testing data
            list_sample_to_read = list(file_feature_selection_top_n_test_sorted.columns.values)

            file_feature_selection_output['sample_id'] = file_feature_selection_output['GEO asscession number'].apply(lambda name: list_sample_to_read.index(name))
            file_feature_selection_output = file_feature_selection_output.sort_values(by = ['sample_id'])
            file_feature_selection_output.drop(columns = 'sample_id', inplace = True)

            # create list of output
            list_desired_output_feature_selection = []
            for element in file_feature_selection_output.loc[:, 'relapse (1=True)']:
                list_desired_output_feature_selection.append(element)
            
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

            print(" # Process : Sequential Forward Selection (SFS)")
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

                        input_no_relapse = []
                        for sample_index in range(0, len(list_top_n_gene_no_relapse_sorted)):
                            list_each_sample = []
                            for element_id in range(0, len(list_top_n_gene_no_relapse_sorted[sample_index])):
                                if (element_id in gene_order_test):
                                    list_each_sample.append(list_top_n_gene_no_relapse_sorted[sample_index][element_id])
                            input_no_relapse.append(list_each_sample)

                        input_testing_data = []
                        for sample_index in range(0, len(list_feature_selection_top_n_test_sorted)):
                            list_each_sample = []
                            for element_id in range(0, len(list_feature_selection_top_n_test_sorted[sample_index])):
                                if (element_id in gene_order_test):
                                    list_each_sample.append(list_feature_selection_top_n_test_sorted[sample_index][element_id])
                            input_testing_data.append(list_each_sample)

                        list_actual_output = calculate.lda(input_testing_data, input_relapse, input_no_relapse)

                        # calculate AUC score
                        auc_score = roc_auc_score(list_desired_output_feature_selection, list_actual_output)
                        # print("auc_score = " + str(auc_score))

                        if (auc_score > max_auc_score):
                            max_auc_score = auc_score
                            gene_index_in_list = i

                            if max_auc_score not in list_auc:
                                list_auc.append(max_auc_score)

                    # do not add gene that already exists in a feature
                    if (gene_index_in_list not in gene_order):
                        gene_order.extend([gene_index_in_list])                       
                    count_iteration += 1  
   
            gene_order.sort()      

            # get gene_name
            gene_order_name = []
            for element in gene_order:
                gene_order_name.append(top_n_genes_name[element])

            # copy required data to be used in evaluation
            top_n_genes_name_for_eval = deepcopy(top_n_genes_name)
            feature_set  = deepcopy(gene_order)
            feature_set_name = deepcopy(gene_order_name) 
        
        # conducting cross-validation on the second dataset
        # prepare data for cross-validation
        print()
        print(" # Conducting the cross-validation on the second dataset")
        print()
        print(" Process : Cross validation ...")
        # split data into k parts
        chunk_relapse_size_cv = math.ceil(len(list_sample_relapse_second_dataset) / num_of_folds)
        chunk_no_relapse_size_cv = math.ceil(len(list_sample_no_relapse_second_dataset) / num_of_folds)

        chunk_list_relapse_cv = list(calculate.chunks(list_sample_relapse_second_dataset, chunk_relapse_size_cv))
        print(" Number of chunks in class relapse : " + str(len(chunk_list_relapse_cv)))

        chunk_list_no_relapse_cv = list(calculate.chunks(list_sample_no_relapse_second_dataset, chunk_no_relapse_size_cv))
        print(" Number of chunks in class non-relapse : " + str(len(chunk_list_no_relapse_cv)))

        check_valid_cv, num_of_chunks_cv = calculate.checkEqualListSize(chunk_list_relapse_cv, chunk_list_no_relapse_cv)
        
        # list and variable to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        # do only if number of chunks of both datasets are equal
        if (check_valid_cv == True):
            for chunk_test_index in range(0, num_of_chunks_cv):
                start_fold_time = time.time()
                result_file.write("\n#### Fold " + str(chunk_test_index + 1) + " ####\n")

                # separating data into testing and training dataset
                # get testing set in this fold
                chunk_test_relapse = chunk_list_relapse_cv[chunk_test_index]
                chunk_test_no_relapse = chunk_list_no_relapse_cv[chunk_test_index]

                print("\n------------------------------------------ K : " + str(chunk_test_index + 1) + " --------------------------------")
                print(" Samples in class relapse used as testing data : " + str(chunk_test_relapse))
                print()
                print(" Samples in class non-relapse used as testing data : " + str(chunk_test_no_relapse))
                print()

                # get training set in this fold
                # for class "relapse"
                chunk_train_relapse = []
                for chunk_train_relapse_index in range(0, num_of_chunks_cv):
                    if (chunk_list_relapse_cv[chunk_train_relapse_index] is not chunk_test_relapse):
                        chunk_train_relapse.append(chunk_list_relapse_cv[chunk_train_relapse_index])

                # for class "non-relapse"
                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks_cv):
                    if (chunk_list_no_relapse_cv[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse_cv[chunk_train_no_relapse_index])

                # merge training data of each class
                # for class "relapse"
                list_train_relapse = []
                for i in range(0, len(chunk_train_relapse)):
                    list_train_relapse.extend(chunk_train_relapse[i])
                print(" Samples in class relapse used as training set " + "(" + str(len(list_train_relapse)) + " samples) : ")
                print(list_train_relapse)
                print()

                # for class "non-relapse"
                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])
                print(" Samplse in class non-relapse used as training set (" + str(len(list_train_no_relapse)) + " samples) : ")
                print(list_train_no_relapse)
                print()

                # preparing data for evaluation and creating classifier
                # for class 'relapse'
                col_to_read_relapse_for_eval = ["ID_REF"]
                col_to_read_relapse_for_eval.extend(list_train_relapse)
                file_training_input_relapse_for_eval = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset, usecols = col_to_read_relapse_for_eval)
                
                top_n_genes_relapse_for_eval = file_training_input_relapse_for_eval.loc[file_training_input_relapse_for_eval['ID_REF'].isin(feature_set_name)]
                top_n_genes_relapse_for_eval['gene_id'] = top_n_genes_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))
                
                top_n_genes_relapse_sorted_for_eval  = top_n_genes_relapse_for_eval.sort_values(by = ['gene_id'])
                top_n_genes_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)

                top_n_genes_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

                # for class 'non-relapse'
                col_to_read_no_relapse_for_eval = ["ID_REF"]
                col_to_read_no_relapse_for_eval.extend(list_train_no_relapse)
                file_training_input_no_relapse_for_eval = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset, usecols = col_to_read_no_relapse_for_eval)

                top_n_genes_no_relapse_for_eval = file_training_input_no_relapse_for_eval.loc[file_training_input_no_relapse_for_eval['ID_REF'].isin(feature_set_name)]
                top_n_genes_no_relapse_for_eval['gene_id'] = top_n_genes_no_relapse_for_eval['ID_REF'].apply(lambda name: feature_set_name.index(name))

                top_n_genes_no_relapse_sorted_for_eval  = top_n_genes_no_relapse_for_eval.sort_values(by = ['gene_id'])
                top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'gene_id', inplace = True)

                top_n_genes_no_relapse_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

                # get testing data
                list_testing_all = []
                list_testing_all.extend(chunk_test_relapse)
                list_testing_all.extend(chunk_test_no_relapse)  

                col_to_read_file_testing = ["ID_REF"]
                col_to_read_file_testing.extend(list_testing_all)

                file_testing = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset, usecols = col_to_read_file_testing)

                file_testing_top_n_test = file_testing.loc[file_testing['ID_REF'].isin(feature_set_name)]
                file_testing_top_n_test['gene_id'] = file_testing_top_n_test['ID_REF'].apply(lambda name: feature_set_name.index(name))   

                file_testing_top_n_test_sorted = file_testing_top_n_test.sort_values(by = ['gene_id'])
                file_testing_top_n_test_sorted.drop(columns = 'gene_id', inplace = True)

                top_n_test_sorted_for_eval = file_testing_top_n_test_sorted
                top_n_test_sorted_for_eval.drop(columns = 'ID_REF', inplace = True)

                # prepare list for testing
                list_top_n_test_sorted = []
                for column in range(0, len(top_n_test_sorted_for_eval)):
                    list_each_sample = []
                    for element in top_n_test_sorted_for_eval.iloc[column]:
                        list_each_sample.append(element)
                    list_top_n_test_sorted.append(list_each_sample)
                list_top_n_test_sorted = list(np.transpose(list_top_n_test_sorted)) 

                # output for testing 
                file_testing_output = file_output_second_dataset.loc[file_output_second_dataset['GEO asscession number'].isin(list_testing_all)]
                list_sample_to_read_for_eval = list(file_testing_top_n_test_sorted.columns.values)

                file_testing_output['sample_id'] = file_testing_output['GEO asscession number'].apply(lambda name: list_sample_to_read_for_eval.index(name))
                file_testing_output = file_testing_output.sort_values(by = ['sample_id'])
                file_testing_output.drop(columns = 'sample_id', inplace = True)

                # create list of output
                list_desired_output_for_eval = []
                for element in file_testing_output.loc[:, 'relapse (1=True)']:
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
                for sample_index in range(0, len(list_top_n_test_sorted)):
                    list_each_sample = []
                    for element_id in range(0, len(list_top_n_test_sorted[sample_index])):
                        if (element_id in feature_set):
                            list_each_sample.append(list_top_n_test_sorted[sample_index][element_id])
                    input_testing_data_for_eval.append(list_each_sample)
                
                list_actual_output_for_eval = calculate.lda(input_testing_data_for_eval, input_relapse_for_eval, input_no_relapse_for_eval)

                # calculate AUC score
                auc_score_for_eval = roc_auc_score(list_desired_output_for_eval, list_actual_output_for_eval)
                list_auc_score.append(auc_score_for_eval)

                print("#### Evaluation of " + str(chunk_test_index + 1) + " - fold ####")
                print(" Feature Set : " + str(feature_set_name))
                print(" Actual Output : " + str(list_actual_output_for_eval))
                print(" Desired Output : " + str(list_desired_output_for_eval))
                print(" AUC ROC score from testing : " + str(auc_score_for_eval))

                # track feature set which gives maximum auc score
                if (auc_score_for_eval > auc_score_max):
                    list_feature_set_max_auc = deepcopy(feature_set_name)
                    auc_score_max = auc_score_for_eval

                # record ending time of this iteration
                end_epoch_time = time.time()
                time_elapse_epoch_second = end_epoch_time - start_epoch_time
                time_elapse_epoch_minute = time_elapse_epoch_second / 60
                time_elapse_epoch_hour = time_elapse_epoch_minute / 60

                # write output to an output file
                result_file.write("Fold : " + str(chunk_test_index + 1) + "\n")
                result_file.write("Feature Set : " + str(feature_set_name) + "\n")
                result_file.write("Actual Output : " + str(list_actual_output_for_eval) + "\n")
                result_file.write("Desired Output : " + str(list_desired_output_for_eval) + "\n")
                result_file.write("AUC ROC Score from testing : " + str(auc_score_for_eval) +  "\n")
                result_file.write("\n")

                end_fold_time = time.time()
                fold_elapse_time_second = end_fold_time - start_fold_time
                fold_elapse_time_minute = fold_elapse_time_second / 60
                fold_elapse_time_minute = round(fold_elapse_time_minute, 2)
                print(" Fold elapse time : " + str(fold_elapse_time_minute) + " minutes")
                result_file.write("Fold elapse time : " + str(fold_elapse_time_minute) + " minutes \n")
                result_file.write("\n")

        list_avg_auc_each_epoch.append(calculate.mean(list_auc_score))

        # record ending time of this iteration
        end_epoch_time = time.time()
        time_elapse_epoch_second = end_epoch_time - start_epoch_time
        time_elapse_epoch_minute = time_elapse_epoch_second / 60
        time_elapse_epoch_hour = time_elapse_epoch_minute / 60

        time_elapse_epoch_minute = round(time_elapse_epoch_minute, 2)
        time_elapse_epoch_hour = round(time_elapse_epoch_hour, 2)

        print()
        print("#### Summary ####")
        print(" Average AUC score : " + str(calculate.mean(list_auc_score)))
        print(" Maximum AUC score : " + str(auc_score_max))
        print(" Feature set which gives highest AUC score : ")
        print(list_feature_set_max_auc)
        print()
        print(" Time Elapse : " + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours)\n")

        result_file.write("\n#### Summary ####\n")
        result_file.write("Average AUC score : " + str(calculate.mean(list_auc_score)) + "\n")
        result_file.write("Maximum AUC score : " + str(auc_score_max) + "\n")
        result_file.write("Size of feature set which gives the highest AUC score from testing : " + str(len(list_feature_set_max_auc)))
        result_file.write("\n")
        result_file.write("Feature set which gives the highest AUC score from testing : " + "\n")
        result_file.write(str(list_feature_set_max_auc))
        result_file.write("\n")
        result_file.write("Time Elapse : " + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours)\n")
        result_file.write("\n")

    end_time = time.time()
    total_elapse_time_second = end_time - start_time

    total_elapse_time_minute = (total_elapse_time_second / 60)
    total_elapse_time_hour = (total_elapse_time_minute / 60)  

    total_elapse_time_minute = round(total_elapse_time_minute, 2)     
    total_elapse_time_hour = round(total_elapse_time_hour, 2)

    # calculate mean over all epoch
    mean_over_all_epoch = calculate.mean(list_avg_auc_each_epoch)
    print(" Average AUC score over " + str(epoch) + " epoch : " + str(mean_over_all_epoch))
    result_file.write("Average AUC score over " + str(epoch) + " epoch : " + str(mean_over_all_epoch) + "\n")

    print(" Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("\n")

    result_file.close()

if __name__ == "__main__":
    main()