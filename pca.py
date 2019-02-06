from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
import xlsxwriter
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main():
     # record start time
    start_time = time.time()

    # prepare data
    # row_to_read = 22283
    row_to_read = 100
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    rows_to_read_file_pathway = 100
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    
    # get list of pathway name
    list_pathway_name = []
    for i in range(0, rows_to_read_file_pathway):
        pathway_name = []
        pathway_name.append(i)
        pathway_name.append(file_pathway.loc[i, "PATHWAY_NAME"])
        list_pathway_name.append(pathway_name)

    # consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['0'])]

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
            print("WARNING : Input must be numeric")
        elif(int(num_of_folds) > len(list_sample_relapse)):
            print("WARNING : Number of folds exceeds the size of the 1st dataset")
        elif(int(num_of_folds) > len(list_sample_no_relapse)):
            print("WARNING : Number of folds exceeds the size of the 2nd dataset")
        elif(int(num_of_folds) <= 1):
            print("WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)

    # get number of epochs
    while True:
        num_of_epochs = input("Number of epochs: ")
        if (num_of_epochs.isnumeric() == False):
            print("WARNING : Input must be numeric")
        elif(int(num_of_epochs) < 1):
            print("WARNING : Number of folds cannot lower than 1")
        else:
            break
    num_of_epochs = int(num_of_epochs)

    # get number of pathways
    while True:
        num_of_pathways_percentage = input("number of pathways to be used as feature set (%) : ")
        if (num_of_pathways_percentage.isnumeric() == False):
            print("WARNING : Input must be numeric")
        elif (int(num_of_pathways_percentage) <= 0):
            print("WARNING : Percentage must greater than 0.")
        elif (int(num_of_pathways_percentage) > 100):
            print("WARNING : Percentage must lower than or equal to 100.")
        else:
            break
    num_of_pathways_percentage = int(num_of_pathways_percentage)

    # ask for scaling method
    print("\n Scaling methods: [1] z-score // [2] narrow scaling (range [0,1]) // [3] wide scaling (range [-1,1])") 
    method_id = None
    while True:
        method_id = input(" Enter method id : ")
        if (method_id not in ["1", "2", "3"]):
            print(" WARNING : Invalid method id")
        else:
            break

    # get output file's name
    # file_name = input("Name of Output File : ")

    # # prepare text file for results to be written in
    # result_file = open(str(file_name) + ".txt", "w+")

    # calculate number of pathways to be used
    num_of_ranked_pathways = (rows_to_read_file_pathway * (num_of_pathways_percentage / 100))
    num_of_ranked_pathways = math.ceil(num_of_ranked_pathways)

    # create list of all samples in trainig data in this fold
    list_all_samples = []
    list_all_samples.extend(list_sample_relapse)
    list_all_samples.extend(list_sample_no_relapse)

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        # result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        
        # create list of indexes used to indicate the position in the list
        list_index_samples_relapse = []
        list_index_samples_no_relapse = []

        for index in range(0, len(list_sample_relapse)):
            list_index_samples_relapse.append(index)
        
        for index in range(0, len(list_sample_no_relapse)):
            list_index_samples_no_relapse.append(index)
        
        # shuffle it to make it flexible for epoch changed
        random.shuffle(list_index_samples_relapse)
        random.shuffle(list_index_samples_no_relapse)
        
        print()
        print("list_index_samples_relapse size : " + str(len(list_index_samples_relapse)))
        print("list_index_samples_relapse : ")
        print(list_index_samples_relapse)
        print()
        print("list_index_samples_no_relapse size : " + str(len(list_index_samples_no_relapse)))
        print("list_index_samples_no_relapse : ")
        print(list_index_samples_no_relapse)

        # split data into k parts
        chunk_relapse_size = math.ceil(len(list_index_samples_relapse) / num_of_folds)
        chunk_no_relapse_size = math.ceil(len(list_index_samples_no_relapse) / num_of_folds)

        chunk_list_relapse = list(calculate.chunks(list_index_samples_relapse, chunk_relapse_size))
        print("number of chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))
        print("chunk_list_relapse : ")
        print(chunk_list_relapse)
        print()

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse, chunk_no_relapse_size))
        print("number of in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))
        print("chunk_list_no_relapse : ")
        print(chunk_list_no_relapse)

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # list to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        # list to collect maximun AUC in each fold
        list_max_auc = []

        # do only if number of chunks of both datasets are equal
        if (check_valid == True):
            for chunk_test_index in range(0, num_of_chunks):
                feature_set = []
                feature_set_name = []

                strat_fold_time = time.time()
                # result_file.write("\n#### Fold " + str(chunk_test_index + 1) + " ####\n")

                # separating data into testing and training dataset
                # get testing set in this fold
                chunk_test_relapse = chunk_list_relapse[chunk_test_index]
                chunk_test_no_relapse = chunk_list_no_relapse[chunk_test_index]

                print("\n------------------------------------------ K : " + str(chunk_test_index + 1) + " --------------------------------")
                print("chunk_test_relapse : " + str(chunk_test_relapse))
                print()
                print("chunk_test_no_relapse : " + str(chunk_test_no_relapse))
                print()

                # get training set in this fold
                chunk_train_relapse = []
                for chunk_train_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_relapse[chunk_train_relapse_index] is not chunk_test_relapse):
                        chunk_train_relapse.append(chunk_list_relapse[chunk_train_relapse_index])
                print("chunk train relapse size = " + str(len(chunk_train_relapse)))
                # print("chunk train relapse = " + str(chunk_train_relapse))
                # print()
                
                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse[chunk_train_no_relapse_index])
                print("chunk train no relapse size = " + str(len(chunk_train_no_relapse)))
                # print("chunk train no relapse = " + str(chunk_train_no_relapse))
                # print()

                # merge training data of each class
                list_train_relapse = []
                for i in range(0, len(chunk_train_relapse)):
                    list_train_relapse.extend(chunk_train_relapse[i])
                print("size of list_train_relapse : " + str(len(list_train_relapse)))
                print("list_train_relapse : ")
                print(list_train_relapse)
                print()

                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])
                print("size of list_train_no_relapse : " + str(len(list_train_no_relapse)))
                print("list_train_no_relapse : ")
                print(list_train_no_relapse)
                print()


                print("------------------------------------------------------------------------------------------------------------")
                list_train_relapse_name = []
                list_train_no_relapse_name = []

                # get sample name and add to a list to be used as column index
                for i in range(0, len(list_train_relapse)):
                    list_train_relapse_name.append(list_sample_relapse[list_train_relapse[i]])
                
                for i in range(0, len(list_train_no_relapse)):
                    list_train_no_relapse_name.append(list_sample_no_relapse[list_train_no_relapse[i]])
                
                # NORMALIZE HERE !!!
                # create list to collect mean of each gene 
                # calculate mean and sd directly from a file
                list_gene_expression_by_entrez = []
                list_gene_name_probe_id = []
                # row_to_read_file_to_cal_zscore = 22283
                row_to_read_file_to_cal_zscore = 100
                col_to_read_file_to_cal_zscore = ["ID_REF"]
                col_to_read_file_to_cal_zscore.extend(list_train_relapse_name)
                col_to_read_file_to_cal_zscore.extend(list_train_no_relapse_name)
                print("col_to_read_file_to_cal_zscore : ")
                print(col_to_read_file_to_cal_zscore)
                file_to_cal_zscore = pd.read_csv("GSE2034-22071 (edited).csv", usecols = col_to_read_file_to_cal_zscore, nrows = row_to_read)
                # num_of_all_genes = len(samples_relapse[0][1])
                num_of_all_samples  = len(col_to_read_file_to_cal_zscore)
                print("------------------------------------------------------------------------------------------------------------")
                # print("num_of_all_genes : " + str(num_of_all_genes))

                # dictionary contains genes idintified by its probe id which contain all gene expression from all samples
                # {1: [gene_probe_id, [exp1,exp2, ...]]}
                genes_expression = {}
                for line_index in range(0, row_to_read_file_to_cal_zscore):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_cal_zscore.iloc[line_index, 1:-1]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_cal_zscore.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression[line_index] = gene_expression_by_probe_id
                # print("genes_expression : ")
                # print(genes_expression)

                # calculate mean and sd of each gene
                list_mean_sd_gene_expression_by_probe_id = []
                for gene_index in range(0, len(genes_expression)):
                    result = []
                    gene_name = genes_expression[gene_index][0]
                    mean_of_list = calculate.mean(genes_expression[gene_index][1])
                    sd_of_list = calculate.sd(genes_expression[gene_index][1])
                    result.append(gene_name)
                    result.append(mean_of_list)
                    result.append(sd_of_list)
                    list_mean_sd_gene_expression_by_probe_id.append(result)
                # print()
                # print("list_mean_sd_gene_expression_by_probe_id : ")
                # print(list_mean_sd_gene_expression_by_probe_id)
                # print()

                # create samples
                samples_relapse = {}
                samples_no_relapse = {}


                for element_index in range(0, len(list_sample_relapse)):
                # for element_index in range(0, 6):
                    print()
                    print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
                    print(str(len(list_sample_relapse) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_sample_relapse[element_index]
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                                list_mean_sd_gene_expression_by_probe_id)

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_relapse[element_index] = sample       
                # print("samples_relapse : ")
                # print(samples_relapse)

                for element_index in range(0, len(list_sample_no_relapse)):
                # for element_index in range(0, 6):
                    print()
                    print("Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
                    print(str(len(list_sample_no_relapse) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_sample_no_relapse[element_index]
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, \
                                list_mean_sd_gene_expression_by_probe_id)

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_no_relapse[element_index] = sample 

                print("samples_relapse : ")
                print(samples_relapse)
                print()    
                print("samples_no_relapse : ")
                print(samples_no_relapse)
                print()

                # splitting lists to use them as an evaluation set and feature selection set
                # given that we use 3-fold cross validation in this level
                print("\n#### given that we use 3-fold cross validation in this level for evaluation and feature selection ####")
                second_num_of_fold = 3
                second_chunk_relapse_size = math.ceil(len(list_train_relapse) / second_num_of_fold)
                second_chunk_no_relapse_size = math.ceil(len(list_train_no_relapse) / second_num_of_fold)

                second_chunk_list_relapse = list(calculate.chunks(list_train_relapse, second_chunk_relapse_size))
                print("# chunks in second_chunk_list_relapse = " + str(len(second_chunk_list_relapse)))
                second_chunk_list_no_relapse = list(calculate.chunks(list_train_no_relapse, second_chunk_no_relapse_size))
                print("# chunks in second_chunk_list_no_relapse = " + str(len(second_chunk_list_no_relapse)))

                second_check_valid, second_num_of_chunks = calculate.checkEqualListSize(second_chunk_list_relapse, second_chunk_list_no_relapse)

                if (second_check_valid is True):
                    # random a chunk of data to be use as a feature selection set
                    second_layer_test_index = random.randint(0, second_num_of_chunks - 1)

                    # get testing data from each class
                    second_layer_test_relapse =  second_chunk_list_relapse[second_layer_test_index]
                    second_layer_test_no_relapse = second_chunk_list_no_relapse[second_layer_test_index]

                    # separate training dataset from testing dataset to use in t-test ranking
                    second_layer_train_relapse = []
                    for second_layer_train_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_relapse[second_layer_train_index] is not second_layer_test_relapse):
                            second_layer_train_relapse.append(second_chunk_list_relapse[second_layer_train_index])
                    print("second_layer_train_relapse size = " + str(len(second_layer_train_relapse)))
                    print("second_layer_train_relapse = " + str(second_layer_train_relapse))
                    print()
                    
                    second_layer_train_no_relapse = []
                    for second_layer_train_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_no_relapse[second_layer_train_index] is not second_layer_test_no_relapse):
                            second_layer_train_no_relapse.append(second_chunk_list_no_relapse[second_layer_train_index])
                    print("second_layer_train_no_relapse size : " + str(len(second_layer_train_no_relapse)))
                    print("second_layer_train_no_relapse : " + str(second_layer_train_no_relapse))    
                    print()

                    # merge all samples in the same class
                    print("\n#### merge all samples in the same class to be used further ####")
                    list_sample_relapse_train_feature_selection = []
                    # ttest_list_sample_relapse = []
                    for i in range(0, len(second_layer_train_relapse)):
                        list_sample_relapse_train_feature_selection.extend(second_layer_train_relapse[i])
                    # print("size of ttest list sample relapse = " + str(len(list_sample_relapse_train_feature_selection)))
                    print("list_sample_relapse_train_feature_selection : " + str(list_sample_relapse_train_feature_selection))
                    
                    list_sample_no_relapse_train_feature_selection = []
                    # ttest_list_sample_no_relapse = []
                    for i in range(0, len(second_layer_train_no_relapse)):
                        list_sample_no_relapse_train_feature_selection.extend(second_layer_train_no_relapse[i])
                    # print("size of ttest list sample no relapse = " + str(len(list_sample_no_relapse_train_feature_selection)))
                    print("list_sample_no_relapse_train_feature_selection : " + str(list_sample_no_relapse_train_feature_selection))

                    # create collection of samples used in feature selection
                    samples_relapse_feature_selection = {}
                    samples_no_relapse_feature_selection = {}

                    for sample_index in range(0, len(list_sample_relapse_train_feature_selection)):
                        index_samples_relapse = list_sample_relapse_train_feature_selection[sample_index]
                        samples_relapse_feature_selection[sample_index] = samples_relapse[index_samples_relapse]
                    print()
                    print("samples_relapse_feature_selection : ")
                    print(samples_relapse_feature_selection)
                    print()

                    for sample_index in range(0, len(list_sample_no_relapse_train_feature_selection)):
                        index_samples_no_relapse = list_sample_no_relapse_train_feature_selection[sample_index]
                        samples_no_relapse_feature_selection[sample_index] = samples_no_relapse[index_samples_no_relapse]
                    print()
                    print("samples_no_relapse_feature_selection : ")
                    print(samples_no_relapse_feature_selection)
                    print()

                    # find CORG
                    #----> START HERE
                    list_corg_each_pathway = []

                    for pathway_index in range(0, rows_to_read_file_pathway):
                        list_ttest_gene_activity = []

                        list_gene_name_in_pathway = []
                        list_gene_sample_relapse = []
                        list_gene_sample_no_relapse = []

                        num_of_genes_in_pathway = len(samples_relapse_feature_selection[0][1][pathway_index][1])

                        # create list of gene entrez id in this pathway
                        for gene_index in range(0, num_of_genes_in_pathway):
                            gene_entrez_id = samples_relapse_feature_selection[0][1][pathway_index][1][gene_index][0]
                            list_gene_name_in_pathway.append(gene_entrez_id)

                        
                        # create list of gene expression 
                        for gene_index in range(0, len(list_gene_name_in_pathway)):
                            list_gene_expression_with_entrez = []
                            list_gene_expression_from_sample = []
                            gene_entrez_id = list_gene_name_in_pathway[gene_index]
                            for sample_index in range(0, len(samples_relapse_feature_selection)):
                                gene_expression = samples_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                                list_gene_expression_from_sample.append(gene_expression)
                            list_gene_expression_with_entrez.append(gene_entrez_id)    
                            list_gene_expression_with_entrez.append(list_gene_expression_from_sample)
                            list_gene_sample_relapse.append(list_gene_expression_with_entrez)
                        # print("list_gene_sample_relapse : ")
                        # print(list_gene_sample_relapse)

                        for gene_index in range(0, len(list_gene_name_in_pathway)):
                            list_gene_expression_with_entrez = []
                            list_gene_expression_from_sample = []
                            gene_entrez_id = list_gene_name_in_pathway[gene_index]
                            for sample_index in range(0, len(samples_no_relapse_feature_selection)):
                                gene_expression = samples_no_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                                list_gene_expression_from_sample.append(gene_expression)
                            list_gene_expression_with_entrez.append(gene_entrez_id)
                            list_gene_expression_with_entrez.append(list_gene_expression_from_sample)    
                            list_gene_sample_no_relapse.append(list_gene_expression_with_entrez)
                    
                        # find t-test score of each gene
                        list_ttest = []
                        for gene_index in range(0, num_of_genes_in_pathway):
                            gene_ttest = []
                            gene_name = list_gene_name_in_pathway[gene_index]
                            gene_ttest_value = stats.ttest_ind(list_gene_sample_relapse[gene_index][1], list_gene_sample_no_relapse[gene_index][1], equal_var = False)[0]
                            # gene_ttest_value = gene_index
                            gene_ttest.append(gene_name)
                            gene_ttest.append(gene_ttest_value)
                            list_ttest.append(gene_ttest)

                        list_ttest.sort(key = lambda x : x[1], reverse = True)
                        
                        # create corg of each pathway (use top 10% as an initial set)
                        num_of_top_gene = int((num_of_genes_in_pathway * 10) / 100)
                        if (num_of_top_gene < 1):
                            num_of_top_gene = 1

                        list_corg_initial = []
                        for i in range(0, num_of_top_gene):
                            list_corg_initial.append(list_ttest[i][0])

                        check_finish = False
                        check_improve_discriminative = True
                        max_discriminative_over_all_features = 0
                        # list_corg = []

                        while (check_finish == False):
                            if (check_improve_discriminative == True):
                                max_ttest_in_consider = 0
                                list_member_gene = []
                                for gene_index in range(0, num_of_genes_in_pathway):
                                    list_gene_to_consider = deepcopy(list_corg_initial)
                                    gene_entrez_id = list_ttest[gene_index][0]
                                    
                                    if (gene_entrez_id not in list_corg_initial):
                                        list_gene_to_consider.extend([gene_entrez_id])
                                    # print("list_gene_to_consider : " + str(list_gene_to_consider))    

                                    # create list of gene expression of each sample using this member genes
                                    list_sample_relapse_find_discrimination = []
                                    list_sample_no_relapse_find_discrimination = []

                                    if (len(list_gene_to_consider) == 1):
                                        gene_entrez_id = list_gene_to_consider[0]
                                        # find gene expression using gene entrez id
                                        for gene_sample_relapse_index in range(0, len(list_gene_sample_relapse)):
                                            if (list_gene_sample_relapse[gene_sample_relapse_index][0] == gene_entrez_id):
                                                list_sample_relapse_find_discrimination.append(list_gene_sample_relapse[gene_sample_relapse_index][1])
                                    else:
                                        for member_gene_index in range(0, len(list_gene_to_consider)):
                                            gene_entrez_id = list_gene_to_consider[member_gene_index]
                                            # find gene expression using gene entrez id
                                            for gene_sample_relapse_index in range(0, len(list_gene_sample_relapse)):
                                                if (list_gene_sample_relapse[gene_sample_relapse_index][0] == gene_entrez_id):
                                                    list_sample_relapse_find_discrimination.append(list_gene_sample_relapse[gene_sample_relapse_index][1])                              
                                    # print("list_sample_relapse_find_discrimination ")
                                    # print(list_sample_relapse_find_discrimination)
                                    # print()


                                    if (len(list_gene_to_consider) == 1):
                                        gene_entrez_id = list_gene_to_consider[0]
                                        # find gene expression using gene entrez id
                                        for gene_sample_no_relapse_index in range(0, len(list_gene_sample_no_relapse)):
                                            if (list_gene_sample_no_relapse[gene_sample_no_relapse_index][0] == gene_entrez_id):
                                                list_sample_no_relapse_find_discrimination.append(list_gene_sample_no_relapse[gene_sample_no_relapse_index][1])
                                    else:
                                        for member_gene_index in range(0, len(list_gene_to_consider)):
                                            gene_entrez_id = list_gene_to_consider[member_gene_index]
                                            # find gene expression using gene entrez id
                                            for gene_sample_no_relapse_index in range(0, len(list_gene_sample_no_relapse)):
                                                if (list_gene_sample_no_relapse[gene_sample_no_relapse_index][0] == gene_entrez_id):
                                                    list_sample_no_relapse_find_discrimination.append(list_gene_sample_no_relapse[gene_sample_no_relapse_index][1])
                                    # print("list_sample_no_relapse_find_discrimination ")
                                    # print(list_sample_no_relapse_find_discrimination)

                                    # calculate activity score
                                    list_sample_relapse_activity_score = []
                                    list_sample_no_relapse_activity_score = []

                                    for sample_index in range(0, len(samples_relapse_feature_selection)):
                                        sum_gene_expression = 0
                                        for gene_index in range(0, len(list_sample_relapse_find_discrimination)):
                                            sum_gene_expression += list_sample_relapse_find_discrimination[gene_index][sample_index]
                                        activity_score = (sum_gene_expression / math.sqrt(len(list_sample_relapse_find_discrimination)))
                                        list_sample_relapse_activity_score.append(activity_score)
                                    # print("list_sample_relapse_activity_score : " + str(list_sample_relapse_activity_score))

                                    for sample_index in range(0, len(samples_no_relapse_feature_selection)):
                                        sum_gene_expression = 0
                                        for gene_index in range(0, len(list_sample_no_relapse_find_discrimination)):
                                            sum_gene_expression += list_sample_no_relapse_find_discrimination[gene_index][sample_index]
                                        activity_score = (sum_gene_expression / math.sqrt(len(list_sample_no_relapse_find_discrimination)))
                                        list_sample_no_relapse_activity_score.append(activity_score)
                                    # print("list_sample_no_relapse_activity_score : " + str(list_sample_no_relapse_activity_score))
                                    # print()

                                    # calculate ttest score of this gene set as a discriminative score
                                    ttest_member_gene_set = stats.ttest_ind(list_sample_relapse_activity_score, list_sample_no_relapse_activity_score, equal_var = False)[0]
                                    # print("ttest_member_gene_set : " + str(ttest_member_gene_set))

                                    if (ttest_member_gene_set > max_ttest_in_consider):
                                        max_ttest_in_consider = ttest_member_gene_set
                                        list_member_gene = deepcopy(list_gene_to_consider)
                                        
                                if (max_ttest_in_consider > max_discriminative_over_all_features):
                                    max_discriminative_over_all_features = max_ttest_in_consider
                                    list_corg_initial = deepcopy(list_member_gene)
                                else:
                                    check_improve_discriminative = False
                            else:
                                check_finish = True          
                        list_corg_each_pathway.append(list_corg_initial)                     
                        # <----- END HERE
                    print("list_corg_each_pathway : ")
                    print(list_corg_each_pathway)
                    print()

                    # create samples for validation
                    samples_relapse_validation = {}
                    samples_no_relapse_validation = {}

                    for sample_index in range(0, len(second_layer_test_relapse)):
                        index_samples_relapse = second_layer_test_relapse[sample_index]
                        samples_relapse_validation[sample_index] = samples_relapse[index_samples_relapse]

                    
                    for sample_index in range(0, len(second_layer_test_no_relapse)):
                        index_samples_no_relapse = second_layer_test_no_relapse[sample_index]
                        samples_no_relapse_validation[sample_index] = samples_no_relapse[index_samples_no_relapse]
                    
                    print("samples_relapse_validation[0] :")
                    print(samples_relapse_validation[0])
                    
                    # calculate pathway activity 
                    samples_relapse_validation_pathway_activity = {}
                    samples_no_relapse_validation_pathway_activity = {}

                    for sample_index in range(0, len(samples_relapse_validation)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_relapse_validation[sample_index][1])):
                            pathway = []
                            sum_gene_expression = 0
                            num_of_corg = len(samples_relapse_validation[sample_index][1][pathway_index])
                            for gene_index in range(0, len(samples_relapse_validation[sample_index][1][pathway_index])):
                                gene_entrez_id = samples_relapse_validation[sample_index][1][pathway_index][1][gene_index][0]
                                gene_expression =  samples_relapse_validation[sample_index][1][pathway_index][1][gene_index][1]
                                if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                                    sum_gene_expression += gene_expression
                            pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                            pathway_name = samples_relapse_validation[sample_index][1][pathway_index][0]
                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_relapse_validation[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)
                        
                        samples_relapse_validation_pathway_activity[sample_index] = list_sample_with_pathway_activity
                    
                    print()
                    print("samples_relapse_validation_pathway_activity : ")
                    print(samples_relapse_validation_pathway_activity)
                    

                    for sample_index in range(0, len(samples_no_relapse_validation)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_no_relapse_validation[sample_index][1])):
                            pathway = []
                            sum_gene_expression = 0
                            num_of_corg = len(samples_no_relapse_validation[sample_index][1][pathway_index])
                            for gene_index in range(0, len(samples_no_relapse_validation[sample_index][1][pathway_index])):
                                gene_entrez_id = samples_no_relapse_validation[sample_index][1][pathway_index][1][gene_index][0]
                                gene_expression =  samples_no_relapse_validation[sample_index][1][pathway_index][1][gene_index][1]
                                if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                                    sum_gene_expression += gene_expression
                            pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                            pathway_name = samples_no_relapse_validation[sample_index][1][pathway_index][0]
                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_no_relapse_validation[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)
                        
                        samples_no_relapse_validation_pathway_activity[sample_index] = list_sample_with_pathway_activity
                    
                    print()
                    print("samples_no_relapse_validation_pathway_activity : ")
                    print(samples_no_relapse_validation_pathway_activity)

                    # sort pathways in each sample using p-value in ascending order
                    # create list contains pathway activity of each sample preparing for calculating p-value
                    list_relapse_pathway_activity_for_pvalue = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway = []
                        for sample_index in range(0, len(samples_relapse_validation_pathway_activity)):
                            pathway_activity = samples_relapse_validation_pathway_activity[sample_index][1][pathway_index][1]
                            pathway.append(pathway_activity)
                        list_relapse_pathway_activity_for_pvalue.append(pathway)
                    
                    list_no_relapse_pathway_activity_for_pvalue = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway = []
                        for sample_index in range(0, len(samples_no_relapse_validation_pathway_activity)):
                            pathway_activity = samples_no_relapse_validation_pathway_activity[sample_index][1][pathway_index][1]
                            pathway.append(pathway_activity)
                        list_no_relapse_pathway_activity_for_pvalue.append(pathway)
                    
                    # calculate p-value
                    list_pvalue_pathway_activity = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway_pvalue = []
                        pathway_name = list_pathway_name[pathway_index][1]
                        pvalue = stats.ttest_ind(list_sample_relapse_activity_score, list_sample_no_relapse_activity_score, equal_var = False)[1]
                        pathway_pvalue.append(pathway_name)
                        pathway_pvalue.append(pvalue)
                        list_pvalue_pathway_activity.append(pathway_pvalue)
                    
                    # sort pathway using p-value in ascending order
                    list_pvalue_pathway_activity.sort(key = lambda x : x[1], reverse = False)

                    print("list_pvalue_pathway_activity : " + str(list_pvalue_pathway_activity))

                    # reorder pathway in each sample
                    samples_relapse_validation_pathway_activity_sorted = {}
                    samples_no_relapse_validation_pathway_activity_sorted = {}
                    
                    for sample_index in range(0, len(samples_relapse_validation_pathway_activity)):
                        sample_name = samples_relapse_validation_pathway_activity[sample_index][0]
                        list_relapse_validation_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_relapse_validation_pathway_activity[sample_index][1])):
                                pathway_name = samples_relapse_validation_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_relapse_validation_pathway_activity[sample_index][1][pathway_index][1]
                                
                                # print("pathway_name : " + str(pathway_name))
                                # print("list_pvalue_pathway_activity[pvalue_index] : " + str(list_pvalue_pathway_activity[pvalue_index]))
                                # print()

                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_relapse_validation_pathway_activity_sorted.append(sample_name)
                        list_relapse_validation_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_relapse_validation_pathway_activity_sorted[sample_index] = list_relapse_validation_pathway_activity_sorted
                               
                    for sample_index in range(0, len(samples_no_relapse_validation_pathway_activity)):
                        sample_name = samples_no_relapse_validation_pathway_activity[sample_index][0]
                        list_no_relapse_validation_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_no_relapse_validation_pathway_activity[sample_index][1])):
                                pathway_name = samples_no_relapse_validation_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_no_relapse_validation_pathway_activity[sample_index][1][pathway_index][1]
                                
                                # print("pathway_name : " + str(pathway_name))
                                # print("list_pvalue_pathway_activity[pvalue_index] : " + str(list_pvalue_pathway_activity[pvalue_index]))
                                # print()

                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_no_relapse_validation_pathway_activity_sorted.append(sample_name)
                        list_no_relapse_validation_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_no_relapse_validation_pathway_activity_sorted[sample_index] = list_relapse_validation_pathway_activity_sorted
                    print()
                    print("num of pathway before : " + str(len(samples_no_relapse_validation_pathway_activity[0][1])))
                    print("BEFORE : ")
                    print(samples_no_relapse_validation_pathway_activity[0])
                    print()
                    print("AFTER ")
                    print("num of pathway after : " + str(len(samples_no_relapse_validation_pathway_activity[0][1])))
                    print(samples_no_relapse_validation_pathway_activity_sorted[0])
                    
                    



                    # print("list_gene_name_in_pathway : " + str(list_gene_name_in_pathway))
                    # print("list_gene_name_in_pathway size : " + str(len(list_gene_name_in_pathway)))
                    # print("list_gene_sample_relapse : " + str(list_gene_sample_relapse))
                    # print("list_gene_sample_relapse size : " + str(len(list_gene_sample_relapse)))
                    # print("list_gene_sample_no_relapse : " + str(list_gene_sample_no_relapse))
                    # print("list_gene_sample_no_relapse size : " + str(len(list_gene_sample_no_relapse)))
                    # print("list_ttest : " + str(list_ttest))
                    # print("list_corg : " + str(list_corg))

                    print("\n-------------------------------------------------------------------------------------------------------------\n")







if __name__ == "__main__":
    main()