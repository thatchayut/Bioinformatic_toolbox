import numpy as np 
import pandas as pd
import random
import math
import calculate
import time

def main():
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

    # calculate number of pathways to be used
    num_of_ranked_pathways = (rows_to_read_file_pathway * (num_of_pathways_percentage / 100))
    num_of_ranked_pathways = math.ceil(num_of_ranked_pathways)

    # create list of all samples in trainig data in this fold
    list_all_samples = []
    list_all_samples.extend(list_sample_relapse)
    list_all_samples.extend(list_sample_no_relapse)

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        
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
                # genes_expression = []
                for line_index in range(0, row_to_read_file_to_cal_zscore):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_cal_zscore.iloc[line_index, 1:-1]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_cal_zscore.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression[line_index] = gene_expression_by_probe_id
                    # genes_expression.append(gene_expression_by_probe_id)
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

                samples_relapse = {}
                samples_no_relapse = {}

                for element_index in range(0, len(list_sample_relapse)):
                    print()
                    print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
                    print(str(len(list_sample_relapse) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_sample_relapse[element_index]
                    pathways = calculate.vgetPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                                list_mean_sd_gene_expression_by_probe_id)
                   
                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_relapse[element_index] = sample   



if __name__ == "__main__":
    main()