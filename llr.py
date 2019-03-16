from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
from copy import deepcopy
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

def main():
    # record start time
    start_time = time.time()

    # prepare data
    # default row_to_read = 22283
    row_to_read_file_input = 100
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read_file_input)
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
    for i in range(0, row_to_read_file_input):
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

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)

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
        num_of_pathways_percentage = input("Number of pathways to be used as feature set (%) : ")
        if (num_of_pathways_percentage.isnumeric() == False):
            print("WARNING : Input must be numeric")
        elif (int(num_of_pathways_percentage) <= 0):
            print("WARNING : Percentage must greater than 0.")
        elif (int(num_of_pathways_percentage) > 100):
            print("WARNING : Percentage must lower than or equal to 100.")
        else:
            break
    num_of_pathways_percentage = int(num_of_pathways_percentage)

    # # get output file's name
    # file_name = input("Name of output file : ")

    # # prepare text file for results to be written in
    # result_file = open(str(file_name) + ".txt", "w+")

    # calculate number of pathways to be used
    num_of_ranked_pathways = (rows_to_read_file_pathway * (num_of_pathways_percentage / 100))
    num_of_ranked_pathways = math.ceil(num_of_ranked_pathways)

    # create list of all samples 
    list_all_samples = []
    list_all_samples.extend(list_sample_relapse)
    list_all_samples.extend(list_sample_no_relapse)

    # print("Process : Creating collection to collect samples and their genes' expression")
    # # create dictionary used to collect pathways of each sample
    # samples_relapse = {}
    # samples_no_relapse = {}

    # # get all pathways of all samples in class 'relapse'
    # for element_index in range(0, len(list_sample_relapse)):
    #     print()
    #     print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
    #     print(str(len(list_sample_relapse) - (element_index + 1)) + " samples left")
    #     print()

    #     sample = []
    #     sample_name = list_sample_relapse[element_index]
    #     pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

    #     sample.append(sample_name)
    #     sample.append(pathways)
    #     samples_relapse[element_index] = sample
    
    # for element_index in range(0, len(list_sample_no_relapse)):
    #     print()
    #     print("Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
    #     print(str(len(list_sample_no_relapse) - (element_index + 1)) + " samples left")
    #     print()

    #     sample = []
    #     sample_name = list_sample_no_relapse[element_index]
    #     pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

    #     sample.append(sample_name)
    #     sample.append(pathways)
    #     samples_no_relapse[element_index] = sample
    
    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        # result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################\n")
        
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

        # split data into k parts
        chunk_relapse_size = math.ceil(len(list_index_samples_relapse) / num_of_folds)
        chunk_no_relapse_size = math.ceil(len(list_index_samples_no_relapse) / num_of_folds)

        chunk_list_relapse = list(calculate.chunks(list_index_samples_relapse, chunk_relapse_size))
        print("number of chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse, chunk_no_relapse_size))
        print("number of in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

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
                # list of feature set
                feature_set_name = []

                start_fold_time = time.time()
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

                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse[chunk_train_no_relapse_index])
                print("chunk train no relapse size = " + str(len(chunk_train_no_relapse)))

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

                # get sample name and add to a list to be used as column index
                list_train_relapse_name = []
                list_train_no_relapse_name = []
                for i in range(0, len(list_train_relapse)):
                    list_train_relapse_name.append(list_sample_relapse[list_train_relapse[i]])
                
                for i in range(0, len(list_train_no_relapse)):
                    list_train_no_relapse_name.append(list_sample_no_relapse[list_train_no_relapse[i]])

                # find mean and sd of gene expression of each gene in each class
                # prepare file to get gene expression

                # default row_to_read_file_to_cal_mean_sd = 22283
                row_to_read_file_to_cal_mean_sd = 100

                col_to_read_file_to_cal_mean_sd_relapse = ["ID_REF"]
                col_to_read_file_to_cal_mean_sd_relapse.extend(list_train_relapse_name)

                col_to_read_file_to_cal_mean_sd_no_relapse = ["ID_REF"]
                col_to_read_file_to_cal_mean_sd_no_relapse.extend(list_train_no_relapse_name)

                file_to_cal_mean_sd_relapse = pd.read_csv("GSE2034-22071 (edited).csv", usecols = col_to_read_file_to_cal_mean_sd_relapse, nrows = row_to_read_file_to_cal_mean_sd)
                file_to_cal_mean_sd_no_relapse = pd.read_csv("GSE2034-22071 (edited).csv", usecols = col_to_read_file_to_cal_mean_sd_no_relapse, nrows = row_to_read_file_to_cal_mean_sd)
                
                # dictionary contains genes idintified by its probe id which contain all gene expression from samples
                # {1: [gene_probe_id, [exp1,exp2, ...]]}

                # for class "relapse"
                genes_expression_relapse_train = {}
                last_index_to_read_file_to_cal_mean_sd_relapse = len(col_to_read_file_to_cal_mean_sd_relapse)
                for line_index in range(0, row_to_read_file_to_cal_mean_sd):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_cal_mean_sd_relapse.iloc[line_index, 1:last_index_to_read_file_to_cal_mean_sd_relapse]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_cal_mean_sd_relapse.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression_relapse_train[line_index] = gene_expression_by_probe_id
                
                # for class "non-relapse"
                genes_expression_no_relapse_train = {}
                last_index_to_read_file_to_cal_mean_sd_no_relapse = len(col_to_read_file_to_cal_mean_sd_no_relapse)
                for line_index in range(0, row_to_read_file_to_cal_mean_sd):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_cal_mean_sd_no_relapse.iloc[line_index, 1:last_index_to_read_file_to_cal_mean_sd_no_relapse]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_cal_mean_sd_no_relapse.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression_no_relapse_train[line_index] = gene_expression_by_probe_id

                # for class "relapse"
                list_mean_sd_gene_expression_by_probe_id_relapse = []
                for gene_index in range(0, len(genes_expression_relapse_train)):
                    result = []
                    gene_name = genes_expression_relapse_train[gene_index][0]
                    mean_of_list = calculate.mean(genes_expression_relapse_train[gene_index][1])
                    sd_of_list = calculate.sd(genes_expression_relapse_train[gene_index][1])
                    result.append(gene_name)
                    result.append(mean_of_list)
                    result.append(sd_of_list)
                    list_mean_sd_gene_expression_by_probe_id_relapse.append(result)
                
                # for class "non-relapse"
                list_mean_sd_gene_expression_by_probe_id_no_relapse = []
                for gene_index in range(0, len(genes_expression_no_relapse_train)):
                    result = []
                    gene_name = genes_expression_no_relapse_train[gene_index][0]
                    mean_of_list = calculate.mean(genes_expression_no_relapse_train[gene_index][1])
                    sd_of_list = calculate.sd(genes_expression_no_relapse_train[gene_index][1])
                    result.append(gene_name)
                    result.append(mean_of_list)
                    result.append(sd_of_list)
                    list_mean_sd_gene_expression_by_probe_id_no_relapse.append(result)
                

                print(" Process : Calculate Log-likelihood ration of each gene ...")
                # calculate gene lambda
                # find mean and sd of gene expression in each class
                # default row_to_read_file_to_get_lambda = 22283
                row_to_read_file_to_get_lambda = 100

                col_to_read_file_to_get_lambda_relapse = ["ID_REF"]
                col_to_read_file_to_get_lambda_relapse.extend(list_sample_relapse)

                col_to_read_file_to_get_lambda_no_relapse = ["ID_REF"]
                col_to_read_file_to_get_lambda_no_relapse.extend(list_sample_no_relapse)

                file_to_get_lambda_relapse = pd.read_csv("GSE2034-22071 (edited).csv", usecols = col_to_read_file_to_get_lambda_relapse, nrows = row_to_read_file_to_get_lambda)
                file_to_get_lambda_no_relapse = pd.read_csv("GSE2034-22071 (edited).csv", usecols = col_to_read_file_to_get_lambda_no_relapse, nrows = row_to_read_file_to_get_lambda)
    
                # dictionary contains genes idintified by its probe id which contain all gene expression from samples
                # {1: [gene_probe_id, [exp1,exp2, ...]]}

                # for class "relapse"
                genes_expression_relapse = {}
                last_index_to_read_file_to_get_lambda_relapse = len(col_to_read_file_to_get_lambda_relapse)
                for line_index in range(0, row_to_read_file_to_get_lambda):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 
                    for element in file_to_get_lambda_relapse.iloc[line_index, 1:last_index_to_read_file_to_get_lambda_relapse]:
                        list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_get_lambda_relapse.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression_relapse[line_index] = gene_expression_by_probe_id

                # for class "non-relapse"
                genes_expression_no_relapse = {}
                last_index_to_read_file_to_get_lambda_no_relapse = len(col_to_read_file_to_get_lambda_no_relapse)
                for line_index in range(0, row_to_read_file_to_get_lambda):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_get_lambda_no_relapse.iloc[line_index, 1:last_index_to_read_file_to_get_lambda_no_relapse]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_get_lambda_no_relapse.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression_no_relapse[line_index] = gene_expression_by_probe_id
                
                # for class "relapse"
                genes_lambda_relapse = {}
                for gene_index in range(0, len(genes_expression_relapse)):
                    result = []
                    gene_name = genes_expression_relapse[gene_index][0]

                    gene_mean_relapse = list_mean_sd_gene_expression_by_probe_id_relapse[gene_index][1]
                    gene_sd_relapse = list_mean_sd_gene_expression_by_probe_id_relapse[gene_index][2]

                    gene_mean_no_relapse = list_mean_sd_gene_expression_by_probe_id_no_relapse[gene_index][1]
                    gene_sd_no_relapse = list_mean_sd_gene_expression_by_probe_id_no_relapse[gene_index][2]

                    list_gene_lambda = []
                    for gene_expression_index in range(0, len(genes_expression_relapse[gene_index][1])):
                        gene_expression = genes_expression_relapse[gene_index][1][gene_expression_index]

                        gene_pdf_under_relapse = norm.pdf(gene_expression, gene_mean_relapse, gene_sd_relapse)
                        gene_pdf_under_no_relapse = norm.pdf(gene_expression, gene_mean_no_relapse, gene_sd_no_relapse)

                        log_gene_pdf_under_relapse = math.log10(gene_pdf_under_relapse)
                        log_gene_pdf_under_no_relapse = math.log10(gene_pdf_under_no_relapse)

                        lambda_value = log_gene_pdf_under_relapse - log_gene_pdf_under_no_relapse

                        list_gene_lambda.append(lambda_value)
                    
                    result.append(gene_name)
                    result.append(list_gene_lambda)

                    genes_lambda_relapse[gene_index] = result
                
                # for class "non-relapse"
                genes_lambda_no_relapse = {}
                for gene_index in range(0, len(genes_expression_no_relapse)):
                    result = []
                    gene_name = genes_expression_no_relapse[gene_index][0]

                    gene_mean_relapse = list_mean_sd_gene_expression_by_probe_id_relapse[gene_index][1]
                    gene_sd_relapse = list_mean_sd_gene_expression_by_probe_id_relapse[gene_index][2]

                    gene_mean_no_relapse = list_mean_sd_gene_expression_by_probe_id_no_relapse[gene_index][1]
                    gene_sd_no_relapse = list_mean_sd_gene_expression_by_probe_id_no_relapse[gene_index][2]

                    list_gene_lambda = []
                    for gene_expression_index in range(0, len(genes_expression_no_relapse[gene_index][1])):
                        gene_expression = genes_expression_no_relapse[gene_index][1][gene_expression_index]

                        gene_pdf_under_relapse = norm.pdf(gene_expression, gene_mean_relapse, gene_sd_relapse)
                        gene_pdf_under_no_relapse = norm.pdf(gene_expression, gene_mean_no_relapse, gene_sd_no_relapse)

                        log_gene_pdf_under_relapse = math.log10(gene_pdf_under_relapse)
                        log_gene_pdf_under_no_relapse = math.log10(gene_pdf_under_no_relapse)

                        lambda_value = log_gene_pdf_under_relapse - log_gene_pdf_under_no_relapse

                        list_gene_lambda.append(lambda_value)
                    
                    result.append(gene_name)
                    result.append(list_gene_lambda)

                    genes_lambda_no_relapse[gene_index] = result
                
                # normalize lambda value
                list_gene_lambda_mean_sd = []
                for gene_index in range(0, len(genes_lambda_relapse)):
                    gene_mean_sd = []
                    gene_probe_id = genes_lambda_relapse[gene_index][0]

                    all_lambda_this_gene = []
                    all_lambda_this_gene.extend(genes_lambda_relapse[gene_index][1])
                    all_lambda_this_gene.extend(genes_lambda_no_relapse[gene_index][1])
                    
                    mean_lambda_this_gene = calculate.mean(all_lambda_this_gene)
                    sd_lambda_this_gene  = calculate.sd(all_lambda_this_gene)

                    gene_mean_sd.append(gene_probe_id)
                    gene_mean_sd.append(mean_lambda_this_gene)
                    gene_mean_sd.append(sd_lambda_this_gene)

                    list_gene_lambda_mean_sd.append(gene_mean_sd)

                genes_lambda_relapse_normalized = {}
                for gene_index in range(0, len(genes_lambda_relapse)):
                    gene_with_lambda = []
                    gene_probe_id = genes_lambda_relapse[gene_index][0]

                    list_gene_lambda = []
                    for lambda_index in range(0, len(genes_lambda_relapse[gene_index][1])):
                        lambda_value = genes_lambda_relapse[gene_index][1][lambda_index]

                        mean_this_gene = list_gene_lambda_mean_sd[gene_index][1]
                        sd_this_gene = list_gene_lambda_mean_sd[gene_index][2]

                        lambda_value_normalized = ((lambda_value - mean_this_gene) / sd_this_gene)

                        list_gene_lambda.append(lambda_value_normalized)
                    
                    gene_with_lambda.append(gene_probe_id)
                    gene_with_lambda.append(list_gene_lambda)

                    genes_lambda_relapse_normalized[gene_index] = gene_with_lambda
                
                genes_lambda_no_relapse_normalized = {}
                for gene_index in range(0, len(genes_lambda_no_relapse)):
                    gene_with_lambda = []
                    gene_probe_id = genes_lambda_no_relapse[gene_index][0]

                    list_gene_lambda = []
                    for lambda_index in range(0, len(genes_lambda_no_relapse[gene_index][1])):
                        lambda_value = genes_lambda_no_relapse[gene_index][1][lambda_index]

                        mean_this_gene = list_gene_lambda_mean_sd[gene_index][1]
                        sd_this_gene = list_gene_lambda_mean_sd[gene_index][2]

                        lambda_value_normalized = ((lambda_value - mean_this_gene) / sd_this_gene)

                        list_gene_lambda.append(lambda_value_normalized)
                    
                    gene_with_lambda.append(gene_probe_id)
                    gene_with_lambda.append(list_gene_lambda)

                    genes_lambda_no_relapse_normalized[gene_index] = gene_with_lambda

                
                print("Process : Creating collection to collect samples and their genes' expression")

                # create dictionary used to collect pathways of each sample
                samples_relapse = {}
                samples_no_relapse = {}

                # get all pathways of all samples 
                # for class 'relapse'
                for sample_index in range(0, len(list_sample_relapse)):
                    print()
                    print("Creating pathways for sample " + str(sample_index + 1) + " relapse is in progress ...")
                    print(str(len(list_sample_relapse) - (sample_index + 1)) + " samples left")
                    print()

                    # create list of lambda of genes in this sample
                    list_gene_lambda_value = []
                    for gene_index in range(0, len(genes_lambda_relapse_normalized)):
                        gene_with_lambda = []
                        for lambda_index in range(0, len(genes_lambda_relapse_normalized[gene_index][1])):
                            if (lambda_index == sample_index):
                                gene_probe_id = genes_lambda_relapse_normalized[gene_index][0]
                                lambda_value = genes_lambda_relapse_normalized[gene_index][1][lambda_index]

                                gene_with_lambda.append(gene_probe_id)
                                gene_with_lambda.append(lambda_value)

                                list_gene_lambda_value.append(gene_with_lambda)

                    sample = []
                    sample_name = list_sample_relapse[sample_index]
                    pathways = calculate.getPathwayLLR(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                                list_gene_lambda_value)

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_relapse[sample_index] = sample
                
                # for class "non-relapse"
                for sample_index in range(0, len(list_sample_no_relapse)):
                    print()
                    print("Creating pathways for sample " + str(sample_index + 1) + " relapse is in progress ...")
                    print(str(len(list_sample_no_relapse) - (sample_index + 1)) + " samples left")
                    print()

                    # create list of lambda of genes in this sample
                    list_gene_lambda_value = []
                    for gene_index in range(0, len(genes_lambda_no_relapse_normalized)):
                        gene_with_lambda = []
                        for lambda_index in range(0, len(genes_lambda_no_relapse_normalized[gene_index][1])):
                            if (lambda_index == sample_index):
                                gene_probe_id = genes_lambda_no_relapse_normalized[gene_index][0]
                                lambda_value = genes_lambda_no_relapse_normalized[gene_index][1][lambda_index]

                                gene_with_lambda.append(gene_probe_id)
                                gene_with_lambda.append(lambda_value)

                                list_gene_lambda_value.append(gene_with_lambda)

                    sample = []
                    sample_name = list_sample_no_relapse[sample_index]
                    pathways = calculate.getPathwayLLR(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                                list_gene_lambda_value)
                    
                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_no_relapse[sample_index] = sample
                


                
                
                    



if __name__ == "__main__":
    main()