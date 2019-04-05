import pandas as pd
import random
import math
import calculate
import time
import add_ons
from scipy import stats
from copy import deepcopy
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

def main():
    # record start time
    start_time = time.time()

    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Method : LLR-Based Classification")
    print(" # Experiment : Within Dataset")
    print(" # This method requires 2 datasets.")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains mapping between gene probe IDs and samples")
    print(" #   [2] Number of rows of the file containing mapping between gene probe IDs and samples to be read")
    print(" #   [3] A file contains mapping between samples and their class")  
    print(" #   [4] A file contains mapping between gene probe IDs and gene entrez IDs")
    print(" #   [5] A file contains pathways and their member genes")
    print(" #   [6] number of rows of the file contaning pathways and their member genes")
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked to provide required information to conduct an experiment including")
    print(" #   [1] Number of epochs")
    print(" #   [2] Number of folds")
    print(" #   [3] Number of pathways to be used as feature set (%)")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()

    # prepare variables
    file_training_input_name = None

    row_to_read_file_input = None
    
    file_training_output_name = None

    rows_to_read_file_pathway = None
    file_ref_name = None
    file_pathway_name = None

    num_of_epochs = None
    num_of_folds = None

    file_name = None

    print(" # Enter required information about the first dataset ")
    print(" 1. Enter name of the file containing mapping between probes IDs and samples ")
    file_training_input_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_input = input(" Number of rows : ")
        if (row_to_read_file_input.isnumeric() == False):
            print(" WARNING : Number of rows must be numeric.")
        elif (int(row_to_read_file_input) < 1):
            print(" WARNING : Number of rows cannot be lower than 1.") 
        else:
            break
    row_to_read_file_input = int(row_to_read_file_input)
    print()

    print(" 3. Enter name of a file containing mapping between samples and their class")
    file_training_output_name = add_ons.checkFileValid()
    print()

    print(" # Enter required information about files related to pathway mapping")
    print(" 1. Enter name of a file containing mapping between gene probe IDs and gene entrez IDs")
    file_ref_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter name of a file containing pathways and their member genes")
    file_pathway_name = add_ons.checkFileValid()
    print()

    print(" 3. Enter number of rows of the file contaning pathways and their member genes")
    while True:
        rows_to_read_file_pathway = input(" Number of rows : ")
        if(rows_to_read_file_pathway.isnumeric() == False):
            print(" WARNING : Number of rows must be numeric.")
        elif (int(rows_to_read_file_pathway) < 1):
            print(" WARNING : Number of rows cannot be lower than 1.") 
        else:
            break
    rows_to_read_file_pathway = int(rows_to_read_file_pathway)
    print()

    # prepare data
    # default row_to_read for "GSE2034-22071 (edited).csv" = 22283
    file_training_input = pd.read_csv(file_training_input_name, nrows = row_to_read_file_input)
    file_training_output= pd.read_csv(file_training_output_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway for "c2.cp.v6.2.entrez.gmt.csv" = 1329
    file_to_convert_name = file_training_input_name
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

    print(" # Enter required information to conduct an experiment")
    print(" 1. Enter number of epochs ")
    while True:
        num_of_epochs = input(" Epochs : ")

        if (num_of_epochs.isnumeric() == False):
            print(" WARNING : Number of epochs must be numeric.")
        elif (int(num_of_epochs) <= 0):
            print(" WARINING : Number of epochs must be greater than 0.")
        else:
            break
    num_of_epochs = int(num_of_epochs)
    print()

    print(" 2. Enter number of folds ")
    while True:
        num_of_folds = input(" Number of folds: ")
        if (num_of_folds.isnumeric() == False):
            print("WARNING : Number of folds must be numeric")
        
        # these conditions are not available in mock-up
        elif(int(num_of_folds) > len(list_sample_relapse)):
            print("WARNING : Number of folds exceeds the size of samples in class relapse.")
        elif(int(num_of_folds) > len(list_sample_no_relapse)):
            print("WARNING : Number of folds exceeds the size of samples in class non-relapse.")

        elif(int(num_of_folds) <= 1):
            print("WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)    
    print()

    file_name = input(" # Enter name of an output file : ")

    # prepare text file for results to be written in
    result_file = open("./result/" +str(file_name) + ".txt", "w+")

    # record dataset
    result_file.write("Dataset : " + str(file_training_input_name) + "\n")
    result_file.write("Pathway reference : " +str(file_pathway_name) + "\n")
    result_file.write("Number of epochs : " + str(num_of_epochs) + "\n")
    result_file.write("Number of folds : " + str(num_of_folds) + "\n")
    result_file.write("\n")

    # create list of all samples 
    list_all_samples = []
    list_all_samples.extend(list_sample_relapse)
    list_all_samples.extend(list_sample_no_relapse)
    
    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    for epoch_count in range(0, num_of_epochs):
        start_epoch_time = time.time()

        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################\n")
        
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
        print(" Number of chunks in class relapse = " + str(len(chunk_list_relapse)))

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse, chunk_no_relapse_size))
        print(" Number of chunks in class non-relapse = " + str(len(chunk_list_no_relapse)))

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # list to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        # list to collect maximun AUC in each fold
        list_max_auc = []

        print()
        print(" # Process : Cross-validation")
        print()
        # do only if number of chunks of both datasets are equal
        if (check_valid == True):
            for chunk_test_index in range(0, num_of_chunks):

                start_fold_time = time.time()
                result_file.write("\n#### Fold " + str(chunk_test_index + 1) + " ####\n")

                # separating data into testing and training dataset
                # get testing set in this fold
                chunk_test_relapse = chunk_list_relapse[chunk_test_index]
                chunk_test_no_relapse = chunk_list_no_relapse[chunk_test_index]

                print("\n------------------------------------------ K : " + str(chunk_test_index + 1) + " --------------------------------")
  
                # get training set in this fold
                chunk_train_relapse = []
                for chunk_train_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_relapse[chunk_train_relapse_index] is not chunk_test_relapse):
                        chunk_train_relapse.append(chunk_list_relapse[chunk_train_relapse_index])
                # print("chunk train relapse size = " + str(len(chunk_train_relapse)))

                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse[chunk_train_no_relapse_index])
                # print("chunk train no relapse size = " + str(len(chunk_train_no_relapse)))

                 # merge training data of each class
                list_train_relapse = []
                for i in range(0, len(chunk_train_relapse)):
                    list_train_relapse.extend(chunk_train_relapse[i])

                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])

                # get sample name and add to a list to be used as column index
                list_train_relapse_name = []
                list_train_no_relapse_name = []
                for i in range(0, len(list_train_relapse)):
                    list_train_relapse_name.append(list_sample_relapse[list_train_relapse[i]])
                print(" Samples in class relapse used as training set : " + str(list_train_relapse_name) + "\n")
                
                for i in range(0, len(list_train_no_relapse)):
                    list_train_no_relapse_name.append(list_sample_no_relapse[list_train_no_relapse[i]])
                print(" Samples in class non-relapsed used as training set : " + str(list_train_no_relapse_name) + "\n")

                # find mean and sd of gene expression of each gene in each class
                # prepare file to get gene expression

                row_to_read_file_to_cal_mean_sd = row_to_read_file_input

                col_to_read_file_to_cal_mean_sd_relapse = ["ID_REF"]
                col_to_read_file_to_cal_mean_sd_relapse.extend(list_train_relapse_name)

                col_to_read_file_to_cal_mean_sd_no_relapse = ["ID_REF"]
                col_to_read_file_to_cal_mean_sd_no_relapse.extend(list_train_no_relapse_name)

                file_to_cal_mean_sd_relapse = pd.read_csv(file_training_input_name, usecols = col_to_read_file_to_cal_mean_sd_relapse, nrows = row_to_read_file_to_cal_mean_sd)
                file_to_cal_mean_sd_no_relapse = pd.read_csv(file_training_input_name, usecols = col_to_read_file_to_cal_mean_sd_no_relapse, nrows = row_to_read_file_to_cal_mean_sd)
                
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
                

                print(" # Process : Calculate Log-likelihood ration of each gene")
                # calculate gene lambda
                # find mean and sd of gene expression in each class
                row_to_read_file_to_get_lambda = row_to_read_file_input

                col_to_read_file_to_get_lambda_relapse = ["ID_REF"]
                col_to_read_file_to_get_lambda_relapse.extend(list_sample_relapse)

                col_to_read_file_to_get_lambda_no_relapse = ["ID_REF"]
                col_to_read_file_to_get_lambda_no_relapse.extend(list_sample_no_relapse)

                file_to_get_lambda_relapse = pd.read_csv(file_training_input_name, usecols = col_to_read_file_to_get_lambda_relapse, nrows = row_to_read_file_to_get_lambda)
                file_to_get_lambda_no_relapse = pd.read_csv(file_training_input_name, usecols = col_to_read_file_to_get_lambda_no_relapse, nrows = row_to_read_file_to_get_lambda)
    
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

                        log_gene_pdf_under_relapse = norm.logpdf(gene_expression, gene_mean_relapse, gene_sd_relapse)
                        log_gene_pdf_under_no_relapse = norm.logpdf(gene_expression, gene_mean_no_relapse, gene_sd_no_relapse)            

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

                        log_gene_pdf_under_relapse = norm.logpdf(gene_expression, gene_mean_relapse, gene_sd_relapse)
                        log_gene_pdf_under_no_relapse = norm.logpdf(gene_expression, gene_mean_no_relapse, gene_sd_no_relapse)

                        lambda_value = log_gene_pdf_under_relapse - log_gene_pdf_under_no_relapse

                        list_gene_lambda.append(lambda_value)
                    
                    result.append(gene_name)
                    result.append(list_gene_lambda)

                    genes_lambda_no_relapse[gene_index] = result
                
                print(" # Process : Normalize log-likelihood ratio of gene in samples")
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

                
                print(" # Process : Creating collection to collect samples and their genes' expression")

                # create dictionary used to collect pathways of each sample
                samples_relapse = {}
                samples_no_relapse = {}

                # get all pathways of all samples 
                # for class 'relapse'
                for sample_index in range(0, len(list_sample_relapse)):
                    print()
                    print(" Creating pathways for sample " + str(sample_index + 1) + " relapse is in progress ...")
                    print(" " + str(len(list_sample_relapse) - (sample_index + 1)) + " samples left")
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
                    print(" Creating pathways for sample " + str(sample_index + 1) + " non-relapse is in progress ...")
                    print(" " + str(len(list_sample_no_relapse) - (sample_index + 1)) + " samples left")
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

                # splitting lists to use them as a marker evaluation set and feature selection set
                # given that we use 3-fold cross validation in this level
                print("\n#### given that we use 3 parts in this level for marker evaluation and feature selection set ####\n")
                second_num_of_fold = 3
                second_chunk_relapse_size = math.ceil(len(list_train_relapse) / second_num_of_fold)
                second_chunk_no_relapse_size = math.ceil(len(list_train_no_relapse) / second_num_of_fold)

                second_chunk_list_relapse = list(calculate.chunks(list_train_relapse, second_chunk_relapse_size))
                print(" Number of chunks in class relapse : " + str(len(second_chunk_list_relapse)) + "\n")
                second_chunk_list_no_relapse = list(calculate.chunks(list_train_no_relapse, second_chunk_no_relapse_size))
                print(" Number of chunks in class non-relapse : " + str(len(second_chunk_list_no_relapse)) + "\n")

                second_check_valid, second_num_of_chunks = calculate.checkEqualListSize(second_chunk_list_relapse, second_chunk_list_no_relapse)

                # variable to collect data from sfs
                feature_set_name = None
                auc_score_feature_selection = None

                print(" # Process : Feature Selection\n")
                if (second_check_valid is True):
                    # random a chunk of data to be use as a feature selection set
                    feature_selection_index = random.randint(0, second_num_of_chunks - 1)

                    # get feature selection set from each class
                    feature_selection_relapse =  second_chunk_list_relapse[feature_selection_index]
                    feature_selection_no_relapse = second_chunk_list_no_relapse[feature_selection_index]

                    # separate marker evaluation dataset from feature selection dataset
                    marker_evaluation_relapse = []
                    for marker_evaluation_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_relapse[marker_evaluation_index] is not feature_selection_relapse):
                            marker_evaluation_relapse.append(second_chunk_list_relapse[marker_evaluation_index])

                    marker_evaluation_no_relapse = []
                    for marker_evaluation_index in range(0, second_num_of_chunks):
                        if (second_chunk_list_no_relapse[marker_evaluation_index] is not feature_selection_no_relapse):
                            marker_evaluation_no_relapse.append(second_chunk_list_no_relapse[marker_evaluation_index])
                
                    # merge all samples in the same class
                    list_sample_relapse_marker_evaluation = []
                    for i in range(0, len(marker_evaluation_relapse)):
                        list_sample_relapse_marker_evaluation.extend(marker_evaluation_relapse[i])

                    list_sample_no_relapse_marker_evaluation = []
                    for i in range(0, len(marker_evaluation_no_relapse)):
                        list_sample_no_relapse_marker_evaluation.extend(marker_evaluation_no_relapse[i])

                    # create collection of samples used in feature selection
                    samples_relapse_marker_evaluation = {}
                    samples_no_relapse_marker_evaluation= {}

                    for sample_index in range(0, len(list_sample_relapse_marker_evaluation)):
                        index_samples_relapse = list_sample_relapse_marker_evaluation[sample_index]
                        samples_relapse_marker_evaluation[sample_index] = samples_relapse[index_samples_relapse]

                    for sample_index in range(0, len(list_sample_no_relapse_marker_evaluation)):
                        index_samples_no_relapse = list_sample_no_relapse_marker_evaluation[sample_index]
                        samples_no_relapse_marker_evaluation[sample_index] = samples_no_relapse[index_samples_no_relapse]

                    # create list contain pathway activity of each samples to be used in sfs
                    samples_relapse_marker_evaluation_pathway_activity = {}
                    samples_no_relapse_marker_evaluation_pathway_activity = {}

                    # for class 'relapse'
                    for sample_index in range(0, len(samples_relapse_marker_evaluation)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_relapse_marker_evaluation[sample_index][1])):
                            pathway = []
                            sum_gene_lambda = 0
                            for gene_index in range(0, len(samples_relapse_marker_evaluation[sample_index][1][pathway_index])):
                                gene_lambda_value =  samples_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
                                sum_gene_lambda += gene_lambda_value
                            pathway_name = samples_relapse_marker_evaluation[sample_index][1][pathway_index][0]

                            pathway.append(pathway_name)
                            pathway.append(sum_gene_lambda)

                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_relapse_marker_evaluation[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)

                        samples_relapse_marker_evaluation_pathway_activity[sample_index] = list_sample_with_pathway_activity

                    # for class 'non-relapse'
                    for sample_index in range(0, len(samples_no_relapse_marker_evaluation)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_no_relapse_marker_evaluation[sample_index][1])):
                            pathway = []
                            sum_gene_lambda = 0
                            for gene_index in range(0, len(samples_no_relapse_marker_evaluation[sample_index][1][pathway_index])):
                                gene_lambda_value =  samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
                                sum_gene_lambda += gene_lambda_value
                            pathway_name = samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][0]

                            pathway.append(pathway_name)
                            pathway.append(sum_gene_lambda)

                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_no_relapse_marker_evaluation[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)

                        samples_no_relapse_marker_evaluation_pathway_activity[sample_index] = list_sample_with_pathway_activity
                    
                    # create samples for feature selection
                    samples_relapse_feature_selection = {}
                    samples_no_relapse_feature_selection = {}

                    # for class 'relapse'
                    for sample_index in range(0, len(feature_selection_relapse)):
                        index_samples_relapse = feature_selection_relapse[sample_index]
                        samples_relapse_feature_selection[sample_index] = samples_relapse[index_samples_relapse]

                    # for class 'non-relapse'
                    for sample_index in range(0, len(feature_selection_no_relapse)):
                        index_samples_no_relapse = feature_selection_no_relapse[sample_index]
                        samples_no_relapse_feature_selection[sample_index] = samples_no_relapse[index_samples_no_relapse]
                    
                    # calculate pathway activity 
                    samples_relapse_feature_selection_pathway_activity = {}
                    samples_no_relapse_feature_selection_pathway_activity = {}

                    # for class 'relapse'
                    for sample_index in range(0, len(samples_relapse_feature_selection)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_relapse_feature_selection[sample_index][1])):
                            pathway = []
                            sum_gene_lambda = 0
                            for gene_index in range(0, len(samples_relapse_feature_selection[sample_index][1][pathway_index])):
                                gene_lambda_value =  samples_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                                sum_gene_lambda += gene_lambda_value
                            pathway_name = samples_relapse_feature_selection[sample_index][1][pathway_index][0]

                            pathway.append(pathway_name)
                            pathway.append(sum_gene_lambda)

                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_relapse_feature_selection[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)

                        samples_relapse_feature_selection_pathway_activity[sample_index] = list_sample_with_pathway_activity
                    
                    # for class 'non-relapse'
                    for sample_index in range(0, len(samples_no_relapse_feature_selection)):
                        list_sample_with_pathway_activity = []
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_no_relapse_feature_selection[sample_index][1])):
                            pathway = []
                            sum_gene_lambda = 0
                            for gene_index in range(0, len(samples_no_relapse_feature_selection[sample_index][1][pathway_index])):
                                gene_lambda_value =  samples_no_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                                sum_gene_lambda += gene_lambda_value
                            pathway_name = samples_no_relapse_feature_selection[sample_index][1][pathway_index][0]

                            pathway.append(pathway_name)
                            pathway.append(sum_gene_lambda)

                            list_pathway_activity.append(pathway)
                        
                        sample_name = samples_no_relapse_feature_selection[sample_index][0]

                        list_sample_with_pathway_activity.append(sample_name)
                        list_sample_with_pathway_activity.append(list_pathway_activity)

                        samples_no_relapse_feature_selection_pathway_activity[sample_index] = list_sample_with_pathway_activity
                    
                    # create list contains pathway activity of each sample preparing for calculating p-value
                    # for class 'relapse'
                    list_relapse_pathway_activity_for_pvalue = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway = []
                        for sample_index in range(0, len(samples_relapse_marker_evaluation_pathway_activity)):
                            pathway_activity = samples_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][1]
                            pathway.append(pathway_activity)
                        list_relapse_pathway_activity_for_pvalue.append(pathway)

                    # for class 'non-relapse'
                    list_no_relapse_pathway_activity_for_pvalue = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway = []
                        for sample_index in range(0, len(samples_no_relapse_marker_evaluation_pathway_activity)):
                            pathway_activity = samples_no_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][1]
                            pathway.append(pathway_activity)
                        list_no_relapse_pathway_activity_for_pvalue.append(pathway)
                    
                    # calculate p-value
                    list_pvalue_pathway_activity = []
                    for pathway_index in range(0, rows_to_read_file_pathway):
                        pathway_pvalue = []
                        pathway_name = list_pathway_name[pathway_index][1]

                        pvalue = stats.ttest_ind(list_relapse_pathway_activity_for_pvalue[pathway_index], list_no_relapse_pathway_activity_for_pvalue[pathway_index], equal_var = False)[1]

                        pathway_pvalue.append(pathway_name)
                        pathway_pvalue.append(pvalue)
                        list_pvalue_pathway_activity.append(pathway_pvalue)
                    
                    # sort pathway using p-value in ascending order
                    list_pvalue_pathway_activity.sort(key = lambda x : x[1], reverse = False)

                    # reorder pathway in each sample
                    samples_relapse_marker_evaluation_pathway_activity_sorted = {}
                    samples_no_relapse_marker_evaluation_pathway_activity_sorted ={}
                    samples_relapse_feature_selection_pathway_activity_sorted = {}
                    samples_no_relapse_feature_selection_pathway_activity_sorted = {}

                    # for marker evaluation class "relapse"
                    for sample_index in range(0, len(samples_relapse_marker_evaluation_pathway_activity)):
                        sample_name = samples_relapse_marker_evaluation_pathway_activity[sample_index][0]
                        list_relapse_train_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_relapse_marker_evaluation_pathway_activity[sample_index][1])):
                                pathway_name = samples_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][1]
                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_relapse_train_pathway_activity_sorted.append(sample_name)
                        list_relapse_train_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_relapse_marker_evaluation_pathway_activity_sorted[sample_index] = list_relapse_train_pathway_activity_sorted
                    
                    # for marker evaluation class "non-relapse"
                    for sample_index in range(0, len(samples_no_relapse_marker_evaluation_pathway_activity)):
                        sample_name = samples_no_relapse_marker_evaluation_pathway_activity[sample_index][0]
                        list_no_relapse_train_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_no_relapse_marker_evaluation_pathway_activity[sample_index][1])):
                                pathway_name = samples_no_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_no_relapse_marker_evaluation_pathway_activity[sample_index][1][pathway_index][1]
                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_no_relapse_train_pathway_activity_sorted.append(sample_name)
                        list_no_relapse_train_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_no_relapse_marker_evaluation_pathway_activity_sorted[sample_index] = list_no_relapse_train_pathway_activity_sorted

                    # for feature selection class "relapse"
                    for sample_index in range(0, len(samples_relapse_feature_selection_pathway_activity)):
                        sample_name = samples_relapse_feature_selection_pathway_activity[sample_index][0]
                        list_relapse_validation_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_relapse_feature_selection_pathway_activity[sample_index][1])):
                                pathway_name = samples_relapse_feature_selection_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_relapse_feature_selection_pathway_activity[sample_index][1][pathway_index][1]                              
                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_relapse_validation_pathway_activity_sorted.append(sample_name)
                        list_relapse_validation_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_relapse_feature_selection_pathway_activity_sorted[sample_index] = list_relapse_validation_pathway_activity_sorted

                    # for feature selection class "non-relapse"
                    for sample_index in range(0, len(samples_no_relapse_feature_selection_pathway_activity)):
                        sample_name = samples_no_relapse_feature_selection_pathway_activity[sample_index][0]
                        list_no_relapse_validation_pathway_activity_sorted = []
                        list_pathway_activity_sorted = []
                        for pvalue_index in range(0, len(list_pvalue_pathway_activity)):
                            pathway = []
                            for pathway_index in range(0, len(samples_no_relapse_feature_selection_pathway_activity[sample_index][1])):
                                pathway_name = samples_no_relapse_feature_selection_pathway_activity[sample_index][1][pathway_index][0]
                                pathway_activity = samples_no_relapse_feature_selection_pathway_activity[sample_index][1][pathway_index][1]                              
                                if (pathway_name == list_pvalue_pathway_activity[pvalue_index][0]):
                                    pathway.append(pathway_name)
                                    pathway.append(pathway_activity)
                                    list_pathway_activity_sorted.append(pathway)
                        list_no_relapse_validation_pathway_activity_sorted.append(sample_name)
                        list_no_relapse_validation_pathway_activity_sorted.append(list_pathway_activity_sorted)
                        samples_no_relapse_feature_selection_pathway_activity_sorted[sample_index] = list_no_relapse_validation_pathway_activity_sorted

                    # get sample name of data in feature selection sample for feature selection
                    # for class 'relapse'
                    feature_selection_relapse_name = []
                    for index in range(0, len(feature_selection_relapse)):
                        index_samples_relapse = feature_selection_relapse[index]
                        feature_selection_relapse_name.append(samples_relapse[index_samples_relapse][0])
                    
                    # for class 'non-relapse'
                    feature_selection_no_relapse_name = []
                    for index in range(0, len(feature_selection_no_relapse)):
                        index_samples_no_relapse = feature_selection_no_relapse[index]
                        feature_selection_no_relapse_name.append(samples_no_relapse[index_samples_no_relapse][0])
                    
                    # merge testing data to be used in lda for feature selection 
                    feature_selection_all_name = []
                    feature_selection_all_name.extend(feature_selection_relapse_name)
                    feature_selection_all_name.extend(feature_selection_no_relapse_name)

                    # create list of desired output
                    file_desired_outputs_feature_selection = file_training_output.loc[file_training_output['GEO asscession number'].isin(feature_selection_all_name)]
                    file_desired_outputs_feature_selection['sample_id'] = file_desired_outputs_feature_selection['GEO asscession number'].apply(lambda name: feature_selection_all_name.index(name)) 
                    file_desired_outputs_feature_selection = file_desired_outputs_feature_selection.sort_values(by = ['sample_id'])
                    file_desired_outputs_feature_selection.drop(columns = 'sample_id', inplace = True)

                    list_desired_outputs_feature_selection = []
                    for element in file_desired_outputs_feature_selection.loc[:, 'relapse (1=True)']:
                        list_desired_outputs_feature_selection.append(element)

                    # create list of pathway name in the same order as in each sample
                    list_pathway_name_feature_selection = []
                    for pathway_index in range(0, len(samples_relapse_feature_selection_pathway_activity_sorted[0][1])):
                        pathway_name = samples_relapse_feature_selection_pathway_activity_sorted[0][1][pathway_index][0]
                        list_pathway_name_feature_selection.append(pathway_name)

                    # create list contain pathway activity of each samples to be used in sfs

                    # for marker evaluation class "relapse"
                    list_sample_relapse_pathway_expression_marker_evaluation = []
                    for sample_index in range(0, len(samples_relapse_marker_evaluation_pathway_activity_sorted)):
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1])):
                            pathway = []

                            pathway_name = samples_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1][pathway_index][0]
                            pathway_activity = samples_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        list_sample_relapse_pathway_expression_marker_evaluation.append(list_pathway_activity)
                    
                    # for marker evaluation class "non-relapse"
                    list_sample_no_relapse_pathway_expression_marker_evaluation = []
                    for sample_index in range(0, len(samples_no_relapse_marker_evaluation_pathway_activity_sorted)):
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_no_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1])):
                            pathway = []

                            pathway_name = samples_no_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1][pathway_index][0]
                            pathway_activity = samples_no_relapse_marker_evaluation_pathway_activity_sorted[sample_index][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        list_sample_no_relapse_pathway_expression_marker_evaluation.append(list_pathway_activity)

                    # for feature selection class "relapse"
                    list_sample_relapse_pathway_expression_feature_selection = []
                    for sample_index in range(0, len(samples_relapse_feature_selection_pathway_activity_sorted)):
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_relapse_feature_selection_pathway_activity_sorted[sample_index][1])):
                            pathway = []

                            pathway_name = samples_relapse_feature_selection_pathway_activity_sorted[sample_index][1][pathway_index][0]
                            pathway_activity = samples_relapse_feature_selection_pathway_activity_sorted[sample_index][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        list_sample_relapse_pathway_expression_feature_selection.append(list_pathway_activity)

                    # for feature selection class "non-relapse"
                    list_sample_no_relapse_pathway_expression_feature_selection = []
                    for sample_index in range(0, len(samples_no_relapse_feature_selection_pathway_activity_sorted)):
                        list_pathway_activity = []
                        for pathway_index in range(0, len(samples_no_relapse_feature_selection_pathway_activity_sorted[sample_index][1])):
                            pathway = []

                            pathway_name = samples_no_relapse_feature_selection_pathway_activity_sorted[sample_index][1][pathway_index][0]
                            pathway_activity = samples_no_relapse_feature_selection_pathway_activity_sorted[sample_index][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)
                            list_pathway_activity.append(pathway)
                        list_sample_no_relapse_pathway_expression_feature_selection.append(list_pathway_activity)
                    
                    # merge testing set for feature selection together
                    list_sample_all_pathway_expression_feature_selection = []
                    list_sample_all_pathway_expression_feature_selection.extend(list_sample_relapse_pathway_expression_feature_selection)
                    list_sample_all_pathway_expression_feature_selection.extend(list_sample_no_relapse_pathway_expression_feature_selection)

                    print(" # Process : Sequential Forward Selection (SFS)")
                    # find feature set using sequential forward selection
                    feature_set_name, auc_score_feature_selection = calculate.sfsAdvance(list_pathway_name_feature_selection, list_desired_outputs_feature_selection, list_sample_relapse_pathway_expression_marker_evaluation, \
                            list_sample_no_relapse_pathway_expression_marker_evaluation, list_sample_all_pathway_expression_feature_selection)

                    list_max_auc.append(auc_score_feature_selection)

                    print(" Feature set: " + str(feature_set_name))
                    print(" AUC score from feature selection : " + str(auc_score_feature_selection))
                    print()
                    result_file.write("Feature set: " + str(feature_set_name))
                    result_file.write("\n")
                    result_file.write("AUC score from feature selection : " + str(auc_score_feature_selection))
                    result_file.write("\n")

                print(" # Process : Evaluation")
                # preparing data for evaluation and creating classifier
                # create classifier
                # for class 'relapse'
                samples_classifier_relapse = {}
                for sample_index in range(0,len(list_train_relapse)):
                    index_samples_relapse = list_train_relapse[sample_index]
                    samples_classifier_relapse[sample_index] = samples_relapse[index_samples_relapse]
                
                # for class 'non-relapse'
                samples_classifier_no_relapse = {}
                for sample_index in range(0,len(list_train_no_relapse)):
                    index_samples_no_relapse = list_train_no_relapse[sample_index]
                    samples_classifier_no_relapse[sample_index] = samples_no_relapse[index_samples_no_relapse]

                # calculate pathway activity of classifier
                # for class 'relapse'
                samples_classifier_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_classifier_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_classifier_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_lambda = 0
                        for gene_index in range(0, len(samples_classifier_relapse[sample_index][1][pathway_index])):
                            gene_lambda_value =  samples_classifier_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            sum_gene_lambda += gene_lambda_value
                        pathway_name = samples_classifier_relapse[sample_index][1][pathway_index][0]

                        pathway.append(pathway_name)
                        pathway.append(sum_gene_lambda)

                        list_pathway_activity.append(pathway)
                    
                    sample_name = samples_classifier_relapse[sample_index][0]

                    list_sample_with_pathway_activity.append(sample_name)
                    list_sample_with_pathway_activity.append(list_pathway_activity)

                    samples_classifier_relapse_pathway_activity[sample_index] = list_sample_with_pathway_activity
                
                # for class 'non-relapse'
                samples_classifier_no_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_classifier_no_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_classifier_no_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_lambda = 0
                        for gene_index in range(0, len(samples_classifier_no_relapse[sample_index][1][pathway_index])):
                            gene_lambda_value =  samples_classifier_no_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            sum_gene_lambda += gene_lambda_value
                        pathway_name = samples_classifier_no_relapse[sample_index][1][pathway_index][0]

                        pathway.append(pathway_name)
                        pathway.append(sum_gene_lambda)

                        list_pathway_activity.append(pathway)
                    
                    sample_name = samples_classifier_no_relapse[sample_index][0]

                    list_sample_with_pathway_activity.append(sample_name)
                    list_sample_with_pathway_activity.append(list_pathway_activity)

                    samples_classifier_no_relapse_pathway_activity[sample_index] = list_sample_with_pathway_activity
                
                # create testing data
                # for class 'relapse'
                samples_testing_relapse = {}
                for sample_index in range(0,len(chunk_test_relapse)):
                    index_samples_relapse = chunk_test_relapse[sample_index]
                    samples_testing_relapse[sample_index] = samples_relapse[index_samples_relapse]
                
                # for class 'non-relapse'
                samples_testing_no_relapse = {}
                for sample_index in range(0,len(chunk_test_no_relapse)):
                    index_samples_no_relapse = chunk_test_no_relapse[sample_index]
                    samples_testing_no_relapse[sample_index] = samples_no_relapse[index_samples_no_relapse]
                
                # calculate pathway activity
                # for class 'relapse'
                samples_testing_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_testing_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_testing_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_lambda = 0
                        for gene_index in range(0, len(samples_testing_relapse[sample_index][1][pathway_index])):
                            gene_lambda_value =  samples_testing_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            sum_gene_lambda += gene_lambda_value
                        pathway_name = samples_testing_relapse[sample_index][1][pathway_index][0]

                        pathway.append(pathway_name)
                        pathway.append(sum_gene_lambda)

                        list_pathway_activity.append(pathway)
                    
                    sample_name = samples_testing_relapse[sample_index][0]

                    list_sample_with_pathway_activity.append(sample_name)
                    list_sample_with_pathway_activity.append(list_pathway_activity)

                    samples_testing_relapse_pathway_activity[sample_index] = list_sample_with_pathway_activity
                
                # for class 'non-relapse'
                samples_testing_no_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_testing_no_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_testing_no_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_lambda = 0
                        for gene_index in range(0, len(samples_testing_no_relapse[sample_index][1][pathway_index])):
                            gene_lambda_value =  samples_testing_no_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            sum_gene_lambda += gene_lambda_value
                        pathway_name = samples_testing_no_relapse[sample_index][1][pathway_index][0]

                        pathway.append(pathway_name)
                        pathway.append(sum_gene_lambda)

                        list_pathway_activity.append(pathway)
                    
                    sample_name = samples_testing_no_relapse[sample_index][0]

                    list_sample_with_pathway_activity.append(sample_name)
                    list_sample_with_pathway_activity.append(list_pathway_activity)

                    samples_testing_no_relapse_pathway_activity[sample_index] = list_sample_with_pathway_activity
                
                # create list to be used in lda
                # for classifier class "relapse"
                list_classifier_relapse_pathway_expression = []
                for sample_index in range(0, len(samples_classifier_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for feature in feature_set_name:
                        for pathway_index in range(0, len(samples_classifier_relapse_pathway_activity[sample_index][1])):
                            pathway_name = samples_classifier_relapse_pathway_activity[sample_index][1][pathway_index][0]
                            pathway_activity = samples_classifier_relapse_pathway_activity[sample_index][1][pathway_index][1]

                            if (pathway_name == feature):
                                list_pathway_activity.append(pathway_activity)

                    list_classifier_relapse_pathway_expression.append(list_pathway_activity)  

                # for classifier class "non-relapse"
                list_classifier_no_relapse_pathway_expression = []
                for sample_index in range(0, len(samples_classifier_no_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for feature in feature_set_name:
                        for pathway_index in range(0, len(samples_classifier_no_relapse_pathway_activity[sample_index][1])):
                            pathway_name = samples_classifier_no_relapse_pathway_activity[sample_index][1][pathway_index][0]
                            pathway_activity = samples_classifier_no_relapse_pathway_activity[sample_index][1][pathway_index][1]

                            if (pathway_name ==  feature):
                                list_pathway_activity.append(pathway_activity)

                    list_classifier_no_relapse_pathway_expression.append(list_pathway_activity)

                # for testing data class "relapse"
                list_testing_relapse_pathway_expression = []
                for sample_index in range(0, len(samples_testing_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for feature in feature_set_name:
                        for pathway_index in range(0, len(samples_testing_relapse_pathway_activity[sample_index][1])):
                            pathway_name = samples_testing_relapse_pathway_activity[sample_index][1][pathway_index][0]
                            pathway_activity = samples_testing_relapse_pathway_activity[sample_index][1][pathway_index][1]

                            if (pathway_name == feature):
                                list_pathway_activity.append(pathway_activity)

                    list_testing_relapse_pathway_expression.append(list_pathway_activity)

                # for testing data class "non-relapse"
                list_testing_no_relapse_pathway_expression = []
                for sample_index in range(0, len(samples_testing_no_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for feature in feature_set_name:
                        for pathway_index in range(0, len(samples_testing_no_relapse_pathway_activity[sample_index][1])):
                            pathway_name = samples_testing_no_relapse_pathway_activity[sample_index][1][pathway_index][0]
                            pathway_activity = samples_testing_no_relapse_pathway_activity[sample_index][1][pathway_index][1]

                            if (pathway_name == feature):
                                list_pathway_activity.append(pathway_activity)

                    list_testing_no_relapse_pathway_expression.append(list_pathway_activity)

                # merge testing data of 2 class together
                list_testing_all_pathway_expression = []
                list_testing_all_pathway_expression.extend(list_testing_relapse_pathway_expression)
                list_testing_all_pathway_expression.extend(list_testing_no_relapse_pathway_expression)

                # create list to contain all sample name used in testing procedure
                list_chunk_test_relapse_name = []
                for index in range(0, len(chunk_test_relapse)):
                    index_samples_relapse = chunk_test_relapse[index]
                    list_chunk_test_relapse_name.append(samples_relapse[index_samples_relapse][0])

                list_chunk_test_no_relapse_name = []
                for index in range(0, len(chunk_test_no_relapse)):
                        index_samples_no_relapse = chunk_test_no_relapse[index]
                        list_chunk_test_no_relapse_name.append(samples_no_relapse[index_samples_no_relapse][0])

                list_samples_name_testing_all = []
                list_samples_name_testing_all.extend(list_chunk_test_relapse_name)
                list_samples_name_testing_all.extend(list_chunk_test_no_relapse_name)

                # create list of desired outputs
                file_desired_outputs = file_training_output.loc[file_training_output['GEO asscession number'].isin(list_samples_name_testing_all)]

                file_desired_outputs['pathway_id'] = file_desired_outputs['GEO asscession number'].apply(lambda name: list_samples_name_testing_all.index(name)) 
                file_desired_outputs = file_desired_outputs.sort_values(by = ['pathway_id'])
                file_desired_outputs.drop(columns = 'pathway_id', inplace = True)

                list_desired_outputs = []
                for element in file_desired_outputs.loc[:, 'relapse (1=True)']:
                    list_desired_outputs.append(element)
                
                # calculate lda 
                list_actual_outputs = calculate.lda(list_testing_all_pathway_expression, list_classifier_relapse_pathway_expression, list_classifier_no_relapse_pathway_expression)

                # calculate AUC score
                auc_score = roc_auc_score(list_desired_outputs, list_actual_outputs)
                list_auc_score.append(auc_score)

                print()
                print("#### Evaluation of " + str(chunk_test_index + 1) + " - fold ####")
                print(" Feature set : " + str(feature_set_name))
                print(" Number of features in feature set: " + str(len(feature_set_name)))
                print(" Number of actual outputs : " + str(len(list_actual_outputs)))
                print(" Actual outputs : ")
                print(list_actual_outputs)
                print()
                print(" Number of desired outputs : " + str(len(list_desired_outputs)))
                print(" Desired outputs : ")
                print(list_desired_outputs)
                print(" AUC score from feature selection : " + str(auc_score_feature_selection))
                print(" AUC score from testing : " + str(auc_score))

                 # track feature set which gives maximum auc score
                if (auc_score > auc_score_max):
                    list_feature_set_max_auc = deepcopy(feature_set_name)
                    auc_score_max = auc_score
                
                result_file.write("Feature set : " + str(feature_set_name) + "\n")
                result_file.write("Number of features in feature set : " + str(len(feature_set_name)) + "\n")
                result_file.write("Number of actual outputs : " + str(len(list_actual_outputs)) + "\n")
                result_file.write(str(list_actual_outputs) + "\n")
                result_file.write("\n")
                result_file.write("Number of desired outputs : " + str(len(list_desired_outputs)) + "\n")
                result_file.write(str(list_desired_outputs) + "\n")
                result_file.write("AUC score from feature selection : " + str(auc_score_feature_selection) + "\n")
                result_file.write("AUC score from testing : " + str(auc_score) + "\n")
                result_file.write("\n")
                

                end_fold_time = time.time()
                fold_elapse_time_second = end_fold_time - start_fold_time
                fold_elapse_time_minute = fold_elapse_time_second / 60
                fold_elapse_time_minute = round(fold_elapse_time_minute, 2)
                print(" Fold elapse time : " + str(fold_elapse_time_minute) + " minutes")
                result_file.write("Fold elapse time : " + str(fold_elapse_time_minute) + " minutes \n")
                result_file.write("\n")
        
        end_epoch_time = time.time()
        time_elapse_epoch_second = end_epoch_time - start_epoch_time
        time_elapse_epoch_minute = time_elapse_epoch_second / 60
        time_elapse_epoch_hour = time_elapse_epoch_minute / 60

        time_elapse_epoch_minute = round(time_elapse_epoch_minute, 2)
        time_elapse_epoch_hour = round(time_elapse_epoch_hour, 2)

        list_avg_auc_each_epoch.append(calculate.mean(list_auc_score))

        print()
        print("#### Summary ####")
        print(" Average AUC score : " + str(calculate.mean(list_auc_score)))
        print(" Maximum AUC score : " + str(auc_score_max))
        print(" Feature set which gives highest AUC score : ")
        print(list_feature_set_max_auc)
        print()
        print(" Time elapse : "  + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours) ")

        result_file.write("\n#### Summary ####\n")

        result_file.write("Maximum AUC ROC score of feature from feature selection in each fold : \n")
        result_file.write(str(list_max_auc))
        result_file.write("\n")

        result_file.write("Average AUC score : " + str(calculate.mean(list_auc_score)) + "\n")
        result_file.write("Maximum AUC score : " + str(auc_score_max) + "\n")
        result_file.write("Size of feature set which gives the highest AUC score from testing : " + str(len(list_feature_set_max_auc)))
        result_file.write("\n")
        result_file.write("Feature set which gives the highest AUC score from testing : " + "\n")
        result_file.write(str(list_feature_set_max_auc))
        result_file.write("\n")
        result_file.write("Time elapse : "  + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours) ")
        result_file.write("\n")
        result_file.write("\n")

    end_time = time.time()
    total_elapse_time_second = end_time - start_time

    total_elapse_time_minute = total_elapse_time_second / 60
    total_elapse_time_hour = total_elapse_time_minute / 60  
    
    total_elapse_time_minute = round(total_elapse_time_minute, 2)
    total_elapse_time_hour = round(total_elapse_time_hour, 2)

    # calculate mean over all epoch
    mean_over_all_epoch = calculate.mean(list_avg_auc_each_epoch)
    print(" Average AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_over_all_epoch))
    result_file.write("Average AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_over_all_epoch) + "\n")

    print(" Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("\n")


    result_file.close()
                
if __name__ == "__main__":
    main()