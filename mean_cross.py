import pandas as pd
import random
import math
import calculate
import time
import xlsxwriter
import add_ons
from scipy import stats
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main():
    # record start time
    start_time = time.time()

    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Method : Mean-Based Classification")
    print(" # Experiment : Cross-Dataset")
    print(" # This method requires 2 datasets.")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains mapping between gene probe IDs and samples of the first dataset")
    print(" #   [2] Number of rows of the file containing mapping between gene probe IDs and samples of the first dataset to be read")
    print(" #   [3] A file contains mapping between samples and their class of the first dataset")  
    print(" #   [4] A file contains mapping between gene probe IDs and samples of the second dataset")
    print(" #   [5] Number of rows of the file contains mapping between gene probe IDs and samples of the second dataset to be read")
    print(" #   [6] A file contains mapping between samples and their class of the second dataset")
    print(" #   [7] A file contains mapping between gene probe IDs and gene entrez IDs")
    print(" #   [8] A file contains pathways and their member genes")
    print(" #   [9] number of rows of the file contaning pathways and their member genes")
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
    file_gene_first_dataset_name = None
    file_gene_second_dataset_name = None

    row_to_read_file_gene_first_dataset = None
    row_to_read_file_gene_second_dataset = None
    
    file_output_first_dataset_name = None
    file_output_second_dataset_name = None

    rows_to_read_file_pathway = None
    file_ref_name = None
    file_pathway_name = None

    num_of_epochs = None
    num_of_folds = None
    num_of_pathways_percentage = None

    file_name = None

    print(" # Enter required information about the first dataset ")
    print(" 1. Enter name of the file containing mapping between probes IDs and samples of the first dataset ")
    file_gene_first_dataset_name = add_ons.getFile()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_gene_first_dataset = input(" Number of rows : ")
        if(row_to_read_file_gene_first_dataset.isnumeric() == False):
            print(" WARNING : Number of rows must be numeric.")
        elif (int(row_to_read_file_gene_first_dataset) < 1):
            print(" WARNING : Number of rows cannot be lower than 1.")
        else:
            break
    row_to_read_file_gene_first_dataset = int(row_to_read_file_gene_first_dataset)     
    print()

    print(" 3. Enter name of a file containing mapping between samples and their class of the first dataset")
    file_output_first_dataset_name = add_ons.getFile()
    print()

    print(" # Enter required information about the second dataset ")
    print(" 1. Enter name of the file containing mapping between probes IDs and samples of the second dataset ")
    file_gene_second_dataset_name = add_ons.getFile()
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

    print(" 3. Enter name of a file containing mapping between samples and their class of the second dataset")
    file_output_second_dataset_name = add_ons.getFile()
    print()

    print(" # Enter required information about files related to pathway mapping")
    print(" 1. Enter name of a file containing mapping between gene probe IDs and gene entrez IDs")
    file_ref_name = add_ons.getFile()
    print()

    print(" 2. Enter name of a file containing pathways and their member genes")
    file_pathway_name = add_ons.getFile()
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
    # for 1st dataset
    # default row_to_read_file_gene_first_dataset for "GSE2034-22071 (edited).csv"= 22283
    file_gene_first_dataset = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset)

    # file_output_first_dataset_name = "mapping_sample_to_class_gse2034.csv"
    file_output_first_dataset = pd.read_csv(file_output_first_dataset_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # for 2nd dataset
    # default row_to_read_file_second_dataset for "GSE3494_GPL96.csv" = 22283
    file_gene_second_dataset = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset) 

    # file_output_second_dataset_name = "mapping_sample_to_class_gse3494.csv"
    file_output_second_dataset = pd.read_csv(file_output_second_dataset_name, usecols = ["GEO asscession number", "relapse (1=True)"])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway for "c2.cp.v6.2.entrez.gmt.csv"= 1329
    file_to_convert_first_dataset_name = file_gene_first_dataset_name
    file_to_convert_second_dataset_name = file_gene_second_dataset_name
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read_file_gene_first_dataset):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_gene_first_dataset.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)

    # get list of pathway name
    list_pathway_name = []
    for i in range(0, rows_to_read_file_pathway):
        pathway_name = []
        pathway_name.append(i)
        pathway_name.append(file_pathway.loc[i, "PATHWAY_NAME"])
        list_pathway_name.append(pathway_name)

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

    print(" 3. Enter number of pathways to be used as feature set (%)")
    while True:
        num_of_pathways_percentage = input(" Number of pathways to be used as feature set (%) : ")
        if (num_of_pathways_percentage.isnumeric() == False):
            print(" WARNING : Input must be numeric")
        elif (int(num_of_pathways_percentage) <= 0):
            print(" WARNING : Percentage must greater than 0.")
        elif (int(num_of_pathways_percentage) > 100):
            print(" WARNING : Percentage must lower than or equal to 100.")
        else:
            break
    num_of_pathways_percentage = int(num_of_pathways_percentage)
    print()

    file_name = input(" # Enter name of an output file : ")

    # prepare text file for results to be written in
    result_file = open("./result/" +str(file_name) + ".txt", "w+")

    # record dataset
    result_file.write("The first dataset : " + str(file_gene_first_dataset_name) + "\n")
    result_file.write("The second dataset : " + str(file_gene_second_dataset_name) + "\n")
    result_file.write("Pathway reference : " + str(file_pathway_name) + "\n")
    result_file.write("Number of epochs : " + str(num_of_epochs) + "\n")
    result_file.write("Number of folds : " + str(num_of_folds) + "\n")
    result_file.write("Number of pathways to be used as feature set (%) : " + str(num_of_pathways_percentage) + "\n")
    result_file.write("\n")

    # calculate number of pathways to be used
    num_of_ranked_pathways = (rows_to_read_file_pathway * (num_of_pathways_percentage / 100))
    num_of_ranked_pathways = math.ceil(num_of_ranked_pathways)

    # create list of all samples in trainig data in this fold
    list_all_samples = []
    list_all_samples.extend(list_sample_relapse_first_dataset)
    list_all_samples.extend(list_sample_no_relapse_first_dataset)

    print(" # Process : Creating collection to collect samples and their genes' expression")
    # create dictionary used to collect pathways of each sample 
    # for class "relapse "
    samples_relapse_first_dataset = {}
    for element_index in range(0, len(list_sample_relapse_first_dataset)):
        print()
        print(" Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
        print(str(len(list_sample_relapse_first_dataset) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_relapse_first_dataset[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_first_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_relapse_first_dataset[element_index] = sample
    
    # for class "non-relapse"
    samples_no_relapse_first_dataset = {}
    for element_index in range(0, len(list_sample_no_relapse_first_dataset)):
        print()
        print(" Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
        print(str(len(list_sample_no_relapse_first_dataset) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_no_relapse_first_dataset[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_first_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_no_relapse_first_dataset[element_index] = sample
    
    print(" # Process : Creating collections of samples with their pathways' activity ...")
    # create collections of samples with their pathways
    # data will be collected in this format
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}} 

    # for class "relapse"
    samples_relapse_pathway_activity_first_dataset = {}
    for samples_index in range(0, len(samples_relapse_first_dataset)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_relapse_first_dataset[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_relapse_first_dataset[samples_index][1][pathway_index][1])):
                gene_expression = samples_relapse_first_dataset[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_relapse_first_dataset[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)
            list_pathway.append(pathway)
        
        sample_name = samples_relapse_first_dataset[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_relapse_pathway_activity_first_dataset[samples_index] = sample
    
    # for class "non-relapse"
    samples_no_relapse_pathway_activity_first_dataset = {}
    for samples_index in range(0, len(samples_no_relapse_first_dataset)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_no_relapse_first_dataset[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_no_relapse_first_dataset[samples_index][1][pathway_index][1])):
                gene_expression = samples_no_relapse_first_dataset[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_no_relapse_first_dataset[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)
            list_pathway.append(pathway)
        
        sample_name = samples_no_relapse_first_dataset[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_no_relapse_pathway_activity_first_dataset[samples_index] = sample

    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    # list to collect feature counter
    list_feature_counter = []

    for epoch_count in range(0, num_of_epochs):
        start_epoch_time = time.time()

        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################")

        # create list of indexes used to indicate the position in the list
        # for 1st dataset
        list_index_samples_relapse_first_dataset = []
        list_index_samples_no_relapse_first_dataset = []

        for index in range(0, len(list_sample_relapse_first_dataset)):
            list_index_samples_relapse_first_dataset.append(index)
        
        for index in range(0, len(list_sample_no_relapse_first_dataset)):
            list_index_samples_no_relapse_first_dataset.append(index)
        
        # for 2nd dataset
        list_index_samples_relapse_second_dataset = []
        list_index_samples_no_relapse_second_dataset = []

        for index in range(0, len(list_sample_relapse_second_dataset)):
            list_index_samples_relapse_second_dataset.append(index)
        
        for index in range(0, len(list_sample_no_relapse_second_dataset)):
            list_index_samples_no_relapse_second_dataset.append(index)
        
        # shuffle it to make it flexible for epoch changed
        # for 1st dataset
        random.shuffle(list_index_samples_relapse_first_dataset)
        random.shuffle(list_index_samples_no_relapse_first_dataset)

        # for 2nd dataset
        random.shuffle(list_index_samples_relapse_second_dataset)
        random.shuffle(list_index_samples_no_relapse_second_dataset)

        # conduct feature selection on 1st dataset
        # split data into 3 parts
        num_of_chunk_feature_selection = 3
        chunk_relapse_size = math.ceil(len(list_index_samples_relapse_first_dataset) / num_of_chunk_feature_selection)
        chunk_no_relapse_size = math.ceil(len(list_index_samples_no_relapse_first_dataset) / num_of_chunk_feature_selection)

        chunk_list_relapse = list(calculate.chunks(list_index_samples_relapse_first_dataset, chunk_relapse_size))
        print("number of chunks in class relapse : " + str(len(chunk_list_relapse)) + "\n")

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse_first_dataset, chunk_no_relapse_size))
        print("number of in class non-relapse : " + str(len(chunk_list_no_relapse)) + "\n")

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # list to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        # variable to collect data from sfs
        feature_set_name = None
        auc_score_feature_selection = None

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

            # for class "non-relapse"
            list_sample_no_relapse_marker_evaluation = []
            for i in range(0, len(marker_evaluation_no_relapse)):
                list_sample_no_relapse_marker_evaluation.extend(marker_evaluation_no_relapse[i])
            
            # get column name
            # for class "relapse"
            list_sample_relapse_marker_evaluation_name = []
            for i in range(0, len(list_sample_relapse_marker_evaluation)):
                sample_index_in_list = list_sample_relapse_marker_evaluation[i]
                sample_name = list_sample_relapse_first_dataset[sample_index_in_list]
                list_sample_relapse_marker_evaluation_name.append(sample_name)
            print(" Samples in class relapse used as marker evaluation set : " + str(list_sample_relapse_marker_evaluation_name) + "\n")
            
            # for class "non-relapse"
            list_sample_no_relapse_marker_evaluation_name = []
            for i in range(0, len(list_sample_no_relapse_marker_evaluation)):
                sample_index_in_list = list_sample_no_relapse_marker_evaluation[i]
                sample_name = list_sample_no_relapse_first_dataset[sample_index_in_list]
                list_sample_no_relapse_marker_evaluation_name.append(sample_name)
            print(" Samples in class non-relapse used as marker evaluation set : " + str(list_sample_no_relapse_marker_evaluation_name) + "\n")
            
            # prepare file used to calculate t-score of each gene
            file_to_cal_gene_tscore_name = file_gene_first_dataset_name
            row_to_read_file_cal_gene_tscore = row_to_read_file_gene_first_dataset

            col_to_read_file_cal_gene_tscore = ["ID_REF"]
            col_to_read_file_cal_gene_tscore.extend(list_sample_relapse_first_dataset)
            col_to_read_file_cal_gene_tscore.extend(list_sample_no_relapse_first_dataset)

            col_to_read_file_cal_gene_tscore_relapse = ["ID_REF"]
            col_to_read_file_cal_gene_tscore_relapse.extend(list_sample_relapse_marker_evaluation_name)

            col_to_read_file_cal_gene_tscore_no_relapse = ["ID_REF"]
            col_to_read_file_cal_gene_tscore_no_relapse.extend(list_sample_no_relapse_marker_evaluation_name)

            file_to_cal_gene_tscore = pd.read_csv(file_to_cal_gene_tscore_name, usecols = col_to_read_file_cal_gene_tscore, nrows = row_to_read_file_cal_gene_tscore)
            file_to_cal_gene_tscore_relapse = pd.read_csv(file_to_cal_gene_tscore_name, usecols = col_to_read_file_cal_gene_tscore_relapse, nrows = row_to_read_file_cal_gene_tscore)
            file_to_cal_gene_tscore_no_relapse = pd.read_csv(file_to_cal_gene_tscore_name, usecols = col_to_read_file_cal_gene_tscore_no_relapse, nrows = row_to_read_file_cal_gene_tscore)

            print()
            print(" # Process : Gathering gene expression from file ...")

            # get gene expression 
            # for class "relapse"
            genes_expression_relapse_first_dataset = {}
            last_index_to_read_file_to_cal_gene_tscore_relapse = len(col_to_read_file_cal_gene_tscore_relapse)
            for line_index in range(0, row_to_read_file_cal_gene_tscore):
                gene_expression_by_probe_id = []              
                list_gene_expression_same_probe_id = [] 

                for element in file_to_cal_gene_tscore_relapse.iloc[line_index, 1:last_index_to_read_file_to_cal_gene_tscore_relapse]:
                        list_gene_expression_same_probe_id.append(element)  

                gene_probe_id = file_to_cal_gene_tscore_relapse.iloc[line_index, 0]
                gene_expression_by_probe_id.append(gene_probe_id)
                gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                genes_expression_relapse_first_dataset[line_index] = gene_expression_by_probe_id
            
            # for class "non-relapse"
            genes_expression_no_relapse_first_dataset= {}
            last_index_to_read_file_to_cal_gene_tscore_no_relapse = len(col_to_read_file_cal_gene_tscore_no_relapse)
            for line_index in range(0, row_to_read_file_cal_gene_tscore):
                gene_expression_by_probe_id = []              
                list_gene_expression_same_probe_id = [] 

                for element in file_to_cal_gene_tscore_no_relapse.iloc[line_index, 1:last_index_to_read_file_to_cal_gene_tscore_no_relapse]:
                        list_gene_expression_same_probe_id.append(element)  

                gene_probe_id = file_to_cal_gene_tscore_no_relapse.iloc[line_index, 0]
                gene_expression_by_probe_id.append(gene_probe_id)
                gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                genes_expression_no_relapse_first_dataset[line_index] = gene_expression_by_probe_id
            
            print(" # Process : Calculating t-score ...")
            # prepare list of gene expression used to calculate t-score

            # for class "relapse"
            list_gene_expression_relapse_first_dataset = []
            for gene_index in range(0, len(genes_expression_relapse_first_dataset)):
                gene_name = genes_expression_relapse_first_dataset[gene_index][0]
                list_gene_expression = genes_expression_relapse_first_dataset[gene_index][1]
                list_gene_expression_relapse_first_dataset.append(list_gene_expression)
            
            # for class "non-relapse"
            list_gene_expression_no_relapse_first_dataset = []
            for gene_index in range(0, len(genes_expression_no_relapse_first_dataset)):
                gene_name = genes_expression_no_relapse_first_dataset[gene_index][0]
                list_gene_expression = genes_expression_no_relapse_first_dataset[gene_index][1]
                list_gene_expression_no_relapse_first_dataset.append(list_gene_expression)
            
            # calculate t-score of each gene
            genes_tscore = {}
            for gene_index in range(0, row_to_read_file_cal_gene_tscore):
                gene_tscore = []
                gene_name = list_gene_name[gene_index][1]
                gene_expression_relapse = list_gene_expression_relapse_first_dataset[gene_index]
                gene_expression_no_relapse = list_gene_expression_no_relapse_first_dataset[gene_index]

                tscore = stats.ttest_ind(gene_expression_relapse, gene_expression_no_relapse, equal_var = False)[0]

                gene_tscore.append(gene_name)
                gene_tscore.append(tscore)

                genes_tscore[gene_index] = gene_tscore
            
            print(" genes's t-score : " + str(genes_tscore) + "\n")

            # create pathways with their activity score
            # prepare files to be used
            file_ref = pd.read_csv(file_ref_name)
            file_to_convert = pd.read_csv(file_to_convert_first_dataset_name)
            file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

            # list all probe id
            print(" # Process : Creating a collection of pathways with their activity score ")
            list_probe_id = []
            for element in file_to_convert.loc[:, 'ID_REF']:
                list_probe_id.append(element)

            num_of_probe_id = len(list_probe_id)

            # scan probe id by its entrez id
            list_entrez_id = []
            for index in range(0, num_of_probe_id):
                entrez_id = file_ref.loc[file_ref['From'].isin([list_probe_id[index]])]
                # add to list if entrez_id is found
                if (entrez_id.empty):
                    list_entrez_id.append("none")
                else:
                    list_entrez_id.append(entrez_id.iloc[0][1]) 
            
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

                                # get t-score of this gene
                                for tscore_index in range(0, len(genes_tscore)):
                                    gene_name = genes_tscore[tscore_index][0]           
                                    if (gene_name == gene_probe_id):
                                        gene_tscore = genes_tscore[tscore_index][1]
                                        list_gene_same_entrez.append(gene_tscore)
                                

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
            
            print(" # Process : Calculating pathway activity by averaging t-score of member genes")
            list_pathway_activity = []
            for pathway_index in range(0, len(pathways)):
                list_gene_tscore_in_pathway = []
                pathway = []
                for gene_index in range(0, len(pathways[pathway_index][1])):
                    gene_tscore = pathways[pathway_index][1][gene_index][1]
                    list_gene_tscore_in_pathway.append(gene_tscore)

                # data to collect as pathway activity
                pathway_name = pathways[pathway_index][0]
                pathway_activity = calculate.mean(list_gene_tscore_in_pathway)

                pathway.append(pathway_name)
                pathway.append(pathway_activity)
                
                list_pathway_activity.append(pathway)

            # sort pathways using their pathway activity
            list_pathway_activity.sort(key = lambda x : x[1], reverse = True)

            # get list of top-rank pathway
            list_top_rank_pathway_activity = []
            for pathway_index in range(0, num_of_ranked_pathways):
                pathway = list_pathway_activity[pathway_index]
                list_top_rank_pathway_activity.append(pathway)
            
            # get list of pathway names 
            list_top_rank_pathway_name = []
            for pathway_index in range(0, len(list_top_rank_pathway_activity)):
                pathway_name = list_top_rank_pathway_activity[pathway_index][0]
                list_top_rank_pathway_name.append(pathway_name)
            
            # prepare data for feature selection
            # create list contain pathway activity of each samples to be used in sfs
            # for marker evaluation set
            # for class 'relapse'
            list_sample_relapse_pathway_activity_marker_evaluation_set = []
            for sample_index in range(0, len(list_sample_relapse_marker_evaluation)):
                list_pathway_activity = []
                sample_index_in_list = list_sample_relapse_marker_evaluation[sample_index]
                for top_ranked_pathway_index in range(0, len(list_top_rank_pathway_name)):
                    for pathway_index in range(0, len(samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1])):
                        pathway = []
                        pathway_name = samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][0]
                            
                        if (pathway_name == list_top_rank_pathway_name[top_ranked_pathway_index]):
                            pathway_activity = samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)

                            list_pathway_activity.append(pathway)

                list_sample_relapse_pathway_activity_marker_evaluation_set.append(list_pathway_activity)
            
            # for class "non-relapse"
            list_sample_no_relapse_pathway_activity_marker_evaluation_set = []
            for sample_index in range(0, len(list_sample_no_relapse_marker_evaluation)):
                    list_pathway_activity = []
                    sample_index_in_list = list_sample_no_relapse_marker_evaluation[sample_index]
                    for top_ranked_pathway_index in range(0, len(list_top_rank_pathway_name)):
                        for pathway_index in range(0, len(samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1])):
                            pathway = []
                            pathway_name = samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][0]
                                                        
                            if (pathway_name == list_top_rank_pathway_name[top_ranked_pathway_index]):
                                pathway_activity = samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][1]

                                pathway.append(pathway_name)
                                pathway.append(pathway_activity)

                                list_pathway_activity.append(pathway)

                    list_sample_no_relapse_pathway_activity_marker_evaluation_set.append(list_pathway_activity)
            
            # for feature selection set           
            # for class "relapse"
            list_sample_relapse_pathway_activity_feature_selection_set = []
            for sample_index in range(0, len(feature_selection_relapse)):
                list_pathway_activity = []
                sample_index_in_list = feature_selection_relapse[sample_index]
                for top_ranked_pathway_index in range(0, len(list_top_rank_pathway_name)):
                    for pathway_index in range(0, len(samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1])):
                        pathway = []
                        pathway_name = samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][0]
                        
                        if (pathway_name == list_top_rank_pathway_name[top_ranked_pathway_index]):
                            pathway_activity = samples_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)

                            list_pathway_activity.append(pathway)

                list_sample_relapse_pathway_activity_feature_selection_set.append(list_pathway_activity)
            
            # for class "non-relapse"
            list_sample_no_relapse_pathway_activity_feature_selection_set = []
            for sample_index in range(0, len(feature_selection_no_relapse)):
                list_pathway_activity = []
                sample_index_in_list = feature_selection_no_relapse[sample_index]
                for top_ranked_pathway_index in range(0, len(list_top_rank_pathway_name)):
                    for pathway_index in range(0, len(samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1])):
                        pathway = []
                        pathway_name = samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][0]
                        
                        if (pathway_name == list_top_rank_pathway_name[top_ranked_pathway_index]):
                            pathway_activity = samples_no_relapse_pathway_activity_first_dataset[sample_index_in_list][1][pathway_index][1]

                            pathway.append(pathway_name)
                            pathway.append(pathway_activity)

                            list_pathway_activity.append(pathway)

                list_sample_no_relapse_pathway_activity_feature_selection_set.append(list_pathway_activity)
            
            # merge both classes together
            list_sample_all_pathway_activity_feature_selection_set = []

            list_sample_all_pathway_activity_feature_selection_set.extend(list_sample_relapse_pathway_activity_feature_selection_set)
            list_sample_all_pathway_activity_feature_selection_set.extend(list_sample_no_relapse_pathway_activity_feature_selection_set)

            # get sample name of samples feature selection set
            # for class "relapse"
            list_sample_relapse_name_feature_selection = []
            for index in range(0, len(feature_selection_relapse)):
                sample_index_in_list = feature_selection_relapse[index]
                list_sample_relapse_name_feature_selection.append(samples_relapse_first_dataset[sample_index_in_list][0])
            
            # for class "non-relapse"
            list_sample_no_relapse_name_feature_selection = []
            for index in range(0, len(feature_selection_no_relapse)):
                sample_index_in_list = feature_selection_no_relapse[index]
                list_sample_no_relapse_name_feature_selection.append(samples_no_relapse_first_dataset[sample_index_in_list][0])
            
            # merge samples' name of both class
            list_sample_name_feature_selection = []
            list_sample_name_feature_selection.extend(list_sample_relapse_name_feature_selection)
            list_sample_name_feature_selection.extend(list_sample_no_relapse_name_feature_selection)

            # create list of desired output
            file_desired_outputs_feature_selection = file_output_first_dataset.loc[file_output_first_dataset['GEO asscession number'].isin(list_sample_name_feature_selection)]
            file_desired_outputs_feature_selection['sample_id'] = file_desired_outputs_feature_selection['GEO asscession number'].apply(lambda name: list_sample_name_feature_selection.index(name)) 
            file_desired_outputs_feature_selection = file_desired_outputs_feature_selection.sort_values(by = ['sample_id'])
            file_desired_outputs_feature_selection.drop(columns = 'sample_id', inplace = True)

            list_desired_outputs_feature_selection = []
            for element in file_desired_outputs_feature_selection.loc[:, 'relapse (1=True)']:
                list_desired_outputs_feature_selection.append(element)
            
            # find feature set using sequential forward selection
            feature_set_name, auc_score_feature_selection = calculate.sfsAdvance(list_top_rank_pathway_name, list_desired_outputs_feature_selection, list_sample_relapse_pathway_activity_marker_evaluation_set, \
                        list_sample_no_relapse_pathway_activity_marker_evaluation_set, list_sample_all_pathway_activity_feature_selection_set)
        
        # record feature frequency only when number of epochs is more than 1
        # (the feature selection is conducted once per epoch)
        if (num_of_epochs > 1):
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

        # conducting cross-validation on the second dataset
        # prepare data for cross-validation

        print(" # Process : Cross-validation")

        # create list of indexes used to indicate the position in the list
        # for class "relapse"
        list_index_samples_relapse_second_dataset = []
        for index in range(0, len(list_sample_relapse_second_dataset)):
            list_index_samples_relapse_second_dataset.append(index)
        
        # for class "non-relapse"
        list_index_samples_no_relapse_second_dataset = []
        for index in range(0, len(list_sample_no_relapse_second_dataset)):
            list_index_samples_no_relapse_second_dataset.append(index)

        # shuffle it to make it flexible for epoch changed
        random.shuffle(list_index_samples_relapse_second_dataset)
        random.shuffle(list_index_samples_no_relapse_second_dataset)

        # split data into k parts
        chunk_relapse_size_cv = math.ceil(len(list_index_samples_relapse_second_dataset) / num_of_folds)
        chunk_no_relapse_size_cv = math.ceil(len(list_index_samples_no_relapse_second_dataset) / num_of_folds)

        chunk_list_relapse_cv = list(calculate.chunks(list_index_samples_relapse_second_dataset, chunk_relapse_size_cv))
        print(" Number of chunks in class relapse : " + str(len(chunk_list_relapse_cv)) + "\n")

        chunk_list_no_relapse_cv = list(calculate.chunks(list_index_samples_no_relapse_second_dataset, chunk_no_relapse_size_cv))
        print(" Number of chunks in class non-relapse : " + str(len(chunk_list_no_relapse_cv)) + "\n")

        check_valid_cv, num_of_chunks_cv = calculate.checkEqualListSize(chunk_list_relapse_cv, chunk_list_no_relapse_cv)

        # create samples with their pathways and member genes
        # for class "relapse"
        samples_relapse_second_dataset = {}
        for element_index in range(0, len(list_sample_relapse_second_dataset)):
            print()
            print(" Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
            print(" " + str(len(list_sample_relapse_second_dataset) - (element_index + 1)) + " samples left")
            print()

            sample = []
            sample_name = list_sample_relapse_second_dataset[element_index]
            pathways = calculate.getPathway(file_ref_name, file_to_convert_second_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

            sample.append(sample_name)
            sample.append(pathways)
            samples_relapse_second_dataset[element_index] = sample
        
        # for class "non-relapse"
        samples_no_relapse_second_dataset = {}
        for element_index in range(0, len(list_sample_no_relapse_second_dataset)):
            print()
            print(" Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
            print(" " + str(len(list_sample_no_relapse_second_dataset) - (element_index + 1)) + " samples left")
            print()

            sample = []
            sample_name = list_sample_no_relapse_second_dataset[element_index]
            pathways = calculate.getPathway(file_ref_name, file_to_convert_second_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

            sample.append(sample_name)
            sample.append(pathways)
            samples_no_relapse_second_dataset[element_index] = sample

        print(" # Process : Creating collections of samples with their pathways' activity ...")
        # create collections of samples with their pathways
        # data will be collected in this format
        # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}} 

        # for class "relapse"
        samples_relapse_pathway_activity_second_dataset = {}
        for samples_index in range(0, len(samples_relapse_second_dataset)):
            sample = []
            list_pathway = []
            for pathway_index in range(0, len(samples_relapse_second_dataset[samples_index][1])):
                list_gene_expression_in_pathway = []
                pathway = []
                for gene_index in range(0, len(samples_relapse_second_dataset[samples_index][1][pathway_index][1])):
                    gene_expression = samples_relapse_second_dataset[samples_index][1][pathway_index][1][gene_index][1]
                    list_gene_expression_in_pathway.append(gene_expression)

                # data to collect as pathway activity
                pathway_name = samples_relapse_second_dataset[samples_index][1][pathway_index][0]
                pathway_activity = calculate.mean(list_gene_expression_in_pathway)

                pathway.append(pathway_name)
                pathway.append(pathway_activity)
                list_pathway.append(pathway)
            
            sample_name = samples_relapse_second_dataset[samples_index][0]
            sample.append(sample_name)
            sample.append(list_pathway)
            samples_relapse_pathway_activity_second_dataset[samples_index] = sample    
        
        # for class "non-relapse"
        samples_no_relapse_pathway_activity_second_dataset = {}
        for samples_index in range(0, len(samples_no_relapse_second_dataset)):
            sample = []
            list_pathway = []
            for pathway_index in range(0, len(samples_no_relapse_second_dataset[samples_index][1])):
                list_gene_expression_in_pathway = []
                pathway = []
                for gene_index in range(0, len(samples_no_relapse_second_dataset[samples_index][1][pathway_index][1])):
                    gene_expression = samples_no_relapse_second_dataset[samples_index][1][pathway_index][1][gene_index][1]
                    list_gene_expression_in_pathway.append(gene_expression)

                # data to collect as pathway activity
                pathway_name = samples_no_relapse_second_dataset[samples_index][1][pathway_index][0]
                pathway_activity = calculate.mean(list_gene_expression_in_pathway)

                pathway.append(pathway_name)
                pathway.append(pathway_activity)
                list_pathway.append(pathway)
            
            sample_name = samples_no_relapse_second_dataset[samples_index][0]
            sample.append(sample_name)
            sample.append(list_pathway)
            samples_no_relapse_pathway_activity_second_dataset[samples_index] = sample

        # list to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_auc_score = []

        # list to collect maximun AUC in each fold
        list_max_auc = []

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
                print(" Samples in class relapse used as testing set : " + str(chunk_test_relapse) + "\n")
                print()
                print(" Samples in class non-relapse used as testing set : " + str(chunk_test_no_relapse) + "\n")
                print()

                # get training set in this fold
                # for class "relapse"
                chunk_train_relapse = []
                for chunk_train_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_relapse_cv[chunk_train_relapse_index] is not chunk_test_relapse):
                        chunk_train_relapse.append(chunk_list_relapse_cv[chunk_train_relapse_index])

                # for class "non-relapse"
                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse_cv[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse_cv[chunk_train_no_relapse_index])

                # merge training data of each class
                # for class "relapse"
                list_train_relapse = []
                for i in range(0, len(chunk_train_relapse)):
                    list_train_relapse.extend(chunk_train_relapse[i])
                print(" Samples in class relapse used as training set : ")
                print(list_train_relapse)
                print()

                # for class "non-relapse"
                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])
                print(" Samples in class non-relapse used as training set : ")
                print(list_train_no_relapse)
                print()

                # preparing data for evaluation and creating classifier
                # prepare data for testing
                # for classifier class 'relapse'
                list_classifier_relapse_pathway_expression = []
                for sample_index in range(0, len(list_train_relapse)):
                    list_pathway_activity = []
                    sample_index_in_list = list_train_relapse[sample_index]
                    for feature_index in range(0, len(feature_set_name)):
                        for pathway_index in range(0, len(samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1])):

                            pathway_name = samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][0]

                            if (pathway_name == feature_set_name[feature_index]):
                                pathway_activity = samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][1]
                                list_pathway_activity.append(pathway_activity)
                                
                    list_classifier_relapse_pathway_expression.append(list_pathway_activity)
                
                # for classifier class "non-relapse"
                list_classifier_no_relapse_pathway_expression = []
                for sample_index in range(0, len(list_train_no_relapse)):
                    list_pathway_activity = []
                    sample_index_in_list = list_train_no_relapse[sample_index]
                    for feature_index in range(0, len(feature_set_name)):
                        for pathway_index in range(0, len(samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1])):
                            # pathway = []
                            pathway_name = samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][0]

                            if (pathway_name == feature_set_name[feature_index]):
                                pathway_activity = samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][1]
                                list_pathway_activity.append(pathway_activity)
                    list_classifier_no_relapse_pathway_expression.append(list_pathway_activity)
                
                # for testing set class 'relapse'
                list_testing_relapse_pathway_expression = []
                for sample_index in range(0, len(chunk_test_relapse)):
                    list_pathway_activity = []
                    sample_index_in_list = chunk_test_relapse[sample_index]
                    for feature_index in range(0, len(feature_set_name)):
                        for pathway_index in range(0, len(samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1])):
   
                            pathway_name = samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][0]

                            if (pathway_name == feature_set_name[feature_index]):
                                pathway_activity = samples_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][1]

                                list_pathway_activity.append(pathway_activity)
                    list_testing_relapse_pathway_expression.append(list_pathway_activity)
                
                # for testing set class 'non-relapse'
                list_testing_no_relapse_pathway_expression = []
                for sample_index in range(0, len(chunk_test_no_relapse)):
                    list_pathway_activity = []
                    sample_index_in_list = chunk_test_no_relapse[sample_index]
                    for feature_index in range(0, len(feature_set_name)):
                        for pathway_index in range(0, len(samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1])):

                            pathway_name = samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][0]

                            if (pathway_name == feature_set_name[feature_index]):
                                pathway_activity = samples_no_relapse_pathway_activity_second_dataset[sample_index_in_list][1][pathway_index][1]

                                list_pathway_activity.append(pathway_activity)
                    list_testing_no_relapse_pathway_expression.append(list_pathway_activity)
                
                # create list to contain all pathway activity of all sample in testing set
                list_testing_all_pathway_expression = []
                list_testing_all_pathway_expression.extend(list_testing_relapse_pathway_expression)
                list_testing_all_pathway_expression.extend(list_testing_no_relapse_pathway_expression)

                # create list to contain all samples' name of testing set
                list_chunk_test_relapse_name = []
                for index in range(0, len(chunk_test_relapse)):
                    sample_index_in_list = chunk_test_relapse[index]
                    list_chunk_test_relapse_name.append(samples_relapse_second_dataset[sample_index_in_list][0])

                list_chunk_test_no_relapse_name = []
                for index in range(0, len(chunk_test_no_relapse)):
                    sample_index_in_list = chunk_test_no_relapse[index]
                    list_chunk_test_no_relapse_name.append(samples_no_relapse_second_dataset[sample_index_in_list][0])
                
                list_samples_name_testing_all = []
                list_samples_name_testing_all.extend(list_chunk_test_relapse_name)
                list_samples_name_testing_all.extend(list_chunk_test_no_relapse_name)

                # create list of desired outputs of samples in testing set
                file_desired_outputs_testing = file_output_second_dataset.loc[file_output_second_dataset['GEO asscession number'].isin(list_samples_name_testing_all)]
                file_desired_outputs_testing['sample_id'] = file_desired_outputs_testing['GEO asscession number'].apply(lambda name: list_samples_name_testing_all.index(name)) 
                file_desired_outputs_testing = file_desired_outputs_testing.sort_values(by = ['sample_id'])
                file_desired_outputs_testing.drop(columns = 'sample_id', inplace = True)

                list_desired_outputs = []
                for element in file_desired_outputs_testing.loc[:, 'relapse (1=True)']:
                    list_desired_outputs.append(element)

                # calculate outputs using lda
                # calculate lda 
                list_actual_outputs = calculate.lda(list_testing_all_pathway_expression, list_classifier_relapse_pathway_expression, list_classifier_no_relapse_pathway_expression)

                # calculate AUC score
                auc_score = roc_auc_score(list_desired_outputs, list_actual_outputs)
                list_auc_score.append(auc_score)

                print()
                print("#### Evaluation of " + str(chunk_test_index + 1) + " - fold ####")
                print(" Feature set : " + str(feature_set_name))
                print(" Number of features in feature set : " + str(len(feature_set_name)))
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
                result_file.write("Actual outputs : " + str(list_actual_outputs) + "\n")
                result_file.write("\n")
                result_file.write("Number of desired outputs : " + str(len(list_desired_outputs)) + "\n")
                result_file.write("Desired outputs : " + str(list_desired_outputs) + "\n")
                result_file.write("AUC score from feature selection : " + str(auc_score_feature_selection) + "\n")
                result_file.write("AUC score from testing : " + str(auc_score) + "\n")
                result_file.write("\n")

                # calculate time used in this fold
                end_fold_time = time.time()
                fold_elapse_time_second = end_fold_time - start_fold_time
                fold_elapse_time_minute = fold_elapse_time_second / 60
                fold_elapse_time_minute = round(fold_elapse_time_minute, 2)

                print(" Fold elapse time : " + str(fold_elapse_time_minute) + " minutes")

                result_file.write("Fold elapse time : " + str(fold_elapse_time_minute) + " minutes \n")
                result_file.write("\n")
        
        # calculate time used in this epoch
        end_epoch_time = time.time()
        time_elapse_epoch_second = end_epoch_time - start_epoch_time
        time_elapse_epoch_minute = time_elapse_epoch_second / 60
        time_elapse_epoch_hour = time_elapse_epoch_minute / 60

        time_elapse_epoch_minute = round(time_elapse_epoch_minute, 2)
        time_elapse_epoch_hour = round(time_elapse_epoch_hour, 2)

        # collect average auc score of this epoch for further use
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
        result_file.write("Average AUC score : " + str(calculate.mean(list_auc_score)) + "\n")
        result_file.write("Maximum AUC score : " + str(auc_score_max) + "\n")
        result_file.write("Feature set which gives highest AUC score : " + "\n")
        result_file.write(str(list_feature_set_max_auc))
        result_file.write("\n")
        result_file.write("Number of features in feature set : ")
        result_file.write(str(len(list_feature_set_max_auc)))
        result_file.write("\n")
        result_file.write("Time elapse : "  + str(time_elapse_epoch_minute) + " minutes (" + str(time_elapse_epoch_hour) + " hours) ")
        result_file.write("\n")

    end_time = time.time()
    total_elapse_time_second = end_time - start_time

    total_elapse_time_minute = total_elapse_time_second / 60
    total_elapse_time_hour = total_elapse_time_minute / 60  

    total_elapse_time_minute = round(total_elapse_time_minute, 2)
    total_elapse_time_hour = round(total_elapse_time_hour, 2)

    # calculate mean over all epoch
    mean_auc_over_all_epoch = calculate.mean(list_avg_auc_each_epoch)

    print(" Average AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_auc_over_all_epoch))

    result_file.write("\nAverage AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_auc_over_all_epoch) + "\n")
    result_file.write("\n")

    if (num_of_epochs > 1):
        # rank feature frequency
        if (len(list_feature_counter) < 10):
            num_of_top_frequent_pathway = len(list_feature_counter)
        else:
            # default number of features to be shown is 10
            num_of_top_frequent_pathway = 10
        
        # rank pathway frequency in descending order
        list_feature_counter.sort(key=lambda x: x[1], reverse=True)

        # add top pathways to a list to be shown
        list_top_pathway_frequency = []
        for top_pathway_index in range(0, num_of_top_frequent_pathway):
            list_top_pathway_frequency.append(list_feature_counter[top_pathway_index])
        
        print(" Feature frequency : ")
        result_file.write("\n")
        result_file.write("Feature frequency :\n")
        for index in range(0, len(list_top_pathway_frequency)):
            feature_name = list_top_pathway_frequency[index][0]
            feature_frequency = list_top_pathway_frequency[index][1]

            print(" " + str(index + 1) + ". " + str(feature_name) + " : " + str(feature_frequency))
            result_file.write(str(index + 1) + ". " + str(feature_name) + " : " + str(feature_frequency) + "\n")

        print()
        result_file.write("\n")

    print(" Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.write("\n")

    result_file.close()   

if __name__ == "__main__":
    main()