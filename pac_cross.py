from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
import add_ons
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main():
    # record start time
    start_time = time.time()

    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Method : CORGs-Based Classification")
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

    file_name = None

    print(" # Enter required information about the first dataset ")
    print(" 1. Enter name of the file containing mapping between probes IDs and samples of the first dataset ")
    file_gene_first_dataset_name = add_ons.checkFileValid()
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
    file_output_first_dataset_name = add_ons.checkFileValid()
    print()

    print(" # Enter required information about the second dataset ")
    print(" 1. Enter name of the file containing mapping between probes IDs and samples of the second dataset ")
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

    print(" 3. Enter name of a file containing mapping between samples and their class of the second dataset")
    file_output_second_dataset_name = add_ons.checkFileValid()
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
    # for 1st dataset
    # default row_to_read_file_gene_first_dataset = 22283
    # row_to_read_file_gene_first_dataset = 22283
    # file_gene_first_dataset_name = "GSE2034-22071 (edited).csv"
    file_gene_first_dataset = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset)

    # file_output_first_dataset_name = "mapping_sample_to_class_gse2034.csv"
    file_output_first_dataset = pd.read_csv(file_output_first_dataset_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # for 2nd dataset
    # default row_to_read_file_second_dataset = 22283
    # row_to_read_file_gene_second_dataset = 22283
    # file_gene_second_dataset_name = "GSE3494_GPL96.csv"
    file_gene_second_dataset = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset) 

    # file_output_second_dataset_name = "mapping_sample_to_class_gse3494.csv"
    file_output_second_dataset = pd.read_csv(file_output_second_dataset_name, usecols = ["GEO asscession number", "relapse (1=True)"])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    # rows_to_read_file_pathway = 1329
    # file_ref_name = "accession_number_to_entrez_id.csv"
    # file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_to_convert_first_dataset_name = file_gene_first_dataset_name
    file_to_convert_second_dataset_name = file_gene_second_dataset_name
    # file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

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

    file_name = input(" # Enter name of an output file : ")

    # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################\n")

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

        print(" # Process : Feature selection")
        # conduct feature selection on 1st dataset
        # split data into 3 parts
        num_of_chunk_feature_selection = 3
        chunk_relapse_size = math.ceil(len(list_index_samples_relapse_first_dataset) / num_of_chunk_feature_selection)
        chunk_no_relapse_size = math.ceil(len(list_index_samples_no_relapse_first_dataset) / num_of_chunk_feature_selection)

        chunk_list_relapse = list(calculate.chunks(list_index_samples_relapse_first_dataset, chunk_relapse_size))
        print(" Number of chunks in class relapse : " + str(len(chunk_list_relapse)) + "\n")

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse_first_dataset, chunk_no_relapse_size))
        print(" Number of chunks in class non-relapse : " + str(len(chunk_list_no_relapse)) + "\n")

        # merge training data of each class
        # for class "relapse"
        list_relapse_first_dataset = []
        for i in range(0, len(chunk_list_relapse)):
            list_relapse_first_dataset.extend(chunk_list_relapse[i])
        # print("size of list_relapse_first_dataset : " + str(len(list_relapse_first_dataset)))
        # print("list_relapse_first_dataset : ")
        # print(list_relapse_first_dataset)
        # print()

        # for class "non-relapse"
        list_no_relapse_first_dataset = []
        for i in range(0, len(chunk_list_no_relapse)):
            list_no_relapse_first_dataset.extend(chunk_list_no_relapse[i])
        # print("size of list_no_relapse_first_dataset : " + str(len(list_no_relapse_first_dataset)))
        # print("list_no_relapse_first_dataset : ")
        # print(list_no_relapse_first_dataset)
        # print()

        # get sample name and add to a list to be used as column index
        # for class "relapse"
        list_relapse_first_dataset_name = []      
        for i in range(0, len(list_relapse_first_dataset)):
            list_relapse_first_dataset_name.append(list_sample_relapse_first_dataset[list_relapse_first_dataset[i]])
        # print(" Samples in class relapse used as marker evaluation set : " + str(list_relapse_first_dataset_name) + "\n")

        # for class "non-relapse"
        list_no_relapse_first_dataset_name = []
        for i in range(0, len(list_no_relapse_first_dataset)):
            list_no_relapse_first_dataset_name.append(list_sample_no_relapse_first_dataset[list_no_relapse_first_dataset[i]])
        # print(" Samples in class non-relapse used as marker evaluation set : " + str(list_no_relapse_first_dataset_name) + "\n")
        
        # prepare file to be used in calculating z-score
        # default row_to_read_file_to_cal_zscore = 22283
        row_to_read_file_to_cal_zscore = row_to_read_file_gene_first_dataset
        col_to_read_file_to_cal_zscore = ["ID_REF"]
        col_to_read_file_to_cal_zscore.extend(list_relapse_first_dataset_name)
        col_to_read_file_to_cal_zscore.extend(list_no_relapse_first_dataset_name)
        # print("col_to_read_file_to_cal_zscore : ")
        # print(col_to_read_file_to_cal_zscore)
        file_to_cal_zscore = pd.read_csv(file_gene_first_dataset_name, usecols = col_to_read_file_to_cal_zscore, nrows = row_to_read_file_to_cal_zscore)
        num_of_all_samples  = len(col_to_read_file_to_cal_zscore)

        # dictionary contains genes idintified by its probe id which contain all gene expression from all samples
        # {1: [gene_probe_id, [exp1,exp2, ...]]}
        genes_expression_first_dataset = {}
        last_index_to_read_file_to_cal_zscore = len(col_to_read_file_to_cal_zscore)
        for line_index in range(0, row_to_read_file_to_cal_zscore):
            gene_expression_by_probe_id = []              
            list_gene_expression_same_probe_id = [] 

            for element in file_to_cal_zscore.iloc[line_index, 1:last_index_to_read_file_to_cal_zscore]:
                    list_gene_expression_same_probe_id.append(element)  

            gene_probe_id = file_to_cal_zscore.iloc[line_index, 0]
            gene_expression_by_probe_id.append(gene_probe_id)
            gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

            genes_expression_first_dataset[line_index] = gene_expression_by_probe_id
        
        # calculate mean and sd of each gene
        list_mean_sd_gene_expression_by_probe_id_first_dataset = []
        for gene_index in range(0, len(genes_expression_first_dataset)):
            result = []
            gene_name = genes_expression_first_dataset[gene_index][0]
            mean_of_list = calculate.mean(genes_expression_first_dataset[gene_index][1])
            sd_of_list = calculate.sd(genes_expression_first_dataset[gene_index][1])
            result.append(gene_name)
            result.append(mean_of_list)
            result.append(sd_of_list)
            list_mean_sd_gene_expression_by_probe_id_first_dataset.append(result)
        
        # create samples with their pathways and member genes
        # for class "relapse"
        samples_relapse_first_dataset = {}
        for element_index in range(0, len(list_sample_relapse_first_dataset)):
            print()
            print(" Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
            print(str(len(list_sample_relapse_first_dataset) - (element_index + 1)) + " samples left")
            print()

            sample = []
            sample_name = list_sample_relapse_first_dataset[element_index]
            pathways = calculate.getPathway(file_ref_name, file_to_convert_first_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                        list_mean_sd_gene_expression_by_probe_id_first_dataset)

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
            pathways = calculate.getPathway(file_ref_name, file_to_convert_first_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, \
                        list_mean_sd_gene_expression_by_probe_id_first_dataset)

            sample.append(sample_name)
            sample.append(pathways)
            samples_no_relapse_first_dataset[element_index] = sample
        
        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        # variable to collect data from sfs
        feature_set_name = None
        list_corg_all_pathway= None
        auc_score_feature_selection = None

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
            # print("marker_evaluation_relapse size = " + str(len(marker_evaluation_relapse)))
            # print("marker_evaluation_relapse = " + str(marker_evaluation_relapse))
            # print()

            # for class "non-relapse"
            marker_evaluation_no_relapse = []
            for marker_evaluation_index in range(0, num_of_chunks):
                if (chunk_list_no_relapse[marker_evaluation_index] is not feature_selection_no_relapse):
                    marker_evaluation_no_relapse.append(chunk_list_no_relapse[marker_evaluation_index])
            print("marker_evaluation_no_relapse size : " + str(len(marker_evaluation_no_relapse)))
            print("marker_evaluation_no_relapse : " + str(marker_evaluation_no_relapse))    
            print()

            # merge all samples in the same class
            print("\n#### merge all samples in the same class to be used later ####")

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

            # create collection of samples used in feature selection
            # for class "relapse"
            samples_relapse_marker_evaluation = {}
            for sample_index in range(0, len(list_sample_relapse_marker_evaluation)):
                index_samples_relapse = list_sample_relapse_marker_evaluation[sample_index]
                samples_relapse_marker_evaluation[sample_index] = samples_relapse_first_dataset[index_samples_relapse]
            # print()
            # print("samples_relapse_marker_evaluation : ")
            # print(samples_relapse_marker_evaluation)
            # print()

            # for class "non-relapse"
            samples_no_relapse_marker_evaluation = {}
            for sample_index in range(0, len(list_sample_no_relapse_marker_evaluation)):
                index_samples_no_relapse = list_sample_no_relapse_marker_evaluation[sample_index]
                samples_no_relapse_marker_evaluation[sample_index] = samples_no_relapse_first_dataset[index_samples_no_relapse]
            # print()
            # print("samples_no_relapse_marker_evaluation : ")
            # print(samples_no_relapse_marker_evaluation)
            # print()

            print(" # Process : Finding Condition-Responsive Genes (CORGs)")
            # find CORG
            # list to collect CORG in each pathway
            list_corg_each_pathway = []

             # find CORG in each pathway 
            for pathway_index in range(0, rows_to_read_file_pathway):
                list_ttest_gene_activity = []

                list_gene_name_in_pathway = []
                list_gene_sample_relapse = []
                list_gene_sample_no_relapse = []

                # get number of genes in this pathway
                num_of_genes_in_pathway = len(samples_relapse_marker_evaluation[0][1][pathway_index][1])

                # create list of gene entrez id in this pathway
                for gene_index in range(0, num_of_genes_in_pathway):
                    gene_entrez_id = samples_relapse_marker_evaluation[0][1][pathway_index][1][gene_index][0]
                    list_gene_name_in_pathway.append(gene_entrez_id)
                
                # create list of gene expression of every samples in the same class to be used in t-test
                for gene_index in range(0, len(list_gene_name_in_pathway)):
                    list_gene_expression_with_entrez = []
                    list_gene_expression_from_sample = []
                    gene_entrez_id = list_gene_name_in_pathway[gene_index]

                    for sample_index in range(0, len(samples_relapse_marker_evaluation)):
                        gene_expression = samples_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
                        list_gene_expression_from_sample.append(gene_expression)

                    list_gene_expression_with_entrez.append(gene_entrez_id)    
                    list_gene_expression_with_entrez.append(list_gene_expression_from_sample)
                    list_gene_sample_relapse.append(list_gene_expression_with_entrez)

                for gene_index in range(0, len(list_gene_name_in_pathway)):
                    list_gene_expression_with_entrez = []
                    list_gene_expression_from_sample = []
                    gene_entrez_id = list_gene_name_in_pathway[gene_index]

                    for sample_index in range(0, len(samples_no_relapse_marker_evaluation)):
                        gene_expression = samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
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

                while (check_finish == False):
                    if (check_improve_discriminative == True):
                        max_ttest_in_consider = 0
                        list_member_gene = []

                        for gene_index in range(0, num_of_genes_in_pathway):
                            list_gene_to_consider = deepcopy(list_corg_initial)
                            gene_entrez_id = list_ttest[gene_index][0]

                            if (gene_entrez_id not in list_corg_initial):
                                list_gene_to_consider.extend([gene_entrez_id])

                                # create list of gene expression of each sample using this member genes
                                list_sample_relapse_find_discrimination = []
                                list_sample_no_relapse_find_discrimination = []

                                # for class 'relapse'
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

                                # for class 'non-relapse'
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
        
                                # calculate activity score
                                list_sample_relapse_activity_score = []
                                list_sample_no_relapse_activity_score = []

                                # for class 'relapse'
                                for sample_index in range(0, len(samples_relapse_marker_evaluation)):
                                    sum_gene_expression = 0
                                    for gene_index in range(0, len(list_sample_relapse_find_discrimination)):
                                        sum_gene_expression += list_sample_relapse_find_discrimination[gene_index][sample_index]
                                        # sum_gene_expression += list_sample_relapse_find_discrimination[sample_index][gene_index]
                                    activity_score = (sum_gene_expression / math.sqrt(len(list_sample_relapse_find_discrimination)))
                                    list_sample_relapse_activity_score.append(activity_score)

                                # for class 'non-relapse'
                                for sample_index in range(0, len(samples_no_relapse_marker_evaluation)):
                                    sum_gene_expression = 0
                                    for gene_index in range(0, len(list_sample_no_relapse_find_discrimination)):
                                        sum_gene_expression += list_sample_no_relapse_find_discrimination[gene_index][sample_index]
                                        # sum_gene_expression += list_sample_no_relapse_find_discrimination[sample_index][gene_index]
                                    activity_score = (sum_gene_expression / math.sqrt(len(list_sample_no_relapse_find_discrimination)))
                                    list_sample_no_relapse_activity_score.append(activity_score)
                                
                                # calculate ttest score of this gene set as a discriminative score
                                ttest_member_gene_set = stats.ttest_ind(list_sample_relapse_activity_score, list_sample_no_relapse_activity_score, equal_var = False)[0]

                                if (ttest_member_gene_set > max_ttest_in_consider):
                                    max_ttest_in_consider = ttest_member_gene_set
                                    list_member_gene = deepcopy(list_gene_to_consider)

                            # tracking member genes set which gives higher discriminative score
                            if (max_ttest_in_consider > max_discriminative_over_all_features):
                                max_discriminative_over_all_features = max_ttest_in_consider
                                list_corg_initial = deepcopy(list_member_gene)
                            else:
                                check_improve_discriminative = False
                    else:
                        check_finish = True 

                list_corg_each_pathway.append(list_corg_initial) 

            print(" CORGs in each pathway : ")
            print(list_corg_each_pathway)
            print()

            # copy list_corg_each_pathway to be used in feature selection process
            list_corg_all_pathway = deepcopy(list_corg_each_pathway)

            # create list contain pathway activity of each samples to be used in sfs
            # create samples which contain only corg in each pathway
            samples_relapse_marker_evaluation_pathway_activity = {}
            samples_no_relapse_marker_evaluation_pathway_activity = {}

            # for class 'relapse'
            for sample_index in range(0, len(samples_relapse_marker_evaluation)):
                list_sample_with_pathway_activity = []
                list_pathway_activity = []
                for pathway_index in range(0, len(samples_relapse_marker_evaluation[sample_index][1])):
                    pathway = []
                    sum_gene_expression = 0
                    num_of_corg = len(list_corg_each_pathway[pathway_index])
                    for gene_index in range(0, len(samples_relapse_marker_evaluation[sample_index][1][pathway_index])):
                        gene_entrez_id = samples_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][0]
                        gene_expression =  samples_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
                        if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                            sum_gene_expression += gene_expression
                    pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                    pathway_name = samples_relapse_marker_evaluation[sample_index][1][pathway_index][0]
                    pathway.append(pathway_name)
                    pathway.append(pathway_activity)
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
                    sum_gene_expression = 0
                    num_of_corg = len(list_corg_each_pathway[pathway_index])
                    for gene_index in range(0, len(samples_no_relapse_marker_evaluation[sample_index][1][pathway_index])):
                        gene_entrez_id = samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][0]
                        gene_expression =  samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][1][gene_index][1]
                        if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                            sum_gene_expression += gene_expression
                    pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                    pathway_name = samples_no_relapse_marker_evaluation[sample_index][1][pathway_index][0]
                    pathway.append(pathway_name)
                    pathway.append(pathway_activity)
                    list_pathway_activity.append(pathway)
                
                sample_name = samples_no_relapse_marker_evaluation[sample_index][0]

                list_sample_with_pathway_activity.append(sample_name)
                list_sample_with_pathway_activity.append(list_pathway_activity)
                
                samples_no_relapse_marker_evaluation_pathway_activity[sample_index] = list_sample_with_pathway_activity
            
            # create samples for feature selection
            # for class 'relapse'
            samples_relapse_feature_selection = {}
            for sample_index in range(0, len(feature_selection_relapse)):
                index_samples_relapse = feature_selection_relapse[sample_index]
                samples_relapse_feature_selection[sample_index] = samples_relapse_first_dataset[index_samples_relapse]
            
            # for class 'non-relapse'
            samples_no_relapse_feature_selection = {}
            for sample_index in range(0, len(feature_selection_no_relapse)):
                index_samples_no_relapse = feature_selection_no_relapse[sample_index]
                samples_no_relapse_feature_selection[sample_index] = samples_no_relapse_first_dataset[index_samples_no_relapse]
            
            # calculate pathway activity 
            
            # for class 'relapse'
            samples_relapse_feature_selection_pathway_activity = {}
            for sample_index in range(0, len(samples_relapse_feature_selection)):
                list_sample_with_pathway_activity = []
                list_pathway_activity = []
                for pathway_index in range(0, len(samples_relapse_feature_selection[sample_index][1])):
                    pathway = []
                    sum_gene_expression = 0
                    num_of_corg = len(list_corg_each_pathway[pathway_index])
                    for gene_index in range(0, len(samples_relapse_feature_selection[sample_index][1][pathway_index])):
                        gene_entrez_id = samples_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][0]
                        gene_expression =  samples_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                        if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                            sum_gene_expression += gene_expression
                    pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                    pathway_name = samples_relapse_feature_selection[sample_index][1][pathway_index][0]
                    pathway.append(pathway_name)
                    pathway.append(pathway_activity)
                    list_pathway_activity.append(pathway)
                
                sample_name = samples_relapse_feature_selection[sample_index][0]

                list_sample_with_pathway_activity.append(sample_name)
                list_sample_with_pathway_activity.append(list_pathway_activity)
                
                samples_relapse_feature_selection_pathway_activity[sample_index] = list_sample_with_pathway_activity
            
            # for class 'non-relapse'
            samples_no_relapse_feature_selection_pathway_activity = {}
            for sample_index in range(0, len(samples_no_relapse_feature_selection)):
                list_sample_with_pathway_activity = []
                list_pathway_activity = []
                for pathway_index in range(0, len(samples_no_relapse_feature_selection[sample_index][1])):
                    pathway = []
                    sum_gene_expression = 0
                    num_of_corg = len(list_corg_each_pathway[pathway_index])
                    for gene_index in range(0, len(samples_no_relapse_feature_selection[sample_index][1][pathway_index])):
                        gene_entrez_id = samples_no_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][0]
                        gene_expression =  samples_no_relapse_feature_selection[sample_index][1][pathway_index][1][gene_index][1]
                        if (gene_entrez_id in list_corg_each_pathway[pathway_index]):
                            sum_gene_expression += gene_expression
                    pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                    pathway_name = samples_no_relapse_feature_selection[sample_index][1][pathway_index][0]
                    pathway.append(pathway_name)
                    pathway.append(pathway_activity)
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
                # pvalue = stats.ttest_ind(list_sample_relapse_activity_score, list_sample_no_relapse_activity_score, equal_var = False)[1]

                pvalue = stats.ttest_ind(list_relapse_pathway_activity_for_pvalue[pathway_index], list_no_relapse_pathway_activity_for_pvalue[pathway_index], equal_var = False)[1]

                pathway_pvalue.append(pathway_name)
                pathway_pvalue.append(pvalue)
                list_pvalue_pathway_activity.append(pathway_pvalue)
            
            # sort pathway using p-value in ascending order
            list_pvalue_pathway_activity.sort(key = lambda x : x[1], reverse = False)

            # reorder pathway in each sample        

            # for marker evaluation class "relapse"
            samples_relapse_marker_evaluation_pathway_activity_sorted = {}
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
            samples_no_relapse_marker_evaluation_pathway_activity_sorted ={}
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
            samples_relapse_feature_selection_pathway_activity_sorted = {}
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
            samples_no_relapse_feature_selection_pathway_activity_sorted = {}
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
                feature_selection_relapse_name.append(samples_relapse_first_dataset[index_samples_relapse][0])
            
            # for class 'non-relapse'
            feature_selection_no_relapse_name = []
            for index in range(0, len(feature_selection_no_relapse)):
                index_samples_no_relapse = feature_selection_no_relapse[index]
                feature_selection_no_relapse_name.append(samples_no_relapse_first_dataset[index_samples_no_relapse][0])
            
            # merge testing data to be used in lda for feature selection 
            feature_selection_all_name = []
            feature_selection_all_name.extend(feature_selection_relapse_name)
            feature_selection_all_name.extend(feature_selection_no_relapse_name)

            # create list of desired output
            file_desired_outputs_feature_selection = file_output_first_dataset.loc[file_output_first_dataset['GEO asscession number'].isin(feature_selection_all_name)]
            file_desired_outputs_feature_selection['sample_id'] = file_desired_outputs_feature_selection['GEO asscession number'].apply(lambda name: feature_selection_all_name.index(name)) 
            file_desired_outputs_feature_selection = file_desired_outputs_feature_selection.sort_values(by = ['sample_id'])
            file_desired_outputs_feature_selection.drop(columns = 'sample_id', inplace = True)

            list_desired_outputs_feature_selection = []
            for element in file_desired_outputs_feature_selection.loc[:, 'relapse (1=True)']:
                list_desired_outputs_feature_selection.append(element)
            # print("list_desired_outputs_feature_selection : " + str(list_desired_outputs_feature_selection))
            # print()

            # create list of pathway name in the same order as in each sample
            list_pathway_name_feature_selection = []
            for pathway_index in range(0, len(samples_relapse_feature_selection_pathway_activity_sorted[0][1])):
                pathway_name = samples_relapse_feature_selection_pathway_activity_sorted[0][1][pathway_index][0]
                list_pathway_name_feature_selection.append(pathway_name)
            
            # create list contain pathway activity of each samples to be used in sfs
            # for marker evaluation set
            # for class "relapse"
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
            
            # print("list_sample_relapse_pathway_expression_marker_evaluation : ")
            # print(list_sample_relapse_pathway_expression_marker_evaluation)
            # print()

            # for class "non-relapse"
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
            
            # print("list_sample_no_relapse_pathway_expression_marker_evaluation : ")
            # print(list_sample_no_relapse_pathway_expression_marker_evaluation)
            # print()

            # for feature selection set
            # for class "relapse"
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
            
            # print("list_sample_relapse_pathway_expression_feature_selection : ")
            # print(list_sample_relapse_pathway_expression_feature_selection)
            # print()

            # for class "non-relapse"
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
            
            # print("list_sample_no_relapse_pathway_expression_feature_selection : ")
            # print(list_sample_no_relapse_pathway_expression_feature_selection)
            # print()

            # merge testing set for feature selection together
            list_sample_all_pathway_expression_feature_selection = []
            list_sample_all_pathway_expression_feature_selection.extend(list_sample_relapse_pathway_expression_feature_selection)
            list_sample_all_pathway_expression_feature_selection.extend(list_sample_no_relapse_pathway_expression_feature_selection)
            # print("list_sample_all_pathway_expression_feature_selection size : " + str(len(list_sample_all_pathway_expression_feature_selection)))
            # print(list_sample_all_pathway_expression_feature_selection)
            # print()

            print(" # Process : Sequential Forward Selection (SFS)")

            # find feature set using sequential forward selection
            feature_set_name, auc_score_feature_selection = calculate.sfsAdvance(list_pathway_name_feature_selection, list_desired_outputs_feature_selection, list_sample_relapse_pathway_expression_marker_evaluation, \
                    list_sample_no_relapse_pathway_expression_marker_evaluation, list_sample_all_pathway_expression_feature_selection)

            # list_max_auc.append(auc_score_feature_selection)

            print(" Feature set : " + str(feature_set_name))
            print(" AUC score from feature selection : " + str(auc_score_feature_selection))
            print()

            result_file.write("Feature set : " + str(feature_set_name))
            result_file.write("\n")
            result_file.write("AUC score from feature selection : " + str(auc_score_feature_selection))
            result_file.write("\n")
            # print("\n-------------------------------------------------------------------------------------------------------------\n")

        # conducting cross-validation on the second dataset
        # prepare data for cross-validation

        print(" Process : Cross-validation")

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
        print(" Number of chunks in class non-relapse  = " + str(len(chunk_list_no_relapse_cv)) + "\n")

        check_valid_cv, num_of_chunks_cv = calculate.checkEqualListSize(chunk_list_relapse_cv, chunk_list_no_relapse_cv)

        # list and variable to track feature set that has the best auc score
        auc_score_max = 0
        list_feature_set_max_auc = []
        list_corg_in_feature_set_max_auc = []

        # list to collect auc score from testing which will be used to calculate an average auc score
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
                print(" Samples in class relapse used as testing set : " + str(chunk_test_relapse))
                print()
                print(" Samples in class non-relapse used as testing set : " + str(chunk_test_no_relapse))
                print()

                # get training set in this fold
                # for class "relapse"
                chunk_train_relapse = []
                for chunk_train_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_relapse_cv[chunk_train_relapse_index] is not chunk_test_relapse):
                        chunk_train_relapse.append(chunk_list_relapse_cv[chunk_train_relapse_index])
                # print("chunk train relapse size = " + str(len(chunk_train_relapse)))

                # for class "non-relapse"
                chunk_train_no_relapse = []
                for chunk_train_no_relapse_index in range(0, num_of_chunks):
                    if (chunk_list_no_relapse_cv[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                        chunk_train_no_relapse.append(chunk_list_no_relapse_cv[chunk_train_no_relapse_index])
                # print("chunk train no relapse size = " + str(len(chunk_train_no_relapse)))

                # merge training data of each class
                # for class "relapse"
                list_train_relapse = []
                for i in range(0, len(chunk_train_relapse)):
                    list_train_relapse.extend(chunk_train_relapse[i])
                # print("size of list_train_relapse : " + str(len(list_train_relapse)))
                # print("list_train_relapse : ")
                # print(list_train_relapse)
                # print()

                # for class "non-relapse"
                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])
                # print("size of list_train_no_relapse : " + str(len(list_train_no_relapse)))
                # print("list_train_no_relapse : ")
                # print(list_train_no_relapse)
                # print()

                # get sample name and add to a list to be used as column index
                # for class "relapse"
                list_train_relapse_name = []
                for i in range(0, len(list_train_relapse)):
                    list_train_relapse_name.append(list_sample_relapse_second_dataset[list_train_relapse[i]])
                print(" Samples in class relapse used as training set : " + str(list_train_relapse_name) + "\n")

                # for class "non-relapse"
                list_train_no_relapse_name = []
                for i in range(0, len(list_train_no_relapse)):
                    list_train_no_relapse_name.append(list_sample_no_relapse_second_dataset[list_train_no_relapse[i]])
                print(" Samples in class non-relapse used as training set : " + str(list_train_no_relapse_name) + "\n")

                # default row_to_read_file_to_cal_zscore = 22283
                row_to_read_file_to_cal_zscore_cv = row_to_read_file_gene_second_dataset

                col_to_read_file_to_cal_zscore_cv = ["ID_REF"]
                col_to_read_file_to_cal_zscore_cv.extend(list_train_relapse_name)
                col_to_read_file_to_cal_zscore_cv.extend(list_train_no_relapse_name)
                print("col_to_read_file_to_cal_zscore : ")
                print(col_to_read_file_to_cal_zscore_cv)
                
                print(file_gene_second_dataset_name)
                file_to_cal_zscore_cv = pd.read_csv(file_gene_second_dataset_name, usecols = col_to_read_file_to_cal_zscore_cv, nrows = row_to_read_file_to_cal_zscore_cv)
                num_of_all_samples_cv  = len(col_to_read_file_to_cal_zscore_cv)

                # dictionary contains genes idintified by its probe id which contain all gene expression from all samples
                # {1: [gene_probe_id, [exp1,exp2, ...]]}
                genes_expression_second_dataset = {}
                last_index_to_read_file_to_cal_zscore_cv = len(col_to_read_file_to_cal_zscore_cv)
                for line_index in range(0, row_to_read_file_to_cal_zscore_cv):
                    gene_expression_by_probe_id = []              
                    list_gene_expression_same_probe_id = [] 

                    for element in file_to_cal_zscore_cv.iloc[line_index, 1:last_index_to_read_file_to_cal_zscore_cv]:
                         list_gene_expression_same_probe_id.append(element)  

                    gene_probe_id = file_to_cal_zscore_cv.iloc[line_index, 0]
                    gene_expression_by_probe_id.append(gene_probe_id)
                    gene_expression_by_probe_id.append(list_gene_expression_same_probe_id)

                    genes_expression_second_dataset[line_index] = gene_expression_by_probe_id
                
                # calculate mean and sd of each gene
                list_mean_sd_gene_expression_by_probe_id_cv = []
                for gene_index in range(0, len(genes_expression_second_dataset)):
                    result = []
                    gene_name = genes_expression_second_dataset[gene_index][0]
                    mean_of_list = calculate.mean(genes_expression_second_dataset[gene_index][1])
                    sd_of_list = calculate.sd(genes_expression_second_dataset[gene_index][1])
                    result.append(gene_name)
                    result.append(mean_of_list)
                    result.append(sd_of_list)
                    list_mean_sd_gene_expression_by_probe_id_cv.append(result)
                
                # create samples with their pathways and member genes
                # for class "relapse"
                samples_relapse_second_dataset = {}
                for element_index in range(0, len(list_sample_relapse_second_dataset)):
                    print()
                    print(" Creating pathways for sample " + str(element_index + 1) + " relapse is in progress for cross-validation...")
                    print(str(len(list_sample_relapse_second_dataset) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_sample_relapse_second_dataset[element_index]
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_second_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway,\
                                list_mean_sd_gene_expression_by_probe_id_cv)

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_relapse_second_dataset[element_index] = sample    
                
                # for class "non-relapse"
                samples_no_relapse_second_dataset = {}
                for element_index in range(0, len(list_sample_no_relapse_second_dataset)):
                    print()
                    print(" Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress for cross-validation...")
                    print(str(len(list_sample_no_relapse_second_dataset) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_sample_no_relapse_second_dataset[element_index]
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_second_dataset_name, file_pathway_name, sample_name, rows_to_read_file_pathway, \
                                list_mean_sd_gene_expression_by_probe_id_cv)

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_no_relapse_second_dataset[element_index] = sample 
                
                # preparing data for evaluation and creating classifier
                # create classifier

                # for class 'relapse'
                samples_classifier_relapse = {}
                for sample_index in range(0,len(list_train_relapse)):
                    index_samples_relapse = list_train_relapse[sample_index]
                    samples_classifier_relapse[sample_index] = samples_relapse_second_dataset[index_samples_relapse]
                
                # for class 'non-relapse'
                samples_classifier_no_relapse = {}
                for sample_index in range(0,len(list_train_no_relapse)):
                    index_samples_no_relapse = list_train_no_relapse[sample_index]
                    samples_classifier_no_relapse[sample_index] = samples_no_relapse_second_dataset[index_samples_no_relapse]
                
                # calculate pathway activity of classifier
                # for class 'relapse'
                samples_classifier_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_classifier_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_classifier_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_expression = 0
                        num_of_corg = len(list_corg_all_pathway[pathway_index])
                        for gene_index in range(0, len(samples_classifier_relapse[sample_index][1][pathway_index])):
                            gene_entrez_id = samples_classifier_relapse[sample_index][1][pathway_index][1][gene_index][0]
                            gene_expression =  samples_classifier_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            if (gene_entrez_id in list_corg_all_pathway[pathway_index]):
                                sum_gene_expression += gene_expression
                        pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                        pathway_name = samples_classifier_relapse[sample_index][1][pathway_index][0]
                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
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
                        sum_gene_expression = 0
                        num_of_corg = len(list_corg_all_pathway[pathway_index])
                        for gene_index in range(0, len(samples_classifier_no_relapse[sample_index][1][pathway_index])):
                            gene_entrez_id = samples_classifier_no_relapse[sample_index][1][pathway_index][1][gene_index][0]
                            gene_expression =  samples_classifier_no_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            if (gene_entrez_id in list_corg_all_pathway[pathway_index]):
                                sum_gene_expression += gene_expression
                        pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                        pathway_name = samples_classifier_no_relapse[sample_index][1][pathway_index][0]
                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
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
                    samples_testing_relapse[sample_index] = samples_relapse_second_dataset[index_samples_relapse]
                
                # for class 'non-relapse'
                samples_testing_no_relapse = {}
                for sample_index in range(0,len(chunk_test_no_relapse)):
                    index_samples_no_relapse = chunk_test_no_relapse[sample_index]
                    samples_testing_no_relapse[sample_index] = samples_no_relapse_second_dataset[index_samples_no_relapse]
                
                # calculate pathway activity
                # for class 'relapse'
                samples_testing_relapse_pathway_activity = {}
                for sample_index in range(0, len(samples_testing_relapse)):
                    list_sample_with_pathway_activity = []
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_testing_relapse[sample_index][1])):
                        pathway = []
                        sum_gene_expression = 0
                        num_of_corg = len(list_corg_all_pathway[pathway_index])
                        for gene_index in range(0, len(samples_testing_relapse[sample_index][1][pathway_index])):
                            gene_entrez_id = samples_testing_relapse[sample_index][1][pathway_index][1][gene_index][0]
                            gene_expression =  samples_testing_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            if (gene_entrez_id in list_corg_all_pathway[pathway_index]):
                                sum_gene_expression += gene_expression
                        pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                        pathway_name = samples_testing_relapse[sample_index][1][pathway_index][0]
                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
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
                        sum_gene_expression = 0
                        num_of_corg = len(list_corg_all_pathway[pathway_index])
                        for gene_index in range(0, len(samples_testing_no_relapse[sample_index][1][pathway_index])):
                            gene_entrez_id = samples_testing_no_relapse[sample_index][1][pathway_index][1][gene_index][0]
                            gene_expression =  samples_testing_no_relapse[sample_index][1][pathway_index][1][gene_index][1]
                            if (gene_entrez_id in list_corg_all_pathway[pathway_index]):
                                sum_gene_expression += gene_expression
                        pathway_activity = (sum_gene_expression / math.sqrt(num_of_corg))
                        pathway_name = samples_testing_no_relapse[sample_index][1][pathway_index][0]
                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
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
                # print("list_classifier_relapse_pathway_expression : ")
                # print(list_classifier_relapse_pathway_expression)
                # print() 

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
                # print("list_classifier_no_relapse_pathway_expression : ")
                # print(list_classifier_no_relapse_pathway_expression)
                # print()   

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
                # print("list_testing_relapse_pathway_expression : ")
                # print(list_testing_relapse_pathway_expression)
                # print()     
                
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
                # print("list_testing_no_relapse_pathway_expression : ")
                # print(list_testing_no_relapse_pathway_expression)
                # print() 

                # merge testing data of 2 class together
                list_testing_all_pathway_expression = []
                list_testing_all_pathway_expression.extend(list_testing_relapse_pathway_expression)
                list_testing_all_pathway_expression.extend(list_testing_no_relapse_pathway_expression)

                # create list to contain all sample name used in testing procedure
                # for class "relapse"
                list_chunk_test_relapse_name = []
                for index in range(0, len(chunk_test_relapse)):
                    index_samples_relapse = chunk_test_relapse[index]
                    list_chunk_test_relapse_name.append(samples_relapse_second_dataset[index_samples_relapse][0])
                print(" Samples in class relapse used as testing set : ")
                print(list_chunk_test_relapse_name)
                print()

                # for class "non-relapse"
                list_chunk_test_no_relapse_name = []
                for index in range(0, len(chunk_test_no_relapse)):
                        index_samples_no_relapse = chunk_test_no_relapse[index]
                        list_chunk_test_no_relapse_name.append(samples_no_relapse_second_dataset[index_samples_no_relapse][0])
                print(" Samples in class non-relapse used as testing site : ")
                print(list_chunk_test_no_relapse_name)
                print()

                list_samples_name_testing_all = []
                list_samples_name_testing_all.extend(list_chunk_test_relapse_name)
                list_samples_name_testing_all.extend(list_chunk_test_no_relapse_name)
                # print("list_samples_name_testing_all : ")
                # print(list_samples_name_testing_all)
                # print()

                # create list of desired outputs
                file_desired_outputs = file_output_second_dataset.loc[file_output_second_dataset['GEO asscession number'].isin(list_samples_name_testing_all)]

                file_desired_outputs['sample_id'] = file_desired_outputs['GEO asscession number'].apply(lambda name: list_samples_name_testing_all.index(name)) 
                file_desired_outputs = file_desired_outputs.sort_values(by = ['sample_id'])
                file_desired_outputs.drop(columns = 'sample_id', inplace = True)

                list_desired_outputs = []
                for element in file_desired_outputs.loc[:, 'relapse (1=True)']:
                    list_desired_outputs.append(element)
                
                # calculate lda 
                list_actual_outputs = calculate.lda(list_testing_all_pathway_expression, list_classifier_relapse_pathway_expression, list_classifier_no_relapse_pathway_expression)

                # calculate AUC score
                auc_score = roc_auc_score(list_desired_outputs, list_actual_outputs)
                list_auc_score.append(auc_score)

                # get corg of each feature in feature set
                # get index of members in feature set in relation to th eorder of pathway in an input file
                list_index_feature_set = []
                for feature_index in range(0, len(feature_set_name)):
                    for pathway_index in range(0, len(list_pathway_name)):
                        feature = feature_set_name[feature_index]
                        pathway = list_pathway_name[pathway_index][1]

                        if (feature == pathway):
                            index_feature_in_list_pathway = list_pathway_name[pathway_index][0]
                            list_index_feature_set.append(index_feature_in_list_pathway)

                # map index of members in feature set to index of CORG 
                list_corg_feature_set = []
                for feature_index in range(0, len(list_index_feature_set)):
                    index_in_list = list_index_feature_set[feature_index]

                    corg_this_feature = list_corg_all_pathway[index_in_list]
                    list_corg_feature_set.append(corg_this_feature)
                
                print()
                print("#### Evaluation of " + str(chunk_test_index + 1) + " - fold ####")
                # print("Feature set : " + str(list_top_ranked_pathways))
                print(" Feature set : " + str(feature_set_name))
                print(" CORG of each feature in feature set : ")
                print(list_corg_feature_set)
                print(" Number of features in feature set : " + str(len(feature_set_name)))
                print(" Number of actual outputs : " + str(len(list_actual_outputs)))
                print(" Actual outputs : ")
                print(list_actual_outputs)
                print()
                print(" Number of desired output : " + str(len(list_desired_outputs)))
                print(" Desired outputs : ")
                print(list_desired_outputs)
                print(" AUC score from feature selection : " + str(auc_score_feature_selection))
                print(" AUC score from testing : " + str(auc_score))

                 # track feature set which gives maximum auc score
                if (auc_score > auc_score_max):
                    list_feature_set_max_auc = deepcopy(feature_set_name)
                    list_corg_in_feature_set_max_auc = deepcopy(list_corg_feature_set)
                    auc_score_max = auc_score
                
                result_file.write("Feature set : " + str(feature_set_name) + "\n")
                result_file.write("CORG of each feature in feature set : \n")
                result_file.write(str(list_corg_feature_set) + "\n")
                result_file.write("Number of features in feature set : " + str(len(feature_set_name)) + "\n")
                result_file.write("Number of actual outputs : " + str(len(list_actual_outputs)) + "\n")
                result_file.write(str(list_actual_outputs) + "\n")
                result_file.write("\n")
                result_file.write("Number of desired outputs : " + str(len(list_desired_outputs)) + "\n")
                result_file.write(str(list_desired_outputs) + "\n")
                result_file.write(" AUC score from feature selection : " + str(auc_score_feature_selection) + "\n")
                result_file.write(" AUC score from testing : " + str(auc_score) + "\n")
                result_file.write("\n")
                

                end_fold_time = time.time()
                fold_elapse_time_second = end_fold_time - start_fold_time
                fold_elapse_time_minute = fold_elapse_time_second / 60
                fold_elapse_time_minute = round(fold_elapse_time_minute, 2)
                print(" Fold elapse time : " + str(fold_elapse_time_minute) + " minutes")
                result_file.write("Fold elapse time : " + str(fold_elapse_time_minute) + " minutes \n")
                result_file.write("\n")

        end_time = time.time()
        total_elapse_time_second = end_time - start_time
        total_elapse_time_minute = total_elapse_time_second / 60
        total_elapse_time_minute = round(total_elapse_time_minute, 2)
        total_elapse_time_hour = total_elapse_time_minute / 60  
        total_elapse_time_hour = round(total_elapse_time_minute / 60)

        list_avg_auc_each_epoch.append(calculate.mean(list_auc_score))

        print()
        print("#### Summary ####")
        print(" Average AUC score : " + str(calculate.mean(list_auc_score)))
        print(" Maximum AUC score : " + str(auc_score_max))
        print(" Feature set which gives highest AUC score : ")
        print(list_feature_set_max_auc)
        print(" CORG of each feature in the feature set which gives the highest AUC score : ")
        print(list_corg_in_feature_set_max_auc)
        print()
        print(" Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")

        result_file.write("\n#### Summary ####\n")

        # result_file.write("Maximum AUC ROC score of feature from feature selection in each fold : \n")
        # result_file.write(str(list_max_auc))
        # result_file.write("\n")

        result_file.write("Average AUC score : " + str(calculate.mean(list_auc_score)) + "\n")
        result_file.write("Maximum AUC score : " + str(auc_score_max) + "\n")
        result_file.write("Size of feature set which gives the highest AUC score from testing : " + str(len(list_feature_set_max_auc)))
        result_file.write("\n")
        result_file.write("Feature set which gives the highest AUC score from testing : " + "\n")
        result_file.write(str(list_feature_set_max_auc))
        result_file.write("\n")
        result_file.write("CORG of each feature in the feature set which gives the highest AUC score : \n")
        result_file.write(str(list_corg_in_feature_set_max_auc))
        result_file.write("\n")
        result_file.write("Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
        result_file.write("\n")
        result_file.write("\n")

        print("----------------------------------------------------------------------------------------------------")
    
    # calculate mean over all epoch
    mean_over_all_epoch = calculate.mean(list_avg_auc_each_epoch)
    print(" Average AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_over_all_epoch))
    result_file.write("Average AUC score over " + str(num_of_epochs) + " epoch : " + str(mean_over_all_epoch) + "\n")

    result_file.close()




if __name__ == "__main__":
    main()