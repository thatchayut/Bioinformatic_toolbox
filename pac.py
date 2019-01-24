from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
import copy
from sklearn.metrics import roc_auc_score

def main():
    # record start time
    start_time = time.time()

    # prepare data
    # row_to_read = 22283
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    rows_to_read_file_pathway = 1329
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
    file_name = input("Name of Output File : ")

    # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    # calculate number of pathways to be used
    num_of_ranked_pathways = (rows_to_read_file_pathway * (num_of_pathways_percentage / 100))
    num_of_ranked_pathways = math.ceil(num_of_ranked_pathways)

    # split data into k parts
    chunk_relapse_size = math.ceil(len(list_sample_relapse) / num_of_folds)
    chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse) / num_of_folds)

    chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_relapse_size))
    print("# chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

    chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_no_relapse_size))
    print("# chunks in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

    check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

    # list to track feature set that has the best auc score
    auc_score_max = 0
    list_feature_set_max_auc = []
    list_auc_score = []

    # do only if number of chunks of both datasets are equal
    if (check_valid == True):
        for chunk_test_index in range(0, num_of_chunks):

            strat_fold_time = time.time()
            result_file.write(" #### Fold " + str(chunk_test_index + 1) + " ####\n")
            # separating data into testing and training dataset
            chunk_test_relapse = chunk_list_relapse[chunk_test_index]
            chunk_test_no_relapse = chunk_list_no_relapse[chunk_test_index]

            print("\n------------------------------------------ K : " + str(chunk_test_index + 1) + " --------------------------------")
            print("test relapse =" + str(chunk_test_relapse))
            print("test no relapse = " + str(chunk_test_no_relapse))

            chunk_train_relapse = []
            for chunk_train_relapse_index in range(0, num_of_chunks):
                if (chunk_list_relapse[chunk_train_relapse_index] is not chunk_test_relapse):
                    chunk_train_relapse.append(chunk_list_relapse[chunk_train_relapse_index])
            print("chunk train relapse size = " + str(len(chunk_train_relapse)))
            print("chunk train relapse = " + str(chunk_train_relapse))
            
            chunk_train_no_relapse = []
            for chunk_train_no_relapse_index in range(0, num_of_chunks):
                if (chunk_list_no_relapse[chunk_train_no_relapse_index] is not chunk_test_no_relapse):
                    chunk_train_no_relapse.append(chunk_list_no_relapse[chunk_train_no_relapse_index])
            print("chunk train no relapse size = " + str(len(chunk_train_no_relapse)))
            print("chunk train no relapse = " + str(chunk_train_no_relapse))
            
            # merge training data of each class
            list_train_relapse = []
            for i in range(0, len(chunk_train_relapse)):
                list_train_relapse.extend(chunk_train_relapse[i])
            print("size of list_train_relapse : " + str(len(list_train_relapse)))
            print(list_train_relapse)

            list_train_no_relapse = []
            for i in range(0, len(chunk_train_no_relapse)):
                list_train_no_relapse.extend(chunk_train_no_relapse[i])
            print("size of list_train_no_relapse : " + str(len(list_train_no_relapse)))

            # create list of all samples in trainig data in this fold
            list_train_all_samples = []
            list_train_all_samples.extend(list_train_relapse)
            list_train_all_samples.extend(list_train_no_relapse)

            # list all gene expression to calculate mean and sd
            print("\n Gathering gene expressions are in progress ... ")
            list_train_all_gene_expression = []
            for line_index in range(0, row_to_read):
                for column in file_training_input.loc[line_index, list_train_all_samples]:
                    list_train_all_gene_expression.append(column)

            mean_all_gene_expression_train = calculate.mean(list_train_all_gene_expression)
            sd_all_gene_expression_train = calculate.sd(list_train_all_gene_expression)
            max_all_gene_expression_train = max(list_train_all_gene_expression)
            min_all_gene_expression_train = min(list_train_all_gene_expression)

            print()
            print(" Mean of all gene expression : " + str(mean_all_gene_expression_train))
            print(" SD of all gene expression : " + str(sd_all_gene_expression_train))
            print(" Max of all gene expression : " + str(max_all_gene_expression_train))
            print(" Min of all gene expression : " + str(min_all_gene_expression_train))

            # this is used to collect feature set
            list_top_ranked_pathways = []

            # this is used to collect data used in calculating lda
            list_testing_relapse_pathway_expression = []
            list_testing_no_relapse_pathway_expression = []

            check_train_valid, num_of_chunks_train = calculate.checkEqualListSize(chunk_train_relapse, chunk_train_no_relapse)

            # get pathway activity of each sample
            # get column name to be read
            if (check_train_valid is True):
                # for chunk_train_index in range(0, len(chunk_train_relapse)):
                # collection of samples containing pathways of each sample
                samples_relapse = {}
                samples_no_relapse = {}

                # identify columns to be read in each chunk in training data
                # for training data with relapse
                for element_index in range(0, len(list_train_relapse)):
                    print()
                    print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
                    print(str(len(list_train_relapse) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_train_relapse[element_index]
                    # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

                    if (method_id is "1"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                                    method = "z_score")
                    elif (method_id is "2"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                    method = "narrow_scaling")            
                    elif (method_id is "3"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                    method = "wide_scaling")    

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_relapse[element_index] = sample
                    
                    # print("Sample " + str(element_index + 1) + " name : " + str(samples_relapse[element_index][0]))
                    # print("Total Pathways : " + str(len(samples_relapse[element_index][1])))
                    # print("1st Pathway : " + str(samples_relapse[element_index][1][0]))
                    # print("1st Pathway name : " + str(samples_relapse[element_index][1][0][0]))
                    # print("Gene expression of 1st Pathway : " + str(samples_relapse[element_index][1][0][1]))
                
                print()
                print("Total number of samples relapse : " + str(len(samples_relapse)))
                print()

                for element_index in range(0, len(list_train_no_relapse)):
                    print()
                    print("Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
                    print(str(len(list_train_no_relapse) - (element_index + 1)) + " samples left")
                    print()

                    sample = []
                    sample_name = list_train_no_relapse[element_index]
                    # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

                    if (method_id is "1"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                                    method = "z_score")
                    elif (method_id is "2"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                    method = "narrow_scaling")            
                    elif (method_id is "3"):
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                    method = "wide_scaling")    

                    sample.append(sample_name)
                    sample.append(pathways)
                    samples_no_relapse[element_index] = sample
                    
                    # print("Sample " + str(element_index + 1) + " name : " + str(samples_no_relapse[element_index][0]))
                    # print("Total Pathways : " + str(len(samples_no_relapse[element_index][1])))
                    # print("1st Pathway : " + str(samples_no_relapse[element_index][1][0]))
                    # print("1st Pathway name : " + str(samples_no_relapse[element_index][1][0][0]))
                    # print("Gene expression of 1st Pathway : " + str(samples_no_relapse[element_index][1][0][1]))
                
                print()
                print("Total number of samples non-relapse : " + str(len(samples_no_relapse)))
                print()


                # version 1: this will use in finding CORGs
                # for samples_index in range(0, len(samples_no_relapse)):
                # calculate z-score for each gene expression in each pathway
                # prepare data to calculate mean and sd
                # list_gene_expression_in_pathway = []

                # get gene expression from samples with relapse
                # for samples_index in range(0, len(samples_relapse)):
                #     # assume using only 1st pahtway
                #     for gene_index in range(0, len(samples_relapse[samples_index][1][0][1])):
                #         # print(samples_relapse[samples_index][1][0][1][gene_index][1])
                #         # get gene expression in 1st pathway of this sample and add to a list 
                #         list_gene_expression_in_pathway.append(samples_relapse[samples_index][1][0][1][gene_index][1])
                
                #     # assume using only 1st pahtway
                #     for gene_index in range(0, len(samples_no_relapse[samples_index][1][0][1])):
                #         # print(samples_relapse[samples_index][1][0][1][gene_index][1])
                #         # get gene expression in 1st pathway of this sample and add to a list 
                #         list_gene_expression_in_pathway.append(samples_no_relapse[samples_index][1][0][1][gene_index][1])
                # print("Total number of genes : " + str(len(list_gene_expression_in_pathway)))

                # mean_list_gene_expression_in_pathway = calculate.mean(list_gene_expression_in_pathway)
                # sd_list_gene_expression_in_pathway = calculate.sd(list_gene_expression_in_pathway)
                
                # print("mean : " + str(mean_list_gene_expression_in_pathway))
                # print("sd : " + str(sd_list_gene_expression_in_pathway))

                # NEXT : calculate z-score

                # version 2: use 'Mean' of each pathway in each sample as a pathway activity
                samples_training_relapse_pathway_activity = {}
                # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
                for samples_index in range(0, len(samples_relapse)):
                    sample = []
                    list_pathway = []
                    for pathway_index in range(0, len(samples_relapse[samples_index][1])):
                        list_gene_expression_in_pathway = []
                        pathway = []
                        for gene_index in range(0, len(samples_relapse[samples_index][1][pathway_index][1])):
                            # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                            # print(gene_index)
                            gene_expression = samples_relapse[samples_index][1][pathway_index][1][gene_index][1]
                            list_gene_expression_in_pathway.append(gene_expression)

                        # data to collect as pathway activity
                        pathway_name = samples_relapse[samples_index][1][pathway_index][0]
                        pathway_activity = calculate.mean(list_gene_expression_in_pathway)

                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
                        list_pathway.append(pathway)
                    
                    sample_name = samples_relapse[samples_index][0]
                    sample.append(sample_name)
                    sample.append(list_pathway)
                    samples_training_relapse_pathway_activity[samples_index] = sample
                # result_file.write("samples_training_relapse_pathway_activity : \n")
                # result_file.write(" " + str(samples_training_relapse_pathway_activity) + "\n")
                # result_file.write("\n")
                
                samples_training_no_relapse_pathway_activity = {}
                # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
                for samples_index in range(0, len(samples_no_relapse)):
                    sample = []
                    list_pathway = []
                    for pathway_index in range(0, len(samples_no_relapse[samples_index][1])):
                        list_gene_expression_in_pathway = []
                        pathway = []
                        for gene_index in range(0, len(samples_no_relapse[samples_index][1][pathway_index][1])):
                            # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                            # print(gene_index)
                            gene_expression = samples_no_relapse[samples_index][1][pathway_index][1][gene_index][1]
                            list_gene_expression_in_pathway.append(gene_expression)

                        # data to collect as pathway activity
                        pathway_name = samples_no_relapse[samples_index][1][pathway_index][0]
                        pathway_activity = calculate.mean(list_gene_expression_in_pathway)

                        pathway.append(pathway_name)
                        pathway.append(pathway_activity)
                        list_pathway.append(pathway)
                    
                    sample_name = samples_no_relapse[samples_index][0]
                    sample.append(sample_name)
                    sample.append(list_pathway)
                    samples_training_no_relapse_pathway_activity[samples_index] = sample
                # result_file.write(" samples_training_no_relapse_pathway_activity : \n")
                # result_file.write(" " + str(samples_training_no_relapse_pathway_activity) + "\n")
                # result_file.write("\n")
                # print(samples_training_relapse_pathway_activity)
                # print()
                # print(samples_training_no_relapse_pathway_activity)

                # calculate t-test score for each pathway
                num_of_all_pathways = rows_to_read_file_pathway

                list_ttest_score = []
                for pathway_index in range(0, num_of_all_pathways):
                    ttest_pathway = []
                    list_pathway_activity_relapse = []
                    for samples_index in range(0, len(samples_training_relapse_pathway_activity)):
                        # print(samples_relapse[samples_index][1][pathway_index][1])
                        list_pathway_activity_relapse.append(samples_training_relapse_pathway_activity[samples_index][1][pathway_index][1])
                    list_pathway_activity_no_relapse = []
                    for samples_index in range(0, len(samples_training_no_relapse_pathway_activity)):
                        # print(samples_relapse[samples_index][1][pathway_index][1])
                        list_pathway_activity_no_relapse.append(samples_training_no_relapse_pathway_activity[samples_index][1][pathway_index][1])
                    
                    # calculate t-test score of this pathway
                    # original : abs_ttest_value = math.fabs(stats.ttest_ind(list_pathway_activity_relapse, list_pathway_activity_no_relapse, equal_var = False)[0])
                    abs_ttest_value = stats.ttest_ind(list_pathway_activity_relapse, list_pathway_activity_no_relapse, equal_var = False)[0]
                    # print(stats.ttest_ind(list_pathway_activity_relapse, list_pathway_activity_no_relapse, equal_var = False))

                    # collect data as [pathway_name, abs_ttest_value]
                    pathway_name = list_pathway_name[pathway_index]
                    ttest_pathway.append(pathway_name)
                    ttest_pathway.append(abs_ttest_value)

                    list_ttest_score.append(ttest_pathway)
                
                # sort t-test score in descending order
                list_ttest_score.sort(key = lambda x : x[1], reverse = True)
                
                print()
                print("list_ttest_score : ")
                print(list_ttest_score)
                print()

                # use number of pathways equals to user's input
                # list_top_ranked_pathway = []
                for i in range(0, num_of_ranked_pathways):
                    list_top_ranked_pathways.append(list_ttest_score[i][0][1])
                
                print("Number of ranked pathways : " + str(num_of_ranked_pathways))
                print("list_top_ranked_pathways : " + str(list_top_ranked_pathways))
                # result_file.write(" Number of ranked pathways : " + str(num_of_ranked_pathways) +"\n")
                # result_file.write(" list_top_ranked_pathways : " + str(list_top_ranked_pathway) + "\n")
                # result_file.write("\n")
                
                # prepare data for testing
                for sample_index in range(0, len(samples_training_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_training_relapse_pathway_activity[sample_index][1])):
                        pathway_name = samples_training_relapse_pathway_activity[sample_index][1][pathway_index][0]
                        # print(pathway_name)
                        pathway_activity = samples_training_relapse_pathway_activity[sample_index][1][pathway_index][1]
                        if (pathway_name in list_top_ranked_pathways):
                            # print("HEYYYYYYYYYYYYYYYYY")
                            list_pathway_activity.append(pathway_activity)
                    list_testing_relapse_pathway_expression.append(list_pathway_activity)
                
                for sample_index in range(0, len(samples_training_no_relapse_pathway_activity)):
                    list_pathway_activity = []
                    for pathway_index in range(0, len(samples_training_no_relapse_pathway_activity[sample_index][1])):
                        pathway_name = samples_training_no_relapse_pathway_activity[sample_index][1][pathway_index][0]
                        # print(pathway_name)
                        pathway_activity = samples_training_no_relapse_pathway_activity[sample_index][1][pathway_index][1]
                        if (pathway_name in list_top_ranked_pathways):
                            # print("HEYYYYYYYYYYYYYYYYY")
                            list_pathway_activity.append(pathway_activity)
                    list_testing_no_relapse_pathway_expression.append(list_pathway_activity)
            
            print("######### Testing #########")
            # use feature set in testing procedure

            # merge samples of 2 classes together
            list_samples_name_testing_all = []
            list_samples_name_testing_all.extend(chunk_test_relapse)
            list_samples_name_testing_all.extend(chunk_test_no_relapse)
            
            # get gene expression of each pathway of each sample in testing set
            samples_testing_all = {}
            for element_index in range(0, len(list_samples_name_testing_all)):
                print()
                print("Creating pathways for samples_testing_all #" + str(element_index + 1) + "  is in progress ...")
                print(str(len(list_samples_name_testing_all) - (element_index + 1)) + " samples left")
                print()

                sample = []
                sample_name = list_samples_name_testing_all[element_index]
                # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

                # Testing set uses the same mean and sd as training set since it's treated as unknown data
                if (method_id is "1"):
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                                method = "z_score")
                elif (method_id is "2"):
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                method = "narrow_scaling")            
                elif (method_id is "3"):
                    pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                                method = "wide_scaling")    

                sample.append(sample_name)
                sample.append(pathways)
                samples_testing_all[element_index] = sample 
        
            # calculate pathway activity of each pathway using 'Mean'
            samples_testing_all_pathway_activity = {}
            # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
            for samples_index in range(0, len(samples_testing_all)):
                sample = []
                list_pathway = []
                for pathway_index in range(0, len(samples_testing_all[samples_index][1])):
                    list_gene_expression_in_pathway = []
                    pathway = []
                    for gene_index in range(0, len(samples_testing_all[samples_index][1][pathway_index][1])):
                        # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                        # print(gene_index)
                        gene_expression = samples_testing_all[samples_index][1][pathway_index][1][gene_index][1]
                        list_gene_expression_in_pathway.append(gene_expression)

                    # data to collect as pathway activity
                    pathway_name = samples_testing_all[samples_index][1][pathway_index][0]
                    pathway_activity = calculate.mean(list_gene_expression_in_pathway)

                    pathway.append(pathway_name)
                    pathway.append(pathway_activity)

                    # add to list only if it is in feature set
                    if (pathway_name in list_top_ranked_pathways):
                        list_pathway.append(pathway)
                
                sample_name = samples_testing_all[samples_index][0]
                sample.append(sample_name)
                sample.append(list_pathway)
                samples_testing_all_pathway_activity[samples_index] = sample

            # write to file
            result_file.write("Mean of each pathway in each samples : \n")
            for sample_index in range(0, len(samples_testing_all_pathway_activity)):
                result_file.write(str(samples_testing_all_pathway_activity[sample_index]))
                result_file.write("\n")

            # get gene expression of each pathway of each sample in testing set
            # samples_testing_relapse = {}
            # for element_index in range(0, len(chunk_test_relapse)):
            #     print()
            #     print("Creating pathways for sample testing " + str(element_index + 1) + " relapse is in progress ...")
            #     print(str(len(chunk_test_relapse) - (element_index + 1)) + " samples left")
            #     print()

            #     sample = []
            #     sample_name = chunk_test_relapse[element_index]
            #     pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

            #     sample.append(sample_name)
            #     sample.append(pathways)
            #     samples_testing_relapse[element_index] = sample
                    
            # samples_testing_no_relapse = {}
            # for element_index in range(0, len(chunk_test_no_relapse)):
            #     print()
            #     print("Creating pathways for sample testing " + str(element_index + 1) + " relapse is in progress ...")
            #     print(str(len(chunk_test_no_relapse) - (element_index + 1)) + " samples left")
            #     print()

            #     sample = []
            #     sample_name = chunk_test_no_relapse[element_index]
            #     pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

            #     sample.append(sample_name)
            #     sample.append(pathways)
            #     samples_testing_no_relapse[element_index] = sample
            # print("samples_testing_relapse : ")
            # print(samples_testing_relapse)
            # print("samples_testing_no_relapse : ")
            # print(samples_testing_no_relapse)

            # calculate pathway activity of each pathway using 'Mean'
            # samples_testing_relapse_pathway_activity = {}
            # # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
            # for samples_index in range(0, len(samples_testing_relapse)):
            #     sample = []
            #     list_pathway = []
            #     for pathway_index in range(0, len(samples_testing_relapse[samples_index][1])):
            #         list_gene_expression_in_pathway = []
            #         pathway = []
            #         for gene_index in range(0, len(samples_testing_relapse[samples_index][1][pathway_index][1])):
            #             # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
            #             # print(gene_index)
            #             gene_expression = samples_testing_relapse[samples_index][1][pathway_index][1][gene_index][1]
            #             list_gene_expression_in_pathway.append(gene_expression)

            #         # data to collect as pathway activity
            #         pathway_name = samples_testing_relapse[samples_index][1][pathway_index][0]
            #         pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            #         pathway.append(pathway_name)

            #         # add to list only if it is in feature set
            #         if (pathway_name in list_top_ranked_pathways):
            #             list_pathway.append(pathway)
                
            #     sample_name = samples_testing_relapse[samples_index][0]
            #     sample.append(sample_name)
            #     sample.append(list_pathway)
            #     samples_testing_relapse_pathway_activity[samples_index] = sample

            # samples_testing_no_relapse_pathway_activity = {}
            # # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
            # for samples_index in range(0, len(samples_testing_no_relapse)):
            #     sample = []
            #     list_pathway = []
            #     for pathway_index in range(0, len(samples_testing_no_relapse[samples_index][1])):
            #         list_gene_expression_in_pathway = []
            #         pathway = []
            #         for gene_index in range(0, len(samples_testing_no_relapse[samples_index][1][pathway_index][1])):
            #             # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
            #             # print(gene_index)
            #             gene_expression = samples_testing_no_relapse[samples_index][1][pathway_index][1][gene_index][1]
            #             list_gene_expression_in_pathway.append(gene_expression)

            #         # data to collect as pathway activity
            #         pathway_name = samples_testing_no_relapse[samples_index][1][pathway_index][0]
            #         pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            #         pathway.append(pathway_name)
            #         pathway.append(pathway_activity)

            #         # add to list only if it is in feature set
            #         if (pathway_name in list_top_ranked_pathways):
            #             list_pathway.append(pathway)
                
            #     sample_name = samples_testing_no_relapse[samples_index][0]
            #     sample.append(sample_name)
            #     sample.append(list_pathway)
            #     samples_testing_no_relapse_pathway_activity[samples_index] = sample
            
            print("samples_testing_all_pathway_activity : ")
            print(samples_testing_all_pathway_activity)
            print("size of samples_testing_all_pathway_activity : " + str(len(samples_testing_all_pathway_activity)))

            # create list of testing dataset to be used in lda calculation
            list_testing_all_pathway_expression = []
            for sample_index in range(0, len(samples_testing_all_pathway_activity)):
                list_pathway_activity = []
                for pathway_index in range(0, len(samples_testing_all_pathway_activity[sample_index][1])):
                    pathway_name = samples_testing_all_pathway_activity[sample_index][1][pathway_index][0]
                    # print(pathway_name)
                    pathway_activity = samples_testing_all_pathway_activity[sample_index][1][pathway_index][1]
                    if (pathway_name in list_top_ranked_pathways):
                        # print("HEYYYYYYYYYYYYYYYYY")
                        list_pathway_activity.append(pathway_activity)
                list_testing_all_pathway_expression.append(list_pathway_activity)
            print()
            print("list_testing_all_pathway_expression : " + str(list_testing_all_pathway_expression))
            print()
            # print("list_testing_all_pathway_expression : ")
            # print(list_testing_all_pathway_expression)
            # print()
            print("list_testing_relapse_pathway_expression : ")
            print(list_testing_relapse_pathway_expression)
            print()
            print("list_testing_no_relapse_pathway_expression : ")
            print(list_testing_no_relapse_pathway_expression)

            # create list of desired outputs
            file_desired_outputs = file_training_output.loc[file_training_output['GEO asscession number'].isin(list_samples_name_testing_all)]
            
            list_desired_outputs = []
            for element in file_desired_outputs.loc[:, 'relapse (1=True)']:
                list_desired_outputs.append(element)

            print()

            # calculate lda 
            list_actual_outputs = calculate.lda(list_testing_all_pathway_expression, list_testing_relapse_pathway_expression, list_testing_no_relapse_pathway_expression)
           
            # calculate AUC score
            auc_score = roc_auc_score(list_desired_outputs, list_actual_outputs)
            list_auc_score.append(auc_score)

            print()
            print("Feature set : " + str(list_top_ranked_pathways))
            print("size of list_actual_outputs : " + str(len(list_actual_outputs)))
            print("list_actual_outputs : ")
            print(list_actual_outputs)
            print()
            print("size of list_desired_outputs : " + str(len(list_desired_outputs)))
            print("list_desired_outputs : ")
            print(list_desired_outputs)
            print("AUC score : " + str(auc_score))

            # track feature set which gives maximum auc score
            if (auc_score > auc_score_max):
                list_feature_set_max_auc = copy.deepcopy(list_top_ranked_pathways)
                auc_score_max = auc_score

            result_file.write("Feature set : " + str(list_top_ranked_pathways) + "\n")
            result_file.write("size of list_actual_outputs : " + str(len(list_actual_outputs)) + "\n")
            result_file.write(str(list_actual_outputs) + "\n")
            result_file.write("\n")
            result_file.write("size of list_desired_outputs : " + str(len(list_desired_outputs)) + "\n")
            result_file.write(str(list_desired_outputs) + "\n")
            result_file.write("AUC score : " + str(auc_score) + "\n")
            result_file.write("\n")
            

            end_fold_time = time.time()
            fold_elapse_time_second = end_fold_time - strat_fold_time
            fold_elapse_time_minute = fold_elapse_time_second / 60
            fold_elapse_time_minute = round(fold_elapse_time_minute, 2)
            print("fold_elapse_time : " + str(fold_elapse_time_minute) + " minutes")
            result_file.write("fold elapse time : " + str(fold_elapse_time_minute) + " minutes \n")
            result_file.write("\n")

    end_time = time.time()
    total_elapse_time_second = end_time - start_time
    total_elapse_time_minute = total_elapse_time_second / 60
    total_elapse_time_minute = round(total_elapse_time_minute, 2)
    total_elapse_time_hour = total_elapse_time_minute / 60  
    total_elapse_time_hour = round(total_elapse_time_minute / 60)
    print("#### Summary ####")
    print(" Average AUC score : " + str(calculate.mean(list_auc_score)))
    print(" Maximum AUC score : " + str(auc_score_max))
    print(" Feature set which gives highest AUC score : ")
    print(list_feature_set_max_auc)
    print()
    print(" Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")

    result_file.write("#### Summary ####\n")
    result_file.write(" Average AUC score : " + str(calculate.mean(list_auc_score)) + "\n")
    result_file.write(" Maximum AUC score : " + str(auc_score_max) + "\n")
    result_file.write(" Feature set which gives highest AUC score : " + "\n")
    result_file.write(str(list_feature_set_max_auc))
    result_file.write("\n")
    result_file.write("Total elapse time : "  + str(total_elapse_time_minute) + " minutes (" + str(total_elapse_time_hour) + " hours) ")
    result_file.close()         

if __name__ == '__main__':
    main()