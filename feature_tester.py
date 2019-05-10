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
    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Feature Tester")
    print(" # Purpose : This program is used to test generated feature set.")
    print(" # You have to manually configure in feature_tester.py as follows")
    print(" #   [1] list_feature : Feature to be tested")
    print(" #   [2] row_to_read_file_input : Number of rows in the file mapping between samples and their gene expression to be read")  
    print(" #   [3] file_training_input : A file contains mapping between samples and their gene expression") 
    print(" #   [4] file_training_output : A file contains mapping between samples and their healt status")  
    print(" #   [5] rows_to_read_file_pathway : Number of rows in the file mapping between pathways and their member genes to be read")  
    print(" #   [6] file_ref_name : A file mapping between gene probe id and gene entrez id")  
    print(" #   [7] file_to_convert_name : A file contains mapping between samples and their gene expression")  
    print(" #   [8] file_pathway_name : A file mapping between pathways and their member genes") 
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print(" # An output file will be created in directory 'result'.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()

    # list of feature to be tested 
    # example : epoch 7 in mean_no_normalize_10_10_10
    list_feature = ['BIOCARTA_INTRINSIC_PATHWAY', 'REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS_THAT_BIND_AU_RICH_ELEMENTS', 'PID_SMAD2_3NUCLEAR_PATHWAY', 'REACTOME_G1_S_SPECIFIC_TRANSCRIPTION', 'REACTOME_RESOLUTION_OF_AP_SITES_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY', 'KEGG_BASAL_TRANSCRIPTION_FACTORS', 'REACTOME_EXTENSION_OF_TELOMERES', 'PID_A6B1_A6B4_INTEGRIN_PATHWAY', 'REACTOME_LIPID_DIGESTION_MOBILIZATION_AND_TRANSPORT', 'REACTOME_BASE_FREE_SUGAR_PHOSPHATE_REMOVAL_VIA_THE_SINGLE_NUCLEOTIDE_REPLACEMENT_PATHWAY', 'REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT', 'PID_MET_PATHWAY', 'KEGG_SPLICEOSOME', 'BIOCARTA_TOLL_PATHWAY', 'PID_AVB3_OPN_PATHWAY', 'REACTOME_CELL_CYCLE_MITOTIC', 'REACTOME_FORMATION_OF_THE_HIV1_EARLY_ELONGATION_COMPLEX', 'REACTOME_DNA_STRAND_ELONGATION', 'REACTOME_CYCLIN_E_ASSOCIATED_EVENTS_DURING_G1_S_TRANSITION_', 'BIOCARTA_SPPA_PATHWAY', 'REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A', 'REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS', 'PID_HIF1A_PATHWAY', 'BIOCARTA_PTEN_PATHWAY', 'REACTOME_GRB2_SOS_PROVIDES_LINKAGE_TO_MAPK_SIGNALING_FOR_INTERGRINS_', 'PID_RETINOIC_ACID_PATHWAY']

    # prepare data
    # default row_to_read = 22283
    row_to_read_file_input = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read_file_input)
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

    # get output file's name
    file_name = input("Name of output file : ")

    # # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    print("Process : Creating collection to collect samples and their genes' expression")
    # create dictionary used to collect pathways of each sample
    samples_relapse = {}
    samples_no_relapse = {}

    # get all pathways of all samples in class 'relapse'
    for element_index in range(0, len(list_sample_relapse)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
        print(str(len(list_sample_relapse) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_relapse[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_relapse[element_index] = sample

    for element_index in range(0, len(list_sample_no_relapse)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
        print(str(len(list_sample_no_relapse) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_no_relapse[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_no_relapse[element_index] = sample
    
    print("Process : Creating collections of samples with their pathways' activity ...")
    
    # create collections of samples with their pathways
    # data will be collected in this format
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    samples_relapse_pathway_activity = {}
    samples_no_relapse_pathway_activity = {}
    for samples_index in range(0, len(samples_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_relapse[samples_index][1][pathway_index][1])):
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
        samples_relapse_pathway_activity[samples_index] = sample

    for samples_index in range(0, len(samples_no_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_no_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_no_relapse[samples_index][1][pathway_index][1])):
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
        samples_no_relapse_pathway_activity[samples_index] = sample
    
    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################")

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

    if (check_valid == True):
        # random index the chunk to be tested
        chunk_test_index = random.randint(0, num_of_chunks - 1)

        # separating data into testing and training dataset
        # get testing set
        chunk_test_relapse = chunk_list_relapse[chunk_test_index]
        chunk_test_no_relapse = chunk_list_no_relapse[chunk_test_index]

        # get training set of this fold
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

        list_train_no_relapse = []
        for i in range(0, len(chunk_train_no_relapse)):
            list_train_no_relapse.extend(chunk_train_no_relapse[i])
        print("size of list_train_no_relapse : " + str(len(list_train_no_relapse)))

        # making classifier
        # get pathways' activity of members in the feature set
        # for class 'relapse'
        list_sample_relapse_pathway_activity_classifier = []
        list_pathway_name_classifier_relapse = []
        for sample_index in range(0, len(list_train_relapse)):
            list_pathway_activity = []
            sample_index_in_list = list_train_relapse[sample_index]
            for feature in list_feature:
                for pathway_index in range(0, len(samples_relapse_pathway_activity[sample_index_in_list][1])):
                    pathway_name = samples_relapse_pathway_activity[sample_index_in_list][1][pathway_index][0]
                        
                    if (pathway_name == feature):
                        pathway_activity = samples_relapse_pathway_activity[sample_index_in_list][1][pathway_index][1]
                        list_pathway_activity.append(pathway_activity)
                        if(pathway_name not in list_pathway_name_classifier_relapse):
                            list_pathway_name_classifier_relapse.append(pathway_name)

            list_sample_relapse_pathway_activity_classifier.append(list_pathway_activity)
        result_file.write("feature set (" + str(len(list_feature)) + ") : \n")
        result_file.write(str(list_feature))
        result_file.write("\n")
        result_file.write("pathway name in class 'relapse' (" + str(len(list_pathway_name_classifier_relapse))+ ") : ")
        result_file.write(str(list_pathway_name_classifier_relapse))
        result_file.write("\n")

        # for class 'non-relapse'
        list_sample_no_relapse_pathway_activity_classifier = []
        list_pathway_name_classifier_no_relapse = []
        for sample_index in range(0, len(list_train_no_relapse)):
            list_pathway_activity = []
            sample_index_in_list = list_train_no_relapse[sample_index]
            for feature in list_feature:
                for pathway_index in range(0, len(samples_no_relapse_pathway_activity[sample_index_in_list][1])):
                    pathway_name = samples_no_relapse_pathway_activity[sample_index_in_list][1][pathway_index][0]
                                                
                    if (pathway_name == feature):
                        pathway_activity = samples_no_relapse_pathway_activity[sample_index_in_list][1][pathway_index][1]
                        list_pathway_activity.append(pathway_activity)
                        if(pathway_name not in list_pathway_name_classifier_no_relapse):
                            list_pathway_name_classifier_no_relapse.append(pathway_name)

            list_sample_no_relapse_pathway_activity_classifier.append(list_pathway_activity)
        result_file.write("pathway name in class 'non-relapse' (" + str(len(list_pathway_name_classifier_no_relapse)) + ") : ")
        result_file.write(str(list_pathway_name_classifier_no_relapse))
        result_file.write("\n")

        # prepare tesing set
        # each sample contains only pathway in feature set
        # for class 'relapse'
        list_sample_relapse_pathway_activity_testing_set = []
        for sample_index in range(0, len(chunk_test_relapse)):
            list_pathway_activity = []
            sample_index_in_list = chunk_test_relapse[sample_index]
            for feature in list_feature:
                for pathway_index in range(0, len(samples_relapse_pathway_activity[sample_index_in_list][1])):
                    pathway_name = samples_relapse_pathway_activity[sample_index_in_list][1][pathway_index][0]
                    
                    if (pathway_name == feature):
                        pathway_activity = samples_relapse_pathway_activity[sample_index_in_list][1][pathway_index][1]
                        list_pathway_activity.append(pathway_activity)

            list_sample_relapse_pathway_activity_testing_set.append(list_pathway_activity)
        
        # for class 'non-relapse'
        list_sample_no_relapse_pathway_activity_testing_set = []
        for sample_index in range(0, len(chunk_test_no_relapse)):
            list_pathway_activity = []
            sample_index_in_list = chunk_test_no_relapse[sample_index]
            for feature in list_feature:
                for pathway_index in range(0, len(samples_no_relapse_pathway_activity[sample_index_in_list][1])):
                    pathway_name = samples_no_relapse_pathway_activity[sample_index_in_list][1][pathway_index][0]
                    
                    if (pathway_name == feature):
                        pathway_activity = samples_no_relapse_pathway_activity[sample_index_in_list][1][pathway_index][1]
                        list_pathway_activity.append(pathway_activity)

            list_sample_no_relapse_pathway_activity_testing_set.append(list_pathway_activity)
        
        # merge testing data to be used in lda for feature selection 
        list_sample_all_pathway_activity_testing_set = []
        list_sample_all_pathway_activity_testing_set.extend(list_sample_relapse_pathway_activity_testing_set)
        list_sample_all_pathway_activity_testing_set.extend(list_sample_no_relapse_pathway_activity_testing_set)

        # get sample name of samples feature selection set
        list_sample_relapse_name_testing_set = []
        for index in range(0, len(chunk_test_relapse)):
            sample_index_in_list = chunk_test_relapse[index]
            list_sample_relapse_name_testing_set.append(samples_relapse[sample_index_in_list][0])

        list_sample_no_relapse_name_testing_set = []
        for index in range(0, len(chunk_test_no_relapse)):
            sample_index_in_list = chunk_test_no_relapse[index]
            list_sample_no_relapse_name_testing_set.append(samples_no_relapse[sample_index_in_list][0])
        
        # merge samples' name of both class
        list_sample_name_testing_set = []
        list_sample_name_testing_set.extend(list_sample_relapse_name_testing_set)
        list_sample_name_testing_set.extend(list_sample_no_relapse_name_testing_set)

        # create list of desired output
        file_desired_outputs_testing = file_training_output.loc[file_training_output['GEO asscession number'].isin(list_sample_name_testing_set)]
        file_desired_outputs_testing['sample_id'] = file_desired_outputs_testing['GEO asscession number'].apply(lambda name: list_sample_name_testing_set.index(name)) 
        file_desired_outputs_testing = file_desired_outputs_testing.sort_values(by = ['sample_id'])
        file_desired_outputs_testing.drop(columns = 'sample_id', inplace = True)

        list_desired_outputs_testing = []
        for element in file_desired_outputs_testing.loc[:, 'relapse (1=True)']:
            list_desired_outputs_testing.append(element)
        
        # linear discrimination analysis
        list_actual_outputs_testing = calculate.lda(list_sample_all_pathway_activity_testing_set, list_sample_relapse_pathway_activity_classifier, list_sample_no_relapse_pathway_activity_classifier)

        # calculate rAUC score
        auc_score = roc_auc_score(list_desired_outputs_testing, list_actual_outputs_testing)

        result_file.write("list_sample_name_testing_set (" + str(len(list_sample_name_testing_set)) + ") : " + str(list_sample_name_testing_set) +"\n")
        result_file.write("list_desired_outputs_testing (" + str(len(list_desired_outputs_testing)) + ") : \n")
        result_file.write(str(list_desired_outputs_testing))
        result_file.write("\n")
        result_file.write("list_actual_outputs_testing (" + str(len(list_actual_outputs_testing)) + ") : \n")
        result_file.write(str(list_actual_outputs_testing))
        result_file.write("\n")
        result_file.write("AUC score : " + str(auc_score) + "\n")
    result_file.close()   


if __name__ == "__main__":
    main()