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
    # list of feature to be tested 
    list_feature = ['REACTOME_PROCESSIVE_SYNTHESIS_ON_THE_LAGGING_STRAND', 'REACTOME_LAGGING_STRAND_SYNTHESIS', 'REACTOME_DNA_STRAND_ELONGATION', 'REACTOME_E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION', 'KEGG_SPLICEOSOME', 'REACTOME_CELL_CYCLE', 'REACTOME_CELL_CYCLE_MITOTIC', 'REACTOME_RORA_ACTIVATES_CIRCADIAN_EXPRESSION', 'REACTOME_COPI_MEDIATED_TRANSPORT', 'REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT', 'BIOCARTA_INTRINSIC_PATHWAY', 'REACTOME_REGULATED_PROTEOLYSIS_OF_P75NTR', 'BIOCARTA_PTDINS_PATHWAY', 'PID_EPHA2_FWD_PATHWAY', 'KEGG_RETINOL_METABOLISM', 'REACTOME_PI3K_CASCADE', 'KEGG_OXIDATIVE_PHOSPHORYLATION', 'REACTOME_MEIOTIC_RECOMBINATION', 'REACTOME_P38MAPK_EVENTS', 'BIOCARTA_LEPTIN_PATHWAY', 'REACTOME_MITOTIC_G1_G1_S_PHASES', 'PID_ILK_PATHWAY', 'KEGG_NEUROTROPHIN_SIGNALING_PATHWAY', 'REACTOME_CYTOCHROME_P450_ARRANGED_BY_SUBSTRATE_TYPE', 'REACTOME_CIRCADIAN_REPRESSION_OF_EXPRESSION_BY_REV_ERBA', 'PID_BMP_PATHWAY', 'REACTOME_SIGNALING_BY_BMP', 'REACTOME_G1_PHASE', 'PID_INTEGRIN_A9B1_PATHWAY', 'REACTOME_AKT_PHOSPHORYLATES_TARGETS_IN_THE_CYTOSOL', 'PID_ARF_3PATHWAY', 'REACTOME_SIGNALING_BY_ERBB2', 'KEGG_THYROID_CANCER', 'REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION', 'REACTOME_PRE_NOTCH_TRANSCRIPTION_AND_TRANSLATION', 'REACTOME_POTASSIUM_CHANNELS', 'REACTOME_METABOLISM_OF_CARBOHYDRATES', 'BIOCARTA_CELLCYCLE_PATHWAY', 'REACTOME_SHC_MEDIATED_SIGNALLING', 'REACTOME_TELOMERE_MAINTENANCE', 'PID_AURORA_B_PATHWAY', 'REACTOME_INSULIN_SYNTHESIS_AND_PROCESSING', 'BIOCARTA_BARRESTIN_PATHWAY', 'BIOCARTA_PITX2_PATHWAY', 'PID_P38_MK2_PATHWAY', 'REACTOME_MEMBRANE_BINDING_AND_TARGETTING_OF_GAG_PROTEINS', 'BIOCARTA_SPPA_PATHWAY', 'BIOCARTA_FIBRINOLYSIS_PATHWAY', 'REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION', 'REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION', 'REACTOME_G1_S_SPECIFIC_TRANSCRIPTION', 'REACTOME_REMOVAL_OF_THE_FLAP_INTERMEDIATE_FROM_THE_C_STRAND', 'REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE', 'REACTOME_XENOBIOTICS']

    file_pathway_to_sample_name = "pathway_map_by_mapper.csv"
    file_sample_to_class_name = "mapping_sample_to_class_full.csv"

    # prepare a file which has pathways mapped to samples
    row_to_read_file_pathway_to_sample = 1329
    file_pathway_to_sample = pd.read_csv(file_pathway_to_sample_name, nrows = row_to_read_file_pathway_to_sample)

    # prepare a file which has samples mapped to their class
    row_to_read_file_sample_to_class = 22283
    file_sample_to_class = pd.read_csv(file_sample_to_class_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # get list of pathway name
    list_pathway_name = []
    for i in range(0, row_to_read_file_pathway_to_sample):
        pathway_name = []
        pathway_name.append(i)
        pathway_name.append(file_pathway_to_sample.loc[i, "PATHWAY_NAME"])
        list_pathway_name.append(pathway_name)

    # consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_sample_to_class.loc[file_sample_to_class['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_sample_to_class.loc[file_sample_to_class['relapse (1=True)'].isin(['0'])]

    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)
    
    # shuffle data to make each chunk does not depend on sample order
    # random.shuffle(list_sample_relapse)
    # print("list_sample_relapse SIZE = " + str(len(list_sample_relapse)))
    # random.shuffle(list_sample_no_relapse)
    # print("list_sample_no_relapse SIZE = " + str(len(list_sample_no_relapse)))

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

    # get output file's name
    file_name = input("Name of output file : ")

    # prepare text file for results to be written in
    result_file = open(str(file_name) + ".txt", "w+")

    # write feature set to an output file
    result_file.write("feature set : " + str(list_feature) + "\n")
    result_file.write("\n")

    print("Process : Creating collections of samples with their pathways' activity ...")
    # create collections of samples with their pathways
    # data will be collected in this format
    # {0:  [GSM1234, [['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', 100]. [..., ...]], ... ]

    samples_relapse_pathway_activity = {}
    samples_no_relapse_pathway_activity = {}

    # for class 'relapse'
    for sample_index in range(0, len(list_sample_relapse)):
        sample = []

        sample_name = list_sample_relapse[sample_index]
        col_to_read = ["PATHWAY_NAME", sample_name]

        pathway_this_sample = pd.read_csv(file_pathway_to_sample_name, nrows = row_to_read_file_pathway_to_sample, usecols = col_to_read)
        
        # get pathway name and pathway activity
        list_pathway_name = []
        for element in pathway_this_sample.loc[:, "PATHWAY_NAME"]:
            list_pathway_name.append(element)
        
        list_pathway_activity = []
        for element in pathway_this_sample.loc[:, sample_name]:
            list_pathway_activity.append(element)
        
        list_pathway = []
        for pathway_index in range(0, row_to_read_file_pathway_to_sample):
            pathway = []
            pathway_name = list_pathway_name[pathway_index]
            pathway_activity = list_pathway_activity[pathway_index]
            
            pathway.append(pathway_name)
            pathway.append(pathway_activity)

            list_pathway.append(pathway)
        
        sample.append(sample_name)
        sample.append(list_pathway)

        samples_relapse_pathway_activity[sample_index] = sample
    
    # for class 'no-relapse'
    for sample_index in range(0, len(list_sample_no_relapse)):
        sample = []

        sample_name = list_sample_no_relapse[sample_index]
        col_to_read = ["PATHWAY_NAME", sample_name]

        pathway_this_sample = pd.read_csv(file_pathway_to_sample_name, nrows = row_to_read_file_pathway_to_sample, usecols = col_to_read)

        # get pathway name and pathway activity
        list_pathway_name = []
        for element in pathway_this_sample.loc[:, "PATHWAY_NAME"]:
            list_pathway_name.append(element)
        
        list_pathway_activity = []
        for element in pathway_this_sample.loc[:, sample_name]:
            list_pathway_activity.append(element)

        list_pathway = []
        for pathway_index in range(0, row_to_read_file_pathway_to_sample):
            pathway = []
            pathway_name = list_pathway_name[pathway_index]
            pathway_activity = list_pathway_activity[pathway_index]
            
            pathway.append(pathway_name)
            pathway.append(pathway_activity)

            list_pathway.append(pathway)
        
        sample.append(sample_name)
        sample.append(list_pathway)

        samples_no_relapse_pathway_activity[sample_index] = sample

    # list to collect auc score of each epoch
    list_auc_all_epoch = []

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################")

        # list to collect auc score of this ecpoch
        list_auc_all_fold = []

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
        # print("number of chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

        chunk_list_no_relapse = list(calculate.chunks(list_index_samples_no_relapse, chunk_no_relapse_size))
        # print("number of in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

        check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

        if (check_valid == True):
            for chunk_test_index in range(0, num_of_chunks):
                result_file.write("\n#### Fold " + str(chunk_test_index + 1) + " ####\n")

                # separating data into testing and training dataset
                # get testing set
                chunk_test_relapse = chunk_list_relapse[chunk_test_index]
                chunk_test_no_relapse = chunk_list_no_relapse[chunk_test_index]

                # get training set of this fold
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
                # print("size of list_train_relapse : " + str(len(list_train_relapse)))

                list_train_no_relapse = []
                for i in range(0, len(chunk_train_no_relapse)):
                    list_train_no_relapse.extend(chunk_train_no_relapse[i])
                # print("size of list_train_no_relapse : " + str(len(list_train_no_relapse)))

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

                                if (pathway_name not in list_pathway_name_classifier_relapse):
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

                                if (pathway_name not in list_pathway_name_classifier_no_relapse):
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
                    list_sample_relapse_name_testing_set.append(samples_relapse_pathway_activity[sample_index_in_list][0])

                list_sample_no_relapse_name_testing_set = []
                for index in range(0, len(chunk_test_no_relapse)):
                    sample_index_in_list = chunk_test_no_relapse[index]
                    list_sample_no_relapse_name_testing_set.append(samples_no_relapse_pathway_activity[sample_index_in_list][0])

                # merge samples' name of both class
                list_sample_name_testing_set = []
                list_sample_name_testing_set.extend(list_sample_relapse_name_testing_set)
                list_sample_name_testing_set.extend(list_sample_no_relapse_name_testing_set)

                # create list of desired output
                file_desired_outputs_testing = file_sample_to_class.loc[file_sample_to_class['GEO asscession number'].isin(list_sample_name_testing_set)]
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

                list_auc_all_fold.append(auc_score)

                result_file.write("list_sample_name_testing_set (" + str(len(list_sample_name_testing_set)) + ") : " + str(list_sample_name_testing_set) +"\n")
                result_file.write("list_desired_outputs_testing (" + str(len(list_desired_outputs_testing)) + ") : \n")
                result_file.write(str(list_desired_outputs_testing))
                result_file.write("\n")
                result_file.write("list_actual_outputs_testing (" + str(len(list_actual_outputs_testing)) + ") : \n")
                result_file.write(str(list_actual_outputs_testing))
                result_file.write("\n")
                result_file.write("AUC score : " + str(auc_score) + "\n")

        mean_auc_all_fold = calculate.mean(list_auc_all_fold)
        list_auc_all_epoch.append(mean_auc_all_fold)
        
        result_file.write("#### Summary of epoch " + str(epoch_count + 1) + " ####\n")
        result_file.write("AUC score of each fold : " + str(list_auc_all_fold) + "\n")
        result_file.write("Average AUC score of this epoch : " + str(mean_auc_all_fold) + "\n")
        result_file.write("\n")

    mean_auc_all_epoch = calculate.mean(list_auc_all_epoch)

    result_file.write("#### Summary ####\n")
    result_file.write("AUC score of each epoch : " + str(list_auc_all_epoch) + "\n")
    result_file.write("Average AUC score over all epoch : " + str(mean_auc_all_epoch) + "\n")

    result_file.close()   

if __name__ == "__main__":
    main()