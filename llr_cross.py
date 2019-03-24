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
    # for 1st dataset
    # default row_to_read_file_gene_first_dataset = 22283
    row_to_read_file_gene_first_dataset = 22283
    file_gene_first_dataset_name = "GSE2034-22071 (edited).csv"
    file_gene_first_dataset = pd.read_csv(file_gene_first_dataset_name, nrows = row_to_read_file_gene_first_dataset)

    file_output_first_dataset_name = "mapping_sample_to_class_gse2034.csv"
    file_output_first_dataset = pd.read_csv(file_output_first_dataset_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # for 2nd dataset
    # default row_to_read_file_second_dataset = 22283
    row_to_read_file_gene_second_dataset = 22283
    file_gene_second_dataset_name = "GSE3494_GPL96.csv"
    file_gene_second_dataset = pd.read_csv(file_gene_second_dataset_name, nrows = row_to_read_file_gene_second_dataset) 

    file_output_second_dataset_name = "mapping_sample_to_class_gse3494.csv"
    file_output_second_dataset = pd.read_csv(file_output_second_dataset_name, usecols = ["SAMPLE_ID", "relapse (1=True)"])
    
    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    rows_to_read_file_pathway = 1329
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
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
    for element in sample_relapse_second_dataset.loc[:, 'SAMPLE_ID']:
        list_sample_relapse_second_dataset.append(element)
    
    list_sample_no_relapse_second_dataset = []
    for element in sample_no_relapse_second_dataset.loc[:, 'SAMPLE_ID']:
        list_sample_no_relapse_second_dataset.append(element)

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

    # # get output file's name
    # file_name = input("Name of output file : ")

    # # prepare text file for results to be written in
    # result_file = open(str(file_name) + ".txt", "w+")

    # list used to collect average auc score of each epoch
    list_avg_auc_each_epoch = []

    for epoch_count in range(0, num_of_epochs):
        print("######################################### epoch : " + str(epoch_count + 1) + "#########################################")
        # result_file.write("######################################### epoch : " + str(epoch_count + 1) + "#########################################\n")

        # conduct feature selection on 1st dataset
        # create list of indexes used to indicate the position in the list
        list_index_samples_relapse_first_dataset = []
        list_index_samples_no_relapse_first_dataset = []

        for index in range(0, len(list_sample_relapse_first_dataset)):
            list_index_samples_relapse_first_dataset.append(index)
        
        for index in range(0, len(list_sample_no_relapse_first_dataset)):
            list_index_samples_no_relapse_first_dataset.append(index)

    

if __name__ == "__main__":
    main()
