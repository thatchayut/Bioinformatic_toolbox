from scipy import stats
from sklearn.metrics import roc_auc_score
from copy import deepcopy

import math
# import calculate
import random
import time
import data_handler
import pandas as pd
import numpy as np
import sys

import config as cf

def test():
    print("testestest")

def gene_based():
    # record start time
    start_time = time.time()

    # acquire required files
    # a file containing mapping between probes IDs and sample
    file_training_input_name = cf.gene_based["file_training_input_name"]
    # check whether this file is valid
    if data_handler.validateFile(file_training_input_name) is False :
        sys.exit(0)

    # number of rows of this file to be read
    row_to_read = cf.gene_based["row_to_read"]
    # check whether number of rows is valid
    if data_handler.validateRowToRead(row_to_read) is False :
        sys.exit(0)

    # name of a file containing mapping between samples and their class
    file_training_output_name = cf.gene_based["file_training_output_name"]
    # check whether this file is valid
    if data_handler.validateFile(file_training_output_name) is False :
        sys.exit(0)

    # acquire required information to conduct an experiment
    # number of epochs
    epoch = cf.gene_based["epoch"]

    # number of folds
    num_of_folds = cf.gene_based["num_of_folds"]

    # number of top-ranked features
    number_of_ranked_gene = cf.gene_based["number_of_ranked_gene"]

    # prepare data

    # row_to_read = 22283
    # file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_input = pd.read_csv(file_training_input_name, nrows = row_to_read)

    # consider non-relapse and relapse (not in specific period of time)
    file_training_output = pd.read_csv(file_training_output_name, usecols = ['GEO asscession number', 'relapse (1=True)'])

    # this will be used in calculating lda
    training_input = file_training_input
    training_output = file_training_output

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    
    # separate data into 2 classes
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
    random.shuffle(list_sample_no_relapse)

