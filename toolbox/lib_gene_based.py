from scipy import stats
from sklearn.metrics import roc_auc_score
from copy import deepcopy

import math
# import calculate
import random
import time
import files_handler
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
    if files_handler.validateFile(file_training_input_name) is False :
        sys.exit(0)

    # number of rows of this file to be read
    row_to_read = cf.gene_based["row_to_read"]

    # name of a file containing mapping between samples and their class
    file_training_output_name = cf.gene_based["file_training_output_name"]
    # check whether this file is valid
    if files_handler.validateFile(file_training_output_name) is False :
        sys.exit(0)

    # acquire required information to conduct an experiment
    # number of epochs
    epoch = cf.gene_based["epoch"]

    # number of folds
    num_of_folds = cf.gene_based["num_of_folds"]

    # number of top-ranked features
    number_of_ranked_gene = cf.gene_based["number_of_ranked_gene"]

