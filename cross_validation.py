#!/usr/bin/python
import pandas as pd
import math
import calculate
import random

def main():
    # prepare data
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output = pd.read_csv("mapping_sample_to_class.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])

    # separate data into 2 classes
    sample_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['0'])]
    # print(sample_no_relapse)
    
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
            print("WARNING : Invalid input must be numeric")
        elif(int(num_of_folds) > len(list_sample_relapse)):
            print("WARNING : Number of folds exceeds the size of the 1st dataset")
        elif(int(num_of_folds) > len(list_sample_no_relapse)):
            print("WARNING : Number of folds exceeds the size of the 2nd dataset")
        elif(int(num_of_folds) < 1):
            print("WARNING : Number of folds cannot lower than 1")
        else:
            break
    num_of_folds = int(num_of_folds)

    # split data into k parts
    chunk_size = math.ceil(len(list_sample_relapse) / num_of_folds)
    # print(chunk_size)    

    chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_size))
    print("chunk_list_relapse SIZE = " + str(len(chunk_list_relapse)))

    chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_size))
    print("chunk_list_no_relapse SIZE = " + str(len(chunk_list_no_relapse)))

    # print(len(file_training_input.columns))
    # print(len(file_training_input.columns))

if __name__ == '__main__':
    main()