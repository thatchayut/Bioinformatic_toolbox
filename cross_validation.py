#!/usr/bin/python
import pandas as pd
import math
import calculate
import random

def main():
    # prepare data
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)

    # list of all sample names to be used in k-fold crossvalidation
    list_samples = file_training_input.columns.tolist()

    # shuffle data to make each chunk does not depend on sample order
    random.shuffle(list_samples)
    print(len(list_samples))

    # get number of folds
    while True:
        num_of_folds = input("Number of folds: ")
        if ((num_of_folds.isnumeric() == False) or (int(num_of_folds) > len(list_samples)) or (int(num_of_folds) < 1)):
            print("Invalid input...")
        else:
            break
    num_of_folds = int(num_of_folds)

    chunk_size = math.ceil(len(list_samples) / num_of_folds)
    print(chunk_size)

    # split data into k parts
    chunk_input = list(calculate.chunks(list_samples, chunk_size))
    print(len(chunk_input))
    
    # print(len(file_training_input.columns))
    # print(len(file_training_input.columns))

if __name__ == '__main__':
    main()