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
    chunk_relapse_size = math.ceil(len(list_sample_relapse) / num_of_folds)
    chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse) / num_of_folds)
    # print(chunk_size)    

    chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_relapse_size))
    print("# chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

    chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_no_relapse_size))
    print("# chunks in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

    check_valid = False
    num_of_chunks = None
    if (len(chunk_list_relapse) == len(chunk_list_no_relapse)):
        check_valid = True
        num_of_chunks = len(chunk_list_relapse)
    else:
        print("# chunks in 1 st set is not equal to # chunks in 2nd")

    print(chunk_list_relapse)
    # print(len(file_training_input.columns))

    # do only if number of chunks is equal
    if (check_valid == True):
        for first_layer_test_index in range(0, num_of_chunks):
            # keep testing data from each class
            first_layer_test_relapse = chunk_list_relapse[first_layer_test_index]
            first_layer_test_no_relapse = chunk_list_no_relapse[first_layer_test_index]
            print("\nINDEX : " + str(first_layer_test_index))
            print("test relapse =" + str(first_layer_test_relapse))
            print("test no relapse = " + str(first_layer_test_no_relapse))
            print()
            # find training data
            # first layer 
            first_layer_train_relapse = []
            for first_layer_train_index in range(0, num_of_chunks):
                if (chunk_list_relapse[first_layer_train_index] is not first_layer_test_relapse):
                    first_layer_train_relapse.append(chunk_list_relapse[first_layer_train_index])
            print("1st layer train relapse size = " + str(len(first_layer_train_relapse)))
            print("1st layer train relapse = " + str(first_layer_train_relapse))

            first_layer_train_no_relapse = []
            for first_layer_train_index in range(0, num_of_chunks):
                if (chunk_list_no_relapse[first_layer_train_index] is not first_layer_test_no_relapse):
                    first_layer_train_no_relapse.append(chunk_list_no_relapse[first_layer_train_index])
            print("1st layer train no relapse size = " + str(len(first_layer_train_no_relapse)))
            print("1st layer train no relapse = " + str(first_layer_train_no_relapse))

            # merge all element in each list to be used in second layer
            print("\nmerge element in each list to be used in the next step")
            second_list_sample_relapse = []
            for i in range(0, len(first_layer_train_relapse)):
                second_list_sample_relapse.extend(first_layer_train_relapse[i])
            print("size of list sample relapse = " + str(len(second_list_sample_relapse)))
            print("list sample relapse for next step = " + str(second_list_sample_relapse))

            second_list_sample_no_relapse = []
            for i in range(0, len(first_layer_train_no_relapse)):
                second_list_sample_no_relapse.extend(first_layer_train_no_relapse[i])
            print("size of list sample no relapse  = " + str(len(second_list_sample_no_relapse)))
            print("list sample no relapse for next step = " + str(second_list_sample_no_relapse))

           

if __name__ == '__main__':
    main()