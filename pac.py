import pandas as pd
import random
import math
import calculate

def main():
    # prepare data
    # row_to_read = 22283
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    rows_to_read_file_pathway = 1329

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)

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
            print("WARNING : Invalid input must be numeric")
        elif(int(num_of_folds) > len(list_sample_relapse)):
            print("WARNING : Number of folds exceeds the size of the 1st dataset")
        elif(int(num_of_folds) > len(list_sample_no_relapse)):
            print("WARNING : Number of folds exceeds the size of the 2nd dataset")
        elif(int(num_of_folds) <= 1):
            print("WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)

    # split data into k parts
    chunk_relapse_size = math.ceil(len(list_sample_relapse) / num_of_folds)
    chunk_no_relapse_size = math.ceil(len(list_sample_no_relapse) / num_of_folds)

    chunk_list_relapse = list(calculate.chunks(list_sample_relapse, chunk_relapse_size))
    print("# chunks in chunk_list_relapse = " + str(len(chunk_list_relapse)))

    chunk_list_no_relapse = list(calculate.chunks(list_sample_no_relapse, chunk_no_relapse_size))
    print("# chunks in chunk_list_no_relapse  = " + str(len(chunk_list_no_relapse)))

    check_valid, num_of_chunks = calculate.checkEqualListSize(chunk_list_relapse, chunk_list_no_relapse)

    # do only if number of chunks of both datasets are equal
    if (check_valid == True):
        for chunk_test_index in range(0, num_of_chunks):
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
            
            check_train_valid, num_of_chunks_train = calculate.checkEqualListSize(chunk_train_relapse, chunk_train_no_relapse)

            # get pathway activity of each sample
            # get column name to be read
            if (check_train_valid is True):
                for chunk_train_index in range(0, len(chunk_train_relapse)):
                    # collection of samples containing pathways of each sample
                    samples = {}

                    # identify columns to be read in each chunk in training data
                    for element_index in range(0, len(chunk_train_relapse[chunk_train_index])):
                        sample = []
                        sample_name = chunk_train_relapse[chunk_train_index][element_index]
                        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

                        sample.append(sample_name)
                        sample.append(pathways)
                        samples[element_index] = sample

if __name__ == '__main__':
    main()