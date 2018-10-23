#!/usr/bin/python
import pandas as pd
import numpy as np
import calculate 

def main():
    # these column names can be retrieved by <code> list(my_dataframe.columns.values) </code>
    cols_to_read = ['GSM36778','GSM36784', 'GSM36789', 'GSM36792', 'GSM36797', 'GSM36811', 'GSM36814']
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", usecols = cols_to_read, nrows = 20)
    file_training_output = pd.read_csv("mapping_sample_to_class.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])
    
    # file_training_output.set_index("GEO asscession number", inplace=True)
    # data = data.loc['GSM36778']
    training_input = file_training_input
    training_output = file_training_output.loc[file_training_output['GEO asscession number'].isin(cols_to_read)]
    print("-------------------------Data to be used------------------------")
    print(training_input)
    print(training_output)
    print()

    # create list preparing for creating matrix
    # create list of output classes
    list_training_output = []
    for element in training_output.loc[:, 'relapses within 5 years (1 = yes, 0=no)']:
        list_training_output.append(element)

    # create list of input classes
    list_training_relapse = []
    list_training_no_relapse = []
    list_training_input = []

    # list containing all input values
    for column in training_input:
        list_each_sample = []
        for element in training_input[column]:
            list_each_sample.append(element)
        list_training_input.append(list_each_sample)

    # list input for each class
    count = 0
    for column in training_input:
        list_each_sample = []
        if (list_training_output[count] == 0):
            for element in training_input[column]:
                list_each_sample.append(element)
            list_training_no_relapse.append(list_each_sample) 
        elif(list_training_output[count] == 1):
            for element in training_input[column]:
                list_each_sample.append(element)
            list_training_relapse.append(list_each_sample) 
        count += 1
    print("-------------- List of Gene expression of all input -------------")
    print(list_training_input)
    print("------------- List of Gene expression of each class -------------")
    print("Relapse within 5 years ...")
    print(list_training_relapse)
    print("NO relapse within 5 years")
    print(list_training_no_relapse)
    print("--------------------- List of output class ----------------------")
    print(list_training_output)
    print("-----------------------------------------------------------------")
    print()


    # create matrix using for calculating 
    matrix_training_input = np.matrix(list_training_input)
    matrix_training_relapse = np.matrix(list_training_relapse)
    matrix_training_no_relapse = np.matrix(list_training_no_relapse)
    matrix_training_output = np.matrix(list_training_output).transpose()
    # print(matrix_training_output.transpose(1,0))
    print("---------------------- Matrix for all input ---------------------")
    print(matrix_training_input)
    print("-------------------- Matrix for each feature --------------------")
    print("Relapse within 5 years ... ")
    print(matrix_training_relapse)
    print("NO relapse within 5 years")
    print(matrix_training_no_relapse)
    print("-------------------- Matrix for output class --------------------")
    print(matrix_training_output)
    print("-----------------------------------------------------------------")
    print()

    # calculate average of each feature
    avg_training_input = calculate.avgFromList(matrix_training_input)
    avg_matrix_relapse = calculate.avgFromList(matrix_training_relapse)
    avg_matrix_no_relapse = calculate.avgFromList(matrix_training_no_relapse)
    print("average matrix training input : ")
    print(avg_training_input)
    print("average matrix relapse : ")
    print(avg_matrix_relapse)
    print("average matrix no relapse : ")
    print(avg_matrix_no_relapse)

if __name__ == '__main__':
    main()